"""Signed original core m15 approximations with full source C0/C1 bounds.

Finite Volterra forcing polynomials are integrated analytically. The source
and all-word errors are integrated separately; no n1 radial Taylor theorem
or enclosure midpoint is substituted for the original solution.
"""
import hashlib
import json

from .original_core_numeric import ROOT, PREFIX, _box, _bounds
from .original_n1_core_volterra_numeric import OriginalN1CoreVolterraNumeric


def c1_product_rows(left, right, ctx):
    """Exact ordinary-C1 product of two finite radial polynomials."""
    out = [[ctx.mpf(0), ctx.mpf(0)] for _ in range(len(left)+len(right)-1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i+j][0] += a[0]*b[0]
            out[i+j][1] += a[1]*b[0]+a[0]*b[1]
    return out


def signed_prefix_rows(leading_U, n1_u, n1_p, radius, ctx):
    """Return 2 int U0*u drho and int p drho-rho*p at the true endpoint."""
    product = c1_product_rows(leading_U, n1_u, ctx)
    axial = [2*sum((row[k]*radius**(n+1)/(n+1) for n,row in enumerate(product)), ctx.mpf(0)) for k in (0,1)]
    # Combine the exact finite endpoint with its integral coefficient first.
    # The constant pressure row cancels, without interval subtraction loss.
    pressure = [-sum((row[k]*radius**(n+1)*n/(n+1) for n,row in enumerate(n1_p)), ctx.mpf(0)) for k in (0,1)]
    return axial, pressure


class OriginalN1CoreM15Numeric:
    def __init__(self, *, solver=None, **kwargs):
        if solver is not None and kwargs:
            raise ValueError('Supply an existing core owner or constructor options')
        self.n1 = solver if solver is not None else OriginalN1CoreVolterraNumeric(**kwargs)
        self.leading = self.n1.leading
        self.ctx,self.iv,self.family = self.n1.ctx,self.n1.iv,self.n1.family
        self.directory = self.leading.directory
        self.hashes,self.code_hashes = dict(self.n1.hashes),dict(self.n1.code_hashes)
        stem = PREFIX+'compliant_current_original_n1_cumulative_moments'
        receipt_path = self.directory/(stem+'_check.json')
        receipt = json.loads(receipt_path.read_text(encoding='utf-8'))
        if (not receipt['all_passed'] or receipt['source_family'] != self.family
                or not receipt['current_original_n1_common_inner_five_cumulative_moment_enclosures_installed']):
            raise ValueError('Accepted same-source original cumulative moments required')
        for name,digest in receipt['input_hashes'].items():
            if self.hashes.get(name,digest) != digest:
                raise ValueError('Conflicting signed core moment source: '+name)
            if name not in self.hashes and hashlib.sha256((self.directory/name).read_bytes()).hexdigest() != digest:
                raise ValueError('Changed signed core moment source: '+name)
            self.hashes[name]=digest
        self.hashes[receipt_path.name]=hashlib.sha256(receipt_path.read_bytes()).hexdigest()
        for name in ('src/openai_ns_reconstruction/original_n1_core_volterra_numeric.py',
                     'src/openai_ns_reconstruction/original_n1_core_m15_numeric.py'):
            self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        self._packets={}

    def _rows(self, point, truth):
        out=[]
        for p,t in zip(point,truth):
            row=self.n1._summary(p,t)
            row['source_enclosure_exact_mpf_tuples']=row['original_source_enclosure_exact_mpf_tuples']
            out.append(row)
        return out

    def _symmetric(self, value):
        hi=_bounds(self.ctx,value)[1]
        if hi<0:raise ArithmeticError('A moment error majorant must be nonnegative')
        return self.iv.mpf([-hi,hi])

    def _uniform_errors(self, data, Z, radius, degree):
        """Pointwise |u-u_poly|<=rho*eu, |p-p_poly|<=rho*ep on[0,r]."""
        iv=self.iv
        errors,_=self.n1._source_tail(data,Z,radius,degree)
        forcing_u=[v/4 for v in errors['u'][:2]]
        forcing_p=errors['p'][:2]
        feedback=[_box(iv,self.n1.tail[key]) for key in ('value_remainder_upper','ordinary_Z1_remainder_upper')]
        h=_box(iv,self.leading.tube['Xh_parameter'])
        correction=_box(iv,self.leading.major['scaled_map_size_upper'])/(1-_box(iv,self.leading.major['scaled_map_Lipschitz_upper']))
        leading=[]
        for k in (0,1):
            # Actual leading U=axis+epsilon*Psi. The analytic correction's
            # tail, not its whole ball, is multiplied by epsilon once.
            leading.append(data['e']*correction*self.leading.tail_factor(iv,degree=degree,
                radial_order=0,axial_order=k,radius=iv.mpf(radius),h=h)['tail_per_Xh_norm'])
        return forcing_u,forcing_p,feedback,leading

    def evaluate(self,Z,rho,*,degree=16,kernel_degree=24,feedback_iterations=1):
        c,iv=self.ctx,self.iv
        if (type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1<c.mpf(Z)<1
                or type(rho) is not str or not c.isfinite(c.mpf(rho)) or not 0<=c.mpf(rho)<=4
                or type(degree) is not int or not 4<=degree<=21
                or type(kernel_degree) is not int or not 4<=kernel_degree<=32
                or type(feedback_iterations) is not int or feedback_iterations not in (0,1)):
            raise ValueError('Original core Z(-1,1), rho[0,4], degree4..21, kernel4..32 and iterations0/1 required')
        key=Z,rho,degree,kernel_degree,feedback_iterations
        if key in self._packets:return self._packets[key]
        rp,rt=c.mpf(rho),iv.mpf(rho)
        ledger={}
        if rp==0:
            point=([c.mpf(0)]*2,[c.mpf(0)]*2)
            truth=([iv.mpf(0)]*2,[iv.mpf(0)]*2)
        else:
            total,_,data,_,principal=self.n1._prepare(Z,degree,kernel_degree,feedback_iterations)
            _,nominal=self.n1._forcing_polynomial(Z,degree,False)
            point=signed_prefix_rows(nominal['rows']['Uz'],total['u'].rows,total['p'].rows,rp,c)
            base=signed_prefix_rows(data['rows']['Uz'],principal['u'].rows,principal['p'].rows,rt,iv)
            fu,fp,feedback,lead_error=self._uniform_errors(data,Z,rho,degree)
            if any(_bounds(c,v)!=(0,0) for v in principal['u'].rows[0][:2]):
                raise ValueError('The original n1 core must have its regular zero axis value')
            magnitude=lambda v:iv.mpf(max(abs(t) for t in _bounds(c,v)))
            U=[sum((magnitude(row[k])*rt**n for n,row in enumerate(data['rows']['Uz'])),iv.mpf(0)) for k in (0,1)]
            u_over_r=[sum((magnitude(row[k])*rt**(n-1) for n,row in enumerate(principal['u'].rows) if n),iv.mpf(0)) for k in (0,1)]

            def product_error(eu,include_leading):
                # Integration of rho gives r^2/2 and the density has factor2.
                out=[U[0]*eu[0]+lead_error[0]*eu[0],
                     U[1]*eu[0]+lead_error[1]*eu[0]+U[0]*eu[1]+lead_error[0]*eu[1]]
                if include_leading:
                    out[0]+=u_over_r[0]*lead_error[0]
                    out[1]+=u_over_r[0]*lead_error[1]+u_over_r[1]*lead_error[0]
                return [rt*rt*v for v in out]

            forcing_axial=product_error(fu,False)
            feedback_axial=product_error(feedback,False)
            leading_axial=product_error([iv.mpf(0)]*2,True)
            forcing_pressure=[iv.mpf('1.5')*rt*rt*v for v in fp]
            feedback_pressure=[iv.mpf('1.5')*rt*rt*v for v in feedback]
            ea=[forcing_axial[k]+feedback_axial[k]+leading_axial[k] for k in (0,1)]
            ep=[forcing_pressure[k]+feedback_pressure[k] for k in (0,1)]
            truth=([base[0][k]+self._symmetric(ea[k]) for k in (0,1)],
                   [base[1][k]+self._symmetric(ep[k]) for k in (0,1)])
            for name,values in dict(axial_forcing=forcing_axial,axial_all_feedback=feedback_axial,
                    axial_leading_analytic_tail=leading_axial,pressure_forcing_and_finite_endpoint=forcing_pressure,
                    pressure_all_feedback_and_finite_endpoint=feedback_pressure).items():
                ledger[name]=self._rows([c.mpf(0)]*2,[self._symmetric(v) for v in values])
        lp,lpiv=self.leading.logP,self.leading.lp_iv
        ll,lliv=self.leading.logLambda,self.leading.ll_iv
        # R=epsilon*rho, U1=epsilon*u, P1=epsilon*Pstar^2*p.
        scale_e=self._rows([-2*ll],[-2*lliv])
        scale_p=scale_e+self._rows([2*lp],[2*lpiv])
        terms=[dict(name=name,coefficient_C0_C1=self._rows(point[j],truth[j]),log_scale_terms=atoms)
            for j,(name,atoms) in enumerate((('core_axial_coupling',scale_e),
                ('core_pressure_integral_and_finite_boundary',scale_p)))]
        out=dict(Z=Z,rho=rho,source_family=self.family,m15_core_C0_C1_terms=terms,
            error_ledger=ledger,source_degree=degree,potential_kernel_degree=kernel_degree,
            nominal_feedback_iterations=feedback_iterations,
            source_formula=dict(source=PREFIX+'compliant_current_original_n1_cumulative_moments.py',
                functions=['rho_integrands','endpoint_cumulants'],
                actual_m15='epsilon^2*(2*integral U_lead*u drho)+epsilon^2*Pstar^2*(integral p drho-rho*p)',
                leading_axial='U_lead=U0+epsilon*Psi; both source sectors retained',
                finite_endpoint='-R*P1=-epsilon^2*Pstar^2*rho*p',
                coefficient_C1_convention='ordinary Z0 and Z1 at fixed original rho'),
            contract=dict(actual_leading_and_original_n1_Volterra_sources=True,
                finite_forcing_polynomials_integrated_analytically=True,
                leading_analytic_source_tail_integrated=True,
                full_all_word_feedback_tail_integrated=True,
                both_pressure_integral_and_finite_endpoint_errors_retained=True,
                ordinary_Z_product_derivatives_retained=True,
                nominal_feedback_distance_charged_against_principal_enclosure=True,
                serialized_decimal_approximation_has_complete_absolute_error=True,
                source_enclosure_midpoints_not_selected_as_values=True,
                physical_epsilon_and_Pstar_powers_applied_once=True),
            limits=dict(signed_core_m15_approximation_with_C0_C1_error=True,
                common_inner_rho_up_to_4=True,inner_cutoff_first_collar=False,lower_collar=False,
                complete_signed_inner_m15=False,active_compact_signed_pressure=False,
                independent_preRp_signed_numeric_seam=False,global_five_moment_numeric_join=False,
                common_global_N=False,exact_heat_exterior=False,stress_cone=False,
                genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out
        return out
