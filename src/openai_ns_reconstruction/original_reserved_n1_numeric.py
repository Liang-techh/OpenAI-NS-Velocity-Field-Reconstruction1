"""Actual reserved n1 velocities and own mass in fixed original source units.

Coefficient enclosures bound the original signed functions; their numerical
centers are proxies. Radius cancellation precedes every axial mass integral.
"""
from fractions import Fraction
import gzip
import hashlib
import json
import math

import mpmath as mp

from .original_core_numeric import ROOT,PREFIX,_bounds


class OriginalReservedN1Numeric:
    names=('alpha1','alpha2','beta1','beta2','beta3')

    def __init__(self,background,*,cells=256):
        if type(cells) is not int or not 32<=cells<=2048:raise ValueError('Reserved primitive quadrature cells32..2048 required')
        self.background=background;self.ctx,self.iv,self.J=background.ctx,background.iv,background.J
        self.family=background.family;self.directory=background.directory
        self.hashes=dict(background.hashes);self.code_hashes=dict(background.code_hashes)
        self.cells=cells;self._integrals={};self._packets={}
        from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_jets
        self.beta_jets=beta_jets
        specifications=(('current_original_n1_reserved_coefficient_enclosures',('current_original_n1_actual_reserved_coefficient_C1_numeric_enclosures_installed',)),
            ('current_original_n1_reserved_linear_operator',('current_original_n1_reserved_normalized_linear_repair_operator_installed',)),
            ('current_original_n1_terminal_integral_identities',('current_original_n1_actual_terminal_integrals_reduce_to_reserved_response','n1_repaired_own_terminal_five_function_zeros')))
        raw={}
        for stem,gates in specifications:
            receipt_path=self.directory/(PREFIX+'compliant_'+stem+'_check.json');receipt=json.loads(receipt_path.read_text(encoding='utf-8'))
            if not receipt['all_passed'] or receipt['source_family']!=self.family or not all(receipt[k] for k in gates):raise ValueError('Checked actual n1 coefficient/support/terminal source required: '+stem)
            for n,h in receipt['input_hashes'].items():
                if self.hashes.get(n,h)!=h or hashlib.sha256((self.directory/n).read_bytes()).hexdigest()!=h:raise ValueError('Changed/conflicting reserved n1 input: '+n)
                self.hashes[n]=h
            self.hashes[receipt_path.name]=hashlib.sha256(receipt_path.read_bytes()).hexdigest()
            p=self.directory/(PREFIX+'compliant_'+stem+'.json.gz');self.hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
            raw[stem]=json.loads(gzip.decompress(p.read_bytes()))
            if raw[stem]['source_family']!=self.family:raise ValueError('Same actual n1 family required')
        coefficients=raw['current_original_n1_reserved_coefficient_enclosures'];operator=raw['current_original_n1_reserved_linear_operator']
        if not coefficients['source_identity']['passed'] or not coefficients['source_identity']['actual_reserved_radius_identity']['passed']:raise ValueError('Actual reserved scale recipe identity required')
        self.coefficients=coefficients['samples']['whole_Z'];self.bindings=coefficients['true_function_bindings']
        self.support_certificate=operator['support_certificate'];self.matrix=operator['actual_normalized_matrix_enclosures']
        self.normal=self._decode(self.matrix['base_normalizer']);self.L=self.iv.log(2);self.ell=self.L/400
        self.centers=[self.L*j/40 for j in (12,13,14,15,16)]
        self.H0=self._decode(self.matrix['axial_divided_rows'][0][0])/self.iv.exp(self.centers[0]/2)
        self.mass_terminal_certificate=dict(archive=PREFIX+'compliant_current_original_n1_terminal_integral_identities.json.gz',
            gates=list(specifications[-1][1]),C1_mass_zero_at_2Rc=True,terminal_integrals_reduced_to_same_selected_response=True)
        self.tail=self.iv.exp(-64*background.patch.leading.defect.logP)
        p=self.directory/(PREFIX+'compliant_current_original_Rp_selected_pulse.json');selected=json.loads(p.read_text(encoding='utf-8'))
        self.logUp=self.iv.log(self._decode(selected['actual_current_interval_constants']['U']))
        self.hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
        name='src/openai_ns_reconstruction/original_reserved_n1_numeric.py';self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _decode(self,row):
        return self.iv.mpf([mp.make_mpf(tuple(row[k+'_exact_mpf_tuple'])) for k in ('lower','upper')])

    def _source_row(self,row):return self.iv.mpf([mp.make_mpf(tuple(v)) for v in row['source_enclosure_exact_mpf_tuples']])

    def _mid(self,v):return sum(_bounds(self.ctx,v),self.ctx.mpf(0))/2

    def _rows(self,values):
        ref=self.background.patch.leading.reference
        return [ref._summary(self._mid(v),v) for v in values]

    def _alpha(self,cells,name):
        c,iv=self.ctx,self.iv;pieces=[];proof=[];lm=self.background.outer.logmui
        for cell in cells:
            values=[]
            for k in (0,1):
                total=iv.mpf(0)
                for term in cell['coefficients'][name]['signed_original_source_terms']:
                    row=term['rows'][k]
                    if row['exact_zero']:continue
                    value=self._decode(term['transformed_ordinary_coefficients'][k])
                    log=self._decode(term['original_absolute_source_log'])+self._decode(term['exact_physical_log_shift'])-lm
                    bound=self._decode(row['absolute_log_upper'])-lm
                    if _bounds(c,bound)[1]<=_bounds(c,-64*self.background.patch.leading.defect.logP)[0]:
                        lo,hi=_bounds(c,value);cap=_bounds(c,self.tail)[1]
                        total+=iv.mpf([0,cap]) if lo>=0 else iv.mpf([-cap,0]) if hi<=0 else iv.mpf([-cap,cap])
                        proof.append(dict(coefficient=name,Z_order=k,sector=term['source_exponents'],admitted_log_cap='-64logPstar',actual_log_bound=bound))
                    else:
                        # Only representable original inner sectors reach this
                        # branch. Hyper-exponential exponent integers are not made.
                        lo,hi=_bounds(c,log)
                        if max(abs(lo),abs(hi))>c.mpf('1e30'):raise ArithmeticError('Unmaterializable n1 coefficient sector requires a proved small cap')
                        total+=value*iv.exp(log)
                values.append(total)
            pieces.append(values)
        return [iv.mpf([min(_bounds(c,v[k])[0] for v in pieces),max(_bounds(c,v[k])[1] for v in pieces)]) for k in (0,1)],proof

    def _integral_to(self,r,sign):
        """int_-1^r beta(t)*exp(sign*ell*t)dt, Taylor3/fourth remainder."""
        c,iv,J=self.ctx,self.iv,self.J
        r=max(c.mpf(-1),min(c.mpf(1),r));key=tuple(r._mpf_),sign
        if key in self._integrals:return self._integrals[key]
        if r==-1:return iv.mpf(0)
        width=(iv.mpf(r)+1)/self.cells;total=iv.mpf(0);k=sign*self.ell
        for i in range(self.cells):
            left=-1+width*i;right=-1+width*(i+1);mid=-1+width*(iv.mpf(i)+iv.mpf('.5'))
            box=iv.mpf([max(c.mpf(-1),_bounds(c,left)[0]),min(c.mpf(1),_bounds(c,right)[1])])
            def jet(q):
                beta=J(iv,self.beta_jets(iv,q));ex=iv.exp(k*q)
                return beta*J(iv,[ex*k**j/math.factorial(j) for j in range(5)])
            at,whole=jet(mid),jet(box)
            total+=at[0]*width+at[2]*width**3/12+whole[4]*width**5/80
        lo,hi=_bounds(c,total);result=iv.mpf([max(0,lo),max(0,hi)])
        self._integrals[key]=result;return result

    def _shape_and_mass(self,s,j):
        c,iv=self.ctx,self.iv;r=(s-self.centers[j])/self.ell;lo,hi=_bounds(c,r)
        shape=self.beta_jets(iv,r)[0]/(self.ell*iv.exp(self.centers[j]/2)*self.normal)
        full=iv.exp(self.centers[j]/2)*self.H0
        if hi<=-1:return iv.mpf(0),iv.mpf(0),full
        if lo>=1:return iv.mpf(0),full,iv.mpf(0)
        factor=iv.exp(self.centers[j]/2)/self.normal
        partial=iv.mpf([_bounds(c,self._integral_to(lo,1))[0],_bounds(c,self._integral_to(hi,1))[1]])*factor
        # beta is even; the suffix is the same primitive with opposite
        # exponential weight and reflected endpoint, without full-partial loss.
        suffix=iv.mpf([_bounds(c,self._integral_to(-hi,-1))[0],_bounds(c,self._integral_to(-lo,-1))[1]])*factor
        return shape,partial,suffix

    def _intersection(self,a,b):
        c,iv=self.ctx,self.iv;out=[]
        for x,y in zip(a,b):
            lo,hi=_bounds(c,x);ll,hh=_bounds(c,y)
            if max(lo,ll)>min(hi,hh):raise ArithmeticError('Actual reserved n1 mass sources have disjoint enclosures')
            out.append(iv.mpf([max(lo,ll),min(hi,hh)]))
        return out

    def beta_unit_atoms(self,packet):
        """log(Utheta1/normalized_shape), from exact pressure-unit recipe.

        Unit_beta + log(sqrt(2/Rc)) = logPressure15-logP-logSigma-logRc.
        Its inverse-mu powers are combined before any interval boxing.
        """
        c,iv=self.ctx,self.iv;ref=self.background.patch.leading.reference
        lp=self.background.patch.leading.defect.logP;lm=self.background.outer.logmui
        phi=iv.mpf('10.99')*iv.sqrt(2*iv.pi)*(iv.exp(5)-iv.exp(3))
        atoms=[self._source_row(v) for v in packet['coordinate']['global_y_log_terms'][:4]]
        # Those four terms are log(Rc/110); local log(x) is deliberately absent.
        atoms.extend([iv.log(110),iv.log(phi)+self.logUp-iv.mpf('17.5'),lp/2,
            iv.mpf('18.5')*iv.exp(-lm),-3*iv.exp(-iv.mpf(2)/3*(iv.log(4)+lm)),
            -iv.mpf(364)/3*lm,iv.mpf(2)/3*iv.log(4)-iv.log(6)/2])
        return [ref._summary(self._mid(v),v) for v in atoms]

    def evaluate(self,Z,x,packet=None):
        c,iv=self.ctx,self.iv
        if type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1<c.mpf(Z)<1:raise ValueError('Original interior decimal Z required')
        if type(x) is not str or not c.isfinite(c.mpf(x)) or not 1<=Fraction(x)<=2:raise ValueError('Original reserved n1 native1<=x<=2 required')
        key=Z,x
        if key in self._packets:return self._packets[key]
        packet=packet or self.background.evaluate(Z,x);z=iv.mpf(Z);zl,zh=_bounds(c,z)
        cells=[v for v in self.coefficients['cells'] if not (_bounds(c,self._decode(v['clipped_Z']))[1]<zl or _bounds(c,self._decode(v['clipped_Z']))[0]>zh)]
        if not cells:raise ValueError('Actual original n1 coefficient cover required')
        alpha=[];proof=[]
        for name in self.names[:2]:
            values,p=self._alpha(cells,name);alpha.append(values);proof.extend(p)
        beta=[]
        for name in self.names[2:]:
            beta.append([iv.mpf([min(_bounds(c,self._decode(v['coefficients'][name]['normalized_C0_C1_enclosures'][k]))[0] for v in cells),
                max(_bounds(c,self._decode(v['coefficients'][name]['normalized_C0_C1_enclosures'][k]))[1] for v in cells)]) for k in (0,1)])
        s=iv.log(iv.mpf(x));shapes=[self._shape_and_mass(s,j) for j in range(5)]
        axial=[sum((a[k]*shapes[j][0] for a,j in zip(alpha,(0,4))),iv.mpf(0)) for k in (0,1)]
        swirl=[sum((b[k]*shapes[j][0] for b,j in zip(beta,(1,2,3))),iv.mpf(0)) for k in (0,1)]
        history=packet['incoming_n1_history_not_extended'];primitive=history['actual_terminal_primitive']['actual_primitive']
        ll=iv.mpf(self.background.owner.amplitude.ll_iv)
        incoming=[self._source_row(v)*iv.exp(-2*ll) for v in primitive[:2]]
        forward=[incoming[k]+sum((a[k]*shapes[j][1] for a,j in zip(alpha,(0,4))),iv.mpf(0)) for k in (0,1)]
        backward=[-sum((a[k]*shapes[j][2] for a,j in zip(alpha,(0,4))),iv.mpf(0)) for k in (0,1)]
        mass=self._intersection(forward,backward)
        RUz=[v*iv.exp(s) for v in axial];delta=iv.mpf(self.background.owner.amplitude.delta_iv)
        RQ=(2*z*RUz[0]-(1+delta)*z*mass[0]-(1-z*z)*mass[1])/(1-delta*z*z)
        active=[j for j in range(5) if _bounds(c,shapes[j][0])[1]>0]
        grouped={}
        for v in proof:grouped.setdefault((v['coefficient'],v['Z_order']),[]).append(v['actual_log_bound'])
        summaries=[]
        for (name,k),values in grouped.items():
            worst=iv.mpf(max(_bounds(c,v)[1] for v in values))
            summaries.append(dict(coefficient=name,Z_order=k,source_sectors=len(values),admitted_log_cap='-64logPstar',
                largest_actual_log_bound=self.background.patch.leading.reference._summary(self._mid(worst),worst)))
        out=dict(Z=Z,native_compact_x=x,source_family=self.family,
            alpha_over_mu_C0_C1=[self._rows(v) for v in alpha],beta_normalized_C0_C1=[self._rows(v) for v in beta],
            Uz1_times_Rc_C0_C1=self._rows(axial),Utheta1_normalized_C0_C1=self._rows(swirl),M1_C0_C1=self._rows(mass),
            R_times_Q1=self._rows([RQ]),R_times_Uz1_C0_C1=self._rows(RUz),incoming_M1_C0_C1=self._rows(incoming),
            forward_M1_C0_C1=self._rows(forward),backward_terminal_M1_C0_C1=self._rows(backward),
            mass_increment_factors=[dict(shape=self._rows([v[0]])[0],partial=self._rows([v[1]])[0],suffix=self._rows([v[2]])[0]) for v in shapes],
            Utheta1_log_unit_atoms=self.beta_unit_atoms(packet),actual_coefficient_function_bindings=self.bindings,
            support_certificate=self.support_certificate,active_bumps=active,terminal_mass_source_identity=self.mass_terminal_certificate,
            small_sector_proofs=summaries,
            contract=dict(actual_original_five_selected_coefficient_enclosures=True,inner_mass_J_over_Lambda_squared_retained=True,
                Rc_cancels_before_axial_mass_integration=True,coefficients_and_units_Z_independent_units_respected=True,
                physical_swirl_already_Utheta1_not_F1=True,positive_order_recovery_weight_one_plus_delta=True,
                same_actual_source_terminal_mass_identity=True,normalized_source_errors_not_selected_as_exact_functions=True),
            limits=dict(velocity_values=True,mass_C1=True,Q1_Z_derivative=False,n1_pressure=False,global_n1_joins=False,
                resolved_signed_coefficient_functions=False,beyond_2Rc=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out;return out
