"""Original leading implicit Rm patch with source-owned Z2 control jets.

The leading five equations use unmodified defects. Finite-N loop histories
remain distinct incoming data for their subsequent original dx/x transport.
"""
from fractions import Fraction
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_Rm_defect_numeric import OriginalRmDefectNumeric


def encode(value):
    if hasattr(value,'coefficients'):return dict(Taylor_coefficients=[encode(v) for v in value.coefficients])
    if hasattr(value,'_mpi_'):return dict(lower_exact_mpf_tuple=list(value._mpi_[0]),upper_exact_mpf_tuple=list(value._mpi_[1]))
    if hasattr(value,'_mpf_'):return dict(exact_mpf_tuple=list(value._mpf_))
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [encode(v) for v in value]
    return value


def implicit_Z2(ctx,J,W,d,alpha):
    """Strict C0 contraction, its C1 inverse and the same-Jacobian Z2 row."""
    from lei_ren_part1_paper_interval_five_bump_inverse import (
        certify,nonlinear_jacobian,identity_error,weighted_norm,dot,mag,intersection)
    from lei_ren_part1_paper_compliant_actual_moment_patch import nonlinear_series
    if any(v.order!=2 for v in d) or alpha.order!=2:raise ValueError('Actual Taylor Z2 source rows required')
    first=certify(ctx,W,[v.truncate(1) for v in d],alpha.truncate(1),scales=[ctx.mpf(1)]*5)
    if not first['certified']:raise ArithmeticError('Original leading patch strict self-map failed')
    mid=mp.matrix([[(mp.make_mpf(v._mpi_[0])+mp.make_mpf(v._mpi_[1]))/2 for v in row] for row in W['L']])
    inverse=mp.inverse(mid);R=[[ctx.mpf(inverse[i,j]) for j in range(5)] for i in range(5)]
    h0=[v[0] for v in first['controls']]
    jq=nonlinear_jacobian(ctx,W,h0,alpha[0])
    jac=[[W['L'][i][j]+jq[i][j] for j in range(5)] for i in range(5)]
    K=identity_error(ctx,R,jac);q=weighted_norm(ctx,K,[ctx.mpf(1)]*5)
    if mp.make_mpf(q._mpi_[1])>=1:raise ArithmeticError('Same-Jacobian second-row inverse failed')
    known=[J(ctx,list(v.coefficients)+[ctx.mpf(0)]) for v in first['controls']]
    Q=nonlinear_series(ctx,W,known,alpha)
    rhs=[d[i][2]+Q[i][2] for i in range(5)];b=[-dot(row,rhs) for row in R]
    radius=max((mag(ctx,v) for v in b),key=lambda v:mp.make_mpf(v._mpi_[1]))/(1-q)
    hi=mp.make_mpf(radius._mpi_[1]);h2=[ctx.mpf([-hi,hi]) for _ in range(5)]
    for _ in range(12):h2=[intersection(ctx,a,b[i]+dot(K[i],h2)) for i,a in enumerate(h2)]
    controls=[J(ctx,list(v.coefficients)+[h2[i]]) for i,v in enumerate(first['controls'])]
    nonlinear=nonlinear_series(ctx,W,controls,alpha)
    residual=[d[i]+nonlinear[i]+sum((controls[j]*W['L'][i][j] for j in range(5)),controls[0]*0)
              for i in range(5)]
    if any(not mp.make_mpf(v._mpi_[0])<=0<=mp.make_mpf(v._mpi_[1]) for row in residual for v in row.coefficients):
        raise ArithmeticError('Same implicit-map diagnostic excludes zero')
    return controls,dict(first=first,second_row_same_Jacobian=jac,second_row_contraction=q,
        second_row_known_rhs=rhs,second_row_inverse_radius=radius,unrefined_map_residual=residual,
        closure_from_strict_contraction_and_same_implicit_function=True,
        residual_zero_containment_is_only_diagnostic=True)


class _Bump:
    """Pinned numerical leaf provider; no old datum/ancestor is constructed."""
    def __init__(self,ctx,raw,constants):
        from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value
        self.ctx=ctx;self.normalization=restore_value(ctx,raw['normalization'])
        self.beta_sup=restore_value(ctx,constants['bump_sup_bound'])
        self.beta_deriv_sup=restore_value(ctx,constants['bump_first_derivative_sup_bound'])
        self.full_weights={(Fraction(v['center']),Fraction(v['power']),v['multiplicity']):
            restore_value(ctx,v['weight_interval']) for v in raw['weight_records'].values()}

    def beta(self,offset):
        from lei_ren_part1_paper_compliant_five_moment_repair import SharedFiveMomentRepair
        return SharedFiveMomentRepair.beta(self,offset)

    def partial_weight(self,x,center,power,multiplicity=1,cells=128):
        from lei_ren_part1_paper_compliant_five_moment_repair import SharedFiveMomentRepair
        return SharedFiveMomentRepair.partial_weight(self,x,center,power,multiplicity,cells)


def weight_matrix(ctx,weight):
    """Original five-bump map, including partial mass weights."""
    L=[[ctx.mpf(0) for _ in range(5)] for _ in range(5)]
    for i,k in enumerate((0,2)):L[0][i]=weight(k,'0');L[1][i]=weight(k,'.6')
    for i in range(3):
        L[2][i+2]=weight(i,'.5');L[3][i+2]=-weight(i,'.1');L[4][i+2]=weight(i,'-.9')
    return dict(L=L,fg=[weight(i,'.5',2) for i in (0,2)],gg=[weight(i,'0',2) for i in (0,2)],
        ff=[weight(i,'0',2) for i in range(3)],ff_over_x=[weight(i,'-1',2) for i in range(3)])


class OriginalRmPatchNumeric:
    centers=(Fraction(5,4),Fraction(3,2),Fraction(7,4))
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        self.defect=OriginalRmDefectNumeric(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        self.reference=self.defect.reference;self.ctx,self.iv,self.J=self.defect.ctx,self.defect.iv,self.defect.J
        self.family=self.defect.family;self.hashes=dict(self.defect.hashes);self.code_hashes=dict(self.defect.code_hashes)
        self.directory=self.reference.leading.core.directory;self.precision=precision
        def read(stem):
            path=self.directory/(PREFIX+stem+'.json');r=json.loads(path.read_text(encoding='utf-8'))
            digest=hashlib.sha256(path.read_bytes()).hexdigest()
            if self.hashes.get(path.name,digest)!=digest:raise ValueError('Changed actual patch source: '+stem)
            self.hashes[path.name]=digest
            return r
        self.accepted=read('compliant_current_original_Rm_defect_numeric')
        if not self.accepted['checks']['passed'] or self.accepted['source_family']!=self.family:
            raise ValueError('Accepted actual signed Rm source required')
        for name,digest in self.accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:raise ValueError('Changed signed inlet implementation')
            self.code_hashes[name]=digest
        self._inlets={v['Z']:v for v in self.accepted['samples']}
        repair=read('compliant_five_moment_repair');check=read('compliant_five_moment_repair_check')
        if not check['actual_implicit_functional_five_moment_closure_independently_checked']:
            raise ValueError('Original independent map/source admission required')
        branch=read('compliant_current_original_true_Rm_patch_C6_functions_check')
        if not branch['all_passed'] or branch['source_family']!=self.family or not branch['same_exact_original_unique_branch_theorem_bound']:
            raise ValueError('Same original leading branch required')
        raw=read('bump_integral_enclosures_check');constants=read('shared_bump_constants')
        if not raw['directed_integrals_certified']:raise ValueError('Original directed bump integrals required')
        scale=repair['axial_box_scale']
        if scale['lower_exact_mpf_tuple']!=scale['upper_exact_mpf_tuple']:raise ValueError('Fixed normalization coordinate required')
        t=mp.make_mpf(tuple(scale['lower_exact_mpf_tuple']));self.t=self.iv.mpf(t);self.tp=self.ctx.mpf(t)
        self.angular=self.iv.exp(-2*self.defect.logP)*(10**10*self.t*self.t)
        self.angularp=self.ctx.exp(-2*self.defect.logPp)*(10**10*self.tp*self.tp)
        self.units=[self.t]*2+[self.angular]*3;self.unitsp=[self.tp]*2+[self.angularp]*3
        if _bounds(self.ctx,self.angular)[0]<=0:raise ValueError('Strictly positive original angular coordinate required')
        from lei_ren_part1_paper_interval_five_bump_inverse import weights
        self.W=weights(self.iv,raw);self.bump=_Bump(self.iv,raw,constants)
        self._weights={};self._point_weights={};self._coefficients={};self._packets={}
        self.normp=self.ctx.quad(lambda s:self.ctx.exp(-1/(1-s*s)) if abs(s)<1 else self.ctx.mpf(0),[-1,0,1])
        for stem in ('compliant_five_moment_repair','compliant_actual_moment_patch','interval_five_bump_inverse',
                     'bump_integral_enclosures','compliant_current_original_true_Rm_patch_C6_functions'):
            path=self.directory/(PREFIX+stem+'.py');self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        name='src/openai_ns_reconstruction/original_Rm_patch_numeric.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _restore(self,row):return self.reference.previous._restore(row)
    def _jet(self,rows,directed):
        ctx=self.iv if directed else self.ctx
        return self.J(ctx,[self._restore(r) if directed else ctx.mpf(r['approximation']) for r in rows])

    def _inlet(self,Z):
        if Z not in self._inlets:self._inlets[Z]=self.defect.evaluate(Z)
        return self._inlets[Z]

    def _point_weight(self,x,i,power,m):
        c=self.ctx;cc=c.mpf(self.centers[i].numerator)/self.centers[i].denominator;r=c.mpf('.025')
        end=min(c.mpf(1),(x-cc)/r)
        if end<=-1:return c.mpf(0)
        if end==1 and m==1 and c.mpf(power)==0:return c.mpf(1)
        key=(str(end),i,power,m)
        if key not in self._point_weights:
            self._point_weights[key]=c.quad(lambda s:(cc+r*s)**c.mpf(power)*
                (c.exp(-1/(1-s*s))/(r*self.normp))**m*r if abs(s)<1 else c.mpf(0),[-1,min(c.mpf(0),end),end])
        return self._point_weights[key]

    def _W(self,x,directed,terminal):
        ctx=self.iv if directed else self.ctx;key=(str(x),directed,terminal)
        if key not in self._weights:
            def w(i,p,m=1):
                if directed:
                    if terminal:return self.bump.full_weights[self.centers[i],Fraction(p),m]
                    return self.bump.partial_weight(x,self.centers[i],p,m)
                return self._point_weight(x,i,p,m)
            self._weights[key]=weight_matrix(ctx,w)
        return self._weights[key]

    def coefficients(self,Z):
        if Z in self._coefficients:return self._coefficients[Z]
        inlet=self._inlet(Z);c,iv,J=self.ctx,self.iv,self.J
        d=[self._jet(r,True)/self.units[i] for i,r in enumerate(inlet['signed_unmodified_defects'])]
        dp=[self._jet(r,False)/self.unitsp[i] for i,r in enumerate(inlet['signed_unmodified_defects'])]
        z=J.variable(iv,iv.mpf(Z),2);zp=J.variable(c,c.mpf(Z),2)
        alpha=(1+z*z)*(1+z*z)*iv.exp(iv.mpf('1.2'))/(10**10)
        alphap=(1+zp*zp)*(1+zp*zp)*c.exp(c.mpf('1.2'))/(10**10)
        Wn={n:v if n in ('L','gg') else [q*self.angular for q in v] for n,v in self.W.items()}
        Wp=self._W(c.mpf(2),False,True)
        Wpn={n:v if n in ('L','gg') else [q*self.angularp for q in v] for n,v in Wp.items()}
        with mp.workdps(self.precision+50):controls,proof=implicit_Z2(iv,J,Wn,d,alpha)
        from lei_ren_part1_paper_interval_five_bump_inverse import nonlinear,nonlinear_jacobian
        from lei_ren_part1_paper_compliant_actual_moment_patch import nonlinear_series
        L=c.matrix(Wpn['L']);guess=c.lu_solve(L,c.matrix([-v[0] for v in dp]))
        def residual(*h):return [sum((Wpn['L'][i][j]*h[j] for j in range(5)),c.mpf(0))+
            nonlinear(c,Wpn,h,alphap[0])[i]+dp[i][0] for i in range(5)]
        h0=list(c.findroot(residual,tuple(guess),tol=c.mpf(10)**(-self.precision+12),maxsteps=30))
        jq=nonlinear_jacobian(c,Wpn,h0,alphap[0]);jac=c.matrix([[Wpn['L'][i][j]+jq[i][j] for j in range(5)] for i in range(5)])
        explicit=[c.mpf(0)]*5;explicit[3]=alphap[1]*sum((Wpn['gg'][i]*h0[i]*h0[i] for i in range(2)),c.mpf(0))
        h1=list(c.lu_solve(jac,c.matrix([-dp[i][1]-explicit[i] for i in range(5)])))
        known=[J(c,[h0[i],h1[i],0]) for i in range(5)];Q=nonlinear_series(c,Wpn,known,alphap)
        h2=list(c.lu_solve(jac,c.matrix([-dp[i][2]-Q[i][2] for i in range(5)])))
        point=[J(c,[h0[i],h1[i],h2[i]]) for i in range(5)]
        data=dict(inlet=inlet,point=[v*self.unitsp[i] for i,v in enumerate(point)],
            truth=[v*self.units[i] for i,v in enumerate(controls)],normalized_point=point,normalized_truth=controls,
            proof=proof,defects_point=[self._jet(r,False) for r in inlet['signed_unmodified_defects']],
            defects_truth=[self._jet(r,True) for r in inlet['signed_unmodified_defects']])
        self._coefficients[Z]=data
        return data

    def evaluate(self,Z,x):
        if type(x) is not str:raise ValueError('Native x decimal string or Rh required')
        c,iv,J=self.ctx,self.iv,self.J
        if x=='Rh':xp,xt=c.exp(1),iv.exp(1);terminal=True
        else:
            xp,xt=c.mpf(x),iv.mpf(x)
            if not c.isfinite(xp) or xp<1 or xp>c.exp(1):raise ValueError('Original 1<=x=R/Rm<=e required')
            terminal=Fraction(x)>=Fraction(71,40)
        key=(Z,x)
        if key in self._packets:return self._packets[key]
        data=self.coefficients(Z);sides=[]
        from lei_ren_part1_paper_compliant_actual_moment_patch import nonlinear_series
        for directed,ctx,xx in ((False,c,xp),(True,iv,xt)):
            h=data['truth' if directed else 'point'];d=data['defects_truth' if directed else 'defects_point']
            z=J.variable(ctx,ctx.mpf(Z),2);zero=J.constant(ctx,0,2)
            invP=ctx.exp(-ctx.mpf(self.defect.logP if directed else self.defect.logPp))
            am=(1+z*z).reciprocal()*ctx.exp(ctx.mpf('-.6'));alpha=(am*am).reciprocal()*(invP*invP)
            W=self.W if directed and terminal else self._W(xx,directed,terminal)
            Q=nonlinear_series(ctx,W,h,alpha)
            unrefined=[d[i]+Q[i]+sum((h[j]*W['L'][i][j] for j in range(5)),zero) for i in range(5)]
            if directed and terminal and any(not _bounds(c,v)[0]<=0<=_bounds(c,v)[1] for row in unrefined for v in row.coefficients):
                raise ArithmeticError('Actual full-weight implicit-map residual excludes zero')
            residual=[zero]*5 if terminal else unrefined
            beta=[]
            for center in self.centers:
                cc=ctx.mpf(center.numerator)/center.denominator
                if directed:beta.append(self.bump.beta(xx-cc))
                else:
                    u=(xx-cc)/c.mpf('.025')
                    v=c.mpf(0) if abs(u)>=1 else c.exp(-1/(1-u*u))/(c.mpf('.025')*self.normp)
                    v1=c.mpf(0) if not v else -2*u*v/(c.mpf('.025')*(1-u*u)**2)
                    beta.append((v,v1))
            f=sum((h[i+2]*beta[i][0] for i in range(3)),zero)
            fy=sum((h[i+2]*beta[i][1]*xx for i in range(3)),zero)
            g=sum((h[i]*beta[k][0] for i,k in enumerate((0,2))),zero)
            gy=sum((h[i]*beta[k][1]*xx for i,k in enumerate((0,2))),zero)
            E=am*(f+xx**ctx.mpf('.1'));Ey=am*(fy+xx**ctx.mpf('.1')*ctx.mpf('.1'))
            U=z*4+g;Uy=gy;M=z*4+residual[0]/xx
            Mz=J(ctx,[(i+1)*M[i+1] for i in range(2)])
            delta=ctx.mpf(self.reference.leading.core.delta_iv if directed else self.reference.leading.core.delta)
            radial=((z*U*2-z*M*(1-delta)-(1-z*z)*Mz)/(1-z*z*delta)).truncate(1)
            theta=residual[2]+xx**ctx.mpf('1.6')*(ctx.mpf(5)/8)
            mixed=z*theta*4+residual[1]
            energy=residual[3]-xx**ctx.mpf('1.2')*(ctx.mpf(5)/12)
            pressure=residual[4]+xx**ctx.mpf('.2')*(ctx.mpf(5)/2)
            raw=dict(E=E,V=U*invP,m=M*invP,h=am*theta*xx**ctx.mpf('-1.5'),
                k=am*mixed*xx**ctx.mpf('-1.5')*invP,
                e=(z*M*8-z*z*16)*(invP*invP)+(am*am)*energy/xx,p=(am*am)*pressure)
            p0=self._jet(data['inlet']['independent_P0_over_Pstar_squared'],directed).truncate(2)
            raw['pressure']=p0+raw['p']
            sides.append(dict(fields=raw,physical_similarity=dict(Uz=U,Ur_over_sqrt_R_over2=radial,
                Utheta_over_Pstar=E,Uz_y=Uy,Utheta_y_over_Pstar=Ey,Mz_over_R=M),
                residual=residual,unrefined=unrefined))
        p,t=sides;rows=self.reference._rows
        packet=dict(Z=Z,chart='original_leading_Rm_patch',native_x=x,source_family=self.family,
            selected_integer=data['inlet']['selected_integer'],
            controls_Z2=[rows(a,b) for a,b in zip(data['point'],data['truth'])],
            normalized_controls_Z2=[rows(a,b) for a,b in zip(data['normalized_point'],data['normalized_truth'])],
            implicit_certificate=encode(data['proof']),
            generic_leading_fields={n:rows(p['fields'][n],t['fields'][n]) for n in t['fields']},
            physical_similarity={n:rows(p['physical_similarity'][n],t['physical_similarity'][n]) for n in t['physical_similarity']},
            actual_partial_defects=[rows(a,b) for a,b in zip(p['residual'],t['residual'])],
            unrefined_full_or_partial_defects=[rows(a,b) for a,b in zip(p['unrefined'],t['unrefined'])],
            independent_P0_over_Pstar_squared=data['inlet']['independent_P0_over_Pstar_squared'],
            finite_N_incoming_history_errors=data['inlet']['finite_N_signed_unknown_history_errors'],
            positive_order_n1_own_history=data['inlet']['positive_order_n1_own_history'],
            contract=dict(map='original LEADING L*h+Q(h,invAm2)+unmodified_d=0',
                coefficient_order=['c1','c2','xi1','xi2','xi3'],normalization='diag(t,t,1e10*t²/Pstar²,1e10*t²/Pstar²,1e10*t²/Pstar²)',
                normalization_coordinates_are_not_field_parameters=True,point_root_is_proxy_with_full_directed_error=True,
                missing_Z3_Z5_not_padded=True,pressure_P0_unchanged_and_added_once=True,
                terminal_leading_five_function_identities=terminal,
                terminal_refinement_from_same_unique_implicit_map=terminal,
                physical_measure='dx/x; source derivative y=log(x)',
                finite_N_not_used_to_resolve_leading_controls=True,nonzero_finite_N_inlet_preserved=True),
            limits=dict(leading_control_and_field_Z_order=2,radial_recovery_Z_order=1,
                leading_original_branch_certified_for_query=True,new_whole_axis_numeric_certificate=False,
                finite_N_patch_corrections_transported=False,modified_patch_fields=False,
                global_physical_field=False,physical_gradients=False,global_N_admission=False,
                genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=packet
        return packet
