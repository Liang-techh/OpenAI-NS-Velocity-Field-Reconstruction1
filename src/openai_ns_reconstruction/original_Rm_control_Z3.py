"""Genuine Rm source/control Z3 from private, untruncated history transport.

The public Z2 modules remain unchanged. All new third rows are evaluated
from actual histories; the checked axial5 archive is independent evidence.
Finite-N perturbations remain separate and are not leading controls.
"""
import gzip
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_Rm_patch_numeric import OriginalRmPatchNumeric,implicit_Z2,encode


def implicit_Z3(iv,J,W,d,alpha,*,precision=135):
    """Enclose the genuine third coefficient using the same full Jacobian."""
    from lei_ren_part1_paper_interval_five_bump_inverse import identity_error,weighted_norm,dot,mag,intersection
    from lei_ren_part1_paper_compliant_actual_moment_patch import nonlinear_series
    if any(row.order!=3 for row in d) or alpha.order!=3:
        raise ValueError('Genuine source third-order jets required')
    with mp.workdps(precision):
        lower,proof=implicit_Z2(iv,J,W,[row.truncate(2) for row in d],alpha.truncate(2))
        jac=proof['second_row_same_Jacobian']
        midpoint=mp.matrix([[sum(_bounds(mp.mp,v))/2 for v in row] for row in W['L']]);ri=mp.inverse(midpoint)
        R=[[iv.mpf(ri[i,j]) for j in range(5)] for i in range(5)]
        K=identity_error(iv,R,jac);q=weighted_norm(iv,K,[iv.mpf(1)]*5)
        if _bounds(mp.mp,q)[1]>=1:raise ArithmeticError('Original same-Jacobian third inverse not contractive')
        known=[J(iv,list(row.coefficients)+[iv.mpf(0)]) for row in lower];Q=nonlinear_series(iv,W,known,alpha)
        rhs=[d[i][3]+Q[i][3] for i in range(5)];b=[-dot(row,rhs) for row in R]
        radius=max((mag(iv,v) for v in b),key=lambda v:_bounds(mp.mp,v)[1])/(1-q)
        hi=_bounds(mp.mp,radius)[1];third=[iv.mpf([-hi,hi]) for _ in range(5)]
        for _ in range(12):third=[intersection(iv,a,b[i]+dot(K[i],third)) for i,a in enumerate(third)]
        truth=[J(iv,list(row.coefficients)+[third[i]]) for i,row in enumerate(lower)]
    final=nonlinear_series(iv,W,truth,alpha)
    diagnostic=[d[i]+final[i]+sum((truth[j]*W['L'][i][j] for j in range(5)),J.constant(iv,0,3)) for i in range(5)]
    if any(not _bounds(mp.mp,v)[0]<=0<=_bounds(mp.mp,v)[1] for row in diagnostic for v in row.coefficients):
        raise ArithmeticError('Same original third-order implicit diagnostic excludes zero')
    return truth,dict(lower_Z2=proof,third_same_Jacobian=jac,third_preconditioner=R,third_inverse_contraction=q,
        third_known_rhs=rhs,third_preconditioned_rhs=b,third_containment_radius=radius,
        residual_only_diagnostic=diagnostic,third_existence_from_same_branch_and_inverse=True)


class OriginalRmControlZ3:
    def __init__(self, *, patch=None, **kwargs):
        if patch is not None and kwargs:raise ValueError('Supply existing patch or options')
        self.patch=patch if patch is not None else OriginalRmPatchNumeric(**kwargs)
        self.defect=self.patch.defect;self.reference=self.patch.reference
        for name in ('ctx','iv','J','family','directory'):
            setattr(self,name,getattr(self.patch,name))
        self.hashes,self.code_hashes=dict(self.patch.hashes),dict(self.patch.code_hashes)
        stem=PREFIX+'compliant_current_original_whole_Z_Rm_centered_inverse'
        receipt_path=self.directory/(stem+'_check.json');archive_path=self.directory/(stem+'.json.gz')
        raw=archive_path.read_bytes();receipt=json.loads(receipt_path.read_bytes());self.archive=json.loads(gzip.decompress(raw))
        if (not receipt['all_passed'] or receipt['source_family']!=self.family
                or self.archive['source_family']!=self.family
                or not receipt['actual_leading_controls_and_same_Jacobian_axial5_checked']
                or not receipt['same_common_unit_box_proves_unique_whole_Z_smooth_leading_family']):
            raise ValueError('Same original source-centered axial5 inverse receipt required')
        for name,digest in receipt['input_hashes'].items():
            if hashlib.sha256((self.directory/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed actual axial5 source: '+name)
            self.hashes[name]=digest
        for path in (receipt_path,archive_path):self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        if self.archive['input_hashes']!= {name:self.hashes[name] for name in self.archive['input_hashes']}:
            raise ValueError('Axial5 archive/source receipt mismatch')
        self._sources={};self._coefficients={};self._packets={}
        name='src/openai_ns_reconstruction/original_Rm_control_Z3.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def source(self,Z):
        if Z in self._sources:return self._sources[Z]
        ref=self.reference;long=ref.previous;c,iv,J=self.ctx,self.iv,self.J;o=ref.owner
        inp=long.inputs(Z);p,t=inp;proof=[];pp,tt={},{}
        for values in (p,t):
            if any(row.order<3 for row in values['six'].values()):
                raise ValueError('Actual internal R110 six histories through Z3 required')
        for rate,m,names in (('1.6',1,('H_over_Phi','Kpsi_over_Phi')),
                             ('1.2',2,('B_over_Phi2',)),('.2',2,('C_over_Phi2',))):
            radius=8*m*ref.leading._abs(t['B'][0])/long.T
            low,high=iv.mpf(rate)-radius,iv.mpf(rate)+radius
            if _bounds(c,low)[0]<=0:raise ArithmeticError('Original positive perturbed rate required')
            kernel=long.kernels.backward_kernel(iv,t['B'].truncate(3),long.T,long.T,m,low,high,proof)
            kp=J.constant(c,1/c.mpf(rate),3)
            for name in names:
                key=dict(H_over_Phi='H',Kpsi_over_Phi='Kpsi',B_over_Phi2='B',C_over_Phi2='C')[name]
                initial=t['six'][key]/(t['values']['Phi']**m)
                logs=t['B'].truncate(3)*m-long.T*iv.mpf(rate)
                memory=long.kernels.normalized_amplitude_decay(iv,logs,low,high,long.T,proof,initial=initial.truncate(3))
                sp=2*p['values']['Psi'].truncate(3) if name=='Kpsi_over_Phi' else J.constant(c,2 if name=='H_over_Phi' else 1,3)
                st=2*t['values']['Psi'].truncate(3) if name=='Kpsi_over_Phi' else J.constant(iv,2 if name=='H_over_Phi' else 1,3)
                pp[name]=sp*kp;tt[name]=memory+st*kernel
        for name,power in (('Mpsi',1),('Apsi',2)):
            targetp=(p['values']['Psi']**power).truncate(3);targett=(t['values']['Psi']**power).truncate(3)
            dp,dt=ref._decay(p['six'][name].truncate(3)-targetp,t['six'][name].truncate(3)-targett,
                1,long.Tp,long.T,proof,name+'_genuine_long_Z3_memory')
            pp[name]=targetp+dp;tt[name]=targett+dt
        long_six=dict(point=pp,truth=tt)
        psi=p['values']['Psi'].truncate(3),t['values']['Psi'].truncate(3)
        initial={};const=lambda x:(J.constant(c,c.mpf(x),3),J.constant(iv,iv.mpf(x),3))
        for name,target,rate in (('H_over_Phi',const('1.25'),'1.6'),
                ('Kpsi_over_Phi',(psi[0]*c.mpf('1.25'),psi[1]*iv.mpf('1.25')),'1.6'),
                ('B_over_Phi2',(J.constant(c,c.mpf(5)/6,3),J.constant(iv,iv.mpf(5)/6,3)),'1.2'),
                ('C_over_Phi2',const(5),'.2'),('Mpsi',psi,1),('Apsi',(psi[0]*psi[0],psi[1]*psi[1]),1)):
            initial[name]=ref._relax((pp[name],tt[name]),target,rate,ref.lengthp,ref.length,proof,name+'_reference_Z3')
        kernels=ref.kernels('1');e,ei=c.exp(-1),iv.exp(-1);six={}
        for name,rate in (('H_over_Phi','1.6'),('B_over_Phi2','1.2'),('C_over_Phi2','.2')):
            target=const('1.25') if name=='H_over_Phi' else (const(5) if name=='C_over_Phi2'
                else (J.constant(c,c.mpf(5)/6,3),J.constant(iv,iv.mpf(5)/6,3)))
            pair=ref._relax(initial[name],target,rate,c.mpf(2),iv.mpf(2),proof,name+'_restore_post_Z3')
            six[name]=pair
        for name,kind,rate,forcing in (('Mpsi','mean',1,psi),('Kpsi_over_Phi','mixed','1.6',(psi[0]*2,psi[1]*2)),
                                       ('Apsi','square',1,(psi[0]*psi[0],psi[1]*psi[1]))):
            kp,kt=kernels[kind]
            a,b=ref._decay(*initial[name],rate,c.mpf(1),iv.mpf(1),proof,name+'_restore_Z3')
            six[name]=ref._decay(a+forcing[0]*kp,b+forcing[1]*kt,rate,c.mpf(1),iv.mpf(1),proof,name+'_post_Z3')
        km,ks=kernels['mixed'],kernels['square']
        hp,ht=initial['H_over_Phi'];dp,dt=initial['Mpsi']
        baseline=dict(gamma1=((e+kernels['mean'][0])*e,(ei+kernels['mean'][1])*ei),
            gamma2=((e+ks[0])*e,(ei+ks[1])*ei),
            Kalpha_over_Phi=((hp*c.exp(-c.mpf('1.6'))+2*km[0])*c.exp(-c.mpf('1.6')),
                            (ht*iv.exp(-iv.mpf('1.6'))+2*km[1])*iv.exp(-iv.mpf('1.6'))),
            Dpsi=((dp*e+psi[0]*ks[0])*e,(dt*ei+psi[1]*ks[1])*ei))
        # Center angular errors before division by the extremely small
        # angular unit. A generic1e-1000 cap would be far too broad here.
        ap,at={},{}
        for short,name,rate in (('H','H_over_Phi','1.6'),('B','B_over_Phi2','1.2'),('C','C_over_Phi2','.2')):
            target=iv.mpf('1.25') if short=='H' else (iv.mpf(5)/6 if short=='B' else iv.mpf(5))
            delta=tt[name]-J.constant(iv,target,3);errors=[]
            for k in range(4):
                mag=ref.leading._abs(delta[k]);log=None if _bounds(c,mag)[1]==0 else iv.log(mag)-iv.mpf(rate)*(ref.length+2)
                errors.append(ref.leading._symmetric(self.defect._small(log,proof,'genuine_'+short+'_Z'+str(k))))
            ap[short]=J.constant(c,0,3);at[short]=J(iv,errors)
        sides=[]
        for directed,ctx in ((0,c),(1,iv)):
            z=J.variable(ctx,ctx.mpf(Z),3);j=ctx.mpf(o.amplitude.j_iv if directed else o.amplitude.j)
            eps=ctx.mpf(ref.leading.core.epsilon_iv if directed else ref.leading.core.epsilon)
            invP=ctx.exp(-ctx.mpf(ref.logP if directed else ref.logPp));am=(1+z*z).reciprocal()*ctx.exp(ctx.mpf('-.6'))
            ss={name:pair[directed] for name,pair in six.items()}
            bb={name:(J.constant(ctx,pair[directed],3) if name in ('gamma1','gamma2') else pair[directed])
                for name,pair in baseline.items()}
            angular=at if directed else ap
            mass=bb['gamma1']*j+ss['Mpsi']*eps
            mixed=(bb['Kalpha_over_Phi']*j+ss['Kpsi_over_Phi']*eps)/2
            square=bb['gamma2']*(j*j)+bb['Dpsi']*(2*j*eps)+ss['Apsi']*(eps*eps)
            defects=[mass,mixed,angular['H']/2,square*(invP*invP)/(am*am)-angular['B']/2,angular['C']/2]
            sides.append(dict(defects=defects,physical_mean=4*z+mass,six=ss,baseline=bb))
        out=dict(point=sides[0],truth=sides[1],long_six=long_six,proofs=proof)
        self._sources[Z]=out;return out

    def coefficients(self,Z):
        if Z in self._coefficients:return self._coefficients[Z]
        c,iv,J=self.ctx,self.iv,self.J;p=self.patch;source=self.source(Z)
        d=[row/p.units[i] for i,row in enumerate(source['truth']['defects'])]
        dp=[row/p.unitsp[i] for i,row in enumerate(source['point']['defects'])]
        z=J.variable(iv,iv.mpf(Z),3);zp=J.variable(c,c.mpf(Z),3)
        alpha=(1+z*z)**2*iv.exp(iv.mpf('1.2'))/(10**10)
        alphap=(1+zp*zp)**2*c.exp(c.mpf('1.2'))/(10**10)
        W={n:v if n in ('L','gg') else [q*p.angular for q in v] for n,v in p.W.items()}
        Wp=p._W(c.mpf(2),False,True);Wp={n:v if n in ('L','gg') else [q*p.angularp for q in v] for n,v in Wp.items()}
        from lei_ren_part1_paper_interval_five_bump_inverse import nonlinear,nonlinear_jacobian
        from lei_ren_part1_paper_compliant_actual_moment_patch import nonlinear_series
        truth,proof=implicit_Z3(iv,J,W,d,alpha,precision=p.precision+50)
        L=c.matrix(Wp['L']);guess=c.lu_solve(L,c.matrix([-row[0] for row in dp]))
        def residual(*h):return [sum((Wp['L'][i][j]*h[j] for j in range(5)),c.mpf(0))+nonlinear(c,Wp,h,alphap[0])[i]+dp[i][0] for i in range(5)]
        h0=list(c.findroot(residual,tuple(guess),tol=c.mpf(10)**(-p.precision+12),maxsteps=30))
        jq=nonlinear_jacobian(c,Wp,h0,alphap[0]);jp=c.matrix([[Wp['L'][i][j]+jq[i][j] for j in range(5)] for i in range(5)])
        point=[J(c,[v]) for v in h0]
        for n in range(1,4):
            known=[J(c,list(row.coefficients)+[c.mpf(0)]) for row in point];Q=nonlinear_series(c,Wp,known,J(c,list(alphap.coefficients[:n+1])))
            solved=c.lu_solve(jp,c.matrix([-dp[i][n]-Q[i][n] for i in range(5)]))
            point=[J(c,list(row.coefficients)+[solved[i]]) for i,row in enumerate(point)]
        out=dict(source=source,defects_point=dp,defects_truth=d,normalized_point=point,normalized_truth=truth,
            point=[row*p.unitsp[i] for i,row in enumerate(point)],truth=[row*p.units[i] for i,row in enumerate(truth)],
            proof=proof)
        self._coefficients[Z]=out;return out

    def evaluate(self,Z):
        if Z in self._packets:return self._packets[Z]
        data=self.coefficients(Z);rows=self.reference._rows
        out=dict(Z=Z,source_family=self.family,
            normalized_defects_Z3=[rows(a,b) for a,b in zip(data['defects_point'],data['defects_truth'])],
            normalized_controls_Z3=[rows(a,b) for a,b in zip(data['normalized_point'],data['normalized_truth'])],
            physical_controls_Z3=[rows(a,b) for a,b in zip(data['point'],data['truth'])],
            physical_mean_at_Rm_Z3=rows(data['source']['point']['physical_mean'],data['source']['truth']['physical_mean']),
            genuine_six_Rm_Z3={name:rows(data['source']['point']['six'][name],data['source']['truth']['six'][name]) for name in data['source']['point']['six']},
            implicit_certificate=encode(data['proof']),source_decay_proofs=self.reference.previous._serializable_proof(data['source']['proofs']),
            convention='Taylor coefficient D_Z^k/k!; all third rows evaluated before serialization',
            contract=dict(actual_R110_six_Z3_retained=True,finite_long_variable_source_kernel_Z3=True,
                reference_restore_postrestore_Z3=True,angular_centered_before_normalization=True,
                original_small_tail_not_signed_value=True,all_original_incoming_memories=True,
                genuine_third_defects_not_zero_padded=True,original_normalization_units_unchanged=True,
                same_original_nonlinear_map_and_Jacobian=True,third_directed_inverse_containment=True,
                point_root_and_recurrence_only_proxies=True,finite_N_not_used_to_resolve_leading_controls=True),
            limits=dict(original_leading_source_and_controls_Z3_with_absolute_errors=True,
                Rm_patch_pressure_rows_Z2=False,finite_N_patch_C3_transported=False,
                Rm_Rh_signed_m15_integrated=False,new_whole_axis_Z3_certificate=False,
                source_graph_C6_numeric_certificate=False,global_five_moment_numeric_join=False,
                common_global_N=False,exact_heat_exterior=False,stress_cone=False,
                genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[Z]=out;return out
