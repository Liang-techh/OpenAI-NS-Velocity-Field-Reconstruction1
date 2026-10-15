"""Full original n1 first-half collar axial Taylor order two and q_Z.

Uses genuine Z2 inlet data and Z5 leading source rows. The local regular
remainder includes W_ZZZ by a Cauchy bound, without inventing a numerical
third derivative. Original inverse-width cancellation precedes evaluation.
"""
import hashlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _box, _bounds
from .original_first_collar_transport import OriginalFirstCollarTransport, WidthChange, axial_second
from .original_n1_core_numeric import linear_rhs
from .original_n1_core_axial2 import OriginalN1CoreAxial2


class OriginalFirstCollarAxial2(OriginalFirstCollarTransport):
    def __init__(self,*,precision=85,directory=None):
        super().__init__(precision=precision,directory=directory)
        self.n1=OriginalN1CoreAxial2(precision=precision,directory=directory)
        if self.n1.family!=self.family:
            raise ValueError('Same original leading and axial2 source family required')
        for target,extra in ((self.hashes,self.n1.hashes),(self.code_hashes,self.n1.code_hashes)):
            for name,digest in extra.items():
                if name in target and target[name]!=digest:
                    raise ValueError('Conflicting axial2 prerequisite: '+name)
                target[name]=digest
        for stem in ('compliant_current_original_n1_core_axial2',
                     'compliant_current_original_first_collar_transport'):
            path=self.leading.core.directory/(PREFIX+stem+'.json')
            report=json.loads(path.read_text(encoding='utf-8'))
            if not report['checks']['passed'] or report['source_family']!=self.family:
                raise ValueError('Accepted same-family Z2 inlet and collar prerequisites required')
            for name,digest in report['implementation_hashes'].items():
                if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed accepted collar prerequisite: '+name)
                self.code_hashes[name]=digest
            self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        relative='src/openai_ns_reconstruction/original_first_collar_axial2.py'
        self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        self._core={}

    def _forcing(self,v):
        z,r,x,de,e,P2=v['z'],v['rho'],v['x'],v['de'],v['e'],v['P2']
        d,L=1-z*z,1-z*z*de
        phi,U,Q=v['Phi'],v['U'],v['Q']
        grad=v['g'];m=d*grad/L;a=-2-de
        f,fr=phi[:2]
        cross=(f*(z*(a-3))+f.dz()*(2*d)-fr*(z*(4*r)))*m/L+f*(d*self.leading._derivative(m))/L
        gf=f*(-2*m*m)+cross*(2*e)-axial_second(phi,z,r,de,a)*(2*e*e)
        gu=-axial_second(U,z,r,de,-1-de)*2
        q,qr,qrr=Q;qrho=qr*r
        flow=(q+q.dz()*(z*((1-de)/2))+qrho)/L+q*(q/2+qrho)
        flow+=U[0]*(q*(z*(-2))+q.dz()*d-qrho*(z*2))/L
        Np=-flow/(2*P2)+(qr*2+qrr*r)/(e*P2)
        gu-=Np*(z*(4*e*P2)*r/L)
        out=dict(p=Np*(2*x),fx=gf,ux=gu)
        if any(min(q.rows)<-1 for q in out.values()):raise ArithmeticError('Unexpected inverse-width power in original first forcing')
        return {name:value.cut(2) for name,value in out.items()}

    def _regular_operator_bound(self,v):
        o,iv,J=self.leading,self.iv,self.Jet
        cut=lambda q:q.truncate(2)
        z,r,x=map(cut,(v['z'],v['rho'],v['x']))
        de,e,P2=v['de'],v['e'],v['P2'];d,L=1-z*z,1-z*z*de
        bound=9*_box(iv,self.n1.tail['regular_solution_bound'])
        gap=_box(iv,self.n1.tail['derivative_Cauchy_gap'])/2
        state={name:J(iv,[o._symmetric(bound/gap**k) for k in range(4)]) for name in self.fields}
        V={n:cut(q) for n,q in state.items()}
        V.update({n+'Z':cut(o._derivative(state[n])) for n in ('f','u','k','p')})
        phi,pr=[cut(q.rows[0]) for q in v['Phi'][:2]]
        U,Ur=[cut(q.rows[0]) for q in v['U'][:2]]
        Uz=cut(o._derivative(v['U'][0].rows[0]));phiz=cut(o._derivative(v['Phi'][0].rows[0]))
        Q=cut(v['Q'][0].rows[0]);EU=Ur*r
        V.update(z=z,rho=r,delta=J.constant(iv,de,2),epsilon=J.constant(iv,e,2),P2=J.constant(iv,P2,2),
            F2=cut(v['S'])/(e*e),L=L,d=d,A=phi,C=phi+pr*r,g=cut(v['g']),
            E=z*U-iv.mpf('.5'),H=z*((1-de)/2)+d*U,B=Q+(1-z*U*2)/L,Q=Q,
            Zf_without_gauge=(z*phi*(-2-de)+d*phiz-z*pr*(2*r))/L,
            Zw=(z*U*(-1-de)+d*Uz-z*EU*2)/L,Euler_U=EU,
            Euler_f=V['fx']*x/2,Euler_u=V['ux']*x/2)
        rf,ru,rp=linear_rhs(V)
        return dict(f=V['fx'],u=V['ux'],k=-V['ux']-V['k']*2/x,
                    p=rp*(2*x),fx=rf-V['fx']*3/x,ux=ru-V['ux']/x)

    def _jump(self,v,s):
        o,c,iv,J=self.leading,self.ctx,self.iv,self.Jet
        point,truth={},{}
        for directed,label in ((False,'point'),(True,'truth')):
            ctx=iv if directed else c;z=J.variable(ctx,ctx.mpf(v['Z']),2)
            de=o.core.delta_iv if directed else o.core.delta;P2=iv.exp(2*o.core.lp_iv) if directed else c.exp(2*o.core.logP)
            L=1-z*z*de;radial=v['inlet'][label]['radial']
            ph=radial['Phi'][1].truncate(2);ps=radial['Psi'][1].truncate(2)/P2
            rows=dict(p=z*ps*(-8)/L,fx=z*z*ph*32/(L*L),ux=z*z*ps*64/(L*L))
            (truth if directed else point).update(rows)
        # Finite-radius/ratio variations of the complete singular measure.
        C=lambda b,err=0:WidthChange(o,b,err)
        z=C(J.variable(iv,iv.mpf(v['Z']),2));de=v['de'];L=1-z*z*de
        rho=C(J.constant(iv,4,2),4*iv.mpf(s)*iv.exp(o.hcap*iv.mpf(s)))
        x=C(J.constant(iv,2,2),iv.mpf(s)*iv.exp(o.hcap*iv.mpf(s)/2))
        one=C(J.constant(iv,1,2),iv.mpf(1))
        a=v['inlet']['truth'];phi=C(a['values']['Phi'].truncate(2),v['delta']['Phi'])
        comp=C(a['radial']['Phi'][0].truncate(2),self._core_comparison_changes('Phi',0,s,2))
        ratio=phi/comp
        ph=C(a['radial']['Phi'][1].truncate(2),self._core_comparison_changes('Phi',1,s,2))
        ps=C(a['radial']['Psi'][1].truncate(2)/v['P2'],self._core_comparison_changes('Psi',1,s,2)/v['P2'])
        varying=dict(p=one*rho*z*ratio*ps*(-2)/L,
                     fx=one*rho*x*z*z*ratio*ph*4/(L*L),
                     ux=one*rho*x*z*z*ratio*ps*8/(L*L))
        sp=c.mpf(s)
        sigma=c.mpf(0) if not sp else 1/(1+c.exp(1/sp**2-1/(1-sp)**2))
        with mp.workdps(280):sigbox=self.sigma_jets(iv,iv.mpf(s))[0]
        return point,truth,{name:q.error for name,q in varying.items()},sigma,sigbox

    def evaluate(self,Z,s):
        c,iv,o=self.ctx,self.iv,self.leading
        if type(s) is not str or not c.isfinite(c.mpf(s)) or not 0<=c.mpf(s)<=c.mpf('.5'):
            raise ValueError('Admitted original first-inner phase is [0,1/2]')
        v=self._inputs(Z,s,order=5);v['Z']=Z
        forcing=self._forcing(v);linear=self._regular_operator_bound(v)
        jp,jt,variation,sigma,sigbox=self._jump(v,s)
        if Z not in self._core:self._core[Z]=self.n1.evaluate(Z,'4',degree=20,kernel_degree=28)
        core=self._core[Z]
        zero=self.Jet.constant(iv,0,2)
        scales=dict(f=iv.mpf(1),u=iv.mpf(1),k=iv.mpf(1),p=iv.mpf(1),fx=v['e']**2,ux=v['e']*v['P2'])
        points,truths,jumps,proofs={},{},{},{}
        for name in self.fields:
            source=forcing[name].rows.get(0,zero) if name in forcing else zero
            regular=o._width_error(iv.mpf(s)*iv.mpf('1.5')*(o._norm(linear[name])+o._norm(source))/scales[name])
            change=o._width_error(o._abs(sigbox)*variation[name]) if name in variation else iv.mpf(0)
            error=o._norm_error_jet(regular+change,2)
            p=jp[name]*sigma if name in jp else self.Jet.constant(c,0,2)
            t=jt[name]*sigbox+error if name in jt else error
            raw=core['normalized_original_n1'][name]
            cp=self.Jet(c,[c.mpf(q['approximation']) for q in raw])
            ct=self.Jet(iv,[iv.mpf([mp.make_mpf(t) for t in q['original_source_enclosure_exact_mpf_tuples']]) for q in raw])
            scale_point=(o.core.epsilon**2 if name=='fx' else o.core.epsilon*c.exp(2*o.core.logP) if name=='ux' else c.mpf(1))
            points[name],truths[name]=cp+p*scale_point,ct+t*scales[name]
            jumps[name]=dict(epsilon_power=2 if name=='fx' else 1 if name=='ux' else 0,
                Pstar_squared_power=1 if name=='ux' else 0,
                coefficient=o._jet_summary(p,t),finite_width_error_is_included=True)
            proofs[name]=dict(regular_scaled_error_upper=o._summary(c.mpf(0),o._symmetric(regular)),
                             singular_variation_scaled_error_upper=o._summary(c.mpf(0),o._symmetric(change)))
        def recover(rows,ctx,de):
            z=self.Jet.variable(ctx,ctx.mpf(Z),2)
            mean=rows['u']+rows['k']
            q=(2*z*rows['u']-z*(1+de)*mean-(1-z*z)*o._derivative(mean))/(1-z*z*de)
            return mean,q.truncate(1)
        pm,pq=recover(points,c,o.core.delta);tm,tq=recover(truths,iv,o.core.delta_iv)
        out={name:o._jet_summary(points[name],truths[name]) for name in self.fields}
        out.update(mean_u=o._jet_summary(pm,tm),q=o._jet_summary(pq,tq))
        def radial(rows,ctx,de,rho,x):
            z=self.Jet.variable(ctx,ctx.mpf(Z),2);L=1-z*z*de;d=1-z*z
            return dict(f_rho=rows['fx']/(2*x),u_rho=rows['ux']/(2*x),
                mean_u_rho=-rows['k']/rho,
                q_rho=((z*rows['ux']/x+z*(1+de)*rows['k']/rho
                    +d*o._derivative(rows['k'])/rho)/L).truncate(1))
        rp=radial(points,c,o.core.delta,c.mpf(4),c.mpf(2))
        rt=radial(truths,iv,o.core.delta_iv,v['rho'].truncate(2),v['x'].truncate(2))
        radial_out={name:o._jet_summary(rp[name],rt[name]) for name in rp}
        return dict(Z=Z,phase=s,source_family=self.family,normalized_original_n1=out,
            native_first_radial_derivatives=radial_out,
            incoming_core=core['normalized_original_n1'],factored_first_collar_increments=jumps,
            finite_width_remainder_budgets=proofs,
            inverse_width_forcing={name:{str(p):[o._summary(c.mpf(0),q[k]) for k in range(3)]
                for p,q in row.rows.items()} for name,row in forcing.items()},
            exact_geometry=dict(rho='4exp(hb*s)',R='epsilon*rho',dx='hb*sqrt(rho)*ds/2'),
            sigma_approximation=c.nstr(sigma,self.precision),
            source_width_cancelled_before_evaluation=True,original_nonzero_core_history_retained=True,
            full_original_finite_width_n1_value_Z2_errors=True,
            axial_coefficient_convention="D_Z^k/k!; multiply order2 values and errors by 2 for ordinary derivatives",
            solution_norm_source='accepted original regular series bound on the same common-inner natural frame',
            limits=dict(admitted_first_phase_max='.5',exact_original_point_oracle=False,
                q_ZZ_available=False,complete_Cartesian_C2_join=False,
                remaining_first_second_lower_transport=False,global_moment_repair=False,
                physical_collar_velocity=False,temporal_n_ge_2_recursion=False,full_corrected_NS=False))
