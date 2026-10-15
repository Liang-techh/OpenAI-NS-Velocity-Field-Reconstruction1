"""Finite-width original n1 transport through the admitted first half collar.

The inverse-width forcing is extracted before applying dx=hb*x*ds/2.
Incoming values are the complete-error core Volterra values. The remaining
ODE integral is bounded with the admitted common-inner full solution norm.
"""
import gzip
import hashlib
import importlib
import json
import math

import mpmath as mp

from .original_core_numeric import PREFIX, ROOT, _box, _bounds
from .original_first_collar_jets import OriginalFirstCollarJets
from .original_n1_core_numeric import linear_rhs
from .original_n1_core_volterra_numeric import OriginalN1CoreVolterraNumeric


class Laurent:
    """Finite original hb powers, with native axial Taylor-jet coefficients."""
    def __init__(self, rows):
        self.rows = dict(rows)
        self.zero = next(iter(rows.values()))*0
    def cast(self, value):
        return value if isinstance(value, Laurent) else Laurent({0:self.zero+value})
    def __add__(self, other):
        other = self.cast(other)
        return Laurent({p:self.rows.get(p,self.zero)+other.rows.get(p,self.zero) for p in self.rows.keys()|other.rows.keys()})
    __radd__ = __add__
    def __neg__(self):
        return Laurent({p:-v for p,v in self.rows.items()})
    def __sub__(self, other):
        return self+-self.cast(other)
    def __rsub__(self, other):
        return self.cast(other)+-self
    def __mul__(self, other):
        other=self.cast(other);rows={}
        for p,a in self.rows.items():
            for q,b in other.rows.items():
                rows[p+q]=rows.get(p+q,self.zero)+a*b
        return Laurent(rows)
    __rmul__=__mul__
    def __truediv__(self, other):
        if isinstance(other,Laurent):
            if set(other.rows)!={0}:raise ValueError('Only width-free denominators may be inverted')
            other=other.rows[0]
        return Laurent({p:v/other for p,v in self.rows.items()})
    def dz(self):
        return Laurent({p:type(v)(v.ctx,[(k+1)*v[k+1] for k in range(v.order)]) for p,v in self.rows.items()})
    def cut(self, order):
        return Laurent({p:v.truncate(order) for p,v in self.rows.items()})


def axial_second(rows,z,rho,de,weight):
    f,fr,frr=rows;d,L=1-z*z,1-z*z*de
    fz,frz=f.dz(),fr.dz()
    A=f*(z*weight)+fz*d-fr*(z*(2*rho))
    Az=f*weight+fz*(z*(weight-2))+fz.dz()*d-fr*(2*rho)-frz*(z*(2*rho))
    Ar=fr*(z*(weight-2))+frz*d-frr*(z*(2*rho))
    return (A*(z*(weight-1+de))/L+(Az/L+A*(z*(2*de))/(L*L))*d-Ar*(z*(2*rho))/L)/L


class WidthChange:
    """a(h)=a(0)+h*e; a jet and a uniform norm bound for e."""
    def __init__(self, owner, base, error=0):
        self.o,self.base,self.error=owner,base,owner.iv.mpf(error)
    def cast(self,v):
        return v if isinstance(v,WidthChange) else WidthChange(self.o,self.base*0+v)
    def __add__(self,v):
        v=self.cast(v);return WidthChange(self.o,self.base+v.base,self.error+v.error)
    __radd__=__add__
    def __neg__(self):return WidthChange(self.o,-self.base,self.error)
    def __sub__(self,v):return self+-self.cast(v)
    def __rsub__(self,v):return self.cast(v)+-self
    def __mul__(self,v):
        v=self.cast(v);o=self.o
        return WidthChange(o,self.base*v.base,o._norm(self.base)*v.error+self.error*o._norm(v.base)+o._width_error(self.error*v.error))
    __rmul__=__mul__
    def inverse(self):
        o=self.o;b=self.base.reciprocal();n=o._norm(b);r=o._width_error(n*self.error)
        if _bounds(o.ctx,r)[1]>=o.ctx.mpf('.5'):raise ArithmeticError('Original collar denominator variation is not admitted')
        return WidthChange(o,b,n*n*self.error/(1-r))
    def __truediv__(self,v):return self*self.cast(v).inverse()


class OriginalFirstCollarTransport:
    fields=('f','u','k','p','fx','ux')
    def __init__(self,*,precision=80,directory=None):
        self.leading=OriginalFirstCollarJets(precision=precision,directory=directory)
        self.n1=OriginalN1CoreVolterraNumeric(precision=precision,directory=directory)
        self.ctx,self.iv,self.Jet=self.leading.ctx,self.leading.iv,self.leading.Jet
        self.family,self.precision=self.leading.family,precision
        if self.n1.family!=self.family:raise ValueError('Common original leading/n1 family required')
        self.hashes={**self.leading.hashes,**self.n1.hashes}
        self.code_hashes={**self.leading.code_hashes,**self.n1.code_hashes}
        for relative in ('src/openai_ns_reconstruction/original_n1_core_volterra_numeric.py',
                         'src/openai_ns_reconstruction/original_first_collar_transport.py'):
            self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        name=PREFIX+'compliant_current_original_n1_inner_analytic_domain.json.gz'
        if self.hashes.get(name)!=hashlib.sha256((self.leading.core.directory/name).read_bytes()).hexdigest():
            raise ValueError('Admitted actual first-inner analytic bounds must be pinned by all-word receipt')
        self.domain=json.loads(gzip.decompress((self.leading.core.directory/name).read_bytes()))
        if self.domain['source_family']!=self.family:raise ValueError('Actual inner analytic family changed')
        self._core={}
        self.sigma_jets=importlib.import_module(PREFIX+'compliant_flat_pulse_derivatives').sigma_jets

    def _core_comparison_changes(self, name, i, s, order):
        o,iv=self.leading,self.iv
        d=self.domain['joint_core_domain']
        bound=_box(iv,d['Xh_norm_majorants'][name])*iv.mpf(d['core_value_factor'])
        rr,rz=iv.mpf(d['Cauchy_rho_radius']),_box(iv,d['Cauchy_Z_radius'])
        change=4*iv.mpf(s)*iv.exp(o.hcap*iv.mpf(s))*bound*math.factorial(i+1)/rr**(i+1)
        return sum((change*o.weight**k/rz**k for k in range(order+1)),iv.mpf(0))

    def _inputs(self,Z,s,order=4):
        o,c,iv,J=self.leading,self.ctx,self.iv,self.Jet
        inlet=o.inlet_jets(Z,degree=24,order=order)
        if Z not in o._prepared:
            with mp.workdps(280):o._prepared[Z]=o.bridge.prepare(Z)
        _,Bell,_,Bpsi=o._source_jet_increments(o._prepared[Z],s,order)
        a=inlet['truth']['values']
        Bphi=o._norm(a['Phi'])*Bell*iv.exp(o._width_error(Bell))
        Bmean=iv.mpf(s)*o._norm(a['Psi']-a['Mpsi'])+o._width_error(iv.mpf(s)*Bpsi)
        delta=dict(Phi=Bphi,Psi=Bpsi,Mpsi=Bmean)
        actual={n:a[n]+o._norm_error_jet(o._width_error(delta[n]),order) for n in delta}
        comparison={}
        for name in ('Phi','Psi'):
            comparison[name]=[inlet['truth']['radial'][name][i]+o._norm_error_jet(
                o._width_error(self._core_comparison_changes(name,i,s,order)),order) for i in range(3)]
        with mp.workdps(280):sigma=self.sigma_jets(iv,iv.mpf([0,iv.mpf(s).b]))
        # Whole-prefix controls; hb is enclosed only after its source log proof.
        hv=iv.mpf([0,_bounds(c,o._width_error(iv.mpf(1)))[1]])
        sig,sp=sigma[0],sigma[1]
        rho=J.constant(iv,iv.mpf([4,_bounds(c,iv.mpf(4)+o._width_error(4*iv.mpf(s)*iv.exp(o.hcap*iv.mpf(s))))[1]]),order)
        x=J.constant(iv,iv.sqrt(rho[0]),order)
        chi=J.constant(iv,1-sig+hv*sig,order)
        chir=rho.reciprocal()*(-(1-hv)*sp)
        z=J.variable(iv,iv.mpf(Z),order);de,e=o.core.delta_iv,o.core.epsilon_iv
        U0=4*z+o.core.j_iv;d,L=1-z*z,1-z*z*de
        A=comparison['Phi'][1]/comparison['Phi'][0]
        Ar=comparison['Phi'][2]/comparison['Phi'][0]-A*A
        ratio=actual['Phi']/comparison['Phi'][0]
        phi=[Laurent({0:actual['Phi']}),Laurent({0:actual['Phi']*chi*A}),
             Laurent({0:actual['Phi']*(chi*chi*A*A+chi*Ar),-1:actual['Phi']*chir*A})]
        psir=ratio*chi*comparison['Psi'][1]
        psi=[Laurent({0:actual['Psi']}),Laurent({0:psir}),
             Laurent({0:ratio*(chi*(chi-1)*A*comparison['Psi'][1]+chi*comparison['Psi'][2]),
                      -1:ratio*chir*comparison['Psi'][1]})]
        mr=(actual['Psi']-actual['Mpsi'])/rho
        mean=[Laurent({0:actual['Mpsi']}),Laurent({0:mr}),Laurent({0:(psir-2*mr)/rho})]
        U=[psi[0]*e+U0,psi[1]*e,psi[2]*e]
        M=[mean[0]*e+U0,mean[1]*e,mean[2]*e]
        Q=[(w*(2*z)-m*(z*(1-de))-m.dz()*d)/L for w,m in zip(U,M)]
        fixed,_,_,_=o.core._seed(Z,24,order,directed=True)
        g=J(iv,[-v for v in fixed['ell_Z_taylor'][:order+1]])
        S=J(iv,fixed['S_Z_taylor'][:order+1])
        return dict(inlet=inlet,actual=actual,comparison=comparison,delta=delta,
                    Phi=phi,Psi=psi,U=U,M=M,Q=Q,z=z,rho=rho,x=x,de=de,e=e,
                    P2=iv.exp(2*o.core.lp_iv),g=g,S=S)

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
        return {name:value.cut(1) for name,value in out.items()}

    def _regular_operator_bound(self,v):
        o,iv,J=self.leading,self.iv,self.Jet
        cut=lambda q:q.truncate(1)
        z,r,x=map(cut,(v['z'],v['rho'],v['x']))
        de,e,P2=v['de'],v['e'],v['P2'];d,L=1-z*z,1-z*z*de
        bound=9*_box(iv,self.n1.tail['regular_solution_bound'])
        gap=_box(iv,self.n1.tail['derivative_Cauchy_gap'])/2
        state={name:J(iv,[o._symmetric(bound/gap**k) for k in range(3)]) for name in self.fields}
        V={n:cut(q) for n,q in state.items()}
        V.update({n+'Z':cut(o._derivative(state[n])) for n in ('f','u','k','p')})
        phi,pr=[cut(q.rows[0]) for q in v['Phi'][:2]]
        U,Ur=[cut(q.rows[0]) for q in v['U'][:2]]
        Uz=cut(o._derivative(v['U'][0].rows[0]));phiz=cut(o._derivative(v['Phi'][0].rows[0]))
        Q=cut(v['Q'][0].rows[0]);EU=Ur*r
        V.update(z=z,rho=r,delta=J.constant(iv,de,1),epsilon=J.constant(iv,e,1),P2=J.constant(iv,P2,1),
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
            ctx=iv if directed else c;z=J.variable(ctx,ctx.mpf(v['Z']),1)
            de=o.core.delta_iv if directed else o.core.delta;P2=iv.exp(2*o.core.lp_iv) if directed else c.exp(2*o.core.logP)
            L=1-z*z*de;radial=v['inlet'][label]['radial']
            ph=radial['Phi'][1].truncate(1);ps=radial['Psi'][1].truncate(1)/P2
            rows=dict(p=z*ps*(-8)/L,fx=z*z*ph*32/(L*L),ux=z*z*ps*64/(L*L))
            (truth if directed else point).update(rows)
        # Finite-radius/ratio variations of the complete singular measure.
        C=lambda b,err=0:WidthChange(o,b,err)
        z=C(J.variable(iv,iv.mpf(v['Z']),1));de=v['de'];L=1-z*z*de
        rho=C(J.constant(iv,4,1),4*iv.mpf(s)*iv.exp(o.hcap*iv.mpf(s)))
        x=C(J.constant(iv,2,1),iv.mpf(s)*iv.exp(o.hcap*iv.mpf(s)/2))
        one=C(J.constant(iv,1,1),iv.mpf(1))
        a=v['inlet']['truth'];phi=C(a['values']['Phi'].truncate(1),v['delta']['Phi'])
        comp=C(a['radial']['Phi'][0].truncate(1),self._core_comparison_changes('Phi',0,s,1))
        ratio=phi/comp
        ph=C(a['radial']['Phi'][1].truncate(1),self._core_comparison_changes('Phi',1,s,1))
        ps=C(a['radial']['Psi'][1].truncate(1)/v['P2'],self._core_comparison_changes('Psi',1,s,1)/v['P2'])
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
        v=self._inputs(Z,s);v['Z']=Z
        forcing=self._forcing(v);linear=self._regular_operator_bound(v)
        jp,jt,variation,sigma,sigbox=self._jump(v,s)
        if Z not in self._core:self._core[Z]=self.n1.evaluate(Z,'4',degree=18,kernel_degree=24)
        core=self._core[Z]
        zero=self.Jet.constant(iv,0,1)
        scales=dict(f=iv.mpf(1),u=iv.mpf(1),k=iv.mpf(1),p=iv.mpf(1),fx=v['e']**2,ux=v['e']*v['P2'])
        points,truths,jumps,proofs={},{},{},{}
        for name in self.fields:
            source=forcing[name].rows.get(0,zero) if name in forcing else zero
            regular=o._width_error(iv.mpf(s)*iv.mpf('1.5')*(o._norm(linear[name])+o._norm(source))/scales[name])
            change=o._width_error(o._abs(sigbox)*variation[name]) if name in variation else iv.mpf(0)
            error=o._norm_error_jet(regular+change,1)
            p=jp[name]*sigma if name in jp else self.Jet.constant(c,0,1)
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
        z,zi=c.mpf(Z),iv.mpf(Z)
        def recover(rows,zz,de):
            mean=rows['u'][0]+rows['k'][0];meanZ=rows['u'][1]+rows['k'][1]
            return mean,(2*zz*rows['u'][0]-(1+de)*zz*mean-(1-zz*zz)*meanZ)/(1-de*zz*zz)
        pm,pq=recover(points,z,o.core.delta);tm,tq=recover(truths,zi,o.core.delta_iv)
        out={name:o._jet_summary(points[name],truths[name]) for name in self.fields}
        out.update(mean_u=o._summary(pm,tm),q=o._summary(pq,tq))
        return dict(Z=Z,phase=s,source_family=self.family,normalized_original_n1=out,
            incoming_core=core['normalized_original_n1'],factored_first_collar_increments=jumps,
            finite_width_remainder_budgets=proofs,
            inverse_width_forcing={name:{str(p):[o._summary(c.mpf(0),q[k]) for k in range(2)]
                for p,q in row.rows.items()} for name,row in forcing.items()},
            exact_geometry=dict(rho='4exp(hb*s)',R='epsilon*rho',dx='hb*sqrt(rho)*ds/2'),
            sigma_approximation=c.nstr(sigma,self.precision),
            source_width_cancelled_before_evaluation=True,original_nonzero_core_history_retained=True,
            full_original_finite_width_n1_value_Z1_errors=True,
            solution_norm_source='accepted original regular series bound on the same common-inner natural frame',
            limits=dict(admitted_first_phase_max='.5',exact_original_point_oracle=False,
                remaining_first_second_lower_transport=False,global_moment_repair=False,
                physical_collar_velocity=False,temporal_n_ge_2_recursion=False,full_corrected_NS=False))
