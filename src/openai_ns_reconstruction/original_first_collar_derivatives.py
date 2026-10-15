"""Actual radial derivative values, retaining inverse-width collar sectors.

No finite differences or inverse hb exponential are used. Second n1 radial
derivatives come from the original full local equation and its enclosed
solution, not from differentiating a nominal radial Taylor polynomial.
"""
import hashlib
import json

import mpmath as mp

from .original_core_numeric import PREFIX, ROOT, _bounds
from .original_first_collar_transport import OriginalFirstCollarTransport, Laurent
from .original_n1_core_numeric import linear_rhs


class OriginalFirstCollarDerivatives:
    def __init__(self,*,precision=85,directory=None):
        self.transport=OriginalFirstCollarTransport(precision=precision,directory=directory)
        self.leading=self.transport.leading
        self.ctx,self.iv,self.Jet=self.transport.ctx,self.transport.iv,self.transport.Jet
        self.family,self.precision=self.transport.family,precision
        self.hashes,self.code_hashes=dict(self.transport.hashes),dict(self.transport.code_hashes)
        name=PREFIX+'compliant_current_original_first_collar_transport.json'
        path=self.leading.core.directory/name;receipt=json.loads(path.read_text(encoding='utf-8'))
        if not receipt['checks']['passed'] or receipt['source_family']!=self.family:
            raise ValueError('Accepted same-family original finite-width n1 transport required')
        for name,digest in receipt['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed original transport implementation: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        relative='src/openai_ns_reconstruction/original_first_collar_derivatives.py'
        self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        self._cache={}

    def _restore(self,row):
        return self.iv.mpf([mp.make_mpf(tuple(t)) for t in row['source_enclosure_exact_mpf_tuples']])

    def _point_sources(self,Z,s,native,directed):
        o,J=self.leading,self.Jet;c=self.iv if directed else self.ctx
        label='truth' if directed else 'point';inlet=native['inlet'][label]
        a=native['actual'] if directed else inlet['values']
        comparison=native['comparison'] if directed else inlet['radial']
        order=4;z=J.variable(c,c.mpf(Z),order)
        de=o.core.delta_iv if directed else o.core.delta
        e=o.core.epsilon_iv if directed else o.core.epsilon
        r=native['rho'] if directed else J.constant(c,4,order)
        x=J(c,[c.sqrt(r[0])]+[c.mpf(0)]*order)
        if directed:
            with mp.workdps(280):sigma=self.transport.sigma_jets(c,c.mpf(s))
            sig,sp=sigma[0],sigma[1]
            h=c.mpf([0,_bounds(self.ctx,o._width_error(c.mpf(1)))[1]])
            g,S=native['g'],native['S']
        else:
            t=c.mpf(s)
            sig=c.mpf(0) if not t else 1/(1+c.exp(1/t**2-1/(1-t)**2))
            sp=c.mpf(0) if not t else sig*(1-sig)*(2/t**3+2/(1-t)**3)
            h=c.mpf(0)  # Numerical proxy; the directed branch retains original hb.
            fixed,_,_,_=o.core._seed(Z,24,order,directed=False)
            g=J(c,[-v for v in fixed['ell_Z_taylor'][:order+1]])
            S=J(c,fixed['S_Z_taylor'][:order+1])
        chi=J.constant(c,1-sig+h*sig,order)
        slope=r.reciprocal()*sp
        A=comparison['Phi'][1]/comparison['Phi'][0]
        Ar=comparison['Phi'][2]/comparison['Phi'][0]-A*A
        ratio=a['Phi']/comparison['Phi'][0]
        phi=[Laurent({0:a['Phi']}),Laurent({0:a['Phi']*chi*A}),
             Laurent({0:a['Phi']*(chi*chi*A*A+chi*Ar+slope*A),
                      -1:-a['Phi']*slope*A})]
        psir=ratio*chi*comparison['Psi'][1]
        psi=[Laurent({0:a['Psi']}),Laurent({0:psir}),
             Laurent({0:ratio*(chi*(chi-1)*A*comparison['Psi'][1]
                             +chi*comparison['Psi'][2]+slope*comparison['Psi'][1]),
                      -1:-ratio*slope*comparison['Psi'][1]})]
        mr=(a['Psi']-a['Mpsi'])/r
        mean=[Laurent({0:a['Mpsi']}),Laurent({0:mr}),Laurent({0:(psir-2*mr)/r})]
        U0=4*z+(o.core.j_iv if directed else o.core.j)
        U=[psi[0]*e+U0,psi[1]*e,psi[2]*e]
        M=[mean[0]*e+U0,mean[1]*e,mean[2]*e]
        d,L=1-z*z,1-z*z*de
        Q=[(u*(2*z)-m*(z*(1-de))-m.dz()*d)/L for u,m in zip(U,M)]
        return dict(Phi=phi,Psi=psi,Mpsi=mean,U=U,M=M,Q=Q,z=z,rho=r,x=x,
                    de=de,e=e,P2=c.exp(2*(o.core.lp_iv if directed else o.core.logP)),g=g,S=S)

    def _state(self,packet,directed):
        c=self.iv if directed else self.ctx
        return {name:self.Jet(c,[self._restore(row) if directed else c.mpf(row['approximation'])
            for row in packet['normalized_original_n1'][name]]) for name in self.transport.fields}

    def _linear(self,v,state):
        """Value of the exact local operator, using the available W and W_Z."""
        J=self.Jet;cut=lambda a:a.truncate(0)
        z,r,x=map(cut,(v['z'],v['rho'],v['x']));c=z.ctx
        de,e,P2=v['de'],v['e'],v['P2'];d,L=1-z*z,1-z*z*de
        V={name:cut(q) for name,q in state.items()}
        V.update({name+'Z':cut(self.leading._derivative(state[name])) for name in ('f','u','k','p')})
        phi,pr=[cut(q.rows[0]) for q in v['Phi'][:2]]
        U,Ur=[cut(q.rows[0]) for q in v['U'][:2]]
        phiz=cut(self.leading._derivative(v['Phi'][0].rows[0]))
        Uz=cut(self.leading._derivative(v['U'][0].rows[0]));EU=Ur*r
        Q=cut(v['Q'][0].rows[0])
        V.update(z=z,rho=r,delta=J.constant(c,de,0),epsilon=J.constant(c,e,0),P2=J.constant(c,P2,0),
            F2=cut(v['S'])/(e*e),L=L,d=d,A=phi,C=phi+pr*r,g=cut(v['g']),
            E=z*U-c.mpf('.5'),H=z*((1-de)/2)+d*U,B=Q+(1-z*U*2)/L,Q=Q,
            Zf_without_gauge=(z*phi*(-2-de)+d*phiz-z*pr*(2*r))/L,
            Zw=(z*U*(-1-de)+d*Uz-z*EU*2)/L,Euler_U=EU,
            Euler_f=V['fx']*x/2,Euler_u=V['ux']*x/2)
        return linear_rhs(V)

    def _derivatives(self,v,state):
        r,x,z=v['rho'],v['x'],v['z'];de=v['de'];d,L=1-z*z,1-z*z*de
        rf,ru,rp=self._linear(v,state);forcing=self.transport._forcing(v)
        f1=state['fx']/(2*x);u1=state['ux']/(2*x)
        k1=-u1-state['k']/r;mean1=-state['k']/r
        f2=(forcing['fx']+rf-state['fx']*(4/x))/(4*r)
        u2=(forcing['ux']+ru-state['ux']*(2/x))/(4*r)
        k2=-u2+u1/r+state['k']*(2/(r*r))
        p1=forcing['p']/(2*x)+rp
        q1=(state['ux']*z/x+state['k']*(z*(1+de))/r
             +self.leading._derivative(state['k'])*d/r)/L
        leading={}
        for name in ('Phi','Psi','Mpsi'):
            for i in (1,2):leading[name+('_rho' if i==1 else '_rhorho')]=v[name][i].cut(1)
        leading['Q_rho']=v['Q'][1].cut(1)
        return leading,dict(f_rho=Laurent({0:f1}),u_rho=Laurent({0:u1}),k_rho=Laurent({0:k1}),
            mean_u_rho=Laurent({0:mean1}),q_rho=Laurent({0:q1.truncate(0)}),
            f_rhorho=f2.cut(0),u_rhorho=u2.cut(0),k_rhorho=k2.cut(0),p_rho=p1.cut(0)),forcing,(rf,ru,rp)

    def _prepare(self,Z,s):
        key=Z,s
        if key not in self._cache:
            packet=self.transport.evaluate(Z,s)
            native=self.transport._inputs(Z,s)
            result={}
            for directed,label in ((False,'point'),(True,'truth')):
                v=self._point_sources(Z,s,native,directed);state=self._state(packet,directed)
                result[label]=dict(inputs=v,state=state,derivatives=self._derivatives(v,state))
            self._cache[key]=packet,result
        return self._cache[key]

    def _export(self,point,truth):
        if set(point.rows)!=set(truth.rows):raise ArithmeticError('Radial width sectors changed')
        return dict(original_hb_powers={str(power):self.leading._jet_summary(point.rows[power],truth.rows[power])
                                      for power in sorted(truth.rows)},
                    coefficient_semantics='derivative = sum hb^power * enclosed coefficient(Z,s)',
                    inverse_original_width_materialized=False)

    def evaluate(self,Z,s):
        c=self.ctx
        if (type(Z) is not str or type(s) is not str
                or not all(c.isfinite(c.mpf(v)) for v in (Z,s))
                or not -1<c.mpf(Z)<1 or not 0<=c.mpf(s)<=c.mpf('.5')):
            raise ValueError('Interior decimal Z and original first phase0..1/2 required')
        packet,data=self._prepare(Z,s)
        pl,pn,_,_=data['point']['derivatives'];tl,tn,_,_=data['truth']['derivatives']
        return dict(Z=Z,phase=s,source_family=self.family,
            leading_radial_derivatives={name:self._export(value,tl[name]) for name,value in pl.items()},
            n1_radial_derivatives={name:self._export(value,tn[name]) for name,value in pn.items()},
            original_width=dict(definition='hb=cstar*actual_K_9_16^-100',
                log_enclosure_exact_mpf_tuples=[list(t) for t in self.leading.logh._mpi_],
                positive_width_not_selected_from_enclosure=True),
            exact_chi_rho_split='-sigma_prime/(hb*rho) + sigma_prime/rho',
            source_order_convention='Z coefficient k = D_Z^k/k!, radial derivatives are ordinary',
            all_second_radial_values_from_full_local_equation=True,
            nonzero_core_state_and_full_feedback_retained=True,
            accepted_n1_radial_increment_factors={name:packet['factored_first_collar_increments'][name]
                                                 for name in ('fx','ux')},
            limitations=dict(first_phase_max='.5',leading_second_rho_Z_order=1,
                n1_first_rho_f_u_k_mean_Z_order=1,n1_second_rho_Z_order=0,
                n1_pressure_first_rho_Z_order=0,n1_q_first_rho_Z_order=0,
                n1_q_Z2_or_pressure_second_rho=False,
                complete_Cartesian_C2_velocity_join=False,remaining_collars=False,
                global_moment_repair=False,temporal_n_ge_2_recursion=False,full_corrected_NS=False))
