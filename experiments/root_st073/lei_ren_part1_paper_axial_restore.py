"""Actual Section 9.38 axial restoration; accumulated moments reach Rh."""
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_primitive import _sigma_mp
from lei_ren_part1_paper_long_reshape import LongReshape, _nodes
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_exit_field import cone_receipt


class AxialRestore:
    def __init__(self,reshape,*,order=64):
        if int(order)!=order or order<8: raise ValueError('Require integer order >=8')
        self.reshape=reshape; self.core=reshape.core
        self.precision=reshape.work_precision; self.order=order

    @lru_cache(maxsize=16)
    def _start(self,Zkey):
        with mp.workdps(self.precision):
            return self.reshape.evaluate_log_offset(self.reshape.reference_end,mp.mpf(Zkey))

    @lru_cache(maxsize=32)
    def _integrals(self,tkey):
        with mp.workdps(self.precision):
            t=mp.mpf(tkey); out=[mp.mpf(0)]*3
            if t==0:return out
            nodes,weights=_nodes(self.order,self.precision)
            # Flat endpoints at 0 and 1; split keeps all intervals bounded.
            for left,right in ((mp.mpf(0),t/2),(t/2,t)):
                mid=(left+right)/2; half=(right-left)/2
                for node,weight in zip(nodes,weights):
                    s=mid+half*node; sig=_sigma_mp(s)
                    values=(mp.exp(-(t-s))*sig,
                            mp.exp(-mp.mpf('1.6')*(t-s))*sig,
                            mp.exp(-(t-s))*sig*sig)
                    out=[a+half*weight*b for a,b in zip(out,values)]
            return out

    def _output(self,t,z,v,V,VZ,Vy,VyZ,dt,dtZ,mass,massZ,mixed,mixedZ,axial,axialZ):
        R=v['R']*mp.exp(t); root=mp.sqrt(2*R)
        u=mp.sqrt(2*v['R'])*v['F']*mp.exp(t/10)
        F=u/root; zeta=-2*z/(1+z*z)
        dp=u*u/2*(-mp.expm1(-mp.mpf('.2')*t))/mp.mpf('.2')
        ds=R*u*u/2*(-mp.expm1(-mp.mpf('1.2')*t))/mp.mpf('1.2')
        m=dict(v['moments']); mz=dict(v['momentsZ'])
        m['theta']+=dt; mz['theta']+=dtZ
        m['z']+=mass; mz['z']+=massZ
        m['theta_z']+=mixed; mz['theta_z']+=mixedZ
        m['p']+=dp; mz['p']+=2*zeta*dp
        raw=v['raw_quadratic_integrals']
        axial+=raw['axial']; axialZ+=raw['axial_Z']
        swirl=raw['swirl']+ds; swirlZ=raw['swirl_Z']+2*zeta*ds
        m['z_theta']=axial-swirl; mz['z_theta']=axialZ-swirlZ
        P=v['P']+dp; PZ=v['PZ']+2*zeta*dp
        a=mp.mpf('.8'); b=2*Vy/u
        stress=evaluate_mp_stress(v['logR']+t,z,self.core.delta,
            Utheta=u,Uz=V,Utheta_y=u/10,Utheta_Z=u*zeta,
            Uz_y=Vy,Uz_Z=VZ,moments=m,moments_Z=mz,P=P,P_Z=PZ,
            shear_theta=-a*F,shear_z=root*Vy/R,precision=self.precision)
        return dict(R=R,logR=v['logR']+t,Z=z,F=F,FZ=F*zeta,Uz=V,UZ=VZ,
            Ur=stress['U_r'],P=P,PZ=PZ,moments=m,momentsZ=mz,stress=stress,
            a=a,b=b,g_y=-a/2,u_y=Vy,u_y_Z=VyZ,F_R=-a*F/(2*R),Uz_R=Vy/R,
            raw_quadratic_integrals=dict(axial=axial,axial_Z=axialZ,swirl=swirl,swirl_Z=swirlZ),
            region='axial_restoration' if t<=1 else 'reference_before_moment_repair',
            axial_restore_phase=t,source_constants_certified=False,
            moment_repair_complete=False,outer_matching_complete=False)

    def evaluate_phase(self,t,Z):
        """t=log(R/Rz): restoration [0,1], reference continuation [1,3]."""
        with mp.workdps(self.precision):
            t=mp.mpf(t); z=mp.mpf(Z)
            if not 0<=t<=3 or not abs(z)<1:
                raise ValueError('Require 0<=log(R/Rz)<=3 and |Z|<1')
            if t>1:return self._continue(t,z)
            v=self._start(mp.nstr(z,self.reshape.precision)); V0=v['Uz']; V0Z=v['UZ']
            delta=4*z-V0; deltaZ=4-V0Z
            sig=_sigma_mp(t)
            sigprime=sig*(1-sig)*(2/t**3+2/(1-t)**3) if 0<t<1 else mp.mpf(0)
            V=V0+delta*sig; VZ=V0Z+deltaZ*sig
            R=v['R']*mp.exp(t); deltaR=v['R']*mp.expm1(t)
            u=mp.sqrt(2*v['R'])*v['F']*mp.exp(t/10)
            scale=mp.sqrt(2)*R**mp.mpf('1.5')*u; zeta=-2*z/(1+z*z)
            dt=scale*(-mp.expm1(-mp.mpf('1.6')*t))/mp.mpf('1.6')
            j1,j16,j2=self._integrals(mp.nstr(t,self.precision))
            mass=V0*deltaR+delta*R*j1
            massZ=V0Z*deltaR+deltaZ*R*j1
            mixed=V0*dt+delta*scale*j16
            mixedZ=(V0Z+V0*zeta)*dt+(deltaZ+delta*zeta)*scale*j16
            axial=V0*V0*deltaR+2*V0*delta*R*j1+delta*delta*R*j2
            axialZ=2*V0*V0Z*deltaR+2*(V0Z*delta+V0*deltaZ)*R*j1+2*delta*deltaZ*R*j2
            return self._output(t,z,v,V,VZ,delta*sigprime,deltaZ*sigprime,
                dt,zeta*dt,mass,massZ,mixed,mixedZ,axial,axialZ)

    @lru_cache(maxsize=16)
    def _restored(self,Zkey):return self.evaluate_phase(1,mp.mpf(Zkey))

    def _continue(self,t,z):
        v=self._restored(mp.nstr(z,self.reshape.precision)); x=t-1
        R=v['R']*mp.exp(x); deltaR=v['R']*mp.expm1(x)
        u=mp.sqrt(2*v['R'])*v['F']*mp.exp(x/10)
        dt=mp.sqrt(2)*R**mp.mpf('1.5')*u*(-mp.expm1(-mp.mpf('1.6')*x))/mp.mpf('1.6')
        zeta=-2*z/(1+z*z); V=4*z; VZ=mp.mpf(4)
        out=self._output(x,z,v,V,VZ,mp.mpf(0),mp.mpf(0),dt,zeta*dt,
            V*deltaR,VZ*deltaR,V*dt,(VZ+V*zeta)*dt,V*V*deltaR,2*V*VZ*deltaR)
        out['axial_restore_phase']=t; out['region']='reference_before_moment_repair'
        return out


def reference_moments(R,u,z):
    """Supplied Section 9.1 outer target; pressure offset stays separate."""
    theta=mp.mpf('5')/8*R*mp.sqrt(2*R)*u
    return dict(theta=theta,z=4*z*R,theta_z=4*z*theta,
                z_theta=16*z*z*R-mp.mpf('5')/12*R*u*u,p=mp.mpf('2.5')*u*u)


def defect_receipt(provider,Z):
    """Expose numerical subtraction losses; never certify zero from cancellation."""
    with mp.workdps(provider.precision):
        z=mp.mpf(Z); v=provider.evaluate_phase(3,z)
        R=v['R']; u=mp.sqrt(2*R)*v['F']; Rm=R/mp.e; Am=u*mp.exp(-mp.mpf('.1'))
        outer=reference_moments(R,u,z)
        delta={k:v['moments'][k]-outer[k] for k in outer}
        zeta=-2*z/(1+z*z)
        outerZ=dict(theta=zeta*outer['theta'],z=4*R,
            theta_z=4*outer['theta']+4*z*zeta*outer['theta'],
            z_theta=32*z*R-mp.mpf('5')/6*R*u*u*zeta,p=5*u*u*zeta)
        deltaZ={k:v['momentsZ'][k]-outerZ[k] for k in outerZ}
        centered=[delta['z']/Rm,(delta['theta_z']-4*z*delta['theta'])/(mp.sqrt(2)*Rm**mp.mpf('1.5')*Am),
            delta['theta']/(mp.sqrt(2)*Rm**mp.mpf('1.5')*Am),
            (delta['z_theta']-8*z*delta['z'])/(Rm*Am*Am),delta['p']/(Am*Am)]
        centeredZ=[deltaZ['z']/Rm,
            (deltaZ['theta_z']-4*delta['theta']-4*z*deltaZ['theta'])/(mp.sqrt(2)*Rm**mp.mpf('1.5')*Am)-zeta*centered[1],
            deltaZ['theta']/(mp.sqrt(2)*Rm**mp.mpf('1.5')*Am)-zeta*centered[2],
            (deltaZ['z_theta']-8*delta['z']-8*z*deltaZ['z'])/(Rm*Am*Am)-2*zeta*centered[3],
            deltaZ['p']/(Am*Am)-2*zeta*centered[4]]
        # Angular differences can be far below the subtraction precision.
        # Source 10.18 bounds their changes independently of that loss.
        reshape=provider.reshape
        old=reshape.switches.evaluate_R('110',z)
        uold=mp.exp(-reshape.logC)/(1+z*z)
        targets0=reference_moments(mp.mpf(110),uold,z)
        delta0={k:old['moments'][k]-targets0[k] for k in targets0}
        scales=(Rm,mp.sqrt(2)*Rm**mp.mpf('1.5')*Am,
                mp.sqrt(2)*Rm**mp.mpf('1.5')*Am,Rm*Am*Am,Am*Am)
        d0=[delta0['z']/scales[0],(delta0['theta_z']-4*z*delta0['theta'])/scales[1],
            delta0['theta']/scales[2],(delta0['z_theta']-8*z*delta0['z'])/scales[3],delta0['p']/scales[4]]
        alpha=mp.exp(mp.log(110)+reshape.T-mp.log(Rm)); s=mp.exp(-1)
        entry=provider._start(mp.nstr(z,reshape.precision))
        eta=abs(entry['Uz']-4*z)+abs(entry['UZ']-4)
        logP=mp.mpf(reshape.switches.comparison.bundle['shared_parameters']['logPstar'])
        changes=[eta*s,2*eta*s**mp.mpf('1.6'),2*alpha**mp.mpf('1.6'),
                 alpha**mp.mpf('1.2')+40*eta*eta*s*mp.exp(-2*logP),10*alpha**mp.mpf('.2')]
        bounds=[abs(x)+y for x,y in zip(d0,changes)]
        # A subtraction can return a nonzero roundoff artifact larger than
        # the conditional physical bound. Such entries are unresolved too.
        unresolved=[i+1 for i,(value,bound) in enumerate(zip(centered,bounds))
                    if value==0 or abs(value)>bound]
        n=lambda value:mp.nstr(value,60)
        return dict(Z=n(z),centered_defect_order=['z','theta_z-4Ztheta','theta','z_theta-8Zz','p'],
            centered_defects_from_finite_subtraction=[n(x) for x in centered],
            centered_Z_derivatives_from_finite_subtraction=[n(x) for x in centeredZ],
            cancellation_unresolved_entries=unresolved,
            finite_subtraction_is_not_an_exact_defect_certificate=True,
            zeros_are_not_certified=True,
            conditional_centered_value_bounds=[n(x) for x in bounds],
            bound_source='Eq10.18 with pointwise eta; uniform C1 assumptions unproved',
            sampled_eta=n(eta),input_e_star_certified=False,
            pressure_axis_offset_preserved=True,moment_repair_complete=False)


def run():
    from lei_ren_part1_paper_core_adapter import build_source_core
    from lei_ren_part1_paper_exit_comparison import Section923Comparison
    from lei_ren_part1_paper_exit_tangents import ExitTangents
    from lei_ren_part1_paper_exit_continuation import ExitContinuation
    from lei_ren_part1_paper_exit_switches import ExitSwitches
    with mp.workdps(260):
        bundle=build_source_core(precision=260,degree=18,j='1e-14',Lambda='1e36',
            logC='5e151',logPstar='14',delta='1e-200')
        hb=mp.exp(-100-100*mp.mpf('1e152'))
        comparison=Section923Comparison(bundle,h_b=hb,transition_steps=16)
        switches=ExitSwitches(ExitContinuation(ExitTangents(comparison,epsilon=hb,steps=16,
            derivative_step='1e-45')),steps=16)
        reshape=LongReshape(switches,order=32)
        coarse=AxialRestore(reshape,order=32); fine=AxialRestore(reshape,order=64)
        receipt=dict(shared_parameters=bundle['shared_parameters'],rows=[],
            source='arXiv:2609.35406v1 Eq9.38',finite_energy_certified=False,
            scale_recursion_established=False,moment_repair_complete=False)
        path=Path(__file__).with_suffix('.json')
        with mp.workdps(fine.precision):
            z=mp.mpf('.3')
            for t in ('0','.5','1','2','3'):
                print('computing restore '+t,flush=True)
                lo=coarse.evaluate_phase(t,z); hi=fine.evaluate_phase(t,z)
                rel=lambda a,b:abs(a-b)/abs(b) if b else abs(a)
                if t=='0':
                    old=fine._start(mp.nstr(z,reshape.precision))
                    receipt['entry_matching']={k:mp.nstr(rel(hi[k],old[k]),45)
                        for k in ('F','FZ','Uz','UZ','Ur','P','PZ','u_y','g_y')}
                    receipt['entry_moment_matching']={kind:{k:mp.nstr(rel(hi[kind][k],old[kind][k]),45)
                        for k in hi[kind]} for kind in ('moments','momentsZ','raw_quadratic_integrals')}
                row=dict(t=t,Uz=mp.nstr(hi['Uz'],60),UZ=mp.nstr(hi['UZ'],60),
                    u_y=mp.nstr(hi['u_y'],60),b=mp.nstr(hi['b'],60),cone=cone_receipt(hi),
                    refinement={k:mp.nstr(rel(lo[k],hi[k]),45) for k in ('Ur','P','PZ')},
                    moment_refinement={kind:{k:mp.nstr(rel(lo[kind][k],hi[kind][k]),45)
                        for k in hi[kind]} for kind in ('moments','momentsZ')})
                receipt['rows'].append(row)
                path.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
                print(json.dumps({k:row[k] for k in ('t','Uz','u_y','b','cone')}),flush=True)
            receipt['defects']=defect_receipt(fine,z)
            # Independent radial finite differences check the new transported
            # primitives in the mapped continuity equation at an interior point.
            t=mp.mpf('.5'); center=fine.evaluate_phase(t,z); divergence=[]
            for h in (mp.mpf('1e-5'),mp.mpf('5e-6')):
                values={i:fine.evaluate_phase(t+i*h,z)['Ur'] for i in (-2,-1,1,2)}
                derivative=(values[-2]-8*values[-1]+8*values[1]-values[2])/(12*h)
                R=center['R']; delta=fine.core.delta; L=1-delta*z*z
                terms=[mp.sqrt(2*R)*derivative/R,center['Ur']/mp.sqrt(2*R),
                    ((1-z*z)*center['UZ']-2*z*center['u_y']-(1+delta)*z*center['Uz'])/L]
                residual=sum(terms)
                divergence.append(dict(h=mp.nstr(h,30),mapped_q_divergence=mp.nstr(residual,45),
                    relative_cancellation=mp.nstr(abs(residual)/sum(abs(x) for x in terms),45)))
            receipt['radial_finite_difference_divergence']=divergence
            path.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
        return receipt


if __name__=='__main__':run()
