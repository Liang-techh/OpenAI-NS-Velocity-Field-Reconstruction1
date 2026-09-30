"""Section 9.30 long swirl reshape with actual endpoint-normalized moments.

Finite quadrature and explicit conditional tail envelopes, not interval
certification. No target moments replace the inherited inner primitives.
"""
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_primitive import _sigma_mp
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_exit_field import cone_receipt


@lru_cache(maxsize=8)
def _nodes(order,precision):
    with mp.workdps(precision):
        return mp.gauss_quadrature(order,'legendre')


class LongReshape:
    def __init__(self,switches,*,A='1e150',logC='5e151',order=24,tail_digits=60):
        if int(order)!=order or order<8 or tail_digits<1:
            raise ValueError('Require integer quadrature order >=8 and positive tail digits')
        self.switches=switches; self.core=switches.core
        self.precision=switches.precision
        with mp.workdps(self.precision):
            self.A=mp.mpf(A); self.T=400*self.A; self.logC=mp.mpf(logC)
            if self.A<=0: raise ValueError('A must be positive')
            if self.logC!=mp.mpf(switches.comparison.bundle['shared_parameters']['logCstar']):
                raise ValueError('Reshape logCstar must match the inherited shared candidate')
            self.work_precision=self.precision+max(0,int(mp.floor(mp.log10(self.T)))+1)+30
            self.cutoff=(tail_digits+20)*mp.log(10)/mp.mpf('.1')
            self.reference_end=mp.mpf(switches.comparison.bundle['shared_parameters']['logRref'])-mp.log(110)-8
            if self.reference_end<=self.T:
                raise ValueError('Rsh must precede the axial restoration radius Rz')
        self.order=order; self.tail_digits=tail_digits

    @lru_cache(maxsize=16)
    def _start(self,Zkey):
        with mp.workdps(self.work_precision):
            z=mp.mpf(Zkey); v=self.switches.evaluate_R('110',z)
            logu=mp.log(v['F'])+mp.log(220)/2
            zeta=v['FZ']/v['F']; target=-2*z/(1+z*z)
            B=logu+self.logC+mp.log(1+z*z)
            BZ=zeta-target
            # This guard verifies the value and first Z derivative at this Z.
            # It does not assert the source mixed C3 bound globally.
            if max(abs(B),abs(BZ))>2*self.A:
                raise ArithmeticError('Local B/BZ exceeds provisional source budget')
            return dict(v=v,logu=logu,zeta=zeta,target=target,B=B,BZ=BZ,
                        zeta_bound=max(abs(zeta),abs(target)))

    def _shape(self,y,start):
        s=y/self.T; sig=_sigma_mp(s)
        if 0<s<1:
            ds=sig*(1-sig)*(2/s**3+2/(1-s)**3)
        else: ds=mp.mpf(0)
        logu=start['logu']+y/10-sig*start['B']
        zeta=start['zeta']-sig*start['BZ']
        a=mp.mpf('.8')+2*ds*start['B']/self.T
        return logu,zeta,a,sig

    def evaluate_log_offset(self,y,Z):
        """Evaluate y=log(R/110), preserving enormous log-radius offsets."""
        with mp.workdps(self.work_precision):
            y=mp.mpf(y); z=mp.mpf(Z)
            if not 0<=y<=self.reference_end or not abs(z)<1:
                raise ValueError('Require 110<=R<=Rz and |Z|<1')
            if y>self.T:
                return self._reference(y,z)
            start=self._start(mp.nstr(z,self.precision)); v=start['v']
            logu,zeta,a,sig=self._shape(y,start)
            R=110*mp.exp(y); u=mp.exp(logu); F=u/mp.sqrt(2*R)
            V=v['Uz']; VZ=v['UZ']; end=min(y,self.cutoff)
            # Integrate in ell=log(R/current integration radius). The
            # normalized integrands decay even when y is 4e152.
            cuts=[mp.mpf(0)]; nextcut=mp.mpf(1)
            while nextcut<end:
                cuts.append(nextcut); nextcut*=2
            if end>0: cuts.append(end)
            totals=[mp.mpf(0) for _ in range(6)]
            nodes,weights=_nodes(self.order,self.work_precision)
            for left,right in zip(cuts,cuts[1:]):
                mid=(left+right)/2; width=(right-left)/2
                for node,weight in zip(nodes,weights):
                    ell=mid+width*node
                    psig=_sigma_mp((y-ell)/self.T)
                    ratio=mp.exp(-ell/10-(psig-sig)*start['B'])
                    pzeta=start['zeta']-psig*start['BZ']
                    theta=mp.exp(-mp.mpf('1.5')*ell)*ratio
                    pressure=ratio**2; swirl=mp.exp(-ell)*pressure
                    for i,value in enumerate((theta,theta*pzeta,pressure,
                                               2*pressure*pzeta,swirl,2*swirl*pzeta)):
                        totals[i]+=width*weight*value
            t,tz,p,pz,s,sz=totals
            theta_scale=mp.sqrt(2)*R**mp.mpf('1.5')*u
            pressure_scale=u*u/2; swirl_scale=R*pressure_scale
            dt=theta_scale*t; dtz=theta_scale*tz
            dp=pressure_scale*p; dpz=pressure_scale*pz
            ds=swirl_scale*s; dsz=swirl_scale*sz
            m=dict(v['moments']); mz=dict(v['momentsZ'])
            m['theta']+=dt; mz['theta']+=dtz
            m['theta_z']+=V*dt; mz['theta_z']+=VZ*dt+V*dtz
            m['z']+=V*(R-110); mz['z']+=VZ*(R-110)
            m['p']+=dp; mz['p']+=dpz
            raw=v['raw_quadratic_integrals']
            axial=raw['axial']+V*V*(R-110)
            axialZ=raw['axial_Z']+2*V*VZ*(R-110)
            swirl=raw['swirl']+ds; swirlZ=raw['swirl_Z']+dsz
            m['z_theta']=axial-swirl; mz['z_theta']=axialZ-swirlZ
            P=v['P']+dp; PZ=v['PZ']+dpz
            stress=evaluate_mp_stress(mp.log(110)+y,z,self.core.delta,
                Utheta=u,Uz=V,Utheta_y=u*(1-a)/2,Utheta_Z=u*zeta,
                Uz_y=0,Uz_Z=VZ,moments=m,moments_Z=mz,P=P,P_Z=PZ,
                shear_theta=-a*F,shear_z=0,precision=self.work_precision)
            # From |B|<=2A and sup sigma'<=8: d_y logu >= .06.
            # Use the weaker source .05 to expose conservative tails.
            omitted=y>end
            tails={key:mp.exp(-rate*end)/rate if omitted else mp.mpf(0)
                   for key,rate in [('theta',mp.mpf('1.55')),
                                    ('p',mp.mpf('.1')),('swirl',mp.mpf('1.1'))]}
            result=dict(R=R,logR=mp.log(110)+y,log_radius_offset=y,Z=z,
                F=F,FZ=F*zeta,Uz=V,UZ=VZ,P=P,PZ=PZ,Ur=stress['U_r'],
                moments=m,momentsZ=mz,stress=stress,a=a,b=mp.mpf(0),
                g_y=-a/2,u_y=mp.mpf(0),F_R=-a*F/(2*R),Uz_R=mp.mpf(0),
                normalized_integrals=dict(theta=t,p=p,swirl=s),
                normalized_integrals_scope='increments from R110, local endpoint normalization',
                normalized_tail_bounds=tails,
                normalized_Z_tail_bounds={k:tails[k]*start['zeta_bound']*(1 if k=='theta' else 2) for k in tails},
                raw_quadratic_integrals=dict(axial=axial,axial_Z=axialZ,swirl=swirl,swirl_Z=swirlZ),
                quadrature_order=self.order,quadrature_error_certified=False,
                tail_assumptions_locally_checked=True,source_constants_certified=False,
                region='long_swirl_reshape',outer_matching_complete=False)
            return result

    @lru_cache(maxsize=16)
    def _shape_endpoint(self,Zkey):
        return self.evaluate_log_offset(self.T,mp.mpf(Zkey))

    def _reference(self,y,z):
        """Exact power-law primitives from Rsh through Rz; no new tails."""
        v=self._shape_endpoint(mp.nstr(z,self.precision)); x=y-self.T
        R=v['R']*mp.exp(x); deltaR=v['R']*mp.expm1(x)
        u=mp.sqrt(2*v['R'])*v['F']*mp.exp(x/10)
        F=u/mp.sqrt(2*R); zeta=-2*z/(1+z*z)
        V=v['Uz']; VZ=v['UZ']
        dt=mp.sqrt(2)*R**mp.mpf('1.5')*u*(-mp.expm1(-mp.mpf('1.6')*x))/mp.mpf('1.6')
        dp=u*u/2*(-mp.expm1(-mp.mpf('.2')*x))/mp.mpf('.2')
        ds=R*u*u/2*(-mp.expm1(-mp.mpf('1.2')*x))/mp.mpf('1.2')
        m=dict(v['moments']); mz=dict(v['momentsZ'])
        m['theta']+=dt; mz['theta']+=zeta*dt
        m['theta_z']+=V*dt; mz['theta_z']+=(VZ+V*zeta)*dt
        m['z']+=V*deltaR; mz['z']+=VZ*deltaR
        m['p']+=dp; mz['p']+=2*zeta*dp
        raw=v['raw_quadratic_integrals']
        axial=raw['axial']+V*V*deltaR; axialZ=raw['axial_Z']+2*V*VZ*deltaR
        swirl=raw['swirl']+ds; swirlZ=raw['swirl_Z']+2*zeta*ds
        m['z_theta']=axial-swirl; mz['z_theta']=axialZ-swirlZ
        P=v['P']+dp; PZ=v['PZ']+2*zeta*dp
        a=mp.mpf('.8')
        stress=evaluate_mp_stress(mp.log(110)+y,z,self.core.delta,
            Utheta=u,Uz=V,Utheta_y=u/10,Utheta_Z=u*zeta,Uz_y=0,Uz_Z=VZ,
            moments=m,moments_Z=mz,P=P,P_Z=PZ,shear_theta=-a*F,shear_z=0,
            precision=self.work_precision)
        return dict(R=R,logR=mp.log(110)+y,log_radius_offset=y,Z=z,
            F=F,FZ=F*zeta,Uz=V,UZ=VZ,P=P,PZ=PZ,Ur=stress['U_r'],
            moments=m,momentsZ=mz,stress=stress,a=a,b=mp.mpf(0),
            g_y=-a/2,u_y=mp.mpf(0),F_R=-a*F/(2*R),Uz_R=mp.mpf(0),
            raw_quadratic_integrals=dict(axial=axial,axial_Z=axialZ,swirl=swirl,swirl_Z=swirlZ),
            normalized_integrals=v['normalized_integrals'],
            normalized_integrals_scope='inherited shaping increments at Rsh; later primitives exact',
            normalized_tail_bounds=v['normalized_tail_bounds'],
            normalized_Z_tail_bounds=v['normalized_Z_tail_bounds'],
            tail_bounds_scope='inherited Rsh normalization; reference primitives exact',
            quadrature_error_certified=False,source_constants_certified=False,
            region='reference_swirl_before_axial_restore',outer_matching_complete=False)

    def evaluate_phase(self,phase,Z):
        with mp.workdps(self.work_precision):
            phase=mp.mpf(phase)
            if not 0<=phase<=1: raise ValueError('Shaping phase must lie in [0,1]')
            return self.evaluate_log_offset(phase*self.T,Z)


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
        switches=ExitSwitches(ExitContinuation(ExitTangents(comparison,epsilon=hb,
                              steps=16,derivative_step='1e-45')),steps=16)
        coarse=LongReshape(switches,order=16); fine=LongReshape(switches,order=32)
        z=mp.mpf('.3'); receipt=dict(source='arXiv:2609.35406v1 Eq9.30',
            shared_parameters=bundle['shared_parameters'],rows=[],
            source_constants_certified=False,finite_energy_certified=False,
            scale_recursion_established=False)
        path=Path(__file__).with_suffix('.json')
        for name,y in [('start',mp.mpf(0)),('one_log_unit',mp.mpf(1)),
                       ('midpoint',fine.T/2),('end',fine.T),
                       ('axial_restore_entry',fine.reference_end)]:
            print('computing '+name,flush=True)
            lo=coarse.evaluate_log_offset(y,z); hi=fine.evaluate_log_offset(y,z)
            with mp.workdps(fine.work_precision):
                errors={k:abs(lo[k]-hi[k])/max(abs(hi[k]),mp.mpf('1e-100'))
                        for k in ('Uz','Ur','P','PZ')}
                errors.update({k:abs(lo['normalized_integrals'][k]-hi['normalized_integrals'][k])
                    /max(abs(hi['normalized_integrals'][k]),mp.mpf('1e-100'))
                    for k in ('theta','p','swirl')})
                if name=='start':
                    old=switches.evaluate_R('110',z)
                    relative=lambda left,right: abs(left-right)/abs(right) if right else abs(left)
                    receipt['initial_match_errors']={k:mp.nstr(relative(hi[k],old[k]),45)
                        for k in ('F','FZ','Uz','UZ','Ur','P','PZ','g_y','u_y')}
                    receipt['initial_moment_match_errors']={kind:{k:mp.nstr(relative(hi[kind][k],old[kind][k]),45)
                        for k in hi[kind]} for kind in ('moments','momentsZ','raw_quadratic_integrals')}
                row=dict(name=name,y=mp.nstr(y,260),a=mp.nstr(hi['a'],45),
                    region=hi['region'],
                    normalized_integrals_scope=hi['normalized_integrals_scope'],
                    F_positive=bool(hi['F']>0),cone=cone_receipt(hi),
                    normalized_integrals={k:mp.nstr(v,45) for k,v in hi['normalized_integrals'].items()},
                    normalized_tail_bounds={k:mp.nstr(v,45) for k,v in hi['normalized_tail_bounds'].items()},
                    normalized_Z_tail_bounds={k:mp.nstr(v,45) for k,v in hi['normalized_Z_tail_bounds'].items()},
                    tail_bounds_scope=hi.get('tail_bounds_scope','local endpoint normalization'),
                    quadrature_refinement={k:mp.nstr(v,45) for k,v in errors.items()})
                if name in ('end','axial_restore_entry'):
                    target=-fine.logC-mp.log(1+z*z)+y/10
                    row['reference_logu_error']=mp.nstr(mp.log(hi['F'])+mp.log(2*hi['R'])/2-target,45)
                receipt['rows'].append(row)
                path.write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
                print(json.dumps(row),flush=True)
        return receipt


if __name__=='__main__':run()
