"""One finite core polynomial supplies velocity, pressure, moments and stress.

The coefficient factory owns axis data and its uncertainty. There is no
implicit join to the source outer field. Queries outside the core fail.
"""
from functools import lru_cache
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_core_recursion import evaluate_core_jets
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress


class CorePolynomial:
    def __init__(self,coefficient_factory,*,Lambda,delta='1e-32',precision=160):
        self.factory=coefficient_factory; self.precision=precision
        with mp.workdps(precision):
            self.Lambda=mp.mpf(str(Lambda)); self.delta=mp.mpf(str(delta))
        if self.Lambda<=0:raise ValueError('Lambda must be positive')

    @lru_cache(maxsize=32)
    def coefficients(self,Z):
        return self.factory(Z)

    def evaluate(self,R,Z):
        with mp.workdps(self.precision):
            r=mp.mpf(str(R)); z=mp.mpf(str(Z))
            if not 0<=r<=mp.mpf('4.1')/self.Lambda or abs(z)>1:
                raise ValueError('Query must remain in 0<=Lambda R<=4.1, |Z|<=1')
            coefficients=self.coefficients(Z)
            a=evaluate_core_jets(coefficients,r)
            f=coefficients['F']; u=coefficients['Uz']
            m={k:mp.mpf(0) for k in ('theta','z','theta_z','z_theta','p')}
            mz=dict(m)
            for n,row in enumerate(f):
                m['theta']+=2*row[0]*r**(n+2)/(n+2)
                mz['theta']+=2*row[1]*r**(n+2)/(n+2)
            for n,row in enumerate(u):
                m['z']+=row[0]*r**(n+1)/(n+1)
                mz['z']+=row[1]*r**(n+1)/(n+1)
            for i,fi in enumerate(f):
                for j,uj in enumerate(u):
                    power=i+j+2
                    m['theta_z']+=2*fi[0]*uj[0]*r**power/power
                    mz['theta_z']+=2*(fi[1]*uj[0]+fi[0]*uj[1])*r**power/power
                for j,fj in enumerate(f):
                    power=i+j+1
                    value=fi[0]*fj[0]; derivative=fi[1]*fj[0]+fi[0]*fj[1]
                    m['p']+=value*r**power/power
                    mz['p']+=derivative*r**power/power
                    m['z_theta']-=value*r**(power+1)/(power+1)
                    mz['z_theta']-=derivative*r**(power+1)/(power+1)
            for i,ui in enumerate(u):
                for j,uj in enumerate(u):
                    power=i+j+1
                    m['z_theta']+=ui[0]*uj[0]*r**power/power
                    mz['z_theta']+=(ui[1]*uj[0]+ui[0]*uj[1])*r**power/power
            # Full FÃ‚Â² integral, rather than its truncated Taylor pressure.
            a['P']=coefficients['P'][0][0]+m['p']
            a['P_Z']=coefficients['P'][0][1]+mz['p']
            a['P_R']=a['F']**2
            root=mp.sqrt(2*r)
            a['Utheta']=root*a['F']; a['moments']=m; a['moments_Z']=mz
            a['Ur']=((2*z*r*a['Uz']-(1-self.delta)*z*m['z']-(1-z*z)*mz['z'])/
                     ((1-self.delta*z*z)*root)) if r else mp.mpf(0)
            return a

    def stress(self,R,Z):
        with mp.workdps(self.precision):
            r=mp.mpf(str(R)); a=self.evaluate(r,Z)
            if r<=0:raise ValueError('Stress formula is evaluated at R>0')
            root=mp.sqrt(2*r)
            return evaluate_mp_stress(mp.log(r),Z,self.delta,
                Utheta=a['Utheta'],Uz=a['Uz'],Utheta_y=root*(r*a['F_R']+a['F']/2),
                Utheta_Z=root*a['F_Z'],Uz_y=r*a['Uz_R'],Uz_Z=a['Uz_Z'],
                moments=a['moments'],moments_Z=a['moments_Z'],P=a['P'],P_Z=a['P_Z'],
                precision=self.precision)

    def physical_chart(self,R,Z,logq,*,nu='.01',phi=0):
        """Cartesian coordinates and velocity of this same regular core."""
        with mp.workdps(self.precision):
            r=mp.mpf(str(R)); z=mp.mpf(str(Z)); q=mp.exp(mp.mpf(str(logq)))
            if abs(z)>=1:raise ValueError('Physical chart requires |Z|<1')
            viscosity=mp.mpf(str(nu)); angle=mp.mpf(str(phi))
            if viscosity<=0:raise ValueError('nu must be positive')
            a=self.evaluate(r,Z); radius=mp.sqrt(2*viscosity*q*r)
            axial=mp.sqrt(viscosity)*q**((1-self.delta)/2)*z
            ur=mp.sqrt(viscosity/q)*a['Ur']
            scale=mp.sqrt(viscosity)*q**(-(1+self.delta)/2)
            ut=scale*a['Utheta']; uz=scale*a['Uz']
            return {'xyz':(radius*mp.cos(angle),radius*mp.sin(angle),axial),
                    'uvw':(ur*mp.cos(angle)-ut*mp.sin(angle),ur*mp.sin(angle)+ut*mp.cos(angle),uz),
                    'tau':q*(1-z*z),'pressure':viscosity*q**(-1-self.delta)*a['P'],
                    'scope':'Unlocalized finite core; outer connection and pressure-jet uncertainties remain.'}


def build_source_core(precision=160,degree=18,Lambda='1e36',*,
                      j='.02',logC=None,logPstar='14',delta='1e-32',
                      continuous_pressure=False,pressure_order=192):
    from lei_ren_part1_paper_outer import PaperOuterSchedule
    from lei_ren_part1_paper_corrected_profile import CorrectedSourceProfile
    from lei_ren_part1_paper_axis_pressure_jets import AxisPressureJets
    from lei_ren_part1_paper_axis_jets import RegularCoreAxisJets
    from lei_ren_part1_paper_core_ra_experiment import build_coefficients
    from lei_ren_part1_paper_axial_correction import signed_log
    with mp.workdps(precision):
        lam=mp.mpf(str(Lambda))
        chosen_logC=2*mp.log(lam) if logC is None else mp.mpf(str(logC))
        logP=mp.mpf(str(logPstar)); chosen_delta=mp.mpf(str(delta))
        logR=mp.log(110)+10*(chosen_logC+logP)
        schedule=PaperOuterSchedule(logPstar=mp.nstr(logP,precision),logRref=mp.nstr(logR,precision),
            delta=mp.nstr(chosen_delta,precision),Md='.5',c_mu='.001',c_delta='.001',c_epsilon='.01')
        profile=CorrectedSourceProfile(schedule=schedule,match_waiting=True,precision=precision)
        if continuous_pressure:
            from lei_ren_part1_paper_continuous_angular_schedule import install_continuous_angular_schedule
            from lei_ren_part1_paper_continuous_axis_pressure import ContinuousAxisPressureJets
            # Build the nonlinear inner core with this datum from the start;
            # never shift a completed field's Z-dependent pressure afterward.
            # CorrectedSourceProfile may replace the schedule when solving
            # its waiting length. Install on the actual profile-owned object.
            install_continuous_angular_schedule(profile.schedule,precision=precision,primitive_precision=100)
            pressure=ContinuousAxisPressureJets(profile,R_a=4/lam,quadrature_order=pressure_order)
        else:
            pressure=AxisPressureJets(profile,R_a=4/lam)
        anchor=pressure.dominant_taylor(0)
        # AxisPressureJets separates P0/PstarÃ‚Â². Restore physical profile
        # units before using the Section 8 nonlinear equations.
        K=mp.exp(2*logP)*mp.mpf(anchor['anchor_K_actual_Z0'])
        axis=RegularCoreAxisJets(j=j,Lambda=lam,logC=chosen_logC,delta=chosen_delta,precision=precision)
        def factory(z):
            center=mp.mpf(str(z))
            jets=mp.taylor(lambda w:K/(1+w*w)**2,center,degree+1)
            return build_coefficients(center,degree,axis=axis,
                pressure_taylor_coefficients=jets,precision=precision)
        core=CorePolynomial(factory,Lambda=lam,delta=chosen_delta,precision=precision)
        return {'profile':profile,'axis':axis,'pressure':pressure,'core':core,'K':K,
                'Lambda':lam,'radial_degree':degree,'precision':precision,
                'shared_parameters':{'j':mp.nstr(axis.j,precision),
                    'logCstar':mp.nstr(chosen_logC,precision),'logPstar':mp.nstr(logP,precision),
                    'delta':mp.nstr(chosen_delta,precision),'logRref':mp.nstr(logR,precision)},
                'pressure_recomputed_for_shared_parameters':True,
                'continuous_preflatten_pressure_anchor':continuous_pressure,
                'post_Rv_pressure_tail_in_axis_jet':False,
                'source_parameter_regime_certified':False}


def run(precision=160,degree=18,Lambda='1e36'):
    from lei_ren_part1_paper_axial_correction import signed_log
    with mp.workdps(precision):
        bundle=build_source_core(precision,degree,Lambda)
        core=bundle['core']; axis=bundle['axis']; pressure=bundle['pressure']
        K=bundle['K']; lam=bundle['Lambda']
        z='.3'; rows=[]
        center=axis.Z0; d=1-center*center; L=1-axis.delta*center*center
        U0=4*center+axis.j; H=axis.H0(center)
        P0=K/(1+center*center)**2; P0Z=-4*K*center/(1+center*center)**3
        U1=((1+axis.delta)/2*(1-2*center*U0)*U0+4*H+d*P0Z
            -2*(1+axis.delta)*center*P0)/(2*L)
        width=abs(U1)/(lam*abs(axis.H0_Z(center)))
        probes=(z,mp.nstr(center,precision),'-0.00444443468970541')+tuple(
            mp.nstr(center+mp.mpf(c)*width,precision) for c in (-16,-4,-1,1,4,16))
        for point in probes:
            for s in ('.5','1','2','4','4.1'):
                r=mp.mpf(s)/lam; a=core.evaluate(r,point); stress=core.stress(r,point)
                kappa=-((2*r*a['F_R'])**2+2*r*a['Uz_R']**2)/(a['F']*2*r*a['F_R'])
                rows.append({'Z':point,'s':s,'F_positive':bool(a['F']>0),'F_R_negative':bool(a['F_R']<0),
                    'F_R_over_F':mp.nstr(a['F_R']/a['F'],40),
                    'kappa':signed_log(kappa,precision),
                    'T_theta_over_S_theta':mp.nstr(abs(stress['T_theta']/stress['S_theta']),30),
                    'T_z_over_S_z':mp.nstr(abs(stress['T_z']/stress['S_z']),30),
                    'Uz':mp.nstr(a['Uz'],precision)})
        field=core.physical_chart(2/lam,z,'-4',phi='.4')
        axisfield=core.physical_chart(0,z,'-4')
        signs_passed=all(row['F_positive'] and row['F_R_negative'] for row in rows)
        return {'radial_degree':degree,'Lambda':mp.nstr(lam,precision),'Z':z,
            'samples':rows,
            'all_probe_signs_passed':signs_passed,
            'narrow_pressure_window':mp.nstr(width,precision),
            'pressure_to_Lambda_ratio_U1_squared_over_Lambda_sigma_squared':mp.nstr(U1**2/(lam*axis.sigma0**2),40),
            'physical_core_sample':{k:[signed_log(v,precision) for v in field[k]] for k in ('xyz','uvw')},
            'axis_velocity': [signed_log(v,precision) for v in axisfield['uvw']],
            'actual_axis_pressure_decomposition':pressure.pressure_receipt(z),
            'all_moments_from_same_polynomial':True,'pressure_radial_identity_exact_for_polynomial':True,
            'future_pressure_tail_derivatives_bounded':False,'outer_connection_complete':False,
            'global_finite_energy_certified':False,'scale_recursion_established':False}


if __name__=='__main__':
    result=run()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'all_probe_signs_passed':result['all_probe_signs_passed'],
                      'sample_count':len(result['samples']),
                      'pressure_to_Lambda_ratio':result['pressure_to_Lambda_ratio_U1_squared_over_Lambda_sigma_squared'],
                      'max_theta_stress_ratio':max(float(row['T_theta_over_S_theta']) for row in result['samples'])}))
