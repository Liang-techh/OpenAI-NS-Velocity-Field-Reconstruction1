"""Independent O.2 turnoff fixture, including retained mass and radial recovery."""
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from scipy.integrate import solve_ivp
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_outer_slope_field import evaluate_transition
from lei_ren_part1_paper_interval_outer_axial_turnoff import IntervalOuterAxialTurnoff,cutoff_integrals
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress


def run():
    c=MPIntervalContext();c.dps=90;count=0
    def B(y):
        s=math.log(y)/1.1
        if s<=0:return 1.
        if s>=1:return 0.
        exponent=-1/(1-s)**2+1/s**2
        if exponent>700:return 1.
        if exponent < -700:return 0.
        return 1-1/(1+math.exp(exponent))
    ode=solve_ivp(lambda y,k:[math.exp(y)*B(y),math.exp(y)*B(y)**2],
        (1,math.exp(1.1)),[0.,0.],method='DOP853',rtol=2e-13,atol=2e-14,dense_output=True)
    if not ode.success:raise AssertionError(ode.message)
    with mp.workdps(100):
        z=IntervalTaylor(c,[c.mpf('.5'),c.mpf(1)])
        A=(1+z*z).reciprocal()*2;P0=z*z+3
        inlet=evaluate_transition(c,z,c.mpf('.01'),c.mpf(7),A,P0,1)
        inlet.update(pressure_schedule_Md='1.1',pressure_datum_kind='synthetic_fixture')
        field=IntervalOuterAxialTurnoff(c,z,c.mpf('.01'),inlet,'1.1','1.1')
        R1=7*mp.e;u1=mp.mpf('1.6')*mp.exp(-mp.mpf('.2'))
        u1z=mp.mpf('-1.28')*mp.exp(-mp.mpf('.2'))
        def contains(box,v,label):
            nonlocal count
            lo,hi=endpoints(box)
            if not lo<=v<=hi:raise AssertionError((label,v,lo,hi))
            count+=1
        for phase in ('.25','.5','.75','1'):
            p=field.evaluate_phase(phase)
            y=mp.exp(mp.mpf('1.1')*mp.mpf(phase));t=y-1;R=R1*mp.exp(t)
            K1,K2=map(lambda v:mp.mpf(str(v)),ode.sol(float(y)))
            contains(p['cutoff_integrals']['B_mass'],K1,'cutoff mass')
            contains(p['cutoff_integrals']['B_squared_mass'],K2,'cutoff squared mass')
            dr=R1*mp.exp(-1)*K1;dr2=R1*mp.exp(-1)*K2
            dtheta_factor=mp.sqrt(2)*R1**mp.mpf('1.5')*(mp.exp(t)-1)
            mixed_factor=mp.sqrt(2)*R1**mp.mpf('1.5')*mp.exp(-1)*K1
            increment=dict(z=(2*dr,4*dr),theta=(u1*dtheta_factor,u1z*dtheta_factor),
                theta_z=(2*u1*mixed_factor,4*(u1+mp.mpf('.5')*u1z)*mixed_factor),
                z_theta=(4*dr2-u1*u1*R1*t/2,16*dr2-u1*u1z*R1*t),
                p=(u1*u1*(1-mp.exp(-t))/2,u1*u1z*(1-mp.exp(-t))))
            expected={}
            for name,values in increment.items():
                expected[name]=[]
                for k,v in enumerate(values):
                    # Initial slope masses are interval-valued; compare the
                    # independent increment plus a contained inlet center.
                    lo,hi=endpoints(inlet['physical_moments'][name][k])
                    expected[name].append((lo+hi)/2+v)
                    contains(p['physical_moments'][name][k],expected[name][-1],(name,k))
            contains(p['Utheta'][0],u1*mp.exp(-t/2),'swirl')
            contains(p['Utheta'][1],u1z*mp.exp(-t/2),'swirl Z')
            q=mp.mpf(phase)
            if q==1:sig=mp.mpf(1);ds=mp.mpf(0)
            else:
                sig=1/(1+mp.exp(-1/(1-q)**2+1/q**2))
                ds=sig*(1-sig)*(2/q**3+2/(1-q)**3)
            bv=1-sig;by=-ds/(mp.mpf('1.1')*y)
            u=u1*mp.exp(-t/2);uz=u1z*mp.exp(-t/2)
            stress=evaluate_mp_stress(mp.log(R),mp.mpf('.5'),mp.mpf('.01'),
                Utheta=u,Uz=2*bv,Utheta_y=-u/2,Utheta_Z=uz,Uz_y=2*by,Uz_Z=4*bv,
                moments={k:v[0] for k,v in expected.items()},
                moments_Z={k:v[1] for k,v in expected.items()},
                P=mp.mpf('3.25')+expected['p'][0],P_Z=1+expected['p'][1],precision=100)
            for key in ('I_theta','I_z','S_theta','S_z','T_theta','T_z'):
                contains(p['stress'][key],stress[key],key)
            contains(p['Ur_value'],stress['U_r'],'radial recovery')
        start=field.evaluate_phase(0)
        for name in inlet['physical_moments']:
            for k in (0,1):
                if start['physical_moments'][name][k]._mpi_!=inlet['physical_moments'][name][k]._mpi_:
                    raise AssertionError(('inlet moment changed',name,k))
        end=field.evaluate_phase(1);tail=field.evaluate_buffer(11)
        for packet in (end,tail):
            for key in ('Uz','Uz_y'):
                if any(v._mpi_!=c.mpf(0)._mpi_ for v in packet[key].coefficients):
                    raise AssertionError('axial cutoff not exactly flat at endpoint')
        for k in (0,1):
            if end['physical_moments']['z'][k]._mpi_!=tail['physical_moments']['z'][k]._mpi_:
                raise AssertionError('retained axial mass reset in zero-Uz buffer')
        if endpoints(tail['physical_moments']['z'][0])[0]<=0:
            raise AssertionError('accumulated Mz disappeared with Uz')
        for Md,datum in (('.5','.5'),('1.1','.5')):
            try:IntervalOuterAxialTurnoff(c,z,c.mpf('.01'),inlet,Md,datum)
            except ValueError:pass
            else:raise AssertionError('invalid paper Md or mixed pressure schedule accepted')
        for bad_inlet in ({k:v for k,v in inlet.items() if k!='pressure_schedule_Md'},
                          dict(inlet,pressure_schedule_Md='.5')):
            try:IntervalOuterAxialTurnoff(c,z,c.mpf('.01'),bad_inlet,'1.1','1.1')
            except ValueError:pass
            else:raise AssertionError('unrecorded or mismatched inlet pressure Md accepted')
        coarse=cutoff_integrals(c,'1.1',1,64);fine=cutoff_integrals(c,'1.1',1,256)
        for key in coarse:
            a,b=endpoints(coarse[key]);d,e=endpoints(fine[key])
            if not a<=d<=e<=b:raise AssertionError('cutoff refinement not nested')
    here=Path(__file__).parent
    result=dict(passed=True,independent_fixture_containments=count,
        exact_C1_inlet_moments_preserved=True,terminal_Uz_and_Uz_y_exact_zero=True,
        retained_Mz_preserved_through_zero_axial_buffer=True,
        incompatible_pressure_Md_and_Md_le_one_rejected=True,
        new_production_pressure_schedule_regenerated=False,
        analytic_preheat_pressure_for_new_Md_certified=False,
        fixture_Md='1.1',fixture_scope='moderate synthetic field; numerical ODE is not a rigorous oracle',
        input_hashes={n:hashlib.sha256((here/n).read_bytes()).hexdigest() for n in
            (Path(__file__).name,'lei_ren_part1_paper_interval_outer_axial_turnoff.py',
            'lei_ren_part1_paper_interval_outer_slope_field.py','lei_ren_part1_paper_mp_stress.py')})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Axial turnoff checks passed:',count,'independent quantities; flat endpoint; retained mass; schedule guard',flush=True)
    return result

if __name__=='__main__':run()
