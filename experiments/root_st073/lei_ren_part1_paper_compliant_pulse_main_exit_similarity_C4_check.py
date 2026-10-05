"""Independent original-integral and full-stress checks for main/exit."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_main_exit_similarity_C4 import (
    CompliantPulseMainExitSimilarityC4,uncapped_selected_pulse_source,
    DOMAIN,FALSE_FLAGS,PREFIX,source_precision)
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_compliant_pulse_main_exit_fixture_integrals import run as load_fixture_integrals
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def independent_main_exit_fixture():
    """Actual gp/beta integrals; independent forward energy and original stress."""
    with mp.workdps(100):
        c=MPIntervalContext();c.dps=130
        mu,delta,Rp,Pstar,D0,Pin,Xp,ell=map(mp.mpf,('.2','.12','2.8','.45','.07','.12','.3','.15'))
        p=1+2*mu;r=1-mu;bp=mp.mpf('.5')+mu;logG=mp.log(mu*D0)
        sigma=lambda x:mp.mpf(0) if x<=0 else mp.mpf(1) if x>=1 else 1/(1+mp.exp(1/x**2-1/(1-x)**2))
        def gp(x):
            if x<=0 or x>=11:return mp.mpf(0)
            primitive=x-mp.mpf('.01') if x>=mp.mpf('.02') else mp.quad(
                lambda v:sigma(50*v),[0,x/2,x],method='gauss-legendre')
            return primitive*sigma(11-x)
        integral_fixture=load_fixture_integrals()
        values=integral_fixture["values"]
        scalar=lambda key:mp.mpf(values[key])
        array=lambda key:list(map(mp.mpf,values[key]))
        normal=scalar("normal");W=array("W");gram=scalar("gram")
        A=array("A");D=array("D");weights=array("weights")
        row2=[a+mu*d for a,d in zip(A,D)]
        entrance_linear=array("entrance_linear");entrance_energy=scalar("entrance_energy")
        exit_energy=scalar("exit_energy");Kpulse=scalar("Kpulse");full_rows=array("full_rows")
        def linear_integral(x,i):
            k=(mp.mpf(".5")-i*mu)/mu
            prim=lambda a:mp.exp(k*a)*((a-mp.mpf(".01"))/k-1/k**2)
            total=entrance_linear[i-1]+prim(min(x,mp.mpf(10)))-prim(mp.mpf(".02"))
            if x>10:total+=mp.quad(lambda a:mp.exp(k*a)*gp(a),[10,(10+x)/2,x])
            return total/mu
        energy_primitive=lambda a:-mp.exp(-2*a)*((a-mp.mpf(".01"))**2/2+(a-mp.mpf(".01"))/2+mp.mpf(".25"))
        def partial_energy(x):
            result=entrance_energy+energy_primitive(min(x,mp.mpf(10)))-energy_primitive(mp.mpf('.02'))
            if x>10:result+=mp.quad(lambda a:mp.exp(-2*a)*gp(a)**2,[10,(10+x)/2,x])
            return result
        C=lambda z:1/(1+z*z)
        incoming=(lambda z:mp.mpf('.1')*(z+z**3),lambda z:mp.mpf('.07')*(z+z**3))
        energy_in=lambda z:-mp.mpf('.015')+mp.mpf('.002')*(z*z+2*z**4+z**6)
        future=lambda z:mp.mpf('.8')+mp.mpf('.1')*z*z
        P0=lambda z:-mp.mpf('.35')*C(z)**2+mp.mpf('.02')*z
        def select(z):
            return uncapped_selected_pulse_source(mp,mu,Kpulse,[fn(z) for fn in incoming],energy_in(z),
                mu*mp.exp(-26)*future(z),A,D,full_rows,weights,logG)
        ap=lambda z:select(z)['ap']
        controls=[lambda z,j=j:select(z)['controls'][j] for j in range(2)]
        J0=lambda z:sum(fn(z)**2*mp.exp(-2*mu*center)*gram for fn,center in zip(controls,(-3,-1)))
        tol=mp.mpf('1e-50');worst=mp.mpf(0);counts={};nonzero_shear=0
        def compare(label,bound,value,unit=1):
            nonlocal worst
            lo,hi=endpoints(bound);miss=max((lo-value)/unit,(value-hi)/unit,mp.mpf(0));worst=max(worst,miss)
            if miss>tol*max(1,abs(value/unit)):raise ArithmeticError('Independent main/exit fixture differs: '+label)
            counts[label]=counts.get(label,0)+1
        jet=lambda fn,z:IntervalTaylor(c,[c.mpf(v) for v in mp.taylor(fn,z,5)])
        zj=mp.mpf('.31')
        selected=uncapped_selected_pulse_source(c,c.mpf(mu),c.mpf(Kpulse),
            [jet(fn,zj) for fn in incoming],jet(energy_in,zj),
            jet(lambda z:mu*mp.exp(-26)*future(z),zj),
            list(map(c.mpf,A)),list(map(c.mpf,D)),list(map(c.mpf,full_rows)),list(map(c.mpf,weights)),c.mpf(logG))
        for j in range(6):
            compare('uncapped_selected_C5',selected['ap'][j],mp.diff(ap,zj,j)/math.factorial(j))
        root=select(zj)
        compare('selected_quadratic',c.mpf(0),root['A2']*root['ap']**2+root['A1']*root['ap']+root['A0'])
        for i in (1,2):
            actual=sum((row2 if i==2 else A)[j]*root['controls'][j] for j in range(2))
            expected=-mu*(incoming[i-1](zj)*mp.exp(-13*(mp.mpf('.5')-i*mu)/mu-logG)+root['ap']*full_rows[i-1])
            compare('selected_linear_rows',c.mpf(actual),expected)
        try:
            uncapped_selected_pulse_source(mp,mu,Kpulse,[mp.mpf(0)]*2,mp.mpf(10),mp.mpf(0),
                A,D,full_rows,weights,logG)
        except ArithmeticError:
            counts['invalid_positive_branch_rejected']=1
        else:
            raise ArithmeticError('Invalid selected positive branch accepted')
        B=lambda y:Pstar*mp.exp(-bp*y)
        R=lambda y:Rp*mp.exp(y)
        def reference(x,z):
            y=x/mu;cc=C(z);aa=ap(z);ut=B(y)*cc;gz=gp(x)
            mi=[mp.exp(-(mp.mpf('.5')-i*mu)*y)*(incoming[i-1](z)+aa*linear_integral(x,i)) for i in (1,2)]
            ev=mp.exp(2*x)*(energy_in(z)+aa*aa*partial_energy(x)/mu-(1-mp.exp(-2*x))/(4*mu))
            press=-cc*cc/(2*p)+mp.exp(-p*(13-x)/mu)*(P0(z)+cc*cc/(2*p))
            uz=ut*aa*gz
            transport=(1-delta)*z*cc*mi[0]+(1-z*z)*mp.diff(
                lambda zz:C(zz)*mp.exp(-(mp.mpf('.5')-mu)*y)*(incoming[0](zz)+ap(zz)*linear_integral(x,1)),z)
            ur=mp.sqrt(R(y)/2)*B(y)*(2*z*cc*aa*gz-transport)/(1-delta*z*z)
            moments=dict(theta=mp.sqrt(2)*R(y)**mp.mpf('1.5')*ut*(1/r+(Xp-1/r)*mp.exp(-r*y)),
                z=R(y)*ut*mi[0],theta_z=mp.sqrt(2)*R(y)**mp.mpf('1.5')*ut**2*mi[1],
                z_theta=R(y)*ut**2*ev,p=Pstar**2*cc**2*(Pin+(1-mp.exp(-p*y))/(2*p)))
            return dict(Ur=ur,Ut=ut,Uz=uz,P=B(y)**2*press,moments=moments)
        # Freeze each scalar integral only for axial differentiation. This
        # is exact since xi is held fixed and every integral is Z-independent.
        cache={}
        original_linear=linear_integral;original_energy=partial_energy
        def linear_cached(x,i):
            key=(mp.nstr(x,110),i)
            if key not in cache:cache[key]=original_linear(x,i)
            return cache[key]
        energy_cache={}
        def energy_cached(x):
            key=mp.nstr(x,110)
            if key not in energy_cache:energy_cache[key]=original_energy(x)
            return energy_cache[key]
        linear_integral=linear_cached;partial_energy=energy_cached
        def direct_stress(x,z):
            data=reference(x,z);y=x/mu;ut=data['Ut'];aa=ap(z)
            uz_y=B(y)*C(z)*aa*(mu*mp.diff(gp,x)-bp*gp(x))
            return evaluate_mp_stress(mp.log(R(y)),z,delta,Utheta=ut,Uz=data['Uz'],Utheta_y=-bp*ut,
                Utheta_Z=B(y)*mp.diff(C,z),Uz_y=uz_y,Uz_Z=B(y)*mp.diff(lambda zz:C(zz)*ap(zz),z)*gp(x),
                moments=data['moments'],moments_Z={name:mp.diff(lambda zz:reference(x,zz)['moments'][name],z)
                    for name in data['moments']},P=data['P'],P_Z=mp.diff(lambda zz:reference(x,zz)['P'],z),
                precision=mp.mp.dps,radius_override=R(y),axial_override=z)
        for x,z in ((mp.mpf(1),zj),(mp.mpf('10.4'),mp.mpf('-.4')),(mp.mpf(11),mp.mpf('.5'))):
            fixture=object.__new__(CompliantPulseMainExitSimilarityC4)
            fixture.ctx=c;fixture.mu=c.mpf(mu);fixture.delta=c.mpf(delta);fixture.Xp=c.mpf(Xp);fixture.cells=128
            fixture.ap=jet(ap,z);fixture.incoming=[jet(fn,z) for fn in incoming]
            fixture.future=jet(future,z);fixture.J0=jet(J0,z);fixture.P0=jet(P0,z)
            fixture.Pin=c.mpf(Pin);fixture.logP=c.mpf(mp.log(Pstar));fixture.logU=c.mpf(0)
            fixture.logRp=c.mpf(mp.log(Rp));fixture.finite=c.mpf(mp.log(D0)+1/mu)
            packet=fixture.main_exit(c.mpf(z),c.mpf(x));data=reference(x,z);y=x/mu
            def aggregate(group,rowkey,label,j,n):
                sectors=packet[group][label]
                if group=='full_absolute_pressure_log_sectors':sectors={label:sectors}
                total=c.mpf(0)
                for sector in sectors.values():
                    total+=sector[rowkey]['y'+str(j)+'_Z'+str(n)]*c.exp(sum(sector['exact_source_log_parts'].values(),c.mpf(0)))
                return total
            for label,key in (('theta','T_theta'),('axial','T_z')):
                for n in range(4):
                    bound=aggregate('full_meridional_stress_log_sectors','full_stress_mixed3_coefficient_enclosures',label,0,n)
                    expected=mp.diff(lambda zz:direct_stress(x,zz)[key],z,n)
                    compare('original_full_stress_axial_C3',bound,expected,mp.sqrt(R(y)/2)*B(y))
            for label,key,unit in (('radial','Ur',mp.sqrt(R(y)/2)*B(y)),('theta','Ut',B(y)),('axial','Uz',B(y))):
                for n in range(5):
                    compare('velocity_axial_C4',aggregate('full_velocity_log_sectors','full_velocity_mixed4_coefficient_enclosures',label,0,n),
                            mp.diff(lambda zz:reference(x,zz)[key],z,n),unit)
            for name in data['moments']:
                unit=R(y)*B(y) if name=='z' else mp.sqrt(2)*R(y)**mp.mpf('1.5')*B(y) if name=='theta' else (
                    mp.sqrt(2)*R(y)**mp.mpf('1.5')*B(y)**2 if name=='theta_z' else R(y)*B(y)**2 if name=='z_theta' else Pstar**2)
                for n in range(5):
                    compare('all_five_raw_moments_axial_C4',aggregate('five_raw_cumulative_moment_log_sectors',
                        'full_moment_mixed4_coefficient_enclosures',name,0,n),
                        mp.diff(lambda zz:reference(x,zz)['moments'][name],z,n),unit)
            for n in range(5):
                bound=sum(aggregate('full_absolute_pressure_log_sectors','full_pressure_mixed4_coefficient_enclosures',label,0,n)
                    for label in packet['full_absolute_pressure_log_sectors'])
                compare('same_absolute_pressure_axial_C4',bound,mp.diff(lambda zz:reference(x,zz)['P'],z,n),B(y)**2)
            if x<11 and abs(mu*mp.diff(gp,x)-bp*gp(x))>0:nonzero_shear+=1
            print('Independent original full-stress fixture xi='+str(x)+' compared',flush=True)
        return dict(comparisons=counts,maximum_normalized_positive_enclosure_miss=mp.nstr(worst,30),
            normalized_tolerance=str(tol),nonzero_axial_shear_fixtures=nonzero_shear,
            original_integral_fixture_generator_sha256=integral_fixture['generator_sha256'],
            original_startup_exit_gp_and_full_beta_integrals_used=True,
            exact_uncapped_selected_quadratic_and_C5_recurrence_checked=True,
            forward_energy_independent_of_companion_backward_energy=True,
            original_full_meridional_stress_evaluator_used=True,
            zero_output_cannot_pass_by_tiny_physical_amplitude=True,
            source_cone_and_corrected_residual_not_tested=True)


@source_precision
def run():
    companion=CompliantPulseMainExitSimilarityC4();actual=companion.report()
    name=PREFIX+'pulse_main_exit_similarity_C4.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
    for path,digest in record['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Main/exit source changed: '+path)
    if record!=encode(pack(actual)):raise ValueError('Current main/exit report differs')
    if record['domain']!=DOMAIN:raise ValueError('Original main/exit domain changed')
    identities=companion.mainproof['identities']
    if not all(identities.values()):raise ValueError('Source identity failed')
    counts={}
    for key in ('whole_original_main','whole_original_exit','original_main_left','common_main_exit','original_exit_gap'):
        packet=actual[key]
        for group,rowkey in (('full_meridional_stress_log_sectors','full_stress_mixed3_coefficient_enclosures'),
            ('full_velocity_log_sectors','full_velocity_mixed4_coefficient_enclosures'),
            ('five_raw_cumulative_moment_log_sectors','full_moment_mixed4_coefficient_enclosures')):
            for sectors in packet[group].values():
                for sector in sectors.values():
                    order=3 if group=='full_meridional_stress_log_sectors' else 4
                    required={'y'+str(j)+'_Z'+str(n) for j in range(order+1) for n in range(order+1-j)}
                    if set(sector[rowkey])!=required:raise ValueError('Incomplete required mixed source grid: '+group)
                    for value in sector[rowkey].values():
                        if not all(mp.isfinite(v) for v in endpoints(value)):raise ArithmeticError('Nonfinite original source row')
                        counts[group]=counts.get(group,0)+1
        for sector in packet['full_absolute_pressure_log_sectors'].values():
            required={'y'+str(j)+'_Z'+str(n) for j in range(5) for n in range(5-j)}
            if set(sector['full_pressure_mixed4_coefficient_enclosures'])!=required:
                raise ValueError('Incomplete required absolute-pressure mixed4 grid')
            for value in sector['full_pressure_mixed4_coefficient_enclosures'].values():
                if not all(mp.isfinite(v) for v in endpoints(value)):raise ArithmeticError('Nonfinite pressure row')
                counts['full_absolute_pressure']=counts.get('full_absolute_pressure',0)+1
        for flag in ('selected_source_is_uncapped_positive_quadratic','coefficients_are_only_enclosures_of_defining_sources',
            'original_full_gp_and_axial_shear_retained','original_forward_linear_histories_not_replaced_by_end_only_histories',
            'same_absolute_pressure_getter_and_raw_Mp_kept_separate'):
            if not packet[flag]:raise ValueError('Main/exit source scope missing')
        for flag in FALSE_FLAGS:
            if packet[flag] or record[flag]:raise ValueError('Main/exit scope overclaimed: '+flag)
    fixture=independent_main_exit_fixture()
    hashes=dict(actual['input_hashes']);hashes[name]=hashlib.sha256(raw).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for suffix in ('.py','.json'):
        fixture_name=PREFIX+'pulse_main_exit_fixture_integrals'+suffix
        hashes[fixture_name]=hashlib.sha256((HERE/fixture_name).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=companion.family,implicit_source_sha256=companion.source,
        domain=DOMAIN,source_function_identities_checked=len(identities),finite_signed_rows_checked=counts,
        current_original_main_exit_report_recomputed=True,independent_original_integral_fixture=fixture,
        actual_original_whole_main_exit_similarity_companion_constructed=True,
        main_exit_and_exit_gap_similarity_source_functional_joins_consumed=True,
        source_caps_used_as_defining_field_values=False,input_hashes=hashes,**{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('PASS full original main/exit similarity companion; physical/cone/global/recursion pending',flush=True)
    return result


if __name__=='__main__':run()
