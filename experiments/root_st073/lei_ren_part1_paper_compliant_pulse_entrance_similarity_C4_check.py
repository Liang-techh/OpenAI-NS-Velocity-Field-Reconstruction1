"""Focused whole-entrance source and independent original startup integrals."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_entrance_similarity_C4 import (
    CompliantPulseEntranceSimilarityC4,DOMAIN,FALSE_FLAGS,PREFIX,source_precision)
from lei_ren_part1_paper_compliant_pulse_main_exit_similarity_C4 import uncapped_selected_pulse_source
from lei_ren_part1_paper_compliant_pulse_main_exit_fixture_integrals import run as fixture_integrals
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def independent_original_entrance_fixture():
    """Independent forward primitives at inlet/interior/main boundary."""
    with mp.workdps(75):
        c=MPIntervalContext();c.dps=100
        mu,delta,Rp,Pstar,D0,Pin,Xp=map(mp.mpf,('.2','.12','2.8','.45','.07','.12','.3'))
        p=1+2*mu;r=1-mu;bp=mp.mpf('.5')+mu
        fi=fixture_integrals();vals=fi['values']
        scalar=lambda key:mp.mpf(vals[key])
        array=lambda key:list(map(mp.mpf,vals[key]))
        A,D,weights,full_rows=[array(key) for key in ('A','D','weights','full_rows')]
        Kpulse,gram=scalar('Kpulse'),scalar('gram')
        C=lambda z:1/(1+z*z)
        incoming=(lambda z:mp.mpf('.1')*(z+z**3),lambda z:mp.mpf('.07')*(z+z**3))
        ein=lambda z:-mp.mpf('.015')+mp.mpf('.002')*(z*z+2*z**4+z**6)
        future=lambda z:mp.mpf('.8')+mp.mpf('.1')*z*z
        P0=lambda z:-mp.mpf('.35')*C(z)**2+mp.mpf('.02')*z
        def select(z):
            return uncapped_selected_pulse_source(mp,mu,Kpulse,[f(z) for f in incoming],
                ein(z),mu*mp.exp(-26)*future(z),A,D,full_rows,weights,mp.log(mu*D0))
        ap=lambda z:select(z)['ap']
        loss=lambda z:sum(v*v*mp.exp(-2*mu*center)*gram
            for v,center in zip(select(z)['controls'],(-3,-1)))
        def sigma(x):
            if x<=0:return mp.mpf(0)
            if x>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/x**2-1/(1-x)**2))
        gp_cache={}
        def gp(x):
            key=mp.nstr(x,85)
            if key not in gp_cache:
                gp_cache[key]=(mp.mpf(0) if x<=0 else mp.mpf('.01') if x==mp.mpf('.02')
                    else mp.quad(lambda b:sigma(50*b),[0,x/2,x],method='gauss-legendre'))
            return gp_cache[key]
        def original_kernel(x,i):
            if x==0:return mp.mpf(0)
            lam=mp.mpf('.5')-i*mu;k=lam/mu
            # Fubini exchanges the exact gp primitive with the forward kernel.
            value=mp.quad(lambda b:mp.exp(k*b)*sigma(50*b),[0,x/2,x])
            return (gp(x)-mp.exp(-k*x)*value)/lam
        def energy(x):
            if x==0:return mp.mpf(0)
            if x==mp.mpf('.02'):return scalar('entrance_energy')
            return mp.quad(lambda a:mp.exp(-2*a)*gp(a)**2,[0,x/2,x],method='gauss-legendre')
        tol=mp.mpf('1e-40');worst=mp.mpf(0);counts={};history_nonzero=radial_nonzero=axial_nonzero=0
        def compare(label,bound,expected,unit=1):
            nonlocal worst
            lo,hi=endpoints(bound)
            miss=max((lo-expected)/unit,(expected-hi)/unit,mp.mpf(0))
            worst=max(worst,miss)
            if miss>tol*max(1,abs(expected/unit)):
                raise ArithmeticError('Independent entrance source differs: '+label)
            counts[label]=counts.get(label,0)+1
        jet=lambda fn,z:IntervalTaylor(c,[c.mpf(v) for v in mp.taylor(fn,z,5)])
        for x,z in ((mp.mpf(0),mp.mpf('.31')),(mp.mpf('.01'),mp.mpf('-.4')),
                    (mp.mpf('.02'),mp.mpf('.5'))):
            y=x/mu;R=Rp*mp.exp(y);B=Pstar*mp.exp(-bp*y)
            g=gp(x);kernels=[original_kernel(x,i) for i in (1,2)];en=energy(x)
            def reference(zz):
                cc=C(zz);aa=ap(zz);ut=B*cc;uz=ut*aa*g
                m=[mp.exp(-(mp.mpf('.5')-i*mu)*y)*incoming[i-1](zz)+aa*kernels[i-1] for i in (1,2)]
                ev=mp.exp(2*x)*(ein(zz)+aa*aa*en/mu-(1-mp.exp(-2*x))/(4*mu))
                pressure=B*B*(-cc*cc/(2*p)+mp.exp(-p*(13-x)/mu)*(P0(zz)+cc*cc/(2*p)))
                transport=(1-delta)*zz*cc*m[0]+(1-zz*zz)*mp.diff(
                    lambda w:C(w)*(mp.exp(-(mp.mpf('.5')-mu)*y)*incoming[0](w)+ap(w)*kernels[0]),zz)
                ur=mp.sqrt(R/2)*B*(2*zz*cc*aa*g-transport)/(1-delta*zz*zz)
                moments=dict(theta=mp.sqrt(2)*R**mp.mpf('1.5')*ut*(1/r+(Xp-1/r)*mp.exp(-r*y)),
                    z=R*ut*m[0],theta_z=mp.sqrt(2)*R**mp.mpf('1.5')*ut*ut*m[1],
                    z_theta=R*ut*ut*ev,p=Pstar*Pstar*cc*cc*(Pin+(1-mp.exp(-p*y))/(2*p)))
                return dict(Ur=ur,Ut=ut,Uz=uz,P=pressure,moments=moments)
            def direct_stress(zz):
                d=reference(zz);ut=d['Ut']
                return evaluate_mp_stress(mp.log(R),zz,delta,Utheta=ut,Uz=d['Uz'],
                    Utheta_y=-bp*ut,Utheta_Z=B*mp.diff(C,zz),
                    Uz_y=B*C(zz)*ap(zz)*(mu*sigma(50*x)-bp*g),
                    Uz_Z=B*mp.diff(lambda w:C(w)*ap(w),zz)*g,
                    moments=d['moments'],moments_Z={key:mp.diff(lambda w:reference(w)['moments'][key],zz)
                        for key in d['moments']},P=d['P'],P_Z=mp.diff(lambda w:reference(w)['P'],zz),
                    precision=mp.mp.dps,radius_override=R,axial_override=zz)
            obj=object.__new__(CompliantPulseEntranceSimilarityC4)
            obj.ctx=c;obj.mu=c.mpf(mu);obj.delta=c.mpf(delta);obj.Xp=c.mpf(Xp);obj.cells=128
            obj.ap=jet(ap,z);obj.incoming=[jet(f,z) for f in incoming];obj.incoming_energy=jet(ein,z)
            obj.future=jet(future,z);obj.J0=jet(loss,z);obj.P0=jet(P0,z);obj.Pin=c.mpf(Pin)
            obj.logRp=c.mpf(mp.log(Rp));obj.logP=c.mpf(mp.log(Pstar));obj.logU=c.mpf(0)
            obj.finite=c.mpf(mp.log(D0)+1/mu)
            from lei_ren_part1_paper_compliant_pulse_entrance_similarity_C4 import compiled_entrance_exporter
            obj.exporter,_=compiled_entrance_exporter()
            packet=obj.entrance(c.mpf(z),c.mpf(x))
            def aggregate(group,rowkey,label,n):
                parts=packet[group][label]
                if group=='full_absolute_pressure_log_sectors':parts={label:parts}
                return sum((v[rowkey]['y0_Z'+str(n)]*c.exp(sum(v['exact_source_log_parts'].values(),c.mpf(0)))
                    for v in parts.values()),c.mpf(0))
            data=reference(z)
            for label,key in (('theta','T_theta'),('axial','T_z')):
                for n in range(4):
                    compare('original_full_stress_axial_C3',
                        aggregate('full_meridional_stress_log_sectors','full_stress_mixed3_coefficient_enclosures',label,n),
                        mp.diff(lambda zz:direct_stress(zz)[key],z,n),mp.sqrt(R/2)*B)
            for label,key,unit in (('radial','Ur',mp.sqrt(R/2)*B),('theta','Ut',B),('axial','Uz',B)):
                for n in range(5):
                    compare('original_velocity_axial_C4',
                        aggregate('full_velocity_log_sectors','full_velocity_mixed4_coefficient_enclosures',label,n),
                        mp.diff(lambda zz:reference(zz)[key],z,n),unit)
            for name in data['moments']:
                unit=R*B if name=='z' else mp.sqrt(2)*R**mp.mpf('1.5')*B if name=='theta' else (
                    mp.sqrt(2)*R**mp.mpf('1.5')*B*B if name=='theta_z' else R*B*B if name=='z_theta' else Pstar*Pstar)
                for n in range(5):
                    compare('all_five_raw_inlet_and_entrance_moments_C4',
                        aggregate('five_raw_cumulative_moment_log_sectors','full_moment_mixed4_coefficient_enclosures',name,n),
                        mp.diff(lambda zz:reference(zz)['moments'][name],z,n),unit)
            for n in range(5):
                total=sum(aggregate('full_absolute_pressure_log_sectors','full_pressure_mixed4_coefficient_enclosures',name,n)
                    for name in packet['full_absolute_pressure_log_sectors'])
                compare('same_absolute_pressure_C4',total,mp.diff(lambda zz:reference(zz)['P'],z,n),B*B)
                forward=packet['forward_full_energy_raw_moment']
                bound=forward['full_moment_mixed4_coefficient_enclosures']['y0_Z'+str(n)]*c.exp(
                    sum(forward['exact_source_log_parts'].values(),c.mpf(0)))
                compare('directly_anchored_forward_energy_C4',bound,
                    mp.diff(lambda zz:reference(zz)['moments']['z_theta'],z,n),R*B*B)
            if x==0 and all(abs(data['moments'][k])>mp.mpf('1e-10') for k in ('z','theta_z','z_theta')):
                history_nonzero+=1
            if abs(data['Ur'])>mp.mpf('1e-10'):radial_nonzero+=1
            if abs(data['Uz'])>mp.mpf('1e-10'):axial_nonzero+=1
            print('Independent original entrance fixture xi='+str(x)+' checked',flush=True)
        if history_nonzero!=1 or radial_nonzero!=3 or axial_nonzero!=2:
            raise ArithmeticError('Required nonzero entrance histories were lost')
        return dict(comparisons=counts,maximum_normalized_positive_enclosure_miss=mp.nstr(worst,30),
            normalized_tolerance=str(tol),nonzero_inlet_history_fixtures=history_nonzero,
            nonzero_radial_velocity_fixtures=radial_nonzero,nonzero_axial_input_fixtures=axial_nonzero,
            independent_original_sigma_primitive_and_Fubini_forward_kernel_used=True,
            forward_incoming_energy_independent_of_companion_backward_energy=True,
            actual_original_full_stress_evaluator_used=True,fixture_is_not_a_cone_or_corrected_residual_test=True)


@source_precision
def run():
    companion=CompliantPulseEntranceSimilarityC4();actual=companion.report()
    name=PREFIX+'pulse_entrance_similarity_C4.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
    for path,digest in record['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Entrance source changed: '+path)
    if record!=encode(pack(actual)):raise ValueError('Current whole entrance companion differs')
    if record['domain']!=DOMAIN:raise ValueError('Original entrance domain changed')
    counts={}
    for key in ('whole_original_entrance','original_inlet','early_ordinary_y_chart','common_entrance_main'):
        packet=actual[key]
        for group,rowkey,order in (('full_meridional_stress_log_sectors','full_stress_mixed3_coefficient_enclosures',3),
            ('full_velocity_log_sectors','full_velocity_mixed4_coefficient_enclosures',4),
            ('five_raw_cumulative_moment_log_sectors','full_moment_mixed4_coefficient_enclosures',4),
            ('full_absolute_pressure_log_sectors','full_pressure_mixed4_coefficient_enclosures',4)):
            groups=({label:{label:part} for label,part in packet[group].items()} if group=='full_absolute_pressure_log_sectors'
                else packet[group])
            required={'y'+str(j)+'_Z'+str(n) for j in range(order+1) for n in range(order+1-j)}
            for parts in groups.values():
                for sector in parts.values():
                    if set(sector[rowkey])!=required:raise ValueError('Entrance mixed derivative grid incomplete')
                    for value in sector[rowkey].values():
                        if not all(mp.isfinite(v) for v in endpoints(value)):raise ValueError('Nonfinite entrance source bound')
                        counts[group]=counts.get(group,0)+1
        for flag in FALSE_FLAGS:
            if packet[flag] or record[flag]:raise ValueError('Entrance scope overclaimed: '+flag)
    # The exact common function is proved by the unchanged exporter AST.
    # At the common point, its current directed outputs must also match.
    main_left=companion.main_exit([-1,1],'.02')
    for key in ('full_meridional_stress_log_sectors','full_velocity_log_sectors',
        'five_raw_cumulative_moment_log_sectors','full_absolute_pressure_log_sectors',
        'main_exit_source_rows','exact_source_logs','original_partial_linear_kernel_bounds'):
        if encode(pack(main_left[key]))!=encode(pack(actual['common_entrance_main'][key])):
            raise ValueError('Same-source xi=.02 exporter output differs: '+key)
    source=actual['original_inlet']['main_exit_source_rows']
    for key in ('Bh','ml','nl'):
        if any(endpoints(v)!=(0,0) for row in source[key] for v in row.coefficients):
            raise ValueError('Only local inlet input/moment jets must be zero')
    for j,v in enumerate(actual['original_inlet']['forward_full_energy_source_rows'][0].coefficients):
        if endpoints(v)!=endpoints(companion.incoming_energy[j]):
            raise ValueError('Forward energy inlet is not the actual incoming datum')
    for key,rowkey,order in (('forward_full_energy_raw_moment','full_moment_mixed4_coefficient_enclosures',4),
        ('forward_full_energy_only_axial_stress','full_stress_mixed3_coefficient_enclosures',3)):
        for packetkey in ('whole_original_entrance','original_inlet','early_ordinary_y_chart','common_entrance_main'):
            grid=actual[packetkey][key][rowkey]
            required={'y'+str(j)+'_Z'+str(n) for j in range(order+1) for n in range(order+1-j)}
            if set(grid)!=required:raise ValueError('Forward-energy derivative grid incomplete')
            for value in grid.values():
                if not all(mp.isfinite(v) for v in endpoints(value)):raise ValueError('Nonfinite forward-energy row')
                counts[key]=counts.get(key,0)+1
    oracle=independent_original_entrance_fixture()
    hashes=dict(actual['input_hashes']);hashes[name]=hashlib.sha256(raw).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for suffix in ('.py','.json'):
        path=PREFIX+'pulse_main_exit_fixture_integrals'+suffix
        hashes[path]=hashlib.sha256((HERE/path).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=companion.family,implicit_source_sha256=companion.source,
        domain=DOMAIN,input_hashes=hashes,source_function_identities_checked=len(actual['source_function_proof']['identities']),
        finite_signed_rows_checked=counts,independent_original_entrance_fixture=oracle,
        current_whole_entrance_source_report_recomputed=True,
        original_inlet_and_entrance_main_similarity_functional_joins_verified=True,
        nonzero_incoming_moments_radial_velocity_and_energy_preserved=True,
        actual_original_whole_entrance_similarity_companion_constructed=True,
        source_caps_used_as_defining_field_values=False,**{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('PASS whole original entrance five-moment/pressure/velocity/stress source; physical/cone pending',flush=True)
    return result


if __name__=='__main__':run()
