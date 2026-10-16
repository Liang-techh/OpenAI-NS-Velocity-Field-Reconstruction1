"""Actual original selected-N Rh-reference numerical/physical continuation."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[1];sys.path[:0]=[str(ROOT/'src'),str(ROOT/'experiments/root_st073')]
DIRECTORY=ROOT/'experiments/root_st073'
from openai_ns_reconstruction.original_core_numeric import PREFIX,_bounds
from openai_ns_reconstruction.original_outer_reference_velocity import OriginalOuterReferenceVelocity
from run_original_Rm_patch_velocity import inlet_physical_seam
from run_original_second_micro_numeric import overlap

NAME=PREFIX+'compliant_current_original_outer_reference_numeric.json'


def independent_own_ODEs(f,packet):
    """Integrate all five own density equations from the actual Rh inlet."""
    # Use the original amplitude context for this independent normalization.
    # Re-exponentiating its logPstar after rounding to the lower field context
    # introduces a relative error of order1e-69 into both inlet and local source.
    c=f.amplitude.ctx;Z=packet['Z'];end=c.mpf(packet['native_offset']);z=c.mpf(Z)
    inlet=f.numeric._inlets[Z];P=c.exp(f.amplitude.lp)
    def density(q):
        E=c.exp(q/10)/(1+z*z);V=4*z/P
        # Normalize the tiny m/k amplitudes before ordinary quadrature so
        # its absolute tolerance cannot discard their local source.
        return dict(m=4*z,h=E,k=E*(4*z),e=V*V-E*E/2,p=E*E/2)
    count=0
    for name,rate in dict(m=1,h=c.mpf('1.5'),k=c.mpf('1.5'),e=1,p=0).items():
        scale=P if name in ('m','k') else c.mpf(1)
        inlet_row=inlet['generic_modified_fields'][name][0]
        initial=c.mpf(inlet_row['approximation'])*scale
        inlet_lo,inlet_hi=_bounds(c,f._restore(inlet_row)*scale)
        inlet_radius=max(abs(inlet_lo-initial),abs(inlet_hi-initial))
        value=initial*c.exp(-rate*(end+5))+c.quad(lambda q:c.exp(-rate*(end-q))*density(q)[name],[-5,end])
        truth=f._restore(packet['generic_leading_fields'][name][0])*scale;lo,hi=_bounds(c,truth)
        # The accepted inlet point is a serialized proxy. Its known directed
        # uncertainty must follow the same homogeneous kernel; normalizing an
        # exp(-logPstar) value can expose lost digits of the huge logPstar.
        tol=inlet_radius*c.exp(-rate*(end+5))+c.mpf('1e-70')*max(abs(value),c.mpf('1e-100'))
        assert lo-tol<=value<=hi+tol,(Z,end,name,c.nstr(value,35),c.nstr(lo,35),c.nstr(hi,35),
            c.nstr(value-lo,20),c.nstr(tol,10));count+=1
    return count


def run():
    started=time.perf_counter();f=OriginalOuterReferenceVelocity();samples=[]
    p=DIRECTORY/(PREFIX+'compliant_current_original_Rm_patch_velocity.json')
    old=json.loads(p.read_text(encoding='utf-8'))
    assert old['checks']['passed'] and old['source_family']==f.family
    f.hashes[p.name]=hashlib.sha256(p.read_bytes()).hexdigest()
    for name,digest in old['implementation_hashes'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
        f.code_hashes[name]=digest
    fieldseams=physical=history=measure=pressure=ode=0
    for Z in ('-.3','0','.3','.17','.23'):
        for offset in ('-5','-4.5','-2.5','0'):
            packet=f.evaluate(Z,offset,log_tau='-2',azimuth='.25');samples.append(packet)
            assert packet['selected_integer']==f.numeric._inlets[Z]['selected_integer'];history+=1
            assert packet['positive_order_n1_own_history']==f.numeric._inlets[Z]['positive_order_n1_own_history'];history+=1
            assert packet['independent_P0_over_Pstar_squared']==f.numeric._inlets[Z]['independent_P0_over_Pstar_squared'];pressure+=1
            for k in range(3):
                # Same independently normalized datum, added exactly once.
                p0=f._restore(packet['independent_P0_over_Pstar_squared'][k])
                pp=f._restore(packet['generic_leading_fields']['p'][k]);P=f._restore(packet['generic_leading_fields']['pressure'][k])
                lo,hi=_bounds(f.ctx,P-p0-pp);assert lo<=0<=hi;pressure+=1
            width=sum((Fraction(*s['exact_end'])-Fraction(*s['exact_left']) for s in packet['physical_transport_steps']),Fraction(0))
            assert width==Fraction(offset)+5;measure+=1
            for step in packet['physical_transport_steps']:
                weight=step['own_weights']['p']
                assert weight['decay']['lower_exact_mpf_tuple']==weight['decay']['upper_exact_mpf_tuple']==[0,1,0,1];measure+=1
            if offset=='-5':
                inlet=f.numeric._inlets[Z]
                for group in ('generic_modified_fields','physical_similarity','finite_N_history_errors'):
                    for name,rows in packet[group].items():
                        for a,b in zip(rows,inlet[group][name]):overlap(a,b);fieldseams+=1
                previous=next(a for a in old['samples'] if a['Z']==Z and a['native_x']=='Rh' and a['physical']['time']['kind']=='log_tau')
                physical+=inlet_physical_seam(f,previous['physical'],packet['physical'])
            if Z in ('-.3','0','.3') and offset in ('-4.5','0'):ode+=independent_own_ODEs(f,packet)
        print('Actual original Rh-to-Rref values and histories: Z='+Z,flush=True)
    scaled=0
    for Z in ('-.3','0','.3'):
        p=f.evaluate(Z,'0',log_tau=None,scaled_time='1e18',azimuth='1.1');samples.append(p)
        assert p['physical']['time']['delta_cancelled_before_time_comparison'];scaled+=1
        normal=next(a for a in samples if a['Z']==Z and a['native_offset']=='0')
        for name,terms in p['physical']['cylindrical'].items():
            for a,b in zip(terms,normal['physical']['cylindrical'][name]):overlap(a['coefficient'],b['coefficient']);scaled+=1
    boundary=0
    for s,expected in [('-4.'+'0'*139+'1',1),('-3.'+'9'*140,2)]:
        p=f.numeric.evaluate('.3',s)
        assert len(p['physical_transport_steps'])==expected;boundary+=1
    rejected=0
    for offset,kwargs in (('-5.1',{}),('.1',{}),('nan',{}),('-2',dict(log_tau='1')),('-2',dict(azimuth='nan'))):
        try:f.evaluate('.3',offset,**kwargs)
        except ValueError:rejected+=1
        else:raise AssertionError('Invalid original reference request admitted')
    code=dict(f.code_hashes);code[Path(__file__).relative_to(ROOT).as_posix()]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    return dict(source_family=f.family,samples=samples,input_hashes=f.hashes,implementation_hashes=code,
        checks=dict(passed=True,actual_original_reference_physical_packets=len(samples),Rh_numeric_seam_coefficients=fieldseams,
            Rh_physical_sum_seams=physical,N_and_nonzero_n1_history_checks=history,
            physical_measure_and_pressure_memory_checks=measure,P0_unit_and_once_only_checks=pressure,
            independent_actual_own_density_ODE_values=ode,scaled_time_checks=scaled,exact_close_boundary_checks=boundary,rejected_requests=rejected),
        execution_seconds=round(time.perf_counter()-started,3),
        limits=dict(through_original_Rref=True,beyond_Rref=False,physical_Cartesian_values=True,native_coordinates_required=True,
            finite_N_values_are_point_plus_directed_error=True,resolved_signed_finite_N_centers=False,
            higher_Z_or_radial_fast_phase_derivatives=False,full_n1_pressure=False,reduced_m13_m15_integrals_continued=False,
            physical_gradients=False,direct_physical_inverse=False,
            global_n1_high_order_joins=False,global_N_admission=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))


def inspect(path):
    result=json.loads(path.read_text(encoding='utf-8'));assert result['checks']['passed']
    for name,digest in result['input_hashes'].items():assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
    for name,digest in result['implementation_hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-report',action='store_true');parser.add_argument('--inspect-report',action='store_true')
    args=parser.parse_args();path=DIRECTORY/NAME
    result=inspect(path) if args.inspect_report else run()
    if args.write_report:path.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=result['checks'],execution_seconds=result['execution_seconds']),indent=2),flush=True)
