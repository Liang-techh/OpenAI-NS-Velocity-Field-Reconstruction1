"""Compute original logarithmic amplitude, resolved peak and leading velocity terms."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from openai_ns_reconstruction.original_core_log_amplitude import OriginalCoreLogAmplitude
from openai_ns_reconstruction.original_core_numeric import DIRECTORY, _box, _bounds

NAME = "lei_ren_part1_paper_compliant_current_original_core_log_amplitude.json"


def real_integral(solver,Z):
    """Independent nominal sinh-coordinate quadrature; not an error certificate."""
    c = solver.ctx
    anchor,hp,sigma = solver.anchor,solver._hp(solver.anchor),solver.sigma
    end = c.asinh((c.mpf(Z)-anchor)*hp/sigma)
    sign = 1 if end >= 0 else -1
    points = [c.mpf(0)]
    for value in (1,2,4,8,16,32,64,128):
        q = c.mpf(sign*value)
        if abs(q) >= abs(end):
            break
        points.append(q)
    points.append(end)
    def integrand(q):
        z = anchor+sigma*c.sinh(q)/hp
        H = solver._H(z)
        return solver._L(z)*H/(H*H+sigma*sigma)*sigma*c.cosh(q)/hp
    return c.quad(integrand,points)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-report",action="store_true")
    parser.add_argument("--inspect-report",action="store_true")
    args=parser.parse_args()
    path=DIRECTORY/NAME
    started=time.monotonic()
    if args.inspect_report:
        if args.write_report:
            parser.error("Inspect separately from construction")
        report=json.loads(path.read_text(encoding="utf-8"))
        for name,digest in report['implementation_hashes'].items():
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
        for name,digest in report['input_hashes'].items():
            assert hashlib.sha256((DIRECTORY/name).read_bytes()).hexdigest()==digest,name
        assert report['checks']['passed']
        print(json.dumps(dict(saved_report_restored=True,full_owner_restored=False,
            checks=report['checks'],execution_seconds=round(time.monotonic()-started,3)),indent=2))
        return
    levels=[]
    for precision in (80,110):
        solver=OriginalCoreLogAmplitude(precision=precision)
        global_values=[solver.global_amplitude(Z) for Z in ('-1','-0.3','0','0.3','1','anchor')]
        peaks=[solver.resolved_peak(xi) for xi in ('-4','-2','-1','-0.5','0','0.5','1','2','4')]
        levels.append(dict(precision=precision,global_amplitude=global_values,resolved_peak=peaks))
        print(f"Completed original log-amplitude and peak at {precision} digits",flush=True)
    c=mp.mp.clone();c.dps=160
    def consistent(left,right):
        delta=abs(c.mpf(left['approximation'])-c.mpf(right['approximation']))
        return delta<=c.mpf(left['absolute_error_upper'])+c.mpf(right['absolute_error_upper'])
    comparisons=0
    for left,right in zip(levels[0]['global_amplitude'],levels[1]['global_amplitude']):
        for key in ('G','log_F0_relative_to_anchor','anchor'):
            assert consistent(left[key],right[key]),(left['Z'],key)
            comparisons+=1
    for left,right in zip(levels[0]['resolved_peak'],levels[1]['resolved_peak']):
        for key in ('scaled_G_Lambda_times_G','F0_relative_to_exact_anchor','width_b'):
            assert consistent(left[key],right[key]),(left['xi'],key)
            comparisons+=1
    integrals=[]
    for Z in ('0','0.3'):
        quadrature=real_integral(solver,Z)
        point=c.mpf(solver.global_amplitude(Z)['G']['approximation'])
        difference=abs(c.mpf(solver.ctx.nstr(quadrature,110))-point)
        assert difference<c.mpf('1e-85')
        integrals.append(dict(Z=Z,nominal_integral=solver.ctx.nstr(quadrature,110),
            difference=c.nstr(difference,40),quadrature_error_certified=False,
            source_error_provenance="archived certified poles and guarded log branches"))
    # Compare local calculations to the existing same-source root chart.
    peak_path=DIRECTORY/'lei_ren_part1_paper_compliant_root_centered_peak.json'
    peak_archive=json.loads(peak_path.read_text(encoding='utf-8'))
    for key,value in solver.family.items():assert peak_archive[key]==value
    aliases={'-0.5':'-.5','0.5':'.5'}
    old_peak_pairs=0
    for packet in levels[-1]['resolved_peak']:
        old=peak_archive['peak_packets'][aliases.get(packet['xi'],packet['xi'])]
        for new_key,old_key in (('scaled_G_Lambda_times_G','actual_scaled_G'),
                                ('F0_relative_to_exact_anchor','normalized_F0_relative_to_shared_anchor')):
            new=packet[new_key]
            box=_box(solver.iv,old[old_key])
            point=solver.iv.mpf(new['approximation'])
            error=solver.iv.mpf(new['absolute_error_upper'])
            low,high=_bounds(c,point-error);nlow,nhigh=_bounds(c,point+error)
            olow,ohigh=_bounds(c,box)
            assert low<=ohigh and nhigh>=olow
            old_peak_pairs+=1
    velocities=[]
    positive_swirl=0
    for sample in solver.numeric['levels'][-1]['samples']:
        for log_tau in ('-1','-4'):
            packet=solver.leading_velocity_terms(sample['Z'],sample['rho'],log_tau=log_tau,azimuth='0.5')
            assert packet['coefficient_orders_included']==[0]
            assert not packet['absolute_exponentials_materialized']
            assert not packet['numerical_float_velocity_vector_available']
            swirl=packet['cylindrical']['utheta']
            if sample['rho']=='0':
                assert c.mpf(swirl['coefficient']['approximation'])==0
                assert not swirl['log_scale_terms']
            else:
                assert c.mpf(swirl['coefficient']['approximation'])>0
                assert len(swirl['log_scale_terms'])==3
                assert swirl['log_scale_terms'][0]['negative_of_exact_mpf_tuple']==solver.selected_logC_tuple
                positive_swirl+=1
            velocities.append(packet)
    offgrid=solver.leading_velocity_terms('0.17','1.5',log_tau='-2',azimuth='0.25')
    assert offgrid['numerical_core_origin']=='on_demand_same_source_production_recurrence'
    assert set(offgrid['coordinates'])=={'r','x','y','z','t'}
    assert c.mpf(offgrid['cylindrical']['utheta']['coefficient']['approximation'])>0
    velocities.append(offgrid)
    positive_swirl+=1
    extra_hashes={peak_path.name:hashlib.sha256(peak_path.read_bytes()).hexdigest()}
    checks=dict(passed=True,precision_comparisons=comparisons,independent_nominal_real_integrals=len(integrals),
        existing_original_peak_source_pairs=old_peak_pairs,leading_physical_velocity_term_packets=len(velocities),
        positive_nonaxis_swirl_coefficient_packets=positive_swirl,
        offgrid_core_point_computed_from_original_recurrence=True,
        original_exact_logC_singleton=True,source_midpoint_selection=False,
        narrow_root_offsets_kept_separate=True,all_original_local_numerator_terms_retained=True,
        absolute_F0_materialized=False,full_float_velocity_field_available=False,temporal_recursion=False)
    paths=[Path(__file__),ROOT/'src/openai_ns_reconstruction/original_core_log_amplitude.py']
    report=dict(source_family=solver.family,levels=levels,independent_integral_diagnostics=integrals,
        leading_physical_velocity_terms=velocities,checks=checks,input_hashes={**solver.hashes,**extra_hashes},
        implementation_hashes={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        limitations=dict(global_numerical_F0_relative_error_control=False,absolute_F0_materialized=False,
            full_float_velocity_vector_available=False,leading_n0_only=True,n1_numerical_layer_included=False,
            global_field_available=False,measured_blowup_dynamics=False,temporal_recursion=False,
            full_corrected_NS_solution=False),execution_seconds=round(time.monotonic()-started,3))
    if args.write_report:path.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    selected=next(p for p in levels[-1]['resolved_peak'] if p['xi']=='1')
    print(json.dumps(dict(checks=checks,relative_peak_at_xi1=selected['F0_relative_to_exact_anchor'],
        execution_seconds=report['execution_seconds']),indent=2))


if __name__=='__main__':main()
