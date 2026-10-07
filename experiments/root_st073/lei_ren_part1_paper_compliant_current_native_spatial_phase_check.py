"""Original spatial phase replay, modular arithmetic, radius seams and scope."""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_native_spatial_phase as current
import lei_ren_part1_paper_compliant_current_native_conditioned_phase_check as checks

phase=current.phase;density=current.density;native=current.native;packets=current.packets
HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;ep=current.ep
require=checks.require;contains=checks.contains


def modular_checks(c):
    p=mp.mp.clone();p.dps=c.dps+20;exponent=10**40
    selected=c.mpf(p.make_mpf((0,3,exponent,2)));mult=Fraction(-137,10000)
    value,proof=current.binary_mod_one(c,selected,mult)
    # Independent CRT: 10000=16*625; lambda(625)=500. The binary
    # exponent is huge, but neither reference nor producer creates 2**exponent.
    r625=(-411*pow(2,exponent%500,625))%625
    expected=Fraction((r625+625*((-r625)%16))%10000,10000)
    require(proof['exact_fraction']==current.fractional_record(expected),'Huge signed binary modulus differs from CRT')
    require(contains(value,c.mpf(expected.numerator)/expected.denominator),'Directed huge binary fraction missing')
    negative=c.mpf(p.make_mpf((0,13,-7,4)))
    value,proof=current.binary_mod_one(c,negative,Fraction(7,11))
    require(proof['exact_fraction']==current.fractional_record(Fraction(91,1408)),'Negative binary exponent modulus differs')
    for selected,mult in ((c.mpf((1,2)),Fraction(1,3)),(c.mpf(1),c.mpf('.3'))):
        try:current.binary_mod_one(c,selected,mult)
        except ValueError:continue
        raise ArithmeticError('Non-selected constant or non-exact multiplier accepted')
    require(current.exact_coordinate(c.mpf('.1337')) is None,'Interval coordinate was selected as an exact field value')
    wrapped=current.ordinary_mod_one(c,c.mpf(('.9','1.1')))
    require(not wrapped['full_period'] and len(wrapped['boxes'])==2,'Endpoint-crossing union lost')
    require(contains(wrapped['boxes'][0],c.mpf('.95')) and contains(wrapped['boxes'][1],c.mpf('.05')),'Wrapped phase cells incorrect')
    require(current.ordinary_mod_one(c,c.mpf((-2,3)))['full_period'],'Wide analytic phase was narrowed')
    return dict(passed=True,huge_binary_exponent_test=str(exponent),independent_signed_CRT_reference=True,
        negative_exponent_fraction_reference=True,nonselected_constants_and_interval_midpoints_rejected=True,
        wrapped_union_and_full_period_cover_checked=True,huge_integer_not_materialized=True)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_chart_count']==17,'All17 original radius charts required')
    require(not any(saved.get(k) for k in packets.OPEN),'Spatial candidate cannot complete global stages')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed spatial phase prerequisite: '+name)
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=current.NativeSpatialPhase(native.NativeGenericSourcePackets(bridge));N=saved['candidate_N'];c=owner.ctx
        independent=modular_checks(c);identity=owner.identity
        require(len(identity['exact_source_radius_minus_same_left_radius_identities'])==17,'Missing affine radius identities')
        require(len(identity['original_native_coordinate_Jacobian_identities'])==17,'Missing original coordinate Jacobians')
        require(len(identity['original_same_radius_periodic_phase_seam_identities'])==16,'Missing radius/phase seams')
        require(owner.theorem['original_axial_and_buffer_radius_coordinates_bound'],'Original actual_y source not bound')
        charts={}
        for chart,old in saved['current_original_radius_and_spatial_phase_records'].items():
            got=owner.query(chart,('.5','.5'),old['query_coordinate'],N)
            require(packets.encode(got['record'])==old['result'],'Original spatial phase source changed: '+chart)
            require(not got['record']['periodic_projection']['full_period'],'Explicit point phase unexpectedly full: '+chart)
            require(not any(got['record'].get(k) for k in packets.OPEN),'Local phase cannot admit a global stage')
            require(got['record']['source_radius_caps_not_consumed'] and got['record']['native_width_conversion_not_reapplied'],'Original affine source/radial conversion required')
            if chart in ('bridge_first','bridge_second'):
                require(not got['offset'].zero and ep(got['offset'].coefficient)[0]>0,'Positive microscopic radius offset dropped')
            charts[chart]=dict(phase_cell_count=len(got['phase_boxes']),bounded_actual_source_phase=True)
        inlet=owner.query('bridge_first',(-1,1),{'selected_sc_multiple':'1/2'},N)
        require(inlet['offset'].zero and all(ep(v)==(0,0) for v in inlet['phase_boxes']),'Same source left radius must give exact zero phase')
        require(packets.encode(inlet['record'])==saved['actual_same_left_inlet'],'Original inlet source changed')
        cover=owner.query('inner_reference',('.5','.5'),c.mpf('.1337'),N)
        require(cover['record']['periodic_projection']['full_period'],'Analytic coordinate enclosure was silently replaced by exact literal')
        candidate=density.NativeCandidateDensities(phase.NativeConditionedPhase(current.current.NativeCorrelatedShearQ(current.prior.NativeSignedInputEnclosures(owner.native))))
        cells=0;active=0
        for chart,old in saved['actual_spatial_candidate_velocity_and_density_records'].items():
            coordinate=saved['current_original_radius_and_spatial_phase_records'][chart]['query_coordinate']
            got=current.spatial_candidate_functions(candidate,owner,chart,('.5','.5'),coordinate,N)
            require(packets.encode(got)==old,'Spatial candidate functions changed: '+chart)
            for cell in got['actual_spatial_candidate_cells']:
                require(cell['status']=='enclosed' and cell['spatial_phase_binding_installed'],'Actual derived phase unavailable')
                require(cell['phase_box_is_derived_not_independently_selected'],'Spatial phase substituted by free parameter')
                require(not cell['original_phase_inverse_and_primitives']['free_phase_parameter_not_spatial_phase'],'Free phase flag not rebound')
                require(set(cell['five_signed_dimensionless_Duhamel_density_kernels'])==set(density.RATES),'Missing original signed kernel')
                require(cell['candidate_N_is_not_global_common_N_admission'] and not cell['density_integrals_or_new_moment_histories_constructed'],'Candidate phase overstated as a global/integral admission')
                require(not any(cell.get(k) for k in packets.OPEN),'Local spatial candidate cannot admit global stages')
                active+=any(not v['exact_zero'] for v in cell['five_signed_dimensionless_Duhamel_density_kernels'].values())
                cells+=1
            print('Actual spatial candidate functions checked:',chart,flush=True)
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},
        native_radius_phase_queries_checked=len(charts),native_Jacobian_identities_checked=17,
        original_radius_phase_seam_identities_checked=16,exact_same_left_inlet_phase_checked=True,
        exact_literal_and_analytic_coordinate_cover_distinguished=True,microscopic_radius_source_retained=True,
        spatial_candidate_source_box_count=5,spatial_candidate_cell_queries_checked=cells,
        signed_spatial_density_kernel_enclosures_checked=cells*5,nonzero_spatial_candidate_cells=active,
        candidate_N=N,independent_modular_arithmetic_checks=independent,charts=charts,
        global_frequency_or_whole_velocity_coverage_or_integral_or_repair_admission=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original native spatial phase PASS:',len(charts),'charts;',cells,'spatial candidate cells',flush=True)
    return result


if __name__=='__main__':run()
