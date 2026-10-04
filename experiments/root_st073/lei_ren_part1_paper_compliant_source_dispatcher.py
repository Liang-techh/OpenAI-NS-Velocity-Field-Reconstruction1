"""Callable same-source leading profile charts, analytic core to heat exterior.

This is a chart dispatcher, NOT a pointwise physical Cartesian field.
It exposes each chart's actual derivative coordinates, normalization and
formal source scales instead of silently equating their numerical boxes.
Usage: field.evaluate('O2_axial', Z='.5', coordinate='.5'). Charts are
selected explicitly because the original huge radii remain formal.
"""
import hashlib
import importlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
# provider, class, method, coordinate description, legal domain, receipt,
# extra method argument (if any). Tuple order is the original radial chain.
ROUTES={
    'core':('core_physical_field','CompliantCorePhysicalField','profiles','rho=R/epsilon_core','[0,4]','core_physical_field_check',None),
    'bridge_first':('bridge_mixed_C4','CompliantBridgeMixedC4','evaluate','s=log(R/Ra)/hb','[0,1]','bridge_mixed_C4_check','first'),
    'bridge_second':('bridge_mixed_C4','CompliantBridgeMixedC4','evaluate','s=log(R/Ra)/hb','[1,2]','bridge_mixed_C4_check','second'),
    'bridge_macro':('bridge_mixed_C4','CompliantBridgeMixedC4','evaluate','fraction of [2hb,log(100/Ra)]','[0,1]','bridge_mixed_C4_check','macro'),
    'switch_first':('microswitch_mixed_C4','CompliantMicroswitchMixedC4','evaluate','s=log(R/100)/hb','[0,1]','microswitch_mixed_C4_check','first'),
    'switch_second':('microswitch_mixed_C4','CompliantMicroswitchMixedC4','evaluate','s=log(R/100)/hb','[1,2]','microswitch_mixed_C4_check','second'),
    'switch_power':('microswitch_mixed_C4','CompliantMicroswitchMixedC4','postpower','fraction from R2 to110','[0,1]','microswitch_mixed_C4_check',None),
    'reshape':('long_reshape_mixed_C4','CompliantLongReshapeMixedC4','evaluate','log(R/110)/T, original T=400 Abar','[0,1]','long_reshape_mixed_C4_check',None),
    'inner_reference':('reference_restore_mixed_C4','CompliantReferenceRestoreMixedC4','reference_branch','fraction from Rsh toRz','[0,1]','reference_restore_mixed_C4_check',None),
    'axial_restore':('reference_restore_mixed_C4','CompliantReferenceRestoreMixedC4','restoration','log(R/Rz)','[0,1]','reference_restore_mixed_C4_check',None),
    'restore_buffer':('reference_restore_mixed_C4','CompliantReferenceRestoreMixedC4','postrestore','log(R/Rref)','[-7,-6]','reference_restore_mixed_C4_check',None),
    'actual_patch':('actual_patch_mixed_C4','CompliantActualPatchMixedC4','evaluate','x=R/Rm','[1,e]','actual_patch_mixed_C4_check',None),
    'Rh_reference':('pre_pulse_mixed_C4','CompliantPrePulseMixedC4','reference','log(R/Rref)','[-5,0]','pre_pulse_mixed_C4_check',None),
    'O2_slope':('pre_pulse_mixed_C4','CompliantPrePulseMixedC4','slope','y=log(R/Rref)','[0,1]','pre_pulse_mixed_C4_check',None),
    'O2_axial':('pre_pulse_mixed_C4','CompliantPrePulseMixedC4','axial','phase=log(y)/Md','[0,1]','pre_pulse_mixed_C4_check',None),
    'O2_buffer':('pre_pulse_mixed_C4','CompliantPrePulseMixedC4','axial','offset=y-exp(Md)','[0,11]','pre_pulse_mixed_C4_check','buffer_offset'),
    'O3_slope_mu':('pre_pulse_mixed_C4','CompliantPrePulseMixedC4','slope_mu','log(R/Rd)','[0,1]','pre_pulse_mixed_C4_check',None),
    'O3_power':('pre_pulse_mixed_C4','CompliantPrePulseMixedC4','power','log(R/Rw)/Tw','[0,1]','pre_pulse_mixed_C4_check',None),
    'pulse_entrance':('pulse_mixed_C4','CompliantPulseMixedC4','entrance','y=log(R/Rp), original smooth entrance','[0,.02/mu]','pulse_mixed_C4_check',None),
    'pulse_main':('pulse_mixed_C4','CompliantPulseMixedC4','main','xi=mu*log(R/Rp)','[.02,10]','pulse_mixed_C4_check',None),
    'pulse_exit':('pulse_mixed_C4','CompliantPulseMixedC4','main','xi=mu*log(R/Rp)','[10,11]','pulse_mixed_C4_check',None),
    'pulse_gap':('pulse_mixed_C4','CompliantPulseMixedC4','gap','xi=mu*log(R/Rp)','[11,12]','pulse_mixed_C4_check',None),
    'pulse_gap_end':('pulse_mixed_C4','CompliantPulseMixedC4','gap_from_end','s=log(R/Rv)','[-1/mu,-4]','pulse_mixed_C4_check',None),
    'pulse_end':('pulse_mixed_C4','CompliantPulseMixedC4','end','s=log(R/Rv)','[-4,0]','pulse_mixed_C4_check',None),
    'flatten':('flatten_mixed_C4','CompliantFlattenMixedC4','flatten','log(R/Rv)','[0,100]','flatten_mixed_C4_check',None),
    'outer_power':('power_angular_C4','CompliantPowerAngularC4','power','postflatten power fraction','[0,1]','power_angular_C4_check',None),
    'outer_angular':('power_angular_C4','CompliantPowerAngularC4','angular','log(R/Rrel)','[-4,0]','power_angular_C4_check',None),
    'steep_entry':('steep_waiting_C4','CompliantSteepWaitingC4','steep_in','log(R/Rrel)','[0,1]','steep_waiting_C4_check',None),
    'steep_power':('steep_waiting_C4','CompliantSteepWaitingC4','steep_power','fraction of original Ts steep power','[0,1]','steep_waiting_C4_check',None),
    'steep_exit':('steep_waiting_C4','CompliantSteepWaitingC4','steep_out','original one-unit exit offset','[0,1]','steep_waiting_C4_check',None),
    'waiting':('steep_waiting_C4','CompliantSteepWaitingC4','waiting','fraction of selected waiting interval','[0,1]','steep_waiting_C4_check',None),
 'heat_collar':('collar_pressure_C4','CompliantCollarPressureC4','collar','t=log(R/Rtail)','[0,3]','collar_pressure_C4_check',None),
 'heat_exterior':('heat_stress_C4','CompliantHeatStressC4','exterior','t=log(R/Rtail)','[3,infinity)','heat_stress_C4_check',None),
}


class CompliantSourceDispatcher:
    def __init__(self):
        self.providers={};self.receipts={};self.hashes={}
        self.family=None;self.source=None

    def provider(self,chart):
        if chart not in ROUTES:raise ValueError('Unknown source chart: '+str(chart))
        stem,cls,method,coordinate,domain,receipt,extra=ROUTES[chart]
        if receipt not in self.receipts:
            name=PREFIX+receipt+'.json';record=json.loads((HERE/name).read_bytes())
            if not record.get('all_passed'):raise ValueError('Source dispatcher requires accepted chart: '+receipt)
            family=record['actual_five_defect_family_sha256'];source=record['implicit_source_sha256']
            if self.family is None:self.family,self.source=family,source
            if (family,source)!=(self.family,self.source):raise ValueError('Mixed source families in dispatcher')
            for path,digest in record['input_hashes'].items():
                if path not in self.hashes:
                    if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Dispatcher source changed: '+path)
                    self.hashes[path]=digest
                elif self.hashes[path]!=digest:raise ValueError('Inconsistent chart source hash: '+path)
            self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest();self.receipts[receipt]=record
        if stem not in self.providers:self.providers[stem]=getattr(importlib.import_module(PREFIX+stem),cls)()
        return self.providers[stem]

    def evaluate(self,chart,Z,coordinate):
        provider=self.provider(chart);stem,cls,method,description,domain,receipt,extra=ROUTES[chart]
        fn=getattr(provider,method)
        with mp.workdps(280):
            if chart in ('core','actual_patch'):packet=fn(coordinate,Z)
            elif extra=='buffer_offset':packet=fn(Z,buffer_offset=coordinate)
            elif extra is not None:packet=fn(Z,coordinate,extra)
            else:packet=fn(Z,coordinate)
        grids={key:value for key,value in packet.items() if key in (
            'ordinary_mixed_profile_grids','physical_velocity_pressure_phase_Z_mixed4',
            'physical_velocity_pressure_y_Z_mixed4','physical_velocity_pressure_x_Z_mixed4',
            'physical_five_primitive_phase_Z_mixed4','physical_five_primitive_y_Z_mixed4',
            'physical_five_primitive_x_Z_mixed4','physical_mixed_derivatives_total_order_le4')}
        if not grids:raise ValueError('Chart did not return its actual mixed4 source packet: '+chart)
        return dict(chart=chart,actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            source_coverage_coordinate=description,source_coordinate_domain=domain,
            derivative_coordinate=('rho,Z' if chart=='core' else 'original phase s,Z with factored logR conversion' if chart in (
                'bridge_first','bridge_second','switch_first','switch_second') else 'y=logR,Z; selector fraction is not derivative coordinate'),
            physical_mixed_grids=grids,
            original_scale_metadata=packet.get('exact_formal_prefactors',dict(
                rule='Use component normalization named in the source grid plus original positive amplitude/radius sources in source_packet; no common numerical rescaling has been applied')),
            source_packet=packet,acceptance_receipt=PREFIX+receipt+'.json',
            heat_absolute_pressure_companion_used=chart in ('heat_collar','heat_exterior') and packet.get('absolute_pressure_same_source_mixed4_available',False),
            heat_exterior_similarity_stress_identity_available=chart=='heat_exterior' and packet.get('heat_exterior_stress_identity_certified',False),
            output_kind='directed source-field enclosures with formal scales; not point coefficients or a physical Cartesian evaluation',
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,temporal_recursion=False)

    def manifest(self):
        for chart in ROUTES:self.provider(chart)
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            ordered_chart_registry={name:dict(provider=PREFIX+spec[0],method=spec[2],coverage_coordinate=spec[3],domain=spec[4],acceptance_receipt=PREFIX+spec[5]+'.json') for name,spec in ROUTES.items()},
            all_profile_source_charts_callable=True,core_to_Rp_pre_pulse_source_chain_available=True,
            same_actual_inner_five_moment_source_and_pressure_family=True,
            full_collar_and_exterior_absolute_pressure_companions_adopted=True,
            source_bound_exterior_similarity_stress_companion_adopted=True,
            actual_physical_heat_stress_transfer_certified=False,
            all_absolute_radii_automatically_dispatched=False,uniform_physical_units_assembled=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,
            original_unlocalized_whole_space_energy_infinite=True,temporal_recursion=False,input_hashes=dict(self.hashes))


def run():
    field=CompliantSourceDispatcher();result=field.manifest();samples={}
    values={'core':'2','bridge_first':'.5','bridge_second':'1.5','bridge_macro':'.5',
        'switch_first':'.5','switch_second':'1.5','switch_power':'.5','reshape':'.5','inner_reference':'.5',
        'axial_restore':'.5','restore_buffer':'-6.5','actual_patch':'1.5','Rh_reference':'-2',
        'O2_slope':'.5','O2_axial':'.5','O2_buffer':'5','O3_slope_mu':'.5','O3_power':'.5',
        'pulse_entrance':'0','pulse_main':'5','pulse_exit':'10.5','pulse_gap':'11.5','pulse_gap_end':'-5',
        'pulse_end':'-2','flatten':'50','outer_power':'.5','outer_angular':'-2','steep_entry':'.5',
        'steep_power':'.5','steep_exit':'.5','waiting':'.5','heat_collar':'1.5','heat_exterior':'4'}
    count=0
    with mp.workdps(280):
        for chart in ROUTES:
            packet=field.evaluate(chart,'.5',values[chart]);rows=0
            for group in packet['physical_mixed_grids'].values():
                for grid in group.values():
                    for value in grid.values():
                        if any(not mp.isfinite(v) for v in endpoints(value)):raise ArithmeticError('Nonfinite dispatched source: '+chart)
                        rows+=1
            samples[chart]=dict(coverage_sample=values[chart],source_grids=packet['physical_mixed_grids'],finite_derivative_rows=rows,
                derivative_coordinate=packet['derivative_coordinate'],output_kind=packet['output_kind'],
                heat_absolute_pressure_companion_used=packet['heat_absolute_pressure_companion_used'],
                heat_exterior_similarity_stress_identity_available=packet['heat_exterior_similarity_stress_identity_available'])
            count+=rows;print('Dispatched actual source chart: '+chart,flush=True)
    result.update(actual_chart_evaluations=samples,all_routes_exercised=True,total_finite_derivative_rows=count,
        all_passed=True,input_hashes=dict(field.hashes))
    result['input_hashes'][Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Same-source core-to-heat dispatcher PASS: '+str(len(samples))+' charts',flush=True)
    return result


if __name__=='__main__':run()
