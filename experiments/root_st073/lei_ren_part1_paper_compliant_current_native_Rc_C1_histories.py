"""Actual original inlet-to-Rc C0/Z history covers, with original quiet power.

Rc is original power offset2 (phase2/Tw), not the power chart endpoint.
Conservative function-domain covers do not establish terminal closure.
"""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_O2_C1_histories as preceding
import lei_ren_part1_paper_compliant_current_native_serial_C1_cell as serial

p=preceding.p;transfer=preceding.transfer;history=preceding.history;density=preceding.density
prior=preceding.prior;native=preceding.native;packets=preceding.packets
HERE,PREFIX,sha=preceding.HERE,preceding.PREFIX,preceding.sha;ep=preceding.ep;RATES=preceding.RATES
ZERO=serial.ZERO;DZ=serial.DZ;whole_period_C1=serial.whole_period_C1
NAME=PREFIX+'current_native_Rc_C1_histories.json';RECEIPT=PREFIX+'current_native_Rc_C1_histories_check.json'
GATE='current_original_inlet_to_actual_Rc_C1_history_covers_executed'
ROUTE=(('transition','O3_slope_mu',0,1),('power_to_r_plus','O3_power',0,1),('power_to_Rc','O3_power',1,2))

def original_O3_cell(owner,*,chart,geometry,endpoint,Z,N,incoming,left_record,right_record,seam,join):
    """Same original density/mass/recovery, explicit inherited correction."""
    if chart not in ('O3_slope_mu','O3_power'):raise ValueError('Original O3 source chart required')
    owner.coordinates.require_family(owner.family)
    if N<160:raise ValueError('Whole-period C0 cap requires candidate N>=160')
    if not geometry['record']['width_and_endpoints_independent_of_Z']:
        raise ValueError('Original fixed-Z integration geometry required')
    c=owner.ctx;source=owner.q_owner.query(chart,Z,geometry['coordinate'])
    roots=source['source']['roots'];root_owner=owner.q_owner.owner.owner
    positive=root_owner.decode(root_owner.inventory[chart]['actual_positive_denominator_theorem'])
    if positive.get('whole_actual_source_positive_not_inferred_from_saved_denominator_box') is not True:
        raise ValueError('Original O3 chart-uniform source positive a theorem required')
    eta=packets.interval(c,root_owner.scales['selected_positive_eta_log'])
    dstar=packets.interval(c,root_owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
    primitives=whole_period_C1(roots,eta,positive['log_actual_a_positive_lower'],dstar)
    E,E_Z=roots['E'][ZERO],roots['E'][DZ];packet=source['source']['packet']
    signed=owner.transfer.owner.signed_owner
    def axial(k):
        row=prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
        return signed.leaf(prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),E.scale.bases,E.ledger)
    got=density.density_Z_kernels(E,E_Z,axial(0),axial(1),primitives['values'],N)
    factors={k:transfer.true_width_kernel(owner.coordinates,geometry,rate) for k,rate in RATES.items()}
    values={k:owner.coordinates.rebase(got['kernels'][k],owner.family)*factors[k]['mass'] for k in RATES}
    jets={k:owner.coordinates.rebase(got['Z_derivatives'][k],owner.family)*factors[k]['mass'] for k in RATES}
    operator=history.C1DuhamelOperator(owner.coordinates)
    transfer.append_true_cell(operator,geometry,values,jets,owner.family)
    correction=operator.apply(incoming['values'],incoming['Z_derivatives'],owner.family)
    background=history.packet_history_functions(owner.transfer.owner.native.query(chart,Z,endpoint),owner.coordinates,signed)
    own={k:background['originals'][k]+correction['values'][k] for k in RATES}
    own_Z={k:background['Z_derivatives'][k]+correction['Z_derivatives'][k] for k in RATES}
    record=dict(source_family=owner.family,chart=chart,Z_box=c.mpf(Z),candidate_N=N,
        original_left_endpoint=left_record,original_right_endpoint=right_record,
        original_same_radius_phase_or_same_chart_seam=seam,original_background_history_P0_join_receipt=join,
        actual_true_cell_geometry=geometry['record'],original_chart_uniform_a_positive_certificate=positive,
        original_full_box_signed_source_and_old_q_slow_status=source['record'],
        original_periodic_parameter_C1_cover=primitives['record'],
        original_true_width_kernel_factors={k:dict(branch=v['branch'],decay=v['decay'].record(),mass=v['mass'].record()) for k,v in factors.items()},
        actual_signed_density_C0_covers={k:v.record() for k,v in got['kernels'].items()},
        actual_signed_density_Z_covers={k:v.record() for k,v in got['Z_derivatives'].items()},
        actual_true_width_C0_contributions={k:v.record() for k,v in values.items()},
        actual_true_width_Z_contributions={k:v.record() for k,v in jets.items()},
        actual_inherited_correction_C0={k:v.record() for k,v in incoming['values'].items()},
        actual_inherited_correction_Z={k:v.record() for k,v in incoming['Z_derivatives'].items()},
        actual_right_correction_C0={k:v.record() for k,v in correction['values'].items()},
        actual_right_correction_Z={k:v.record() for k,v in correction['Z_derivatives'].items()},
        original_right_background_and_separate_P0_Z=background['record'],
        actual_right_own_history_C0={k:v.record() for k,v in own.items()},
        actual_right_own_history_Z={k:v.record() for k,v in own_Z.items()},
        actual_cell_C1_operator=operator.record(),
        old_unresolved_q_rows_not_used_for_first_Z_primitive_bounds=True,
        original_true_width_and_native_ordinary_y_conversion_applied_once=True,
        actual_incoming_correction_not_reset_or_background_double_added=True,
        global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False))
    return dict(record=record,geometry=geometry,source=source,primitives=primitives,kernels=got,
        factors=factors,operator=operator,incoming=incoming,contributions=values,Z_derivatives=jets,
        correction=correction,background=background,own=own,own_Z=own_Z)


def original_O3_join(family):
    name=PREFIX+'current_O3_transition_background_tensor_check.json';r=json.loads((HERE/name).read_bytes())
    flags=('current_actual_O3_transition_full_tensor_available','current_actual_O3_transition_power_completed_tensor_join_certified')
    if not r['all_passed'] or any(r.get(k)!=v for k,v in family.items()) or not all(r.get(k) for k in flags):
        raise ValueError('Same accepted original O3 source join required')
    proof=r['current_actual_O3_transition_source_and_tensor_theorem']
    endpoint=proof['current_actual_transition_source_and_power_endpoint']
    identities=endpoint['exact_transition_power_source_endpoint']['exact_five_history_interface_identities']['Rw_slope_mu_power']
    required=('actual_transition_power_same_Rw','actual_same_parent_all_five_histories_absolute_P0_and_flat_logU_jets')
    if not proof['source_function_equality_precedes_common_tensor_bounds'] or not all(endpoint[k] for k in required) or not all(identities[k] for k in RATES):
        raise ValueError('Original source-function/five-history/P0 seam proof required')
    return dict(receipt=name,sha256=sha(name),source_family=family,original_interface='Rw_slope_mu_power',
        original_source_function_equality_precedes_tensor_bounds=True,
        actual_same_parent_all_five_histories_absolute_P0_and_flat_logU_jets=True,
        exact_five_history_interface_identities=identities)


class NativeRcC1Histories:
    def __init__(self,owner):
        if type(owner) is not preceding.NativeO2C1Histories:raise ValueError('Accepted actual O3-inlet owner required')
        r=json.loads((HERE/preceding.RECEIPT).read_bytes())
        if not r['all_passed'] or not r[preceding.GATE] or r['source_family']!=owner.family:
            raise ValueError('Same checked original inlet-to-O2-exit C1 source required')
        self.owner=owner;self.ctx=owner.ctx;self.family=owner.family;self.coordinates=owner.coordinates
        self.transfer=owner.transfer;self.service=owner.service;self.q_owner=owner.q_owner
        self.service.bind_hashes(r['input_hashes'])
        self.service.bind_hashes({name:sha(name) for name in (preceding.RECEIPT,Path(__file__).name,
            PREFIX+'current_O3_transition_background_tensor_check.json',PREFIX+'current_O3_power_cone_check.json',
            PREFIX+'current_generic_shear_O3_sources_check.json',PREFIX+'current_generic_shear_loop_domain.json',
            PREFIX+'current_generic_shear_loop_domain_check.json')})
        self.join=original_O3_join(self.family)
        domain=json.loads((HERE/(PREFIX+'current_generic_shear_loop_domain.json')).read_bytes())
        receipt=json.loads((HERE/(PREFIX+'current_generic_shear_loop_domain_check.json')).read_bytes())
        if not receipt['all_passed'] or not receipt['current_original_generic_loop_two_sided_collars_and_repair_geometry_certified'] or domain['source_family']!=self.family:
            raise ValueError('Original same-family loop-edge/Rc geometry required')
        self.domain=domain['current_original_generic_loop_domain'];reserve=self.domain['original_r_plus_through_twice_Rc_admission_and_reservation']
        if not reserve['original_whole_power_cone_consumed_on_r_plus_through_twice_Rc'] or not reserve['strict_excess_not_computed_by_subtracting_rounded_a_minus2']:
            raise ValueError('Original canonical power and Rc reservation required')
        offsets=reserve['log_radius_offsets_from_same_Rw']
        if ep(packets.interval(self.ctx,offsets['r_plus']))!=(1,1) or ep(packets.interval(self.ctx,offsets['Rc']))!=(2,2):
            raise ValueError('Original r_plus/Rc offsets1/2 required')
        self.reservation=reserve

    @native.inlet.source_precision
    def route(self,Z=(-1,1),N=1024):
        N=density.spatial.density.candidate_integer(N)
        if N<160:raise ValueError('Original whole-period cover requires N>=160')
        initial=self.owner.route(Z,N);incoming=initial['correction'];cells={}
        cumulative=history.C1DuhamelOperator(self.coordinates)
        previous_chart,previous_endpoint='O3_slope_mu',0
        Tw=self.transfer.geometry.binder.fixed['Tw']
        for label,chart,left,right in ROUTE:
            if chart=='O3_power':
                lo=dict(original_power_offset=left);hi=dict(original_power_offset=right)
                geometry=self.transfer.geometry.cell(chart,lo,hi);endpoint=self.ctx.mpf(right)/Tw
                left_record=dict(chart=previous_chart,coordinate=previous_endpoint) if left==0 else dict(chart=chart,original_power_offset=left)
                right_record=dict(chart=chart,original_power_offset=right,original_phase_expression=str(right)+'/Tw')
            else:
                geometry=self.transfer.geometry.cell(chart,left,right);endpoint=self.ctx.mpf(right)
                left_record=dict(chart=chart,coordinate=left);right_record=dict(chart=chart,coordinate=right)
            if previous_chart==chart:
                seam='same original '+chart+' source'
                join=dict(same_original_chart_and_source_function=True,original_power_offset=left) if chart=='O3_power' else dict(same_original_chart_and_source_function=True,coordinate=left)
            else:
                seam=previous_chart+' -> '+chart
                if seam not in self.transfer.geometry.binder.identity['original_same_radius_periodic_phase_seam_identities']:
                    raise ValueError('Exact original O3 radius/phase seam required')
                join=self.join
            cell=original_O3_cell(self,chart=chart,geometry=geometry,endpoint=endpoint,Z=Z,N=N,incoming=incoming,
                left_record=left_record,right_record=right_record,seam=seam,join=join)
            transfer.append_true_cell(cumulative,geometry,cell['contributions'],cell['Z_derivatives'],self.family)
            cell['record']['O3_inlet_to_current_endpoint_C1_operator']=cumulative.record()
            if chart=='O3_power':
                if not all(v.zero for v in cell['primitives']['values'].values()) or not all(v.zero for v in [*cell['kernels']['kernels'].values(),*cell['kernels']['Z_derivatives'].values()]):
                    raise ValueError('Original complete power q flat proof and exact zero increments required')
                cell['record']['original_whole_power_q_q_Z_and_density_C0_Z_exact_zero']=True
                cell['record']['original_canonical_power_excess_source']='2*mu; all Z derivatives exact zero'
                cell['record']['original_actual_power_offsets']=[left,right]
                cell['record']['quiet_local_density_does_not_reset_inherited_histories']=True
            cells[label]=cell;incoming=cell['correction'];previous_chart,previous_endpoint=chart,right
        last=cells['power_to_Rc'];background=last['background']
        record=dict(source_family=self.family,Z_box=self.ctx.mpf(Z),candidate_N=N,
            known_actual_original_O3_inlet_C1_history=initial['record'],
            ordered_original_O3_Rc_cells=[dict(label=label,chart=chart,left=left,right=right,
                coordinate_kind='original_power_offset' if chart=='O3_power' else 'original_transition_offset') for label,chart,left,right in ROUTE],
            actual_O3_Rc_serial_C1_records={label:cell['record'] for label,cell in cells.items()},
            O3_inlet_to_actual_Rc_C1_operator=cumulative.record(),
            original_Rc_endpoint=dict(radius_expression='Rw*exp(2)',original_power_offset=2,
                original_power_phase_expression='2/Tw',actual_direct_endpoint_phase_cover=2/Tw,
                original_Rc_endpoint_not_whole_power_phase1=True),
            original_right_collar_inside_power=dict(original_power_offset_box=self.ctx.mpf(('.75','1.25')),
                original_r_plus_offset=1,whole_original_power_q_flat_source_proof_covers_both_sides=True),
            original_Rc_background_and_separate_P0_Z=background['record'],
            actual_Rc_correction_C0={k:v.record() for k,v in incoming['values'].items()},
            actual_Rc_correction_Z={k:v.record() for k,v in incoming['Z_derivatives'].items()},
            actual_Rc_own_history_C0={k:v.record() for k,v in last['own'].items()},
            actual_Rc_own_history_Z={k:v.record() for k,v in last['own_Z'].items()},
            no_original_interval_skipped_from_true_inlet_to_actual_Rc=True,
            original_whole_power_flat_support_preserves_inherited_pressure_memory=True,
            actual_original_Rc_conservative_function_domain_C0_Z_covers_available=True,
            quantitative_terminal_five_Z_identities_not_established=True,
            original_outer_2Rc_to_Rb_admission_still_missing=True,
            global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,initial=initial,cells=cells,cumulative=cumulative,correction=incoming,
            own=last['own'],own_Z=last['own_Z'],background=background)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=NativeRcC1Histories(preceding.NativeO2C1Histories(preceding.make_middle_owner(bridge)));records={}
        for name,Z in (('whole_Z',(-1,1)),('Z_interval',('.49','.51'))):
            records[name]=owner.route(Z,N=1024)['record'];print('Actual original Rc C1 covers:',name,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},candidate_N=1024,native_Z_query_count=2,
        actual_added_O3_Rc_cell_count=3,actual_original_inlet_to_Rc_C1_records=records,
        original_O3_transition_power_source_join=owner.join,original_Rc_reservation=owner.reservation,
        actual_original_Rc_conservative_function_domain_C0_Z_covers_available=True,
        original_entire_inlet_to_Rc_radial_route_covered=True,
        quantitative_terminal_five_Z_identities_not_established=True,
        global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Actual original inlet sc/2 through every inner/middle/O2 interval, full O3 transition, power offsets0->1->2 to directly queried Rc phase2/Tw; five correction/own C0/Z function-domain covers and exact power flat increments. Conservative enclosures; quantitative terminal identities, full loop/outer admission, global N/cone, higher jets/heat/energy/recursion/full NS remain open.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
