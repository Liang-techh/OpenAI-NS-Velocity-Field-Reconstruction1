"""Actual five-moment C0/Z histories through the complete original O2 route."""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_middle_O2_inlet_C1_histories as preceding
import lei_ren_part1_paper_compliant_current_native_serial_C1_cell as serial

p=preceding.p;transfer=preceding.transfer;history=preceding.history;density=preceding.density
prior=preceding.prior;native=preceding.native;packets=preceding.packets
HERE,PREFIX,sha=preceding.HERE,preceding.PREFIX,preceding.sha;ep=preceding.ep;RATES=preceding.RATES
NAME=PREFIX+'current_native_O2_C1_histories.json'
RECEIPT=PREFIX+'current_native_O2_C1_histories_check.json'
GATE='current_original_inlet_to_O2_exit_C1_history_covers_and_parameter_bounds_executed'
ROUTE=(('slope_0_12','O2_slope',0,'.12'),('slope_12_13','O2_slope','.12','.13'),
    ('slope_13_14','O2_slope','.13','.14'),('slope_14_15','O2_slope','.14','.15'),
    ('slope_15_1','O2_slope','.15',1),('axial','O2_axial',0,1),
    ('buffer_0_9','O2_buffer',0,9),('buffer_9_11','O2_buffer',9,11))
LOCAL=('slope_12_13','slope_13_14','slope_14_15')


def source_joins(family):
    name=PREFIX+'current_O2_background_tensor_check.json';r=json.loads((HERE/name).read_bytes())
    flags=('current_actual_O2_reference_slope_axial_buffer_full_tensors_available',
        'current_actual_O2_full_meridional_decomposition_available','current_actual_four_O2_completed_tensor_joins_certified')
    if not r['all_passed'] or any(r.get(k)!=v for k,v in family.items()) or not all(r.get(k) for k in flags):
        raise ValueError('Same accepted original O2 source/history/P0 joins required')
    if r['current_actual_O2_source_and_tensor_theorem']['source_function_equality_precedes_common_tensor_bounds'] is not True:
        raise ValueError('Original O2 source equality proof required')
    result={}
    for interface in ('reference_slope','slope_axial','axial_buffer','buffer_transition'):
        if interface not in r['current_actual_four_O2_tensor_interface_rows']:
            raise ValueError('Original O2 source interface absent: '+interface)
        result[interface]=dict(receipt=name,sha256=sha(name),source_family=family,
            required_original_join_flags=list(flags),original_interface=interface,
            original_source_function_equality_precedes_tensor_bounds=True)
    return result


def make_middle_owner(bridge):
    c1=density.NativeDensityC1LocalIntegrals(p.first.NativePhaseFirstJets(p.slow.NativeQSlowJets(
        p.current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge))))))
    second=preceding.preceding.preceding.NativeSecondBridgeC1Histories(p.NativeFirstBridgeC1Histories(
        transfer.NativeTrueChartC1Transfer(history.NativeC1HistoryTransfer(c1))))
    return preceding.NativeMiddleO2InletC1Histories(preceding.preceding.NativeBridgeSwitchC1Histories(second))


class NativeO2C1Histories:
    def __init__(self,owner):
        if type(owner) is not preceding.NativeMiddleO2InletC1Histories:raise ValueError('Accepted actual O2-inlet owner required')
        r=json.loads((HERE/preceding.RECEIPT).read_bytes())
        if not r['all_passed'] or not r[preceding.GATE] or r['source_family']!=owner.family:
            raise ValueError('Same original checked inlet-to-O2-inlet C1 source required')
        self.owner=owner;self.ctx=owner.ctx;self.family=owner.family;self.coordinates=owner.coordinates
        self.transfer=owner.transfer;self.service=owner.service;self.q_owner=owner.q_owner
        self.service.bind_hashes(r['input_hashes'])
        self.service.bind_hashes({name:sha(name) for name in (preceding.RECEIPT,Path(__file__).name,
            Path(serial.__file__).name,PREFIX+'current_O2_background_tensor_check.json',
            PREFIX+'current_O2_relaxed_buffer_cone_check.json')})
        self.joins=source_joins(self.family)

    @native.inlet.source_precision
    def route(self,Z=(-1,1),N=1024):
        N=density.spatial.density.candidate_integer(N)
        if N<160:raise ValueError('Original whole-period cover requires candidate N>=160')
        initial=self.owner.route(Z,N=N);incoming=initial['correction']
        previous_chart,previous_endpoint='O2_slope',0
        cells={};cumulative=history.C1DuhamelOperator(self.coordinates)
        local=history.C1DuhamelOperator(self.coordinates);local_incoming=None
        for label,chart,left,right in ROUTE:
            geometry=self.transfer.geometry.cell(chart,left,right)
            if chart==previous_chart:
                if density.spatial.exact_coordinate(left)!=density.spatial.exact_coordinate(previous_endpoint):
                    raise ValueError('Original adjacent same-chart endpoints required')
                seam='same '+chart+' chart at '+str(left)
                join=dict(same_original_chart_and_source_function=True,coordinate=str(left))
            else:
                seam=previous_chart+' -> '+chart
                if seam not in self.transfer.geometry.binder.identity['original_same_radius_periodic_phase_seam_identities']:
                    raise ValueError('Exact original O2 radius/phase seam required')
                join=self.joins['slope_axial' if chart=='O2_axial' else 'axial_buffer']
            if label==LOCAL[0]:local_incoming=incoming
            cell=serial.serial_cell(self,chart=chart,geometry=geometry,endpoint=self.ctx.mpf(right),Z=Z,N=N,
                incoming=incoming,left_record=dict(chart=previous_chart,coordinate=previous_endpoint),
                right_record=dict(chart=chart,coordinate=right),seam=seam,join=join)
            transfer.append_true_cell(cumulative,geometry,cell['contributions'],cell['Z_derivatives'],self.family)
            if label in LOCAL:transfer.append_true_cell(local,geometry,cell['contributions'],cell['Z_derivatives'],self.family)
            cell['record']['O2_inlet_to_current_endpoint_C1_operator']=cumulative.record()
            if chart=='O2_buffer':
                cell['record']['buffer_selector_to_shared_offset']='shared_offset=selector-11'
                cell['record']['smaller_relaxed_buffer_admission_scope']='selector[0,9] = shared_offset[-11,-2]'
                cell['record']['last_two_units_source_covered_without_extended_cone_admission']=left==9
                cell['record']['global_strict_or_relaxed_cone_admission_from_this_integral']=False
            cells[label]=cell;incoming=cell['correction'];previous_chart,previous_endpoint=chart,right
        seam='O2_buffer -> O3_slope_mu'
        if seam not in self.transfer.geometry.binder.identity['original_same_radius_periodic_phase_seam_identities']:
            raise ValueError('Original terminal buffer/transition radius and phase seam required')
        background=history.packet_history_functions(self.transfer.owner.native.query('O3_slope_mu',Z,0),
            self.coordinates,self.transfer.owner.signed_owner)
        own={k:background['originals'][k]+incoming['values'][k] for k in RATES}
        own_Z={k:background['Z_derivatives'][k]+incoming['Z_derivatives'][k] for k in RATES}
        record=dict(source_family=self.family,Z_box=self.ctx.mpf(Z),candidate_N=N,
            known_actual_original_O2_inlet_C1_history=initial['record'],
            ordered_original_O2_cells=[dict(label=label,chart=chart,left=left,right=right) for label,chart,left,right in ROUTE],
            actual_O2_serial_cell_C1_records={label:cell['record'] for label,cell in cells.items()},
            actual_local_O2_12_incoming_C0={k:v.record() for k,v in local_incoming['values'].items()},
            actual_local_O2_12_incoming_Z={k:v.record() for k,v in local_incoming['Z_derivatives'].items()},
            actual_local_O2_12_15_three_cell_C1_operator=local.record(),
            O2_inlet_to_O2_exit_C1_operator=cumulative.record(),
            original_terminal_O2_O3_radius_phase_seam=seam,
            original_terminal_O2_O3_background_history_P0_join=self.joins['buffer_transition'],
            actual_O3_transition0_background_and_separate_P0_Z=background['record'],
            actual_O2_exit_correction_C0={k:v.record() for k,v in incoming['values'].items()},
            actual_O2_exit_correction_Z={k:v.record() for k,v in incoming['Z_derivatives'].items()},
            actual_O3_transition0_own_history_C0={k:v.record() for k,v in own.items()},
            actual_O3_transition0_own_history_Z={k:v.record() for k,v in own_Z.items()},
            no_original_interval_skipped_from_true_inlet_to_O2_exit=True,
            actual_O2_12_incoming_then_original_three_adjacent_cells_applied=True,
            original_buffer_last_two_units_not_given_missing_cone_admission=True,
            new_periodic_derivative_bounds_preserve_original_loop_and_source_parameters=True,
            original_O3_to_Rc_route_still_missing=True,
            conservative_covers_not_quantitative_terminal_closure=True,
            global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,initial=initial,cells=cells,cumulative=cumulative,local=local,
            local_incoming=local_incoming,correction=incoming,own=own,own_Z=own_Z,background=background)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        owner=NativeO2C1Histories(make_middle_owner(bridge));records={}
        for name,Z in (('whole_Z',(-1,1)),('Z_interval',('.49','.51'))):
            records[name]=owner.route(Z,N=1024)['record'];print('Actual complete O2 C1 exit:',name,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},candidate_N=1024,native_Z_query_count=2,
        complete_added_original_chart_count=3,actual_original_O2_cell_count=len(ROUTE),
        actual_original_inlet_to_O2_exit_C1_records=records,
        original_periodic_parameter_C1_theorem=serial.periodic_parameter_theorem(),
        actual_O2_12_incoming_and_adjacent_three_cell_C1_transfer_installed=True,
        actual_original_inlet_to_O3_transition0_C1_covers_installed=True,
        original_O3_to_Rc_route_still_missing=True,global_inlet_to_Rc_histories_admitted=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Actual original inlet -> complete inner/middle route -> full O2 slope/axial/buffer -> directly queried O3 transition0; five C0/Z correction/own histories on whole Z[-1,1] and[.49,.51], with actual O2 .12 incoming and .12-.15 three-cell transfer. Parameter-uniform original primitive derivative bounds, exact widths and separate backgrounds/P0_Z. Conservative bounds; O3/Rc/terminal/common N/cone/recursion/full NS remain open.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
