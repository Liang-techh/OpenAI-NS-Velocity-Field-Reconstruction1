"""Actual serial C1 histories through middle charts to the original O2 inlet.

The source chart order, positive denominator theorems and exact seams are
original. Whole-period primitive bounds and true log-radius masses act on
signed density functions, retaining actual R110 correction memory.
"""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_bridge_switch_C1_histories as preceding

p=preceding.p;transfer=preceding.transfer;history=preceding.history
density=preceding.density;prior=preceding.prior;native=preceding.native;packets=preceding.packets
HERE,PREFIX,sha=preceding.HERE,preceding.PREFIX,preceding.sha
RATES=preceding.RATES;ep=preceding.ep;ZERO=preceding.ZERO;DZ=preceding.DZ
NAME=PREFIX+'current_native_middle_O2_inlet_C1_histories.json'
RECEIPT=PREFIX+'current_native_middle_O2_inlet_C1_histories_check.json'
GATE='current_original_inlet_to_O2_inlet_middle_C1_history_covers_executed'
ROUTE=(('reshape',0,1),('inner_reference',0,1),('axial_restore',0,1),
    ('restore_buffer',-7,-6),('actual_patch',1,'e'),('Rh_reference',-5,0))


JOIN_SPECS={
    'reshape':('current_switch_power_background_tensor',
        ('current_actual_switch_power_full_tensor_available','current_actual_switch_power_full_meridional_decomposition_available','current_actual_R110_completed_tensor_join_certified'),
        'current_actual_R110_tensor_interface_rows','power_reshape',
        ('current_actual_switch_power_source_and_tensor_theorem','current_actual_switch_power_source_pressure_and_R110_theorem','same_full_current_R2_histories_V110_F0_and_P0_retained')),
    'inner_reference':('current_reshape_background_tensor',
        ('current_actual_long_reshape_full_tensor_available','current_actual_long_reshape_full_meridional_decomposition_available','current_actual_Rsh_completed_tensor_join_certified'),
        'current_actual_Rsh_tensor_interface_rows','reshape_reference',
        ('current_actual_reshape_source_units_and_tensor_theorem','current_actual_variable_reshape_source_pressure_and_Rsh_theorem','source_function_identity_precedes_completed_tensor_triangle_bounds')),
    'axial_restore':('current_restore_background_tensor',
        ('current_actual_inner_reference_axial_restore_buffer_full_tensors_available','current_actual_restore_full_meridional_decomposition_available','current_actual_three_restore_completed_tensor_joins_certified'),
        'current_actual_three_restore_tensor_interface_rows','reference_restore',
        ('current_actual_restore_source_units_and_tensor_theorem','current_actual_restore_source_pressure_and_three_endpoint_theorem','actual_current_Rm_open_neighborhood_function_theorem','same_P0_and_source_functions_arbitrary_smooth_Z')),
    'Rh_reference':('current_actual_patch_background_tensor',
        ('current_actual_five_moment_patch_full_tensor_available','current_actual_five_moment_patch_full_meridional_decomposition_available','current_actual_patch_Rh_completed_tensor_join_certified','current_actual_patch_six_internal_support_tensor_traces_certified'),
        'current_actual_patch_Rh_and_six_support_tensor_interface_rows','patch_Rh',
        ('current_actual_patch_source_units_and_tensor_theorem','current_actual_patch_source_pressure_support_and_Rh_theorem','same_analytic_P0_function_and_first_six_projection_not_hash_only'))}
JOIN_SPECS['restore_buffer']=(*JOIN_SPECS['axial_restore'][:3],'restore_buffer',JOIN_SPECS['axial_restore'][4])
JOIN_SPECS['actual_patch']=(*JOIN_SPECS['axial_restore'][:3],'buffer_patch',JOIN_SPECS['axial_restore'][4])
JOIN_SPECS['O2_inlet']=('current_O2_background_tensor',
    ('current_actual_O2_reference_slope_axial_buffer_full_tensors_available','current_actual_O2_full_meridional_decomposition_available','current_actual_four_O2_completed_tensor_joins_certified'),
    'current_actual_four_O2_tensor_interface_rows','reference_slope',
    ('current_actual_O2_source_and_tensor_theorem','source_function_equality_precedes_common_tensor_bounds'))


def original_join_receipts(family):
    result={}
    for chart,(stem,flags,rows,interface,path) in JOIN_SPECS.items():
        name=PREFIX+stem+'_check.json';receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or any(receipt.get(k)!=v for k,v in family.items()) or not all(receipt.get(k) for k in flags):
            raise ValueError('Same accepted original background/history/P0 join receipt required: '+chart)
        if interface not in receipt[rows]:raise ValueError('Original completed source interface absent: '+interface)
        proof=receipt
        for key in path:proof=proof[key]
        if proof is not True:raise ValueError('Same original source function proof required: '+chart)
        result[chart]=dict(receipt=name,sha256=sha(name),source_family=family,
            required_original_completed_join_flags=list(flags),original_interface_row=rows+'.'+interface,
            original_source_function_identity_flag_path=list(path),source_identity_flag_is_true=True,
            background_history_and_analytic_P0_identity_not_inferred_from_radius_only=True)
    return result


class NativeMiddleO2InletC1Histories:
    def __init__(self,owner):
        if type(owner) is not preceding.NativeBridgeSwitchC1Histories:
            raise ValueError('Accepted actual inlet-to-R110 C1 owner required')
        receipt=json.loads((HERE/preceding.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[preceding.GATE] or receipt['source_family']!=owner.family:
            raise ValueError('Same checked inlet-to-R110 C1 source required')
        self.owner=owner;self.transfer=owner.transfer;self.ctx=owner.ctx;self.family=owner.family
        self.coordinates=owner.coordinates;self.service=owner.service;self.q_owner=owner.q_owner
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({name:sha(name) for name in (preceding.RECEIPT,Path(__file__).name,
            PREFIX+'current_reshape_background_tensor.py',PREFIX+'current_restore_background_tensor.py',
            PREFIX+'current_actual_patch_background_tensor.py',PREFIX+'current_O2_background_tensor.py')})
        self.joins=original_join_receipts(self.family)
        self.service.bind_hashes({v['receipt']:v['sha256'] for v in self.joins.values()})

    def _cell(self,chart,left,right,Z,N,incoming,previous_chart,previous_endpoint):
        if (chart,left,right) not in ROUTE:raise ValueError('Original complete middle/reference chart required')
        c=self.ctx
        if chart=='actual_patch':
            endpoint=c.exp(1)
            coordinate=c.mpf((1,ep(endpoint)[1]))
            geometry=self.transfer.geometry.build(chart,1,{},coordinate,['1','original exp(1)'])
        else:
            endpoint=c.mpf(right);geometry=self.transfer.geometry.cell(chart,left,right)
        seam=previous_chart+' -> '+chart
        identity=geometry['record']['exact_original_radius_Jacobian_identities']
        if not identity['passed'] or seam not in identity['original_same_radius_periodic_phase_seam_identities']:
            raise ValueError('Accepted exact neighboring original radius/phase seam required')
        source=self.q_owner.query(chart,Z,geometry['coordinate']);roots=source['source']['roots']
        root_owner=self.q_owner.owner.owner;signed=self.transfer.owner.signed_owner
        positive=root_owner.decode(root_owner.inventory[chart]['actual_positive_denominator_theorem'])
        if not positive['source_function_positivity_not_inferred_from_saved_box']:
            raise ValueError('Original analytic whole-chart a-positive lower required')
        eta=packets.interval(c,root_owner.scales['selected_positive_eta_log'])
        dstar=packets.interval(c,root_owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
        primitives=p.whole_period_C1(roots,eta,positive['log_actual_a_positive_lower'],dstar)
        E,E_Z=roots['E'][ZERO],roots['E'][DZ];packet=source['source']['packet']
        def axial(k):
            row=prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
            return signed.leaf(prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),E.scale.bases,E.ledger)
        kernels=density.density_Z_kernels(E,E_Z,axial(0),axial(1),primitives['values'],N)
        factors={key:transfer.true_width_kernel(self.coordinates,geometry,rate) for key,rate in RATES.items()}
        values={key:self.coordinates.rebase(kernels['kernels'][key],self.family)*factors[key]['mass'] for key in RATES}
        jets={key:self.coordinates.rebase(kernels['Z_derivatives'][key],self.family)*factors[key]['mass'] for key in RATES}
        operator=history.C1DuhamelOperator(self.coordinates)
        transfer.append_true_cell(operator,geometry,values,jets,self.family)
        correction=operator.apply(incoming['values'],incoming['Z_derivatives'],self.family)
        background=history.packet_history_functions(self.transfer.owner.native.query(chart,Z,endpoint),self.coordinates,signed)
        own={key:background['originals'][key]+correction['values'][key] for key in RATES}
        own_Z={key:background['Z_derivatives'][key]+correction['Z_derivatives'][key] for key in RATES}
        record=dict(source_family=self.family,chart=chart,Z_box=c.mpf(Z),candidate_N=N,
            original_left_endpoint=dict(chart=previous_chart,coordinate=previous_endpoint),
            original_right_endpoint=dict(chart=chart,coordinate=right,analytic_source_expression='exp(1)' if right=='e' else str(right)),
            exact_original_neighbor_radius_and_periodic_phase_seam=seam,
            original_background_history_P0_source_join_receipt=self.joins[chart],
            original_actual_patch_direct_source_receipt=self.joins['Rh_reference'] if chart=='actual_patch' else None,
            actual_whole_native_chart_geometry=geometry['record'],
            original_chart_uniform_a_positive_certificate=positive,
            original_full_box_signed_C1_source=source['record'],original_whole_period_C1_cover=primitives['record'],
            original_true_width_kernel_factors={k:dict(branch=v['branch'],decay=v['decay'].record(),mass=v['mass'].record()) for k,v in factors.items()},
            true_log_radius_signed_C0_contributions={k:v.record() for k,v in values.items()},
            true_log_radius_signed_Z_contributions={k:v.record() for k,v in jets.items()},
            actual_inherited_correction_C0={k:v.record() for k,v in incoming['values'].items()},
            actual_inherited_correction_Z={k:v.record() for k,v in incoming['Z_derivatives'].items()},
            actual_right_correction_C0={k:v.record() for k,v in correction['values'].items()},
            actual_right_correction_Z={k:v.record() for k,v in correction['Z_derivatives'].items()},
            actual_right_original_background_and_separate_P0_Z=background['record'],
            actual_right_own_history_C0={k:v.record() for k,v in own.items()},
            actual_right_own_history_Z={k:v.record() for k,v in own_Z.items()},
            actual_cell_C1_operator=operator.record(),
            original_source_rows_are_ordinary_log_radius_and_Z=True,
            original_true_width_and_source_Jacobian_applied_once=True,
            actual_incoming_correction_not_reset_or_background_double_added=True,
            source_cutoff_or_phase_crossings_not_deleted=True,
            conservative_C1_covers_do_not_prove_tight_signed_error=True,
            global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,geometry=geometry,source=source,primitives=primitives,factors=factors,
            incoming=incoming,operator=operator,contributions=values,Z_derivatives=jets,
            correction=correction,background=background,own=own,own_Z=own_Z)

    @native.inlet.source_precision
    def route(self,Z=(-1,1),N=1024):
        N=density.spatial.density.candidate_integer(N)
        if N<160:raise ValueError('Accepted whole-period primitive cover requires candidate N>=160')
        initial=self.owner.route(Z,N=N);incoming=initial['correction']
        previous_chart,previous_endpoint='switch_power',1
        cells={};cumulative=history.C1DuhamelOperator(self.coordinates)
        for chart,left,right in ROUTE:
            cell=self._cell(chart,left,right,Z,N,incoming,previous_chart,previous_endpoint)
            transfer.append_true_cell(cumulative,cell['geometry'],cell['contributions'],cell['Z_derivatives'],self.family)
            cell['record']['R110_to_current_endpoint_C1_operator']=cumulative.record()
            cells[chart]=cell;incoming=cell['correction'];previous_chart,previous_endpoint=chart,right
        final=cells[ROUTE[-1][0]]
        terminal_seam='Rh_reference -> O2_slope'
        if terminal_seam not in self.transfer.geometry.binder.identity['original_same_radius_periodic_phase_seam_identities']:
            raise ValueError('Original reference-to-O2-inlet radius/phase seam required')
        background=history.packet_history_functions(self.transfer.owner.native.query('O2_slope',Z,0),
            self.coordinates,self.transfer.owner.signed_owner)
        own={key:background['originals'][key]+incoming['values'][key] for key in RATES}
        own_Z={key:background['Z_derivatives'][key]+incoming['Z_derivatives'][key] for key in RATES}
        record=dict(source_family=self.family,Z_box=self.ctx.mpf(Z),candidate_N=N,
            known_actual_inlet_to_R110_C1_history=initial['record'],
            ordered_original_native_route=[dict(chart=chart,left=left,right=right) for chart,left,right in ROUTE],
            actual_serial_chart_C1_history_records={chart:cell['record'] for chart,cell in cells.items()},
            R110_to_O2_inlet_C1_operator=cumulative.record(),
            actual_original_inlet_to_O2_inlet_correction_C0={k:v.record() for k,v in incoming['values'].items()},
            actual_original_inlet_to_O2_inlet_correction_Z={k:v.record() for k,v in incoming['Z_derivatives'].items()},
            original_terminal_reference_slope_radius_and_phase_seam=terminal_seam,
            original_terminal_reference_slope_background_history_P0_join_receipt=self.joins['O2_inlet'],
            actual_O2_slope0_background_directly_queried=True,
            actual_O2_inlet_original_background_and_separate_P0_Z=background['record'],
            actual_O2_inlet_own_history_C0={k:v.record() for k,v in own.items()},
            actual_O2_inlet_own_history_Z={k:v.record() for k,v in own_Z.items()},
            no_original_interval_skipped_between_true_inlet_and_O2_inlet=True,
            actual_R110_incoming_C0_Z_functions_carried_through_every_chart=True,
            actual_patch_exp1_endpoint_and_exact_unit_log_width_retained=True,
            reusable_R110_to_O2_inlet_affine_operator_installed=True,
            original_O2_inlet_to_Rc_route_still_missing=True,
            conservative_covers_do_not_prove_terminal_control_or_common_N=True,
            global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,initial=initial,cells=cells,cumulative=cumulative,correction=incoming,
            own=own,own_Z=own_Z,background=background)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        c1=density.NativeDensityC1LocalIntegrals(p.first.NativePhaseFirstJets(p.slow.NativeQSlowJets(
            p.current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge))))))
        second=preceding.preceding.NativeSecondBridgeC1Histories(p.NativeFirstBridgeC1Histories(
            transfer.NativeTrueChartC1Transfer(history.NativeC1HistoryTransfer(c1))))
        owner=NativeMiddleO2InletC1Histories(preceding.NativeBridgeSwitchC1Histories(second))
        records={}
        for name,Z in (('whole_Z',(-1,1)),('Z_interval',('.49','.51'))):
            result=owner.route(Z,N=1024);records[name]=result['record']
            print('Actual inlet-to-O2_inlet middle/reference C1 route:',name,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},candidate_N=1024,native_Z_query_count=2,
        complete_added_original_chart_count=len(ROUTE),actual_inlet_to_O2_inlet_C1_history_records=records,
        actual_original_inlet_to_O2_inlet_C1_histories_installed=True,
        R110_correction_memory_retained_across_all_six_new_charts=True,
        wide_signed_derivative_covers_need_quantitative_refinement=True,
        global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Actual original inlet sc/2 -> accepted complete inner R110 route -> reshape/reference/axial restoration/buffer/actual patch/Rh_reference -> original O2 inlet (reference coordinate0), with five correction/own C0/Z functions on whole Z[-1,1] and[.49,.51]. Exact original widths/seams, analytic exp(1) patch endpoint and separate backgrounds/P0_Z. Conservative covers; O2/O3-to-Rc/terminal/common N/cone/recursion/full NS remain incomplete.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
