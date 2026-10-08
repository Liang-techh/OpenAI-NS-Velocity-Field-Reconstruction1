"""Actual full-predicate O2 correction continues through six cells to Rc.

Accepted inherited source covers are restored losslessly. The unchanged
actual downstream sources execute afresh; historical route entrypoints do
not run. Direct and spatial inlet covers remain separate alternatives.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_full_predicate_incoming_driver as driver
import lei_ren_part1_paper_compliant_current_native_Rc_C1_histories as downstream
import lei_ren_part1_paper_compliant_current_original_O2_source_parameter_frame as parameters

base,prior,ep=driver.base,driver.prior,driver.ep
encode=downstream.packets.encode
native=downstream.native;serial=downstream.serial;history=downstream.history
HERE,PREFIX,sha=driver.HERE,driver.PREFIX,driver.sha
NAME=PREFIX+'current_original_O2_full_predicate_Rc_continuation.json.gz'
RECEIPT=PREFIX+'current_original_O2_full_predicate_Rc_continuation_check.json'
GATE='actual_full_predicate_O2_N1024_correction_continued_through_six_original_cells_to_Rc'
KEYS=driver.KEYS;N=1024
ROUTE=(('axial','O2_axial',0,1),('buffer_0_9','O2_buffer',0,9),('buffer_9_11','O2_buffer',9,11),
    ('transition','O3_slope_mu',0,1),('power_to_r_plus','O3_power',0,1),('power_to_Rc','O3_power',1,2))


class OriginalO2FullPredicateRcContinuation:
    def __init__(self,bridge=None):
        if bridge is None:bridge,_=native.inlet.native_bridge_owner()
        self.shell=downstream.NativeRcC1Histories(downstream.preceding.NativeO2C1Histories(
            downstream.preceding.make_middle_owner(bridge)))
        self.ctx=c=self.shell.ctx;self.family=self.shell.family;self.coordinates=self.shell.coordinates
        self.hashes=dict(self.shell.service.hashes)
        self.manifest=driver.accepted(self.hashes,driver,self.family)
        self.saved=self.manifest['actual_full_predicate_original_O2_source_incoming_driver']
        if self.manifest['candidate_N']!=N or self.saved['explicit_candidate_N']!=N:
            raise ValueError('Same accepted actual O2 candidate N1024 required')
        if self.saved['exact_y_window']!=['0','1'] or self.saved['exact_Z_range']!=['-1','1']:
            raise ValueError('Accepted complete O2 slope and whole-Z inlet required')
        if self.saved['normalized_own_units']!=driver.integrals.five.UNITS or self.saved['original_P0_datum_sha256']!=self.family['datum_enclosure_sha256']:
            raise ValueError('Same original five units and analytic P0 required')
        if not self.saved['original_background_added_once_and_analytic_P0_separate'] or not self.saved['exact_rate_zero_pressure_keeps_all_upstream_memory']:
            raise ValueError('Accepted actual correction/background/pressure split required')
        frame=parameters.OriginalO2SourceParameterFrame()
        self.atlas=driver.spatial.source.O2MixedAtlas(frame,lower=-1,upper=1,logq=c.mpf(0))
        actual=self.atlas.record();old=self.saved['actual_upstream_source_binding']['output_original_source_atlas']
        if actual['source_family']!=self.family or actual['source_basis_order']!=old['source_basis_order'] or actual['exact_Z_bounds']!=old['exact_Z_bounds']:
            raise ValueError('Same original O2 correction basis contract required')
        for live,saved in zip(self.atlas.bases,old['defining_basis'],strict=True):
            if ep(live)!=ep(driver.integrals.interval(self.atlas.ctx,saved)):
                raise ValueError('Accepted original correction defining basis changed')
        if s.simplify(frame.definitions['logPstar']-(s.exp(40)+11))!=0:
            raise ValueError('Exact original shared Pstar definition required')
        seed=self.shell.transfer.owner.native.seed
        seed_family=dict(zip(downstream.packets.FAMILY_KEYS,(seed.family,seed.source,seed.datum_sha)))
        if seed_family!=self.family or ep(c.mpf(seed.params.Md))!=(40,40) or seed.params is not seed.pre.datum.parameters:
            raise ValueError('Same actual original pressure-parameter owner required')
        if self.saved['actual_upstream_source_binding']['original_Pstar_definition']!='logPstar=exp(40)+11':
            raise ValueError('Accepted incoming exact original Pstar contract required')
        if ep(self.coordinates.logP_squared)!=ep(c.mpf(ep(seed.logP*2))):
            raise ValueError('Common factor must be twice the live original parameter-owner logP')
        nativeP2=2*self.atlas.bases[0];commonP2=self.coordinates.logP_squared
        if max(ep(nativeP2)[0],ep(commonP2)[0])>min(ep(nativeP2)[1],ep(commonP2)[1]):
            raise ValueError('Original shared Pstar-squared source disagrees')
        self.initial_native=[];self.initial=[]
        for name in ('propagated_actual_direct_C0_Z_corrections','propagated_actual_spatial_C0_Z_alternative_corrections'):
            groups=self.saved[name]
            if len(groups)!=2 or any(set(group)!=set(KEYS) for group in groups):
                raise ValueError('All five actual correction C0/Z covers required')
            restored=[{key:driver.integrals.restore_value(self.atlas,group[key]) for key in KEYS} for group in groups]
            self.initial_native.append(restored)
            self.initial.append(dict(values={key:self.native_to_common(value) for key,value in restored[0].items()},
                Z_derivatives={key:self.native_to_common(value) for key,value in restored[1].items()}))
        self.inlet_binding=dict(source_manifest=driver.NAME,source_manifest_sha256=sha(driver.NAME),
            source_receipt=driver.RECEIPT,source_receipt_sha256=sha(driver.RECEIPT),
            source_record_keys=['propagated_actual_direct_C0_Z_corrections','propagated_actual_spatial_C0_Z_alternative_corrections'],
            accepted_entire_source_function_covers_restored_not_selected_field_points=True,
            inherited_corrections_not_original_plus_correction_histories=True,
            original_source_family=self.family,candidate_N=N,exact_Z_range=['-1','1'],
            original_upstream_and_O2_producers_not_replayed=True,
            native_factor_identity='p*logP+c*logC+l*logL+offset = (p/2)*log(Pstar^2)+(offset+c*logC+l*logL)',
            exact_Pstar_definition='logPstar=exp(40)+11',
            live_common_factor_definition='CommonSourceCoordinates(c, actual_original_seed.logP*2, same family)',
            same_live_original_pressure_parameter_owner=True,
            native_common_source_constructor=PREFIX+'current_native_C1_history_transfer.py',
            native_common_source_constructor_sha256=sha(PREFIX+'current_native_C1_history_transfer.py'),
            original_native_atlas=self.atlas.record(),common_coordinates=self.coordinates.record(),
            source_C_L_cover_dependence_can_widen_common_offset_not_sharpen_source_errors=True,
            ordinary_Z_source_rows_rebased_not_differentiated_again=True)
        self.hashes.update(frame.hashes)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def native_to_common(self,value):
        if type(value) is not prior.ScaledEnclosure or value.scale.bases is not self.atlas.bases or value.ledger is not self.atlas.ledger:
            raise ValueError('Live restored same-source O2 correction cover required')
        p=value.scale.powers
        if p[3]!=0 or p[4]!=0 or p[0]!=int(p[0]):
            raise ValueError('Collected original whole-Z correction factors required')
        c=self.ctx
        offset=c.mpf(ep(value.scale.offset))
        for index in (1,2):
            if p[index]:offset+=c.mpf(ep(self.atlas.bases[index]))*p[index]
        return prior.ScaledEnclosure(prior.FormalScale(self.coordinates.bases,(0,p[0]/2,0,0,0),offset),
            c.mpf(ep(value.coefficient)),self.coordinates.ledger)

    @native.inlet.source_precision
    def integrate(self,*,N=1024):
        if type(N) is not int or N!=1024:raise ValueError('Accepted inherited corrections require unchanged candidate N1024')
        shell=self.shell;c=self.ctx;incoming=self.initial[0];alternative=self.initial[1]
        cumulative=history.C1DuhamelOperator(self.coordinates);cells=[];altcells=[]
        previous_chart,previous_endpoint='O2_slope',1
        joins=downstream.preceding.source_joins(self.family);Tw=shell.transfer.geometry.binder.fixed['Tw']
        if ep(Tw-2)[0]<=0:raise ValueError('Actual original Rc must lie inside the power source chart')
        with mp.workdps(c.dps+40):
            for index,(label,chart,left,right) in enumerate(ROUTE,1):
                if chart=='O3_power':
                    geometry=shell.transfer.geometry.cell(chart,dict(original_power_offset=left),dict(original_power_offset=right))
                    endpoint=c.mpf(right)/Tw
                    left_record=dict(chart=previous_chart,coordinate=previous_endpoint) if left==0 else dict(chart=chart,original_power_offset=left)
                    right_record=dict(chart=chart,original_power_offset=right,original_phase_expression=str(right)+'/Tw')
                else:
                    geometry=shell.transfer.geometry.cell(chart,left,right);endpoint=c.mpf(right)
                    left_record=dict(chart=previous_chart,coordinate=previous_endpoint);right_record=dict(chart=chart,coordinate=right)
                if chart==previous_chart:
                    seam='same original '+chart+' source'
                    join=dict(same_original_chart_and_source_function=True,
                        original_power_offset=left) if chart=='O3_power' else dict(same_original_chart_and_source_function=True,coordinate=left)
                else:
                    seam=previous_chart+' -> '+chart
                    if seam not in shell.transfer.geometry.binder.identity['original_same_radius_periodic_phase_seam_identities']:
                        raise ValueError('Actual neighboring original radius/common-N phase seam required')
                    join=joins[dict(O2_axial='slope_axial',O2_buffer='axial_buffer',O3_slope_mu='buffer_transition').get(chart,'buffer_transition')] if chart!='O3_power' else shell.join
                api=downstream.original_O3_cell if chart.startswith('O3_') else serial.serial_cell
                cell=api(shell,chart=chart,geometry=geometry,endpoint=endpoint,Z=(-1,1),N=N,incoming=incoming,
                    left_record=left_record,right_record=right_record,seam=seam,join=join)
                downstream.transfer.append_true_cell(cumulative,geometry,cell['contributions'],cell['Z_derivatives'],self.family)
                if chart=='O3_power':
                    if set(cell['primitives']['values'])!={'A','A_Z','B_over_Pstar','B_Z_over_Pstar'}:
                        raise ValueError('Both original primitive ordinary-Z rows must be retained')
                    if not all(value.zero for value in cell['primitives']['values'].values()) or not all(value.zero for value in
                        (*cell['kernels']['kernels'].values(),*cell['kernels']['Z_derivatives'].values(),*cell['contributions'].values(),*cell['Z_derivatives'].values())):
                        raise ValueError('Original whole power flat primitive/density C0/Z increments must be exact zero')
                    cell['record'].update(original_whole_power_q_q_Z_density_C0_Z_exact_zero=True,
                        original_canonical_power_excess_source='2*mu; original positive source not rounded a-minus2',
                        original_actual_power_offsets=[left,right],quiet_local_density_preserves_actual_inherited_correction=True)
                alt=cell['operator'].apply(alternative['values'],alternative['Z_derivatives'],self.family)
                altown={key:cell['background']['originals'][key]+alt['values'][key] for key in KEYS}
                altownZ={key:cell['background']['Z_derivatives'][key]+alt['Z_derivatives'][key] for key in KEYS}
                encode=lambda group:{key:value.record() for key,value in group.items()}
                cell['record'].update(ordered_continuation_label=label,ordered_continuation_step=index,
                    O2_slope_exit_to_current_endpoint_C1_operator=cumulative.record(),
                    new_full_predicate_original_O2_inherited_correction_used=True)
                cells.append(cell['record'])
                altcells.append(dict(label=label,candidate_N=N,source_family=self.family,
                    inherited_C0_Z_alternative=[encode(alternative['values']),encode(alternative['Z_derivatives'])],
                    corrected_C0_Z_alternative=[encode(alt['values']),encode(alt['Z_derivatives'])],
                    own_C0_Z_alternative=[encode(altown),encode(altownZ)],
                    identical_actual_source_cell_operator_used=True,alternative_not_added_to_direct_history=True))
                incoming=cell['correction'];alternative=alt;previous_chart,previous_endpoint=chart,right
                print('Actual full-predicate O2 correction continuation:',label,'->',right,flush=True)
            background=cell['background'];P0=background['P0'];P0Z=background['P0_Z']
            pressure=P0+cell['own']['p'];pressureZ=P0Z+cell['own_Z']['p']
            altpressure=P0+altown['p'];altpressureZ=P0Z+altownZ['p']
            composite=[cumulative.apply(start['values'],start['Z_derivatives'],self.family) for start in self.initial]
            result=dict(source_family=self.family,candidate_N=N,exact_Z_range=['-1','1'],
                accepted_full_predicate_O2_correction_binding=self.inlet_binding,
                actual_rehydrated_native_C0_Z_correction_alternatives=[[encode(group) for group in groups] for groups in self.initial_native],
                actual_rebased_common_C0_Z_correction_alternatives=[[encode(group['values']),encode(group['Z_derivatives'])] for group in self.initial],
                ordered_original_six_cell_route=[dict(label=l,chart=ch,left=a,right=b) for l,ch,a,b in ROUTE],
                actual_six_downstream_source_C0_Z_cell_records=cells,spatial_inlet_alternative_cell_records=altcells,
                O2_slope_exit_to_actual_Rc_C1_operator=cumulative.record(),
                direct_and_spatial_alternative_composite_corrections=[[encode(row['values']),encode(row['Z_derivatives'])] for row in composite],
                actual_Rc_direct_correction_C0_Z=[encode(incoming['values']),encode(incoming['Z_derivatives'])],
                actual_Rc_spatial_alternative_correction_C0_Z=[encode(alternative['values']),encode(alternative['Z_derivatives'])],
                actual_Rc_direct_own_history_C0_Z=[encode(cell['own']),encode(cell['own_Z'])],
                actual_Rc_spatial_alternative_own_history_C0_Z=[encode(altown),encode(altownZ)],
                original_Rc_background_and_separate_P0_Z=background['record'],
                absolute_Rc_direct_pressure_C0_Z=[pressure.record(),pressureZ.record()],
                absolute_Rc_spatial_alternative_pressure_C0_Z=[altpressure.record(),altpressureZ.record()],
                original_Rc_endpoint=dict(radius_expression='Rw*exp(2)',original_power_offset=2,original_power_phase_expression='2/Tw',
                    actual_direct_endpoint_phase_cover=2/Tw,original_Rc_not_power_phase1=True),
                original_Rc_reservation=shell.reservation,original_transition_power_source_join=shell.join,
                same_candidate_N1024_source_five_units_and_P0_through_every_cell=True,
                no_old_upstream_or_old_O2_slope_route_recomputed=True,
                no_gap_from_new_full_predicate_O2_slope_exit_to_original_Rc=True,
                true_radius_width_and_ordinary_log_radius_density_applied_once=True,
                all_original_periodic_source_parameters_and_whole_Z_error_covers_retained=True,
                exact_quiet_power_local_increments_do_not_reset_upstream_pressure=True,
                direct_and_spatial_histories_are_alternatives_not_added=True,
                conservative_C0_Z_covers_not_quantitative_terminal_closure=True,
                global_inlet_to_Rc_histories_admitted=False,functional_terminal_identity_solved=False,
                current_whole_N_selected=False)
        self.hashes.update(shell.service.hashes)
        return dict(report=result,correction=incoming,spatial_alternative_correction=alternative,
            own=cell['own'],own_Z=cell['own_Z'],background=background,cumulative=cumulative)


@native.inlet.source_precision
def run():
    began=time.monotonic();bridge,construction=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime() as runtime:
        owner=OriginalO2FullPredicateRcContinuation(bridge);result=owner.integrate()
    report=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
        actual_full_predicate_O2_correction_to_Rc_continuation=result['report'],
        original_source_runtime=runtime.record(),original_bridge_construction=construction,
        actual_source_owned_full_predicate_to_Rc_C0_Z_covers_installed=True,
        actual_all_route_incoming_histories_installed=False,global_inlet_to_Rc_histories_admitted=False,
        functional_terminal_identity_solved=False,current_whole_N_selected=False,
        all_17_chart_or_24_cell_oracle_installed=False,actual_five_controls_installed=False,
        **dict.fromkeys(downstream.packets.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,
        scope='Actual accepted full-y/full-Z original O2 N1024 source correction propagated through axial, both buffer segments, original O3 transition and two quiet power offset cells to Rc=Rw*exp(2). Fresh original downstream source C0/Z queries, native/common basis adapter, actual inherited errors/pressure memory and separate P0. Direct/spatial alternatives. Conservative covers only; terminal controls, global N, full source oracle, higher jets/cone/heat/energy/recursion/corrected NS remain open.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(encode(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Actual original full-predicate O2 correction propagated to Rc',flush=True);return report


if __name__=='__main__':run()
