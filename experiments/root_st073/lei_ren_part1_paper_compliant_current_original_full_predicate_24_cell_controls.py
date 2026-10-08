"""Actual13+five exact full-predicate O2+six Rc candidate control ranges.

Accepted source-function operator covers are restored into one live target
owner and applied in original order. Finite Picard and residual enclosures
do not imply a fixed point, global frequency, or functional terminal closure.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rc_full_predicate_joint_targets as joint
import lei_ren_part1_paper_compliant_current_original_O2_five_slope_partitions as partitions
import lei_ren_part1_paper_compliant_current_native_Rc_functional_controls as controls

preceding=joint.preceding;driver=preceding.driver;native=joint.native;history=preceding.history
parameters=joint.parameters;packets=joint.packets;fixed=joint.fixed
HERE,PREFIX,sha=joint.HERE,joint.PREFIX,joint.sha;ep=joint.ep;encode=joint.encode;KEYS=joint.KEYS;N=1024
NAME=PREFIX+'current_original_full_predicate_24_cell_controls.json.gz'
RECEIPT=PREFIX+'current_original_full_predicate_24_cell_controls_check.json'
GATE='actual_original_full_predicate_24_cell_N1024_histories_and_finite_control_residual_ranges_installed'


def records(group):return {key:value.record() for key,value in group.items()}


def c1records(rows,labels):return {key:dict(C0=row.value.record(),Z=row.Z.record()) for key,row in zip(labels,rows,strict=True)}


class OriginalFullPredicate24CellControls:
    def __init__(self,bridge=None):
        self.joint=joint.OriginalRcFullPredicateJointTargets(bridge)
        self.target_owner=self.joint.target_owner;self.ctx=c=self.target_owner.ctx
        self.coordinates=self.target_owner.coordinates;self.family=self.target_owner.family;self.hashes=dict(self.joint.hashes)
        self.five_manifest=driver.accepted(self.hashes,partitions,self.family)
        self.five=self.five_manifest['actual_original_O2_five_slope_partitions']
        self.rc_manifest=self.joint.manifest
        self.rc=self.rc_manifest['actual_full_predicate_O2_correction_to_Rc_continuation']
        middle=driver.middle
        self.middle_manifest=driver.accepted(self.hashes,middle,self.family)
        for module in (middle.preceding,middle.preceding.preceding,middle.p):
            driver.accepted(self.hashes,module,self.family)
        common=self.five['accepted_actual_N1024_density_and_upstream_binding']['original_common_coordinates']
        if common['source_family']!=self.family or self.five['candidate_N']!=N or self.five['exact_Z_range']!=['-1','1']:
            raise ValueError('Same actual original whole-Z N1024 five-slope source required')
        for live,saved in zip(self.coordinates.bases,common['common_log_bases'],strict=True):
            if ep(live)!=ep(packets.interval(c,saved)):raise ValueError('Actual five-slope and target factors differ')
        w=self.middle_manifest['actual_inlet_to_O2_inlet_C1_history_records']['whole_Z']
        r110=w['known_actual_inlet_to_R110_C1_history'];phase2=r110['known_actual_inlet_to_phase2_C1_history']
        first=phase2['known_original_first_bridge_phase1_history'];initial=first['known_original_initial_collar_history']
        if initial['exact_collar_fraction_endpoints']!=['1/2','3/4'] or not initial['original_correction_initial_condition_is_the_checked_sc_half_inlet']:
            raise ValueError('Actual checked sc-half source inlet and flat collar required')
        self.prefix=[]
        def add(label,chart,raw,geometry_key,operator_key,increment_keys,background_key,output_keys,path):
            self.require_identity(raw)
            if raw[geometry_key]['chart']!=chart:raise ValueError('Original source cell chart changed')
            operator=raw[operator_key] if operator_key else None
            if increment_keys is None:
                increments=[operator['cumulative_signed_increment_C0_enclosures'],operator['cumulative_signed_increment_Z_enclosures']]
            else:increments=[raw[key] for key in increment_keys]
            self.prefix.append(dict(label=label,chart=chart,raw=raw,geometry=raw[geometry_key],operator=operator,
                increments=increments,background=raw[background_key],saved_out=[raw[key] for key in output_keys],
                source_binding=dict(manifest=middle.NAME,manifest_sha256=sha(middle.NAME),record_path=path)))
        basepath=['actual_inlet_to_O2_inlet_C1_history_records','whole_Z']
        rpath=basepath+['known_actual_inlet_to_R110_C1_history'];spath=rpath+['known_actual_inlet_to_phase2_C1_history']
        fpath=spath+['known_original_first_bridge_phase1_history'];ipath=fpath+['known_original_initial_collar_history']
        add('initial_flat_collar','bridge_first',initial,'actual_true_chart_cell_geometry','true_width_affine_operator',None,
            'original_right_background_and_separate_P0_Z',('actual_correction_C0_enclosures','actual_correction_Z_enclosures'),ipath)
        add('active_first_bridge','bridge_first',first,'actual_entire_remaining_first_bridge_geometry',None,
            ('true_log_radius_signed_C0_contributions','true_log_radius_signed_Z_contributions'),
            'original_phase1_background_and_separate_P0_Z',('actual_first_bridge_exit_correction_C0','actual_first_bridge_exit_correction_Z'),fpath)
        add('second_bridge','bridge_second',phase2,'actual_whole_second_bridge_geometry','actual_second_bridge_C1_operator',
            ('true_log_radius_signed_C0_contributions','true_log_radius_signed_Z_contributions'),
            'original_phase2_background_and_separate_P0_Z',('actual_phase2_correction_C0','actual_phase2_correction_Z'),spath)
        for parent,path in ((r110,rpath),(w,basepath)):
            for chart,raw in parent['actual_serial_chart_C1_history_records'].items():
                if raw['chart']!=chart:raise ValueError('Original upstream chart key/source mismatch')
                add(chart,chart,raw,'actual_whole_native_chart_geometry','actual_cell_C1_operator',
                    ('true_log_radius_signed_C0_contributions','true_log_radius_signed_Z_contributions'),
                    'actual_right_original_background_and_separate_P0_Z',('actual_right_correction_C0','actual_right_correction_Z'),
                    path+['actual_serial_chart_C1_history_records',chart])
        if [(row['label'],row['chart']) for row in self.prefix]!=[(row[0],row[1]) for row in parameters.ROUTE[:13]]:
            raise ValueError('Actual original13 prefix order required')
        self.upstream_saved=[w['actual_original_inlet_to_O2_inlet_correction_C0'],w['actual_original_inlet_to_O2_inlet_correction_Z']]
        self.five_rows=self.five['actual_original_five_slope_cell_records']
        self.six_rows=self.rc['actual_six_downstream_source_C0_Z_cell_records']
        if [(row['label'],row['chart']) for row in self.five_rows]!=[(row[0],row[1]) for row in parameters.ROUTE[13:18]]:
            raise ValueError('Actual original five-slope source order required')
        if [(row['ordered_continuation_label'],row['chart']) for row in self.six_rows]!=[(row[0],row[1]) for row in parameters.ROUTE[18:]]:
            raise ValueError('Actual original six downstream source order required')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.source_ranges=None

    def require_identity(self,record):
        if record['source_family']!=self.family or record['candidate_N']!=N or ep(packets.interval(self.ctx,record['Z_box']))!=(-1,1):
            raise ValueError('Same actual source family,N1024 and whole-Z cell required')

    def restore(self,record):return driver.upstream.restore_common_source(record,self.coordinates)

    def restore_group(self,group):
        if set(group)!=set(KEYS):raise ValueError('All five original source rows required')
        return {key:self.restore(value) for key,value in group.items()}

    def restore_operator(self,record):
        if record['steps']!=1 or record['original_rates']!={key:str(rate) for key,rate in history.RATES.items()}:
            raise ValueError('One source-owned original own-rate cell operator required')
        operator=history.C1DuhamelOperator(self.coordinates)
        operator.coefficients=self.restore_group(record['incoming_C0_Z_decay_coefficients'])
        operator.increments=self.restore_group(record['cumulative_signed_increment_C0_enclosures'])
        operator.Z_increments=self.restore_group(record['cumulative_signed_increment_Z_enclosures']);operator.steps=1
        return operator

    def restored_background(self,record):
        return dict(values=self.restore_group(record['original_normalized_history_C0_enclosures']),
            Z_derivatives=self.restore_group(record['original_normalized_history_Z_enclosures']),
            P0=self.restore(record['original_separate_P0_over_Pstar_squared']),P0_Z=self.restore(record['original_separate_P0_Z_over_Pstar_squared']))

    @native.inlet.source_precision
    def assemble(self,*,N=1024):
        if type(N) is not int or N!=1024:raise ValueError('Accepted actual24 source ranges require unchanged integer N1024')
        c=self.ctx;coords=self.coordinates;cells=[]
        # This is the checked original source inlet, not a new zero fallback.
        initial=self.prefix[0]
        if any(not self.restore(value).zero for group in initial['saved_out'] for value in group.values()):
            raise ValueError('Original exact-flat inlet theorem must retain zero initial correction')
        incoming=dict(values={key:coords.scalar(0) for key in KEYS},Z_derivatives={key:coords.scalar(0) for key in KEYS})
        with mp.workdps(c.dps+40):
            for index,spec in enumerate(parameters.ROUTE):
                label,chart,left,right=spec
                if index<13:
                    source=self.prefix[index];raw=source['raw'];binding=source['source_binding']
                    increments=[self.restore_group(group) for group in source['increments']]
                    if source['operator'] is None:
                        geometry=self.target_owner.geometry(label,chart,left,right)
                        operator=history.C1DuhamelOperator(coords)
                        operator.coefficients={key:preceding.downstream.transfer.true_width_kernel(coords,geometry,rate)['decay'] for key,rate in history.RATES.items()}
                        operator.increments,operator.Z_increments=increments;operator.steps=1
                        if encode(geometry['record'])!=source['geometry']:
                            raise ValueError('Actual active-first-bridge source geometry changed')
                    else:operator=self.restore_operator(source['operator'])
                    background=self.restored_background(source['background'])
                    source_geometry=source['geometry'];saved_out=source['saved_out']
                elif index<18:
                    raw=self.five_rows[index-13];self.require_identity(raw)
                    operator=self.restore_operator(raw['actual_cell_C1_operator'])
                    increments=[self.restore_group(raw['changed_C0_contributions']),self.restore_group(raw['changed_Z_contributions'])]
                    background=dict(values=self.restore_group(raw['original_right_background_C0']),
                        Z_derivatives=self.restore_group(raw['original_right_background_Z']),
                        P0=self.restore(raw['original_separate_P0']),P0_Z=self.restore(raw['original_separate_P0_Z']))
                    binding=dict(manifest=partitions.NAME,manifest_sha256=sha(partitions.NAME),
                        record_path=['actual_original_O2_five_slope_partitions','actual_original_five_slope_cell_records',index-13])
                    source_geometry=raw['actual_true_cell_geometry'];saved_out=[raw['actual_right_correction_C0'],raw['actual_right_correction_Z']]
                else:
                    raw=self.six_rows[index-18];self.require_identity(raw)
                    operator=self.restore_operator(raw['actual_cell_C1_operator'])
                    increments=[self.restore_group(raw['actual_true_width_C0_contributions']),self.restore_group(raw['actual_true_width_Z_contributions'])]
                    background=self.restored_background(raw['original_right_background_and_separate_P0_Z'])
                    binding=dict(manifest=preceding.NAME,manifest_sha256=sha(preceding.NAME),
                        record_path=['actual_full_predicate_O2_correction_to_Rc_continuation','actual_six_downstream_source_C0_Z_cell_records',index-18])
                    source_geometry=raw['actual_true_cell_geometry'];saved_out=[raw['actual_right_correction_C0'],raw['actual_right_correction_Z']]
                correction=operator.apply(incoming['values'],incoming['Z_derivatives'],self.family)
                own={key:background['values'][key]+correction['values'][key] for key in KEYS}
                ownZ={key:background['Z_derivatives'][key]+correction['Z_derivatives'][key] for key in KEYS}
                record=dict(label=label,chart=chart,candidate_N=N,source_family=self.family,Z_box=c.mpf((-1,1)),
                    actual_original_source_cell_binding=binding,actual_original_geometry=source_geometry,
                    actual_cell_C1_operator=operator.record(),actual_local_C0_contributions=records(increments[0]),actual_local_Z_contributions=records(increments[1]),
                    actual_inherited_correction_C0=records(incoming['values']),actual_inherited_correction_Z=records(incoming['Z_derivatives']),
                    actual_right_correction_C0=records(correction['values']),actual_right_correction_Z=records(correction['Z_derivatives']),
                    actual_right_background_C0=records(background['values']),actual_right_background_Z=records(background['Z_derivatives']),
                    actual_right_own_history_C0=records(own),actual_right_own_history_Z=records(ownZ),
                    original_separate_P0=background['P0'].record(),original_separate_P0_Z=background['P0_Z'].record(),
                    absolute_pressure_C0=(background['P0']+own['p']).record(),absolute_pressure_Z=(background['P0_Z']+ownZ['p']).record(),
                    accepted_entire_source_operator_covers_not_numerical_source_owners_restored=True,
                    actual_new_incoming_applied_not_old_outgoing_hydrated=True,
                    no_integrated_increment_mass_or_original_background_double_addition=True,
                    original_candidate_N_and_true_geometry_and_phase_unchanged=True,
                    quiet_power_local_zero_retains_actual_upstream_correction=True)
                cells.append(dict(record=record,operator=operator,values=increments[0],Z_derivatives=increments[1],
                    incoming=incoming,correction=correction,background=background,saved_ancestor_out=saved_out))
                incoming=correction
                print('Actual original24 source operator:',index+1,label,flush=True)
            if len(cells)!=24 or [(row['record']['label'],row['record']['chart']) for row in cells]!=[(row[0],row[1]) for row in parameters.ROUTE]:
                raise ValueError('Actual complete original24 source range route required')
            amplitude=self.joint.amplitude_source()
            target=dict(**fixed.fixed_N_target_rows(incoming['values'],incoming['Z_derivatives'],amplitude['A'],amplitude['AZ'],
                amplitude['logA'],amplitude['mu_source'],amplitude['logmu']),
                **{key:amplitude[key] for key in ('A','AZ','mu','logA','logmu')})
            record=dict(source_family=self.family,candidate_N=N,Z_box=c.mpf((-1,1)),
                original_exact_function_route_cells=24,actual_source_prefix_slope_downstream_counts=[13,5,6],
                actual_original24_source_cell_records=[row['record'] for row in cells],original_common_coordinates=coords.record(),
                actual_Rc_correction_C0_Z=[records(incoming['values']),records(incoming['Z_derivatives'])],
                actual_Rc_joint_target_C0_Z=[records(target['values']),records(target['Z_derivatives'])],
                actual_Rc_joint_numerator_C0_Z=[target['joint_numerator'].record(),target['joint_numerator_Z'].record()],
                original_Rc_amplitude=amplitude['record'],
                full24_original_C1_integral_range_transport_enclosed=True,
                same_original24_function_definitions_and_source_radius_phase_not_fixtures=True,
                source_operator_covers_reapplied_without_upstream_or_O2_or_Rc_source_query_replay=True,
                original_P0_and_background_outside_correction_target=True,
                spatial_O2_global_alternative_not_summed_into_direct24_history=True,
                actual_five_controls_installed=False,functional_terminal_identity_solved=False,current_whole_N_selected=False)
            self.source_ranges=dict(record=record,cells=cells,history=incoming['values'],Z_derivatives=incoming['Z_derivatives'],target=target,amplitude=amplitude)
        self.hashes.update(self.joint.hashes);self.hashes.update(self.target_owner.service.hashes)
        return self.source_ranges

    @native.inlet.source_precision
    def finite_controls(self,source_ranges,*,iterations=3):
        if source_ranges is not self.source_ranges or source_ranges is None:
            raise ValueError('Issued actual original24 source function ranges required')
        role=controls.roles.RoleBoundNativeRcFunctions(controls.source.NativeRcFunctionTransport(self.target_owner))
        owner=controls.NativeRcFunctionalControls(role);range_owner=controls.paired.NativePairedC1Transport(role)
        live=owner.controls((-1,1),N,range_owner=range_owner,source_ranges=source_ranges,iterations=iterations)
        built=live['built']
        if built['source_family']!=self.family or [(row['label'],row['chart']) for row in built['cells']]!=[(row[0],row[1]) for row in parameters.ROUTE]:
            raise ValueError('Control exact integral functions must have the same actual original24 source route')
        self.hashes.update(owner.hashes);self.hashes.update(range_owner.hashes);self.hashes.update(self.target_owner.service.hashes)
        diagnostics=[];c=self.ctx;matrix=live['exact_matrix_enclosures'];target=source_ranges['target']
        with mp.workdps(c.dps+40):
            for depth,row in enumerate(live['ranges']):
                Q=controls.range_quadratic(c,target['mu'],target['logmu'],matrix,row)
                zero=row[0].value.scalar(0);residual=[]
                for j,d in enumerate(live['N_scaled_target_ranges']):
                    residual.append(controls.C1Enclosure(
                        d.value+sum((q.value*b for q,b in zip(row,matrix['linear_enclosure'][j],strict=True)),zero)+Q[j].value*(c.mpf(1)/N),
                        d.Z+sum((q.Z*b for q,b in zip(row,matrix['linear_enclosure'][j],strict=True)),zero)+Q[j].Z*(c.mpf(1)/N)))
                cap=lambda pairs,labels:{key:parameters.repair.LogUpper.add(c,[parameters.magnitude(q.value),parameters.magnitude(q.Z)]).record()
                    for key,q in zip(labels,pairs,strict=True)}
                diagnostics.append(dict(finite_depth=depth,actual_control_C0_Z_ranges=c1records(row,controls.CONTROLS),
                    actual_repair_equation_residual_C0_Z_ranges=c1records(residual,controls.ROWS),
                    actual_control_C1_upper_caps=cap(row,controls.CONTROLS),actual_residual_C1_upper_caps=cap(residual,controls.ROWS)))
            contraction=controls.fixed_N_contraction_diagnostic(live)
        report=dict(candidate_N=N,source_family=self.family,exact_Z_range=['-1','1'],finite_picard_depth=iterations,
            actual_finite_control_and_residual_ranges=diagnostics,
            actual_N_scaled_target_C0_Z_ranges=c1records(live['N_scaled_target_ranges'],controls.ROWS),
            original_exact_weight_enclosures=live['exact_weight_enclosures'],original_exact_matrix_enclosures=matrix,
            fixed_N_sufficient_contraction_diagnostic=contraction,
            exact_control_graph_nodes=built['graph'].nodes,exact_original_source_roles=built['source_roles'],
            exact_function_source_graph_sha256=built['source_graph_sha256'],
            exact_finite_control_function_roots=[controls.pair_roots(row,controls.CONTROLS) for row in built['finite_picard_sequence']],
            exact_repair_residual_function_roots=controls.pair_roots(built['control_residual'],controls.ROWS),
            numerical_residual_definition='N*r + B_exact*h + Q_exact(h)/N; same exact first-Z product rules',
            same_live_role_target_paired_range_owner=True,source_input_is_actual_original24_continuous_integral_covers=True,
            exact_source_functions_not_defined_by_enclosure_endpoints=True,
            finite_iterates_not_solved_controls=True,certified_fixed_point_tail_installed=False,
            actual_five_controls_installed=False,functional_terminal_identity_solved=False,current_whole_N_selected=False)
        return dict(report=report,live=live,control_owner=owner,role_owner=role,range_owner=range_owner)


@native.inlet.source_precision
def run():
    began=time.monotonic();bridge,construction=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime() as runtime:
        owner=OriginalFullPredicate24CellControls(bridge);source=owner.assemble();result=owner.finite_controls(source)
    report=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
        actual_original24_full_predicate_source_range_bridge=source['record'],actual_original24_finite_control_diagnostics=result['report'],
        original_bridge_construction=construction,original_source_runtime=runtime.record(),
        full24_original_C1_integral_range_transport_enclosed=True,actual_finite_picard_and_residual_C0_Z_ranges_installed=True,
        actual_five_controls_installed=False,functional_terminal_identity_solved=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual original13+five exact full-predicate O2+six downstream source operators compose atN1024 over fullZ and bind the same original Rc joint target to finite control/residual enclosures. Separate backgrounds/P0 and full source errors. Fixed-point tail, terminal/global N/sharp source/heat/cone/recursion/full corrected NS remain open.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(encode(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Actual original full-predicate24 source histories and finite control residuals generated',flush=True);return report


if __name__=='__main__':run()
