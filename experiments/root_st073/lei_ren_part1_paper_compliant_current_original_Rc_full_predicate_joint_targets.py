"""Actual fixed-N full-predicate Rc correction targets on the whole Z interval.

Accepted source-function covers are restored losslessly into one live target
owner. The actual original Rc amplitude is queried afresh. These finite-N
targets do not supply the missing original five-slope-cell control records.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_full_predicate_Rc_continuation as preceding
import lei_ren_part1_paper_compliant_current_native_Rc_parameter_targets as parameters
import lei_ren_part1_paper_compliant_current_native_cutoff_range_transport as fixed

native=preceding.native;downstream=preceding.downstream;packets=downstream.packets
HERE,PREFIX,sha=preceding.HERE,preceding.PREFIX,preceding.sha
ep=preceding.ep;encode=preceding.encode;KEYS=preceding.KEYS;N=1024
NAME=PREFIX+'current_original_Rc_full_predicate_joint_targets.json.gz'
RECEIPT=PREFIX+'current_original_Rc_full_predicate_joint_targets_check.json'
GATE='actual_original_full_predicate_Rc_N1024_whole_Z_five_joint_targets_installed'
FIELDS=('actual_Rc_direct_correction_C0_Z','actual_Rc_spatial_alternative_correction_C0_Z')


class OriginalRcFullPredicateJointTargets:
    def __init__(self,bridge=None):
        if bridge is None:bridge,_=native.inlet.native_bridge_owner()
        shell=downstream.NativeRcC1Histories(downstream.preceding.NativeO2C1Histories(
            downstream.preceding.make_middle_owner(bridge)))
        self.target_owner=parameters.NativeRcParameterTargets(shell)
        self.ctx=c=self.target_owner.ctx;self.coordinates=self.target_owner.coordinates
        self.family=self.target_owner.family;self.hashes=dict(shell.service.hashes)
        self.manifest=preceding.driver.accepted(self.hashes,preceding,self.family)
        saved=self.manifest['actual_full_predicate_O2_correction_to_Rc_continuation']
        if self.manifest['candidate_N']!=N or saved['candidate_N']!=N or saved['exact_Z_range']!=['-1','1']:
            raise ValueError('Same accepted actual full-Z candidate N1024 Rc correction required')
        endpoint=saved['original_Rc_endpoint']
        if endpoint['original_power_offset']!=2 or endpoint['original_power_phase_expression']!='2/Tw' or not endpoint['original_Rc_not_power_phase1']:
            raise ValueError('Actual original Rc at power offset2 required')
        binding=saved['accepted_full_predicate_O2_correction_binding']
        if not binding['same_live_original_pressure_parameter_owner'] or binding['exact_Pstar_definition']!='logPstar=exp(40)+11':
            raise ValueError('Accepted correction must retain the actual pressure-parameter factor')
        common=binding['common_coordinates']
        if common['source_family']!=self.family:
            raise ValueError('Same accepted correction common source family required')
        for actual,record in zip(self.coordinates.bases,common['common_log_bases'],strict=True):
            if ep(actual)!=ep(packets.interval(c,record)):
                raise ValueError('Accepted correction and live target defining factors differ')
        self.corrections=[]
        for field in FIELDS:
            groups=saved[field]
            if len(groups)!=2 or any(set(group)!=set(KEYS) for group in groups):
                raise ValueError('All five accepted actual Rc correction C0/Z rows required')
            self.corrections.append(dict(values={key:preceding.driver.upstream.restore_common_source(groups[0][key],self.coordinates) for key in KEYS},
                Z_derivatives={key:preceding.driver.upstream.restore_common_source(groups[1][key],self.coordinates) for key in KEYS}))
        self.binding=dict(source_manifest=preceding.NAME,source_manifest_sha256=sha(preceding.NAME),
            source_receipt=preceding.RECEIPT,source_receipt_sha256=sha(preceding.RECEIPT),source_record_fields=list(FIELDS),
            source_family=self.family,candidate_N=N,exact_Z_range=['-1','1'],
            common_coordinates=self.coordinates.record(),same_actual_original_pressure_parameter_owner=True,
            complete_signed_source_function_covers_restored_without_source_point_selection=True,
            correction_errors_not_original_background_or_absolute_P0_pressure=True,
            upstream_O2_and_six_downstream_source_suites_not_replayed=True,
            direct_and_spatial_correction_covers_are_alternatives_not_added=True)
        self.hashes.update({Path(__file__).name:sha(Path(__file__).name),
            Path(parameters.__file__).name:sha(Path(parameters.__file__).name),
            Path(fixed.__file__).name:sha(Path(fixed.__file__).name)})
        self.amplitude=None

    @native.inlet.source_precision
    def amplitude_source(self):
        target=self.target_owner;c=self.ctx;coords=self.coordinates
        Tw=target.transfer.geometry.binder.fixed['Tw']
        if ep(Tw-2)[0]<=0:raise ValueError('Original Rc must be inside the actual power source chart')
        endpoint=2/Tw
        if ep(endpoint)!=ep(packets.interval(c,self.manifest['actual_full_predicate_O2_correction_to_Rc_continuation']['original_Rc_endpoint']['actual_direct_endpoint_phase_cover'])):
            raise ValueError('Actual Rc source endpoint changed')
        source=target.q_owner.query('O3_power',(-1,1),endpoint)
        roots=source['source']['roots']
        logA=packets.interval(c,target.owner.reservation['positive_Ac_over_S_log_lower'])
        rawA=coords.rebase(roots['E'][parameters.ZERO],self.family)
        AZ=coords.rebase(roots['E'][parameters.DZ],self.family)
        A=rawA.positive_intersection(logA)
        mu=packets.interval(c,target.owner.domain['right_collar_mu'])
        if mu._mpi_!=target.repair_mu._mpi_ or ep(mu)[0]<=0:
            raise ValueError('Same actual strictly positive source and repair mu required')
        cover=target.repair_manifest['new_repair_geometry']['actual_mu_cover_for_constants_only']
        if ep(cover)[0]>ep(mu)[0] or ep(cover)[1]<ep(mu)[1]:
            raise ValueError('Actual source mu outside accepted repair inverse cover')
        logmu=c.mpf(ep(c.ln(mu))[0])
        background=downstream.history.packet_history_functions(source['source']['packet'],coords,target.transfer.owner.signed_owner)
        record=dict(source_family=self.family,candidate_N=N,chart='O3_power',exact_Z_range=['-1','1'],
            original_Rc_power_offset=2,original_Rc_phase_expression='2/Tw',actual_source_endpoint=endpoint,
            actual_original_source_query=source['record'],original_raw_amplitude_C0=rawA.record(),
            original_positive_reserved_amplitude_C0=A.record(),original_amplitude_Z=AZ.record(),
            original_positive_amplitude_log_lower=logA,original_source_mu=mu,original_source_logmu_lower=logmu,
            original_separate_P0=background['P0'].record(),original_separate_P0_Z=background['P0_Z'].record(),
            actual_phase_independent_Rc_E_E_Z_source_roots_used=True,source_E_Z_not_derivative_of_positive_cap=True,
            unmodulated_original_amplitude_not_candidate_changed_E=True,
            exact_source_positive_Ac_over_S_reservation_not_selected_amplitude=True,
            original_source_and_repair_mu_same_interval=True,original_P0_not_in_correction_targets=True)
        self.amplitude=dict(A=A,AZ=AZ,mu=mu,mu_source=coords.scalar(mu),logA=logA,logmu=logmu,record=record)
        self.hashes.update(target.service.hashes)
        return self.amplitude

    def normalise(self,correction,amplitude):
        if not any(correction is issued for issued in self.corrections) or amplitude is not self.amplitude or amplitude is None:
            raise ValueError('Issued actual Rc source correction and amplitude required')
        for group in correction.values():
            if set(group)!=set(KEYS):raise ValueError('All five issued correction rows required')
            for value in group.values():
                if value.ctx is not self.ctx or value.scale.bases is not self.coordinates.bases or value.ledger is not self.coordinates.ledger:
                    raise ValueError('Live same-target source coordinates and ledger required')
        got=fixed.fixed_N_target_rows(correction['values'],correction['Z_derivatives'],amplitude['A'],amplitude['AZ'],
            amplitude['logA'],amplitude['mu_source'],amplitude['logmu'])
        return dict(**got,**{key:amplitude[key] for key in ('A','AZ','mu','logA','logmu')})

    @native.inlet.source_precision
    def targets(self,*,N=1024):
        if type(N) is not int or N!=1024:raise ValueError('Accepted actual correction targets require unchanged integer N1024')
        c=self.ctx
        with mp.workdps(c.dps+40):
            amplitude=self.amplitude_source();targets=[self.normalise(row,amplitude) for row in self.corrections]
            rows=[]
            for label,target in zip(('direct','spatial_alternative'),targets,strict=True):
                caps={key:parameters.repair.LogUpper.add(c,[parameters.magnitude(target['values'][key]*N),
                    parameters.magnitude(target['Z_derivatives'][key]*N)]) for key in parameters.repair.ROWS}
                rows.append(dict(alternative=label,actual_target_C0={key:value.record() for key,value in target['values'].items()},
                    actual_target_Z={key:value.record() for key,value in target['Z_derivatives'].items()},
                    actual_joint_numerator=target['joint_numerator'].record(),actual_joint_numerator_Z=target['joint_numerator_Z'].record(),
                    actual_N_scaled_target_C1_upper_caps={key:value.record() for key,value in caps.items()},
                    norm='sup_Z|N*r_i| + sup_Z|N*r_i_Z|; conservative upper bound',
                    signed_joint_before_division_not_independent_normalized_row_subtraction=True,
                    interval_hulls_do_not_prove_quantitative_joint_cancellation=True))
            report=dict(source_family=self.family,candidate_N=N,exact_Z_range=['-1','1'],
                actual_accepted_Rc_correction_binding=self.binding,actual_original_Rc_amplitude_source=amplitude['record'],
                actual_correction_C0_Z_alternatives=[[{key:value.record() for key,value in row[field].items()}
                    for field in ('values','Z_derivatives')] for row in self.corrections],
                actual_five_joint_target_alternatives=rows,original_target_row_order=parameters.repair.ROWS,
                fixed_N_target_definition='r=(delta_m/A,(delta_k-A*delta_m)/(mu*A^2),delta_h/A,delta_e/A^2,delta_p/A^2)',
                joint_Z_definition='C_Z=delta_k_Z-A_Z*delta_m-A*delta_m_Z',
                full_predicate_axis_and_both_Z_signs_retained=True,
                exact_original_24_cell_control_range_precondition_met=False,
                missing_control_record_dependency='Reintegrate full-predicate O2 at exact original five slope boundaries 0,.12,.13,.14,.15,1;64 uniform bins cross these seams.',
                full24_original_C1_integral_range_transport_enclosed=False,
                actual_five_controls_installed=False,functional_terminal_identity_solved=False,current_whole_N_selected=False)
        return dict(report=report,targets=targets,amplitude=amplitude,corrections=self.corrections)


@native.inlet.source_precision
def run():
    began=time.monotonic();bridge,construction=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime() as runtime:
        owner=OriginalRcFullPredicateJointTargets(bridge);live=owner.targets()
    report=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
        actual_full_predicate_Rc_joint_targets=live['report'],original_bridge_construction=construction,
        original_source_runtime=runtime.record(),actual_original_Rc_joint_target_C0_Z_covers_installed=True,
        full24_original_C1_integral_range_transport_enclosed=False,actual_five_controls_installed=False,
        functional_terminal_identity_solved=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual accepted full-predicate N1024 Rc correction functions and fresh original source amplitude/ordinary-Z produce all five joint whole-Z finite targets. Separate P0/background and direct/spatial alternatives; no selected field points. Quantitative joint cancellation, original24-cell control record adapter, fixed point, terminal/global N/recursion/full NS remain open.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(encode(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Actual original full-predicate Rc five C0/Z joint targets generated',flush=True);return report


if __name__=='__main__':run()
