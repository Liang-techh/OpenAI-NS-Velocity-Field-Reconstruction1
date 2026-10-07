"""Carry actual first-bridge C1 correction functions through bridge_second.

The unchanged original loop is bounded with the accepted whole-period C1
backend. The true microscopic width is applied once. Correction memory,
original endpoint backgrounds and separate P0/P0_Z stay distinct.
"""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_first_bridge_C1_histories as prior_stage

transfer=prior_stage.transfer;history=prior_stage.history;density=prior_stage.density
prior=prior_stage.prior;native=prior_stage.native;packets=prior_stage.packets
HERE,PREFIX,sha=prior_stage.HERE,prior_stage.PREFIX,prior_stage.sha
RATES=prior_stage.RATES;ep=prior_stage.ep;ZERO=prior_stage.ZERO;DZ=prior_stage.DZ
NAME=PREFIX+'current_native_second_bridge_C1_histories.json'
RECEIPT=PREFIX+'current_native_second_bridge_C1_histories_check.json'
GATE='current_original_inlet_to_phase2_second_bridge_C1_history_covers_executed'


class NativeSecondBridgeC1Histories:
    def __init__(self,owner):
        if type(owner) is not prior_stage.NativeFirstBridgeC1Histories:raise ValueError('Accepted original first-bridge C1 owner required')
        receipt=json.loads((HERE/prior_stage.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[prior_stage.GATE] or receipt['source_family']!=owner.family:
            raise ValueError('Same accepted original whole first-bridge C1 stage required')
        self.owner=owner;self.transfer=owner.owner;self.ctx=owner.ctx;self.family=owner.family
        self.coordinates=owner.coordinates;self.service=owner.service
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({name:sha(name) for name in (prior_stage.RECEIPT,Path(__file__).name,
            PREFIX+'current_inner_relaxed_inputs.py',PREFIX+'current_inner_relaxed_inputs_check.json')})
    @native.inlet.source_precision
    def second_bridge(self,Z=(-1,1),N=1024):
        N=density.spatial.density.candidate_integer(N)
        if N<160:raise ValueError('Accepted whole-period cover requires candidate N>=160')
        incoming=self.owner.first_bridge(Z,N=N);c=self.ctx
        geometry=self.transfer.geometry.cell('bridge_second',1,2)
        source=self.owner.q_owner.query('bridge_second',Z,geometry['coordinate'])
        roots=source['source']['roots'];signed=self.transfer.owner.signed_owner
        root_owner=self.owner.q_owner.owner.owner
        positive=root_owner.decode(root_owner.inventory['bridge_second']['actual_positive_denominator_theorem'])
        if not positive['source_function_positivity_not_inferred_from_saved_box']:
            raise ValueError('Original whole-domain analytic a lower required')
        eta=packets.interval(c,root_owner.scales['selected_positive_eta_log'])
        dstar=packets.interval(c,root_owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
        primitives=prior_stage.whole_period_C1(roots,eta,positive['log_actual_a_positive_lower'],dstar)
        E=roots['E'][ZERO];E_Z=roots['E'][DZ];packet=source['source']['packet']
        def axial(k):
            row=prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
            return signed.leaf(prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),E.scale.bases,E.ledger)
        got=density.density_Z_kernels(E,E_Z,axial(0),axial(1),primitives['values'],N)
        factors={key:transfer.true_width_kernel(self.coordinates,geometry,rate) for key,rate in RATES.items()}
        values={key:self.coordinates.rebase(got['kernels'][key],self.family)*factors[key]['mass'] for key in RATES}
        jets={key:self.coordinates.rebase(got['Z_derivatives'][key],self.family)*factors[key]['mass'] for key in RATES}
        operator=history.C1DuhamelOperator(self.coordinates)
        transfer.append_true_cell(operator,geometry,values,jets,self.family)
        correction=operator.apply(incoming['correction']['values'],incoming['correction']['Z_derivatives'],self.family)
        background=history.packet_history_functions(self.transfer.owner.native.query('bridge_second',Z,2),self.coordinates,signed)
        own={key:background['originals'][key]+correction['values'][key] for key in RATES}
        own_Z={key:background['Z_derivatives'][key]+correction['Z_derivatives'][key] for key in RATES}
        record=dict(source_family=self.family,Z_box=c.mpf(Z),candidate_N=N,
            known_original_first_bridge_phase1_history=incoming['record'],
            actual_whole_second_bridge_geometry=geometry['record'],
            original_whole_domain_a_positive_certificate=positive,
            original_whole_second_bridge_source_C1=source['record'],
            original_whole_period_C1_cover=primitives['record'],
            true_log_radius_signed_C0_contributions={k:v.record() for k,v in values.items()},
            true_log_radius_signed_Z_contributions={k:v.record() for k,v in jets.items()},
            actual_second_bridge_C1_operator=operator.record(),
            actual_inherited_phase1_correction_C0={k:v.record() for k,v in incoming['correction']['values'].items()},
            actual_inherited_phase1_correction_Z={k:v.record() for k,v in incoming['correction']['Z_derivatives'].items()},
            actual_phase2_correction_C0={k:v.record() for k,v in correction['values'].items()},
            actual_phase2_correction_Z={k:v.record() for k,v in correction['Z_derivatives'].items()},
            original_phase2_background_and_separate_P0_Z=background['record'],
            actual_phase2_own_history_C0={k:v.record() for k,v in own.items()},
            actual_phase2_own_history_Z={k:v.record() for k,v in own_Z.items()},
            original_exact_radius_seam=geometry['record']['exact_original_radius_Jacobian_identities'],
            correction_incoming_is_actual_phase1_not_own_history_or_zero_reset=True,
            no_gap_between_original_inlet_and_phase2=True,
            analytic_a_lower_uniform_on_original_bridge_second_and_whole_Z=True,
            source_C0_Z_queried_on_full_native_coordinate_box=True,
            original_background_and_P0_not_reset_or_double_added=True,
            microscopic_width_applied_once_in_original_log_radius_integral=True,
            first_and_second_bridge_covers_are_conservative_not_tight_terminal_error=True,
            global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,incoming=incoming,geometry=geometry,source=source,primitives=primitives,
            factors=factors,operator=operator,contributions=values,Z_derivatives=jets,
            correction=correction,background=background,own=own,own_Z=own_Z)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        c1=density.NativeDensityC1LocalIntegrals(prior_stage.first.NativePhaseFirstJets(prior_stage.slow.NativeQSlowJets(
            prior_stage.current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge))))))
        owner=NativeSecondBridgeC1Histories(prior_stage.NativeFirstBridgeC1Histories(
            transfer.NativeTrueChartC1Transfer(history.NativeC1HistoryTransfer(c1))))
        records={}
        for name,Z in (('whole_Z',(-1,1)),('Z_interval',('.49','.51'))):
            records[name]=owner.second_bridge(Z,N=1024)['record'];print('Actual whole second-bridge C1 exit:',name,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},candidate_N=1024,native_Z_query_count=2,
        actual_whole_second_bridge_C1_history_records=records,actual_original_inlet_to_phase2_C1_histories_installed=True,
        actual_phase1_memory_not_reset_or_double_added=True,whole_Z_second_bridge_positive_a_lower_is_analytic=True,
        covers_can_be_very_wide_tight_signed_oscillatory_error_still_required=True,
        global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Actual original inlet sc/2 to bridge_second phase2 correction and own C1 history covers on whole Z[-1,1] and[.49,.51]. Same original incoming phase1 correction, true h_bridge width and separate original endpoint background/P0_Z. Conservative bounds; downstream macro route, terminal repair/common N/cone/recursion/full NS remain open.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
