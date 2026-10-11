"""Selected actual n1 swirl enclosures including all residual feedback.

The signed zeroth source integral is retained as center. Its errors are
augmented by the uniform all-word remainder; no unknown feedback term is
set to zero. Axial/pressure accuracy and full solved-field contracts stay open.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_weighted_feedback_bound as feedback
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

op=feedback.op;weighted=feedback.weighted;core=feedback.core;ep=feedback.ep;domain=feedback.domain
HERE,PREFIX,sha=feedback.HERE,feedback.PREFIX,feedback.sha
NAME=PREFIX+'current_original_n1_coupled_swirl_enclosure.json.gz'
RECEIPT=PREFIX+'current_original_n1_coupled_swirl_enclosure_check.json'
GATE='current_original_n1_selected_actual_swirl_enclosures_include_full_residual_feedback'
OPEN=feedback.OPEN
FIELDS=('F1','d_x_F1')


def coupled_budget(value,error,k):
    result=core.integral_budget(value,error,k)
    for old,new in (('relative_initial_integral_width_certified','relative_coupled_swirl_width_certified'),
                    ('relative_initial_integral_width_target','relative_coupled_swirl_width_target'),
                    ('initial_source_integral_normalized_enclosure','coupled_swirl_normalized_enclosure'),
                    ('signed_model_integral_normalized_enclosure','signed_zeroth_source_integral_normalized_enclosure'),
                    ('initial_integral_exact_zero','coupled_swirl_exact_zero')):
        if old in result:result[new]=result.pop(old)
    result['remaining_B0_B1_feedback_included']=True
    return result


@dataclass(frozen=True,eq=False)
class OriginalN1CoupledSwirlEnclosurePacket:
    temporal_hierarchy_order:int=1


class CurrentOriginalN1CoupledSwirlEnclosure:
    @source_precision
    def __init__(self,accepted_feedback,require_checked=True):
        if type(accepted_feedback) is not feedback.CurrentOriginalN1WeightedFeedbackBound or not accepted_feedback.acceptance_loaded:
            raise ValueError('Accepted original all-word feedback bound required')
        accepted_feedback.assert_graph();self.feedback=accepted_feedback;self.source=accepted_feedback.axial.inner
        self.inputs=accepted_feedback.inputs;self.local=accepted_feedback.local;self.c=accepted_feedback.c
        self.bound_packet=accepted_feedback.domain(0);self.bound=accepted_feedback.report(self.bound_packet)
        self.family=copy.deepcopy(accepted_feedback.family);self.hashes=dict(accepted_feedback.hashes)
        for name in (feedback.NAME,feedback.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self._feedback=accepted_feedback;self._identities=(self.source,self.inputs,self.local,self.c)
        self._bound=op.snapshot(self.bound);self.acceptance_loaded=False;self._packets={}
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not receipt[GATE] or receipt['source_family']!=self.family or any(receipt[k] for k in OPEN):
                raise ValueError('Coupled swirl receipt differs')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Coupled swirl source changed '+name)
                domain.parent.parent.original.bind(self.hashes,name,digest)
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)

    def assert_graph(self):
        if (self.feedback is not self._feedback or (self.source,self.inputs,self.local,self.c)!=self._identities
            or self.source is not self.feedback.axial.inner or self.inputs is not self.feedback.inputs
            or self.local is not self.feedback.local or self.c is not self.feedback.c
            or not self.feedback.acceptance_loaded or self.hashes!=self._hashes or self.family!=self.feedback.family
            or op.snapshot(self.bound)!=self._bound):raise ValueError('Coupled swirl source graph changed')
        self.feedback.report(self.bound_packet);return self.feedback.assert_graph()

    @source_precision
    def evaluate(self,Z,coordinate):
        self.assert_graph();began=time.monotonic();c=self.c
        # The inherited original dispatch admits only exact endpoints in the
        # selected common interval. Integrating cells inside remain intervals.
        p=self.source.evaluate(Z,coordinate);sv=self.source.report(p);entry=self.source._packets[id(p)]
        values,initial_errors=entry[6:8];algebra=next(iter(values.values())).algebra
        tail=self.bound['all_word_feedback_tail']
        error=algebra.lift(op.Jet(c,[tail['value_remainder_upper'],tail['ordinary_Z1_remainder_upper']]))
        errors={field:initial_errors[field]+error for field in FIELDS}
        gradient=op.Jet(c,self.feedback.weighted.source.amplitude.gradient_jets(sv['Z']))
        physical={field:weighted.transform(v,gradient,(0,1,0,1),1) for field,v in values.items()}
        factor=weighted.absolute_jet(weighted.relative_amplitude(algebra,gradient,1,1))
        physical_errors={field:algebra.shift(factor*e,(0,1,0,1)) for field,e in errors.items()}
        physical_feedback=algebra.shift(factor*error,(0,1,0,1))
        budgets={field:{'Z'+str(k):coupled_budget(values[field],errors[field],k) for k in (0,1)} for field in FIELDS}
        physical_budgets={field:{'Z'+str(k):coupled_budget(physical[field],physical_errors[field],k) for k in (0,1)} for field in FIELDS}
        result=dict(source_family=self.family,Z=sv['Z'],query=sv['query'],actual_exit_geometry=sv['actual_exit_geometry'],
            normalized_signed_source_integral_center={field:op.export(v,1) for field,v in values.items()},
            normalized_initial_integral_error={field:op.export(v,1) for field,v in initial_errors.items()},
            normalized_all_feedback_error=op.export(error,1),normalized_total_coupled_error={field:op.export(v,1) for field,v in errors.items()},
            physical_signed_source_integral_center={field:op.export(v,1) for field,v in physical.items()},
            physical_all_feedback_error=op.export(physical_feedback,1),physical_total_coupled_error={field:op.export(v,1) for field,v in physical_errors.items()},
            normalized_coupled_swirl_budgets=budgets,physical_coupled_swirl_budgets=physical_budgets,
            passing_normalized_coupled_width_rows=sum(v['relative_coupled_swirl_width_certified'] for b in budgets.values() for v in b.values()),
            passing_physical_coupled_width_rows=sum(v['relative_coupled_swirl_width_certified'] for b in physical_budgets.values() for v in b.values()),
            budgets_now_include_all_residual_feedback=True,
            full_six_row_remaining_feedback_uniform_bound_included=True,feedback_tail_retained_indices=[0,0],
            ordinary_Z1_tail_tube_radius=tail['ordinary_Z1_tube_radius'],original_uniform_feedback_gate=self.bound[feedback.GATE],
            actual_Phi_collar_and_own_mean_pressure_inputs_retained=True,
            unknown_feedback_not_set_to_zero=True,physical_F0_Lambda_Z_derivative_gauge_added_once=True,
            only_F1_and_d_x_F1_components_enclosed_here=True,axial_pressure_mean_Q_accuracy_still_open=True,
            solved_coefficient_function_contracts_still_open=True,n1_own_repair_and_n2_still_open=True,
            execution_seconds=time.monotonic()-began,**{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        token=OriginalN1CoupledSwirlEnclosurePacket();state=(result,values,errors,physical,physical_errors,error)
        self._packets[id(token)]=(token,result,p,state,op.snapshot(state))
        return token

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalN1CoupledSwirlEnclosurePacket or entry is None or entry[0] is not packet:raise ValueError('Live coupled swirl packet required')
        self.source.report(entry[2])
        if op.snapshot(entry[3])!=entry[4]:raise ValueError('Coupled swirl packet changed')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_feedback):
    owner=CurrentOriginalN1CoupledSwirlEnclosure(accepted_feedback,require_checked=False)
    packets={name:owner.evaluate('.371',coordinate) for name,coordinate in (('actual_first_inlet','0'),('selected_actual_inner_exit','1'))}
    raw=dict(source_family=owner.family,samples={name:owner.report(p) for name,p in packets.items()},
             input_hashes=owner.hashes,**{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    print('N1_COUPLED_SWIRL_ENCLOSURES_COMPUTED',len(packets),flush=True)
    return owner,packets
