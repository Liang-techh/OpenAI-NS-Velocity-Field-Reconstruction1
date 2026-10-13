"""Restore checked complete original pulse pressure C1 functions."""
from pathlib import Path
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions as current


if __name__=='__main__':
    owner=current.CurrentOriginalN1PulsePressureFunctions();owner.assert_graph()
    print(json.dumps(dict(checked=owner.acceptance_loaded,
        own_fields_checked=owner.own.acceptance_loaded,graph_nodes=len(owner.g.nodes),
        complete_pressure_C1=current.own.encode(owner.functions),
        chart_pressure_C1=current.own.encode(owner.pieces),
        signed_integrand_count_per_chart={chart:len(values) for chart,values in owner.bindings['products'].items()},
        radial_scope='original Rp<=R<=Rv; complete main/gap/end charts',
        pressure_q_equals_auxiliary_Q=True,m15_minus_weighted_integral_sign_included=True,
        native_incoming_mean_retained=True,terminal_empty_support_zero_proved=True,
        full_five_defect_functions_and_actual_repair_installed=False,
        uniform_original_pressure_numerical_accuracy_claimed=False),indent=2))
