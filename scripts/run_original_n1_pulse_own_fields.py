"""Restore checked original C3 coefficients and own main/gap/end fields."""
from pathlib import Path
import json
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_pulse_own_fields as current


if __name__=='__main__':
    owner=current.CurrentOriginalN1PulseOwnFields();owner.assert_graph()
    charts={key:dict(query=field['domain']['coordinate'],
        physical_Utheta_U_Uy_m_my_C3={name:current.encode(field[name])
            for name in ('Utheta','U','U_y','m','m_y')},
        auxiliary_Q_Qy_and_physical_Ur_C2={name:current.encode(field[name])
            for name in ('Q','Q_y','Ur')},
        cumulative_swirl_base_pressure_C3=current.encode(field['absolute_swirl_base_pressure']))
        for key,field in owner.fields.items()}
    print(json.dumps(dict(checked=owner.acceptance_loaded,
        actual_C3_coefficients_checked=owner.jets.acceptance_loaded,
        graph_nodes=len(owner.g.nodes),charts=charts,
        physical_Pstar_applied_once=True,own_gap_history_retained=True,
        Q_definition='Ur/sqrt(R/2)',
        full_pulse_pressure_and_five_moment_repair_closed=False,
        uniform_original_numerical_accuracy_claimed=False),indent=2))
