"""Check independent source germs and reject asymmetric patch changes."""
import json
from pathlib import Path
import sympy as s
import lei_ren_part1_paper_compliant_current_patch_support_velocity_pressure as source


def run():
    data = json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name) != digest:
            raise ValueError('Accepted patch support source changed: '+name)
    theorem = source.exact_patch_support_theorem()
    if theorem != data['exact_actual_patch_support_velocity_pressure_theorem']:
        raise ValueError('Actual two-sided source theorem differs')
    edges = theorem['exact_support_source_theorems']
    if tuple(map(int,edges)) != source.EDGES:
        raise ValueError('Actual six support inventory differs')
    if not edges['71']['last_edge_uses_current_five_functional_implicit_equations']:
        raise ValueError('Last support needs the actual unique implicit closure')
    count = sum(len(value['identities']) for value in edges.values())
    negatives = {}
    for perturb,label in (('H','Utheta_over_Pstar'),('mass','Ur_over_sqrt_Rm_over_2'),('P0','P_over_Pstar2')):
        _,packets,_ = source.replayed_trace_packets(source.SourceAST(),perturb)
        left = packets[0]['physical_velocity_pressure_y_Z_mixed4'][label]['y0_Z0']
        right = packets[1]['physical_velocity_pressure_y_Z_mixed4'][label]['y0_Z0']
        if s.cancel(left-right) == 0:
            raise ArithmeticError('A changed independent right germ was hidden: '+perturb)
        negatives[perturb] = dict(right_source_changed=True,source_row_mismatch_detected=True)
    if any(data.get(flag) is not False for flag in (
        'mapped_physical_spatial4_time1_trace_views_available',
        'current_modified_global_velocity_interfaces_certified',
        'common_N_modified_cones_energy_recursion_full_NS_certified')):
        raise ValueError('Local source theorem exceeds completed scope')
    hashes = dict(data['input_hashes'])
    for name in (source.NAME,Path(__file__).name):
        hashes[name] = source.sha(name)
    receipt = dict(all_passed=True,source_family=theorem['source_family'],
        exact_support_source_identities=count,
        generic_mixed4_and_unit_source_identities=len(theorem['generic_two_germ_mixed4_and_current_unit_identities']),
        six_actual_source_edges_certified=list(source.EDGES),
        independently_changed_right_germ_negative_controls=negatives,
        support_source_identity_not_numeric_zero_containment_or_self_comparison=True,
        same_unique_current_implicit_family_and_original_P0_consumed=True,
        current_six_patch_support_velocity_pressure_mixed4_source_joins_certified=True,
        mapped_physical_spatial4_time1_trace_views_available=False,
        current_modified_global_velocity_interfaces_certified=False,
        common_N_modified_cones_energy_recursion_full_NS_certified=False,input_hashes=hashes)
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Six actual patch support source joins PASS:',count,'support identities;',
        receipt['generic_mixed4_and_unit_source_identities'],'mixed/source/unit identities; 3 asymmetric changes detected',flush=True)
    return receipt


if __name__=='__main__':
    run()
