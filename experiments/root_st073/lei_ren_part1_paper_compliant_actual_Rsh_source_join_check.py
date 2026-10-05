"""Focused acceptance of the new current-source Rsh mixed4 interface."""
import json
from pathlib import Path
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import (
    build, sha, HERE, PREFIX, SCOPES)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes

NAME=PREFIX+"actual_Rsh_source_join.json"


def run():
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    expected=build()
    if raw!=expected:
        raise ValueError("Current exact-source Rsh admission differs from fresh source proof")
    proof=raw["boundary_source_proof"]
    grid={"y"+str(k)+"_Z"+str(n) for k in range(5) for n in range(5-k)}
    rows=proof["physical_mixed4_rows_implied_by_exact_function_identities"]
    if len(rows)!=9 or any(set(keys)!=grid for keys in rows.values()):
        raise ValueError("All nine physical functions and all mixed4 rows required")
    if (not proof["passed"] or proof["exact_six_centered_history_function_identities"]!=6
            or proof["exact_original_flat_log_mixed4_boundary_rows"]!=15
            or proof["total_physical_mixed4_rows_implied"]!=135
            or proof["finite_left_neighborhood_replaced"]):
        raise ValueError("Current Rsh boundary theorem scope changed")
    if any(raw[k] for k in SCOPES) or any(raw["native_historical_Rsh_flags_preserved"].values()):
        raise ValueError("Source join incorrectly promotes point/global/downstream scope")
    result=dict(actual_five_defect_family_sha256=raw["actual_five_defect_family_sha256"],
        implicit_source_sha256=raw["implicit_source_sha256"],
        datum_enclosure_sha256=raw["datum_enclosure_sha256"],
        current_input_hashes_checked=len(raw["input_hashes"]),
        current_source_AST_bindings=raw["current_Rsh_source_bindings"],
        exact_six_centered_history_function_identities=6,
        exact_flat_log_boundary_mixed4_rows=15,
        physical_functions_with_exact_boundary_source_identity=9,
        physical_mixed4_rows_implied_by_exact_source_identities=135,
        independently_reused_physical_fixture_receipts=raw["generic_fixture_receipts"],
        current_Rsh_parent_P0_V_E_and_physical_normalization_verified=True,
        source_equality_proof_not_enclosure_overlap=True,
        original_flat_cutoff_endpoint_derivatives_checked=5,
        constant_power_boundary_extension_does_not_replace_left_field=True,
        current_Rsh_source_functional_join_certified=True,
        **dict.fromkeys(SCOPES,False),all_passed=True,
        input_hashes={**raw["input_hashes"],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix(".json").write_bytes((json.dumps(result,indent=2)+"\n").encode("utf8"))
    print("PASS current Rsh functional mixed4 join: 6 exact histories and 135 implied physical rows",flush=True)
    return result


if __name__=="__main__":
    run()
