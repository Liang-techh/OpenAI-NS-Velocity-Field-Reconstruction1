"""Run and inspect the proved restoration/postrestore viscous contraction.

``--write-report`` constructs the new owner, evaluates only the three saved
axial queries, and runs the independent active-seam fixture.  Inspection reads
the saved packet and verifies hashes and serialized invariants without
constructing the owner or replaying an ancestor.
"""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path
import sys
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
DIRECTORY = ROOT / "experiments" / "root_st073"

from openai_ns_reconstruction.original_core_numeric import PREFIX  # noqa: E402
from openai_ns_reconstruction.original_active_viscous_m15 import (  # noqa: E402
    OriginalActiveViscousM15,
)


NAME = PREFIX + "compliant_current_original_active_viscous_m15_numeric.json"
PARENT_NAME = PREFIX + "compliant_current_original_micro_viscous_m15_numeric.json"
FIXTURE_NAME = "scripts/check_original_active_viscous_m15.py"
RUNNER_NAME = "scripts/run_original_active_viscous_m15.py"
SAVED_Z = ("-.3", "0", ".3")
MISSING_Z = "0.1"
RESTORATION = "restoration_finite_N_viscous_error"
POSTRESTORE = "postrestore_finite_N_viscous_error"
REPLACEMENT = "restoration_postrestore_finite_N_viscous_external_endpoints"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def check_source_seam_proof(proof):
    if not isinstance(proof, dict) or proof.get("passed") is not True:
        raise AssertionError(("source seam proof", proof))
    if (proof.get("left"), proof.get("right")) != ("restoration", "postrestore"):
        raise AssertionError(("active seam identity", proof))
    if proof.get("global_phase_not_reset") is not True:
        raise AssertionError("global phase reset at active seam")
    if proof.get("actual_modified_and_original_leading_traces_equal") is not True:
        raise AssertionError("modified/leading trace equality missing")
    if proof.get("independent_of_interval_caps") is not True:
        raise AssertionError("interval caps used as source identities")
    if proof.get("actual_signed_centers_evaluated") is not False:
        raise AssertionError("source cap was presented as a signed value")
    if proof.get("native_Jacobians") != [1, 1]:
        raise AssertionError(("native Jacobian identity", proof.get("native_Jacobians")))
    radius_flag_names = (
        "common_ordinary_Z_independent_radius",
        "radius_ordinary_Z_independent",
        "radius_Z_independent",
        "ordinary_Z_independent_radius",
    )
    if not any(proof.get(name) is True for name in radius_flag_names):
        raise AssertionError("explicit radius ordinary-Z independence proof missing")
    if proof.get("native_Jacobians_proved_equal") is not True or proof.get("physical_R_factor_in_both_ordinary_Z_rows") is not True:
        raise AssertionError("actual physical measure or R factor proof missing")

    rows = proof.get("identity_rows")
    if not isinstance(rows, list) or len(rows) != 26:
        raise AssertionError(("active seam identity row count", None if rows is None else len(rows)))
    expected = {
        ("geometry", "radius", 0, 83014, 84533),
        ("geometry", "global_offset", 0, 116796, 119414),
        ("geometry", "actual_phase", 0, 116798, 119416),
        ("geometry", "physical_measure", 0, 1, 1),
        ("modified", "U", 0, 177665, 183268),
        ("modified", "U", 1, 177666, 183269),
        ("modified", "U_y", 0, 177671, 183274),
        ("modified", "U_y", 1, 177672, 183275),
        ("modified", "m", 0, 177677, 183280),
        ("modified", "m", 1, 177678, 183281),
        ("modified", "m_y", 0, 177689, 183292),
        ("modified", "m_y", 1, 177690, 183293),
        ("modified", "Q", 0, 177760, 183363),
        ("modified", "Q", 1, 177764, 183367),
        ("modified", "Q_y", 0, 177848, 183451),
        ("modified", "Q_y", 1, 177852, 183455),
        ("modified", "R_Q_y", 0, 178007, 183610),
        ("modified", "R_Q_y", 1, 178008, 183611),
        ("original_leading", "V", 0, 95414, 84945),
        ("original_leading", "V", 1, 95415, 84946),
        ("original_leading", "V", 2, 95416, 0),
        ("original_leading", "m", 0, 95421, 96436),
        ("original_leading", "m", 1, 95422, 96437),
        ("original_leading", "m", 2, 95423, 96438),
        ("original_leading", "V_y", 0, 95606, 0),
        ("original_leading", "V_y", 1, 95607, 0),
    }
    seen = set()
    for row in rows:
        required = {"group", "quantity", "ordinary_Z_order", "left_node", "right_node", "canonical_endpoint_sha256"}
        if set(row) != required:
            raise AssertionError(("source seam row shape", row))
        digest = row["canonical_endpoint_sha256"]
        if not isinstance(digest, str) or len(digest) != 64:
            raise AssertionError(("source seam row digest", row))
        key = (row["group"], row["quantity"], row["ordinary_Z_order"], row["left_node"], row["right_node"])
        seen.add(key)
    if seen != expected:
        raise AssertionError(("actual source endpoint IDs", sorted(seen ^ expected)))
    if len(proof.get("empty_integrals_eliminated", [])) != 19:
        raise AssertionError("active seam empty-integral count changed")
    return dict(rows=len(rows), actual_source_endpoint_ids=True, empty_integrals=19, radius_ordinary_Z_independent=True)


def check_trace_chain(packet):
    trace = packet.get("signed_trace_chain")
    if not isinstance(trace, dict):
        raise AssertionError("signed trace chain missing")
    if trace.get("signed_coefficients") != {"restoration_native0": 1, "Rm": -1}:
        raise AssertionError(("signed trace coefficients", trace))
    if trace.get("cancelled_interfaces") != ["restoration_end"]:
        raise AssertionError(("cancelled interface", trace))
    return dict(
        signed_coefficients=copy.deepcopy(trace["signed_coefficients"]),
        cancelled_interfaces=copy.deepcopy(trace["cancelled_interfaces"]),
    )


def check_optional_external_inlet_proof(packet):
    """Keep the unresolved reference inlet diagnostic outside the proof aliases."""
    proof = packet.get("external_inlet_trace_probe")
    if not isinstance(proof, dict) or proof.get("passed") is not False:
        raise AssertionError(("unresolved external inlet diagnostic", proof))
    if (proof.get("left"), proof.get("right")) != ("reference", "restoration"):
        raise AssertionError(("external inlet seam identity", proof))
    required_true = (
        "global_phase_not_reset",
        "independent_of_interval_caps",
        "physical_R_factor_in_both_ordinary_Z_rows",
    )
    if any(proof.get(key) is not True for key in required_true):
        raise AssertionError(("external inlet proof flags", proof))
    if proof.get("actual_signed_centers_evaluated") is not False:
        raise AssertionError(("external inlet cap interpretation", proof))
    # The inlet proof may compare native charts with different Jacobian nodes;
    # its role is the physical trace identity.  Preserve the explicit result
    # when the producer supplies it without requiring equal native measures.
    if proof.get("native_Jacobians_proved_equal") is True:
        raise AssertionError(("external inlet native measure was overclaimed", proof))
    rows = proof.get("identity_rows")
    if not isinstance(rows, list) or not rows:
        raise AssertionError(("external inlet identity rows", proof))
    for row in rows:
        if not isinstance(row, dict) or not isinstance(row.get("canonical_endpoint_sha256"), str):
            raise AssertionError(("external inlet identity row", row))
    first = proof.get("unproved_rows", [None])[0]
    expected = ("modified", "U_y", 0, 171446, 177671)
    if not isinstance(first, dict) or tuple(first.get(k) for k in ("group", "quantity", "ordinary_Z_order", "left_node", "right_node")) != expected:
        raise AssertionError(("first unresolved physical inlet derivative", first))
    verified = {(r["group"], r["quantity"], r["ordinary_Z_order"]) for r in rows}
    if not {("geometry", k, 0) for k in ("radius", "global_offset", "actual_phase")} <= verified:
        raise AssertionError("exact inlet geometry diagnostic missing")
    if not {("modified", "U", k) for k in (0, 1)} <= verified:
        raise AssertionError("exact inlet velocity value diagnostic missing")
    return dict(present=True, identity_rows=len(rows), unproved_rows=len(proof["unproved_rows"]),
        native_measures_not_identified=True,source_inlet_identity_not_claimed=True)


def check_contract(packet):
    required = {
        "raw_213_verbatim": True,
        "only_proved_viscous_seam_contracted": True,
        "modified_and_leading_traces_both_verified": True,
        "both_external_traces_retained": True,
        "original_both_full_positive_aggregates_retained": True,
        "reference_inlet_reinterpretation_not_used": True,
        "ordinary_Z1_and_physical_R_once": True,
        "all_other_source_volume_width_errors_retained": True,
        "independent_target_unchanged": True,
        "ancestors_not_reexecuted": True,
    }
    contract = packet.get("contract")
    if not isinstance(contract, dict) or any(contract.get(k) is not v for k, v in required.items()):
        raise AssertionError(("contract flags", contract))
    limits = packet.get("limits")
    if not isinstance(limits, dict) or limits.get("restoration_postrestore_viscous_internal_trace_resolved") is not True:
        raise AssertionError(("resolved active seam limit", limits))
    for key, value in limits.items():
        if key != "restoration_postrestore_viscous_internal_trace_resolved" and value is not False:
            raise AssertionError(("saved open limit overclaimed", key, value))
    return dict(contract_flags=copy.deepcopy(required), internal_trace_resolved=True, all_other_limits_open=True)


def check_external_budget(packet, raw):
    witness = packet.get("external_endpoint_budget_witness")
    if not isinstance(witness, dict):
        raise AssertionError("external endpoint budget witness missing")
    source = witness.get("source")
    if not isinstance(source, dict) or source.get("passed") is not True:
        raise AssertionError(("external budget source", source))
    if source.get("source") != "src/openai_ns_reconstruction/original_restore_m15.py":
        raise AssertionError(("external budget source path", source))
    if source.get("method") != "OriginalRestoreM15.finite_errors":
        raise AssertionError(("external budget source method", source))
    if source.get("physical_radius_power_for_viscous") != 1:
        raise AssertionError(("physical R power", source))
    if source.get("positive_upper_lower_union_not_difference") is not True:
        raise AssertionError("external endpoint budgets were subtracted")
    if source.get("source_cap_not_function_value") is not True:
        raise AssertionError("external budget cap treated as a function value")
    aggregates=[]
    for name,index in (("restoration",139),("postrestore",155)):
        retained = witness.get(name+"_positive_aggregate_retained")
        if retained != raw[index]["log_scale_terms"][0]:
            raise AssertionError(("original full aggregate changed", name))
        if retained.get("positive_error_bound_atom_only") is not True or retained.get("source_value_selected") is not False:
            raise AssertionError(("aggregate interpretation", name))
        if any(r["source_enclosure_exact_mpf_tuples"] != [[1,1,0,1],[0,1,0,1]] for r in raw[index]["coefficient_C0_C1"]):
            raise AssertionError(("original symmetric coefficient budget", name))
        aggregates.append(retained)
    for key in ("reference_lower_packet_not_reinterpreted_as_restoration_inlet","no_new_separate_external_endpoint_bound"):
        if witness.get(key) is not True:raise AssertionError(("unproved inlet budget reused",key))
    union = witness.get("positive_union")
    if not isinstance(union, dict) or union.get("positive_error_bound_atom_only") is not True or union.get("source_value_selected") is not False:
        raise AssertionError(("positive union interpretation", union))
    tuples = union.get("source_enclosure_exact_mpf_tuples")
    if not isinstance(tuples, list) or len(tuples) != 2 or any(len(value) != 4 for value in tuples):
        raise AssertionError(("positive union exact tuple", union))
    ctx=mp.mp.clone();ctx.dps=160;iv=MPIntervalContext();iv.dps=160
    old_upper=max(ctx.make_mpf(tuple(r["source_enclosure_exact_mpf_tuples"][1])) for r in aggregates)
    required=iv.mpf(old_upper)+iv.log(2)
    if ctx.make_mpf(tuple(tuples[1])) < ctx.make_mpf(required._mpi_[1]):
        raise AssertionError("serialized union does not retain both complete budgets")
    return dict(source_ast_passed=True, full_original_aggregates=2,original_budget_sum_retained=True,
        no_new_separate_inlet_bound=True, inherited_budget_source_identity_still_open=True,
        physical_R_once=True,cap_as_value=False)


def check_packet(packet, parent_packet, Z):
    if packet.get("Z") != Z or packet.get("source_family") != parent_packet.get("source_family"):
        raise AssertionError(("packet identity", Z))
    raw = packet.get("raw_m15_C0_C1_terms")
    parent_raw = parent_packet.get("raw_m15_C0_C1_terms")
    if not isinstance(raw, list) or len(raw) != 213 or raw != parent_raw:
        raise AssertionError(("raw 213 terms differ from parent", Z))
    if packet.get("original_raw_terms") != 213 or packet.get("previous_reduced_terms") != 200:
        raise AssertionError(("raw/previous reduced counts", Z))
    reduced = packet.get("correlated_reduced_m15_C0_C1_terms")
    old_reduced = parent_packet.get("correlated_reduced_m15_C0_C1_terms")
    if not isinstance(reduced, list) or len(reduced) != 199 or not isinstance(old_reduced, list) or len(old_reduced) != 200:
        raise AssertionError(("active reduced count", Z))
    if packet.get("correlated_reduced_terms") != 199:
        raise AssertionError(("serialized active reduced count", Z))
    old_names = [term["name"] for term in old_reduced]
    new_names = [term["name"] for term in reduced]
    if old_names.count(RESTORATION) != 1 or old_names.count(POSTRESTORE) != 1:
        raise AssertionError(("external viscous sectors", Z))
    if RESTORATION in new_names or POSTRESTORE in new_names or new_names.count(REPLACEMENT) != 1:
        raise AssertionError(("replacement endpoint name", Z))
    expected_names = [
        REPLACEMENT if term["name"] == POSTRESTORE else term["name"]
        for term in old_reduced
        if term["name"] != RESTORATION
    ]
    if new_names != expected_names:
        raise AssertionError(("active reduced term order", Z))
    unchanged_old = [term for term in old_reduced if term["name"] not in {RESTORATION, POSTRESTORE}]
    unchanged_new = [term for term in reduced if term["name"] != REPLACEMENT]
    if len(unchanged_old) != 198 or len(unchanged_new) != 198:
        raise AssertionError(("unchanged active terms", Z))
    if any(canonical(left) != canonical(right) for left, right in zip(unchanged_old, unchanged_new)):
        raise AssertionError(("non-viscous active terms changed", Z))
    replacement = packet.get("replacement_outer_endpoints")
    replacements = [term for term in reduced if term["name"] == REPLACEMENT]
    if len(replacements) != 1 or replacements[0] != replacement:
        raise AssertionError(("replacement serialization", Z))
    binding = replacement.get("actual_function_binding", {})
    required_binding = {
        "error_enclosure_only": True,
        "actual_internal_trace_proved_equal": True,
        "both_external_traces_retained_without_zero_assumption": True,
        "ordinary_Z1": True,
        "no_error_cap_subtraction": True,
    }
    if any(binding.get(key) is not value for key, value in required_binding.items()):
        raise AssertionError(("replacement binding", Z, binding))
    if binding.get("function") != "delta(RQy)(restoration native0)-delta(RQy)(postrestore native1)":
        raise AssertionError(("replacement physical endpoint", Z))
    old_budget = packet.get("old_budget_terms_retained_as_provenance")
    if old_budget != [raw[139], raw[155]]:
        raise AssertionError(("external budget provenance", Z))
    if packet.get("independent_target") != parent_packet.get("independent_target"):
        raise AssertionError(("independent target changed", Z))
    seam = packet.get("source_seam_proofs")
    if not isinstance(seam, list) or len(seam) != 1:
        raise AssertionError(("source seam proof list", Z))
    seam_check = check_source_seam_proof(seam[0])
    chain_check = check_trace_chain(packet)
    inlet_check = check_optional_external_inlet_proof(packet)
    contract_check = check_contract(packet)
    budget_check = check_external_budget(packet, raw)
    return dict(
        Z=Z,
        raw_terms=213,
        raw_same_parent=True,
        previous_reduced_terms=200,
        correlated_reduced_terms=199,
        unchanged_reduced_terms=198,
        external_budget_indices=[139, 155],
        both_external_budgets_retained=True,
        replacement_name=REPLACEMENT,
        source_seam=seam_check,
        external_inlet_source=inlet_check,
        signed_trace_chain=chain_check,
        external_budget=budget_check,
        contract=contract_check,
        independent_target_unchanged=True,
    )


def run():
    started = time.perf_counter()
    from check_original_active_viscous_m15 import run as run_fixture

    fixture = run_fixture()
    if not fixture.get("passed"):
        raise AssertionError("independent active viscous fixture failed")
    owner = OriginalActiveViscousM15(precision=100)
    parent_by_Z = {row["Z"]: row for row in owner.parent["samples"]}
    packets = [owner.evaluate(Z) for Z in SAVED_Z]
    packet_checks = [check_packet(packet, parent_by_Z[Z], Z) for packet, Z in zip(packets, SAVED_Z)]
    missing_rejected = 0
    try:
        owner.evaluate(MISSING_Z)
    except (ValueError, ArithmeticError):
        missing_rejected = 1
    else:
        raise AssertionError(("unsupported saved Z admitted", MISSING_Z))
    hashes = dict(owner.hashes)
    code = dict(owner.code_hashes)
    code[FIXTURE_NAME] = sha(ROOT / FIXTURE_NAME)
    code[RUNNER_NAME] = sha(ROOT / RUNNER_NAME)
    return dict(
        source_family=owner.family,
        samples=packets,
        input_hashes=hashes,
        implementation_hashes=code,
        checks=dict(
            passed=True,
            independent_fixture=fixture,
            saved_Z=SAVED_Z,
            active_viscous_packets=len(packets),
            raw_213_same_parent=True,
            previous_correlated_200=True,
            correlated_reduced_199=True,
            unchanged_reduced_198=True,
            replaced_external_indices=[139, 155],
            both_external_endpoint_budgets_retained=True,
            physical_R_once=True,
            ordinary_Z1=True,
            source_identity_rows_per_packet=26,
            source_endpoint_ids_actual=True,
            radius_ordinary_Z_independent=True,
            cap_as_function_value=False,
            independent_target_unchanged=True,
            missing_Z_rejected=missing_rejected,
            production_is_new_calculation_only=True,
        ),
        packet_checks=packet_checks,
        execution_seconds=round(time.perf_counter() - started, 3),
        limits=dict(
            restoration_postrestore_viscous_internal_trace_resolved=True,
            remaining_external_trace_signed_values_resolved=False,
            other_internal_traces_resolved=False,
            independent_numeric_Rp_seam=False,
            global_five_moment_numeric_join=False,
            exact_heat_exterior=False,
            stress_cone=False,
            flat_remainder=False,
            genuine_temporal_n_ge_2=False,
            full_corrected_NS=False,
        ),
    )


def inspect(path):
    """Verify hashes and saved invariants without constructing the owner."""
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("checks", {}).get("passed") is not True or data.get("source_family") is None:
        raise AssertionError("saved active viscous report is not passed")
    for base, group in ((DIRECTORY, data["input_hashes"]), (ROOT, data["implementation_hashes"])):
        for name, digest in group.items():
            if sha(base / name) != digest:
                raise AssertionError(("hash changed", name))
    parent = json.loads((DIRECTORY / PARENT_NAME).read_text(encoding="utf-8"))
    parent_by_Z = {row["Z"]: row for row in parent["samples"]}
    samples = data.get("samples")
    if not isinstance(samples, list) or tuple(sample["Z"] for sample in samples) != SAVED_Z:
        raise AssertionError("saved active Z packet set changed")
    for sample in samples:
        check_packet(sample, parent_by_Z[sample["Z"]], sample["Z"])
    fixture = data["checks"].get("independent_fixture")
    if not isinstance(fixture, dict) or fixture.get("passed") is not True:
        raise AssertionError("independent active fixture evidence missing")
    if fixture.get("physical", {}).get("checks", {}).get("both_external_endpoints_nonzero") is not True:
        raise AssertionError("external endpoint fixture evidence missing")
    active_fingerprint = fixture.get("active_fingerprint")
    required_active_fingerprint = (
        "exact_additive_normalization",
        "nonrational_negative_cancellation",
        "opaque_source_argument_distinction",
        "radius_endpoint_mismatch",
        "phase_endpoint_mismatch",
    )
    if not isinstance(active_fingerprint, dict) or any(
        active_fingerprint.get(key) is not True for key in required_active_fingerprint
    ):
        raise AssertionError("active fingerprint adversarial evidence missing")
    if data["checks"].get("missing_Z_rejected") != 1:
        raise AssertionError("unsupported-Z rejection evidence changed")
    if data["checks"].get("source_identity_rows_per_packet") != 26:
        raise AssertionError("active source row count changed")
    if data["checks"].get("unchanged_reduced_198") is not True:
        raise AssertionError("unchanged reduced-term evidence missing")
    if data["checks"].get("both_external_endpoint_budgets_retained") is not True:
        raise AssertionError("external endpoint retention evidence missing")
    if data["checks"].get("cap_as_function_value") is not False:
        raise AssertionError("cap interpretation overclaimed")
    for key, value in data["limits"].items():
        expected = key == "restoration_postrestore_viscous_internal_trace_resolved"
        if value is not expected:
            raise AssertionError(("saved limit changed", key, value))
    return data


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-report", action="store_true")
    parser.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = DIRECTORY / NAME
    data = inspect(path) if args.inspect_report else run()
    if args.write_report:
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8", newline="\n")
        temporary.replace(path)
    print(json.dumps(dict(checks=data["checks"], execution_seconds=data["execution_seconds"]), indent=2), flush=True)
