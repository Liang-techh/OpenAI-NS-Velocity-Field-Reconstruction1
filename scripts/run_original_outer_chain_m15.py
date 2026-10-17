"""Write or inspect the source-bound six-seam outer-chain m15 packet.

``--write-report`` constructs the new owner once and runs the independent
local-germ fixture once.  ``--inspect-report`` validates only the serialized
packet and its pinned files; it never constructs an owner or imports the
fixture.  The two modes are mutually exclusive.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
from typing import Any, Mapping

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
DIRECTORY = ROOT / "experiments" / "root_st073"
for path in (str(ROOT), str(SRC), str(SCRIPTS)):
    if path not in sys.path:
        sys.path.insert(0, path)

from openai_ns_reconstruction.original_core_numeric import PREFIX  # noqa: E402
from openai_ns_reconstruction.original_outer_chain_numeric_m15 import (  # noqa: E402
    NEW_EXTERNAL_TERM,
    NEW_REDUCED_COUNT,
    OLD_EXTERNAL_TERM,
    OUTER_M15_REPORT,
    PACKET_SPECS,
    PREVIOUS_REDUCED_COUNT,
    RAW_COUNT,
    REMOVED_OUTER_TERMS,
    REPORT_SHA256,
    SAVED_Z,
    UNCHANGED_COUNT,
)


NAME = PREFIX + "compliant_current_original_outer_chain_m15_numeric.json"
FIXTURE_NAME = "scripts/check_original_outer_chain_m15.py"
RUNNER_NAME = "scripts/run_original_outer_chain_m15.py"
PARENT = PREFIX + "compliant_current_original_Rm_inlet_m15_numeric.json"
PARENT_SHA = "3f936cbb5f91fd05677f0c292ed783ad72788b0971ed869298f43be9dfc0850d"

SEAMS = (
    "patch_reference",
    "reference_slope",
    "slope_axial",
    "axial_buffer",
    "buffer_transition",
    "transition_power",
)
JOIN_EXPECTED = {
    "patch_reference": ("Rm_patch", "Rh_reference", 2948, 3348),
    "reference_slope": ("Rh_reference", "O2_slope", 0, 0),
    "slope_axial": ("O2_slope", "O2_axial", 1, 0),
    "axial_buffer": ("O2_axial", "O2_buffer", 1, 0),
    "buffer_transition": ("O2_buffer", "O3_transition", 4, 0),
    "transition_power": ("O3_transition", "O3_power", 1, 0),
}


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def digest(value: Any) -> str:
    return hashlib.sha256(canonical(value).encode("utf-8")).hexdigest()


def require(condition: bool, message: Any) -> None:
    if not condition:
        raise AssertionError(message)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _resolve(name: str, base: Path | None = None) -> Path:
    """Resolve a saved hash key without allowing a missing file to pass."""
    raw = Path(name)
    candidates = []
    if raw.is_absolute():
        candidates.append(raw)
    else:
        candidates.extend((ROOT / raw, (base or DIRECTORY) / raw,
                           ROOT / "experiments" / "root_st073" / raw))
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0] if candidates else raw


def check_hash_group(group: Any, label: str, base: Path | None = None) -> None:
    require(isinstance(group, Mapping) and group, (label, "hash map missing"))
    for name, expected in group.items():
        require(isinstance(name, str) and isinstance(expected, str) and len(expected) == 64,
                (label, name, "malformed digest"))
        path = _resolve(name, base)
        require(path.is_file(), (label, name, "pinned file missing"))
        require(sha(path) == expected, (label, name, "pinned file changed"))


def _ctx_pair():
    ctx = mp.mp.clone()
    ctx.dps = 150
    iv = MPIntervalContext()
    iv.dps = 180
    return ctx, iv


def _record_tuples(row: Any):
    require(isinstance(row, Mapping), ("record", row))
    values = row.get("source_enclosure_exact_mpf_tuples")
    if values is None and "lower_exact_mpf_tuple" in row:
        values = [row["lower_exact_mpf_tuple"], row["upper_exact_mpf_tuple"]]
    require(isinstance(values, list) and len(values) == 2 and
            all(isinstance(value, list) and len(value) == 4 for value in values),
            ("record tuples", row))
    return values


def _decode_record(row: Mapping[str, Any], ctx, iv):
    return iv.mpf([ctx.make_mpf(tuple(value)) for value in _record_tuples(row)])


def _bounds(ctx, value):
    return tuple(ctx.make_mpf(row) for row in value._mpi_)


def _check_record(row: Any, label: str) -> None:
    require(isinstance(row, Mapping), (label, "missing record"))
    _record_tuples(row)
    require(row.get("positive_error_bound_atom_only") is True and
            row.get("source_value_selected") is False,
            (label, "record is not an error-only positive atom"))


def _check_positive_union(values: Any, record: Any, label: str) -> dict[str, Any]:
    require(isinstance(values, list), (label, "union inputs missing"))
    _check_record(record, label)
    ctx, iv = _ctx_pair()
    decoded = [_decode_record(value, ctx, iv) for value in values]
    require(decoded, (label, "empty positive union"))
    top = max(_bounds(ctx, value)[1] for value in decoded)
    expected = iv.mpf(top) + iv.log(len(decoded))
    require(_record_tuples(record) == [list(value) for value in expected._mpi_],
            (label, "serialized MPF union changed"))
    return {"components": len(decoded), "exact_mpf_union": True, "cap_as_value": False}


def _find(rows: Any, **criteria: Any) -> Mapping[str, Any]:
    require(isinstance(rows, list), ("rows", criteria))
    for row in rows:
        if isinstance(row, Mapping) and all(row.get(key) == value for key, value in criteria.items()):
            return row
    raise AssertionError(("missing row", criteria))


def _check_fixture_shape(fixture: Any) -> dict[str, Any]:
    require(isinstance(fixture, Mapping) and fixture.get("passed") is True,
            "independent outer-chain fixture failed")
    expected = {
        "synthetic_local_germs_only": True,
        "original_source_values_evaluated": False,
        "original_Poisson_inverse_evaluated": False,
        "NS_residual": False,
        "finite_N": 17,
        "interface_count": 6,
        "power_transport_with_nonzero_mean": True,
    }
    for key, value in expected.items():
        require(fixture.get(key) == value, ("fixture flag", key, fixture.get(key)))
    controls = fixture.get("controls")
    expected_controls = {
        "extra_N_on_fast_term_rejected",
        "frozen_J_second_derivative_rejected",
        "phase_reset_rejected",
        "nonzero_incoming_mean",
        "power_mean_not_zero",
        "power_D_not_zero",
        "missing_physical_R_rejected",
        "wrong_telescoping_sign_rejected",
        "inlet_retained",
        "outlet_retained",
        "leading_subtraction_retained",
        "genuine_ordinary_Z1",
    }
    require(isinstance(controls, Mapping) and set(controls) == expected_controls and
            all(value is True for value in controls.values()), ("fixture controls", controls))
    require(fixture.get("control_count") == len(expected_controls), "fixture control count changed")
    require(mp.mpf(str(fixture.get("maximum_supporting_difference"))) < mp.mpf("1e-38"),
            "fixture tolerance changed")
    chain = fixture.get("signed_chain_C0_Z1_differences")
    require(isinstance(chain, list) and len(chain) == 2 and
            all(mp.mpf(str(value)) < mp.mpf("1e-38") for value in chain),
            "fixture signed chain differences changed")
    return {
        "passed": True,
        "finite_N": fixture.get("finite_N"),
        "interface_count": fixture.get("interface_count"),
        "control_count": fixture.get("control_count"),
        "maximum_supporting_difference": fixture.get("maximum_supporting_difference"),
        "power_transport_with_nonzero_mean": True,
    }


def _check_source(source: Any, family: Any) -> dict[str, Any]:
    require(isinstance(source, Mapping) and source.get("passed") is True,
            "six actual outer source joins are not passed")
    require(source.get("source_only") is True and source.get("source_values_selected") is False and
            source.get("numeric_caps_used") is False, "source/value boundary changed")
    require(source.get("source_family", {}).get("matches") is True,
            "outer source family mismatch")
    require(source.get("source_prefix", {}).get("passed") is True,
            "outer source graph prefix changed")
    seams = source.get("seams")
    require(isinstance(seams, Mapping) and set(SEAMS) <= set(seams), "six source seams missing")
    rows = {"U": [], "Uy": [], "Q": [], "Qy": []}
    seam_checks = {}
    for seam in SEAMS:
        row = seams[seam]
        require(isinstance(row, Mapping) and row.get("passed") is True, (seam, "source seam"))
        left, right, left_query, right_query = JOIN_EXPECTED[seam]
        require(row.get("left_chart") == left and row.get("right_chart") == right and
                row.get("left_native_query") == left_query and row.get("right_native_query") == right_query,
                (seam, "endpoint scope"))
        expected_lengths = {"U_rows": 3, "Uy_rows": 2, "Q_rows": 3, "Qy_rows": 2}
        for key, count in expected_lengths.items():
            values = row.get(key)
            require(isinstance(values, list) and len(values) == count and
                    all(isinstance(value, Mapping) and value.get("passed") is True for value in values),
                    (seam, key, "actual source rows"))
            rows[key.removesuffix("_rows")].extend(values)
        seam_checks[seam] = {key: len(row[key]) for key in expected_lengths}
    require(len(rows["U"]) == 18 and len(rows["Uy"]) == 12 and
            len(rows["Q"]) == 18 and len(rows["Qy"]) == 12,
            "source row coverage changed")
    charts = source.get("charts")
    require(isinstance(charts, Mapping), "source chart scopes missing")
    for name in ("Rm_patch", "Rh_reference", "O2_slope", "O2_axial", "O2_buffer", "O3_transition", "O3_power"):
        require(isinstance(charts.get(name), Mapping) and charts[name].get("passed") is True,
                (name, "source chart"))
    require(source.get("limits", {}).get("source_values_selected") is False and
            source.get("limits", {}).get("opaque_global_alias") is False,
            "source limits changed")
    return {"passed": True, "seams": seam_checks, "U_rows": 18, "Uy_rows": 12,
            "Q_rows": 18, "Qy_rows": 12, "cap_as_value": False}


def _check_native(native: Any) -> dict[str, Any]:
    require(isinstance(native, Mapping) and native.get("passed") is True,
            "native source guards are not passed")
    require(native.get("source_only") is True and native.get("numeric_caps_used") is False and
            native.get("source_values_selected") is False and native.get("ancestor_replayed") is False,
            "native source/value boundary changed")
    required = native.get("checks", {}).get("required")
    require(isinstance(required, Mapping) and all(value is True for value in required.values()),
            "native required source guards incomplete")
    premises = native.get("premises", {})
    for key in ("d_native_log_R_equals_J", "physical_phase_uses_common_N",
                "physical_U_is_Pstar_V_N", "physical_Uy_uses_slow_fast_once",
                "physical_Qy_R_once", "variable_J_not_replaced_by_constant"):
        require(premises.get(key) is True, ("native premise", key))
    common = native.get("checks", {}).get("common_phase", {})
    require(common.get("passed") is True and common.get("N_node") is not None,
            "common N/phase source guard missing")
    charts = native.get("chartguards", {})
    require(isinstance(charts, Mapping), "native chart guards missing")
    for name in ("Rm_patch", "Rh_reference", "O2_slope", "O2_axial", "O2_buffer", "O3_transition", "O3_power"):
        row = charts.get(name)
        require(isinstance(row, Mapping) and row.get("passed") is True, (name, "native chart guard"))
    for name in ("Rm_patch", "O2_axial"):
        require(charts[name].get("native_J_is_varying") is True and
                charts[name].get("d_native_log_R_equals_J") is True,
                (name, "variable J source dependency"))
    seams = native.get("seamgeometry")
    require(isinstance(seams, Mapping) and set(SEAMS) <= set(seams), "native seam guards missing")
    for name in SEAMS:
        require(seams[name].get("passed") is True and seams[name].get("source_values_selected") is False,
                (name, "native seam guard"))
    power = native.get("checks", {}).get("power_transport_certificate", {})
    require(isinstance(power, Mapping), "power transport source guard missing")
    return {"passed": True, "common_N_node": common.get("N_node"),
            "variable_J": True, "power_transport_passed": power.get("passed") is True,
            "power_transport": power}


def _check_mean(mean: Any) -> dict[str, Any]:
    require(isinstance(mean, Mapping) and mean.get("passed") is True,
            "mean/history source certificate is not passed")
    require(mean.get("source_only") is True and mean.get("source_values_selected") is False and
            mean.get("numeric_caps_used") is False and mean.get("ancestor_replayed") is False,
            "mean source/value boundary changed")
    modified = mean.get("modified_physical_mean_rows")
    histories = mean.get("finite_N_history_rows")
    require(isinstance(modified, Mapping) and isinstance(histories, Mapping),
            "mean/history rows missing")
    for seam in SEAMS:
        require(modified.get(seam, {}).get("passed") is True and
                histories.get(seam, {}).get("passed") is True,
                (seam, "mean/history seam"))
    leading = mean.get("leading_mean_rows", {})
    require(leading.get("passed") is True and leading.get("original_leading_output_alias_issued") is False and
            leading.get("original_leading_controls_opaque") is True,
            "leading mean alias boundary changed")
    leading_qy = mean.get("original_leading_Qy_input_rows")
    require(isinstance(leading_qy, Mapping) and leading_qy.get("passed") is True and
            leading_qy.get("row_count") == 96 and leading_qy.get("radial_row_count") == 24 and
            leading_qy.get("source_node_alias_issued") is False and
            leading_qy.get("source_premises_only") is True,
            "independent leading Qy input provenance changed")
    leading_qy_inputs = mean.get("leading_Qy_inputs")
    require(isinstance(leading_qy_inputs, Mapping) and set(SEAMS) <= set(leading_qy_inputs) and
            all(isinstance(leading_qy_inputs[seam], Mapping) and
                leading_qy_inputs[seam].get("passed") is True for seam in SEAMS),
            "per-seam leading Qy inputs are not source-bound")
    trace = mean.get("trace_contract", {})
    require(trace.get("common_N_and_rate_source_retained") is True and
            trace.get("variable_native_Jacobians_retained") is True and
            trace.get("transition_power_nonzero_mean_retained") is True,
            "mean trace contract changed")
    return {"passed": True, "seams": list(SEAMS), "leading_Qy_common": True,
            "nonzero_transition_power_mean": True, "cap_as_value": False}


def _check_recovery(recovery: Any) -> dict[str, Any]:
    require(isinstance(recovery, Mapping) and recovery.get("passed") is True,
            "strict six-seam recovery proof is not passed")
    for key in ("source_values_selected", "numeric_caps_used"):
        if key in recovery:
            require(recovery.get(key) is False, ("recovery source/value boundary", key))
    candidate = recovery.get("joins") if isinstance(recovery.get("joins"), Mapping) else recovery.get("seams")
    require(isinstance(candidate, Mapping), "per-seam recovery rows missing")
    rows = {}
    for seam in SEAMS:
        row = candidate.get(seam)
        require(isinstance(row, Mapping) and row.get("passed") is True, (seam, "recovery"))
        trace_rows = row.get("trace_rows")
        require(isinstance(trace_rows, list) and len(trace_rows) == 15,
                (seam, "recovery trace row count"))
        quantities = [trace.get("quantity") for trace in trace_rows if isinstance(trace, Mapping)]
        require(len(quantities) == 15 and
                quantities.count("U") == 3 and quantities.count("Uy") == 2 and
                quantities.count("m") == 3 and quantities.count("Q") == 3 and
                quantities.count("Qy") == 2 and quantities.count("RQy") == 2,
                (seam, "recovery trace quantities"))
        leading_qy = row.get("leading_Qy_common")
        require(isinstance(leading_qy, Mapping) and leading_qy.get("passed") is True,
                (seam, "leading Qy common recovery"))
        d_rows = row.get("D_equality_rows")
        require(isinstance(d_rows, list) and len(d_rows) == 2 and
                sorted(d.get("ordinary_Z_order") for d in d_rows if isinstance(d, Mapping)) == [0, 1] and
                all(isinstance(d, Mapping) and d.get("passed") is True for d in d_rows),
                (seam, "physical D equality recovery"))
        require(row.get("endpoint_wrapper_binding_passed") is True,
                (seam, "scoped endpoint wrapper recovery"))
        rows[seam] = {key: row.get(key) for key in ("left_chart", "right_chart", "left_endpoint", "right_endpoint")
                      if key in row}
    return {"passed": True, "seams": rows, "strict_per_seam": True, "cap_as_value": False}


def _check_budget(sample: Mapping[str, Any], *, power_transport_used: bool) -> dict[str, Any]:
    budget = sample.get("external_endpoint_budget_witness")
    require(isinstance(budget, Mapping), "outer endpoint budget missing")
    own = budget.get("own_longreshape0")
    require(isinstance(own, Mapping), "own long native0 budget missing")
    own_rows = own.get("radius_weighted_C0_C1_log_bounds")
    require(isinstance(own_rows, list) and len(own_rows) == 2, "own long0 C0/Z1 rows changed")
    for row in own_rows:
        _check_record(row, "own long0 row")
    aggregate = budget.get("outer_chain_upper_aggregate")
    require(isinstance(aggregate, Mapping), "outer-chain upper aggregate missing")
    atom = aggregate.get("aggregate_log_atom")
    _check_record(atom, "outer-chain aggregate atom")
    aggregate_join = "buffer_transition" if power_transport_used else "transition_power"
    aggregate_spec = PACKET_SPECS[aggregate_join]
    expected_index = aggregate_spec["raw_index"]
    aggregate_term = aggregate.get("aggregate", {})
    require(aggregate.get("raw_index") == expected_index and
            aggregate.get("chart") == aggregate_spec["chart"] and
            aggregate.get("aggregate_coordinate") == aggregate_spec["aggregate_coordinate"] and
            aggregate.get("term_name") == aggregate_spec["term_name"] and
            aggregate.get("source_report") == aggregate_spec["report"] and
            isinstance(aggregate_term, Mapping) and
            aggregate_term.get("name") == aggregate_spec["term_name"] and
            aggregate_term.get("log_scale_terms") == [atom] and
            atom.get("positive_error_bound_atom_only") is True and
            atom.get("source_value_selected") is False,
            "outer aggregate scope changed")
    inputs = budget.get("positive_union_inputs")
    require(isinstance(inputs, list) and len(inputs) == 3 and inputs[:2] == own_rows and inputs[2] == atom,
            "own long0 plus outer aggregate inputs changed")
    union = _check_positive_union(inputs, budget.get("positive_union"), "long native0 plus outer aggregate")
    require(budget.get("no_cap_subtraction_value_midpoint_or_overlap") is True and
            budget.get("power_transport_certificate_used") is power_transport_used,
            "outer positive budget contract changed")
    return {"passed": True, "aggregate_raw_index": expected_index, **union, "cap_as_value": False}


def _check_samples(report: Mapping[str, Any], parent: Mapping[str, Any], *, power_transport_used: bool) -> list[dict[str, Any]]:
    samples = report.get("samples")
    parent_samples = {row.get("Z"): row for row in parent.get("samples", []) if isinstance(row, Mapping)}
    require(isinstance(samples, list) and tuple(row.get("Z") for row in samples) == SAVED_Z,
            "saved outer-chain Z set changed")
    checks = []
    for sample in samples:
        require(isinstance(sample, Mapping) and sample.get("Z") in parent_samples, ("sample", sample.get("Z")))
        base = parent_samples[sample["Z"]]
        raw = sample.get("raw_m15_C0_C1_terms")
        previous = base.get("correlated_reduced_m15_C0_C1_terms")
        reduced = sample.get("correlated_reduced_m15_C0_C1_terms")
        require(isinstance(raw, list) and len(raw) == RAW_COUNT and raw == base.get("raw_m15_C0_C1_terms"),
                (sample.get("Z"), "raw 213 changed"))
        require(isinstance(previous, list) and len(previous) == PREVIOUS_REDUCED_COUNT and
                isinstance(reduced, list) and len(reduced) == NEW_REDUCED_COUNT,
                (sample.get("Z"), "reduced counts changed"))
        old_external = [row for row in previous if row.get("name") == OLD_EXTERNAL_TERM]
        new_external = [row for row in reduced if row.get("name") == NEW_EXTERNAL_TERM]
        removed = [row for row in previous if row.get("name") in REMOVED_OUTER_TERMS]
        require(len(old_external) == 1 and len(new_external) == 1 and
                len(removed) == len(REMOVED_OUTER_TERMS) and
                {row.get("name") for row in removed} == set(REMOVED_OUTER_TERMS),
                (sample.get("Z"), "outer term identities changed"))
        old_keep = [row for row in previous if row.get("name") not in set(REMOVED_OUTER_TERMS) | {OLD_EXTERNAL_TERM}]
        new_keep = [row for row in reduced if row.get("name") != NEW_EXTERNAL_TERM]
        require(len(old_keep) == UNCHANGED_COUNT and new_keep == old_keep,
                (sample.get("Z"), "other 189 sectors changed"))
        require(sum(row.get("name") == "reserved_n1_own_swirl_once" for row in reduced) == 1,
                (sample.get("Z"), "reserved term count"))
        require(sample.get("independent_target") == base.get("independent_target"),
                (sample.get("Z"), "independent target changed"))
        contract = sample.get("contract", {})
        for key in (
            "raw_213_verbatim", "other_189_verbatim", "independent_target_unchanged",
            "reserved_retained_once", "ordinary_Z1_and_physical_R_once",
            "leading_subtraction_retained", "all_six_outer_source_joins_proved",
            "six_signed_outer_viscous_sectors_removed", "outer_upper_aggregate_selected_as_error_only",
            "flat_U_does_not_zero_history_mean_or_Qy", "no_error_cap_as_value",
            "all_other_volume_width_numerical_errors_retained", "ancestors_not_reexecuted",
            "source_bound_signed_trace_contraction",
        ):
            require(contract.get(key) is True, (sample.get("Z"), "sample contract", key))
        sample_limits = sample.get("limits", {})
        require(isinstance(sample_limits, Mapping) and
                sample_limits.get("outer_chain_internal_trace_resolved") is True and
                all(value is False for key, value in sample_limits.items()
                    if key != "outer_chain_internal_trace_resolved"),
                (sample.get("Z"), "sample limits"))
        checks.append({"Z": sample.get("Z"), "raw_213": True, "reduced_190": True,
                       "unchanged_189": True, "budget": _check_budget(sample, power_transport_used=power_transport_used),
                       "reserved_once": True, "target_unchanged": True})
    return checks


def _check_packet(report: Mapping[str, Any], parent: Mapping[str, Any], fixture_summary: Mapping[str, Any] | None = None) -> dict[str, Any]:
    require(report.get("source_family") == parent.get("source_family"), "outer source family changed")
    require(report.get("parent_report") == PARENT and report.get("parent_report_sha256") == PARENT_SHA,
            "outer parent identity changed")
    source_check = _check_source(report.get("source_derivative_theorem"), report.get("source_family"))
    native_check = _check_native(report.get("native_source_guards"))
    mean_check = _check_mean(report.get("mean_source_certificate"))
    recovery_check = _check_recovery(report.get("source_seam_proofs"))
    power = report.get("power_transport_certificate")
    power_used = isinstance(power, Mapping) and power.get("passed") is True
    if power_used:
        power_flags = [
            key for key, value in power.items()
            if value is True and any(token in str(key).lower().replace("_", "")
                                     for token in ("constant", "equal", "identical", "lowerupper", "dpower"))
        ]
        require(power_flags, "power 202 aggregate lacks an explicit transport identity")
    require(report.get("numerical_source_binding", {}).get("power_transport_certificate_used") is power_used,
            "power transport usage changed")
    samples = _check_samples(report, parent, power_transport_used=power_used)
    contract = report.get("contract", {})
    require(contract == {
        "accepted_ancestor_constructors_not_run": True,
        "source_functions_not_cap_values": True,
        "raw_213_and_independent_target_retained": True,
    }, "top contract changed")
    limits = report.get("limits", {})
    require(isinstance(limits, Mapping) and all(value is False for value in limits.values()),
            "global completion limits changed")
    if fixture_summary is not None:
        require(fixture_summary.get("passed") is True, "fixture summary not passed")
    return {"source": source_check, "native": native_check, "mean": mean_check,
            "recovery": recovery_check, "samples": samples,
            "power_transport_used": power_used, "global_limits_all_false": True,
            "cap_as_value": False}


def _load_parent() -> Mapping[str, Any]:
    path = DIRECTORY / PARENT
    require(path.is_file() and sha(path) == PARENT_SHA, "accepted Rm inlet parent changed")
    parent = read_json(path)
    require(isinstance(parent, Mapping) and parent.get("checks", {}).get("passed") is True,
            "accepted Rm inlet parent is not passed")
    return parent


def run() -> dict[str, Any]:
    started = time.perf_counter()
    # This is the sole write-mode fixture invocation.
    from check_original_outer_chain_m15 import check_fixture
    fixture = check_fixture()
    fixture_summary = _check_fixture_shape(fixture)
    from openai_ns_reconstruction.original_outer_chain_m15 import OriginalOuterChainM15
    owner = OriginalOuterChainM15(precision=150)
    packet = owner.build()
    code = dict(owner.code_hashes)
    code[FIXTURE_NAME] = sha(ROOT / FIXTURE_NAME)
    code[RUNNER_NAME] = sha(ROOT / RUNNER_NAME)
    report = dict(packet, input_hashes=dict(owner.hashes), implementation_hashes=code)
    parent = _load_parent()
    packet_check = _check_packet(report, parent, fixture_summary)
    report["checks"] = {
        "passed": True,
        "independent_fixture": fixture,
        "fixture_invariants": fixture_summary,
        "packet_invariants": packet_check,
        "saved_Z": list(SAVED_Z),
        "raw_213_verbatim": True,
        "previous_reduced_196_to_new_190": True,
        "unchanged_other_sectors_189": True,
        "six_actual_source_joins_U3_Uy2_Q3_Qy2": True,
        "native_variable_J_common_phase_and_R_once": True,
        "mean_and_leading_Qy_common_six": True,
        "strict_per_seam_recovery": True,
        "own_long0_plus_outer_aggregate_positive_union": True,
        "power_202_requires_explicit_transport_certificate": True,
        "flat_U_does_not_zero_mean": True,
        "source_functions_not_caps": True,
        "production_is_new_owner_and_fixture_once": True,
    }
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path: Path) -> dict[str, Any]:
    """Inspect only serialized report shape and pinned files."""
    data = read_json(path)
    require(isinstance(data, Mapping) and data.get("checks", {}).get("passed") is True,
            "saved outer-chain report is not passed")
    check_hash_group(data.get("input_hashes"), "input", DIRECTORY)
    check_hash_group(data.get("implementation_hashes"), "implementation", ROOT)
    parent = _load_parent()
    fixture = data.get("checks", {}).get("independent_fixture")
    fixture_summary = _check_fixture_shape(fixture)
    packet_check = _check_packet(data, parent, fixture_summary)
    expected = {
        "saved_Z": list(SAVED_Z),
        "raw_213_verbatim": True,
        "previous_reduced_196_to_new_190": True,
        "unchanged_other_sectors_189": True,
        "six_actual_source_joins_U3_Uy2_Q3_Qy2": True,
        "native_variable_J_common_phase_and_R_once": True,
        "mean_and_leading_Qy_common_six": True,
        "strict_per_seam_recovery": True,
        "own_long0_plus_outer_aggregate_positive_union": True,
        "power_202_requires_explicit_transport_certificate": True,
        "flat_U_does_not_zero_mean": True,
        "source_functions_not_caps": True,
    }
    for key, value in expected.items():
        require(data.get("checks", {}).get(key) == value, ("saved check", key))
    return {"checks": data["checks"], "fixture_invariants": fixture_summary,
            "packet_invariants": packet_check, "execution_seconds": data.get("execution_seconds")}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group(required=True)
    modes.add_argument("--write-report", action="store_true")
    modes.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = DIRECTORY / NAME
    if args.write_report:
        data = run()
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n",
                             encoding="utf-8", newline="\n")
        temporary.replace(path)
    else:
        require(path.is_file(), ("missing report", str(path)))
        data = inspect(path)
    print(json.dumps({"checks": data.get("checks", {}),
                      "execution_seconds": data.get("execution_seconds")}, indent=2), flush=True)
