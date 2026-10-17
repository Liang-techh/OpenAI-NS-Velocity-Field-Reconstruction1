"""Write or inspect the source-bound compact/quiet tail m15 packet.

Write mode constructs the new tail owner and the independent local-germ
fixture once.  Inspect mode reads only the serialized report and pinned hash
maps; it never constructs an owner, imports a fixture, or replays an ancestor.
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


NAME = PREFIX + "compliant_current_original_tail_m15_numeric.json"
FIXTURE_NAME = "scripts/check_original_tail_m15.py"
RUNNER_NAME = "scripts/run_original_tail_m15.py"
PARENT = PREFIX + "compliant_current_original_outer_chain_m15_numeric.json"
PARENT_SHA = "cf185c36affc4183984aff5ef9caf119c68ea90d0b8e95535b1f8698e2161e53"
SAVED_Z = ("-.3", "0", ".3")
RAW_COUNT = 213
PARENT_REDUCED_COUNT = 190
NEW_REDUCED_COUNT = 189
UNCHANGED_COUNT = 188
COMPACT_RAW_INDEX = 210
# This is the external endpoint term produced by the accepted current
# outer-chain parent.  The parent's provenance field intentionally names an
# older ancestor term and must not be used as the replacement identity here.
CURRENT_PARENT_EXTERNAL_TERM = "long_outer_chain_finite_N_viscous_external_endpoints"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value: Any) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def require(condition: bool, message: Any) -> None:
    if not condition:
        raise AssertionError(message)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _resolve(name: str, base: Path | None = None) -> Path:
    raw = Path(name)
    if raw.is_absolute():
        return raw
    candidates = (ROOT / raw, (base or DIRECTORY) / raw,
                  ROOT / "experiments" / "root_st073" / raw)
    for candidate in candidates:
        if candidate.is_file():
            return candidate
    return candidates[0]


def check_hash_group(group: Any, label: str, base: Path | None = None) -> None:
    require(isinstance(group, Mapping) and group, (label, "hash map missing"))
    for name, expected in group.items():
        require(isinstance(name, str) and isinstance(expected, str) and len(expected) == 64,
                (label, name, "malformed digest"))
        path = _resolve(name, base)
        require(path.is_file() and sha(path) == expected, (label, name, "pinned file changed"))


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
            (label, "not an error-only positive atom"))


def _check_positive_union(values: Any, record: Any, label: str) -> dict[str, Any]:
    require(isinstance(values, list) and values, (label, "union inputs missing"))
    _check_record(record, label)
    ctx, iv = _ctx_pair()
    decoded = [_decode_record(value, ctx, iv) for value in values]
    top = max(_bounds(ctx, value)[1] for value in decoded)
    expected = iv.mpf(top) + iv.log(len(decoded))
    require(_record_tuples(record) == [list(value) for value in expected._mpi_],
            (label, "serialized MPF union changed"))
    return {"components": len(decoded), "exact_mpf_union": True, "cap_as_value": False}


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _find_flag(value: Any, names: set[str]) -> bool:
    if isinstance(value, Mapping):
        for key, item in value.items():
            if key in names and item is True:
                return True
            if isinstance(item, (Mapping, list, tuple)) and _find_flag(item, names):
                return True
    elif isinstance(value, (list, tuple)):
        return any(_find_flag(item, names) for item in value)
    return False


def _check_fixture_shape(fixture: Any) -> dict[str, Any]:
    require(isinstance(fixture, Mapping) and fixture.get("passed") is True,
            "independent tail fixture failed")
    expected = {
        "synthetic_local_germs_only": True,
        "original_source_values_evaluated": False,
        "original_Poisson_inverse_evaluated": False,
        "NS_residual": False,
        "finite_N": 17,
        "interface_count": 2,
        "quiet_constant_D": True,
        "quiet_mean_nonzero_decay": True,
    }
    for key, value in expected.items():
        require(fixture.get(key) == value, ("fixture flag", key, fixture.get(key)))
    controls = fixture.get("controls")
    expected_controls = {
        "extra_N_on_fast_term_rejected",
        "varying_J_first_derivative_required",
        "varying_J_second_derivative_required",
        "phase_reset_rejected",
        "nonzero_incoming_mean",
        "quiet_mean_not_zeroed",
        "quiet_D_constant_while_M_decays",
        "leading_subtraction_retained",
        "physical_R_required",
        "ordinary_Z1_nonzero",
        "compact_endpoint_scope_retained",
        "flat_beta_cutoff_derivatives",
    }
    require(isinstance(controls, Mapping) and set(controls) == expected_controls and
            all(value is True for value in controls.values()), ("fixture controls", controls))
    require(fixture.get("control_count") == len(expected_controls), "fixture control count changed")
    terminal_case = fixture.get("terminal_relative_zero_case")
    require(isinstance(terminal_case, Mapping) and terminal_case.get("passed") is True and
            terminal_case.get("relative_D_zero") is True and
            terminal_case.get("absolute_leading_Qy_nonzero") is True and
            terminal_case.get("absolute_mean_not_zeroed") is True,
            "terminal relative-zero fixture collapsed absolute history")
    require(mp.mpf(str(fixture.get("maximum_supporting_difference"))) < mp.mpf("1e-70"),
            "fixture tolerance changed")
    chain = fixture.get("signed_chain_C0_Z1_differences")
    require(isinstance(chain, list) and len(chain) == 2 and
            all(mp.mpf(str(value)) < mp.mpf("1e-70") for value in chain),
            "fixture chain differences changed")
    return {
        "passed": True,
        "finite_N": fixture.get("finite_N"),
        "interface_count": fixture.get("interface_count"),
        "control_count": fixture.get("control_count"),
        "maximum_supporting_difference": fixture.get("maximum_supporting_difference"),
    }


def _check_source_certificates(report: Mapping[str, Any]) -> dict[str, Any]:
    source = report.get("source_derivative_theorem")
    native = report.get("native_source_guards")
    seams = report.get("source_seam_proofs")
    for label, cert in (("source derivative", source), ("native", native), ("seams", seams)):
        require(isinstance(cert, Mapping) and cert.get("passed") is True, f"{label} certificate failed")
        for key in ("source_values_selected", "numeric_caps_used", "ancestor_replayed"):
            if key in cert:
                require(cert.get(key) is False, (label, key, "must remain false"))
    seam_map = seams.get("joins") if isinstance(seams.get("joins"), Mapping) else seams.get("seams")
    require(isinstance(seam_map, Mapping), "tail source seam map missing")
    selected = {}
    for name, row in seam_map.items():
        if not isinstance(row, Mapping) or row.get("passed") is not True:
            continue
        selected[name] = row
    require(len(selected) >= 2, ("tail source interfaces", list(selected)))
    for name in ("power_compact", "compact_quiet"):
        require(name in selected, (name, "required tail source seam missing"))
    # Require the native helper's concrete schema.  Recursive flag searches
    # previously accepted unrelated diagnostic booleans and also looked for
    # names the helper never emits (for example ``common_phase``).
    native_checks = _mapping(native.get("checks"))
    required = native_checks.get("required")
    require(isinstance(required, Mapping) and required and
            all(value is True for value in required.values()),
            "native required source guards incomplete")
    compact_radius = _mapping(_mapping(native.get("radius_guards")).get("compact"))
    require(compact_radius.get("native_J_is_varying") is True and
            compact_radius.get("native_J_is_constant") is False and
            compact_radius.get("d_native_log_R_equals_J") is True,
            "varying compact Jacobian source claim missing")
    common_n = _mapping(native_checks.get("common_N"))
    require(common_n.get("passed") is True and
            common_n.get("outer") == common_n.get("compact") == common_n.get("quiet") == 18,
            "common N source claim missing")
    premises = _mapping(native.get("premises"))
    require(premises.get("common_N_and_global_phase") is True and
            premises.get("physical_R_once") is True,
            "native phase/R source premises missing")
    return {"passed": True, "interfaces": sorted(selected), "native_source_only": True}


def _check_recovery(report: Mapping[str, Any]) -> dict[str, Any]:
    recovery = report.get("source_seam_proofs")
    seam_map = recovery.get("joins") if isinstance(recovery, Mapping) and isinstance(recovery.get("joins"), Mapping) else _mapping(recovery).get("seams")
    require(isinstance(seam_map, Mapping), "tail recovery seams missing")
    rows = {}
    for name in ("power_compact", "compact_quiet"):
        row = seam_map.get(name)
        require(isinstance(row, Mapping) and row.get("passed") is True, (name, "recovery"))
        trace = row.get("trace_rows")
        require(isinstance(trace, list) and len(trace) == 15, (name, "trace row count"))
        quantities = [item.get("quantity") for item in trace if isinstance(item, Mapping)]
        require(len(quantities) == 15 and quantities.count("U") == 3 and quantities.count("Uy") == 2 and
                quantities.count("m") == 3 and quantities.count("Q") == 3 and
                quantities.count("Qy") == 2 and quantities.count("RQy") == 2,
                (name, "trace quantities"))
        require(_mapping(row.get("leading_Qy_common")).get("passed") is True,
                (name, "leading Qy common"))
        d_rows = row.get("D_equality_rows")
        require(isinstance(d_rows, list) and len(d_rows) == 2 and
                sorted(item.get("ordinary_Z_order") for item in d_rows if isinstance(item, Mapping)) == [0, 1] and
                all(isinstance(item, Mapping) and item.get("passed") is True for item in d_rows),
                (name, "D equality rows"))
        require(row.get("endpoint_wrapper_binding_passed") is True,
                (name, "endpoint wrapper binding"))
        rows[name] = {"trace_rows": 15, "D_orders": [0, 1]}
    return {"passed": True, "interfaces": rows, "strict_per_seam": True}


def _terminal_certificate(report: Mapping[str, Any]) -> tuple[bool, Mapping[str, Any] | None]:
    value = report.get("terminal_D_zero_certificate")
    if value is None:
        return False, None
    require(isinstance(value, Mapping) and value.get("passed") is True,
            "terminal D-zero certificate failed")
    for key in ("source_only", "D_upper_zero", "ordinary_Z1", "physical_R_once", "absolute_mean_not_zeroed"):
        require(value.get(key) is True, ("terminal D-zero certificate", key))
    for key in ("source_values_selected", "numeric_caps_used", "ancestor_replayed"):
        if key in value:
            require(value.get(key) is False, ("terminal D-zero certificate", key))
    mass = report.get("compact_mass_row_identity")
    require(isinstance(mass, Mapping) and mass.get("passed") is True,
            "compact terminal mass identity missing")
    require(mass.get("m_relative_zero_source_identity") is True and
            mass.get("absolute_modified_m_equals_leading_m") is True and
            mass.get("absolute_leading_histories_not_zeroed") is True,
            "compact terminal mass source identity incomplete")
    for key in ("numeric_value_selected", "source_values_selected"):
        if key in mass:
            require(mass.get(key) is False, ("compact terminal mass identity", key))
    quiet = report.get("quiet_transport_certificate")
    if quiet is not None:
        require(isinstance(quiet, Mapping) and quiet.get("passed") is True,
                "quiet transport certificate failed")
        require(_find_flag(quiet, {"absolute_mean_not_zeroed", "quiet_constant_D", "nonzero_mean_retained"}),
                "quiet transport nonzero-mean source claim missing")
    return True, value


def _aggregate_from_budget(budget: Mapping[str, Any]) -> Mapping[str, Any] | None:
    # Only the canonical active packet is an admissible budget component.
    # Terminal-zero packets retain raw-index-210 data as provenance, which
    # must never be mistaken for an active positive-union aggregate.
    value = budget.get("compact_upper_aggregate")
    return value if isinstance(value, Mapping) else None


def _check_budget(sample: Mapping[str, Any], terminal_zero: bool, compact_name: str) -> dict[str, Any]:
    budget = sample.get("external_endpoint_budget_witness")
    require(isinstance(budget, Mapping), "tail external budget missing")
    own = budget.get("own_longreshape0") or budget.get("own_long0")
    require(isinstance(own, Mapping), "own long0 budget missing")
    own_rows = own.get("radius_weighted_C0_C1_log_bounds")
    require(isinstance(own_rows, list) and len(own_rows) == 2, "own long0 C0/Z1 rows changed")
    for row in own_rows:
        _check_record(row, "own long0")
    require(budget.get("compact210_aggregate_used") is (not terminal_zero),
            "compact aggregate usage mode changed")
    aggregate = _aggregate_from_budget(budget)
    if terminal_zero:
        require(aggregate is None, "terminal D-zero budget must not retain compact aggregate")
        terminal_budget_cert = budget.get("terminal_D_zero_certificate")
        require(isinstance(terminal_budget_cert, Mapping) and
                terminal_budget_cert.get("passed") is True,
                "terminal D-zero budget certificate missing")
        inputs = budget.get("positive_union_inputs")
        require(isinstance(inputs, list) and len(inputs) == 2 and inputs == own_rows,
                "terminal D-zero union must contain only own long0 rows")
    else:
        require(isinstance(aggregate, Mapping), "compact upper aggregate missing")
        atom = aggregate.get("aggregate_log_atom")
        _check_record(atom, "compact aggregate")
        term = aggregate.get("aggregate", {})
        require(aggregate.get("raw_index") == COMPACT_RAW_INDEX and
                compact_name == aggregate.get("term_name") and
                isinstance(term, Mapping) and term.get("name") == compact_name and
                term.get("log_scale_terms") == [atom],
                "compact aggregate packet scope changed")
        inputs = budget.get("positive_union_inputs")
        require(isinstance(inputs, list) and len(inputs) == 3 and inputs[:2] == own_rows and inputs[2] == atom,
                "long0 plus compact aggregate inputs changed")
    union = _check_positive_union(inputs, budget.get("positive_union"), "tail positive union")
    require(budget.get("no_cap_subtraction_value_midpoint_or_overlap") is True,
            "tail budget subtraction contract changed")
    return {"passed": True, "components": union["components"], "terminal_D_zero": terminal_zero,
            "aggregate_raw_index": None if terminal_zero else COMPACT_RAW_INDEX, "cap_as_value": False}


def _check_samples(report: Mapping[str, Any], parent: Mapping[str, Any], terminal_zero: bool) -> list[dict[str, Any]]:
    samples = report.get("samples")
    parent_samples = {row.get("Z"): row for row in parent.get("samples", []) if isinstance(row, Mapping)}
    require(isinstance(samples, list) and tuple(row.get("Z") for row in samples) == SAVED_Z,
            "tail saved Z set changed")
    out = []
    compact_names = set()
    for sample in samples:
        require(isinstance(sample, Mapping) and sample.get("Z") in parent_samples, ("sample", sample.get("Z")))
        base = parent_samples[sample["Z"]]
        raw = sample.get("raw_m15_C0_C1_terms")
        previous = base.get("correlated_reduced_m15_C0_C1_terms")
        reduced = sample.get("correlated_reduced_m15_C0_C1_terms")
        require(isinstance(raw, list) and len(raw) == RAW_COUNT and raw == base.get("raw_m15_C0_C1_terms"),
                (sample.get("Z"), "raw 213 changed"))
        require(isinstance(previous, list) and len(previous) == PARENT_REDUCED_COUNT and
                isinstance(reduced, list) and len(reduced) == NEW_REDUCED_COUNT,
                (sample.get("Z"), "tail term counts changed"))
        compact_name = raw[COMPACT_RAW_INDEX].get("name")
        require(isinstance(compact_name, str) and "compact" in compact_name.lower() and
                "viscous" in compact_name.lower(), (sample.get("Z"), "compact viscous raw term"))
        compact_names.add(compact_name)
        # The accepted parent has exactly one current external replacement.
        # Do not read `previous_external_replacement_retained_as_provenance`:
        # that field deliberately records the older Rm/Rh ancestor term, not
        # the term present in this parent's 190-term reduced view.
        old_external_name = CURRENT_PARENT_EXTERNAL_TERM
        require(sum(row.get("name") == old_external_name for row in previous) == 1,
                (sample.get("Z"), "current parent external term"))
        old_keep = [row for row in previous if row.get("name") not in {old_external_name, compact_name}]
        require(len(old_keep) == UNCHANGED_COUNT, (sample.get("Z"), "parent unchanged sector count"))
        # Derive the sole replacement by exact list subtraction rather than
        # assuming a new product name.
        old_keep_names = {row.get("name") for row in old_keep}
        new_external = [row for row in reduced if row.get("name") not in old_keep_names]
        old_external = next(row for row in previous if row.get("name") == old_external_name)
        require(len(new_external) == 1 and new_external[0].get("name") != compact_name and
                (new_external[0].get("name") != old_external_name or new_external[0] != old_external),
                (sample.get("Z"), "new tail external replacement"))
        require([row for row in previous if row.get("name") == old_external_name].__len__() == 1 and
                [row for row in previous if row.get("name") == compact_name].__len__() == 1,
                (sample.get("Z"), "removed term identities"))
        new_keep = [row for row in reduced if row.get("name") in old_keep_names]
        require(new_keep == old_keep, (sample.get("Z"), "other 188 sectors changed"))
        require(len(new_keep) == UNCHANGED_COUNT, (sample.get("Z"), "tail unchanged sector count"))
        require(sum(row.get("name") == "reserved_n1_own_swirl_once" for row in reduced) == 1,
                (sample.get("Z"), "reserved term count"))
        require(sample.get("independent_target") == base.get("independent_target"),
                (sample.get("Z"), "independent target changed"))
        contract = sample.get("contract", {})
        required_contract = {
            "raw_213_verbatim", "previous_190_to_new_189", "other_188_verbatim",
            "six_signed_inner_lower_retained", "independent_target_unchanged",
            "reserved_retained_once", "quiet_mean_square_retained_once",
            "ordinary_Z1_and_physical_R_once", "leading_subtraction_retained",
            "power_compact_source_join_proved", "compact_quiet_source_join_proved",
            "compact210_aggregate_retained", "source_functions_not_cap_values",
            "ancestors_not_reexecuted",
        }
        for key in required_contract:
            require(contract.get(key) is True, (sample.get("Z"), "sample contract", key))
        require(contract.get("compact210_aggregate_used_as_error_only") is (not terminal_zero),
                (sample.get("Z"), "sample compact aggregate mode"))
        require(contract.get("terminal_D_zero_certificate_used") is terminal_zero,
                (sample.get("Z"), "sample terminal certificate mode"))
        quiet_used = isinstance(report.get("quiet_transport_certificate"), Mapping) and \
            report.get("quiet_transport_certificate", {}).get("passed") is True
        require(contract.get("quiet_constant_D_transport_used") is quiet_used,
                (sample.get("Z"), "sample quiet transport mode"))
        require(contract.get("old_outer_full_202_aggregate_used") is False,
                (sample.get("Z"), "sample old aggregate"))
        sample_limits = sample.get("limits", {})
        if isinstance(sample_limits, Mapping):
            require(all(value is False or ("internal" in str(key).lower() and "trace" in str(key).lower())
                       for key, value in sample_limits.items()),
                    (sample.get("Z"), "sample limits"))
        out.append({"Z": sample.get("Z"), "raw_213": True, "reduced_189": True,
                    "unchanged_188": True, "budget": _check_budget(sample, terminal_zero, compact_name),
                    "compact_term": compact_name, "target_unchanged": True})
    require(len(compact_names) == 1, "compact raw term name changed across Z")
    return out


def _check_packet(report: Mapping[str, Any], parent: Mapping[str, Any], fixture_summary: Mapping[str, Any] | None = None) -> dict[str, Any]:
    require(report.get("source_family") == parent.get("source_family"), "tail source family changed")
    source_check = _check_source_certificates(report)
    recovery_check = _check_recovery(report)
    terminal_zero, terminal = _terminal_certificate(report)
    samples = _check_samples(report, parent, terminal_zero)
    binding = _mapping(report.get("numerical_source_binding"))
    for key in ("physical_radius_power", "physical_R_once"):
        if key in binding:
            require(binding.get(key) in (1, True), ("numeric source binding", key))
    if "ordinary_Z_orders" in binding:
        require(binding.get("ordinary_Z_orders") == [0, 1], "numeric ordinary-Z scope changed")
    require(report.get("contract") == {
        "accepted_ancestor_constructors_not_run": True,
        "source_functions_not_cap_values": True,
        "raw_213_and_independent_target_retained": True,
    }, "top contract changed")
    limits = report.get("limits")
    require(isinstance(limits, Mapping) and all(value is False for value in limits.values()),
            "global completion limits changed")
    if fixture_summary is not None:
        require(fixture_summary.get("passed") is True, "fixture summary failed")
    return {"source": source_check, "recovery": recovery_check, "terminal_D_zero": terminal_zero,
            "samples": samples, "global_limits_all_false": True, "cap_as_value": False}


def _load_parent() -> Mapping[str, Any]:
    path = DIRECTORY / PARENT
    require(path.is_file() and sha(path) == PARENT_SHA, "accepted outer-chain parent changed")
    parent = read_json(path)
    require(isinstance(parent, Mapping) and parent.get("checks", {}).get("passed") is True,
            "accepted outer-chain parent is not passed")
    return parent


def run() -> dict[str, Any]:
    started = time.perf_counter()
    from check_original_tail_m15 import check_fixture
    fixture = check_fixture()
    fixture_summary = _check_fixture_shape(fixture)
    from openai_ns_reconstruction.original_tail_m15 import OriginalTailM15
    owner = OriginalTailM15(precision=150)
    packet = owner.build()
    code = dict(getattr(owner, "code_hashes", {}))
    code[FIXTURE_NAME] = sha(ROOT / FIXTURE_NAME)
    code[RUNNER_NAME] = sha(ROOT / RUNNER_NAME)
    report = dict(packet, input_hashes=dict(getattr(owner, "hashes", {})), implementation_hashes=code)
    parent = _load_parent()
    packet_check = _check_packet(report, parent, fixture_summary)
    report["checks"] = {
        "passed": True,
        "independent_fixture": fixture,
        "fixture_invariants": fixture_summary,
        "packet_invariants": packet_check,
        "saved_Z": list(SAVED_Z),
        "raw_213_verbatim": True,
        "parent_reduced_190_to_new_189": True,
        "unchanged_other_sectors_188": True,
        "two_actual_source_interfaces": True,
        "quiet_constant_D_nonzero_mean_retained": True,
        "terminal_D_zero_certificate_scoped": packet_check["terminal_D_zero"],
        "source_functions_not_caps": True,
        "production_is_new_owner_and_fixture_once": True,
    }
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path: Path) -> dict[str, Any]:
    data = read_json(path)
    require(isinstance(data, Mapping) and data.get("checks", {}).get("passed") is True,
            "saved tail report is not passed")
    check_hash_group(data.get("input_hashes"), "input", DIRECTORY)
    check_hash_group(data.get("implementation_hashes"), "implementation", ROOT)
    parent = _load_parent()
    fixture_summary = _check_fixture_shape(data.get("checks", {}).get("independent_fixture"))
    packet_check = _check_packet(data, parent, fixture_summary)
    for key, value in {
        "saved_Z": list(SAVED_Z),
        "raw_213_verbatim": True,
        "parent_reduced_190_to_new_189": True,
        "unchanged_other_sectors_188": True,
        "two_actual_source_interfaces": True,
        "quiet_constant_D_nonzero_mean_retained": True,
        "source_functions_not_caps": True,
    }.items():
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
