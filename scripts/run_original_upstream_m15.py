"""Write or inspect the source-bound upstream m15 contraction.

The accepted tail packet is the only numerical parent.  Write mode constructs
the new owner and the independent local-germ fixture once.  Inspect mode reads
the serialized packet and pinned hashes only; it never imports the fixture or
constructs an owner.
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


NAME = PREFIX + "compliant_current_original_upstream_m15_numeric.json"
FIXTURE_NAME = "scripts/check_original_upstream_m15.py"
RUNNER_NAME = "scripts/run_original_upstream_m15.py"
PARENT = PREFIX + "compliant_current_original_tail_m15_numeric.json"
PARENT_SHA = "a394d3b538e420936e382b36b63dc7405b7a098b57a557c668deba7a2cc60c9a"
SAVED_Z = ("-.3", "0", ".3")
RAW_COUNT = 213
PARENT_REDUCED_COUNT = 189
NEW_REDUCED_COUNT = 187
UNCHANGED_COUNT = 186
SECOND_SWITCH_RAW_INDEX = 68
POST_POWER_RAW_INDEX = 87
OLD_EXTERNAL_TERM = "long_compact2Rc_finite_N_viscous_external_endpoints"
REMOVED_PARENT_TERMS = (
    "second_switch_finite_N_viscous_error",
    "post_power_finite_N_viscous_error",
    OLD_EXTERNAL_TERM,
)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


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


def _find_flag(value: Any, names: set[str]) -> bool:
    if isinstance(value, Mapping):
        return any((key in names and item is True) or
                   (isinstance(item, (Mapping, list, tuple)) and _find_flag(item, names))
                   for key, item in value.items())
    if isinstance(value, (list, tuple)):
        return any(_find_flag(item, names) for item in value)
    return False


def _rows(value: Any, names: tuple[str, ...]) -> Any:
    if not isinstance(value, Mapping):
        return None
    for name in names:
        candidate = value.get(name)
        if isinstance(candidate, list):
            return candidate
    return None


def _quantity(row: Mapping[str, Any]) -> str:
    raw = str(row.get("quantity", row.get("field", row.get("name", ""))))
    return {"U_y": "Uy", "V_y": "Vy", "Q_y": "Qy", "R_Q_y": "RQy",
            "UY": "Uy", "VY": "Vy", "QY": "Qy"}.get(raw, raw)


def _check_source_rows(group: Any, expected: Mapping[str, int], label: str) -> dict[str, Any]:
    require(isinstance(group, list), (label, "source rows missing"))
    require(all(isinstance(row, Mapping) and row.get("passed") is True for row in group),
            (label, "source row failed"))
    counts = {name: sum(_quantity(row) == name for row in group) for name in expected}
    require(counts == dict(expected), (label, counts, expected))
    for row in group:
        if "ordinary_Z_order" in row:
            require(isinstance(row["ordinary_Z_order"], int) and row["ordinary_Z_order"] >= 0,
                    (label, "bad ordinary-Z order"))
    return {"rows": len(group), "quantities": counts}


def _check_pressure_algebra_witness(join: Mapping[str, Any], label: str) -> dict[str, Any]:
    """Validate the serialized exact pressure-jet congruence witness.

    This is deliberately a schema check over source-function identities.  The
    runner never treats the witness as a numerical cap or evaluates its input
    digests.
    """
    witness = join.get("pressure_algebra_witness")
    require(isinstance(witness, Mapping) and witness.get("passed") is True,
            (label, "pressure algebra witness"))
    rows = witness.get("rows")
    require(isinstance(rows, list) and len(rows) == 5,
            (label, "pressure witness row count"))
    expected_orders = {"Q": [0, 1, 2], "Qy": [0, 1]}
    seen: dict[str, list[int]] = {key: [] for key in expected_orders}
    for index, row in enumerate(rows):
        require(isinstance(row, Mapping) and row.get("passed") is True,
                (label, "pressure witness row", index))
        quantity = _quantity(row)
        require(quantity in expected_orders,
                (label, "pressure witness quantity", index, quantity))
        order = row.get("ordinary_Z_order")
        require(isinstance(order, int) and order in expected_orders[quantity],
                (label, "pressure witness ordinary-Z order", index, order))
        require(order not in seen[quantity],
                (label, "duplicate pressure witness order", quantity, order))
        seen[quantity].append(order)
        require(row.get("reduced_difference_monomials") == 0,
                (label, "nonzero reduced pressure difference", index))
        require(isinstance(row.get("unreduced_difference_monomials"), int) and
                row.get("unreduced_difference_monomials") >= 0,
                (label, "pressure witness unreduced shape", index))
        require(row.get("ordinary_derivative_product_rule") is True and
                row.get("positive_original_L_shared") is True,
                (label, "pressure witness source rules", index))
        identities = row.get("actual_input_identity_sha256")
        require(isinstance(identities, Mapping) and identities,
                (label, "pressure witness input identities", index))
        require(all(isinstance(name, str) and isinstance(value, str) and len(value) == 64
                    for name, value in identities.items()),
                (label, "pressure witness input digest", index))
    require(all(sorted(seen[key]) == expected_orders[key] for key in expected_orders),
            (label, "pressure witness order coverage", seen))
    leading = witness.get('independent_leading_rows')
    require(isinstance(leading,list) and len(leading)==2 and
            {row.get('ordinary_Z_order') for row in leading}=={0,1} and
            all(row.get('passed') is True and row.get('quantity')=='leading_Qy' and
                row.get('reduced_difference_monomials')==0 and row.get('common_axially_constant_Pstar') is True
                for row in leading), (label,'independent leading Qy algebra'))
    formula_hash = witness.get("formula_AST_sha256")
    require(isinstance(formula_hash, str) and len(formula_hash) == 64,
            (label, "pressure witness formula AST hash"))
    require(witness.get("source_values_selected") is False and
            witness.get("numeric_caps_used") is False,
            (label, "pressure witness value provenance"))
    return {
        "passed": True,
        "rows": len(rows),
        "quantities": {key: len(value) for key, value in seen.items()},
        "orders": {key: sorted(value) for key, value in seen.items()},
        "all_reduced_difference_monomials_zero": True,
        "formula_AST_sha256": formula_hash,
        "source_values_selected": False,
        "numeric_caps_used": False,
    }


def _check_actual_function_identities(source: Mapping[str, Any]) -> dict[str, Any]:
    """Require explicit source-function identity evidence at both joins."""
    seam_map = source.get("seams")
    require(isinstance(seam_map, Mapping), "source derivative seam map missing")
    output = {}
    for name in ("second_power", "power_long"):
        seam = seam_map.get(name)
        require(isinstance(seam, Mapping), (name, "source derivative seam missing"))
        identity = seam.get("actual_function_identity")
        require(isinstance(identity, Mapping) and identity.get("passed") is True,
                (name, "actual function identity evidence missing"))
        for key in ("source_values_selected", "numeric_caps_used", "ancestor_replayed"):
            if key in identity:
                require(identity.get(key) is False, (name, key, "identity must remain source-only"))
        require(source.get('actual_function_congruence_completed') is True and
                identity.get('native_derivative_congruence_consumed_by_downstream_adapter') is True and
                identity.get('native_derivative_congruence_unproved_here') is False,
                (name, 'actual function congruence is incomplete'))
        loop = _mapping(source.get('actual_loop_function_identity'))
        derivative = _mapping(source.get('native_derivative_congruence'))
        require(loop.get('passed') is True and derivative.get('passed') is True,
                (name, 'actual loop and native derivative identities required'))
        loop_seam = _mapping(_mapping(loop.get('seams')).get(name))
        require(loop_seam.get('passed') is True and _mapping(loop_seam.get('scoped_fingerprint')).get('passed') is True,
                (name, 'scoped source authorizations failed'))
        for key, quantities in (('root_rows', ('E','V','a','b','p1','p2','t0','Delta','kappa')),
                                ('loop_rows', ('q','B_over_Pstar','B_over_Pstar_phi'))):
            rows = loop_seam.get(key)
            require(isinstance(rows,list) and len(rows)==3*len(quantities), (name,key,'row count'))
            require({(r.get('quantity'),r.get('ordinary_Z_order')) for r in rows} ==
                    {(q,k) for q in quantities for k in (0,1,2)}, (name,key,'coverage'))
            require(all(r.get('passed') is True and r.get('left_endpoint_sha256')==r.get('right_endpoint_sha256')
                        for r in rows), (name,key,'actual fingerprints differ'))
        require(_mapping(loop_seam.get('physical_y_derivative_congruence')).get('passed') is True,
                (name,'physical source derivative congruence failed'))
        def value_digest(value):
            return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()
        require(identity.get('derivative_certificate_sha256')==value_digest(derivative), (name,'foreign native derivative certificate'))
        fields=identity.get('rows')
        require(isinstance(fields,list) and len(fields)==9, (name,'actual field identity rows'))
        for quantity,key,count in (('U','modified_U_rows',3),('Uy','modified_Uy_rows',2),('m','modified_m_rows',4)):
            for k in range(count):
                row=next((v for v in fields if v.get('quantity')==quantity and v.get('ordinary_Z_order')==k),None)
                original=seam[key][k]
                require(isinstance(row,Mapping) and row.get('passed') is True and
                        row.get('left_node')==original.get('left_node') and row.get('right_node')==original.get('right_node') and
                        row==original.get('actual_function_identity'), (name,quantity,k,'actual input pair'))
                if quantity=='Uy':
                    require(row.get('issued_derivative_wiring_sha256')==value_digest(derivative['charts']), (name,quantity,k,'foreign derivative wiring'))
                if row.get('exact_endpoint_fingerprints') is True:
                    require(row.get('left_endpoint_sha256')==row.get('right_endpoint_sha256'), (name,quantity,k,'direct identity'))
                else:
                    require(row.get('actual_loop_operator_identity_proved') is True and
                            row.get('actual_loop_identity_sha256')==value_digest(loop), (name,quantity,k,'operator identity'))
        output[name] = {
            "passed": True,
            "source_values_selected": identity.get("source_values_selected", False),
            "numeric_caps_used": identity.get("numeric_caps_used", False),
        }
    return {"passed": True, "joins": output}


def _check_fixture_shape(fixture: Any) -> dict[str, Any]:
    require(isinstance(fixture, Mapping) and fixture.get("passed") is True,
            "independent upstream fixture failed")
    expected = {
        "synthetic_local_germs_only": True,
        "original_source_values_evaluated": False,
        "original_Poisson_inverse_evaluated": False,
        "NS_residual": False,
        "finite_N": 17,
        "interface_count": 2,
        "physical_R_once": True,
        "ordinary_Z1_nonzero": True,
        "nonzero_incoming_mean": True,
        "leading_subtraction_nonzero": True,
    }
    for key, value in expected.items():
        require(fixture.get(key) == value, ("fixture flag", key, fixture.get(key)))
    controls = fixture.get("controls")
    expected_controls = {
        "shared_finite_N_once", "phase_reset_rejected", "second_positive_J",
        "post_positive_varying_J", "long_positive_varying_J",
        "second_physical_y_derivative_required", "post_physical_y_derivative_required",
        "long_physical_y_derivative_required", "nonzero_incoming_mean",
        "leading_subtraction_retained", "physical_R_required", "ordinary_Z1_nonzero",
        "post_endpoint_scope_retained", "long_endpoint_scope_retained",
        "flat_width_cutoff_derivatives", "finite_N_error_atoms_are_not_values",
    }
    require(isinstance(controls, Mapping) and set(controls) == expected_controls and
            all(value is True for value in controls.values()), ("fixture controls", controls))
    require(fixture.get("control_count") == len(expected_controls), "fixture control count changed")
    require(mp.mpf(str(fixture.get("maximum_supporting_difference"))) < mp.mpf("1e-70"),
            "fixture tolerance changed")
    chain = fixture.get("signed_chain_C0_Z1_differences")
    require(isinstance(chain, list) and len(chain) == 2 and
            all(mp.mpf(str(value)) < mp.mpf("1e-70") for value in chain),
            "fixture chain differences changed")
    finite_rows = fixture.get("finite_N_error_rows")
    require(isinstance(finite_rows, list) and len(finite_rows) == 2 and
            all(row.get("positive_error_bound_atom_only") is True and
                row.get("source_value_selected") is False and
                row.get("ordinary_Z_orders") == [0, 1] for row in finite_rows if isinstance(row, Mapping)),
            "fixture finite-N rows changed")
    return {"passed": True, "finite_N": fixture.get("finite_N"),
            "interface_count": fixture.get("interface_count"),
            "control_count": fixture.get("control_count"),
            "maximum_supporting_difference": fixture.get("maximum_supporting_difference")}


def _check_source_certificates(report: Mapping[str, Any]) -> dict[str, Any]:
    source = report.get("source_derivative_theorem")
    native = report.get("native_source_guards")
    seams = report.get("source_seam_proofs")
    for label, cert in (("source derivative", source), ("native", native), ("seams", seams)):
        require(isinstance(cert, Mapping) and cert.get("passed") is True,
                f"{label} certificate failed")
        for key in ("source_values_selected", "numeric_caps_used", "ancestor_replayed"):
            if key in cert:
                require(cert.get(key) is False, (label, key, "must remain false"))
    seam_map = seams.get("joins") if isinstance(seams.get("joins"), Mapping) else seams.get("seams")
    require(isinstance(seam_map, Mapping), "upstream source seam map missing")
    selected = {}
    pressure_witnesses = {}
    for name in ("second_power", "power_long"):
        row = seam_map.get(name)
        require(isinstance(row, Mapping) and row.get("passed") is True,
                (name, "required source seam missing"))
        selected[name] = row
    # The pressure trace is the existing 15-row physical trace.  The upstream
    # adapter additionally exposes separate modified/leading source rows;
    # accepting only a generic flag here would permit a cap-only alias.
    for name, row in selected.items():
        trace = row.get("trace_rows")
        require(isinstance(trace, list) and len(trace) == 15 and
                all(isinstance(item, Mapping) and item.get("passed") is True for item in trace),
                (name, "pressure trace"))
        modified = _rows(row, ("modified_source_rows", "modified_trace_rows", "actual_modified_rows"))
        leading = _rows(row, ("leading_source_rows", "leading_trace_rows", "actual_leading_rows"))
        _check_source_rows(modified, {"U": 3, "Uy": 2, "m": 4}, (name, "modified"))
        _check_source_rows(leading, {"V": 3, "m": 3, "Vy": 2}, (name, "leading"))
        d_rows = row.get("D_equality_rows")
        require(isinstance(d_rows, list) and len(d_rows) == 2 and
                sorted(item.get("ordinary_Z_order") for item in d_rows if isinstance(item, Mapping)) == [0, 1] and
                all(isinstance(item, Mapping) and item.get("passed") is True for item in d_rows),
                (name, "D ordinary-Z rows"))
        require(row.get("endpoint_wrapper_binding_passed") is True,
                (name, "endpoint wrapper binding"))
        pressure_witnesses[name] = _check_pressure_algebra_witness(row, name)
    function_identities = _check_actual_function_identities(source)
    require(_find_flag(native, {"physical_R_once", "physical_radius_once", "R_factor_once"}),
            "physical R source premise missing")
    require(_find_flag(native, {"ordinary_Z1", "ordinary_Z_orders_0_1", "ordinary_Z1_available"}) or
            _mapping(native.get("scope")).get("ordinary_Z_orders") == [0, 1],
            "ordinary Z1 source premise missing")
    require(_find_flag(native, {"second_native_J_positive", "post_power_native_J_positive",
                                "long_native_J_positive", "variable_J_not_replaced_by_constant",
                                "native_J_is_varying", "nonconstant_J_source_dependency"}),
            "native Jacobian source premise missing")
    return {"passed": True, "interfaces": ["second_power", "power_long"],
            "native_source_only": True, "trace_rows_each": 15,
            "pressure_algebra_witness": pressure_witnesses,
            "actual_function_identity": function_identities}


def _check_recovery(report: Mapping[str, Any]) -> dict[str, Any]:
    seams = _mapping(report.get("source_seam_proofs"))
    seam_map = seams.get("joins") if isinstance(seams.get("joins"), Mapping) else seams.get("seams")
    rows = {}
    for name in ("second_power", "power_long"):
        row = _mapping(seam_map.get(name))
        trace = row.get("trace_rows")
        quantities = [_quantity(item) for item in trace if isinstance(item, Mapping)]
        require(quantities.count("U") == 3 and quantities.count("Uy") == 2 and
                quantities.count("m") == 3 and quantities.count("Q") == 3 and
                quantities.count("Qy") == 2 and quantities.count("RQy") == 2,
                (name, "pressure trace quantities"))
        rows[name] = {"trace_rows": len(trace), "pressure_quantities": quantities}
    return {"passed": True, "interfaces": rows, "strict_per_seam": True}


def _check_budget(sample: Mapping[str, Any]) -> dict[str, Any]:
    budget = sample.get("external_endpoint_budget_witness")
    require(isinstance(budget, Mapping), "upstream external budget missing")
    own = None
    for key in ("own_second_switch0", "own_second0", "own_second_switch_native0"):
        if isinstance(budget.get(key), Mapping):
            own = budget[key]
            break
    own_rows = own.get("radius_weighted_C0_C1_log_bounds") if own is not None else None
    mode = str(budget.get("mode", budget.get("selected_mode", "")))
    fallback = None
    for key in ("full_second_switch_68_fallback", "second_switch_68_fallback",
                "full_raw68_fallback", "second_switch_upper_aggregate"):
        if isinstance(budget.get(key), Mapping):
            fallback = budget[key]
            break
    if own is not None:
        require(isinstance(own_rows, list) and len(own_rows) == 2,
                "own second-switch C0/Z1 rows changed")
        for row in own_rows:
            _check_record(row, "own second-switch")
        require(mode in ("", "own_second_switch0", "own_second0") or
                "second" in mode.lower(), "unexpected upstream budget mode")
        active = own_rows
        active_kind = "own_second_switch0"
    else:
        require(isinstance(fallback, Mapping), "neither own second-switch nor raw68 fallback exists")
        atom = fallback.get("aggregate_log_atom", fallback.get("log_atom", fallback.get("atom")))
        _check_record(atom, "full raw68 fallback")
        active = [atom]
        active_kind = "full_raw68_fallback"
    inputs = budget.get("positive_union_inputs")
    require(isinstance(inputs, list) and inputs == active,
            "active upstream budget inputs changed")
    union = _check_positive_union(inputs, budget.get("positive_union"), "upstream positive union")
    require(budget.get("no_cap_subtraction_value_midpoint_or_overlap") is True,
            "upstream budget subtraction contract changed")
    require(_find_flag(budget, {"physical_R_once", "physical_radius_once", "ordinary_Z1"}),
            "upstream budget physical source scope missing")
    return {"passed": True, "mode": active_kind, "components": union["components"],
            "cap_as_value": False}


def _check_samples(report: Mapping[str, Any], parent: Mapping[str, Any]) -> list[dict[str, Any]]:
    samples = report.get("samples")
    parent_samples = {row.get("Z"): row for row in parent.get("samples", []) if isinstance(row, Mapping)}
    require(isinstance(samples, list) and tuple(row.get("Z") for row in samples) == SAVED_Z,
            "upstream saved Z set changed")
    output = []
    new_names = set()
    for sample in samples:
        require(isinstance(sample, Mapping) and sample.get("Z") in parent_samples,
                ("sample", sample.get("Z")))
        base = parent_samples[sample["Z"]]
        raw = sample.get("raw_m15_C0_C1_terms")
        previous = base.get("correlated_reduced_m15_C0_C1_terms")
        reduced = sample.get("correlated_reduced_m15_C0_C1_terms")
        require(isinstance(raw, list) and len(raw) == RAW_COUNT and raw == base.get("raw_m15_C0_C1_terms"),
                (sample.get("Z"), "raw 213 changed"))
        require(raw[SECOND_SWITCH_RAW_INDEX].get("name") == REMOVED_PARENT_TERMS[0] and
                raw[POST_POWER_RAW_INDEX].get("name") == REMOVED_PARENT_TERMS[1],
                (sample.get("Z"), "upstream raw finite-N identities"))
        require(isinstance(previous, list) and len(previous) == PARENT_REDUCED_COUNT and
                isinstance(reduced, list) and len(reduced) == NEW_REDUCED_COUNT,
                (sample.get("Z"), "upstream term counts changed"))
        removed = [row for row in previous if row.get("name") in REMOVED_PARENT_TERMS]
        require(len(removed) == 3 and {row.get("name") for row in removed} == set(REMOVED_PARENT_TERMS),
                (sample.get("Z"), "removed upstream terms"))
        keep = [row for row in previous if row.get("name") not in set(REMOVED_PARENT_TERMS)]
        require(len(keep) == UNCHANGED_COUNT, (sample.get("Z"), "parent unchanged count"))
        keep_names = {row.get("name") for row in keep}
        replacement = [row for row in reduced if row.get("name") not in keep_names]
        require(len(replacement) == 1, (sample.get("Z"), "single upstream replacement"))
        replacement = replacement[0]
        replacement_name = replacement.get("name")
        require(isinstance(replacement_name, str) and "second" in replacement_name.lower() and
                "external" in replacement_name.lower() and "viscous" in replacement_name.lower(),
                (sample.get("Z"), "upstream replacement identity", replacement_name))
        new_names.add(replacement_name)
        new_keep = [row for row in reduced if row.get("name") in keep_names]
        require(new_keep == keep and len(new_keep) == UNCHANGED_COUNT,
                (sample.get("Z"), "other 186 sectors changed"))
        binding = _mapping(replacement.get("actual_function_binding"))
        require(binding and _find_flag(binding, {"source_functions_not_cap_values", "error_enclosure_only",
                                                   "source_bound", "signed_source_function"}) or
                isinstance(replacement.get("log_scale_terms"), list),
                (sample.get("Z"), "replacement source binding"))
        require(sum(row.get("name") == "reserved_n1_own_swirl_once" for row in reduced) == 1,
                (sample.get("Z"), "reserved term count"))
        require(sample.get("independent_target") == base.get("independent_target"),
                (sample.get("Z"), "independent target changed"))
        contract = sample.get("contract")
        require(isinstance(contract, Mapping), (sample.get("Z"), "sample contract missing"))
        required_contract = {
            "raw_213_verbatim", "previous_189_to_new_187", "other_186_verbatim",
            "independent_target_unchanged", "reserved_retained_once",
            "ordinary_Z1_and_physical_R_once", "leading_subtraction_retained",
            "second_power_source_join_proved", "power_long_source_join_proved",
            "no_error_cap_as_value", "all_other_volume_width_numerical_errors_retained",
            "source_functions_not_cap_values", "ancestors_not_reexecuted",
        }
        for key in required_contract:
            require(contract.get(key) is True, (sample.get("Z"), "sample contract", key))
        require(contract.get("old_compact_external_removed") is True and
                contract.get("new_second_switch_external_once") is True,
                (sample.get("Z"), "sample external replacement contract"))
        limits = sample.get("limits")
        if isinstance(limits, Mapping):
            require(all(value is False or ("internal" in str(key).lower() and "trace" in str(key).lower())
                        for key, value in limits.items()), (sample.get("Z"), "sample limits"))
        output.append({"Z": sample.get("Z"), "raw_213": True, "reduced_187": True,
                       "unchanged_186": True, "budget": _check_budget(sample),
                       "replacement": replacement_name, "target_unchanged": True})
    require(len(new_names) == 1, "upstream replacement name changed across Z")
    return output


def _check_packet(report: Mapping[str, Any], parent: Mapping[str, Any],
                  fixture_summary: Mapping[str, Any] | None = None) -> dict[str, Any]:
    require(report.get("source_family") == parent.get("source_family"), "upstream source family changed")
    guard_controls = report.get('pressure_input_identity_guard_controls')
    require(isinstance(guard_controls, Mapping) and guard_controls.get('passed') is True and
            guard_controls.get('control_count') == 3 and
            guard_controls.get('controls') == {'missing_identity': True, 'foreign_node_pair': True, 'unequal_fingerprints': True} and
            guard_controls.get('source_values_evaluated') is False, 'pressure input identity negative controls failed')
    source_check = _check_source_certificates(report)
    recovery_check = _check_recovery(report)
    samples = _check_samples(report, parent)
    binding = _mapping(report.get("numerical_source_binding"))
    require(binding, "upstream numeric source binding missing")
    require(binding.get("physical_radius_power", 1) in (1, True) and
            binding.get("physical_R_once", True) in (1, True),
            "upstream physical R scope changed")
    require(binding.get("ordinary_Z_orders", [0, 1]) == [0, 1],
            "upstream ordinary-Z scope changed")
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
    return {"source": source_check, "recovery": recovery_check, "samples": samples,
            "global_limits_all_false": True, "cap_as_value": False}


def _load_parent() -> Mapping[str, Any]:
    path = DIRECTORY / PARENT
    require(path.is_file() and sha(path) == PARENT_SHA, "accepted tail parent changed")
    parent = read_json(path)
    require(isinstance(parent, Mapping) and _mapping(parent.get("checks")).get("passed") is True,
            "accepted tail parent is not passed")
    return parent


def run() -> dict[str, Any]:
    started = time.perf_counter()
    from check_original_upstream_m15 import check_fixture
    fixture = check_fixture()
    fixture_summary = _check_fixture_shape(fixture)
    from openai_ns_reconstruction.original_upstream_m15 import OriginalUpstreamM15
    owner = OriginalUpstreamM15(precision=150)
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
        "parent_reduced_189_to_new_187": True,
        "unchanged_other_sectors_186": True,
        "two_actual_source_interfaces": True,
        "second_power_and_power_long_source_proved": True,
        "pressure_algebra_witness_checked": True,
        "actual_function_identity_checked": True,
        "source_functions_not_caps": True,
        "production_is_new_owner_and_fixture_once": True,
    }
    report["execution_seconds"] = round(time.perf_counter() - started, 3)
    return report


def inspect(path: Path) -> dict[str, Any]:
    data = read_json(path)
    require(isinstance(data, Mapping) and _mapping(data.get("checks")).get("passed") is True,
            "saved upstream report is not passed")
    check_hash_group(data.get("input_hashes"), "input", DIRECTORY)
    check_hash_group(data.get("implementation_hashes"), "implementation", ROOT)
    parent = _load_parent()
    fixture_summary = _check_fixture_shape(_mapping(data.get("checks")).get("independent_fixture"))
    packet_check = _check_packet(data, parent, fixture_summary)
    checks = _mapping(data.get("checks"))
    for key, value in {
        "saved_Z": list(SAVED_Z),
        "raw_213_verbatim": True,
        "parent_reduced_189_to_new_187": True,
        "unchanged_other_sectors_186": True,
        "two_actual_source_interfaces": True,
        "second_power_and_power_long_source_proved": True,
        "pressure_algebra_witness_checked": True,
        "actual_function_identity_checked": True,
        "source_functions_not_caps": True,
    }.items():
        require(checks.get(key) == value, ("saved check", key))
    return {"checks": checks, "fixture_invariants": fixture_summary,
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
