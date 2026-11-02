"""One-entry, no-addition consumer for the saved W0-chi records.

This v2 leaf is intentionally different from ``n1_w0_chi_attach_v1``.  The
parent formal rows and the raw W0-chi signed records are retained as separate
formal-expression segments.  The four error-only records are retained as a
separate reference stream.  No powers grouping, interval addition, source
reconstruction, scale application, endpoint claim, or ledger write is done.

The runner authenticates all frozen bytes before decoding either successful
receipt.  The archived v1 failure is authenticated as raw gzip bytes only and
is never decoded.  The synthetic fixture runs before any saved receipt JSON is
decoded and uses coefficients that raise on arithmetic.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "n1_w0_chi_ungrouped_v3"
PREVIOUS_COMPONENT = "n1_w0_chi_ungrouped_v2"
DIR = ROOT / "tmp" / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entered.json", "result.json"))

PRIOR = "experiments/root_st073/n1_segmented_inner_v2_receipt.json.gz"
CORRECTION = "experiments/root_st073/n1_core_kchi_feedback_v2_receipt.json.gz"
FAILURE = "experiments/root_st073/n1_w0_chi_attach_v1_failure.json.gz"
FAILURE_V2 = "experiments/root_st073/n1_w0_chi_ungrouped_v2_failure.json.gz"
CORE_W0 = "experiments/root_st073/n1_core_w0_m15_v1_receipt.json.gz"
LEAF = "src/openai_ns_reconstruction/n1_w0_chi_ungrouped_v2.py"
FIXTURE = "src/openai_ns_reconstruction/n1_w0_chi_ungrouped_fixture_v3.py"
RUNNER = "scripts/run_n1_w0_chi_ungrouped_v3_once.py"

REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR, CORRECTION, FAILURE, FAILURE_V2)
RECEIPTS = frozenset((PRIOR, CORRECTION))
RAW_ONLY = frozenset((FAILURE, FAILURE_V2))
PINNED_RECEIPTS = RECEIPTS | RAW_ONLY

RAW_SHA256 = {
    PRIOR: "c5b719e18ecabeeaed82bc496d5885178bd4bd7f0e0b29e9430b016ea8287c29",
    CORRECTION: "0ff1992f08754f2e15c8ec8e65f3ed396b38dc3101821c35eb90635664a753aa",
    FAILURE: "c42ca1a83bb4af4ca3b6122651c907da1f8e5b0f719522afd0539f33478ebedf",
    FAILURE_V2: "61892de3b51ee4b16e252e1823b0c7fe1a2e98aee0b2766be207535909e17549",
}
RESULT_SHA256 = {
    PRIOR: "03118b6f84254f902b92367dcfa5d6f3917f48aae70766727f2651daf75bb893",
    CORRECTION: "30926458ca6ec25a05ae6cde00d846adb15f9fe50361e46ec5e80deb7c6abae1",
}
POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDERS = (0, 1)
SIGNED_COUNTS = {0: 1, 1: 3}
ERROR_COUNTS = {0: 1, 1: 3}
PRIOR_COUNTS = {0: 615, 1: 863}


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def sha_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(value: object, message: str) -> None:
    if not value:
        raise ValueError(message)


def save(path: Path, value: object) -> str:
    raw = (
        json.dumps(
            value,
            sort_keys=True,
            ensure_ascii=False,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    with path.open("xb") as handle:
        handle.write(raw)
    return sha(raw)


def _reject_float_tree(value: Any, label: str) -> None:
    if isinstance(value, float):
        raise TypeError(label + " contains a host floating-point value")
    if isinstance(value, Mapping):
        for key, item in value.items():
            _reject_float_tree(key, label + ".key")
            _reject_float_tree(item, label)
    elif isinstance(value, (list, tuple, set, frozenset)):
        for item in value:
            _reject_float_tree(item, label)


def encode(value: Any) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return {
            "lower_exact_mpf_tuple": list(raw[0]),
            "upper_exact_mpf_tuple": list(raw[1]),
            "directed": True,
        }
    if isinstance(value, Mapping):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return [encode(item) for item in sorted(value, key=repr)]
    if isinstance(value, float):
        raise TypeError("cannot serialize a host floating-point value")
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    raise TypeError("cannot serialize result type " + type(value).__name__)


def decode(value: Any, label: str = "value") -> Any:
    """Hydrate exact interval endpoint dictionaries without midpointing."""

    _reject_float_tree(value, label)
    if isinstance(value, Mapping):
        has_lower = "lower_exact_mpf_tuple" in value
        has_upper = "upper_exact_mpf_tuple" in value
        if has_lower or has_upper:
            require(has_lower and has_upper, label + " has incomplete exact endpoints")
            lower = value["lower_exact_mpf_tuple"]
            upper = value["upper_exact_mpf_tuple"]
            require(
                isinstance(lower, (list, tuple))
                and isinstance(upper, (list, tuple))
                and len(lower) == len(upper) == 4
                and all(type(item) is int for item in list(lower) + list(upper)),
                label + " endpoint tuple",
            )
            # Keep receipt intervals opaque until the caller-owned mpmath
            # context is available.  The actual hydration occurs in
            # ``hydrate`` below, never in the synthetic fixture.
            return {"__exact_interval__": (tuple(lower), tuple(upper))}
        return {key: decode(item, label + "." + str(key)) for key, item in value.items()}
    if isinstance(value, list):
        return [decode(item, label + "[]") for item in value]
    if isinstance(value, tuple):
        return [decode(item, label + "[]") for item in value]
    return value


def hydrate(context: Any, value: Any, label: str = "value") -> Any:
    if isinstance(value, Mapping):
        if set(value) == {"__exact_interval__"}:
            lower, upper = value["__exact_interval__"]
            result = context.make_mpf((lower, upper))
            require(result.a <= result.b, label + " reversed interval")
            return result
        return {key: hydrate(context, item, label + "." + str(key)) for key, item in value.items()}
    if isinstance(value, list):
        return [hydrate(context, item, label + "[]") for item in value]
    return value


def _lookup(mapping: Mapping[Any, Any], key: int | str, label: str) -> Any:
    if key in mapping:
        return mapping[key]
    alternate: int | str = str(key) if isinstance(key, int) else int(key)
    if alternate in mapping:
        return mapping[alternate]
    raise KeyError(label + " missing " + str(key))


def _orders(value: Any, label: str) -> dict[int, Any]:
    require(isinstance(value, Mapping), label + " must be an ordinary-Z mapping")
    result = {order: _lookup(value, order, label) for order in ORDERS}
    require(
        len(value) == 2 or all(key in value for key in (0, 1, "0", "1")),
        label + " must expose exactly Z0 and Z1",
    )
    return result


def _manifest_hashes(receipt: Mapping[str, Any], label: str) -> dict[str, str]:
    source = receipt.get("source")
    require(isinstance(source, Mapping), label + " source manifest missing")
    rows = source.get("files")
    require(isinstance(rows, list), label + " source manifest files missing")
    result: dict[str, str] = {}
    for row in rows:
        require(isinstance(row, Mapping), label + " source manifest row")
        path, digest = row.get("path"), row.get("sha256")
        require(isinstance(path, str) and isinstance(digest, str), label + " source manifest path")
        require(path not in result, label + " duplicate source manifest path " + path)
        result[path] = digest
    return result


def _canonical_result_sha(result: Mapping[str, Any]) -> str:
    raw = (
        json.dumps(result, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")
    return sha(raw)


def _receipt(raw: Mapping[str, bytes], path: str, component: str, expected_sha: str) -> dict[str, Any]:
    value = json.loads(gzip.decompress(raw[path]))
    require(isinstance(value, Mapping), path + " receipt envelope")
    result = value.get("result")
    require(isinstance(result, Mapping), path + " result missing")
    require(result.get("status") == "returned" and result.get("component") == component, path + " component/status lock")
    require(_canonical_result_sha(result) == expected_sha, path + " saved result SHA")
    return dict(value)


def module(path: str, name: str) -> Any:
    require(path in (LEAF, FIXTURE), "unexpected import path: " + path)
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    require(spec is not None and spec.loader is not None, "module unavailable: " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def _previous_pending_guard(value: Any) -> bool:
    if not isinstance(value, Mapping):
        return False
    component = value.get("component")
    phase = value.get("phase")
    return component in (PREVIOUS_COMPONENT, PREVIOUS_COMPONENT + "_pending") and phase in (
        "component_pending",
        PREVIOUS_COMPONENT + "_pending",
        "component_published",
    )


def preflight() -> tuple[dict[str, str], str, str, str]:
    """Freeze source and raw receipt bytes before any gzip decode."""

    if not __debug__:
        raise RuntimeError("optimized mode disables one-entry guards")
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")

    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    require(
        isinstance(brief, Mapping)
        and brief.get("component") == COMPONENT
        and brief.get("phase") in (COMPONENT + "_pending", "component_pending")
        and brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("one_entry_only") is True
        and brief.get("do_not_reenter") is True,
        "brief pending guard",
    )
    checkpoint = brief.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "source checkpoint missing")
    require(Path(str(checkpoint.get("path", ""))).resolve() == SOURCE.resolve(), "source checkpoint path")
    source_raw = SOURCE.read_bytes()
    require(checkpoint.get("sha256") == sha(source_raw), "source checkpoint SHA")

    preimage_value = checkpoint.get("preimage_path")
    preimage_sha = checkpoint.get("preimage_sha256")
    if preimage_value is None:
        preimage_value = brief.get("previous_pending_preimage")
        if isinstance(preimage_value, Mapping):
            preimage_sha = preimage_value.get("sha256")
            preimage_value = preimage_value.get("path")
    require(isinstance(preimage_value, str) and isinstance(preimage_sha, str), "pending preimage pointer")
    preimage_path = Path(preimage_value)
    preimage_raw = preimage_path.read_bytes()
    require(sha(preimage_raw) == preimage_sha, "pending preimage SHA")
    preimage = json.loads(preimage_raw)
    require(_previous_pending_guard(preimage), "previous attach pending guard")

    source = json.loads(source_raw)
    require(
        isinstance(source, Mapping)
        and source.get("component") == COMPONENT
        and source.get("one_entry_only") is True
        and source.get("previous_component") == PREVIOUS_COMPONENT,
        "source manifest guard",
    )
    policy = source.get("source_policy", {})
    require(isinstance(policy, Mapping), "source policy missing")
    for key in (
        "source_bind_after_freeze",
        "authenticate_all_raw_gzip_bytes_before_decode",
        "previous_receipts_decode_only_after_entry",
        "fixture_before_saved_receipt_consumption",
        "complete_result_and_failures_are_saved",
        "no_ancestor_import_or_reconstruction",
        "no_numeric_addition",
        "formal_expression_segments_are_separate",
        "failure_archive_raw_only",
    ):
        require(policy.get(key) is True, "source policy guard " + key)

    rows = source.get("files")
    require(isinstance(rows, list) and len(rows) == len(REQUIRED), "source manifest file count")
    hashes: dict[str, str] = {}
    for row in rows:
        require(isinstance(row, Mapping), "source manifest row")
        path = row.get("path")
        require(isinstance(path, str) and path in REQUIRED and path not in hashes, "unexpected source file")
        digest = sha_path(ROOT / path)
        require(digest == row.get("sha256"), "changed frozen source file " + path)
        if path in PINNED_RECEIPTS:
            require(row.get("raw_sha256_before_decode") is True, "receipt lacks raw preflight marker " + path)
            require(digest == RAW_SHA256[path], "changed pinned raw receipt " + path)
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), "incomplete source manifest")
    return hashes, sha(source_raw), sha(brief_raw), sha(preimage_raw)


def _term_ids(term: Mapping[str, Any]) -> list[str]:
    values: list[Any] = []
    for container in (term.get("origin"), term):
        if not isinstance(container, Mapping):
            continue
        for key in ("source_ids", "source_id", "source_paths", "source_path"):
            value = container.get(key)
            if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
                values.extend(value)
            elif value is not None:
                values.append(value)
    return [str(value) for value in values if isinstance(value, (str, int))]


def _term_shape(term: Any, label: str, *, stream: str) -> dict[str, Any]:
    require(isinstance(term, Mapping), label + " term mapping")
    require("powers" in term and isinstance(term["powers"], Sequence), label + " powers missing")
    require(len(term["powers"]) == 4 and all(type(power) is int for power in term["powers"]), label + " powers")
    require(hasattr(term.get("coefficient"), "_mpi_"), label + " coefficient is not a directed interval")
    origin = term.get("origin")
    require(isinstance(origin, Mapping), label + " origin")
    require(_term_ids(term), label + " lineage/source reference missing")
    if stream == "error":
        require(origin.get("error_only") is True, label + " origin.error_only")
        require(term.get("error_only", True) is True, label + " error_only")
        require(term.get("signed_field") is False and origin.get("signed_field") is False, label + " signed marker")
    elif stream == "signed":
        require(term.get("signed_field") is True and origin.get("signed_field") is True, label + " signed marker")
        require(term.get("error_only", False) is False and origin.get("error_only", False) is False, label + " error marker")
    elif stream != "parent":
        raise ValueError(label + " unknown stream " + stream)
    return dict(term)


def _branch_row(branch: Mapping[str, Any], order: int) -> Mapping[str, Any]:
    row = _lookup(branch, order, "W0_chi branch")
    require(isinstance(row, Mapping), "W0_chi Z" + str(order) + " row")
    return row


def _validate_correction(context: Any, output: Mapping[str, Any]) -> tuple[dict[int, dict[str, list[dict[str, Any]]]], dict[str, Any]]:
    m15 = output.get("m15")
    require(isinstance(m15, Mapping), "Kchi correction m15 mapping")
    basis = m15.get("source_power_basis", output.get("source_power_basis"))
    require(tuple(basis or ()) == POWER_BASIS, "Kchi source power basis mismatch")
    branches = m15.get("branches")
    require(isinstance(branches, Mapping), "Kchi correction branches mapping")
    branch = branches.get("W0_chi")
    require(isinstance(branch, Mapping), "Kchi W0_chi branch missing")
    packet: dict[int, dict[str, list[dict[str, Any]]]] = {}
    refs: dict[str, Any] = {}
    seen_paths: set[tuple[Any, ...]] = set()
    for order in ORDERS:
        row = _branch_row(branch, order)
        require(row.get("axial", []) == [] and row.get("axial_error", []) == [], "W0_chi axial channels must be empty")
        signed_raw = row.get("swirl")
        error_raw = row.get("swirl_error")
        require(isinstance(signed_raw, list) and isinstance(error_raw, list), "W0_chi channels must be lists")
        require(len(signed_raw) == SIGNED_COUNTS[order], "W0_chi signed count Z" + str(order))
        require(len(error_raw) == ERROR_COUNTS[order], "W0_chi error count Z" + str(order))
        signed: list[dict[str, Any]] = []
        errors: list[dict[str, Any]] = []
        for channel, raw_rows, is_error in (("swirl", signed_raw, False), ("swirl_error", error_raw, True)):
            for index, raw in enumerate(raw_rows):
                term = hydrate(context, decode(raw, "correction.W0_chi"), "correction.W0_chi")
                term = _term_shape(term, "correction.W0_chi.Z" + str(order) + "." + channel, stream="error" if is_error else "signed")
                path = ("result", "output", "m15", "branches", "W0_chi", str(order), channel, index)
                require(path not in seen_paths, "duplicate generated correction path")
                seen_paths.add(path)
                generated = "W0_chi:" + str(order) + ":" + channel + ":" + str(index)
                refs[generated] = {
                    "receipt": CORRECTION,
                    "result_sha256": RESULT_SHA256[CORRECTION],
                    "path": list(path),
                    "error_only": is_error,
                    "powers": list(term["powers"]),
                    # Full upstream lineage stays in the raw term.  The
                    # source index may use a different namespace; absence of a
                    # same-name upstream id is not a rejection condition.
                    "source_refs": _term_ids(term),
                }
                (errors if is_error else signed).append(term)
        packet[order] = {"axial": [], "swirl": signed, "axial_error": [], "swirl_error": errors}
    require(len(refs) == 8 and len(seen_paths) == 8, "W0_chi generated record references")
    return packet, {
        "generated_source_index": refs,
        "m15_path": ["result", "output", "m15"],
        "receipt": CORRECTION,
        "result_sha256": RESULT_SHA256[CORRECTION],
        "generated_paths": [list(path) for path in sorted(seen_paths, key=repr)],
    }


def _flags_false(value: Mapping[str, Any], keys: Sequence[str], label: str) -> None:
    for key in keys:
        if key in value:
            require(value[key] is False, label + " flag " + key)


def _channel_rows(value: Any, channel: str) -> dict[int, list[Any]]:
    if not isinstance(value, Mapping):
        return {0: [], 1: []}
    candidates: list[Any] = []
    for key in (
        channel,
        "signed" if channel == "signed" else "error_only",
        "signed_rows" if channel == "signed" else "error_rows",
        "errors" if channel != "signed" else "signed_records",
        "raw_signed" if channel == "signed" else "raw_errors",
    ):
        if key in value:
            candidates.append(value[key])
    for container_key in ("W0_chi", "error_stream", "error_tail", "formal_expression"):
        container = value.get(container_key)
        if isinstance(container, Mapping):
            candidates.append(container.get(channel))
            candidates.append(container.get("signed" if channel == "signed" else "error"))
            candidates.append(container.get("error_only" if channel != "signed" else "signed_rows"))
    for candidate in candidates:
        if isinstance(candidate, Mapping):
            try:
                rows = {order: list(_orders(candidate, "leaf " + channel)[order]) for order in ORDERS}
            except (KeyError, TypeError, ValueError):
                continue
            if all(isinstance(rows[order], list) for order in ORDERS):
                return rows
    return {0: [], 1: []}


def _expression_rows(actual: Mapping[str, Any], order: int) -> list[Any]:
    expression = actual.get("formal_expression")
    if expression is None:
        expression = actual.get("formal_expression_by_order")
    require(isinstance(expression, Mapping), "leaf formal_expression missing")
    selected = _lookup(expression, order, "leaf formal_expression")
    if isinstance(selected, Mapping):
        prior = selected.get("prior")
        signed = selected.get("signed")
        if isinstance(prior, Sequence) and not isinstance(prior, (str, bytes, bytearray)) and isinstance(signed, Sequence) and not isinstance(signed, (str, bytes, bytearray)):
            return list(prior) + list(signed)
        for key in ("terms", "children", "subterms", "items"):
            if isinstance(selected.get(key), list):
                return selected[key]
    if isinstance(selected, Sequence) and not isinstance(selected, (str, bytes, bytearray)):
        return list(selected)
    raise TypeError("leaf formal_expression Z" + str(order) + " is not a segment list")


def _entry_term(value: Any) -> Any:
    if isinstance(value, Mapping):
        for key in ("term", "raw_term", "input_term", "record", "child"):
            if key in value:
                return value[key]
    return value


def _validate_leaf_result(
    actual: Mapping[str, Any],
    prior_rows: Mapping[int, Sequence[Any]],
    correction: Mapping[int, Mapping[str, list[dict[str, Any]]]],
    source_family: Any,
    source_logs: Any,
) -> dict[str, Any]:
    _flags_false(
        actual,
        (
            "physical_scale_reapplied",
            "physical_factors_applied",
            "physical_weight_reapplied",
            "radial_weight_reapplied",
            "geometry_reapplied",
            "endpoint_cancellation_proved",
            "terminal_closure_proved",
            "ledger_212_changed",
            "ledger_appended",
            "ledger_written",
            "writes_ledger",
            "numeric_values_computed",
            "midpoint_selected",
            "point_value_selected",
            "coefficients_operated_on",
            "numeric_addition_performed",
            "grouped_by_powers",
            "power_net_computed",
            "formal_sum_recomputed",
            "error_stream_summed",
            "error_only_added_to_signed",
        ),
        "leaf",
    )
    require(actual.get("owner_identity_proved", False) is False, "leaf owner identity must remain unproved")
    require(actual.get("correction_attached", True) is True, "leaf correction marker")
    require(actual.get("prior_and_signed_segments_separate", True) is True, "leaf segments not separate")
    require(actual.get("input_terms_retained_by_reference", True) is True, "leaf did not retain input terms")
    require(actual.get("error_terms_in_signed_expression", False) is False, "leaf inserted error-only terms")
    require(actual.get("numeric_addition_performed", False) is False, "leaf performed numeric addition")
    require(actual.get("coefficients_operated_on", False) is False, "leaf operated on coefficients")
    require(actual.get("error_terms_are_independent_references", False) is False, "leaf claimed error independence")
    require(actual.get("error_terms_are_separate_references", True) is True, "leaf did not keep error references")
    if "source_family" in actual:
        require(actual["source_family"] == source_family, "leaf source family changed")
    if "source_logs" in actual:
        require(actual["source_logs"] == source_logs, "leaf source logs changed")

    signed = _channel_rows(actual, "signed")
    errors = _channel_rows(actual, "error_only")
    for order in ORDERS:
        require(len(signed[order]) == SIGNED_COUNTS[order], "leaf raw signed count Z" + str(order))
        require(len(errors[order]) == ERROR_COUNTS[order], "leaf raw error count Z" + str(order))
        for index, term in enumerate(signed[order]):
            require(term is correction[order]["swirl"][index], "leaf signed raw identity changed")
            _term_shape(term, "leaf signed Z" + str(order), stream="signed")
        for index, term in enumerate(errors[order]):
            require(term is correction[order]["swirl_error"][index], "leaf error raw identity changed")
            _term_shape(term, "leaf error Z" + str(order), stream="error")

        expression = [_entry_term(item) for item in _expression_rows(actual, order)]
        expected = list(prior_rows[order]) + list(correction[order]["swirl"])
        require(len(expression) == len(expected), "leaf segmented expression count Z" + str(order))
        require(all(expression[index] is expected[index] for index in range(len(expected))), "leaf parent/signed order changed")
        require(not any(item is error for item in expression for error in correction[order]["swirl_error"]), "leaf error entered expression")

    return {
        "signed_counts": {str(order): len(signed[order]) for order in ORDERS},
        "error_only_counts": {str(order): len(errors[order]) for order in ORDERS},
        "prior_formal_counts": {str(order): len(prior_rows[order]) for order in ORDERS},
        "parent_order_retained": True,
        "raw_signed_records_retained": True,
        "raw_error_records_retained": True,
        "formal_expression_segmented": True,
        "error_only_not_in_signed": True,
        "error_independence_not_claimed": True,
        "physical_scale_reapplied": False,
        "endpoint_cancellation_proved": False,
        "terminal_closure_proved": False,
        "ledger_212_changed": False,
        "owner_identity_proved": False,
    }


def _record_refs(packet: Mapping[int, Mapping[str, list[dict[str, Any]]]], channel: str) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = {}
    for order in ORDERS:
        rows = packet[order][channel]
        result[str(order)] = [
            {
                "receipt": CORRECTION,
                "result_sha256": RESULT_SHA256[CORRECTION],
                "path": ["result", "output", "m15", "branches", "W0_chi", str(order), channel, index],
                "record": encode(term),
            }
            for index, term in enumerate(rows)
        ]
    return result


def consume(hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    # Re-authenticate all frozen bytes after exclusive entry creation.  This
    # reads the v1 failure archive only as bytes; it is never decompressed.
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "changed after entry: " + path)
        raw[path] = data
    for path in PINNED_RECEIPTS:
        require(sha(raw[path]) == RAW_SHA256[path], "changed pinned raw receipt: " + path)
    progress["raw_bytes_authenticated"] = True
    progress["failure_archive_raw_authenticated"] = True
    progress["failure_archive_decoded"] = False

    leaf = module(LEAF, COMPONENT + "_leaf")
    fixture = module(FIXTURE, COMPONENT + "_fixture")
    fixture_runner = getattr(fixture, "run_fixture", getattr(fixture, "run", None))
    require(callable(fixture_runner), "fixture must export run_fixture")
    progress["fixture_started"] = True
    try:
        # The fixture owns its non-numeric probe context and runs before the
        # first successful receipt is decompressed or parsed.
        proof = fixture_runner(leaf, None)
        progress["fixture_leaf_calls"] = int(proof.get("leaf_calls", 0)) if isinstance(proof, Mapping) else 0
        progress["fixture"] = encode(proof)
        require(isinstance(proof, Mapping) and proof.get("passed") is True, "synthetic fixture failed")
        checks = proof.get("checks")
        if isinstance(checks, Mapping):
            require(all(value is True for value in checks.values()), "synthetic fixture has a failed check")
    except BaseException as exc:
        progress["fixture_error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        raise
    progress["fixture_completed_before_input_json_decode"] = True
    progress["input_json_decode_started"] = True

    # Only these two successful receipts are decoded.  FAILURE is deliberately
    # absent from both calls below.
    prior_receipt = _receipt(raw, PRIOR, "n1_segmented_inner_v2", RESULT_SHA256[PRIOR])
    progress["prior_receipt_decoded"] = True
    correction_receipt = _receipt(raw, CORRECTION, "n1_core_kchi_feedback_v2", RESULT_SHA256[CORRECTION])
    progress["correction_receipt_decoded"] = True

    context_module = importlib.util.find_spec("mpmath")
    require(context_module is not None, "mpmath is required only for saved interval hydration")
    from mpmath.ctx_iv import MPIntervalContext

    context = MPIntervalContext()
    context.prec = 1800
    prior_output = prior_receipt["result"].get("output")
    correction_output = correction_receipt["result"].get("output")
    require(isinstance(prior_output, Mapping), "prior output mapping")
    require(isinstance(correction_output, Mapping), "correction output mapping")
    prior_result = prior_output.get("result")
    require(isinstance(prior_result, Mapping), "prior nested result mapping")
    prior_formal_raw = decode(prior_result.get("formal_sum"), "prior.result.formal_sum")
    prior_formal = hydrate(context, prior_formal_raw, "prior.result.formal_sum")
    prior_rows = _orders(prior_formal, "prior.result.formal_sum")
    for order in ORDERS:
        require(isinstance(prior_rows[order], list), "prior formal Z" + str(order) + " row")
        require(len(prior_rows[order]) == PRIOR_COUNTS[order], "prior formal count changed at Z" + str(order))
        prior_rows[order] = [
            _term_shape(term, "prior.result.formal_sum.Z" + str(order), stream="parent")
            for term in prior_rows[order]
        ]
    _flags_false(
        prior_result,
        (
            "physical_scale_reapplied",
            "physical_factors_applied",
            "physical_weight_reapplied",
            "radial_weight_reapplied",
            "geometry_reapplied",
            "endpoint_cancellation_proved",
            "terminal_closure_proved",
            "ledger_212_changed",
            "ledger_appended",
            "ledger_written",
            "writes_ledger",
        ),
        "prior",
    )
    require(prior_result.get("owner_identity_proved", prior_result.get("source_owner_identity_proved", False)) is False, "prior owner identity")

    prior_family_raw = prior_output.get("source_family", prior_result.get("source_family"))
    prior_logs_raw = prior_output.get("source_logs", prior_result.get("source_logs"))
    correction_family_raw = correction_output.get("source_family")
    correction_logs_raw = correction_output.get("source_logs")
    source_family = hydrate(context, decode(prior_family_raw, "prior.source_family"), "prior.source_family")
    source_logs = hydrate(context, decode(prior_logs_raw, "prior.source_logs"), "prior.source_logs")
    require(source_family == hydrate(context, decode(correction_family_raw, "correction.source_family"), "correction.source_family"), "source family mismatch")
    require(source_logs == hydrate(context, decode(correction_logs_raw, "correction.source_logs"), "correction.source_logs"), "source logs mismatch")
    source_power_basis = prior_result.get("source_power_basis", prior_output.get("source_power_basis"))
    require(tuple(source_power_basis or ()) == POWER_BASIS, "source power basis mismatch")

    prior_manifest = _manifest_hashes(prior_receipt, "prior segmented receipt")
    correction_manifest = _manifest_hashes(correction_receipt, "Kchi correction receipt")
    require(CORE_W0 in prior_manifest and CORE_W0 in correction_manifest, "common core W0 source-manifest pin missing")
    require(prior_manifest[CORE_W0] == correction_manifest[CORE_W0], "common core W0 source-manifest hash mismatch")
    correction_packet, correction_meta = _validate_correction(context, correction_output)
    progress["correction_shape_validated"] = True
    prior_source_index = hydrate(context, decode(prior_result.get("source_index", {}), "prior.result.source_index"), "prior.result.source_index")
    correction_source_index = hydrate(context, decode(correction_output.get("input_index", {}), "correction.output.input_index"), "correction.output.input_index")
    require(isinstance(prior_source_index, Mapping), "prior source index mapping")
    require(isinstance(correction_source_index, Mapping), "correction source index mapping")

    progress["actual_leaf_build_entered"] = True
    try:
        # Canonical v2 API: exactly three positional packets.  This is the
        # sole actual leaf call and cannot fall back to an older leaf.
        actual = leaf.build(
            prior_rows,
            {order: correction_packet[order]["swirl"] for order in ORDERS},
            {order: correction_packet[order]["swirl_error"] for order in ORDERS},
        )
    finally:
        progress["actual_leaf_build_calls"] += 1
    require(isinstance(actual, Mapping), "ungrouped leaf result mapping")
    progress["leaf_result_validated"] = True
    validation = _validate_leaf_result(actual, prior_rows, correction_packet, source_family, source_logs)

    return {
        "schema": COMPONENT,
        "source_power_basis": list(POWER_BASIS),
        "ordinary_Z_orders": [0, 1],
        "source_family": encode(source_family),
        "source_logs": encode(source_logs),
        "prior_receipt_reference": {
            "path": PRIOR,
            "raw_sha256": RAW_SHA256[PRIOR],
            "result_sha256": RESULT_SHA256[PRIOR],
            "output_path": ["result", "output"],
            "formal_sum_path": ["result", "output", "result", "formal_sum"],
        },
        "correction_receipt_reference": {
            "path": CORRECTION,
            "raw_sha256": RAW_SHA256[CORRECTION],
            "result_sha256": RESULT_SHA256[CORRECTION],
            "W0_chi_path": ["result", "output", "m15", "branches", "W0_chi"],
        },
        "failure_archive_reference": {
            "path": FAILURE,
            "raw_sha256": RAW_SHA256[FAILURE],
            "decoded": False,
            "used_as_success": False,
        },
        "parent_source_index_reference": {
            "receipt": PRIOR,
            "result_sha256": RESULT_SHA256[PRIOR],
            "path": ["result", "output", "result", "source_index"],
        },
        "correction_source_index_reference": {
            "receipt": CORRECTION,
            "result_sha256": RESULT_SHA256[CORRECTION],
            "path": ["result", "output", "input_index"],
            "complete_index_retained": True,
            "missing_same_name_upstream_id_is_not_rejection": True,
        },
        "common_core_W0_receipt_pin": {
            "path": CORE_W0,
            "prior_manifest_sha256": prior_manifest[CORE_W0],
            "correction_manifest_sha256": correction_manifest[CORE_W0],
            "equal": True,
        },
        "prior_formal_row_counts": {str(order): len(prior_rows[order]) for order in ORDERS},
        "raw_signed_records": _record_refs(correction_packet, "swirl"),
        "raw_error_records": _record_refs(correction_packet, "swirl_error"),
        "generated_source_paths": correction_meta["generated_paths"],
        "validation": validation,
        "formal_expression_compact": {
            "parent_rows_referenced": True,
            "parent_receipt": PRIOR,
            "signed_segment_counts": {str(order): SIGNED_COUNTS[order] for order in ORDERS},
            "error_segment_counts": {str(order): ERROR_COUNTS[order] for order in ORDERS},
            "expanded_parent_serialized": False,
            "power_grouping_performed": False,
            "coefficient_addition_performed": False,
        },
        "prior_errors_preserved": True,
        "prior_ledger_212_preserved": True,
        "ledger_212_changed": False,
        "ledger_written": False,
        "physical_scale_reapplied": False,
        "endpoint_cancellation_proved": False,
        "terminal_closure_proved": False,
        "source_owner_identity_proved": False,
        "formal_expression_segmented": True,
        "error_tail_reference_only": True,
        "error_tail_added_to_signed": False,
        "error_tail_is_independent_reference": False,
        "error_tail_independence_not_proved": True,
        "failure_archive_decoded": False,
        "actual_leaf_build_calls": progress["actual_leaf_build_calls"],
        "fixture_leaf_build_calls": progress["fixture_leaf_calls"],
    }


def main() -> int:
    hashes, source_sha256, brief_sha256, preimage_sha256 = preflight()
    DIR.mkdir(parents=True, exist_ok=True)
    entry_sha256 = save(
        ENTRY,
        {
            "component": COMPONENT,
            "previous_component": PREVIOUS_COMPONENT,
            "source_sha256": source_sha256,
            "brief_sha256": brief_sha256,
            "preimage_sha256": preimage_sha256,
            "input_hashes": hashes,
            "entered": True,
            "exclusive_entry_and_result": True,
            "raw_gzip_authenticated_before_decode": True,
            "failure_archive_raw_only": True,
            "fixture_before_saved_receipt_decode": True,
            "no_numeric_addition": True,
        },
    )
    progress: dict[str, Any] = {
        "raw_bytes_authenticated": False,
        "failure_archive_raw_authenticated": False,
        "failure_archive_decoded": False,
        "fixture_started": False,
        "fixture_completed_before_input_json_decode": False,
        "input_json_decode_started": False,
        "prior_receipt_decoded": False,
        "correction_receipt_decoded": False,
        "correction_shape_validated": False,
        "fixture_leaf_calls": 0,
        "actual_leaf_build_entered": False,
        "actual_leaf_build_calls": 0,
        "leaf_result_validated": False,
    }
    outcome: dict[str, Any] = {
        "component": COMPONENT,
        "previous_component": PREVIOUS_COMPONENT,
        "source_sha256": source_sha256,
        "entry_sha256": entry_sha256,
        "status": "failed",
        "raw_gzip_authenticated_before_decode": True,
        "failure_archive_raw_only": True,
        "fixture_before_saved_receipt_decode": True,
        "progress": progress,
        "output": None,
    }
    try:
        outcome["output"] = consume(hashes, progress)
        outcome["status"] = "returned"
        outcome["error"] = None
    except BaseException as exc:
        outcome["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    result_sha256 = save(RESULT, outcome)
    print(
        json.dumps(
            {
                "component": COMPONENT,
                "status": outcome["status"],
                "entry_sha256": entry_sha256,
                "result_sha256": result_sha256,
                "error": (outcome.get("error") or {}).get("message"),
            },
            sort_keys=True,
        )
    )
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
