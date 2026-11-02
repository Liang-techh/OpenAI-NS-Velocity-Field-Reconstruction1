"""Freeze-checked, one-entry attachment of the saved W0-chi correction.

Only the two authenticated receipt blobs named below are consumed. Their raw
bytes are hashed before either gzip blob is decoded; the synthetic fixture is
run first. The prior segmented result is used only through its saved
``formal_sum`` rows, while the correction contributes the already aggregated
``W0_chi`` swirl and its separate error stream. No old producer, cell
reduction, physical rescaling, endpoint term, or ledger write is performed.
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

from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "n1_w0_chi_attach_v1"
PREVIOUS_COMPONENT = "n1_segmented_inner_v2"
DIR = ROOT / "tmp" / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entered.json", "result.json"))

PRIOR = "experiments/root_st073/n1_segmented_inner_v2_receipt.json.gz"
CORRECTION = "experiments/root_st073/n1_core_kchi_feedback_v2_receipt.json.gz"
CORE_W0 = "experiments/root_st073/n1_core_w0_m15_v1_receipt.json.gz"
LEAF = "src/openai_ns_reconstruction/n1_w0_chi_attach.py"
FIXTURE = "src/openai_ns_reconstruction/n1_w0_chi_attach_fixture_v1.py"
RUNNER = "scripts/run_n1_w0_chi_attach_v1_once.py"
REQUIRED = (LEAF, FIXTURE, RUNNER, PRIOR, CORRECTION)
RECEIPTS = frozenset((PRIOR, CORRECTION))

RAW_SHA256 = {
    PRIOR: "c5b719e18ecabeeaed82bc496d5885178bd4bd7f0e0b29e9430b016ea8287c29",
    CORRECTION: "0ff1992f08754f2e15c8ec8e65f3ed396b38dc3101821c35eb90635664a753aa",
}
RESULT_SHA256 = {
    PRIOR: "03118b6f84254f902b92367dcfa5d6f3917f48aae70766727f2651daf75bb893",
    CORRECTION: "30926458ca6ec25a05ae6cde00d846adb15f9fe50361e46ec5e80deb7c6abae1",
}
POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
ORDERS = (0, 1)
SIGNED_COUNTS = {0: 1, 1: 3}
SIGNED_NET_COUNTS = {0: 1, 1: 2}
ERROR_COUNTS = {0: 1, 1: 3}


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


def decode(c: MPIntervalContext, value: Any, label: str = "value") -> Any:
    """Hydrate all endpoint dictionaries in one caller-owned context."""

    _reject_float_tree(value, label)
    if isinstance(value, Mapping):
        has_lower = "lower_exact_mpf_tuple" in value
        has_upper = "upper_exact_mpf_tuple" in value
        if has_lower or has_upper:
            require(has_lower and has_upper, label + " has incomplete exact endpoints")
            if "directed" in value:
                require(value["directed"] is True, label + " explicitly disables directed endpoints")
            lower = value["lower_exact_mpf_tuple"]
            upper = value["upper_exact_mpf_tuple"]
            require(
                isinstance(lower, (list, tuple))
                and isinstance(upper, (list, tuple))
                and len(lower) == len(upper) == 4
                and all(type(item) is int for item in list(lower) + list(upper)),
                label + " endpoint tuple",
            )
            result = c.make_mpf((tuple(lower), tuple(upper)))
            require(result.a <= result.b, label + " reversed interval")
            return result
        return {key: decode(c, item, label + "." + str(key)) for key, item in value.items()}
    if isinstance(value, list):
        return [decode(c, item, label + "[]") for item in value]
    if isinstance(value, tuple):
        return [decode(c, item, label + "[]") for item in value]
    if isinstance(value, float):
        raise TypeError(label + " contains a host floating-point value")
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
    require(len(value) == 2 or all(key in value for key in (0, 1, "0", "1")), label + " must expose exactly Z0 and Z1")
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


def preflight() -> tuple[dict[str, str], str, str, str]:
    """Validate the frozen pending brief and hash every byte before entry."""

    if not __debug__:
        raise RuntimeError("optimized mode disables one-entry guards")
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")

    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    require(
        isinstance(brief, Mapping)
        and brief.get("component") == COMPONENT
        and brief.get("phase") == COMPONENT + "_pending"
        and brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("one_entry_only") is True
        and brief.get("do_not_reenter") is True,
        "brief pending guard",
    )
    checkpoint = brief.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "source checkpoint missing")
    require(Path(str(checkpoint.get("path", ""))).resolve() == SOURCE.resolve(), "source checkpoint path")
    require(checkpoint.get("sha256") == sha(SOURCE.read_bytes()), "source checkpoint SHA")

    preimage_path = Path(str(checkpoint.get("preimage_path", "")))
    preimage_raw = preimage_path.read_bytes()
    require(sha(preimage_raw) == checkpoint.get("preimage_sha256"), "brief preimage SHA")
    preimage = json.loads(preimage_raw)
    require(
        isinstance(preimage, Mapping)
        and preimage.get("component") == PREVIOUS_COMPONENT
        and preimage.get("phase") == "component_published",
        "previous published preimage guard",
    )

    source_raw = SOURCE.read_bytes()
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
        if path in RECEIPTS:
            require(row.get("raw_sha256_before_decode") is True, "receipt lacks raw preflight marker " + path)
            require(digest == RAW_SHA256[path], "changed pinned receipt " + path)
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), "incomplete source manifest")
    return hashes, sha(source_raw), sha(brief_raw), sha(preimage_raw)


def _term_ids(term: Mapping[str, Any]) -> list[str]:
    values: list[Any] = []
    origin = term.get("origin")
    if isinstance(origin, Mapping):
        for key in ("source_ids", "source_id", "source_paths", "source_path"):
            value = origin.get(key)
            if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
                values.extend(value)
            elif value is not None:
                values.append(value)
    for key in ("source_ids", "source_id", "source_paths", "source_path"):
        value = term.get(key)
        if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
            values.extend(value)
        elif value is not None:
            values.append(value)
    return [str(value) for value in values if isinstance(value, (str, int))]


def _validate_powers(term: Mapping[str, Any], label: str) -> tuple[int, int, int, int]:
    powers = term.get("powers")
    require(
        isinstance(powers, Sequence)
        and not isinstance(powers, (str, bytes, bytearray))
        and len(powers) == 4
        and all(type(power) is int for power in powers),
        label + " formal powers",
    )
    return tuple(powers)  # type: ignore[return-value]


def _validate_input_term(term: Any, label: str) -> dict[str, Any]:
    require(isinstance(term, Mapping), label + " term mapping")
    require(hasattr(term.get("coefficient"), "_mpi_"), label + " coefficient was not hydrated")
    _validate_powers(term, label)
    require(_term_ids(term), label + " source path/id missing")
    return dict(term)


def _branch_order(branch: Any, order: int, label: str) -> Mapping[str, Any]:
    rows = _orders(branch, label)
    value = rows[order]
    require(isinstance(value, Mapping), label + " Z" + str(order) + " row")
    return value


def _marker(term: Mapping[str, Any], *, error_only: bool, label: str) -> None:
    origin = term.get("origin")
    require(isinstance(origin, Mapping), label + " origin")
    require(term.get("signed_field") is (not error_only), label + " signed_field marker")
    require(origin.get("signed_field") is (not error_only), label + " origin signed_field marker")
    require(bool(term.get("uncertainty", origin.get("uncertainty", False))), label + " uncertainty marker")
    if error_only:
        require(term.get("error_only") is True or origin.get("error_only") is True, label + " error_only marker")
        require(
            term.get("symmetric_error") is True
            or origin.get("symmetric_error") is True
            or origin.get("bound_type") in ("symmetric_error", "tail_bound", "kernel_tail"),
            label + " symmetric error marker",
        )
    else:
        require(term.get("error_only", False) is False and origin.get("error_only", False) is False, label + " signed/error separation")


def _validate_correction(
    c: MPIntervalContext, output: Mapping[str, Any]
) -> tuple[dict[int, dict[str, list[dict[str, Any]]]], dict[str, Any]]:
    m15 = output.get("m15")
    require(isinstance(m15, Mapping), "Kchi correction m15 mapping")
    correction_basis = m15.get("source_power_basis", output.get("source_power_basis"))
    require(tuple(correction_basis or ()) == POWER_BASIS, "Kchi correction source power basis mismatch")
    branches = m15.get("branches")
    require(isinstance(branches, Mapping), "Kchi correction branches mapping")
    branch = branches.get("W0_chi")
    require(isinstance(branch, Mapping), "Kchi correction W0_chi branch missing")

    packet: dict[int, dict[str, list[dict[str, Any]]]] = {}
    index: dict[str, Any] = {}
    for order in ORDERS:
        row = _branch_order(branch, order, "Kchi W0_chi")
        for field in ("axial", "axial_error", "swirl", "swirl_error"):
            raw = row.get(field, [])
            require(isinstance(raw, list), "Kchi W0_chi " + field + " must be a list")
        require(not row["axial"] and not row["axial_error"], "Kchi W0_chi axial channels must be empty")
        signed = [dict(decode(c, term, "correction.W0_chi.Z" + str(order) + ".swirl")) for term in row["swirl"]]
        errors = [dict(decode(c, term, "correction.W0_chi.Z" + str(order) + ".swirl_error")) for term in row["swirl_error"]]
        require(len(signed) == SIGNED_COUNTS[order], "Kchi W0_chi signed count Z" + str(order))
        require(len(errors) == ERROR_COUNTS[order], "Kchi W0_chi error count Z" + str(order))
        for channel, terms, is_error in (("swirl", signed, False), ("swirl_error", errors, True)):
            for number, term in enumerate(terms):
                label = "correction.W0_chi.Z" + str(order) + "." + channel + "[" + str(number) + "]"
                _validate_input_term(term, label)
                _marker(term, error_only=is_error, label=label)
                generated = "W0_chi:" + str(order) + ":" + channel + ":" + str(number)
                origin = dict(term.get("origin", {}))
                parents = _term_ids(term)
                origin["source_ids"] = [generated] + [item for item in origin.get("source_ids", ()) if str(item) != generated]
                origin["parent_source_ids"] = parents
                origin["receipt"] = CORRECTION
                origin["path"] = ["result", "output", "m15", "branches", "W0_chi", str(order), channel, number]
                term["origin"] = origin
                term["attach_source_id"] = generated
                term["receipt"] = CORRECTION
                term["path"] = list(origin["path"])
                index[generated] = {
                    "receipt": CORRECTION,
                    "result_sha256": RESULT_SHA256[CORRECTION],
                    "path": list(origin["path"]),
                    "parent_source_ids": parents,
                    "error_only": is_error,
                    "powers": list(term["powers"]),
                }
        packet[order] = {"axial": [], "swirl": signed, "axial_error": [], "swirl_error": errors}
        require(
            sorted(tuple(term["powers"]) for term in signed) == sorted(tuple(term["powers"]) for term in errors),
            "Kchi W0_chi signed/error powers differ at Z" + str(order),
        )
    return packet, {
        "generated_source_index": index,
        "m15_path": ["result", "output", "m15"],
        "receipt": CORRECTION,
        "result_sha256": RESULT_SHA256[CORRECTION],
    }


def _flags_false(value: Mapping[str, Any], keys: Sequence[str], label: str) -> None:
    for key in keys:
        if key in value:
            require(value[key] is False, label + " flag " + key)


def _channel_rows(result: Mapping[str, Any], channel: str) -> dict[int, list[Any]]:
    """Read a leaf's explicit signed/error channel without grouping by powers."""

    candidates: list[Any] = []
    for key in (
        channel,
        "signed" if channel == "signed" else "error_only",
        "signed_correction" if channel == "signed" else "error_correction",
    ):
        if key in result:
            candidates.append(result[key])
    for container_key in ("correction", "attached", "W0_chi", "m15"):
        container = result.get(container_key)
        if not isinstance(container, Mapping):
            continue
        candidates.append(container.get(channel))
        branches = container.get("branches")
        if isinstance(branches, Mapping):
            branch = branches.get("W0_chi")
            if isinstance(branch, Mapping):
                candidates.append(branch.get(channel))
                rows: dict[int, Any] = {}
                for order in ORDERS:
                    try:
                        rows[order] = _branch_order(branch, order, "leaf W0_chi").get(channel, [])
                    except (KeyError, TypeError, ValueError):
                        rows[order] = []
                candidates.append(rows)
    for candidate in candidates:
        if isinstance(candidate, Mapping):
            try:
                return {order: list(_orders(candidate, "leaf " + channel)[order]) for order in ORDERS}
            except (KeyError, TypeError, ValueError):
                continue
    return {order: [] for order in ORDERS}


def _all_term_ids(rows: Mapping[int, Sequence[Any]]) -> list[str]:
    result: list[str] = []
    for order in ORDERS:
        for term in rows.get(order, ()):
            if isinstance(term, Mapping):
                result.extend(_term_ids(term))
    return result


def _validate_leaf_result(
    actual: Mapping[str, Any],
    prior_formal: Mapping[int, Any],
    correction: Mapping[int, Mapping[str, list[dict[str, Any]]]],
    correction_index: Mapping[str, Any],
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
        ),
        "leaf",
    )
    require(actual.get("owner_identity_proved") is False, "leaf owner identity must remain unproved")
    if "source_family" in actual:
        require(actual["source_family"] == source_family, "leaf source family changed")
    if "source_logs" in actual:
        require(actual["source_logs"] == source_logs, "leaf source logs changed")

    raw_signed = actual.get("signed_correction")
    require(isinstance(raw_signed, Mapping), "leaf raw signed_correction missing")
    raw_signed_rows = {
        order: list(_orders(raw_signed, "leaf signed_correction")[order])
        for order in ORDERS
    }
    require(
        len(raw_signed_rows[0]) == SIGNED_COUNTS[0]
        and len(raw_signed_rows[1]) == SIGNED_COUNTS[1],
        "leaf raw signed W0_chi counts",
    )
    for order in ORDERS:
        require(
            all(
                any(term is original for original in correction[order]["swirl"])
                for term in raw_signed_rows[order]
            ),
            "leaf raw signed rows changed or duplicated",
        )
    signed_net_raw = actual.get("signed_correction_net")
    require(isinstance(signed_net_raw, Mapping), "leaf signed_correction_net missing")
    signed = {
        order: list(_orders(signed_net_raw, "leaf signed_correction_net")[order])
        for order in ORDERS
    }
    require(
        len(signed[0]) == SIGNED_NET_COUNTS[0]
        and len(signed[1]) == SIGNED_NET_COUNTS[1],
        "leaf signed W0_chi net counts",
    )
    for order in ORDERS:
        require(
            all(
                raw_signed_rows[order][index] is correction[order]["swirl"][index]
                for index in range(len(raw_signed_rows[order]))
            ),
            "leaf raw signed W0_chi order/identity changed at Z" + str(order),
        )
    signed_alias = actual.get("signed")
    if isinstance(signed_alias, Mapping):
        signed_alias_rows = _orders(signed_alias, "leaf signed alias")
        for order in ORDERS:
            require(
                len(signed_alias_rows[order]) == SIGNED_NET_COUNTS[order]
                and all(
                    signed_alias_rows[order][index] is signed[order][index]
                    for index in range(SIGNED_NET_COUNTS[order])
                ),
                "leaf signed alias is not the powers-net stream at Z" + str(order),
            )
    errors = _channel_rows(actual, "error_only")
    require(len(errors[0]) == ERROR_COUNTS[0] and len(errors[1]) == ERROR_COUNTS[1], "leaf error-only W0_chi counts")
    for order in ORDERS:
        require(
            all(
                errors[order][index] is correction[order]["swirl_error"][index]
                for index in range(len(errors[order]))
            ),
            "leaf error-only W0_chi order/identity changed at Z" + str(order),
        )
    signed_ids = set(_all_term_ids(signed))
    error_ids = set(_all_term_ids(errors))
    generated_signed = {
        source_id
        for source_id, record in correction_index.items()
        if isinstance(record, Mapping) and record.get("error_only") is False
    }
    generated_error = {
        source_id
        for source_id, record in correction_index.items()
        if isinstance(record, Mapping) and record.get("error_only") is True
    }
    require(generated_signed.issubset(signed_ids), "leaf dropped W0_chi signed source IDs")
    require(generated_error.issubset(error_ids), "leaf dropped W0_chi error source IDs")
    require(not signed_ids.intersection(generated_error), "leaf inserted error-only source into signed stream")

    consumed_ids: list[int] = []
    for order in ORDERS:
        for number, net_term in enumerate(signed[order]):
            input_terms = net_term.get("input_terms") if isinstance(net_term, Mapping) else None
            require(isinstance(input_terms, Sequence) and not isinstance(input_terms, (str, bytes, bytearray)), "leaf net input_terms")
            for raw_term in input_terms:
                require(any(raw_term is candidate for candidate in raw_signed_rows[order]), "leaf net input is not a raw signed term")
                require(not bool(raw_term.get("error_only", False)) if isinstance(raw_term, Mapping) else False, "leaf net consumed an error term")
                consumed_ids.append(id(raw_term))
            require(input_terms, "leaf net group has no raw signed input")
    require(
        sorted(consumed_ids) == sorted(id(term) for order in ORDERS for term in raw_signed_rows[order]),
        "leaf net did not consume each raw signed term exactly once",
    )

    for order in ORDERS:
        for number, term in enumerate(signed[order]):
            require(isinstance(term, Mapping), "leaf signed term mapping")
            _validate_powers(term, "leaf signed term")
            _marker(term, error_only=False, label="leaf signed Z" + str(order) + "[" + str(number) + "]")
        for number, term in enumerate(errors[order]):
            require(isinstance(term, Mapping), "leaf error term mapping")
            _validate_powers(term, "leaf error term")
            _marker(term, error_only=True, label="leaf error Z" + str(order) + "[" + str(number) + "]")
        require(
            {tuple(term["powers"]) for term in signed[order]}
            == {tuple(term["powers"]) for term in errors[order]},
            "leaf mixed signed/error powers at Z" + str(order),
        )

    prior_output = actual.get("prior_formal_sum")
    formal_output = actual.get("formal_sum")
    require(isinstance(prior_output, Mapping), "leaf prior_formal_sum missing")
    require(isinstance(formal_output, Mapping), "leaf formal_sum missing")
    prior_rows_out = _orders(prior_output, "leaf prior_formal_sum")
    formal_rows_out = _orders(formal_output, "leaf formal_sum")
    net_rows_out = _orders(signed_net_raw, "leaf signed_correction_net")
    for order in ORDERS:
        expected_prior = list(prior_formal[order])
        actual_prior = list(prior_rows_out[order])
        actual_formal = list(formal_rows_out[order])
        actual_net = list(net_rows_out[order])
        require(len(actual_prior) == len(expected_prior), "leaf prior row count changed at Z" + str(order))
        require(
            all(actual_prior[index] is expected for index, expected in enumerate(expected_prior)),
            "leaf prior row identity changed at Z" + str(order),
        )
        require(len(actual_formal) == len(expected_prior) + len(actual_net), "leaf formal prefix/suffix count at Z" + str(order))
        require(
            all(actual_formal[index] is expected for index, expected in enumerate(expected_prior)),
            "leaf formal prior prefix changed at Z" + str(order),
        )
        require(
            all(actual_formal[len(expected_prior) + index] is term for index, term in enumerate(actual_net)),
            "leaf formal signed suffix changed at Z" + str(order),
        )
    require(actual.get("cells_resummed", False) is False, "leaf re-summed grouped cells")
    require(actual.get("reaggregated_cells", False) is False, "leaf re-aggregated grouped cells")
    return {
        "signed_counts": {str(order): len(signed[order]) for order in ORDERS},
        "error_only_counts": {str(order): len(errors[order]) for order in ORDERS},
        "signed_source_id_count": len(signed_ids),
        "error_only_source_id_count": len(error_ids),
        "prior_source_id_count": sum(len(prior_formal[order]) for order in ORDERS),
        "prior_prefix_retained_once": True,
        "error_only_not_in_signed": True,
        "powers_kept_separate": True,
        "physical_scale_reapplied": False,
        "endpoint_cancellation_proved": False,
        "ledger_212_changed": False,
        "owner_identity_proved": False,
    }


def _compact_leaf_result(actual: Mapping[str, Any]) -> dict[str, Any]:
    """Keep the result structure without serializing parent aliases repeatedly."""

    compact: dict[str, Any] = {}
    omitted_parent_aliases = {
        "parent_input",
        "parent_formal_sum_input",
        "prior_formal_sum",
        "prior",
        "parent",
        "prior_terms",
        "signed_correction_input",
        "error_only_input",
        "correction_input",
        "selected_correction_input",
        "w0_chi_input",
        "prior_source_index",
        "correction_source_index",
        "source_index",
        "parent_error_stream",
        "segmented_candidate",
        "candidate",
        "prior_plus_correction",
        "incremental",
        "formal_sum",
        "formal_expression",
        "incremental_formal_expression",
        "formal_expression_by_order",
    }
    for key, value in actual.items():
        if key in omitted_parent_aliases:
            continue
        if key in {"prior_output", "prior_receipt", "prior_result"}:
            compact[key + "_reference"] = {
                "receipt": PRIOR,
                "result_sha256": RESULT_SHA256[PRIOR],
                "path": [
                    "result",
                    "output",
                    "result",
                    "formal_sum" if key != "prior_source_index" else "source_index",
                ],
            }
            continue
        compact[str(key)] = encode(value)
    compact["prior_formal_sum_reference"] = {
        "receipt": PRIOR,
        "result_sha256": RESULT_SHA256[PRIOR],
        "path": ["result", "output", "result", "formal_sum"],
    }
    compact["formal_sum_reference"] = {
        "prior_prefix_receipt": PRIOR,
        "prior_prefix_result_sha256": RESULT_SHA256[PRIOR],
        "prior_prefix_path": ["result", "output", "result", "formal_sum"],
        "signed_suffix_path": ["result_structure", "signed_correction_net"],
        "expanded_formal_sum_serialized": False,
    }
    compact["prior_source_index_reference"] = {
        "receipt": PRIOR,
        "result_sha256": RESULT_SHA256[PRIOR],
        "path": ["result", "output", "result", "source_index"],
    }
    return compact


def consume(hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    # Authenticate every frozen byte again after the exclusive entry was saved.
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "changed after entry: " + path)
        raw[path] = data
    for path in RECEIPTS:
        require(sha(raw[path]) == RAW_SHA256[path], "changed pinned raw receipt: " + path)
    progress["raw_bytes_authenticated"] = True

    leaf = module(LEAF, COMPONENT + "_leaf")
    fixture = module(FIXTURE, COMPONENT + "_fixture")
    context = MPIntervalContext()
    context.prec = 1800
    fixture_runner = getattr(fixture, "run_fixture", getattr(fixture, "run", None))
    require(callable(fixture_runner), "fixture must export run_fixture")
    progress["fixture_started"] = True
    try:
        proof = fixture_runner(leaf, context)
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

    # No saved receipt is decoded until the fixture has passed.
    progress["input_json_decode_started"] = True
    prior_receipt = _receipt(raw, PRIOR, PREVIOUS_COMPONENT, RESULT_SHA256[PRIOR])
    progress["prior_receipt_decoded"] = True
    correction_receipt = _receipt(raw, CORRECTION, "n1_core_kchi_feedback_v2", RESULT_SHA256[CORRECTION])
    progress["correction_receipt_decoded"] = True

    prior_output = prior_receipt["result"].get("output")
    correction_output = correction_receipt["result"].get("output")
    require(isinstance(prior_output, Mapping), "prior output mapping")
    require(isinstance(correction_output, Mapping), "correction output mapping")
    prior_result = prior_output.get("result")
    require(isinstance(prior_result, Mapping), "prior nested result mapping")
    prior_formal = decode(context, prior_result.get("formal_sum"), "prior.result.formal_sum")
    prior_rows = _orders(prior_formal, "prior.result.formal_sum")
    for order in ORDERS:
        require(isinstance(prior_rows[order], list), "prior formal Z" + str(order) + " row")
        prior_rows[order] = [
            _validate_input_term(term, "prior.result.formal_sum.Z" + str(order) + "[" + str(i) + "]")
            for i, term in enumerate(prior_rows[order])
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
    prior_owner = prior_result.get("source_owner_identity_proved", prior_result.get("owner_identity_proved"))
    require(prior_owner is False, "prior owner identity must remain unproved")

    prior_family_raw = prior_output.get("source_family", prior_result.get("source_family"))
    prior_logs_raw = prior_output.get("source_logs", prior_result.get("source_logs"))
    correction_family_raw = correction_output.get("source_family")
    correction_logs_raw = correction_output.get("source_logs")
    _reject_float_tree(prior_family_raw, "prior.source_family")
    _reject_float_tree(prior_logs_raw, "prior.source_logs")
    _reject_float_tree(correction_family_raw, "correction.source_family")
    _reject_float_tree(correction_logs_raw, "correction.source_logs")
    source_family = decode(context, prior_family_raw, "prior.source_family")
    source_logs = decode(context, prior_logs_raw, "prior.source_logs")
    require(source_family == decode(context, correction_family_raw, "correction.source_family"), "source family mismatch")
    require(source_logs == decode(context, correction_logs_raw, "correction.source_logs"), "source logs mismatch")
    source_power_basis = prior_result.get("source_power_basis", prior_output.get("source_power_basis"))
    require(tuple(source_power_basis or ()) == POWER_BASIS, "source power basis mismatch")

    prior_manifest = _manifest_hashes(prior_receipt, "prior segmented receipt")
    correction_manifest = _manifest_hashes(correction_receipt, "Kchi correction receipt")
    require(CORE_W0 in prior_manifest and CORE_W0 in correction_manifest, "core W0 source-manifest pin missing")
    require(prior_manifest[CORE_W0] == correction_manifest[CORE_W0], "cross-receipt core W0 source-manifest pin mismatch")
    prior_output_pins = prior_output.get("input_receipt_sha256", {})
    if isinstance(prior_output_pins, Mapping) and CORE_W0 in prior_output_pins:
        require(prior_output_pins[CORE_W0] == prior_manifest[CORE_W0], "prior core W0 output pin mismatch")

    correction_packet, correction_meta = _validate_correction(context, correction_output)
    progress["correction_shape_validated"] = True
    prior_source_index = decode(context, prior_result.get("source_index", {}), "prior.result.source_index")
    correction_source_index = decode(context, correction_output.get("input_index", {}), "correction.output.input_index")
    source_index = {
        "prior": prior_source_index,
        "W0_chi": correction_meta["generated_source_index"],
    }
    owner_certificate = {
        "source_family_exact": True,
        "source_logs_exact": True,
        "owner_identity_proved": False,
        "cross_receipt_core_W0_manifest_pin": True,
    }
    scale_certificate = {
        "prior_records_already_physically_scaled": True,
        "correction_records_already_physically_scaled": True,
        "physical_scale_reapplied": False,
        "endpoint_added": False,
        "ledger_written": False,
    }

    progress["actual_leaf_build_entered"] = True
    try:
        actual = leaf.build(
            context,
            prior_formal_sum=prior_rows,
            correction={
                "branches": {"W0_chi": correction_packet},
                "source_index": correction_meta["generated_source_index"],
                "receipt": CORRECTION,
                "result_sha256": RESULT_SHA256[CORRECTION],
                "aggregated_eight_cells": True,
                "cells_resummed": False,
            },
            w0_chi=correction_packet,
            source_family=source_family,
            source_logs=source_logs,
            source_atoms=source_logs,
            source_index=source_index,
            prior_source_index=prior_source_index,
            correction_source_index=correction_source_index,
            owner_certificate=owner_certificate,
            scale_certificate=scale_certificate,
            certificate=owner_certificate,
        )
    finally:
        progress["actual_leaf_build_calls"] += 1
    require(isinstance(actual, Mapping), "W0 chi attach leaf result mapping")
    progress["leaf_result_validated"] = True
    validation = _validate_leaf_result(
        actual,
        prior_rows,
        correction_packet,
        correction_meta["generated_source_index"],
        source_family,
        source_logs,
    )

    correction_compact = {
        str(order): {
            "swirl": encode(correction_packet[order]["swirl"]),
            "swirl_error": encode(correction_packet[order]["swirl_error"]),
            "axial": [],
            "axial_error": [],
        }
        for order in ORDERS
    }
    return {
        "schema": COMPONENT,
        "source_family": encode(source_family),
        "source_logs": encode(source_logs),
        "source_power_basis": list(POWER_BASIS),
        "ordinary_Z_orders": [0, 1],
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
        "prior_source_index_reference": {
            "receipt": PRIOR,
            "path": ["result", "output", "result", "source_index"],
        },
        "correction_source_index_reference": {
            "receipt": CORRECTION,
            "result_sha256": RESULT_SHA256[CORRECTION],
            "path": ["result", "output", "input_index"],
        },
        "parent_ledger_formula_reference": {
            "receipt": PRIOR,
            "result_sha256": RESULT_SHA256[PRIOR],
            "formal_sum_path": ["result", "output", "result", "formal_sum"],
            "source_power_basis_path": ["result", "output", "result", "source_power_basis"],
            "ledger_212_changed_path": ["result", "output", "result", "ledger_212_changed"],
            "ledger_212_changed": False,
        },
        "W0_chi_source_index": encode(correction_meta["generated_source_index"]),
        "cross_receipt_core_W0_pin": {
            "path": CORE_W0,
            "prior_manifest_sha256": prior_manifest[CORE_W0],
            "correction_manifest_sha256": correction_manifest[CORE_W0],
            "equal": True,
        },
        "correction": correction_compact,
        "result_structure": _compact_leaf_result(actual),
        "validation": validation,
        "prior_errors_preserved": True,
        "prior_ledger_212_preserved": True,
        "ledger_212_changed": False,
        "ledger_written": False,
        "endpoint_cancellation_proved": False,
        "physical_scale_reapplied": False,
        "same_interval_context_validated": False,
        "same_interval_context_is_caller_obligation": True,
        "old_prefix_repeated": False,
        "grouped_eight_cells_resummed": False,
        "source_owner_identity_proved": False,
        "tail_source_reference": {
            "receipt": CORRECTION,
            "path": ["result", "output", "m15", "branches", "W0_chi", "swirl_error"],
            "candidate_error_not_summed": True,
            "not_a_signed_field": True,
            "error_tail_is_independent_reference": False,
            "error_tail_kept_as_separate_reference": True,
            "error_tail_independence_not_proved": True,
            "upstream_receipts": [
                "experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz",
                "experiments/root_st073/n1_kchi_inlet_tail_bind_v1_receipt.json.gz",
            ],
        },
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
            "fixture_before_saved_receipt_decode": True,
        },
    )
    progress: dict[str, Any] = {
        "raw_bytes_authenticated": False,
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




