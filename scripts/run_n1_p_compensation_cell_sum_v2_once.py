"""Run the fresh P-only physical m15 cell-sum consumer exactly once.

Only the two named successful gzip receipts are consumed.  Their compressed
bytes are authenticated before the exclusive entry receipt, the independent
fixture runs before either receipt is decoded, and the new leaf is called once
for the actual packet.  No producer, ancestor module, prior formal parent, or
continuation input is re-entered.
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
from typing import Any, Callable


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "n1_p_compensation_cell_sum_v2"
PREVIOUS_COMPONENT = "n1_p_compensation_cell_sum_v1"
DIR = ROOT / "tmp" / COMPONENT
SOURCE = DIR / "source.json"
ENTRY = DIR / "entered.json"
RESULT = DIR / "result.json"

P_RECEIPT = "experiments/root_st073/n1_core_kchi_compensation_v2_receipt.json.gz"
PRIOR_RECEIPT = "experiments/root_st073/n1_w0_chi_ungrouped_v3_receipt.json.gz"
FAILURE = "experiments/root_st073/n1_p_compensation_cell_sum_v1_failure.json.gz"
FEEDBACK_RECEIPT = "experiments/root_st073/n1_core_kchi_feedback_v2_receipt.json.gz"
LEAF = "src/openai_ns_reconstruction/n1_p_compensation_cell_sum_v2.py"
FIXTURE = "src/openai_ns_reconstruction/n1_p_compensation_cell_fixture_v2.py"
RUNNER = "scripts/run_n1_p_compensation_cell_sum_v2_once.py"

REQUIRED = (LEAF, FIXTURE, RUNNER, P_RECEIPT, PRIOR_RECEIPT, FAILURE)
GZIP_INPUTS = frozenset((P_RECEIPT, PRIOR_RECEIPT, FAILURE))
RAW_SHA256 = {
    P_RECEIPT: "c862f0560e8cb4b053abbbfcc982e4b4c1b602a806628e2d709399878e1c3aaf",
    PRIOR_RECEIPT: "3123a436957522fd3de3926dccacf6b463c28587a50951df6377337ffbab0821",
    FAILURE: "012e2d410a74ee16eb6eaf6859bac00884a29514c55a9330006ce09d54e421a1",
}
RESULT_SHA256 = {
    P_RECEIPT: "341a09693c0fc46682a6561c37df468ca536e62c59ee7dc944ad9841e16a53b5",
    PRIOR_RECEIPT: "d0f9f8d466828ea9cb02505f99d21500c6cab1966eeb113af22ea6941e5b4ee3",
}
FEEDBACK_RAW_SHA256 = "0ff1992f08754f2e15c8ec8e65f3ed396b38dc3101821c35eb90635664a753aa"
BRANCHES = ("G0_P_W0", "G0_P_DeltaW0", "D_P_W0", "D_P_DeltaW0")
ORDERS = (0, 1)
CELLS = tuple(range(8))
POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FORBIDDEN_NAMES = {
    "rm_loop_density_identity.json",
    "lei_ren_part1_paper_compliant_current_original_rm_loop_density_identity.json",
}


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


def _decode_exact(value: Any, label: str = "value") -> Any:
    """Keep receipt intervals as endpoint tokens until one context exists."""

    _reject_float_tree(value, label)
    if isinstance(value, Mapping):
        has_lower = "lower_exact_mpf_tuple" in value
        has_upper = "upper_exact_mpf_tuple" in value
        if has_lower or has_upper:
            require(has_lower and has_upper, label + " has incomplete exact interval")
            lower = value["lower_exact_mpf_tuple"]
            upper = value["upper_exact_mpf_tuple"]
            require(
                isinstance(lower, (list, tuple))
                and isinstance(upper, (list, tuple))
                and len(lower) == len(upper) == 4
                and all(type(item) is int for item in list(lower) + list(upper)),
                label + " exact endpoint tuple",
            )
            return {"__exact_interval__": (tuple(lower), tuple(upper))}
        return {
            key: _decode_exact(item, label + "." + str(key))
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_decode_exact(item, label + "[]") for item in value]
    return value


def _hydrate(context: Any, value: Any, label: str = "value") -> Any:
    if isinstance(value, Mapping):
        if set(value) == {"__exact_interval__"}:
            lower, upper = value["__exact_interval__"]
            result = context.make_mpf((lower, upper))
            require(result.a <= result.b, label + " reversed interval")
            return result
        return {
            key: _hydrate(context, item, label + "." + str(key))
            for key, item in value.items()
        }
    if isinstance(value, list):
        return [_hydrate(context, item, label + "[]") for item in value]
    return value


def _canonical_result_sha(result: Mapping[str, Any]) -> str:
    raw = (
        json.dumps(
            result,
            sort_keys=True,
            ensure_ascii=False,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    return sha(raw)


def _forbidden(path: Path) -> bool:
    name = path.name.lower()
    return name in FORBIDDEN_NAMES or "continuation" in name


def _source_rows(receipt: Mapping[str, Any], label: str) -> dict[str, str]:
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


def preflight() -> tuple[dict[str, str], str, str, str | None]:
    """Authenticate source, brief, preimage, and raw gzip bytes before entry."""

    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    source_raw = SOURCE.read_bytes()
    checkpoint = brief.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "source checkpoint missing")
    require(
        Path(str(checkpoint.get("path", ""))).resolve() == SOURCE.resolve()
        and checkpoint.get("sha256") == sha(source_raw),
        "source checkpoint",
    )
    require(
        brief.get("component") == COMPONENT
        and brief.get("phase") in (COMPONENT + "_pending", "component_pending")
        and brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("one_entry_only") is True
        and brief.get("do_not_reenter") is True,
        "brief pending guard",
    )

    preimage_path_value = checkpoint.get("preimage_path")
    preimage_sha = checkpoint.get("preimage_sha256")
    if preimage_path_value is None:
        previous = brief.get("previous_pending_preimage")
        if isinstance(previous, Mapping):
            preimage_path_value = previous.get("path")
            preimage_sha = previous.get("sha256")
    require(isinstance(preimage_path_value, str) and isinstance(preimage_sha, str), "pending preimage pointer")
    preimage_path = Path(preimage_path_value)
    require(not _forbidden(preimage_path), "forbidden preimage path")
    preimage_raw = preimage_path.read_bytes()
    preimage_sha256 = sha(preimage_raw)
    require(preimage_sha256 == preimage_sha, "pending preimage SHA")
    preimage = json.loads(preimage_raw)
    require(preimage.get("component") == PREVIOUS_COMPONENT and preimage.get("phase") == PREVIOUS_COMPONENT + "_pending", "published predecessor guard")

    manifest = json.loads(source_raw)
    require(
        isinstance(manifest, Mapping)
        and manifest.get("component") == COMPONENT
        and manifest.get("one_entry_only") is True
        and manifest.get("previous_component") == PREVIOUS_COMPONENT,
        "source manifest guard",
    )
    policy = manifest.get("source_policy", {})
    require(isinstance(policy, Mapping), "source policy mapping")
    for key in (
            "authenticate_all_raw_gzip_bytes_before_decode",
            "fixture_before_saved_receipt_consumption",
            "previous_receipts_decode_only_after_entry",
            "complete_result_and_failures_are_saved",
            "no_ancestor_import_or_reconstruction",
            "source_bind_after_freeze",
    ):
        require(policy.get(key) is True, "source policy guard " + key)

    rows = manifest.get("files")
    require(isinstance(rows, list) and len(rows) == len(REQUIRED), "source manifest file count")
    hashes: dict[str, str] = {}
    for row in rows:
        require(isinstance(row, Mapping), "source manifest row")
        path = row.get("path")
        require(isinstance(path, str) and path in REQUIRED and path not in hashes, "unexpected source file")
        require(not _forbidden(Path(path)), "forbidden source path")
        digest = sha_path(ROOT / path)
        require(digest == row.get("sha256"), "changed frozen source file " + path)
        if path in GZIP_INPUTS:
            require(row.get("raw_sha256_before_decode") is True, "receipt lacks raw preflight marker " + path)
            require(digest == RAW_SHA256[path], "changed pinned raw receipt " + path)
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), "incomplete source manifest")
    return hashes, sha(source_raw), sha(brief_raw), preimage_sha256


def _module(path: str, name: str) -> Any:
    require(path in (LEAF, FIXTURE), "unexpected import path: " + path)
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    require(spec is not None and spec.loader is not None, "module unavailable: " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def _leaf_callable(module: Any) -> Callable[..., Any]:
    candidate = getattr(module, "build", None)
    require(callable(candidate), "new P cell-sum leaf must expose build")
    return candidate


def _as_mapping(value: Any, label: str) -> Mapping[str, Any]:
    require(isinstance(value, Mapping), label + " must be a mapping")
    return value


def _term_ids(value: Mapping[str, Any]) -> list[str]:
    result: list[str] = []
    for holder in (value, value.get("origin")):
        if not isinstance(holder, Mapping):
            continue
        raw = holder.get("source_ids", holder.get("source_id", ()))
        if isinstance(raw, str):
            result.append(raw)
        elif isinstance(raw, Sequence) and not isinstance(raw, (str, bytes, bytearray)):
            result.extend(str(item) for item in raw)
    return result


def _term_classes(value: Mapping[str, Any]) -> list[Any]:
    origin = value.get("origin")
    if isinstance(origin, Mapping):
        classes = origin.get("error_classes", ())
    else:
        classes = value.get("error_classes", ())
    if isinstance(classes, str):
        return [classes]
    return list(classes) if isinstance(classes, Sequence) else []


def _input_index(output: Mapping[str, Any], label: str) -> Mapping[str, Any]:
    for key in ("input_index", "source_index"):
        value = output.get(key)
        if isinstance(value, Mapping):
            return value
    nested = output.get("result")
    if isinstance(nested, Mapping):
        for key in ("input_index", "source_index"):
            value = nested.get(key)
            if isinstance(value, Mapping):
                return value
    return {}


def _potential_terms(context: Any, output: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    potential = output.get("potential_compensation_m15")
    require(isinstance(potential, Mapping), "P receipt potential_compensation_m15 missing")
    require(set(potential) == set(BRANCHES), "P branch roster")
    packet: dict[str, Any] = {}
    paths: dict[str, Any] = {}
    for branch in BRANCHES:
        branch_raw = potential[branch]
        require(isinstance(branch_raw, Mapping), branch + " branch mapping")
        packet[branch] = {}
        for order in ORDERS:
            order_key = str(order)
            row = branch_raw.get(order_key, branch_raw.get(order))
            require(isinstance(row, Mapping), branch + " missing ordinary Z" + str(order))
            require(set(str(key) for key in row) == {str(cell) for cell in CELLS}, branch + " cell roster")
            packet[branch][order_key] = {}
            for cell in CELLS:
                cell_key = str(cell)
                raw_terms = row.get(cell_key, row.get(cell))
                require(isinstance(raw_terms, list), branch + " Z" + str(order) + " cell" + str(cell) + " terms")
                hydrated_terms: list[dict[str, Any]] = []
                for index, raw_term in enumerate(raw_terms):
                    term_token = _decode_exact(raw_term, f"P.{branch}.Z{order}.cell{cell}[{index}]")
                    term = _hydrate(context, term_token, f"P.{branch}.Z{order}.cell{cell}[{index}]")
                    require(isinstance(term, Mapping), "P term mapping")
                    powers = term.get("powers")
                    require(
                        isinstance(powers, Sequence)
                        and len(powers) == 4
                        and all(type(power) is int for power in powers),
                        "P term powers",
                    )
                    require(hasattr(term.get("coefficient"), "_mpi_"), "P coefficient interval")
                    origin = term.get("origin")
                    require(isinstance(origin, Mapping), "P term origin")
                    require(
                        isinstance(origin.get("source_ids"), Sequence)
                        and not isinstance(origin.get("source_ids"), (str, bytes, bytearray))
                        and bool(origin.get("source_ids")),
                        "P origin source_ids",
                    )
                    require(
                        isinstance(origin.get("error_classes"), Sequence)
                        and not isinstance(origin.get("error_classes"), (str, bytes, bytearray)),
                        "P origin error_classes",
                    )
                    require(isinstance(origin.get("uncertainty"), bool), "P origin uncertainty")
                    require(origin.get("target_cell") == cell, "P target_cell mismatch")
                    require(type(origin.get("lead_order")) is int and type(origin.get("word_order")) is int, "P derivative pair")
                    require(origin["lead_order"] in (0, 1) and origin["word_order"] in (0, 1), "P derivative pair range")
                    require(origin["lead_order"] + origin["word_order"] == order, "P derivative pair/order mismatch")
                    require(term.get("error_only", False) is not True and origin.get("error_only", False) is not True, "error-only term in finite P stream")
                    term_copy = dict(term)
                    term_copy["source_path"] = ["result", "output", "potential_compensation_m15", branch, order_key, cell_key, index]
                    term_copy["receipt"] = P_RECEIPT
                    term_copy["result_sha256"] = RESULT_SHA256[P_RECEIPT]
                    hydrated_terms.append(term_copy)
                    paths[f"{branch}:Z{order}:cell{cell}:{index}"] = {
                        "receipt": P_RECEIPT,
                        "result_sha256": RESULT_SHA256[P_RECEIPT],
                        "path": list(term_copy["source_path"]),
                        "source_ids": list(origin["source_ids"]),
                        "error_classes": list(origin["error_classes"]),
                        "uncertainty": origin["uncertainty"],
                        "target_cell": cell,
                        "lead_order": origin["lead_order"],
                        "word_order": origin["word_order"],
                    }
                packet[branch][order_key][cell_key] = hydrated_terms
    return packet, paths


def _receipt(raw: Mapping[str, bytes], path: str, component: str) -> dict[str, Any]:
    require(sha(raw[path]) == RAW_SHA256[path], "pinned raw receipt changed " + path)
    value = json.loads(gzip.decompress(raw[path]))
    require(isinstance(value, Mapping), path + " receipt envelope")
    result = value.get("result")
    require(isinstance(result, Mapping), path + " result missing")
    require(result.get("status") == "returned" and result.get("component") == component, path + " status/component")
    require(_canonical_result_sha(result) == RESULT_SHA256[path], path + " saved result SHA")
    return dict(value)


def _source_scope(output: Mapping[str, Any], label: str) -> tuple[Any, Any, Any, Mapping[str, Any]]:
    family = output.get("source_family")
    logs = output.get("source_logs")
    basis = output.get("source_power_basis")
    require(family is not None and logs is not None, label + " source family/logs")
    require(isinstance(logs, Mapping) and set(logs) == set(POWER_BASIS), label + " source atom names")
    require(tuple(basis or POWER_BASIS) == POWER_BASIS, label + " source power basis")
    index = _input_index(output, label)
    return family, logs, tuple(basis or POWER_BASIS), index


def _validate_members(actual: Mapping[str, Any], packet: Mapping[str, Any]) -> dict[str, int]:
    """Audit the canonical branch/order/cell member stream exactly once."""

    branches = actual.get("branches")
    require(isinstance(branches, Mapping), "leaf branches missing for member audit")
    expected: dict[tuple[str, int, int, int], Mapping[str, Any]] = {}
    for branch in BRANCHES:
        rows = packet[branch]
        for order in ORDERS:
            for cell in CELLS:
                terms = rows[str(order)][str(cell)]
                for index, term in enumerate(terms):
                    expected[(branch, order, cell, index)] = term
    seen: set[tuple[str, int, int, int]] = set()
    group_counts: dict[str, int] = {}
    member_count = 0
    for branch in BRANCHES:
        branch_value = branches.get(branch)
        require(isinstance(branch_value, Mapping), "leaf branch missing " + branch)
        orders = branch_value.get("ordinary_Z", branch_value.get("orders"))
        require(isinstance(orders, Mapping), "leaf ordinary-Z rows missing " + branch)
        for order in ORDERS:
            row = orders.get(order, orders.get(str(order)))
            require(isinstance(row, Mapping), "leaf Z row missing " + branch + " Z" + str(order))
            cells = row.get("cells", row.get("target_cells"))
            require(isinstance(cells, Mapping), "leaf cells missing " + branch + " Z" + str(order))
            for cell in CELLS:
                cell_value = cells.get(cell, cells.get(str(cell)))
                require(isinstance(cell_value, Mapping), "leaf cell missing")
                groups = cell_value.get("groups", cell_value.get("aggregates"))
                require(isinstance(groups, list), "leaf cell groups missing")
                group_counts[f"{branch}:Z{order}:cell{cell}"] = len(groups)
                for group in groups:
                    require(isinstance(group, Mapping), "leaf group mapping")
                    members = group.get("members")
                    require(isinstance(members, list), "leaf group members missing")
                    for member in members:
                        require(isinstance(member, Mapping), "leaf member mapping")
                        key = (
                            member.get("branch"),
                            member.get("ordinary_Z_order"),
                            member.get("target_cell"),
                            member.get("input_index", member.get("index")),
                        )
                        require(key in expected, "leaf member is not an input term")
                        require(key not in seen, "leaf input member duplicated")
                        require(member.get("record") is expected[key], "leaf member record not retained by reference")
                        require(member.get("source_ids") == expected[key]["origin"]["source_ids"], "member source ids changed")
                        require(member.get("error_classes") == expected[key]["origin"]["error_classes"], "member error classes changed")
                        require(member.get("source_path") == expected[key].get("source_path"), "member source path changed")
                        seen.add(key)
                        member_count += 1
    require(seen == set(expected), "not every P term appears exactly once")
    return {
        "expected_term_count": len(expected),
        "member_count": member_count,
        "group_count": sum(group_counts.values()),
        "group_counts": group_counts,
    }


def consume(hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    # Re-authenticate all frozen bytes after entry.  No gzip stream is opened
    # until after the independent fixture returns successfully.
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "changed after entry " + path)
        raw[path] = data
    for path in GZIP_INPUTS:
        require(sha(raw[path]) == RAW_SHA256[path], "changed pinned raw receipt " + path)
    progress["raw_gzip_authenticated_before_decode"] = True

    leaf_module = _module(LEAF, COMPONENT + "_leaf")
    fixture_module = _module(FIXTURE, COMPONENT + "_fixture")
    context_type = __import__("mpmath.ctx_iv", fromlist=["MPIntervalContext"]).MPIntervalContext
    context = context_type()
    context.prec = 1800
    fixture_runner = getattr(fixture_module, "run_fixture", getattr(fixture_module, "run", None))
    require(callable(fixture_runner), "fixture must export run_fixture")
    progress["fixture_started"] = True
    proof = fixture_runner(leaf_module, context)
    progress["fixture"] = encode(proof)
    require(isinstance(proof, Mapping) and proof.get("passed") is True, "synthetic fixture failed")
    checks = proof.get("checks")
    if isinstance(checks, Mapping):
        require(all(value is True for value in checks.values()), "synthetic fixture has failed checks")
    progress["fixture_completed_before_input_json_decode"] = True
    progress["input_json_decode_started"] = True

    p_receipt = _receipt(raw, P_RECEIPT, "n1_core_kchi_compensation_v2")
    prior_receipt = _receipt(raw, PRIOR_RECEIPT, "n1_w0_chi_ungrouped_v3")
    progress["p_receipt_decoded"] = True
    progress["prior_receipt_decoded"] = True
    p_output = _as_mapping(p_receipt["result"].get("output"), "P output")
    prior_output = _as_mapping(prior_receipt["result"].get("output"), "prior output")
    p_family, p_logs, p_basis, p_index = _source_scope(p_output, "P")
    prior_family, prior_logs, prior_basis, prior_index = _source_scope(prior_output, "prior")
    p_index = _hydrate(context, _decode_exact(p_index, "P.input_index"), "P.input_index")
    prior_index = _hydrate(context, _decode_exact(prior_index, "prior.input_index"), "prior.input_index")
    require(p_family == prior_family, "source family mismatch")
    for atom in POWER_BASIS:
        left, right = p_logs[atom], prior_logs[atom]
        require(isinstance(left, Mapping) and isinstance(right, Mapping), "source log packet")
        for endpoint in ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple"):
            require(endpoint in left and endpoint in right and left[endpoint] == right[endpoint], "source log exact endpoint mismatch: " + atom)
    require(tuple(p_basis) == tuple(prior_basis) == POWER_BASIS, "source power basis mismatch")
    require(isinstance(p_index, Mapping), "P input_index missing")
    for key in (
        "physical_scale_reapplied",
        "physical_factors_applied",
        "physical_weight_reapplied",
        "radial_weight_reapplied",
        "geometry_reapplied",
        "ledger_212_changed",
        "ledger_appended",
        "ledger_written",
    ):
        if key in prior_output:
            require(prior_output[key] is False, "prior scope changed: " + key)
    p_manifest = _source_rows(p_receipt, "P")
    prior_manifest = _source_rows(prior_receipt, "prior")
    require(
        p_manifest.get(FEEDBACK_RECEIPT) == FEEDBACK_RAW_SHA256
        and prior_manifest.get(FEEDBACK_RECEIPT) == FEEDBACK_RAW_SHA256,
        "common FEEDBACK source-manifest pin",
    )

    # The prior receipt is a parent reference only.  Preserve its complete
    # input-index object if present, but never expand/replay its formal rows.
    prior_parent_reference = {
        "receipt": PRIOR_RECEIPT,
        "raw_sha256": RAW_SHA256[PRIOR_RECEIPT],
        "result_sha256": RESULT_SHA256[PRIOR_RECEIPT],
        "output_path": ["result", "output"],
        "input_index_path": ["result", "output", "input_index"] if prior_index else None,
        "input_index_retained": bool(prior_index),
        "formal_parent_expanded": False,
    }
    potential, p_paths = _potential_terms(context, p_output)
    require(len(p_paths) == sum(len(potential[b][str(o)][str(c)]) for b in BRANCHES for o in ORDERS for c in CELLS), "P term path count")

    source_params = {
        key: p_output[key]
        for key in ("sc", "domain", "source_power_basis", "source_family", "source_logs")
        if key in p_output
    }
    caller_obligation = {
        "authenticated_branches": list(BRANCHES),
        "upstream_source_authenticated": True,
        "same_interval_context": True,
        "condition": "P finite branches are raw-authenticated and hydrated in one MPIntervalContext.prec=1800",
    }
    leaf = _leaf_callable(leaf_module)
    progress["actual_leaf_build_entered"] = True
    # Canonical API of n1_p_compensation_cell_sum: context first, branch packet
    # second, with one explicit caller-authentication obligation.  This is the
    # sole actual call; no retry/fallback is permitted.
    actual = leaf(context, potential, caller_obligation=caller_obligation)
    progress["actual_leaf_build_calls"] = 1
    require(isinstance(actual, Mapping), "P cell-sum leaf result mapping")
    require(actual.get("schema") == COMPONENT, "P cell-sum schema")
    require(actual.get("cross_cell_grouping") is False, "cross-cell grouping claimed")
    require(actual.get("cross_branch_grouping") is False, "cross-branch grouping claimed")
    require(actual.get("error_only_as_signed") is False, "error-only mixed into signed stream")
    require(actual.get("physical_scale_reapplied") is False, "physical scale reapplied")
    require(actual.get("physical_measure_reapplied") is False, "physical measure reapplied")
    require(actual.get("ledger_appended") is False, "ledger appended")
    require(actual.get("same_source_claimed") is False, "same source claimed")
    require(actual.get("same_gauge_claimed") is False, "same gauge claimed")
    require(actual.get("physical_field_claimed") is False, "physical field claimed")
    require(actual.get("tail_added") is False, "tail included")
    require(actual.get("input_records_retained_by_reference") is True, "input records not retained")
    require(actual.get("member_records_retained_by_reference") is True, "member records not retained")
    require(actual.get("source_paths_retained") is True, "source paths not retained")
    require(actual.get("provenance_retained") is True, "provenance not retained")
    require(actual.get("outward_interval_addition_only") is True, "non-directed operation")
    member_audit = _validate_members(actual, potential)
    output = dict(actual)
    output["input_index"] = encode(p_index)
    output["parent_input_index"] = encode(prior_index)
    output["potential_compensation_m15_input"] = encode(potential)
    output["potential_compensation_m15_input_paths"] = p_paths
    output["member_audit"] = member_audit
    output["source_family"] = encode(p_family)
    output["source_logs"] = encode(p_logs)
    output["source_power_basis"] = list(POWER_BASIS)
    output["retained_source_physical_params"] = encode(source_params)
    output["prior_source_physical_params"] = encode(
        {
            key: prior_output[key]
            for key in (
                "sc",
                "domain",
                "source_power_basis",
                "source_family",
                "source_logs",
                "physical_params",
                "physical_parameters",
            )
            if key in prior_output
        }
    )
    output["parent_receipt_reference"] = prior_parent_reference
    output["p_receipt_reference"] = {
        "path": P_RECEIPT,
        "raw_sha256": RAW_SHA256[P_RECEIPT],
        "result_sha256": RESULT_SHA256[P_RECEIPT],
        "potential_compensation_m15_path": ["result", "output", "potential_compensation_m15"],
    }
    output["grouping_key"] = "(branch, ordinary_Z_order, target_cell, lead_order, word_order, powers, uncertainty, ordered error_classes)"
    output["prior_direct_m15_weights_reapplied"] = False
    output["prior_direct_m15_weights_preserved"] = True
    output["prior_212_terms_added"] = False
    output["gauge_common_owner_proved"] = False
    output["tail_included"] = False
    output["P_only_candidate_non_point_intervals"] = True
    output["complete_m15_computed"] = False
    output["terminal_closure_proved"] = False
    output["ledger_212_changed"] = False
    output["ancestor_reexecuted"] = False
    output["actual_leaf_build_calls"] = progress["actual_leaf_build_calls"]
    output["fixture_leaf_build_calls"] = 1
    return output


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
            "ancestor_reexecuted": False,
        },
    )
    progress: dict[str, Any] = {
        "raw_gzip_authenticated_before_decode": False,
        "fixture_started": False,
        "fixture_completed_before_input_json_decode": False,
        "input_json_decode_started": False,
        "p_receipt_decoded": False,
        "prior_receipt_decoded": False,
        "actual_leaf_build_entered": False,
        "actual_leaf_build_calls": 0,
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

