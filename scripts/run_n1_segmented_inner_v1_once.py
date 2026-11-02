"""One-entry segmented inner consumer for the saved W0, nested K0-W1, and collar packets.

Every frozen byte is authenticated before a receipt is decoded.  The synthetic
fixture runs after the raw-byte check and before any gzip decode.  The leaf is
called once for the actual saved packet; fixture and actual-call progress,
including complete tracebacks on failure, are written to the exclusive result.
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
COMPONENT = "n1_segmented_inner_v1"
DIR = ROOT / "tmp" / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entered.json", "result.json"))

CORE = "experiments/root_st073/n1_core_w0_m15_v1_receipt.json.gz"
W1 = "experiments/root_st073/n1_core_k0_feedback_m15_v2_receipt.json.gz"
COLLAR = "experiments/root_st073/first_collar_m15_group_v2_receipt.json.gz"
OLD_COLLAR = "experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz"

LEAF = "src/openai_ns_reconstruction/n1_segmented_inner_candidate.py"
FIXTURE = "src/openai_ns_reconstruction/n1_segmented_inner_fixture_v1.py"
RUNNER = "scripts/run_n1_segmented_inner_v1_once.py"
REQUIRED = (LEAF, FIXTURE, RUNNER, CORE, W1, COLLAR)
RECEIPTS = frozenset((CORE, W1, COLLAR))

RAW_SHA256 = {
    CORE: "157f5677742b71a469c67e0efd21de0d94d230b4d3c2968d82acd2972413bda7",
    W1: "ca947df9fdf7765d9211f428375d81956af783d89074f9e782ba56aeff91dad0",
    COLLAR: "86eb1b46c0fc965145ab596242c55b39498507db39b230022bdca0b264ebca60",
}
RESULT_SHA256 = {
    CORE: "ed1cc60175f34625142e7e7b62a93514be93cae8eb70ba0634ed67b152decfd8",
    W1: "3f7f3398bcad74e43acfeebe6eef64371e5156e77240c2cbe34742b776450952",
    COLLAR: "a9c98be65d77208f778e00ccb7f701120ab414d157b140e9eaf927c3f9f1640f",
}
OLD_COLLAR_SHA256 = "c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660"

FORBIDDEN_BASENAMES = {
    "resume_first_switch_2026_10_17.json",
    "rm_loop_density_identity.json",
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
        json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")
    with path.open("xb") as handle:
        handle.write(raw)
    return sha(raw)


def _forbidden(path: Path) -> bool:
    name = path.name.lower()
    return name in FORBIDDEN_BASENAMES or name.endswith("rm_loop_density_identity.json")


def _reject_float_tree(value: Any, label: str) -> None:
    if isinstance(value, float):
        raise TypeError(label + " must not contain Python floating-point values")
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
        raise TypeError("cannot serialize a Python floating-point value")
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    raise TypeError("cannot serialize result type " + type(value).__name__)


def decode(c: MPIntervalContext, value: Any, label: str = "value") -> Any:
    """Rebuild every saved directed endpoint in the one 1800-bit context."""

    if isinstance(value, Mapping):
        has_lower = "lower_exact_mpf_tuple" in value
        has_upper = "upper_exact_mpf_tuple" in value
        if has_lower or has_upper:
            require(has_lower and has_upper, label + " has incomplete endpoints")
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
        raise TypeError(label + " contains a Python float")
    return value


def _lookup(mapping: Mapping[Any, Any], key: int | str, label: str) -> Any:
    if key in mapping:
        return mapping[key]
    alternate = str(key) if isinstance(key, int) else int(key)
    if alternate in mapping:
        return mapping[alternate]
    raise KeyError(label + " missing key " + str(key))


def _orders(mapping: Any, label: str) -> dict[int, Any]:
    require(isinstance(mapping, Mapping), label + " must be a mapping")
    result = {order: _lookup(mapping, order, label) for order in (0, 1)}
    require(len(mapping) == 2, label + " must contain exactly Z0 and Z1")
    return result


def _term_ids(term: Mapping[str, Any], label: str) -> list[str]:
    origin = term.get("origin")
    require(isinstance(origin, Mapping), label + ".origin missing")
    ids = origin.get("source_ids", origin.get("sources"))
    require(isinstance(ids, Sequence) and not isinstance(ids, (str, bytes, bytearray)) and ids, label + " source ids missing")
    require(isinstance(origin.get("uncertainty"), bool), label + " uncertainty missing")
    errors = origin.get("error_classes")
    require(isinstance(errors, Sequence) and not isinstance(errors, (str, bytes, bytearray)), label + " error classes missing")
    return [item for item in ids if isinstance(item, str) and item]


def _source_index(output: Mapping[str, Any], label: str, expected_count: int | None = None) -> Mapping[str, Any]:
    index = output.get("source_index", output.get("input_index"))
    require(isinstance(index, Mapping), label + " source index missing")
    if expected_count is not None:
        require(len(index) == expected_count, label + " source index count")
    return index


def module(path: str, name: str) -> Any:
    require(path in (LEAF, FIXTURE), "unexpected import path: " + path)
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    require(spec is not None and spec.loader is not None, "module unavailable: " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def _manifest_hashes(manifest: Any, label: str) -> dict[str, str]:
    require(isinstance(manifest, Mapping), label + " source manifest missing")
    rows = manifest.get("files")
    require(isinstance(rows, list), label + " source manifest files missing")
    result: dict[str, str] = {}
    for row in rows:
        require(isinstance(row, Mapping), label + " source manifest row")
        path, digest = row.get("path"), row.get("sha256")
        require(isinstance(path, str) and isinstance(digest, str) and path not in result, label + " source manifest path")
        result[path] = digest
    return result


def _canonical_result_sha(result: Mapping[str, Any]) -> str:
    raw = (
        json.dumps(result, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False)
        + "\n"
    ).encode("utf-8")
    return sha(raw)


def _receipt(
    raw: Mapping[str, bytes],
    path: str,
    component: str,
    result_sha: str,
    *,
    frozen_key: str,
) -> dict[str, Any]:
    # This helper is called only after the caller has authenticated the raw
    # gzip bytes.  It deliberately does not accept a filesystem path.
    value = json.loads(gzip.decompress(raw[path]))
    require(isinstance(value, Mapping), path + " receipt envelope")
    require(all(key in value for key in ("source", "entry", "result", frozen_key)), path + " receipt envelope fields")
    result = value.get("result")
    require(isinstance(result, Mapping), path + " result missing")
    require(result.get("status") == "returned" and result.get("component") == component, path + " published component lock")
    require(_canonical_result_sha(result) == result_sha, path + " saved result SHA")
    return dict(value)


def preflight() -> tuple[dict[str, str], str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    brief_raw = BRIEF.read_bytes()
    source_raw = SOURCE.read_bytes()
    brief = json.loads(brief_raw)
    source = json.loads(source_raw)
    require(
        brief.get("component") == COMPONENT
        and brief.get("phase") == COMPONENT + "_pending"
        and brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("one_entry_only") is True
        and brief.get("do_not_reenter") is True,
        "brief pending guard",
    )
    checkpoint = brief.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "source checkpoint missing")
    require(Path(checkpoint.get("path", "")).resolve() == SOURCE.resolve(), "source checkpoint path")
    require(checkpoint.get("sha256") == sha(source_raw), "source checkpoint SHA")
    preimage_path = Path(checkpoint.get("preimage_path", ""))
    require(not _forbidden(preimage_path), "forbidden source preimage")
    preimage_raw = preimage_path.read_bytes()
    require(sha(preimage_raw) == checkpoint.get("preimage_sha256"), "source preimage SHA")
    preimage = json.loads(preimage_raw)
    require(
        isinstance(preimage, Mapping)
        and preimage.get("component") == "first_collar_m15_group_v2"
        and preimage.get("phase") == "component_published",
        "published predecessor preimage guard",
    )
    require(
        source.get("component") == COMPONENT
        and source.get("one_entry_only") is True
        and source.get("previous_component") == "first_collar_m15_group_v2",
        "source manifest lock",
    )
    policy = source.get("source_policy", {})
    require(isinstance(policy, Mapping), "source policy missing")
    require(policy.get("actual_leaf_calls") == 1, "source policy leaf call count")
    require(policy.get("exclusive_entry_and_result") is True, "source policy exclusive outputs")
    require(policy.get("preserve_full_traceback") is True, "source policy traceback")
    require(policy.get("raw_hash_before_gzip_decode") is True, "source policy raw hash")
    require(policy.get("synthetic_fixture_before_gzip_decode") is True, "source policy fixture order")

    rows = source.get("files")
    require(isinstance(rows, list) and len(rows) == len(REQUIRED), "source manifest file count")
    hashes: dict[str, str] = {}
    for row in rows:
        require(isinstance(row, Mapping), "source manifest row")
        path = row.get("path")
        require(isinstance(path, str) and path in REQUIRED and path not in hashes, "unexpected source file")
        require(not _forbidden(Path(path)), "forbidden source file")
        digest = sha_path(ROOT / path)
        require(digest == row.get("sha256"), "changed source file " + path)
        if path in RECEIPTS:
            require(row.get("raw_sha256_before_decode") is True, "receipt not raw-pinned " + path)
            require(digest == RAW_SHA256[path], "changed pinned receipt " + path)
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), "incomplete source manifest")
    return hashes, sha(source_raw), sha(brief_raw)


def _hydrate_term(c: MPIntervalContext, raw_term: Any, label: str) -> dict[str, Any]:
    require(isinstance(raw_term, Mapping), label + " term mapping")
    term = decode(c, raw_term, label)
    require(isinstance(term, Mapping), label + " term mapping")
    require("coefficient" in term, label + " coefficient missing")
    require(hasattr(term["coefficient"], "_mpi_"), label + " coefficient is not directed")
    _term_ids(term, label)
    powers = term.get("powers")
    require(isinstance(powers, Sequence) and not isinstance(powers, (str, bytes, bytearray)) and len(powers) == 4, label + " powers")
    require(all(type(power) is int for power in powers), label + " powers are not integers")
    # F0_power/Lambda_power are part of the saved core contract and must not
    # be reconstructed from a shortened powers tuple.
    require("F0_power" in term and "Lambda_power" in term, label + " legacy formal powers missing")
    return dict(term)


def _term_rows(c: MPIntervalContext, raw_rows: Any, label: str) -> dict[int, list[dict[str, Any]]]:
    rows = _orders(raw_rows, label)
    result: dict[int, list[dict[str, Any]]] = {}
    for order in (0, 1):
        values = rows[order]
        require(isinstance(values, list), label + " Z" + str(order) + " must be a list")
        result[order] = [_hydrate_term(c, value, f"{label}[{order}][{index}]") for index, value in enumerate(values)]
    return result


def _core_packet(c: MPIntervalContext, receipt: Mapping[str, Any]) -> tuple[dict[str, Any], set[str]]:
    output = receipt["result"]["output"]
    require(isinstance(output, Mapping), "core output mapping")
    saved = output.get("result")
    require(isinstance(saved, Mapping), "core nested result mapping")
    require(saved.get("rho_formula") == "4*t^2", "core rho formula")
    require(output.get("domain") == "rho=4, real Z in [-1,1]", "core domain")
    require(output.get("core_W1_radial_history_missing") is True, "core full-W1 scope flag")
    require(output.get("pressure_direct_density_no_endpoint_added") is True, "core pressure scope flag")
    require(output.get("F0_connection_reused_without_recomputation") is True, "core F0 reuse flag")
    require(saved.get("W1_included") is False and saved.get("pressure_endpoint_added") is False, "core nested scope flags")
    index = _source_index(output, "core", expected_count=57)

    axial_raw = _orders(saved.get("axial"), "core result axial")
    swirl_raw = _orders(saved.get("swirl"), "core result swirl")
    axial: dict[int, dict[int, list[dict[str, Any]]]] = {0: {}, 1: {}}
    swirl: dict[int, dict[int, list[dict[str, Any]]]] = {0: {}, 1: {}}
    cells: list[dict[str, Any]] = []
    for cell in range(8):
        channels: dict[str, dict[int, list[dict[str, Any]]]] = {"axial": {0: [], 1: []}, "swirl": {0: [], 1: []}}
        for channel, raw_rows, target in (("axial", axial_raw, axial), ("swirl", swirl_raw, swirl)):
            for order in (0, 1):
                row = raw_rows[order]
                require(isinstance(row, Mapping) or isinstance(row, list), "core " + channel + " row")
                raw_terms = _lookup(row, cell, "core " + channel + " cell") if isinstance(row, Mapping) else row[cell]
                require(isinstance(raw_terms, list), "core " + channel + " cell terms")
                terms = [_hydrate_term(c, term, f"core.{channel}[{order}][{cell}][{i}]") for i, term in enumerate(raw_terms)]
                for term_index, term in enumerate(terms):
                    term["word"] = "W0"
                    term["segment"] = "core_w0"
                    term["source_receipt"] = CORE
                    term["receipt_input_path"] = [CORE, "result", "output", "result", channel, str(order), cell, term_index]
                target[order][cell] = terms
                channels[channel][order] = terms
        cells.append({"cell": cell, "channels": channels, "word": "W0", "segment": "core_w0"})
    pressure = _term_rows(c, saved.get("pressure"), "core result pressure")
    for order in (0, 1):
        for term_index, term in enumerate(pressure[order]):
            term["word"] = "W0"
            term["segment"] = "core_w0_pressure"
            term["source_receipt"] = CORE
            term["receipt_input_path"] = [CORE, "result", "output", "result", "pressure", str(order), term_index]
    packet = {
        "schema": "n1_segmented_inner_core_w0",
        "cells": cells,
        "axial": axial,
        "swirl": swirl,
        "pressure": pressure,
        "source_index": decode(c, index, "core input index"),
        "source_receipt": CORE,
        "source_result_sha256": RESULT_SHA256[CORE],
        "rho_formula": saved["rho_formula"],
        "domain": output["domain"],
        "physically_scaled": True,
        "full_original_core_W1_missing": True,
        "nested_K0_W1_saved_not_full_W1": True,
        "allword_tail_added": output.get("allword_tail_added") is True,
    }
    used = {sid for channel in ("axial", "swirl") for rows in packet[channel].values() for terms in rows.values() for term in terms for sid in _term_ids(term, "core term")}
    used.update(sid for terms in pressure.values() for term in terms for sid in _term_ids(term, "core pressure term"))
    require(used.issubset(set(index)), "core source index lost a term")
    packet["actual_source_ids"] = used
    return packet, used


def _w1_packet(c: MPIntervalContext, receipt: Mapping[str, Any]) -> tuple[dict[str, Any], set[str]]:
    output = receipt["result"]["output"]
    require(isinstance(output, Mapping), "nested W1 output mapping")
    m15 = output.get("m15")
    require(isinstance(m15, Mapping), "nested W1 m15 mapping")
    require("axial" in m15 and "swirl" in m15 and "pressure" not in m15, "nested W1 pressure exclusion")
    require(output.get("pressure_added", output.get("g4_pressure_added", False)) is False, "nested W1 pressure flag")
    index = _source_index(output, "nested W1", expected_count=461)
    axial_raw = _orders(m15["axial"], "nested W1 axial")
    swirl_raw = _orders(m15["swirl"], "nested W1 swirl")
    axial: dict[int, dict[int, list[dict[str, Any]]]] = {0: {}, 1: {}}
    swirl: dict[int, dict[int, list[dict[str, Any]]]] = {0: {}, 1: {}}
    cells: list[dict[str, Any]] = []
    for cell in range(8):
        channels: dict[str, dict[int, list[dict[str, Any]]]] = {"axial": {0: [], 1: []}, "swirl": {0: [], 1: []}}
        for channel, raw_rows, target in (("axial", axial_raw, axial), ("swirl", swirl_raw, swirl)):
            for order in (0, 1):
                row = raw_rows[order]
                raw_terms = _lookup(row, cell, "nested W1 " + channel + " cell") if isinstance(row, Mapping) else row[cell]
                require(isinstance(raw_terms, list), "nested W1 cell term list")
                terms = [_hydrate_term(c, term, f"nested_w1.{channel}[{order}][{cell}][{i}]") for i, term in enumerate(raw_terms)]
                for term_index, term in enumerate(terms):
                    term["word"] = "W1-K0"
                    term["segment"] = "core_w1_k0"
                    term["source_receipt"] = W1
                    term["receipt_input_path"] = [W1, "result", "output", "m15", channel, str(order), cell, term_index]
                target[order][cell] = terms
                channels[channel][order] = terms
        cells.append({"cell": cell, "channels": channels, "word": "W1-K0", "segment": "core_w1_k0"})
    packet = {
        "schema": "n1_segmented_inner_nested_k0_w1",
        "cells": cells,
        "axial": axial,
        "swirl": swirl,
        "pressure": {},
        "source_index": decode(c, index, "nested W1 input index"),
        "source_receipt": W1,
        "source_result_sha256": RESULT_SHA256[W1],
        "ordinary_Z_orders": [0, 1],
        "pressure_added": False,
        "full_W1": False,
    }
    used = {sid for channel in ("axial", "swirl") for rows in packet[channel].values() for terms in rows.values() for term in terms for sid in _term_ids(term, "nested W1 term")}
    require(used.issubset(set(index)), "nested W1 source index lost a term")
    packet["actual_source_ids"] = used
    return packet, used


def _validate_prior_groups(groups: Mapping[int, list[Any]], *, expected_members: int) -> tuple[set[str], set[tuple[str, int, int]]]:
    source_ids: set[str] = set()
    member_keys: set[tuple[str, int, int]] = set()
    count = 0
    for order in (0, 1):
        for group in groups[order]:
            require(isinstance(group, Mapping), "collar prior group mapping")
            origin = group.get("origin")
            members = group.get("members")
            require(isinstance(origin, Mapping) and isinstance(members, list), "collar group provenance")
            require(isinstance(origin.get("source_ids"), Sequence) and isinstance(origin.get("error_classes"), Sequence), "collar group sources/errors")
            require(origin.get("member_count") == len(members), "collar group member count")
            for index, member in enumerate(members):
                require(isinstance(member, Mapping), "collar group member")
                channel = member.get("channel")
                member_order = member.get("ordinary_Z_order")
                member_index = member.get("index")
                key = (str(channel), int(member_order), int(member_index))
                require(key not in member_keys, "duplicate collar group member")
                member_keys.add(key)
                record = member.get("record")
                provenance = member.get("provenance")
                require(isinstance(record, Mapping) and isinstance(provenance, Mapping), "collar member record/provenance")
                require(provenance == record.get("origin"), "collar member provenance changed")
                ids = _term_ids(record, "collar member")
                source_ids.update(ids)
                count += 1
    require(count == expected_members, "collar prior group member total")
    return source_ids, member_keys


def _collar_packet(c: MPIntervalContext, receipt: Mapping[str, Any]) -> tuple[dict[str, Any], set[str]]:
    output = receipt["result"]["output"]
    require(isinstance(output, Mapping), "collar output mapping")
    groups_raw = _orders(output.get("groups"), "collar output groups")
    groups: dict[int, list[Any]] = {0: [], 1: []}
    for order in (0, 1):
        require(isinstance(groups_raw[order], list), "collar group order list")
        groups[order] = [decode(c, group, f"collar groups[{order}][{i}]") for i, group in enumerate(groups_raw[order])]
    # The published v2 output has all three indexes.  They are preserved as
    # a bundle, rather than reduced to only the IDs consumed by the leaf.
    for key in ("input_index", "prior_input_index", "group_input_index"):
        require(isinstance(output.get(key), Mapping), "collar output index " + key)
    require(len(output["input_index"]) == 1975 and len(output["prior_input_index"]) == 1975, "collar index count")
    require(sum(len(groups[order]) for order in (0, 1)) == 1461, "collar prior group count")
    member_source_ids, _ = _validate_prior_groups(groups, expected_members=1975)
    index = decode(c, output["input_index"], "collar input index")
    require(member_source_ids.issubset(set(output["group_input_index"])), "collar group source index lost a member")
    packet = {
        "schema": "n1_segmented_inner_prior_collar_groups",
        "groups": groups,
        "source_index": index,
        "input_index": index,
        "prior_input_index": decode(c, output["prior_input_index"], "collar prior index"),
        "group_input_index": decode(c, output["group_input_index"], "collar group index"),
        "source_receipt": COLLAR,
        "source_result_sha256": RESULT_SHA256[COLLAR],
        "domain": output.get("domain"),
        "rho_formula": "4*exp(hb*sc*q/4)",
        "physically_scaled": output.get("already_physically_scaled_channels") is True,
        "physical_factors_reapplied": output.get("physical_factors_reapplied") is True,
        "allword_tail_added": output.get("allword_tail_added") is True,
    }
    packet["actual_source_ids"] = member_source_ids
    return packet, member_source_ids


def _flat_core_rows(core: Mapping[str, Any], nested_w1: Mapping[str, Any]) -> dict[int, list[dict[str, Any]]]:
    """Flatten saved cell/channel hierarchy only at the leaf boundary."""

    rows: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    for channel in ("axial", "swirl"):
        for order in (0, 1):
            for cell in range(8):
                for term in core[channel][order][cell]:
                    record = dict(term)
                    record.update({
                        "segment": "core",
                        "word": "W0",
                        "receipt": CORE,
                        "path": list(term["receipt_input_path"]),
                        "channel": channel,
                        "cell": cell,
                    })
                    rows[order].append(record)
                for term in nested_w1[channel][order][cell]:
                    record = dict(term)
                    record.update({
                        "segment": "core",
                        "word": "W1-K0",
                        "receipt": W1,
                        "path": list(term["receipt_input_path"]),
                        "channel": channel,
                        "cell": cell,
                    })
                    rows[order].append(record)
    for order in (0, 1):
        for term in core["pressure"][order]:
            record = dict(term)
            record.update({
                "segment": "core",
                "word": "W0",
                "receipt": CORE,
                "path": list(term["receipt_input_path"]),
                "channel": "pressure",
                "cell": None,
            })
            rows[order].append(record)
    return rows


def _flat_collar_rows(collar: Mapping[str, Any]) -> dict[int, list[dict[str, Any]]]:
    rows: dict[int, list[dict[str, Any]]] = {0: [], 1: []}
    for order in (0, 1):
        for group_index, group in enumerate(collar["groups"][order]):
            # Each published group is one new leaf input.  Its already-summed
            # coefficient is consumed directly; the complete prior member
            # list, paths, and provenance remain nested in this record.
            record = dict(group)
            record.update({
                "segment": "collar",
                "word": "collar-prior-group",
                "receipt": COLLAR,
                "path": [COLLAR, "result", "output", "groups", str(order), group_index],
                "channel": "collar_group",
                "cell": None,
            })
            rows[order].append(record)
    return rows


def _source_ids(value: Any) -> set[str]:
    found: set[str] = set()
    if isinstance(value, Mapping):
        for key, item in value.items():
            if key in ("source_ids", "sources") and isinstance(item, Sequence) and not isinstance(item, (str, bytes, bytearray)):
                found.update(x for x in item if isinstance(x, str))
            found.update(_source_ids(item))
    elif isinstance(value, (list, tuple)):
        for item in value:
            found.update(_source_ids(item))
    return found


def _segments(result: Mapping[str, Any]) -> tuple[Any, Any]:
    segments = result.get("segments")
    if isinstance(segments, Mapping):
        return segments.get("C0"), segments.get("Z1")
    return result.get("C0"), result.get("Z1")


def _path_key(value: Any) -> Any:
    if isinstance(value, Mapping):
        return tuple(sorted((str(key), _path_key(item)) for key, item in value.items()))
    if isinstance(value, (list, tuple)):
        return tuple(_path_key(item) for item in value)
    return (type(value).__name__, value)


def _row_key(segment: str, order: int, record: Mapping[str, Any]) -> tuple[Any, ...]:
    return (segment, order, str(record.get("receipt")), _path_key(record.get("path")))


def _member_key(segment: str, order: int, member: Mapping[str, Any]) -> tuple[Any, ...]:
    return (segment, order, str(member.get("receipt")), _path_key(member.get("path")))


def _namespace_for_record(segment: str, record: Mapping[str, Any]) -> str:
    if segment == "collar":
        return "collar"
    word = record.get("word", "")
    return "core_w1" if isinstance(word, str) and "w1" in word.lower() else "core_w0"


def _validate_result_members(
    result: Mapping[str, Any],
    core_rows: Mapping[int, Sequence[Any]],
    collar_rows: Mapping[int, Sequence[Any]],
) -> tuple[dict[str, set[str]], dict[str, Any]]:
    """Check that flattened input members survive exactly once in groups."""

    expected: dict[tuple[Any, ...], Mapping[str, Any]] = {}
    input_counts: dict[str, dict[str, int]] = {
        "core_w0": {"0": 0, "1": 0},
        "core_w1": {"0": 0, "1": 0},
        "collar": {"0": 0, "1": 0},
    }
    prior_group_count = 0
    prior_member_count = 0
    for segment, rows in (("core", core_rows), ("collar", collar_rows)):
        for order in (0, 1):
            for record in rows[order]:
                require(isinstance(record, Mapping), "flattened input record mapping")
                key = _row_key(segment, order, record)
                require(key not in expected, "duplicate flattened input member path")
                expected[key] = record
                input_counts[_namespace_for_record(segment, record)][str(order)] += 1
                if segment == "collar":
                    members = record.get("members")
                    require(isinstance(members, list), "flattened collar prior group members")
                    prior_group_count += 1
                    prior_member_count += len(members)
    require(prior_group_count == 1461, "flattened collar prior group count")
    require(prior_member_count == 1975, "flattened collar prior member count")

    grouped = result.get("segment_groups", result.get("groups_by_segment"))
    require(isinstance(grouped, Mapping), "candidate segment groups missing")
    actual: dict[tuple[Any, ...], Mapping[str, Any]] = {}
    actual_ids: dict[str, set[str]] = {"core_w0": set(), "core_w1": set(), "collar": set()}
    group_counts: dict[str, dict[str, int]] = {"core": {"0": 0, "1": 0}, "collar": {"0": 0, "1": 0}}
    member_counts: dict[str, dict[str, int]] = {
        "core_w0": {"0": 0, "1": 0},
        "core_w1": {"0": 0, "1": 0},
        "collar": {"0": 0, "1": 0},
    }
    for segment in ("core", "collar"):
        rows = grouped.get(segment)
        require(isinstance(rows, Mapping), "candidate " + segment + " groups mapping")
        for order in (0, 1):
            row = _lookup(rows, order, "candidate " + segment + " groups")
            require(isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray)), "candidate group row list")
            group_counts[segment][str(order)] = len(row)
            for group in row:
                require(isinstance(group, Mapping), "candidate grouped record mapping")
                members = group.get("members")
                require(isinstance(members, Sequence) and not isinstance(members, (str, bytes, bytearray)), "candidate group members list")
                for member in members:
                    require(isinstance(member, Mapping), "candidate member mapping")
                    key = _member_key(segment, order, member)
                    require(key not in actual, "candidate member path repeated")
                    original = member.get("record")
                    require(isinstance(original, Mapping), "candidate original member record")
                    expected_record = expected.get(key)
                    require(expected_record is not None, "candidate member path not in flattened input")
                    require(original is expected_record, "candidate member record was copied or replaced")
                    require(member.get("origin") is expected_record.get("origin"), "candidate member origin changed")
                    require(member.get("provenance") is expected_record.get("origin"), "candidate member provenance changed")
                    require(original.get("coefficient") is expected_record.get("coefficient"), "candidate member coefficient changed")
                    actual[key] = member
                    namespace = _namespace_for_record(segment, expected_record)
                    member_counts[namespace][str(order)] += 1
                    actual_ids[namespace].update(_source_ids(expected_record.get("origin")))
    require(set(actual) == set(expected), "candidate member path set differs from flattened input")
    return actual_ids, {
        "input_counts": input_counts,
        "candidate_member_counts": member_counts,
        "candidate_group_counts": group_counts,
        "prior_group_count": prior_group_count,
        "prior_member_count": prior_member_count,
        "expected_member_count": len(expected),
        "candidate_member_count": len(actual),
        "member_paths_exactly_once": True,
        "origin_and_coefficient_identity_retained": True,
    }


def _validate_leaf_result(
    result: Any,
    namespace_indexes: Mapping[str, Mapping[str, Any]],
    family: Any,
    logs: Any,
    core_rows: Mapping[int, Sequence[Any]],
    collar_rows: Mapping[int, Sequence[Any]],
    progress: dict[str, Any],
) -> None:
    require(isinstance(result, Mapping), "segmented leaf result mapping")
    require(result.get("schema") == COMPONENT, "segmented leaf schema")
    require(result.get("source_family") == family, "segmented source family changed")
    require(result.get("source_logs") == logs, "segmented source logs changed")
    c0, z1 = _segments(result)
    require(c0 is not None and z1 is not None and bool(c0) and bool(z1), "partial C0/Z1 segments missing")
    result_ids, member_stats = _validate_result_members(result, core_rows, collar_rows)
    result_index = result.get("source_index", result.get("input_index"))
    require(isinstance(result_index, Mapping), "segmented source index incomplete")
    core_index = result_index.get("core")
    collar_index = result_index.get("collar")
    require(isinstance(core_index, Mapping) and isinstance(collar_index, Mapping), "segmented source index namespaces incomplete")
    core_actual = set(namespace_indexes["core_w0"]["actual_ids"]) | set(namespace_indexes["core_w1"]["actual_ids"])
    require(core_actual.issubset(set(core_index)), "segmented core source index incomplete")
    require(set(namespace_indexes["collar"]["actual_ids"]).issubset(set(collar_index)), "segmented collar source index incomplete")
    for namespace, details in namespace_indexes.items():
        actual_ids = set(details["actual_ids"])
        prior_index = details["prior_index"]
        require(actual_ids.issubset(set(prior_index)), namespace + " actual IDs are not covered by its prior source index")
        require(actual_ids.issubset(result_ids[namespace]), namespace + " member IDs missing from candidate groups")
    progress["member_validation"] = member_stats
    require(result.get("physical_scale_reapplied", False) is False, "physical scale reapplied")
    require(result.get("endpoint_cancellation_proved", False) is False, "endpoint cancellation claimed")
    require(result.get("terminal_closure_proved", False) is False, "terminal closure claimed")
    require(result.get("ledger_212_changed", False) is False and result.get("ledger_212_consulted", False) is False, "212 ledger changed or consulted")
    require(result.get("P1_is_moment", False) is False and result.get("pressure_is_moment", False) is False, "P1/moment alias")
    require(result.get("allword_tail_added", False) is False, "all-word tail added")
    progress["leaf_result_validated"] = True


def _compact_prior_group_record(record: Mapping[str, Any]) -> dict[str, Any]:
    """Retain a prior collar group without expanding alias fields."""

    return {
        "segment": record.get("segment"),
        "word": record.get("word"),
        "channel": record.get("channel"),
        "cell": record.get("cell"),
        "receipt": record.get("receipt"),
        "path": record.get("path"),
        "ordinary_Z_order": record.get("ordinary_Z_order"),
        "powers": record.get("powers"),
        "coefficient": encode(record.get("coefficient")),
        "group_key": encode(record.get("group_key")),
        "origin": encode(record.get("origin")),
        "members": [
            _compact_member(member)
            for member in record.get("members", ())
            if isinstance(member, Mapping)
        ],
    }


def _compact_member(member: Mapping[str, Any]) -> dict[str, Any]:
    # Keep the original record and both provenance spellings, but do not walk
    # any candidate alias that points back to a group.
    return {
        "segment": member.get("segment"),
        "word": member.get("word"),
        "channel": member.get("channel"),
        "cell": member.get("cell"),
        "receipt": member.get("receipt"),
        "path": member.get("path"),
        "ordinary_Z_order": member.get("ordinary_Z_order"),
        "index": member.get("index"),
        "record_index": member.get("record_index"),
        "record": (
            _compact_prior_group_record(member["record"])
            if isinstance(member.get("record"), Mapping) and isinstance(member["record"].get("members"), list)
            else encode(member.get("record"))
        ),
        "origin": encode(member.get("origin")),
        "provenance": encode(member.get("provenance")),
    }


def _compact_group(group: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "segment": group.get("segment"),
        "ordinary_Z_order": group.get("ordinary_Z_order"),
        "powers": group.get("powers"),
        "coefficient": encode(group.get("coefficient")),
        "group_key": encode(group.get("group_key")),
        "origin": encode(group.get("origin")),
        "members": [_compact_member(member) for member in group.get("members", ()) if isinstance(member, Mapping)],
    }


def _compact_groups(grouped: Any) -> dict[str, dict[str, list[dict[str, Any]]]]:
    result: dict[str, dict[str, list[dict[str, Any]]]] = {}
    if not isinstance(grouped, Mapping):
        return result
    for segment, rows in grouped.items():
        if not isinstance(rows, Mapping):
            continue
        result[str(segment)] = {}
        for order in (0, 1):
            row = rows.get(order, rows.get(str(order), ()))
            result[str(segment)][str(order)] = [
                _compact_group(group)
                for group in row
                if isinstance(group, Mapping)
            ]
    return result


def _compact_formal(formal: Any) -> dict[str, list[dict[str, Any]]]:
    result: dict[str, list[dict[str, Any]]] = {"0": [], "1": []}
    if not isinstance(formal, Mapping):
        return result
    # group_refs are retained as segment/order/index references, never as the
    # aliased group dictionaries themselves.
    group_indices: dict[tuple[str, int], dict[int, int]] = {}
    for segment in ("core", "collar"):
        rows = formal.get("_segment_groups")
        if isinstance(rows, Mapping):
            segment_rows = rows.get(segment, {})
            group_indices[(segment, 0)] = {id(group): i for i, group in enumerate(segment_rows.get(0, segment_rows.get("0", ())))}
            group_indices[(segment, 1)] = {id(group): i for i, group in enumerate(segment_rows.get(1, segment_rows.get("1", ())))}
    for order in (0, 1):
        row = formal.get(order, formal.get(str(order), ()))
        for term in row:
            if not isinstance(term, Mapping):
                continue
            refs: dict[str, Any] = {}
            raw_refs = term.get("group_refs", {})
            for segment in ("core", "collar"):
                ref = raw_refs.get(segment) if isinstance(raw_refs, Mapping) else None
                refs[segment] = None if ref is None else {
                    "segment": segment,
                    "ordinary_Z_order": order,
                    "group_index": group_indices.get((segment, order), {}).get(id(ref)),
                }
            origin = term.get("origin", {})
            if isinstance(origin, Mapping):
                origin = {
                    key: origin.get(key)
                    for key in (
                        "operation", "segments", "ordinary_Z_order", "source_ids",
                        "source_ids_by_segment", "uncertainty", "uncertainty_by_segment",
                        "error_classes", "error_classes_by_segment", "source_group_count",
                        "member_count", "provenance_retained", "source_equality_proved",
                        "endpoint_cancellation_claimed",
                    )
                    if key in origin
                }
            result[str(order)].append(
                {
                    "segment": "segmented",
                    "segments": list(term.get("segments", ())),
                    "ordinary_Z_order": order,
                    "powers": term.get("powers"),
                    "coefficient": encode(term.get("coefficient")),
                    "word_by_segment": encode(term.get("word_by_segment", {})),
                    "channel_by_segment": encode(term.get("channel_by_segment", {})),
                    "cell_by_segment": encode(term.get("cell_by_segment", {})),
                    "receipt_by_segment": encode(term.get("receipt_by_segment", {})),
                    "path_by_segment": encode(term.get("path_by_segment", {})),
                    "group_refs": refs,
                    "group_key": encode(term.get("group_key")),
                    "origin": encode(origin),
                }
            )
    return result


def _compact_leaf_result(actual: Mapping[str, Any], enriched: Mapping[str, Any], core: Mapping[str, Any], nested_w1: Mapping[str, Any], collar: Mapping[str, Any]) -> dict[str, Any]:
    grouped = actual.get("segment_groups", actual.get("groups_by_segment", {}))
    # Attach grouped data transiently so _compact_formal can resolve refs by
    # identity, then remove the helper before serializing.
    formal = actual.get("formal_sum", {})
    if isinstance(formal, Mapping):
        formal_for_compact = dict(formal)
        formal_for_compact["_segment_groups"] = grouped
    else:
        formal_for_compact = {}
    return {
        "schema": COMPONENT,
        "source_family": encode(enriched.get("source_family")),
        "source_logs": encode(enriched.get("source_logs")),
        "source_atoms": encode(actual.get("source_atoms")),
        "source_power_basis": encode(actual.get("source_power_basis")),
        "ordinary_Z_orders": [0, 1],
        "segment_groups": _compact_groups(grouped),
        "formal_sum": _compact_formal(formal_for_compact),
        "partial_segment_labels": {"C0": "formal_sum[0]", "Z1": "formal_sum[1]"},
        "source_index": {
            "core_w0": encode(core.get("source_index")),
            "core_w1": encode(nested_w1.get("source_index")),
            "collar": encode(collar.get("source_index")),
        },
        "scope": encode(actual.get("scope")),
        "geometry_certificate": encode(actual.get("geometry_certificate")),
        "owner_certificate": encode(actual.get("owner_certificate")),
        "scale_certificate": encode(actual.get("scale_certificate")),
        "records_already_physically_weighted": None,
        "physical_weighting_is_caller_obligation": True,
        "parent_receipts_describe_inputs_as_physically_weighted": True,
        "physical_factors_applied": False,
        "physical_scale_reapplied": False,
        "pressure_endpoint_added": False,
        "P1_is_moment": False,
        "source_owner_identity_proved": False,
        "full_original_core_W1_missing": True,
        "nested_K0_W1_saved_not_full_W1": True,
        "Kchi_minus_K0_missing": True,
        "feedback_missing": True,
        "allword_tail_missing": True,
        "endpoint_cancellation_proved": False,
        "terminal_closure_proved": False,
        "allword_tail_added": False,
        "ledger_212_changed": False,
    }


def _compact_cell_packet(packet: Mapping[str, Any]) -> dict[str, Any]:
    cells: list[dict[str, Any]] = []
    for cell in packet.get("cells", ()):
        if not isinstance(cell, Mapping):
            continue
        channels: dict[str, dict[str, list[Any]]] = {}
        for channel, rows in (cell.get("channels", {}) or {}).items():
            if not isinstance(rows, Mapping):
                continue
            channels[str(channel)] = {
                str(order): [encode(term) for term in terms]
                for order, terms in rows.items()
            }
        cells.append({
            "cell": cell.get("cell"),
            "word": cell.get("word"),
            "segment": cell.get("segment"),
            "channels": channels,
        })
    pressure = {
        str(order): [encode(term) for term in terms]
        for order, terms in packet.get("pressure", {}).items()
    } if isinstance(packet.get("pressure"), Mapping) else {}
    return {
        "schema": packet.get("schema"),
        "source_receipt": packet.get("source_receipt"),
        "source_result_sha256": packet.get("source_result_sha256"),
        "rho_formula": packet.get("rho_formula"),
        "domain": packet.get("domain"),
        "cells": cells,
        "pressure": pressure,
        "source_index": encode(packet.get("source_index")),
        "full_original_core_W1_missing": packet.get("full_original_core_W1_missing") is True,
        "nested_K0_W1_saved_not_full_W1": packet.get("nested_K0_W1_saved_not_full_W1") is True,
    }


def consume(hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    # Re-read every frozen byte and authenticate its raw hash before importing
    # a new module or decoding any saved receipt.
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
    fixture_runner = getattr(fixture, "run_fixture", None)
    require(callable(fixture_runner), "fixture must export run_fixture")
    try:
        proof = fixture_runner(leaf, context)
        progress["fixture_leaf_calls"] = int(proof.get("leaf_calls", 1)) if isinstance(proof, Mapping) else 1
        progress["fixture"] = encode(proof)
        require(isinstance(proof, Mapping) and proof.get("passed") is True, "synthetic fixture failed")
        checks = proof.get("checks")
        if isinstance(checks, Mapping):
            require(all(checks.values()), "synthetic fixture has failed checks")
    except BaseException as exc:
        progress["fixture_error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
        raise
    progress["fixture_completed_before_input_json_decode"] = True

    # Only now are the three pinned gzip receipts decoded.
    progress["input_json_decode_started"] = True
    core_receipt = _receipt(raw, CORE, "n1_core_w0_m15_v1", RESULT_SHA256[CORE], frozen_key="frozen_byte_sha256")
    progress["core_receipt_decoded"] = True
    w1_receipt = _receipt(raw, W1, "n1_core_k0_feedback_m15_v2", RESULT_SHA256[W1], frozen_key="frozen_byte_sha256")
    progress["nested_w1_receipt_decoded"] = True
    collar_receipt = _receipt(raw, COLLAR, "first_collar_m15_group_v2", RESULT_SHA256[COLLAR], frozen_key="frozen_code_base64")
    progress["collar_receipt_decoded"] = True

    core_result = core_receipt["result"]
    w1_result = w1_receipt["result"]
    collar_result = collar_receipt["result"]
    core_output = core_result["output"]
    w1_output = w1_result["output"]
    collar_output = collar_result["output"]
    require(isinstance(core_output, Mapping) and isinstance(w1_output, Mapping) and isinstance(collar_output, Mapping), "prior outputs mapping")

    # Exact family/log equality is a compatibility condition only.  It is not
    # treated as proof that this new leaf owns any source atom.
    family = core_output.get("source_family")
    logs = core_output.get("source_logs")
    require(family == w1_output.get("source_family") == collar_output.get("source_family"), "source family mismatch")
    require(logs == w1_output.get("source_logs") == collar_output.get("source_logs"), "source logs mismatch")

    core_manifest = _manifest_hashes(core_receipt.get("source"), "core receipt")
    w1_manifest = _manifest_hashes(w1_receipt.get("source"), "nested W1 receipt")
    collar_manifest = _manifest_hashes(collar_receipt.get("source"), "collar receipt")
    require(core_manifest.get(OLD_COLLAR) == OLD_COLLAR_SHA256, "core old collar linkage")
    require(collar_manifest.get(OLD_COLLAR) == OLD_COLLAR_SHA256, "collar source manifest old collar linkage")
    require(w1_manifest.get(CORE) == RAW_SHA256[CORE], "nested W1 source manifest W0 linkage")
    require(not any("n1_first_collar_g4_pressure_v1_receipt" in path for path in w1_manifest), "nested W1 source manifest contains pressure receipt")
    require(core_output.get("previous_collar_receipt") == OLD_COLLAR, "core previous collar path")

    # Hydration and source-shape checks all use the same 1800-bit context.
    core_packet, core_ids = _core_packet(context, core_receipt)
    w1_packet, w1_ids = _w1_packet(context, w1_receipt)
    collar_packet, collar_ids = _collar_packet(context, collar_receipt)
    namespace_indexes = {
        "core_w0": {"actual_ids": core_ids, "prior_index": core_packet["source_index"]},
        "core_w1": {"actual_ids": w1_ids, "prior_index": w1_packet["source_index"]},
        "collar": {"actual_ids": collar_ids, "prior_index": collar_packet["group_input_index"]},
    }
    source_logs = decode(context, logs, "source_logs")
    source_family = decode(context, family, "source_family")
    geometry = {
        "core_rho_formula": core_packet["rho_formula"],
        "collar_rho_formula": "4*exp(hb*sc*q/4)",
        "core_domain": core_packet["domain"],
        "collar_domain": collar_output.get("domain"),
        "rho_seam": True,
        "seam_scope": "geometric interface only; no endpoint cancellation",
        "core_endpoint_rho": "rho=4 at t=1",
        "collar_endpoint_rho": "rho=4 at q=0",
        "owner_identity_proved": False,
    }
    # The candidate's public API is row-major, while the receipt evidence is
    # retained in the rich packets below.  Flattening only creates references
    # to the already hydrated records; it performs no physical rescaling.
    core_rows = _flat_core_rows(core_packet, w1_packet)
    collar_rows = _flat_collar_rows(collar_packet)
    owner_certificate = {
        "source_family_exact": True,
        "source_logs_exact": True,
        "rho_seam_only": True,
        "owner_identity_proved": False,
        "source_index_members_complete": True,
    }
    scale_certificate = {
        "core_already_physically_scaled": True,
        "collar_already_physically_scaled": True,
        "physical_scale_reapplied": False,
        "pressure_endpoint_added": False,
    }

    progress["actual_leaf_build_entered"] = True
    try:
        actual = leaf.build(
            context,
            core=core_rows,
            collar=collar_rows,
            source_family=source_family,
            source_atoms=source_logs,
            source_logs=source_logs,
            geometry_certificate=geometry,
            owner_certificate=owner_certificate,
            scale_certificate=scale_certificate,
        )
    finally:
        progress["actual_leaf_build_calls"] += 1
    require(isinstance(actual, Mapping), "segmented leaf actual result mapping")
    enriched = dict(actual)
    enriched["schema"] = COMPONENT
    _validate_leaf_result(enriched, namespace_indexes, source_family, source_logs, core_rows, collar_rows, progress)

    return {
        "source_family": encode(source_family),
        "source_logs": encode(source_logs),
        "geometry": geometry,
        "fixture": progress.get("fixture"),
        "core_w0_by_cell_channel": _compact_cell_packet(core_packet),
        "nested_w1_k0_by_cell_channel": _compact_cell_packet(w1_packet),
        "prior_collar_groups": _compact_groups({"collar": collar_packet["groups"]}).get("collar", {}),
        "prior_collar_index_bundle": {
            "input_index": encode(collar_packet["input_index"]),
            "prior_input_index": encode(collar_packet["prior_input_index"]),
            "group_input_index": encode(collar_packet["group_input_index"]),
        },
        "source_index_namespaces": {
            "core_w0": {
                "prior_index": encode(core_packet["source_index"]),
                "actual_ids": sorted(core_ids),
            },
            "core_w1": {
                "prior_index": encode(w1_packet["source_index"]),
                "actual_ids": sorted(w1_ids),
            },
            "collar": {
                "prior_index": encode(collar_packet["group_input_index"]),
                "actual_ids": sorted(collar_ids),
            },
        },
        "member_validation": progress.get("member_validation"),
        "result": _compact_leaf_result(actual, enriched, core_packet, w1_packet, collar_packet),
        "input_receipt_sha256": {CORE: RAW_SHA256[CORE], W1: RAW_SHA256[W1], COLLAR: RAW_SHA256[COLLAR]},
        "input_result_sha256": {CORE: RESULT_SHA256[CORE], W1: RESULT_SHA256[W1], COLLAR: RESULT_SHA256[COLLAR]},
        "published_source_linkage": {
            "old_collar_receipt": OLD_COLLAR,
            "old_collar_raw_sha256": OLD_COLLAR_SHA256,
            "core_manifest_old_collar_sha256": core_manifest.get(OLD_COLLAR),
            "collar_manifest_old_collar_sha256": collar_manifest.get(OLD_COLLAR),
            "nested_w1_manifest_core_sha256": w1_manifest.get(CORE),
        },
        "source_owner_identity_proved": False,
        "physical_scale_reapplied": False,
        "tail_interpreted_as_signed_field": False,
        "full_original_core_W1_missing": True,
        "nested_K0_W1_saved_not_full_W1": True,
        "Kchi_minus_K0_missing": True,
        "feedback_missing": True,
        "allword_tail_missing": True,
        "nested_W1_pressure_added": False,
        "Kchi_K0_included": False,
        "allword_tail_included": False,
        "endpoint_cancellation_proved": False,
        "terminal_closure_proved": False,
        "ledger_212_changed": False,
        "actual_leaf_build_calls": progress["actual_leaf_build_calls"],
        "fixture_leaf_build_calls": progress["fixture_leaf_calls"],
    }


def main() -> int:
    hashes, source_sha256, brief_sha256 = preflight()
    entry_sha256 = save(
        ENTRY,
        {
            "component": COMPONENT,
            "source_sha256": source_sha256,
            "brief_sha256": brief_sha256,
            "input_hashes": hashes,
            "entered": True,
            "raw_gzip_authenticated_before_decode": True,
            "exclusive_entry_and_result": True,
        },
    )
    progress: dict[str, Any] = {
        "raw_bytes_authenticated": False,
        "fixture_completed_before_input_json_decode": False,
        "input_json_decode_started": False,
        "core_receipt_decoded": False,
        "nested_w1_receipt_decoded": False,
        "collar_receipt_decoded": False,
        "fixture_leaf_calls": 0,
        "actual_leaf_build_entered": False,
        "actual_leaf_build_calls": 0,
        "leaf_result_validated": False,
    }
    outcome: dict[str, Any] = {
        "component": COMPONENT,
        "source_sha256": source_sha256,
        "entry_sha256": entry_sha256,
        "status": "failed",
        "raw_gzip_authenticated_before_decode": True,
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
    print(json.dumps({"component": COMPONENT, "status": outcome["status"], "entry_sha256": entry_sha256, "result_sha256": result_sha256, "error": (outcome.get("error") or {}).get("message")}))
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())

