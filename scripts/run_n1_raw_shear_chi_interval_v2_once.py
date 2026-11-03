"""One-shot v2 receipt for the raw lower-collar shear chi consumer.

This runner keeps the v1 leaf frozen and changes only the fixture/runner
boundary.  The complete v1 failed receipt is authenticated and retained as
evidence; the v1 runner and v1 fixture are never re-entered.  The source
family receipt is fully hashed and checked, but only its structural
``source_family`` field is passed to the v1 leaf.  This prevents unrelated
archival metadata (notably execution-time floats) from entering the leaf's
structural comparison.

Only graph nodes 35951 and 35954 perform new Fraction arithmetic.  Nodes
35949 and 35953 are consumed as the authenticated saved sigma and sigma_s
intervals.  No physical m15, Kchi, live-owner, or complete g5 claim is made.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.machinery
import importlib.util
import json
import math
import os
import sys
import traceback
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path
from types import ModuleType
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"

COMPONENT = "n1_raw_shear_chi_interval_v2"
PHASE = COMPONENT + "_pending"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_MODULE = "openai_ns_reconstruction.n1_raw_shear_chi_interval_v1"
FIXTURE_MODULE = "openai_ns_reconstruction.n1_raw_shear_chi_interval_v2_fixture"
LEAF_REL = "src/openai_ns_reconstruction/n1_raw_shear_chi_interval_v1.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_raw_shear_chi_interval_v2_fixture.py"
RUNNER_REL = "scripts/run_n1_raw_shear_chi_interval_v2_once.py"
V1_FAILURE_REL = "experiments/root_st073/n1_raw_shear_chi_interval_v1_failure.json.gz"
GRAPH_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz"
)
SIGMA_REL = "experiments/root_st073/n1_raw_sigma_scoped_oracle_v2_receipt.json.gz"
HB_REL = "experiments/root_st073/n1_collar_local_g5_j2_tail_v1_receipt.json.gz"
SOURCE_FAMILY_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions_check.json"
)

LEAF_SCHEMA = "n1_raw_shear_chi_interval_v1"
FIXTURE_SCHEMA = COMPONENT + "_fixture"
RESULT_SCHEMA = COMPONENT + "_result"

GRAPH_SHA256 = "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe"
SIGMA_SHA256 = "698f558d6ae5d179ea7700a48fc0af2d31d630be64fe1c24b639cb04e68d507a"
HB_SHA256 = "67894bac58f4fbbcd5437bf970dadafeb918fe2d9638425f8e5165d21af4530b"
SOURCE_FAMILY_SHA256 = "f3c4d1b230b3835da6f4f270926fa85f74ca488f8331e05fbbb894a7c22a0158"
LEAF_SHA256 = "71e0bc342776cff2a6b0cdd434d22060fe0c72f63ec1c96ad249dbaf52462028"
V1_FAILURE_SHA256 = "3df35ac8157b1cd1352eefb74e620764829e8817ade47bf7925c8891e2eb4af6"
HB_CERTIFICATE_SOURCE_SHA256 = "0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e"

FROZEN_INPUTS = (
    LEAF_REL,
    FIXTURE_REL,
    RUNNER_REL,
    V1_FAILURE_REL,
    GRAPH_REL,
    SIGMA_REL,
    HB_REL,
    SOURCE_FAMILY_REL,
)
FROZEN_SET = frozenset(FROZEN_INPUTS)
GZIP_INPUTS = frozenset((GRAPH_REL, SIGMA_REL, HB_REL, V1_FAILURE_REL))
EXPECTED_HASHES = {
    LEAF_REL: LEAF_SHA256,
    V1_FAILURE_REL: V1_FAILURE_SHA256,
    GRAPH_REL: GRAPH_SHA256,
    SIGMA_REL: SIGMA_SHA256,
    HB_REL: HB_SHA256,
    SOURCE_FAMILY_REL: SOURCE_FAMILY_SHA256,
}

CONSTANT_ONE_NODE = 1
CONSTANT_MINUS_ONE_NODE = 25085
HB_NODE = 8
HB_MINUS_ONE_NODE = 35926
SIGMA_PRODUCT_NODE = 35950
BINDER_NODE = 35848
SIGMA_NODE = 35949
SIGMA_S_NODE = 35953
CHI_NODE = 35951
CHI_S_NODE = 35954
PROGRAM_NODES = (
    CONSTANT_ONE_NODE,
    CONSTANT_MINUS_ONE_NODE,
    HB_NODE,
    HB_MINUS_ONE_NODE,
    SIGMA_PRODUCT_NODE,
    CHI_NODE,
    CHI_S_NODE,
)
SELECTED_NODES = PROGRAM_NODES + (SIGMA_NODE, SIGMA_S_NODE)
REAL_NODES = (CHI_NODE, CHI_S_NODE)
REUSED_NODES = (SIGMA_NODE, SIGMA_S_NODE)
SHA256_HEX = frozenset("0123456789abcdefABCDEF")
MISSING = object()


class RunnerInputError(ValueError):
    """Raised when a frozen v2 one-entry contract is not satisfied."""


def _require(condition: object, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def _sha_bytes(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _sha_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _json_safe(value: Any, seen: set[int] | None = None) -> Any:
    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, Fraction):
        return {"type": "Fraction", "numerator": value.numerator, "denominator": value.denominator}
    if isinstance(value, float):
        _require(math.isfinite(value), "non-finite float in receipt")
        return value
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, bytes):
        return {"type": "bytes", "length": len(value), "sha256": _sha_bytes(value)}
    if hasattr(value, "lower") and hasattr(value, "upper") and type(value).__name__ != "str":
        try:
            return {
                "type": type(value).__name__,
                "lower": _json_safe(value.lower, seen),
                "upper": _json_safe(value.upper, seen),
            }
        except Exception:
            pass
    active = seen if seen is not None else set()
    identity = id(value)
    if identity in active:
        return {"type": "cycle"}
    active.add(identity)
    try:
        if isinstance(value, Mapping):
            return {str(key): _json_safe(item, active) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [_json_safe(item, active) for item in value]
        if isinstance(value, (set, frozenset)):
            return sorted((_json_safe(item, active) for item in value), key=repr)
        if hasattr(value, "__dict__"):
            return {str(key): _json_safe(item, active) for key, item in vars(value).items()}
        return repr(value)
    finally:
        active.discard(identity)


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(_json_safe(value), ensure_ascii=False, sort_keys=True, separators=(",", ":"), allow_nan=False)
        .encode("utf-8")
        + b"\n"
    )


def _write_exclusive(path: Path, value: Any) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = _json_bytes(value)
    flags = os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_BINARY", 0)
    descriptor = os.open(str(path), flags, 0o600)
    try:
        with os.fdopen(descriptor, "wb") as stream:
            descriptor = -1
            stream.write(raw)
            stream.flush()
            os.fsync(stream.fileno())
    finally:
        if descriptor != -1:
            os.close(descriptor)
    return _sha_bytes(raw)


def _decode_object(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise RunnerInputError(label + " is not ordinary UTF-8 JSON") from exc
    _require(isinstance(value, Mapping), label + " is not a JSON object")
    return dict(value)


def _decode_gzip_object(raw: bytes, label: str) -> dict[str, Any]:
    try:
        decoded = gzip.decompress(raw)
    except (OSError, gzip.BadGzipFile) as exc:
        raise RunnerInputError(label + " is not authenticated gzip") from exc
    return _decode_object(decoded, label)


def _sha_text(value: Any, label: str) -> str:
    _require(isinstance(value, str) and len(value) == 64 and all(c in SHA256_HEX for c in value), "invalid SHA-256 for " + label)
    return value.lower()


def _repo_file(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    relative = Path(raw_path.replace("\\", "/"))
    _require(not relative.is_absolute() and ".." not in relative.parts, label + " path must be repository-relative")
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(path.is_relative_to(ROOT.resolve()) and path.is_file(), label + " path is missing or escapes repository")
    return key, path


def _manifest_hashes(manifest: Mapping[str, Any]) -> dict[str, str]:
    _require(manifest.get("schema") == COMPONENT + "_manifest", "manifest schema mismatch")
    _require(manifest.get("component") == COMPONENT and manifest.get("phase") == PHASE, "manifest component/phase mismatch")
    for key in ("one_entry_only", "do_not_reenter", "source_manifest_self_hash_forbidden"):
        _require(manifest.get(key) is True, "manifest guard missing: " + key)
    _require(manifest.get("runner_path") == RUNNER_REL, "manifest runner path mismatch")
    contract = manifest.get("input_contract")
    _require(isinstance(contract, Mapping), "manifest input_contract missing")
    expected = {
        "leaf": LEAF_REL,
        "leaf_schema": LEAF_SCHEMA,
        "fixture": FIXTURE_REL,
        "runner": RUNNER_REL,
        "v1_failure_evidence": V1_FAILURE_REL,
        "v1_failure_archive": True,
        "graph": GRAPH_REL,
        "sigma_receipt": SIGMA_REL,
        "hb_receipt": HB_REL,
        "source_family_receipt": SOURCE_FAMILY_REL,
        "fixture_api": "run_fixture(leaf)",
        "consumer_api": "compute(selected_rows, sigma_interval, sigma_s_interval, hb_evidence, graph_sha256, source_family, source_family_receipt, z_condition, rows_authenticated, sigma_bounds_authenticated)",
        "consumer_source_family_receipt": {"source_family": "source_family_receipt.source_family"},
        "reused_nodes": list(REUSED_NODES),
        "real_nodes": list(REAL_NODES),
        "fraction_only_new_nodes": True,
    }
    for key, value in expected.items():
        _require(contract.get(key) == value, "manifest input_contract." + key + " mismatch")
    scope = manifest.get("scope")
    _require(isinstance(scope, Mapping), "manifest scope missing")
    for key, value in {
        "real_nodes": list(REAL_NODES),
        "reused_nodes": list(REUSED_NODES),
        "only_new_nodes": True,
        "complete_g5": False,
        "physical_m15": False,
        "P1_is_m15": False,
        "Kchi": False,
        "live_owner": False,
        "source_identity_verified": False,
        "v1_failure_reentered": False,
    }.items():
        _require(scope.get(key) == value, "manifest scope mismatch: " + key)
    files = manifest.get("files")
    _require(isinstance(files, list) and len(files) == len(FROZEN_INPUTS), "manifest frozen file count mismatch")
    hashes: dict[str, str] = {}
    for index, item in enumerate(files):
        _require(isinstance(item, Mapping), f"manifest files[{index}] is not an object")
        key, path = _repo_file(item.get("path"), f"manifest files[{index}]")
        _require(key in FROZEN_SET and key not in hashes, "manifest file is out of scope or duplicated: " + key)
        expected_hash = _sha_text(item.get("sha256"), key)
        actual = _sha_file(path)
        _require(actual == expected_hash, "source hash mismatch: " + key)
        if key in GZIP_INPUTS:
            _require(item.get("raw_sha256_before_decode") is True, "gzip raw hash pin missing: " + key)
        hashes[key] = actual
    _require(frozenset(hashes) == FROZEN_SET, "manifest does not freeze exactly the required inputs")
    for key, expected_hash in EXPECTED_HASHES.items():
        _require(hashes[key] == expected_hash, "accepted raw SHA mismatch: " + key)
    return {key: hashes[key] for key in FROZEN_INPUTS}


def _brief_record(brief: Mapping[str, Any]) -> Mapping[str, Any]:
    if brief.get("component") == COMPONENT:
        return brief
    nested = brief.get(COMPONENT)
    if isinstance(nested, Mapping):
        return nested
    continuation = brief.get("continuation")
    if isinstance(continuation, Mapping) and continuation.get("component") == COMPONENT:
        return continuation
    return brief


def _validate_pending_brief(source_sha: str) -> str:
    _require(BRIEF.is_file(), "pending resume brief is missing")
    raw = BRIEF.read_bytes()
    record = _brief_record(_decode_object(raw, "pending resume brief"))
    _require(record.get("component") == COMPONENT and record.get("phase") == PHASE, "pending brief component/phase mismatch")
    _require(record.get("one_entry_only") is True and record.get("do_not_reenter") is True, "pending brief guard missing")
    checkpoint = record.get("source_checkpoint")
    if isinstance(checkpoint, Mapping):
        path_value, expected = checkpoint.get("path"), checkpoint.get("sha256")
    else:
        path_value, expected = checkpoint, None
    _require(isinstance(path_value, str), "pending brief source checkpoint path missing")
    path = Path(path_value.replace("\\", "/"))
    if not path.is_absolute():
        path = (ROOT / path).resolve()
    _require(path == SOURCE.resolve(), "pending brief points to another source checkpoint")
    _require(_sha_text(expected, "pending source checkpoint") == source_sha, "pending brief checkpoint SHA mismatch")
    return _sha_bytes(raw)


def _preflight() -> tuple[dict[str, str], str, str, dict[str, Any]]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered or has a result")
    _require(SOURCE.is_file(), "source manifest checkpoint is missing")
    source_raw = SOURCE.read_bytes()
    source_sha = _sha_bytes(source_raw)
    manifest = _decode_object(source_raw, "source manifest")
    brief_sha = _validate_pending_brief(source_sha)
    return _manifest_hashes(manifest), source_sha, brief_sha, manifest


def _read_authenticated(hashes: Mapping[str, str]) -> dict[str, bytes]:
    raw: dict[str, bytes] = {}
    for relative in FROZEN_INPUTS:
        data = (ROOT / relative).read_bytes()
        _require(_sha_bytes(data) == hashes[relative], "frozen input changed after entry: " + relative)
        raw[relative] = data
    for relative in GZIP_INPUTS:
        _require(_sha_bytes(raw[relative]) == EXPECTED_HASHES[relative], "gzip raw SHA changed before decode: " + relative)
    return raw


def _load_module(relative: str, module_name: str, allowed: set[str]) -> ModuleType:
    _require(relative in (LEAF_REL, FIXTURE_REL), "unexpected module path")
    before = set(sys.modules)
    package_name = "openai_ns_reconstruction"
    package_dir = ROOT / "src" / package_name
    previous_package = sys.modules.get(package_name)
    namespace = ModuleType(package_name)
    namespace.__path__ = [str(package_dir)]
    namespace.__package__ = package_name
    namespace.__spec__ = importlib.machinery.ModuleSpec(package_name, loader=None, is_package=True)
    sys.modules[package_name] = namespace
    old_path = list(sys.path)
    source_root = str(ROOT / "src")
    if source_root not in sys.path:
        sys.path.insert(0, source_root)
    try:
        spec = importlib.util.spec_from_file_location(module_name, ROOT / relative)
        _require(spec is not None and spec.loader is not None, "module loader missing: " + relative)
        module = importlib.util.module_from_spec(spec)
        sys.modules[module_name] = module
        spec.loader.exec_module(module)
    finally:
        sys.path[:] = old_path
        if previous_package is None:
            sys.modules.pop(package_name, None)
        else:
            sys.modules[package_name] = previous_package
    for name in set(sys.modules) - before:
        if name == package_name or name in allowed:
            continue
        if name.startswith(package_name + "."):
            raise RunnerInputError("unexpected project module imported: " + name)
    return module


def _validate_fixture(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "fixture result is not an object")
    result = dict(value)
    _require(result.get("schema") == FIXTURE_SCHEMA, "v2 fixture schema mismatch")
    _require(result.get("passed") is True, "v2 fixture did not pass")
    checks = result.get("checks")
    _require(isinstance(checks, Mapping) and checks and all(item is True for item in checks.values()), "v2 fixture checks failed")
    _require(result.get("synthetic_only") is True, "v2 fixture is not synthetic-only")
    for key in ("used_real_graph", "used_prior_interval", "used_saved_result", "used_old_fixture", "used_old_module", "old_fixture_called", "old_producer_called", "historical_modules_called"):
        if key in result:
            _require(result.get(key) is False, "v2 fixture consumed forbidden input: " + key)
    return result


def _fraction(value: Any, label: str) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if type(value) is int:
        return Fraction(value)
    _require(isinstance(value, Mapping), label + " is not an exact Fraction record")
    if "type" in value:
        _require(value.get("type") == "Fraction", label + " is not a Fraction record")
    numerator, denominator = value.get("numerator"), value.get("denominator")
    _require(type(numerator) is int and type(denominator) is int and denominator > 0, label + " is not exact integer data")
    return Fraction(numerator, denominator)


def _fractionize(value: Any, label: str) -> Any:
    if isinstance(value, Mapping):
        if value.get("type") == "Fraction":
            return _fraction(value, label)
        return {key: _fractionize(item, label + "." + str(key)) for key, item in value.items()}
    if isinstance(value, list):
        return [_fractionize(item, label + "[]") for item in value]
    if isinstance(value, tuple):
        return tuple(_fractionize(item, label + "[]") for item in value)
    if isinstance(value, float):
        raise RunnerInputError(label + " contains a float; exact Fraction input required")
    return value


def _interval(value: Any, label: str) -> tuple[Fraction, Fraction]:
    if isinstance(value, Mapping):
        nested = value.get("interval", MISSING)
        if nested is not MISSING:
            return _interval(nested, label + ".interval")
        lower, upper = value.get("lower", MISSING), value.get("upper", MISSING)
    else:
        _require(isinstance(value, Sequence) and not isinstance(value, (str, bytes)) and len(value) == 2, label + " is not an interval")
        lower, upper = value[0], value[1]
    lo, hi = _fraction(lower, label + ".lower"), _fraction(upper, label + ".upper")
    _require(lo <= hi, label + " lower exceeds upper")
    return lo, hi


def _named_interval(mapping: Mapping[str, Any], names: tuple[str, ...], label: str) -> tuple[Fraction, Fraction]:
    found = [(name, _interval(mapping[name], label + "." + name)) for name in names if name in mapping]
    _require(found, label + " is missing")
    first_name, first = found[0]
    for name, value in found[1:]:
        _require(value == first, label + " aliases disagree: " + first_name + "/" + name)
    return first


def _saved_sigma(receipt: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any], tuple[Fraction, Fraction], tuple[Fraction, Fraction]]:
    candidate = receipt.get("result")
    _require(isinstance(candidate, Mapping), "saved sigma receipt result missing")
    _require(candidate.get("status") == "returned", "saved sigma receipt is not returned")
    _require(candidate.get("real_roots") == [SIGMA_NODE, SIGMA_S_NODE], "saved sigma real-root scope mismatch")
    _require(candidate.get("real_root_evaluation_calls") == 2, "saved sigma root evaluation count mismatch")
    sidecar = candidate.get("graph_sidecar")
    _require(isinstance(sidecar, Mapping) and sidecar.get("graph_sha256") == GRAPH_SHA256, "saved sigma graph sidecar SHA mismatch")
    direct = candidate.get("prior_direct_interval")
    _require(isinstance(direct, Mapping), "saved sigma prior_direct_interval missing")
    _require(direct.get("scope") == "direct lower-collar chi sigma(s) only", "saved sigma scope mismatch")
    _require(direct.get("fixed_s_in_z") is True and direct.get("fixed_s_over_sc_in_z") is True, "saved sigma Z premise missing")
    _require(direct.get("sigma_upper_is_not_a_value") is True and direct.get("sigma_value") is None, "saved sigma cap was materialized")
    sigma = _named_interval(direct, ("sigma", "sigma_interval"), "saved sigma interval")
    sigma_s = _named_interval(direct, ("d_sigma_ds", "sigma_s_interval", "derivative"), "saved sigma_s interval")
    _require(_interval(candidate.get("prior_sigma_interval"), "prior_sigma_interval") == sigma, "prior sigma aliases disagree")
    _require(_interval(candidate.get("prior_sigma_s_interval"), "prior_sigma_s_interval") == sigma_s, "prior sigma_s aliases disagree")
    real_results = candidate.get("real_results")
    _require(isinstance(real_results, list) and len(real_results) == 2, "saved sigma real_results are missing")
    seen: set[int] = set()
    for index, row in enumerate(real_results):
        _require(isinstance(row, Mapping), f"saved sigma real_results[{index}] is not an object")
        root_id = row.get("root")
        _require(root_id in REUSED_NODES and root_id not in seen, f"saved sigma real_results[{index}] root mismatch")
        seen.add(root_id)
        value = row.get("value")
        _require(isinstance(value, Mapping) and value.get("node_id") == root_id, f"saved sigma real_results[{index}] value mismatch")
        _require(value.get("graph_sha256") == GRAPH_SHA256, f"saved sigma real root {root_id} graph SHA mismatch")
        expected = sigma if root_id == SIGMA_NODE else sigma_s
        _require(_interval(value, f"saved sigma root {root_id}") == expected, f"saved sigma root {root_id} interval mismatch")
    _require(seen == set(REUSED_NODES), "saved sigma roots 35949/35953 were not both returned")
    return dict(candidate), dict(direct), sigma, sigma_s


def _op(row: Mapping[str, Any]) -> Any:
    operation, shorthand = row.get("operation", MISSING), row.get("op", MISSING)
    if operation is not MISSING and shorthand is not MISSING:
        _require(operation == shorthand, "graph operation/op disagreement")
    return operation if operation is not MISSING else shorthand


def _arguments(row: Mapping[str, Any]) -> list[Any]:
    value = row.get("arguments", row.get("operands", row.get("argument", MISSING)))
    if "argument" in row and not isinstance(value, Sequence):
        value = [value]
    _require(isinstance(value, Sequence) and not isinstance(value, (str, bytes)), "graph arguments are missing")
    return list(value)


def _validate_graph(graph: Mapping[str, Any]) -> dict[str, Any]:
    nodes = graph.get("exact_graph_nodes")
    _require(isinstance(nodes, list) and len(nodes) > max(SELECTED_NODES), "graph selected rows are missing")
    selected: dict[str, dict[str, Any]] = {}
    for node_id in SELECTED_NODES:
        row = nodes[node_id]
        _require(isinstance(row, Mapping), f"graph row {node_id} is not an object")
        selected[str(node_id)] = dict(row)
    one = selected[str(CONSTANT_ONE_NODE)]
    minus_one = selected[str(CONSTANT_MINUS_ONE_NODE)]
    hb = selected[str(HB_NODE)]
    _require(_op(one) == "exact_rational", "node 1 operation mismatch")
    _require(_fraction(one.get("value", one.get("constant", one)), "node 1 value") == Fraction(1), "node 1 is not exact one")
    _require(_op(minus_one) == "negative", "node 25085 must be negative")
    minus_args = _arguments(minus_one)
    _require(len(minus_args) == 1 and minus_args[0] == CONSTANT_ONE_NODE, "node 25085 must negate node 1")
    _require(_op(hb) in {"current_original_source_parameter", "original_source_parameter", "source_parameter"}, "node 8 operation mismatch")
    _require(hb.get("name", hb.get("parameter", hb.get("symbol"))) == "hbB", "node 8 hbB name mismatch")
    expected = {
        HB_MINUS_ONE_NODE: ("sum", [HB_NODE, CONSTANT_MINUS_ONE_NODE]),
        SIGMA_PRODUCT_NODE: ("product", [HB_MINUS_ONE_NODE, SIGMA_NODE]),
        CHI_NODE: ("sum", [CONSTANT_ONE_NODE, SIGMA_PRODUCT_NODE]),
        CHI_S_NODE: ("product", [HB_MINUS_ONE_NODE, SIGMA_S_NODE]),
    }
    for node_id, (operation, args) in expected.items():
        _require(_op(selected[str(node_id)]) == operation, f"node {node_id} operation mismatch")
        _require(_arguments(selected[str(node_id)]) == args, f"node {node_id} arguments mismatch")
    family = hb.get("source_family", MISSING)
    if family is MISSING and isinstance(hb.get("source"), Mapping):
        family = hb["source"].get("source_family", hb["source"].get("family", MISSING))
    _require(isinstance(family, Mapping), "graph node 8 source_family missing")
    current_receipt = hb.get("current_source_receipt", hb.get("source_receipt", MISSING))
    if current_receipt is MISSING and isinstance(hb.get("source"), Mapping):
        current_receipt = hb["source"].get("current_source_receipt", hb["source"].get("source_receipt", MISSING))
    _require(current_receipt == SOURCE_FAMILY_REL or current_receipt == Path(SOURCE_FAMILY_REL).name, "node 8 source receipt mismatch")
    return {
        "node_count": len(nodes),
        "selected_nodes": selected,
        "node8_source_family": dict(family),
        "node8_current_source_receipt": current_receipt,
        "source_identity_verified": False,
        "live_owner_binding": False,
    }


def _locate_hb(receipt: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any], str]:
    result = receipt.get("result")
    _require(isinstance(result, Mapping), "hb receipt result missing")
    output = result.get("output")
    _require(isinstance(output, Mapping), "hb receipt result.output missing")
    actual = output.get("actual_local_j2_tail")
    _require(isinstance(actual, Mapping), "hb receipt actual_local_j2_tail missing")
    geometry, certificate, family = actual.get("geometry"), actual.get("certificate_passthrough"), output.get("source_family")
    _require(isinstance(geometry, Mapping), "hb geometry missing")
    _require(isinstance(certificate, Mapping), "hb certificate_passthrough missing")
    _require(isinstance(family, Mapping), "hb output source_family missing")
    return dict(geometry), dict(certificate), dict(family), "result.output.actual_local_j2_tail.geometry"


def _validate_hb(geometry: Mapping[str, Any], certificate: Mapping[str, Any], path: str) -> dict[str, Any]:
    _require(geometry.get("strict_hb_admitted") is True, "strict hb certificate missing")
    _require(geometry.get("hb_sc_Z_independent") is True, "hb Z-independence certificate missing")
    interval = geometry.get("hb_interval")
    _require(isinstance(interval, Mapping), "hb_interval missing")
    _require(interval.get("lower_exact_mpf_tuple") == [0, 0, 0, 0], "hb lower tuple mismatch")
    _require(interval.get("upper_exact_mpf_tuple") == [0, 1, -1, 1], "hb upper tuple mismatch")
    _require(certificate.get("source_sha256") == HB_CERTIFICATE_SOURCE_SHA256, "hb certificate source SHA mismatch")
    return {
        "geometry_path": path,
        "strict_hb_admitted": True,
        "hb_Z_independent": True,
        "condition": "0<hb<1/2",
        "source_sha256": HB_CERTIFICATE_SOURCE_SHA256,
        "hb_interval_encoding": {
            "lower_exact_mpf_tuple": [0, 0, 0, 0],
            "upper_exact_mpf_tuple": [0, 1, -1, 1],
        },
        "interval_bounds": (Fraction(0), Fraction(1, 2)),
        "source_identity_verified": False,
        "live_owner_binding": False,
    }


def _validate_source_family(graph_scope: Mapping[str, Any], source_receipt: Mapping[str, Any], hb_family: Mapping[str, Any]) -> dict[str, Any]:
    graph_family = graph_scope["node8_source_family"]
    receipt_family = source_receipt.get("source_family")
    _require(isinstance(receipt_family, Mapping), "source receipt source_family missing")
    _require(_json_bytes(graph_family) == _json_bytes(receipt_family), "graph/source receipt source_family mismatch")
    _require(_json_bytes(graph_family) == _json_bytes(hb_family), "graph/hb source_family mismatch")
    return {
        "same_source_family": True,
        "graph_node8_source_family": graph_family,
        "source_receipt_source_family": receipt_family,
        "hb_geometry_source_family": hb_family,
        "source_receipt_sha256": SOURCE_FAMILY_SHA256,
        "source_identity_verified": False,
        "live_owner_binding": False,
    }


def _validate_v1_failure(value: Mapping[str, Any]) -> dict[str, Any]:
    _require(value.get("component") == "n1_raw_shear_chi_interval_v1", "v1 failure component mismatch")
    _require(value.get("status") == "failed", "v1 evidence is not a failed receipt")
    error = value.get("error")
    _require(isinstance(error, Mapping), "v1 failure error missing")
    _require(error.get("type") == "GraphRecordError", "v1 failure type mismatch")
    _require("execution_seconds" in str(error.get("message")), "v1 failure does not preserve pre-math float evidence")
    progress = value.get("progress")
    _require(isinstance(progress, Mapping) and progress.get("consumer_calls") == 1, "v1 failure consumer entry count missing")
    _require("consumer_result" not in value, "v1 failure unexpectedly contains a math result")
    return {
        "path": V1_FAILURE_REL,
        "sha256": V1_FAILURE_SHA256,
        "status": "failed",
        "error_type": error.get("type"),
        "error_message": error.get("message"),
        "consumer_calls": progress.get("consumer_calls"),
        "math_result_present": False,
        "reentered": False,
        "receipt": dict(value),
    }


def _product(a: tuple[Fraction, Fraction], b: tuple[Fraction, Fraction]) -> tuple[Fraction, Fraction]:
    values = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return min(values), max(values)


def _node_interval(result: Mapping[str, Any], aliases: tuple[str, ...], node_id: int) -> tuple[Fraction, Fraction]:
    for key in aliases:
        if key in result:
            return _interval(result[key], "consumer node " + str(node_id))
    node_results = result.get("node_results")
    if isinstance(node_results, Mapping):
        row = node_results.get(str(node_id), node_results.get(node_id))
        if isinstance(row, Mapping):
            return _interval(row.get("interval", row), "consumer node " + str(node_id))
    raise RunnerInputError("consumer interval missing for node " + str(node_id))


def _validate_consumer(value: Any, hb_scope: Mapping[str, Any], sigma: tuple[Fraction, Fraction], sigma_s: tuple[Fraction, Fraction]) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "consumer result is not an object")
    output = dict(value)
    _require(output.get("schema") == LEAF_SCHEMA + "_result", "consumer schema mismatch")
    _require(output.get("status") == "conditional_raw_g5_lower_collar_shear_chi_interval", "consumer status mismatch")
    _require(output.get("coordinate") == "raw_s", "consumer coordinate mismatch")
    for key, expected in {
        "s_derivative_is_Z_derivative": False,
        "s_derivative_is_q_derivative": False,
        "z_derivative": False,
        "q_derivative": False,
        "chi_s_is_Z_derivative": False,
        "chi_s_is_q_derivative": False,
        "ordinary_s_derivative": True,
        "fraction_only_new_nodes": True,
        "reused_sigma_intervals": True,
        "hb_Z_independent": True,
        "source_family_receipt_checked": True,
        "not_physical_m15": True,
        "not_Kchi_axial_potential": True,
    }.items():
        _require(output.get(key) is expected, "consumer flag mismatch: " + key)
    for key in ("physical_m15", "full_g5", "source_identity_verified", "owner_binding", "live_owner_verified", "sigma_recomputed", "sigma_s_recomputed", "prior_sigma_intervals_reused"):
        if key in output and key not in {"prior_sigma_intervals_reused"}:
            _require(output.get(key) is False, "consumer forbidden claim: " + key)
    _require(output.get("prior_sigma_intervals_reused") is True, "consumer prior sigma reuse missing")
    selected = output.get("selected_node_ids")
    _require(isinstance(selected, Sequence) and not isinstance(selected, (str, bytes)) and set(selected) == set(PROGRAM_NODES), "consumer node scope mismatch")
    chi = _node_interval(output, ("chi_interval", "chi"), CHI_NODE)
    chi_s = _node_interval(output, ("chi_s_interval", "chi_s"), CHI_S_NODE)
    hb = hb_scope["interval_bounds"]
    expected_chi_product = _product((hb[0] - 1, hb[1] - 1), sigma)
    expected_chi = (Fraction(1) + expected_chi_product[0], Fraction(1) + expected_chi_product[1])
    expected_chi_s = _product((hb[0] - 1, hb[1] - 1), sigma_s)
    _require(chi == expected_chi, "node 35951 Fraction interval mismatch")
    _require(chi_s == expected_chi_s, "node 35954 Fraction interval mismatch")
    return {
        "result": output,
        "chi_35951": chi,
        "chi_s_35954": chi_s,
        "fraction_endpoints": True,
        "arithmetic_checked": True,
        "scope": {
            "new_nodes": list(REAL_NODES),
            "reused_nodes": list(REUSED_NODES),
            "complete_g5": False,
            "physical_m15": False,
            "Kchi": False,
            "live_owner": False,
            "source_identity_verified": False,
        },
    }


def _consume(hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    allowed = {"openai_ns_reconstruction", LEAF_MODULE, FIXTURE_MODULE}
    progress["stage"] = "module_import"
    try:
        leaf = _load_module(LEAF_REL, LEAF_MODULE, allowed)
        fixture = _load_module(FIXTURE_REL, FIXTURE_MODULE, allowed)
    except BaseException as exc:
        progress["module_import_exception"] = {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}
        raise
    progress["modules_imported_after_entry"] = True
    _require(getattr(leaf, "SCHEMA", None) == LEAF_SCHEMA, "frozen v1 leaf schema mismatch")
    _require(getattr(fixture, "SCHEMA", None) == FIXTURE_SCHEMA, "v2 fixture schema mismatch")

    fixture_runner = getattr(fixture, "run_fixture", None)
    _require(callable(fixture_runner), "v2 fixture has no run_fixture(leaf) API")
    progress["fixture_calls"] = 1
    progress["stage"] = "fixture_execution"
    try:
        fixture_raw = fixture_runner(leaf)
    except BaseException as exc:
        progress["fixture_exception"] = {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}
        raise
    progress["fixture_result_raw"] = _json_safe(fixture_raw)
    progress["stage"] = "fixture_validation"
    try:
        fixture_result = _validate_fixture(fixture_raw)
    except BaseException as exc:
        progress["fixture_validation_exception"] = {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}
        raise
    progress["fixture"] = fixture_result

    progress["stage"] = "authenticated_raw_read"
    raw = _read_authenticated(hashes)
    progress["raw_input_lengths"] = {path: len(data) for path, data in raw.items()}
    progress["stage"] = "gzip_decode"
    graph = _decode_gzip_object(raw[GRAPH_REL], GRAPH_REL)
    sigma_receipt = _decode_gzip_object(raw[SIGMA_REL], SIGMA_REL)
    hb_receipt = _decode_gzip_object(raw[HB_REL], HB_REL)
    source_family_receipt = _decode_object(raw[SOURCE_FAMILY_REL], SOURCE_FAMILY_REL)
    v1_failure_archive = _decode_gzip_object(raw[V1_FAILURE_REL], V1_FAILURE_REL)
    v1_failure = v1_failure_archive.get("result")
    _require(isinstance(v1_failure, Mapping), "v1 failure archive top-level result missing")
    progress["raw_gzip_authenticated_before_decode"] = True
    progress["gzip_decoded_after_fixture"] = True

    progress["stage"] = "v1_failure_validation"
    failure_scope = _validate_v1_failure(v1_failure)
    progress["v1_failure_scope"] = {key: value for key, value in failure_scope.items() if key != "receipt"}
    progress["stage"] = "graph_validation"
    graph_scope = _validate_graph(graph)
    progress["stage"] = "hb_validation"
    hb_geometry, hb_certificate, hb_family, hb_path = _locate_hb(hb_receipt)
    hb_scope = _validate_hb(hb_geometry, hb_certificate, hb_path)
    progress["stage"] = "source_family_validation"
    family_scope = _validate_source_family(graph_scope, source_family_receipt, hb_family)
    progress["stage"] = "saved_sigma_validation"
    prior_result, prior_direct, sigma, sigma_s = _saved_sigma(sigma_receipt)
    saved_direct = _fractionize(prior_direct, "prior_direct_interval")

    # Preserve the complete receipt hash/family check, while passing only the
    # structural field to the leaf so execution_seconds and similar metadata
    # cannot enter its exact structural comparator.
    reduced_family_receipt = {"source_family": source_family_receipt["source_family"]}
    hb_evidence = {
        "strict_hb_admitted": True,
        "hb_sc_Z_independent": True,
        "hb_interval": {"lower": Fraction(0), "upper": Fraction(1, 2)},
        "source_family": hb_family,
        "fixed_s": True,
        "fixed_hb": True,
        "fixed_s_hb": True,
        "ordinary_Z1": Fraction(0),
        "s_derivative_not_Z_derivative": True,
        "s_derivative_is_q_derivative": False,
        "source_certificate_sha256": HB_CERTIFICATE_SOURCE_SHA256,
    }
    consumer = getattr(leaf, "compute", None)
    _require(callable(consumer), "frozen v1 leaf compute API missing")
    selected_rows = {int(key): value for key, value in graph_scope["selected_nodes"].items()}
    progress["stage"] = "consumer_call"
    progress["consumer_calls"] = 1
    try:
        consumer_raw = consumer(
            selected_rows,
            {"lower": sigma[0], "upper": sigma[1]},
            {"lower": sigma_s[0], "upper": sigma_s[1]},
            hb_evidence,
            graph_sha256=GRAPH_SHA256,
            source_family=graph_scope["node8_source_family"],
            source_family_receipt=reduced_family_receipt,
            z_condition={
                "fixed_s": True,
                "fixed_hb": True,
                "fixed_s_hb": True,
                "ordinary_Z1": Fraction(0),
                "s_derivative_not_Z_derivative": True,
                "s_derivative_is_q_derivative": False,
            },
            rows_authenticated=True,
            sigma_bounds_authenticated=True,
        )
    except BaseException as exc:
        progress["consumer_exception"] = {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()}
        raise
    progress["consumer_result_raw"] = _json_safe(consumer_raw)
    progress["stage"] = "consumer_validation"
    consumer_scope = _validate_consumer(consumer_raw, hb_scope, sigma, sigma_s)
    progress["consumer_result"] = consumer_scope
    return {
        "fixture": fixture_result,
        "v1_failure_evidence": failure_scope,
        "graph_scope": graph_scope,
        "hb_certificate_scope": hb_scope,
        "hb_geometry_path": hb_path,
        "source_family_scope": family_scope,
        "source_family_receipt_sha256": SOURCE_FAMILY_SHA256,
        "source_family_receipt_full_sha256": SOURCE_FAMILY_SHA256,
        "source_family_receipt_full_authenticated": True,
        "source_family_receipt_projection_keys": ["source_family"],
        "source_family_receipt_passed_to_consumer": reduced_family_receipt,
        "source_family_receipt_reduced_for_structural_leaf": True,
        "prior_result_status": prior_result.get("status"),
        "prior_direct_interval": prior_direct,
        "prior_direct_interval_fraction_adapter": saved_direct,
        "prior_sigma_interval": {"lower": sigma[0], "upper": sigma[1]},
        "prior_sigma_s_interval": {"lower": sigma_s[0], "upper": sigma_s[1]},
        "consumer_api": "compute",
        "consumer_scope": consumer_scope,
        "consumer_calls": 1,
        "fraction_only_new_nodes": True,
        "ancestor_reexecuted": False,
        "v1_failure_reentered": False,
        "old_fixture_called": False,
        "old_producer_called": False,
        "complete_g5": False,
        "physical_m15": False,
        "P1_is_m15": False,
        "Kchi": False,
        "source_identity_verified": False,
        "live_owner_binding": False,
    }


def main() -> int:
    try:
        hashes, source_sha, brief_sha, manifest = _preflight()
    except BaseException:
        print("status=preflight_failed")
        return 1
    try:
        entry_sha = _write_exclusive(
            ENTRY,
            {
                "schema": COMPONENT + "_entry",
                "component": COMPONENT,
                "phase": PHASE,
                "manifest_sha256": source_sha,
                "brief_sha256": brief_sha,
                "input_hashes": hashes,
                "raw_gzip_authenticated_before_decode": True,
                "entry_written_before_import": True,
                "fixture_before_decode": True,
                "v1_failure_evidence_reentered": False,
                "entered": True,
            },
        )
    except BaseException:
        print("status=entry_failed")
        return 1
    progress: dict[str, Any] = {
        "stage": "entered",
        "modules_imported_after_entry": False,
        "fixture_calls": 0,
        "consumer_calls": 0,
        "raw_gzip_authenticated_before_decode": False,
        "gzip_decoded_after_fixture": False,
    }
    outcome: dict[str, Any] = {
        "schema": RESULT_SCHEMA,
        "component": COMPONENT,
        "phase": PHASE,
        "manifest_sha256": source_sha,
        "brief_sha256": brief_sha,
        "entry_sha256": entry_sha,
        "input_hashes": hashes,
        "status": "failed",
        "progress": progress,
        "error": None,
        "traceback": None,
        "source_manifest_self_hash_forbidden": True,
        "raw_gzip_authenticated_before_decode": True,
        "complete_g5": False,
        "physical_m15": False,
        "P1_is_m15": False,
        "Kchi": False,
        "source_identity_verified": False,
        "live_owner_binding": False,
        "v1_failure_reentered": False,
    }
    try:
        outcome.update(_consume(hashes, progress))
        outcome["status"] = "returned"
    except BaseException as exc:
        outcome["error"] = {"type": type(exc).__name__, "message": str(exc)}
        outcome["traceback"] = traceback.format_exc()
    try:
        _write_exclusive(RESULT, outcome)
    except BaseException:
        print("status=result_write_failed")
        return 1
    print("status=" + str(outcome["status"]))
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
