"""Run the raw sigma scoped-oracle boundary exactly once.

This runner authenticates the new oracle, its fixture, the existing scoped
evaluator, and the two saved gzip inputs before any archive is decoded. The
exclusive entry is written before imports. The fixture is run once, then the
authenticated graph and prior direct-interval receipt are decoded, and only
roots 35949 and 35953 are evaluated once each.

The roots are conditional graph rows. The oracle returns the saved direct
lower-collar sigma interval; it does not integrate, rebuild, or import an old
producer. The result never claims source identity, a live owner, physical
m 15, or a complete g 5 evaluation. The source manifest hashes this runner
as a frozen input; it does not hash the manifest checkpoint itself.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import math
import os
import sys
import traceback
from collections.abc import Mapping
from fractions import Fraction
from pathlib import Path
from types import ModuleType
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"

COMPONENT = "n1_raw_sigma_scoped_oracle_v2"
PHASE = COMPONENT + "_pending"
CHECKPOINT_DIR = ROOT / "tmp" / COMPONENT
SOURCE = CHECKPOINT_DIR / "source.json"
ENTRY = CHECKPOINT_DIR / "entered.json"
RESULT = CHECKPOINT_DIR / "result.json"

LEAF_REL = "src/openai_ns_reconstruction/n1_raw_sigma_scoped_oracle_v1.py"
FIXTURE_REL = "src/openai_ns_reconstruction/n1_raw_sigma_scoped_oracle_v2_fixture.py"
EVALUATOR_REL = "src/openai_ns_reconstruction/n1_scoped_graph_evaluator_v1.py"
RUNNER_REL = "scripts/run_n1_raw_sigma_scoped_oracle_v2_once.py"
GRAPH_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz"
)
PRIOR_REL = "experiments/root_st073/n1_collar_sigma_actual_domain_v1_receipt.json.gz"

SIGMA_SOURCE_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_outer_initial.py"
)
SIGMA_PRIME_SOURCE_REL = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_flat_pulse_derivatives.py"
)

# The runner is included in the frozen source set, but the source manifest
# itself is not. This is a direct runner hash pin, not a self-hash constant in
# the runner, so it does not form a self-referential digest.
FROZEN_INPUTS = (
    LEAF_REL,
    FIXTURE_REL,
    EVALUATOR_REL,
    RUNNER_REL,
    GRAPH_REL,
    PRIOR_REL,
    SIGMA_SOURCE_REL,
    SIGMA_PRIME_SOURCE_REL,
)
FROZEN_SET = frozenset(FROZEN_INPUTS)
GZIP_INPUTS = frozenset({GRAPH_REL, PRIOR_REL})

GRAPH_SHA256 = "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe"
PRIOR_SHA256 = "b42ad7a699de3d9bc501b6f8f9d8040919079cafd214ce7bb9274c639b90d209"
SIGMA_SOURCE_SHA256 = "a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657"
SIGMA_PRIME_SOURCE_SHA256 = "32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7"

MANIFEST_SCHEMA = COMPONENT + "_manifest"
LEAF_SCHEMA = "n1_raw_sigma_scoped_oracle_v1"
FIXTURE_SCHEMA = COMPONENT + "_fixture"
EVALUATOR_SCHEMA = "n1_scoped_graph_evaluator_v1"
EVALUATOR_RESULT_SCHEMA = "n1_scoped_graph_evaluator_result_v1"

ARGUMENT_NODE_ID = 35848
SIGMA_NODE_ID = 35949
SIGMA_PRIME_NODE_ID = 35953
ROOTS = (SIGMA_NODE_ID, SIGMA_PRIME_NODE_ID)
EXPECTED_BINDER_INTERVAL = (Fraction(0), Fraction(1, 2**128))

SHA256_HEX = frozenset("0123456789abcdefABCDEF")
_MISSING = object()


class RunnerInputError(ValueError):
    """Raised for a fail-closed manifest, provenance, or scope violation."""


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


class _OpaqueScopedBackend:
    """Validate exact values while forbidding extra graph arithmetic."""

    def __init__(self, interval_types: tuple[type[Any], ...]) -> None:
        self.interval_types = interval_types

    def validate_scalar(self, value: Any) -> bool:
        return isinstance(value, (Fraction, *self.interval_types))

    def _forbidden(self, *_args: Any, **_kwargs: Any) -> Any:
        raise RunnerInputError("unexpected non-target graph arithmetic")

    rational = _forbidden
    add = _forbidden
    mul = _forbidden
    divide = _forbidden
    positive_divide = _forbidden
    unary = _forbidden
    integral = _forbidden

    def is_positive(self, _value: Any) -> bool:
        raise RunnerInputError("unexpected positivity query outside scoped roots")


def _json_safe(value: Any, seen: set[int] | None = None) -> Any:
    """Serialize exact Fraction interval objects without float coercion."""

    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, Fraction):
        return {
            "type": "Fraction",
            "numerator": value.numerator,
            "denominator": value.denominator,
        }
    if type(value).__name__ == "FractionInterval" and hasattr(value, "lower") and hasattr(value, "upper"):
        return {
            "type": "FractionInterval",
            "lower": _json_safe(value.lower),
            "upper": _json_safe(value.upper),
        }
    if type(value).__name__ == "ScopedSigmaInterval" and hasattr(value, "lower") and hasattr(value, "upper"):
        fields = (
            "lower",
            "upper",
            "node_id",
            "argument_node_id",
            "input_interval",
            "quantity",
            "graph_sha256",
            "source_result_sha256",
            "receipt_sha256",
            "certified",
            "conditional",
            "source_identity_verified",
            "owner_binding",
            "normalized",
            "full_g5",
            "physical_m15",
            "coordinate",
            "source_form",
        )
        return {
            "type": "ScopedSigmaInterval",
            **{
                field: _json_safe(getattr(value, field))
                for field in fields
                if hasattr(value, field)
            },
        }
    if isinstance(value, float):
        return value if math.isfinite(value) else repr(value)
    if isinstance(value, Path):
        return value.as_posix()
    if isinstance(value, bytes):
        return {"type": "opaque_bytes", "length": len(value), "sha256": _sha_bytes(value)}

    active = seen if seen is not None else set()
    identity = id(value)
    if identity in active:
        return {"type": "opaque_cycle", "repr": repr(value)}
    active.add(identity)
    try:
        if isinstance(value, Mapping):
            return {str(key): _json_safe(item, active) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [_json_safe(item, active) for item in value]
        if isinstance(value, (set, frozenset)):
            return sorted((_json_safe(item, active) for item in value), key=repr)
        return {"type": "opaque_value", "class": type(value).__name__, "repr": repr(value)}
    finally:
        active.discard(identity)


def _json_bytes(value: Any) -> bytes:
    return (
        json.dumps(
            _json_safe(value),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")


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
    _require(
        isinstance(value, str)
        and len(value) == 64
        and all(character in SHA256_HEX for character in value),
        "invalid SHA-256 for " + label,
    )
    return value.lower()


def _repo_file(raw_path: Any, label: str) -> tuple[str, Path]:
    _require(isinstance(raw_path, str) and raw_path, label + " path is missing")
    relative = Path(raw_path.replace("\\", "/"))
    _require(
        not relative.is_absolute() and ".." not in relative.parts,
        label + " path must be repository-relative",
    )
    key = relative.as_posix()
    path = (ROOT / relative).resolve()
    _require(
        path.is_relative_to(ROOT.resolve()) and path.is_file(),
        label + " path is missing or escapes repository: " + key,
    )
    return key, path


def _manifest_hashes(manifest: Mapping[str, Any]) -> dict[str, str]:
    _require(manifest.get("schema") == MANIFEST_SCHEMA, "manifest schema mismatch")
    _require(manifest.get("component") == COMPONENT, "manifest component mismatch")
    _require(manifest.get("phase") == PHASE, "manifest phase mismatch")
    _require(manifest.get("one_entry_only") is True, "manifest one-entry guard missing")
    _require(manifest.get("do_not_reenter") is True, "manifest re-entry guard missing")
    _require(
        manifest.get("source_manifest_self_hash_forbidden") is True,
        "manifest must forbid self-hashing",
    )
    _require(manifest.get("runner_path") == RUNNER_REL, "manifest runner path mismatch")

    contract = manifest.get("input_contract")
    _require(isinstance(contract, Mapping), "manifest input_contract missing")
    expected_contract = {
        "leaf": LEAF_REL,
        "fixture": FIXTURE_REL,
        "scoped_evaluator": EVALUATOR_REL,
        "runner": RUNNER_REL,
        "graph": GRAPH_REL,
        "prior_receipt": PRIOR_REL,
        "sigma_source": SIGMA_SOURCE_REL,
        "sigma_prime_source": SIGMA_PRIME_SOURCE_REL,
        "fixture_api": "run_fixture(leaf)",
        "oracle_api": "build_oracle(saved_result, graph_sidecar, saved_result_authenticated=True)",
        "evaluator_api": "ScopedEvaluator(graph, backend, free_variables, oracle).evaluate(root)",
        "real_roots": list(ROOTS),
        "binder_node": ARGUMENT_NODE_ID,
        "binder_interval": "[0,2^-128]",
    }
    for key, expected in expected_contract.items():
        _require(contract.get(key) == expected, "manifest input_contract." + key + " mismatch")

    scope = manifest.get("scope")
    _require(isinstance(scope, Mapping), "manifest scope missing")
    _require(scope.get("real_roots") == list(ROOTS), "manifest real roots mismatch")
    _require(scope.get("only_each_root_once") is True, "manifest root call guard missing")
    _require(scope.get("complete_g5") is False, "manifest complete g5 claim is forbidden")
    _require(scope.get("physical_m15") is False, "manifest physical m15 claim is forbidden")
    _require(scope.get("live_owner") is False, "manifest live owner claim is forbidden")

    files = manifest.get("files")
    _require(isinstance(files, list) and len(files) == len(FROZEN_INPUTS), "manifest file count")
    hashes: dict[str, str] = {}
    for index, row in enumerate(files):
        _require(isinstance(row, Mapping), f"manifest files[{index}] is not an object")
        key, path = _repo_file(row.get("path"), f"manifest files[{index}]")
        _require(key in FROZEN_SET, "manifest contains out-of-scope file: " + key)
        _require(key not in hashes, "manifest contains duplicate file: " + key)
        expected = _sha_text(row.get("sha256"), key)
        actual = _sha_file(path)
        _require(actual == expected, "source hash mismatch: " + key)
        if key in GZIP_INPUTS:
            _require(row.get("raw_sha256_before_decode") is True, "gzip raw pin missing: " + key)
        hashes[key] = actual
    _require(frozenset(hashes) == FROZEN_SET, "manifest does not contain exactly the frozen inputs")
    _require(hashes[GRAPH_REL] == GRAPH_SHA256, "real graph raw SHA is not the accepted pin")
    _require(hashes[PRIOR_REL] == PRIOR_SHA256, "prior receipt raw SHA is not the accepted pin")
    _require(
        hashes[SIGMA_SOURCE_REL] == SIGMA_SOURCE_SHA256,
        "sigma source metadata file SHA is not the accepted pin",
    )
    _require(
        hashes[SIGMA_PRIME_SOURCE_REL] == SIGMA_PRIME_SOURCE_SHA256,
        "sigma-prime source metadata file SHA is not the accepted pin",
    )
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


def _checkpoint_path(value: Any) -> Path:
    _require(isinstance(value, str) and value, "pending brief source checkpoint is missing")
    candidate = Path(value.replace("\\", "/"))
    if not candidate.is_absolute():
        candidate = (ROOT / candidate).resolve()
    else:
        candidate = candidate.resolve()
    _require(candidate == SOURCE.resolve(), "pending brief points to another source checkpoint")
    return candidate


def _validate_pending_brief(source_sha: str) -> str:
    _require(BRIEF.is_file(), "pending resume brief is missing")
    raw = BRIEF.read_bytes()
    brief = _decode_object(raw, "pending resume brief")
    record = _brief_record(brief)
    _require(record.get("component") == COMPONENT, "pending brief component mismatch")
    _require(record.get("phase") == PHASE, "pending brief phase mismatch")
    _require(record.get("one_entry_only") is True, "pending brief one-entry guard missing")
    _require(record.get("do_not_reenter") is True, "pending brief re-entry guard missing")
    checkpoint = record.get("source_checkpoint")
    _checkpoint_path(checkpoint.get("path") if isinstance(checkpoint, Mapping) else checkpoint)
    checkpoint_sha = checkpoint.get("sha256") if isinstance(checkpoint, Mapping) else None
    _require(
        _sha_text(checkpoint_sha, "pending source checkpoint") == source_sha,
        "pending brief checkpoint SHA mismatch",
    )
    return _sha_bytes(raw)


def _preflight() -> tuple[dict[str, str], str, str, dict[str, Any]]:
    _require(not ENTRY.exists() and not RESULT.exists(), "component already entered or has a result")
    _require(SOURCE.is_file(), "source manifest checkpoint is missing")
    source_raw = SOURCE.read_bytes()
    source_sha = _sha_bytes(source_raw)
    manifest = _decode_object(source_raw, "source manifest")
    brief_sha = _validate_pending_brief(source_sha)
    hashes = _manifest_hashes(manifest)
    return hashes, source_sha, brief_sha, manifest


def _read_authenticated(hashes: Mapping[str, str]) -> dict[str, bytes]:
    raw: dict[str, bytes] = {}
    for relative in FROZEN_INPUTS:
        data = (ROOT / relative).read_bytes()
        _require(_sha_bytes(data) == hashes[relative], "frozen input changed after entry: " + relative)
        raw[relative] = data
    _require(_sha_bytes(raw[GRAPH_REL]) == GRAPH_SHA256, "graph raw SHA changed before decode")
    _require(_sha_bytes(raw[PRIOR_REL]) == PRIOR_SHA256, "prior raw SHA changed before decode")
    _require(
        _sha_bytes(raw[SIGMA_SOURCE_REL]) == SIGMA_SOURCE_SHA256,
        "sigma source metadata SHA changed after entry",
    )
    _require(
        _sha_bytes(raw[SIGMA_PRIME_SOURCE_REL]) == SIGMA_PRIME_SOURCE_SHA256,
        "sigma-prime source metadata SHA changed after entry",
    )
    return raw


def _load_module(relative: str, module_name: str, allowed: set[str]) -> ModuleType:
    _require(relative in {LEAF_REL, FIXTURE_REL, EVALUATOR_REL}, "unexpected module path")
    before = set(sys.modules)
    spec = importlib.util.spec_from_file_location(module_name, ROOT / relative)
    _require(spec is not None and spec.loader is not None, "module loader missing: " + relative)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    for name in set(sys.modules) - before:
        if name == "openai_ns_reconstruction" or name in allowed:
            continue
        if name.startswith("openai_ns_reconstruction."):
            raise RunnerInputError("unexpected project module imported: " + name)
    return module


def _validate_fixture(value: Any) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "fixture result is not an object")
    result = dict(value)
    _require(result.get("schema") == FIXTURE_SCHEMA, "fixture schema mismatch")
    _require(result.get("passed") is True, "new oracle fixture did not pass")
    checks = result.get("checks")
    _require(
        isinstance(checks, Mapping) and all(item is True for item in checks.values()),
        "fixture checks failed",
    )
    _require(result.get("synthetic_only") is True, "fixture is not synthetic-only")
    for key, message in (
        ("used_prior_interval", "fixture consumed the saved interval"),
        ("used_saved_result", "fixture consumed the saved result"),
        ("used_real_graph", "fixture consumed the real graph"),
        ("used_graph_sidecar", "fixture consumed the real graph sidecar"),
        ("used_old_fixture", "fixture consumed an old fixture"),
        ("used_old_module", "fixture consumed an old module"),
        ("old_fixture_called", "fixture called an old fixture"),
        ("old_producer_called", "fixture called an old producer"),
    ):
        if key in result:
            _require(result.get(key) is False, message)
    return result


def _fraction(value: Any, label: str) -> Fraction:
    if isinstance(value, Fraction):
        return value
    if isinstance(value, Mapping):
        _require(value.get("type") == "Fraction", label + " is not a Fraction record")
        numerator = value.get("numerator")
        denominator = value.get("denominator")
    else:
        numerator = value
        denominator = 1
    _require(
        type(numerator) is int and type(denominator) is int,
        label + " is not exact integer data",
    )
    _require(denominator > 0, label + " denominator is not positive")
    return Fraction(numerator, denominator)


def _fraction_pair(value: Any, label: str) -> tuple[Fraction, Fraction]:
    _require(
        isinstance(value, (list, tuple)) and len(value) == 2,
        label + " is not a two-entry interval",
    )
    lower = _fraction(value[0], label + ".lower")
    upper = _fraction(value[1], label + ".upper")
    _require(lower <= upper, label + " lower exceeds upper")
    return lower, upper


def _fractionize(value: Any, label: str) -> Any:
    """Convert only explicit Fraction JSON records to exact Fraction objects."""

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


def _interval_bounds(value: Any, label: str) -> tuple[Fraction, Fraction]:
    if isinstance(value, Mapping):
        lower = _fraction(value.get("lower"), label + ".lower")
        upper = _fraction(value.get("upper"), label + ".upper")
    else:
        _require(hasattr(value, "lower") and hasattr(value, "upper"), label + " has no exact endpoints")
        lower = _fraction(getattr(value, "lower"), label + ".lower")
        upper = _fraction(getattr(value, "upper"), label + ".upper")
    _require(lower <= upper, label + " lower exceeds upper")
    return lower, upper


def _interval_alias_bounds(value: Any, label: str) -> tuple[Fraction, Fraction]:
    """Decode either saved JSON interval alias without losing exactness."""

    if isinstance(value, Mapping):
        nested = value.get("interval", _MISSING)
        if nested is not _MISSING:
            return _interval_alias_bounds(nested, label + ".interval")
        return _interval_bounds(value, label)
    return _fraction_pair(value, label)


def _named_interval_alias_bounds(
    mapping: Mapping[str, Any], names: tuple[str, ...], label: str
) -> tuple[Fraction, Fraction]:
    """Normalize every present exact interval alias before comparing aliases."""

    found = [
        (name, _interval_alias_bounds(mapping[name], label + "." + name))
        for name in names
        if name in mapping
    ]
    _require(found, label + " is missing")
    first_name, first = found[0]
    for name, value in found[1:]:
        _require(value == first, label + " aliases disagree: " + first_name + "/" + name)
    return first


def _saved_prior_result(receipt: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    candidate = receipt.get("result")
    if not isinstance(candidate, Mapping):
        candidate = receipt
    _require(candidate.get("status") == "returned", "prior sigma receipt is not returned")
    actual = candidate.get("actual_result")
    _require(isinstance(actual, Mapping), "prior sigma actual_result is missing")
    consumer = actual.get("consumer_result")
    direct = actual.get("direct_lower_collar")
    nested_direct = consumer.get("direct_lower_collar") if isinstance(consumer, Mapping) else None
    if not isinstance(direct, Mapping):
        direct = nested_direct
    elif isinstance(nested_direct, Mapping):
        _require(
            _interval_alias_bounds(
                direct.get("interval_leaf_input"), "prior direct interval input"
            )
            == _interval_alias_bounds(
                nested_direct.get("interval_leaf_input"), "prior nested direct interval input"
            ),
            "prior direct interval aliases disagree",
        )
        if "evaluated_argument_interval" in direct and "evaluated_argument_interval" in nested_direct:
            _require(
                _interval_alias_bounds(
                    direct["evaluated_argument_interval"],
                    "prior direct evaluated argument interval",
                )
                == _interval_alias_bounds(
                    nested_direct["evaluated_argument_interval"],
                    "prior nested evaluated argument interval",
                ),
                "prior evaluated argument aliases disagree",
            )
    if not isinstance(consumer, Mapping):
        consumer = actual
    _require(isinstance(direct, Mapping), "prior direct lower-collar result is missing")
    direct = dict(direct)
    _require(consumer.get("status") == "bounded", "prior sigma consumer is not bounded")
    _require(direct.get("scope") == "direct lower-collar chi sigma(s) only", "prior scope mismatch")
    _require(direct.get("fixed_s_in_z") is True, "prior fixed-s-in-Z condition missing")
    _require(direct.get("fixed_s_over_sc_in_z") is True, "prior fixed-s/sc-in-Z condition missing")
    _require(direct.get("interval_leaf_call_count") == 1, "prior direct interval call count mismatch")
    _require(direct.get("sigma_upper_is_not_a_value") is True, "prior sigma cap was selected as a value")
    _require(direct.get("sigma_value") is None, "prior sigma value was materialized")
    _require(consumer.get("source_identity_verified") is False, "prior source identity claim is forbidden")
    _require(consumer.get("owner_binding") is False, "prior owner claim is forbidden")
    _require(consumer.get("live_owner_binding") is False, "prior live-owner claim is forbidden")
    _require(consumer.get("physical_m15") is False, "prior physical m 15 claim is forbidden")
    direct_input = _interval_alias_bounds(
        direct.get("interval_leaf_input"), "prior direct interval input"
    )
    _require(direct_input == EXPECTED_BINDER_INTERVAL, "prior direct input is not [0, 2^-128]")
    evaluated_argument = direct.get("evaluated_argument_interval")
    _require(isinstance(evaluated_argument, Mapping), "prior evaluated argument interval is missing")
    evaluated_input = _interval_alias_bounds(
        evaluated_argument, "prior evaluated argument interval"
    )
    _require(evaluated_input == direct_input, "prior evaluated argument interval alias mismatch")
    _named_interval_alias_bounds(
        direct,
        ("sigma", "sigma_interval"),
        "prior direct sigma interval",
    )
    _named_interval_alias_bounds(
        direct,
        ("d_sigma_ds", "sigma_s_interval", "derivative"),
        "prior direct sigma_s interval",
    )
    return dict(candidate), direct


def _operation(record: Mapping[str, Any]) -> Any:
    operation = record.get("operation", _MISSING)
    shorthand = record.get("op", _MISSING)
    if operation is not _MISSING and shorthand is not _MISSING:
        _require(operation == shorthand, "graph operation/op disagreement")
    return operation if operation is not _MISSING else shorthand


def _validate_graph_scope(graph: Mapping[str, Any]) -> dict[str, Any]:
    nodes = graph.get("exact_graph_nodes")
    _require(
        isinstance(nodes, list) and len(nodes) > max(ROOTS),
        "exact_graph_nodes is missing target roots",
    )
    argument = nodes[ARGUMENT_NODE_ID]
    sigma = nodes[SIGMA_NODE_ID]
    sigma_prime = nodes[SIGMA_PRIME_NODE_ID]
    _require(isinstance(argument, Mapping), "binder graph row is not an object")
    _require(isinstance(sigma, Mapping), "sigma graph row is not an object")
    _require(isinstance(sigma_prime, Mapping), "sigma-prime graph row is not an object")
    _require(_operation(argument) == "bound_variable", "binder graph row is not bound_variable")
    _require(argument.get("name") == "original_lower_collar_s", "binder graph owner name mismatch")
    _require(_operation(sigma) == "exact_original_source_flat_sigma", "sigma graph operation mismatch")
    _require(sigma.get("argument") == ARGUMENT_NODE_ID, "sigma graph binder mismatch")
    _require(_operation(sigma_prime) == "analytic_unary", "sigma-prime graph operation mismatch")
    prime_name = sigma_prime.get("name", _MISSING)
    prime_function = sigma_prime.get("function", _MISSING)
    if prime_name is not _MISSING and prime_function is not _MISSING:
        _require(prime_name == prime_function, "sigma-prime name/function disagreement")
    prime_owner = prime_name if prime_name is not _MISSING else prime_function
    _require(prime_owner == "original_flat_sigma_prime", "sigma-prime owner mismatch")
    _require(sigma_prime.get("argument") == ARGUMENT_NODE_ID, "sigma-prime graph binder mismatch")
    return {
        "node_count": len(nodes),
        "selected_nodes": {
            "35848": dict(argument),
            "35949": dict(sigma),
            "35953": dict(sigma_prime),
        },
        "owner_check": "conditional graph-row metadata only",
        "graph_binding": False,
        "source_identity_verified": False,
        "live_owner_binding": False,
    }


def _validate_source_hash_bindings(
    graph_scope: Mapping[str, Any], hashes: Mapping[str, str]
) -> dict[str, Any]:
    """Check recorded module hashes without importing either source module."""

    selected = graph_scope.get("selected_nodes")
    _require(isinstance(selected, Mapping), "graph scope selected nodes are missing")
    sigma = selected.get(str(SIGMA_NODE_ID))
    sigma_prime = selected.get(str(SIGMA_PRIME_NODE_ID))
    _require(isinstance(sigma, Mapping), "sigma selected row is missing")
    _require(isinstance(sigma_prime, Mapping), "sigma-prime selected row is missing")
    sigma_binding = sigma.get("source_binding")
    prime_binding = sigma_prime.get("derivative_source")
    _require(isinstance(sigma_binding, Mapping), "sigma source binding is missing")
    _require(isinstance(prime_binding, Mapping), "sigma-prime source binding is missing")
    _require(
        sigma_binding.get("module") == Path(SIGMA_SOURCE_REL).name,
        "sigma source module name does not match frozen path",
    )
    _require(
        prime_binding.get("module") == Path(SIGMA_PRIME_SOURCE_REL).name,
        "sigma-prime source module name does not match frozen path",
    )
    _require(
        sigma_binding.get("module_sha256") == hashes[SIGMA_SOURCE_REL],
        "sigma graph module SHA does not match frozen source",
    )
    _require(
        prime_binding.get("module_sha256") == hashes[SIGMA_PRIME_SOURCE_REL],
        "sigma-prime graph module SHA does not match frozen source",
    )
    return {
        "sigma_module": sigma_binding.get("module"),
        "sigma_module_sha256": sigma_binding.get("module_sha256"),
        "sigma_prime_module": prime_binding.get("module"),
        "sigma_prime_module_sha256": prime_binding.get("module_sha256"),
        "source_identity_verified": False,
        "live_owner_binding": False,
    }


def _build_oracle(
    leaf: ModuleType,
    saved_result: Mapping[str, Any],
    graph_sidecar: Mapping[str, Any],
) -> Any:
    factory = getattr(leaf, "build_oracle", None)
    _require(callable(factory), "new oracle leaf has no build_oracle API")
    oracle = factory(saved_result, graph_sidecar, saved_result_authenticated=True)
    _require(callable(oracle), "new oracle factory did not return a callable")
    return oracle


def _validate_evaluation(
    value: Any,
    root: int,
    expected_bounds: tuple[Fraction, Fraction],
) -> dict[str, Any]:
    _require(isinstance(value, Mapping), "scoped evaluator result is not an object")
    result = dict(value)
    _require(result.get("schema") == EVALUATOR_RESULT_SCHEMA, "scoped result schema mismatch")
    _require(result.get("status") == "conditional", "scoped result is not conditional")
    _require(result.get("conditional") is True, "scoped result conditional marker missing")
    _require(result.get("certified") is False, "scoped result is certified unexpectedly")
    _require(result.get("root") == root, "scoped result root mismatch")
    returned = result.get("value")
    _require(
        type(returned).__name__ == "ScopedSigmaInterval",
        "oracle did not return a ScopedSigmaInterval",
    )
    _require(
        _interval_bounds(returned, "scoped oracle result") == expected_bounds,
        "oracle did not return the saved interval for this root",
    )
    claims = result.get("claims")
    _require(isinstance(claims, Mapping), "scoped result claims missing")
    for key in ("function_equality", "physical_m15", "Kchi", "Navier_Stokes", "certified"):
        _require(claims.get(key) is False, "forbidden scoped claim: " + key)
    provenance = result.get("provenance")
    _require(isinstance(provenance, Mapping), "scoped result provenance missing")
    for key in ("source_identity_verified", "backend_verified", "oracle_verified", "graph_hash_authenticated"):
        _require(provenance.get(key) is False, "forbidden scoped provenance: " + key)
    return result


def _consume(hashes: Mapping[str, str], progress: dict[str, Any]) -> dict[str, Any]:
    # Nothing below this line runs before the entry and successful fixture.
    allowed = {
        "openai_ns_reconstruction",
        "openai_ns_reconstruction." + COMPONENT,
        "openai_ns_reconstruction." + COMPONENT + "_fixture",
        "openai_ns_reconstruction.n1_scoped_graph_evaluator_v1",
    }
    leaf = _load_module(LEAF_REL, "openai_ns_reconstruction." + COMPONENT, allowed)
    fixture = _load_module(
        FIXTURE_REL,
        "openai_ns_reconstruction." + COMPONENT + "_fixture",
        allowed,
    )
    evaluator_module = _load_module(
        EVALUATOR_REL,
        "openai_ns_reconstruction.n1_scoped_graph_evaluator_v1",
        allowed,
    )
    progress["modules_imported_after_entry"] = True
    _require(getattr(leaf, "SCHEMA", None) == LEAF_SCHEMA, "new oracle leaf schema mismatch")
    _require(getattr(fixture, "SCHEMA", None) == FIXTURE_SCHEMA, "new oracle fixture schema mismatch")
    _require(
        getattr(evaluator_module, "SCHEMA", None) == EVALUATOR_SCHEMA,
        "scoped evaluator schema mismatch",
    )

    fixture_runner = getattr(fixture, "run_fixture", None)
    _require(callable(fixture_runner), "new oracle fixture has no run_fixture(leaf) API")
    progress["fixture_calls"] = 1
    progress["fixture_result_raw"] = None
    fixture_result_raw = fixture_runner(leaf)
    progress["fixture_result_raw"] = fixture_result_raw
    fixture_result = _validate_fixture(fixture_result_raw)
    progress["fixture_result"] = fixture_result

    # Only after the fixture has completed are the authenticated gzip bytes
    # decompressed. The raw hash is checked again immediately before decode.
    raw = _read_authenticated(hashes)
    graph = _decode_gzip_object(raw[GRAPH_REL], GRAPH_REL)
    receipt = _decode_gzip_object(raw[PRIOR_REL], PRIOR_REL)
    progress["raw_hash_authenticated_before_decode"] = True
    progress["gzip_decoded_after_entry_and_fixture"] = True

    prior_result, prior_direct = _saved_prior_result(receipt)
    prior_bounds_by_node = {
        SIGMA_NODE_ID: _named_interval_alias_bounds(
            prior_direct,
            ("sigma", "sigma_interval"),
            "prior direct sigma interval",
        ),
        SIGMA_PRIME_NODE_ID: _named_interval_alias_bounds(
            prior_direct,
            ("d_sigma_ds", "sigma_s_interval", "derivative"),
            "prior direct sigma_s interval",
        ),
    }
    saved_result_adapter = _fractionize(prior_direct, "prior direct lower-collar")
    graph_scope = _validate_graph_scope(graph)
    source_hash_bindings = _validate_source_hash_bindings(graph_scope, hashes)
    progress["prior_direct_interval"] = prior_direct
    progress["prior_bounds_by_node"] = prior_bounds_by_node
    progress["saved_result_fraction_adapter"] = saved_result_adapter
    progress["graph_scope"] = graph_scope
    progress["source_hash_bindings"] = source_hash_bindings

    graph_sidecar = {
        "authenticated": True,
        "graph_sha256": GRAPH_SHA256,
        "nodes": {
            ARGUMENT_NODE_ID: graph_scope["selected_nodes"]["35848"],
            SIGMA_NODE_ID: graph_scope["selected_nodes"]["35949"],
            SIGMA_PRIME_NODE_ID: graph_scope["selected_nodes"]["35953"],
        },
    }
    oracle = _build_oracle(leaf, saved_result_adapter, graph_sidecar)
    fraction_interval = getattr(leaf, "FractionInterval", None)
    sigma_interval = getattr(leaf, "ScopedSigmaInterval", None)
    _require(fraction_interval is not None, "new oracle leaf FractionInterval is missing")
    _require(sigma_interval is not None, "new oracle leaf ScopedSigmaInterval is missing")
    binder_interval = fraction_interval(*EXPECTED_BINDER_INTERVAL)
    oracle_calls: list[dict[str, Any]] = []
    progress["oracle_calls"] = oracle_calls
    progress["oracle_call_count"] = 0

    def counting_oracle(node_id: int, record: Any, category: str, env: Any) -> Any:
        _require(node_id in ROOTS, "oracle reached a non-target root")
        lookup = getattr(env, "lookup", None)
        _require(callable(lookup), "scoped evaluator did not provide lexical lookup")
        _require(
            lookup(ARGUMENT_NODE_ID) == binder_interval,
            "binder 35848 did not receive [0, 2^-128]",
        )
        expected_bounds = prior_bounds_by_node[node_id]
        oracle_calls.append(
            {
                "node_id": node_id,
                "category": category,
                "record": record,
                "expected_bounds": expected_bounds,
            }
        )
        progress["oracle_call_count"] = len(oracle_calls)
        value = oracle(node_id, record, category, env)
        _require(type(value).__name__ == "ScopedSigmaInterval", "oracle returned the wrong interval type")
        _require(
            _interval_bounds(value, "oracle value") == expected_bounds,
            "oracle returned a value different from this root's saved interval",
        )
        return value

    scoped_evaluator = getattr(evaluator_module, "ScopedEvaluator", None)
    _require(scoped_evaluator is not None, "ScopedEvaluator class is missing")
    evaluator = scoped_evaluator(
        graph,
        _OpaqueScopedBackend((fraction_interval, sigma_interval)),
        free_variables={ARGUMENT_NODE_ID: binder_interval},
        oracle=counting_oracle,
    )
    real_results: list[dict[str, Any]] = []
    real_results_raw: list[Any] = []
    progress["real_results"] = real_results
    progress["real_results_raw"] = real_results_raw
    progress["real_root_evaluation_calls"] = 0
    progress["current_root"] = None
    for root in ROOTS:
        progress["current_root"] = root
        progress["real_root_evaluation_calls"] = len(real_results) + 1
        evaluated = evaluator.evaluate(root)
        real_results_raw.append(evaluated)
        progress["real_results_raw"] = real_results_raw
        validated = _validate_evaluation(
            evaluated,
            root,
            prior_bounds_by_node[root],
        )
        real_results.append(validated)
        progress["real_results"] = real_results
    progress["current_root"] = None
    _require(
        [row["node_id"] for row in oracle_calls] == list(ROOTS),
        "real roots were not each called once",
    )
    _require(len(oracle_calls) == 2, "oracle call count is not exactly two")

    progress["scoped_evaluator_instances"] = 1
    progress["real_root_evaluation_calls"] = len(real_results)
    progress["oracle_calls"] = oracle_calls
    progress["real_results"] = real_results
    return {
        "fixture": fixture_result,
        "prior_result_status": prior_result.get("status"),
        "prior_direct_interval": prior_direct,
        "prior_bounds_by_node": prior_bounds_by_node,
        "prior_sigma_interval": {
            "lower": prior_bounds_by_node[SIGMA_NODE_ID][0],
            "upper": prior_bounds_by_node[SIGMA_NODE_ID][1],
        },
        "prior_sigma_s_interval": {
            "lower": prior_bounds_by_node[SIGMA_PRIME_NODE_ID][0],
            "upper": prior_bounds_by_node[SIGMA_PRIME_NODE_ID][1],
        },
        "graph_scope": graph_scope,
        "source_hash_bindings": source_hash_bindings,
        "graph_sidecar": graph_sidecar,
        "binder_node": ARGUMENT_NODE_ID,
        "binder_interval": binder_interval,
        "real_roots": list(ROOTS),
        "oracle_calls": oracle_calls,
        "real_results": real_results,
        "scoped_evaluator_instances": 1,
        "real_root_evaluation_calls": len(real_results),
        "complete_g5": False,
        "ancestor_reexecuted": False,
        "old_fixture_called": False,
        "old_producer_called": False,
        "physical_m15": False,
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
                "manifest_schema": manifest.get("schema"),
                "manifest_sha256": source_sha,
                "brief_sha256": brief_sha,
                "input_hashes": hashes,
                "runner_path": RUNNER_REL,
                "source_manifest_self_hash_forbidden": True,
                "raw_gzip_authenticated_before_decode": True,
                "entry_written_before_import": True,
                "fixture_must_complete_before_decode": True,
                "entered": True,
            },
        )
    except BaseException:
        print("status=entry_failed")
        return 1

    progress: dict[str, Any] = {
        "modules_imported_after_entry": False,
        "fixture_calls": 0,
        "raw_hash_authenticated_before_decode": False,
        "gzip_decoded_after_entry_and_fixture": False,
        "scoped_evaluator_instances": 0,
        "real_root_evaluation_calls": 0,
        "oracle_calls": [],
        "real_results": [],
    }
    outcome: dict[str, Any] = {
        "schema": COMPONENT + "_result",
        "component": COMPONENT,
        "phase": PHASE,
        "manifest_sha256": source_sha,
        "brief_sha256": brief_sha,
        "entry_sha256": entry_sha,
        "input_hashes": hashes,
        "status": "failed",
        "progress": progress,
        "error": None,
        "source_manifest_self_hash_forbidden": True,
        "raw_gzip_authenticated_before_decode": True,
        "complete_g5": False,
        "ancestor_reexecuted": False,
        "old_fixture_called": False,
        "old_producer_called": False,
        "physical_m15": False,
        "source_identity_verified": False,
        "live_owner_binding": False,
    }
    try:
        outcome.update(_consume(hashes, progress))
        outcome["status"] = "returned"
    except BaseException as exc:
        outcome["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }

    try:
        _write_exclusive(RESULT, outcome)
    except BaseException:
        print("status=result_write_failed")
        return 1
    print("status=" + str(outcome["status"]))
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
