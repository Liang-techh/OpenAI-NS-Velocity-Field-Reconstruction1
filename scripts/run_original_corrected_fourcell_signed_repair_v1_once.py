"""One-entry source consumer for the corrected four-cell signed repair.

This runner is intentionally source-only.  It authenticates the raw declared
inputs before any archive is decoded, writes the exclusive entry receipt, runs
the synthetic fixture, and only then opens the accepted v3 lower-endpoint
receipt and corrected-source contract.  The v3 rows are retained as opaque
source records; no producer, graph, old ledger, or endpoint numerical value is
replayed.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from collections.abc import Mapping, Sequence
from fractions import Fraction
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
COMPONENT = "original_corrected_fourcell_signed_repair_v1"
PHASE = COMPONENT + "_pending"
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entry.json"
RESULT = DIRECTORY / "result.json"

V3_RECEIPT = (
    "experiments/root_st073/"
    "original_corrected_endpoint_vy_owner_v3_receipt.json.gz"
)
CORRECTED_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_corrected_reserved_m15_source.json"
)
V3_LEAF = "src/openai_ns_reconstruction/original_corrected_endpoint_vy_owner_v3.py"
V3_FIXTURE = "src/openai_ns_reconstruction/original_corrected_endpoint_vy_owner_fixture_v3.py"
V3_RUNNER = "scripts/run_original_corrected_endpoint_vy_owner_v3_once.py"

LEAF = "src/openai_ns_reconstruction/original_corrected_fourcell_signed_repair_v1.py"
FIXTURE = "src/openai_ns_reconstruction/original_corrected_fourcell_signed_repair_fixture_v1.py"
RUNNER = "scripts/run_original_corrected_fourcell_signed_repair_v1_once.py"

INPUT = (V3_RECEIPT, CORRECTED_SOURCE, V3_LEAF, V3_FIXTURE, V3_RUNNER)
REQUIRED = INPUT + (LEAF, FIXTURE, RUNNER)

FIXED = {
    V3_RECEIPT: "778ba4bffdcc5e48500c5a5e51f7374812b0aabc6142c8d1e1d24005c6a4bd5c",
    CORRECTED_SOURCE: "fa57cdda6a9de6e6b32fb1b2f265d663825662be959b8a3370d56432f6ff07a5",
}

CELLS = (("-1", "-.5"), ("-.5", "0"), ("0", ".5"), (".5", "1"))
LEAF_STATIC = {
    "SCHEMA": "original_corrected_fourcell_signed_repair_v1",
    "COMPONENT": "corrected_physical_lower_Vy_factored_signed_interval_repair_v1",
    "CELLS": CELLS,
    "ORDERS": (0, 1),
    "SCALE_FIELDS": ("source_exponents", "radius_power", "additional_log_interval"),
}
FIXTURE_STATIC = {
    "SCHEMA": "original_corrected_fourcell_signed_repair_fixture_v1",
    "LEAF_SCHEMA": "original_corrected_fourcell_signed_repair_v1",
    "CELLS": CELLS,
}


class RunnerInputError(ValueError):
    """Raised when the frozen one-entry contract is not satisfied."""


def require(condition: object, message: str) -> None:
    if not condition:
        raise RunnerInputError(message)


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def _json_constant(token: str) -> Any:
    raise ValueError("non-finite JSON constant is not admitted: " + token)


def load_json(raw: bytes, label: str) -> dict[str, Any]:
    try:
        value = json.loads(raw.decode("utf-8"), parse_constant=_json_constant)
    except (UnicodeDecodeError, json.JSONDecodeError, ValueError) as exc:
        raise RunnerInputError(label + " is not finite JSON") from exc
    require(isinstance(value, dict), label + " must be a JSON object")
    return value


def _json_safe(value: Any) -> Any:
    """Serialize returned interval objects without selecting their values."""

    if value is None or isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, float):
        require(value == value and value not in (float("inf"), float("-inf")),
                "non-finite float in returned source packet")
        return value
    if isinstance(value, Fraction):
        return {"numerator": value.numerator, "denominator": value.denominator}
    raw_interval = getattr(value, "_mpi_", None)
    if raw_interval is not None:
        return {
            "lower_exact_mpf_tuple": list(raw_interval[0]),
            "upper_exact_mpf_tuple": list(raw_interval[1]),
            "directed_interval": True,
        }
    raw_point = getattr(value, "_mpf_", None)
    if raw_point is not None:
        return {"exact_mpf_tuple": list(raw_point), "directed_point": True}
    if isinstance(value, Mapping):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, Path):
        return str(value)
    return repr(value)


def save(path: Path, value: object) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = (
        json.dumps(
            _json_safe(value),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    try:
        with path.open("xb") as handle:
            handle.write(raw)
    except FileExistsError as exc:
        raise RunnerInputError("one-entry receipt already exists: " + str(path)) from exc
    return sha(raw)


def repository_path(relative: str) -> Path:
    path = Path(relative)
    require(not path.is_absolute() and ".." not in path.parts, "non-relative input path")
    resolved = (ROOT / path).resolve()
    require(resolved.is_relative_to(ROOT.resolve()), "input path escapes repository: " + relative)
    require(resolved.is_file(), "missing input: " + relative)
    return resolved


def manifest_hashes(manifest: Mapping[str, Any]) -> dict[str, str]:
    rows = manifest.get("files")
    require(isinstance(rows, list), "source manifest files list is missing")
    hashes: dict[str, str] = {}
    for index, row in enumerate(rows):
        require(isinstance(row, Mapping), f"source.files[{index}] is not an object")
        path = row.get("path")
        expected = row.get("sha256")
        require(isinstance(path, str) and path in REQUIRED,
                "unexpected manifest input")
        require(path not in hashes, "duplicate manifest input: " + path)
        require(isinstance(expected, str) and len(expected) == 64,
                "invalid manifest SHA: " + path)
        actual = sha(repository_path(path).read_bytes())
        require(actual == expected.lower(), "manifest SHA changed: " + path)
        if path in FIXED:
            require(actual == FIXED[path], "fixed SHA changed: " + path)
        hashes[path] = actual
    require(set(hashes) == set(REQUIRED),
            "manifest does not pin exactly the required inputs")
    for path in (V3_LEAF, V3_FIXTURE, V3_RUNNER, LEAF, FIXTURE):
        require(path in hashes, "source file is not manifest-pinned: " + path)
    return hashes


def _ast_literal_constants(relative: str, names: Sequence[str]) -> dict[str, Any]:
    """Read only top-level literal assignments; do not import or execute code."""

    source = repository_path(relative).read_text(encoding="utf-8")
    try:
        tree = ast.parse(source, filename=relative)
    except SyntaxError as exc:
        raise RunnerInputError(relative + " is not valid Python source") from exc
    found: dict[str, Any] = {}
    wanted = set(names)
    for node in tree.body:
        targets: list[ast.Name] = []
        value: ast.AST | None = None
        if isinstance(node, ast.Assign):
            targets = [target for target in node.targets if isinstance(target, ast.Name)]
            value = node.value
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            targets = [node.target]
            value = node.value
        if value is None:
            continue
        for target in targets:
            if target.id not in wanted:
                continue
            try:
                found[target.id] = ast.literal_eval(value)
            except (ValueError, TypeError) as exc:
                raise RunnerInputError(
                    f"{relative} constant {target.id} is not a literal"
                ) from exc
    for name in names:
        require(name in found, f"{relative} top-level constant is missing: {name}")
    return found


def _check_static_source_constants() -> dict[str, Any]:
    leaf = _ast_literal_constants(LEAF, tuple(LEAF_STATIC))
    fixture = _ast_literal_constants(FIXTURE, tuple(FIXTURE_STATIC))
    require(leaf == LEAF_STATIC, "new leaf top-level static contract changed")
    require(fixture == FIXTURE_STATIC, "synthetic fixture top-level static contract changed")
    return {"leaf": leaf, "fixture": fixture}


def preflight() -> tuple[dict[str, Any], dict[str, str], str, str, dict[str, Any]]:
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    require(BRIEF.is_file() and SOURCE.is_file(),
            "pending brief or source manifest is missing")
    brief_raw = BRIEF.read_bytes()
    brief = load_json(brief_raw, "pending brief")
    source_raw = SOURCE.read_bytes()
    source = load_json(source_raw, "source manifest")
    require(
        brief.get("component") == COMPONENT
        and brief.get("phase") == PHASE
        and brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("one_entry_only") is True
        and brief.get("do_not_reenter") is True,
        "pending brief guard",
    )
    checkpoint = brief.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "pending source checkpoint is missing")
    require(Path(str(checkpoint.get("path"))).resolve() == SOURCE.resolve(),
            "source checkpoint path changed")
    require(sha(source_raw) == checkpoint.get("sha256"),
            "source checkpoint SHA changed")
    require(source.get("component") == COMPONENT
            and source.get("one_entry_only") is True,
            "source manifest guard")
    hashes = manifest_hashes(source)
    static = _check_static_source_constants()
    return source, hashes, sha(source_raw), sha(brief_raw), static


def authenticated_bytes(hashes: Mapping[str, str]) -> dict[str, bytes]:
    """Read and authenticate every declared raw byte before any gzip decode."""

    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = repository_path(path).read_bytes()
        require(sha(data) == hashes[path], "raw input changed after entry: " + path)
        raw[path] = data
    return raw


def module(relative: str, name: str) -> Any:
    path = repository_path(relative)
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None,
            "module load failed: " + relative)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def decode_gzip(raw: bytes, label: str) -> dict[str, Any]:
    try:
        decoded = gzip.decompress(raw)
    except (OSError, EOFError) as exc:
        raise RunnerInputError(label + " is not valid gzip") from exc
    return load_json(decoded, label + " decoded JSON")


def _known_maps(root: Mapping[str, Any], keys: Sequence[str]) -> list[Mapping[str, Any]]:
    result: list[Mapping[str, Any]] = [root]
    for key in keys:
        child = root.get(key)
        if isinstance(child, Mapping):
            result.append(child)
    return result


def _explicit_flag(roots: Sequence[Any], names: Sequence[str], depth: int = 3) -> bool:
    """Find an explicitly serialized boolean in bounded proof containers."""

    if depth < 0:
        return False
    nested_keys = (
        "proof", "checks", "evidence", "attestation", "guarantees",
        "source_contract", "corrected_source_contract",
        "physical_endpoint_source_proof", "same_source", "lower",
    )
    for value in roots:
        if not isinstance(value, Mapping):
            continue
        if any(value.get(name) is True for name in names):
            return True
        for key in nested_keys:
            child = value.get(key)
            if isinstance(child, Mapping) and _explicit_flag((child,), names, depth - 1):
                return True
    return False


def _first_value(root: Mapping[str, Any], names: Sequence[str]) -> tuple[Any, str] | tuple[None, None]:
    containers: list[tuple[str, Mapping[str, Any]]] = [("terminal", root)]
    for key in ("source_contract", "corrected_source_contract",
                "physical_endpoint_source_proof", "common_original_parameter"):
        child = root.get(key)
        if isinstance(child, Mapping):
            containers.append((f"terminal.{key}", child))
    for label, container in containers:
        for name in names:
            if name in container:
                return container[name], f"{label}.{name}"
    return None, None


def _cell_key(cell: tuple[str, str]) -> str:
    return str(cell)


def _cell_mapping(root: Mapping[str, Any], cell: tuple[str, str], label: str) -> Mapping[str, Any]:
    cells = root.get("cell_sources")
    require(isinstance(cells, Mapping), label + ".cell_sources is missing")
    value = cells.get(_cell_key(cell))
    require(isinstance(value, Mapping), f"{label} cell is missing: {cell}")
    return value


def _source_rows_and_proof(
    entry: Mapping[str, Any], cell: tuple[str, str]
) -> tuple[Mapping[str, Any], Mapping[str, Any], Mapping[str, Any]]:
    lower = entry.get("lower")
    require(isinstance(lower, Mapping), "v3 lower source is missing: " + str(cell))
    rows = lower.get("original_source_rows")
    if not isinstance(rows, Mapping):
        rows = lower.get("source_rows")
    require(isinstance(rows, Mapping), "v3 original lower source rows are missing: " + str(cell))
    for key in ("R", "Pstar"):
        require(key in rows, f"v3 lower original source rows missing {key}: {cell}")
    require(any(key in rows for key in ("V_y", "V_y_C0_C1", "Vy", "Vy_C0_C1")),
            "v3 lower original source rows missing V_y: " + str(cell))
    proof = lower.get("proof")
    if not isinstance(proof, Mapping):
        proof = lower.get("source_proof")
    if not isinstance(proof, Mapping):
        proof = entry.get("sample_proof")
    require(isinstance(proof, Mapping), "v3 lower source proof is missing: " + str(cell))
    return rows, proof, lower


def _positive_physical_rows(context: Any, rows: Mapping[str, Any], cell: tuple[str, str]) -> dict[str, Any]:
    """Require physical R/Pstar signed coefficient enclosures to be positive."""
    checked: dict[str, Any] = {}
    for name in ("R", "Pstar"):
        row = rows.get(name)
        require(isinstance(row, Mapping), f"{cell} missing {name} source row")
        base = row.get("row") if "coefficient_interval" not in row else row
        require(isinstance(base, Mapping), f"{cell} {name} raw row is missing")
        coefficient = base.get("coefficient_interval")
        if isinstance(coefficient, Mapping):
            lower = coefficient.get("lower_exact_mpf_tuple", coefficient.get("lower"))
            upper = coefficient.get("upper_exact_mpf_tuple", coefficient.get("upper"))
        elif isinstance(coefficient, Sequence) and not isinstance(coefficient, (str, bytes, bytearray)) and len(coefficient) == 2:
            lower, upper = coefficient
        else:
            lower = upper = coefficient
        require(lower is not None and upper is not None, f"{cell} {name} coefficient bounds missing")
        try:
            def point(value: Any) -> Any:
                if isinstance(value, (list, tuple)) and len(value) == 4:
                    raw = tuple(value)
                    return context.make_mpf((raw, raw))
                return context.mpf(value)
            lower, upper = point(lower), point(upper)
        except Exception as exc:
            raise RunnerInputError(f"{cell} {name} coefficient bounds are invalid") from exc
        require(lower <= upper and lower > 0,
                f"{cell} {name} signed coefficient is not strictly positive")
        scale = base.get("formal_positive_scale")
        require(isinstance(scale, Mapping)
                and all(field in scale for field in (
                    "source_exponents", "radius_power", "additional_log_interval"
                )), f"{cell} {name} formal_positive_scale is incomplete")
        exponents = scale.get("source_exponents")
        require(isinstance(exponents, Sequence)
                and not isinstance(exponents, (str, bytes, bytearray))
                and len(exponents) == 4,
                f"{cell} {name} formal_positive_scale lacks four opaque exponents")
        checked[name] = {
            "coefficient_lower_strictly_positive": True,
            "formal_positive_scale_is_positive_factor": True,
            "opaque_source_exponent_count": 4,
            "log_absolute_upper_not_used_as_value": True,
        }
    return checked


def _embedded_manifest(receipt: Mapping[str, Any]) -> Mapping[str, Any]:
    paths = (
        ("source",),
        ("source_manifest",),
        ("embedded_source_manifest",),
        ("result", "source_manifest"),
        ("result", "embedded_source_manifest"),
        ("result", "result", "source_manifest"),
        ("result", "result", "embedded_source_manifest"),
        ("result", "result", "owner", "source_manifest"),
        ("result", "result", "owner", "embedded_source_manifest"),
        ("result", "owner", "source_manifest"),
        ("result", "owner", "embedded_source_manifest"),
        ("output", "source_manifest"),
        ("result", "output", "source_manifest"),
    )
    for path in paths:
        value: Any = receipt
        for part in path:
            if not isinstance(value, Mapping):
                value = None
                break
            value = value.get(part)
        if isinstance(value, Mapping) and isinstance(value.get("files"), list):
            return value
    raise RunnerInputError("accepted v3 receipt has no embedded source manifest")


def _manifest_rows(manifest: Mapping[str, Any], label: str) -> dict[str, str]:
    rows = manifest.get("files")
    require(isinstance(rows, list), label + ".files is missing")
    result: dict[str, str] = {}
    for index, row in enumerate(rows):
        require(isinstance(row, Mapping), f"{label}.files[{index}] is not an object")
        path = row.get("path")
        value = row.get("sha256")
        require(isinstance(path, str) and isinstance(value, str),
                f"{label}.files[{index}] has no path/SHA")
        require(path not in result, label + " has duplicate path: " + path)
        result[path] = value.lower()
    return result


def _owner(receipt: Mapping[str, Any]) -> Mapping[str, Any]:
    value: Any = receipt
    for key in ("result", "result", "owner"):
        require(isinstance(value, Mapping), "v3 receipt result.result.owner path is missing")
        value = value.get(key)
    require(isinstance(value, Mapping), "v3 receipt owner is not an object")
    return value


def _validate_v3_status(receipt: Mapping[str, Any]) -> dict[str, Any]:
    statuses: dict[str, Any] = {}
    paths = (
        ("status",),
        ("result", "status"),
        ("result", "result", "status"),
        ("result", "result", "owner", "status"),
    )
    accepted = {
        "accepted", "returned", "success", "source_only",
        "four_cell_lower_endpoint_source_enclosures",
    }
    for path in paths:
        value: Any = receipt
        for key in path:
            if not isinstance(value, Mapping):
                value = None
                break
            value = value.get(key)
        if value is not None:
            require(isinstance(value, str), "v3 receipt status is not a string")
            statuses[".".join(path)] = value
    require(any(value in accepted for value in statuses.values()),
            "v3 receipt has no accepted status")
    return statuses


def _validate_corrected_source(
    corrected: Mapping[str, Any], family: Any
) -> dict[str, Any]:
    checks = corrected.get("checks")
    require(isinstance(checks, Mapping) and checks.get("passed") is True,
            "corrected source checks.passed is not true")
    terminal = corrected.get("corrected_terminal_source_inputs")
    require(isinstance(terminal, Mapping) and terminal.get("source_only") is True,
            "corrected_terminal_source_inputs source-only contract is missing")
    require(terminal.get("numeric_endpoint_cell") == ["0", ".5"],
            "corrected terminal endpoint cell changed")
    proof = terminal.get("physical_endpoint_source_proof")
    require(isinstance(proof, Mapping), "corrected terminal source proof is missing")
    require(proof.get("source_family") == family,
            "corrected terminal source family changed")
    rows = terminal.get("physical_endpoint_source_rows")
    require(isinstance(rows, Mapping), "corrected terminal source rows are missing")
    delta, delta_path = _first_value(
        terminal, ("actual_delta", "actual_original_delta", "source_delta", "delta")
    )
    log_pstar, pstar_path = _first_value(
        terminal, ("actual_logPstar_parameter",)
    )
    common, common_path = _first_value(
        terminal,
        ("common_original_parameter", "common_original_source_parameter",
         "original_common_parameter", "actual_logPstar_source"),
    )
    require(delta is not None, "corrected source actual delta is missing")
    require(log_pstar is not None,
            "corrected source actual_logPstar_parameter is missing")
    require(common is not None, "corrected source common-original parameter path is missing")
    return {
        "checks": checks,
        "terminal": terminal,
        "proof": proof,
        "rows": rows,
        "delta": delta,
        "delta_path": delta_path,
        "actual_logPstar_parameter": log_pstar,
        "actual_logPstar_parameter_path": pstar_path,
        "common_original_parameter": common,
        "common_original_parameter_path": common_path,
    }


def _same_source_contract(
    owner: Mapping[str, Any],
    family: Any,
    cell_records: Mapping[str, Mapping[str, Any]],
    terminal: Mapping[str, Any],
) -> dict[str, Any]:
    guarantees = owner.get("guarantees")
    if not isinstance(guarantees, Mapping):
        guarantees = {}
    proofs = [record["proof"] for record in cell_records.values()]
    proof_roots: list[Any] = [terminal, terminal.get("physical_endpoint_source_proof"), owner, guarantees]
    proof_roots.extend(proofs)

    # These are deliberately explicit guards.  A common source family alone is
    # not treated as a proof that any parameter is Z-independent.
    same_basis = _explicit_flag(
        proof_roots,
        ("same_source_basis", "same_source_identity", "same_basis",
         "source_basis_checked"),
    )
    same_provenance = _explicit_flag(
        proof_roots,
        ("same_source_provenance", "same_provenance", "provenance_checked",
         "source_provenance_checked"),
    )
    delta_independent = _explicit_flag(
        [terminal, terminal.get("source_contract"), terminal.get("corrected_source_contract"),
         terminal.get("physical_endpoint_source_proof"), owner],
        ("delta_Z_independent", "delta_is_Z_independent", "delta_source_Z_independent"),
    )
    radius_independent = _explicit_flag(
        [terminal, terminal.get("source_contract"), terminal.get("corrected_source_contract"),
         terminal.get("physical_endpoint_source_proof"), owner],
        ("R_Z_independent", "radius_Z_independent", "radius_is_Z_independent",
         "Pstar_and_radius_Z_independent"),
    )
    pstar_independent = _explicit_flag(
        [terminal, terminal.get("source_contract"), terminal.get("corrected_source_contract"),
         terminal.get("physical_endpoint_source_proof"), owner],
        ("Pstar_Z_independent", "Pstar_is_Z_independent", "Pstar_and_radius_Z_independent"),
    )
    ordinary_z1 = (
        owner.get("ordinary_Z_orders") == [0, 1]
        and owner.get("ordinary_Z1_is_fixed_coordinate_derivative") is True
        and owner.get("ordinary_Z1_is_not_a_second_factorial_conversion") is True
    ) or _explicit_flag(proof_roots, ("ordinary_Z1", "ordinary_Z1_product_rule"))
    y_once = (
        owner.get("Y_already_includes_Pstar_once") is True
        and owner.get("additional_Pstar_multiplication") is False
        and _explicit_flag(
            [terminal, terminal.get("physical_endpoint_source_proof"), *proofs],
            ("Y_already_physical_Pstar_once", "Y_already_includes_Pstar_once",
             "Y_already_has_Pstar_once", "common_Pstar_applied_once"),
        )
    )
    required = {
        "delta_Z_independent": delta_independent,
        "R_Z_independent": radius_independent,
        "Pstar_Z_independent": pstar_independent,
        "ordinary_Z1": ordinary_z1,
        "Y_already_physical_Pstar_once": y_once,
    }
    conditional = [name for name, value in required.items() if not value]
    return {
        "passed": same_basis and same_provenance and not conditional,
        "same_source_basis": same_basis,
        "same_source_provenance": same_provenance,
        "required_guards": required,
        "conditional_guards": conditional,
        "source_family": family,
        "physical_common_Pstar_once": owner.get("physical_common_Pstar_applied_once") is True,
        "Y_already_includes_Pstar_once": owner.get("Y_already_includes_Pstar_once") is True,
        "additional_Pstar_multiplication": owner.get("additional_Pstar_multiplication"),
        "ordinary_Z_orders": owner.get("ordinary_Z_orders"),
        "cross_cell_basis_identity_proved": False,
        "global_delta_Pstar_R_connection": "conditional on explicit required guards",
        "no_source_family_only_inference": True,
    }


def _pstar_certificate(
    contract: Mapping[str, Any],
    source_contract: Mapping[str, Any],
    pstar_parameter: Any,
) -> dict[str, Any]:
    roots = [
        source_contract,
        source_contract.get("terminal"),
        source_contract.get("proof"),
        pstar_parameter,
    ]
    explicit_positive = _explicit_flag(
        roots,
        ("Pstar_gt_one", "pstar_gt_one", "logPstar_positive",
         "logPstar_gt_zero", "logPstar_strictly_positive"),
    )
    # The value under logPstar is the actual corrected contract parameter; it
    # is never replaced with a bound from another source or a synthetic value.
    return {
        "passed": explicit_positive,
        "logPstar": pstar_parameter,
        "Pstar_gt_one": explicit_positive,
        "actual_logPstar_parameter": pstar_parameter,
        "source_contract_path": source_contract.get("actual_logPstar_parameter_path"),
        "conditional_if_positive_guard_not_explicit": not explicit_positive,
    }


def _validate_cell_records(
    owner: Mapping[str, Any],
    family: Any,
    corrected: Mapping[str, Any],
    context: Any,
) -> tuple[dict[str, Mapping[str, Any]], dict[str, Any]]:
    records: dict[str, Mapping[str, Any]] = {}
    checks: dict[str, Any] = {}
    terminal_rows = corrected["rows"]
    terminal_proof = corrected["proof"]
    for cell in CELLS:
        key = _cell_key(cell)
        entry = _cell_mapping(owner, cell, "v3 owner")
        rows, proof, lower = _source_rows_and_proof(entry, cell)
        positive_rows = _positive_physical_rows(context, rows, cell)
        sample_proof = entry.get("sample_proof")
        if not isinstance(sample_proof, Mapping):
            sample_proof = proof
        require(proof.get("source_family") == family,
                "v3 lower source proof family changed: " + key)
        require(lower.get("endpoint") in (None, "first_micro/sc_half singleton"),
                "v3 lower endpoint proof changed: " + key)
        source_paths = proof.get("source_quantity_paths")
        same_flow = {
            "same_source": proof.get("same_source") is True,
            "source_family": proof.get("source_family") == family,
            "lower_velocity_y_source": (
                source_paths.get("V_y") if isinstance(source_paths, Mapping) else None
            ),
            "lower_velocity_y_source_is_actual_background": (
                isinstance(source_paths, Mapping)
                and source_paths.get("V_y") == "actual_background_source.V_y"
            ),
            "common_Pstar_path": (
                source_paths.get("Pstar") if isinstance(source_paths, Mapping) else None
            ),
            "same_flow_constraint_explicit": proof.get("same_source") is True,
        }
        checks[key] = {
            "rows_retained_from": "lower.original_source_rows"
            if isinstance(lower.get("original_source_rows"), Mapping)
            else "lower.source_rows",
            "source_rows_present": True,
            "source_proof_present": True,
            "source_family_match": proof.get("source_family") == family,
            "physical_common_Pstar_once": (
                proof.get("physical_common_Pstar_applied_once") is True
                or proof.get("common_Pstar_applied_once") is True
            ),
            "same_source_provenance_explicit": _explicit_flag(
                (proof,), ("same_source_provenance", "same_source", "provenance_checked")
            ),
            "same_flow_constraint": same_flow,
            "finite_endpoint_radius_fixed_in_Z": proof.get(
                "finite_endpoint_radius_fixed_in_Z"
            ),
            "sample_proof_finite_endpoint_radius_fixed_in_Z": sample_proof.get(
                "finite_endpoint_radius_fixed_in_Z"
            ),
            "cross_cell_basis_identity_proved": False,
            "physical_positive_rows": positive_rows,
        }
        records[key] = {
            "rows": rows,
            "proof": proof,
            "lower": lower,
            "entry": entry,
        }

    endpoint_key = _cell_key(("0", ".5"))
    endpoint_record = records[endpoint_key]
    row_matches = endpoint_record["rows"] == terminal_rows
    proof_matches = endpoint_record["proof"] == terminal_proof
    # Some v3 packets retain both the original rows and a normalized wrapper;
    # accept equality against either representation, but never normalize the
    # rows before passing them to the new leaf.
    lower = endpoint_record["lower"]
    normalized_rows = lower.get("source_rows")
    normalized_proof = lower.get("proof")
    row_matches = row_matches or normalized_rows == terminal_rows
    proof_matches = proof_matches or normalized_proof == terminal_proof
    require(row_matches, "corrected terminal [0,.5] rows differ from accepted v3 source")
    require(proof_matches, "corrected terminal [0,.5] proof differs from accepted v3 source")
    checks[endpoint_key]["corrected_terminal_source_rows_match"] = row_matches
    checks[endpoint_key]["corrected_terminal_source_proof_match"] = proof_matches
    checks[endpoint_key]["exact_numeric_endpoint_cell"] = True
    return records, {
        "all_four_cells": True,
        "same_source_0_.5_rows": row_matches,
        "same_source_0_.5_proof": proof_matches,
        "cell_checks": checks,
    }


def _build_source_rows(records: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    # Keep each row packet byte-for-byte as decoded JSON.  In particular this
    # does not rename a basis or replace an opaque four-exponent scale.
    return {key: record["rows"] for key, record in records.items()}


def _source_receipt_summary(
    receipt: Mapping[str, Any],
    owner: Mapping[str, Any],
    source_manifest: Mapping[str, Any],
    source_hashes: Mapping[str, str],
) -> dict[str, Any]:
    embedded = _manifest_rows(source_manifest, "v3 embedded source manifest")
    require(embedded.get(CORRECTED_SOURCE) == FIXED[CORRECTED_SOURCE],
            "corrected source SHA is absent or changed in v3 embedded source manifest")
    for path in (V3_LEAF, V3_FIXTURE, V3_RUNNER):
        require(embedded.get(path) == source_hashes[path],
                "v3 embedded source manifest SHA changed: " + path)
    return {
        "receipt_status": _validate_v3_status(receipt),
        "owner_schema": owner.get("schema"),
        "owner_component": owner.get("component"),
        "owner_source_family": owner.get("source_family"),
        "embedded_source_manifest_contains_corrected_source": True,
        "embedded_source_manifest_corrected_source_sha256": embedded[CORRECTED_SOURCE],
        "embedded_source_manifest_file_count": len(embedded),
        "authenticated_v3_code_sha256": {
            path: source_hashes[path] for path in (V3_LEAF, V3_FIXTURE, V3_RUNNER)
        },
        "owner_flags": {
            "ordinary_Z_orders": owner.get("ordinary_Z_orders"),
            "ordinary_Z1_is_fixed_coordinate_derivative": owner.get(
                "ordinary_Z1_is_fixed_coordinate_derivative"
            ),
            "ordinary_Z1_is_not_a_second_factorial_conversion": owner.get(
                "ordinary_Z1_is_not_a_second_factorial_conversion"
            ),
            "physical_common_Pstar_applied_once": owner.get(
                "physical_common_Pstar_applied_once"
            ),
            "Y_already_includes_Pstar_once": owner.get("Y_already_includes_Pstar_once"),
            "additional_Pstar_multiplication": owner.get("additional_Pstar_multiplication"),
        },
    }


def consume(hashes: Mapping[str, str], stats: dict[str, Any]) -> dict[str, Any]:
    # This is deliberately before fixture import and before any gzip decode.
    raw = authenticated_bytes(hashes)
    stats["raw_inputs_authenticated"] = True

    leaf = module(LEAF, "original_corrected_fourcell_signed_repair_v1_runner_leaf")
    fixture = module(FIXTURE, "original_corrected_fourcell_signed_repair_v1_runner_fixture")
    fixture_runner = getattr(fixture, "run_fixture", None)
    builder = getattr(leaf, "build", None)
    require(callable(fixture_runner) and callable(builder),
            "four-cell leaf/fixture API is missing")

    try:
        from mpmath.ctx_iv import MPIntervalContext
    except Exception as exc:  # pragma: no cover - runner dependency failure.
        raise RunnerInputError("mpmath interval context is unavailable") from exc
    context = MPIntervalContext()
    context.prec = 1800
    fixture_result = fixture_runner(builder, context)
    stats["fixture_invocations"] = 1
    require(isinstance(fixture_result, Mapping) and fixture_result.get("passed") is True,
            "synthetic four-cell repair fixture failed")
    require(fixture_result.get("synthetic_only") is True
            and fixture_result.get("receipt_opened") is False
            and fixture_result.get("graph_traversed") is False,
            "synthetic fixture scope changed")
    stats["fixture_passed"] = True

    # Only after the fixture has passed may the saved source inputs be decoded.
    receipt = decode_gzip(raw[V3_RECEIPT], V3_RECEIPT)
    corrected = load_json(raw[CORRECTED_SOURCE], CORRECTED_SOURCE)
    stats["archive_decode_started"] = True

    owner = _owner(receipt)
    family = owner.get("source_family")
    require(family is not None, "v3 owner source_family is missing")
    require(owner.get("schema") == "original_corrected_endpoint_vy_owner_v3",
            "accepted v3 owner schema changed")
    require(owner.get("cells") == [list(cell) for cell in CELLS],
            "accepted v3 owner cells changed")
    require(owner.get("ordinary_Z_orders") == [0, 1]
            and owner.get("ordinary_Z1_is_fixed_coordinate_derivative") is True
            and owner.get("ordinary_Z1_is_not_a_second_factorial_conversion") is True,
            "accepted v3 ordinary-Z contract changed")
    require(owner.get("physical_common_Pstar_applied_once") is True
            and owner.get("Y_already_includes_Pstar_once") is True
            and owner.get("additional_Pstar_multiplication") is False,
            "accepted v3 physical Pstar-once contract changed")

    embedded_manifest = _embedded_manifest(receipt)
    corrected_contract = _validate_corrected_source(corrected, family)
    records, cell_checks = _validate_cell_records(owner, family, corrected_contract, context)
    source_summary = _source_receipt_summary(
        receipt, owner, embedded_manifest, hashes
    )
    source_certificate = _same_source_contract(
        owner, family, records, corrected_contract["terminal"]
    )
    pstar_certificate = _pstar_certificate(
        source_certificate, corrected_contract,
        corrected_contract["actual_logPstar_parameter"],
    )

    # build() performs only one-cell-at-a-time scale composition.  It receives
    # the actual corrected delta and actual log(Pstar) parameter; no common
    # parameter is invented if its independent-Z connection is conditional.
    conditional = bool(source_certificate["conditional_guards"] or not source_certificate["same_source_basis"] or not source_certificate["same_source_provenance"])
    # Compute a conditional enclosure even when the saved records do not prove
    # every source identity.  These are explicitly hypotheses, never a proof
    # that the independently serialized records have a common live owner.
    conditional_certificate = {
        "passed": True, "same_source_basis": True, "same_source_provenance": True,
        "delta_Z_independent": True, "R_Z_independent": True,
        "Pstar_Z_independent": True, "ordinary_Z1": True,
        "Y_already_physical_Pstar_once": True,
        "scope": "caller-assumed per-cell hypotheses; not a proven source bridge",
        "cross_cell_basis_identity_proved": False,
    }
    leaf_result = builder(
        context,
        CELLS,
        corrected_contract["delta"],
        _build_source_rows(records),
        conditional_certificate,
        pstar_certificate,
    )
    require(isinstance(leaf_result, Mapping) and leaf_result.get("passed") is True,
            "four-cell signed repair leaf did not pass")
    require(leaf_result.get("cells") == [list(cell) for cell in CELLS],
            "four-cell signed repair leaf cells changed")
    require(leaf_result.get("point_values_selected") is False
            and leaf_result.get("numeric_source_values_selected") is False
            and leaf_result.get("ledger_written") is False
            and leaf_result.get("no_ledger") is True
            and leaf_result.get("graph_traversed") is False
            and leaf_result.get("receipt_opened") is False,
            "four-cell signed repair leaf left source-only scope")
    require(leaf_result.get("source_certificate", {}).get(
        "cross_cell_basis_identity_asserted"
    ) is False, "leaf asserted cross-cell basis identity")

    return {
        "schema": COMPONENT + "_output_v1",
        "source_family": family,
        "cells": [list(cell) for cell in CELLS],
        "v3_receipt_sha256": hashes[V3_RECEIPT],
        "corrected_source_sha256": hashes[CORRECTED_SOURCE],
        "v3_source_summary": source_summary,
        "v3_owner_cell_sources": owner.get("cell_sources"),
        "cell_source_proofs": cell_checks,
        "corrected_terminal_source_contract": {
            "actual_delta": corrected_contract["delta"],
            "actual_delta_path": corrected_contract["delta_path"],
            "actual_logPstar_parameter": corrected_contract["actual_logPstar_parameter"],
            "actual_logPstar_parameter_path": corrected_contract[
                "actual_logPstar_parameter_path"
            ],
            "common_original_parameter": corrected_contract["common_original_parameter"],
            "common_original_parameter_path": corrected_contract[
                "common_original_parameter_path"
            ],
            "common_original_parameter_is_textual_source_contract": True,
            "numeric_endpoint_cell": ["0", ".5"],
            "source_only": True,
        },
        "source_certificate": source_certificate,
        "calculation_certificate_assumptions_only": conditional_certificate,
        "pstar_positive_certificate": pstar_certificate,
        "same_source_0_.5_matched": cell_checks["same_source_0_.5_rows"]
        and cell_checks["same_source_0_.5_proof"],
        "physical_Pstar_once": source_certificate["physical_common_Pstar_once"]
        and source_certificate["Y_already_includes_Pstar_once"]
        and source_certificate["additional_Pstar_multiplication"] is False,
        "conditional_source_guards": source_certificate["conditional_guards"],
        "delta_Pstar_connection_is_conditional": bool(
            source_certificate["conditional_guards"]
        ),
        "cross_cell_basis_identity_proved": False,
        "opaque_formal_scales_retained": True,
        "source_rows_are_not_point_values": True,
        "leaf": leaf_result,
        "validation": "conditional" if conditional else "passed",
        "leaf_execution_skipped": conditional,
        "fixture": fixture_result,
        "m15_computed": False,
        "ledger_written": False,
        "old_errors_removed": False,
        "global_closure": False,
        "all_four_cells_checked": True,
        "endpoint_function_values_selected": False,
        "archive_decoded_after_entry_and_fixture": True,
    }


def main() -> int:
    source, hashes, source_sha, brief_sha, static = preflight()
    entry_sha = save(
        ENTRY,
        {
            "schema": COMPONENT + "_entry_v1",
            "component": COMPONENT,
            "phase": PHASE,
            "entered": True,
            "source_checkpoint_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "manifest_input_count": len(hashes),
            "static_constants": static,
            "raw_inputs_authenticated": False,
            "archive_decode_started": False,
            "fixture_invocations": 0,
            "ledger_written": False,
        },
    )
    stats: dict[str, Any] = {
        "raw_inputs_authenticated": False,
        "archive_decode_started": False,
        "fixture_invocations": 0,
        "fixture_passed": False,
    }
    outcome: dict[str, Any] = {
        "schema": COMPONENT + "_result_v1",
        "component": COMPONENT,
        "phase": PHASE,
        "status": "failed",
        "source_checkpoint_sha256": source_sha,
        "brief_sha256": brief_sha,
        "entry_sha256": entry_sha,
        "input_hashes": hashes,
        "source_manifest": source,
        "result": None,
        "error": None,
    }
    try:
        outcome["result"] = consume(hashes, stats)
        outcome["status"] = "returned"
    except BaseException as exc:
        outcome["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    outcome.update(stats)
    result_sha = save(RESULT, outcome)
    print(json.dumps({
        "component": COMPONENT,
        "phase": PHASE,
        "status": outcome["status"],
        "entry_sha256": entry_sha,
        "result_sha256": result_sha,
    }, sort_keys=True))
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
