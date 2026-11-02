"""Run the frozen Cartesian fixed-Z pullback leaf exactly once.

The prior row-5 receipt is authenticated by raw bytes only.  After the
exclusive entry is written, this runner reads the two original source files
for static coordinate labels, imports only the new leaf and fixture, and
invokes the fixture-owned evaluation once.  No old source consumer or
receipt decoder is called.
"""

from __future__ import annotations

import hashlib
import importlib.util
import inspect
import json
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from typing import Any, Callable

from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "cartesian_z_pullback_v1"
DIR = ROOT / "tmp" / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entry.json", "result.json"))

LEAF = "src/openai_ns_reconstruction/cartesian_z_pullback_v1.py"
FIXTURE = "src/openai_ns_reconstruction/cartesian_z_pullback_fixture_v1.py"
RUNNER = "scripts/run_cartesian_z_pullback_v1_once.py"
CARTESIAN_SOURCE = "src/openai_ns_reconstruction/original_n1_cartesian_core_source.py"
PHYSICAL_SOURCE = "src/openai_ns_reconstruction/original_n1_physical_source.py"
PREVIOUS_RECEIPT = "experiments/root_st073/row5_firstword_product_v4_receipt.json.gz"
REQUIRED = (
    LEAF,
    FIXTURE,
    RUNNER,
    CARTESIAN_SOURCE,
    PHYSICAL_SOURCE,
    PREVIOUS_RECEIPT,
)

PIN = {
    CARTESIAN_SOURCE: "bc65a55d79df9812c404a856c8afd0d3a250c0ec6511ea87f50dc80f4546199e",
    PHYSICAL_SOURCE: "542178f6a71dd20e7091dc5f76f115c78274717fd4cce5be391e4c2dbd1158bb",
    PREVIOUS_RECEIPT: "1db10d5abf6b619297d1dd0be282e7dad845b9b5fdf0ecdff000b15693fccfc1",
}

COORDINATE_FORMULA = {
    "lambda_squared": "tau/(1-Z^2)",
    "R": "(X^2+Y^2)/(2*lambda^2)",
    "rho": "R/epsilon",
    "z_phys": "lambda^(1-delta)*Z",
    "row5_current_x": "sqrt(rho), normalized; not Cartesian X",
}

STATIC_LABELS = {
    CARTESIAN_SOURCE: (
        "class OriginalN1CartesianCoreSource",
        "rho=(x*x+y*y)/(2*lambda**2*epsilon)<=4",
        "g.add(1, g.neg(g.mul(zeta, zeta)))",
        'g.unary("sqrt", g.quotient(tau_node, d, "tau>0 and |Z|<1"))',
        "R = g.quotient(g.add(g.mul(x, x), g.mul(y, y)), g.mul(g.constant(2), lam, lam), \"lambda>0\")",
        "rho = g.quotient(R, self.epsilon, \"same original epsilon>0\")",
        "z=g.mul(power(g.add(1, g.neg(de))), zeta)",
        "Z_independent_substitution=False",
        "physical_coordinate_pullback=True",
        "physical_xy_coordinates=True",
        "physical_z_inverse_available=False",
    ),
    PHYSICAL_SOURCE: (
        "class OriginalN1PhysicalSource",
        "lambda_source=lam",
        "coordinates=dict(x=g.mul(r,cos),y=g.mul(r,sin),z=z",
        "parametric_physical_coordinates=True",
        "inverse_Cartesian_chart_selection_available=False",
    ),
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
    raw = (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")
    with path.open("xb") as handle:
        handle.write(raw)
    return sha(raw)


def preflight() -> tuple[dict[str, str], str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    require(
        brief.get("component") == COMPONENT
        and brief.get("phase") == COMPONENT + "_pending"
        and brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False
        and brief.get("one_entry_only") is True
        and brief.get("do_not_reenter") is True,
        "brief guard",
    )
    source_raw = SOURCE.read_bytes()
    checkpoint = brief["source_checkpoint"]
    require(Path(checkpoint["path"]).resolve() == SOURCE.resolve(), "source checkpoint path")
    require(checkpoint["sha256"] == sha(source_raw), "source checkpoint hash")
    require(
        Path(checkpoint["preimage_path"]).resolve()
        == Path("C:/Users/z5242/Documents/Codex/cartesian_z_pre.json").resolve(),
        "preimage path",
    )
    require(sha_path(Path(checkpoint["preimage_path"])) == checkpoint["preimage_sha256"], "preimage hash")
    manifest = json.loads(source_raw)
    require(manifest.get("component") == COMPONENT and manifest.get("one_entry_only") is True, "manifest guard")
    rows = manifest.get("files")
    require(isinstance(rows, list) and len(rows) == len(REQUIRED), "manifest file count")
    hashes: dict[str, str] = {}
    for row in rows:
        require(isinstance(row, Mapping), "manifest row")
        path = row.get("path")
        require(isinstance(path, str) and path in REQUIRED and path not in hashes, "unexpected manifest input")
        digest = sha_path(ROOT / path)
        require(digest == row.get("sha256"), "changed input: " + path)
        if path in PIN:
            require(digest == PIN[path], "changed pinned input: " + path)
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), "missing manifest input")
    return hashes, sha(source_raw), sha(brief_raw)


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
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return {"unserializable_type": type(value).__name__, "repr": repr(value)}


def module(path: str, name: str) -> Any:
    require(path in (LEAF, FIXTURE), "unexpected import path: " + path)
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    require(spec is not None and spec.loader is not None, "new component unavailable: " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def bind_source_labels(raw: Mapping[str, bytes]) -> dict[str, Any]:
    """Read only frozen source text and verify the exact coordinate labels."""

    matched: dict[str, list[str]] = {}
    for path, labels in STATIC_LABELS.items():
        text = raw[path].decode("utf-8")
        missing = [label for label in labels if label not in text]
        require(not missing, "source coordinate label missing: " + path + ": " + repr(missing))
        matched[path] = list(labels)
    return {
        "status": "matched",
        "read_after_exclusive_entry": True,
        "static_text_only": True,
        "labels": matched,
        "old_source_not_imported": True,
    }


def _weight_power_value(owner: Any) -> tuple[bool, Any]:
    for name in ("WEIGHT_POWER", "DEFAULT_WEIGHT_POWER", "weight_power"):
        if hasattr(owner, name):
            return True, getattr(owner, name)
    return False, None


def _invoke_with_optional_weight_power(function: Callable[..., Any], args: tuple[Any, ...], owner: Any) -> tuple[Any, Any]:
    """Call once, adding weight_power only when the inspected signature requires it."""

    try:
        parameters = inspect.signature(function).parameters
    except (TypeError, ValueError):
        return function(*args), "signature_unavailable_default_call"
    parameter = parameters.get("weight_power")
    if parameter is None or parameter.default is not inspect.Parameter.empty:
        return function(*args), "fixture_default"
    found, value = _weight_power_value(owner)
    require(found, "required weight_power has no frozen owner value")
    return function(*args, weight_power=value), value


def run_fixture_once(leaf: Any, fixture: Any, context: MPIntervalContext) -> tuple[Any, Any]:
    run = getattr(fixture, "run_fixture", None)
    require(callable(run), "fixture.run_fixture is required")
    return _invoke_with_optional_weight_power(run, (leaf.evaluate, context), fixture)


def run_optional_stencil(leaf: Any, fixture: Any, context: MPIntervalContext) -> dict[str, Any]:
    stencil = getattr(leaf, "stencil", None)
    stencil_runner = getattr(fixture, "run_stencil", None)
    if not callable(stencil) and not callable(stencil_runner):
        return {"present": False, "called": False, "status": "not_present"}
    require(callable(stencil), "fixture exposes stencil runner but leaf.stencil is absent")
    if callable(stencil_runner):
        value, weight_power = _invoke_with_optional_weight_power(stencil_runner, (stencil, context), fixture)
        return {
            "present": True,
            "called": True,
            "runner": "fixture.run_stencil",
            "weight_power": weight_power,
            "result": encode(value),
        }
    try:
        parameters = inspect.signature(stencil).parameters
    except (TypeError, ValueError) as exc:
        raise TypeError("leaf.stencil requires fixture.run_stencil") from exc
    required = [
        parameter
        for parameter in parameters.values()
        if parameter.kind in (parameter.POSITIONAL_ONLY, parameter.POSITIONAL_OR_KEYWORD, parameter.KEYWORD_ONLY)
        and parameter.default is inspect.Parameter.empty
    ]
    require(len(required) == 1, "leaf.stencil requires fixture.run_stencil")
    params = {"Z": "symbolic_Z", "rho": "symbolic_rho", "lambda": "symbolic_lambda",
              "delta": "symbolic_delta", "fixed_conditions": {
                  "X": True, "Y": True, "tau": True, "epsilon": True, "delta": True}}
    return {"present": True, "called": True, "runner": "leaf.stencil",
            "symbolic_only": True, "result": encode(stencil(params))}


def consume(hashes: Mapping[str, str]) -> dict[str, Any]:
    output: dict[str, Any] = {
        "coordinate_formula": COORDINATE_FORMULA,
        "source_pins": {path: hashes[path] for path in REQUIRED},
        "raw_receipt_sha256_only": {PREVIOUS_RECEIPT: hashes[PREVIOUS_RECEIPT]},
        "source_binding": None,
        "fixture": None,
        "stencil": None,
        "missing_derivative_inputs": [
            "f_rho",
            "f_rho_Z",
            "f_rho_rho",
            "physical_h5",
            "Cartesian X/Y derivative jets unless explicitly supplied by a caller",
        ],
        "owner_status": {
            "callable_owner": "unknown",
            "source_expression_owner": "unknown",
            "frame_owner": "unknown",
            "gauge_owner": "unknown",
            "state_owner": "unknown",
            "physical_h5_owner": "unknown",
        },
        "claims": {
            "row5_x_is_cartesian_X": False,
            "navier_stokes_solution": False,
            "m15": False,
            "ledger_appended": False,
            "unproved_cancellation_used": False,
        },
    }
    try:
        raw_sources: dict[str, bytes] = {}
        for path in REQUIRED:
            input_path = ROOT / path
            if path == PREVIOUS_RECEIPT:
                require(sha_path(input_path) == hashes[path], "changed after entry: " + path)
                continue
            if path in (CARTESIAN_SOURCE, PHYSICAL_SOURCE):
                data = input_path.read_bytes()
                require(sha(data) == hashes[path], "changed after entry: " + path)
                raw_sources[path] = data
            else:
                require(sha_path(input_path) == hashes[path], "changed after entry: " + path)

        output["source_binding"] = bind_source_labels(raw_sources)
        leaf = module(LEAF, "cartesian_z_pullback_v1_new_leaf")
        fixture = module(FIXTURE, "cartesian_z_pullback_fixture_v1_new_fixture")
        require(callable(getattr(leaf, "evaluate", None)), "leaf.evaluate is required")
        context = MPIntervalContext()
        context.prec = 1800
        proof, weight_power = run_fixture_once(leaf, fixture, context)
        encoded_proof = encode(proof)
        output["fixture"] = {
            "call": "fixture.run_fixture(leaf.evaluate, c)",
            "weight_power": weight_power,
            "complete_return": encoded_proof,
            "checks": encoded_proof.get("checks") if isinstance(encoded_proof, Mapping) else None,
            "interval_precision": 1800,
        }
        require(isinstance(proof, Mapping), "fixture result must be a mapping")
        require(proof.get("passed") is True, "fixture checks failed")

        output["stencil"] = run_optional_stencil(leaf, fixture, context)
        if output["stencil"].get("called"):
            stencil_result = output["stencil"].get("result")
            if isinstance(stencil_result, Mapping) and "passed" in stencil_result:
                require(stencil_result["passed"] is True, "stencil checks failed")
        output["status"] = "returned"
        output["result_status"] = "fixture_and_optional_stencil_complete"
    except BaseException as exc:
        output["status"] = "failed"
        output["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    return output


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
        },
    )
    output = consume(hashes)
    result = {
        "component": COMPONENT,
        "source_sha256": source_sha256,
        "entry_sha256": entry_sha256,
        "status": output["status"],
        "output": output,
    }
    result_sha256 = save(RESULT, result)
    print(json.dumps({"component": COMPONENT, "status": result["status"], "result_sha256": result_sha256}))
    return 0 if result["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
