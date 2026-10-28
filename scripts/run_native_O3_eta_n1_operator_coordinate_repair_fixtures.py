"""Single-use coordinate-repair fixture; source-only until the root runs it.

This leaf is deliberately independent of the entered Green leaf.  It never
imports the old runner, the old operator, or any Green export.  After the
exclusive ``entered.json`` write it imports only the coordinate-repair
operator, the already-published response interface for ``bind``/``evaluate``
only, and the not-yet-entered finite-word module.  All W/W_Z values used by
the m15 check are explicit synthetic raw tuples, not Green output.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import traceback
import uuid
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any, Callable, Mapping


SCHEMA = "native_O3_eta_n1_operator_coordinate_repair_fixture_v1"
RUNTIME_MODEL_METADATA = None
FIXTURE_KIND = "synthetic_coordinate_repair_fixture_only"
OLD212_UNTOUCHED = True
REQUIRED_COMPONENT = "native_O3_eta_n1_operator_coordinate_repair"
REQUIRED_PHASE = "prepared_unexecuted_operator_coordinate_repair"
SOURCE_REF_KEY = "operator_coordinate_repair_source_checkpoint"
PENDING_B_KEY = "pending_B"

R_ROOT = Path(r"C:/Users/z5242/Documents/Codex/2026-09-13/before-doing-any-project-work-verify/work/st073-transition-next")
W_ROOT = Path(r"C:/Users/z5242/Documents/Codex/2026-10-16/https-github-com-liang-techh-openai/work")
CHECKPOINT_DIR = R_ROOT / "tmp" / "native_O3_eta_n1_operator_coordinate_repair_v1"
SOURCE_MANIFEST = CHECKPOINT_DIR / "source.json"
BRIEF_DEFAULT = CHECKPOINT_DIR / "brief_prepared.json"
ENTERED_RECORD = CHECKPOINT_DIR / "entered.json"
RESULT_RECORD = CHECKPOINT_DIR / "result.json"

EXPECTED_FILES = {
    "actual_operator_coordinate_repair": R_ROOT / "src" / "openai_ns_reconstruction" / "native_O3_eta_n1_actual_operator_coordinate_repair.py",
    "response_interface_repaired": R_ROOT / "src" / "openai_ns_reconstruction" / "native_O3_eta_n1_response_interface_repaired.py",
    "volterra_words": R_ROOT / "src" / "openai_ns_reconstruction" / "native_O3_eta_n1_volterra_words.py",
    "runner": Path(__file__).resolve(),
}

PRIOR_FAILURES = {
    "operator_bridge_v3": {
        "path": R_ROOT / "tmp" / "native_O3_eta_n1_operator_bridge_v3" / "result.json",
        "sha256": "03687d2c5bd9e5e3d9689e30d46b3a12ce6304243d24016b0a6db9f3e18d5048",
    },
    "operator_bridge_repair_v1": {
        "path": R_ROOT / "tmp" / "native_O3_eta_n1_operator_bridge_repair_v1" / "result.json",
        "sha256": "b0a5194d8cfa78b015e2fee22fc289ea3c30831396c7892d34eb234ac9ed6ff1",
    },
}

ASSUMED_API = {
    "actual_operator": "build_actual_six_row_operator",
    "response_bind": "bind_same_source_n1_response",
    "response_evaluate": "evaluate_n1_m15",
    "volterra_words": "prepare_finite_neumann_words",
    "green_calls": False,
    "old_runner_calls": False,
    "word_contract": {
        "green_scope": "common_core_fixed_x",
        "green_kernel_Z_independent": True,
        "preconditioned_green": True,
    },
}


class FixtureRefusal(RuntimeError):
    pass


class FixtureFailure(RuntimeError):
    pass


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _read(path: Path) -> Any:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def _write_exclusive(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _resolve_r_path(raw: Any) -> Path:
    path = Path(str(raw))
    return (path if path.is_absolute() else R_ROOT / path).resolve()


def _path_keys(path: Path) -> tuple[str, ...]:
    keys: list[str] = []
    for root in (R_ROOT, W_ROOT):
        try:
            keys.append(path.resolve().relative_to(root.resolve()).as_posix())
        except ValueError:
            pass
    keys.extend((path.name, path.resolve().as_posix()))
    return tuple(dict.fromkeys(keys))


def _sha_record(value: Any) -> str | None:
    if isinstance(value, str):
        return value
    if isinstance(value, Mapping) and isinstance(value.get("sha256"), str):
        return value["sha256"]
    return None


def _immutable_reference(references: Any, label: str, path: Path, expected: str) -> dict[str, Any]:
    if not isinstance(references, list):
        raise FixtureRefusal("source.json immutable_references must be a list")
    target = path.resolve()
    for item in references:
        if not isinstance(item, Mapping) or item.get("sha256") != expected:
            continue
        item_path = _resolve_r_path(item.get("path")) if item.get("path") else None
        if item_path != target:
            continue
        if not target.exists() or _sha256(target) != expected:
            raise FixtureRefusal("immutable reference changed: " + target.as_posix())
        return {"label": label, "path": target.as_posix(), "sha256": expected}
    raise FixtureRefusal("missing immutable reference: " + label)


def _pending_b_is_true(manifest: Mapping[str, Any], brief: Mapping[str, Any]) -> bool:
    values = (manifest.get(PENDING_B_KEY), brief.get(PENDING_B_KEY))
    for value in values:
        if value is True:
            return True
        if isinstance(value, Mapping) and value.get("prepared") is True:
            return True
    return False


def _preflight() -> dict[str, Any]:
    if ENTERED_RECORD.exists() or RESULT_RECORD.exists():
        raise FixtureRefusal("single-use refusal: coordinate repair entered/result exists")
    if not SOURCE_MANIFEST.exists():
        raise FixtureRefusal("missing coordinate repair source.json")
    manifest = _read(SOURCE_MANIFEST)
    if not isinstance(manifest, Mapping):
        raise FixtureRefusal("source.json must be an object")
    if manifest.get("prepared") is not True or manifest.get("brief_prepared") is not True:
        raise FixtureRefusal("source.json prepared and brief_prepared must be true")
    brief_raw = manifest.get("brief_path")
    if brief_raw is None and isinstance(manifest.get("brief"), Mapping):
        brief_raw = manifest["brief"].get("path")
    brief = BRIEF_DEFAULT if brief_raw is None else _resolve_r_path(brief_raw)
    if not brief.exists():
        raise FixtureRefusal("missing prepared coordinate repair brief: " + brief.as_posix())
    brief_data = _read(brief)
    if not isinstance(brief_data, Mapping) or brief_data.get("prepared") is not True:
        raise FixtureRefusal("prepared coordinate brief must contain prepared=true")
    if brief_data.get("component") != REQUIRED_COMPONENT:
        raise FixtureRefusal("coordinate repair component mismatch")
    if brief_data.get("phase") not in (REQUIRED_PHASE, "pending_B", "prepared_pending_B"):
        raise FixtureRefusal("coordinate repair phase mismatch")
    if brief_data.get("goal_active") is not True or brief_data.get("full_goal_complete") is not False or brief_data.get("do_not_reenter") is not True:
        raise FixtureRefusal("coordinate repair lifecycle state mismatch")
    if not _pending_b_is_true(manifest, brief_data):
        raise FixtureRefusal("pending_B must be prepared before entry")
    source_ref = brief_data.get(SOURCE_REF_KEY)
    if not isinstance(source_ref, Mapping) or source_ref.get("sha256") != _sha256(SOURCE_MANIFEST):
        raise FixtureRefusal("coordinate repair brief/source SHA mismatch")

    immutable = {
        label: _immutable_reference(manifest.get("immutable_references"), label, data["path"], data["sha256"])
        for label, data in PRIOR_FAILURES.items()
    }
    records = manifest.get("files", manifest.get("code_sha256"))
    if not isinstance(records, Mapping):
        raise FixtureRefusal("source.json needs files/code_sha256 mapping")
    hashes: dict[str, str] = {}
    for label, path in EXPECTED_FILES.items():
        if not path.exists() or not path.is_file():
            raise FixtureRefusal(f"missing {label}: {path}")
        expected = None
        for key in _path_keys(path):
            if key in records:
                expected = _sha_record(records[key])
                if expected is not None:
                    break
        if expected is None:
            raise FixtureRefusal("source.json has no SHA record for " + path.as_posix())
        actual = _sha256(path)
        hashes[label] = actual
        if expected.lower() != actual:
            raise FixtureRefusal(f"SHA mismatch for {label}: expected {expected}, actual {actual}")
    return {
        "manifest_path": SOURCE_MANIFEST.as_posix(),
        "brief_path": brief.as_posix(),
        "source_manifest_schema": manifest.get("schema"),
        "source_hashes": hashes,
        "pending_B": True,
        "prepared_brief": dict(brief_data),
        "prior_failure_immutable_references": immutable,
        "old_runner_executed": False,
        "green_functions_called": False,
    }


def _enter(preflight: Mapping[str, Any]) -> dict[str, Any]:
    entered = {
        "schema": SCHEMA,
        "entered": True,
        "entered_at_utc": _now(),
        "run_id": str(uuid.uuid4()),
        "runtime_model_metadata": RUNTIME_MODEL_METADATA,
        "fixture_kind": FIXTURE_KIND,
        "synthetic_fixture": True,
        "synthetic_raw_m15": True,
        "synthetic_mock_green": True,
        "accepted_as_proof": False,
        "actual_source_identity_bound": False,
        "old212_untouched": OLD212_UNTOUCHED,
        "old_runner_executed": False,
        "green_functions_called": False,
        "preflight": dict(preflight),
        "assumed_api": ASSUMED_API,
    }
    try:
        _write_exclusive(ENTERED_RECORD, entered)
    except FileExistsError as exc:
        raise FixtureRefusal("entered.json appeared during coordinate repair admission") from exc
    return entered


def _load_new(label: str, path: Path) -> ModuleType:
    spec = importlib.util.spec_from_file_location("_native_o3_eta_n1_coordinate_repair_" + label, path)
    if spec is None or spec.loader is None:
        raise FixtureFailure("cannot import coordinate-repair " + label)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _export(module: ModuleType, name: str, label: str) -> Callable[..., Any]:
    value = getattr(module, name, None)
    if not callable(value):
        raise FixtureFailure(f"missing assumed export {label}.{name}")
    return value


def _part(value: Any, name: str) -> Any:
    return value.get(name) if isinstance(value, Mapping) else getattr(value, name, None)


def _ser(value: Any, mp: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(k): _ser(v, mp) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [_ser(v, mp) for v in value]
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    try:
        return mp.nstr(value, 60)
    except Exception:
        return repr(value)


def _check_scalar(mp: Any, actual: Any, expected: Any, tol: Any) -> dict[str, Any]:
    error = abs(actual - expected)
    return {"actual": _ser(actual, mp), "expected": _ser(expected, mp), "absolute_error": _ser(error, mp), "passed": bool(error <= tol)}


def _matrix_expectations(mp: Any, values: Mapping[str, Any]) -> dict[str, Any]:
    F, Fr, U, Ur = values["F"], values["Fr"], values["U"], values["Ur"]
    Q, eps, delta, rho, Z, x = (values[k] for k in ("Q", "epsilon", "delta", "rho", "Z", "x"))
    d = 1 - Z * Z
    L = 1 - delta * Z * Z
    half = mp.mpf(1) / 2
    H = Z * (1 - delta) * half + d * U
    Cf = F + rho * Fr
    return {
        "B0[4,1]": eps * 4 * x * F,
        "B0[6,1]": -eps**2 * 8 * Z * rho * F / L,
        "B1[5,1]": eps * 2 * H / L,
        "B1[5,2]": -eps * 2 * d * Cf / L,
        "B1[6,4]": eps * 2 * d / L,
    }


def _m15_closed(mp: Any, epsilon: Any, rho: Any, Z: Any, ordinary_z1: bool = False) -> Any:
    if ordinary_z1:
        bracket = -(1 + 2 * Z) * rho**2 - 2 * rho**3 / 3 - rho**2 / 4
    else:
        bracket = (2 + Z) * (1 - Z) * rho**2 + 2 * (1 - Z) * rho**3 / 3 - (1 + Z) * rho**2 / 4
    return epsilon * bracket


def _execute_after_entry(entered: Mapping[str, Any]) -> dict[str, Any]:
    import mpmath as mp

    mp.mp.dps = 100
    operator = _load_new("actual_operator_coordinate_repair", EXPECTED_FILES["actual_operator_coordinate_repair"])
    response = _load_new("response_interface_repaired", EXPECTED_FILES["response_interface_repaired"])
    words = _load_new("volterra_words", EXPECTED_FILES["volterra_words"])
    build_operator = _export(operator, ASSUMED_API["actual_operator"], "actual_operator_coordinate_repair")
    bind = _export(response, ASSUMED_API["response_bind"], "response_interface_repaired")
    evaluate = _export(response, ASSUMED_API["response_evaluate"], "response_interface_repaired")
    prepare_words = _export(words, ASSUMED_API["volterra_words"], "volterra_words")
    tol = mp.mpf("1e-80")

    # Independently rounded mpmath coordinate values: do not assert an exact identity.
    x = mp.mpf(3) / mp.mpf(5)
    rho = mp.mpf(9) / mp.mpf(25)
    Z = mp.mpf(1) / mp.mpf(4)
    values = {
        "F": mp.mpf(2), "Fr": mp.mpf(3) / 10, "FZ": mp.mpf(1) / 6,
        "U": mp.mpf(5) / 4, "Ur": -mp.mpf(1) / 8, "UZ": mp.mpf(2) / 9,
        "Q": mp.mpf(7) / 10, "epsilon": mp.mpf(1) / 7, "delta": mp.mpf(1) / 5,
        "rho": rho, "Z": Z, "x": x,
    }
    matrix_value = build_operator(**values)
    coordinate_check = _part(matrix_value, "coordinate_check") or {}
    coordinate_result = {
        "x": _ser(x, mp), "rho": _ser(rho, mp),
        "independently_rounded_mpmath": True,
        "accepted_without_error": True,
        "coordinate_check": _ser(coordinate_check, mp),
        "exact_identity_proof_claimed": any(coordinate_check.get(k) is True for k in ("rho_equals_x_squared_exact", "exact_identity", "exact_proof")),
        "not_exact_proof": not any(coordinate_check.get(k) is True for k in ("rho_equals_x_squared_exact", "exact_identity", "exact_proof")),
    }
    expected = _matrix_expectations(mp, values)
    b0, b1 = _part(matrix_value, "B0"), _part(matrix_value, "B1")
    if not isinstance(b0, Mapping) or not isinstance(b1, Mapping):
        raise FixtureFailure("coordinate-repair operator did not return B0/B1 mappings")
    matrix_checks = {}
    for label, expected_value in expected.items():
        row, key = label.split("[")
        key = key[:-1]
        actual = (b0 if row == "B0" else b1).get(key)
        if actual is None:
            raise FixtureFailure("missing sparse matrix entry " + label)
        matrix_checks[label] = _check_scalar(mp, actual, expected_value, tol)
    matrix_report = {
        "inputs": _ser(values, mp),
        "independent_fixed_rational_data": {key: str(value) for key, value in {"x": "3/5", "rho": "9/25", "Z": "1/4", "epsilon": "1/7", "delta": "1/5", "F": "2", "Fr": "3/10", "U": "5/4", "Ur": "-1/8"}.items()},
        "coordinate": coordinate_result,
        "checks": matrix_checks,
        "passed": coordinate_result["accepted_without_error"] and coordinate_result["not_exact_proof"] and all(row["passed"] for row in matrix_checks.values()),
    }

    # Raw synthetic m15 only: response.bind/evaluate are called, never Green.
    epsilon, Pstar = mp.mpf(1) / 7, mp.mpf(11) / 3
    def raw_W(radial: Any, ordinary_z: Any) -> tuple[Any, ...]:
        return (0, (1 - ordinary_z) * radial, 0, (1 + ordinary_z) * radial / 2, 0, 0)
    def raw_W_Z(radial: Any, ordinary_z: Any) -> tuple[Any, ...]:
        return (0, -radial, 0, radial / 2, 0, 0)
    def integrate(fn: Callable[[Any], Any], lo: Any, hi: Any) -> Any:
        return mp.quad(fn, [lo, hi])
    response_binding = bind(W=raw_W, W_Z=raw_W_Z, Ulead=lambda r, z: 2 + z + r, Ulead_Z=lambda r, z: 1, integrate=integrate, epsilon=epsilon, Pstar=Pstar, scales_Z_independent=True, source_id="coordinate_repair_synthetic_raw_W_not_green", same_source_evidence={"synthetic_fixture": True, "actual_source_identity_bound": False}, axis_conditions={"W_axis_zero": True, "Uz1_axis_zero": True, "P1_axis_normalized": True}, error_data={"synthetic_raw_m15": True, "actual_source_identity_bound": False})
    m15_value = evaluate(response_binding, rho, Z)
    actual_m15, actual_z1 = _part(m15_value, "m15"), _part(m15_value, "ordinary_Z1")
    expected_m15, expected_z1 = _m15_closed(mp, epsilon, rho, Z), _m15_closed(mp, epsilon, rho, Z, True)
    m15_report = {
        "synthetic_raw_W_not_green": True, "not_actual_source": True,
        "W": "(0,(1-Z)rho,0,(1+Z)rho/2,0,0)", "W_Z": "(0,-rho,0,rho/2,0,0)",
        "inputs": {"rho": _ser(rho, mp), "Z": _ser(Z, mp), "epsilon": _ser(epsilon, mp), "Pstar": _ser(Pstar, mp), "Ulead": "2+Z+rho", "Ulead_Z": "1", "measure": "d rho"},
        "m15": _check_scalar(mp, actual_m15, expected_m15, tol),
        "ordinary_Z1": _check_scalar(mp, actual_z1, expected_z1, tol),
    }
    m15_report["passed"] = bool(m15_report["m15"]["passed"] and m15_report["ordinary_Z1"]["passed"])

    word_Z = mp.mpf(2) / 5
    zero = [lambda s, z: 0, lambda s, z: 0, lambda s, z: 0]
    g6 = [lambda s, z: z * z, lambda s, z: 2 * z, lambda s, z: 2]
    c0 = [lambda s, z: {"6,2": 1}, lambda s, z: {}, lambda s, z: {}]
    b1 = [lambda s, z: {"6,2": z}, lambda s, z: {"6,2": 1}, lambda s, z: {}]
    def mock_green(*, x: Any, Z: Any, g4: Callable[..., Any], g5: Callable[..., Any], g6: Callable[..., Any], integrate: Callable[..., Any]) -> tuple[Any, ...]:
        return (0, g6(x, Z), 0, 0, 0, 0)
    word_kwargs = {"x": mp.mpf(1) / 2, "Z": word_Z, "K": 1, "d": 1, "g4_jets": zero, "g5_jets": zero, "g6_jets": g6, "C0_jets": c0, "B1_jets": b1, "principal_green_response": mock_green, "integrate": lambda fn, lo, hi: 0, **ASSUMED_API["word_contract"]}
    word_eval = prepare_words(**word_kwargs).evaluate()
    terms, partial = _part(word_eval, "terms"), _part(word_eval, "partial_sum")
    actual_words = ((_part(terms[0], "jet_values")[0][1], _part(terms[0], "jet_values")[1][1]), (_part(terms[1], "jet_values")[0][1], _part(terms[1], "jet_values")[1][1]), (partial[0][1], partial[1][1]))
    expected_words = ((word_Z**2, 2 * word_Z), (3 * word_Z**2, 6 * word_Z), (4 * word_Z**2, 8 * word_Z))
    word_checks = {name: {"actual": _ser(got, mp), "expected": _ser(want, mp), "passed": all(bool(abs(a - b) <= tol) for a, b in zip(got, want))} for name, got, want in zip(("word0", "word1", "partial_sum"), actual_words, expected_words)}
    missing_rejected, missing_error = False, None
    try:
        prepare_words(**dict(word_kwargs, g6_jets=g6[:2], C0_jets=c0[:2], B1_jets=b1[:2]))
    except Exception as exc:
        missing_rejected, missing_error = True, {"type": type(exc).__name__, "message": str(exc)}
    word_report = {"synthetic_mock_green": True, "not_actual_physical_green": True, "K": 1, "d": 1, "required_jet_order": 2, "g6_jets": ["Z^2", "2Z", "2"], "C0_6_2": "1", "B1_6_2": "Z", "B1_Z1_6_2": "1", "checks": word_checks, "missing_Z2_rejected_before_evaluate": missing_rejected, "missing_input_error": missing_error, "passed": all(item["passed"] for item in word_checks.values()) and missing_rejected}

    passed = bool(matrix_report["passed"] and m15_report["passed"] and word_report["passed"])
    return {"schema": SCHEMA, "runtime_model_metadata": RUNTIME_MODEL_METADATA, "fixture_kind": FIXTURE_KIND, "synthetic_fixture": True, "synthetic_raw_m15": True, "synthetic_mock_green": True, "accepted_as_proof": False, "actual_source_identity_bound": False, "old212_untouched": OLD212_UNTOUCHED, "old_runner_executed": False, "green_functions_called": False, "full_goal_complete": False, "passed": passed, "operator_coordinate_matrix": matrix_report, "m15": m15_report, "finite_volterra_words": word_report, "state": {"coordinate_operator_called_after_entry": True, "response_bind_evaluate_only": True, "green_functions_called": False, "old_runner_executed": False, "old_helper_called": False, "actual_source_value_selected": False, "actual_source_identity_bound": False}, "prior_failure_immutable_references": entered["preflight"]["prior_failure_immutable_references"], "entered": dict(entered)}


def _failure(entered: Mapping[str, Any], exc: BaseException) -> dict[str, Any]:
    return {"schema": SCHEMA, "runtime_model_metadata": RUNTIME_MODEL_METADATA, "fixture_kind": FIXTURE_KIND, "synthetic_fixture": True, "synthetic_raw_m15": True, "synthetic_mock_green": True, "accepted_as_proof": False, "actual_source_identity_bound": False, "old212_untouched": OLD212_UNTOUCHED, "old_runner_executed": False, "green_functions_called": False, "full_goal_complete": False, "passed": False, "failure_stage": "execution", "error_type": type(exc).__name__, "error": str(exc), "traceback": traceback.format_exc(), "state": {"green_functions_called": False, "old_runner_executed": False, "old_helper_called": False, "actual_source_identity_bound": False}, "prior_failure_immutable_references": entered.get("preflight", {}).get("prior_failure_immutable_references", {}), "entered": dict(entered)}


def main() -> int:
    try:
        preflight = _preflight()
        entered = _enter(preflight)
    except Exception as exc:
        print(json.dumps({"passed": False, "stage": "preflight_or_enter", "error": str(exc), "runtime_model_metadata": None}, ensure_ascii=False))
        return 2
    try:
        result = _execute_after_entry(entered)
    except Exception as exc:
        result = _failure(entered, exc)
    try:
        _write_exclusive(RESULT_RECORD, result)
    except Exception as exc:
        print(json.dumps({"passed": False, "stage": "result_write", "error": str(exc), "entered_path": ENTERED_RECORD.as_posix()}, ensure_ascii=False))
        return 3
    print(json.dumps({"passed": bool(result.get("passed")), "stage": "complete", "result_path": RESULT_RECORD.as_posix(), "fixture_kind": FIXTURE_KIND, "green_functions_called": False, "old_runner_executed": False, "runtime_model_metadata": None}, ensure_ascii=False))
    return 0 if result.get("passed") else 1


if __name__ == "__main__":
    raise SystemExit(main())
