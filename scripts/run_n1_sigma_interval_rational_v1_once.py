"""Enter the frozen ``n1_sigma_interval_rational_v1`` leaf exactly once.

This launcher is deliberately source-local.  Its only project files are the
five paths frozen by ``tmp/n1_sigma_interval_rational_v1/source.json``:

* the new leaf module;
* its independent synthetic fixture;
* this launcher; and
* the outer-initial and original cutoff-definition sources, which are
  authenticated but never imported or evaluated here.

The continuation brief and the source manifest are checked before an entry is
created.  The entry is then created with exclusive mode before either new
module is imported.  The fixture is the sole caller and is invoked once as
``run_fixture(module)``.  A returned value or a complete traceback is written
to ``result.json``.  A fixture result is synthetic evidence only; it does not
claim a live source owner or a physical result.

The runner never reads a receipt or gzip artifact, calls an ancestor project
calculation, or accesses a live graph.
"""

from __future__ import annotations

import contextlib
from fractions import Fraction
import hashlib
import importlib.util
import io
import json
import sys
import traceback
from collections.abc import Mapping
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
NAME = "n1_sigma_interval_rational_v1"
PHASE = NAME + "_pending"
DIR = ROOT / "tmp" / NAME
SOURCE = DIR / "source.json"
ENTRY = DIR / "entry.json"
RESULT = DIR / "result.json"

MODULE = "src/openai_ns_reconstruction/n1_sigma_interval_rational_v1.py"
FIXTURE = "src/openai_ns_reconstruction/n1_sigma_interval_rational_v1_fixture.py"
RUNNER = "scripts/run_n1_sigma_interval_rational_v1_once.py"
CUTOFF = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions.py"
)
OUTER_INITIAL = "experiments/root_st073/lei_ren_part1_paper_compliant_outer_initial.py"
FROZEN_FILES = (MODULE, FIXTURE, RUNNER, OUTER_INITIAL, CUTOFF)
SOURCE_DEFINITIONS = (OUTER_INITIAL, CUTOFF)

FIXTURE_API = "run_fixture(module)"
SYNTHETIC_ONLY = {
    "synthetic_fixture_only": True,
    "real_source_owner_claimed": False,
    "physical_result_claimed": False,
    "live_graph_called": False,
    "ancestor_project_calculation_called": False,
}


class ExistingArtifactError(RuntimeError):
    """The one-entry contract was already consumed."""


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def require(condition: object, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _resolve_pointer(value: object, *, base: Path) -> Path:
    require(isinstance(value, str) and value, "path pointer missing")
    candidate = Path(value)
    if not candidate.is_absolute():
        candidate = base / candidate
    return candidate.resolve()


def _root_file(path: str) -> Path:
    """Resolve a manifest path while preventing an escape from this tree."""

    candidate = Path(path)
    require(not candidate.is_absolute(), "manifest path must be relative: " + path)
    resolved = (ROOT / candidate).resolve()
    try:
        resolved.relative_to(ROOT.resolve())
    except ValueError as exc:
        raise ValueError("manifest path escapes repository: " + path) from exc
    return resolved


def _json(raw: bytes, label: str) -> Any:
    try:
        return json.loads(raw)
    except BaseException as exc:
        raise ValueError("invalid JSON: " + label) from exc


def _load_brief(source_raw: bytes) -> tuple[dict[str, Any], bytes]:
    brief_raw = BRIEF.read_bytes()
    brief = _json(brief_raw, "continuation brief")
    require(isinstance(brief, dict), "continuation brief must be an object")
    require(brief.get("component") == NAME, "brief component guard")
    require(brief.get("phase") == PHASE, "brief pending phase guard")
    require(
        brief.get("goal_active") is True
        and brief.get("full_goal_complete") is False,
        "brief goal guard",
    )
    require(brief.get("one_entry_only") is True, "brief one-entry guard")
    require(brief.get("do_not_reenter") is True, "brief do-not-reenter guard")

    checkpoint = brief.get("source_checkpoint")
    require(isinstance(checkpoint, Mapping), "brief source checkpoint missing")
    checkpoint_path = _resolve_pointer(checkpoint.get("path"), base=ROOT)
    require(checkpoint_path == SOURCE.resolve(), "brief source checkpoint path mismatch")
    require(
        checkpoint.get("sha256") == sha(source_raw),
        "brief source checkpoint hash mismatch",
    )
    return brief, brief_raw


def _manifest(source_raw: bytes) -> tuple[dict[str, str], bytes]:
    source = _json(source_raw, "source manifest")
    require(isinstance(source, dict), "source manifest must be an object")
    require(source.get("component") == NAME, "manifest component guard")
    require(source.get("one_entry_only") is True, "manifest one-entry guard")
    require(source.get("do_not_reenter") is True, "manifest do-not-reenter guard")
    if "phase" in source:
        require(source.get("phase") == PHASE, "manifest pending phase guard")

    contract = source.get("input_contract")
    require(isinstance(contract, Mapping), "manifest input contract missing")
    require(contract.get("module") == MODULE, "module contract mismatch")
    require(contract.get("fixture") == FIXTURE, "fixture contract mismatch")
    require(contract.get("runner") == RUNNER, "runner contract mismatch")
    require(contract.get("outer_initial_source") == OUTER_INITIAL, "outer-initial source contract mismatch")
    require(contract.get("original_cutoff_source") == CUTOFF, "cutoff source contract mismatch")
    require(contract.get("fixture_api") == FIXTURE_API, "fixture API contract mismatch")

    rows = source.get("files")
    require(isinstance(rows, list) and len(rows) == len(FROZEN_FILES), "manifest must freeze exactly five files")
    required = set(FROZEN_FILES)
    hashes: dict[str, str] = {}
    for row in rows:
        require(isinstance(row, Mapping), "manifest file row must be an object")
        path = row.get("path")
        require(isinstance(path, str) and path, "manifest file path missing")
        require(path in required and path not in hashes, "unexpected or duplicate frozen file: " + path)
        require(not path.lower().endswith((".gz", ".gzip")), "gzip input is forbidden")
        digest = sha(_root_file(path).read_bytes())
        require(digest == row.get("sha256"), "manifest hash mismatch: " + path)
        hashes[path] = digest
    require(set(hashes) == required, "manifest file set mismatch")

    return hashes, source_raw


def preflight() -> tuple[dict[str, str], str, str]:
    """Read only the brief, manifest and five frozen file bytes."""

    if ENTRY.exists() or RESULT.exists():
        raise ExistingArtifactError("entry or result already exists")
    require(SOURCE.exists(), "source manifest missing")
    source_raw = SOURCE.read_bytes()
    _, brief_raw = _load_brief(source_raw)
    hashes, _ = _manifest(source_raw)
    return hashes, sha(source_raw), sha(brief_raw)


def _authenticate(hashes: Mapping[str, str]) -> dict[str, bytes]:
    raw: dict[str, bytes] = {}
    for path, expected in hashes.items():
        data = _root_file(path).read_bytes()
        require(sha(data) == expected, "frozen input changed after entry: " + path)
        raw[path] = data
    return raw


def _import_new(path: str, name: str) -> Any:
    """Import exactly one of the two new modules named by the manifest."""

    require(path in (MODULE, FIXTURE), "unexpected module import: " + path)
    spec = importlib.util.spec_from_file_location(name, _root_file(path))
    require(spec is not None and spec.loader is not None, "module spec missing: " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def _json_value(value: Any, seen: set[int] | None = None) -> Any:
    """Serialize fixture output without evaluating project objects."""

    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, Fraction):
        return {
            "__fraction__": True,
            "numerator": str(value.numerator),
            "denominator": str(value.denominator),
        }
    if seen is None:
        seen = set()
    identity = id(value)
    if identity in seen:
        return {"__cycle__": True, "type": type(value).__name__}
    seen.add(identity)
    try:
        if isinstance(value, Mapping):
            return {str(key): _json_value(item, seen) for key, item in value.items()}
        if isinstance(value, (list, tuple)):
            return [_json_value(item, seen) for item in value]
        if isinstance(value, (set, frozenset)):
            return [_json_value(item, seen) for item in sorted(value, key=repr)]
        return {
            "__unsupported_type__": type(value).__name__,
            "repr": repr(value),
        }
    finally:
        seen.discard(identity)


def _error(exc: BaseException) -> dict[str, str]:
    return {
        "type": type(exc).__name__,
        "message": str(exc),
        "traceback": traceback.format_exc(),
    }


def _write_exclusive(path: Path, value: object) -> str:
    raw = (
        json.dumps(
            _json_value(value),
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
        + "\n"
    ).encode("utf-8")
    with path.open("xb") as handle:
        handle.write(raw)
    return sha(raw)


def consume(hashes: Mapping[str, str], progress: dict[str, Any]) -> None:
    """Authenticate, import the two new files, and invoke the fixture once."""

    _authenticate(hashes)
    module = _import_new(MODULE, NAME + "_module")
    progress["imported_modules"] = [MODULE]
    fixture = _import_new(FIXTURE, NAME + "_fixture")
    progress["imported_modules"] = [MODULE, FIXTURE]

    runner = getattr(fixture, "run_fixture", None)
    require(callable(runner), "fixture must export run_fixture(module)")
    progress["fixture_call_count"] = 1
    # This is the sole fixture invocation.  The runner intentionally does not
    # inspect a passed flag as a physical or source-owner admission.
    value = runner(module)
    progress["fixture_result"] = _json_value(value)
    progress["fixture_returned"] = True

    # Re-authenticate all frozen bytes after the call.  Both source-definition
    # files are still only byte/hash inputs; neither was imported or evaluated.
    _authenticate(hashes)


def main() -> int:
    if not __debug__:
        raise RuntimeError("optimized mode disables one-entry guards")

    try:
        hashes, source_sha256, brief_sha256 = preflight()
    except ExistingArtifactError:
        print("status=rejected_existing_entry_or_result")
        return 1
    except BaseException:
        print("status=preflight_failed")
        return 1

    try:
        entry_sha256 = _write_exclusive(
            ENTRY,
            {
                "component": NAME,
                "source_sha256": source_sha256,
                "brief_sha256": brief_sha256,
                "input_hashes": hashes,
                "source_definition_paths": SOURCE_DEFINITIONS,
                "source_definition_imported": False,
                "fixture_api": FIXTURE_API,
                "entered": True,
                "exclusive": True,
                **SYNTHETIC_ONLY,
            },
        )
    except FileExistsError:
        print("status=rejected_entry_race")
        return 1
    except BaseException:
        print("status=entry_write_failed")
        return 1

    out: dict[str, Any] = {
        "component": NAME,
        "source_sha256": source_sha256,
        "brief_sha256": brief_sha256,
        "entry_sha256": entry_sha256,
        "input_hashes": hashes,
        "source_definition_paths": SOURCE_DEFINITIONS,
        "source_definition_imported": False,
        "status": "failed",
        "fixture_api": FIXTURE_API,
        "fixture_call_count": 0,
        "imported_modules": [],
        "fixture_returned": False,
        **SYNTHETIC_ONLY,
    }
    captured_stdout = io.StringIO()
    captured_stderr = io.StringIO()
    try:
        with contextlib.redirect_stdout(captured_stdout), contextlib.redirect_stderr(captured_stderr):
            consume(hashes, out)
    except BaseException as exc:
        out["error"] = _error(exc)
    else:
        out["status"] = "returned"
    finally:
        out["captured_module_stdout"] = captured_stdout.getvalue()
        out["captured_module_stderr"] = captured_stderr.getvalue()

    try:
        _write_exclusive(RESULT, out)
    except FileExistsError:
        print("status=rejected_result_race")
        return 1
    except BaseException:
        print("status=result_write_failed")
        return 1
    out_status = out["status"]
    # The launcher emits status only.  Fixture output and tracebacks live in
    # result.json so a leaf cannot pollute the one-line launcher protocol.
    print("status=" + str(out_status))
    return 0 if out_status == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
