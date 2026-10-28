"""Enter and invoke the frozen M-history source binding exactly once."""

from pathlib import Path
import hashlib
import importlib.util
import json
import sys
import traceback


R = Path(__file__).resolve().parents[1]
B = R.parent.parent / "RESUME_BRIEF_2026_10_16.json"
L = R / "tmp/compact_M_history_binding_v1"
SOURCE, CAPSULE, ENTRY, RESULT = (
    L / name for name in ("source.json", "capsule.json", "entered.json", "result.json")
)
CODE = R / "src/openai_ns_reconstruction/compact_M_history_binding.py"
COMPONENT = "compact_M_history_binding_v1"
PENDING_PHASE = "compact_M_history_binding_pending"


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def save_new(path: Path, value: dict) -> str:
    raw = (
        json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)
        + "\n"
    ).encode("utf-8")
    with path.open("xb") as stream:
        stream.write(raw)
    return hashlib.sha256(raw).hexdigest()


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def _path_from_manifest(raw_path: str) -> Path:
    path = Path(raw_path)
    return path if path.is_absolute() else R / path


def _validate_checkpoint_hashes(source: dict, capsule: dict, brief: dict) -> None:
    """Validate every declared source/capsule input before writing ``entered``."""
    _require(source.get("component") == COMPONENT, "source component mismatch")
    _require(capsule.get("component", COMPONENT) == COMPONENT,
             "capsule component mismatch")
    _require(source.get("capsule_sha256") == sha(CAPSULE),
             "capsule checkpoint SHA mismatch")
    checkpoint = brief.get("source_checkpoint")
    _require(isinstance(checkpoint, dict),
             "M history source checkpoint is absent from pending brief")
    _require(checkpoint.get("sha256") == sha(SOURCE),
             "source checkpoint SHA mismatch")

    declared = source.get("files")
    _require(isinstance(declared, list) and declared,
             "source file hash manifest is absent")
    for item in declared:
        _require(isinstance(item, dict), "invalid source file hash entry")
        path = _path_from_manifest(item["path"])
        _require(path.is_file(), f"missing frozen file: {path}")
        _require(sha(path) == item["sha256"],
                 f"frozen file changed: {item['path']}")

    input_hashes = capsule.get("input_sha256", {})
    _require(isinstance(input_hashes, dict), "capsule input hash map is absent")
    for raw_path, expected in input_hashes.items():
        path = _path_from_manifest(raw_path)
        _require(path.is_file(), f"missing capsule input: {path}")
        _require(sha(path) == expected,
                 f"capsule input changed: {raw_path}")


def main() -> int:
    if ENTRY.exists() or RESULT.exists():
        raise RuntimeError("already entered; never reenter")

    source = json.loads(SOURCE.read_bytes())
    capsule = json.loads(CAPSULE.read_bytes())
    brief = json.loads(B.read_bytes())
    _require(brief.get("component") == COMPONENT, "pending brief component mismatch")
    _require(brief.get("phase") == PENDING_PHASE, "pending brief phase mismatch")
    _require(brief.get("goal_active") is True, "pending brief is not active")
    _require(brief.get("full_goal_complete") is False,
             "pending brief already marks the full goal complete")
    _require(brief.get("do_not_reenter") is True,
             "pending brief lacks the no-reentry guard")
    _validate_checkpoint_hashes(source, capsule, brief)

    entry = save_new(ENTRY, {
        "component": COMPONENT,
        "source_sha256": sha(SOURCE),
        "capsule_sha256": sha(CAPSULE),
        "entered": True,
        "consumer_calls_started": 0,
        "ancestor_math_calls": 0,
        "runtime_model_metadata": None,
    })

    try:
        spec = importlib.util.spec_from_file_location("_compact_M_history_binding", CODE)
        if spec is None or spec.loader is None:
            raise RuntimeError("kernel load failed")
        module = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = module
        spec.loader.exec_module(module)
        result = module.bind_M_history_binding(R)
        outcome = {
            "component": COMPONENT,
            "status": "returned",
            "source_sha256": sha(SOURCE),
            "entry_sha256": entry,
            "capsule_sha256": sha(CAPSULE),
            "consumer_calls": 1,
            "ancestor_math_calls": 0,
            "runtime_model_metadata": None,
            "result": result,
        }
    except Exception as exc:
        outcome = {
            "component": COMPONENT,
            "status": "failed",
            "source_sha256": sha(SOURCE),
            "entry_sha256": entry,
            "capsule_sha256": sha(CAPSULE),
            "consumer_calls": 1,
            "ancestor_math_calls": 0,
            "runtime_model_metadata": None,
            "exception": f"{type(exc).__name__}: {exc}",
            "traceback": traceback.format_exc(),
        }

    digest = save_new(RESULT, outcome)
    print(json.dumps({
        "component": COMPONENT,
        "status": outcome["status"],
        "result_sha256": digest,
        "consumer_calls": 1,
        "ancestor_math_calls": 0,
        "runtime_model_metadata": None,
    }))
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
