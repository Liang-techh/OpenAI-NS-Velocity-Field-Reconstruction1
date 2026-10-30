"""Freeze and enter the new non-common Gamma Z1/admission component once.

No historical producer or evaluator is imported. Every declared input is
authenticated as raw bytes before entry and before any archive decompression.
An exclusive entry receipt precedes all imports and source admission, and a
success or failure result is persisted before it is printed.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "gamma_noncommon_admission_v1"
PHASE = "gamma_noncommon_admission_pending"
DIRECTORY = ROOT / "tmp" / COMPONENT
SOURCE = DIRECTORY / "source.json"
ENTRY = DIRECTORY / "entered.json"
RESULT = DIRECTORY / "result.json"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def write_once(path: Path, obj: dict) -> str:
    data = json.dumps(obj, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode()
    with path.open("xb") as f:
        f.write(data)
    return hashlib.sha256(data).hexdigest()


def load_new_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, "cannot load new leaf")
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def preflight() -> tuple[dict, dict[str, str], str]:
    require(not ENTRY.exists() and not RESULT.exists(), "one-entry component already entered")
    brief = json.loads(BRIEF.read_bytes())
    source = json.loads(SOURCE.read_bytes())
    require(brief.get("component") == COMPONENT and brief.get("phase") == PHASE,
            "brief not pending for this component")
    require(brief.get("goal_active") is True and brief.get("full_goal_complete") is False
            and brief.get("do_not_reenter") is True, "brief goal flags mismatch")
    checkpoint = brief.get("source_checkpoint", {})
    require(Path(checkpoint.get("path", "")).resolve() == SOURCE.resolve()
            and checkpoint.get("sha256") == sha(SOURCE), "checkpoint mismatch")
    require(source.get("component") == COMPONENT and source.get("one_entry_only") is True
            and source.get("no_source_graph_traversal") is True
            and source.get("all_raw_inputs_authenticated_before_decompression") is True,
            "checkpoint safeguards missing")
    files = source.get("files")
    require(isinstance(files, list) and files, "empty input manifest")
    hashes: dict[str, str] = {}
    for item in files:
        relative = item["path"]
        candidate = Path(relative)
        require(not candidate.is_absolute() and ".." not in candidate.parts,
                "manifest path must stay inside repository")
        path = (ROOT / candidate).resolve()
        require(path.is_relative_to(ROOT.resolve()) and path.is_file(), "invalid input path")
        require(relative not in hashes and item["sha256"] == sha(path), "raw SHA mismatch")
        hashes[relative] = item["sha256"]
    required = {
        "src/openai_ns_reconstruction/gamma_noncommon_z1.py",
        "src/openai_ns_reconstruction/gamma_source_admission.py",
        "scripts/run_gamma_noncommon_admission_once.py",
        "tmp/gamma_weighted_energy_delta_v1/result.json",
    }
    for stem in (
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions",
        "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_complete_future_functions",
    ):
        required.update({stem + ".py", stem + "_check.json", stem + ".json.gz"})
    require(required.issubset(hashes), "required frozen inputs missing")
    return source, hashes, sha(SOURCE)


def main() -> int:
    source, hashes, checkpoint_sha = preflight()
    entry_sha = write_once(ENTRY, {"component": COMPONENT, "phase": PHASE,
        "source_checkpoint_sha256": checkpoint_sha, "input_hashes": hashes,
        "consumer_calls_started": 0, "ancestor_math_calls": 0,
        "source_values_computed": False, "runtime_model_metadata": None})
    calls = []
    try:
        algebra = load_new_module(ROOT / "src/openai_ns_reconstruction/gamma_noncommon_z1.py",
                                  "gamma_noncommon_z1_new_leaf")
        admission = load_new_module(ROOT / "src/openai_ns_reconstruction/gamma_source_admission.py",
                                    "gamma_source_admission_new_leaf")
        calls.append("noncommon_synthetic_fixture")
        algebra_fixture = algebra.verify_synthetic_fixture()
        require(algebra_fixture.get("checks_passed") is True,
                "noncommon synthetic fixture failed")
        calls.append("admission_synthetic_negative_fixture")
        negative_fixture = admission.synthetic_negative_fixture()
        require(negative_fixture.get("checks_passed") is True,
                "admission negative fixture failed")
        calls.append("canonical_source_admission")
        accepted = admission.admit_canonical_gamma_sources(ROOT)
        require(accepted.get("admitted") is True and
                accepted.get("P_F_owner_correspondence_proved") is False and
                accepted.get("full_future_P_F_proved") is False,
                "canonical source admission flags invalid")
        output = {"component": COMPONENT, "phase": PHASE, "status": "returned",
                  "source_checkpoint_sha256": checkpoint_sha, "entry_sha256": entry_sha,
                  "input_hashes": hashes, "consumer_calls": calls,
                  "ancestor_math_calls": 0, "source_values_computed": False,
                  "certified": False, "runtime_model_metadata": None,
                  "result": {"algebra_synthetic_fixture": algebra_fixture,
                             "admission_negative_fixture": negative_fixture,
                             "canonical_source_admission": accepted}, "error": None}
    except Exception as exc:
        output = {"component": COMPONENT, "phase": PHASE, "status": "failed",
                  "source_checkpoint_sha256": checkpoint_sha, "entry_sha256": entry_sha,
                  "input_hashes": hashes, "consumer_calls": calls,
                  "ancestor_math_calls": 0, "source_values_computed": False,
                  "certified": False, "runtime_model_metadata": None,
                  "error": {"type": type(exc).__name__, "message": str(exc)},
                  "traceback": traceback.format_exc()}
    result_sha = write_once(RESULT, output)
    print(json.dumps({"component": COMPONENT, "status": output["status"],
                      "result_sha256": result_sha, "consumer_calls": calls}, sort_keys=True))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
