"""One-entry source-only W0 ordinary-Z2 four-cell consumer.

The runner authenticates every raw input before decoding the three gzip
receipts.  It passes only the selected ordinary-Z2 core and forcing rows to
the new leaf; the older W0 Z0/Z1 receipt is checked as a separate provenance
record and is never rebuilt or added to the new result.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from pathlib import Path

from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
COMPONENT = "w0_z2_four_cell_v2"
DIR = ROOT / "tmp" / COMPONENT
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entered.json", "result.json"))

GREEN = "experiments/root_st073/n1_actual_collar_green_v3_receipt.json.gz"
FORCING = "experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz"
W0_Z0Z1 = "experiments/root_st073/n1_w0_collar_sectors_v1_receipt.json.gz"
SELECTED = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json"
LEAF = "src/openai_ns_reconstruction/w0_z2_four_cell_v2.py"
FIXTURE = "src/openai_ns_reconstruction/w0_z2_four_cell_fixture_v2.py"
RUNNER = "scripts/run_w0_z2_four_cell_v2_once.py"

FAILURE = "experiments/root_st073/w0_z2_four_cell_v1_failure.json.gz"
REQUIRED = (LEAF, FIXTURE, RUNNER, GREEN, FORCING, W0_Z0Z1, SELECTED, FAILURE)
FIXED = {
    FAILURE: "35004fd38c4106bbe8d4c48dd45fbd7e56982458740b29e047565ce066197468",
    GREEN: "fad35594ccf5ef81392f94f97ae72b3ffb94fcaa554438fd5702ef2c4343c4b5",
    FORCING: "96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2",
    W0_Z0Z1: "4229c29a2a21207f8b8da3ba53d70bbf6f50e09eb67110bd3b7cb30fa6577c06",
    SELECTED: "0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e",
}

FIELDS = ("F1", "Uz1", "K1")
ROWS = ("5", "6")
ATOMS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def require(condition: object, message: str) -> None:
    if not condition:
        raise ValueError(message)


def save(path: Path, value: object) -> str:
    raw = (json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
    with path.open("xb") as handle:
        handle.write(raw)
    return sha(raw)


def interval(c: MPIntervalContext, record: object, label: str = "interval") -> object:
    require(isinstance(record, dict), f"{label} must be an exact interval record")
    lower = record.get("lower_exact_mpf_tuple")
    upper = record.get("upper_exact_mpf_tuple")
    require(isinstance(lower, (list, tuple)) and isinstance(upper, (list, tuple)), f"{label} endpoints missing")
    require(len(lower) == len(upper) == 4, f"{label} endpoints must be four-tuples")
    value = c.make_mpf((tuple(lower), tuple(upper)))
    require(value.a <= value.b, f"{label} interval is reversed")
    return value


def encode(value: object) -> object:
    """Encode every interval by its exact directed endpoint tuples."""

    if hasattr(value, "_mpi_"):
        lower, upper = value._mpi_
        return {
            "lower_exact_mpf_tuple": list(lower),
            "upper_exact_mpf_tuple": list(upper),
            "directed": True,
        }
    if isinstance(value, dict):
        return {str(key): encode(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError(f"unsupported output value {type(value).__name__}")


def module(path: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    require(spec is not None and spec.loader is not None, f"cannot load {path}")
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def preflight() -> tuple[dict[str, str], str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    require(SOURCE.exists(), "source manifest missing")
    brief_raw = BRIEF.read_bytes()
    brief = json.loads(brief_raw)
    source_raw = SOURCE.read_bytes()
    checkpoint = brief.get("source_checkpoint")
    require(brief.get("component") == COMPONENT, "brief component guard")
    require(brief.get("phase") == COMPONENT + "_pending", "brief phase guard")
    require(brief.get("goal_active") is True and brief.get("full_goal_complete") is False, "goal guard")
    require(brief.get("one_entry_only") is True and brief.get("do_not_reenter") is True, "one-entry guard")
    require(isinstance(checkpoint, dict), "source checkpoint missing")
    require(Path(checkpoint["path"]).resolve() == SOURCE.resolve(), "source checkpoint path mismatch")
    require(checkpoint["sha256"] == sha(source_raw), "source checkpoint hash mismatch")

    manifest = json.loads(source_raw)
    require(manifest.get("component") == COMPONENT and manifest.get("one_entry_only") is True, "manifest guard")
    hashes: dict[str, str] = {}
    for row in manifest.get("files", ()):
        path = row.get("path")
        require(path in REQUIRED and path not in hashes, f"unexpected or duplicate manifest path: {path}")
        digest = sha((ROOT / path).read_bytes())
        require(digest == row.get("sha256"), f"manifest hash mismatch: {path}")
        if path in FIXED:
            require(digest == FIXED[path], f"pinned input changed: {path}")
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), "manifest input set mismatch")
    return hashes, sha(source_raw), sha(brief_raw)


def _bind_terms(
    c: MPIntervalContext,
    terms: object,
    *,
    receipt: str,
    path_prefix: list[object],
    stable_prefix: str,
    source_index: dict[str, dict[str, object]],
) -> list[dict[str, object]]:
    require(isinstance(terms, (list, tuple)), f"term list missing at {path_prefix}")
    bound: list[dict[str, object]] = []
    for index, raw in enumerate(terms):
        require(isinstance(raw, dict), f"term {path_prefix + [index]} is not a mapping")
        require("coefficient" in raw and "powers" in raw and "origin" in raw, f"incomplete term {path_prefix + [index]}")
        powers = raw["powers"]
        require(isinstance(powers, (list, tuple)) and len(powers) == len(ATOMS), f"formal powers missing at {path_prefix + [index]}")
        require(all(isinstance(power, int) and not isinstance(power, bool) for power in powers), f"noninteger formal power at {path_prefix + [index]}")
        origin = raw["origin"]
        require(isinstance(origin, dict), f"origin metadata missing at {path_prefix + [index]}")
        stable_id = f"{stable_prefix}:{index}"
        source_ids = origin.get("source_ids", [])
        source_paths = origin.get("source_paths", [])
        source_index[stable_id] = {
            "receipt": receipt,
            "path": path_prefix + [index],
            "original_source_ids": copy.deepcopy(source_ids),
            "original_source_paths": copy.deepcopy(source_paths),
            "origin": copy.deepcopy(origin),
            "powers": list(powers),
        }
        term = copy.deepcopy(raw)
        term["coefficient"] = interval(c, raw["coefficient"], f"coefficient {path_prefix + [index]}")
        # These stable identifiers are additive metadata.  The original
        # origin/source IDs and error classes above are not rewritten.
        term["input_id"] = stable_id
        term["stable_id"] = stable_id
        term["receipt_source_path"] = [receipt, *path_prefix, index]
        bound.append(term)
    return bound


def _collect_stable_ids(value: object, known: set[str], seen: set[str]) -> None:
    if isinstance(value, dict):
        for key in ("input_id", "stable_id"):
            item = value.get(key)
            if isinstance(item, str) and item in known:
                seen.add(item)
        origin = value.get("origin")
        if isinstance(origin, dict):
            for key in ("input_id", "stable_id"):
                item = origin.get(key)
                if isinstance(item, str) and item in known:
                    seen.add(item)
        for item in value.values():
            _collect_stable_ids(item, known, seen)
    elif isinstance(value, (list, tuple)):
        for item in value:
            _collect_stable_ids(item, known, seen)


def _validate_q_cells(c: MPIntervalContext, cells: object, label: str) -> None:
    require(isinstance(cells, (list, tuple)) and len(cells) == 4, f"{label} must contain four q cells")
    for index, cell in enumerate(cells):
        require(isinstance(cell, dict), f"{label}[{index}] must be a mapping")
        for key in ("q_lower", "q_upper", "s_lower", "s_upper"):
            require(key in cell, f"{label}[{index}] missing {key}")
        q_lower = c.mpf(index) / 4
        q_upper = c.mpf(index + 1) / 4
        def same_endpoint(value: object, expected: object) -> bool:
            if hasattr(value, "_mpi_"):
                return value._mpi_ == expected._mpi_
            if isinstance(value, dict):
                return interval(c, value, f"{label}[{index}] endpoint")._mpi_ == expected._mpi_
            return False

        require(same_endpoint(cell["q_lower"], q_lower), f"{label}[{index}] q lower mismatch")
        require(same_endpoint(cell["q_upper"], q_upper), f"{label}[{index}] q upper mismatch")


def consume(hashes: dict[str, str], progress: dict[str, object]) -> dict[str, object]:
    # Authenticate every file first.  Only after this loop may receipt bytes
    # be decoded; the gzip payloads themselves are never used for hashing.
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], f"input changed after entry: {path}")
        raw[path] = data

    c = MPIntervalContext()
    c.prec = 1800
    leaf = module(LEAF, "w0_z2_four_cell_v2_leaf")
    fixture = module(FIXTURE, "w0_z2_four_cell_v2_fixture")
    fixture_result = fixture.run_fixture(leaf.build, c)
    progress["fixture"] = encode(fixture_result)
    require(isinstance(fixture_result, dict) and fixture_result.get("passed") is True, "synthetic fixture failed")
    green = json.loads(gzip.decompress(raw[GREEN]))
    forcing = json.loads(gzip.decompress(raw[FORCING]))
    old_w0 = json.loads(gzip.decompress(raw[W0_Z0Z1]))
    selected = json.loads(raw[SELECTED].decode("utf-8"))
    for name, receipt in ((GREEN, green), (FORCING, forcing), (W0_Z0Z1, old_w0)):
        require(receipt.get("result", {}).get("status") == "returned", f"{name} is not returned")
    green_output = green["result"]["output"]
    forcing_output = forcing["result"]["output"]
    old_output = old_w0["result"]["output"]
    green_inputs = {row["path"]: row["sha256"] for row in green["source"]["files"]}
    old_inputs = {row["path"]: row["sha256"] for row in old_w0["source"]["files"]}
    require(green_inputs.get(FORCING) == hashes[FORCING], "green forcing chain")
    require(old_inputs.get(GREEN) == hashes[GREEN] and old_inputs.get(FORCING) == hashes[FORCING], "old W0 chain")
    metadata = forcing_output.get("calculation_metadata", {})
    require(metadata.get("F0_relative_derivative_applied_once") is True, "forcing F0 derivative")
    require(metadata.get("source_derivatives_are_ordinary") is True, "ordinary source derivatives")
    family = green_output["source_family"]
    source_logs = green_output["source_logs"]
    require(forcing_output["source_family"] == family and old_output["source_family"] == family, "source family mismatch")
    require(forcing_output["source_logs"] == source_logs and old_output["source_logs"] == source_logs, "source logs mismatch")
    require(selected.get("source_family") == family, "selected source mismatch")
    require(forcing_output.get("domain") == green_output.get("domain") == old_output.get("domain"), "source domain mismatch")
    measure = selected["actual_cutoff_complete_support_measure"]
    require(measure.get("passed") is True and measure.get("source_family") == family, "support measure")
    require(measure.get("measure_is_Z_independent") is True and "0<hb<1/2" in measure.get("inequality_proof", ""), "geometry premise")
    sc = interval(c, green_output["sc"], "green sc")
    exact_sc = tuple(measure["source_sc_exact_mpf_tuple"])
    require(sc._mpi_ == (exact_sc, exact_sc), "selected exact sc mismatch")
    if "sc" in old_output:
        require(interval(c, old_output["sc"], "old W0 sc")._mpi_ == sc._mpi_, "old W0 sc mismatch")
    old_result = old_output["result"]
    require(old_output.get("saved_core_transport_receipt") == GREEN and old_output.get("source_forcing_receipt") == FORCING, "old W0 pins")
    require(tuple(old_result.get("ordinary_Z_orders", ())) == (0, 1), "old W0 orders")
    _validate_q_cells(c, old_result.get("q_cells"), "old W0 q_cells")
    require(old_output.get("allword_tail_added") is False and old_output.get("m15_computed") is False and old_output.get("appended_to_cumulative_m15") is False, "old W0 scope")
    core_saved = green_output["result"]["core_W0_transport"]["2"]
    forcing_saved = forcing_output["g"]["2"]
    require(set(FIELDS).issubset(core_saved) and set(ROWS).issubset(forcing_saved), "source rows missing")
    source_index: dict[str, dict[str, object]] = {}
    core = {field: _bind_terms(c, core_saved[field], receipt=GREEN,
        path_prefix=["result", "output", "result", "core_W0_transport", "2", field],
        stable_prefix=f"green:W0_Z2:{field}", source_index=source_index) for field in FIELDS}
    g = {row: _bind_terms(c, forcing_saved[row], receipt=FORCING,
        path_prefix=["result", "output", "g", "2", row],
        stable_prefix=f"forcing:g_Z2:{row}", source_index=source_index) for row in ROWS}
    # Keep the complete synthetic return in the outcome before checking its
    # status, so a failed fixture remains reviewable in the saved result.
    outcome: dict[str, object] = {
        "source_family": family,
        "source_logs": source_logs,
        "domain": green_output["domain"],
        "sc": green_output["sc"],
        "fixture": fixture_result,
        "input_index": source_index,
        "green_receipt": GREEN,
        "forcing_receipt": FORCING,
        "prior_W0_Z0_Z1_receipt": W0_Z0Z1,
        "selected_words": SELECTED,
        "actual_order": 2,
        "old_orders_remain_separate": [0, 1],
        "old_errors_and_tails_remain_separate": True,
    }
    require(isinstance(fixture_result, dict) and fixture_result.get("passed") is True, "synthetic fixture failed")

    # This is the only actual invocation of the new leaf.  The fixture above
    # receives the callable for synthetic checks, but its calls are not actual
    # receipt evaluation.
    result = leaf.build(c, sc=sc, core=core, g=g)
    require(isinstance(result, dict), "new leaf result must be a mapping")
    require(result.get("ordinary_Z_order") == 2 and tuple(result.get("ordinary_Z_orders", ())) == (2,), "new result is not ordinary Z2")
    require(result.get("source_z2_is_ordinary_derivative") is True, "ordinary Z2 marker missing")
    require(result.get("ordinary_Z_factorial_divided") is False, "ordinary Z2 factorial semantics changed")
    require(tuple(result.get("source_power_basis", ())) == ATOMS, "formal atom basis changed")
    require(set(result.get("core_transport", {})) == set(FIELDS), "new core transport fields mismatch")
    _validate_q_cells(c, result.get("q_cells"), "new q_cells")
    require(result.get("allword_tail_added") is False and result.get("m15_computed") is False, "new scope added tail/m15")
    require(result.get("old_W1_consumed") is False and result.get("same_source_cancellation_claimed") is False, "new scope changed")

    seen: set[str] = set()
    _collect_stable_ids(result, set(source_index), seen)
    require(seen == set(source_index), "new result lost a stable input ID")
    outcome.update(
        {
            "result": encode(result),
            "distinct_input_ids": len(source_index),
            "distinct_output_input_ids": len(seen),
            "source_terms_not_merged": True,
            "formal_powers_preserved": True,
            "exact_interval_endpoints_encoded": True,
            "allword_tail_added": False,
            "m15_computed": False,
            "cartesian_owner_proved": False,
            "ledger_appended": False,
        }
    )
    return outcome


def main() -> int:
    hashes, source_sha, brief_sha = preflight()
    entry_sha = save(
        ENTRY,
        {
            "component": COMPONENT,
            "source_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "entered": True,
        },
    )
    output: dict[str, object] = {
        "component": COMPONENT,
        "source_sha256": source_sha,
        "entry_sha256": entry_sha,
        "status": "failed",
    }
    try:
        output.update(status="returned", output=consume(hashes, output), error=None)
    except BaseException as exc:
        output["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    result_sha = save(RESULT, output)
    print(json.dumps({"component": COMPONENT, "status": output["status"], "result_sha256": result_sha}))
    return 0 if output["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
