"""One-entry sealed runner for the source-bound local second-potential m15 leaf.

The runner authenticates the frozen manifest and every raw input before any
receipt is decompressed.  It only imports the requested leaf and its
independent fixture; all ancestor results are consumed as sealed receipts.
"""

from __future__ import annotations

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
COMPONENT = "n1_local_second_potential_m15_v2"
DIR = ROOT / "tmp" / COMPONENT
SOURCE = DIR / "source.json"
ENTRY = DIR / "entry.json"
RESULT = DIR / "result.json"

FORCING = "experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz"
KCHI = "experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz"
BIND = "experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz"
PRESS = "experiments/root_st073/n1_first_collar_g4_pressure_v1_receipt.json.gz"
SWIRL = "experiments/root_st073/n1_first_word_swirl_m15_v1_receipt.json.gz"
CORRECTION = "experiments/root_st073/n1_first_word_swirl_m15_v1_provenance_correction.json"
SELECT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json"
PREVIOUS = "experiments/root_st073/n1_collar_qp_k0_m15_cells_v1_receipt.json.gz"
OLD_TAIL = "experiments/root_st073/n1_collar_local_g5_j2_tail_v1_receipt.json.gz"
PAIR = "experiments/root_st073/n1_local_p_dg_correlated_pair_v1_receipt.json.gz"
COMP = "experiments/root_st073/n1_core_kchi_compensation_v2_receipt.json.gz"

LEAF = "src/openai_ns_reconstruction/n1_local_second_potential_m15_v2.py"
FIXTURE = "src/openai_ns_reconstruction/n1_local_second_potential_m15_fixture_v2.py"
RUNNER = "scripts/run_n1_local_second_potential_m15_v2_once.py"

REQUIRED = (
    FORCING,
    KCHI,
    BIND,
    PRESS,
    SWIRL,
    CORRECTION,
    SELECT,
    PREVIOUS,
    OLD_TAIL,
    PAIR,
    COMP,
    LEAF,
    FIXTURE,
    RUNNER,
)

FIXED = {
    FORCING: "96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2",
    KCHI: "b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0",
    BIND: "5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8",
    PRESS: "c229cf188d306413f6edd0641a729321ad4567c2587ebf4701ada47720853660",
    SWIRL: "a626e71fa2b0ddf40b79ce992016da2af8ac9955db3f17337029f3c86c2513b4",
    CORRECTION: "f2133bf5824b2f9636b430a33a2ba1d19fbf391a7f90a81de0c8269d1fe638c4",
    SELECT: "0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e",
    PREVIOUS: "48409ed314e95fe0f3117f6c27ed74b5aa8b56822c810ba5050d2b000b08e5c2",
    OLD_TAIL: "67894bac58f4fbbcd5437bf970dadafeb918fe2d9638425f8e5165d21af4530b",
    PAIR: "7fabbd9cd30ccdf6e64288ebdd56bc7e744f3023a2ff1966f3e9e96de2651b54",
    COMP: "c862f0560e8cb4b053abbbfcc982e4b4c1b602a806628e2d709399878e1c3aaf",
}


def sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def require(value: object, message: str) -> None:
    if not value:
        raise ValueError(message)


def save(path: Path, value: object) -> str:
    data = (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()
    with path.open("xb") as handle:
        handle.write(data)
    return sha(data)


def preflight() -> tuple[dict[str, str], str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "already entered")
    brief_bytes = BRIEF.read_bytes()
    brief = json.loads(brief_bytes)
    source_bytes = SOURCE.read_bytes()
    checkpoint = brief["source_checkpoint"]
    require(
        brief["component"] == COMPONENT
        and brief["phase"] == COMPONENT + "_pending"
        and brief["goal_active"] is True
        and brief["full_goal_complete"] is False
        and brief["one_entry_only"] is True
        and brief["do_not_reenter"] is True,
        "brief guard",
    )
    require(
        Path(checkpoint["path"]).resolve() == SOURCE.resolve()
        and checkpoint["sha256"] == sha(source_bytes),
        "source checkpoint",
    )
    manifest = json.loads(source_bytes)
    require(manifest["component"] == COMPONENT and manifest["one_entry_only"] is True, "manifest")

    hashes: dict[str, str] = {}
    for row in manifest["files"]:
        path = row["path"]
        require(path in REQUIRED and path not in hashes, "unexpected input")
        digest = sha((ROOT / path).read_bytes())
        require(digest == row["sha256"] and (path not in FIXED or digest == FIXED[path]), "changed " + path)
        hashes[path] = digest
    require(set(hashes) == set(REQUIRED), "missing input")
    return hashes, sha(source_bytes), sha(brief_bytes)


def iv(context: MPIntervalContext, record: dict) -> object:
    lower = tuple(record["lower_exact_mpf_tuple"])
    upper = tuple(record["upper_exact_mpf_tuple"])
    require(len(lower) == len(upper) == 4, "interval endpoints")
    value = context.make_mpf((lower, upper))
    require(value.a <= value.b, "reversed interval")
    return value


def enc(value: object) -> object:
    if hasattr(value, "_mpi_"):
        return {
            "lower_exact_mpf_tuple": list(value._mpi_[0]),
            "upper_exact_mpf_tuple": list(value._mpi_[1]),
            "directed": True,
        }
    if isinstance(value, dict):
        return {str(key): enc(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [enc(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError(type(value).__name__)


def module(path: str, name: str):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None:
        raise ImportError(path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def _source_pins(receipt: dict) -> dict[str, str]:
    return {item["path"]: item["sha256"] for item in receipt["source"]["files"]}


def consume(hashes: dict[str, str]) -> dict:
    # Read and authenticate every raw byte before decompression or receipt inspection.
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "changed after entry " + path)
        raw[path] = data

    context = MPIntervalContext()
    context.prec = 1800
    leaf = module(LEAF, "n1_local_second_potential_m15_v2")
    fixture = module(FIXTURE, "n1_local_second_potential_m15_fixture_v2")
    proof = fixture.run_fixture(leaf.build, context)
    require(proof["passed"], "fixture failed " + json.dumps(proof, sort_keys=True, separators=(",", ":")))

    receipts = {
        path: json.loads(gzip.decompress(raw[path]))
        for path in (FORCING, KCHI, BIND, PRESS, SWIRL, PREVIOUS, OLD_TAIL, PAIR, COMP)
    }
    require(all(receipt["result"]["status"] == "returned" for receipt in receipts.values()), "receipt status")

    forcing = receipts[FORCING]["result"]["output"]
    kchi = receipts[KCHI]["result"]["output"]
    binding = receipts[BIND]["result"]["output"]
    pressure = receipts[PRESS]["result"]["output"]
    swirl = receipts[SWIRL]["result"]["output"]
    previous = receipts[PREVIOUS]["result"]["output"]
    old_receipt = receipts[OLD_TAIL]["result"]["output"]
    pair = receipts[PAIR]["result"]["output"]
    compensation = receipts[COMP]["result"]["output"]
    selected = json.loads(raw[SELECT])
    correction = json.loads(raw[CORRECTION])

    family = forcing["source_family"]
    atoms = forcing["source_logs"]
    for name, value in (
        ("kchi", kchi),
        ("binding", binding),
        ("pressure", pressure),
        ("swirl", swirl),
        ("previous", previous),
        ("old_tail", old_receipt),
        ("pair", pair),
        ("compensation", compensation),
    ):
        require(value["source_family"] == family and value["source_logs"] == atoms, "common source family/logs: " + name)
    require(selected["source_family"] == family, "selected source family")

    # Verify the direct receipt-to-receipt source chain without importing any ancestor.
    old_pins = _source_pins(receipts[OLD_TAIL])
    previous_pins = _source_pins(receipts[PREVIOUS])
    binding_pins = _source_pins(receipts[BIND])
    compensation_pins = _source_pins(receipts[COMP])
    pair_pins = _source_pins(receipts[PAIR])
    require(
        old_pins[FORCING] == hashes[FORCING]
        and all(old_pins[path] == hashes[path] for path in (KCHI, BIND, PRESS, SELECT))
        and binding_pins[KCHI] == hashes[KCHI]
        and all(compensation_pins[path] == hashes[path] for path in (FORCING, KCHI, PRESS, SELECT))
        and all(pair_pins[path] == hashes[path] for path in (FORCING, KCHI, BIND, COMP, PRESS, SELECT))
        and all(previous_pins[path] == hashes[path] for path in (SWIRL, CORRECTION, SELECT)),
        "authenticated source chain",
    )
    require(binding['source_checks']['one_query_jet_shared_in_integration'] is True
            and binding['source_checks']['native_axis_jet_uses_Z_only'] is True
            and binding['phase_shear_cutoff_identified_with_axis_potential'] is False,
            'Z-only shared axis-potential owner')
    require(compensation['normalization_identity_scope'] == 'S-conjugated Achi extraction implies physical +2chi E51; remaining gauge/source identity remains conditional', 'conditional physical P convention')
    require(pair['source_paths_aligned'] is True and len(pair['source_matches']) == 32
            and old_receipt['local_j2_tail_only'] is True and old_receipt['m15_computed'] is False,
            'paired source scope and archived tail')
    require(
        correction["receipt_sha256"] == hashes[SWIRL]
        and correction["numerical_result_changed"] is False
        and correction["input_id"] == "F0_connection",
        "leading amplitude correction",
    )

    contract = selected["actual_cutoff_complete_support_measure"]
    sc = iv(context, pressure["sc"])
    exact_sc = tuple(contract["source_sc_exact_mpf_tuple"])
    require(
        sc._mpi_ == (exact_sc, exact_sc)
        and sc.a > 0
        and pressure["sc"] == forcing.get("sc", pressure["sc"])
        and swirl['sc'] == pressure['sc']
        and contract["passed"]
        and contract["measure_is_Z_independent"]
        and "0<hb<1/2" in contract["inequality_proof"],
        "physical geometry",
    )

    source_index: dict[str, dict] = {}
    g5: dict[str, list[dict]] = {"Z0": [], "Z1": []}
    for order in (0, 1):
        label = "Z" + str(order)
        terms: list[dict] = []
        for index, term in enumerate(forcing["g"][str(order)]["5"]):
            sid = f"collar_g5:Z{order}:{index}"
            origin = term.get("origin", {})
            parents = list(origin.get("source_ids", term.get("source_ids", [])))
            require(parents, "forcing parent source ids")
            pointer = {
                "receipt": FORCING,
                "path": ["result", "output", "g", str(order), "5", index],
                "parent_source_ids": parents,
            }
            require(old_receipt["source_index"].get(sid) == pointer, "old tail source path " + sid)
            matches = [m for m in pair['source_matches'] if m['order'] == order and m['index'] == index]
            require(len(matches) == 1 and matches[0]['first_id'] == sid and matches[0]['source_path'] == pointer['path'], 'paired source path ' + sid)
            source_index[sid] = pointer
            row = {
                "coefficient": iv(context, term["coefficient"]),
                "powers": term["powers"],
                "source_ids": [sid],
                "source_paths": ["source_index:" + sid],
                "parent_source_ids": parents,
                "error_classes": origin.get("error_classes", term.get("error_classes", [])),
                "uncertainty": origin.get("uncertainty", term.get("uncertainty", False)),
                "error_only": bool(term.get("error_only", False)),
            }
            terms.append(row)
        require(len(terms) == 16, "collar g5 sector count")
        g5[label] = terms

    chi0_id = "chi:Z0"
    chi1_id = "chi:Z1"
    chi0 = kchi["chi"]["chi"][0]
    chi1 = kchi["chi"]["chi"][1]
    for order, saved, sid in ((0, chi0, chi0_id), (1, chi1, chi1_id)):
        bound = binding["saved_chi_Z0_Z1"][str(order)]["coefficient"]
        require(saved == bound, "axis chi binding mismatch Z" + str(order))
        source_index[sid] = {
            "receipt": KCHI,
            "path": ["result", "output", "chi", "chi", order],
            "binding_receipt": BIND,
            "axis_potential_Z_only": True,
        }
    chi = {
        "coefficient": iv(context, chi0),
        "powers": [0, 0, 0, 0],
        "source_ids": [chi0_id],
        "source_paths": ["source_index:" + chi0_id],
        "error_classes": ["whole_Z_axis_potential_interval"],
        "uncertainty": True,
    }
    chi_Z = {
        "coefficient": iv(context, chi1),
        "powers": [0, 0, 0, 0],
        "source_ids": [chi1_id],
        "source_paths": ["source_index:" + chi1_id],
        "error_classes": ["whole_Z_axis_potential_interval"],
        "uncertainty": True,
    }

    leading: dict[str, list[dict]] = {"Z0": [], "Z1": []}
    for order in (0, 1):
        label = "Z" + str(order)
        terms = []
        for index, term in enumerate(swirl["F0_connection"]["leading"][str(order)]):
            sid = f"lead_F:{order}:{index}"
            origin = term["origin"]
            pointer = {
                "receipt": SWIRL,
                "path": ["result", "output", "F0_connection", "leading", str(order), index],
                "parent_source_ids": origin["source_ids"],
                "connection_path_correction": CORRECTION,
            }
            require(previous["source_index"].get(sid) == pointer, "leading source path " + sid)
            source_index[sid] = pointer
            terms.append(
                {
                    "coefficient": iv(context, term["coefficient"]),
                    "powers": term["powers"],
                    "source_ids": [sid],
                    "source_paths": ["source_index:" + sid],
                    "error_classes": origin.get("error_classes", []),
                    "uncertainty": origin.get("uncertainty", False),
                    "error_only": False,
                }
            )
        require(terms, "missing leading F sectors")
        leading[label] = terms

    certificate = {
        "same_source": True,
        "chi_is_axis_potential": True,
        "geometry_Z_independent": True,
        "g5_whole_phase_enclosure": True,
        "P_sign": 2,
        "eta_mass": context.mpf(5) / context.mpf(8),
        "q_domain": context.mpf([0, 1]),
    }
    hb = context.mpf([0, "0.5"])
    old_tail = old_receipt["actual_local_j2_tail"]
    result = leaf.build(
        context,
        hb,
        sc,
        chi,
        chi_Z,
        g5,
        leading,
        certificate,
        old_j2_tail=old_tail,
        source_family=family,
        source_atoms=["hb", "F0_anchor", "Pstar_squared", "Lambda"],
    )
    require(
        result["schema"] == COMPONENT
        and result["ledger_written"] is False
        and result["old_j2_tail_cancelled"] is False
        and result["old_errors_preserved"] is True,
        "leaf scope",
    )

    counts = {
        group + ":" + label: len(result[group][label])
        for group in ("finite_j2", "P_correlated_negative", "j_ge_3")
        for label in ("Z0", "Z1")
    }
    for group in ("finite_j2", "P_correlated_negative", "j_ge_3"):
        for label in ("Z0", "Z1"):
            for term in result[group][label]:
                require(all(sid in source_index for sid in term.get("source_ids", [])), "unindexed output provenance")
    for item in result['P_correlated_pairs']:
        parent, negative = item['D2g'], item['P']
        require(negative['coefficient']._mpi_ == (-parent['coefficient'])._mpi_
                and parent['powers'] == negative['powers']
                and parent['source_ids'] == negative['source_ids']
                and parent['error_classes'] == negative['error_classes']
                and item['symbolic_sum']['independent_boxed_interval_sum_performed'] is False,
                'correlated pair identity and provenance')
    require(all(t['error_only'] is True and t['symmetric_error'] is True
                and t['assigned_tail_value'] is False
                for label in ('Z0', 'Z1') for t in result['j_ge_3'][label]), 'tail is error only')

    return {
        "source_family": family,
        "source_logs": atoms,
        "source_index": source_index,
        "fixture": proof,
        "counts": counts,
        "input_counts": {label: len(rows) for label, rows in g5.items()},
        "result": enc(
            {
                "schema": result["schema"],
                "geometry": result["geometry"],
                "chi": result["chi"],
                "finite_j2": result["finite_j2"],
                "P_correlated_negative": result["P_correlated_negative"],
                "P_correlated_pairs": result["P_correlated_pairs"],
                "j_ge_3": result["j_ge_3"],
                "old_j2_tail": result["old_j2_tail"],
            }
        ),
        "old_tail_reference": {"receipt": OLD_TAIL, "sha256": hashes[OLD_TAIL]},
        "pair_reference": {"receipt": PAIR, "sha256": hashes[PAIR]},
        "compensation_reference": {"receipt": COMP, "sha256": hashes[COMP]},
        "old_errors_preserved": True,
        "old_errors_removed": False,
        "appended_to_ledger": False,
        "ledger_written": False,
        "complete_m15_computed": False,
        "terminal_closure_proved": False,
    }


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
    outcome = {
        "component": COMPONENT,
        "source_sha256": source_sha,
        "entry_sha256": entry_sha,
        "status": "failed",
    }
    try:
        outcome.update(status="returned", output=consume(hashes), error=None)
    except BaseException as exc:
        outcome["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    # Persist the complete sealed outcome exactly once before emitting status.
    result_sha = save(RESULT, outcome)
    print(json.dumps({"component": COMPONENT, "status": outcome["status"], "result_sha256": result_sha}))
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
