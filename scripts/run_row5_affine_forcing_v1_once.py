"""Run the source-bound row-5 affine forcing adapter exactly once.

The actual forcing receipt is used only as an interval/provenance packet.
The new leaf is imported after exclusive entry and after every input byte has
been authenticated.  The fixture is the only caller of ``evaluate``; the
actual packet is passed to ``build`` once and is never evaluated with W data.
"""

from __future__ import annotations

import gzip
import hashlib
import importlib.util
import json
import sys
import traceback
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from mpmath.ctx_iv import MPIntervalContext


ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
NAME = "row5_affine_forcing_v1"
DIR = ROOT / "tmp" / NAME
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entry.json", "result.json"))

LEAF = "src/openai_ns_reconstruction/row5_affine_forcing_v1.py"
FIXTURE = "src/openai_ns_reconstruction/row5_affine_forcing_fixture_v1.py"
RUNNER = "scripts/run_row5_affine_forcing_v1_once.py"
G5 = "experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz"
KCHI = "experiments/root_st073/n1_core_kchi_first_v1_receipt.json.gz"
BIND = "experiments/root_st073/n1_collar_axis_potential_binding_v1_receipt.json.gz"
SELECT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_selected_words.json"

REQUIRED = (LEAF, FIXTURE, RUNNER, G5, KCHI, BIND, SELECT)
PIN = {
    G5: "96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2",
    KCHI: "b0d654203e8c4acabdd0b61ca7b79c1b7050314b0c003a0c2f8757bde94973a0",
    BIND: "5859d5cb4bb82d372d8ac2f3ea1629e07c9787ffbaca149c99d1aba5edb8c6e8",
    SELECT: "0e04c87a011bfec74e0ceb0f822774a03495a08a839aaf4e378d21949736cf1e",
}

SOURCE_POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
MATRIX_KEYS = {
    "C0": ("5,1", "5,2", "5,3", "5,5"),
    "B1": ("5,1", "5,2", "5,3"),
}
COEFFICIENT_KEYS = (
    "coefficient",
    "coefficient_interval",
    "raw_interval",
    "directed_interval",
    "interval",
    "enclosure",
    "box",
    "value",
)


def sha(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def require(value: object, message: str) -> None:
    if not value:
        raise ValueError(message)


def save(path: Path, value: object) -> str:
    raw = (json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False) + "\n").encode("utf-8")
    with path.open("xb") as handle:
        handle.write(raw)
    return sha(raw)


def clone(value: Any) -> Any:
    """Copy receipt-shaped JSON without changing its exact interval tuples."""

    if isinstance(value, Mapping):
        return {key: clone(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [clone(item) for item in value]
    return value


def encode_intervals(value: Any) -> Any:
    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        return {
            "lower_exact_mpf_tuple": list(raw[0]),
            "upper_exact_mpf_tuple": list(raw[1]),
            "directed": True,
        }
    if isinstance(value, Mapping):
        return {str(key): encode_intervals(item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [encode_intervals(item) for item in value]
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    raise TypeError("cannot serialize result type " + type(value).__name__)


def interval(c: MPIntervalContext, value: Any, label: str = "interval") -> Any:
    """Decode an exact directed interval without midpoint conversion."""

    raw = getattr(value, "_mpi_", None)
    if raw is not None:
        result = c.make_mpf((raw[0], raw[1]))
    elif isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            require(
                "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value,
                label + " requires both exact endpoints",
            )
            lower = tuple(value["lower_exact_mpf_tuple"])
            upper = tuple(value["upper_exact_mpf_tuple"])
            require(len(lower) == len(upper) == 4, label + " endpoint tuples")
            require(all(type(item) is int for item in lower + upper), label + " endpoint integers")
            result = c.make_mpf((lower, upper))
        else:
            nested = next((value[key] for key in COEFFICIENT_KEYS if key in value), None)
            if nested is not None and nested is not value:
                return interval(c, nested, label + ".coefficient")
            lower = next((value[key] for key in ("lower", "lo", "left") if key in value), None)
            upper = next((value[key] for key in ("upper", "hi", "right") if key in value), None)
            require(lower is not None and upper is not None, label + " is not interval-compatible")
            result = c.make_mpf((lower, upper))
    elif isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)) and len(value) == 2:
        result = c.make_mpf((value[0], value[1]))
    else:
        result = c.mpf(value)
    require(getattr(result, "_mpi_", None) is not None, label + " is not directed")
    require(result.a <= result.b, label + " has reversed endpoints")
    return result


def module(path: str, name: str) -> Any:
    require(path in (LEAF, FIXTURE), "unexpected import path")
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    if spec is None or spec.loader is None:
        raise RuntimeError("new component unavailable: " + path)
    loaded = importlib.util.module_from_spec(spec)
    sys.modules[name] = loaded
    spec.loader.exec_module(loaded)
    return loaded


def preflight() -> tuple[dict[str, str], str, str]:
    require(not ENTRY.exists() and not RESULT.exists(), "component already entered")
    brief_raw = BRIEF.read_bytes()
    source_raw = SOURCE.read_bytes()
    brief = json.loads(brief_raw)
    manifest = json.loads(source_raw)
    require(brief["component"] == NAME and brief["phase"] == NAME + "_pending", "brief guard")
    require(brief["goal_active"] is True and brief["full_goal_complete"] is False, "goal guard")
    require(brief["one_entry_only"] is True and brief["do_not_reenter"] is True, "exclusive-entry guard")
    checkpoint = brief["source_checkpoint"]
    require(checkpoint["sha256"] == sha(source_raw), "source checkpoint hash")
    require(Path(checkpoint["path"]).resolve() == SOURCE.resolve(), "source checkpoint path")
    require(Path(checkpoint["preimage_path"]).resolve() == Path("C:/Users/z5242/Documents/Codex/row5_affine_pre.json").resolve(), "preimage path")
    require(sha(Path(checkpoint["preimage_path"]).read_bytes()) == checkpoint["preimage_sha256"], "preimage hash")
    require(manifest["component"] == NAME and manifest["one_entry_only"] is True, "manifest guard")
    require(len(manifest["files"]) == len(REQUIRED), "manifest file count")
    require({row["path"] for row in manifest["files"]} == set(REQUIRED), "manifest file set")
    hashes: dict[str, str] = {}
    for row in manifest["files"]:
        path = row["path"]
        require(path not in hashes, "duplicate manifest path: " + path)
        digest = sha((ROOT / path).read_bytes())
        require(digest == row["sha256"], "changed input: " + path)
        if path in PIN:
            require(digest == PIN[path], "changed pinned input: " + path)
        hashes[path] = digest
    return hashes, sha(source_raw), sha(brief_raw)


def _receipt_output(receipt: Mapping[str, Any], path: str) -> dict[str, Any]:
    result = receipt.get("result")
    require(isinstance(result, Mapping), "missing result: " + path)
    require(result.get("status") == "returned", "receipt not returned: " + path)
    output = result.get("output")
    require(isinstance(output, Mapping), "missing output: " + path)
    return dict(output)


def _as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return list(value)
    return [value]


def _unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    for value in values:
        if value not in result:
            result.append(value)
    return result


def _metadata_values(raw: Mapping[str, Any], origin: Mapping[str, Any], names: Sequence[str]) -> list[Any]:
    values: list[Any] = []
    for container in (raw, origin):
        for name in names:
            if name in container:
                values.extend(_as_list(container[name]))
    return _unique(values)


def _metadata_flag(raw: Mapping[str, Any], origin: Mapping[str, Any], name: str, default: bool = False) -> bool:
    if name in raw:
        value = raw[name]
    elif name in origin:
        value = origin[name]
    else:
        value = default
    require(type(value) is bool, "sector " + name + " must be boolean")
    return value


def _term_coefficient(raw: Mapping[str, Any], label: str) -> Any:
    for key in COEFFICIENT_KEYS:
        if key in raw:
            return clone(raw[key])
    raise KeyError(label + " is missing coefficient")


def _term_powers(raw: Mapping[str, Any], label: str) -> Any:
    for key in ("powers", "formal_powers", "source_powers"):
        if key in raw:
            return clone(raw[key])
    origin = raw.get("origin")
    if isinstance(origin, Mapping):
        for key in ("powers", "formal_powers", "source_powers"):
            if key in origin:
                return clone(origin[key])
    raise KeyError(label + ".powers is required")


def _order_value(mapping: Mapping[Any, Any], order: int, label: str) -> Any:
    for key in (str(order), order, "Z" + str(order), "ordinary_Z" + str(order)):
        if key in mapping:
            return mapping[key]
    raise KeyError(label + " missing order " + str(order))


def _adapt_term(
    raw: Any,
    *,
    family: str,
    z_label: str,
    key: str,
    index: int,
    receipt_path: list[Any],
    source_index: dict[str, Any],
) -> dict[str, Any]:
    require(isinstance(raw, Mapping), "actual sector must be a mapping")
    original = clone(dict(raw))
    origin_value = raw.get("origin", {})
    origin = clone(origin_value)
    origin_mapping = origin if isinstance(origin, Mapping) else {}
    parents = _metadata_values(
        raw,
        origin_mapping,
        ("source_ids", "source_id", "sources", "source", "input_id"),
    )
    uncertainty = _metadata_flag(raw, origin_mapping, "uncertainty", False)
    error_classes = _metadata_values(raw, origin_mapping, ("error_classes", "error_class", "uncertainty_class"))
    source_error_classes = _metadata_values(raw, origin_mapping, ("source_error_classes", "source_error_class"))
    if uncertainty:
        require(bool(error_classes), "uncertain actual sector has no error classes")
    powers = _term_powers(raw, family + "." + z_label + "." + key + "." + str(index))
    require(isinstance(powers, (list, tuple)) and len(powers) == 4
            and all(type(power) is int for power in powers), "actual formal sector powers")
    if family == "g5":
        require(list(powers[1:]) == [1, 0, -1], "raw g5 F0/Lambda basis")
    _term_coefficient(raw, family + "." + z_label + "." + key + "." + str(index))
    stable_id = "row5_affine:" + family + ":" + z_label + ":" + key + ":" + str(index)
    require(stable_id not in source_index, "duplicate stable source id: " + stable_id)
    row = clone(dict(raw))
    if "powers" not in row:
        row["powers"] = clone(powers)
    if "formal_powers" not in row:
        row["formal_powers"] = clone(powers)
    row["source_ids"] = _unique(list(parents) + [stable_id])
    row["parent_source_ids"] = list(parents)
    row["receipt_path"] = clone(receipt_path)
    old_paths = _as_list(raw.get("source_paths", []))
    if "source_path" in raw:
        old_paths.append(clone(raw["source_path"]))
    row["source_paths"] = _unique(old_paths + [clone(receipt_path)])
    if "source_path" not in row:
        row["source_path"] = clone(receipt_path)
    row["original_origin"] = clone(origin_value)
    row["original_source_term"] = original
    row["source_term"] = original
    row["uncertainty"] = uncertainty
    row["error_classes"] = list(error_classes)
    if source_error_classes:
        row["source_error_classes"] = list(source_error_classes)
    source_index[stable_id] = {
        "receipt": G5,
        "path": clone(receipt_path),
        "family": family,
        "source_ids": list(parents) + [stable_id],
        "parent_source_ids": list(parents),
        "origin": clone(origin_value),
        "powers": clone(powers),
        "coefficient": _term_coefficient(raw, family + "." + z_label + "." + key + "." + str(index)),
        "uncertainty": uncertainty,
        "error_classes": list(error_classes),
        "source_error_classes": list(source_error_classes),
    }
    return row


def _build_packet(forcing: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    g = forcing.get("g")
    matrix = forcing.get("matrix")
    require(isinstance(g, Mapping), "actual forcing output g is missing")
    require(isinstance(matrix, Mapping), "actual forcing output matrix is missing")
    source_index: dict[str, Any] = {}
    packet: dict[str, Any] = {"g5": {"Z0": [], "Z1": []}, "C0": {"Z0": {}, "Z1": {}}, "B1": {"Z0": {}, "Z1": {}}}
    for order in (0, 1):
        label = "Z" + str(order)
        saved = _order_value(g, order, "actual forcing g")
        require(isinstance(saved, Mapping), "actual forcing g order is not a mapping")
        terms = saved.get("5") if "5" in saved else saved.get(5)
        require(isinstance(terms, Sequence) and not isinstance(terms, (str, bytes, bytearray)), "actual g5 row is not a sequence")
        for index, raw in enumerate(terms):
            path = ["result", "output", "g", str(order), "5", index]
            packet["g5"][label].append(
                _adapt_term(raw, family="g5", z_label=label, key="g5", index=index, receipt_path=path, source_index=source_index)
            )
        require(packet["g5"][label], "actual g5 row is empty: " + label)
    for name, keys in MATRIX_KEYS.items():
        section = matrix.get(name)
        require(isinstance(section, Mapping), "actual forcing matrix section missing: " + name)
        for order in (0, 1):
            label = "Z" + str(order)
            saved = _order_value(section, order, "actual forcing matrix " + name)
            require(isinstance(saved, Mapping), "actual matrix order is not a mapping: " + name + "/" + label)
            for key in keys:
                terms = saved.get(key)
                if terms is None and isinstance(key, str) and "," in key:
                    terms = saved.get(tuple(int(item) for item in key.split(",")))
                require(isinstance(terms, Sequence) and not isinstance(terms, (str, bytes, bytearray)), "actual matrix row is not a sequence: " + name + "/" + label + "/" + key)
                packet[name][label][key] = []
                for index, raw in enumerate(terms):
                    path = ["result", "output", "matrix", name, str(order), key, index]
                    packet[name][label][key].append(
                        _adapt_term(raw, family=name, z_label=label, key=key, index=index, receipt_path=path, source_index=source_index)
                    )
                require(packet[name][label][key], "actual matrix row is empty: " + name + "/" + label + "/" + key)
    return packet, source_index, {"g": clone(g), "matrix": clone(matrix)}


def _axis_term(
    c: MPIntervalContext,
    kchi: Mapping[str, Any],
    binding: Mapping[str, Any],
    order: int,
    source_index: dict[str, Any],
) -> dict[str, Any]:
    chi_section = kchi.get("chi")
    require(isinstance(chi_section, Mapping), "axis KCHI chi section missing")
    chi_values = chi_section.get("chi")
    require(isinstance(chi_values, Sequence) and not isinstance(chi_values, (str, bytes, bytearray)), "axis KCHI chi rows missing")
    raw_chi = chi_values[order]
    binding_rows = binding.get("saved_chi_Z0_Z1")
    require(isinstance(binding_rows, Mapping), "axis BIND saved jet missing")
    bound = binding_rows.get(str(order))
    require(isinstance(bound, Mapping) and "coefficient" in bound, "axis BIND coefficient missing")
    chi_iv = interval(c, raw_chi, "axis KCHI Z" + str(order))
    bound_iv = interval(c, bound["coefficient"], "axis BIND Z" + str(order))
    require(chi_iv._mpi_ == bound_iv._mpi_, "axis exact interval mismatch Z" + str(order))
    label = "Z" + str(order)
    stable_id = "axis_chi:" + label
    common_id = "axis_chi:whole_Z"
    path = ["result", "output", "chi", "chi", order]
    error_classes = ["whole_Z_axis_potential_interval"]
    if common_id not in source_index:
        source_index[common_id] = {
            "receipt": KCHI,
            "binding_receipt": BIND,
            "scope": "whole-Z axis source jet",
            "source_ids": [common_id],
            "error_classes": list(error_classes),
            "uncertainty": True,
            "whole_Z": True,
        }
    term = {
        "coefficient": clone(raw_chi),
        "powers": [0, 0, 0, 0],
        "formal_powers": [0, 0, 0, 0],
        "power_basis": list(SOURCE_POWER_BASIS),
        "source_ids": [common_id, stable_id],
        "source_paths": [clone(path)],
        "source_path": clone(path),
        "receipt_path": clone(path),
        "origin": {
            "source_ids": [common_id, stable_id],
            "error_classes": list(error_classes),
            "uncertainty": True,
            "whole_Z": True,
        },
        "original_origin": {"source_ids": [common_id, stable_id], "error_classes": list(error_classes), "uncertainty": True, "whole_Z": True},
        "original_source_term": {"coefficient": clone(raw_chi), "powers": [0, 0, 0, 0]},
        "uncertainty": True,
        "error_classes": list(error_classes),
        "whole_Z": True,
        "midpoint_values_used": False,
    }
    require(stable_id not in source_index, "duplicate axis source id: " + stable_id)
    source_index[stable_id] = {
        "receipt": KCHI,
        "binding_receipt": BIND,
        "path": clone(path),
        "source_ids": [common_id, stable_id],
        "coefficient": clone(raw_chi),
        "binding_coefficient": clone(bound["coefficient"]),
        "error_classes": list(error_classes),
        "uncertainty": True,
        "whole_Z": True,
    }
    return term


def consume(hashes: dict[str, str]) -> dict[str, Any]:
    # Authenticate every byte before any gzip decompression or new-module import.
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "changed after entry: " + path)
        raw[path] = data

    leaf = module(LEAF, NAME + "_leaf")
    fixture = module(FIXTURE, NAME + "_fixture")
    c = MPIntervalContext()
    c.prec = 1800
    proof = fixture.run_fixture(leaf.build, leaf.evaluate, c)
    require(isinstance(proof, Mapping), "fixture result must be a mapping")
    require(proof.get("passed") is True, "fresh synthetic fixture did not pass")
    if isinstance(proof.get("checks"), Mapping):
        require(all(value is True for value in proof["checks"].values()), "fresh synthetic fixture has failed checks")

    # Only now decode the authenticated gzip receipts.  No ancestor module or
    # ancestor fixture is imported or called.
    forcing_receipt = json.loads(gzip.decompress(raw[G5]))
    kchi_receipt = json.loads(gzip.decompress(raw[KCHI]))
    bind_receipt = json.loads(gzip.decompress(raw[BIND]))
    selected = json.loads(raw[SELECT])
    forcing = _receipt_output(forcing_receipt, G5)
    kchi = _receipt_output(kchi_receipt, KCHI)
    binding = _receipt_output(bind_receipt, BIND)

    family = forcing.get("source_family")
    source_logs = forcing.get("source_logs")
    require(family is not None and source_logs is not None, "actual forcing source family/logs missing")
    for label, output in (("KCHI", kchi), ("BIND", binding)):
        require(output.get("source_family") == family, "source family mismatch: " + label)
        require(output.get("source_logs") == source_logs, "source logs mismatch: " + label)
    require(selected.get("source_family") == family, "selected source family mismatch")
    if "source_logs" in selected:
        require(selected["source_logs"] == source_logs, "selected source logs mismatch")

    metadata = forcing.get("calculation_metadata")
    require(isinstance(metadata, Mapping), "actual forcing calculation metadata missing")
    require(metadata.get("same_source_dynamic_C0_and_B1_Z0_to_Z2") is True, "same-source C0/B1 metadata guard")
    require(metadata.get("F0_relative_derivative_applied_once") is True, "F0 derivative metadata guard")
    require(metadata.get("g5_formula") == "2*epsilon*(-axial_second(F))", "historical raw g5 formula guard")
    require(metadata.get("g5_anchor_rule") == "F(Z+h)/F0(Z) jet at each current Z; no extra D_F", "historical g5 derivative guard")
    checks = binding.get("source_checks")
    require(isinstance(checks, Mapping) and checks.get("one_query_jet_shared_in_integration") is True, "axis source-jet binding")
    require(binding.get("phase_shear_cutoff_identified_with_axis_potential") is False, "phase shear substituted for axis potential")

    source_index: dict[str, Any] = {}
    packet, packet_index, raw_selected = _build_packet(forcing)
    source_index.update(packet_index)
    chi0 = _axis_term(c, kchi, binding, 0, source_index)
    chi1 = _axis_term(c, kchi, binding, 1, source_index)
    expression = leaf.build(c, packet, chi0, chi1)
    require(isinstance(expression, Mapping), "row-5 build result must be a mapping")

    # The adapter is source-bound and deliberately does not assert physical
    # x/rho ownership, a live W owner, a physical h5 value, a cancellation, or
    # any ledger/m15 result.
    require(expression.get("physical_h5_materialized") is False, "physical h5 was materialized")
    require(expression.get("h5_is_physical_function") is False, "physical h5 function was asserted")
    require(expression.get("source_owner_proved") is False, "source owner was asserted")
    require(expression.get("same_state_owner_proved") is False, "state owner was asserted")
    require(expression.get("axis_same_owner_proved") is False, "axis live owner was asserted")
    require(expression.get("physical_dx_applied") is False, "physical dx was applied")
    require(expression.get("rho_x2_chain_rule_applied") is False, "rho=x^2 chain rule was applied")
    require(expression.get("interval_sum_in_build") is False, "source intervals were summed in build")
    scope = expression.get("scope")
    require(isinstance(scope, Mapping), "row-5 scope metadata missing")
    for key in ("ledger_written", "m15", "physical_h5", "state_materialized", "ns_validation"):
        require(scope.get(key) is False, "forbidden scope flag asserted: " + key)

    actual_inputs = {
        "receipt": G5,
        "receipt_sha256": hashes[G5],
        "source_family": clone(family),
        "source_logs": clone(source_logs),
        "g": clone(forcing["g"]),
        "matrix": clone(forcing["matrix"]),
        "calculation_metadata": clone(metadata),
        "selected_packet_sections": raw_selected,
        "axis": {
            "KCHI": {
                "Z0": clone(kchi["chi"]["chi"][0]),
                "Z1": clone(kchi["chi"]["chi"][1]),
            },
            "BIND": {
                "Z0": clone(binding["saved_chi_Z0_Z1"]["0"]),
                "Z1": clone(binding["saved_chi_Z0_Z1"]["1"]),
            },
        },
    }
    return {
        "source_family": clone(family),
        "source_logs": clone(source_logs),
        "source_index": source_index,
        "fixture": encode_intervals(proof),
        "actual_source_inputs": actual_inputs,
        "packet": packet,
        "expression": encode_intervals(expression),
        "input_receipt_sha256": {path: hashes[path] for path in (G5, KCHI, BIND, SELECT)},
        "axis_source_ids": {
            "common": "axis_chi:whole_Z",
            "Z0": list(chi0["source_ids"]),
            "Z1": list(chi1["source_ids"]),
        },
        "scope": {
            "synthetic_fixture_only": True,
            "actual_evaluate_called": False,
            "actual_build_calls": 1,
            "W_values_supplied": False,
            "physical_h5_values": False,
            "physical_h5_owner_proved": False,
            "W_owner_proved": False,
            "fixed_x_rho_x2_relationship": "caller_condition",
            "geometry": "caller_condition",
            "live_owner": False,
            "ledger_written": False,
            "m15": False,
            "cancellation": False,
            "ancestor_called": False,
            "ancestor_imports": False,
            "midpoint_values_used": False,
        },
    }


def main() -> int:
    if not __debug__:
        raise RuntimeError("optimized mode disables required preflight checks")
    hashes, source_sha, brief_sha = preflight()
    entry_sha = save(
        ENTRY,
        {
            "component": NAME,
            "source_sha256": source_sha,
            "brief_sha256": brief_sha,
            "input_hashes": hashes,
            "entered": True,
            "ancestor_called": False,
            "raw_gzip_authenticated_before_decompression": True,
        },
    )
    outcome: dict[str, Any] = {
        "component": NAME,
        "source_sha256": source_sha,
        "entry_sha256": entry_sha,
        "status": "failed",
        "raw_gzip_authenticated_before_decompression": True,
        "ancestor_called": False,
        "actual_evaluate_called": False,
        "fixture": None,
        "expression": None,
    }
    try:
        payload = consume(hashes)
        outcome.update(
            status="returned",
            output=payload,
            fixture=payload.get("fixture"),
            expression=payload.get("expression"),
            error=None,
        )
    except BaseException as exc:
        outcome["error"] = {
            "type": type(exc).__name__,
            "message": str(exc),
            "traceback": traceback.format_exc(),
        }
    result_sha = save(RESULT, outcome)
    print(
        json.dumps(
            {
                "component": NAME,
                "status": outcome["status"],
                "entry_sha256": entry_sha,
                "result_sha256": result_sha,
                "error": (outcome.get("error") or {}).get("message"),
            }
        )
    )
    return 0 if outcome["status"] == "returned" else 1


if __name__ == "__main__":
    raise SystemExit(main())
