"""Run the source-bound row-5 first-word product consumer exactly once.

All input bytes are authenticated before any receipt is decompressed or new
code is imported.  The consumer receives the saved affine expression and a
finite first-word state assembled from the authenticated join/Green outputs.
It produces a product ledger only; formal atoms, physical ``dx``, physical
``h5``, a live owner, and the residual ``R = W - W_finite`` are not supplied
or proved here.  The v1 failure archive is authenticated for provenance by
raw SHA only; it is never decompressed, imported, or re-entered.
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

ROOT = Path(__file__).resolve().parents[1]
BRIEF = ROOT.parent.parent / "RESUME_BRIEF_2026_10_16.json"
NAME = "row5_firstword_product_v4"
DIR = ROOT / "tmp" / NAME
SOURCE, ENTRY, RESULT = (DIR / name for name in ("source.json", "entry.json", "result.json"))

LEAF = "src/openai_ns_reconstruction/row5_firstword_product_v4.py"
FIXTURE = "src/openai_ns_reconstruction/row5_firstword_product_fixture_v4.py"
RUNNER = "scripts/run_row5_firstword_product_v4_once.py"
AFFINE = "experiments/root_st073/row5_affine_forcing_v1_receipt.json.gz"
JOIN = "experiments/root_st073/n1_core_collar_low_join_v2_receipt.json.gz"
GREEN = "experiments/root_st073/n1_actual_collar_green_v3_receipt.json.gz"
FORCING = "experiments/root_st073/n1_actual_collar_forcing_v1_receipt.json.gz"
FAILURE_V1 = "experiments/root_st073/row5_firstword_product_v1_failure.json.gz"
FAILURE_V2 = "experiments/root_st073/row5_firstword_product_v2_failure.json.gz"
FAILURE_V3_PREFLIGHT = "experiments/root_st073/row5_firstword_product_v3_preflight_failure.json.gz"
REQUIRED = (LEAF, FIXTURE, RUNNER, AFFINE, JOIN, GREEN, FORCING, FAILURE_V1, FAILURE_V2, FAILURE_V3_PREFLIGHT)
PIN = {
    AFFINE: "b64d14b59b178502e562fe593d79419cc2e21200f80c9bf0006905d1b6e72515",
    JOIN: "a2e25ccec22f97dc99d1eb80609f0d806d4276e7a75e763ca882d90532e7ce67",
    GREEN: "fad35594ccf5ef81392f94f97ae72b3ffb94fcaa554438fd5702ef2c4343c4b5",
    FORCING: "96e8a3b64140fb2c111e8117951e96674dfa3c2bc7ff41ebe34448dccbbd8be2",
    FAILURE_V1: "331d752260480dcd8f9b80cddc697b3b8c8b85c9e0bdd625bf1cbd20966e0a97",
    FAILURE_V2: "b729f77899bdea1e4583feef6126d98f024048f1b291a06b6a71807aa1088c67",
    FAILURE_V3_PREFLIGHT: "45af50e1cf504c61fb065a53dcc170000f3399b2bc91f62a69efc28954a9c462",
}

POWER_BASIS = ("hb", "F0_anchor", "Pstar_squared", "Lambda")
FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
STATE_COMPONENTS = (1, 2, 3, 5)
ROWS = ("h5", "h5_Z1")
COEFFICIENT_KEYS = (
    "coefficient_interval",
    "coefficient",
    "directed_interval",
    "enclosure",
    "interval",
    "raw_interval",
    "value",
)
POWER_KEYS = ("formal_powers", "powers", "source_powers")
ID_KEYS = ("source_ids", "source_id", "state_ids", "state_id", "input_id", "term_id")
PATH_KEYS = (
    "source_paths",
    "source_path",
    "state_paths",
    "state_path",
    "input_paths",
    "input_path",
    "paths",
    "path",
    "parent_receipt_path",
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


def decode_intervals(c: MPIntervalContext, value: Any) -> Any:
    """Decode only exact directed-MPF records; never use a midpoint/float."""

    if isinstance(value, Mapping):
        if "lower_exact_mpf_tuple" in value or "upper_exact_mpf_tuple" in value:
            require(
                "lower_exact_mpf_tuple" in value and "upper_exact_mpf_tuple" in value,
                "directed interval requires both exact endpoints",
            )
            lower = tuple(value["lower_exact_mpf_tuple"])
            upper = tuple(value["upper_exact_mpf_tuple"])
            require(len(lower) == len(upper) == 4, "directed endpoint tuple length")
            require(all(type(item) is int for item in lower + upper), "directed endpoint tuple integers")
            interval = c.make_mpf((lower, upper))
            require(interval.a <= interval.b, "directed interval has reversed endpoints")
            return interval
        return {key: decode_intervals(c, item) for key, item in value.items()}
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return [decode_intervals(c, item) for item in value]
    return value


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
    require(
        brief["component"] == NAME
        and brief["phase"] == NAME + "_pending"
        and brief["goal_active"] is True
        and brief["full_goal_complete"] is False,
        "brief guard",
    )
    require(brief["one_entry_only"] is True and brief["do_not_reenter"] is True, "exclusive-entry guard")
    checkpoint = brief["source_checkpoint"]
    require(Path(checkpoint["path"]).resolve() == SOURCE.resolve(), "source checkpoint path")
    require(checkpoint["sha256"] == sha(source_raw), "source checkpoint hash")
    require(Path(checkpoint["preimage_path"]).resolve() == Path("C:/Users/z5242/Documents/Codex/row5_firstword_v4_pre.json").resolve(), "preimage path")
    require(sha(Path(checkpoint["preimage_path"]).read_bytes()) == checkpoint["preimage_sha256"], "preimage hash")
    require(manifest["component"] == NAME and manifest["one_entry_only"] is True, "manifest guard")
    require(manifest.get("supersedes") == "row5_firstword_product_v3", "v3 supersession guard")
    failure_archives = manifest.get("failure_archives")
    require(isinstance(failure_archives, Sequence) and len(failure_archives) == 3, "failure archive manifest guard")
    failure_manifest = {item.get("path"): item for item in failure_archives if isinstance(item, Mapping)}
    require(set(failure_manifest) == {FAILURE_V1, FAILURE_V2, FAILURE_V3_PREFLIGHT}, "failure archive path guard")
    require(all(failure_manifest[path].get("sha256") == PIN[path] for path in (FAILURE_V1, FAILURE_V2, FAILURE_V3_PREFLIGHT)), "failure archive pin guard")
    require(all(failure_manifest[path].get("raw_sha256_only") is True for path in (FAILURE_V1, FAILURE_V2, FAILURE_V3_PREFLIGHT)), "failure archive raw-only guard")
    previous_checkpoint = manifest.get("previous_source_checkpoint")
    require(isinstance(previous_checkpoint, Mapping), "previous source checkpoint missing")
    require(previous_checkpoint.get("sha256") == "d587382bd496638123740e42a3bc4fd6b07a4f86d0d846091ac718d7b26db3f1", "previous source SHA mismatch")
    previous_source = Path(previous_checkpoint["path"])
    previous_preimage = Path(previous_checkpoint["preimage_path"])
    require(previous_preimage.resolve() == Path("C:/Users/z5242/Documents/Codex/row5_firstword_v3_pre.json").resolve(), "unexpected previous preimage path")
    previous_source_raw = previous_source.read_bytes()
    require(sha(previous_source_raw) == previous_checkpoint["sha256"], "previous source checkpoint changed")
    require(sha(previous_preimage.read_bytes()) == previous_checkpoint["preimage_sha256"], "previous preimage changed")
    previous_manifest = json.loads(previous_source_raw)
    require(previous_manifest.get("component") == "row5_firstword_product_v3", "unexpected superseded source")
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


def receipt_output(receipt: Mapping[str, Any], path: str) -> dict[str, Any]:
    result = receipt.get("result")
    require(isinstance(result, Mapping), "missing result: " + path)
    require(result.get("status") == "returned", "receipt not returned: " + path)
    output = result.get("output")
    require(isinstance(output, Mapping), "missing output: " + path)
    return dict(output)


def as_list(value: Any) -> list[Any]:
    if value is None:
        return []
    if isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)):
        return list(value)
    return [value]


def unique(values: Sequence[Any]) -> list[Any]:
    result: list[Any] = []
    for value in values:
        if value not in result:
            result.append(value)
    return result


def metadata_values(raw: Mapping[str, Any], names: Sequence[str]) -> list[Any]:
    values: list[Any] = []
    for name in names:
        if name in raw:
            values.extend(clone(item) for item in as_list(raw[name]))
    for wrapper in ("origin", "metadata", "provenance", "original_term"):
        nested = raw.get(wrapper)
        if isinstance(nested, Mapping):
            for name in names:
                if name in nested:
                    values.extend(clone(item) for item in as_list(nested[name]))
    return unique(values)


def path_values(raw: Mapping[str, Any]) -> list[Any]:
    values: list[Any] = []
    singular = ("source_path", "state_path", "input_path", "path", "parent_receipt_path")
    plural = ("source_paths", "state_paths", "input_paths", "paths")
    for key in singular:
        if key in raw:
            values.append(clone(raw[key]))
    for key in plural:
        if key in raw:
            values.extend(clone(item) for item in as_list(raw[key]))
    for wrapper in ("origin", "metadata", "provenance", "original_term"):
        nested = raw.get(wrapper)
        if not isinstance(nested, Mapping):
            continue
        for key in singular:
            if key in nested:
                values.append(clone(nested[key]))
        for key in plural:
            if key in nested:
                values.extend(clone(item) for item in as_list(nested[key]))
    return unique(values)


def first_value(raw: Mapping[str, Any], keys: Sequence[str], label: str) -> Any:
    for key in keys:
        if key in raw:
            return raw[key]
    raise KeyError(label)


def powers(raw: Mapping[str, Any], label: str) -> Any:
    value = first_value(raw, POWER_KEYS, label + ".formal_powers")
    if isinstance(value, Mapping):
        value = [value.get(name) for name in POWER_BASIS]
    require(isinstance(value, Sequence) and not isinstance(value, (str, bytes, bytearray)), label + " powers")
    require(len(value) == 4 and all(type(item) is int for item in value), label + " powers must be four integers")
    return clone(value)


def coefficient(raw: Mapping[str, Any], label: str) -> Any:
    return clone(first_value(raw, COEFFICIENT_KEYS, label + ".coefficient"))


def _state_path(receipt: str, node_section: str, order: int, field: str, index: int) -> list[Any]:
    if node_section == "finite_W1_Z2":
        require(order == 2, "finite_W1_Z2 is an ordinary Z2 section")
        return ["result", "output", "result", node_section, field, index]
    return ["result", "output", "result", node_section, str(order), field, index]


def adapt_state_term(
    c: MPIntervalContext,
    raw: Any,
    *,
    receipt: str,
    node_section: str,
    component: int,
    order: int,
    field: str,
    index: int,
    state_index: dict[str, Any],
) -> dict[str, Any]:
    require(isinstance(raw, Mapping), "finite state term must be a mapping")
    label = f"{receipt}:{node_section}:W{component}:Z{order}:{field}:{index}"
    original = clone(dict(raw))
    original_coefficient = coefficient(raw, label)
    original_powers = powers(raw, label)
    parent_ids = metadata_values(raw, ID_KEYS)
    require(parent_ids, label + " has no parent source IDs")
    receipt_token = Path(receipt).name.removesuffix(".json.gz")
    stable_id = f"row5_firstword:state:{receipt_token}:W{component}:Z{order}:{field}:{index}"
    require(stable_id not in state_index, "duplicate stable state source ID: " + stable_id)
    json_path = _state_path(receipt, node_section, order, field, index)
    row = clone(dict(raw))
    # Keep the exact receipt term and its error metadata nested before the
    # directed interval objects are decoded for the consumer.
    row["original_receipt_term"] = original
    row["original_source_term"] = original
    row["original_coefficient"] = original_coefficient
    row["original_formal_powers"] = original_powers
    row["parent_source_ids"] = list(parent_ids)
    row["state_source_id"] = stable_id
    row["receipt_path"] = clone(json_path)
    old_paths = path_values(raw)
    row["source_paths"] = unique(old_paths + [clone(json_path)])
    row["source_path"] = clone(json_path) if "source_path" not in row else row["source_path"]
    row["source_ids"] = unique(list(parent_ids) + [stable_id])
    row["parent_receipt_path"] = clone(json_path)
    # Leaf._interval accepts an MP interval; decoding these top-level fields
    # leaves original_receipt_term/original_source_term untouched.
    for key in COEFFICIENT_KEYS:
        if key in row:
            row[key] = decode_intervals(c, row[key])
    state_index[stable_id] = {
        "receipt": receipt,
        "json_path": clone(json_path),
        "state_component": component,
        "ordinary_Z_order": order,
        "field": field,
        "source_id": stable_id,
        "parent_source_ids": list(parent_ids),
        "source_ids": list(parent_ids) + [stable_id],
        "origin": clone(raw.get("origin")),
        "powers": original_powers,
        "coefficient": original_coefficient,
        "error_classes": clone(raw.get("error_classes", raw.get("origin", {}).get("error_classes", []))) if isinstance(raw.get("origin", {}), Mapping) else clone(raw.get("error_classes", [])),
        "uncertainty": clone(raw.get("uncertainty", raw.get("origin", {}).get("uncertainty"))) if isinstance(raw.get("origin", {}), Mapping) else clone(raw.get("uncertainty")),
        "original_term": original,
    }
    return row


def source_index_for_expression(expression: Mapping[str, Any]) -> dict[str, Any]:
    rows = expression.get("rows")
    require(isinstance(rows, Mapping), "affine expression rows missing")
    result: dict[str, Any] = {}
    for row_name in ROWS:
        terms = rows.get(row_name)
        require(isinstance(terms, Sequence) and not isinstance(terms, (str, bytes, bytearray)), "affine expression row missing: " + row_name)
        for index, raw in enumerate(terms):
            require(isinstance(raw, Mapping), "affine expression term is not a mapping")
            term_powers = powers(raw, "expression." + row_name + "." + str(index))
            term_coefficient = coefficient(raw, "expression." + row_name + "." + str(index))
            ids = metadata_values(raw, ID_KEYS)
            paths = path_values(raw)
            require(ids, "affine expression term has no source IDs")
            require(paths, "affine expression term has no source paths")
            stable_id = f"row5_firstword:expression:{row_name}:{index}"
            result[stable_id] = {
                "source_id": stable_id,
                "row": row_name,
                "index": index,
                "parent_source_ids": ids,
                "source_paths": paths,
                "powers": term_powers,
                "coefficient": term_coefficient,
                "origin": clone(raw.get("origin")),
                "original_term": clone(dict(raw)),
            }
    return result


def state_terms(
    c: MPIntervalContext,
    source: Mapping[str, Any],
    *,
    receipt: str,
    node_section: str,
    component: int,
    order: int,
    field: str,
    state_index: dict[str, Any],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    section = source.get(node_section)
    require(isinstance(section, Mapping), "finite state section missing: " + node_section)
    if node_section == "finite_W1_Z2":
        require(order == 2, "finite_W1_Z2 is an ordinary Z2 section")
        order_section = section
    else:
        order_section = section.get(str(order), section.get(order))
    require(isinstance(order_section, Mapping), f"{node_section} Z{order} missing")
    raw_terms = order_section.get(field)
    require(isinstance(raw_terms, Sequence) and not isinstance(raw_terms, (str, bytes, bytearray)), f"{node_section} Z{order} {field} missing")
    require(raw_terms, f"{node_section} Z{order} {field} empty")
    adapted: list[dict[str, Any]] = []
    raw_index: list[dict[str, Any]] = []
    for index, raw in enumerate(raw_terms):
        adapted.append(
            adapt_state_term(
                c,
                raw,
                receipt=receipt,
                node_section=node_section,
                component=component,
                order=order,
                field=field,
                index=index,
                state_index=state_index,
            )
        )
        raw_index.append(clone(raw))
    return adapted, raw_index


def build_state(c: MPIntervalContext, join: Mapping[str, Any], green: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    join_result = join.get("result")
    green_result = green.get("result")
    require(isinstance(join_result, Mapping) and isinstance(green_result, Mapping), "finite producer result missing")
    state: dict[str, dict[str, list[dict[str, Any]]]] = {str(component): {} for component in STATE_COMPONENTS}
    index: dict[str, Any] = {}
    counts: dict[str, dict[str, dict[str, int]]] = {str(component): {} for component in STATE_COMPONENTS}

    source_for = {
        1: "F1",
        2: "Uz1",
        3: "K1",
        5: "d_x_F1",
    }
    for component in (1, 2, 3):
        field = source_for[component]
        for order in (0, 1):
            adapted, raw = state_terms(
                c,
                join_result,
                receipt=JOIN,
                node_section="finite_sum",
                component=component,
                order=order,
                field=field,
                state_index=index,
            )
            state[str(component)][f"Z{order}"] = adapted
            counts[str(component)][f"Z{order}"] = {field: len(adapted)}
        adapted, raw = state_terms(
            c,
            green_result,
            receipt=GREEN,
            node_section="finite_W1_Z2",
            component=component,
            order=2,
            field=field,
            state_index=index,
        )
        state[str(component)]["Z2"] = adapted
        counts[str(component)]["Z2"] = {field: len(adapted)}

    # Green publishes the complete six-field W1/Z2 finite transport; the
    # Uz1/K1 fields are the corresponding W2/W3 Z2 state components.  W5/Z2
    # is intentionally absent because the affine row-5 expression does not
    # require that derivative.
    for order in (0, 1):
        field = source_for[5]
        adapted, raw = state_terms(
            c,
            join_result,
            receipt=JOIN,
            node_section="finite_sum",
            component=5,
            order=order,
            field=field,
            state_index=index,
        )
        state["5"][f"Z{order}"] = adapted
        counts["5"][f"Z{order}"] = {field: len(adapted)}

    return state, {"terms": index, "counts": counts}, counts


def validate_chain(
    affine: Mapping[str, Any],
    forcing: Mapping[str, Any],
    join: Mapping[str, Any],
    green: Mapping[str, Any],
    hashes: Mapping[str, str],
) -> tuple[Any, Any, dict[str, Any]]:
    family = forcing.get("source_family")
    logs = forcing.get("source_logs")
    require(family is not None and logs is not None, "actual forcing source family/logs missing")
    require(forcing.get("same_source_radial_and_core_receipts_verified") is True, "actual forcing source-chain verification missing")
    require(forcing.get("allword_tail_added") is False, "actual forcing added all-word tail")
    forcing_metadata = forcing.get("calculation_metadata")
    require(isinstance(forcing_metadata, Mapping), "actual forcing calculation metadata missing")
    require(forcing_metadata.get("same_source_dynamic_C0_and_B1_Z0_to_Z2") is True, "forcing same-source C0/B1 metadata guard")
    require(forcing_metadata.get("F0_relative_derivative_applied_once") is True, "forcing F0 derivative metadata guard")
    require(forcing_metadata.get("g5_formula") == "2*epsilon*(-axial_second(F))", "forcing g5 formula guard")
    require(forcing_metadata.get("g5_anchor_rule") == "F(Z+h)/F0(Z) jet at each current Z; no extra D_F", "forcing g5 derivative guard")
    for label, output in (("affine", affine), ("join", join), ("green", green)):
        require(output.get("source_family") == family, "source family mismatch: " + label)
        require(output.get("source_logs") == logs, "source logs mismatch: " + label)
    input_sha = affine.get("input_receipt_sha256")
    require(isinstance(input_sha, Mapping) and input_sha.get(FORCING) == hashes[FORCING], "affine forcing receipt pin missing")
    actual_inputs = affine.get("actual_source_inputs")
    require(isinstance(actual_inputs, Mapping), "affine actual source inputs missing")
    require(actual_inputs.get("receipt") == FORCING and actual_inputs.get("receipt_sha256") == hashes[FORCING], "affine actual forcing source mismatch")
    expression = affine.get("expression")
    require(isinstance(expression, Mapping), "affine expression missing")
    require(expression.get("schema") == "row5_affine_forcing_v1", "unexpected affine expression schema")
    require(tuple(expression.get("source_power_basis", ())) == POWER_BASIS, "affine power basis mismatch")
    for key in (
        "physical_h5_materialized",
        "h5_is_physical_function",
        "physical_dx_applied",
        "rho_x2_chain_rule_applied",
        "interval_sum_in_build",
        "formal_powers_materialized",
        "axis_same_owner_proved",
        "same_state_owner_proved",
        "source_owner_proved",
    ):
        require(expression.get(key) is False, "affine forbidden/proof flag changed: " + key)
    expression_scope = expression.get("scope")
    require(isinstance(expression_scope, Mapping), "affine scope missing")
    for key in ("state_materialized", "physical_h5", "physical_dx", "ledger_written", "m15", "ns_validation"):
        require(expression_scope.get(key) is False, "affine scope flag changed: " + key)

    join_result = join.get("result")
    green_result = green.get("result")
    require(isinstance(join_result, Mapping) and isinstance(green_result, Mapping), "producer result missing")
    require(join.get("domain") == green.get("domain"), "join/green domain mismatch")
    require(join.get("domain") == "rho=4 exp(hb*s), 0<=s<=sc/4, real Z in [-1,1]", "unexpected full-collar domain")
    require(join.get("sc") == green.get("sc"), "join/green sc mismatch")
    require(join.get("phase") == green.get("phase"), "join/green phase mismatch")
    require(join.get("allword_tail_added") is False and join_result.get("allword_tail_added") is False, "join added all-word tail")
    require(green.get("allword_tail_added") is False and green_result.get("allword_tail_added") is False, "green added all-word tail")
    require(join.get("prior_Z2_green_receipt") == GREEN, "join did not pin green receipt")
    require(join.get("local_receipt") == "experiments/root_st073/n1_actual_collar_low_word_v1_receipt.json.gz", "join local collar receipt changed")
    require(join_result.get("ordinary_Z_orders") == [0, 1], "join ordinary Z orders changed")
    require(join_result.get("source_power_basis") == list(POWER_BASIS), "join power basis mismatch")
    require(join_result.get("source_family") == family, "join inner source family mismatch")
    require(join_result.get("finite_sum_is_unmerged_sector_list") is True, "join finite sectors were merged")
    require(join_result.get("local_added_once") is True and join_result.get("core_history_reset") is False, "join history/local contract changed")
    require(join_result.get("source_callback_identity_proved") is False, "join source callback identity was asserted")
    require(join_result.get("same_source_cancellation_claimed") is False, "join cancellation was asserted")
    require(green_result.get("source_power_basis") == list(POWER_BASIS), "green power basis mismatch")
    require(green_result.get("source_family") == family, "green inner source family mismatch")
    require(green_result.get("ordinary_Z_fixed_phase_and_R") is True, "green ordinary-Z contract missing")
    require(green.get("core_source_index_in_parent_receipt") is True, "green parent source-index chain missing")
    require(green.get("prior_forcing_source_index_in_parent_receipt") is True, "green forcing source-index chain missing")
    require(green_result.get("core_history_reset") is False, "green history was reset")
    require(green_result.get("source_callback_identity_proved") is False, "green source callback identity was asserted")
    require(green_result.get("same_source_cancellation_claimed") is False, "green cancellation was asserted")
    for key in ("error_classes_preserved", "uncertainty_preserved", "source_ids_union_retained", "source_count_retained", "aggregation_is_outward_interval_addition_only"):
        require(green_result.get(key) is True, "green provenance contract changed: " + key)
    local_index = green.get("local_source_index")
    require(isinstance(local_index, Mapping) and local_index, "green local forcing source index missing")
    require(all(isinstance(item, Mapping) and item.get("receipt") == FORCING for item in local_index.values()), "green local source is not forcing-pinned")
    source_pin_checks = {
        "actual_forcing_sha256": hashes[FORCING],
        "affine_actual_source_receipt": actual_inputs.get("receipt"),
        "green_local_source_index_receipts": sorted({item.get("receipt") for item in local_index.values()}),
        "green_prior_forcing_source_index_in_parent_receipt": True,
        "join_prior_Z2_green_receipt": join.get("prior_Z2_green_receipt"),
        "join_local_chain_receipt": join.get("local_receipt"),
        "join_direct_forcing_pin": False,
    }
    return family, logs, source_pin_checks


def consume(hashes: dict[str, str], outcome: dict[str, Any]) -> dict[str, Any]:
    # Authenticate every raw byte before any gzip decompression or new-module
    # import.  In particular, no producer module is imported or called.
    raw: dict[str, bytes] = {}
    for path in REQUIRED:
        data = (ROOT / path).read_bytes()
        require(sha(data) == hashes[path], "changed after entry: " + path)
        raw[path] = data

    from mpmath.ctx_iv import MPIntervalContext

    leaf = module(LEAF, NAME + "_leaf")
    fixture = module(FIXTURE, NAME + "_fixture")
    c = MPIntervalContext()
    c.prec = 1800
    proof = fixture.run_fixture(leaf.build, c)
    require(isinstance(proof, Mapping), "fixture result must be a mapping")
    # Preserve every fixture check in the eventual receipt before fail-closed validation.
    outcome["fixture"] = encode_intervals(proof)
    require(proof.get("schema") == "row5_firstword_product_fixture_v4", "unexpected v4 fixture schema")
    require(proof.get("passed") is True, "fresh synthetic fixture did not pass")
    if isinstance(proof.get("checks"), Mapping):
        require(all(value is True for value in proof["checks"].values()), "fresh synthetic fixture has failed checks")

    affine_receipt = json.loads(gzip.decompress(raw[AFFINE]))
    join_receipt = json.loads(gzip.decompress(raw[JOIN]))
    green_receipt = json.loads(gzip.decompress(raw[GREEN]))
    forcing_receipt = json.loads(gzip.decompress(raw[FORCING]))
    # All failure archives are authenticated and retained only as raw-SHA provenance; none is decompressed or re-entered.
    require(sha(raw[FAILURE_V1]) == PIN[FAILURE_V1], "v1 failure archive changed after entry")
    require(sha(raw[FAILURE_V2]) == PIN[FAILURE_V2], "v2 failure archive changed after entry")
    require(sha(raw[FAILURE_V3_PREFLIGHT]) == PIN[FAILURE_V3_PREFLIGHT], "v3 preflight archive changed after entry")
    green_pins = {item["path"]: item["sha256"] for item in green_receipt["source"]["files"]}
    join_pins = {item["path"]: item["sha256"] for item in join_receipt["source"]["files"]}
    affine_pins = {item["path"]: item["sha256"] for item in affine_receipt["source"]["files"]}
    require(green_pins[FORCING] == hashes[FORCING], "green forcing manifest pin mismatch")
    require(join_pins[GREEN] == hashes[GREEN], "join green manifest pin mismatch")
    require(affine_pins[FORCING] == hashes[FORCING], "affine forcing manifest pin mismatch")
    affine = receipt_output(affine_receipt, AFFINE)
    join = receipt_output(join_receipt, JOIN)
    green = receipt_output(green_receipt, GREEN)
    forcing = receipt_output(forcing_receipt, FORCING)
    family, logs, source_pin_checks = validate_chain(affine, forcing, join, green, hashes)

    expression_raw = clone(affine["expression"])
    expression = decode_intervals(c, expression_raw)
    expression_index = source_index_for_expression(expression_raw)
    state, state_index, state_counts = build_state(c, join, green)

    # This is the sole actual receipt-backed consumer call.  The fixture's
    # synthetic calls above are deliberately reported separately from it.
    outcome["actual_build_calls"] = 1
    result = leaf.build(c, expression, state)
    require(isinstance(result, Mapping), "product result must be a mapping")
    require(result.get("schema") == NAME, "product result schema mismatch")
    require(result.get("status") == "conditional_firstword_substituted_expression", "product status mismatch")
    require(result.get("formal_atoms_materialized") is False, "formal atoms were materialized")
    require(result.get("physical_h5_materialized") is False and result.get("full_h5") is False, "physical h5 was materialized")
    require(result.get("interval_sum_performed") is False, "interval sum was performed")
    require(result.get("unknown_residual", {}).get("symbol") == "R", "unknown residual symbol missing")
    require(result["unknown_residual"].get("preserved") is True, "unknown residual was not preserved")
    require(result["unknown_residual"].get("included_in_products") is False, "unknown residual entered products")
    require(result["unknown_residual"].get("added_to_ledger") is False, "unknown residual entered ledger")
    conditions = result.get("conditions")
    require(isinstance(conditions, Mapping), "product conditions missing")
    require(conditions.get("physical_dx_applied") is False, "physical dx was applied")
    require(conditions.get("same_gauge") == "caller condition; not proved here", "gauge proof was asserted")
    require(conditions.get("live_owner") == "caller condition; not proved here", "live owner proof was asserted")
    require(conditions.get("same_state_owner") == "caller condition; not proved here", "state owner proof was asserted")
    require(result.get("scope", {}).get("no_error_cancellation") is True, "error cancellation scope changed")

    result_counts = {
        row: len(result.get("rows", {}).get(row, []))
        for row in ROWS
    }
    require(all(count > 0 for count in result_counts.values()), "product rows are empty")
    return {
        "source_family": clone(family),
        "source_logs": clone(logs),
        "domain": join["domain"],
        "same_k0_firstword_construction": {
            "source_family_equal": True,
            "source_logs_equal": True,
            "sc_equal": True,
            "phase_equal": True,
            "full_collar_domain_equal": True,
            "domain": join["domain"],
            "join_finite_source": "K0 first-word finite_sum Z0/Z1",
            "green_finite_source": "K0 first-word finite_W1_Z2",
            "prior_Z2_green_receipt": join["prior_Z2_green_receipt"],
            "same_live_owner_proved": False,
            "same_gauge_proved": False,
            "same_physical_h5_proved": False,
        },
        "input_receipt_sha256": {path: hashes[path] for path in (AFFINE, JOIN, GREEN, FORCING)},
        "failure_archives": [
            {"path": FAILURE_V1, "sha256": hashes[FAILURE_V1], "raw_sha256_only": True},
            {"path": FAILURE_V2, "sha256": hashes[FAILURE_V2], "raw_sha256_only": True},
            {"path": FAILURE_V3_PREFLIGHT, "sha256": hashes[FAILURE_V3_PREFLIGHT], "raw_sha256_only": True},
        ],
        "source_pin_checks": source_pin_checks,
        "source_index": expression_index,
        "expression_source_index": expression_index,
        "state_index": state_index,
        "state_counts": state_counts,
        "state_availability": {
            "W1": ["Z0", "Z1", "Z2"],
            "W2": ["Z0", "Z1", "Z2"],
            "W3": ["Z0", "Z1", "Z2"],
            "W5": ["Z0", "Z1"],
            "missing": ["W5_Z2"],
        },
        "expression_counts": {row: len(expression_raw["rows"][row]) for row in ROWS},
        "result_counts": result_counts,
        "counts": {
            "expression": {row: len(expression_raw["rows"][row]) for row in ROWS},
            "state": state_counts,
            "products": result_counts,
        },
        "fixture": encode_intervals(proof),
        "expression": expression_raw,
        "result": encode_intervals(result),
        "scope": {
            "actual_build_calls": 1,
            "synthetic_fixture_build_calls_not_receipt_backed": True,
            "actual_evaluate_called": False,
            "formal_atoms_materialized": False,
            "physical_dx_applied": False,
            "physical_h5_materialized": False,
            "full_state": False,
            "unknown_residual_preserved": True,
            "ledger_written": False,
            "same_live_owner_proved": False,
            "same_gauge_proved": False,
            "same_physical_h5_proved": False,
            "ancestor_called": False,
            "ancestor_imports": False,
            "raw_gzip_authenticated_before_decompression": True,
            "directed_mpf_raw_tuples_retained": True,
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
            "actual_build_calls": 0,
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
        "actual_build_calls": 0,
        "fixture": None,
        "output": None,
    }
    try:
        payload = consume(hashes, outcome)
        outcome.update(
            status="returned",
            output=payload,
            fixture=payload.get("fixture"),
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






