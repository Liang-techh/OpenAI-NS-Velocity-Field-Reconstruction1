"""Exact signed pressure correlation after source-owned seam proofs.

The saved 185690-node source is restored from serialized rows, without
executing any of its owners.  Only a proved adjacent trace may replace a
lower wrapper.  Standalone nodes and every old cap remain in the prefix.
"""
from __future__ import annotations

from collections import deque
from copy import deepcopy
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from .original_signed_graph_pairs import GraphPairs


PARENT = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_long_reference_seam_prerequisites.json"
PARENT_SHA256 = "0a91144aca1caaf0dd2b11170b1aee4339f676eb5f213eea236ae8c814cea0c4"
Q4 = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.json.gz"
HELPER = "src/openai_ns_reconstruction/original_signed_graph_pairs.py"
MODULE = "src/openai_ns_reconstruction/original_long_reference_correlated_pressure.py"
PROOF_MODULE = "src/openai_ns_reconstruction/original_long_reference_modified_seams.py"
SEAMS = ("post_power_upper_to_long_lower", "long_upper_to_reference_lower")
BASE_COUNT = 185690
BASE_DIGEST = "8422242e6c2bce6750323cf2aafd2a5a4bce2bf406f033ddd6a18ad61ae33d28"


def require(value: Any, message: str) -> None:
    if not value:
        raise ValueError(message)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=True).encode()).hexdigest()


def expected_pin(parent: Mapping[str, Any], path: str) -> str:
    if path.startswith("experiments/root_st073/"):
        value = parent["input_hashes"][Path(path).name]
    else:
        value = parent["implementation_hashes"][path]
    require(isinstance(value, str), "Missing accepted source pin: " + path)
    return value


def source_pins(root: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    result = {}
    for path in (PARENT, Q4, HELPER, MODULE, PROOF_MODULE):
        actual = sha(root / path)
        expected = (PARENT_SHA256 if path == PARENT else actual
                    if path in (MODULE, PROOF_MODULE) else expected_pin(parent, path))
        require(actual == expected, "Source bytes changed: " + path)
        result[path] = dict(path=path, expected_sha256=expected,
                            actual_sha256=actual, passed=True, source_executed=False)
    return result


def restored_nodes(root: Path, parent: Mapping[str, Any]) -> list[dict[str, Any]]:
    with gzip.open(root / Q4, "rt", encoding="utf-8") as stream:
        archive = json.load(stream)
    nodes = deepcopy(archive["exact_graph_nodes"])
    require(len(nodes) == 183700, "Q4 serialized node count changed")
    pieces = (
        ("actual_micro_macro_correlated_m15", "merged_namespace", 183700, 987),
        ("actual_macro_switch_RQy_seam", "graph_patch", 184687, 16),
        ("actual_macro_first_switch_pressure_source", "graph_patch", 184703, 195),
        ("actual_switch_post_pressure_source", "graph_patch", 184898, 410),
        ("actual_long_reference_pressure_source", "graph_patch", 185308, 382),
    )
    for key, field, start, length in pieces:
        packet = parent[key][field]
        require(packet["base_node_count"] == start and len(nodes) == start,
                "Serialized segment starts changed: " + key)
        suffix = packet["appended_nodes"]
        require(len(suffix) == length, "Serialized segment length changed: " + key)
        nodes.extend(deepcopy(suffix))
    require(len(nodes) == BASE_COUNT and digest(nodes) == BASE_DIGEST,
            "Full accepted pressure graph prefix changed")
    require(parent["actual_long_reference_pressure_source"]["graph_patch"][
                "full_canonical_nodes_sha256"] == BASE_DIGEST,
            "Parent graph digest changed")
    return nodes


def linear_coefficients(nodes: list[dict[str, Any]], root: int) -> dict[int, int]:
    """Collect only sum/negative edges; every other operation is an atom."""
    out: dict[int, int] = {}
    pending = [(root, 1)]
    while pending:
        node, sign = pending.pop()
        row = nodes[node]
        if row["operation"] == "sum":
            pending.extend((child, sign) for child in row["arguments"])
        elif row["operation"] == "negative":
            pending.append((row["argument"], -sign))
        elif row["operation"] == "exact_rational" and int(row["numerator"]) == 0:
            continue
        else:
            out[node] = out.get(node, 0) + sign
            if out[node] == 0:
                del out[node]
    return out


def integer_references(value: Any, limit: int) -> set[int]:
    """Conservative reference set: metadata integers may add extra edges.

    Every graph reference in these pinned archives is an integer.  Following
    additional integer metadata can reject a cancellation, but cannot hide
    a surviving dependency.  Numerical strings are not interpreted.
    """
    if type(value) is int:
        return {value} if 0 <= value < limit else set()
    if isinstance(value, Mapping):
        return set().union(*(integer_references(item, limit) for item in value.values()))
    if isinstance(value, (list, tuple)):
        return set().union(*(integer_references(item, limit) for item in value))
    return set()


def dependency_masks(nodes: list[dict[str, Any]], targets: list[int]) -> dict[int, int]:
    """Propagate target bits once over the reversed conservative graph."""
    parents: list[list[int]] = [[] for _ in nodes]
    for index, row in enumerate(nodes):
        for child in integer_references(row, len(nodes)):
            parents[child].append(index)
    masks: dict[int, int] = {}
    pending: deque[int] = deque()
    for bit, target in enumerate(targets):
        masks[target] = masks.get(target, 0) | (1 << bit)
        pending.append(target)
    while pending:
        child = pending.popleft()
        for parent in parents[child]:
            mask = masks.get(parent, 0) | masks[child]
            if mask != masks.get(parent, 0):
                masks[parent] = mask
                pending.append(parent)
    return masks


def reduce_pair(graph: GraphPairs, pair: list[int], targets: list[int]) -> tuple[list[int], list[dict[str, Any]]]:
    reduced, rows = [], []
    coefficient_rows = []
    for order, root in enumerate(pair):
        coefficients = linear_coefficients(graph.nodes, root)
        cancelled = {str(node): coefficients.get(node, 0) for node in targets}
        require(all(value == 0 for value in cancelled.values()),
                "A proved internal endpoint did not cancel by its signed coefficient")
        terms = []
        for atom, coefficient in sorted(coefficients.items()):
            # Keep integer multiplicity as repeated signed atoms.  No product
            # needs to be opened by the sum/negative verification below.
            terms.extend([atom] * coefficient if coefficient > 0
                         else [graph.neg(atom)] * (-coefficient))
        result = graph.add(*terms)
        require(linear_coefficients(graph.nodes, result) == coefficients,
                "Exact signed-atom rebuilding changed a coefficient")
        reduced.append(result)
        coefficient_rows.append(coefficients)
        rows.append(dict(ordinary_Z_order=order, input_handle=root, reduced_handle=result,
                         signed_atoms={str(k): v for k, v in sorted(coefficients.items())},
                         target_coefficients=cancelled,
                         exact_integer_arithmetic=True))
    masks = dependency_masks(graph.nodes, targets)
    for coefficients, result, receipt in zip(coefficient_rows, reduced, rows):
        surviving = {str(atom): [target for bit, target in enumerate(targets)
                                 if masks.get(atom, 0) & (1 << bit)]
                     for atom in coefficients}
        require(not any(surviving.values()),
                "A surviving pressure atom retains an internal endpoint dependency")
        require(masks.get(result, 0) == 0,
                "Rebuilt pressure source retains an internal endpoint dependency")
        receipt["surviving_target_dependencies"] = surviving
    return reduced, rows


def build(root: str | Path, parent: Mapping[str, Any], proof: Mapping[str, Any]) -> dict[str, Any]:
    root = Path(root)
    require(sha(root / PARENT) == PARENT_SHA256
            and json.loads((root / PARENT).read_bytes()) == parent,
            "Exact accepted prerequisite parent is required")
    pins = source_pins(root, parent)
    require(proof["passed"] is True and proof["parent_report_sha256"] == PARENT_SHA256
            and proof["source_family"] == parent["source_family"]
            and proof["module_sha256"] == pins[PROOF_MODULE]["actual_sha256"]
            and proof["source_executed"] is False,
            "Modified seam proof provenance is missing")
    source = parent["actual_long_reference_pressure_source"]
    old = source["joined_pressure_source"]
    require(old["joined_contracted"] is None
            and old["joined_uncontracted"] == [185688, 185689],
            "Accepted raw joined source changed")
    flags = [proof["seams"][name]["finite_N_Qy_identity_proven"] is True for name in SEAMS]
    common = dict(passed=True, parent_report=PARENT, parent_report_sha256=PARENT_SHA256,
                  source_family=parent["source_family"], source_pins=pins,
                  source_executed=False, accepted_ancestors_replayed=False,
                  ordinary_Z_orders=[0, 1], proved_seams=dict(zip(SEAMS, flags)),
                  module_sha256=sha(root / MODULE))
    if not any(flags):
        return common | dict(graph_patch=None, joined_contracted=None,
            limits=dict(full_modified_seams=False, internal_endpoints_removed=False,
                        full_pressure_numeric=False, global_five_moment_closure=False,
                        standalone_caps_retained=True, all_old_data_retained=True))

    nodes = restored_nodes(root, parent)
    graph = GraphPairs(nodes)
    require(graph.nodes[:BASE_COUNT] == nodes, "Private graph copy changed the prefix")
    previous_uppers = (old["post_power_upper_RQy"], old["long_reshape_upper_RQy"])
    raw_lowers = (old["long_reshape_lower_RQy"], old["reference_lower_RQy"])
    uppers = (old["long_reshape_upper_RQy"], old["reference_upper_RQy"])
    halves = (old["long_reshape_half_S1"], old["reference_half_S1"])
    current = list(old["parent_joined_contracted"])
    targets, substituted, stages = [], [], []
    for index, seam in enumerate(SEAMS):
        if flags[index]:
            require(proof["seams"][seam]["left_RQy_upper"] == previous_uppers[index]
                    and proof["seams"][seam]["right_RQy_lower_wrapper"] == raw_lowers[index],
                    "Proved seam traces do not match the actual pressure graph wrappers")
            require(proof["seams"][seam]["physical_Qy_C1_equal"] is True,
                    "RQy contraction needs an ordinary-C1 physical Qy identity")
            require(proof["seams"][seam]["physical_RQy_C1_equal"] is True,
                    "RQy contraction also needs the same physical radius")
            lower = previous_uppers[index]
            targets.extend(previous_uppers[index])
            substituted.append(dict(seam=seam, old_lower_wrapper=list(raw_lowers[index]),
                                    equal_adjacent_upper=list(lower)))
        else:
            lower = raw_lowers[index]
        current = [graph.add(current[k], halves[index][k], lower[k], graph.neg(uppers[index][k]))
                   for k in range(2)]
        stages.append(list(current))
    before = list(current)
    joined, rows = reduce_pair(graph, current, targets)
    require(all(linear_coefficients(graph.nodes, value).get(uppers[-1][k]) == -1
                for k, value in enumerate(joined)),
            "The external reference RQy trace lost its signed coefficient")
    require(graph.nodes[:BASE_COUNT] == nodes, "A standalone accepted source node was changed")
    suffix = deepcopy(graph.nodes[BASE_COUNT:])
    return common | dict(
        graph_patch=dict(base_node_count=BASE_COUNT, base_canonical_nodes_sha256=BASE_DIGEST,
                         base_prefix_unchanged=True, appended_nodes=suffix,
                         appended_node_count=len(suffix), full_node_count=len(graph.nodes),
                         full_canonical_nodes_sha256=digest(graph.nodes)),
        joined_uncontracted=list(old["joined_uncontracted"]),
        source_substitutions=substituted, joined_before_linearization=before,
        staged_join_handles=stages, joined_contracted=joined,
        exact_linearization=dict(passed=True, targets=targets, rows=rows,
                                 only_sum_negative_expanded=True,
                                 conservative_integer_dependency_check=True),
        outer_reference_RQy_retained=list(old["reference_upper_RQy"]),
        local_S1_added_once=True, cumulative_S1_not_added=True,
        limits=dict(full_modified_seams=all(flags), internal_endpoints_removed=True,
                    full_pressure_numeric=False, global_five_moment_closure=False,
                    standalone_caps_retained=True, all_old_data_retained=True))
