"""Build and inspect the source-bound reference physical-germ packet.

``--write-report`` constructs only :class:`OriginalReferenceGerm` and the
independent source fixture.  It does not evaluate an axial production query
or replay an accepted ancestor.  ``--inspect-report`` verifies the saved
packet, its file hashes, and its serialized graph invariants without creating
the owner.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
SCRIPTS = ROOT / "scripts"
DIRECTORY = ROOT / "experiments" / "root_st073"
for path in (str(ROOT), str(SRC), str(SCRIPTS)):
    if path not in sys.path:
        sys.path.insert(0, path)

from openai_ns_reconstruction.original_core_numeric import PREFIX  # noqa: E402
from openai_ns_reconstruction.original_reference_viscous_m15 import (  # noqa: E402
    ACTIVE,
    ACTIVE_SHA256,
    OriginalReferenceGerm,
)


NAME = PREFIX + "compliant_current_original_reference_physical_germ.json"
FIXTURE_NAME = "scripts/check_original_reference_viscous_m15.py"
RUNNER_NAME = "scripts/run_original_reference_germ.py"
MISSING_Z = "0.1"
SOURCE_PORT_COORDINATE = 81520
SOURCE_PORT_ENDPOINT = 1


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def _require(condition, message):
    if not condition:
        raise AssertionError(message)


def _check_hash_group(base: Path, group: dict, label: str) -> None:
    _require(isinstance(group, dict) and group, f"{label} hash map missing")
    for name, digest in group.items():
        _require(isinstance(name, str) and isinstance(digest, str) and len(digest) == 64,
                 (label, name, "malformed digest"))
        target = base / name
        _require(target.is_file(), (label, name, "missing pinned file"))
        _require(sha(target) == digest, (label, name, "hash changed"))


def _expected_ledger(parent: dict) -> list[dict]:
    samples = parent.get("samples")
    _require(isinstance(samples, list) and len(samples) == 3, "accepted parent samples changed")

    def digest(value) -> str:
        return hashlib.sha256(canonical(value).encode("utf-8")).hexdigest()

    result = []
    for sample in samples:
        raw = sample.get("raw_m15_C0_C1_terms")
        reduced = sample.get("correlated_reduced_m15_C0_C1_terms")
        target = sample.get("independent_target")
        _require(isinstance(raw, list) and len(raw) == 213, (sample.get("Z"), "parent raw count"))
        _require(isinstance(reduced, list) and len(reduced) == 199,
                 (sample.get("Z"), "parent reduced count"))
        _require(target is not None, (sample.get("Z"), "parent independent target"))
        result.append({
            "Z": sample.get("Z"),
            "raw_terms": 213,
            "reduced_terms": 199,
            "raw_ledger_sha256": digest(raw),
            "reduced_ledger_sha256": digest(reduced),
            "independent_target_sha256": digest(target),
        })
    return result


def _check_common_rows(rows) -> dict:
    _require(isinstance(rows, list) and len(rows) == 94,
             ("common source identity rows", None if rows is None else len(rows)))
    allowed_orders = {0, 1, 2}
    expected=set()
    for group,keys,orders in (
        ('geometry',('radius','global_offset','actual_phase'),(0,)),
        ('fields',('E','V','m','h','k','e','p','P0','pressure'),(0,1,2)),
        ('radial',('E_y','V_y'),(0,1,2)),
        ('actual',('phi','U','phi_y','U_y','physical_Utheta','log_Utheta'),(0,1,2)),
        ('root_values',('E','V','a','b','p1','p2','t0','Delta','kappa'),(0,1)),
        ('loop_values',('B_over_Pstar','B_over_Pstar_phi'),(0,1)),
        ('retained_history_inlet',('H','M','K','A','B','C'),(0,1,2))):
        expected.update((group,key,k) for key in keys for k in orders)
    seen=[]
    for row in rows:
        _require(isinstance(row, dict), ("common source row", row))
        required = {"group", "quantity", "ordinary_Z_order", "left_node", "right_node",
                    "canonical_endpoint_sha256"}
        _require(set(row) == required, ("common source row shape", row))
        _require(row["ordinary_Z_order"] in allowed_orders, ("ordinary Z order", row))
        _require(isinstance(row["left_node"], int) and row["left_node"] >= 0 and
                 isinstance(row["right_node"], int) and row["right_node"] >= 0,
                 ("source endpoint IDs", row))
        _require(isinstance(row["canonical_endpoint_sha256"], str) and
                 len(row["canonical_endpoint_sha256"]) == 64, ("source digest", row))
        seen.append((row['group'],row['quantity'],row['ordinary_Z_order']))
    _require(set(seen)==expected and len(set(seen))==len(seen),'complete unique source input grids required')
    return {"rows": len(rows), "ordinary_Z_orders": sorted(allowed_orders), "source_endpoint_ids_actual": True}


def _check_radial_laws(laws) -> dict:
    _require(isinstance(laws, dict) and laws.get("passed") is True, "radial source law failed")
    rows = laws.get("rows")
    _require(isinstance(rows, list) and len(rows) == 6, "whole-chart radial row count changed")
    seen = set()
    for row in rows:
        _require(isinstance(row, dict), ("radial row", row))
        key = (row.get("chart"), row.get("ordinary_Z_order"))
        seen.add(key)
        _require(row.get("chart") in {"reference", "restoration"}, ("radial chart", row))
        _require(row.get("ordinary_Z_order") in {0, 1, 2}, ("radial Z order", row))
        _require(isinstance(row.get("E_node"), int) and isinstance(row.get("E_y_node"), int),
                 ("radial source IDs", row))
        _require(row.get("endpoint_sampling_used") is False and
                 row.get("whole_native_chart_function_identity") is True,
                 ("radial whole-chart proof", row))
        _require(isinstance(row.get("exact_zero_sha256"), str) and
                 len(row["exact_zero_sha256"]) == 64, ("radial exact zero", row))
    _require(seen == {(chart, order) for chart in ("reference", "restoration") for order in range(3)},
             ("radial row coverage", sorted(seen)))
    _require(laws.get("source_values_selected") is False, "radial source values selected")
    _require(laws.get("physical_y2_E_law_by_differentiation") is True,
             "E second physical-y law missing")
    return {"passed": True, "rows": len(rows), "whole_chart_E_y_equals_E_over_10": True}


def _local_children(row: dict) -> list[int]:
    operation = row.get("operation")
    if operation in {"sum", "product"}:
        return row.get("arguments", [])
    if operation == "positive_quotient":
        return [row.get("numerator"), row.get("denominator")]
    return []


def _check_recipe_graph(recipes) -> dict:
    _require(isinstance(recipes, dict), "physical recipe graph missing")
    nodes = recipes.get("recipe_nodes")
    _require(isinstance(nodes, list) and len(nodes) == 489,
             ("recipe node count", None if nodes is None else len(nodes)))
    for index, row in enumerate(nodes):
        _require(isinstance(row, dict) and isinstance(row.get("operation"), str),
                 ("recipe node", index))
        for child in _local_children(row):
            _require(isinstance(child, int) and 0 <= child < index,
                     ("recipe graph edge is not backward", index, child))
        if row["operation"] == "original_endpoint_function":
            # This is an original source-graph port, not a recipe-local edge.
            _require(row.get("original_coordinate") == SOURCE_PORT_COORDINATE and
                     row.get("original_endpoint") == SOURCE_PORT_ENDPOINT,
                     ("source port scope", index, row))
            _require(row.get("ordinary_Z_order") in {0, 1, 2} and
                     isinstance(row.get("original_node"), int) and row["original_node"] >= 0,
                     ("source port ID/order", index, row))
            _require(row.get("function_not_numeric_atom") is True,
                     ("source port treated as numeric", index))
    _require(recipes.get("ordinary_Z_orders") == [0, 1, 2] and
             recipes.get("physical_y_orders") == [0, 1, 2], "recipe jet orders changed")
    _require(recipes.get("source_functions_not_evaluated") is True and
             recipes.get("no_error_caps_used") is True and
             recipes.get("physical_jet_only") is True,
             "recipe source/cap contract changed")

    common = recipes.get("common_fields")
    fields = {"E", "V", "m", "h", "k", "e", "p", "P0", "pressure"}
    _require(isinstance(common, dict) and set(common) == {"y0", "y1", "y2"},
             "common physical field jets missing")
    for level in ("y0", "y1", "y2"):
        table = common[level]
        _require(isinstance(table, dict) and set(table) == fields, ("common field names", level))
        for name, jet in table.items():
            _require(isinstance(jet, list) and len(jet) == 3, ("common field jet", level, name))
            _require(all(isinstance(node, int) and 0 <= node < len(nodes) for node in jet),
                     ("common field node bounds", level, name))

    root_p2 = recipes.get("root_p2")
    _require(isinstance(root_p2, dict), "root p2 recipe missing")
    for key in ("value_Z0_Z1", "physical_y_Z0_Z1", "Iz_Z0_Z1", "Iz_y_Z0_Z1"):
        jet = root_p2.get(key)
        _require(isinstance(jet, list) and len(jet) == 2 and
                 all(isinstance(node, int) and 0 <= node < len(nodes) for node in jet),
                 ("root p2 node bounds", key))
    _require(root_p2.get("identity") == "p2_y=(9/10)*p2+(R/E)*Iz_y", "root p2 identity changed")
    ports = [row for row in nodes if row.get("operation") == "original_endpoint_function"]
    _require(ports and {row["ordinary_Z_order"] for row in ports} == {0, 1, 2},
             "source port ordinary-Z coverage missing")
    return {
        "recipe_nodes": len(nodes),
        "source_port_coordinate": SOURCE_PORT_COORDINATE,
        "source_port_endpoint": SOURCE_PORT_ENDPOINT,
        "source_port_orders": [0, 1, 2],
        "bounded_graph_backward_indices_only": True,
        "source_ports_are_definitions": True,
        "common_fields_y0_y1_y2": True,
        "root_p2_identity": True,
    }


def _check_guards(guards) -> dict:
    _require(isinstance(guards, dict) and guards.get("passed") is True, "source guards failed")
    _require(guards.get("global_periodic_lift_proved") is False,
             "global phase lift was overclaimed")
    checks = guards.get("checks")
    _require(isinstance(checks, dict) and isinstance(checks.get("global_periodic_lift"), dict),
             "source guard diagnostic missing")
    _require(checks["global_periodic_lift"].get("passed") is False,
             "global periodic lift diagnostic changed")
    for name, check in checks.items():
        _require(isinstance(check, dict), ("source guard shape", name))
        if name != "global_periodic_lift":
            _require(check.get("passed") is True, ("source guard failed", name))
    return {"passed": True, "global_periodic_lift_proved": False,
            "other_source_premises_passed": True}


def _check_certificate(certificate) -> dict:
    _require(isinstance(certificate, dict) and certificate.get("passed") is True,
             "physical input two-jet certificate missing")
    _require(certificate.get("direct_endpoint_fingerprint_for_modified_Uy") is False,
             "modified Uy direct fingerprint was overclaimed")
    _require(certificate.get("universal_loop_native_derivative_naturality_not_yet_certified") is True,
             "universal native derivative was overclaimed")
    _require(certificate.get("source_values_selected") is False,
             "physical certificate selected source values")
    return {"passed": True, "direct_modified_Uy_fingerprint": False,
            "universal_native_derivative_unproved": True}


def _check_contract(contract) -> dict:
    required = {
        "accepted_ancestors_not_reexecuted": True,
        "raw_213_and_reduced_199_unchanged": True,
        "all_old_error_budgets_and_independent_target_unchanged": True,
        "derived_recipe_nodes_are_not_claimed_as_original_graph_nodes": True,
        "source_function_definitions_not_numeric_field_evaluation": True,
        "no_new_trace_alias": True,
        "no_error_cap_subtraction": True,
    }
    _require(isinstance(contract, dict), "reference contract missing")
    for name, expected in required.items():
        _require(contract.get(name) is expected, ("contract flag", name, contract.get(name)))
    return {"raw_213_reduced_199_unchanged": True, "no_new_trace_alias": True,
            "no_error_cap_subtraction": True, "source_definitions_only": True}


def _check_limits(limits) -> dict:
    _require(isinstance(limits, dict), "limits missing")
    _require(limits.get("physical_input_two_jet_functions_implemented") is True,
             "physical input two-jet implementation missing")
    for name, value in limits.items():
        if name != "physical_input_two_jet_functions_implemented":
            _require(value is False, ("open limit overclaimed", name, value))
    return {"physical_input_two_jet_functions_implemented": True,
            "all_other_limits_open": True}


def check_packet(packet: dict, parent: dict) -> dict:
    _require(isinstance(packet, dict), "reference packet missing")
    _require(packet.get("source_family") == parent.get("source_family"), "source family changed")
    _require(packet.get("parent_report") == ACTIVE and
             packet.get("parent_report_sha256") == ACTIVE_SHA256,
             "accepted parent identity changed")
    common = _check_common_rows(packet.get("common_source_identity_rows"))
    radial = _check_radial_laws(packet.get("leading_radial_laws"))
    guards = _check_guards(packet.get("source_guards"))
    recipes = _check_recipe_graph(packet.get("physical_germ_recipes"))
    graph=packet['physical_germ_recipes']['recipe_nodes']
    source_fields={(r['quantity'],r['ordinary_Z_order']):r['left_node'] for r in packet['common_source_identity_rows'] if r['group']=='fields'}
    for name,jet in packet['physical_germ_recipes']['common_fields']['y0'].items():
        for k,i in enumerate(jet):
            _require(graph[i].get('operation')=='original_endpoint_function' and
                     graph[i].get('original_node')==source_fields[name,k],('wrong original field port',name,k))
    certificate = _check_certificate(packet.get("physical_two_jet_certificate"))
    contract = _check_contract(packet.get("contract"))
    limits = _check_limits(packet.get("limits"))
    ledger = packet.get("unchanged_m15_ledger")
    expected = _expected_ledger(parent)
    _require(ledger == expected, "raw/reduced/target ledger changed")
    return {
        "common_source_identity_rows": common,
        "leading_radial_laws": radial,
        "source_guards": guards,
        "physical_germ_recipes": recipes,
        "physical_two_jet_certificate": certificate,
        "contract": contract,
        "limits": limits,
        "unchanged_m15_ledger": {"samples": len(ledger), "raw_terms": 213,
                                  "reduced_terms": 199, "independent_target_unchanged": True},
    }


def _check_fixture(fixture: dict) -> dict:
    _require(isinstance(fixture, dict) and fixture.get("passed") is True,
             "independent reference fixture failed")
    sections = ("symbolic", "quadrature", "pressure_jets", "poisson_loop", "recipe_numeric_fixture")
    _require(all(isinstance(fixture.get(name), dict) for name in sections),
             "independent fixture sections missing")
    return {"passed": True, "sections": list(sections), "source_ports_are_synthetic": True}


def _check_missing_coordinate(owner) -> dict:
    try:owner.evaluate(MISSING_Z)
    except ValueError:return {"missing_Z":MISSING_Z,"rejected":1,"reason":"source-only owner explicitly rejects numerical evaluation"}
    raise AssertionError('Source germ was silently evaluated as original axial data')


def run() -> dict:
    started = time.perf_counter()
    from check_original_reference_viscous_m15 import run as run_fixture

    fixture = run_fixture()
    fixture_check = _check_fixture(fixture)
    owner = OriginalReferenceGerm(precision=100)
    packet = owner.build()
    parent = owner.parent
    packet_check = check_packet(packet, parent)
    missing = _check_missing_coordinate(owner)

    owner.code_hashes[FIXTURE_NAME] = sha(ROOT / FIXTURE_NAME)
    owner.code_hashes[RUNNER_NAME] = sha(ROOT / RUNNER_NAME)
    code = dict(owner.code_hashes)
    hashes = dict(owner.hashes)
    checks = {
        "passed": True,
        "independent_fixture": fixture,
        "packet_invariants": packet_check,
        "fixture_invariants": fixture_check,
        "source_port_coordinate": SOURCE_PORT_COORDINATE,
        "source_port_endpoint": SOURCE_PORT_ENDPOINT,
        "recipe_graph_backward_indices_only": True,
        "raw_213_reduced_199_ledger_unchanged": True,
        "all_old_budgets_unchanged": True,
        "independent_target_unchanged": True,
        "direct_modified_Uy_fingerprint": False,
        "universal_native_derivative_unproved": True,
        "global_periodic_lift_proved": False,
        "no_inlet_B_y_or_RQy_alias": True,
        "missing_Z_rejected": missing["rejected"],
        "production_is_new_source_packet_only": True,
    }
    return dict(packet, input_hashes=hashes, implementation_hashes=code,
                checks=checks, execution_seconds=round(time.perf_counter() - started, 3))


def inspect(path: Path) -> dict:
    """Verify a saved packet using hashes and serialized data only."""
    data = json.loads(path.read_text(encoding="utf-8"))
    _require(data.get("checks", {}).get("passed") is True, "saved reference packet is not passed")
    _check_hash_group(DIRECTORY, data.get("input_hashes"), "input")
    _check_hash_group(ROOT, data.get("implementation_hashes"), "implementation")
    parent_path = DIRECTORY / ACTIVE
    parent_raw = parent_path.read_bytes()
    _require(hashlib.sha256(parent_raw).hexdigest() == ACTIVE_SHA256, "accepted active parent changed")
    parent = json.loads(parent_raw)
    packet_check = check_packet(data, parent)
    _require(data.get("checks", {}).get("source_port_coordinate") == SOURCE_PORT_COORDINATE,
             "serialized source coordinate changed")
    _require(data.get("checks", {}).get("source_port_endpoint") == SOURCE_PORT_ENDPOINT,
             "serialized source endpoint changed")
    _require(data["checks"].get("recipe_graph_backward_indices_only") is True,
             "serialized graph bound missing")
    _require(data["checks"].get("raw_213_reduced_199_ledger_unchanged") is True and
             data["checks"].get("all_old_budgets_unchanged") is True,
             "serialized old-budget ledger claim missing")
    _require(data["checks"].get("independent_target_unchanged") is True,
             "serialized target preservation missing")
    _require(data["checks"].get("global_periodic_lift_proved") is False and
             data["checks"].get("direct_modified_Uy_fingerprint") is False and
             data["checks"].get("universal_native_derivative_unproved") is True,
             "unresolved phase/native derivative premise was overclaimed")
    _require(data["checks"].get("no_inlet_B_y_or_RQy_alias") is True,
             "unproved inlet alias was serialized")
    _require(data["checks"].get("missing_Z_rejected") == 1, "unsupported-Z evidence changed")
    fixture_check = _check_fixture(data["checks"].get("independent_fixture"))
    return dict(packet_invariants=packet_check, fixture_invariants=fixture_check,
                checks=data["checks"], execution_seconds=data.get("execution_seconds"))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-report", action="store_true")
    parser.add_argument("--inspect-report", action="store_true")
    args = parser.parse_args()
    path = DIRECTORY / NAME
    if args.inspect_report:
        data = inspect(path)
    else:
        data = run()
    if args.write_report:
        temporary = path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8", newline="\n")
        temporary.replace(path)
    print(json.dumps({"checks": data.get("checks", {}),
                      "execution_seconds": data.get("execution_seconds")}, indent=2), flush=True)
