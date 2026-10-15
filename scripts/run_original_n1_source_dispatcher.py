"""Restore/query regional n1 source expressions; never replay ancestor owners."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from openai_ns_reconstruction.original_n1_source_dispatcher import OriginalN1SourceDispatcher

OUTPUT = ROOT / "experiments" / "root_st073" / "lei_ren_part1_paper_compliant_current_original_n1_source_dispatcher"


def check(owner):
    requests = 0
    slots = 0
    for chart, packet in owner.charts.items():
        fractions = ("1/4", "1/2", "1") if packet["lower_open"] else ("0", "1/2", "1")
        parameters = [dict(tail_log_offset=q) if chart == "tail" else dict(fraction=q) for q in fractions]
        for group, fields in owner.contract(chart).items():
            for field, contract in fields.items():
                for y in range(contract["physical_y_order"] + 1):
                    for z in range(contract["ordinary_Z_order"] + 1):
                        slots += 1
                        for parameter in parameters:
                            result = owner.query(chart, field, group=group, physical_y_order=y,
                                                 ordinary_Z_order=z, **parameter)
                            expected = packet[group][field][y][z]
                            assert result.source_node == expected
                            assert owner.graph.nodes[result.bound_node] == dict(
                                operation="function_substitution", expression=expected,
                                variable=packet["coordinate"], value=result.native_query_node,
                                Z_independent_substitution=True)
                            assert owner.graph.nodes[result.radius_node] == dict(
                                operation="function_substitution", expression=packet["radius"],
                                variable=packet["coordinate"], value=result.native_query_node,
                                Z_independent_substitution=True)
                            assert not result.numerical_oracle_available
                            requests += 1
    bad = [
        ("missing", "P1", dict(fraction="1/2")),
        ("inner_core", "P1", dict(fraction="0")),
        ("inner_core", "P1", dict(fraction="1/2", physical_y_order=2)),
        ("inner_first", "F1", dict(fraction="1/2", physical_y_order=3)),
        ("lower", "M1", dict(fraction="1/2", physical_y_order=4)),
        ("O2_axial", "P1", dict(fraction="1/2", physical_y_order=4)),
        ("reserved", "Q1", dict(fraction="1/2", ordinary_Z_order=3)),
        ("tail", "P1", dict(fraction="1/2")),
        ("tail", "P1", dict(tail_log_offset="-1")),
        ("quiet", "P1", dict(fraction="3/2")),
        ("quiet", "P1", dict(fraction="-1/2")),
        ("quiet", "P1", dict(fraction=0.5)),
        ("quiet", "P1", dict(fraction=True)),
        ("quiet", "P1", dict(fraction="1/2", physical_y_order=True)),
        ("quiet", "P1", dict(fraction="1/2", ordinary_Z_order=-1)),
        ("quiet", "not_a_field", dict(fraction="1/2")),
        ("quiet", "P1", dict(fraction="1/2", tail_log_offset="0")),
        ("quiet", "P1", dict(fraction="1/2", group="leading")),
    ]
    for chart, field, kwargs in bad:
        before = len(owner.graph.nodes)
        try:
            owner.query(chart, field, **kwargs)
        except (ValueError, TypeError):
            pass
        else:
            raise AssertionError("Unsupported query accepted")
        assert len(owner.graph.nodes) == before, "Rejected queries must not mutate the graph"
    # Independent operator identities for the only newly constructed core rows.
    # Historical source/NS proofs are not rerun here.
    import sympy as s
    rho, eps = s.symbols("rho eps", positive=True)
    f = s.Function("f")(rho)
    D = lambda v: rho * s.diff(v, rho)
    identities = [D(D(f)) - (rho*s.diff(f, rho) + rho**2*s.diff(f, rho, 2))]
    for exponent, weight in ((s.Rational(1, 2), s.sqrt(2*eps*rho)), (s.Integer(1), eps*rho)):
        identities += [D(weight*f) - weight*(exponent*f + rho*s.diff(f, rho)),
            D(D(weight*f)) - weight*(exponent**2*f + (2*exponent+1)*rho*s.diff(f, rho)
                                      + rho**2*s.diff(f, rho, 2))]
    assert all(s.simplify(value) == 0 for value in identities)
    assert len(owner.charts) == 27 and slots == 6318 and requests == 18954
    return dict(passed=True, regional_charts=27, finite_mixed_scalar_slots=slots,
                source_queries=requests, rejected_unsupported_queries=len(bad),
                independent_core_physical_operator_identities=len(identities))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--write-report", action="store_true")
    parser.add_argument("--query", nargs=2, metavar=("CHART", "FIELD"))
    parser.add_argument("--fraction")
    parser.add_argument("--tail-log-offset")
    parser.add_argument("--group", choices=("fields", "cumulative"), default="fields")
    parser.add_argument("--dy", type=int, default=0)
    parser.add_argument("--dZ", type=int, default=0)
    args = parser.parse_args()
    if args.write_report and not args.check:
        parser.error("--write-report requires --check")
    began = time.monotonic()
    owner = OriginalN1SourceDispatcher()
    result = dict(source_only_restored=True, full_owner_restored=False,
                  archived_graph_nodes=owner.graph.archived_nodes,
                  regional_charts=list(owner.charts), coverage=owner.coverage)
    if args.check:
        result["checks"] = check(owner)
    if args.query:
        result["query"] = owner.query(*args.query, fraction=args.fraction,
            tail_log_offset=args.tail_log_offset, group=args.group,
            physical_y_order=args.dy, ordinary_Z_order=args.dZ).as_dict()
        result["query_graph_suffix"] = owner.graph.nodes[owner.graph.archived_nodes:]
    result["execution_seconds"] = round(time.monotonic()-began, 3)
    if args.write_report:
        manifest = owner.manifest()
        manifest_path = OUTPUT.with_suffix(".json")
        manifest_path.write_text(json.dumps(manifest, indent=2)+"\n", encoding="utf-8", newline="\n")
        paths = [Path(__file__), ROOT/"src/openai_ns_reconstruction/original_n1_source_dispatcher.py", manifest_path]
        result["implementation_hashes"] = {str(p.relative_to(ROOT)).replace("\\", "/"):
            hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
        OUTPUT.with_name(OUTPUT.name+"_check.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8", newline="\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
