"""Construct/check original inner y4 sources from saved six-row equations."""
import argparse
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/"src"))
sys.setrecursionlimit(max(20000, sys.getrecursionlimit()))
from openai_ns_reconstruction.original_n1_source_dispatcher import OriginalN1SourceDispatcher
from openai_ns_reconstruction.original_n1_inner_high_sources import DIAGONAL, FIELDS

STEM = "lei_ren_part1_paper_compliant_current_original_n1_inner_high_source_functions"
DIRECTORY = ROOT/"experiments"/"root_st073"


class Algebra:
    """Independent scalar interpretation, holding actual matrix entries free."""
    def __init__(self, owner, held=None, analytic=True):
        import sympy as s
        self.s, self.g = s, owner.graph
        self.cache = dict(held or {})
        self.analytic = analytic

    def at(self, i):
        if i in self.cache:
            return self.cache[i]
        s, g = self.s, self.g
        row = g.nodes[i]
        op = row["operation"]
        if op == "exact_rational":
            out = s.Rational(row["numerator"], row["denominator"])
        elif op == "bound_variable":
            out = s.Symbol(row["name"], real=True)
        elif op == "sum":
            out = s.Add(*(self.at(j) for j in row["arguments"]))
        elif op == "product":
            out = s.Mul(*(self.at(j) for j in row["arguments"]))
        elif op == "negative":
            out = -self.at(row["argument"])
        elif op == "positive_quotient":
            out = self.at(row["numerator"])/self.at(row["denominator"])
        elif op == "analytic_unary" and row["name"] in ("exp", "log", "sqrt", "positive_sqrt", "log1p"):
            x = self.at(row["argument"])
            exp = s.exp if self.analytic else s.Function("source_exp")
            out = {"exp": exp, "log": s.log, "sqrt": s.sqrt, "positive_sqrt": s.sqrt,
                   "log1p": lambda x: s.log(1+x)}[row["name"]](x)
        elif op == "exact_original_source_flat_sigma" or op == "analytic_unary" and row["name"].startswith("original_flat_sigma"):
            name = row.get("name", "original_flat_sigma")
            order = {"original_flat_sigma": 0, "original_flat_sigma_prime": 1,
                     "original_flat_sigma_second": 2}.get(name, row.get("derivative_order"))
            out = s.Function("sigma"+str(order))(self.at(row["argument"]))
        elif op == "current_original_unique_regular_n1_component":
            out = s.Symbol("W_"+row["chart"]+"_"+row["quantity"]+"_Z"+str(row["ordinary_slow_Z_derivative_order"]))
        elif op == "regular_even_radial_axis_continuation":
            out = self.at(row["interior_expression"])
        elif op == "function_substitution":
            out = self.at(row["expression"]).subs(self.at(row["variable"]), self.at(row["value"]))
        elif op in ("current_original_source_parameter", "original_source_parameter"):
            out = s.Symbol(row["name"], positive=True)
        else:
            raise ValueError("Unbound independent checker operation: "+op)
        self.cache[i] = out
        return out


def check(owner):
    import sympy as s
    c, g = owner.inner_high_calculus, owner.graph
    retained = 0
    for chart, groups in owner.inner_high_lower_rows.items():
        for group, fields in groups.items():
            for field, jets in fields.items():
                assert owner.charts[chart][group][field][:len(jets)] == jets
                retained += sum(map(len, jets))
    assert retained == 414
    print("INNER_HIGH_LOWER_ROWS_PRESERVED", retained, flush=True)
    rhs_pairs = 0
    lower_second_pairs = 0
    for name, system in c.inner["systems"].items():
        r = s.Symbol("rho", positive=True)
        held = {system["rho"]: r, system["x"]: s.sqrt(r), c.parameters["hb"]: s.Symbol("hb", positive=True)}
        held.update({i:s.Symbol("parameter_"+key,positive=True) for key,i in c.parameters.items() if i not in (0,1)})
        held[system["rho"]], held[system["x"]] = r, s.sqrt(r)
        for label in ("B0", "B1", "g"):
            for key, ids in system[label].items():
                for k, i in enumerate(ids):
                    if i not in (0, 1):
                        held[i] = s.Symbol(label+key.replace(",", "_")+"_Z"+str(k))
        a = Algebra(owner, held, analytic=False)
        state = c.inner["functions"][name]
        W = s.Matrix([a.at(state[field][0]) for field in FIELDS])
        WZ = s.Matrix([a.at(state[field][1]) for field in FIELDS])
        B = {}
        for label in ("B0", "B1"):
            matrix = s.zeros(6)
            for key, ids in system[label].items():
                row, col = map(int, key.split(","))
                matrix[row-1,col-1] = a.at(ids[0])
            B[label] = matrix
        forcing = s.Matrix([a.at(system["g"][str(i+1)][0]) if str(i+1) in system["g"] else 0 for i in range(6)])
        expected = s.sqrt(r)/2*(B["B0"]*W+B["B1"]*WZ+forcing)-s.diag(*DIAGONAL)*W/2
        for i, field in enumerate(FIELDS):
            assert s.cancel(a.at(c.equations[name][field]["physical_rhs"])-expected[i]) == 0
            rhs_pairs += 1
        # Compare an actual differentiated lower row with the retained rho2
        # row, using arbitrary original B/g entries and regular state values.
        chart = "inner_core" if name == "core" else "inner_first"
        for field in ("F1", "Uz1"):
            jets = owner.inner_high_lower_rows[chart]["fields"][field]
            for z in range(4):
                actual = c.dy(chart, jets[1][z])
                assert s.cancel(a.at(actual)-a.at(jets[2][z])) == 0, (chart, field, z)
                lower_second_pairs += 1
        print("INNER_HIGH_ACTUAL_EQUATION_AND_LOW_RADIAL", name, rhs_pairs, lower_second_pairs, flush=True)
    # Independent analytic fixtures exercise full moving-endpoint Leibniz
    # and capture-aware composition, including a bound variable collision.
    x = g.node("bound_variable", name="inner_high_check_x")
    t = g.node("bound_variable", name="inner_high_check_t")
    c.variables["inner_high_check_x"], c.variables["inner_high_check_t"] = x, t
    power = lambda v, n: g.mul(*([v]*n))
    fixtures = [g.mul(g.unary("exp", x), power(x, 3)), g.quotient(g.add(1, power(x, 2)), g.add(1, x)),
        g.node("function_substitution", expression=g.mul(power(t, 3), x), variable=t, value=power(x, 2),
               Z_independent_substitution=True),
        g.node("function_substitution", expression=power(x, 3), variable=x, value=power(x, 2),
               Z_independent_substitution=True)]
    for f in fixtures:
        a = Algebra(owner)
        sx = a.at(x)
        assert s.simplify(a.at(c.partial(f,x))-s.diff(a.at(f),sx)) == 0
    body = g.add(power(t, 2), g.mul(t,x))
    integral = g.node("definite_integral", variable="inner_high_check_t", integrand=body,
                      lower=x, upper=power(x,2))
    derivative = c.partial(integral,x)
    # Evaluate integrals independently rather than copying the source AD rule.
    class IntegralAlgebra(Algebra):
        def at(self, i):
            row = self.g.nodes[i]
            if row["operation"] == "definite_integral":
                return s.integrate(self.at(row["integrand"]),
                    (self.at(t),self.at(row["lower"]),self.at(row["upper"])))
            return super().at(i)
    a = IntegralAlgebra(owner)
    assert s.simplify(a.at(derivative)-s.diff(a.at(integral),a.at(x))) == 0
    # Actual flat-source formulas agree with derivatives of the original
    # logistic function on its open support. Endpoint flatness is inherited.
    a = Algebra(owner)
    sx, sig = a.at(x), s.Symbol("sigma_value")
    class SigmaAlgebra(Algebra):
        def at(self, i):
            row = self.g.nodes[i]
            if row["operation"] == "exact_original_source_flat_sigma" and row["argument"] == x:
                return sig
            return super().at(i)
    a = SigmaAlgebra(owner)
    # Differentiate the original logistic by an independent rational
    # derivation with d(sigma)/dx=sigma*(1-sigma)*d(logit)/dx.
    logit = 1/(1-sx)**2-1/sx**2
    expected = sig
    for order in range(1,5):
        node = c.sigma(x,order)
        lhs = a.at(g.nodes[node]["interior_body"])
        expected = s.diff(expected,sx)+s.diff(expected,sig)*sig*(1-sig)*s.diff(logit,sx)
        assert s.cancel(lhs-expected) == 0, ("sigma",order)
    # All newly exposed finite mixed rows address the same original domains.
    queries = 0
    for chart in ("inner_core", "inner_first", "lower"):
        for group, fields in owner.contract(chart).items():
            for field, contract in fields.items():
                assert contract["physical_y_order"] == (4 if group == "fields" and field != "P1" else 3)
                for y in range(contract["physical_y_order"]+1):
                    for z in range(contract["ordinary_Z_order"]+1):
                        for fraction in ("1/4", "1/2", "1"):
                            result = owner.query(chart, field, group=group, fraction=fraction,
                                                 physical_y_order=y, ordinary_Z_order=z)
                            assert result.source_node == owner.charts[chart][group][field][y][z]
                            assert result.source_node < owner.source_construction_end
                            queries += 1
    for (node, variable), derivative in c.regular_requests.items():
        row = g.nodes[node]
        assert "radial_order" not in row and row["all_words_retained_in_analytic_limit"]
        assert g.nodes[derivative]["operation"] != "current_original_unique_regular_n1_component"
    constructed = g.nodes[g.archived_nodes:owner.source_construction_end]
    regular = [r for r in constructed if r["operation"] == "current_original_unique_regular_n1_component"]
    assert max(r["ordinary_slow_Z_derivative_order"] for r in regular) == 5
    # A changed definition, an invented regular radial tag and a derivative
    # beyond the recorded finite analytic orders must each fail closed.
    template = dict(g.nodes[c.inner["functions"]["core"]["F1"][0]])
    bad = []
    for update in ({"radial_order":1},{"exact_function_definition":"changed"},
                   {"ordinary_slow_Z_derivative_order":5}):
        attrs = dict(template)
        attrs.update(update)
        op = attrs.pop("operation")
        bad.append(g.node(op,**attrs))
    for i in bad:
        try:
            c.partial(i,c.z)
        except ValueError:
            pass
        else:
            raise AssertionError("Unadmitted regular source derivative accepted")
    return dict(passed=True, unchanged_lower_scalar_rows=retained, actual_six_row_equation_identities=rhs_pairs,
        actual_retained_second_derivative_comparisons=lower_second_pairs, independent_partial_Leibniz_fixtures=5,
        original_sigma_interior_derivative_identities=4, inner_mixed_source_queries=queries,
        regular_radial_derivatives_expanded_from_equation=len(c.regular_requests),
        regular_source_max_Z=5, unadmitted_source_requests_rejected=len(bad),
        analytic_core_max_radial=max(g.nodes[i]["radial_order"] for i in c.core_requests),
        analytic_core_max_Z=max(g.nodes[i]["ordinary_slow_Z_derivative_order"] for i in c.core_requests))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--check", action="store_true")
    p.add_argument("--write-report", action="store_true")
    p.add_argument("--restore-source", action="store_true")
    args = p.parse_args()
    if args.write_report and not args.check:
        p.error("--write-report requires --check")
    began = time.monotonic()
    o = OriginalN1SourceDispatcher().extend_inner_high_sources(from_archive=args.restore_source)
    if args.restore_source:
        if args.check or args.write_report:
            p.error("Archive-only recovery is separate from construction/checking")
        queries = 0
        for chart in ("inner_core","inner_first","lower"):
            for group, fields in o.contract(chart).items():
                for field, contract in fields.items():
                    for q in ("1/4","1"):
                        o.query(chart,field,group=group,fraction=q,physical_y_order=contract["physical_y_order"],
                                ordinary_Z_order=contract["ordinary_Z_order"])
                        queries += 1
        print(json.dumps(dict(source_only_restored=True,full_owner_restored=False,
            inner_high_source_restored=o.inner_high_source_restored,
            actual_restored_source_queries=queries,constructed_graph_nodes=o.source_construction_end,
            coverage=o.coverage,execution_seconds=round(time.monotonic()-began,3)),indent=2))
        return
    construction_end = o.source_construction_end
    report = dict(source_family=o.source_family, graph_prefix_nodes=o.graph.archived_nodes,
                  constructed_graph_nodes=construction_end, coverage=o.coverage)
    if args.check:
        report["checks"] = check(o)
    if args.write_report:
        paths = [Path(__file__),ROOT/"src/openai_ns_reconstruction/original_n1_inner_high_sources.py",
                 ROOT/"src/openai_ns_reconstruction/original_n1_source_dispatcher.py"]
        pins = {str(path.relative_to(ROOT)).replace("\\","/"):hashlib.sha256(path.read_bytes()).hexdigest() for path in paths}
        archive = dict(source_family=o.source_family, base=o.source_archives["Rm_outer_mixed_source_joins"],
            graph_prefix_nodes=o.graph.archived_nodes, graph_suffix=o.graph.nodes[o.graph.archived_nodes:construction_end],
            charts={name:o.charts[name] for name in ("inner_core","inner_first","lower")},
            lower_rows=o.inner_high_lower_rows, original_equations=o.inner_high_calculus.equations,
            analytic_source_extension_nodes=[i for i in range(o.graph.archived_nodes,construction_end)
                if o.graph.nodes[i]["operation"] in ("current_original_unique_regular_n1_component",
                                                    "current_original_unique_analytic_core_component")],
            regular_radial_equation_requests=[dict(source=i,variable=v,result=j)
                for (i,v),j in o.inner_high_calculus.regular_requests.items()],
            sigma_interior_source_pin=o.inner_high_calculus.sigma_pin,
            original_input_hashes=o.input_hashes, implementation_hashes=pins, coverage=o.coverage)
        path = DIRECTORY/(STEM+".json.gz")
        path.write_bytes(gzip.compress((json.dumps(archive,separators=(",",":"))+"\n").encode(),mtime=0))
        # Independent serialization recovery, without a new owner or replay.
        with gzip.open(path,"rt",encoding="utf-8") as f:
            restored = json.load(f)
        assert restored == archive
        report.update(implementation_hashes=pins, archive_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                      source_archive_roundtrip=True)
    report["execution_seconds"] = round(time.monotonic()-began,3)
    if args.write_report:
        (DIRECTORY/(STEM+"_check.json")).write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8",newline="\n")
    print(json.dumps(report,indent=2))


if __name__ == "__main__":
    main()
