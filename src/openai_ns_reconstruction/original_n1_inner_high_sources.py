"""Positive-radius inner high sources from the original six-row equation.

No radial derivative label is added to the regular n1 solution. Its partial
derivatives are expanded through its pinned differential equation. Analytic
nonlinear core leaves retain their original admitted component definitions.
"""
import json
import math
from fractions import Fraction
import gzip
import hashlib
from pathlib import Path

FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")
DIAGONAL = (0, 0, 2, 0, 3, 1)


def _definition(row, *orders):
    return json.dumps({k: v for k, v in row.items() if k not in orders}, sort_keys=True, separators=(",", ":"))


class InnerSourceCalculus:
    """Exact partial chain/Leibniz rules and the actual regular n1 equation."""
    def __init__(self, owner):
        self.owner = owner
        self.g = owner.graph
        self.rin = owner._inner_inputs["rin"]
        self.high = owner._inner_inputs["high"]
        self.inner = self.rin["inputs"]["inner_Volterra_bindings"]
        self.parameters = self.inner["parameters"]
        self.variables = {r["name"]: i for i, r in enumerate(self.g.nodes) if r["operation"] == "bound_variable"}
        self.z = self.variables["Z"]
        self.cache = {}
        self.equations = {}
        self.regular_requests = {}
        self.core_requests = set()
        self.sigma_requests = {}
        self.regular_templates = { _definition(self.g.nodes[ids[0]], "ordinary_slow_Z_derivative_order")
            for values in self.inner["functions"].values() for ids in values.values() }
        self.core_templates = { _definition(row, "radial_order", "ordinary_slow_Z_derivative_order")
                                for row in self.high["core_templates"].values() }
        self.sigma_templates = {row["argument"]: dict(row) for row in self.g.nodes
                                if row["operation"] == "exact_original_source_flat_sigma"}
        self.sigma_template = dict(self.g.nodes[self.high["sigma"][0]])
        if self.sigma_template["operation"] != "exact_original_source_flat_sigma":
            raise ValueError("Actual original sigma source template required")
        self.sigma_pin = self.g.nodes[self.high["sigma"][3]]["derivative_source"]
        for chart, system in self.inner["systems"].items():
            v = self.inner["coordinates"][chart]
            r, x = system["rho"], system["x"]
            if chart == "first":
                self.cache[(v, r)] = self.g.mul(self.parameters["hb"], r)
                self.cache[(v, x)] = self.g.mul(self.g.constant("1/2"), self.parameters["hb"], x)
        for group in self.inner["systems"].values():
            for label in ("B0", "B1", "g"):
                for ids in group[label].values():
                    self._seed_Z(ids)
        old = self.rin["inputs"]["core_lower_C5_C3_bindings"]
        for ids in old["core_functions"].values():
            self._seed_Z(ids)
        for rows in old["core_radial_rows"].values():
            for ids in rows:
                self._seed_Z(ids)
        for rows in old["radial_fields"].values():
            for ids in rows:
                self._seed_Z(ids)
        for ids in old["fields"].values():
            self._seed_Z(ids)
        # Each equation is d_y W = x/2 * (B0 W + B1 dZ W + g) - D W/2.
        # The full six-row equation includes its original kinematic entries.
        for chart, system in self.inner["systems"].items():
            state = self.inner["functions"][chart]
            x = system["x"]
            rows = {}
            for i, field in enumerate(FIELDS):
                terms = [system["g"][str(i+1)][0]] if str(i+1) in system["g"] else []
                for label, zorder in (("B0", 0), ("B1", 1)):
                    terms += [self.g.mul(ids[0], state[FIELDS[int(key.split(",")[1])-1]][zorder])
                              for key, ids in system[label].items() if key.split(",")[0] == str(i+1)]
                H = self.g.add(*terms)
                rhs = self.g.add(self.g.mul(self.g.constant("1/2"), x, H),
                                 self.g.mul(self.g.constant(Fraction(-DIAGONAL[i], 2)), state[field][0]))
                rows[field] = dict(H=H, physical_rhs=rhs)
            self.equations[chart] = rows

    def _seed_Z(self, ids):
        for left, right in zip(ids, ids[1:]):
            key = (self.z, left)
            if key in self.cache and self.cache[key] != right:
                raise ValueError("Conflicting accepted ordinary-Z source rows")
            self.cache[key] = right

    def _sub(self, body, variable, value):
        return self.g.node("function_substitution", expression=body, variable=variable, value=value,
                           Z_independent_substitution=self.partial(value, self.z) == 0,
                           exact_partial_chain_rule_substitution=True)

    def sigma(self, argument, order):
        if not 1 <= order <= 4:
            raise ValueError("Only the required actual sigma orders 1..4 are admitted")
        key = (argument, order)
        if key in self.sigma_requests:
            return self.sigma_requests[key]
        g = self.g
        row = dict(self.sigma_templates.get(argument, self.sigma_template))
        row.pop("operation")
        row["argument"] = argument
        sig = g.node("exact_original_source_flat_sigma", **row)
        complement = g.add(1, g.neg(argument))
        power = lambda v, n: g.mul(*([v]*n))
        inv = lambda n, sign: g.add(g.quotient(g.constant(sign), power(argument, n)),
                                  g.quotient(g.constant(abs(sign)), power(complement, n)))
        L1, L2, L3, L4 = inv(3, 2), inv(4, -6), inv(5, 24), inv(6, -120)
        q = g.mul(sig, g.add(1, g.neg(sig)))
        a = g.add(1, g.mul(g.constant(-2), sig))
        if order == 1:
            body = g.mul(q, L1)
        elif order == 2:
            body = g.mul(q, g.add(L2, g.mul(a, L1, L1)))
        elif order == 3:
            body = g.mul(q, g.add(L3, g.mul(g.constant(3), a, L1, L2),
                                 g.mul(g.add(1, g.mul(g.constant(-6), q)), power(L1, 3))))
        else:
            body = g.mul(q, g.add(L4, g.mul(a, g.add(g.mul(g.constant(4), L1, L3), g.mul(g.constant(3), L2, L2))),
                g.mul(g.constant(6), g.add(1, g.mul(g.constant(-6), q)), L1, L1, L2),
                g.mul(a, g.add(1, g.mul(g.constant(-12), q)), power(L1, 4))))
        name = "original_flat_sigma_" + ({1: "prime", 2: "second"}.get(order, "order"+str(order)))
        out = g.node("analytic_unary", name=name, argument=argument, interior_body=body,
            active_interval=[0, 1], flat_value=0, derivative_order=order,
            active_body_not_evaluated_on_flat_branch=True, ordinary_from_Taylor_factor=math.factorial(order),
            derivative_source=self.sigma_pin, exact_function_not_a_numerical_bound=True)
        self.sigma_requests[key] = out
        return out

    def _regular(self, index, variable):
        g = self.g
        row = g.nodes[index]
        if ("radial_order" in row or _definition(row, "ordinary_slow_Z_derivative_order") not in self.regular_templates
                or row["source_family"] != self.owner.source_family):
            raise ValueError("Changed original regular n1 function definition")
        k = row["ordinary_slow_Z_derivative_order"]
        if variable == row["Z_variable"]:
            if k >= 5:
                raise ValueError("Requested regular Z order outside this finite construction")
            attrs = dict(row)
            attrs.pop("operation")
            attrs["ordinary_slow_Z_derivative_order"] = k+1
            return g.node(row["operation"], **attrs)
        if variable != row["coordinate"]:
            return 0
        chart, field = row["chart"], row["quantity"]
        result = self.equations[chart][field]["physical_rhs"]
        for _ in range(k):
            result = self.partial(result, self.z)
        result = (g.quotient(result, self.inner["systems"][chart]["rho"])
                  if chart == "core" else g.mul(self.parameters["hb"], result))
        self.regular_requests[(index, variable)] = result
        return result

    def partial(self, index, variable):
        key = (variable, index)
        if key in self.cache:
            return self.cache[key]
        g = self.g
        row = g.nodes[index]
        op = row["operation"]
        D = lambda i: self.partial(i, variable)
        if op in ("exact_rational", "mathematical_pi", "mathematical_positive_infinity"):
            out = 0
        elif op == "bound_variable":
            out = 1 if index == variable else 0
        elif op in ("current_original_source_parameter", "original_source_parameter", "exact_positive_integer_expression"):
            out = 0
        elif op == "current_original_unique_regular_n1_component":
            out = self._regular(index, variable)
        elif op == "current_original_unique_analytic_core_component":
            if (_definition(row, "radial_order", "ordinary_slow_Z_derivative_order") not in self.core_templates
                    or row["source_family"] != self.owner.source_family):
                raise ValueError("Changed admitted analytic core definition")
            terms = []
            for coordinate, order in (("rho", "radial_order"), ("Z_variable", "ordinary_slow_Z_derivative_order")):
                factor = D(row[coordinate])
                if factor:
                    attrs = dict(row)
                    attrs.pop("operation")
                    attrs[order] += 1
                    if attrs["radial_order"] > 4 or attrs["ordinary_slow_Z_derivative_order"] > 6:
                        raise ValueError("Analytic core request exceeds finite required source orders")
                    source = g.node(op, **attrs)
                    self.core_requests.add(source)
                    terms.append(g.mul(factor, source))
            out = g.add(*terms)
        elif op == "regular_even_radial_axis_continuation":
            # This construction is scoped to rho>0; no new axis limit is claimed.
            out = D(row["interior_expression"])
        elif op == "sum":
            out = g.add(*(D(i) for i in row["arguments"]))
        elif op == "product":
            ids = row["arguments"]
            out = g.add(*(g.mul(D(i), *(j for pos, j in enumerate(ids) if pos != k))
                          for k, i in enumerate(ids) if D(i)))
        elif op == "negative":
            out = g.neg(D(row["argument"]))
        elif op == "positive_quotient":
            n, d = row["numerator"], row["denominator"]
            out = g.quotient(g.add(g.mul(D(n), d), g.neg(g.mul(n, D(d)))), g.mul(d, d))
        elif op in ("analytic_unary", "exact_original_source_flat_sigma"):
            x = row["argument"]
            dx = D(x)
            name = row.get("name", "original_flat_sigma")
            if not dx:
                out = 0
            elif name in ("exp", "expm1"):
                out = g.mul(g.unary("exp", x), dx)
            elif name == "log":
                out = g.quotient(dx, x)
            elif name == "log1p":
                out = g.quotient(dx, g.add(1, x))
            elif name in ("sqrt", "positive_sqrt"):
                out = g.quotient(dx, g.mul(g.constant(2), index))
            elif name.startswith("original_flat_sigma"):
                order = {"original_flat_sigma": 0, "original_flat_sigma_prime": 1,
                         "original_flat_sigma_second": 2}.get(name)
                if order is None:
                    order = int(name.removeprefix("original_flat_sigma_order"))
                out = g.mul(self.sigma(x, order+1), dx)
            else:
                raise ValueError("Unissued analytic inner partial: " + name)
        elif op == "function_substitution":
            body, v, arg = (row[k] for k in ("expression", "variable", "value"))
            direct = 0 if v == variable else D(body)
            darg = D(arg)
            out = g.add(self._sub(direct, v, arg) if direct else 0,
                        g.mul(self._sub(self.partial(body, v), v, arg), darg) if darg else 0)
        elif op == "definite_integral":
            body, lo, hi = (row[k] for k in ("integrand", "lower", "upper"))
            v = self.variables[row["variable"]]
            partial = 0 if v == variable else D(body)
            integral = 0
            if partial:
                attrs = dict(row)
                attrs.pop("operation")
                attrs["integrand"] = partial
                attrs["exact_partial_of_same_integrand"] = True
                integral = g.node(op, **attrs)
            out = g.add(integral, g.mul(self._sub(body, v, hi), D(hi)) if D(hi) else 0,
                        g.neg(g.mul(self._sub(body, v, lo), D(lo))) if D(lo) else 0)
        else:
            raise ValueError("Unissued actual inner derivative operation: " + op)
        self.cache[key] = out
        return out

    def dy(self, chart, node):
        variable = self.inner["coordinates"]["core" if chart == "inner_core" else "first"]
        partial = self.partial(node, variable)
        return (self.g.mul(variable, partial) if chart == "inner_core"
                else self.g.quotient(partial, self.parameters["hb"]))


def install_inner_high_sources(owner):
    """Extend interval functions while preserving every existing lower row."""
    calculus = InnerSourceCalculus(owner)
    before = {}
    for name in ("inner_core", "inner_first", "lower"):
        chart = owner.charts[name]
        before[name] = {group: {key: [list(j) for j in jets] for key, jets in chart[group].items()}
                        for group in ("fields", "cumulative")}
        for group in ("fields", "cumulative"):
            for field, jets in chart[group].items():
                # Work on new lists: archived binding payloads stay unchanged.
                extended = [list(j) for j in jets]
                target = 4 if group == "fields" and field != "P1" else 3
                while len(extended) <= target:
                    extended.append([calculus.dy(name, i) for i in extended[-1]])
                chart[group][field] = extended
        chart["high_source_method"] = "differentiate same six-row Volterra equation and retained signed histories"
    owner.inner_high_calculus = calculus
    owner.inner_high_lower_rows = before
    owner.coverage["positive_radius_inner_high_source_functions"] = True
    owner.coverage["uniform_positive_radius_finite_derivative_slots"] = True
    owner.coverage["inner_high_interface_germs_checked"] = False
    owner.coverage["unresolved"] = ["axis and Cartesian regularity", "higher inner interface germ comparisons",
        "original numerical evaluation/bounds", "physical recovery and genuine n>=2 temporal recursion"]
    owner._validate_rows()
    owner.source_construction_end = len(owner.graph.nodes)
    return calculus


def restore_inner_high_sources(owner):
    """Load the checked small source suffix; no derivative owner is rebuilt."""
    stem = "lei_ren_part1_paper_compliant_current_original_n1_inner_high_source_functions"
    path = owner.directory/(stem+".json.gz")
    receipt = json.loads((owner.directory/(stem+"_check.json")).read_text(encoding="utf-8"))
    if not receipt["checks"]["passed"] or not receipt.get("source_archive_roundtrip"):
        raise ValueError("Checked inner high source receipt required")
    if hashlib.sha256(path.read_bytes()).hexdigest() != receipt["archive_sha256"]:
        raise ValueError("Changed inner high source archive")
    with gzip.open(path,"rt",encoding="utf-8") as stream:
        raw = json.load(stream)
    root = Path(__file__).resolve().parents[2]
    for name, digest in receipt["implementation_hashes"].items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest() != digest:
            raise ValueError("Changed inner high source implementation: "+name)
    if (raw["implementation_hashes"] != receipt["implementation_hashes"]
            or raw["original_input_hashes"] != owner.input_hashes
            or raw["source_family"] != owner.source_family
            or raw["base"] != owner.source_archives["Rm_outer_mixed_source_joins"]
            or raw["graph_prefix_nodes"] != owner.graph.archived_nodes
            or raw["graph_prefix_nodes"]+len(raw["graph_suffix"]) != receipt["constructed_graph_nodes"]):
        raise ValueError("Same accepted original source family and prefix required")
    offset = owner.graph.archived_nodes
    current_suffix = owner.graph.nodes[offset:]
    if current_suffix != raw["graph_suffix"][:len(current_suffix)]:
        raise ValueError("Restore before issuing queries or constructing another source suffix")
    owner.graph.nodes.extend(raw["graph_suffix"][len(current_suffix):])
    for i in range(offset,len(owner.graph.nodes)):
        owner.graph._keys[owner.graph._key(owner.graph.nodes[i])] = i
    for name, chart in raw["charts"].items():
        owner.charts[name] = chart
    owner.inner_high_lower_rows = raw["lower_rows"]
    owner.coverage = raw["coverage"]
    owner.inner_high_source_restored = True
    owner.source_construction_end = len(owner.graph.nodes)
    owner._validate_rows()
    return owner
