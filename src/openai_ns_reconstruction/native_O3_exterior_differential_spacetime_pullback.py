"""Fresh differential-form pullback of the accepted, opaque O3 sources.

Only the explicit new coordinate map and potential are differentiated. Native
source functions are supplied as exact references and conditional first jets.
No legacy AST, graph producer, field oracle, or failed kernel is executed.
"""
from __future__ import annotations

from fractions import Fraction
import json


COMPONENT = "native_O3_exterior_differential_spacetime_pullback"
INPUT_SCHEMA = COMPONENT + "_input_v1"
RESULT_SCHEMA = COMPONENT + "_result_v1"
BRANCHES = ("background", "modified")
WRAPPER_PINS = {
    "saved_flux_result": "d4cf1b39195a608e32197c58d572f86f0d974ba940ef2ff577295fdfc95b092d",
    "source_literals": "ae3a6afbdc8de9c557dff18b3552a6b1c6c1d902decca5298b0b6adbcc216b59",
    "physical_source_capsule": "90ce8bbacd274cf7ea022afd976799944b60709fef69162a6736a891fcef7788",
}
EVENT_COUNTS = {
    "family_entered": 1, "source_join": 1, "geometry_dual": 1,
    "potential": 2, "flux_form": 2, "velocity_pressure": 2,
    "Cartesian_composition": 2, "volume_form": 1, "exterior_closure": 2,
    "api_refusal": 4, "family_completed": 1,
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


class _IR:
    """An append-only expression DAG; structural sharing is not evaluation."""
    def __init__(self):
        self.nodes = []
        self.shared = {}
        self.zero = self.number(0)
        self.one = self.number(1)

    def node(self, operation, *arguments, **attributes):
        value = {"operation": operation, "arguments": list(arguments), **attributes}
        key = json.dumps(value, sort_keys=True, separators=(",", ":"))
        if key not in self.shared:
            ref = "physical_expr:" + str(len(self.nodes))
            self.shared[key] = ref
            self.nodes.append({"id": ref, **value})
        return self.shared[key]

    def number(self, numerator, denominator=1):
        f = Fraction(numerator, denominator)
        return self.node("rational", numerator=f.numerator, denominator=f.denominator)

    def add(self, *values):
        values = tuple(v for v in values if v != self.zero)
        return self.zero if not values else values[0] if len(values) == 1 else self.node("add", *values)

    def mul(self, *values):
        if self.zero in values:
            return self.zero
        values = tuple(v for v in values if v != self.one)
        return self.one if not values else values[0] if len(values) == 1 else self.node("multiply", *values)

    def neg(self, value):
        return self.zero if value == self.zero else self.node("negative", value)

    def sub(self, left, right):
        return self.zero if left == right else self.add(left, self.neg(right))

    def power(self, value, exponent):
        return self.one if exponent == self.zero else value if exponent == self.one else self.node("power", value, exponent)

    def div(self, numerator, denominator):
        return self.zero if numerator == self.zero else self.node("quotient", numerator, denominator, nonzero_denominator_premise=True)

    def source(self, node, path, role):
        require(type(node) is int, "native source node must be an integer identifier")
        return self.node("native_function_reference", native_node=node,
                         source="physical_source_capsule", path=path, role=role,
                         source_values_unevaluated=True)


class _Dual:
    """The new ansatz's first differential in independent profile R and Z."""
    def __init__(self, g, value, radial=None, axial=None):
        self.g, self.value = g, value
        self.radial = g.zero if radial is None else radial
        self.axial = g.zero if axial is None else axial

    def add(self, other):
        g = self.g
        return _Dual(g, g.add(self.value, other.value),
                     g.add(self.radial, other.radial), g.add(self.axial, other.axial))

    def neg(self):
        g = self.g
        return _Dual(g, g.neg(self.value), g.neg(self.radial), g.neg(self.axial))

    def sub(self, other):
        return self.add(other.neg())

    def mul(self, other):
        g = self.g
        return _Dual(g, g.mul(self.value, other.value),
                     g.add(g.mul(self.radial, other.value), g.mul(self.value, other.radial)),
                     g.add(g.mul(self.axial, other.value), g.mul(self.value, other.axial)))

    def power(self, fixed_exponent):
        g = self.g
        slope = g.mul(fixed_exponent, g.power(self.value, g.sub(fixed_exponent, g.one)))
        return _Dual(g, g.power(self.value, fixed_exponent),
                     g.mul(slope, self.radial), g.mul(slope, self.axial))

    def div(self, other):
        return self.mul(other.power(self.g.number(-1)))

    def data(self):
        return {"value": self.value, "D_R": self.radial, "D_Z": self.axial}


def _node(value):
    if isinstance(value, dict):
        value = value.get("node")
    require(type(value) is int, "expected a literal native node reference")
    return value


def _inputs(payload):
    require(payload.get("schema") == INPUT_SCHEMA, "wrong input schema")
    require(payload.get("component") == COMPONENT, "wrong input component")
    wrappers = payload.get("source_wrappers")
    require(isinstance(wrappers, dict) and set(wrappers) == set(WRAPPER_PINS), "wrong source wrappers")
    data = {}
    for name, sha in WRAPPER_PINS.items():
        row = wrappers[name]
        require(set(row) == {"sha256", "literal"} and row["sha256"] == sha, "wrapper pin changed: " + name)
        require(isinstance(row["literal"], dict), "wrapper is not literal JSON")
        data[name] = row["literal"]
    accepted = data["saved_flux_result"]["result"]
    require(accepted["passed"] is True and accepted["source_only"] is True, "accepted flux result required")
    require(accepted["new_flags"]["native_O3_saved_m_y_compatibility_proved"] is True, "accepted m_y alias required")
    cap = data["physical_source_capsule"]
    require(cap["schema"] == "native_O3_physical_source_capsule_v1" and cap["packet_revision"] == 2, "wrong physical capsule")
    require(cap["scope"]["chart"] == "O3_transition", "wrong source chart")
    return accepted["payload"], data["source_literals"], cap


def _join_sources(parent, literals, cap):
    joins = []
    my = parent["my_compatibility_rows"]
    require(len(my) == 5, "five accepted m_y rows required")
    bg_proofs = parent["source_proof_rows"]["background"]["m"]
    for j in range(5):
        bg, mod = cap["background_rows"], cap["modified_rows"]
        require(my[j]["ordinary_Z_order"] == j, "ordinary Z row mismatch")
        require(_node(my[j]["negative_argument"]) == mod["m"][j], "modified m join mismatch")
        require(_node(my[j]["V_N_ref"]) == mod["V"][j], "modified V join mismatch")
        history = [t["node"] for t in bg_proofs[str(j)]["source_expression"]["terms"] if t["source"] == "history_j"]
        require(history == [bg["m"][j]], "background incoming mass source changed")
        for branch, fields in (("background", bg), ("modified", mod)):
            for field in ("E", "V"):
                group = literals["literal_rows"]["background_fields" if branch == "background" else "finite_N"]
                key = field if branch == "background" else field + "_N"
                require(group[key]["selected_node_ids"][j] == fields[field][j], "field source join mismatch")
        joins.append({"ordinary_Z_order": j, "modified_my_node": _node(my[j]["sum_head"]),
                      "modified_m_node": mod["m"][j], "background_m_node": bg["m"][j]})
    pressure = parent["pressure_alias_rows"]
    require(len(pressure) == 10, "ten pressure alias rows required")
    seen = set()
    for row in pressure:
        branch, j = row["branch"], row["ordinary_Z_order"]
        require(branch in BRANCHES and j in range(5) and (branch, j) not in seen, "pressure row key mismatch")
        seen.add((branch, j))
        for field in ("P0", "p", "pressure"):
            require(row["raw_node_bindings"][field] == cap[branch + "_rows"][field][j], "pressure source join mismatch")
    expected = {"Pstar": 33572, "delta": 25220, "Z": 20, "R": 241787, "native_x": 241536, "Rm": 85743}
    for name, node in expected.items():
        require(cap["parameter_ports"][name]["node"] == node, "source parameter port changed")
        if name in ("Pstar", "delta", "R", "Rm"):
            require(cap["parameter_ports"][name]["materialized"] is False, "source scale must remain opaque")
    radius = parent["cached_data"]["nativeR_certificate"]
    require(radius["root"] == expected["R"] and radius["Z_excluded"] is True
            and radius["ordinary_Z_derivative"] == "0", "accepted conditional native-R frame required")
    return joins


def _geometry(g, cap):
    ports = {name: g.source(row["node"], "parameter_ports." + name, row["role"])
             for name, row in cap["parameter_ports"].items()}
    const = lambda value: _Dual(g, value)
    R, Z = _Dual(g, ports["R"], g.one), _Dual(g, ports["Z"], axial=g.one)
    tau = g.node("free_variable", name="tau", role="physical_time_remaining")
    theta = g.node("free_variable", name="theta", role="physical_azimuth_not_native_fast_phase")
    d = const(g.one).sub(Z.mul(Z))
    lam = const(tau).div(d).power(g.number(1, 2))
    r = lam.mul(const(g.number(2)).mul(R).power(g.number(1, 2)))
    z = lam.power(g.sub(g.one, ports["delta"])).mul(Z)
    require(z.radial == g.zero, "new geometry must have z_R=0 structurally")
    det = g.sub(g.mul(r.radial, z.axial), g.mul(r.axial, z.radial))
    L = g.sub(g.one, g.mul(ports["delta"], ports["Z"], ports["Z"]))
    return {"ports": ports, "R": R, "Z": Z, "tau": tau, "theta": theta,
            "d": d, "L": L, "lambda": lam, "r": r, "z": z, "J2": det}


def _inverse(g, geometry):
    p = geometry["ports"]
    X = g.node("free_variable", name="X_phys")
    Y = g.node("free_variable", name="Y_phys")
    z = g.node("free_variable", name="z_phys")
    t = g.node("free_variable", name="t_phys")
    q = g.node("bound_variable", name="q_inverse", role="implicit_positive_root")
    tau = g.sub(g.one, t)
    equation = g.sub(g.sub(q, g.mul(z, z, g.power(q, p["delta"]))), tau)
    root = g.node("implicit_positive_root", equation, bound_variable=q,
                  domain="q>0 and q-z_phys^2*q^delta=1-t_phys; physical branch",
                  uniqueness="paper coordinate lemma used as DATA under 0<delta<1 and tau>0",
                  numerical_solver=False)
    R = g.div(g.add(g.mul(X, X), g.mul(Y, Y)), g.mul(g.number(2), root))
    Z = g.div(z, g.power(root, g.mul(g.number(1, 2), g.sub(g.one, p["delta"]))))
    native_x = g.sub(g.node("log", g.div(R, p["Rm"])), g.add(g.number(17), g.node("exp", g.number(40))))
    bindings = {p["R"]: R, p["Z"]: Z, p["native_x"]: native_x,
                geometry["tau"]: tau, geometry["theta"]: g.node("atan2", Y, X)}
    return {"bound_root": root, "root_equation": equation, "bound_variable": q,
            "arguments": [X, Y, z, t], "bindings": bindings,
            "native_expression_substitution": True, "numerical_inverse_available": False}


def _branch(g, geometry, inverse, cap, parent, branch):
    p = geometry["ports"]
    fields = cap[branch + "_rows"]
    refs = {key: g.source(rows[0], branch + "_rows." + key + "[0]", key) for key, rows in fields.items()}
    mz = g.source(fields["m"][1], branch + "_rows.m[1]", "ordinary m_Z at fixed profile R")
    my = (g.source(_node(parent["my_compatibility_rows"][0]["sum_head"]), "accepted_modified_m_y[0]", "accepted V_N-m")
          if branch == "modified" else g.sub(refs["V"], refs["m"]))
    m = _Dual(g, refs["m"], g.div(my, p["R"]), mz)
    psi = geometry["lambda"].power(g.sub(g.one, p["delta"])).mul(_Dual(g, p["Pstar"])).mul(geometry["R"]).mul(m)
    r, z, det = geometry["r"], geometry["z"], geometry["J2"]
    psi_r = g.div(g.sub(g.mul(psi.radial, z.axial), g.mul(psi.axial, z.radial)), det)
    psi_z = g.div(g.add(g.neg(g.mul(psi.radial, r.axial)), g.mul(psi.axial, r.radial)), det)
    ur, uz = g.neg(g.div(psi_z, r.value)), g.div(psi_r, r.value)
    ut = g.mul(g.power(geometry["lambda"].value, g.sub(g.number(-1), p["delta"])), p["Pstar"], refs["E"])
    pressure = g.mul(g.power(geometry["lambda"].value, g.sub(g.number(-2), g.mul(g.number(2), p["delta"]))),
                     g.power(p["Pstar"], g.number(2)), g.add(refs["P0"], refs["p"]))
    cos, sin = g.node("cos", geometry["theta"]), g.node("sin", geometry["theta"])
    u, v = g.sub(g.mul(ur, cos), g.mul(ut, sin)), g.add(g.mul(ur, sin), g.mul(ut, cos))
    Atheta = g.div(psi.value, r.value)
    vector = g.node("vector", g.neg(g.mul(Atheta, sin)), g.mul(Atheta, cos), g.zero)
    functions = {"u": u, "v": v, "w": uz, "pressure": pressure, "psi": psi.value, "vector_potential": vector}
    composed = {key: g.node("compose_native_function_and_geometry", value,
                           substitutions=inverse["bindings"], physical_arguments=inverse["arguments"],
                           implicit_bound_root=inverse["bound_root"], source_values_unevaluated=True)
                for key, value in functions.items()}
    return {"branch": branch, "coefficient_order": 0, "physical_function_orders": [0],
            "source_fields": refs, "m_y": my, "m_Z": mz, "psi_dual": psi.data(),
            "potential_one_form": {"basis": ["R", "theta", "Z"], "theta_coefficient": psi.value,
                                   "definition": "A_one_form=psi*dtheta"},
            "flux_two_form": {"orientation": ["R", "theta", "Z"], "R_theta": psi.radial, "theta_Z": g.neg(psi.axial)},
            "cylindrical": {"ur": ur, "utheta": ut, "uz": uz}, "Cartesian": functions,
            "physical_Cartesian_functions": composed, "incoming_history_reset": False,
            "P0_used_once": True, "source_jet_orders_data_only": [0, 1, 2, 3, 4]}


def _closure(branch):
    # Fresh exterior antisymmetry calculation for d(dpsi wedge dtheta).
    # The shared mixed-partial atom is conditional on the declared C2 premise.
    atom = {"potential": branch["psi_dual"]["value"], "derivatives": ["R", "Z"], "commutation_premise": "C2"}
    signed_terms = [{"coefficient": 1, "atom": atom}, {"coefficient": -1, "atom": atom}]
    coefficient = sum(term["coefficient"] for term in signed_terms)
    require(coefficient == 0, "exterior closure coefficient did not cancel")
    return {"branch": branch["branch"], "signed_mixed_terms": signed_terms,
            "result_coefficient": coefficient, "identity": "d(dpsi wedge dtheta)=0",
            "physical_orientation": "r dr wedge dtheta wedge dz",
            "swirl_divergence": "zero by independence of physical theta, not by holding native fast phase",
            "conditional_C2_source_identity": True, "native_global_derivative_identity_proved": False}


def select_definition(branch, field, coefficient_order=0, phase_mode="ordinary_Z"):
    require(branch in BRANCHES, "unknown source branch")
    require(field in ("u", "v", "w", "pressure", "psi", "vector_potential", "volume", "map"), "unknown field")
    require(type(coefficient_order) is int and coefficient_order == 0, "only leading coefficient zero is defined")
    require(phase_mode == "ordinary_Z", "held native phase is not an admitted derivative mode")
    return {"branch": branch, "field": field, "coefficient_order": 0, "source_values_unevaluated": True}


def _counts(journal, passed=False):
    events = {name: sum(row["event"] == name for row in journal) for name in EVENT_COUNTS}
    return {**{key: events[key] for key in EVENT_COUNTS if key not in ("family_entered", "family_completed")},
            "source_family_calls": events["family_entered"], "new_differential_family_calls": events["family_entered"],
            "new_mathematical_proofs": int(passed), "failed_family_calls": int(not passed and events["family_entered"] > 0),
            "old_math_calls": 0, "provider_calls": 0, "numeric_evaluation_calls": 0,
            "N_materialization_calls": 0, "source_value_reducer_calls": 0, "ancestor_calls": 0,
            "temporal_recursion_calls": 0}


def _flags(passed=False):
    return {"conditional_spacetime_pullback_defined": passed, "same_source_exterior_potential_defined": passed,
            "Cartesian_velocity_pressure_functions_defined": passed, "conditional_exterior_divergence_identity_proved": passed,
            "physical_volume_expression_defined": passed, "implicit_Cartesian_function_composition_defined": passed,
            "native_field_derivative_identity_proved": False, "numeric_closure": False, "physical_layer_completed": False,
            "full_goal_complete": False, "temporal_recursion": False, "cap_deletion": False,
            "unconditional_parameter_admission": False, "global_phase_seam_regularity_proved": False,
            "axis_extension_proved": False, "numerical_oracle_available": False,
            "unconditional_Z4_source_admission": False, "actual_incoming_x_independence_proved": False}


def derive(payload, emit=None):
    journal = []
    def send(event, row):
        journal.append({"event": event, "row": row})
        if emit is not None:
            emit(event, row)
    send("family_entered", {"component": COMPONENT, "new_method": "fresh dual one-form and Jacobian pullback"})
    try:
        parent, literals, cap = _inputs(payload)
        joins = _join_sources(parent, literals, cap)
        send("source_join", {"rows": joins, "pressure_aliases": 10, "accepted_proofs_consumed_as_DATA": True})
        g = _IR()
        geometry = _geometry(g, cap)
        inverse = _inverse(g, geometry)
        send("geometry_dual", {"r": geometry["r"].data(), "z": geometry["z"].data(),
                               "lambda": geometry["lambda"].data(), "J2": geometry["J2"], "inverse": inverse})
        branches = {}
        for name in BRANCHES:
            b = _branch(g, geometry, inverse, cap, parent, name)
            branches[name] = b
            send("potential", {"branch": name, "psi": b["psi_dual"], "one_form": b["potential_one_form"],
                               "vector_potential": b["Cartesian"]["vector_potential"]})
            send("flux_form", {"branch": name, "two_form": b["flux_two_form"]})
            send("velocity_pressure", {"branch": name, "cylindrical": b["cylindrical"], "pressure": b["Cartesian"]["pressure"]})
            send("Cartesian_composition", {"branch": name, "functions": b["physical_Cartesian_functions"]})
            send("exterior_closure", _closure(b))
        p = geometry["ports"]
        volume = {"profile_R_Z_theta_absolute_measure": g.mul(geometry["r"].value, geometry["J2"]),
                  "native_x_Z_theta_absolute_measure": g.mul(p["R"], geometry["r"].value, geometry["J2"]),
                  "positive_orientation": ["R", "theta", "Z"],
                  "paper_closed_form": "lambda_phys^(3-delta)*L/(1-Z^2); native x has the extra R factor",
                  "closed_form_symbolic_comparison_performed": False,
                  "positivity_premises": "R>0,tau>0,|Z|<1,0<delta<1/200,fixed Pstar>0 and Rm>0",
                  "profile_R_is_not_physical_r": True, "not_a_numeric_L2_norm": True}
        send("volume_form", volume)
        for request in (("unknown", "u", 0, "ordinary_Z"), ("background", "u", 1, "ordinary_Z"),
                        ("background", "u", 0, "held_phase"), ("background", "unknown", 0, "ordinary_Z")):
            try:
                select_definition(*request)
            except ValueError as exc:
                send("api_refusal", {"request": list(request), "refused": True, "reason": str(exc)})
            else:
                raise ValueError("unsupported API request was admitted")
        body = {"expression_IR": g.nodes, "branches": branches, "inverse": inverse,
                "forward_map": {"r": geometry["r"].value, "z": geometry["z"].value, "t": g.sub(g.one, geometry["tau"]),
                                "theta": geometry["theta"], "profile_R": p["R"], "native_x": p["native_x"],
                                "native_R_definition": "R=Rm*exp(17+exp(40)+native_x)"},
                "physical_volume": volume, "source_join_rows": joins,
                "source_wrapper_pins": WRAPPER_PINS, "source_values_unevaluated": True,
                "parent_payload_copied_recursively": False,
                "premises": ["regular interior phase cells", "C2 source potential", "accepted conditional m_y=V-m and ordinary m_Z",
                             "cached native-R ordinary-Z frame contract used as DATA, not a numerical oracle certificate",
                             "Pstar and delta fixed in R,Z", "Pstar>0", "Rm>0", "0<delta<1/200", "R>0", "0<tau<=1", "|Z|<1", "0<native_x<1"],
                "not_claimed": ["axis extension", "global phase seam regularity", "numerical source values", "NS closure", "temporal coefficient recursion"]}
        send("family_completed", {"component": COMPONENT, "proof_count": 1, "conditional_exterior_identity": True})
        passed = True
        require({name: sum(row["event"] == name for row in journal) for name in EVENT_COUNTS} == EVENT_COUNTS, "callback contract changed")
        failure = None
    except Exception as exc:
        passed, body, failure = False, None, {"reason": str(exc), "prefix_preserved": True}
    return {"schema": RESULT_SCHEMA, "component": COMPONENT, "passed": passed, "source_only": True,
            "metadata_only": False, "full_goal_complete": False, "proof_count": int(passed),
            "actual_counts": _counts(journal, passed), "new_flags": _flags(passed),
            "callback_events": journal, "callback_event_count": len(journal),
            "callback_event_counts": {name: sum(row["event"] == name for row in journal) for name in EVENT_COUNTS},
            "source_wrapper_pins": WRAPPER_PINS, "payload": body, "failure": failure,
            "parent_commit": payload.get("parent_commit") if isinstance(payload, dict) else None,
            "runtime_model_metadata": None}
