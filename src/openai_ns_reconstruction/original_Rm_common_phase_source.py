"""Exact source binding of the common Rm phase and saved periodic theorem.

This is a new source operator, not a numerical phase evaluation.  It consumes
the accepted radius definition and a saved universal periodic-loop theorem.
No derivative of fractional_part is taken and no inverse conditioning or
native loop derivative is claimed.
"""
from __future__ import annotations

import copy
import hashlib
import json


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                    ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def prove(nodes, left, right, radius_lemma, saved_lift, seed_admission, generic_packet, *, loop_source_sha256):
    """Bind the actual two endpoint phases to the same unfolded phase.

    The radius lemma has already established R_left(1)=Rm and R_right(x)=Rm*x
    with D_y=x*d/dx.  Congruence of quotient/log/sum/product/fractional_part
    gives the common endpoint value.  The *unfolded* phase has D_y theta=N;
    periodic smooth composition, rather than the discontinuous representative,
    supplies the interpretation at integer cuts.
    """
    require(radius_lemma.get("verified") is True and radius_lemma.get("status") == "proved",
            "accepted analytic radius source required")
    descriptors = radius_lemma["exact_source_descriptors"]
    require(all(nodes[int(index)] == row for index, row in descriptors.items()),
            "radius source descriptor changed")
    require(left["coordinate"] == 84528 and right["coordinate"] == 185676
            and left["radius"] == 84533 and right["radius"] == 187214,
            "foreign radius chart")
    require(nodes[187214] == {"operation": "product", "arguments": [85743, 185676]}
            and nodes[85743] == {"operation": "function_substitution", "expression": 84533,
                                "variable": 84528, "value": 1, "Z_independent_substitution": True},
            "endpoint radius identity not bound")
    require(nodes[187963] == {"operation": "function_substitution", "expression": 187214,
                             "variable": 185676, "value": 1, "Z_independent_substitution": True},
            "right closed radius endpoint changed")
    require(nodes[24954] == dict(operation="product", arguments=[8, 10, 28])
            and nodes[65604] == dict(operation="product", arguments=[13, 35565]),
            "common origin/Ra source changed")
    phase_rows = []
    ids = set(int(index) for index in descriptors)
    for chart, offset, logref, quotient, product, phase in (
        (left, 119414, 119413, 119412, 119415, 119416),
        (right, 196044, 196043, 196042, 196045, 196046),
    ):
        require(chart["global_offset"] == offset and chart["actual_phase"] == phase,
                "actual phase binding changed")
        require(nodes[phase] == dict(operation="analytic_unary", name="fractional_part", argument=product)
                and nodes[product] == dict(operation="product", arguments=[18, offset])
                and nodes[offset] == dict(operation="sum", arguments=[logref, 97298])
                and nodes[logref] == dict(operation="analytic_unary", name="log", argument=quotient),
                "actual global phase recipe changed")
        q = nodes[quotient]
        require(q["operation"] == "positive_quotient" and q["numerator"] == chart["radius"]
                and q["denominator"] == 65604 and bool(q.get("source_positive_certificate")),
                "positive R/Ra source changed")
        phase_rows.append(dict(coordinate=chart["coordinate"], endpoint=1, radius=chart["radius"],
                               global_offset=offset, actual_phase=phase, common_N=18,
                               common_Ra=65604, common_phase_origin=24954))
        ids.update((offset, logref, quotient, product, phase))
    require(nodes[97298] == dict(operation="negative", argument=24954), "phase origin sign changed")
    require(nodes[18]["operation"] == "exact_positive_integer_expression"
            and nodes[18].get("Z_independent") is True
            and nodes[18].get("one_common_N_for_all_charts_and_control_denominators") is True,
            "common positive N changed")
    # These are the complete actual phase dependencies, not a metadata claim
    # that a chart's total Z derivative happens to vanish.
    zfree = set()
    def visit(index):
        if index in zfree:
            return
        row = nodes[index]
        operation = row["operation"]
        if operation == "exact_rational":
            children = []
        elif operation in ("current_original_source_parameter", "exact_positive_integer_expression"):
            require(row.get("Z_independent") is True, "phase source leaf is not Z independent")
            children = []
        elif operation == "bound_variable":
            require(index in (84528, 185676), "foreign free variable in actual phase")
            children = []
        elif operation in ("sum", "product"):
            children = row["arguments"]
        elif operation in ("analytic_unary", "negative"):
            children = [row["argument"]]
        elif operation == "positive_quotient":
            children = [row["numerator"], row["denominator"]]
        elif operation == "function_substitution":
            require(row.get("Z_independent_substitution") is True
                    and row["variable"] in (84528, 185676), "unguarded phase substitution")
            children = [row["expression"], row["value"], row["variable"]]
        else:
            raise ValueError("unhandled phase dependency: " + operation)
        for child in children:
            require(isinstance(child, int) and 0 <= child < index, "phase graph is not backward")
            visit(child)
        zfree.add(index)
    visit(119416)
    visit(196046)
    ids.update(zfree)
    cert = saved_lift["phase_lift_certificate"]
    require(cert["passed"] is True and cert["global_periodic_lift_proved"] is True
            and cert["mixed_phase_parameter_jets_through_order2_proved"] is True,
            "saved universal periodic theorem missing")
    pin = cert["source_guard"]["loop_pin"]
    require(pin["passed"] is True and pin["sha256"] == loop_source_sha256
            and cert["source_ast"]["passed"] is True and cert["formal_identities"]["passed"] is True,
            "saved theorem is bound to another loop recipe")
    checks = cert["formal_identities"]["checks"]
    for name in ("phase_period", "K_positive", "A_period_difference", "B_period_difference",
                 "parameter_and_mixed_jet_differences", "phase_derivative_endpoint_formulas",
                 "strict_r_inside_unit_interval", "q_zero_included", "flat_A_is_zero_function",
                 "flat_B_is_zero_function", "flat_Aphi_is_zero_function", "flat_Bphi_is_zero_function",
                 "fixed_q_zero_tangent_mixed_jets", "positive_gamma_at_flat_boundary"):
        require(checks[name] is True, "saved theorem lacks " + name)
    require(cert["source_guard"]["required_checks"]["parameters"] is True,
            "saved common positive parameter admission missing")
    require(seed_admission.get("passed") is True and generic_packet.get("passed") is True,
            "actual seed/generic admissions required for loop domain")
    fields = {(r["field"], r["order"]): r for r in seed_admission["fields"]}
    e, ey = fields["E", 0], fields["E_y", 0]
    positive_monomial = [["exp_rational", -3, 5, 1],
                        ["inverse_expression", [[[["atom", "symbol:Z", 2]], 1, 1], [[], 1, 1]], 1]]
    require(e["status"] == ey["status"] == "proved"
            and e["inlet"] == e["patch"] == dict(terms=[dict(monomial=positive_monomial, coefficient="1")])
            and ey["inlet"] == ey["patch"] == dict(terms=[dict(monomial=positive_monomial, coefficient="1/10")]),
            "saved actual positive E / E_y=E/10 source forms changed")
    a = next(r for r in generic_packet["root_rows"] if r["root"] == "a" and r["order"] == 0)
    require(a["status"] == "proved"
            and a["postrestore"]["root_node"] == left["roots"]["a"][0]
            and a["patch"]["root_node"] == right["roots"]["a"][0]
            and a["postrestore"]["value"]["text"] == "(SEED:E[0] - 2*SEED:E_y[0])/SEED:E[0]"
            and a["postrestore"]["value"] == a["patch"]["value"],
            "actual generic positive a definition changed")
    ids.update((18, 97298, 24954, 65604))
    return dict(schema="actual_Rm_common_phase_source_v1", passed=True, source_only=True,
        actual_common_phase_proved=True, smooth_periodic_lift_bound=True,
        phase_nodes=[119416, 196046], coordinates=[84528, 185676], endpoints=[1, 1],
        phase_rows=phase_rows, exact_source_descriptors={str(i):copy.deepcopy(nodes[i]) for i in sorted(ids)},
        radius_lemma_digest=digest(radius_lemma), saved_lift_certificate_digest=digest(cert),
        source_graph_digest=digest(nodes), seed_admission_digest=digest(seed_admission),
        endpoint_radius_reduction=dict(left={"closed_node":85743, "rule":"guarded R_left(p=1)=Rm"},
            right={"closed_node":187963, "body_node":187214, "rule":"(Rm*x)|x=1=Rm*1=Rm"},
            common_radius_source_node=85743),
        same_loop_source_sha256=loop_source_sha256,
        endpoint_identity="frac(N*(log(Rm/Ra)-phase_origin)) on both sides",
        unfolded_native_identity="D_y[N*(log(R/Ra)-phase_origin)]=N; D_y=d/dp on left and x*d/dx on right",
        proof_rules=["guarded radius endpoint substitution", "congruence of the displayed actual phase constructors",
                     "elementary log chain rule on positive R/Ra", "saved smooth periodic composition theorem"],
        periodic_theorem_scope="universal active loop a,E,eta,d_star>0; common lazy-flat extension; parameter/phase mixed jets through order2",
        actual_loop_domain=dict(E="exp(-3/5)/(1+Z^2)>0 for every real Z", E_y="E/10",
            a="(E-2*E_y)/E=4/5>0", eta_dstar="same positive canonical parameters of the saved theorem, same graph nodes21/22",
            cutoff="if Delta<eta then 2*eta-Delta>eta>0; otherwise same lazy zero branch",
            arbitrary_real_p2_allowed=True, quantitative_inverse_conditioning=False),
        ordinary_orders=[0, 1], fractional_part_differentiated=False,
        actual_phase_Z_independence_proved=True, Z_independent_dependency_nodes=sorted(zfree),
        numerical_phase_value_computed=False, inverse_conditioning_proved=False,
        native_loop_derivative_proved=False, density_seam_proved=False,
        all_212_terms_and_errors_unchanged=True, full_goal_complete=False)
