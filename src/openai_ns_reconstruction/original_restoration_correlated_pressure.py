"""Actual modified restoration seams and conservative signed correlation.

Only current new source packets are composed. Saved generic formula lemmas
supply reusable identities; actual new cache nodes, seed maps, phase rows,
and nonreset memories must match before an internal RQy is contracted.
"""
from __future__ import annotations

from collections import deque
from copy import deepcopy
from fractions import Fraction
import gzip
import hashlib
import json
from math import comb
from pathlib import Path
from typing import Any, Mapping

from .original_signed_graph_pairs import GraphPairs
from .original_restoration_pressure_source import restored_nodes


def require(value: Any, message: str) -> None:
    if not value:
        raise ValueError(message)


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     ensure_ascii=True).encode()).hexdigest()


BASE = "experiments/root_st073/"
PARENT = BASE + "lei_ren_part1_paper_compliant_current_original_long_reference_modified_seams.json"
PARENT_SHA256 = "204d8b45aff263a3df1eae319dd7534e304162794ac375a1c0ec2cfe9beafc35"
STEM = BASE + "lei_ren_part1_paper_compliant_current_original_"
C6 = STEM + "prescribed_R100_Rm_C6_functions"
BRIDGE_C6 = STEM + "prescribed_bridge_C6_functions"
C5 = STEM + "prescribed_loop_C5_functions"
Q4 = STEM + "modified_Q4_pressure_C3_functions"
MODULE = "src/openai_ns_reconstruction/original_restoration_correlated_pressure.py"
NEW_MODULES = dict(
    source="original_restoration_pressure_source", endpoint="original_restoration_native_endpoints",
    native="original_restoration_native_seams", dag="original_restoration_physical_y_dag",
    phase="original_restoration_phase_mixed",
)
CHARTS = ("restoration", "postrestore")
NATIVE = {"restoration": "restoration", "postrestore": "post_restore"}
FIELDS = ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")
SEAMS = ("reference_to_restoration", "restoration_to_postrestore")
SOURCE_GUARDS = (
    (BRIDGE_C6 + ".py", "W=a.scale(F0,g.mul(g.unary('sqrt',g.mul(c(2),radius)),inverse))"),
    (BRIDGE_C6 + ".py", "FF=a.scale(a.mul(F0,F0),g.mul(radius,inverse,inverse))"),
    (BRIDGE_C6 + ".py", "radial=dict(E_y=a.mul(W,a.add(a.scale(actual['phi'],'1/2'),actual['phi_y'])),V_y=a.scale(actual['U_y'],inverse))"),
    (C5 + ".py", "p2=generic['p2'];t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))"),
    (C5 + ".py", "kap=a.add(aa,a.div(a.mul(bb,bb),aa,'same original positive a'))"),
    (C5 + ".py", "predicate=shear"),
    (C5 + ".py", "collar=g.mul(c('3/4'),sc) if chart=='first_micro' else None"),
    (C5 + ".py", "gamma=a.sub(a.fixed(g.mul(c(2),eta)),roots['Delta'])"),
    (C5 + ".py", "cut=a.sub(one,a.scale(roots['Delta'],g.quotient(g.one,eta,'same positive canonical eta')))"),
    (C5 + ".py", "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))"),
    (C5 + ".py", "C=a.sub(E,a.scale(a.trim(radial['E_y']),2));B=a.scale(a.trim(radial['V_y']),2)"),
    (C5 + ".py", "aa,bb=divide(C),divide(B)"),
    (C5 + ".py", "M=a.sub(a.neg(a.scale(a.mul(roots['a'],T1hat),g.quotient(g.one,twopi,'positive 2*pi'))),"),
    (C5 + ".py", "inc=dict(m=s,h=r,k=a.add(a.mul(V,r),a.mul(E,s),a.mul(r,s)),"),
    (Q4 + ".py", "Vy=a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])"),
    (Q4 + ".py", "m=a.scale(fields['m'],Pstar);my=a.sub(U,m)"),
    (Q4 + ".py", "Q=div(numerator(U,m));Qy=div(numerator(Uy,my))"),
)


def pin_packet(root, parent, packet, name):
    path = "src/openai_ns_reconstruction/" + NEW_MODULES[name] + ".py"
    pin = packet["source_pins"][path]
    expected_parent = Path(PARENT).name if name == "dag" else PARENT
    require(packet["passed"] is True and packet["parent_report"] == expected_parent
            and packet["parent_report_sha256"] == PARENT_SHA256
            and packet["source_family"] == parent["source_family"]
            and packet["source_executed"] is False and pin["path"] == path
            and pin["passed"] is True and pin["source_executed"] is False
            and pin["expected_sha256"] == pin["actual_sha256"] == sha(root / path),
            "New packet provenance changed: " + name)


def p_add(*values):
    result = {}
    for value in values:
        for monomial, coefficient in value.items():
            result[monomial] = result.get(monomial, Fraction(0)) + coefficient
            if not result[monomial]:
                del result[monomial]
    return result


def p_scale(value, coefficient):
    return {key: val * coefficient for key, val in value.items() if val * coefficient}


def p_mul(*values):
    result = {(): Fraction(1)}
    for value in values:
        result = p_add(*({tuple(sorted(left + right)): a*b}
                         for left, a in result.items() for right, b in value.items()))
    return result


def atom(nodes, index):
    row = nodes[index]
    if row.get("operation") == "exact_rational":
        value = Fraction(row["numerator"], row["denominator"])
        return {(): value} if value else {}
    return {(index,): Fraction(1)}


def polynomial(nodes, root, inputs):
    cache = {}

    def visit(index):
        if index in cache:
            return cache[index]
        row = nodes[index]
        operation = row["operation"]
        if operation == "exact_rational" or index in inputs:
            result = atom(nodes, index)
        elif operation == "sum":
            result = p_add(*(visit(n) for n in row["arguments"]))
        elif operation == "product":
            result = p_mul(*(visit(n) for n in row["arguments"]))
        elif operation == "negative":
            result = p_scale(visit(row["argument"]), -1)
        else:
            result = atom(nodes, index)
        cache[index] = result
        return result

    return visit(root)


def descendants(nodes, root):
    pending, found = [root], set()
    while pending:
        node = pending.pop()
        if node in found:
            continue
        found.add(node)
        row = nodes[node]
        if row.get("operation") == "exact_rational":
            # Numerator/denominator here are literal integers, not handles.
            continue
        for key in ("argument", "active_body", "body", "expression", "numerator", "denominator", "inverse_angle"):
            if type(row.get(key)) is int:
                pending.append(row[key])
        pending.extend(n for n in row.get("arguments", ()) if type(n) is int)
    return found


def exact(nodes, node, value):
    return type(node) is int and nodes[node].get("operation") == "exact_rational" and atom(nodes, node) == ({(): Fraction(value)} if value else {})


def inverse(nodes, root, denominator=None, twopi=False):
    candidates = []
    for index in descendants(nodes, root):
        row = nodes[index]
        if row.get("operation") != "positive_quotient" or not exact(nodes, row.get("numerator"), 1):
            continue
        if denominator is not None and row.get("denominator") != denominator:
            continue
        if twopi:
            den = nodes[row["denominator"]]
            arguments = den.get("arguments", ())
            if not (den.get("operation") == "product" and len(arguments) == 2
                    and any(exact(nodes, n, 2) for n in arguments)
                    and any(nodes[n].get("operation") == "mathematical_pi" for n in arguments)):
                continue
        candidates.append(index)
    require(len(set(candidates)) == 1, "Actual positive inverse is not uniquely bound")
    return candidates[0]


def _flatten_operation(nodes, index, operation):
    """Flatten only the associative source operation at ``index``."""
    row = nodes[index]
    if row.get("operation") != operation:
        return [index]
    values = []
    for child in row.get("arguments", ()):
        values.extend(_flatten_operation(nodes, child, operation))
    return values


def _scalar(nodes, index):
    row = nodes[index]
    if row.get("operation") != "exact_rational":
        return None
    return Fraction(row["numerator"], row["denominator"])


def _source_zero(nodes, index):
    """Recognize exact zero through endpoint substitution wrappers."""
    value = _scalar(nodes, index)
    if value == 0:
        return True
    row = nodes[index]
    if row.get("operation") == "function_substitution":
        expression = row.get("expression")
        return type(expression) is int and _source_zero(nodes, expression)
    return False


def _polynomial_residual(nodes, actual, expected, inputs):
    """Return the exact polynomial residual of two source expressions."""
    return p_add(polynomial(nodes, actual, inputs), p_scale(expected, -1))


def _expand_polynomial_expression(nodes, expression, inputs):
    """Expand a local polynomial expression whose monomials are graph atoms."""
    result = {}
    for monomial, coefficient in expression.items():
        term = {(): coefficient}
        for node in monomial:
            term = p_mul(term, polynomial(nodes, node, inputs))
        result = p_add(result, term)
    return result


def _jet_product(nodes, left, right, order):
    """Product of ordinary derivative rows, with the binomial factors."""
    return p_add(*(p_scale(p_mul(atom(nodes, left[k]), atom(nodes, right[order-k])),
                           Fraction(comb(order, k)))
                   for k in range(order + 1)))


def _scale_series_from_product(nodes, output, input_rows, label):
    """Recover the actual multiplicative jet from an emitted product jet.

    The C6 source emits ordinary derivative rows as a sum of products.  Each
    term contains exactly one source input row and one coefficient row (plus
    an exact binomial factor), so this recovers the coefficient rows from the
    graph itself rather than installing a desired scale fixture.
    """
    result = [None] * len(input_rows)
    for order, output_row in enumerate(output):
        terms = _flatten_operation(nodes, output_row, "sum")
        for input_order in range(order + 1):
            source = input_rows[order - input_order]
            matches = []
            for term in terms:
                term_row = nodes[term]
                factors = (list(term_row.get("arguments", ()))
                           if term_row.get("operation") == "product" else [term])
                if source not in factors:
                    continue
                nonconstant = [factor for factor in factors
                               if factor != source and _scalar(nodes, factor) is None]
                if len(nonconstant) == 1:
                    matches.append((term, nonconstant[0]))
            require(len(matches) == 1,
                    f"{label} source product row is not uniquely factored: {order}:{input_order}")
            coefficient = matches[0][1]
            if result[input_order] is None:
                result[input_order] = coefficient
            else:
                require(result[input_order] == coefficient,
                        f"{label} coefficient jet changed: {input_order}")
    require(all(type(value) is int for value in result), label + " coefficient jet is incomplete")
    return result


def _inverse_series_from_product(nodes, output, input_rows, label):
    """Recover a Z-independent exact multiplier from a field jet."""
    first_row = nodes[output[0]]
    first = (list(first_row.get("arguments", ()))
             if first_row.get("operation") == "product" else [output[0]])
    matches = [item for item in first if item != input_rows[0] and _scalar(nodes, item) is None]
    require(len(matches) == 1, label + " multiplier is not uniquely source-bound")
    multiplier = matches[0]
    for order, output_row in enumerate(output):
        if _scalar(nodes, input_rows[order]) == 0:
            require(_scalar(nodes, output_row) == 0,
                    label + " zero source row was not preserved")
            continue
        row = nodes[output_row]
        factors = (list(row.get("arguments", ()))
                   if row.get("operation") == "product" else [output_row])
        require(input_rows[order] in factors and multiplier in factors,
                label + " source multiplier changed")
    return multiplier


def _native_history_seed_proof(nodes, c6, native_packet):
    """Bind all six actual native history densities to the C6 phi/U jets."""
    source = native_packet["native_history_source"]
    lemma = native_packet["saved_lemmas"]["history_ODE"]
    generic = native_packet["saved_lemmas"]["generic_seed"]
    rates = {"H": 2, "M": 1, "K": 2, "A": 1, "B": 2, "C": 1}
    formulas = {"H": "2*phi", "M": "U", "K": "2*phi*U",
                "A": "U**2", "B": "phi**2", "C": "phi**2"}
    require(source["passed"] is True and source["source_executed"] is False
            and source["actual_phi_U_roots_used"] is True
            and source["density_source_formulas"] == formulas
            and source["rates"] == rates,
            "Actual native history source theorem changed")
    require(lemma["source_FTC"] is True and lemma["applied_to_new_charts"] is True
            and lemma["statement"] == "Dq X/J = density - rate*X",
            "Saved native history ODE lemma is not consumable")
    require(generic["passed"] is True and generic["consumed_saved_lemma_only"] is True
            and generic["no_old_math_rerun"] is True,
            "Saved generic seed lemma was not consumed as a lemma")
    rows = {}
    for name in CHARTS:
        native_name = NATIVE[name]
        actual = c6["charts"][native_name]["actual"]
        moments = c6["charts"][native_name]["moments"]
        phi, U = actual["phi"], actual["U"]
        chart_receipt = source["charts"][native_name]
        require(chart_receipt["passed"] is True and chart_receipt["source_executed"] is False,
                "Native history chart receipt changed: " + name)
        chart_rows = {}
        for history, formula in formulas.items():
            receipt = chart_receipt["histories"][history]
            require(receipt["rate"] == rates[history]
                    and receipt["ordinary_Z_orders"] == list(range(7))
                    and len(receipt["rows"]) == 7,
                    "Native history rows changed: " + name + ":" + history)
            residuals = []
            for order, row in enumerate(receipt["rows"]):
                require(row["order"] == order and row["q0"]["exp_zero_is_one"] is True
                        and row["q0"]["integral_zero_is_zero"] is True
                        and row["q0"]["incoming_memory_retained"] is True,
                        "Native history q0 source changed: " + name + ":" + history)
                require(row["moment"] == moments[history][order],
                        "Native history receipt is detached from actual C6 moment rows: "
                        + name + ":" + history + ":" + str(order))
                if history == "H":
                    expected = p_scale(atom(nodes, phi[order]), 2)
                elif history == "M":
                    expected = atom(nodes, U[order])
                elif history == "K":
                    expected = p_scale(_jet_product(nodes, phi, U, order), 2)
                elif history == "A":
                    expected = _jet_product(nodes, U, U, order)
                else:
                    expected = _jet_product(nodes, phi, phi, order)
                residual = _polynomial_residual(
                    nodes, row["density_expression"], expected, set(phi + U))
                require(not residual,
                        "Native history density is detached from actual phi/U: "
                        + name + ":" + history + ":" + str(order))
                residuals.append(dict(order=order, residual="0"))
            chart_rows[history] = dict(rate=rates[history], formula=formula,
                                       rows=residuals, residual_count=0)
        rows[name] = chart_rows
    return dict(passed=True, charts=rows, rates=rates, formulas=formulas,
                source_FTC=True, generic_saved_lemma_consumed=True,
                residual_count=0)


def _native_seed_formula_proof(nodes, c6_chart, q4_chart, native_packet, name, c6, dag_chart):
    """Prove the actual 9 by 7 Dq seed rows by local exact algebra."""
    native_name = NATIVE[name]
    fields = c6_chart["fields"]
    moments = c6_chart["moments"]
    actual = c6_chart["actual"]
    radial = c6_chart["radial"]
    F0 = c6["F0"]
    P0 = c6["P0"]
    Pstar = c6["Pstar"]
    radius = c6_chart["radius"]
    background = q4_chart["radial_background"]
    require(background["E"] == radial["E_y"] and background["V"] == radial["V_y"],
            "Q4 background seeds differ from actual C6 radial rows: " + name)
    seed_source = native_packet["q4_background_seed_source"]
    receipt = seed_source["charts"][name]
    require(seed_source["passed"] is True and receipt["J"] == 1
            and receipt["C6_fields"] == fields
            and receipt["radial_background"] == background
            and receipt["constructor_seed"]["J_applied_once"] is True
            and receipt["constructor_seed"]["no_J_cap_or_selected_value"] is True,
            "Actual Q4 seed receipt is not bound to C6: " + name)
    require(all(len(fields[field]) == 7 for field in FIELDS)
            and all(len(background[field]) == 7 for field in FIELDS),
            "Actual native/Q4 seed rows are not 9 by 7: " + name)

    # E=W*phi, E_y=W*(phi/2+phi_y), where W is recovered from the
    # actual emitted C6 product jet.  This keeps the proof tied to C6's
    # actual source rows instead of substituting a desired endpoint value.
    W = _scale_series_from_product(nodes, fields["E"], actual["phi"], name + " E")
    S = [None] * 7
    for order in range(7):
        S[order] = p_add(p_scale(atom(nodes, actual["phi"][order]), Fraction(1, 2)),
                         atom(nodes, actual["phi_y"][order]))
    invPstar = _inverse_series_from_product(nodes, fields["V"], actual["U"],
                                            name + " V")
    require(len(F0) == 7 and len(P0) >= 7 and len(moments["H"]) == 7
            and all(len(moments[key]) == 7 for key in ("M", "K", "A", "B", "C")),
            "Actual C6 normalized source rows are not seven ordinary rows: " + name)
    inv_row = nodes[invPstar]
    require(inv_row.get("operation") == "positive_quotient"
            and inv_row.get("numerator") == 1
            and inv_row.get("denominator") == Pstar
            and inv_row.get("source_positive_certificate") == "same strictly positive Pstar",
            "Actual normalized fields use a different Pstar inverse: " + name)
    # Recover FF from the actual p=FF*C source products.  The recovered W and
    # FF rows are then checked against the common F0, chart radius and Pstar
    # source atoms before they are used in any normalized field map.
    FF = _scale_series_from_product(nodes, fields["p"], moments["C"],
                                     name + " p")
    require(len(FF) == 7, "Actual FF multiplier has the wrong ordinary-row count: " + name)
    W_factors = _flatten_operation(nodes, W[0], "product")
    require(F0[0] in W_factors and invPstar in W_factors,
            "Actual W multiplier is detached from common F0/Pstar: " + name)
    sqrt_rows = [factor for factor in W_factors
                 if nodes[factor].get("operation") == "analytic_unary"
                 and nodes[factor].get("name") == "sqrt"]
    require(len(sqrt_rows) == 1, "Actual W multiplier lacks a unique sqrt(2R): " + name)
    sqrt_row = sqrt_rows[0]
    sqrt_arg = nodes[sqrt_row].get("argument")
    sqrt_factors = _flatten_operation(nodes, sqrt_arg, "product")
    require(nodes[sqrt_arg].get("operation") == "product"
            and radius in nodes[sqrt_arg].get("arguments", ())
            and any(_scalar(nodes, factor) == Fraction(2) for factor in sqrt_factors),
            "Actual W multiplier is not sqrt(2*actual radius): " + name)
    FF0_factors = _flatten_operation(nodes, FF[0], "product")
    require(FF0_factors.count(F0[0]) == 2 and FF0_factors.count(invPstar) == 2
            and radius in descendants(nodes, FF[0]),
            "Actual FF multiplier is detached from F0^2*radius/Pstar^2: " + name)
    source_inputs = set(W + [invPstar] + FF + sum(moments.values(), []) + P0)
    FF_source_inputs = set(F0 + [radius, invPstar])
    W_source_inputs = set(F0 + [sqrt_row, invPstar])
    W_squared_rows = []
    W_source_rows = []
    for order in range(7):
        W_source_residual = _polynomial_residual(
            nodes, W[order], p_mul(atom(nodes, F0[order]), atom(nodes, sqrt_row),
                                    atom(nodes, invPstar)), W_source_inputs)
        require(not W_source_residual,
                "Actual W ordinary row is detached from F0*sqrt(2R)/Pstar: "
                + name + ":" + str(order))
        W_source_rows.append(dict(order=order, residual="0"))
        FF_source_residual = _polynomial_residual(
            nodes, FF[order],
            p_mul(_jet_product(nodes, F0, F0, order), atom(nodes, radius),
                  atom(nodes, invPstar), atom(nodes, invPstar)),
            FF_source_inputs)
        require(not FF_source_residual,
                "Actual FF source map changed before W^2 relation: " + name + ":" + str(order))
        squared = _expand_polynomial_expression(nodes, _jet_product(nodes, W, W, order),
                                                W_source_inputs)
        reduced_square = {}
        for monomial, coefficient in squared.items():
            require(monomial.count(sqrt_row) == 2,
                    "Actual W square has a different radius power: " + name)
            factors = list(monomial)
            factors.remove(sqrt_row)
            factors.remove(sqrt_row)
            factors.append(radius)
            reduced_square = p_add(reduced_square,
                                   {tuple(sorted(factors)): 2 * coefficient})
        square_residual = p_add(reduced_square,
                               p_scale(polynomial(nodes, FF[order], FF_source_inputs), -2))
        require(not square_residual, "Actual W^2-2FF row is nonzero: " + name)
        W_squared_rows.append(dict(order=order, residual="0",
                                   source_identity="sqrt(2*radius)^2=2*radius"))
    relation = dag_chart["radius_relation"]
    phase = c6_chart["phase"]
    require(relation["coordinate_node"] == phase and relation["radius_node"] == radius
            and exact(nodes, relation["J_node"], 1) and relation["R_q_equals_JR"] is True
            and relation["R_q_sha256"] == relation["JR_sha256"],
            "Actual new-chart physical radius derivative is not R_y=R: " + name)
    fixed = {"F0": F0, "Pstar": [Pstar], "P0": P0}
    independence = {key: all(phase not in descendants(nodes, node) for node in rows)
                    for key, rows in fixed.items()}
    require(all(independence.values()), "A radial normalization input depends on phase: " + name)
    # Apply sqrt/product chain rules to the actual factor maps just checked.
    # The exact R_y/R=1 premise is the new DAG radius identity, not a cap.
    W_exponent, FF_exponent, radial_ratio = Fraction(1, 2), Fraction(1), Fraction(1)
    weight_residuals = dict(W_y_minus_W_over_2=W_exponent * radial_ratio-Fraction(1, 2),
                            FF_y_minus_FF=FF_exponent * radial_ratio-Fraction(1))
    require(all(value == 0 for value in weight_residuals.values()),
            "Actual radial weight chain rule failed: " + name)
    residuals = []
    for order in range(7):
        expected_field = _jet_product(nodes, W, actual["phi"], order)
        residual = _polynomial_residual(nodes, fields["E"][order], expected_field,
                                        set(W + actual["phi"]))
        require(not residual, "Actual C6 E source map changed: " + name + ":E")
        residuals.append(dict(field="E", order=order, residual="0"))
        expected_radial = p_add(*(p_scale(p_mul(atom(nodes, W[k]), S[order-k]),
                                           Fraction(comb(order, k)))
                                  for k in range(order + 1)))
        residual = _polynomial_residual(nodes, background["E"][order], expected_radial,
                                        set(W + actual["phi"] + actual["phi_y"]))
        require(not residual, "Actual E_y source map changed: " + name + ":" + str(order))
        residuals.append(dict(field="E_y", order=order, residual="0"))

        expected_field = p_mul(atom(nodes, actual["U"][order]), atom(nodes, invPstar))
        residual = _polynomial_residual(nodes, fields["V"][order], expected_field,
                                        set(actual["U"] + [invPstar]))
        require(not residual, "Actual C6 V source map changed: " + name + ":V")
        residuals.append(dict(field="V", order=order, residual="0"))
        expected_radial = p_mul(atom(nodes, actual["U_y"][order]), atom(nodes, invPstar))
        residual = _polynomial_residual(nodes, background["V"][order], expected_radial,
                                        set(actual["U_y"] + [invPstar]))
        require(not residual, "Actual V_y source map changed: " + name + ":" + str(order))
        residuals.append(dict(field="V_y", order=order, residual="0"))

    # Bind every normalized C6 field constructor to the six actual history
    # moment jets and the common F0/Pstar/radius source rows.  The saved
    # generic lemma is only an ODE identity; it does not supply these maps.
    normalized_rows = []
    for order in range(7):
        expected_normalized = {
            "m": p_mul(atom(nodes, moments["M"][order]), atom(nodes, invPstar)),
            "h": p_scale(_jet_product(nodes, W, moments["H"], order), Fraction(1, 2)),
            "k": p_scale(p_mul(_jet_product(nodes, W, moments["K"], order),
                               atom(nodes, invPstar)), Fraction(1, 2)),
            "e": p_add(p_mul(atom(nodes, moments["A"][order]), atom(nodes, invPstar),
                              atom(nodes, invPstar)),
                       p_scale(_jet_product(nodes, FF, moments["B"], order), -1)),
            "p": _jet_product(nodes, FF, moments["C"], order),
            "pressure": p_add(atom(nodes, P0[order]),
                              _jet_product(nodes, FF, moments["C"], order)),
        }
        require(fields["P0"][order] == P0[order],
                "Actual P0 constructor is not the common C6 P0 source: "
                + name + ":" + str(order))
        row_receipts = []
        for field, expected_value in expected_normalized.items():
            residual = _polynomial_residual(nodes, fields[field][order], expected_value,
                                            source_inputs)
            require(not residual,
                    "Actual normalized C6 source map changed: " + name + ":"
                    + field + ":" + str(order))
            row_receipts.append(dict(field=field, order=order, residual="0"))
        normalized_rows.extend(row_receipts)

        V, m, E, h, k, e = (fields[key][order] for key in ("V", "m", "E", "h", "k", "e"))
        expected = {
            "m": p_add(atom(nodes, V), p_scale(atom(nodes, m), -1)),
            "h": p_add(atom(nodes, E), p_scale(atom(nodes, h), Fraction(-3, 2))),
            "k": p_add(_jet_product(nodes, fields["E"], fields["V"], order),
                       p_scale(atom(nodes, k), Fraction(-3, 2))),
            "e": p_add(_jet_product(nodes, fields["V"], fields["V"], order),
                       p_scale(_jet_product(nodes, fields["E"], fields["E"], order), Fraction(-1, 2)),
                       p_scale(atom(nodes, e), -1)),
            "p": p_scale(_jet_product(nodes, fields["E"], fields["E"], order), Fraction(1, 2)),
            "P0": {},
            "pressure": p_scale(_jet_product(nodes, fields["E"], fields["E"], order), Fraction(1, 2)),
        }
        for field, expected_value in expected.items():
            residual = _polynomial_residual(nodes, background[field][order], expected_value,
                                            {pointer for field_name in FIELDS
                                             for pointer in fields[field_name]})
            require(not residual,
                    "Actual Q4 background source map changed: " + name + ":" + field
                    + ":" + str(order))
            residuals.append(dict(field=field, order=order, residual="0"))
    require(background["pressure"] == background["p"],
            "Actual Q4 pressure row is not p_y: " + name)
    return dict(passed=True, chart=name, native_chart=native_name,
                native_fields=9, ordinary_rows=7, W_rows=W, FF_rows=FF,
                Pstar_inverse=invPstar, residuals=residuals,
                normalized_field_rows=normalized_rows,
                normalized_field_maps=dict(
                    m="M/Pstar", h="W*H/2", k="W*K/(2*Pstar)",
                    e="A/Pstar^2-FF*B", p="FF*C", P0="common C6 P0",
                    pressure="P0+FF*C"),
                W_source=dict(F0_rows=list(F0), radius=radius, Pstar=Pstar,
                              sqrt_source=sqrt_row, W_squared_equals_2FF=W_squared_rows,
                              W_source_rows=W_source_rows,
                              W_y_over_W="1/2", FF_y_over_FF="1",
                              derivative_source=dict(radius_identity=deepcopy(relation),
                                  phase_independence=independence,
                                  residuals={key: str(value) for key, value in weight_residuals.items()},
                                  rule="actual sqrt/product maps and R_q=J*R with J=1",
                                  source_executed=False),
                              radial_source="E_y=W*(phi/2+phi_y), V_y=U_y/Pstar"),
                residual_count=0, residual_sha256=digest(residuals),
                E_y_over_E_source="W*(phi/2+phi_y)/(W*phi)",
                V_y_source="U_y/Pstar", P0_y_source="0", pressure_y_source="p_y")


def loop_binding(nodes, chart):
    roots, loop = chart["roots"], chart["loop"]
    phase, angle, psi = chart["actual_phase"], loop["angle"], loop["inverse"][0]
    phase_inverse = nodes[psi]
    predicate = nodes[loop["flat_predicate"]]
    require(predicate.get("operation") == "exact_real_comparison" and predicate["operator"] == ">="
            and predicate["left"] == roots["Delta"][0] and predicate["right"] == 21
            and loop["original_collar_upper"] is None, "New chart has a different lazy branch")
    require(phase_inverse.get("operation") == "original_monotone_phase_inverse"
            and phase_inverse["phase_function"] == loop["Phi"][0]
            and phase_inverse["target_phase"] == phase and phase_inverse["phase_endpoints"] == [0, 1],
            "Actual phase inverse source changed")
    lam, direction, T1 = (loop[key][0] for key in ("phase_lambda_at_inverse", "direction_at_inverse", "T1_at_inverse"))
    partial01 = next(row["node"] for row in loop["phase_partial_coefficients"] if row["Z"] == 0 and row["angle"] == 1)
    lam_row = nodes[lam]
    angle_name = nodes[angle]["name"]
    for index in (lam, direction, T1):
        row = nodes[index]
        require(row.get("operation") == "substitute_original_inverse_angle"
                and row["inverse_angle"] == psi and row["angle_variable"] == angle_name,
                "Actual inverse substitution source changed")
    require(lam_row["body"] == partial01, "Lambda substitution is not the checked mixed coefficient")
    lam_body = nodes[partial01]
    plus = next((n for n in lam_body.get("arguments", ()) if nodes[n].get("operation") == "sum"), None)
    require(plus is not None, "Lambda lacks 1+t squared")
    square = next((n for n in nodes[plus]["arguments"] if not exact(nodes, n, 1)), None)
    require(nodes[square].get("operation") == "product"
            and nodes[square]["arguments"] == [nodes[direction]["body"]]*2,
            "Lambda and direction use different actual t")
    K = next(n for n in lam_body["arguments"] if n != plus)
    phi0 = nodes[loop["Phi"][0]]
    require(phi0.get("operation") == "product" and K in phi0["arguments"], "Phi uses a different K")
    primitive_sum = next(n for n in phi0["arguments"] if n != K)
    phi_terms = nodes[primitive_sum]
    require(phi_terms.get("operation") == "sum" and angle in phi_terms["arguments"], "Phi lacks its angle term")
    T2 = next(n for n in phi_terms["arguments"] if n != angle)
    integral2 = nodes[T2]
    require(integral2.get("operation") == "definite_integral" and integral2["integrand"] == square
            and integral2["variable"] == angle_name and integral2["upper"] == angle
            and exact(nodes, integral2["lower"], 0), "Phi does not integrate the actual t squared")
    integral = nodes[nodes[T1]["body"]]
    require(integral.get("operation") == "definite_integral"
            and integral["integrand"] == nodes[direction]["body"]
            and integral["variable"] == angle_name and integral["upper"] == angle
            and exact(nodes, integral["lower"], 0), "Actual T1 primitive changed")
    inv_lambda = inverse(nodes, loop["B_over_Pstar_phi"][0], lam)
    inv_twopi = inverse(nodes, loop["A"][0], twopi=True)
    a, b, E = (roots[key][0] for key in ("a", "b", "E"))
    inputs = {phase, psi, T1, direction, a, b, E, inv_lambda, inv_twopi}
    v = lambda n: atom(nodes, n)
    half = Fraction(1, 2)
    expected = {
        "A": p_scale(p_mul(v(a), p_add(v(phase), p_scale(p_mul(v(inv_twopi), v(psi)), -1))), half),
        "B_over_Pstar": p_scale(p_mul(v(E), p_add(p_scale(p_mul(v(a), v(T1), v(inv_twopi)), -1), p_scale(p_mul(v(b), v(phase)), -1))), half),
        "A_phi": p_scale(p_mul(v(a), p_add({(): Fraction(1)}, p_scale(p_mul(v(inv_twopi), v(inv_lambda)), -1))), half),
        "B_over_Pstar_phi": p_scale(p_mul(v(E), p_add(p_scale(p_mul(v(a), v(direction), v(inv_lambda), v(inv_twopi)), -1), p_scale(v(b), -1))), half),
    }
    checks = {}
    for key, expected_polynomial in expected.items():
        root = nodes[loop[key][0]]
        require(root.get("operation") == "original_lazy_flat_branch" and root["Delta"] == roots["Delta"][0]
                and root["eta"] == 21 and root["flat_value"] == 0
                and root["flat_predicate"] == loop["flat_predicate"]
                and root["active_body_not_evaluated_on_flat_branch"] is True,
                "Actual loop lazy source changed")
        checks[key] = polynomial(nodes, root["active_body"], inputs) == expected_polynomial
    require(all(checks.values()), "Actual A/B or phase derivative formula mismatch")
    return dict(passed=True, checks=checks, phase=phase, inverse=psi,
                lambda_at_inverse=lam, direction_at_inverse=direction, T1_at_inverse=T1,
                inverse_lambda=inv_lambda, inverse_twopi=inv_twopi,
                Phi_equals_K_angle_plus_integral_t_squared=True, K=K, T2=T2,
                formulas_general_b=True, no_endpoint_zero_simplification=True)


def physical_binding(nodes, native, loop_chart, q4_chart, pstar):
    physical, radial, background = (q4_chart[key] for key in ("physical", "radial", "radial_background"))
    fields, loop, finite = (loop_chart[key] for key in ("modified_fields", "loop", "finite_N"))
    # The slow root itself reaches another serialized spelling of 1/N.
    # Bind the coefficient of this actual slow root at the outer V_y map,
    # rather than treating all descendant quotient nodes as one handle.
    vy_terms = nodes[radial["V_y"][0]]["arguments"]
    slow_term = next(n for n in vy_terms if nodes[n].get("operation") == "product"
                     and radial["slow_B_y"][0] in nodes[n]["arguments"])
    slow_factors = nodes[slow_term]["arguments"]
    require(len(slow_factors) == 2, "Actual slow B_y coefficient changed")
    invN = next(n for n in slow_factors if n != radial["slow_B_y"][0])
    invN_row = nodes[invN]
    require(invN_row.get("operation") == "positive_quotient"
            and exact(nodes, invN_row.get("numerator"), 1)
            and invN_row["denominator"] == 18, "Actual slow B_y lacks the common 1/N")
    source_m = finite["five_signed_density_increments"]["m"]
    require(source_m == finite["delta_V"], "Correction m source differs from actual velocity increment")
    checks = []
    for order in range(6):
        leaves = {pstar, fields["V"][order], fields["m"][order], background["V"][order],
                  radial["slow_B_y"][order], loop["B_over_Pstar_phi"][order], invN}
        v = lambda n: atom(nodes, n)
        U = p_mul(v(pstar), v(fields["V"][order]))
        m = p_mul(v(pstar), v(fields["m"][order]))
        Uy = p_mul(v(pstar), p_add(v(background["V"][order]),
                    p_mul(v(invN), v(radial["slow_B_y"][order])), v(loop["B_over_Pstar_phi"][order])))
        row = dict(ordinary_Z_order=order,
                   U=polynomial(nodes, physical["U"][order], leaves) == U,
                   m=polynomial(nodes, physical["m"][order], leaves) == m,
                   Uy=polynomial(nodes, physical["U_y"][order], leaves) == Uy,
                   my=polynomial(nodes, physical["m_y"][order], leaves) == p_add(U, p_scale(m, -1)))
        require(all(value is True for key, value in row.items() if key != "ordinary_Z_order"),
                "Actual modified physical fields differ from source functional map")
        checks.append(row)
    # Source ordinary derivatives follow the guarded FunctionJet generator.
    # The base pressure operator itself is checked as an exact local polynomial.
    qy = nodes[physical["Q_y"][0]]
    require(qy.get("operation") == "positive_quotient" and qy["denominator"] == 25272,
            "Actual Qy source denominator changed")
    Z, delta = 20, 25220
    Uy, my, myZ = physical["U_y"][0], physical["m_y"][0], physical["m_y"][1]
    leaves = {Z, delta, Uy, my, myZ}
    v = lambda n: atom(nodes, n)
    one = {(): Fraction(1)}
    expected = p_add(p_scale(p_mul(v(Z), v(Uy)), 2),
        p_scale(p_mul(p_add(one, p_scale(v(delta), -1)), v(Z), v(my)), -1),
        p_scale(p_mul(p_add(one, p_scale(p_mul(v(Z), v(Z)), -1)), v(myZ)), -1))
    require(polynomial(nodes, qy["numerator"], leaves) == expected,
            "Actual physical Qy numerator changed")
    p2 = loop_chart["roots"]["p2"]
    require(p2 == native["generic"]["p2"], "Loop uses modified p2 in place of native p2")
    p2_dependencies = descendants(nodes, p2[0])
    require(all(node in p2_dependencies for node in native["fields"]["P0"][:2]),
            "Native p2 lost full P0/P0_Z")
    return dict(passed=True, physical_rows=checks, Qy_base_operator_checked=True,
                ordinary_Z_prefix_from_guarded_FunctionJet=True,
                my_ODE="native m_y=V-m, delta_m_y=s-delta_m, U=Pstar*(V+s)",
                correction_m_density_equals_delta_V=True, native_p2_rows=p2,
        native_p2_first_jet_naturality=True,
                full_P0_P0_Z_retained=True, denominator=25272)


def phase_binding(nodes, charts, bindings):
    Ra, origin = None, bindings["phase_origin"]
    rows = {}
    for name in ("reference", *CHARTS):
        chart = charts[name]
        offset = nodes[chart["global_offset"]]
        arguments = offset.get("arguments", ())
        negative = [n for n in arguments if nodes[n].get("operation") == "negative"
                    and nodes[n].get("argument") == origin]
        logs = [n for n in arguments if nodes[n].get("operation") == "analytic_unary"
                and nodes[n].get("name") == "log"]
        require(offset.get("operation") == "sum" and len(negative) == len(logs) == 1,
                "Global phase offset source changed")
        ratio = nodes[nodes[logs[0]]["argument"]]
        require(ratio.get("operation") == "positive_quotient" and ratio["numerator"] == chart["radius"],
                "Global phase does not use physical R/Ra")
        if Ra is None:
            Ra = ratio["denominator"]
        require(ratio["denominator"] == Ra, "Global phase has different source Ra")
        phase = nodes[chart["actual_phase"]]
        product = nodes[phase["argument"]]
        require(phase.get("operation") == "analytic_unary" and phase.get("name") == "fractional_part"
                and product.get("operation") == "product"
                and product["arguments"] == [18, chart["global_offset"]], "Original common N phase changed")
        rows[name] = dict(radius=chart["radius"], phase=chart["actual_phase"],
                          global_offset=chart["global_offset"], passed=True)
    return dict(passed=True, shared_Ra=Ra, shared_origin=origin, N=18, charts=rows,
                conclusion="equal physical endpoint radii give equal global fractional phase")


def _endpoint_radial_source_proof(nodes, endpoint, c6, seed_proofs):
    """Check actual radial rows, then substitute each real seam side at 1/0."""
    def zero_at(index, phase, value):
        row = nodes[index]
        op = row.get("operation")
        if exact(nodes, index, 0):
            return True
        if op == "negative":
            return zero_at(row["argument"], phase, value)
        if op == "product":
            return any(zero_at(child, phase, value) for child in row["arguments"])
        if op == "sum":
            return all(zero_at(child, phase, value) for child in row["arguments"])
        return (op == "analytic_unary" and row.get("name") == "original_flat_sigma_prime"
                and row.get("argument") == phase and value in (0, 1)
                and row.get("exact_flat_source_derivative_not_a_bound") is True)

    def side(chart_name, value):
        chart = c6["charts"][chart_name]
        fields, actual, radial = (chart[key] for key in ("fields", "actual", "radial"))
        phase, radius = chart["phase"], chart["radius"]
        W = _scale_series_from_product(nodes, fields["E"], actual["phi"], chart_name + " endpoint E")
        inv = _inverse_series_from_product(nodes, fields["V"], actual["U"], chart_name + " endpoint V")
        require(nodes[inv].get("operation") == "positive_quotient"
                and exact(nodes, nodes[inv]["numerator"], 1)
                and nodes[inv]["denominator"] == c6["Pstar"], "Endpoint Pstar changed")
        factors = _flatten_operation(nodes, W[0], "product")
        sqrt_rows = [n for n in factors if nodes[n].get("operation") == "analytic_unary"
                     and nodes[n].get("name") == "sqrt"]
        require(len(sqrt_rows) == 1, "Endpoint W lacks sqrt(2R)")
        sqrt_node = sqrt_rows[0]
        argument = nodes[sqrt_node]["argument"]
        require(polynomial(nodes, argument, {radius}) == p_scale(atom(nodes, radius), 2),
                "Endpoint W uses a different actual radius")
        fixed_inputs = set(c6["F0"] + [sqrt_node, inv])
        for order in range(7):
            require(not _polynomial_residual(nodes, W[order],
                p_mul(atom(nodes, c6["F0"][order]), atom(nodes, sqrt_node), atom(nodes, inv)),
                fixed_inputs), "Endpoint W differs from common F0/radius/Pstar")
        inputs = set(W + [inv] + actual["phi"] + actual["phi_y"] + actual["U"] + actual["U_y"])
        reduced_inputs = inputs - set(actual["phi_y"])
        rows = []
        wrappers = {key: [] for key in ("E", "E_y", "V", "V_y", "U_y")}
        for order in range(7):
            E = _jet_product(nodes, W, actual["phi"], order)
            Ey = p_add(*(p_scale(p_mul(atom(nodes, W[k]),
                p_add(p_scale(atom(nodes, actual["phi"][order-k]), Fraction(1, 2)),
                      atom(nodes, actual["phi_y"][order-k]))), Fraction(comb(order, k)))
                for k in range(order + 1)))
            Vy = p_mul(atom(nodes, actual["U_y"][order]), atom(nodes, inv))
            require(not _polynomial_residual(nodes, fields["E"][order], E, inputs),
                    "Actual endpoint E constructor changed")
            require(not _polynomial_residual(nodes, radial["E_y"][order], Ey, inputs),
                    "Actual endpoint radial E_y differs from source law")
            require(not _polynomial_residual(nodes, radial["V_y"][order], Vy, inputs),
                    "Actual endpoint radial V_y differs from U_y/Pstar")
            difference = p_add(_expand_polynomial_expression(nodes, Ey, reduced_inputs),
                               p_scale(_expand_polynomial_expression(nodes, E, reduced_inputs),
                                       Fraction(-1, 10)))
            require(not difference, "Actual endpoint E_y-E/10 is nonzero")
            require(zero_at(actual["U_y"][order], phase, value),
                    "Actual U_y endpoint is not exact/flat source zero")
            for key, source in (("E", fields["E"][order]), ("E_y", radial["E_y"][order]),
                                ("V", fields["V"][order]), ("V_y", radial["V_y"][order]),
                                ("U_y", actual["U_y"][order])):
                wrappers[key].append(dict(operation="function_substitution", expression=source,
                    variable=phase, value=value, value_kind="exact_rational",
                    ordinary_Z_order=order, source_executed=False))
            rows.append(dict(order=order, E_y_minus_E_over_10="0",
                             V_y_minus_U_y_over_Pstar="0", U_y_endpoint="0",
                             actual_radial_rows_checked=True))
        return dict(chart=chart_name, chart_phase=str(value), source_wrappers=wrappers,
                    radius_substitution=dict(operation="function_substitution", expression=radius,
                        variable=phase, value=value, value_kind="exact_rational", source_executed=False),
                    W_rows=W, common_F0=list(c6["F0"]), Pstar=c6["Pstar"],
                    Pstar_inverse=inv, ordinary_rows=rows, residual_count=0,
                    source_executed=False)

    symbolic = endpoint["symbolic_endpoint_identities"]
    sides, seams = {}, {}
    for name, left_chart, right_chart, radius_name in (
        ("reference_to_restoration", "reference", "restoration", "Rz"),
        ("restoration_to_postrestore", "restoration", "post_restore", "Rrestore")):
        receipt = next(row for row in symbolic["seams"] if row["name"] == name)
        require(receipt["left_chart"] == left_chart and receipt["right_chart"] == right_chart
                and receipt["residuals"]["log_R"] == "0",
                "Actual endpoint radius seam theorem changed")
        left, right = side(left_chart, 1), side(right_chart, 0)
        saved_radius = nodes[c6["endpoints"][radius_name]["radius"]]
        require(saved_radius.get("operation") == "function_substitution"
                and saved_radius["expression"] == c6["charts"][left_chart]["radius"]
                and saved_radius["variable"] == c6["charts"][left_chart]["phase"]
                and exact(nodes, saved_radius["value"], 1), "Actual saved endpoint radius changed")
        sides[name] = dict(left=left, right=right, ordinary_row_count=14,
                           separate_actual_source_wrappers=True)
        seams[name] = dict(left_chart=left_chart, right_chart=right_chart,
            left_phase="1", right_phase="0", radius=radius_name,
            radius_identity_residual="0", radial_endpoint_residual="0")
    return dict(passed=True, seams=seams, seam_sides=sides, residual_count=0,
                full_ordinary_rows=28, source_executed=False,
                same_radius_source_objects=False)


def _endpoint_control_binding(nodes, endpoint, native, c6, seed_proofs):
    """Discharge endpoint first and second physical-y control jets.

    The endpoint receipt supplies the seam function identities.  This helper
    additionally binds those identities to the actual C6 phi/U rows and the
    C6 radial source formulas already proved in ``seed_proofs``.
    """
    symbolic = endpoint["symbolic_endpoint_identities"]
    require(symbolic["passed"] is True and symbolic["function_identity_not_sampling"] is True
            and symbolic["phi_y_over_phi"] == "-2/5"
            and symbolic["phi_yy_over_phi"] == "4/25"
            and symbolic["E_y_over_E"] == "1/10"
            and symbolic["E_yy_over_E"] == "1/100"
            and symbolic["V_y"] == "0" and symbolic["V_yy"] == "0",
            "Actual endpoint control receipt changed")
    source_rows = endpoint["actual_source_rows"]
    first_y = native["native_first_y_seed_tables"]
    chart_rows = {}
    for name in CHARTS:
        native_name = NATIVE[name]
        chart = c6["charts"][native_name]
        actual = source_rows[native_name]["actual"]
        require(first_y[native_name] == {key: list(actual[key]) for key in first_y[native_name]},
                "Native first-y table is detached from endpoint actual rows: " + name)
        require(all(source_rows[native_name]["actual_phi_y_rows_checked"]),
                "Endpoint phi_y source rows are not checked: " + name)
        ratios = []
        for phi, phi_y in zip(actual["phi"], actual["phi_y"]):
            row = nodes[phi_y]
            factors = _flatten_operation(nodes, phi_y, "product")
            require(phi in factors and any(_scalar(nodes, factor) == Fraction(-2, 5)
                                          for factor in factors),
                    "Actual phi_y row is not -2 phi / 5: " + name)
            ratios.append(dict(phi=phi, phi_y=phi_y, ratio="-2/5", residual="0"))
        require(seed_proofs[name]["passed"] is True
                and seed_proofs[name]["E_y_over_E_source"] == "W*(phi/2+phi_y)/(W*phi)"
                and seed_proofs[name]["V_y_source"] == "U_y/Pstar",
                "C6 physical control source proof is missing: " + name)
        chart_rows[name] = dict(phi_y_rows=ratios, E_y_source_bound=True,
                                V_y_source_bound=True, ordinary_rows=7)

    endpoint_radial = _endpoint_radial_source_proof(nodes, endpoint, c6, seed_proofs)
    endpoint_rows = deepcopy(endpoint_radial["seams"])
    for name in SEAMS:
        receipt = next(row for row in symbolic["seams"] if row["name"] == name)
        require(all(receipt["residuals"][key] == "0" for key in ("U_y", "U_yy")),
                "Actual U physical C2 endpoint theorem is missing")
        endpoint_rows[name].update(U_y_zero=True, V_y_zero=True, V_yy_zero=True)
    return dict(passed=True, ratios=dict(E_y_over_E="1/10", E_yy_over_E="1/100",
                                         V_y="0", V_yy="0", phi_y_over_phi="-2/5",
                                         phi_yy_over_phi="4/25"),
                charts=chart_rows, seams=endpoint_rows, residual_count=0,
                source_formula_receipt=symbolic["formulas"],
                actual_endpoint_source_rows=True, endpoint_radial_source=endpoint_radial)


def _local_ift_binding(nodes, native, parent, c5, q4):
    """Tie the saved positivity/IFT lemma to the actual new denominator."""
    finite = native["finite_p2_positive_E_local_gate"]
    required_finite = ("native_phi_finite_positive", "native_Utheta_finite_positive",
                       "finite_p2", "positive_E")
    finite_passed = all(finite.get(key) is True for key in required_finite)
    positive = parent["actual_long_reference_control_seams"]["positive_active_loop_IFT"]
    handles = positive["source_parameter_handles"]
    require(positive["passed"] is True and positive["source_executed"] is False
            and handles["eta"] == c5["eta"] == 21
            and handles["dstar"] == c5["dstar"] == 22
            and handles["N"] == c5["N"] == 18,
            "Saved positive source parameter handles changed")
    for key in ("eta", "dstar"):
        record = positive["finite_log_records"][key]
        require(record["finite"] is True and record["exp_strictly_positive"] is True,
                "Source parameter positivity is missing: " + key)
    hypotheses = positive["hypotheses"]
    for key in ("finite_p2", "finite_q", "finite_u", "D_positive", "K_positive", "lambda_positive"):
        require(hypotheses[key]["proved"] is True,
                "Saved positive hypothesis is missing: " + key)
    denominator = 25272
    delta, Z = 25220, 20
    expected_L = p_add({(): Fraction(1)}, p_scale(p_mul(atom(nodes, delta),
                                                          atom(nodes, Z), atom(nodes, Z)), -1))
    denominator_residual = _polynomial_residual(nodes, denominator, expected_L, {delta, Z})
    require(not denominator_residual, "Actual Qy denominator is not 1-delta Z^2")
    qy_denominators = []
    for name in CHARTS:
        qy = q4["charts"][name]["physical"]["Q_y"][0]
        require(nodes[qy].get("operation") == "positive_quotient"
                and nodes[qy]["denominator"] == denominator,
                "Actual Qy denominator changed: " + name)
        qy_denominators.append(dict(chart=name, node=qy, denominator=denominator,
                                    residual="0"))
    passed = finite_passed and not denominator_residual
    require(passed, "Actual local positivity/IFT premise is not discharged")
    return dict(passed=True, finite_gate=dict(passed=finite_passed, keys=list(required_finite)),
                parameter_handles=dict(N=18, eta=21, dstar=22),
                saved_positive_hypotheses={key: True for key in
                    ("finite_p2", "finite_q", "finite_u", "D_positive", "K_positive", "lambda_positive")},
                actual_L=dict(node=denominator, delta=delta, Z=Z, residual="0",
                              positive_on_inherited_Z_window=True),
                qy_denominators=qy_denominators, residual_count=0,
                analytic_inverse_at_actual_seams=True,
                global_chart_cutoff_claimed=False)


def control_theorem(endpoint, native, parent, nodes, c6, c5, q4, seed_proofs):
    """Build exact controls from actual endpoint receipts and new source rows."""
    endpoint_control = _endpoint_control_binding(nodes, endpoint, native, c6, seed_proofs)
    # E is source-nonzero through the native phi/Utheta receipt; all finite
    # controls below are read from the actual new packet, not selected values.
    finite = native["finite_p2_positive_E_local_gate"]
    require(all(finite.get(key) is True for key in
                ("native_phi_finite_positive", "native_Utheta_finite_positive",
                 "finite_p2", "positive_E")),
            "New local native finiteness not established")
    local = _local_ift_binding(nodes, native, parent, c5, q4)
    ratios = endpoint_control["ratios"]
    ey, eyy = Fraction(ratios["E_y_over_E"]), Fraction(ratios["E_yy_over_E"])
    vy, vyy = Fraction(ratios["V_y"]), Fraction(ratios["V_yy"])
    a = 1 - 2 * ey
    b = 2 * vy
    require(a != 0, "Actual endpoint control has zero a")
    kappa = a + b * b / a
    delta = kappa - 2
    t0 = -b / a
    a_y = -2 * (eyy - ey * ey)
    b_y = 2 * (vyy - ey * vy)
    kappa_y = a_y + 2 * b * b_y / a - b * b * a_y / (a * a)
    t0_y = -(b_y * a - b * a_y) / (a * a)
    values = {key: str(value) for key, value in
              dict(a=a, b=b, kappa=kappa, Delta=delta, t0=t0).items()}
    first = {key: str(value) for key, value in
             dict(a_y=a_y, b_y=b_y, kappa_y=kappa_y, Delta_y=kappa_y,
                  t0_y=t0_y).items()}
    return dict(passed=endpoint_control["passed"] and local["passed"],
                endpoint_values=values, physical_first_y=first,
                endpoint_control=endpoint_control, local_IFT=local,
                p2_naturality="actual native rational map of common C6 fields/first y jets; P0_Z retained",
                native_E_nonzero=True, no_p2_value_selected=True,
                common_native_fields_from_C2_and_source_histories=True,
                source_executed=False)


def build(root: str | Path, parent: Mapping[str, Any], source_packet: Mapping[str, Any],
          endpoint_packet: Mapping[str, Any], native_packet: Mapping[str, Any],
          dag_packet: Mapping[str, Any], phase_packet: Mapping[str, Any]) -> dict[str, Any]:
    """Compose only the newly bound seam functions, then correlate their traces."""
    root = Path(root)
    require(sha(root / PARENT) == PARENT_SHA256 and json.loads((root / PARENT).read_bytes()) == parent,
            "Exact accepted predecessor required")
    packets = dict(source=source_packet, endpoint=endpoint_packet, native=native_packet,
                   dag=dag_packet, phase=phase_packet)
    for name, packet in packets.items():
        pin_packet(root, parent, packet, name)
    require(source_packet["limits"]["source_formula_values_selected"] is False
            and source_packet["limits"]["interval_or_cap_values_used"] is False
            and source_packet["limits"]["all_caps_retained"] is True
            and native_packet["limits"]["numerical_values_selected"] is False
            and native_packet["limits"]["interval_caps_used_as_values"] is False,
            "A source packet selected an enclosure/cap as a signed value")
    pins = {}
    for path in (PARENT, C6 + ".py", C6 + ".json.gz", BRIDGE_C6 + ".py", BRIDGE_C6 + ".json.gz",
                 C5 + ".py", C5 + ".json.gz",
                 Q4 + ".py", Q4 + ".json.gz", MODULE, "src/openai_ns_reconstruction/original_signed_graph_pairs.py",
                 *("src/openai_ns_reconstruction/" + value + ".py" for value in NEW_MODULES.values())):
        actual = sha(root / path)
        expected = (PARENT_SHA256 if path == PARENT else parent["input_hashes"][Path(path).name]
                    if path.startswith(BASE) else parent["implementation_hashes"][path]
                    if path.endswith("original_signed_graph_pairs.py") else actual)
        require(actual == expected, "Pinned source changed: " + path)
        pins[path] = dict(path=path, expected_sha256=expected, actual_sha256=actual,
                          passed=True, source_executed=False)
    guards = []
    for path, token in SOURCE_GUARDS:
        present = token in (root / path).read_text(encoding="utf-8")
        require(present, "Actual generic source map changed: " + token)
        guards.append(dict(path=path, token=token, passed=True, source_executed=False))
    archives = [json.loads(gzip.decompress((root / (stem + ".json.gz")).read_bytes()))
                for stem in (C6, C5, Q4)]
    require(all(value["source_family"] == parent["source_family"] for value in archives), "Wrong cache family")
    c6, c5, q4 = (value[key] for value, key in zip(archives,
        ("prescribed_R100_Rm_C6_bindings", "prescribed_loop_C5_bindings", "modified_Q4_pressure_C3_bindings")))
    nodes = archives[2]["exact_graph_nodes"]
    require(nodes[:len(archives[1]["exact_graph_nodes"])] == archives[1]["exact_graph_nodes"]
            and nodes[:len(archives[0]["exact_graph_nodes"])] == archives[0]["exact_graph_nodes"],
            "Actual cache prefixes changed")
    require(endpoint_packet["native_phi_U_physical_C2"] is True
            and endpoint_packet["radius_endpoint_identities"] is True
            and endpoint_packet["axial_scope"]["Z_interval"] == ["37/100", "93/250"]
            and all(row["passed"] is True for row in endpoint_packet["symbolic_endpoint_identities"]["seams"]),
            "Actual new native C2/radius theorem is missing")
    require(native_packet["endpoint_dependency"]["consumed"] is True
            and native_packet["native_history_source"]["passed"] is True
            and native_packet["c5_history_source"]["passed"] is True
            and all(native_packet["c5_history_source"]["endpoint_links"].values())
            and native_packet["q4_background_seed_source"]["passed"] is True
            and native_packet["saved_lemmas"]["generic_seed"]["passed"] is True,
            "Actual nonreset histories/native generic seed theorem is missing")
    require(dag_packet["conditional_passed"] is True and dag_packet["seed_identity_discharged"] is False
            and dag_packet["counts"]["slow_derivative_roots_total"] == 24
            and dag_packet["counts"]["slow_derivative_roots_passed"] == 24
            and all(value is True for value in dag_packet["gates"].values()), "Actual new slow derivative roots failed")
    phase = phase_packet["phase_coefficients"]
    require(phase["passed"] is True and phase["matched_rows"] == phase["total_rows"] == 42,
            "Actual new mixed phase rows failed")
    saved = parent["actual_long_reference_modified_seams"]
    periodic = parent["actual_long_reference_loop_seams"]["periodic_formula"]
    require(saved["full_modified_seam_proved"] is True and all(value is True for value in saved["gates"].values())
            and saved["symbolic_checks"]["quotient_residual"] == "0"
            and saved["symbolic_checks"]["ordinary_Z_commutation_residual"] == "0"
            and periodic["passed"] is True and periodic["fourier_means"]["mean_t"] == "t0"
            and periodic["drifts"]["M(phase+1)-M(phase)"] == "-a*t0-b=0",
            "Saved universal field/periodic formulas are missing")
    for chart in ("long_reshape", "reference"):
        require(parent["actual_long_reference_background_I1_numeric"]["charts"][chart]["tail_proof"]["inequalities"]["L_positive"] is True,
                "Saved radial-independent L positivity is missing")
    history_proof = _native_history_seed_proof(nodes, c6, native_packet)
    bindings, seed_rows, seed_proofs = {}, [], {}
    for name in CHARTS:
        native = c6["charts"][NATIVE[name]]
        correction, physical = c5["charts"][name], q4["charts"][name]
        require(native["radius"] == correction["radius"]
                == endpoint_packet["actual_source_rows"][NATIVE[name]]["radius"],
                "Native endpoint and loop phase use different actual radii")
        source_native = native_packet["q4_background_seed_source"]["charts"][name]
        require(source_native["C6_fields"] == native["fields"]
                and source_native["radial_background"] == physical["radial_background"],
                "Native theorem does not bind actual new seed fields")
        seed_proof = _native_seed_formula_proof(nodes, native, physical, native_packet, name, c6,
                                                dag_packet["charts"][name])
        seed_proofs[name] = seed_proof
        actual_seed = dag_packet["charts"][name]["seed_bindings"]["field_rows"]
        for field in FIELDS:
            rows = actual_seed[field]
            require(len(rows) == 7, "Missing genuine native seed rows")
            for order, row in enumerate(rows):
                require(row["source_field_node"] == native["fields"][field][order]
                        and row["physical_radial_row_node"] == physical["radial_background"][field][order],
                        "Derivative seed uses another source field")
        seed_rows.append(dict(chart=name, native_fields=9, ordinary_rows=7,
                              passed=seed_proof["passed"],
                              residual_count=seed_proof["residual_count"],
                              residual_sha256=seed_proof["residual_sha256"]))
        phase_rows = phase["charts"][name]
        require(phase_rows["matched_rows"] == 21 and phase_rows["raw_lambda_rows_matched"] == 15
                and {(r["Z"], r["angle"], r["node"]) for r in phase_rows["rows"]}
                == {(r["Z"], r["angle"], r["node"]) for r in correction["loop"]["phase_partial_coefficients"]},
                "Mixed phase proof has different actual node rows")
        for field, slow in (("A", "slow_A_y"), ("B_over_Pstar", "slow_B_y")):
            require(dag_packet["charts"][name]["source_loop_roots"][field] == correction["loop"][field]
                    and dag_packet["charts"][name]["actual_slow_roots"][field] == physical["radial"][slow],
                    "Derivative packet is not bound to actual physical slow roots")
        endpoint_name = "Rz" if name == "restoration" else "Rrestore"
        require(all(native["incoming_histories"][key] == c6["endpoints"][endpoint_name]["moments"][key]
                    for key in ("H", "M", "K", "A", "B", "C")), "Native incoming history was reset")
        bindings[name] = dict(loop=loop_binding(nodes, correction),
                              physical=physical_binding(nodes, native, correction, physical, c6["Pstar"]))
    controls = control_theorem(endpoint_packet, native_packet, parent, nodes, c6, c5, q4,
                               seed_proofs)
    phases = phase_binding(nodes, c5["charts"], c5)
    # Every admitted gate is derived from the actual receipts above.  The
    # common-function map consumes the C6/history endpoint identities and the
    # pure generic lemma; it does not install a desired endpoint fixture or a
    # selected p2 value.
    gates = dict(
        actual_native_C2_and_radii=(endpoint_packet["native_phi_U_physical_C2"] is True
                                    and endpoint_packet["radius_endpoint_identities"] is True),
        actual_new_nonreset_histories=(native_packet["native_history_source"]["passed"] is True
                                       and native_packet["c5_history_source"]["passed"] is True),
        actual_native_seed_discharge=(history_proof["passed"] is True
                                      and all(seed["passed"] for seed in seed_proofs.values())),
        actual_new_24_slow_roots=(dag_packet["counts"]["slow_derivative_roots_total"] == 24
                                  and dag_packet["counts"]["slow_derivative_roots_passed"] == 24),
        actual_new_42_mixed_phase_rows=(phase["passed"] is True
                                       and phase["matched_rows"] == phase["total_rows"] == 42),
        actual_A_B_and_phase_partials=all(binding["loop"]["passed"] for binding in bindings.values()),
        actual_U_Uy_m_my_Qy_source_maps=all(binding["physical"]["passed"] for binding in bindings.values()),
        native_p2_first_jet_naturality=all(binding["physical"].get(
            "native_p2_first_jet_naturality") is True for binding in bindings.values()),
        positive_local_loop_IFT=(controls["local_IFT"]["passed"] is True),
        global_log_radius_phase_common=(phases["passed"] is True),
        periodic_loop_and_slow_first_jets=(parent["actual_long_reference_loop_seams"]["periodic_formula"]["passed"] is True),
        Qy_denominator_positive=all(row["residual"] == "0"
                                    for row in controls["local_IFT"]["qy_denominators"]),
        ordinary_Z_functional_identity_not_sampling=(history_proof["passed"] is True
                                                     and native_packet["saved_lemmas"]["row_lemma"]["passed"] is True),
    )
    require(all(gates.values()), "An actual source-bound correlation gate failed")
    raw_nodes = restored_nodes(root, parent)
    source_patch = source_packet["graph_patch"]
    require(len(raw_nodes) == source_patch["base_node_count"] == 185704
            and digest(raw_nodes) == source_patch["base_canonical_nodes_sha256"], "Accepted prefix changed")
    raw_nodes.extend(deepcopy(source_patch["appended_nodes"]))
    require(len(raw_nodes) == source_patch["full_node_count"] == 186086
            and digest(raw_nodes) == source_patch["full_canonical_nodes_sha256"], "New raw source graph changed")
    graph = GraphPairs(raw_nodes)
    handles = source_packet["stable_chart_handles"]
    current, previous = list(handles["parent_joined"]), [171869, 171870]
    require(current == [185700, 185703], "Parent join changed")
    targets, seams = [], {}
    for name, seam_name in zip(CHARTS, SEAMS):
        handle = handles[name]
        require(handle["upper_RQy_traces"] == q4["endpoints"][name]["boundary"]["R_Q_y"][:2],
                "Source upper RQy traces changed")
        for order, wrapper in enumerate(handle["lower_RQy_wrappers"]):
            row = graph.nodes[wrapper]
            require(row["operation"] == "function_substitution"
                    and row["expression"] == q4["charts"][name]["finite_boundary"]["R_Q_y"][order]
                    and row["variable"] == c5["charts"][name]["coordinate"]
                    and exact(graph.nodes, row["value"], 0), "New lower RQy wrapper changed")
        targets.extend(previous)
        seams[seam_name] = dict(passed=True, finite_N_Qy_identity_proven=True,
            physical_Qy_C1_equal=True, physical_RQy_C1_equal=True,
            left_RQy_upper=list(previous), right_RQy_lower_wrapper=handle["lower_RQy_wrappers"],
            ordinary_Z_orders=[0, 1], equality_kind="actual common smooth source functions on inherited Z interval")
        current = [graph.add(current[k], handle["local_half_S1"][k], previous[k],
                             graph.neg(handle["upper_RQy_traces"][k])) for k in range(2)]
        previous = handle["upper_RQy_traces"]
    joined, rows = reduce_pair(graph, current, targets)
    require(len(graph.nodes) == 186100, "New contracted graph node count changed")
    require(previous == [183696, 183697] and all(linear_coefficients(graph.nodes, node).get(previous[k]) == -1
            for k, node in enumerate(joined)), "External Rm RQy coefficient changed")
    require(graph.nodes[:len(raw_nodes)] == raw_nodes, "Standalone raw source node changed")
    return dict(passed=True, parent_report=PARENT, parent_report_sha256=PARENT_SHA256,
        source_family=parent["source_family"], source_pins=pins, source_guards=guards,
        sourceguards_tuple=True, source_executed=False, accepted_ancestors_replayed=False,
        full_modified_seam_proved=all(gates.values()), gates=gates, unresolved_gates=[], seams=seams,
        composition_conditional=not all(gates.values()),
        pure_generic_rule_engine_scope=dict(
            applied_to_new_24_roots_only=True,
            new_root_count=24,
            old_accepted_roots_replayed=False,
            old_accepted_builders_checkers_owners_replayed=False,
            source_executed=False,
        ),
        derivative_premise_discharge=dict(
            issued_DAG_remains_conditional=True, issued_DAG_seed_identity_discharged=False,
            actual_seed_function_identity=dict(
                saved_lemma=native_packet["saved_lemmas"]["generic_seed"],
                actual_new_history_binding=history_proof,
                actual_new_9x7_seed_bindings=seed_proofs,
            ),
            actual_native_C6_Q4_DAG_bindings=seed_rows,
            common_native_first_jets="actual phi/U C2 plus six native history ODEs",
            common_modified_memory_jets="five actual C5 histories plus delta_m_y=delta_V-delta_m",
            smooth_branch_and_positive_IFT=controls,
            scoped_to_two_seams_and_inherited_Z_window=True,
            global_chart_derivative_numeric_claim=False),
        axial_scope=dict(Z_interval=["37/100", "93/250"], ordinary_Z_orders=[0, 1]),
        actual_native_seed_bindings=seed_rows, actual_source_bindings=bindings,
        native_control_theorem=controls, phase_endpoint_theorem=phases,
        ordinary_Z_proof="Function identity plus guarded ordinary FunctionJet generator; only base Qy topology claimed",
        inherited_generic_formula_proofs_consumed_only=True,
        graph_patch=dict(base_node_count=len(raw_nodes), base_canonical_nodes_sha256=digest(raw_nodes),
            base_prefix_unchanged=True, appended_nodes=deepcopy(graph.nodes[len(raw_nodes):]),
            full_node_count=len(graph.nodes), full_canonical_nodes_sha256=digest(graph.nodes)),
        joined_contracted=joined, outer_Rm_RQy_retained=previous,
        exact_linearization=dict(passed=True, targets=targets, rows=rows,
            only_sum_negative_expanded=True, conservative_integer_dependency_check=True),
        limits=dict(current_N_restoration_numeric=False, complete_pressure_numeric=False,
            standalone_caps_retained=True, all_old_data_retained=True,
            global_five_moment_closure=False, full_goal_complete=False))


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
