"""Strict source bridge for the three defining-function links of heat C4.

This module closes three producer-level links which were intentionally left
conditional by ``compliant_heat_pressure_source_bridge``:

* the two incoming amplitudes are evaluations of the same
  ``SharedOuterBuffer.power`` method;
* the exact inverse radius and Gamma argument are tied to the selected heat
  radius, while every finite ``S`` box remains an enclosure;
* the C4 Gamma derivative and infinite-tail providers enclose derivatives of
  the positive Gamma expectation used by the exact heat component.

The retained pressure history (``P0``, ``Mp`` and ``Prv``) is deliberately
outside this module.  In particular, this file never promotes the conditional
absolute-pressure candidate to an admitted pressure transfer.
"""

import ast
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_angular_high_jets import integrated_gamma_tails
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import (
    gamma_deficit_mixed,
    positive_moment_derivative,
)
from lei_ren_part1_paper_compliant_angular_high_jets_check import (
    identities as admitted_angular_high_jets_identities,
)
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE = Path(__file__).parent
PREFIX = "lei_ren_part1_paper_"


def _input_hashes(names):
    """Hash every source file whose AST/runtime binding is reported below."""
    return {
        name: hashlib.sha256((HERE / name).read_bytes()).hexdigest()
        for name in names
    }


def _tree(stem):
    return ast.parse((HERE / (PREFIX + stem + ".py")).read_text(encoding="utf8"))


def _class_method(stem, class_name, method_name):
    tree = _tree(stem)
    cls = next(
        (
            node
            for node in tree.body
            if isinstance(node, ast.ClassDef) and node.name == class_name
        ),
        None,
    )
    if cls is None:
        raise ValueError(f"Missing class source: {stem}.{class_name}")
    method = next(
        (
            node
            for node in cls.body
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
            and node.name == method_name
        ),
        None,
    )
    if method is None:
        raise ValueError(f"Missing method source: {stem}.{class_name}.{method_name}")
    return method


def _assignment_nodes(stem, method_name, target):
    method = next(
        node
        for node in ast.walk(_tree(stem))
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == method_name
    )
    return [
        node
        for node in ast.walk(method)
        if isinstance(node, ast.Assign)
        and any(ast.unparse(item) == target for item in node.targets)
    ]


def _assignment_expr(stem, method_name, target, *, index=0):
    nodes = _assignment_nodes(stem, method_name, target)
    if len(nodes) <= index:
        raise ValueError(f"Missing source assignment: {stem}.{method_name}:{target}")
    return ast.unparse(nodes[index].value), nodes[index].value


def _call_keyword_expr(stem, method_name, callee, keyword, *, index=0):
    """Return one keyword value from a concrete call in a method body."""
    method = next(
        node
        for node in ast.walk(_tree(stem))
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == method_name
    )
    values = [
        item.value
        for node in ast.walk(method)
        if isinstance(node, ast.Call) and ast.unparse(node.func) == callee
        for item in node.keywords
        if item.arg == keyword
    ]
    if len(values) <= index:
        raise ValueError(f"Missing call keyword: {stem}.{method_name}:{callee}.{keyword}")
    return ast.unparse(values[index]), values[index]


def _require_expr(stem, method_name, target, expected, *, index=0):
    actual_text, actual = _assignment_expr(stem, method_name, target, index=index)
    expected_node = ast.parse(expected, mode="eval").body
    if ast.dump(actual, include_attributes=False) != ast.dump(
        expected_node, include_attributes=False
    ):
        raise ValueError(
            f"Strict source binding changed: {stem}.{method_name}:{target}: "
            f"{actual_text} != {expected}"
        )
    return actual_text


def _require_import(stem, symbol, module):
    tree = _tree(stem)
    found = []
    for node in tree.body:
        if isinstance(node, ast.ImportFrom) and node.module == module:
            found.extend(alias.name for alias in node.names)
    if symbol not in found:
        raise ValueError(f"Strict import binding missing: {stem}: {module}.{symbol}")
    return f"{module}.{symbol}"


def _require_constructor(stem, method_name, target, constructor):
    text, value = _assignment_expr(stem, method_name, target)
    if not isinstance(value, ast.Call) or ast.unparse(value.func) != constructor:
        raise ValueError(
            f"Constructor origin changed: {stem}.{method_name}:{target}: {text}"
        )
    return text


def _common_buffer_callable_binding():
    """Trace both producer objects to one concrete class method in the AST."""
    imports = {
        "high_to_amplitude": _require_import(
            "compliant_axial_high_jets",
            "CompliantAxialAmplitude",
            "lei_ren_part1_paper_compliant_axial_amplitude_selection",
        ),
        "amplitude_to_pulse": _require_import(
            "compliant_axial_amplitude_selection",
            "SharedOuterPulseMap",
            "lei_ren_part1_paper_compliant_outer_pulse_map",
        ),
        "pulse_to_buffer": _require_import(
            "compliant_outer_pulse_map",
            "SharedOuterBuffer",
            "lei_ren_part1_paper_compliant_outer_buffer",
        ),
        "angular_to_buffer": _require_import(
            "compliant_outer_angular_candidate",
            "SharedOuterBuffer",
            "lei_ren_part1_paper_compliant_outer_buffer",
        ),
    }

    constructors = {
        "axial_base": _require_constructor(
            "compliant_axial_high_jets", "__init__", "self.base", "CompliantAxialAmplitude"
        ),
        "amplitude_pulse": _require_constructor(
            "compliant_axial_amplitude_selection",
            "__init__",
            "self.pulse",
            "SharedOuterPulseMap",
        ),
        "pulse_buffer": _require_constructor(
            "compliant_outer_pulse_map", "__init__", "self.buffer", "SharedOuterBuffer"
        ),
        "angular_buffer": _require_constructor(
            "compliant_outer_angular_candidate",
            "__init__",
            "self.buffer",
            "SharedOuterBuffer",
        ),
    }

    # These are the two production evaluations.  Their arguments differ only
    # in enclosure resolution and literal-vs-local spelling of Z/phase.
    axial_call = _require_expr(
        "compliant_axial_high_jets",
        "_incoming_constants",
        "p0",
        "buffer.power(0, 1, cells=128)",
    )
    angular_call = _require_expr(
        "compliant_outer_angular_candidate",
        "__init__",
        "inlet",
        "self.buffer.power('0', 1, cells)",
    )
    buffer_origin = _require_expr(
        "compliant_axial_high_jets",
        "_incoming_constants",
        "buffer",
        "self.base.pulse.buffer",
    )
    _require_expr(
        "compliant_axial_high_jets",
        "_incoming_constants",
        "U",
        "box(p0['Utheta_over_Pstar'][0])",
    )
    _require_expr(
        "compliant_outer_angular_candidate",
        "__init__",
        "self.loguRp0",
        "c.ln(inlet['Utheta_over_Pstar'][0])",
    )

    power = _class_method("compliant_outer_buffer", "SharedOuterBuffer", "power")
    if not power.args.args or power.args.args[0].arg != "self":
        raise ValueError("SharedOuterBuffer.power is not an instance method")
    power_source = (HERE / (PREFIX + "compliant_outer_buffer.py")).resolve()

    # Resolve the concrete class attribute without recomputing the expensive
    # production packets.  AST imports/constructors/calls above establish the
    # object path; this identity check establishes the callable object.
    from lei_ren_part1_paper_compliant_outer_buffer import SharedOuterBuffer

    method = SharedOuterBuffer.__dict__["power"]
    if method is not SharedOuterBuffer.power:
        raise ValueError("SharedOuterBuffer.power runtime descriptor changed")
    if Path(method.__code__.co_filename).resolve() != power_source:
        raise ValueError("Runtime SharedOuterBuffer.power source differs from AST source")

    return {
        "verified": True,
        "imports": imports,
        "constructor_origins": constructors,
        "production_evaluations": {
            "axial_high_jets_p0": axial_call,
            "axial_high_jets_U": "box(p0['Utheta_over_Pstar'][0])",
            "outer_angular_inlet": angular_call,
            "outer_angular_loguRp0": "c.ln(inlet['Utheta_over_Pstar'][0])",
            "axial_buffer_origin": buffer_origin,
        },
        "resolved_callable": {
            "module": "lei_ren_part1_paper_compliant_outer_buffer",
            "class": "SharedOuterBuffer",
            "method": "power",
            "source": power_source.name,
            "same_unbound_function_object": method is SharedOuterBuffer.__dict__["power"],
        },
        "runtime_callable_identity_recheck": True,
        "proof_does_not_require_interval_box_equality_or_packet_recomputation": True,
    }


def _c4_constants_data_path():
    """Bind the high-jet U constant into the C4 scale path."""
    imports = {
        "fifth_to_high": _require_import(
            "compliant_fifth_axial_jets",
            "CompliantAxialHighJets",
            "lei_ren_part1_paper_compliant_axial_high_jets",
        ),
        "flatten_to_inlet": _require_import(
            "compliant_flatten_mixed_C4",
            "CompliantPowerInletC4",
            "lei_ren_part1_paper_compliant_power_inlet_C4",
        ),
        "power_to_flatten": _require_import(
            "compliant_power_angular_C4",
            "CompliantFlattenMixedC4",
            "lei_ren_part1_paper_compliant_flatten_mixed_C4",
        ),
        "steep_to_power": _require_import(
            "compliant_steep_waiting_C4",
            "CompliantPowerAngularC4",
            "lei_ren_part1_paper_compliant_power_angular_C4",
        ),
        "collar_to_steep": _require_import(
            "compliant_collar_Gamma_C4",
            "CompliantSteepWaitingC4",
            "lei_ren_part1_paper_compliant_steep_waiting_C4",
        ),
    }

    constructors = {
        "fifth_high": _require_constructor(
            "compliant_fifth_axial_jets", "__init__", "self.fourth", "CompliantAxialHighJets"
        ),
        "flatten_inlet": _require_constructor(
            "compliant_flatten_mixed_C4", "__init__", "self.inlet", "CompliantPowerInletC4"
        ),
        "power_flatten": _require_constructor(
            "compliant_power_angular_C4", "__init__", "self.flatten", "CompliantFlattenMixedC4"
        ),
        "steep_power": _require_constructor(
            "compliant_steep_waiting_C4", "__init__", "self.outer", "CompliantPowerAngularC4"
        ),
        "collar_steep": _require_constructor(
            "compliant_collar_Gamma_C4", "__init__", "self.steep", "CompliantSteepWaitingC4"
        ),
    }

    # The fifth-order receipt is produced from the high-jet object and emits
    # its Z-independent constants under the exact key consumed by PowerInlet.
    _require_expr("compliant_fifth_axial_jets", "select", "k", "self.fourth.constants")
    _, incoming_value = _call_keyword_expr(
        "compliant_fifth_axial_jets", "select", "result.update", "incoming"
    )
    if not isinstance(incoming_value, ast.Call) or ast.unparse(incoming_value.func) != "dict":
        raise ValueError("Fifth axial receipt incoming payload is not a dict call")
    keywords = {item.arg: item.value for item in incoming_value.keywords if item.arg}
    if "Z_independent_constant_definitions" not in keywords or not isinstance(
        keywords["Z_independent_constant_definitions"], ast.Name
    ) or keywords["Z_independent_constant_definitions"].id != "k":
        raise ValueError("Fifth axial receipt no longer emits high-jet constants")
    _require_expr(
        "compliant_power_inlet_C4",
        "__init__",
        "self.constants",
        "{k: read_interval(c, v) for k, v in fifth['whole_Z']['selected']['incoming']['Z_independent_constant_definitions'].items() if isinstance(v, dict) and 'lower' in v}",
    )
    _require_expr(
        "compliant_flatten_mixed_C4", "__init__", "self.U", "self.inlet.constants['U']"
    )
    _require_expr(
        "compliant_flatten_mixed_C4",
        "__init__",
        "self.logEv2_parts",
        "dict(inlet_log=2*c.ln(self.U), inverse_mu_term=-13/self.mu, finite_offset=c.mpf(-26))",
    )
    _require_expr(
        "compliant_power_angular_C4", "__init__", "self.S_cap", "read_interval(c, angular['strong_inverse_radius_positive_cap'])"
    )
    _require_expr(
        "compliant_collar_Gamma_C4", "__init__", "self.Ev2", "self.outer.flatten.Ev2"
    )
    _require_expr(
        "compliant_collar_Gamma_C4",
        "__init__",
        "self.pressure_scale",
        "self.Ev2 * self.theta_base ** 2",
    )

    return {
        "verified": True,
        "imports": imports,
        "constructor_path": constructors,
        "receipt_producer": "fifth_axial_jets.select -> self.fourth.constants",
        "receipt_key": "whole_Z.selected.incoming.Z_independent_constant_definitions",
        "c4_scale_path": [
            "CompliantPowerInletC4.constants['U']",
            "CompliantFlattenMixedC4.U",
            "CompliantFlattenMixedC4.logEv2_parts.inlet_log",
            "CompliantFlattenMixedC4.Ev2 (positive cap only)",
            "CompliantPowerAngularC4.flatten.Ev2",
            "CompliantCollarGammaC4.outer.flatten.Ev2",
            "CompliantCollarGammaC4.pressure_scale=Ev2*theta_base**2",
        ],
        "caps_are_data_enclosures": True,
    }


def _radius_and_xi_binding():
    """Tie exact S and xi to the canonical heat radius and all cap routes."""
    logC, logPstar, yd, Tw, logmu, Ts, waiting, mu = s.symbols(
        "logC logPstar yd Tw logmu Ts waiting mu", positive=True
    )
    env = {
        "self.logC": logC,
        "self.params.logPstar": logPstar,
        "self.params.yd": yd,
        "self.params.Tw": Tw,
        "self.params.log_mu": logmu,
        "self.params.Ts": Ts,
        "self.angular.waiting": waiting,
        "self.params.mu": mu,
        "self.mu": mu,
        "self.heat.logRref": s.Symbol("logRref"),
        "self.heat.tail_finite": s.Symbol("tail_finite"),
    }
    _require_expr(
        "compliant_exact_heat_component",
        "__init__",
        "self.logRref",
        "c.ln(110) + 10 * (self.logC + self.params.logPstar)",
    )
    _require_expr(
        "compliant_exact_heat_component",
        "__init__",
        "self.tail_finite",
        "self.params.yd + 1 + self.params.Tw + 100 - 30 * self.params.log_mu + 2 + self.params.Ts + self.angular.waiting",
    )
    exact_log_rref = assignment(
        "compliant_exact_heat_component", "__init__", "self.logRref", env
    )
    exact_tail_finite = assignment(
        "compliant_exact_heat_component", "__init__", "self.tail_finite", env
    )
    log_rref = s.Symbol("logRref")
    tail_finite = s.Symbol("tail_finite")
    exact_log_rref_expected = s.log(110) + 10 * (logC + logPstar)
    exact_tail_expected = yd + 1 + Tw + 100 - 30 * logmu + 2 + Ts + waiting
    if s.simplify(exact_log_rref - exact_log_rref_expected) != 0:
        raise ArithmeticError("Canonical logRref assignment changed")
    if s.simplify(exact_tail_finite - exact_tail_expected) != 0:
        raise ArithmeticError("Canonical finite tail assignment changed")

    _require_expr(
        "compliant_exact_heat_component",
        "__init__",
        "self.logradius_terms",
        "dict(selected_reference=self.logRref, pulse_term=13/self.params.mu, finite_offset=self.tail_finite)",
    )
    _require_expr(
        "compliant_outer_angular_repair",
        "__init__",
        "logRtail",
        "self.heat.logRref + 13 / self.mu + self.heat.tail_finite",
    )
    _require_expr(
        "compliant_outer_angular_repair",
        "__init__",
        "self.inverse_radius_margin",
        "logRtail + self.strong_logS_cap",
    )
    _require_expr(
        "compliant_outer_angular_repair",
        "__init__",
        "self.strong_S_cap",
        "c.exp(self.strong_logS_cap)",
    )
    _require_expr(
        "compliant_exact_heat_component", "__init__", "self.Scap", "c.exp(-1000)"
    )

    # The exact selected radius is a formal expression.  Rtail is its
    # exponential, so the reciprocal identity is proved before any cap is
    # introduced.
    log_tail = log_rref + 13 / mu + tail_finite
    exact_S = s.exp(-log_tail)
    Rtail = s.exp(log_tail)
    reciprocal_residual = s.simplify(s.expand_power_exp(exact_S * Rtail - 1))
    if reciprocal_residual != 0:
        raise ArithmeticError("Exact S=1/Rtail identity failed")
    Z_symbol, t_symbol = s.symbols("Z t", real=True)
    xi_from_S = 2 * (1 - Z_symbol**2) * exact_S * s.exp(-t_symbol)
    xi_from_Rtail = 2 * (1 - Z_symbol**2) * s.exp(-t_symbol) / Rtail
    xi_residual = s.simplify(s.expand_power_exp(xi_from_S - xi_from_Rtail))
    if xi_residual != 0:
        raise ArithmeticError("Exact xi=2*(1-Z^2)*S*exp(-t) identity failed")

    _require_expr(
        "compliant_collar_Gamma_C4",
        "gamma_deficit_mixed",
        "d",
        "IntervalTaylor(c, [1-Z**2, -2*Z, -1, 0, 0, 0])",
    )
    _require_expr(
        "compliant_collar_Gamma_C4", "gamma_deficit_mixed", "base", "d*(2*decay)"
    )
    _require_expr(
        "compliant_collar_Gamma_C4",
        "gamma_deficit_mixed",
        "S",
        "c.mpf([0, endpoints(S_cap)[1]])",
    )
    _require_expr("compliant_collar_Gamma_C4", "gamma_deficit_mixed", "xi", "base*S")
    _require_expr(
        "compliant_collar_Gamma_C4", "local_Gamma", "Sc", "self.S*c.exp(-t)"
    )
    _require_expr(
        "compliant_exact_heat_component", "deficit", "xi_cap", "2*self.Scap*decay"
    )
    xi_definition_text, xi_definition_value = _call_keyword_expr(
        "compliant_exact_heat_component", "deficit", "dict", "xi_definition"
    )
    expected_xi_definition = "xi=2*(1-Z^2)*exp(-offset)/Rtail"
    if not isinstance(xi_definition_value, ast.Constant) or xi_definition_value.value != expected_xi_definition:
        raise ValueError("Canonical exact xi definition changed")
    _require_expr(
        "compliant_corrected_outer_field",
        "heat_local",
        "Sbox",
        "self.repair.strong_S_cap*self.pulse.factor(-t)",
        index=0,
    )
    _require_expr(
        "compliant_corrected_outer_field",
        "heat_local",
        "Sbox",
        "c.mpf([0, endpoints(Sbox)[1]])",
        index=1,
    )
    _require_expr(
        "compliant_corrected_outer_field",
        "collar_shape",
        "S",
        "c.mpf([0, endpoints(self.repair.strong_S_cap)[1]])",
    )
    heat_definition, heat_literal = _assignment_expr(
        "compliant_exact_heat_component", "__init__", "self.heat_definition"
    )
    expected_heat_definition = (
        "H_delta(xi)=Gamma(1+a)^-1 integral_0^infinity "
        "exp(-v)*v^a*(1+xi*v)^-a dv, a=delta/2"
    )
    if not isinstance(heat_literal, ast.Constant) or heat_literal.value != expected_heat_definition:
        raise ValueError("Canonical exact Gamma heat definition changed")

    # The exact component's guard proves its exact S lies below exp(-1000),
    # while the repair guard proves exact S < strong_S_cap.  Neither cap is
    # therefore used as the defining function.
    guard_text = None
    repair_init = _class_method("compliant_outer_angular_repair", "CompliantAngularRepair", "__init__")
    for node in ast.walk(repair_init):
        if isinstance(node, ast.If) and "inverse_radius_margin" in ast.unparse(node.test):
            guard_text = ast.unparse(node.test)
            break
    if guard_text is None:
        raise ValueError("Strong inverse-radius cap guard disappeared")
    exact_init = _class_method("compliant_exact_heat_component", "SharedExactHeatComponent", "__init__")
    exact_cap_guard = None
    for node in ast.walk(exact_init):
        if isinstance(node, ast.If) and "-self.logRref" in ast.unparse(node.test):
            exact_cap_guard = ast.unparse(node.test)
            break
    if exact_cap_guard is None:
        raise ValueError("Canonical exp(-1000) cap guard disappeared")

    return {
        "verified": True,
        "exact_logRref_assignment": str(exact_log_rref),
        "exact_tail_finite_assignment": str(exact_tail_finite),
        "exact_logRtail": "logRref + 13/mu + tail_finite",
        "exact_S_definition": "exp(-(logRref + 13/mu + tail_finite)) = 1/Rtail",
        "reciprocal_identity": True,
        "canonical_heat_definition": expected_heat_definition,
        "xi_definition": "xi=2*(1-Z^2)*S*exp(-t) = 2*(1-Z^2)*exp(-offset)/Rtail",
        "xi_source_literal": xi_definition_text,
        "xi_symbolic_identity": s.sstr(xi_residual),
        "c4_route": {
            "gamma_deficit_input": "S_cap -> S=[0,S_cap]",
            "local_gamma_current_radius": "Sc=S_box*exp(-t)",
            "cap_is_not_definition": True,
        },
        "closure_route": {
            "exact_component": "SharedExactHeatComponent.deficit uses canonical H_delta and xi_cap=2*(1-Z^2)*exp(-offset)",
            "corrected_outer_heat_local": "strong_S_cap*exp(-t), then [0, upper]",
            "corrected_outer_collar_shape": "[0, strong_S_cap upper]",
            "cap_is_not_definition": True,
        },
        "cap_guards": {
            "exact_component_exp_minus_1000_guard": exact_cap_guard,
            "repair_strong_cap_margin_guard": guard_text,
        },
        "exact_radius_bound_by_caps": True,
        "numerical_caps_are_boxes_not_exact_functions": True,
    }


def _poly_mul(left, right, y_order, z_order):
    out = {}
    for (iy, iz), lv in left.items():
        for (jy, jz), rv in right.items():
            if iy + jy <= y_order and iz + jz <= z_order:
                out[(iy + jy, iz + jz)] = out.get((iy + jy, iz + jz), mp.mpf(0)) + lv * rv
    return out


def _poly_pow(base, power, y_order, z_order):
    out = {(0, 0): mp.mpf(1)}
    for _ in range(power):
        out = _poly_mul(out, base, y_order, z_order)
    return out


def _canonical_H_derivative(a, xi, order):
    return (
        (-1) ** order
        * mp.rf(a, order)
        * mp.quad(
            lambda v: mp.exp(-v)
            * v ** (a + order)
            * (1 + xi * v) ** (-a - order),
            [0, 1, mp.inf],
        )
        / mp.gamma(1 + a)
    )


def _canonical_D_derivatives(a, xi, order):
    values = []
    H0 = mp.quad(
        lambda v: mp.exp(-v) * v**a * (1 + xi * v) ** (-a),
        [0, 1, mp.inf],
    ) / mp.gamma(1 + a)
    values.append((1 - H0) / (a * xi))
    for n in range(1, order + 1):
        Hn = _canonical_H_derivative(a, xi, n)
        values.append((-Hn / a - n * values[n - 1]) / xi)
    return values


def _canonical_bivariate_coefficients(a, S, Z, t, y_order, z_order):
    """Taylor coefficients of the actual integral, not an interval proxy."""
    xi0 = 2 * (1 - Z**2) * S * mp.exp(-t)
    d_z = {(0, 0): 1 - Z**2, (0, 1): -2 * Z, (0, 2): -1}
    exp_y = {(k, 0): (-1) ** k / mp.factorial(k) for k in range(y_order + 1)}
    xi = {}
    for (iy, iz), ev in exp_y.items():
        for (jy, jz), dv in d_z.items():
            if iy + jy <= y_order and iz + jz <= z_order:
                xi[(iy + jy, iz + jz)] = (
                    xi.get((iy + jy, iz + jz), mp.mpf(0))
                    + 2 * S * mp.exp(-t) * ev * dv
                )
    delta = dict(xi)
    delta[(0, 0)] -= xi0
    h = {}
    for n in range(y_order + z_order + 1):
        coefficient = _canonical_H_derivative(a, xi0, n) / mp.factorial(n)
        for key, value in _poly_pow(delta, n, y_order, z_order).items():
            h[key] = h.get(key, mp.mpf(0)) + coefficient * value
    h[(0, 0)] = 1 - h.get((0, 0), mp.mpf(0))
    for key in list(h):
        if key != (0, 0):
            h[key] = -h[key]
    return {key: value / (a * S) for key, value in h.items()}


def _moment_and_derivative_proof():
    """Check Gamma moments, derivative formulas, and production bounds."""
    a, xi, v = s.symbols("a xi v", positive=True)
    moment_steps = {}
    for order in range(11):
        lhs = s.gamma(1 + a + order) / s.gamma(1 + a)
        rhs = s.prod(1 + a + j for j in range(order))
        residual = s.simplify(s.expand_func(lhs) - rhs)
        if residual != 0:
            raise ArithmeticError(f"Gamma moment recurrence failed at {order}")
        moment_steps[str(order)] = True
    derivative_steps = {}
    for order in range(10):
        residual = s.simplify(
            s.diff((1 + xi * v) ** (-a), xi, order)
            - (-1) ** order * s.rf(a, order) * v**order * (1 + xi * v) ** (-a - order)
        )
        if residual != 0:
            raise ArithmeticError(f"Gamma integrand derivative failed at {order}")
        derivative_steps[str(order)] = True

    # Match both M/nextM branches exactly to the Gamma moment formulas used
    # by the production interval bounds.
    divided_M = _assignment_expr(
        "compliant_collar_Gamma_C4", "positive_moment_derivative", "M", index=0
    )[0]
    divided_next = _assignment_expr(
        "compliant_collar_Gamma_C4", "positive_moment_derivative", "nextM", index=0
    )[0]
    plain_M = _assignment_expr(
        "compliant_collar_Gamma_C4", "positive_moment_derivative", "M", index=1
    )[0]
    plain_next = _assignment_expr(
        "compliant_collar_Gamma_C4", "positive_moment_derivative", "nextM", index=1
    )[0]
    if divided_M != "pochhammer(c, 1 + a, n) * pochhammer(c, 1 + a, n + 1) / (n + 1)":
        raise ValueError("Divided Gamma M production expression changed")
    if divided_next != "pochhammer(c, 1 + a, n + 1) * pochhammer(c, 1 + a, n + 2) / (n + 2)":
        raise ValueError("Divided Gamma nextM production expression changed")
    if plain_M != "pochhammer(c, 1 + a, n - 1) * pochhammer(c, 1 + a, n)":
        raise ValueError("Plain Gamma M production expression changed")
    if plain_next != "pochhammer(c, 1 + a, n) * pochhammer(c, 1 + a, n + 1)":
        raise ValueError("Plain Gamma nextM production expression changed")

    # The tangent inequality (1+x)^(-b) >= 1-bx gives the production lower
    # bound; the upper bound is positivity plus the zero-x Gamma moment.
    bound_steps = {}
    for order in range(5):
        plain_m = s.prod(1 + a + j for j in range(max(order - 1, 0))) * s.prod(
            1 + a + j for j in range(order)
        ) if order else s.Integer(0)
        if order:
            plain_n = s.prod(1 + a + j for j in range(order)) * s.prod(
                1 + a + j for j in range(order + 1)
            )
            expected_plain_tangent = (a + order) * s.prod(
                1 + a + j for j in range(order - 1)
            ) * s.prod(1 + a + j for j in range(order + 1))
            if s.simplify(plain_n - expected_plain_tangent) != 0:
                raise ArithmeticError("Plain nextM tangent moment mismatch")
        divided_m = s.prod(1 + a + j for j in range(order)) * s.prod(
            1 + a + j for j in range(order + 1)
        ) / (order + 1)
        divided_n = s.prod(1 + a + j for j in range(order + 1)) * s.prod(
            1 + a + j for j in range(order + 2)
        ) / (order + 2)
        if order >= 0:
            tangent = (1 + a + order) * s.prod(1 + a + j for j in range(order)) * s.prod(
                1 + a + j for j in range(order + 2)
            ) / (order + 2)
            if s.simplify(divided_n - tangent) != 0:
                raise ArithmeticError("Divided nextM tangent moment mismatch")
        bound_steps[str(order)] = True

    # Finite positive fixture: evaluate the actual Gamma integral and check
    # inclusion in both production interval bound variants.
    with mp.workdps(45):
        c = MPIntervalContext()
        c.dps = 70
        aval = mp.mpf(".15")
        xval = mp.mpf(".2")
        dvals = _canonical_D_derivatives(aval, xval, 4)
        divided_checks = 0
        for order, value in enumerate(dvals):
            box = positive_moment_derivative(c, c.mpf(".15"), c.mpf(".2"), order, True)
            lo, hi = endpoints(box)
            if not lo <= value <= hi:
                raise ArithmeticError(f"Divided Gamma derivative escaped bound {order}")
            divided_checks += 1
        plain_checks = 0
        for order in range(1, 5):
            value = _canonical_H_derivative(aval, xval, order) / aval
            box = positive_moment_derivative(c, c.mpf(".15"), c.mpf(".2"), order, False)
            lo, hi = endpoints(box)
            if not lo <= value <= hi:
                raise ArithmeticError(f"Plain Gamma derivative escaped bound {order}")
            plain_checks += 1

    return {
        "gamma_moment_recurrence": moment_steps,
        "integrand_derivative_formula": derivative_steps,
        "production_M_nextM_expressions": {
            "divided_M": divided_M,
            "divided_nextM": divided_next,
            "plain_M": plain_M,
            "plain_nextM": plain_next,
        },
        "tangent_and_positivity_bounds": bound_steps,
        "divided_derivative_fixture_checks": divided_checks,
        "plain_derivative_fixture_checks": plain_checks,
        "canonical_H_is_positive_integral": True,
        "canonical_H_definition_used": "Gamma(1+a)^-1 * integral exp(-v)*v^a*(1+xi*v)^(-a) dv",
    }


def _admitted_angular_tail_receipt():
    """Validate the already-admitted independent Gamma tail fixture."""
    name = "lei_ren_part1_paper_compliant_angular_high_jets_check.json"
    receipt = json.loads((HERE / name).read_text(encoding="utf8"))
    gamma = receipt.get("independent_Gamma_fixture", {})
    if not receipt.get("all_passed") or not gamma.get("passed"):
        raise ValueError("Admitted angular Gamma tail receipt is not complete")
    for source, digest in receipt.get("input_hashes", {}).items():
        path = HERE / source
        if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError("Admitted angular Gamma source changed: " + source)
    return {
        "receipt": name,
        "all_passed": True,
        "true_Gamma_value_and_axial_tail_derivative_checks": True,
        "finite_pressure_energy_tail_errors_bounded": True,
        "source_hashes_revalidated": True,
    }


def _gamma_enclosure_proof():
    """Bind C4 composition/tail source and test its enclosure semantics."""
    admitted_angular = admitted_angular_high_jets_identities()
    admitted_tail = _admitted_angular_tail_receipt()
    required_angular = {
        "Gamma_integrand_derivative1",
        "Gamma_integrand_derivative2",
        "Gamma_integrand_derivative3",
        "Gamma_integrand_derivative4",
        "angular_leading_integral_cancellation",
        "square_pressure_deficit",
        "quadratic_axial_composition",
    }
    if not required_angular.issubset(admitted_angular):
        raise ValueError("Admitted angular Gamma source lemma is incomplete")
    _require_import(
        "compliant_collar_Gamma_C4",
        "integrated_gamma_tails",
        "lei_ren_part1_paper_compliant_angular_high_jets",
    )
    _require_expr(
        "compliant_collar_Gamma_C4",
        "gamma_deficit_mixed",
        "F",
        "sum((powers[n] * (positive_moment_derivative(c, a, xi[0], n, True) / math.factorial(n)) for n in range(6)), base * 0)",
    )
    _require_expr(
        "compliant_collar_Gamma_C4",
        "gamma_deficit_mixed",
        "G",
        "{l: sum((powers[n] * (positive_moment_derivative(c, a, xi[0], l + n) / math.factorial(n)) for n in range(6)), base * 0) for l in range(1, yorder + 1)}",
    )
    _require_expr(
        "compliant_collar_Gamma_C4",
        "gamma_deficit_mixed",
        "dx",
        "IntervalTaylor(c, [0] + list(xi.coefficients[1:]))",
    )
    _require_expr(
        "compliant_collar_Gamma_C4", "gamma_deficit_mixed", "rows", "[base * F]"
    )
    _, augmented = _assignment_expr("compliant_collar_Gamma_C4", "gamma_deficit_mixed", "row")
    # ``row += ...`` is an AugAssign, so inspect it directly.
    method = next(
        node
        for node in ast.walk(_tree("compliant_collar_Gamma_C4"))
        if isinstance(node, ast.FunctionDef) and node.name == "gamma_deficit_mixed"
    )
    aug = [node for node in ast.walk(method) if isinstance(node, ast.AugAssign) and ast.unparse(node.target) == "row"]
    if len(aug) != 1 or ast.unparse(aug[0].value) != "base ** l * G[l] * (stirling_second(k, l) * (-1) ** (k + 1) * S ** (l - 1))":
        raise ValueError("Gamma y/chain production recurrence changed")

    # Tail source: derivative coefficients are squared by the explicit
    # convolution and integrated against exact exponential tails.
    tail_method = next(
        node
        for node in ast.walk(_tree("compliant_angular_high_jets"))
        if isinstance(node, ast.FunctionDef) and node.name == "integrated_gamma_tails"
    )
    tail_assignments = [
        ast.unparse(node.value)
        for node in ast.walk(tail_method)
        if isinstance(node, ast.Assign)
        and any(
            isinstance(target, ast.Tuple)
            and [ast.unparse(element) for element in target.elts] == ["A", "errors"]
            for target in node.targets
        )
    ]
    if tail_assignments != ["gamma_leading_error_weights(c, Z, a, S_cap, order)"]:
        raise ValueError("Gamma tail leading/error source changed")
    _require_expr(
        "compliant_angular_high_jets",
        "integrated_gamma_tails",
        "full",
        "[dict(w) for w in errors]",
    )
    _require_expr(
        "compliant_angular_high_jets",
        "integrated_gamma_tails",
        "square",
        "[{} for _ in range(order + 1)]",
    )

    # Product coefficient identity used by the square convolution.
    f = s.symbols("f0:6")
    x = s.symbols("x")
    product = s.expand(sum(f[j] * x**j for j in range(6)) ** 2)
    product_proofs = {}
    for order in range(6):
        rhs = sum(f[i] * f[order - i] for i in range(order + 1))
        if s.expand(product).coeff(x, order) != rhs:
            raise ArithmeticError("Gamma tail product coefficient identity failed")
        product_proofs[str(order)] = True

    # The y chain rule is checked on the canonical positive Gamma integrand
    # itself.  For xi=B*S*exp(-y), d/dy=-xi*d/dxi and the production rows
    # use the Stirling expansion of this operator on the actual integrand.
    xi, v, aa, B, S0 = s.symbols("xi v a B S", positive=True)
    integrand = (1 + xi * v) ** (-aa)
    chain_proofs = {}
    for order in range(5):
        direct = integrand
        for _ in range(order):
            direct = -xi * s.diff(direct, xi)
        operator_sum = (
            integrand
            if order == 0
            else sum(
                s.functions.combinatorial.numbers.stirling(order, ell, kind=2)
                * (-1) ** order
                * xi**ell
                * s.diff(integrand, xi, ell)
                for ell in range(1, order + 1)
            )
        )
        if s.simplify(direct - operator_sum) != 0:
            raise ArithmeticError("Gamma y chain coefficient identity failed")
        # The deficit carries the additional minus sign and 1/(a*S).  After
        # xi=B*S, this is exactly the production factor
        # B^ell*S^(ell-1)*(-1)^(k+1)/a multiplying H^(ell).
        if order:
            deficit_direct = -direct.subs(xi, B * S0) / (aa * S0)
            production_factor = sum(
                s.functions.combinatorial.numbers.stirling(order, ell, kind=2)
                * (-1) ** (order + 1)
                * B**ell
                * S0 ** (ell - 1)
                * s.diff(integrand, xi, ell).subs(xi, B * S0)
                / aa
                for ell in range(1, order + 1)
            )
            if s.simplify(deficit_direct - production_factor) != 0:
                raise ArithmeticError("Gamma deficit chain factor identity failed")
        chain_proofs[str(order)] = True

    # Exercise both production functions on a finite positive cap.  The
    # actual integral coefficients below are generated from the explicit
    # positive Gamma integrand and checked against every requested C4 row.
    with mp.workdps(40):
        c = MPIntervalContext()
        c.dps = 65
        aval = mp.mpf(".15")
        S_exact = mp.mpf(".1")
        cap = c.mpf(".2")
        zval = mp.mpf(".3")
        tval = mp.mpf(".2")
        exact = _canonical_bivariate_coefficients(aval, S_exact, zval, tval, 4, 2)
        rows = gamma_deficit_mixed(
            c,
            c.mpf(str(zval)),
            c.mpf(str(aval)),
            cap,
            c.mpf(str(tval)),
            yorder=4,
        )
        mixed_checks = 0
        for iy in range(5):
            for iz in range(3):
                lo, hi = endpoints(rows[iy][iz])
                value = exact[(iy, iz)]
                if not lo <= value <= hi:
                    raise ArithmeticError(f"Gamma mixed enclosure escaped at {(iy, iz)}")
                mixed_checks += 1
        tails = integrated_gamma_tails(
            c,
            c.mpf(str(zval)),
            c.mpf(str(aval)),
            cap,
            order=2,
            start=3,
        )
        tail_finiteness = 0
        for name in ("theta", "pressure", "energy"):
            for index in range(3):
                lo, hi = endpoints(tails[name][index])
                if not (mp.isfinite(lo) and mp.isfinite(hi) and lo <= hi):
                    raise ArithmeticError(f"Gamma tail enclosure invalid: {name}[{index}]")
                tail_finiteness += 1

    # Exact exponential tail primitives are the final enclosure step.
    q, start, y = s.symbols("q start y", positive=True)
    tail_integral_residual = s.simplify(
        s.integrate(s.exp(-q * y), (y, start, s.oo)) - s.exp(-q * start) / q
    )
    if tail_integral_residual != 0:
        raise ArithmeticError("Gamma exponential tail primitive changed")

    return {
        "verified": True,
        "production_positive_moment_derivative_bound": True,
        "production_gamma_deficit_mixed_chain": True,
        "production_integrated_gamma_tail_convolution": True,
        "product_coefficient_checks": product_proofs,
        "explicit_gamma_y_chain_checks": chain_proofs,
        "exact_exponential_tail_primitive": True,
        "mixed_C4_fixture_checks": mixed_checks,
        "tail_fixture_coefficients_checked": tail_finiteness,
        "finite_fixture_uses_positive_Gamma_integrand": True,
        "finite_fixture_is_not_a_replacement_for_source_binding": True,
        "admitted_angular_high_jets_source_lemma": {
            "source": "lei_ren_part1_paper_compliant_angular_high_jets_check.identities",
            "checks": admitted_angular,
        },
        "admitted_angular_high_jets_tail_fixture": admitted_tail,
    }


def defining_function_bridge():
    """Return explicit evidence for the three defining-function links."""
    common = _common_buffer_callable_binding()
    constants = _c4_constants_data_path()
    radius = _radius_and_xi_binding()
    moments = _moment_and_derivative_proof()
    gamma = _gamma_enclosure_proof()
    bound_sources = [
        "lei_ren_part1_paper_compliant_axial_high_jets.py",
        "lei_ren_part1_paper_compliant_axial_amplitude_selection.py",
        "lei_ren_part1_paper_compliant_outer_pulse_map.py",
        "lei_ren_part1_paper_compliant_outer_buffer.py",
        "lei_ren_part1_paper_compliant_outer_angular_candidate.py",
        "lei_ren_part1_paper_compliant_fifth_axial_jets.py",
        "lei_ren_part1_paper_compliant_power_inlet_C4.py",
        "lei_ren_part1_paper_compliant_flatten_mixed_C4.py",
        "lei_ren_part1_paper_compliant_power_angular_C4.py",
        "lei_ren_part1_paper_compliant_steep_waiting_C4.py",
        "lei_ren_part1_paper_compliant_collar_Gamma_C4.py",
        "lei_ren_part1_paper_compliant_angular_high_jets_check.py",
        "lei_ren_part1_paper_compliant_angular_high_jets_check.json",
        "lei_ren_part1_paper_compliant_exact_heat_component.py",
        "lei_ren_part1_paper_compliant_outer_angular_repair.py",
        "lei_ren_part1_paper_compliant_corrected_outer_field.py",
        "lei_ren_part1_paper_compliant_angular_high_jets.py",
    ]
    result = {
        "common_buffer_callable_origin_verified": common["verified"],
        "c4_constants_data_path_verified": constants["verified"],
        "exact_heat_radius_and_xi_binding_verified": radius["verified"],
        "canonical_gamma_derivative_enclosure_verified": moments is not None and gamma["verified"],
        "three_defining_function_links_verified": True,
        "pressure_history_binding_verified": False,
        "complete_defining_function_history_bridge_verified": False,
        "absolute_pressure_same_source_mixed4_available": False,
        "source_history_transfer_conditional": True,
        "input_hashes": _input_hashes(bound_sources),
        "common_buffer_callable": common,
        "c4_constants_path": constants,
        "exact_radius_and_xi": radius,
        "gamma_derivative_enclosure": {**moments, **gamma},
        "remaining_gap": "Actual C4 Ptail/forward_pressure integrals to retained P0/Mp/Prv histories are intentionally left to the pressure-history bridge.",
        "scope_excludes": [
            "heat-exterior stress identity",
            "global admissible stress/cone lift",
            "positive-order temporal recursion",
            "oscillatory correction",
        ],
    }
    return result


def run():
    result = defining_function_bridge()
    print(json.dumps(result, indent=2, default=str))
    return result


if __name__ == "__main__":
    run()
