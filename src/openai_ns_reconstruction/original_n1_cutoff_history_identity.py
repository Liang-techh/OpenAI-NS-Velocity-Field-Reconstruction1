"""Source-bound first-chart n1 cutoff, own history, and signed m15 identity.

The accepted inner graph supplies the raw regular first-chart state.  This
packet binds the separate own cutoff

    eta(s) = 1 - sigma(4 * (s / (sc/4)) - 2)

to that state, then records the exact history formulas for mean, K1, Q1, V1,
pressure, and signed m15.  The core inlet and the leading N1p pressure source
remain present after eta becomes zero.  The finite numerical evaluator and
its intervals are not promoted to function values; an optional raw_binding
must be supplied by the root consumer for any numerical prefix gate.

The raw Volterra ``system['N1p']`` is the physical source.  The historical
finite normalized pressure routine represents that same source as
``raw_N1p/Pstar^2``; the physical radial pressure row restores ``Pstar^2``
once before the signed m15 weight applies ``-epsilon^2*rho``.

Only saved JSON/gzip graph packets, source ASTs, and fresh symbolic algebra
are inspected.  No accepted producer, owner, checker, recurrence, norm
calculation, or quadrature is executed.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

import sympy as sp


PARENT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_native_collar_transport.json"
)
PARENT_SHA256 = "7c05f7438a27d9caa4378cb89b9ab6d33b44138086e86292b73f2fd8e3f4bf96"

CUTOFF_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py"
)
CUTOFF_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.json.gz"
)
CUTOFF_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions_check.json"
)
CUTOFF_SOURCE_SHA256 = "b9f8b854b51b9eb0eda338503c8d25e6a1b69307da54b75bf40c1bd212ff1b29"
CUTOFF_ARCHIVE_SHA256 = "600f4c6e69cff5b95e43312c1a9da6e77ba5e204d9873eabf27d72ba09100078"
CUTOFF_CHECK_SHA256 = "688fc660eb5d7d38cb3f255c97275818de7451466cd121646803af7ba6c6ba47"

RADIAL_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_radial_functions.py"
)
RADIAL_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_radial_functions.json.gz"
)
RADIAL_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_radial_functions_check.json"
)
RADIAL_SOURCE_SHA256 = "d80c05494bf5901e09235d2757c9da1a7ff103b925e337b61774f738ad3f68b2"
RADIAL_ARCHIVE_SHA256 = "91653044b450157ecb74fd5023bb61ccf1d60a959d8ecd56fbfb869f8d50b5b3"
RADIAL_CHECK_SHA256 = "0d0ac44de54df037fbd4d170a705fb69dd73634f5aa1133e72332d6fbe6ba330"

VOLTERRA_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
VOLTERRA_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz"
)
VOLTERRA_CHECK = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions_check.json"
)
VOLTERRA_SOURCE_SHA256 = "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208"
VOLTERRA_ARCHIVE_SHA256 = "1b44be8c6cbfffb692c7a0430ecbc8a590f46001a233c99c59ca8c5e2c44a6fe"
VOLTERRA_CHECK_SHA256 = "3cfdc3ce6a7e8c2a98a31fae985c02a2acb4e28025d313ee35db4cffc6e7b47c"

CUTOFF_NUMERIC_SOURCE = "src/openai_ns_reconstruction/original_n1_cutoff_numeric.py"
CUTOFF_NUMERIC_SOURCE_SHA256 = "7a2004b378ee495c6c053d04a8394406937f3a6fca235a721138f1ce4661c731"
M15_NUMERIC_SOURCE = "src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py"
M15_NUMERIC_SOURCE_SHA256 = "004341ab4eb1f4693995bc86a748eae4c20ff2b5406fb0d784c4416963f5988b"
M15_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_m15_numeric.json"
)
M15_REPORT_SHA256 = "0615f83f6a0c67ee506ca143c1d4ba698d06a92a63e018992b608c15dc930d95"

FAMILY_KEYS = (
    "actual_five_defect_family_sha256",
    "implicit_source_sha256",
    "datum_enclosure_sha256",
)
ORDINARY_Z_ORDERS = [0, 1, 2, 3]
RAW_FIELDS = ("F1", "Uz1", "K1", "P1", "d_x_F1", "d_x_Uz1")


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load(path: Path) -> Mapping[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    _require(isinstance(value, Mapping), f"Saved mapping required: {path}")
    return value


def _load_gzip(path: Path) -> Mapping[str, Any]:
    value = json.loads(gzip.open(path, "rt", encoding="utf-8").read())
    _require(isinstance(value, Mapping), f"Saved gzip mapping required: {path}")
    return value


def _family(value: Any) -> dict[str, Any]:
    value = value if isinstance(value, Mapping) else {}
    nested = value.get("source_family")
    if isinstance(nested, Mapping):
        value = nested
    return {key: value.get(key) for key in FAMILY_KEYS}


def _pin(root: Path, relative: str, expected: str) -> dict[str, Any]:
    actual = _sha(root / relative)
    return {
        "path": relative,
        "sha256": actual,
        "expected_sha256": expected,
        "passed": actual == expected,
        "source_executed": False,
    }


def _find_function(tree: ast.AST, name: str, cls: str | None = None) -> ast.AST:
    scope: ast.AST = tree
    if cls is not None:
        classes = [node for node in ast.walk(tree)
                   if isinstance(node, ast.ClassDef) and node.name == cls]
        _require(len(classes) == 1, f"Unique class required: {cls}")
        scope = classes[0]
    functions = [node for node in ast.walk(scope)
                 if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
                 and node.name == name]
    _require(len(functions) == 1, f"Unique function required: {name}")
    return functions[0]


def _ast_guard(root: Path, relative: str, name: str,
               fragments: tuple[str, ...], cls: str | None = None) -> dict[str, Any]:
    row: dict[str, Any] = {
        "source_path": str(root / relative),
        "function": name,
        "class": cls,
        "required_fragments": {fragment: False for fragment in fragments},
        "function_AST_sha256": None,
        "source_executed": False,
        "passed": False,
    }
    try:
        source = (root / relative).read_text(encoding="utf-8")
        fn = _find_function(ast.parse(source, filename=relative), name, cls)
        segment = ast.get_source_segment(source, fn) or ""
        compact = "".join(segment.split())
        row["required_fragments"] = {
            fragment: "".join(fragment.split()) in compact
            for fragment in fragments
        }
        row["function_AST_sha256"] = hashlib.sha256(
            ast.dump(fn, include_attributes=False).encode("utf-8")
        ).hexdigest()
        row["passed"] = all(row["required_fragments"].values())
    except (OSError, SyntaxError, UnicodeError, ValueError) as exc:
        row["error"] = str(exc)
    return row


def _source_ast(root: Path) -> tuple[dict[str, Any], bool]:
    guards = {
        "cutoff_geometry": _ast_guard(
            root, CUTOFF_SOURCE, "original_cutoff_rows", (
                "q=g.quotient(s,phase,'positive original common-inner phase')",
                "arg=g.sub(g.mul(c(4),q),c(2));sigma=own.jets.complete.Gamma.sigma(g,arg)",
                "value=g.sub(g.one,sigma)",
                "exact_hb_inverse_and_inverse_squared_retained=True",
                "Z_independent_cutoff=True",
            ),
        ),
        "own_history_and_pressure": _ast_guard(
            root, CUTOFF_SOURCE, "own_cutoff_fields", (
                "history=inner.collar.integral(g,body,v,g.zero,s,actual_cutoff_own_mean=True,history_constant_after_velocity_zero_by_flat_support=True,core_inlet_not_reset=True)",
                "mean=a.div(a.add(a.scale(inlet,4),history),a.fixed(r))",
                "pressure=a.add(inner.collar.at(g,inner_owner.functions['core']['P1'],rcore,c(4)),inner.collar.integral(g,pbody,t,g.zero,s,actual_cutoff_own_pressure=True,leading_N1p_not_cut=True,nonzero_core_pressure_inlet_retained=True))",
                "recovered=radial.recovered_radial_rows(g,cut,r,eps,p['delta'])",
            ),
        ),
        "signed_m15_integrand": _ast_guard(
            root, CUTOFF_SOURCE, "signed_rho_integrands", (
                "pressure=a.add(a.scale(a.mul(f,F),2),system['N1p'])",
                "m13=a.scale(pressure,epsilon)",
                "m15=a.add(a.scale(a.mul(u,U),g.mul(g.constant(2),epsilon)),a.neg(a.scale(pressure,g.mul(eps2,r))))",
            ),
        ),
        "raw_system_physical_N1p": _ast_guard(
            root, VOLTERRA_SOURCE, "raw_system", (
                "Np=a.add(a.scale('-1/2',S),a.scale(parameters['Lambda'],a.add(a.scale(2,c1(Q[1])),a.scale(r,c1(Q[2])))))",
                "g4=g.mul(c(2),eps,x,Np.value)",
                "g6=g.sub(g.mul(c(2),eps,Nz),g.quotient(g.mul(c(4),zv,r,eps,eps,Np.value),Lv,'positive original L'))",
                "N1p=derivative.c3(Np.value)",
            ),
        ),
        "physical_pressure_radial_function": _ast_guard(
            root, RADIAL_SOURCE, "pressure_radial_function", (
                "return a.scale(a.add(a.scale(a.mul(leadingF,F1),2),N1p),epsilon)",
            ),
        ),
        "radial_recovery": _ast_guard(
            root, VOLTERRA_SOURCE, "recovery", (
                "Q=b.div(b.add(b.mul(b.add(b.fixed(1),b.neg(de)),b.mul(z,c2(U))),b.neg(b.mul(b.add(b.fixed(1),de),b.mul(z,c2(K)))),b.neg(b.mul(d,mz))),L)",
                "return dict(mean1=mean,radial_primitive_M1=a.scale(mean,R),Q1=Q,V1=b.scale(Q,R)",
                "temporal_n1_mean_weight_is_1_plus_delta=True",
                "raw_physical_units_retained=True",
            ),
        ),
        "green_history": _ast_guard(
            root, VOLTERRA_SOURCE, "green_apply", (
                "complete_history_not_reset_at_core_seam=True",
                "signed_body_not_magnitude_cap=True",
            ),
        ),
        "numeric_cutoff_recovery": _ast_guard(
            root, CUTOFF_NUMERIC_SOURCE, "_evaluate", (
                "meanp=primitive_p/4;meant=primitive_t/rho",
                "kp,kt=meanp-up,meant-ut",
                "qp,qt=recover(up,meanp,c,o.core.delta),recover(ut,meant,iv,o.core.delta_iv)",
                "pp=cp['p']+jp['p']*sig",
                "leading_N1p_retained_once=True",
            ),
        ),
        "numeric_m15_source": _ast_guard(
            root, M15_NUMERIC_SOURCE, "source_cell", (
                "pressure=(-(Np+coupling)*vt['rho']*vt['rho'])",
                "coupling=(v['S']*vt['Phi'][0].rows[0]*f*(2/(v['e']**3*v['P2'])))",
            ),
        ),
        "numeric_normalized_pressure_source": _ast_guard(
            root, M15_NUMERIC_SOURCE, "normalized_pressure_source", (
                "Original first-collar N1p/Pstar^2",
                "return -flow/(2*P2)+(qr*2+qrr*r)/(e*P2)",
            ),
        ),
        "numeric_m15_history": _ast_guard(
            root, M15_NUMERIC_SOURCE, "evaluate", (
                "terms=copy.deepcopy(inlet['m15_core_C0_C1_terms'])",
                "append('first_lower_axial_width'",
                "append('first_lower_regular_pressure_width'",
                "density_after_inverse_hb_cancellation='A(s,h)*sigma_prime(s) ds'",
                "actual_signed_core_inlet_retained_once=True",
            ),
        ),
    }
    return guards, all(row.get("passed") is True for row in guards.values())


def _raw_rows(nodes: list[Mapping[str, Any]], start: int, chart: str,
              upper: int) -> dict[str, Any]:
    rows: dict[str, Any] = {}
    passed = True
    for field_index, field in enumerate(RAW_FIELDS):
        entries = []
        for order in ORDINARY_Z_ORDERS:
            index = start + field_index * 4 + order
            node = nodes[index]
            expected_power = 1 if field.startswith("d_x_") else 2
            ok = (
                node.get("operation") == "current_original_unique_regular_n1_component"
                and node.get("chart") == chart
                and node.get("quantity") == field
                and node.get("ordinary_slow_Z_derivative_order") == order
                and node.get("coordinate_domain_upper") == upper
                and node.get("all_words_retained_in_analytic_limit") is True
                and node.get("unique_regular_analytic_solution") is True
                and node.get("raw_units_not_normalized_potential_kernel") is True
                and node.get("source_parameter_bounds_not_selected_as_values") is True
                and node.get("zero_axis_regular_power") == expected_power
            )
            passed = passed and ok
            entries.append({
                "node": index,
                "ordinary_Z_order": order,
                "quantity": field,
                "passed": ok,
            })
        rows[field] = entries
    return {"passed": passed, "chart": chart, "coordinate_domain_upper": upper,
            "fields": rows, "ordinary_Z_orders": list(ORDINARY_Z_ORDERS)}


def _graph_checks(archive: Mapping[str, Any], cutoff_check: Mapping[str, Any],
                  radial_archive: Mapping[str, Any], volterra_archive: Mapping[str, Any],
                  radial_check: Mapping[str, Any], volterra_check: Mapping[str, Any],
                  m15_report: Mapping[str, Any], family: Mapping[str, Any]) -> dict[str, Any]:
    nodes = archive.get("exact_graph_nodes", [])
    _require(isinstance(nodes, list) and len(nodes) == 53165,
             "Pinned cutoff graph node count changed")
    core = _raw_rows(nodes, 43300, "core", 13)
    first = _raw_rows(nodes, 43324, "first", 35824)

    eta_nodes = {
        "fraction": 52487,
        "four_q_minus_two": 52489,
        "sigma": 52490,
        "eta": 52494,
    }
    eta_checks = {
        "q_is_s_over_sc_quarter": nodes[52487].get("operation") == "positive_quotient"
        and nodes[52487].get("numerator") == 35848
        and nodes[52487].get("denominator") == 35824,
        "argument_is_four_q_minus_two": nodes[52488].get("arguments") == [13, 52487]
        and nodes[52489].get("arguments") == [52488, 25304],
        "sigma_is_exact_source": nodes[52490].get("operation") == "exact_original_source_flat_sigma"
        and nodes[52490].get("exact_function_not_its_interval_cap") is True,
        "eta_is_one_minus_sigma": nodes[52494].get("operation") == "sum"
        and nodes[52494].get("arguments") == [1, 52493],
    }
    cutoff_rows = archive.get("inner_cutoff_bindings", {}).get("cutoff_rows", {})
    cutoff_rows_ok = cutoff_rows == {
        "core": [[1, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
        "first": [[52494, 0, 0, 0], [52498, 0, 0, 0], [52506, 0, 0, 0]],
    }

    bindings = archive.get("inner_cutoff_bindings", {})
    fields = bindings.get("fields", {})
    first_fields = fields.get("first", {})
    recovery_handles = {
        "first_Q1": first_fields.get("Q1"),
        "first_V1": first_fields.get("V1"),
        "first_pressure_P1": first_fields.get("pressure_P1"),
        "first_mean1": first_fields.get("rows", {}).get("mean1"),
        "first_K1": first_fields.get("rows", {}).get("K1"),
    }
    recovery_ok = all(isinstance(value, list) and len(value) >= 3
                      for value in recovery_handles.values())

    m15_integrands = bindings.get("integrands", {})
    m15_core = m15_integrands.get("core", {}).get("m15")
    m15_first = m15_integrands.get("first", {}).get("m15")
    cumulative = bindings.get("cumulative", {})
    history_nodes = {
        "core_m15": [52601, 52602, 52603, 52604],
        "rho4_port": [53125, 53126],
        "first_increment": [53137, 53138],
        "combined": [53141, 53142],
        "endpoint_C0_C3": [53161, 53162, 53163, 53164],
    }
    history_checks = {
        "core_m15_integrand_handles": m15_core == [52576, 52577, 52578, 52579],
        "first_m15_integrand_handles": m15_first == [53076, 53077, 53078, 53079],
        "core_m15_C0_C3": cumulative.get("core", {}).get("m15") == history_nodes["core_m15"],
        "core_m15_integrals_are_exact": all(
            nodes[index].get("operation") == "definite_integral"
            and nodes[index].get("paper_moment") == "m15"
            and nodes[index].get("exact_function_integral") is True
            and nodes[index].get("numerical_value_not_installed") is True
            for index in history_nodes["core_m15"]
        ),
        "rho4_port_C0_C1": nodes[53125].get("expression") == 52601
        and nodes[53125].get("value") == 13
        and nodes[53126].get("expression") == 52602
        and nodes[53126].get("value") == 13,
        "first_increment_C0_C1": all(
            nodes[index].get("operation") == "definite_integral"
            and nodes[index].get("paper_moment") == "m15"
            and nodes[index].get("full_core_to_first_history_retained") is True
            and nodes[index].get("leading_N1p_not_cut") is True
            and nodes[index].get("numerical_value_not_installed") is True
            for index in history_nodes["first_increment"]
        ),
        "combined_C0_C1": nodes[53141].get("arguments") == [53125, 53137]
        and nodes[53142].get("arguments") == [53126, 53138],
        "endpoint_C0_C3": all(
            nodes[index].get("operation") == "function_substitution"
            and nodes[index].get("Z_independent_substitution") is True
            for index in history_nodes["endpoint_C0_C3"]
        ),
    }
    lower_nodes = volterra_archive.get("exact_graph_nodes", [])
    lower_checks = {
        "leading_lower_C0_C1": len(lower_nodes) > 36593
        and all(
            lower_nodes[index].get("operation") == "definite_integral"
            and lower_nodes[index].get("original_direct_N1p") is True
            and lower_nodes[index].get("original_m15_minus_applied_once") is True
            and lower_nodes[index].get("signed_integrand_not_a_cap") is True
            for index in (36592, 36593)
        ),
    }
    report_checks = {
        "cutoff_report": cutoff_check.get("all_passed") is True
        and cutoff_check.get(
            "current_original_n1_true_inner_cutoff_and_five_endpoint_functions_installed"
        ) is True,
        "cutoff_own_history": cutoff_check.get("checks", {}).get("own_history", {}).get("passed") is True
        and cutoff_check.get("checks", {}).get("own_history", {}).get(
            "surviving_leading_pressure_source_after_velocity_zero") is True,
        "cutoff_source": cutoff_check.get("checks", {}).get("source", {}).get("passed") is True
        and cutoff_check.get("checks", {}).get("source", {}).get(
            "five_signed_inner_endpoint_functions_not_caps") is True,
        "radial_report": radial_check.get("all_passed") is True
        and radial_check.get("current_original_n1_true_inner_even_radial_functions_installed") is True,
        "volterra_report": volterra_check.get("all_passed") is True
        and volterra_check.get("current_original_n1_true_inner_Volterra_function_adapter_installed") is True,
        "m15_report": m15_report.get("checks", {}).get("passed") is True
        and m15_report.get("source_family") == dict(family),
    }
    families = {
        "cutoff": _family(archive) == dict(family),
        "cutoff_check": _family(cutoff_check) == dict(family),
        "radial": _family(radial_archive) == dict(family),
        "radial_check": _family(radial_check) == dict(family),
        "volterra": _family(volterra_archive) == dict(family),
        "volterra_check": _family(volterra_check) == dict(family),
        "m15": _family(m15_report) == dict(family),
    }
    checks = {
        "graph_node_count": len(nodes) == 53165,
        "raw_core_prefix": core["passed"],
        "raw_first_prefix": first["passed"],
        "eta_graph": all(eta_checks.values()),
        "cutoff_rows": cutoff_rows_ok,
        "recovered_Q1_V1_pressure_handles": recovery_ok,
        "m15_history_graph": all(history_checks.values()),
        "leading_lower_graph": all(lower_checks.values()),
        "reports": all(report_checks.values()),
        "families": all(families.values()),
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "families": families,
        "raw_core": core,
        "raw_first": first,
        "eta": {"nodes": eta_nodes, "checks": eta_checks,
                 "formula": "eta(s)=1-sigma(4*(s/(sc/4))-2)"},
        "cutoff_rows": cutoff_rows,
        "recovery_handles": recovery_handles,
        "history_nodes": history_nodes,
        "history_checks": history_checks,
        "leading_lower_checks": lower_checks,
        "report_checks": report_checks,
        "source_graph_node_count": len(nodes),
    }


def _symbolic_history() -> dict[str, Any]:
    s, t, Z, sc, hb, epsilon, delta, Pstar = sp.symbols(
        "s t Z sc hb epsilon delta Pstar", positive=True)
    Sigma = sp.Function("Sigma")
    F1_fn = sp.Function("F1")
    U1_fn = sp.Function("Uz1")
    Flead_fn = sp.Function("F_lead")
    N1p_fn = sp.Function("N1p")
    Ulead_fn = sp.Function("U_lead")
    M0 = sp.Function("M0")(Z)
    P0 = sp.Function("P0")(Z)
    rho_t = 4 * sp.exp(hb * t)
    rho_s = 4 * sp.exp(hb * s)
    eta_t = 1 - Sigma(16 * t / sc - 2)
    eta_s = eta_t.xreplace({t: s})
    U_t = U1_fn(rho_t, Z)
    F_t = F1_fn(rho_t, Z)
    lead_F_t = Flead_fn(rho_t, Z)
    n1p_t = N1p_fn(rho_t, Z)
    mean_integral = sp.Integral(hb * rho_t * eta_t * U_t, (t, 0, s))
    mean = (4 * M0 + mean_integral) / rho_s
    cut_U = eta_s * U1_fn(rho_s, Z)
    cut_F = eta_s * F1_fn(rho_s, Z)
    K = mean - cut_U
    Q = (2 * Z * cut_U - (1 + delta) * Z * mean
         - (1 - Z**2) * sp.diff(mean, Z)) / (1 - delta * Z**2)
    V = epsilon * rho_s * Q
    pressure_source_s = epsilon * hb * rho_s * (
        2 * Flead_fn(rho_s, Z) * cut_F + N1p_fn(rho_s, Z)
    )
    n1p_raw_s = N1p_fn(rho_s, Z)
    n1p_normalized_s = n1p_raw_s / Pstar**2
    pressure_radial_s = epsilon * (
        2 * Flead_fn(rho_s, Z) * cut_F + n1p_raw_s
    )
    pressure_radial_from_normalized_s = epsilon * (
        2 * Flead_fn(rho_s, Z) * cut_F + Pstar**2 * n1p_normalized_s
    )
    pressure_source_from_normalized_s = epsilon * hb * rho_s * (
        2 * Flead_fn(rho_s, Z) * cut_F + Pstar**2 * n1p_normalized_s
    )
    pressure = P0 + sp.Integral(
        epsilon * hb * rho_t * (2 * lead_F_t * eta_t * F_t + n1p_t),
        (t, 0, s),
    )
    m15 = (
        2 * epsilon * Ulead_fn(rho_s, Z) * cut_U
        - epsilon**2 * rho_s * (
            2 * Flead_fn(rho_s, Z) * cut_F + N1p_fn(rho_s, Z)
        )
    )
    m15_from_normalized = (
        2 * epsilon * Ulead_fn(rho_s, Z) * cut_U
        - epsilon**2 * rho_s * (
            2 * Flead_fn(rho_s, Z) * cut_F + Pstar**2 * n1p_normalized_s
        )
    )
    m13 = epsilon * (
        2 * Flead_fn(rho_s, Z) * cut_F + N1p_fn(rho_s, Z)
    )
    m15_increment = sp.Integral(
        m15.xreplace({s: t}) * hb * rho_t, (t, 0, s)
    )
    core_m15 = sp.Function("C_core_m15")(Z)
    inner_m15 = core_m15 + m15_increment
    residuals = {
        "mean_history_ODE": sp.diff(mean, s) - hb * (cut_U - mean),
        "K_definition": K - (mean - cut_U),
        "Q_definition": Q - (
            2 * Z * cut_U - (1 + delta) * Z * mean
            - (1 - Z**2) * sp.diff(mean, Z)
        ) / (1 - delta * Z**2),
        "V_definition": V - epsilon * rho_s * Q,
        "pressure_FTC": sp.diff(pressure, s) - pressure_source_s,
        "pressure_radial_source_scale_restore": (
            pressure_radial_from_normalized_s - pressure_radial_s
        ),
        "pressure_history_source_scale_restore": (
            pressure_source_from_normalized_s - pressure_source_s
        ),
        "raw_N1p_scale_restore": Pstar**2 * n1p_normalized_s - n1p_raw_s,
        "m15_normalized_pressure_scale_restore": m15_from_normalized - m15,
        "m15_history_FTC": sp.diff(m15_increment, s) - m15 * hb * rho_s,
        "inner_m15_join": inner_m15 - core_m15 - m15_increment,
    }
    base_checks = {name: sp.simplify(value) == 0
                   for name, value in residuals.items()}
    ordinary_Z = {
        name: [sp.simplify(sp.diff(value, Z, order)) == 0
               for order in ORDINARY_Z_ORDERS]
        for name, value in residuals.items()
    }
    pressure_after_eta_zero = sp.simplify(pressure_source_s.subs(eta_s, 0))
    m15_after_eta_zero = sp.simplify(m15.subs(eta_s, 0))
    m15_pressure_coefficient = sp.expand(m15).coeff(N1p_fn(rho_s, Z))
    checks = {
        **base_checks,
        "ordinary_Z0_to_Z3_residuals": all(
            all(values) for values in ordinary_Z.values()
        ),
        "leading_N1p_survives_eta_zero": pressure_after_eta_zero.has(
            N1p_fn(rho_s, Z)
        ) and m15_after_eta_zero.has(N1p_fn(rho_s, Z)),
        "m15_negative_pressure_coefficient": sp.simplify(
            m15_pressure_coefficient + epsilon**2 * rho_s
        ) == 0,
        "m15_pressure_weight_is_epsilon2_rho": sp.simplify(
            sp.expand(m15).coeff(N1p_fn(rho_s, Z))
            + epsilon**2 * rho_s
        ) == 0,
        "Pstar_squared_restores_raw_N1p": sp.simplify(
            Pstar**2 * n1p_normalized_s - n1p_raw_s
        ) == 0,
        "m15_Pstar_squared_consumed_once": sp.simplify(
            m15_from_normalized - m15
        ) == 0,
        "m15_distinct_from_m13": sp.simplify(m15 - m13) != 0,
        "epsilon_and_hb_retained": pressure_source_s.has(epsilon)
        and pressure_source_s.has(hb)
        and m15.has(epsilon),
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "ordinary_Z_residuals": ordinary_Z,
        "base_residuals": base_checks,
        "formulas": {
            "q": "q=s/(sc/4)=4*s/sc",
            "eta": "1-sigma(4*q-2)",
            "rho": "4*exp(hb*s)",
            "R": "epsilon*rho",
            "mean": "[4*M_core(4,Z)+Integral_0^s hb*rho(t)*eta(t)*Uz1(rho(t),Z)dt]/rho(s)",
            "K1": "mean-eta*Uz1",
            "Q1": "[2*Z*(eta*Uz1)-(1+delta)*Z*mean-(1-Z^2)*mean_Z]/(1-delta*Z^2)",
            "V1": "epsilon*rho*Q1",
            "pressure": "P_core(4,Z)+Integral_0^s epsilon*hb*rho(t)*(2*F_lead*eta*F1+N1p)dt",
            "N1p_units": "raw_system N1p is physical; numeric normalized source is N1p/Pstar^2",
            "pressure_radial_function": "P1_rho=epsilon*(2*F_lead*cut_F1+N1p)",
            "pressure_radial_function_from_normalized": "P1_rho=epsilon*(2*F_lead*cut_F1+Pstar^2*(N1p/Pstar^2))",
            "m15_density_per_drho": "2*epsilon*U_lead*(eta*Uz1)-epsilon^2*rho*(2*F_lead*eta*F1+N1p)",
            "m15_pressure_from_normalized": "-epsilon^2*rho*(2*F_lead*eta*F1+Pstar^2*(N1p/Pstar^2))",
            "m15_first_increment": "Integral_0^s m15_density(rho(t),Z)*hb*rho(t)dt",
            "m15_inner": "m15_core(Z)+m15_first_increment(s,Z)",
        },
        "leading_pressure_after_eta_zero": str(pressure_after_eta_zero),
        "m15_after_eta_zero": str(m15_after_eta_zero),
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "source_executed": False,
    }


def _parent_pressure_units(saved_parent: Mapping[str, Any]) -> dict[str, Any]:
    """Consume exact parent pressure/source handles without selecting values.

    The raw Volterra ``N1p`` is physical.  The older finite evaluator exposes
    ``N1p/Pstar^2`` only after normalization.  This gate binds both statements
    to the saved source AST and row-4 source handles; it does not infer them
    from matching strings in a numerical packet.
    """

    partition = saved_parent.get("actual_collar_pressure_partition", {})
    core = saved_parent.get("actual_core_pressure_m15_primitive", {})
    partition_checks = partition.get("checks", {})
    core_checks = core.get("checks", {})
    ast_bindings = core.get("source_AST_bindings", {}).get("bindings", {})
    raw_system = ast_bindings.get("raw_system", {})
    raw_matrix = ast_bindings.get("raw_matrix_expressions", {})
    raw_pressure = ast_bindings.get("raw_pressure_source_binding", {})
    row_equation = core.get("actual_raw_pressure_rho_equation", {})
    scales = core.get("normalization_and_scale_binding", {})
    checks = {
        "parent_partition_passed": partition.get("passed") is True,
        "same_leading_fields": (
            partition.get("full_pressure_partition_for_same_leading_fields")
            is True
        ),
        "partition_exact_pressure": (
            partition_checks.get("exact_pressure_partition") is True
        ),
        "physical_Pstar_squared_canceled_once": (
            partition.get("physical_Pstar_squared_canceled_once") is True
        ),
        "core_pressure_packet_passed": core.get("passed") is True,
        "core_pressure_checks_passed": core_checks.get("pressure_equation") is True
        and core_checks.get("exact_scales") is True,
        "raw_system_ast_bound": raw_system.get("passed") is True,
        "raw_matrix_ast_bound": raw_matrix.get("passed") is True,
        "raw_pressure_source_ast_bound": raw_pressure.get("passed") is True,
        "raw_pressure_source_is_N1p": (
            raw_pressure.get("raw_system_field")
            == "fields['F_rho_rows'][0].value"
            and raw_pressure.get("field_meaning") == "f=F_lead=F0*Phi"
        ),
        "row4_raw_source": (
            raw_pressure.get("row4_literal")
            == "B0['4,1']=epsilon*4*x*F_lead"
            and row_equation.get("g4") == "2*epsilon*x*N1p"
        ),
        "physical_pressure_rho_equation": (
            row_equation.get("P1_rho")
            == "epsilon*(2*F_lead*F1+N1p)"
        ),
        "normalization_applied_once": scales.get("powers_applied_once") is True,
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "raw_system_N1p": "physical N1p; no Pstar^2 division",
        "physical_pressure_radial_function": (
            "P1_rho=epsilon*(2*F_lead*cut_F1+N1p)"
        ),
        "numeric_normalized_pressure_source": "N1p/Pstar^2",
        "physical_m15_pressure_weight": "-epsilon^2*rho*(2*F_lead*cut_F1+N1p)",
        "Pstar_squared_consumption": (
            "restore N1p once in normalized source, then apply physical "
            "epsilon^2*rho m15 weight once"
        ),
        "source_executed": False,
    }


def _raw_gate(raw_binding: Mapping[str, Any] | None,
              family: Mapping[str, Any]) -> dict[str, Any]:
    if raw_binding is None:
        return {
            "provided": False,
            "passed": False,
            "external_root_gate": True,
            "missing": "raw_binding with exact source coefficient enclosures",
            "function_source_gate": {
                "passed": False,
                "missing": "actual first-raw function source binding",
            },
            "numerical_prefix_gate": {
                "passed": False,
                "missing": "explicit actual_first_raw_numerical_prefix_bound",
            },
            "raw_prefix_numerical_gate": False,
            "source_executed": False,
        }
    nested = raw_binding.get("actual_forward_prefix_bound")
    candidate = nested if isinstance(nested, Mapping) else raw_binding
    archive_schema = candidate.get("archive_schema", {})
    candidate_family = _family(candidate)
    source_family = (
        candidate_family
        if any(value is not None for value in candidate_family.values())
        else _family(raw_binding)
    )
    family_or_parent_pin = (
        source_family == dict(family)
        or raw_binding.get("parent_mapping_matches_file") is True
    )
    first_raw = raw_binding
    if not any(key in first_raw for key in (
        "actual_first_raw_function_bound",
        "actual_first_raw_all_word_convergence_bound",
        "actual_first_raw_core_entry_bound",
    )) and isinstance(raw_binding.get("actual_first_raw_function"), Mapping):
        first_raw = raw_binding["actual_first_raw_function"]
    function_checks = {
        "actual_first_raw_function_bound": (
            first_raw.get("actual_first_raw_function_bound") is True
        ),
        "actual_first_raw_all_word_convergence_bound": (
            first_raw.get("actual_first_raw_all_word_convergence_bound") is True
        ),
        "actual_first_raw_core_entry_bound": (
            first_raw.get("actual_first_raw_core_entry_bound") is True
        ),
        "source_executed": (
            first_raw.get("source_executed") is False
            and first_raw.get("archive_source_executed") is not True
        ),
        "source_family": (
            _family(first_raw) == dict(family)
            or source_family == dict(family)
        ),
        "ordinary_Z_orders": (
            first_raw.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS
            or candidate.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS
            or archive_schema.get("ordinary_Z_orders") == ORDINARY_Z_ORDERS
        ),
        "source_only": (
            first_raw.get("source_only") is True
            or (
                first_raw.get("accepted_owner_builder_fixture_recurrence_quadrature_replayed")
                is False
                and first_raw.get("point_values_selected") is False
            )
        ),
    }
    # A finite canonical prefix packet is not the first-raw function gate.  It
    # only becomes a numerical gate when the producer explicitly supplies this
    # separate admission flag; shape, ``passed``, and archive schema alone are
    # insufficient.
    numerical_checks = {
        "explicit_actual_first_raw_numerical_prefix_bound": (
            first_raw.get("actual_first_raw_numerical_prefix_bound") is True
            or candidate.get("actual_first_raw_numerical_prefix_bound") is True
        ),
        "prefix_passed": (
            first_raw.get("passed") is True or candidate.get("passed") is True
        ),
        "source_only": raw_binding.get("source_only") is True
        or (
            candidate.get("accepted_producers_replayed") is False
            and (
                candidate.get("point_values_selected") is False
                or candidate.get("coefficient_intervals_are_not_selected_points") is True
                or archive_schema.get("approximation_decimal_used") is False
            )
        ),
        "source_not_executed": (
            raw_binding.get("source_executed") is not True
            and candidate.get("source_executed") is False
        ),
        "family_or_parent_pin": family_or_parent_pin,
        "ordinary_Z_orders": function_checks["ordinary_Z_orders"],
    }
    function_passed = all(function_checks.values())
    numerical_passed = all(numerical_checks.values())
    return {
        "provided": True,
        "passed": function_passed,
        "external_root_gate": True,
        "checks": function_checks,
        "function_source_gate": {
            "passed": function_passed,
            "checks": function_checks,
            "required_contract": [
                "actual_first_raw_function_bound=True",
                "actual_first_raw_all_word_convergence_bound=True",
                "actual_first_raw_core_entry_bound=True",
                "source_executed=False",
                "source_family matches parent",
                "ordinary_Z_orders=[0,1,2,3]",
            ],
        },
        "numerical_prefix_gate": {
            "passed": numerical_passed,
            "checks": numerical_checks,
            "required_contract": [
                "actual_first_raw_numerical_prefix_bound=True",
                "source-only exact coefficient enclosures",
                "ordinary_Z_orders=[0,1,2,3]",
            ],
        },
        "raw_prefix_numerical_gate": numerical_passed,
        "source_executed": False,
        "source_family": source_family,
        "binding_shape": "actual_forward_prefix_bound" if isinstance(nested, Mapping)
        else "direct",
    }


def build(root: str | Path, parent: Mapping[str, Any],
          raw_binding: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Build the source-only cutoff/history identity packet."""

    repo = Path(root)
    parent_path = repo / PARENT
    _require(_sha(parent_path) == PARENT_SHA256, "Pinned native-collar parent changed")
    saved_parent = _load(parent_path)
    _require(dict(parent) == dict(saved_parent),
             "Supplied parent mapping differs from pinned native-collar report")
    family = _family(saved_parent)
    _require(all(value is not None for value in family.values()),
             "Parent source family is incomplete")

    specs = {
        "parent": (PARENT, PARENT_SHA256),
        "cutoff_source": (CUTOFF_SOURCE, CUTOFF_SOURCE_SHA256),
        "cutoff_archive": (CUTOFF_ARCHIVE, CUTOFF_ARCHIVE_SHA256),
        "cutoff_check": (CUTOFF_CHECK, CUTOFF_CHECK_SHA256),
        "radial_source": (RADIAL_SOURCE, RADIAL_SOURCE_SHA256),
        "radial_archive": (RADIAL_ARCHIVE, RADIAL_ARCHIVE_SHA256),
        "radial_check": (RADIAL_CHECK, RADIAL_CHECK_SHA256),
        "volterra_source": (VOLTERRA_SOURCE, VOLTERRA_SOURCE_SHA256),
        "volterra_archive": (VOLTERRA_ARCHIVE, VOLTERRA_ARCHIVE_SHA256),
        "volterra_check": (VOLTERRA_CHECK, VOLTERRA_CHECK_SHA256),
        "cutoff_numeric_source": (CUTOFF_NUMERIC_SOURCE, CUTOFF_NUMERIC_SOURCE_SHA256),
        "m15_numeric_source": (M15_NUMERIC_SOURCE, M15_NUMERIC_SOURCE_SHA256),
        "m15_report": (M15_REPORT, M15_REPORT_SHA256),
    }
    source_pins = {name: _pin(repo, relative, expected)
                   for name, (relative, expected) in specs.items()}
    source_pins_ok = all(pin["passed"] for pin in source_pins.values())

    cutoff_archive = _load_gzip(repo / CUTOFF_ARCHIVE)
    radial_archive = _load_gzip(repo / RADIAL_ARCHIVE)
    volterra_archive = _load_gzip(repo / VOLTERRA_ARCHIVE)
    cutoff_check = _load(repo / CUTOFF_CHECK)
    radial_check = _load(repo / RADIAL_CHECK)
    volterra_check = _load(repo / VOLTERRA_CHECK)
    m15_report = _load(repo / M15_REPORT)
    graph = _graph_checks(
        cutoff_archive, cutoff_check, radial_archive, volterra_archive,
        radial_check, volterra_check, m15_report, family,
    )
    source_ast, source_ast_ok = _source_ast(repo)
    symbolic = _symbolic_history()
    pressure_units = _parent_pressure_units(saved_parent)
    raw_gate = _raw_gate(raw_binding, family)

    parent_checks = saved_parent.get("actual_native_collar_transport", {})
    parent_scope = {
        "passed": saved_parent.get("checks", {}).get("passed") is True,
        "native_transport_passed": parent_checks.get("passed") is True,
        "raw_n1_identity_still_open": parent_checks.get("raw_n1_cut_response_identity") is False,
        "six_terms_not_removed": parent_checks.get("six_signed_source_equality") is False,
        "source_not_executed": parent_checks.get(
            "accepted_ancestor_or_quadrature_replayed") is False,
    }
    parent_scope_ok = all(parent_scope.values())
    checks = {
        "parent_hash_and_mapping": True,
        "parent_scope": parent_scope_ok,
        "source_pins": source_pins_ok,
        "saved_graphs": graph["passed"],
        "source_AST": source_ast_ok,
        "fresh_symbolic_history": symbolic["passed"],
        "N1p_units_and_Pstar_once": pressure_units["passed"],
        "raw_input_gate_external": True,
        "no_numeric_values_promoted": True,
    }
    passed = all(checks.values())
    return {
        "passed": passed,
        "source_only": True,
        "source_executed": False,
        "parent_report": PARENT,
        "parent_report_sha256": PARENT_SHA256,
        "source_family": family,
        "checks": checks,
        "source_pins": source_pins,
        "parent_scope": parent_scope,
        "source_AST": source_ast,
        "graph_binding": graph,
        "fresh_symbolic_history": symbolic,
        "N1p_units": pressure_units,
        "raw_input_gate": raw_gate,
        "raw_input_gate_is_external": True,
        "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        "raw_state": {
            "core_nodes": list(range(43300, 43324)),
            "first_nodes": list(range(43324, 43348)),
            "fields": list(RAW_FIELDS),
            "raw_units_not_normalized_potential_kernel": True,
            "point_values_selected": False,
        },
        "cutoff_identity": {
            "eta": "1-sigma(4*(s/(sc/4))-2)",
            "only_F1_and_Uz1_are_cut": True,
            "K1_Q1_V1_pressure_are_recovered_from_own_history": True,
            "leading_N1p_remains_after_eta_zero": True,
            "core_inlet_not_reset": True,
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
        },
        "history_identity": {
            "core_m15_C0_C3_nodes": [52601, 52602, 52603, 52604],
            "rho4_port_nodes": [53125, 53126],
            "first_increment_C0_C1_nodes": [53137, 53138],
            "combined_C0_C1_nodes": [53141, 53142],
            "C_inner_endpoint_C0_C3_nodes": [53161, 53162, 53163, 53164],
            "leading_lower_m15_C0_C1_nodes": [36592, 36593],
            "m15_distinct_from_m13": True,
            "signed_pressure_term_negative_once": True,
            "six_signed_terms_retained": ["m11", "m12", "m13", "m14", "m15"],
        },
        "limits": {
            "raw_prefix_numerical_gate": raw_gate.get(
                "numerical_prefix_gate", {}
            ).get("passed") is True,
            "full_raw_n1_function_identity": False,
            "full_global_five_moment_closure": False,
            "finite_numeric_values_as_functions": False,
            "m13_substituted_for_m15": False,
            "accepted_owner_or_quadrature_replayed": False,
        },
    }
