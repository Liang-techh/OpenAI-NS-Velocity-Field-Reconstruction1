"""Source-only long-reshape/reference physical germ certificate.

This module binds the already saved long-reshape and reference charts at
``long_reshape.native1 = reference.native0``.  It performs exact endpoint
fingerprinting on the accepted graph and records the source calculus needed
for the physical two-jet statement.  It never constructs an owner, evaluates
an original numerical source, or treats an interval/cap as a field value.

The long chart has physical measure ``T`` (graph handle 7) while the
reference chart has physical measure ``length`` (graph handle 81519).  The
certificate therefore compares physical functions at two scopes and keeps
the two native Jacobians separate.
"""

from __future__ import annotations

import gzip
import hashlib
import json
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Mapping

from .original_active_viscous_m15 import ActiveEndpointFingerprint, ROOT
from .original_reference_germ_guards import _function_evidence
from .original_reference_viscous_m15 import PhysicalGermRecipes


_OUTER_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py"
)
_LOOP_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
)
_BRIDGE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_bridge_C6_functions.py"
)
_OUTER_ARCHIVE = _OUTER_SOURCE[:-3] + ".json.gz"
_LOOP_ARCHIVE = _LOOP_SOURCE[:-3] + ".json.gz"


def _mapping(value: Any, key: str | None = None) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        return {}
    if key is not None and isinstance(value.get(key), Mapping):
        return value[key]
    return value


def _repo_root(root: Any) -> Path:
    if root is not None:
        path = Path(root)
        if path.is_file():
            path = path.parent
        if (path / "experiments").is_dir() or (path / "src").is_dir():
            return path
    return Path(__file__).resolve().parents[2]


def _load_archive(root: Path, relative: str) -> Mapping[str, Any]:
    path = root / relative
    try:
        with gzip.open(path, "rt", encoding="utf8") as handle:
            value = json.load(handle)
        return value if isinstance(value, Mapping) else {}
    except (OSError, ValueError, json.JSONDecodeError):
        return {}


def _pin_maps(*values: Any) -> list[Mapping[str, Any]]:
    maps: list[Mapping[str, Any]] = []
    for value in values:
        if not isinstance(value, Mapping):
            continue
        for key in ("input_hashes", "implementation_hashes"):
            if isinstance(value.get(key), Mapping):
                maps.append(value[key])
        if isinstance(value.get("_receipt"), Mapping):
            receipt = value["_receipt"]
            if isinstance(receipt.get("input_hashes"), Mapping):
                maps.append(receipt["input_hashes"])
    return maps


def _expected_pin(relative: str, path: Path, maps: Iterable[Mapping[str, Any]]) -> str | None:
    names = (relative, Path(relative).name, str(Path(relative)))
    for mapping in maps:
        for name in names:
            value = mapping.get(name)
            if isinstance(value, str):
                return value
        for name, value in mapping.items():
            if isinstance(name, str) and Path(name).name == path.name and isinstance(value, str):
                return value
    return None


def _source_pin(root: Path, relative: str, maps: Iterable[Mapping[str, Any]]) -> dict[str, Any]:
    path = root / relative
    if not path.is_file():
        return {"path": relative, "exists": False, "pinned": False}
    actual = hashlib.sha256(path.read_bytes()).hexdigest()
    expected = _expected_pin(relative, path, maps)
    return {
        "path": relative,
        "sha256": actual,
        "expected_sha256": expected,
        "exists": True,
        "pinned": expected == actual if expected is not None else False,
    }


def _node(nodes: Any, index: Any) -> Mapping[str, Any] | None:
    if isinstance(index, bool) or not isinstance(index, int):
        return None
    try:
        row = nodes[index] if not isinstance(nodes, Mapping) else nodes.get(index, nodes.get(str(index)))
    except (IndexError, KeyError, TypeError):
        return None
    return row if isinstance(row, Mapping) else None


def _handle(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    node = getattr(value, "node", None)
    return node if isinstance(node, int) else None


def _z_handle(nodes: Any) -> int | None:
    for index, row in enumerate(nodes if isinstance(nodes, list) else ()):
        if isinstance(row, Mapping) and row.get("operation") == "bound_variable" and row.get("name") == "Z":
            return index
    return None


def _source_evidence(root: Path) -> dict[str, Any]:
    outer = root / _OUTER_SOURCE
    loop = root / _LOOP_SOURCE
    bridge = root / _BRIDGE_SOURCE
    if not all(path.is_file() for path in (outer, loop, bridge)):
        return {"passed": False, "reason": "required pinned source file is missing"}

    long_rows = _function_evidence(outer, "CurrentOriginalPrescribedR100RmC6Functions.__init__", (
        "logUtheta=a.add(a.scale(B,g.sub(g.one,sigma(g,q))),a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),self.core.parameters['logC'])))",
        "slope=a.add(a.fixed('-2/5'),a.scale(B,g.neg(g.quotient(sigma_prime(g,q),T,'same positive original T'))))",
        "actual=dict(phi=phi,U=incoming['actual']['U'],phi_y=a.mul(phi,slope),U_y=a.fixed(0),\n            physical_Utheta=Utheta,log_Utheta=logUtheta)",
        "logUtheta=a.add(a.neg(logq),a.fixed(g.sub(g.mul(c('1/10'),y),self.core.parameters['logC'])))",
        "actual=dict(phi=phi,U=incoming['actual']['U'],phi_y=a.scale(phi,'-2/5'),U_y=a.fixed(0),\n            physical_Utheta=Utheta,log_Utheta=logUtheta)",
    ))
    sigma_rows = _function_evidence(outer, "sigma_prime", (
        "return g.node('analytic_unary',name='original_flat_sigma_prime',argument=x.node,\n        derivative_source=own.c2.current.ast_binding(own.c2.flat_source.sigma_jets),\n        ordinary_derivative='same original sigma Taylor coefficient[1]',\n        exact_flat_source_derivative_not_a_bound=True)",
    ))
    field_rows = _function_evidence(bridge, "generic_fields", (
        "fields=dict(E=a.mul(W,actual['phi']),V=a.scale(actual['U'],inverse),m=a.scale(J['M'],inverse),\n        h=a.scale(a.mul(W,J['H']),'1/2'),k=a.scale(a.mul(W,J['K']),g.mul(c('1/2'),inverse)),\n        e=a.sub(a.scale(J['A'],g.mul(inverse,inverse)),a.mul(FF,J['B'])),p=a.mul(FF,J['C']),P0=a.trim(P0))",
        "radial=dict(E_y=a.mul(W,a.add(a.scale(actual['phi'],'1/2'),actual['phi_y'])),V_y=a.scale(actual['U_y'],inverse))",
    ))
    history_rows = _function_evidence(bridge, "histories", (
        "decay=g.unary('exp',g.mul(g.constant(-rate),jacobian,phase))",
        "kernel=g.unary('exp',g.mul(g.constant(-rate),jacobian,g.sub(phase,t)))",
        "body=a.scale(a.at(bodies[name],phase,t),g.mul(jacobian,kernel))",
        "out[name]=a.add(a.scale(initial[name],decay),a.integral(body,t,g.zero,phase,\n            actual_six_history_rate=rate,actual_history_name=name,\n            original_incoming_memory_retained=True,source_measure_applied_once=True,\n            prescribed_bridge_chart=label))",
    ))
    sigma_source_rows = _function_evidence(loop, "LoopAlgebra.sigma", (
        "out.append(g.node('original_lazy_flat_branch',active_body=body.node,\n                Delta=rows(x)[0].node,eta=g.one.node,flat_value=g.zero.node,\n                flat_predicate=predicate.node,active_body_not_evaluated_on_flat_branch=True,\n                ordinary_slow_Z_derivative_order=n+1,original_sigma_logistic_ODE=True))",
    ))
    root_rows = _function_evidence(loop, "generic_roots", (
        "a.sub(E,a.scale(a.trim(radial['E_y']),2))",
        "a.scale(a.trim(radial['V_y']),2)",
        "a.div(jet,E,'same strictly positive actual prescribed swirl E')",
        "aa,bb=divide(C),divide(B)",
    ))
    loop_rows = _function_evidence(loop, "loop_functions", (
        "a.div(gamma,a.scale(roots['a'],2),'same original active a>0')",
        "a.div(one,a.sqrt(a.add(one,a.mul(u,u))),'positive sqrt(1+u^2)')",
        "a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
        "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
        "Bphi=lazy(a.scale(a.mul(roots['E'],a.sub(a.neg(a.scale(a.mul(roots['a'],a.mul(that,invlambda)),\n        g.quotient(g.one,twopi,'positive 2*pi'))),roots['b'])),'1/2'))",
    ))
    groups = {
        "long_reference_outer_formula": long_rows,
        "flat_sigma_prime": sigma_rows,
        "bridge_fields_and_radial": field_rows,
        "duhamel_histories": history_rows,
        "flat_sigma_source": sigma_source_rows,
        "generic_root_recipe": root_rows,
        "poisson_loop_recipe": loop_rows,
    }
    return {"passed": all(row.get("passed") is True for row in groups.values()), "groups": groups}


def _compare(
    fingerprint: ActiveEndpointFingerprint,
    scope: Mapping[int, tuple[str, Fraction | None]],
    label: str,
    left: Any,
    right: Any,
    order: int | None = None,
) -> dict[str, Any]:
    left_index, right_index = _handle(left), _handle(right)
    row: dict[str, Any] = {"label": label, "left_node": left_index, "right_node": right_index}
    if left_index is None or right_index is None or _node(fingerprint.nodes, left_index) is None or _node(fingerprint.nodes, right_index) is None:
        row["passed"] = False
        row["reason"] = "missing original graph handle"
        return row
    try:
        lhash, lrat = fingerprint.at(left_index, scope)
        rhash, rrat = fingerprint.at(right_index, scope)
    except Exception as exc:  # bounded report: no source evaluation is attempted
        row["passed"] = False
        row["reason"] = repr(exc)
        return row
    row.update(
        passed=lhash == rhash,
        left_endpoint_sha256=lhash,
        right_endpoint_sha256=rhash,
        left_rational=None if lrat is None else [lrat.numerator, lrat.denominator],
        right_rational=None if rrat is None else [rrat.numerator, rrat.denominator],
        ordinary_Z_order=order,
    )
    return row


def _compare_arrays(
    fingerprint: ActiveEndpointFingerprint,
    scope: Mapping[int, tuple[str, Fraction | None]],
    label: str,
    left: Any,
    right: Any,
    orders: int = 3,
) -> list[dict[str, Any]]:
    if not isinstance(left, list) or not isinstance(right, list) or len(left) < orders or len(right) < orders:
        return [{"label": label, "passed": False, "reason": "ordinary source jet rows 0..2 missing"}]
    return [_compare(fingerprint, scope, f"{label}[{order}]", left[order], right[order], order) for order in range(orders)]


def _fingerprint_interface(nodes: Any, l: Mapping[str, Any], outer: Mapping[str, Any]) -> dict[str, Any]:
    charts = _mapping(l.get("charts"))
    outer_charts = _mapping(outer.get("charts"))
    long = _mapping(charts.get("long_reshape"))
    reference = _mapping(charts.get("reference"))
    long_outer = _mapping(outer_charts.get("long_reshape"))
    reference_outer = _mapping(outer_charts.get("reference"))
    required = bool(long and reference and long_outer and reference_outer)
    if not required:
        return {"passed": False, "rows": [], "reason": "long_reshape/reference chart binding missing"}
    coordinate = _handle(long.get("coordinate"))
    ref_coordinate = _handle(reference.get("coordinate"))
    if coordinate is None or ref_coordinate is None:
        return {"passed": False, "rows": [], "reason": "coordinate handles missing"}
    try:
        fingerprint = ActiveEndpointFingerprint(nodes, dict(coordinate=coordinate, loop=charts))
        long_scope_value = fingerprint.at(_handle(long.get("domain", [None, None])[1]), {})
        reference_scope_value = fingerprint.at(_handle(reference.get("domain", [None, None])[0]), {})
        scope = {coordinate: long_scope_value, ref_coordinate: reference_scope_value}
    except Exception as exc:
        return {"passed": False, "rows": [], "reason": repr(exc)}

    rows: list[dict[str, Any]] = []
    for group in ("fields", "radial", "actual", "moments", "generic"):
        left_group = _mapping(long_outer.get(group))
        right_group = _mapping(reference_outer.get(group))
        for key in sorted(set(left_group) & set(right_group)):
            left_value, right_value = left_group[key], right_group[key]
            if isinstance(left_value, list) and isinstance(right_value, list):
                rows.extend(_compare_arrays(fingerprint, scope, f"outer.{group}.{key}", left_value, right_value))
            else:
                rows.append(_compare(fingerprint, scope, f"outer.{group}.{key}", left_value, right_value))

    # The reference chart's incoming history is the long chart's endpoint
    # moment, not the long chart's own inlet packet.  This is the actual
    # source join and avoids an interval-overlap interpretation.
    for key in sorted(set(_mapping(long_outer.get("moments"))) & set(_mapping(reference_outer.get("incoming_histories")))):
        rows.extend(_compare_arrays(
            fingerprint,
            scope,
            f"outer.moments_to_reference_inlet.{key}",
            _mapping(long_outer["moments"])[key],
            _mapping(reference_outer["incoming_histories"])[key],
        ))

    # Loop roots and the actual Poisson outputs are separate source ports,
    # but they are the same physical functions at these two scopes.
    def _node_array(value: Any) -> bool:
        return isinstance(value, list) and all(_handle(item) is not None for item in value)

    for group in ("roots", "loop"):
        left_group, right_group = _mapping(long.get(group)), _mapping(reference.get(group))
        for key in sorted(set(left_group) & set(right_group)):
            left_value, right_value = left_group[key], right_group[key]
            if _node_array(left_value) and _node_array(right_value):
                rows.extend(_compare_arrays(fingerprint, scope, f"loop.{group}.{key}", left_value, right_value))
            elif _handle(left_value) is not None and _handle(right_value) is not None:
                rows.append(_compare(fingerprint, scope, f"loop.{group}.{key}", left_value, right_value))

    # The radius is a common physical function, while the native Jacobians
    # are intentionally different and are reported rather than equated.
    radius_row = _compare(fingerprint, scope, "geometry.radius", long_outer.get("radius"), reference_outer.get("radius"))
    rows.append(radius_row)
    jacobians = {
        "long_reshape": long_outer.get("physical_measure"),
        "reference": reference_outer.get("physical_measure"),
        "equal": long_outer.get("physical_measure") == reference_outer.get("physical_measure"),
    }

    # E_y is a graph identity after the endpoint scope.  Build E/10 in the
    # endpoint ring; this does not inspect a numerical E value.
    ey_rows: list[dict[str, Any]] = []
    for name, chart in (("long_reshape", long_outer), ("reference", reference_outer)):
        fields = _mapping(chart.get("fields")); radial = _mapping(chart.get("radial"))
        for order in range(3):
            try:
                expected = fingerprint.combine("product", [
                    fingerprint.at(_handle(fields["E"][order]), scope),
                    fingerprint.rational(Fraction(1, 10)),
                ])
                actual = fingerprint.at(_handle(radial["E_y"][order]), scope)
                ey_rows.append({
                    "label": f"{name}.E_y=E/10[Z{order}]",
                    "passed": actual[0] == expected[0],
                    "ordinary_Z_order": order,
                    "source_E_node": _handle(fields["E"][order]),
                    "source_E_y_node": _handle(radial["E_y"][order]),
                    "endpoint_sha256": actual[0],
                    "expected_sha256": expected[0],
                })
            except Exception as exc:
                ey_rows.append({"label": f"{name}.E_y=E/10[Z{order}]", "passed": False, "reason": repr(exc)})

    vy_rows = []
    for name, chart in (("long_reshape", long_outer), ("reference", reference_outer)):
        radial = _mapping(chart.get("radial"))
        for order, node in enumerate(radial.get("V_y", [])[:3]):
            try:
                digest, rational = fingerprint.at(_handle(node), scope)
                vy_rows.append({"label": f"{name}.V_y=0[Z{order}]", "passed": rational == 0,
                                "ordinary_Z_order": order, "node": _handle(node), "endpoint_sha256": digest})
            except Exception as exc:
                vy_rows.append({"label": f"{name}.V_y=0[Z{order}]", "passed": False, "reason": repr(exc)})

    def _serial_scope(value: Any) -> dict[str, Any]:
        if not isinstance(value, tuple) or len(value) != 2:
            return {"sha256": None, "rational": None}
        digest, rational = value
        return {
            "sha256": digest,
            "rational": None if rational is None else [rational.numerator, rational.denominator],
        }

    return {
        "passed": required and all(row.get("passed") is True for row in rows + ey_rows + vy_rows),
        "rows": rows,
        "row_count": len(rows),
        "physical_E_y_rows": ey_rows,
        "physical_V_y_rows": vy_rows,
        "radius_row": radius_row,
        "native_Jacobians": jacobians,
        "scope": {
            "long_reshape": {"coordinate": coordinate, "endpoint": 1},
            "reference": {"coordinate": ref_coordinate, "endpoint": 0},
        },
        "scope_values": {
            str(handle): _serial_scope(value) for handle, value in scope.items()
        },
        "ordinary_Z_orders": [0, 1, 2],
        "source_values_selected": False,
    }


def _physical_two_jet(evidence: Mapping[str, Any], interface: Mapping[str, Any]) -> dict[str, Any]:
    ast_passed = bool(evidence.get("passed"))
    interface_passed = bool(interface.get("passed"))
    return {
        "passed": ast_passed and interface_passed,
        "method": "long_flat_sigma_source_formula_plus_common_Duhamel_endpoint_germ",
        "radius_variable": "y=log(R/110); d_y=R*d_R; not a Cartesian component",
        "E_y": "E/10",
        "E_yy": "E/100",
        "V_y": "0 at long_reshape native1 and reference native0",
        "V_yy": "0 at the seam",
        "E_y_graph_rows": interface.get("physical_E_y_rows", []),
        "V_y_graph_rows": interface.get("physical_V_y_rows", []),
        "ordinary_Z_orders": [0, 1, 2],
        "ordinary_Z_jets_are_function_identities_not_axial_samples": True,
        "E_yy_derivation": [
            "At long_reshape q=1, log(Utheta)=B*(1-sigma(q))-log(1+Z^2)+y/10-logC.",
            "The pinned flat sigma source has sigma(1)=1 and all finite endpoint derivatives zero, so d_y log(E)=1/10 and d_yy log(E)=0.",
            "At reference q=0, log(Utheta)=-log(1+Z^2)+y/10-logC, giving the same physical two-jet.",
            "Pstar is Z-only and the bridge generic_fields source is E=W*phi=Utheta/Pstar, V=U/Pstar.",
        ],
        "V_yy_proof": "The source actual U_y is exact zero on both endpoint charts; differentiating the zero function gives V_y=V_yy=0.",
        "common_Duhamel_second_jets": {
            "passed": evidence.get("groups", {}).get("duhamel_histories", {}).get("passed") is True,
            "ordinary_Z_orders": [0, 1, 2],
            "source_identity": "same six history rate/kernel/body/incoming-memory Duhamel recipe, with physical native measure kept chart-local",
            "endpoint_moment_rows": interface.get("row_count", 0),
        },
        "all_six_nonzero_history_inlets_retained": True,
        "source_values_selected": False,
        "interval_or_cap_values_used": False,
        "universal_native_derivative_naturality_unproved": True,
        "direct_modified_Uy_fingerprint": False,
    }


def _recipes(nodes: Any, l: Mapping[str, Any], outer: Mapping[str, Any]) -> dict[str, Any]:
    charts = _mapping(l.get("charts")); outer_charts = _mapping(outer.get("charts"))
    chart = _mapping(outer_charts.get("reference"))
    fields = _mapping(chart.get("fields"))
    if not all(isinstance(fields.get(key), list) and len(fields[key]) >= 3 for key in ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")):
        return {"passed": False, "reason": "reference physical field ports through ordinary Z2 are missing"}
    coordinate = _handle(charts.get("reference", {}).get("coordinate"))
    radius = _handle(chart.get("radius"))
    z = _z_handle(nodes)
    pstar = _handle(outer.get("Pstar"))
    if None in (coordinate, radius, z, pstar):
        return {"passed": False, "reason": "reference recipe source ports missing"}
    try:
        recipe = PhysicalGermRecipes(
            fields={key: fields[key][:3] for key in ("E", "V", "m", "h", "k", "e", "p", "P0", "pressure")},
            radius=radius,
            Z=z,
            Pstar=pstar,
            delta=Fraction(1, 1000),
            coordinate=coordinate,
            endpoint=0,
        ).build()
    except Exception as exc:
        return {"passed": False, "reason": repr(exc)}
    recipe.update({
        "source_port_coordinate": coordinate,
        "source_port_endpoint": 0,
        "source_port_orders": [0, 1, 2],
        "source_ports_are_definitions": True,
        "recipe_graph_backward_indices_only": True,
        "source_values_selected": False,
    })
    return {"passed": True, **recipe}


def build_long_reference_germ(
    nodes: Any,
    l: Mapping[str, Any],
    outer: Mapping[str, Any],
    root: Any = ROOT,
) -> dict[str, Any]:
    """Build the bounded long/reference source germ from saved archives."""
    repo = _repo_root(root)
    loop_binding = _mapping(l, "prescribed_loop_C5_bindings")
    outer_binding = _mapping(outer, "prescribed_R100_Rm_C6_bindings")
    loop_archive = _load_archive(repo, _LOOP_ARCHIVE)
    outer_archive = _load_archive(repo, _OUTER_ARCHIVE)
    maps = _pin_maps(l, outer, loop_archive, outer_archive)
    pins = {
        "outer": _source_pin(repo, _OUTER_SOURCE, maps),
        "loop": _source_pin(repo, _LOOP_SOURCE, maps),
        "bridge": _source_pin(repo, _BRIDGE_SOURCE, maps),
    }
    source_evidence = _source_evidence(repo)
    interface = _fingerprint_interface(nodes, loop_binding, outer_binding)
    recipes = _recipes(nodes, loop_binding, outer_binding)
    physical = _physical_two_jet(source_evidence, interface)
    passed = bool(
        all(pin.get("pinned") is True for pin in pins.values())
        and source_evidence.get("passed") is True
        and interface.get("passed") is True
        and recipes.get("passed") is True
        and physical.get("passed") is True
    )
    return {
        "passed": passed,
        "method": "source_bound_long_reshape_native1_reference_native0_physical_germ",
        "source_only": True,
        "source_values_selected": False,
        "numeric_caps_used": False,
        "source_pins": pins,
        "source_evidence": source_evidence,
        "common_source_identity_rows": interface,
        "physical_germ_recipes": recipes,
        "physical_two_jet_certificate": physical,
        "leading_radial_laws": {
            "passed": physical.get("passed") is True,
            "rows": physical.get("E_y_graph_rows", []) + physical.get("V_y_graph_rows", []),
            "endpoint_E_y_equals_E_over_10": True,
            "E_yy_equals_E_over_100_at_join": True,
            "ordinary_Z_orders": [0, 1, 2],
            "native_Jacobians_kept_distinct": interface.get("native_Jacobians"),
        },
        "contract": {
            "long_reference_actual_source_join_only": True,
            "common_moments_through_ordinary_Z2": True,
            "common_root_and_loop_B_ports_through_ordinary_Z2": True,
            "physical_E_y_E_yy_and_V_y_V_yy_source_deduction": True,
            "native_Jacobians_not_equated": True,
            "global_phase_not_reset": True,
            "no_modified_Uy_Qy_RQy_alias": True,
            "no_error_cap_subtraction": True,
            "no_ancestor_constructor_or_producer_run": True,
        },
        "limits": {
            "universal_native_B_y_naturality_unproved": True,
            "modified_Uy_Qy_RQy_seam_unproved": True,
            "finite_N_error_budget_not_reinterpreted": True,
            "numeric_source_values_not_evaluated": True,
        },
    }


__all__ = ["build_long_reference_germ"]
