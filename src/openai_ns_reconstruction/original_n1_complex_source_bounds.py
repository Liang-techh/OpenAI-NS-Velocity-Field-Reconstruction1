"""Directed source/operator bounds for the original n=1 complex core.

This is a source-only evaluator.  It reads the accepted n=0 analytic-core
receipts and the literal n=1 raw-system AST; it never imports or runs an
owner, checker, recurrence, fixture, or quadrature.  The X_h coefficient
embedding is evaluated afresh with ``MPIntervalContext`` at ``|rho| <= 9/2``.

The module stops at genuine source and weighted A0/A1 bounds.  It does not
claim the n=1 solution norm or a radial Taylor remainder.  Physical Lambda is
kept out of the arithmetic: the accepted ``Lambda*epsilon=1`` identity is
used to reduce the pressure-head terms before interval evaluation.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Mapping

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext


PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_radial_IBP.json"
PARENT_SHA256 = "67945b4e9a8f2b1cdaa9a255b42a3bd9e5377f514cb12409ba41e9f17f957af0"
LATEST_PARENT_NAME = "lei_ren_part1_paper_compliant_current_original_C5_seam_joined_history.json"
LATEST_PARENT_SHA256 = "35fdc828c112e5fb09bf28a3f709b1d8a9817b54cd3a99fc11f76cfd92251d16"
RHO_RADIUS = "9/2"
RHO_TARGET = "4"
RAW_SOURCE = "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"


def _read(path: Path) -> Any:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as stream:
            return json.load(stream, parse_int=str)
    return json.loads(path.read_text(encoding="utf-8"), parse_int=str)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _mapping_digest(value: Any) -> str:
    """Digest JSON mappings after normalizing parser-dependent integers."""

    def normalize(item: Any) -> Any:
        if isinstance(item, bool) or item is None or isinstance(item, str):
            return item
        if isinstance(item, int):
            return str(item)
        if isinstance(item, float):
            return repr(item)
        if isinstance(item, Mapping):
            return {str(key): normalize(item[key]) for key in sorted(item, key=str)}
        if isinstance(item, (list, tuple)):
            return [normalize(part) for part in item]
        return repr(item)

    payload = json.dumps(normalize(value), sort_keys=True, separators=(",", ":"), ensure_ascii=False)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _endpoints(value: Any) -> tuple[mp.mpf, mp.mpf]:
    return tuple(mp.make_mpf(v) for v in value._mpi_)


def _int_tuple(value: tuple[mp.mpf, mp.mpf]) -> list[list[int]]:
    return [[int(part) for part in endpoint._mpf_] for endpoint in value]


def _interval(ctx: MPIntervalContext, value: Any) -> Any:
    """Restore a saved directed interval without converting through Decimal."""

    if hasattr(value, "_mpi_"):
        return value
    if isinstance(value, Mapping):
        lo = value.get("lower_exact_mpf_tuple")
        hi = value.get("upper_exact_mpf_tuple")
        if lo is not None and hi is not None:
            left = mp.make_mpf(tuple(int(part) for part in lo))
            right = mp.make_mpf(tuple(int(part) for part in hi))
            return ctx.mpf([left, right])
        if "lower" in value and "upper" in value:
            return ctx.mpf([str(value["lower"]), str(value["upper"])])
        if "upper" in value:
            return ctx.mpf(str(value["upper"]))
    return ctx.mpf(str(value))


def _lower_scalar(ctx: MPIntervalContext, value: Any) -> Any:
    lo, _ = _endpoints(_interval(ctx, value))
    return ctx.mpf(str(lo))


def _upper_scalar(ctx: MPIntervalContext, value: Any) -> Any:
    _, hi = _endpoints(_interval(ctx, value))
    return ctx.mpf(str(hi))


def _record(value: Any) -> dict[str, Any]:
    lo, hi = _endpoints(value)
    return {
        "lower_exact_mpf_tuple": list(map(int, lo._mpf_)),
        "upper_exact_mpf_tuple": list(map(int, hi._mpf_)),
        "lower_decimal_display": mp.nstr(lo, 28),
        "upper_decimal_display": mp.nstr(hi, 28),
        "directed_interval": True,
    }


def _serialization_roundtrip(ctx: MPIntervalContext, value: Any) -> int:
    """Ensure every new binary endpoint record is usable by its consumer."""
    if isinstance(value, Mapping):
        if value.get("directed_interval") is True:
            restored = _interval(ctx, value)
            expected = tuple(tuple(int(v) for v in value[key]) for key in
                ("lower_exact_mpf_tuple", "upper_exact_mpf_tuple"))
            if restored._mpi_ != expected:
                raise ValueError("Directed interval serialization did not round trip exactly")
            return 1
        return sum(_serialization_roundtrip(ctx, v) for v in value.values())
    if isinstance(value, (list, tuple)):
        return sum(_serialization_roundtrip(ctx, v) for v in value)
    return 0


def _upper_abs(ctx: MPIntervalContext, value: Any) -> Any:
    lo, hi = _endpoints(value)
    if lo <= 0 <= hi:
        return ctx.mpf([0, max(-lo, hi)])
    return ctx.mpf([0, max(abs(lo), abs(hi))])


def _max_interval(values: list[Any]) -> Any:
    if not values:
        raise ValueError("nonempty interval maximum required")
    return max(values, key=lambda value: _endpoints(value)[1])


def _sum(ctx: MPIntervalContext, values: list[Any]) -> Any:
    result = ctx.mpf(0)
    for value in values:
        result += value
    return result


def _family(report: Mapping[str, Any]) -> Mapping[str, Any] | None:
    value = report.get("source_family")
    return value if isinstance(value, Mapping) else None


def _delta_admission(domain: Mapping[str, Any], tube: Mapping[str, Any]) -> dict[str, Any]:
    """Bind the actual positive delta recipe before using the tube cap.

    The saved tube uses 1e-200 as a coefficient cap.  The actual pressure
    graph supplies ``min(1e-200, exp(-4 logPstar-30))`` with
    ``logPstar=exp(40)+11``.  The elementary inequalities below prove that
    the exponential branch is below the cap without materializing either
    the physical scale or the tiny delta.
    """

    analytic = (
        domain.get("original_pressure", {})
        .get("inherited_exact_P0_binding", {})
        .get("analytic_definition", {})
    )
    family = tube.get("core_family_definition", {})
    tube_delta = tube.get("delta", {})
    formula_ok = analytic.get("delta") == "min(1e-200,exp(-4logPstar-30))"
    logpstar_ok = analytic.get("logPstar") == "exp(Md)+11"
    md_ok = analytic.get("Md") == "40"
    cap_ok = family.get("delta_upper") == "1e-200"
    tube_delta_ok = (
        isinstance(tube_delta, Mapping)
        and tube_delta.get("lower") == "1.0e-200"
        and tube_delta.get("upper") == "1.0e-200"
    )
    L_lower = domain.get("denominator_bounds", {}).get("L_modulus_lower")
    saved_L_recorded = isinstance(L_lower, Mapping)
    saved_L_display_is_one = saved_L_recorded and L_lower.get("lower") == "1.0"
    return {
        "passed": formula_ok and logpstar_ok and md_ok and cap_ok and tube_delta_ok,
        "actual_delta_formula": analytic.get("delta"),
        "actual_logPstar_formula": analytic.get("logPstar"),
        "actual_Md": analytic.get("Md"),
        "tube_delta_cap": family.get("delta_upper"),
        "tube_delta_record": tube_delta,
        "actual_delta_positive": formula_ok and logpstar_ok and md_ok,
        "actual_delta_le_tube_cap": formula_ok and cap_ok and tube_delta_ok,
        "actual_delta_exponential_branch_proved": formula_ok and logpstar_ok and md_ok,
        "saved_L_modulus_lower": L_lower,
        "saved_L_lower_recorded": saved_L_recorded,
        "saved_L_display_is_one": saved_L_display_is_one,
        "saved_L_used_as_quantitative_lower": False,
        "elementary_proof": {
            "exp40_lower": "exp(40) >= 1+40+40^2/2 = 841",
            "logPstar_lower": "logPstar = exp(40)+11 >= 852",
            "negative_exponent_lower": "4*logPstar+30 >= 3438",
            "e3_lower": "e^3 > 1+3+9/2+27/6 = 13 > 10",
            "cap_comparison": "exp(-3438) < 10^-200",
            "conclusion": "0 < actual_delta = exp(-4*logPstar-30) < 1e-200",
        },
    }


def _fresh_L_lower(ctx: MPIntervalContext, domain: Mapping[str, Any], tube: Mapping[str, Any]) -> Any:
    """Derive |1-delta Z^2| from the admitted delta cap and Z tube."""

    outer = _interval(ctx, domain["common_tubes"]["outer_Z_radius"])
    delta_cap = _interval(ctx, tube["delta"])
    z_abs_upper = 1 + _upper_scalar(ctx, outer)
    delta_upper = _upper_scalar(ctx, delta_cap)
    # |1-delta Z^2| >= 1-delta_cap*|Z|^2, with |Z| <= 1+outer.
    return _lower_scalar(ctx, ctx.mpf(1) - delta_upper * z_abs_upper * z_abs_upper)


def _pin_files(root: Path, parent: Mapping[str, Any], paths: Mapping[str, Path]) -> dict[str, Any]:
    """Check actual bytes against the accepted parent input_hashes map."""

    pins = parent.get("input_hashes")
    checked: dict[str, Any] = {}
    missing: list[str] = []
    mismatches: list[str] = []
    if not isinstance(pins, Mapping):
        return {"passed": False, "checked": {}, "missing": ["parent.input_hashes"], "mismatches": []}
    for key, path in paths.items():
        candidates = (path.name, key, path.as_posix())
        selected = next((candidate for candidate in candidates if candidate in pins), None)
        if selected is None:
            missing.append(f"parent.input_hashes[{path.name}]")
            continue
        expected = pins[selected]
        if isinstance(expected, Mapping):
            expected = expected.get("sha256")
        actual = _sha256(path)
        ok = expected == actual
        checked[key] = {"parent_key": selected, "expected": expected, "actual": actual, "match": ok}
        if not ok:
            mismatches.append(f"{path.name}: parent hash differs")
    return {"passed": not missing and not mismatches, "checked": checked, "missing": missing, "mismatches": mismatches}


def _ast_binding(path: Path) -> dict[str, Any]:
    """Inspect the actual raw matrix literals without executing their module."""

    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "raw_matrix_expressions"]
    if len(functions) != 1:
        return {"passed": False, "reason": "unique raw_matrix_expressions AST missing"}
    node = functions[0]
    dictionaries: dict[str, ast.Dict] = {}
    for child in ast.walk(node):
        if isinstance(child, ast.Assign) and isinstance(child.value, ast.Dict):
            for target in child.targets:
                if isinstance(target, ast.Name) and target.id in ("B0", "B1"):
                    dictionaries[target.id] = child.value
    expected_kinematic = {"1,5", "2,6", "3,6"}
    expected_b0 = {"4,1", "5,1", "5,2", "5,3", "5,5", "6,1", "6,2", "6,3", "6,4", "6,6"}
    expected_b1 = {"5,1", "5,2", "5,3", "6,2", "6,3", "6,4"}

    def keys(table: ast.Dict) -> set[str]:
        return {key.value for key in table.keys if isinstance(key, ast.Constant) and isinstance(key.value, str)}

    b0 = keys(dictionaries["B0"]) if "B0" in dictionaries else set()
    b1 = keys(dictionaries["B1"]) if "B1" in dictionaries else set()
    normalized = "".join(ast.unparse(node).split())
    expected_literals = {
        ("B0", "1,5"): "sy.Integer(1)",
        ("B0", "2,6"): "sy.Integer(1)",
        ("B0", "3,6"): "sy.Integer(-1)",
        ("B0", "4,1"): "eps*4*x*f",
        ("B0", "5,1"): "eps*2*(q+E*(-2+de)/L)",
        ("B0", "5,2"): "eps*(2*Zf-2*(-1+de)*z*Cf/L)",
        ("B0", "5,3"): "-eps*2*(1+de)*z*Cf/L",
        ("B0", "5,5"): "eps*x*B",
        ("B0", "6,1"): "-eps**2*8*z*r*f/L",
        ("B0", "6,2"): "eps*2*(E*(-1+de)/L+Zw-wr*z*(-1+de)*r/L)",
        ("B0", "6,3"): "-eps*2*(1+de)*z*r*wr/L",
        ("B0", "6,4"): "-eps*4*z/L",
        ("B0", "6,6"): "eps*x*B",
        ("B1", "5,1"): "eps*2*H/L",
        ("B1", "5,2"): "-eps*2*d*Cf/L",
        ("B1", "5,3"): "-eps*2*d*Cf/L",
        ("B1", "6,2"): "eps*2*(H-wr*d*r)/L",
        ("B1", "6,3"): "-eps*2*d*r*wr/L",
        ("B1", "6,4"): "eps*2*d/L",
    }
    literal_matches: dict[str, bool] = {}
    for (table, key), expected in expected_literals.items():
        table_node = dictionaries.get(table)
        actual = None
        if table_node is not None:
            for literal_key, value in zip(table_node.keys, table_node.values):
                if isinstance(literal_key, ast.Constant) and literal_key.value == key:
                    actual = "".join(ast.unparse(value).split())
                    break
        literal_matches[f"{table}:{key}"] = actual == expected

    def degree(node: ast.AST) -> tuple[int, int]:
        """Minimum x-degree and parity for the polynomial AST (denominators are rho-free)."""

        if isinstance(node, ast.Name):
            return ((1, 1) if node.id == "x" else (2, 0) if node.id == "r" else (0, 0))
        if isinstance(node, ast.Constant) or isinstance(node, ast.Call):
            return (0, 0)
        if isinstance(node, ast.UnaryOp):
            return degree(node.operand)
        if isinstance(node, ast.BinOp):
            left, lp = degree(node.left)
            right, rp = degree(node.right)
            if isinstance(node.op, (ast.Add, ast.Sub)):
                return (min(left, right), lp) if lp == rp else (min(left, right), 0)
            if isinstance(node.op, ast.Mult):
                return (left + right, (lp + rp) % 2)
            if isinstance(node.op, ast.Div):
                return (left, lp) if right == 0 and rp == 0 else (0, 1)
            if isinstance(node.op, ast.Pow) and isinstance(node.right, ast.Constant) and isinstance(node.right.value, int):
                return (left * node.right.value, (lp * node.right.value) % 2)
        return (0, 1)

    axis_weights = {1: 2, 2: 2, 3: 2, 4: 2, 5: 1, 6: 1}
    source_powers = {4: 1, 5: 0, 6: 0}
    weighted_degree_matches: dict[str, bool] = {}
    for table, table_node in dictionaries.items():
        for literal_key, value in zip(table_node.keys, table_node.values):
            if not isinstance(literal_key, ast.Constant) or not isinstance(literal_key.value, str):
                continue
            i, j = (int(part) for part in literal_key.value.split(","))
            if i not in source_powers:
                continue
            raw_degree, parity = degree(value)
            shift = axis_weights[j] - source_powers[i]
            weighted_degree_matches[f"{table}:{literal_key.value}"] = raw_degree + shift >= 2 and (parity + shift) % 2 == 0
    raw_system_nodes = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "raw_system"]
    raw_system_normalized = "".join(ast.unparse(raw_system_nodes[0]).split()) if len(raw_system_nodes) == 1 else ""
    raw_system_fragments = (
        "a.add(a.scale('-1/2',S),a.scale(parameters['Lambda'],a.add(a.scale(2,c1(Q[1])),a.scale(r,c1(Q[2])))))",
        "g.mul(c(2),eps,x,Np.value)",
        "g.mul(c(2),eps,Nt)",
        "g.sub(g.mul(c(2),eps,Nz),g.quotient(g.mul(c(4),zv,r,eps,eps,Np.value),Lv,'positiveoriginalL'))",
    )
    raw_system_matches = {fragment: fragment in raw_system_normalized for fragment in raw_system_fragments}
    fragments = (
        "eps*4*x*f",
        "eps**2*8*z*r*f/L",
        "eps*4*z/L",
        "eps*2*H/L",
        "eps*2*d/L",
    )
    fragments_passed = all(fragment in normalized for fragment in fragments)
    return {
        "passed": b0 == expected_kinematic | expected_b0 and b1 == expected_b1 and fragments_passed and all(literal_matches.values()) and all(weighted_degree_matches.values()) and len(weighted_degree_matches) == 16 and all(raw_system_matches.values()),
        "B0_keys": sorted(b0),
        "B1_keys": sorted(b1),
        "off_kinematic_count": len(b0 - expected_kinematic) + len(b1),
        "expected_off_kinematic_count": 16,
        "fragments": {fragment: fragment in normalized for fragment in fragments},
        "literal_expression_match": literal_matches,
        "weighted_axis_degree_match": weighted_degree_matches,
        "weighted_axis_degree_proof": all(weighted_degree_matches.values()) and len(weighted_degree_matches) == 16,
        "raw_system_fragment_match": raw_system_matches,
        "raw_system_source_binding": all(raw_system_matches.values()),
        "source_executed": False,
    }


def _core_field_bounds(
    ctx: MPIntervalContext,
    domain: Mapping[str, Any],
    transfer: Mapping[str, Any],
    tube: Mapping[str, Any],
    radius: Any,
) -> dict[str, Any]:
    """Bind Phi/Psi X_h norms to actual complex-rho mixed-derivative bounds."""

    joint = domain["joint_core_domain"]
    bphi = _interval(ctx, joint["Xh_norm_majorants"]["Phi"])
    bpsi = _interval(ctx, joint["Xh_norm_majorants"]["Psi"])
    eta = _interval(ctx, domain["common_tubes"]["outer_Z_radius"])
    z_abs = 1 + _upper_scalar(ctx, eta)
    j = _interval(ctx, tube["j"])
    epsilon = _interval(ctx, transfer["epsilon"])
    f0 = _interval(ctx, domain["original_axis_amplitude"]["complex_modulus_upper"])
    delta = _interval(ctx, tube["delta"])
    matrix_loss = _lower_scalar(ctx, domain["common_tubes"]["matrix_Z_radius"])
    gg_loss = _lower_scalar(ctx, domain["common_tubes"]["Gg_Z_radius"])

    # The X_h norm controls the joint Taylor series, not only the axial
    # centre jets.  On |rho|<5 and |Z-centre|<h/2 the coefficient majorant is
    # the binomial series 1/(1-|rho|/20-|Z-centre|/h).  Its boundary value is
    # 1/(1-1/4-1/2)=4.  Cauchy in rho with the admitted radius 1/2 then gives
    # the honest sup bound 4 * B * i! * 2**i for the i-th radial derivative on
    # |rho|<=9/2.  This is the actual X_h definition, rather than the
    # ``core_value_factor`` metadata alone.
    factor4 = ctx.mpf(4)
    f_radial = [f0 * factor4 * bphi * math.factorial(i) * 2**i for i in range(3)]
    u_radial = [epsilon * factor4 * bpsi * math.factorial(i) * 2**i for i in range(3)]
    u_radial[0] += 4 * z_abs + j
    # The accepted domain exports the same positive radial envelope for U and
    # the own mean M; use that saved envelope for the Q recovery bound.
    m_radial = list(u_radial)
    fz = f_radial[0] / matrix_loss
    frz = f_radial[1] / matrix_loss
    fzz = 2 * f_radial[0] / (matrix_loss**2)
    uz = u_radial[0] / matrix_loss
    urz = u_radial[1] / matrix_loss
    uzz = 2 * u_radial[0] / (matrix_loss**2)
    return {
        "Bphi_Xh": bphi,
        "Bpsi_Xh": bpsi,
        "F": f_radial,
        "Fr": f_radial[1],
        "Frr": f_radial[2],
        "FZ": fz,
        "FrZ": frz,
        "FZZ": fzz,
        "U": u_radial,
        "Ur": u_radial[1],
        "Urr": u_radial[2],
        "UZ": uz,
        "UrZ": urz,
        "UZZ": uzz,
        "M": m_radial,
        "MZ": m_radial[0] / matrix_loss,
        "radius": radius,
        "z_abs": z_abs,
        "delta": delta,
        "epsilon": epsilon,
        "f0": f0,
        "matrix_loss": matrix_loss,
        "Gg_loss": gg_loss,
        "L_lower": _fresh_L_lower(ctx, domain, tube),
        "L_lower_formula": "1 - delta_cap*(1+outer_Z_radius)^2",
        "embedding_formula": "joint Xh binomial majorant; sup <= 4*B*i!*2^i on |rho|<=9/2, |Z-centre|<=eta",
        "joint_binomial_majorant": "sum_(n,m)>=0 binomial(n+m,m)(1/4)^n(1/2)^m = 4",
        "rho_cauchy_radius": "1/2",
        "uniform_complex_sup_factor4": True,
    }


def _source_bounds(ctx: MPIntervalContext, fields: Mapping[str, Any], domain: Mapping[str, Any]) -> dict[str, Any]:
    """Bound Q, Ntheta, Nz, S and the pressure-head source by positive majorants."""

    z = fields["z_abs"]
    delta = fields["delta"]
    L = fields["L_lower"]
    R = fields["radius"]
    d = 1 + z * z
    one_plus_delta = 1 + delta
    # The saved tube cap bounds the actual positive graph delta.  For the
    # coefficient 1-delta use its absolute majorant 1+|delta|; using
    # 1-delta here would be an underbound because the saved cap is not the
    # actual value.
    one_minus_delta_abs = 1 + delta
    U, Ur, Urr = fields["U"][0], fields["U"][1], fields["U"][2]
    M, Mr, Mrr = fields["M"][0], fields["M"][1], fields["M"][2]
    matrix_loss, gg_loss = fields["matrix_loss"], fields["Gg_loss"]
    Q = [
        (2 * z * U + one_plus_delta * z * M + d * M / matrix_loss) / L,
        (2 * z * Ur + one_plus_delta * z * Mr + d * Mr / matrix_loss) / L,
        (2 * z * Urr + one_plus_delta * z * Mrr + d * Mrr / matrix_loss) / L,
    ]
    QZ = Q[0] / gg_loss
    f, fr, frr = fields["F"][0], fields["F"][1], fields["F"][2]
    fz, frz, fzz = fields["FZ"], fields["FrZ"], fields["FZZ"]
    uz, urz, uzz = fields["UZ"], fields["UrZ"], fields["UZZ"]

    def second(rows: tuple[Any, Any, Any, Any, Any, Any], weight: Any) -> Any:
        value, radial, second_radial, zderiv, radial_zderiv, second_zderiv = rows
        A = abs(weight) * z * value + d * zderiv + 2 * z * R * radial
        Az = (
            abs(weight) * value
            + abs(weight - 2) * z * zderiv
            + d * second_zderiv
            + 2 * R * radial
            + 2 * z * R * radial_zderiv
        )
        Ar = abs(weight - 2) * z * radial + d * radial_zderiv + 2 * z * R * second_radial
        return (
            A * abs(weight - 1 + delta) * z / L
            + (Az / L + A * 2 * delta * z / (L * L)) * d / L
            + Ar * 2 * z * R / (L * L)
        )

    ntheta = second((f, fr, frr, fz, frz, fzz), -2 - delta)
    nz = second((U, Ur, Urr, uz, urz, uzz), -1 - delta)
    S = (
        (Q[0] + Q[0] * one_minus_delta_abs * z / (2 * gg_loss) + Q[1] * R) / L
        + Q[0] * (Q[0] / 2 + Q[1] * R)
        + U * (2 * z * Q[0] + d * QZ + 2 * z * R * Q[1]) / L
    )
    pressure_head = 2 * Q[1] + R * Q[2]
    eps_np = fields["epsilon"] * S / 2 + pressure_head
    g4_over_x = 2 * eps_np
    g5 = 2 * fields["epsilon"] * ntheta
    g6 = 2 * fields["epsilon"] * nz + 4 * z * R / L * fields["epsilon"] * eps_np
    return {
        "Q": Q,
        "QZ": QZ,
        "S": S,
        "N1theta": ntheta,
        "N1z": nz,
        "pressure_head_without_Lambda": pressure_head,
        "epsilon_Np": eps_np,
        "g4_over_x": g4_over_x,
        "g5": g5,
        "g6": g6,
        "formula": "eps*Np <= eps*S/2 + 2*Qr + rho*Qrr by Lambda*epsilon=1",
    }


def _weighted_operator_bounds(ctx: MPIntervalContext, fields: Mapping[str, Any], source: Mapping[str, Any], domain: Mapping[str, Any]) -> dict[str, Any]:
    """Evaluate the exact weighted A0/A1 row majorants from the raw AST."""

    z, delta, L, R, eps = (
        fields["z_abs"],
        fields["delta"],
        fields["L_lower"],
        fields["radius"],
        fields["epsilon"],
    )
    d = 1 + z * z
    f, fr, fz, w, wr, wz, q = (
        fields["F"][0], fields["F"][1], fields["FZ"], fields["U"][0], fields["U"][1], fields["UZ"], source["Q"][0]
    )
    Cf = f + R * fr
    Zf = ((2 + delta) * z * f + d * fz + 2 * z * R * fr) / L
    Zw = ((1 + delta) * z * w + d * wz + 2 * z * R * wr) / L
    H = z * (1 + delta) / 2 + d * w
    B = q + (1 + 2 * z * w) / L
    E = z * w + ctx.mpf(".5")
    a0 = {
        "4,1": eps * 4 * f,
        "5,1": eps * 2 * (q + E * (2 + delta) / L),
        "5,2": eps * (2 * Zf + 2 * (1 + delta) * z * Cf / L),
        "5,3": eps * (2 * (1 + delta) * z * Cf / L),
        "5,5": eps * B,
        # After x^(p_1-q_6)=rho and division by rho, this entry is still the
        # raw rho factor from B0[6,1].
        "6,1": eps * eps * 8 * z * R * f / L,
        "6,2": eps * 2 * ((1 + delta) * E / L + Zw + (1 + delta) * z * R * wr / L),
        "6,3": eps * 2 * (1 + delta) * z * R * wr / L,
        "6,4": eps * 4 * z / L,
        "6,6": eps * B,
    }
    a1 = {
        "5,1": eps * 2 * H / L,
        "5,2": eps * 2 * d * Cf / L,
        "5,3": eps * 2 * d * Cf / L,
        "6,2": eps * 2 * (H + wr * d * R) / L,
        "6,3": eps * 2 * d * R * wr / L,
        "6,4": eps * 2 * d / L,
    }
    rows0 = {"4": _sum(ctx, [a0["4,1"]]), "5": _sum(ctx, [a0["5,1"], a0["5,2"], a0["5,3"], a0["5,5"]]), "6": _sum(ctx, [a0["6,1"], a0["6,2"], a0["6,3"], a0["6,4"], a0["6,6"]])}
    rows1 = {"5": _sum(ctx, [a1["5,1"], a1["5,2"], a1["5,3"]]), "6": _sum(ctx, [a1["6,2"], a1["6,3"], a1["6,4"]])}
    return {"A0_entries": a0, "A1_entries": a1, "A0_rows": rows0, "A1_rows": rows1, "C0": _max_interval(list(rows0.values())), "C1": _max_interval(list(rows1.values()))}


def build(root: str | Path, parent: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Build a directed source/operator receipt without running producers."""

    repo = Path(root)
    base = repo / "experiments" / "root_st073"
    legacy_parent_path = base / PARENT_NAME
    latest_parent_path = base / LATEST_PARENT_NAME
    parent_candidates = {
        LATEST_PARENT_NAME: (latest_parent_path, LATEST_PARENT_SHA256),
        PARENT_NAME: (legacy_parent_path, PARENT_SHA256),
    }
    parent_mapping_match = True
    if parent is None:
        selected = next(
            ((name, path, digest) for name, (path, digest) in parent_candidates.items() if path.exists() and _sha256(path) == digest),
            None,
        )
        if selected is None:
            selected = (PARENT_NAME, legacy_parent_path, PARENT_SHA256)
        selected_parent_name, selected_parent_path, selected_parent_sha256 = selected
        parent_report = _read(selected_parent_path) if selected_parent_path.exists() else {}
    else:
        parent_report = parent
        selected = next(
            (
                (name, path, digest)
                for name, (path, digest) in parent_candidates.items()
                if path.exists() and _sha256(path) == digest and _mapping_digest(parent_report) == _mapping_digest(_read(path))
            ),
            None,
        )
        if selected is None:
            selected_parent_name, selected_parent_path, selected_parent_sha256 = PARENT_NAME, legacy_parent_path, PARENT_SHA256
            parent_mapping_match = False
        else:
            selected_parent_name, selected_parent_path, selected_parent_sha256 = selected
    names = {
        "parent": selected_parent_name,
        "domain": "lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain.json.gz",
        "domain_check": "lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain_check.json",
        "transfer": "lei_ren_part1_paper_compliant_core_transfer.json",
        "common_check": "lei_ren_part1_paper_compliant_current_core_common_fixed_point_check.json",
        "tube": "lei_ren_part1_paper_shared_analytic_tube.json",
        "volterra": "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.json.gz",
        "volterra_check": "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions_check.json",
        "operator": "lei_ren_part1_paper_compliant_current_original_n1_inner_operator.json.gz",
        "operator_check": "lei_ren_part1_paper_compliant_current_original_n1_inner_operator_check.json",
        "raw_source": RAW_SOURCE,
    }
    paths = {key: base / name for key, name in names.items()}
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        return {"passed": False, "actual_source_bounds_proved": False, "diagnostic_only": True, "missing_input_contract": missing}
    parent_hash_ok = _sha256(paths["parent"]) == selected_parent_sha256 and parent_mapping_match
    parent_file_mapping_digest = _mapping_digest(_read(paths["parent"]))
    # The selected parent report is checked by its immutable file pin above;
    # its own bytes are intentionally not required to occur in its input map.
    pin_paths = {key: path for key, path in paths.items() if key != "parent"}
    pin = _pin_files(repo, parent_report, pin_paths)
    records = {key: _read(path) for key, path in paths.items() if key not in ("raw_source",)}
    radial = records["parent"]
    domain, domain_check = records["domain"], records["domain_check"]
    transfer, common, tube = records["transfer"], records["common_check"], records["tube"]
    volterra_check, operator_check = records["volterra_check"], records["operator_check"]
    delta_admission = _delta_admission(domain, tube)
    family = _family(domain)
    parent_family = _family(radial)
    family_ok = isinstance(family, Mapping) and parent_family == family
    for report in (volterra_check, operator_check):
        family_ok = family_ok and _family(report) == family
    for report in (transfer, common):
        family_ok = family_ok and all(report.get(key) == family.get(key) for key in ("implicit_source_sha256", "datum_enclosure_sha256"))
    source_identity = radial.get("n1_source_vector_identity", {})
    source_ok = source_identity.get("passed") is True and source_identity.get("symbolic", {}).get("checks", {}).get("shared_core_Lambda_epsilon_reciprocal") is True
    core_transfer_ok = transfer.get("contraction_proved") is True and transfer.get("analytic_core_exists_for_compliant_source") is True
    xh = domain.get("joint_core_domain", {}).get("Xh_norm_majorants", {})
    # The accepted transfer stores a directed interval [1, upper], while the
    # domain receipt records the selected upper endpoint as a point interval.
    # Bind the domain majorant to that exact transfer upper endpoint rather
    # than comparing the two differently shaped JSON records.
    def same_upper(field: str, transfer_field: str) -> bool:
        domain_value = xh.get(field)
        transfer_value = transfer.get(transfer_field)
        if not isinstance(domain_value, Mapping) or not isinstance(transfer_value, Mapping):
            return False
        return domain_value.get("upper_exact_mpf_tuple") == transfer_value.get("upper_exact_mpf_tuple")

    xh_match = same_upper("Phi", "Phi_ball_norm_upper") and same_upper("Psi", "Psi_ball_norm_upper")
    core_rows = domain.get("actual_core_radial_bounds", {})
    mean_envelope_match = core_rows.get("Uz") == core_rows.get("mean") and isinstance(core_rows.get("Uz"), list)
    axis_binding = domain.get("original_axis_amplitude", {})
    f0_binding_ok = (
        axis_binding.get("source_formula") == "exp(-selected_logCstar-Lambda*G(Z))"
        and axis_binding.get("no_source_amplitude_materialized") is True
    )
    tube_radii = domain.get("common_tubes", {})
    check_ctx = MPIntervalContext()
    check_ctx.dps = 80
    fresh_L_lower_check = _fresh_L_lower(check_ctx, domain, tube)
    fresh_L_lower_positive = _lower_scalar(check_ctx, fresh_L_lower_check) > 0
    nested_complex_Z_tubes = False
    if all(name in tube_radii for name in ("outer_Z_radius", "matrix_Z_radius", "Gg_Z_radius", "solution_Z_radius")):
        nested_complex_Z_tubes = (
            _lower_scalar(check_ctx, tube_radii["outer_Z_radius"])
            >= _upper_scalar(check_ctx, tube_radii["matrix_Z_radius"])
            and _lower_scalar(check_ctx, tube_radii["matrix_Z_radius"])
            >= _upper_scalar(check_ctx, tube_radii["Gg_Z_radius"])
            and _lower_scalar(check_ctx, tube_radii["Gg_Z_radius"])
            >= _upper_scalar(check_ctx, tube_radii["solution_Z_radius"])
        )
    outer_within_Xh_disk = False
    outer_plus_matrix_within_Xh_disk = False
    matrix_at_least_two_Gg = False
    if "outer_Z_radius" in tube_radii and "Z_disk_radius" in domain.get("joint_core_domain", {}):
        outer_within_Xh_disk = _upper_scalar(check_ctx, tube_radii["outer_Z_radius"]) <= _lower_scalar(check_ctx, domain["joint_core_domain"]["Z_disk_radius"])
    if all(name in tube_radii for name in ("outer_Z_radius", "matrix_Z_radius", "Gg_Z_radius")) and "Z_disk_radius" in domain.get("joint_core_domain", {}):
        outer_plus_matrix_within_Xh_disk = (
            _upper_scalar(check_ctx, _interval(check_ctx, tube_radii["outer_Z_radius"]) + _interval(check_ctx, tube_radii["matrix_Z_radius"]))
            <= _lower_scalar(check_ctx, domain["joint_core_domain"]["Z_disk_radius"])
        )
        matrix_at_least_two_Gg = _lower_scalar(check_ctx, tube_radii["matrix_Z_radius"]) >= 2 * _upper_scalar(check_ctx, tube_radii["Gg_Z_radius"])
    proof = common.get("actual_Banach_and_common_coefficient_proof", {})
    xh_embedding = (
        tube.get("norm_definition") == "paper8.30:20^n h^m(n+1)^2(m+1)^2|d_Z^m f_n|/[m! binomial(n+m,m)]"
        and xh_match
        and proof.get("actual_source_operator_defined_on_same_complete_Xh_pair") is True
        and common.get("same_solution_all_order_tail_proof", {}).get("same_actual_fixed_point_coefficient_norm_used") is True
        and domain_check.get("all_passed") is True
        and core_transfer_ok
        and domain.get("joint_core_domain", {}).get("rho_complex_modulus_less_than") == "5"
        and mean_envelope_match
        and f0_binding_ok
        and nested_complex_Z_tubes
        and outer_within_Xh_disk
        and outer_plus_matrix_within_Xh_disk
        and matrix_at_least_two_Gg
        and delta_admission["passed"]
    )
    domain_ok = (
        xh_embedding
        and domain.get("source_bindings", {}).get("actual_core_collar_is_not_analytic_core_extension") is True
        and domain.get("source_inner_common_holomorphic_bounds_constructed") is True
        and domain.get("all_fixed_radial_derivatives_holomorphic_on_same_domain_by_smooth_ODE") is True
        and fresh_L_lower_positive
        and delta_admission["passed"]
    )
    ast_receipt = _ast_binding(paths["raw_source"])
    structural = (
        parent_hash_ok
        and pin["passed"]
        and family_ok
        and source_ok
        and core_transfer_ok
        and domain_check.get("all_passed") is True
        and domain_ok
        and ast_receipt.get("passed") is True
        and volterra_check.get("all_passed") is True
        and operator_check.get("all_passed") is True
        and volterra_check.get("current_original_n1_true_inner_Volterra_function_adapter_installed") is True
        and operator_check.get("current_original_n1_actual_inner_operator_input_enclosures_installed") is True
    )
    if not structural:
        return {
            "passed": False,
            "actual_source_bounds_proved": False,
            "diagnostic_only": True,
            "parent_report_sha256": _sha256(paths["parent"]),
            "parent_mapping_digest": _mapping_digest(parent_report),
            "parent_file_mapping_digest": parent_file_mapping_digest,
            "parent_mapping_matches_file": parent_mapping_match,
            "parent_pins": pin,
            "source_family_match": family_ok,
            "xh_embedding_proved": xh_embedding,
            "delta_admission": delta_admission,
            "fresh_L_lower_positive": fresh_L_lower_positive,
            "fresh_L_lower": _record(fresh_L_lower_check),
            "nested_complex_Z_tubes": nested_complex_Z_tubes,
            "outer_plus_matrix_within_Xh_disk": outer_plus_matrix_within_Xh_disk,
            "matrix_at_least_two_Gg": matrix_at_least_two_Gg,
            "raw_ast_binding": ast_receipt,
            "missing_input_contract": [
                "accepted parent pins and source-family agreement",
                "actual Xh complex-rho embedding hypotheses",
                "accepted raw n1 source/operator receipts",
            ],
        }

    ctx = MPIntervalContext()
    ctx.dps = 220
    with mp.workdps(ctx.dps + 40):
        R = ctx.mpf(9) / 2
        fields = _core_field_bounds(ctx, domain, transfer, tube, R)
        source = _source_bounds(ctx, fields, domain)
        operator = _weighted_operator_bounds(ctx, fields, source, domain)
        # The word lemma consumes the actual Cauchy gap from the Gg tube to
        # the solution tube.  Do not substitute a radius unless the saved
        # nested domains prove that it is this gap.
        gg_radius = _interval(ctx, domain["common_tubes"]["Gg_Z_radius"])
        solution_radius = _interval(ctx, domain["common_tubes"]["solution_Z_radius"])
        sigma_gap = gg_radius - solution_radius
        sigma = _lower_scalar(ctx, sigma_gap)
        msrc = _max_interval([source["g4_over_x"] / 2, source["g5"] / 4, source["g6"] / 2])
        beta = ctx.e * R / 2 * (operator["C0"] + operator["C1"] / sigma)
        beta_ok = _endpoints(beta)[1] <= mp.mpf("0.5")
        sigma_gap_ok = _endpoints(sigma_gap)[0] > 0
        missing_contract: list[str] = []
        if not sigma_gap_ok:
            missing_contract.append("positive Gg-to-solution Z Cauchy gap")
        if not beta_ok:
            missing_contract.append("directed beta <= 1/2 at R=9/2")
        result = {
            "passed": not missing_contract,
            "actual_source_bounds_proved": True,
            "actual_weighted_operator_bounds_proved": True,
            "actual_solution_norm_proved": False,
            "actual_radial_Taylor_remainder_certified": False,
            "diagnostic_only": False,
            "source_only": True,
            "source_executed": False,
            "accepted_owner_checker_fixture_recurrence_quadrature_replayed": False,
            "parent_report": selected_parent_name,
            "parent_report_sha256": _sha256(paths["parent"]),
            "parent_mapping_digest": _mapping_digest(parent_report),
            "parent_file_mapping_digest": parent_file_mapping_digest,
            "parent_mapping_matches_file": parent_mapping_match,
            "parent_pins": pin,
            "source_family_match": family_ok,
            "raw_ast_binding": ast_receipt,
            "delta_admission": delta_admission,
            "fresh_L_lower_positive": fresh_L_lower_positive,
            "domain": {
                "rho_source_disk": "|rho| < 5",
                "rho_bound_radius": RHO_RADIUS,
                "rho_target_radius": RHO_TARGET,
                "outer_Z_radius": _record(_interval(ctx, domain["common_tubes"]["outer_Z_radius"])),
                "matrix_Z_radius": _record(_interval(ctx, domain["common_tubes"]["matrix_Z_radius"])),
                "Gg_Z_radius": _record(_interval(ctx, domain["common_tubes"]["Gg_Z_radius"])),
                "solution_Z_radius": _record(_interval(ctx, domain["common_tubes"]["solution_Z_radius"])),
                "L_modulus_lower_saved": _record(_interval(ctx, domain["denominator_bounds"]["L_modulus_lower"])),
                "L_modulus_lower_from_delta_cap": _record(fields["L_lower"]),
                "L_modulus_lower_formula": fields["L_lower_formula"],
                "actual_collar_used_as_core_extension": False,
            },
            "xh_embedding": {
                "proved": xh_embedding,
                "norm_definition": tube["norm_definition"],
                "Phi_Xh_norm": _record(fields["Bphi_Xh"]),
                "Psi_Xh_norm": _record(fields["Bpsi_Xh"]),
                "embedding_formula": fields["embedding_formula"],
                "radial_geometric_ratio": "R/20=9/40<1",
                "joint_boundary_rho_ratio": "5/20=1/4",
                "uniform_complex_rho_sup": True,
                "factor4_from_joint_binomial_series": True,
                "mean_envelope_same_as_Uz": mean_envelope_match,
                "F0_complex_axis_binding": f0_binding_ok,
                "nested_complex_Z_tubes": nested_complex_Z_tubes,
                "outer_tube_within_Xh_Z_disk": outer_within_Xh_disk,
                "outer_plus_matrix_within_Xh_Z_disk": outer_plus_matrix_within_Xh_disk,
                "matrix_at_least_two_Gg": matrix_at_least_two_Gg,
            },
            "n0_fields": {key: _record(value) for key, value in {"F": fields["F"][0], "Fr": fields["F"][1], "Frr": fields["F"][2], "FZ": fields["FZ"], "FrZ": fields["FrZ"], "FZZ": fields["FZZ"], "U": fields["U"][0], "Ur": fields["U"][1], "Urr": fields["U"][2], "UZ": fields["UZ"], "UrZ": fields["UrZ"], "UZZ": fields["UZZ"], "M": fields["M"][0], "Mr": fields["M"][1], "Mrr": fields["M"][2], "MZ": fields["MZ"]}.items()},
            "source_bounds": {key: (_record(value) if hasattr(value, "_mpi_") else [_record(v) for v in value]) for key, value in source.items() if key not in ("formula",)},
            "weighted_A0_A1": {
                "A0_entries": {key: _record(value) for key, value in operator["A0_entries"].items()},
                "A1_entries": {key: _record(value) for key, value in operator["A1_entries"].items()},
                "A0_rows": {key: _record(value) for key, value in operator["A0_rows"].items()},
                "A1_rows": {key: _record(value) for key, value in operator["A1_rows"].items()},
                "C0": _record(operator["C0"]),
                "C1": _record(operator["C1"]),
                "weighted_unknown": "(F1/x^2,Uz1/x^2,K1/x^2,P1/x^2,F1x/x,Uz1x/x)",
                "normalized_source": "(g4/x,g5,g6)",
            },
            "green_input": {
                "Msrc": _record(msrc),
                "sigma": _record(sigma),
                "sigma_available_gap": _record(sigma_gap),
                "sigma_definition": "Gg_Z_radius - solution_Z_radius",
                "beta": _record(beta),
                "beta_le_half": beta_ok,
                "sigma_gap_positive": sigma_gap_ok,
                "word_bound": "Msrc*(R/2)^k*(C0+k*C1/sigma)^k/k!",
                "geometric_beta": "e*R/2*(C0+C1/sigma)",
                "physical_Lambda_materialized": False,
            },
            "gates": {
                "parent_hashes": pin["passed"],
                "source_family": family_ok,
                "source_normalization": source_ok,
                "core_transfer": core_transfer_ok,
                "analytic_domain_receipt": domain_check.get("all_passed") is True,
                "xh_complex_rho_embedding": xh_embedding,
                "mean_envelope_binding": mean_envelope_match,
                "F0_complex_axis_binding": f0_binding_ok,
                "nested_complex_Z_tubes": nested_complex_Z_tubes,
                "outer_tube_within_Xh_Z_disk": outer_within_Xh_disk,
                "outer_plus_matrix_within_Xh_Z_disk": outer_plus_matrix_within_Xh_disk,
                "matrix_at_least_two_Gg": matrix_at_least_two_Gg,
                "delta_cap_admission": delta_admission["passed"],
                "actual_L_lower_from_delta_cap": fresh_L_lower_positive,
                "L_holomorphic_lower_positive": domain_ok,
                "actual_raw_source_AST": ast_receipt.get("passed") is True,
                "directed_interval_arithmetic": True,
                "sigma_gap_positive": sigma_gap_ok,
                "beta_le_half": beta_ok,
                "solution_norm": False,
                "radial_Taylor_remainder": False,
            },
            "missing_input_contract": missing_contract,
            "physical_scale": "ordinary d_R^i=exp(i*logLambda)d_rho^i; no Lambda or physical derivative is materialized",
        }
        result["exact_interval_records_roundtripped"] = _serialization_roundtrip(ctx, result)
        result["gates"]["interval_serialization_roundtrip"] = result["exact_interval_records_roundtripped"] > 0
        return result


__all__ = ["build"]
