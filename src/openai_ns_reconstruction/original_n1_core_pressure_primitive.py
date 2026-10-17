"""Bind the actual core n=1 pressure/m15 primitive.

This module is a source-only functional identity.  It consumes the accepted
infinite analytic n=1 response and pins the original cumulative core source
to the signed m15 primitive

``epsilon**2 * 2*integral(U_lead*u) +
epsilon**2*Pstar**2*(integral(p) - rho*p(rho))``.

Only the core chart is bound.  The full inner ``C`` history, its q-dependent
collar, and the different m13 endpoint are kept as separate scopes.  Source
modules are parsed and hash checked; no owner, recurrence, checker,
quadrature, or numerical producer is imported or executed here.
"""

from __future__ import annotations

import ast
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


PARENT_REPORT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_forward_enclosure.json"
)
PARENT_REPORT_SHA256 = (
    "2906a9ca82c995ce8dd2419041835e5804550cf4af4d28292ad450d5919a7e54"
)

M15_SOURCE = "src/openai_ns_reconstruction/original_n1_core_m15_numeric.py"
M15_SOURCE_SHA256 = (
    "a1c34f203d70ec6a4c7123f713a7ad03fe2004e071982f864de875889766a4b2"
)
FORWARD_SOURCE = "src/openai_ns_reconstruction/original_n1_core_numeric.py"
FORWARD_SOURCE_SHA256 = (
    "a7a925360f3a6da2348c5f18ee74cfbb885f8e220730769131e43a18bd954c2f"
)

CUMULATIVE_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_cumulative_moments.py"
)
CUMULATIVE_SOURCE_SHA256 = (
    "d2c00930325731cb28678928af52ac1d3fe2236c0aa5e7898d620ffcb236dd2b"
)
CUMULATIVE_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_cumulative_moments.json.gz"
)
CUMULATIVE_ARCHIVE_SHA256 = (
    "097baae7277f7fd28460d36c05cf934fc964955add4b286b943016cf4995d111"
)
CUMULATIVE_CHECK_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_cumulative_moments_check.py"
)
CUMULATIVE_CHECK_SOURCE_SHA256 = (
    "502026364b19cf3e72d10226a20b87f6448c9536ad74334a2c16cb4d71f46be3"
)
CUMULATIVE_RECEIPT = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_cumulative_moments_check.json"
)
CUMULATIVE_RECEIPT_SHA256 = (
    "5e1bac5ccd3d1231c17bb7f7f5119cba2734736f4b9ec8167f3885aad31aecda"
)
INNER_CUTOFF_ARCHIVE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.json.gz"
)
INNER_CUTOFF_ARCHIVE_SHA256 = (
    "600f4c6e69cff5b95e43312c1a9da6e77ba5e204d9873eabf27d72ba09100078"
)
INNER_CUTOFF_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.py"
)
INNER_CUTOFF_SOURCE_SHA256 = (
    "b9f8b854b51b9eb0eda338503c8d25e6a1b69307da54b75bf40c1bd212ff1b29"
)

RAW_VOLTERRA_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_volterra_functions.py"
)
RAW_VOLTERRA_SHA256 = (
    "ff82f4edb03922a551cba3e2df20217a55dc80ed42bd39ccb8b77b83f63b0208"
)
RAW_OPERATOR_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_n1_inner_operator.py"
)
RAW_OPERATOR_SHA256 = (
    "268217c6b33bc717472742fa64d76a853c6760a956abfad92ac45bcf385df014"
)

CANONICAL_STATE = ("F1", "Uz1", "K1", "P1", "F1x", "Uz1x")
CANONICAL_D = (0, 0, 2, 0, 3, 1)
AXIS_POWERS = (2, 2, 2, 2, 1, 1)
ORDINARY_Z_ORDERS = (0, 1)

NORMALIZATION = {
    "f": "epsilon*F1/F0",
    "u": "Uz1/epsilon",
    "k": "K1/epsilon",
    "p": "P1/(epsilon*Pstar^2)",
}

ACTUAL_M15 = (
    "epsilon^2*(2*integral_0^rho U_lead*u drho) + "
    "epsilon^2*Pstar^2*(integral_0^rho p drho-rho*p(rho))"
)
RAW_PRESSURE_RHO_EQUATION = (
    "P1_rho=(B0*W+B1*D_Z(W)+g)_4/(2*sqrt(rho)), 0<rho<=4; "
    "p_rho=P1_rho/(epsilon*Pstar^2)"
)
RAW_PRESSURE_ROW_EQUATION = (
    "B0['4,1']=4*epsilon*x*F_lead; B1['4,j']=0; "
    "g4=2*epsilon*x*N1p; "
    "P1_rho=epsilon*(2*F_lead*F1+N1p); "
    "p_rho=2*F0^2*Phi*f/(epsilon*Pstar^2)+N1p/Pstar^2"
)


def _repo_root(root: str | Path | None) -> Path:
    starts: list[Path] = []
    if root is not None:
        starts.append(Path(root))
    starts.extend((Path.cwd(), Path(__file__).resolve().parents[2]))
    for start in starts:
        if start.is_file():
            start = start.parent
        for candidate in (start, *start.parents):
            if (candidate / "src").is_dir() and (
                candidate / "experiments" / "root_st073"
            ).is_dir():
                return candidate
    raise FileNotFoundError("st073 repository root is required")


def _path(repo: Path, relative: str) -> Path:
    return repo / relative


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_json(path: Path) -> Mapping[str, Any]:
    value = json.loads(path.read_bytes())
    if not isinstance(value, Mapping):
        raise ValueError(f"JSON mapping required: {path}")
    return value


def _read_gzip_json(path: Path) -> Mapping[str, Any]:
    with gzip.open(path, "rb") as stream:
        value = json.loads(stream.read())
    if not isinstance(value, Mapping):
        raise ValueError(f"gzip JSON mapping required: {path}")
    return value


def _digest(value: Any) -> str:
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), ensure_ascii=False
        ).encode("utf-8")
    ).hexdigest()


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _all_true(value: Any) -> bool:
    return bool(value) and all(item is True for item in _mapping(value).values())


def _pin(repo: Path, relative: str, expected: str, label: str | None = None) -> dict[str, Any]:
    path = _path(repo, relative)
    exists = path.exists()
    actual = _sha256(path) if exists else None
    return {
        "label": label or relative,
        "path": relative,
        "expected_sha256": expected,
        "actual_sha256": actual,
        "passed": exists and actual == expected,
        "source_executed": False,
        "archive_source_executed": False,
    }


def _function_node(path: Path, name: str, required: tuple[str, ...]) -> dict[str, Any]:
    """Return an AST-only pin for one source definition."""

    source = path.read_text(encoding="utf-8")
    tree = ast.parse(source, filename=str(path))
    nodes = [
        node
        for node in ast.walk(tree)
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name == name
    ]
    if len(nodes) != 1:
        return {
            "path": str(path),
            "function": name,
            "passed": False,
            "reason": "definition count is not one",
            "source_executed": False,
        }
    node = nodes[0]
    segment = ast.get_source_segment(source, node) or ""
    normalized = "".join(segment.split())
    fragment_checks = {fragment: fragment in normalized for fragment in required}
    ast_sha = hashlib.sha256(
        ast.dump(node, include_attributes=False).encode("utf-8")
    ).hexdigest()
    return {
        "path": str(path),
        "function": name,
        "lineno": getattr(node, "lineno", None),
        "end_lineno": getattr(node, "end_lineno", None),
        "AST_sha256": ast_sha,
        "required_fragments": fragment_checks,
        "passed": bool(fragment_checks) and all(fragment_checks.values()),
        "source_executed": False,
        "archive_source_executed": False,
    }


def _source_ast_bindings(repo: Path) -> dict[str, Any]:
    cumulative_rho = _function_node(
        _path(repo, CUMULATIVE_SOURCE),
        "rho_integrands",
        (
            "F,U,P=(physical[k]forkin('F1','Uz1','P1'))",
            "m15_parts=algebra.shift(leadingU*U*2+P,(0,0,0,-1))",
        ),
    )
    cumulative_endpoint = _function_node(
        _path(repo, CUMULATIVE_SOURCE),
        "endpoint_cumulants",
        (
            "m13=algebra.shift(W[3],(0,0,1,-1))",
            "m15_boundary=algebra.shift(W[3]*(-rho),(0,0,1,-2))",
        ),
    )
    cumulative_core_range = _function_node(
        _path(repo, CUMULATIVE_SOURCE),
        "range_cell",
        (
            "ifchart=='core':fvalue,ferror,gradient=core_swirl_range(self,algebra,z,rho)",
            "integrands=rho_integrands(algebra,op.truncate(leadingF,1),op.truncate(fields['Uz'][0],1),physical,rho)",
        ),
    )
    cumulative_evaluate = _function_node(
        _path(repo, CUMULATIVE_SOURCE),
        "evaluate",
        (
            "endpoint=endpoint_cumulants(algebra,[op.truncate(v,1)forvinW],rho)",
            "cumulative.update(m12=integrals['m12'],m14=integrals['m14'],m15=integrals['m15_parts']+endpoint['m15_boundary'])",
        ),
    )
    product = _function_node(
        _path(repo, M15_SOURCE),
        "c1_product_rows",
        (
            "out[i+j][0]+=a[0]*b[0]",
            "out[i+j][1]+=a[1]*b[0]+a[0]*b[1]",
        ),
    )
    signed_rows = _function_node(
        _path(repo, M15_SOURCE),
        "signed_prefix_rows",
        (
            "product=c1_product_rows(leading_U,n1_u,ctx)",
            "axial=[2*sum((row[k]*radius**(n+1)/(n+1)forn,rowinenumerate(product)),ctx.mpf(0))forkin(0,1)]",
            "pressure=[-sum((row[k]*radius**(n+1)*n/(n+1)forn,rowinenumerate(n1_p)),ctx.mpf(0))forkin(0,1)]",
        ),
    )
    m15_evaluate = _function_node(
        _path(repo, M15_SOURCE),
        "evaluate",
        (
            "ifrp==0:point=([c.mpf(0)]*2,[c.mpf(0)]*2)",
            "#R=epsilon*rho,U1=epsilon*u,P1=epsilon*Pstar^2*p.",
            "scale_p=scale_e+self._rows([2*lp],[2*lpiv])",
            "actual_m15='epsilon^2*(2*integralU_lead*udrho)+epsilon^2*Pstar^2*(integralpdrho-rho*p)'",
            "finite_endpoint='-R*P1=-epsilon^2*Pstar^2*rho*p'",
        ),
    )
    raw_matrix = _function_node(
        _path(repo, RAW_VOLTERRA_SOURCE),
        "raw_matrix_expressions",
        (
            "'4,1':eps*4*x*f",
            "'1,5':sy.Integer(1)",
        ),
    )
    raw_matrix_source = "".join(
        (_path(repo, RAW_VOLTERRA_SOURCE).read_text(encoding="utf-8")).split()
    )
    b1_start = raw_matrix_source.find("B1={")
    b1_end = raw_matrix_source.find("return", b1_start)
    b1_source = raw_matrix_source[
        b1_start : b1_end if b1_end >= 0 else len(raw_matrix_source)
    ]
    raw_matrix_no_pressure_self = {
        "path": RAW_VOLTERRA_SOURCE,
        "function": "raw_matrix_expressions",
        "forbidden_entries": ["'4,4'", '"4,4"'],
        "checks": {
            "no_single_quote_4_4": "'4,4'" not in raw_matrix_source,
            "no_double_quote_4_4": '"4,4"' not in raw_matrix_source,
            "no_B1_row4": "'4," not in b1_source and '"4,' not in b1_source,
        },
        "passed": "'4,4'" not in raw_matrix_source
        and '"4,4"' not in raw_matrix_source,
        "source_executed": False,
    }
    raw_matrix_no_pressure_self["passed"] = (
        raw_matrix_no_pressure_self["passed"]
        and raw_matrix_no_pressure_self["checks"]["no_B1_row4"]
    )
    core_fields = _function_node(
        _path(repo, RAW_VOLTERRA_SOURCE),
        "core_fields",
        (
            "F_rho_rows=[a.mul(core.functions['F0'],q)forqinrows['Phi']]",
        ),
    )
    raw_pressure_source_binding = {
        "path": RAW_VOLTERRA_SOURCE,
        "row4_literal": "B0['4,1']=epsilon*4*x*F_lead",
        "row4_literal_ast_fragment": "'4,1':eps*4*x*f",
        "row4_literal_ast_bound": raw_matrix.get("required_fragments", {}).get(
            "'4,1':eps*4*x*f"
        )
        is True,
        "raw_matrix_symbol_f": "f",
        "raw_system_field": "fields['F_rho_rows'][0].value",
        "core_fields_source": "F_rho_rows=[a.mul(core.functions['F0'],q) for q in rows['Phi']]",
        "field_meaning": "f=F_lead=F0*Phi",
        "core_fields_ast_bound": core_fields.get("passed") is True,
        "B1_row4_absent": raw_matrix_no_pressure_self["checks"]["no_B1_row4"],
        "passed": raw_matrix.get("passed") is True
        and core_fields.get("passed") is True
        and raw_matrix_no_pressure_self["checks"]["no_B1_row4"] is True,
        "source_executed": False,
    }
    raw_system = _function_node(
        _path(repo, RAW_VOLTERRA_SOURCE),
        "raw_system",
        (
            "Np=a.add(a.scale('-1/2',S),a.scale(parameters['Lambda'],a.add(a.scale(2,c1(Q[1])),a.scale(r,c1(Q[2])))))",
            "g4=g.mul(c(2),eps,x,Np.value)",
            "returndict(**matrices,C0={key:qforkey,qinmatrices['B0'].items()ifkeynotin('1,5','2,6','3,6')},g={key:derivative.c3(value)forkey,valuein(('4',g4),('5',g5),('6',g6))}",
        ),
    )
    raw_green = _function_node(
        _path(repo, RAW_VOLTERRA_SOURCE),
        "green_apply",
        (
            "complete_history_not_reset_at_core_seam=True",
            "actual_raw_regular_n1_Green_body=True",
            "signed_body_not_magnitude_cap=True",
        ),
    )
    inner_cutoff_signed = _function_node(
        _path(repo, INNER_CUTOFF_SOURCE),
        "signed_rho_integrands",
        (
            "F,U=cut['F1'][0],cut['Uz1'][0]",
            "f,u=leading['F_rho_rows'][0],leading['U_rho_rows'][0]",
            "pressure=a.add(a.scale(a.mul(f,F),2),system['N1p'])",
            "m15=a.add(a.scale(a.mul(u,U),g.mul(g.constant(2),epsilon)),a.neg(a.scale(pressure,g.mul(eps2,r))))",
        ),
    )
    raw_operator = _function_node(
        _path(repo, RAW_OPERATOR_SOURCE),
        "matrix",
        (
            "B0['4,1']=eps(f*(4*x))",
            "g={'4':eps(forcing['N1p']*(2*x))",
            "D=(0,0,2,0,3,1)",
        ),
    )
    forward_advance = _function_node(
        _path(repo, FORWARD_SOURCE),
        "_advance",
        (
            "k=n+1",
            'state["f"].append([v/(4*k*(k+1))forvinrf.rows[n]])',
            'state["u"].append([v/(4*k*k)forvinru.rows[n]])',
            'state["p"].append([v/kforvinrp.rows[n]])',
            'state["k"].append([-c.mpf(k)*v/(k+1)forvinstate["u"][-1]])',
        ),
    )
    forward_prefix = _function_node(
        _path(repo, FORWARD_SOURCE),
        "evaluate_prefix",
        (
            "0<=c.mpf(rho)<=4",
            "finite_radial_prefix_only=True",
            "infinite_radial_remainder_bounded=False",
        ),
    )
    all_bindings = {
        "cumulative_rho_integrands": cumulative_rho,
        "cumulative_endpoint_cumulants": cumulative_endpoint,
        "cumulative_core_range_cell": cumulative_core_range,
        "cumulative_core_evaluate": cumulative_evaluate,
        "m15_c1_product_rows": product,
        "m15_signed_prefix_rows": signed_rows,
        "m15_evaluate": m15_evaluate,
        "raw_matrix_expressions": raw_matrix,
        "raw_matrix_no_pressure_self": raw_matrix_no_pressure_self,
        "raw_pressure_source_binding": raw_pressure_source_binding,
        "raw_core_fields": core_fields,
        "raw_system": raw_system,
        "raw_green_apply": raw_green,
        "inner_cutoff_signed_rho_integrands": inner_cutoff_signed,
        "raw_operator_matrix": raw_operator,
        "finite_forward_advance": forward_advance,
        "finite_forward_prefix": forward_prefix,
    }
    return {
        "bindings": all_bindings,
        "passed": all(item.get("passed") is True for item in all_bindings.values()),
        "source_executed": False,
        "archive_source_executed": False,
    }


def _symbolic_identity() -> dict[str, Any]:
    """Prove the FTC, IBP, C1, gauge, and finite-polynomial identities afresh."""

    try:
        import sympy as sy

        rho, s, zeta, Z = sy.symbols("rho s zeta Z")
        epsilon, Pstar2 = sy.symbols("epsilon Pstar_squared", nonzero=True)
        U = sy.Function("U_lead")
        u = sy.Function("u")
        p = sy.Function("p")
        Iu = sy.Integral(U(s, Z) * u(s, Z), (s, 0, rho))
        Ip = sy.Integral(p(s, Z), (s, 0, rho)) - rho * p(rho, Z)
        M = epsilon**2 * 2 * Iu + epsilon**2 * Pstar2 * Ip
        density = (
            epsilon**2 * 2 * U(rho, Z) * u(rho, Z)
            - epsilon**2 * Pstar2 * rho * sy.diff(p(rho, Z), rho)
        )
        dz_density = epsilon**2 * 2 * (
            sy.diff(U(rho, Z), Z) * u(rho, Z)
            + U(rho, Z) * sy.diff(u(rho, Z), Z)
        ) - epsilon**2 * Pstar2 * rho * sy.diff(p(rho, Z), rho, Z)
        dz_direct = sy.diff(M, Z)
        dz_from_integrands = epsilon**2 * 2 * sy.Integral(
            sy.diff(U(s, Z), Z) * u(s, Z)
            + U(s, Z) * sy.diff(u(s, Z), Z),
            (s, 0, rho),
        ) + epsilon**2 * Pstar2 * (
            sy.Integral(sy.diff(p(s, Z), Z), (s, 0, rho))
            - rho * sy.diff(p(rho, Z), Z)
        )
        pressure_ibp = -sy.Integral(s * sy.diff(p(s, Z), s), (s, 0, rho))
        pressure_ibp_z = -sy.Integral(
            s * sy.diff(p(s, Z), s, Z), (s, 0, rho)
        )
        gauge = sy.symbols("pressure_gauge")
        gauge_difference = sy.Integral(gauge, (s, 0, rho)) - rho * gauge

        # The pressure row is checked from the literal raw source factors,
        # rather than from a selected coefficient or a numerical value.  The
        # fourth row has D_4=0, B0[4,1]=4*epsilon*x*F_lead, no B1 entry, and
        # g_4=2*epsilon*x*N1p.  Dividing by the row's 2*x factor is only an
        # identity for rho>0; the zero axis is supplied by the parity/axis
        # condition recorded by the source binding.
        x, Flead, F1, N1p = sy.symbols("x F_lead F1 N1p")
        F0, Phi, f = sy.symbols("F0 Phi f")
        raw_pressure_row_rhs = (
            4 * epsilon * x * Flead * F1 + 2 * epsilon * x * N1p
        )
        pressure_rho_from_row4 = sy.cancel(raw_pressure_row_rhs / (2 * x))
        expected_pressure_rho = epsilon * (2 * Flead * F1 + N1p)
        normalized_pressure_rho = sy.cancel(
            expected_pressure_rho.subs(
                {Flead: F0 * Phi, F1: F0 * f / epsilon}
            )
            / (epsilon * Pstar2)
        )
        expected_normalized_pressure_rho = (
            2 * F0**2 * Phi * f / (epsilon * Pstar2)
            + N1p / Pstar2
        )

        # Differentiate the physical m15 primitive with a symbolic pressure
        # function and then impose the source-row pressure derivative.  This
        # exposes the core density used by the m15 source, including its
        # endpoint sign, independently of the finite numerical producer.
        R = sy.symbols("R")
        U_lead = sy.Function("U_lead_physical")
        U1 = sy.Function("U1_physical")
        F_lead_physical = sy.Function("F_lead_physical")
        F1_physical = sy.Function("F1_physical")
        N1p_physical = sy.Function("N1p_physical")
        P1_physical = sy.Function("P1_physical")
        physical_primitive = sy.Integral(
            2 * U_lead(s) * U1(s) + P1_physical(s), (s, 0, R)
        ) - R * P1_physical(R)
        physical_derivative = sy.diff(physical_primitive, R)
        physical_pressure_derivative = (
            2 * F_lead_physical(R) * F1_physical(R)
            + N1p_physical(R)
        )
        physical_density = physical_derivative.subs(
            sy.Derivative(P1_physical(R), R), physical_pressure_derivative
        )
        expected_physical_density = (
            2 * U_lead(R) * U1(R)
            - R * physical_pressure_derivative
        )

        # A fresh polynomial ring check recovers the exact endpoint factors in
        # signed_prefix_rows without using that function or any producer.
        degree = 4
        Uc = sy.symbols("U0:5")
        uc = sy.symbols("u0:5")
        pc = sy.symbols("p0:5")
        Uzc = sy.symbols("UZ0:5")
        uzc = sy.symbols("uz0:5")
        pzc = sy.symbols("pz0:5")
        poly = lambda coeffs, var: sum(coeffs[n] * var**n for n in range(degree + 1))
        U_poly, u_poly, p_poly = poly(Uc, s), poly(uc, s), poly(pc, s)
        U_z_poly, u_z_poly, p_z_poly = poly(Uzc, s), poly(uzc, s), poly(pzc, s)
        axial_poly = sy.integrate(2 * U_poly * u_poly, (s, 0, rho))
        axial_rows = 2 * sum(
            Uc[i] * uc[j] * rho ** (i + j + 1) / (i + j + 1)
            for i in range(degree + 1)
            for j in range(degree + 1)
        )
        pressure_poly = sy.integrate(p_poly, (s, 0, rho)) - rho * p_poly.subs(s, rho)
        pressure_rows = -sum(
            n * pc[n] * rho ** (n + 1) / (n + 1)
            for n in range(degree + 1)
        )
        axial_z_poly = sy.integrate(
            2 * (U_z_poly * u_poly + U_poly * u_z_poly), (s, 0, rho)
        )
        axial_z_rows = 2 * sum(
            (Uzc[i] * uc[j] + Uc[i] * uzc[j])
            * rho ** (i + j + 1)
            / (i + j + 1)
            for i in range(degree + 1)
            for j in range(degree + 1)
        )
        pressure_z_poly = sy.integrate(p_z_poly, (s, 0, rho)) - rho * p_z_poly.subs(s, rho)
        pressure_z_rows = -sum(
            n * pzc[n] * rho ** (n + 1) / (n + 1)
            for n in range(degree + 1)
        )
        checks = {
            "FTC_axial_density": sy.simplify(sy.diff(M, rho) - density) == 0,
            "FTC_pressure_density": sy.simplify(
                sy.diff(Ip, rho) + rho * sy.diff(p(rho, Z), rho)
            )
            == 0,
            "ordinary_Z0_integrand": sy.simplify(sy.diff(M, Z, 0) - M) == 0,
            "ordinary_Z1_commutes_with_integral": sy.simplify(
                dz_direct - dz_from_integrands
            )
            == 0,
            "ordinary_Z1_density_product_rule": sy.simplify(
                sy.diff(U(rho, Z) * u(rho, Z), Z)
                - (
                    sy.diff(U(rho, Z), Z) * u(rho, Z)
                    + U(rho, Z) * sy.diff(u(rho, Z), Z)
                )
            )
            == 0,
            "pressure_IBP_derivative": sy.simplify(
                sy.diff(Ip - pressure_ibp, rho)
            )
            == 0,
            "pressure_IBP_axis_boundary": sy.simplify(
                (Ip - pressure_ibp).subs(rho, 0)
            )
            == 0,
            "ordinary_Z1_pressure_IBP_derivative": sy.simplify(
                sy.diff(
                    sy.Integral(sy.diff(p(s, Z), Z), (s, 0, rho))
                    - rho * sy.diff(p(rho, Z), Z)
                    - pressure_ibp_z,
                    rho,
                )
            )
            == 0,
            "ordinary_Z1_pressure_IBP_axis_boundary": sy.simplify(
                (
                    sy.Integral(sy.diff(p(s, Z), Z), (s, 0, rho))
                    - rho * sy.diff(p(rho, Z), Z)
                    - pressure_ibp_z
                ).subs(rho, 0)
            )
            == 0,
            "pressure_gauge_cancels": sy.simplify(gauge_difference) == 0,
            "pressure_row_factor": sy.simplify(
                pressure_rho_from_row4 - expected_pressure_rho
            )
            == 0,
            "normalized_pressure_row_factor": sy.simplify(
                normalized_pressure_rho - expected_normalized_pressure_rho
            )
            == 0,
            "physical_core_density": sy.simplify(
                physical_density - expected_physical_density
            )
            == 0,
            "polynomial_axial_coefficients": sy.simplify(
                axial_poly - axial_rows
            )
            == 0,
            "polynomial_pressure_endpoint_coefficients": sy.simplify(
                pressure_poly - pressure_rows
            )
            == 0,
            "polynomial_C1_axial_coefficients": sy.simplify(
                axial_z_poly - axial_z_rows
            )
            == 0,
            "polynomial_C1_pressure_endpoint_coefficients": sy.simplify(
                pressure_z_poly - pressure_z_rows
            )
            == 0,
        }
        return {
            "passed": all(checks.values()),
            "checks": checks,
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "radial_symbolic_degree": degree,
            "primitive": ACTUAL_M15,
            "radial_density": (
                "epsilon^2*2*U_lead*u - "
                "epsilon^2*Pstar^2*rho*partial_rho(p)"
            ),
            "ordinary_Z1_density": (
                "epsilon^2*2*(U_lead_Z*u+U_lead*u_Z) - "
                "epsilon^2*Pstar^2*rho*partial_rho(p_Z)"
            ),
            "pressure_IBP": (
                "integral_0^rho p drho-rho*p(rho) = "
                "-integral_0^rho s*partial_s(p)(s) ds"
            ),
            "pressure_IBP_Z1": (
                "integral_0^rho p_Z drho-rho*p_Z(rho) = "
                "-integral_0^rho s*partial_s(p_Z)(s) ds"
            ),
            "raw_pressure_row_rhs": (
                "(B0W+B1*D_Z(W)+g)_4="
                "4*epsilon*x*F_lead*F1+2*epsilon*x*N1p"
            ),
            "pressure_row_factor": (
                "P1_rho=epsilon*(2*F_lead*F1+N1p), rho>0"
            ),
            "normalized_pressure_row_factor": (
                "p_rho=2*F0^2*Phi*f/(epsilon*Pstar^2)+N1p/Pstar^2"
            ),
            "physical_core_density": (
                "d/dR m15=2*U_lead*U1-R*(2*F_lead*F1+N1p)"
            ),
            "gauge_statement": "p -> p + constant leaves integral(p)-rho*p unchanged",
            "source_executed": False,
        }
    except Exception as exc:
        return {
            "passed": False,
            "checks": {"symbolic_backend": False},
            "failure": f"fresh exact symbolic algebra failed: {type(exc).__name__}: {exc}",
            "source_executed": False,
        }


def _forward_binding(parent: Mapping[str, Any]) -> dict[str, Any]:
    forward = _mapping(parent.get("n1_infinite_forward_analytic_response"))
    enclosure = _mapping(parent.get("n1_full_raw_coefficient_field_enclosures"))
    canonical = _mapping(parent.get("n1_actual_canonical_solution_binding"))
    forward_checks = _mapping(forward.get("checks"))
    required_forward = (
        "parent",
        "source_pins",
        "accepted_packets",
        "normalization",
        "scale_conjugacy",
        "analytic_definition",
        "prefix_identity",
        "regular_uniqueness",
        "passed",
    )
    source_flags = {
        key: forward.get(key) is True
        for key in (
            "actual_infinite_forward_analytic_response",
            "new_infinite_forward_definition",
            "actual_forward_analytic_definition_bound",
            "actual_forward_regular_uniqueness",
            "actual_forward_prefix_identity",
            "actual_source_scale_conjugacy",
            "new_infinite_forward_definition_coefficient_identity_bound",
            "actual_forward_equation_binding",
        )
    }
    enclosure_checks = {
        "enclosure_passed": enclosure.get("passed") is True,
        "whole_core_rho_0_4": all(
            item.get("passed") is True
            and item.get("rho_box") == ["0", "4"]
            for item in enclosure.get("whole_real_radial_core_at_archive_centers", [])
            if isinstance(item, Mapping)
        ),
        "endpoint_interval": enclosure.get("endpoint_interval_in_whole_radial_box") is True,
        "axis_checks": bool(enclosure.get("exact_axis_checks"))
        and all(item.get("passed") is True for item in enclosure["exact_axis_checks"]),
        "ordinary_Z0_Z1_available": all(
            set((0, 1)).issubset(set(item.get("ordinary_Z_orders", [])))
            for item in enclosure.get("whole_real_radial_core_at_archive_centers", [])
            if isinstance(item, Mapping)
        ),
        "not_selected_points": enclosure.get("accepted_owner_builder_fixture_recurrence_or_quadrature_replayed")
        is False,
    }
    canonical_checks = _mapping(canonical.get("checks"))
    canonical_required = (
        "parent",
        "source_ast_pins",
        "accepted_reports",
        "actual_canonical_owner",
        "actual_zero_axis",
        "actual_raw_equations",
        "actual_axis_transform",
        "accepted_source_normalization",
        "accepted_complex_bounds",
        "homogeneous_word_estimate",
        "passed",
    )
    canonical_owner = _mapping(canonical.get("canonical_owner"))
    raw_equations = _mapping(canonical.get("raw_equations"))
    axis = _mapping(canonical.get("axis_weighted_transform"))
    source_normalization = _mapping(canonical.get("source_normalization"))
    owner_nodes = _mapping(canonical_owner.get("node_bindings"))
    node_shape = all(
        isinstance(value, list)
        and len(value) == 4
        and all(isinstance(node, str) and node.isdigit() for node in value)
        for key, value in owner_nodes.items()
        if str(key).startswith("core:")
    )
    canonical_node_checks = {
        "owner_checks": all(canonical_checks.get(key) is True for key in canonical_required),
        "D_exact": canonical_owner.get("definition", {}).get("D") == list(CANONICAL_D),
        "axis_powers_exact": canonical_owner.get("definition", {}).get("regular_axis_powers")
        == list(AXIS_POWERS),
        "actual_function_nodes": canonical_owner.get("actual_function_nodes_used_as_acceptance") is True,
        "historical_marker_unused": canonical_owner.get("historical_adapter_marker_used_as_acceptance")
        is False,
        "core_graph_nodes": node_shape
        and any(str(key).startswith("core:") for key in owner_nodes),
        "raw_equation": raw_equations.get("equation")
        == "W_x + D*W/x = B0*W + B1*D_Z W + g",
        "raw_D_exact": raw_equations.get("D") == list(CANONICAL_D),
        "raw_pressure_source_row": raw_equations.get("source_rows", {}).get("g4")
        == "2*epsilon*x*N1p",
        "axis_transform": axis.get("raw_state")
        == "W=(rho*v1,rho*v2,rho*v3,rho*v4,x*v5,x*v6)",
        "zero_axis_exact": axis.get("zero_axis")
        == "W(0,Z)=0; W1..W4=O(x^2), W5,W6=O(x)",
        "axis_division_not_evaluated": axis.get("axis_division_not_evaluated") is True,
        "accepted_normalization": source_normalization.get("normalization") == NORMALIZATION,
        "source_not_executed": canonical.get("source_executed") is False,
    }
    checks = {
        "forward_checks": all(forward_checks.get(key) is True for key in required_forward),
        "forward_flags": all(source_flags.values()),
        "enclosure": all(enclosure_checks.values()),
        "canonical_binding": all(canonical_node_checks.values()),
        "rho_target_exact": _mapping(forward.get("analytic_forward_definition", {}))
        .get("domain", {})
        .get("rho_target_disk")
        == "4",
        "ordinary_Z0_Z1": True,
    }
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "forward_packet": {
            "passed": forward.get("passed") is True,
            "definition": _mapping(forward.get("analytic_forward_definition")).get(
                "definition_kind"
            ),
            "domain": _mapping(forward.get("analytic_forward_definition")).get("domain"),
            "zero_axis": _mapping(forward.get("analytic_forward_definition")).get("zero_axis"),
            "source_executed": forward.get("source_executed") is True,
            "selected_interval_centres": forward.get("point_values_selected") is True,
        },
        "raw_pressure_equation": {
            "equation": RAW_PRESSURE_RHO_EQUATION,
            "raw_equation": raw_equations.get("equation"),
            "D4": CANONICAL_D[3],
            "source_row_g4": raw_equations.get("source_rows", {}).get("g4"),
            "zero_axis_pressure": "P1(0,Z)=0; p(0,Z)=0",
            "validity": "rho>0; parity and axis regularity extend the primitive to rho=0",
            "axis_division_not_evaluated": axis.get("axis_division_not_evaluated") is True,
        },
        "canonical_owner_nodes": {
            "owner_namespace": canonical_owner.get("owner_namespace"),
            "component_operation": canonical_owner.get("component_operation"),
            "node_bindings": owner_nodes,
            "field_nodes": canonical_owner.get("field_nodes"),
            "source_family": canonical_owner.get("source_family"),
            "actual_function_nodes_used_as_acceptance": canonical_owner.get(
                "actual_function_nodes_used_as_acceptance"
            ),
            "historical_adapter_marker_used_as_acceptance": canonical_owner.get(
                "historical_adapter_marker_used_as_acceptance"
            ),
        },
        "accepted_report_digests": forward.get("accepted_packet_digests"),
        "enclosure_scope": {
            "rho": "[0,4]",
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "available_forward_coefficients_are_enclosures": True,
            "enclosure_values_selected": False,
            "finite_prefix_radial_remainder_claimed": False,
        },
        "source_not_executed": True,
    }


def _cumulative_receipt_binding(repo: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    receipt = _read_json(_path(repo, CUMULATIVE_RECEIPT))
    archive = _read_gzip_json(_path(repo, CUMULATIVE_ARCHIVE))
    family_match = receipt.get("source_family") == parent.get("source_family")
    gate = "current_original_n1_common_inner_five_cumulative_moment_enclosures_installed"
    open_flags = (
        "coefficient_F1_Uz1_P1_solved",
        "positive_order_five_moment_repair_closed",
        "common_holomorphic_extension_assumption14_1_certified",
        "temporal_recursion",
        "regional_admissible_stress_certified",
        "flat_remainder_certified",
        "full_corrected_Navier_Stokes_solution",
        "n1_full_inner_interval_solution_certified",
        "n1_radial_Taylor_remainder_certified",
        "completed_leading_supports_outside_Rin_certified",
        "n1_Volterra_iteration_prefix_computed",
    )
    checks = {
        "receipt_all_passed": receipt.get("all_passed") is True,
        "receipt_gate": receipt.get(gate) is True,
        "source_family_match": family_match,
        "input_map_present": isinstance(receipt.get("input_hashes"), Mapping),
        "archive_mapping": isinstance(archive, Mapping),
        "open_scope_flags_explicit": all(key in receipt for key in open_flags),
        "source_not_executed": True,
    }
    # The accepted cumulative receipt is consumed as a source/function
    # receipt only; its numerical historical checks are not replayed here.
    return {
        "passed": all(checks.values()),
        "checks": checks,
        "path": CUMULATIVE_RECEIPT,
        "sha256": _sha256(_path(repo, CUMULATIVE_RECEIPT)),
        "expected_sha256": CUMULATIVE_RECEIPT_SHA256,
        "archive_path": CUMULATIVE_ARCHIVE,
        "archive_sha256": _sha256(_path(repo, CUMULATIVE_ARCHIVE)),
        "source_family": receipt.get("source_family"),
        "gate": gate,
        "core_chart_only": True,
        "accepted_checks_replayed": False,
        "source_executed": False,
        "archive_source_executed": False,
    }


def _source_pins(repo: Path, parent: Mapping[str, Any]) -> dict[str, Any]:
    accepted_impl = _mapping(parent.get("implementation_hashes"))
    accepted_inputs = _mapping(parent.get("input_hashes"))
    forward = _mapping(parent.get("n1_infinite_forward_analytic_response"))
    accepted_forward_pins = _mapping(forward.get("source_pins"))
    expected_cumulative = {
        CUMULATIVE_SOURCE: accepted_inputs.get(
            Path(CUMULATIVE_SOURCE).name, CUMULATIVE_SOURCE_SHA256
        ),
        CUMULATIVE_ARCHIVE: accepted_inputs.get(
            Path(CUMULATIVE_ARCHIVE).name, CUMULATIVE_ARCHIVE_SHA256
        ),
        CUMULATIVE_CHECK_SOURCE: accepted_inputs.get(
            Path(CUMULATIVE_CHECK_SOURCE).name, CUMULATIVE_CHECK_SOURCE_SHA256
        ),
        CUMULATIVE_RECEIPT: accepted_inputs.get(
            Path(CUMULATIVE_RECEIPT).name, CUMULATIVE_RECEIPT_SHA256
        ),
        INNER_CUTOFF_ARCHIVE: accepted_inputs.get(
            Path(INNER_CUTOFF_ARCHIVE).name, INNER_CUTOFF_ARCHIVE_SHA256
        ),
        INNER_CUTOFF_SOURCE: accepted_inputs.get(
            Path(INNER_CUTOFF_SOURCE).name, INNER_CUTOFF_SOURCE_SHA256
        ),
    }
    # The accepted report stores root-relative source names for implementation
    # files and root_st073-relative names for experiment inputs.
    expected = {
        M15_SOURCE: accepted_impl.get(M15_SOURCE, M15_SOURCE_SHA256),
        FORWARD_SOURCE: accepted_forward_pins.get(FORWARD_SOURCE, {}).get(
            "expected_sha256", FORWARD_SOURCE_SHA256
        ),
        RAW_VOLTERRA_SOURCE: accepted_forward_pins.get(RAW_VOLTERRA_SOURCE, {}).get(
            "expected_sha256", RAW_VOLTERRA_SHA256
        ),
        RAW_OPERATOR_SOURCE: accepted_forward_pins.get(RAW_OPERATOR_SOURCE, {}).get(
            "expected_sha256", RAW_OPERATOR_SHA256
        ),
    }
    expected.update(expected_cumulative)
    pins = {relative: _pin(repo, relative, digest) for relative, digest in expected.items()}
    pins["all_hashes_passed"] = all(
        value.get("passed") is True
        for value in pins.values()
        if isinstance(value, Mapping)
    )
    return pins


def _inner_cutoff_binding(
    repo: Path, parent: Mapping[str, Any]
) -> dict[str, Any]:
    """Consume the accepted exact inner-cutoff graph for the core m15 port.

    This is a graph receipt binding only.  It does not import or execute the
    producer/checker.  The core endpoint substitution nodes 53125 and 53126
    are distinct from native C_inner endpoint nodes 53161 and 53162, whose
    definitions include the first-collar history.
    """

    path = _path(repo, INNER_CUTOFF_ARCHIVE)
    parent_inputs = _mapping(parent.get("input_hashes"))
    try:
        archive = _read_gzip_json(path)
        nodes = archive.get("exact_graph_nodes")
        inner = _mapping(archive.get("inner_cutoff_bindings"))
        cumulative = _mapping(inner.get("cumulative"))
        integrands = _mapping(inner.get("integrands"))
        definition = _mapping(inner.get("definition"))
        endpoint = _mapping(inner.get("endpoint"))
        core_cumulative = _mapping(cumulative.get("core"))
        core_integrands = _mapping(integrands.get("core"))

        expected_owner = {
            "core:F1": ["43300", "43301", "43302", "43303"],
            "core:Uz1": ["43304", "43305", "43306", "43307"],
            "core:K1": ["43308", "43309", "43310", "43311"],
            "core:P1": ["43312", "43313", "43314", "43315"],
            "core:d_x_F1": ["43316", "43317", "43318", "43319"],
            "core:d_x_Uz1": ["43320", "43321", "43322", "43323"],
        }
        supplied_owner = _mapping(
            _mapping(parent.get("n1_actual_canonical_solution_binding")).get(
                "canonical_owner"
            )
        )
        supplied_node_bindings = _mapping(supplied_owner.get("node_bindings"))
        owner_map_match = all(
            supplied_node_bindings.get(key) == value
            for key, value in expected_owner.items()
        )

        core_ids = [index for values in expected_owner.values() for index in values]
        core_nodes = [nodes[int(index)] for index in core_ids] if isinstance(nodes, list) and len(nodes) > 43323 else []
        expected_quantities = (
            ["F1"] * 4
            + ["Uz1"] * 4
            + ["K1"] * 4
            + ["P1"] * 4
            + ["d_x_F1"] * 4
            + ["d_x_Uz1"] * 4
        )
        core_node_shape = bool(core_nodes) and all(
            node.get("chart") == "core"
            and node.get("operation") == "current_original_unique_regular_n1_component"
            and node.get("quantity") == expected_quantities[index]
            for index, node in enumerate(core_nodes)
        )
        green_provider_match = bool(core_nodes) and all(
            _mapping(node.get("defining_Green_provider")).get("function")
            == "green_apply"
            and _mapping(node.get("defining_Green_provider")).get("module_sha256")
            == RAW_VOLTERRA_SHA256
            and _mapping(node.get("defining_feedback_provider")).get("function")
            == "feedback_vectors"
            and _mapping(node.get("defining_feedback_provider")).get("module_sha256")
            == RAW_VOLTERRA_SHA256
            and node.get("exact_function_definition")
            == "W=sum_k>=0 H^k W0; W0=(I+J)Gg; H=(I+J)G(C0+B1*dZ)"
            for node in core_nodes
        )

        raw_operator_bindings = [
            _mapping(node.get("actual_raw_operator_bindings")).get("core")
            for node in core_nodes
        ]
        raw_core = raw_operator_bindings[0] if raw_operator_bindings else {}
        raw_operator_equal = bool(raw_operator_bindings) and all(
            binding == raw_core for binding in raw_operator_bindings
        )
        raw_b0 = _mapping(raw_core.get("B0"))
        raw_b1 = _mapping(raw_core.get("B1"))
        raw_g = _mapping(raw_core.get("g"))
        raw_operator_pressure_shape = (
            raw_core.get("D") == list(CANONICAL_D)
            and "4,1" in raw_b0
            and "4,1" not in raw_b1
            and "4" in raw_g
            and isinstance(raw_b0.get("4,1"), list)
            and len(raw_b0.get("4,1")) == 4
            and isinstance(raw_g.get("4"), list)
            and len(raw_g.get("4")) == 4
        )

        integrand_ids = core_integrands.get("m15")
        integrand_nodes = (
            [nodes[int(index)] for index in integrand_ids]
            if isinstance(nodes, list)
            and isinstance(integrand_ids, list)
            and all(isinstance(index, int) and index < len(nodes) for index in integrand_ids)
            else []
        )
        cumulative_ids = core_cumulative.get("m15")
        cumulative_nodes = (
            [nodes[int(index)] for index in cumulative_ids]
            if isinstance(nodes, list)
            and isinstance(cumulative_ids, list)
            and all(isinstance(index, int) and index < len(nodes) for index in cumulative_ids)
            else []
        )
        core_m15_substitution_nodes = (
            [nodes[index] for index in range(52597, 52601)]
            if isinstance(nodes, list) and len(nodes) > 52600
            else []
        )
        core_m15_substitution_exact = bool(core_m15_substitution_nodes) and all(
            node
            == {
                "operation": "function_substitution",
                "expression": 52576 + offset,
                "variable": 35558,
                "value": 52580,
                "Z_independent_substitution": True,
            }
            for offset, node in enumerate(core_m15_substitution_nodes)
        )

        substitution_53125 = (
            nodes[53125]
            if isinstance(nodes, list) and len(nodes) > 53125
            else {}
        )
        substitution_53126 = (
            nodes[53126]
            if isinstance(nodes, list) and len(nodes) > 53126
            else {}
        )
        node13 = nodes[13] if isinstance(nodes, list) and len(nodes) > 13 else {}
        node35558 = (
            nodes[35558] if isinstance(nodes, list) and len(nodes) > 35558 else {}
        )
        node52601 = (
            nodes[52601] if isinstance(nodes, list) and len(nodes) > 52601 else {}
        )
        node52602 = (
            nodes[52602] if isinstance(nodes, list) and len(nodes) > 52602 else {}
        )
        native_endpoint = (
            [nodes[53161], nodes[53162]]
            if isinstance(nodes, list) and len(nodes) > 53162
            else []
        )
        substitution_exact = (
            substitution_53125
            == {
                "Z_independent_substitution": True,
                "expression": 52601,
                "operation": "function_substitution",
                "value": 13,
                "variable": 35558,
            }
            and substitution_53126
            == {
                "Z_independent_substitution": True,
                "expression": 52602,
                "operation": "function_substitution",
                "value": 13,
                "variable": 35558,
            }
        )
        core_integral_nodes_exact = bool(cumulative_nodes) and all(
            node.get("operation") == "definite_integral"
            and node.get("paper_moment") == "m15"
            and node.get("variable") == "original_n1_cutoff_core_moment_t"
            and node.get("upper") == 35558
            and node.get("lower") == 0
            and node.get("exact_function_integral") is True
            and node.get("actual_signed_n1_inner_moment") is True
            and node.get("axis_integrand_regular") is True
            for node in cumulative_nodes
        )
        direct_density_nodes = bool(integrand_nodes) and all(
            node.get("operation") == "sum"
            and isinstance(node.get("arguments"), list)
            and len(node.get("arguments")) == 2
            for node in integrand_nodes
        )
        # Each accepted core m15 integrand is bound to two exact graph terms:
        # an axial product and a pressure product carrying an explicit minus
        # node.  Checking only the enclosing sum would lose the signed
        # pressure endpoint contribution.
        core_m15_term_specs = (
            (52576, 52563, 52550, 52100, 52572, 52568, 52124, 52567),
            (52577, 52564, 52553, 52100, 52573, 52569, 52125, 52567),
            (52578, 52565, 52557, 52100, 52574, 52570, 52126, 52567),
            (52579, 52566, 52562, 52100, 52575, 52571, 52127, 52567),
        )
        core_m15_term_binding = []
        for (
            integrand_id,
            axial_id,
            axial_left,
            axial_right,
            pressure_negative_id,
            pressure_product_id,
            pressure_left,
            pressure_right,
        ) in core_m15_term_specs:
            integrand = nodes[integrand_id] if isinstance(nodes, list) else {}
            axial = nodes[axial_id] if isinstance(nodes, list) else {}
            pressure_negative = (
                nodes[pressure_negative_id] if isinstance(nodes, list) else {}
            )
            pressure_product = (
                nodes[pressure_product_id] if isinstance(nodes, list) else {}
            )
            exact = (
                integrand
                == {
                    "operation": "sum",
                    "arguments": [axial_id, pressure_negative_id],
                }
                and axial
                == {
                    "operation": "product",
                    "arguments": [axial_left, axial_right],
                }
                and pressure_negative
                == {"operation": "negative", "argument": pressure_product_id}
                and pressure_product
                == {
                    "operation": "product",
                    "arguments": [pressure_left, pressure_right],
                }
            )
            core_m15_term_binding.append(
                {
                    "integrand_node": integrand_id,
                    "axial_product": {
                        "node": axial_id,
                        "arguments": [axial_left, axial_right],
                    },
                    "negative_pressure_product": {
                        "negative_node": pressure_negative_id,
                        "product_node": pressure_product_id,
                        "arguments": [pressure_left, pressure_right],
                    },
                    "exact": exact,
                }
            )
        core_m15_terms_exact = bool(core_m15_term_binding) and all(
            term.get("exact") is True for term in core_m15_term_binding
        )
        seam_history_nodes = (
            [nodes[index] for index in range(43369, 43373)]
            if isinstance(nodes, list) and len(nodes) > 43372
            else []
        )
        complete_history_not_reset_at_core_seam = bool(seam_history_nodes) and all(
            node.get("operation") == "definite_integral"
            and node.get("variable") == "original_n1_Green_t"
            and node.get("lower") == 0
            and node.get("upper") == 1
            and node.get("exact_function_integral") is True
            and node.get("actual_raw_regular_n1_Green_body") is True
            and node.get("complete_history_not_reset_at_core_seam") is True
            and node.get("signed_body_not_magnitude_cap") is True
            and node.get("endpoint_independent_of_Z") is True
            and node.get("zero_endpoint_integrand_continuation") == 0
            for node in seam_history_nodes
        )
        native_endpoint_exact = (
            endpoint.get("m15") == [53161, 53162, 53163, 53164]
            and len(native_endpoint) == 2
            and native_endpoint[0]
            == {
                "Z_independent_substitution": True,
                "expression": 53141,
                "operation": "function_substitution",
                "value": 35824,
                "variable": 35848,
            }
            and native_endpoint[1]
            == {
                "Z_independent_substitution": True,
                "expression": 53142,
                "operation": "function_substitution",
                "value": 35824,
                "variable": 35848,
            }
        )
        checks = {
            "archive_hash": _sha256(path) == INNER_CUTOFF_ARCHIVE_SHA256,
            "parent_input_hash": parent_inputs.get(Path(INNER_CUTOFF_ARCHIVE).name)
            == INNER_CUTOFF_ARCHIVE_SHA256,
            "source_family_match": archive.get("source_family") == parent.get("source_family"),
            "exact_graph_node_count": isinstance(nodes, list) and len(nodes) == 53165,
            "core_cumulative_m15_nodes": core_cumulative.get("m15")
            == [52601, 52602, 52603, 52604],
            "core_m15_integrand_nodes": integrand_ids == [52576, 52577, 52578, 52579],
            "core_integrand_shapes": direct_density_nodes,
            "core_m15_exact_product_sign_terms": core_m15_terms_exact,
            "core_cumulative_integral_shapes": core_integral_nodes_exact,
            "core_m15_integral_substitutions": core_m15_substitution_exact,
            "complete_history_not_reset_at_core_seam": complete_history_not_reset_at_core_seam,
            "core_endpoint_substitutions": substitution_exact,
            "core_endpoint_substitution_axis": node13
            == {"denominator": 1, "numerator": 4, "operation": "exact_rational"}
            and node35558
            == {"name": "original_n1_true_core_rho", "operation": "bound_variable"},
            "extended_prefix_owner_map": owner_map_match,
            "extended_prefix_core_nodes": core_node_shape,
            "extended_prefix_volterra_graph": green_provider_match,
            "extended_prefix_raw_operator_same_source": raw_operator_equal
            and raw_operator_pressure_shape,
            "native_C_inner_endpoint_nodes": native_endpoint_exact,
            "definition_core_pressure_source_retained": definition.get(
                "leading_pressure_forcing_unchanged"
            )
            is True
            and definition.get(
                "pressure_and_two_pressure_moments_keep_leading_source_after_cutoff_zero"
            )
            is True,
            "definition_m15_endpoint_sign": definition.get(
                "m15_includes_negative_R_pressure_derivative"
            )
            is True,
        }
        return {
            "passed": all(checks.values()),
            "checks": checks,
            "path": INNER_CUTOFF_ARCHIVE,
            "sha256": _sha256(path),
            "expected_sha256": INNER_CUTOFF_ARCHIVE_SHA256,
            "graph_node_count": len(nodes) if isinstance(nodes, list) else None,
            "core_cumulative_m15": core_cumulative.get("m15"),
            "core_m15_integrands": integrand_ids,
            "core_endpoint_substitution_nodes": [53125, 53126],
            "core_endpoint_substitution_expression_nodes": [52601, 52602],
            "core_endpoint_rho_node": 13,
            "core_endpoint_rho": "4",
            "core_endpoint_variable": 35558,
            "core_integrand_density_nodes": integrand_ids,
            "core_m15_integral_substitution_nodes": [52597, 52598, 52599, 52600],
            "core_m15_exact_product_sign_terms": core_m15_term_binding,
            "core_seam_history_nodes": [43369, 43370, 43371, 43372],
            "complete_history_not_reset_at_core_seam": complete_history_not_reset_at_core_seam,
            "native_C_inner_endpoint_nodes": [53161, 53162],
            "native_C_inner_includes_first_collar": True,
            "native_C_inner_endpoint_moment_nodes": endpoint.get("m15"),
            "extended_prefix_core_nodes": list(range(43300, 43324)),
            "extended_prefix_quantities": expected_quantities,
            "extended_prefix_raw_operator_source": RAW_VOLTERRA_SOURCE,
            "density_link": {
                "integrand_nodes": integrand_ids,
                "identity": "inner-cutoff m15 rho integrand=2*epsilon*u*U-epsilon^2*r*(2*f*F+N1p), equivalently dR density=2*U_lead*U1-R*(2*F_lead*F1+N1p); cumulative split 2*U_lead*U1+P1 is paired with endpoint -R*P1",
                "source_executed": False,
            },
            "source_executed": False,
            "archive_source_executed": False,
        }
    except Exception as exc:
        return {
            "passed": False,
            "checks": {"archive_binding": False},
            "failure": f"inner-cutoff graph binding failed closed: {type(exc).__name__}: {exc}",
            "path": INNER_CUTOFF_ARCHIVE,
            "expected_sha256": INNER_CUTOFF_ARCHIVE_SHA256,
            "source_executed": False,
            "archive_source_executed": False,
        }


def _endpoint_handle_binding(
    repo: Path,
    parent: Mapping[str, Any],
    inner_cutoff: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Bind requested core endpoint handles only when source data has them.

    The accepted report and source definitions are treated as data and AST
    text.  A decimal substring inside an exact coefficient is not a handle;
    JSON scalar equality or an AST integer constant is required.
    """

    requested = (53125, 53126)
    locations: dict[str, list[str]] = {str(value): [] for value in requested}

    def visit(value: Any, path: str) -> None:
        if isinstance(value, Mapping):
            for key, item in value.items():
                visit(item, f"{path}.{key}")
        elif isinstance(value, list):
            for index, item in enumerate(value):
                visit(item, f"{path}[{index}]")
        elif value in requested:
            locations[str(value)].append(path)

    visit(parent, PARENT_REPORT)
    ast_sources = (RAW_VOLTERRA_SOURCE, RAW_OPERATOR_SOURCE, CUMULATIVE_SOURCE)
    for relative in ast_sources:
        path = _path(repo, relative)
        if not path.exists():
            continue
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
        except (OSError, SyntaxError):
            continue
        for node in ast.walk(tree):
            if isinstance(node, ast.Constant) and node.value in requested:
                locations[str(node.value)].append(f"{relative}:line{node.lineno}")

    json_sources = (CUMULATIVE_RECEIPT,)
    for relative in json_sources:
        path = _path(repo, relative)
        if not path.exists():
            continue
        try:
            visit(_read_json(path), relative)
        except (OSError, ValueError, json.JSONDecodeError):
            continue

    # The exact endpoint substitutions live in the accepted inner-cutoff
    # graph, rather than in the raw source modules or cumulative receipt.
    if _mapping(inner_cutoff).get("checks", {}).get(
        "core_endpoint_substitutions"
    ) is True:
        locations["53125"].append(f"{INNER_CUTOFF_ARCHIVE}.exact_graph_nodes[53125]")
        locations["53126"].append(f"{INNER_CUTOFF_ARCHIVE}.exact_graph_nodes[53126]")

    available = all(locations[str(value)] for value in requested)
    return {
        "requested": list(requested),
        "source_exact_available": available,
        "bound": available,
        "locations": locations,
        "reason": (
            "exact source/report endpoint handles are bound"
            if available
            else "no exact endpoint handles 53125 or 53126 occur in the accepted source/report graph"
        ),
        "optional": True,
        "source_executed": False,
    }


def _failure(reason: str, parent_pin: Mapping[str, Any] | None = None) -> dict[str, Any]:
    result: dict[str, Any] = {
        "passed": False,
        "parent_report": PARENT_REPORT,
        "parent_report_sha256": PARENT_REPORT_SHA256,
        "failure": reason,
        "checks": {"passed": False, "failure": reason},
        "actual_core_pressure_m15_functional_identity_bound": False,
        "actual_core_m15_pressure_primitive_bound": False,
        "actual_raw_pressure_rho_equation_bound": False,
        "actual_pressure_rhs_factor_bound": False,
        "actual_core_differentiated_density_bound": False,
        "actual_canonical_core_endpoint_handle_bound": False,
        "accepted_inner_cutoff_graph": False,
        "actual_canonical_core_cumulative_binding": False,
        "actual_symbolic_primitive_identity": False,
        "actual_zero_axis_pressure_gauge_bound": False,
        "source_executed": False,
        "old_pure_definitions_executed": False,
        "archive_source_executed": False,
        "accepted_owner_builder_fixture_checker_recurrence_norm_quadrature_replayed": False,
        "point_values_selected": False,
        "caps_as_values": False,
        "full_C_inner_identity": False,
        "q_dependent_collar_bridge": False,
        "P1_m13_equals_m15": False,
        "six_term_source_equality": False,
        "global_moment_closure": False,
        "remaining_scope": True,
    }
    if parent_pin is not None:
        result["parent_report_pin"] = dict(parent_pin)
    return result


def build(root: str | Path | None = None, parent: Mapping[str, Any] | None = None) -> dict[str, Any]:
    """Build the source-bound core m15 pressure primitive receipt."""

    parent_pin: dict[str, Any] | None = None
    try:
        repo = _repo_root(root)
        parent_path = _path(repo, PARENT_REPORT)
        actual_parent_sha = _sha256(parent_path)
        if actual_parent_sha != PARENT_REPORT_SHA256:
            return _failure("accepted forward-enclosure report hash changed")
        saved_parent = _read_json(parent_path)
        supplied = saved_parent if parent is None else parent
        supplied_shape = "complete_forward_enclosure_report"
        if "n1_infinite_forward_analytic_response" not in supplied:
            nested = supplied
            accepted_nested = saved_parent.get("n1_infinite_forward_analytic_response")
            if _digest(nested) != _digest(accepted_nested):
                return _failure("supplied parent is not the accepted forward packet")
            supplied = saved_parent
            supplied_shape = "nested_forward_packet"
        elif _digest(supplied) != _digest(saved_parent):
            return _failure("supplied parent mapping does not match accepted report")
        parent_pin = {
            "path": PARENT_REPORT,
            "sha256": actual_parent_sha,
            "expected_sha256": PARENT_REPORT_SHA256,
            "mapping_digest": _digest(supplied),
            "supplied_mapping_shape": supplied_shape,
            "source_executed": False,
        }

        pins = _source_pins(repo, supplied)
        ast_bindings = _source_ast_bindings(repo)
        forward_binding = _forward_binding(supplied)
        cumulative_binding = _cumulative_receipt_binding(repo, supplied)
        symbolic = _symbolic_identity()
        inner_cutoff = _inner_cutoff_binding(repo, supplied)
        endpoint_handle = _endpoint_handle_binding(repo, supplied, inner_cutoff)
        density_link_checks = {
            "source_core_integrands": _mapping(inner_cutoff.get("checks")).get(
                "core_integrand_shapes"
            )
            is True,
            "source_core_cumulative_m15": _mapping(inner_cutoff.get("checks")).get(
                "core_cumulative_integral_shapes"
            )
            is True,
            "source_core_m15_integral_substitutions": _mapping(
                inner_cutoff.get("checks")
            ).get("core_m15_integral_substitutions")
            is True,
            "source_signed_integrand_ast": _mapping(
                ast_bindings["bindings"].get(
                    "inner_cutoff_signed_rho_integrands"
                )
            ).get("passed")
            is True,
            "pressure_rhs_factor": symbolic.get("checks", {}).get(
                "pressure_row_factor"
            )
            is True,
            "FTC_pressure_endpoint": symbolic.get("checks", {}).get(
                "FTC_pressure_density"
            )
            is True,
            "physical_core_density": symbolic.get("checks", {}).get(
                "physical_core_density"
            )
            is True,
        }
        inner_cutoff["density_link"] = {
            **_mapping(inner_cutoff.get("density_link")),
            "pressure_rhs": "P1_rho=epsilon*(2*F_lead*F1+N1p)",
            "direct_inner_cutoff_density": "2*epsilon*u*U-epsilon^2*r*(2*f*F+N1p)",
            "FTC": "d/dR[integral(2*U_lead*U1+P1)dR-R*P1]=2*U_lead*U1-R*P1_rho",
            "physical_density": "2*U_lead*U1-R*(2*F_lead*F1+N1p)",
            "checks": density_link_checks,
            "passed": all(density_link_checks.values()),
            "source_executed": False,
        }
        inner_cutoff["checks"] = {
            **_mapping(inner_cutoff.get("checks")),
            "density_pressure_FTC_link": inner_cutoff["density_link"]["passed"],
        }
        inner_cutoff["passed"] = inner_cutoff.get("passed") is True and inner_cutoff[
            "density_link"
        ]["passed"] is True

        ast_checks = {
            "all_bindings": ast_bindings.get("passed") is True,
            "cumulative_rho_integrands": _mapping(
                ast_bindings["bindings"].get("cumulative_rho_integrands")
            ).get("passed")
            is True,
            "cumulative_endpoint_cumulants": _mapping(
                ast_bindings["bindings"].get("cumulative_endpoint_cumulants")
            ).get("passed")
            is True,
            "cumulative_core_range": _mapping(
                ast_bindings["bindings"].get("cumulative_core_range_cell")
            ).get("passed")
            is True,
            "cumulative_core_evaluate": _mapping(
                ast_bindings["bindings"].get("cumulative_core_evaluate")
            ).get("passed")
            is True,
            "m15_signed_rows": _mapping(
                ast_bindings["bindings"].get("m15_signed_prefix_rows")
            ).get("passed")
            is True,
            "m15_evaluate": _mapping(ast_bindings["bindings"].get("m15_evaluate")).get(
                "passed"
            )
            is True,
            "raw_system": _mapping(ast_bindings["bindings"].get("raw_system")).get(
                "passed"
            )
            is True,
            "raw_green_apply": _mapping(
                ast_bindings["bindings"].get("raw_green_apply")
            ).get("passed")
            is True,
            "inner_cutoff_signed_rho_integrands": _mapping(
                ast_bindings["bindings"].get(
                    "inner_cutoff_signed_rho_integrands"
                )
            ).get("passed")
            is True,
            "raw_operator": _mapping(
                ast_bindings["bindings"].get("raw_operator_matrix")
            ).get("passed")
            is True,
            "raw_matrix_no_pressure_self": _mapping(
                ast_bindings["bindings"].get("raw_matrix_no_pressure_self")
            ).get("passed")
            is True,
            "raw_pressure_source_binding": _mapping(
                ast_bindings["bindings"].get("raw_pressure_source_binding")
            ).get("passed")
            is True,
            "raw_core_fields": _mapping(
                ast_bindings["bindings"].get("raw_core_fields")
            ).get("passed")
            is True,
            "finite_forward_ast": _mapping(
                ast_bindings["bindings"].get("finite_forward_advance")
            ).get("passed")
            is True,
        }
        pressure_row = {
            "equation": RAW_PRESSURE_RHO_EQUATION,
            "canonical_equation": "W_x + D*W/x = B0*W + B1*D_Z W + g",
            "D": list(CANONICAL_D),
            "pressure_row_index": 4,
            "D4": CANONICAL_D[3],
            "B0_row4": "B0['4,1']=4*epsilon*x*F_lead; B0['4,j']=0 for j!=1",
            "B1_row4": "B1['4,j']=0 for every j",
            "g4": "2*epsilon*x*N1p",
            "P1_x": "4*epsilon*x*F_lead*F1+2*epsilon*x*N1p",
            "P1_rho": "epsilon*(2*F_lead*F1+N1p)",
            "p_rho": "2*F0^2*Phi*f/(epsilon*Pstar^2)+N1p/Pstar^2",
            "raw_row_equation": RAW_PRESSURE_ROW_EQUATION,
            "raw_matrix_row4_literal_bound": _mapping(
                ast_bindings["bindings"].get("raw_pressure_source_binding")
            ).get("row4_literal_ast_bound")
            is True,
            "raw_matrix_field_meaning_bound": _mapping(
                ast_bindings["bindings"].get("raw_pressure_source_binding")
            ).get("core_fields_ast_bound")
            is True,
            "raw_matrix_B1_row4_absent": _mapping(
                ast_bindings["bindings"].get("raw_pressure_source_binding")
            ).get("B1_row4_absent")
            is True,
            "symbolic_pressure_row_factor": symbolic.get("checks", {}).get(
                "pressure_row_factor"
            )
            is True,
            "symbolic_normalized_pressure_row_factor": symbolic.get(
                "checks", {}
            ).get("normalized_pressure_row_factor")
            is True,
            "symbolic_physical_core_density": symbolic.get("checks", {}).get(
                "physical_core_density"
            )
            is True,
            "valid_for": "0<rho<=4; ordinary Z0 and Z1",
            "axis": "P1(0,Z)=p(0,Z)=0; use parity extension, never divide at rho=0",
            "pressure_row_has_no_self_W4": _mapping(
                ast_bindings["bindings"].get("raw_matrix_no_pressure_self")
            ).get("passed")
            is True,
        }
        pressure_equation_passed = all(
            (
                pressure_row.get("D4") == 0,
                pressure_row.get("pressure_row_has_no_self_W4") is True,
                pressure_row.get("raw_matrix_row4_literal_bound") is True,
                pressure_row.get("raw_matrix_field_meaning_bound") is True,
                pressure_row.get("raw_matrix_B1_row4_absent") is True,
                pressure_row.get("symbolic_pressure_row_factor") is True,
                pressure_row.get("symbolic_normalized_pressure_row_factor") is True,
                pressure_row.get("symbolic_physical_core_density") is True,
            )
        )
        scale_binding = {
            "normalization": dict(NORMALIZATION),
            "physical_fields": {
                "U1": "epsilon*u",
                "P1": "epsilon*Pstar^2*p",
                "R": "epsilon*rho",
                "dR": "epsilon*drho",
            },
            "axial_factor": "epsilon^2",
            "pressure_factor": "epsilon^2*Pstar^2",
            "pressure_endpoint": "-R*P1=-epsilon^2*Pstar^2*rho*p(rho)",
            "powers_applied_once": ast_checks["m15_evaluate"]
            and ast_checks["m15_signed_rows"],
            "gauge_constant_cancels": symbolic.get("checks", {}).get(
                "pressure_gauge_cancels"
            )
            is True,
        }
        canonical_port = {
            "chart": "core",
            "source_module": CUMULATIVE_SOURCE,
            "functions": ["rho_integrands", "endpoint_cumulants"],
            "rho_integrands": {
                "m15_parts": "dR source maps to 2*U_lead*U1+P1; after U1=epsilon*u and dR=epsilon*drho, the two scales are recorded below",
                "mapped_axial": "epsilon^2*2*U_lead*u",
                "mapped_pressure_integrand": "epsilon^2*Pstar^2*p",
                "inner_cutoff_direct_density": "2*epsilon*u*U-epsilon^2*r*(2*f*F+N1p)",
                "inner_cutoff_integral_substitutions": {
                    "52597": "expression=52576, variable=35558, value=52580",
                    "52598": "expression=52577, variable=35558, value=52580",
                    "52599": "expression=52578, variable=35558, value=52580",
                    "52600": "expression=52579, variable=35558, value=52580",
                },
            },
            "endpoint_cumulants": {
                "m15_boundary": "-rho*P1 in rho source coordinates = -epsilon^2*Pstar^2*rho*p",
                "m13": "W[3]=P1 endpoint only; this is not m15",
            },
            "pressure_rhs": {
                "raw_row": RAW_PRESSURE_ROW_EQUATION,
                "field_meaning": "raw f=F_lead=F0*Phi; canonical F1=F0*f/epsilon",
                "rho_equation": "P1_rho=epsilon*(2*F_lead*F1+N1p)",
                "normalized_equation": "p_rho=2*F0^2*Phi*f/(epsilon*Pstar^2)+N1p/Pstar^2",
                "physical_m15_density": "d/dR m15=2*U_lead*U1-R*(2*F_lead*F1+N1p)",
            },
            "core_evaluate_combination": "integrals['m15_parts']+endpoint['m15_boundary']",
            "inner_cutoff_graph": {
                "path": INNER_CUTOFF_ARCHIVE,
                "core_integrand_nodes": [52576, 52577, 52578, 52579],
                "core_cumulative_m15_nodes": [52601, 52602, 52603, 52604],
                "core_endpoint_nodes": [53125, 53126],
                "native_C_inner_endpoint_nodes": [53161, 53162],
                "native_C_inner_includes_first_collar": True,
                "extended_prefix_core_nodes": list(range(43300, 43324)),
                "binding": inner_cutoff,
            },
            "canonical_actual_core_endpoint_handle": {
                **endpoint_handle,
            },
            "q_collar_included": False,
            "full_C_inner_includes_qcollar": True,
        }
        domain = {
            "rho": "[0,4]",
            "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
            "source_identity_first_domain": "rho>0 and L!=0",
            "axis_extension": "regular parity and zero-axis identity extend to rho=0",
            "no_global_arbitrary_Z_claim": True,
            "no_selected_interval_center": True,
        }
        scope = {
            "core_m15_identity": True,
            "full_C_inner_identity": False,
            "q_dependent_collar_bridge": False,
            "P1_m13_equals_m15": False,
            "six_term_source_equality": False,
            "global_moment_closure": False,
            "temporal_n_ge_2": False,
            "full_forward_response_function_identity": False,
            "remaining_analytic_remainder": True,
            "finite_prefix_radial_remainder_bound": False,
        }
        checks = {
            "parent_report": parent_pin.get("sha256") == PARENT_REPORT_SHA256,
            "source_pins": pins.get("all_hashes_passed") is True,
            "ast_source_bindings": all(ast_checks.values()),
            "accepted_forward_object": forward_binding.get("passed") is True,
            "accepted_cumulative_receipt": cumulative_binding.get("passed") is True,
            "accepted_inner_cutoff_graph": inner_cutoff.get("passed") is True,
            "canonical_core_endpoint_handle": endpoint_handle.get("bound") is True,
            "core_density_pressure_FTC_link": inner_cutoff.get("density_link", {}).get(
                "passed"
            )
            is True,
            "fresh_symbolic_identity": symbolic.get("passed") is True,
            "pressure_equation": pressure_equation_passed,
            "exact_scales": scale_binding.get("powers_applied_once") is True,
            "zero_axis": forward_binding.get("checks", {}).get("canonical_binding") is True
            and _mapping(
                _mapping(supplied.get("n1_actual_canonical_solution_binding")).get(
                    "axis_weighted_transform"
                )
            ).get("zero_axis")
            == "W(0,Z)=0; W1..W4=O(x^2), W5,W6=O(x)",
            "ordinary_Z0_Z1": domain.get("ordinary_Z_orders") == [0, 1],
            "core_scope_only": scope["core_m15_identity"]
            and not scope["full_C_inner_identity"],
            "no_replay": True,
        }
        passed = all(checks.values())
        return {
            "passed": passed,
            "parent_report": PARENT_REPORT,
            "parent_report_sha256": PARENT_REPORT_SHA256,
            "parent_report_pin": parent_pin,
            "checks": {**checks, "passed": passed},
            "actual_core_pressure_m15_functional_identity_bound": passed,
            "actual_core_m15_pressure_primitive_bound": passed,
            "actual_raw_pressure_rho_equation_bound": pressure_equation_passed,
            "actual_pressure_rhs_factor_bound": pressure_equation_passed,
            "actual_core_differentiated_density_bound": symbolic.get("checks", {}).get(
                "physical_core_density"
            )
            is True,
            "actual_canonical_core_cumulative_binding": cumulative_binding.get("passed") is True,
            "actual_symbolic_primitive_identity": symbolic.get("passed") is True,
            "actual_zero_axis_pressure_gauge_bound": pressure_row.get("D4") == 0
            and pressure_row.get("axis", "").startswith("P1(0,Z)=p(0,Z)=0"),
            "source_pins": pins,
            "source_AST_bindings": ast_bindings,
            "accepted_forward_binding": forward_binding,
            "accepted_cumulative_receipt": cumulative_binding,
            "accepted_inner_cutoff_binding": inner_cutoff,
            "actual_raw_pressure_rho_equation": pressure_row,
            "canonical_actual_core_endpoint_handle": endpoint_handle,
            "actual_canonical_core_endpoint_handle_bound": endpoint_handle.get("bound")
            is True,
            "normalization_and_scale_binding": scale_binding,
            "canonical_core_cumulative_port": canonical_port,
            "symbolic_identity": symbolic,
            "domain": domain,
            "forward_first_two_terms": {
                "formula": ACTUAL_M15,
                "axial": "epsilon^2*2*integral_0^rho U_lead*u drho",
                "pressure": "epsilon^2*Pstar^2*(integral_0^rho p drho-rho*p(rho))",
                "pressure_rho_derivative": "p_rho=2*F0^2*Phi*f/(epsilon*Pstar^2)+N1p/Pstar^2",
                "physical_core_density": "d/dR m15=2*U_lead*U1-R*(2*F_lead*F1+N1p)",
                "ordinary_Z_orders": list(ORDINARY_Z_ORDERS),
                "source_enclosure_errors_retained": True,
                "source_values_selected": False,
            },
            "scope": scope,
            "source_executed": False,
            "old_pure_definitions_executed": False,
            "archive_source_executed": False,
            "accepted_owner_builder_fixture_checker_recurrence_norm_quadrature_replayed": False,
            "point_values_selected": False,
            "caps_as_values": False,
            "tool_metadata": {
                "route": "fixed luna_max_worker",
                "initial_fork_turns": "none",
                "requested_model": "GPT-5.6 Luna/max",
                "independent_runtime_model_metadata_exposed": False,
            },
        }
    except Exception as exc:
        return _failure(
            f"source binding failed closed: {type(exc).__name__}: {exc}", parent_pin
        )


__all__ = [
    "PARENT_REPORT",
    "PARENT_REPORT_SHA256",
    "ACTUAL_M15",
    "RAW_PRESSURE_RHO_EQUATION",
    "RAW_PRESSURE_ROW_EQUATION",
    "INNER_CUTOFF_ARCHIVE",
    "build",
]
