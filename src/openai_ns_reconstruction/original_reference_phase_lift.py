"""Source-bound smooth periodic lift certificate for the prescribed loop.

The finite phase graph stores the fundamental inverse on ``[0, 2*pi]`` and
then stores ``fractional_part(N*global_offset)``.  This module supplies the
analytic extension needed at an integer cut: the Poisson loop and its
primitive obey exact period shifts with all finite mixed parameter jets. It
only reads an accepted graph and the pinned source file.  It does not build
a producer, evaluate caps, or identify a physical ``y`` derivative.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import Any, Mapping

from .original_reference_germ_guards import _function_evidence


_LOOP_SOURCE = (
    "experiments/root_st073/"
    "lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"
)
_REQUIRED_SOURCE_CHECKS = (
    "source_pins",
    "accepted_archives",
    "parameters",
    "phase_global_R_Z0",
    "phase_source_recipe",
    "generic_roots",
    "flat_cutoff",
)


def _mapping(value: Any) -> Mapping[str, Any]:
    return value if isinstance(value, Mapping) else {}


def _node(nodes: Any, index: Any) -> Mapping[str, Any] | None:
    if isinstance(index, bool):
        return None
    if not isinstance(index, int):
        return None
    try:
        row = nodes[index] if not isinstance(nodes, Mapping) else nodes.get(index, nodes.get(str(index)))
    except (IndexError, KeyError, TypeError):
        return None
    return row if isinstance(row, Mapping) else None


def _repo_root(root: Any) -> Path:
    if root is not None:
        path = Path(root)
        if path.is_file():
            path = path.parent
        if (path / "experiments").is_dir() or (path / "src").is_dir():
            return path
    return Path(__file__).resolve().parents[2]


def _source_pin(path: Path, pin: Any, root: Path) -> dict[str, Any]:
    """Validate the caller's exact inherited pin; never accept a fresh hash."""
    row = _mapping(pin)
    canonical = Path(_LOOP_SOURCE).as_posix()
    supplied_path = str(row.get("path", "")).replace("\\", "/")
    path_ok = supplied_path == canonical
    exists = path.is_file()
    actual = hashlib.sha256(path.read_bytes()).hexdigest() if exists else None
    expected = row.get("expected_sha256")
    pin_ok = bool(
        path_ok
        and exists
        and row.get("exists") is True
        and row.get("pinned") is True
        and isinstance(expected, str)
        and expected == actual
        and row.get("sha256") == actual
    )
    return {
        "passed": pin_ok,
        "path": supplied_path,
        "canonical_path": canonical,
        "exists": exists,
        "sha256": actual,
        "expected_sha256": expected,
        "inherited_pin": row,
    }


def _phase_rows(nodes: Any, loop: Mapping[str, Any]) -> tuple[bool, dict[str, Any]]:
    charts = _mapping(loop.get("charts"))
    common_n = loop.get("N")
    required_names = ("reference", "restoration")
    required_present = all(name in charts and isinstance(charts[name], Mapping) for name in required_names)
    rows: dict[str, Any] = {}
    passed = required_present and isinstance(common_n, int) and not isinstance(common_n, bool)
    for name, chart_value in charts.items():
        chart = _mapping(chart_value)
        if not chart or "actual_phase" not in chart:
            continue
        phase_handle = chart.get("actual_phase")
        phase = _node(nodes, phase_handle)
        argument = phase.get("argument") if phase else None
        product = _node(nodes, argument)
        args = product.get("arguments", ()) if product else ()
        if not isinstance(args, (list, tuple)):
            args = ()
        n_count = sum(value == common_n for value in args)
        offset = chart.get("global_offset")
        offset_count = sum(value == offset for value in args)
        row_passed = bool(
            phase
            and phase.get("operation") == "analytic_unary"
            and phase.get("name") == "fractional_part"
            and product
            and product.get("operation") == "product"
            and n_count == 1
            and offset_count == 1
            and len(args) == 2
            and phase_handle not in {chart.get("phase_center"), chart.get("phase_origin")}
        )
        rows[name] = {
            "passed": row_passed,
            "coordinate": chart.get("coordinate"),
            "global_offset": offset,
            "actual_phase": phase_handle,
            "phase_node": phase,
            "product_node": product,
            "domain": chart.get("domain"),
        }
        passed = passed and row_passed
    passed = passed and all(name in rows for name in required_names)
    return passed, rows


def _formal_identities() -> dict[str, Any]:
    """Exact algebra plus explicit analytic steps for the smooth global lift.

    q is real, including zero. The positive a,d domain is the active
    universal recipe; the source lazy branch is treated separately below.
    Derivatives are taken on actual shift differences, never on a supplied
    zero or on a numerical endpoint value.
    """
    import sympy as sp

    a, E, d, eta = sp.symbols("a E d eta", positive=True)
    b, q, p2 = sp.symbols("b q p2", real=True)
    u = p2*q/d
    r = u/sp.sqrt(1+u**2)
    alpha = 2*q/sp.sqrt(1+u**2)
    t0 = -b/a
    gap = sp.factor(1-r**2)
    K = 1/(2*sp.pi*(1+t0**2+2*q**2))
    w, rho = sp.symbols("w rho", nonzero=True)
    cosine = (w+1/w)/2
    D = 1-rho*(w+1/w)+rho**2
    series = (w/(1-rho*w)+(1/w)/(1-rho/w))/2
    fourier_difference = sp.factor((cosine-rho)/D-series)
    # For |rho|<1, this Laurent identity is the sum of two convergent
    # geometric series on |w|=1. Integrating their cosine coefficients gives
    # mean zero and mean square 1/(2*(1-rho^2)), by orthogonality.
    period = sp.simplify(K*(2*sp.pi+2*sp.pi*(t0**2+2*q**2)))
    phi, psi, T1 = sp.symbols("phi psi T1", real=True)
    A = a*(phi-psi/(2*sp.pi))/2
    B = E*(-a*T1/(2*sp.pi)-b*phi)/2
    shift = {phi:phi+1,psi:psi+2*sp.pi,T1:T1+2*sp.pi*t0}
    A_difference = A.xreplace(shift)-A
    B_difference = B.xreplace(shift)-B
    params=(a,b,E,q,p2,d)
    # Smooth quasiperiodic inverse and T1 identities hold on an open
    # parameter domain, so their substitution into these finite differences
    # is valid before taking parameter or phase derivatives.
    variables=params+(phi,psi,T1)
    derivative_rows=[]
    for name,expression in (("A",A_difference),("B",B_difference)):
        for i,x in enumerate(variables):
            derivative_rows.append(dict(function=name,variables=[str(x)],
                exact_zero=sp.simplify(sp.diff(expression,x))==0))
            for z in variables[i:]:
                derivative_rows.append(dict(function=name,variables=[str(x),str(z)],
                    exact_zero=sp.simplify(sp.diff(expression,x,z))==0))
    theta, lam = sp.symbols("theta lam", real=True, nonzero=True)
    Aphi = a*(1-1/(2*sp.pi*lam))/2
    Bphi = E*(-a*theta/(2*sp.pi*lam)-b)/2
    # theta=t(psi) and lam=lambda(psi) are 2*pi-periodic. These actual
    # phase-derivative formulas therefore have equal endpoint functions.
    Aphi_chain=sp.diff(A,phi)+sp.diff(A,psi)/lam
    Bphi_chain=sp.diff(B,phi)+sp.diff(B,T1)*theta/lam
    phase_shift_checks=[sp.simplify(Aphi_chain-Aphi)==0,sp.simplify(Bphi_chain-Bphi)==0]
    flat_psi=2*sp.pi*phi
    flat_T1=t0*flat_psi
    A_flat_raw=A.subs(psi,flat_psi)
    B_flat_raw=B.subs(T1,flat_T1)
    A_flat=sp.simplify(A_flat_raw)
    B_flat=sp.simplify(B_flat_raw)
    Aphi_flat=sp.simplify(Aphi.subs({theta:t0,lam:1/(2*sp.pi)}))
    Bphi_flat=sp.simplify(Bphi.subs({theta:t0,lam:1/(2*sp.pi)}))
    flat_tangent_checks=[]
    for expression in (A_flat_raw,B_flat_raw):
        for x in (a,b,E,p2,d,phi):
            flat_tangent_checks.append(sp.simplify(sp.diff(expression,x))==0)
            for z in (a,b,E,p2,d,phi):
                flat_tangent_checks.append(sp.simplify(sp.diff(expression,x,z))==0)
    Delta=sp.symbols("Delta",real=True)
    gamma_at_cut=sp.simplify((2*eta-Delta).subs(Delta,eta))
    checks={
        "one_minus_r2":sp.simplify(gap-1/(1+u**2))==0,
        "strict_r_inside_unit_interval":bool(gap>0),
        "q_zero_included":q.is_positive is None and sp.simplify(r.subs(q,0))==0,
        "alpha_squared_identity":sp.simplify(alpha**2-4*q**2*gap)==0,
        "poisson_geometric_series":fourier_difference==0,
        "phase_period":period==1,
        "K_positive":bool(K>0),
        "A_period_difference":sp.simplify(A_difference)==0,
        "B_period_difference":sp.simplify(B_difference)==0,
        "parameter_and_mixed_jet_differences":all(v['exact_zero'] for v in derivative_rows),
        "phase_derivative_endpoint_formulas":all(phase_shift_checks),
        "flat_A_is_zero_function":A_flat==0,
        "flat_B_is_zero_function":B_flat==0,
        "flat_Aphi_is_zero_function":Aphi_flat==0,
        "flat_Bphi_is_zero_function":Bphi_flat==0,
        "fixed_q_zero_tangent_mixed_jets":all(flat_tangent_checks),
        "positive_gamma_at_flat_boundary":gamma_at_cut==eta and bool(gamma_at_cut>0),
    }
    checks={key:bool(value) for key,value in checks.items()}
    return dict(passed=all(checks.values()),checks=checks,
        parameter_domain='universal loop: a,E,d_star,eta>0; b,p2,q real with q independent; t0=-b/a; u=p2*q/d_star',
        source_sqrt_domain='active cutoff requires a>0 and 2*eta-Delta>0; proved with strict gap at the physical seam',
        derivative_rows=derivative_rows,derivative_row_count=len(derivative_rows),
        exact_shift_differences={'A':str(sp.factor(A_difference)),'B':str(sp.factor(B_difference))},
        identities=dict(r='u/sqrt(1+u^2)',alpha='2*q/sqrt(1+u^2)',
            fourier='(cos(psi)-r)/D=sum_{k>=1}r^(k-1)cos(k*psi)',
            integral_t='2*pi*t0',integral_t_squared='2*pi*(t0^2+2*q^2)',
            phase_period='Phi(psi+2*pi)=Phi(psi)+1',
            inverse_period='psi(phi+1)=psi(phi)+2*pi',
            primitive_period='T1(psi+2*pi)=T1(psi)+2*pi*t0',
            loop_period='A(phi+1)=A(phi), B(phi+1)=B(phi)',
            phase_derivatives='A_phi=a*(1-1/(2*pi*lambda))/2; B_phi=E*(-a*t/(2*pi*lambda)-b)/2'),
        analytic_proof=[
            'For finite real u, |r|<1. D=|1-r*exp(i*psi)|^2>0 at every angle.',
            'The geometric cosine series converges absolutely, uniformly with all finite parameter/angle derivatives on compact subsets of |r|<1. Cosine orthogonality gives integral(t)=2*pi*t0 and integral(t^2)=2*pi*(t0^2+2*q^2).',
            'K>0 and lambda=K*(1+t^2)>0. Periodic t gives Phi(psi+2*pi)=Phi(psi)+1, so Phi is a smooth increasing bijection of R. Its unique inverse is smooth in phase and finite parameters by the implicit function theorem.',
            'The inverse and T1 shift identities inserted into the actual A/B formulas give the exact differences checked here. These are function identities on an open parameter domain, so every finite parameter/phase mixed derivative is periodic.',
            'Composition with fractional_part equals this smooth periodic extension at integer phase cuts. In particular it admits ordinary Z1 of a physical-y derivative (second mixed jets); no derivative of a discontinuous bare fractional_part is used.',
            'At q=0, r=alpha=0, Phi=psi/(2*pi), T1=t0*psi; the actual A/B formulas vanish identically. At Delta=eta, gamma=eta>0 and a>0, so the square-root factor is smooth while sigma(1-Delta/eta) is flat. Every q jet vanishes there. The zero lazy extension therefore has the same finite jets, uniformly in phase.'
        ],
        numeric_parameter_values_used=False,quantitative_inverse_conditioning_claimed=False)


def _seam_domain(nodes: Any, loop: Mapping[str, Any]) -> dict[str, Any]:
    """Prove the active domain from actual root functions, for all slow Z."""
    import sys
    from fractions import Fraction
    from .original_reference_viscous_m15 import PhysicalEndpointFingerprint

    charts=loop['charts']
    sys.setrecursionlimit(max(sys.getrecursionlimit(),12000))
    f=PhysicalEndpointFingerprint(nodes,dict(coordinate=charts['reference']['coordinate']))
    scope={charts['reference']['coordinate']:f.one,charts['restoration']['coordinate']:f.zero}
    rows=[]
    expected=dict(a=Fraction(4,5),b=Fraction(0),t0=Fraction(0),Delta=Fraction(-6,5))
    for name in ('reference','restoration'):
        for key,value in expected.items():
            node=charts[name]['roots'][key][0]
            digest,rational=f.at(node,scope)
            rows.append(dict(chart=name,quantity=key,original_node=node,
                exact_value=[value.numerator,value.denominator],
                canonical_endpoint_sha256=digest,passed=rational==value))
    return dict(passed=all(row['passed'] for row in rows),rows=rows,
        method='exact_positive_quotient_polynomial_endpoint_identity',
        all_slow_Z_function_identity=True,numeric_Z_values_used=False,
        ordinary_Z1_domain_justification='These rational equalities hold as functions before any Z substitution; source C5 regularity permits differentiation, giving constant root Z1=0 and a strict active gap for every |Z|<=1',
        native_endpoints={'reference':1,'restoration':0},
        active_branch='a=4/5>0, Delta=-6/5<eta, gamma=2*eta+6/5>0; eta,d_star,E positive by pinned source',
        q_positive='sigma(1-Delta/eta)=1, q^2=(5*eta+3)/4>0',
        exact_polynomial_maximum_terms=f.maximum_terms)


def build_phase_lift_certificate(
    nodes: Any,
    l: Mapping[str, Any],
    source_guards: Mapping[str, Any] | None,
    root: Any = None,
) -> dict[str, Any]:
    """Return a strict, source-bound certificate for the fractional phase lift."""
    guards = _mapping(source_guards)
    loop = _mapping(l.get("prescribed_loop_C5_bindings")) if isinstance(l, Mapping) and "prescribed_loop_C5_bindings" in l else _mapping(l)
    repo = _repo_root(root)
    loop_path = repo / _LOOP_SOURCE
    pin = _mapping(_mapping(guards.get("source_pins")).get("loop"))
    pin_evidence = _source_pin(loop_path, pin, repo)

    guard_checks = _mapping(guards.get("checks"))
    source_guard_rows = {
        name: bool(_mapping(guard_checks.get(name)).get("passed") is True)
        for name in _REQUIRED_SOURCE_CHECKS
    }
    source_guard_passed = bool(guards.get("passed") is True and all(source_guard_rows.values()))

    loop_evidence = _function_evidence(loop_path, "loop_functions", (
        "gamma=a.sub(a.fixed(g.mul(c(2),eta)),roots['Delta'])",
        "cut=a.sub(one,a.scale(roots['Delta'],g.quotient(g.one,eta,'same positive canonical eta')))",
        "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))",
        "u=a.scale(a.mul(roots['p2'],q),g.quotient(g.one,dstar,'same positive canonical d_star'))",
        "hinv=a.div(one,a.sqrt(a.add(one,a.mul(u,u))),'positive sqrt(1+u^2)')",
        "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)",
        "numerator=m.add(cc,m.scale(rr,-1))",
        "den=m.add(m.constant(1),m.scale(m.mul(rr,cc),-2),m.mul(rr,rr))",
        "t=m.add(m.axial(roots['t0']),m.mul(m.axial(alpha),m.div(numerator,den)))",
        "ratio=a.add(one,a.mul(roots['t0'],roots['t0']),a.scale(a.mul(q,q),2))",
        "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
        "t2=m.mul(t,t)",
        "integrate=lambda ref:controls.integral(g,ref,g.nodes[angle.node]['name'],g.zero,angle,measure='original angle dpsi with Z and native radius held fixed',true_loop_C5_original_integrand=True,angle_endpoints_independent_of_Z=True)",
        "T1={(k,0):integrate(t[k,0]) for k in range(6)}",
        "T2=FunctionJet(tuple(g.mul(c(math.factorial(k)),integrate(t2[k,0])) for k in range(6)))",
        "P=a.add(a.fixed(angle),T2);Phi=a.mul(K,P)",
        "lam=m.mul(m.axial(K),m.add(m.constant(1),t2))",
        "chi=a.sub(a.fixed(phase),a.scale(psi,g.quotient(g.one,twopi,'positive 2*pi')))",
        "M=a.sub(a.neg(a.scale(a.mul(roots['a'],T1hat),g.quotient(g.one,twopi,'positive 2*pi'))),a.scale(roots['b'],phase))",
        "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))",
        "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
        "Aphi=lazy(a.scale(a.mul(roots['a'],a.sub(one,a.scale(invlambda,g.quotient(g.one,twopi,'positive 2*pi')))),'1/2'))",
        "Bphi=lazy(a.scale(a.mul(roots['E'],a.sub(a.neg(a.scale(a.mul(roots['a'],a.mul(that,invlambda)),g.quotient(g.one,twopi,'positive 2*pi'))),roots['b'])),'1/2'))",
    )) if loop_path.is_file() else {"passed": False, "scope": "loop_functions"}
    roots_evidence = _function_evidence(loop_path, "generic_roots", (
        "C=a.sub(E,a.scale(a.trim(radial['E_y']),2));B=a.scale(a.trim(radial['V_y']),2)",
        "divide=lambda jet:a.div(jet,E,'same strictly positive actual prescribed swirl E')",
        "aa,bb=divide(C),divide(B)",
        "t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))",
    )) if loop_path.is_file() else {"passed": False, "scope": "generic_roots"}
    source_ast_passed = bool(loop_evidence.get("passed") is True and roots_evidence.get("passed") is True)

    phase_passed, phase_rows = _phase_rows(nodes, loop)
    from .original_phase_graph_guards import build_phase_graph_guards
    graph_guards = build_phase_graph_guards(nodes, loop)
    seam_domain = _seam_domain(nodes, loop)
    identities = _formal_identities()
    charts = _mapping(loop.get("charts"))
    source_handles = {
        name: {
            key: _mapping(chart).get(key)
            for key in ("actual_phase", "global_offset", "coordinate", "radius", "physical_measure")
        }
        for name, chart in charts.items()
        if isinstance(chart, Mapping)
    }
    loop_outputs = {
        name: {
            key: _mapping(_mapping(chart).get("loop")).get(key)
            for key in ("Phi", "phase_lambda_at_inverse", "A", "B_over_Pstar", "A_phi", "B_over_Pstar_phi")
        }
        for name, chart in charts.items()
        if isinstance(chart, Mapping)
    }

    passed = bool(
        source_guard_passed
        and pin_evidence["passed"]
        and source_ast_passed
        and phase_passed
        and graph_guards.get("passed") is True
        and seam_domain.get("passed") is True
        and identities.get("passed") is True
    )
    return {
        "passed": passed,
        "method": "pinned_prescribed_loop_poisson_smooth_periodic_lift",
        "scope": "reference native1/restoration native0; universal active recipe and q-flat boundary",
        "global_periodic_lift_proved": passed,
        "mixed_phase_parameter_jets_through_order2_proved": passed,
        "graph_guards": graph_guards,
        "physical_seam_domain": seam_domain,
        "source_only": True,
        "numeric_caps_used": False,
        "B_y_alias_authorized": False,
        "source_guard": {
            "passed": source_guard_passed,
            "required_checks": source_guard_rows,
            "source_guards_passed": guards.get("passed") is True,
            "loop_pin": pin_evidence,
        },
        "source_ast": {
            "passed": source_ast_passed,
            "loop_functions": loop_evidence,
            "generic_roots": roots_evidence,
        },
        "phase_graph": {
            "passed": phase_passed,
            "common_N": loop.get("N"),
            "charts": phase_rows,
            "handles": source_handles,
            "loop_outputs": loop_outputs,
            "uses_fractional_part_of_N_global_offset": phase_passed,
            "uses_phase_center_alias": False,
        },
        "formal_identities": identities,
        "theorem": {
            "fundamental_inverse_domain": "[0,1] phase corresponds to [0,2*pi] angle",
            "global_lift": "Phi(psi+2*pi)=Phi(psi)+1",
            "endpoint_values": "A(0)=A(1)=B(0)=B(1)=0",
            "endpoint_parameter_mixed_jets": "actual A/B shift differences vanish on an open parameter domain; all finite mixed derivatives agree",
            "ordinary_Z1_of_physical_y": "smooth periodic composition admits mixed derivatives; actual native-derivative graph naturality remains a separate obligation",
        },
        "limits": {
            "full_modified_Uy_Qy_RQy_identity": False,
            "B_y_alias_authorized": False,
            "numeric_phase_value_alias": False,
            "phase_center_alias": False,
        },
    }


__all__ = ["build_phase_lift_certificate"]
