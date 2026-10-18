"""Same-source fixed-phase inverse and loop mixed logarithmic bounds.

All rows are ordinary derivatives.  The y derivative is physical slow y at
fixed target phase; the N times phase chain term is kept separate.  No
unknown inverse angle, signed source value, or finite-N phase is selected.
"""
from __future__ import annotations

from pathlib import Path

from .original_core_signed_direction import require, sha, read, record, source_guard
from .original_micro_direction_C0 import load, PREFIX, LOOP
from .original_log_mixed_jets import LogMixed
from .original_log_phase_jets import LogPhase, inverse_coefficients, compose_inverse

PARENT = PREFIX + "compliant_current_original_micro_q_mixed_bounds.json"
PARENT_SHA256 = "a28ae027174ceefa92c80c7f7437689ae964c83ca5c116c4c088a2ec156d2c39"
HELPER = "src/openai_ns_reconstruction/original_log_phase_jets.py"
HELPER_SHA256 = "703b2bbd1ba12fc2e82486fe5244a7cfee95749d2f486d35cfa6772f34958401"


def pack(jet):
    return [dict(physical_y_order=i, ordinary_Z_order=k, exact_zero=v is None,
                 log_absolute_upper_bound=None if v is None else record(v))
            for (i, k), v in jet.items()]


def phase_pack(jet):
    return [dict(physical_y_order=i, ordinary_Z_order=k, free_angle_order=j,
                 exact_zero=v is None,
                 log_absolute_upper_bound=None if v is None else record(v))
            for (i, k, j), v in jet.items()]


def unpack(c, rows):
    values = {(r["physical_y_order"], r["ordinary_Z_order"]):
              None if r["exact_zero"] else read(c, r["log_absolute_upper_bound"]) for r in rows}
    return LogMixed.from_mapping(c, values)


def roots(c, parent, chart, loghb):
    """Consume saved derivative rows; only active value rows are tightened."""
    ea = parent["actual_micro_Ea_mixed_bounds"]["actual_Ea_mixed_bounds"]
    bp = parent["actual_micro_b_p2_mixed_bounds"]["actual_b_t0_p2_mixed_bounds"]
    aa, ee = {}, {}
    for k in range(4):
        cap = read(c, ea["a_Z_absolute_factored_chi_R"][k]["absolute_upper_enclosure"])
        aa[0, k] = c.ln(cap.b) + (loghb if chart == "second" else 0)
        aa[1, k] = read(c, ea["a_yZ_log_absolute_upper_" + chart][k]["log_absolute_upper_bound"])
        for i, key in ((0, "E_Z_over_E"), (1, "E_yZ_over_E")):
            ee[i, k] = c.ln(read(c, ea[key][k]["absolute_upper_enclosure"]).b)
    aa[0, 0] = c.mpf(min(c.ln(3).b, aa[0, 0].b))
    def channel(key, ykey):
        return LogMixed.from_mapping(c, {
            (i, k): read(c, bp[(key if i == 0 else ykey)][k]["log_absolute_upper_bound"])
            for i in (0, 1) for k in range(4)})
    b = channel("b_Z_" + chart, "b_yZ_" + chart)
    bb = b.as_dict(); bb[0, 0] = c.mpf(min(c.ln(2).b, bb[0, 0].b))
    return dict(a=LogMixed.from_mapping(c, aa), b=LogMixed.from_mapping(c, bb),
                t0=channel("t0_Z_both", "t0_yZ_both"),
                E_over_frozen_E0=LogMixed.from_mapping(c, ee))


def primitive(engine, jet):
    """Angle primitive from zero; ordinary rows shift with no 1/j factor."""
    c = engine.c
    return LogPhase.from_mapping(c, {
        (i, k, j): (None if jet[i, k, 0] is None else jet[i, k, 0] + c.ln(2*c.pi))
        if j == 0 else jet[i, k, j-1] for i, k, j in engine.indices})


def endpoint_controls():
    """Exact Fourier normalization, independent of numerical angle sampling."""
    import sympy as s
    r,z,q,h,a,b=s.symbols("r z q h a b", nonzero=True)
    cosine=(z+1/z)/2
    H=(cosine-r)/(1-2*r*cosine+r*r)
    generating=(z/(1-r*z)+1/(z-r))/2
    require(s.cancel(H-generating)==0,"Poisson Fourier generating identity failed")
    # H=sum_{n>=1} r^(n-1)cos(n*psi), absolutely convergent for |r|<1.
    # Orthogonality gives mean H=0 and mean H^2=1/[2(1-r^2)].
    require(s.cancel((2*q*h)**2/(2*h*h)-2*q*q)==0,"Poisson variance normalization failed")
    require(s.simplify(-a*(-b/a)-b)==0,"Original t0=-b/a endpoint gauge failed")
    return dict(passed=True,cosine_generating_identity_exact=True,
        Fourier_series="H=(cos(psi)-r)/D=sum(n>=1,r^(n-1)*cos(n*psi))",
        convergence="absolute for |r|<1, including r=0 by continuation; accepted 1-r^2>0",
        mean_H="0",mean_H_squared="1/(2*(1-r^2))",
        alpha="2*q*hinv",gap="1-r^2=hinv^2",mean_t="t0",mean_t_squared="t0^2+2*q^2",
        T1_end="2*pi*t0",T2_end="2*pi*(t0^2+2*q^2)",
        Phi_zero="0",Phi_twopi="1",Phi_strictly_increasing=True,
        inverse_exists_uniquely_for_every_target_phase=True,
        exact_t0_identity="t0=-b/a",A_B_both_phase_endpoints_exact_zero=True,
        identities_hold_as_parameter_functions_not_at_sampled_points=True,
        normalized_ordinary_derivative_conversion_not_used_in_endpoint_proof=True)


def chart_bounds(c, parent, chart, C, loghb):
    m, p = LogMixed(c), LogPhase(c)
    saved = parent["actual_micro_q_mixed_bounds"]["charts"][chart]
    rr = {key: unpack(c, saved[key + "_mixed_log_caps"]) for key in ("r", "alpha", "K")}
    source = roots(c, parent, chart, loghb)
    r, alpha, K = [p.lift_mixed(rr[key]) for key in ("r", "alpha", "K")]
    cosine = LogPhase.from_mapping(c, {idx: c.mpf(0) if idx[:2] == (0, 0) else None
                                       for idx in p.indices})
    one = p.constant(c.mpf(0))
    denominator = p.add(one, p.scale(p.mul(r, cosine), c.ln(2)), p.mul(r, r))
    cond = saved["scaled_Poisson_conditioning"]
    logDlower = read(c, cond["full_free_angle_D_log_lower"])
    loglamlower = read(c, cond["full_free_angle_lambda_log_lower"])
    require(logDlower.b < denominator[0,0,0].a,"Positive D source floor inconsistent with norm cap")
    invD = p.reciprocal(denominator, logDlower)
    t = p.add(p.lift_mixed(source["t0"]), p.mul(alpha, p.mul(p.add(cosine, r), invD)))
    t2 = p.mul(t, t)
    lam = p.mul(K, p.add(one, t2))
    require(loglamlower.b < lam[0,0,0].a,"Positive lambda source floor inconsistent with norm cap")
    T1, T2 = primitive(p, t), primitive(p, t2)
    angle = LogPhase.from_mapping(c, {idx: c.ln(2*c.pi) if idx == (0, 0, 0)
                                      else c.mpf(0) if idx == (0, 0, 1) else None for idx in p.indices})
    rawPhi = p.mul(K, p.add(angle, T2))
    phi_rows = {idx: rawPhi[idx] if idx[2] == 0 else lam[idx[0], idx[1], idx[2]-1]
                for idx in p.indices}
    # The accepted endpoint normalization and monotonicity imply 0<=Phi<=1.
    phi_rows[0, 0, 0] = c.mpf(0)
    Phi = LogPhase.from_mapping(c, phi_rows)
    inverse, increment, proof = inverse_coefficients(p, Phi, loglamlower)
    T1hat, that, lamhat = [compose_inverse(p, jet, increment) for jet in (T1, t, lam)]
    target = m.constant(c.mpf(0))  # |target phase|<=1, all slow derivatives exactly zero.
    invlamhat = m.reciprocal(lamhat, loglamlower)
    a, b = source["a"], source["b"]
    chi = m.add(target, m.scale(inverse, -c.ln(2*c.pi)))
    A = m.scale(m.mul(a, chi), -c.ln(2))
    Brel = m.scale(m.add(m.scale(m.mul(a, T1hat), -c.ln(2*c.pi)), m.mul(b, target)), -c.ln(2))
    Aphi = m.scale(m.mul(a, m.add(m.constant(c.mpf(0)), m.scale(invlamhat, -c.ln(2*c.pi)))), -c.ln(2))
    Bphirel = m.scale(m.add(m.scale(m.mul(a, m.mul(that, invlamhat)), -c.ln(2*c.pi)), b), -c.ln(2))
    # E0 is the actual positive E at the differentiation point, frozen as a
    # reference only.  This multiplication retains every E_y/E_Z cross term.
    B = m.mul(source["E_over_frozen_E0"], Brel)
    Bphi = m.mul(source["E_over_frozen_E0"], Bphirel)
    jets = dict(inverse=inverse, T1_at_inverse=T1hat, direction_at_inverse=that,
                lambda_at_inverse=lamhat, A=A, B_over_Pstar_over_frozen_E0=B,
                A_phi=Aphi, B_over_Pstar_phi_over_frozen_E0=Bphi)
    require(all(v is not None and v.b<(9000*C).a for jet in jets.values() for _,v in jet.items()),
            "Fixed-phase mixed envelope failed")
    require(all(jet[0,k].b<(5600*C).a for jet in (A,B) for k in range(4)),
            "Fixed-phase ordinary Z loop envelope failed")
    ratios = {name: [dict(physical_y_order=i, ordinary_Z_order=k,
                         log_cap_over_selected_logC=record(value/C))
                     for (i,k),value in jet.items() if value is not None]
              for name,jet in jets.items()}
    return dict(passed=True, ordinary_Z_orders=[0,1,2,3], physical_y_orders=[0,1],
        active_branch_only=True, free_angle=["0", "2*pi"], target_phase=["0", "1"],
        positive_D_log_lower=record(logDlower), positive_lambda_log_lower=record(loglamlower),
        free_angle_jets={name:phase_pack(jet) for name,jet in dict(D=denominator,t=t,T1=T1,T2=T2,Phi=Phi,lambda_value=lam).items()},
        inverse_recurrence=proof, mixed_log_caps={name:pack(jet) for name,jet in jets.items()},
        log_cap_ratios=ratios, E_product_rule_retained=True,
        strict_log_envelopes=dict(all_eight_loop_inverse_fields=9000,
            A_and_B_over_Pstar_over_frozen_E0_Z0_Z3=5600,reference="selected positive logC",passed=True),
        E_reference_is_point_frozen_not_function_division=True,
        bound_signs_are_not_actual_source_signs=True,
        exact_source_formulas=dict(A="a/2*(phase-psi/(2*pi))",
            B_over_Pstar="E/2*(-a*T1hat/(2*pi)-b*phase)",
            A_phi="a/2*(1-1/(2*pi*lambda_hat))",
            B_over_Pstar_phi="E/2*(-a*t_hat/(2*pi*lambda_hat)-b)"))


def build(root, parent):
    root = Path(root)
    require(sha(root/PARENT) == PARENT_SHA256 and load(root/PARENT) == parent, "Accepted q parent required")
    require(sha(root/HELPER) == HELPER_SHA256, "Reviewed implicit jet helper changed")
    old = parent["actual_micro_q_mixed_bounds"]
    require(old["passed"] and old["actual_free_angle_positive_D_and_lambda_log_bounds"]
            and old["natural_collar_flat_proof"]["passed"] and old["shear_flat_glue"]["passed"],
            "Accepted active conditioning and flat glue required")
    from mpmath.ctx_iv import MPIntervalContext
    from .original_log_phase_jets import symbolic_controls
    c = MPIntervalContext(); c.prec = 600
    actual = old["parameter_envelopes"]["actual_logs"]
    C, loghb = read(c, actual["C"]), read(c, actual["hb"])
    pins = {PARENT:dict(path=PARENT,expected_sha256=PARENT_SHA256),
            HELPER:dict(path=HELPER,expected_sha256=HELPER_SHA256)}
    loop_path=LOOP+".py"; loop_hash=parent["input_hashes"][Path(loop_path).name]
    require(sha(root/loop_path)==loop_hash,"Actual C5 phase source changed")
    pins[loop_path]=dict(path=loop_path,expected_sha256=loop_hash)
    ast = source_guard(root,loop_path,"loop_functions",(
        "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)",
        "den=m.add(m.constant(1),m.scale(m.mul(rr,cc),-2),m.mul(rr,rr))",
        "t=m.add(m.axial(roots['t0']),m.mul(m.axial(alpha),m.div(numerator,den)))",
        "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
        "T1={(k,0):integrate(t[k,0]) for k in range(6)}",
        "P=a.add(a.fixed(angle),T2);Phi=a.mul(K,P)",
        "lam=m.mul(m.axial(K),m.add(m.constant(1),t2))",
        "phase_endpoints=[0,1]", "derivative_positive_certificate='K*(1+t^2)>0; same original a/q-flat construction'",
        "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))",
        "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))"))
    root_ast=source_guard(root,loop_path,"generic_roots",(
        "aa,bb=divide(C),divide(B)",
        "p2=generic['p2'];t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))"))
    return dict(passed=True,source_family=parent["source_family"],source_pins=pins,
        source_AST=ast,root_source_AST=root_ast,endpoint_normalization=endpoint_controls(),
        new_implicit_jet_controls=symbolic_controls(),
        charts={chart:chart_bounds(c,parent,chart,C,loghb) for chart in ("first","second")},
        domain=old["domain"],ordinary_Z_orders=[0,1,2,3],physical_y_orders=[0,1],
        derivative_scope="physical slow y and ordinary Z at fixed target phase; not total y",
        flat_branch=dict(passed=True,apply_before_active_formula=True,q_A_B_Aphi_Bphi_exact_zero=True,
            inverse="psi=2*pi*phase",slow_inverse_y_Z_derivatives_exact_zero=True,
            phase_inverse_phi_derivative="2*pi",raw_flat_a_b_not_bounded_by_active_caps=True,
            glue="accepted q C4 flatness, smooth positive-lambda implicit inverse, exact lazy zero loops",
            q_zero_identities="t=t0; lambda=1/(2*pi); Phi=psi/(2*pi); T1=t0*psi; a*t0+b=0",
            loop_jets_zero_by_composition_and_exact_q_zero_identity=True,
            local_positive_D_lambda_smoothness_not_uniform_floor_on_flat_claim=True,
            natural_collar_receipt_consumed=True,shear_flat_receipt_consumed=True),
        total_y_chain_rule="d_total_y = d_slow_y + N*d_phase; phase_Z=0",
        B_scale_convention="source B_over_Pstar=E*M/2; physical B=Pstar*B_over_Pstar",
        phase_or_inverse_mixed_provider_complete=True,loop_AB_fixed_phase_mixed_provider_complete=True,
        finite_N_derivative_errors=False,numerical_Cy_Z3=False,numerical_mean_Z3=False,numerical_QZZ=False,
        full_C5_numeric_containment=False,global_five_moment_closure=False,full_goal_complete=False,
        inherited_signed_sources_and_all_errors_unchanged=True,accepted_ancestors_replayed=False)
