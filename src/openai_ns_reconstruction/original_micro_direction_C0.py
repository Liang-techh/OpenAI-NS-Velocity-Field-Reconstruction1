"""Actual two-micro comparison direction and finite-N C0 error radii.

The positive comparison direction is bounded from its own six histories.
The second comparison is never replaced by the continued core function.
"""
from __future__ import annotations
import gzip
import json
from pathlib import Path

from .original_core_signed_direction import require, sha, digest, record, read, source_guard

BASE = "experiments/root_st073/"
PREFIX = BASE + "lei_ren_part1_paper_"
PARENT = PREFIX + "compliant_current_original_core_signed_direction.json"
PARENT_SHA256 = "18968d9762341c53e0f161193d2e998126abb4cd22d73e3df393a62d1cb4d12e"
HELPER = "src/openai_ns_reconstruction/original_loop_C0_kernel.py"
HELPER_SHA256 = "d7d3c0d3eb28b400e401e63ff788e11038dd5bc25bfa0f1fca355c53d44ffc17"
BRIDGE = PREFIX + "compliant_current_original_prescribed_bridge_C6_functions"
INLET = PREFIX + "compliant_current_original_core_inlet_axial6_generic_source"
LOOP = PREFIX + "compliant_current_original_prescribed_loop_C5_functions"
CORE = PREFIX + "compliant_current_original_n1_core_fixed_point_functions.py"
CORE_CHECK = CORE[:-3]+"_check.json"
UNIFORM = PREFIX + "compliant_core_uniform_bounds.json"
NORMS = PREFIX + "compliant_physical_norm_family.json"
SIGMA = PREFIX + "compliant_current_original_n1_Gamma_future_functions.py"
RATES_SOURCE = PREFIX + "compliant_actual_bridge_integrals.py"


def load(path):
    b = Path(path).read_bytes()
    return json.loads(gzip.decompress(b) if str(path).endswith(".gz") else b)


def comparison_algebra():
    import sympy as s
    z, de, rho, eps, lam, g, ph, f, fz = s.symbols("Z delta rho epsilon Lambda g phi F0 F0Z", nonzero=True)
    H, HZ, K, KZ, M, MZ = s.symbols("JH JHZ JK JKZ M MZ")
    d, L = 1-z*z, 1-de*z*z
    W = 1-(1-de)*z*M-d*MZ
    angular = (1-de/2)*f*H-(1-de)*z*(f*HZ+fz*H)/2-d*(f*KZ+fz*K)+(2*de-1)*z*f*K
    A0 = (1-de/2)*H-(1-de)*z*HZ/2-d*KZ+(2*de-1)*z*K
    B0 = (1-de)*z*H/2+d*K
    original = eps*rho/L*(-W+angular/(2*f*ph))
    reduced = rho/L*(g*B0/(2*ph)+eps*(-W+A0/(2*ph)))
    residual = s.cancel((original.subs(fz, -lam*g*f)-reduced).subs(lam, 1/eps))
    require(residual == 0, "Comparison F0 dressing or physical radius differs")
    chi, Dbar, w = s.symbols("chi Dbar w", nonzero=True)
    Ey = w*(ph/2-chi*Dbar*ph/2)
    require(s.cancel((w*ph-2*Ey)/(w*ph)-chi*Dbar) == 0,
            "Actual root a must use the physical log-radius derivative")
    return dict(passed=True, F0_dressing_and_epsilon_Lambda_identity=True,
                exact_reduced_formula="R*Dcmp=rho/L*(gaxis*B0/(2*barphi)+epsilon*(-W+A0/(2*barphi)))",
                Pstar_absent=True, actual_root_a_equals_chi_R_Dbar=True, second_comparison_equals_core_assumed=False,
                history_bound_proof="positive exp(-rate*y) memory plus positive kernel; forcing/rate and inlet in [lo,hi] imply history in [lo,hi]",
                ordinary_Z_history_bound="Z-independent rates, endpoints and hb; differentiate source integrands, then use the same positive kernel")


def comparison_bounds(c, parent, uniform):
    a = parent["actual_core_signed_direction"]
    require(a["passed"] and a["signed_actual_Phi_rho_negative"] and a["micro_core_argument_domain"]["passed"],
            "Accepted continuous core sign and both-micro radial coverage required")
    p = a["parameters"]; j, de, logP = [read(c, p[k]) for k in ("j", "delta", "logPstar")]
    eps = c.exp(-(4*logP+1000)); z = c.mpf([c.mpf(37)/100, c.mpf(93)/250])
    rho = c.mpf([4, c.mpf(41)/10]); hbcap = c.mpf(1)/1000
    tube = a["cases"]["core_tube"]
    ph = read(c, tube["Phi"]); pr = read(c, tube["Phi_rho"])
    require(ph.a > c.mpf(1)/4 and ph.b < 1 and pr.b < 0 and a["micro_core_argument_domain"]["hb_strict_upper"] == "1/1000",
            "Actual core monotonicity, positive floor and source width required")
    rows = {(r["scaled_radial_order"], r["axial_order"]): r for r in uniform["derivative_bounds"]}
    bound = lambda key, field: read(c, rows[key][field]).b
    err = lambda key: bound(key, "Uz_minus_U0_absolute_bound")
    Ucore = 4*z+j+c.mpf([-1,1])*err((0,0))
    variation = 2*hbcap*rho.b*err((1,0))
    Ubar = Ucore+c.mpf([-1,1])*variation
    UZ = 4+err((0,1))+2*hbcap*rho.b*err((1,1))
    require(Ubar.a > c.mpf(7)/5 and Ubar.b < 2 and UZ.b < 5,
            "Positive comparison axial field and ordinary derivative bound unresolved")
    global_lo, global_hi = [read(c, uniform[k]) for k in ("actual_Phi_lower", "actual_Phi_upper")]
    require(global_lo.a > c.mpf(1)/4 and global_hi.b < 2 and Ucore.a > c.mpf(7)/5 and Ucore.b < 2,
            "Same-source core moment inlets must obey the history invariants")
    # alpha in [0,1] and Phi_rho<0 imply corePhi(rho(s)) <= barphi(s) <= corePhi(4).
    # Its ordinary Z derivative is bounded from its defining integral, not
    # by differentiating a cap or replacing the comparison by the core.
    C, C1 = bound((0,1), "Phi_absolute_bound"), bound((1,1), "Phi_absolute_bound")
    logphiZ = C/ph.a+2*hbcap*rho.b*(C1/ph.a+c.fabs(pr).b*C/ph.a**2)
    barphiZ = ph.b*logphiZ
    HZ = max(C, barphiZ.b)
    KZ = max((2*C+10).b, (2*barphiZ+5*ph.b).b)
    H = c.mpf([c.mpf(1)/4,2]); K = c.mpf([c.mpf(7)/20,4]); M = c.mpf([c.mpf(7)/5,2])
    require(min(global_lo.a*Ucore.a, ph.a*Ubar.a) > K.a and
            max(global_hi.b*Ucore.b, ph.b*Ubar.b) < K.b,
            "K comparison moment positivity cannot be inferred without its positive densities")
    d, L = 1-z**2, 1-de*z**2
    H0 = (1-de)*z/2+(1-z**2)*(4*z+j); sig = j/500
    gaxis = L*H0/(H0**2+sig**2)
    B0 = (1-de)*z*H/2+d*K
    Wabs = 1+c.fabs(1-de).b*c.fabs(z).b*M.b+c.fabs(d).b*5
    Aabs = c.fabs(1-de/2).b*H.b+c.fabs((1-de)*z/2).b*HZ+c.fabs(d).b*KZ+c.fabs((2*de-1)*z).b*K.b
    small = eps*(Wabs+Aabs/(2*ph.a))
    leading = gaxis*B0/(2*ph)
    direction = rho/L*(leading+c.mpf([-1,1])*small.b)
    require(gaxis.a > 0 and B0.a > 0 and L.a > 0 and direction.a > 1,
            "Actual comparison direction positivity failed")
    values = dict(epsilon=eps, Phi_comparison=ph, U_comparison=Ubar, U_comparison_Z_absolute=UZ,
                  global_core_Phi_lower=global_lo, global_core_Phi_upper=global_hi,
                  core_Phi_Z_absolute=C, core_Phi_rhoZ_absolute=C1,
                  comparison_logPhi_Z_absolute=logphiZ, comparison_Phi_Z_absolute=barphiZ,
                  JH=H, JK=K, M=M, JH_Z_absolute=HZ, JK_Z_absolute=KZ, M_Z_absolute=c.mpf(5),
                  gaxis=gaxis, B0=B0, A0_absolute=Aabs, W_absolute=Wabs,
                  epsilon_error_absolute=small, positive_leading=leading, R_D_comparison=direction)
    return dict(passed=True, domain=dict(s=["0","2"], Z=["37/100","93/250"]),
                bounds={k:record(v) for k,v in values.items()},
                two_micro_comparison_positive=True, source_bound_from_own_comparison_histories=True,
                core_and_comparison_fields_identified_on_second_micro=False,
                Phi_comparison_ordering="Phi_core(rho(s),Z)<=barphi(s,Z)<=Phi_core(4,Z)",
                actual_native_phi_positive_and_nonincreasing=True,
                actual_root_a_formula="a=chi_R*(R*D_comparison), hb<=chi_R<=1",
                core_perturbation_errors_are_bounds_not_signed_values=True)


def amplitude_and_errors(c, root, parent, norms, nodes):
    p = parent["actual_core_signed_direction"]["parameters"]
    j, de, logP = [read(c,p[k]) for k in ("j","delta","logPstar")]
    logLambda = 4*logP+1000
    t = c.mpf([-j.b,c.mpf(93)/250])
    Hp = (1-de)/2+4-12*t**2-2*j*t
    logC = read(c,norms["selected_logCstar"])
    identity = parent["actual_native_core_angular_direction"]["angular_direction_identity"]["nonvanishing_domain"]
    require(Hp.a > 0 and logC.a > 0 and identity["F0_positive_bound"],
            "Anchored gauge sign or actual selected positive logC missing")
    # The exact source anchor is the H0 root in [-j,0].  H0 increases
    # from that root to every Z of this tube, so the oriented G integral
    # is nonnegative and F0=exp(-logC-Lambda*G)<=1.
    core_receipt = load(root/CORE_CHECK)
    scalar = core_receipt["checks"]["actual_oracle"]["exact_parameter_and_native_P0_identity"]["actual_source_scalar_and_unique_anchor_identity"]
    require(scalar["unique_H_root_in_original_bracket"] and scalar["logC_matches_original_selected_source"] and
            scalar["anchor_leaf_uses_this_exact_source"], "Consume the same actual anchor/logC identity")
    logeta = read(c,nodes[21]["canonical_parameter_proof"]["original_canonical_eta_log"])
    logdstar = read(c,nodes[22]["canonical_parameter_proof"]["original_canonical_dstar_log"])
    require(logeta.b < 0 and nodes[21]["name"] == "eta" and nodes[22]["name"] == "d_star",
            "Actual canonical eta/dstar logs required")
    require(nodes[21]["canonical_parameter_proof"]["exact_original_parameter_tuples_equal"] and
            nodes[22]["canonical_parameter_proof"]["exact_original_parameter_tuples_equal"],
            "The eta/dstar logs must belong to the actual unchanged parameters")
    logN = read(c,parent["actual_micro_phase_mean"]["directed_physical_prefactor"]["log_common_N"])
    require(parent["actual_signed_lower_m15"]["source_log_basis_names"][0] == "hb", "Joint mean needs the actual hb log")
    loghb = read(c,parent["actual_signed_lower_m15"]["source_log_bases"][0])
    rho_upper = read(c,parent["actual_core_signed_direction"]["micro_core_argument_domain"]["core_argument_upper"])
    ph_upper = read(c,parent["actual_core_signed_direction"]["cases"]["core_tube"]["Phi"]).b
    logPE = c.ln(ph_upper)+(c.ln(2)+c.ln(rho_upper.b)-logLambda)/2
    logs = dict(relative_delta_E=c.ln(4)-logN, relative_delta_V=c.ln(2)-logN,
                physical_Pstar_delta_E=c.ln(4)+logPE-logN,
                physical_Pstar_delta_V=c.ln(2)+logPE-logN,
                physical_joint_delta_m=c.ln(4)+loghb+logPE-logN)
    threshold = -1024*c.ln(2)
    require(all(v.b < threshold.a for v in logs.values()), "Requested C0 roundoff radius not certified")
    return dict(passed=True, H0_derivative_on_anchor_path=record(Hp), actual_positive_logC=record(logC),
                F0_positive_and_at_most_one=True, exact_anchored_G_nonnegative=True,
                log_eta=record(logeta), eta_positive_and_below_one=True,
                log_dstar=record(logdstar), dstar_strict_positive=True,
                log_common_N=record(logN), log_Pstar_E_upper=record(logPE),
                logarithmic_error_radii={k:record(v) for k,v in logs.items()},
                common_coarse_error_radius="2^-1024", coarse_radius_positive=True,
                ordinary_Z_order=0, phase_free_uniform_bounds=True,
                exact_graphs_retained=True, zero_is_error_box_center_not_function_value=True,
                phase_N_hb_F0_not_materialized=True,
                joint_mean_bound="abs(delta_m)<=4*Pstar*hb*E_upper/N",
                joint_mean_incoming_history_retained=True, background_mean_not_zeroed=True,
                numerical_Z1_or_Z3_error_enclosure=False)


def build(root, parent):
    from mpmath.ctx_iv import MPIntervalContext
    root = Path(root)
    require(sha(root/PARENT) == PARENT_SHA256 and load(root/PARENT) == parent, "Accepted signed core parent required")
    require(sha(root/HELPER) == HELPER_SHA256, "Reviewed new loop C0 helper changed")
    paths = [BRIDGE+s for s in (".py",".json.gz","_check.json")]+[INLET+s for s in (".py",".json.gz","_check.json")]+[
        LOOP+s for s in (".py",".json.gz","_check.json")]+[CORE,CORE_CHECK,UNIFORM,NORMS,SIGMA,RATES_SOURCE]
    pins = {PARENT:dict(path=PARENT,expected_sha256=PARENT_SHA256),
            HELPER:dict(path=HELPER,expected_sha256=HELPER_SHA256)}
    for path in paths:
        h = parent["input_hashes"][Path(path).name]
        require(sha(root/path) == h, "Changed actual micro source: "+path)
        pins[path] = dict(path=path,expected_sha256=h)
    import ast
    rate_assignment = next(n.value for n in ast.parse((root/RATES_SOURCE).read_text()).body
                           if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id == "RATES" for t in n.targets))
    require(isinstance(rate_assignment,ast.Call) and isinstance(rate_assignment.func,ast.Name) and rate_assignment.func.id == "dict" and
            {k.arg:ast.literal_eval(k.value) for k in rate_assignment.keywords} == dict(H=2,M=1,K=2,A=1,B=2,C=1) and
            "RATES=bridge.RATES" in "".join((root/(BRIDGE+".py")).read_text().split()),
            "Actual comparison history rates changed")
    bridge, inlet, loop = [load(root/(s+".json.gz")) for s in (BRIDGE,INLET,LOOP)]
    for s in (BRIDGE,INLET,LOOP):
        r = load(root/(s+"_check.json")); require(r["all_passed"] and r["source_family"] == parent["source_family"], "Checked same source required")
    bn, inn, nodes = bridge["exact_graph_nodes"], inlet["exact_graph_nodes"], loop["exact_graph_nodes"]
    require(nodes[:len(bn)] == bn and bn[:len(inn)] == inn, "One actual graph prefix, not just source-family equality")
    bb, ib, lb = bridge["prescribed_bridge_C6_bindings"], inlet["core_inlet_axial6_bindings"], loop["prescribed_loop_C5_bindings"]
    require(all(bb["core7"][k][:7] == ib["core_functions"][k] for k in ib["core_functions"]) and
            all(bb["initial7"][k][:7] == ib["inlet_moments"][k] for k in ib["inlet_moments"]), "Unreset actual core inlets required")
    require(lb["charts"]["first_micro"]["roots"] == lb["charts"]["second_micro"]["roots"] and
            lb["N"] == 18 and lb["eta"] == 21 and lb["dstar"] == 22 and
            nodes[18] == parent["actual_micro_correlated_m15"]["common_N"], "Same roots and common actual parameters required")
    first, second = [lb["charts"][k] for k in ("first_micro","second_micro")]
    require(bb["micro"]["domain"] == [0,2] and bb["micro"]["subcharts"] == dict(first=[0,1],second=[1,2]) and
            bb["micro"]["physical_measure"] == first["physical_measure"] == second["physical_measure"] == 8 and
            bb["micro"]["radius"] == first["radius"] == second["radius"] and ib["Pstar"] == 33572 and
            nodes[33572] == dict(operation="analytic_unary",name="exp",argument=5) and
            second["incoming"] == lb["endpoints"]["first_micro"] and
            all(row == 0 for row in first["incoming"]["m"]),
            "Physical Pstar/hb and original nonreset correction histories must be the actual source")
    shear = second["loop"]["flat_predicate"]
    require(nodes[shear] == dict(operation="exact_real_comparison",operator=">=",left=second["roots"]["Delta"][0],right=21),
            "The actual shear-flat branch must be Delta>=eta")
    for chart in (first,second):
        predicate = chart["loop"]["flat_predicate"]
        require(predicate == shear or (nodes[predicate]["operation"] == "logical_or" and shear in nodes[predicate]["arguments"]),
                "Every active micro branch must imply Delta<eta")
        for field in ("A","B_over_Pstar","A_phi","B_over_Pstar_phi"):
            row = nodes[chart["loop"][field][0]]
            require(row["operation"] == "original_lazy_flat_branch" and row["flat_value"] == 0 and
                    row["flat_predicate"] == predicate and row["active_body_not_evaluated_on_flat_branch"],
                    "Actual A/B flat branches must be exactly zero")
    require(nodes[35563] == dict(operation="analytic_unary",name="exp",argument=35562) and
            nodes[35565] == dict(operation="analytic_unary",name="exp",argument=35564) and
            nodes[35564] == dict(operation="negative",argument=35562) and
            nodes[35653] == dict(operation="negative",argument=35649) and
            nodes[35649]["operation"] == "product" and set(nodes[35649]["arguments"]) == {35563,bb["core7"]["gradient"][0]} and
            nodes[35658]["operation"] == "product" and set(nodes[35658]["arguments"]) == {35657,35653} and
            bb["core7"]["logF0"][1] == 35653 and bb["core7"]["F0"][:2] == [35657,35658] and
            nodes[ib["radius"]] == dict(operation="product",arguments=[13,35565]) and
            nodes[13] == dict(operation="exact_rational",numerator=4,denominator=1),
            "Actual F0_Z/F0=-Lambda*gaxis and R=epsilon*rho source nodes differ")
    guards = source_guards(root)
    c = MPIntervalContext(); c.prec = 600
    uniform = load(root/UNIFORM); cb = comparison_bounds(c,parent,uniform)
    errors = amplitude_and_errors(c,root,parent,load(root/NORMS),nodes)
    from .original_loop_C0_kernel import symbolic_controls, directed_log_bounds
    controls = symbolic_controls()
    logs = directed_log_bounds(c,read(c,errors["log_common_N"]))
    return dict(passed=True, source_family=parent["source_family"], source_pins=pins, source_AST=guards,
                comparison_algebra=comparison_algebra(), comparison_direction=cb, actual_C0_errors=errors,
                universal_loop_controls=controls, directed_relative_errors=logs,
                loop_cut_variable_binding=dict(helper_cut="sigma_flat(1-Delta/eta)",
                    source_cut_argument="1-Delta/eta",helper_cut_range=[0,1],
                    source_argument_itself_not_assumed_below_one=True),
                actual_root_domain_gates=dict(passed=True,s=[0,2],hb_strictly_positive_and_below_one=True,
                    sigma_flat_range=[0,1],chi_R_formula="1+(hb-1)*sigma_flat(s)",
                    chi_R_in_positive_interval="[hb,1]",positive_comparison_implies_native_phi_nonincreasing=True,
                    active_implies_Delta_less_eta=True,flat_implies_A_and_B_exact_zero=True),
                actual_micro_source_handles=dict(comparison=bb["comparison_micro"], native=bb["micro"],
                    first=lb["charts"]["first_micro"], second=lb["charts"]["second_micro"]),
                both_micro_actual_a_and_E_positive=True, actual_finiteN_and_joint_mean_C0_enclosed=True,
                ordinary_Z_orders=[0], full_C5_numeric_containment=False, numerical_Cy_Z3=False,
                numerical_QZZ=False, forward_Rp_discrepancy_enclosed=False, global_five_moment_closure=False,
                old_source_packets_and_all_errors_unchanged=True, accepted_ancestors_replayed=False, full_goal_complete=False)


def source_guards(root):
    return dict(
        comparison=source_guard(root,BRIDGE+".py","__init__",(
            "alpha=g.sub(g.one,sigma(g.sub(u,g.one)))", "rho=lambda q:g.mul(c(4),g.unary('exp',g.mul(hb,q)))",
            "barJ=histories(a,barphi,barU,self.initial7,s,hb,'comparison_micro')",
            "chi=g.add(g.one,g.mul(g.sub(hb,g.one),sigma(s)))",
            "R=g.mul(Ra,g.unary('exp',g.mul(hb,s)))")),
        comparison_densities=source_guard(root,BRIDGE+".py","densities",(
            "return dict(H=a.scale(phi,2),M=U,K=a.scale(a.mul(phi,U),2)" ,)),
        histories=source_guard(root,BRIDGE+".py","histories",(
            "out[name]=a.add(a.scale(initial[name],decay),a.integral(body,t,g.zero,phase,",
            "original_incoming_memory_retained=True,source_measure_applied_once=True")),
        directions=source_guard(root,BRIDGE+".py","comparison_directions",(
            "f,h,k=(b.mul(F0,jet) for jet in (phi,J['H'],J['K']))", "W=a.sub(one,a.add(a.mul(a.sub(one,de),a.mul(z,M)),a.mul(d,Mz)))",
            "D=a.div(a.add(a.neg(W),a.div(angular,a.scale(a.trim(f),2),")),
        core_inlets=source_guard(root,INLET+".py","__init__",(
            "densities=dict(H=a.scale(phi,g.mul(c(2),t)),M=u,K=a.scale(a.mul(phi,u),g.mul(c(2),t))",
            "self.inlet_moments={key:a.at(jet,rho,c(4)) for key,jet in self.moment_functions.items()}")),
        core_gauge=source_guard(root,CORE,"__init__",(
            "gradient=a.div(a.mul(L,H),a.add(a.mul(H,H),a.fixed(g.mul(sigma,sigma))))",
            "logF=a.add(a.fixed(g.neg(self.parameters['logC'])),a.neg(a.scale(Gjet,Lambda)))",
            "Gjet=c3.C3Function(G,gradient.value,gradient.Z,gradient.ZZ)",
            "G=controls.integral(g,gv,g.nodes[v.node]['name'],self.parameters['anchor'],z,")),
        native_E=source_guard(root,BRIDGE+".py","generic_fields",(
            "W=a.scale(F0,g.mul(g.unary('sqrt',g.mul(c(2),radius)),inverse))", "fields=dict(E=a.mul(W,actual['phi'])",
            "radial=dict(E_y=a.mul(W,a.add(a.scale(actual['phi'],'1/2'),actual['phi_y']))")),
        native_phi=source_guard(root,BRIDGE+".py","actual_fields",(
            "phi=a.mul(initial['phi'],a.exp(logratio))",
            "phiy=a.scale(a.mul(known['D_over_R'],phi),g.mul(chi,radius,g.constant('-1/2')))")),
        sigma=source_guard(root,SIGMA,"sigma",(
            "definition='0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside'",)),
        loop_roots=source_guard(root,LOOP+".py","generic_roots",(
            "C=a.sub(E,a.scale(a.trim(radial['E_y']),2))", "aa,bb=divide(C),divide(B)",
            "kap=a.add(aa,a.div(a.mul(bb,bb),aa,'same original positive a'))")),
        loop_bodies=source_guard(root,LOOP+".py","loop_functions",(
            "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))",
            "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))", "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
            "r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)",
            "K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')",
            "Aphi=lazy(a.scale(a.mul(roots['a'],a.sub(one,a.scale(invlambda,g.quotient(g.one,twopi,'positive 2*pi')))),'1/2'))",
            "Bphi=lazy(a.scale(a.mul(roots['E'],a.sub(a.neg(a.scale(a.mul(roots['a'],a.mul(that,invlambda)),")),
        loop_sigma=source_guard(root,LOOP+".py","sigma",("value=previous.sigma(g,rows(x)[0])",)),
        lazy_flat=source_guard(root,LOOP+".py","lazy",(
            "flat_value=g.zero.node", "active_body_not_evaluated_on_flat_branch=True")),
        finiteN=source_guard(root,LOOP+".py","finite_N_functions",(
            "r=a.mul(roots['E'],exm1);s=a.scale(loop['B_over_Pstar'],invN)",)),
        joint_mean=source_guard(root,LOOP+".py","__init__",(
            "incoming={key:a.fixed(0) for key in RATES}", "incoming=self.endpoints[name]",
            "lower=g.mul(self.sc,g.constant('1/2')) if name=='first_micro' else g.one if name=='second_micro' else g.zero")))
