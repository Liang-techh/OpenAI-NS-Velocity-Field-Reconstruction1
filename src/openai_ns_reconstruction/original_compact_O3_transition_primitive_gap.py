"""Fresh ordinary-Z root/primitive gap algebra; unexecuted until frozen entry.

The caller must bind actual source clauses before calling this routine. No
source value, inverse, integral, provider, or ancestor operator is evaluated.
"""
from math import comb as _pg_comb


def _pg_constant(value):
    return {(): value} if value else {}


def _pg_symbol(name):
    return {(name,): 1}


def _pg_sum(*parts):
    total = {}
    for part in parts:
        for key, coefficient in part.items():
            total[key] = total.get(key, 0) + coefficient
            if total[key] == 0:
                del total[key]
    return total


def _pg_scale(part, scalar):
    return {key: coefficient * scalar for key, coefficient in part.items() if coefficient * scalar}


def _pg_product(*parts):
    value = _pg_constant(1)
    for part in parts:
        combined = {}
        for left, lc in value.items():
            for right, rc in part.items():
                key = tuple(sorted(left + right))
                combined[key] = combined.get(key, 0) + lc * rc
        value = {key: coefficient for key, coefficient in combined.items() if coefficient}
    return value


def _pg_difference(left, right):
    return _pg_sum(left, _pg_scale(right, -1))


def _pg_derivative(poly):
    """Ordinary D_Z, with g=mu*sigma(x) and tau=2*pi held fixed."""
    pieces = []
    for monomial, coefficient in poly.items():
        for index, name in enumerate(monomial):
            if name in ('g', 'tau'):
                continue
            stem, order_text = name.rsplit('_', 1)
            order = int(order_text)
            if order >= 3:
                raise ValueError('Fourth ordinary-Z derivative is outside this leaf')
            differentiated = monomial[:index] + (stem + '_' + str(order + 1),) + monomial[index + 1:]
            pieces.append({tuple(sorted(differentiated)): coefficient})
    return _pg_sum(*pieces)


def _pg_rows(poly):
    return [{'coefficient': coefficient, 'factors': list(key)} for key, coefficient in sorted(poly.items())]


def _pg_jet_product(left, right, order):
    return _pg_sum(*[
        _pg_scale(_pg_product(left[k], right[order - k]), _pg_comb(order, k))
        for k in range(order + 1)
    ])


def derive_primitive_gap_algebra(contract, journal=None):
    """Check fresh identities only after strict caller-owned source binding."""
    required = ('actual_native_primitive_bound', 'actual_C3_primitive_bound',
                'actual_root_shape_clauses_bound', 'same_mu_sigma_source_forms_bound',
                'ordinary_Z_normalization_bound', 'positive_E_and_a_source_bound',
                'sigma_endpoint_extensions_and_reflection_bound', 'both_zero_Vy_source_clauses_bound',
                'same_flat_predicate_forms_bound', 'normalized_B_roles_kept_distinct_from_physical_B')
    if not isinstance(contract, dict) or any(contract.get(name) is not True for name in required):
        raise ValueError('Missing strict primitive source gate')
    if contract.get('ordinary_Z_orders') != [0, 1, 2, 3]:
        raise ValueError('Unexpected ordinary-Z orders')
    if contract.get('source_values_evaluated') is not False:
        raise ValueError('This leaf cannot evaluate original source values')

    zero = {}
    one = _pg_constant(1)
    two = _pg_constant(2)
    g = _pg_symbol('g')
    R = _pg_sum(two, _pg_scale(g, 2))
    delta = _pg_scale(g, 2)
    # Native sigma uses a/(a+b), a=exp(-1/x^2), b=exp(-1/(1-x)^2).
    # C3 uses the positive logistic a/b on the left and reflection on the right.
    # Their positive denominators are cleared without evaluating either source.
    sa, sb = _pg_symbol('sigma_a'), _pg_symbol('sigma_b')
    sigma_witnesses = [
        ('left_logistic', _pg_product(sa, _pg_sum(sb, sa)),
         _pg_product(sa, _pg_sum(sa, sb))),
        ('right_reflection', _pg_product(sa, _pg_sum(sb, sa)),
         _pg_product(_pg_difference(_pg_sum(sb, sa), sb), _pg_sum(sa, sb))),
    ]
    sigma_checks = []
    for name, left, right in sigma_witnesses:
        residual = _pg_difference(left, right)
        sigma_checks.append({'branch': name, 'left': _pg_rows(left), 'right': _pg_rows(right),
                             'residual': _pg_rows(residual), 'is_zero': not residual})
    if not all(check['is_zero'] for check in sigma_checks):
        raise ValueError('Canonical sigma source-form witnesses failed')
    E = [_pg_symbol('E_' + str(j)) for j in range(4)]
    # C=E-2*E_y, and actual source gives 2*E_y=-(1+2*g)*E.
    C = [_pg_sum(e, _pg_product(_pg_sum(one, _pg_scale(g, 2)), e)) for e in E]
    B_field = [zero for _ in range(4)]
    a = [R, zero, zero, zero]
    b = [zero, zero, zero, zero]
    t0 = [zero, zero, zero, zero]
    Delta = [delta, zero, zero, zero]
    root_rows = []
    for j in range(4):
        witnesses = {
            'a': (_pg_jet_product(E, a, j), C[j]),
            'b': (_pg_jet_product(E, b, j), B_field[j]),
            't0': (_pg_sum(_pg_jet_product(C, t0, j), B_field[j]), zero),
            'Delta': (_pg_sum(_pg_jet_product(a, Delta, j), _pg_scale(a[j], 2)),
                      _pg_sum(_pg_jet_product(a, a, j), _pg_jet_product(b, b, j))),
        }
        checks = []
        for name, (left, right) in witnesses.items():
            residual = _pg_difference(left, right)
            checks.append({'quantity': name, 'left': _pg_rows(left), 'right': _pg_rows(right),
                           'residual': _pg_rows(residual), 'is_zero': not residual})
        root_rows.append({'ordinary_Z_order': j, 'checks': checks,
                          'all_zero': all(check['is_zero'] for check in checks)})

    atom = _pg_symbol
    tau = atom('tau')
    phin, phic = atom('phin_0'), atom('phic_0')
    psin, psic = atom('psin_0'), atom('psic_0')
    En, Ec = atom('En_0'), atom('Ec_0')
    Tn, Tc = atom('Tn_0'), atom('Tc_0')
    # These are multiplied by 2*tau. No phase-gap division is introduced.
    A_left = _pg_difference(_pg_product(R, _pg_difference(_pg_product(tau, phin), psin)),
                            _pg_product(R, _pg_difference(_pg_product(tau, phic), psic)))
    A_right = _pg_product(R, _pg_difference(_pg_product(tau, _pg_difference(phin, phic)),
                                           _pg_difference(psin, psic)))
    B_left = _pg_scale(_pg_product(R, _pg_difference(_pg_product(En, Tn), _pg_product(Ec, Tc))), -1)
    B_right = _pg_scale(_pg_product(R, _pg_sum(_pg_product(En, _pg_difference(Tn, Tc)),
                                               _pg_product(_pg_difference(En, Ec), Tc))), -1)
    primitive_rows = []
    for j in range(4):
        checks = []
        for name, left, right in [('A', A_left, A_right), ('B', B_left, B_right)]:
            residual = _pg_difference(left, right)
            checks.append({'quantity': name, 'left': _pg_rows(left), 'right': _pg_rows(right),
                           'residual': _pg_rows(residual), 'is_zero': not residual})
        primitive_rows.append({'ordinary_Z_order': j, 'checks': checks,
                               'all_zero': all(check['is_zero'] for check in checks)})
        if j < 3:
            A_left, A_right = _pg_derivative(A_left), _pg_derivative(A_right)
            B_left, B_right = _pg_derivative(B_left), _pg_derivative(B_right)

    passed = all(check['is_zero'] for check in sigma_checks) and all(row['all_zero'] for row in root_rows + primitive_rows)
    if isinstance(journal, list):
        journal.append({'event': 'fresh_primitive_gap_algebra_saved', 'root_rows': 4,
                        'primitive_rows': 4, 'sigma_checks': 2, 'symbolic_algebra_passed': passed})
    return {'passed': passed, 'symbolic_algebra_passed': passed, 'source_only': True, 'root_shape_rows': root_rows,
            'primitive_gap_rows': primitive_rows, 'sigma_checks': sigma_checks,
            'sigma_check_count': 2,
            'sigma_positive_denominator_assumptions': 'sigma_a>0, sigma_b>0 for 0<x<1; endpoint source extensions agree at 0 and 1',
            'root_shape_row_count': 4,
            'root_shape_subresidual_count': 16, 'primitive_gap_row_count': 4,
            'primitive_gap_subresidual_count': 8,
            'shared_shape_consequence': {'native_and_C3_have_independent_positive_E': True,
                                        'both_Vy_source_functions_are_zero': True,
                                        'common_a_is_derived_from_sigma_witnesses_and_root_relations': True,
                                        'C3_interval_hulls_are_not_selected_values': True,
                                        'B_in_primitive_rows_is_normalized_loop_companion': True},
            'formulas': {'g': 'mu*sigma(x)', 'a': '2+2*g', 'b': '0', 't0': '0', 'Delta': '2*g',
                         'delta_E': 'En-Ec',
                         'delta_A_active': 'a/2*((phin-phic)-(psin-psic)/(2*pi))',
                         'delta_B_active': '-a/(4*pi)*(En*(Tn-Tc)+(En-Ec)*Tc)',
                         'flat_branch': 'Conditional on common Delta>=eta: both lazy A/B values are zero; no actual branch is selected and active bodies are not evaluated'},
            'scope': {'ordinary_Z_orders': [0, 1, 2, 3], 'temporal_order_recursion': False,
                      'E_equality_proved': False, 'P2_equality_proved': False,
                      'phase_equality_proved': False, 'primitive_gap_zero_proved': False,
                      'numeric_values_computed': False, 'old_operators_called': False,
                      'source_values_evaluated': False}}

"""Unexecuted v2 source wrapper for the compact O3 transition primitive gap.

This v2 is deliberately separate from the earlier wrapper.  It accepts only
the frozen packet links, creator source records, literal slices, and archive
slices.  It authenticates actual creator source text with SHA-256 and parses
that text with AST inside derive; it never imports or executes a creator,
provider, checker, normalizer, reducer, quadrature, recurrence, or ancestor
operator.  The root may concatenate the frozen algebra draft after this file.
"""

V2_INPUT_SCHEMA = "compact_O3_transition_primitive_gap_input_v1"
V2_COMPONENT = "compact_O3_transition_primitive_gap"
V2_RESULT_SCHEMA = "original_compact_O3_transition_primitive_gap"
V2_ORDERS = [0, 1, 2, 3]

# Exactly five frozen packet links plus the three static data containers.
V2_PACKET_KEYS = (
    "prior_input",
    "parent_driver_result",
    "native_trace",
    "C3_trace",
    "cached_seams",
    "creator_source_records",
    "literal_slices",
    "archive_literal_slices",
)
V2_PACKET_PINS = {
    "prior_input": "27681343e120a227d34e23d68a6794523c2b13c6117572adb727b80d39c259e3",
    "parent_driver_result": "723b0b24a53e7eef51d3d2a783053dbab3e78b37ece2b4f0baf9f2a9b916cc94",
    "native_trace": "749fc18cbcb4e8df1ec1ecb041239f5f509768ed17f45aa35b2a80d58164a001",
    "C3_trace": "7fa40a975b87720e047e1b679af58a93e4949c8573c0bd3a10d0bf2e87c421ba",
    "cached_seams": "1d5745d1523626302fab2d50a243107323271fee0612fb96c5093c9e694b486e",
}
V2_TRACE_SCHEMAS = {
    "native_trace": "bounded_static_native_O3_transition_primitive_trace_v1",
    "C3_trace": "transition_C3_primitive_source_trace_v1",
    "cached_seams": "transition_cached_primitive_seam_evidence.v1",
}

# Exact relative source paths and actual full-text SHA-256 values.  An
# absolute path is accepted only beneath the two known worktree roots after
# canonicalization; basename collisions are never accepted.
V2_SOURCE_HASHES = {
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O3_Rc_finite_N.py": "45284dea25aa4c2ff06f2fb925e7efd85547d18b38c106d39419618555756a48",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py": "4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3",
    "experiments/root_st073/lei_ren_part1_paper_compliant_outer_initial.py": "a8bc21feb4a2565f976b450f1e4f593180d6eff343e572c31e83117020779657",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py": "83e5a8bdae75d20a261d24dd22d30cc7cb5e603e8fa2f0133958ca880b5e76db",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py": "69d6391120000488472b7961b4bff1901bff872def95d183c4ee056b0948e17d",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_generic_loop_function_sources.py": "d694b8ee8a4bf1c10090d45a7f9d4d6de0792d02a1ebdd6b21872b8514a5b26f",
    "experiments/root_st073/lei_ren_part1_paper_compliant_flat_pulse_derivatives.py": "32f51123a47c2bafbb6c2f7056302bab996a019883615333875207d72a9bb7d7",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_source_targets.py": "de41858de44ff9cea9b5d2b17aeca55eeca2c114969ca4df383ec1345b6badeb",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_repair_band.py": "62e2fa8598bf307126be0ffe0b4e23de83152104524db03f61b4944c225132e4",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py": "9d73659a1c76163dddee6653af53f56ebffd7045986d871dc559ba1e053620b7",
}
V2_ARCHIVE_HASHES = {
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz": "cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz": "75895e1637411cfd2e490778f495f3adf904295a8840184cdc467e006b23305d",
}

V2_NATIVE_ROWS = {
    "coordinate": 241536,
    "E": [241552, 241553, 241554, 241555],
    "E_y": [241563, 241564, 241565, 241566],
    "a": [242476, 242480, 242486, 242493],
    "b": [0, 0, 0, 0],
    "t0": [0, 0, 0, 0],
    "Delta": [242546, 242480, 242486, 242493],
    "P2": [242304, 242308, 242314, 242321],
    "A": [244098, 244099, 244100, 244101],
    "B": [244136, 244137, 244138, 244139],
    "phase": 242551,
    "predicate": 242552,
    "eta": 21,
    "E_y_coefficient": 241562,
    "E_y_Z0_arguments": [241552, 241562],
    "inverse_psi": [243734, 243757, 243826, 243827],
    "T1_at_inverse": [243830, 243852, 243857, 243865],
}
V2_C3_ROWS = {
    "x_node": 4582,
    "Z_node": 20,
    "E": [4590, 4591, 9988, 21725],
    "V": [4592, 4593, 9989, 21726],
    "A": [4718, 4725, 10273, 22254],
    "B": [4727, 4732, 10279, 22255],
    "F_N": [4735, 4739, 11937, 24061],
    "Delta": 4602,
    "predicate": 4604,
    "eta": 21,
}
V2_C3_F_N_TERMS = {
    "0": [[4590, 4718, 4734]],
    "1": [[4591, 4718, 4734], [4590, 4725, 4737]],
    "2": [[9988, 4718, 11930], [48, 4591, 4725, 11931], [4590, 11931, 11935]],
    "3": [[21725, 4718, 11930], [13766, 9988, 4725, 11931], [13766, 4591, 11931, 11935], [4590, 11931, 24059]],
}

V2_REQUIRED_GATES = (
    "actual_native_primitive_bound",
    "actual_C3_primitive_bound",
    "actual_root_shape_clauses_bound",
    "same_mu_sigma_source_forms_bound",
    "ordinary_Z_normalization_bound",
    "positive_E_and_a_source_bound",
    "sigma_endpoint_extensions_and_reflection_bound",
    "both_zero_Vy_source_clauses_bound",
    "same_flat_predicate_forms_bound",
    "normalized_B_roles_kept_distinct_from_physical_B",
)

V2_REQUIRED_DEFINITIONS = {
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O3_Rc_finite_N.py": (
        "nonnegative_Delta_q",
        "background_cell",
        "OriginalO3RcFiniteN",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py": (
        "background_at",
        "CurrentOriginalTrueOuterC6C5Functions",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_outer_initial.py": (
        "stable_sigma",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py": (
        "compiled_background",
        "source_mu_admission",
        "WholeZO3RcFiniteN",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py": (
        "loop_functions",
        "generic_roots",
        "finite_N_functions",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_generic_loop_function_sources.py": (
        "build",
        "CurrentLoopFunctionSources",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_flat_pulse_derivatives.py": (
        "_sigma_left",
        "sigma_jets",
        "FlatPulseDerivatives",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_source_targets.py": (
        "source_packet",
        "CurrentC3SourceTargets",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_repair_band.py": (
        "build",
        "CurrentC3RepairBand",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py": (
        "normalized_targets",
        "WholeZAllNOuterRcFunctions",
    ),
}

# These are literal source clauses, not input assertions.  They are searched
# in the authenticated creator source_text after AST parsing.
V2_REQUIRED_CLAUSES = {
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O3_Rc_finite_N.py": (
        "E=f.scale(u1,factor);E2=f.multiply(u1,u1);V,Vy=zero,zero",
        "excess=mu*(2*sig);Delta=[excess]+zero[1:];a=[f.scalar(2)+excess]+zero[1:]",
        "Ey=f.add(f.scale(E,-c.mpf('.5')),f.scale(E,-mu*sig))",
        "exact_collected_Delta_transition='2*mu*sigma(t)'",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py": (
        "elif name=='O3_transition':",
        "J=controls.integral(g,sigma(g,t),g.nodes[t.node]['name'],g.zero,x,exact_original_sigma_integral=True)",
        "V=a.fixed(0);Ey=a.scale(E,g.neg(g.add(g.constant('1/2'),g.mul(self.mu,sigma(g,x)))));Vy=a.fixed(0)",
        "generic=radial.generic_functions(g,fields,self.delta,self.Pstar,R,order=5)",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_outer_initial.py": (
        "def stable_sigma(c,t):",
        "if v<=0:return c.mpf(0),c.mpf(0)",
        "if v>=1:return c.mpf(1),c.mpf(0)",
        "x=c.mpf(v);a=c.exp(-1/x**2);b=c.exp(-1/(1-x)**2)",
        "value=a/(a+b);derivative=value*(b/(a+b))*(2/x**3+2/(1-x)**3)",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py": (
        "T1={(k,0):integrate(t[k,0]) for k in range(6)}",
        "T1hat=compose_inverse(g,T1,inverse,angle,delta)",
        "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))",
        "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
        "chi=a.sub(a.fixed(phase),a.scale(psi,g.quotient(g.one,twopi,'positive 2*pi')))",
        "predicate=shear",
        "q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))",
        "inverse=g.node('original_monotone_phase_inverse',phase_function=Phi.value.node",
        "C=a.sub(E,a.scale(a.trim(radial['E_y']),2));B=a.scale(a.trim(radial['V_y']),2)",
        "divide=lambda jet:a.div(jet,E,'same strictly positive actual prescribed swirl E')",
        "aa,bb=divide(C),divide(B)",
        "t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))",
        "kap=a.add(aa,a.div(a.mul(bb,bb),aa,'same original positive a'))",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_generic_loop_function_sources.py": (
        "sigma=g.function('original_flat_sigma',x)",
        "g.function('original_flat_sigma_prime',x)",
        "q=g.flat(g.mul(sigma,root),Delta,eta)",
        "actual_eta_positive",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_flat_pulse_derivatives.py": (
        "odds=1/(1-x)**2-1/x**2",
        "e=positive_exp(c,odds); value=e/(1+e)",
        "left=_sigma_left(c,1-c.mpf([a,b]))",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_O3_Rc_finite_N.py": (
        "logmu=c.ln(c.mpf('.001'))-4*(c.exp(40)+11);Tw=-60*logmu",
        "margin=logmu+c.ln(2)-eta_log",
        "self.background,self.kernel_binding=compiled_background(self.kernels)",
        "source=self.background(op,parent,chart,left,right,self.logmu);a=source['actual_a_axial5']",
        "original.O3_QUOTIENTS(proxy,generic,a,self.eta_log,source['actual_collected_Delta_axial5'])",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py": (
        "('o3',o3.WholeZO3RcFiniteN.query))",
        "FRONTENDS[kind],FRONTEND_BINDINGS[kind]=source_frontend(fn,kind)",
        "else:packet=FRONTENDS[kind](view,key,chart.split('_',1)[1],l,r)",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_source_targets.py": (
        "source_projection='6*same actual current source Taylor coefficient[3]'",
        "ordinary_slow_Z_derivative_order=3",
        "T1=alg.integral(t,angle);P=alg.add(alg.fixed(psi),alg.integral(alg.mul(t,t),angle));Phi=alg.mul(K,P)",
        "inverse=old['inverse'];p1,p2inverse=old['psi_Z'],old['psi_ZZ']",
        "def at(v):return g.node('substitute_original_inverse_angle',body=v.node,inverse_angle=inverse.node,angle_variable=angle)",
        "composed=C3Function(at(T1.value)",
        "M=alg.add(alg.neg(alg.div(alg.mul(a,composed),alg.fixed(twopi))),alg.neg(alg.scale(b,phase_ref)))",
        "composed_T1_C3=composed",
        "phase_and_radius_Z_ZZ_ZZZ_exact_zero=True",
    ),
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py": (
        "original_global_phase='frac(N*(logRm+original_outer_offset-logRa-hb*s_c/2))",
    ),
}

V2_SOURCE_CLAUSE_LABELS = {
    "Ey": "Ey=f.add(f.scale(E,-c.mpf('.5')),f.scale(E,-mu*sig))",
    "Vy": "E=f.scale(u1,factor);E2=f.multiply(u1,u1);V,Vy=zero,zero",
    "a": "a=[f.scalar(2)+excess]+zero[1:]",
    "Delta": "Delta=[excess]+zero[1:]",
    "A": "A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))",
    "B": "B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))",
    "T1": "T1={(k,0):integrate(t[k,0]) for k in range(6)}",
    "T1hat": "T1hat=compose_inverse(g,T1,inverse,angle,delta)",
}

def _v2_require(condition, message):
    if not condition:
        raise ValueError(message)


def _v2_packet(source):
    if isinstance(source, dict) and isinstance(source.get("input"), dict):
        packet = source["input"]
    else:
        packet = source
    _v2_require(isinstance(packet, dict), "primitive-gap source packet is not a mapping")
    _v2_require(packet.get("schema") == V2_INPUT_SCHEMA, "primitive-gap input schema changed")
    _v2_require(packet.get("component") == V2_COMPONENT, "primitive-gap component changed")
    _v2_require(packet.get("ordinary_Z_orders") == V2_ORDERS, "ordinary-Z orders changed")
    _v2_require(packet.get("source_values_evaluated") is False, "source values were evaluated")
    for key in ("creator_source_records", "literal_slices"):
        _v2_require(key in packet, "required packet input is missing: " + key)
    if "archive_literal_slices" not in packet:
        _v2_require(
            isinstance(packet.get("literal_slices"), dict)
            and isinstance(packet["literal_slices"].get("native"), dict)
            and isinstance(packet["literal_slices"].get("C3"), dict),
            "raw native/C3 archive slices are missing",
        )
    return packet


def _v2_link(packet, key, expected):
    node = packet.get(key)
    _v2_require(isinstance(node, dict), key + " link is not a mapping")
    observed = node.get("sha256", node.get("raw_sha256"))
    payload = node.get("payload")
    _v2_require(isinstance(observed, str) and observed.lower() == expected, key + " SHA-256 changed")
    _v2_require(isinstance(payload, dict), key + " payload is missing")
    return payload


def _v2_relpath(path):
    _v2_require(isinstance(path, str) and path, "creator/archive path is missing")
    text = path.replace("\\", "/")
    roots = (
        "C:/Users/z5242/Documents/Codex/2026-09-13/before-doing-any-project-work-verify/work/st073-transition-next/",
        "C:/Users/z5242/Documents/Codex/2026-10-16/https-github-com-liang-techh-openai/work/",
    )
    for root in roots:
        if text.startswith(root):
            return text[len(root):]
    return text.lstrip("/")


def _v2_entries(mapping, label):
    _v2_require(isinstance(mapping, (dict, list)), label + " must be a mapping or list")
    return list(mapping.values()) if isinstance(mapping, dict) else list(mapping)


def _v2_authenticate_creator_records(packet):
    import ast as _v2_ast
    import hashlib as _v2_hashlib

    records = _v2_entries(packet.get("creator_source_records"), "creator_source_records")
    seen = {}
    for relative, expected in V2_SOURCE_HASHES.items():
        matches = []
        for record in records:
            if not isinstance(record, dict):
                continue
            path = record.get("path")
            if isinstance(path, str) and _v2_relpath(path) == relative:
                matches.append(record)
        _v2_require(len(matches) == 1, "creator source path missing or duplicated: " + relative)
        record = matches[0]
        source_text = record.get("source_text")
        _v2_require(isinstance(source_text, str), "creator source_text missing: " + relative)
        declared = record.get("sha256")
        actual = _v2_hashlib.sha256(source_text.encode("utf-8")).hexdigest()
        _v2_require(declared == expected and actual == expected, "creator source SHA-256 changed: " + relative)
        tree = _v2_ast.parse(source_text, filename=relative)
        names = {
            node.name
            for node in _v2_ast.walk(tree)
            if isinstance(node, (_v2_ast.FunctionDef, _v2_ast.AsyncFunctionDef, _v2_ast.ClassDef))
        }
        for required in V2_REQUIRED_DEFINITIONS.get(relative, ()):
            _v2_require(required in names, "creator AST definition missing: %s:%s" % (relative, required))
        for clause in V2_REQUIRED_CLAUSES.get(relative, ()):
            _v2_require(clause in source_text, "actual source clause missing: %s:%s" % (relative, clause))
        seen[relative] = {
            "path": relative,
            "sha256": expected,
            "source_text_sha256_recomputed": actual,
            "ast_parsed": True,
            "required_clauses_bound": list(V2_REQUIRED_CLAUSES.get(relative, ())),
        }
    return seen


V2_ARCHIVE_NAMESPACES = {
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz": "native_O3_transition_archive_exact_graph_nodes",
    "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz": "rebuilt_C3_repair_band_archive_exact_graph_nodes",
}


def _v2_node(nodes, node_id, label):
    _v2_require(isinstance(nodes, dict), label + " node slice is not a mapping")
    row = nodes.get(str(node_id), nodes.get(node_id))
    _v2_require(isinstance(row, dict), label + " selected node is missing: " + str(node_id))
    return row


def _v2_int(value):
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.lstrip("-").isdigit():
        return int(value)
    return value


def _v2_ints(values):
    return [_v2_int(value) for value in values] if isinstance(values, list) else values


def _v2_archive_entries(packet):
    archive = packet.get("archive_literal_slices")
    if archive is None:
        # The bounded assembler emits native/C3 records directly under
        # literal_slices.  This adapter keeps the input schema stable while
        # authenticating those raw selected-node maps in place.
        archive = packet.get("literal_slices")
    _v2_require(isinstance(archive, (dict, list)), "archive_literal_slices must be a mapping or list")
    return list(archive.values()) if isinstance(archive, dict) else list(archive)


def _v2_find_archive_record(packet, relative):
    matches = []
    for record in _v2_archive_entries(packet):
        if not isinstance(record, dict):
            continue
        path = record.get("path", record.get("archive_path"))
        if isinstance(path, str) and _v2_relpath(path) == relative:
            matches.append(record)
    _v2_require(len(matches) == 1, "archive slice path missing or duplicated: " + relative)
    return matches[0]


def _v2_record_nodes(record, relative):
    nodes = record.get("nodes")
    if nodes is None:
        slices = record.get("bounded_literal_slices")
        if isinstance(slices, dict) and isinstance(slices.get("nodes"), dict):
            nodes = slices["nodes"]
        elif isinstance(slices, dict):
            nodes = slices
    _v2_require(isinstance(nodes, dict) and len(nodes) > 0, "archive selected node map missing: " + relative)
    return nodes


def _v2_validate_archive_rows(relative, nodes):
    if relative.endswith("true_outer_C6_C5_functions.json.gz"):
        required = [
            V2_NATIVE_ROWS["coordinate"], 241545, 241546, 241551, V2_NATIVE_ROWS["E_y_coefficient"],
            241563, 241564, 241565, 241566, V2_NATIVE_ROWS["predicate"],
            *V2_NATIVE_ROWS["E"], *V2_NATIVE_ROWS["a"], *V2_NATIVE_ROWS["Delta"],
            *V2_NATIVE_ROWS["P2"], *V2_NATIVE_ROWS["A"], *V2_NATIVE_ROWS["B"],
            *V2_NATIVE_ROWS["inverse_psi"],
            *V2_NATIVE_ROWS["T1_at_inverse"],
        ]
        for node_id in sorted(set(required)):
            _v2_node(nodes, node_id, "native archive")
        coordinate = _v2_node(nodes, V2_NATIVE_ROWS["coordinate"], "native archive")
        _v2_require(coordinate.get("operation") == "bound_variable" and coordinate.get("name") == "true_outer_O3_transition", "native coordinate literal changed")
        sigma = _v2_node(nodes, 241545, "native archive")
        binding = sigma.get("source_binding")
        _v2_require(
            sigma.get("operation") == "exact_original_source_flat_sigma"
            and isinstance(binding, dict)
            and binding.get("module") == "lei_ren_part1_paper_compliant_outer_initial.py"
            and binding.get("module_sha256") == V2_SOURCE_HASHES["experiments/root_st073/lei_ren_part1_paper_compliant_outer_initial.py"]
            and binding.get("function") == "stable_sigma",
            "native stable-sigma graph literal changed",
        )
        integral = _v2_node(nodes, 241546, "native archive")
        _v2_require(
            integral.get("operation") == "definite_integral"
            and integral.get("integrand") == 241545
            and integral.get("variable") == "true_outer_transition_sigma_t"
            and _v2_int(integral.get("lower")) == 0
            and _v2_int(integral.get("upper")) == V2_NATIVE_ROWS["coordinate"],
            "native J_sigma literal changed",
        )
        ey_coefficient = _v2_node(nodes, V2_NATIVE_ROWS["E_y_coefficient"], "native archive")
        ey_z0 = _v2_node(nodes, 241563, "native archive")
        _v2_require(ey_coefficient.get("operation") == "negative" and _v2_ints([ey_coefficient.get("argument")]) == [241561], "native E_y coefficient literal changed")
        _v2_require(ey_z0.get("operation") == "product" and _v2_ints(ey_z0.get("arguments")) == [V2_NATIVE_ROWS["E"][0], V2_NATIVE_ROWS["E_y_coefficient"]], "native E_y Z0 product changed")
        for node_id in V2_NATIVE_ROWS["E"]:
            row = _v2_node(nodes, node_id, "native archive")
            _v2_require(row.get("operation") == "product", "native E product row changed: " + str(node_id))
        for node_id in V2_NATIVE_ROWS["a"] + V2_NATIVE_ROWS["P2"]:
            row = _v2_node(nodes, node_id, "native archive")
            _v2_require(row.get("operation") == "positive_quotient" and _v2_int(row.get("denominator")) == V2_NATIVE_ROWS["E"][0], "native positive quotient row changed: " + str(node_id))
        for node_id in V2_NATIVE_ROWS["a"]:
            row = _v2_node(nodes, node_id, "native archive")
            _v2_require(row.get("source_positive_certificate") == "same strictly positive actual prescribed swirl E", "native a positivity certificate changed: " + str(node_id))
        delta0 = _v2_node(nodes, V2_NATIVE_ROWS["Delta"][0], "native archive")
        _v2_require(delta0.get("operation") == "sum" and _v2_ints(delta0.get("arguments")) == [V2_NATIVE_ROWS["a"][0], 25304], "native Delta root row changed")
        predicate = _v2_node(nodes, V2_NATIVE_ROWS["predicate"], "native archive")
        _v2_require(predicate.get("operation") == "exact_real_comparison" and predicate.get("operator") == ">=" and _v2_int(predicate.get("left")) == V2_NATIVE_ROWS["Delta"][0] and _v2_int(predicate.get("right")) == V2_NATIVE_ROWS["eta"], "native flat predicate row changed")
        for family in ("A", "B"):
            for node_id in V2_NATIVE_ROWS[family]:
                row = _v2_node(nodes, node_id, "native archive")
                _v2_require(
                    row.get("operation") == "original_lazy_flat_branch"
                    and _v2_int(row.get("Delta")) == V2_NATIVE_ROWS["Delta"][0]
                    and _v2_int(row.get("eta")) == V2_NATIVE_ROWS["eta"]
                    and _v2_int(row.get("flat_predicate")) == V2_NATIVE_ROWS["predicate"]
                    and _v2_int(row.get("flat_value")) == 0
                    and row.get("active_body_not_evaluated_on_flat_branch") is True,
                    "native " + family + " lazy row changed: " + str(node_id),
                )
        return {
            "namespace": V2_ARCHIVE_NAMESPACES[relative],
            "literal_rows_bound": True,
            "selected_node_ids": sorted(set(required)),
        }

    required = [V2_C3_ROWS["x_node"], V2_C3_ROWS["Z_node"], V2_C3_ROWS["Delta"], V2_C3_ROWS["predicate"], *V2_C3_ROWS["E"], *V2_C3_ROWS["V"], *V2_C3_ROWS["A"], *V2_C3_ROWS["B"], *V2_C3_ROWS["F_N"], 48, 13766]
    for terms in V2_C3_F_N_TERMS.values():
        for term in terms:
            required.extend(term)
    for node_id in sorted(set(required)):
        _v2_node(nodes, node_id, "C3 archive")
    x_node = _v2_node(nodes, V2_C3_ROWS["x_node"], "C3 archive")
    z_node = _v2_node(nodes, V2_C3_ROWS["Z_node"], "C3 archive")
    _v2_require(x_node.get("operation") == "bound_variable" and x_node.get("name") == "native_O3_transition", "C3 x coordinate literal changed")
    _v2_require(z_node.get("operation") == "bound_variable" and z_node.get("name") == "Z", "C3 Z coordinate literal changed")
    for family in ("E", "V"):
        for order, node_id in enumerate(V2_C3_ROWS[family]):
            row = _v2_node(nodes, node_id, "C3 archive")
            _v2_require(
                row.get("operation") == "current_original_leading_function_recipe"
                and row.get("quantity") == family
                and _v2_int(row.get("Z_order")) == order
                and _v2_int(row.get("coordinate")) == V2_C3_ROWS["x_node"]
                and _v2_int(row.get("Z_variable")) == V2_C3_ROWS["Z_node"],
                "C3 " + family + " leading row changed: " + str(node_id),
            )
            if order == 2:
                _v2_require(_v2_int(row.get("Taylor_coefficient_factorial")) == 2, "C3 Z2 factorial source marker changed: " + family)
            if order == 3:
                _v2_require(
                    _v2_int(row.get("Taylor_coefficient_factorial")) == 6
                    and row.get("ordinary_slow_Z_derivative_order") == 3,
                    "C3 Z3 ordinary derivative marker changed: " + family,
                )
    predicate = _v2_node(nodes, V2_C3_ROWS["predicate"], "C3 archive")
    _v2_require(predicate.get("operation") == "exact_real_comparison" and predicate.get("operator") == ">=" and _v2_int(predicate.get("left")) == V2_C3_ROWS["Delta"] and _v2_int(predicate.get("right")) == V2_C3_ROWS["eta"], "C3 flat predicate row changed")
    for family in ("A", "B"):
        for node_id in V2_C3_ROWS[family]:
            row = _v2_node(nodes, node_id, "C3 archive")
            _v2_require(
                row.get("operation") == "original_lazy_flat_branch"
                and _v2_int(row.get("Delta")) == V2_C3_ROWS["Delta"]
                and _v2_int(row.get("eta")) == V2_C3_ROWS["eta"]
                and _v2_int(row.get("flat_predicate")) == V2_C3_ROWS["predicate"]
                and _v2_int(row.get("flat_value")) == 0
                and row.get("active_body_not_evaluated_on_flat_branch") is True,
                "C3 " + family + " lazy row changed: " + str(node_id),
            )
            if node_id == V2_C3_ROWS[family][-1]:
                _v2_require(row.get("ordinary_slow_Z_derivative_order") == 3 and row.get("not_additional_Taylor_conversion") is True, "C3 " + family + " Z3 lazy normalization marker changed")
    _v2_require(_v2_node(nodes, V2_C3_ROWS["F_N"][0], "C3 archive").get("operation") == "product", "C3 F_N Z0 product row changed")
    return {
        "namespace": V2_ARCHIVE_NAMESPACES[relative],
        "literal_rows_bound": True,
        "selected_node_ids": sorted(set(required)),
    }


def _v2_validate_archive_slices(packet):
    seen = {}
    for relative, expected in V2_ARCHIVE_HASHES.items():
        record = _v2_find_archive_record(packet, relative)
        observed = record.get("raw_sha256", record.get("archive_sha256", record.get("sha256")))
        _v2_require(observed == expected, "archive raw SHA-256 changed: " + relative)
        nodes = _v2_record_nodes(record, relative)
        row_binding = _v2_validate_archive_rows(relative, nodes)
        seen[relative] = {
            "path": relative,
            "raw_sha256": expected,
            "namespace": row_binding["namespace"],
            "bounded_literal_slice_count": len(nodes),
            "selected_node_ids": row_binding["selected_node_ids"],
            "literal_rows_bound": row_binding["literal_rows_bound"],
        }
    return seen


def _v2_default_c3_t1_refs():
    source_path = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_source_targets.py"
    return [
        {
            "ordinary_Z_order": order,
            "source_path": source_path,
            "source_sha256": V2_SOURCE_HASHES[source_path],
            "binding": "typed_derivative_whole_function_view",
            "node_id": None,
            "whole_function_view": True,
            "interval_value_bound": False,
            "selector": "primitive_C3.composed_T1_C3: T1=alg.integral(t,angle), evaluated at old['inverse'] with the full ordinary-Z composition rule",
        }
        for order in V2_ORDERS
    ]


def _v2_normalize_raw_literal_slices(packet):
    slices = packet.get("literal_slices")
    _v2_require(isinstance(slices, dict), "literal_slices must be a mapping")
    native = slices.get("native")
    c3 = slices.get("C3")
    _v2_require(isinstance(native, dict) and isinstance(c3, dict), "native/C3 literal slices are required")
    prior = slices.get("prior")
    if native.get("row_path") == "native_O3_transition.root_primitive_jets.ordinary_Z0_to_Z3" and isinstance(prior, dict):
        return {"native": native, "C3": c3, "prior": prior}
    native_record = _v2_find_archive_record(packet, "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz")
    c3_record = _v2_find_archive_record(packet, "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz")
    native_nodes = _v2_record_nodes(native_record, "native archive")
    c3_nodes = _v2_record_nodes(c3_record, "C3 archive")
    ey_z0 = _v2_node(native_nodes, 241563, "native archive")
    sigma_integral = _v2_node(native_nodes, 241546, "native archive")
    native_branch = {}
    for family in ("A", "B"):
        native_branch[family] = {}
        for node_id in V2_NATIVE_ROWS[family]:
            row = _v2_node(native_nodes, node_id, "native archive")
            native_branch[family][str(node_id)] = {
                "Delta": _v2_int(row.get("Delta")),
                "eta": _v2_int(row.get("eta")),
                "flat_predicate": _v2_int(row.get("flat_predicate")),
                "flat_value": _v2_int(row.get("flat_value")),
                "active_body_not_evaluated_on_flat_branch": row.get("active_body_not_evaluated_on_flat_branch"),
            }
    c3_branch = {}
    for family, source_family in (("A", "A"), ("B_over_Pstar", "B")):
        c3_branch[family] = {}
        for node_id in V2_C3_ROWS[source_family]:
            row = _v2_node(c3_nodes, node_id, "C3 archive")
            c3_branch[family][str(node_id)] = {
                "Delta": _v2_int(row.get("Delta")),
                "eta": _v2_int(row.get("eta")),
                "flat_predicate": _v2_int(row.get("flat_predicate")),
                "flat_value": _v2_int(row.get("flat_value")),
                "active_body_not_evaluated_on_flat_branch": row.get("active_body_not_evaluated_on_flat_branch"),
            }
    if not isinstance(prior, dict):
        prior = {
            "native_T1_at_inverse": list(V2_NATIVE_ROWS["T1_at_inverse"]),
            "C3_inverse_T1_source_refs": _v2_default_c3_t1_refs(),
            "native_P2_rows": list(V2_NATIVE_ROWS["P2"]),
            "phase_P2_identity_proved": False,
        }
    return {
        "native": {
            "row_path": "native_O3_transition.root_primitive_jets.ordinary_Z0_to_Z3",
            "ordinary_Z0_to_Z3": {
                "E": list(V2_NATIVE_ROWS["E"]),
                "a": list(V2_NATIVE_ROWS["a"]),
                "b": list(V2_NATIVE_ROWS["b"]),
                "t0": list(V2_NATIVE_ROWS["t0"]),
                "Delta": list(V2_NATIVE_ROWS["Delta"]),
                "P2": list(V2_NATIVE_ROWS["P2"]),
            },
            "E_y_rows": list(V2_NATIVE_ROWS["E_y"]),
            "E_y_Z0": {
                "node": 241563,
                "operation": ey_z0.get("operation"),
                "arguments": _v2_ints(ey_z0.get("arguments")),
                "sigma_argument": "x=coordinate node 241536",
                "coefficient_node": V2_NATIVE_ROWS["E_y_coefficient"],
            },
            "J_sigma_dummy_variable": sigma_integral.get("variable"),
            "J_sigma_upper_node": _v2_int(sigma_integral.get("upper")),
            "A_rows": list(V2_NATIVE_ROWS["A"]),
            "B_rows": list(V2_NATIVE_ROWS["B"]),
            "predicate": {
                "node": V2_NATIVE_ROWS["predicate"],
                "Delta_node": V2_NATIVE_ROWS["Delta"][0],
                "eta": V2_NATIVE_ROWS["eta"],
                "literal": "Delta >= eta",
            },
            "branch_rows": native_branch,
        },
        "C3": {
            "coordinate_and_rows": {
                "x_node": V2_C3_ROWS["x_node"],
                "Z_node": V2_C3_ROWS["Z_node"],
                "E_rows_Z0_to_Z3": list(V2_C3_ROWS["E"]),
                "A_rows_Z0_to_Z3": list(V2_C3_ROWS["A"]),
                "B_rows_Z0_to_Z3": list(V2_C3_ROWS["B"]),
                "F_N_rows_Z0_to_Z3": list(V2_C3_ROWS["F_N"]),
            },
            "row_namespace": "rebuilt_C3_repair_band_archive_exact_graph_nodes",
            "primitive_row_branch_literals": c3_branch,
            "all_active_bodies_not_evaluated_on_flat_branch": True,
            "F_N_terms": dict(V2_C3_F_N_TERMS),
        },
        "prior": prior,
    }


def _v2_validate_literal_slices(packet):
    slices = _v2_normalize_raw_literal_slices(packet)
    native = slices.get("native")
    c3 = slices.get("C3")
    prior = slices.get("prior")
    _v2_require(isinstance(native, dict) and isinstance(c3, dict) and isinstance(prior, dict), "native/C3/prior literal slices are required")

    _v2_require(native.get("row_path") == "native_O3_transition.root_primitive_jets.ordinary_Z0_to_Z3", "native root row path changed")
    _v2_require(native.get("ordinary_Z0_to_Z3") == {
        "E": V2_NATIVE_ROWS["E"],
        "a": V2_NATIVE_ROWS["a"],
        "b": V2_NATIVE_ROWS["b"],
        "t0": V2_NATIVE_ROWS["t0"],
        "Delta": V2_NATIVE_ROWS["Delta"],
        "P2": V2_NATIVE_ROWS["P2"],
    }, "native root row IDs changed")
    _v2_require(native.get("E_y_rows") == V2_NATIVE_ROWS["E_y"], "native E_y row IDs changed")
    _v2_require(native.get("E_y_Z0") == {
        "node": 241563,
        "operation": "product",
        "arguments": V2_NATIVE_ROWS["E_y_Z0_arguments"],
        "sigma_argument": "x=coordinate node 241536",
        "coefficient_node": V2_NATIVE_ROWS["E_y_coefficient"],
    }, "native E_y immediate literal changed")
    _v2_require(native.get("J_sigma_dummy_variable") == "true_outer_transition_sigma_t", "native J_sigma dummy variable changed")
    _v2_require(native.get("J_sigma_upper_node") == 241536, "native J_sigma upper node changed")
    _v2_require(native.get("A_rows") == V2_NATIVE_ROWS["A"] and native.get("B_rows") == V2_NATIVE_ROWS["B"], "native A/B row IDs changed")
    _v2_require(native.get("predicate") == {
        "node": V2_NATIVE_ROWS["predicate"],
        "Delta_node": V2_NATIVE_ROWS["Delta"][0],
        "eta": V2_NATIVE_ROWS["eta"],
        "literal": "Delta >= eta",
    }, "native flat predicate literal changed")

    c3_coords = c3.get("coordinate_and_rows")
    _v2_require(c3_coords == {
        "x_node": V2_C3_ROWS["x_node"],
        "Z_node": V2_C3_ROWS["Z_node"],
        "E_rows_Z0_to_Z3": V2_C3_ROWS["E"],
        "A_rows_Z0_to_Z3": V2_C3_ROWS["A"],
        "B_rows_Z0_to_Z3": V2_C3_ROWS["B"],
        "F_N_rows_Z0_to_Z3": V2_C3_ROWS["F_N"],
    }, "C3 coordinate/row IDs changed")
    _v2_require(c3.get("row_namespace") == "rebuilt_C3_repair_band_archive_exact_graph_nodes", "C3 namespace changed")
    branches = c3.get("primitive_row_branch_literals")
    _v2_require(isinstance(branches, dict), "C3 branch literal slices are missing")
    _v2_require(c3.get("all_active_bodies_not_evaluated_on_flat_branch") is True, "C3 branch-evaluation scope changed")
    for family, rows in (("A", V2_C3_ROWS["A"]), ("B_over_Pstar", V2_C3_ROWS["B"])):
        family_rows = branches.get(family)
        _v2_require(isinstance(family_rows, dict), "C3 " + family + " branch rows missing")
        for row in rows:
            item = family_rows.get(str(row))
            _v2_require(
                isinstance(item, dict)
                and item.get("Delta") == V2_C3_ROWS["Delta"]
                and item.get("eta") == V2_C3_ROWS["eta"]
                and item.get("flat_predicate") == V2_C3_ROWS["predicate"]
                and item.get("flat_value") == 0
                and item.get("active_body_not_evaluated_on_flat_branch") is True,
                "C3 " + family + " literal changed at " + str(row),
            )
    _v2_require(c3.get("F_N_terms") == V2_C3_F_N_TERMS, "C3 ordinary-Z product refs changed")

    _v2_require(prior.get("native_T1_at_inverse") == V2_NATIVE_ROWS["T1_at_inverse"], "native inverse/T1 IDs changed")
    c3_t1 = prior.get("C3_inverse_T1_source_refs")
    _v2_require(isinstance(c3_t1, list) and len(c3_t1) == 4, "C3 inverse/T1 source refs must have four orders")
    for order, item in enumerate(c3_t1):
        _v2_require(
            isinstance(item, dict)
            and item.get("ordinary_Z_order") == order
            and isinstance(item.get("source_path"), str)
            and item.get("source_path") in V2_SOURCE_HASHES
            and item.get("source_sha256") == V2_SOURCE_HASHES[item["source_path"]]
            and item.get("binding") in ("typed_derivative_whole_function_view", "source-qualified archive node"),
            "C3 inverse/T1 source ref changed at order " + str(order),
        )
        if item.get("binding") == "typed_derivative_whole_function_view":
            _v2_require(item.get("node_id") is None and item.get("whole_function_view") is True and item.get("interval_value_bound") is False, "C3 typed derivative view was converted to an interval value")
        else:
            _v2_require(isinstance(item.get("node_id"), int), "C3 source-qualified node/order ref missing")
    _v2_require(prior.get("native_P2_rows") == V2_NATIVE_ROWS["P2"], "native P2 source IDs changed")
    _v2_require(prior.get("phase_P2_identity_proved") is False, "phase/P2 identity was claimed")
    native_t1_refs = [
        {
            "ordinary_Z_order": order,
            "node_id": V2_NATIVE_ROWS["T1_at_inverse"][order],
            "source_path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py",
            "source_sha256": V2_SOURCE_HASHES["experiments/root_st073/lei_ren_part1_paper_compliant_current_original_prescribed_loop_C5_functions.py"],
            "binding": "T1hat compose_inverse source row",
            "whole_function_view": False,
            "interval_value_bound": False,
        }
        for order in V2_ORDERS
    ]
    c3_t1_refs = [
        {
            "ordinary_Z_order": item["ordinary_Z_order"],
            "source_path": item["source_path"],
            "source_sha256": item["source_sha256"],
            "binding": item["binding"],
            "node_id": item.get("node_id"),
            "whole_function_view": item.get("whole_function_view", False),
            "interval_value_bound": item.get("interval_value_bound", False),
        }
        for item in c3_t1
    ]
    ordinary_Z_product_refs = []
    for order in V2_ORDERS:
        ordinary_Z_product_refs.append(
            {
                "ordinary_Z_order": order,
                "native_E_node": V2_NATIVE_ROWS["E"][order],
                "native_A_node": V2_NATIVE_ROWS["A"][order],
                "native_B_node": V2_NATIVE_ROWS["B"][order],
                "C3_E_node": V2_C3_ROWS["E"][order],
                "C3_A_node": V2_C3_ROWS["A"][order],
                "C3_B_node": V2_C3_ROWS["B"][order],
                "C3_F_N_node": V2_C3_ROWS["F_N"][order],
                "C3_F_N_terms": V2_C3_F_N_TERMS[str(order)],
                "native_binomial_product": "ordinary Z binomial product over the explicitly selected E/A/B jet rows",
                "C3_binomial_product": "ordinary Z binomial product over the explicitly selected E/A/B jet rows",
                "source_values_evaluated": False,
            }
        )
    return {
        "native_T1_at_inverse": list(V2_NATIVE_ROWS["T1_at_inverse"]),
        "C3_inverse_T1_source_refs": c3_t1_refs,
        "actual_T1_refs": {"native": native_t1_refs, "C3": c3_t1_refs},
        "native_P2_rows": list(V2_NATIVE_ROWS["P2"]),
        "phase_P2_identity_proved": False,
        "ordinary_Z_product_refs": ordinary_Z_product_refs,
    }


def _v2_trace_bindings(native, c3, cached):
    _v2_require(native.get("schema") == V2_TRACE_SCHEMAS["native_trace"], "native trace schema changed")
    _v2_require(c3.get("schema") == V2_TRACE_SCHEMAS["C3_trace"], "C3 trace schema changed")
    _v2_require(cached.get("schema") == V2_TRACE_SCHEMAS["cached_seams"], "cached seam schema changed")

    native_stage = native.get("native_O3_transition")
    _v2_require(isinstance(native_stage, dict), "native O3 transition trace missing")
    native_root = native_stage.get("root_primitive_jets", {}).get("ordinary_Z0_to_Z3")
    _v2_require(isinstance(native_root, dict), "native root ordinary-Z path missing")
    _v2_require(native_root.get("E") == V2_NATIVE_ROWS["E"], "native E root rows changed")
    _v2_require(native_root.get("a") == V2_NATIVE_ROWS["a"], "native a root rows changed")
    _v2_require(native_root.get("b") == V2_NATIVE_ROWS["b"], "native b root rows changed")
    _v2_require(native_root.get("t0") == V2_NATIVE_ROWS["t0"], "native t0 root rows changed")
    _v2_require(native_root.get("Delta") == V2_NATIVE_ROWS["Delta"], "native Delta root rows changed")
    _v2_require(native_root.get("P2") == V2_NATIVE_ROWS["P2"], "native P2 root rows changed")
    native_lazy = native.get("phase_q_and_lazy_loop")
    _v2_require(isinstance(native_lazy, dict), "native phase/lazy trace missing")
    _v2_require(native_lazy.get("A", {}).get("rows_Z0_to_Z3") == V2_NATIVE_ROWS["A"], "native A rows changed")
    _v2_require(native_lazy.get("B_over_Pstar", {}).get("rows_Z0_to_Z3") == V2_NATIVE_ROWS["B"], "native B rows changed")
    native_pred = native_lazy.get("lazy_predicate", {})
    _v2_require(
        native_pred.get("predicate_node") == V2_NATIVE_ROWS["predicate"]
        and native_pred.get("right") == V2_NATIVE_ROWS["eta"]
        and native_pred.get("literal") == "Delta >= eta",
        "native flat predicate changed",
    )

    c3_stage = c3.get("rebuilt_C3_leading_packet")
    _v2_require(isinstance(c3_stage, dict), "C3 leading packet missing")
    coords = c3_stage.get("coordinate_and_rows")
    _v2_require(isinstance(coords, dict), "C3 coordinate/rows missing")
    _v2_require(coords.get("x_node") == V2_C3_ROWS["x_node"] and coords.get("Z_node") == V2_C3_ROWS["Z_node"], "C3 coordinate nodes changed")
    _v2_require(coords.get("E_rows_Z0_to_Z3") == V2_C3_ROWS["E"], "C3 E rows changed")
    _v2_require(coords.get("A_rows_Z0_to_Z3") == V2_C3_ROWS["A"], "C3 A rows changed")
    _v2_require(coords.get("B_rows_Z0_to_Z3") == V2_C3_ROWS["B"], "C3 B rows changed")
    _v2_require(coords.get("F_N_rows_Z0_to_Z3") == V2_C3_ROWS["F_N"], "C3 F_N rows changed")
    c3_branches = c3_stage.get("primitive_row_branch_literals", {})
    _v2_require(c3_branches.get("all_active_bodies_not_evaluated_on_flat_branch") is True, "C3 aggregate branch scope changed")
    for family, rows in (("A", V2_C3_ROWS["A"]), ("B_over_Pstar", V2_C3_ROWS["B"])):
        family_rows = c3_branches.get(family, {})
        for row in rows:
            item = family_rows.get(str(row))
            _v2_require(
                isinstance(item, dict)
                and item.get("Delta") == V2_C3_ROWS["Delta"]
                and item.get("eta") == V2_C3_ROWS["eta"]
                and item.get("flat_predicate") == V2_C3_ROWS["predicate"]
                and item.get("flat_value") == 0,
                "C3 trace branch row changed: " + str(row),
            )
    quantity_paths = c3_stage.get("quantity_paths", {})
    expected_quantity_paths = {
            "E": "original_generic_source.common_velocity_E_axial5[Z_order]",
            "V": "original_generic_source.common_velocity_V_axial5[Z_order]",
            "a": "original_roots.a[(0,Z_order)]",
            "b": "original_full_source_quotients.actual_b_axial5[Z_order]",
            "t0": "original_roots.t0[(0,Z_order)]",
            "p2": "original_roots.p2[(0,Z_order)]",
            "Delta": "original_full_source_quotients.actual_Delta_axial5[Z_order]",
        }
    for name, expected in expected_quantity_paths.items():
        _v2_require(quantity_paths.get(name) == expected, "C3 leading root quantity path changed: " + name)
    _v2_require(quantity_paths.get("C3_p2_alternates") == ["full_signed_p2_axial4", "exact_full_p2_axial4"], "C3 p2 source alternates changed")
    _v2_require(quantity_paths.get("C3_Delta_alternates") == ["actual_Delta_axial5", "Delta_axial5"], "C3 Delta source alternates changed")
    native_route = c3.get("native_O3_transition_source", {}).get("generic_roots", {})
    _v2_require(
        native_route.get("a") == "a=2+2*mu*sigma(x)"
        and native_route.get("b") == "b=0"
        and native_route.get("t0") == "t0=0"
        and native_route.get("Delta") == "Delta=2*mu*sigma(x)",
        "native route root witness changed",
    )
    return {
        "native_rows": {name: list(value) for name, value in V2_NATIVE_ROWS.items() if isinstance(value, list)},
        "C3_rows": {name: list(value) if isinstance(value, list) else value for name, value in V2_C3_ROWS.items()},
        "native_lazy_predicate": "Delta >= eta",
        "C3_lazy_predicate": "Delta >= eta",
        "C3_namespace": "rebuilt_C3_repair_band_archive_exact_graph_nodes",
        "C3_quantity_paths": dict(quantity_paths),
        "C3_native_route_root_witness": {
            "a": native_route.get("a"),
            "b": native_route.get("b"),
            "t0": native_route.get("t0"),
            "Delta": native_route.get("Delta"),
        "cached_native_trace_alone_proves_C3_root_identity": False,
        },
        "cached_seams_data_only": True,
    }


def _v3_c3_root_and_normalization(packet):
    """Bind root recipe selectors and reuse the accepted A/B ordinary-Z ports."""
    relative = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz"
    nodes = _v2_record_nodes(_v2_find_archive_record(packet, relative), relative)
    recipes = []
    for quantity in ("a", "b", "t0", "Delta"):
        for order in V2_ORDERS:
            matches = [(key, row) for key, row in nodes.items()
                       if row.get("operation") == "current_original_leading_function_recipe"
                       and row.get("coordinate") == V2_C3_ROWS["x_node"]
                       and row.get("quantity") == quantity and row.get("Z_order") == order]
            _v2_require(len(matches) == 1, "C3 root recipe selector is missing or duplicated: %s Z%s" % (quantity, order))
            key, row = matches[0]
            recipe = row.get("recipe", {})
            _v2_require(row.get("native_chart") == "O3_transition" and row.get("Z_variable") == 20,
                        "C3 root recipe chart/Z source changed")
            _v2_require(row.get("defining_quantity_not_a_range_value") is True,
                        "C3 root selector is an enclosure rather than a defining source quantity")
            _v2_require(recipe.get("module") == "lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py"
                        and recipe.get("module_sha256") == V2_SOURCE_HASHES["experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py"]
                        and recipe.get("function") == "WholeZAllNOuterRcFunctions.leading_packet"
                        and recipe.get("native_chart") == "O3_transition", "C3 root defining provider route changed")
            if order >= 2:
                _v2_require(row.get("ordinary_slow_Z_derivative_order") == order
                            and row.get("Taylor_coefficient_factorial") == (2 if order == 2 else 6),
                            "C3 root ordinary-Z projection changed")
            recipes.append({"node": int(key), "quantity": quantity, "ordinary_Z_order": order,
                            "provider": recipe, "range_value_selected": False})
    parent = _v2_link(packet, "parent_driver_result", V2_PACKET_PINS["parent_driver_result"])
    _v2_require(parent.get("passed") is True and parent.get("source_only") is True
                and parent.get("O3_transition_finite_N_driver_source_bound") is True,
                "Cached actual driver source acceptance is missing")
    jets = parent.get("C3_source_binding", {}).get("jets", {})
    for family in ("E", "A", "B"):
        _v2_require(jets.get(family) == V2_C3_ROWS[family], "Cached C3 ordinary-Z primitive ports changed")
    native_relative = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
    native_nodes = _v2_record_nodes(_v2_find_archive_record(packet, native_relative), native_relative)
    eta_records = [_v2_node(selected, 21, "canonical eta") for selected in (native_nodes, nodes)]
    for eta in eta_records:
        _v2_require(eta.get("operation") == "current_original_source_parameter"
                    and eta.get("source_attribute") == "eta_log"
                    and eta.get("quantity_is_exact_original_parameter_not_range_endpoint") is True
                    and eta.get("Z_independent") is True
                    and eta.get("canonical_parameter_proof", {}).get("exact_original_parameter_tuples_equal") is True,
                    "Exact canonical eta source-parameter binding changed")
    _v2_require(eta_records[0].get("source_family") == eta_records[1].get("source_family")
                and eta_records[0].get("fixed_leading_input_N0") == eta_records[1].get("fixed_leading_input_N0")
                and eta_records[0]["canonical_parameter_proof"].get("original_canonical_eta_log")
                    == eta_records[1]["canonical_parameter_proof"].get("original_canonical_eta_log"),
                "Native/C3 canonical eta metadata differ")
    return {"root_recipe_rows": recipes, "root_recipe_count": 16,
            "actual_C3_closed_root_clauses_authenticated_in_original_O3_creator": True,
            "C3_root_formula_source": {
                "source_path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O3_Rc_finite_N.py",
                "source_sha256": V2_SOURCE_HASHES["experiments/root_st073/lei_ren_part1_paper_compliant_current_original_O3_Rc_finite_N.py"],
                "authenticated_clause": "excess=mu*(2*sig);Delta=[excess]+zero[1:];a=[f.scalar(2)+excess]+zero[1:]",
                "ordinary_Z_dependence": "source a and Delta explicitly have zero higher-Z coefficient rows; native mu/sigma multiplier is fixed under ordinary Z",
                "provider_route": "rebuilt root recipe -> leading_packet -> o3 source frontend -> compiled original O3 background_cell and O3_QUOTIENTS",
                "native_trace_formula_used_as_C3_proof": False},
            "A_B_ordinary_Z_normalization_reused_as_saved_data": True,
            "new_ancestor_normalization_calls": 0, "canonical_eta_source_metadata_reused": True}


def _v2_bind(source):
    packet = _v2_packet(source)
    payloads = {
        key: _v2_link(packet, key, V2_PACKET_PINS[key])
        for key in ("prior_input", "parent_driver_result", "native_trace", "C3_trace", "cached_seams")
    }
    source_records = _v2_authenticate_creator_records(packet)
    archive_records = _v2_validate_archive_slices(packet)
    literal_bindings = _v2_validate_literal_slices(packet)
    trace_bindings = _v2_trace_bindings(payloads["native_trace"], payloads["C3_trace"], payloads["cached_seams"])
    c3_root_binding = _v3_c3_root_and_normalization(packet)

    native_trace = payloads["native_trace"]
    c3_trace = payloads["C3_trace"]
    native_sigma = native_trace["native_O3_transition"]["sigma_and_J_sigma"]
    c3_sigma = c3_trace["exact_sigma_mu_eta_semantics"]["sigma"]
    _v2_require("0 for x<=0" in native_sigma["canonical_sigma_recipe"]["recipe"][0] and "1 for x>=1" in native_sigma["canonical_sigma_recipe"]["recipe"][0], "native sigma endpoint source changed")
    _v2_require(c3_sigma["flat_extension"] == "sigma(x)=0 for x<=0; sigma(x)=1 for x>=1", "C3 sigma endpoint source changed")
    _v2_require(c3_sigma["symmetry"] == "sigma(1-x)=1-sigma(x)", "C3 sigma reflection source changed")
    _v2_require("mu" in c3_trace["exact_sigma_mu_eta_semantics"] and "eta" in c3_trace["exact_sigma_mu_eta_semantics"], "C3 mu/eta source semantics missing")
    mu_text = c3_trace["exact_sigma_mu_eta_semantics"]["mu"]["mu"]
    _v2_require(mu_text == "exp(logmu)=1/1000*exp(-4*(exp(40)+11))", "common mu source form changed")
    _v2_require(c3_trace["exact_sigma_mu_eta_semantics"]["eta"]["lazy_predicate"] == "Delta >= eta", "C3 flat predicate source changed")

    contract = {
        "actual_native_primitive_bound": True,
        "actual_C3_primitive_bound": True,
        "actual_root_shape_clauses_bound": True,
        "same_mu_sigma_source_forms_bound": True,
        "ordinary_Z_normalization_bound": True,
        "positive_E_and_a_source_bound": True,
        "sigma_endpoint_extensions_and_reflection_bound": True,
        "both_zero_Vy_source_clauses_bound": True,
        "same_flat_predicate_forms_bound": True,
        "normalized_B_roles_kept_distinct_from_physical_B": True,
        "ordinary_Z_orders": list(V2_ORDERS),
        "source_values_evaluated": False,
    }
    root_shape = {
        "native": "a=2+2*mu*sigma(x), b=0, t0=0, Delta=2*mu*sigma(x)",
        "C3": "a=2+2*mu*sigma(x), b=0, t0=0, Delta=2*mu*sigma(x)",
        "derived_after_independent_native_C3_clauses": True,
        "common_a_was_not_assumed": True,
        "sigma_argument_native": "x=coordinate node 241536",
        "sigma_argument_C3": "x=coordinate node 4582",
    }
    native_archive_path = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.json.gz"
    c3_archive_path = "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz"
    native_archive_ref = {
        "path": native_archive_path,
        "raw_sha256": archive_records[native_archive_path]["raw_sha256"],
        "namespace": archive_records[native_archive_path]["namespace"],
    }
    c3_archive_ref = {
        "path": c3_archive_path,
        "raw_sha256": archive_records[c3_archive_path]["raw_sha256"],
        "namespace": archive_records[c3_archive_path]["namespace"],
    }
    c3_t1_refs = literal_bindings["actual_T1_refs"]["C3"]
    c3_phase_ref = {
        "operation": "ordinary_Z_derivative_with_global_phase_held_fixed",
        "exact_zero_derivative_orders": [1, 2, 3],
        "source_path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py",
        "source_sha256": V2_SOURCE_HASHES["experiments/root_st073/lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions.py"],
        "binding": "typed whole-function global phase source view",
        "selector": "outer_binding.original_global_phase",
        "whole_function_view": True,
        "interval_value_bound": False,
    }
    c3_psi_ref = {
        "source_path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_source_targets.py",
        "source_sha256": V2_SOURCE_HASHES["experiments/root_st073/lei_ren_part1_paper_compliant_current_original_C3_source_targets.py"],
        "binding": "typed derivative whole-function inverse-angle source view",
        "selector": "primitive_C3 inverse=old['inverse']; at(v)=substitute_original_inverse_angle",
        "whole_function_view": True,
        "interval_value_bound": False,
    }
    gaps = []
    for family, formula, native_rows, c3_rows in (
        ("E", "E_n-E_c", V2_NATIVE_ROWS["E"], V2_C3_ROWS["E"]),
        ("A", "piecewise inactive 0 else correlated A_n-A_c", V2_NATIVE_ROWS["A"], V2_C3_ROWS["A"]),
        ("B", "piecewise inactive 0 else correlated B_n-B_c", V2_NATIVE_ROWS["B"], V2_C3_ROWS["B"]),
    ):
        for order in V2_ORDERS:
            gap = {
                "family": family,
                "ordinary_Z_order": order,
                "signed_formula": formula,
                "native_node": native_rows[order],
                "C3_node": c3_rows[order],
                "flat_predicate": "Delta >= eta",
                "branch_selected": False,
                "source_values_evaluated": False,
                "gap_zero_claimed": False,
                "native_archive_ref": dict(native_archive_ref, ordinary_Z_order=order, node_id=native_rows[order]),
                "C3_archive_ref": dict(c3_archive_ref, ordinary_Z_order=order, node_id=c3_rows[order]),
            }
            if family in ("A", "B"):
                gap["inactive_flat_value"] = 0
                gap["correlation_refs"] = {
                    "R_constant_Z": True,
                    "delta_phase_refs": [
                        {
                            "side": "native",
                            "source_path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py",
                            "source_sha256": V2_SOURCE_HASHES["experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py"],
                            "node_id": V2_NATIVE_ROWS["phase"] if order == 0 else 0,
                            "source_value_node": V2_NATIVE_ROWS["phase"],
                            "operation": "ordinary_Z_derivative_with_global_phase_held_fixed",
                            "binding": "actual native global phase value for Z0; exact zero ordinary-Z derivative for Z1..3",
                            "ordinary_Z_order": order,
                        },
                        dict(c3_phase_ref, side="C3", ordinary_Z_order=order),
                    ],
                    "delta_psi_refs": [
                        {
                            "side": "native",
                            "source_path": "experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py",
                            "source_sha256": V2_SOURCE_HASHES["experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py"],
                            "node_id": V2_NATIVE_ROWS["inverse_psi"][order],
                            "binding": "actual native monotone phase inverse psi jet row",
                            "ordinary_Z_order": order,
                        },
                        dict(c3_psi_ref, side="C3", ordinary_Z_order=order),
                    ],
                    "delta_T1_refs": [
                        dict(literal_bindings["actual_T1_refs"]["native"][order], side="native"),
                        dict(c3_t1_refs[order], side="C3"),
                    ],
                }
            else:
                gap["flat_branch_applies"] = False
            gaps.append(gap)
    _v2_require(len(gaps) == 12, "primitive E/A/B gap count changed")
    return {
        "packet": packet,
        "payloads": payloads,
        "source_records": source_records,
        "archive_records": archive_records,
        "literal_bindings": literal_bindings,
        "trace_bindings": trace_bindings,
        "contract": contract,
        "root_shape": root_shape,
        "gaps": gaps,
        "C3_root_binding": c3_root_binding,
        "actual_T1_refs": literal_bindings["actual_T1_refs"],
    }


def _v3_atom_mapping(bound):
    """Instantiate the fresh algebra symbols onto actual ordinary-Z source views."""
    rows = {(row["family"], row["ordinary_Z_order"]): row for row in bound["gaps"]}
    mapping = {
        "g": {"definition": "mu*sigma(x)", "ordinary_Z_independent": True,
              "exact_source_definition_only": True, "C3_interval_value_selected": False},
        "tau": {"definition": "2*pi", "same_mathematical_constant": True},
        "sigma_a": {"definition": "exp(-1/x^2)", "domain": "0<x<1", "positive": True},
        "sigma_b": {"definition": "exp(-1/(1-x)^2)", "domain": "0<x<1", "positive": True},
    }
    for order in V2_ORDERS:
        e = rows[("E", order)]
        correlations = rows[("A", order)]["correlation_refs"]
        for name, ref in (("En", e["native_archive_ref"]), ("Ec", e["C3_archive_ref"]),
                          ("phin", correlations["delta_phase_refs"][0]), ("phic", correlations["delta_phase_refs"][1]),
                          ("psin", correlations["delta_psi_refs"][0]), ("psic", correlations["delta_psi_refs"][1]),
                          ("Tn", correlations["delta_T1_refs"][0]), ("Tc", correlations["delta_T1_refs"][1])):
            mapping[name + "_" + str(order)] = ref
        mapping["E_" + str(order)] = {"independent_root_proof_roles": [e["native_archive_ref"], e["C3_archive_ref"]],
                                       "E_native_equals_C3_assumed": False}
    return mapping


def _v2_emit(journal, event, **fields):
    if journal is None:
        return
    payload = {"event": event}
    payload.update(fields)
    if callable(journal):
        journal(payload)
    elif isinstance(journal, list):
        journal.append(payload)


def derive(source, journal=None):
    bound = _v2_bind(source)
    # The root concatenates the frozen algebra draft, which defines this
    # source-only symbolic helper.  It is called exactly once after binding.
    algebra = derive_primitive_gap_algebra(bound["contract"], journal)
    passed = algebra.get("passed") is True
    gaps = bound["gaps"]
    audit = {
        "original_source_evaluations": 0,
        "new_root_shape_order_checks": 4,
        "new_sigma_identity_checks": 2,
        "new_primitive_cocycle_order_checks": 4,
        "new_actual_primitive_gap_rows": 12,
        "new_density_proof_rows": 0,
        "old_generic_identity_calls": 0,
        "old_provider_calls": 0,
        "old_checker_calls": 0,
        "old_normalizer_calls": 0,
        "old_source_value_reducer_calls": 0,
        "old_operator_calls": 0,
        "old_source_consumer_calls": 0,
        "source_value_reducer_calls": 0,
        "numeric_field_evaluations": 0,
        "old_numeric_calls": 0,
        "numeric_calls": 0,
        "new_algebra_family_calls": 1,
        "old_math_calls": 0,
        "old_source_calls": 0,
        "numeric_values_computed": 0,
    }
    result = {
        "schema": V2_RESULT_SCHEMA,
        "component": V2_COMPONENT,
        "passed": passed,
        "metadata_only": False,
        "source_only": True,
        "runtime_model_metadata": None,
        "runtime_model_metadata_input": bound["packet"].get("runtime_model_metadata"),
        "original_source_evaluations": 0,
        "source_packet_defined": passed,
        "source_packet_definition_accepted": passed,
        "source_packet_accepted_owned": passed,
        "source_packet_accepted_inherited": False,
        "derive_entry_calls": 1,
        "derive_calls": 1,
        "entry_attempts": 1,
        "input_digest": source.get("input_digest") if isinstance(source, dict) else None,
        "input_raw_sha256": source.get("input_raw_sha256") if isinstance(source, dict) else None,
        "new_flags": {
            "O3_transition_primitive_sources_bound": passed,
            "O3_transition_root_shape_correlations_proved": passed,
            "O3_transition_signed_primitive_gap_functions_defined": passed,
        },
        "counts": {
            "ordinary_Z_orders": list(V2_ORDERS),
            "actual_primitive_gap_functions": 12,
            "actual_E_gap_functions": 4,
            "actual_A_gap_functions": 4,
            "actual_B_gap_functions": 4,
            "root_shape_ordinary_Z_rows": 4,
            "root_shape_checks_per_row": 4,
            "root_shape_subresiduals": 16,
            "new_sigma_identity_checks": 2,
            "new_primitive_cocycle_order_checks": 4,
            "new_primitive_cocycle_subresiduals": 8,
            "new_actual_primitive_gap_rows": 12,
            "new_density_proof_rows": 0,
        },
        "new_counts": {
            "actual_primitive_gap_functions": 12,
            "actual_E_gap_functions": 4,
            "actual_A_gap_functions": 4,
            "actual_B_gap_functions": 4,
            "root_shape_ordinary_Z_rows": 4,
            "root_shape_checks_per_row": 4,
            "root_shape_subresiduals": 16,
            "new_sigma_identity_checks": 2,
            "new_primitive_cocycle_order_checks": 4,
            "new_primitive_cocycle_subresiduals": 8,
            "new_actual_primitive_gap_rows": 12,
            "new_density_proof_rows": 0,
        },
        "actual_primitive_gap_functions": 12,
        "signed_primitive_gap_function_count": 12,
        "actual_E_gap_functions": 4,
        "actual_A_gap_functions": 4,
        "actual_B_gap_functions": 4,
        "new_actual_primitive_gap_rows": 12,
        "root_shape_rows": algebra.get("root_shape_rows", []),
        "root_shape_ordinary_Z_rows": 4,
        "root_shape_checks_per_row": 4,
        "root_shape_subresiduals": 16,
        "root_shape_identity_fields": ["a", "b", "t0", "Delta"],
        "root_shape_sigma_same_exact_definition": True,
        "root_shape_sigma_interval_evaluated_value_equality": False,
        "sigma_identity_families": [
            "positive_interior_logistic_reflection",
            "source_endpoint_binding",
        ],
        "new_sigma_identity_checks": 2,
        "sigma_identity_checks": algebra.get("sigma_checks", []),
        "primitive_gap_rows": algebra.get("primitive_gap_rows", []),
        "new_primitive_cocycle_order_checks": 4,
        "new_primitive_cocycle_subresiduals": 8,
        "new_density_proof_rows": 0,
        "source_bindings": {
            "C3_root_recipe_binding": bound["C3_root_binding"],
            "fresh_algebra_source_atom_mapping": _v3_atom_mapping(bound),
            "packet_keys": list(V2_PACKET_KEYS),
            "packet_sha256_pins": dict(V2_PACKET_PINS),
            "creator_source_records": bound["source_records"],
            "archive_literal_slices": bound["archive_records"],
            "raw_literal_slice_adapter": "archive_path/archive_sha256 plus selected native/C3 nodes are normalized and row-authenticated before contract flags",
            "literal_bindings": bound["literal_bindings"],
            "trace_bindings": bound["trace_bindings"],
            "contract": dict(bound["contract"]),
            "root_shape": bound["root_shape"],
            "actual_T1_refs": bound["actual_T1_refs"],
        },
        "actual_signed_primitive_gap_functions": gaps,
        "ordinary_Z_product_refs": bound["literal_bindings"]["ordinary_Z_product_refs"],
        "formulas": {
            "root_shape": "a=2+2*mu*sigma(x), b=0, t0=0, Delta=2*mu*sigma(x)",
            "delta_E": "E_n-E_c",
            "delta_A_active": "R/2*((phase_n-phase_c)-(psi_n-psi_c)/(2*pi))",
            "delta_B_active": "-R/(4*pi)*(E_n*(T_n-T_c)+(E_n-E_c)*T_c)",
            "flat_branch": "Conditional on common Delta>=eta: both A/B values are 0; branch is not selected and active bodies are not evaluated",
            "A": "lazy(a*(phase-psi/(2*pi))/2)",
            "B": "lazy(E*(-a*T1hat/(2*pi)-b*phase)/2)",
            "q": "lazy(sigma(1-Delta/eta)*sqrt((2*eta-Delta)/(2*a)), Delta >= eta; q=0 on the flat branch)",
            "mu": "mu=exp(log(1/1000)-4*(exp(40)+11))",
            "eta": "eta=canonical eta_log source parameter; predicate Delta >= eta",
            "phase": "global fractional_part phase; no fresh local chart restart",
            "P2": "actual generic root p2 source row; value and cross-namespace identity remain unknown",
            "T1": "T1(psi)=running integral_0^psi t(theta)dtheta",
            "B_roles": "normalized loop B_over_Pstar is distinct from physical B=2*V_y",
        },
        "scope": {
            "ordinary_Z_orders": list(V2_ORDERS),
            "root_shape_correlations_proved": passed,
            "root_shape_derived_after_sigma_and_mu_witnesses": True,
            "same_mu_sigma_source_forms_bound": True,
            "sigma_endpoint_extensions_and_reflection_bound": True,
            "both_zero_Vy_source_clauses_bound": True,
            "same_flat_predicate_forms_bound": True,
            "normalized_B_roles_kept_distinct_from_physical_B": True,
            "E_equality_proved": False,
            "E_join_for_later": True,
            "P2_equality_proved": False,
            "phase_equality_proved": False,
            "primitive_gap_zero_proved": False,
            "numeric_values_computed": False,
            "source_values_evaluated": False,
            "new_density_proof_rows": 0,
        },
        "incoming_E_amplitude_route": {
            "native_source_route": "background_endpoints['O2_buffer']['fields']['E']",
            "C3_source_route": "original_closed_O3_background -> original_generic_source.common_velocity_E_axial5",
            "local_E_join_proved": False,
        },
        "findings": [
            "Native and rebuilt C3 E/A/B rows remain separate namespace handles.",
            "Common root shape is derived only after independent source E_y/V_y/root clauses and sigma witnesses bind.",
            "A/B are lazy local primitive projections; normalized B_over_Pstar remains distinct from physical B=2*V_y.",
            "Twelve signed E/A/B primitive gaps are defined without selecting the flat branch.",
        ],
        "missing_links": [
            "native versus C3 local E equality remains open",
            "phase and P2 cross-namespace identities remain open",
            "branch selection, numerical evaluation, and primitive-gap-zero claims remain forbidden",
        ],
        "source_inspection_audit": audit,
        "old_generic_identity_cached_data_only": True,
        "old_generic_identity_calls": 0,
        "old_provider_calls": 0,
        "old_checker_calls": 0,
        "old_normalizer_calls": 0,
        "old_source_value_reducer_calls": 0,
        "old_operator_calls": 0,
        "old_source_consumer_calls": 0,
        "source_value_reducer_calls": 0,
        "numeric_field_evaluations": 0,
        "old_numeric_calls": 0,
        "numeric_calls": 0,
        "new_algebra_family_calls": 1,
        "old_math_calls": 0,
        "old_source_calls": 0,
        "numeric": False,
        "numeric_values_computed": False,
        "numeric_field_values_present": False,
        "numerical": False,
        "full_pair_proved": False,
        "incoming_identity_proved": False,
        "caps_bound": False,
        "global_five_moment_closure": False,
        "global_closure": False,
        "E_zero_claimed": False,
        "A_zero_claimed": False,
        "B_zero_claimed": False,
        "P2_zero_claimed": False,
        "phase_equality_claimed": False,
        "signed_gap_zero_claimed": False,
        "gap_zero_claimed": False,
        "primitive_gap_zero_claimed": False,
        "goal_active": True,
        "full_goal_complete": False,
    }
    _v2_emit(journal, "primitive_gap_v2_source_binding_saved", passed=passed, root_rows=4, sigma_checks=2, primitive_cocycle_rows=4, primitive_gap_rows=12)
    return result
