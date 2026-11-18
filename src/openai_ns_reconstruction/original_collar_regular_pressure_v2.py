"""Canonical regular pressure and actual cut swirl at the saved q=2 inlet."""
from __future__ import annotations

import hashlib

import sympy as sp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext

from .original_signed_density_expression import SignedDensityExpressionGraph
from .original_collar_regular_pressure_point import enclose

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
LOWER = PREFIX + 'n1_lower_collar_functions.json.gz'
CUTOFF = PREFIX + 'n1_inner_cutoff_functions.json.gz'
FORWARD_CODE = 'src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py'
FTC = 'original_collar_signed_ftc_v2_receipt.json.gz'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def pair_mul(a, b):
    return a[0]*b[0], a[1]*b[0]+a[0]*b[1]


def pair_div(a, b):
    value = a[0]/b[0]
    return value, (a[1]-value*b[1])/b[0]


def _proof():
    r, s, z, h, eps, delta = sp.symbols('rho s Z hb epsilon delta')
    p, mean = sp.Function('Psi')(r), sp.Function('meanPsi')(r)
    # The true own mean satisfies (rho*meanPsi)_rho=Psi.
    H_rho = sp.expand(sp.diff(r*(p-mean), r)).subs(sp.diff(mean, r), (p-mean)/r)
    require(sp.simplify(H_rho-r*sp.diff(p,r)) == 0, 'own mean FTC failed')
    T, K, M = (sp.Function(name)(s,z) for name in ('T','K','D'))
    sigma = sp.Function('sigma')(s)
    L, d = 1-delta*z*z, 1-z*z
    chi = 1-(1-h)*sigma
    bracket = (1-delta)*z*T+d*sp.diff(T,z)
    H = sp.Function('rho_Psi_minus_mean')(s,z)
    D = ((1-delta)*z*H+d*sp.diff(H,z))/L
    D_s = sp.diff(D,s).subs({sp.diff(H,s,z): h*chi*sp.diff(T,z), sp.diff(H,s): h*chi*T})
    require(sp.cancel(D_s-h*chi*bracket/L) == 0, 'ordinary axial derivative of mean history changed')
    B = 2*z*chi*T/L-M
    direct = -eps**2*sp.diff(B,s).subs({sp.diff(T,s): h*K, sp.diff(M,s): h*chi*bracket/L})
    visc = eps**2*h*chi*(bracket-2*z*K)/L
    singular = eps**2*2*z*(1-h)*T*sp.diff(sigma,s)/L
    residual = sp.expand(direct-visc-singular)
    require(sp.cancel(residual) == 0 and sp.cancel(sp.diff(residual,z)) == 0,
            'regular and singular ordinary Z0/Z1 identity failed')
    q1, q2, q3 = sp.symbols('b_r b_rr raw_N1p')
    # Q_y=epsilon*r*b_r, and N1p uses (Q_yy+Q_y)/R.
    Qy, Qyy = eps*r*q1, eps*(r*q1+r*r*q2)
    require(sp.cancel(-eps**2*h*r*r*(Qyy+Qy)/(eps*r)+eps**2*h*r*(2*r*q1+r*r*q2)) == 0,
            'physical R*dR Jacobian changed')
    F0, phi, Fcut, P2 = sp.symbols('F0 Phi F1_cut P2')
    fcut, S0 = eps*Fcut/F0, eps**2*F0**2
    old = -eps**2*h*r*r*P2*(q3/P2+2*S0*phi*fcut/(eps**3*P2))
    expect = -eps**2*h*r*r*(q3+2*F0*phi*Fcut)
    require(sp.cancel(old-expect) == 0, 'cut swirl unit conversion changed')
    return dict(own_mean_FTC=True, ordinary_axial_derivative_of_mean_FTC=True, physical_weight=True,
                singular_partition_C0_C1=True,
                cut_swirl_F1_conversion=True, ordinary_Z_orders=[0,1])


def _source_graph(inputs, pins):
    low, cut, check, ftc, q2 = (inputs[name] for name in ('lower','cutoff','cut_check','ftc','q2'))
    native = inputs['native']['actual_native_collar_transport']
    family = low['source_family']
    require(all(p['source_family']==family for p in (cut,check,q2,inputs['native']))
            and check['all_passed'] and check['checks']['source']['true_regular_solution_and_axis_jets_used']
            and check['input_hashes'][CUTOFF] == pins['cutoff'] and native['passed']
            and native['actual_quotient_equals_ratio'] and native['domain']['canonical_q_interval']==[0,2]
            and native['native_defining_leading_function_equals_canonical_collar'], 'different actual collar source family')
    old_partition = q2['actual_collar_pressure_partition']
    accepted_cut = q2['actual_complete_first_collar_m15_increment']
    require(old_partition['full_pressure_partition_for_same_leading_fields']
            and accepted_cut['passed'] and accepted_cut['raw_allword_tail_added_once_per_field_at_evaluation_point'],
            'actual pressure and cutoff source prerequisites missing')
    lo_nodes, nodes = low['exact_graph_nodes'],cut['exact_graph_nodes']
    require(nodes[:len(lo_nodes)] == lo_nodes, 'cutoff archive no longer extends the exact lower-collar graph')
    lower = low['lower_collar_bindings']; f,p = lower['fields'],lower['parameters']
    saved = cut['inner_cutoff_bindings']; kfirst = saved['fields']['first']; radial = lower['core_radial_rows']['Psi'][1]
    require(saved['geometry']['keep_fraction']=='1/2' and saved['geometry']['zero_fraction']=='3/4'
            and saved['geometry']['Z_independent_cutoff'] and saved['definition']['leading_pressure_forcing_unchanged']
            and f['domain']=={'lower':p['sc_quarter'],'upper':p['sc_half']}, 'original cutoff support or collar changed')
    eta = saved['cutoff_rows']['first'][0][0]
    cutF = kfirst['rows']['F1'][0][:2]
    rawF = []
    for order, index in enumerate(cutF):
        product = nodes[index]
        require(product['operation']=='product' and len(product['arguments'])==2
                and eta in product['arguments'], 'cut F1 is not the actual product with eta')
        raw = next(x for x in product['arguments'] if x!=eta)
        row = nodes[raw]
        require(row['operation']=='current_original_unique_regular_n1_component' and row['quantity']=='F1'
                and row['chart']=='first' and row['coordinate']==f['query']
                and row['coordinate_domain_upper']==p['sc_quarter']
                and row['ordinary_slow_Z_derivative_order']==order
                and row['all_words_retained_in_analytic_limit'], 'cut product lost an actual raw all-word F1 port')
        rawF.append(raw)
    previous = ftc['actual_FTC']
    require(ftc['passed'] and ftc['raw_input_sha256'][f'experiments/root_st073/{LOWER}']==pins['lower']
            and previous['radial_FTC']['density_C0_C1']==[131,134]
            and previous['signed_primitive']['integral_C0_C1']==[157,158]
            and previous['limits']['upper']==p['sc_half'], 'different K and signed A-A0 source')

    g=SignedDensityExpressionGraph()
    def ref(index):
        return g.reference(CUTOFF,pins['cutoff'],index,source_nodes=nodes)

    s, rho, rho_core = ref(f['query']), ref(f['rho']), ref(35558)
    lo, upper = g.zero, ref(p['sc_half'])
    epsilon, hb, delta, z = (ref(index) for index in (p['epsilon'],p['hb'],p['delta'],20))
    chi, mass = ref(f['chi']), epsilon**2*hb*rho**2
    L, L1, d, d1 = 1-delta*z*z, -2*delta*z, 1-z*z, -2*z
    def at(expression, variable, where):
        return g.insert(dict(operation='function_substitution', expression=expression.node,
                             variable=variable.node, value=where.node))

    coordinate = g.insert(dict(operation='bound_variable',name='regular_pressure_integration_s'))
    def integral(pair, end, label):
        return tuple(g.insert(dict(operation='definite_integral', variable=coordinate.node, lower=lo.node,
                    upper=end.node, integrand=at(expr,s,coordinate).node,
                    ordinary_Z_order=i, endpoint_independent_of_Z=True, label=label)) for i,expr in enumerate(pair))

    def ftc_ref(i, ordinary, endpoint=None):
        row = dict(operation='saved_exact_source_reference', report=FTC, report_sha256=pins['ftc'],
                   local_node=i, ordinary_Z_order=ordinary)
        if endpoint is not None:
            row.update(bind_local_variable=previous['source_factor']['variable'], endpoint_node=endpoint.node)
        else:
            row.update(bind_local_variable=previous['source_factor']['variable'], endpoint_node=s.node)
        return g.insert(row)

    q, q1, q2 = (ref(i) for i in f['ratio'][:3])
    a, a1, a2 = (at(ref(i),rho_core,rho) for i in radial[:3])
    T, T1, T2 = rho**2*q*a, rho**2*(q1*a+q*a1), rho**2*(q2*a+2*q1*a1+q*a2)
    K,K1=(ftc_ref(i,k) for k,i in enumerate(previous['radial_FTC']['density_C0_C1']))
    b = (1-delta)*z*T+d*T1-2*z*K
    b1 = (1-delta)*(T+z*T1)+d1*T1+d*T2-2*(K+z*K1)
    visc = tuple(epsilon**2*hb*chi*v for v in pair_div((b,b1),(L,L1)))
    inertia = tuple(mass*value/2 for value in map(ref,lower['pressure_bindings']['S']))
    swirl = tuple(-2*mass*v for v in pair_mul(tuple(map(ref,f['F_rho_rows'][0][:2])),tuple(map(ref,cutF))))
    regularVisc=integral(visc,upper,'regular_viscous'); flow=integral(inertia,upper,'leading_inertial')
    cutoffPiece=integral(swirl,ref(p['sc_quarter'])*g.constant('3/4'),'cut_swirl_only_through_support')
    row4=tuple(a+b+c for a,b,c in zip(regularVisc,flow,cutoffPiece,strict=True))
    principal=tuple(ftc_ref(i,k,upper) for k,i in enumerate(previous['signed_primitive']['principal_epsilon2_A0_sigma_C0_C1']))
    remainder=tuple(ftc_ref(i,k,upper) for k,i in enumerate(previous['signed_primitive']['integral_C0_C1']))
    direct=integral(tuple(map(ref,lower['pressure_bindings']['signed_physical_integrands']['m15'])),upper,'original_direct_N1p')
    properties=dict(row4_regular_plus_cut=[v.node for v in row4], row5_principal=[v.node for v in principal],
                    row6_signed_width=[v.node for v in remainder],
                    source_sum_4_5_6=[(a+b+c).node for a,b,c in zip(row4,principal,remainder,strict=True)],
                    direct_N1p_plus_cut=[(p+c).node for p,c in zip(direct,cutoffPiece,strict=True)],
                    regular_visc=[v.node for v in regularVisc],leading_inertial=[v.node for v in flow],
                    cut_swirl=[v.node for v in cutoffPiece],cut_F1_source_nodes=cutF, raw_F1_source_nodes=rawF,
                    leading_F_source_nodes=f['F_rho_rows'][0][:2], cutoff_eta_source_node=eta,
                    exact_graph_prefix_rows=len(lo_nodes), terms_are_function_definitions=True,
                    legacy_forward_finite_response_equals_all_word=False)
    return dict(ports=properties, **g.packet())


def _diagnose(forward):
    require(forward['passed'] and forward['source_principal_recipe_linked'], 'missing accepted singular comparison')
    iv,mp=MPIntervalContext(),MPContext()
    iv.prec=1800
    def read(row):
        return iv.make_mpf(tuple(tuple(row[f'{side}_exact_mpf_tuple']) for side in ('lower','upper')))
    rows={}
    for z in ('-.3','.3'):
        pair=forward['samples'][z]['ordinary_Z'][1]
        exact,old=map(read,(pair['canonical_principal'],pair['archived_principal_coefficient']))
        lo,hi=max(exact.a,old.a),min(exact.b,old.b)
        edge=lambda a:mp.nstr(mp.make_mpf(a._mpi_[1]),20)
        rows[z]=dict(exact_enclosure_overlap=lo<=hi, old_width_upper=edge((old.b-old.a).b),
                     canonical_width_upper=edge((exact.b-exact.a).b),
                     missing_lower_canonical_slice=edge(max(iv.mpf(0),(old.a-exact.a).b)),
                     missing_upper_canonical_slice=edge(max(iv.mpf(0),(exact.b-old.b).b)),
                     disjoint_gap_lower=edge(max(iv.mpf(0),(lo-hi).a)),
                     actual_value_outside_old_proved=lo>hi)
    return rows


def evaluate(*, inputs, pins):
    q2=inputs['q2']
    recipe=inputs['forward_code']
    require(hashlib.sha256(recipe.encode()).hexdigest()==q2['implementation_hashes'][FORWARD_CODE]
            and "pressure=(-(Np+coupling)*vt['rho']*vt['rho']).truncate(1)" in recipe
            and "v['e']**3*v['P2']" in recipe and "leading_N1p_after_velocity_support_retained=True" in recipe,
            'old forward fourth recipe changed')
    graph=_source_graph(inputs,pins)
    proof=_proof()
    point=enclose(q2)
    return dict(passed=True, symbolic_source_proof=proof, canonical_function_bridge=graph, selected_point=point,
                fifth_row_Z1_interval_diagnostics=_diagnose(inputs['forward']),
                original_measure='epsilon^2*hb*rho^2 ds', cutoff_support_q='[0,3/4]',
                full_forward_six_equal_canonical=False, forward212_terms_added_or_removed=0,
                old_error_caps_removed=False, independent_terminal_matching_closed=False,
                ancestor_mathematics_replayed=False)
