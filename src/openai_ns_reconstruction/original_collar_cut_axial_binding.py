"""Bind only the actual all-word axial cutoff contribution to m15."""
from __future__ import annotations

from .original_signed_density_expression import SignedDensityExpressionGraph


CUTOFF = 'lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions.json.gz'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def bind(*, cutoff, lower, q2, cutoff_sha):
    family = q2['source_family']
    history = q2['actual_n1_cutoff_history_identity']
    native = q2['actual_native_collar_transport']
    velocity = q2['actual_complete_first_collar_m15_increment']
    require(cutoff['source_family']==lower['source_family']==family==history['source_family']==velocity['source_family']
            and q2['input_hashes'][CUTOFF]==cutoff_sha and q2['checks']['passed']
            and history['passed'] and history['graph_binding']['checks']['eta_graph']
            and native['passed'] and native['native_defining_leading_function_equals_canonical_collar']
            and velocity['passed'] and velocity['raw_allword_tail_added_once_per_field_at_evaluation_point'],
            'current actual cutoff and leading source packets required')

    original, lo_nodes = cutoff['exact_graph_nodes'], lower['exact_graph_nodes']
    require(original[:len(lo_nodes)]==lo_nodes, 'different native leading collar graph')
    low = lower['lower_collar_bindings']
    f, p = low['fields'], low['parameters']
    cut = cutoff['inner_cutoff_bindings']
    eta = cut['cutoff_rows']['first'][0][0]
    c0c1 = cut['fields']['first']['rows']['Uz1'][0][:2]
    raw = [row['node'] for row in history['graph_binding']['raw_first']['fields']['Uz1'][:2]]
    leading = f['U_rho_rows'][0][:2]
    require(cut['geometry']['keep_fraction']=='1/2'
            and cut['geometry']['zero_fraction']=='3/4'
            and cut['geometry']['phase_is_sc_quarter']
            and cut['geometry']['Z_independent_cutoff']
            and cut['cutoff_rows']['first'][0]==[history['graph_binding']['eta']['nodes']['eta'],0,0,0]
            and leading==native['canonical_function_handles']['U_C0_C3'][:2]
            and history['graph_binding']['raw_first']['coordinate_domain_upper']==p['sc_quarter'],
            'different axial source, phase or support')
    for order,(j,k) in enumerate(zip(c0c1,raw,strict=True)):
        product, state = original[j],original[k]
        require(product['operation']=='product' and sorted(product['arguments'])==sorted((eta,k))
                and state['operation']=='current_original_unique_regular_n1_component'
                and state['quantity']=='Uz1' and state['chart']=='first'
                and state['coordinate']==f['query'] and state['coordinate_domain_upper']==p['sc_quarter']
                and state['ordinary_slow_Z_derivative_order']==order
                and state['all_words_retained_in_analytic_limit'],
                'cut axial port must multiply the actual raw all-word response')

    q2_axial = q2['actual_correlated_q2_inlet']
    terms = {}
    for order in ('0','1'):
        part = q2_axial['ordinary_Z'][order]['source_leaves']['cut_axial']
        source = velocity['ordinary_Z'][order]['cut_axial_velocity']
        require(part['signed_approximation']==source['signed_inlet_product_integral']
                and part['error_radius']==source['error_radius'], 'selected cut-axial packet substituted another error')
        terms[order] = part
    require(q2['actual_cutoff_complete_support_measure']['exact_unweighted_q_mass']=='1/2+(1/4)*(1/2)=5/8',
            'unknown complete cutoff mass')

    g = SignedDensityExpressionGraph()
    def ref(node):
        return g.reference(CUTOFF,cutoff_sha,node,source_nodes=original)

    U,UZ = map(ref,leading)
    u,uZ = map(ref,c0c1)
    epsilon,hb,rho = map(ref,(p['epsilon'],p['hb'],f['rho']))
    jacobian = hb*rho
    density = (2*epsilon*jacobian*U*u, 2*epsilon*jacobian*(UZ*u+U*uZ))
    s = ref(f['query'])
    bvar = g.insert(dict(operation='bound_variable',name='m15_cut_axial_s'))
    upper = g.constant('3/4')*ref(p['sc_quarter'])
    integral = []
    for k,body in enumerate(density):
        integrand = g.insert(dict(operation='function_substitution',expression=body.node,
                                  variable=s.node,value=bvar.node))
        integral.append(g.insert(dict(operation='definite_integral',integrand=integrand.node,
            variable=bvar.node,lower=g.zero.node,upper=upper.node,
            ordinary_Z_order=k,endpoint_independent_of_Z=True)).node)
    return dict(passed=True, **g.packet(), raw_all_word_Uz1_C0_C1=raw,
                cut_Uz1_C0_C1=c0c1, leading_U_C0_C1=leading,
                actual_axial_integral_C0_C1=integral, radial_density_nodes=[r.node for r in density],
                bounds_of_support_in_s=['0','3*sc/16'],
                physical_radial_density='2*epsilon*hb*rho*Ulead*Uz1cut',
                legacy_unit='epsilon^2*hb*(sc/4)',
                same_unit_recipe='u_cut=epsilon*Uz1cut; ds=(sc/4)*dq',
                selected_q2_signed_and_error_terms=terms, selected_packet_only_at_saved_Z=True,
                canonical_whole_Z_function_defined=True, all_Z_numeric_coverage=False,
                legacy_finite_Uz1_equals_raw_all_word=False, prior_ledger_mutated=False,
                independent_terminal_matching_closed=False, ancestor_calculations_replayed=False)
