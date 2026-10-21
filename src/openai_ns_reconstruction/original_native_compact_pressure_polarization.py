"""Actual compact pressure difference with shared physical measure and signs.

This leaf constructs only new expressions. Archived graphs and source theorems
are read as data; no archived transducer, normalizer or numerical owner runs.
The two compact fields remain independent source functions throughout.
"""
from fractions import Fraction
import ast
import hashlib
import json

PREFIX = 'lei_ren_part1_paper_compliant_current_original_'
PARENT = PREFIX + 'native_physical_scale_source_identity.json'
NATIVE = PREFIX + 'true_compact_C5_functions.json.gz'
OLD = PREFIX + 'C3_repair_band.json.gz'
PRE = PREFIX + 'n1_preRp_moment_functions.json.gz'
QUIET = PREFIX + 'true_quiet_Rp_C5_pressure_C3_functions.json.gz'
NATIVE_CHECK = PREFIX + 'true_compact_C5_functions_check.json'
OLD_CHECK = PREFIX + 'C3_repair_band_check.json'
JOIN = PREFIX + 'native_Rc_2Rc_pressure_join_lift.json'
NATIVE_SOURCE = PREFIX + 'true_compact_C5_functions.py'
OLD_SOURCE = PREFIX + 'C3_repair_band.py'
PRE_SOURCE = PREFIX + 'n1_preRp_moment_functions.py'
PHYSICAL_SOURCE = PREFIX + 'modified_Q4_pressure_C3_functions.py'
GEOMETRY_SOURCE = PREFIX + 'true_outer_C6_C5_functions.py'
PARTIAL_SOURCE = PREFIX + 'n1_partial_mean_Q.py'
NONRADIAL_SOURCE = PREFIX + 'n1_physical_nonradial_pressure_fields.py'
SCHEMA = 'original_native_compact_pressure_polarization_v1'
RECORD_HASHES = {'lei_ren_part1_paper_compliant_current_original_native_physical_scale_source_identity.json': 'b122d23659279f7f1f906756e44bbe2f189df8e31fc214c50c994f97f8d71760', 'lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz': '3593dd1ba3e510081090dd315d8581c1253afd5192233a284079c0646a556269', 'lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz': 'cf51e2131c3d1cb76f93a2deb76af12421e2deb75f7888bd310084b3b3eec193', 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz': '3921b728c3782d42802fef9e5a6f01bc244f65441811756d89dc98754aac0cce', 'lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions.json.gz': 'bbc8ae0361bdfbb75255707e7a234eb39f531533e6f39ad718cf262e27ad42c4', 'lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.py': 'f0b0333b19a0e2504d023b97e841033d3cbcc23427c57ad9672264866c01ec7b', 'lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions_check.json': '6c933d077f9c151eda246acd4d5f694d32dc8272fc3fb65a031b9b9cf28c8e7b', 'lei_ren_part1_paper_compliant_current_original_C3_repair_band.py': '62e2fa8598bf307126be0ffe0b4e23de83152104524db03f61b4944c225132e4', 'lei_ren_part1_paper_compliant_current_original_C3_repair_band_check.json': '00149597ab6206348a62dc39037ff53b7c5bcfec5d7eec6ae2591444f0f91a23', 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.py': '661933739091ce43f055813510aad016756c2094adf6a6ea2e652def7f129a2e', 'lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py': '7398b8cbba75453c566b7b7281d98827e128e89441124ea67b96d9493ca092a7', 'lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions.py': '4464a934d303a95a0da157e9e289f6ed8a04d5ff0adcb0d663857f0d9c57dcd3', 'lei_ren_part1_paper_compliant_current_original_native_Rc_2Rc_pressure_join_lift.json': '1f4b6aedf6683299117fb3bc2cddb90241a216f2b83296dddba04f60b53df39f', 'lei_ren_part1_paper_compliant_current_original_n1_partial_mean_Q.py': 'd826c150d021fa318201ce9b393b8efb92b459546e4ce8377c2e852f2d65cbc6', 'lei_ren_part1_paper_compliant_current_original_n1_physical_nonradial_pressure_fields.py': 'c2ae34e166a85542efd78c8768ed26191da90a945c17df1ea6b086d326cad3c7'}
DECODED_RECORD_HASHES = {'lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.json.gz': 'b6a7227da5e2427561cd22fcff614003cc14f3182537e1dbe278a6de84322510', 'lei_ren_part1_paper_compliant_current_original_C3_repair_band.json.gz': '845bade42d86d80e0ef82964e449c9167b18fc26c10f3849ae88768204f8698a', 'lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions.json.gz': '348df68bb426a2dfda6012e89f678b1e9c4ceb68b12f65f36853d5579977bc64', 'lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions.json.gz': '3bc0f9984a0aeb6f447d177c40c25116b8a59c3b90ee09b8661f029cbe032a40'}
STAGE_FLAGS = (
    'actual_two_compact_pressure_sources_bound',
    'common_physical_R_squared_log_measure_proved',
    'exact_compact_density_polarization_C1_proved',
    'exact_compact_m15_defect_representation_proved',
    'signed_terminal_active_plus_compact_density_defect_proved',
)
FALSE_FLAGS = (
    'native_compact_m15_source_identity_proved',
    'native_compact_V_m_source_identity_proved',
    'native_active_m15_source_identity_proved',
    'terminal_difference_zero_proved', 'old_error_terms_removed',
    'same_density_certificate', 'finite_N_history_identity', 'physical_Q_y_Z1',
    'global_five_moment_closure', 'legacy_feedback_reference_alias_proved',
    'prior_adapter_geometry_alias_proved', 'prior_adapter_P0_alias_proved',
    'time_order_recursion_claimed', 'numeric_values_computed', 'old_operators_replayed',
    'full_goal_complete', 'old212_adapter_source_identity_proved',
    'independent_terminal_function_match_complete',
    'native_vs_registered_preRp_m15_identity_proved',
)


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
        ensure_ascii=True, allow_nan=False).encode()).hexdigest()


def function(tree, name):
    return next((n for n in ast.walk(tree)
                 if isinstance(n, ast.FunctionDef) and n.name == name), None)


def syntax(tree, name, *parts):
    node = function(tree, name)
    text = ast.unparse(node) if node is not None else ''
    return all(part in text for part in parts)


def archived_dependency_ids(nodes, root):
    """Structural references only; no source operation is interpreted."""
    seen, pending = set(), [root]
    scalar_keys = {
        'negative': ('argument',), 'analytic_unary': ('argument',),
        'positive_quotient': ('numerator', 'denominator'),
        'function_substitution': ('expression', 'variable', 'value'),
        'definite_integral': ('integrand', 'lower', 'upper'),
    }
    while pending:
        ref = pending.pop()
        if ref in seen:
            continue
        seen.add(ref)
        row = nodes[ref]
        if row['operation'] in ('sum', 'product'):
            pending.extend(row['arguments'])
        else:
            pending.extend(row[key] for key in scalar_keys.get(row['operation'], ())
                           if isinstance(row.get(key), int))
    return seen


class NewExpressionGraph:
    """New exact C1 expressions; old source functions are opaque atoms."""

    def __init__(self):
        self.nodes, self.cache, self.poly_cache = [], {}, {}

    def node(self, operation, arguments=(), **data):
        row = dict(operation=operation, arguments=list(arguments), **data)
        key = json.dumps(row, sort_keys=True, separators=(',', ':'))
        if key not in self.cache:
            self.cache[key] = len(self.nodes)
            self.nodes.append(row)
        return self.cache[key]

    def c(self, value):
        q = Fraction(value)
        return self.node('exact_rational', numerator=q.numerator, denominator=q.denominator)

    def atom(self, name, **data):
        return self.node('new_formal_or_source_atom', name=name, **data)

    def add(self, *args):
        return self.node('sum', args)

    def mul(self, *args):
        return self.node('product', args)

    def neg(self, a):
        return self.mul(self.c(-1), a)

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def quotient(self, a, b):
        return self.node('positive_quotient', (a, b),
            certificate='same original L=1-delta*Z^2 positive on admitted axial domain')

    def jet(self, value, derivative=None):
        return (value, self.c(0) if derivative is None else derivative)

    def ja(self, *args):
        return tuple(self.add(*(p[k] for p in args)) for k in range(2))

    def jm(self, a, b):
        return (self.mul(a[0], b[0]),
                self.add(self.mul(a[1], b[0]), self.mul(a[0], b[1])))

    def jn(self, a):
        return tuple(self.neg(x) for x in a)

    def js(self, a, b):
        return self.ja(a, self.jn(b))

    def scale(self, value, a):
        return self.jm(self.jet(self.c(value)), a)

    def polynomial(self, ref):
        """Normalize only this fresh graph's ring operations, never old DAGs."""
        if ref in self.poly_cache:
            return self.poly_cache[ref]
        row = self.nodes[ref]
        if row['operation'] == 'exact_rational':
            result = {(): Fraction(row['numerator'], row['denominator'])}
        elif row['operation'] == 'sum':
            result = {}
            for arg in row['arguments']:
                for monomial, coefficient in self.polynomial(arg).items():
                    result[monomial] = result.get(monomial, Fraction(0)) + coefficient
        elif row['operation'] == 'product':
            result = {(): Fraction(1)}
            for arg in row['arguments']:
                out = {}
                for left, a in result.items():
                    for right, b in self.polynomial(arg).items():
                        key = tuple(sorted(left + right))
                        out[key] = out.get(key, Fraction(0)) + a * b
                result = out
        else:
            result = {(ref,): Fraction(1)}
        result = {m: c for m, c in result.items() if c}
        self.poly_cache[ref] = result
        return result


def derive(source_records, *, expected_input_digest):
    checks = {}

    def check(name, condition):
        checks[name] = bool(condition)

    check('frozen_input_digest', digest(source_records) == expected_input_digest)
    data, trees = {}, {}
    for name, pin in RECORD_HASHES.items():
        row = source_records[name]
        raw_ok = row['raw_sha256'] == pin and row['raw_authenticated_before_decode'] is True
        check('authenticated:' + name, raw_ok and
            row['decoded_sha256'] == DECODED_RECORD_HASHES.get(name, pin))
        if name.endswith('.py'):
            trees[name] = ast.parse(row['source_text'])
        else:
            data[name] = json.loads(row['source_text'])
    parent = data[PARENT]['component_result']
    native, old, pre, quiet = (data[k] for k in (NATIVE, OLD, PRE, QUIET))
    n = native['true_compact_C5_bindings']
    o = old['actual_C3_repair_band_functions']
    p = pre['defining_bindings']['compact']
    nc, oc = data[NATIVE_CHECK], data[OLD_CHECK]
    family = parent['source_family']
    check('same_source_family', all(v['source_family'] == family for v in (native, old, pre, quiet)))
    check('accepted_parent_common_geometry_and_scale', parent['passed'] is True
        and parent['native_original_absolute_Rc_Rp_identity_proved'] is True
        and parent['old_preRp_Pstar_source_identity_proved'] is True
        and parent['terminal_difference_reduced_to_active_compact'] is True)
    check('accepted_actual_native_compact_C5_and_old_C3', nc['all_passed'] is True
        and oc['all_passed'] is True
        and nc['checks']['actual_source']['passed'] is True
        and nc['checks']['actual_source']['actual_Rc_inlets_P0_N_mu_and_nonzero_pressure_history_retained'] is True)
    check('actual_old_compact_archive_consumed', pre['input_hashes'].get(OLD) == RECORD_HASHES[OLD])
    old_q_dependencies = archived_dependency_ids(pre['exact_graph_nodes'], p['physical_integrals']['I1_q'][0])
    check('actual_old_compact_Q_uses_same_delta_and_original_scale_port',
          {25220, 25203, o['original_band_variable']} <= old_q_dependencies)
    check('exact_old_band_pre_native_prefix_and_delta',
        pre['exact_graph_nodes'][:len(old['exact_graph_nodes'])] == old['exact_graph_nodes']
        and native['exact_graph_nodes'][:len(pre['exact_graph_nodes'])] == pre['exact_graph_nodes']
        and n['delta'] == 25220 and n['Pstar'] == 33572
        and p['common_Pstar_applied_once'] is True
        and p['actual_C3_original_band_variable'] == o['original_band_variable'])
    check('actual_common_x_domain_and_radius_pullback',
        o['actual_reserved_band'] == [1, 2]
        and n['background']['phase'] == 260181
        and native['exact_graph_nodes'][n['background']['radius']] == {
            'operation': 'product', 'arguments': [n['Rc'], n['background']['phase']]}
        and syntax(trees[PRE_SOURCE], 'compact_functions',
            'R = g.mul(Rc, x)', 'for edge in (g.one, g.constant(2))')
        and syntax(trees[NATIVE_SOURCE], '__init__', 'domain=[g.one, g.constant(2)]'))
    check('actual_native_prefix_in_quiet', quiet['exact_graph_nodes'][:len(native['exact_graph_nodes'])]
        == native['exact_graph_nodes'])
    check('native_physical_joins_accepted', data[JOIN]['component_result']['passed'] is True
        and data[JOIN]['component_result']['native_Rc_pressure_join_C3_proved'] is True
        and data[JOIN]['component_result']['native_2Rc_pressure_join_C3_proved'] is True)
    check('original_old_compact_pressure_source_syntax', syntax(trees[PRE_SOURCE], 'compact_functions',
        "raw['actual_C3_repair_band_functions']", "profile[0]['V']",
        "f['complete_histories']['m']", "profile[1]['V']", "f['complete_history_y']['m']",
        'recover_Q_increment(g, U, m, delta)', 'recover_Q_increment(g, Uy, my, delta)',
        "measure='physical dR=Rc*dx'", "'physical R*dR=Rc^2*x*dx'"))
    check('original_old_reduced_pressure_and_endpoints_syntax', syntax(trees[PRE_SOURCE], 'pressure_functions',
        "boundary('R_squared_q')", "i['I1_q']", "i['I1_U_q_Z']",
        "i['I1_U_squared']", "i['I1_U_U_Z']", "boundary('R_q_y')",
        'm15=a.sub(half(S1), boundary(\'R_q_y\'))'))
    check('original_old_Q_and_nonradial_product_source_syntax',
        syntax(trees[PARTIAL_SOURCE], 'recover_Q_increment',
            'C2Function(m.Z, m.ZZ, m.ZZZ)', 'alg.scale(alg.mul(z, low(V)), 2)',
            'alg.fixed(delta)', 'alg.neg(alg.mul(d, mz))', 'return alg.div(numerator, L)')
        and syntax(trees[NONRADIAL_SOURCE], 'nonradial_integrands',
            'I1_q=qr', 'I1_q_Z=qz', 'I1_q_squared=_z1_mul(g, qr, qr)',
            'I1_U_q_Z=uqz', 'I1_U_squared=uu', 'I1_U_U_Z=uuz'))
    check('native_physical_field_and_density_syntax', syntax(trees[PHYSICAL_SOURCE], 'physical_pressure_fields',
        "a.scale(fields['V'], Pstar)", "a.scale(fields['m'], Pstar)", 'my = a.sub(U, m)',
        'Qy = div(numerator(Uy, my))')
        and syntax(trees[PHYSICAL_SOURCE], 'reduced_pressure_S',
            'Qz = a.shift', "a.scale(Q, '1/2')", 'a.mul(U,', 'Qy'))
    check('native_weight_and_cumulative_pressure_syntax', syntax(trees[GEOMETRY_SOURCE], 'pressure_geometry',
        "(('S0', 1), ('S1', 2))", 'weight = g.mul(jac, *[R] * power)',
        "m15=a.sub(a.scale(cumulative['S1'], '1/2')", "initial['R_Q_y']"))
    check('native_compact_log_coordinate_and_full_history_syntax', syntax(trees[NATIVE_SOURCE], '__init__',
        "x = g.symbol('true_compact_repair_x')", "offset = g.unary('log', x)",
        "g.quotient(g.one, x", "self.modified_fields[k] = a.add(a.trim(fields[k]), self.history[k])",
        "self.radial = dict(E_y=f[1]['E'], V_y=f[1]['V'], m_y=a.sub(self.modified_fields['V'], self.modified_fields['m'])"))
    check('old_mass_radial_source_syntax', syntax(trees[OLD_SOURCE], 'build',
        "m=alg.add(alg.scale(complete['m'], -1), V)",
        'complete_history_y=dy', 'axial_Z_orders=[0, 1, 2, 3]'))
    check('actual_saved_compact_pressure_ports',
        pre['signed_pieces']['compact']['m15'] == [31155, 31156]
        and n['terminal']['pressure']['m15'][:2] == [262906, 262907]
        and n['inlet']['pressure']['m15'][:2] == [263147, 263148])

    g = NewExpressionGraph()
    zero, one = g.c(0), g.c(1)
    z = g.atom('Z', ordinary_Z_order=0)
    delta = g.atom('same_original_delta', Z_independent=True, native_source_node=n['delta'])
    rc = g.atom('common_positive_Rc', Z_independent=True, native_source_node=n['Rc'],
        alias_parent=RECORD_HASHES[PARENT])
    scale = g.atom('common_positive_Pstar', Z_independent=True, native_source_node=n['Pstar'],
        alias_parent=RECORD_HASHES[PARENT])
    x = g.atom('common_compact_x', domain=['1', '2'], Z_independent=True)
    Z, D = g.jet(z, one), g.jet(delta)
    d = g.js(g.jet(one), g.jm(Z, Z))
    L = g.js(g.jet(one), g.jm(D, g.jm(Z, Z)))
    invL = (g.quotient(one, L[0]),
            g.quotient(g.neg(L[1]), g.mul(L[0], L[0])))
    a = g.scale('1/2', g.js(g.jet(one), D))
    graph_digests = {name: digest(data[name]['exact_graph_nodes']) for name in (NATIVE, OLD, PRE, QUIET)}
    source_bindings = {}

    def ports(label, archive, ids, coordinate):
        ids = list(ids)
        source_bindings[label] = dict(archive=archive, graph_digest=graph_digests[archive],
            source_nodes=ids, native_coordinate=coordinate, common_coordinate=x,
            exact_pullback='same dimensionless x=R/Rc in [1,2]', ordinary_Z_rows=True)
        return tuple(g.atom(label + '_Z' + str(k), source_archive=archive,
            source_graph_digest=graph_digests[archive], source_node=node,
            ordinary_Z_order=k, coordinate_pullback=x, bound_or_midpoint=False)
            for k, node in enumerate(ids))

    physical = n['physical']
    native_x = n['background']['phase']
    Nu = ports('native_U', NATIVE, physical['U'][:2], native_x)
    Nq = ports('native_Q', NATIVE, physical['Q'][:3], native_x)
    Nqy = ports('native_Q_y', NATIVE, physical['Q_y'][:2], native_x)
    old_x = o['original_band_variable']
    old_v = ports('old_normalized_V', OLD, o['profiles_at_band_x']['radial_y_rows'][0]['V'][:3], old_x)
    old_m = ports('old_normalized_complete_m', OLD, o['complete_histories']['m'][:4], old_x)
    old_vy = ports('old_normalized_V_y', OLD, o['profiles_at_band_x']['radial_y_rows'][1]['V'][:2], old_x)
    old_my = ports('old_normalized_complete_m_y', OLD, o['complete_history_y']['m'][:3], old_x)
    Ou = tuple(g.mul(scale, q) for q in old_v)
    Om = tuple(g.mul(scale, q) for q in old_m)
    Ouy = tuple(g.mul(scale, q) for q in old_vy)
    Omy = tuple(g.mul(scale, q) for q in old_my)

    def recover_q(U, M, order):
        # Genuine ordinary derivatives, obtained from L*Q=numerator.
        zj = (z, one, zero)
        dj = (d[0], d[1], g.c(-2))
        lj = (L[0], L[1], g.mul(g.c(-2), delta))
        def product(A, B, k):
            coefficients = (1, 2, 1) if k == 2 else (1,) * (k + 1)
            return g.add(*(g.mul(g.c(coefficients[j]), A[j], B[k-j]) for j in range(k+1)))
        numerator = []
        for k in range(order+1):
            numerator.append(g.sub(g.sub(g.mul(g.c(2), product(zj, U, k)),
                g.mul(g.sub(one, delta), product(zj, M, k))), product(dj, M[1:], k)))
        result = []
        for k in range(order+1):
            value = numerator[k]
            for j in range(1, k+1):
                value = g.sub(value, g.mul(g.c(2 if k == 2 and j == 1 else 1), lj[j], result[k-j]))
            result.append(g.quotient(value, lj[0]))
        return tuple(result)

    Oq, Oqy = recover_q(Ou, Om, 2), recover_q(Ouy, Omy, 1)
    check('actual_Qzz_Qyz_Uz_rows_bound', len(Nq) == 3 and len(Nqy) == 2
        and len(Nu) == 2 and len(old_m) == 4 and len(old_my) == 3)

    def S(U, Q, Qy, Qz):
        A = g.ja(Q, Qy, g.jm(g.jm(a, Z), Qz))
        B = g.ja(g.scale(-2, g.jm(Z, Q)), g.jm(d, Qz),
                 g.scale(-2, g.jm(Z, Qy)))
        return g.ja(g.jm(invL, A), g.scale('1/2', g.jm(Q, Q)),
                    g.jm(Q, Qy), g.jm(invL, g.jm(U, B)))

    sn, so = S(Nu, Nq[:2], Nqy, Nq[1:3]), S(Ou[:2], Oq[:2], Oqy, Oq[1:3])
    du, dq, dqy, dqz = (g.js(v, w) for v, w in
        ((Nu, Ou[:2]), (Nq[:2], Oq[:2]), (Nqy, Oqy), (Nq[1:3], Oq[1:3])))
    Bn = g.ja(g.scale(-2, g.jm(Z, Nq[:2])), g.jm(d, Nq[1:3]),
               g.scale(-2, g.jm(Z, Nqy)))
    Bd = g.ja(g.scale(-2, g.jm(Z, dq)), g.jm(d, dqz),
               g.scale(-2, g.jm(Z, dqy)))
    polarized = g.ja(
        g.jm(invL, g.ja(dq, dqy, g.jm(g.jm(a, Z), dqz))),
        g.scale('1/2', g.jm(g.ja(Nq[:2], Oq[:2]), dq)),
        g.jm(Nq[:2], dqy), g.jm(Oqy, dq),
        g.jm(invL, g.jm(du, Bn)), g.jm(invL, g.jm(Ou[:2], Bd)))
    direct = g.js(sn, so)
    witnesses = [g.sub(direct[k], polarized[k]) for k in range(2)]
    check('new_exact_polarization_polynomials_zero_Z0_Z1',
          all(not g.polynomial(ref) for ref in witnesses))
    # The axial IBP coefficient is -2(1+delta)+4=2(1-delta).
    ibp = g.sub(g.add(g.mul(g.c(-2), g.add(one, delta)), g.c(4)),
                g.mul(g.c(2), g.sub(one, delta)))
    check('new_exact_axial_IBP_coefficient_identity', not g.polynomial(ibp))
    formal = {key: g.atom('universal_' + key) for key in ('U', 'Uy', 'Uz', 'm', 'mz')}
    def numerator(U, M, Mz):
        return g.sub(g.sub(g.mul(g.c(2), z, U),
            g.mul(g.sub(one, delta), z, M)), g.mul(d[0], Mz))
    continuity_left = g.add(numerator(formal['U'], formal['m'], formal['mz']),
        numerator(formal['Uy'], g.sub(formal['U'], formal['m']),
                  g.sub(formal['Uz'], formal['mz'])))
    continuity_right = g.sub(g.add(g.mul(g.add(one, delta), z, formal['U']),
        g.mul(g.c(2), z, formal['Uy'])), g.mul(d[0], formal['Uz']))
    continuity_witness = g.sub(continuity_left, continuity_right)
    check('new_exact_Q_plus_Qy_continuity_numerator_identity',
          not g.polynomial(continuity_witness))

    outputs = []
    parent_graph = parent['definition_graph']
    parent_ports = parent_graph['nodes']
    for k in range(2):
        weighted = g.mul(rc, rc, x, polarized[k])
        integral = g.node('new_exact_source_definite_integral', (weighted,),
            variable=x, lower='1', upper='2', ordinary_Z_order=k,
            measure='R^2 dy = R dR = Rc^2*x dx', geometry_Z_independent=True,
            genuine_compact_C1_sources_not_numeric_bounds=True)
        ends = [g.node('new_exact_source_endpoint_pullback', (dqy[k],),
            variable=x, value=edge, ordinary_Z_order=k) for edge in ('1', '2')]
        boundary = g.mul(rc, g.sub(g.mul(g.c(2), ends[1]), ends[0]))
        defect = g.sub(g.mul(g.c('1/2'), integral), boundary)
        current_ports = {}
        expected = {'native_quiet_Rc': (258845, 258846),
                    'native_quiet_2Rc': (262906, 262907),
                    'corrected_active_physical': (53466, 53470),
                    'old_compact_physical': (31155, 31156)}
        for role, handles in expected.items():
            rows = [row for row in parent_ports if row.get('namespace') == role
                    and row.get('ordinary_Z_order') == k]
            check('parent_source_port:' + role + ':Z' + str(k), len(rows) == 1
                and rows[0]['source_node'] == handles[k]
                and rows[0]['operation'] == 'source_function_port'
                and rows[0]['bound_or_midpoint'] is False)
            row = rows[0] if rows else {}
            current_ports[role] = g.atom(role + '_Z' + str(k), source_parent=PARENT,
                parent_graph_digest=parent_graph['nodes_digest'],
                parent_port=row.get('id'), source_node=handles[k], ordinary_Z_order=k)
        active = g.sub(current_ports['native_quiet_Rc'], current_ports['corrected_active_physical'])
        archived_compact_defect = g.sub(g.sub(current_ports['native_quiet_2Rc'],
            current_ports['native_quiet_Rc']), current_ports['old_compact_physical'])
        representation = g.node('new_qualified_exact_compact_pressure_representation',
            (archived_compact_defect, defect), admitted=all(checks.values()),
            premises=['actual full source fields and fixed Rc/Pstar pullback',
                'continuity numerator identity', 'two exact integration by parts identities',
                'accepted native pressure inlet join', 'new density polarization identity'],
            no_source_function_alias_assumed=True, ordinary_Z_order=k)
        terminal = g.add(active, defect)
        outputs.append(dict(ordinary_Z_order=k, density_difference=direct[k],
            polarized_density_difference=polarized[k], polynomial_zero_witness=witnesses[k],
            signed_compact_m15_defect=defect, weighted_integral=integral,
            endpoint_flux_difference=boundary, endpoint_coefficients=['1', '-2'],
            source_representation_certificate=representation,
            saved_compact_difference=archived_compact_defect, active_difference=active,
            terminal_with_explicit_compact_density_defect=terminal,
            compact_difference_zero_claimed=False,
            identity='native_2Rc-native_Rc-old_compact = Rc^2/2*int_1^2 x*DeltaS dx - Rc*(2*DeltaQy(2)-DeltaQy(1))',
            terminal_identity='F-T=(native_Rc-corrected_active)+compact_density_and_endpoint_defect'))

    binding_ok = all(checks.values())
    stages = dict(actual_two_compact_pressure_sources_bound=binding_ok,
        common_physical_R_squared_log_measure_proved=binding_ok,
        exact_compact_density_polarization_C1_proved=binding_ok,
        exact_compact_m15_defect_representation_proved=binding_ok,
        signed_terminal_active_plus_compact_density_defect_proved=binding_ok)
    result = dict(schema=SCHEMA, passed=binding_ok, source_only=True,
        ordinary_Z_orders=[0, 1], ordinary_Z_not_temporal_n=True,
        source_family=family, source_prerequisite_audit=checks,
        definition_graph=dict(namespace=SCHEMA, nodes=g.nodes, node_count=len(g.nodes),
            nodes_digest=digest(g.nodes), outputs=outputs, proof_graph_not_numeric_evaluator=True,
            archived_graphs_rebuilt=False, source_graph_digests=graph_digests),
        source_interface=dict(source_bindings=source_bindings,
            native_compact_physical=n['physical'], old_compact_pressure=p,
            source_field_differences=dict(U=du, Q=dq, Q_y=dqy, Q_Z=dqz),
            unaliased_incoming_background=n['incoming_background'],
            unaliased_incoming_correction=n['incoming_correction'],
            old_complete_history_retained=o['complete_histories'],
            old_compact_signed_m15=pre['signed_pieces']['compact']['m15'],
            native_pressure_local=n['pressure_local'],
            continuity_polynomial_zero_witness=continuity_witness,
            axial_IBP_polynomial_zero_witness=ibp,
            regularity_scope='actual C3/C5 functions on the original compact x and axial domain; mixed radial derivatives and continuous ordinary Z1 justify fixed-endpoint differentiation; no quantitative high-order norm claimed',
            integration_by_parts_identity='I1(U*U_y)=1/2*[R^2*U^2]-I1(U^2)',
            physical_Q_plus_Q_y_identity='((1+delta)*Z*U+2*Z*U_y-(1-Z^2)*U_Z)/(1-delta*Z^2)',
            measure='R^2 dy = R dR = Rc^2*x dx',
            pressure_not_P1=True, source_fields_not_identified=True,
            unchanged_parent_terminal_graph=data[PARENT]['component_result']['definition_graph']['nodes_digest']),
        retained_construction_data=parent['retained_construction_data'],
        unresolved_gates=[key for key, value in checks.items() if not value],
        remaining_construction_gates=['native_compact_V_m_incoming_and_control_function_identity',
            'native_Rc_corrected_active_physical_density_and_lower_endpoint_identity',
            'old212_same_source_packet_error_identity', 'global_five_moment_closure'],
        native_compact_m15_source_identity_proved=False,
        native_compact_V_m_source_identity_proved=False,
        native_active_m15_source_identity_proved=False,
        terminal_difference_zero_proved=False, old_error_terms_removed=False,
        legacy_numeric_chain_reconfigured=False, new_unentered_expression_simulations=0,
        runtime_model_metadata=None,
        source_inspection_audit=dict(new_source_consumer_calls=1, accepted_ancestor_operator_calls=0,
            old_operators_replayed=False, numeric_values_computed=False,
            archived_definition_graphs_rebuilt=False, old_pressure_IBP_reused_as_source_syntax=True,
            new_polynomial_certificates_only=True), **stages)
    for flag in FALSE_FLAGS:
        result[flag] = False
    return result
