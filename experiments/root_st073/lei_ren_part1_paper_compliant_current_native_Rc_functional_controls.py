"""Source-bound five-control Picard functions and live C1 iterate enclosures.

Exact bump integrals define the repair map. Directed integral/matrix ranges
enclose its finite iterates; no endpoint becomes a coefficient function.
Finite iterates are not fixed-point controls without a compatible N/tail.
"""
from dataclasses import dataclass
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_Rc_factored_source_oracle as roles
import lei_ren_part1_paper_compliant_current_native_paired_C1_transport as paired

source = roles.transport
repair = source.current.repair
HERE, PREFIX, sha = source.HERE, source.PREFIX, source.sha
packets = source.current.packets
ROWS, CONTROLS = repair.ROWS, repair.CONTROLS
NAME = PREFIX+'current_native_Rc_functional_controls.json'
RECEIPT = PREFIX+'current_native_Rc_functional_controls_check.json'
GATE = 'current_original_source_bound_five_control_Picard_C1_functions_and_finite_ranges_defined'
ep = packets.recovery.endpoints


def integral(g, integrand, variable, lower, upper, **metadata):
    g.ids([integrand, lower, upper])
    return g.node('definite_integral', integrand=integrand.node, variable=variable,
        lower=lower.node, upper=upper.node, exact_function_integral=True,
        numerical_value_not_installed=True, **metadata)


def exact_weights(g, mu):
    """Exact translated original bump integrals in t in[-1,1]."""
    L = g.unary('log', g.constant(2)); ell = g.mul(L, g.constant('1/40'))
    centers = [g.mul(L, g.constant(q)) for q in ('1/5', '1/2', '4/5')]
    variable = 'repair_normalized_bump_t'; t = g.symbol(variable)
    beta = g.node('compact_raw_beta', argument=t.node,
        definition='exp(-1/(1-t^2)) on abs(t)<1, zero otherwise',
        defining_module=PREFIX+'outer_pulse_map.py',
        defining_module_sha256=sha(PREFIX+'outer_pulse_map.py'))
    integrate = lambda value: integral(g, value, variable, g.constant(-1), g.one,
        measure='normalized original log-bump coordinate dt')
    normal = integrate(beta)
    ratio = lambda value: g.quotient(value, normal, 'original positive raw bump integral J0')
    w = g.mul(ell, t)
    powers = [g.constant('1/2'), g.neg(g.add(g.constant('1/2'), mu)),
        g.neg(g.add(g.constant('3/2'), mu))]
    H = [ratio(integrate(g.mul(beta, g.unary('exp', g.mul(p, w))))) for p in powers]
    Hax = ratio(integrate(g.mul(beta, g.unary('exp', g.neg(g.mul(mu, w))))))
    D = []
    for i in (0, 2):
        logx = g.add(centers[i], w)
        # (x^-mu-1)/mu = -log(x)*exprel(-mu*log(x)).
        value = g.mul(beta, g.neg(logx), g.unary('exprel', g.neg(g.mul(mu, logx))))
        D.append(ratio(integrate(value)))
    gram_den = g.mul(ell, normal, normal)
    base_gram = [g.quotient(integrate(g.mul(beta, beta, g.unary('exp', g.mul(g.constant(p), w)))),
        gram_den, 'original positive ell and J0^2') for p in ('-1/2', -1, -2)]
    gram = [[g.mul(base, g.unary('exp', g.mul(g.constant(p), ci))) for ci in centers]
        for base, p in zip(base_gram, ('-1/2', -1, -2))]
    E = [[g.mul(g.constant(sign), H[r], g.unary('exp', g.mul(powers[r], ci)))
        for ci in centers] for r, sign in enumerate((1, -1, 1))]
    gap = g.sub(centers[2], centers[0]); step = g.sub(centers[1], centers[0])
    axial_positive = g.mul(gap, Hax, g.unary('exp', g.neg(g.mul(mu, centers[0]))),
        g.unary('exprel', g.neg(g.mul(mu, gap))))
    # Three negative power differences and the negative S row cancel.
    swirl_positive = g.mul(*H, g.unary('exp', g.mul(g.add(*powers), centers[0])))
    for i, j in ((0, 1), (0, 2), (1, 2)):
        delta = g.sub(powers[i], powers[j])
        swirl_positive = g.mul(swirl_positive, g.unary('exp', g.mul(powers[i], step)),
            delta, step, g.unary('exprel', g.neg(g.mul(delta, step))))
    inverse_swirl = []
    for i in range(3):
        row = []
        for j in range(3):
            ii = [k for k in range(3) if k != j]; jj = [k for k in range(3) if k != i]
            cofactor = g.mul(g.constant((-1)**(i+j)), g.sub(
                g.mul(E[ii[0]][jj[0]], E[ii[1]][jj[1]]),
                g.mul(E[ii[0]][jj[1]], E[ii[1]][jj[0]])))
            row.append(g.quotient(cofactor, swirl_positive, 'original positive swirl Vandermonde determinant'))
        inverse_swirl.append(row)
    B = [[g.one, g.one, g.zero, g.zero, g.zero], [*D, g.zero, g.zero, g.zero]]
    B += [[g.zero, g.zero, *row] for row in E]
    return dict(mu=mu, log_length=L, ell=ell, centers=centers, normal=normal,
        powers=powers, H=H, Hax=Hax, D=D, cross=[gram[0][i] for i in (0, 2)],
        energy=gram[1], pressure=gram[2], B=B, inverse_swirl=inverse_swirl,
        axial_positive_determinant=axial_positive, swirl_positive_determinant=swirl_positive)


def inverse_action(g, W, vector):
    """Exact B inverse action with only theorem-positive denominators."""
    yM, yD, *swirl = vector; D0, D2 = W['D']; den = W['axial_positive_determinant']
    result = [g.quotient(g.sub(yD, g.mul(D2, yM)), den, 'original positive D0-D2'),
        g.quotient(g.sub(g.mul(D0, yM), yD), den, 'original positive D0-D2')]
    result += [g.add(*(g.mul(v, q) for v, q in zip(row, swirl))) for row in W['inverse_swirl']]
    return result


def quadratic(g, W, h):
    a0, a2, e0, e1, e2 = h; C = W['cross']; E = W['energy']; P = W['pressure']
    return [g.zero, g.quotient(g.add(g.mul(C[0], a0, e0), g.mul(C[1], a2, e2)),
        W['mu'], 'same original strictly positive mu'), g.zero,
        g.sub(g.add(g.mul(E[0], a0, a0), g.mul(E[2], a2, a2)),
            g.mul(g.constant('1/2'), g.add(*(g.mul(w, e, e) for w, e in zip(E, (e0, e1, e2)))))),
        g.mul(g.constant('1/2'), g.add(*(g.mul(w, e, e) for w, e in zip(P, (e0, e1, e2)))))]


def quadratic_action(g, W, h, v):
    a0, a2, e0, e1, e2 = h; b0, b2, f0, f1, f2 = v
    C = W['cross']; E = W['energy']; P = W['pressure']
    return [g.zero, g.quotient(g.add(g.mul(C[0], g.add(g.mul(e0, b0), g.mul(a0, f0))),
        g.mul(C[1], g.add(g.mul(e2, b2), g.mul(a2, f2)))), W['mu'], 'same original strictly positive mu'),
        g.zero, g.sub(g.mul(g.constant(2), g.add(g.mul(E[0], a0, b0), g.mul(E[2], a2, b2))),
            g.add(*(g.mul(w, e, f) for w, e, f in zip(E, (e0, e1, e2), (f0, f1, f2))))),
        g.add(*(g.mul(w, e, f) for w, e, f in zip(P, (e0, e1, e2), (f0, f1, f2))))]


def c1_vectors(g, values, jets):
    return [source.C1Function(v, j) for v, j in zip(values, jets)]


def exact_control_graph(built, *, iterations=3):
    if type(iterations) is not int or not 1 <= iterations <= 16:
        raise ValueError('Explicit finite Picard depth in[1,16] required')
    g = built['graph']; W = exact_weights(g, built['parameters']['mu']); N = built['N']
    d = [built['N_scaled_targets'][key] for key in ROWS]
    seq = [c1_vectors(g, [g.zero]*5, [g.zero]*5)]
    overN = lambda q: g.quotient(q, N, 'one common original positive integer N')
    for unused in range(iterations):
        h = seq[-1]; Q = quadratic(g, W, [q.value for q in h])
        QZ = quadratic_action(g, W, [q.value for q in h], [q.Z for q in h])
        next_h = inverse_action(g, W, [g.add(q.value, overN(v)) for q, v in zip(d, Q)])
        next_Z = inverse_action(g, W, [g.add(q.Z, overN(v)) for q, v in zip(d, QZ)])
        seq.append(c1_vectors(g, list(map(g.neg, next_h)), list(map(g.neg, next_Z))))
    final = seq[-1]; Q = quadratic(g, W, [q.value for q in final])
    QZ = quadratic_action(g, W, [q.value for q in final], [q.Z for q in final])
    residual = c1_vectors(g,
        [g.add(q.value, g.add(*(g.mul(b, v.value) for b, v in zip(row, final))), overN(r))
            for q, row, r in zip(d, W['B'], Q)],
        [g.add(q.Z, g.add(*(g.mul(b, v.Z) for b, v in zip(row, final))), overN(r))
            for q, row, r in zip(d, W['B'], QZ)])
    columns = [quadratic_action(g, W, [q.value for q in final],
        [g.one if j == i else g.zero for j in range(5)]) for i in range(5)]
    jacobian = [[g.add(W['B'][i][j], overN(columns[j][i])) for j in range(5)] for i in range(5)]
    # Physical profile handles use exactly the same coefficients; no seam or
    # recovered pressure is admitted merely by having these band recipes.
    variable = 'repair_x'; x = g.symbol(variable); logx = g.unary('log', x)
    bumps = []
    for center in W['centers']:
        t = g.quotient(g.sub(logx, center), W['ell'], 'original positive log bump width')
        raw = g.node('compact_raw_beta', argument=t.node,
            definition='exp(-1/(1-t^2)) on abs(t)<1, zero otherwise',
            defining_module=PREFIX+'outer_pulse_map.py',
            defining_module_sha256=sha(PREFIX+'outer_pulse_map.py'))
        bumps.append(g.quotient(raw, g.mul(W['ell'], W['normal'], x), 'original positive ell J0 x'))
    F = source.C1Function(overN(g.add(*(g.mul(b, final[i+2].value) for i, b in enumerate(bumps)))),
        overN(g.add(*(g.mul(b, final[i+2].Z) for i, b in enumerate(bumps)))))
    G = source.C1Function(overN(g.add(g.mul(bumps[0], final[0].value), g.mul(bumps[2], final[1].value))),
        overN(g.add(g.mul(bumps[0], final[0].Z), g.mul(bumps[2], final[1].Z))))
    power = g.unary('exp', g.neg(g.mul(g.add(g.constant('1/2'), W['mu']), logx)))
    Eprofile = g.c1mul(built['amplitude'], g.c1add(source.C1Function(power, g.zero), F))
    Vprofile = g.c1mul(built['amplitude'], G)
    return dict(**built, repair_weights=W, finite_picard_sequence=seq,
        control_residual=residual, implicit_Z_matrix=jacobian,
        implicit_Z_rhs=[g.neg(q.Z) for q in d], band_profiles=dict(F=F, G=G, E=Eprofile, V=Vprofile),
        band_variable=variable, band_bumps=bumps)


class FunctionEvaluator:
    """Explicit approximate scalar/range evaluator, with oracle-owned quadrature.

    Interval exprel returns a positive outer enclosure. Neither this routine
    nor an arbitrary injected quadrature certifies a control or residual.
    """
    def __init__(self, built, *, oracle, Z, N, ctx):
        if type(N) is not int or N < source.current.MIN_N:
            raise ValueError('One common integer N>=160 required')
        if oracle is None or not all(callable(getattr(oracle, key, None)) for key in ('parameter', 'source', 'integrate')):
            raise source.SourceOracleRequired('Explicit original parameter/source/integral oracle required')
        if getattr(oracle, 'source_family', None) != built['source_family']:
            raise source.SourceOracleRequired('Same original source family required')
        self.mode = getattr(oracle, 'mode', None)
        if self.mode not in ('original_source', 'synthetic_reference'):
            raise source.SourceOracleRequired('Explicit original_source or synthetic_reference mode required')
        if self.mode == 'original_source' and getattr(oracle, 'source_graph_sha256', None) != built['source_graph_sha256']:
            raise source.SourceOracleRequired('Accepted original graph hash required')
        self.built, self.oracle, self.Z, self.N, self.ctx = built, oracle, Z, N, ctx
        self.cache = {}

    def __call__(self, root, variables=None):
        g = self.built['graph']; g.ids([root]); variables = {} if variables is None else variables
        return self.walk(root.node, variables)

    def walk(self, index, variables):
        def strictly_positive(v):
            return ep(v)[0] > 0 if hasattr(v, '_mpi_') else v > 0
        def frozen(v):
            return getattr(v, '_mpi_', getattr(v, '_mpf_', v))
        key = index, tuple((k, frozen(v)) for k, v in sorted(variables.items()))
        if key in self.cache: return self.cache[key]
        row = self.built['graph'].nodes[index]; op = row['operation']; c = self.ctx
        child = lambda i: self.walk(i, variables)
        if op == 'exact_rational': value = c.mpf(row['numerator'])/row['denominator']
        elif op == 'bound_variable': value = variables[row['name']]
        elif op == 'original_source_parameter':
            value = self.oracle.parameter(row['name'])
            if row['name'] == 'mu' and not strictly_positive(value):
                raise ArithmeticError('Actual repair mu must be strictly positive')
        elif op == 'shared_positive_integer': value = c.mpf(self.N)
        elif op == 'sum': value = sum((child(i) for i in row['arguments']), c.mpf(0))
        elif op == 'product':
            value = c.mpf(1)
            for i in row['arguments']: value *= child(i)
        elif op == 'negative': value = -child(row['argument'])
        elif op == 'positive_quotient':
            denominator = child(row['denominator'])
            if not strictly_positive(denominator): raise ArithmeticError('Original positive denominator violated')
            value = child(row['numerator'])/denominator
        elif op == 'compact_raw_beta':
            t = child(row['argument'])
            value = repair.raw_beta(c, t) if hasattr(t, '_mpi_') else c.exp(-1/(1-t*t)) if abs(t) < 1 else c.mpf(0)
        elif op == 'analytic_unary':
            t = child(row['argument']); name = row['name']
            if name == 'exprel':
                value = repair.exp_average(c, t) if hasattr(t, '_mpi_') else c.mpf(1) if t == 0 else c.expm1(t)/t
            elif name == 'fractional_part': value = t-c.floor(t)
            else: value = getattr(c, name)(t)
        elif op == 'original_function_graph':
            value = self.oracle.source(row, coordinate=child(row['coordinate']), phase=child(row['phase']), Z=self.Z, N=self.N)
        elif op == 'definite_integral':
            value = self.oracle.integrate(lambda t: self.walk(row['integrand'], {**variables, row['variable']: t}),
                child(row['lower']), child(row['upper']))
        else: raise ValueError('Unknown exact function operation: '+op)
        if isinstance(value, (dict, repair.LogUpper)) or hasattr(value, 'scale'):
            raise TypeError('Source oracle returned a cap/factored range instead of an explicit function value')
        self.cache[key] = value
        return value


@dataclass(frozen=True)
class C1Enclosure:
    value: object
    Z: object


def range_quadratic(c, mu, logmu, matrix, h):
    """Enclose Q and DQ*h_Z using original signed ranges, not bound values."""
    a, b, e, f, k = h; C, E, P = (matrix[key] for key in ('cross_weights', 'energy_weights', 'pressure_weights'))
    zero = a.value.scalar(0); m = zero.scalar(mu)
    div = lambda x: x.positive_divide(m, logmu)
    sq = lambda v: v.value*v.value
    crossZ = lambda v, w: v.Z*w.value+v.value*w.Z
    values = [zero, div(a.value*e.value*C[0]+b.value*k.value*C[1]), zero,
        sq(a)*E[0]+sq(b)*E[2]-(sq(e)*E[0]+sq(f)*E[1]+sq(k)*E[2])*c.mpf('0.5'),
        (sq(e)*P[0]+sq(f)*P[1]+sq(k)*P[2])*c.mpf('0.5')]
    jets = [zero, div(crossZ(a, e)*C[0]+crossZ(b, k)*C[1]), zero,
        (a.value*a.Z*E[0]+b.value*b.Z*E[2])*2-(e.value*e.Z*E[0]+f.value*f.Z*E[1]+k.value*k.Z*E[2]),
        e.value*e.Z*P[0]+f.value*f.Z*P[1]+k.value*k.Z*P[2]]
    return [C1Enclosure(v, j) for v, j in zip(values, jets)]


class NativeRcFunctionalControls:
    def __init__(self, role_owner):
        if type(role_owner) is not roles.RoleBoundNativeRcFunctions:
            raise TypeError('Accepted original role-bound function owner required')
        self.role_owner, self.family = role_owner, role_owner.family
        self.target = role_owner.owner.owner; self.ctx = self.target.ctx
        checked = json.loads((HERE/roles.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(roles.GATE) or checked['source_family'] != self.family:
            raise ValueError('Checked same-family source-role binder required')
        self.hashes = {**role_owner.hashes, **checked['input_hashes'], roles.RECEIPT:sha(roles.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)}
        self.target.service.bind_hashes(self.hashes)

    def build(self, *, iterations=3):
        return exact_control_graph(self.role_owner.build(), iterations=iterations)

    @paired.native.inlet.source_precision
    def controls(self, Z_box, N, *, range_owner, source_ranges, iterations=3):
        """Return exact finite functions plus enclosures on the supplied live domain."""
        if type(iterations) is not int or not 1 <= iterations <= 16:
            raise ValueError('Explicit finite Picard depth in[1,16] required')
        if type(range_owner) is not paired.NativePairedC1Transport or range_owner.target is not self.target:
            raise TypeError('Same live accepted paired original source range owner required')
        got = source_ranges; c = self.ctx; rec = got['record']
        if type(N) is not int or N < source.current.MIN_N or rec['candidate_N'] != N:
            raise ValueError('Live range frequency must equal exact graph frequency')
        if c.mpf(Z_box)._mpi_ != c.mpf(rec['Z_box'])._mpi_ or rec['source_family'] != self.family:
            raise ValueError('Exact live source family and Z domain required')
        if got['target'] is None or len(got['cells']) != 24 or not rec['full24_original_C1_integral_range_transport_enclosed']:
            raise ValueError('Actual whole24 source histories required')
        checked = json.loads((HERE/paired.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(paired.GATE) or checked['source_family'] != self.family:
            raise ValueError('Checked paired whole-source prerequisite required')
        self.hashes.update({**range_owner.hashes, paired.NAME:sha(paired.NAME), paired.RECEIPT:sha(paired.RECEIPT)})
        self.target.service.bind_hashes(self.hashes)
        target = got['target']; mu = target['mu']
        if mu._mpi_ != self.target.repair_mu._mpi_ or ep(mu)[0] <= 0:
            raise ValueError('Same actual strictly positive source/repair mu required')
        W = repair.fresh_weights(c, mu, cells=512); B = repair.fresh_linear_inverse(c, mu, W)
        d = [C1Enclosure(target['values'][key]*N, target['Z_derivatives'][key]*N) for key in ROWS]
        if any(q.value.ctx is not c or q.Z.ctx is not c for q in d):
            raise ValueError('Original arithmetic context required')
        zero = d[0].value.scalar(0); seq = [[C1Enclosure(zero, zero) for unused in range(5)]]
        for unused in range(iterations):
            Q = range_quadratic(c, mu, target['logmu'], B, seq[-1])
            rhs = [C1Enclosure(q.value+r.value*(c.mpf(1)/N), q.Z+r.Z*(c.mpf(1)/N)) for q, r in zip(d, Q)]
            seq.append([C1Enclosure(-sum((q.value*b for q, b in zip(rhs, row)), zero),
                -sum((q.Z*b for q, b in zip(rhs, row)), zero)) for row in B['inverse_enclosure']])
        built = self.build(iterations=iterations)
        return dict(built=built, ranges=seq, exact_weight_enclosures=W, exact_matrix_enclosures=B,
            N_scaled_target_ranges=d, original_source_ranges=got, candidate_N=N, Z_box=c.mpf(Z_box))


def pair_roots(pairs, labels):
    return {key:dict(value=q.value.node, Z=q.Z.node) for key, q in zip(labels, pairs)}


def fixed_N_contraction_diagnostic(live):
    """A repair-only test at this N, never an all-N source certificate."""
    c = live['N_scaled_target_ranges'][0].value.ctx; matrix = live['exact_matrix_enclosures']
    target = live['original_source_ranges']['target']; const = lambda v: repair.LogUpper.constant(c, v)
    caps = [repair.LogUpper.add(c, [source.current.magnitude(q.value), source.current.magnitude(q.Z)])
        for q in live['N_scaled_target_ranges']]
    D = repair.LogUpper(c, c.mpf(max(ep(q.log)[1] for q in caps if q.log is not None)))
    CA = const(matrix['inverse_infinity_norm_upper'])
    rows = [const(sum(matrix['cross_weights'])).divide_positive(target['logmu']),
        const(matrix['energy_weights'][0]+matrix['energy_weights'][2]+sum(matrix['energy_weights'])/2),
        const(sum(matrix['pressure_weights'])/2)]
    CQ = repair.LogUpper(c, c.mpf(max(ep(q.log)[1] for q in rows)))
    rho = D*CA*const(2); Lipschitz = CA*CQ*rho*const(c.mpf(2)/live['candidate_N'])
    half = c.ln(c.mpf('0.5'))
    return dict(fixed_N_only=True, candidate_N=live['candidate_N'], C1_target_bound=D.record(),
        exact_inverse_bound=CA.record(), quadratic_C1_bound=CQ.record(), proposed_ball_radius=rho.record(),
        Lipschitz_upper=Lipschitz.record(), half_contraction_sufficient_bound_passed=ep(Lipschitz.log)[1] <= ep(half)[0],
        failed_sufficient_bound_does_not_prove_actual_map_diverges=True,
        target_bound_cannot_be_rescaled_to_a_different_N=True, globally_compatible_N_admitted=False)


@paired.native.inlet.source_precision
def run(role_owner, *, range_owner, source_ranges, iterations=3, return_live=False):
    began = time.monotonic(); owner = NativeRcFunctionalControls(role_owner)
    original = source_ranges['record']; live = owner.controls(original['Z_box'], original['candidate_N'],
        range_owner=range_owner, source_ranges=source_ranges, iterations=iterations)
    built = live['built']; W = built['repair_weights']; c = owner.ctx
    norms = []
    for row in live['ranges']:
        caps = [repair.LogUpper.add(c, [source.current.magnitude(q.value), source.current.magnitude(q.Z)]) for q in row]
        logs = [ep(q.log)[1] for q in caps if q.log is not None]
        norms.append(dict(control_C1_ranges={key:dict(value=q.value.record(), Z=q.Z.record()) for key, q in zip(CONTROLS, row)},
            C1_norm_log_upper=None if not logs else c.mpf(max(logs))))
    result = dict(source_family=owner.family, **{GATE:True}, candidate_N=live['candidate_N'], Z_box=live['Z_box'],
        exact_function_graph_nodes=built['graph'].nodes, original_source_roles=built['source_roles'],
        exact_N_scaled_target_roots=pair_roots([built['N_scaled_targets'][key] for key in ROWS], ROWS),
        exact_finite_Picard_function_roots=[pair_roots(row, CONTROLS) for row in built['finite_picard_sequence']],
        exact_residual_roots=pair_roots(built['control_residual'], ROWS),
        exact_implicit_Z_matrix_roots=[[v.node for v in row] for row in built['implicit_Z_matrix']],
        exact_implicit_Z_rhs_roots=[v.node for v in built['implicit_Z_rhs']],
        exact_repaired_band_profile_roots=pair_roots(list(built['band_profiles'].values()), list(built['band_profiles'])),
        exact_weight_roots=dict(mu=W['mu'].node, J0=W['normal'].node, D=[q.node for q in W['D']],
            H=[q.node for q in W['H']], cross=[q.node for q in W['cross']], energy=[q.node for q in W['energy']],
            pressure=[q.node for q in W['pressure']], axial_positive=W['axial_positive_determinant'].node,
            swirl_positive=W['swirl_positive_determinant'].node),
        defining_equation='B_exact(mu)*h+N*r(N,Z)+Q_exact(mu,h)/N=0',
        differentiated_iterate='h_next_Z=-B_exact^-1*(N*r_Z+DQ(h)*h_Z/N)',
        correct_cross_weight='integral sqrt(x)*g_i(x)^2 dx; log-coordinate exponent -1/2',
        row_order=ROWS, control_order=CONTROLS, finite_Picard_depth=iterations,
        live_original_C1_iterate_range_records=norms, exact_weight_enclosures=live['exact_weight_enclosures'],
        exact_matrix_enclosures=live['exact_matrix_enclosures'],
        fixed_N_repair_contraction_diagnostic=fixed_N_contraction_diagnostic(live),
        original_source_C1_target_ranges_used_not_log_caps=True, inherited_original_continuous_cells=24,
        function_coefficients_never_selected_from_enclosure_endpoints=True,
        original_P0_and_P0_Z_unchanged=True, actual_source_graph_ancestor_constructors_called=False,
        numerical_original_source_point_oracle_installed=False, finite_iterates_are_not_solved_controls=True,
        certified_fixed_point_tail_installed=False, actual_five_controls_installed=False,
        actual_terminal_Z_function_closure_installed=False, current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN, False), input_hashes=owner.hashes, execution_seconds=time.monotonic()-began,
        scope='Exact original N*r-bound five-control integral matrix/inverse, nonlinear action, C1 Picard functions, band profile recipes and actual same-source finite iterate enclosures. No contraction/tail, solved controls, terminal closure, global N or recursion admitted.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result), indent=2)+'\n', encoding='utf8')
    print('Source-bound five-control C1 Picard map and live finite ranges:', iterations, 'iterations; fixed point open', flush=True)
    return (result, owner, live) if return_live else result
