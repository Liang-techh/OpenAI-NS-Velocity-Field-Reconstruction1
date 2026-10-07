"""Independent original bump/physical-moment checks and live Picard bindings.

Manufactured coefficient/reference tests validate the operator. They do not
evaluate the current original source field or establish terminal closure.
"""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_native_Rc_functional_controls as current
import lei_ren_part1_paper_compliant_current_native_signed_input_enclosures as signed

HERE, PREFIX, sha = current.HERE, current.PREFIX, current.sha
source, repair = current.source, current.repair
ep = current.ep


def near(c, got, expected, tol='1e-43'):
    if abs(got-expected) > c.mpf(tol)*(1+abs(expected)):
        raise AssertionError('Independent integral/action mismatch')


def scalar_fixture(owner, mu, c):
    g = source.FunctionTransportGraph()
    mu_root = g.node('original_source_parameter', name='mu', exact_source='manufactured positive mu reference only')
    built = dict(graph=g, source_family=owner.family, source_graph_sha256='manufactured_operator_reference_only')
    class ReferenceOracle:
        source_family = owner.family
        mode = 'synthetic_reference'
        def parameter(self, name):
            if name != 'mu': raise AssertionError('Unexpected original-field request in bump-only fixture')
            return c.mpf(mu)
        def source(self, *args, **kwargs): raise AssertionError('Bump checker must not manufacture current densities')
        def integrate(self, function, left, right):
            return c.quad(function, [left, (left+right)/2, right])
    W = current.exact_weights(g, mu_root)
    evaluate = current.FunctionEvaluator(built, oracle=ReferenceOracle(), Z=c.mpf(0), N=2048, ctx=c)
    return g, W, evaluate


def direct_physical_reference(owner):
    """Direct x-coordinate integrals of original physical increment densities."""
    c = mp.mp.clone(); c.dps = 75; L = c.log(2); ell = L/40
    centers = [L/5, L/2, 4*L/5]
    raw = lambda t: c.exp(-1/(1-t*t)) if abs(t) < 1 else c.mpf(0)
    normal = c.quad(raw, [-1, 0, 1])
    def bump(i, x): return raw((c.log(x)-centers[i])/ell)/(ell*normal*x)
    def integrate(function):
        return sum(c.quad(function, [c.exp(ci-ell), c.exp(ci), c.exp(ci+ell)]) for ci in centers)
    h = list(map(c.mpf, ('.031', '-.019', '.023', '-.017', '.011')))
    hZ = list(map(c.mpf, ('-.021', '.014', '.007', '.009', '-.012')))
    checks = weights_checked = inverse_checked = range_checked = 0
    for mu_text in ('.04', '1e-40', '.166666666666'):
        mu = c.mpf(mu_text); g, W, value = scalar_fixture(owner, mu_text, c)
        B = c.matrix([[value(v) for v in row] for row in W['B']])
        axial_det = value(W['axial_positive_determinant']); swirl_det = value(W['swirl_positive_determinant'])
        near(c, axial_det, B[1,0]-B[1,1]); near(c, swirl_det, c.det(B[2:5,2:5]))
        if not axial_det > 0 or not swirl_det > 0: raise AssertionError('Exact determinant positivity lost')
        for i in (0, 2):
            direct = integrate(lambda x: c.expm1(-mu*c.log(x))/mu*bump(i,x))
            near(c, value(W['D'][0 if i == 0 else 1]), direct); weights_checked += 1
        for group, power in (('cross', '.5'), ('energy', '0'), ('pressure', '-1')):
            indices = (0, 2) if group == 'cross' else (0, 1, 2)
            for j, i in enumerate(indices):
                direct = integrate(lambda x: x**c.mpf(power)*bump(i,x)**2)
                near(c, value(W[group][j]), direct); weights_checked += 1
        rhs = list(map(c.mpf, ('-.23', '.17', '.09', '-.08', '.11')))
        roots = current.inverse_action(g, W, [g.constant(str(v)) for v in rhs])
        independent = c.lu_solve(B, c.matrix(rhs))
        for root, expected in zip(roots, independent): near(c, value(root), expected); inverse_checked += 1
        iv=MPIntervalContext();iv.dps=80
        mu_box=iv.mpf([str(mu*(1-c.mpf('1e-12'))),str(mu*(1+c.mpf('1e-12')))])
        ranged_weights=repair.fresh_weights(iv,mu_box,cells=512)
        ranged_matrix=repair.fresh_linear_inverse(iv,mu_box,ranged_weights)
        basis=tuple(iv.mpf(0) for unused in range(5));ledger=dict(directed_small_exponential_tails=0,
            positive_function_denominator_intersections=0,positive_function_root_intersections=0,directed_independent_log_rescalings=0)
        def box(v):
            return signed.ScaledEnclosure(signed.FormalScale(basis),
                iv.mpf([str(v-c.mpf('1e-12')),str(v+c.mpf('1e-12'))]),ledger)
        ranged_h=[current.C1Enclosure(box(v),box(j)) for v,j in zip(h,hZ)]
        ranged_Q=current.range_quadratic(iv,mu_box,iv.ln(mu_box),ranged_matrix,ranged_h)
        for N in (160, 2048):
            F = lambda x: sum(h[i+2]*bump(i,x) for i in range(3))/N
            G = lambda x: (h[0]*bump(0,x)+h[1]*bump(2,x))/N
            FZ = lambda x: sum(hZ[i+2]*bump(i,x) for i in range(3))/N
            GZ = lambda x: (hZ[0]*bump(0,x)+hZ[1]*bump(2,x))/N
            power = lambda x: x**(-c.mpf('.5')-mu)
            # Use the accepted original physical increment definitions at
            # E0=power,V0=0,dE=F,dV=G, and apply each cumulative rate once.
            density = lambda x: repair.packets.recovery.increment_densities(power(x), c.mpf(0), F(x), G(x))
            M = integrate(lambda x: density(x)['m'])
            J = integrate(lambda x: c.sqrt(x)*density(x)['k'])
            moments = [M, (J-M)/mu, integrate(lambda x:c.sqrt(x)*density(x)['h']),
                integrate(lambda x:density(x)['e']), integrate(lambda x:density(x)['p']/x)]
            MZ = integrate(GZ)
            JZ = integrate(lambda x:c.sqrt(x)*((power(x)+F(x))*GZ(x)+FZ(x)*G(x)))
            jets = [MZ, (JZ-MZ)/mu, integrate(lambda x:c.sqrt(x)*FZ(x)),
                integrate(lambda x:2*G(x)*GZ(x)-power(x)*FZ(x)-F(x)*FZ(x)),
                integrate(lambda x:(power(x)*FZ(x)+F(x)*FZ(x))/x)]
            hr = [g.constant(str(q)) for q in h]; jr = [g.constant(str(q)) for q in hZ]
            Q = [value(q) for q in current.quadratic(g,W,hr)]
            QZ = [value(q) for q in current.quadratic_action(g,W,hr,jr)]
            linear, linearZ = B*c.matrix(h), B*c.matrix(hZ)
            for i in range(5):
                near(c, linear[i]+Q[i]/N, N*moments[i]); checks += 1
                near(c, linearZ[i]+QZ[i]/N, N*jets[i]); checks += 1
                for enclosed, physical in ((ranged_Q[i].value,N*(N*moments[i]-linear[i])),
                    (ranged_Q[i].Z,N*(N*jets[i]-linearZ[i]))):
                    if i in (0,2):
                        near(c,physical,c.mpf(0)); continue
                    lo,hi=ep(enclosed.finite_interval())
                    if not lo <= physical <= hi: raise AssertionError('Physical quadratic value outside directed action range')
                    range_checked += 1
    return dict(passed=True, direct_physical_moment_and_Z_comparisons=checks,
        independent_x_coordinate_weight_integrals=weights_checked, independent_dense_inverse_comparisons=inverse_checked,
        direct_physical_quadratic_values_inside_fresh_C1_ranges=range_checked,
        positive_determinant_comparisons=6, reference_mu_values=['.04','1e-40','.166666666666'], reference_N_values=[160,2048],
        fresh_raw_bump_quadrature_not_saved_matrix_midpoints=True, cross_integrates_sqrt_x_not_inverse_sqrt_x=True,
        manufactured_operator_coefficients_only=True, current_original_field_or_terminal_closure_evaluated=False)


def physical_amplitude_normalization():
    """Original physical A/A_Z units cancel only with the full quotient rule."""
    Z = sy.Symbol('Z',real=True); x, mu = sy.symbols('x mu',positive=True)
    A,F,G = [sy.Function(name)(Z) for name in ('A','F','G')]
    power = x**(-sy.Rational(1,2)-mu)
    original = repair.packets.recovery.increment_densities(A*power,sy.Integer(0),A*F,A*G)
    normalized = dict(m=G,h=sy.sqrt(x)*F,k=sy.sqrt(x)*(power+F)*G,
        e=G*G-power*F-F*F/2,p=(power*F+F*F/2)/x)
    count = 0
    for key in normalized:
        degree = 1 if key in ('m','h') else 2
        rate = {'m':1,'h':sy.Rational(3,2),'k':sy.Rational(3,2),'e':1,'p':0}[key]
        raw = x**(rate-1)*original[key]
        if sy.simplify(raw/A**degree-normalized[key]) != 0: raise AssertionError('Physical normalization differs')
        jet = sy.diff(raw,Z)/A**degree-degree*sy.diff(A,Z)*raw/A**(degree+1)
        if sy.simplify(jet-sy.diff(normalized[key],Z)) != 0: raise AssertionError('A_Z unit rule lost')
        count += 2
    return dict(passed=True, independent_physical_A_and_A_Z_identities=count,
        original_increment_density_callable=repair.packets.recovery.__name__+'.increment_densities')


def symbolic_function_map(built):
    """Differentiate each exact local Picard step, cutting prior iterates as functions."""
    g = built['graph']; Z = sy.Symbol('Z',real=True); checks = 0
    d_atoms = {}
    for key, pair in built['N_scaled_targets'].items():
        f = sy.Function('actual_N_scaled_'+key.replace('=','_').replace('/','_'))(Z)
        d_atoms[pair.value.node], d_atoms[pair.Z.node] = f, sy.diff(f,Z)
    def walker(atoms):
        cache = {}
        def walk(i):
            if i in atoms: return atoms[i]
            if i in cache: return cache[i]
            row=g.nodes[i]; op=row['operation']
            if op=='exact_rational': q=sy.Rational(row['numerator'],row['denominator'])
            elif op=='original_source_parameter': q=sy.Symbol('parameter_'+row['name'],positive=True)
            elif op=='shared_positive_integer': q=sy.Symbol('N',positive=True,integer=True)
            elif op=='sum': q=sum(walk(j) for j in row['arguments'])
            elif op=='product': q=sy.prod(walk(j) for j in row['arguments'])
            elif op=='negative': q=-walk(row['argument'])
            elif op=='positive_quotient': q=walk(row['numerator'])/walk(row['denominator'])
            elif op=='definite_integral': q=sy.Symbol('exact_Z_independent_weight_integral_'+str(i),positive=True)
            elif op=='analytic_unary':
                arg=walk(row['argument'])
                q=sy.Function('exprel')(arg) if row['name']=='exprel' else getattr(sy,row['name'])(arg)
            elif op=='bound_variable': q=sy.Symbol(row['name'],positive=True)
            elif op=='compact_raw_beta': q=sy.Function('raw_beta')(walk(row['argument']))
            else: raise AssertionError('Unexpected uncut original source in repair action: '+op)
            cache[i]=q; return q
        return walk
    for step in range(1,len(built['finite_picard_sequence'])):
        atoms = dict(d_atoms)
        for i,pair in enumerate(built['finite_picard_sequence'][step-1]):
            if pair.value == g.zero: continue
            f=sy.Function('prior_control_'+str(i))(Z)
            atoms[pair.value.node],atoms[pair.Z.node]=f,sy.diff(f,Z)
        walk=walker(atoms)
        for pair in built['finite_picard_sequence'][step]:
            if sy.simplify(sy.diff(walk(pair.value.node),Z)-walk(pair.Z.node)) != 0:
                raise AssertionError('Finite Picard derivative differs')
            checks += 1
    atoms=dict(d_atoms)
    final=built['finite_picard_sequence'][-1]
    for i,pair in enumerate(final):
        f=sy.Function('final_control_'+str(i))(Z); atoms[pair.value.node],atoms[pair.Z.node]=f,sy.diff(f,Z)
    A=sy.Function('actual_A')(Z)
    atoms[built['amplitude'].value.node],atoms[built['amplitude'].Z.node]=A,sy.diff(A,Z)
    walk=walker(atoms)
    for pair in [*built['control_residual'],*built['band_profiles'].values()]:
        if sy.simplify(sy.diff(walk(pair.value.node),Z)-walk(pair.Z.node)) != 0:
            raise AssertionError('Residual/band profile derivative differs')
        checks += 1
    h_symbols=sy.symbols('a0 a2 e0 e1 e2',real=True)
    matrix_atoms=dict(d_atoms)
    for pair,h_symbol in zip(final,h_symbols): matrix_atoms[pair.value.node]=h_symbol
    matrix_walk=walker(matrix_atoms)
    Q=current.quadratic(g,built['repair_weights'],[pair.value for pair in final])
    N=matrix_walk(built['N'].node);jacobian_checks=0
    for i in range(5):
        for j in range(5):
            expected=matrix_walk(built['repair_weights']['B'][i][j].node)+sy.diff(matrix_walk(Q[i].node),h_symbols[j])/N
            if sy.simplify(matrix_walk(built['implicit_Z_matrix'][i][j].node)-expected) != 0:
                raise AssertionError('Implicit Z matrix column/normalization differs')
            jacobian_checks += 1
    return dict(passed=True, exact_Picard_residual_and_band_C1_derivatives=checks,
        independent_symbolic_implicit_Jacobian_entries=jacobian_checks,
        genuine_N_scaled_target_functions_kept=True, finite_depth_does_not_prove_convergence=True)


def live_source_bindings(owner,live,range_owner):
    built=live['built']; g=built['graph']; original=live['original_source_ranges']
    if len(built['source_roles']) != 212 or len(built['cells']) != 24: raise AssertionError('Source family route changed')
    if len(live['ranges']) != 4: raise AssertionError('Expected three live Picard iterations')
    for pair in live['ranges'][0]:
        if not pair.value.zero or not pair.Z.zero: raise AssertionError('Nonzero initial control iterate')
    for key,pair in zip(current.ROWS,live['N_scaled_target_ranges']):
        if pair.value.scale.bases is not original['target']['values'][key].scale.bases:
            raise AssertionError('Original fixed factor basis changed')
        if pair.value.ledger is not original['target']['values'][key].ledger:
            raise AssertionError('Directed arithmetic ledger changed')
    N_nodes=[i for i,row in enumerate(g.nodes) if row['operation']=='shared_positive_integer']
    if N_nodes != [built['N'].node]: raise AssertionError('Independent repair frequency admitted')
    diagnostic=current.fixed_N_contraction_diagnostic(live)
    if diagnostic['half_contraction_sufficient_bound_passed']:
        raise AssertionError('Current enormous fixed-N bound unexpectedly certifies contraction')
    failed=0
    for args in ((original['record']['Z_box'],live['candidate_N']+1), ((0,0),live['candidate_N'])):
        try: owner.controls(*args,range_owner=range_owner,source_ranges=original,iterations=3)
        except ValueError: failed += 1
    if failed != 2: raise AssertionError('Range domain/frequency guard missing')
    try: current.FunctionEvaluator(built,oracle=None,Z=0,N=2048,ctx=mp.mp)
    except source.SourceOracleRequired: failed += 1
    else: raise AssertionError('Missing numerical source oracle accepted')
    try: g.mul(g.one,{'log_absolute_upper':'0'})
    except TypeError: failed += 1
    else: raise AssertionError('Cap accepted as function coefficient')
    # Test mu=0 rejection on the isolated exact weight graph only.
    c=mp.mp.clone();c.dps=40;unused,W,value=scalar_fixture(owner,'0',c)
    try: value(W['D'][0])
    except ArithmeticError: failed += 1
    else: raise AssertionError('Generic divided repair silently continued to mu=0')
    return dict(passed=True, original_signed_density_and_amplitude_roles=212,
        original_continuous_radial_cells=24, live_C0_and_Z_control_ranges=3*5*2,
        original_target_domain=live['Z_box'], candidate_N=live['candidate_N'], fail_closed_guards=failed,
        original_context_basis_ledger_and_actual_mu_preserved=True,
        fixed_N_sufficient_contraction_test_passed=False, actual_fixed_point_or_terminal_closure_admitted=False)


@current.paired.native.inlet.source_precision
def run(owner,live,*,range_owner):
    began=time.monotonic()
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
            independent_physical_integrals=direct_physical_reference(owner),
            original_physical_amplitude_normalization=physical_amplitude_normalization(),
            exact_C1_function_rules=symbolic_function_map(live['built']),
            actual_source_and_live_iterate_binding=live_source_bindings(owner,live,range_owner),
            read_only_review_model=dict(model='gpt-5.6-luna',reasoning_effort='max'),
            current_original_source_points_or_actual_terminal_integrals_evaluated=False,
            certified_fixed_point_tail_installed=False,actual_five_controls_installed=False,
            actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
            **dict.fromkeys(current.packets.OPEN,False),
            input_hashes={**owner.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name)},
            scope='Independent original physical density/x-integral and derivative normalization, exact bump/inverse/Picard C1 identities and actual full-source finite range provenance. Manufactured operator references do not establish current original-field closure.')
    result['execution_seconds']=time.monotonic()-began
    (HERE/current.RECEIPT).write_text(json.dumps(current.packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Source-bound control-map independent physical/C1 check PASS; fixed point remains open',flush=True)
    return result
