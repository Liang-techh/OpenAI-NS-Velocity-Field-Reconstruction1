"""Sharp same-function magnitude coordinates and collected Poisson Z jets.

Retain the original source, cutoff, inverse, signed densities and fixed N.
Only interval evaluation coordinates and two dependent factors change.
"""
import ast
import gzip
import json
from pathlib import Path
import time
from types import FunctionType
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_correlated_functions as previous
import lei_ren_part1_paper_compliant_flat_pulse_derivatives as cutoff_derivatives

graph,bridge,backend = previous.graph,previous.bridge,previous.backend
HERE,PREFIX,sha,bind,ep = previous.HERE,previous.PREFIX,previous.sha,previous.bind,previous.ep
CELLS,RATES,C0,Z,OPEN = previous.CELLS,previous.RATES,previous.C0,previous.Z,previous.OPEN
serialized,encode = previous.serialized,previous.encode
NAME = PREFIX+'current_original_whole_Z_bridge_sharp_functions.json.gz'
RECEIPT = PREFIX+'current_original_whole_Z_bridge_sharp_functions_check.json'
GATE = 'current_original_whole_Z_sharp_magnitude_and_collected_Poisson_Z_bridge_installed'


def original_q_squared_Z(source,qrows):
    """Differentiate the original q squared, collecting its first eta ratio.

    The original q C0 and linear q_Z remain untouched. A mixed Delta box
    uses the union of negative, transition and flat derivative covers.
    It never applies an active root or positive gamma to the entire box.
    """
    q = qrows[C0]; c = q.ctx
    if source['q'] is not q:raise ValueError('Same original q C0 object required')
    if q.zero:
        if not qrows[Z].zero:raise ValueError('Original exact flat q requires zero q_Z')
        return {C0:q,Z:q},dict(exact_original_flat_or_source_collar=True)
    proof = source['original_full_source_quotients']
    a = source['roots']['a']
    if proof['actual_a_axial5'][0] is not a[C0] or proof['actual_a_axial5'][1] is not a[Z]:
        raise ValueError('Same original a C0/Z objects required')
    Delta,dz = proof['actual_Delta_axial5'][:2]
    for row in (qrows[Z],a[C0],a[Z],Delta,dz):q.coerce(row)
    loga = proof['actual_positive_a']['source_log_lower']
    if ep(a[C0].coefficient)[0] <= 0:raise ValueError('Actual a positivity required')
    square = previous.switch.parameters.original.square(q)
    sigma_cap = cutoff_derivatives.sigma_tail_bound(c,1,c.mpf('.5'))
    # g'= -1 on Delta<=0, -sigma^2-2sigma*sigma'*(2-Delta/eta)
    # on [0,eta], and 0 on [eta,infinity). The independently checked
    # original global sigma cap (about 32) bounds the transition term.
    negative = ep(Delta.coefficient)[1] <= 0
    bound = c.mpf(1) if negative else 1+4*sigma_cap
    derivative = q.scalar(-1 if negative else c.mpf([-ep(bound)[1],0]))
    square_Z = (derivative*dz*c.mpf('.5')-square*a[Z]).positive_divide(a[C0],loga)
    return {C0:square,Z:square_Z},dict(
        original_q_C0_object_preserved=True,original_linear_q_Z_preserved=True,
        original_q_squared_C0_from_nonnegative_square=True,
        source_formula='q^2=g(Delta)/(2a)',
        derivative_formula='g_prime*Delta_Z/(2a)-q^2*a_Z/a',
        g_prime_signed_union=derivative.record(),
        sigma_prime_original_checked_global_bound=sigma_cap,
        same_source_negative_transition_flat_union_not_field_selection=not negative,
        first_transition_eta_over_eta_cancelled_before_interval_evaluation=True,
        no_y_derivative_rows_constructed=True,derivative_of_cap_not_used=True)


def bound_scale_intersection(value,bound):
    """Intersect signed endpoints in the theorem-bound scale.

    Ratio logs retain the shared formal bases. No large positive exponential
    is formed, and a negative tail remains in the bounding source's units.
    This returns an enclosure of the same evaluated function, not a value
    chosen from the theorem's magnitude cap.
    """
    value.coerce(bound)
    if value.zero:return value
    c = value.ctx
    radius = max(abs(v) for v in ep(bound.coefficient))
    if not radius:raise ArithmeticError('Nonzero function incompatible with exact zero magnitude bound')
    r = c.mpf(radius)
    logR = c.ln(r)
    dlo,dhi = ep((value.scale-bound.scale).evaluate())
    lo,hi = ep(value.coefficient)
    def endpoint(coefficient,side):
        if not coefficient:return mp.mpf(0)
        positive = coefficient > 0
        exponent = dlo if (side == 'lower') == positive else dhi
        logratio = c.mpf(exponent)+c.ln(c.mpf(abs(coefficient)))-logR
        l,u = ep(logratio)
        # Intersection clips only the outward endpoint. Inward endpoints
        # beyond the opposite edge certify disjoint source/theorem covers.
        if (side == 'lower' and positive or side == 'upper' and not positive) and l > 0:
            raise ArithmeticError('Signed source interval and same-function magnitude theorem are disjoint')
        chosen = l if (side == 'lower') == positive else u
        chosen = min(chosen,mp.mpf(0))
        mag = r*value.bounded_exp(c.mpf(chosen))
        signed = mag if positive else -mag
        return ep(signed)[0 if side == 'lower' else 1]
    lower = max(endpoint(lo,'lower'),-radius)
    upper = min(endpoint(hi,'upper'),radius)
    if lower > upper:raise ArithmeticError('Empty signed magnitude intersection')
    prior = backend.signed.density.local.prior
    return prior.ScaledEnclosure(bound.scale,c.mpf((lower,upper)),value.ledger)


GEOMETRY_PATTERNS = (
    ('-loop.u * loop.hinv * loop.hinv * loop.hinv * u.d[key]', '-loop.uhinvcube * u.d[key]'),
    ('loop.hinv * loop.hinv * loop.hinv * u.d[key]', 'loop.hinvcube * u.d[key]'))


def collect_geometry_assignments(target):
    patterns = [ast.dump(ast.parse(before,mode='eval').body) for before,_ in GEOMETRY_PATTERNS]
    counts = [0,0]
    class Collect(ast.NodeTransformer):
        def visit(self,node):
            for i,wanted in enumerate(patterns):
                if ast.dump(node) == wanted:
                    counts[i] += 1
                    return ast.copy_location(ast.parse(GEOMETRY_PATTERNS[i][1],mode='eval').body,node)
            return super().visit(node)
    Collect().visit(target)
    if counts != [1,1]:raise ValueError('Exactly two original dependent Poisson Z products required')
    quadratic = ast.dump(ast.parse('q.square()',mode='eval').body)
    qcount = 0
    class Quadratic(ast.NodeTransformer):
        def visit_Call(self,node):
            nonlocal qcount
            if ast.dump(node) == quadratic:
                qcount += 1
                return ast.copy_location(ast.Name(id='original_q2_dual',ctx=ast.Load()),node)
            return self.generic_visit(node)
    Quadratic().visit(target)
    if qcount != 5:raise ValueError('Exactly five original quadratic q sites required')
    indices = [i for i,n in enumerate(target.body) if isinstance(n,ast.Assign)
        and ast.unparse(n.targets[0]) == 'nu']
    if len(indices) != 1:raise ValueError('Original phase nu injection site required')
    target.body[indices[0]:indices[0]] = ast.parse(
        "q2rows=source['original_q_squared_C0_Z']\n"
        "original_q2_dual=source_dual(q2rows[ZERO],q2rows)").body
    return target


COMPILE_INVERSE,BINDER_REWRITE = previous.replacement_function(previous.compile_inverse_identity_Z,
    [('changed = ast.dump(target)', 'collect_geometry_assignments(target)\nchanged = ast.dump(target)')],
    'compile_collected_inverse_identity_Z')
COMPILE_INVERSE = FunctionType(COMPILE_INVERSE.__code__,
    {**COMPILE_INVERSE.__globals__,'collect_geometry_assignments':collect_geometry_assignments},COMPILE_INVERSE.__name__)
SHARP_Z_FIRST,INVERSE_BINDING = COMPILE_INVERSE()
INVERSE_BINDING.update(collected_geometry_replacements=[dict(before=a,after=b) for a,b in GEOMETRY_PATTERNS],
    exact_collected_geometry_replacement_counts=[1,1],
    exact_original_q_square_replacement_count=5,
    direct_q_squared_Z_from_same_original_cutoff_source=True,
    original_root_q_p2_derivatives_and_implicit_inverse_rules_unchanged=True,
    source_bound_uhinvcube_and_hinvcube_not_independent_products=True)

SHARP_BRANCH_BODY = FunctionType(previous.correlated_inverse_branch.__code__,
    {**previous.correlated_inverse_branch.__globals__,'CORRELATED_Z_FIRST':SHARP_Z_FIRST,
     'INVERSE_BINDING':INVERSE_BINDING,'intersect_same_function_magnitude':bound_scale_intersection},
    'sharp_same_source_inverse_Z',previous.correlated_inverse_branch.__defaults__)


def sharp_inverse_branch(source,qrows,dstar_log,phi,branch):
    if 'original_q_squared_C0_Z' not in source:
        q2rows,q2proof = original_q_squared_Z(source,qrows)
        source = dict(source,original_q_squared_C0_Z=q2rows,original_q_squared_Z_proof=q2proof)
    got = SHARP_BRANCH_BODY(source,qrows,dstar_log,phi,branch)
    got['record']['original_direct_q_squared_Z_cover'] = source['original_q_squared_Z_proof']
    if got['values'] is not None and got['record']['geometry'] != 'flat':
        got['record']['same_function_active_flat_magnitude_theorem'].update(
            intersection_coordinate='actual bounding source formal scale',
            signed_endpoints_clipped_by_log_comparison_before_exponentiation=True,
            tiny_tails_not_moved_to_a_dominating_unrelated_scale=True)
        got['record']['collected_original_Poisson_geometry_Z_factors'] = dict(
            u_hinverse_cubed=got['loop'].uhinvcube.record(),
            hinverse_cubed=got['loop'].hinvcube.record(),
            original_q_p2_Z_chain_retained=True)
    return got


QUERY_BODY,QUERY_BINDING = previous.replacement_function(backend.WholeZSignedLoopFunctions.query,
    [("source=dict(q=qrows[C0],roots=roots)",
      "source=dict(q=qrows[C0],roots=roots,original_C_B_E_C0_Z=packet['original_C_B_E_C0_Z'],"
      "original_eta_log=packet['original_eta_log'],original_full_source_quotients=packet['original_full_source_quotients'])")],
    'sharp_same_source_bridge_query')
SHARP_QUERY = FunctionType(QUERY_BODY.__code__,
    {**bridge.BRIDGE_QUERY.__globals__,'branch_Z_functions':sharp_inverse_branch},
    QUERY_BODY.__name__,QUERY_BODY.__defaults__)


class WholeZSharpBridgeFunctions(previous.WholeZCorrelatedBridgeFunctions):
    def __init__(self,dps=500):
        super().__init__(dps)
        receipt = json.loads((HERE/previous.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(previous.GATE) \
                or receipt['source_family'] != self.identity or receipt['candidate_N'] != self.N:
            raise ValueError('Accepted actual correlated bridge arithmetic required')
        for name,digest in receipt['input_hashes'].items():bind(self.hashes,name,digest)
        bind(self.hashes,previous.RECEIPT,sha(previous.RECEIPT))
        for module in (previous,backend.signed):
            name = Path(module.__file__).name
            bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        flat = PREFIX+'flat_pulse_derivatives'
        checked = json.loads((HERE/(flat+'_check.json')).read_bytes())
        if not checked.get('all_passed') or not checked.get('original_radial_shape_derivatives_C4_available') \
                or checked['actual_five_defect_family_sha256'] != self.identity['actual_five_defect_family_sha256']:
            raise ValueError('Same-family checked original sigma derivative theorem required')
        for name,digest in checked['input_hashes'].items():
            if sha(name) != digest:raise ValueError('Changed cutoff derivative prerequisite: '+name)
            bind(self.hashes,name,digest)
        for name in (flat+'.py',flat+'.json',flat+'_check.json'):
            bind(self.hashes,name,sha(name))
        original = PREFIX+'current_generic_shear_loop.py'
        tree = ast.parse((HERE/original).read_text(encoding='utf8'))
        step = next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name == 'flat_step')
        odds = next(n.value for n in ast.walk(step) if isinstance(n,ast.Assign)
            and any(isinstance(t,ast.Name) and t.id == 'odds' for t in n.targets))
        if ast.dump(odds) != ast.dump(ast.parse('1/(1-x)**2-1/x**2',mode='eval').body):
            raise ValueError('Same original flat cutoff odds AST required')
        bind(self.hashes,original,sha(original))

    def query(self,ends,chart,left,right=None):
        if chart not in bridge.CHARTS:return super().query(ends,chart,left,right)
        l = self.coordinate(chart,left)
        r = self.coordinate(chart,left if right is None else right)
        value = SHARP_QUERY(self,ends,chart,l,r)
        value['actual_source_bounded_exponent_callback_binding'] = bridge.DENSITY_BINDING
        value['correlated_inverse_query_binding'] = QUERY_BINDING
        value['actual_sharp_magnitude_and_collected_geometry_binding'] = INVERSE_BINDING
        return value


def run():
    began = time.monotonic()
    with mp.workdps(540):
        owner = WholeZSharpBridgeFunctions()
        rows = []
        for ends in CELLS:
            rows.append(serialized(owner.retransport(ends)))
            print('Sharp bridge and actual Rc retransport: '+str(ends),flush=True)
        report = dict(**{GATE:True},source_family=owner.identity,candidate_N=owner.N,
            exact_Z_partition=CELLS,actual_sharp_bridge_and_Rc_cells=rows,
            inverse_binding=INVERSE_BINDING,compiler_injection_binding=BINDER_REWRITE,query_binding=QUERY_BINDING,
            actual_terminal_controls_installed=False,actual_global_frequency_admitted=False,
            physical_original_exterior_five_targets_closed=False,**dict.fromkeys(OPEN,False),
            input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Original actual bridge functions with signed logarithmic intersection in bound units '
                'two analytically collected Poisson Z geometry factors and direct original q squared Z cover; live bridge integration and '
                'exact unchanged downstream local-driver retransport through Rc. Original q, cutoff, '
                'phase, widths, P0 and N remain. No controls, global N, exterior/heat/cone or recursion admission.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(serialized(report)),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__ == '__main__':run()
