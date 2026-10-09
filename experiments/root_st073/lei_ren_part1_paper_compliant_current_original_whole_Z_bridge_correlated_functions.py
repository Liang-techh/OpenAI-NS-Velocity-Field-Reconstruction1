"""Same bridge functions evaluated through correlated source identities.

The defining phase, cutoff, densities, radial measures and fixed N remain
unchanged. C/B/E quotients and inverse primitives are rewritten before
interval evaluation; derivatives are of these identities, not of caps.
"""
import ast
import copy
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
import time
from types import FunctionType, SimpleNamespace
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rc_relative_defects as relative

graph = relative.graph
bridge, backend = graph.previous, graph.backend
switch, phase = bridge.switch, backend.phase
HERE, PREFIX, sha, bind, ep = graph.HERE, graph.PREFIX, graph.sha, graph.bind, graph.ep
CELLS, RATES, C0, Z, OPEN = graph.CELLS, graph.RATES, graph.C0, graph.Z, graph.OPEN
serialized, encode = graph.serialized, graph.encode
NAME = PREFIX+'current_original_whole_Z_bridge_correlated_functions.json.gz'
RECEIPT = PREFIX+'current_original_whole_Z_bridge_correlated_functions_check.json'
GATE = 'current_original_whole_Z_correlated_bridge_C1_functions_and_Rc_retransport_installed'


def replacement_function(fn, replacements, name):
    """Replace exactly declared statements, with reversible AST evidence."""
    node = ast.parse(textwrap.dedent(inspect.getsource(fn))).body[0]
    original = ast.dump(node)
    edits = []
    for before, after in replacements:
        wanted = ast.dump(ast.parse(before).body[0])
        hits = [i for i, row in enumerate(node.body) if ast.dump(row) == wanted]
        if len(hits) != 1:
            raise ValueError('Exactly one defining statement required: '+before)
        index = hits[0]
        old = node.body[index]
        new = ast.parse(after).body
        node.body[index:index+1] = new
        edits.append((index, old, new, before, after))
    changed = ast.dump(node)
    for index, old, new, _, _ in reversed(edits):
        node.body[index:index+len(new)] = [old]
    if ast.dump(node) != original:
        raise ValueError('Undeclared defining function change')
    for index, _, new, _, _ in edits:
        node.body[index:index+1] = new
    node.name = name
    env = dict(fn.__globals__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[node], type_ignores=[])),
        '<same function; correlated source identities>', 'exec'), env)
    return env[name], dict(original_module=Path(inspect.getfile(fn)).name,
        original_module_sha256=sha(Path(inspect.getfile(fn)).name),
        original_AST_sha256=hashlib.sha256(original.encode()).hexdigest(),
        rewritten_AST_sha256=hashlib.sha256(changed.encode()).hexdigest(),
        reversible_declared_statement_rewrites_only=True,
        exact_replacements=[dict(before=e[3], after=e[4]) for e in edits])


CORRELATED_QUOTIENTS, QUOTIENT_BINDING = replacement_function(switch.general_quotients, [
    ("t0=parameters.quotient(f,parameters.scale(b,-1),a,pa)",
     "C=n['C']; pc=positive_jet(f,C,'actual_same_source_C'); "
     "t0=parameters.quotient(f,parameters.scale(n['B'],-1),C,pc)"),
    ("kappa=parameters.add(f,a,parameters.quotient(f,square,a,pa))",
     "C2=parameters.multiply(f,C,C); C2[0]=parameters.original.square(C[0]); "
     "B2=parameters.multiply(f,n['B'],n['B']); B2[0]=parameters.original.square(n['B'][0]); "
     "CE=parameters.multiply(f,C,E); pce=positive_jet(f,CE,'actual_same_source_CE'); "
     "joint=parameters.add(f,C2,B2); kappa=parameters.quotient(f,joint,CE,pce)"),
    ("Delta=parameters.add(f,kappa,[f.scalar(-2)]+[f.scalar(0)]*5)",
     "jointDelta=parameters.add(f,joint,parameters.scale(CE,-2)); "
     "Delta=parameters.quotient(f,jointDelta,CE,pce)")], 'correlated_general_quotients')


def compile_inverse_identity_Z():
    """Keep F/Fx and all tails; rewrite A/B before first differentiation."""
    tree = ast.parse(Path(phase.first.__file__).read_text(encoding='utf8'))
    target = next(n for n in tree.body if getattr(n,'name',None) == 'conditioned_first_jets')
    original = ast.dump(target)
    edits = []
    rules = {
        "A=(a*numerator).divide(nu,0)*c.mpf('.5')":
            "A=a*(First(loop.scalar(c.mpf(phi)))-psi)*c.mpf('.5')",
        "B=E*(t0*A-divr(a*q*hinv*difference)*c.mpf('.5'))":
            "B=-Bsource*(First(loop.scalar(c.mpf(phi)))-psi)*c.mpf('.5')-divr(Csource*q*hinv*difference)*c.mpf('.5')",
        "A=(a*numerator).divide(nu,0)*(1/(4*c.pi))":
            "A=a*(First(loop.scalar(c.mpf(phi)))-psi)*c.mpf('.5')",
        "B=E*(t0*A-a*q*hinv*W1*(1/(2*c.pi)))":
            "B=-Bsource*(First(loop.scalar(c.mpf(phi)))-psi)*c.mpf('.5')-Csource*q*hinv*W1*(1/(2*c.pi))"}
    mapping = {ast.dump(ast.parse(k).body[0]):(k,v) for k,v in rules.items()}
    class Rewrite(ast.NodeTransformer):
        def visit_Assign(self, node):
            key = ast.dump(node)
            if key in mapping:
                before, after = mapping[key]
                edits.append((before,after))
                return ast.copy_location(ast.parse(after).body[0], node)
            return self.generic_visit(node)
    Rewrite().visit(target)
    if set(before for before,_ in edits) != set(rules) or len(edits) != 4:
        raise ValueError('Original Mobius and small-r primitive assignments changed')
    variable = "variables={name:source_dual(rows[ZERO],rows) for name,rows in source['roots'].items()}"
    idx = [i for i,n in enumerate(target.body) if ast.dump(n) == ast.dump(ast.parse(variable).body[0])]
    if len(idx) != 1:
        raise ValueError('Same source-dual injection site required')
    target.body[idx[0]+1:idx[0]+1] = ast.parse(
        "numerators=source['original_C_B_E_C0_Z']\n"
        "Csource=source_dual(numerators['C'][ZERO],numerators['C'])\n"
        "Bsource=source_dual(numerators['B'][ZERO],numerators['B'])").body
    idx = [i for i,n in enumerate(target.body) if isinstance(n,ast.Assign)
        and ast.unparse(n.targets[0]) == '(philo, phihi)']
    if len(idx) != 1:
        raise ValueError('Explicit phi derivative must precede unchanged symmetry handling')
    target.body[idx[0]:idx[0]] = ast.parse(
        "outputs['A_phi'] += a.v*c.mpf('.5')\n"
        "outputs['B_phi_over_Pstar'] -= Bsource.v*c.mpf('.5')").body
    changed = ast.dump(target)
    # Use the existing exact six-site restriction, with the altered tree as
    # its only input. All First, source_dual and tail bodies remain original.
    real_ast = phase.compile_Z_only.__globals__['ast']
    proxy_ast = SimpleNamespace(**{**vars(real_ast), 'parse': lambda *args,**kwargs:copy.deepcopy(tree)})
    compiler = FunctionType(phase.compile_Z_only.__code__,
        {**phase.compile_Z_only.__globals__, 'ast':proxy_ast}, 'compile_correlated_Z_only')
    fn, axes = compiler()
    return fn, dict(original_module=Path(phase.first.__file__).name,
        original_module_sha256=sha(Path(phase.first.__file__).name),
        original_first_dual_AST_sha256=hashlib.sha256(original.encode()).hexdigest(),
        rewritten_first_dual_AST_sha256=hashlib.sha256(changed.encode()).hexdigest(),
        exact_primitive_rewrites=[dict(before=a,after=b) for a,b in edits],
        original_F_Fx_inverse_selection_and_series_tails_unchanged=True,
        same_source_C_B_rows_injected=True,
        explicit_phi_terms=dict(A_phi='a/2',B_phi_over_Pstar='-Bsource/2'),
        actual_A_Z_and_B_Z_from_identity_not_retained_expansion=True,
        exact_Z_only_axis_restrictions=axes['exact_axis_restrictions'],
        derivative_of_support_cap_not_used=True)


CORRELATED_Z_FIRST, INVERSE_BINDING = compile_inverse_identity_Z()


def intersect_same_function_magnitude(value, bound):
    """Intersect two proved enclosures, never select the bounding function."""
    if value.zero:
        return value
    value.coerce(bound)
    reference = graph.switch_backend.directed_dominating_union([value,bound]).scale
    c = value.ctx
    ordinary = lambda row:row.coefficient*row.bounded_exp((row.scale-reference).evaluate())
    lo,hi = ep(ordinary(value))
    radius = max(abs(v) for v in ep(ordinary(bound)))
    lower,upper = max(lo,-radius),min(hi,radius)
    if lower > upper:
        raise ArithmeticError('Same-function primitive enclosure and theorem are disjoint')
    prior = backend.signed.density.local.prior
    return prior.ScaledEnclosure(reference,c.mpf((lower,upper)),value.ledger)


def correlated_inverse_branch(source, qrows, dstar_log, phi, branch):
    namespace = dict(CORRELATED_Z_FIRST.__globals__)
    original_phase = namespace['phase']
    namespace['phase'] = SimpleNamespace(**{**vars(original_phase),
        'ConditionedPhase':lambda src,log:backend.signed.BranchConditionedPhase(src,log,branch)})
    fn = FunctionType(CORRELATED_Z_FIRST.__code__,namespace,
        'same_source_correlated_inverse_Z',CORRELATED_Z_FIRST.__defaults__)
    got = fn(source,qrows,dstar_log,phi)
    got['record'].update(correlated_inverse_identity_binding=INVERSE_BINDING,
        conditional_signed_u_geometry=branch['name'],
        source_roots_q_Z_and_original_geometry_not_selected=True)
    if got['values'] is not None and set(got['values']) != set(backend.OUTPUTS):
        raise ValueError('Only genuine Z/phi primitive slots may be exported')
    if got['values'] is not None and got['record']['geometry'] != 'flat':
        if ep(source['original_eta_log'])[1] > 0:
            raise ArithmeticError('Primitive magnitude theorem requires actual eta<=1')
        # On the active source Delta<eta, kappa=a*(1+t0^2), and
        # a*nu=kappa+sigma^2*(2*eta-Delta)<=2+2*eta. On the
        # complementary source the actual primitives are exactly flat.
        # Thus a<=3, |A|<=a/2<=1.5. Also |b_L|<=a*nu/2,
        # |b|<=kappa/2, so |(B/E)_phi|<=1+3*eta/4<=1.75.
        # Integrate over a period of length one from the zero origin.
        # These bounds apply to the same functions, including mixed cells.
        before = dict(got['values'])
        E = source['roots']['E'][C0]
        for key,cap in (('A',E.scalar(2)),('B_over_Pstar',E*2),
                ('A_phi',E.scalar(4)),('B_phi_over_Pstar',E*2)):
            got['values'][key] = intersect_same_function_magnitude(before[key],cap)
        got['record'].update(actual_identity_primitive_before_magnitude_intersection=
            {key:v.record() for key,v in before.items()},
            same_function_active_flat_magnitude_theorem=dict(actual_eta_log=source['original_eta_log'],
                active_domain='Delta<eta',flat_domain='Delta>=eta',
                a_nu_identity='kappa+sigma^2*(2eta-Delta)',a_nu_upper='2+2eta',
                A_magnitude_upper=2,B_over_E_magnitude_upper=2,
                A_phi_magnitude_upper=4,B_phi_over_E_magnitude_upper=2,
                actual_C0_phi_functions_intersected_not_replaced=True,
                Z_derivative_bounds_unchanged_by_magnitude_theorem=True,
                derivative_of_cap_not_used=True),
            original_A_B_first_derivative_enclosures={key:v.record() for key,v in got['values'].items()})
    return got


# Query keeps the same bounded exponent/expm1 density namespace; only the
# declared source bundle and inverse callback change.
QUERY_BODY, QUERY_BINDING = replacement_function(backend.WholeZSignedLoopFunctions.query,
    [("source=dict(q=qrows[C0],roots=roots)",
      "source=dict(q=qrows[C0],roots=roots,original_C_B_E_C0_Z=packet['original_C_B_E_C0_Z'],original_eta_log=packet['original_eta_log'])")],
    'correlated_bridge_query')
CORRELATED_QUERY = FunctionType(QUERY_BODY.__code__,
    {**bridge.BRIDGE_QUERY.__globals__, 'branch_Z_functions':correlated_inverse_branch},
    QUERY_BODY.__name__,QUERY_BODY.__defaults__)
SOURCE_QUERY = FunctionType(bridge.WholeZBridgeSignedFunctions.source_query.__code__,
    {**bridge.WholeZBridgeSignedFunctions.source_query.__globals__,
     'switch':SimpleNamespace(**{**vars(switch),'general_quotients':CORRELATED_QUOTIENTS})},
    'same_bridge_source_with_correlated_quotients',bridge.WholeZBridgeSignedFunctions.source_query.__defaults__,
    bridge.WholeZBridgeSignedFunctions.source_query.__closure__)


class WholeZCorrelatedBridgeFunctions(graph.WholeZFullSignedTransport):
    def __init__(self,dps=500):
        super().__init__(dps)
        receipt = json.loads((HERE/graph.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(graph.GATE) \
                or receipt['source_family'] != self.identity or receipt['candidate_N'] != self.N:
            raise ValueError('Accepted same-source full signed graph required')
        for name,digest in receipt['input_hashes'].items():
            bind(self.hashes,name,digest)
        bind(self.hashes,graph.RECEIPT,sha(graph.RECEIPT))
        self.accepted = json.loads(gzip.decompress((HERE/graph.NAME).read_bytes()))
        if self.accepted['source_family'] != self.identity or self.accepted['candidate_N'] != self.N:
            raise ValueError('Cached independent local functions have different source/N')
        self.cells = {tuple(row['exact_Z_cell']):row for row in self.accepted['actual_signed_full_transports']}
        for module in (relative,phase.first,switch):
            name = Path(module.__file__).name
            bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def source_query(self,ends,chart,left,right=None):
        if chart not in bridge.CHARTS:
            return super().source_query(ends,chart,left,right)
        packet = SOURCE_QUERY(self,ends,chart,left,right)
        n = packet['original_generic_source']['actual_generic_source_numerators']
        packet['original_C_B_E_C0_Z'] = {name:{C0:n[name][0],Z:n[name][1]} for name in ('C','B','E')}
        packet['original_eta_log'] = self.switch_source.eta_log
        packet['correlated_quotient_binding'] = QUOTIENT_BINDING
        return packet

    def query(self,ends,chart,left,right=None):
        if chart not in bridge.CHARTS:
            return super().query(ends,chart,left,right)
        l = self.coordinate(chart,left)
        r = self.coordinate(chart,left if right is None else right)
        value = CORRELATED_QUERY(self,ends,chart,l,r)
        value['actual_source_bounded_exponent_callback_binding'] = bridge.DENSITY_BINDING
        value['correlated_inverse_query_binding'] = QUERY_BINDING
        return value

    def retransport(self,ends):
        """Recompute bridge; hydrate only unchanged independent local drivers.

        Every new incoming history is propagated with its original memory.
        Old incoming/final intervals are never subtracted or reused as new
        histories. The live full_transport method remains available too.
        """
        key = tuple(ends)
        if key not in self.cells:
            raise ValueError('Admitted whole-Z cell required')
        f = self.owner(key).flow
        stored = self.cells[key]
        decode = lambda v:relative.restore_row(f,v)
        if any(not relative.exact_row(a,b) for a,b in zip(self.owner(key).P0,
                [decode(v) for v in stored['exact_common_P0_axial5']])):
            raise ValueError('Unchanged independent source P0 required')
        pairs = lambda rows:{name:[decode(v) for v in pair] for name,pair in rows.items()}
        fresh_bridge = self.bridge_transport(key)
        incoming = fresh_bridge['actual_R100_signed_correction_C0_Z']
        stages, trace = stored['actual_source_owned_zero_inlet_to_Rc_signed_stages'], []
        for stage,field in (('switch','actual_signed_switch_windows'),('long','actual_signed_long_windows'),
                ('patch',None),('outer','actual_signed_outer_windows')):
            windows = stages[stage][field] if field else [dict(actual_chart='actual_patch',
                original_incoming_own_rate_memory=stages['patch']['original_incoming_own_rate_memory'],
                actual_signed_local_driver_C0_Z=stages['patch']['actual_signed_local_driver_at_Rh_C0_Z'])]
            for old in windows:
                memory = {name:decode(v) for name,v in old['original_incoming_own_rate_memory'].items()}
                local = pairs(old['actual_signed_local_driver_C0_Z'])
                outgoing = {name:[memory[name]*incoming[name][i]+local[name][i] for i in range(2)] for name in RATES}
                trace.append(dict(stage=stage,actual_chart=old['actual_chart'],
                    actual_incoming_correction_C0_Z=incoming,original_incoming_own_rate_memory=memory,
                    unchanged_independent_local_driver_C0_Z=local,actual_exit_signed_correction_C0_Z=outgoing))
                incoming = outgoing
        outer = stages['outer']
        A = [decode(v) for v in outer['actual_positive_Rc_normalized_amplitude_C0_Z']]
        logmu = self.source.logmu
        target = relative.normalized_targets(f,incoming,A,f.factor((0,0,0,0,0),logmu),logmu,self.N)
        weights = relative.repair.fresh_weights(self.c,self.c.mpf([0,ep(self.c.mpf(1)/6)[1]]))
        matrix = relative.repair.fresh_linear_inverse(self.c,self.c.mpf([0,ep(self.c.mpf(1)/6)[1]]),weights)
        conditions = relative.repair.contraction_log_conditions(self.c,target['target_C1_caps'],matrix,weights,logmu,self.c.mpf(0))
        sufficient = ep(self.c.ln(self.N))[0] >= ep(conditions['repair_sufficient_common_log_N_lower'])[1]
        conditions.update(actual_fixed_log_N=self.c.ln(self.N),
            fixed_candidate_N_satisfies_sufficient_contraction_bound=sufficient,
            sufficient_test_failure_is_not_a_nonexistence_proof=True,
            fixed_N_cap_does_not_admit_a_changed_N=True)
        return dict(source_identity=self.identity,candidate_N=self.N,exact_Z_cell=list(key),
            exact_common_P0_axial5=self.owner(key).P0,actual_correlated_bridge_transport=fresh_bridge,
            actual_downstream_retransport_windows=trace,actual_Rc_signed_correction_C0_Z=incoming,
            **target,actual_fixed_N_repair_log_conditions=conditions,
            evaluation_mode='live_correlated_bridge_then_hash_bound_unchanged_local_driver_retransport',
            old_incoming_or_final_histories_not_reused=True,source_local_drivers_independent_of_correction=True,
            original_density_radial_measures_phase_N_and_pressure_memory_unchanged=True,
            actual_terminal_controls_installed=False,actual_global_frequency_admitted=False,
            physical_original_exterior_five_targets_closed=False,**dict.fromkeys(OPEN,False))


def run():
    began = time.monotonic()
    with mp.workdps(540):
        owner = WholeZCorrelatedBridgeFunctions()
        rows = []
        for ends in CELLS:
            rows.append(serialized(owner.retransport(ends)))
            print('Correlated bridge and actual Rc retransport: '+str(ends),flush=True)
        report = dict(**{GATE:True},source_family=owner.identity,candidate_N=owner.N,
            exact_Z_partition=CELLS,actual_correlated_bridge_and_Rc_cells=rows,
            quotient_binding=QUOTIENT_BINDING,inverse_binding=INVERSE_BINDING,query_binding=QUERY_BINDING,
            original_bounded_density_binding=bridge.DENSITY_BINDING,
            actual_terminal_controls_installed=False,actual_global_frequency_admitted=False,
            physical_original_exterior_five_targets_closed=False,**dict.fromkeys(OPEN,False),
            input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual same-function C/B/E and inverse-identity C1 bridge arithmetic; '
                'live zero-inlet bridge integrals, then exact accepted independent downstream local-driver '
                'hydration and new incoming propagation through Rc. No control, exterior closure, '
                'global N, higher finite-N jets or temporal recursion admission.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(serialized(report)),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__ == '__main__':run()
