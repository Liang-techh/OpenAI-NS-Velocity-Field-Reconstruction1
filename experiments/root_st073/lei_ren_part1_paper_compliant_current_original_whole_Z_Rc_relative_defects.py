"""Actual fixed-N Rc defects relative to the unchanged leading input.

Appendix C.2 / Lei-Ren Section 11 restores the moments of that same input,
not a newly chosen particular power-law history. The signed integral graph
therefore defines the relative defect. This does not prove that the input
itself already satisfies every original exterior terminal condition.
"""
import ast
from fractions import Fraction
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_whole_Z_full_signed_transport as graph
import lei_ren_part1_paper_compliant_current_generic_moment_repair_operator as repair
from lei_ren_part1_paper_compliant_current_original_R100_endpoint import restore_row

HERE, PREFIX, sha, bind, ep = graph.HERE, graph.PREFIX, graph.sha, graph.bind, graph.ep
CELLS, RATES, OPEN = graph.CELLS, graph.RATES, graph.OPEN
serialized, encode = graph.serialized, graph.encode
NAME = PREFIX+'current_original_whole_Z_Rc_relative_defects.json.gz'
RECEIPT = PREFIX+'current_original_whole_Z_Rc_relative_defects_check.json'
GATE = 'current_original_whole_Z_fixed_N_relative_Rc_C1_defects_and_targets_installed'
REFERENCE = 'the identical unchanged leading input of Appendix C.2 / Section 11'


def exact_row(a, b):
    return a.ctx is b.ctx and a.scale.bases is b.scale.bases and a.ledger is b.ledger \
        and a.coefficient._mpi_ == b.coefficient._mpi_ and a.scale.powers == b.scale.powers \
        and a.scale.offset._mpi_ == b.scale.offset._mpi_ and a.zero == b.zero


def same_source(f, rows):
    if any(v.ctx is not f.c or v.scale.bases is not f.logs or v.ledger is not f.ledger for v in rows):
        raise ValueError('One actual source context, formal basis and ledger required')


def reference_binding():
    """Bind the fixed-input subtraction, without an exterior-target claim."""
    fn = graph.WholeZFullSignedTransport.outer_transport
    node = ast.parse(textwrap.dedent(inspect.getsource(fn))).body[0]
    statement = "complete = {name: [background[name][i]+outgoing[name][i] for i in range(2)] for name in RATES}"
    wanted = ast.dump(ast.parse(statement).body[0])
    if sum(ast.dump(n) == wanted for n in ast.walk(node) if isinstance(n, ast.Assign)) != 1:
        raise ValueError('Same leading input plus current signed correction identity changed')
    recovery = repair.packets.recovery
    if {k:Fraction(str(v)) for k,v in recovery.RATES.items()} != RATES:
        raise ValueError('Original normalized moment rates changed')
    return dict(reference=REFERENCE, original_paper='Appendix C.2, equations C.15-C.17',
        independent_reconstruction='Lei-Ren Part I v2, Section 11',
        fixed_input_not_a_quiet_power_particular_solution=True,
        exact_complete_history_assignment=statement,
        original_transport_AST_sha256=hashlib.sha256(ast.dump(node).encode()).hexdigest(),
        original_normalized_units=recovery.UNITS,
        own_history_difference='(H_original+D)-H_original=D as functions, before interval evaluation',
        pressure_difference='(P0+p_original+D_p)-(P0+p_original)=D_p',
        same_axis_P0_cancels_only_in_relative_pressure_difference=True,
        original_exterior_target_and_heat_pressure_identity_not_proved_here=True)


def normalized_targets(f, defects, amplitude, mu, logmu, N):
    """Normalize actual fixed-N totals, with the joint numerator first."""
    c = f.c
    if set(defects) != set(RATES) or any(len(row) != 2 for row in defects.values()) or len(amplitude) != 2:
        raise ValueError('Five actual C0/Z pairs and source amplitude C0/Z required')
    same_source(f, [v for pair in defects.values() for v in pair]+list(amplitude)+[mu])
    A, AZ = amplitude
    if ep(A.coefficient)[0] <= 0 or ep(mu.coefficient)[0] <= 0:
        raise ValueError('The actual source amplitude and mu must be strictly positive')
    logA = c.mpf(ep(A.scale.evaluate()+c.ln(A.coefficient))[0])
    ratio = AZ.positive_divide(A, logA)
    den = A*A
    values, jets = {}, {}
    for out, key, degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2)):
        D, DZ = defects[key]
        divisor = A if degree == 1 else den
        values[out] = D.positive_divide(divisor, degree*logA)
        jets[out] = (DZ-D*ratio*degree).positive_divide(divisor, degree*logA)
    numerator = defects['k'][0]-A*defects['m'][0]
    numeratorZ = defects['k'][1]-AZ*defects['m'][0]-A*defects['m'][1]
    divided = repair.ROWS[1]
    values[divided] = numerator.positive_divide(den*mu, 2*logA+logmu)
    jets[divided] = (numeratorZ-numerator*ratio*2).positive_divide(den*mu, 2*logA+logmu)
    targets = {key:[values[key], jets[key]] for key in repair.ROWS}
    scaled = {key:[v*c.mpf(N) for v in pair] for key,pair in targets.items()}
    caps = {}
    for key, pair in scaled.items():
        terms = [repair.LogUpper(c, None if v.zero else v.record()['log_absolute_upper']) for v in pair]
        caps[key] = repair.LogUpper.add(c, terms)
    nonzero = [cap for cap in caps.values() if cap.log is not None]
    maximum = repair.LogUpper(c, c.mpf(max(ep(cap.log)[1] for cap in nonzero))) if nonzero else repair.LogUpper(c,None)
    return dict(actual_fixed_N_normalized_target_C0_Z=targets,
        actual_fixed_N_N_scaled_target_C0_Z=scaled,
        actual_joint_axial_numerator_C0_Z=[numerator, numeratorZ],
        actual_positive_A_rc_C0_Z=amplitude, actual_positive_mu=mu,
        actual_positive_A_rc_log_lower=logA, actual_positive_mu_log=logmu,
        actual_A_rc_Z_over_A_rc=ratio,
        target_C1_caps=dict(transformed_N_scaled_target_C1_caps={key:cap.record() for key,cap in caps.items()},
            whole_target_C1_cap=maximum.record(), norm='max_i(sup|N*r_i|+sup|N*r_i_Z|)',
            fixed_candidate_N=N, uniform_across_changed_N=False),
        source_joint_numerator_formed_before_division_by_actual_mu=True,
        interval_hulls_do_not_prove_joint_cancellation=True)


def physical_differences(f, radius, defects):
    S = f.factor((0,.5,0,0,0))
    # Rm contains half-powers, so R^(3/2) can contain quarter-powers.
    # The existing factor algebra intentionally accepts only half-powers.
    # Collect this positive Z-independent physical unit in a formal log
    # offset instead of expanding an astronomical physical exponential.
    radial32 = f.factor((0,0,0,0,0),
        (radius.scale.evaluate()+f.c.ln(radius.coefficient))*f.c.mpf('1.5'))
    units = dict(m=radius*S, h=radial32*S*f.c.sqrt(2),
        k=radial32*S*S*f.c.sqrt(2), e=radius*S*S, p=S*S)
    return {key:[v*units[key] for v in pair] for key,pair in defects.items()}


class WholeZRcRelativeDefects:
    def __init__(self, dps=500):
        self.graph = graph.WholeZFullSignedTransport(dps)
        self.c, self.identity, self.N = self.graph.c, self.graph.identity, self.graph.N
        self.hashes = dict(self.graph.hashes)
        receipt = json.loads((HERE/graph.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(graph.GATE) \
                or receipt['source_family'] != self.identity or receipt['candidate_N'] != self.N:
            raise ValueError('Accepted actual same-source signed full graph required')
        for name,digest in receipt['input_hashes'].items():
            bind(self.hashes,name,digest)
        bind(self.hashes,graph.RECEIPT,sha(graph.RECEIPT))
        self.accepted = json.loads(gzip.decompress((HERE/graph.NAME).read_bytes()))
        if not self.accepted.get(graph.GATE) or self.accepted['source_family'] != self.identity \
                or self.accepted['candidate_N'] != self.N:
            raise ValueError('Accepted graph function evaluation has different source/N')
        self.cells = {tuple(row['exact_Z_cell']):row for row in self.accepted['actual_signed_full_transports']}
        if set(self.cells) != set(CELLS):
            raise ValueError('Complete accepted whole-Z evaluation partition required')
        self.binding = reference_binding()
        for module in (repair, repair.packets.recovery):
            name = Path(module.__file__).name
            bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        self.weights = repair.fresh_weights(self.c,self.c.mpf([0,ep(self.c.mpf(1)/6)[1]]))
        self.matrix = repair.fresh_linear_inverse(self.c,self.c.mpf([0,ep(self.c.mpf(1)/6)[1]]),self.weights)

    def _domain(self, ends):
        key = tuple(ends)
        if key not in self.cells:
            raise ValueError('One admitted current whole-Z cell required; no saved point-label fallback')
        return key

    def _consume(self, ends, transport, mode):
        op = self.graph.owner(ends)
        f, c = op.flow, op.c
        if transport['source_identity'] != self.identity or transport['candidate_N'] != self.N \
                or tuple(transport['exact_Z_cell']) != tuple(ends):
            raise ValueError('Identical source family/N/Z cell required')
        if not transport.get('actual_same_source_signed_C0_Z_integral_graph_through_Rc'):
            raise ValueError('Actual defining zero-inlet through Rc graph required')
        P0 = transport['exact_common_P0_axial5']
        if len(P0) != 6 or any(not exact_row(a,b) for a,b in zip(P0,op.P0)):
            raise ValueError('Same independent source P0 exact tuples required')
        outer = transport['actual_source_owned_zero_inlet_to_Rc_signed_stages']['outer']
        radius = op.Rm_factor*f.factor((0,.5,0,0,0),9)
        if not exact_row(transport['exact_Rc_radius'],radius):
            raise ValueError('Actual Rc=Rm*Pstar*exp(9) required')
        D = transport['actual_Rc_signed_correction_C0_Z']
        if any(not exact_row(a,b) for name in RATES for a,b in zip(D[name],outer['actual_Rc_signed_correction_C0_Z'][name])):
            raise ValueError('Current full graph output required, not a saved target')
        A = outer['actual_positive_Rc_normalized_amplitude_C0_Z']
        logmu = self.graph.source.logmu
        mu = f.factor((0,0,0,0,0),logmu)
        target = normalized_targets(f,D,A,mu,logmu,self.N)
        # The current signed source has already passed its actual-N
        # exponent budget. Test the additional repair condition separately.
        conditions = repair.contraction_log_conditions(c,target['target_C1_caps'],self.matrix,self.weights,logmu,c.mpf(0))
        sufficient = ep(c.ln(self.N))[0] >= ep(conditions['repair_sufficient_common_log_N_lower'])[1]
        # The existing theorem is conditional. Do not claim a control
        # solution when this actual fixed-N source cap misses its threshold.
        conditions.update(actual_fixed_candidate_N=self.N, actual_fixed_log_N=c.ln(self.N),
            same_fixed_N_signed_source_budget_already_accepted=True,
            fixed_candidate_N_satisfies_sufficient_contraction_bound=sufficient,
            sufficient_test_failure_is_not_a_nonexistence_proof=True,
            threshold_from_this_fixed_N_cap_cannot_select_a_changed_N=True,
            actual_signed_defect_functions_and_control_functions_not_evaluated=False,
            actual_signed_defect_functions_evaluated=True, actual_control_functions_evaluated=False,
            unique_C1_controls_for_this_fixed_candidate_N_certified=sufficient)
        return dict(source_identity=self.identity,exact_Z_cell=list(ends),candidate_N=self.N,
            reference_binding=self.binding, evaluation_mode=mode,
            exact_common_P0_axial5=op.P0,exact_Rc_radius=radius,
            actual_relative_Rc_normalized_five_defects_C0_Z=D,
            actual_relative_Rc_physical_five_moment_differences_C0_Z=physical_differences(f,radius,D),
            **target, actual_fixed_N_repair_log_conditions=conditions,
            complete_minus_identical_input_cancelled_as_functions_not_independent_interval_hulls=True,
            actual_Rc_relative_defect_functions_installed=True,
            physical_original_exterior_five_targets_closed=False,
            actual_terminal_controls_installed=False, actual_global_frequency_admitted=False,
            **dict.fromkeys(OPEN,False))

    def query(self, ends):
        """Live defining function: evaluate the signed graph, then normalize."""
        key = self._domain(ends)
        return self._consume(key,self.graph.full_transport(key),'live_defining_signed_integral_graph')

    def accepted_evaluation(self, ends):
        """Reuse a hash-bound, freshly checked function evaluation explicitly.

        This method never turns support caps or chosen enclosure endpoints
        into source functions. query() remains the actual defining callable.
        """
        key = self._domain(ends)
        stored = self.cells[key]
        f = self.graph.owner(key).flow
        decode = lambda row: restore_row(f,row)
        pairs = lambda rows:{name:[decode(v) for v in pair] for name,pair in rows.items()}
        outer = stored['actual_source_owned_zero_inlet_to_Rc_signed_stages']['outer']
        D = pairs(stored['actual_Rc_signed_correction_C0_Z'])
        packet = dict(source_identity=stored['source_identity'],candidate_N=stored['candidate_N'],
            exact_Z_cell=stored['exact_Z_cell'],actual_same_source_signed_C0_Z_integral_graph_through_Rc=
                stored['actual_same_source_signed_C0_Z_integral_graph_through_Rc'],
            exact_common_P0_axial5=[decode(v) for v in stored['exact_common_P0_axial5']],
            exact_Rc_radius=decode(stored['exact_Rc_radius']),actual_Rc_signed_correction_C0_Z=D,
            actual_source_owned_zero_inlet_to_Rc_signed_stages=dict(outer=dict(
                actual_Rc_signed_correction_C0_Z=pairs(outer['actual_Rc_signed_correction_C0_Z']),
                actual_positive_Rc_normalized_amplitude_C0_Z=[decode(v) for v in outer['actual_positive_Rc_normalized_amplitude_C0_Z']])))
        return self._consume(key,packet,'accepted_current_full_graph_evaluation_exact_MPI_hydration')


def exact_quotient_theorem():
    z = sy.symbols('z')
    A,D,m,k = [sy.Function(name)(z) for name in ('A','D','m','k')]
    mu = sy.symbols('mu',positive=True)
    identities = {}
    for degree in (1,2):
        expression = sy.diff(D/A**degree,z)-(sy.diff(D,z)-degree*sy.diff(A,z)*D/A)/A**degree
        if sy.simplify(expression) != 0: raise ArithmeticError('Amplitude quotient Z rule differs')
        identities['degree_'+str(degree)] = True
    C = k-A*m
    if sy.simplify(sy.diff(C/(mu*A*A),z)-(sy.diff(k,z)-sy.diff(A,z)*m-A*sy.diff(m,z)-2*sy.diff(A,z)*C/A)/(mu*A*A)) != 0:
        raise ArithmeticError('Joint divided-row Z rule differs')
    identities['joint_axial_divided_Z'] = True
    return dict(passed=True,identities=identities,mu_and_radial_limits_Z_independent=True)


def run():
    began = time.monotonic()
    with mp.workdps(540):
        owner = WholeZRcRelativeDefects()
        rows = []
        for ends in CELLS:
            rows.append(serialized(owner.accepted_evaluation(ends)))
            print('Actual same-input Rc relative defects/targets: '+str(ends),flush=True)
        report = dict(**{GATE:True},source_family=owner.identity,candidate_N=owner.N,
            exact_Z_partition=CELLS,actual_Rc_relative_defect_cells=rows,
            exact_quotient_Z_theorem=exact_quotient_theorem(),reference_binding=owner.binding,
            original_exact_bump_weight_enclosures=owner.weights,original_exact_divided_inverse_enclosures=owner.matrix,
            actual_Rc_relative_defect_functions_installed=True,actual_terminal_controls_installed=False,
            actual_global_frequency_admitted=False,physical_original_exterior_five_targets_closed=False,
            **dict.fromkeys(OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Callable actual fixed-N relative Rc C1 defects from the current zero-inlet graph; '
                'exact accepted function-evaluation hydration; live source A_rc/mu/P0/radius guards; '
                'physical differences, joint normalization and actual fixed-N contraction diagnostic. '
                'Controls, input exterior closure, higher finite-N jets, global N/heat/cone and temporal recursion remain open.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(serialized(report)),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__ == '__main__': run()
