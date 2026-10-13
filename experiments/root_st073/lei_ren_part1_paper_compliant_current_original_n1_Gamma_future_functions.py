"""Actual Gamma heat deficits and complete infinite-tail C1 functions.

The original native angular history determines waiting and Rtail.  S is
exp(-logRtail), never a cap, a midpoint, or zero.  Positive Gamma integrals
define the functions; existing directed bounds remain separate evidence.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_pulse_kernel_functions as pulse
import lei_ren_part1_paper_compliant_exact_heat_component as heat_source
import lei_ren_part1_paper_compliant_outer_angular_candidate as angular_source
import lei_ren_part1_paper_compliant_angular_high_jets as high_jets
import lei_ren_part1_paper_compliant_outer_initial as flat_source

parent, source, controls = pulse.parent, pulse.source, pulse.controls
HERE, PREFIX, sha = pulse.HERE, pulse.PREFIX, pulse.sha
NAME = PREFIX + 'current_original_n1_Gamma_future_functions.json.gz'
RECEIPT = PREFIX + 'current_original_n1_Gamma_future_functions_check.json'
GATE = 'current_original_n1_exact_Gamma_three_future_C1_functions_installed'
OPEN = ('current_original_n1_complete_future_energy_C1_function_installed',) + pulse.OPEN


class C1Algebra(parent.C1Algebra):
    def exp(self, pair):
        value = self.g.unary('exp', pair.value)
        return source.C1Function(value, self.g.mul(value, pair.Z))

    def log(self, pair, positive='positive original logarithm argument'):
        return source.C1Function(self.g.unary('log', pair.value),
            self.g.quotient(pair.Z, pair.value, positive))

    def integral(self, pair, variable, lower, upper, **metadata):
        self.g.ids([pair.value, pair.Z, variable, lower, upper])
        name = self.g.nodes[variable.node]['name']
        return source.C1Function(*(controls.integral(self.g, ref, name, lower, upper,
            endpoint_independent_of_Z=True, ordinary_Z_derivative_order=order,
            same_defining_integrand_before_enclosure=True, **metadata)
            for order, ref in enumerate((pair.value, pair.Z))))


def sigma(g, argument):
    return g.node('exact_original_source_flat_sigma', argument=argument.node,
        definition='0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside',
        source_binding=parent.active.current.ast_binding(flat_source.stable_sigma),
        exact_function_not_its_interval_cap=True)


def phi(g, argument):
    return g.node('exact_original_source_heat_phi', argument=argument.node,
        definition='exp(-4/(3-t)^2) for t<3; zero for t>=3',
        source_binding=parent.active.current.ast_binding(high_jets.heat_future_jets),
        source_nested_function='phi',
        exact_function_not_its_interval_cap=True)


def angular_history_functions(g, mu, delta, epsilon, Xp, logRp):
    """Actual pulse memory, flatten/steep histories and scalar waiting root."""
    if type(Xp) is not source.C1Function:
        raise TypeError('Actual native angular C1 function required')
    g.ids([mu, delta, epsilon, Xp.value, Xp.Z, logRp])
    a = C1Algebra(g)
    f, c = a.fixed, g.constant
    z = source.C1Function(g.symbol('Z'), g.one)
    q = a.add(f(1), a.mul(z, z))
    r, k = g.sub(g.one, mu), g.sub(g.one, g.mul(c('1/2'), delta))
    eq = g.quotient(g.one, r, 'original 0<mu<1')
    Xv = a.add(f(eq), a.scale(g.unary('exp',
        g.quotient(g.mul(c(-13), r), mu, 'positive original mu')), a.sub(Xp, f(eq))))
    t = g.symbol('original_heat_history_t')
    v = g.symbol('original_heat_history_primitive_v')
    logratio = a.log(a.scale('1/2', q), 'q/2>=1/2')
    F = lambda x: a.exp(a.scale(sigma(g, g.mul(c('1/100'), x)), logratio))
    flatten_primitive = a.integral(a.scale(g.unary('exp', g.mul(r, v)), F(v)),
        v, g.zero, c(100), original_angular_transport_integral=True)
    Xf = a.div(a.scale(g.unary('exp', g.mul(c(-100), r)),
        a.add(Xv, flatten_primitive)), F(c(100)), 'positive original flatten factor')
    logmu = g.unary('log', mu)
    Lrel = g.mul(c(-30), logmu)
    Xrel = a.add(f(eq), a.scale(g.unary('exp', g.neg(g.mul(r, Lrel))), a.sub(Xf, f(eq))))
    J = controls.integral(g, sigma(g, v), g.nodes[v.node]['name'], g.zero, t,
        original_unit_steep_primitive=True, endpoint_independent_of_Z=True)
    Iin = controls.integral(g, g.unary('exp', g.mul(r, g.sub(t, J))),
        g.nodes[t.node]['name'], g.zero, g.one, original_unit_steep_integral=True)
    Iout = controls.integral(g, g.unary('exp', g.mul(k, J)),
        g.nodes[t.node]['name'], g.zero, g.one, original_unit_steep_integral=True)
    Ts = g.mul(c(4), g.sub(g.unary('log', c(2)), g.unary('log', delta)))
    Xin = a.scale(g.unary('exp', g.mul(c('-1/2'), r)), a.add(Xrel, f(Iin)))
    Xs = a.add(Xin, f(Ts))
    Xt = a.scale(g.unary('exp', g.mul(c('-1/2'), k)), a.add(Xs, f(Iout)))
    sig, ph = sigma(g, t), phi(g, t)
    W = g.add(g.sub(g.one, sig), g.mul(sig, ph))
    Jcollar = controls.integral(g, g.mul(g.unary('exp', g.mul(k, t)), W),
        g.nodes[t.node]['name'], g.zero, c(3), original_preheat_collar_atom=True)
    Xt0 = g.node('function_substitution', expression=Xt.value.node,
        variable=z.value.node, value=g.zero.node, Z_independent_substitution=True)
    eqk = g.quotient(g.one, k, 'original 0<delta<=1/2')
    logone = g.unary('log1p', g.neg(epsilon))
    waiting = g.quotient(g.add(g.unary('log', g.sub(Xt0, eqk)), logone,
        g.neg(g.unary('log', epsilon)), g.neg(g.unary('log', g.add(eqk, Jcollar)))),
        k, 'positive waiting rate; Xt(0)-1/k>0 from checked original candidate')
    logRtail = g.add(logRp, g.quotient(c(13), mu, 'positive original mu'),
        c(102), Lrel, Ts, waiting)
    S = g.unary('exp', g.neg(logRtail))
    atoms = {}
    for name, rate in (('energy', delta), ('pressure', g.add(g.one, delta))):
        weight = g.unary('exp', g.neg(g.mul(rate, t)))
        atoms[name] = {key: controls.integral(g, g.mul(weight, body),
            g.nodes[t.node]['name'], g.zero, c(3), original_preheat_collar_atom=True)
            for key, body in (('W', W), ('W_squared', g.mul(W, W)))}
    return dict(mu=mu, delta=delta, epsilon=epsilon, r=r, k=k, Lrel=Lrel, Ts=Ts,
        waiting=waiting, log_one_minus_epsilon=logone, Jcollar=Jcollar,
        logRp=logRp, logRtail=logRtail, S=S, preheat_atoms=atoms), dict(
        Xp=Xp, Xv=Xv, Xf=Xf, Xrel=Xrel, steep_in=Xin, steep_power=Xs, Xt=Xt,
        unit_sigma_primitive=J, steep_in_integral=Iin, steep_out_integral=Iout)


def Gamma_future_functions(g, delta, epsilon, S):
    """Full positive Gamma expectation and all three exact C1 tails."""
    g.ids([delta, epsilon, S])
    c, a = g.constant, C1Algebra(g)
    fixed = a.fixed
    z = source.C1Function(g.symbol('Z'), g.one)
    d = a.sub(fixed(1), a.mul(z, z))
    aa = g.mul(c('1/2'), delta)
    t = g.symbol('original_Gamma_log_radius_offset')
    v = g.symbol('original_Gamma_expectation_v')
    infinity = g.node('mathematical_positive_infinity', exact_improper_endpoint=True)
    decay = g.unary('exp', g.neg(t))
    xi = g.mul(c(2), d.value, S, decay)
    x = g.mul(xi, v)
    log1px = g.unary('log1p', x)
    weight = g.unary('exp', g.add(g.neg(v),
        g.mul(g.add(aa, g.one), g.unary('log', v))))
    normalization = g.unary('gamma', g.add(g.one, aa))
    expectation = lambda body: controls.integral(g, body, g.nodes[v.node]['name'],
        g.zero, infinity, exact_positive_Gamma_expectation=True,
        full_improper_integral=True, Gamma_normalization=normalization.node,
        zero_endpoint_integrand_continuation=0)
    # 1-(1+x)^(-a) = a*x*log1prel(x)*exprel(-a*log1p(x)).
    # Both divided functions have their exact continuation 1 at x=0.
    Dvalue = g.quotient(g.mul(c(2), d.value, decay, expectation(g.mul(weight,
        g.unary('log1prel', x), g.unary('exprel', g.neg(g.mul(aa, log1px)))))),
        normalization, 'positive Gamma(1+a)')
    DZ = g.quotient(g.mul(c(-4), z.value, decay, expectation(g.mul(weight,
        g.unary('exp', g.neg(g.mul(g.add(aa, g.one), log1px)))))),
        normalization, 'positive Gamma(1+a)')
    Dhat = source.C1Function(Dvalue, DZ)
    sig, ph = sigma(g, t), phi(g, t)
    C = g.sub(g.one, g.mul(epsilon, ph))
    W = g.add(g.sub(g.one, sig), g.mul(sig, ph))
    pre = g.sub(g.one, g.mul(epsilon, W))
    actualD = a.scale(g.mul(sig, C), Dhat)
    K = a.sub(fixed(pre), a.scale(g.mul(aa, S), actualD))
    square = a.mul(actualD, a.sub(fixed(g.mul(c(2), pre)), a.scale(g.mul(aa, S), actualD)))
    exterior_square = a.mul(Dhat, a.sub(fixed(2), a.scale(g.mul(aa, S), Dhat)))
    weights = dict(theta=g.mul(aa, g.unary('exp', g.mul(g.sub(g.one, aa), t))),
        pressure=g.mul(c('1/2'), g.unary('exp', g.neg(g.mul(g.add(g.one, delta), t)))),
        energy=g.unary('exp', g.neg(g.mul(delta, t))))
    collar, exterior, total = {}, {}, {}
    for name, weight_ref in weights.items():
        first, last = (actualD, Dhat) if name == 'theta' else (square, exterior_square)
        metadata = dict(full_exact_Gamma_future_defect=name,
            scaled_heat_deficit_not_its_cap=True, differentiation_under_integral=True)
        collar[name] = a.integral(a.scale(weight_ref, first), t, g.zero, c(3), **metadata)
        exterior[name] = a.integral(a.scale(weight_ref, last), t, c(3), infinity,
            full_improper_integral=True, dominated_by_exponentially_integrable_Gamma_bound=True,
            strict_delta_positive_required=True, **metadata)
        total[name] = a.add(collar[name], exterior[name])
    return dict(Dhat=Dhat, theta=total['theta'], pressure=total['pressure'], energy=total['energy']), dict(
        collar=collar, exterior=exterior, actual_scaled_deficit=actualD,
        actual_bracket=K, preheat_bracket=pre, squared_difference_scaled=square,
        exterior_squared_difference_scaled=exterior_square, offset_variable=t,
        Gamma_variable=v, xi=xi, a=aa, infinity=infinity)


def checked_sources(family, hashes):
    gates = {
        'exact_heat_component': ('exact_selected_heat_component_specified',),
        'outer_angular_candidate': ('actual_waiting_root_in_original_same_source_enclosure',
            'original_continuous_preheat_waiting_log_equation_encloses_zero'),
        'current_heat_source': ('all_passed', 'current_heat_source_ownership_certified')}
    records = {}
    for stem, required in gates.items():
        name, receipt_name = PREFIX + stem + '.json', PREFIX + stem + '_check.json'
        raw, receipt = parent.linear.read(name), parent.linear.read(receipt_name)
        if not all(receipt.get(key) for key in required):
            raise ValueError('Original Gamma/waiting acceptance required: ' + receipt_name)
        for key, value in family.items():
            if raw.get(key) != value or (key in receipt and receipt[key] != value):
                raise ValueError('Original Gamma/waiting source family differs: ' + key)
        incoming = {**raw['input_hashes'], **receipt['input_hashes'],
            name: sha(name), receipt_name: sha(receipt_name)}
        for path, digest in incoming.items():
            if sha(path) != digest or (path in hashes and hashes[path] != digest):
                raise ValueError('Changed Gamma/waiting source: ' + path)
            hashes[path] = digest
        records[stem] = dict(raw=name, receipt=receipt_name, required_gates=required)
    return records


class CurrentOriginalN1GammaFutureFunctions:
    def __init__(self, pulse_owner=None, require_checked=True):
        self.pulse = pulse_owner or pulse.CurrentOriginalN1PulseKernelFunctions()
        self.pulse.assert_graph()
        self.g, self.family = self.pulse.g, copy.deepcopy(self.pulse.family)
        self.hashes = dict(self.pulse.hashes)
        self.source_receipts = checked_sources(self.family, self.hashes)
        raw = parent.linear.read(parent.quiet.FRAME_NAME)
        original = parent.linear.realization.restore_graph(raw['exact_graph_nodes'])
        importer = parent.GraphImport(self.g, original)
        ref = lambda i: importer(source.FunctionRef(original, i))
        frame = raw['actual_Rp_native_frame_functions']
        native = frame['actual_identified_Rp_native_frame_C3']
        Xp = source.C1Function(*map(ref, native['X'][:2]))
        logRp = ref(frame['original_logRp'])
        logP = source.FunctionRef(self.g, self.pulse.preRp.operator.background['parameters']['logP'])
        # Checked Md=40 selects the strict exponential branch of min(1e-200,...).
        delta = self.g.unary('exp', self.g.neg(self.g.add(self.g.mul(self.g.constant(4), logP),
            self.g.constant(30))))
        epsilon = self.g.mul(self.g.constant('1/1000'), delta)
        self.schedule, self.histories = angular_history_functions(self.g, self.pulse.mu,
            delta, epsilon, Xp, logRp)
        self.functions, self.pieces = Gamma_future_functions(self.g, delta, epsilon, self.schedule['S'])
        self.native_binding = dict(frame_archive=parent.quiet.FRAME_NAME,
            original_X_C1=native['X'][:2], original_logRp=frame['original_logRp'],
            imported_X_C1=parent.encode(Xp), imported_logRp=logRp.node,
            delta_definition='exp(-4*(exp(40)+11)-30)',
            epsilon_definition='delta/1000', original_mu_preserved=True)
        for module in (heat_source, angular_source, high_jets, flat_source, pulse, parent):
            self.hashes[Path(module.__file__).name] = sha(Path(module.__file__).name)
        self.hashes[Path(__file__).name] = sha(Path(__file__).name)
        self.prefix = copy.deepcopy(self.g.nodes)
        self.snapshot = self._snapshot()
        self.acceptance_loaded = False
        if require_checked:
            receipt, saved = parent.linear.read(RECEIPT), parent.linear.read(NAME)
            if (not receipt.get('all_passed') or not receipt.get(GATE)
                    or receipt['source_family'] != self.family
                    or any(receipt.get(key) for key in OPEN)
                    or saved['exact_graph_nodes'] != self.g.nodes
                    or saved['function_bindings'] != self.snapshot):
                raise ValueError('Checked original exact Gamma future graph required')
            for name, digest in receipt['input_hashes'].items():
                if sha(name) != digest:
                    raise ValueError('Changed Gamma future source: ' + name)
            self.acceptance_loaded = True

    def _snapshot(self):
        return parent.encode(dict(functions=self.functions, pieces=self.pieces,
            schedule=self.schedule, histories=self.histories, native_binding=self.native_binding,
            source_receipts=self.source_receipts, input_hashes=self.hashes))

    def assert_graph(self):
        self.pulse.assert_graph()
        if self.g.nodes[:len(self.prefix)] != self.prefix or self._snapshot() != self.snapshot:
            raise ValueError('Original exact Gamma/future defining functions changed')


def run():
    began = time.monotonic()
    owner = CurrentOriginalN1GammaFutureFunctions(require_checked=False)
    raw = dict(source_family=owner.family, exact_graph_nodes=owner.g.nodes,
        function_bindings=owner.snapshot, exact_inverse_radius_not_zero=True,
        native_nonzero_pulse_angular_memory_retained=True, actual_waiting_root_defined=True,
        all_three_complete_infinite_exterior_C1_integrals_defined=True,
        interval_caps_not_selected_as_defining_values=True, finite_Gamma_series_not_used=True,
        full_future_energy_and_angular_branch_still_required=True,
        **dict.fromkeys((GATE,) + OPEN, False), input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw, separators=(',', ':'))+'\n').encode(), mtime=0))
    print('ORIGINAL_N1_EXACT_GAMMA_FUTURE_FUNCTIONS_DEFINED', len(owner.g.nodes), flush=True)
    return owner, raw


if __name__ == '__main__':
    run()
