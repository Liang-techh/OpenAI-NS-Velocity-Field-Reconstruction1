"""Exact original finite-mu pulse kernels in the common n1 repair graph.

Fixed source integrals remain defining integrals, never interval midpoints.
The positive selector accepts true same-graph incoming/future C1 functions;
it does not silently substitute a future energy bound for that function.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions as parent
import lei_ren_part1_paper_compliant_current_original_Rp_pulse_basis_energy as basis
import lei_ren_part1_paper_compliant_current_original_Rp_pulse_coefficients_refined as coefficients
import lei_ren_part1_paper_compliant_axial_pulse_field as field

source, controls = parent.source, parent.controls
HERE, PREFIX, sha = parent.HERE, parent.PREFIX, parent.sha
NAME = PREFIX + 'current_original_n1_pulse_kernel_functions.json.gz'
RECEIPT = PREFIX + 'current_original_n1_pulse_kernel_functions_check.json'
GATE = 'current_original_n1_exact_pulse_kernel_functions_installed'
OPEN = ('current_original_n1_full_pulse_C1_function_adapter_installed',
        'n1_complete_signed_defect_C1_function_adapter_installed',
        'n1_selected_linear_coefficient_functions_installed',
        'n1_repaired_own_M_V_P_installed', 'positive_order_five_moment_repair_closed',
        'temporal_recursion', 'full_corrected_Navier_Stokes_solution')


def checked(module, family, hashes):
    raw, receipt = parent.linear.read(module.NAME), parent.linear.read(module.RECEIPT)
    if (not receipt.get('all_passed') or not receipt.get(module.GATE)
            or receipt['source_family'] != family or raw['source_family'] != family):
        raise ValueError('Checked original pulse source required ' + module.NAME)
    incoming = dict(receipt['input_hashes'])
    for name in (module.NAME, module.RECEIPT, Path(module.__file__).name):
        incoming[name] = sha(name)
    for name, digest in incoming.items():
        if sha(name) != digest or (name in hashes and hashes[name] != digest):
            raise ValueError('Changed or conflicting pulse source ' + name)
        hashes[name] = digest
    return raw


def kernel_functions(g, mu):
    """Original gp, row integrals, divided matrix and nonzero end energy."""
    g.ids([mu])
    one, two = g.one, g.constant(2)
    ell = g.constant('3/20')
    t = g.symbol('original_pulse_beta_coordinate')
    beta = g.node('compact_raw_beta', argument=t.node,
        definition='exp(-1/(1-t^2)) on abs(t)<1, zero otherwise',
        defining_module=parent.linear.PREFIX+'outer_pulse_map.py',
        defining_module_sha256=sha(parent.linear.PREFIX+'outer_pulse_map.py'))
    integral = lambda body, var, lo, hi: controls.integral(g, body, var, lo, hi,
        endpoint_independent_of_Z=True, original_function_integral=True,
        numerical_enclosure_not_selected_value=True)
    norm = integral(beta, g.nodes[t.node]['name'], g.constant(-1), one)
    divide = lambda n, d, proof: g.quotient(n, d, proof)
    J = lambda k: divide(integral(g.mul(beta, g.unary('exp', g.mul(k, ell, t))),
        g.nodes[t.node]['name'], g.constant(-1), one), norm, 'positive original beta normalization')
    lam = [g.sub(g.constant('1/2'), g.mul(g.constant(i), mu)) for i in (1, 2)]
    j1, j2 = map(J, lam)
    A, D, centers = [], [], [g.constant(-3), g.constant(-1)]
    for center in centers:
        s = g.add(center, g.mul(ell, t))
        A.append(g.mul(g.unary('exp', g.mul(lam[0], center)), j1))
        # (exp(-mu*s)-1)/mu = -s*exprel(-mu*s), including finite mu.
        body = g.mul(beta, g.unary('exp', g.mul(lam[0], s)), g.neg(s),
                     g.unary('exprel', g.neg(g.mul(mu, s))))
        D.append(divide(integral(body, g.nodes[t.node]['name'], g.constant(-1), one),
                        norm, 'positive original beta normalization'))
    determinant_abs = g.mul(two, g.unary('exp', g.add(g.constant(-2), g.mul(g.constant(5), mu))),
                            j1, j2, g.unary('exprel', g.mul(two, mu)))
    gram = divide(integral(g.mul(beta, beta,
        g.unary('exp', g.neg(g.mul(two, mu, ell, t)))),
        g.nodes[t.node]['name'], g.constant(-1), one),
        g.mul(ell, norm, norm), 'positive ell and beta normalization squared')
    Kj = [g.mul(g.unary('exp', g.sub(g.constant(-26), g.mul(two, mu, c))), gram)
          for c in centers]
    xi = g.symbol('original_pulse_xi')
    v = g.symbol('original_pulse_gp_primitive_v')
    sigma = lambda x: g.node('exact_original_source_flat_sigma', argument=x.node,
        definition='0 for x<=0; 1 for x>=1; exp(-1/x^2)/(exp(-1/x^2)+exp(-1/(1-x)^2)) inside',
        source_binding=parent.active.current.ast_binding(field.sigma_enclosure),
        exact_function_not_its_interval_cap=True)
    primitive = integral(sigma(g.mul(g.constant(50), v)), g.nodes[v.node]['name'], g.zero, xi)
    gp = g.mul(primitive, sigma(g.sub(g.constant(11), xi)))
    energy = integral(g.mul(g.unary('exp', g.neg(g.mul(two, xi))), gp, gp),
                      g.nodes[xi.node]['name'], g.zero, g.constant(11))
    # The source normalization is an explicit choice based on u0=(4mu)^(1/3),
    # not a numerical minimizer/row enclosure used as a constant.
    logmu = g.unary('log', mu)
    logu0 = g.mul(g.constant('1/3'), g.unary('log', g.mul(g.constant(4), mu)))
    saddleL = g.unary('exp', g.neg(g.mul(two, logu0)))
    logpref = g.add(g.neg(divide(one, mu, 'positive original mu')),
        g.neg(g.mul(g.constant(3), saddleL)), g.mul(two, logu0),
        g.neg(g.mul(g.constant('1/2'), g.unary('log', g.constant(6)))), g.neg(logmu))
    logE = g.sub(logpref, logmu)
    Pi = []
    for l in lam:
        logweight = g.sub(divide(g.mul(l, g.sub(xi, g.constant(13))), mu,
                               'positive original mu'), logpref)
        Pi.append(divide(integral(g.mul(g.unary('exp', logweight), gp),
            g.nodes[xi.node]['name'], g.zero, g.constant(11)), mu, 'positive original mu'))
    nu = [g.mul(mu, k, g.unary('exp', g.mul(two, logE))) for k in Kj]
    factors = [g.unary('exp', g.sub(divide(g.mul(g.constant(-13), l), mu,
                                         'positive original mu'), logpref)) for l in lam]
    return dict(mu=mu,ell=ell,normalization=norm,lambda_rows=lam,centers=centers,
        first_row=A,divided_row=D,negative_divided_determinant=g.neg(determinant_abs),
        positive_determinant_absolute=determinant_abs,energy_gram=gram,
        original_end_energy_weights=Kj,log_common_prefactor=logpref,log_end_scale=logE,
        incoming_factors=factors,scaled_main_rows=Pi,fixed_pulse_energy=energy,
        exact_nonzero_end_energy=nu,original_gp=gp,original_gp_variable=xi)


def original_incoming_functions(g):
    """The actual source-bound native inlet, including its own histories."""
    raw=parent.linear.read(parent.quiet.FRAME_NAME)
    original=parent.linear.realization.restore_graph(raw['exact_graph_nodes'])
    importer=parent.GraphImport(g,original)
    packet=raw['actual_Rp_native_frame_functions']['actual_identified_Rp_native_frame_C3']
    if not packet['exact_common_m_k_already_normalized_by_Pstar']:
        raise ValueError('Actual native incoming normalization required')
    ref=lambda i:importer(source.FunctionRef(original,i))
    out={key:source.C1Function(*map(ref,packet[key][:2])) for key in ('u','m1','m2','energy')}
    return out


def select_functions(g, W, incoming, incoming_energy, future):
    """Original divided inverse and positive energy root, with true C1 rows."""
    pairs = list(incoming)+[incoming_energy,future]
    if len(incoming)!=2 or any(type(p) is not source.C1Function for p in pairs):
        raise TypeError('Two genuine incoming and two energy C1 functions required')
    for p in pairs:
        g.ids([p.value,p.Z])
    a = parent.C1Algebra(g)
    fixed = a.fixed
    A,D = W['first_row'],W['divided_row']
    inv = lambda r1,r2: [a.div(a.sub(a.scale(A[1],r2),a.scale(D[1],r1)),
        fixed(W['positive_determinant_absolute']),'positive absolute pulse determinant'),
        a.div(a.sub(a.scale(D[0],r1),a.scale(A[0],r2)),
        fixed(W['positive_determinant_absolute']),'positive absolute pulse determinant')]
    rows = [a.scale(k,p) for k,p in zip(W['incoming_factors'],incoming)]
    u = inv(a.neg(a.scale(W['mu'],rows[0])),a.neg(a.sub(rows[1],rows[0])))
    v = inv(a.neg(fixed(g.mul(W['mu'],W['scaled_main_rows'][0]))),
            a.neg(fixed(g.sub(*reversed(W['scaled_main_rows'])))))
    base = g.mul(g.constant('1/4'),g.sub(g.one,g.unary('exp',g.constant(-26))))
    target = a.add(future,a.neg(a.scale(W['mu'],incoming_energy)),fixed(base))
    A2 = a.add(fixed(W['fixed_pulse_energy']),
        *(a.scale(w,a.mul(p,p)) for w,p in zip(W['exact_nonzero_end_energy'],v)))
    A1 = a.add(*(a.scale(g.mul(g.constant(2),w),a.mul(p,q))
                for w,p,q in zip(W['exact_nonzero_end_energy'],u,v)))
    A0 = a.add(a.neg(target),*(a.scale(w,a.mul(p,p))
                             for w,p in zip(W['exact_nonzero_end_energy'],u)))
    half = a.div(A1,a.scale(2,A2),'positive original quadratic A2')
    disc = a.sub(a.mul(half,half),a.div(A0,A2,'positive original quadratic A2'))
    root = g.unary('sqrt',disc.value)
    sqrt = source.C1Function(root,g.quotient(disc.Z,g.mul(g.constant(2),root),
                                           'positive original quadratic discriminant square root'))
    ap = a.sub(sqrt,half)
    C = [a.add(p,a.mul(ap,q)) for p,q in zip(u,v)]
    return dict(ap=ap,C1=C[0],C2=C[1]),dict(A2=A2,A1=A1,A0=A0,
        discriminant=disc,incoming_rows=rows,target=target)


class CurrentOriginalN1PulseKernelFunctions:
    def __init__(self, preRp=None, require_checked=True):
        self.preRp = preRp or parent.SavedCurrentOriginalN1PreRpMomentFunctions()
        self.preRp.assert_graph()
        self.g,self.family = self.preRp.g,copy.deepcopy(self.preRp.family)
        self.hashes = dict(self.preRp.hashes)
        self.basis_raw = checked(basis,self.family,self.hashes)
        self.coefficients_raw = checked(coefficients,self.family,self.hashes)
        self.mu = source.FunctionRef(self.g,self.preRp.operator.background['parameters']['mu'])
        self.functions = kernel_functions(self.g,self.mu)
        self.incoming = original_incoming_functions(self.g)
        for module in (field,parent):
            self.hashes[Path(module.__file__).name] = sha(Path(module.__file__).name)
        self.hashes[Path(__file__).name] = sha(Path(__file__).name)
        self.prefix = copy.deepcopy(self.g.nodes)
        self._functions = parent.encode(self.functions)
        self._incoming = parent.encode(self.incoming)
        self._hashes = copy.deepcopy(self.hashes)
        self.acceptance_loaded = False
        if require_checked:
            receipt,raw = parent.linear.read(RECEIPT),parent.linear.read(NAME)
            if (not receipt.get('all_passed') or not receipt.get(GATE)
                    or receipt['source_family']!=self.family
                    or any(receipt.get(k) for k in OPEN)
                    or raw['exact_graph_nodes']!=self.g.nodes
                    or raw['exact_pulse_kernel_functions']!=self._functions
                    or raw['actual_native_incoming_C1_functions']!=self._incoming):
                raise ValueError('Checked exact original pulse kernel graph required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed pulse kernel source '+name)
            self.acceptance_loaded = True

    def assert_graph(self):
        self.preRp.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix
                or parent.encode(self.functions)!=self._functions
                or parent.encode(self.incoming)!=self._incoming or self.hashes!=self._hashes):
            raise ValueError('Original pulse kernel bodies changed')

    def coefficient_functions(self,incoming,incoming_energy,future):
        self.assert_graph()
        return select_functions(self.g,self.functions,incoming,incoming_energy,future)

    def original_coefficient_functions(self,future):
        self.assert_graph()
        return select_functions(self.g,self.functions,
            [self.incoming['m1'],self.incoming['m2']],self.incoming['energy'],future)


def run():
    began=time.monotonic()
    owner=CurrentOriginalN1PulseKernelFunctions(require_checked=False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,
        exact_pulse_kernel_functions=parent.encode(owner.functions),
        actual_native_incoming_C1_functions=parent.encode(owner.incoming),
        source_bound_basis_archive=basis.NAME,source_bound_coefficients_archive=coefficients.NAME,
        gp_definition='integral_0^xi sigma(50v) dv times sigma(11-xi)',
        determinant_identity='det=-2*exp(-2+5mu)*J(.5-mu)*J(.5-2mu)*exprel(2mu)',
        source_mu_nonzero=True,end_energy_integrals_not_zeroed=True,
        enclosures_not_selected_constants=True,full_future_C1_function_still_required=True,
        full_pulse_pressure_functions_not_installed=True,input_hashes=owner.hashes,
        **dict.fromkeys((GATE,)+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_EXACT_PULSE_KERNEL_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw


if __name__=='__main__':run()
