"""Correlated actual quiet pressure moments with a directed scalar kernel.

The physical Pstar cancels the exact original exp(-logP) before evaluation.
The absolute original Rc remains an exact source scale; no radius is
materialized. The original O2 turnoff mean kernel is evaluated around its
actual saddle, retaining positive tails and a quadrature error enclosure.
"""
import copy
import gzip
import json
import math
from pathlib import Path
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_primitives as prior
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE, PREFIX, sha = prior.HERE, prior.PREFIX, prior.sha
NAME = PREFIX + 'current_original_n1_quiet_pressure_source.json.gz'
RECEIPT = PREFIX + 'current_original_n1_quiet_pressure_source_check.json'
CONSTANTS = PREFIX + 'current_original_Rp_native_constants.json'
CONSTANTS_RECEIPT = PREFIX + 'current_original_Rp_native_constants_check.json'
CONSTANTS_SOURCE = PREFIX + 'current_original_Rp_native_constants.py'
GATE = 'current_original_n1_quiet_source_scaled_signed_pressure_moments_installed'
OPEN = ('complete_signed_pressure_tails_certified', 'positive_order_five_moment_repair_closed',
        'temporal_recursion', 'full_corrected_Navier_Stokes_solution',
        'absolute_physical_radial_scale_materialization_installed')


def ep(value):
    return tuple(mp.make_mpf(item) for item in value._mpi_)


def record(value):
    lo, hi = ep(value)
    return dict(lower=str(lo), upper=str(hi),
                lower_exact_mpf_tuple=list(value._mpi_[0]),
                upper_exact_mpf_tuple=list(value._mpi_[1]))


def read_interval(c, value):
    return c.make_mpf((tuple(value['lower_exact_mpf_tuple']),
                       tuple(value['upper_exact_mpf_tuple'])))


def _positive(value, name):
    if ep(value)[0] <= 0:
        raise ValueError('Positive directed source enclosure required: ' + name)
    return value


def kernel_profile(c, Y, point):
    """H=s+D, for sigma(x)=exp(-D)/(1+exp(-D)), x=-log(1-s/Y)/40."""
    x = -c.ln(1 - point / Y) / 40
    _positive(x, 'turnoff x')
    D = 1 / x**2 - 1 / (1 - x)**2
    dx = 1 / (40 * (Y - point))
    Dx = -2 / x**3 - 2 / (1 - x)**3
    Dxx = 6 / x**4 - 6 / (1 - x)**4
    return dict(x=x, D=D, H=point + D,
                H1=1 + Dx * dx,
                H2=(Dxx / 1600 + Dx / 40) / (Y - point)**2)


@source_precision
def original_turnoff_kernel(c=None, cells=2048):
    """Directed full original kernel at Y=exp(40), including both tails.

    On x<=1/100 the explicit H'' is positive. Monotonic endpoint bounds
    control each normalized f=e^(H(center)-H(s)). The local trapezoid error
    is bounded by h^3 sup|f''|/12, with f''=f*(H'^2-H'').
    """
    if type(cells) is not int or cells < 128:
        raise ValueError('At least 128 exact uniform quadrature cells required')
    c = c if c is not None else MPIntervalContext()
    c.dps = max(c.dps, 140)
    Y = c.exp(40)
    with mp.workdps(100):
        y = mp.exp(40)
        def derivative(u):
            x = -mp.log1p(-u / y) / 40
            return 1 - (2 / x**3 + 2 / (1 - x)**3) / (40 * (y - u))
        seed = (2 * (40*y)**2)**(mp.mpf(1)/3)
        point = mp.findroot(derivative, (seed, seed*(1+mp.mpf('1e-4'))))
        center_decimal = mp.nstr(point, 75)
    center = c.mpf(center_decimal)
    epsilon = c.mpf('1e-40')
    root_box = center + c.mpf([-1, 1]) * epsilon
    left_root, right_root = center-epsilon, center+epsilon
    if not (ep(kernel_profile(c,Y,left_root)['H1'])[1] < 0 <
            ep(kernel_profile(c,Y,right_root)['H1'])[0]):
        raise ValueError('Actual source saddle bracket not certified')
    width = c.sqrt(center)
    a, b = center-12*width, center+12*width
    xa, xb = kernel_profile(c,Y,a)['x'], kernel_profile(c,Y,b)['x']
    if not (ep(a)[0] > 0 and ep(xb)[1] < mp.mpf('0.01')):
        raise ValueError('Original positive convex-kernel window required')
    reference = kernel_profile(c,Y,center)['H']
    saddle_H = kernel_profile(c,Y,root_box)['H']
    step = (b-a)/cells
    total, error = c.mpf(0), c.mpf(0)
    previous = kernel_profile(c,Y,a)
    previous_f = c.exp(reference-previous['H'])
    for index in range(cells):
        left, right = a+index*step, a+(index+1)*step
        endpoint = kernel_profile(c,Y,right)
        endpoint_f = c.exp(reference-endpoint['H'])
        if ep(right)[1] <= ep(left_root)[0]:
            maximum = endpoint_f
        elif ep(left)[0] >= ep(right_root)[1]:
            maximum = previous_f
        else:
            maximum = c.exp(reference-saddle_H)
        box = c.mpf([ep(left)[0],ep(right)[1]])
        derivatives = kernel_profile(c,Y,box)
        _positive(derivatives['H2'], 'local H convexity')
        M2 = maximum*(derivatives['H1']**2+derivatives['H2'])
        total += (previous_f+endpoint_f)*step/2
        error += step**3*M2/12
        previous, previous_f = endpoint, endpoint_f
    lo, hi = ep(total)
    E = ep(error)[1]
    core = _positive(c.mpf([lo-E,hi+E]), 'normalized trapezoid enclosure')
    left_profile, right_profile = kernel_profile(c,Y,a), kernel_profile(c,Y,b)
    left_tail = c.exp(reference-left_profile['H'])/_positive(-left_profile['H1'],'left tail slope')
    right_tail = c.exp(reference-right_profile['H'])/_positive(right_profile['H1'],'right tail slope')
    cutoff = Y*(1-c.exp(c.mpf('-0.4')))
    far_tail = c.exp(reference-cutoff)
    sigma_error = c.exp(-right_profile['D'])
    lower = ep(core/(1+sigma_error))[0]
    upper = ep(core+left_tail+right_tail+far_tail)[1]
    normalized = _positive(c.mpf([lower,upper]), 'full positive normalized turnoff integral')
    log_kernel = c.ln(normalized)-reference
    # The original inlet-memory term is exp(1-Y), retained by log-sum-exp.
    factor = 1+c.exp(1-Y-log_kernel)
    log_F = log_kernel+c.ln(factor)
    relative = c.exp(c.mpf(ep(log_F)[1]-ep(log_F)[0]))-1
    return dict(log_kernel=log_kernel, log_F=log_F, reference_H=reference,
                normalized_integral=normalized, saddle_bracket=root_box,
                center_decimal=center_decimal, window=[a,b], original_Y=Y,
                trapezoid_error=error, normalized_positive_tails=[left_tail,right_tail,far_tail],
                sigma_denominator_error=sigma_error, kernel_relative_width=relative,
                cells=cells, coefficient_precision=c.dps,
                convexity_domain='0<x<=1/100, Md=40',
                original_lower_memory_exp_1_minus_Y_retained=True)


def _graph_expression(nodes, index, bindings):
    if index in bindings:
        return bindings[index]
    data=nodes[index]; op=data['operation']; at=lambda i:_graph_expression(nodes,i,bindings)
    if op=='exact_rational': return s.Rational(data['numerator'],data['denominator'])
    if op=='sum': return s.Add(*(at(i) for i in data['arguments']))
    if op=='product': return s.Mul(*(at(i) for i in data['arguments']))
    if op=='negative': return -at(data['argument'])
    if op=='positive_quotient': return at(data['numerator'])/at(data['denominator'])
    if op=='analytic_unary' and data['name'] in ('exp','log'):
        return getattr(s,data['name'])(at(data['argument']))
    if op=='exact_original_source_flat_sigma':
        return s.Function('original_sigma')(at(data['argument']))
    raise ValueError('Unexpected actual turnoff source node '+str(index))


class CurrentOriginalN1QuietPressureSource:
    @source_precision
    def __init__(self, require_checked=True):
        self.primitives=prior.CurrentOriginalN1QuietPressurePrimitives()
        self.family=copy.deepcopy(self.primitives.family)
        self.hashes=dict(self.primitives.hashes)
        constants=prior.read(CONSTANTS); receipt=prior.read(CONSTANTS_RECEIPT)
        if (not receipt['all_passed'] or not receipt['current_original_Rp_native_cached_constant_function_identity_installed']
                or constants['source_family']!=self.family or receipt['source_family']!=self.family):
            raise ValueError('Actual original native physical-unit identity required')
        prior._merge_hashes(self.hashes,receipt['input_hashes'])
        for name in (prior.NAME,prior.RECEIPT,CONSTANTS,CONSTANTS_RECEIPT,CONSTANTS_SOURCE,Path(__file__).name):
            self.hashes[name]=sha(name)
        self.nodes=copy.deepcopy(self.primitives.frame_raw['exact_graph_nodes'])
        self.functions=copy.deepcopy(self.primitives.frame_raw['actual_Rp_native_frame_functions'])
        kappa=self.primitives.kappa
        symbols={str(q):q for q in kappa.free_symbols}
        Tw,logP=symbols['Tw'],symbols['logP']
        Y=s.exp(40)
        kernel=s.cancel(kappa*s.exp(Tw+logP+12)/4-s.exp(1-Y))
        if not (kernel.is_Function and str(kernel.func).startswith('original_kernel_') and kernel.args==(Y,)):
            raise ValueError('Actual original scalar turnoff kernel not identified')
        self.kernel_node=int(str(kernel.func).rsplit('_',1)[1])
        node=self.nodes[self.kernel_node]
        if node['operation']!='definite_integral' or not node['coordinate_and_kernel_Z_independent']:
            raise ValueError('Original Z-independent mean kernel integral required')
        bound_nodes={i for i,row in enumerate(self.nodes) if row['operation']=='bound_variable'}
        names={self.nodes[i]['name']:i for i in bound_nodes}
        integration=names[node['variable']]
        y,u=s.symbols('y u',positive=True)
        actual=_graph_expression(self.nodes,node['integrand'],{integration:y,node['upper']:u})
        expected=s.exp(y-u)*s.Function('original_sigma')(1-s.log(y)/40)
        if s.simplify(actual-expected)!=0 or _graph_expression(self.nodes,node['lower'],{})!=1:
            raise ValueError('Actual turnoff integrand/endpoint definition differs')
        M=s.sympify(constants['actual_current_native_constants']['M'])
        # Accepted native physical history M=Pstar*kappa, Pstar=exp(logP).
        if s.simplify(M-s.exp(logP)*kappa)!=0:
            raise ValueError('Actual physical Pstar normalization differs')
        self.identity=dict(actual_kernel_node=self.kernel_node,
            original_defining_integrand=s.sstr(actual),exact_cutoff='exp(40)',
            physical_Rc_node=self.functions['original_logRc'],
            physical_K='4*Rc*exp(-14)*(kernel_at_exp40+exp(1-exp40))',
            Pstar_equals_exp_logP_cancelled_before_evaluation=True,
            actual_compact_repaired_quiet_source_consumed=True,
            physical_Rc_is_absolute_not_original_phase_offset=True,
            exact_original_quiet_length_node=self.functions['quiet_length'])
        self._fingerprint=prior._fingerprint((self.nodes,self.functions,self.identity,self.hashes))
        self.acceptance_loaded=False
        if require_checked:
            checked=prior.read(RECEIPT)
            if not checked['all_passed'] or not checked[GATE] or checked['source_family']!=self.family or any(checked[k] for k in OPEN):
                raise ValueError('Actual quiet source-scaled pressure receipt/scope differs')
            prior._merge_hashes(self.hashes,checked['input_hashes'])
            self.acceptance_loaded=True
        self._hashes=copy.deepcopy(self.hashes)
        self._fingerprint=prior._fingerprint((self.nodes,self.functions,self.identity,self.hashes))

    def assert_source(self):
        self.primitives.assert_source()
        if self.hashes!=self._hashes or prior._fingerprint((self.nodes,self.functions,self.identity,self.hashes))!=self._fingerprint:
            raise ValueError('Actual quiet pressure source graph changed')

    @source_precision
    def evaluate(self, cells=2048):
        self.assert_source();c=MPIntervalContext();c.dps=140
        kernel=original_turnoff_kernel(c,cells)
        rows=self.primitives.quiet_raw['actual_four_Z_C3_power_ranges']
        L=read_interval(c,rows[0]['actual_positive_quiet_length'])
        _positive(L,'actual quiet length')
        if any(read_interval(c,row['actual_positive_quiet_length'])._mpi_!=L._mpi_ for row in rows[1:]):
            raise ValueError('One original Z-independent quiet length required')
        logF=kernel['log_F']
        m13=c.ln(2)-28+2*logF+c.ln(-c.expm1(-L))
        m15=c.ln(4)-28+2*logF+c.ln(L)
        log_width=max(ep(m13)[1]-ep(m13)[0],ep(m15)[1]-ep(m15)[0])
        relative=c.exp(c.mpf(log_width))-1
        view=dict(source_family=self.family,source_identity=self.identity,
            original_log_Rc_source_node=self.functions['original_logRc'],
            quiet_segment='actual_2Rc_to_Rp',
            m13=dict(sign='positive',radial_scale_power=1,log_coefficient=record(m13)),
            m15=dict(sign='negative',radial_scale_power=2,log_absolute_coefficient=record(m15)),
            K=dict(sign='positive',radial_scale_power=1,log_coefficient=record(c.ln(4)-14+logF)),
            original_quiet_length=record(L),
            source_scaled_Z1_exact_zero=True,
            normalized_quiet_coefficient_relative_width=record(relative),
            normalized_quiet_coefficient_relative_width_target='0.001',
            normalized_quiet_coefficient_relative_width_passed=ep(relative)[1]<mp.mpf('.001'),
            kernel={k:record(v) if hasattr(v,'_mpi_') else
                    [record(q) for q in v] if isinstance(v,list) else v for k,v in kernel.items()},
            physical_radius_scale_retained_as_original_function=True,
            signed_coefficients_not_magnitude_caps=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        return view


def run(cells=2048):
    began=time.monotonic();owner=CurrentOriginalN1QuietPressureSource(False)
    view=owner.evaluate(cells)
    raw=dict(source_family=owner.family,actual_quiet_source=view,input_hashes=owner.hashes,
             source_identity=owner.identity,execution_seconds=time.monotonic()-began,
             **{GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('QUIET_ACTUAL_SOURCE_SCALED_PRESSURE_COMPUTED',cells,view['normalized_quiet_coefficient_relative_width_passed'],flush=True)
    return owner,view


if __name__=='__main__':run()
