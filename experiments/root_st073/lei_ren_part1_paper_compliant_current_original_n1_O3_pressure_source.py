"""Actual signed O3 pressure values, retaining the finite-N incoming mean.

U is exactly zero on the original O3 power segment. Its incoming mean is
nonzero and decays with the original own rate one. Exact radial integration
and the accepted scalar baseline kernel give signed source-scaled values;
the finite-N perturbation is enclosed rather than silently discarded.
The third mean derivative remains coarse, so Z1 accuracy is not installed.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s

import lei_ren_part1_paper_compliant_current_original_n1_quiet_pressure_source as quiet
import lei_ren_part1_paper_compliant_current_original_n1_all_chart_physical_pressure_fields as physical

atlas, partial, rg = physical.all_chart, physical.partial, physical.rg
HERE, PREFIX, sha = quiet.HERE, quiet.PREFIX, quiet.sha
source_precision, record, read_interval, ep = quiet.source_precision, quiet.record, quiet.read_interval, quiet.ep
NAME = PREFIX + 'current_original_n1_O3_pressure_source.json.gz'
RECEIPT = PREFIX + 'current_original_n1_O3_pressure_source_check.json'
TARGET_NAME = PREFIX + 'current_original_C3_source_targets.json.gz'
TARGET_RECEIPT = PREFIX + 'current_original_C3_source_targets_check.json'
GATE = 'current_original_n1_O3_source_scaled_signed_pressure_values_installed'
DERIVATIVE_GATE = 'current_original_n1_O3_source_scaled_Z1_accuracy_installed'
OPEN = ('complete_signed_pressure_tails_certified', 'positive_order_five_moment_repair_closed',
        'temporal_recursion', 'regional_admissible_stress_certified',
        'flat_remainder_certified', 'full_corrected_Navier_Stokes_solution',
        'absolute_physical_radial_scale_materialization_installed')


def general_zero_axial_primitives(B, B_Z, Z, delta, a, b):
    """Exact physical dR / R*dR integrals for q=B(Z)/R, U=0."""
    B, B_Z, Z, delta, a, b = map(s.sympify, (B, B_Z, Z, delta, a, b))
    if a.is_number and b.is_number and not (bool(a > 0) and bool(b > a)):
        raise ValueError('Positive ordered physical radial endpoints required')
    g = (1-delta)*Z/(4*(1-delta*Z**2))
    return dict(m13=B**2*(1/a-1/b)/4-g*B_Z*s.log(b/a),
                m15=-B**2*s.log(b/a)/4+g*B_Z*(b-a))


def normalized_pressure_functions(Z, delta, btilde):
    """Actual O3 segment [Rc*exp(-2),Rc], normalized by Rc and Rc²."""
    Z, delta, btilde = map(s.sympify, (Z, delta, btilde))
    g = (1-delta)*Z/(4*(1-delta*Z**2))
    B_Z = s.diff(btilde, Z)
    return dict(m13=(s.exp(2)-1)*btilde**2/4-2*g*B_Z,
                m15=-btilde**2/2+(1-s.exp(-2))*g*B_Z)


def _read(name):
    return quiet.prior.read(name)


def _cap_record(bound):
    return dict(exact_zero=bound.log is None,
                log_absolute_upper=None if bound.log is None else record(bound.log),
                source_value_or_exponential_not_materialized=True)


def _read_cap(c, row):
    return rg.LogUpper(c) if row['exact_zero'] else rg.LogUpper(c, read_interval(c, row['log_absolute_upper']))


def _graph_expression(nodes, index, bindings=None):
    bindings = dict(bindings or {})
    for i, row in enumerate(nodes):
        if row['operation'] in ('bound_variable', 'current_original_source_parameter'):
            bindings.setdefault(i, s.Symbol(row['name'], real=True))
    return quiet._graph_expression(nodes, index, bindings)


def leading_baseline_expression(nodes, index):
    """Interpret the actual leading mean, keeping original kernels exact."""
    memo = {}
    def at(i):
        if i in memo:
            return memo[i]
        row = nodes[i]; operation = row['operation']
        if operation == 'exact_rational':
            value = s.Rational(row['numerator'],row['denominator'])
        elif operation == 'bound_variable':
            value = s.Symbol(row['name'],real=True)
        elif operation == 'sum':
            value = s.Add(*(at(j) for j in row['arguments']))
        elif operation == 'product':
            value = s.Mul(*(at(j) for j in row['arguments']))
        elif operation == 'negative':
            value = -at(row['argument'])
        elif operation == 'positive_quotient':
            value = at(row['numerator'])/at(row['denominator'])
        elif operation == 'analytic_unary' and row['name'] in ('exp','log'):
            value = getattr(s,row['name'])(at(row['argument']))
        elif operation == 'function_substitution' and row['Z_independent_substitution']:
            value = at(row['expression']).subs(at(row['variable']),at(row['value']))
        elif operation == 'definite_integral' and row.get('original_kernel_function') and row.get('coordinate_and_kernel_Z_independent'):
            upper = at(row['upper'])
            value = s.Integer(0) if at(row['lower'])==upper else s.Function('original_kernel_'+str(i))(upper)
        else:
            raise ValueError('Unsupported exact O3 leading mean source '+str(i))
        memo[i] = value
        return value
    return at(index)


class CurrentOriginalN1O3PressureSource:
    @source_precision
    def __init__(self, require_checked=True):
        self.hashes = {}
        specs = (
            ('quiet', quiet.NAME, quiet.RECEIPT, (quiet.GATE,)),
            ('targets', TARGET_NAME, TARGET_RECEIPT,
             ('current_original_actual_17_chart_C3_source_and_same_inverse_functions_installed',
              'current_original_actual_17_chart_C3_density_transport_and_five_target_functions_installed')),
            ('atlas', atlas.NAME, atlas.RECEIPT, (atlas.GATE,)),
            ('physical', physical.NAME, physical.RECEIPT, (physical.GATE,)),
            ('ranges', rg.NAME, rg.RECEIPT, rg.GATES),
            ('frame', quiet.prior.FRAME_NAME, quiet.prior.FRAME_RECEIPT, quiet.prior.FRAME_GATES),
        )
        self.sources = {}
        self.direct_hashes = {}
        for key, name, receipt_name, gates in specs:
            raw, checked = _read(name), _read(receipt_name)
            if not checked['all_passed'] or not all(checked[gate] for gate in gates):
                raise ValueError('Accepted original O3 dependency required: '+receipt_name)
            family = checked['source_family']
            if raw['source_family'] != family or (self.sources and family != self.family):
                raise ValueError('One actual original source family required')
            self.family = copy.deepcopy(family)
            for data in (raw, checked):
                for item, digest in data['input_hashes'].items():
                    if item in self.hashes and self.hashes[item] != digest:
                        raise ValueError('Conflicting original O3 dependency '+item)
                    self.hashes[item] = digest
            for item in (name, receipt_name):
                self.hashes[item] = self.direct_hashes[item] = sha(item)
            self.sources[key] = raw
        self.hashes[Path(__file__).name] = sha(Path(__file__).name)
        for item, digest in self.hashes.items():
            if sha(item) != digest:
                raise ValueError('Changed accepted original O3 input '+item)

        targets = self.sources['targets']
        nodes = targets['exact_graph_nodes']
        frame = self.sources['frame']
        frame_nodes = frame['exact_graph_nodes']
        if frame_nodes[:len(nodes)] != nodes:
            raise ValueError('Same original C3 graph prefix required')
        self.window = next(row for row in targets['actual_C3_density_transport_targets']['windows']
                           if row['chart'] == 'O3_power')
        self.incoming_nodes = self.window['incoming']['m']
        if self.window['density']['m'] != [0,0,0,0] or self.window['contributions']['m'] != [0,0,0,0]:
            raise ValueError('Exact original zero O3 local mean density required')
        if any(i == 0 for i in self.incoming_nodes):
            raise ValueError('Retain the actual nonzero predecessor mean')
        cells = self.window['original_physical_measure_and_partitions']
        if [(row['exact_left'], row['exact_right']) for row in cells] != [([0,1],[1,1]),([1,1],[2,1])]:
            raise ValueError('Original complete O3 native union required')
        left_node, right_node = cells[0]['left_offset'], cells[-1]['right_offset']
        native = next(i for i,row in enumerate(nodes) if row.get('name') == 'native_O3_power')
        offset_node = next(i for i,row in enumerate(nodes)
                           if row['operation'] == 'sum' and native in row['arguments'])
        t = s.Symbol('t', real=True)
        left, right = _graph_expression(nodes,left_node), _graph_expression(nodes,right_node)
        offset = _graph_expression(nodes,offset_node,{native:t})
        if s.simplify(right-left) != 2 or s.simplify(offset-right-(t-2)) != 0:
            raise ValueError('Actual R/Rc=exp(t-2) geometry required')
        native_functions = frame['actual_Rp_native_frame_functions']
        if native_functions['original_Rc_phase_offset'] != right_node:
            raise ValueError('Actual absolute Rc must use the original O3 right offset')
        Rc_node = native_functions['original_logRc']
        origin_node = native_functions['original_log_r_minus']
        if set(frame_nodes[Rc_node]['arguments']) != {right_node,origin_node}:
            raise ValueError('Absolute original Rc origin must be retained')
        for j, (before, after) in enumerate(zip(self.incoming_nodes,self.window['outgoing']['m'])):
            H = s.Symbol('incoming_%d'%j)
            if s.simplify(_graph_expression(nodes,after,{before:H})-s.exp(-2)*H) != 0:
                raise ValueError('Original rate-one predecessor memory required')

        sample = self.sources['physical']['samples']['O3_power']
        if sample['exact_native_cell'] != [[0,1],[2,1]] or sample['exact_Z_cell'] != ['0','.5']:
            raise ValueError('Accepted whole O3 U=0 query required')
        if not all(row['exact_zero'] for row in sample['whole_cell_physical_U_Z3_magnitude']):
            raise ValueError('Actual physical U through Z3 must be exactly zero')
        self.inlet = self.sources['atlas']['samples']['O3_power']
        if self.inlet['exact_Z_cell'] != ['0','.5'] or not self.inlet['same_selected_N_division_once_no_integer_materialization']:
            raise ValueError('Same original finite-N incoming cap required')
        self.whole_Z_transports = {}
        for report in self.sources['ranges']['actual_four_Z_C3_target_transports']:
            ends = tuple(report['exact_Z_cell'])
            transport = next(row for row in report['actual_17_chart_C3_transports'] if row['chart']=='O3_power')
            if report['source_identity'] != self.family:
                raise ValueError('Actual whole-Z O3 ranges use a different source family')
            if [(row['exact_left'],row['exact_right']) for row in transport['cells']] != [([0,1],[1,1]),([1,1],[2,1])]:
                raise ValueError('Whole-Z O3 ranges must cover both original native cells')
            for row in transport['cells']:
                quantities = row['actual_source_and_primitive_range']
                if not all(value['exact_zero'] for name in ('V','B') for value in quantities[name]):
                    raise ValueError('Original U and local density must vanish on the entire O3 cover')
            self.whole_Z_transports[ends] = transport
        cover = [('-1','-.5'),('-.5','0'),('0','.5'),('.5','1')]
        if list(self.whole_Z_transports) != cover:
            raise ValueError('Exact original four axial cells required')
        if self.whole_Z_transports[('0','.5')]['actual_N_scaled_incoming_C3']['m'] != self.inlet['original_incoming_N_scaled_mean_Z3']:
            raise ValueError('All-chart incoming history must equal its original C3 source ranges')
        N_node, N_recipe = next((i,row) for i,row in enumerate(nodes) if row.get('type')=='DoubleDyadicInteger')
        if not N_recipe['one_common_N_for_all_charts_and_control_denominators']:
            raise ValueError('Same selected original N required')
        Z, delta = s.symbols('Z delta', real=True)
        N = s.Symbol('original_selected_N', positive=True, integer=True)
        H = s.Function('original_O3_incoming_mean_%d'%self.incoming_nodes[0])(Z)
        source_kernel_node = self.sources['quiet']['source_identity']['actual_kernel_node']
        F = s.Function('original_kernel_%d'%source_kernel_node)(s.exp(40))+s.exp(1-s.exp(40))
        k = 4*s.exp(-14)*F
        leading2_node = native_functions['actual_leading_2Rc_histories_C3']['m'][0]
        leading2_mean = leading_baseline_expression(frame_nodes,leading2_node)
        # The original leading continuation has zero U from Rc through 2Rc.
        # Thus 2*m(2Rc)=m(Rc). Check its actual defining graph against k,
        # rather than assuming a repaired quiet coefficient is the inlet.
        leadingRc_mean = 2*leading2_mean
        if s.simplify(s.exp(s.exp(40)+11)*leadingRc_mean-k*Z) != 0:
            raise ValueError('Actual unmodified O3 mean does not match the original scalar kernel')
        D = s.exp(s.exp(40)+9)*H/N
        delta_b = (-(1-delta)*Z*D-(1-Z**2)*s.diff(D,Z))/(1-delta*Z**2)
        btilde = -k+delta_b
        normalized = normalized_pressure_functions(Z,delta,btilde)
        self.functions = dict(btilde_Z0_to_Z2=[s.sstr(s.diff(btilde,Z,j)) for j in range(3)],
            normalized_m13_Z0_Z1=[s.sstr(normalized['m13']),s.sstr(s.diff(normalized['m13'],Z))],
            normalized_m15_Z0_Z1=[s.sstr(normalized['m15']),s.sstr(s.diff(normalized['m15'],Z))],
            incoming_mean_symbol=s.sstr(H),incoming_ordinary_Z3_source_nodes=self.incoming_nodes,
            incoming_source_graph_namespace=TARGET_NAME,baseline_kernel_graph_namespace=quiet.prior.FRAME_NAME,
            baseline_kernel_source_node=source_kernel_node,original_selected_N_source_node=N_node,
            original_delta_recipe='exp(-4*(exp(40)+11)-30)',
            physical_Pstar_recipe='exp(exp(40)+11)',actual_incoming_coefficient='exp(-2)*Pstar*H_in/N',
            actual_signed_functions_not_caps_or_selected_centers=True)
        self.identity = dict(chart='O3_power',exact_Z_cell=['0','.5'],exact_native_union=[[0,1],[2,1]],
            exact_Z_cover=[list(ends) for ends in cover],whole_Z_zero_U_and_local_density_from_complete_original_cells=True,
            incoming_node_rows=list(self.incoming_nodes),local_density_node_rows=[0,0,0,0],
            outgoing_node_rows=self.window['outgoing']['m'],geometry_dlogR_width='2',
            original_Rc_phase_offset_shared=True,original_Rc_phase_offset_node=right_node,
            original_absolute_logRc_node=Rc_node,original_log_r_minus_node=origin_node,
            radial_segment='[Rc*exp(-2),Rc]',own_mean_rate_one_and_zero_local_source=True,
            original_U_Z3_exact_zero=True,original_selected_N_division_once=True,
            original_selected_N_node=N_node,original_selected_N_recipe=copy.deepcopy(N_recipe),
            physical_Pstar_applied_once=True,finite_N_incoming_mean_not_reset=True,
            baseline_K_from_same_original_leading_power_and_quiet_source=True,
            actual_leading_2Rc_mean_source_node=leading2_node,
            exact_leading_Rc_common_mean=s.sstr(leadingRc_mean),
            physical_leading_Rc_mean_equals_k_Z=True,
            baseline_source_identity=copy.deepcopy(self.sources['quiet']['source_identity']),
            physical_radius_and_N_not_materialized=True)
        self.acceptance_loaded = False
        if require_checked:
            checked = _read(RECEIPT)
            if not checked['all_passed'] or not checked[GATE] or checked[DERIVATIVE_GATE] or checked['source_family'] != self.family or any(checked[k] for k in OPEN):
                raise ValueError('Original O3 pressure receipt/scope differs')
            for item,digest in checked['input_hashes'].items():
                if sha(item) != digest:
                    raise ValueError('Changed original O3 pressure input '+item)
                self.hashes[item] = digest
            self.acceptance_loaded = True
        self._hashes = copy.deepcopy(self.hashes)
        self._fingerprint = quiet.prior._fingerprint((self.sources,self.identity,self.functions,self.family))

    def assert_source(self):
        if self.hashes != self._hashes or quiet.prior._fingerprint((self.sources,self.identity,self.functions,self.family)) != self._fingerprint:
            raise ValueError('Original O3 pressure source changed')
        for item,digest in self.direct_hashes.items():
            if sha(item) != digest:
                raise ValueError('Changed direct original O3 source '+item)

    @source_precision
    def evaluate(self, ends=('0','.5')):
        self.assert_source()
        ends = tuple(ends)
        if ends not in self.whole_Z_transports:
            raise ValueError('Exact accepted original O3 axial cell required')
        c = MPIntervalContext(); c.dps = 180
        bd = rg.Bounds(c)
        incoming_rows = self.whole_Z_transports[ends]['actual_N_scaled_incoming_C3']['m']
        incoming = rg.JetBound(*(_read_cap(c,row) for row in incoming_rows))
        epsilon_N = _read_cap(c,self.inlet['epsilon_N_log_bound'])
        logN = read_interval(c,self.inlet['selected_N_log'])
        # Check the exact stored binary tuples. Negating a gigantic mpf after
        # lowering working precision can otherwise falsely reject equality.
        expected_negative_lower = list(self.inlet['selected_N_log']['lower_exact_mpf_tuple'])
        expected_negative_lower[0] = 1-expected_negative_lower[0]
        saved_epsilon = self.inlet['epsilon_N_log_bound']['log_absolute_upper']
        if (saved_epsilon['lower_exact_mpf_tuple'] != expected_negative_lower or
                saved_epsilon['upper_exact_mpf_tuple'] != expected_negative_lower):
            raise ValueError('Directed epsilon_N bound must retain actual selected N')
        selected = partial.current.DoubleDyadicInteger(self.identity['original_selected_N_recipe']['J']).log_interval(c)
        if not (ep(selected)[0] <= ep(logN)[1] and ep(logN)[0] <= ep(selected)[1]):
            raise ValueError('Saved incoming cap has a different selected N')
        D = bd.scaled(incoming,epsilon_N*rg.LogUpper(c,c.exp(40)+9))
        zero = bd.fixed(bd.zero)
        delta = c.exp(-4*(c.exp(40)+11)-30)
        b0,b1,b2 = partial.coefficient_Q_bounds(c,ends,delta,zero,D)
        logk = read_interval(c,self.sources['quiet']['actual_quiet_source']['K']['log_coefficient'])
        k = rg.LogUpper(c,logk)
        from lei_ren_part1_paper_interval_taylor import IntervalTaylor
        z = IntervalTaylor.variable(c,c.mpf(list(ends)),1)
        g = z*(1-delta)/((1-z*z*delta)*4)
        g0,g1 = (bd.constant(ep(abs(g[j]))[1]) for j in range(2))
        quadratic = bd.sum(bd.scale(k*b0,2),b0*b0)
        A, f = bd.constant((c.exp(2)-1)/4),bd.constant(1-c.exp(-2))
        E13 = bd.sum(A*quadratic,bd.scale(g0*b1,2))
        E15 = bd.sum(bd.scale(quadratic,c.mpf('.5')),f*g0*b1)
        linear = bd.sum(k,b0)*b1
        gderivative = bd.sum(g1*b1,g0*b2)
        E13Z = bd.sum(bd.scale(A*linear,2),bd.scale(gderivative,2))
        E15Z = bd.sum(linear,f*gderivative)
        baseline = dict(m13=2*logk+c.ln((c.exp(2)-1)/4),m15=2*logk-c.ln(2))
        errors = dict(m13=E13,m15=E15)
        relative_tail = c.mpf('1e-60')
        values = {}
        proofs = {}
        for name, logarithm in baseline.items():
            error = errors[name]
            difference = error.log-c.mpf(ep(logarithm)[0])
            if ep(difference)[1] >= ep(c.ln(relative_tail))[0]:
                raise ValueError('Actual finite-N pressure value correction is not controlled')
            lo,hi = ep(logarithm)
            directed = c.mpf([lo,hi])+c.mpf([ep(c.ln(1-relative_tail))[0],ep(c.ln(1+relative_tail))[1]])
            values[name] = dict(sign='positive' if name=='m13' else 'negative',
                radial_scale_power=1 if name=='m13' else 2,
                log_coefficient=record(directed),log_absolute_coefficient=record(directed))
            proofs[name] = dict(true_log_error=_cap_record(error),
                log_error_relative_to_baseline_lower=record(difference),
                nonzero_directed_relative_tail_upper='1e-60',
                exact_zero_not_substituted_for_small_source=True)
        coefficient_intervals = [read_interval(c,row['log_coefficient']) for row in values.values()]
        width = max(ep(value)[1]-ep(value)[0] for value in coefficient_intervals)
        relative_width = c.expm1(c.mpf(width))
        identity = {**self.identity,'exact_Z_cell':list(ends)}
        view = dict(source_family=self.family,source_identity=identity,functions=self.functions,exact_Z_cell=list(ends),
            m13=values['m13'],m15=values['m15'],baseline_log_coefficients={k:record(v) for k,v in baseline.items()},
            original_incoming_N_scaled_mean_Z3=[_cap_record(v) for v in rg.rows(incoming)],
            physical_incoming_D_Z3_magnitudes=[_cap_record(v) for v in rg.rows(D)],
            source_scaled_delta_b_Z2_magnitudes=[_cap_record(v) for v in (b0,b1,b2)],
            finite_N_incoming_mean_retained=True,value_error_proofs=proofs,
            value_error_relative_to_baseline_upper='1e-60',
            source_scaled_Z1_error_magnitudes=dict(m13=_cap_record(E13Z),m15=_cap_record(E15Z)),
            source_scaled_Z1_accuracy=False,
            normalized_coefficient_relative_width=record(relative_width),
            normalized_coefficient_relative_width_target='0.001',
            normalized_coefficient_relative_width_passed=ep(relative_width)[1]<mp.mpf('.001'),
            interval_scope='all Z in [%s,%s], full actual O3_power segment'%ends,
            physical_radius_scale_retained_as_original_function=True,
            **{GATE:self.acceptance_loaded,DERIVATIVE_GATE:False},**dict.fromkeys(OPEN,False))
        if not view['normalized_coefficient_relative_width_passed']:
            raise ValueError('Actual O3 normalized signed values must meet stated width')
        return view

    def evaluate_all(self):
        return [self.evaluate(ends) for ends in self.whole_Z_transports]


def run():
    began = time.monotonic()
    owner = CurrentOriginalN1O3PressureSource(False)
    owner.evaluated_cells = owner.evaluate_all()
    view = next(row for row in owner.evaluated_cells if row['exact_Z_cell']==['0','.5'])
    raw = dict(source_family=owner.family,actual_O3_pressure_source=view,
               actual_four_Z_O3_pressure_values=owner.evaluated_cells,
               source_identity=owner.identity,functions=owner.functions,input_hashes=owner.hashes,
               execution_seconds=time.monotonic()-began,
               **{GATE:False,DERIVATIVE_GATE:False},**dict.fromkeys(OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_O3_SIGNED_PRESSURE_VALUES_COMPUTED',len(owner.evaluated_cells),
          all(row['normalized_coefficient_relative_width_passed'] for row in owner.evaluated_cells),flush=True)
    return owner,view


if __name__ == '__main__':
    run()
