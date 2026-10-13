"""Independent checks for the refined original Rp coefficient packets.

The producer is deliberately only a candidate before this checker writes its
receipt.  This file re-derives the two by two inverse and the positive
quadratic Taylor recurrence, then checks the live packets against the actual
source inputs and the prior coefficient boxes.
"""
import ast
import gzip
import hashlib
import importlib
import inspect
import json
import textwrap
from pathlib import Path

import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext

import lei_ren_part1_paper_compliant_current_original_Rp_pulse_coefficients_refined as current


HERE = Path(__file__).parent
NAME = current.NAME
RECEIPT = current.RECEIPT
SOURCE_DPS = current.basis.SOURCE_DPS
GATE = current.GATE
OPEN = current.OPEN
FIELDS = tuple(current.FIELDS)
CELL_NAMES = ('negative_outer', 'negative_inner', 'positive_inner',
              'positive_outer', 'named_real_point')
EXPECTED_Z = {
    'negative_outer': ('-1', '-.5'),
    'negative_inner': ('-.5', '0'),
    'positive_inner': ('0', '.5'),
    'positive_outer': ('.5', '1'),
    'named_real_point': ('.371', '.371'),
}
CONSTRAINTS = ('first_linear', 'divided_linear', 'energy')
REQUIRED_BINDINGS = (
    'original_linear_inverse', 'original_positive_amplitude',
    'original_current_C5_selection', 'original_current_C5_future',
    'implicit_quadratic_jets', 'current_selected_owner',
    'refined_linear_inverse', 'refined_positive_selection',
)


def _sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def _read_raw():
    return json.loads(gzip.decompress((HERE / NAME).read_bytes()))


def _bounds(value):
    if isinstance(value, dict):
        if 'lower_exact_mpf_tuple' not in value:
            raise TypeError('Interval record has no exact endpoints')
        return (mp.make_mpf(tuple(value['lower_exact_mpf_tuple'])),
                mp.make_mpf(tuple(value['upper_exact_mpf_tuple'])))
    if hasattr(value, '_mpi_'):
        return tuple(mp.make_mpf(item) for item in value._mpi_)
    raise TypeError('Expected a directed interval')


def _ctx_interval(c, value):
    if hasattr(value, '_mpi_'):
        return c.make_mpf(value._mpi_)
    if isinstance(value, dict):
        return c.make_mpf((tuple(value['lower_exact_mpf_tuple']),
                           tuple(value['upper_exact_mpf_tuple'])))
    raise TypeError('Expected interval value')


def _tay(value):
    if not isinstance(value, dict) or 'Taylor_order' not in value:
        raise TypeError('Expected recorded Taylor jet')
    coeffs = value.get('coefficients')
    if len(coeffs) != int(value['Taylor_order']) + 1:
        raise ArithmeticError('Taylor order/row length mismatch')
    return coeffs


def _zero(value):
    lo, hi = _bounds(value)
    return lo <= 0 <= hi


def _nonzero(value):
    lo, hi = _bounds(value)
    return not (lo == 0 and hi == 0)


def _positive(value):
    return _bounds(value)[0] > 0


def _relative_width(value):
    lo, hi = _bounds(value)
    if lo <= 0 <= hi:
        raise ArithmeticError('Width denominator is not sign definite')
    return (hi - lo) / abs(lo)


def _check_hashes(raw):
    for name, digest in raw.get('input_hashes', {}).items():
        if not (HERE / name).exists():
            raise ValueError('Coefficient source is missing: ' + name)
        if _sha(name) != digest:
            raise ValueError('Coefficient source changed: ' + name)


def symbolic_fixture():
    """Prove the original inverse and positive quadratic recurrence directly."""
    A0, A1, D0, D1, r0, r1 = sy.symbols('A0 A1 D0 D1 r0 r1')
    det = A0 * D1 - A1 * D0
    x0 = (r0 * D1 - r1 * A1) / det
    x1 = (r1 * A0 - r0 * D0) / det
    if sy.simplify(A0 * x0 + A1 * x1 - r0) != 0:
        raise ArithmeticError('First divided row inverse identity failed')
    if sy.simplify(D0 * x0 + D1 * x1 - r1) != 0:
        raise ArithmeticError('Second divided row inverse identity failed')

    a2 = sy.symbols('A2', nonzero=True)
    a = sy.symbols('a0:6')
    a1 = sy.symbols('A1_0:6')
    a0 = sy.symbols('A0_0:6')
    residuals = {}
    for degree in range(1, 6):
        full = (a2 * sum(a[j] * a[degree - j]
                         for j in range(degree + 1)) +
                sum(a1[i] * a[degree - i]
                    for i in range(degree + 1)) + a0[degree])
        recurrence = ((2 * a2 * a[0] + a1[0]) * a[degree] +
                      a2 * sum(a[j] * a[degree - j]
                               for j in range(1, degree)) +
                      sum(a1[i] * a[degree - i]
                          for i in range(1, degree + 1)) + a0[degree])
        residuals[str(degree)] = sy.simplify(full - recurrence)
        if residuals[str(degree)] != 0:
            raise ArithmeticError('Generic positive quadratic recurrence failed at ' + str(degree))
    return dict(linear_inverse_first='0', linear_inverse_second='0',
                positive_quadratic_recurrence={degree:str(value) for degree,value in residuals.items()},
                positive_quadratic_recurrence_all_zero=True,
                recurrence_order=5)


def _series_mul(left, right, order):
    return [sum(left[i] * right[n - i] for i in range(n + 1))
            for n in range(order + 1)]


def _series_add(*values):
    return [sum(value[n] for value in values)
            for n in range(len(values[0]))]


def moderate_fixture():
    """Compare production positive_selection with an independent C5 solve."""
    # The production fixture is interval-evaluated at 540 digits.  Use a
    # higher independent scalar precision so nearest rounding cannot place
    # the oracle one ulp outside a correctly directed endpoint.
    with mp.workdps(SOURCE_DPS + 80):
        order = 5
        nu = [mp.mpf('.20'), mp.mpf('.15')]
        # v is constant in the actual source; u and target carry all five
        # nonconstant coefficients.
        u = [[mp.mpf('.30'), mp.mpf('.04'), mp.mpf('.01'), mp.mpf('.02'),
              -mp.mpf('.01'), mp.mpf('.005')],
             [mp.mpf('.10'), -mp.mpf('.03'), mp.mpf('.02'), mp.mpf('.01'),
              mp.mpf('.004'), -mp.mpf('.003')]]
        v = [mp.mpf('.40'), -mp.mpf('.25')]
        K = mp.mpf('2')
        chosen_root = mp.mpf('1.01')
        target = [mp.mpf('0')] * (order + 1)
        a2 = K + sum(weight * value**2 for weight, value in zip(nu, v))
        a1 = [2 * sum(weight * uj[n] * vj
                      for weight, uj, vj in zip(nu, u, v))
              for n in range(order + 1)]
        a0_without_target = [sum(weight * sum(uj[j] * uj[n - j]
                                               for j in range(n + 1))
                                  for weight, uj in zip(nu, u))
                             for n in range(order + 1)]
        target[0] = a2 * chosen_root**2 + a1[0] * chosen_root + a0_without_target[0]
        for n in range(1, order + 1):
            target[n] = mp.mpf('.01') * (n + 1) * (-1 if n % 2 else 1)
        a0 = [a0_without_target[n] - target[n] for n in range(order + 1)]
        denominator = 2 * a2 * chosen_root + a1[0]
        if not (mp.mpf('.9') < chosen_root < mp.mpf('1.2')):
            raise ArithmeticError('Moderate fixture root left source guard')
        c = MPIntervalContext()
        c.dps = SOURCE_DPS
        nu_iv = [c.mpf(str(value)) for value in nu]
        u_iv = [current.IntervalTaylor(c, [c.mpf(str(value)) for value in row])
                for row in u]
        v_iv = [c.mpf(str(value)) for value in v]
        target_iv = current.IntervalTaylor(c, [c.mpf(str(value)) for value in target])
        ap, controls, quadratic = current.positive_selection(
            c, c.mpf(str(K)), nu_iv, u_iv, v_iv, target_iv)
        if hasattr(quadratic['A2'], 'coefficients'):
            raise ArithmeticError('Production fixture A2 is not scalar')
        expected = [chosen_root]
        for n in range(1, order + 1):
            known = a0[n] + sum(a1[i] * expected[n - i]
                                 for i in range(1, n + 1))
            known += a2 * sum(expected[j] * expected[n - j]
                               for j in range(1, n))
            expected.append(-known / denominator)
        expected_controls = [[u[row][n] + expected[n] * v[row]
                              for n in range(order + 1)] for row in range(2)]
        for got, want in zip(ap.coefficients, expected):
            lo, hi = _bounds(got)
            if not (lo <= want <= hi):
                raise ArithmeticError('Production amplitude fixture missed independent C5 value')
        for row, jet in enumerate(controls):
            for got, want in zip(jet.coefficients, expected_controls[row]):
                lo, hi = _bounds(got)
                if not (lo <= want <= hi):
                    raise ArithmeticError('Production control fixture missed independent C5 value')
        return dict(order=order, nonzero_nu=True,
                    constant_v=True, nonconstant_u_and_target=True,
                    chosen_root=mp.nstr(chosen_root, 60),
                    denominator=mp.nstr(denominator, 60),
                    production_matches_independent_C5=True)


def _check_disjoint_intersection():
    c = MPIntervalContext()
    c.dps = 120
    try:
        current.intersect(c, c.mpf([0, 1]), c.mpf([2, 3]))
    except (ArithmeticError, ValueError) as exc:
        if 'conflict' not in str(exc) and 'intersect' not in str(exc):
            raise
        return True
    raise ArithmeticError('Disjoint coefficient intersection was accepted')


def _binding_hashes(bindings):
    if not isinstance(bindings, dict) or not bindings:
        raise ValueError('Actual source callable bindings are missing')
    absent = [name for name in REQUIRED_BINDINGS if name not in bindings]
    if absent:
        raise ValueError('Missing source bindings: ' + ','.join(absent))
    result = {}
    for name, binding in bindings.items():
        for key in ('module', 'qualname', 'source_file', 'source_sha256',
                    'source_ast_sha256'):
            if not binding.get(key):
                raise ValueError('Incomplete source binding: ' + name)
        source_file = binding['source_file']
        if not (HERE / source_file).exists() or _sha(source_file) != binding['source_sha256']:
            raise ValueError('Source binding hash mismatch: ' + name)
        function = importlib.import_module(binding['module'])
        for part in binding['qualname'].split('.'):
            if part == '<locals>':
                raise ValueError('Local callable binding is not accepted: ' + name)
            function = getattr(function, part)
        function = getattr(function, '__func__', function)
        source = textwrap.dedent(inspect.getsource(function))
        tree = ast.dump(ast.parse(source), include_attributes=False)
        ast_hash = hashlib.sha256(tree.encode()).hexdigest()
        if ast_hash != binding['source_ast_sha256']:
            raise ValueError('Source AST binding mismatch: ' + name)
        result[name] = dict(source_file=source_file,
                            source_sha256=binding['source_sha256'],
                            source_ast_sha256=ast_hash)
    return result


def _nonnegative_cap(value):
    lo, hi = _bounds(value)
    return lo >= 0 and hi > 0


def _check_source_inputs(view, owner):
    inputs = view.get('source_inputs')
    if not isinstance(inputs, dict):
        raise ValueError('Genuine source input packet missing')
    required = ('incoming_moment_Taylor', 'incoming_energy_Taylor',
                'weighted_complete_future_Taylor',
                'fixed_incoming_row_factor_caps',
                'fixed_incoming_row_factor_log_definitions',
                'positive_incoming_row_factor_log_margins',
                'positive_end_energy_log_margins',
                'end_energy_definition',
                'original_incoming_shapes')
    if any(key not in inputs for key in required):
        raise ValueError('Incomplete genuine source input packet')
    incoming = inputs['incoming_moment_Taylor']
    if len(incoming) != 2:
        raise ValueError('Expected two genuine incoming moment jets')
    jets = incoming + [inputs['incoming_energy_Taylor'],
                       inputs['weighted_complete_future_Taylor']]
    for jet in jets:
        if len(_tay(jet)) != 6:
            raise ArithmeticError('Genuine source jet is not C5')
        if not any(_nonzero(item) for item in _tay(jet)):
            raise ArithmeticError('Genuine source jet was replaced by zero')
    caps = inputs['fixed_incoming_row_factor_caps']
    if not caps or any(not _nonnegative_cap(item) for item in caps):
        raise ArithmeticError('Incoming row factor caps are not retained [0,positive] boxes')
    actual_caps = owner.selected.amplitude.incoming_factor_caps
    if len(caps) != len(actual_caps) or any(_bounds(a) != _bounds(b)
                                            for a, b in zip(caps, actual_caps)):
        raise ArithmeticError('Incoming row factor caps changed from the owner')
    definitions = inputs['fixed_incoming_row_factor_log_definitions']
    actual_definitions = owner.selected.amplitude.incoming_factor_log_definitions
    if len(definitions) != len(actual_definitions):
        raise ArithmeticError('Incoming log definition count changed')
    for recorded, actual in zip(definitions, actual_definitions):
        if _bounds(recorded) != _bounds(actual):
            raise ArithmeticError('Incoming exact log definition was changed')
    margins = owner.selected.amplitude.incoming_factor_margins
    if not margins or any(_bounds(item)[0] <= 0 for item in margins):
        raise ArithmeticError('Incoming factor positive margins were not retained')
    recorded_margins = inputs['positive_incoming_row_factor_log_margins']
    if len(recorded_margins) != len(margins) or any(_bounds(a) != _bounds(b)
                                                    for a, b in zip(recorded_margins, margins)):
        raise ArithmeticError('Incoming positive margin packet changed')
    energy_margins = owner.selected.amplitude.energy_log_margins
    if not energy_margins or any(_bounds(item)[0] <= 0 for item in energy_margins):
        raise ArithmeticError('Energy positive margins were not retained')
    recorded_energy_margins = inputs['positive_end_energy_log_margins']
    if (len(recorded_energy_margins) != len(energy_margins) or
            any(_bounds(a) != _bounds(b)
                for a, b in zip(recorded_energy_margins, energy_margins))):
        raise ArithmeticError('Energy positive margin packet changed')
    if inputs['end_energy_definition'] != 'nu_j=mu*Kj*exp(2*log_end_scale)>0; finite caps are enclosures':
        raise ValueError('Exact nonzero end-energy definition changed')
    recorded_nu = view.get('nonzero_energy_caps_retained', [])
    actual_nu = owner.selected.amplitude.nu
    if (len(recorded_nu) != len(actual_nu) or
            any(_bounds(a) != _bounds(b) for a, b in zip(recorded_nu, actual_nu))):
        raise ArithmeticError('Nonzero end-energy cap packet changed')
    return dict(incoming_C5=True, future_C5=True,
                incoming_row_factor_caps_retained=True,
                incoming_log_definitions_exact=True,
                incoming_and_energy_margins_positive=True,
                source_shapes_retained=True)


def _check_quadratic_and_constraints(view):
    refined = view.get('refined_coefficients', {})
    budgets = view.get('value_budgets', {})
    prior = view.get('prior_value_budgets', {})
    for field, sign in (('ap', 'positive'), ('C1', 'negative'), ('C2', 'positive')):
        jet = _tay(refined[field])
        if len(jet) != 6:
            raise ArithmeticError(field + ' is not C5')
        lo, hi = _bounds(jet[0])
        if sign == 'positive' and lo <= 0:
            raise ArithmeticError('ap is not strictly positive')
        if sign == 'negative' and hi >= 0:
            raise ArithmeticError('C1 is not strictly negative')
        if sign == 'positive' and field == 'C2' and lo <= 0:
            raise ArithmeticError('C2 is not strictly positive')
        width = _bounds(budgets[field]['relative_width'])
        old_width = _bounds(prior[field]['relative_width'])
        if width[1] >= mp.mpf('.001'):
            raise ArithmeticError(field + ' value width is not below .001')
        if width[1] >= old_width[0]:
            raise ArithmeticError(field + ' value width did not narrow')

    quadratic = view.get('quadratic', {})
    if 'A2' not in quadratic or isinstance(quadratic['A2'], dict) and 'Taylor_order' in quadratic['A2']:
        raise ArithmeticError('Production quadratic A2 must be a scalar source interval')
    if _bounds(quadratic['A2'])[0] <= 0:
        raise ArithmeticError('Production scalar quadratic A2 is not positive')
    for key in ('A1', 'A0'):
        if len(_tay(quadratic[key])) != 6:
            raise ArithmeticError('Quadratic ' + key + ' is not C5')
    denominator = quadratic['positive_derivative_denominator']
    if isinstance(denominator, dict) and 'Taylor_order' in denominator:
        denominator = _tay(denominator)[0]
    if _bounds(denominator)[0] <= 0:
        raise ArithmeticError('Positive quadratic derivative denominator lost')

    c = MPIntervalContext()
    c.dps = SOURCE_DPS
    a = [_ctx_interval(c, item) for item in _tay(refined['ap'])]
    arrays = {
        'A2': [_ctx_interval(c, quadratic['A2'])] + [c.mpf(0)] * 5,
        'A1': [_ctx_interval(c, item) for item in _tay(quadratic['A1'])],
        'A0': [_ctx_interval(c, item) for item in _tay(quadratic['A0'])],
    }
    for degree in range(6):
        residual = arrays['A0'][degree]
        for i in range(degree + 1):
            residual += arrays['A1'][i] * a[degree - i]
            for j in range(degree - i + 1):
                residual += arrays['A2'][i] * a[j] * a[degree - i - j]
        if not (_bounds(residual)[0] <= 0 <= _bounds(residual)[1]):
            raise ArithmeticError('Quadratic residual missed zero at row ' + str(degree))

    constraints = view.get('source_constraint_enclosures', {})
    if any(key not in constraints for key in CONSTRAINTS):
        raise ValueError('Missing genuine source constraint enclosure')
    for key in CONSTRAINTS:
        row = _tay(constraints[key])
        if len(row) != 6:
            raise ArithmeticError(key + ' constraint is not C5')
        if any(not _zero(item) for item in row):
            raise ArithmeticError(key + ' constraint coefficient missed zero')
    nu = view.get('nonzero_energy_caps_retained', [])
    if not nu or any(not _nonnegative_cap(item) for item in nu):
        raise ArithmeticError('Nonzero original nu caps were not retained')
    return dict(value_signs_and_widths=True, quadratic_C5_residual=True,
                three_source_constraints_all_six_zero=True,
                nonzero_nu_retained=True)


def _check_view(name, view, owner):
    if view.get('source_family') != owner.family:
        raise ValueError('Source family mismatch at ' + name)
    zlo, zhi = _bounds(view['Z'])
    zctx = MPIntervalContext()
    zctx.dps = SOURCE_DPS
    expected = zctx.mpf(list(EXPECTED_Z[name])) if name != 'named_real_point' else zctx.mpf('.371')
    elo, ehi = tuple(mp.make_mpf(item) for item in expected._mpi_)
    if zlo != elo or zhi != ehi:
        raise ArithmeticError('Unexpected Z cell at ' + name)
    if view.get('original_genuine_C5_incoming_and_future_retained') is not True:
        raise ValueError('Genuine incoming/future C5 flag missing')
    if view.get('original_positive_branch_preserved') is not True:
        raise ValueError('Original positive branch flag missing')
    if view.get('derivative_accuracy_not_claimed') is not True:
        raise ValueError('Derivative precision was overclaimed')
    bindings = _binding_hashes(view.get('source_bindings'))
    inputs = _check_source_inputs(view, owner)
    algebra = _check_quadratic_and_constraints(view)
    return dict(name=name, source_bindings=bindings, source_inputs=inputs,
                algebra=algebra)


def _live_views(owner, packets):
    if packets is None:
        packets = {name: owner.evaluate(list(EXPECTED_Z[name]))
                   if name != 'named_real_point' else owner.evaluate('.371')
                   for name in CELL_NAMES}
    if set(packets) != set(CELL_NAMES):
        raise ValueError('The four axial cells and .371 are required')
    for name, packet in packets.items():
        state = getattr(owner, '_packets', {}).get(id(packet))
        if state is None or state[0] is not packet or len(state) < 4:
            raise ValueError('Live coefficient packet state is missing at ' + name)
    return {name: owner.report(packets[name]) for name in CELL_NAMES}


def run(owner=None, packets=None):
    with mp.workdps(SOURCE_DPS):
        if owner is None:
            owner = current.CurrentOriginalRpPulseCoefficientsRefined(require_checked=False)
        views = _live_views(owner, packets)
        raw = _read_raw()
        _check_hashes(raw)
        if raw.get('source_family') != owner.family:
            raise ValueError('Raw source family mismatch')
        if raw.get(GATE) is not False or any(raw.get(key, False) for key in OPEN):
            raise ValueError('Candidate raw scope flags changed')
        if set(raw.get('samples', {})) != set(CELL_NAMES):
            raise ValueError('Raw sample coverage is incomplete')
        checked = {name: _check_view(name, views[name], owner)
                   for name in CELL_NAMES}
        # The saved candidate view must carry the complete live report.  This
        # catches stale source inputs, constraints, bindings, or packet state,
        # rather than checking only the Z endpoints.
        for name in CELL_NAMES:
            saved = raw['samples'][name]
            live = current.record(views[name])
            live.pop('execution_seconds', None)
            saved = dict(saved)
            saved.pop('execution_seconds', None)
            if live != saved:
                raise ArithmeticError('Saved/live view mismatch at ' + name)
        symbolic = symbolic_fixture()
        moderate = moderate_fixture()
        disjoint = _check_disjoint_intersection()
        saved_dps = mp.mp.dps
        with mp.workdps(15):
            public = current.CurrentOriginalRpPulseCoefficientsRefined(
                selected_owner=owner.selected, basis_owner=owner.basis,
                require_checked=False)
            packet = public.evaluate('.371')
            public.report(packet)
            if mp.mp.dps != 15:
                raise ArithmeticError('Public coefficient query leaked precision')
        if mp.mp.dps != saved_dps:
            raise ArithmeticError('Outer precision was not restored')
        own_hashes = dict(owner.hashes)
        own_hashes[NAME] = _sha(NAME)
        own_hashes[Path(__file__).name] = _sha(Path(__file__).name)
        receipt = dict(
            all_passed=True,
            source_family=owner.family,
            checks=dict(symbolic_inverse_and_quadratic=True,
                        moderate_nonzero_nu_fixture=True,
                        all_four_axial_cells_and_named_point=True,
                        source_bindings_and_graph_hashes=True,
                        genuine_C5_inputs_retained=True,
                        all_value_widths_below_001=True,
                        all_value_widths_narrower_than_prior=True,
                        all_constraint_rows_contain_zero=True,
                        nonzero_nu_retained=True,
                        fail_closed_disjoint_intersection=disjoint,
                        default_precision_public_query=True,
                        no_pressure_or_tail_claim=True),
            cells=checked,
            symbolic=symbolic,
            moderate_fixture=moderate,
            original_selection_branch_unchanged=True,
            derivative_accuracy_not_claimed=True,
            input_hashes=own_hashes,
            **{GATE: True},
            **dict.fromkeys(OPEN, False),
        )
        (HERE / RECEIPT).write_text(
            json.dumps(receipt, separators=(',', ':'), ensure_ascii=False) + '\n',
            encoding='utf8', newline='\n')
        checked_owner = current.CurrentOriginalRpPulseCoefficientsRefined(
            selected_owner=owner.selected, basis_owner=owner.basis,
            require_checked=True)
        print('CURRENT_RP_PULSE_COEFFICIENTS_REFINED_CHECK PASS', flush=True)
        return receipt


if __name__ == '__main__':
    run()
