"""Exact signed pressure primitives on the actual quiet 2 Rc to Rp segment.

The accepted native-frame and repaired quiet artifacts identify one
Z-independent ``kappa``.  This adapter keeps ``K = Rp * Pstar * kappa``
formal and derives the physical quiet fields and their two pressure moments;
it does not evaluate a radius, Pstar, N, or a numerical defect.
"""
import copy
import gzip
import hashlib
import json
from pathlib import Path

import sympy as s


HERE = Path(__file__).resolve().parent
PREFIX = 'lei_ren_part1_paper_compliant_'
FRAME_NAME = PREFIX + 'current_original_C3_Rp_native_frame.json.gz'
FRAME_RECEIPT = PREFIX + 'current_original_C3_Rp_native_frame_check.json'
FRAME_SOURCE = PREFIX + 'current_original_C3_Rp_native_frame.py'
FRAME_CHECK = PREFIX + 'current_original_C3_Rp_native_frame_check.py'
QUIET_NAME = PREFIX + 'current_original_C3_power_to_Rp.json.gz'
QUIET_RECEIPT = PREFIX + 'current_original_C3_power_to_Rp_check.json'
QUIET_SOURCE = PREFIX + 'current_original_C3_power_to_Rp.py'
QUIET_CHECK = PREFIX + 'current_original_C3_power_to_Rp_check.py'

NAME = PREFIX + 'current_original_n1_quiet_pressure_primitives.json'
RECEIPT = PREFIX + 'current_original_n1_quiet_pressure_primitives_check.json'
GATE = 'current_original_n1_quiet_2Rc_to_Rp_pressure_primitives_installed'
OPEN = ('actual_numeric_defect_evaluation_not_installed',
        'complete_pressure_tails_and_n1_repair_not_installed',
        'global_physical_solution_not_installed')
FRAME_GATES = (
    'current_original_complete_C3_chain_to_Rp_native_frame_identity_installed',
    'current_original_Rp_native_true_Z0_to_Z5_inlet_functions_installed',
    'current_original_Rp_native_exact_radius_identity_installed')
QUIET_GATES = (
    'current_original_actual_repaired_exit_C3_quiet_power_and_pulse_frame_installed',
    'current_original_quiet_y4_Z3_profile_history_functions_installed',
    'current_original_quiet_C3_directed_magnitude_and_positive_swirl_ranges_installed')


def sha(name):
    path = Path(name)
    if not path.is_absolute():
        path = HERE / path
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(name):
    data = (HERE / name).read_bytes()
    return json.loads(gzip.decompress(data) if name.endswith('.gz') else data)


def _merge_hashes(out, incoming):
    for name, digest in incoming.items():
        if name in out and out[name] != digest:
            raise ValueError('Conflicting accepted quiet-pressure input ' + name)
        if sha(name) != digest:
            raise ValueError('Changed accepted quiet-pressure input ' + name)
        out[name] = digest


def _fingerprint(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'))


def _formal_symbols():
    return dict(Z=s.Symbol('Z', real=True), R=s.Symbol('R', positive=True),
                Rp=s.Symbol('R_p', positive=True), Pstar=s.Symbol('Pstar', positive=True),
                delta=s.Symbol('delta', real=True))


def _sympy_frame(frame):
    z = _formal_symbols()['Z']
    return {key: s.sympify(value, locals={'Z': z}) for key, value in frame.items()}


class CurrentOriginalN1QuietPressurePrimitives:
    def __init__(self, require_checked=True):
        self.frame_raw = read(FRAME_NAME)
        self.frame_receipt = read(FRAME_RECEIPT)
        self.quiet_raw = read(QUIET_NAME)
        self.quiet_receipt = read(QUIET_RECEIPT)
        self.family = copy.deepcopy(self.frame_raw['source_family'])
        if (self.frame_receipt['source_family'] != self.family or
                self.quiet_raw['source_family'] != self.family or
                self.quiet_receipt['source_family'] != self.family):
            raise ValueError('Accepted quiet/native artifacts have different source families')
        if (not self.frame_receipt['all_passed'] or
                not all(self.frame_receipt[key] for key in FRAME_GATES) or
                not self.quiet_receipt['all_passed'] or
                not all(self.quiet_receipt[key] for key in QUIET_GATES)):
            raise ValueError('Accepted quiet/native receipts required')
        self.hashes = {}
        for data in (self.frame_raw, self.frame_receipt,
                     self.quiet_raw, self.quiet_receipt):
            _merge_hashes(self.hashes, data.get('input_hashes', {}))
        for name in (FRAME_NAME, FRAME_RECEIPT, FRAME_SOURCE, FRAME_CHECK,
                     QUIET_NAME, QUIET_RECEIPT, QUIET_SOURCE, QUIET_CHECK,
                     Path(__file__).name):
            self.hashes[name] = sha(name)
        self.frame_strings = copy.deepcopy(
            self.frame_raw['actual_analytic_native_inlet_functions'])
        self.endpoint_nodes = dict(
            quiet=copy.deepcopy(self.quiet_raw['actual_C3_repaired_exit_power_and_pulse_functions']),
            native=copy.deepcopy(self.frame_raw['actual_Rp_native_frame_functions']))
        self.symbols = _formal_symbols()
        self.Z, self.R = self.symbols['Z'], self.symbols['R']
        self.Rp, self.Pstar = self.symbols['Rp'], self.symbols['Pstar']
        self.delta = self.symbols['delta']
        parsed = _sympy_frame(self.frame_strings)
        self.mean_common = s.cancel(parsed['u'] * parsed['m1'])
        self.kappa = s.cancel(self.mean_common / self.Z)
        self.q_common_U0 = s.cancel(
            (-(1 - self.delta) * self.Z * self.mean_common -
             (1 - self.Z**2) * s.diff(self.mean_common, self.Z)) /
            (1 - self.delta * self.Z**2))
        self.K_definition = s.Mul(self.Rp, self.Pstar, self.kappa)
        self._raw_fingerprint = _fingerprint((self.frame_raw, self.quiet_raw))
        self.acceptance_loaded = False
        if require_checked:
            checked = read(RECEIPT)
            if (not checked['all_passed'] or not checked[GATE] or
                    checked['source_family'] != self.family or
                    any(checked.get(key, False) for key in OPEN)):
                raise ValueError('Quiet pressure primitive receipt/scope differs')
            _merge_hashes(self.hashes, checked['input_hashes'])
            self.acceptance_loaded = True
        self._hashes = copy.deepcopy(self.hashes)
        self._source_fingerprint = _fingerprint((self.frame_strings, self.endpoint_nodes,
                                                  self.family, self.hashes))

    def assert_source(self):
        if self.hashes != self._hashes:
            raise ValueError('Quiet pressure source hashes changed')
        for name, digest in self.hashes.items():
            if sha(name) != digest:
                raise ValueError('Quiet pressure source changed ' + name)
        if _fingerprint((self.frame_raw, self.quiet_raw)) != self._raw_fingerprint:
            raise ValueError('Accepted quiet/native raw artifact changed')
        if _fingerprint((self.frame_strings, self.endpoint_nodes,
                         self.family, self.hashes)) != self._source_fingerprint:
            raise ValueError('Quiet pressure endpoint source changed')

    def actual_source_identity(self):
        self.assert_source()
        dz_kappa = s.simplify(s.diff(self.kappa, self.Z))
        q_plus_kappa = s.simplify(self.q_common_U0 + self.kappa)
        return dict(mean_common=s.sstr(self.mean_common), kappa=s.sstr(self.kappa),
                    q_common_U0=s.sstr(self.q_common_U0),
                    kappa_Z_derivative=s.sstr(dz_kappa),
                    q_common_plus_kappa=s.sstr(q_plus_kappa),
                    actual_frame_strings=copy.deepcopy(self.frame_strings),
                    same_original_endpoint_function_nodes=copy.deepcopy(self.endpoint_nodes),
                    exact_Z_independent_kappa=(dz_kappa == 0),
                    exact_q_common_at_Rp_is_minus_kappa=(q_plus_kappa == 0))

    def formal_K(self):
        self.assert_source()
        return self.K_definition

    def recover_physical_q(self, K=None, R=None, U=0):
        self.assert_source()
        K = self.formal_K() if K is None else s.sympify(K)
        R = self.R if R is None else s.sympify(R)
        m = K * self.Z / R
        d = 1 - self.Z**2
        L = 1 - self.delta * self.Z**2
        return s.cancel((2 * self.Z * s.sympify(U) -
                         (1 - self.delta) * self.Z * m - d * s.diff(m, self.Z)) / L)

    def primitive(self, K, a, b):
        """Signed m13/m15 pressure primitives on 2 Rc <= R <= Rp."""
        self.assert_source()
        K, a, b = map(s.sympify, (K, a, b))
        if a.is_number and b.is_number and not (bool(a > 0) and bool(b > a)):
            raise ValueError('Positive ordered physical R endpoints required')
        R = self.R
        q = -K / R
        S = -q**2 / 2
        N1p = s.cancel(-S / 2 + 2 * s.diff(q, R) + R * s.diff(q, R, 2))
        i0 = s.integrate(N1p, (R, a, b))
        i1 = s.integrate(R * N1p, (R, a, b))
        return dict(K=K, m=K * self.Z / R, U=s.Integer(0), q=q,
                    q_Z=s.Integer(0), q_y=-q, q_R=s.diff(q, R),
                    q_RR=s.diff(q, R, 2), S=s.cancel(S), N1p=N1p,
                    N1p_integral=s.simplify(i0), R_N1p_integral=s.simplify(i1),
                    m13=s.simplify(i0), m15=s.simplify(-i1),
                    physical_R=True, quiet_segment='actual_2Rc_to_Rp_only',
                    m13_sign='+', m15_sign='-')


def run(owner=None):
    owner = owner if owner is not None else CurrentOriginalN1QuietPressurePrimitives(False)
    proof = owner.actual_source_identity()
    raw = dict(candidate_quiet_pressure_primitives_constructed=True,
               source_family=owner.family,
               actual_analytic_native_inlet_functions=owner.frame_strings,
               same_original_endpoint_function_nodes=owner.endpoint_nodes,
               source_K_definition=dict(K='Rp*Pstar*kappa',
                                        kappa=proof['kappa'],
                                        Rp_source='actual_Rp_native_frame_functions.original_logRp',
                                        Pstar_source='outer.owner.owner(ends).Pstar'),
               actual_source_identity=proof,
               exact_quiet_segment_only=True,
               input_hashes=owner.hashes, **{GATE: False},
               **dict.fromkeys(OPEN, False))
    (HERE / NAME).write_text(json.dumps(raw, separators=(',', ':')) + '\n', encoding='utf8', newline='\n')
    return owner


if __name__ == '__main__':
    run()
