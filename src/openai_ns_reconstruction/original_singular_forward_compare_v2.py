"""Relate the singular canonical slices to the saved 212-row forward ledger.

Point diagnostics compare serialized finite coefficients in fixed exact source
units. The complete raw rows and their log-scale errors remain authoritative.
"""
from __future__ import annotations

import ast
import hashlib
import json
from decimal import Decimal

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext


INNER = 'src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py'
MICRO = 'src/openai_ns_reconstruction/original_first_micro_m15.py'
FULL = 'src/openai_ns_reconstruction/original_first_micro_full_m15.py'
COLLAR = 'src/openai_ns_reconstruction/original_first_collar_transport.py'
JETS = 'src/openai_ns_reconstruction/original_first_collar_jets.py'
Z_KEYS = ('-.3', '0', '.3')


def require(value, message):
    if not value:
        raise ValueError(message)


def read(iv, row):
    return iv.make_mpf(tuple(tuple(row[f'{side}_exact_mpf_tuple']) for side in ('lower', 'upper')))


def old_box(iv, row):
    return iv.make_mpf(tuple(tuple(bits) for bits in row['source_enclosure_exact_mpf_tuples']))


def digest(data):
    return hashlib.sha256(json.dumps(data, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def sign(value):
    if value.a > 0:
        return 'positive'
    if value.b < 0:
        return 'negative'
    return 'exact_zero' if value.a == value.b == 0 else 'unresolved'


def record(mp, value):
    return dict(lower_exact_mpf_tuple=list(value._mpi_[0]), upper_exact_mpf_tuple=list(value._mpi_[1]),
                display=[mp.nstr(mp.make_mpf(bits), 16) for bits in value._mpi_], sign=sign(value))


def contains(container, value):
    return bool(container.a <= value.a and value.b <= container.b)


def _negative(bits):
    return [1-bits[0], *bits[1:]] if bits[1] else bits


def _method(text, owner, name):
    cls = next(node for node in ast.parse(text).body if isinstance(node, ast.ClassDef) and node.name == owner)
    return next(node for node in cls.body if isinstance(node, ast.FunctionDef) and node.name == name)


def _assignment(method, target, expression):
    nodes = [node.value for node in ast.walk(method) if isinstance(node, ast.Assign)
             and any(isinstance(key, ast.Name) and key.id == target for key in node.targets)]
    wanted = ast.dump(ast.parse(expression, mode='eval').body)
    require(len(nodes) == 1 and ast.dump(nodes[0]) == wanted, f'changed archived source recipe {target}')


def _source_recipes(inner, micro, quiet, codes):
    for receipt, sources in ((inner, (INNER, COLLAR, JETS)), (micro, tuple(codes)), (quiet, tuple(codes))):
        require(receipt['checks']['passed'] and all(receipt['implementation_hashes'][name]
                == hashlib.sha256(codes[name].encode()).hexdigest() for name in sources), 'changed source implementation')
    a = _method(codes[INNER], 'OriginalN1InnerM15Numeric', '_prepare')
    b = _method(codes[MICRO], 'OriginalFirstMicroM15', 'prepare')
    lower_sigma = _method(codes[FULL], 'OriginalFirstMicroFullM15', 'sigma_atom')
    require('return super().sigma_atom(s)' in ast.unparse(lower_sigma),
            'full first-micro lower cutoff no longer uses the inherited log-odds endpoint')
    transport = _method(codes[COLLAR], 'OriginalFirstCollarTransport', '_inputs')
    jets = _method(codes[JETS], 'OriginalFirstCollarJets', 'inlet_jets')
    _assignment(a, 'singular', '2*rr*rr*zi*one*(phi/comp)*psir/L')
    _assignment(a, 'basep', "inp['radial']['Psi'][1].truncate(2)*z.truncate(2)*32/(1-z.truncate(2)**2*de)")
    _assignment(b, 'singular', '2*r*r*z*one*(phi/comp)*psir/L')
    _assignment(b, 'ap', "inp['radial']['Psi'][1].truncate(2)*zp*32/(1-zp*zp*self.o.core.delta)")
    _assignment(transport, 'inlet', 'o.inlet_jets(Z,degree=24,order=order)')
    text = ast.unparse(jets)
    require('J, r = (self.Jet, c.mpf(4))' in text and 'for n in range(degree):' in text
            and 'radial[name][i] +=' in text, 'finite radius-4 inlet or tail changed')


def _endpoint_rows(z, inner, micro, quiet):
    final = [p for p in quiet['samples'] if p['Z'] == z and p['chart'] == 'quiet' and p['native_coordinate'] == '1']
    require(len(final) == 1 and len(final[0]['m15_C0_C1_terms']) == 212, 'missing unique forward Rp ledger')
    terms = final[0]['m15_C0_C1_terms']
    parents = [p for p in inner['samples'] if p['Z'] == z and p['q'] == '2' and p['m15_C0_C1_terms'] == terms[:6]]
    following = [p for p in micro['samples'] if p['Z'] == z and p['incoming_q2_m15_terms'] == terms[:6]
                 and p['coordinate'] == {'kind':'exact_decimal_phase','value':'1'}
                 and p['m15_C0_C1_terms'] == terms[:13]]
    require(parents and following and all(p['original_sc_tuple'] == parents[0]['exact_sc_tuple'] for p in following),
            'terminal ledger has no retained same-source q2/micro prefix')
    a, e, b = (terms[k] for k in (4, 5, 8))
    require(a['name'] == 'first_lower_singular_pressure' and e['name'] == 'first_lower_singular_width_error'
            and b['name'] == 'first_micro_singular_lower'
            and a['actual_function_binding']['A0'] == '32Z*Psi_core_rho(4)/L'
            and e['actual_function_binding']['proof'] == '|A-A0|_norm<=hb*E, integrate against positive d_sigma'
            and b['actual_function_binding']['endpoint_sign'] == -1, 'incorrect endpoint or error sector')
    require([x['source_enclosure_exact_mpf_tuples'] for x in a['log_scale_terms']]
            == [x['source_enclosure_exact_mpf_tuples'] for x in b['log_scale_terms']]
            and len(a['log_scale_terms']) == 2 and len(e['log_scale_terms']) == 3
            and a['log_scale_terms'][0] == e['log_scale_terms'][0]
            and a['log_scale_terms'][1] == e['log_scale_terms'][2], 'source scale or d-sigma endpoint changed')
    for ca, cb in zip(a['coefficient_C0_C1'], b['coefficient_C0_C1'], strict=True):
        require(ca['source_enclosure_exact_mpf_tuples']
                == [_negative(v) for v in cb['source_enclosure_exact_mpf_tuples'][::-1]]
                and Decimal(ca['approximation']) == Decimal(cb['approximation']).copy_negate(),
                'archived endpoint opposites changed')
    return terms, parents[0], a, e


def _inv_p2(iv, log_p2):
    # Peel off the huge binary exponent exactly, and exponentiate only the
    # remaining bounded log factor. The physical field scale stays formal.
    bits = (-log_p2/iv.ln(2))._mpi_[0]
    s, man, exponent, _ = bits
    require(s == 1 and exponent < 0, 'unexpected Pstar log scale')
    shift = (-man) >> (-exponent)
    residual = -log_p2-shift*iv.ln(2)
    require(residual.a > -1 and residual.b < 2, 'unresolved binary log split')
    scale = iv.make_mpf(((0, 1, shift, 1),)*2)
    return scale*iv.exp(residual)


def _canonical_principal(iv, sample, Z, delta):
    if Z.a == Z.b == 0:
        H = read(iv, sample['ordinary_C0_C1']['inlet_hydro_and_nonlinear'][0])
        return iv.mpf(0), 32*H
    W, W_Z = (read(iv, row) for row in sample['ordinary_C0_C1']['inlet'])
    L = 1-delta*Z**2
    return 32*Z*W/L, 32*((W+Z*W_Z)/L+2*delta*Z**2*W/L**2)


def evaluate(*, whole, q2, inner, micro, quiet, native, inner_code, micro_code, collar_code, jets_code, full_code):
    codes = {INNER: inner_code, MICRO: micro_code, COLLAR: collar_code, JETS: jets_code, FULL: full_code}
    _source_recipes(inner, micro, quiet, codes)
    family = q2['source_family']
    leading = native['actual_native_collar_transport']
    require(all(r['source_family'] == family for r in (inner, micro, quiet, native))
            and leading['passed'] and leading['actual_quotient_equals_ratio']
            and leading['native_defining_leading_function_equals_canonical_collar']
            and leading['domain']['canonical_q_interval'] == [0, 2]
            and whole['passed'] and whole['covers_every_real_Z'], 'unbound canonical collar family')
    packet = q2['actual_correlated_q2_inlet']
    iv, mp = MPIntervalContext(), MPContext()
    iv.prec = 1800
    S = read(iv, packet['selected_phase'])
    require(S.a == S.b and S.a > 0, 'expected exact positive sc/2')
    odds = 1/(1-S)**2-1/S**2
    require(odds.b < -1024, 'small logistic odds cap unavailable')
    kappa = 1/(1+iv.mpf([0, iv.mpf(2)**(-1024)]))
    delta = read(iv, q2['actual_selected_leading_pressure_primitive']['delta_enclosure'])
    p_inv = _inv_p2(iv, read(iv, packet['source_log_bases'][2]))
    rows = {}
    for key in Z_KEYS:
        ledger, parent, principal, error = _endpoint_rows(key, inner, micro, quiet)
        step = tuple(parent['exact_phase_tuple'])
        require(tuple(parent['exact_sc_tuple']) == (step[0],step[1],step[2]+2,step[3])
                and S._mpi_[0] == (step[0],step[1],step[2]+1,step[3])
                and parent['exact_sc_tuple'] == leading['domain']['sc_exact_mpf_tuple'], 'incorrect q2/sc coordinate')
        old_log = principal['log_scale_terms'][1]
        require(contains(old_box(iv, old_log), odds), 'exact common logistic log is outside archived enclosure')
        sample = whole['exact_samples']['0' if key == '0' else ('-0.3' if key[0] == '-' else '0.3')]
        factor = iv.mpf(1) if key == '0' else p_inv
        A = _canonical_principal(iv, sample, iv.mpf(key), delta)
        result = []
        for k, exact in enumerate(A):
            approx = principal['coefficient_C0_C1'][k]
            cap = error['coefficient_C0_C1'][k]
            require(cap['approximation'] == '0.0', 'term6 has an unexpected signed approximation')
            canonical5 = kappa*exact
            canonical6 = kappa*read(iv, sample['ordinary_C0_C1']['primitive'][k])
            archived5 = old_box(iv, approx)*factor
            archived6 = old_box(iv, cap)*p_inv
            result.append(dict(ordinary_Z_order=k, canonical_principal=record(mp, canonical5),
                archived_principal_coefficient=record(mp, archived5),
                principal_diagnostic_difference=record(mp, canonical5-iv.mpf(approx['approximation'])*factor),
                principal_full_containment_diagnostic=contains(archived5, canonical5),
                signed_remainder_vs_zero_approximation=record(mp, canonical6),
                archived_remainder_coefficient_cap=record(mp, archived6),
                signed_remainder_inside_archived_cap=contains(archived6, canonical6)))
        rows[key] = dict(ledger_count=len(ledger), all_forward_terms_digest=digest(ledger),
            raw_singular_rows={str(i): ledger[i] for i in (4, 5, 8, 9)},
            reference_pair='canonical +A0 at slot 4 and -A0 at slot 8, same s=sc/2',
            shared_principal_units='epsilon^2*exp(odds)' + ('' if key == '0' else '*Pstar^2'),
            remainder_units='epsilon^2*hb*exp(odds)*Pstar^2',
            exact_log_minus_serialized_log=record(mp, odds-iv.mpf(old_log['approximation'])), ordinary_Z=result)
    return dict(passed=True, source_principal_recipe_linked=True, same_endpoint_reference=[4,8],
                source_reference_contract='epsilon^2*kappa*exp(odds)*(A0-A0)=0 on the fixed interface before boxing',
                exact_source_scale_atoms_materialized=False, numerical_diagnostics_use_exact_source_scale=True,
                saved_approximation_equals_exact_function=False, signed_remainder_source_view=True,
                distinct_width_caps_preserved=[5,9], samples=rows, forward_numeric_ledger_mutated=False,
                forward212_full_source_equivalence=False, terminal_comparison_closed=False,
                interval_overlap_used_as_proof=False, ancestors_executed=False)
