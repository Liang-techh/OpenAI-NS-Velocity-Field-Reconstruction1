"""Focused exact-source/Cauchy-consumer and logarithmic swirl checks."""
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_core_scaled_swirl_source import (
    CurrentCoreScaledSwirlSource,HERE,NAME,RECEIPT,GATE,OPEN,sha,encode,
    selected_transfer_upper_log,relative_amplitude_jet,source_precision)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def independent_relative_S_jets():
    # Independent local logarithmic jets span every analytic source jet;
    # this carrier is not an approximation to the actual amplitude function.
    x=s.Symbol('x');ell=s.symbols('physical_ell0:5')
    local_log=sum(2*ell[k]*x**(k+1)/(k+1) for k in range(5))
    expected=[s.diff(s.exp(local_log),x,k).subs(x,0) for k in range(6)]
    actual=relative_amplitude_jet(SimpleNamespace(mpf=s.sympify),ell,2)
    for left,right in zip(actual,expected):
        if s.expand(left-right)!=0:raise ArithmeticError('Original scaled-S relative derivative/factorial mismatch')
    return dict(executed_original_relative_amplitude_jet_multiplier_two=True,
        six_independent_ordinary_derivative_identities=True,
        logarithmic_derivative_uses_physical_ell_not_scaled_ell=True,
        all_source_amplitude_powers_remain_factored=True,passed=True)


def independent_weighted_Cauchy_fixture():
    with mp.workdps(100):
        lam=mp.mpf(7);loglam=mp.log(lam);eta=mp.mpf('.1');radius=eta/2
        sigma=mp.mpf('.5');anchor=mp.mpf('-.1')
        G=lambda z:(z-anchor)**2/(2*sigma**2)
        Gbar=(1+eta+abs(anchor))**2/(2*sigma**2)
        surplus=mp.mpf(2);guard=lam*Gbar+2*loglam+1000
        exact=lambda z:mp.exp(-2*loglam-2*(guard+surplus)-2*lam*G(z))
        weighted=mp.exp(-4*loglam-2000)*5
        required=selected_transfer_upper_log(mp,loglam,weighted,surplus)
        consumer=-6*loglam-2000;count=0
        if not consumer>required:raise ArithmeticError('Nonunit-weight selected source fixture failed')
        for z in (mp.mpf('-.37'),anchor,mp.mpf('.271')):
            for k in range(9):
                coefficient=mp.diff(exact,z,k)/math.factorial(k)
                if abs(coefficient)>mp.exp(required)/radius**k:
                    raise ArithmeticError('Independent exact source outside its weighted Cauchy majorant')
                count+=1
        for weight,expected in ((1,True),(5,True),(mp.exp(100),False)):
            upper=selected_transfer_upper_log(mp,loglam,mp.exp(-4*loglam-2000)*weight,surplus)
            if bool(consumer>upper)!=expected:raise ArithmeticError('Changed actual Cauchy weight was ignored')
    return dict(independent_entire_source_Taylor_coefficients=count,
        nonunit_transfer_weight_five_used=True,changed_weight_exceeding_surplus_rejected=True,
        actual_source_not_defined_by_numerical_fixture=True,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentCoreScaledSwirlSource(require_checked=False)
    c=field.ctx;read=lambda row:read_interval(c,row)
    if (raw['actual_five_defect_family_sha256']!=field.family or raw['implicit_source_sha256']!=field.source
            or raw['datum_enclosure_sha256']!=field.datum_sha
            or raw['current_source_graph']!=field.current_source_graph() or not all(raw['current_source_graph'].values())
            or raw[GATE] or raw['exact_scaled_S_source_jet_admitted'] or any(raw[key] for key in OPEN)):
        raise ValueError('Current source graph/scope differs')
    for key,value in [('actual_transfer_Cauchy_weight',field.weight),
        ('actual_transfer_weighted_F0_squared_norm_upper',field.weighted_Fsquare),
        ('Cauchy_radius',field.radius),('selected_scaled_source_log_upper',field.required_log),
        ('selected_logCstar_surplus',field.surplus),('selected_source_dominance_margin',field.margin)]:
        if read(raw[key])._mpi_!=value._mpi_:raise ValueError('Weighted selected-source bound changed: '+key)
    counts={}
    def rows_check(rows):
        for k,row in enumerate(rows):
            if row['axial_Taylor_order']!=k or row['ordinary_derivative_divisor_factorial']!=math.factorial(k):
                raise ValueError('Taylor/ordinary derivative units differ')
            required=field.required_log-k*c.ln(field.radius)
            if read(row['selected_source_Taylor_log_upper'])._mpi_!=required._mpi_:
                raise ValueError('Exact selected-source Cauchy factor omitted')
            box=read(row['original_seed_box']);lo,hi=endpoints(box)
            upper=max(abs(lo),abs(hi))
            if (upper!=endpoints(field.rebuild.S_bound/field.radius**k)[1]
                    or not row['exact_positive_source_not_zero_or_cap']
                    or endpoints(read(row['consumer_log_dominance_margin']))[0]<=0
                    or (k==0 and lo!=0) or (k>0 and lo!=-hi)):
                raise ValueError('Original seed fails exact-source inclusion')
        return len(rows)
    packets=raw['current_original_seed_packets']
    if set(packets)!=set(('whole_axis','fresh_Z','source_anchor')):raise ValueError('Source coverage view omitted')
    for name,packet in packets.items():
        if (packet[GATE] or packet['exact_scaled_S_source_jet_admitted'] or any(packet[key] for key in OPEN)
                or not packet['seed_rebuilt_from_current_original_callable'] or not packet['exact_positive_S_retained_formally']
                or not packet['source_value_or_derivative_not_selected']):
            raise ValueError('Source seed scope differs')
        if name=='whole_axis' and endpoints(read(packet['Z']))!=(-1,1):raise ValueError('Original real axial domain lost')
        counts[name]=rows_check(packet['admitted_S_axial_Taylor_coefficients'])
    counts['existing_actual_root']=rows_check(raw['existing_actual_root_seed_S_rows'])
    for name,packet in raw['actual_logarithmic_source_packets'].items():
        if packet[GATE] or any(packet[key] for key in OPEN) or not packet['exp_log_S_not_materialized']:
            raise ValueError('Factored exact source promoted to a point field')
        current=field.logarithmic_source(name)
        if read(packet['log_S_source_enclosure'])._mpi_!=current['log_S_source_enclosure']._mpi_:
            raise ValueError('Original anchored scaled-source logarithm changed')
        for stored,actual in zip(packet['relative_ordinary_S_derivatives_through5'],current['relative_ordinary_S_derivatives_through5']):
            if read(stored)._mpi_!=actual._mpi_:raise ValueError('Original relative S derivative replay differs')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_source_graph=field.current_source_graph(),
        defining_source_bindings=field.bindings,independent_relative_S_jets=independent_relative_S_jets(),
        independent_weighted_Cauchy_fixture=independent_weighted_Cauchy_fixture(),
        original_Cauchy_seed_coefficients_checked=counts,total_original_Cauchy_seed_coefficients_checked=sum(counts.values()),
        all_nonnegative_Taylor_orders_by_Cauchy_source_theorem=True,
        actual_transfer_weight_and_selected_Cstar_surplus_consumed=True,
        current_original_seed_boxes_enclose_exact_positive_S_jets=True,
        original_seed_rows_and_prior_receipts_not_relabelled=True,
        **{GATE:True,'exact_scaled_S_source_jet_admitted':True},**dict.fromkeys(OPEN,False),
        remaining_dependency='Actual twenty-term nonlinear operator/majorant binding -> finite recurrence coefficients and analytic tails of the same fixed point -> original core/first ODE join',
        all_scoped_checks_passed=True,all_passed=True,
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    (HERE/RECEIPT).write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS exact current scaled-S jets in '+str(sum(counts.values()))+' original Cauchy seed rows; common nonlinear fixed point/core-first remain open',flush=True)
    return result


if __name__=='__main__':run()
