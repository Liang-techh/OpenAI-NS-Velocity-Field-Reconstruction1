"""Accept the explicit raw pressure generator, leaving P0 identity open."""
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_raw_preheat_pressure_operator import (
    CurrentRawPreheatPressureOperator,native_flatten_pressure_integral,
    raw_pressure_operator_source_proof,HERE,NAME,RECEIPT,GATES,OPEN,VIEWS,
    BETA2,BETA0,sha,pack,encode,endpoints,IntervalTaylor)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def independent_flatten_density_fixture():
    """Original cutoff and binomial recurrence, independent of Taylor helper."""
    with mp.workdps(80):
        c=MPIntervalContext();c.dps=80
        exact=SimpleNamespace(flatten=SimpleNamespace(ctx=c,cells=64,prate=c.mpf('1.04')))
        def sigma(x):
            if x<=0:return mp.mpf(0)
            if x>=1:return mp.mpf(1)
            A=mp.exp(-1/(x*x));B=mp.exp(-1/((1-x)*(1-x)))
            return A/(A+B)
        rows=0
        for Z in ('0','.43'):
            z=mp.mpf(Z);q0=1+z*z;q1=2*z;q2=mp.mpf(1)
            def density(v,n):
                sg=sigma(v/100);alpha=2*sg-2
                coeff=[q0**alpha*2**(-2*sg)]
                for k in range(n):
                    value=q1*(alpha-k)*coeff[k]
                    if k:value+=q2*(2*alpha-k+1)*coeff[k-1]
                    coeff.append(value/(q0*(k+1)))
                return mp.exp(-mp.mpf('1.04')*v)*coeff[n]/2
            enclosure=native_flatten_pressure_integral(exact,c.mpf(Z))
            for n in range(6):
                value=mp.quad(lambda v:density(v,n),[0,10,30,50,70,100])
                lo,hi=endpoints(enclosure[n])
                if not lo<=value<=hi:raise ArithmeticError('Native flatten pressure density/axial coefficient fails independent integral')
                rows+=1
            if Z!='0':
                missing_q=mp.quad(lambda v:mp.exp(-mp.mpf('1.04')*v)*2**(-2*sigma(v/100))/2,[0,10,30,50,70,100])
                if abs(missing_q-mp.quad(lambda v:density(v,0),[0,10,30,50,70,100]))<mp.mpf('1e-3'):
                    raise ArithmeticError('Fixture does not detect missing axial q normalization')
        return dict(independent_original_cutoff_density_axial5_coefficients=rows,
            exact_binomial_coefficient_recurrence_used=True,missing_q_normalization_detected=True,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentRawPreheatPressureOperator(require_checked=False)
    proof=raw_pressure_operator_source_proof(field.angular)
    if encode(pack(proof))!=raw['raw_pressure_operator_source_proof'] or not proof['passed']:
        raise ValueError('Explicit raw pressure generator proof changed')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Raw pressure producer claims acceptance')
    if (raw['actual_five_defect_family_sha256'],raw['implicit_source_sha256'],raw['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
        raise ValueError('Raw pressure source/family/datum differs')
    fixture=independent_flatten_density_fixture();count=0
    for name,Z in VIEWS.items():
        packet=field.evaluate(Z)
        if encode(pack(packet))!=raw['current_raw_pressure_operator_views'][name]:raise ValueError('Raw pressure packet changed')
        if set(packet['raw_preheat_pressure_atoms'])!=set(BETA2+BETA0+('z_flatten',)) or any(packet[k] for k in OPEN):
            raise ValueError('Raw pressure fourteen-stage scope differs')
        if not packet['original_P0_not_reset_or_replaced'] or not packet['full_raw_tail_uses_H_one_and_no_angular_bumps']:
            raise ValueError('Raw pressure source changed or P0 reset')
        incoming=field.flat.inlet.incoming(field.ctx.mpf(Z))
        if encode(pack(packet['original_analytic_P0_Taylor_retained']))!=encode(pack(incoming['P0'])):
            raise ValueError('Original analytic P0 was changed')
        total=packet['correlated_actual_native_Rp_to_Rv_prefix']+packet['analytic_flatten_Fflat_Taylor']
        total+=sum((packet['raw_preheat_pressure_atoms'][k] for k in BETA0),IntervalTaylor.constant(field.ctx,0,5))
        if encode(pack(total))!=encode(pack(packet['raw_complete_pressure_integral_Taylor'])):
            raise ValueError('Complete raw pressure operator omits a stage')
        # Source density/FTC proves this equality; interval overlap is only
        # a subsequent diagnostic of units and stage propagation.
        native=packet['raw_forward_terminal_plus_full_preheat_pressure_enclosure']
        for left,right in zip((incoming['P0']+total).coefficients,native.coefficients):
            al,ah=endpoints(left);bl,bh=endpoints(right)
            if max(al,bl)>min(ah,bh):raise ArithmeticError('Raw and native forward pressure bounds contradict source identity')
        count+=6
    hashes=dict(raw['input_hashes']);hashes[NAME]=sha(NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    receipt=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,raw_pressure_operator_source_proof=proof,
        independent_flatten_density_fixture=fixture,raw_pressure_operator_axial5_coefficients=count,
        original_implicit_Fflat_identification_still_required=True,all_passed=True,
        input_hashes=hashes,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(receipt)),indent=2)+'\n',encoding='utf8')
    print('Explicit fourteen-stage native raw pressure generator PASS; original P0 identity and Cp zero remain open',flush=True)
    return receipt


if __name__=='__main__':run()
