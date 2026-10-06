"""Current pressure quadratic/full-future balance, without zeroing raw Cp."""
import json
from pathlib import Path
import mpmath as mp

from lei_ren_part1_paper_compliant_current_pressure_terminal_balance import (
    CurrentPressureTerminalBalance,current_pressure_balance_proof,HERE,PREFIX,
    NAME,RECEIPT,GATES,SCOPES,VIEWS,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def independent_signed_bump_pressure_fixture():
    """Measure signed, quadratic swirl changes on both actual bump supports."""
    with mp.workdps(85):
        mu=mp.mpf('.02');p=1+2*mu;ell=mp.mpf('.15')
        raw=lambda x:mp.exp(-1/(1-x*x)) if abs(x)<1 else mp.mpf(0)
        norm=mp.quad(raw,[-1,0,1]);beta=lambda v:raw(v/ell)/(ell*norm)
        B=mp.quad(lambda v:mp.exp(-p*v)*beta(v),[-ell,0,ell])
        D=mp.quad(lambda v:mp.exp(-p*v)*beta(v)**2,[-ell,0,ell])
        breaks=list(map(mp.mpf,('-4','-3.15','-3','-2.85','-1.15','-1','-.85','0')))
        rows=0
        for d1,d2 in (('.0013','-.0022'),('.0018','-.003'),('.0007','-.0012')):
            d1=mp.mpf(d1);d2=mp.mpf(d2)
            h=lambda v:d1*beta(v+3)+d2*beta(v+1)
            direct=mp.quad(lambda v:mp.exp(-p*v)*((1+h(v))**2-1)/2,breaks)
            weight=(d1*B+d1*d1*D/2)*mp.exp(3*p)+(d2*B+d2*d2*D/2)*mp.exp(p)
            if abs(direct-weight)>mp.mpf('1e-65'):raise ArithmeticError('Direct original signed pressure density disagrees with both bump weights')
            # A missing square term must fail this independent density check.
            linear=B*(d1*mp.exp(3*p)+d2*mp.exp(p))
            if abs(direct-linear)<mp.mpf('1e-8'):raise ArithmeticError('Fixture does not detect missing bump-square pressure')
            rows+=1
        return dict(direct_supported_swirl_pressure_integral_rows=rows,
            both_signed_bumps_and_quadratic_terms_measured=True,
            dropping_bump_square_term_fails=True,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentPressureTerminalBalance(require_checked=False)
    proof=current_pressure_balance_proof(field.angular)
    if encode(pack(proof))!=raw['current_pressure_balance_source_proof'] or not all(proof['identities'].values()):raise ValueError('Current pressure function proof changed')
    if any(raw[k] for k in GATES+SCOPES):raise ValueError('Pressure producer claims acceptance')
    if (raw['actual_five_defect_family_sha256'],raw['implicit_source_sha256'],raw['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):raise ValueError('Current pressure source/family/datum differs')
    fixture=independent_signed_bump_pressure_fixture();count=0
    for name,Z in VIEWS.items():
        packet=field.evaluate(Z)
        if encode(pack(packet))!=raw['current_pressure_balance_views'][name]:raise ValueError('Current pressure packet changed')
        if not packet['raw_pressure_constant_retained_not_zeroed'] or any(packet[k] for k in SCOPES):raise ValueError('Current pressure balance overclaims raw datum or global scope')
        if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in packet['source_proved_bump_minus_heat_loss_Taylor'].coefficients):raise ValueError('Pressure balance function jet not zero through5')
        actual=field.angular.terminal_constants(Z)['pressure_infinity']
        if encode(pack(actual))!=encode(pack(packet['actual_pressure_infinity_retained'])):raise ValueError('Actual Cp was reset')
        # Overlap is a diagnostic only; the source equations prove balance.
        for a,b in zip(actual.coefficients,packet['same_raw_preheat_pressure_constant_enclosure'].coefficients):
            al,ah=endpoints(a);bl,bh=endpoints(b)
            if max(al,bl)>min(ah,bh):raise ArithmeticError('Forward and raw pressure enclosures contradict the source identity')
        count+=6
    hashes=dict(raw['input_hashes']);hashes[NAME]=sha(NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_pressure_balance_source_proof=proof,
        independent_signed_pressure_fixture=fixture,source_balance_axial5_coefficients=count,
        raw_preheat_constant_kept_pending_function_identification=True,
        all_passed=True,input_hashes=hashes,**dict.fromkeys(GATES,True),**dict.fromkeys(SCOPES,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current quadratic/full-Gamma pressure balance PASS; original raw pressure datum bridge remains open',flush=True)
    return result


if __name__=='__main__':run()
