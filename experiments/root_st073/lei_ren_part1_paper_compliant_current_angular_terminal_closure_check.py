"""Current angular terminal cancellation and retained-pressure scope checks."""
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_angular_terminal_closure import (
    CurrentAngularTerminalClosure,HERE,PREFIX,NAME,RECEIPT,GATES,SCOPES,VIEWS,
    sha,pack,encode,endpoints,current_angular_source_proof)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def independent_corrected_terminal_fixture():
    """Finite parameters distinguish raw waiting reference from corrected XT."""
    with mp.workdps(95):
        mu=mp.mpf('.02');r=1-mu;k=mp.mpf('.9');eps=mp.mpf('.01')
        L=mp.mpf(4);Ts=mp.mpf('.4');S=mp.mpf('.001')
        Iin=mp.mpf('1.3');Iout=mp.mpf('.7');xf0=mp.mpf('2.1')
        def sigma(v):
            if v<=0:return mp.mpf(0)
            if v>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/v**2-1/(1-v)**2))
        def phi(v):return mp.exp(-4/(3-v)**2) if v<3 else mp.mpf(0)
        J=mp.quad(lambda v:mp.exp(k*v)*(1-sigma(v)+sigma(v)*phi(v)),[0,.2,.5,1,2,2.8,3])
        raw=lambda v:mp.exp(-1/(1-v*v)) if abs(v)<1 else mp.mpf(0)
        norm=mp.quad(raw,[-1,0,1]);ell=mp.mpf('.15')
        Aweight=mp.quad(lambda v:mp.exp(r*ell*v)*raw(v)/norm,[-1,0,1])
        def transfer(xf,A):
            XR=1/r+(xf-1/r)*mp.exp(-r*L)+A
            XS=(XR+Iin)*mp.exp(-r/2);XQ=XS+Ts;XT=(XQ+Iout)*mp.exp(-k/2)
            return XT
        XT0=transfer(xf0,0)
        W=(mp.log(XT0-1/k)+mp.log(1-eps)-mp.log(eps)-mp.log(1/k+J))/k
        if W<=0:raise ArithmeticError('Moderate fixture raw waiting root is not positive')
        logmult=(r+k)/2+k*W-mp.log(1-eps);checks=0;axis_increment=None
        for value in ('-1','-.47','0','.47','1'):
            z=mp.mpf(value);xf=xf0-mp.mpf('.3')*z*z;Theta=mp.mpf('.2')*(1-z*z)
            rhs=mp.exp(-r*L)*(xf0-xf)+mp.exp(logmult)*S*Theta
            # Two distinct nonzero supported bumps realize the exact equation.
            d1=rhs*mp.exp(3*r)/(3*Aweight);d2=2*rhs*mp.exp(r)/(3*Aweight)
            pastA=Aweight*(d1*mp.exp(-3*r)+d2*mp.exp(-r))
            XT=transfer(xf,pastA);tail=1/k+(XT-1/k)*mp.exp(-k*W)
            target=1/k+eps*J+S*Theta
            if abs((1-eps)*tail-target)>mp.mpf('1e-80'):
                raise ArithmeticError('Independent corrected terminal misses full heat angular future')
            checks+=1
            if z==0:
                rawtail=1/k+(transfer(xf,0)-1/k)*mp.exp(-k*W)
                axis_increment=(1-eps)*(tail-rawtail)
                if abs(axis_increment-S*Theta)>mp.mpf('1e-80') or axis_increment<=0:
                    raise ArithmeticError('Nonzero Gamma repair at the axis was discarded')
                if abs((1-eps)*rawtail-target)<mp.mpf('1e-6'):
                    raise ArithmeticError('Fixture fails to distinguish raw and corrected terminal')
        return dict(independent_corrected_terminal_rows=checks,
            original_sigma_phi_J_integral_evaluated=True,normalized_original_beta_A_integral_evaluated=True,
            prescribed_waiting_computed_from_raw_XT0=True,axis_nonzero_heat_correction=axis_increment,
            omitting_actual_bump_at_axis_fails_terminal_identity=True,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentAngularTerminalClosure(require_checked=False)
    if (raw['actual_five_defect_family_sha256'],raw['implicit_source_sha256'],raw['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
        raise ValueError('Current angular closure source/family/datum differs')
    if any(raw[k] for k in GATES+SCOPES):raise ValueError('Angular producer claims acceptance')
    proof=current_angular_source_proof(field.exact)
    if encode(pack(proof))!=raw['current_actual_angular_source_identity_proof'] or not all(proof['identities'].values()):
        raise ValueError('Current actual angular source proof differs')
    fixture=independent_corrected_terminal_fixture();theta_count=angular_count=0
    for name,(chart,Z,t) in VIEWS.items():
        packet=field.evaluate(chart,Z,t)
        if encode(pack(packet))!=raw['current_angular_terminal_views'][name]:raise ValueError('Current angular terminal packet differs')
        if any(packet[k] for k in SCOPES) or not packet['actual_pressure_constant_not_eliminated']:
            raise ValueError('Current angular closure overstates pressure/physical/global scope')
        if any(endpoints(value)!=(mp.mpf(0),mp.mpf(0)) for value in packet['actual_Dtheta_Taylor'].coefficients):
            raise ValueError('Actual angular constant is not source-proved zero through5')
        constants=field.terminal_constants(Z);terminal=constants['current_repaired_forward_terminal']
        # Consistency with the rebound forward enclosures is checked AFTER
        # the exact source identity. It never establishes the zero constant.
        if any(not endpoints(v)[0]<=0<=endpoints(v)[1] for v in packet['original_forward_subtraction_enclosure'].coefficients):
            raise ArithmeticError('Rebound forward enclosure contradicts the proved current angular identity')
        expected_Cp=terminal['Ptail']+packet['current_full_collar_future_pressure']*terminal['pressure_scale']
        if encode(pack(packet['actual_pressure_infinity_retained']))!=encode(pack(expected_Cp)):
            raise ValueError('Current pressure constant does not use rebound forward and full future sources')
        if not packet['current_exact_repair_runtime_used'] or not terminal['current_exact_repair_directly_consumed']:
            raise ValueError('Current runtime fails to consume the new exact repair')
        if field.heat.exact_logS_terms is not field.exact.repair.heat.logS_terms or any(hasattr(field.heat,v) for v in ('outer','steep','data','collar','exterior')):
            raise ValueError('Restricted shape view exposes a stale terminal source')
        rows=packet['actual_stress_factored_mixed4']['theta_Qtheta']
        if len(rows)!=15 or packet['actual_angular_Taylor_after_source_closure'].order!=5:
            raise ValueError('Current angular mixed4 or axial5 rows missing')
        if chart=='heat_exterior' and any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in rows.values()):
            raise ValueError('Full Gamma angular stress is not zero after actual closure')
        theta_count+=len(rows);angular_count+=6
    hashes=dict(raw['input_hashes']);hashes[NAME]=sha(NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_actual_angular_source_identity_proof=proof,
        independent_corrected_terminal_fixture=fixture,current_angular_stress_mixed4_rows=theta_count,
        current_stable_angular_axial5_coefficients=angular_count,current_actual_Dtheta_zero_axial_order=5,
        pressure_infinity_and_axial_pressure_stress_retained=True,
        current_runtime_directly_consumes_exact_repair_coefficients_waiting_radius=True,
        pressure_constant_recomputed_from_same_rebound_forward_and_full_future=True,
        full_Gamma_angular_stress_zero_but_full_stress_zero_not_claimed=True,
        all_passed=True,input_hashes=hashes,**dict.fromkeys(GATES,True),**dict.fromkeys(SCOPES,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual angular terminal source identity PASS: '+str(theta_count)+' stress4 rows; Cp/pressure closure pending',flush=True)
    return result


if __name__=='__main__':run()
