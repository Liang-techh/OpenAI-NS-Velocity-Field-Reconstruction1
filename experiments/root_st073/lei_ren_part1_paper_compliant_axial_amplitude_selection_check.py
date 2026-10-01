"""Independent affine/positive-root/C1 identities and full-energy source gates."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_axial_amplitude_selection.json'


def run():
    r=json.loads((HERE/NAME).read_bytes());hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Selected amplitude dependency changed: '+name)
    K,nu1,nu2,u1,u2,v1,v2,a,B=s.symbols('K nu1 nu2 u1 u2 v1 v2 a B',real=True)
    A2=K+nu1*v1*v1+nu2*v2*v2
    A1=2*(nu1*u1*v1+nu2*u2*v2);A0=nu1*u1*u1+nu2*u2*u2-B
    identities={}
    def zero(name,expression):
        if s.simplify(expression)!=0:raise ArithmeticError('Selected amplitude identity failed: '+name)
        identities[name]=True
    zero('actual_affine_energy_quadratic',K*a*a+nu1*(u1+v1*a)**2+nu2*(u2+v2*a)**2-B-(A2*a*a+A1*a+A0))
    C2,C1,C0=s.symbols('C2 C1 C0',real=True);D=s.symbols('D',positive=True)
    root=(-C1+D)/(2*C2)
    zero('positive_root_equation',s.expand(C2*root**2+C1*root+C0).subs(D**2,C1*C1-4*C2*C0))
    zero('positive_root_derivative_denominator',2*C2*root+C1-D)
    h=s.symbols('h',real=True)
    zero('completed_square_equivalence',C2*(-C1/(2*C2)+h)**2+C1*(-C1/(2*C2)+h)+C0-(C2*h*h+C0-C1*C1/(4*C2)))
    u1z,u2z,Bz,az=s.symbols('u1z u2z Bz az',real=True)
    F=K*a*a+nu1*(u1+v1*a)**2+nu2*(u2+v2*a)**2-B
    derivative=s.diff(F,a)*az+s.diff(F,u1)*u1z+s.diff(F,u2)*u2z-Bz
    zero('actual_implicit_C1_recovery',derivative.subs(az,(Bz-2*nu1*(u1+v1*a)*u1z-2*nu2*(u2+v2*a)*u2z)/s.diff(F,a)))
    det,A11,A12,D11,D12,P1,P2=s.symbols('det A11 A12 D11 D12 P1 P2')
    mu=s.symbols('mu',positive=True)
    V1=(-mu*P1*D12+(P2-P1)*A12)/det
    V2=(-(P2-P1)*A11+mu*P1*D11)/det
    zero('affine_row1_slope',(A11*V1+A12*V2+mu*P1).subs(det,A11*D12-A12*D11))
    zero('affine_divided_difference_row_slope',(D11*V1+D12*V2+P2-P1).subs(det,A11*D12-A12*D11))
    # Explicit selected-root gates on the whole continuum, not only fixtures.
    c=MPIntervalContext();c.dps=160
    get=lambda v:read_interval(c,v)
    whole=r['whole_Z_C1_selected_amplitude'];ap=whole['selected_ap_Taylor']['coefficients']
    lo,hi=endpoints(get(ap[0]))
    if not mp.mpf('.9')<lo<=hi<mp.mpf('1.2'):raise ArithmeticError('Whole-Z positive branch failed')
    if not (endpoints(get(whole['low_endpoint_energy_residual']))[1]<0
        and endpoints(get(whole['high_endpoint_energy_residual']))[0]>0
        and endpoints(get(whole['positive_root_derivative_denominator']))[0]>0
        and endpoints(get(whole['positive_discriminant']))[0]>0):
        raise ArithmeticError('Unique positive branch/derivative gates failed')
    for value in r['exact_end_energy_log_margin_checks']:
        if endpoints(get(value))[0]<=0:raise ArithmeticError('Formal end-scale cap failed')
    for value in r['actual_incoming_row_factor_log_margins']:
        if endpoints(get(value))[0]<=0:raise ArithmeticError('Actual incoming row-scale cap failed')
    if not r['actual_incoming_jets_and_first_derivatives_explicitly_scaled']:
        raise ValueError('Incoming derivative source linkage missing')
    if endpoints(get(r['positive_end_energy_cap']))[0]<=0:raise ArithmeticError('End tail zeroed')
    for sample in [whole]+r['samples']:
        coefficients=sample['selected_scaled_end_coefficient_Taylor']
        if len(sample['selected_ap_Taylor']['coefficients'])!=2:
            raise ArithmeticError('Selected amplitude lost C1')
        if not endpoints(get(coefficients[0]['coefficients'][0]))[1]<0<endpoints(get(coefficients[1]['coefficients'][0]))[0]:
            raise ArithmeticError('Selected end coefficient signs failed')
        if not sample['actual_end_linear_and_energy_moment_closure']:
            raise ValueError('Actual selected moment closure missing')
    precise=MPIntervalContext();precise.dps=240
    for sample in [whole]+r['samples']:
        for row,(source,scaled) in enumerate(zip(sample['actual_incoming_moment_Taylor_enclosures'],sample['actual_row_normalized_incoming_Taylor_enclosures'])):
            factor=read_interval(precise,r['actual_incoming_row_factor_positive_caps'][row])
            for raw_value,scaled_value in zip(source['coefficients'],scaled['coefficients']):
                expected=read_interval(precise,raw_value)*factor
                actual=read_interval(precise,scaled_value)
                elo,ehi=endpoints(expected);alo,ahi=endpoints(actual)
                if not alo<=elo<=ehi<=ahi:
                    raise ArithmeticError('Actual incoming value/derivative was not scaled from source jets')
    if not r['full_future_corrected_swirl_energy_used'] or not r['actual_ap_selected']:
        raise ValueError('Energy-selected amplitude prerequisite missing')
    if r['actual_O4_partial_pulse_field_installed'] or r['full_outer_five_moment_match'] or r['whole_outer_cone_certified']:
        raise ValueError('Selected-amplitude scope overclaimed')
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],
        implicit_source_sha256=r['implicit_source_sha256'],independent_selected_amplitude_identities=identities,
        whole_Z_positive_root_and_C1_inverse_checked=True,formal_positive_end_energy_cap_checked=True,
        actual_incoming_value_and_derivative_row_scaling_independently_checked=True,
        selected_c1_negative_c2_positive_checked=True,complete_future_energy_source_bound=True,all_passed=True,
        actual_ap_selected=True,actual_O4_partial_pulse_field_installed=False,
        full_outer_five_moment_match=False,whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Selected actual ap:',len(identities),'independent affine/quadratic/C1 identities and whole-Z positive-root/end-tail/sign gates PASS',flush=True)
    return result


if __name__=='__main__':run()
