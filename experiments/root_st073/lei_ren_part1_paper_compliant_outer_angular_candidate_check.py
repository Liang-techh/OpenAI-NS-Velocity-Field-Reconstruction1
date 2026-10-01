"""Independent angular moment ODEs and raw waiting-source checks."""
# Recomputed for the distinct compliant pressure/moment family.
# Formula origin: lei_ren_part1_paper_shared_outer_angular_candidate_check.py; old .01 receipts remain unchanged.

import hashlib
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_outer_angular_candidate import SharedOuterAngularCandidate
from lei_ren_part1_paper_compliant_outer_initial_check import overlaps
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def symbolic():
    t,mu,k,X0,q=s.symbols('t mu k X0 q',positive=True);sig=s.Function('sigma')(t)
    J,I=s.Function('J')(t),s.Function('I')(t);a=1-mu
    F=(q/2)**sig;X=(X0+I)*s.exp(-a*t)/F
    flatten=s.simplify((s.diff(X,t)-(1-(a+s.diff(sig,t)*s.log(q/2))*X)).subs(s.diff(I,t),s.exp(a*t)*F))
    X=(X0+I)*s.exp(-a*(t-J))
    steep_in=s.simplify((s.diff(X,t)-(1-a*(1-sig)*X)).subs({s.diff(J,t):sig,s.diff(I,t):s.exp(a*(t-J))}))
    X=(X0+I)*s.exp(-k*J)
    steep_out=s.simplify((s.diff(X,t)-(1-k*sig*X)).subs({s.diff(J,t):sig,s.diff(I,t):s.exp(k*J)}))
    power=s.simplify(s.diff(X0+t,t)-1)
    X=1/k+(X0-1/k)*s.exp(-k*t)
    waiting=s.simplify(s.diff(X,t)-(1-k*X))
    # X is derived from the physical angular primitive, including its
    # local logarithmic slope, rather than an independent freely fit field.
    R,U,M,slope=s.symbols('R U M slope',positive=True)
    X=M/(s.sqrt(2)*R**s.Rational(3,2)*U)
    physical=s.simplify(s.diff(X,M)*(R*s.sqrt(2*R)*U)+s.diff(X,R)*R+s.diff(X,U)*U*slope-(1-(s.Rational(3,2)+slope)*X))
    eps,Jc,Xt=s.symbols('eps Jc Xt',positive=True)
    exp_wait=(Xt-1/k)*(1-eps)/(eps*(1/k+Jc))
    tail_X=1/k+(Xt-1/k)/exp_wait
    source=s.simplify(tail_X-(1/k+eps/(1-eps)*(1/k+Jc)))
    residuals=[flatten,steep_in,steep_out,power,waiting,physical,source]
    if any(v!=0 for v in residuals):raise ArithmeticError('Actual angular/waiting identity failed '+str(residuals))
    return dict(angular_stage_ODE_residuals=[str(v) for v in residuals[:5]],
        physical_angular_primitive_residual=str(physical),continuous_preheat_waiting_source_residual=str(source),
        total_symbolic_identities=7)


def run():
    result=symbolic();field=SharedOuterAngularCandidate();c=field.ctx
    name=PREFIX+'compliant_outer_angular_candidate.json';raw=json.loads((HERE/name).read_bytes())
    for source,digest in raw['input_hashes'].items():
        if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Angular candidate changed: '+source)
    interfaces=[]
    with mp.workdps(210):
        for Z in ('-1','0','.5','1'):
            pairs=[('f',field.flatten(Z,1),field.reference_buffer(Z,0)),
                   ('rel',field.reference_buffer(Z),field.steep_in(Z,0)),
                   ('s',field.steep_in(Z),field.steep_power(Z,0)),
                   ('q',field.steep_power(Z),field.before_waiting(Z,0)),
                   ('t',field.before_waiting(Z),field.waiting_field(Z,0))]
            for label,left,right in pairs:
                checks={key:all(overlaps(a,b) for a,b in zip(left[key].coefficients,right[key].coefficients))
                        for key in ('X','relative_log_profile')}
                if not all(checks.values()):raise ArithmeticError('Angular interface mismatch '+label+' '+Z)
                interfaces.append(dict(Z=Z,interface=label,checks=checks))
        Xt=field.before_waiting('0')['X'][0];k=field.restore_rate
        residual=k*field.waiting-c.ln(Xt-1/k)-field.waiting_logone+field.params.log_epsilon+c.ln(1/k+field.collarJ)
        if not endpoints(residual)[0]<=0<=endpoints(residual)[1]:raise ArithmeticError('Waiting log-equation enclosure excludes zero')
        actual_width=endpoints(field.waiting)[1]-endpoints(field.waiting)[0]
        old_width=endpoints(field.params.waiting['root_interval'])[1]-endpoints(field.params.waiting['root_interval'])[0]
        if actual_width>=old_width:raise ArithmeticError('Waiting interval did not improve')
        for api,args in ((field.flatten,([-1,1],'.5')),(field.waiting_field,([-1,1],1))):
            record=api(*args)
            if record['X'].order!=1 or record['relative_log_profile'].order!=1:
                raise ArithmeticError('Whole-axis angular jet lost')
        departure=field.waiting_field('0',1)['X_departure'][0]
        if endpoints(departure)[0]<=0:raise ArithmeticError('Small positive waiting departure discarded')
    result.update(interface_checks=interfaces,first_axial_derivative_enclosures_preserved=True,
        actual_waiting_root_in_original_same_source_enclosure=True,
        actual_waiting_root_strictly_more_tightly_enclosed=True,
        original_continuous_preheat_waiting_log_equation_encloses_zero=True,
        positive_waiting_departure_separately_retained=True,whole_axis_C1_angular_calls_checked=True,
        actual_five_defect_family_sha256=field.initial.family,
        angular_pressure_corrections_installed=False,actual_ap_selected=False,
        exact_heat_exterior_installed=False,temporal_recursion=False,input_hashes=dict(raw['input_hashes']))
    for source in (Path(__file__).name,name):
        result['input_hashes'][source]=hashlib.sha256((HERE/source).read_bytes()).hexdigest()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Angular candidate: 7 independent identities, 20 interfaces and refined same-source waiting PASS; heat/corrections pending',flush=True)
    return result


if __name__=='__main__':run()
