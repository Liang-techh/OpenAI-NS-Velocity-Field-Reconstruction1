"""Independent physical-row, scaling and saddle checks for F33."""
# Recomputed for the distinct compliant pressure/moment family.
# Formula origin: lei_ren_part1_paper_shared_outer_pulse_map_check.py; old .01 receipts remain unchanged.

import hashlib
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_outer_pulse_map import SharedOuterPulseMap
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def symbolic():
    t,mu,Rp,E0=s.symbols('t mu Rp E0',positive=True);B=s.symbols('B')
    R=Rp*s.exp(t);E=E0*s.exp(-(s.Rational(1,2)+mu)*t);V=E*B
    rows=[s.simplify(R*V/(Rp*E0)-s.exp((s.Rational(1,2)-mu)*t)*B),
          s.simplify(R*s.sqrt(2*R)*E*V/(s.sqrt(2)*Rp**s.Rational(3,2)*E0**2)-s.exp((s.Rational(1,2)-2*mu)*t)*B),
          s.simplify(R*(V**2-E**2/2)/(Rp*E0**2)-s.exp(-2*mu*t)*(B**2-s.Rational(1,2)))]
    u0,x=s.symbols('u0 x',positive=True);v=1+x*u0/s.sqrt(6);L=1/u0**2;u=u0*v;k0=2/u0**3
    phi=x**2*(2*v+1)/(6*v**2)
    saddle=[s.simplify(k0*u+1/u**2-3*L-phi),s.simplify(s.diff(phi,x,2)-1/v**4),
            s.simplify(s.diff(phi,x)-x*(v*v+v+1)/(3*v**3)),
            s.simplify(s.diff(u,x)-u0**2/s.sqrt(6))]
    A11,A12,D1,D2,R1,R2=s.symbols('A11 A12 D1 D2 R1 R2')
    determinant=A11*D2-A12*D1
    C1=(-mu*R1*D2+(R2-R1)*A12)/determinant
    C2=(-(R2-R1)*A11+mu*R1*D1)/determinant
    # Check the ORIGINAL two rows after restoring common exp(logpref).
    correction=[s.simplify(A11*C1+A12*C2+mu*R1),
                s.simplify((A11+mu*D1)*C1+(A12+mu*D2)*C2+mu*R2)]
    x=s.symbols('x',positive=True)
    quotient=s.simplify((s.exp((-s.Rational(1,2)+2*mu)*x)-s.exp((-s.Rational(1,2)+mu)*x))/mu
                       -s.exp((-s.Rational(1,2)+mu)*x)*(s.exp(mu*x)-1)/mu)
    b1,b2=s.symbols('b1 b2',positive=True);l1=s.Rational(1,2)-mu;l2=s.Rational(1,2)-2*mu
    det=s.exp(-3*l1-l2)*b1*b2-s.exp(-l1-3*l2)*b1*b2
    detcheck=s.simplify(s.expand((det+2*s.exp(-2+6*mu)*s.sinh(mu)*b1*b2).rewrite(s.exp)))
    # Disjoint supports remove all cross terms. Re-derive Section7.34
    # from the physical energy, including the FUTURE corrected swirl.
    Kp,K1,K2,a,c1,c2,incoming,tail=s.symbols('Kp K1 K2 a c1 c2 incoming tail')
    accumulated=incoming+a*a*Kp/mu+K1*c1*c1+K2*c2*c2-(1-s.exp(-26))/(4*mu)-tail/2
    energy=s.simplify(mu*accumulated-(Kp*a*a+mu*(K1*c1*c1+K2*c2*c2)-(1-s.exp(-26))/4+mu*incoming-mu*tail/2))
    residuals=rows+saddle+correction+[quotient,detcheck,energy]
    if any(v!=0 for v in residuals):raise ArithmeticError('Independent pulse identity failed: '+str(residuals))
    return dict(physical_linear_and_energy_weight_residuals=[str(v) for v in rows],
        saddle_and_convexity_residuals=[str(v) for v in saddle],
        original_two_linear_row_residuals=[str(v) for v in correction],
        exact_divided_difference_residual=str(quotient),determinant_residual=str(detcheck),
        full_future_tail_energy_normalization_residual=str(energy),total_symbolic_identities=len(residuals))


def run():
    result=symbolic();field=SharedOuterPulseMap();c=field.ctx
    name=PREFIX+'compliant_outer_pulse_map.json';raw=json.loads((HERE/name).read_bytes())
    for source,digest in raw['input_hashes'].items():
        if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Pulse map changed: '+source)
    if not raw['near_equal_rows_recovered_by_analytic_divided_difference'] or raw['actual_ap_selected']:
        raise ValueError('Pulse map scope or divided-difference provenance invalid')
    with mp.workdps(210):
        det=field.basis['divided_difference_determinant']
        if endpoints(det)[1]>=0:raise ArithmeticError('Divided-difference inverse lacks determinant gap')
        if not endpoints(field.Kpulse)[0]>mp.mpf('.24') or not endpoints(field.Kpulse)[1]<mp.mpf('.246'):
            raise ArithmeticError('Directed Kpulse range mismatch')
        diagnostics=[]
        for Z,a in [('-1','.9'),('0','1'),('.5','1.2'),('1','1'),([-1,1],'1')]:
            trial=field.coefficients(Z,a);C=trial['scaled_coefficients'];cap=field.incoming_cap
            rho=[c.mpf([-endpoints(cap)[1],endpoints(cap)[1]])+field.rows['scaled_full_rows'][i]*c.mpf(a) for i in (0,1)]
            A=field.basis['first_row'];D=field.basis['exact_divided_difference_row']
            r1=C[0]*A[0]+C[1]*A[1]+rho[0]*field.mu
            r2=C[0]*(A[0]+field.mu*D[0])+C[1]*(A[1]+field.mu*D[1])+rho[1]*field.mu
            for residual in (r1,r2):
                if not endpoints(residual[0])[0]<=0<=endpoints(residual[0])[1]:
                    raise ArithmeticError('Original row diagnostic does not enclose zero')
            if not endpoints(C[0][0])[1]<0<endpoints(C[1][0])[0]:raise ArithmeticError('Coefficient sign gate failed')
            diagnostics.append(dict(Z=Z,trial_a=a,original_row_value_enclosures_contain_zero=True,
                c1_negative_c2_positive=True,derivatives_at_fixed_trial_a_enclosed=True))
        for row in (0,1):
            if endpoints(field.rows['positive_scaled_tail_bounds'][row])[0]<=0:
                raise ArithmeticError('True positive pulse tail discarded')
        if endpoints(field.logscale)[1]>=-1000:raise ArithmeticError('Formal tiny coefficient scale invalid')
    result.update(actual_trial_source_and_interval_inverse_checked=True,diagnostics=diagnostics,
        directed_original_row_residuals_are_diagnostics_not_closure_proof=True,
        positive_unresolved_pulse_and_incoming_tails_preserved=True,
        actual_five_defect_family_sha256=field.initial.family,
        actual_ap_selected=False,corrected_swirl_energy_tail_available=False,
        complete_O4_corrected_field_built=False,heat_exterior_matched=False,temporal_recursion=False,
        input_hashes=dict(raw['input_hashes']))
    for source in (Path(__file__).name,name):
        result['input_hashes'][source]=hashlib.sha256((HERE/source).read_bytes()).hexdigest()
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Independent actual trial pulse: 12 symbolic identities, five source/inverse examples PASS; ap remains pending heat tail',flush=True)
    return result


if __name__=='__main__':run()
