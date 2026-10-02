"""Independent fifth-order branch, Gamma, flatten and selected-root checks."""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_fifth_axial_jets import angular_fifth, selected_fifth, preheat_fifth
from lei_ren_part1_paper_compliant_angular_high_jets import gamma_deficit_jets, integrated_gamma_tails, log_taylor
from lei_ren_part1_paper_compliant_outer_initial import stable_sigma
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_fifth_axial_jets.json'


def contains(label, interval, value, error=0):
    lo,hi=endpoints(interval)
    if not lo <= value-error <= value+error <= hi:
        raise ArithmeticError('Independent fifth coefficient outside enclosure: '+label)


def identities():
    out={}
    def zero(name,expr):
        if s.simplify(expr)!=0: raise ArithmeticError('Fifth-order identity failed: '+name)
        out[name]=True
    h=s.symbols('h'); x=s.symbols('x0:6'); y=s.symbols('y0:6'); p,q,K=s.symbols('p q K')
    b1,b2=s.symbols('b1 b2'); J21=1+2*K*x[0]; J22=q*(1+2*K*y[0]); det=p*J22-J21
    known=K*sum(x[i]*x[5-i]+q*y[i]*y[5-i] for i in range(1,5))
    xn=(J22*b1-(b2-known))/det; yn=(p*(b2-known)-J21*b1)/det
    fx=sum(x[i]*h**i for i in range(6)); fy=sum(y[i]*h**i for i in range(6))
    zero('angular_fifth_linear_equation',(s.expand(p*fx+fy).coeff(h,5)-b1).subs({x[5]:xn,y[5]:yn}))
    zero('angular_fifth_quadratic_equation',(s.expand(fx+q*fy+K*(fx**2+q*fy**2)).coeff(h,5)-b2).subs({x[5]:xn,y[5]:yn}))
    A2=s.symbols('A2'); A1=s.symbols('A1_0:6'); A0=s.symbols('A0_0:6'); a=s.symbols('a0:6')
    D=2*A2*a[0]+A1[0]
    a5=-(A2*sum(a[i]*a[5-i] for i in range(1,5))+sum(A1[i]*a[5-i] for i in range(1,6))+A0[5])/D
    fa=sum(a[i]*h**i for i in range(6)); f1=sum(A1[i]*h**i for i in range(6)); f0=sum(A0[i]*h**i for i in range(6))
    zero('selected_positive_root_fifth_equation',s.expand(A2*fa**2+f1*fa+f0).coeff(h,5).subs(a[5],a5))
    Z=s.symbols('Z'); zero('incoming_moment_fifth',s.diff(Z+Z**3,Z,5))
    zero('incoming_energy_fifth',s.diff(Z**2+2*Z**4+Z**6,Z,5)/s.factorial(5)-6*Z)
    alpha,v,xi=s.symbols('alpha v xi',positive=True)
    zero('Gamma_fifth_integrand',s.diff((1+xi*v)**(-alpha),xi,5)+s.rf(alpha,5)*v**5*(1+xi*v)**(-alpha-5))
    return out


def branch_fixtures():
    with mp.workdps(85):
        c=MPIntervalContext(); c.dps=75; radius=mp.mpf('1e-65'); checks={}
        box=lambda v:c.mpf([v-radius,v+radius])
        p,q,K=map(mp.mpf,('.14','.12','.02'))
        b1=lambda z:mp.exp(z/3)+mp.mpf('.1')*z*z
        b2=lambda z:mp.mpf('.03')/(1+z*z)
        def branch(z):
            b,d=b1(z),b2(z); A=K*(1+q*p*p); B=1-q*p-2*K*q*p*b; C=q*b+K*q*b*b-d
            x=-2*C/(B+mp.sqrt(B*B-4*A*C)); return x,b-p*x
        A2=mp.mpf('.75'); A1=lambda z:mp.mpf('.04')*(z+z**3)
        A0=lambda z:-mp.mpf('1.2')-mp.exp(-z*z)+mp.mpf('.03')*(z+z**3)**2
        root=lambda z:(-A1(z)+mp.sqrt(A1(z)**2-4*A2*A0(z)))/(2*A2)
        for Z in map(mp.mpf,('0','.5')):
            prior=[IntervalTaylor(c,[box(v) for v in mp.taylor(lambda z:branch(z)[j],Z,4)]) for j in (0,1)]
            out,_=angular_fifth(c,prior,c.mpf(mp.taylor(b1,Z,5)[5]),c.mpf(mp.taylor(b2,Z,5)[5]),c.mpf(p),c.mpf(q),c.mpf(K))
            for j in (0,1):
                label='angular_branch'+str(j)+'_Z'+str(Z)
                contains(label,out[j][5],mp.taylor(lambda z:branch(z)[j],Z,5)[5]); checks[label]=True
            ap4=IntervalTaylor(c,[box(v) for v in mp.taylor(root,Z,4)])
            f1=IntervalTaylor(c,[box(v) for v in mp.taylor(A1,Z,5)])
            f0=IntervalTaylor(c,[box(v) for v in mp.taylor(A0,Z,5)])
            ap=selected_fifth(c,ap4,c.mpf(A2),f1,f0,2*c.mpf(A2)*ap4[0]+f1[0])
            label='selected_positive_root_Z'+str(Z)
            contains(label,ap[5],mp.taylor(root,Z,5)[5]); checks[label]=True
        return dict(checks=checks,actual_Md40_source_admission=False,passed=True)


def flatten_fixtures():
    with mp.workdps(75):
        c=MPIntervalContext(); c.dps=65; Z=mp.mpf('.5'); Q=1+Z*Z; u=2*Z/Q; v=1/Q
        def sigma(t):
            if t<=0:return mp.mpf(0)
            if t>=1:return mp.mpf(1)
            A,B=mp.exp(-1/t**2),mp.exp(-1/(1-t)**2); return A/(A+B)
        inverse5=lambda p:(2/Q)**p*(-mp.rf(p,3)*u*v*v/2+mp.rf(p,4)*u**3*v/6-mp.rf(p,5)*u**5/120)
        falling=lambda p,n:mp.fprod(p-i for i in range(n))
        power5=lambda p:(Q/2)**p*(falling(p,3)*u*v*v/2+falling(p,4)*u**3*v/6+falling(p,5)*u**5/120)
        r=SimpleNamespace(rate=c.mpf('.975'),cells=256,angular=SimpleNamespace(Xv=c.mpf('1.125')))
        expected=-mp.mpf('1.125')*mp.exp(-mp.mpf('97.5'))*inverse5(1)
        expected-=mp.quad(lambda t:mp.exp(-mp.mpf('.975')*(100-t))*inverse5(1-sigma(t/100)),[0,25,50,75,90,100])
        contains('preheat_fifth',preheat_fifth(c,c.mpf(Z),r),expected)
        q=IntervalTaylor(c,[c.mpf(Q),c.mpf(2*Z),1,0,0,0]); logq=log_taylor(q)-c.ln(2)
        flat=q*0; mu=mp.mpf('.025'); cells=512
        for i in range(cells):
            left=c.mpf(100)*i/cells; right=c.mpf(100)*(i+1)/cells
            t=c.mpf([endpoints(left)[0],endpoints(right)[1]])
            flat+=(logq*(2*stable_sigma(c,t/100)[0])).exp()*((right-left)*c.exp(-2*c.mpf(mu)*t))
        expected=mp.quad(lambda t:mp.exp(-2*mu*t)*power5(2*sigma(t/100)),[0,25,50,75,90,100])
        contains('energy_flatten_fifth',flat[5],expected)
        contains('log_Q_fifth',log_taylor(q)[5],mp.taylor(lambda z:mp.log(1+z*z),Z,5)[5])
        return dict(preheat_fifth=True,energy_flatten_fifth=True,log_Q_fifth=True,
            actual_Md40_source_admission=False,passed=True)


def gamma_fixtures():
    with mp.workdps(55):
        c=MPIntervalContext(); c.dps=70; a,S,Z,t=map(mp.mpf,('.1','.002','.5','3'))
        H=lambda xi:xi**(-1-a)*mp.hyperu(1+a,2,1/xi)
        D=lambda z:(1-H(2*(1-z*z)*S*mp.exp(-t)))/(a*S)
        actual=gamma_deficit_jets(c,c.mpf(Z),c.mpf(a),c.mpf(S),5,c.mpf(t))
        contains('true_Gamma_deficit_fifth',actual[5],mp.taylor(D,Z,5)[5])
        tails=integrated_gamma_tails(c,c.mpf(Z),c.mpf(a),c.mpf(S),5)
        def theta(z):
            xi=2*(1-z*z)*S*mp.exp(-3)
            return mp.exp(3*(1-a))*(xi**(-1-a)*mp.hyperu(1+a,3,1/xi)/(1-a)-1/(1-a))/S
        contains('true_Gamma_theta_tail_fifth',tails['theta'][5],mp.taylor(theta,Z,5)[5])
        cache={}
        def square5(t):
            if t in cache:return cache[t]
            xi=2*(1-Z*Z)*S*mp.exp(-t)
            delta=[mp.mpf(0),-4*Z*S*mp.exp(-t),-2*S*mp.exp(-t),0,0,0]
            coeffs=[mp.mpf(0)]*6; power=[mp.mpf(1),0,0,0,0,0]
            for m in range(6):
                hm=(-1)**m*mp.rf(a,m)*mp.rf(1+a,m)*xi**(-1-a-m)*mp.hyperu(1+a+m,2,1/xi)/math.factorial(m)
                coeffs=[coeffs[j]+hm*power[j] for j in range(6)]
                power=[sum(power[i]*delta[j-i] for i in range(j+1)) for j in range(6)]
            cache[t]=-sum(coeffs[i]*coeffs[5-i] for i in range(6))/(2*a*S)
            return cache[t]
        # Independently bound |J_5-A_5 exp(-t)| <= C*S exp(-2t).
        # Derive composition coefficients from a symbolic polynomial,
        # rather than using the production Gamma error-weight routine.
        h=s.symbols('h'); dz=-2*s.Rational(1,2)*h-h*h; d=1-Z*Z
        B=(1+a)**2*(2+a); C=[2*B*d*d,4*B*d*2*abs(Z),4*B*d,0,0,0]
        for m in range(2,6):
            bound=mp.rf(1+a,m-1)*mp.rf(1+a,m)*2**m*S**(m-2)/math.factorial(m)
            polynomial=s.Poly(s.expand(dz**m),h)
            for j in range(6): C[j]+=bound*abs(mp.mpf(str(polynomial.nth(j))))
        lead=[2*(1+a)*d,-4*(1+a)*Z,-2*(1+a),0,0,0]
        constant=C[5]+a/2*sum((abs(lead[i])+S*C[i])*(abs(lead[5-i])+S*C[5-i]) for i in range(6))
        if constant>=mp.mpf('1e6'):raise ArithmeticError('Independent omitted-tail constant not bounded')
        checks=dict(true_Gamma_deficit_fifth=True,true_Gamma_theta_tail_fifth=True)
        limit=mp.mpf(30)
        for name,rate,factor in (('pressure',1+2*a,1),('energy',2*a,2)):
            ref=mp.quad(lambda tt:mp.exp(-rate*tt)*square5(tt),[3,8,16,limit])
            error=mp.mpf('1e6')*S*mp.exp(-(rate+2)*limit)/(rate+2)
            contains('true_Gamma_'+name+'_tail_fifth',tails[name][5],factor*ref,factor*error)
            checks['true_Gamma_'+name+'_tail_fifth']=True
        return dict(checks=checks,independent_infinite_tail_remainder_constant=str(constant),
            actual_Md40_source_admission=False,passed=True)


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Fifth-order dependency changed: '+name)
    c=MPIntervalContext(); c.dps=240; read=lambda v:read_interval(c,v)
    angular=json.loads((HERE/'lei_ren_part1_paper_compliant_angular_high_jets.json').read_bytes())
    energy=json.loads((HERE/'lei_ren_part1_paper_compliant_future_energy_high_jets.json').read_bytes())
    axial=json.loads((HERE/'lei_ren_part1_paper_compliant_axial_high_jets.json').read_bytes())
    mu_name='lei_ren_part1_paper_compliant_outer_pulse_map.json'
    if mu_name not in hashes:raise ValueError('Actual source mu receipt not bound')
    mu=read(json.loads((HERE/mu_name).read_bytes())['mu'])
    prefix_checks=0
    for point,pa,pe,ps in zip(r['samples']+[r['whole_Z']],angular['samples']+[angular['whole_Z_C4_coefficients']],energy['samples']+[energy['whole_Z_C4']],axial['samples']+[axial['whole_Z_C4']]):
        a,e,j=point['angular'],point['energy'],point['selected']
        pairs=list(zip(a['scaled_coefficient_Taylor'],pa['scaled_coefficient_Taylor']))
        pairs+=list(zip(a['physical_coefficient_Taylor'],pa['physical_coefficient_Taylor']))
        pairs+=[(a['preheat_difference_scaled_Taylor'],pa['preheat_difference_scaled_Taylor']),
                (e['complete_future_energy_Taylor'],pe['complete_future_energy_Taylor']),
                (j['selected_ap_Taylor'],ps['selected_ap_Taylor'])]
        pairs+=list(zip(j['selected_scaled_end_coefficient_Taylor'],ps['selected_scaled_end_coefficient_Taylor']))
        for current,prior in pairs:
            if len(current['coefficients'])!=6:raise ValueError('Fifth-order function lost')
            for n in range(5):
                lo,hi=endpoints(read(current['coefficients'][n])); pl,ph=endpoints(read(prior['coefficients'][n]))
                if not lo <= pl <= ph <= hi:raise ArithmeticError('Admitted C4 prefix changed')
                prefix_checks+=1
        det=read(a['derivative_jacobian_determinant'])
        if endpoints(det)[0]<=0<=endpoints(det)[1]:raise ArithmeticError('Actual fifth angular inverse lost')
        if endpoints(read(a['fifth_physical_derivative_sum_bound'])-mu/100)[1]>=0:
            raise ArithmeticError('Actual uniform fifth angular derivative bound lost')
        if endpoints(read(j['positive_root_derivative_denominator']))[0]<=0:raise ArithmeticError('Actual fifth selected inverse lost')
        controls=j['selected_scaled_end_coefficient_Taylor']
        if not endpoints(read(controls[0]['coefficients'][0]))[1]<0<endpoints(read(controls[1]['coefficients'][0]))[0]:
            raise ArithmeticError('Actual selected fifth branch end signs lost')
        for key in ('full_pulse_C4_installed','full_outer_C4_certified','whole_outer_cone_certified','temporal_recursion'):
            if j[key]:raise ValueError('Fifth source layer overclaimed')
    for jet in r['samples'][1]['angular']['physical_coefficient_Taylor']:
        if endpoints(read(jet['coefficients'][5]))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Fifth angular odd parity lost')
    if endpoints(read(r['samples'][1]['energy']['complete_future_energy_Taylor']['coefficients'][5]))!=(mp.mpf(0),mp.mpf(0)):
        raise ArithmeticError('Fifth future energy parity lost')
    proof=identities(); branch=branch_fixtures(); flatten=flatten_fixtures()
    print('Fifth implicit/root/flatten derivatives and preserved actual C4 prefixes PASS',flush=True)
    gamma=gamma_fixtures()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(); hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        fifth_equation_identities=proof,independent_branch_fixtures=branch,independent_flatten_fixtures=flatten,
        independent_Gamma_fixtures=gamma,actual_C4_prefix_comparisons=prefix_checks,
        angular_coefficient_C5_available=True,complete_future_corrected_energy_C5_available=True,
        actual_selected_ap_c1_c2_C5_available=True,all_passed=True,full_pulse_C4_installed=False,
        full_outer_C4_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('Separate fifth-order same-source angular/energy/selected functions: independent true Gamma fifth derivatives PASS',flush=True)
    return out


if __name__=='__main__':run()
