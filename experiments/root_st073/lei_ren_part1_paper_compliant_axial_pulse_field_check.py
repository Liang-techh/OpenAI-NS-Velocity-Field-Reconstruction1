"""Independent paper recovery/ODE identities and selected pulse interfaces.

Numerical shape and finite-mu fixtures supplement the exact identities;
neither fixture claims a global cone or a radial Z derivative.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_axial_pulse_field import (
    gp, decay_integral, backward_bump_weights)

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_axial_pulse_field.json'


def exact_identities():
    passed={}
    def zero(name,expression):
        if s.simplify(expression)!=0:raise ArithmeticError('Pulse identity failed: '+name)
        passed[name]=True
    t,R,Z,mu,delta=s.symbols('t R Z mu delta',real=True)
    a=s.Rational(1,2)+mu;U=s.exp(-a*t);B=s.Function('B')(t)
    Mz=s.Function('Mz')(t);Mm=s.Function('Mm')(t);Mt=s.Function('Mt')(t);Me=s.Function('Me')(t)
    n1=s.exp(t)*U;n2=s.sqrt(2)*s.exp(s.Rational(3,2)*t)*U**2
    nt=s.sqrt(2)*s.exp(s.Rational(3,2)*t)*U;ne=s.exp(t)*U**2
    zero('normalized_Mz_ODE',s.diff(Mz/n1,t).subs(s.diff(Mz,t),n1*B)-(B-(s.Rational(1,2)-mu)*Mz/n1))
    zero('normalized_mixed_ODE',s.diff(Mm/n2,t).subs(s.diff(Mm,t),n2*B)-(B-(s.Rational(1,2)-2*mu)*Mm/n2))
    zero('normalized_angular_ODE',s.diff(Mt/nt,t).subs(s.diff(Mt,t),nt)-(1-(1-mu)*Mt/nt))
    zero('normalized_energy_ODE',s.diff(Me/ne,t).subs(s.diff(Me,t),ne*(B**2-s.Rational(1,2)))-(B**2-s.Rational(1,2)+2*mu*Me/ne))
    # The paper R variable is not the physical cylindrical radius.
    M=s.Function('M')(R,Z);Uz=s.Function('Uz')(R,Z);L=1-delta*Z**2;d=1-Z**2
    N=(2*Z*R*Uz-(1-delta)*Z*M-d*s.diff(M,Z))/L
    differentiated=s.diff(N,R).subs(s.diff(M,R,Z),s.diff(Uz,Z)).subs(s.diff(M,R),Uz)
    zero('paper_3_8_from_integrated_3_9',L*differentiated-((1+delta)*Z*Uz-d*s.diff(Uz,Z)+2*Z*R*s.diff(Uz,R)))
    m,mz,b,theta=s.symbols('m mz b theta',real=True);q=1+Z**2
    normalized=(2*Z*b-(1-delta)*Z*m-d*(mz-2*Z*m/q))/L
    recovered=(2*Z*R*(theta*b)-(1-delta)*Z*(R*theta*m)-d*(R*theta*(mz-2*Z*m/q)))/(L*s.sqrt(2*R))
    zero('normalized_radial_q_recovery',recovered-normalized*s.sqrt(R/2)*theta)
    P=s.Function('P')(t)
    zero('pressure_radial_balance',2*s.diff(P,t).subs(s.diff(P,t),U**2/2)-U**2)
    z,v,lam,E=s.symbols('z v lam E',real=True);G=s.Function('G')
    mb=-E*s.Integral(s.exp(lam*(v-z))*G(v),(v,z,0))
    zero('backward_moment_FTC',s.diff(mb,z)-(E*G(z)-lam*mb))
    ev=s.symbols('ev',real=True)
    eb=s.exp(2*mu*z)*ev+s.Integral(s.exp(2*mu*(z-v))/2,(v,z,0))-E**2*s.Integral(s.exp(2*mu*(z-v))*G(v)**2,(v,z,0))
    zero('backward_energy_FTC',s.diff(eb,z)-(E**2*G(z)**2-s.Rational(1,2)+2*mu*eb))
    D,K,cj=s.symbols('D K cj',real=True)
    eg=s.exp(-2*D)*ev+(1-s.exp(-2*D))/(4*mu)-s.exp(-2*D)*E**2*K*cj**2
    zero('inactive_gap_energy_ODE',-mu*s.diff(eg,D)-(-s.Rational(1,2)+2*mu*eg))
    x,row,LS,lu,lmu=s.symbols('xi row LS logu logmu',real=True)
    logE=-1/mu-3*LS+2*lu-s.log(6)/2-2*lmu
    reduced=(D/2-1)/mu-row*D-3*LS+2*lu-s.log(6)/2-2*lmu
    zero('gap_reduced_log_no_inverse_mu_cancellation',logE+(s.Rational(1,2)-row*mu)*D/mu-reduced)
    center,gram=s.symbols('center gram',real=True)
    zero('end_energy_Rp_to_Rv_units',s.exp(26)*s.exp(-26-2*mu*center)*gram-s.exp(-2*mu*center)*gram)
    zero('gap_chart_overlap_coordinate',13+mu*(-s.Rational(9,10)/mu)-s.Rational(121,10))
    pr=1+2*mu
    zero('end_pressure_segmented_log',-pr*(13/mu+z)-(-13*pr/mu-pr*z))
    # Forward/main and backward/gap energy agree by the selected equation.
    ap,Kmain,e0,endenergy=s.symbols('ap Kmain e0 endenergy',real=True)
    chosen=(1-s.exp(-26))/4-mu*e0+mu*s.exp(-26)*ev-mu*s.exp(-26)*endenergy
    forward=s.exp(22)*(e0+ap**2*Kmain/mu-(1-s.exp(-22))/(4*mu))
    backward=s.exp(-4)*ev+(1-s.exp(-4))/(4*mu)-s.exp(-4)*endenergy
    zero('main_gap_energy_selected_identity',(forward-backward).subs(ap**2*Kmain,chosen))
    return passed


def kernel_fixtures(normalization):
    # Reference quadrature runs at 90 digits, above the directed fixture
    # context so rounding of an exact value such as .49 is also enclosed.
    c=MPIntervalContext();c.dps=70;results={}
    with mp.workdps(90):
        norm=mp.quad(lambda r:mp.exp(-1/(1-r*r)) if abs(r)<1 else 0,[-1,0,1])
        def sig(x):
            if x<=0:return mp.mpf(0)
            if x>=1:return mp.mpf(1)
            odds=-1/x**2+1/(1-x)**2
            if odds>=0:return 1/(1+mp.exp(-odds))
            e=mp.exp(odds);return e/(1+e)
        def shape(x):
            if x<=0 or x>=11:return mp.mpf(0)
            P=x-mp.mpf('.01') if x>=mp.mpf('.02') else mp.quad(lambda v:sig(50*v),[0,x])
            return P*sig(11-x)
        def contains(interval,value,label):
            lo,hi=endpoints(interval)
            if not lo<=value<=hi:raise ArithmeticError('Kernel fixture outside enclosure: '+label)
            results[label]=True
        for x in (mp.mpf('.01'),mp.mpf('.5'),mp.mpf('10.5'),mp.mpf(11)):
            contains(gp(c,c.mpf(x))['value'],shape(x),'gp_'+str(x))
        # Tiny positive lengths must not be evaluated by 1-exp(-tiny).
        length=c.mpf('1e-1000');tiny=decay_integral(c,2,length)
        if endpoints(tiny)[0]<=0 or not endpoints(tiny/length)[0]>.99:
            raise ArithmeticError('Sub-precision positive integral lost')
        results['positive_sub_precision_decay_integral']=True
        normal=read_interval(c,normalization);mu=mp.mpf('.04');ell=mp.mpf('.15')
        for start in (mp.mpf('-.4'),mp.mpf(0),mp.mpf('.075'),mp.mpf('.2')):
            actual=backward_bump_weights(c,c.mpf(mu),normal,c.mpf(start),cells=256)
            lo=max(-ell,min(ell,start))
            beta=lambda v:mp.exp(-1/(1-(v/ell)**2))/(ell*norm) if abs(v)<ell else mp.mpf(0)
            points=[lo]+([0] if lo<0 else [])+[ell] if lo<ell else None
            refs=[mp.quad(lambda v:mp.exp((mp.mpf('.5')-row*mu)*v)*beta(v),points) if points else mp.mpf(0) for row in (1,2)]
            refs.append(mp.quad(lambda v:mp.exp(-2*mu*v)*beta(v)**2,points) if points else mp.mpf(0))
            for row,(interval,value) in enumerate(zip(actual,refs)):
                contains(interval,value,'future_beta_'+str(start)+'_'+str(row))
    return results


def finite_mu_interface_fixture():
    """Independent integrals of the actual gp and normalized end bumps."""
    with mp.workdps(65):
        mu=mp.mpf('.04');tv=13/mu;ap=mp.mpf('1.01');mi=(mp.mpf('.03'),mp.mpf('-.02'))
        ell=mp.mpf('.15');centers=(-3,-1);ev=mp.mpf('.4')
        normal=mp.quad(lambda r:mp.exp(-1/(1-r*r)) if abs(r)<1 else 0,[-1,0,1])
        def sigma(x):
            if x<=0:return mp.mpf(0)
            if x>=1:return mp.mpf(1)
            odds=-1/x**2+1/(1-x)**2
            if odds>=0:return 1/(1+mp.exp(-odds))
            e=mp.exp(odds);return e/(1+e)
        def g(x):
            if x<=0 or x>=11:return mp.mpf(0)
            P=x-mp.mpf('.01') if x>=mp.mpf('.02') else mp.quad(lambda v:sigma(50*v),[0,x])
            return P*sigma(11-x)
        beta=lambda v:mp.exp(-1/(1-(v/ell)**2))/(ell*normal) if abs(v)<ell else mp.mpf(0)
        weights=lambda lam:mp.quad(lambda v:mp.exp(lam*v)*beta(v),[-ell,0,ell])
        lambdas=[mp.mpf('.5')-row*mu for row in (1,2)]
        Pi=[mp.quad(lambda x:mp.exp(lam*(x-13)/mu)*g(x),[0,.02,1,10,10.5,11])/mu for lam in lambdas]
        A=mp.matrix([[mp.exp(lam*center)*weights(lam) for center in centers] for lam in lambdas])
        rhs=mp.matrix([-mi[row]*mp.exp(-lambdas[row]*tv)-ap*Pi[row] for row in range(2)])
        coeff=mp.lu_solve(A,rhs)
        errors=[]
        for offset in (-2/mu,-4,-2,mp.mpf(0)):
            for row,lam in enumerate(lambdas):
                base=(mi[row]*mp.exp(-lam*tv)+ap*Pi[row])*mp.exp(-lam*offset)
                past=mp.mpf(0);future=mp.mpf(0)
                for cj,center in zip(coeff,centers):
                    cut=max(-ell,min(ell,offset-center))
                    if cut>-ell:
                        past+=cj*mp.exp(lam*(center-offset))*mp.quad(lambda v:mp.exp(lam*v)*beta(v),[-ell]+([0] if cut>0 else [])+[cut])
                    if cut<ell:
                        future-=cj*mp.exp(lam*(center-offset))*mp.quad(lambda v:mp.exp(lam*v)*beta(v),[cut]+([0] if cut<0 else [])+[ell])
                error=abs(base+past-future)/(1+abs(future));errors.append(error)
        Kmain=mp.quad(lambda x:mp.exp(-2*x)*g(x)**2,[0,.02,1,10,10.5,11])
        gram=mp.quad(lambda v:mp.exp(-2*mu*v)*beta(v)**2,[-ell,0,ell])
        endenergy=sum(cj**2*mp.exp(-2*mu*center)*gram for cj,center in zip(coeff,centers))
        e0=mp.exp(-26)*ev+(1-mp.exp(-26))/(4*mu)-ap**2*Kmain/mu-mp.exp(-26)*endenergy
        forward=mp.exp(22)*(e0+ap**2*Kmain/mu-(1-mp.exp(-22))/(4*mu))
        backward=mp.exp(-4)*ev+(1-mp.exp(-4))/(4*mu)-mp.exp(-4)*endenergy
        errors.append(abs(forward-backward)/(1+abs(backward)))
        if max(errors)>mp.mpf('1e-35'):raise ArithmeticError('Finite-mu directed interface normalization fixture failed')
        return dict(mu=str(mu),moment_interfaces=8,main_gap_energy_interfaces=1,
                    max_relative_error=mp.nstr(max(errors),20),passed=True)


def run():
    r=json.loads((HERE/NAME).read_bytes());hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Selected pulse dependency changed: '+name)
    identities=exact_identities();c=MPIntervalContext();c.dps=240;get=lambda v:read_interval(c,v)
    whole=r['whole_Z_terminal'];jetkeys=('Uz_over_Utheta','Mz_over_R_Utheta',
        'Mtheta_z_over_sqrt2_R_3half_Utheta_squared','Mztheta_over_R_Utheta_squared')
    for key in jetkeys:
        if len(whole[key]['coefficients'])!=2:raise ArithmeticError('Terminal axial C1 lost')
    for key in jetkeys[:3]:
        for v in whole[key]['coefficients']:
            if endpoints(get(v))!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Whole-Z terminal linear closure failed')
    if endpoints(get(whole['Ur_over_sqrt_R_over_2_Utheta']))!=(mp.mpf(0),mp.mpf(0)):
        raise ArithmeticError('Terminal radial recovery not zero')
    energy=whole[jetkeys[3]]['coefficients'];target=whole['positive_terminal_future_energy_Taylor']['coefficients']
    if endpoints(get(energy[0]))[0]<=0:raise ArithmeticError('Positive terminal half-future-energy target lost')
    for actual,expected in zip(energy,target):
        alo,ahi=endpoints(get(actual));elo,ehi=endpoints(get(expected))
        if not alo<=elo<=ehi<=ahi:raise ArithmeticError('Terminal energy no longer encloses exact selected target')
    # Directed overlap is an additional enclosure check. Exact chart identity
    # follows from the selected moment equations, not from zero-box overlap.
    interfaces={}
    def overlap(left,right,label):
        for key in jetkeys:
            for v,w in zip(left[key]['coefficients'],right[key]['coefficients']):
                alo,ahi=endpoints(get(v));blo,bhi=endpoints(get(w))
                if max(alo,blo)>min(ahi,bhi):raise ArithmeticError('Selected pulse interface is disjoint: '+label+' '+key)
        interfaces[label]=True
    overlap(r['main_samples'][-1],r['gap_samples'][0],'main_gap_xi11')
    overlap(r['gap_samples'][2],r['gap_end_samples'][0],'gap_chart_overlap_xi12_1')
    overlap(r['gap_end_samples'][-1],r['end_samples'][0],'gap_end_s_minus4')
    for rv,rp in zip(r['full_end_energy_weights_in_Rv_units'],r['selector_end_energy_weights_in_Rp_units']):
        a=get(rv);b=get(rp)*c.exp(26)
        if max(endpoints(a)[0],endpoints(b)[0])>min(endpoints(a)[1],endpoints(b)[1]):
            raise ArithmeticError('End-energy physical normalization mismatch')
    samples=r['main_samples']+r['entrance_samples']+r['gap_samples']+r['gap_end_samples']+r['end_samples']+[whole]
    for point in samples:
        if point['Ur_Z_available'] or point['whole_outer_cone_certified'] or point['temporal_recursion']:
            raise ValueError('Partial pulse scope overclaimed')
        for key in jetkeys:
            if len(point[key]['coefficients'])!=2:raise ArithmeticError('Partial moment C1 lost')
        if not point['pressure']['exact_pressure_source_unchanged']:raise ValueError('Preheat P0 replaced')
        p=point['pressure']
        for mpv,p0v,pv in zip(p['Mp_over_Pstar_squared']['coefficients'],p['P0_over_Pstar_squared']['coefficients'],p['P_over_Pstar_squared']['coefficients']):
            expected=get(mpv)+get(p0v);actual=get(pv)
            if not endpoints(actual)[0]<=endpoints(expected)[0]<=endpoints(expected)[1]<=endpoints(actual)[1]:
                raise ArithmeticError('Pressure restoration is not P0+Mp')
    for point in [r['gap_end_samples'][-1],r['end_samples'][0]]:
        if endpoints(get(point['pressure']['retained_swirl_pressure_decay_log_parts']['finite_offset']))[0]<=0:
            raise ArithmeticError('Finite pressure offset lost in absolute inverse-mu log')
    if not r['full_O4_interval_composed'] or not r['inactive_gap_partial_fields_callable']:
        raise ValueError('O.4 gap composition missing')
    if r['full_outer_five_moment_match'] or r['whole_outer_cone_certified'] or r['global_admissible_stress_lift_constructed'] or r['temporal_recursion']:
        raise ValueError('Global completion overclaimed')
    kernels=kernel_fixtures(r['bump_normalization']);fixture=finite_mu_interface_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],
        implicit_source_sha256=r['implicit_source_sha256'],independent_pulse_identities=identities,
        directed_chart_interfaces=interfaces,independent_kernel_fixtures=kernels,finite_mu_fixture=fixture,
        whole_Z_terminal_linear_zero_and_positive_energy_checked=True,unchanged_preheat_pressure_checked=True,
        radial_recovery_paper_3_8_3_9_checked=True,radial_Z_derivative_certified=False,
        actual_O4_partial_pulse_field_installed=True,all_passed=True,
        full_outer_five_moment_match=False,whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Selected O.4 pulse:',len(identities),'independent identities;',len(interfaces),'directed interfaces; whole-Z terminal, kernel and finite-mu gates PASS',flush=True)
    return result


if __name__=='__main__':run()
