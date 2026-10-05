"""Focused original entrance cone checks and independent startup integrals."""
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_entrance_cone import (
    current_sources,entrance_cone_identities,whole_entrance_bounds,entrance_log_envelopes,
    full_stress_ratio,DOMAIN,PREFIX,TAIL_DOMAIN,ADMISSIONS,FALSE_FLAGS,source_precision)
from lei_ren_part1_paper_compliant_pulse_entrance_similarity_C4 import (
    CompliantPulseEntranceSimilarityC4,compiled_entrance_exporter)
from lei_ren_part1_paper_compliant_pulse_main_exit_fixture_integrals import run as fixture_integrals
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE=Path(__file__).parent


def independent_original_entrance_fixture():
    """Full original stresses, exact power memory and two-IBP startup identity."""
    with mp.workdps(90):
        c=MPIntervalContext();c.dps=120
        mu,delta=mp.mpf('.001'),mp.mpf('1e-8');r=1-mu;lam=mp.mpf('.5')-mu;p=1+2*mu
        Rp,Pstar,D0,Tw=mp.exp(200),mp.exp(-120),mp.mpf('.07'),mp.mpf(80)
        bp=mp.mpf('.5')+mu;pre=mp.mpf('.17');Xp=1/r+pre*mp.exp(-r*Tw)
        C=lambda z:1/(1+z*z)
        ap=lambda z:mp.mpf('1.02')+mu*mp.mpf('.02')*z*z
        incoming=(lambda z:mp.mpf('.0002')*(z+z**3),lambda z:mp.mpf('.00014')*(z+z**3))
        ein=lambda z:-mp.mpf('.015')+mp.mpf('.002')*(z*z+2*z**4+z**6)
        J0=lambda z:mp.mpf('.2')+mp.mpf('.05')*z*z
        P0=lambda z:-mp.mpf('.35')*C(z)**2+mp.mpf('.02')*z
        vals=fixture_integrals()['values'];Kpulse=mp.mpf(vals['Kpulse'])
        future=lambda z:mp.exp(26)/mu*(ap(z)**2*Kpulse-(1-mp.exp(-26))/4+mu*ein(z))+D0**2*J0(z)
        def sigma(t):
            if t<=0:return mp.mpf(0)
            if t>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/t**2-1/(1-t)**2))
        def sig1(t):
            if t<=0 or t>=1:return mp.mpf(0)
            v=sigma(t);return v*(1-v)*(2/t**3+2/(1-t)**3)
        cache={}
        def gp(x):
            key=mp.nstr(x,100)
            if key not in cache:
                cache[key]=mp.mpf(0) if x<=0 else mp.mpf('.01') if x==mp.mpf('.02') else mp.quad(
                    lambda a:sigma(50*a),[0,x/2,x],method='gauss-legendre')
            return cache[key]
        def kernel(x,i):
            if x==0:return mp.mpf(0)
            rate=mp.mpf('.5')-i*mu
            forward=mp.quad(lambda a:mp.exp(rate*a/mu)*sigma(50*a),[0,x/2,x])
            return (gp(x)-mp.exp(-rate*x/mu)*forward)/rate
        def energy(x):
            if x==0:return mp.mpf(0)
            if x==mp.mpf('.02'):return mp.mpf(vals['entrance_energy'])
            return mp.quad(lambda a:mp.exp(-2*a)*gp(a)**2,[0,x/2,x],method='gauss-legendre')
        checks={};worst=mp.mpf(0);tol=mp.mpf('1e-55');measured=[];nonzero=0
        def compare(name,a,b):
            nonlocal worst
            error=abs(a-b)/max(1,abs(a),abs(b));worst=max(worst,error)
            if error>tol:raise ArithmeticError('Independent original entrance mismatch '+name)
            checks[name]=True
        def enclose(name,bound,expected):
            nonlocal worst
            lo,hi=endpoints(bound);error=max(lo-expected,expected-hi,mp.mpf(0))/max(1,abs(expected))
            worst=max(worst,error)
            if error>tol:raise ArithmeticError('Independent source enclosure misses '+name)
            checks[name]=True
        jet=lambda fn,z:IntervalTaylor(c,[c.mpf(v) for v in mp.taylor(fn,z,5)])
        for x,z in ((mp.mpf(0),mp.mpf('.31')),(mp.mpf('.01'),mp.mpf(0)),(mp.mpf('.02'),mp.mpf('-.4'))):
            y=x/mu;R=Rp*mp.exp(y);B=Pstar*mp.exp(-bp*y);H=mp.exp(-r*y)
            g,gd=gp(x),sigma(50*x);kk=[kernel(x,i) for i in (1,2)];en=energy(x)
            for i in (1,2):
                rate=mp.mpf('.5')-i*mu
                rho=mp.mpf(0) if x==0 else mu/rate**2*mp.quad(
                    lambda a:mp.exp(-rate*(x-a)/mu)*50*sig1(50*a),[0,x/2,x])
                compare('actual_two_IBP/'+str(x)+'/'+str(i),kk[i-1],g/rate-mu*gd/rate**2+rho)
            compare('exact_pre_Tw_angular_memory/'+str(x),(Xp-1/r)*H,pre*mp.exp(-r*(Tw+y)))
            def data(zz):
                cc=C(zz);aa=ap(zz);ut=B*cc;uz=ut*aa*g
                m=[aa*kk[i-1]+mp.exp(-(mp.mpf('.5')-i*mu)*y)*incoming[i-1](zz) for i in (1,2)]
                ev=mp.exp(2*x)*(ein(zz)+aa*aa*en/mu-(1-mp.exp(-2*x))/(4*mu))
                pressure=B*B*(-cc*cc/(2*p)+mp.exp(-p*(13-x)/mu)*(P0(zz)+cc*cc/(2*p)))
                moments=dict(theta=mp.sqrt(2)*R**mp.mpf('1.5')*ut*(1/r+pre*mp.exp(-r*(Tw+y))),
                    z=R*ut*m[0],theta_z=mp.sqrt(2)*R**mp.mpf('1.5')*ut*ut*m[1],
                    z_theta=R*ut*ut*ev,p=mp.mpf(0))
                return dict(ut=ut,uz=uz,pressure=pressure,moments=moments)
            d=data(z)
            stress=evaluate_mp_stress(mp.log(R),z,delta,Utheta=d['ut'],Uz=d['uz'],
                Utheta_y=-bp*d['ut'],Utheta_Z=B*mp.diff(C,z),
                Uz_y=B*C(z)*ap(z)*(mu*gd-bp*g),Uz_Z=B*mp.diff(lambda zz:C(zz)*ap(zz),z)*g,
                moments=d['moments'],moments_Z={key:mp.diff(lambda zz:data(zz)['moments'][key],z) for key in d['moments']},
                P=d['pressure'],P_Z=mp.diff(lambda zz:data(zz)['pressure'],z),precision=90,radius_override=R,axial_override=z)
            obj=object.__new__(CompliantPulseEntranceSimilarityC4)
            obj.ctx=c;obj.mu=c.mpf(mu);obj.delta=c.mpf(delta);obj.Xp=c.mpf(Xp);obj.cells=128
            obj.ap=jet(ap,z);obj.incoming=[jet(f,z) for f in incoming];obj.incoming_energy=jet(ein,z)
            obj.future=jet(future,z);obj.J0=jet(J0,z);obj.P0=jet(P0,z);obj.Pin=c.mpf('.12')
            obj.logRp=c.mpf(mp.log(Rp));obj.logP=c.mpf(mp.log(Pstar));obj.logU=c.mpf(0)
            obj.finite=c.mpf(mp.log(D0)+1/mu);obj.exporter,_=compiled_entrance_exporter()
            def exact_kernel(ctx,mm,xx,rate,cells):
                index=0 if endpoints(rate)[0]>mp.mpf('.4985') else 1
                return dict(enclosure=c.mpf(kk[index]))
            obj.exporter.__globals__['partial_linear_kernel']=exact_kernel
            obj.exporter.__globals__['partial_future_energy']=lambda *args:c.mpf(Kpulse-en)
            # Exact named endpoints stay decimal points across MP contexts.
            packet=obj.entrance(c.mpf(z),c.mpf(mp.nstr(x,20)))
            normalize=mp.sqrt(R/2)*B
            totals={}
            for label,key in (('theta','T_theta'),('axial','T_z')):
                total=sum((part['full_stress_mixed3_coefficient_enclosures']['y0_Z0']*c.exp(
                    sum(part['exact_source_log_parts'].values(),c.mpf(0))-c.mpf(mp.log(normalize)))
                    for part in packet['full_meridional_stress_log_sectors'][label].values()),c.mpf(0))
                totals[label]=total
                enclose('full_original_stress/'+label+'/'+str(x),total,stress[key]/normalize)
            w=stress['T_z']/stress['T_theta'];a=2+2*mu;b=2*mu*ap(z)*gd-(1+2*mu)*ap(z)*g
            enclose('total_w_ratio/'+str(x),totals['axial']/totals['theta'],w)
            second=2*b*w+b*b/a+(a-2)*w*w
            if stress['T_theta']<=0 or a-b*w<=0 or 2-second<=0:raise ArithmeticError('Original fixture cone fails')
            F=d['ut']/mp.sqrt(2*R);dot=stress['T_theta']*F*(-a+b*w);cross=stress['T_theta']*F*(-b-a*w)
            kap=a+b*b/a
            compare('full_directional_factorization/'+str(x),
                (2*dot*dot-(kap-2)*cross*cross)/(stress['T_theta']*F)**2,(a*a+b*b)*(2-second))
            for nu in (mp.mpf('.01'),mp.mpf('.7')):
                positive=nu*mp.mpf('.83')**(-mp.mpf('2.5'))
                compare('fixed_nu_total_stress_ratio/'+str(x)+'/'+str(nu),
                    (positive*stress['T_z'])/(positive*stress['T_theta']),w)
            if x>0 and abs(d['uz'])>0:nonzero+=1
            measured.append(dict(xi=str(x),Z=str(z),Ttheta_normalized=mp.nstr(stress['T_theta']/normalize,20),
                bw=mp.nstr(b*w,20),directional_bracket=mp.nstr(2-second,20)))
            print('Independent original entrance cone fixture xi='+str(x)+' checked',flush=True)
        if nonzero!=2:raise ArithmeticError('Fixture axial input removed')
        return dict(checks=len(checks),identities=checks,tolerance=str(tol),maximum_normalized_error=mp.nstr(worst,30),
            nonzero_axial_input_fixtures=nonzero,measured_fixture_cones=measured,
            original_integral_and_forward_energy_independent_of_backward_exporter=True,
            exact_pre_Tw_memory_and_total_full_shear_stress_retained=True,
            fixtures_are_not_continuous_actual_source_proof=True)


@source_precision
def run():
    records,hashes,family=current_sources();name=PREFIX+'pulse_entrance_cone.json'
    raw=(HERE/name).read_bytes();record=json.loads(raw)
    proof=entrance_cone_identities(records);bounds=whole_entrance_bounds(records)
    if record['exact_source_cone_proof']!=encode(pack(proof)) or record['bounds']!=encode(pack(bounds)):
        raise ValueError('Current entrance cone source or bounds differ')
    if record['domain']!=DOMAIN or record['tail_domain']!=TAIL_DOMAIN:raise ValueError('Entrance domain differs')
    for path,digest in record['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Cone source changed: '+path)
    for flag in ADMISSIONS:
        if not record[flag]:raise ValueError('Missing entrance cone admission')
    for flag in FALSE_FLAGS:
        if record[flag]:raise ValueError('Cone scope overclaimed')
    c=MPIntervalContext();c.dps=300;count=0
    for category in ('positive_parameter_margins','positive_source_log_margins','positive_monotonicity_margins','positive_algebraic_margins'):
        for label,value in record['bounds'][category].items():
            lo,hi=endpoints(read_interval(c,value))
            if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Nonpositive whole entrance margin '+label)
            count+=1
    print('Continuous entrance source, every full stress sector and strict margins admitted',flush=True)
    oracle=independent_original_entrance_fixture()
    hashes.update(record['input_hashes']);hashes[name]=hashlib.sha256(raw).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    for stem in ('pulse_main_exit_fixture_integrals.py','pulse_main_exit_fixture_integrals.json'):
        path=PREFIX+stem;hashes[path]=hashlib.sha256((HERE/path).read_bytes()).hexdigest()
    result=dict(all_passed=True,actual_five_defect_family_sha256=family[0],implicit_source_sha256=family[1],
        domain=DOMAIN,tail_domain=TAIL_DOMAIN,input_hashes=hashes,positive_continuous_margins=count,
        source_cone_identities=len(proof['identities']),independent_original_entrance_fixture=oracle,
        all_original_stress_sectors_and_pre_Tw_memory_preserved=True,
        entrance_physical_Cartesian_receipt_consumed_without_rerun=True,
        **{flag:True for flag in ADMISSIONS},**{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('PASS continuous original entrance full-shear two-vector cone and main/tail composition; global/recursion pending',flush=True)
    return result


if __name__=='__main__':run()
