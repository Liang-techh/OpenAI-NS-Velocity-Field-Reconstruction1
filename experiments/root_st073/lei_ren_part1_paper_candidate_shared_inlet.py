"""One-source finite candidate moments, recovered pressure and inlet jets.

Amplitude factors cancel algebraically before normalized inlet division.
These are center jets; no terminal functional moment closure is asserted.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_candidate_combined_core_budget import DEFAULT_STATE
from lei_ren_part1_paper_interval_taylor import IntervalTaylor, constant
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def radial_product(a,b):
    zero=a[0]*0
    return [sum((a[i]*b[n-i] for i in range(len(a))
                 if 0<=n-i<len(b)),zero) for n in range(len(a)+len(b)-1)]


def integral(poly,s,weight=0):
    return sum((v*s**(n+weight+1)/(n+weight+1)
                for n,v in enumerate(poly)),poly[0]*0)


def diff(jet):
    return IntervalTaylor(jet.ctx,[(k+1)*jet[k+1] for k in range(jet.order)])


def finite_moments(phi,u,s):
    pp=radial_product(phi,phi)
    return dict(theta=2*integral(phi,s,1),z=integral(u,s),
                theta_z=2*integral(radial_product(phi,u),s,1),
                p=integral(pp,s),u_squared=integral(radial_product(u,u),s),
                weighted_phi_squared=integral(pp,s,1))


def normalized_inlet(phi,u,m,S,ell,p0,lam,s,z,dt,endpoints_override=None):
    ctx=z.ctx;eps=1/lam;L=1-z*z*dt;d=1-z*z
    phi_exit=sum((v*s**n for n,v in enumerate(phi)),phi[0]*0)
    u_exit=sum((v*s**n for n,v in enumerate(u)),u[0]*0)
    phi_s=sum((n*v*s**(n-1) for n,v in enumerate(phi) if n),phi[0]*0)
    u_s=sum((n*v*s**(n-1) for n,v in enumerate(u) if n),u[0]*0)
    if endpoints_override:
        phi_exit=endpoints_override['phi_exit'];u_exit=endpoints_override['u_exit']
        phi_s=endpoints_override['phi_s'];u_s=endpoints_override['u_s']
    pressure=p0+S*m['p']*eps
    mzt=m['u_squared']-S*m['weighted_phi_squared']*eps
    V=z*m['z']*(1-dt)+d*diff(m['z'])-s
    C=m['theta']*(1-dt/2)-z*(diff(m['theta'])+ell*m['theta'])*((1-dt)/2)
    C-=d*(diff(m['theta_z'])+ell*m['theta_z'])
    C+=z*m['theta_z']*(2*dt-1)
    itheta=(phi_exit*V/L+C/(2*L*s))*eps
    iz=(V*u_exit+(m['z']-z*diff(m['z']))*((1-dt)/2)
        +z*mzt*(2*dt)-d*diff(mzt)
        +(z*pressure*(2*(1+dt))-d*diff(pressure))*s)/(L*ctx.sqrt(2*s))
    ur=(2*z*s*u_exit-z*m['z']*(1-dt)-d*diff(m['z']))/(L*ctx.sqrt(2*s))
    return dict(phi_exit=phi_exit,u_exit=u_exit,pressure=pressure,mzt=mzt,
        itheta=itheta,iz=iz,ur=ur,ratio=itheta/phi_exit,
        ttheta=itheta+phi_s*(2*s),tz=iz+u_s*(ctx.sqrt(2*s)*lam))


def fixture():
    ctx=MPIntervalContext();ctx.dps=80
    with mp.workdps(120):
        phi=[constant(ctx,1,3),constant(ctx,1,3)]
        u=[constant(ctx,2,3),constant(ctx,-1,3)]
        s=ctx.mpf('.5');m=finite_moments(phi,u,s)
        x=mp.mpf('.5')
        expected=dict(theta=x*x+2*x**3/3,z=2*x-x*x/2,
            theta_z=2*x*x+2*x**3/3-x**4/2,p=x+x*x+x**3/3,
            u_squared=4*x-2*x*x+x**3/3,
            weighted_phi_squared=x*x/2+2*x**3/3+x**4/4)
        for name,value in expected.items():
            lo,hi=endpoints(m[name][0])
            if not lo<=value<=hi:raise AssertionError(('radial moment fixture',name))
        # Compare the normalized transfer with the established original
        # physical stress evaluator on an unrelated polynomial field.
        from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
        lam=ctx.mpf(10);eps=ctx.mpf('.1');dt=ctx.mpf('.01')
        z=IntervalTaylor.variable(ctx,ctx.mpf('.3'),3)
        out=normalized_inlet(phi,u,m,constant(ctx,1,3),constant(ctx,0,3),
            constant(ctx,-3,3),lam,s,z,dt)
        r=x/10;root=mp.sqrt(2*r);ph=1+x;uv=2-x
        moments={name:expected[name]*(mp.mpf('.01') if name in ('theta','theta_z') else mp.mpf('.1'))
                 for name in ('theta','z','theta_z','p')}
        moments['z_theta']=mp.mpf('.1')*expected['u_squared']-mp.mpf('.01')*expected['weighted_phi_squared']
        direct=evaluate_mp_stress(mp.log(r),'.3','.01',Utheta=root*ph,Uz=uv,
            Utheta_y=root*(ph/2+x),Utheta_Z=0,Uz_y=-x,Uz_Z=0,
            moments=moments,moments_Z={name:0 for name in moments},
            P=-3+moments['p'],P_Z=0,precision=120)
        for name,reference,factor in (('itheta','I_theta',1),('iz','I_z',mp.sqrt(mp.mpf('.1'))),
            ('ur','U_r',mp.sqrt(mp.mpf('.1'))),('ttheta','T_theta',1),
            ('tz','T_z',mp.sqrt(mp.mpf('.1')))):
            lo,hi=endpoints(out[name][0]);target=direct[reference]/factor
            if not lo<=target<=hi:raise AssertionError(('physical inlet fixture',name))
        return dict(independent_polynomial_integral_checks=6,original_physical_stress_checks=5,passed=True)


def run(state_name=DEFAULT_STATE,s_exit='4',budget_name='lei_ren_part1_paper_candidate_combined_core_budget.json',trace_name='lei_ren_part1_paper_candidate_inlet_trace.json',output_name=None):
    checks=fixture();base=Path(__file__).parent;path=base/state_name
    raw=path.read_bytes();state=json.loads(raw)
    for name,digest in state['source_hashes'].items():
        name=({'driver':state['target'].get('driver_file','lei_ren_part1_paper_candidate_gauge_core.py'),
               'pressure_input':state['target']['pressure_file']}).get(name,name)
        if hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest:
            raise AssertionError('Candidate dependency changed:'+name)
    ctx=MPIntervalContext();ctx.dps=state['target']['precision']
    with mp.workdps(ctx.dps+40):
        def read(v):return ctx.mpf([mp.make_mpf(tuple(v['lower_exact_mpf_tuple'])),
                                   mp.make_mpf(tuple(v['upper_exact_mpf_tuple']))])
        def jet(row):return IntervalTaylor(ctx,[read(v) for v in row[:4]])
        lam=ctx.mpf(state['target']['Lambda']);eps=1/lam;s=ctx.mpf(s_exit)
        if endpoints(s)[0]<=0 or endpoints(s)[1]>mp.mpf('4.1'):
            raise ValueError('Exit outside (0,4.1]')
        z=IntervalTaylor.variable(ctx,ctx.mpf(state['target']['Z']),3)
        dt=ctx.mpf(state['target']['delta']);L=1-z*z*dt;d=1-z*z
        phi=[jet(row)/lam**n for n,row in enumerate(state['A_rows'])]
        u=[jet(row)/lam**n for n,row in enumerate(state['Uz_rows'])]
        if min(p.order for p in phi+u)<3:raise ValueError('Need three axial derivatives')
        ell=jet(state['fixed_jets']['ell_Z_taylor'])
        S=jet(state['fixed_jets']['S_Z_taylor'])
        F0=read(state['fixed_jets']['F0_interval'])
        m=finite_moments(phi,u,s)
        # Completion is the prescribed radial pressure primitive of F².
        # No independent pressure fit or anchor change is used.
        inlet=normalized_inlet(phi,u,m,S,ell,jet(state['fixed_jets']['P0_Z_taylor']),lam,s,z,dt)
        # Integrate the proved uniform radial remainders, then transfer the
        # resulting analytic core enclosures through the same inlet formulas.
        budget_path=base/budget_name
        budget=json.loads(budget_path.read_text())
        if budget['state_sha256']!=hashlib.sha256(raw).hexdigest():
            raise AssertionError('Combined tail budget belongs to a different state')
        if not budget['exact_gauge_recursion_to_analytic_fixed_point_identity_audited']:
            raise AssertionError('Analytic coefficient identity is missing')
        for name,digest in budget['input_hashes'].items():
            if hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest:
                raise AssertionError('Combined budget dependency changed:'+name)
        tails={(r['component'],r['scaled_radial_order'],r['axial_order']):
               read(r['conditional_infinite_radial_tail']) for r in budget['derivative_budgets']}
        # The proof radius is4.1; smaller exits inherit the same bounds.
        Ephi=[tails['Phi',0,k] for k in range(4)]
        Eu=[tails['Psi',0,k]*eps for k in range(4)]
        def polynomial_bounds(poly):
            return [sum((ctx.mpf(max(abs(v) for v in endpoints(p[k])))*
                         mp.factorial(k)*s**n for n,p in enumerate(poly)),ctx.mpf(0))
                    for k in range(4)]
        Bphi=polynomial_bounds(phi);Bu=polynomial_bounds(u)
        def product_errors(Bf,Bg,Ef,Eg):
            return [sum(((Ef[j]*Bg[k-j]+Bf[j]*Eg[k-j]+Ef[j]*Eg[k-j])*math.comb(k,j)
                         for j in range(k+1)),ctx.mpf(0)) for k in range(4)]
        Epu=product_errors(Bphi,Bu,Ephi,Eu)
        Epp=product_errors(Bphi,Bphi,Ephi,Ephi)
        Euu=product_errors(Bu,Bu,Eu,Eu)
        moment_errors=dict(theta=[e*s*s for e in Ephi],z=[e*s for e in Eu],
            theta_z=[e*s*s for e in Epu],p=[e*s for e in Epp],
            u_squared=[e*s for e in Euu],weighted_phi_squared=[e*s*s/2 for e in Epp])
        def inflate(p,errors):
            coeffs=[]
            for k,error in enumerate(errors):
                upper=endpoints(error/math.factorial(k))[1]
                coeffs.append(p[k]+ctx.mpf([-upper,upper]))
            return IntervalTaylor(ctx,coeffs)
        enclosed_m={name:inflate(value,moment_errors[name]) for name,value in m.items()}
        # Exact same physical axis data; only the omitted radial orders differ.
        radial_phi=sum((n*v*s**(n-1) for n,v in enumerate(phi) if n),phi[0]*0)
        radial_u=sum((n*v*s**(n-1) for n,v in enumerate(u) if n),u[0]*0)
        overrides=dict(phi_exit=inflate(inlet['phi_exit'],Ephi),
            u_exit=inflate(inlet['u_exit'],Eu),
            phi_s=inflate(radial_phi,[tails['Phi',1,k] for k in range(3)]),
            u_s=inflate(radial_u,[tails['Psi',1,k]*eps for k in range(3)]))
        enclosed=normalized_inlet(phi,u,enclosed_m,S,ell,
            jet(state['fixed_jets']['P0_Z_taylor']),lam,s,z,dt,overrides)
        for name in ('ttheta','tz'):
            if any(not endpoints(v)[0]<=0<=endpoints(v)[1]
                   for v in enclosed[name].coefficients):
                raise AssertionError('Analytic core stress-zero consistency failed:'+name)
        trace_path=base/trace_name
        trace=json.loads(trace_path.read_text())
        if trace['state_sha256']!=hashlib.sha256(raw).hexdigest():
            raise AssertionError('Independent angular trace has different core data')
        direct=read(trace['finite_Ttheta_over_F0_integral'])
        lo,hi=endpoints(inlet['ttheta'][0]);dl,dh=endpoints(direct)
        if max(lo,dl)>min(hi,dh):raise AssertionError('Five-moment inlet and angular integral disagree')
        def pack(p):return [encode(v) for v in p.coefficients]
        report=dict(state_file=path.name,state_sha256=hashlib.sha256(raw).hexdigest(),
            accepted_schedule_sha256=state['accepted_schedule_sha256'],
            Lambda=state['target']['Lambda'],Z=state['target']['Z'],scaled_exit=s,
            radial_degree=state['completed_radial_order'],ordinary_axial_Taylor_order=3,
            inlet_axial_Taylor_order=2,
            moment_scales={'theta':'F0*epsilon²','z':'epsilon',
                'theta_z':'F0*epsilon²','z_theta':'epsilon','p':'F0²*epsilon'},
            normalized_moment_jets={name:pack(inlet['mzt'] if name=='z_theta' else m[name])
                                    for name in ('theta','z','theta_z','z_theta','p')},
            normalized_Phi_exit=pack(inlet['phi_exit']),Uz_exit=pack(inlet['u_exit']),
            recovered_physical_pressure_exit=pack(inlet['pressure']),
            I_theta_over_F0=pack(inlet['itheta']),D_Itheta_over_F=pack(inlet['ratio']),
            I_z_over_sqrt_epsilon=pack(inlet['iz']),Ur_over_sqrt_epsilon=pack(inlet['ur']),
            finite_Ttheta_over_F0=pack(inlet['ttheta']),finite_Tz_over_sqrt_epsilon=pack(inlet['tz']),
            analytic_normalized_moment_enclosures={name:pack(enclosed['mzt'] if name=='z_theta' else enclosed_m[name])
                for name in ('theta','z','theta_z','z_theta','p')},
            analytic_core_inlet_enclosures={name:pack(enclosed[name]) for name in
                ('phi_exit','u_exit','pressure','itheta','iz','ur','ratio','ttheta','tz')},
            normalized_moment_radial_remainder_derivative_bounds=moment_errors,
            combined_budget_sha256=hashlib.sha256(budget_path.read_bytes()).hexdigest(),
            input_hashes={name:hashlib.sha256((base/name).read_bytes()).hexdigest()
                for name in (Path(__file__).name,'lei_ren_part1_paper_interval_taylor.py',
                             'lei_ren_part1_paper_mp_stress.py',trace_path.name,budget_path.name)},
            independent_moment_and_angular_trace_enclosures_overlap=True,
            analytic_core_stress_jets_include_exact_zero=True,
            zero_containment_is_not_terminal_matching_certificate=True,
            F0_at_center=F0,fixture=checks,
            physical_amplitude_factor_cancelled_before_inlet_division=True,
            prescribed_pressure_primitive_P_R_equals_F_squared=True,
            pressure_anchor_changed=False,independent_pressure_fit_used=False,
            finite_and_analytic_core_enclosures_separately_reported=True,
            infinite_core_moment_remainders_enclosed=True,
            five_terminal_functional_moments_closed=False,
            whole_axis_inlet_certified=False,shared_collar_regenerated=False,
            temporal_recursion=False)
        output=base/output_name if output_name else Path(__file__).with_suffix('.json')
        output.write_text(json.dumps(encode(report),indent=2)+'\n',encoding='utf-8')
        print('Shared candidate five core moments and inlet:',state['completed_radial_order'],
              'radial orders; independent integral checks6',flush=True)
        print('finite normalized Tz interval',*[mp.nstr(v,18) for v in endpoints(inlet['tz'][0])],flush=True)
        return report


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--state-file',default=DEFAULT_STATE)
    parser.add_argument('--scaled-exit',default='4')
    parser.add_argument('--combined-budget-file',default='lei_ren_part1_paper_candidate_combined_core_budget.json')
    parser.add_argument('--trace-file',default='lei_ren_part1_paper_candidate_inlet_trace.json')
    parser.add_argument('--output-name')
    args=parser.parse_args()
    run(args.state_file,args.scaled_exit,args.combined_budget_file,args.trace_file,args.output_name)
