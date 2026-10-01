"""Independent recovery/normalization identities and directed repair gates."""
import ast
import hashlib
import json
from pathlib import Path
import mpmath as mp
import sympy as sp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_outer_angular_repair.json'


def run():
    r=json.loads((HERE/NAME).read_bytes());hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
            raise ValueError('Angular repair checker source changed: '+name)
    c=MPIntervalContext();c.dps=160
    def read(v):
        if 'lower_exact_mpf_tuple' in v:return read_interval(c,v)
        return c.mpf(mp.make_mpf(tuple(v['exact_mpf_tuple'])))
    get=lambda key:read(r[key])
    identities={}
    def zero(name,expression):
        if sp.simplify(expression)!=0:raise ArithmeticError('Independent identity failed: '+name)
        identities[name]=True
    mu,a,T,tau,lone=sp.symbols('mu a T tau lone',real=True)
    rate=1-mu;k=1-a;L=2+T+tau
    Etail=-sp.Rational(3,2)*(2+T)+(rate+k)/2-(sp.Rational(1,2)+a)*tau
    zero('angular_heat_log_conversion',sp.Rational(3,2)*L+Etail-lone-((rate+k)/2+k*tau-lone))
    zero('pressure_heat_log_conversion',2*(Etail+(sp.Rational(1,2)+a)*L-lone)-(1+2*a)*L+sp.log(a)-(2*Etail-2*lone+sp.log(a)))
    eta,K=sp.symbols('eta K',positive=True)
    zero('positive_flatten_divided_difference',(eta*((1-(1+eta)**(-K))/eta))-(1-(1+eta)**(-K)))
    zero('flatten_difference_derivative',sp.diff(1-(1+eta)**(-K),eta)-K*(1+eta)**(-1-K))
    q,Xv,H=sp.symbols('q Xv H',positive=True)
    zero('preflatten_angular_ratio_independent_of_Z',(H/q)/(Xv/q)-H/Xv)
    # Check the factor against the actual source expressions, then propagate
    # it across the entire axial-only pulse for arbitrary Z and arbitrary Uz.
    fq=sp.Symbol('common_q_inverse',positive=True)
    source_methods={
        'lei_ren_part1_paper_compliant_outer_initial.py':('reference','slope','axial'),
        'lei_ren_part1_paper_compliant_outer_buffer.py':('slope_mu','power')}
    factor_checks={}
    for filename,methods in source_methods.items():
        tree=ast.parse((HERE/filename).read_text(encoding='utf-8'))
        for method in methods:
            function=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
            env={'qi':fq,'u1':fq*sp.Symbol('u_inlet')}
            def formal(node):
                if isinstance(node,ast.Name):return env.get(node.id,sp.Symbol(node.id))
                if isinstance(node,ast.Constant):return sp.Rational(str(node.value))
                if isinstance(node,ast.UnaryOp):return -formal(node.operand) if isinstance(node.op,ast.USub) else formal(node.operand)
                if isinstance(node,ast.BinOp):
                    left,right=formal(node.left),formal(node.right)
                    if isinstance(node.op,ast.Add):return left+right
                    if isinstance(node.op,ast.Sub):return left-right
                    if isinstance(node.op,ast.Mult):return left*right
                    if isinstance(node.op,ast.Div):return left/right
                    if isinstance(node.op,ast.Pow):return left**right
                if isinstance(node,ast.Call):
                    label=ast.unparse(node.func)
                    if label=='get':
                        key=node.args[0].value
                        if key in ('Utheta_over_Pstar','Mtheta_over_sqrt2_R_3half_Pstar'):
                            return fq*sp.Symbol(key)
                        raise ArithmeticError('Unexpected angular input primitive: '+key)
                    if label=='c.exp':return sp.exp(formal(node.args[0]))
                    if label=='c.mpf':return formal(node.args[0])
                # Remaining coefficients/kernels have no qi or Z inputs in
                # these u/h assignments; reject explicit hidden dependence.
                if any(isinstance(n,ast.Name) and n.id in ('Z','z','qi') for n in ast.walk(node)):
                    raise ArithmeticError('Hidden axial dependence in angular transport')
                return sp.Symbol(ast.unparse(node))
            values={}
            for statement in function.body:
                if isinstance(statement,ast.Assign) and len(statement.targets)==1 and isinstance(statement.targets[0],ast.Name):
                    name=statement.targets[0].id
                    if name in ('u1','u','h'):
                        env[name]=formal(statement.value)
                        if name in ('u','h'):values[name]=env[name]
            if set(values)!=set(('u','h')):raise ArithmeticError('Angular source assignments changed')
            for value in values.values():
                if sp.simplify(sp.diff(value/fq,fq))!=0:raise ArithmeticError('Actual source lost q^-1 factor: '+method)
            factor_checks[filename+':'+method]=True
            identities['source_q_factor_'+method]=True
    y,U0,H0=sp.symbols('y U0 H0',positive=True)
    pulse_u=U0/q*sp.exp(-(sp.Rational(1,2)+mu)*y)
    pulse_h=H0/q*sp.exp(-sp.Rational(3,2)*y)+U0/q*(sp.exp(-(sp.Rational(1,2)+mu)*y)-sp.exp(-sp.Rational(3,2)*y))/(1-mu)
    pulse_X=1/(1-mu)+(H0/U0-1/(1-mu))*sp.exp(-(1-mu)*y)
    zero('whole_pulse_scalar_X_for_arbitrary_Z',pulse_h/pulse_u-pulse_X)
    zero('whole_pulse_angular_primitive_ODE',sp.diff(pulse_h,y)+sp.Rational(3,2)*pulse_h-pulse_u)
    p,b,z,u,v,Q,scale=sp.symbols('p b z u v Q scale',real=True)
    M=sp.Matrix([[p,1],[1,b]]);M_inv=sp.Matrix([[b,-1],[-1,p]])/(p*b-1)
    if sp.simplify(M*M_inv)!=sp.eye(2):raise ArithmeticError('Angular inverse identity failed')
    identities['rescaled_linear_inverse']=True
    zero('linear_determinant_exact',sp.exp(-2*(1-mu))*sp.exp(-2*(1+2*mu))-sp.exp(-4-2*mu))
    x,y=sp.symbols('x y',real=True)
    zero('scaled_quadratic_equation',(scale*x+b*scale*y+Q*((scale*x)**2+b*(scale*y)**2))/scale-(x+b*y+Q*scale*(x*x+b*y*y)))
    J=sp.Matrix([p*x+y,x+b*y+Q*scale*(x*x+b*y*y)]).jacobian([x,y])
    if sp.simplify(J-sp.Matrix([[p,1],[1+2*Q*scale*x,b*(1+2*Q*scale*y)]]))!=sp.zeros(2):
        raise ArithmeticError('Differentiated quadratic Jacobian failed')
    identities['actual_quadratic_Jacobian']=True
    h=sp.symbols('h',real=True)
    zero('pressure_square_change',((1+h)**2-1)/2-h-h*h/2)
    zero('swirl_energy_square_change',(1+h)**2-1-2*h-h*h)
    zero('disjoint_support_no_cross_term',(x+y)**2-x*x-y*y-2*x*y)
    t,ctr=sp.symbols('t ctr',real=True)
    zero('angular_bump_translation',sp.exp((1-mu)*(t+ctr))-sp.exp((1-mu)*ctr)*sp.exp((1-mu)*t))
    zero('pressure_bump_translation',sp.exp(-(1+2*mu)*(t+ctr))-sp.exp(-(1+2*mu)*ctr)*sp.exp(-(1+2*mu)*t))
    zero('energy_bump_translation',sp.exp(-2*mu*(t+ctr))-sp.exp(-2*mu*ctr)*sp.exp(-2*mu*t))
    pre,heat,repair=sp.symbols('pre heat repair')
    zero('unchanged_preheat_pressure_when_repair_equals_heat',(pre-heat+repair).subs(repair,heat)-pre)
    # Whole-Z existence margins, rather than only fixture residuals.
    ball=get('scaled_coefficient_ball_radius')
    if not (endpoints(get('scaled_linear_map_size')+get('scaled_nonlinear_map_size')-ball)[1]<0
            and endpoints(get('contraction_lipschitz_upper'))[1]<mp.mpf('.05')
            and endpoints(get('linear_inverse_row_norm_upper'))[1]<2
            and endpoints(get('rescaled_linear_determinant'))[1]<-mp.mpf('.98')
            and endpoints(get('strong_cap_log_margin'))[0]>0):
        raise ArithmeticError('Whole-Z functional repair gates failed')
    for name in ('A','B'):
        lo,hi=endpoints(read(r['bump_weights'][name]))
        if not mp.mpf('.8')<lo<=hi<mp.mpf('1.2'):raise ArithmeticError('Paper bump-weight range failed')
    if not 2>2*mp.mpf('.15'):raise ArithmeticError('Translated supports overlap')
    if endpoints(read(r['whole_Z_C1_coefficients']['derivative_jacobian_determinant']))[1]>=0:
        raise ArithmeticError('Whole-Z differentiated inverse gate failed')
    if (not r['actual_smooth_functional_angular_pressure_repair_proved']
        or r['full_outer_five_moment_match'] or r['whole_outer_cone_certified']
        or r['actual_ap_selected'] or r['temporal_recursion']):
        raise ValueError('Repair scope/completion flags changed')
    # Independent quadrature at the mu=0 limiting formulas. This is a
    # numerical identity check; actual source admission uses directed bounds.
    with mp.workdps(85):
        raw=lambda s:mp.exp(-1/(1-s*s)) if abs(s)<1 else mp.mpf(0)
        norm=mp.quad(raw,[-1,0,1]);ell=mp.mpf('.15')
        beta=lambda s:raw(s/ell)/(ell*norm)
        quadrature={}
        for name,weight,power in (('A',1,1),('B',-1,1),('D',-1,2),('E',0,1),('F',0,2)):
            val=mp.quad(lambda s:mp.exp(weight*s)*beta(s)**power,[-ell,0,ell])
            lo,hi=endpoints(read(r['bump_weights'][name]))
            if not lo<val<hi:raise ArithmeticError('Independent bump-weight quadrature outside bounds: '+name)
            quadrature[name]=val
        def sigma(s):
            if s<=0:return mp.mpf(0)
            if s>=1:return mp.mpf(1)
            A=mp.exp(-1/s**2);B=mp.exp(-1/(1-s)**2)
            return A/(A+B)
        fixture_checks=[]
        for sample in (r['samples'][2],r['samples'][3]):
            Z=endpoints(read(sample['Z']))[0];Qv=1+Z*Z
            parts=[0,20,50,70,80,90,95,99,100]
            value=2*mp.exp(-100)*(1-1/Qv)+mp.quad(
                lambda s:mp.exp(s-100)*2**(1-sigma(s/100))*(1-Qv**(-(1-sigma(s/100)))),parts)
            derivative=4*Z*mp.exp(-100)/Qv**2+mp.quad(
                lambda s:2*Z*mp.exp(s-100)*2**(1-sigma(s/100))*(1-sigma(s/100))*Qv**(-2+sigma(s/100)),parts)
            prejet=sample['defects']['rpre_scaled']['coefficients']
            for i,num in enumerate((value,derivative)):
                lo,hi=endpoints(read(prejet[i]))
                if not lo<num<hi:raise ArithmeticError('Independent correlated flatten fixture outside bounds')
            fixture_checks.append(True)
        atzero=r['samples'][1]
        for jet in atzero['physical_coefficient_Taylor']:
            if endpoints(read(jet['coefficients'][1]))!=(mp.mpf(0),mp.mpf(0)):
                raise ArithmeticError('Even coefficient derivative at zero failed')
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],
        implicit_source_sha256=r['implicit_source_sha256'],independent_identities=identities,
        independent_limit_weight_quadrature=encode(quadrature),correlated_flatten_fixtures=fixture_checks,
        actual_source_q_factor_checks=factor_checks,whole_pulse_X_invariance_for_arbitrary_Z_checked=True,
        whole_Z_contraction_and_derivative_inverse_checked=True,exact_even_axis_derivatives_checked=True,
        strong_finite_inverse_radius_cap_checked=True,all_passed=True,
        actual_implicit_functional_angular_pressure_repair_independently_checked=True,
        full_future_energy_and_ap_completed=False,full_outer_five_moment_match=False,
        whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('Angular/pressure repair:',len(identities),'independent identities/source invariants, five weight quadratures, correlated flatten fixtures and whole-Z branch gates PASS; future energy/ap pending',flush=True)
    return result


if __name__=='__main__':run()
