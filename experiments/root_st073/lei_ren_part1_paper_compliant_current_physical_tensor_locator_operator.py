"""Original radius-source geometry and implicit finite-point map identities."""
import ast
import copy
import functools
import inspect
import mpmath as mp
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_current_tensor_registry import ROUTES,ADJACENT
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_waiting_stress_C3 import source_precision

V=s.Symbol('native_coordinate',real=True)
NAMES=('log_epsilon','log_Ra','h_bridge','h_switch','log_Rref','log_Rp','log_Rv',
       'mu','T','log_gap','Md','Tw','Lrel','Ts','wait','log_C','log_P')
SYMBOLS={name:s.Symbol(name,real=True) for name in NAMES}
DOMAINS={row[0]:(s.Integer(0),s.Integer(1)) for row in ROUTES}
DOMAINS.update(core_positive_radius=(0,4),bridge_second=(1,2),switch_second=(1,2),
    restore_buffer=(-7,-6),actual_patch=(1,s.E),Rh_reference=(-5,0),O2_buffer=(0,11),
    pulse_entrance=(0,s.Rational(1,50)),pulse_main=(s.Rational(1,50),10),pulse_exit=(10,11),
    pulse_gap=(11,12),pulse_end=(-4,0),flatten=(0,100),outer_angular=(-4,0),
    heat_collar=(0,3),heat_exterior=(3,s.oo))

@functools.lru_cache(maxsize=1)
def exact_original_radius_program():
    """Keep the exact source expressions before numerical radius caps.

    Remove only the two rounded logR clamps and two singleton endpoint
    shortcuts. The original full source formulas then work symbolically.
    """
    asts=SourceAST();fn=copy.deepcopy(asts.method('global_physical_assembly','radius'))
    fn.name='_original_unclipped_radius_source';fn.decorator_list=[]
    removed=dict(numerical_logR_clamps=0,rounded_singleton_shortcuts=0)
    class SourceOnly(ast.NodeTransformer):
        def visit_Assign(self,node):
            if len(node.targets)==1 and ast.unparse(node.targets[0])=='logR' and any(
                isinstance(part,ast.Call) and ast.unparse(part.func)=='endpoints' for part in ast.walk(node.value)):
                removed['numerical_logR_clamps']+=1;return None
            return self.generic_visit(node)
        def visit_If(self,node):
            if any(isinstance(part,ast.Call) and ast.unparse(part.func)=='endpoints' for part in ast.walk(node.test)):
                removed['rounded_singleton_shortcuts']+=1;return None
            return self.generic_visit(node)
        def visit_Return(self,node):
            if not isinstance(node.value,ast.Tuple) or len(node.value.elts)!=2:raise ValueError('Original radius/source return changed')
            node.value=node.value.elts[0];return node
    fn=SourceOnly().visit(fn)
    if removed!=dict(numerical_logR_clamps=2,rounded_singleton_shortcuts=2):
        raise ValueError('Original exact radius/cap distinction changed')
    if any(isinstance(node,ast.Call) and ast.unparse(node.func)=='endpoints' for node in ast.walk(fn)):
        raise ValueError('Rounded source selection remains in exact radius geometry')
    env=dict(CompliantGlobalPhysicalAssembly.radius.__globals__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original-exact-radius-source-before-caps>','exec'),env)
    for stem,method in (('pre_pulse_mixed_C4','axial'),('current_pulse_gap_background_tensor','chart'),
        ('current_pulse_entrance_incoming_background_tensor','chart'),('current_core_stress_operator','core_raw_source')):
        asts.method(stem,method)
    return env[fn.name],dict(only_numerical_caps_and_rounded_shortcuts_removed=removed,
        all_original_exact_radius_expressions_and_chart_branching_retained=True,input_hashes=asts.hashes,passed=True)

@functools.lru_cache(maxsize=1)
def original_native_radius_recipes():
    program,proof=exact_original_radius_program();a=SYMBOLS
    c=SimpleNamespace(mpf=s.sympify,ln=s.log)
    params=SimpleNamespace(mu=a['mu'],Tw=a['Tw'])
    owner=SimpleNamespace(ctx=c,logRref=a['log_Rref'],logRp=a['log_Rp'],logP=a['log_P'],params=params)
    steep=SimpleNamespace(outer=SimpleNamespace(Lrel=a['Lrel']),Ts=a['Ts'],wait=a['wait'])
    provider=SimpleNamespace(bridge=SimpleNamespace(r=s.exp(a['log_Ra'])),
        reshape=SimpleNamespace(T=a['T']),reference=SimpleNamespace(reshape=SimpleNamespace(T=a['T']),loggap=a['log_gap']),
        Lrel=a['Lrel'],outer=steep.outer,Ts=a['Ts'],wait=a['wait'],heat=SimpleNamespace(steep=steep))
    result={}
    for region,*_ in ROUTES:
        if region=='core_positive_radius':value=a['log_epsilon']+s.log(V)
        else:
            h=a['h_bridge'] if region.startswith('bridge_') else a['h_switch']
            packet=dict(width_enclosure_is_not_source=h,source_radius_tree='same original formal radius tree',
                source_width_log='same exact positive h log',exact_positive_width_log='same exact positive h log',
                formal_radius_tree='same original formal radius tree')
            coordinate=V
            if region=='O2_axial':packet.update(actual_y=s.exp(a['Md']*V),exact_radius_source='original exp(Md*phase)')
            if region=='O2_buffer':packet.update(actual_y=s.exp(a['Md'])+V,exact_radius_source='original exp(Md)+buffer_offset')
            if region=='pulse_entrance':coordinate=V/a['mu']
            if region=='pulse_gap_end':coordinate=-(4*a['mu']+(1-4*a['mu'])*(1-V))/a['mu']
            value=program(owner,region,coordinate,packet,provider)
        result[region]=s.expand(value)
    if list(result)!=[row[0] for row in ROUTES]:raise ValueError('All original native radius recipes required')
    return result,proof

def evaluate_expression(ctx,expression,values):
    if expression in values:return values[expression]
    if expression is s.E:return ctx.exp(ctx.mpf(1))
    if expression is s.oo:return ctx.mpf('inf')
    if expression.is_Number:return ctx.mpf(str(expression.p))/int(expression.q) if expression.is_Rational else ctx.mpf(str(expression))
    if isinstance(expression,s.Add):return sum((evaluate_expression(ctx,arg,values) for arg in expression.args),ctx.mpf(0))
    if isinstance(expression,s.Mul):
        value=ctx.mpf(1)
        for arg in expression.args:value*=evaluate_expression(ctx,arg,values)
        return value
    if isinstance(expression,s.Pow):return evaluate_expression(ctx,expression.base,values)**evaluate_expression(ctx,expression.exp,values)
    if expression.func is s.log:return ctx.ln(evaluate_expression(ctx,expression.args[0],values))
    if expression.func is s.exp:return ctx.exp(evaluate_expression(ctx,expression.args[0],values))
    raise ValueError('Unsupported original radius source expression: '+str(expression))

@functools.lru_cache(maxsize=1)
def original_radial_cover_theorem():
    recipes,proof=original_native_radius_recipes();a=SYMBOLS
    substitutions={a['log_Ra']:s.log(4)+a['log_epsilon'],a['h_switch']:a['h_bridge'],
        a['log_Rref']:s.log(110)+10*(a['log_C']+a['log_P']),
        a['log_gap']:10*(a['log_C']+a['log_P'])-a['T']}
    substitutions[a['log_Rp']]=substitutions[a['log_Rref']]+a['log_P']+1+a['Tw']
    substitutions[a['log_Rv']]=substitutions[a['log_Rp']]+13/a['mu']
    regions=list(recipes);joins={}
    for j,(_,name,_,_) in enumerate(ADJACENT):
        left=recipes[regions[j]].subs(V,DOMAINS[regions[j]][1])
        right=recipes[regions[j+1]].subs(V,DOMAINS[regions[j+1]][0])
        residual=s.simplify((left-right).subs(substitutions).subs(a['log_P'],s.exp(a['Md'])+11))
        if residual!=0:raise ValueError('Original exact radial source boundary differs: '+name+' '+str(residual))
        joins[name]=dict(left_region=regions[j],right_region=regions[j+1],identity=True,
            left_native_endpoint=str(DOMAINS[regions[j]][1]),right_native_endpoint=str(DOMAINS[regions[j+1]][0]))
    slopes={region:str(s.diff(expression,V)) for region,expression in recipes.items()}
    return dict(original_exact_radius_source_program=proof,all_32_adjacent_exact_radius_source_identities=joins,
        original_native_log_radius_derivatives=slopes,
        source_parameter_equalities_required={str(k):str(v) for k,v in substitutions.items()},
        original_log_Pstar_source='exp(Md)+11',positive_source_requirements=[
            'epsilon>0; hb>0, same original hb at bridge and switches',
            'log100-logRa-2hb>0; log(110/100)-2hb>0',
            'T>0; log_gap-8>0; Md>0; Tw>0; 0<mu<1/4',
            'Lrel-4>0; Ts>0; wait>0'],
        radial_cover_endpoints='R=0 is the nonsingular core axis; heat_exterior offset grows to +infinity',
        formal_source_equalities_not_rounded_interval_overlap=True,input_hashes=proof['input_hashes'],passed=True)

@functools.lru_cache(maxsize=1)
def implicit_physical_source_theorem():
    asts=SourceAST()
    for stem,method in (('global_physical_assembly','evaluate'),('global_physical_assembly_check','implicit_physical_fixture'),
        ('collar_physical_C2','physical_source_row'),('collar_physical_C2','collar')):
        asts.method(stem,method)
    lam,nu,tau=s.symbols('lambda nu tau',positive=True);zeta,delta=s.symbols('Z delta',real=True)
    zphys=s.sqrt(nu)*zeta*lam**(1-delta)
    residual=s.simplify(s.expand_power_base(lam**2-lam**(2*delta)*zphys**2/nu-tau,force=True)
        -(lam**2*(1-zeta**2)-tau))
    if residual!=0:raise ValueError('Original viscosity/implicit source normalization differs')
    A,B=s.symbols('positive_tau positive_axial_term',positive=True)
    derivative=2-2*delta*B/(A+B)
    low=s.simplify(derivative-2*(1-delta));high=s.simplify(2-derivative)
    if low!=2*delta*A/(A+B) or high!=2*delta*B/(A+B):raise ValueError('Log lambda inverse monotonicity differs')
    rho,angle=s.symbols('rho angle',real=True)
    polar=s.trigsimp((s.sqrt(2*rho)*s.cos(angle))**2+(s.sqrt(2*rho)*s.sin(angle))**2-2*rho)
    if polar!=0:raise ValueError('Constrained Cartesian core polar identity differs')
    return dict(original_implicit_viscosity_coordinate_identity=True,
        defining_equation='lambda^2-(z_phys^2/nu)*lambda^(2delta)=tau',
        log_root_equation='F(q)=2q-logaddexp(log_tau,2log(abs(z_phys)/sqrt(nu))+2delta*q)=0',
        exact_global_derivative_bounds='2(1-delta)<=F_prime<=2 for 0<=delta<1',
        derivative_minus_lower=str(low),upper_minus_derivative=str(high),
        radial_source='logR=2log(r_phys)-log2-lognu-2q',
        exact_finite_point_relation='log(1-Z^2)=log_tau-2q; Z=(z_phys/sqrt(nu))*exp(-(1-delta)*q)',
        constrained_polar_Cartesian_core_rho_identity=True,
        actual_log_tau_is_not_replaced_by_2log_lambda=True,
        actual_lambda_factor_must_be_separate_from_native_sqrt_tau_upper=True,
        input_hashes=asts.hashes,passed=True)


def intersect(ctx,a,b):
    al,ah=endpoints(a);bl,bh=endpoints(b);lo=max(al,bl);hi=min(ah,bh)
    return None if lo>hi else ctx.mpf([lo,hi])


def nonpositive_exp(ctx,value):
    """Sound enclosure without constructing an astronomical MPF exponent.

    This is an enclosure of the original exp, never a replacement source.
    The exact logarithmic expression must remain attached to its caller.
    """
    lo,hi=endpoints(value)
    if hi>0:raise ValueError('Nonpositive exponential argument required')
    if hi<-1024:return ctx.mpf([0,endpoints(ctx.exp(ctx.mpf(-1024)))[1]])
    upper=endpoints(ctx.exp(ctx.mpf(hi)))[1]
    lower=mp.mpf(0) if lo<-1024 else endpoints(ctx.exp(ctx.mpf(lo)))[0]
    return ctx.mpf([lower,upper])


def logaddexp(ctx,a,b):
    shift=max(endpoints(a)[1],endpoints(b)[1])
    return ctx.mpf(shift)+ctx.ln(nonpositive_exp(ctx,a-shift)+nonpositive_exp(ctx,b-shift))


@source_precision
def implicit_log_coordinate_map(ctx,z_phys,log_tau,viscosity,delta,relative_tolerance='1e-40',max_steps=256):
    """Enclose all exact finite physical roots under the supplied parameter box.

    Midpoints are trial evaluation coordinates only. No parameter or defining
    field is replaced by a midpoint. Interval Newton uses the global positive
    derivative bound and never selects an arbitrary root from a parameter box.
    """
    z=ctx.mpf(z_phys);lt=ctx.mpf(log_tau);nu=ctx.mpf(viscosity);d=ctx.mpf(delta)
    if not all(mp.isfinite(v) for v in endpoints(z)+endpoints(lt)+endpoints(nu)+endpoints(d)):
        raise ValueError('Finite physical z, log tau, viscosity and delta required')
    dl,dh=endpoints(d)
    if endpoints(nu)[0]<=0 or dl<0 or dh>=1:raise ValueError('nu>0 and 0<=delta<1 required')
    tolerance=mp.mpf(relative_tolerance)
    if tolerance<=0 or not mp.isfinite(tolerance) or max_steps<1:raise ValueError('Positive finite locator tolerance/step count required')
    zl,zh=endpoints(z);absolute=abs(z)
    if zl==zh==0:
        q=lt/2;Z=ctx.mpf(0);log_complement=ctx.mpf(0);iterations=0;status='exact_axial_zero_source'
    else:
        abslo,abshi=endpoints(absolute)
        logabs=ctx.mpf([-mp.inf,endpoints(ctx.ln(ctx.mpf(abshi)))[1]]) if abslo==0 else ctx.ln(absolute)
        axial_log=2*logabs-ctx.ln(nu)
        lo=endpoints(lt/2)[0]
        # Both terms are at most half exp(2q) at this bound, for every
        # admissible delta; the interval division also handles negative q.
        hi=max(endpoints((lt+ctx.ln(2))/2)[1],endpoints((axial_log+ctx.ln(2))/(2*(1-d)))[1])
        if hi<lo:raise ValueError('Implicit source bracket inverted')
        derivative=ctx.mpf([endpoints(2*(1-d))[0],2])
        status='step_limit';iterations=0
        for iterations in range(1,max_steps+1):
            if hi-lo<=tolerance*(1+max(abs(lo),abs(hi))):status='requested_log_root_width';break
            mid=(lo+hi)/2
            f=2*ctx.mpf(mid)-logaddexp(ctx,lt,axial_log+2*d*ctx.mpf(mid))
            fl,fh=endpoints(f)
            if fh<0:lo=mid
            elif fl>0:hi=mid
            else:
                contracted=intersect(ctx,ctx.mpf([lo,hi]),ctx.mpf(mid)-f/derivative)
                if contracted is None:raise ValueError('Directed implicit source root enclosure lost')
                nl,nh=endpoints(contracted)
                if nh-nl>=(hi-lo)*mp.mpf('.999'):status='parameter_or_precision_limited';break
                lo,hi=nl,nh
        q=ctx.mpf([lo,hi]);log_complement=lt-2*q
        # The exact implicit relation proves these two log quantities <=0;
        # their interval upper ends can overshoot by outward rounding.
        log_complement=intersect(ctx,log_complement,ctx.mpf([-mp.inf,0]))
        log_Z=intersect(ctx,logabs-ctx.ln(nu)/2-(1-d)*q,ctx.mpf([-mp.inf,0]))
        if log_complement is None or log_Z is None:raise ValueError('Finite-point source relation inconsistent')
        magnitude=nonpositive_exp(ctx,log_Z)
        squared=intersect(ctx,1-nonpositive_exp(ctx,log_complement),ctx.mpf([0,1]))
        if squared is None:raise ValueError('Exact finite-point complementary Z source inconsistent')
        magnitude=intersect(ctx,magnitude,ctx.sqrt(squared))
        if magnitude is None:raise ValueError('Two exact physical Z source relations differ')
        Z=magnitude if zl>=0 else -magnitude if zh<=0 else ctx.mpf([-endpoints(magnitude)[1],endpoints(magnitude)[1]])
    return dict(physical_z=z,requested_log_tau=lt,physical_viscosity=nu,source_delta=d,
        actual_log_lambda=q,Z=Z,exact_log_one_minus_Z_squared_source=log_complement,
        exact_relation='lambda^2-(z_phys^2/nu)*lambda^(2delta)=tau; lambda^2*(1-Z^2)=tau',
        true_finite_point_has_abs_Z_strictly_below_one=True,
        Z_enclosure_can_touch_infinity_limit_due_to_precision=True,
        log_root_enclosure_width=ctx.mpf(endpoints(q)[1]-endpoints(q)[0]),
        solver_status=status,iterations=iterations,relative_log_root_tolerance=str(relative_tolerance),
        midpoint_used_only_as_root_iteration_trial=True,actual_log_tau_retained=True,
        exp_truncation_is_enclosure_only_not_source=True)


@functools.lru_cache(maxsize=1)
def constrained_core_program():
    """Original Cartesian source, with rho enclosed under its joint constraint.

    Interval X/Y boxes alone overestimate rho at the disk boundary. The only
    change is accepting an independently enclosed rho of the same exact point.
    Callers must bind X=sqrt(2rho)cos(theta), Y=sqrt(2rho)sin(theta).
    """
    from lei_ren_part1_paper_compliant_current_core_axis_background_tensor import CurrentCoreAxisBackgroundTensor
    asts=SourceAST();fn=copy.deepcopy(asts.method('current_core_axis_background_tensor','cartesian'))
    fn.name='_joint_polar_cartesian_source';fn.decorator_list=[]
    fn.args.kwonlyargs.append(ast.arg(arg='rho_source'));fn.args.kw_defaults.append(None)
    replacements=0
    for node in ast.walk(fn):
        if isinstance(node,ast.Assign) and len(node.targets)==1 and ast.unparse(node.targets[0])=='rho':
            if ast.unparse(node.value)!='(x ** 2 + y ** 2) / 2':raise ValueError('Original core radius assignment changed')
            node.value=ast.parse('c.mpf(rho_source)',mode='eval').body;replacements+=1
    if replacements!=1:raise ValueError('Exactly one joint radius enclosure assignment required')
    env=dict(inspect.unwrap(CurrentCoreAxisBackgroundTensor.cartesian).__globals__)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original-Cartesian-core-with-exact-joint-polar-constraint>','exec'),env)
    return env[fn.name],dict(original_full_Cartesian_source_AST_retained=True,
        only_joint_radius_enclosure_assignment_changed=True,
        required_exact_joint_constraint='X=sqrt(2rho)cos(theta); Y=sqrt(2rho)sin(theta); rho=(X^2+Y^2)/2',
        pointwise_identity=implicit_physical_source_theorem()['constrained_polar_Cartesian_core_rho_identity'],
        original_full_mixed_derivative_operators_and_six_source_sectors_retained=True,
        input_hashes=asts.hashes,passed=True)
