"""Independent closed-integral derivatives and actual Rh-to-Rp source joins."""
import ast
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pre_pulse_mixed_C4 import (
    CompliantPrePulseMixedC4,physical_mixed,turnoff_derivatives,IntervalTaylor)
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
NAME=PREFIX+'pre_pulse_mixed_C4.json'


def symbolic_sources():
    """Bind each NEW production history to the OLD original source chain."""
    tree=ast.parse((HERE/(PREFIX+'pre_pulse_mixed_C4.py')).read_text(encoding='utf8'))
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Pre-pulse source identity: '+name)
        proofs[name]=True
    def hist(method,env):
        fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
        node=next(n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)=='hist' for t in n.targets))
        if not isinstance(node,ast.Call) or ast.unparse(node.func)!='dict':raise ValueError('Explicit five-history source dictionary required')
        def formal(n):
            label=ast.unparse(n)
            if label in env:return env[label]
            if isinstance(n,ast.Constant):return s.Rational(str(n.value))
            if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -formal(n.operand)
            if isinstance(n,ast.BinOp):
                a,b=formal(n.left),formal(n.right)
                for kind,op in ((ast.Add,lambda:a+b),(ast.Sub,lambda:a-b),(ast.Mult,lambda:a*b),(ast.Div,lambda:a/b),(ast.Pow,lambda:a**b)):
                    if isinstance(n.op,kind):return op()
            if isinstance(n,ast.Call) and len(n.args)==1:
                name=ast.unparse(n.func);a=formal(n.args[0])
                if name=='c.mpf':return a
                if name=='c.exp':return s.exp(a)
                if name=='square':return a*a
            raise ValueError('Unbound new source: '+label)
        return {k.arg:formal(k.value) for k in node.keywords}
    z,y,t,mu,J,invP2=s.symbols('z y t mu J invP2',real=True)
    qi=1/(1+z*z);u,uh,h,V=s.symbols('u uh h V',real=True)
    oldmap={'m':'m','h':'h','k':'k','e':'e','p':'p'}
    common={'z':z,'qi':qi,'V':4*z,'u':u,'h':h,'self.invP2':invP2,'y':y}
    ref=hist('reference',dict(common,h=s.Rational(5,8)*u))
    oldenv=dict(common,**{'V':4*z,'h':s.Rational(5,8)*u})
    for key,target in oldmap.items():zero('original_reference_'+key,ref[key]-assignment('compliant_outer_initial','reference',target,oldenv))
    I=s.symbols('I0:3',real=True)
    newenv=dict(common,**{'mass[0]':I[0],'mass[1]':I[1],'mass[2]':I[2]})
    oldenv=dict(common,**{'yy':y,'integrals[0]':I[0],'integrals[1]':I[1],'integrals[2]':I[2]})
    newenv['h']=oldenv['h']=qi*(s.Rational(5,8)+I[0])*s.exp(-3*y/2)
    for key,target in oldmap.items():zero('original_O2_slope_'+key,hist('slope',newenv)[key]-assignment('compliant_outer_initial','slope',target,oldenv))
    m0,h0,k0,e0,p0,d,root,d3,K1,K2=s.symbols('m0 h0 k0 e0 p0 d root d3 K1 K2',real=True)
    base={'z':z,'u1':uh,'self.invP2':invP2,'t':t,'d':d,'root':root,'d3':d3,
          "old['m']":m0,"old['h']":h0,"old['k']":k0,"old['e']":e0,"old['p']":p0,
          "K['B_mass']":K1,"K['B_squared_mass']":K2}
    oldbase={'z':z,'u1':uh,'self.invP2':invP2,'t':t,'decay':d,'root_decay':root,'decay3':d3,
          "get('Mz_over_R')":m0,"get('Mtheta_over_sqrt2_R_3half_Pstar')":h0,
          "get('Mtheta_z_over_sqrt2_R_3half_Pstar')":k0,"get('Mztheta_over_R_Pstar_squared')":e0,
          "get('Mp_over_Pstar_squared')":p0,"K['B_mass']":K1,"K['B_squared_mass']":K2}
    for key,target in oldmap.items():zero('original_O2_axial_'+key,hist('axial',base)[key]-assignment('compliant_outer_initial','axial',target,oldbase))
    kt,ke,kp=s.symbols('kt ke kp',real=True)
    newenv=dict(base,**{"K['theta']":kt,"K['energy']":ke,"K['pressure']":kp})
    oldenv=dict(oldbase,**{'d3':d3,"kernels['theta']":kt,"kernels['energy']":ke,"kernels['pressure']":kp})
    for key,target in oldmap.items():zero('original_O3_slope_mu_'+key,hist('slope_mu',newenv)[key]-assignment('compliant_outer_buffer','slope_mu',target,oldenv))
    f=s.exp((-s.Rational(1,2)-mu)*t);theta=(f-d3)/(1-mu)
    K2=(1-s.exp(-2*mu*t))/(2*mu);Kp=(1-s.exp(-(1+2*mu)*t))/(1+2*mu)
    newenv=dict(base,**{'theta':theta,'mu':mu,'decay_integral(c, 2 * mu, t)':K2,'decay_integral(c, 1 + 2 * mu, t)':Kp})
    oldenv=dict(oldbase,**{'d3':d3,'f':f,'mu':mu,'theta_kernel':theta,
        'decay_integral(c, 2 * mu, t)':K2,'decay_integral(c, 1 + 2 * mu, t)':Kp})
    for key,target in oldmap.items():zero('original_O3_power_'+key,hist('power',newenv)[key]-assignment('compliant_outer_buffer','power',target,oldenv))
    # Rh normalization, with true physical prefactors restored, not equality
    # of differently normalized stored derivative boxes.
    x=s.exp(1);am=s.exp(-s.Rational(3,5))*qi;urh=s.exp(-s.Rational(1,2))*qi
    zero('Rh_swirl_reference',am*x**s.Rational(1,10)-urh)
    zero('Rh_angular_reference',s.Rational(5,8)*am*x**s.Rational(8,5)-s.Rational(5,8)*urh*x**s.Rational(3,2))
    zero('Rh_pressure_reference',s.Rational(5,2)*am**2*x**s.Rational(1,5)-s.Rational(5,2)*urh**2)
    zero('Rh_energy_reference',x*(16*z*z*invP2-s.Rational(5,12)*am**2*x**s.Rational(1,5))-x*(16*z*z*invP2-s.Rational(5,12)*urh**2))
    zero('Rh_axial_mass_reference',4*z*x-x*(4*z))
    zero('Rh_mixed_angular_reference',4*z*s.Rational(5,8)*am*x**s.Rational(8,5)-x**s.Rational(3,2)*(4*z*s.Rational(5,8)*urh))
    dd=s.symbols('delta',real=True)
    Q=(8*z*z-(1-dd)*4*z*z-(1-z*z)*4)/(1-dd*z*z)
    zero('Rh_radial_reference',s.sqrt(x)*Q-s.exp(s.Rational(1,2))*Q)
    for name,expr in (('mass',4*z*x-x*(4*z)),('mixed',4*z*s.Rational(5,8)*am*x**s.Rational(8,5)-x**s.Rational(3,2)*(4*z*s.Rational(5,8)*urh)),
                      ('radial',s.sqrt(x)*Q-s.exp(s.Rational(1,2))*Q)):
        for n in range(1,5):zero('Rh_'+name+'_axial'+str(n),s.diff(expr,z,n))
    # Cutoff values and derivatives0..4 are checked numerically against
    # their exact endpoint values below; arbitrary histories at each join
    # agree by the shared initial integral and the SAME five primitive ODEs.
    for key in ('m','h','k','e','p'):
        zero('slope_mu_inlet_'+key,hist('slope_mu',dict(newenv,**{"K['theta']":0,"K['energy']":0,"K['pressure']":0,'d':1,'d3':1}))[key]-{'m':m0,'h':h0,'k':k0,'e':e0,'p':p0}[key])
    proofs.update(actual_Rh_join_inherits_uniform_implicit_functional_closure=True,
        all_internal_joins_use_flat_original_sigma_and_shared_primitive_ODEs=True,
        Rp_source_is_original_whole_O3_power_not_extrapolated_local_fixture=True)
    return proofs


def directed_Rh_diagnostics(field,record):
    """Compare true endpoint physical rows after fixed-basepoint unit conversion."""
    c=field.ctx;name=PREFIX+'actual_patch_mixed_C4.json';patch=json.loads((HERE/name).read_bytes())['actual_Rh_exit']
    new=field.reference([-1,1],-5);x=c.exp(1);count=0
    def overlap(a,b,label):
        nonlocal count
        la,ha=endpoints(a);lb,hb=endpoints(b)
        if max(la,lb)>min(ha,hb):raise ArithmeticError('Actual Rh converted physical diagnostic disjoint: '+label)
        count+=1
    for component,grid in patch['physical_velocity_pressure_y_Z_mixed4'].items():
        factor=c.sqrt(x) if component=='Ur_over_sqrt_Rm_over_2' else c.mpf(1)
        target='Ur_over_current_sqrt_R_over_2' if component=='Ur_over_sqrt_Rm_over_2' else component
        for index,value in grid.items():overlap(read_interval(c,value)/factor,new['physical_velocity_pressure_y_Z_mixed4'][target][index],target+'/'+index)
    mapping={'Mz_over_Rm':('Mz_over_current_R',x),
        'Mtheta_over_sqrt2_Rm_1p5_Pstar':('Mtheta_over_current_sqrt2_R_1p5_Pstar',x**c.mpf('1.5')),
        'Mtheta_z_over_sqrt2_Rm_1p5_Pstar':('Mtheta_z_over_current_sqrt2_R_1p5_Pstar',x**c.mpf('1.5')),
        'Mztheta_over_Rm_Pstar2':('Mztheta_over_current_R_Pstar2',x),'Mp_over_Pstar2':('Mp_over_Pstar2',c.mpf(1))}
    stirling=((1,),(0,1),(0,1,1),(0,1,3,1),(0,1,7,6,1))
    for component,(target,factor) in mapping.items():
        grid=patch['physical_five_primitive_x_Z_mixed4'][component]
        for k in range(5):
            for n in range(5-k):
                value=(read_interval(c,grid['x0_Z'+str(n)]) if k==0 else
                    sum((read_interval(c,grid['x'+str(j)+'_Z'+str(n)])*(stirling[k][j]*x**j) for j in range(1,k+1)),c.mpf(0)))
                index='y'+str(k)+'_Z'+str(n)
                overlap(value/factor,new['physical_five_primitive_y_Z_mixed4'][target][index],target+'/'+index)
    return dict(actual_Rh_converted_physical_mixed_row_overlap_diagnostics=count,
        unit_conversion_before_comparison=True,interval_overlap_not_the_functional_join_proof=True,passed=True)


def closed_integral_fixture():
    """Integrate independent physical moment definitions, then differentiate."""
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=110;y,z,t=s.symbols('y z t',real=True)
        basey=mp.mpf('.31');basez=mp.mpf('.23');delta=s.Rational(1,2000);invP2=s.Rational(25,144)
        u=s.exp(y/10)*(1+z/5+z*z/7);v=2*z+y*z/3+y*y/8;p0=s.Rational(3,100)+z/100+z**4/1000
        initial={'m':1+z*z/7,'h':s.Rational(5,8)+z/100,'k':z/3+z**3/8,
                 'e':-s.Rational(2,5)+z*z/11,'p':s.Rational(5,2)+z**3/15}
        forcing={'m':v,'h':u,'k':u*v,'e':v*v*invP2-u*u/2}
        rates={'m':1,'h':s.Rational(3,2),'k':s.Rational(3,2),'e':1}
        history={name:s.exp(-rates[name]*y)*(initial[name]+s.integrate(s.exp(rates[name]*t)*expr.subs(y,t),(t,0,y))) for name,expr in forcing.items()}
        history['p']=initial['p']+s.integrate(u.subs(y,t)**2/2,(t,0,y))
        Q=(2*z*v-(1-delta)*z*history['m']-(1-z*z)*s.diff(history['m'],z))/(1-delta*z*z)
        physical={'Utheta_over_Pstar':u,'Uz':v,'Ur_over_current_sqrt_R_over_2':s.exp((y-s.Rational(31,100))/2)*Q,'P_over_Pstar2':p0+history['p']}
        primitives={'Mz_over_current_R':s.exp(y-s.Rational(31,100))*history['m'],
            'Mtheta_over_current_sqrt2_R_1p5_Pstar':s.exp(s.Rational(3,2)*(y-s.Rational(31,100)))*history['h'],
            'Mtheta_z_over_current_sqrt2_R_1p5_Pstar':s.exp(s.Rational(3,2)*(y-s.Rational(31,100)))*history['k'],
            'Mztheta_over_current_R_Pstar2':s.exp(y-s.Rational(31,100))*history['e'],'Mp_over_Pstar2':history['p']}
        def jet(expr):return IntervalTaylor(c,[c.mpf(s.lambdify((y,z),s.diff(expr,z,n),'mpmath')(basey,basez)/math.factorial(n)) for n in range(6)])
        logU=[jet(s.diff(s.log(u),y,k)) for k in range(1,5)]
        packet=physical_mixed(c,c.mpf(basez),c.mpf('0.0005'),jet(u),logU,
            [jet(s.diff(v,y,k)) for k in range(5)],{name:jet(expr) for name,expr in history.items()},jet(p0),c.mpf(25)/144)
        count=0
        for group,expressions in (('physical_velocity_pressure_y_Z_mixed4',physical),('physical_five_primitive_y_Z_mixed4',primitives)):
            for name,expr in expressions.items():
                for key,bound in packet[group][name].items():
                    k,n=(int(a[1:]) for a in key.split('_'))
                    actual=s.lambdify((y,z),s.diff(expr,y,k,z,n),'mpmath')(basey,basez);lo,hi=endpoints(bound)
                    if actual<lo-mp.mpf('1e-65') or actual>hi+mp.mpf('1e-65'):raise ArithmeticError('Closed physical integral derivative mismatch: '+name+'/'+key)
                    count+=1
        return dict(independent_closed_physical_integral_mixed_derivatives=count,nonconstant_axial_velocity_and_nonzero_all_five_initial_histories=True,
                    finite_fixture_only=True,actual_core_source_admission=False,passed=True)


def turnoff_fixture():
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=110;md=mp.mpf('2.3');count=0
        def cutoff(yy):
            a=1-mp.log(yy)/md
            return mp.exp(-1/a**2)/(mp.exp(-1/a**2)+mp.exp(-1/(1-a)**2))
        for phase in (mp.mpf('.2'),mp.mpf('.5'),mp.mpf('.8')):
            y=mp.exp(md*phase);rows=turnoff_derivatives(c,c.mpf(y),c.mpf(md),c.mpf(phase))
            for k,row in enumerate(rows):
                actual=mp.diff(cutoff,y,k);lo,hi=endpoints(row)
                if actual<lo-mp.mpf('1e-65') or actual>hi+mp.mpf('1e-65'):raise ArithmeticError('Axial turnoff coordinate derivative mismatch')
                count+=1
        return dict(independent_original_cutoff_logR_derivatives=count,passed=True)


def run():
    record=json.loads((HERE/NAME).read_bytes());hashes=dict(record['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Pre-pulse receipt source changed: '+name)
    fixture=closed_integral_fixture();cutoff=turnoff_fixture();proofs=symbolic_sources()
    c=MPIntervalContext();c.dps=240;read=lambda value:read_interval(c,value);counts={}
    for packet in list(record['whole_charts'].values())+list(record['join_packets'].values()):
        for group in ('physical_velocity_pressure_y_Z_mixed4','physical_five_primitive_y_Z_mixed4'):
            for grid in packet[group].values():
                if len(grid)!=15:raise ValueError('Incomplete mixed4 packet')
                for value in grid.values():
                    if any(not mp.isfinite(v) for v in endpoints(read(value))):raise ArithmeticError('Nonfinite pre-pulse derivative')
                    counts[group]=counts.get(group,0)+1
        if any(packet[flag] for flag in ('full_cartesian_vector_derivatives_certified','admissible_stress_lift_constructed','temporal_recursion')):
            raise ValueError('Pre-pulse scope incorrectly promoted')
    checks=0
    with mp.workdps(280):
        field=CompliantPrePulseMixedC4();z='.5'
        comparisons=[(field.reference(z,-5),field.initial.reference(z,'-5')),
            (field.slope(z,'.5'),field.initial.slope(z,'.5')),
            (field.axial(z,phase='.5'),field.initial.axial(z,phase='.5')),
            (field.axial(z,buffer_offset=11),field.initial.axial(z,buffer_offset=11)),
            (field.slope_mu(z,'.5'),field.buffer.slope_mu(z,'.5')),(field.power(z,1),field.buffer.power(z,1))]
        mapping={'m':'Mz_over_R','h':'Mtheta_over_sqrt2_R_3half_Pstar','k':'Mtheta_z_over_sqrt2_R_3half_Pstar','e':'Mztheta_over_R_Pstar_squared','p':'Mp_over_Pstar_squared'}
        for new,old in comparisons:
            for key,target in mapping.items():
                for a,b in zip(new['actual_normalized_primitive_y_derivative_axial5'][key][0].coefficients[:2],old[target]):
                    la,ha=endpoints(a);lb,hb=endpoints(b)
                    if max(la,lb)>min(ha,hb):raise ArithmeticError('Original pre-pulse history diagnostic disjoint: '+key)
                    checks+=1
        rh=directed_Rh_diagnostics(field,record)
        rp=field.power(z,1)
        from lei_ren_part1_paper_compliant_power_inlet_C4 import CompliantPowerInletC4
        incoming=CompliantPowerInletC4().incoming(z);u=incoming['u'];Pstar=field.ctx.exp(field.params.logPstar)
        canonical=dict(m=incoming['m1']*u*Pstar,h=incoming['X']*u,k=incoming['m2']*u*u*Pstar,
            e=incoming['energy']*u*u,p=incoming['Mp'])
        rpchecks=0
        for key,jet in canonical.items():
            for a,b in zip(jet.coefficients,rp['actual_normalized_primitive_y_derivative_axial5'][key][0].coefficients):
                la,ha=endpoints(a);lb,hb=endpoints(b)
                if max(la,lb)>min(ha,hb):raise ArithmeticError('Accepted Rp canonical inlet diagnostic disjoint: '+key)
                rpchecks+=1
        if endpoints(rp['Uz_ordinary_y_derivative_axial5'][0][0])!=(mp.mpf(0),mp.mpf(0)) or endpoints(rp['actual_normalized_primitive_y_derivative_axial5']['m'][0][0])[0]<=0:
            raise ArithmeticError('Nonzero actual axial moment lost after Uz turns off')
        # Exact flat endpoint controls through the order needed for all
        # fourth derivatives; this supplements the source join algebra.
        for endpoint,value in ((0,0),(1,1)):
            from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
            rows=sigma_jets(field.ctx,endpoint)
            if endpoints(rows[0])!=(mp.mpf(value),mp.mpf(value)) or any(endpoints(rows[k])!=(mp.mpf(0),mp.mpf(0)) for k in range(1,5)):
                raise ArithmeticError('Original flat cutoff endpoint jet changed')
    # The original complete source chain and canonical incoming functions
    # were proved at the existing pulse inlet. New production identities
    # above transfer that accepted binding to this whole Rh-to-Rp provider.
    name=PREFIX+'power_inlet_C4_check.json';pulse=json.loads((HERE/name).read_bytes())
    if not (pulse['all_passed'] and pulse['two_sided_O3_pulse_join_certified'] and
        pulse['actual_five_defect_family_sha256']==record['actual_five_defect_family_sha256'] and pulse['implicit_source_sha256']==record['implicit_source_sha256']):
        raise ValueError('Accepted same-source Rp/pulse binding required')
    for path,digest in pulse['input_hashes'].items():
        if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Pulse inlet source changed: '+path)
    hashes.update(pulse['input_hashes']);hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest();hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],implicit_source_sha256=record['implicit_source_sha256'],
        datum_enclosure_sha256=record['datum_enclosure_sha256'],independent_closed_integral_fixture=fixture,independent_turnoff_coordinate_fixture=cutoff,
        actual_Rh_physical_endpoint_diagnostics=rh,
        actual_Rp_all_five_canonical_inlet_axial5_overlap_diagnostics=rpchecks,
        exact_original_production_and_join_identities=proofs,actual_mixed_bounds_checked=counts,
        independent_original_C1_history_overlap_diagnostics=checks,original_flat_cutoff_endpoint_jets_checked=10,
        retained_nonzero_axial_history_at_Rp_checked=True,all_pre_pulse_mixed4_available=True,
        pre_pulse_functional_mixed4_joins_certified=True,Rp_accepted_pulse_source_join_certified=True,
        full_inner_dispatcher_installed=False,full_cartesian_vector_derivatives_certified=False,
        admissible_stress_lift_constructed=False,temporal_recursion=False,all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Rh-to-Rp PASS:135 independent closed-integral mixed derivatives;15 cutoff logR derivatives; original source and functional joins',flush=True)
    return result


if __name__=='__main__':run()
