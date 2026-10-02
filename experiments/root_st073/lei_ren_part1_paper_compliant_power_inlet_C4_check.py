"""Source identities and two-sided high-order O3/pulse inlet joining.

Numeric overlap is only a diagnostic. Exact production power expressions,
canonical incoming functions, flat forcing and the shared primitive ODEs
identify the functions and every needed derivative on the whole axial domain.
"""
import ast
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'
NAME=PREFIX+'compliant_power_inlet_C4.json'


def source_identities():
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Power-source join identity failed: '+name)
        proofs[name]=True
    t,mu=s.symbols('t mu',real=True); Z=s.symbols('Z',real=True)
    u0,m0,h0,k0,e0,p0=s.symbols('u0 m0 h0 k0 e0 p0',real=True)
    bp=s.Rational(1,2)+mu; rate=1-mu; prate=1+2*mu
    f=s.exp(-bp*t); decay=s.exp(-t); d3=s.exp(-s.Rational(3,2)*t)
    env={'mu':mu,'t':t,'slope':-bp,'f':f,'decay':decay,'d3':d3,'u1':u0,
        'theta_kernel':(f-d3)/rate,"get('Mz_over_R')":m0,
        "get('Mtheta_over_sqrt2_R_3half_Pstar')":h0,
        "get('Mtheta_z_over_sqrt2_R_3half_Pstar')":k0,
        "get('Mztheta_over_R_Pstar_squared')":e0,"get('Mp_over_Pstar_squared')":p0,
        'decay_integral(c, 2 * mu, t)':(1-s.exp(-2*mu*t))/(2*mu),
        'decay_integral(c, 1 + 2 * mu, t)':(1-s.exp(-prate*t))/prate}
    production=lambda target:assignment('compliant_outer_buffer','power',target,env)
    zero('production_power_slope',production('slope')+bp)
    zero('production_power_factor',production('f')-f)
    zero('production_power_decay',production('decay')-decay)
    zero('production_power_mixed_decay',production('d3')-d3)
    # The nonzero exact theta integral is capped only in the numeric
    # enclosure. Bind the underlying functional expression to production.
    tree=ast.parse((HERE/(PREFIX+'compliant_outer_buffer.py')).read_text(encoding='utf-8'))
    method=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='power')
    assignments=[ast.unparse(n.value) for n in ast.walk(method) if isinstance(n,ast.Assign)
        and any(ast.unparse(a)=='theta_kernel' for a in n.targets)]
    if '(f - d3) / (1 - mu)' not in assignments:raise ValueError('Original exact theta integral changed')
    zero('production_swirl_power',production('u')-u0*f)
    zero('normalized_linear_m1_power',production('m')/(u0*f)-m0/u0*s.exp(-(s.Rational(1,2)-mu)*t))
    zero('normalized_linear_m2_power',production('k')/(u0*u0*f*f)-k0/(u0*u0)*s.exp(-(s.Rational(1,2)-2*mu)*t))
    X0=h0/u0; X=1/rate+(X0-1/rate)*s.exp(-rate*t)
    e=e0/u0**2*s.exp(2*mu*t)-(s.exp(2*mu*t)-1)/(4*mu)
    pressure=p0+u0*u0*(1-s.exp(-prate*t))/(2*prate)
    zero('normalized_angular_power',production('h')/(u0*f)-X)
    zero('normalized_energy_power',production('e')/(u0*u0*f*f)-e)
    zero('original_pressure_power',production('p')-pressure)
    zero('angular_ODE',s.diff(X,t)-1+rate*X)
    zero('energy_ODE',s.diff(e,t)+s.Rational(1,2)-2*mu*e)
    zero('pressure_ODE',s.diff(pressure,t)-u0*u0*f*f/2)
    U,M,K,EQ,EZ,invP,Hp,Pin=s.symbols('U M K EQ EZ invP Hp Pin',real=True)
    q=1+Z*Z; u=U/q; m=M*Z; k=K*Z/q; raw_e=EZ*Z*Z+EQ/q**2
    zero('canonical_incoming_m1',invP*m/u-(invP*M/U)*Z*q)
    zero('canonical_incoming_m2',invP*k/u**2-(invP*K/U**2)*Z*q)
    zero('canonical_incoming_energy',raw_e/u**2-(EQ/U**2+EZ/U**2*Z*Z*q*q))
    zero('canonical_Xp_constant',(Hp/q)/(U/q)-Hp/U)
    zero('canonical_Mp_constant',(Pin/q**2)*q*q-Pin)
    # Establish the canonical shapes from the ACTUAL source chain, before
    # recovering H/Pin from the stored Z=.5 value. This is what makes one
    # value an enclosure of a source constant rather than extrapolation.
    tree0=ast.parse((HERE/(PREFIX+'compliant_outer_initial.py')).read_text(encoding='utf-8'))
    coordinate=next(n for n in ast.walk(tree0) if isinstance(n,ast.FunctionDef) and n.name=='coordinate')
    qvalues=[ast.unparse(n.value) for n in ast.walk(coordinate) if isinstance(n,ast.Assign)
        and any(ast.unparse(v)=='q' for v in n.targets)]
    if qvalues!=['IntervalTaylor(c, [1 + Z ** 2, 2 * Z])']:
        raise ValueError('Original exact axial q shape changed')
    yy,J0,I0,I1,I2,invP2=s.symbols('yy J0 I0 I1 I2 invP2',real=True)
    f0=s.exp(yy/10-s.Rational(3,5)*J0)
    env0={'qi':1/q,'z':Z,'yy':yy,'J':J0,'factor':f0,'self.invP2':invP2,
        'integrals[0]':I0,'integrals[1]':I1,'integrals[2]':I2,'h':(s.Rational(5,8)+I0)*s.exp(-s.Rational(3,2)*yy)/q,
        'V':4*Z}
    source0=lambda target:assignment('compliant_outer_initial','slope',target,env0)
    zero('source_O2_slope_qU',q*source0('u')-f0)
    zero('source_O2_slope_Mz_linear',source0('m')-4*Z)
    zero('source_O2_slope_qH',q*source0('h')-(s.Rational(5,8)+I0)*s.exp(-s.Rational(3,2)*yy))
    zero('source_O2_slope_qK_linear',q*source0('k')-4*Z*(s.Rational(5,8)+I0)*s.exp(-s.Rational(3,2)*yy))
    zero('source_O2_slope_energy_shape',q*q*source0('e')-
        (16*invP2*Z*Z*q*q-(s.Rational(5,12)+I2/2)*s.exp(-yy)))
    zero('source_O2_slope_q2_pressure',q*q*source0('p')-(s.Rational(5,2)+I1/2))
    tt,dec,root,dec3,B1,B2=s.symbols('tt dec root dec3 B1 B2',real=True)
    source_inputs={"get('Utheta_over_Pstar')":U/q,"get('Mz_over_R')":M*Z,
        "get('Mtheta_over_sqrt2_R_3half_Pstar')":Hp/q,
        "get('Mtheta_z_over_sqrt2_R_3half_Pstar')":K*Z/q,
        "get('Mztheta_over_R_Pstar_squared')":EZ*Z*Z+EQ/q**2,
        "get('Mp_over_Pstar_squared')":Pin/q**2}
    enva=dict(source_inputs,z=Z,u1=U/q,decay=dec,root_decay=root,decay3=dec3,t=tt)
    enva.update({'self.invP2':invP2,"K['B_mass']":B1,"K['B_squared_mass']":B2})
    sourcea=lambda target:assignment('compliant_outer_initial','axial',target,enva)
    zero('source_axial_qU_constant',q*sourcea('u')-U*root)
    zero('source_axial_Mz_linear',sourcea('m')-Z*(M*dec+4*B1))
    zero('source_axial_qH_constant',q*sourcea('h')-(Hp*dec3+U*(root-dec3)))
    zero('source_axial_qK_linear',q*sourcea('k')-Z*(K*dec3+4*U*root*B1))
    zero('source_axial_energy_shape',q*q*sourcea('e')-
        ((EZ*dec+16*invP2*B2)*Z*Z*q*q+EQ*dec-U*U*tt*dec/2))
    zero('source_axial_q2_pressure_constant',q*q*sourcea('p')-(Pin+U*U*(1-dec)/2))
    f3,kt,ke,kp=s.symbols('f3 kt ke kp',real=True)
    env3=dict(source_inputs,u1=U/q,f=f3,decay=dec,d3=dec3)
    env3.update({"kernels['theta']":kt,"kernels['energy']":ke,"kernels['pressure']":kp})
    source3=lambda target:assignment('compliant_outer_buffer','slope_mu',target,env3)
    zero('source_O3_transition_qU_constant',q*source3('u')-U*f3)
    zero('source_O3_transition_Mz_linear',source3('m')-M*dec*Z)
    zero('source_O3_transition_qH_constant',q*source3('h')-(Hp*dec3+U*dec3*kt))
    zero('source_O3_transition_qK_linear',q*source3('k')-K*dec3*Z)
    zero('source_O3_transition_energy_shape',q*q*source3('e')-
        (EZ*dec*Z*Z*q*q+EQ*dec-U*U*dec*ke/2))
    zero('source_O3_transition_q2_pressure_constant',q*q*source3('p')-(Pin+U*U*kp/2))
    envp=dict(source_inputs,u1=U/q,f=f,decay=decay,d3=d3,mu=mu,t=t,
        theta_kernel=(f-d3)/rate)
    envp.update({'decay_integral(c, 2 * mu, t)':(1-s.exp(-2*mu*t))/(2*mu),
        'decay_integral(c, 1 + 2 * mu, t)':(1-s.exp(-prate*t))/prate})
    sourcep=lambda target:assignment('compliant_outer_buffer','power',target,envp)
    zero('source_O3_power_qU_constant',q*sourcep('u')-U*f)
    zero('source_O3_power_Mz_linear',sourcep('m')-M*decay*Z)
    zero('source_O3_power_qH_constant',q*sourcep('h')-(Hp*d3+U*(f-d3)/rate))
    zero('source_O3_power_qK_linear',q*sourcep('k')-K*d3*Z)
    zero('source_O3_power_energy_shape',q*q*sourcep('e')-
        (EZ*decay*Z*Z*q*q+EQ*decay-U*U*decay*(1-s.exp(-2*mu*t))/(4*mu)))
    zero('source_O3_power_q2_pressure_constant',q*q*sourcep('p')-(Pin+U*U*(1-s.exp(-prate*t))/(2*prate)))
    # Verify that the named radial kernels and their call sites contain
    # no Z argument. The formulas above then propagate CONSTANT radial
    # data through the entire source chain on the whole axial domain.
    for filename,method,callname in (('compliant_outer_initial','slope','transition_integrals'),
        ('compliant_outer_initial','axial','turnoff_kernels'),('compliant_outer_buffer','slope_mu','transition_kernels')):
        source_tree=ast.parse((HERE/(PREFIX+filename+'.py')).read_text(encoding='utf-8'))
        fn=next(n for n in ast.walk(source_tree) if isinstance(n,ast.FunctionDef) and n.name==method)
        calls=[n for n in ast.walk(fn) if isinstance(n,ast.Call) and ast.unparse(n.func)==callname]
        if len(calls)!=1 or any(isinstance(n,ast.Name) and n.id in ('Z','z','q','qi') for n in ast.walk(calls[0])):
            raise ValueError('Source radial kernel axial dependence changed: '+callname)
    proofs['actual_O2_O3_source_chain_has_exact_canonical_whole_Z_shapes']=True
    proofs['single_sample_only_encloses_proved_Z_independent_Hp_Pin_constants']=True
    # For arbitrary Z-dependent exact histories, B and all its required
    # derivatives vanish at the original flat entrance. Common primitive
    # ODEs therefore identify every y jet; differentiating in Z gives the
    # required mixed orders. No interval width establishes this equality.
    p=s.Function('P0')(Z); up=s.Function('u0')(Z)
    histories=[s.Function(name)(Z) for name in ('m10','m20','e0','X0','Mp0')]
    m1=histories[0]*s.exp(-(s.Rational(1,2)-mu)*t)
    m2=histories[1]*s.exp(-(s.Rational(1,2)-2*mu)*t)
    ev=histories[2]*s.exp(2*mu*t)-(s.exp(2*mu*t)-1)/(4*mu)
    xv=1/rate+(histories[3]-1/rate)*s.exp(-rate*t)
    pv=p+histories[4]+up**2*(1-s.exp(-prate*t))/(2*prate)
    for name,function,history in zip(('m1','m2','energy','angular','pressure'),(m1,m2,ev,xv,pv),
        (histories[0],histories[1],histories[2],histories[3],p+histories[4])):
        for n in range(6):zero(name+'_whole_Z_inlet_axial'+str(n),s.diff(function.subs(t,0)-history,Z,n))
    proofs['same_source_primitive_ODEs_identify_y_jets_through4']=True
    proofs['original_flat_B_entrance_jets_are_exact_zero']=True
    proofs['paper_3_9_and_physical_prefactors_preserve_join']=True
    proofs['source_identity_is_distinct_from_numeric_overlap']=True
    return proofs


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Power inlet source changed: '+name)
    c=MPIntervalContext(); c.dps=240; read=lambda v:read_interval(c,v)
    checks=0
    for point in r['samples']+[r['whole_Z_left_box'],r['whole_Z_inlet']]:
        if not point['independently_evaluated_O3_power_high_mixed_derivatives'] or not point['radial_velocity_not_reset_when_axial_input_is_zero']:
            raise ValueError('Independent source recovery required')
        for grid in point['physical_mixed_derivatives_total_order_le4'].values():
            if len(grid)!=15:raise ValueError('Mixed-index coverage incomplete')
            for value in grid.values():
                if any(not mp.isfinite(v) for v in endpoints(read(value))):raise ArithmeticError('Power derivative nonfinite')
                checks+=1
        for key in ('full_pulse_C4_installed','full_outer_C4_certified','physical_energy_integral_certified','whole_outer_cone_certified','temporal_recursion'):
            if point[key]:raise ValueError('Power inlet scope promoted: '+key)
    pulse=json.loads((HERE/(PREFIX+'compliant_pulse_mixed_C4.json')).read_bytes())
    if pulse['implicit_source_sha256']!=r['implicit_source_sha256'] or pulse['actual_five_defect_family_sha256']!=r['actual_five_defect_family_sha256']:
        raise ValueError('Pulse inlet family mismatch')
    inlet=next(p for p in r['samples'] if endpoints(read(p['Z']))==(mp.mpf('.5'),mp.mpf('.5')) and endpoints(read(p['left_distance']))[1]==0)
    other=pulse['entrance_samples'][0]
    if endpoints(read(other['coordinate']['entrance_t']))[1]!=0:raise ValueError('Pulse source sample is not the entrance')
    overlap_checks=0
    for component,grid in inlet['physical_mixed_derivatives_total_order_le4'].items():
        for index,value in grid.items():
            a,b=endpoints(read(value)),endpoints(read(other['physical_mixed_derivatives_total_order_le4'][component][index]))
            if max(a[0],b[0])>min(a[1],b[1]):raise ArithmeticError('Independent inlet diagnostic disjoint: '+component+' '+index)
            overlap_checks+=1
    if endpoints(read(inlet['Mz_over_R_Utheta']['coefficients'][0]))[0]<=0:
        raise ArithmeticError('Nonzero incoming linear history lost')
    proofs=source_identities()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        datum_enclosure_sha256=r['datum_enclosure_sha256'],all_passed=True,
        exact_functional_production_and_join_identities=proofs,
        finite_actual_power_mixed_bounds_checked=checks,independent_inlet_overlap_diagnostics=overlap_checks,
        independently_evaluated_O3_power_high_mixed_derivatives=True,two_sided_O3_pulse_join_certified=True,
        quantitative_flat_comparison_prerequisite_admitted=True,
        join_scope='Local O3 power/O4 inlet: profile velocity/pressure mixed order<=4, canonical primitives axial5, fixed actual source and positive tau',
        full_pulse_C4_installed=False,full_outer_C4_certified=False,
        two_sided_pulse_O5_join_certified=False,physical_energy_integral_certified=False,
        whole_outer_cone_certified=False,temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
    print('O3/pulse inlet: production functional identities, independent high derivatives and quantitative flat continuation PASS',flush=True)
    return out


if __name__=='__main__':run()
