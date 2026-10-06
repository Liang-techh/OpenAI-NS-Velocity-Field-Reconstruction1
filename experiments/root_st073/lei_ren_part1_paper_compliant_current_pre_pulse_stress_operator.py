"""Full paper stress on actual raw pre-pulse histories and variable jets.

Only R powers and the constant Pstar are factored out. No pure-power
amplitude rate or vanishing meridional history is substituted.
"""
import ast
import copy
import hashlib
import math
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_current_pulse_entrance_incoming_background_tensor import (
    HERE,PREFIX,SourceAST,sha,copy_jet,IntervalTaylor,endpoints)
from lei_ren_part1_paper_compliant_collar_stress_C3 import axial_derivative,shifted_rows
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import exponential_derivatives
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import (
    original_full_stress,axial_operator_rows)


def raw_pre_stress_rows(c,delta,z,u,V,history,P):
    """Raw m,h,k,e,p units; ordinary rows include their full y derivatives."""
    m,h,k,e=(history[name] for name in ('m','h','k','e'))
    d=1-z*z;L=1-z*z*delta;b=(1-delta)/2;kk=1-delta/2
    transport=[z*m[j]*(1-delta)+d*axial_derivative(m[j]) for j in range(5)]
    nonlinear_theta=product_rows(u,transport)
    nonlinear_z=product_rows(V,transport)
    theta={
        'local_transport':dict(mode=(.5,1,0,0),shape=[-u[j]/L for j in range(4)]),
        'retained_angular_moment':dict(mode=(.5,1,0,0),shape=[(h[j]*kk-z*b*axial_derivative(h[j]))/L for j in range(4)]),
        'retained_mixed_moment':dict(mode=(.5,1,0,0),shape=[(z*(2*delta-1)*k[j]-d*axial_derivative(k[j]))/L for j in range(4)]),
        'meridional_transport':dict(mode=(.5,1,0,0),shape=[nonlinear_theta[j]/L for j in range(4)]),
        'variable_radial_shear':dict(mode=(-.5,1,0,0),shape=[2*u[j+1]-u[j] for j in range(4)])}
    axial={
        'local_axial_transport':dict(mode=(.5,0,0,0),shape=[-V[j]/L for j in range(4)]),
        'nonlinear_meridional_transport':dict(mode=(.5,0,0,0),shape=[nonlinear_z[j]/L for j in range(4)]),
        'retained_linear_axial_moment':dict(mode=(.5,0,0,0),shape=[(m[j]-z*axial_derivative(m[j]))*(1-delta)/(2*L) for j in range(4)]),
        'retained_full_energy':dict(mode=(.5,2,0,0),shape=[(z*(2*delta)*e[j]-d*axial_derivative(e[j]))/L for j in range(4)]),
        'actual_absolute_pressure':dict(mode=(.5,2,0,0),shape=[(z*(2*(1+delta))*P[j]-d*axial_derivative(P[j]))/L for j in range(4)]),
        'axial_radial_shear':dict(mode=(-.5,0,0,0),shape=[2*V[j+1] for j in range(4)])}
    for parts in (theta,axial):
        for part in parts.values():
            part['full_derivative_rows']=shifted_rows(part['shape'],c.mpf(str(part['mode'][0])))
    return dict(theta=theta,axial=axial)


def raw_pre_velocity_rows(c,pre):
    raw=pre['actual_normalized_primitive_y_derivative_axial5']
    u0=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in pre['Utheta_over_Pstar_axial5_coefficients']])
    logjet=[copy_jet(c,value) for value in pre['log_Utheta_ordinary_y_derivatives']]
    u=[u0*value for value in exponential_derivatives(logjet)]
    V=[copy_jet(c,value) for value in pre['Uz_ordinary_y_derivative_axial5']]
    Q=[copy_jet(c,value) for value in pre['actual_Q_y_derivative_axial4']]
    radial=shifted_rows(Q,c.mpf('.5'),4)
    return dict(theta=u,axial=V,radial=radial)


def raw_pre_remainder_sectors(c,delta,z,velocity):
    """Same full NS remainder, with the actual raw amplitude units."""
    Ur,Ut,Uz=(velocity[key] for key in ('radial','theta','axial'))
    L=1-z*z*delta
    time=[(Ur[j]/2+z*axial_derivative(Ur[j])*((1-delta)/2)+Ur[j+1])/L for j in range(3)]
    radial_viscosity=shifted_rows([-(Ur[j+2]-Ur[j]/4)*2 for j in range(3)],-1,2)
    dzUr=axial_operator_rows(c,Ur,-1,z,delta,count=2,order=1)
    radial_transport=shifted_rows(product_rows(Ur[:3],Ur[1:4]),-c.mpf('.5'),2)
    axial_transport=product_rows(Uz[:3],dzUr)
    nonlinear=[radial_transport[j]+axial_transport[j] for j in range(3)]
    return dict(radial=dict(
        time=dict(mode=(.5,0,0,0),normalization_half=True,beta=-3,rows=time),
        radial_viscosity=dict(mode=(-.5,0,0,0),normalization_half=True,beta=-3,rows=radial_viscosity),
        nonlinear_transport=dict(mode=(.5,0,0,0),normalization_half=True,beta=-3,rows=nonlinear),
        axial_viscosity=dict(mode=(.5,0,0,0),normalization_half=True,beta=-3+2*delta,
            rows=[-value for value in axial_operator_rows(c,Ur,-1,z,delta)])),
        theta=dict(axial_viscosity=dict(mode=(0,1,0,0),normalization_half=False,beta=-3+delta,
            rows=[-value for value in axial_operator_rows(c,Ut,-1-delta,z,delta)])),
        axial=dict(axial_viscosity=dict(mode=(0,0,0,0),normalization_half=False,beta=-3+delta,
            rows=[-value for value in axial_operator_rows(c,Uz,-1-delta,z,delta)])))


def raw_pre_formula_theorem():
    """Replay original (3.16)-(3.18) with arbitrary variable source functions."""
    y,z,delta=s.symbols('y Z delta',real=True);R0,Ps=s.symbols('R0 Pstar',positive=True)
    R=R0*s.exp(y);d=1-z*z;L=1-delta*z*z
    names=('u','V','m','h','k','e','p')
    f={name:s.Function('actual_raw_'+name)(y,z) for name in names}
    rows={key:[s.diff(value,y,j) for j in range(5)] for key,value in f.items()}
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)))
    asts=SourceAST();env=dict(math=math,axial_derivative=lambda value:s.diff(value,z),product_rows=product_rows,shifted_rows=shifted_rows)
    operator=asts.replay('current_pre_pulse_stress_operator','raw_pre_stress_rows',env)
    result=operator(c,delta,z,rows['u'],rows['V'],rows,rows['p'])
    moments=dict(theta=s.sqrt(2)*R**s.Rational(3,2)*Ps*f['h'],z=R*f['m'],
        theta_z=s.sqrt(2)*R**s.Rational(3,2)*Ps*f['k'],z_theta=R*Ps**2*f['e'],p=Ps**2*f['p'])
    actual,hashes=original_full_stress(dict(a=Ps*f['u'],b=f['V'],ay=Ps*rows['u'][1],by=rows['V'][1],
        az=Ps*s.diff(f['u'],z),bz=s.diff(f['V'],z),dt=delta,z=z,R=R,root=s.sqrt(2*R),L=L,d=d,
        m=moments,mz={key:s.diff(value,z) for key,value in moments.items()},p=Ps**2*f['p'],pz=Ps**2*s.diff(f['p'],z)))
    checks={}
    for label,paper in dict(theta=actual['Itheta']+actual['Stheta'],axial=actual['Iz']+actual['Sz']).items():
        for j in range(4):
            lifted=sum(R**s.Rational(str(part['mode'][0]))*Ps**part['mode'][1]*part['full_derivative_rows'][j]/s.sqrt(2) for part in result[label].values())
            difference=s.simplify(s.expand_power_exp(lifted-s.diff(paper,y,j)))
            if difference!=0:raise ArithmeticError('Actual variable raw paper stress row differs: '+label+str(j))
            for n in range(4-j):
                if s.diff(difference,z,n)!=0:raise ArithmeticError('Variable raw stress mixed derivative differs')
                checks[label+'_y'+str(j)+'_Z'+str(n)]=True
    return dict(identities=checks,original_full_paper_stress_AST_replayed=True,
        arbitrary_variable_log_amplitude_and_all_five_histories_supported=True,
        only_exact_radius_powers_have_constant_shift_rates=True,
        no_positive_amplitude_or_radius_materialized=True,
        input_hashes={**asts.hashes,**hashes},passed=True)


def compiled_raw_pre_physical_lift():
    """Reuse the full physical mapper; change only the remainder units."""
    import lei_ren_part1_paper_compliant_pulse_end_physical_C2 as original
    asts=SourceAST();fn=copy.deepcopy(asts.method('pulse_end_physical_C2','lift_physical_packet'))
    fn.decorator_list=[];changes=[]
    for node in ast.walk(fn):
        if isinstance(node,ast.Assign) and ast.unparse(node)=='source_errors = pulse_remainder_sectors(c, delta, None, z, velocity)':
            node.value=ast.parse('raw_pre_remainder_sectors(c,delta,z,velocity)',mode='eval').body;changes.append('source_errors')
    if changes!=['source_errors']:raise ValueError('Unreviewed raw pre physical lift adaptation')
    env=dict(vars(original));env['raw_pre_remainder_sectors']=raw_pre_remainder_sectors
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<original full physical lift; actual raw amplitude units>','exec'),env)
    name=PREFIX+'current_pre_pulse_stress_operator.py';asts.hashes[name]=sha(name)
    return env[fn.name],dict(original_stress_divergence_completion_cartesian_and_physical_operators_unchanged=True,
        only_source_remainder_amplitude_modes_adapted=True,input_hashes=asts.hashes)


def raw_pre_velocity_remainder_theorem():
    """Replay variable-jet recovery and unchanged remainder coefficient operators."""
    from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z
    from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_operators
    class Projection(FunctionJet):
        def __pow__(self,n):return self.function(self.ctx,self.expr**n,self.order)
    c=SimpleNamespace(mpf=lambda v:v[0] if isinstance(v,(tuple,list)) else s.Rational(str(v)) if isinstance(v,(str,int,float)) else v)
    jet=lambda value,order=5:Projection.function(c,value,order)
    u=jet(s.Function('same_actual_u')(Z));logjet=[jet(s.Function('same_actual_logjet'+str(j))(Z)) for j in range(4)]
    V=[jet(s.Function('same_actual_Vjet'+str(j))(Z)) for j in range(5)]
    histories={name:jet(s.Function('same_actual_'+name)(Z)) for name in ('m','h','k','e','p')}
    delta,invP2=s.symbols('delta invPstar_squared',real=True)
    asts=SourceAST();env=dict(math=math,IntervalTaylor=Projection,
        square=lambda v:v*v,derivative=lambda v:jet(s.diff(v.expr,Z),v.order-1),
        axial_derivative=lambda v:jet(s.diff(v.expr,Z),v.order-1),
        copy_jet=lambda ctx,v:v,endpoints=lambda v:(v,v),physical_operators=physical_operators)
    for stem,name in (('pre_pulse_mixed_C4','rate_rows'),('pre_pulse_mixed_C4','product_rows'),
            ('long_reshape_mixed_C4','exponential_derivatives')):asts.replay(stem,name,env)
    source=asts.replay('pre_pulse_mixed_C4','physical_mixed',env)(c,Z,delta,u,logjet,V,histories,jet(s.Function('same_actual_P0')(Z)),invP2)
    packet=dict(source,Utheta_over_Pstar_axial5_coefficients=u.coefficients,
        log_Utheta_ordinary_y_derivatives=logjet,Uz_ordinary_y_derivative_axial5=V)
    env['shifted_rows']=shifted_rows
    velocity=asts.replay('current_pre_pulse_stress_operator','raw_pre_velocity_rows',env)(c,packet)
    checks={}
    for label,key in dict(theta='Utheta_over_Pstar',axial='Uz',radial='Ur_over_current_sqrt_R_over_2').items():
        source_rows=source['physical_velocity_pressure_y_derivative_Taylor'][key]
        for j in range(5):
            difference=s.cancel(velocity[label][j].expr-source_rows[j].expr)
            if difference!=0:raise ArithmeticError('Variable actual source velocity recovery differs: '+label+str(j))
            for n in range(5-j):checks[label+'_y'+str(j)+'_Z'+str(n)]=s.diff(difference,Z,n)==0
    # Compare the actual coefficient program with the original full
    # remainder. Only amplitude modes change: radial has no Pstar, theta
    # has one Pstar, and the local axial source has none.
    env['product_rows']=product_rows
    env['axial_n']=lambda value,order:jet(s.diff(value.expr,Z,order),value.order-order)
    asts.replay('pulse_end_physical_C2','axial_operator_rows',env)
    left=asts.replay('current_pre_pulse_stress_operator','raw_pre_remainder_sectors',env)(c,delta,jet(Z),velocity)
    right=asts.replay('pulse_end_physical_C2','pulse_remainder_sectors',env)(c,delta,None,jet(Z),velocity)
    for label,parts in left.items():
        for name,part in parts.items():
            original=right[label][name]
            if (part['mode'][0],part['normalization_half'],part['beta'])!=(original['mode'][0],original['normalization_half'],original['beta']):raise ArithmeticError('Original remainder radius/physical units differ')
            for j,(a,b) in enumerate(zip(part['rows'],original['rows'])):
                difference=s.cancel(a.expr-b.expr)
                if difference!=0:raise ArithmeticError('Original full remainder coefficient changed')
                for n in range(3-j):checks[label+'_'+name+'_y'+str(j)+'_Z'+str(n)]=s.diff(difference,Z,n)==0
    return dict(identities=checks,original_source_physical_mixed_and_remainder_operators_replayed=True,
        arbitrary_log_amplitude_derivatives_and_local_axial_velocity_jets_retained=True,
        radial_sqrt_R_over_2_theta_Pstar_axial_one_amplitude_modes=True,
        input_hashes=asts.hashes,passed=True)
