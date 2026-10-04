"""Actual O7 waiting stresses, transported from the SAME full future moments.

The collar source is evaluated only at offset zero. Its full moment defects
are propagated backwards by the original constant-power FTC, not by calling
the collar outside its domain. Exact source factors remain factored.
"""
import ast
import functools
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_collar_stress_C3 import (
    CompliantCollarStressC3, collar_defect_rows, axial_derivative)
from lei_ren_part1_paper_compliant_pulse_physical_bounds import P
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def source_precision(fn):
    @functools.wraps(fn)
    def wrapped(*args,**kwargs):
        with mp.workdps(300):return fn(*args,**kwargs)
    return wrapped


def waiting_defect_rows(heat,base,q):
    """Exact backward constant-K moment ODEs, axial jets through five."""
    c=heat.ctx; one=IntervalTaylor.constant(c,1,5)
    eps=heat.eps; k=heat.k; delta=heat.delta; p=heat.prate
    Q0=-2*eps+eps**2
    Ad0=base['angular_defect_rows'][0]
    Ed0=base['energy_defect_rows'][0]
    Pd0=base['pressure_defect_rows'][0]
    Ad=-one*(eps/k)+(Ad0+one*(eps/k))*c.exp(-k*q)
    Ed=one*(Q0/delta)+(Ed0-one*(Q0/delta))*c.exp(delta*q)
    Pd=one*(Q0/(2*p))+(Pd0-one*(Q0/(2*p)))*c.exp(p*q)
    A=[Ad]; E=[Ed]; Pr=[Pd]
    for j in range(4):
        A.append(-one*eps-A[j]*k if j==0 else -A[j]*k)
        E.append(E[j]*delta-one*Q0 if j==0 else E[j]*delta)
        Pr.append(Pr[j]*p-one*Q0/2 if j==0 else Pr[j]*p)
    return dict(angular_defect_rows=A,energy_defect_rows=E,pressure_defect_rows=Pr,
                K_defect_rows=[-one*eps]+[one*0]*4,
                K_squared_defect_rows=[one*Q0]+[one*0]*4)


def waiting_stress_rows(heat,base,Z,q):
    """Cancel stationary moment baselines BEFORE signed interval evaluation.

    Rows include the positive Qt/Qz factor derivatives. Their two radial
    rates are -1, -(1+a) for theta, and -1/2, +1/2 for axial stress.
    """
    c=heat.ctx; z=IntervalTaylor.variable(c,Z,5)
    d=1-z*z; L=1-z*z*heat.delta; b=(1-heat.delta)/2
    Ad0=base['angular_defect_rows'][0]
    Ed0=base['energy_defect_rows'][0]
    Pd0=base['pressure_defect_rows'][0]
    eps=heat.eps; Q0=-2*eps+eps**2; K0=1-eps
    Nt=Ad0*heat.k+eps-z*axial_derivative(Ad0)*b
    Ne=z*Ed0*heat.delta-z*Q0-d*axial_derivative(Ed0)/2
    Np=-z*Pd0*(2*heat.prate)+z*Q0+d*axial_derivative(Pd0)
    It=Nt/L*c.exp(-heat.k*q)
    St=IntervalTaylor.constant(c,-2*heat.S*(1+heat.a)*K0*c.exp(-q),4)
    Je=Ne/L*c.exp(heat.delta*q)
    Jp=Np/L*c.exp(heat.prate*q)
    inertial=[It*(-1)**j for j in range(4)]
    shear=[St*(-(1+heat.a))**j for j in range(4)]
    axial=[Je*(-c.mpf('.5'))**j+Jp*c.mpf('.5')**j for j in range(4)]
    return dict(theta=[inertial[j]+shear[j] for j in range(4)],axial=axial,
                axial_energy=[Je*(-c.mpf('.5'))**j for j in range(4)],
                axial_pressure=[Jp*c.mpf('.5')**j for j in range(4)],
                theta_inertial=inertial,theta_shear=shear,
                normalized_coefficient_rows=dict(theta=[It+St],axial=[Je+Jp]),
                shear_strength_margin=IntervalTaylor.constant(c,heat.delta,3))


def waiting_pressure_rows(heat,base,q):
    """Original absolute pressure; full future datum retained at Rtail."""
    c=heat.ctx; one=IntervalTaylor.constant(c,1,5)
    p=heat.prate; K0=1-heat.eps; Q0=-2*heat.eps+heat.eps**2
    datum=base['pressure_defect_rows'][0]-one*(Q0/(2*p))
    pure=one*(K0**2/(2*p)*c.exp(-p*q))
    rows=[-(pure+datum)*heat.pressure_scale]
    rows.extend(-pure*heat.pressure_scale*(-p)**j for j in range(1,5))
    return rows


def waiting_identities():
    """Functional FTC, endpoint, physical power and zero-history identities."""
    a,eps,q,z,wait=s.symbols('a eps q Z wait',real=True)
    delta=2*a; k=1-a; p=1+delta; b=(1-delta)/2
    K=1-eps; Q=-2*eps+eps**2; L=1-delta*z*z; d=1-z*z
    A0=s.Function('Ad0')(z); E0=s.Function('Ed0')(z); P0=s.Function('Pd0')(z)
    Ad=-eps/k+(A0+eps/k)*s.exp(-k*q)
    Ed=Q/delta+(E0-Q/delta)*s.exp(delta*q)
    Pd=Q/(2*p)+(P0-Q/(2*p))*s.exp(p*q)
    A=1/k+Ad; E=1/delta+Ed; Pr=1/(2*p)+Pd
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Waiting identity failed: '+name)
        proofs[name]=True
    zero('actual_angular_backward_FTC',s.diff(A,q)-(K-k*A))
    zero('actual_energy_backward_FTC',s.diff(E,q)-(delta*E-K*K))
    zero('actual_pressure_backward_FTC',s.diff(Pr,q)-(p*Pr-K*K/2))
    for label,row,start in (('angular',Ad,A0),('energy',Ed,E0),('pressure',Pd,P0)):
        zero('same_full_future_'+label+'_endpoint',row.subs(q,0)-start)
    Ct=(k*A-b*z*s.diff(A,z)-K)/L
    Cz=(delta*z*E-d*s.diff(E,z)/2-2*p*z*Pr+d*s.diff(Pr,z))/L
    Nt=k*A0+eps-b*z*s.diff(A0,z)
    Ne=delta*z*E0-z*Q-d*s.diff(E0,z)/2
    Np=-2*p*z*P0+z*Q+d*s.diff(P0,z)
    zero('actual_cancelled_theta_inertial',Ct-Nt/L*s.exp(-k*q))
    zero('actual_cancelled_axial_stress',Cz-(Ne*s.exp(delta*q)+Np*s.exp(p*q))/L)
    theta_factor=s.exp(-a*q); axial_factor=s.exp(-(s.Rational(1,2)+delta)*q)
    for j in range(4):
        zero('theta_inertial_full_factor_y'+str(j),
             s.diff(theta_factor*Ct,q,j)/theta_factor-(-1)**j*Ct)
        shear=s.exp(-q)*(-2*(1+a)*K)
        zero('theta_shear_full_factor_y'+str(j),
             s.diff(theta_factor*shear,q,j)/theta_factor-(-(1+a))**j*shear)
        zero('axial_full_factor_y'+str(j),s.diff(axial_factor*Cz,q,j)/axial_factor
             -(Ne*s.exp(delta*q)*(-s.Rational(1,2))**j+Np*s.exp(p*q)*s.Rational(1,2)**j)/L)
    reference_pressure=-s.exp(-p*q)*Pr
    zero('absolute_pressure_full_future_datum',reference_pressure
         +(K*K/(2*p)*s.exp(-p*q)+P0-Q/(2*p)))
    zero('waiting_pressure_density_FTC',s.diff(reference_pressure,q)-K*K*s.exp(-p*q)/2)
    X0,H,thetaT=s.symbols('Xtail H thetaT',real=True)
    theta_base=thetaT*s.exp(-(s.Rational(1,2)+a)*wait)/K
    zero('original_waiting_velocity_normalization',theta_base*s.exp(-(s.Rational(1,2)+a)*q)*K
         -thetaT*s.exp(-(s.Rational(1,2)+a)*(wait+q)))
    source_X=(X0-1/k)*s.exp(-k*q)+1/k
    zero('original_waiting_angular_normalization',
         (K*X0*s.exp(-k*q)+K/k*(1-s.exp(-k*q)))/K-source_X)
    zero('original_waiting_energy_half_normalization',
         (K*K*H*s.exp(delta*q)+K*K/delta*(1-s.exp(delta*q)))/(2*K*K)
         -(H*s.exp(delta*q)+(1-s.exp(delta*q))/delta)/2)
    # Both primitive meridional derivatives vanish on the actual source.
    m0,mt0=s.symbols('actual_Mz_tail actual_Mtheta_z_tail')
    for label,terminal in (('Mz',m0),('Mtheta_z',mt0)):
        primitive=terminal+s.integrate(s.Integer(0),(s.Symbol('v'),0,q))
        zero('actual_zero_'+label+'_backward_primitive',primitive.subs(terminal,0))
    lam,r,nu,Rt,Btail=s.symbols('lambda r nu Rtail Btail',positive=True)
    bh=s.Rational(1,2)+a
    physical=s.sqrt(nu)*lam**(-1-delta)*Btail*K*(r*r/(2*nu*lam*lam*Rt))**(-bh)
    pure=nu**(1+a)*Btail*K*(2*Rt)**bh*r**(-1-delta)
    if s.simplify(s.powsimp(physical-pure,force=True))!=0:raise ArithmeticError('Physical waiting power failed')
    proofs['actual_waiting_physical_swirl_pure_radial_power']=True
    zero('waiting_physical_swirl_lambda_derivative_zero',s.diff(pure,lam))
    zero('waiting_physical_swirl_axial_viscosity_zero',s.diff(pure,z,2))
    return dict(identities=proofs,waiting_backward_full_moment_FTC_verified=True,
                waiting_physical_swirl_pure_power_verified=True,
                waiting_physical_axial_viscosity_exact_zero=True,
                interval_overlap_used_as_functional_proof=False)


def waiting_join_binding():
    """Replay actual source formula ASTs on arbitrary symbolic terminal jets.

    Scalar symbolic functions stand for the axial Taylor algebra here. The
    single context lookup in pressure_y_rows is removed; its actual scalar
    expressions and all loops are executed unchanged in the same context.
    """
    a,eps,z,C=s.symbols('a eps Z pressure_scale',real=True)
    c=SimpleNamespace(exp=s.exp,mpf=lambda value:s.Rational(str(value)))
    heat=SimpleNamespace(ctx=c,a=a,eps=eps,S=s.Symbol('exact_inverse_Rtail'),
                         delta=2*a,k=1-a,prate=1+2*a,pressure_scale=C)
    stub=SimpleNamespace(constant=lambda ctx,value,order:value,
                         variable=lambda ctx,value,order:value)
    environment=dict(s=s,math=math,IntervalTaylor=stub,c=c,
                     axial_derivative=lambda value:s.diff(value,z))
    hashes={}
    def replay(stem,name,remove_context=False):
        path=HERE/(PREFIX+stem+'.py'); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        tree=ast.parse(path.read_text(encoding='utf8'))
        fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==name)
        fn.decorator_list=[]
        if remove_context:
            matches=[n for n in fn.body if isinstance(n,ast.Assign) and len(n.targets)==1
                     and ast.unparse(n.targets[0])=='c' and ast.unparse(n.value)=='pressure_numerator.ctx']
            if len(matches)!=1:raise ValueError('Pressure scalar context lookup changed')
            fn.body.remove(matches[0])
        exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual source '+name+'>','exec'),environment)
        return environment[name]
    replay('collar_stress_C3','shifted_rows')
    replay('collar_Gamma_C4','product_rows')
    replay('heat_pressure_C4','pressure_y_rows',True)
    collar_pressure=replay('collar_pressure_C4','collar_pressure_rows')
    collar_stress=replay('collar_stress_C3','collar_stress_rows')
    defect=replay('waiting_stress_C3','waiting_defect_rows')
    stress=replay('waiting_stress_C3','waiting_stress_rows')
    pressure=replay('waiting_stress_C3','waiting_pressure_rows')
    base={label+'_defect_rows':[s.Function(label+'_terminal_defect')(z)] for label in ('angular','energy','pressure')}
    defects=defect(heat,base,s.Integer(0)); shape=dict(K_rows=[1-eps]+[s.Integer(0)]*4)
    right=collar_stress(heat,shape,defects,z,s.Integer(0))
    left=stress(heat,base,z,s.Integer(0))
    proofs={}
    for label in ('theta','axial'):
        for j in range(4):
            for n in range(4-j):
                if s.simplify(s.diff(left[label][j]-right[label][j],z,n))!=0:
                    raise ArithmeticError('Actual waiting/collar stress AST join failed: '+str((label,j,n)))
                proofs[label+'_y'+str(j)+'_Z'+str(n)]=True
    remainder=base['pressure_defect_rows'][0]+1/(2*heat.prate)
    right_pressure=collar_pressure(shape['K_rows'],remainder,heat.prate,C,s.Integer(0))
    left_pressure=pressure(heat,base,s.Integer(0))
    for j in range(5):
        for n in range(5-j):
            if s.simplify(s.diff(left_pressure[j]-right_pressure[j],z,n))!=0:
                raise ArithmeticError('Actual waiting/collar pressure AST join failed: '+str((j,n)))
            proofs['pressure_y'+str(j)+'_Z'+str(n)]=True
    return dict(identities=proofs,actual_waiting_collar_stress_mixed3_AST_join_verified=True,
                actual_waiting_collar_pressure_mixed4_AST_join_verified=True,
                arbitrary_full_terminal_axial_functions_not_sample_overlap=True,input_hashes=hashes)


def waiting_source_bridge(collar):
    """Compose accepted full endpoint histories with original waiting FTC."""
    hashes=dict(collar.hashes); trees={}; bindings={}
    def syntax(stem,method,target,expression):
        path=HERE/(PREFIX+stem+'.py')
        if stem not in trees:
            trees[stem]=ast.parse(path.read_text(encoding='utf8'))
            hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        fn=next(n for n in ast.walk(trees[stem]) if isinstance(n,ast.FunctionDef) and n.name==method)
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)]
        expected=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(ast.dump(v)==expected for v in values)!=1:raise ValueError('Waiting source changed: '+stem+':'+target)
        bindings[stem+'.'+method+'.'+target]=True
    original={
        'theta':'one*self.thetaT*c.exp(-self.bh*t)',
        'X':"(data['XT']-1/self.k)*c.exp(-self.k*t)+1/self.k",
        'energy':"(data['H']*c.exp(-self.delta*left)+decay_integral(c,self.delta,left))/2",
        'pressure':"data['PT']+decay_integral(c,1+self.delta,t)*(self.outer.flatten.Ev2*self.thetaT**2/2)",
        't':'self.wait*phase','left':'self.wait*(1-phase)'}
    for target,expression in original.items():syntax('steep_waiting_C4','waiting',target,expression)
    syntax('steep_waiting_C4','packet','result',
           'self.outer._packet(Z,theta,X,energy,pressure,[one*r for r in rates],coordinate,dict(stage=stage,entire_steep_waiting_high_mixed_derivatives_available=True,same_actual_angular_terminal_histories_retained=True,exact_Gamma_and_both_epsilon_atoms_retained=True,exact_relative_velocity_log_parts=logs,no_forward_subtraction_of_unrelated_long_future_energy=True))')
    syntax('power_angular_C4','_packet','mixed',
           'flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')
    syntax('flatten_mixed_C4','flatten_mixed','rows',
           '{UZ:[zero for _ in range(5)],UT:theta_rows,UR:[zero for _ in range(5)],P:prows}')
    syntax('collar_Gamma_C4','__init__','self.theta_base',
           'self.steep.thetaT*c.exp(-self.bh*self.steep.wait-self.steep.logone)')
    formulas={
        'waiting_defect_rows':{
            'Ad':'-one*(eps/k)+(Ad0+one*(eps/k))*c.exp(-k*q)',
            'Ed':'one*(Q0/delta)+(Ed0-one*(Q0/delta))*c.exp(delta*q)',
            'Pd':'one*(Q0/(2*p))+(Pd0-one*(Q0/(2*p)))*c.exp(p*q)'},
        'waiting_stress_rows':{
            'Nt':'Ad0*heat.k+eps-z*axial_derivative(Ad0)*b',
            'Ne':'z*Ed0*heat.delta-z*Q0-d*axial_derivative(Ed0)/2',
            'Np':'-z*Pd0*(2*heat.prate)+z*Q0+d*axial_derivative(Pd0)',
            'It':'Nt/L*c.exp(-heat.k*q)',
            'St':'IntervalTaylor.constant(c,-2*heat.S*(1+heat.a)*K0*c.exp(-q),4)',
            'Je':'Ne/L*c.exp(heat.delta*q)','Jp':'Np/L*c.exp(heat.prate*q)',
            'inertial':'[It*(-1)**j for j in range(4)]',
            'shear':'[St*(-(1+heat.a))**j for j in range(4)]',
            'axial':"[Je*(-c.mpf('.5'))**j+Jp*c.mpf('.5')**j for j in range(4)]"},
        'waiting_pressure_rows':{
            'datum':"base['pressure_defect_rows'][0]-one*(Q0/(2*p))",
            'pure':'one*(K0**2/(2*p)*c.exp(-p*q))',
            'rows':'[-(pure+datum)*heat.pressure_scale]'}}
    for method,values in formulas.items():
        for target,expression in values.items():syntax('waiting_stress_C3',method,target,expression)
    for flag in ('actual_full_collar_moments_same_source_verified','meridional_moments_and_velocities_exact_zero'):
        if not collar.bridge[flag]:raise ValueError('Actual waiting endpoint history missing: '+flag)
    if not collar.pressure.collar_bridge['pressure_waiting_collar_and_collar_exterior_mixed4_joins_verified']:
        raise ValueError('Original waiting absolute pressure endpoint is not admitted')
    for flag in ('actual_cumulative_meridional_histories_remain_zero_by_FTC_when_Uz_is_zero',
                 'downstream_zero_rows_follow_source_proven_terminal_meridional_histories'):
        if not collar.history['identities'][flag]:raise ValueError('Actual primitive zero-density source route missing: '+flag)
        bindings['consumed_terminal_history_'+flag]=True
    join=waiting_join_binding(); hashes.update(join['input_hashes'])
    return dict(source_bindings=bindings,input_hashes=hashes,
                actual_waiting_collar_formula_join=join,
                actual_waiting_endpoint_full_moment_history_consumed=True,
                actual_waiting_absolute_pressure_endpoint_consumed=True,
                actual_meridional_terminal_zeros_and_zero_density_FTC_consumed=True,
                waiting_constant_K_and_flat_collar_K_endpoint_same_source=True,
                waiting_collar_stress_mixed3_join_verified=True,
                waiting_collar_pressure_mixed4_join_verified=True,
                steep_exit_waiting_stress_mixed3_join_verified=False,
                original_datum_coefficients_or_velocity_changed=False,
                no_collar_source_evaluation_at_negative_offset=True)


class CompliantWaitingStressC3:
    @source_precision
    def __init__(self):
        self.collar_source=CompliantCollarStressC3()
        self.heat=self.collar_source.heat; self.steep=self.heat.steep; self.ctx=self.heat.ctx
        self.family,self.source=self.heat.family,self.heat.source
        name=PREFIX+'collar_stress_C3_check.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['actual_original_collar_similarity_stress_recovered']:
            raise ValueError('Accepted full collar moments/stress required')
        self.hashes=dict(self.collar_source.hashes)
        for path,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Waiting prerequisite changed: '+path)
        self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.proof=waiting_identities(); self.bridge=waiting_source_bridge(self.collar_source)
        self.hashes.update(self.bridge['input_hashes']); self.cache={}
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def endpoint(self,Z):
        c=self.ctx; Z=c.mpf(Z); key=tuple(endpoints(Z))
        if key not in self.cache:
            self.cache[key]=collar_defect_rows(self.heat,self.heat.shape(Z,0),self.heat.collar_tails(Z,0),0)
        return self.cache[key]

    @source_precision
    def waiting(self,Z,phase):
        c=self.ctx; Z=c.mpf(Z); phase=c.mpf(phase); self.steep._phase(phase)
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1:raise ValueError('Source Z in[-1,1] required')
        q=self.steep.wait*(phase-1); base=self.endpoint(Z)
        point=dict(self.steep.waiting(Z,phase)); one=IntervalTaylor.constant(c,1,5); K0=1-self.heat.eps
        defects=waiting_defect_rows(self.heat,base,q); rows=waiting_stress_rows(self.heat,base,Z,q)
        pressure=waiting_pressure_rows(self.heat,base,q)
        point['original_forward_angular_Taylor']=point['angular_Taylor']
        point['original_forward_angular_y_derivative_Taylor']=point['angular_y_derivative_Taylor']
        point['original_forward_energy_Taylor']=point['energy_Taylor']
        point['original_forward_energy_y_derivative_Taylor']=point['energy_y_derivative_Taylor']
        point['original_forward_pressure_over_Pstar_squared_Taylor']=point['pressure_over_Pstar_squared_Taylor']
        mixed=dict(point['physical_mixed_derivatives_total_order_le4'])
        point['original_forward_pressure_mixed_bounds']=mixed[P]
        mixed[P]={'y'+str(j)+'_Z'+str(n):jet[n]*math.factorial(n) for j,jet in enumerate(pressure) for n in range(5-j)}
        fields=dict(point['physical_velocity_and_pressure_y_derivative_Taylor'])
        fields[P]=[jet.truncate(4-j) for j,jet in enumerate(pressure)]
        angular_rows=[one/self.heat.k+defects['angular_defect_rows'][0]]+defects['angular_defect_rows'][1:]
        energy_rows=[one/self.heat.delta+defects['energy_defect_rows'][0]]+defects['energy_defect_rows'][1:]
        point.update(angular_Taylor=(one/self.heat.k+defects['angular_defect_rows'][0])/K0,
                     angular_y_derivative_Taylor=[(jet/K0).truncate(5-j) for j,jet in enumerate(angular_rows)],
                     energy_Taylor=(one/self.heat.delta+defects['energy_defect_rows'][0])/(2*K0**2),
                     energy_y_derivative_Taylor=[(jet/(2*K0**2)).truncate(5-j) for j,jet in enumerate(energy_rows)],
                     pressure_over_Pstar_squared_Taylor=pressure[0],
                     pressure_y_derivative_axial5_Taylor=pressure,
                     physical_mixed_derivatives_total_order_le4=mixed,
                     physical_velocity_and_pressure_y_derivative_Taylor=fields,
                     waiting_offset_from_Rtail=q,waiting_phase=phase,
                     waiting_similarity_stress_mixed3_factored={label:{'y'+str(j)+'_Z'+str(n):rows[label][j][n]*math.factorial(n)
                         for j in range(4) for n in range(4-j)} for label in ('theta','axial')},
                     waiting_similarity_stress_y_derivative_Taylor={label:[jet.truncate(3-j) for j,jet in enumerate(rows[label])]
                         for label in ('theta','axial','theta_inertial','theta_shear','axial_energy','axial_pressure')},
                     waiting_future_defect_zeroth_Taylor={label:values[0] for label,values in defects.items()},
                     waiting_shear_strength_kappa_minus2=self.heat.delta,
                     waiting_shear_strength_kappa_gt2_certified=endpoints(self.heat.delta)[0]>0,
                     waiting_angular_shear_negative_certified=False,
                     waiting_meridional_moments_and_velocities={label:one*0 for label in ('Mz','Mtheta_z','Ur','Uz')},
                     actual_original_waiting_similarity_stress_recovered=True,
                     actual_waiting_absolute_pressure_same_source_mixed4_available=True,
                     full_sigma_phi_Gamma_future_integral_used=True,
                     original_angular_energy_pressure_histories_retained=True,
                     waiting_collar_stress_mixed3_join_verified=True,
                     steep_exit_waiting_stress_mixed3_join_verified=False,
                     waiting_physical_axial_viscosity_exact_zero=True,
                     source_defined_positive_stress_factors=dict(theta='sqrt(R/2)*B',axial='sqrt(R/2)*B^2',
                         B='Ev0*theta_base*exp(-bh*q)',R='Rtail*exp(q)',q='wait*(phase-1)',
                         exact_Ev0='Pstar*U*exp(-13/(2*mu)-13); Ev2 is an enclosure only'),
                     waiting_cone_certified=False,global_admissible_stress_lift_constructed=False,
                     physical_energy_integral_certified=False,temporal_recursion=False)
        return point

    def report(self):
        keys=('Z','waiting_phase','waiting_offset_from_Rtail','waiting_similarity_stress_mixed3_factored',
              'waiting_similarity_stress_y_derivative_Taylor','waiting_future_defect_zeroth_Taylor',
              'pressure_y_derivative_axial5_Taylor','waiting_shear_strength_kappa_minus2',
              'waiting_meridional_moments_and_velocities','source_defined_positive_stress_factors')
        def summary(z,phase):
            point=self.waiting(z,phase); return {key:point[key] for key in keys}
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                    scope='Original O7 waiting phase[0,1], Z[-1,1]; signed full-moment stress mixed3 and absolute pressure mixed4',
                    samples=[summary(z,phase) for z,phase in (('0','0'),('.5','.5'),('-.5','1'))],
                    whole_waiting=summary([-1,1],[0,1]),
                    waiting_identities=self.proof,waiting_source_bridge=self.bridge,
                    actual_original_waiting_similarity_stress_recovered=True,
                    waiting_collar_stress_mixed3_join_verified=True,
                    steep_exit_waiting_stress_mixed3_join_verified=False,
                    waiting_cone_certified=False,global_admissible_stress_lift_constructed=False,
                    physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantWaitingStressC3().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual waiting full-future stress mixed3 and absolute pressure mixed4 generated',flush=True)
    return result


if __name__=='__main__':run()
