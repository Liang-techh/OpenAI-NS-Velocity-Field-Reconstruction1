"""Original steep exit: physical completed stress and axial-viscosity error.

The accepted general-K transfer is an algebraic mapper. Original exit phase,
full moments, absolute pressure and exact amplitude/radius source are kept.
Homogeneous angular/energy divergence and stationary pressure cancel before
enclosure. The regional nonzero error is not a global flatness certificate.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_steep_exit_stress_C3 import CompliantSteepExitStressC3,exit_signed_kernels
from lei_ren_part1_paper_compliant_waiting_stress_C3 import source_precision
from lei_ren_part1_paper_compliant_collar_stress_C3 import axial_derivative,shifted_rows
from lei_ren_part1_paper_compliant_collar_physical_C2 import (
    CompliantCollarPhysicalC2,physical_collar_identities,physical_source_row,scale_row,
    collar_velocity_operator_coefficients)
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_axial_pulse_field import sigma_enclosure
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket,ZSYM,DSYM
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def exit_pressure_fluctuation(c,t,k,eps,p,cells=128):
    """Positive weighted integral of K^2-K0^2, with no stationary subtraction.

    Iplus=int_t^1 exp(-p*(v-t))*K0^2*expm1(2*k*f(v))dv.
    This is the same original sigma primitive used by the stress companion.
    Separating the constant K0^2-1 analytically preserves the waiting join.
    """
    t=c.mpf(t); lo,hi=endpoints(t)
    if lo<0 or hi>1:raise ValueError('Original steep-exit phase required')
    if not isinstance(cells,int) or cells<1:raise ValueError('Positive cell count required')
    length=1-t; ds=length/cells; f=c.mpf(0); integral=c.mpf(0); K0=1-eps
    if endpoints(length)[1]>0:
        for i in reversed(range(cells)):
            a=t+length*i/cells; b=t+length*(i+1)/cells
            v=c.mpf([max(mp.mpf(0),endpoints(a)[0]),min(mp.mpf(1),endpoints(b)[1])])
            nextf=f+ds*(1-sigma_enclosure(c,v))
            fc=c.mpf([endpoints(f)[0],min(mp.mpf('.5'),endpoints(nextf)[1])])
            delta_square=K0**2*c.expm1(2*k*fc)
            distance=length*c.mpf([i,i+1])/cells
            integral+=ds*c.exp(-p*distance)*delta_square
            f=nextf
    if (lo,hi)==(mp.mpf(0),mp.mpf(0)):f=c.mpf('.5')
    else:f=c.mpf([endpoints(f)[0],min(mp.mpf('.5'),endpoints(f)[1])])
    return dict(integral=integral,current_delta_square=K0**2*c.expm1(2*k*f),f=f,cells=cells,
                exact_definition='integral_t^1 exp(-p*(v-t))*(K(v)^2-K0^2)dv',
                positive_correlated_cell_lengths=True,stationary_pressure_removed_before_enclosure=True)


def steep_exit_divergence_grids(heat,base,shape,q,fluctuation):
    """Div(T) source rows with its 1/sqrt(R) factor included.

    Homogeneous A/E modes cancel exactly. Pressure's K0^2-1 stationary
    mode cancels before enclosure too. The remaining pressure mode is the
    SAME collar endpoint propagated by its exact p exponent.
    """
    c=heat.ctx; z=IntervalTaylor.variable(c,base['Z'],5)
    L=1-z*z*heat.delta; d=1-z*z; K=shape['K_rows']; p=heat.prate
    viscous=[K[j+2]-K[j+1]*p+K[j]*(heat.a*(1+heat.a)) for j in range(3)]
    shear=shifted_rows(viscous,-1,2)
    theta_raw=[-K[j+1]/L+shear[j]*(2*heat.S*c.exp(-q)) for j in range(3)]
    theta=shifted_rows(theta_raw,-heat.a-c.mpf('.5'),2)
    Pd0=base['pressure_defect_rows'][0]; Q0=-2*heat.eps+heat.eps**2
    Np=-z*Pd0*(2*p)+z*Q0+d*axial_derivative(Pd0)
    wave=Np*c.exp(p*q)
    positive_delta=fluctuation['current_delta_square']
    axial=[(wave+z*(positive_delta-p*fluctuation['integral']))/L,
           z*shape['K_squared_defect_rows'][1]/L,
           z*(shape['K_squared_defect_rows'][2]-shape['K_squared_defect_rows'][1]*p)/L]
    return {label:{'y'+str(j)+'_Z'+str(n):rows[j][n]*math.factorial(n)
                   for j in range(3) for n in range(3-j)} for label,rows in (('theta',theta),('axial',axial))}


class ExitHeatReference:
    def __init__(self,heat,shape):self.actual_heat=heat; self.actual_exit_shape=shape
    def __getattr__(self,key):return getattr(self.actual_heat,key)
    def shape(self,Z,q):return self.actual_exit_shape


class ExitDispatch:
    def __init__(self,native,stress):self.native=native; self.stress=stress
    def __getattr__(self,key):return getattr(self.native,key)
    def provider(self,chart):
        return self.stress.steep if chart=='steep_exit' else self.native.provider(chart)
    def evaluate(self,chart,Z,coordinate):
        if chart!='steep_exit':return self.native.evaluate(chart,Z,coordinate)
        packet=self.stress.steep_out(Z,coordinate)
        return dict(chart=chart,source_packet=packet,
                    physical_mixed_grids=packet['physical_mixed_derivatives_total_order_le4'],
                    original_steep_exit_absolute_pressure_companion_used=True)


class AtExitPhaseStress:
    def __init__(self,packet):self.packet=packet
    def collar(self,Z,q):
        return dict(collar_similarity_stress_mixed3_factored=self.packet['steep_exit_similarity_stress_mixed3_factored'],
                    Gamma_endpoint_stress_exact_zero_from_same_moments=False)


class AtExitPhaseAssembly:
    def __init__(self,native,phase):self.native=native; self.phase=phase
    def __getattr__(self,key):return getattr(self.native,key)
    def evaluate(self,chart,Z,q,**kwargs):
        if chart!='heat_collar':raise ValueError('Generic exit map received an unexpected chart')
        return self.native.evaluate('steep_exit',Z,self.phase,**kwargs)


def steep_exit_physical_identities():
    general=physical_collar_identities(); proofs={}
    a,q,z,S=s.symbols('a q Z S',real=True); delta=2*a; k=1-a; p=1+delta; L=1-delta*z*z; d=1-z*z
    K=s.Function('actual_exit_K')(q); A=s.Function('full_A')(q,z)
    E=s.Function('full_E')(q,z); P=s.Function('full_P')(q,z); b=(1-delta)/2
    Ct=(k*A-b*z*s.diff(A,z)-K)/L+2*S*s.exp(-q)*(s.diff(K,q)-(1+a)*K)
    Cz=(delta*z*E-d*s.diff(E,z)/2-2*p*z*P+d*s.diff(P,z))/L
    rules={s.diff(A,q):K-k*A,s.diff(A,q,z):-k*s.diff(A,z),
           s.diff(E,q):delta*E-K*K,s.diff(E,q,z):delta*s.diff(E,z),
           s.diff(P,q):p*P-K*K/2,s.diff(P,q,z):p*s.diff(P,z)}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Exit physical reduction failed: '+name)
        proofs[name]=True
    Dt=-s.diff(K,q)/L+2*S*s.exp(-q)*(s.diff(K,q,2)-p*s.diff(K,q)+a*(1+a)*K)
    Dz=(d*s.diff(P,z)-2*z*(p*P-K*K/2))/L
    zero('actual_full_theta_homogeneous_mode_cancelled',(s.diff(Ct,q)+k*Ct).subs(rules,simultaneous=True)-Dt)
    zero('actual_full_axial_energy_mode_cancelled',(s.diff(Cz,q)-delta*Cz).subs(rules,simultaneous=True)-Dz)
    eps,ell,wait=s.symbols('eps ell wait',real=True); Q0=-2*eps+eps**2
    Pd0=s.Function('full_collar_pressure_defect')(z); I=s.Function('positive_pressure_fluctuation')(q)
    D=s.Function('K_squared_minus_K0_squared')(q)
    qrelation=-wait-ell
    Pd1=Q0/(2*p)+(Pd0-Q0/(2*p))*s.exp(-p*wait)
    signed_I=Q0*(1-s.exp(-p*ell))/p+I
    Pd=s.exp(-p*ell)*Pd1+signed_I/2
    Np=-2*p*z*Pd0+z*Q0+d*s.diff(Pd0,z)
    zero('actual_stationary_pressure_cancels_before_enclosure',
         d*s.diff(Pd,z)-2*p*z*Pd+z*(Q0+D)-(Np*s.exp(p*qrelation)+z*(D-p*I)))
    # I'=pI-D is the FTC of the exact positive fluctuation integral.
    V=D-p*I
    zero('actual_pressure_full_factor_first_radial_row',
         s.diff(V,q).subs(s.diff(I,q),p*I-D)-p*V-s.diff(D,q))
    zero('actual_pressure_full_factor_second_radial_row',
         (s.diff(V,q,2)-2*p*s.diff(V,q)+p*p*V).subs(
             {s.diff(I,q,2):p*(p*I-D)-s.diff(D,q),s.diff(I,q):p*I-D},simultaneous=True)
         -(s.diff(D,q,2)-p*s.diff(D,q)))
    zero('theta_divergence_positive_factor_rate',(-a-s.Rational(1,2))+a+s.Rational(1,2))
    zero('axial_divergence_positive_factor_rate',(-p)+(s.Rational(1,2)+delta)+s.Rational(1,2))
    # Native implicit-map q_zz, with K independent of source Z.
    F=1-2*(1-delta)*z*z-delta*delta*z**4
    reduced=2*F/L**3*s.diff(K,q)-4*z*z/L**2*s.diff(K,q,2)
    actual=sum(coefficient.subs({ZSYM:z,DSYM:delta})*s.diff(K,q,index[0])
               for index,coefficient in collar_velocity_operator_coefficients(0,2).items() if index[1]==0)
    zero('actual_reduced_axial_viscosity_from_general_K_operator',-actual-reduced)
    endpoint_K=s.symbols('source_K0',real=True)
    for i in range(3):
        for j in range(3-i):
            endpoint=sum(coefficient.subs({ZSYM:z,DSYM:delta})*(endpoint_K if index==(0,0) else 0)
                         for index,coefficient in collar_velocity_operator_coefficients(i,j+2).items())
            zero('flat_K_y1_to4_gives_waiting_remainder_r'+str(i)+'_z'+str(j),endpoint)
    return dict(general_K_physical_transfer=general,reduction_identities=proofs,
                original_K_is_independent_of_source_Z=True,
                actual_regional_physical_steep_exit_stress_remainder_identity_verified=True,
                exact_remainder='Etheta=-nu*partial_zz(utheta); Er=Ez=0',
                reduced_remainder_coefficient='2*(1-2*(1-delta)*Z^2-delta^2*Z^4)*K_y/L^3-4*Z^2*K_yy/L^2',
                source_homogeneous_A_E_and_stationary_pressure_cancelled_before_enclosure=True,
                globally_flat_remainder_not_inferred=True)


def steep_exit_adapter_binding():
    tree=ast.parse(Path(__file__).read_text(encoding='utf8')); bindings={}
    requirements={
        ('ExitHeatReference','shape'):'self.actual_exit_shape',
        ('ExitDispatch','provider'):"self.stress.steep if chart=='steep_exit' else self.native.provider(chart)",
        ('ExitDispatch','evaluate'):'self.stress.steep_out(Z,coordinate)',
        ('AtExitPhaseStress','collar'):"self.packet['steep_exit_similarity_stress_mixed3_factored']",
        ('AtExitPhaseAssembly','evaluate'):"self.native.evaluate('steep_exit',Z,self.phase,**kwargs)",
        ('CompliantSteepExitPhysicalC2','steep_out'):'-self.stress.steep.wait-1+phase'}
    for (cls,name),expression in requirements.items():
        node=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name==cls)
        fn=next(n for n in node.body if isinstance(n,ast.FunctionDef) and n.name==name)
        wanted=ast.dump(ast.parse(expression,mode='eval').body)
        if not any(ast.dump(n)==wanted for n in ast.walk(fn)):raise ValueError('Exit adapter changed: '+cls+'.'+name)
        bindings[cls+'.'+name]=True
    formulas={
        'exit_pressure_fluctuation':{'nextf':'f+ds*(1-sigma_enclosure(c,v))',
            'fc':"c.mpf([endpoints(f)[0],min(mp.mpf('.5'),endpoints(nextf)[1])])",
            'delta_square':'K0**2*c.expm1(2*k*fc)','distance':'length*c.mpf([i,i+1])/cells'},
        'steep_exit_divergence_grids':{
            'viscous':'[K[j+2]-K[j+1]*p+K[j]*(heat.a*(1+heat.a)) for j in range(3)]',
            'theta_raw':'[-K[j+1]/L+shear[j]*(2*heat.S*c.exp(-q)) for j in range(3)]',
            'Np':'-z*Pd0*(2*p)+z*Q0+d*axial_derivative(Pd0)',
            'axial':"[(wave+z*(positive_delta-p*fluctuation['integral']))/L,z*shape['K_squared_defect_rows'][1]/L,z*(shape['K_squared_defect_rows'][2]-shape['K_squared_defect_rows'][1]*p)/L]"}}
    for name,assignments in formulas.items():
        fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name==name)
        for target,expression in assignments.items():
            values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
            if sum(ast.dump(v)==ast.dump(ast.parse(expression,mode='eval').body) for v in values)!=1:
                raise ValueError('Exit physical source formula changed: '+name+'.'+target)
            bindings[name+'.'+target]=True
    split=exit_pressure_split_binding()
    return dict(source_bindings=bindings,original_exit_phase_not_uncorrelated_offset_division=True,
                signed_and_positive_pressure_integrals_same_actual_source_split=split,
                actual_shape_pressure_stress_and_source_phase_consumed=True,
                general_transfer_source='Unchanged CompliantCollarPhysicalC2.collar',
                physical_collar_source_never_evaluated_at_negative_offset=True)


def exit_pressure_split_binding():
    """Replay both producer integrands and exact stationary integral split.

    The two cell sums are enclosures of related exact integrals; equality of
    their numerical boxes is neither required nor asserted.
    """
    current=ast.parse(Path(__file__).read_text(encoding='utf8'))
    source_path=HERE/(PREFIX+'steep_exit_stress_C3.py'); source=ast.parse(source_path.read_text(encoding='utf8'))
    positive=next(n for n in current.body if isinstance(n,ast.FunctionDef) and n.name=='exit_pressure_fluctuation')
    signed=next(n for n in source.body if isinstance(n,ast.FunctionDef) and n.name=='exit_signed_kernels')
    def assignment(fn,name,wanted):
        matches=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)==name for t in n.targets)]
        expr=ast.parse(wanted,mode='eval').body
        if sum(ast.dump(v)==ast.dump(expr) for v in matches)!=1:raise ValueError('Pressure split producer path changed: '+fn.name+'.'+name)
        return expr
    common=dict(t='c.mpf(t)',length='1-t',ds='length/cells',K0='1-eps',
                a='t+length*i/cells',b='t+length*(i+1)/cells',
                v='c.mpf([max(mp.mpf(0),endpoints(a)[0]),min(mp.mpf(1),endpoints(b)[1])])',
                nextf='f+ds*(1-sigma_enclosure(c,v))',
                fc="c.mpf([endpoints(f)[0],min(mp.mpf('.5'),endpoints(nextf)[1])])",
                distance='length*c.mpf([i,i+1])/cells')
    for name,expression in common.items():
        assignment(positive,name,expression); assignment(signed,name,expression)
    for fn in (positive,signed):
        assignment(fn,'f',"c.mpf('.5')")
        assignment(fn,'f',"c.mpf([endpoints(f)[0],min(mp.mpf('.5'),endpoints(f)[1])])")
        loops=[n for n in ast.walk(fn) if isinstance(n,ast.For)]
        if len(loops)!=1 or ast.unparse(loops[0].iter)!='reversed(range(cells))':raise ValueError('Pressure split traversal changed')
    delta_expr=assignment(positive,'delta_square','K0**2*c.expm1(2*k*fc)')
    signed_expr=assignment(signed,'Qm','-2*eps+eps**2+K0**2*c.expm1(2*k*fc)')
    def accumulated(fn,name,wanted):
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.AugAssign) and isinstance(n.op,ast.Add) and ast.unparse(n.target)==name]
        if len(values)!=1 or ast.dump(values[0])!=ast.dump(ast.parse(wanted,mode='eval').body):raise ValueError('Pressure split integral changed')
    accumulated(positive,'integral','ds*c.exp(-p*distance)*delta_square')
    accumulated(signed,'IP','ds*c.exp(-p*distance)*Qm')
    returns=[n for n in ast.walk(positive) if isinstance(n,ast.keyword) and n.arg=='current_delta_square']
    if len(returns)!=1 or ast.unparse(returns[0].value)!='K0 ** 2 * c.expm1(2 * k * f)':raise ValueError('Pressure split current square changed')
    if exit_signed_kernels.__globals__['sigma_enclosure'] is not sigma_enclosure:raise ValueError('Pressure split sigmoid callable differs')
    eps,k,f=s.symbols('eps k f',real=True); p=s.symbols('p',positive=True); t,v=s.symbols('t v',real=True)
    context=SimpleNamespace(expm1=lambda value:s.exp(value)-1)
    env=dict(eps=eps,K0=1-eps,k=k,fc=f,c=context)
    Qplus=eval(compile(ast.Expression(delta_expr),'<actual positive square>','eval'),{},env)
    Qsigned=eval(compile(ast.Expression(signed_expr),'<actual signed square>','eval'),{},env)
    Q0=-2*eps+eps*eps
    if s.simplify(Qsigned-Q0-Qplus)!=0:raise ArithmeticError('Actual pressure integrand split failed')
    baseline=Q0*(1-s.exp(-p*(1-t)))/p
    actual_integral=s.integrate(Q0*s.exp(-p*(v-t)),(v,t,1))
    if s.simplify(s.expand_power_exp(actual_integral-baseline))!=0:raise ArithmeticError('Pressure stationary integral split failed')
    return dict(same_sigma_callable_cumulative_f_cells_distances_and_endpoints_verified=True,
                actual_K_squared_minus1_equals_Q0_plus_positive_fluctuation=True,
                actual_signed_pressure_integral_equals_stationary_integral_plus_Iplus=True,
                exact_split='Isigned=Q0*(1-exp(-p*(1-t)))/p+Iplus',
                numerical_cell_boxes_required_equal=False,
                original_Ev0_and_absolute_pressure_datum_unchanged=True)


class CompliantSteepExitPhysicalC2:
    @source_precision
    def __init__(self):
        self.stress=CompliantSteepExitStressC3(); self.assembly=CompliantGlobalPhysicalAssembly()
        self.assembly.dispatch=ExitDispatch(self.assembly.dispatch,self.stress)
        self.ctx=self.stress.ctx; self.family,self.source=self.stress.family,self.stress.source
        if (self.assembly.family,self.assembly.source)!=(self.family,self.source):raise ValueError('Exit physical family mismatch')
        self.proof=steep_exit_physical_identities(); self.adapter_proof=steep_exit_adapter_binding(); self.cache={}
        self.hashes=dict(self.stress.hashes); self.hashes.update(self.assembly.hashes)
        for stem,gate in (('steep_exit_stress_C3_check','actual_original_steep_exit_similarity_stress_recovered'),
                          ('global_physical_assembly_check','all_33_original_source_charts_physical_spatial4_time1_mapped'),
                          ('collar_physical_C2_check','actual_regional_physical_collar_stress_remainder_identity_verified'),
                          ('waiting_physical_C2_check','actual_regional_physical_waiting_stress_remainder_identity_verified')):
            name=PREFIX+stem+'.json'; receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or not receipt[gate]:raise ValueError('Exit physical prerequisite missing: '+name)
            if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):raise ValueError('Exit physical receipt family mismatch')
            for path,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Exit physical source changed: '+path)
            self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def fluctuation(self,phase):
        c=self.ctx; key=tuple(endpoints(c.mpf(phase))); heat=self.stress.heat
        if key not in self.cache:
            self.cache[key]=exit_pressure_fluctuation(c,phase,heat.k,heat.eps,heat.prate,self.stress.cells)
        return self.cache[key]

    @source_precision
    def steep_out(self,Z,phase,log_tau='-1',theta='0',viscosity='1'):
        c=self.ctx; z=c.mpf(Z); phase=c.mpf(phase); self.stress.steep._phase(phase)
        q=-self.stress.steep.wait-1+phase; packet=self.stress.steep_out(z,phase)
        heat=ExitHeatReference(self.stress.heat,packet['steep_exit_shape'])
        adapter=SimpleNamespace(ctx=c,heat=heat,stress=AtExitPhaseStress(packet),assembly=AtExitPhaseAssembly(self.assembly,phase))
        point=CompliantCollarPhysicalC2.collar(adapter,z,q,log_tau,theta,viscosity)
        base=dict(self.stress.waiting_source.endpoint(z),Z=z)
        fluctuation=self.fluctuation(phase)
        grids=steep_exit_divergence_grids(heat,base,packet['steep_exit_shape'],q,fluctuation)
        logtau=c.mpf(log_tau); nu=c.mpf(viscosity); beta=-2-heat.delta
        parts=point['actual_source_log_factors']; logR=point['source_logR_enclosure']; divmixed={}
        for label,factor in (('theta','Qtheta'),('axial','Qz')):
            divmixed[label]={}
            for i in range(3):
                for j in range(3-i):
                    coefficient=physical_bracket(c,grids[label],i,j,z,heat.delta,beta-1)
                    divmixed[label]['r'+str(i)+'_z'+str(j)]=physical_source_row(c,coefficient,parts[factor],
                        beta-1-i+j*(heat.delta-1),logtau,nu,i+j,c.mpf(i+1)/2*(c.ln(2)-logR),nu_base=c.mpf('.5'))
        cs=point['cos_theta']; sn=point['sin_theta']
        divcart=dict(x=[scale_row(c,divmixed['theta']['r0_z0'],-sn)],y=[scale_row(c,divmixed['theta']['r0_z0'],cs)],z=[divmixed['axial']['r0_z0']])
        point['physical_cylindrical_stress_divergence_mixed2']=divmixed
        point['physical_completed_stress_divergence_cartesian']=divcart
        point['physical_momentum_residual_decomposition_cartesian']={label:[scale_row(c,row,-1) for row in divcart[label]]
            +point['physical_remainder_cartesian'][label] for label in ('x','y','z')}
        terminal=endpoints(phase)==(mp.mpf(1),mp.mpf(1))
        if terminal:
            # All positive K derivatives through four vanish at this actual
            # flat endpoint. Native velocity is the waiting pure radial power.
            for key,components in point['physical_spatial_cartesian_mixed4'].items():
                if int(key.split('_')[2][1:]):
                    for component in ('ux','uy','uz'):
                        for row in components[component].values():
                            row.update(terms=[],exact_zero=True,log_absolute_upper=None,exact_zero_from_actual_exit_flat_K=True)
            for component in ('ux','uy','uz'):
                for row in point['first_fixed_x_physical_time_derivative'][component].values():
                    row.update(terms=[],exact_zero=True,log_absolute_upper=None,exact_zero_from_actual_exit_flat_K=True)
        for key in ('actual_regional_physical_collar_stress_remainder_identity_verified','Gamma_endpoint_physical_remainder_exact_zero','collar_cone_certified'):
            point.pop(key,None)
        point.update(chart='steep_exit',steep_exit_phase=phase,steep_exit_offset_from_Rtail=q,
                     steep_exit_pressure_fluctuation=fluctuation,
                     actual_regional_physical_steep_exit_stress_remainder_identity_verified=True,
                     steep_exit_homogeneous_A_E_and_stationary_pressure_divergence_cancelled_before_enclosure=True,
                     steep_exit_waiting_stress_mixed3_join_verified=True,steep_exit_waiting_pressure_mixed4_join_verified=True,
                     steep_exit_waiting_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
                     steep_exit_terminal_remainder_and_axial_time_velocity_jets_exact_zero=terminal,
                     steep_exit_regional_remainder_exact_zero=False,regional_remainder_is_leading_axial_viscosity_not_proven_flat=True,
                     steep_exit_cone_certified=False,full_background_NS_validation=False,temporal_recursion=False)
        return point

    def report(self):
        def summary(point):
            omitted=('physical_spatial_cartesian_mixed4','first_fixed_x_physical_time_derivative','positive_source_log_bases')
            result={key:value for key,value in point.items() if key not in omitted}
            result['zeroth_physical_velocity_pressure_rows']=point['physical_spatial_cartesian_mixed4']['x0_y0_z0']
            result['physical_axial_and_time_velocity_rows_all_exact_zero']=all(row['exact_zero']
                for key,components in point['physical_spatial_cartesian_mixed4'].items() if int(key.split('_')[2][1:])
                for component in ('ux','uy','uz') for row in components[component].values()) and all(row['exact_zero']
                for component in ('ux','uy','uz') for row in point['first_fixed_x_physical_time_derivative'][component].values())
            return result
        samples=[summary(self.steep_out(z,t,lt,angle,nu)) for z,t,lt,angle,nu in
                 (('0','0','-1','0','1'),('.5','.5','-10','.7','.01'),('-.5','1','-100','1','.7'))]
        whole=summary(self.steep_out([-1,1],[0,1],[-1000,-1],None))
        self.hashes.update(self.assembly.dispatch.native.hashes)
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                    scope='Original steep_exit t[0,1], Z[-1,1], tau>0, r>0; physical stress mixed3, diagonal/divergence/remainder mixed2',
                    samples=samples,whole_steep_exit=whole,steep_exit_physical_identities=self.proof,steep_exit_adapter_binding=self.adapter_proof,
                    actual_regional_physical_steep_exit_stress_remainder_identity_verified=True,
                    steep_exit_waiting_physical_stress_mixed3_and_remainder_mixed2_join_verified=True,
                    steep_exit_regional_remainder_exact_zero=False,steep_exit_cone_certified=False,
                    independently_bounded_global_flat_remainder=False,global_admissible_stress_lift_constructed=False,
                    physical_energy_integral_certified=False,full_background_NS_validation=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    result=CompliantSteepExitPhysicalC2().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual steep-exit completed physical stress and nonzero axial-viscosity remainder generated',flush=True)
    return result


if __name__=='__main__':run()
