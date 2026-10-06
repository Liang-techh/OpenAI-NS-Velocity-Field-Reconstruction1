"""Original full raw core source units, exact integrated stress and nonzero remainder."""
import ast
import copy
import math
from types import SimpleNamespace

import sympy as s
from lei_ren_part1_paper_compliant_current_microswitch_stress_operator import (
    compile_function,expanded_stress_sectors,factored_rows_record,shifted_rows)
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST
from lei_ren_part1_paper_compliant_current_core_interior_moments import stirling,key
from lei_ren_part1_paper_compliant_microswitch_mixed_C4 import FactoredAlgebra
from lei_ren_part1_paper_compliant_frozen_comparison_field import dress
from lei_ren_part1_paper_interval_taylor import IntervalTaylor


def assemble_core_raw_rows(c,algebra,phi,V,moments,Q,p0,ratios,ratios2):
    """Direct original raw-pre units; no microscopic phase inversion."""
    phi=[dress(row,ratios) for row in phi]
    H=[dress(row,ratios) for row in moments['H']]
    K=[dress(row,ratios) for row in moments['K']]
    B=[dress(row,ratios2) for row in moments['B']]
    C=[dress(row,ratios2) for row in moments['C']]
    unit=lambda row:algebra.shift(row,(0,0,0,.5))
    square_unit=lambda row:algebra.shift(row,(0,0,0,1))
    theta=[unit(row) for row in shifted_rows(phi,c.mpf('.5'),4)]
    h=[unit(row)/2 for row in shifted_rows(H,c.mpf('.5'),4)]
    k=[unit(row)/2 for row in shifted_rows(K,c.mpf('.5'),4)]
    pp=[square_unit(row)/2 for row in shifted_rows(C,1,4)]
    bb=[square_unit(row)/2 for row in shifted_rows(B,1,4)]
    histories=dict(m=[algebra.lift(row) for row in moments['M']],h=h,k=k,
        e=[algebra.shift(moments['A'][j],(0,-1,0,0))-bb[j] for j in range(5)],p=pp)
    return dict(algebra=algebra,velocity=dict(theta=theta,axial=[algebra.lift(row) for row in V],
        radial=[algebra.lift(row) for row in shifted_rows(Q,c.mpf('.5'),4)]),
        histories=histories,absolute_pressure=[algebra.lift(p0)+pp[0]]+pp[1:],
        original_core_h_and_k_unit_half_retained=True,
        original_composite_energy_A_minus_swirl_B_retained=True,
        original_P0_separate_from_cumulative_C=True,ordinary_logR_rows_before_any_source_resolution=True)


def core_raw_source(field,Z,rho):
    c=field.ctx;mom=field.interior.evaluate(Z,rho);packet,_=field.interior.rows(Z);r=c.mpf(rho)
    phi={};vel={}
    for i in range(5):
        for k in range(6):
            original=field.interior.rebuild.profile(packet,r,radial_order=i,axial_order=k)
            phi[i,k]=original['source_profile_enclosures']['Phi']
            vel[i,k]=original['source_profile_enclosures']['Uz']
    def Euler(grids):
        return [IntervalTaylor(c,[sum((stirling(j,i)*r**i*grids[i,k] for i in range(j+1)),c.mpf(0))/math.factorial(k)
                 for k in range(6)]) for j in range(5)]
    original_phi=Euler(phi);V=Euler(vel)
    moments={name:[IntervalTaylor(c,[mom['ordinary_logR4_axial6_moment_grids'][name]['y'+str(j)+'_Z'+str(k)]/math.factorial(k)
        for k in range(6)]) for j in range(5)] for name in ('H','M','K','A','B','C')}
    Q=Euler({(i,k):mom['normalized_radial_Q_rho4_axial5'][key(i,k)] for i in range(5) for k in range(6)})
    p0=IntervalTaylor(c,mom['original_analytic_axis_pressure_axial6_Taylor_coefficients'][:6])
    ratios=mom['original_F0_relative_axial6_ordinary_derivatives']
    ratios2=mom['original_F0_squared_relative_axial6_ordinary_derivatives']
    logR=c.ln(field.core.epsilon)+c.ln(r)
    # Use the original admitted bound directly, without endpoint arithmetic choosing a source value.
    from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
    logF0=c.mpf([endpoints(-field.core.logC-field.core.Lambda*field.core.Gbar)[0],endpoints(-field.core.logC)[1]])
    logs=(c.mpf(0),2*field.physical.logP,2*logF0,c.ln(2)+logR+2*logF0-2*field.physical.logP)
    algebra=FactoredAlgebra(c,logs,[])
    raw=assemble_core_raw_rows(c,algebra,original_phi,V,moments,Q,p0,ratios,ratios2)
    return raw,dict(original_same_fixed_point_density_and_moment_packet=mom,
        original_ordinary_logR_phi_and_Uz_rows=dict(Phi=original_phi,Uz=V),
        logR=logR,formal_radius='epsilon_core*rho',formal_unit='sqrt(2*R)*F0base/Pstar',
        original_positive_epsilon_log=c.ln(field.core.epsilon),
        no_phase_width_source_or_inverse_width_in_core=True,
        fixed_basepoint_F0_bound_encloses_source_not_selects_value=True)


def core_raw_unit_theorem():
    from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z
    class Projection(FunctionJet):
        def __pow__(self,n):return self.function(self.ctx,self.expr**n,self.order)
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)) if isinstance(value,(str,int,float)) else value)
    y=s.Symbol('y',real=True);unit,Ps=s.symbols('core_unit Pstar',positive=True)
    amp=s.Function('relative_F0')(Z);J=lambda value:Projection.function(c,value,5)
    names=('Phi','V','H','M','K','A','B','C','Q')
    functions={name:s.Function('core_'+name)(y,Z) for name in names}
    rows={name:[J(s.diff(value,y,j).subs(y,0)) for j in range(5)] for name,value in functions.items()}
    p0=s.Function('normalized_P0')(Z)
    class Algebra:
        lift=staticmethod(lambda row:row if isinstance(row,Projection) else J(row))
        shift=staticmethod(lambda row,powers:Algebra.lift(row)*unit**(2*s.Rational(str(powers[3])))*Ps**(2*powers[1]))
    asts=SourceAST();env=dict(math=math,shifted_rows=shifted_rows,
        dress=lambda row,ratios:row*J(amp if ratios=='first' else amp**2))
    adapter=compile_function(asts.method('current_core_stress_operator','assemble_core_raw_rows'),env,'<original direct core units symbolic replay>')
    raw=adapter(c,Algebra(),rows['Phi'],rows['V'],{n:rows[n] for n in ('H','M','K','A','B','C')},rows['Q'],J(p0),'first','second')
    expected=dict(theta=unit*s.exp(y/2)*amp*functions['Phi'],axial=functions['V'],radial=s.exp(y/2)*functions['Q'],
        m=functions['M'],h=unit*s.exp(y/2)*amp*functions['H']/2,k=unit*s.exp(y/2)*amp*functions['K']/2,
        e=functions['A']/Ps**2-unit**2*s.exp(y)*amp**2*functions['B']/2,
        p=unit**2*s.exp(y)*amp**2*functions['C']/2,
        pressure=p0+unit**2*s.exp(y)*amp**2*functions['C']/2)
    checks={}
    for label,source in expected.items():
        actual=(raw['absolute_pressure'] if label=='pressure' else raw['velocity'][label] if label in raw['velocity'] else raw['histories'][label])
        for j,row in enumerate(actual):
            difference=s.cancel(row.expr-s.diff(source,y,j).subs(y,0))
            if difference!=0:raise ValueError('Original core raw source unit mismatch: '+label+str(j))
            for k in range(5-j):checks[label+'_y'+str(j)+'_Z'+str(k)]=s.diff(difference,Z,k)==0
    return dict(original_core_raw_unit_identities=checks,total_identities=len(checks),
        core_swirl_unit_has_no_Phi_or_hb_factor=True,
        original_h_k_half_and_composite_e_and_absolute_P0_retained=True,input_hashes=asts.hashes,passed=True)


def core_integrated_stress_theorem():
    from lei_ren_part1_paper_compliant_current_core_first_interface_check import integrated_core_equations_proof
    from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
    integrated=integrated_core_equations_proof();asts=SourceAST()
    R,z,delta=s.symbols('R Z delta',real=True,positive=True);Ps=s.Symbol('Pstar',positive=True)
    F=s.Function('actual_F')(R,z);V=s.Function('actual_V')(R,z);P0=s.Function('actual_P0')(z)
    IV,IA,IH,IK,IB,IC=[s.Function(n)(R,z) for n in ('I_V','I_V2','I_RF','I_RFV','I_RF2','I_F2')]
    densities={IV:V,IA:V**2,IH:R*F,IK:R*F*V,IB:R*F**2,IC:F**2}
    def ftc(expression):
        def selected(node):return isinstance(node,s.Derivative) and node.expr in densities and R in node.variables
        def replace(node):
            order=dict(node.variable_count)
            return s.diff(densities[node.expr],R,order.get(R,0)-1,z,order.get(z,0))
        return expression.replace(selected,replace)
    Dy=lambda value:R*s.diff(value,R)
    def rows(value):
        result=[value]
        for _ in range(4):result.append(ftc(Dy(result[-1])))
        return result
    u=rows(s.sqrt(2*R)*F/Ps);vv=rows(V)
    history=dict(m=rows(IV/R),h=rows(2*IH/(s.sqrt(2)*R**s.Rational(3,2)*Ps)),
        k=rows(2*IK/(s.sqrt(2)*R**s.Rational(3,2)*Ps)),e=rows((IA-IB)/(R*Ps**2)),p=rows(IC/Ps**2))
    pressure=rows((P0+IC)/Ps**2)
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)))
    env=dict(math=math,axial_derivative=lambda value:s.diff(value,z),product_rows=product_rows,shifted_rows=shifted_rows)
    raw=asts.replay('current_pre_pulse_stress_operator','raw_pre_stress_rows',env)(c,delta,z,u,vv,history,pressure)
    d=1-z*z;L=1-delta*z*z;M=IV/R;W=1-(1-delta)*z*M-d*s.diff(M,z)
    At=(2*L*R**2*s.diff(F,R)-W*R**2*F+(1-delta/2)*IH
        -(1-delta)*z*s.diff(IH,z)/2-d*s.diff(IK,z)+(2*delta-1)*z*IK)
    Az=(2*L*R*s.diff(V,R)-W*R*V+(1-delta)*IV/2-(1-delta)*z*s.diff(IV,z)/2
        +2*delta*z*IA-d*s.diff(IA,z)-R*(d*s.diff(P0,z)-2*(1+delta)*z*P0)
        -d*(R*s.diff(IC,z)-s.diff(IB,z))+2*(1+delta)*z*(R*IC-IB)+2*z*IB)
    checks={}
    for label,expected in (('theta',At/(L*R)),('axial',Az/(L*s.sqrt(2*R)))):
        total=sum((R**s.Rational(str(part['mode'][0]))*Ps**part['mode'][1]*part['full_derivative_rows'][0]/s.sqrt(2)
            for part in raw[label].values()),s.Integer(0))
        difference=s.cancel(ftc(total-expected))
        if difference!=0:raise ValueError('Full original core stress differs from integrated ODE: '+label)
        for j in range(4):
            for k in range(4-j):checks[label+'_y'+str(j)+'_Z'+str(k)]=s.diff(difference,R,j,z,k)==0
    return dict(full_original_raw_stress_to_integrated_core_ODE_identities=checks,
        consumed_original_integrated_core_equations_and_axis_constants=integrated,
        exact_common_source_core_stress_zero_before_physical_bounds=True,
        zero_stress_does_not_zero_any_NS_remainder=True,input_hashes=asts.hashes,passed=True)


def reduce_core_stress_from_source(rows,theorem):
    if not theorem['passed'] or not theorem['exact_common_source_core_stress_zero_before_physical_bounds']:
        raise ValueError('Exact actual core ODE theorem required')
    # Only the TOTAL vanishes. Original signed sectors are exported separately.
    # A separate aggregate is passed to the unchanged linear physical mapper.
    return {label:{'integrated_core_total_stress_zero':dict(mode=(0,0,0,0),
        full_derivative_rows=[row*0 for row in next(iter(parts.values()))['full_derivative_rows']])}
        for label,parts in rows.items()}
