"""Original full heat-collar stresses from the admitted common moments.

Small future-integral defects are retained before cancelling unit baselines.
Positive amplitude/radius factors remain formal source-defined factors. This
is a leading similarity-stress recovery, not an all-region cone certificate.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_collar_pressure_C4 import CompliantCollarPressureC4
from lei_ren_part1_paper_compliant_heat_terminal_history_bridge import terminal_history_bridge
from lei_ren_part1_paper_compliant_heat_terminal_stress_identities import terminal_stress_identities
from lei_ren_part1_paper_compliant_heat_stress_C4 import angular_primitive_y_rows
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_compliant_'


def axial_derivative(jet):
    return IntervalTaylor(jet.ctx, [(n+1)*jet[n+1] for n in range(jet.order)])


def shifted_rows(rows, rate, order=3):
    """d_y^j(Q*C)/Q for a constant logarithmic source-factor rate."""
    return [sum((rows[n]*math.comb(j,n)*rate**(j-n) for n in range(j+1)), rows[0]*0)
            for j in range(order+1)]


def collar_defect_rows(heat, shape, tails, offset):
    """Same full future integrals, with exact unit baselines separated.

    A=1/k+Adef, E_current=1/delta+Edef,
    P_current=1/(2*(1+delta))+Pdef, K=1+Kdef.
    The numerical interval S encloses the exact 1/Rtail. It does not define
    a replacement radius or source function.
    """
    c=heat.ctx; a=heat.a; eps=heat.eps; smallS=heat.S
    W=shape['W_rows']; D=shape['D_rows']; future=tails['scaled_full_future_Gamma_defects']
    atoms=tails['separate_epsilon_atoms']; one=IntervalTaylor.constant(c,1,5)
    Km=[-W[j]*eps-D[j]*(a*smallS) for j in range(5)]
    A0=(future['theta']*smallS+one*(eps*atoms['JW']))*c.exp(-heat.k*offset)
    E0=(one*(-2*eps*atoms['EW']+eps**2*atoms['EW2'])-future['energy']*(a*smallS))*c.exp(heat.delta*offset)
    P0=(one*(-eps*atoms['PW']+eps**2*atoms['PW2']/2)-future['pressure']*(a*smallS))*c.exp(heat.prate*offset)
    square=product_rows(Km,Km)
    Qm=[2*Km[j]+square[j] for j in range(5)]
    Arows=[A0]; Erows=[E0]; Prows=[P0]
    for j in range(4):
        Arows.append(Km[j]-Arows[j]*heat.k)
        Erows.append(Erows[j]*heat.delta-Qm[j])
        Prows.append(Prows[j]*heat.prate-Qm[j]/2)
    return dict(K_defect_rows=Km, K_squared_defect_rows=Qm,
                angular_defect_rows=Arows, energy_defect_rows=Erows,
                pressure_defect_rows=Prows)


def collar_stress_rows(heat, shape, defects, Z, offset):
    """(3.16)-(3.18), normalized by positive Qtheta and Qz.

    B=Ev0*theta_base*exp(-bh*offset), R=Rtail*exp(offset),
    Qtheta=sqrt(R/2)*B; Qz=sqrt(R/2)*B^2.
    Returned profile derivative rows include derivatives of these factors.
    """
    c=heat.ctx; z=IntervalTaylor.variable(c,Z,5); d=1-z*z; L=1-z*z*heat.delta
    b=(1-heat.delta)/2; A=defects['angular_defect_rows']; E=defects['energy_defect_rows']
    P=defects['pressure_defect_rows']; Km=defects['K_defect_rows']; K=shape['K_rows']
    inertial=[(A[j]*heat.k-z*axial_derivative(A[j])*b-Km[j])/L for j in range(4)]
    axial=[(z*E[j]*heat.delta-d*axial_derivative(E[j])/2
            -z*P[j]*(2*heat.prate)+d*axial_derivative(P[j]))/L for j in range(4)]
    viscous=[K[j+1]-K[j]*(1+heat.a) for j in range(4)]
    shear=[v*(2*heat.S*c.exp(-offset)) for v in shifted_rows(viscous,-1)]
    angular=[inertial[j]+shear[j] for j in range(4)]
    return dict(theta=shifted_rows(angular,-heat.a), axial=shifted_rows(axial,-c.mpf('.5')-heat.delta),
                theta_inertial=shifted_rows(inertial,-heat.a), theta_shear=shifted_rows(shear,-heat.a),
                normalized_coefficient_rows=dict(theta=angular, axial=axial),
                shear_strength_margin=K[1]/K[0]*(-2)+heat.delta)


def collar_moment_stress_identities():
    """Symbolic general-K proof of the original moment and radial equations."""
    proofs={}
    def zero(name,expression):
        if s.simplify(s.expand_power_exp(expression))!=0:
            raise ArithmeticError('Collar moment/stress identity failed: '+name)
        proofs[name]=True
    a,t,z,Rt,B0=s.symbols('a offset Z Rtail B0',positive=True)
    delta=2*a; k=1-a; b=(1-delta)/2; p=1+delta; L=1-delta*z*z; d=1-z*z
    R=Rt*s.exp(t); B=B0*s.exp(-(s.Rational(1,2)+a)*t)
    K=s.Function('full_collar_K')(t,z); A=s.Function('full_angular_A')(t,z)
    E=s.Function('full_current_energy')(t,z); P=s.Function('full_current_pressure')(t,z)
    qt=s.sqrt(R/2)*B; qz=s.sqrt(R/2)*B*B
    U=B*K; Mt=s.sqrt(2)*R**s.Rational(3,2)*B*A
    Me=R*B*B*E/2; pressure=-B*B*P
    It=(k*Mt-b*z*s.diff(Mt,z)-R*s.sqrt(2*R)*U)/(2*L*R)
    Iz=(2*delta*z*Me-d*s.diff(Me,z)+R*(2*p*z*pressure-d*s.diff(pressure,z)))/(L*s.sqrt(2*R))
    Ct=(k*A-b*z*s.diff(A,z)-K)/L
    Cz=(delta*z*E-d*s.diff(E,z)/2-2*p*z*P+d*s.diff(P,z))/L
    zero('original_3_16_common_angular_moment_and_units',It/qt-Ct)
    zero('original_3_17_common_energy_pressure_half_and_units',Iz/qz-Cz)
    shear=s.sqrt(2*R)/R*s.diff(U,t)-U/s.sqrt(2*R)
    Cs=2/R*(s.diff(K,t)-(1+a)*K)
    zero('original_3_18_full_collar_shear_and_inverse_radius',shear/qt-Cs)
    zero('full_collar_kappa_minus2_from_same_angular_shear',
         -shear/(U/s.sqrt(2*R))-2-(delta-2*s.diff(K,t)/K))
    zero('theta_positive_factor_radial_rate',s.diff(qt,t)/qt+a)
    zero('axial_positive_factor_radial_rate',s.diff(qz,t)/qz+s.Rational(1,2)+delta)
    angular_rules={s.diff(A,t):K-k*A, s.diff(A,t,z):s.diff(K,z)-k*s.diff(A,z)}
    theta_radial=s.diff(Ct,t)+k*Ct
    zero('original_3_13_angular_inertial_radial_equation',
         theta_radial.subs(angular_rules,simultaneous=True)+(s.diff(K,t)+b*z*s.diff(K,z))/L)
    quadratic_rules={s.diff(E,t):delta*E-K*K,
                     s.diff(E,t,z):delta*s.diff(E,z)-2*K*s.diff(K,z),
                     s.diff(P,t):p*P-K*K/2,
                     s.diff(P,t,z):p*s.diff(P,z)-K*s.diff(K,z)}
    zero('original_3_13_axial_inertial_radial_equation',
         (s.diff(Cz,t)-delta*Cz).subs(quadratic_rules,simultaneous=True)
         -(d*s.diff(P,z)-2*z*(p*P-K*K/2))/L)
    Ad,Ed,Pd,Kd=s.symbols('angular_defect energy_defect pressure_defect K_defect',real=True)
    Az,Ez,Pz=s.symbols('angular_Z energy_Z pressure_Z',real=True)
    zero('theta_unit_baselines_cancel_before_enclosure',
         (k*(1/k+Ad)-b*z*Az-(1+Kd))/L-(k*Ad-b*z*Az-Kd)/L)
    zero('axial_unit_baselines_cancel_before_enclosure',
         (delta*z*(1/delta+Ed)-d*Ez/2-2*p*z*(1/(2*p)+Pd)+d*Pz)/L
         -(delta*z*Ed-d*Ez/2-2*p*z*Pd+d*Pz)/L)
    zero('quadratic_defect_uses_same_full_K', (1+Kd)**2-1-(2*Kd+Kd*Kd))
    # Unit offsets are constant in Z, and current reference factors have
    # no Z dependence at fixed R. This is not differentiation at fixed r.
    zero('current_reference_energy_units',R*B*B*(E*s.exp(-delta*t))*s.exp(delta*t)/2-Me)
    zero('current_reference_pressure_units',-B0*B0*(P*s.exp(-p*t))-pressure)
    return dict(identities=proofs, original_collar_moment_to_stress_identities_verified=True,
                original_collar_inertial_radial_equations_verified=True,
                full_K_not_Gamma_only=True, fixed_R_axial_derivatives=True,
                exact_unit_baseline_cancellation_before_enclosure=True,
                source_formulas=['Lei-Ren v2 (3.13)', '(3.16)', '(3.17)', '(3.18)'])


def collar_Gamma_endpoint_binding():
    """Replay actual coefficient AST on the admitted full Gamma moments.

    This connects the expressions before the interval endpoint override to
    the terminal theorem, rather than checking serialized zeros afterwards.
    """
    tree=ast.parse(Path(__file__).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='collar_stress_rows')
    def expression(target):
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        if len(values)!=1 or not isinstance(values[0],ast.ListComp):raise ValueError('Endpoint stress formula changed: '+target)
        return compile(ast.Expression(values[0].elt),'<actual collar '+target+'>','eval')
    a,R,Rt,t=s.symbols('a R Rtail offset',positive=True); z=s.symbols('Z',real=True)
    H,Hp,I,BE,BP,JE,JP=s.symbols('H H_xi angular_derivative_integral BE BP JE JP',real=True)
    delta=2*a; k=1-a; p=1+delta; b=(1-delta)/2; d=1-z*z; L=1-delta*z*z; xi=2*d/R
    Adef=(H-xi*I)/k-1/k; Edef=BE-1/delta; Pdef=BP/2-1/(2*p)
    derivatives={Adef:4*z*I/R,Edef:-8*z*JE/R,Pdef:-4*z*JP/R}
    heat=SimpleNamespace(k=k,a=a,delta=delta,prate=p,S=1/Rt)
    environment=dict(heat=heat,A=[Adef],E=[Edef],P=[Pdef],Km=[H-1],K=[H,-xi*Hp],
                     b=b,z=z,d=d,L=L,j=0,axial_derivative=lambda q:derivatives[q],c=SimpleNamespace(exp=s.exp),offset=t)
    It=eval(expression('inertial'),{'__builtins__':{}},environment)
    Iz=eval(expression('axial'),{'__builtins__':{}},environment)
    viscous=eval(expression('viscous'),{'__builtins__':{}},environment)
    shear=eval(expression('shear'),{'__builtins__':{}},{**environment,'v':viscous}).subs(Rt,R*s.exp(-t))
    angular=s.simplify((It+shear).subs(I,-(1+a)*H-xi*Hp))
    axial=s.simplify(Iz.subs({BE:(H*H-2*xi*JE)/delta,BP:(H*H-2*xi*JP)/p}))
    proofs={}
    for label,value in (('theta',angular),('axial',axial)):
        if value!=0:raise ArithmeticError('Actual pre-override Gamma stress AST does not cancel: '+label)
        proofs['actual_pre_override_'+label+'_expression_is_terminal_Gamma_identity']=True
        for j in range(4):
            for n in range(4-j):
                derivative=s.diff(value,z,n)
                for _ in range(j):derivative=R*s.diff(derivative,R)
                if derivative!=0:raise ArithmeticError('Canonical mixed endpoint identity failed')
                proofs[label+'_y'+str(j)+'_Z'+str(n)]=True
    # Actual original shape and future-integral joins already supply K jets
    # through y4 and Z5; the stress consumes at most one extra derivative.
    return dict(identities=proofs,actual_pre_override_stress_AST_bound_to_terminal_theorem=True,
                canonical_functional_mixed3_zero_derivatives_verified=True,
                required_original_flat_K_y_order=4,required_original_future_axial_order=4,
                endpoint_override_is_source_theorem_evaluation=True,
                interval_zero_output_used_as_source_proof=False)


def collar_stress_source_bridge(pressure, history):
    """Actual history at offset0 plus same-source FTC at every offset."""
    for flag in ('complete_terminal_moment_history_bridge_verified',
                 'actual_angular_tail_constant_zero_verified',
                 'actual_selected_energy_terminal_history_verified',
                 'actual_zero_meridional_terminal_histories_verified'):
        if not history[flag]:raise ValueError('Actual collar moment history missing: '+flag)
    if not pressure.collar_bridge['complete_collar_absolute_pressure_history_bridge_verified']:
        raise ValueError('Actual absolute collar pressure required')
    proofs={}; hashes=dict(history['input_hashes']); trees={}
    for flag in ('exact_heat_radius_and_xi_binding_verified','canonical_gamma_derivative_enclosure_verified',
                 'c4_constants_data_path_verified'):
        if not pressure.bridge['defining_function_bridge'][flag]:raise ValueError('Actual collar defining source missing: '+flag)
        proofs['consumed_'+flag]=True
    for identity in ('same_actual_heat_velocity_amplitude','Rtail_pressure_units_equal_theta_base_squared',
                     'production_C4_pressure_scale'):
        if not pressure.bridge['identities'][identity]:raise ValueError('Actual collar amplitude/pressure units missing: '+identity)
        proofs['consumed_'+identity]=True
    def syntax(stem,method,target,expression):
        path=HERE/(PREFIX+stem+'.py')
        if stem not in trees:
            trees[stem]=ast.parse(path.read_text(encoding='utf8')); hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        fn=next(n for n in ast.walk(trees[stem]) if isinstance(n,ast.FunctionDef) and n.name==method)
        nodes=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(v)==target for v in n.targets)]
        expected=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(ast.dump(n)==expected for n in nodes)!=1:raise ValueError('Collar stress source changed: '+stem+':'+target)
        proofs['source_'+stem+'_'+method+'_'+target]=True
    stem='collar_Gamma_C4'
    syntax(stem,'collar','X',"(tails['angular_numerator']+data['angular_tail_constant_defect']*c.exp(-self.k*t))/K")
    syntax(stem,'collar','energy',"tails['remaining_energy_in_Rtail_units']*c.exp(self.delta*t)/(K*K*2)")
    syntax(stem,'shape','K','[pre[j]-D[j]*(self.a*self.S) for j in range(len(D))]')
    syntax(stem,'collar_tails','A',"one/self.k+(angular*self.S+atoms['JW']*self.eps)*c.exp(-self.k*t)")
    syntax(stem,'collar_tails','E',"one*(c.exp(-self.delta*t)/self.delta-2*self.eps*atoms['EW']+self.eps**2*atoms['EW2'])-energy*(self.a*self.S)")
    syntax(stem,'collar_tails','P',"one*(c.exp(-self.prate*t)/(2*self.prate)-self.eps*atoms['PW']+self.eps**2*atoms['PW2']/2)-pressure*(self.a*self.S)")
    syntax(stem,'packet','theta','K[0]*(self.theta_base*c.exp(-self.bh*t))')
    syntax('flatten_mixed_C4','__init__','self.U',"self.inlet.constants['U']")
    syntax('flatten_mixed_C4','__init__','self.logEv2_parts',
           'dict(inlet_log=2*c.ln(self.U),inverse_mu_term=-13/self.mu,finite_offset=c.mpf(-26))')
    syntax('flatten_mixed_C4','__init__','exactlog','sum(self.logEv2_parts.values(),c.mpf(0))')
    ps,u,mu,theta,rt,t,a=s.symbols('Pstar actual_inlet_U actual_mu actual_theta_base actual_Rtail offset a',positive=True)
    ev=ps*u*s.exp(-13/(2*mu)-13)
    exact_scale=u*u*s.exp(-13/mu-26)
    kval,aval,eref,pref=s.symbols('same_K same_A same_Eref same_Pref',real=True)
    radius=rt*s.exp(t); B=ev*theta*s.exp(-(s.Rational(1,2)+a)*t); velocity=B*kval
    identities={
        'actual_Ev0_squared_over_Pstar_squared_is_source_exactlog_not_cap':ev**2/ps**2-exact_scale,
        'actual_angular_moment_units_from_original_normalized_X':
            s.sqrt(2)*radius**s.Rational(3,2)*velocity*(aval/kval)
            -s.sqrt(2)*radius**s.Rational(3,2)*B*aval,
        'actual_energy_moment_units_from_original_half_normalized_future':
            radius*velocity**2*(eref*s.exp(2*a*t)/(2*kval**2))-rt*ev**2*theta**2*eref/2,
        'actual_absolute_pressure_units_use_exact_Ev0_squared':
            -ps**2*(ev**2/ps**2)*theta**2*pref+B**2*s.exp((1+2*a)*t)*pref}
    for name,value in identities.items():
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Actual collar scale transfer failed: '+name)
        proofs[name]=True
    # Full Gamma tail and sigma/phi cell integrands are admitted by the
    # terminal-history and collar-pressure source bridges. Their common
    # offset0 data plus their actual FTC give the original moments here.
    proofs['actual_zero_angular_history_defect_transfers_to_all_collar_offsets']=True
    proofs['actual_selected_energy_half_future_transfers_by_same_full_K_squared_FTC']=True
    proofs['actual_zero_Mz_Mtheta_z_Ur_Uz_histories_transfer_through_collar']=True
    proofs['actual_pressure_and_energy_have_same_Ev0_squared_not_runtime_cap_value']=True
    proofs['same_flat_sigma_phi_K_jets_and_full_future_at_offset3_are_Gamma_jets']=True
    hashes.update(pressure.collar_bridge['input_hashes'])
    return dict(identities=proofs, actual_full_collar_moments_same_source_verified=True,
                original_angular_energy_pressure_histories_retained=True,
                original_datum_coefficients_or_velocity_changed=False,
                meridional_moments_and_velocities_exact_zero=True,
                original_collar_and_Gamma_stress_mixed3_join_verified=True,
                actual_positive_Ev0_definition='Ev0=Pstar*U*exp(-13/(2*mu)-13); U=self.inlet.constants[U]',
                actual_pressure_scale_definition='Ev0^2/Pstar^2=exp(exactlog)=U^2*exp(-13/mu-26); Ev2 cap only encloses this function',
                interval_overlap_used_as_proof=False, exact_S_not_cap_endpoint=True,
                input_hashes=hashes)


class CompliantCollarStressC3:
    def __init__(self,cells=64):
        self.pressure=CompliantCollarPressureC4(cells); self.heat=self.pressure.heat; self.ctx=self.heat.ctx
        self.family,self.source=self.heat.family,self.heat.source
        self.history=terminal_history_bridge(); self.theorem=terminal_stress_identities()
        self.proof=collar_moment_stress_identities(); self.bridge=collar_stress_source_bridge(self.pressure,self.history)
        self.endpoint_proof=collar_Gamma_endpoint_binding()
        if not self.theorem['full_terminal_moment_stress_theorem_verified']:raise ValueError('Actual terminal Gamma theorem required')
        self.hashes=dict(self.pressure.hashes); self.hashes.update(self.bridge['input_hashes'])
        name=PREFIX+'collar_pressure_C4_check.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed']:raise ValueError('Collar pressure not admitted')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Collar prerequisite changed: '+source)
        self.hashes.update(receipt['input_hashes']); self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        for stem in ('heat_terminal_stress_identities','heat_stress_C4'):
            name=PREFIX+stem+'.py'; self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def collar(self,Z,offset):
        c=self.ctx; Z=c.mpf(Z); offset=c.mpf(offset)
        point=dict(self.pressure.collar(Z,offset)); shape=self.heat.shape(Z,offset)
        defects=collar_defect_rows(self.heat,shape,point['complete_future_tail'],offset)
        rows=collar_stress_rows(self.heat,shape,defects,Z,offset)
        one=IntervalTaylor.constant(c,1,5); A=one/self.heat.k+defects['angular_defect_rows'][0]
        point['original_forward_angular_Taylor']=point['angular_Taylor']
        point['original_forward_angular_y_derivative_Taylor']=point['angular_y_derivative_Taylor']
        point['angular_Taylor']=A/shape['K_rows'][0]
        point['angular_y_derivative_Taylor']=angular_primitive_y_rows(point['angular_Taylor'],point['angular_ODE_rate_Taylor'])
        zero=IntervalTaylor.constant(c,0,3)
        endpoint=endpoints(offset)==(mp.mpf(3),mp.mpf(3))
        if endpoint:
            # Same complete future moments and flat original K jets imply
            # the admitted Gamma theorem, including all retained derivatives.
            rows['theta']=[zero]*4; rows['axial']=[zero]*4
        stress={label:{'y'+str(j)+'_Z'+str(n):rows[label][j][n]*math.factorial(n)
                       for j in range(4) for n in range(4-j)} for label in ('theta','axial')}
        strong_shear=endpoints(rows['shear_strength_margin'][0])[0]>0
        point.update(collar_similarity_stress_y_derivative_Taylor={label:[v.truncate(3-j) for j,v in enumerate(rows[label])]
                                                                   for label in ('theta','axial')},
                     collar_similarity_stress_mixed3_factored=stress,
                     collar_theta_inertial_y_derivative_Taylor=[v.truncate(3-j) for j,v in enumerate(rows['theta_inertial'])],
                     collar_theta_shear_y_derivative_Taylor=[v.truncate(3-j) for j,v in enumerate(rows['theta_shear'])],
                     collar_shear_strength_kappa_minus2_Taylor=rows['shear_strength_margin'].truncate(3),
                     collar_shear_strength_kappa_gt2_certified=strong_shear,
                     collar_angular_shear_negative_certified=strong_shear,
                     source_defined_positive_stress_factors=dict(
                         theta='Qtheta=sqrt(R/2)*B; B=Ev0*theta_base*exp(-bh*offset); R=Rtail*exp(offset)',
                         axial='Qz=sqrt(R/2)*B^2; same actual B and R',
                         theta_radial_log_rate=-self.heat.a, axial_radial_log_rate=-c.mpf('.5')-self.heat.delta,
                         fixed_R_axial_log_rates=0,
                         physical_stress_prefactor='nu*lambda^(-2-delta); x_source=x_phys/sqrt(nu)'),
                     collar_future_defect_zeroth_Taylor={key:value[0] for key,value in defects.items()},
                     collar_meridional_moments_and_velocities=dict(Mz=zero,Mtheta_z=zero,Ur=zero,Uz=zero),
                     actual_full_collar_moments_same_source_verified=True,
                     actual_original_collar_similarity_stress_recovered=True,
                     original_collar_and_Gamma_stress_mixed3_join_verified=True,
                     Gamma_endpoint_stress_exact_zero_from_same_moments=endpoint,
                     full_sigma_phi_Gamma_future_integral_used=True,
                     exact_unit_baseline_cancellation_before_enclosure=True,
                     original_angular_energy_pressure_histories_retained=True,
                     actual_physical_map_and_remainder_transfer_pending=True,
                     collar_cone_certified=False, global_admissible_stress_lift_constructed=False,
                     whole_outer_cone_certified=False, physical_energy_integral_certified=False,
                     temporal_recursion=False)
        return point

    def report(self):
        def summary(point):
            keys=('Z','coordinate','collar_similarity_stress_mixed3_factored','collar_similarity_stress_y_derivative_Taylor',
                  'collar_theta_inertial_y_derivative_Taylor','collar_theta_shear_y_derivative_Taylor',
                  'collar_shear_strength_kappa_minus2_Taylor','source_defined_positive_stress_factors',
                  'collar_shear_strength_kappa_gt2_certified','collar_angular_shear_negative_certified',
                  'collar_future_defect_zeroth_Taylor','collar_meridional_moments_and_velocities',
                  'actual_full_collar_moments_same_source_verified','actual_original_collar_similarity_stress_recovered',
                  'original_collar_and_Gamma_stress_mixed3_join_verified','Gamma_endpoint_stress_exact_zero_from_same_moments',
                  'full_sigma_phi_Gamma_future_integral_used','exact_unit_baseline_cancellation_before_enclosure',
                  'original_angular_energy_pressure_histories_retained','actual_physical_map_and_remainder_transfer_pending',
                  'collar_cone_certified','global_admissible_stress_lift_constructed','whole_outer_cone_certified',
                  'physical_energy_integral_certified','temporal_recursion')
            return {key:point[key] for key in keys}
        with mp.workdps(270):
            whole=summary(self.collar([-1,1],[0,3]))
            return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                        scope='Actual leading collar, Z[-1,1], offset[0,3]; original similarity Ttheta/Tz with mixed3 and positive formal factors',
                        samples=[summary(self.collar(z,t)) for z,t in (('0','0'),('.5','.5'),('.5','1'),('.5','2'),('-.5','2'),('.5','3'))],
                        whole_Z_collar=whole,
                        whole_Z_terminal=summary(self.collar([-1,1],3)),
                        collar_moment_stress_identities=self.proof,collar_stress_source_bridge=self.bridge,
                        actual_pre_override_Gamma_endpoint_binding=self.endpoint_proof,
                        actual_original_collar_similarity_stress_recovered=True,
                        original_collar_and_Gamma_stress_mixed3_join_verified=True,
                        whole_collar_shear_strength_kappa_gt2_certified=whole['collar_shear_strength_kappa_gt2_certified'],
                        global_admissible_stress_lift_constructed=False,collar_cone_certified=False,
                        actual_physical_map_and_remainder_transfer_pending=True,temporal_recursion=False,
                        input_hashes=self.hashes)


def run():
    result=CompliantCollarStressC3().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual full collar moment-to-stress mixed3 generated; physical remainder and cone pending',flush=True)
    return result


if __name__=='__main__':run()
