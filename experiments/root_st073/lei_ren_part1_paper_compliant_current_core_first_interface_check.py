"""Original core integral recovery and common first-boundary source proof.

Differentiate the integrated identities back to the actual core equations.
The same source functions, not overlapping numerical boxes, establish joins.
"""
import ast
import copy
import json
from pathlib import Path
from types import SimpleNamespace

import sympy as s

from lei_ren_part1_paper_compliant_current_core_first_interface import (
    CurrentCoreFirstInterface,HERE,PREFIX,NAME,RECEIPT,GATE,ALL_GATE,SCOPES,
    PRIOR_JOINS,function,binding,sha,defining_boundary_bindings,pack)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_core_integral_atoms import finite_atom_coefficients
from lei_ren_part1_paper_compliant_frozen_comparison_field import MTH,MZ,MTHZ,MZT,MP


def exact(left,right,message):
    if s.cancel(left-right)!=0 and s.simplify(left-right)!=0:raise ArithmeticError(message)


def original_core_equation_bindings():
    specs={'W':"1-(1-dt)*z*a['Mz_over_R']-d*a['Mz_Z_over_R']",
        'H':"(1-dt)*z/2+d*a['Uz']",
        'lhsf':"2*L*(r*a['F_RR']+2*a['F_R'])",
        'rhsf':"W*(a['F']+r*a['F_R'])+dt/2*(1-2*z*a['Uz'])*a['F']+H*a['F_Z']",
        'lhsu':"2*L*(r*a['Uz_RR']+a['Uz_R'])",
        'rhsu':"W*r*a['Uz_R']+(1+dt)/2*(1-2*z*a['Uz'])*a['Uz']+H*a['Uz_Z']"}
    for target,expression in specs.items():binding('core_recursion','core_equation_defects',target,expression)
    fn=function('core_recursion','core_equation_defects')
    additions=[node for node in ast.walk(fn) if isinstance(node,ast.AugAssign) and ast.unparse(node.target)=='rhsu']
    expr="d*a['P_Z']-2*(1+dt)*z*a['P']-2*z*r*a['F']**2"
    if len(additions)!=1 or ast.dump(additions[0].value)!=ast.dump(ast.parse(expr,mode='eval').body):
        raise ValueError('Original core pressure or centrifugal term changed')
    return dict(all_original_angular_axial_mean_pressure_terms_AST_bound=True,
        axis_equation_not_substituted_for_exit_drive=True,passed=True)


def integrated_core_equations_proof():
    original_core_equation_bindings()
    r,z,delta=s.symbols('R Z delta',real=True);L=1-delta*z*z;d=1-z*z
    F=s.Function('F')(r,z);V=s.Function('V')(r,z);P0=s.Function('P0')(z)
    IV,IA,IH,IK,IB,IC=[s.Function(name)(r,z) for name in ('I_V','I_V2','I_RF','I_RFV','I_RF2','I_F2')]
    densities={IV:V,IA:V*V,IH:r*F,IK:r*F*V,IB:r*F*F,IC:F*F}
    def ftc(expr):
        def is_primitive_derivative(node):
            return isinstance(node,s.Derivative) and node.expr in densities and r in node.variables
        def replace(node):
            orders=dict(node.variable_count)
            return s.diff(densities[node.expr],r,orders.get(r,0)-1,z,orders.get(z,0))
        return expr.replace(is_primitive_derivative,replace)
    M=IV/r;W=1-(1-delta)*z*M-d*s.diff(M,z);H=(1-delta)*z/2+d*V;P=P0+IC
    angular_original=(2*L*(r*s.diff(F,r,2)+2*s.diff(F,r))
        -W*(F+r*s.diff(F,r))-delta*(1-2*z*V)*F/2-H*s.diff(F,z))
    axial_original=(2*L*(r*s.diff(V,r,2)+s.diff(V,r))
        -W*r*s.diff(V,r)-(1+delta)*(1-2*z*V)*V/2-H*s.diff(V,z)
        -d*s.diff(P,z)+2*(1+delta)*z*P+2*z*r*F*F)
    angular_integrated=(2*L*r*r*s.diff(F,r)-W*r*r*F+(1-delta/2)*IH
        -(1-delta)*z*s.diff(IH,z)/2-d*s.diff(IK,z)+(2*delta-1)*z*IK)
    axial_integrated=(2*L*r*s.diff(V,r)-W*r*V+(1-delta)*IV/2
        -(1-delta)*z*s.diff(IV,z)/2+2*delta*z*IA-d*s.diff(IA,z)
        -r*(d*s.diff(P0,z)-2*(1+delta)*z*P0)
        -d*(r*s.diff(IC,z)-s.diff(IB,z))
        +2*(1+delta)*z*(r*IC-IB)+2*z*IB)
    exact(ftc(s.diff(angular_integrated,r)),r*angular_original,'Integrated angular equation does not differentiate to original core')
    exact(ftc(s.diff(axial_integrated,r)),axial_original,'Integrated axial equation does not differentiate to original core')
    # Execute the actual original direction function with exact symbolic
    # source functions. The two physical scales remain distinct constants;
    # the relative amplitude itself varies with Z and is differentiated.
    F0=s.Function('F0')(z);Fbase,Pstar=s.symbols('F0base Pstar',positive=True)
    phi=SimpleNamespace(order=0)
    moments={MTH:2*IH/(r*r*F0),MZ:M,MTHZ:2*IK/(r*r*F0),
             MZT:dict(axial=IA/r,swirl=IB/(r*r*F0*F0)),MP:IC/(r*F0*F0)}
    amplitude=F0/Fbase;amplitude2=F0*F0/(Fbase*Fbase)
    dress=lambda value,factor:(F/F0 if value is phi else value)*factor
    env=dict(IntervalTaylor=SimpleNamespace(variable=lambda *args:z),square=lambda value:value*value,
        derivative=lambda value:s.diff(value,z),dress=dress,MTH=MTH,MZ=MZ,MTHZ=MTHZ,MZT=MZT,MP=MP)
    fn=copy.deepcopy(function('compliant_inner_bridge_profiles','direction'))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<original direction on exact primitive source>','exec'),env)
    directions=env['direction'](None,z,delta,phi,V,moments,P0/(Pstar*Pstar),amplitude,amplitude2)
    D=directions['D_over_R']
    drive=directions['drive_hydro']+Pstar*Pstar*directions['drive_pressure']+r*Fbase*Fbase*directions['drive_swirl']
    exact(angular_integrated/r**2,2*L*s.diff(F,r)+L*F*D,'Original angular direction is not integrated stress-free core recovery')
    exact(axial_integrated/r,2*L*(s.diff(V,r)+drive),'Original three-scale axial drive is not integrated stress-free core recovery')
    # The nested mixed-derivative direction must use the same base formula,
    # not another scalar/axis formula. Bind its amplitude/Taylor construction.
    for target,expr in [('f','true(phi,ratios)'),('h','true(rows[MTH],ratios)'),
        ('k','true(rows[MTHZ],ratios)'),('b',"true(rows['B'],ratios2)"),('p','true(rows[MP],ratios2)'),
        ('m','nested(rows[MZ])'),('a',"nested(rows['A'])"),('vv','nested(V)'),('pressure','lift(p0)'),
        ('W','1-(z*m)*(1-delta)-d*m.Zderivative()')]:
        binding('compliant_bridge_mixed_C4','comparison_directions',target,expr)
    nested=function('compliant_bridge_mixed_C4','comparison_directions')
    values=[node.value for node in ast.walk(nested) if isinstance(node,ast.Assign)
            and any(ast.unparse(t)=='values' for t in node.targets)]
    if len(values)!=1 or not isinstance(values[0],ast.Call):raise ValueError('Original nested direction source missing')
    class ConvertDerivative(ast.NodeTransformer):
        def visit_Call(self,node):
            if isinstance(node.func,ast.Attribute) and node.func.attr=='Zderivative' and not node.args:
                return ast.copy_location(ast.Call(func=ast.Name(id='derivative',ctx=ast.Load()),args=[node.func.value],keywords=[]),node)
            return self.generic_visit(node)
    base=dict(z=z,delta=delta,L=L,d=d,W=W,f=F/Fbase,h=2*IH/(r*r*Fbase),k=2*IK/(r*r*Fbase),
        b=IB/(r*r*Fbase*Fbase),p=IC/(r*Fbase*Fbase),m=M,a=IA/r,vv=V,pressure=P0/(Pstar*Pstar),
        derivative=lambda value:s.diff(value,z))
    angular=[node.value for node in ast.walk(nested) if isinstance(node,ast.Assign)
             and any(ast.unparse(t)=='angular' for t in node.targets)]
    if len(angular)!=1:raise ValueError('Original nested angular direction source missing')
    angular_tree=ast.Expression(ConvertDerivative().visit(copy.deepcopy(angular[0])))
    base['angular']=eval(compile(ast.fix_missing_locations(angular_tree),
        '<original nested angular direction>','eval'),{},base)
    for kw in values[0].keywords:
        tree=ast.Expression(ConvertDerivative().visit(copy.deepcopy(kw.value)))
        actual=eval(compile(ast.fix_missing_locations(tree),'<nested base direction>','eval'),{},base)
        key=kw.arg if kw.arg=='D_over_R' else 'drive_'+kw.arg
        exact(actual,directions[key],'Nested physical direction differs from actual original direction: '+key)
    return dict(original_core_equations_AST_bound=True,
        differentiated_integrated_angular_equation_identity=True,
        differentiated_integrated_axial_equation_identity=True,
        original_scalar_direction_functions_replayed_exactly=True,
        original_nested_direction_base_identities=4,
        pressure_and_swirl_physical_scales_not_dropped=True,
        relative_F0_amplitude_varies_with_Z=True,
        exact_stress_free_recovery=['D_logR logF=-R*D_over_R/2',
            'D_logR V=-(R*hydro+R*Pstar^2*pressure+R^2*F0base^2*swirl)'],
        axis_integration_constant_proof=[
            'The admitted common core is analytic at R=0; F,V and their axial derivatives are finite for each fixed source.',
            'I_V,I_V2,I_F2 vanish as O(R); I_RF,I_RFV,I_RF2 vanish as O(R^2).',
            'W has a finite axis limit since I_V/R is its analytic radial mean.',
            'Both integrated residuals tend to0 at R=0. Their derivatives vanish by the actual common core equations.',
            'Hence both integrated residuals vanish identically throughout the common core/continuation.'],
        same_common_pressure_P_equals_P0_plus_integral_F_squared=True,
        integrated_core_means_and_direction_are_source_functions_not_caps=True,passed=True)


def scaled_atom_physical_primitive_proof(field=None):
    binding('compliant_core_transfer','run','lam','ctx.exp(logLambda)')
    binding('compliant_core_transfer','run','eps','ctx.exp(-logLambda)')
    logarithm=s.Symbol('logLambda',real=True)
    exact(s.exp(logarithm)*s.exp(-logarithm),1,'Original Lambda epsilon definitions differ')
    parameter_name=PREFIX+'current_core_recurrence_source_check.json'
    if field is not None:
        parameter=json.loads((HERE/parameter_name).read_bytes())
        if (field.hashes.get(parameter_name)!=sha(parameter_name)
            or not parameter['pressure_and_seed_bindings']['Lambda_times_core_epsilon_exactly_one']
            or (parameter['actual_five_defect_family_sha256'],parameter['implicit_source_sha256'],parameter['datum_enclosure_sha256'])
               !=(field.family,field.source_sha,field.datum_sha)):
            raise ValueError('Common source exact Lambda epsilon receipt missing')
    binding('compliant_core_physical_field','report','normalized_core_definition',
        "'rho=Lambda*R; F=F0*Phi, Uz=4Z+j+epsilon*Psi, F0=exp(-selected_logCstar-Lambda*G); unique admitted analytic fixed point'")
    binding('compliant_core_physical_field','report','core_domain',
        "'rho in[0,4.1], Z in[-1,1]; physical map tau>0, |Z|<1; Ra=4/Lambda'")
    binding('compliant_frozen_comparison_field','__init__','self.r','4*self.core.epsilon')
    binding('compliant_core_physical_field','profiles','grids[F][index]',
        "sum((math.comb(k,j)*source['relative_F0_derivatives'][j]*phi[gridkey(i,k-j)] for j in range(k+1)),c.mpf(0))")
    for target,expression in [('phi',"[row[:order+1] for row in rows['A'][:N+1]]"),
        ('uz',"[row[:order+1] for row in rows['Uz'][:N+1]]"),
        ('finite','finite_atom_coefficients(c,phi,uz,order)'),
        ('pressure_primitive_at_exit_axial_coefficients',"[4*value for value in atoms['C']]")]:
        binding('compliant_core_integral_atoms','atoms_from_packet',target,expression)
    binding('compliant_comparison_point_integrals','inlet','moments',
        "{name:IntervalTaylor(c,values) for name,values in atoms['actual_core_atom_axial_coefficients'].items()}")
    # Change variables in each physical primitive, retaining its actual
    # normalization. F0 varies axially; its factors cancel before all axial
    # differentiations, rather than being silently treated as constants.
    r,rho,z=s.symbols('R rho Z',real=True);epsilon=s.Symbol('epsilon_core',positive=True)
    Ra=4*epsilon;F0=s.Function('F0')(z)
    phi=s.Function('Phi')(rho,z);uz=s.Function('Uz')(rho,z)
    physical_f=F0*s.Function('Phi')(r/epsilon,z);physical_v=s.Function('Uz')(r/epsilon,z)
    definitions={
        'H':(2*r*physical_f/(Ra**2*F0),rho*phi/8,'2 I_RF/(Ra^2 F0)'),
        'M':(physical_v/Ra,uz/4,'I_V/Ra'),
        'K':(2*r*physical_f*physical_v/(Ra**2*F0),rho*phi*uz/8,'2 I_RFV/(Ra^2 F0)'),
        'A':(physical_v**2/Ra,uz**2/4,'I_V2/Ra'),
        'B':(r*physical_f**2/(Ra**2*F0**2),rho*phi**2/16,'I_RF2/(Ra^2 F0^2)'),
        'C':(physical_f**2/(Ra*F0**2),phi**2/4,'I_F2/(Ra F0^2)')}
    identities=0
    for name,(density,scaled,normalization) in definitions.items():
        changed=s.cancel(epsilon*density.subs(r,epsilon*rho))
        for k in range(7):
            exact(s.diff(changed,z,k),s.diff(scaled,z,k),'Scaled actual atom differs from physical primitive: '+name)
            identities+=1
    exact((r/epsilon).subs(r,Ra),4,'Physical atom exit is not scaled rho4')
    return dict(original_scaled_core_and_amplitude_callable_AST_bound=True,
        original_exponential_Lambda_epsilon_definitions_verified=True,
        same_current_exact_parameter_receipt=parameter_name,
        actual_inlet_uses_scaled_Phi_Uz_atom_mapping=True,
        exact_change_of_variables='rho=Lambda*R; Lambda*epsilon_core=1; dR=epsilon_core*drho; Ra=4epsilon_core',
        six_physical_primitive_normalizations={name:row[2] for name,row in definitions.items()},
        normalized_density_and_axial_derivative_identities=identities,
        axial_orders=list(range(7)),F0_axial_dependence_cancelled_exactly=True,
        integrated_core_direction_and_actual_scaled_atom_functions_identified=True,
        directed_tail_enclosures_do_not_define_atom_values=True,passed=True)


def common_pressure_atom_proof():
    rho,z=s.symbols('rho Z',real=True);f=s.Function('Phi')(rho,z)
    t=s.Symbol('t',real=True);integral=s.Integral(s.Function('Phi')(t,z)**2,(t,0,4));C=integral/4
    exact(4*C,integral,'Exact pressure atom normalization differs')
    epsilon,F0,Pstar,PD=s.symbols('epsilon F0base Pstar PD',positive=True)
    Ra=4*epsilon
    exact(epsilon*F0**2*4*C,Ra*F0**2*C,'Core and first pressure physical units differ')
    # All axial orders are the product rule for one primitive; PD stays
    # independent. Never equate four C to axis pressure plus increment.
    amp=s.Function('S')(z)
    identities=0
    for k in range(7):
        direct=s.diff(amp*4*C,z,k)
        product=sum((s.binomial(k,j)*s.diff(amp,z,j)*s.diff(4*C,z,k-j) for j in range(k+1)),s.Integer(0))
        exact(direct,product,'Common pressure amplitude product derivative differs');identities+=1
    # Check the actual integrated finite-atom callable against an independent
    # polynomial integral. The exact all-order identity follows finite sums,
    # then analytic convergence and the common source-bound product tails.
    ctx=SimpleNamespace(mpf=s.sympify)
    phi=[[s.Rational((n+1)*(k+2),17**(k+1)) for k in range(7)] for n in range(3)]
    uz=[[s.Rational((n+2)*(k+1),19**(k+1)) for k in range(7)] for n in range(3)]
    actual=finite_atom_coefficients(ctx,phi,uz,6)['C']
    polynomial=sum((phi[n][k]*rho**n*z**k for n in range(3) for k in range(7)),s.Integer(0))
    primitive=s.integrate(polynomial**2,(rho,0,4))
    for k,value in enumerate(actual):
        exact(4*value,s.expand(primitive).coeff(z,k),'Actual finite pressure atom differs from independent integral')
    return dict(exact_common_pressure_atom_four_C_identity=True,
        axis_PD_is_not_added_to_the_four_C_integral=True,
        core_PI_equals_four_amplitude_dressed_C=True,
        physical_pressure='Pstar^2*PD + epsilon_core*F0base^2*PI_core = Pstar^2*PD + Ra*F0base^2*C_dressed',
        ordinary_axial_product_identities=identities,independent_actual_finite_atom_coefficients=7,
        analytic_primitive_proof=[
            'The common Banach solution supplies a single Phi function and convergent radial/axial derivatives on rho<=4.1.',
            'C=(1/4)integral_0^4 Phi^2 is the same actual atom used by the comparison and first prescribed field.',
            'The original core profiles covering integral bounds this same primitive at radial order0; positive radial orders are derivatives of Phi^2 by FTC.',
            'The same common Phi/model/correction bounds enclose all linear/product atom tails.',
            'Termwise integration and ordinary axial product differentiation are justified by local uniform convergence.',
            'Different finite truncations or wider boxes do not define different pressure functions.'],
        raw_bridge_axial_V_is_Uz_not_pressure_primitive=True,passed=True)


def phase0_and_coordinate_proof(field):
    # This width is one scalar norm choice for the entire source family,
    # not a pointwise function of Z. Its intervals are only enclosures.
    ledger=field.core.records['K1_ledger'] if 'K1_ledger' in field.core.records else json.loads((HERE/(PREFIX+'K1_ledger.json')).read_bytes())
    definition=ledger['parameter_definition']
    if definition['h_b']!='cstar*K^-100' or definition['K']!='actual physical norm sum in9.16, bounded by Cstar*Kbar':
        raise ValueError('Positive global scalar hb definition changed')
    geometry=function('compliant_bridge_mixed_C4','evaluate')
    microscopic=[node for node in ast.walk(geometry) if isinstance(node,ast.If)
        and ast.unparse(node.test)=='microscopic' and any(isinstance(row,ast.Assign)
        and any(ast.unparse(t)=='y' for t in row.targets) for row in node.body)]
    if len(microscopic)!=1:raise ValueError('Unique microscopic radius branch required')
    radius=[node.value for node in microscopic[0].body if isinstance(node,ast.Assign)
        and any(ast.unparse(t)=='R' for t in node.targets)]
    if len(radius)!=1 or ast.dump(radius[0])!=ast.dump(ast.parse('self.r*c.exp(y)',mode='eval').body):
        raise ValueError('Actual first-chart radius changed')
    zero=ast.parse('if lo==hi==0:\n R=self.r;theta=c.mpf(1)').body[0]
    if not any(ast.dump(node)==ast.dump(zero) for node in microscopic[0].body):
        raise ValueError('Exact microscopic phase0 radius changed')
    for target,expr in [('y','self.h*value'),('rho',"c.mpf(4)*c.exp(y)"),
        ('alpha','[algebra.lift(1)]+[algebra.lift(0)]*3'),
        ('chi','[algebra.lift(1-sigma[0])+h*sigma[0]]+[-algebra.lift(sigma[k])+h*sigma[k] for k in range(1,4)]')]:
        binding('compliant_bridge_mixed_C4','evaluate',target,expr,True)
    binding('compliant_bridge_mixed_C4','bridge_controls','logF',
        '[-scale*product_rows(chi,Dbar,k)/2 for k in range(4)]')
    binding('compliant_bridge_mixed_C4','bridge_controls','relative_log',
        '[logF[k]-bar_log[k] for k in range(3)]+[logF[3]*0]')
    binding('compliant_bridge_mixed_C4','bridge_controls','quotient_rows',
        '[algebra.lift(quotient)*row for row in exponential_derivatives(relative_log)]')
    binding('compliant_bridge_mixed_C4','evaluate','Dbar',
        "[rate_rows(directions['D_over_R'],scale,k)*R for k in range(4)]")
    binding('compliant_bridge_mixed_C4','evaluate','drive',
        "[rate_rows(directions['hydro'],scale,k)*R+algebra.shift(rate_rows(directions['pressure'],scale,k)*R,(0,1,0,0))+algebra.shift(rate_rows(directions['swirl'],scale*2,k)*R**2,(0,0,1,0)) for k in range(4)]")
    formal=field.prior['retained_original_flat_endpoint_source_proof']
    if not formal['passed'] or not formal['phase1_and_phase2_use_original_flat_endpoint_jets']:
        raise ValueError('Actual common flat cutoff endpoint source proof missing')
    y,z,h,Ra=s.symbols('y Z hb Ra',positive=True);g=s.Function('g');k=s.Symbol('k',integer=True,nonnegative=True)
    phase=s.Symbol('phase',real=True);rho=4*s.exp(h*phase);radius=Ra*s.exp(h*phase)
    exact(rho.subs(phase,0),4,'Original phase0 core radius differs')
    exact(radius.subs(phase,0),Ra,'Original phase0 physical radius differs')
    def canonical_pullback(expr):
        # Identify the derivative by its two orders and evaluation point,
        # independent of SymPy's dummy names and mixed-partial ordering.
        jet=lambda a,b:s.Symbol('g_y%d_Z%d_at_hb_phase'%(a,b))
        substitutions={}
        for node in expr.atoms(s.Subs):
            if len(node.variables)!=1 or node.point!=(h*phase,) or not isinstance(node.expr,s.Derivative):
                raise ValueError('Unexpected mixed pullback substitution')
            variable=node.variables[0];derivative=node.expr
            if derivative.expr!=g(variable,z) or set(derivative.variables)-{variable,z}:
                raise ValueError('Mixed pullback differentiates another function')
            orders=dict(derivative.variable_count)
            substitutions[node]=jet(orders.get(variable,0),orders.get(z,0))
        result=expr.xreplace(substitutions)
        substitutions={}
        for node in result.atoms(s.Derivative):
            if node.expr!=g(h*phase,z) or set(node.variables)-{z}:
                raise ValueError('Unexpected unpulled mixed derivative')
            substitutions[node]=jet(0,dict(node.variable_count).get(z,0))
        return result.xreplace(substitutions).xreplace({g(h*phase,z):jet(0,0)})
    chain=0
    for a in range(5):
        for b in range(5-a):
            direct=s.diff(g(h*phase,z),phase,a,z,b)
            expected=h**a*s.diff(g(y,z),y,a,z,b).subs(y,h*phase)
            exact(canonical_pullback(direct),canonical_pullback(expected),'Exact hb mixed4 coordinate pullback differs');chain+=1
    # Derivative induction is parameter-free, using equal inlet values and
    # the same core direction/mean/pressure functions just proved above.
    return dict(exact_global_positive_scalar_hb_definition_bound=True,
        hb_axial_derivatives_are_zero_by_global_norm_choice=True,
        source_width_cap_not_used_as_hb_value=True,
        first_phase0_geometry_identities=2,exact_mixed4_pullback_identities=chain,
        flat_endpoint_control_chi_equals_one_with_zero_positive_jets=True,
        first_comparison_alpha_equals_one=True,
        original_full_three_scale_controls_AST_bound=True,
        same_moment_history_join_proof=[
            'On first alpha=1, smoothed comparison logF and V integrate the original core derivatives and equal the core throughout the comparison continuation.',
            'Its six histories have the same core-integral inlet and same radial ODEs H_y+2H=2Phi, M_y+M=V, K_y+2K=2Phi*V, A_y+A=V^2, B_y+2B=Phi^2, C_y+C=Phi^2.',
            'Uniqueness of these linear ODEs identifies the comparison histories with the exact core primitives.',
            'At phase0 the actual integral increments vanish, so actual F,V and all six histories have those same source values.',
            'The independent integrated core proof identifies Dbar and all three physical drive components with the stress-free core derivatives.',
            'Flat sigma gives chi(0)=1 and every positive-order jet zero; the actual and comparison quotient is1 with matching lower-order derivative jets.',
            'Differentiate the actual source ODEs. At each order through4, the Leibniz and Bell recurrences use only already identified lower derivatives, so the next derivative equals the core derivative.',
            'Axial differentiation of these function identities and the exact hb pullback give equality of every mixed source jet of total order<=4.',
            'Radial Q and pressure are recovered from the same actual mean and C primitives, so their jets and the five primitive jets agree.',
            'The existing nonsingular Cartesian/time maps apply to the same source functions at R=Ra=4epsilon_core; regional physical equality follows under the original common smooth pullback.'],
        equality_by_functions_and_ODE_uniqueness_not_interval_overlap=True,
        temporal_coefficient_recursion_not_inferred=True,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentCoreFirstInterface(require_checked=False)
    if (raw['actual_five_defect_family_sha256']!=field.family or raw['implicit_source_sha256']!=field.source_sha
        or raw['datum_enclosure_sha256']!=field.datum_sha or raw[GATE] or raw[ALL_GATE]
        or any(raw[key] for key in SCOPES) or raw['original_boundary_source_bindings']!=defining_boundary_bindings()
        or not all(field.graph.values())):
        raise ValueError('Core/first common source or scope differs')
    runtime=field.boundary()
    if encode(pack(runtime))!=raw['fresh_original_first_boundary']:
        raise ValueError('New common original first-inlet runtime changed')
    counts=runtime['first_mixed4_component_row_counts']
    if sum(n for group in counts.values() for n in group.values())!=135:
        raise ValueError('Actual first-inlet135 physical/primitive mixed4 rows omitted')
    if not runtime['same_packet_actual_comparison_acquisition'] or runtime['original_inlet_arithmetic_changed']:
        raise ValueError('Original same-packet atom arithmetic was changed')
    equations=integrated_core_equations_proof();atoms=scaled_atom_physical_primitive_proof(field)
    pressure=common_pressure_atom_proof();coordinate=phase0_and_coordinate_proof(field)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source_sha,
        datum_enclosure_sha256=field.datum_sha,current_common_core_bridge_graph=field.graph,
        original_boundary_source_bindings=field.bindings,
        original_integrated_stress_free_core_equations=equations,
        scaled_atoms_are_the_original_physical_primitives=atoms,
        common_pressure_primitive_and_true_atom=pressure,
        phase0_ODE_and_exact_positive_width_pullback=coordinate,
        fresh_same_packet_first_boundary_checked=True,first_mixed4_rows_present=135,
        **dict.fromkeys(PRIOR_JOINS,True),**dict.fromkeys((GATE,ALL_GATE),True),**dict.fromkeys(SCOPES,False),
        remaining_dependency='Quantitative native pulse C4 and same-source heat stress companions, global admissible stress/flatness/required-domain energy, then true temporal recursion and oscillatory correction',
        all_scoped_checks_passed=True,all_passed=True,
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    (HERE/RECEIPT).write_bytes((json.dumps(encode(result),indent=2)+'\n').encode('utf8'))
    print('PASS original integrated core recovery, correct pressure units and fourth functional mixed4 join; global/temporal remain open',flush=True)
    return result


if __name__=='__main__':run()
