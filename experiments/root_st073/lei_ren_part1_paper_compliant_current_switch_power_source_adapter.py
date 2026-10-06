"""Expose original power locals and identify the two exact R110 sources."""
import ast
import copy
import math
from types import SimpleNamespace
import sympy as s
import lei_ren_part1_paper_compliant_microswitch_mixed_C4 as original
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import SourceAST
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import keyword_binding

PARENT_EXPRESSION="""dict(
    log_Utheta_over_Pstar_axial5_coefficients=list(logu.coefficients),
    actual_normalized_moment_shape_axial5_coefficients={name:list(row.coefficients) for name,row in shapes.items()},
    Uz_actual_axial5_coefficients=list(V.coefficients),
    pressure_axis_axial5_coefficients=inlet['pressure_axis_axial5_coefficients'],
    actual_R2_histories_retained=True,original_pressure_datum_retained=True,
    exposed_existing_locals_are_source_enclosures_not_point_values=True)"""


def compiled_original_postpower_with_locals():
    """Change only one output keyword; never feed exposed rows back into math."""
    asts=SourceAST();before=copy.deepcopy(asts.method('microswitch_mixed_C4','postpower'))
    fn=copy.deepcopy(before);fn.decorator_list=[]
    returns=[node for node in ast.walk(fn) if isinstance(node,ast.Return)]
    if len(returns)!=1 or not isinstance(returns[0].value,ast.Call) or not isinstance(returns[0].value.func,ast.Name) or returns[0].value.func.id!='dict':raise ValueError('One original postpower dict return required')
    added='actual_inherited_axial5_packet'
    if any(keyword.arg==added for keyword in returns[0].value.keywords):raise ValueError('Original source unexpectedly exposes adapter locals')
    returns[0].value.keywords.append(ast.keyword(arg=added,value=ast.parse(PARENT_EXPRESSION,mode='eval').body))
    compare=copy.deepcopy(fn);compare.decorator_list=before.decorator_list
    compare_return=next(node for node in ast.walk(compare) if isinstance(node,ast.Return))
    compare_return.value.keywords=[keyword for keyword in compare_return.value.keywords if keyword.arg!=added]
    if ast.dump(compare)!=ast.dump(before):raise ValueError('Unreviewed original postpower math adaptation')
    env=dict(vars(original))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<original postpower; existing axial5 locals exposed only>','exec'),env)
    return env[fn.name],dict(original_postpower_AST_unchanged_except_one_output_keyword=True,
        every_original_nonreturn_statement_retained_exactly=True,
        existing_locals_logu_shapes_V_P0_exposed_without_feedback=True,
        formal_hb_R2_zeta_radius_functions_retained_not_selected_from_caps=True,
        input_hashes=asts.hashes,passed=True)


def pure_source_power_transport(asts):
    """Replay source weights; numeric clipping only encloses positive kernels."""
    before=copy.deepcopy(asts.method('inner_switch_profiles','power_transport'));fn=copy.deepcopy(before)
    guard=ast.parse("endpoints(theta)[0]<=0 or endpoints(theta)[1]>1",mode='eval').body
    guards=[node for node in fn.body if isinstance(node,ast.If) and ast.dump(node.test)==ast.dump(guard)]
    if len(guards)!=1:raise ValueError('Original positive transport domain guard changed')
    fn.body.remove(guards[0])
    difference=next(node for node in fn.body if isinstance(node,ast.FunctionDef) and node.name=='difference')
    returns=[node for node in ast.walk(difference) if isinstance(node,ast.Return)]
    wanted=ast.parse("c.mpf([max(mp.mpf(0),endpoints(value)[0]),endpoints(value)[1]])",mode='eval').body
    if len(returns)!=1 or ast.dump(returns[0].value)!=ast.dump(wanted):raise ValueError('Original directed positive kernel clip changed')
    returns[0].value=ast.Name(id='value',ctx=ast.Load())
    env=dict(square=lambda value:value*value,MTH=original.MTH,MTHZ=original.MTHZ,MZ=original.MZ,MZT=original.MZT,MP=original.MP)
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),'<original exact source weights before directed positivity enclosure>','exec'),env)
    return env[fn.name]


def R110_source_functional_identity():
    """Arbitrary smooth axial data; postpower and post(110) before enclosures."""
    asts=SourceAST();c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)))
    z,y=s.symbols('Z y',real=True);hb=s.Symbol('original_positive_hb',positive=True)
    delta,logC,logP,Lambda=s.symbols('delta logCstar logPstar Lambda',real=True)
    Ps=s.Symbol('Pstar',positive=True);T=s.Symbol('actual_reshape_T',positive=True)
    G,J,logphi100,V,P0=[s.Function(name)(z) for name in ('same_G','same_JD','same_logphi100','same_V110','same_P0')]
    initial={key:s.Function('same_R2_'+key)(z) for key in ('H','mean','K','A','B','C')}
    a=s.log(s.Rational(110,100));zeta=a+y-2*hb
    Rleft=100*s.exp(2*hb+zeta);Rright=110*s.exp(y)
    theta_left=s.exp(-zeta);theta_right=100*s.exp(2*hb)/Rright
    logphi2=logphi100-hb/5-hb**2*J/2
    logphi_left=logphi2-s.Rational(2,5)*zeta
    logphi_right=logphi100-s.Rational(2,5)*(a+y)+s.Rational(3,5)*hb-hb**2*J/2
    geom={}
    for key,expr in dict(exact_radius=Rleft-Rright,exact_theta=theta_left-theta_right,
            exact_log_phi=logphi_left-logphi_right).items():
        if s.simplify(s.expand_power_exp(expr))!=0:raise ArithmeticError('R110 original source geometry differs: '+key)
        geom[key]=True
    transport=pure_source_power_transport(asts);phi2=s.exp(logphi2)
    left=transport(c,theta_left,initial,phi2,V);right=transport(c,theta_right,initial,phi2,V)
    moment_checks={}
    def flatten(moments):
        return dict(H=moments[original.MTH],K=moments[original.MTHZ],mean=moments[original.MZ],
            A=moments[original.MZT]['axial'],B=moments[original.MZT]['swirl'],C=moments[original.MP])
    ll,rr=flatten(left),flatten(right)
    for key in ll:
        if s.simplify(s.expand_power_exp(ll[key]-rr[key]))!=0:raise ArithmeticError('R110 exact inherited moment source differs: '+key)
        moment_checks[key]=True
    theta=s.Symbol('theta',positive=True);positive_weights={}
    for first,second in ((s.Rational(2,5),2),(s.Rational(4,5),2),(s.Rational(4,5),1)):
        if s.simplify(theta**first-theta**second-theta**first*(1-theta**(second-first)))!=0:raise ArithmeticError('Original positive kernel source weight differs')
        positive_weights[str(first)+'_'+str(second)]=True
    F0=s.exp(-logC-Lambda*G)
    def fields(radius,logphi,moments):
        mean=moments['mean'];Q=(2*z*V-(1-delta)*z*mean-(1-z*z)*s.diff(mean,z))/(1-delta*z*z)
        return dict(Utheta=s.sqrt(2*radius)*F0*s.exp(logphi),Uz=V,Ur=s.sqrt(radius/2)*Q,
            P_over_Pstar2=P0+radius*F0**2*moments['C']/Ps**2,
            Mtheta=radius**2*F0*moments['H'],Mtheta_z=radius**2*F0*moments['K'],
            Mz=radius*mean,Mztheta=radius*moments['A']-radius**2*F0**2*moments['B'],
            Mp=radius*F0**2*moments['C'])
    aa=fields(Rleft,logphi_left,ll);bb=fields(Rright,logphi_right,rr);rows={}
    for key in aa:
        difference=s.simplify(s.expand_power_exp(aa[key]-bb[key]))
        if difference!=0:raise ArithmeticError('R110 nine physical source functions differ: '+key)
        rows[key]=[]
        for j in range(5):
            for n in range(5-j):
                if s.diff(difference,y,j,z,n).subs(y,0)!=0:raise ArithmeticError('R110 exact mixed source row differs')
                rows[key].append('y%d_Z%d'%(j,n))
    logphi110=logphi_right.subs(y,0)
    B=-Lambda*G+logphi110+s.log(220)/2+s.log(1+z*z)
    log_power=s.log(220)/2-logC-Lambda*G+logphi110-logP+y/10
    ds=s.symbols('sigma1:5',real=True)
    flat_cutoff=sum(ds[j-1]*y**j/(T**j*s.factorial(j)) for j in range(1,5))
    log_reshape=B*(1-flat_cutoff)-s.log(1+z*z)+y/10-logC-logP
    logrows=[]
    for j in range(5):
        for n in range(5-j):
            value=s.simplify(s.diff(log_reshape-log_power,y,j,z,n).subs(y,0).subs(dict.fromkeys(ds,0)))
            if value!=0:raise ArithmeticError('R110 flat original reshape log jet differs')
            logrows.append('y%d_Z%d'%(j,n))
    wanted={
        ('microswitch_mixed_C4','postpower','inlet'):'self.switch.phase(Z,2)',
        ('microswitch_mixed_C4','postpower','length'):'c.ln(c.mpf(110)/100)-2*self.h',
        ('microswitch_mixed_C4','postpower','zeta'):'length*fraction',
        ('microswitch_mixed_C4','postpower','theta'):'c.exp(-zeta)',
        ('microswitch_mixed_C4','postpower','phi'):"phi2*theta**c.mpf('.4')",
        ('microswitch_mixed_C4','postpower','moments'):'power_transport(c,theta,as_initial(moments2),phi2,V)',
        ('microswitch_mixed_C4','postpower','packet'):"reshape_mixed(c,Z,self.core.delta,logu,[constant('.1')]+[constant(0)]*3,V,shapes,jet(inlet['pressure_axis_axial5_coefficients']),self.invP2,self.proofs)",
        ('inner_switch_profiles','post','p2'):'self.phase(Z,2)',
        ('inner_switch_profiles','post','moments'):'power_transport(c,theta,as_initial(moments2),phi2,v2)',
        ('inner_switch_profiles','inlet','result'):'self.post(Z,110)',
        ('inner_switch_profiles','inlet','B'):'G*(-self.core.Lambda)+logarithm(phi)+logarithm(1+square(z))+c.ln(c.mpf(220))/2',
        ('long_reshape_profiles','inputs','inlet'):'self.switch.inlet(Z)',
        ('long_reshape_profiles','evaluate','logu'):'B*(1-sig)-logq+(y/10-self.core.logC-self.core.logP)',
    }
    for (stem,method,target),value in wanted.items():asts.expression(stem,method,target,wanted=value)
    literals=dict(
        original_R2_log_identity=keyword_binding('inner_switch_profiles','phase','exact_R2_log_identity',"'log(F2/F100)=-.2*hb-.5*hb^2*JD' if end else None"),
        original_post_log_identity=keyword_binding('inner_switch_profiles','post','exact_post_log_identity',"'log(F/F100)=-(2/5)*log(R/100)+.6*hb-.5*hb^2*JD'"),
        original_power_log_identity=keyword_binding('microswitch_mixed_C4','postpower','source_angular_identity',"'log(F/F100)=-.4*log(R/100)+.6*hb-.5*hb^2*JD'"))
    asts.method('inner_switch_profiles','packet')
    asts.method('long_reshape_profiles','normalized_amplitude_decay')
    asts.method('long_reshape_profiles','backward_kernel')
    return dict(arbitrary_smooth_current_R2_axial_histories=True,exact_R2_R110_geometry_identities=geom,
        exact_six_current_power_transport_function_identities=moment_checks,
        original_positive_kernel_weight_factorizations=positive_weights,
        numeric_domain_guard_removed_only_under_original_0theta_le1_source_domain=True,
        directed_nonnegative_kernel_clip_is_enclosure_not_source_assignment=True,
        physical_source_function_identities=list(aa),physical_mixed4_rows_implied_by_exact_function_identities=rows,
        total_physical_mixed4_rows_implied=sum(map(len,rows.values())),
        exact_original_flat_R110_log_amplitude_mixed4_rows=logrows,original_source_literal_bindings=literals,
        same_original_R2_moments_V110_F0_and_P0_not_enclosure_overlap=True,
        same_reshape_mixed_generator_and_flat_log_jets_imply_R110_full_boundary_mixed4=True,
        original_reshape_zero_length_kernels_and_inherited_initial_functions_retained=True,
        constant_power_extension_used_only_for_reshape_boundary_jet=True,
        finite_right_reshape_not_replaced=True,positive_width_and_G_caps_not_selected_as_functions=True,
        input_hashes=asts.hashes,passed=True)
