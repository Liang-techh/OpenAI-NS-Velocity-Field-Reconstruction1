"""Current native angular repair -> actual whole-Z heat terminal identity.

The prescribed waiting root belongs to the UNCORRECTED reference field.
The live angular terminal includes the two-bump correction. Their exact
source equations cancel the current Gamma angular defect. Pressure and
energy/selected-pulse source installation remain separate obligations.
"""
import ast
import json
from pathlib import Path
from types import MethodType,SimpleNamespace

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_exact_repair_branch import (
    CurrentExactRepairBranch,HERE,PREFIX,sha,binding,function,pack)
from lei_ren_part1_paper_compliant_current_collar_stress_mixed_C4 import (
    CurrentCollarStressMixedC4,collar_stress_mixed4,shape_radial5)
from lei_ren_part1_paper_compliant_current_heat_pressure_stress import mixed,constant_stress_rows
from lei_ren_part1_paper_compliant_collar_pressure_C4 import collar_pressure_rows
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import CompliantCollarGammaC4
from lei_ren_part1_paper_compliant_future_energy_high_jets import copy_jet
from lei_ren_part1_paper_compliant_axial_pulse_field import decay_integral
from lei_ren_part1_paper_compliant_outer_angular_repair import intersect
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_defect_rows
from lei_ren_part1_paper_compliant_current_steep_waiting_source import actual_angular_zero_future_bindings
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_angular_terminal_closure.json'
RECEIPT=PREFIX+'current_angular_terminal_closure_check.json'
GATES=('current_native_angular_repair_function_identification_certified',
       'current_heat_angular_terminal_constant_eliminated',
       'current_heat_angular_stress_mixed4_after_terminal_closure_recovered')
SCOPES=('current_exact_repair_installed_in_all_physical_charts',
    'current_heat_terminal_constants_eliminated','current_heat_pressure_terminal_constant_eliminated',
    'heat_exterior_stress_identity_certified','global_completed_tensor_admissibility',
    'admissible_stress_lift_constructed','full_background_NS_validation','physical_energy_integral_certified',
    'independently_bounded_flat_remainder','full_point_physical_field_evaluation',
    'full_cartesian_vector_derivatives_certified','temporal_recursion')
VIEWS={'whole_collar':('heat_collar',[-1,1],[0,3]),
    'collar_exit':('heat_collar',[-1,1],3),
    'whole_exterior':('heat_exterior',[-1,1],[3,mp.inf]),
    'fresh_exterior':('heat_exterior','.467',4)}


@source_precision
def current_forward_terminal(exact,Z):
    """Original forward angular/absolute-pressure formulas on the NEW repair.

    The native flatten angular and pressure data are independent of the
    selected future-energy input. No old outer.angular/steep.data is read.
    """
    native=exact.companion.heat;c=native.ctx;z=c.mpf(Z)
    box=lambda v:c.mpf(endpoints(v))
    repair=exact.repair;steep=native.steep;flatten=exact.flatten
    if exact.future.params is not repair.params or steep.future.params is not exact.future.params:
        raise ValueError('Parameter-only transition kernels must share the rebound exact parameter source')
    source=exact.angular5.angular(z)
    coeff=[copy_jet(c,v) for v in source['physical_coefficient_Taylor']]
    weights={name:box(v) for name,v in repair.weights.items()}
    f=flatten.flatten(z,100);r=1-native.mu;prate=1+2*native.mu
    L=box(-30*repair.params.log_mu);k=1-native.delta/2
    W=box(repair.angular.waiting);logone=box(repair.angular.waiting_logone)
    Ts=box(repair.params.Ts)
    bp=c.mpf('.5')+native.mu;bh=c.mpf('.5')+native.delta/2
    thetaR=c.exp(-100*bp-bp*L)/2
    thetaS=thetaR*c.exp(-bp-r/2)
    thetaQ=thetaS*c.exp(-c.mpf('1.5')*Ts)
    thetaT=thetaQ*c.exp(-c.mpf('1.5')+k/2)
    pastA=IntervalTaylor.constant(c,0,5);pastP=pastA
    for dj,center in zip(coeff,(-3,-1)):
        pastA+=dj*(c.exp(r*center)*weights['A'])
        pastP+=(dj*weights['B']+dj*dj*weights['D']/2)*c.exp(-prate*center)
    one=IntervalTaylor.constant(c,1,5)
    XR=one/r+(f['angular_Taylor']-1/r)*c.exp(-r*L)+pastA
    PR=f['pressure']['P_over_Pstar_squared']+decay_integral(c,prate,L)*(flatten.Ev2*c.exp(-100*prate)/8)
    PR+=pastP*(flatten.Ev2*c.exp(-prate*(100+L))/4)
    XS=(XR+steep.infull['angular'])*c.exp(-r/2)
    XQ=XS+Ts;XT=(XQ+steep.outfull['angular'])*c.exp(-k/2)
    PS=PR+steep.infull['pressure']*(flatten.Ev2*thetaR**2)
    PQ=PS+decay_integral(c,3,Ts)*(flatten.Ev2*thetaS**2/2)
    PT=PQ+steep.outfull['pressure']*(flatten.Ev2*thetaQ**2)
    Xtail=one/k+(XT-1/k)*c.exp(-k*W)
    Ptail=PT+decay_integral(c,1+native.delta,W)*(flatten.Ev2*thetaT**2/2)
    theta_base=thetaT*c.exp(-bh*W-logone)
    pressure_scale=flatten.Ev2*theta_base**2
    return dict(XR=XR,XS=XS,XQ=XQ,XT=XT,Xtail=Xtail,
        PR=PR,PS=PS,PQ=PQ,PT=PT,Ptail=Ptail,pastA=pastA,pastP=pastP,
        current_exact_angular_C5_coefficients=coeff,waiting=W,logone=logone,
        theta_base=theta_base,pressure_scale=pressure_scale,
        exact_current_inverse_radius_log_terms=repair.heat.logS_terms,
        current_exact_repair_directly_consumed=True,old_outer_terminal_not_consumed=True,
        transition_kernels_share_same_exact_parameter_source=True,
        selected_future_energy_not_installed=True)


def current_heat_shape_view(exact):
    """Restricted canonical heat methods with rebound exact radius source."""
    original=exact.companion.heat;c=original.ctx
    cap=c.mpf(endpoints(exact.repair.strong_S_cap))
    logS=sum((c.mpf(endpoints(v)) for v in exact.repair.heat.logS_terms.values()),c.mpf(0))
    native_S_enclosure=exact.pulse.factor(logS)
    source_S_enclosure=intersect(c,native_S_enclosure,c.mpf([0,endpoints(cap)[1]]))
    view=SimpleNamespace(ctx=c,mu=original.mu,delta=original.delta,a=original.a,
        eps=original.eps,k=original.k,bh=original.bh,prate=original.prate,
        Scap=cap,S=source_S_enclosure,cells=original.cells,
        shape_cache={},tail_cache={},gamma_cache={},
        exact_logRtail_terms=exact.repair.heat.logradius_terms,
        exact_logS_terms=exact.repair.heat.logS_terms,
        native_inverse_radius_factor_consumed=True,
        exact_inverse_radius_definition='S=exp(sum(exact_logS_terms)); Scap is only an enclosure')
    for method in ('shape','local_Gamma','collar_tails'):
        setattr(view,method,MethodType(getattr(CompliantCollarGammaC4,method),view))
    return view


def return_binding(stem,name,expression,allow_zero=False):
    fn=function(stem,name);returns=[n.value for n in ast.walk(fn) if isinstance(n,ast.Return)]
    expected=ast.dump(ast.parse(expression,mode='eval').body)
    other=[v for v in returns if ast.dump(v)!=expected]
    zero_return=ast.dump(ast.parse('c.mpf(0)',mode='eval').body)
    if sum(ast.dump(v)==expected for v in returns)!=1 or (other and (not allow_zero or any(ast.dump(v)!=zero_return for v in other))):
        raise ValueError('Angular terminal production return changed: '+stem+'.'+name)


def current_angular_source_proof(exact):
    """Translate live source formulas, including repair and raw waiting root."""
    identities={};bindings={}
    def zero(name,left,right=0):
        if s.simplify(s.expand_power_exp(left-right))!=0:raise ArithmeticError('Current angular source identity differs: '+name)
        identities[name]=True
    specs={
        ('compliant_current_power_angular_source','data'):{
            'source':'self.fifth.angular(Z)',
            'coeff':"[copy_jet(c,v) for v in source['physical_coefficient_Taylor']]",
            'f':'self.flatten.flatten(Z,100)'},
        ('compliant_power_angular_C4','angular'):{
            'data':'self.data(Z)','y':'self.Lrel+s',
            'Xbase':"eq+(data['flatten_exit_X']-1/self.rate)*c.exp(-self.rate*y)",
            'X':'(Xbase+pastA*c.exp(-self.rate*s))/F[0]'},
        ('compliant_current_steep_waiting_source','data'):{
            'terminal':'self.outer.angular(Z,0)','XR':"terminal['angular_Taylor']",
            'XS':"(XR+self.infull['angular'])*c.exp(-self.rate/2)",
            'XQ':'XS+self.Ts','XT':"(XQ+self.outfull['angular'])*c.exp(-self.k/2)"},
        ('compliant_steep_waiting_C4','waiting'):{
            'data':'self.data(Z)','t':'self.wait*phase',
            'X':"(data['XT']-1/self.k)*c.exp(-self.k*t)+1/self.k"},
        ('compliant_outer_angular_repair','preheat_difference'):{'start':'2*self.angular.Xv*c.exp(-100*self.rate)'},
        ('compliant_outer_angular_repair','__init__'):{
            'self.p':'c.exp(-2*self.rate)',
            'self.logscale':'30*self.rate*self.params.log_mu','self.scale':'c.exp(self.logscale)',
            'self.log_theta_multiplier':'(self.rate+self.angular.restore_rate)/2+self.angular.restore_rate*self.angular.waiting-self.angular.waiting_logone'},
        ('compliant_outer_angular_repair','coefficients'):{
            'b1':"r*(c.exp(self.rate)/self.weights['A'])"},
        ('compliant_outer_angular_repair','defects'):{
            'exact_r_definition':"'r=scale*(Xf(0)-Xf(Z))+exp(log_theta_multiplier)*S*Theta_hat(Z)'"},
        ('compliant_outer_angular_candidate','reference_buffer'):{
            'departure':"(inlet['X']-eq)*c.exp(-self.rate*t)",'X':'departure+eq'},
        ('compliant_outer_angular_candidate','steep_in'):{'X':"(inlet['X']+I)*c.exp(-self.rate*(t-J))"},
        ('compliant_outer_angular_candidate','before_waiting'):{'X':"(inlet['X']+I)*c.exp(-self.restore_rate*J)"},
        ('compliant_outer_angular_candidate','collar_preheat_integral'):{'bracket':'1-sig+sig*fv'},
        ('compliant_outer_angular_candidate','f'):{'y':'(3-t)/2'},
        ('compliant_collar_Gamma_C4','shape'):{
            'sig':'sigma_jets(c,t)','phi':'phi_jets(c,t)',
            'W':'[unity[j]-sr[j]+v for j,v in enumerate(product_rows(sr,fr))]',
            'C':'[unity[j]-fr[j]*self.eps for j in range(5)]'},
        ('compliant_collar_Gamma_C4','collar_tails'):{
            'A':'one/self.k+(angular*self.S+atoms[\'JW\']*self.eps)*c.exp(-self.k*t)'},
        ('compliant_collar_Gamma_C4','data'):{
            'terminal':'self.steep.waiting(Z,1)','heat0':'self.collar_tails(Z,0)',
            'defect':"Xtail*(1-self.eps)-heat0['angular_numerator']"},
        ('compliant_current_angular_terminal_closure','current_forward_terminal'):{
            'repair':'exact.repair','source':'exact.angular5.angular(z)',
            'coeff':"[copy_jet(c,v) for v in source['physical_coefficient_Taylor']]",
            'weights':'{name:box(v) for name,v in repair.weights.items()}',
            'f':'flatten.flatten(z,100)','L':'box(-30*repair.params.log_mu)',
            'W':'box(repair.angular.waiting)','logone':'box(repair.angular.waiting_logone)',
            'Ts':'box(repair.params.Ts)',
            'XR':"one/r+(f['angular_Taylor']-1/r)*c.exp(-r*L)+pastA",
            'XS':"(XR+steep.infull['angular'])*c.exp(-r/2)",
            'XQ':'XS+Ts','XT':"(XQ+steep.outfull['angular'])*c.exp(-k/2)",
            'Xtail':'one/k+(XT-1/k)*c.exp(-k*W)',
            'PR':"f['pressure']['P_over_Pstar_squared']+decay_integral(c,prate,L)*(flatten.Ev2*c.exp(-100*prate)/8)",
            'thetaR':'c.exp(-100*bp-bp*L)/2','thetaS':'thetaR*c.exp(-bp-r/2)',
            'thetaQ':"thetaS*c.exp(-c.mpf('1.5')*Ts)",
            'thetaT':"thetaQ*c.exp(-c.mpf('1.5')+k/2)",
            'PS':"PR+steep.infull['pressure']*(flatten.Ev2*thetaR**2)",
            'PQ':'PS+decay_integral(c,3,Ts)*(flatten.Ev2*thetaS**2/2)',
            'PT':"PQ+steep.outfull['pressure']*(flatten.Ev2*thetaQ**2)",
            'Ptail':'PT+decay_integral(c,1+native.delta,W)*(flatten.Ev2*thetaT**2/2)',
            'theta_base':'thetaT*c.exp(-bh*W-logone)',
            'pressure_scale':'flatten.Ev2*theta_base**2'},
        ('compliant_steep_waiting_C4','kernels'):{
            'self.kernel_cache[key]':'transition_kernels(self.ctx,t,self.mu,self.delta,kind,self.cells)'}}
    for (stem,method),values in specs.items():
        for target,expression in values.items():
            binding(stem,method,target,expression,allow_other_assignments=(target=='bracket'))
            bindings[stem+'.'+method+':'+target]=True
    return_binding('compliant_outer_angular_candidate','f','c.exp(-1/y**2)',allow_zero=True)
    bindings['preheat_phi_original_return']=True
    class_assignment('current_steep_waiting_source','CurrentSteepWaitingC4','__init__','self.wait','box(f.angular.waiting)')
    class_assignment('current_power_angular_source','CurrentPowerAngularC4','__init__','self.weights',
        '{k:box(v) for k,v in future.repair.weights.items()}')
    class_assignment('current_power_angular_source','CurrentPowerAngularC4','__init__','self.normalization','self.flat.normalization')
    class_assignment('future_swirl_energy','CompliantFutureSwirlEnergy','__init__','self.Lrel','-30*self.params.log_mu')
    for target,expression in (('self.Ts','box(f.params.Ts)'),
        ('self.infull',"self.kernels(1,'in')"),('self.outfull',"self.kernels(1,'out')")):
        class_assignment('current_steep_waiting_source','CurrentSteepWaitingC4','__init__',target,expression)
    if exact.future.params is not exact.repair.params or exact.companion.heat.steep.future.params is not exact.future.params:
        raise ValueError('Transition parameter source alias differs')
    return_binding('compliant_outer_angular_repair','inverse','((u*self.q-v)/self.det,(v*self.p-u)/self.det)')
    # Actual support endpoint theorem identifies F(0)=1 and full past A;
    # it does not use numerical overlap or discard a nonzero past history.
    endpoint=actual_angular_zero_future_bindings()
    if not endpoint['passed']:raise ValueError('Current full-support angular endpoint theorem missing')
    proof=json.loads((HERE/(PREFIX+'current_exact_repair_branch_check.json')).read_bytes())
    if not proof['current_exact_native_repair_source_proof']['same_equations_plus_common_ball_uniqueness_bind_both_paths']:
        raise ValueError('Current exact affine input/unique branch theorem required')
    heat_binding=exact.companion.current_bindings['shared_exact_Gamma_future_binding']
    if not heat_binding['passed'] or not all(heat_binding['exact_full_Gamma_future_identities'].values()):
        raise ValueError('Same current full Gamma angular future functions required')
    kernel_binding=exact.companion.source_owner.before.bindings['exact_current_kernel_and_parameter_source_bindings']
    if not kernel_binding['passed'] or not all(kernel_binding['actual_shared_kernel_functional_identities'].values()):
        raise ValueError('Current original sigma/transition/J(1)=1/2 source theorem required')
    stress_receipt=json.loads((HERE/(PREFIX+'current_heat_pressure_stress_check.json')).read_bytes())
    if not stress_receipt['canonical_projected_full_Gamma_moment_stress_identities']['full_terminal_moment_stress_theorem_verified']:
        raise ValueError('Canonical full Gamma angular stress source theorem required')
    r,k,L,Ts,W,lone,eps,S=s.symbols('rate k L Ts W logone eps S',real=True)
    xf,xf0,A,Iin,Iout,Theta,J=s.symbols('Xf Xf0 A Iin Iout Theta J',real=True)
    scale=s.exp(-r*L);logtransfer=-(r+k)/2-k*W
    logmult=(r+k)/2+k*W-lone
    past=assignment('compliant_power_angular_C4','angular','pastA',{
        'dj':s.Symbol('dj'), 'self.rate':r,'center':s.Symbol('center'),"past['A']":s.Symbol('weightA')},augmented=True)
    dj,center,weightA=s.symbols('dj center weightA')
    zero('actual_current_full_bump_integrand',past,dj*s.exp(r*center)*weightA)
    # The two support centers and normalized integrals give the exact first
    # repair equation. The current source theorem already binds that equation.
    zero('two_bump_angular_equation',sum(d*s.exp(r*j)*weightA for d,j in zip(s.symbols('d1 d2'),(-3,-1))),
        weightA*(s.exp(-3*r)*s.Symbol('d1')+s.exp(-r)*s.Symbol('d2')))
    # Replay the first equation of the actual inverse/quadratic system.
    # The nonlinear pressure term changes only its second right-hand side.
    p,q,u,v=s.symbols('p q u v');det=p*q-1
    nx=(u*q-v)/det;ny=(v*p-u)/det
    zero('actual_repair_inverse_preserves_first_equation',p*nx+ny,u)
    x,y,rhs_scaled,coefficient_scale=s.symbols('x y rhs_scaled coefficient_scale')
    bumped=weightA*coefficient_scale*(s.exp(-3*r)*x+s.exp(-r)*y)
    zero('actual_scaled_quadratic_first_equation_to_pastA',
        bumped.subs(y,rhs_scaled*s.exp(r)/weightA-s.exp(-2*r)*x),coefficient_scale*rhs_scaled)
    current_angular=assignment('compliant_current_angular_terminal_closure','current_forward_terminal','pastA',{
        'dj':dj,'r':r,'center':center,"weights['A']":weightA},augmented=True)
    zero('new_runtime_pastA_is_original_full_bump_integrand',current_angular,past)
    pp,B,D=s.symbols('pressure_rate weightB weightD')
    old_pastP=assignment('compliant_power_angular_C4','angular','pastP',{
        'dj':dj,"past['B']":B,"past['D']":D,'self.prate':pp,'center':center},augmented=True)
    new_pastP=assignment('compliant_current_angular_terminal_closure','current_forward_terminal','pastP',{
        'dj':dj,"weights['B']":B,"weights['D']":D,'prate':pp,'center':center},augmented=True)
    zero('new_runtime_pastP_is_original_full_bump_integrand',new_pastP,old_pastP)
    Ev2,fullP=s.symbols('Ev2 fullP')
    old_pressure_change=assignment('compliant_power_angular_C4','angular','pressure',{
        'pastP':fullP,'self.flatten.Ev2':Ev2,'self.prate':pp,'self.Lrel':L},augmented=True)
    new_pressure_change=assignment('compliant_current_angular_terminal_closure','current_forward_terminal','PR',{
        'pastP':fullP,'flatten.Ev2':Ev2,'prate':pp,'L':L},augmented=True)
    zero('new_runtime_pressure_bump_uses_original_absolute_units',new_pressure_change,old_pressure_change)
    XR=assignment('compliant_power_angular_C4','angular','Xbase',{
        'eq':1/r,"data['flatten_exit_X']":xf,'self.rate':r,'y':L})+A
    env={'XR':XR,"self.infull['angular']":Iin,'self.rate':r,'self.k':k,'self.Ts':Ts,
        "self.outfull['angular']":Iout}
    XS=assignment('compliant_current_steep_waiting_source','data','XS',env);env['XS']=XS
    XQ=assignment('compliant_current_steep_waiting_source','data','XQ',env);env['XQ']=XQ
    XT=assignment('compliant_current_steep_waiting_source','data','XT',env)
    Xtail=assignment('compliant_steep_waiting_C4','waiting','X',{"data['XT']":XT,'self.k':k,'t':W})
    XT0=XT.subs({xf:xf0,A:0});rawtail0=Xtail.subs({xf:xf0,A:0})
    zero('actual_current_XR_full_history',XR,1/r+(xf-1/r)*scale+A)
    zero('actual_corrected_vs_raw_waiting_reference',Xtail-rawtail0,s.exp(logtransfer)*(scale*(xf-xf0)+A))
    # Repair is evaluated on the exact current Xv* source, including at Z=0.
    repaired_A=scale*(xf0-xf)+s.exp(logmult)*S*Theta
    zero('unique_runtime_coefficients_realize_corrected_A',
        (coefficient_scale*rhs_scaled).subs(rhs_scaled,repaired_A/coefficient_scale),repaired_A)
    zero('angular_heat_multiplier_inverse_current_transfer',(1-eps)*s.exp(logtransfer+logmult).subs(lone,s.log(1-eps)),1)
    zero('current_corrected_vs_raw0_heat_increment',
        ((1-eps)*(Xtail.subs(A,repaired_A)-rawtail0)).subs(lone,s.log(1-eps)),S*Theta)
    decay=eps*(1/k+J)/((1-eps)*(XT0-1/k))
    rawroot=1/k+(XT0-1/k)*decay
    zero('raw_prescribed_waiting_target',(1-eps)*rawroot,1/k+eps*J)
    corrected=rawroot+s.exp(logtransfer)*(scale*(xf-xf0)+repaired_A)
    heat0=assignment('compliant_collar_Gamma_C4','collar_tails','A',{
        'one':1,'self.k':k,'angular':Theta,'self.S':S,"atoms['JW']":J,'self.eps':eps,'t':0})
    residual=((1-eps)*corrected-heat0).subs(lone,s.log(1-eps))
    zero('current_actual_Dtheta_identically_zero',residual)
    # Native vs repair transition kernels are the same exact integrands;
    # the directed quadrature decompositions merely bound them differently.
    v,aa,bb,dt,jv=s.symbols('v left right dt Jv',real=True)
    zero('unit_in_exponential_weight_primitive',s.diff(s.exp(r*v)/r,v),s.exp(r*v))
    zero('current_in_angular_kernel',s.exp(r*(v-jv)),s.exp(r*v)*s.exp(-r*jv))
    zero('current_out_angular_kernel',s.exp(k*jv),s.exp(k*jv))
    zero('same_preheat_phi_function',-1/((3-v)/2)**2,-4/(3-v)**2)
    sig,phi=s.symbols('sigma phi')
    # The second assignment intersects a directed enclosure with[0,1];
    # it does not change this AST-bound defining integrand.
    oldJ=1-sig+sig*phi
    zero('same_current_collar_W',oldJ,1-sig+sig*phi)
    weighted=assignment('compliant_outer_angular_candidate','collar_preheat_integral','result',
        {'k':k,'b':bb,'a':aa,'bracket':oldJ},augmented=True)
    zero('old_J_weight_is_same_exact_measure',weighted,oldJ*(s.exp(k*bb)-s.exp(k*aa))/k)
    currentJ=assignment('compliant_collar_Gamma_C4','collar_tails',"atoms['JW']",{'ds':dt,'self.k':k,'v':v,'W':oldJ},augmented=True)
    zero('current_J_same_integrand',currentJ,dt*s.exp(k*v)*oldJ)
    # Integrating identical functions over[0,3] identifies J, not its boxes.
    z=s.symbols('Z',real=True)
    substituted=residual.subs({xf:s.Function('Xf')(z),Theta:s.Function('Theta')(z)})
    for n in range(6):zero('current_Dtheta_axial_derivative_'+str(n),s.diff(substituted,z,n))
    return dict(actual_current_production_AST_bindings=bindings,identities=identities,
        current_exact_native_angular_branch_identified_by_uniform_uniqueness=True,
        native_live_C4_C5_functions_not_serialized_jet_equality=True,
        same_exact_preheat_collarJ_and_current_JW_function=True,
        same_current_exact_inverse_radius_and_full_Gamma_Theta=True,
        admitted_current_sigma_and_transition_kernel_identities_consumed=len(kernel_binding['actual_shared_kernel_functional_identities']),
        admitted_canonical_full_Gamma_angular_stress_theorem_consumed=True,
        actual_quadratic_first_equation_to_runtime_full_bump_replayed=True,
        rebound_exact_C5_coefficients_waiting_radius_and_pressure_adapter_AST_bound=True,
        parameter_only_transition_kernels_bound_to_same_rebound_parameter_source=True,
        original_pressure_bump_integrand_and_absolute_unit_conversion_replayed=True,
        old_companion_evaluate_not_used_for_recovered_runtime=True,
        support_endpoint_proof=endpoint,
        raw_waiting_root_not_confused_with_corrected_actual_XT=True,
        nonzero_heat_repair_at_Z0_retained=True,
        exact_function_identity_implies_whole_Z_axial5=True,
        interval_overlap_not_used_to_zero_constant=True,
        angular_transfer_log=str(logtransfer),angular_heat_multiplier_log=str(logmult),
        full_future_or_selected_pulse_installation_claim=False,passed=True)


class CurrentAngularTerminalClosure:
    @source_precision
    def __init__(self,exact=None,collar=None,require_checked=True):
        self.exact=exact if exact is not None else CurrentExactRepairBranch()
        if not self.exact.acceptance_loaded:raise ValueError('Checked current exact repair required')
        self.companion=self.exact.companion;self.heat=current_heat_shape_view(self.exact);self.ctx=self.heat.ctx
        self.collar=collar if collar is not None else CurrentCollarStressMixedC4(companion=self.companion)
        if not self.collar.acceptance_loaded or self.collar.companion is not self.companion:
            raise ValueError('Checked collar stress and exact repair need the same current source graph')
        self.family=self.exact.family;self.source=self.exact.source;self.datum_sha=self.exact.datum_sha
        self.proof=current_angular_source_proof(self.exact)
        self.hashes=dict(self.exact.hashes);self.hashes.update(self.collar.hashes)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.cache={};self.terminals={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[1]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in SCOPES):
                raise ValueError('Current angular closure admission source/scope differs')
            if receipt['current_actual_angular_source_identity_proof']!=encode(pack(self.proof)):
                raise ValueError('Current angular source proof changed')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def terminal_constants(self,Z):
        c=self.ctx;z=c.mpf(Z);key=z._mpi_
        if key not in self.terminals:
            terminal=current_forward_terminal(self.exact,z);tail0=self.heat.collar_tails(z,0)
            pressure_infinity=terminal['Ptail']+tail0['remaining_pressure_in_Rtail_units']*terminal['pressure_scale']
            self.terminals[key]=dict(current_repaired_forward_terminal=terminal,
                current_forward_angular_subtraction=terminal['Xtail']*(1-self.heat.eps)-tail0['angular_numerator'],
                current_full_collar_future_pressure=tail0['remaining_pressure_in_Rtail_units'],
                pressure_infinity=pressure_infinity,
                exact_repair_coefficients_waiting_and_radius_used=True)
        return self.terminals[key]

    @source_precision
    def evaluate(self,chart,Z,t):
        if chart not in ('heat_collar','heat_exterior'):raise ValueError('Current original heat chart required')
        c=self.ctx;z=c.mpf(Z);t=c.mpf(t);key=(chart,z._mpi_,t._mpi_)
        if endpoints(z)[0]<-1 or endpoints(z)[1]>1 or endpoints(t)[0]<0 or (chart=='heat_collar' and endpoints(t)[1]>3) or (chart=='heat_exterior' and endpoints(t)[0]<3):
            raise ValueError('Original heat collar/exterior domain required')
        if key in self.cache:return self.cache[key]
        constants=self.terminal_constants(z);terminal=constants['current_repaired_forward_terminal']
        Cp=constants['pressure_infinity'];zero=IntervalTaylor.constant(c,0,5)
        if chart=='heat_collar':
            shape=shape_radial5(self.heat,z,t);tails=self.heat.collar_tails(z,t)
            defects=collar_defect_rows(self.heat,shape,tails,t)
            canonical=collar_stress_mixed4(self.heat,shape,defects,z,t)
            angular=tails['angular_numerator']/shape['K_rows'][0]
            remaining=tails['remaining_pressure_in_Rtail_units']
        else:
            shape=self.heat.local_Gamma(z,t);angular=shape['angular_numerator']/shape['K_rows'][0]
            remaining=shape['pressure_numerator']*c.exp(-self.heat.prate*t)
            canonical=dict(theta=[zero]*5,axial=[zero]*5)
        extra=constant_stress_rows(self.heat,z,t,zero,Cp,4)
        pressure=collar_pressure_rows(shape['K_rows'][:5],remaining,self.heat.prate,terminal['pressure_scale'],t)
        pressure[0]+=Cp
        out=dict(chart=chart,Z=z,coordinate=t,actual_angular_Taylor_after_source_closure=angular,
            actual_Dtheta_Taylor=zero,source_proof_eliminates_only_angular_constant=True,
            original_forward_subtraction_enclosure=constants['current_forward_angular_subtraction'],
            current_repaired_forward_terminal=terminal,
            current_full_collar_future_pressure=constants['current_full_collar_future_pressure'],
            actual_pressure_infinity_retained=Cp,
            actual_stress_factored_mixed4=dict(theta_Qtheta=mixed(canonical['theta'],4),
                axial_Qz=mixed(canonical['axial'],4),axial_pressure_Qpressure=mixed(extra['axial_pressure_constant'],4)),
            stable_current_absolute_pressure_mixed4=mixed(pressure,4),
            positive_source_stress_factors=dict(theta_Qtheta='sqrt(R/2)*B',axial_Qz='sqrt(R/2)*B^2',
                axial_pressure_Qpressure='sqrt(R/2)*Pstar^2',B='Ev0*theta_base*exp(-(1+delta)*t/2)',
                Ev0='Pstar*U*exp(-13/(2mu)-13)',R='Rtail*exp(t)',
                exact_logRtail_terms=self.exact.repair.heat.logradius_terms,
                exact_Ev0_squared_over_Pstar_squared_logs=self.exact.flatten.logEv2_parts,
                theta_base_from_current_exact_waiting=True,
                physical_tensor_prefactor='nu*lambda^(-2-delta); global transfer remains open'),
            current_exact_repair_runtime_used=True,
            canonical_Gamma_angular_stress_zero=chart=='heat_exterior',
            actual_pressure_constant_not_eliminated=True)
        out.update({gate:self.acceptance_loaded for gate in GATES});out.update({gate:False for gate in SCOPES})
        self.cache[key]=out;return out


@source_precision
def run(field=None):
    field=field if field is not None else CurrentAngularTerminalClosure(require_checked=False)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_actual_angular_source_identity_proof=field.proof,
        current_angular_terminal_views={name:field.evaluate(chart,Z,t) for name,(chart,Z,t) in VIEWS.items()},
        **dict.fromkeys(GATES+SCOPES,False),input_hashes=field.hashes)
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current corrected angular history closes Dtheta; pressure constant retained',flush=True)
    return result


if __name__=='__main__':run()
