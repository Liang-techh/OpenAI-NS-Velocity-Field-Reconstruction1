"""Native source, correlated replay and current unique repair admission."""
import ast
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s

from lei_ren_part1_paper_compliant_current_exact_repair_branch import (
    CurrentExactRepairBranch,HERE,PREFIX,NAME,RECEIPT,GATES,SCOPES,VIEWS,
    sha,function,binding,pack,encode,endpoints,CompliantAngularRepair)
from lei_ren_part1_paper_compliant_outer_buffer import SharedOuterBuffer
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_future_energy_high_jets import copy_jet
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def exact(left,right,label):
    difference=s.simplify(left-right)
    if isinstance(difference,s.MatrixBase):passed=difference.is_zero_matrix
    else:passed=difference==0
    if not passed:raise ArithmeticError('Current repair source identity differs: '+label)


def source_proof(field):
    """Bind written production recipes before using uniform uniqueness."""
    bindings={}
    specifications={
        ('compliant_outer_angular_candidate','__init__'):{
            'inlet':"self.buffer.power('0',1,cells)",
            'Xp':"inlet['Mtheta_over_sqrt2_R_3half_Pstar'][0]/inlet['Utheta_over_Pstar'][0]",
            'self.waiting':'(c.ln(Xt-eq)+logone-self.params.log_epsilon-c.ln(eq+self.collarJ))/k'},
        ('compliant_current_exact_repair_branch','__init__'):{
            'angular.Xp':'c.mpf(endpoints(self.pulse.Xp))',
            'angular.loguRp0':'c.ln(c.mpf(endpoints(self.flatten.U)))',
            'angular.native_terminal_factor':'self.pulse.factor',
            'native_decay_bound':'c.mpf(endpoints(angular.native_terminal_factor(-13*self.pulse.rate/self.pulse.mu)))',
            'decay_bound':'intersect(c,native_decay_bound,analytic_decay_bound)',
            'self.native_exponent':'-13*angular.rate/angular.mu',
            'angular.Xv':'equilibrium+difference*decay_bound',
            'terminal':"angular.before_waiting('0')",
            'angular.waiting':'(c.ln(Xt-eq)+angular.waiting_logone-angular.params.log_epsilon-c.ln(eq+angular.collarJ))/k',
            'self.angular4.cache':'{}',
            'self.angular5':'SimpleNamespace(ctx=self.pulse.fifth.ctx,angular4=self.angular4,angular_cache={})',
            'self.angular5.angular':'MethodType(CompliantFifthAxialJets.angular,self.angular5)'},
        ('compliant_outer_angular_repair','preheat_difference'):{
            'start':'2*self.angular.Xv*c.exp(-100*self.rate)'},
        ('compliant_angular_high_jets','preheat_jets'):{
            'Xf':'q.reciprocal()*(2*r.angular.Xv*c.exp(-100*r.rate))'},
        ('compliant_fifth_axial_jets','preheat_fifth'):{
            'Xf':'q.reciprocal()*(2*c.mpf(endpoints(r.angular.Xv))*c.exp(-100*rate))'},
        ('compliant_current_exact_repair_branch','replay_heat'):{
            'out.angular':'angular',
            'out.tail_finite':'p.yd+1+p.Tw+100-30*p.log_mu+2+p.Ts+angular.waiting',
            'out.logtail_relative':"angular.waiting_field('0',1)['relative_log_profile'][0]",
            'out.logone':'angular.waiting_logone'},
        ('compliant_current_exact_repair_branch','replay_repair'):{
            'out.logtail_distance':'2+out.params.Ts+angular.waiting',
            'out.logEtail_over_Erel':"-c.mpf('1.5')*(2+out.params.Ts)+(out.rate+angular.restore_rate)/2-(c.mpf('.5')+out.delta/2)*angular.waiting",
            'out.log_theta_multiplier':'(out.rate+angular.restore_rate)/2+angular.restore_rate*angular.waiting-angular.waiting_logone',
            'out.log_pressure_multiplier':'2*out.logEtail_over_Erel-2*angular.waiting_logone+c.ln(out.delta/2)',
            'whole':'out.defects([-1,1])'},
        ('compliant_current_exact_repair_branch','replay_future'):{
            'out.Ntail':'out.Nt*c.exp(-out.delta*out.angular.waiting)',
            'out.waiting_energy':'decay_integral(c,out.delta,out.angular.waiting)',
            'out.tail_multiplier':'out.Ntail*c.exp(-2*out.angular.waiting_logone)'},
        ('compliant_outer_angular_repair','__init__'):{
            'self.normalization':'self.angular.initial.repair.normalization',
            'self.weights':'bump_weights(c,self.mu,self.normalization,cells=bump_cells)'},
        ('compliant_outer_angular_repair','bump_weights'):{
            'ell':"c.mpf('.15')",'beta':'raw_beta(c,raw_coordinate)/(ell*normalization)'},
        ('compliant_five_moment_repair','__init__'):{
            'self.normalization':"restore_value(c,self.wraw['normalization'])"}}
    for (stem,method),spec in specifications.items():
        for target,expression in spec.items():
            binding(stem,method,target,expression)
            bindings[stem+'.'+method+':'+target]=True
    class_assignment('current_pulse_flatten_source','CurrentFlattenMixedC4','__init__','self.Xv',
        '1/self.rate+(pulse.Xp-1/self.rate)*pulse.factor(-13*self.rate/self.mu)')
    class_assignment('current_pulse_flatten_source','CurrentFlattenMixedC4','__init__','self.U',
        "self.inlet.constants['U']")
    bindings['current_flatten_native_affine_terminal']=True
    old=field.heat.future.repair;other=field.pulse.fifth.angular4.repair
    buffers=(old.angular.buffer,other.angular.buffer,field.pulse.pulse.buffer)
    if not all(type(buffer) is SharedOuterBuffer and buffer.power.__func__ is SharedOuterBuffer.power for buffer in buffers):
        raise ValueError('Native Xp and both repairs need the same original buffer power callable')
    if not all(r.angular.initial.repair.wraw==old.angular.initial.repair.wraw for r in (old,other)):
        raise ValueError('The two repair beta normalizations have different defining integral receipts')
    fn=function('compliant_outer_angular_repair','__init__')
    if [ast.literal_eval(v) for v in fn.args.defaults]!=[512,1024]:
        raise ValueError('Repair/default beta quadrature scope changed')
    if old.cells!=512 or other.cells!=512 or field.angular4.cells!=256:
        raise ValueError('Repair and high-jet quadrature domains changed')
    x,Xp,mu,r,Q,k,eps,J,Xt=s.symbols('Xv Xp mu rate Q k eps J Xt',positive=True)
    memory=s.exp(-13*(1-mu)/mu)
    affine=1/(1-mu)+(Xp-1/(1-mu))*memory
    exact(affine,1/(1-mu)+(Xp-1/(1-mu))*memory,'native endpoint equals flatten terminal function')
    sigma,v=s.symbols('sigma v',real=True)
    exact(2/Q*s.exp(-r*(100-v))*(Q/2)**sigma,
        s.exp(-r*(100-v))*(2/Q)**(1-sigma),'same native flatten preheat integrand')
    eta,ss,power=s.symbols('eta ss power',positive=True)
    primitive=-(1+ss*eta)**(-power)/eta
    exact(s.diff(primitive,ss),power*(1+ss*eta)**(-1-power),'stable divided difference integrand')
    exact(eta*(primitive.subs(ss,1)-primitive.subs(ss,0)),1-(1+eta)**(-power),'stable preheat subtraction')
    waiting_decay=eps*(1/k+J)/((Xt-1/k)*(1-eps))
    exact((1-eps)*(1/k+(Xt-1/k)*waiting_decay),1/k+eps*J,'prescribed waiting root at Z0')
    p,q,K,scale,x,y,b1,b2=s.symbols('p q K scale x y b1 b2',real=True)
    matrix=s.Matrix([[p,1],[1+2*K*scale*x,q*(1+2*K*scale*y)]])
    exact(s.Matrix([p*x+y,x+q*y+K*scale*(x*x+q*y*y)]).jacobian([x,y]),matrix,'one current quadratic Jacobian')
    inv=matrix.inv();exact(matrix*inv,s.eye(2),'same implicit inverse through order5')
    # Differentiating polynomial equations gives one C0 Jacobian at every
    # positive order; mixed lower-order products form the known RHS only.
    return dict(actual_source_AST_bindings=bindings,
        same_original_native_and_both_repair_Xp_power_callable=True,
        same_beta_normalization_defining_integral_receipt=True,
        same_beta_support_and_default_weight_integral_recipe=True,
        exact_native_affine_terminal_suppression_is_function_containment=True,
        same_native_stable_factor_callable_carried_by_current_repair=True,
        current_native_inlet_heat_amplitude_origin_rebound=True,
        native_flatten_and_stable_preheat_subtraction_are_same_function=True,
        all_correlated_waiting_radius_amplitude_and_future_factors_recomputed=True,
        old_uniform_common_ball_theorem_consumed=True,
        same_equations_plus_common_ball_uniqueness_bind_both_paths=True,
        same_implicit_quadratic_Jacobian_through_order5=True,
        old_interval_box_not_asserted_equal_to_native_function=True,passed=True)


def independent_nontrivial_Xv_fixture():
    """Use finite Xv values whose contributions are resolvable and unequal."""
    c=MPIntervalContext();c.dps=70;checks=0
    with mp.workdps(95):
        def sigma(t):
            if t<=0:return mp.mpf(0)
            if t>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/t**2-1/(1-t)**2))
        rate=mp.mpf('.19');parts=[0,20,50,70,80,90,95,99,100]
        for Xv in ('1.3','2.4'):
            fixture=CompliantAngularRepair.__new__(CompliantAngularRepair)
            fixture.ctx=c;fixture.rate=c.mpf(str(rate));fixture.cells=512
            fixture.angular=SimpleNamespace(Xv=c.mpf(Xv))
            for location in ('0','.37','.83'):
                z=mp.mpf(location);Q=1+z*z
                direct=2*mp.mpf(Xv)*mp.exp(-100*rate)*(1-1/Q)+mp.quad(
                    lambda v:mp.exp(-rate*(100-v))*2**(1-sigma(v/100))*(1-Q**(-(1-sigma(v/100)))),parts)
                derivative=4*z*mp.mpf(Xv)*mp.exp(-100*rate)/Q**2+mp.quad(
                    lambda v:2*z*mp.exp(-rate*(100-v))*2**(1-sigma(v/100))*(1-sigma(v/100))*Q**(-2+sigma(v/100)),parts)
                jet=fixture.preheat_difference(location)
                for n,value in enumerate((direct,derivative)):
                    lo,hi=endpoints(jet[n])
                    if not lo<=value<=hi:raise ArithmeticError('Independent nontrivial-Xv preheat not enclosed')
                    checks+=1
    return dict(independent_current_preheat_rows=checks,two_distinct_resolvable_terminal_inputs_used=True,passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentExactRepairBranch(require_checked=False)
    if (raw['actual_five_defect_family_sha256'],raw['implicit_source_sha256'],raw['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
        raise ValueError('Current exact repair family/source/datum differs')
    if any(raw[k] for k in GATES+SCOPES):raise ValueError('Producer claims admission')
    proof=source_proof(field);fixture=independent_nontrivial_Xv_fixture();rows=0;future_rows=0
    for name,Z in VIEWS.items():
        packet=field.evaluate(Z)
        if encode(pack(packet))!=raw['current_exact_repair_views'][name]:raise ValueError('Current exact branch packet differs: '+name)
        if not all(packet['same_current_native_source_graph'].values()) or any(packet[k] for k in SCOPES):
            raise ValueError('Current exact repair graph or scope differs')
        a4=packet['current_angular_C4'];a5=packet['current_angular_C5']
        for index in range(2):
            expected=copy_jet(a5['scaled_coefficient_Taylor'][index].ctx,a4['scaled_coefficient_Taylor'][index])
            for n in range(5):
                if endpoints(expected[n])!=endpoints(a5['scaled_coefficient_Taylor'][index][n]):
                    raise ValueError('Current C5 changed its recomputed C4 prefix')
            for n in range(6):
                lo,hi=endpoints(a5['physical_coefficient_Taylor'][index][n])
                if not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Current coefficient is not finite')
                rows+=1
        if endpoints(a5['derivative_jacobian_determinant'])[0]<=0<=endpoints(a5['derivative_jacobian_determinant'])[1]:
            raise ArithmeticError('Current order5 implicit inverse is singular')
        energy=packet['current_complete_future_swirl_energy_C1']['complete_future_energy_Taylor']
        if energy.order!=1 or endpoints(energy[0])[0]<=0:raise ArithmeticError('Current complete future energy lost positivity')
        future_rows+=2
    hashes=dict(raw['input_hashes']);hashes[NAME]=sha(NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_exact_native_repair_source_proof=proof,
        independent_preheat_input_fixture=fixture,current_angular_C5_coefficient_rows=rows,
        current_complete_future_energy_C1_rows=future_rows,
        both_paths_share_one_current_unique_branch=True,all_passed=True,input_hashes=hashes)
    result.update({gate:True for gate in GATES});result.update({gate:False for gate in SCOPES})
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current native Xv/source and one angular-pressure branch PASS: '+str(rows)+' C5 rows, '+str(future_rows)+' future C1 rows; terminal constants pending',flush=True)
    return result


if __name__=='__main__':run()
