"""Current second quadratic moment equation cancels the actual heat loss.

This replays the pressure units on the checked exact-repair heat runtime.
It reduces Cp to the ORIGINAL raw-preheat pressure constant. Identification
of that raw cumulative history with the analytic datum is a separate input,
and zero Cp is deliberately not admitted by this balance alone.
"""
import ast,json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_angular_terminal_closure import (
    CurrentAngularTerminalClosure,HERE,PREFIX,sha,function,binding,pack,encode,
    endpoints,IntervalTaylor,current_forward_terminal,return_binding)
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_pressure_terminal_balance.json'
RECEIPT=PREFIX+'current_pressure_terminal_balance_check.json'
GATES=('current_pressure_quadratic_to_actual_bump_function_certified',
       'current_pressure_bump_full_Gamma_heat_loss_balance_certified',
       'current_pressure_terminal_reduced_to_raw_preheat_constant')
SCOPES=('current_heat_pressure_terminal_constant_eliminated',
    'current_heat_terminal_constants_eliminated','heat_exterior_stress_identity_certified',
    'current_exact_repair_installed_in_all_physical_charts','global_completed_tensor_admissibility',
    'admissible_stress_lift_constructed','full_background_NS_validation',
    'physical_energy_integral_certified','independently_bounded_flat_remainder',
    'full_point_physical_field_evaluation','full_cartesian_vector_derivatives_certified','temporal_recursion')
VIEWS={'whole_Z':[-1,1],'axis':0,'fresh':'.523','endpoint':1}


def current_pressure_balance_proof(angular):
    """Source algebra on arbitrary smooth Z functions, never box equality."""
    if not angular.acceptance_loaded:raise ValueError('Checked current angular runtime source required')
    exact=angular.exact;proofs={};bindings={}
    def zero(name,left,right=0):
        if s.simplify(s.expand_power_exp(left-right))!=0:raise ArithmeticError('Current pressure identity differs: '+name)
        proofs[name]=True
    specs={
        ('compliant_outer_angular_repair','__init__'):{
            'self.q':'c.exp(-2*self.prate)','self.k':"self.weights['D']/(2*self.weights['B'])"},
        ('compliant_outer_angular_repair','coefficients'):{
            's':"rhs['s_scaled']",
            'b2':"s*(c.exp(-3*self.prate)/self.weights['B'])",
            'nonlinear':'self.k*self.scale*(x*x+self.q*y*y)',
            '(nx, ny)':'self.inverse(b1[0],b2[0]-nonlinear)'},
        ('compliant_current_exact_repair_branch','replay_repair'):{
            'out.logEtail_over_Erel':"-c.mpf('1.5')*(2+out.params.Ts)+(out.rate+angular.restore_rate)/2-(c.mpf('.5')+out.delta/2)*angular.waiting",
            'out.log_pressure_multiplier':'2*out.logEtail_over_Erel-2*angular.waiting_logone+c.ln(out.delta/2)',
            'out.pressure_heat_over_scale_cap':'c.exp(out.log_pressure_multiplier+out.strong_logS_cap-out.logscale)'},
        ('compliant_outer_angular_repair','defects'):{
            'exact_s_definition':"'sH=exp(log_pressure_multiplier)*S*Pressure_hat(Z)'"},
        ('compliant_angular_high_jets','coefficients'):{
            'r':'self.repair',
            'oldrhs':"old['defects']",
            'sh':"scaled_heat(heat['pressure'], r.pressure_heat_over_scale_cap, oldrhs['s_scaled'])",
            'b2':"sh*(c.exp(-3*r.prate)/r.weights['B'])",
            'physical':"[v*r.scale for v in out['scaled_coefficient_Taylor']]"},
        ('compliant_fifth_axial_jets','angular'):{
            'r':'self.angular4.repair','old':'self.angular4.coefficients(endpoints(Z))',
            'sh5':"heat['pressure'][5]*c.mpf([0,endpoints(r.pressure_heat_over_scale_cap)[1]])",
            'sh':"append_fifth(c,old['s_scaled_Taylor'],sh5)",
            'prior':"[copy_jet(c,j) for j in old['scaled_coefficient_Taylor']]",
            '(coeffs, det)':"angular_fifth(c,prior,rhs[5]*(c.exp(box(r.rate))/box(r.weights['A'])),sh[5]*(c.exp(-3*box(r.prate))/box(r.weights['B'])),box(r.p),box(r.q),box(r.k)*box(r.scale))",
            'physical':'[j*box(r.scale) for j in coeffs]'},
        ('compliant_current_angular_terminal_closure','current_forward_terminal'):{
            'source':'exact.angular5.angular(z)',
            'weights':'{name:box(v) for name,v in repair.weights.items()}',
            'pressure_scale':'flatten.Ev2*theta_base**2',
            'theta_base':'thetaT*c.exp(-bh*W-logone)'},
        ('compliant_collar_Gamma_C4','collar_tails'):{
            'square':'D*(pre*2-D*(self.a*self.S))',
            'P':"one*(c.exp(-self.prate*t)/(2*self.prate)-self.eps*atoms['PW']+self.eps**2*atoms['PW2']/2)-pressure*(self.a*self.S)"}}
    for (stem,method),rows in specs.items():
        for target,expression in rows.items():
            binding(stem,method,target,expression)
            bindings[stem+'.'+method+':'+target]=True
    return_binding('compliant_outer_angular_repair','inverse','((u*self.q-v)/self.det,(v*self.p-u)/self.det)')
    # The checked native branch and C5 prefix prove one smooth solution of
    # precisely these equations, not six independently selectable numbers.
    receipt=json.loads((HERE/(PREFIX+'current_exact_repair_branch_check.json')).read_bytes())
    shared=receipt['current_exact_native_repair_source_proof']
    if not shared['same_equations_plus_common_ball_uniqueness_bind_both_paths']:
        raise ValueError('Current exact quadratic branch identification required')
    if exact.angular4.repair is not exact.future.repair:raise ValueError('Pressure and angular branches differ')
    heat=exact.companion.current_bindings['shared_exact_Gamma_future_binding']
    if not heat['passed'] or not all(heat['exact_full_Gamma_future_identities'].values()):
        raise ValueError('Current pressure repair must use the same complete Gamma future')
    x,y,p,q,b1,b2,K,scale=s.symbols('x y p q b1 b2 K scale',real=True)
    v=b2-K*scale*(x*x+q*y*y);det=p*q-1
    nx=(b1*q-v)/det;ny=(v*p-b1)/det
    zero('actual_inverse_second_equation',nx+q*ny,v)
    B,D,pp,rhs=s.symbols('B D pressure_rate rhs_pressure_scaled',real=True)
    bump=B*(s.exp(3*pp)*scale*x+s.exp(pp)*scale*y)+D/2*(s.exp(3*pp)*scale**2*x*x+s.exp(pp)*scale**2*y*y)
    normalized=B*s.exp(3*pp)*scale*(x+q*y+K*scale*(x*x+q*y*y))
    zero('original_pressure_bump_is_scaled_second_quadratic_equation',
        normalized.subs({q:s.exp(-2*pp),K:D/(2*B)}),bump)
    zero('actual_unique_quadratic_branch_realizes_pressure_rhs',
        (B*s.exp(3*pp)*scale*b2).subs(b2,rhs*s.exp(-3*pp)/B),scale*rhs)
    dj,center=s.symbols('dj center',real=True)
    old=assignment('compliant_power_angular_C4','angular','pastP',{
        'dj':dj,"past['B']":B,"past['D']":D,'self.prate':pp,'center':center},augmented=True)
    new=assignment('compliant_current_angular_terminal_closure','current_forward_terminal','pastP',{
        'dj':dj,"weights['B']":B,"weights['D']":D,'prate':pp,'center':center},augmented=True)
    zero('new_current_pastP_uses_original_supported_pressure_density',new,old)
    mu,a,Ts,W,lone,L,Ev2=s.symbols('mu a Ts W logone L Ev2',real=True)
    r=1-mu;k=1-a;bp=s.Rational(1,2)+mu;bh=s.Rational(1,2)+a
    logR=-100*bp-bp*L-s.log(2)
    logT=logR-bp-r/2-s.Rational(3,2)*Ts-s.Rational(3,2)+k/2
    logbase=logT-bh*W-lone
    ell=-s.Rational(3,2)*(2+Ts)+(r+k)/2-bh*W
    zero('recomputed_theta_base_over_Erel_log_units',logbase-logR,ell-lone)
    Erel2=Ev2*s.exp(-(1+2*mu)*(100+L))/4
    pressure_scale=Ev2*s.exp(2*logbase)
    zero('new_pressure_scale_exact_Erel_units',pressure_scale,Erel2*s.exp(2*ell-2*lone))
    S,Ph,eps,PW,PW2=s.symbols('S Ph epsilon PW PW2',real=True)
    sh=s.exp(2*ell-2*lone)*a*S*Ph
    logscale,logS,logmultiplier=s.symbols('logscale logS log_pressure_multiplier',real=True)
    # The original defect definition specifies one exact function. The
    # cap recipes above bound its C4/C5 coefficients; they never define
    # its value. All orders inherit the same factor and unique branch.
    scaled_rhs_function=s.exp(logmultiplier+logS-logscale)*Ph
    zero('current_exact_scaled_pressure_rhs_times_scale_is_original_sH',
        s.exp(logscale)*scaled_rhs_function,s.exp(logmultiplier)*s.exp(logS)*Ph)
    zero('current_pressure_rhs_exact_log_multiplier_and_inverse_radius',
        (s.exp(logmultiplier)*s.exp(logS)*Ph).subs({logmultiplier:2*ell-2*lone+s.log(a),logS:s.log(S)}),sh)
    # Substitute the ORIGINAL second equation of the unique branch, not
    # a numerically chosen pair. This explicitly connects runtime pastP
    # to the heat RHS before using the cancellation below.
    # Replay its named RHS after the preceding exact factorization.
    zero('runtime_pastP_equals_pressure_heat_rhs_from_unique_second_equation',
        (B*s.exp(3*pp)*scale*b2).subs(b2,sh*s.exp(-3*pp)/(B*scale)),sh)
    loss=pressure_scale*a*S*Ph
    zero('actual_pressure_bump_equals_same_complete_Gamma_pressure_loss',Erel2*sh,loss)
    prate=1+2*a
    pre=1/(2*prate)-eps*PW+eps**2*PW2/2
    actual=assignment('compliant_collar_Gamma_C4','collar_tails','P',{
        'one':1,'t':0,'self.prate':prate,'self.eps':eps,
        "atoms['PW']":PW,"atoms['PW2']":PW2,'pressure':Ph,'self.a':a,'self.S':S})
    zero('same_full_future_pressure_preheat_minus_Gamma_loss',actual,pre-a*S*Ph)
    rawtail,p0,rawmass=s.symbols('raw_Ptail P0 raw_forward_mass',real=True)
    Cp=rawtail+Erel2*sh+pressure_scale*actual
    raw_constant=rawtail+pressure_scale*pre
    zero('current_actual_Cp_equals_original_raw_preheat_constant',Cp,raw_constant)
    z=s.symbols('Z',real=True)
    functional=(Cp-raw_constant).subs(Ph,s.Function('Pressure_hat')(z))
    for n in range(6):zero('pressure_balance_axial_derivative_'+str(n),s.diff(functional,z,n))
    return dict(actual_current_production_AST_bindings=bindings,identities=proofs,
        same_current_unique_quadratic_C5_function_consumed=True,
        runtime_pastP_to_exact_pressure_heat_rhs_replayed=True,
        current_C4_prefix_and_C5_pressure_RHS_same_exact_function_bound=True,
        scaled_rhs_function_definition='exp(log_pressure_multiplier + logS_exact - logscale) * Pressure_hat(Z)',
        RHS_coefficient_caps_only_enclose_this_function=True,
        same_complete_pressure_Gamma_future_integrand_consumed=True,
        actual_Pstar_squared_units_and_both_half_factors_retained=True,
        pressure_multiplier_log=str(2*ell-2*lone+s.log(a)),
        Cp_reduced_to='P0 + complete_original_raw_preheat_pressure_integral',
        raw_preheat_datum_to_forward_integral_identification_not_inferred=True,
        no_post_propagation_pressure_patch=True,zero_Cp_not_admitted=True,passed=True)


class CurrentPressureTerminalBalance:
    @source_precision
    def __init__(self,angular=None,require_checked=True):
        self.angular=angular if angular is not None else CurrentAngularTerminalClosure()
        self.proof=current_pressure_balance_proof(self.angular)
        self.exact=self.angular.exact;self.heat=self.angular.heat;self.ctx=self.angular.ctx
        self.family=self.angular.family;self.source=self.angular.source;self.datum_sha=self.angular.datum_sha
        self.hashes=dict(self.angular.hashes);self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.cache={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in SCOPES):
                raise ValueError('Current pressure balance admission source/scope differs')
            if receipt['current_pressure_balance_source_proof']!=encode(pack(self.proof)):
                raise ValueError('Current pressure source proof changed')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def evaluate(self,Z):
        c=self.ctx;z=c.mpf(Z);key=z._mpi_
        if key in self.cache:return self.cache[key]
        constants=self.angular.terminal_constants(z);forward=constants['current_repaired_forward_terminal']
        future=self.heat.collar_tails(z,0);atoms=future['separate_epsilon_atoms']
        pre=IntervalTaylor.constant(c,1,5)*(1/(2*self.heat.prate)-self.heat.eps*atoms['PW']+self.heat.eps**2*atoms['PW2']/2)
        Erel2=self.exact.flatten.Ev2*c.exp(-(1+2*self.heat.mu)*(100+c.mpf(endpoints(-30*self.exact.repair.params.log_mu))))/4
        bump=forward['pastP']*Erel2
        loss=future['scaled_full_future_Gamma_defects']['pressure']*(self.heat.a*self.heat.S*forward['pressure_scale'])
        rawtail=forward['Ptail']-bump
        rawconstant=rawtail+pre*forward['pressure_scale']
        out=dict(Z=z,current_repaired_forward_terminal=forward,
            actual_pressure_infinity_retained=constants['pressure_infinity'],
            same_raw_preheat_pressure_constant_enclosure=rawconstant,
            actual_absolute_pressure_bump_enclosure=bump,actual_absolute_Gamma_pressure_loss_enclosure=loss,
            source_proved_bump_minus_heat_loss_Taylor=IntervalTaylor.constant(c,0,5),
            raw_preheat_remaining_pressure_in_Rtail_units=pre,
            raw_pressure_constant_retained_not_zeroed=True,
            exact_positive_Erel2_source_log_terms=dict(
                native_Ev2=self.exact.flatten.logEv2_parts,
                relative_factor=-(1+2*self.heat.mu)*(100+c.mpf(endpoints(-30*self.exact.repair.params.log_mu)))-c.ln(4)),
            output_enclosures_not_selected_values=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(SCOPES,False))
        self.cache[key]=out;return out


@source_precision
def run(field=None):
    field=field if field is not None else CurrentPressureTerminalBalance(require_checked=False)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_pressure_balance_source_proof=field.proof,
        current_pressure_balance_views={name:field.evaluate(Z) for name,Z in VIEWS.items()},
        input_hashes=field.hashes,**dict.fromkeys(GATES+SCOPES,False))
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current pressure quadratic bump cancels complete Gamma loss; raw-preheat datum bridge still required',flush=True)
    return result


if __name__=='__main__':run()
