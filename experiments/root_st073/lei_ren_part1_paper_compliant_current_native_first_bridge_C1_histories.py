"""Whole first-bridge C1 history covers across genuine cutoff/phase crossings.

The original loop and its unique periodic inverse define the functions.
Analytic active-branch and whole-period bounds enclose them without choosing
an inverse phase or deleting unresolved active intervals. These conservative
covers supply the real incoming correction at the phase1 bridge seam.
"""
from fractions import Fraction
import json
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_true_chart_C1_transfer as transfer

history=transfer.history;density=transfer.density;first=density.first;slow=density.slow
current=density.current;prior=transfer.prior;packets=transfer.packets;native=transfer.native
HERE,PREFIX,sha=transfer.HERE,transfer.PREFIX,transfer.sha;ep=transfer.ep;RATES=transfer.RATES
NAME=PREFIX+'current_native_first_bridge_C1_histories.json'
RECEIPT=PREFIX+'current_native_first_bridge_C1_histories_check.json'
GATE='current_original_whole_first_bridge_cutoff_phase_C1_history_covers_executed'
ZERO=(0,0);DZ=(0,1)


def absolute_upper(value):
    """A directed magnitude bound, explicitly not a defining field value."""
    if value.zero:return value
    c=value.ctx;magnitude=max(abs(v) for v in ep(value.coefficient))
    upper=ep(value.scale.evaluate()+c.ln(c.mpf(magnitude)))[1]
    return prior.ScaledEnclosure(prior.FormalScale(value.scale.bases,offset=c.mpf(upper)),1,value.ledger)


def symmetric_bound(value):
    cap=absolute_upper(value)
    return cap if cap.zero else prior.ScaledEnclosure(cap.scale,(-1,1),cap.ledger)


def positive_restriction(value,lower_log,upper_log):
    """Intersect a positive FUNCTION branch with its proved log bounds."""
    c=value.ctx;lo,hi=ep(value.coefficient)
    if hi<=0:return None
    lower=ep(c.mpf(lower_log))[0];upper=min(ep(c.mpf(upper_log))[1],ep(value.scale.evaluate()+c.ln(c.mpf(hi)))[1])
    if lo>0:lower=max(lower,ep(value.scale.evaluate()+c.ln(c.mpf(lo)))[0])
    if lower>upper:return None
    return prior.ScaledEnclosure(prior.FormalScale(value.scale.bases,offset=c.mpf((lower,upper))),1,value.ledger)


def whole_cutoff_C1(roots,eta_log,log_a_lower):
    """Bound body, transition and flat branches separately, then their union."""
    original_a=roots['a'][ZERO];c=original_a.ctx;scalar=original_a.scalar
    eta_log=c.mpf(eta_log)
    if ep(eta_log)[1]>ep(c.ln(c.mpf('.5')))[0]:raise ValueError('Original eta<=1/2 required for active constraints')
    a=positive_restriction(original_a,log_a_lower,c.ln(3))
    if a is None:
        # kappa>=a>3>=2+eta excludes the entire original active support.
        return dict(q=scalar(0),q_Z=scalar(0),a=None,t0=None,t0_Z=None,
            record=dict(active_support_empty_by_original_kappa_ge_a=True,branches=['exact_flat']))
    b=scalar(c.mpf(('-1.5','1.5')));a_Z=roots['a'][DZ];b_Z=roots['b'][DZ]
    t0=(-b).positive_divide(a,log_a_lower)
    t0_Z=(-b_Z-t0*a_Z).positive_divide(a,log_a_lower)
    eta=prior.ScaledEnclosure(prior.FormalScale(a.scale.bases,offset=eta_log),1,a.ledger)
    Delta=roots['kappa_minus2'][ZERO];gamma=eta*2-Delta;gamma_Z=-roots['kappa_minus2'][DZ]
    branches=[];values=[];derivatives=[]
    for name,lower,upper in (('body_Delta_le_zero',eta_log+c.ln(2),c.ln(3)),
        ('transition_zero_lt_Delta_lt_eta',eta_log,eta_log+c.ln(2))):
        active_gamma=positive_restriction(gamma,lower,upper)
        if active_gamma is None:continue
        ratio=active_gamma.positive_divide(a*2,log_a_lower+c.ln(2))
        root_lower=(lower-c.ln(6))/2
        root=current.nonnegative_sqrt(ratio).positive_intersection(root_lower)
        ratio_Z=(gamma_Z-ratio*a_Z*2).positive_divide(a*2,log_a_lower+c.ln(2))
        root_Z=ratio_Z.positive_divide(root*2,root_lower+c.ln(2))
        if name.startswith('body'):
            cutoff=scalar(1);cutoff_Z=scalar(0);argument=None
        else:
            # Only the transition is evaluated on[0,1]. At both flat ends,
            # original sigma derivatives vanish, so no seam distribution or
            # derivative boundary term is introduced.
            argument=c.mpf((0,1));coefficients=prior.sigma_jets(c,argument)
            cutoff=scalar(coefficients[0])
            cutoff_Z=scalar(coefficients[1])*(-roots['kappa_minus2'][DZ]).positive_divide(eta,eta_log)
        q=root*cutoff;q_Z=root_Z*cutoff+root*cutoff_Z
        values.append(q);derivatives.append(q_Z)
        branches.append(dict(branch=name,gamma_C0=active_gamma.record(),
            original_cutoff_argument=argument,q_C0=q.record(),q_Z=q_Z.record(),
            sigma_rows_are_Taylor_coefficients_first_order_factorial_is_one=True))
    values.append(scalar(0));derivatives.append(scalar(0))
    # Common positive arithmetic coordinates give a valid union including
    # exact flat zero. The endpoints bound functions; they never define q.
    q_upper=max((absolute_upper(v) for v in values if not v.zero),key=lambda v:ep(v.scale.evaluate())[1],default=scalar(0))
    q=prior.ScaledEnclosure(q_upper.scale,(0,1),q_upper.ledger) if not q_upper.zero else q_upper
    q_Z_upper=max((absolute_upper(v) for v in derivatives if not v.zero),key=lambda v:ep(v.scale.evaluate())[1],default=scalar(0))
    q_Z=symmetric_bound(q_Z_upper)
    return dict(q=q,q_Z=q_Z,a=a,t0=t0,t0_Z=t0_Z,
        record=dict(original_active_constraint='kappa=a+b^2/a<2+eta; eta<=1/2 =>0<a<=3,|b|<=3/2,a*q<=sqrt(a*(3-a)/2)',
            constrained_active_a_C0=a.record(),constrained_active_b_C0=b.record(),
            active_constraints_not_applied_to_original_field_outside_support=True,
            separately_bounded_original_active_cutoff_branches=branches,
            flat_Delta_ge_eta_q_and_q_Z_exact_zero=True,
            entire_original_branch_union_q_C0=q.record(),entire_original_branch_union_q_Z=q_Z.record(),
            ordinary_first_Z_only=True,source_derivative_rows_not_clipped_or_zeroed=True,
            numerical_upper_coordinates_are_bounds_not_source_function_values=True))


def original_periodic_theorem():
    a,t0,q,h,W1,W2,s,psi=sy.symbols('a t0 q h W1 W2 s psi',real=True)
    nu=1+t0*t0+2*q*q;T2=t0*t0*psi+4*t0*q/h*W1+4*q*q*s*W2
    F=(psi+T2)/(2*sy.pi*nu);J=t0*q/h*W1+q*q*(s*W2-psi/2)
    if sy.cancel(a*J/(sy.pi*nu)-a*(F-psi/(2*sy.pi))/2)!=0:raise ArithmeticError('Original primitive/phase identity failed')
    r,cosine=sy.symbols('r cosine',real=True);D=1-2*r*cosine+r*r;w=(cosine-r)/D
    if sy.cancel(sy.diff(w,r)-(-1/D+2*(cosine-r)**2/D**2))!=0:raise ArithmeticError('Original signed-r derivative identity failed')
    return dict(passed=True,original_A_equals_a_half_phase_minus_uniform_angle_identity=True,
        original_smooth_signed_r_w_derivative_identity=True,
        W1_bound_derivation='For |r|<=1/2, original Fourier series gives |W1|<=sum |r|^(k-1)/k<=2. For |r|>=1/2, W1=(angle-psi)/(2r), both angles in[0,2pi], giving |W1|<=2pi. Use safe4pi.',
        sW2_bound_derivation='For |r|<=1/2, the original Fourier series gives sW2<=pi+|r|/(1-|r|)+1/[2(1-|r|)]<=pi+2. For |r|>=1/2, original Mobius formula sW2=[(2-3s)*angle+s*psi+2r*sin(angle)]/(4r^2)<=6pi+2. Use safe100pi.',
        original_active_A_bound='|t0*q|/nu<=1/2, q^2/nu<=1/2, h^-1<=1 imply |A|<=52.5*a<53*a<=159',
        original_active_B_over_E_bound='t0*a=-b, |b|<=3/2, a*q<=sqrt(a*(3-a)/2)<=sqrt(9/8) imply |B/E|<=53*3/2+2sqrt(9/8)<100',
        original_inverse_derivative='F_psi=(1+t^2)/(2pi*nu)>0; fixed fractional phi in[0,1] gives |psi_Z|<=2pi*|nu_Z|+|T2_Z|_fixed_psi',
        cutoff_smooth_seams='sigma and all positive derivatives are flat at0 and1; the active/flat union defines the same original smooth q/A/B functions',
        phase_is_one_period_not_unwrapped_N_y=True)


def whole_period_C1(roots,eta_log,log_a_lower,dstar_log):
    cutoff=whole_cutoff_C1(roots,eta_log,log_a_lower);c=roots['a'][ZERO].ctx;scalar=roots['a'][ZERO].scalar
    if cutoff['a'] is None:
        zero=scalar(0);return dict(values=dict(A=zero,A_Z=zero,B_over_Pstar=zero,B_Z_over_Pstar=zero),record=cutoff['record'])
    a=cutoff['a'];q=cutoff['q'];q_Z=cutoff['q_Z'];t0=cutoff['t0'];t0_Z=cutoff['t0_Z']
    dstar=prior.ScaledEnclosure(prior.FormalScale(a.scale.bases,offset=c.mpf(dstar_log)),1,a.ledger)
    u=(roots['p2'][ZERO]*q).positive_divide(dstar,dstar_log)
    u_Z=(roots['p2'][DZ]*q+roots['p2'][ZERO]*q_Z).positive_divide(dstar,dstar_log)
    U,UZ,T,TZ,Q,QZ,AU,AZ=[absolute_upper(v) for v in (u,u_Z,t0,t0_Z,q,q_Z,a,roots['a'][DZ])]
    one=scalar(1);pi=c.pi;oneplus=one+U
    invD=current.square(current.square(oneplus))*4
    wcap=invD*2;wr=invD+current.square(invD)*8
    W1cap=scalar(4*pi);W2cap=current.square(wcap)*(2*pi);sW2cap=scalar(100*pi);scap=scalar(101*pi)
    # r_Z=h^-3*u_Z; s_Z=-2u*h^-4*u_Z; (h^-1)_Z=-u*h^-3*u_Z.
    # h>=1 and rho>=1/[2(1+|u|)^2] give safe formal bounds without
    # subtracting r from1 or choosing a signed Mobius chart.
    hZ=U*UZ;sZ=hZ*2;W1Z=wr*UZ*(2*pi);W2Z=wcap*wr*UZ*(4*pi)
    firstZ=(TZ*Q+T*QZ+T*Q*hZ)*W1cap+T*Q*W1Z
    weightedW2Z=sZ*W2cap+W2Z
    JZ=firstZ+Q*QZ*scap*2+current.square(Q)*weightedW2Z
    nuZ=T*TZ*2+Q*QZ*4
    T2Z=T*TZ*(4*pi)+firstZ*4+Q*QZ*sW2cap*8+current.square(Q)*weightedW2Z*4
    psiZ=nuZ*(2*pi)+T2Z
    A=scalar(c.mpf((-159,159)))
    fixedAZ=(AZ.positive_divide(a,log_a_lower)+nuZ)*159+AU*JZ*(1/pi)
    Apsi=AU*(T*Q*wcap+current.square(Q)*(current.square(wcap)+c.mpf('.5')))*(1/pi)
    AZtotal=absolute_upper(fixedAZ+Apsi*psiZ)
    E=roots['E'][ZERO];EZ=absolute_upper(roots['E'][DZ]);EU=absolute_upper(E)
    B=E*c.mpf((-100,100))
    BZtotal=EZ*100+EU*(TZ*159+T*AZtotal+(AZ*Q+AU*QZ+AU*Q*hZ)*2
        +AU*Q*(W1Z+wcap*psiZ)*(1/(2*pi)))
    values=dict(A=A,A_Z=symmetric_bound(AZtotal),B_over_Pstar=B,B_Z_over_Pstar=symmetric_bound(BZtotal))
    record=dict(original_three_branch_cutoff_C1_cover=cutoff['record'],original_whole_period_theorem=original_periodic_theorem(),
        original_whole_period_primitive_C0_Z_covers={k:v.record() for k,v in values.items()},
        implicit_actual_inverse_Z_derivative_upper=absolute_upper(psiZ).record(),
        formal_original_Poisson_denominator_inverse_upper=absolute_upper(invD).record(),
        original_entire_period_fractional_phase_cover=c.mpf((0,1)),
        actual_periodic_phase_covered_without_sample_or_selected_inverse=True,
        active_constraint_bounds_do_not_redefine_original_velocity_or_loop=True,
        derivative_covers_can_be_very_wide_and_do_not_prove_small_global_error=True,
        q_and_phase_crossings_have_genuine_C1_covers=True,first_Z_only=True)
    return dict(values=values,record=record)


class NativeFirstBridgeC1Histories:
    def __init__(self,owner):
        if type(owner) is not transfer.NativeTrueChartC1Transfer:raise ValueError('Same true-chart C1 transfer owner required')
        receipt=json.loads((HERE/transfer.RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt[transfer.GATE] or receipt['source_family']!=owner.family:raise ValueError('Accepted true-chart C1 stage required')
        self.owner=owner;self.ctx=owner.ctx;self.family=owner.family;self.service=owner.service;self.coordinates=owner.coordinates
        self.q_owner=owner.owner.owner.owner.owner
        self.service.bind_hashes(receipt['input_hashes'])
        self.service.bind_hashes({name:sha(name) for name in (transfer.RECEIPT,Path(__file__).name,
            PREFIX+'current_generic_shear_loop.py',PREFIX+'flat_pulse_derivatives.py')})
    @native.inlet.source_precision
    def first_bridge(self,Z=(-1,1),N=1024):
        N=density.spatial.density.candidate_integer(N)
        if N<160:raise ValueError('This analytic whole-period primitive cover requires candidate N>=160')
        initial=self.owner.initial_collar(Z,N=N);c=self.ctx;sc=self.owner.owner.binder.sc
        left=sc*c.mpf(3)/4;coordinate=c.mpf((ep(left)[0],1));d=1-left
        geometry=self.owner.geometry.build('bridge_first',0,{'bridge':d},coordinate,
            [{'selected_sc_multiple':'3/4'},'1'])
        source=self.q_owner.query('bridge_first',Z,coordinate);roots=source['source']['roots']
        positive=self.q_owner.owner.owner.decode(self.q_owner.owner.owner.inventory['bridge_first']['actual_positive_denominator_theorem'])
        eta=packets.interval(c,self.q_owner.owner.owner.scales['selected_positive_eta_log'])
        dstar=packets.interval(c,self.q_owner.owner.owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
        primitives=whole_period_C1(roots,eta,positive['log_actual_a_positive_lower'],dstar)
        E=roots['E'][ZERO];E_Z=roots['E'][DZ];packet=source['source']['packet']
        signed=self.owner.owner.signed_owner
        def axial(k):
            row=prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
            return signed.leaf(prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),E.scale.bases,E.ledger)
        got=density.density_Z_kernels(E,E_Z,axial(0),axial(1),primitives['values'],N)
        factors={key:transfer.true_width_kernel(self.coordinates,geometry,rate) for key,rate in RATES.items()}
        contributions={key:self.coordinates.rebase(got['kernels'][key],self.family)*factors[key]['mass'] for key in RATES}
        jets={key:self.coordinates.rebase(got['Z_derivatives'][key],self.family)*factors[key]['mass'] for key in RATES}
        operator=history.C1DuhamelOperator(self.coordinates)
        transfer.append_true_cell(operator,geometry,contributions,jets,self.family)
        correction=operator.apply(initial['correction']['values'],initial['correction']['Z_derivatives'],self.family)
        background_packet=self.owner.owner.native.query('bridge_first',Z,1)
        background=history.packet_history_functions(background_packet,self.coordinates,signed)
        own={key:background['originals'][key]+correction['values'][key] for key in RATES}
        own_Z={key:background['Z_derivatives'][key]+correction['Z_derivatives'][key] for key in RATES}
        record=dict(source_family=self.family,Z_box=c.mpf(Z),candidate_N=N,
            known_original_initial_collar_history=initial['record'],actual_entire_remaining_first_bridge_geometry=geometry['record'],
            original_whole_cell_source_and_cutoff_C1=source['record'],original_whole_period_C1_cover=primitives['record'],
            original_whole_bridge_signed_density_C0_covers={k:v.record() for k,v in got['kernels'].items()},
            original_whole_bridge_signed_density_Z_covers={k:v.record() for k,v in got['Z_derivatives'].items()},
            true_log_radius_signed_C0_contributions={k:v.record() for k,v in contributions.items()},
            true_log_radius_signed_Z_contributions={k:v.record() for k,v in jets.items()},
            actual_first_bridge_exit_correction_C0={k:v.record() for k,v in correction['values'].items()},
            actual_first_bridge_exit_correction_Z={k:v.record() for k,v in correction['Z_derivatives'].items()},
            original_phase1_background_and_separate_P0_Z=background['record'],
            actual_first_bridge_exit_own_history_C0={k:v.record() for k,v in own.items()},
            actual_first_bridge_exit_own_history_Z={k:v.record() for k,v in own_Z.items()},
            no_interval_skipped_between_true_inlet_and_first_bridge_exit=True,
            actual_active_loop_retained_across_cutoff_crossings=True,
            actual_whole_first_bridge_C1_incoming_correction_functions_installed=True,
            first_bridge_covers_are_conservative_not_small_error_or_terminal_closure=True,
            full_active_phase_inverse_point_evaluator_installed=False,
            global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return dict(record=record,correction=correction,background=background,own=own,own_Z=own_Z,
            geometry=geometry,contributions=contributions,Z_derivatives=jets,primitives=primitives,source=source)


@native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic()
    if bridge is None:bridge,_=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime():
        c1=density.NativeDensityC1LocalIntegrals(first.NativePhaseFirstJets(slow.NativeQSlowJets(current.NativeCorrelatedShearQ(prior.NativeSignedInputEnclosures(native.NativeGenericSourcePackets(bridge))))))
        owner=NativeFirstBridgeC1Histories(transfer.NativeTrueChartC1Transfer(history.NativeC1HistoryTransfer(c1)))
        records={}
        for name,Z in (('whole_Z',(-1,1)),('Z_interval',('.49','.51'))):
            records[name]=owner.first_bridge(Z,N=1024)['record'];print('Actual whole first-bridge C1 exit histories:',name,flush=True)
    result=dict(source_family=owner.family,**{GATE:True},candidate_N=1024,
        actual_whole_first_bridge_C1_history_records=records,native_Z_query_count=2,
        actual_original_inlet_to_phase1_C1_histories_installed=True,
        active_bridge_cutoff_or_phase_crossings_not_deleted=True,
        covers_can_be_very_wide_tight_signed_oscillatory_error_still_required=True,
        global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=dict(owner.service.hashes),
        scope='Actual original inlet-to-phase1 first-bridge correction and own C1 history covers on whole Z[-1,1] and[.49,.51], using separately bounded original cutoff branches and analytic whole-period inverse/primitive derivative covers. Conservative bounds; no downstream bridge/global Rc/repair/common N/cone/recursion/full NS admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
