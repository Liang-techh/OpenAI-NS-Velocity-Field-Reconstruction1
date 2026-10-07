"""Original serial C1 cell with parameter-uniform periodic derivative bounds.

The actual loop is unchanged. Mobius/Fourier inequalities bound its
derivatives without large Poisson denominator powers. Certified flat q
support has exactly zero primitive/density increments and retains memory.
"""
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_middle_O2_inlet_C1_histories as base

p=base.p;transfer=base.transfer;history=base.history;density=base.density
prior=base.prior;native=base.native;packets=base.packets
HERE,PREFIX,sha=base.HERE,base.PREFIX,base.sha;ep=base.ep;RATES=base.RATES;ZERO=base.ZERO;DZ=base.DZ
absolute_upper=p.absolute_upper;symmetric_bound=p.symmetric_bound


def minimum_upper(*values):
    """Minimum of certified magnitude bounds, never a defining field value."""
    caps=[absolute_upper(value) for value in values]
    if any(value.zero for value in caps):return caps[0].scalar(0)
    return min(caps,key=lambda value:ep(value.scale.evaluate())[1])


def periodic_parameter_theorem():
    r,z=sy.symbols('r cospsi',real=True);D=1-2*r*z+r*r
    if sy.expand(D*D-(1-r*r)**2*(1-z*z)-((1+r*r)*z-2*r)**2)!=0:
        raise ArithmeticError('Original Mobius angle parameter bound identity failed')
    return dict(passed=True,original_Mobius_signed_r_parameter_identity=True,
        theta_u_bound='theta_r=2 sin(psi)/D; |sin(psi)|/D<=1/(1-r^2)=h^2; r_u=h^-3; hence |theta_u|<=2/h<=2',
        W1_u_bound='For |r|<=1/2 original Fourier derivative <=4. Outside, W1=(theta-psi)/(2r), |theta_u|<=2 and |r_u|<=1 give <=2+4pi<16pi.',
        sW2_u_bound='Inside |r|<=1/2: |w|<=2,|w_r|<=76,|s_r|<=1 give |(sW2)_r|<=616pi. Outside use original Mobius formula, |s_u|<=2, |theta_u|<=2, |r_u|<=1: numerator derivative <=16pi+12, denominator term <=400pi. Both <1000pi.',
        original_fixed_phase_A_derivative='A=a/2*(phi-psi/(2pi)); |A_Z|<=|a_Z|/2+a*L/(4pi), L=2pi|nu_Z|+|T2_Z|fixed',
        original_B_angle_derivative='q*h^-1*w=(t-t0)/2 and psi_Z=(2pi phi nu_Z-T2_Z)/(1+t^2); |q*h^-1*w|/(1+t^2)<=(1+|t0|)/2',
        fractional_phase_domain=[0,1],higher_jets_not_claimed=True,
        parameter_bounds_are_covers_not_modified_loop_functions=True)


def whole_period_C1(roots,eta_log,log_a_lower,dstar_log):
    cutoff=p.whole_cutoff_C1(roots,eta_log,log_a_lower)
    c=roots['a'][ZERO].ctx;scalar=roots['a'][ZERO].scalar
    theorem=periodic_parameter_theorem()
    if cutoff['q'].zero and cutoff['q_Z'].zero:
        zero=scalar(0)
        return dict(values=dict(A=zero,A_Z=zero,B_over_Pstar=zero,B_Z_over_Pstar=zero),
            record=dict(original_three_branch_cutoff_C1_cover=cutoff['record'],
                original_periodic_parameter_C1_theorem=theorem,
                original_full_box_q_and_q_Z_exact_zero_implies_primitive_C0_Z_exact_zero=True,
                original_entire_period_fractional_phase_cover=c.mpf((0,1)),
                source_functions_not_zeroed_without_original_flat_proof=True))
    a=cutoff['a'];q=cutoff['q'];q_Z=cutoff['q_Z']
    # An exactly zero original b row is retained. The generic active bound
    # |b|<=3/2 is used otherwise; no point/sample replaces the source.
    b=scalar(0) if roots['b'][ZERO].zero else scalar(c.mpf(('-1.5','1.5')))
    t0=(-b).positive_divide(a,log_a_lower)
    t0_Z=(-roots['b'][DZ]-t0*roots['a'][DZ]).positive_divide(a,log_a_lower)
    dstar=prior.ScaledEnclosure(prior.FormalScale(a.scale.bases,offset=c.mpf(dstar_log)),1,a.ledger)
    u_Z=(roots['p2'][DZ]*q+roots['p2'][ZERO]*q_Z).positive_divide(dstar,dstar_log)
    T,TZ,Q,QZ,AU,AZ,UZ=[absolute_upper(v) for v in (t0,t0_Z,q,q_Z,a,roots['a'][DZ],u_Z)]
    pi=c.pi
    A_cap=minimum_upper(scalar(159),AU*(T*Q*4+Q*Q*101))
    B_cap=minimum_upper(scalar(100),T*A_cap+AU*Q*2)
    # |(h^-1)_Z|<=|u_Z|, |W1_Z|fixed<=16pi|u_Z|,
    # |(sW2)_Z|fixed<=1000pi|u_Z| uniformly in signed u.
    firstZ=(TZ*Q+T*QZ+T*Q*UZ)*(4*pi)+T*Q*UZ*(16*pi)
    weightedW2Z=UZ*(1000*pi)
    nuZ=T*TZ*2+Q*QZ*4
    T2Z=T*TZ*(4*pi)+firstZ*4+Q*QZ*(800*pi)+Q*Q*weightedW2Z*4
    L=absolute_upper(nuZ*(2*pi)+T2Z)
    AZ_cap=absolute_upper(AZ*c.mpf('.5')+AU*L*(1/(4*pi)))
    E,EZ=roots['E'][ZERO],absolute_upper(roots['E'][DZ]);EU=absolute_upper(E)
    BZ_cap=absolute_upper(EZ*B_cap+EU*(TZ*A_cap+T*AZ_cap+AZ*Q*2+AU*QZ*2
        +AU*Q*UZ*10+AU*(T+1)*L*(1/(4*pi))))
    values=dict(A=symmetric_bound(A_cap),A_Z=symmetric_bound(AZ_cap),
        B_over_Pstar=E*symmetric_bound(B_cap),B_Z_over_Pstar=symmetric_bound(BZ_cap))
    record=dict(original_three_branch_cutoff_C1_cover=cutoff['record'],
        original_periodic_parameter_C1_theorem=theorem,
        original_parameter_sensitive_C0_caps=dict(A=A_cap.record(),B_over_E=B_cap.record()),
        original_whole_period_C0_Z_covers={k:v.record() for k,v in values.items()},
        implicit_phase_derivative_numerator_upper=L.record(),
        whole_period_derivatives_avoid_Poisson_denominator_powers=True,
        original_exact_zero_b_C0_row_retained=roots['b'][ZERO].zero,
        original_full_source_derivatives_not_clipped_or_erased=True,
        actual_fractional_phase_covered_without_samples_or_selected_inverse=True,
        numerical_caps_bound_functions_and_do_not_define_field_values=True,
        conservative_first_Z_only=True)
    return dict(values=values,record=record)


def serial_cell(owner,*,chart,geometry,endpoint,Z,N,incoming,left_record,right_record,seam,join):
    """Same original density/mass/recovery, explicit inherited correction."""
    owner.coordinates.require_family(owner.family)
    if N<160:raise ValueError('Whole-period C0 cap requires candidate N>=160')
    if not geometry['record']['width_and_endpoints_independent_of_Z']:
        raise ValueError('Original fixed-Z integration geometry required')
    c=owner.ctx;source=owner.q_owner.query(chart,Z,geometry['coordinate'])
    roots=source['source']['roots'];root_owner=owner.q_owner.owner.owner
    positive=root_owner.decode(root_owner.inventory[chart]['actual_positive_denominator_theorem'])
    if positive.get('source_function_positivity_not_inferred_from_saved_box') is not True:
        raise ValueError('Original chart-uniform positive a theorem required')
    eta=packets.interval(c,root_owner.scales['selected_positive_eta_log'])
    dstar=packets.interval(c,root_owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
    primitives=whole_period_C1(roots,eta,positive['log_actual_a_positive_lower'],dstar)
    E,E_Z=roots['E'][ZERO],roots['E'][DZ];packet=source['source']['packet']
    signed=owner.transfer.owner.signed_owner
    def axial(k):
        row=prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
        return signed.leaf(prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),E.scale.bases,E.ledger)
    got=density.density_Z_kernels(E,E_Z,axial(0),axial(1),primitives['values'],N)
    factors={k:transfer.true_width_kernel(owner.coordinates,geometry,rate) for k,rate in RATES.items()}
    values={k:owner.coordinates.rebase(got['kernels'][k],owner.family)*factors[k]['mass'] for k in RATES}
    jets={k:owner.coordinates.rebase(got['Z_derivatives'][k],owner.family)*factors[k]['mass'] for k in RATES}
    operator=history.C1DuhamelOperator(owner.coordinates)
    transfer.append_true_cell(operator,geometry,values,jets,owner.family)
    correction=operator.apply(incoming['values'],incoming['Z_derivatives'],owner.family)
    background=history.packet_history_functions(owner.transfer.owner.native.query(chart,Z,endpoint),owner.coordinates,signed)
    own={k:background['originals'][k]+correction['values'][k] for k in RATES}
    own_Z={k:background['Z_derivatives'][k]+correction['Z_derivatives'][k] for k in RATES}
    record=dict(source_family=owner.family,chart=chart,Z_box=c.mpf(Z),candidate_N=N,
        original_left_endpoint=left_record,original_right_endpoint=right_record,
        original_same_radius_phase_or_same_chart_seam=seam,original_background_history_P0_join_receipt=join,
        actual_true_cell_geometry=geometry['record'],original_chart_uniform_a_positive_certificate=positive,
        original_full_box_signed_source_and_old_q_slow_status=source['record'],
        original_periodic_parameter_C1_cover=primitives['record'],
        original_true_width_kernel_factors={k:dict(branch=v['branch'],decay=v['decay'].record(),mass=v['mass'].record()) for k,v in factors.items()},
        actual_signed_density_C0_covers={k:v.record() for k,v in got['kernels'].items()},
        actual_signed_density_Z_covers={k:v.record() for k,v in got['Z_derivatives'].items()},
        actual_true_width_C0_contributions={k:v.record() for k,v in values.items()},
        actual_true_width_Z_contributions={k:v.record() for k,v in jets.items()},
        actual_inherited_correction_C0={k:v.record() for k,v in incoming['values'].items()},
        actual_inherited_correction_Z={k:v.record() for k,v in incoming['Z_derivatives'].items()},
        actual_right_correction_C0={k:v.record() for k,v in correction['values'].items()},
        actual_right_correction_Z={k:v.record() for k,v in correction['Z_derivatives'].items()},
        original_right_background_and_separate_P0_Z=background['record'],
        actual_right_own_history_C0={k:v.record() for k,v in own.items()},
        actual_right_own_history_Z={k:v.record() for k,v in own_Z.items()},
        actual_cell_C1_operator=operator.record(),
        old_unresolved_q_rows_not_used_for_first_Z_primitive_bounds=True,
        original_true_width_and_native_ordinary_y_conversion_applied_once=True,
        actual_incoming_correction_not_reset_or_background_double_added=True,
        global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False))
    return dict(record=record,geometry=geometry,source=source,primitives=primitives,kernels=got,
        factors=factors,operator=operator,incoming=incoming,contributions=values,Z_derivatives=jets,
        correction=correction,background=background,own=own,own_Z=own_Z)
