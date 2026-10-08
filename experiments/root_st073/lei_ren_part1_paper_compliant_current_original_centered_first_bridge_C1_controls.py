"""Centered original first-bridge Z derivatives and actual24 control ranges.

The actual source functions are unchanged. Branchwise a*nu=v and paired
weighted derivatives retain correlations before taking magnitude covers.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_full_Z_paired_C1_controls as accepted
import lei_ren_part1_paper_compliant_current_native_centered_phase_conditioning as centered

base=accepted.base;prior=accepted.prior;first=accepted.first;serial=accepted.serial
HERE,PREFIX,sha=accepted.HERE,accepted.PREFIX,accepted.sha;ep=accepted.ep;encode=accepted.encode
ZERO,DZ,KEYS,N=accepted.ZERO,accepted.DZ,accepted.KEYS,1024
NAME=PREFIX+'current_original_centered_first_bridge_C1_controls.json.gz'
RECEIPT=PREFIX+'current_original_centered_first_bridge_C1_controls_check.json'
GATE='actual_original_full_Z_centered_first_bridge_derivatives_and_24_cell_finite_control_ranges_installed'


def centered_branch_C1(roots,qrows,rrows,branch,eta_log,logamin,dstar_log):
    """Fixed-free-angle weighted derivatives, then the inverse correction once."""
    a=roots['a'][ZERO];c=a.ctx;scalar=a.scalar;absup=first.absolute_upper;minimum=serial.minimum_upper
    if branch['name']=='flat':
        if any(not v.zero for v in (*qrows.values(),*rrows.values())):raise ValueError('Exact flat q/q2 first rows required')
        zero=scalar(0)
        return dict(values=dict(A_Z=zero,B_Z_over_Pstar=zero),record=dict(original_flat_branch_exact_primitive_Z_zero=True))
    if branch['name'] not in ('negative','transition'):raise ValueError('Original cutoff branch required')
    eta=prior.ScaledEnclosure(prior.FormalScale(a.scale.bases,offset=eta_log),1,a.ledger)
    dstar=prior.ScaledEnclosure(prior.FormalScale(a.scale.bases,offset=dstar_log),1,a.ledger)
    DeltaZ=roots['kappa_minus2'][DZ];gamma_lower=eta_log+(c.ln(2) if branch['name']=='negative' else 0)
    gamma=(eta*2-branch['Delta']).positive_intersection(gamma_lower)
    logg=(gamma_lower-c.ln(2))/2
    g=first.current.nonnegative_sqrt(gamma*c.mpf('.5')).positive_intersection(logg)
    gZ=(-DeltaZ).positive_divide(g*4,logg+c.ln(4))
    if branch['name']=='negative':
        s,sp=c.mpf(1),c.mpf(0);sigmaZ=scalar(0);alpha,alphaZ=g,gZ
        M2Z_direct=-DeltaZ*c.mpf('.5');vZ_direct=scalar(0);Hprime=c.mpf(-1)
    else:
        theta=branch['theta'];s,sp=prior.sigma_jets(c,1-theta)[:2]
        sigmaZ=(-DeltaZ).positive_divide(eta,eta_log)*sp
        alpha,alphaZ=g*s,gZ*s+g*sigmaZ
        Hprime=-2*s*sp*(2-theta)-s*s
        M2Z_direct=DeltaZ*(Hprime*c.mpf('.5'));vZ_direct=DeltaZ*(1+Hprime)
    Q,QZ,R,RZ=[absup(v) for v in (qrows[ZERO],qrows[DZ],rrows[ZERO],rrows[DZ])]
    AU,AZ,BU,BZ,EU,EZ=[absup(v) for v in (a,roots['a'][DZ],roots['b'][ZERO],roots['b'][DZ],roots['E'][ZERO],roots['E'][DZ])]
    alpha0,alpha1=absup(alpha),absup(alphaZ)
    # sqrt(a)*q_Z=alpha_Z-alpha*a_Z/(2a), before b or a weighting.
    sqrt_a_qZ=alpha1+(alpha0*AZ).positive_divide(a*2,logamin+c.ln(2))
    limit=2 if branch['name']=='negative' else c.mpf('2.5')
    sqrt_a_cap=minimum(absup(first.current.nonnegative_sqrt(a)),scalar(c.sqrt(limit)))
    rho=scalar(c.sqrt(limit)) # |b|/sqrt(a)<=sqrt(kappa)<=sqrt(limit).
    b_qZ=minimum(BU*QZ,rho*sqrt_a_qZ)
    a_qZ=minimum(AU*QZ,sqrt_a_cap*sqrt_a_qZ)
    b_R=minimum(BU*R,rho*alpha0*Q)
    raw_M2=a*rrows[ZERO];M2=minimum(raw_M2,alpha0*alpha0,scalar(c.mpf('1.5')))
    raw_M2Z=roots['a'][DZ]*rrows[ZERO]+a*rrows[DZ]
    M2Z=minimum(raw_M2Z,M2Z_direct)
    raw_vZ=DeltaZ+raw_M2Z*2;vZ=minimum(raw_vZ,vZ_direct)
    cZ=absup(roots['p2'][DZ].positive_divide(dstar,dstar_log));pi=c.pi
    # P+uP_u<=1, |P_u|<=20pi, H+uH_u/2<=pi+3/2.
    # The accepted H<=pi also gives |uH_u|/2<=2pi+3/2.
    bqP_Z=BZ*Q*(4*pi)+b_qZ+b_R*cZ*(20*pi)
    aqP_Z=AZ*Q*(4*pi)+a_qZ+M2*cZ*(20*pi)
    M2H_Z=M2Z*(pi+c.mpf('1.5'))+AZ*R*(2*pi+c.mpf('1.5'))+M2*Q*cZ*(1000*pi)
    Jtilde_Z=bqP_Z+M2H_Z+M2Z*pi
    # |K|<=pi/2 from original Phi=(psi+4K)/(2pi); v>=2.
    KZ=(Jtilde_Z+vZ*(pi/2))*c.mpf('.5')
    inverse=KZ*(3/pi)
    AZcap=AZ*c.mpf('.5')+inverse
    # Exact M inverse term is 2v*t*K_Z/[pi*(1+t²)], magnitude<=3|K_Z|/pi.
    MZcap=BZ+aqP_Z*(1/pi)+inverse
    BZcap=EZ*c.mpf('1.5')+EU*MZcap*c.mpf('.5')
    values=dict(A_Z=first.symmetric_bound(AZcap),B_Z_over_Pstar=first.symmetric_bound(BZcap))
    record=dict(original_branch=branch['name'],fixed_free_angle_derivatives_not_total_inverse_rows=True,
        original_v='a*nu=kappa+2*a*q²',v_positive_lower='2',v_upper='3',
        original_branch_v_Z_exact='0' if branch['name']=='negative' else '(1+Hcut_prime(theta))*Delta_Z',
        Hcut='sigma(1-theta)^2*(2-theta)',Hcut_prime_wrt_theta=Hprime,
        original_eta_Z_exact_zero=True,original_gamma_C0=gamma.record(),original_unmodulated_g_C0_Z=[g.record(),gZ.record()],
        original_sigma_C0=s,original_sigma_Z=sigmaZ.record(),original_alpha_C0_Z=[alpha.record(),alphaZ.record()],
        original_sqrt_a_times_q_Z_upper=sqrt_a_qZ.record(),original_b_times_q_Z_upper=b_qZ.record(),original_a_times_q_Z_upper=a_qZ.record(),
        original_b_times_q_squared_upper=b_R.record(),original_aq2_raw=raw_M2.record(),original_aq2_upper=M2.record(),
        original_aq2_Z_raw=raw_M2Z.record(),original_aq2_Z_direct=M2Z_direct.record(),original_aq2_Z_upper=M2Z.record(),
        original_v_Z_raw=raw_vZ.record(),original_v_Z_direct=vZ_direct.record(),original_v_Z_upper=vZ.record(),
        original_chi_Z_upper=cZ.record(),original_fixed_bqP_Z_upper=bqP_Z.record(),original_fixed_aqP_Z_upper=aqP_Z.record(),
        original_fixed_M2H_Z_upper=M2H_Z.record(),original_fixed_centered_Jtilde_Z_upper=Jtilde_Z.record(),
        original_fixed_centered_K_Z_upper=KZ.record(),original_inverse_contribution_upper=inverse.record(),
        original_total_A_Z_upper=AZcap.record(),original_total_M_Z_upper=MZcap.record(),original_total_B_Z_upper=BZcap.record(),
        original_total_primitive_Z_covers=base.records(values),
        exact_weighted_H_rule='(a*q²*H)_Z=M2_Z*(H+u*H_u/2)-a_Z*q²*u*H_u/2+M2*q*chi_Z*H_u',
        original_free_M_Z='-b_Z*(phi-psi/(2pi))-(a*qP)_Z/pi',
        original_total_M_inverse='2*v*t*K_Z/(pi*(1+t²)); 2|t|/(1+t²)<=1',
        original_paired_weights_before_magnitude_bounds=True,only_original_a_g_and_dstar_positive_denominators=True,
        no_cutoff_q_floor_or_source_derivative_clipping=True)
    return dict(values=values,record=record)


def whole_centered_C1(roots,eta_log,logamin,dstar_log):
    c=roots['a'][ZERO].ctx
    if ep(eta_log)[1]>ep(-c.ln(2))[0]:raise ValueError('Original eta<=1/2 required')
    candidates,empty=accepted.cutoff.conditional_cutoff_branches(roots['kappa_minus2'][ZERO],eta_log)
    rows=[];covers=[]
    for branch in candidates:
        restricted,correlation=accepted.correlated_active_roots(roots,branch,logamin)
        if restricted is None:empty.append(dict(name=branch['name'],proof=correlation));continue
        q=accepted.cutoff.conditional_q_jet(restricted,eta_log,logamin,branch)['rows']
        square=accepted.q2.conditional_q2_jet(restricted,eta_log,logamin,branch)
        got=centered_branch_C1(restricted,q,square,branch,eta_log,logamin,dstar_log)
        rows.append(dict(name=branch['name'],condition=branch['condition'],conditional_Delta=branch['Delta'].record(),theta=branch.get('theta'),
            original_active_correlation=correlation,original_q_C0_Z=[q[ZERO].record(),q[DZ].record()],
            original_q_squared_C0_Z=[square[ZERO].record(),square[DZ].record()],original_centered_first_Z=got['record']))
        covers.append(got['values'])
    if not covers:raise ValueError('Complete nonempty original cutoff union required')
    hull={key:accepted.density.local.same_source_union([row[key] for row in covers]) for key in ('A_Z','B_Z_over_Pstar')}
    return dict(values=hull,record=dict(original_native_log_bases=roots['a'][ZERO].scale.bases,
        original_source_root_C0_Z={name:{str(order):row[order].record() for order in (ZERO,DZ)} for name,row in roots.items()},
        original_eta_log=eta_log,original_log_a_positive_lower=logamin,original_dstar_log=dstar_log,
        complete_original_cutoff_branches=rows,branches_proved_empty=empty,whole_original_centered_primitive_Z_hulls=base.records(hull),
        original_centered_source_theorem=centered.exact_kernel_theorem(),original_paired_Poisson_theorem=accepted.paired.paired_parameter_theorem(),
        all_original_body_transition_flat_branches_hulled_before_primitive_selection=True,first_Z_only=True))


class OriginalCenteredFirstBridgeC1Controls:
    def __init__(self,bridge=None):
        self.previous=accepted.OriginalFullZPairedC1Controls(bridge);self.base=self.previous.base
        self.ctx=self.base.ctx;self.coordinates=self.base.coordinates;self.family=self.base.family;self.target_owner=self.base.target_owner
        self.hashes=dict(self.previous.hashes);self.manifest=base.driver.accepted(self.hashes,accepted,self.family)
        base.driver.accepted(self.hashes,centered,self.family)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name);self.source=None

    @base.native.inlet.source_precision
    def refine(self,*,N=1024):
        if type(N) is not int or N!=1024:raise ValueError('Accepted original N1024 required')
        source=self.base.assemble();c=self.ctx;coords=self.coordinates
        prior_source=self.manifest['actual_full_Z_paired_C1_original24_source_ranges']
        with mp.workdps(c.dps+40):
            # Reuse accepted integrated derivative covers, not numerical source owners.
            for cell,row in zip(source['cells'],prior_source['actual_original24_source_cell_records'],strict=True):
                for key in ('actual_original_geometry','actual_original_source_cell_binding','actual_local_C0_contributions','actual_right_background_C0','actual_right_background_Z','original_separate_P0','original_separate_P0_Z'):
                    if encode(cell['record'][key])!=row[key]:raise ValueError('Original source/C0/geometry/background/P0 changed: '+key)
                jets=self.base.restore_group(row['actual_local_Z_contributions'])
                cell['Z_derivatives']=jets;cell['operator'].Z_increments=jets
            geom=self.target_owner.geometry(*base.parameters.ROUTE[1]);query=self.target_owner.q_owner.query('bridge_first',(-1,1),geom['coordinate'])
            roots=query['source']['roots'];root_owner=self.target_owner.q_owner.owner.owner
            positive=root_owner.decode(root_owner.inventory['bridge_first']['actual_positive_denominator_theorem'])
            if not positive['source_function_positivity_not_inferred_from_saved_box']:raise ValueError('Whole original first bridge positive a theorem required')
            eta=base.packets.interval(c,root_owner.scales['selected_positive_eta_log']);logamin=positive['log_actual_a_positive_lower']
            dstar=base.packets.interval(c,root_owner.scales['logarithmic_selected_positive_lower_constants']['d_star'])
            got=whole_centered_C1(roots,eta,logamin,dstar);oldproof=prior_source['actual_full_Z_paired_C1_refinement_records'][0]
            bases=roots['a'][ZERO].scale.bases;ledger=roots['a'][ZERO].ledger
            def restore(record):
                s=record['formal_positive_scale']
                return prior.ScaledEnclosure(prior.FormalScale(bases,tuple(s['source_exponents'])+(s['radius_power'],),base.packets.interval(c,s['additional_log_interval'])),
                    base.packets.interval(c,record['coefficient_interval']),ledger)
            primitives={key:restore(row) for key,row in oldproof['selected_original_primitive_C0_Z'].items()}
            comparisons={}
            for key in ('A_Z','B_Z_over_Pstar'):primitives[key],comparisons[key]=accepted.smaller_cover(primitives[key],got['values'][key])
            # All E/EZ/V/VZ stay genuine rows of the freshly queried original packet.
            signed=self.target_owner.transfer.owner.signed_owner;packet=query['source']['packet']
            def axial(k):
                row=prior.signed.ordinary_axial_coefficient(packet.velocity['axial'][0],k)
                return signed.leaf(prior.signed.expressions.RadiusPolynomial(packet.algebra,{0:row}),bases,ledger)
            density=accepted.density.density_Z_kernels(roots['E'][ZERO],roots['E'][DZ],axial(0),axial(1),primitives,N)
            cell=source['cells'][1];oldjets=cell['Z_derivatives'];newjets={};decisions={};factors={}
            for key,rate in base.history.RATES.items():
                factors[key]=base.preceding.downstream.transfer.true_width_kernel(coords,geom,rate)
                new=coords.rebase(density['Z_derivatives'][key],self.family)*factors[key]['mass']
                newjets[key],decisions[key]=accepted.smaller_cover(oldjets[key],new)
            cell['Z_derivatives']=newjets;cell['operator'].Z_increments=newjets
            proof=dict(original_cell_index=1,label='active_first_bridge',chart='bridge_first',source_family=self.family,candidate_N=N,Z_box=c.mpf((-1,1)),
                original_geometry=geom['record'],fresh_original_source_query=query['record'],original_positive_a_theorem=positive,
                original_complete_centered_first_Z=got['record'],previous_accepted_primitive_C0_Z=oldproof['selected_original_primitive_C0_Z'],
                selected_original_primitive_C0_Z=base.records(primitives),primitive_Z_comparisons=comparisons,
                original_density_velocity_C0_Z=base.records(density['velocities']),original_signed_density_Z=base.records(density['Z_derivatives']),
                original_mass_and_decay={key:{field:value[field].record() for field in ('mass','decay')} for key,value in factors.items()},
                previous_accepted_local_Z=base.records(oldjets),selected_local_Z=base.records(newjets),local_Z_comparisons=decisions,
                same_original_C0_geometry_background_P0_and_source_functions=True)
            incoming=dict(values={key:coords.scalar(0) for key in KEYS},Z_derivatives={key:coords.scalar(0) for key in KEYS});budgets=[]
            for index,cell in enumerate(source['cells']):
                row=cell['record'];out=cell['operator'].apply(incoming['values'],incoming['Z_derivatives'],self.family);bg=cell['background']
                for key in KEYS:
                    if encode(out['values'][key].record())!=prior_source['actual_original24_source_cell_records'][index]['actual_right_correction_C0'][key]:
                        raise ValueError('Original C0 affine history changed')
                ownZ={key:bg['Z_derivatives'][key]+out['Z_derivatives'][key] for key in KEYS}
                cell['incoming']=incoming;cell['correction']=out
                row.update(actual_cell_C1_operator=cell['operator'].record(),actual_local_Z_contributions=base.records(cell['Z_derivatives']),
                    actual_inherited_correction_Z=base.records(incoming['Z_derivatives']),actual_right_correction_Z=base.records(out['Z_derivatives']),
                    actual_right_own_history_Z=base.records(ownZ),absolute_pressure_Z=(bg['P0_Z']+ownZ['p']).record(),
                    original_C0_and_background_P0_unchanged=True)
                if index==1:row['centered_first_Z_refinement_binding']=dict(module=Path(__file__).name)
                budgets.append(dict(label=row['label'],outgoing_Z_log_caps={key:accepted.magnitude_log(v) for key,v in out['Z_derivatives'].items()}));incoming=out
            a=source['amplitude'];target=dict(**base.fixed.fixed_N_target_rows(incoming['values'],incoming['Z_derivatives'],a['A'],a['AZ'],a['logA'],a['mu_source'],a['logmu']),
                **{key:a[key] for key in ('A','AZ','mu','logA','logmu')})
            source.update(history=incoming['values'],Z_derivatives=incoming['Z_derivatives'],target=target)
            source['record'].update(actual_original24_source_cell_records=[cell['record'] for cell in source['cells']],
                actual_Rc_correction_C0_Z=[base.records(incoming['values']),base.records(incoming['Z_derivatives'])],
                actual_Rc_joint_target_C0_Z=[base.records(target['values']),base.records(target['Z_derivatives'])],
                actual_Rc_joint_numerator_C0_Z=[target['joint_numerator'].record(),target['joint_numerator_Z'].record()],
                actual_centered_first_bridge_C1_refinement=proof,actual_updated_24_cell_Z_error_budget=budgets,
                accepted_other23_local_source_Z_covers_retained=True,original_C0_histories_background_P0_geometry_and_phase_unchanged=True)
        self.source=source;self.base.hashes.update(self.hashes);self.hashes.update(self.target_owner.service.hashes)
        print('Original centered first-bridge local Z reductions',sum(row['strict_absolute_upper_reduction'] for row in decisions.values()),flush=True)
        return source

    def finite_controls(self,source,*,iterations=3):
        if source is not self.source or source is None:raise ValueError('Issued original centered source ranges required')
        result=self.base.finite_controls(source,iterations=iterations);self.hashes.update(self.base.hashes);return result


@base.native.inlet.source_precision
def run():
    began=time.monotonic();bridge,construction=base.native.inlet.native_bridge_owner()
    with base.native.inlet.CheckedSourceRuntime() as runtime:
        owner=OriginalCenteredFirstBridgeC1Controls(bridge);source=owner.refine();result=owner.finite_controls(source)
    old=owner.manifest['actual_refined_finite_control_diagnostics']
    report=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
        actual_original_centered_first_bridge24_source_ranges=source['record'],actual_centered_finite_control_diagnostics=result['report'],
        previous_accepted_report=dict(filename=accepted.NAME,sha256=sha(accepted.NAME)),previous_N_scaled_target_C0_Z=old['actual_N_scaled_target_C0_Z_ranges'],
        previous_fixed_N_contraction_diagnostic=old['fixed_N_sufficient_contraction_diagnostic'],
        original_source_construction=construction,original_source_runtime=runtime.record(),
        full24_original_C1_integral_range_transport_enclosed=True,actual_finite_picard_and_residual_C0_Z_ranges_installed=True,
        actual_five_controls_installed=False,functional_terminal_identity_solved=False,current_whole_N_selected=False,**dict.fromkeys(base.packets.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Original full-Z first bridge centered source derivatives propagate through the actual24 N1024 chain and finite control/residual ranges. Other23 accepted source integrals, C0/geometry/P0/phase and exact functions unchanged. No global N, fixed point, terminal, heat/cone/recursion/pulse/full corrected NS completion.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(encode(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Actual centered first-bridge original24 control diagnostics generated',flush=True);return report


if __name__=='__main__':run()
