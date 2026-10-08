"""Independent centered first-Z algebra and original24 successor acceptance."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_centered_first_bridge_C1_controls as current
import lei_ren_part1_paper_compliant_current_original_full_Z_paired_C1_controls_check as preceding

base=current.base;accepted=current.accepted;checked=preceding.checked;ep=current.ep;iv=base.packets.interval
encoded=preceding.encoded;KEYS=current.KEYS;ZERO,DZ=current.ZERO,current.DZ


def identities():
    z,psi,phi=sy.symbols('Z psi phi',real=True)
    a,b,q,x=[sy.Function(name)(z) for name in ('a','b','q','chi')]
    u=x*q;P=sy.Function('P')(u);H=sy.Function('H')(u);R=q*q;M2=a*R
    Pu=sy.Subs(sy.Derivative(sy.Function('P')(sy.Symbol('s')),sy.Symbol('s')),sy.Symbol('s'),u)
    Hu=sy.Subs(sy.Derivative(sy.Function('H')(sy.Symbol('s')),sy.Symbol('s')),sy.Symbol('s'),u)
    qz,xz,az,bz=[sy.diff(value,z) for value in (q,x,a,b)]
    qP_Z=qz*(P+u*Pu)+R*xz*Pu
    M2H_Z=sy.diff(M2,z)*(H+u*Hu/2)-az*R*u*Hu/2+M2*q*xz*Hu
    assert sy.simplify(sy.diff(q*P,z)-qP_Z)==0
    assert sy.simplify(sy.diff(M2*H,z)-M2H_Z)==0
    J=-b*q*P+M2*(H-psi/2)
    JZ=-bz*q*P-b*qP_Z+M2H_Z-sy.diff(M2,z)*psi/2
    assert sy.simplify(sy.diff(J,z)-JZ)==0
    alpha=sy.sqrt(a)*q
    assert sy.simplify(sy.diff(alpha,z)-alpha*az/(2*a)-sy.sqrt(a)*qz)==0
    eta=sy.Symbol('eta',positive=True);Delta=sy.Function('Delta')(z)
    bodyR=(2*eta-Delta)/(2*a);kappa=2+Delta
    assert sy.simplify(kappa+2*a*bodyR-(2+2*eta))==0
    assert sy.simplify(sy.diff(a*bodyR,z)+sy.diff(Delta,z)/2)==0
    theta=sy.Function('theta')(z);s=sy.Function('sigma')(1-theta)
    sp=sy.Subs(sy.Derivative(sy.Function('sigma')(sy.Symbol('s')),sy.Symbol('s')),sy.Symbol('s'),1-theta)
    Hcut=s*s*(2-theta);Hprime=-2*s*sp*(2-theta)-s*s
    assert sy.simplify(sy.diff(Hcut,z)-Hprime*sy.diff(theta,z))==0
    # Ordinary Delta_Z/eta cancels the eta in M2; eta has no slow derivative.
    transitionM2=eta*Hcut/2
    assert sy.simplify(sy.diff(transitionM2,z)-Hprime*eta*sy.diff(theta,z)/2)==0
    v=sy.Function('v')(z);K=J/v
    assert sy.simplify(sy.diff(K,z)-(JZ-K*sy.diff(v,z))/v)==0
    t,kz,nu,p_psi=sy.symbols('t K_Z nu P_psi',real=True)
    psiz=-4*nu*kz/(1+t*t)
    A_inverse=-a*psiz/(4*sy.pi)
    assert sy.simplify(A_inverse.subs(nu,v/a)-v*kz/(sy.pi*(1+t*t)))==0
    M_psi=b/(2*sy.pi)-a*q*p_psi/sy.pi
    assert sy.simplify(M_psi+a*(-b/a+2*q*p_psi)/(2*sy.pi))==0
    M_inverse=(-a*t/(2*sy.pi))*psiz
    assert sy.simplify(M_inverse.subs(nu,v/a)-2*v*t*kz/(sy.pi*(1+t*t)))==0
    assert sy.expand((1+t*t)**2-4*t*t-(t*t-1)**2)==0
    theorem=current.centered.exact_kernel_theorem()
    bounds=theorem['weighted_partial_angle_kernel_bounds']['uniform_all_finite_signed_u_and_partial_psi']
    assert bounds['H']=='pi'
    assert theorem['centered_original_phase_theorem']['active_v_range']=='2<=v<=2+2eta<=3'
    return dict(passed=True,independent_centered_weighted_and_total_inverse_first_Z_identities=13,
        body_v_Z_exact_zero=True,transition_Hcut_prime_is_theta_derivative=True,
        original_H_full_period_nonnegative_integral_bound='pi; accepted periodic mixed-conditioning proof',
        inverse_phase_contribution_applied_once=True,exact_M_inverse='2*v*t*K_Z/(pi*(1+t²))',
        inverse_M_magnitude_inequality='(1+t²)²-4t²=(t²-1)²>=0',no_q_positive_floor=True)


def branch_and_primitive_checks(owner,source,saved,previous):
    c=owner.ctx;proof=saved['actual_centered_first_bridge_C1_refinement'];old=previous['actual_full_Z_paired_C1_refinement_records'][0]
    body=proof['original_complete_centered_first_Z'];oldbody=old['actual_complete_paired_cutoff_C1']
    assert proof['source_family']==owner.family and proof['candidate_N']==1024 and ep(iv(c,proof['Z_box']))==(-1,1)
    assert body['original_native_log_bases']==oldbody['original_native_log_bases']
    assert body['original_source_root_C0_Z']==oldbody['original_source_root_C0_Z']
    assert encoded(source['record'])==saved,checked.first_difference(encoded(source['record']),saved)
    assert body['all_original_body_transition_flat_branches_hulled_before_primitive_selection']
    assert body['branches_proved_empty']==oldbody['branches_proved_empty']
    rows=body['complete_original_cutoff_branches'];oldrows=oldbody['actual_conditional_branches']
    assert [r['name'] for r in rows]==[r['name'] for r in oldrows]==['negative','transition','flat']
    bases=tuple(iv(c,v) for v in body['original_native_log_bases'])
    ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
        positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    restore=lambda r:preceding.native_restore(r,bases,ledger)
    branches=0;comparisons=0;vchecks=0
    for row,oldrow in zip(rows,oldrows,strict=True):
        assert row['condition']==oldrow['condition']
        assert row['conditional_Delta']==oldrow['conditional_Delta_C0']
        assert row['original_q_C0_Z']==[oldrow['q_C0'],oldrow['q_Z']]
        assert row['original_q_squared_C0_Z']==[oldrow['direct_q_squared_C0'],oldrow['direct_q_squared_Z']]
        record=row['original_centered_first_Z'];branches+=1
        if row['name']=='flat':
            assert record['original_flat_branch_exact_primitive_Z_zero'];continue
        assert record['fixed_free_angle_derivatives_not_total_inverse_rows']
        assert record['original_eta_Z_exact_zero'] and record['v_positive_lower']=='2' and record['v_upper']=='3'
        M2Z=restore(record['original_aq2_Z_upper']);vZ=restore(record['original_v_Z_upper'])
        if row['name']=='negative':assert restore(record['original_v_Z_direct']).zero and vZ.zero
        elif row['name']=='transition':
            theta=iv(c,row['theta']);s,sp=current.prior.sigma_jets(c,1-theta)[:2]
            hp=-2*s*sp*(2-theta)-s*s
            assert ep(hp)==ep(iv(c,record['Hcut_prime_wrt_theta']))
        vchecks+=1
        bqP,aqP,M2H,JZ,KZ,inverse,A,M,B=[restore(record[key]) for key in (
            'original_fixed_bqP_Z_upper','original_fixed_aqP_Z_upper','original_fixed_M2H_Z_upper',
            'original_fixed_centered_Jtilde_Z_upper','original_fixed_centered_K_Z_upper','original_inverse_contribution_upper',
            'original_total_A_Z_upper','original_total_M_Z_upper','original_total_B_Z_upper')]
        assert checked.overlaps(JZ,bqP+M2Z*c.pi+M2H)
        assert checked.overlaps(KZ,(vZ*(c.pi/2)+JZ)*c.mpf('.5'))
        assert checked.overlaps(inverse,KZ*(3/c.pi))
        roots=body['original_source_root_C0_Z'];AZ=current.first.absolute_upper(restore(roots['a'][str(DZ)]))
        BZ=current.first.absolute_upper(restore(roots['b'][str(DZ)]));E=current.first.absolute_upper(restore(roots['E'][str(ZERO)]))
        EZ=current.first.absolute_upper(restore(roots['E'][str(DZ)]))
        assert checked.overlaps(A,inverse+AZ*c.mpf('.5'))
        assert checked.overlaps(M,inverse+BZ+aqP*(1/c.pi))
        assert checked.overlaps(B,EZ*c.mpf('1.5')+E*M*c.mpf('.5'))
        for key in ('A_Z','B_Z_over_Pstar'):preceding.symmetric(restore(record['original_total_primitive_Z_covers'][key]))
        comparisons+=6
    for key in ('A','B_over_Pstar'):
        assert proof['selected_original_primitive_C0_Z'][key]==old['selected_original_primitive_C0_Z'][key]
    for key in KEYS:
        assert proof['previous_accepted_local_Z'][key]==previous['actual_original24_source_cell_records'][1]['actual_local_Z_contributions'][key]
    for index,(row,oldrow) in enumerate(zip(saved['actual_original24_source_cell_records'],previous['actual_original24_source_cell_records'],strict=True)):
        if index!=1:assert row['actual_local_Z_contributions']==oldrow['actual_local_Z_contributions']
    return dict(passed=True,same_fresh_original_first_bridge_source_basis_roots_and_q_first_rows=True,
        complete_original_cutoff_branches=branches,independent_centered_v_derivative_and_denominator_checks=vchecks,
        rearranged_centered_and_total_inverse_cap_overlaps=comparisons,
        original_C0_primitives_and_other23_local_Z_covers_unchanged=True,
        primitive_Z_strict_reductions=sum(r['strict_absolute_upper_reduction'] for r in proof['primitive_Z_comparisons'].values()),
        local_Z_strict_reductions=sum(r['strict_absolute_upper_reduction'] for r in proof['local_Z_comparisons'].values()))


@base.native.inlet.source_precision
def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    flags=('actual_five_controls_installed','functional_terminal_identity_solved','current_whole_N_selected',*base.packets.OPEN)
    assert manifest[current.GATE] and manifest['candidate_N']==1024 and all(manifest[key] is False for key in flags)
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    previous=json.loads(gzip.decompress((current.HERE/accepted.NAME).read_bytes()))
    assert manifest['previous_accepted_report']['sha256']==current.sha(accepted.NAME)
    bridge,_=base.native.inlet.native_bridge_owner();checks={'independent_centered_first_Z_identities':identities()}
    with base.native.inlet.CheckedSourceRuntime():
        service=current.OriginalCenteredFirstBridgeC1Controls(bridge);source=service.refine();owner=service.base
        with mp.workdps(owner.ctx.dps+40):
            saved=manifest['actual_original_centered_first_bridge24_source_ranges']
            checks['original_C1_seam_binding']=preceding.derivative_identities(owner.ctx)
            checks['complete_original_centered_first_bridge']=branch_and_primitive_checks(owner,source,saved,previous['actual_full_Z_paired_C1_original24_source_ranges'])
            checks['actual24_affine_and_new_joint_target']=preceding.propagate(owner,source,saved,previous['actual_full_Z_paired_C1_original24_source_ranges'])
            print('Centered first-Z algebra, complete original branches and actual24 propagation PASS',flush=True)
        result=service.finite_controls(source)
        with mp.workdps(owner.ctx.dps+40):
            checks['same_exact_functions_and_actual_finite_controls']=preceding.finite_controls(owner,source,result,
                manifest['actual_centered_finite_control_diagnostics'],previous['actual_refined_finite_control_diagnostics'])
            count=0
            for callback in (lambda:service.refine(N=True),lambda:service.refine(N=2048),lambda:service.finite_controls(None),lambda:service.finite_controls(dict(source))):
                try:callback()
                except (ValueError,TypeError):count+=1
                else:raise AssertionError('Invalid frequency or unissued centered source admitted')
            assert count==4;checks['issued_original_source_and_frequency_guards']=dict(passed=True,rejections=count)
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,candidate_N=1024,**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw)),
        full24_original_C1_integral_range_transport_enclosed=True,actual_finite_picard_and_residual_C0_Z_ranges_installed=True,
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Full-Z original first bridge centered weighted first derivatives, complete source/seam coverage, actual24 affine/target and same exact finite control/residual functions. Ancestor suites reused. No global N, fixed point, terminal, heat/stress/flat/recursion/pulse/corrected NS admission.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.encode(receipt),indent=2).encode()+b'\n')
    print('Actual original centered first-bridge24 controls focused acceptance PASS',flush=True);return receipt


if __name__=='__main__':run()
