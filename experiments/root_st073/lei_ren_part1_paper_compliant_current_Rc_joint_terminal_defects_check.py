"""Independent joint-density/transport/normalization math and saved replay.

Only the new algebra and source composition are checked. Accepted upstream
owners and phase inverse solvers are not reconstructed or rerun.
"""
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_Rc_joint_terminal_defects as source
import lei_ren_part1_paper_compliant_current_transition_normalized_q_source_check as compare
from lei_ren_part1_paper_compliant_current_transition_complete_prefix_check import exact_replay_equal

HERE,sha,ep,iv=source.HERE,source.sha,source.ep,source.iv


def independent_joint_checks(family):
    c=MPIntervalContext();c.dps=100;p=mp.mp.clone();p.dps=350
    coordinates=source.native.HalfPstarCoordinates(c,c.mpf(0),family);t=coordinates.scalar(1)
    checks=0;reflection=0;transport=0;normalization=0
    for mu0 in ('.003','1e-300'):
        mu=p.mpf(mu0);epsilon=source.rc.density.density.factored_expm1(t.scalar(c.mpf(mu0))*c.mpf('-2.5'))
        assert not epsilon.zero and ep(epsilon.coefficient)[1]<0
        for Z0 in ('-.37','.37'):
            Z=p.mpf(Z0)
            for s0 in ('.2','10.5'):
                s=p.mpf(s0)
                for phi0 in ('.17','.63'):
                    sn=p.sin(2*p.pi*p.mpf(phi0))
                    E=lambda y,z:p.exp(p.mpf('.3')-y/2)/(1+z*z)
                    A=lambda y,z:(p.mpf('.3')+p.mpf('.002')*y+p.mpf('.02')*z)*sn
                    B=lambda y,z:(p.mpf('.2')+p.mpf('.003')*y-p.mpf('.01')*z)*sn
                    a=lambda value:t.scalar(c.mpf(str(value)))
                    primitive=dict(A=a(A(s,Z)),A_Z=a(p.mpf('.02')*sn),A_y=a(p.mpf('.002')*sn),
                        B_over_Pstar=a(B(s,Z)),B_Z_over_Pstar=a(-p.mpf('.01')*sn),B_y_over_Pstar=a(p.mpf('.003')*sn))
                    for nonzero_V in (False,True):
                        V=lambda y,z:(p.mpf('.1')+p.mpf('.03')*y)*z if nonzero_V else p.mpf(0)
                        # Independently differentiate physical k and effective-amplitude*m.
                        def direct(y,z):
                            e=E(y,z);dE=e*p.expm1(A(y,z)/source.N);dV=B(y,z)/source.N
                            k=V(y,z)*dE+e*dV+dE*dV;m=dV
                            return k-e*p.exp(-p.mpf('2.5')*mu)*m
                        got=source.joint_density(a(E(s,Z)),a(p.diff(lambda z:E(s,z),Z)),a(V(s,Z)),
                            a(p.diff(lambda z:V(s,z),Z)),primitive,epsilon,Ey=None if nonzero_V else a(-E(s,Z)/2))
                        compare.enclosed(p,c,got['C'],direct(s,Z));checks+=1
                        compare.enclosed(p,c,got['C_Z'],p.diff(lambda z:direct(s,z),Z));checks+=1
                        if not nonzero_V:
                            compare.enclosed(p,c,got['C_y'],p.diff(lambda y:direct(y,Z),s));checks+=1
                            mean=source.buffer.reflected_density_mean(a(E(s,Z)),primitive,source.N)['k']
                            x=A(s,Z)/source.N;z=B(s,Z)/source.N;e=E(s,Z);f=p.exp(-p.mpf('2.5')*mu)
                            pair=(e*z*(p.exp(x)-f)-e*z*(p.exp(-x)-f))/2
                            compare.enclosed(p,c,mean,pair);reflection+=1
        # Exact constant-primitive integrals independently verify different
        # rates and the sign/exponent of the effective amplitude.
        for L0 in ('.02','11','40'):
            L=p.mpf(L0);D=p.mpf(3);E0=p.mpf('1.3');x=p.mpf('.0003');z=p.mpf('.0007');f=p.exp(-p.mpf('2.5')*mu)
            Aend=E0*p.exp(-(L+D)/2)*f
            k=E0*z*p.exp(x)*p.exp(-p.mpf('1.5')*(L+D))*p.expm1(L)
            m=z*(-p.expm1(-L))*p.exp(-D)
            joint=E0*z*(p.exp(x)-f)*p.exp(-p.mpf('1.5')*(L+D))*p.expm1(L)
            assert p.almosteq(k-Aend*m,joint,rel_eps=p.mpf('1e-320'));transport+=1
    # Independent quotient differentiation checks all five rows at both Z signs.
    for Z0 in ('-.37','.37'):
        Z=p.mpf(Z0);mu=p.mpf('.003');A=lambda z:p.mpf('1.3')/(1+z*z)
        functions=dict(m=lambda z:p.mpf('.1')+z*z,h=lambda z:p.mpf('.2')-z,
            k=lambda z:p.mpf('.3')+z**3,e=lambda z:p.mpf('.4')+z,p=lambda z:p.mpf('.5')-z*z)
        a=lambda value:t.scalar(c.mpf(str(value)))
        values={key:a(fn(Z)) for key,fn in functions.items()};jets={key:a(p.diff(fn,Z)) for key,fn in functions.items()}
        C=lambda z:functions['k'](z)-A(z)*functions['m'](z)
        amp=dict(A=a(A(Z)),A_Z=a(p.diff(A,Z)),mu=a(mu),ratio=a(p.diff(A,Z)/A(Z)),
            log_A_lower=c.ln(c.mpf(str(A(Z)))),log_mu_lower=c.ln(c.mpf('.003')))
        targets=source.normalize(values,jets,a(C(Z)),a(p.diff(C,Z)),amp)
        reference=dict(M=lambda z:functions['m'](z)/A(z),I=lambda z:functions['h'](z)/A(z),
            S=lambda z:functions['e'](z)/A(z)**2,Cp=lambda z:functions['p'](z)/A(z)**2)
        reference[source.DROW]=lambda z:C(z)/(mu*A(z)**2)
        for key,fn in reference.items():
            compare.enclosed(p,c,targets['values'][key],fn(Z));normalization+=1
            compare.enclosed(p,c,targets['Z_derivatives'][key],p.diff(fn,Z));normalization+=1
    return dict(passed=True,independent_original_density_joint_C0_Z_slow_y_comparisons=checks,
        independent_frozen_reflection_joint_mean_comparisons=reflection,
        independent_exact_different_rate_Duhamel_identity_checks=transport,
        independent_five_normalized_target_quotient_derivative_comparisons=normalization,
        tiny_positive_mu_nonzero_joint_mismatch_checked=True,fixtures_are_not_actual_field_values=True)


def run():
    began=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    assert not manifest['actual_five_controls_or_functional_terminal_closure_admitted']
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    assert all(manifest[key] is False for key in source.common.current.FLAGS)
    c=MPIntervalContext();c.dps=240;tiles=0;axial=0;buffer=0;transition=0
    with mp.workdps(300):
        independent=independent_joint_checks(manifest['source_family'])
        for archive in manifest['actual_original_Rc_joint_terminal_defect_archives']:
            data=source.load_archive(archive);bd=source.load_archive(data['accepted_original_buffer_archive'])
            ad=source.load_archive(bd['accepted_original_axial_archive']);td=source.load_archive(bd['accepted_original_transition_archive'])
            assert data['source_family']==bd['source_family']==ad['source_family']==td['source_family']==manifest['source_family']
            assert data['exact_Z_range']==bd['exact_Z_range']==ad['exact_Z_range']==td['exact_Z_range']
            coordinates=source.native.HalfPstarCoordinates(c,iv(c,data['common_directed_coordinate_theorem']['common_log_bases'][1]),manifest['source_family'])
            replay=source.execute_tile(coordinates,bd,ad,td)
            expected={key:value for key,value in data.items() if key!='accepted_original_buffer_archive'}
            exact_replay_equal(json.loads(json.dumps(source.encode(replay))),expected,'new-joint-terminal-defect')
            assert data['original_Rc_amplitude_and_parameter']['original_mu']['sign']=='positive'
            assert not data['original_Rc_amplitude_and_parameter']['original_nonzero_relative_mu_increment']['exact_zero']
            assert data['joint_density_formed_before_axial_buffer_transition_branch_hulls']
            assert data['source_correlated_slope_inlet_and_full_function_oracle_still_missing']
            assert not data['actual_five_controls_or_functional_terminal_closure_admitted']
            controls=data['original_five_control_linear_response_enclosures']
            assert ep(iv(c,controls['original_exact_inverse_enclosure']['divided_axial_determinant']))[1]<0
            assert ep(iv(c,controls['original_exact_inverse_enclosure']['swirl_determinant']))[0]>0
            assert not controls['actual_nonlinear_control_functions_or_repair_admitted']
            axial+=data['actual_original_axial_joint_transport']['actual_original_source_joint_density_queries']
            buffer+=data['actual_original_buffer_joint_transport']['actual_joint_phase_queries']
            transition+=data['actual_original_transition_joint_transport']['actual_joint_phase_queries'];tiles+=1
            print('Actual source-correlated five Rc terminal defects replay:',data['exact_Z_range'],'PASS',flush=True)
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=sha(source.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        independent_joint_density_transport_and_normalization_checks=independent,
        actual_two_tile_target_and_original_linear_control_enclosure_replays=tiles,
        actual_axial_joint_density_queries_replayed=axial,actual_buffer_joint_density_queries_replayed=buffer,
        actual_transition_joint_density_queries_replayed=transition,
        original_source_family_units_true_width_phase_and_positive_mu_bound=True,
        actual_joint_transport_and_error_source_decomposition_replayed=True,
        original_inverse_solvers_accepted_upstream_or_density_producers_not_rerun=True,
        slope_component_correlation_and_original_function_control_oracle_still_open=True,
        actual_nonlinear_controls_functional_closure_whole_Z_axis_global_N_heat_stress_recursion_admitted=False,
        **dict.fromkeys(source.common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
