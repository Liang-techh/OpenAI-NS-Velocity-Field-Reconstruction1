"""Focused new slope algebra, full saved sum and sampled source replay.

No accepted producer or inverse solver is repeated. Independent finite
fixtures test the effective amplitude and different-rate transport; actual
source cells are restored with their original bases and saved inverse boxes.
"""
import json
from pathlib import Path
import time
from types import SimpleNamespace
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_slope_joint_source as source
import lei_ren_part1_paper_compliant_current_transition_normalized_q_source_check as compare
from lei_ren_part1_paper_compliant_current_transition_complete_prefix_check import exact_replay_equal

HERE,sha,ep,iv=source.HERE,source.sha,source.ep,source.iv


def independent_checks(family):
    c=MPIntervalContext();c.dps=100;p=mp.mp.clone();p.dps=350
    coordinates=source.native.HalfPstarCoordinates(c,c.mpf(0),family);t=coordinates.scalar(1)
    a=lambda x:t.scalar(c.mpf(str(x)))
    count=0;transport=0;basis=0
    for mu0 in ('.003','1e-300'):
        mu=p.mpf(mu0)
        for z0 in ('-.37','.37'):
            z=p.mpf(z0)
            for y0 in ('.1','.7','1'):
                y=p.mpf(y0)
                for phi0 in ('.17','.63'):
                    sn=p.sin(2*p.pi*p.mpf(phi0));J=lambda v:v*v/2
                    E=lambda v,w:p.exp(v/10-p.mpf('.6')*J(v))/(1+w*w)
                    V=lambda v,w:(p.mpf('.1')+p.mpf('.02')*v)*w
                    A=lambda v,w:(p.mpf('.02')+p.mpf('.003')*v+p.mpf('.04')*w)*sn
                    B=lambda v,w:(p.mpf('.05')+p.mpf('.005')*v-p.mpf('.03')*w)*sn
                    # Derive a_eff from A_Rc and exact remaining distance,
                    # independently of the producer's simplified formula.
                    logP=p.mpf(11)
                    Aend=lambda w:p.exp(-logP/2-p.mpf('1.2')-p.mpf('2.5')*mu)/(1+w*w)
                    ae=lambda v,w:Aend(w)*p.exp((logP+3-v)/2)
                    def direct(v,w):
                        de=E(v,w)*p.expm1(A(v,w)/source.N);dv=B(v,w)/source.N
                        return V(v,w)*de+E(v,w)*dv+de*dv-ae(v,w)*dv
                    primitives=dict(A=a(A(y,z)),B_over_Pstar=a(B(y,z)))
                    kernel=SimpleNamespace(c=c,scalar=t.scalar,q=t,
                        primitives=lambda coordinate,chart:primitives)
                    roots=dict(E={(0,0):a(E(y,z)),(0,1):a(p.diff(lambda w:E(y,w),z))},
                        V={(0,0):a(V(y,z)),(0,1):a(p.diff(lambda w:V(y,w),z))})
                    phase=dict(original_inverse=dict(coordinate_interval=source.encode(c.mpf('.2')),chart='psi'),
                        native_original_primitive_Z_enclosures=source.encode(dict(A_Z_slow=a(p.mpf('.04')*sn).record(),B_Z_slow=a(-p.mpf('.03')*sn).record())))
                    got=source.original_joint_density(kernel,roots,c.mpf(y0),c.mpf(z0),phase,c.mpf(mu0))
                    compare.enclosed(p,c,got['C'],direct(y,z));count+=1
                    compare.enclosed(p,c,got['C_Z'],p.diff(lambda w:direct(y,w),z));count+=1
                    eps=got['record']['original_nonzero_relative_mu_increment']
                    assert not eps['exact_zero'] and ep(eps['coefficient_interval'])[1]<0
    # Direct closed-form integrals with distinct recovery rates.
    for L0 in ('.1','1','3'):
        L=p.mpf(L0);D=p.mpf(2);gamma=p.mpf('.1');E0=p.mpf('1.3');V=p.mpf('.2');de0=E0*p.expm1(p.mpf('.003'));dv=p.mpf('.0007');ae0=p.mpf('1.4')
        Aend=ae0*p.exp(-(L+D)/2)
        ownk=(V*de0+(E0+de0)*dv)*p.exp(-p.mpf('1.5')*(L+D))*p.expm1((gamma+p.mpf('1.5'))*L)/(gamma+p.mpf('1.5'))
        ownm=dv*(-p.expm1(-L))*p.exp(-D)
        joint=p.exp(-p.mpf('1.5')*(L+D))*((V*de0+(E0+de0)*dv)*p.expm1((gamma+p.mpf('1.5'))*L)/(gamma+p.mpf('1.5'))-ae0*dv*p.expm1(L))
        assert p.almosteq(ownk-Aend*ownm,joint,rel_eps=p.mpf('1e-320'));transport+=1
    # Exact native-to-half power collection, including delta, R and u slots.
    for powers in ((1,2,3,1,-1),(-2,0,1,-2,2),(0,-1,0,0,1)):
        logP=c.mpf(7);logC=c.mpf(2);y=c.mpf('.3');logL=c.mpf('-.2');logu=c.mpf(-9)
        bases=(logP,-4*logP-30,logL,logu,c.ln(110)+10*(logC+logP)+y)
        old=source.prior.ScaledEnclosure(source.prior.FormalScale(bases,powers,c.mpf('.17')),c.mpf(('-2','3')),coordinates.ledger)
        target=source.native.HalfPstarCoordinates(c,logP*2,family)
        got=source.native_to_half(old,target,logC,y)
        assert got.scale.powers==(0,0,0,2*(powers[0]-4*powers[1]+10*powers[4]),0)
        lo,hi=ep(old.scale.evaluate());gl,gh=ep(got.scale.evaluate());assert max(lo,gl)<=min(hi,gh)
        assert got.coefficient._mpi_==old.coefficient._mpi_;basis+=1
    # Microscopic local log-|u| factors must survive expm1 as nonzero.
    tiny=source.prior.ScaledEnclosure(source.prior.FormalScale(t.scale.bases,offset=c.mpf(-2000)),c.mpf(('-1','1')),t.ledger)
    result=source.native_factored_expm1(tiny)
    assert not result.zero and result.scale.powers==tiny.scale.powers and result.scale.offset._mpi_==tiny.scale.offset._mpi_
    return dict(passed=True,independent_joint_C0_Z_and_effective_amplitude_comparisons=count,
        independent_different_rate_closed_form_integrals=transport,native_to_half_power_collection_checks=basis,
        microscopic_native_expm1_factor_and_original_positive_mu_retained=True,fixtures_are_not_actual_field_values=True)


@source.rc.native.inlet.source_precision
def run():
    began=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    assert all(manifest[key] is False for key in source.common.current.FLAGS)
    original_manifest=json.loads((HERE/source.source.NAME).read_bytes())
    weighted_manifest=json.loads((HERE/source.weighted.NAME).read_bytes())
    c=MPIntervalContext();c.dps=240;tiles=0;cells=0;sampled=0;phase_queries=0
    independent=independent_checks(manifest['source_family'])
    for archive in manifest['actual_original_slope_joint_source_archives']:
        data=source.terminal.load_archive(archive)
        old=source.terminal.load_archive(data['accepted_original_terminal_archive'])
        bd=source.terminal.load_archive(old['accepted_original_buffer_archive']);ad=source.terminal.load_archive(bd['accepted_original_axial_archive'])
        assert data['source_family']==old['source_family']==manifest['source_family'] and data['candidate_N']==source.N
        assert data['exact_Z_range']==ad['exact_Z_range']
        coordinates=source.native.HalfPstarCoordinates(c,iv(c,data['common_directed_coordinate_theorem']['common_log_bases'][1]),manifest['source_family'])
        amp=source.terminal.amplitude(coordinates,ad)
        C=coordinates.scalar(0);CZ=coordinates.scalar(0);rows=data['original_complete_2048_cell_joint_source_rows']
        assert len(rows)==2048
        for index,row in enumerate(rows):
            assert row['original_source_cell_index']==index
            assert source.Fraction(row['exact_y_cell'][0])==source.Fraction(index,2048) and source.Fraction(row['exact_y_cell'][1])==source.Fraction(index+1,2048)
            C=C+source.terminal.complete.restore_half_source(row['actual_joint_C_contribution'],coordinates)
            CZ=CZ+source.terminal.complete.restore_half_source(row['actual_joint_C_Z_contribution'],coordinates)
            phase_queries+=len(row['actual_original_phase_joint_sources']);cells+=1
        exact_replay_equal(source.encode(C.record()),data['slope_own_joint_C_at_slope_exit'],'complete-slope-C-sum')
        exact_replay_equal(source.encode(CZ.record()),data['slope_own_joint_C_Z_at_slope_exit'],'complete-slope-CZ-sum')
        report=next(row for row in original_manifest['original_complete_same_N_source_incoming_O2_integral_refinements'] if row['ordered_source_cells']==2048 and source.weighted.same_Z(row['exact_Z_range'],data['exact_Z_range']))
        assert report['whole_original_density_and_density_Z_source_archives']==data['accepted_original_slope_source_archives']
        for original_archive in report['whole_original_density_and_density_Z_source_archives']:
            original=source.terminal.load_archive(original_archive);left,right=original_archive['source_cell_index_range']
            for index in (left,right-1):
                got=source.source_row(coordinates,original['records'][index-left],index,iv(c,ad['original_logC']),iv(c,ad['original_dstar_log']),amp['mu_interval'])
                exact_replay_equal(source.encode(got['record']),rows[index],'saved-original-slope-source');sampled+=1
        pressure=next(row for row in weighted_manifest['actual_weighted_pressure_tiles'] if source.weighted.same_Z(row['exact_Z_range'],data['exact_Z_range']))
        primary=next(row for row in weighted_manifest['genuine_original_O2_weighted_pressure_replays'] if row['ordered_source_cells']==2048 and source.weighted.same_Z(row['exact_Z_range'],data['exact_Z_range']))
        trace=source.pre_slope_attribution(coordinates,source.terminal.load_archive(pressure['original_source_archive']),pressure)
        exact_replay_equal(source.encode(trace),data['actual_pre_slope_C1_contribution_attribution'],'pre-slope-source-attribution')
        assert all(row['label']=='active_first_bridge' for row in trace['dominant_inherited_Z_source'].values())
        composed=source.terminal_composition(coordinates,old,bd,ad,report,primary,C,CZ,amp,trace)
        exact_replay_equal(source.encode(composed),data['actual_updated_Rc_terminal_defects_and_source_ledger'],'new-Rc-target-composition')
        assert composed['actual_pre_slope_root_correlation_not_recovered']
        assert not composed['actual_nonlinear_controls_or_terminal_function_identity_admitted']
        assert data['slope_own_is_separate_from_actual_upstream_incoming']
        assert all(data[key] is False for key in source.common.current.FLAGS);tiles+=1
        print('Original slope joint source, separated memory and updated Rc targets:',data['exact_Z_range'],'PASS',flush=True)
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=sha(source.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        independent_slope_algebra_basis_and_transport_checks=independent,actual_complete_source_cells_summed=cells,
        actual_saved_original_source_cells_replayed=sampled,actual_source_phase_queries_bound=phase_queries,
        actual_terminal_composition_and_upstream_attribution_replays=tiles,
        accepted_inverse_solvers_and_upstream_producers_not_rerun=True,
        actual_pre_slope_function_correlation_nonlinear_controls_whole_Z_axis_global_N_heat_stress_recursion_admitted=False,
        **dict.fromkeys(source.common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8');return result


if __name__=='__main__':run()
