"""Conditional source restriction, shared q/u basis and full-domain evidence."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_predicate_source_frames as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

ep=current.ep;C0,Y,Z,YZ=current.C0,current.Y,current.Z,current.YZ


def source_restriction_calculus():
    p,q,d,py,pz,pyz,qy,gamma=s.symbols('p q d py pz pyz qy gamma',nonzero=True)
    u=p*q/d
    assert s.cancel(d*u/q-p)==0
    assert s.cancel((py*q+p*qy)/d-u*(gamma+qy/q)).subs(py,p*gamma)==0
    assert s.cancel(s.diff((p+s.Symbol('z')*pz)*q/d,s.Symbol('z'))-pz*q/d)==0
    yy,zz=s.symbols('y z')
    pp=p+yy*py+zz*pz+yy*zz*pyz;qq=q+yy*qy
    assert s.expand(s.diff(pp*qq/d,yy,zz).subs({yy:0,zz:0})-(pyz*q+pz*qy)/d)==0
    # q>0 and the accepted g never vanishes: these two identities extend
    # through p2=u=0 and do not use a reciprocal p2 or reciprocal u.
    return dict(passed=True,original_source_restriction_identities=4,
        only_positive_q_dstar_and_nonzero_carrier_divided=True,
        regular_axis_identity_extends_through_u_zero=True,
        predicate_union_proof='Every real u has abs(u)<=1/4 or u>=3/16 or u<=-3/16; overlap3/16<=abs(u)<=1/4',
        signed_q_lower_proof='|u|=|p2|*q/dstar<=Lambda0*max|g/Lambda0|*q/dstar, using |Z|<=1; hence signed abs(u)>=3/16 gives the recorded q lower',
        predicates_restrict_actual_source_function_not_entire_outer_rectangle=True)


def separate_log_arithmetic(owner):
    raw=owner.source.source_frame(64,0,Z_lower=-1,Z_upper=0);c=owner.ctx;p=mp.mp.clone();p.dps=100
    a=current.O2QUAtlas(owner.owner.inputs.frame,lower=-1,upper=0,
        logq=c.mpf((-10,-1)),logu=c.ln(c.mpf(['.1875','20'])))
    comparisons=0
    for qp in (-2,0,2):
        for up in (-3,0,3):
            left=current.prior.ScaledEnclosure(current.prior.FormalScale(a.bases,(0,0,0,qp,up)),c.mpf(['-.7','1.3']),a.ledger)
            right=current.prior.ScaledEnclosure(current.prior.FormalScale(a.bases,(0,0,0,-qp,-up)),c.mpf(['-.2','.8']),a.ledger)
            result=a.add(left,right).finite_interval()
            for lq in (-10,-1):
                for u in (p.mpf('.1875'),p.mpf(20)):
                    for lc,rc in ((p.mpf('-.7'),p.mpf('-.2')),(p.mpf('1.3'),p.mpf('.8'))):
                        target=lc*p.exp(lq*qp)*u**up+rc*p.exp(-lq*qp)*u**(-up)
                        saved.contains(result,target,p.mpf('1e-35'));comparisons+=1
    rejected=0
    other=current.O2QUAtlas(owner.owner.inputs.frame,lower=-1,upper=0,logq=c.mpf((-10,-1)),logu=c.mpf(0))
    for call in (lambda:a.add(a.scalar(1),other.scalar(1)),
            lambda:current.rebase_original_root(other,raw.roots['q'][C0],raw_frame=raw,source_owner=owner.source)):
        try:call()
        except ValueError:rejected+=1
    assert rejected==2
    return dict(passed=True,separate_q_u_formal_addition_endpoint_comparisons=comparisons,
        foreign_basis_and_invalid_q_rebase_guards=rejected,
        fifth_slot_is_defined_logabsu_not_radius=True)


def whole_original_predicate_contract(owner,manifest):
    c=owner.ctx;count=min(owner.source.parent.parent.levels);rows=manifest['whole_original_O2_predicate_source_cells']
    assert len(rows)==count and manifest['exact_outer_original_domain']==dict(y=['0','1'],Z=['-1','1'])
    assert manifest['all_real_original_u_including_zero_covered_by_predicates']
    rootrows=0;frames=0;signed=0
    for index,item in enumerate(rows):
        assert set(item)=={'regular','positive','negative'}
        for branch,record in item.items():
            assert record['source_family']==owner.family and record['branch']==branch
            assert record['exact_y_cell']==[str(s.Rational(index,count)),str(s.Rational(index+1,count))]
            assert record['exact_outer_Z_bounds']==list(map(str,dict(regular=(-1,1),positive=(-1,0),negative=(0,1))[branch]))
            assert record['conditional_domain_not_claimed_entire_outer_rectangle']
            assert record['same_bases_context_and_ledger_for_all_rebuilt_roots']
            assert record['original_pressure_family_and_incoming_histories_unchanged']
            assert record['original_p2_Z_yZ_template_and_pressure_error_rows_retained']
            atlas=record['atlas'];assert atlas['source_basis_order'][-2:]==['original_logq','defined_original_logabsu_or_zero']
            assert atlas['actual_radius_power_is_exact_zero_and_not_in_fifth_slot']
            oq=saved.interval(c,record['original_positive_logq_before_predicate']);pq=saved.interval(c,record['predicate_logq'])
            assert ep(oq)[0]<=ep(pq)[0]<=ep(pq)[1]<=ep(oq)[1]
            for name,row in record['actual_root_jets'].items():
                assert set(row)==set(map(str,(C0,Y,Z,YZ)))
                for value in row.values():
                    assert value['encloses_original_source_function'] and not value['point_value_selected'];rootrows+=1
                if name in ('a','q','nu'):assert row[str(Z)]['exact_zero'] and row[str(YZ)]['exact_zero']
            assert not record['actual_root_jets']['q'][str(C0)]['exact_zero']
            assert not record['actual_root_jets']['q'][str(Y)]['exact_zero']
            assert not record['actual_root_jets']['nu'][str(Y)]['exact_zero']
            if branch=='regular':
                assert ep(oq)==ep(pq)
                u=record['actual_u_mixed'][str(C0)]
                saved.contains(saved.saved_value(c,u,bases=tuple(saved.interval(c,v) for v in atlas['defining_basis'])),c.mpf(['-.25','.25']))
                assert record['kernel_geometry']['branch']=='small_r_series'
            else:
                signed+=1
                implication=saved.interval(c,record['signed_predicate_q_lower_implication'])
                assert ep(pq)[0]>=ep(implication)[0]
                logu=saved.interval(c,atlas['defining_basis'][4])
                assert ep(logu)[0]>=ep(c.ln(c.mpf(3)/16))[0]
                assert record['actual_u_mixed'][str(C0)]['formal_positive_scale']['radius_power']==1
                assert record['kernel_geometry']['branch']=='signed_Mobius'
            frames+=1
    # Fresh issued identity/closed branch guards, not an ancestor replay.
    f=owner.frame(count,0,branch='regular');a=f.roots['q'].atlas;rejected=0
    for call in (lambda:owner.describe(replace(f)),lambda:owner.frame(count,0,branch='unsigned'),
            lambda:current.O2QUAtlas(owner.owner.inputs.frame,lower=-1,upper=0,logq=a.bases[3],logu=c.mpf((-3,0)))):
        try:call()
        except ValueError:rejected+=1
    assert rejected==3 and owner.describe(f) is f.record
    return dict(passed=True,continuous_original_y_cells=count,issued_conditional_source_frames=frames,
        strict_signed_source_frames=signed,actual_mixed_root_rows=rootrows,issued_identity_and_predicate_guards=rejected,
        all_real_u_and_both_outer_Z_sides_including_axis_covered=True,
        q_nu_a_varying_y_rows_and_original_transverse_pressure_rows_retained=True,
        conditional_restriction_not_selected_source_and_not_full_rectangle_geometry_claim=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('actual_mixed_primitives_installed','actual_changed_five_integrals_installed','all_17_chart_or_24_cell_oracle_installed',
        'actual_five_controls_installed','current_whole_N_selected',*current.source.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2PredicateSources();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (('conditional_original_source_calculus',source_restriction_calculus),
                ('separate_original_q_u_basis_arithmetic',lambda:separate_log_arithmetic(owner)),
                ('whole_original_O2_issued_predicate_domain',lambda:whole_original_predicate_contract(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            uncompressed_bytes=len(raw),compressed_bytes=path.stat().st_size),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Issued conditional full original O2 source frame interface only: original q/u log roles, all mixed roots, sign-correlated y rows, full named predicate cover. Mixed primitives/integrals, global N and full reconstruction remain open.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('Full original O2 conditional source interface PASS',flush=True);return report


if __name__=='__main__':run()
