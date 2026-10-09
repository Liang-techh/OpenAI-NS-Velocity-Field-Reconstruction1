"""Independent active/join integral references and live whole-Rm source replay."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_Rm_whole_density_integrals as current
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow

fields,ep=current.fields,current.ep


def finite_fixture():
    c=MPIntervalContext();c.dps=200;p=mp.mp.clone();p.dps=240
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    P0=[f.scalar(1),f.scalar('.2')];Rm=f.scalar(7)
    active=((1,1),(5,4),current.JOIN);terminal=(current.JOIN,(2,1),'Rh');comparisons=0
    def pc(value):return p.e if value=='Rh' else p.mpf(value[0])/value[1]
    def reference(left,right,target,rate,Z,sign,derivative=False):
        yl,yr=p.log(pc(left)),p.log(pc(right));w=yr-yl
        suffix=p.log(pc(target))-yr;r=p.mpf(rate.numerator)/rate.denominator
        mass=w if not rate else -p.expm1(-r*w)/r
        first_u=w*w/2 if not rate else (1-(1+r*w)*p.exp(-r*w))/(r*r)
        first_s=w*mass-first_u
        constant=3+2*yl if derivative else 2+yl+Z*(3+2*yl)
        slope=2 if derivative else 1+2*Z
        return sign*p.exp(-r*suffix)*(constant*mass+slope*first_s)
    def packet(left,right,target,Z,sign):
        xl,xr=current.coordinate(c,left),current.coordinate(c,right)
        y=c.mpf([ep(c.ln(xl))[0],ep(c.ln(xr))[1]])
        kernels={key:f.scalar(sign*(2+y+c.mpf(str(Z))*(3+2*y))) for key in current.RATES}
        jets={key:f.scalar(sign*(3+2*y)) for key in current.RATES}
        geometry=current.previous.exact_geometry(left,right)
        source=dict(source_family='finite affine diagnostic',actual_closed_radial_source_cell=True,
            actual_original_spatial_Z_density_interface_installed=True,
            actual_original_Rm_radius_phase=dict(source_geometry=geometry,candidate_N=257,
                exact_source_Rm_factor=Rm,actual_Rm_radius_phase_Z_independent=True,
                source_radial_logarithmic_cell_width=c.ln(xr/xl)),
            actual_source_bound_phase_density_cells=[dict(source_family='finite affine diagnostic',
                source_geometry=geometry,original_common_P0_axial5=P0,
                exact_same_shared_radius_factor=Rm*c.mpf([ep(xl)[0],ep(xr)[1]]),
                actual_original_five_signed_density_C0_Z=dict(kernels=kernels,Z_derivatives=jets))])
        return current.integrate_source_cell(f,source,left,right,target,
            source_family='finite affine diagnostic',P0=P0,Rm_factor=Rm)
    for Z in (p.mpf('.3'),p.mpf('-.4')):
        for sign in (-1,1):
            totals=[]
            for partition,target in ((active,current.JOIN),(terminal,'Rh')):
                packets=[]
                for left,right in zip(partition,partition[1:]):
                    got=packet(left,right,target,Z,sign);packets.append(got)
                    for key,rate in current.RATES.items():
                        for n,row in enumerate(got['actual_cell_to_target_integral_C0_Z'][key]):
                            expected=reference(left,right,target,rate,Z,sign,n==1)
                            lo,hi=ep(row.finite_interval());assert lo<=expected<=hi;comparisons+=1
                totals.append({key:[sum((row['actual_cell_to_target_integral_C0_Z'][key][n] for row in packets),
                                       f.scalar(0)) for n in range(2)] for key in current.RATES})
            tail_width=c.mpf(1)-c.ln(current.coordinate(c,current.JOIN))
            whole=current.previous.affine_transport(f,totals[0],totals[1],tail_width)
            incoming={key:[f.scalar((i+1)*sign),f.scalar((i+2)*sign)] for i,key in enumerate(current.RATES)}
            carried=current.previous.affine_transport(f,incoming,whole,c.mpf(1))
            for i,(key,rate) in enumerate(current.RATES.items()):
                r=p.mpf(rate.numerator)/rate.denominator
                for n in (0,1):
                    expected=sum(reference(a,b,'Rh',rate,Z,sign,n==1)
                                 for partition in (active,terminal) for a,b in zip(partition,partition[1:]))
                    lo,hi=ep(whole[key][n].finite_interval());assert lo<=expected<=hi;comparisons+=1
                    expected+=(i+1+n)*sign*p.exp(-r)
                    lo,hi=ep(carried[key][n].finite_interval());assert lo<=expected<=hi;comparisons+=1
    current.validate_active_partition(c,current.ACTIVE_PARTITION)
    rejected=0
    for partition in (((1,1),current.JOIN),current.ACTIVE_PARTITION[:-1],
                      ((1,1),(51,40),(49,40),(59,40),(61,40),(69,40),current.JOIN)):
        try:current.validate_active_partition(c,partition)
        except ValueError:rejected+=1
        else:raise AssertionError('Incomplete or nonordered source atlas accepted')
    for value in ((9,10),(3,1)):
        try:current.coordinate(c,value)
        except ValueError:rejected+=1
        else:raise AssertionError('Outside actual Rm chart accepted')
    return dict(passed=True,independent_active_join_terminal_and_nonzero_incoming_comparisons=comparisons,
        source_atlas_endpoint_and_edge_guard_rejections=rejected,
        original_own_five_rates_and_zero_pressure_decay_checked=True,
        affine_composition_compared_to_direct_exact_integrals=True,
        finite_manufactured_sources_are_diagnostic_only=True)


def native(owner,report):
    cells=cell_rows=active_rows=whole_rows=memory=unknown=0;signs={};refinements={};statuses={};unknown_geometry={}
    for label in ('0','.5'):
        packet=owner.contribution(label)
        assert current.base.encoded(fields.serialized(packet))==report['frames'][label]
        op=owner.phase.upstream.upstream.upstream.owner(label).op;f=op.flow;c=op.c
        assert packet['exact_common_P0_axial5'] is op.P0
        partition=packet['exact_actual_active_partition']
        current.validate_active_partition(c,partition)
        atlas=packet['complete_active_source_atlas'];resolved=packet['actual_active_cells_to_join']
        unresolved=packet['actual_unresolved_active_cells']
        assert len(partition)==len(atlas)+1==len(resolved)+len(unresolved)+1
        seen_resolved=[];seen_unresolved=[]
        for left,right,entry in zip(partition,partition[1:],atlas):
            assert entry['source_geometry']==current.previous.exact_geometry(left,right)
            if entry['density_integral_enclosed']:
                index=entry['resolved_cell_index'];seen_resolved.append(index)
                assert entry['source_geometry']==resolved[index]['source_geometry']
            else:
                index=entry['unresolved_cell_index'];seen_unresolved.append(index)
                assert entry['source_geometry']==unresolved[index]['source_geometry']
        assert sorted(seen_resolved)==list(range(len(resolved)))
        assert sorted(seen_unresolved)==list(range(len(unresolved)))
        for cell in resolved:
            assert cell['target_coordinate']==current.JOIN and cell['source_family']==owner.family
            assert cell['actual_phase_source']['exact_source_Rm_factor'] is op.Rm_factor
            assert cell['actual_phase_source']['actual_Rm_radius_phase_Z_independent']
            assert ep(cell['source_radial_logarithmic_cell_width'])[0]>0
            assert ep(cell['source_logarithmic_suffix_to_target'])[0]>=0
            assert cell['no_extra_R_Pstar_or_N_factor'] and cell['source_cells_and_phase_unions_not_samples']
            for key,rate in current.RATES.items():
                assert ep(cell['original_positive_own_rate_masses'][key])[0]>0
                if not rate:assert ep(cell['original_own_rate_suffix_decays'][key])==(1,1)
                for row in cell['actual_cell_to_target_integral_C0_Z'][key]:
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger
                    assert all(mp.isfinite(v) for v in ep(row.coefficient));cell_rows+=1
            assert all(cell[key] is False for key in fields.previous.OPEN);cells+=1
        unknown_geometry[label]=[cell['source_geometry'] for cell in unresolved]
        for cell in unresolved:
            assert cell['exact_common_P0_axial5'] is op.P0 and cell['source_family']==owner.family
            assert cell['unresolved_integral_not_zero_or_replaced_by_a_sample']
            source=cell['actual_source_bound_unresolved_phase_density']
            assert source['source_family']==owner.family and source['actual_closed_radial_source_cell']
            assert not source['actual_original_spatial_Z_density_interface_installed']
            assert source['actual_original_Rm_radius_phase']['source_geometry']==cell['source_geometry']
            assert source['actual_original_Rm_radius_phase']['exact_source_Rm_factor'] is op.Rm_factor
            assert source['actual_source_bound_phase_density_cells']
            assert any(row['actual_original_five_signed_density_C0_Z'] is None
                       for row in source['actual_source_bound_phase_density_cells'])
            for row in source['actual_source_bound_phase_density_cells']:
                assert row['original_common_P0_axial5'] is op.P0 and row['source_family']==owner.family
                status=row['original_phase_Z_only_result']['status']
                assert status in {'requires_signed_source_refinement'},'New unresolved source status needs explicit scope: '+status
                statuses[status]=statuses.get(status,0)+1
                if status=='requires_signed_source_refinement':
                    u=row['original_phase_Z_only_result']['geometry']['original_u']
                    assert u['sign']=='undetermined' and not u['point_value_selected']
                    lo,hi=ep(u['coefficient_interval']);assert lo<=0<=hi
            unknown+=1
        assert packet['unknown_source_integral_count']==len(unresolved)
        assert packet['complete_whole_patch_density_integrals_installed']==(not unresolved)
        assert packet['complete_active_density_integrals_installed']==(not unresolved)
        assert packet['unresolved_integrals_are_explicit_unknown_source_terms']
        active=packet['actual_active_resolved_cell_integral_C0_Z']
        tail=packet['actual_terminal_local_defect_integral_C0_Z']
        whole=packet['actual_whole_patch_known_cell_contribution_C0_Z']
        assert ep(packet['full_original_logarithmic_interval_width'])==(1,1)
        expected=current.previous.affine_transport(f,active,tail,packet['actual_terminal_logarithmic_interval_width'])
        assert current.base.encoded(fields.serialized(expected))==current.base.encoded(fields.serialized(whole))
        assert ep(packet['original_terminal_incoming_to_Rh_decays']['p'])==(1,1)
        for key in current.RATES:
            for n in (0,1):
                value=sum((cell['actual_cell_to_target_integral_C0_Z'][key][n]
                           for cell in packet['actual_active_cells_to_join']),f.scalar(0))
                assert current.base.encoded(fields.serialized(value))==current.base.encoded(fields.serialized(active[key][n]))
                for row in (active[key][n],whole[key][n]):
                    assert row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger
                active_rows+=1;whole_rows+=1
        for name in ('actual_leading_Rm_memory','actual_leading_join_memory','actual_leading_Rh_memory'):
            histories=packet[name];assert set(histories)==set(current.RATES)
            for rows in histories.values():
                assert len(rows)==6
                assert all(row.ctx is c and row.scale.bases is f.logs and row.ledger is f.ledger for row in rows)
            memory+=1
        assert packet['finite_N_Rm_incoming_defect_is_unsupplied_affine_argument']
        assert packet['finite_N_complete_prefix_corrected_Rh_or_moment_closure_not_claimed']
        assert all(packet[key] is False for key in fields.previous.OPEN)
        try:owner.transport_supplied_incoming(label,None)
        except ValueError:pass
        else:raise AssertionError('Unspecified Rm correction silently reset')
        if unresolved:
            incoming={key:[f.scalar(1),f.scalar(1)] for key in current.RATES}
            try:owner.transport_supplied_incoming(label,incoming)
            except ValueError:pass
            else:raise AssertionError('Unknown actual source integrals silently omitted')
        signs[label]={key:[row.record()['sign'] for row in rows] for key,rows in whole.items()}
        refinements[label]=len(packet['actual_source_refinements'])
    return dict(passed=True,actual_resolved_active_closed_source_cells=cells,
        actual_active_cell_to_join_C0_Z_rows=cell_rows,actual_active_resolved_C0_Z_integral_rows=active_rows,
        actual_whole_patch_known_contribution_C0_Z_rows=whole_rows,actual_leading_source_memory_packets=memory,
        unresolved_actual_active_source_cells=unknown,unresolved_original_phase_status_counts=statuses,
        unresolved_exact_source_geometries=unknown_geometry,
        signed_source_interval_crossing_not_a_certified_exact_root_or_singularity=True,
        actual_source_refinement_counts=refinements,conditional_whole_patch_integral_signs=signs,
        exact_support_edge_atlas_and_active_terminal_suffix_composition_checked=True,
        real_Rm_incoming_prefix_whole_Z_moment_closure_repair_and_global_N_not_admitted=True)


def run():
    began=time.monotonic()
    with mp.workdps(260):fixture=finite_fixture()
    print('Independent whole-Rm dx/x and nonzero inlet composition PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalRmWholeDensityIntegrals(require_checked=False);evidence=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    assert report['original_integral_chart_binding']==current.INTEGRATION_BINDING
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_exact_active_terminal_affine_integral_fixture=fixture,actual_live_source=evidence,
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual Rm active atlas, known density integrals and explicit unknown cells PASS',flush=True);return result


if __name__=='__main__':run()
