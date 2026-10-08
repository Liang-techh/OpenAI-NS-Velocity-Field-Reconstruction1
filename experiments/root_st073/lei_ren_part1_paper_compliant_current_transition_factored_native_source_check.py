"""Focused original source arithmetic, saved native replay and local integrals.

The finite fixtures execute the unmodified original operator at independently
chosen physical s points. Actual saved native parents are replayed without
constructing native owners or repeating any accepted upstream integrations.
"""
import gzip
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_transition_factored_native_source as source
import lei_ren_part1_paper_compliant_current_generic_tail_phase_integrals_check as saved
import lei_ren_part1_paper_compliant_current_transition_normalized_q_source_check as qcheck

HERE,ep,iv,sha=source.HERE,source.ep,source.iv,source.sha
same=saved.same


def check_scaled_tree(c,record,value,path='source'):
    if isinstance(value,source.prior.ScaledEnclosure):
        try:same(c,record,value)
        except AssertionError as exc:raise AssertionError(path) from exc
        return 1
    if isinstance(value,source.ScaledTaylor):return sum(check_scaled_tree(c,r,v,path+'/'+str(i)) for i,(r,v) in enumerate(zip(record,value.coefficients,strict=True)))
    if isinstance(value,dict):
        assert set(record)==set(value)
        return sum(check_scaled_tree(c,record[k],v,path+'/'+str(k)) for k,v in value.items())
    if isinstance(value,(list,tuple)):return sum(check_scaled_tree(c,r,v,path+'/'+str(i)) for i,(r,v) in enumerate(zip(record,value,strict=True)))
    assert record==value
    return 0


def restore_parent(c,record):
    histories={key:[IntervalTaylor(c,[iv(c,v) for v in rows]) for rows in jets]
        for key,jets in record['original_history_axial_Taylor_coefficients'].items()}
    return dict(actual_normalized_primitive_y_derivative_axial5=histories,
        original_P0_axial5_coefficients=[iv(c,v) for v in record['original_P0_axial_Taylor_coefficients']],
        log_Utheta_over_Pstar_base_source=iv(c,record['original_log_Utheta_over_Pstar_base_source']))


def record_backed_adapter(c,coordinates,record,parent):
    adapter=object.__new__(source.FactoredOriginalTransition)
    adapter.c=c;adapter.coordinates=coordinates;adapter.template=coordinates.scalar(1)
    adapter.mu=iv(c,record['original_mu']);adapter.eta=iv(c,record['original_eta_log'])
    adapter.owner=SimpleNamespace(family=coordinates.family)
    adapter.proof=source.canonical_amplitude_proof();adapter.parents={tuple(record['exact_Z_range']):parent}
    adapter.pre=SimpleNamespace(params=SimpleNamespace(mu=adapter.mu),delta=iv(c,record['original_delta']))
    return adapter


def original_finite_checks(family):
    c=MPIntervalContext();c.dps=85
    direct=MPIntervalContext();direct.dps=150
    p=mp.mp.clone();p.dps=350
    count=0;sigma_checks=0;kernel_checks=0
    for D in (64,128):
        coordinates=source.HalfPstarCoordinates(c,c.mpf(32),family)
        mu=c.exp(c.ln(c.mpf('.001'))-D);eta=c.ln(mu)-D
        adapter=object.__new__(source.FactoredOriginalTransition)
        adapter.c=c;adapter.coordinates=coordinates;adapter.template=coordinates.scalar(1)
        adapter.mu=mu;adapter.eta=eta;adapter.owner=SimpleNamespace(family=family)
        adapter.proof=source.canonical_amplitude_proof();adapter.parents={}
        def parent(ctx,Z):
            z=IntervalTaylor.variable(ctx,ctx.mpf(Z),5);qi=(1+z*z).reciprocal()
            ud=qi*ctx.exp(ctx.mpf('.3')-8)
            history={k:[IntervalTaylor(ctx,[ctx.mpf('.1'),ctx.mpf('.01')]+[ctx.mpf(0)]*4)]*5 for k in source.rc.RATES}
            return dict(actual_normalized_primitive_y_derivative_axial5=history,
                Utheta_over_Pstar_axial5_coefficients=ud.coefficients,
                log_Utheta_over_Pstar_base_source=ctx.mpf('.3')-8-ctx.ln(1+z[0]**2),
                original_P0_axial5_coefficients=[ctx.mpf('.3'),ctx.mpf('.02')]+[ctx.mpf(0)]*4)
        adapter.pre=SimpleNamespace(params=SimpleNamespace(mu=mu),delta=c.mpf('.2'),
            axial=lambda Z,buffer_offset:parent(c,Z))
        coord=source.previous.OriginalTransitionCoordinates(coordinates,mu,eta)
        geometry=coord.geometry('xi','1/2','3/4');got=adapter.source(('.37','.37'),geometry,8)
        for numerator in (4,5,6):
            s=p.mpf(numerator)/(8*p.sqrt(D))
            sig=lambda x:1/(1+p.exp(1/x**2-1/(1-x)**2))
            rows=source.original_sigma_rows(adapter.template,c.mpf(p.nstr(s,340)))
            for order,value in enumerate(rows):
                qcheck.enclosed(p,c,value,p.diff(sig,s,order));sigma_checks+=1
            # Execute the original operator with ordinary interval arithmetic,
            # an independent precision and exact physical source point.
            Z=('.37','.37');rawparent=parent(direct,Z)
            pre=SimpleNamespace(ctx=direct,params=SimpleNamespace(mu=direct.exp(direct.ln(direct.mpf('.001'))-D)),
                delta=direct.mpf('.2'),cells=8,axial=lambda ZZ,buffer_offset:rawparent,
                coordinates=lambda ZZ:(IntervalTaylor.variable(direct,direct.mpf(ZZ),5),None))
            emitted={}
            def packet(ZZ,chart,offset,u,logu,logU,V,history,extra):
                p0=IntervalTaylor(direct,rawparent['original_P0_axial5_coefficients'])
                emitted.update(source.original.physical_mixed(direct,ZZ,pre.delta,u,logU,V,history,p0,direct.exp(-32)))
                return emitted
            pre.packet=packet
            source.original.CompliantPrePulseMixedC4.slope_mu(pre,Z,direct.mpf(p.nstr(s,340)))
            physical=got['values']['physical']
            for name in ('physical_velocity_pressure_y_Z_mixed4','physical_five_primitive_y_Z_mixed4'):
                for key,grid in emitted[name].items():
                    for order,value in grid.items():
                        for endpoint in ep(value):qcheck.enclosed(p,c,physical[name][key][order],p.mpf(endpoint));count+=1
            reference=source.original.transition_kernels(direct,direct.mpf(p.nstr(s,340)),pre.params.mu,8)
            kernel=got['record']['original_safe_transition_kernels']['original_kernel_values']
            for key in ('J','theta','energy','pressure'):
                value=saved.native_restore(source.encode(kernel[key]),coordinates.bases,coordinates.ledger)
                for endpoint in ep(reference[key]):qcheck.enclosed(p,c,value,p.mpf(endpoint));kernel_checks+=1
    return dict(passed=True,independent_physical_s_original_sigma_y0_through_y4_comparisons=sigma_checks,
        independent_unmodified_original_physical_velocity_and_five_history_mixed4_endpoint_comparisons=count,
        independent_original_transition_kernel_endpoint_comparisons=kernel_checks,
        original_unmodified_operator_and_independent_precision_used=True,
        fixtures_not_substituted_for_actual_native_source=True)


def run():
    began=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    assert manifest['original_typed_local_five_signed_C0_and_Z_integral_contributions_executed']
    assert not manifest['full_prefix_density_integrals_or_Rc_targets_admitted']
    assert all(manifest[k] is False for k in source.common.current.FLAGS)
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    c=MPIntervalContext();c.dps=240
    count=0;phase_queries=0;integrals=0;query_summaries=[]
    with mp.workdps(300):
        independent=original_finite_checks(manifest['source_family'])
        inventory=json.loads((HERE/(source.PREFIX+'current_generic_shear_loop_jet_bounds.json')).read_bytes())
        positive=inventory['current_actual_loop_jet_log_bounds_by_chart']['O3_slope_mu']['actual_positive_denominator_theorem']
        scales=json.loads((HERE/(source.PREFIX+'current_generic_shear_uniform_inputs.json')).read_bytes())['current_actual_logarithmic_loop_scales']
        norm_family=json.loads((HERE/(source.PREFIX+'physical_norm_family.json')).read_bytes())
        pressure=json.loads((HERE/(source.PREFIX+'pressure_source.json')).read_bytes())['compliant_source']
        assert pressure['implicit_source_definition']['delta']=='min(1e-200,exp(-4logPstar-30))'
        assert norm_family['implicit_source_sha256']==manifest['source_family']['implicit_source_sha256']
        assert norm_family['datum_enclosure_sha256']==manifest['source_family']['datum_enclosure_sha256']
        for archive in manifest['original_factored_native_source_archives']:
            compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
            assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
            assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
            data=json.loads(raw)
            assert data['source_family']==manifest['source_family'] and data['candidate_N']==source.N
            Z=data['exact_Z_range'];basis=data['common_directed_coordinate_theorem']['common_log_bases']
            coordinates=source.HalfPstarCoordinates(c,iv(c,basis[1]),manifest['source_family'])
            assert source.encode(coordinates.record())==data['common_directed_coordinate_theorem']
            pl,ph=ep(coordinates.logP_squared);fl,fh=ep((c.exp(40)+11)*2)
            assert pl<=fl<=fh<=ph
            assert all(ep(iv(c,basis[j]))==(0,0) for j in (0,2,4))
            assert ep(iv(c,basis[3]))==ep(coordinates.logP_squared/4)
            assert ep(iv(c,data['original_logC']))==ep(iv(c,norm_family['selected_logCstar']))
            assert ep(iv(c,data['original_dstar_log']))==ep(iv(c,scales['logarithmic_selected_positive_lower_constants']['d_star']))
            assert ep(iv(c,data['original_positive_a_lower_log']))==ep(iv(c,positive['log_actual_a_positive_lower']))
            summary=[]
            for row in data['original_source_rows']:
                saved_source=row['original_typed_native_source'];inputs=saved_source['lossless_actual_native_parent_inputs']
                assert inputs['exact_Z_range']==Z
                parent=restore_parent(c,inputs);adapter=record_backed_adapter(c,coordinates,inputs,parent)
                assert ep(adapter.mu)==ep(iv(c,positive['actual_positive_mu']))
                assert ep(adapter.eta)==ep(iv(c,scales['selected_positive_eta_log']))
                delta=iv(c,inputs['original_delta']);expected=c.exp(-4*(c.exp(40)+11)-30)
                assert ep(delta)[0]<=ep(expected)[0]<=ep(expected)[1]<=ep(delta)[1]
                assert ep(expected)[1]<ep(c.mpf('1e-200'))[0]
                parameter=saved_source['actual_original_typed_geometry']['original_parameter_source']
                assert ep(adapter.mu)==ep(iv(c,parameter['original_mu'])) and ep(adapter.eta)==ep(iv(c,parameter['original_eta_log']))
                coord=source.previous.OriginalTransitionCoordinates(coordinates,adapter.mu,adapter.eta)
                g=saved_source['actual_original_typed_geometry']
                geometry=coord.geometry(g['typed_coordinate_kind'],*g['exact_normalized_endpoints'])
                assert source.encode(geometry['record'])==g
                cells=saved_source['original_safe_transition_kernels']['cells']
                got=adapter.source(Z,geometry,cells)
                assert source.encode(got['record']['same_original_source_program'])==saved_source['same_original_source_program']
                assert source.encode(got['record']['same_original_physical_mixed_program'])==saved_source['same_original_physical_mixed_program']
                assert source.encode(adapter.proof)==saved_source['canonical_original_amplitude_source']
                count+=check_scaled_tree(c,saved_source['actual_original_ordinary_velocity_history_and_pressure_rows'],got['values']['physical'])
                count+=check_scaled_tree(c,saved_source['original_current_amplitude'],got['values']['u'])
                assert source.encode(got['record']['original_safe_transition_kernels'])==saved_source['original_safe_transition_kernels']
                assert source.encode(got['record']['original_profile_exponential_corrections'])==saved_source['original_profile_exponential_corrections']
                for key in source.rc.RATES:
                    assert source.encode(got['record']['original_inherited_Rd_history_rows'][key])==saved_source['original_inherited_Rd_history_rows'][key]
                norm=coord.normalized_excess(adapter.template,geometry)
                seed=SimpleNamespace(core=SimpleNamespace(logC=iv(c,data['original_logC'])))
                roots=source.signed_root_source(adapter,got,norm['rows'],adapter.eta,iv(c,data['original_positive_a_lower_log']),seed)
                assert source.encode(roots['record'])==row['original_signed_root_source']
                assert source.encode(source.previous.correlated_original_q(norm['rows'],adapter.eta,iv(c,data['original_positive_a_lower_log']))['record'])==row['original_correlated_q_source']
                phase=row['actual_original_typed_phase']
                assert phase['actual_phase_is_N_times_original_s_not_N_times_xi_or_k']
                assert phase['nonzero_microscopic_phase_component_retained']
                assert phase['same_original_affine_radius_identity']['passed']
                same(c,phase['original_microscopic_s_component'],got['values']['original_s_source'])
                same(c,phase['original_microscopic_N_s_component'],got['values']['original_s_source']*source.N)
                assert phase['original_base_source_phase']['original_spatial_phase_bound_at_candidate_N']
                decoded_phase=dict(phase,actual_phase_boxes=[iv(c,box) for box in phase['actual_phase_boxes']])
                local=source.original_local_density_integrals(adapter,got,roots,geometry,decoded_phase,iv(c,data['original_dstar_log']))
                assert source.encode(local['inverses'])==row['actual_original_conditioned_first_jet_queries']
                assert source.encode(local['record'])==row['actual_original_local_five_C0_Z_density_and_integral']
                assert local['record']['status']=='enclosed'
                assert not local['record']['full_original_prefix_incoming_or_Rc_targets_admitted']
                assert row['original_full_prefix_integral_not_claimed']
                phase_queries+=len(local['inverses']);integrals+=2*len(source.rc.RATES)
                summary.append(dict(label=row['label'],p2_sign=roots['roots']['p2'][source.ZERO].record()['sign'],
                    inverse_statuses=[x['status'] for x in local['inverses']],local_five_C0_Z_integral_status=local['record']['status']))
            query_summaries.append(dict(exact_Z_range=Z,actual_source_queries=summary))
        assert query_summaries==manifest['actual_original_query_summaries']
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=sha(source.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        independent_original_arithmetic_checks=independent,
        actual_native_source_parent_physical_rows_replayed=count,
        actual_original_phase_first_jet_queries_replayed=phase_queries,
        original_typed_local_five_C0_Z_integral_contribution_rows_replayed=integrals,
        actual_parent_inputs_kept_lossless=True,accepted_upstream_native_integration_not_rerun=True,
        same_actual_mu_eta_dstar_pressure_delta_and_selected_logCstar_bound=True,
        full_original_prefix_axial_targets_controls_recursion_and_NS_not_claimed=True,
        **dict.fromkeys(source.common.current.FLAGS,False),input_hashes=hashes,
        execution_seconds=time.monotonic()-began)
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original factored native source/local five C0-Z integrals PASS; full prefix remains open',flush=True)
    return result


if __name__=='__main__':run()
