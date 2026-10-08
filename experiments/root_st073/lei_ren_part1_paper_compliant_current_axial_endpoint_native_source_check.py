"""Focused independent endpoint ODE checks and lossless whole-axial replay.

No upstream owners or accepted quadrature producers are rerun. Finite M=1
fixtures are separate from the actual M=40/N1024 strict-sign source records.
"""
from fractions import Fraction
import gzip
import hashlib
import json
import math
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_axial_endpoint_native_source as source
import lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator as stress
import lei_ren_part1_paper_compliant_current_transition_normalized_q_source_check as compare
from lei_ren_part1_paper_compliant_current_transition_complete_prefix_check import exact_replay_equal

HERE,sha,ep,iv=source.HERE,source.sha,source.ep,source.iv


def independent_original_endpoint_checks(family):
    c=MPIntervalContext();c.dps=100;direct=MPIntervalContext();direct.dps=160
    p=mp.mp.clone();p.dps=350
    coordinates=source.native.HalfPstarCoordinates(c,2*(c.exp(1)+11),family)
    coord=source.axial.AxialCoordinates(coordinates,-c.mpf(512)-coordinates.logP_squared,1)
    Z=['.37','.37'];z=IntervalTaylor.variable(direct,direct.mpf(Z),5);qi=(1+z*z).reciprocal()
    logP=direct.exp(1)+11;P=direct.exp(logP);EY=qi*direct.exp(direct.mpf('5.8')-logP/2)
    endpoint=dict(m=1+z,h=qi*direct.mpf('.03'),k=qi*z*direct.mpf('.02'),
        e=z*z*direct.mpf('.001'),p=qi*qi*direct.mpf('.1'))
    raw=dict(m=endpoint['m']*direct.exp(-11),
        h=endpoint['h']*direct.exp(-direct.mpf('16.5'))+EY*(direct.exp(-direct.mpf('5.5'))-direct.exp(-direct.mpf('16.5'))),
        k=endpoint['k']*direct.exp(-direct.mpf('16.5')),
        e=(endpoint['e']-EY*EY*direct.mpf('5.5'))*direct.exp(-11),
        p=endpoint['p']+EY*EY*((1-direct.exp(-11))/2))
    P0=z*z*direct.mpf('.03')+direct.mpf('.2')
    inputs=dict(original_delta=c.mpf('.001'),original_P0_axial_Taylor_coefficients=P0.coefficients,
        original_history_axial_Taylor_coefficients={k:[v.coefficients] for k,v in raw.items()})
    data=dict(exact_Z_range=Z,original_source_rows=[dict(original_typed_native_source=dict(lossless_actual_native_parent_inputs=inputs))])
    adapter=source.OriginalAxialEndpoint(coordinates,coord,data)
    _,mass=source.native.original.slope_masses(direct,direct.mpf(1),2048)
    left_h=qi*((direct.mpf(5)/8+mass[0])*direct.exp(-direct.mpf('1.5')))
    left=dict(m=z*4,h=left_h,k=left_h*z*4,
        e=z*z*(16/P**2)-qi*qi*((direct.mpf(5)/12+mass[2]/2)*direct.exp(-1)),
        p=qi*qi*(direct.mpf('2.5')+mass[1]/2))
    for key in source.KEYS:
        for n in range(6):
            for x in ep(endpoint[key][n]):compare.enclosed(p,c,adapter.right[key][n],p.mpf(x))
    count=60;root_checks=0;width_checks=0
    def sigma(r):
        if r<=0:return p.mpf(0)
        if r>=1:return p.mpf(1)
        return 1/(1+p.exp(1/r**2-1/(1-r)**2))
    for side in ('left','right'):
        for kind,a,b in (('xi','0','3/4'),('mixed','3/4','8'),('k','0','-3')):
            geometry=coord.geometry(side,kind,a,b);got=adapter.source(geometry)
            norm=source.axial.normalized_excess(coord,Z,geometry);q=source.axial.original_axial_q(coord,norm['rows'])
            roots=source.original_signed_roots(adapter,got,norm,q,geometry,c.mpf('.1'))['roots']
            ra,rb=(p.mpf(x) for x in ep(geometry['cover']))
            ya,yb=(p.exp(ra),p.exp(rb)) if side=='left' else (p.exp(1-rb),p.exp(1-ra))
            compare.enclosed(p,c,geometry['physical_width'],yb-ya);width_checks+=1
            for r in ((ra+rb)/2,rb):
                y=p.exp(r) if side=='left' else p.exp(1-r)
                distance=y-1 if side=='left' else p.exp(1)-y
                # Independent quadrature of the original weighted B source,
                # followed by exact forward/backward variation of constants.
                if distance:
                    points=[0,distance/2,9*distance/10,distance]
                    phase=lambda s:p.log(1+s) if side=='left' else -p.log(1-s/p.exp(1))
                    J=p.quad(lambda s:p.exp(s if side=='left' else -s)*sigma(phase(s)),points)
                    J2=p.quad(lambda s:p.exp(s if side=='left' else -s)*(2*sigma(phase(s))-sigma(phase(s))**2 if side=='left' else sigma(phase(s))**2),points)
                else:J=J2=p.mpf(0)
                d=direct.mpf(p.nstr(distance,340));j=direct.mpf(p.nstr(J,340));j2=direct.mpf(p.nstr(J2,340))
                EL=qi*direct.exp(-direct.mpf('.2'))
                if side=='left':
                    dec=direct.exp(-d);d3=direct.exp(-direct.mpf('1.5')*d)
                    hist=dict(m=z*4-z*(4*dec*j),
                        h=left['h']*d3+EL*(direct.exp(-d/2)-d3),
                        k=left['k']*d3+EL*z*4*d3*(direct.expm1(d)-j),
                        e=left['e']*dec+z*z*(16/P**2)*dec*(direct.expm1(d)-j2)-EL*EL*d*dec/2,
                        p=left['p']+EL*EL*(1-dec)/2)
                    E=EL*direct.exp(-d/2)
                else:
                    dec=direct.exp(d);d3=direct.exp(direct.mpf('1.5')*d)
                    hist=dict(m=(endpoint['m']-z*4*j)*dec,
                        h=(endpoint['h']-EY*(1-direct.exp(-d)))*d3,
                        k=(endpoint['k']-EY*z*4*j)*d3,
                        e=(endpoint['e']-z*z*(16/P**2)*j2+EY*EY*d/2)*dec,
                        p=endpoint['p']-EY*EY*direct.expm1(d)/2)
                    E=EY*direct.exp(d/2)
                yy=direct.mpf(p.nstr(y,340))
                # Reflection keeps derivatives of the tiny sigma tail
                # separate from its baseline 1 at the left endpoint.
                if side=='left':
                    B=[1-sigma(p.log(y))]+[-p.diff(lambda yy:sigma(p.log(yy)),y,j) for j in range(1,5)]
                else:B=[p.diff(lambda yy:sigma(1-p.log(yy)),y,j) for j in range(5)]
                V=[z*(4*direct.mpf(p.nstr(v,340))) for v in B];zero=z*0
                physical=source.native.original.physical_mixed(direct,Z,direct.mpf('.001'),E,[zero-direct.mpf('.5')]+[zero]*3,V,hist,P0,1/P**2)
                for name in ('physical_velocity_pressure_y_Z_mixed4','physical_five_primitive_y_Z_mixed4'):
                    for key,grid in physical[name].items():
                        for order,value in grid.items():
                            for x in ep(value):
                                try:compare.enclosed(p,c,got['values']['physical'][name][key][order],p.mpf(x))
                                except AssertionError as exc:raise AssertionError((side,kind,a,b,p.nstr(r,12),name,key,order)) from exc
                                count+=1
                U=physical['physical_velocity_pressure_y_derivative_Taylor']['Utheta_over_Pstar']
                Hrows=physical['actual_normalized_primitive_y_derivative_axial5']
                press=[P0+Hrows['p'][0]]+Hrows['p'][1:]
                parts=stress.raw_pre_stress_rows(direct,direct.mpf('.001'),z,U,V,Hrows,press)
                radius=direct.exp(10*direct.mpf('.1')+direct.ln(110)+10*logP+yy)
                for label,omit in (('theta','variable_radial_shear'),('axial','axial_radial_shear')):
                    rawnum=zero
                    for key,part in parts[label].items():
                        if key==omit:continue
                        factor=1 if label=='theta' else 1/P if part['mode'][1]==0 else P
                        rawnum=rawnum+part['shape'][0]*factor
                    ref=rawnum*(radius)/U[0]
                    for n in (0,1):
                        for x in ep(ref[n]):compare.enclosed(p,c,roots['p1' if label=='theta' else 'p2'][(0,n)],p.mpf(x));root_checks+=1
    # Analytic mass tests exercise nonzero microscopic corrections separately
    # from the baseline 1 in exp(x); no repeated accepted q derivative suite.
    mass_checks=0
    for d in ('1e-40','.001','.2'):
        distance=coordinates.scalar(c.mpf(d));sig=coordinates.scalar(c.mpf((0,c.mpf(d)**2)))
        for rate in (0,1):
            value=source.positive_prefix_integral(sig,distance,rate)
            ref=p.quad(lambda x:x*x*p.exp(-rate*x),[0,p.mpf(d)])
            compare.enclosed(p,c,value,ref);mass_checks+=1
    return dict(passed=True,independent_original_physical_and_primitive_mixed4_endpoint_comparisons=count,
        independent_original_signed_stress_Pstar_and_radius_comparisons=root_checks,
        independent_physical_width_checks=width_checks,independent_weighted_prefix_mass_checks=mass_checks,
        finer_original_slope_mass_enclosures_used_only_in_finite_fixtures=True,
        original_unmodified_ordinary_physical_and_stress_operators_used=True,
        finite_M1_fixture_not_substituted_for_actual_M40_source=True)


def tile_limits_and_memory(c,data):
    assert data['complete_original_axial_operator_on_two_strict_sign_Z_tiles_executed']
    assert not data['full_buffer_global_N_targets_controls_recursion_admitted']
    rows=data['actual_original_endpoint_source_queries'];assert len(rows)==12
    for row in rows:
        g=row['original_typed_native_source'];local=row['actual_original_density_and_local_integrals']
        assert g['ordinary_y_derivatives_and_all_five_history_ODEs_unmodified']
        assert local['status']=='enclosed' and local['original_nonlinear_density_precedes_overlapping_branch_hull']
        assert row['actual_original_phase']['actual_phase_is_N_times_original_log_radius']
        assert ep(iv(c,local['original_true_ordinary_y_width']['coefficient_interval']))[0]>0
        assert all(x['actual_original_inverse_first_jet']['status']=='enclosed' for x in local['original_conditional_phase_inverse_and_density_queries'])
    transport=data['original_whole_axial_transport_with_genuine_incoming']
    steps=transport['original_increasing_y_source_steps'];labels=[(s['endpoint'],s['label']) for s in steps]
    expected=[('left',r[0]) for r in source.axial.PLAN]+[('middle','strict_flat_middle')]+[('right',r[0]) for r in reversed(source.axial.PLAN)]
    assert labels==expected and len(steps)==15
    assert transport['original_exact_endpoint_telescope']==source.original_axial_endpoint_telescope()
    incoming=transport['actual_genuine_slope_incoming']['original_genuine_slope_exit_record']
    assert incoming['ordered_source_cells']==2048 and incoming['candidate_N']==source.N and incoming['exact_y_window']==['0','1']
    assert transport['P0_separate_and_added_once'] and transport['entire_original_axial_phase_zero_to_one_covered']
    assert not transport['original_full_buffer_Rc_targets_controls_global_N_axis_heat_stress_recursion_admitted']
    assert ep(iv(c,transport['original_whole_axial_C0_Z_operator']['incoming_C0_Z_decay_coefficients']['p']['coefficient_interval']))==(1,1)
    for step in steps:
        if step['zero_own_source_keeps_incoming_and_pressure_memory']:
            assert all(v['coefficient_interval']['lower']=='0.0' and v['coefficient_interval']['upper']=='0.0' for v in step['local_signed_C0_increment'].values())


def run():
    began=time.monotonic();manifest=json.loads((HERE/source.NAME).read_bytes())
    assert manifest[source.GATE] and manifest['candidate_N']==source.N
    assert manifest['complete_original_axial_C0_Z_operator_with_genuine_slope_incoming_executed']
    assert not manifest['full_buffer_global_N_targets_controls_recursion_admitted']
    for name,digest in manifest['input_hashes'].items():assert sha(name)==digest,name
    assert all(manifest[k] is False for k in source.common.current.FLAGS)
    c=MPIntervalContext();c.dps=240;queries=0;rows=0
    tail=json.loads((HERE/source.complete.tail.NAME).read_bytes())
    with mp.workdps(300):
        independent=independent_original_endpoint_checks(manifest['source_family'])
        for archive in manifest['actual_original_axial_endpoint_native_source_archives']:
            compressed=(HERE/archive['filename']).read_bytes();raw=gzip.decompress(compressed)
            assert len(compressed)==archive['compressed_bytes'] and len(raw)==archive['uncompressed_bytes']
            assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
            data=json.loads(raw);parent=data['accepted_original_native_parent_archive']
            original=gzip.decompress((HERE/parent['filename']).read_bytes())
            assert hashlib.sha256(original).hexdigest()==parent['lossless_original_json_sha256']
            coordinates=source.native.HalfPstarCoordinates(c,iv(c,data['common_directed_coordinate_theorem']['common_log_bases'][1]),manifest['source_family'])
            coord=source.axial.AxialCoordinates(coordinates,iv(c,data['original_parameter_sources']['original_eta_log']))
            replay=source.execute_tile(coord,json.loads(original),tail)
            for key,value in replay.items():exact_replay_equal(json.loads(json.dumps(source.encode(value))),data[key],key)
            tile_limits_and_memory(c,data);queries+=12;rows+=120
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=sha(source.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{source.GATE:True},source_family=manifest['source_family'],candidate_N=source.N,
        independent_original_endpoint_checks=independent,actual_original_axial_endpoint_queries_replayed=queries,
        actual_original_signed_C0_Z_integral_rows_replayed=rows,whole_axial_operator_cells_replayed=30,
        original_increasing_y_order_exact_telescope_and_quiet_pressure_memory_checked=True,
        genuine_2048_cell_slope_incoming_and_axial_exit_correction_own_history_pressure_replayed=True,
        accepted_upstream_owners_and_producers_not_rerun=True,
        full_buffer_Rc_update_global_Z_N_targets_controls_heat_stress_recursion_not_claimed=True,
        **dict.fromkeys(source.common.current.FLAGS,False),input_hashes=hashes,execution_seconds=time.monotonic()-began)
    (HERE/source.RECEIPT).write_text(json.dumps(source.encode(result),indent=2)+'\n',encoding='utf8')
    print('Whole original axial C0/Z operator and genuine incoming PASS; downstream buffer remains open',flush=True)
    return result


if __name__=='__main__':run()
