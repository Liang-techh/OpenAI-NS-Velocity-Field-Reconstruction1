"""Whole-cell source checks and independent defining-ODE integral reference.

The producer's complete JSON is shipped losslessly compressed. The independent
floating ODE computes only the eta=0 limiting reference at the same native
phase, with an explicit bound on the native positive-eta difference. The
certified result remains the directed actual-source whole-cell integral.
"""
import gzip
import hashlib
import json
import math
from pathlib import Path
import time
import mpmath as mp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
import lei_ren_part1_paper_compliant_current_original_O2_pressure_own_integral as current

base=current.base;ep=current.ep
COMPACT=current.NAME+'.gz'


def interval(c,value):return base.conditioned.packets.interval(c,value)


def independent_limit_reference(origin,N):
    """Defining J/P ODE plus scalar phase inversion; no native eta replacement."""
    p=mp.mp.clone();p.dps=100
    phi0=float(current.density.slow.read_scalar(p,origin['approximate_original_fractional_phase']))
    loop=base.point.source.inertial.profiles.loop
    def rhs(y,values):
        sigma=float(loop.flat_step(mp.fp,y));a=.8+1.2*sigma
        q2=max(0.,(2-a)/(2*a));nu=1+2*q2;phi=(phi0+N*y)%1
        target=lambda x:x+q2*math.sin(4*math.pi*x)/(2*math.pi*nu)-phi
        psi_fraction=brentq(target,0.,1.,xtol=5e-15)
        A=a*(phi-psi_fraction)/2
        E=math.exp(y/10-.6*values[0]);EN=E*math.exp(A/N)
        before=current.density.recovery.history_densities(E,0.)['p']
        after=current.density.recovery.history_densities(EN,0.)['p']
        return [sigma,after-before]
    values=[]
    for tol in (1e-10,5e-12):
        got=solve_ivp(rhs,(0.,1.),(0.,0.),method='DOP853',rtol=tol,atol=tol/100,max_step=1/(64*N))
        assert got.success and abs(got.y[0,-1]-.5)<1e-11
        values.append(float(got.y[1,-1]))
    assert abs(values[0]-values[1])<1e-10
    return dict(passed=True,independent_defining_J_and_pressure_ODE_values=values,
        numerical_refinement_difference=abs(values[0]-values[1]),
        same_original_native_phase_origin_and_candidate_N_used=True,
        eta_zero_is_diagnostic_limit_only_not_selected_native_field=True,
        limit_error_proof='At fixed native phase |dA/d(q²)|<=a/(4*pi*nu)<=1/(2*pi), d(q²)/deta=1/a<=5/4, |dp/dA|<=exp(2.2). Thus |integral_native-integral_eta0|<=2*eta.',
        floating_ODE_not_used_as_certified_native_integral_proof=True)


def native_whole_cells(manifest):
    c=base.MPIntervalContext();c.dps=260;p=mp.mp.clone();p.dps=100
    widths=[];integrals=[];source_checks=0;wraps=[];eta_cap=None
    point_data=json.loads((current.HERE/current.density.NAME).read_bytes())
    anchor=next(row for row in point_data['actual_original_O2_signed_density_queries']
        if row['original_Z_exact']=='0' and row['explicit_candidate_N']==7)
    anchor_pressure=interval(c,anchor['actual_signed_density_point_queries'][0]['five_signed_own_rate_density_enclosures']['p']['coefficient_interval'])
    # Midplane normalized pressure has no physical source exponent here.
    record=anchor['actual_signed_density_point_queries'][0]['five_signed_own_rate_density_enclosures']['p']
    assert not any(record['formal_positive_scale']['source_exponents']) and record['formal_positive_scale']['radius_power']==0
    anchor_pressure*=c.exp(interval(c,record['formal_positive_scale']['additional_log_interval']))
    with mp.workdps(300):
        for refinement in manifest['actual_original_pressure_integral_refinements']:
            count=refinement['ordered_source_cells'];cells=refinement['whole_source_cells']
            assert len(cells)==count and refinement['exact_y_window']==['0','1']
            assert refinement['original_Z_exact']=='0' and refinement['explicit_candidate_N']==7
            assert refinement['original_own_rate']==0 and refinement['no_additional_R_or_Jacobian_applied_to_normalized_density']
            assert refinement['original_inlet_pressure_and_P0_not_reset'] and refinement['contribution_only_requires_adding_incoming_pressure_history']
            assert not refinement['full_five_controls_or_Z_functional_terminal_identity_installed']
            total=c.mpf(0);dy=c.mpf(1)/count
            for i,row in enumerate(cells):
                assert row['exact_y_cell']==[str(i)+'/'+str(count),str(i+1)+'/'+str(count)]
                assert row['entire_cell_source_and_inverse_enclosure_not_point_quadrature']
                assert not row['source_caps_or_midpoints_selected_as_field_values']
                assert not row['q_source_cell']['exact_zero'] and not row['positive_eta_source']['exact_zero']
                eta_log=interval(c,row['positive_eta_source']['log_absolute_upper'])
                assert ep(eta_log)[1]<ep(-2*c.dps*c.ln(10))[0]
                eta_cap=ep(c.exp(-2*c.dps*c.ln(10)))[1]
                hull=interval(c,row['signed_pressure_density_whole_cell']);contribution=interval(c,row['signed_rate0_integral_contribution'])
                actual=dy*hull
                assert ep(actual)==ep(contribution)
                total+=contribution
                for phase in row['true_common_N_phase_boxes']:
                    pl,ph=ep(interval(c,phase));assert 0<=pl<=ph<=1
            saved=interval(c,refinement['directed_signed_pressure_own_rate0_integral'])
            assert ep(total)==ep(saved)
            widths.append(ep(saved)[1]-ep(saved)[0]);integrals.append(saved);wraps.append(refinement['phase_wrap_cell_count'])
            index=int(p.mpf('.53')*count);lo,hi=ep(interval(c,cells[index]['signed_pressure_density_whole_cell']))
            al,ah=ep(anchor_pressure);assert lo<=al<=ah<=hi
            for y in (p.mpf('.23'),p.mpf('.53'),p.mpf('.79'),p.mpf('.97')):
                cell=cells[min(int(y*count),count-1)]
                J=p.quad(lambda x:base.point.source.inertial.profiles.loop.flat_step(p,x),[0,y/2,y])
                jl,jh=ep(interval(c,cell['J_source_cell']));assert jl<=J<=jh
                f=p.exp(y/10-p.mpf(3)*J/5);fl,fh=ep(interval(c,cell['f_source_cell']));assert fl<=f<=fh
                source_checks+=2
        assert all(widths[i+1]<widths[i]/2 for i in range(len(widths)-1))
        lower=max(ep(v)[0] for v in integrals);upper=min(ep(v)[1] for v in integrals);assert lower<=upper
        finest=integrals[-1];assert ep(finest)[1]<0
    return dict(passed=True,actual_original_whole_cell_integral_levels=len(integrals),
        ordered_source_cell_counts=[r['ordered_source_cells'] for r in manifest['actual_original_pressure_integral_refinements']],
        integral_interval_widths=widths,finest_directed_signed_integral=finest,
        finest_integral_strictly_negative=True,whole_refinement_intervals_have_common_intersection=True,
        independent_defining_J_f_source_checks=source_checks,accepted_actual_source_point_pressure_anchor_enclosed=True,
        exact_phase_wrap_cell_counts=wraps,native_eta_upper_for_independent_limit_error=eta_cap,
        native_eta_to_zero_limit_pressure_integral_difference_upper=ep(c.mpf(eta_cap)*2)[1],
        incoming_pressure_memory_and_original_P0_not_reset=True,
        partial_pressure_contribution_not_claimed_as_full_five_control_closure=True)


def run():
    began=time.monotonic();raw=current.HERE/current.NAME;packed=current.HERE/COMPACT
    content=raw.read_bytes() if raw.exists() else gzip.decompress(packed.read_bytes())
    manifest=json.loads(content);assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    packed.write_bytes(gzip.compress(content,compresslevel=9,mtime=0))
    assert gzip.decompress(packed.read_bytes())==content
    flags=('actual_changed_five_moment_integral_evaluated','numerical_original_source_point_or_integral_oracle_installed',
        'actual_five_controls_installed','current_whole_N_selected',*base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    contracts=native_whole_cells(manifest)
    ref=independent_limit_reference(manifest['actual_original_pressure_integral_refinements'][-1]['original_true_radius_phase_at_y0'],7)
    lo,hi=ep(contracts['finest_directed_signed_integral'])
    assert all(lo<value<hi for value in ref['independent_defining_J_and_pressure_ODE_values'])
    result=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],
        actual_original_source_whole_cell_integral_contracts=contracts,
        independent_defining_ODE_pressure_limit_reference=ref,
        one_actual_original_O2_pressure_rate0_integral_enclosed=True,
        compressed_producer_report=dict(filename=COMPACT,lossless_original_json_sha256=hashlib.sha256(content).hexdigest(),
            uncompressed_bytes=len(content),compressed_bytes=packed.stat().st_size),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],COMPACT:current.sha(COMPACT),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began,
        scope='One actual original full-window O2 pressure own-rate0 contribution at Z0/N7, from whole-cell source/phase/inverse enclosures. Independent ODE and positive-eta limit bound corroborate it; full five histories/controls, Z-functional terminal matching, global N and recursion remain open.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual original O2 pressure own integral: source cells, refinement, independent ODE and signed contribution PASS',flush=True)
    return result


if __name__=='__main__':run()
