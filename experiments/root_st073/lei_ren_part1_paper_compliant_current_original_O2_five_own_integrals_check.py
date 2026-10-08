"""Five directed own integrals, defining Duhamel ODE and inlet-memory checks.

The independent floating calculation is a diagnostic eta=0 reference at the
fixed native phase, with separate O(sqrt(eta))/O(eta) comparison bounds. It
does not replace the directed positive-eta source cell integral certificate.
No ancestor producer or original nested mass quadrature is reexecuted.
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
import lei_ren_part1_paper_compliant_current_original_O2_five_own_integrals as current

base=current.base;ep=current.ep;KEYS=tuple(current.RATES)
FLAGS=('C1_or_Z_functional_transport_installed','actual_changed_five_moment_integral_evaluated',
    'numerical_original_source_point_or_integral_oracle_installed','actual_five_controls_installed',
    'current_whole_N_selected',*base.point.source.inertial.profiles.loop.OPEN)


def interval(c,value):return current.interval(c,value)


def overlaps(a,b):
    return max(ep(a)[0],ep(b)[0])<=min(ep(a)[1],ep(b)[1])


def original_mass(c,rate,left,right):
    """Independent closed-form final-endpoint positive exponential integral."""
    r=c.mpf(rate)
    return right-left if r==0 else (c.exp(-r*(1-right))-c.exp(-r*(1-left)))/r


def independent_references(origin,N):
    p=mp.mp.clone();p.dps=100
    phi0=float(current.density.slow.read_scalar(p,origin['approximate_original_fractional_phase']))
    rates=[float(current.RATES[key]) for key in KEYS]
    loop=base.point.source.inertial.profiles.loop
    def rhs(y,values):
        sigma=float(loop.flat_step(mp.fp,y));a=.8+1.2*sigma
        q2=max(0.,(2-a)/(2*a));nu=1+2*q2;phi=(phi0+N*y)%1
        target=lambda x:x+q2*math.sin(4*math.pi*x)/(2*math.pi*nu)-phi
        fraction=brentq(target,0.,1.,xtol=5e-15);angle=2*math.pi*fraction
        A=a*(phi-fraction)/2;E=math.exp(y/10-.6*values[0])
        B=-a*E*math.sqrt(q2)*math.sin(angle)/(2*math.pi)
        EN=E*math.exp(A/N);VN=B/N
        before=current.density.recovery.history_densities(E,0.)
        after=current.density.recovery.history_densities(EN,VN)
        return [sigma]+[after[k]-before[k]-rates[i]*values[1+i] for i,k in enumerate(KEYS)]+[
            before[k]-rates[i]*values[6+i] for i,k in enumerate(KEYS)]
    results=[]
    for tol in (1e-10,5e-12):
        got=solve_ivp(rhs,(0.,1.),[0.,*([0.]*5),0.,5/8,0.,-5/12,5/2],
            method='DOP853',rtol=tol,atol=tol/100,max_step=1/(64*N),dense_output=True)
        assert got.success and abs(got.y[0,-1]-.5)<1e-11
        values=got.sol([.25,.5,.75,1.])
        results.append(dict(tolerance=tol,terminal_window_contributions={k:float(values[1+i,-1]) for i,k in enumerate(KEYS)},
            terminal_original_histories={k:float(values[6+i,-1]) for i,k in enumerate(KEYS)},
            prefixes=[dict(y=float(y),window_contributions={k:float(values[1+i,j]) for i,k in enumerate(KEYS)},
                original_histories={k:float(values[6+i,j]) for i,k in enumerate(KEYS)})
                for j,y in enumerate((.25,.5,.75,1.))]))
    differences={key:max(abs(results[0][field][key]-results[1][field][key])
        for field in ('terminal_window_contributions','terminal_original_histories')) for key in KEYS}
    assert max(differences.values())<1e-9
    return dict(passed=True,independent_defining_J_and_five_Duhamel_ODE_runs=results,
        terminal_numerical_refinement_differences=differences,
        original_nonzero_inlets=[0.,5/8,0.,-5/12,5/2],
        original_history_densities_evaluated_before_and_after_modulation=True,
        original_B_and_axial_increment_retained=True,
        same_original_native_phase_origin_and_candidate_N_used=True,
        eta_zero_is_diagnostic_limit_only_not_selected_native_field=True,
        floating_ODE_not_used_as_certified_native_integral_proof=True)


def native_limit_error(c,eta):
    """Fixed phi/N7 error; constants refer to differences, not magnitudes."""
    root=ep(c.sqrt(c.mpf(eta)))[1]
    errors={key:ep(c.mpf(eta if key in ('h','p') else root)*factor)[1]
        for key,factor in dict(m=1,h=1,k=4,e=6,p=2).items()}
    return dict(passed=True,native_positive_eta_upper=eta,
        normalized_five_contribution_native_minus_eta0_absolute_error_upper=errors,
        fixed_native_phase_comparison=True,candidate_N=7,
        inverse_derivative_proof='For s=q^2, |dpsi/ds|=|sin(2psi)|/[nu*(nu+2*s*cos(2psi))] <=1/[nu*sqrt(1+4*s)]<=1. s_eta-s_0=eta/a<=5*eta/4.',
        primitive_error_proof='|Delta A|<=5*eta/(8*pi); |Delta B|<=exp(.1)/pi*(sqrt(5*eta/4)+sqrt(.75)*5*eta/4)<sqrt(eta) for eta<=.5.',
        density_error_proof='At N7, |u|<.5, |Delta u|<=sqrt(eta), |Delta delta_E|<=eta, |E_N|<=exp(1.1). Five density differences are bounded by sqrt(eta), eta, 4sqrt(eta), 6sqrt(eta), 2eta.',
        error_not_absolute_density_bound=True,
        nonnegative_rates_and_unit_window_preserve_sup_density_error=True,
        physical_Pstar_R_normalizations_not_materialized=True)


def native_cells(manifest,parent):
    c=base.MPIntervalContext();c.dps=260;p=mp.mp.clone();p.dps=c.dps+60
    all_integrals={key:[] for key in KEYS};widths={key:[] for key in KEYS}
    mass_checks=0;pressure_checks=0;prefix_checks=0
    with mp.workdps(300):
        for got,saved in zip(manifest['actual_original_five_own_integral_refinements'],
                             parent['actual_original_pressure_integral_refinements'],strict=True):
            count=got['ordered_source_cells'];cells=got['whole_source_cell_density_and_mass_records']
            assert count==saved['ordered_source_cells'] and len(cells)==count
            assert got['source_family']==manifest['source_family']==saved['source_family']
            assert got['original_Z_exact']=='0' and got['exact_y_window']==['0','1'] and got['explicit_candidate_N']==7
            assert got['own_rates']==current.RATES and got['normalized_own_units']==current.UNITS
            assert got['original_P0_datum_sha256']==manifest['source_family']['datum_enclosure_sha256']
            assert got['original_P0_not_reset'] and got['exact_positive_Duhamel_mass_not_an_extra_R_Jacobian']
            assert got['original_reference_inlet_binding']['original_nonzero_h_e_p_preserved']
            assert got['C0_midplane_five_transport_only'] and not got['C1_or_Z_functional_transport_installed']
            assert all(not got[k] for k in ('actual_five_controls_installed','current_whole_N_selected',
                'numerical_original_source_point_or_integral_oracle_installed'))
            assert got['accepted_whole_source_cache']==dict(filename=current.pressure.NAME+'.gz',source_cell_level=count)
            total={key:c.mpf(0) for key in KEYS};old_total={key:c.mpf(0) for key in KEYS}
            prefix={key:c.mpf(0) for key in KEYS};old_prefix=current.original_inlet(c)
            width=c.mpf(1)/count;local,unused,decay=current.masses(c,width,c.mpf(0),width)
            for i,(row,source) in enumerate(zip(cells,saved['whole_source_cells'],strict=True)):
                assert row['original_source_cache_cell_index']==i and row['exact_y_cell']==source['exact_y_cell']
                assert row['B_and_delta_V_not_zeroed_with_original_V'] and row['source_ranges_not_selected_as_points']
                assert row['true_common_N_phase_piece_count']==len(source['true_common_N_phase_boxes'])
                left,right=c.mpf(i)/count,c.mpf(i+1)/count
                lo,hi=p.mpf(i)/count,p.mpf(i+1)/count
                f=interval(c,source['f_source_cell']);old_rhs=dict(m=c.mpf(0),h=f,k=c.mpf(0),e=-f*f/2,p=f*f/2)
                for key in KEYS:
                    hull=interval(c,row['five_signed_density_whole_cell_hulls'][key])
                    weight=interval(c,row['positive_original_final_endpoint_masses'][key])
                    contribution=interval(c,row['five_signed_final_endpoint_contributions'][key])
                    wl,wh=ep(weight);assert 0<wl<=wh
                    exact=original_mass(p,current.RATES[key],lo,hi);assert wl<=exact<=wh
                    assert ep(hull*weight)==ep(contribution)
                    total[key]+=contribution;old_total[key]+=old_rhs[key]*weight
                    prefix[key]=prefix[key]*decay[key]+hull*local[key]
                    old_prefix[key]=old_prefix[key]*decay[key]+old_rhs[key]*local[key]
                    mass_checks+=1
                assert ep(interval(c,row['five_signed_density_whole_cell_hulls']['p']))==ep(interval(c,source['signed_pressure_density_whole_cell']))
                pressure_checks+=1
                if (i+1)%(count//4)==0:
                    stored=got['actual_cumulative_window_contribution_prefixes'][(i+1)//(count//4)-1]
                    assert stored['exact_y']==str(i+1)+'/'+str(count) and stored['physical_incoming_defects_not_set_to_zero']
                    for key in KEYS:
                        assert ep(prefix[key])==ep(interval(c,stored['window_contributions'][key]))
                        assert overlaps(old_prefix[key],interval(c,stored['original_own_histories'][key]))
                        prefix_checks+=1
            for key in KEYS:
                stored=interval(c,got['five_own_rate_integral_contributions'][key])
                assert ep(total[key])==ep(stored)
                decay_total=interval(c,got['incoming_defect_decay_coefficients'][key])
                assert ep(decay_total)[0]<=p.exp(-p.mpf(current.RATES[key]))<=ep(decay_total)[1]
                old=current.original_inlet(c)[key]*decay_total+old_total[key]
                assert overlaps(old,interval(c,got['original_five_histories_at_y1'][key]))
                assert overlaps(prefix[key],stored)
                all_integrals[key].append(stored);widths[key].append(ep(stored)[1]-ep(stored)[0])
            assert ep(interval(c,got['five_own_rate_integral_contributions']['p']))==ep(interval(c,saved['directed_signed_pressure_own_rate0_integral']))
        for key in KEYS:
            assert all(widths[key][i+1]<widths[key][i]/2 for i in range(len(widths[key])-1))
            assert max(ep(v)[0] for v in all_integrals[key])<=min(ep(v)[1] for v in all_integrals[key])
        finest={key:values[-1] for key,values in all_integrals.items()}
        assert ep(finest['m'])[1]<0 and ep(finest['h'])[1]<0 and ep(finest['e'])[0]>0 and ep(finest['p'])[1]<0
        # k is deliberately not assigned a sign: the finest enclosure crosses zero.
        assert ep(finest['k'])[0]<0<ep(finest['k'])[1]
        eta=ep(c.exp(-2*c.dps*c.ln(10)))[1]
        assert all(ep(interval(c,row['positive_eta_source']['log_absolute_upper']))[1]<ep(-2*c.dps*c.ln(10))[0]
            for row in parent['actual_original_pressure_integral_refinements'][-1]['whole_source_cells'])
    return dict(passed=True,ordered_source_cell_counts=[r['ordered_source_cells'] for r in manifest['actual_original_five_own_integral_refinements']],
        independent_closed_form_positive_Duhamel_cell_mass_checks=mass_checks,
        accepted_pressure_density_cell_equalities=pressure_checks,cumulative_prefix_row_checks=prefix_checks,
        refinement_interval_widths=widths,finest_five_normalized_contribution_intervals=finest,
        all_five_refinement_intervals_have_common_intersection=True,
        finest_proved_signs=dict(m='negative',h='negative',k='unresolved',e='positive',p='negative'),
        incoming_decay_factors_and_original_nonzero_inlets_preserved=True,
        accepted_pressure_integral_exactly_preserved=True),native_limit_error(c,eta)


def inlet_contract(record):
    c=base.MPIntervalContext();c.dps=260;family=record['source_family'];datum=family['datum_enclosure_sha256']
    incoming={k:c.mpf([str(i+1),str(i+1)+'.25']) for i,k in enumerate(KEYS)}
    with mp.workdps(300):
        got=current.apply_incoming(c,record,incoming,source_family=family,original_P0_datum_sha256=datum)
        for key in KEYS:
            base_value=interval(c,record['original_five_histories_at_y1'][key])+interval(c,record['five_own_rate_integral_contributions'][key])
            factor=interval(c,record['incoming_defect_decay_coefficients'][key])
            assert ep(got['own_five_histories_at_y1'][key])==ep(base_value+factor*incoming[key])
        assert ep(interval(c,record['incoming_defect_decay_coefficients']['p']))==ep(c.mpf(1))
    guards=0
    for inlet,fam,P0 in (({'p':1},family,datum),(incoming,dict(family,implicit_source_sha256='wrong'),datum),
                         (incoming,family,'wrong'),(dict(incoming,m=c.mpf('inf')),family,datum)):
        try:current.apply_incoming(c,record,inlet,source_family=fam,original_P0_datum_sha256=P0)
        except ValueError:guards+=1
        else:raise AssertionError('Invalid incoming/source/P0 was admitted')
    assert guards==4
    return dict(passed=True,five_explicit_nonzero_diagnostic_inlet_enclosures_tested=True,
        pressure_rate0_preserves_full_incoming_memory=True,missing_inlet_family_datum_and_nonfinite_guards=guards,
        diagnostic_inlets_not_claimed_as_actual_native_preceding_chart_defects=True,
        original_P0_not_added_to_cumulative_pressure=True,C1_or_Z_functional_transport_installed=False)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;content=gzip.decompress(path.read_bytes());manifest=json.loads(content)
    assert manifest[current.GATE] and manifest['all_five_actual_original_midplane_own_rate_contributions_enclosed']
    assert manifest['no_original_defining_quadratures_or_pressure_cell_producer_reexecuted']
    assert all(manifest[k] is False for k in FLAGS)
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    receipt=json.loads((current.HERE/current.pressure.RECEIPT).read_bytes())
    archive=receipt['compressed_producer_report'];packed=current.HERE/archive['filename']
    raw=gzip.decompress(packed.read_bytes());assert hashlib.sha256(raw).hexdigest()==archive['lossless_original_json_sha256']
    parent=json.loads(raw)
    contracts,eta=native_cells(manifest,parent)
    reference=independent_references(parent['actual_original_pressure_integral_refinements'][-1]['original_true_radius_phase_at_y0'],7)
    c=base.MPIntervalContext();c.dps=260
    with mp.workdps(300):
        for run in reference['independent_defining_J_and_five_Duhamel_ODE_runs']:
            for result in manifest['actual_original_five_own_integral_refinements']:
                for key in KEYS:
                    lo,hi=ep(interval(c,result['five_own_rate_integral_contributions'][key]))
                    error=eta['normalized_five_contribution_native_minus_eta0_absolute_error_upper'][key]
                    value=c.mpf(run['terminal_window_contributions'][key]);vl,vh=ep(value+c.mpf((-error,error)))
                    assert lo<=vl<=vh<=hi,(key,result['ordered_source_cells'])
                    old=run['terminal_original_histories'][key];ol,oh=ep(interval(c,result['original_five_histories_at_y1'][key]))
                    assert ol<=old<=oh
                for stored,ref in zip(result['actual_cumulative_window_contribution_prefixes'],run['prefixes'],strict=True):
                    for key in KEYS:
                        lo,hi=ep(interval(c,stored['window_contributions'][key]));assert lo<=ref['window_contributions'][key]<=hi
                        lo,hi=ep(interval(c,stored['original_own_histories'][key]));assert lo<=ref['original_histories'][key]<=hi
    inlet=inlet_contract(manifest['actual_original_five_own_integral_refinements'][-1])
    result=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],
        actual_original_whole_cell_five_integral_contracts=contracts,
        independent_defining_J_and_five_Duhamel_reference=reference,fixed_native_phase_eta_limit_error=eta,
        explicit_incoming_five_history_affine_contract=inlet,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(content).hexdigest(),
            uncompressed_bytes=len(content),compressed_bytes=path.stat().st_size),
        all_five_actual_original_midplane_own_rate_contributions_enclosed=True,
        **dict.fromkeys(FLAGS,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began,
        scope='Five actual original O2 C0 midplane own-rate contributions at candidate N7 and prescribed original nonzero histories, directed whole-source-cell covers with exact exponential masses and explicit incoming memory. Independent defining ODE and native eta error corroborate. No Z-functional matching, all-chart oracle, controls, global N, C1 stress or coefficient recursion.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual original O2 five own integrals, defining Duhamel ODE and incoming memory PASS',flush=True)
    return result


if __name__=='__main__':run()
