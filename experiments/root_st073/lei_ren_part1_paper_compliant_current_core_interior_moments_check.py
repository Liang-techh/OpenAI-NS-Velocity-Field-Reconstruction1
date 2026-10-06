"""Admit six same-fixed-point interior core moments, axis shapes and Euler jets."""
import copy
import gzip
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_core_interior_moments import (
    CurrentCoreInteriorMoments,HERE,NAME,RECEIPT,GATES,OPEN,SPEC,INDICES,VIEWS,
    sha,pack,encode,endpoints,source_precision,key,stirling,integrated_polynomial,mixed_product_tail,_verify_hashes)
from lei_ren_part1_paper_compliant_current_core_recurrence_source import binding


def finite(value):
    lo,hi=endpoints(value)
    if not mp.isfinite(lo) or not mp.isfinite(hi) or lo>hi:raise ValueError('Finite directed core row required')


def independent_operator_proof():
    # Every possible finite density degree (0..48), not just random examples.
    r=s.Symbol('rho',nonnegative=True);rows_count=0
    class Exact:
        @staticmethod
        def mpf(v):return s.Rational(v)
    c=Exact()
    for n in range(49):
        rows=[[s.Integer(int(m==n))] for m in range(n+1)]
        for i in range(5):
            for w,norm in ((0,1),(1,1),(1,2)):
                produced=integrated_polynomial(c,rows,r,i,0,w,norm)
                expected=s.diff(s.Rational(norm,n+w+1)*r**n,r,i)
                if s.expand(produced-expected)!=0:raise ValueError('Normalized finite integration derivative changed')
                rows_count+=1
    t=s.Symbol('t',nonnegative=True)
    weights={str((i,w,norm)):s.integrate(norm*t**(i+w),(t,0,1))==s.Rational(norm,i+w+1)
        for i in range(5) for w,norm in ((0,1),(1,1),(1,2))}
    # Compare mixed product tail with independently expanded Leibniz products.
    tail_rows=0
    PN={(i,k):s.Symbol('p_'+str(i)+'_'+str(k)) for i,k in INDICES}
    QN={(i,k):s.Symbol('q_'+str(i)+'_'+str(k)) for i,k in INDICES}
    PE={(i,k):s.Symbol('e_'+str(i)+'_'+str(k)) for i,k in INDICES}
    QE={(i,k):s.Symbol('f_'+str(i)+'_'+str(k)) for i,k in INDICES}
    for i,k in INDICES:
        expected=s.expand(sum(math.comb(i,a)*((PN[a,l]+PE[a,l])*(QN[i-a,k-l]+QE[i-a,k-l])
            -PN[a,l]*QN[i-a,k-l]) for a in range(i+1) for l in range(k+1)))
        if s.expand(mixed_product_tail(c,PN,QN,PE,QE,i,k)-expected)!=0:raise ValueError('Mixed density product tail identity changed')
        tail_rows+=1
    # Explicit original scaled log derivative: seed ell=-gradient and F0=exp(-C-Lambda*G).
    binding('compliant_core_coefficient_rebuild','seed','fixed',
        "dict(ell_Z_taylor=[-v for v in gradient.coefficients],S_Z_taylor=S,U0_Z_taylor=list(u.coefficients),P0_Z_taylor=[scale*v for v in pressure['normalized_pressure_coefficients']])")
    binding('compliant_current_core_interior_moments','evaluate','ell',
        "[self.core.Lambda*v for v in packet['fixed']['ell_Z_taylor'][:6]]")
    # Bell/Taylor recurrence proven for arbitrary sixth-order log-source jet.
    x=s.Symbol('x');e=s.symbols('e0:6');log=sum(e[j]*x**(j+1)/(j+1) for j in range(6))
    from lei_ren_part1_paper_compliant_current_core_interior_moments import relative_ratios
    ratio={}
    for multiplier in (1,2):
        rows=relative_ratios(c,e,multiplier)
        # A truncated exponential power series is exact through degree6.
        series=s.series(s.exp(multiplier*log),x,0,7).removeO()
        ratio[str(multiplier)]=[s.expand(rows[k]-series.coeff(x,k)*math.factorial(k))==0 for k in range(7)]
    if not all(weights.values()) or not all(all(v) for v in ratio.values()):raise ValueError('Integral/Bell source operator proof failed')
    return dict(all_density_degree0_to48_integrated_derivative_identities=rows_count,
        normalized_integral_tail_weight_identities=weights,mixed_radial4_axial6_product_tail_identities=tail_rows,
        original_seed_scaled_log_derivative_bound=True,relative_F0_and_F0_squared_axial6_identities=ratio,passed=True)


def check_view(view,c):
    if any(view[k] for k in OPEN) or not all(view[k] for k in (
            'original_pressure_axis_and_increment_retained_separately','actual_nonlinear_fixed_point_and_all_density_product_tails_retained',
            'original_mean_axial6_retained_before_Q_axial5','normalized_axis_uses_no_division_by_rho',
            'physical_axis_tensor_limits_not_admitted','source_enclosures_not_resolved_point_values',
            'finite_polynomial_is_not_complete_solution')):raise ValueError('Interior source/scope lost')
    if view['radial_degree']!=24 or view['axial_depth']!=6:raise ValueError('Original coupled radial depth lost')
    layouts=('finite_integrated_axial_Taylor_coefficient_grids','directed_integrated_axial_Taylor_tail_bounds',
        'ordinary_rho4_axial6_moment_grids')
    expected={key(i,k) for i,k in INDICES};count=0
    for layout in layouts:
        if set(view[layout])!=set(SPEC):raise ValueError('Six full core moment shapes required')
        for name,rows in view[layout].items():
            if set(rows)!=expected:raise ValueError('Full rho4/axial6 grid required')
            for value in rows.values():finite(value);count+=1
            if layout.startswith('directed') and any(endpoints(v)[0]<0 for v in rows.values()):raise ValueError('Negative core tail bound')
    log=view['ordinary_logR4_axial6_moment_grids']
    for name,rows in log.items():
        if set(rows)!={'y'+str(q)+'_Z'+str(k) for q in range(5) for k in range(7)}:raise ValueError('Core Euler4/axial6 grid required')
        for q in range(5):
            for k in range(7):
                expected_value=sum((stirling(q,j)*view['rho']**j*view[layouts[-1]][name][key(j,k)] for j in range(q+1)),c.mpf(0))
                if endpoints(rows['y'+str(q)+'_Z'+str(k)])!=endpoints(expected_value):raise ValueError('Core rho-to-logR conversion changed')
                finite(rows['y'+str(q)+'_Z'+str(k)]);count+=1
    for label in ('original_analytic_axis_pressure_axial6_Taylor_coefficients','original_F0_relative_axial6_ordinary_derivatives',
            'original_F0_squared_relative_axial6_ordinary_derivatives'):
        if len(view[label])!=7:raise ValueError('Original full axial6 source required')
        for value in view[label]:finite(value);count+=1
    for label,required in (('normalized_pressure_increment_rho4_axial6',expected),
            ('normalized_radial_Q_rho4_axial5',{key(i,k) for i in range(5) for k in range(6)})):
        if set(view[label])!=required:raise ValueError('Pressure/Q core derivative inventory lost')
        for value in view[label].values():finite(value);count+=1
    axis=endpoints(view['rho'])==(mp.mpf(0),mp.mpf(0))
    if axis:
        for name,expected_value in (('H',1),('B',mp.mpf('.5')),('C',1)):
            if endpoints(view[layouts[-1]][name][key(0,0)])!=(expected_value,expected_value):raise ValueError('Exact nonsingular core axis shape lost')
        for rows in log.values():
            for q in range(1,5):
                for k in range(7):
                    if endpoints(rows['y'+str(q)+'_Z'+str(k)])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Positive Euler derivative nonzero on core axis')
        for k in range(7):
            if endpoints(view['normalized_pressure_increment_rho4_axial6'][key(0,k)])!=(mp.mpf(0),mp.mpf(0)):raise ValueError('Axis pressure increment nonzero')
    return dict(derivative_and_source_rows_checked=count,full_six_rho4_axial6_shapes=True,axis=axis)


@source_precision
def run(field=None):
    raw=json.loads(gzip.decompress((HERE/NAME).read_bytes()));_verify_hashes(raw)
    field=field if field is not None else CurrentCoreInteriorMoments(require_checked=False);field.assert_graph()
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k!='current_core_interior_views'}:raise ValueError('Same core source manifest differs')
    if any(raw[k] for k in GATES+OPEN) or set(raw['current_core_interior_views'])!=set(VIEWS):raise ValueError('Interior producer scope/view inventory differs')
    field.cache={};field.rows_cache={};field.tail_cache={}
    proof=independent_operator_proof();counts={}
    for name,args in VIEWS.items():
        value=field.evaluate(*args)
        if encode(pack(value))!=raw['current_core_interior_views'][name]:raise ValueError('Current full interior moment replay differs: '+name)
        counts[name]=check_view(value,field.ctx)
        print('Check six current interior core moments: '+name,flush=True)
    rejected=[]
    for label,args in (('negative_rho',('.2','-.001')),('past_core',('.2','4.001')),
            ('outside_Z',(2,1)),('nonfinite_Z',('inf',1)),('nonfinite_rho',('.2','inf'))):
        try:field.evaluate(*args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Invalid current core source domain admitted')
    clone=copy.copy(field);clone.rebuild=object()
    try:clone.assert_graph()
    except ValueError:rejected.append('foreign_fixed_point_rebuild')
    else:raise ValueError('Foreign nonlinear source accepted')
    clone=copy.copy(field);clone.common=copy.copy(field.common);clone.common.acceptance_loaded=False
    try:clone.assert_graph()
    except ValueError:rejected.append('unchecked_common_fixed_point')
    else:raise ValueError('Unchecked fixed-point prerequisite accepted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_core_interior_view_count=len(VIEWS),
        current_core_interior_view_counts=counts,independent_current_core_integral_product_Bell_operator_proof=proof,
        current_core_interior_rows_checked=sum(v['derivative_and_source_rows_checked'] for v in counts.values()),
        current_core_interior_source_theorem=field.proof,invalid_domains_and_foreign_sources_rejected=rejected,
        existing_completed_tensor_inventory=raw['existing_completed_tensor_inventory'],
        current_velocity_pressure_atlas_inventory=raw['current_velocity_pressure_atlas_inventory'],
        scope='Six same-fixed-point normalized interior moment functions with full rho4/axial6/Euler4, original P0 plus dressed increment and Q axial5. Nonsingular moment-axis values and source inlet normalization. Core full tensor, completed core/bridge tensor join, axis tensor/remainder, global/cone/time/energy/points/n-dependent recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current six full interior core moment functions/rho4/axial6/Euler4 PASS; core full tensor remains open',flush=True)
    return result


if __name__=='__main__':run()
