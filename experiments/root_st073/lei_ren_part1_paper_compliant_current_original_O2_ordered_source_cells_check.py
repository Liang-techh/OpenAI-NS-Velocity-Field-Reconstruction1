"""Ordered source certificates and independent defining mass ODE checks.

Independent point checks corroborate continuous source covers. Their floating
reference has an explicit diagnostic tolerance and is not a functional proof.
Exact coefficient/pressure/basis identities and directed whole-cell ranges
provide the production certificate. No ancestor producer is reexecuted.
"""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from scipy.integrate import solve_ivp
import lei_ren_part1_paper_compliant_current_original_O2_ordered_source_cells as current

base=current.base;ep=current.ep
FLAGS=('numerical_original_source_point_or_integral_oracle_installed','inverse_or_modulated_nonmidplane_integral_installed',
    'actual_changed_five_moment_integral_evaluated','actual_five_controls_installed','current_whole_N_selected',
    *base.point.source.inertial.profiles.loop.OPEN)


def interval(c,value):return current.interval(c,value)


def overlaps(a,b):return max(ep(a)[0],ep(b)[0])<=min(ep(a)[1],ep(b)[1])


def defining_ODE_reference():
    loop=base.point.source.inertial.profiles.loop;solutions=[]
    def rhs(y,v):
        sigma=float(loop.flat_step(mp.fp,y));J=v[0]
        return [sigma,*[float(mp.fp.exp(float(rate)*y-.6*power*J)) for rate,power in current.MASS_SPECS]]
    for tol in (1e-11,1e-13):
        got=solve_ivp(rhs,(0.,1.),[0.]*4,method='DOP853',rtol=tol,atol=tol/100,max_step=1/64,dense_output=True)
        assert got.success and abs(got.y[0,-1]-.5)<1e-12
        solutions.append(got)
    difference=max(abs(solutions[0].y[i,-1]-solutions[1].y[i,-1]) for i in range(4))
    assert difference<1e-11
    return solutions,dict(passed=True,independent_defining_J_M0_M1_M2_terminal_values=[list(map(float,v.y[:,-1])) for v in solutions],
        numerical_refinement_difference=float(difference),
        defining_mass_rates=['8/5','1/5','6/5'],defining_mass_powers=[1,2,2],
        floating_solution_not_the_directed_source_function_proof=True)


def ordered_cells(manifest,pressure_manifest,five_manifest,solutions):
    c=base.MPIntervalContext();c.dps=manifest['effective_numerical_precisions']['directed_defining_mass_digits']
    p=mp.mp.clone();p.dps=c.dps+60
    widths=[[] for _ in range(3)];mass_check_count=0;profile_check_count=0;reference_check_count=0
    with mp.workdps(c.dps+40):
        for level,source,own in zip(manifest['actual_original_ordered_source_levels'],
                pressure_manifest['actual_original_pressure_integral_refinements'],
                five_manifest['actual_original_five_own_integral_refinements'],strict=True):
            count=level['ordered_source_cells'];nodes=level['ordered_nodes'];cells=level['whole_source_cells']
            assert count==source['ordered_source_cells']==own['ordered_source_cells'] and len(nodes)==count+1 and len(cells)==count
            assert level['source_family']==manifest['source_family']==own['source_family']
            assert level['original_alpha_binding']==current.source_alpha_binding()
            assert level['original_pressure_cancellation_identity']['passed']
            assert level['no_nested_original_quadrature_or_ancestor_producer_reexecuted']
            assert level['effective_numerical_precisions']==manifest['effective_numerical_precisions']
            total=[c.mpf(0)]*3
            for i,(row,old) in enumerate(zip(cells,source['whole_source_cells'],strict=True)):
                assert row['exact_y_cell']==old['exact_y_cell'] and row['original_source_cache_cell_index']==i
                assert row['source_ranges_not_selected_as_field_points'] and row['nonmonotone_H_D_f_enclosed_on_entire_y_mass_cell']
                J=interval(c,row['original_J_source_cover']);assert ep(J)==ep(interval(c,old['J_source_cell']))
                left,right=c.mpf(i)/count,c.mpf(i+1)/count;ys=current.hull(c,left,right)
                node_left,node_right=nodes[i],nodes[i+1]
                for j,(rate,power) in enumerate(current.MASS_SPECS):
                    r=current.rational(c,rate);weight=interval(c,row['exact_positive_original_mass_weights'][j])
                    expected=(c.exp(r*right)-c.exp(r*left))/r;assert ep(expected)==ep(weight)
                    pr=p.mpf(int(rate.p))/int(rate.q);exact=(p.exp(pr*p.mpf(i+1)/count)-p.exp(pr*p.mpf(i)/count))/pr
                    assert ep(weight)[0]<=exact<=ep(weight)[1] and ep(weight)[0]>0
                    contribution=interval(c,row['directed_original_mass_cell_contributions'][j])
                    assert ep(weight*c.exp(-3*power*J/5))==ep(contribution)
                    assert ep(total[j])==ep(interval(c,node_left['ordered_original_mass_prefixes'][j]))
                    total[j]+=contribution
                    assert ep(total[j])==ep(interval(c,node_right['ordered_original_mass_prefixes'][j]))
                    mass_cell=current.hull(c,interval(c,node_left['ordered_original_mass_prefixes'][j]),
                        interval(c,node_right['ordered_original_mass_prefixes'][j]))
                    assert ep(mass_cell)==ep(interval(c,row['whole_cell_mass_covers'][j]));mass_check_count+=1
                M=[interval(c,v) for v in row['whole_cell_mass_covers']]
                profiles={k:interval(c,v) for k,v in row['whole_original_radial_profile_covers'].items()}
                assert ep(profiles['f'])==ep(interval(c,old['f_source_cell']))
                assert ep(profiles['a'])==ep(interval(c,old['a_source_cell']))
                assert ep(profiles['H'])==ep((c.mpf(5)/8+M[0])*c.exp(-3*ys/2))
                assert ep(profiles['D'])==ep((c.mpf(5)/12+M[2]/2)*c.exp(-ys))
                assert ep(profiles['P'])==ep(c.mpf(5)/2+M[1]/2)
                expected_W=current.hull(c,interval(c,node_right['ordered_original_mass_suffixes'][1]),
                    interval(c,node_left['ordered_original_mass_suffixes'][1]))/2+c.exp(-c.mpf(2)/5)/2
                assert ep(expected_W)==ep(profiles['remaining_pressure_mass'])
                assert all(ep(v)[0]>0 for v in profiles.values());profile_check_count+=6
            for j in range(3):
                replay=c.mpf(0)
                for i in reversed(range(count)):
                    replay+=interval(c,cells[i]['directed_original_mass_cell_contributions'][j])
                    assert ep(replay)==ep(interval(c,nodes[i]['ordered_original_mass_suffixes'][j]))
                for index in (0,count//4,count//2,3*count//4,count):
                    node=nodes[index]
                    assert overlaps(interval(c,node['ordered_original_mass_prefixes'][j])+interval(c,node['ordered_original_mass_suffixes'][j]),total[j])
                widths[j].append(ep(total[j])[1]-ep(total[j])[0])
                for solution in solutions:
                    assert ep(total[j])[0]<=solution.y[j+1,-1]<=ep(total[j])[1];reference_check_count+=1
            alpha=interval(c,level['original_alpha_directed_source_enclosure'])
            assert ep(alpha)==ep(c.mpf(5)/2+total[1]/2+c.exp(-c.mpf(2)/5)/2)
            for key in ('h','e','p'):
                assert overlaps(interval(c,level['original_normalized_midplane_reference_histories_at_y1'][key]),
                    interval(c,own['original_five_histories_at_y1'][key]))
            for index in (count//4,count//2,3*count//4):
                node=nodes[index];y=index/count
                for solution in solutions:
                    J,M0,M1,M2=map(float,solution.sol(y))
                    references=dict(H=(5/8+M0)*mp.fp.exp(-1.5*y),D=(5/12+M2/2)*mp.fp.exp(-y),P=5/2+M1/2)
                    for key,value in references.items():
                        lo,hi=ep(interval(c,node[key]));assert lo<=value<=hi;reference_check_count+=1
            assert ep(interval(c,nodes[0]['H']))==ep(c.mpf(5)/8) and ep(interval(c,nodes[0]['D']))==ep(c.mpf(5)/12)
            assert ep(interval(c,nodes[0]['P']))==ep(c.mpf(5)/2)
        for values in widths:assert all(values[i+1]<values[i]/2 for i in range(len(values)-1))
    return dict(passed=True,exact_original_mass_cell_contribution_checks=mass_check_count,
        full_y_mass_cell_profile_cover_checks=profile_check_count,independent_defining_ODE_terminal_and_prefix_checks=reference_check_count,
        defining_mass_interval_widths_by_refinement=widths,source_mass_rates_not_confused_with_own_history_decay=True,
        prefix_suffix_and_pressure_alpha_relations_preserved=True,
        original_midplane_H_minus_D_P_agree_with_accepted_own_histories=True)


def coefficient_checks(manifest,solutions):
    c=base.MPIntervalContext();c.dps=260;p=mp.mp.clone();p.dps=100
    # The checked original templates are rebuilt symbolically, not the
    # producer's substituted/coefficient lambdas or selected point fields.
    frame=base.point.source.OriginalO2SourceParameterFrame(50)
    templates=base.point.finite_coefficient_templates(frame)
    functions={key:tuple((powers,sy.lambdify(templates['inputs'],expr,modules=[{'mpf':p.mpf},'mpmath'])) for powers,expr in terms)
        for key,terms in templates['rows'].items()}
    counts=0;late_counts=0;unresolved=[];geometries=[];signs=set();diagnostic_tolerance=2e-12
    alpha=mp.fp.mpf(5/2)+solutions[-1].y[2,-1]/2+mp.fp.exp(-.4)/2
    with mp.workdps(300):
        for query in manifest['whole_cell_C0_and_ordinary_Z_coefficient_queries']:
            assert query['source_family']==manifest['source_family'] and query['whole_source_ranges_not_field_points']
            assert query['ordinary_Z_values_not_Taylor_coefficients']
            assert query['original_P0_datum_sha256']==manifest['source_family']['datum_enclosure_sha256']
            basis=query['original_basis_contract'];assert basis['order']==['logPstar','logdelta','logL','zero','logR']
            assert basis['semantic_consumer_binding']==current.source_basis_binding() and basis['delta_and_R_Z_independent']
            assert ep(interval(c,basis['positive_L_interval']))[0]>0
            assert all(not query[k] for k in ('inverse_or_modulated_nonmidplane_integral_installed',
                'numerical_original_source_point_or_integral_oracle_installed','actual_five_controls_installed','current_whole_N_selected'))
            zl,zh=[p.mpf(sy.Rational(v).p)/int(sy.Rational(v).q) for v in query['exact_Z_range']]
            yl,yh=[p.mpf(sy.Rational(v).p)/int(sy.Rational(v).q) for v in query['exact_y_cell']]
            for key,pair in query['full_original_factored_coefficient_function_ranges'].items():
                for order,rows in enumerate(pair):
                    for row in rows:
                        r,ps,d,L=row['original_R_Pstar_delta_L_powers']
                        assert row['no_field_point_or_coefficient_midpoint_selected']
                        if key in ('p1','p2'):
                            assert r==1 and L==-1-order and ps in ((0,) if key=='p1' else (-1,1))
                            if key=='p2':signs.add(ps)
                        if row['late_pressure_coefficient_error_log_upper'] is not None:late_counts+=1
            for y in (yl,(yl+yh)/2,yh):
                J,M0,M1,M2=map(float,solutions[-1].sol(float(y)))
                f=p.exp(y/10-p.mpf(3)*J/5);H=(p.mpf(5)/8+M0)*p.exp(-3*y/2)
                D=(p.mpf(5)/12+p.mpf(M2)/2)*p.exp(-y);P=p.mpf(5)/2+p.mpf(M1)/2
                for z in (zl,(zl+zh)/2,zh):
                    q=1+z*z;values=(z,f,H,D,P,-alpha/q**2,4*alpha*z/q**3,(4-20*z*z)*alpha/q**4)
                    for key,terms in functions.items():
                        rows={tuple(row['original_R_Pstar_delta_L_powers']):row for row in query['full_original_factored_coefficient_function_ranges'][key[0]][key[1]]}
                        for powers,fn in terms:
                            reference=p.mpf(fn(*values));cover=interval(c,rows[powers]['finite_coefficient_function_cover']) if powers in rows else c.mpf(0)
                            lo,hi=ep(cover)
                            assert float(lo)-diagnostic_tolerance<=float(reference)<=float(hi)+diagnostic_tolerance,(key,powers,query['exact_Z_range'])
                            counts+=1
            geometry=query['conditioned_geometry'];branch=geometry['branch'];geometries.append(branch)
            assert not geometry['flat'] and not geometry['q']['exact_zero']
            if branch=='requires_signed_source_refinement':unresolved.append(dict(y=query['exact_y_cell'],Z=query['exact_Z_range']))
            if query['exact_Z0_parity_and_nonzero_p2_Z_preserved']:
                roots=query['native_source_root_enclosures']
                assert roots['p2']['y0_Z0']['exact_zero'] and roots['V']['y0_Z0']['exact_zero']
                assert not roots['p2']['y0_Z1']['exact_zero'] and not roots['V']['y0_Z1']['exact_zero']
    assert late_counts>0 and signs=={-1,1} and 'signed_Mobius' in geometries and 'small_r_series' in geometries
    return dict(passed=True,independent_original_unsubstituted_coefficient_reference_checks=counts,
        floating_reference_diagnostic_absolute_tolerance=diagnostic_tolerance,
        floating_point_checks_not_whole_function_proof=True,
        late_pressure_coefficient_errors_retained=late_counts,p2_original_Pstar_sectors=[-1,1],
        original_conditioned_consumer_semantic_basis_mapping_checked=True,
        exact_midplane_parity_and_nonzero_p2_Z_checked=True,conditioned_geometry_branches=geometries,
        unresolved_signed_source_boxes_explicitly_retained=unresolved,
        nonmidplane_inverse_or_integral_not_claimed=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;content=gzip.decompress(path.read_bytes());manifest=json.loads(content)
    assert manifest[current.GATE] and all(manifest[k] is False for k in FLAGS)
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    p_receipt=json.loads((current.HERE/current.five.pressure.RECEIPT).read_bytes())
    f_receipt=json.loads((current.HERE/current.five.RECEIPT).read_bytes())
    def load_archive(receipt):
        data=receipt['compressed_producer_report'];raw=gzip.decompress((current.HERE/data['filename']).read_bytes())
        assert hashlib.sha256(raw).hexdigest()==data['lossless_original_json_sha256'];return json.loads(raw)
    parent=load_archive(p_receipt);own=load_archive(f_receipt)
    solutions,reference=defining_ODE_reference();cells=ordered_cells(manifest,parent,own,solutions)
    coefficients=coefficient_checks(manifest,solutions)
    result=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],
        actual_original_ordered_source_integral_contracts=cells,independent_defining_mass_ODE_reference=reference,
        full_C0_ordinary_Z_coefficient_source_contracts=coefficients,effective_numerical_precisions=manifest['effective_numerical_precisions'],
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(content).hexdigest(),
            uncompressed_bytes=len(content),compressed_bytes=path.stat().st_size),
        **dict.fromkeys(FLAGS,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began,
        scope='Original ordered defining mass prefix/suffix source cells and full fixed/nonzero/Z-range C0/ordinary-Z coefficient covers, correlated remaining-pressure mass and complete formal source/late-error sectors. Continuous covers and exact symbolic identities certify source ranges; independent floating original ODE/template checks corroborate. Some source geometry boxes remain unresolved. No nonmidplane inverse/integral, functional matching, all-chart oracle, controls, global N, stress or recursion.')
    (current.HERE/current.RECEIPT).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Original O2 ordered source masses, pressure suffix and C0/Z coefficient covers PASS',flush=True)
    return result


if __name__=='__main__':run()
