"""Focused same-current pulse seam and exact boundary/coverage checks."""
import copy
import gzip
import json
import math
from pathlib import Path
import time
from types import SimpleNamespace
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_pulse_mixed_seams as current
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import gp_jets
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def raw_recurrence_transfer():
    """Independent raw recovery from the actual normalized primitive ODEs."""
    z=s.Symbol('Z',real=True);mu=s.Symbol('mu',positive=True)
    u=s.Function('current_inlet_u')(z);bp=s.Rational(1,2)+mu
    B=[s.Function('arbitrary_B_y'+str(k))(z) for k in range(5)]
    initial=[s.Function('incoming_'+key)(z) for key in ('m1','X','m2','energy')]
    normalized={key:[v] for key,v in zip(('Mz','Mtheta','Mtheta_z','Mztheta'),initial)}
    rates={'Mz':s.Rational(1,2)-mu,'Mtheta':1-mu,
        'Mtheta_z':s.Rational(1,2)-2*mu,'Mztheta':-2*mu}
    shift=lambda rows,rate,k:sum(s.binomial(k,j)*rate**(k-j)*rows[j] for j in range(k+1))
    product=lambda a,b,k:sum(s.binomial(k,j)*a[j]*b[k-j] for j in range(k+1))
    E=[u*(-bp)**k for k in range(5)]
    V=[u*shift(B,-bp,k) for k in range(5)]
    EV=[product(E,V,k) for k in range(5)]
    signed=[product(V,V,k)-product(E,E,k)/2 for k in range(5)]
    prefactor={'Mz':u,'Mtheta':s.sqrt(2)*u,'Mtheta_z':s.sqrt(2)*u*u,'Mztheta':u*u}
    count=0
    for k in range(1,5):
        square=product(B,B,k-1)
        density={'Mz':B[k-1],'Mtheta':s.Integer(1) if k==1 else s.Integer(0),
            'Mtheta_z':B[k-1],'Mztheta':square-(s.Rational(1,2) if k==1 else 0)}
        for key in normalized:normalized[key].append(density[key]-rates[key]*normalized[key][-1])
        wanted={'Mz':shift(V,1,k-1),'Mtheta':s.sqrt(2)*shift(E,s.Rational(3,2),k-1),
            'Mtheta_z':s.sqrt(2)*shift(EV,s.Rational(3,2),k-1),'Mztheta':shift(signed,1,k-1)}
        for key in normalized:
            difference=s.expand(prefactor[key]*shift(normalized[key],rates[key],k)-wanted[key])
            for n in range(5-k):assert s.diff(difference,z,n)==0,(key,k,n);count+=1
        # Absolute pressure and Mp use a distinct true theta^2 scale for
        # positive radial orders. P0 is an arbitrary independent Z function.
        p0=s.Function('independent_P0')(z)
        for n in range(5-k):
            assert s.diff(s.expand(product(E,E,k-1)/2-u*u*(-2*bp)**(k-1)/2),z,n)==0
            assert s.diff(p0,z,n).diff(s.Symbol('y'))==0
            count+=2
    return dict(passed=True,independent_raw_primitive_density_mixed_identities=count,
        nonzero_active_B_and_arbitrary_incoming_histories_retained=True,
        absolute_pressure_P0_and_radial_theta_squared_units_separate=True)


def interpretation(owner):
    radius=owner.radius
    field=SimpleNamespace(graph=owner.graph,parameters=radius.parameters,
        contracts=radius.frame.bridge.leading.contracts)
    mu=s.Symbol('actual_positive_mu',positive=True)
    return current.pulse.radius.RadiusInterpreter(field,False,{
        radius.logRp.node:s.Symbol('actual_absolute_logRp',real=True),
        radius.functions['mu'].node:mu,
        radius.functions['logP'].node:s.Symbol('actual_logPstar',real=True)}),mu


def exact_geometry_and_scales(owner,views):
    reader,mu=interpretation(owner);scales=0;native=0
    for name,view in views.items():
        left,right=view['left'],view['right']
        assert s.simplify(reader.at(left['geometry']['logR'])-reader.at(right['geometry']['logR']))==0,name
        JL=reader.at(left['geometry']['native_to_log_radius_jacobian'])
        JR=reader.at(right['geometry']['native_to_log_radius_jacobian'])
        assert s.simplify(JL-(1 if name in ('entrance_main','gap_end') else 1/mu))==0
        assert s.simplify(JR-(1 if name in ('gap_coordinate','gap_end') else 1/mu))==0
        for key,grid in left['log_radius_mixed_rows'].items():
            assert len(grid)==15
            for label,a in grid.items():
                b=right['log_radius_mixed_rows'][key][label]
                assert a.powers==b.powers and a.source_units==b.source_units
                alog=sum(reader.at(v.node) for _,v in a.log_scale_parts)
                blog=sum(reader.at(v.node) for _,v in b.log_scale_parts)
                assert s.simplify(alog-blog)==0,(name,key,label,'exact common scale')
                scales+=1
                k,n=a.derivative
                for side,base,J in ((left,a,JL),(right,b,JR)):
                    row=side['native_coordinate_mixed_rows'][key]['n%d_Z%d'%(k,n)]
                    assert row.coefficients is base.coefficients and row.powers[-1]==k
                    total=sum(reader.at(v.node) for _,v in row.log_scale_parts)
                    expected=sum(reader.at(v.node) for _,v in base.log_scale_parts)+k*s.log(J)
                    assert s.simplify(total-expected)==0,(name,key,label,'exact J power')
                    native+=1
    entrance=views['entrance_main']['left']['geometry'];gap=views['gap_coordinate']['right']['geometry']
    assert s.simplify(reader.at(entrance['exact_native_coordinate_function'])-s.Rational(1,50)/mu)==0
    assert s.simplify(reader.at(gap['exact_native_coordinate_function'])+1/mu)==0
    return dict(passed=True,exact_lazy_radius_seams=5,exact_common_scale_mixed_identities=scales,
        exact_native_Jacobian_mixed_identities=native,
        formal_boundaries_are_function_nodes_not_interval_endpoints=True)


def diagnostic_rows(owner,views):
    ends=current.pulse.radius.post.selected.inlet.endpoints
    count=0;differences={}
    for name,view in views.items():
        differences[name]={}
        for key,grid in view['left']['log_radius_mixed_rows'].items():
            differences[name][key]={}
            for label,a in grid.items():
                b=view['right']['log_radius_mixed_rows'][key][label]
                al,ah=ends(a.coefficients[0]);bl,bh=ends(b.coefficients[0])
                assert max(al,bl)<=min(ah,bh),(name,key,label,'directed source inconsistency')
                delta=a.coefficients[0]-b.coefficients[0];dl,dh=ends(delta)
                assert dl<=0<=dh,(name,key,label)
                differences[name][key][label]=dict(signed_common_scale_difference_enclosure=delta,
                    directed_difference_width_bound=owner.ctx.mpf(dh)-owner.ctx.mpf(dl))
                count+=1
    return count,differences


@source_precision
def run(before=None,observed_owner=None,observed_views=None,observed_strips=None):
    began=time.monotonic();candidate=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    for name,digest in candidate['input_hashes'].items():assert current.sha(name)==digest,name
    owner=current.CurrentOriginalRpPulseMixedSeams(before,require_checked=False)
    assert not owner.acceptance_loaded and not any(candidate[k] for k in current.GATES+current.OPEN)
    assert candidate['source_family']==owner.family_record
    assert candidate['actual_source_graph']==owner.assert_graph()
    assert candidate['arbitrary_current_source_mixed4_theorem']==owner.theorem
    assert candidate['source_bound_canonical_theorem']==owner.canonical
    assert candidate['actual_selected_energy_branch_theorem']==owner.energy_branch
    assert candidate['accepted_uncapped_source_certificate_flags']==owner.uncapped_certificate_flags
    assert candidate['actual_current_positive_selection_terminal_equation_proof']==owner.selection_equation_proof
    assert candidate['current_production_bindings']==owner.production_bindings
    if observed_owner is not None:
        # Reuse only typed observations of this exact unchanged source
        # instance. No scalar rows are reconstructed from archived receipts.
        assert type(observed_owner) is type(owner) and observed_owner.before is before
        assert observed_owner.graph is owner.graph and observed_owner.hashes==owner.hashes
        assert all(observed_owner.assert_graph().values())
        assert set(observed_views)==set(current.SEAMS) and set(observed_strips)==set(current.OVERLAPS)
        views=observed_views;strips=observed_strips
    else:
        views={name:owner.evaluate(name,candidate['fresh_Z']) for name in current.SEAMS}
        strips={name:owner.supplemental(name,candidate['fresh_Z'],(a+b)/2)
            for name,(_,a,b) in current.OVERLAPS.items()}
    for name,view in views.items():
        assert current.pulse.raw.packed(current.report(view))==candidate['actual_five_pulse_mixed_seam_views'][name],name
        assert not any(view[k] for k in current.GATES+current.OPEN)
        for side in ('left','right'):
            assert all(v.graph is owner.graph for grid in view[side]['log_radius_mixed_rows'].values()
                for row in grid.values() for _,v in row.log_scale_parts)
    for name,view in strips.items():
        assert current.pulse.raw.packed(current.mixed.view_report(view))==candidate['actual_supplemental_coverage_views'][name]
        reader,_=interpretation(owner)
        chart,a,b=current.OVERLAPS[name];v=(a+b)/2
        expected=owner.radius._map(chart,owner.graph.constant(v))
        assert s.simplify(reader.at(view['geometry']['logR'])-reader.at(expected['logR'].node))==0
    assert current.pulse.raw.packed(owner.coverage())==candidate['rounded_coverage_ledger']
    independent=raw_recurrence_transfer();exact=exact_geometry_and_scales(owner,views)
    count,differences=diagnostic_rows(owner,views)
    for value in gp_jets(owner.ctx,11):assert current.pulse.radius.post.selected.inlet.endpoints(value)==(0,0)
    for center in (-3,-1):
        for value in owner.pulse.flat.beta(-4-center):assert current.pulse.radius.post.selected.inlet.endpoints(value)==(0,0)
    rejected=[]
    for name in ('missing',):
        try:owner.evaluate(name,'.5')
        except ValueError:rejected.append('unknown_seam')
        else:raise AssertionError('Unknown seam accepted')
    for name,v in (('entrance_main_overlap',0),('gap_coordinate_overlap',13)):
        try:owner.supplemental(name,'.5',v)
        except ValueError:rejected.append(name+'_invalid_domain')
        else:raise AssertionError('Out-of-domain supplemental source accepted')
    bad=copy.copy(owner);bad.pulse=copy.copy(owner.pulse);bad.pulse.main=lambda *a,**k:None
    try:bad.assert_graph()
    except (ValueError,AttributeError):rejected.append('changed_original_source_method')
    else:raise AssertionError('Changed original source method accepted')
    bad=copy.copy(owner);bad.uncapped_certificate_flags=dict(owner.uncapped_certificate_flags)
    bad.uncapped_certificate_flags['exact_uncapped_selected_sources_used']=False
    try:bad.assert_graph()
    except ValueError:rejected.append('missing_uncapped_source_certificate')
    else:raise AssertionError('Missing uncapped selection certificate accepted')
    bad=copy.copy(owner);bad.selection_equation_proof=dict(owner.selection_equation_proof)
    bad.selection_equation_proof['same_current_positive_quadratic_gives_forward_full_future_half']=False
    try:bad.assert_graph()
    except ValueError:rejected.append('changed_current_positive_selection_equation')
    else:raise AssertionError('Changed current selection equation accepted')
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        arbitrary_current_source_theorem_counts=dict(base=len(owner.theorem['boundary_five_primitive_function_identities']),
            primitive_mixed=len(owner.theorem['current_instantiated_five_primitive_mixed4_identities']),
            velocity_pressure_mixed=len(owner.theorem['current_instantiated_velocity_absolute_pressure_mixed4_identities']),
            first_linear_axial5=len(owner.theorem['first_linear_moment_axial5_identities'])),
        independent_raw_density_recurrence_transfer=independent,
        exact_five_seam_geometry_scale_and_native_transfer=exact,
        actual_common_scale_overlap_diagnostics=count,
        same_current_signed_scaled_difference_diagnostics=differences,
        actual_supplemental_source_views_checked=len(strips),
        source_flat_gp11_and_beta_minus4_jets_checked=True,
        exact_energy_branch_identity_precedes_any_interval_comparison=True,
        accepted_uncapped_current_positive_selection_and_terminal_energy_equations_required=True,
        actual_same_source_typed_observations_not_receipt_rows_used=True,
        uniform_global_physical_time_and_point_admission_remain_open=True,
        rejected_sources_and_domains=rejected,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.pulse.raw.packed(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_PULSE_MIXED_SEAMS five current joins, formal boundaries and supplemental coverage',flush=True)
    return result


if __name__=='__main__':run()
