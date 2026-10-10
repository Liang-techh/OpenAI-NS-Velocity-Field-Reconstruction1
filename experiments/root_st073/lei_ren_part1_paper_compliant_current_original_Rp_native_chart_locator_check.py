"""Focused current affine inverse, source-bound radius and cover checks."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_native_chart_locator as current
from lei_ren_part1_paper_compliant_current_original_Rp_segmented_radius_check import interpreter
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def symbolic_checks(owner,anchors):
    r=owner.radius;g=owner.graph;mu=s.Symbol('mu',positive=True)
    C,P=s.symbols('logC logP',positive=True);L,T,W=s.symbols('L T W',positive=True)
    bindings={ref.node:s.Symbol('current_'+name,positive=True) for name,ref in r.parameters.items()}
    bindings.update({r.parameters['logC'].node:C,r.parameters['logP'].node:P,r.functions['mu'].node:mu,
        r.functions['Lrel'].node:L,r.functions['Ts'].node:T,r.functions['waiting'].node:W})
    q=interpreter(r,False,bindings)
    expected=10*C+10*P+s.log(110)-60*s.log(mu)+12+s.exp(40)
    assert s.simplify(s.expand_log(q.at(r.logRp)-expected,force=True))==0
    identities=0
    for chart,value in current.ANCHORS.items():
        native=s.Rational(value);row=anchors[chart]['chart_rows'][chart]
        got=q.at(row['exact_native_coordinate_function'])
        assert s.simplify(s.expand(got-native))==0,(chart,'exact affine inverse')
        assert row['exact_native_Jacobian']==r.maps[chart]['native_to_log_radius_jacobian']
        assert row['exact_source_log_radius_origin']==r.logRp
        assert row['exact_source_pulse_term']==r.maps[chart]['pulse_term']
        identities+=1
    # This proves the independent polynomial normalizer removes the same
    # microscopic-origin monomial, rather than numerically underflowing it.
    reader=owner.reader();poly=reader.polynomial(r.logRp)
    forbidden={r.parameters[key].node for key in ('hbB','hbS','sc','T')}
    assert all(not forbidden.intersection(key) for key in poly)
    logC=owner.closed.exact.repair.heat.logC
    assert logC._mpi_[0]==logC._mpi_[1]
    c=owner.ctx;lp=c.exp(40)+11;logmu=c.ln(c.mpf('1/1000'))-4*lp
    independent=10*logC+11*lp+1+c.ln(110)-60*logmu
    actual=reader.at(r.logRp)
    a,b=current.ends(actual);x,y=current.ends(independent)
    assert max(a,x)<=min(b,y)
    return dict(passed=True,exact_current_affine_inverse_identities=identities,
        original_absolute_logRp_identity=True,current_logC_exact_dyadic_retained=True,
        microscopic_origin_cancels_before_any_hbB_sc_numeric_evaluation=True,
        independent_directed_logRp_formula_overlap=True)


def native_and_boundary_checks(owner,anchors):
    c=owner.ctx;g=owner.graph;r=owner.radius
    for chart,value in current.ANCHORS.items():
        view=anchors[chart];row=view['chart_rows'][chart];point=Fraction(value)
        lo,hi=current.ends(row['directed_native_coordinate'])
        exact=current.inverse.box(c,point)
        assert lo<=exact.a and exact.b<=hi,(chart,'anchor enclosed')
        assert chart in view['possible_charts'] and chart in view['contained_charts']
        assert row['entire_coordinate_enclosure_in_chart'] and row['interior_resolved']
        assert not view['source_point_value_installed'] and not view['native_source_interval_callback_available']
    domains=owner.maps
    for chart in current.CHARTS:
        lower,upper=owner.domain(chart)
        assert (lower,upper)==(domains[chart]['lower'],domains[chart]['upper'])
    seams=[]
    for left,right in zip(current.CHARTS,current.CHARTS[1:]):
        endpoint=domains[left]['upper'];assert endpoint is not None
        radius=owner.source_endpoint_log_radius(left,'right');view=owner.locate_function(radius)
        assert left in view['possible_charts'] and right in view['possible_charts'],(left,right)
        assert not view['exact_seam_source_call_performed']
        q=interpreter(r,False,{r.functions['waiting'].node:s.Symbol('W',positive=True),
            **{ref.node:s.Symbol(key,positive=True) for key,ref in r.parameters.items()}})
        other=r._map(right,domains[right]['lower'])['logR']
        assert s.simplify(q.at(radius)-q.at(other))==0,(left,right,'exact seam radius')
        seams.append(dict(left=left,right=right,both_native_candidates_retained=True,
            source_packet_not_selected=True))
    # A box which straddles the flatten/power seam must not become a single
    # selected source chart. This is a numerical cover fixture only.
    logR=r._map('flatten',g.constant(100))['logR']
    reader=owner.reader();base=r._map('flatten',g.zero)['logR']
    fixture=copy.copy(reader);fixture.polys={};fixture.memo={}
    # Place the interval on a fresh atomic exact variable only inside this
    # private fixture, leaving every current source bound unchanged.
    offset=g.symbol('native_locator_boundary_cover_fixture')
    fixture.bindings={**reader.bindings,offset.node:c.mpf([-1,1])}
    crossed=owner._locate(g.add(logR,offset),fixture,None)
    assert {'flatten','outer_power'}<=set(crossed['possible_charts'])
    assert not crossed['source_point_value_installed']
    return dict(passed=True,actual_current_interior_affine_anchors=15,
        actual_current_directed_boundary_locations=seams,boundary_overlap_not_discarded=True,
        private_interval_cover_fixture_not_a_physical_source_value=True)


@source_precision
def run(before,observed_owner,anchors,points,locations):
    began=time.monotonic();candidate=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    for name,digest in candidate['input_hashes'].items():assert current.sha(name)==digest,name
    owner=current.CurrentOriginalRpNativeChartLocator(before,require_checked=False)
    assert type(observed_owner) is type(owner) and observed_owner.before is before
    assert observed_owner.hashes==owner.hashes and all(observed_owner.assert_graph().values())
    assert not any(candidate[key] for key in current.GATES+current.OPEN)
    assert candidate['source_family']==owner.family_record and candidate['actual_current_graph']==owner.assert_graph()
    assert set(anchors)==set(current.ANCHORS)
    graph=candidate['exact_current_native_locator_graph'];assert graph==owner.graph.nodes[:len(graph)]
    for key,view in anchors.items():
        assert current.report(view)==candidate['actual_fifteen_source_affine_anchor_locations'][key]
        assert not any(view[k] for k in current.GATES+current.OPEN)
    assert set(points)==set(locations)==set(current.inverse.CASES)|set(current.inverse.LOG_CASES)
    for key,view in locations.items():
        assert current.report(view)==candidate['actual_nine_physical_inverse_locations'][key]
        fresh=owner.locate_inverse(points[key]);assert current.report(fresh)==current.report(view)
        assert fresh['definitely_before_current_Rp'] and not fresh['possible_charts'],key
        assert not fresh['source_point_value_installed']
    exact=symbolic_checks(owner,anchors);cover=native_and_boundary_checks(owner,anchors)
    rejected=[]
    for label,fn in (
        ('numeric_radius_box_is_not_exact_source_function',lambda:owner.locate_function(owner.ctx.mpf(1))),
        ('unregistered_same_graph_constant',lambda:owner.locate_function(owner.graph.constant(1))),
        ('unregistered_moderate_exponential',lambda:owner.locate_function(owner.graph.unary('exp',owner.graph.constant(1)))),
        ('unknown_source_variable',lambda:owner.locate_function(owner.graph.symbol('unadmitted_locator_input'))),
        ('absolute_radius_exponential_not_log_radius',lambda:owner.locate_function(owner.graph.unary('exp',owner.radius.logRp))),
        ('foreign_log_radius_graph',lambda:owner.locate_function(current.radius_source.FunctionRef(copy.deepcopy(owner.graph),0))),
        ('inverse_from_foreign_owner',lambda:owner.locate_inverse(copy.deepcopy(points['ordinary'])))):
        try:fn()
        except (ValueError,TypeError):rejected.append(label)
        else:raise AssertionError('Unadmitted locator source accepted: '+label)
    changed=dict(points['ordinary']);changed['source_logR_enclosure']=owner.ctx.mpf(0)
    try:owner.locate_inverse(changed)
    except ValueError:rejected.append('mutated_inverse_result')
    else:raise AssertionError('Changed inverse result accepted')
    bad=copy.copy(owner);bad.bindings=dict(owner.bindings)
    bad.bindings[owner.radius.parameters['logC'].node]=owner.ctx.mpf(0)
    try:bad.assert_graph()
    except ValueError:rejected.append('substituted_logC_source_bound')
    else:raise AssertionError('Changed exact source parameter accepted')
    bad=copy.copy(owner);bad.maps={key:dict(value) for key,value in owner.maps.items()}
    bad.maps['flatten']['upper']=owner.graph.constant(101)
    try:bad.assert_graph()
    except ValueError:rejected.append('changed_native_source_domain')
    else:raise AssertionError('Changed native source domain accepted')
    bad=copy.copy(owner);bad.source_exponentials=owner.source_exponentials|{owner.graph.unary('exp',owner.graph.constant(1)).node}
    try:bad.assert_graph()
    except ValueError:rejected.append('expanded_parameter_exponential_whitelist')
    else:raise AssertionError('Expanded parameter exponential whitelist accepted')
    try:owner.reader().at(owner.graph.unary('exp',owner.graph.constant(1)))
    except ValueError:rejected.append('nonparameter_exponential_oracle')
    else:raise AssertionError('Unadmitted moderate exponential oracle accepted')
    assert all(row['no_source_packet_or_velocity_value_called'] for row in candidate['actual_source_call_trace'])
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        original_affine_inverse_and_source_logRp_checks=exact,actual_native_coordinate_cover_checks=cover,
        actual_live_physical_inverse_inputs_classified_before_Rp=9,
        original_N_definition=owner.before.original_N_definition,
        all_overlapping_chart_candidates_retained=True,current_exact_parameter_functions_not_replaced_with_bounds=True,
        source_native_interval_and_point_callbacks_left_open=True,
        rejected_unadmitted_inputs=rejected,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            current.PREFIX+'current_original_Rp_segmented_radius_check.py':current.sha(current.PREFIX+'current_original_Rp_segmented_radius_check.py')},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.report(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_NATIVE_LOCATOR current affine inverse and directed chart cover',flush=True)
    return result
