"""Original-scale inverse identities, source admission and physical rows."""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_correlated_physical_source as current
from lei_ren_part1_paper_compliant_current_original_Rp_segmented_radius_check import interpreter
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def canonical(value):return json.loads(json.dumps(current.report(value)))


def independent_root_algebra():
    Q,d=s.symbols('q delta',real=True);a=s.Symbol('abs_Z',positive=True)
    lt=2*Q+s.log(1-a*a);lz=s.log(a)+(1-d)*Q
    factor=s.exp(lt-2*Q)+s.exp(2*lz+2*d*Q-2*Q)
    assert s.simplify(s.expand(factor-1))==0
    A,B=s.symbols('time_term axial_term',positive=True)
    derivative=2-2*d*B/(A+B)
    assert s.simplify(derivative-2*(1-d)-2*d*A/(A+B))==0
    return dict(passed=True,original_implicit_equation_forward_factorization=True,
        original_strict_monotonicity_margin_identity=True,
        domain_assumptions='0<delta<1; 0<=abs_Z<1; tau>0; signed Z recovered separately')


def source_graph_identities(owner,requests,inverses,delivered):
    r=owner.radius;g=owner.graph;identities=0
    bindings={ref.node:s.Symbol('original_'+name,positive=True) for name,ref in r.parameters.items()}
    bindings.update({ref.node:s.Symbol('source_'+name,positive=True) for name,ref in r.functions.items()})
    bindings[owner.before.actual.owner.raw.logU0.node]=s.Symbol('source_inlet_logU0',real=True)
    for name,request in requests.items():
        q0=interpreter(r,False,dict(bindings));f=request.forward_coordinates;refs=inverses[name]['coordinate_functions']
        qf=q0.at(f['loglambda']);bound=dict(bindings)
        bound[refs['log_lambda'].node]=qf
        q=interpreter(r,False,bound)
        assert s.cancel(s.expand(q.at(refs['logR'])-q0.at(f['logR'])))==0,(name,'exact logR inverse')
        assert s.simplify(q.at(refs['Z'])-s.Rational(request.Z.numerator,request.Z.denominator))==0,(name,'signed Z inverse')
        identities+=2
        finite=s.Rational(request.finite_log_time_offset.numerator,request.finite_log_time_offset.denominator)
        z=s.Rational(request.Z.numerator,request.Z.denominator)
        assert s.simplify(q0.at(f['log_r'])-(finite+s.log(2)-s.log(1-z*z))/2)==0,(name,'physical log_r cancellation')
        identities+=1
        location=owner.locate(inverses[name]);native=location['chart_rows'][request.chart]['exact_native_coordinate_function']
        assert s.cancel(s.expand(q.at(native)-s.Rational(request.native.numerator,request.native.denominator)))==0,(name,'native inverse')
        identities+=1
        if name in delivered:
            data=delivered[name];token=owner._deliveries[id(data)][2];reader,_,_=owner._validate_delivery(data)
            exact=owner.radius.geometry(request.chart,request.native)
            geometry=data['actual_source_view']['geometry']
            # Request geometry and inverse-token geometry represent the same
            # exact logR; retain this identity when evaluating physical scales.
            q.bindings[token.function.node]=s.Rational(request.native.numerator,request.native.denominator)
            for key in ('logR','offset','origin','pulse_term','native_to_log_radius_jacobian'):
                assert s.cancel(s.expand(q.at(current.box.pulse.radius.FunctionRef(g,geometry[key]))-
                    q.at(current.box.pulse.radius.FunctionRef(g,exact[key]))))==0,(name,key)
                identities+=1
            assert reader.aliases[token.function.node]==g.constant(request.native).node
    return dict(passed=True,independent_symbolic_current_source_coordinate_identities=identities,
        original_log_time_native_and_axial_correlation_preserved=True,
        request_and_inverse_box_geometry_exact_logR_identity_retained=True)


@source_precision
def run(before,observed_owner,requests,inverses,delivered,mapped):
    began=time.monotonic();candidate=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.CurrentOriginalRpCorrelatedPhysicalSource(before,require_checked=False)
    assert type(observed_owner) is type(owner) and observed_owner.before is before
    assert observed_owner.hashes==owner.hashes and all(observed_owner.assert_graph().values())
    for name,digest in candidate['input_hashes'].items():assert current.sha(name)==digest,name
    assert candidate['source_family']==owner.family_record and not any(candidate[k] for k in current.GATES+current.OPEN)
    assert candidate['actual_current_graph']==owner.assert_graph()
    graph=candidate['exact_current_correlated_physical_graph'];assert graph==owner.graph.nodes[:len(graph)]
    assert set(requests)==set(inverses)==set(current.CASES)|set(current.INVERSE_ONLY)
    assert set(delivered)==set(mapped)==set(current.CASES)
    rows=0;times=0;groups=0;summaries={}
    for name,view in inverses.items():
        assert canonical(view)==candidate['actual_correlated_physical_inverses'][name]
        assert view['strict_numerical_Z_interior_resolved'] and view['exact_root_node_retained_with_separate_equality_proof']
        assert not view['numerical_absolute_radius_lambda_time_or_z_materialized']
        assert not any(view[k] for k in current.GATES+current.OPEN)
        located=observed_owner.locate(view);request=requests[name]
        assert request.chart in located['contained_charts'] and not located['definitely_before_current_Rp']
        lr=view['physical_log_radius_enclosure'];lo,hi=current.ends(lr)
        # Independent finite expression, with no original logC arithmetic.
        c=owner.ctx;z=current.inverse.box(c,request.Z);off=current.inverse.box(c,request.finite_log_time_offset)
        finite=(off+c.ln(2)-c.ln(1-z*z))/2;fl,fh=current.ends(finite)
        assert max(lo,fl)<=min(hi,fh) and -100<lo<=hi<100
        assert current.ends(view['directed_inverse_mapping']['requested_log_tau'])[1]<0
    for name,data in delivered.items():
        assert canonical(data)==candidate['actual_correlated_source_deliveries'][name]
        assert data['correlated_physical_input_source_delivery'] and data['native_enclosure_refined_by_exact_affine_identity']
        assert not data['physical_values_materialized']
        source=data['actual_source_view'];request=requests[name]
        assert source['whole_native_and_axial_boxes_passed_to_original_algorithms']
        assert source['source_input_provenance']['kind']=='live_proved_correlated_original_physical_inverse'
        assert source['source_input_provenance']['exact_affine_native_identity']==data['exact_native_inverse_identity'].node
        point=current.inverse.box(owner.ctx,request.native);native=data['source_native_coordinate_enclosure']
        nl,nh=current.ends(native);pl,ph=current.ends(point)
        assert pl<=nl<=nh<=ph
        physical=mapped[name]
        assert json.loads(json.dumps(current.locator.report(current.physical.report(physical))))==candidate['actual_correlated_physical_operator_rows'][name]
        assert physical['actual_correlated_physical_point_source_consumed']
        for component,values in physical['Cartesian_spatial_rows'].items():
            assert len(values)==35
            for label,row in values.items():
                rows+=1
                assert row.component==component and row.terms
                for group in row.groups():
                    groups+=1;assert group['signed_coefficient'].ctx is owner.ctx and group['signed_coefficient'].order==0
                    assert all(ref.graph is owner.graph for _,ref in group['log_scale_parts'])
        assert set(physical['fixed_x_time_rows'])==set(current.physical.COMPONENTS)
        times+=len(physical['fixed_x_time_rows'])
        summaries[name]=dict(chart=request.chart,physical_log_r=inverses[name]['physical_log_radius_enclosure'],
            directed_native_box=native,directed_axial_box=source['directed_axial_coordinate'],
            actual_spatial_rows=140,actual_fixed_x_time_rows=4,actual_source_called=True,
            absolute_physical_values_materialized=False)
    independent=independent_root_algebra()
    exact=source_graph_identities(observed_owner,requests,inverses,delivered)
    # Pure original numerical roots on finite logs corroborate the equality
    # proof. These fixtures do not define original scales or current delta.
    numeric=[]
    for Z in ('521/1000','-521/1000','0'):
        c=owner.ctx;z=current.inverse.box(c,Fraction(Z));q=c.mpf(-7);lt=2*q+c.ln(1-z*z)
        sign=1 if Fraction(Z)>0 else -1 if Fraction(Z)<0 else 0
        lz=None if not sign else c.ln(abs(z))+(1-owner.inverse.delta)*q
        solved=current.inverse.log_coordinate_map(c,lz,lt,owner.inverse.delta,sign)
        qlo,qhi=current.ends(solved['actual_log_lambda']);zl,zh=current.ends(solved['Z'])
        assert qlo<=-7<=qhi and max(zl,current.ends(z)[0])<=min(zh,current.ends(z)[1])
        numeric.append(dict(Z=Z,original_pure_log_root_overlap=True,current_original_delta_unchanged=True))
    rejected=[]
    tests=(
        ('unregistered_copied_input',lambda:observed_owner.invert(copy.copy(requests['exterior']))),
        ('foreign_owner_input',lambda:owner.invert(requests['exterior'])),
        ('foreign_owner_inverse',lambda:owner.locate(inverses['exterior'])),
        ('copied_inverse',lambda:observed_owner.locate(dict(inverses['exterior']))),
        ('copied_source_delivery',lambda:observed_owner.physical_rows(dict(delivered['exterior']))),
        ('foreign_owner_source_delivery',lambda:owner.physical_rows(delivered['exterior'])),
        ('axial_endpoint',lambda:owner.physical_input('heat_exterior','4','1')),
        ('unknown_chart',lambda:owner.physical_input('old_practical_heat','4','.5')),
        ('imported_numeric_native_box',lambda:owner.physical_input('heat_exterior',owner.ctx.mpf([4,5]),'.5')),
        ('unknown_graph_log_oracle',lambda:owner.physical_input('heat_exterior','4','.5',owner.graph.symbol('unadmitted_log_time'))),
        ('different_source_chart',lambda:observed_owner.source(inverses['exterior'],'heat_collar')))
    for name,call in tests:
        try:call()
        except (TypeError,ValueError):rejected.append(name)
        else:raise AssertionError('Unadmitted correlated physical source accepted: '+name)
    changed=owner.physical_input('heat_exterior','4','.521');changed.geometry['logR']=0
    try:owner.invert(changed)
    except ValueError:rejected.append('mutated_live_source_geometry')
    else:raise AssertionError('Mutated source input accepted')
    modified=dict(inverses['exterior']);modified['directed_inverse_mapping']=dict(modified['directed_inverse_mapping'])
    modified['directed_inverse_mapping']['Z']=owner.ctx.mpf(0)
    try:observed_owner.locate(modified)
    except ValueError:rejected.append('mutated_imported_inverse_Z')
    else:raise AssertionError('Mutated inverse accepted')
    assert rows==840 and times==24
    result=dict(all_passed=True,source_family=owner.family_record,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        independent_original_inverse_algebra=independent,current_exact_source_graph_identities=exact,
        independent_finite_log_original_root_fixtures=numeric,
        actual_correlated_physical_source_summaries=summaries,
        actual_correlated_physical_spatial_rows=rows,actual_correlated_fixed_x_time_rows=times,
        actual_signed_exact_physical_scale_groups=groups,
        unchanged_original_physical_mapper_AST_binding=owner.inverse.forward_unit_binding,
        original_physical_mapper_function_reused_unchanged=True,
        physical_input_family_is_restricted_not_an_arbitrary_x_y_z_t_oracle=True,
        physical_absolute_relative_accuracy_materialization_global_axis_stress_and_recursion_open=True,
        rejected_unadmitted_inputs=rejected,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            current.PREFIX+'current_original_Rp_segmented_radius_check.py':current.sha(current.PREFIX+'current_original_Rp_segmented_radius_check.py')},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.report(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_CORRELATED_PHYSICAL_SOURCE actual original-scale inputs and physical rows',flush=True)
    return result
