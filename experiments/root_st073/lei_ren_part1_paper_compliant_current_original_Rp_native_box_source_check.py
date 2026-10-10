"""Focused whole-coordinate-box forwarding and original unit checks."""
import copy
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_native_box_source as current
from lei_ren_part1_paper_compliant_current_original_Rp_segmented_radius_check import interpreter
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def _overlap(a,b):
    lo,hi=current.ends(a);l,h=current.ends(b)
    return max(lo,l)<=min(hi,h)


def exact_scale_checks(owner,views):
    g=owner.graph;r=owner.before.radius;identities=0
    for chart,view in views.items():
        coordinate=current.pulse.radius.FunctionRef(g,view['geometry']['exact_native_coordinate_function'])
        anchor=Fraction(current.locator.ANCHORS[chart])
        bindings={ref.node:s.Symbol('original_'+name,positive=True) for name,ref in r.parameters.items()}
        bindings.update({ref.node:s.Symbol('source_'+name,positive=True) for name,ref in r.functions.items()})
        bindings[owner.actual.owner.raw.logU0.node]=s.Symbol('source_inlet_logU0',real=True)
        bindings[coordinate.node]=s.Rational(anchor.numerator,anchor.denominator)
        q=interpreter(r,False,bindings)
        exact_geometry=r.geometry(chart,anchor)
        for key in ('origin','pulse_term','offset','logR','native_to_log_radius_jacobian'):
            assert s.cancel(s.expand(q.at(current.pulse.radius.FunctionRef(g,view['geometry'][key]))-
                q.at(current.pulse.radius.FunctionRef(g,exact_geometry[key]))))==0,(chart,key)
            identities+=1
        assert view['geometry']['native_to_log_radius_jacobian']==r.maps[chart]['native_to_log_radius_jacobian'].node
        for name,rows in view['log_radius_mixed_rows'].items():
            assert len(rows)==15
            for label,row in rows.items():
                k,j=row.derivative
                expected=(owner.actual.owner.scale(chart,anchor,row.powers[:3]) if chart in current.pulse.PULSE else
                    owner.actual.owner.raw.scale(chart,exact_geometry,row.powers[:3]))
                got=dict(row.log_scale_parts)
                for part,value in expected.items():
                    assert s.cancel(s.expand(q.at(got[part])-q.at(value)))==0,(chart,name,label,part)
                    identities+=1
                assert row.source_units[-1]=='native_Jacobian' and row.powers[-1]==0
                native=view['native_coordinate_mixed_rows'][name]['n'+str(k)+'_Z'+str(j)]
                assert native.derivative==(k,j) and native.powers[-1]==k
                assert native.coefficients[0]._mpi_==row.coefficients[0]._mpi_
                J=current.pulse.radius.FunctionRef(g,view['geometry']['native_to_log_radius_jacobian'])
                wanted=g.mul(g.constant(k),g.unary('log',J)) if k else g.zero
                assert dict(native.log_scale_parts)['native_coordinate_Jacobian']==wanted
                identities+=1
    return dict(passed=True,current_exact_geometry_scale_and_native_Jacobian_identities=identities,
        symbolic_cell_variable_specialization_agrees_with_unchanged_scalar_source=True,
        native_derivative_units_retain_exact_J_power=True)


@source_precision
def run(before,observed_owner,views):
    began=time.monotonic();candidate=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.CurrentOriginalRpNativeBoxSource(before,require_checked=False)
    assert type(observed_owner) is type(owner) and observed_owner.before is before
    assert observed_owner.hashes==owner.hashes
    assert all(observed_owner.assert_graph().values())
    for name,digest in candidate['input_hashes'].items():assert current.sha(name)==digest,name
    assert candidate['source_family']==owner.family_record and not any(candidate[key] for key in current.GATES+current.OPEN)
    assert candidate['current_scalar_wrapper_adapter_proof']==owner.adapter_proof
    assert candidate['actual_current_graph']==owner.assert_graph()
    graph=candidate['exact_current_box_source_graph'];assert graph==owner.graph.nodes[:len(graph)]
    assert set(views)==set(current.CHARTS)
    summaries={};row_count=0;scalar_count=0;overlap_count=0
    for chart,view in views.items():
        assert current.report(view)==candidate['actual_fifteen_whole_box_source_views'][chart]
        assert view['whole_native_and_axial_boxes_passed_to_original_algorithms']
        assert view['coordinate_error_is_propagated_by_directed_source_arithmetic']
        assert not view['uniform_source_accuracy_or_global_physical_field_certified']
        assert not any(view[key] for key in current.GATES+current.OPEN)
        nl,nh=current.ends(view['geometry']['directed_native_coordinate'])
        zl,zh=current.ends(view['directed_axial_coordinate'])
        assert nl<nh and -1<zl<zh<1
        # Calls on a genuine exact source coordinate serve only as a
        # consistency comparison. They never define the whole-box result.
        scalar=(owner.actual.evaluate(chart,'521/1000',current.locator.ANCHORS[chart])
            if chart in ('pulse_exit','steep_power','heat_exterior') else None)
        if scalar is not None:scalar_count+=1
        widths=[];max_row=0
        if scalar is not None:assert set(scalar['log_radius_mixed_rows'])==set(view['log_radius_mixed_rows'])
        for name,rows in view['log_radius_mixed_rows'].items():
            for key,row in rows.items():
                box=row.coefficients[0];lo,hi=current.ends(box)
                assert mp.isfinite(lo) and mp.isfinite(hi) and lo<=hi
                if scalar is not None:
                    assert _overlap(box,scalar['log_radius_mixed_rows'][name][key].coefficients[0]),(chart,name,key)
                    overlap_count+=1
                row_count+=1
                widths.append(owner.ctx.mpf(hi)-owner.ctx.mpf(lo));max_row=max(max_row,abs(lo),abs(hi))
        packet=view['original_forward_source_packet']
        assert current.ends(packet['Z'])==(zl,zh)
        assert view['geometry'].get('exact_native_coordinate') is None
        summaries[chart]=dict(native_box_width=owner.ctx.mpf(nh)-owner.ctx.mpf(nl),
            axial_box_width=owner.ctx.mpf(zh)-owner.ctx.mpf(zl),
            ordinary_logR_Z_rows=len(widths),maximum_signed_scaled_row_width=max(widths),
            maximum_signed_scaled_magnitude=owner.ctx.mpf(max_row),
            source_packet_uses_whole_Z_box=True,unchanged_scalar_source_rows_overlap=scalar is not None)
        print('Checked current whole-box source units/rows:',chart,flush=True)
    exact=exact_scale_checks(owner,views)
    rejected=[]
    tests=(
        ('native_chart_crossing',lambda:owner.source_cell('flatten','.5','.6','99','101')),
        ('axial_domain_crossing',lambda:owner.source_cell('flatten','.9','1.01','49','51')),
        ('reversed_cell',lambda:owner.source_cell('flatten','.6','.5','49','51')),
        ('unknown_chart',lambda:owner.source_cell('old_practical_heat','.5','.6','4','5')),
        ('angular_support_crossing',lambda:owner.source_cell('outer_angular','.5','.6','-16/5','-31/10')),
        ('collar_tail_crossing',lambda:owner.source_cell('heat_collar','.5','.6','29/10','31/10')),
        ('imported_native_interval',lambda:owner.source_cell('flatten','.5','.6',owner.ctx.mpf([49,50]),51)),
        ('unregistered_token',lambda:owner.box_radius.geometry('flatten',current.NativeBoxCoordinate(
            'flatten',owner.graph.constant(50),owner.ctx.mpf(50)))),
        ('before_Rp_inverse_not_an_outer_source',lambda:owner.from_inverse(owner.before.before.cartesian('.8','-.3','.2','.95'),'heat_exterior')))
    for name,call in tests:
        try:call()
        except (ValueError,TypeError):rejected.append(name)
        else:raise AssertionError('Unadmitted source box accepted: '+name)
    token_row=next(iter(observed_owner._tokens.values()));token=token_row[0]
    try:owner.box_radius.geometry(token.chart,token)
    except ValueError:rejected.append('foreign_owner_live_token')
    else:raise AssertionError('Foreign box owner token accepted')
    changed=copy.copy(observed_owner);changed.transport=current._Proxy(observed_owner.transport.actual)
    changed.transport.__dict__.update(observed_owner.transport.__dict__)
    changed.transport.owner=owner.pulse
    try:changed.assert_graph()
    except ValueError:rejected.append('substituted_proxy_source_owner')
    else:raise AssertionError('Substituted proxy source accepted')
    assert len(candidate['actual_whole_box_source_call_trace'])==15
    assert all(row['current_original_whole_box_source_called'] and not row['absolute_scale_materialized']
        for row in candidate['actual_whole_box_source_call_trace'])
    result=dict(all_passed=True,source_family=owner.family_record,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        current_scalar_wrapper_adapter_proof=owner.adapter_proof,
        current_exact_source_units_and_geometry_checks=exact,
        actual_whole_box_source_summaries=summaries,
        actual_whole_box_logR_Z_ordinary_rows=row_count,
        actual_whole_box_native_Z_ordinary_rows=row_count,
        unchanged_original_scalar_source_calls=scalar_count,
        unchanged_original_scalar_row_overlap_checks=overlap_count,
        scalar_comparisons_are_not_uniform_error_proofs=True,
        absolute_or_relative_physical_accuracy_not_certified=True,
        actual_original_scale_physical_inverse_source_delivery_still_open=True,
        rejected_unadmitted_inputs=rejected,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            current.PREFIX+'current_original_Rp_segmented_radius_check.py':current.sha(current.PREFIX+'current_original_Rp_segmented_radius_check.py')},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.locator.report(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_NATIVE_BOX_SOURCE fifteen unchanged original whole-box callbacks',flush=True)
    return result
