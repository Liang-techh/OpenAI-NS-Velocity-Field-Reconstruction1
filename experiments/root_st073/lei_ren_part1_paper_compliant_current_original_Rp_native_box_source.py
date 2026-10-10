"""Whole-box source callbacks on the accepted original fifteen-chart graph.

The original coefficient algorithms are unchanged. Narrow, checked AST
adapters remove only wrapper scalar coercions and pass the entire directed
box to the already interval-native providers. Exact coordinate functions,
source units and Jacobians stay separate from their numeric enclosures.
This is a similarity-source callback, not a global physical point field.
"""
import ast
import copy
from dataclasses import dataclass
from fractions import Fraction
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
import time
from types import MethodType

import lei_ren_part1_paper_compliant_current_original_Rp_native_chart_locator as locator
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

mixed=locator.physical.mixed
pulse=mixed.pulse
HERE,PREFIX,sha=locator.HERE,locator.PREFIX,locator.sha
NAME=PREFIX+'current_original_Rp_native_box_source.json.gz'
RECEIPT=PREFIX+'current_original_Rp_native_box_source_check.json'
GATES=('current_original_Rp_fifteen_chart_whole_box_source_callbacks_installed',
       'current_original_Rp_native_box_exact_function_scales_preserved')
OPEN=tuple(key for key in locator.OPEN if key not in GATES)
CHARTS=locator.CHARTS
ends=locator.ends


def _dump(node):return ast.dump(node,include_attributes=False)


def adapter(function,replacements):
    """Fail closed unless every exact, scalar-only replacement matches once."""
    original=inspect.unwrap(function)
    tree=ast.parse(textwrap.dedent(inspect.getsource(original)))
    before=_dump(tree);program=copy.deepcopy(tree);changes=[]
    for old,new in replacements:
        wanted=ast.parse(old).body
        if len(wanted)!=1:raise ValueError('One wrapper statement per adapter change required')
        needle=_dump(wanted[0]);found=[]
        for parent in ast.walk(program):
            for field,value in ast.iter_fields(parent):
                if isinstance(value,list):
                    for index,item in enumerate(value):
                        if isinstance(item,ast.stmt) and _dump(item)==needle:
                            found.append((value,index))
        if len(found)!=1:raise ValueError('Current wrapper source changed: '+old)
        body,index=found[0];body[index:index+1]=ast.parse(new).body if new else []
        changes.append(dict(original_statement=old,box_statement=new))
    namespace=dict(original.__globals__)
    exec(compile(ast.fix_missing_locations(program),inspect.getsourcefile(original),'exec'),namespace)
    return namespace[original.__name__],dict(
        current_wrapper=original.__module__+'.'+original.__qualname__,
        current_wrapper_file_sha256=sha(Path(inspect.getsourcefile(original)).name),
        original_AST_sha256=hashlib.sha256(before.encode()).hexdigest(),
        adapter_AST_sha256=hashlib.sha256(_dump(program).encode()).hexdigest(),
        only_listed_scalar_coercions_or_coordinate_forwarding_changed=changes)


def programs():
    specifications={
        'pulse_scale':(pulse.CurrentOriginalRpRawPulseTransport.scale,[
            ('value = radius.exact_coordinate(coordinate)','v = coordinate.function'),
            ('v = g.constant(value)','')]),
        'pulse_view':(pulse.CurrentOriginalRpRawPulseTransport.pulse_view,[
            ('value = radius.exact_coordinate(coordinate)','v = coordinate.enclosure'),
            ('v = c.mpf(value.numerator) / value.denominator','')]),
        'raw_evaluate':(pulse.raw.CurrentOriginalRpRawHistoryTransport.evaluate,[
            ("packet = self.post.evaluate(chart, Z, coordinate)['source_packet']",
             "packet = self.post.evaluate(chart, Z, coordinate.enclosure)['source_packet']"),
            ('hats = self.post.histories(chart, Z, coordinate, packet)',
             'hats = self.post.histories(chart, Z, coordinate.enclosure, packet)')]),
        'heat_evaluate':(pulse.closed.CurrentOriginalRpSameRepairHeatClosure.evaluate,[
            ('packet = self.pressure.evaluate(chart, Z, coordinate)',
             'packet = self.pressure.evaluate(chart, Z, coordinate.enclosure)'),
            ('value = raw_source.radius.exact_coordinate(coordinate)','t = coordinate.enclosure'),
            ('t = c.mpf(value.numerator) / value.denominator','')]),
        'heat_velocity_rows':(mixed.CurrentOriginalRpMixedTransport.velocity_rows,[
            ('value = pulse.radius.exact_coordinate(coordinate)','t = coordinate.enclosure'),
            ('t = c.mpf(value.numerator) / value.denominator','')])}
    functions={};proof={}
    for name,(function,changes) in specifications.items():
        functions[name],proof[name]=adapter(function,changes)
    return functions,proof


@dataclass(frozen=True)
class NativeBoxCoordinate:
    chart: str
    function: object
    enclosure: object


class _Proxy:
    def __init__(self,actual):self.actual=actual;self.call_trace=[]
    def __getattr__(self,name):return getattr(self.actual,name)
    def bind(self,name,function):setattr(self,name,MethodType(function,self))


class _BoxRadius(_Proxy):
    def __init__(self,owner):
        super().__init__(owner.before.radius);self.box_owner=owner

    def geometry(self,chart,coordinate):
        self.box_owner._validate_token(chart,coordinate)
        actual=self.actual._map(chart,coordinate.function)
        return dict(chart=chart,exact_native_coordinate_function=coordinate.function.node,
            directed_native_coordinate=coordinate.enclosure,
            **{key:ref.node for key,ref in actual.items()},
            native_coordinate_scope='source function variable or admitted live physical inverse; whole box enclosure',
            absolute_radius_is_a_lazy_function=True,short_offset_not_added_to_huge_absolute_value=True,
            original_waiting_integrals_not_replaced_by_root_box=True,
            numerical_absolute_radius_point_value_installed=False)

    def local_step(self,chart,left,right):
        # Existing raw wrappers use local_step(v,v) only to request the
        # constant affine Jacobian. An interval difference is not substituted.
        self.box_owner._validate_token(chart,left)
        if right is not left:raise ValueError('Box adapter supplies only the same-coordinate Jacobian query')
        J=self.actual._map(chart,left.function)['native_to_log_radius_jacobian']
        bound=self.box_owner.before.reader().at(J)
        if ends(bound)[0]<=0:raise ValueError('Positive exact current native Jacobian required')
        return dict(exact_log_radius_increment=self.actual.graph.zero.node,
            native_to_log_radius_jacobian_bound=bound,log_radius_increment_bound=self.actual.ctx.mpf(0),
            constant_origin_cancelled_before_arithmetic=True,enclosure_not_used_as_defining_function=True)


class CurrentOriginalRpNativeBoxSource:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else locator.CurrentOriginalRpNativeChartLocator()
        if type(self.before) is not locator.CurrentOriginalRpNativeChartLocator or not self.before.acceptance_loaded:
            raise ValueError('Accepted current original native locator required')
        self.actual=self.before.before.before.transport
        if type(self.actual) is not mixed.CurrentOriginalRpMixedTransport or not self.actual.acceptance_loaded:
            raise ValueError('Accepted same current mixed source required')
        self.graph=self.before.graph;self.ctx=self.before.ctx;self.family_record=self.before.family_record
        self.hashes=dict(self.before.hashes);self.acceptance_loaded=False;self.call_trace=[];self._tokens={}
        for name in (locator.NAME,locator.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        functions,self.adapter_proof=programs();self._adapter_functions=functions
        self.box_radius=_BoxRadius(self)
        self.raw=_Proxy(self.actual.owner.raw);self.raw.radius=self.box_radius
        self.raw.bind('evaluate',functions['raw_evaluate'])
        self.heat=_Proxy(self.actual.owner.closed);self.heat.before=self.raw
        self.heat.bind('evaluate',functions['heat_evaluate'])
        self.pulse=_Proxy(self.actual.owner);self.pulse.radius=self.box_radius
        self.pulse.raw=self.raw;self.pulse.closed=self.heat
        self.pulse.bind('scale',functions['pulse_scale']);self.pulse.bind('pulse_view',functions['pulse_view'])
        for method in ('factor','expression_owner','evaluate'):
            self.pulse.bind(method,getattr(pulse.CurrentOriginalRpRawPulseTransport,method))
        self.transport=_Proxy(self.actual);self.transport.owner=self.pulse
        self.transport.bind('velocity_rows',functions['heat_velocity_rows'])
        for method in ('factor','evaluate'):
            self.transport.bind(method,getattr(mixed.CurrentOriginalRpMixedTransport,method))
        self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or any(receipt[key] for key in OPEN):
                raise ValueError('Current native box source receipt or scope differs')
            if receipt['source_family']!=self.family_record or receipt['current_scalar_wrapper_adapter_proof']!=self.adapter_proof:
                raise ValueError('Current source family or scalar wrapper binding differs')
            pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        result=dict(accepted_current_native_locator=self.before.acceptance_loaded and all(self.before.assert_graph().values()),
            accepted_current_mixed_source=self.actual.acceptance_loaded and all(self.actual.assert_graph().values()),
            same_current_expression_graph=self.graph is self.actual.graph is self.before.graph,
            same_current_interval_context=self.ctx is self.actual.ctx is self.before.ctx,
            same_original_fifteen_chart_provider=self.actual is self.before.before.before.transport,
            same_selected_repair_and_closed_heat=self.heat.actual is self.actual.owner.closed
                and self.pulse.actual is self.actual.owner and self.raw.actual is self.actual.owner.raw,
            current_box_proxy_wiring=self.transport.owner is self.pulse and self.pulse.raw is self.raw
                and self.pulse.closed is self.heat and self.heat.before is self.raw
                and self.pulse.radius is self.raw.radius is self.box_radius and self.box_radius.box_owner is self,
            current_scalar_adapter_methods_unchanged=all(getattr(proxy,name).__func__ is self._adapter_functions[key]
                for proxy,name,key in ((self.pulse,'scale','pulse_scale'),(self.pulse,'pulse_view','pulse_view'),
                    (self.raw,'evaluate','raw_evaluate'),(self.heat,'evaluate','heat_evaluate'),
                    (self.transport,'velocity_rows','heat_velocity_rows'))))
        if not all(result.values()):raise ValueError('Current box source graph differs: '+str(result))
        return result

    def _admit(self,chart,Z,native):
        if chart not in CHARTS:raise ValueError('Actual current source chart required')
        zl,zh=ends(Z);nl,nh=ends(native)
        if not (-1<zl<=zh<1):raise ValueError('Strict original axial domain required for the entire Z box')
        reader=self.before.reader();lower,upper=self.before.domain(chart)
        ll,lh=ends(reader.at(lower));rl=None if upper is None else ends(reader.at(upper))[0]
        if nl<lh or rl is not None and nh>rl:
            raise ValueError('Entire source box must lie in one current native chart; split a crossing box')
        # The angular provider has different bump-support formulas on each
        # interval. Do not silently select a branch for a support crossing.
        if chart=='outer_angular':
            for center in (Fraction(-3),Fraction(-1)):
                for endpoint in (center-Fraction(3,20),center+Fraction(3,20)):
                    point=locator.inverse.box(self.ctx,endpoint)
                    if nl<ends(point)[1] and nh>ends(point)[0]:
                        raise ValueError('Split an angular box at its exact bump-support boundary')

    def _register(self,chart,function,enclosure,Z,Zfunction,provenance):
        self._admit(chart,Z,enclosure)
        token=NativeBoxCoordinate(chart,function,enclosure)
        self._tokens[id(token)]=(token,function.node,enclosure._mpi_,Z._mpi_,Zfunction.node,copy.deepcopy(provenance))
        return token

    def _validate_token(self,chart,coordinate):
        if type(coordinate) is not NativeBoxCoordinate or id(coordinate) not in self._tokens:
            raise ValueError('Live coordinate token issued by this current box owner required')
        row=self._tokens[id(coordinate)]
        if row[0] is not coordinate or coordinate.chart!=chart or coordinate.function.graph is not self.graph or \
                (coordinate.function.node,coordinate.enclosure._mpi_)!=(row[1],row[2]):
            raise ValueError('Unchanged same-chart current coordinate token required')

    def _evaluate(self,token,Z,Zfunction,provenance):
        self.assert_graph();self._validate_token(token.chart,token)
        record=self._tokens[id(token)]
        if (Z._mpi_,Zfunction.node,provenance)!=(record[3],record[4],record[5]):
            raise ValueError('Unchanged current axial coordinate and source provenance required')
        view=self.transport.evaluate(token.chart,Z,token)
        view.update(exact_axial_coordinate_function=Zfunction.node,directed_axial_coordinate=Z,
            source_input_provenance=provenance,whole_native_and_axial_boxes_passed_to_original_algorithms=True,
            scalar_midpoint_or_endpoint_not_used_as_source_value=True,
            box_geometry_retains_exact_source_functions=True,
            coordinate_error_is_propagated_by_directed_source_arithmetic=True,
            uniform_source_accuracy_or_global_physical_field_certified=False,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        self.call_trace.append(dict(chart=token.chart,native_box=token.enclosure,Z_box=Z,
            exact_native_function=token.function.node,exact_Z_function=Zfunction.node,
            current_original_whole_box_source_called=True,absolute_scale_materialized=False))
        return view

    @source_precision
    def source_cell(self,chart,Z_left,Z_right,native_left,native_right):
        """An independent original-coordinate cell, not a physical point query."""
        values=tuple(locator.inverse.exact_scalar(v) for v in (Z_left,Z_right,native_left,native_right))
        a,b,l,r=values
        if a>b or l>r:raise ValueError('Ordered exact source cell endpoints required')
        c=self.ctx
        Z=c.mpf([locator.inverse.box(c,a).a,locator.inverse.box(c,b).b])
        native=c.mpf([locator.inverse.box(c,l).a,locator.inverse.box(c,r).b])
        self._admit(chart,Z,native)
        exact=dict(chart=chart,axial_bounds=[str(a),str(b)],native_bounds=[str(l),str(r)],
            kind='independent_original_similarity_coordinate_cell')
        key=hashlib.sha256(json.dumps(exact,sort_keys=True).encode()).hexdigest()
        function=self.graph.symbol('current_original_native_cell_'+key)
        Zfunction=self.graph.symbol('current_original_axial_cell_'+key)
        token=self._register(chart,function,native,Z,Zfunction,exact)
        return self._evaluate(token,Z,Zfunction,exact)

    @source_precision
    def from_inverse(self,view,chart):
        """Consume a witnessed whole inverse box only after full chart admission."""
        location=self.before.locate_inverse(view)
        if chart not in location['contained_charts'] or not location['physical_Z_numeric_interior_resolved']:
            raise ValueError('Whole inverse enclosure must be contained in one current chart and strict Z domain')
        row=location['chart_rows'][chart];Z=self.ctx.mpf(view['directed_inverse_mapping']['Z'])
        Zfunction=view['coordinate_functions']['Z']
        provenance=dict(kind='unchanged_live_current_physical_inverse',
            exact_input_record=copy.deepcopy(view['exact_input_record']),exact_logR=location['exact_input_logR'].node)
        token=self._register(chart,row['exact_native_coordinate_function'],row['directed_native_coordinate'],Z,Zfunction,provenance)
        return self._evaluate(token,Z,Zfunction,provenance)


def report(view):return locator.report(mixed.view_report(view))


@source_precision
def run(before=None):
    began=time.monotonic();owner=CurrentOriginalRpNativeBoxSource(before,require_checked=False);views={}
    # Nonzero box widths exercise both coordinate dimensions. Tiny native
    # widths retain usefulness even for the very large exact source Jacobians.
    dz=Fraction(1,10**10);dv=Fraction(1,10**30);z=Fraction(521,1000)
    for chart,anchor in locator.ANCHORS.items():
        n=Fraction(anchor);views[chart]=owner.source_cell(chart,z-dz,z+dz,n-dv,n+dv)
        print('Actual current whole native/Z box source:',chart,flush=True)
    result=dict(source_family=owner.family_record,actual_current_graph=owner.assert_graph(),
        current_scalar_wrapper_adapter_proof=owner.adapter_proof,
        actual_fifteen_whole_box_source_views={chart:report(view) for chart,view in views.items()},
        actual_whole_box_source_call_trace=owner.call_trace,
        exact_current_box_source_graph=owner.graph.nodes,
        source_scope='Directed similarity source coefficients and total-order-four mixed rows; exact log scale factors. Correlated original-scale physical query and materialization remain open',
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(locator.report(result),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,views
