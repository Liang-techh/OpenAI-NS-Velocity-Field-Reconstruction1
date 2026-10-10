"""Source-bound original-scale physical log inputs and exact inverse delivery.

This restricted input family preserves log_tau=-logR_source+finite_offset.
The original monotone implicit inverse remains an exact root node. Its
proved equality to the original forward coordinate is a separate alias,
used before directed arithmetic so astronomical common terms cancel.
Actual source callbacks consume the resulting live coordinate enclosures.
Unrestricted physical inputs and physical value materialization stay open.
"""
import copy
from dataclasses import dataclass
from fractions import Fraction
import gzip
import hashlib
import json
from pathlib import Path
import time
from types import MethodType

import lei_ren_part1_paper_compliant_current_original_Rp_native_box_source as box
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

locator=box.locator
inverse=locator.inverse
physical=locator.physical
HERE,PREFIX,sha=box.HERE,box.PREFIX,box.sha
NAME=PREFIX+'current_original_Rp_correlated_physical_source.json.gz'
RECEIPT=PREFIX+'current_original_Rp_correlated_physical_source_check.json'
GATES=('current_original_Rp_source_bound_correlated_physical_log_input_family_installed',
       'current_original_Rp_correlated_unit_inverse_identity_and_cancellation_installed',
       'current_original_Rp_correlated_physical_inverse_actual_source_delivery_installed',
       'current_original_Rp_correlated_physical_scaled_spatial_time_rows_installed')
OPEN=tuple(key for key in box.OPEN if key not in GATES)
ends=locator.ends


class CorrelatedGraphBounds(locator.CancelledGraphBounds):
    """Only this live caller's exact equality proofs supply the aliases."""
    def __init__(self,reader,aliases):
        super().__init__(reader.graph,reader.ctx,reader.bindings,reader.allowed_exponentials)
        self.aliases=dict(aliases)

    def polynomial(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.aliases:return super().polynomial(self.aliases[i])
        if i in self.bindings:return {(i,):Fraction(1)}
        node=self.graph.nodes[i]
        if node['operation']=='positive_quotient':
            # Normalize a/positive_d to a times one common reciprocal.
            # Original source maps use both v/mu and v*(1/mu).
            reciprocal=self.graph.quotient(self.graph.one,box.pulse.radius.FunctionRef(self.graph,node['denominator']),
                'current original positive denominator; exact reciprocal normalization')
            if i==reciprocal.node:return {(i,):Fraction(1)}
            return {tuple(sorted(key+(reciprocal.node,))):value
                for key,value in self.polynomial(node['numerator']).items()}
        return super().polynomial(i)

    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.aliases:return super().at(self.aliases[i])
        return super().at(i)


@dataclass(frozen=True)
class CorrelatedPhysicalInput:
    chart: str
    native: Fraction
    Z: Fraction
    finite_log_time_offset: Fraction
    theta: Fraction
    geometry: dict
    forward_coordinates: dict
    physical_log_abs_z: object


class CurrentOriginalRpCorrelatedPhysicalSource:
    @source_precision
    def __init__(self,before=None,require_checked=True):
        self.before=before if before is not None else box.CurrentOriginalRpNativeBoxSource()
        if type(self.before) is not box.CurrentOriginalRpNativeBoxSource or not self.before.acceptance_loaded:
            raise ValueError('Accepted current whole-box source owner required')
        self.locator=self.before.before;self.inverse=self.locator.before;self.physical=self.inverse.before
        self.graph=self.before.graph;self.ctx=self.before.ctx;self.radius=self.locator.radius
        self.family_record=self.before.family_record;self.original_N_definition=self.inverse.original_N_definition
        self.hashes=dict(self.before.hashes);self.acceptance_loaded=False;self.call_trace=[]
        self._inputs={};self._inverse_views={};self._deliveries={}
        for name in (box.NAME,box.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[key] for key in GATES) or any(receipt[key] for key in OPEN):
                raise ValueError('Current correlated physical receipt or scope differs')
            if receipt['source_family']!=self.family_record:raise ValueError('Current correlated source family differs')
            box.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    def assert_graph(self):
        result=dict(accepted_current_whole_box_source=self.before.acceptance_loaded and all(self.before.assert_graph().values()),
            same_current_original_inverse_and_physical_map=self.inverse is self.locator.before
                and self.physical is self.inverse.before and all(self.inverse.assert_graph().values()),
            same_original_graph_and_context=self.graph is self.inverse.graph is self.physical.graph
                and self.ctx is self.before.ctx is self.inverse.ctx,
            unchanged_original_N_and_delta=self.original_N_definition==self.physical.original_N_definition
                and self.inverse.delta_function==self.radius.functions['delta'],
            original_unit_viscosity_forward_binding=self.inverse.forward_unit_binding==inverse.current_forward_unit_binding())
        if not all(result.values()):raise ValueError('Current correlated physical graph differs: '+str(result))
        return result

    def _input_record(self,request):
        if any(value.graph is not self.graph for value in request.forward_coordinates.values()) or \
                request.physical_log_abs_z is not None and request.physical_log_abs_z.graph is not self.graph:
            raise ValueError('Same original physical graph inputs required')
        return dict(chart=request.chart,native=str(request.native),Z=str(request.Z),
            finite_log_time_offset=str(request.finite_log_time_offset),theta=str(request.theta),
            source_geometry=request.geometry,
            exact_forward_coordinate_nodes={key:value.node for key,value in request.forward_coordinates.items()},
            exact_physical_log_abs_z=None if request.physical_log_abs_z is None else request.physical_log_abs_z.node)

    def _fingerprint(self,value):
        return hashlib.sha256(json.dumps(locator.report(value),sort_keys=True,separators=(',',':')).encode()).hexdigest()

    def _validate_input(self,request):
        self.assert_graph();record=self._inputs.get(id(request))
        if type(request) is not CorrelatedPhysicalInput or record is None or record[0] is not request or \
                record[1]!=self._fingerprint(self._input_record(request)):
            raise ValueError('Unchanged live correlated physical input from this current owner required')
        expected=self.radius.geometry(request.chart,request.native)
        if expected!=request.geometry:raise ValueError('Original source geometry changed')
        g=self.graph;logR=box.pulse.radius.FunctionRef(g,expected['logR'])
        lt=g.add(g.neg(logR),g.constant(request.finite_log_time_offset))
        actual=self.physical.coordinates(expected,request.Z,lt,request.theta)
        if request.forward_coordinates!=actual:raise ValueError('Original forward coordinate input changed')
        logz=None if request.Z==0 else g.add(g.unary('log',g.constant(abs(request.Z))),
            g.mul(g.sub(g.one,self.inverse.delta_function),actual['loglambda']))
        if request.physical_log_abs_z!=logz:raise ValueError('Exact physical axial logarithm changed')

    @source_precision
    def physical_input(self,chart,native,Z,finite_log_time_offset='0',theta='7/10'):
        """A restricted exact physical input family, with original source scale.

        The input's time and axial log may be far outside decimal or binary
        physical-value materialization. Their exact functions remain finite.
        """
        self.assert_graph();g=self.graph
        n,z,offset,angle=map(inverse.exact_scalar,(native,Z,finite_log_time_offset,theta))
        if not -1<z<1:raise ValueError('Strict original axial coordinate required')
        geometry=self.radius.geometry(chart,n)
        logR=box.pulse.radius.FunctionRef(g,geometry['logR'])
        lt=g.add(g.neg(logR),g.constant(offset))
        coords=self.physical.coordinates(geometry,z,lt,angle)
        reader=self.locator.reader({self.inverse.delta_function.node:self.inverse.delta})
        if ends(reader.at(lt))[1]>=0:raise ValueError('This near-critical family requires exact 0<tau<1')
        logz=None if z==0 else g.add(g.unary('log',g.constant(abs(z))),
            g.mul(g.sub(g.one,self.inverse.delta_function),coords['loglambda']))
        request=CorrelatedPhysicalInput(chart,n,z,offset,angle,geometry,coords,logz)
        self._inputs[id(request)]=(request,self._fingerprint(self._input_record(request)))
        return request

    def _identity(self,request,refs):
        """Check source algebra before registering any root-coordinate alias."""
        self._validate_input(request);g=self.graph;c=request.forward_coordinates
        q=c['loglambda'];z=c['Z'];d=c['delta'];complement=g.sub(g.one,g.mul(z,z))
        reader=self.locator.reader()
        required=(g.add(g.mul(g.constant(2),q),g.neg(c['log_tau']),g.unary('log',complement)),
            g.sub(c['log_r'],g.add(q,g.mul(g.constant('1/2'),g.add(g.unary('log',g.constant(2)),c['logR'])))))
        if any(reader.polynomial(value) for value in required):raise ValueError('Original forward algebra differs')
        if request.Z:
            axial=g.sub(request.physical_log_abs_z,g.add(g.unary('log',g.constant(abs(request.Z))),
                g.mul(g.sub(g.one,d),q)))
            if reader.polynomial(axial):raise ValueError('Original physical axial algebra differs')
            node=g.nodes[refs['log_lambda'].node]
            if node['operation']!='exact_original_physical_log_lambda_inverse' or \
                    (node['log_tau'],node['log_abs_z'],node['delta'],node['viscosity'])!=(
                    c['log_tau'].node,request.physical_log_abs_z.node,d.node,g.one.node):
                raise ValueError('Original implicit inverse root node required')
        if not 1-request.Z**2>0 or (1-request.Z**2)+request.Z**2!=1:
            raise ValueError('Positive original root factorization required')
        # The displayed factorization follows the checked forward exponents:
        # exp(log_tau)=(1-Z^2)exp(2q),
        # exp(2log_abs_z+2delta*q)=Z^2 exp(2q).
        # Their sum is exp(2q). The admitted original theorem gives uniqueness.
        proof=g.node('exact_original_correlated_physical_inverse_identity',
            original_inverse_root=refs['log_lambda'].node,forward_log_lambda=q.node,
            exact_source_Z=z.node,physical_log_r=c['log_r'].node,physical_log_tau=c['log_tau'].node,
            physical_log_abs_z=None if request.physical_log_abs_z is None else request.physical_log_abs_z.node,
            source_logR=c['logR'].node,original_delta=d.node,viscosity=g.one.node,
            positive_time_factor='1-Z^2>0; exp(log_tau)=(1-Z^2)exp(2q)',
            axial_factor='Z^2 exp(2q), or zero on exact Z=0 branch',
            root_factorization='(1-Z^2)+Z^2=1',
            uniqueness='original implicit F_prime>=2(1-delta)>0; 0<delta<1',
            source_family=self.family_record,current_forward_AST_sha256=self.inverse.forward_unit_binding[
                'current_forward_coordinates_AST_sha256'])
        aliases={refs['log_lambda'].node:q.node,refs['Z'].node:z.node}
        aliases={key:value for key,value in aliases.items() if key!=value}
        return proof,aliases

    @source_precision
    def invert(self,request):
        self._validate_input(request);g=self.graph;c=self.ctx;f=request.forward_coordinates
        sign=1 if request.Z>0 else -1 if request.Z<0 else 0
        refs=self.inverse.inverse_functions(f['log_r'],request.physical_log_abs_z,f['log_tau'],f['theta'],sign)
        proof,aliases=self._identity(request,refs)
        reader=CorrelatedGraphBounds(self.locator.reader({self.inverse.delta_function.node:self.inverse.delta}),aliases)
        q=reader.at(refs['log_lambda']);z=reader.at(refs['Z']);lr=reader.at(refs['logR'])
        physical_lr=reader.at(refs['log_r']);theta=reader.at(refs['theta'])
        result=dict(input_kind='live_source_bound_correlated_original_physical_logs',
            exact_physical_input_record=self._input_record(request),coordinate_functions=refs,
            exact_original_inverse_identity=proof,exact_root_aliases=aliases,
            directed_inverse_mapping=dict(actual_log_lambda=q,Z=z,solver_status='exact_forward_identity_and_original_root_uniqueness',
                iterations=0,source_delta=self.inverse.delta,physical_viscosity=c.mpf(1),
                requested_log_tau=reader.at(refs['log_tau']),
                requested_log_abs_z=None if request.physical_log_abs_z is None else reader.at(request.physical_log_abs_z)),
            physical_log_radius_enclosure=physical_lr,theta_enclosure=theta,source_logR_enclosure=lr,
            strict_numerical_Z_interior_resolved=-1<ends(z)[0]<=ends(z)[1]<1,
            exact_root_node_retained_with_separate_equality_proof=True,
            common_original_scales_cancelled_before_numeric_coordinate_arithmetic=True,
            finite_offset_not_lost_beside_original_logC=True,
            numerical_absolute_radius_lambda_time_or_z_materialized=False,
            supported_input_family='physical forward logs with original log_tau=-source_logR+exact finite offset',
            arbitrary_physical_graph_input_or_global_field_installed=False,
            original_N_definition=self.original_N_definition,source_family=self.family_record,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        self._inverse_views[id(result)]=(result,self._fingerprint(result),request,dict(aliases))
        return result

    def _validated_reader(self,view):
        record=self._inverse_views.get(id(view))
        if record is None or record[0] is not view or record[1]!=self._fingerprint(view):
            raise ValueError('Unchanged live correlated inverse from this current owner required')
        request=record[2];proof,aliases=self._identity(request,view['coordinate_functions'])
        if proof!=view['exact_original_inverse_identity'] or aliases!=record[3] or aliases!=view['exact_root_aliases']:
            raise ValueError('Current original inverse equality proof changed')
        return CorrelatedGraphBounds(self.locator.reader({self.inverse.delta_function.node:self.inverse.delta}),aliases)

    @source_precision
    def locate(self,view):
        reader=self._validated_reader(view)
        return self.locator._locate(view['coordinate_functions']['logR'],reader,view['strict_numerical_Z_interior_resolved'])

    @source_precision
    def source(self,view,chart=None):
        reader=self._validated_reader(view);location=self.locate(view)
        request=self._inverse_views[id(view)][2]
        chart=request.chart if chart is None else chart
        if chart!=request.chart:raise ValueError('This source-bound physical family uses its original chart; boundary protocols remain separate')
        if chart not in location['contained_charts'] or not location['physical_Z_numeric_interior_resolved']:
            raise ValueError('Entire correlated inverse must be admitted to one current source chart')
        row=location['chart_rows'][chart];refs=view['coordinate_functions'];Z=reader.at(refs['Z'])
        provenance=dict(kind='live_proved_correlated_original_physical_inverse',
            exact_physical_input_record=view['exact_physical_input_record'],
            original_inverse_root=refs['log_lambda'].node,
            original_inverse_equality_proof=view['exact_original_inverse_identity'].node,
            exact_logR=refs['logR'].node)
        route=self.locator.maps[chart];g=self.graph
        difference=g.sub(refs['logR'],route['base'])
        if reader.polynomial(g.sub(difference,g.mul(route['jacobian'],g.constant(request.native)))):
            raise ValueError('Original exact affine native-coordinate identity required')
        native_proof=g.node('exact_original_correlated_native_inverse_identity',
            original_root_equality_proof=view['exact_original_inverse_identity'].node,
            native_inverse_function=row['exact_native_coordinate_function'].node,original_native_coordinate=g.constant(request.native).node,
            original_source_chart=chart,original_source_base=route['base'].node,
            original_source_Jacobian=route['jacobian'].node,
            identity='logR-base=J*native; J>0; exact affine inverse equals original native coordinate')
        # Both enclosures bound the same function after this exact equality
        # proof. Their intersection refines rounding/dependency error; it
        # does not select a midpoint, endpoint or replacement coordinate.
        current_lo,current_hi=ends(row['directed_native_coordinate'])
        exact_lo,exact_hi=ends(inverse.box(self.ctx,request.native))
        lo,hi=max(current_lo,exact_lo),min(current_hi,exact_hi)
        if lo>hi:raise ArithmeticError('Proved native coordinate enclosures disagree')
        native_box=self.ctx.mpf([lo,hi]);provenance['exact_affine_native_identity']=native_proof.node
        token=self.before._register(chart,row['exact_native_coordinate_function'],native_box,Z,refs['Z'],provenance)
        source=self.before._evaluate(token,Z,refs['Z'],provenance)
        self.call_trace.append(dict(chart=chart,actual_original_scale_correlated_physical_inverse_consumed=True,
            original_root_equality_proof=view['exact_original_inverse_identity'].node,
            source_provider_called=True,conditional_clipped_chart_intersection_not_consumed=True))
        result=dict(physical_inverse=view,native_location=location,actual_source_view=source,
            exact_native_inverse_identity=native_proof,
            source_native_coordinate_enclosure=native_box,native_enclosure_refined_by_exact_affine_identity=True,
            correlated_physical_input_source_delivery=True,physical_values_materialized=False,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        self._deliveries[id(result)]=(result,self._fingerprint(report(result)),token,view)
        return result

    def _validate_delivery(self,delivery):
        record=self._deliveries.get(id(delivery))
        if record is None or record[0] is not delivery or record[1]!=self._fingerprint(report(delivery)):
            raise ValueError('Unchanged live correlated source delivery from this current owner required')
        reader=self._validated_reader(record[3]);token=record[2]
        self.before._validate_token(token.chart,token)
        request=self._inverse_views[id(record[3])][2]
        proof=self.graph.nodes[delivery['exact_native_inverse_identity'].node]
        if (proof['native_inverse_function'],proof['original_native_coordinate'])!=(token.function.node,self.graph.constant(request.native).node):
            raise ValueError('Current native inverse identity changed')
        reader.aliases[token.function.node]=self.graph.constant(request.native).node
        return reader,token,request

    @source_precision
    def physical_rows(self,delivery):
        """Map the live source to original physical operators at this point.

        Exact point coordinates follow the proved inverse identities. The
        coefficient boxes still carry every directed source/rounding error.
        The existing exact-rational point guard is not weakened globally.
        """
        reader,token,request=self._validate_delivery(delivery)
        view=delivery['actual_source_view'];Zbox=reader.at(delivery['physical_inverse']['coordinate_functions']['Z'])
        proxy=box._Proxy(self.physical)
        owner=self
        def source_view(proxy,observed,Z):
            if observed is not view or inverse.exact_scalar(Z)!=request.Z:
                raise ValueError('This physical map requires its unchanged correlated source point')
            owner._validate_delivery(delivery)
            geometry=owner.before.box_radius.geometry(token.chart,token)
            if owner._fingerprint(geometry)!=owner._fingerprint(view['geometry']):
                raise ValueError('Current source token geometry required')
            for name in physical.SOURCE.values():
                rows=view['log_radius_mixed_rows'][name];base=view['original_factorized_values'][name]
                if len(rows)!=15:raise ValueError('Full ordinary similarity mixed4 grid required')
                for k in range(5):
                    for j in range(5-k):
                        row=rows['y%d_Z%d'%(k,j)]
                        powers=(tuple(map(Fraction,(0,2,2 if token.chart in box.pulse.PULSE else 0)))
                            if name=='pressure' and k else base.powers)
                        expected=owner.before.transport.factor(name,token.chart,token,geometry,powers,row.coefficients,k,j,False)
                        if type(row) is not physical.mixed.FactorizedMixedSourceRow or row.derivative!=(k,j) or \
                                row.coefficients.order!=0 or row.coefficients.ctx is not owner.ctx or \
                                (row.powers,row.source_units,row.log_scale_parts)!=(expected.powers,expected.source_units,expected.log_scale_parts):
                            raise ValueError('Original box source units and coefficient context required')
            return request.native,request.Z,Zbox
        def coordinates(proxy,geometry,Z,log_tau,theta):
            if inverse.exact_scalar(Z)!=request.Z or log_tau!=request.forward_coordinates['log_tau'] or \
                    inverse.exact_scalar(theta)!=request.theta or owner._fingerprint(geometry)!=owner._fingerprint(view['geometry']):
                raise ValueError('Proved original correlated physical coordinate functions required')
            return request.forward_coordinates
        proxy.source_view=MethodType(source_view,proxy);proxy.coordinates=MethodType(coordinates,proxy)
        proxy.bind('map_source_view',physical.CurrentOriginalRpPhysicalSourceMap.map_source_view)
        rows=proxy.map_source_view(view,request.Z,request.forward_coordinates['log_tau'],request.theta)
        rows.update(correlated_original_physical_input_identity=delivery['physical_inverse']['exact_original_inverse_identity'],
            correlated_original_native_inverse_identity=delivery['exact_native_inverse_identity'],
            actual_correlated_physical_point_source_consumed=True,
            source_coefficient_and_coordinate_enclosures_retained=True,
            unrestricted_physical_input_and_materialization_still_open=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        return rows


def report(result):
    if 'actual_source_view' in result:
        return locator.report({**{key:value for key,value in result.items() if key!='actual_source_view'},
            'actual_source_view':box.report(result['actual_source_view'])})
    return locator.report(result)


CASES={'pulse':('pulse_exit','21/2','.521','0','7/10'),
    'flatten':('flatten','50','.521','1/3','7/10'),
    'angular':('outer_angular','-2','.521','-1/5','7/10'),
    'waiting':('waiting','1/2','.521','2/7','7/10'),
    'collar':('heat_collar','3/2','.521','0','7/10'),
    'exterior':('heat_exterior','4','.521','0','7/10')}
INVERSE_ONLY={'negative_axial':('heat_exterior','9/2','-.521','0','7/10'),
    'zero_axial':('heat_exterior','9/2','0','0','7/10')}


@source_precision
def run(before=None):
    began=time.monotonic();owner=CurrentOriginalRpCorrelatedPhysicalSource(before,require_checked=False)
    requests={};inverses={};delivered={};mapped={}
    for name,args in {**CASES,**INVERSE_ONLY}.items():
        requests[name]=owner.physical_input(*args);inverses[name]=owner.invert(requests[name])
        if name in CASES:
            delivered[name]=owner.source(inverses[name])
            mapped[name]=owner.physical_rows(delivered[name])
            print('Actual original-scale correlated physical input delivered:',name,flush=True)
    result=dict(source_family=owner.family_record,actual_current_graph=owner.assert_graph(),
        actual_correlated_physical_inverses={key:report(value) for key,value in inverses.items()},
        actual_correlated_source_deliveries={key:report(value) for key,value in delivered.items()},
        actual_correlated_physical_operator_rows={key:locator.report(physical.report(value)) for key,value in mapped.items()},
        exact_current_correlated_physical_graph=owner.graph.nodes,
        actual_correlated_physical_source_call_trace=owner.call_trace,
        restricted_source_bound_input_family=True,unrestricted_physical_input_and_materializer_still_open=True,
        **dict.fromkeys(GATES+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(report(result),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,requests,inverses,delivered,mapped
