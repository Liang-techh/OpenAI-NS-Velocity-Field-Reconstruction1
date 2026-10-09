"""Defining reference/O2 point leaves for the original shared-N source graph.

Actual point coefficients, quadrature errors and inverse/slow-Z objects feed
the accepted function DAG. Ordinary interval callbacks are available only
when the resulting source enclosure is materializable with directed tails.
Other charts, parameter/integral callbacks and full control evaluation remain
explicitly unavailable. Formal point values are never cast or hidden.
"""
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_all_chart_point_phase as precise
import lei_ren_part1_paper_compliant_current_original_reference_point_oracle as reference
import lei_ren_part1_paper_compliant_current_original_O2_signed_densities as density
import lei_ren_part1_paper_compliant_current_native_local_signed_integrals as unions

base=density.base;slow=density.slow;prior=base.prior;ep=density.ep
transport=precise.functions;HERE,PREFIX,sha=precise.HERE,precise.PREFIX,precise.sha
NAME=PREFIX+'current_original_point_source_leaves.json'
RECEIPT=PREFIX+'current_original_point_source_leaves_check.json'
GATE='actual_original_reference_and_O2_graph_bound_C0_Z_point_source_callbacks_installed'
SUPPORTED=('Rh_reference','O2_slope')


def digest_rows(rows):
    import hashlib
    return hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':')).encode()).hexdigest()


@dataclass(frozen=True)
class OriginalPointSourceFrame:
    chart: str
    coordinate: object
    Z: object
    N: int
    ctx: object
    values: dict
    record: dict
    query: dict
    pieces: tuple


class OriginalPointSourceLeaves:
    # A partial provider is not admitted by the full scalar control evaluator.
    mode='original_point_source_partial'
    def __init__(self,phase_service=None):
        if phase_service is None:
            bridge,_=precise.phase.native.inlet.native_bridge_owner()
            with precise.phase.native.inlet.CheckedSourceRuntime():
                phase_service=precise.OriginalAllChartPointPhase(precise.phase.NativeSpatialPhase(
                    precise.phase.native.NativeGenericSourcePackets(bridge)))
        if type(phase_service) is not precise.OriginalAllChartPointPhase:
            raise TypeError('Genuine accepted original point phase service required')
        self.phase=phase_service;self.family=self.source_family=phase_service.family
        self.source_graph_sha256=phase_service.source_graph_sha256;self.hashes={}
        self._frames={};self._issued={};self.materialization_records=[]
        for module in (precise,reference,density,transport):self.bind_receipt(module)
        self.O2=density.OriginalO2SignedDensities();self.reference=reference.OriginalReferencePointOracle()
        if self.O2.family!=self.family or self.reference.family!=self.family:
            raise ValueError('One original pressure/source family required')
        for closure in (phase_service.hashes,self.O2.hashes,self.reference.hashes):self.bind_hashes(closure)
        self.views=json.loads(gzip.decompress((HERE/transport.sources.VIEWS).read_bytes()))
        if any(self.views[chart]['source_family']!=self.family for chart in SUPPORTED):
            raise ValueError('Original point graph chart family differs')
        chosen=ep(self.phase.binder.fixed['logC'])[0]._mpf_
        for owner in (self.O2.owner,self.reference.owner):
            if tuple(owner.inputs.frame.selected_logCstar_mpf_tuple)!=tuple(chosen):
                raise ValueError('Point coefficients and actual radius use different selected logCstar')
        self.built=self.exact_transport_rows()
        self.graph_digest=digest_rows(self.built['graph'].nodes)
        self.rows={(row['chart'],row['function_role']):row for row in self.built['graph'].nodes
            if row['operation']=='original_function_graph'}
        self.bind_hashes({name:sha(name) for name in (transport.NAME,transport.sources.VIEWS,
            Path(unions.__file__).name,Path(__file__).name)})

    def bind_hashes(self,closure):
        for name,digest in closure.items():
            if sha(name)!=digest:raise ValueError('Original point dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original point closures disagree: '+name)
            self.hashes[name]=digest

    def bind_receipt(self,module):
        receipt=json.loads((HERE/module.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(module.GATE) or receipt['source_family']!=self.family:
            raise ValueError('Accepted same-source point prerequisite required: '+module.RECEIPT)
        self.bind_hashes({**receipt['input_hashes'],module.RECEIPT:sha(module.RECEIPT)})

    def exact_transport_rows(self):
        """Reuse checked exact expression syntax only, never numerical covers."""
        saved=json.loads((HERE/transport.NAME).read_bytes())
        if not saved[transport.GATE] or saved['source_family']!=self.family:raise ValueError('Original function graph required')
        g=transport.FunctionTransportGraph()
        for index,row in enumerate(saved['function_graph_nodes']):
            if g.node(**row).node!=index:raise ValueError('Accepted transport expression interning changed')
        roles={}
        def annotate(index,namespace,role):
            row=g.nodes[index]
            if row['operation']!='original_function_graph':raise ValueError('Original source row required')
            if index in roles and roles[index]!=(namespace,role):raise ValueError('Ambiguous original source role')
            row['source_graph_namespace']=namespace;row['function_role']=role;roles[index]=(namespace,role)
        for cell in saved['exact_original_cells']:
            if cell['source_flat_exact_zero']:continue
            for key,pair in cell['contributions'].items():
                for order in ('value','Z'):
                    integral=g.nodes[pair[order]];integrand=g.nodes[integral['integrand']]
                    nodes=[i for i in integrand['arguments'] if g.nodes[i]['operation']=='original_function_graph']
                    if len(nodes)!=1:raise ValueError('One original signed density source per integral required')
                    annotate(nodes[0],'function_graph_nodes','density_'+key+('_Z' if order=='Z' else '_C0'))
        for order,role in (('value','Rc_E_C0'),('Z','Rc_E_Z')):
            annotate(saved['exact_Rc_amplitude_roots'][order],'original_signed_input_graph.jet_expression_dag',role)
        if len(roles)!=212:raise ValueError('All212 original density/amplitude source roles required')
        g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(g.nodes)}
        Nnodes=[i for i,row in enumerate(g.nodes) if row['operation']=='shared_positive_integer']
        if len(Nnodes)!=1:raise ValueError('One shared original N required')
        parameters={row['name']:transport.FunctionRef(g,i) for i,row in enumerate(g.nodes)
            if row['operation']=='original_source_parameter'}
        return dict(graph=g,N=transport.FunctionRef(g,Nnodes[0]),parameters=parameters,
            source_family=self.family,source_graph_sha256=self.source_graph_sha256,source_roles=roles)

    def require_row(self,row):
        g=self.built['graph']
        if digest_rows(g.nodes)!=self.graph_digest or not any(row is item for item in g.nodes):
            raise ValueError('Unchanged issued original source row required')
        if row.get('operation')!='original_function_graph' or row.get('graph_sha256')!=self.source_graph_sha256:
            raise ValueError('Same defining original source graph required')
        chart=row['chart'];role=row.get('function_role')
        if chart not in SUPPORTED:raise NotImplementedError('Defining point leaves not installed on '+chart)
        if row.get('source_graph_namespace')!='function_graph_nodes' or not role or not role.startswith('density_'):
            raise ValueError('Original full signed density namespace required')
        key=role.split('_')[1]
        if key not in transport.RATES or role not in ('density_'+key+'_C0','density_'+key+'_Z'):
            raise ValueError('Original C0/Z density role required')
        view=self.views[chart]
        expected=view['five_signed_increment_rate_roots'][key] if role.endswith('_C0') else view['five_signed_increment_rate_first_derivatives']['Z'][key]
        if row['source_node']!=expected or row['shared_N']!=self.built['N'].node:
            raise ValueError('Original root/role/shared N binding differs')
        return key,'C0' if role.endswith('_C0') else 'Z'

    def frame(self,*,chart,coordinate,Z,N,bits=80):
        if chart not in SUPPORTED:raise NotImplementedError('Defining point chart not installed: '+chart)
        q=precise.rational(coordinate);z=base.point.pressure.exact_Z(Z)
        if type(N) is not int or N<160 or N.bit_length()>4096:raise ValueError('Explicit original graph candidate integer N>=160 required')
        if type(bits) is not int or not 4<=bits<=256:raise ValueError('Inverse bits in[4,256] required')
        if chart=='Rh_reference':reference.reference_coordinate(q)
        else:base.point.source.exact_coordinate(q)
        cache_key=(chart,q,z,N,bits)
        if cache_key in self._frames:return self._frames[cache_key]
        probe=self.rows[chart,'density_m_C0'];self.require_row(probe)
        geometry=self.phase.phase_for_source(probe,self.built,coordinate=str(q),N=N)
        owner=self.reference.owner if chart=='Rh_reference' else self.O2.owner;c=owner.ctx
        reused=(q,z) in owner.inputs.point_cache
        with mp.workdps(c.dps+40):
            query=reference.reference_query(owner,y=q,Z=z) if chart=='Rh_reference' else owner.query(y=q,Z=z)
            kernel=query['kernel'];pieces=[];view=self.views[chart]
            a=base.conditioned.bounded_value(query['roots']['a'][(0,0)])
            if not 0<=ep(a)[0]<=ep(a)[1]<=2:raise ValueError('Original point modulation |A/N|<=1 proof requires0<a<=2')
            for box in geometry['actual_phase_directed_boxes']:
                inverse=kernel.evaluate(c.mpf(box),bits=bits)
                if inverse['status']!='enclosed':raise ArithmeticError('Original point source inverse needs refinement')
                selected=inverse['selected_inverse']
                primitives=kernel.primitives(selected['coordinate_interval'],selected['chart'])
                jets,derivative=slow.slow_values(kernel,query['roots'],selected['coordinate_interval'],selected['chart'])
                interpreter=density.BoundDensityGraph(view,query,primitives,jets,N)
                values=interpreter.outputs()
                inverse.update(free_phase_parameter_not_spatial_phase=False,original_common_N_and_radius_phase_bound=True)
                contract=interpreter.record()
                contract.update(point_chart=chart,actual_point_a_enclosure=a,
                    modulation_argument_proof='actual defining Rh_reference a=4/5 or O2 a=.8+1.2*sigma in[.8,2]; same A=a*(phi-psi/(2*pi))/2',
                    B_over_Pstar_already_normalized_not_multiplied_again=True)
                pieces.append(dict(values=values,primitives=primitives,jets=jets,inverse=inverse,
                    derivative=derivative,interpreter=interpreter,execution_contract=contract))
            union=lambda values:unions.same_source_union(values)
            values={order:{key:union([piece['values'][name][key] for piece in pieces]) for key in transport.RATES}
                for order,name in (('C0','densities'),('Z','density_Z'))}
            if any(value.ctx is not c or value.ledger is not query['ledger'] or value.scale.bases is not kernel.q.scale.bases
                for row in values.values() for value in row.values()):raise ValueError('One source point context/basis/ledger required')
            record=dict(source_family=self.family,chart=chart,original_coordinate_exact=str(q),original_Z_exact=str(z),explicit_candidate_N=N,
                source_graph_sha256=self.source_graph_sha256,actual_original_point_phase=geometry,
                actual_defining_point_inputs_and_errors=base.point.record_query(query['point']),
                defining_point_coefficient_cache_reused=reused,source_factor_basis=query['basis_contract'],
                native_point_inverse_and_derivative_records=[dict(inverse=p['inverse'],derivative=p['derivative'],graph=p['execution_contract']) for p in pieces],
                actual_original_five_C0_Z_point_source_enclosures={order:{key:value.record() for key,value in row.items()} for order,row in values.items()},
                numerical_arithmetic_ledger=query['ledger'],original_source_quadrature_errors_and_pressure_tail_errors_retained=True,
                actual_graph_nodes_executed_not_cover_values_selected=True,phase_union_is_one_point_error_cover_not_duplicate_source_mass=True,
                ordinary_Z_at_fixed_true_radius_phase=True,original_P0_P0_Z_not_modified=True,
                per_log_radius_density_no_extra_Jacobian_applied=True,point_enclosure_not_a_continuous_integral=True,
                full_17_chart_source_or_integral_oracle_installed=False,actual_five_controls_installed=False,current_whole_N_selected=False)
            frame=OriginalPointSourceFrame(chart,q,z,N,c,values,record,query,tuple(pieces))
        self._frames[cache_key]=frame;self._issued[id(frame)]=frame;return frame

    def dispatch(self,row,frame):
        key,order=self.require_row(row)
        if type(frame) is not OriginalPointSourceFrame or self._issued.get(id(frame)) is not frame or frame.chart!=row['chart']:
            raise ValueError('Issued matching original point frame required')
        self.phase.phase_for_source(row,self.built,coordinate=str(frame.coordinate),N=frame.N)
        return frame.values[order][key]

    def source_factored(self,row,*,coordinate,Z,N,bits=80):
        self.require_row(row)
        return self.dispatch(row,self.frame(chart=row['chart'],coordinate=coordinate,Z=Z,N=N,bits=bits))

    def source(self,row,*,coordinate,Z,N,phase=None,bits=80):
        """A real ordinary interval, only after certified point/tail arithmetic."""
        if phase is not None:raise ValueError('Point source derives its true phase; external free/naive phase is forbidden')
        frame=self.frame(chart=row['chart'],coordinate=coordinate,Z=Z,N=N,bits=bits)
        value=self.dispatch(row,frame)
        with mp.workdps(frame.ctx.dps+40):
            result=base.conditioned.bounded_value(value)
        if not hasattr(result,'_mpi_') or hasattr(result,'scale'):raise TypeError('Directed ordinary source point interval required')
        self.materialization_records.append(dict(chart=frame.chart,coordinate=str(frame.coordinate),Z=str(frame.Z),N=N,
            role=row['function_role'],source_node=row['source_node'],ordinary_interval=result,
            original_point_factor_enclosure=value.record(),directed_small_exponential_tail_enclosure_retained=True,
            field_value_not_selected_from_range_cap=True))
        return result

    def parameter(self,name):raise NotImplementedError('Full original parameter evaluator not installed by this point-only provider')
    def integrate(self,*args,**kwargs):raise NotImplementedError('Certified original definite-integral callback still required')


@precise.phase.native.inlet.source_precision
def run(*,return_live=False):
    began=time.monotonic();owner=OriginalPointSourceLeaves();records=[];frames=[];dispatches=0
    samples=(('Rh_reference','-2.337','-.37',257),('Rh_reference','-2.337','0',257),
        ('Rh_reference','-1.713','.37',257),('Rh_reference','0','.37',257),
        ('O2_slope','0','.37',257),('O2_slope','.53','.37',257),
        ('O2_slope','.731','-.37',257),('O2_slope','.731','0',160),('O2_slope','1','.37',512))
    for chart,coordinate,Z,N in samples:
        frame=owner.frame(chart=chart,coordinate=coordinate,Z=Z,N=N);records.append(frame.record);frames.append(frame)
        for key in transport.RATES:
            for order in ('C0','Z'):
                row=owner.rows[chart,'density_'+key+'_'+order];assert owner.dispatch(row,frame) is frame.values[order][key];dispatches+=1
        print('Original graph point source leaves:',chart,coordinate,Z,N,flush=True)
    ordinary=0;formal=[]
    for frame in frames:
        for key in transport.RATES:
            for order in ('C0','Z'):
                row=owner.rows[frame.chart,'density_'+key+'_'+order]
                try:owner.source(row,coordinate=str(frame.coordinate),Z=str(frame.Z),N=frame.N);ordinary+=1
                except ArithmeticError:formal.append(dict(chart=frame.chart,coordinate=str(frame.coordinate),Z=str(frame.Z),N=frame.N,role=row['function_role'],
                    original_point_source_still_retained_as_factored_value=True,not_replaced_with_zero_or_finite_cap=True))
    result=dict(**{GATE:True},source_family=owner.family,source_graph_sha256=owner.source_graph_sha256,
        exact_original_transport_graph_sha256=owner.graph_digest,mode=owner.mode,supported_defining_point_charts=list(SUPPORTED),
        actual_original_point_source_frames=records,original_role_bound_point_dispatches=dispatches,
        actual_ordinary_interval_point_callbacks=owner.materialization_records,ordinary_interval_callbacks=ordinary,
        unresolved_ordinary_materializations_with_factored_source_retained=formal,
        full_17_chart_source_or_integral_oracle_installed=False,full_scalar_control_evaluator_compatibility_installed=False,
        actual_original_numerical_integrals_evaluated=False,actual_five_controls_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(precise.phase.packets.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual reference/O2 defining point C0/Z density callbacks bound to original graph roots and precise radius phase. Directed ordinary interval results only where bounded; huge derivatives remain factored. All-chart leaves, parameter/integral callbacks, actual controls/global N, recursion and corrected NS remain open.')
    (HERE/NAME).write_bytes(json.dumps(precise.encode(result),indent=2).encode()+b'\n')
    return (result,owner,frames) if return_live else result


if __name__=='__main__':run()
