"""Role-bound original density/Rc-amplitude range oracle in the live context.

Original function ranges remain ScaledEnclosure objects. This companion to
the scalar oracle API does not cast formal scales to ordinary scalar values,
read saved coefficient midpoints, or admit a full original numerical route.
"""
from dataclasses import dataclass
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_Rc_function_transport as transport
import lei_ren_part1_paper_compliant_current_native_density_C1_local_integrals as density

HERE, PREFIX, sha = transport.HERE, transport.PREFIX, transport.sha
NAME = PREFIX+'current_native_Rc_factored_source_oracle.json'
RECEIPT = PREFIX+'current_native_Rc_factored_source_oracle_check.json'
GATE = 'current_original_role_bound_Rc_factored_source_oracle_selected_queries_executed'
current = transport.current
packets = current.packets


class RoleBoundNativeRcFunctions:
    """Add explicit graph namespace/role to the immutable accepted graph."""
    def __init__(self, owner):
        if type(owner) is not transport.NativeRcFunctionTransport:
            raise TypeError('Same accepted original function-transport owner required')
        checked = json.loads((HERE/transport.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(transport.GATE) or checked['source_family'] != owner.family:
            raise ValueError('Checked original function transport required')
        self.owner, self.family, self.service = owner, owner.family, owner.service
        self.hashes = {**checked['input_hashes'], transport.RECEIPT:sha(transport.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)}

    def build(self):
        built = self.owner.build()
        g = built['graph'];roles = {}
        def annotate(index, namespace, role):
            row = g.nodes[index]
            if row['operation'] != 'original_function_graph':
                raise ValueError('Original function source reference required')
            if index in roles and roles[index] != (namespace,role):
                raise ValueError('Ambiguous shared function source node')
            row['source_graph_namespace'], row['function_role'] = namespace, role
            roles[index] = (namespace,role)
        for cell in built['cells']:
            if cell['source_flat_exact_zero']:
                continue
            for key,pair in cell['contributions'].items():
                for order in ('value','Z'):
                    integral = g.nodes[pair[order]]
                    integrand = g.nodes[integral['integrand']]
                    indices = [i for i in integrand['arguments'] if g.nodes[i]['operation']=='original_function_graph']
                    if len(indices)!=1: raise ValueError('One full signed density function per integral required')
                    annotate(indices[0], 'function_graph_nodes', 'density_'+key+('_Z' if order=='Z' else '_C0'))
        annotate(built['amplitude'].value.node,'original_signed_input_graph.jet_expression_dag','Rc_E_C0')
        annotate(built['amplitude'].Z.node,'original_signed_input_graph.jet_expression_dag','Rc_E_Z')
        if len(roles)!=212: raise ValueError('All 210 density and two Rc source references required')
        # Annotations change interning keys. Keep the graph usable for later
        # control-map extensions without retaining stale untagged keys.
        g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(g.nodes)}
        built['source_roles']=roles
        return built


@dataclass
class FactoredSourceFrame:
    chart: str
    N: int
    values: object
    record: dict


class NativeFactoredSourceOracle:
    mode = 'original_factored_range'
    def __init__(self, role_owner, built=None):
        if type(role_owner) is not RoleBoundNativeRcFunctions:
            raise TypeError('Explicit original role-bound functions required')
        self.role_owner = role_owner
        self.target_owner = role_owner.owner.owner
        self.owner = self.target_owner.transfer.owner.owner
        if type(self.owner) is not density.NativeDensityC1LocalIntegrals:
            raise TypeError('Existing accepted original density/first-jet owner required')
        if self.owner.ctx is not self.target_owner.ctx or self.owner.family!=role_owner.family:
            raise ValueError('Same live arithmetic context and original family required')
        self.ctx, self.source_family, self.service = self.owner.ctx, self.owner.family, self.owner.service
        self.built = role_owner.build() if built is None else built
        self.source_graph_sha256 = self.built['source_graph_sha256']
        receipt = json.loads((HERE/density.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(density.GATE) or receipt['source_family']!=self.source_family:
            raise ValueError('Accepted original density C0/Z backend required')
        self.hashes={**role_owner.hashes, **receipt['input_hashes'], density.RECEIPT:sha(density.RECEIPT)}
        self.service.bind_hashes(self.hashes)

    def require_row(self,row):
        if row.get('operation')!='original_function_graph' or row.get('graph_sha256')!=self.source_graph_sha256:
            raise ValueError('Same accepted original source graph required')
        chart,role,namespace = row['chart'],row.get('function_role'),row.get('source_graph_namespace')
        view=self.role_owner.owner.views[chart]
        if role in ('Rc_E_C0','Rc_E_Z'):
            if chart!='O3_power' or namespace!='original_signed_input_graph.jet_expression_dag':
                raise ValueError('Exact original signed Rc E namespace required')
            key='y0_Z0' if role=='Rc_E_C0' else 'y0_Z1'
            expected=view['original_signed_input_graph']['jet_expression_dag']['roots']['E'][key]
        elif role and role.startswith('density_'):
            if namespace!='function_graph_nodes': raise ValueError('Original full loop-density namespace required')
            key=role.split('_')[1]
            if key not in transport.RATES or role not in ('density_'+key+'_C0','density_'+key+'_Z'):
                raise ValueError('Unknown signed original density role')
            expected=(view['five_signed_increment_rate_roots'][key] if role.endswith('_C0') else
                view['five_signed_increment_rate_first_derivatives']['Z'][key])
        else:
            raise ValueError('Explicit original function role required; integer node ID is insufficient')
        if row['source_node']!=expected: raise ValueError('Role does not match original defining source node')
        return role

    @density.native.inlet.source_precision
    def density_frame(self,*,chart,Z,coordinate,N):
        N=density.density.candidate_integer(N)
        if N<transport.current.MIN_N: raise ValueError('Original function route requires integer N>=160')
        got=self.owner.spatial_query(chart,Z,coordinate,N)
        enclosed=bool(got['cells']) and all(cell['values'] is not None for cell in got['cells'])
        values=None
        if enclosed:
            union=density.local.same_source_union
            values=dict(C0={key:union([cell['values']['kernels'][key] for cell in got['cells']]) for key in transport.RATES},
                Z={key:union([cell['values']['Z_derivatives'][key] for cell in got['cells']]) for key in transport.RATES})
            if any(value.ctx is not self.ctx for rows in values.values() for value in rows.values()):
                raise ValueError('Original factored ranges changed arithmetic context')
        packet=got['source']['packet']
        background=current.history.packet_history_functions(packet,self.target_owner.coordinates,
            self.target_owner.transfer.owner.signed_owner)
        record=dict(status='enclosed' if enclosed else 'requires_source_or_phase_refinement',
            chart=chart,source_family=self.source_family,candidate_N=N,
            actual_original_spatial_source=got['record'],
            derived_original_global_phase=got['geometry']['record'],
            separate_original_P0=background['P0'].record(),separate_original_P0_Z=background['P0_Z'].record(),
            signed_density_C0_ranges=None if values is None else {key:v.record() for key,v in values['C0'].items()},
            signed_density_Z_ranges=None if values is None else {key:v.record() for key,v in values['Z'].items()},
            arithmetic='same live ScaledEnclosure context, formal source/radius factors and ledger',
            saved_cover_midpoints_or_plain_scalar_cast_used=False,
            phase_derived_from_original_radius_not_supplied_as_free_angle=True,
            all_phase_cells_unioned_with_source_functions_retained=True,
            original_nonzero_V_and_E_and_all_cross_terms_retained=True,
            phase_solver_backend='accepted NativePhaseFirstJets, original bits60 solver; unresolved status retained',
            numerical_integral_or_full_route_or_controls_admitted=False,
            **dict.fromkeys(packets.OPEN,False))
        return FactoredSourceFrame(chart,N,values,record)

    @density.native.inlet.source_precision
    def amplitude_frame(self,*,Z,N):
        N=density.density.candidate_integer(N)
        if N<transport.current.MIN_N: raise ValueError('Original function route requires integer N>=160')
        binder=self.target_owner.transfer.geometry.binder
        geometry=binder.query('O3_power',Z,{'original_power_offset':'2'},N)
        got=self.target_owner.q_owner.query('O3_power',Z,geometry['raw']['coordinate'])
        roots=got['source']['roots'];values={key:roots['E'][order] for key,order in (('Rc_E_C0',(0,0)),('Rc_E_Z',(0,1)))}
        if any(value.ctx is not self.ctx for value in values.values()):
            raise ValueError('Original Rc amplitude context changed')
        packet=got['source']['packet']
        background=current.history.packet_history_functions(packet,self.target_owner.coordinates,
            self.target_owner.transfer.owner.signed_owner)
        record=dict(status='enclosed',source_family=self.source_family,chart='O3_power',candidate_N=N,
            original_endpoint='original_power_offset=2; coordinate=2/same Tw',
            actual_original_spatial_geometry=geometry['record'],
            original_source=got['record'],ranges={key:value.record() for key,value in values.items()},
            separate_original_P0=background['P0'].record(),separate_original_P0_Z=background['P0_Z'].record(),
            Rc_amplitude_not_candidate_modulated_E=True,phase_independent_amplitude_from_original_signed_roots=True,
            source_coefficients_or_range_midpoints_not_selected=True,
            numerical_integral_or_full_route_or_controls_admitted=False,**dict.fromkeys(packets.OPEN,False))
        return FactoredSourceFrame('O3_power',N,values,record)

    def dispatch_function_range(self,row,frame):
        role=self.require_row(row)
        if type(frame) is not FactoredSourceFrame or frame.chart!=row['chart']:
            raise ValueError('Matching live original chart frame required')
        if frame.values is None: raise ArithmeticError('Refine unresolved original source/phase; no fallback zero')
        if role in ('Rc_E_C0','Rc_E_Z'):
            if role not in frame.values: raise ValueError('Original Rc amplitude frame required')
            result=frame.values[role]
        else:
            key=role.split('_')[1];order='Z' if role.endswith('_Z') else 'C0'
            if order not in frame.values: raise ValueError('Original density frame required')
            result=frame.values[order][key]
        if result.ctx is not self.ctx: raise ValueError('Factored function range must stay in original context')
        return result


@density.native.inlet.source_precision
def run(role_owner):
    began=time.monotonic();built=role_owner.build();oracle=NativeFactoredSourceOracle(role_owner,built)
    points=dict(inner_reference='.1337',O2_slope='.1337',O2_buffer='5.337',
        O3_slope_mu='.537',O3_power={'original_power_offset':'.537'})
    records={};dispatches=0;enclosed=0
    rows=[row for row in built['graph'].nodes if row['operation']=='original_function_graph']
    byrole={(row['chart'],row['function_role']):row for row in rows}
    # O3 power is exactly quiet, so it has no density graph references.
    # Its explicit frame still checks that original E/V/P0 are preserved.
    for N in (1024,2048):
        for chart,coordinate in points.items():
            frame=oracle.density_frame(chart=chart,Z=('.5','.5'),coordinate=coordinate,N=N)
            records[chart+'_N'+str(N)]=frame.record
            if frame.values is not None:
                enclosed+=1
                for key in transport.RATES:
                    for order in ('C0','Z'):
                        row=byrole.get((chart,'density_'+key+'_'+order))
                        if row is not None:
                            oracle.dispatch_function_range(row,frame);dispatches+=1
            print('Original factored source oracle:',chart,N,frame.record['status'],flush=True)
    amplitude=oracle.amplitude_frame(Z=('.49','.51'),N=2048)
    for role in ('Rc_E_C0','Rc_E_Z'):
        oracle.dispatch_function_range(byrole['O3_power',role],amplitude);dispatches+=1
    if enclosed==0: raise ArithmeticError('No actual original source query resolved')
    result=dict(source_family=oracle.source_family,**{GATE:True},
        original_function_role_annotations={str(i):dict(namespace=namespace,role=role) for i,(namespace,role) in built['source_roles'].items()},
        original_factored_source_query_records=records,original_Rc_amplitude_query_record=amplitude.record,
        declared_density_queries=10,enclosed_density_queries=enclosed,
        actual_role_dispatched_factored_ranges=dispatches,candidate_N_queries=[1024,2048],
        numerical_source_ranges_from_live_original_queries=True,
        original_parameter_or_function_point_scalar_values_installed=False,
        all_17_chart_original_range_oracle_admitted=False,
        full_factored_function_graph_evaluator_installed=False,
        actual_original_numerical_integrals_evaluated=False,actual_five_controls_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=oracle.hashes,
        scope='Role-bound original function graph dispatch to selected actual live density C0/Z ranges at two distinct candidate N values and original Rc E/E_Z, preserving factored units/context/global phase/P0. Selected source-range oracle only; scalar values, all-chart/full-graph integration, controls/global N/closure/recursion remain open.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    return result
