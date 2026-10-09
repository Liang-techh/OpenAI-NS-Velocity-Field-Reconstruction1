"""Genuine original O2 five-route C0/Z graph integral callbacks.

Accepted defining source-cell ranges retain their quadrature and pressure
errors. New exact shared-N phase covers replace the archived N7 phases.
The actual signed density DAG is executed before branch union and positive
own-rate integration. Fixed nonzero Z, partial evaluator contract only.
"""
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_reference_integral_callback as reference
import lei_ren_part1_paper_compliant_current_original_O2_density_Z_integrals as cells

leaves=reference.leaves;precise=leaves.precise;transport=leaves.transport
positive=cells.current;HERE,PREFIX,sha=leaves.HERE,leaves.PREFIX,leaves.sha;ep=leaves.ep
NAME=PREFIX+'current_original_O2_integral_callback.json.gz'
RECEIPT=PREFIX+'current_original_O2_integral_callback_check.json'
GATE='actual_original_O2_five_route_graph_bound_C0_Z_integral_callbacks_installed'


@dataclass(frozen=True)
class OriginalO2IntegralFrame:
    Z: object
    N: int
    count: int
    bits: int
    values: dict
    whole_values: dict
    record: dict


class OriginalO2IntegralCallback:
    mode='original_O2_integral_partial'
    def __init__(self,provider=None):
        self.provider=leaves.OriginalPointSourceLeaves() if provider is None else provider
        if type(self.provider) is not leaves.OriginalPointSourceLeaves:
            raise TypeError('Genuine issued point-source provider required')
        self.family=self.source_family=self.provider.family
        self.source_graph_sha256=self.provider.source_graph_sha256
        self.hashes=dict(self.provider.hashes);self.cache={};self.issued={}
        for module in (leaves,reference,cells):
            accepted=json.loads((HERE/module.RECEIPT).read_bytes())
            if not accepted.get('all_passed') or not accepted.get(module.GATE) or accepted['source_family']!=self.family:
                raise ValueError('Accepted same-source integral prerequisite required')
            self.bind_hashes({**accepted['input_hashes'],module.RECEIPT:sha(module.RECEIPT)})
        self.cells=cells.OriginalO2DensityZIntegrals();self.c=self.cells.c
        self.bind_hashes(self.cells.hashes)
        frame=self.cells.parent.owner.inputs.frame;current=self.provider.O2.owner.inputs.frame
        if self.cells.family!=self.family or frame.selected_logCstar_mpf_tuple!=current.selected_logCstar_mpf_tuple or any(
            frame.definitions[key]!=current.definitions[key] for key in ('logPstar','log_delta','logRref','L')):
            raise ValueError('Same defining O2 parameter/pressure/source family required')
        self.view=self.provider.views['O2_slope']
        if self.view!=self.cells.parent.owner.scales.graph:
            raise ValueError('Identical original full signed O2 density DAG required')
        saved=json.loads((HERE/transport.NAME).read_bytes())
        self.route=[cell for cell in saved['exact_original_cells'] if cell['chart']=='O2_slope']
        if len(self.route)!=5:raise ValueError('Five unchanged original O2 route cells required')
        self.roles={};self.rows={};self.node_ids={};self.sources={};self.bounds=[];self.contracts=[]
        g=self.provider.built['graph']
        for index,cell in enumerate(self.route):
            left,_=reference.symbolic_node(self.provider,cell['lower']);right,_=reference.symbolic_node(self.provider,cell['upper'])
            if not (left.is_Rational and right.is_Rational and 0<=left<right<=1):
                raise ValueError('Exact original slope subwindow required')
            if index and self.bounds[-1][1]!=left:raise ValueError('Original O2 route gap')
            self.bounds.append((left,right))
            for key,pair in cell['contributions'].items():
                for jet,label in (('C0','value'),('Z','Z')):
                    node=pair[label];row=g.nodes[node];integrand=g.nodes[row['integrand']]
                    source=[g.nodes[i] for i in integrand['arguments'] if g.nodes[i]['operation']=='original_function_graph']
                    if len(source)!=1:raise ValueError('One original unweighted density per integral required')
                    self.roles[node]=(index,key,jet);self.rows[index,key,jet]=row;self.node_ids[index,key,jet]=node;self.sources[node]=source[0]
                    self.contracts.append(self.require_integral(row))
        if self.bounds[0][0]!=0 or self.bounds[-1][1]!=1:raise ValueError('Complete original O2 window[0,1] required')
        x=self.provider.phase.x;offset=self.provider.phase.maps['O2_slope']
        if sy.diff(offset,x)!=1 or sy.simplify(offset-self.provider.phase.maps['Rh_reference'])!=0:
            raise ValueError('Original O2 affine radius/shared phase origin required')
        self.parameter_phase_binding=dict(passed=True,same_selected_logCstar_and_original_pressure_parameters=True,
            exact_O2_and_reference_radius_offset_formulas_equal=True,original_coordinate_Jacobian=1,
            same_precise_original_source_phase_recipe=True,original_r_minus_and_microscopic_origin_retained=True)
        self.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})

    def bind_hashes(self,closure):
        for name,digest in closure.items():
            if sha(name)!=digest:raise ValueError('Original O2 integral dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original O2 closures disagree: '+name)
            self.hashes[name]=digest

    def require_integral(self,row):
        g=self.provider.built['graph']
        if leaves.digest_rows(g.nodes)!=self.provider.graph_digest:raise ValueError('Unchanged original graph required')
        issued=[i for i,item in enumerate(g.nodes) if row is item]
        if len(issued)!=1 or issued[0] not in self.roles:raise ValueError('Issued original O2 integral row required')
        node=issued[0];index,key,jet=self.roles[node];cell=self.route[index];left,right=self.bounds[index]
        if row.get('operation')!='definite_integral' or not row.get('exact_function_integral'):
            raise ValueError('Exact original definite integral required')
        if (row['lower'],row['upper'])!=(cell['lower'],cell['upper']):raise ValueError('Unchanged original O2 endpoints required')
        source=self.sources[node];actual_key,actual_jet=self.provider.require_row(source)
        if source['chart']!='O2_slope' or (actual_key,actual_jet)!=(key,jet):raise ValueError('Original O2 C0/Z source role differs')
        if g.nodes[source['coordinate']]!={'operation':'bound_variable','name':row['variable']}:
            raise ValueError('Same original O2 coordinate required')
        x=sy.Symbol(row['variable'],real=True);rate=transport.RATES[key]
        expression,density=reference.symbolic_node(self.provider,row['integrand'],source)
        expected=sy.exp(-sy.Rational(rate.numerator,rate.denominator)*(right-x))*density
        if sy.simplify(expression-expected)!=0:raise ValueError('Original endpoint/kernel/source/Jacobian differs')
        if row['measure']!='native coordinate; original dy/dcoordinate applied exactly once':
            raise ValueError('Original O2 native measure required')
        return dict(integral_node=node,route_index=index,route_label=cell['label'],key=key,ordinary_order=jet,
            source_node=source['source_node'],source_role=source['function_role'],exact_endpoints=[str(left),str(right)],
            exact_reduced_integrand='exp(-rate*(original_route_right-coordinate))*issued_original_signed_density',
            rate=str(rate),original_coordinate_Jacobian=1,own_rate_kernel_endpoint_not_replaced_by_one=True,
            actual_issued_full_signed_DAG_executed_directly=True)

    def phase_boxes(self,index,left,right,N):
        source=self.sources[self.node_ids[index,'m','C0']]
        phase=self.provider.phase.phase_for_source(source,self.provider.built,coordinate=str(left),N=N)
        c=self.c;shift=c.mpf((0,ep(c.mpf(N)*(c.mpf(int((right-left).p))/int((right-left).q)))[1]))
        projection=precise.periodic_add(c,[c.mpf(ep(v)) for v in phase['actual_phase_directed_boxes']],shift)
        return projection['boxes'],dict(phase,exact_source_cell_phase_increment='N*(right-left)',
            archived_N7_phase_rows_not_used=True,ordinary_Z_radius_phase_derivative_exact_zero=True)

    def frame(self,*,Z,N,count=64,bits=24):
        z=leaves.base.point.pressure.exact_Z(Z)
        if z==0:raise ValueError('Nonzero O2 Z required by this signed integral callback')
        if type(N) is not int or N<160 or N.bit_length()>4096:raise ValueError('Explicit original N>=160, at most4096 bits required')
        if type(count) is not int or count not in self.cells.parent.parent.levels:raise ValueError('Accepted original defining-source grid level required')
        if type(bits) is not int or not 8<=bits<=256:raise ValueError('Inverse bits in[8,256] required')
        self.require_integral(self.rows[0,'m','C0']);key=(z,N,count,bits)
        if key in self.cache:return self.cache[key]
        c=self.c;records=[];began=time.monotonic();values={node:c.mpf(0) for node in self.roles}
        with mp.workdps(c.dps+40):
            for index,(start,end) in enumerate(self.bounds):
                for i in range(int(sy.floor(start*count)),int(sy.ceiling(end*count))):
                    # Restrict a certified whole source-cell range to its
                    # overlap with the exact original decimal route interval.
                    # No point interpolation or new source quadrature occurs.
                    left=max(sy.Rational(i,count),start);right=min(sy.Rational(i+1,count),end)
                    if not left<right:raise ValueError('Positive exact source/route intersection required')
                    query=self.cells.parent.query(count,i,Z_lower=str(z),Z_upper=str(z))
                    if not query['pieces']:raise ArithmeticError('Original O2 signed source requires refinement')
                    boxes,phase=self.phase_boxes(index,left,right,N);ranges={jet:{k:[] for k in transport.RATES} for jet in ('C0','Z')};proofs=[]
                    for source_index,original in enumerate(query['pieces']):
                        part=cells.correlated_p2_Z_source(original,query['record']);kernel=part['kernel']
                        for phi in boxes:
                            if ep(phi)==(0,1):
                                chart='E' if kernel.geometry=='signed_Mobius' else 'psi';coordinate=c.mpf((0,1))
                                inverse=dict(status='enclosed',chart=chart,coordinate_interval=coordinate,
                                    actual_source_phase=phi,proof='original strict inverse maps a complete period onto[0,1]')
                                derivative=inverse
                            else:
                                inverse0=kernel.evaluate(phi,bits=bits)
                                if inverse0['status']!='enclosed':raise ArithmeticError('Original O2 source inverse requires refinement')
                                inverse=inverse0['selected_inverse'];coordinate=inverse['coordinate_interval'];chart=inverse['chart']
                                derivative=cells.derivative_inverse(kernel,phi,bits)
                            primitive=kernel.primitives(coordinate,chart)
                            jets,derivative_contract=cells.slow_values(part,derivative['coordinate_interval'],derivative['chart'])
                            graph=positive.CellDensityGraph(self.view,part,primitive,jets,N);got=graph.outputs()
                            for jet,label in (('C0','densities'),('Z','density_Z')):
                                for density_key,value in got[label].items():
                                    source=self.sources[self.node_ids[index,density_key,jet]]
                                    if graph.evaluate(source['source_node']) is not value:
                                        raise ValueError('Actual issued full signed density root was not executed')
                                    if value.ctx is not c or value.ledger is not part['ledger'] or value.scale.bases is not part['q'].scale.bases:
                                        raise ValueError('One local source context/basis/ledger required')
                                    ranges[jet][density_key].append(positive.bounded(value))
                            proofs.append(dict(source_piece_index=source_index,actual_true_phase_cover=phi,
                                original_inverse=inverse,original_derivative_inverse=derivative,
                                original_slow_Z_contract=derivative_contract,
                                correlated_original_p2_Z_carrier=part['p2_Z_common_carrier_binding'],
                                actual_density_graph_root_contract=graph.record(),
                                local_source_factors_canceled_before_finite_export=True,
                                nonlinear_signed_DAG_before_source_and_phase_union=True))
                    hull=lambda rows:c.mpf((min(ep(v)[0] for v in rows),max(ep(v)[1] for v in rows)))
                    densities={jet:{k:hull(v) for k,v in rows.items()} for jet,rows in ranges.items()};masses={};contributions={}
                    for density_key,rate in transport.RATES.items():
                        rr=c.mpf(rate.numerator)/rate.denominator;L=c.mpf(int(left.p))/int(left.q);R=c.mpf(int(right.p))/int(right.q);B=c.mpf(int(end.p))/int(end.q)
                        mass=R-L if not rate else (c.exp(-rr*(B-R))-c.exp(-rr*(B-L)))/rr
                        if ep(mass)[0]<=0:raise ArithmeticError('Positive original O2 subwindow mass lost')
                        masses[density_key]=mass;contributions[density_key]={}
                        for jet in ('C0','Z'):
                            node=self.node_ids[index,density_key,jet]
                            contribution=densities[jet][density_key]*mass;values[node]+=contribution
                            contributions[density_key][jet]=contribution
                    records.append(dict(route_index=index,exact_y_cell=[str(left),str(right)],source=query['record'],
                        accepted_whole_source_cell_restricted_as_function_range=True,
                        actual_original_phase_origin_and_increment=phase,actual_phase_boxes=boxes,
                        actual_source_piece_and_density_graph_records=proofs,five_signed_C0_Z_density_hulls=densities,
                        exact_positive_own_rate_masses=masses,C0_Z_integral_contributions=contributions,
                        original_route_kernel_right_endpoint=str(end),original_coordinate_Jacobian_applied_once=1,
                        alternatives_unioned_not_added_as_source_masses=True))
                print('Actual original O2 route integrals:',count,'N',N,self.route[index]['label'],flush=True)
            whole={k:{jet:c.mpf(0) for jet in ('C0','Z')} for k in transport.RATES}
            for index,(_,end) in enumerate(self.bounds):
                B=c.mpf(int(end.p))/int(end.q)
                for density_key,rate in transport.RATES.items():
                    decay=c.exp(-(c.mpf(rate.numerator)/rate.denominator)*(1-B))
                    for jet in ('C0','Z'):
                        node=self.node_ids[index,density_key,jet]
                        whole[density_key][jet]+=decay*values[node]
        if leaves.digest_rows(self.provider.built['graph'].nodes)!=self.provider.graph_digest:
            raise ValueError('Original graph changed during integral evaluation')
        record=dict(source_family=self.family,source_graph_sha256=self.source_graph_sha256,original_Z_exact=str(z),
            explicit_candidate_N=N,accepted_defining_source_cell_level=count,inverse_bits=bits,
            exact_original_integral_bindings=self.contracts,original_parameter_and_phase_identity=self.parameter_phase_binding,
            actual_C0_Z_original_graph_integral_intervals={str(node):value for node,value in values.items()},
            original_O2_full_window_contribution_intervals=whole,genuine_source_cell_phase_integral_records=records,
            actual_full_signed_density_DAG_not_cap_values_executed=True,
            archived_N7_phase_or_integral_results_not_rescaled_or_evaluated=True,
            original_source_quad_and_pressure_late_errors_retained=True,
            positive_q_and_correlated_p2_Z_and_source_sign_branches_retained=True,
            normalized_intervals_exported_before_cross_cell_sum=True,
            original_incoming_corrections_and_separate_P0_not_reset=True,
            error_contract='Accepted entire defining-source cells and pressure errors; current actual N/radius phase covers; certified inverse and genuine ordinary-Z formulas; nonlinear signed graph before branch union; positive own-rate masses to each original route endpoint.',
            fixed_nonzero_Z_partial_callback_not_whole_Z_terminal_proof=True,
            full_scalar_control_evaluator_compatibility_installed=False,
            full_17_chart_source_or_24_cell_integral_oracle_installed=False,
            actual_five_controls_installed=False,current_whole_N_selected=False,execution_seconds=time.monotonic()-began)
        frame=OriginalO2IntegralFrame(z,N,count,bits,values,whole,record);self.cache[key]=frame;self.issued[id(frame)]=frame;return frame

    def dispatch(self,row,frame):
        contract=self.require_integral(row)
        if type(frame) is not OriginalO2IntegralFrame or self.issued.get(id(frame)) is not frame:
            raise ValueError('Issued genuine original O2 integral frame required')
        value=frame.values[contract['integral_node']]
        if not hasattr(value,'_mpi_') or hasattr(value,'scale') or value.ctx is not self.c:
            raise TypeError('Ordinary directed same-source integral interval required')
        return value

    def integrate(self,row,*,Z,N,count=64,bits=24):
        self.require_integral(row);return self.dispatch(row,self.frame(Z=Z,N=N,count=count,bits=bits))
    def parameter(self,name):raise NotImplementedError('Full original parameter evaluator remains open')
    def source(self,*args,**kwargs):raise NotImplementedError('Complete source/evaluator adapter remains open')


@precise.phase.native.inlet.source_precision
def run(*,return_live=False):
    began=time.monotonic();owner=OriginalO2IntegralCallback();frames=[]
    for count,N in ((64,160),(256,160),(64,257)):
        frame=owner.frame(Z='.37',N=N,count=count);frames.append(frame)
        for row in owner.rows.values():owner.integrate(row,Z='.37',N=N,count=count)
    result=dict(**{GATE:True},source_family=owner.family,source_graph_sha256=owner.source_graph_sha256,
        mode=owner.mode,exact_original_integral_bindings=owner.contracts,
        actual_original_O2_integral_frames=[frame.record for frame in frames],
        actual_original_graph_integral_callbacks=50*len(frames),
        original_source_integral_graph_unchanged=leaves.digest_rows(owner.provider.built['graph'].nodes)==owner.provider.graph_digest,
        full_17_chart_source_or_24_cell_integral_oracle_installed=False,
        full_scalar_control_evaluator_compatibility_installed=False,
        actual_five_controls_installed=False,actual_terminal_Z_function_closure_installed=False,
        current_whole_N_selected=False,**dict.fromkeys(precise.phase.packets.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Genuine issued O2 five-route C0/Z integral callbacks at fixed nonzero Z and actual common candidate N>=160, using original source-cell errors, new precise phase and directly executed original signed density roots. Not an N7 rescale, full-Z terminal closure, full scalar evaluator, global N, recursion or corrected NS.')
    raw=json.dumps(precise.encode(result),indent=2).encode()+b'\n'
    (HERE/NAME).write_bytes(gzip.compress(raw,compresslevel=9,mtime=0))
    return (result,owner,frames) if return_live else result


if __name__=='__main__':run()
