"""Portable actual Rh continuation: signed C3 correction and complete histories.

Restore accepted exact function handles; do not replay expensive numerical
owners or choose enclosure endpoints as values. Rh offset a is in [-5,0].
"""
from fractions import Fraction
import copy
import gzip
import json
import math
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_C3_source_targets as target
import lei_ren_part1_paper_compliant_current_original_C3_target_ranges as ranges
import lei_ren_part1_paper_compliant_current_original_C3_limit_controls as controls
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rh_functional_join as join

phase,current,source=target.phase,target.current,target.source
HERE,PREFIX,sha,ep=target.HERE,target.PREFIX,target.sha,target.ep
NAME=PREFIX+'current_original_Rh_C3_continuation.json.gz'
RECEIPT=PREFIX+'current_original_Rh_C3_continuation_check.json'
GATES=('current_original_Rh_C3_partial_correction_and_complete_histories_installed',
    'current_original_Rh_C3_whole_reference_history_magnitude_ranges_installed')
OPEN=target.OPEN


def read(name):
    data=(HERE/name).read_bytes()
    return json.loads(gzip.decompress(data) if name.endswith('.gz') else data)


def accepted_inputs():
    hashes={};reports={}
    for module in (target,ranges,controls,join):
        report,receipt=read(module.NAME),read(module.RECEIPT)
        if not receipt['all_passed'] or not all(receipt[k] for k in module.GATES):
            raise ValueError('Accepted actual source required '+module.RECEIPT)
        if report['source_family']!=receipt['source_family']:raise ValueError('Same report/receipt source required')
        for name,digest in receipt['input_hashes'].items():
            if name in hashes and hashes[name]!=digest:raise ValueError('Conflicting source provenance '+name)
            hashes[name]=digest
        hashes[module.NAME]=sha(module.NAME);hashes[module.RECEIPT]=sha(module.RECEIPT)
        reports[module.NAME]=report
    connection=read(current.NAME);reports[current.NAME]=connection;hashes[current.NAME]=sha(current.NAME)
    hashes[Path(source.__file__).name]=sha(Path(source.__file__).name);hashes[Path(__file__).name]=sha(Path(__file__).name)
    identity=reports[target.NAME]['source_family']
    if any(q['source_family']!=identity for q in reports.values()):raise ValueError('One actual source family required')
    for name,digest in hashes.items():
        if sha(name)!=digest:raise ValueError('Changed accepted source '+name)
    return reports,hashes


def restore_graph(report):
    g=source.FunctionTransportGraph();nodes=report['exact_graph_nodes']
    if nodes[:2]!=g.nodes:raise ValueError('Exact original zero/one prefix required')
    g.nodes=copy.deepcopy(nodes)
    g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(g.nodes)}
    if len(g.keys)!=len(g.nodes):raise ValueError('Duplicate accepted graph identity')
    return g


def build(field):
    g=field.graph;alg=target.C3Algebra(g);ref=lambda i:source.FunctionRef(g,i)
    jet=lambda rows:target.C3Function(*map(ref,rows));w=field.window;native=field.native
    a=g.symbol('Rh_reference_endpoint');left=g.constant(-5);distance=g.add(a,g.constant(5))
    kernel=lambda rate:g.unary('exp',g.neg(g.mul(g.constant(rate),g.sub(a,native))))
    incoming={key:jet(row) for key,row in w['incoming'].items()}
    density={key:jet(row) for key,row in w['density'].items()}
    local={};partial={}
    for key,rate in current.RATES.items():
        memory=g.one if rate==0 else g.unary('exp',g.neg(g.mul(g.constant(rate),distance)))
        local[key]=target.C3Function(*[current.controls.integral(g,g.mul(kernel(rate),q),
            'native_Rh_reference',left,a,measure='da=dlogR; Rh Jacobian one, applied once',
            actual_current_chart='Rh_reference',N_scaled_density=True,own_rate=str(rate),
            ordinary_slow_Z_derivative_order=j,endpoint_and_kernel_Z_independent=True)
            for j,q in enumerate(target.rows(density[key]))])
        partial[key]=alg.add(alg.scale(incoming[key],memory),local[key])
    substitute=lambda q:target.C3Function(*[g.node('function_substitution',expression=v.node,
        variable=native.node,value=a.node,Z_independent_substitution=True) for v in target.rows(q)])
    primitive=field.raw['actual_C3_primitive_functions']['Rh_reference']
    leadingE,leadingV=(substitute(jet(primitive['source_roots'][key])) for key in ('E','V'))
    F,B=substitute(jet(w['F_N'])),substitute(jet(primitive['B']))
    invN=g.quotient(g.one,field.N,'same exact positive selected repair-only integer')
    E,V=alg.add(leadingE,alg.scale(F,invN)),alg.add(leadingV,alg.scale(B,invN))
    recipe=copy.deepcopy(field.connection['exact_source_recipe_manifest']['Rh_reference'])
    recipe['quantity_paths'].update({**{('history_'+key):'original_generic_source.common_own_five_histories_axial5.'+key+'[Z_order]'
        for key in current.RATES},'P0':'original_generic_source.common_original_P0_axial5[Z_order]'})
    def leaf(quantity,coordinate):
        return target.C3Function(*[g.node('current_original_leading_function_recipe',source_family=field.identity,
            native_chart='Rh_reference',recipe=recipe,quantity=quantity,Z_order=j,coordinate=coordinate.node,
            Z_variable=g.symbol('Z').node,ordinary_slow_Z_derivative_order=j,Taylor_coefficient_factorial=math.factorial(j),
            source_projection_binding=current.ast_binding(CurrentRhC3Continuation.project_leading_packet),
            quantity_path=recipe['quantity_paths'][quantity],defining_quantity_not_a_range_value=True,
            leading_history_normalization_already_applied=True,phase_independent_leading_source=True,
            coordinate_independent=(quantity=='P0'),original_P0_source_registry_guard=field.P0_guards,
            exact_source_op_P0_function_not_enclosure_endpoint=(quantity=='P0')) for j in range(4)])
    leading={key:leaf('history_'+key,a) for key in current.RATES}
    P0=leaf('P0',left)
    complete={key:alg.add(leading[key],alg.scale(partial[key],invN)) for key in current.RATES}
    EE=alg.mul(E,E)
    first=dict(m=alg.add(V,alg.neg(complete['m'])),h=alg.add(E,alg.scale(complete['h'],'-3/2')),
        k=alg.add(alg.mul(E,V),alg.scale(complete['k'],'-3/2')),
        e=alg.add(alg.mul(V,V),alg.scale(EE,'-1/2'),alg.neg(complete['e'])),p=alg.scale(EE,'1/2'))
    seam=g.node('proved_current_Rh_same_source_C3_history_continuation_identity',source_family=field.identity,
        incoming_function_rows=target.encoded(incoming),actual_patch_outgoing_function_rows=w['incoming'],
        accepted_leading_join_receipt=join.RECEIPT,accepted_leading_join_receipt_sha256=sha(join.RECEIPT),
        original_complete_support_function_proof=field.join_report['exact_current_function_proof'],
        coordinate_identity='actual_patch x=e, Rh_reference a=-5; logR=logRm+1',
        correction_identity='partial at a=-5 is same predecessor H; complete=leading+H/N exactly once',
        global_phase_identity='same selected N and same original absolute log-radius offset',
        range_overlap_not_function_identity=True)
    return dict(endpoint=a,native_coordinate=native,reference_domain=[-5,0],
        original_N=field.N,inverse_N=invN,source_incoming_N_scaled_C3=incoming,source_density_N_scaled_C3=density,
        local_partial_N_scaled_C3=local,partial_N_scaled_C3=partial,
        original_full_outgoing_N_scaled_C3={key:jet(row) for key,row in w['outgoing'].items()},
        actual_leading_history_C3=leading,actual_corrected_E_C3=E,actual_corrected_V_C3=V,
        actual_independent_P0_C3=P0,actual_complete_history_C3=complete,actual_complete_history_y_C3=first,
        actual_absolute_pressure_C3=alg.add(P0,complete['p']),Rh_seam_identity=seam,
        actual_geometry_source_binding=field.geometry_binding,original_independent_P0_source_registry=field.P0_guards,
        axial_orders=[0,1,2,3],complete_history_y_orders=[0,1],
        selected_repair_integer_and_nonzero_incoming_retained=True,
        relative_zero_does_not_reset_absolute_pressure_or_memory=True,
        global_physical_time_Cartesian_heat_cone_and_temporal_recursion_not_admitted=True)


class CurrentRhC3Continuation:
    def __init__(self,require_checked=True):
        reports,self.hashes=accepted_inputs();self.raw=reports[target.NAME];self.identity=self.raw['source_family']
        self.connection=reports[current.NAME];self.range_report=reports[ranges.NAME]
        self.join_report=reports[join.NAME];self.control_report=reports[controls.NAME]
        self.graph=restore_graph(self.raw);self.prefix=copy.deepcopy(self.graph.nodes)
        raw=self.connection['exact_current_integral_and_control_graph'];self.N=source.FunctionRef(self.graph,raw['N'])
        if self.prefix[:len(self.connection['exact_function_graph_nodes'])]!=self.connection['exact_function_graph_nodes']:
            raise ValueError('Exact shared original graph prefix required')
        self.windows=self.raw['actual_C3_density_transport_targets']['windows']
        index=next(i for i,row in enumerate(self.windows) if row['chart']=='Rh_reference')
        self.window=self.windows[index];self.predecessor=self.windows[index-1]
        if self.predecessor['chart']!='actual_patch' or self.predecessor['outgoing']!=self.window['incoming']:
            raise ValueError('Actual patch nonzero predecessor function memory required')
        self.geometry=raw['windows'][index];self.native=self.graph.symbol('native_Rh_reference')
        self.geometry_binding=self.bind_geometry(raw)
        self.c=MPIntervalContext();self.c.dps=540
        self.rows={tuple(row['exact_Z_cell']):next(w for w in row['actual_17_chart_C3_transports'] if w['chart']=='Rh_reference')
            for row in self.range_report['actual_four_Z_C3_target_transports']}
        self.join_rows={tuple(row['join']['exact_Z_cell']):row['join'] for row in self.join_report['source_cells']}
        selected=self.connection['actual_whole_Z_frequency_connection']['selected_integer']
        if self.control_report['actual_quantitative_C3_bounds']['actual_same_repair_integer']!=selected:
            raise ValueError('Same actual selected repair integer required')
        self.P0_guards=self.bind_P0_registry()
        self.functions=build(self);self.acceptance_loaded=False
        if require_checked:
            receipt=read(RECEIPT)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual Rh continuation required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed Rh continuation source '+name)
            self.acceptance_loaded=True

    def bind_geometry(self,raw):
        g=self.graph;w=self.geometry;offsets,domains,x=current.exact_geometry()
        parameters={key:source.FunctionRef(g,node) for key,node in raw['parameters'].items()}
        expr=lambda q:source.expression(g,q,parameters,self.native)
        expected=dict(offset=expr(offsets['Rh_reference']).node,
            left_offset=expr(offsets['Rh_reference'].subs(x,-5)).node,
            right_offset=expr(offsets['Rh_reference'].subs(x,0)).node,
            Jacobian=g.one.node,domain=[g.constant(-5).node,g.zero.node])
        if w['chart']!='Rh_reference' or domains['Rh_reference']!=(-5,0):raise ValueError('Actual Rh reference chart/domain required')
        for key,value in expected.items():
            if w[key]!=value:raise ValueError('Actual Rh geometry changed '+key)
        if w['width']!=g.sub(source.FunctionRef(g,w['right_offset']),source.FunctionRef(g,w['left_offset'])).node:
            raise ValueError('Actual original reference width required')
        if w['phase']!=g.unary('fractional_part',g.mul(self.N,source.FunctionRef(g,w['offset']))).node:
            raise ValueError('Same selected N and original radius phase required')
        primitive=self.raw['actual_C3_primitive_functions']['Rh_reference']
        if any(rows[:2]!=w['source_roots'][key] for key,rows in primitive['source_roots'].items()):
            raise ValueError('Same original leading C3 root prefix required')
        if self.window['original_own_rate_memory']!=w['memory']:
            raise ValueError('Same original own-rate window memories required')
        return dict(**expected,width=w['width'],phase=w['phase'],selected_N=self.N.node,
            original_current_offset_and_phase_Z_independent=True,source_roots_C3_prefix_equals_geometry_roots=True)

    def bind_P0_registry(self):
        packets={tuple(row['exact_Z_cell']):row for row in self.raw['actual_four_Z_17_chart_third_source_packets']
            if row['chart']=='Rh_reference'}
        if packets.keys()!=self.join_rows.keys() or packets.keys()!=self.rows.keys():
            raise ValueError('All four same actual source P0 cells required')
        guards={}
        for ends,packet in packets.items():
            seam=self.join_rows[ends]
            if (packet['source_family']!=self.identity or seam['source_identity']!=self.identity
                    or not seam['shared_original_P0'] or not seam['current_unique_implicit_owner_shared_by_both_branches']
                    or not packet['same_live_basis_context_and_ledger'] or not packet['live_original_phase_Z_exact_zero']
                    or len(packet['actual_independent_P0_ordinary_Z3'])!=4
                    or any(not cell['actual_source_and_primitive_range']['proof']['exact_same_source_P0'] for cell in self.rows[ends]['cells'])):
                raise ValueError('Original op.P0 source registry identity required')
            guards[str(ends)]=dict(actual_C3_source_packet=packet,
                original_op_P0_identity_from_shared_owner=seam['shared_original_P0'],
                actual_source_registry_guard=current.ast_binding(target.CurrentC3SourceTargets.source_packet),
                accepted_same_owner_P0_guard=current.ast_binding(join.WholeZRhFunctionalJoin.seam),
                exact_function_is_same_op_P0_at_every_coordinate=True,
                numerical_enclosure_equality_not_used_as_function_identity=True)
        return guards

    @staticmethod
    def project_leading_packet(packet):
        generic=packet['original_generic_source']
        if packet['exact_common_P0_axial5'] is not generic['common_original_P0_axial5']:
            if not current.current.previous.equivalent_rows(packet['exact_common_P0_axial5'],generic['common_original_P0_axial5']):
                raise ValueError('Same independent source P0 rows required')
        def project(row):
            if len(row)<4:raise ValueError('Actual leading Taylor rows zero through three required')
            return [math.factorial(j)*row[j] for j in range(4)]
        return dict(histories={key:project(row) for key,row in generic['common_own_five_histories_axial5'].items()},
            P0=project(generic['common_original_P0_axial5']))

    def correction_functions(self,offset=None):
        if offset is None:return self.functions['partial_N_scaled_C3']
        q=Fraction(*offset) if isinstance(offset,tuple) else Fraction(offset)
        if not -5<=q<=0:raise ValueError('Rh_reference offset in [-5,0]; use next native chart beyond zero')
        if q==-5:return self.functions['source_incoming_N_scaled_C3']
        if q==0:return self.functions['original_full_outgoing_N_scaled_C3']
        g=self.graph;a=self.functions['endpoint'];value=g.constant(q)
        return {key:target.C3Function(*[g.node('function_substitution',expression=v.node,variable=a.node,
            value=value.node,Z_independent_substitution=True) for v in target.rows(row)])
            for key,row in self.functions['partial_N_scaled_C3'].items()}

    def quantitative_range(self,ends):
        ends=tuple(ends)
        if ends not in self.rows:raise ValueError('Admitted actual axial cell required')
        c=self.c;bd=ranges.Bounds(c);row=self.rows[ends];seam=self.join_rows[ends]
        readcap=lambda rec:ranges.LogUpper(c,None if rec['exact_zero'] else current.packets.interval(c,rec['log_absolute_upper']))
        jet=lambda rows:ranges.JetBound(*(readcap(rec) for rec in rows[:4]))
        incoming={key:jet(q) for key,q in row['actual_N_scaled_incoming_C3'].items()}
        leading={key:ranges.JetBound(*(bd.scale(readcap(rec),math.factorial(j)) for j,rec in enumerate(rows[:4])))
            for key,rows in seam['actual_reference_histories'].items()}
        correction=dict(incoming);Ecap=Vcap=bd.fixed(bd.zero);cells=[]
        for cell in row['cells']:
            left,right=(Fraction(*cell[key]) for key in ('exact_left','exact_right'))
            width=c.mpf((right-left).numerator)/(right-left).denominator
            E,V=(jet(cell['actual_source_and_primitive_range'][key]) for key in ('E','V'))
            Ecap,Vcap=bd.add(Ecap,E),bd.add(Vcap,V)
            EE=bd.product(E,E);ld=dict(m=V,h=E,k=bd.product(E,V),
                e=bd.add(bd.product(V,V),bd.scaled(EE,bd.constant(c.mpf('.5')))),p=bd.scaled(EE,bd.constant(c.mpf('.5'))))
            masses={}
            for key,rate in current.RATES.items():
                r=c.mpf(rate.numerator)/rate.denominator;mass=width if rate==0 else -c.expm1(-r*width)/r
                if ep(mass)[0]<=0:raise ValueError('Original positive reference own-rate mass required')
                masses[key]=mass
                leading[key]=bd.add(leading[key],bd.scaled(ld[key],bd.constant(mass)))
                correction[key]=bd.add(correction[key],bd.scaled(jet(cell['actual_density_ranges']['N_scaled'][key]),bd.constant(mass)))
            cells.append(dict(exact_left=cell['exact_left'],exact_right=cell['exact_right'],original_positive_own_rate_masses=masses))
        epsilon=readcap(self.control_report['actual_quantitative_C3_bounds']['actual_inverse_epsilon_upper'])
        complete={key:bd.add(leading[key],bd.scaled(correction[key],epsilon)) for key in current.RATES}
        # Source E correction F_N/N and V correction B/N; full-window sums are conservative sup covers.
        F=B=bd.fixed(bd.zero)
        for cell in row['cells']:
            F=bd.add(F,jet(cell['actual_density_ranges']['F']))
            B=bd.add(B,jet(cell['actual_source_and_primitive_range']['B']))
        E,V=bd.add(Ecap,bd.scaled(F,epsilon)),bd.add(Vcap,bd.scaled(B,epsilon));EE=bd.product(E,E)
        first=dict(m=bd.add(V,complete['m']),h=bd.add(E,bd.scaled(complete['h'],bd.constant(c.mpf('1.5')))),
            k=bd.add(bd.product(E,V),bd.scaled(complete['k'],bd.constant(c.mpf('1.5')))),
            e=bd.add(bd.product(V,V),bd.scaled(EE,bd.constant(c.mpf('.5'))),complete['e']),p=bd.scaled(EE,bd.constant(c.mpf('.5'))))
        P0=ranges.JetBound(*(bd.scale(readcap(rec),math.factorial(j)) for j,rec in enumerate(seam['exact_common_P0_axial5'][:4])))
        return dict(source_family=self.identity,exact_Z_cell=ends,whole_reference_offset=[-5,0],
            actual_N_scaled_partial_history_C3_bounds=ranges.record(correction),actual_leading_partial_history_C3_bounds=ranges.record(leading),
            actual_complete_partial_history_C3_bounds=ranges.record(complete),actual_complete_history_y_C3_bounds=ranges.record(first),
            actual_corrected_E_V_C3_bounds=ranges.record(dict(E=E,V=V)),actual_independent_P0_C3_bounds=ranges.record(P0),
            actual_absolute_pressure_C3_bound=ranges.record(bd.add(P0,complete['p'])),original_native_cells=cells,
            actual_selected_inverse_N_upper=epsilon.record(),predecessor_decay_and_suffix_at_most_one=True,
            range_covers_all_partial_endpoints=True,range_caps_not_function_values=True,
            original_nonzero_memory_and_independent_P0_retained=True)


def source_bindings():
    return dict(accepted_exact_graph_restore=current.ast_binding(restore_graph),build=current.ast_binding(build),
        actual_leading_Taylor_projection=current.ast_binding(CurrentRhC3Continuation.project_leading_packet),
        partial_endpoint_domain_and_aliases=current.ast_binding(CurrentRhC3Continuation.correction_functions),
        actual_Rh_geometry_binding=current.ast_binding(CurrentRhC3Continuation.bind_geometry),
        original_exact_P0_registry_binding=current.ast_binding(CurrentRhC3Continuation.bind_P0_registry),
        quantitative_range=current.ast_binding(CurrentRhC3Continuation.quantitative_range),
        original_C3_density_transport=current.ast_binding(target.build_transport),
        original_leading_Rh_function_identity=current.ast_binding(join.exact_join_proof),
        original_source_geometry=current.ast_binding(current.exact_geometry))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        field=CurrentRhC3Continuation(require_checked=False)
        rows=[field.quantitative_range(ends) for ends in field.rows]
        report=dict(candidate_actual_Rh_C3_continuation_constructed=True,source_family=field.identity,
            accepted_C3_source_graph_prefix_length=len(field.prefix),exact_graph_nodes=field.graph.nodes,
            actual_Rh_C3_continuation_functions=target.encoded(field.functions),actual_four_Z_Rh_C3_partial_ranges=ranges.record(rows),
            source_bindings=source_bindings(),**dict.fromkeys(GATES+OPEN,False),
            current_numeric_point_field_oracle_installed=False,global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed=False,
            input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(target.encoded(report),separators=(',',':'))+'\n').encode(),mtime=0))
    print('Actual source-owned Rh C3 partial correction, complete histories and whole-reference bounds constructed',flush=True)
    return field


if __name__=='__main__':run()
