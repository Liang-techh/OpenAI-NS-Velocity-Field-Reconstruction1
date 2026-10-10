"""Portable source-owned O2/O3 partial histories, with the true dlogR measure.

This stops at the original O3_power endpoint. The selected pulse/collar/heat
registry and physical-time reconstruction are separate unfinished layers.
"""
from fractions import Fraction
import copy
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_Rh_C3_continuation as rh

target,ranges,current,source=rh.target,rh.ranges,rh.current,rh.source
HERE,PREFIX,sha,ep=rh.HERE,rh.PREFIX,rh.sha,rh.ep
NAME=PREFIX+'current_original_outer_C3_continuation.json.gz'
RECEIPT=PREFIX+'current_original_outer_C3_continuation_check.json'
GATES=('current_original_outer_C3_partial_and_complete_history_chain_installed',
    'current_original_outer_C3_whole_native_history_ranges_installed')
OPEN=rh.OPEN
CHARTS=('O2_slope','O2_axial','O2_buffer','O3_transition','O3_power')


def substitute(g,jet,variable,value):
    return target.C3Function(*[g.node('function_substitution',expression=q.node,
        variable=variable.node,value=value.node,Z_independent_substitution=True)
        for q in target.rows(jet)])


def interval_constant(c,expr):
    """Evaluate only a parameter-free exact physical width, without float casts."""
    if expr.is_Rational:return c.mpf(int(expr.p))/int(expr.q)
    if expr.func==sy.exp:return c.exp(interval_constant(c,expr.args[0]))
    if expr.is_Add:return sum((interval_constant(c,q) for q in expr.args),c.mpf(0))
    if expr.is_Mul:
        result=c.mpf(1)
        for q in expr.args:result*=interval_constant(c,q)
        return result
    if expr.is_Pow and expr.args[1].is_Integer:return interval_constant(c,expr.args[0])**int(expr.args[1])
    raise ValueError('Parameter-free exact outer radius width required '+str(expr))


def build_chart(field,chart,leading_inlet):
    g=field.graph;alg=target.C3Algebra(g);ref=lambda i:source.FunctionRef(g,i)
    jet=lambda row:target.C3Function(*map(ref,row));w=field.windows[chart];geo=field.geometry[chart]
    native=g.symbol('native_'+chart);endpoint=g.symbol(chart+'_history_endpoint')
    left,right=map(ref,geo['domain']);offset=ref(geo['offset']);Jacobian=ref(geo['Jacobian'])
    at_endpoint=g.node('function_substitution',expression=offset.node,variable=native.node,
        value=endpoint.node,Z_independent_substitution=True)
    distance=g.sub(at_endpoint,ref(geo['left_offset']))
    primitive=field.raw['actual_C3_primitive_functions'][chart]
    E0,V0=(jet(primitive['source_roots'][key]) for key in ('E','V'))
    EE=alg.mul(E0,E0)
    leading_density=dict(m=V0,h=E0,k=alg.mul(E0,V0),
        e=alg.add(alg.mul(V0,V0),alg.scale(EE,'-1/2')),p=alg.scale(EE,'1/2'))
    incoming={key:jet(row) for key,row in w['incoming'].items()}
    density={key:jet(row) for key,row in w['density'].items()}
    def duhamel(inlets,densities,kind):
        local,partial={},{}
        for key,rate in current.RATES.items():
            memory=g.one if rate==0 else g.unary('exp',g.neg(g.mul(g.constant(rate),distance)))
            kernel=g.one if rate==0 else g.unary('exp',g.neg(g.mul(g.constant(rate),g.sub(at_endpoint,offset))))
            def integrate(q,j):
                if q==g.zero:return g.zero
                return current.controls.integral(g,g.mul(kernel,q,Jacobian),'native_'+chart,left,endpoint,
                    measure='original dlogR: native Jacobian exactly once',actual_current_chart=chart,
                    history_kind=kind,own_rate=str(rate),ordinary_slow_Z_derivative_order=j,
                    endpoints_kernel_radius_Jacobian_and_global_phase_Z_independent=True)
            local[key]=target.C3Function(*[integrate(q,j) for j,q in enumerate(target.rows(densities[key]))])
            partial[key]=alg.add(alg.scale(inlets[key],memory),local[key])
        return local,partial
    local,partial=duhamel(incoming,density,'actual N-scaled correction')
    leading_local,leading=duhamel(leading_inlet,leading_density,'actual original leading history')
    invN=field.rh.functions['inverse_N']
    F=alg.fixed(0) if w['F_N'] is None else jet(w['F_N'])
    B=jet(primitive['B'])
    E,V=(substitute(g,q,native,endpoint) for q in
        (alg.add(E0,alg.scale(F,invN)),alg.add(V0,alg.scale(B,invN))))
    complete={key:alg.add(leading[key],alg.scale(partial[key],invN)) for key in current.RATES}
    EE=alg.mul(E,E)
    first=dict(m=alg.add(V,alg.neg(complete['m'])),h=alg.add(E,alg.scale(complete['h'],'-3/2')),
        k=alg.add(alg.mul(E,V),alg.scale(complete['k'],'-3/2')),
        e=alg.add(alg.mul(V,V),alg.scale(EE,'-1/2'),alg.neg(complete['e'])),p=alg.scale(EE,'1/2'))
    leading_out={key:substitute(g,q,endpoint,right) for key,q in leading.items()}
    P0=field.rh.functions['actual_independent_P0_C3']
    return dict(chart=chart,endpoint=endpoint,native_coordinate=native,domain=geo['domain'],
        exact_native_domain=field.domains[chart],actual_offset_at_endpoint=at_endpoint,actual_native_Jacobian=Jacobian,
        source_incoming_N_scaled_C3=incoming,source_density_N_scaled_C3=density,
        local_partial_N_scaled_C3=local,partial_N_scaled_C3=partial,
        original_full_outgoing_N_scaled_C3={key:jet(row) for key,row in w['outgoing'].items()},
        actual_leading_incoming_C3=leading_inlet,actual_leading_density_C3=leading_density,
        actual_leading_local_C3=leading_local,actual_leading_history_C3=leading,actual_leading_outgoing_C3=leading_out,
        actual_corrected_E_C3=E,actual_corrected_V_C3=V,actual_complete_history_C3=complete,
        actual_complete_history_y_C3=first,actual_independent_P0_C3=P0,actual_absolute_pressure_C3=alg.add(P0,complete['p']),
        pressure_datum_same_exact_Rh_function=True,quiet_source_does_not_reset_incoming=True)


class CurrentOuterC3Continuation:
    def __init__(self,require_checked=True):
        self.rh=rh.CurrentRhC3Continuation();self.graph=self.rh.graph
        self.hashes={**self.rh.hashes,rh.NAME:sha(rh.NAME),rh.RECEIPT:sha(rh.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)}
        self.identity,self.raw,self.connection=self.rh.identity,self.rh.raw,self.rh.connection
        self.prefix=copy.deepcopy(self.graph.nodes);self.N=self.rh.N;self.c=self.rh.c
        raw=self.connection['exact_current_integral_and_control_graph']
        self.windows={w['chart']:w for w in self.rh.windows if w['chart'] in CHARTS}
        self.geometry={w['chart']:w for w in raw['windows'] if w['chart'] in CHARTS}
        offsets,domains,x=current.exact_geometry();self.offsets,self.x=offsets,x
        self.domains={key:[int(q) for q in domains[key]] for key in CHARTS}
        self.geometry_bindings=self.bind_geometry(raw)
        self.P0_registry=self.bind_P0_registry()
        self.rows={tuple(row['exact_Z_cell']):{w['chart']:w for w in row['actual_17_chart_C3_transports'] if w['chart'] in CHARTS}
            for row in self.rh.range_report['actual_four_Z_C3_target_transports']}
        g=self.graph;f=self.rh.functions
        leading={key:substitute(g,q,f['endpoint'],g.zero) for key,q in f['actual_leading_history_C3'].items()}
        self.initial_leading=leading;self.functions={}
        for chart in CHARTS:
            self.functions[chart]=build_chart(self,chart,leading)
            leading=self.functions[chart]['actual_leading_outgoing_C3']
        self.acceptance_loaded=False
        if require_checked:
            receipt=rh.read(RECEIPT)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual outer C3 history chain required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed actual outer continuation source '+name)
            self.acceptance_loaded=True

    def bind_geometry(self,raw):
        g=self.graph;parameters={key:source.FunctionRef(g,node) for key,node in raw['parameters'].items()}
        bindings={};previous=self.rh.geometry;before=self.rh.window
        for chart in CHARTS:
            w=self.geometry[chart];native=g.symbol('native_'+chart);left,right=self.domains[chart]
            expr=lambda q:source.expression(g,q,parameters,native)
            expected=dict(offset=expr(self.offsets[chart]).node,left_offset=expr(self.offsets[chart].subs(self.x,left)).node,
                right_offset=expr(self.offsets[chart].subs(self.x,right)).node,
                Jacobian=expr(sy.diff(self.offsets[chart],self.x)).node,domain=[g.constant(left).node,g.constant(right).node])
            for key,value in expected.items():
                if w[key]!=value:raise ValueError('Actual outer geometry changed '+chart+' '+key)
            if w['left_offset']!=previous['right_offset']:raise ValueError('Exact original radius seam required '+chart)
            if w['phase']!=g.unary('fractional_part',g.mul(self.N,source.FunctionRef(g,w['offset']))).node:
                raise ValueError('Same original absolute phase and selected N required '+chart)
            if self.windows[chart]['incoming']!=before['outgoing']:raise ValueError('Nonzero predecessor correction required '+chart)
            if self.windows[chart]['original_own_rate_memory']!=w['memory']:raise ValueError('Original own-rate memory required '+chart)
            primitive=self.raw['actual_C3_primitive_functions'][chart]
            if any(rows[:2]!=w['source_roots'][key] for key,rows in primitive['source_roots'].items()):
                raise ValueError('Original C3 root prefix required '+chart)
            bindings[chart]=dict(**expected,phase=w['phase'],selected_N=self.N.node,
                previous_chart=before['chart'],exact_predecessor_radius_and_C3_correction_memory_retained=True)
            previous,before=w,self.windows[chart]
        return bindings

    def bind_P0_registry(self):
        guards={};rh_packets={tuple(p['exact_Z_cell']):p for p in self.raw['actual_four_Z_17_chart_third_source_packets'] if p['chart']=='Rh_reference'}
        source_rows=self.rh.range_report['actual_four_Z_C3_target_transports']
        for chart in CHARTS:
            packets=[p for p in self.raw['actual_four_Z_17_chart_third_source_packets'] if p['chart']==chart]
            if len(packets)!=4:raise ValueError('Four actual P0 source cells required '+chart)
            guards[chart]=[]
            for p in packets:
                ends=tuple(p['exact_Z_cell']);window=next(w for row in source_rows if tuple(row['exact_Z_cell'])==ends
                    for w in row['actual_17_chart_C3_transports'] if w['chart']==chart)
                if (p['source_family']!=self.identity or not p['same_live_basis_context_and_ledger']
                    or not p['live_original_phase_Z_exact_zero']
                    or p['actual_independent_P0_ordinary_Z3']!=rh_packets[ends]['actual_independent_P0_ordinary_Z3']
                    or any(not cell['actual_source_and_primitive_range']['proof']['exact_same_source_P0'] for cell in window['cells'])):
                    raise ValueError('Same original P0 source registry required '+chart)
                guards[chart].append(dict(exact_Z_cell=ends,actual_source_packet=p,
                    original_source_P0_identity_guard=current.ast_binding(target.CurrentC3SourceTargets.source_packet),
                    exact_original_Rh_P0_function_reused=True,range_equality_is_only_secondary_registry_check=True))
        return guards

    def correction_functions(self,chart,endpoint=None):
        if chart not in CHARTS:raise ValueError('Admitted original O2/O3 native chart required')
        f=self.functions[chart]
        if endpoint is None:return f['partial_N_scaled_C3']
        q=Fraction(*endpoint) if isinstance(endpoint,tuple) else Fraction(endpoint)
        left,right=self.domains[chart]
        if not left<=q<=right:raise ValueError('Use next original native chart outside '+chart)
        if q==left:return f['source_incoming_N_scaled_C3']
        if q==right:return f['original_full_outgoing_N_scaled_C3']
        return {key:substitute(self.graph,row,f['endpoint'],self.graph.constant(q)) for key,row in f['partial_N_scaled_C3'].items()}

    def quantitative_chain(self,ends):
        ends=tuple(ends)
        if ends not in self.rows:raise ValueError('Admitted actual axial cell required')
        c=self.c;bd=ranges.Bounds(c)
        readcap=lambda rec:ranges.LogUpper(c,None if rec['exact_zero'] else current.packets.interval(c,rec['log_absolute_upper']))
        jet=lambda rows:ranges.JetBound(*(readcap(rec) for rec in rows[:4]))
        inlet=self.rh.quantitative_range(ends)
        leading={key:jet(q) for key,q in inlet['actual_leading_partial_history_C3_bounds'].items()}
        correction={key:jet(q) for key,q in inlet['actual_N_scaled_partial_history_C3_bounds'].items()}
        P0=jet(inlet['actual_independent_P0_C3_bounds']);epsilon=readcap(inlet['actual_selected_inverse_N_upper']);windows=[]
        for chart in CHARTS:
            before_leading,before_correction=dict(leading),dict(correction)
            Ecap=Vcap=Fcap=Bcap=bd.fixed(bd.zero);cells=[]
            for cell in self.rows[ends][chart]['cells']:
                left,right=(sy.Rational(*cell[key]) for key in ('exact_left','exact_right'))
                width_expr=sy.simplify(self.offsets[chart].subs(self.x,right)-self.offsets[chart].subs(self.x,left))
                width=interval_constant(c,width_expr)
                if ep(width)[0]<=0:raise ValueError('Strictly positive physical dlogR width required')
                E,V=(jet(cell['actual_source_and_primitive_range'][key]) for key in ('E','V'))
                Ecap,Vcap=bd.add(Ecap,E),bd.add(Vcap,V)
                Fcap=bd.add(Fcap,jet(cell['actual_density_ranges']['F']))
                Bcap=bd.add(Bcap,jet(cell['actual_source_and_primitive_range']['B']))
                EE=bd.product(E,E);ld=dict(m=V,h=E,k=bd.product(E,V),
                    e=bd.add(bd.product(V,V),bd.scaled(EE,bd.constant(c.mpf('.5')))),p=bd.scaled(EE,bd.constant(c.mpf('.5'))))
                masses={}
                for key,rate in current.RATES.items():
                    r=c.mpf(rate.numerator)/rate.denominator;mass=width if rate==0 else -c.expm1(-r*width)/r
                    if ep(mass)[0]<=0:raise ValueError('Strictly positive original own-rate mass required')
                    masses[key]=mass
                    leading[key]=bd.add(leading[key],bd.scaled(ld[key],bd.constant(mass)))
                    correction[key]=bd.add(correction[key],bd.scaled(jet(cell['actual_density_ranges']['N_scaled'][key]),bd.constant(mass)))
                cells.append(dict(exact_left=cell['exact_left'],exact_right=cell['exact_right'],
                    exact_physical_width=str(width_expr),physical_width=width,original_positive_own_rate_masses=masses,
                    nonunit_Jacobian_in_radius_width_once=True))
            complete={key:bd.add(leading[key],bd.scaled(correction[key],epsilon)) for key in current.RATES}
            E,V=bd.add(Ecap,bd.scaled(Fcap,epsilon)),bd.add(Vcap,bd.scaled(Bcap,epsilon));EE=bd.product(E,E)
            first=dict(m=bd.add(V,complete['m']),h=bd.add(E,bd.scaled(complete['h'],bd.constant(c.mpf('1.5')))),
                k=bd.add(bd.product(E,V),bd.scaled(complete['k'],bd.constant(c.mpf('1.5')))),
                e=bd.add(bd.product(V,V),bd.scaled(EE,bd.constant(c.mpf('.5'))),complete['e']),p=bd.scaled(EE,bd.constant(c.mpf('.5'))))
            windows.append(dict(chart=chart,exact_native_domain=self.domains[chart],original_native_cells=cells,
                actual_leading_incoming_C3_bounds=ranges.record(before_leading),actual_N_scaled_incoming_C3_bounds=ranges.record(before_correction),
                actual_leading_partial_history_C3_bounds=ranges.record(leading),actual_N_scaled_partial_history_C3_bounds=ranges.record(correction),
                actual_complete_partial_history_C3_bounds=ranges.record(complete),actual_complete_history_y_C3_bounds=ranges.record(first),
                actual_corrected_E_V_C3_bounds=ranges.record(dict(E=E,V=V)),actual_independent_P0_C3_bounds=ranges.record(P0),
                actual_absolute_pressure_C3_bound=ranges.record(bd.add(P0,complete['p'])),
                predecessor_decay_and_suffix_at_most_one=True,range_covers_all_native_partial_endpoints=True,
                actual_incoming_majorant_from_previous_whole_partial_cover=True,range_caps_not_function_values=True))
        return dict(source_family=self.identity,exact_Z_cell=ends,actual_five_chart_C3_partial_ranges=windows,
            original_Rh_inlet_cover_retained=True,final_selected_pulse_collar_heat_registry_not_admitted=True)


def source_bindings():
    return dict(build_chart=current.ast_binding(build_chart),exact_endpoint_substitution=current.ast_binding(substitute),
        exact_physical_width_interval=current.ast_binding(interval_constant),
        actual_geometry_binding=current.ast_binding(CurrentOuterC3Continuation.bind_geometry),
        original_P0_registry_binding=current.ast_binding(CurrentOuterC3Continuation.bind_P0_registry),
        correction_aliases=current.ast_binding(CurrentOuterC3Continuation.correction_functions),
        quantitative_chain=current.ast_binding(CurrentOuterC3Continuation.quantitative_chain),
        original_density_transport=current.ast_binding(target.build_transport),original_geometry=current.ast_binding(current.exact_geometry))


def run():
    began=time.monotonic()
    with mp.workdps(540):
        field=CurrentOuterC3Continuation(require_checked=False)
        rows=[field.quantitative_chain(ends) for ends in field.rows]
        report=dict(candidate_actual_outer_C3_continuation_constructed=True,source_family=field.identity,
            accepted_Rh_C3_graph_prefix_length=len(field.prefix),exact_graph_nodes=field.graph.nodes,
            actual_five_chart_C3_continuation_functions=target.encoded(field.functions),
            actual_geometry_bindings=field.geometry_bindings,original_same_P0_source_registry=field.P0_registry,
            actual_four_Z_outer_C3_partial_ranges=ranges.record(rows),source_bindings=source_bindings(),
            **dict.fromkeys(GATES+OPEN,False),current_numeric_point_field_oracle_installed=False,
            global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed=False,
            final_selected_native_pulse_postpulse_collar_heat_registry_installed=False,
            input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(target.encoded(report),separators=(',',':'))+'\n').encode(),mtime=0))
    print('Actual O2/O3 C3 partial and complete histories with native dlogR weights constructed',flush=True)
    return field


if __name__=='__main__':run()
