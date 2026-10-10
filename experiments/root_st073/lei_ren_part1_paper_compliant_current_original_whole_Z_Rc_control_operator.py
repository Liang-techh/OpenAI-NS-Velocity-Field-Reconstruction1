"""Actual whole-Z Rc target functions in the original five-control map.

Exact expression handles define finite Picard functions. Directed source,
weight and iterate ranges enclose them; endpoints never become functions.
Nonzero incoming defects persist. Finite iterates are not solved controls.
"""
import ast
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_sharp_functions as current
import lei_ren_part1_paper_compliant_current_native_Rc_functional_controls as algebra

relative = current.previous.relative
repair = relative.repair
HERE,PREFIX,sha,bind,ep = current.HERE,current.PREFIX,current.sha,current.bind,current.ep
CELLS,OPEN = current.CELLS,current.OPEN
ROWS,CONTROLS = repair.ROWS,repair.CONTROLS
Pair = algebra.C1Enclosure
NAME = PREFIX+'current_original_whole_Z_Rc_control_operator.json.gz'
RECEIPT = PREFIX+'current_original_whole_Z_Rc_control_operator_check.json'
GATE = 'current_original_whole_Z_actual_Rc_control_map_and_finite_C1_functions_installed'


def add(a,b):return Pair(a.value+b.value,a.Z+b.Z)
def negative(a):return Pair(-a.value,-a.Z)
def scale(a,b):return Pair(a.value*b,a.Z*b)
def multiply(a,b):return Pair(a.value*b.value,a.Z*b.value+a.value*b.Z)


def inverse_action(matrix,vector):
    zero = vector[0].value.scalar(0)
    return [Pair(sum((q.value*b for q,b in zip(vector,row)),zero),
        sum((q.Z*b for q,b in zip(vector,row)),zero)) for row in matrix['inverse_enclosure']]


def control_step(c,mu,logmu,N,matrix,d,h):
    if type(N) is not int or N < 1:raise ValueError('One actual positive integer N required')
    if len(d) != 5 or len(h) != 5:raise ValueError('Five source rows and five controls required')
    reference = d[0].value
    for q in (*d,*h):
        reference.coerce(q.value); reference.coerce(q.Z)
    Q = algebra.range_quadratic(c,mu,logmu,matrix,h)
    return [negative(q) for q in inverse_action(matrix,
        [add(a,scale(b,c.mpf(1)/N)) for a,b in zip(d,Q)])]


def picard_ranges(c,mu,logmu,N,matrix,d,iterations):
    if type(iterations) is not int or not 1 <= iterations <= 16:
        raise ValueError('Explicit finite Picard depth in [1,16] required')
    zero = d[0].value.scalar(0)
    sequence = [[Pair(zero,zero) for unused in range(5)]]
    for unused in range(iterations):sequence.append(control_step(c,mu,logmu,N,matrix,d,sequence[-1]))
    # Exact function identity B*h_j+d=-Q(h_(j-1))/N removes the
    # linear cancellation before enclosure. No independent endpoints
    # are subtracted to assert an exact-zero linear residual.
    last = algebra.range_quadratic(c,mu,logmu,matrix,sequence[-1])
    prior = algebra.range_quadratic(c,mu,logmu,matrix,sequence[-2])
    residual = [scale(add(a,negative(b)),c.mpf(1)/N) for a,b in zip(last,prior)]
    return sequence,residual


def pair_record(p):return [p.value.record(),p.Z.record()]
def pairs_record(rows,labels):return {key:pair_record(q) for key,q in zip(labels,rows)}


def power_binding():
    name = PREFIX+'current_original_O3_Rc_finite_N.py'
    tree = ast.parse((HERE/name).read_text(encoding='utf8'))
    fn = next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name == 'background_cell')
    statements = ('sig=c.mpf(1);factor=c.exp(-t/2-mu_cover*t)',
        'E=f.scale(u1,factor)', 'V,Vy=zero,zero')
    for text in statements:
        for wanted in ast.parse(text).body:
            if sum(ast.dump(n)==ast.dump(wanted) for n in ast.walk(fn)) != 1:
                raise ValueError('Defining original power-band velocity identity changed')
    reserve = PREFIX+'current_generic_shear_O3_sources.py'
    text = (HERE/reserve).read_text(encoding='utf8')
    if 'Utheta=Ac_theta(Z)*(R/Rc)^(-1/2-mu), Uz=0' not in text \
            or 'Ac_theta(Z)=Utheta(Rw,Z)*exp(-2*(1/2+mu))' not in text:
        raise ValueError('Original same-amplitude Rc..2Rc reservation required')
    return dict(source_modules={name:sha(name),reserve:sha(reserve)},
        original_power_velocity='E_Rw*exp(-(1/2+mu)*t), V=0',
        same_Rc_amplitude='E_Rw*exp(-2*(1/2+mu))',
        exact_band_change_of_variables='t=2+log(x), x=R/Rc',
        original_band_velocity='E=A_rc*x^(-1/2-mu), V=0',
        original_power_driver_exact_flat=True,
        nonzero_predecessor_moment_defects_retained=True,
        reservation_does_not_admit_repaired_terminal_or_heat_exterior=True)


class WholeZRcControlOperator:
    def __init__(self,dps=500):
        self.owner = current.WholeZSharpBridgeFunctions(dps)
        self.c,self.identity,self.N = self.owner.c,self.owner.identity,self.owner.N
        checked = json.loads((HERE/current.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(current.GATE) \
                or checked['source_family'] != self.identity or checked['candidate_N'] != self.N:
            raise ValueError('Checked actual sharp whole-Z target functions required')
        self.hashes = dict(self.owner.hashes)
        for name,digest in checked['input_hashes'].items():
            if sha(name) != digest:raise ValueError('Changed actual target prerequisite: '+name)
            bind(self.hashes,name,digest)
        bind(self.hashes,current.RECEIPT,sha(current.RECEIPT))
        self.saved = json.loads(gzip.decompress((HERE/current.NAME).read_bytes()))
        if not self.saved[current.GATE] or self.saved['source_family'] != self.identity \
                or self.saved['candidate_N'] != self.N:
            raise ValueError('Same actual original source and N required')
        self.cells = {tuple(row['exact_Z_cell']):row for row in self.saved['actual_sharp_bridge_and_Rc_cells']}
        if set(self.cells) != set(CELLS):raise ValueError('All four actual source domains required')
        self.binding = power_binding()
        for name,digest in self.binding['source_modules'].items():bind(self.hashes,name,digest)
        for module in (algebra,algebra.source,repair):
            name = Path(module.__file__).name; bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        self.logmu = self.owner.source.logmu
        self.mu = self.c.exp(self.logmu)
        if ep(self.mu)[0] <= 0:raise ValueError('Original actual positive mu required')
        self.weights = repair.fresh_weights(self.c,self.mu,cells=512)
        self.matrix = repair.fresh_linear_inverse(self.c,self.mu,self.weights)
        self.graph_binding = dict(provider_module=Path(current.__file__).name,
            provider_module_sha256=sha(Path(current.__file__).name),provider_gate=current.GATE,
            provider_report=current.NAME,provider_report_sha256=sha(current.NAME),
            provider_receipt=current.RECEIPT,provider_receipt_sha256=sha(current.RECEIPT),
            defining_live_method='WholeZSharpBridgeFunctions.retransport',
            source_family=self.identity,candidate_N=self.N,
            target_definition=relative.REFERENCE,reference_binding=relative.reference_binding(),
            fixed_input_pressure_P0_cancels_only_in_relative_difference=True,
            source_functions_are_integrals_not_saved_interval_endpoints=True)
        self.graph_digest = hashlib.sha256(json.dumps(self.graph_binding,sort_keys=True).encode()).hexdigest()

    def build(self,iterations=3):
        g = algebra.source.FunctionTransportGraph()
        mu = g.node('actual_whole_Z_source_parameter',name='mu',source_binding_sha256=self.graph_digest)
        N = g.node('shared_positive_integer',name='N',actual_integer=self.N)
        def source_node(quantity,key,order):
            return g.node('actual_whole_Z_Rc_source_function',quantity=quantity,row=key,Z_order=order,
                source_binding_sha256=self.graph_digest,Z_variable='Z',shared_N=N.node,
                defining_provider='WholeZSharpBridgeFunctions.retransport',
                numerical_value_not_installed=True)
        targets = {key:algebra.source.C1Function(source_node('N_scaled_relative_target',key,0),
            source_node('N_scaled_relative_target',key,1)) for key in ROWS}
        amplitude = algebra.source.C1Function(source_node('original_Rc_amplitude','E',0),
            source_node('original_Rc_amplitude','E',1))
        built = dict(graph=g,parameters=dict(mu=mu),N=N,N_scaled_targets=targets,
            amplitude=amplitude,source_family=self.identity,source_graph_sha256=self.graph_digest,
            source_roles=self.graph_binding)
        return algebra.exact_control_graph(built,iterations=iterations)

    def target_functions(self,ends,mode='live'):
        key = tuple(ends)
        if key not in self.cells:raise ValueError('Admitted actual whole-Z source cell required')
        f = self.owner.owner(key).flow
        if mode == 'live':
            row = self.owner.retransport(key)
            targets = row['actual_fixed_N_N_scaled_target_C0_Z']
            amplitude = row['actual_positive_A_rc_C0_Z']
        elif mode == 'accepted':
            row = self.cells[key]
            decode = lambda v:relative.restore_row(f,v)
            if any(not relative.exact_row(v,w) for v,w in zip(self.owner.owner(key).P0,
                    [decode(v) for v in row['exact_common_P0_axial5']])):
                raise ValueError('Same actual P0 and Z derivatives required')
            targets = {k:[decode(v) for v in pair] for k,pair in row['actual_fixed_N_N_scaled_target_C0_Z'].items()}
            amplitude = [decode(v) for v in row['actual_positive_A_rc_C0_Z']]
            outer = self.owner.cells[key]['actual_source_owned_zero_inlet_to_Rc_signed_stages']['outer']
            original = outer['actual_original_Rc_signed_source_function']['actual_original_source_packet']
            if not all(v['exact_zero'] for v in original['original_q_C0_Z'].values()):
                raise ValueError('Same original Rc power velocity must be exactly flat')
            expected = [decode(v) for v in outer['actual_positive_Rc_normalized_amplitude_C0_Z']]
            if any(not relative.exact_row(a,b) for a,b in zip(amplitude,expected)):
                raise ValueError('Actual same original endpoint amplitude required')
        else:raise ValueError('Explicit live or accepted actual source evaluation required')
        relative.same_source(f,[v for pair in targets.values() for v in pair]+amplitude)
        return dict(flow=f,d=[Pair(*targets[key]) for key in ROWS],amplitude=Pair(*amplitude),
            exact_common_P0_axial5=self.owner.owner(key).P0,evaluation_mode=mode,
            source_binding_sha256=self.graph_digest)

    def controls(self,ends,iterations=3,mode='live'):
        target = self.target_functions(ends,mode)
        sequence,residual = picard_ranges(self.c,self.mu,self.logmu,self.N,self.matrix,target['d'],iterations)
        return dict(**target,sequence=sequence,residual=residual,built=self.build(iterations),
            exact_Z_cell=list(ends),candidate_N=self.N,iterations=iterations,
            finite_iterates_are_not_solved_controls=True)

    def band(self,live,x):
        """Original supported finite-iterate trial velocities in x=R/Rc."""
        c = self.c; x = c.mpf(x)
        if not 1 <= ep(x)[0] <= ep(x)[1] <= 2:raise ValueError('Original Rc..2Rc band required')
        if live['source_binding_sha256'] != self.graph_digest or live['candidate_N'] != self.N:
            raise ValueError('Same original source-bound finite iterate required')
        w = self.weights; logx = c.ln(x)
        bumps = [repair.raw_beta(c,(logx-center)/w['radius'])/(w['radius']*w['raw_normalization']*x)
            for center in w['centers']]
        h = live['sequence'][-1]; zero = h[0].value.scalar(0); zpair = Pair(zero,zero)
        F = zpair; G = zpair
        for b,q in zip(bumps,h[2:]):F=add(F,scale(q,b/self.N))
        for b,q in zip((bumps[0],bumps[2]),h[:2]):G=add(G,scale(q,b/self.N))
        power = c.exp(-(c.mpf('.5')+self.mu)*logx)
        A = live['amplitude']
        E = add(scale(A,power),multiply(A,F)); V = multiply(A,G)
        return dict(F=F,G=G,E=E,V=V,
            leading_E=scale(A,power),leading_V=zpair,
            actual_control_solution_certified=False,complete_band_terminal_closure_certified=False)


def run():
    began = time.monotonic()
    with mp.workdps(540):
        owner = WholeZRcControlOperator(); built=owner.build(); rows=[]
        for ends in CELLS:
            live = owner.controls(ends,mode='accepted')
            bands = {x:{key:pair_record(value) for key,value in owner.band(live,x).items()
                if isinstance(value,Pair)} for x in ('1','1.5','2')}
            rows.append(dict(exact_Z_cell=list(ends),candidate_N=owner.N,
                exact_common_P0_axial5=live['exact_common_P0_axial5'],
                source_evaluation_mode=live['evaluation_mode'],source_binding_sha256=owner.graph_digest,
                actual_N_scaled_target_C0_Z=pairs_record(live['d'],ROWS),
                finite_control_sequence_C0_Z=[pairs_record(h,CONTROLS) for h in live['sequence']],
                actual_finite_iterate_residual_C0_Z=pairs_record(live['residual'],ROWS),
                original_supported_trial_band_C0_Z=bands,
                actual_fixed_N_solution_certified=False,actual_terminal_controls_installed=False))
            print('Actual Rc control map/three finite iterates: '+str(ends),flush=True)
        roots=lambda pairs,labels:{key:dict(value=q.value.node,Z=q.Z.node) for key,q in zip(labels,pairs)}
        report=dict(**{GATE:True},source_family=owner.identity,candidate_N=owner.N,
            source_function_binding=owner.graph_binding,source_binding_sha256=owner.graph_digest,
            original_fixed_leading_band_binding=owner.binding,exact_function_graph_nodes=built['graph'].nodes,
            exact_finite_control_roots=[roots(h,CONTROLS) for h in built['finite_picard_sequence']],
            exact_control_residual_roots=roots(built['control_residual'],ROWS),
            exact_trial_band_roots={k:dict(value=v.value.node,Z=v.Z.node) for k,v in built['band_profiles'].items()},
            exact_weight_enclosures=owner.weights,exact_matrix_enclosures=owner.matrix,
            actual_four_Z_control_map_cells=rows,finite_iterate_depth=3,
            actual_source_bound_control_map_callable=True,actual_scalar_point_control_oracle_installed=False,
            actual_terminal_controls_installed=False,actual_fixed_N_solution_certified=False,
            actual_global_frequency_admitted=False,physical_original_exterior_five_targets_closed=False,
            **dict.fromkeys(OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
            scope='Actual same-fixed-input whole-Z Rc target oracle, exact original five-control map '
                'and three C1 finite functions with source-bound enclosures and supported trial band. '
                'Nonzero incoming/pressure memory retained. No solved control, contraction tail, '
                'terminal zero, global N, heat/stress or temporal recursion admission.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(current.encode(current.serialized(report)),
            separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__ == '__main__':run()
