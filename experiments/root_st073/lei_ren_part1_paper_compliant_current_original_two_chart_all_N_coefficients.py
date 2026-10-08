"""Live two-chart actual source roles and original all-N coefficient values.

True point/error inputs, actual shared-N phase and inverse are consumed.
The original exact coefficient program is interpreted with directed formal
factors, separately on every phase piece before nonlinear coefficient unions.
The legacy scalar evaluator and its guards are unchanged. This is a local
coefficient runtime, not the 24-cell integral/five-control/global-N oracle.
"""
from dataclasses import dataclass
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
from types import MappingProxyType
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_reference_point_oracle as reference
import lei_ren_part1_paper_compliant_current_native_Rc_all_N_function_controls as exact

base=reference.base;slow=reference.slow;point=base.point;prior=base.prior
HERE,PREFIX,sha=reference.HERE,reference.PREFIX,reference.sha
NAME=PREFIX+'current_original_two_chart_all_N_coefficients.json'
RECEIPT=PREFIX+'current_original_two_chart_all_N_coefficients_check.json'
GATE='original_two_chart_live_role_dispatch_and_exact_N_dependent_C0_Z_coefficients_evaluated'
ep=reference.ep
ROLES={
    'all_N_original_E_C0':('E',(0,0)),
    'all_N_original_E_Z':('E',(0,1)),
    'all_N_original_V_C0':('V',(0,0)),
    'all_N_original_V_Z':('V',(0,1)),
    'all_N_periodic_A_C0':('A',(0,0)),
    'all_N_periodic_A_Z':('A',(0,1)),
    'all_N_periodic_B_C0':('B',(0,0)),
    'all_N_periodic_B_Z':('B',(0,1)),
}
CHARTS=('Rh_reference','O2_slope')


def candidate_N(N):
    if type(N) is not int or N<160:raise ValueError('One explicit integer candidate N>=160 required')
    return N


def same_interval(a,b):return a._mpi_==b._mpi_


def directed_exprel(c,argument):
    """Exact defining exponential average with a directed Taylor remainder."""
    radius=max(abs(x) for x in ep(argument))
    if radius>1:raise ValueError('Exprel source argument must be bounded by1')
    if not radius:return c.mpf(1)
    mean=c.mpf(1);power=c.mpf(1);order=64
    for k in range(1,order+1):
        power*=argument;mean+=power/c.factorial(k+1)
    # Integrate the exponential Taylor remainder over t in[0,1].
    tail=ep(c.exp(1)*c.mpf(radius)**(order+1)/c.factorial(order+2))[1]
    return mean+c.mpf((-tail,tail))


def union(values):
    if not values:raise ValueError('Nonempty original phase pieces required')
    first=values[0];result=first
    # Union after the original nonlinear expression is evaluated per piece.
    for value in values[1:]:
        first.scale.pair(value.scale)
        if value.ledger is not first.ledger:raise ValueError('Same point/source ledger required')
        reference_log=max(ep(result.scale.evaluate())[1],ep(value.scale.evaluate())[1])
        scale=prior.FormalScale(first.scale.bases,offset=first.ctx.mpf(reference_log))
        left=result.coefficient*result.bounded_exp((result.scale-scale).evaluate())
        right=value.coefficient*value.bounded_exp((value.scale-scale).evaluate())
        result=prior.ScaledEnclosure(scale,first.ctx.mpf((min(ep(left)[0],ep(right)[0]),max(ep(left)[1],ep(right)[1]))),first.ledger)
    return result


@dataclass(frozen=True)
class OriginalPointPiece:
    chart: str
    coordinate: object
    Z: object
    N: int
    family: dict
    graph_sha256: str
    phase: object
    values: dict
    record: dict

    @property
    def ctx(self):return self.values['all_N_original_E_C0'].ctx


class TwoChartOriginalPointSources:
    mode='original_point_factored_enclosures'
    def __init__(self):
        self.reference=reference.OriginalReferencePointOracle()
        self.O2=slow.OriginalO2ConditionedSlowZ()
        self.source_family=self.reference.family
        if self.O2.family!=self.source_family:raise ValueError('One original point family required')
        self.source_graph_sha256=sha(point.source.inertial.profiles.VIEWS)
        self.views=json.loads(gzip.decompress((HERE/point.source.inertial.profiles.VIEWS).read_bytes()))
        self.hashes={}
        for stem,gate in (('current_original_reference_point_oracle',reference.GATE),
            ('current_original_O2_conditioned_slow_Z',slow.GATE),
            ('current_native_Rc_all_N_function_controls',exact.GATE)):
            name=PREFIX+stem+'_check.json';receipt=json.loads((HERE/name).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(gate) or receipt['source_family']!=self.source_family:
                raise ValueError('Checked same-family original point/coefficient program required')
            for filename,digest in {**receipt['input_hashes'],name:sha(name)}.items():
                if sha(filename)!=digest:raise ValueError('Actual point source changed: '+filename)
                if filename in self.hashes and self.hashes[filename]!=digest:raise ValueError('Point hash closures disagree')
                self.hashes[filename]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.cache={};self.issued_pieces={};self.issued_graphs={}

    def source_identity(self,chart,role):
        if chart not in CHARTS:raise ValueError('Unsupported actual point chart: '+str(chart))
        if role not in ROLES:raise ValueError('Explicit original all-N E/V/A/B C0/Z role required')
        view=self.views[chart];name,order=ROLES[role]
        if view['source_family']!=self.source_family:raise ValueError('Original view family differs')
        if name=='E':
            namespace='original_signed_input_graph.jet_expression_dag'
            node=view['original_signed_input_graph']['jet_expression_dag']['roots']['E']['y0_Z'+str(order[1])]
        elif name=='V':
            namespace='function_graph_nodes';label='V_y0_Z'+str(order[1])
            nodes=[i for i,r in enumerate(view['function_graph_nodes']) if r['operation']=='source_derivative' and r['name']==label]
            if len(nodes)!=1:raise ValueError('Unique original axial source row required')
            node=nodes[0]
        else:
            namespace='function_graph_nodes'
            key=('A' if order==(0,0) else 'A_Z_slow') if name=='A' else ('B_over_Pstar' if order==(0,0) else 'B_Z_slow')
            node=view['roots'][key]
        return namespace,node

    def require_row(self,row,piece):
        if type(piece) is not OriginalPointPiece or self.issued_pieces.get(id(piece)) is not piece or piece.family!=self.source_family or piece.graph_sha256!=self.source_graph_sha256:
            raise ValueError('Same typed original point/family/graph frame required')
        if row.get('operation')!='original_function_graph' or row.get('graph_sha256')!=self.source_graph_sha256:
            raise ValueError('Original function graph source required')
        if row.get('chart')!=piece.chart or row.get('graph_file')!=point.source.inertial.profiles.VIEWS:
            raise ValueError('Same original source chart and graph file required')
        namespace,node=self.source_identity(piece.chart,row.get('function_role'))
        if row.get('source_graph_namespace')!=namespace or row.get('source_node')!=node:
            raise ValueError('Original role namespace/node mismatch')
        if row.get('Z_variable')!='Z' or not row.get('source_coefficients_are_function_recipes_not_cover_values'):
            raise ValueError('Original source-function contract required')
        return row['function_role']

    def source(self,row,*,piece,coordinate,phase,Z,N):
        role=self.require_row(row,piece);c=piece.ctx
        if candidate_N(N)!=piece.N or point.pressure.exact_Z(Z)!=piece.Z:
            raise ValueError('Same actual N/Z as point frame required')
        expected=c.mpf(int(piece.coordinate.p))/int(piece.coordinate.q)
        if not same_interval(c.mpf(coordinate),expected) or not same_interval(c.mpf(phase),piece.phase):
            raise ValueError('Exact bound coordinate and owned actual phase piece required')
        value=piece.values[role]
        if type(value) is not prior.ScaledEnclosure or value.ctx is not c:
            raise TypeError('Original live directed factored function value required')
        return value

    def pieces(self,*,chart,coordinate,Z,N,bits=80):
        N=candidate_N(N);self.source_identity(chart,'all_N_original_E_C0')
        coordinate=reference.reference_coordinate(coordinate) if chart=='Rh_reference' else point.source.exact_coordinate(coordinate)
        Z=point.pressure.exact_Z(Z)
        if type(bits) is not int or not 4<=bits<=256:raise ValueError('Inverse bits in[4,256] required')
        key=(chart,coordinate,Z,N,bits)
        if key in self.cache:return self.cache[key]
        owner=self.reference.owner if chart=='Rh_reference' else self.O2.owner
        phase=owner.radius.evaluate(y=coordinate,N=N);c=owner.ctx;pieces=[]
        with mp.workdps(c.dps+40):
            query=(reference.reference_query(owner,y=coordinate,Z=Z) if chart=='Rh_reference' else owner.query(y=coordinate,Z=Z))
            kernel=query['kernel']
            for box in phase['true_original_phase_directed_boxes']:
                target=c.mpf([box['lower'],box['upper']]);inverse=kernel.evaluate(target,bits=bits)
                if inverse['status']!='enclosed':raise ArithmeticError('Refine actual original point source/phase; no fallback')
                selected=inverse['selected_inverse'];primitive=kernel.primitives(selected['coordinate_interval'],selected['chart'])
                jets,derivative=slow.slow_values(kernel,query['roots'],selected['coordinate_interval'],selected['chart'])
                values={}
                for role,(name,order) in ROLES.items():
                    if name in ('E','V'):value=query['roots'][name][order]
                    elif name=='A':value=primitive['A'] if order==(0,0) else jets['A_Z_slow']
                    else:value=primitive['B_over_Pstar'] if order==(0,0) else jets['B_Z_slow']
                    values[role]=value
                anchor=values['all_N_original_E_C0']
                for value in values.values():
                    anchor.scale.pair(value.scale)
                    if value.ledger is not anchor.ledger:raise ValueError('One common point basis/ledger required')
                record=dict(chart=chart,source_family=self.source_family,original_coordinate_exact=str(coordinate),
                    original_Z_exact=str(Z),candidate_N=N,actual_owned_phase=target,
                    original_radius_phase=phase,actual_inverse=inverse,ordinary_Z_contract=derivative,
                    original_point_coefficients_and_errors=point.record_query(query['point']),
                    all_original_role_values={r:v.record() for r,v in values.items()},
                    true_phase_Z_exact_zero=True,one_basis_and_ledger_for_all_roles=True,
                    B_units='B_over_Pstar',range_endpoint_or_midpoint_selected=False)
                piece=OriginalPointPiece(chart,coordinate,Z,N,self.source_family,self.source_graph_sha256,target,MappingProxyType(values),record)
                self.issued_pieces[id(piece)]=piece;pieces.append(piece)
        self.cache[key]=pieces;return pieces


def coefficient_graph(dispatcher,chart):
    """Build the original exact coefficient program with bound source roles."""
    g=exact.source.FunctionTransportGraph();coordinate=g.symbol('original_coordinate')
    phase=g.symbol('owned_original_phase');N=g.node('shared_positive_integer',name='N')
    refs={}
    for role in ROLES:
        namespace,node=dispatcher.source_identity(chart,role)
        refs[role]=g.node('original_function_graph',chart=chart,graph_file=point.source.inertial.profiles.VIEWS,
            graph_sha256=dispatcher.source_graph_sha256,source_node=node,coordinate=coordinate.node,
            phase=phase.node,shared_N=N.node,Z_variable='Z',source_graph_namespace=namespace,function_role=role,
            source_coefficients_are_function_recipes_not_cover_values=True)
    pairs={name:exact.source.C1Function(refs['all_N_'+prefix+'_C0'],refs['all_N_'+prefix+'_Z'])
        for name,prefix in (('E','original_E'),('V','original_V'),('A','periodic_A'),('B','periodic_B'))}
    coefficients,F=exact.coefficient_pairs(g,pairs['E'],pairs['V'],pairs['A'],pairs['B'],N)
    built=dict(graph=g,coefficients=coefficients,F=F,N=N,refs=refs,source_family=dispatcher.source_family,
        source_graph_sha256=dispatcher.source_graph_sha256)
    dispatcher.issued_graphs[id(built)]=built
    return built


class FactoredPointCoefficientEvaluator:
    """Evaluate a bounded modulation graph without weakening scalar guards."""
    def __init__(self,built,dispatcher,piece):
        if type(dispatcher) is not TwoChartOriginalPointSources or type(piece) is not OriginalPointPiece:
            raise TypeError('Typed live original dispatcher and phase piece required')
        if dispatcher.issued_graphs.get(id(built)) is not built or dispatcher.issued_pieces.get(id(piece)) is not piece or built['source_family']!=piece.family or built['source_graph_sha256']!=piece.graph_sha256:
            raise ValueError('Same original exact coefficient graph/family required')
        self.built=built;self.dispatcher=dispatcher;self.piece=piece;self.ctx=piece.ctx
        self.scalar=piece.values['all_N_original_E_C0'].scalar;self.cache={}
    def __call__(self,root):
        self.built['graph'].ids([root]);return self.walk(root.node)
    def walk(self,index):
        if index in self.cache:return self.cache[index]
        g=self.built['graph'];row=g.nodes[index];op=row['operation'];c=self.ctx;child=self.walk
        if op=='exact_rational':value=self.scalar(c.mpf(row['numerator'])/row['denominator'])
        elif op=='bound_variable':
            if row['name']=='original_coordinate':value=self.scalar(c.mpf(int(self.piece.coordinate.p))/int(self.piece.coordinate.q))
            elif row['name']=='owned_original_phase':value=self.scalar(self.piece.phase)
            else:raise ValueError('Only source-owned point variables supported')
        elif op=='shared_positive_integer':value=self.scalar(candidate_N(self.piece.N))
        elif op=='original_function_graph':
            if row['shared_N']!=self.built['N'].node:raise ValueError('One common original N node required')
            value=self.dispatcher.source(row,piece=self.piece,coordinate=child(row['coordinate']).finite_interval(),
                phase=child(row['phase']).finite_interval(),Z=self.piece.Z,N=self.piece.N)
        elif op=='sum':value=sum((child(i) for i in row['arguments']),self.scalar(0))
        elif op=='negative':value=-child(row['argument'])
        elif op=='product':
            args=row['arguments']
            if len(args)==2 and args[0]==args[1]:value=base.current.square(child(args[0]))
            else:
                value=self.scalar(1)
                for i in args:value=value*child(i)
        elif op=='positive_quotient':
            denominator=child(row['denominator'])
            if ep(denominator.coefficient)[0]<=0:raise ArithmeticError('Actual strictly positive denominator required')
            lower=ep(denominator.scale.evaluate()+c.ln(c.mpf(ep(denominator.coefficient)[0])))[0]
            value=child(row['numerator']).positive_divide(denominator,lower)
        elif op=='analytic_unary':
            argument=base.conditioned.bounded_value(child(row['argument']))
            if max(abs(x) for x in ep(argument))>1:raise ArithmeticError('Bounded original modulation argument required')
            if row['name']=='exprel':
                value=self.scalar(directed_exprel(c,argument))
                value.ledger['directed_exprel_order64_remainders']=value.ledger.get('directed_exprel_order64_remainders',0)+1
            elif row['name']=='exp':value=self.scalar(c.exp(argument))
            else:raise ValueError('Unsupported factored analytic operation')
        else:raise ValueError('Coefficient-only evaluator rejects unsupported integral/control operation: '+op)
        if type(value) is not prior.ScaledEnclosure:raise TypeError('Directed factored values required throughout')
        self.cache[index]=value;return value


def evaluate_coefficients(dispatcher,*,chart,coordinate,Z,N,bits=80):
    built=coefficient_graph(dispatcher,chart);pieces=dispatcher.pieces(chart=chart,coordinate=coordinate,Z=Z,N=N,bits=bits)
    evaluated=[]
    for piece in pieces:
        evaluator=FactoredPointCoefficientEvaluator(built,dispatcher,piece)
        values={order:{key:dict(C0=evaluator(pair.value),Z=evaluator(pair.Z))
            for key,pair in rows.items()} for order,rows in built['coefficients'].items()}
        F=dict(C0=evaluator(built['F'].value),Z=evaluator(built['F'].Z))
        evaluated.append(dict(piece=piece,values=values,F=F))
    values={order:{key:{jet:union([p['values'][order][key][jet] for p in evaluated]) for jet in ('C0','Z')}
        for key in exact.RATES} for order in exact.ORDERS}
    record=dict(source_family=dispatcher.source_family,chart=chart,original_coordinate_exact=str(pieces[0].coordinate),
        original_Z_exact=str(pieces[0].Z),candidate_N=N,source_graph_sha256=dispatcher.source_graph_sha256,
        source_role_count=len(ROLES),exact_coefficient_program=Path(exact.__file__).name,
        exact_coefficient_expression_nodes=built['graph'].nodes,
        actual_phase_pieces=[dict(source=p['piece'].record,F_N={k:v.record() for k,v in p['F'].items()},
            coefficients={str(o):{k:{j:v.record() for j,v in pair.items()} for k,pair in rows.items()} for o,rows in p['values'].items()})
            for p in evaluated],
        actual_N_dependent_C0_Z_coefficients={str(o):{k:{j:v.record() for j,v in pair.items()} for k,pair in rows.items()} for o,rows in values.items()},
        phase_piece_nonlinear_coefficients_evaluated_before_union=True,
        original_E_V_B_and_all_cross_terms_retained=True,exact_order_minus1_functions_not_zeroed=True,
        original_scalar_evaluator_guards_unchanged=True,source_caps_or_midpoints_used_as_values=False,
        all_17_chart_or_24_cell_integral_oracle_installed=False,actual_five_controls_installed=False,
        current_whole_N_selected=False,**dict.fromkeys(point.source.inertial.profiles.loop.OPEN,False))
    return dict(record=record,values=values,evaluated=evaluated,built=built)


def run():
    began=time.monotonic();dispatcher=TwoChartOriginalPointSources();records=[]
    for chart,coordinate,Z,N in (('Rh_reference','-2.337','.37',160),('Rh_reference','-2.337','0',160),
        ('Rh_reference','-2.337','.37',256),('O2_slope','.53','.37',160),('O2_slope','.53','0',160)):
        got=evaluate_coefficients(dispatcher,chart=chart,coordinate=coordinate,Z=Z,N=N);records.append(got['record'])
        print('Live original all-N coefficient functions:',chart,coordinate,Z,N,flush=True)
    result=dict(**{GATE:True},source_family=dispatcher.source_family,mode=dispatcher.mode,
        actual_original_point_coefficient_queries=records,supported_original_charts=list(CHARTS),
        implemented_original_source_roles=list(ROLES),original_exact_coefficient_program_reused=True,
        partial_live_role_dispatch_and_factored_coefficient_runtime_installed=True,
        original_scalar_evaluator_guards_unchanged=True,no_ancestor_source_producer_reexecuted=True,
        all_17_chart_or_24_cell_integral_oracle_installed=False,actual_five_controls_installed=False,
        current_whole_N_selected=False,**dict.fromkeys(point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=dispatcher.hashes,execution_seconds=time.monotonic()-began,
        scope='Live actual Rh_reference/O2_slope E/V/A/B C0/Z roles and original exact N-dependent -1/-2 coefficient program at common candidate N>=160. Directed factored runtime evaluates nonlinear coefficients before phase-piece union. Not a full scalar/integral/five-control/global-N/recursion oracle.')
    (HERE/NAME).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
