"""Resolve original leading source recipes and the actual Rc leading interface.

The source replay interprets actual defining assignments in exact point/Taylor
algebra. Scalar enclosure helpers denote their already checked mathematical
integrals; no archived interval endpoint or magnitude is a function value.
"""
import ast
import copy
import gzip
import inspect
import json
import math
from pathlib import Path
import textwrap
import time
from types import SimpleNamespace
import sympy as s
import lei_ren_part1_paper_compliant_current_original_outer_leading_C3_identity as leading
import lei_ren_part1_paper_compliant_current_original_outer_leading_C3_identity_check as leading_check

outer,band,target,current,source=leading.outer,leading.band,leading.target,leading.current,leading.source
HERE,PREFIX,sha=leading.HERE,leading.PREFIX,leading.sha
NAME=PREFIX+'current_original_outer_leading_C3_bridge.json.gz'
RECEIPT=PREFIX+'current_original_outer_leading_C3_bridge_check.json'
GATES=('current_original_outer_leading_C3_provider_functions_resolved',
    'current_original_outer_leading_C3_Duhamel_functions_identified',
    'current_outer_leading_endpoint_band_seed_function_identity_installed',
    'current_outer_complete_history_to_repair_inlet_function_identity_installed')
OPEN=leading.OPEN


class ExactTaylor:
    """Ordinary Taylor coefficients, with exact order-three algebra."""
    def __init__(self,c,coefficients):self.coefficients=list(coefficients)[:4]
    def __getitem__(self,j):return self.coefficients[j]
    def coerce(self,q):return q if isinstance(q,ExactTaylor) else ExactTaylor(None,[s.sympify(q),0,0,0])
    def __add__(self,q):
        q=self.coerce(q);return ExactTaylor(None,[self[j]+q[j] for j in range(4)])
    __radd__=__add__
    def __neg__(self):return ExactTaylor(None,[-q for q in self.coefficients])
    def __sub__(self,q):return self+-self.coerce(q)
    def __rsub__(self,q):return self.coerce(q)+-self
    def __mul__(self,q):
        q=self.coerce(q);return ExactTaylor(None,[sum(self[k]*q[n-k] for k in range(n+1)) for n in range(4)])
    __rmul__=__mul__
    def reciprocal(self):
        out=[1/self[0]]
        for n in range(1,4):out.append(-sum(self[k]*out[n-k] for k in range(1,n+1))/self[0])
        return ExactTaylor(None,out)
    def __truediv__(self,q):return self*self.coerce(q).reciprocal()
    def exp(self):
        out=[s.exp(self[0])]
        for n in range(1,4):out.append(sum(k*self[k]*out[n-k] for k in range(1,n+1))/n)
        return ExactTaylor(None,out)
    def log(self):
        out=[s.log(self[0])]
        for n in range(1,4):out.append((n*self[n]-sum(k*out[k]*self[n-k] for k in range(1,n)))/(n*self[0]))
        return ExactTaylor(None,out)


class ExactPointFlow:
    def __init__(self,logP):self.logP=logP
    def scalar(self,q):return s.Rational(str(q)) if isinstance(q,float) else s.sympify(q)
    def factor(self,powers,offset=0):
        if any(q for j,q in enumerate(powers) if j!=1):raise ValueError('Unexpected leading source logarithmic basis')
        return s.exp(2*s.Rational(str(powers[1]))*self.logP+offset)
    def jet(self,q):return list(q.coefficients) if isinstance(q,ExactTaylor) else [self.scalar(q),0,0,0]
    def add(self,*rows):return [sum(row[n] for row in rows) for n in range(min(map(len,rows)))]
    def multiply(self,a,b):return [sum(a[k]*b[n-k] for k in range(n+1)) for n in range(min(len(a),len(b)))]
    def scale(self,row,q):return [v*q for v in row]
    def ordinary_cover(self,q):return q  # exact parameter denoted by its enclosure


def defining_slice(fn,names,env,branches=None):
    """Execute actual pure assignments; unsupported defining paths fail closed."""
    names=set(names);branches={} if branches is None else branches
    tree=ast.parse(textwrap.dedent(inspect.getsource(fn))).body[0]
    picked=[];seen=set()
    def visit(statements):
        for n in statements:
            if isinstance(n,ast.Assign):
                assigned={v.id for t in n.targets for v in ast.walk(t) if isinstance(v,ast.Name)}
                if assigned & names:
                    if not assigned<=names:raise ValueError('Partially selected defining assignment '+ast.unparse(n))
                    if assigned & seen:raise ValueError('Duplicate defining assignment '+ast.unparse(n))
                    seen.update(assigned);picked.append(copy.deepcopy(n))
            elif isinstance(n,ast.If) and ast.unparse(n.test) in branches:
                visit(n.body if branches[ast.unparse(n.test)] else n.orelse)
            elif isinstance(n,ast.If):
                assigned={t.id for child in ast.walk(n) if isinstance(child,ast.Assign)
                    for dst in child.targets for t in ast.walk(dst) if isinstance(t,ast.Name)}
                if assigned & names:raise ValueError('Unselected defining branch '+ast.unparse(n.test))
    visit(tree.body)
    if seen!=names:raise ValueError('Missing defining source assignments '+str(names-seen))
    allowed={'dict':dict,'list':list,'sum':sum,'range':range,'min':min,'len':len,'enumerate':enumerate}
    scope={'__builtins__':allowed,**env}
    exec(compile(ast.fix_missing_locations(ast.Module(body=picked,type_ignores=[])),
        '<exact mathematical defining slice '+fn.__name__+'>','exec'),scope)
    return scope,dict(function_binding=current.ast_binding(fn),
        actual_defining_assignments=[ast.unparse(n) for n in picked],selected_branches=branches,
        exact_scalar_primitives_not_enclosure_values=True)


class OriginalLeadingSourceReplay:
    """Replay original raw histories and their unchanged common normalization."""
    def __init__(self,field):
        self.field=field;self.q=leading_check.Interpreter(field.leading)
        self.z,self.logP,self.mu=self.q.z,self.q.logP,self.q.mu
        self.flow=ExactPointFlow(self.logP)
        self.c=SimpleNamespace(mpf=lambda q:s.Rational(str(q)),exp=s.exp,ln=s.log)
        self.ref=SimpleNamespace(z=ExactTaylor(None,[self.z,1,0,0]),
            zrows=[self.z,1,0,0],Vref=[4*self.z,4,0,0])
        self.P0=[s.Function('original_P0')(self.z).diff(self.z,j)/math.factorial(j) for j in range(4)]
        self.op=SimpleNamespace(flow=self.flow,c=self.c,reference=self.ref,zrows=self.ref.zrows,
            P0=self.P0,Pstar=s.exp(self.logP))
        self.raw={};self.normalized={};self.bindings={}
        self.build()

    def endpoint(self,raw,x,t):
        return {key:[s.sympify(v).subs(x,t) for v in value] for key,value in raw.items()}

    def normalize(self,E,V,H):
        recover=leading.providers.o3.reference.long.generic.recover_inputs
        raw=dict(velocity=dict(theta=[E],axial=[V]),histories={key:[value] for key,value in H.items()})
        env=dict(op=self.op,packet={'raw_current_radius_y_derivative_axial_coefficients':raw,
            'original_P0_normalized_axial5':self.op.P0},f=self.flow,c=self.c,add=self.flow.add)
        out,binding=defining_slice(recover,('scale','invS','S','raw','E','V','histories','m','h','k','e','p','pressure'),env)
        return dict(E=out['E'],V=out['V'],histories=out['histories'],
            pressure=out['pressure'],P0=self.op.P0),binding

    def turnoff(self,y,phase):
        fn=leading.providers.slope.original.original.turnoff_derivatives
        node=ast.parse(textwrap.dedent(inspect.getsource(fn))).body[0]
        u=s.Symbol('exact_sigma_argument',real=True)
        def sigma_jet(ctx,value):
            return [s.diff(leading_check.Sigma(u),u,j).subs(u,value)/math.factorial(j) for j in range(5)]
        scope={**fn.__globals__,'sigma_jets':sigma_jet}
        exec(compile(ast.fix_missing_locations(ast.Module(body=[node],type_ignores=[])),
            '<unchanged original physical-y turnoff derivatives>','exec'),scope)
        return scope[fn.__name__](self.c,y,40,phase)

    def build(self):
        charts=self.field.leading.functions['charts'];providers=leading.providers
        for name,chart in charts.items():
            f,c=self.flow,self.c;x=self.q.at(chart['native_coordinate']);y=self.q.at(chart['physical_coordinate'])
            env=dict(op=self.op,ref=self.ref,f=f,c=c,x=x,y=y,t=x,selector=x,chart=name.split('_',1)[-1])
            stages=[]
            if name=='Rh_reference':
                env.update(unit=lambda v:[v,0,0,0],zero=[0]*6,
                    reference=SimpleNamespace(background=SimpleNamespace(logarithm=lambda q:q.log())))
                src,proof=defining_slice(providers.rh.original.background_cell,('shapes','logE'),env);stages.append(proof)
                src['source']=dict(actual_normalized_six_history_shapes=src['shapes'],original_log_E_axial5=src['logE'],
                    actual_velocity_V=self.ref.Vref,actual_velocity_V_y=[0]*6)
                src['long']=SimpleNamespace(IntervalTaylor=ExactTaylor)
                out,proof=defining_slice(providers.rh.reference.recover_cell,
                    ('shapes','logu','relative','E','E2','V','Vy','invP2','hist'),src);stages.append(proof)
            elif name=='O2_slope':
                K=chart['kernels'];env.update(J=self.q.at(K['J']),mass=[self.q.at(K[k]) for k in ('theta','pressure','energy')])
                out,proof=defining_slice(providers.slope.original.background_cell,
                    ('factor','qi','qi2','E','V','zero','h','hist'),env);stages.append(proof)
            elif name in ('O2_axial','O2_buffer'):
                parent=charts['O2_slope'];u=self.q.at(parent['native_coordinate'])
                old=self.endpoint(self.raw['O2_slope']['histories'],u,1)
                E1=[v.subs(u,1) for v in self.raw['O2_slope']['E']]
                env['parent']=dict(original_closed_O2_slope_background={'raw':{'raw_current_radius_y_derivative_axial_coefficients':
                    dict(velocity={'theta':[E1]},histories={key:[value] for key,value in old.items()})}})
                env.update(B=self.turnoff(y,x) if name=='O2_axial' else [0]*5,
                    K={key:self.q.at(chart['kernels'][key]) for key in ('B_mass','B_squared_mass')})
                out,proof=defining_slice(providers.axial.original.background_cell,
                    ('raw','u1','old','t','formal','d','root','d3','E','E2','z','z2','V','Vy','zero','hist'),env);stages.append(proof)
            else:
                pname='O2_buffer' if name=='O3_transition' else 'O3_transition';parent=charts[pname]
                u=self.q.at(parent['native_coordinate']);end=11 if pname=='O2_buffer' else 1
                old=self.endpoint(self.raw[pname]['histories'],u,end);E1=[v.subs(u,end) for v in self.raw[pname]['E']]
                key='original_closed_O2_axial_buffer_background' if name=='O3_transition' else 'original_closed_O3_background'
                env['parent']={key:{'raw':{'raw_current_radius_y_derivative_axial_coefficients':
                    dict(velocity={'theta':[E1]},histories={k:[value] for k,value in old.items()})}}}
                K=charts['O3_transition']['kernels']
                exactK={key:self.q.at(value) for key,value in K.items()}
                env.update(logmu=s.log(self.mu),transition_kernel_hull=lambda *args:exactK,
                    previous=SimpleNamespace(previous=SimpleNamespace(original=SimpleNamespace(sigma_jets=lambda *args:[leading_check.Sigma(x)]))),
                    incoming=SimpleNamespace(weighted=SimpleNamespace(terminal=SimpleNamespace(local=SimpleNamespace(
                        packets=SimpleNamespace(recovery=SimpleNamespace(exp_average=lambda ctx,v:leading_check.Exprel(v))))))),
                    positive_decay_mass=lambda ctx,k,t:t*leading_check.Exprel(-k*t))
                names=('key','raw','u1','old','mu','mu_cover','zero','d','d3','K','sig','factor','theta','energy','pressure','E','E2','V','Vy','hist')
                if name=='O3_power':names=names+('average',)
                out,proof=defining_slice(providers.o3.original.background_cell,names,env,
                    {"chart == 'transition'":name=='O3_transition'});stages.append(proof)
            self.raw[name]=dict(E=out['E'],V=out['V'],histories=out['hist'])
            if name in ('O2_axial','O2_buffer'):self.raw[name].update(Vy=out['Vy'],B=out['B'])
            common,normproof=self.normalize(out['E'],out['V'],out['hist'])
            self.normalized[name]=common;self.bindings[name]=dict(original_source_replay=stages,unchanged_normalization=normproof)


def build_aliases(field,nodes=None):
    g=field.graph;charts=field.leading.functions['charts'];aliases={};quantities={'E','V'}|{'history_'+k for k in current.RATES}
    nodes=field.prefix if nodes is None else nodes
    for i,n in enumerate(nodes):
        if n['operation']!='current_original_leading_function_recipe' or n.get('native_chart') not in charts or n.get('quantity') not in quantities:continue
        chart=charts[n['native_chart']];quantity=n['quantity'];j=n['Z_order']
        if n['source_family']!=field.identity or not n['defining_quantity_not_a_range_value']:raise ValueError('Actual original leading recipe required')
        if not 0<=j<=3:raise ValueError('Ordinary leading order zero through three required')
        packet=chart[quantity] if quantity in ('E','V') else chart['histories'][quantity[8:]]
        ref=target.rows(packet)[j];c=nodes[n['coordinate']]
        if c['operation']=='exact_rational':coordinate=g.constant(s.Rational(c['numerator'],c['denominator']))
        elif c['operation']=='bound_variable':coordinate=g.symbol(c['name'])
        else:raise ValueError('Unbound original source recipe coordinate '+str(c))
        value=leading.scalar_at(g,ref,chart['native_coordinate'],coordinate)
        aliases[str(i)]=dict(original_leaf=i,explicit_original_source_row=value.node,native_chart=n['native_chart'],
            quantity=quantity,ordinary_Z_order=j,actual_source_recipe=n['recipe'],
            source_graph='accepted_outer' if nodes is field.prefix else 'accepted_repair_band',
            resolved_graph='canonical_bridge_graph',same_original_source_family_and_common_units=True)
    if not aliases:raise ValueError('Actual original leading source leaf aliases required')
    return aliases


def build_interface(field):
    g=field.graph;alg=target.C3Algebra(g);b=field.leading.bridge
    ref=lambda i:source.FunctionRef(g,i);jet=lambda q:target.C3Function(*map(ref,q))
    source_endpoint=field.leading.functions['actual_source_Rc_leading_C3']
    invN=g.quotient(g.one,ref(b.outer_N),'same original positive selected repair integer')
    correction={key:alg.scale(jet(q),invN) for key,q in b.outer_terminal.items()}
    complete={key:alg.add(source_endpoint[key],correction[key]) for key in current.RATES}
    rh=outer.rh.read(outer.rh.NAME)['actual_Rh_C3_continuation_functions']
    P0=jet(rh['actual_independent_P0_C3']);c=field.leading.functions['charts']['O3_power']
    amplitude=outer.substitute(g,c['E'],c['native_coordinate'],g.constant(2))
    original=b.outer_report['actual_five_chart_C3_continuation_functions'];identified={}
    for name,packet in original.items():
        explicit=field.leading.functions['charts'][name];endpoint=ref(packet['endpoint'])
        identified[name]={key:outer.substitute(g,q,explicit['native_coordinate'],endpoint)
            for key,q in explicit['histories'].items()}
    return dict(actual_outer_leading_Duhamel_as_original_source_C3=identified,
        actual_outer_Rc_leading_C3=source_endpoint,actual_outer_Rc_correction_C3=correction,
        actual_outer_Rc_complete_C3=complete,actual_outer_Rc_amplitude_C3=amplitude,
        actual_independent_P0_C3=P0,actual_outer_Rc_absolute_pressure_C3=alg.add(P0,complete['p']),
        same_original_N_divided_once=True,complete_history_not_relative_zero=True)


class CurrentOuterLeadingC3Bridge:
    def __init__(self,require_checked=True):
        self.leading=leading.CurrentOuterLeadingC3Identity();self.identity=self.leading.identity
        self.graph=self.leading.graph;self.prefix=copy.deepcopy(self.graph.nodes)
        self.hashes={**self.leading.hashes,leading.NAME:sha(leading.NAME),leading.RECEIPT:sha(leading.RECEIPT),
            Path(leading_check.__file__).name:sha(Path(leading_check.__file__).name),Path(__file__).name:sha(Path(__file__).name)}
        self.aliases=build_aliases(self)
        self.band_aliases=build_aliases(self,self.leading.bridge.band_nodes)
        self.functions=build_interface(self);self.acceptance_loaded=False
        if require_checked:
            receipt=outer.rh.read(RECEIPT)
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES):raise ValueError('Checked actual leading provider resolver required')
            for name,value in receipt['input_hashes'].items():
                if sha(name)!=value:raise ValueError('Changed original leading source bridge '+name)
            self.acceptance_loaded=True

    def resolve_leading_leaf(self,node,consumer='outer'):
        """Return a canonical bridge-graph ref; consumer selects the source graph.

        Both graphs are re-expressed using their source recipe coordinates and
        the same original normalized mathematical functions. Returned node IDs
        always belong to self.graph, never to the caller's source graph. This is
        an exact function map, not a live numerical flow/ledger conversion.
        """
        key=str(node.node if hasattr(node,'node') else node)
        aliases=self.aliases if consumer=='outer' else self.band_aliases if consumer=='band' else None
        if aliases is None or key not in aliases:raise ValueError('Unadmitted source recipe; use its original provider')
        return source.FunctionRef(self.graph,aliases[key]['explicit_original_source_row'])

    def complete_Rc_functions(self):return self.functions['actual_outer_Rc_complete_C3']


def run():
    began=time.monotonic();field=CurrentOuterLeadingC3Bridge(require_checked=False)
    report=dict(candidate_actual_original_leading_C3_provider_bridge_constructed=True,source_family=field.identity,
        exact_graph_nodes=field.graph.nodes,accepted_original_leading_graph_prefix_length=len(field.prefix),
        actual_original_leading_provider_aliases=field.aliases,
        actual_band_leading_provider_aliases=field.band_aliases,
        resolved_function_graph_contract=dict(returned_refs_owner='canonical_bridge_graph',
            consumer_selects_original_source_graph_not_returned_owner=True,
            same_original_source_family_and_common_V_m_k_Pstar_normalization=True,
            numerical_flow_ledger_or_point_oracle_conversion=False),
        actual_original_leading_and_complete_Rc_functions=target.encoded(field.functions),
        source_replay_binding=current.ast_binding(OriginalLeadingSourceReplay),
        selected_native_pulse_constructor_consumes_current_C3_frame=False,
        current_numeric_point_field_oracle_installed=False,
        global_physical_time_Cartesian_heat_cone_and_temporal_recursion_installed=False,
        **dict.fromkeys(GATES+OPEN,False),input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(report,separators=(',',':'))+'\n').encode(),mtime=0))
    print('Original leading C3 source recipe resolver constructed',flush=True)
    return field


if __name__=='__main__':run()
