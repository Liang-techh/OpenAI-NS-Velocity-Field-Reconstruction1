"""Node-aware genuine reference/O2 evaluation and unclosed affine histories.

The unchanged original graph is evaluated with issued integral callbacks.
True source phase is derived by its owner, never rounded from raw radii.
Unknown preceding histories remain explicit boundary functions. This is a
partial original evaluator and does not satisfy the full control oracle.
"""
from dataclasses import dataclass
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_O2_integral_callback as O2

reference=O2.reference;leaves=O2.leaves;precise=O2.precise;transport=O2.transport
HERE,PREFIX,sha,ep=O2.HERE,O2.PREFIX,O2.sha,O2.ep
NAME=PREFIX+'current_original_route_function_evaluator.json'
RECEIPT=PREFIX+'current_original_route_function_evaluator_check.json'
GATE='actual_original_accepted_six_route_node_aware_evaluator_and_affine_history_operator_installed'


@dataclass(frozen=True)
class OriginalLocalHistoryFrame:
    Z: object
    N: int
    forcing: dict
    decay: dict
    record: dict


class OriginalRouteFunctionEvaluator:
    mode='original_accepted_route_partial'
    def __init__(self,*,reference_count=4,O2_count=64,reference_bits=32,O2_bits=24):
        self.provider=leaves.OriginalPointSourceLeaves();self.built=self.provider.built
        self.reference=reference.OriginalReferenceIntegralCallback(self.provider)
        self.O2=O2.OriginalO2IntegralCallback(self.provider);self.c=self.O2.c
        if type(reference_count) is not int or not 1<=reference_count<=256:
            raise ValueError('Reference partition in[1,256] required')
        if type(O2_count) is not int or O2_count not in self.O2.cells.parent.parent.levels:
            raise ValueError('Accepted O2 source partition required')
        if type(reference_bits) is not int or not 4<=reference_bits<=256 or type(O2_bits) is not int or not 8<=O2_bits<=256:
            raise ValueError('Accepted inverse-bit settings required')
        self.settings=dict(reference_count=reference_count,O2_count=O2_count,
            reference_bits=reference_bits,O2_bits=O2_bits)
        self.family=self.source_family=self.provider.family
        self.source_graph_sha256=self.provider.source_graph_sha256
        self.hashes=dict(self.provider.hashes);self.cache={};self.frames={};self.issued={}
        for owner in (self.reference,self.O2):self.bind_hashes(owner.hashes)
        for module in (reference,O2):
            receipt=json.loads((HERE/module.RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(module.GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted same-family graph integral prerequisite required')
            self.bind_hashes({**receipt['input_hashes'],module.RECEIPT:sha(module.RECEIPT)})
        self.integrals={node:self.reference for node in self.reference.roles}
        self.integrals.update({node:self.O2 for node in self.O2.roles})
        if len(self.integrals)!=60:raise ValueError('Exactly60 genuine accepted-route integral nodes required')
        saved=json.loads((HERE/transport.NAME).read_bytes())
        self.route=[row for row in saved['exact_original_cells'] if row['chart'] in leaves.SUPPORTED]
        positions=[i for i,row in enumerate(saved['exact_original_cells']) if row['chart'] in leaves.SUPPORTED]
        if len(self.route)!=6 or positions!=list(range(positions[0],positions[0]+6)):
            raise ValueError('Six contiguous unchanged original route cells required')
        if self.route[0]['chart']!='Rh_reference' or any(row['chart']!='O2_slope' for row in self.route[1:]):
            raise ValueError('Reference then five O2 cells required')
        self.boundaries={};self.boundary_symbols={};self.history_symbols={}
        for key,pair in self.route[0]['incoming'].items():
            for jet,label in (('C0','value'),('Z','Z')):
                node=pair[label];self.boundaries[node]=(key,jet)
                self.boundary_symbols[node]=sy.Symbol('unknown_original_incoming_'+key+'_'+jet,real=True)
        if len(self.boundaries)!=10:raise ValueError('Ten distinct original upstream C0/Z functions required')
        self.integral_symbols={node:sy.Symbol('issued_integral_'+str(node),real=True) for node in self.integrals}
        self.identities=self.affine_identities()
        self.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})

    def bind_hashes(self,closure):
        for name,digest in closure.items():
            if sha(name)!=digest:raise ValueError('Original evaluator dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original evaluator closures disagree: '+name)
            self.hashes[name]=digest

    def unchanged(self):
        if leaves.digest_rows(self.built['graph'].nodes)!=self.provider.graph_digest:
            raise ValueError('Unchanged issued original function graph required')

    def ordinary(self,value):
        if not hasattr(value,'_mpi_') or hasattr(value,'scale') or isinstance(value,dict):
            raise TypeError('Directed ordinary interval required; no cap, scalar selection or formal factors')
        # Endpoints are converted outward into one integration/evaluation
        # context. No local formal basis is transported across owners.
        result=self.c.mpf(ep(value))
        if result.ctx is not self.c:raise TypeError('One directed numeric evaluation context required')
        return result

    def request(self,Z,N):
        z=leaves.base.point.pressure.exact_Z(Z)
        if z==0:raise ValueError('Accepted integral/history evaluator requires fixed nonzero Z')
        if type(N) is not int or N<160 or N.bit_length()>4096:
            raise ValueError('Explicit common integer N>=160, at most4096 bits required')
        self.unchanged();return z

    def parameter(self,name):
        """Original source recipes; tiny widths retain a nonzero upper tail."""
        self.unchanged();phase=self.provider.phase;b=phase.binder;c=self.c
        if name not in self.built['parameters']:raise ValueError('Issued original source parameter required')
        with mp.workdps(c.dps+40):
            if name in ('logP','Tw','Md','sc'):value=phase.analytic(phase.symbols[name],c)
            elif name in ('logC','T'):value=c.mpf(b.fixed[name])
            elif name in ('hbB','hbS'):
                value=O2.positive.small_exp(c,c.mpf(b.loghB if name=='hbB' else b.loghS))
            elif name=='mu':
                logmu=c.ln(c.mpf(1)/1000)-4*phase.analytic(phase.symbols['logP'],c)
                value=O2.positive.small_exp(c,logmu)
            else:raise NotImplementedError('Original parameter recipe unavailable')
            return self.ordinary(value)

    def integrate(self,row,*,Z,N):
        z=self.request(Z,N);g=self.built['graph']
        nodes=[i for i,item in enumerate(g.nodes) if row is item]
        if len(nodes)!=1 or nodes[0] not in self.integrals:
            raise NotImplementedError('Only issued reference/O2 integral nodes are installed')
        owner=self.integrals[nodes[0]];s=self.settings
        options=dict(count=s['reference_count'],bits=s['reference_bits']) if owner is self.reference else dict(count=s['O2_count'],bits=s['O2_bits'])
        with mp.workdps(self.c.dps+40):return self.ordinary(owner.integrate(row,Z=str(z),N=N,**options))

    def source(self,row,*,coordinate,Z,N,phase=None):
        if phase is not None:raise ValueError('Native source owner derives phase; arbitrary supplied phase is forbidden')
        self.request(Z,N);self.provider.require_row(row)
        with mp.workdps(self.c.dps+40):
            return self.ordinary(self.provider.source(row,coordinate=coordinate,Z=Z,N=N))

    def symbolic(self,index,memo=None):
        """Collect original arithmetic with ten true upstream boundary atoms."""
        memo={} if memo is None else memo
        if index in memo:return memo[index]
        if index in self.boundary_symbols:return self.boundary_symbols[index]
        if index in self.integral_symbols:return self.integral_symbols[index]
        row=self.built['graph'].nodes[index];op=row['operation'];child=lambda i:self.symbolic(i,memo)
        if op=='exact_rational':value=sy.Rational(row['numerator'],row['denominator'])
        elif op=='original_source_parameter':value=sy.Symbol(row['name'],real=True)
        elif op=='shared_positive_integer':value=sy.Symbol('N',positive=True,integer=True)
        elif op=='bound_variable':value=sy.Symbol(row['name'],real=True)
        elif op=='sum':value=sy.Add(*(child(i) for i in row['arguments']))
        elif op=='product':value=sy.Mul(*(child(i) for i in row['arguments']))
        elif op=='negative':value=-child(row['argument'])
        elif op=='positive_quotient':value=child(row['numerator'])/child(row['denominator'])
        elif op=='analytic_unary' and row['name'] in ('exp','log'):
            value=(sy.exp if row['name']=='exp' else sy.log)(child(row['argument']))
        else:raise NotImplementedError('Unsupported source is not replaced at affine boundary: '+op)
        memo[index]=value;return value

    def affine_identities(self):
        """Verify actual outgoing roots, not a newly invented history DAG."""
        records=[];g=self.built['graph'];memo={}
        for route_index,cell in enumerate(self.route):
            right,_=reference.symbolic_node(self.provider,cell['upper'])
            width=sy.Integer(5) if route_index==0 else sy.Integer(5)+right
            if not width.is_Rational:raise ValueError('Exact accepted-route relative width required')
            for key,pair in cell['outgoing'].items():
                rate=transport.RATES[key];r=sy.Rational(rate.numerator,rate.denominator)
                for jet,label in (('C0','value'),('Z','Z')):
                    boundary=self.route[0]['incoming'][key][label];node=pair[label]
                    expected=sy.exp(-r*width)*self.boundary_symbols[boundary]
                    for i,previous in enumerate(self.route[:route_index+1]):
                        endpoint=sy.Integer(0) if i==0 else reference.symbolic_node(self.provider,previous['upper'])[0]
                        current_endpoint=sy.Integer(0) if route_index==0 else right
                        integral=previous['contributions'][key][label]
                        expected+=sy.exp(-r*(current_endpoint-endpoint))*self.integral_symbols[integral]
                    actual=self.symbolic(node,memo)
                    if sy.simplify(actual-expected)!=0:raise ValueError('Original affine outgoing/history identity differs')
                    self.history_symbols[route_index,key,jet]=(node,width)
                    records.append(dict(route_index=route_index,route_label=cell['label'],key=key,ordinary_order=jet,
                        actual_outgoing_node=node,unknown_original_upstream_node=boundary,exact_relative_width=str(width),
                        exact_original_affine_history_identity_passed=True,
                        incoming_not_reset_or_selected=True,ordinary_Z_geometry_has_no_boundary_terms=True))
        return records

    def walk(self,index,variables,z,N,*,forcing=False,cache=None):
        cache={} if cache is None else cache
        key=(index,tuple(sorted((name,str(q)) for name,q in variables.items())))
        if key in cache:return cache[key]
        c=self.c;row=self.built['graph'].nodes[index];op=row['operation']
        child=lambda i:self.walk(i,variables,z,N,forcing=forcing,cache=cache)
        if forcing and index in self.boundaries:
            # Evaluate the additive component of a proved affine operator.
            # This does not evaluate or assert zero original boundary data.
            value=c.mpf(0)
        elif op=='exact_rational':value=c.mpf(row['numerator'])/row['denominator']
        elif op=='bound_variable':
            if row['name'] not in variables:raise ValueError('Explicit exact source coordinate binding required')
            q=variables[row['name']];value=c.mpf(int(q.p))/int(q.q)
        elif op=='original_source_parameter':value=self.parameter(row['name'])
        elif op=='shared_positive_integer':value=c.mpf(N)
        elif op=='sum':value=sum((child(i) for i in row['arguments']),c.mpf(0))
        elif op=='product':
            value=c.mpf(1)
            for i in row['arguments']:value*=child(i)
        elif op=='negative':value=-child(row['argument'])
        elif op=='positive_quotient':
            denominator=child(row['denominator'])
            if ep(denominator)[0]<=0:raise ArithmeticError('Directed original positive denominator not proved')
            value=child(row['numerator'])/denominator
        elif op=='analytic_unary':
            name=row['name']
            if name=='fractional_part':raise NotImplementedError('Native phase must be derived inside its issued source callback')
            # Exact radius differences and widths are collected before
            # numeric arithmetic; huge common radii cannot erase a width.
            expr,_=reference.symbolic_node(self.provider,row['argument'])
            expr=sy.simplify(expr.subs({sy.Symbol(k,real=True):v for k,v in variables.items()}))
            if not expr.is_Rational:raise NotImplementedError('Only exact rational reduced analytic arguments are installed')
            arg=c.mpf(int(expr.p))/int(expr.q)
            if name=='exp':
                if ep(arg)[1]>ep(2*c.dps*c.ln(10))[0]:raise ArithmeticError('Unbounded ordinary exponential not materialized')
                value=O2.positive.small_exp(c,arg)
            elif name=='log':
                if ep(arg)[0]<=0:raise ArithmeticError('Positive log argument required')
                value=c.ln(arg)
            else:raise NotImplementedError('Original analytic operation unavailable')
        elif op=='original_function_graph':
            coordinate=grow=self.built['graph'].nodes[row['coordinate']]
            if coordinate['operation']!='bound_variable' or coordinate['name'] not in variables:
                raise ValueError('Exact native source coordinate binding required')
            value=self.source(row,coordinate=str(variables[grow['name']]),Z=str(z),N=N)
        elif op=='definite_integral':value=self.integrate(row,Z=str(z),N=N)
        else:raise NotImplementedError('Original function operation unavailable: '+op)
        value=self.ordinary(value);cache[key]=value;return value

    def evaluate(self,root,*,Z,N,variables=None):
        z=self.request(Z,N);self.built['graph'].ids([root])
        values={} if variables is None else {name:precise.rational(value) for name,value in variables.items()}
        key=(root.node,z,N,tuple(sorted((name,str(q)) for name,q in values.items())))
        if key not in self.cache:
            with mp.workdps(self.c.dps+40):self.cache[key]=self.walk(root.node,values,z,N)
        self.unchanged();return self.cache[key]

    def history_frame(self,*,Z,N):
        z=self.request(Z,N);key=(z,N)
        if key in self.frames:return self.frames[key]
        c=self.c;forcing={};decay={};cache={}
        with mp.workdps(c.dps+40):
            for role,(node,width) in self.history_symbols.items():
                route_index,density_key,jet=role;rate=transport.RATES[density_key]
                rr=c.mpf(rate.numerator)/rate.denominator
                decay[role]=c.exp(-rr*(c.mpf(int(width.p))/int(width.q)))
                forcing[role]=self.walk(node,{},z,N,forcing=True,cache=cache)
        record=dict(source_family=self.family,source_graph_sha256=self.source_graph_sha256,
            original_Z_exact=str(z),explicit_candidate_N=N,settings=self.settings,
            exact_original_affine_history_identities=self.identities,
            actual_six_route_forcing_and_homogeneous_intervals=[dict(route_index=i,key=k,ordinary_order=j,
                additive_original_integral_contribution=forcing[i,k,j],original_upstream_multiplier=decay[i,k,j])
                for i,k,j in forcing],
            exact_unknown_upstream_function_nodes={str(node):dict(key=k,ordinary_order=j) for node,(k,j) in self.boundaries.items()},
            unknown_incoming_remains_explicit_not_zeroed=True,separate_original_P0_and_P0_Z_unchanged=True,
            additive_operator_component_not_actual_full_history=True,
            original_Z_independent_geometry_no_extra_Z_boundary_terms=True,
            genuine_reference_and_O2_source_integral_errors_retained=True,
            full_scalar_control_evaluator_compatibility_installed=False,actual_original_upstream_history_installed=False)
        frame=OriginalLocalHistoryFrame(z,N,forcing,decay,record)
        self.frames[key]=frame;self.issued[id(frame)]=frame;self.unchanged();return frame

    def apply_assumed_boundary(self,frame,incoming):
        """Conditional operator action; no assertion of actual upstream data."""
        self.unchanged()
        if type(frame) is not OriginalLocalHistoryFrame or self.issued.get(id(frame)) is not frame:
            raise ValueError('Issued genuine affine history frame required')
        roles=set(self.boundaries.values())
        if set(incoming)!=roles:raise ValueError('All ten explicitly assumed upstream C0/Z enclosures required')
        if any(not hasattr(value,'_mpi_') or value.ctx is not self.c or hasattr(value,'scale') for value in incoming.values()):
            raise TypeError('Directed same-context assumed boundary enclosures required')
        with mp.workdps(self.c.dps+40):
            return {role:frame.decay[role]*incoming[role[1:]]+value for role,value in frame.forcing.items()}


@precise.phase.native.inlet.source_precision
def run(*,return_live=False):
    began=time.monotonic();owner=OriginalRouteFunctionEvaluator();g=owner.built['graph'];frames=[];integrals=[]
    for N in (160,257):
        frame=owner.history_frame(Z='.37',N=N);frames.append(frame)
        rows={str(node):owner.evaluate(transport.FunctionRef(g,node),Z='.37',N=N) for node in owner.integrals}
        integrals.append(dict(explicit_candidate_N=N,original_Z_exact='37/100',actual_issued_integral_values=rows))
        print('Original node-aware six-route affine history:',N,'60 issued integrals /60 outgoing operators',flush=True)
    points=[]
    for chart,q in (('Rh_reference','-2.337'),('O2_slope','.537')):
        for key in transport.RATES:
            for jet in ('C0','Z'):
                row=owner.provider.rows[chart,'density_'+key+'_'+jet];node=next(i for i,r in enumerate(g.nodes) if r is row)
                variable=g.nodes[row['coordinate']]['name']
                value=owner.evaluate(transport.FunctionRef(g,node),Z='.37',N=160,variables={variable:q})
                points.append(dict(node=node,chart=chart,coordinate=q,key=key,ordinary_order=jet,value=value))
    result=dict(**{GATE:True},source_family=owner.family,source_graph_sha256=owner.source_graph_sha256,mode=owner.mode,
        actual_original_integral_expression_frames=integrals,actual_original_source_expression_points=points,
        actual_original_affine_history_operator_frames=[frame.record for frame in frames],
        accepted_route_source_bound_phase_and_integral_dispatch_installed=True,
        actual_original_upstream_history_installed=False,full_scalar_control_evaluator_compatibility_installed=False,
        full_17_chart_source_or_24_cell_integral_oracle_installed=False,actual_five_controls_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(precise.phase.packets.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,
        scope='Genuine issued reference/O2 node-aware integral and source evaluation and exact original six-route affine history operators at fixed nonzero Z/current shared N. Unknown upstream functions and original separate P0 retained. No full history values, scalar control compatibility, whole-Z closure, global N, recursion or corrected NS.')
    (HERE/NAME).write_bytes(json.dumps(precise.encode(result),indent=2).encode()+b'\n')
    return (result,owner,frames) if return_live else result


if __name__=='__main__':run()
