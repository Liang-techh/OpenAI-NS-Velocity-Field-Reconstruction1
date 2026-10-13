"""True ordinary Z2/Z3 functions of the original future and pulse branch.

Accepted value/first-Z handles stay unchanged. Native derivatives come
from the checked C3 frame. Gamma derivatives retain the positive improper
integrals and actual nonzero S; numerical ranges are never differentiated.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_complete_future_functions as complete

parent, source, controls = complete.parent, complete.source, complete.controls
c3 = parent.active.c3
HERE, PREFIX, sha = complete.HERE, complete.PREFIX, complete.sha
NAME = PREFIX+'current_original_n1_pulse_C3_functions.json.gz'
RECEIPT = PREFIX+'current_original_n1_pulse_C3_functions_check.json'
GATE = 'current_original_n1_actual_future_and_pulse_C3_functions_installed'
OPEN = complete.OPEN


def encode(value):
    if type(value) is c3.C3Function:
        return [q.node for q in (value.value,value.Z,value.ZZ,value.ZZZ)]
    if isinstance(value,dict):return {key:encode(item) for key,item in value.items()}
    if isinstance(value,(tuple,list)):return [encode(item) for item in value]
    return parent.encode(value)


def Gamma_higher_functions(g,owner):
    """Exact Dhat_ZZ/ZZZ; no divided expression singular at Z=+/-1."""
    a=owner.pieces['a'];S=owner.schedule['S']
    t,v=owner.pieces['offset_variable'],owner.pieces['Gamma_variable']
    xi=owner.pieces['xi'];z=g.symbol('Z');c=g.constant
    normal=g.unary('gamma',g.add(g.one,a))
    log=g.unary('log1p',g.mul(xi,v))
    moments=[]
    for k in (1,2,3):
        ak=g.add(a,c(k))
        body=g.unary('exp',g.add(g.neg(v),g.mul(ak,g.unary('log',v)),g.neg(g.mul(ak,log))))
        integral=controls.integral(g,body,g.nodes[v.node]['name'],g.zero,owner.pieces['infinity'],
            exact_positive_Gamma_expectation=True,full_improper_integral=True,
            Gamma_normalization=normal.node,Gamma_moment_order=k,
            zero_endpoint_integrand_continuation=0)
        moments.append(g.quotient(integral,normal,'positive Gamma(1+a)'))
    M1,M2,M3=moments
    e1,e2,e3=[g.unary('exp',g.mul(c(-k),t)) for k in (1,2,3)]
    ap1,ap2=g.add(a,g.one),g.add(a,c(2))
    second=g.add(g.mul(c(-4),e1,M1),g.mul(c(-16),ap1,z,z,S,e2,M2))
    third=g.add(g.mul(c(-48),ap1,z,S,e2,M2),
        g.mul(c(-64),ap1,ap2,z,z,z,S,S,e3,M3))
    return c3.C3Function(owner.functions['Dhat'].value,owner.functions['Dhat'].Z,second,third),dict(
        positive_Gamma_moments=moments,normalization=normal,
        exact_nonzero_S=S,offset_variable=t)


class OrdinaryZDerivative:
    """Exact DAG chain/product rules, fail closed on unknown sources."""
    def __init__(self,g):
        self.g=g;self.z=g.symbol('Z')
        self.derivatives={g.zero.node:g.zero,g.one.node:g.zero,self.z.node:g.one}

    def seed(self,left,right):
        self.g.ids([left,right])
        prior=self.derivatives.get(left.node)
        if prior is not None and prior!=right:
            raise ValueError('Conflicting actual derivative seed')
        self.derivatives[left.node]=right

    def seed_pair(self,pair):
        if type(pair) is not source.C1Function:raise TypeError('Genuine C1 seed required')
        self.seed(pair.value,pair.Z)

    def seed_jet(self,jet):
        if type(jet) is not c3.C3Function:raise TypeError('Genuine C3 seed required')
        rows=(jet.value,jet.Z,jet.ZZ,jet.ZZZ)
        for left,right in zip(rows,rows[1:]):self.seed(left,right)

    def operands(self,i):
        row=self.g.nodes[i];op=row['operation']
        if op in ('exact_rational','mathematical_positive_infinity','mathematical_pi','bound_variable'):
            return []
        if op in ('current_original_source_parameter','original_source_parameter'):
            if not row.get('Z_independent'):raise ValueError('Unproved scalar source independence')
            return []
        if op in ('sum','product'):return row['arguments']
        if op in ('negative','analytic_unary','exact_original_source_flat_sigma',
                  'exact_original_source_heat_phi','compact_raw_beta'):return [row['argument']]
        if op=='positive_quotient':return [row['numerator'],row['denominator']]
        if op=='definite_integral':
            if row['variable']=='Z':raise ValueError('Axial bound variable requires separate capture handling')
            return [row['integrand'],row['lower'],row['upper']]
        if op=='function_substitution':
            return ([row['value']] if row['variable']==self.z.node else [row['expression'],row['value']])
        raise ValueError('Actual derivative source not seeded: '+op)

    def derivative(self,value):
        self.g.ids([value]);stack=[(value.node,False)];visiting=set()
        while stack:
            i,ready=stack.pop()
            if i in self.derivatives:continue
            if not ready:
                if i in visiting:raise ValueError('Cyclic defining derivative graph')
                visiting.add(i);stack.append((i,True))
                stack.extend((j,False) for j in reversed(self.operands(i)) if j not in self.derivatives)
            else:
                self.derivatives[i]=self.rule(i);visiting.remove(i)
        return self.derivatives[value.node]

    def rule(self,i):
        g=self.g;row=g.nodes[i];op=row['operation'];ref=lambda j:source.FunctionRef(g,j)
        d=lambda j:self.derivatives[j];c=g.constant
        if not self.operands(i):return g.zero
        if op=='sum':return g.add(*(d(j) for j in row['arguments']))
        if op=='product':
            args=row['arguments']
            return g.add(*(g.mul(d(j),*(ref(k) for n,k in enumerate(args) if n!=p))
                for p,j in enumerate(args) if d(j)!=g.zero))
        if op=='negative':return g.neg(d(row['argument']))
        if op=='positive_quotient':
            n,v=row['numerator'],row['denominator']
            return g.quotient(g.sub(g.mul(d(n),ref(v)),g.mul(ref(n),d(v))),
                g.mul(ref(v),ref(v)),'square of the same checked nonzero denominator')
        if op in ('analytic_unary','exact_original_source_flat_sigma','exact_original_source_heat_phi','compact_raw_beta'):
            arg=row['argument'];dx=d(arg)
            if dx==g.zero:return g.zero
            if op!='analytic_unary':raise ValueError('Nonconstant axial flat/bump derivative must be source-bound')
            x=ref(arg);name=row['name']
            if name=='exp':return g.mul(ref(i),dx)
            if name=='log':return g.quotient(dx,x,'same positive original log argument')
            if name=='log1p':return g.quotient(dx,g.add(g.one,x),'same positive 1+x')
            if name in ('sqrt','positive_sqrt'):
                return g.quotient(dx,g.mul(c(2),ref(i)),'same positive selected square root')
            raise ValueError('Nonconstant analytic derivative not defined: '+name)
        if op=='definite_integral':
            if d(row['lower'])!=g.zero or d(row['upper'])!=g.zero:
                raise ValueError('Axial moving endpoint needs full Leibniz law')
            body=d(row['integrand'])
            if body==g.zero:return g.zero
            attrs=copy.deepcopy(row);attrs.pop('operation');attrs['integrand']=body.node
            attrs['ordinary_Z_derivative_order']=attrs.get('ordinary_Z_derivative_order',0)+1
            attrs['derivative_of_same_defining_integral']=True
            attrs['differentiation_under_integral']=True
            attrs['numerical_enclosure_not_differentiated']=True
            return g.node('definite_integral',**attrs)
        if op=='function_substitution':
            if d(row['value'])!=g.zero:raise ValueError('Axial-dependent substitution needs full chain rule')
            if row['variable']==self.z.node:return g.zero
            attrs=copy.deepcopy(row);attrs.pop('operation');attrs['expression']=d(row['expression']).node
            return g.node('function_substitution',**attrs)
        raise ValueError('No ordinary derivative rule for '+op)

    def lift(self,pair):
        self.seed_pair(pair)
        zz=self.derivative(pair.Z);zzz=self.derivative(zz)
        return c3.C3Function(pair.value,pair.Z,zz,zzz)


def native_C3_functions(g):
    raw=parent.linear.read(parent.quiet.FRAME_NAME)
    original=parent.linear.realization.restore_graph(raw['exact_graph_nodes'])
    importer=parent.GraphImport(g,original)
    native=raw['actual_Rp_native_frame_functions']['actual_identified_Rp_native_frame_C3']
    ref=lambda i:importer(source.FunctionRef(original,i))
    return {key:c3.C3Function(*map(ref,native[key])) for key in ('u','m1','m2','energy','X','Mp','P0','pressure')}


class CurrentOriginalN1PulseC3Functions:
    def __init__(self,complete_owner=None,require_checked=True):
        self.complete=complete_owner or complete.CurrentOriginalN1CompleteFutureFunctions()
        if not self.complete.acceptance_loaded:raise ValueError('Checked actual complete future C1 owner required')
        self.complete.assert_graph();self.g=self.complete.g;self.family=copy.deepcopy(self.complete.family)
        self.hashes=dict(self.complete.hashes)
        receipt=parent.linear.read(complete.RECEIPT)
        for name,digest in receipt['input_hashes'].items():
            if sha(name)!=digest or (name in self.hashes and self.hashes[name]!=digest):
                raise ValueError('Changed complete future ancestor '+name)
            self.hashes[name]=digest
        self.hashes[complete.NAME]=sha(complete.NAME);self.hashes[complete.RECEIPT]=sha(complete.RECEIPT)
        self.native=native_C3_functions(self.g)
        self.Gamma_Dhat,self.Gamma_evidence=Gamma_higher_functions(self.g,self.complete.Gamma)
        self.differentiator=OrdinaryZDerivative(self.g)
        for jet in self.native.values():self.differentiator.seed_jet(jet)
        self.differentiator.seed_jet(self.Gamma_Dhat)
        groups=dict(heat={key:value for key,value in self.complete.Gamma.functions.items() if key!='Dhat'},
            histories={key:value for key,value in self.complete.Gamma.histories.items() if type(value) is source.C1Function},
            angular=self.complete.angular,future=self.complete.functions,pulse=self.complete.pulse_functions)
        for group in groups.values():
            for pair in group.values():self.differentiator.seed_pair(pair)
        self.functions={name:{key:self.differentiator.lift(pair) for key,pair in group.items()}
            for name,group in groups.items()}
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(self.g.nodes);self.snapshot=self._snapshot();self.acceptance_loaded=False
        if require_checked:
            receipt,raw=parent.linear.read(RECEIPT),parent.linear.read(NAME)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family
                    or any(receipt.get(key) for key in OPEN) or raw['exact_graph_nodes']!=self.g.nodes
                    or raw['function_bindings']!=self.snapshot):
                raise ValueError('Checked actual C3 pulse/future functions required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed C3 pulse/future source '+name)
            self.acceptance_loaded=True

    def _snapshot(self):
        return encode(dict(functions=self.functions,native=self.native,Gamma_Dhat=self.Gamma_Dhat,
            Gamma_evidence=self.Gamma_evidence,input_hashes=self.hashes))

    def assert_graph(self):
        self.complete.assert_graph()
        if self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot:
            raise ValueError('Actual C3 pulse/future function bodies changed')


def run():
    began=time.monotonic();owner=CurrentOriginalN1PulseC3Functions(require_checked=False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,function_bindings=owner.snapshot,
        accepted_value_and_first_Z_handles_preserved=True,actual_Z2_Z3_not_zero_padded=True,
        exact_improper_Gamma_integrals_retained=True,source_ranges_not_differentiated=True,
        **dict.fromkeys((GATE,)+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_ACTUAL_PULSE_C3_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw


if __name__=='__main__':run()
