"""Actual O3 power radial extension and own Rc mixed source join.

The O3 local loop is proved flat before differentiating its actual source.
All incoming background/correction/pressure histories are retained.
"""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_n1_combined_outer_mixed_source_functions as previous
import lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions as leading_source
import lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions as compact_source

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer=previous.outer
compact_high=previous.previous
high=compact_high.previous
radial=high.radial
NAME=PREFIX+'current_original_n1_Rc_mixed_source_join.json.gz'
RECEIPT=PREFIX+'current_original_n1_Rc_mixed_source_join_check.json'
GATES=('current_original_actual_O3_power_flat_loop_source_identity_installed',
    'current_original_actual_O3_power_y4_C5_Q4_C4_pressure_y3_C3_source_functions_installed',
    'current_original_actual_Rc_compact_inlet_history_pressure_source_identity_installed',
    'current_original_n1_Rc_own_y4_pressure_cumulative_y3_mixed_source_join_installed')
OPEN=tuple(dict.fromkeys(tuple(k for k in previous.OPEN if k not in GATES)+(
    'current_original_n1_Rin_rminus_own_mixed_source_joins_installed',)))


class CurrentOriginalN1RcMixedSourceJoin:
    def __init__(self,combined=None,require_checked=True):
        self.combined=combined if combined is not None else previous.CurrentOriginalN1CombinedOuterMixedSourceFunctions()
        if not self.combined.acceptance_loaded:raise ValueError('Accepted combined outer source required')
        self.combined.assert_graph();self.c=self.combined.compact;self.own=self.c.own_outer
        self.g=g=self.combined.g;self.family=copy.deepcopy(self.combined.family);self.hashes=dict(self.combined.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):
            self.hashes[name]=sha(name)
        raw=parent.linear.read(leading_source.NAME)
        if raw['source_family']!=self.family or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]:
            raise ValueError('Actual original outer source prefix required')
        data=raw['true_outer_C6_C5_bindings']
        self.leading=copy.deepcopy(dict(chart=data['charts']['O3_power'],background=data['background']['O3_power'],
            endpoint=data['endpoints']['O3_power'],background_endpoint=data['background_endpoints']['O3_power'],
            initial=data['initial'],mu_proof=data['mu_proof'],mu=data['mu'],eta=data['eta'],N=data['N'],bindings=data['bindings']))
        d=self.leading['chart'];self.coordinate=self.ref(d['coordinate']);self.radius=self.ref(d['radius'])
        self.mu=self.ref(data['mu']);self.delta=self.c.delta;self.Pstar=self.c.Pstar
        a=outer.Algebra(g,5);b=outer.Algebra(g,3);v=outer.Algebra(g,2)
        self.profiles=[{key:self.jet(d['modified_fields'][key]) for key in ('E','V')},
            {key:self.jet(d['radial'][key+'_y']) for key in ('E','V')}]
        bp=g.add(g.constant('1/2'),self.mu)
        for j in range(2,5):self.profiles.append(dict(E=a.scale(self.profiles[0]['E'],g.mul(g.constant((-1)**j),high.power(g,bp,j))),V=a.fixed(0)))
        self.histories=[{key:self.jet(d['modified_fields'][key]) for key in compact_source.RATES}]
        E,V=([row[key] for row in self.profiles] for key in ('E','V'))
        for n in range(4):
            EE=compact_high.product_row(a,E,E,n)
            density=dict(m=V[n],h=E[n],k=compact_high.product_row(a,E,V,n),
                e=a.sub(compact_high.product_row(a,V,V,n),a.scale(EE,'1/2')),p=a.scale(EE,'1/2'))
            row={key:a.sub(density[key],a.scale(self.histories[-1][key],rate)) for key,rate in compact_source.RATES.items()}
            if n==0:row['m']=self.jet(d['radial']['m_y'])
            self.histories.append(row)
        self.physical={}
        for j,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy')):
            if j<2:
                for key in ('U','m','Q'):self.physical[key+suffix]=self.jet(d['physical'][key+suffix])
            else:
                self.physical['U'+suffix]=a.scale(V[j],self.Pstar)
                self.physical['m'+suffix]=a.scale(self.histories[j]['m'],self.Pstar)
                self.physical['Q'+suffix]=radial.leading_Q(g,self.physical['U'+suffix],self.physical['m'+suffix],self.delta)
        self.density=radial.pressure_density(g,self.physical,self.delta,self.radius)
        self.higher_density=high.higher_pressure_density(g,self.physical,self.delta,self.radius,self.density)
        self.leading_left=self.at_leading(self.ref(d['domain'][1]))
        self.leading_right=self.c.at_compact(g.one)
        self.compact_inlet=copy.deepcopy(self.c.c['inlet'])
        intermediate=self.own.inputs[outer.intermediate.NAME]
        self.intermediate=copy.deepcopy(intermediate['charts']['O3_power'])
        self.intermediate_endpoint=copy.deepcopy(intermediate['endpoints']['O3_power']['right'])
        f=self.own.jets(self.intermediate['fields']);cumulative=self.own.jets(self.intermediate['cumulative'])
        fields={key:[f[key]]+[(v if key in ('Q1','V1') else b).fixed(0) for _ in range(4)] for key in previous.STATE}
        for key in ('mean1','K1','Q1'):
            alg=v if key=='Q1' else b;fields[key]=[alg.scale(f[key],(-1)**j) for j in range(5)]
        pressure=[f['P1'],self.density['P1_y'],self.higher_density['P1_yy'],self.higher_density['P1_yyy']]
        fields['P1']=pressure
        moment_rows={key:[cumulative[key]]+[b.fixed(0) for _ in range(3)] for key in previous.MOMENTS}
        moment_rows['m11']=fields['M1'][:4];moment_rows['m13']=pressure
        moment_rows['m15']=[cumulative['m15'],b.scale(pressure[1],g.neg(self.radius)),
            b.scale(b.add(pressure[1],pressure[2]),g.neg(self.radius)),
            b.scale(b.add(pressure[1],b.scale(pressure[2],2),pressure[3]),g.neg(self.radius))]
        self.left=dict(coordinate=self.coordinate,domain=[self.ref(i) for i in d['domain']],radius=self.radius,
            physical_y_Jacobian=g.one,fields=fields,cumulative=moment_rows)
        self.join=dict(left=self.at_left(self.left['domain'][1]),right=self.combined.at('reserved',g.zero))
        self.contract=copy.deepcopy(previous.CONTRACT)
        self.definition=dict(actual_O3_selected_finite_N_source_and_original_mu_eta_guard_required=True,
            flat_local_loop_does_not_zero_any_incoming_histories=True,
            actual_compact_x1_empty_integrals_and_shared_boundary_data_required=True,
            total_own_pressure_FTC_and_once_R_weighted_m15_sign_retained=True,
            finite_own_y4_C3_QV_C2_pressure_cumulative_y3_C3_contract=True,
            other_inner_joins_numeric_global_n1_recursion_corrected_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self._owners=(self.combined,self.c,self.own,g);self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def ref(self,i):return source.FunctionRef(self.g,i)
    def jet(self,ids):return outer.FunctionJet(tuple(self.ref(i) for i in ids))
    def at_leading(self,x):
        at=lambda packet:{key:outer.Algebra(self.g,len(outer.rows(jet))-1).at(jet,self.coordinate,x) for key,jet in packet.items()}
        return dict(profiles=[at(row) for row in self.profiles],histories=[at(row) for row in self.histories],
            physical=at(self.physical),density=at(self.density),higher_density=at(self.higher_density))
    def at_left(self,x):
        at=lambda packet:{key:[outer.Algebra(self.g,len(outer.rows(jet))-1).at(jet,self.coordinate,x) for jet in jets] for key,jets in packet.items()}
        return dict(fields=at(self.left['fields']),cumulative=at(self.left['cumulative']),
            radius=outer.Algebra(self.g,3).at(outer.Algebra(self.g,3).fixed(self.radius),self.coordinate,x).value)
    def _snapshot(self):
        return outer.encode(dict(leading=self.leading,coordinate=self.coordinate,radius=self.radius,
            profiles=self.profiles,histories=self.histories,physical=self.physical,density=self.density,higher_density=self.higher_density,
            leading_left=self.leading_left,leading_right=self.leading_right,compact_inlet=self.compact_inlet,
            intermediate=self.intermediate,intermediate_endpoint=self.intermediate_endpoint,left=self.left,join=self.join,
            contract=self.contract,definition=self.definition))
    def assert_graph(self):
        self.combined.assert_graph()
        if (self._owners!=(self.combined,self.c,self.own,self.g) or self.hashes!=self._hashes
            or self.snapshot!=self._snapshot() or self.g.nodes[:len(self.prefix)]!=self.prefix
            or self.numerical_oracle_available):raise ValueError('Actual own Rc source changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed Rc source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        if (not rec['all_passed'] or not all(rec[k] for k in GATES) or any(rec[k] for k in OPEN)
            or pins.pop(NAME)!=sha(NAME) or pins.pop(Path(__file__).stem+'_check.py')!=sha(Path(__file__).stem+'_check.py')
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['Rc_mixed_source_join_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):raise ValueError('Checked Rc source required')
        self.acceptance_loaded=True
    def at_left_fraction(self,fraction):
        self.assert_graph()
        if not self.acceptance_loaded:raise ValueError('Checked Rc source required')
        if not isinstance(fraction,(int,str,Fraction)) or isinstance(fraction,bool):raise TypeError('Exact rational fraction required')
        q=Fraction(fraction)
        if not 0<=q<=1:raise ValueError('Chart fraction must lie in [0,1]')
        lo,hi=self.left['domain'];x=self.g.add(lo,self.g.mul(self.g.constant(q),self.g.sub(hi,lo)))
        return dict(chart='O3_power',exact_native_fraction=str(q),**self.at_left(x),global_n1=False,numerical_oracle_available=False)


def run(combined=None):
    began=time.monotonic();owner=CurrentOriginalN1RcMixedSourceJoin(combined,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        Rc_mixed_source_join_bindings=owner.snapshot,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_O3_HIGHER_AND_OWN_RC_SOURCE_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw
