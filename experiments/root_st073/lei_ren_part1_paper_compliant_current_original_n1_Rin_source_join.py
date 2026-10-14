"""Actual cutoff/lower Rin join with retained pressure and moment histories.

State rows extend through physical y2, pressure/moments through y1.
Higher radial pressure, the full mixed contract and global n1 remain open.
"""
import copy
import gzip
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_n1_Rc_mixed_source_join as previous
import lei_ren_part1_paper_compliant_current_original_n1_inner_cutoff_functions as cutoff_source
import lei_ren_part1_paper_compliant_current_original_n1_analytic_core_lower_pressure_C3_functions as lower_high

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer=previous.outer
inner_source=cutoff_source.radial.inner
lower_source=inner_source.collar
NAME=PREFIX+'current_original_n1_Rin_source_join.json.gz'
RECEIPT=PREFIX+'current_original_n1_Rin_source_join_check.json'
GATES=('current_original_actual_Rin_inner_lower_pressure_source_identity_installed',
    'current_original_n1_Rin_own_y2_pressure_moment_y1_source_join_installed')
OPEN=tuple(dict.fromkeys(previous.OPEN+(
    'current_original_n1_Rin_full_y4_pressure_y3_mixed_source_join_installed',)))
CONTRACT=[dict(group=row['group'],field=row['field'],physical_y_order=(1 if row['group']=='cumulative' or row['field']=='P1' else 2),
    ordinary_Z_order=row['ordinary_Z_order']) for row in previous.previous.CONTRACT]


class CurrentOriginalN1RinSourceJoin:
    def __init__(self,rc=None,require_checked=True):
        self.rc=rc if rc is not None else previous.CurrentOriginalN1RcMixedSourceJoin()
        if not self.rc.acceptance_loaded:raise ValueError('Accepted actual Rc source required')
        self.rc.assert_graph();self.g=g=self.rc.g
        self.family=copy.deepcopy(self.rc.family);self.hashes=dict(self.rc.hashes)
        self.inputs={}
        specifications=((cutoff_source,'inner_cutoff_bindings'),(inner_source,'inner_Volterra_bindings'),
            (lower_source,'lower_collar_bindings'),(lower_high,'core_lower_C5_C3_bindings'))
        for module,key in specifications:
            raw=parent.linear.read(module.NAME)
            if raw['source_family']!=self.family or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]:
                raise ValueError('Same actual inner/lower source prefix required')
            self.inputs[key]=copy.deepcopy(raw[key])
            receipt=parent.linear.read(module.RECEIPT)
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed Rin ancestor '+name)
                self.hashes[name]=digest
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self.hashes[name]=sha(name)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        self.cut=self.inputs['inner_cutoff_bindings'];self.inner=self.inputs['inner_Volterra_bindings']
        self.lower=self.inputs['lower_collar_bindings'];self.leading=self.inputs['core_lower_C5_C3_bindings']
        self.parameters={key:self.ref(node) for key,node in self.lower['parameters'].items()}
        self.coordinate=self.ref(self.inner['coordinates']['first']);self.endpoint=self.parameters['sc_quarter']
        self.rho=self.ref(self.lower['fields']['rho']);self.radius=self.ref(self.lower['fields']['radius'])
        eps=self.parameters['epsilon'];g2=g.mul(self.rho,self.rho)
        a=outer.Algebra(g,3);b=outer.Algebra(g,2)
        cut=self.cut['fields']['first'];native=cut['rows']
        self.native_fields=copy.deepcopy(cut)
        fields={key:[] for key in previous.previous.STATE}
        for key in ('F1','Uz1','mean1','K1'):
            v,r,rr=map(self.jet,native[key])
            fields[key]=[v,a.scale(r,self.rho),a.add(a.scale(r,self.rho),a.scale(rr,g2))]
        mean=list(map(self.jet,native['mean1']))
        fields['M1']=[self.jet(cut['radial_primitive_M1']),
            a.scale(a.add(mean[0],a.scale(mean[1],self.rho)),self.radius),
            a.scale(a.add(mean[0],a.scale(mean[1],g.mul(g.constant(3),self.rho)),a.scale(mean[2],g2)),self.radius)]
        F=list(map(self.jet,native['F1']));sqrt2R=g.unary('sqrt',g.mul(g.constant(2),self.radius))
        fields['Utheta1']=[a.scale(F[0],sqrt2R),
            a.scale(a.add(a.scale(F[0],'1/2'),a.scale(F[1],self.rho)),sqrt2R),
            a.scale(a.add(a.scale(F[0],'1/4'),a.scale(F[1],g.mul(g.constant(2),self.rho)),a.scale(F[2],g2)),sqrt2R)]
        for key in ('Q1','V1'):
            v,r,rr=map(self.jet,cut[key]);fields[key]=[v,b.scale(r,self.rho),b.add(b.scale(r,self.rho),b.scale(rr,g2))]
        fields['P1']=[self.jet(cut['pressure_P1']),a.scale(self.jet(cut['pressure_rho']),self.rho)]
        moments={key:[self.jet(value),a.scale(self.jet(self.cut['integrands']['first'][key]),self.rho)]
            for key,value in self.cut['cumulative']['first'].items()}
        self.left=dict(fields=fields,cumulative=moments)
        self.intermediate=copy.deepcopy(self.rc.own.inputs[outer.intermediate.NAME]['charts']['lower'])
        f={key:self.jet(value) for key,value in self.intermediate['fields'].items()}
        right_fields={key:[f[key]]+[(b if key in ('Q1','V1') else a).fixed(0) for _ in range(2)]
            for key in previous.previous.STATE}
        for key in ('mean1','K1','Q1'):
            alg=b if key=='Q1' else a
            right_fields[key]=[alg.scale(f[key],(-1)**j) for j in range(3)]
        N=self.jet(self.leading['N1p'])
        right_fields['P1']=[f['P1'],a.scale(N,self.radius)]
        right_moments={key:[self.jet(value),a.fixed(0)] for key,value in self.intermediate['cumulative'].items()}
        right_moments['m13'][1]=right_fields['P1'][1]
        right_moments['m15'][1]=a.scale(N,g.neg(g.mul(self.radius,self.radius)))
        self.right=dict(fields=right_fields,cumulative=right_moments)
        self.join={key:self.at(packet,self.endpoint) for key,packet in (('left',self.left),('right',self.right))}
        self.contract=copy.deepcopy(CONTRACT)
        self.definition=dict(actual_inner_leading_first_is_lower_source=True,
            Lambda_epsilon_one_and_physical_y_equals_rho_d_rho=True,
            exact_original_flat_cutoff_argument_at_Rin_is_two=True,
            own_velocity_zero_does_not_zero_mean_pressure_or_moment_history=True,
            positive_order_mean_weight='1+delta',leading_pressure_weight='1-delta',
            actual_native_rho0_rho1_rho2_state_and_own_recovery_retained=True,
            actual_signed_cutoff_FTC_and_empty_lower_integrals_required=True,
            state_y0_y2_C3_QV_C2_pressure_moments_y0_y1_C3=True,
            higher_radial_pressure_full_mixed_contract_global_n1_numeric_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def ref(self,node):return source.FunctionRef(self.g,node)
    def jet(self,ids):return outer.FunctionJet(tuple(map(self.ref,ids)))
    def at(self,packet,value):
        return {group:{key:[outer.Algebra(self.g,len(outer.rows(jet))-1).at(jet,self.coordinate,value) for jet in jets]
            for key,jets in data.items()} for group,data in packet.items()}
    def _snapshot(self):return outer.encode(dict(inputs=self.inputs,coordinate=self.coordinate,endpoint=self.endpoint,
        rho=self.rho,radius=self.radius,native_fields=self.native_fields,intermediate=self.intermediate,
        left=self.left,right=self.right,join=self.join,contract=self.contract,definition=self.definition))
    def assert_graph(self):
        self.rc.assert_graph()
        if (self.snapshot!=self._snapshot() or self.g.nodes[:len(self.prefix)]!=self.prefix
            or self.hashes!=self._hashes or self.numerical_oracle_available):raise ValueError('Actual Rin source changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed Rin source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        if (not rec['all_passed'] or not all(rec[key] for key in GATES) or any(rec[key] for key in OPEN)
            or pins.pop(NAME)!=sha(NAME) or pins.pop(Path(__file__).stem+'_check.py')!=sha(Path(__file__).stem+'_check.py')
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['Rin_source_join_bindings']!=self.snapshot or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual Rin source required')
        self.acceptance_loaded=True


def run(rc=None):
    began=time.monotonic();owner=CurrentOriginalN1RinSourceJoin(rc,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        Rin_source_join_bindings=owner.snapshot,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_RIN_SOURCE_JOIN_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw
