"""Actual own rminus join and genuine high first-micro flat-collar source.

The first-micro derivatives use its actual ratio/U/mean primitives and
the original common-core IVP; lower output derivatives are not copied.
"""
import copy
import gzip
import json
import math
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_lower_higher_radial_Rin_join as previous
import lei_ren_part1_paper_compliant_current_original_lower_active_common_source_join_functions as common_source
import lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions as modified_source

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer,high,rc_source=previous.outer,previous.high,previous.rc_source
bridge_source=common_source.bridge_source
loop_source=modified_source.previous
NAME=PREFIX+'current_original_n1_rminus_mixed_source_join.json.gz'
RECEIPT=PREFIX+'current_original_n1_rminus_mixed_source_join_check.json'
GATES=('current_original_actual_first_micro_flat_collar_high_radial_C5_Q4_source_installed',
    'current_original_n1_rminus_own_y4_pressure_cumulative_y3_mixed_source_join_installed',
    'current_original_n1_Rin_rminus_own_mixed_source_joins_installed')
OPEN=tuple(key for key in previous.OPEN if key not in GATES)
CONTRACT=copy.deepcopy(previous.CONTRACT)
STIRLING=((1,),(0,1),(0,1,1),(0,1,3,1),(0,1,7,6,1))


def y_row(a,g,rho,native,n):
    if n==0:return native[0]
    return a.add(*(a.scale(native[k],g.mul(g.constant(STIRLING[n][k]),previous.high.power(g,rho,k))) for k in range(1,n+1)))


class CurrentOriginalN1RminusMixedSourceJoin:
    def __init__(self,lower=None,require_checked=True):
        self.lower=lower if lower is not None else previous.CurrentOriginalLowerHigherRadialRinJoin()
        if not self.lower.acceptance_loaded:raise ValueError('Accepted actual high lower/Rin source required')
        self.lower.assert_graph();self.g=g=self.lower.g;self.family=copy.deepcopy(self.lower.family);self.hashes=dict(self.lower.hashes)
        self.inputs={}
        for module,key in ((common_source,'lower_active_common_source_bindings'),(bridge_source,'prescribed_bridge_C6_bindings'),
            (modified_source,'modified_Q4_pressure_C3_bindings'),(loop_source,'prescribed_loop_C5_bindings')):
            raw=parent.linear.read(module.NAME);rec=parent.linear.read(module.RECEIPT)
            if raw['source_family']!=self.family or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]:
                raise ValueError('Actual same-family first-micro/common source prefix required')
            for name,digest in rec['input_hashes'].items():
                if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed rminus ancestor '+name)
                self.hashes[name]=digest
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self.hashes[name]=sha(name)
            self.inputs[key]=copy.deepcopy(raw[key])
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        self.common=self.inputs['lower_active_common_source_bindings'];self.bridge=self.inputs['prescribed_bridge_C6_bindings']
        self.loop=self.inputs['prescribed_loop_C5_bindings']['charts']['first_micro']
        self.modified=self.inputs['modified_Q4_pressure_C3_bindings']['charts']['first_micro']
        self.micro=self.bridge['micro'];self.parameters=dict(self.lower.parameters)
        self.coordinate=self.ref(self.loop['coordinate']);self.radius=self.ref(self.loop['radius'])
        self.jacobian=self.ref(self.loop['physical_measure']);self.active_sc=self.ref(self.common['sc_binding']['active_sc_node'])
        self.left_endpoint=self.parameters['sc_half'];self.right_endpoint=self.ref(self.loop['domain'][0])
        self.collar_domain=[self.right_endpoint,g.mul(g.constant('3/4'),self.active_sc)]
        self.rho=g.mul(g.constant(4),g.unary('exp',g.mul(self.parameters['hb'],self.coordinate)))
        a=outer.Algebra(g,5);b=outer.Algebra(g,3);v=outer.Algebra(g,2)
        eps,delta=self.parameters['epsilon'],self.parameters['delta']
        self.collar_core={key:[a.at(jet,self.lower.core_rho,self.rho) for jet in jets] for key,jets in self.lower.core_radial.items()}
        # These scalar chi derivatives come from the identical original
        # sigma/hb source; actual field primitives below remain distinct.
        self.chi=[a.fixed(self.ref(self.micro['chi']))]+[a.at(jet,self.lower.coordinate,self.coordinate) for jet in self.lower.chi[1:]]
        cp=self.collar_core['Phi'];self.A=[a.div(cp[1],cp[0],'same positive actual comparison Phi')]
        self.A.append(a.div(a.sub(cp[2],a.mul(self.A[0],cp[1])),cp[0],'same positive actual comparison Phi'))
        self.A.append(a.div(a.sub(cp[3],a.add(a.mul(self.A[0],cp[2]),a.scale(a.mul(self.A[1],cp[1]),2))),cp[0],'same positive actual comparison Phi'))
        cm1=[a.sub(self.chi[0],a.fixed(1)),*self.chi[1:]]
        self.g_rows=[previous.product_row(a,cm1,self.A,n) for n in range(3)]
        self.ratio=[a.trim(self.jet(self.micro['actual']['ratio']))]
        for n in range(1,4):self.ratio.append(previous.product_row(a,self.ratio,self.g_rows,n-1))
        self.B=[previous.product_row(a,self.ratio,self.chi,n) for n in range(4)]
        physical=self.modified['physical'];axis=a.trim(self.jet(self.bridge['core7']['U0']))
        self.native=dict(Phi=[a.trim(self.jet(self.micro['actual']['phi']))],
            Psi=[a.scale(a.sub(self.jet(physical['U']),axis),g.quotient(g.one,eps,'same positive original epsilon'))],
            mean_Psi=[a.scale(a.sub(self.jet(physical['m']),axis),g.quotient(g.one,eps,'same positive original epsilon'))])
        for key in ('Phi','Psi'):
            for n in range(1,5):self.native[key].append(a.add(*(a.scale(a.mul(self.B[k],self.collar_core[key][n-k]),math.comb(n-1,k)) for k in range(n))))
        for n in range(1,5):self.native['mean_Psi'].append(a.scale(a.sub(self.native['Psi'][n-1],a.scale(self.native['mean_Psi'][-1],n)),g.quotient(g.one,self.rho,'same positive original rho')))
        self.native['U']=[self.jet(physical['U'])]+[a.scale(jet,eps) for jet in self.native['Psi'][1:]]
        self.native['m']=[self.jet(physical['m'])]+[a.scale(jet,eps) for jet in self.native['mean_Psi'][1:]]
        self.physical={};self.reconstructed_y={}
        for key in ('U','m'):
            candidate=[y_row(a,g,self.rho,self.native[key],n) for n in range(5)]
            self.reconstructed_y[key]=candidate
            self.physical[key]=[self.jet(physical[key]),self.jet(physical[key+'_y']),*candidate[2:]]
        self.physical['Q']=[self.jet(physical['Q']),self.jet(physical['Q_y'])]+[
            rc_source.radial.leading_Q(g,self.physical['U'][n],self.physical['m'][n],delta) for n in range(2,5)]
        self.recovered_Q=[rc_source.radial.leading_Q(g,U,m,delta) for U,m in zip(self.physical['U'],self.physical['m'])]
        f={key+suffix:self.physical[key][n] for key in ('U','m','Q') for n,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy'))}
        self.density=rc_source.radial.pressure_density(g,f,delta,self.radius)
        self.higher_density=high.higher_pressure_density(g,f,delta,self.radius,self.density)
        self.raw_pressure_y=[self.density['P1_y'],self.higher_density['P1_yy'],self.higher_density['P1_yyy']]
        self.raw_m15_y=[b.scale(self.raw_pressure_y[0],g.neg(self.radius)),self.higher_density['m15_yy'],self.higher_density['m15_yyy']]
        self.N_y=[self.density['N1p'],self.higher_density['N1p_y'],self.higher_density['N1p_yy']]
        N,Ny,Nyy=self.N_y;negative_R2=g.neg(g.mul(self.radius,self.radius))
        self.pressure_y=[b.scale(N,self.radius),b.scale(b.add(N,Ny),self.radius),b.scale(b.add(N,b.scale(Ny,2),Nyy),self.radius)]
        self.m15_y=[b.scale(N,negative_R2),b.scale(b.add(b.scale(N,2),Ny),negative_R2),b.scale(b.add(b.scale(N,4),b.scale(Ny,4),Nyy),negative_R2)]
        self.intermediate=copy.deepcopy(self.lower.rin.rc.own.inputs[common_source.previous.NAME]['charts']['first_micro'])
        values={key:self.jet(ids) for key,ids in self.intermediate['fields'].items()}
        fields={}
        for key in rc_source.previous.STATE:
            alg=v if key in ('Q1','V1') else b
            fields[key]=[values[key]]+[(alg.scale(values[key],(-1)**n) if key in ('mean1','K1','Q1') else alg.fixed(0)) for n in range(1,5)]
        fields['P1']=[values['P1'],*self.pressure_y]
        cumulative={key:[self.jet(ids)]+[b.fixed(0) for _ in range(3)] for key,ids in self.intermediate['cumulative'].items()}
        cumulative['m13']=list(fields['P1']);cumulative['m15'][1:]=self.m15_y
        self.right=dict(fields=fields,cumulative=cumulative)
        self.join=dict(left=self.lower.rin.at(self.lower.right,self.left_endpoint),right=self.at(self.right,self.right_endpoint))
        self.contract=copy.deepcopy(CONTRACT)
        self.definition=dict(actual_first_micro_primitives_ratio_Phi_U_mean_and_all_background_histories_retained=True,
            first_micro_finite_N_flat_source_not_positive_order_own_incoming_data=True,
            genuine_active_radial4_from_actual_common_core_IVP_not_copied_lower_output=True,
            active_pressure_y1_y2_y3_from_actual_reduced_pressure_primitives_and_FTC=True,
            physical_y_is_native_partial_divided_by_original_hb_once=True,
            distinct_exact_sc_nodes_resolved_only_by_accepted_original_selector_binding=True,
            nonzero_own_moment_pressure_histories_and_one_plus_delta_recovery_retained=True,
            active_high_pressure_scope_is_first_micro_original_flat_collar=True,
            full_finite_rminus_y4_pressure_cumulative_y3_C3_QV_C2_contract=True,
            remaining_intermediate_joins_global_n1_numeric_all_order_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def ref(self,i):return source.FunctionRef(self.g,i)
    def jet(self,ids):return outer.FunctionJet(tuple(map(self.ref,ids)))
    def at(self,packet,value):
        return {group:{key:[outer.Algebra(self.g,len(outer.rows(jet))-1).at(jet,self.coordinate,value) for jet in jets]
            for key,jets in data.items()} for group,data in packet.items()}
    def _snapshot(self):return outer.encode(dict(inputs=self.inputs,coordinate=self.coordinate,radius=self.radius,jacobian=self.jacobian,
        active_sc=self.active_sc,left_endpoint=self.left_endpoint,right_endpoint=self.right_endpoint,collar_domain=self.collar_domain,
        rho=self.rho,collar_core=self.collar_core,chi=self.chi,A=self.A,g_rows=self.g_rows,ratio=self.ratio,B=self.B,
        native=self.native,physical=self.physical,reconstructed_y=self.reconstructed_y,recovered_Q=self.recovered_Q,
        density=self.density,higher_density=self.higher_density,N_y=self.N_y,raw_pressure_y=self.raw_pressure_y,
        raw_m15_y=self.raw_m15_y,pressure_y=self.pressure_y,m15_y=self.m15_y,
        intermediate=self.intermediate,right=self.right,join=self.join,contract=self.contract,definition=self.definition))
    def assert_graph(self):
        self.lower.assert_graph()
        if (self.snapshot!=self._snapshot() or self.g.nodes[:len(self.prefix)]!=self.prefix or self.hashes!=self._hashes
            or self.numerical_oracle_available):raise ValueError('Actual own rminus source changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed rminus source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        if (not rec['all_passed'] or not all(rec[key] for key in GATES) or any(rec[key] for key in OPEN)
            or pins.pop(NAME)!=sha(NAME) or pins.pop(Path(__file__).stem+'_check.py')!=sha(Path(__file__).stem+'_check.py')
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['rminus_mixed_source_join_bindings']!=self.snapshot or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual own rminus source required')
        self.acceptance_loaded=True


def run(lower=None):
    began=time.monotonic();owner=CurrentOriginalN1RminusMixedSourceJoin(lower,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        rminus_mixed_source_join_bindings=owner.snapshot,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_FIRST_MICRO_HIGH_SOURCE_AND_OWN_RMINUS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw
