"""Actual frozen-macro high sources and the coordinate-changing own seam."""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_micro_active_radial4_pressure_C3_functions as previous

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer,high,rc_source=previous.outer,previous.high,previous.rc_source
NAME=PREFIX+'current_original_micro_macro_mixed_source_join.json.gz'
RECEIPT=PREFIX+'current_original_micro_macro_mixed_source_join_check.json'
GATES=('current_original_actual_frozen_macro_U_mean_y4_C5_Q_y4_C4_pressure_y3_C3_sources_installed',
    'current_original_n1_second_micro_frozen_macro_full_mixed_source_join_installed')
OPEN=tuple(key for key in previous.OPEN if key not in GATES)
CHARTS=('second_micro','frozen_macro')


class CurrentOriginalMicroMacroMixedSourceJoin:
    def __init__(self,base=None,require_checked=True):
        self.base=base if base is not None else previous.CurrentOriginalMicroActiveRadial4PressureC3Functions()
        if not self.base.acceptance_loaded:raise ValueError('Accepted actual active micro source required')
        self.base.assert_graph();self.g=g=self.base.g;self.family=copy.deepcopy(self.base.family)
        self.hashes=dict(self.base.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        self.rminus=self.base.base;self.lower=self.base.lower;self.parameters=dict(self.base.parameters)
        self.bridge=copy.deepcopy(self.rminus.bridge);self.loop_input=copy.deepcopy(self.base.loop_input)
        self.modified_input=copy.deepcopy(self.base.modified_input)
        self.loop={name:self.loop_input['charts'][name] for name in CHARTS}
        self.modified={name:self.modified_input['charts'][name] for name in CHARTS}
        self.background=self.bridge['macro'];self.coordinate=self.ref(self.background['phase'])
        self.radius=self.ref(self.background['radius']);self.jacobian=self.ref(self.background['physical_measure'])
        self.phase=self.ref(self.loop['frozen_macro']['actual_phase']);self.N=self.base.N
        self.background_y={key:self.jet(ids) for key,ids in self.modified['frozen_macro']['radial_background'].items()}
        D=previous.PeriodicTotalNativeDerivative(g,self.coordinate,self.phase,self.coordinate,
            self.jacobian,self.N,self.base.guard_core)
        D.shared['family']=self.family;self.factory=D;self.requests=[]
        for key,ids in self.background['fields'].items():
            for node,dy in zip(ids,outer.rows(self.background_y[key])):D.seed(self.ref(node),g.mul(self.jacobian,dy))
        a=outer.Algebra(g,5);b=outer.Algebra(g,3);old=self.modified['frozen_macro']
        invJ=g.quotient(g.one,self.jacobian,'same positive original frozen-macro physical Jacobian')
        physical={key:[self.jet(old['physical'][key]),self.jet(old['physical'][key+'_y'])] for key in ('U','m','Q')}
        for n in range(2,5):
            refs=[]
            for ref in outer.rows(physical['U'][-1]):
                native=D.derivative(ref);refs.append(g.mul(invJ,native))
                self.requests.append(dict(input=ref.node,native_output=native.node,physical_output=refs[-1].node,physical_y_order=n))
            physical['U'].append(outer.FunctionJet(tuple(refs)))
            physical['m'].append(a.sub(physical['U'][n-1],physical['m'][n-1]))
            physical['Q'].append(rc_source.radial.leading_Q(g,physical['U'][n],physical['m'][n],self.parameters['delta']))
        f={key+suffix:physical[key][n] for key in ('U','m','Q') for n,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy'))}
        density=rc_source.radial.pressure_density(g,f,self.parameters['delta'],self.radius)
        higher=high.higher_pressure_density(g,f,self.parameters['delta'],self.radius,density)
        Ny=[density['N1p'],higher['N1p_y'],higher['N1p_yy']];n0,n1,n2=Ny
        raw_P=[density['P1_y'],higher['P1_yy'],higher['P1_yyy']]
        raw_m15=[b.scale(raw_P[0],g.neg(self.radius)),higher['m15_yy'],higher['m15_yyy']]
        negative_R2=g.neg(g.mul(self.radius,self.radius))
        P=[b.scale(n0,self.radius),b.scale(b.add(n0,n1),self.radius),b.scale(b.add(n0,b.scale(n1,2),n2),self.radius)]
        m15=[b.scale(n0,negative_R2),b.scale(b.add(b.scale(n0,2),n1),negative_R2),b.scale(b.add(b.scale(n0,4),b.scale(n1,4),n2),negative_R2)]
        before=previous.previous.common_source.previous
        intermediate=copy.deepcopy(self.lower.rin.rc.own.inputs[before.NAME]['charts']['frozen_macro'])
        values={key:self.jet(ids) for key,ids in intermediate['fields'].items()};fields={}
        for key in rc_source.previous.STATE:
            alg=outer.Algebra(g,2 if key in ('Q1','V1') else 3)
            fields[key]=[values[key]]+[alg.scale(values[key],(-1)**n) if key in ('mean1','K1','Q1') else alg.fixed(0) for n in range(1,5)]
        fields['P1']=[values['P1'],*P]
        cumulative={key:[self.jet(ids)]+[b.fixed(0) for _ in range(3)] for key,ids in intermediate['cumulative'].items()}
        cumulative['m13']=list(fields['P1']);cumulative['m15'][1:]=m15
        self.charts={'second_micro':self._second_micro_packet(),
            'frozen_macro':dict(domain=[g.zero,g.constant('1/8')],physical=physical,density=density,higher_density=higher,
                N_y=Ny,pressure_y=P,m15_y=m15,raw_pressure_y=raw_P,raw_m15_y=raw_m15,
                intermediate=intermediate,own=dict(fields=fields,cumulative=cumulative))}
        self.contract=copy.deepcopy(self.base.contract)
        self.join=dict(left=self.at('second_micro',self.charts['second_micro']['own'],g.constant(2)),
            right=self.at('frozen_macro',self.charts['frozen_macro']['own'],g.zero))
        self.definition=dict(actual_macro_has_fresh_native_factory_and_N_times_length_phase_lift=True,
            actual_second_micro_endpoint_packet_rebuilt_from_original_full_chart_sources=True,
            micro_hb_and_macro_length_native_Jacobians_remain_distinct=True,
            actual_comparison_frozen_profile_and_nonzero_six_incoming_histories_retained=True,
            actual_modified_U_y_differentiated_and_mean_ODE_Q_recovery_retained=True,
            actual_pressure_operator_signed_primitives_and_full_own_moment_inputs_retained=True,
            finite_source_seam_not_global_n1_or_numerical_mixed_norms=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def _second_micro_packet(self):
        """Derive endpoint sources anew; do not relabel the near-one packet."""
        g=self.g;old=self.modified['second_micro'];bg=self.bridge['micro']
        coordinate=self.ref(bg['phase']);radius=self.ref(bg['radius']);J=self.ref(bg['physical_measure'])
        phase=self.ref(self.loop['second_micro']['actual_phase'])
        self.second_background_y={key:self.jet(ids) for key,ids in old['radial_background'].items()}
        D=previous.PeriodicTotalNativeDerivative(g,coordinate,phase,coordinate,J,self.N,self.base.guard_core)
        D.shared['family']=self.family;self.second_factory=D;self.second_requests=[]
        for key,ids in bg['fields'].items():
            for node,dy in zip(ids,outer.rows(self.second_background_y[key])):D.seed(self.ref(node),g.mul(J,dy))
        a=outer.Algebra(g,5);b=outer.Algebra(g,3)
        invJ=g.quotient(g.one,J,'same positive original micro physical Jacobian')
        physical={key:[self.jet(old['physical'][key]),self.jet(old['physical'][key+'_y'])] for key in ('U','m','Q')}
        for n in range(2,5):
            refs=[]
            for ref in outer.rows(physical['U'][-1]):
                native=D.derivative(ref);refs.append(g.mul(invJ,native))
                self.second_requests.append(dict(input=ref.node,native_output=native.node,
                    physical_output=refs[-1].node,physical_y_order=n))
            physical['U'].append(outer.FunctionJet(tuple(refs)))
            physical['m'].append(a.sub(physical['U'][n-1],physical['m'][n-1]))
            physical['Q'].append(rc_source.radial.leading_Q(g,physical['U'][n],physical['m'][n],self.parameters['delta']))
        f={key+suffix:physical[key][n] for key in ('U','m','Q')
            for n,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy'))}
        density=rc_source.radial.pressure_density(g,f,self.parameters['delta'],radius)
        higher=high.higher_pressure_density(g,f,self.parameters['delta'],radius,density)
        Ny=[density['N1p'],higher['N1p_y'],higher['N1p_yy']];n0,n1,n2=Ny
        P=[b.scale(n0,radius),b.scale(b.add(n0,n1),radius),b.scale(b.add(n0,b.scale(n1,2),n2),radius)]
        negative_R2=g.neg(g.mul(radius,radius))
        m15=[b.scale(n0,negative_R2),b.scale(b.add(b.scale(n0,2),n1),negative_R2),
            b.scale(b.add(b.scale(n0,4),b.scale(n1,4),n2),negative_R2)]
        before=previous.previous.common_source.previous
        intermediate=copy.deepcopy(self.lower.rin.rc.own.inputs[before.NAME]['charts']['second_micro'])
        values={key:self.jet(ids) for key,ids in intermediate['fields'].items()};fields={}
        for key in rc_source.previous.STATE:
            alg=outer.Algebra(g,2 if key in ('Q1','V1') else 3)
            fields[key]=[values[key]]+[alg.scale(values[key],(-1)**n)
                if key in ('mean1','K1','Q1') else alg.fixed(0) for n in range(1,5)]
        fields['P1']=[values['P1'],*P]
        cumulative={key:[self.jet(ids)]+[b.fixed(0) for _ in range(3)] for key,ids in intermediate['cumulative'].items()}
        cumulative['m13']=list(fields['P1']);cumulative['m15'][1:]=m15
        return dict(domain=[g.constant('15/8'),g.constant(2)],physical=physical,density=density,higher_density=higher,
            N_y=Ny,pressure_y=P,m15_y=m15,raw_pressure_y=[density['P1_y'],higher['P1_yy'],higher['P1_yyy']],
            raw_m15_y=[b.scale(density['P1_y'],g.neg(radius)),higher['m15_yy'],higher['m15_yyy']],
            intermediate=intermediate,own=dict(fields=fields,cumulative=cumulative))
    def ref(self,node):return source.FunctionRef(self.g,node)
    def jet(self,ids):return outer.FunctionJet(tuple(map(self.ref,ids)))
    def at(self,name,packet,value):
        coordinate=self.ref(self.loop[name]['coordinate'])
        return {group:{key:[outer.Algebra(self.g,len(outer.rows(jet))-1).at(jet,coordinate,value)
            for jet in jets] for key,jets in fields.items()} for group,fields in packet.items()}
    def _snapshot(self):
        return outer.encode(dict(coordinate=self.coordinate,radius=self.radius,jacobian=self.jacobian,phase=self.phase,N=self.N,
            background=self.background,background_y=self.background_y,loop=self.loop,modified=self.modified,
            charts=self.charts,requests=self.requests,second_requests=self.second_requests,
            second_background_y=self.second_background_y,contract=self.contract,join=self.join,definition=self.definition))
    def assert_graph(self):
        self.base.assert_graph()
        if (self.snapshot!=self._snapshot() or self.g.nodes[:len(self.prefix)]!=self.prefix
            or self.hashes!=self._hashes or self.numerical_oracle_available):raise ValueError('Actual micro/macro source changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed micro/macro source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        if (not rec['all_passed'] or len(rec.get('checks',{}))!=7
            or not all(rec[key] for key in GATES) or any(rec[key] for key in OPEN)
            or pins.pop(NAME)!=sha(NAME) or pins.pop(Path(__file__).stem+'_check.py')!=sha(Path(__file__).stem+'_check.py')
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['micro_macro_mixed_source_join_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):raise ValueError('Checked actual micro/macro source required')
        self.acceptance_loaded=True


def run(base=None):
    began=time.monotonic();owner=CurrentOriginalMicroMacroMixedSourceJoin(base,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        micro_macro_mixed_source_join_bindings=owner.snapshot,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_FROZEN_MACRO_HIGH_SOURCE_DEFINED',len(owner.g.nodes),flush=True);return owner,raw
