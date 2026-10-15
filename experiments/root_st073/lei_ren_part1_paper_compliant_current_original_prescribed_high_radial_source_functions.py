"""Actual high radial source functions on the ten prescribed loop charts.

Each packet is freshly derived from its full original chart, with its own
physical Jacobian and periodic lift. Source coverage does not admit global
mixed norms, numerical values, global n1, or temporal coefficient recursion.
"""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_micro_macro_mixed_source_join as previous
import lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions as downstream_source

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer,high,rc_source=previous.outer,previous.high,previous.rc_source
NAME=PREFIX+'current_original_prescribed_high_radial_source_functions.json.gz'
RECEIPT=PREFIX+'current_original_prescribed_high_radial_source_functions_check.json'
GATES=('current_original_ten_prescribed_full_chart_U_mean_y4_C5_Q_y4_C4_sources_installed',
    'current_original_ten_prescribed_full_chart_N_y2_pressure_y3_C3_sources_installed',
    'current_original_ten_prescribed_beforeRc_own_high_source_functions_installed',
    'current_original_n1_frozen_macro_first_switch_full_mixed_source_join_installed')
OPEN=tuple(key for key in previous.OPEN if key not in GATES)
CHARTS=('first_micro','second_micro','frozen_macro','first_switch','second_switch',
    'post_power','long_reshape','reference','restoration','postrestore')
BACKGROUND_NAMES=dict(post_power='post_switch',postrestore='post_restore')


class CurrentOriginalPrescribedHighRadialSourceFunctions:
    def __init__(self,base=None,require_checked=True):
        self.base=base if base is not None else previous.CurrentOriginalMicroMacroMixedSourceJoin()
        if not self.base.acceptance_loaded:raise ValueError('Accepted actual micro/macro seam required')
        self.base.assert_graph();self.g=g=self.base.g;self.family=copy.deepcopy(self.base.family)
        self.hashes=dict(self.base.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        raw=parent.linear.read(downstream_source.NAME);rec=parent.linear.read(downstream_source.RECEIPT)
        if not rec['all_passed'] or raw['source_family']!=self.family or not all(rec[key] for key in downstream_source.GATES):
            raise ValueError('Accepted original full downstream background required')
        if raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]:raise ValueError('Same downstream source graph required')
        for name,digest in rec['input_hashes'].items():
            if self.hashes.get(name)!=digest or sha(name)!=digest:raise ValueError('Changed downstream input '+name)
        self.downstream=copy.deepcopy(raw['prescribed_R100_Rm_C6_bindings'])
        self.active=self.base.base;self.rminus=self.base.rminus;self.lower=self.base.lower
        self.parameters=dict(self.base.parameters);self.N=self.base.N
        self.bridge=copy.deepcopy(self.base.bridge);self.loop_input=copy.deepcopy(self.base.loop_input)
        self.modified_input=copy.deepcopy(self.base.modified_input)
        self.loop={name:self.loop_input['charts'][name] for name in CHARTS}
        self.modified={name:self.modified_input['charts'][name] for name in CHARTS}
        self.backgrounds={name:copy.deepcopy(self.bridge['micro' if name in CHARTS[:2] else 'macro']
            if name in CHARTS[:3] else self.downstream['charts'][BACKGROUND_NAMES.get(name,name)]) for name in CHARTS}
        before=previous.previous.previous.common_source.previous
        self.intermediate=copy.deepcopy(self.lower.rin.rc.own.inputs[before.NAME])
        self.factories={};self.requests={};self.background_y={};self.charts={}
        for name in CHARTS:
            self.charts[name]=self.build_chart(name)
            print('ACTUAL_FULL_CHART_HIGH_SOURCE_DEFINED',name,len(g.nodes),flush=True)
        self.contract=copy.deepcopy(self.base.contract)
        self.macro_first_join=dict(left=self.at_own('frozen_macro',g.one),right=self.at_own('first_switch',g.zero))
        self.definition=dict(ten_actual_original_domains_not_relabelled_seam_neighborhoods=True,
            every_chart_fresh_total_native_factory_uses_its_actual_N_times_J=True,
            full_actual_U_y_differentiation_mean_ODE_and_one_minus_delta_Q=True,
            original_signed_pressure_and_nonzero_predecessor_histories_retained=True,
            own_one_plus_delta_and_actual_beforeRc_support_retained=True,
            source_coverage_not_global_interface_closure_or_numerical_bounds=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()

    def build_chart(self,name):
        g=self.g;old=self.modified[name];bg=self.backgrounds[name];loop=self.loop[name]
        coordinate=self.ref(loop['coordinate']);radius=self.ref(loop['radius']);J=self.ref(loop['physical_measure'])
        phase=self.ref(loop['actual_phase']);self.background_y[name]={key:self.jet(ids) for key,ids in old['radial_background'].items()}
        D=previous.previous.PeriodicTotalNativeDerivative(g,coordinate,phase,coordinate,J,self.N,self.active.guard_core)
        D.shared['family']=self.family;self.factories[name]=D;self.requests[name]=[]
        for key,ids in bg['fields'].items():
            for node,dy in zip(ids,outer.rows(self.background_y[name][key])):D.seed(self.ref(node),g.mul(J,dy))
        a=outer.Algebra(g,5);b=outer.Algebra(g,3)
        certificate=('same positive original micro physical Jacobian' if name in CHARTS[:2] else
            'same positive original frozen-macro physical Jacobian' if name=='frozen_macro' else
            'same positive original prescribed chart physical Jacobian')
        invJ=g.quotient(g.one,J,certificate)
        physical={key:[self.jet(old['physical'][key]),self.jet(old['physical'][key+'_y'])] for key in ('U','m','Q')}
        for n in range(2,5):
            refs=[]
            for ref in outer.rows(physical['U'][-1]):
                native=D.derivative(ref);refs.append(g.mul(invJ,native))
                self.requests[name].append(dict(input=ref.node,native_output=native.node,
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
        packet=copy.deepcopy(self.intermediate['charts'][name])
        values={key:self.jet(ids) for key,ids in packet['fields'].items()};fields={}
        for key in rc_source.previous.STATE:
            alg=outer.Algebra(g,2 if key in ('Q1','V1') else 3)
            fields[key]=[values[key]]+[alg.scale(values[key],(-1)**n)
                if key in ('mean1','K1','Q1') else alg.fixed(0) for n in range(1,5)]
        fields['P1']=[values['P1'],*P]
        cumulative={key:[self.jet(ids)]+[b.fixed(0) for _ in range(3)] for key,ids in packet['cumulative'].items()}
        cumulative['m13']=list(fields['P1']);cumulative['m15'][1:]=m15
        return dict(domain=[self.ref(node) for node in loop['domain']],coordinate=coordinate,radius=radius,
            jacobian=J,phase=phase,physical=physical,density=density,higher_density=higher,N_y=Ny,
            pressure_y=P,m15_y=m15,raw_pressure_y=[density['P1_y'],higher['P1_yy'],higher['P1_yyy']],
            raw_m15_y=[b.scale(density['P1_y'],g.neg(radius)),higher['m15_yy'],higher['m15_yyy']],
            intermediate=packet,own=dict(fields=fields,cumulative=cumulative))

    def ref(self,node):return source.FunctionRef(self.g,node)
    def jet(self,ids):return outer.FunctionJet(tuple(map(self.ref,ids)))
    def at_own(self,name,value):
        chart=self.charts[name];coordinate=chart['coordinate']
        return {group:{key:[outer.Algebra(self.g,len(outer.rows(jet))-1).at(jet,coordinate,value)
            for jet in jets] for key,jets in fields.items()} for group,fields in chart['own'].items()}
    def at_fraction(self,name,fraction):
        if name not in CHARTS:raise ValueError('Actual prescribed chart required')
        q=Fraction(fraction)
        if not 0<=q<=1:raise ValueError('Native fraction must lie in the original complete chart')
        chart=self.charts[name];lo,hi=chart['domain'];g=self.g
        value=g.add(lo,g.mul(g.constant(q),g.sub(hi,lo)));coordinate=chart['coordinate']
        at=lambda jets:[outer.Algebra(g,len(outer.rows(jet))-1).at(jet,coordinate,value) for jet in jets]
        return dict(chart=name,native_query=value,radius=outer.Algebra(g,0).at(outer.FunctionJet((chart['radius'],)),coordinate,value),
            physical={key:at(jets) for key,jets in chart['physical'].items()},
            N_y=at(chart['N_y']),pressure_y=at(chart['pressure_y']),m15_y=at(chart['m15_y']),
            own=self.at_own(name,value),numerical_oracle_available=False)
    def _snapshot(self):
        return outer.encode(dict(backgrounds=self.backgrounds,downstream=self.downstream,loop=self.loop,modified=self.modified,
            N=self.N,background_y=self.background_y,charts=self.charts,requests=self.requests,
            contract=self.contract,macro_first_join=self.macro_first_join,definition=self.definition))
    def assert_graph(self):
        self.base.assert_graph()
        if (self.snapshot!=self._snapshot() or self.g.nodes[:len(self.prefix)]!=self.prefix
            or self.hashes!=self._hashes or self.numerical_oracle_available):raise ValueError('Actual prescribed high sources changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed prescribed high input '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        if (not rec['all_passed'] or len(rec.get('checks',{}))!=7
            or not all(rec[key] for key in GATES) or any(rec[key] for key in OPEN)
            or pins.pop(NAME)!=sha(NAME) or pins.pop(Path(__file__).stem+'_check.py')!=sha(Path(__file__).stem+'_check.py')
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['prescribed_high_radial_source_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual prescribed high sources required')
        self.acceptance_loaded=True


def run(base=None):
    began=time.monotonic();owner=CurrentOriginalPrescribedHighRadialSourceFunctions(base,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        prescribed_high_radial_source_bindings=owner.snapshot,**dict.fromkeys(GATES+OPEN,False),
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_TEN_PRESCRIBED_FULL_CHART_HIGH_SOURCES_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw
