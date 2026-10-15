"""Seven original own-field interfaces from postrestore through O3 entry.

Reuse accepted complete high-radial charts, including the existing O3/Rc
owner. Endpoint packets contain actual source handles and nonzero histories.
"""
import copy
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rm_outer_high_radial_source_functions as previous

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer=previous.outer
NAME=PREFIX+'current_original_Rm_outer_mixed_source_joins.json.gz'
RECEIPT=PREFIX+'current_original_Rm_outer_mixed_source_joins_check.json'
JOINS=(('postrestore_patch','postrestore','Rm_patch'),
       ('patch_reference','Rm_patch','Rh_reference'),
       ('reference_slope','Rh_reference','O2_slope'),
       ('slope_axial','O2_slope','O2_axial'),
       ('axial_buffer','O2_axial','O2_buffer'),
       ('buffer_transition','O2_buffer','O3_transition'),
       ('transition_power','O3_transition','O3_power'))
CHARTS=('postrestore',)+previous.CHARTS+('O3_power',)
GATES=tuple('current_original_n1_'+key+'_full_mixed_source_join_installed' for key,*_ in JOINS)
OPEN=tuple(key for key in previous.OPEN if key not in GATES)


class CurrentOriginalRmOuterMixedSourceJoins:
    def __init__(self,base=None,require_checked=True):
        fresh=base is None
        self.base=base if base is not None else previous.CurrentOriginalRmOuterHighRadialSourceFunctions()
        if not self.base.acceptance_loaded:raise ValueError('Accepted complete Rm/outer high sources required')
        if fresh and require_checked:previous.previous._restore_checked_layout(self.base,NAME,RECEIPT)
        self.base.assert_graph();self.g=g=self.base.g;self.family=copy.deepcopy(self.base.family)
        self.hashes=dict(self.base.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        self.prescribed=self.base.base.base.base
        self.rc=self.base.lower.rin.rc
        self.charts=dict(postrestore=self.prescribed.charts['postrestore'],**self.base.charts)
        r=self.rc;d=r.leading['chart'];n0,n1,n2=r.density['N1p'],r.higher_density['N1p_y'],r.higher_density['N1p_yy']
        self.charts['O3_power']=dict(domain=r.left['domain'],coordinate=r.coordinate,radius=r.radius,
            jacobian=g.one,phase=r.ref(d['actual_phase']),
            physical={key:[r.physical[key+suffix] for suffix in ('','_y','_yy','_yyy','_yyyy')] for key in ('U','m','Q')},
            N_y=[n0,n1,n2],pressure_y=r.left['fields']['P1'][1:],m15_y=r.left['cumulative']['m15'][1:],
            intermediate=r.intermediate,own=dict(fields=r.left['fields'],cumulative=r.left['cumulative']))
        self.loop=dict(self.prescribed.loop,**self.base.loop,O3_power=d)
        self.modified=dict(self.prescribed.modified,**self.base.modified,O3_power=d)
        self.backgrounds=dict(self.prescribed.backgrounds,**self.base.backgrounds,O3_power=r.leading['background'])
        self.joins={}
        for key,left,right in JOINS:
            self.joins[key]=dict(left_chart=left,right_chart=right,
                left=self.at_endpoint(left,'right'),right=self.at_endpoint(right,'left'))
        # Differentiate the actual flat logistic source, not a zero-padded
        # endpoint packet. O2 axial derivatives here are native; the checked
        # nonconstant-J lift supplies the corresponding physical derivatives.
        logistic=previous.high_source.downstream_source
        self.flat_controls={}
        for name,arg in (('O2_slope',self.charts['O2_slope']['coordinate']),
                         ('O2_axial',g.sub(g.one,self.charts['O2_axial']['coordinate'])),
                         ('O3_transition',self.charts['O3_transition']['coordinate'])):
            D=self.base.factories[name];sigma=[logistic.sigma(g,arg)];prime=[logistic.sigma_prime(g,arg)]
            for _ in range(5):sigma.append(D.derivative(sigma[-1]))
            for _ in range(4):prime.append(D.derivative(prime[-1]))
            self.flat_controls[name]=dict(sigma=sigma,prime=prime)
        self.patch_bumps=[];D=self.base.factories['Rm_patch'];x=self.charts['Rm_patch']['coordinate']
        for gamma in self.base.gamma:
            refs=[gamma[0]]
            for _ in range(5):refs.append(g.mul(x,D.derivative(refs[-1])))
            self.patch_bumps.append(refs)
        self.contract=copy.deepcopy(self.base.contract)
        self.definition=dict(seven_actual_original_interfaces_and_1722_own_mixed_pairs=True,
            accepted_six_complete_high_sources_and_original_O3_power_Rc_owner_reused=True,
            unique_Rm_leading_branch_and_Rh_normalization_not_midpoint_data=True,
            complete_nonzero_correction_and_signed_pressure_memories_retained=True,
            original_variable_native_Jacobians_and_global_periodic_phase_retained=True,
            finite_source_joins_not_global_n1_numeric_bounds_or_temporal_recursion=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()

    def __getattr__(self,name):return getattr(self.base,name)
    def at_endpoint(self,name,edge):
        if name not in CHARTS or edge not in ('left','right'):raise ValueError('Actual original endpoint required')
        chart=self.charts[name];value=chart['domain'][edge=='right'];x=chart['coordinate'];g=self.g
        at=lambda jets:[outer.Algebra(g,len(outer.rows(jet))-1).at(jet,x,value) for jet in jets]
        return dict(native_query=value,radius=outer.Algebra(g,0).at(outer.FunctionJet((chart['radius'],)),x,value),
            physical={field:at(jets) for field,jets in chart['physical'].items()},N_y=at(chart['N_y']),
            own={group:{field:at(jets) for field,jets in fields.items()} for group,fields in chart['own'].items()},
            numerical_oracle_available=False)
    def _snapshot(self):
        return outer.encode(dict(charts=self.charts,joins=self.joins,flat_controls=self.flat_controls,patch_bumps=self.patch_bumps,
            contract=self.contract,definition=self.definition))
    def assert_graph(self):
        self.base.assert_graph()
        if (self._snapshot()!=self.snapshot or self.g.nodes[:len(self.prefix)]!=self.prefix
            or self.hashes!=self._hashes or self.numerical_oracle_available):raise ValueError('Actual Rm/outer interfaces changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed Rm/outer interface input '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        if (not rec['all_passed'] or len(rec.get('checks',{}))!=6 or not all(rec[key] for key in GATES)
            or any(rec[key] for key in OPEN) or rec['actual_own_mixed_rows']!=1722
            or pins.pop(NAME)!=sha(NAME) or pins.pop(Path(__file__).stem+'_check.py')!=sha(Path(__file__).stem+'_check.py')
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['Rm_outer_mixed_source_join_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual Rm/outer mixed source interfaces required')
        self.acceptance_loaded=True


def run(base=None):
    began=time.monotonic();owner=CurrentOriginalRmOuterMixedSourceJoins(base,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        Rm_outer_mixed_source_join_bindings=owner.snapshot,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_SEVEN_RM_OUTER_MIXED_SOURCE_INTERFACES_DEFINED',len(owner.prefix),flush=True)
    return owner,raw
