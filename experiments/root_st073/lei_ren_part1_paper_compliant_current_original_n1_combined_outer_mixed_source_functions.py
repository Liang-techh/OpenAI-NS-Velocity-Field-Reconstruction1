"""Finite actual own mixed source across the original 2Rc/Rp/Rv seams.

This lifts the retained pressure primitives by their proved FTC, rather
than replacing positive-order pressure with an unrelated leading output.
Numerical bounds, all-order smoothness and global n1 are separate gates.
"""
import copy
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_compact_higher_radial_pressure_C3_functions as previous

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer=previous.outer
NAME=PREFIX+'current_original_n1_combined_outer_mixed_source_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_combined_outer_mixed_source_functions_check.json'
GATES=('current_original_n1_finite_outer_mixed_derivative_contract_installed',
    'current_original_n1_2Rc_Rp_Rv_actual_own_C3_QV_C2_trace_joins_installed',
    'current_original_n1_2Rc_Rp_Rv_radial_pressure_join_C3_installed',
    'current_original_n1_finite_outer_cumulative_y3_C3_source_joins_installed',
    'current_original_n1_internal_pulse_own_y4_pressure_cumulative_y3_mixed_source_joins_installed',
    'current_original_n1_reserved_to_pulse_own_source_accessor_installed',
    'current_original_n1_full_outer_radial_mixed_source_joins_installed')
OPEN=tuple(dict.fromkeys(tuple(k for k in previous.OPEN if k not in GATES)+(
    'current_original_n1_outer_all_order_smoothness_quantitative_bounds_installed',
    'current_original_n1_outer_own_source_numeric_bounds_installed',
    'current_original_n1_Rin_rminus_Rc_own_field_joins_installed',
    'current_original_n1_global_field_join_C3_installed')))
STATE=('Uz1','Utheta1','F1','M1','mean1','K1','Q1','V1')
MOMENTS=('m11','m12','m13','m14','m15')
CONTRACT=tuple(dict(group='fields',field=key,physical_y_order=4,ordinary_Z_order=2 if key in ('Q1','V1') else 3)
    for key in STATE)+(dict(group='fields',field='P1',physical_y_order=3,ordinary_Z_order=3),)+tuple(
    dict(group='cumulative',field=key,physical_y_order=3,ordinary_Z_order=3) for key in MOMENTS)
ALIASES={'M1_y':'M1','mean1_y':'mean1','Q1_y':'Q1','V1_y':'V1'}
ROWS_PER_SEAM=sum((row['physical_y_order']+1)*(row['ordinary_Z_order']+1) for row in CONTRACT)
SEAMS=('2Rc','Rp','main_gap','gap_end','Rv')


def weighted_product(g,rows,weight,rate,n):
    """D_y^n(weight * product), with weight_y=rate*weight."""
    a=outer.Algebra(g,3)
    product=[previous.product_row(a,rows[0],rows[1],j) for j in range(n+1)]
    return a.scale(a.add(*(a.scale(product[j],g.mul(g.constant(math.comb(n,j)),
        previous.previous.power(g,rate,n-j))) for j in range(n+1))),weight)


class CurrentOriginalN1CombinedOuterMixedSourceFunctions:
    def __init__(self,compact=None,require_checked=True):
        self.compact=compact if compact is not None else previous.CurrentOriginalCompactHigherRadialPressureC3Functions()
        if not self.compact.acceptance_loaded:raise ValueError('Accepted actual compact/own higher source required')
        self.compact.assert_graph();self.high=self.compact.high;self.radial=self.high.joins
        self.own=self.compact.own_outer;self.g=g=self.compact.g
        self.family=copy.deepcopy(self.compact.family);self.hashes=dict(self.compact.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):
            self.hashes[name]=sha(name)
        a=outer.Algebra(g,3);v=outer.Algebra(g,2);self.charts={}
        c=self.compact;b=c.repaired
        cumulative={key:[self.own.jet(b['cumulative'][key])] for key in MOMENTS}
        cumulative['m11']=c.own_reserved['M1'][:4]
        cumulative['m13']=list(c.own_pressure.values())
        cumulative['m15'] += [a.scale(c.own_pressure['P1_y'],g.neg(c.own_R)),
            c.own_m15_higher['m15_yy'],c.own_m15_higher['m15_yyy']]
        one=[a.fixed(1),a.fixed(0),a.fixed(0)]
        weight12=g.mul(g.constant(2),c.own_R,c.own_R)
        weight14=g.mul(c.own_R,g.unary('sqrt',g.mul(g.constant(2),c.own_R)))
        for n in range(3):
            cumulative['m12'].append(weighted_product(g,[c.own_reserved['F1'],one],weight12,g.constant(2),n))
            cumulative['m14'].append(weighted_product(g,[c.theta0_radial,c.own_reserved['Uz1']],weight14,g.constant('3/2'),n))
        self.charts['reserved']=dict(coordinate=c.phase,domain=[g.zero,g.unary('log',g.constant(2))],
            radius=c.own_R,physical_y_Jacobian=g.one,fields={**c.own_reserved,'P1':list(c.own_pressure.values())},
            cumulative=cumulative,pressure_lift='actual retained compact primitive plus signed patch FTC')
        for name,states in c.own_outside.items():
            if name=='quiet':
                packet=self.own.charts[name];coordinate=packet['coordinate'];domain=packet['domain'];R=packet['radius']
                Py=self.own.quiet_pressure_y;high=self.high.quiet['higher_density'];jac=g.one
            else:
                packet=self.radial.current_own[name];meta=self.high.charts[name]
                coordinate=packet['coordinate'];domain=packet['domain'];R=packet['radius'];jac=meta['physical_y_Jacobian']
                Py=packet['fields']['P1_y'];high=meta['higher_density']
            pressure=[packet['fields']['P1'],Py,high['P1_yy'],high['P1_yyy']]
            cumulative={key:[packet['cumulative'][key]]+[a.fixed(0) for _ in range(3)] for key in MOMENTS}
            cumulative['m11']=states['M1'][:4];cumulative['m13']=pressure
            cumulative['m15']=[packet['cumulative']['m15'],a.scale(Py,g.neg(R)),
                a.scale(a.add(Py,pressure[2]),g.neg(R)),
                a.scale(a.add(Py,a.scale(pressure[2],2),pressure[3]),g.neg(R))]
            self.charts[name]=dict(coordinate=coordinate,domain=domain,radius=R,physical_y_Jacobian=jac,
                fields={**states,'P1':pressure},cumulative=cumulative,
                pressure_lift='same retained own cumulative plus actual fixed-pressure FTC; no patch outside support')
        # The tail is the actual saved mass tail at its translated coordinate.
        # Its nonzero terminal functions are kept even when a separate repair
        # theorem ultimately identifies them with zero.
        tail=self.own.inputs[outer.intermediate.previous.NAME]
        self.tail_coordinate=self.own.ref(tail['pressure_coordinate'])
        translated=g.add(self.tail_coordinate,self.own.ref(tail['tail_origins']['pressure_log_radius']),
            g.neg(self.own.ref(tail['tail_origins']['mass_log_radius'])))
        self.tail_mass_coordinate=translated
        mass={key:a.at(self.own.jet(ids),self.own.ref(tail['mass_coordinate']),translated)
            if len(ids)==4 else v.at(self.own.jet(ids),self.own.ref(tail['mass_coordinate']),translated)
            for key,ids in tail['mass_tail'].items()}
        R=g.unary('exp',g.add(self.own.ref(tail['tail_origins']['pressure_log_radius']),self.tail_coordinate))
        fields={key:[mass[key]]+[(v if key in ('Q1','V1') else a).fixed(0) for _ in range(4)] for key in STATE}
        for key in ('mean1','K1','Q1'):
            alg=v if key=='Q1' else a;fields[key]=[alg.scale(mass[key],(-1)**j) for j in range(5)]
        terminal=self.own.jets(tail['terminal']);fields['P1']=[terminal['m13']]+[a.fixed(0) for _ in range(3)]
        self.charts['tail']=dict(coordinate=self.tail_coordinate,domain=[g.zero,None],radius=R,physical_y_Jacobian=g.one,
            fields=fields,cumulative={key:[terminal[key]]+[a.fixed(0) for _ in range(3)] for key in MOMENTS},
            pressure_lift='same saved terminal pressure constant; same saved mass tail coordinate translation')
        self.contract=copy.deepcopy(CONTRACT);self.aliases=dict(ALIASES)
        self.joins={seam:dict(left=self.at(left,lv),right=self.at(right,rv)) for seam,left,lv,right,rv in (
            ('2Rc','reserved',self.charts['reserved']['domain'][1],'quiet',g.zero),
            ('Rp','quiet',self.charts['quiet']['domain'][1],'main',self.charts['main']['domain'][0]),
            ('main_gap','main',self.charts['main']['domain'][1],'gap',self.charts['gap']['domain'][0]),
            ('gap_end','gap',self.charts['gap']['domain'][1],'end',self.charts['end']['domain'][0]),
            ('Rv','end',self.charts['end']['domain'][1],'tail',g.zero))}
        self.definition=dict(finite_physical_y4_state_y3_pressure_cumulative_contract=True,
            ordinary_Z3_state_pressure_cumulative_QV_Z2=True,original_both_sides_and_nonzero_histories_preserved=True,
            actual_own_pressure_FTC_lift_required_before_higher_leading_derivative_composition=True,
            physical_y_main_gap_mu_end_quiet_reserved_tail_one=True,
            numerical_bounds_arbitrary_order_smoothness_inner_joins_global_n1_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self._owners=(self.compact,self.high,self.radial,self.own,g)
        self.acceptance_loaded=False;self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def at(self,name,value):
        chart=self.charts[name];x=chart['coordinate']
        group=lambda packet:{key:[outer.Algebra(self.g,len(outer.rows(jet))-1).at(jet,x,value) for jet in jets]
            for key,jets in packet.items()}
        return dict(fields=group(chart['fields']),cumulative=group(chart['cumulative']),
            radius=outer.Algebra(self.g,3).at(outer.Algebra(self.g,3).fixed(chart['radius']),x,value).value)
    def _snapshot(self):
        return outer.encode(dict(charts=self.charts,joins=self.joins,contract=self.contract,aliases=self.aliases,
            tail_coordinate=self.tail_coordinate,tail_mass_coordinate=self.tail_mass_coordinate,definition=self.definition))
    def assert_graph(self):
        self.compact.assert_graph()
        if (self._owners!=(self.compact,self.high,self.radial,self.own,self.g) or self.hashes!=self._hashes
            or self.snapshot!=self._snapshot() or self.g.nodes[:len(self.prefix)]!=self.prefix
            or self.numerical_oracle_available):raise ValueError('Actual combined outer source changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed combined outer source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        if (not rec['all_passed'] or not all(rec[key] for key in GATES) or any(rec[key] for key in OPEN)
            or pins.pop(NAME)!=sha(NAME) or pins.pop(Path(__file__).stem+'_check.py')!=sha(Path(__file__).stem+'_check.py')
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['combined_outer_mixed_source_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked combined outer mixed source required')
        self.acceptance_loaded=True
    def at_fraction(self,name,fraction):
        self.assert_graph()
        if not self.acceptance_loaded:raise ValueError('Checked combined outer source required')
        if name not in self.charts or name=='tail':raise ValueError('Choose a finite original own chart')
        if not isinstance(fraction,(int,str,Fraction)) or isinstance(fraction,bool):raise TypeError('Exact rational fraction required')
        q=Fraction(fraction)
        if not 0<=q<=1:raise ValueError('Chart fraction must lie in [0,1]')
        lo,hi=self.charts[name]['domain'];value=self.g.add(lo,self.g.mul(self.g.constant(q),self.g.sub(hi,lo)))
        return dict(chart=name,exact_native_fraction=str(q),**self.at(name,value),
            finite_outer_source_accepted=True,numerical_oracle_available=False,global_n1=False)
    def at_tail_log_offset(self,offset):
        self.assert_graph()
        if not self.acceptance_loaded:raise ValueError('Checked combined outer source required')
        if not isinstance(offset,(int,str,Fraction)) or isinstance(offset,bool):raise TypeError('Exact rational offset required')
        q=Fraction(offset)
        if q<0:raise ValueError('Tail log(R/Rv) must be nonnegative')
        return dict(chart='tail',exact_log_R_over_Rv=str(q),**self.at('tail',self.g.constant(q)),
            finite_outer_source_accepted=True,numerical_oracle_available=False,global_n1=False)


def run(compact=None):
    began=time.monotonic();owner=CurrentOriginalN1CombinedOuterMixedSourceFunctions(compact,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,input_hashes=owner.hashes,
        combined_outer_mixed_source_bindings=owner.snapshot,**dict.fromkeys(GATES+OPEN,False),
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_COMBINED_OUTER_MIXED_SOURCE_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw
