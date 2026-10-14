"""Actual own field source assembly across 2Rc, Rp and Rv.

Keep both traces and all pressure histories. Ordinary trace equality and
radial pressure/mixed matching are separate acceptance requirements.
"""
import copy
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_Rp_native_C5_pulse_source as previous
import lei_ren_part1_paper_compliant_current_original_n1_intermediate_own_C3_chart_functions as intermediate

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
FunctionJet,rows,encode,Algebra=previous.previous.FunctionJet,previous.previous.rows,previous.previous.encode,previous.previous.Algebra
NAME=PREFIX+'current_original_n1_outer_own_source_joins.json.gz'
RECEIPT=PREFIX+'current_original_n1_outer_own_source_joins_check.json'
CONSTRUCTION_RECEIPT=PREFIX+'current_original_n1_outer_own_source_joins_construction_check.json'
GATES=('current_original_n1_2Rc_Rp_Rv_actual_own_C3_QV_C2_trace_joins_installed',
       'current_original_n1_reserved_to_pulse_own_source_accessor_installed')
OPEN=tuple(dict.fromkeys(tuple(k for k in previous.OPEN if k not in GATES)+(
    'current_original_n1_2Rc_Rp_Rv_radial_pressure_join_C3_installed',
    'current_original_n1_Rin_rminus_Rc_own_field_joins_installed',
    'current_original_n1_outer_own_source_numeric_bounds_installed')))
FIELD_KEYS=('Uz1','Utheta1','F1','M1','mean1','K1','P1','Q1','V1','M1_y','mean1_y','Q1_y','V1_y')


def source_operands(row):
    """Closed references, including the actual implicit/Volterra input maps."""
    op=row['operation']
    if op in previous.SCHEMA:
        out=[]
        for key in previous.SCHEMA[op]:
            value=row.get(key)
            if value is not None:out.extend(value if isinstance(value,list) else [value])
        return out
    def packed(value):
        if isinstance(value,dict):return [i for v in value.values() for i in packed(v)]
        if isinstance(value,list):return [i for v in value for i in packed(v)]
        if type(value) is not int:raise ValueError('Non-reference in actual source input map')
        return [value]
    simple={'current_original_axis_pressure_function_derivative':('original_function','Z_variable'),
        'current_original_unique_analytic_core_component':('rho','Z_variable'),
        'original_Rm_unique_leading_implicit_control':('system','Z_variable'),
        'original_compact_unique_Banach_limit_component':('system','Z_variable')}
    if op in simple:return [row[key] for key in simple[op]]
    if op=='current_original_unique_regular_n1_component':
        out=[row[key] for key in ('coordinate','Z_variable','coordinate_domain_lower','coordinate_domain_upper')]
        raw=row['actual_raw_operator_bindings'];assert set(raw)=={'core','first'}
        for branch in raw.values():
            assert set(branch)=={'B0','B1','C0','g','N1theta','N1z','N1p','x','rho','D'}
            assert branch['D']==[0,0,2,0,3,1]
            for key in ('B0','B1','C0','g','N1theta','N1z','N1p','x','rho'):out.extend(packed(branch[key]))
        return out
    if op=='original_Rm_unique_leading_implicit_system':
        return [i for key in ('exact_weight_nodes','exact_defect_value_nodes','invAm2',
            'exact_Z_independent_row_and_control_units') for i in packed(row[key])]
    if op=='original_compact_unique_Banach_limit_system':
        return [i for key in ('variables','exact_map','actual_target_values','actual_selected_N',
            'original_weights','initial_vector') for i in packed(row[key])]
    raise ValueError('Unproved own coefficient operand '+op)


class ReservedRadialDerivative(previous.pulse3.OrdinaryZDerivative):
    """Actual first radial derivative of the original selected bump source."""
    def __init__(self,g,phase,coefficients):
        super().__init__(g)
        self.derivatives[self.z.node]=g.zero;self.z=phase;self.derivatives[phase.node]=g.one
        for jet in coefficients.values():
            for ref in rows(jet):self.seed(ref,g.zero)
    def rule(self,i):
        row=self.g.nodes[i]
        if row['operation']=='compact_raw_beta':
            dx=self.derivatives[row['argument']]
            if dx==self.g.zero:return self.g.zero
            first,_=previous.compact.leading.beta_first_two(self.g,
                source.FunctionRef(self.g,row['argument']),source.FunctionRef(self.g,i))
            return self.g.mul(first,dx)
        return super().rule(i)


def own_Q(g,U,mean,delta):
    a=Algebra(g,2);z=FunctionJet((g.symbol('Z'),g.one,g.zero));one=a.fixed(1)
    d=a.fixed(delta);L=a.sub(one,a.mul(d,a.mul(z,z)))
    numerator=a.sub(a.scale(a.mul(z,a.trim(U)),2),
        a.add(a.mul(a.add(one,d),a.mul(z,a.trim(mean))),a.mul(a.sub(one,a.mul(z,z)),a.shift(mean))))
    return a.div(numerator,L,'original positive L=1-delta*Z^2 on common axial domain')


class CurrentOriginalN1OuterOwnSourceJoins:
    def __init__(self,inlet=None,require_checked=True):
        global FunctionJet,rows,encode,Algebra
        self.inlet=inlet if inlet is not None else previous.CurrentOriginalRpNativeC5PulseSource()
        if not self.inlet.acceptance_loaded:raise ValueError('Checked full actual Rp source required')
        FunctionJet,rows,encode,Algebra=previous.previous.FunctionJet,previous.previous.rows,previous.previous.encode,previous.previous.Algebra
        self.inlet.assert_graph();self.g=g=self.inlet.g;self.family=copy.deepcopy(self.inlet.family)
        self.hashes=dict(self.inlet.hashes);self.inputs={}
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):
            self.hashes[name]=sha(name)
        for module,key in ((intermediate,'intermediate_own_C3_chart_bindings'),
            (intermediate.previous,'repaired_own_C3_terminal_bindings'),
            (intermediate.previous.previous,'complete_C3_defect_repair_bindings')):
            self.inputs[module.NAME]=self.pin(module)[key]
        b=self.inputs[intermediate.previous.NAME];i=self.inputs[intermediate.NAME]
        c=self.inputs[intermediate.previous.previous.NAME];q=self.inlet.inputs[previous.quiet.NAME]
        self.delta=self.inlet.geometry['delta'];self.reserved_phase=self.ref(c['repair_velocity']['log_R_over_Rc'])
        self.reserved_radius=self.ref(c['repair_velocity']['physical_radius'])
        self.coefficients=self.jets(c['coefficients']);a=Algebra(g,3);v=Algebra(g,2)
        for key in ('Uz1','Utheta1','F1'):
            assert b['fields'][key]==c['repair_velocity'][key]
        assert b['patch_bodies']['physical_F1']==c['repair_velocity']['F1']
        self.coefficient_phase_independence=[]
        for key,jet in self.coefficients.items():
            for n,ref in enumerate(rows(jet)):
                closure=self.operand_closure(ref.node)
                if self.reserved_phase.node in closure or self.reserved_radius.node in closure:
                    raise ValueError('Reserved coefficient depends on the own query phase/radius')
                self.coefficient_phase_independence.append(dict(coefficient=key,order=n,
                    actual_source_node=ref.node,closed_source_nodes=len(closure),no_reserved_query_phase_or_radius=True))
        assert self.reserved_radius.node in self.operand_closure(b['fields']['mean1'][0])
        assert self.reserved_phase.node in self.operand_closure(c['repair_velocity']['Uz1'][0])
        derivative=ReservedRadialDerivative(g,self.reserved_phase,self.coefficients)
        self.reserved_fields=self.jets(b['fields'])
        self.reserved_velocity_y={key:FunctionJet(tuple(derivative.derivative(ref) for ref in rows(self.reserved_fields[key])))
            for key in ('Uz1','Utheta1','F1')}
        Qy=own_Q(g,self.reserved_velocity_y['Uz1'],self.reserved_fields['mean1_y'],self.delta)
        self.reserved_fields['Q1_y']=Qy
        self.reserved_fields['V1_y']=v.scale(v.add(self.reserved_fields['Q1'],Qy),self.reserved_radius)
        self.charts={'reserved':dict(coordinate=self.reserved_phase,domain=[g.zero,g.unary('log',g.constant(2))],
            radius=self.reserved_radius,fields={k:self.reserved_fields[k] for k in FIELD_KEYS},
            cumulative=self.jets(b['cumulative']),physical_log_radius_jacobian=g.one)}
        for name in ('quiet','pulse_main','pulse_gap','pulse_end'):
            packet=i['charts'][name]
            self.charts[name]=dict(coordinate=self.ref(packet['coordinate']),domain=[self.ref(n) for n in packet['domain']],
                radius=self.ref(packet['radius']),fields=self.jets(packet['fields']),cumulative=self.jets(packet['cumulative']),
                physical_log_radius_jacobian=self.inlet.geometry['mu'] if name in ('pulse_main','pulse_gap') else g.one)
        self.quiet_pressure_y=a.scale(self.jet(q['N1p']),self.ref(q['R']))
        self.pressure_radial_traces={'2Rc':dict(left=a.at(self.jet(b['total_pressure_y']),self.reserved_phase,g.unary('log',g.constant(2))),
            right=a.at(self.quiet_pressure_y,self.ref(q['s']),g.zero),matching_not_yet_accepted=True)}
        self.joins={
            '2Rc':dict(left=self.at('reserved',self.charts['reserved']['domain'][1]),right=self.decode_trace(i['endpoints']['quiet']['left'])),
            'Rp':dict(left=self.decode_trace(i['endpoints']['quiet']['right']),right=self.decode_trace(i['endpoints']['pulse_main']['left']))}
        tail_radius=g.unary('exp',self.ref(b['tail_origins']['pressure_log_radius']))
        mass_y=g.sub(self.ref(b['tail_origins']['pressure_log_radius']),self.ref(b['tail_origins']['mass_log_radius']))
        mass={key:Algebra(g,len(rows(jet))-1).at(jet,self.ref(b['mass_coordinate']),mass_y)
            for key,jet in self.jets(b['mass_tail']).items()}
        mass.update(P1=self.jet(b['pressure_tail']['P1']),M1_y=a.fixed(0),mean1_y=a.neg(mass['mean1']),
            Q1_y=v.neg(mass['Q1']),V1_y=v.fixed(0))
        self.tail=dict(fields=mass,cumulative=self.jets(b['terminal']),radius=tail_radius)
        assert set(FIELD_KEYS)<=set(mass)
        self.joins['Rv']=dict(left=self.decode_trace(i['endpoints']['pulse_end']['right']),right=self.tail)
        self.definition=dict(actual_saved_both_sides_not_replaced_by_common_zero=True,
            own_recovery_one_plus_delta=True,fixed_leading_N1p_retained_separately=True,
            existing_C3_QV_C2_source_functions_and_exact_physical_chart_domains=True,
            original_reserved_first_radial_bump_derivatives_constructed=True,
            radial_pressure_matching_mixed_bounds_other_joins_global_n1_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self._family=copy.deepcopy(self.family);self._owners=(self.inlet,g)
        self.numerical_oracle_available=False;self.acceptance_loaded=False
        if require_checked:self.load_acceptance()
    def ref(self,i):return source.FunctionRef(self.g,i)
    def jet(self,ids):return FunctionJet(tuple(self.ref(n) for n in ids))
    def jets(self,packet):return {k:self.jet(v) for k,v in packet.items()}
    def decode_trace(self,packet):
        return dict(fields=self.jets(packet['fields']),cumulative=self.jets(packet['cumulative']),radius=self.ref(packet['radius']))
    def operand_closure(self,node):
        seen=set();stack=[node]
        while stack:
            i=stack.pop()
            if i in seen:continue
            if type(i) is not int or not 0<=i<len(self.g.nodes):raise ValueError('Invalid own source operand')
            seen.add(i);stack.extend(source_operands(self.g.nodes[i]))
        return seen
    def pin(self,module):
        data=parent.linear.read(module.NAME);rec=parent.linear.read(module.RECEIPT)
        if (not rec['all_passed'] or not all(rec[k] for k in module.GATES) or any(rec[k] for k in module.OPEN)
            or rec['source_family']!=self.family or data['source_family']!=self.family
            or data['exact_graph_nodes']!=self.g.nodes[:len(data['exact_graph_nodes'])]):
            raise ValueError('Checked original own source required '+module.NAME)
        pins=dict(rec['input_hashes']);archive=pins.pop(module.NAME);checker=Path(module.__file__).stem+'_check.py'
        if archive!=sha(module.NAME) or pins.pop(checker)!=sha(checker) or pins!=data['input_hashes']:
            raise ValueError('Exact own source archive/checker pins required')
        for name,digest in rec['input_hashes'].items():
            if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed own source '+name)
            self.hashes[name]=digest
        for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self.hashes[name]=sha(name)
        return data
    def at(self,name,value):
        chart=self.charts[name];x=chart['coordinate']
        group=lambda packet:{k:Algebra(self.g,len(rows(q))-1).at(q,x,value) for k,q in packet.items()}
        return dict(fields=group(chart['fields']),cumulative=group(chart['cumulative']),
            radius=Algebra(self.g,3).at(Algebra(self.g,3).fixed(chart['radius']),x,value).value)
    def _snapshot(self):
        return encode(dict(inputs=self.inputs,coefficients=self.coefficients,delta=self.delta,
            coefficient_phase_independence=self.coefficient_phase_independence,
            reserved_phase=self.reserved_phase,reserved_radius=self.reserved_radius,
            reserved_velocity_y=self.reserved_velocity_y,charts=self.charts,joins=self.joins,tail=self.tail,
            quiet_pressure_y=self.quiet_pressure_y,pressure_radial_traces=self.pressure_radial_traces,definition=self.definition))
    def assert_graph(self):
        self.inlet.assert_graph()
        if (self.hashes!=self._hashes or self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot
            or self.family!=self._family or self._owners!=(self.inlet,self.g) or self.numerical_oracle_available is not False):
            raise ValueError('Actual own seam functions changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed own seam source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT)
        pins=dict(rec['input_hashes']);archive=pins.pop(NAME);checker=Path(__file__).stem+'_check.py'
        if (not rec['all_passed'] or not all(rec[k] for k in GATES) or any(rec[k] for k in OPEN)
            or archive!=sha(NAME) or pins.pop(checker)!=sha(checker) or pins!=self.hashes
            or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['outer_own_source_join_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual own seam source required')
        self.acceptance_loaded=True
    def at_fraction(self,name,fraction):
        self.assert_graph()
        if not self.acceptance_loaded:raise ValueError('Checked seam source required')
        if name not in self.charts:raise ValueError('Choose a saved reserved/quiet/pulse chart')
        if not isinstance(fraction,(int,str,Fraction)) or isinstance(fraction,bool):raise TypeError('Exact rational fraction required')
        q=Fraction(fraction)
        if not 0<=q<=1:raise ValueError('Chart fraction must lie in [0,1]')
        lo,hi=self.charts[name]['domain'];value=self.g.add(lo,self.g.mul(self.g.constant(q),self.g.sub(hi,lo)))
        return dict(chart=name,exact_fraction=str(q),**self.at(name,value),numerical_oracle_available=False,global_n1=False)


def run(inlet=None):
    began=time.monotonic();owner=CurrentOriginalN1OuterOwnSourceJoins(inlet,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,outer_own_source_join_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_OUTER_OWN_SOURCE_JOINS_DEFINED',len(owner.g.nodes),flush=True);return owner,raw
