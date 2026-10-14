"""Actual current-source Rp/Rv pressure-y joins and own pulse transport.

Ordinary outer traces remain accepted separately. Higher mixed matching,
original-scale numerical bounds and the global n1 accessor remain open.
"""
import copy
import gzip
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions as previous

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
outer=previous.previous
rp=outer.previous
NAME=PREFIX+'current_original_n1_outer_radial_pressure_source_joins.json.gz'
RECEIPT=PREFIX+'current_original_n1_outer_radial_pressure_source_joins_check.json'
GATES=('current_original_n1_Rp_Rv_total_radial_pressure_C3_source_joins_installed',
    'current_original_n1_Rp_Rv_actual_leading_radial_recovery_source_joins_installed',
    'current_original_n1_pulse_canonical_current_partial_pressure_source_identity_installed')
OPEN=tuple(dict.fromkeys(tuple(k for k in previous.OPEN if k not in GATES)+(
    'current_original_n1_full_outer_radial_mixed_source_joins_installed',)))


def own_zero_support_fields(g,cumulative,R,delta):
    a=outer.Algebra(g,3);b=outer.Algebra(g,2);zero=a.fixed(0)
    mean=a.scale(cumulative['m11'],g.quotient(g.one,R,'positive actual own pulse radius'))
    Q=outer.own_Q(g,zero,mean,delta)
    return dict(Uz1=zero,Utheta1=zero,F1=zero,M1=cumulative['m11'],mean1=mean,K1=mean,
        P1=cumulative['m13'],Q1=Q,V1=b.scale(Q,R),M1_y=zero,mean1_y=a.neg(mean),
        Q1_y=b.neg(Q),V1_y=b.fixed(0))


class CurrentOriginalN1OuterRadialPressureSourceJoins:
    def __init__(self,pulse=None,require_checked=True):
        self.pulse=pulse if pulse is not None else previous.CurrentOriginalPulseRadialPressureC3Functions()
        if not self.pulse.acceptance_loaded:raise ValueError('Checked actual pulse radial/FTC source required')
        self.pulse.assert_graph();self.outer=self.pulse.outer;self.inlet=self.pulse.inlet
        self.g=g=self.pulse.g;self.family=copy.deepcopy(self.pulse.family);self.hashes=dict(self.pulse.hashes)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).stem+'_check.py'):
            self.hashes[name]=sha(name)
        self.quiet=copy.deepcopy(self.inlet.inputs[rp.quiet.NAME])
        self.intermediate=copy.deepcopy(self.outer.inputs[outer.intermediate.NAME])
        self.canonical=copy.deepcopy(self.inlet.inputs[rp.pulse5.NAME])
        rec=parent.linear.read(outer.CONSTRUCTION_RECEIPT)
        if (not rec['all_three_ordinary_axial_trace_identities_passed'] or rec['ordinary_own_field_and_cumulative_trace_rows']!=204
            or rec['source_family']!=self.family):raise ValueError('Accepted actual 204 ordinary trace rows required')
        self.radius_proof=copy.deepcopy(rec['checks']['actual_old_current_radius_sources'])
        if not self.radius_proof['passed']:raise ValueError('Accepted actual old/current radius proof required')
        a=outer.Algebra(g,3);q=self.quiet
        self.quiet_pressure_y=self.outer.quiet_pressure_y
        quiet_right=a.at(self.quiet_pressure_y,self.ref(q['s']),self.ref(q['length']))
        self.quiet_N1p_right=a.at(self.jet(q['N1p']),self.ref(q['s']),self.ref(q['length']))
        self.quiet_physical_right=self.jets(q['Rp_endpoint']['physical'])
        self.pressure_traces=dict(Rp=dict(left=quiet_right,right=self.pulse.endpoints['main']['left']['density']['P1_y']),
            Rv=dict(left=self.pulse.endpoints['end']['right']['density']['P1_y'],right=a.fixed(0)))
        self.N1p_traces=dict(Rp=dict(left=self.quiet_N1p_right,right=self.pulse.endpoints['main']['left']['density']['N1p']),
            Rv=dict(left=self.pulse.endpoints['end']['right']['density']['N1p'],right=a.fixed(0)))
        self.Rp_cumulative=dict(self.outer.joins['Rp']['left']['cumulative'])
        self.current_own={}
        for name,data in self.pulse.charts.items():
            cumulative=dict(self.Rp_cumulative)
            for key in ('m13','m15'):cumulative[key]=a.add(cumulative[key],data['signed_pressure_increment'][key])
            fields=own_zero_support_fields(g,cumulative,data['radius'],self.pulse.delta)
            fields['P1_y']=data['density']['P1_y']
            self.current_own[name]=dict(coordinate=data['coordinate'],radius=data['radius'],domain=data['domain'],
                fields=fields,cumulative=cumulative,physical_y_Jacobian=data['physical_y_Jacobian'])
        # Preserve the actual terminal pressure function as a constant
        # physical-y extension. Its derivative is zero; its value is not reset.
        b=self.outer.inputs[outer.intermediate.previous.NAME]
        self.tail_coordinate=self.ref(b['pressure_coordinate'])
        self.tail_pressure=self.outer.tail['fields']['P1']
        self.tail_pressure_definition=parent.active.current.ast_binding(outer.intermediate.previous.CurrentOriginalN1RepairedOwnC3TerminalFunctions.__init__)
        # Source projection maps are drafts until the checker establishes
        # actual input, radius, integral and finite-tensor equivalences.
        mapping=dict(self.inlet.input_substitution)
        for name,data in self.canonical['leading_fields'].items():
            for key,ids in data.items():
                for old,new in zip(ids,outer.rows(self.inlet.leading_fields[name][key])):
                    if old not in (0,1):mapping[old]=new
        for name,data in self.canonical['selected_pulse'].items():
            for old,new in zip(data,outer.rows(self.inlet.selected_pulse[name])):
                if old not in (0,1):mapping[old]=new
        self.field_projection=copy.copy(mapping);self.radius_projection={}
        main_radius=self.intermediate['charts']['pulse_main']['radius']
        radius_row=g.nodes[main_radius];assert radius_row['operation']=='analytic_unary' and radius_row['name']=='exp'
        log_row=g.nodes[radius_row['argument']];assert log_row['operation']=='sum'
        query=self.pulse.charts['main']['coordinate'].node
        import lei_ren_part1_paper_compliant_current_original_n1_outer_own_source_joins_check as source_check
        fixed=[i for i in log_row['arguments'] if not source_check.free_query_dependency(self.outer,i,query)]
        if len(fixed)!=1:raise ValueError('Exact original common pulse logRp source required')
        self.old_logRp=self.ref(fixed[0]);self.current_logRp=g.unary('log',self.inlet.geometry['Rp'])
        mapping[self.old_logRp.node]=self.current_logRp
        for name,data in self.pulse.charts.items():
            old=self.intermediate['pulse_bindings'][name]['radius']
            self.radius_projection[old]=data['radius'];mapping[old]=data['radius']
            raw_radius=self.intermediate['charts']['pulse_'+name]['radius']
            self.radius_projection[raw_radius]=data['radius'];mapping[raw_radius]=data['radius']
        rewrite=rp.PulseSourceSubstitution(g,mapping)
        self.projected_old={}
        for name,binding in self.intermediate['pulse_bindings'].items():
            self.projected_old[name]=dict(
                products={k:self.jet([rewrite(self.ref(i)).node for i in ids]) for k,ids in binding['products'].items()},
                integrands={k:self.jet([rewrite(self.ref(i)).node for i in ids]) for k,ids in binding['current_integrands'].items()},
                partial_integrals={k:self.jet([rewrite(self.ref(i)).node for i in ids]) for k,ids in binding['current_integrals'].items()},
                complete_integrals={k:self.jet([rewrite(self.ref(i)).node for i in ids]) for k,ids in self.canonical['integrals'][name].items()},
                upper={k:self.jet([rewrite(self.ref(i)).node for i in ids]) for k,ids in binding['upper'].items()},
                lower={k:self.jet([rewrite(self.ref(i)).node for i in ids]) for k,ids in binding['lower'].items()})
        # Complete old pressure map inputs are projected only after their
        # equality is checked; no whole pressure output is replaced.
        pressure_mapping=dict(mapping)
        for name,binding in self.intermediate['pulse_bindings'].items():
            data=self.pulse.charts[name]
            packets=((binding['current_integrals'],data['partial_integrals']),
                (self.canonical['integrals'][name],self.pulse.complete_integrals[name]),
                (binding['upper'],data['upper_tensors']),(binding['lower'],self.pulse.global_lower_tensors))
            for old,new in packets:
                for key,ids in old.items():
                    for i,ref in zip(ids,outer.rows(new[key])):
                        if i not in (0,1):pressure_mapping[i]=ref
        pressure_rewrite=rp.PulseSourceSubstitution(g,pressure_mapping)
        self.projected_old_pressure={name:{key:self.jet([pressure_rewrite(self.ref(i)).node for i in ids])
            for key,ids in self.intermediate['charts']['pulse_'+name]['fixed_pressure'].items()}
            for name in self.pulse.charts}
        self.definition=dict(actual_quiet_Rp_native_endpoint_mean_projection_required=True,
            actual_current_pulse_leading_Uyy_Qyy_N1p_C3_consumed=True,
            existing_204_ordinary_trace_rows_preserved=True,
            current_own_pressure_cumulative_uses_actual_partial_primitives=True,
            nonzero_Rp_own_mass_and_pressure_history_retained=True,
            constant_tail_pressure_value_not_reset=True,
            old_current_partial_source_identity_requires_inputs_radii_integrands_and_tensors=True,
            own_one_plus_delta_and_leading_one_minus_delta_separate=True,
            higher_mixed_numeric_full_outer_global_n1_recursion_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self._owners=(self.pulse,self.outer,self.inlet,g);self.acceptance_loaded=False
        self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def ref(self,i):return source.FunctionRef(self.g,i)
    def jet(self,ids):return outer.FunctionJet(tuple(self.ref(i) for i in ids))
    def jets(self,data):return {k:self.jet(v) for k,v in data.items()}
    def _snapshot(self):
        return outer.encode(dict(quiet=self.quiet,intermediate=self.intermediate,canonical=self.canonical,
            radius_proof=self.radius_proof,quiet_pressure_y=self.quiet_pressure_y,quiet_N1p_right=self.quiet_N1p_right,
            quiet_physical_right=self.quiet_physical_right,pressure_traces=self.pressure_traces,N1p_traces=self.N1p_traces,
            Rp_cumulative=self.Rp_cumulative,current_own=self.current_own,tail_coordinate=self.tail_coordinate,
            tail_pressure=self.tail_pressure,tail_pressure_definition=self.tail_pressure_definition,
            field_projection={str(k):v for k,v in self.field_projection.items()},
            radius_projection={str(k):v for k,v in self.radius_projection.items()},
            old_logRp=self.old_logRp,current_logRp=self.current_logRp,
            projected_old=self.projected_old,projected_old_pressure=self.projected_old_pressure,definition=self.definition))
    def assert_graph(self):
        self.pulse.assert_graph()
        if (self.hashes!=self._hashes or self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot
            or self._owners!=(self.pulse,self.outer,self.inlet,self.g) or self.numerical_oracle_available):
            raise ValueError('Actual radial pressure source join changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed radial source '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT);pins=dict(rec['input_hashes'])
        archive=pins.pop(NAME);checker=Path(__file__).stem+'_check.py'
        if (not rec['all_passed'] or not all(rec[k] for k in GATES) or any(rec[k] for k in OPEN)
            or archive!=sha(NAME) or pins.pop(checker)!=sha(checker) or pins!=self.hashes
            or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or raw['outer_radial_pressure_source_join_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]):
            raise ValueError('Checked actual radial-pressure source joins required')
        self.acceptance_loaded=True


def run(pulse=None):
    began=time.monotonic();owner=CurrentOriginalN1OuterRadialPressureSourceJoins(pulse,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,
        outer_radial_pressure_source_join_bindings=owner.snapshot,input_hashes=owner.hashes,
        **dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_RP_RV_RADIAL_PRESSURE_SOURCE_JOINS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw
