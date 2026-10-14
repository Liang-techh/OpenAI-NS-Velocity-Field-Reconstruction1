"""Actual current Rp C5 source callback and selected pulse graph consumer.

The callback returns the eight archived current functions, without choosing
interval centers or re-instantiating numerical ancestors. Source equality to
the canonical native inlet must be checked before accepting the consumer.
"""
import copy
import gzip
import importlib
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_lower_active_common_source_join_functions as previous
import lei_ren_part1_paper_compliant_current_original_lower_active_archive_recovery as recovery
import lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions as quiet
import lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions as compact
import lei_ren_part1_paper_compliant_current_original_true_outer_C6_C5_functions as outer
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions as pulse5
import lei_ren_part1_paper_compliant_current_original_n1_pulse_C3_functions as pulse3
import lei_ren_part1_paper_compliant_current_original_C3_Rp_native_frame as old_frame
import lei_ren_part1_paper_compliant_current_original_outer_leading_C3_identity as old_leading
import lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions as importer

source,parent,HERE,PREFIX,sha=previous.source,previous.parent,previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_Rp_native_C5_pulse_source.json.gz'
RECEIPT=PREFIX+'current_original_Rp_native_C5_pulse_source_check.json'
FRAME_KEYS=('u','m1','m2','X','energy','Mp','P0','pressure')
GATES=('current_original_true_Rp_C5_selected_native_pulse_inlet_bridge_installed',
    'current_original_actual_Rp_native_C5_all_eight_source_functions_consumed',
    'current_original_actual_Rp_C5_frame_equals_canonical_native_smooth_functions')
OPEN=tuple(dict.fromkeys(tuple(k for k in previous.OPEN if k not in GATES)+(
    'current_original_Rp_native_C5_interval_selected_callback_installed',
    'current_original_Rp_native_C5_original_scale_numeric_bounds_installed')))
SCHEMA={**importer.REFERENCE_KEYS,'mathematical_positive_infinity':(),
    'exact_original_source_heat_phi':('argument',)}


class PulseSourceSubstitution:
    """Exact expression substitution with a closed graph operand schema."""
    def __init__(self,g,bindings):
        self.g=g;self.bindings=dict(bindings);self.mapping=dict(bindings)
        if 0 in bindings or 1 in bindings:raise ValueError('Global constants cannot be inlet placeholders')
    def __call__(self,ref):
        g=self.g;g.ids([ref]);stack=[(ref.node,False)];visiting=set()
        while stack:
            i,ready=stack.pop()
            if i in self.mapping:continue
            row=g.nodes[i];op=row['operation']
            if op not in SCHEMA:raise ValueError('Unsupported selected pulse expression '+op)
            children=[]
            for key in SCHEMA[op]:
                value=row.get(key)
                if value is not None:children.extend(value if isinstance(value,list) else [value])
            if not ready:
                if i in visiting:raise ValueError('Cyclic selected pulse source')
                visiting.add(i);stack.append((i,True))
                stack.extend((n,False) for n in reversed(children) if n not in self.mapping)
                continue
            attrs=copy.deepcopy(row);attrs.pop('operation')
            for key in SCHEMA[op]:
                value=attrs.get(key)
                if value is not None:
                    attrs[key]=[self.mapping[n].node for n in value] if isinstance(value,list) else self.mapping[value].node
            self.mapping[i]=g.node(op,**attrs);visiting.remove(i)
        return self.mapping[ref.node]


class CurrentOriginalRpNativeC5PulseSource:
    def __init__(self,seam=None,require_checked=True):
        self.seam=seam if seam is not None else recovery.restore(True)
        if not self.seam.acceptance_loaded:raise ValueError('Accepted actual common source graph required')
        self.seam.assert_graph();self.g=g=self.seam.g
        self.family=copy.deepcopy(self.seam.family);self.hashes=dict(self.seam.hashes)
        self.inputs={}
        for module,key in ((quiet,'true_quiet_Rp_C5_C3_bindings'),(compact,'true_compact_C5_bindings'),
            (outer,'true_outer_C6_C5_bindings'),(pulse5,'C5_pulse_C3_pressure_bindings'),
            (pulse3,'function_bindings')):
            data=self.pin(module,True);self.inputs[module.NAME]=copy.deepcopy(data[key])
        for module in (old_frame,old_leading,old_frame.quiet):self.inputs[module.NAME]=self.pin(module,False)
        for name in (previous.NAME,previous.RECEIPT,Path(previous.__file__).name,previous.RECOVERY):self._pin(name,sha(name))
        helper=PREFIX+'current_original_Rp_native_C5_terminal_projection.py'
        self._pin(helper,sha(helper))
        q=self.inputs[quiet.NAME];self.current_frame={key:self.jet(q['Rp_frame'][key]) for key in FRAME_KEYS}
        if set(q['Rp_frame'])!=set(FRAME_KEYS):raise ValueError('All eight actual current Rp functions required')
        self.geometry={key:self.ref(q[key]) for key in ('N','mu','Pstar','delta','s','R','Rp','R2','Tw','length')}
        self.current_P0=self.jet(q['P0'])
        p=self.inputs[pulse3.NAME];h=self.inputs[pulse5.NAME]
        self.old_native={key:self.jet(p['native'][key]) for key in FRAME_KEYS}
        bindings={}
        for packet in (p['native'],h['native']):
            for key,old in packet.items():
                for n,i in enumerate(old):
                    if i in (0,1):continue
                    value=previous.rows(self.current_frame[key])[n]
                    if i in bindings and bindings[i]!=value:raise ValueError('Ambiguous native source replacement')
                    bindings[i]=value
        self.input_substitution=bindings;rewrite=PulseSourceSubstitution(g,bindings)
        self.selected_pulse={key:self.jet([rewrite(self.ref(i)).node for i in ids]) for key,ids in h['selected_pulse'].items()}
        self.leading_fields={chart:{key:self.jet([rewrite(self.ref(i)).node for i in ids])
            for key,ids in fields.items()} for chart,fields in h['leading_fields'].items()}
        self.pressure={key:self.jet([rewrite(self.ref(i)).node for i in ids]) for key,ids in h['pressure'].items()}
        self.pulse_coordinates={};self.absolute_pressure={}
        a=previous.Algebra(g,3);mu=self.geometry['mu'];P=self.geometry['Pstar']
        rate=g.add(g.one,g.mul(g.constant(2),mu))
        for chart in ('main','gap','end'):
            query=g.symbol('original_own_pulse_'+('end_s' if chart=='end' else chart+'_xi'))
            y=g.add(g.quotient(g.constant(13),mu,'positive original mu'),query) if chart=='end' else g.quotient(query,mu,'positive original mu')
            radius=g.mul(self.geometry['Rp'],g.unary('exp',y))
            mass=g.mul(y,g.unary('exprel',g.neg(g.mul(rate,y))))
            self.pulse_coordinates[chart]=dict(query=query,y=y,radius=radius,
                coordinate='s=log(R/Rv)' if chart=='end' else 'xi=mu*log(R/Rp)',
                original_mu_and_actual_current_Rp_used=True)
            self.absolute_pressure[chart]=a.scale(a.add(a.trim(self.current_frame['pressure']),
                a.scale(a.mul(a.trim(self.current_frame['u']),a.trim(self.current_frame['u'])),g.mul(g.constant('1/2'),mass))),g.mul(P,P))
        self.rewritten_nodes={str(i):ref.node for i,ref in rewrite.mapping.items()}
        self.definition=dict(actual_current_Rp_eight_C5_jets_are_direct_callback_inputs=True,
            selected_three_C5_coefficients_and_main_gap_end_C5_Q4_pressure_C3_are_source_substituted=True,
            canonical_native_function_identity_required_before_acceptance=True,
            normalized_frame_not_multiplied_by_Pstar=True,
            physical_velocity_Pstar_and_pressure_Pstar_squared_once=True,
            original_mu_phase_Gamma_kernels_and_nonzero_histories_retained=True,
            numerical_interval_selected_owner_not_created_or_modified=True,
            original_scale_bounds_global_joins_recursion_and_corrected_NS_open=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.snapshot=self._snapshot();self.prefix=copy.deepcopy(g.nodes);self._hashes=dict(self.hashes)
        self._family=copy.deepcopy(self.family);self.acceptance_loaded=False
        self.numerical_oracle_available=False
        if require_checked:self.load_acceptance()
    def ref(self,i):return source.FunctionRef(self.g,i)
    def jet(self,ids):return previous.FunctionJet(tuple(self.ref(i) for i in ids))
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed Rp source '+name)
        self.hashes[name]=digest
    def pin(self,module,common_graph):
        data=parent.linear.read(module.NAME);rec=parent.linear.read(module.RECEIPT)
        gates=getattr(module,'GATES',(getattr(module,'GATE',None),))
        if not rec['all_passed'] or not all(rec[k] for k in gates) or any(rec[k] for k in module.OPEN):
            raise ValueError('Accepted Rp dependency scope required '+module.NAME)
        if data['source_family']!=self.family or rec['source_family']!=self.family:raise ValueError('Same original source family required')
        if common_graph and data['exact_graph_nodes']!=self.g.nodes[:len(data['exact_graph_nodes'])]:
            raise ValueError('Exact accepted Rp graph prefix required '+module.NAME)
        pins=dict(rec['input_hashes']);archive=pins.pop(module.NAME,None)
        checker=Path(module.__file__).stem+'_check.py';checksha=pins.pop(checker,None)
        if archive!=sha(module.NAME) or checksha!=sha(checker) or pins!=data['input_hashes']:
            raise ValueError('Exact Rp dependency archive and checker pins required '+module.NAME)
        for name,digest in rec['input_hashes'].items():self._pin(name,digest)
        for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self._pin(name,sha(name))
        return data
    def _snapshot(self):
        return previous.encode(dict(current_frame=self.current_frame,current_P0=self.current_P0,
            geometry=self.geometry,old_native=self.old_native,input_substitution={str(i):ref for i,ref in self.input_substitution.items()},
            selected_pulse=self.selected_pulse,leading_fields=self.leading_fields,pressure=self.pressure,
            pulse_coordinates=self.pulse_coordinates,absolute_pressure=self.absolute_pressure,
            rewritten_nodes=self.rewritten_nodes,inputs=self.inputs,definition=self.definition))
    def assert_graph(self):
        self.seam.assert_graph()
        if (self.g is not self.seam.g or self.family!=self._family or self.hashes!=self._hashes
            or self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot
            or self.numerical_oracle_available is not False):raise ValueError('Actual Rp pulse source changed')
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed Rp dependency '+name)
    def load_acceptance(self):
        raw,rec=parent.linear.read(NAME),parent.linear.read(RECEIPT)
        pins=dict(rec['input_hashes']);archive=pins.pop(NAME);checker=Path(__file__).stem+'_check.py';checksha=pins.pop(checker)
        if (not rec['all_passed'] or not all(rec[k] for k in GATES) or any(rec[k] for k in OPEN)
            or pins!=self.hashes or raw['input_hashes']!=self.hashes or rec['source_family']!=self.family
            or archive!=sha(NAME) or checksha!=sha(checker) or raw['Rp_native_C5_pulse_source_bindings']!=self.snapshot
            or raw['exact_graph_nodes']!=self.g.nodes[:len(raw['exact_graph_nodes'])]
            or set(rec['checks'])!={'source_inputs','outer_common_sources','compact_quiet','callback_and_consumer'}
            or not all(check['passed'] for check in rec['checks'].values())
            or len(rec['checks']['compact_quiet']['actual_eight_Rp_frame_C5_rows'])!=48
            or len(rec['checks']['callback_and_consumer']['actual_absolute_background_pressure_C3_rows'])!=12
            or len(rec['checks']['compact_quiet']['actual_terminal_integral_and_unique_branch_source_proof']['actual_integral_to_residual_C5_rows'])!=30):
            raise ValueError('Checked actual native C5 pulse source required')
        self.acceptance_loaded=True
    def incoming_source_rows(self):
        self.assert_graph();return dict(self.current_frame)
    def pulse_source_functions(self):
        self.assert_graph()
        if not self.acceptance_loaded:raise ValueError('Source equality must be accepted before consuming current pulse functions')
        return dict(selected=self.selected_pulse,leading_fields=self.leading_fields,pressure=self.pressure,
            coordinates=self.pulse_coordinates,absolute_pressure=self.absolute_pressure,
            native=self.current_frame,numerical_oracle_available=False,global_n1=False)


def run(seam=None):
    began=time.monotonic();owner=CurrentOriginalRpNativeC5PulseSource(seam,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,Rp_native_C5_pulse_source_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ACTUAL_CURRENT_RP_NATIVE_C5_PULSE_SOURCE_DEFINED',len(owner.g.nodes),flush=True);return owner,raw


if __name__=='__main__':raise SystemExit('run() restores checked sources without numerical ancestor replay')
