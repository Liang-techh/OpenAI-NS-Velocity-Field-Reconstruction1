"""Live-source numerical containment for the true active/compact/quiet C1 graph.

Signed functions are retained graph integrals. Absolute source bounds are
nonzero symmetric enclosures, with the quiet signs and exact Z1 zero kept.
The live all-chart field owner is required: old magnitude archives alone
do not retain the function references needed for this attachment.
"""
import copy
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_n1_lower_numeric_bridge as lower
import lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions as true
import lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_C1 as caps
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

HERE,PREFIX,sha=lower.HERE,lower.PREFIX,lower.sha
op,ep=lower.op,lower.ep
NAME=PREFIX+'current_original_n1_preRp_numeric_bridge.json.gz'
RECEIPT=PREFIX+'current_original_n1_preRp_numeric_bridge_check.json'
GATE='current_original_true_preRp_pressure_C1_numeric_bridge_installed'
OPEN=lower.OPEN


@dataclass(frozen=True,eq=False)
class OriginalTruePreRpNumericPacket:
    scope:str='original rminus to Rp pressure value/Z1 enclosures'


class CurrentOriginalN1PreRpNumericBridge:
    @source_precision
    def __init__(self,accepted_lower,fields,numerical,require_checked=True):
        if type(accepted_lower) is not lower.CurrentOriginalN1LowerNumericBridge or not accepted_lower.acceptance_loaded:
            raise ValueError('Checked original lower numerical bridge required')
        if type(fields) is not true.active.CurrentOriginalN1AllChartPhysicalPressureFields or not fields.acceptance_loaded:
            raise ValueError('Live checked original all-chart physical field required')
        if type(numerical) is not caps.CurrentOriginalN1PreRpPressureC1 or not numerical.acceptance_loaded:
            raise ValueError('Checked same-source pre-Rp numerical C1 bounds required')
        self.lower,self.fields,self.numeric=accepted_lower,fields,numerical
        self.preRp=accepted_lower.inner.full.preRp;self.g=self.preRp.g;self.c=numerical.c
        self.family=copy.deepcopy(self.preRp.family)
        self.functions,self.pieces=self.preRp.functions,self.preRp.pieces
        self.quiet=true.quiet.CurrentOriginalN1QuietPressurePrimitives()
        self.bindings=true.encode(dict(functions=self.functions,pieces=self.pieces))
        self._bindings=copy.deepcopy(self.bindings)
        self.source_identity=self.assert_source_identity()
        self._identity=op.snapshot(self.source_identity)
        self._value_sources=op.snapshot((numerical.parents['old']['source_data']['cells'],
            numerical.parents['old']['original_physical_units']))
        self._source_parents=op.snapshot((numerical.parents['sources'],numerical.parents['endpoints']))
        self._owners=(self.lower,self.fields,self.numeric,self.preRp,self.quiet,self.g,self.c)
        self.hashes={}
        for owner in (accepted_lower,fields,numerical,self.preRp):
            for name,digest in owner.hashes.items():
                if name in self.hashes and self.hashes[name]!=digest:
                    raise ValueError('Conflicting original pre-Rp source '+name)
                self.hashes[name]=digest
        for module in (lower,true,caps):
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self.hashes[name]=sha(name)
        for name in (Path(__file__).name,Path(__file__).stem+'_check.py'):self.hashes[name]=sha(name)
        for name,digest in self.hashes.items():
            if sha(name)!=digest:raise ValueError('Changed pre-Rp numerical source '+name)
        self._hashes=dict(self.hashes);self._packets={};self.acceptance_loaded=False
        if require_checked:
            receipt=lower.full.parent.linear.read(RECEIPT);pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
                or pins!=self.hashes or archive!=sha(NAME) or receipt.get('source_family')!=self.family
                or receipt.get('source_identity')!=op.snapshot(self.source_identity)):
                raise ValueError('Checked true pre-Rp numerical bridge required')
            self.acceptance_loaded=True

    def assert_source_identity(self):
        self.lower.assert_graph();self.fields.assert_graph();self.preRp.assert_graph();self.numeric.assert_graph()
        if self.fields.family!=self.family or self.numeric.family!=self.family:
            raise ValueError('Pre-Rp true functions and numerical source family differ')
        selected=self.fields.g.nodes[self.fields.N.node]
        repair=self.preRp.operator.reports[true.linear.quiet.band.NAME]['actual_selected_repair_integer']
        if selected['operation']!='exact_positive_integer_expression' or any(selected[k]!=repair[k] for k in ('type','J','definition')):
            raise ValueError('Live active field and actual compact repair N differ')
        before=len(self.g.nodes)
        replay,bindings=true.active_functions(self.preRp.operator,self.fields)
        if (true.encode(replay)!=true.encode(self.pieces['active']) or len(self.g.nodes)!=before
            or bindings['exact_signed_integrals']!=self.preRp.raw['defining_bindings']['active']['exact_signed_integrals']
            or bindings['finite_endpoints']!=self.preRp.raw['defining_bindings']['active']['finite_endpoints']):
            raise ValueError('Live active function and endpoint bodies differ from retained true functions')
        importer=true.GraphImport(self.g,self.fields.g)
        scale=importer(self.fields.functions['first_micro']['scale'].value)
        delta=importer(self.fields.delta_node)
        compact,compact_bindings=true.compact_functions(self.preRp.operator,scale,delta)
        quiet_replay,quiet_bindings=true.quiet_functions(self.preRp.operator,scale,self.quiet)
        if (true.encode(compact)!=true.encode(self.pieces['compact'])
            or true.encode(quiet_replay)!=true.encode(self.pieces['quiet']) or len(self.g.nodes)!=before
            or compact_bindings!=self.preRp.raw['defining_bindings']['compact']
            or quiet_bindings['constant_Z_identity']!=self.preRp.raw['defining_bindings']['quiet']['constant_Z_identity']):
            raise ValueError('Actual compact C3 or quiet native function recipes differ')
        data=self.numeric.data;old=self.numeric.parents['old']['source_data'];sources=self.numeric.parents['sources']
        control=self.numeric.parents['old_control']['actual_quantitative_C3_bounds']
        if (self.numeric.parents['old_band']['actual_selected_repair_integer']!=repair
            or control['actual_same_repair_integer']!=repair):
            raise ValueError('Numerical primitive/control and true compact selected N differ')
        logN=caps.interval(self.c,control['actual_log_N'])
        if (len(self.fields.charts)!=17 or len(old['actual_chart_views'])!=68
            or len(data['upgraded_chart_views'])!=68 or len(data['native_primitive_bounds'])!=228
            or not data['exact_shared_selected_N_checked'] or not data['original_functions_and_selected_N_unchanged']
            or not data['original_target_and_control_C3_update']['same_selected_N']
            or not old['geometry']['original_Rc_ge_one']
            or not old['scope']['original_shared_Rc_cancelled_before_boxing']):
            raise ValueError('Complete original pre-Rp source, history and scale coverage required')
        expected={(tuple(ends),chart) for ends in caps.CELLS for chart in self.fields.charts}
        if {(tuple(v['ends']),v['chart']) for v in data['upgraded_chart_views']}!=expected:
            raise ValueError('Missing or different original axial/chart source coverage')
        for view in data['upgraded_chart_views']:
            if not all(view[k] for k in ('original_leading_history_retained','all_preceding_charts_in_mean',
                'common_Pstar_applied_once','original_B_function_unchanged')) or (
                view['q_recovery_source']!='physical._physical_q_bounds'
                or caps.interval(self.c,view['actual_selected_N_log'])._mpi_!=logN._mpi_):
                raise ValueError('Actual pre-Rp mean history was reset or changed')
        for view in data['native_primitive_bounds']:
            primitive=true.active.physical._z1_encode(self.fields.modified.target.primitives[view['chart']])
            if (view['original_function_handles']!=primitive or not view['predecessor_never_reset']
                or view['rate_one_kernel_upper']!=1
                or not view['actual_whole_chart_rate_zero_width_used_as_subcell_upper']):
                raise ValueError('Original finite-N primitive chain changed')
        quiet=sources['quiet_source'];identity=quiet['source_identity']
        qbind=self.preRp.raw['defining_bindings']['quiet'];qzero=qbind['constant_Z_identity']
        if (not quiet['source_scaled_Z1_exact_zero'] or not quiet['signed_coefficients_not_magnitude_caps']
            or not identity['Pstar_equals_exp_logP_cancelled_before_evaluation']
            or not identity['actual_compact_repaired_quiet_source_consumed']
            or not identity['physical_Rc_is_absolute_not_original_phase_offset']
            or not qzero['exact_Z_independent_kappa'] or not qzero['exact_q_common_at_Rp_is_minus_kappa']
            or any(pair.Z!=self.g.zero for pair in self.pieces['quiet'].values())):
            raise ValueError('Actual quiet native kernel, signed source or exact Z1 zero differs')
        for cell in data['cells']:
            for key,row in cell['fields'].items():
                if not all(row[k] for k in ('quiet_Z_exact_zero','counted_once','original_pressure_units_fixed',
                    'no_signed_zero_assigned_to_bounded_terms')):
                    raise ValueError('Pressure derivative units or nonzero errors changed')
                if any(value is not None and ep(caps.interval(self.c,value))[1]>=-1000
                    for value in row['normalized_log_caps'].values()):
                    raise ValueError('Pre-Rp derivative bounds exceed original enclosure cover')
        return dict(passed=True,live_actual_field_owner_used_not_archive_family_only=True,
            exact_replayed_active_signed_C1_integrals_and_finite_endpoints=True,
            exact_replayed_compact_and_quiet_physical_C1_function_bodies=True,
            active_graph_and_reference_bindings=true.encode(replay),
            original_selected_N_and_compact_repair_integer=repair,
            all_17_charts_68_views_228_native_cells_retained=True,
            all_228_C3_primitive_handles_identical_to_actual_live_source_functions=True,
            all_68_derivative_views_use_the_same_actual_N_log_and_q_recovery_recipe=True,
            actual_full_mean_history_and_Pstar_once_preserved=True,
            quiet_native_source_identity=identity,quiet_true_function_bindings=qbind,
            same_original_radius_units_shared_Rc_cancelled_before_boxing=True,
            checked_C1_input_and_endpoint_receipt_chains_are_authoritative=True,
            magnitude_bounds_attached_as_symmetric_enclosures_not_signed_centers=True,
            original_physical_units=copy.deepcopy(self.numeric.parents['old']['original_physical_units']),
            complete_pulse_pressure_and_reserved_repair_not_admitted=True)

    def assert_graph(self):
        if ((self.lower,self.fields,self.numeric,self.preRp,self.quiet,self.g,self.c)!=self._owners
            or self.hashes!=self._hashes or op.snapshot(self.source_identity)!=self._identity
            or true.encode(dict(functions=self.functions,pieces=self.pieces))!=self._bindings
            or op.snapshot((self.numeric.parents['old']['source_data']['cells'],
                self.numeric.parents['old']['original_physical_units']))!=self._value_sources
            or op.snapshot((self.numeric.parents['sources'],self.numeric.parents['endpoints']))!=self._source_parents):
            raise ValueError('True pre-Rp numerical attachment changed')
        self.lower.assert_graph();self.fields.assert_graph();self.preRp.assert_graph();self.numeric.assert_graph()
        self.quiet.assert_source()

    @source_precision
    def evaluate(self,Z):
        self.assert_graph();c=self.c;z=c.mpf(Z)
        if ep(z)[0]<-1 or ep(z)[1]>1:raise ValueError('Original Z domain [-1,1] required')
        active=[];beta=c.exp(-1000);negative=c.mpf([-1,1])*beta
        values={key:[] for key in ('m13','m15')};piece_rows=[]
        old=self.numeric.parents['old']['source_data']
        for index,(cell,oldcell) in enumerate(zip(self.numeric.data['cells'],old['cells'])):
            left,right=map(c.mpf,cell['ends'])
            if ep(z)[1]<ep(left)[0] or ep(z)[0]>ep(right)[1]:continue
            if tuple(cell['ends'])!=tuple(oldcell['ends']):raise ValueError('Original value/Z1 axial cells differ')
            active.append(index);pieces={}
            for key in values:
                v=oldcell['fields'][key]
                if not v['all_pieces_counted_once'] or not v['no_signed_zero_assigned']:
                    raise ValueError('Pre-Rp value source omitted a nonzero piece')
                if any(ep(caps.interval(c,q))[1]>=-1000 for q in v['normalized_log_caps'].values()):
                    raise ValueError('Pre-Rp value bound exceeds original enclosure cover')
                quiet=c.mpf([0,ep(beta)[1]]) if key=='m13' else c.mpf([-ep(beta)[1],0])
                pieces[key]=dict(active=[negative,negative],compact=[negative,negative],quiet=[quiet,c.mpf(0)])
                values[key].append([2*negative+quiet,2*negative])
            piece_rows.append(dict(exact_Z_cell=cell['ends'],piece_enclosures=pieces))
        if not active:raise ValueError('Missing original real-Z cell coverage')
        hull=lambda rows:c.mpf([min(ep(v)[0] for v in rows),max(ep(v)[1] for v in rows)])
        result=dict(source_family=self.family,Z=z,true_function_bindings=self.bindings,
            radial_scope='original rminus to Rp',active_cells=active,
            normalized_true_preRp_C1_enclosures={k:[hull([row[j] for row in rows]) for j in (0,1)] for k,rows in values.items()},
            original_positive_cover_per_nonzero_piece=beta,cell_piece_enclosures=piece_rows,
            normalization=self.source_identity['original_physical_units'],
            ordinary_Z0_Z1_enclosures_not_defining_function_values=True,
            quiet_m13_positive_m15_negative_and_Z1_exact_zero=True,
            inner_lower_and_pulse_pressure_not_added=True,full_pressure_derivative_accuracy_not_claimed=True,
            no_original_radius_or_hyper_exponential_cap_materialized=True,
            **{GATE:self.acceptance_loaded},**dict.fromkeys(OPEN,False))
        packet=OriginalTruePreRpNumericPacket()
        self._packets[id(packet)]=(packet,result,op.snapshot(result));return packet

    @source_precision
    def report(self,packet):
        self.assert_graph();entry=self._packets.get(id(packet))
        if type(packet) is not OriginalTruePreRpNumericPacket or entry is None or entry[0] is not packet or op.snapshot(entry[1])!=entry[2]:
            raise ValueError('Live true pre-Rp C1 enclosure packet required')
        return {**copy.deepcopy(entry[1]),GATE:self.acceptance_loaded}


@source_precision
def run(accepted_lower,fields,numerical):
    began=time.monotonic();owner=CurrentOriginalN1PreRpNumericBridge(accepted_lower,fields,numerical,False)
    packets={name:owner.evaluate(Z) for name,Z in (('point','.371'),('boundary','0'),('whole_Z',['-1','1']))}
    raw=dict(source_family=owner.family,source_identity=owner.source_identity,true_function_bindings=owner.bindings,
        samples={n:owner.report(p) for n,p in packets.items()},input_hashes=owner.hashes,
        **{GATE:False},**dict.fromkeys(OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(op.snapshot(raw),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,raw,packets
