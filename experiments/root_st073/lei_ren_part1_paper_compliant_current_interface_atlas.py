"""Compose current 14 adjacent / 8 internal traces and compact-time bounds.

The admitted source theorems establish equality before interval bounds.
Physical envelopes enclose the actual traces, not selected point values.
This atlas covers the affected interfaces, not every join in the 33 charts.
"""
import copy
import importlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_angular_support_differences import (
    CurrentAngularSupportDifferences,exact_current_KR_source,HERE,PREFIX,OPEN,
    sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_angular_support_interfaces import EDGES as ANGULAR_EDGES
from lei_ren_part1_paper_compliant_current_pulse_support_interfaces import physical_summary
from lei_ren_part1_paper_compliant_pulse_end_support_interfaces import EDGES as PULSE_EDGES
from lei_ren_part1_paper_compliant_current_physical_interfaces import SEAMS,interface_ledger,supported_beta_edges
from lei_ren_part1_paper_compliant_current_postpulse_interfaces import ordinary_coordinate_proof,statement
from lei_ren_part1_paper_compliant_current_pulse_interfaces import PULSESEAMS
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_interface_atlas.json'
RECEIPT=PREFIX+'current_interface_atlas_check.json'
GATES=('current_22_interface_common_source_composition_certified',
    'current_22_interface_compact_time_physical_spatial4_time1_bounds_available')
STEMS=('current_physical_interfaces','current_postpulse_interfaces','current_pulse_interfaces',
    'current_pulse_support_interfaces','current_angular_support_interfaces','current_angular_support_differences')
PROOF_KEYS=('current_end_flatten_native_source_proof','current_postpulse_source_trace_proof',
    'current_pulse_source_trace_proof','current_pulse_support_source_proof','current_angular_support_source_proof')
BASE_LOG_TAU='-1'
DEFAULT_LOG_TAU=('-3','-1')


def contribution_rows(packet):
    return physical_summary(packet)


def validate_rows(c,rows):
    if len(rows)!=216:raise ValueError('Complete physical spatial4/time1 contribution layout required')
    for row in rows.values():
        gamma=restore_value(c,row['physical_lambda_exponent'])
        if not all(mp.isfinite(v) for v in endpoints(gamma)) or endpoints(gamma)[1]>0:
            raise ValueError('The uniform whole-Z transfer requires finite nonpositive lambda powers')
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None:raise ValueError('Exact zero cannot carry a log bound')
        elif not all(mp.isfinite(v) for v in endpoints(restore_value(c,row['log_absolute_upper']))):
            raise ValueError('Finite positive-time source log bound required')


def merge_trace_rows(c,parts):
    """A common upper for two source-equal traces; overlap proves nothing."""
    if not parts or any(set(v)!=set(parts[0]) for v in parts):raise ValueError('Common trace contribution layout differs')
    result={}
    for key in parts[0]:
        rows=[part[key] for part in parts]
        gamma=restore_value(c,rows[0]['physical_lambda_exponent'])
        if any(endpoints(restore_value(c,row['physical_lambda_exponent']))!=endpoints(gamma) for row in rows):
            raise ValueError('Common physical trace lambda powers differ')
        # A nonzero enclosure can bound an exact zero; require zero on both
        # sides before emitting an exact-zero common contribution.
        values=[endpoints(restore_value(c,row['log_absolute_upper']))[1] for row in rows if not row['exact_zero']]
        result[key]=dict(exact_zero=not values,log_absolute_upper=c.mpf(max(values)) if values else None,
            physical_lambda_exponent=gamma)
    validate_rows(c,result)
    return result


def time_sector_theorem():
    statement('compliant_global_physical_assembly','evaluate','loglambda_bound=logtau/2')
    statement('compliant_global_physical_assembly','log_row',
        'logupper=combined+lambda_exponent*loglambda_bound+c.ln(c.mpf(magnitude))')
    statement('compliant_global_physical_assembly','log_row',
        'bound=c.mpf(max(upper))+c.ln(len(upper)) if upper else None')
    gamma,t,t0,C=s.symbols('gamma log_tau baseline_log_tau source_log_upper',real=True)
    if s.expand((C+gamma*t/2)-(C+gamma*t0/2)-gamma*(t-t0)/2)!=0:
        raise ArithmeticError('Physical source time transfer identity failed')
    return dict(actual_original_log_row_and_lambda_bound_AST_bound=True,
        identity='logUpper(t)=logUpper(t0)+gamma*(t-t0)/2 before outward enclosure',
        source_relation='lambda^2*(1-Z^2)=tau; lambda>=sqrt(tau), gamma<=0',
        conclusion='baseline common upper + sup(gamma*(log_tau-baseline_log_tau)/2) encloses every contribution over a finite compact log(tau) sector',
        all_source_factors_and_ordinary_derivative_lambda_powers_retained=True,
        no_uniform_bound_at_tau_zero_or_global_temporal_flatness_inferred=True,passed=True)


class CurrentInterfaceAtlas:
    @source_precision
    def __init__(self,differences=None,require_checked=True):
        self.differences=differences if differences is not None else CurrentAngularSupportDifferences()
        self.angular=self.differences.support;self.pulse_support=self.angular.pulse_support
        self.pulse_interfaces=self.pulse_support.pulse_interfaces
        self.postpulse=self.pulse_interfaces.postpulse;self.interfaces=self.postpulse.interfaces
        self.owners=(self.interfaces,self.postpulse,self.pulse_interfaces,self.pulse_support,self.angular,self.differences)
        self.physical=self.differences.physical;self.ctx=self.physical.ctx
        self.family=self.physical.family;self.source=self.physical.source;self.datum_sha=self.physical.datum_sha
        self.records={};self.receipts={};self.hashes=dict(self.differences.hashes)
        for stem,owner in zip(STEMS,self.owners):
            if owner.physical is not self.physical or not owner.acceptance_loaded:
                raise ValueError('All current interface layers must retain one checked physical owner')
            name=PREFIX+stem+'.json';receipt_name=PREFIX+stem+'_check.json'
            raw=json.loads((HERE/name).read_bytes());receipt=json.loads((HERE/receipt_name).read_bytes())
            _verify_hashes(raw);_verify_hashes(receipt)
            scoped_gates=importlib.import_module('lei_ren_part1_paper_compliant_'+stem).GATES
            if not receipt['all_passed'] or not all(receipt[key] for key in scoped_gates) or receipt['input_hashes'].get(name)!=sha(name):
                raise ValueError('Each source producer must be consumed by its checked receipt')
            for record in (raw,receipt):
                if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'],record['datum_enclosure_sha256'])!=(self.family,self.source,self.datum_sha):
                    raise ValueError('Different current family/implicit source/analytic pressure datum')
                if any(record[k] for k in OPEN):raise ValueError('An input claims a broader global gate')
                for path,digest in record['input_hashes'].items():
                    if path in self.hashes and self.hashes[path]!=digest:raise ValueError('Conflicting current source hash: '+path)
                    self.hashes[path]=digest
            self.records[stem]=raw;self.receipts[stem]=receipt
            self.hashes[name]=sha(name);self.hashes[receipt_name]=sha(receipt_name)
        for owner,key,stem in zip(self.owners[:5],PROOF_KEYS,STEMS):
            if encode(pack(owner.proof))!=self.records[stem][key]:raise ValueError('Live admitted source theorem differs: '+stem)
        self.graph=self.assert_graph();self.coordinates=ordinary_coordinate_proof(self.physical)
        self.time_proof=time_sector_theorem();self.KR=exact_current_KR_source(self.physical)
        self.local_angular_companion={
            'original_angular_full_moment_transport_theorem':self.differences.transport,
            'original_paper_stress_units':self.differences.units,
            'current_angular_local_difference_source_theorem':self.differences.proof,
            'current_original_stress_and_viscosity_difference_operator_theorem':self.differences.operator_proof}
        angular_record=self.records['current_angular_support_differences']
        if any(encode(pack(value))!=angular_record[key] for key,value in self.local_angular_companion.items()):
            raise ValueError('Live angular local companion theorem differs')
        self.adjacent=copy.deepcopy(self.receipts['current_pulse_interfaces']['current_source_identified_interface_ledger'])
        self.internal=copy.deepcopy(self.receipts['current_angular_support_interfaces']['current_internal_support_ledger'])
        self.ledger_proof=self.assert_ledgers()
        self.baseline={};self.pulse_full_packets={};self.acceptance_loaded=False
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Composed interface receipt exceeds scope')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.acceptance_loaded=True

    def assert_graph(self):
        graph=self.physical.assert_graph()
        links=dict(one_physical_owner=all(owner.physical is self.physical for owner in self.owners),
            same_selected_pulse=self.pulse_support.pulse is self.physical.pulse,
            same_current_angular=self.angular.outer is self.physical.history.outer,
            checked_six_interface_layers=all(owner.acceptance_loaded for owner in self.owners),
            same_family_source_and_pressure_datum=all((owner.family,owner.source,owner.datum_sha)==(self.family,self.source,self.datum_sha) for owner in self.owners))
        if not all(links.values()):raise ValueError('Current common interface owner graph differs')
        return dict(current_complete_physical_graph=graph,current_interface_links=links,passed=True)

    def assert_ledgers(self):
        """Receipt ledger data must match the live owner/coordinate recipes."""
        expected_adjacent=interface_ledger(self.physical)
        for row in expected_adjacent.values():row['current_functional_mixed4_trace_admitted']=True
        if encode(pack(expected_adjacent))!=self.adjacent:
            raise ValueError('Adjacent ledger owner/coordinate/Jacobian/order/history/unit source differs')
        expected_internal=supported_beta_edges(self.physical)
        for row in expected_internal.values():row['current_full_field_interface_trace_admitted']=True
        if encode(pack(expected_internal))!=self.internal:
            raise ValueError('Internal ledger owner/chart/beta/order/exact-edge source differs')
        return dict(adjacent_live_owner_coordinate_Jacobian_history_and_units_bound_count=14,
            internal_live_owner_chart_beta_order_and_edge_bound_count=8,
            exact_live_ledger_recipes_compared_after_source_admission=True,passed=True)

    @source_precision
    def baseline_rows(self,name):
        self.assert_graph()
        self.assert_ledgers()
        if name in self.baseline:return self.baseline[name]
        if name in self.adjacent:
            if name=='end_flatten':point=self.records['current_physical_interfaces']['current_end_flatten_views']['whole_Z']
            elif name in {row[0] for row in PULSESEAMS}:point=self.records['current_pulse_interfaces']['current_whole_Z_pulse_trace_views'][name]
            else:point=self.records['current_postpulse_interfaces']['current_whole_Z_postpulse_trace_views'][name]
            traces=[restore_value(self.ctx,point[key]) for key in ('left_physical_trace','right_physical_trace')]
            for trace in traces:
                if endpoints(trace['Z'])!=(-1,1) or endpoints(trace['requested_log_tau'])!=(-1,-1):
                    raise ValueError('Whole-Z source trace must use the pinned baseline time')
                if endpoints(trace['cos_theta'])!=(-1,1) or endpoints(trace['sin_theta'])!=(-1,1):
                    raise ValueError('All-angle source bound required')
            rows=merge_trace_rows(self.ctx,[contribution_rows(trace) for trace in traces])
        elif name in self.internal:
            edge=self.internal[name]
            if edge['chart']=='outer_angular':
                views=self.records['current_angular_support_interfaces']['current_whole_Z_angular_support_endpoint_views']
                point=next(v for v in views if v['exact_source_edge']==edge['exact_source_edge'])
                rows=restore_value(self.ctx,point['original_physical_spatial4_time1_endpoint_log_bounds'])
            else:
                # Source equality comes from the checked flat-input/FTC theorem.
                # The interval evaluation is solely a bound on the actual full
                # nonzero field at that rational edge; it does not prove equality.
                q=s.Rational(edge['exact_source_edge']);coordinate=self.ctx.mpf(str(q.p))/q.q
                trace=self.physical.evaluate('pulse_end',[-1,1],coordinate,log_tau=BASE_LOG_TAU,theta=None)
                self.pulse_full_packets[name]=trace
                rows=contribution_rows(trace)
            validate_rows(self.ctx,rows)
        else:raise ValueError('Unknown affected current interface: '+name)
        self.baseline[name]=rows
        return rows

    @source_precision
    def bounds(self,name,log_tau=DEFAULT_LOG_TAU,Z=(-1,1)):
        c=self.ctx;sector=c.mpf(log_tau);Z=c.mpf(Z)
        if not all(mp.isfinite(v) for v in endpoints(sector)):raise ValueError('Finite compact log(tau) sector required')
        if endpoints(Z)[0]<-1 or endpoints(Z)[1]>1:raise ValueError('Source Z must lie in [-1,1]')
        rows={};groups={}
        for key,anchor in self.baseline_rows(name).items():
            gamma=restore_value(c,anchor['physical_lambda_exponent'])
            upper=(None if anchor['exact_zero'] else c.mpf(endpoints(restore_value(c,anchor['log_absolute_upper']))[1])+gamma*(sector-c.mpf(BASE_LOG_TAU))/2)
            rows[key]=dict(exact_zero=anchor['exact_zero'],log_absolute_upper=upper,physical_lambda_exponent=gamma)
            groups.setdefault(key.rsplit(':',1)[0],[]).append(upper)
        summed={key:(None if not any(v is not None for v in values) else
            c.mpf(max(endpoints(v)[1] for v in values if v is not None))+c.ln(sum(v is not None for v in values))) for key,values in groups.items()}
        if len(summed)!=144:raise ValueError('Four Cartesian fields spatial4 and time1 aggregation incomplete')
        return dict(interface=name,interface_kind='adjacent' if name in self.adjacent else 'internal_support',
            requested_log_tau=sector,requested_Z=Z,coefficient_bound_Z=c.mpf([-1,1]),all_angles_covered=True,
            physical_contribution_log_bounds=rows,cartesian_velocity_pressure_total_log_bounds=summed,
            actual_full_trace_not_local_difference=True,exact_source_equality_from_checked_function_theorem=True,
            output_kind='upper bounds on actual physical velocity/absolute-pressure derivatives, not resolved point values',
            **dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_common_interface_graph=self.graph,
            current_live_ledger_source_proof=self.ledger_proof,
            retained_current_angular_local_stress_error_companion_theorems=self.local_angular_companion,
            local_companion_theorems_not_actual_completed_tensor=True,
            exact_current_coordinate_source_proof=self.coordinates,compact_time_bound_source_theorem=self.time_proof,
            current_exact_KR_source=self.KR,current_adjacent_source_ledger=self.adjacent,current_internal_source_ledger=self.internal,
            current_source_identified_adjacent_interface_count=14,current_source_identified_internal_support_count=8,
            source_receipt_composition={stem:dict(producer=PREFIX+stem+'.json',receipt=PREFIX+stem+'_check.json',
                producer_sha256=sha(PREFIX+stem+'.json'),receipt_sha256=sha(PREFIX+stem+'_check.json')) for stem in STEMS},
            source_domain='Z=[-1,1]; finite physical traces R>0, |Z|<1, tau>0; Z=+/-1 are source limiting/infinite-space sectors',
            quantitative_domain='all angles, whole source Z, any requested finite compact log(tau) sector; no uniform tau->0 upper',
            uncovered_obligations=['other joins among retained 18 charts','actual completed stress tensor and tensor interfaces',
                'global residual/remainder/cone/lift','core-axis attachment (these 22 outer interfaces have R>0)',
                'resolved coefficient points','prescribed-domain energy','n-dependent temporal recursion'],
            inherited_nonzero_histories_and_analytic_absolute_pressure_not_reset=True,
            interval_overlap_not_used_to_prove_trace_equality=True,input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentInterfaceAtlas(require_checked=False)
    result=field.manifest();result['current_22_compact_time_interface_bounds']={}
    for name in (*field.adjacent,*field.internal):
        result['current_22_compact_time_interface_bounds'][name]=field.bounds(name)
    result['current_four_pulse_full_endpoint_baseline_log_bounds']={name:field.baseline_rows(name) for name,row in field.internal.items() if row['chart']=='pulse_end'}
    result['input_hashes']=field.hashes
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current 22-interface source atlas and actual compact-time bounds generated',flush=True)
    return result


if __name__=='__main__':run()
