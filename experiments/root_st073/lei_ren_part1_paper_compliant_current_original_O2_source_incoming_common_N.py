"""Same-N actual upstream correction covers and genuine continuous O2 integrals.

The checked upstream [-1,1] N1024 history covers remain factored. Fresh N1024
true-radius phase unions drive the original O2 value/ordinary-Z density graph.
No N7 contribution, source-cover midpoint, or zero downstream inlet is used.
"""
import gzip
import hashlib
import json
from pathlib import Path
from fractions import Fraction
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_continuous_Z_transport as current
import lei_ren_part1_paper_compliant_current_native_middle_O2_inlet_C1_histories as middle
import lei_ren_part1_paper_compliant_current_native_Rc_function_transport as functions

base=current.base;ep=current.ep;previous=current.current
HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
NAME=PREFIX+'current_original_O2_source_incoming_common_N.json'
RECEIPT=PREFIX+'current_original_O2_source_incoming_common_N_check.json'
GATE='original_O2_same_N1024_actual_upstream_correction_and_genuine_continuous_C1_integrals_connected'
KEYS=current.KEYS;N=1024


def interval(c,value):
    return c.mpf(value) if hasattr(value,'_mpi_') else previous.interval(c,value)


def attach_receipt(hashes,module,family):
    receipt=json.loads((HERE/module.RECEIPT).read_bytes())
    report=json.loads((HERE/module.NAME).read_bytes())
    if not receipt.get('all_passed') or not receipt.get(module.GATE) or report['source_family']!=family:
        raise ValueError('Accepted same-family source required: '+module.RECEIPT)
    for name,digest in {**receipt['input_hashes'],module.NAME:sha(module.NAME),module.RECEIPT:sha(module.RECEIPT)}.items():
        if sha(name)!=digest:raise ValueError('Same-N source changed: '+name)
        if name in hashes and hashes[name]!=digest:raise ValueError('Source closures disagree: '+name)
        hashes[name]=digest
    return report,receipt


def source_function_prefix_binding(report):
    cells=report['exact_original_cells']
    indices=[i for i,row in enumerate(cells) if row['chart']=='O2_slope']
    if not indices or indices[0]!=13:raise ValueError('Complete 13-cell original pre-O2 route required')
    index=indices[0];first=cells[index]
    expected=['initial_flat_collar','active_first_bridge','second_bridge','bridge_macro',
        'switch_first','switch_second','switch_power','reshape','inner_reference','axial_restore',
        'restore_buffer','actual_patch','Rh_reference']
    if [row['label'] for row in cells[:index]]!=expected:
        raise ValueError('Original upstream chart sequence changed')
    if first['incoming']!=cells[index-1]['outgoing']:raise ValueError('O2 incoming function nodes lost upstream memory')
    nodes=report['function_graph_nodes'];zero=next(i for i,n in enumerate(nodes)
        if n['operation']=='exact_rational' and n['numerator']==0)
    if any(pair!=dict(value=zero,Z=zero) for pair in cells[0]['incoming'].values()):
        raise ValueError('Actual original initial correction must be source-defined zero')
    if not cells[0]['source_flat_exact_zero'] or any(row['source_flat_exact_zero'] for row in cells[1:index]):
        raise ValueError('Active upstream charts cannot be declared quiet')
    return dict(exact_function_source_manifest=functions.NAME,source_manifest_sha256=sha(functions.NAME),
        common_N_lower=report['common_N_lower'],original_true_inlet='r_minus at sc/2',
        source_defined_initial_correction_zero=True,original_background_and_P0_not_zeroed=True,
        ordered_pre_O2_function_cells=expected,O2_first_cell_index=index,
        exact_O2_incoming_correction_function_roots=first['incoming'],
        exact_O2_source_cell_contribution_function_roots=first['contributions'],
        O2_incoming_is_previous_Rh_reference_outgoing=True,
        graph_roots_define_functions_not_numeric_values=True)


def restore_common_source(record,coordinates):
    """Restore a complete signed function cover, never a selected point."""
    c=coordinates.ctx;scale=record['formal_positive_scale']
    powers=tuple(scale['source_exponents'])+(scale['radius_power'],)
    if any(powers[i] for i in (0,2,3,4)):
        raise ValueError('Expected accepted common Pstar-squared source coordinate')
    value=base.prior.ScaledEnclosure(base.prior.FormalScale(coordinates.bases,powers,
        interval(c,scale['additional_log_interval'])),interval(c,record['coefficient_interval']),coordinates.ledger)
    if value.zero!=record['exact_zero'] or record['point_value_selected'] or not record['encloses_original_source_function']:
        raise ValueError('Actual incoming cover cannot become a point or fallback zero')
    return value


def absolute_log_upper(value):
    magnitude=max(abs(x) for x in ep(value.coefficient))
    return None if not magnitude else value.scale.evaluate()+value.ctx.ln(value.ctx.mpf(magnitude))


class OriginalO2SourceIncomingCommonN:
    def __init__(self):
        self.parent=current.OriginalO2ContinuousZTransport();self.c=self.parent.c;self.family=self.parent.family
        self.hashes=dict(self.parent.hashes)
        attach_receipt(self.hashes,current,self.family)
        incoming,receipt=attach_receipt(self.hashes,middle,self.family)
        exact,unused=attach_receipt(self.hashes,functions,self.family)
        self.exact_binding=source_function_prefix_binding(exact)
        if incoming['candidate_N']!=N or self.exact_binding['common_N_lower']>N:
            raise ValueError('Same checked candidate N1024 required')
        self.incoming=incoming['actual_inlet_to_O2_inlet_C1_history_records']['whole_Z']
        if self.incoming['candidate_N']!=N or self.incoming['source_family']!=self.family:
            raise ValueError('Same original upstream source/N required')
        if not self.incoming['no_original_interval_skipped_between_true_inlet_and_O2_inlet']:
            raise ValueError('Complete actual upstream route required')
        self.ordered=self.parent.parent.parent.parent;self.original=self.ordered.owner
        self.pressure=previous.five.pressure
        background=self.incoming['actual_O2_inlet_original_background_and_separate_P0_Z']
        old_bases=background['common_directed_coordinate_theorem']['common_log_bases']
        self.coordinates=middle.history.CommonSourceCoordinates(self.c,interval(self.c,old_bases[1]),self.family)
        if any(ep(interval(self.c,old_bases[i]))!=(0,0) for i in (0,2,3,4)):
            raise ValueError('Original common source basis changed')
        if middle.prior.ScaledEnclosure is not base.prior.ScaledEnclosure or (
            {k:Fraction(str(v)) for k,v in middle.RATES.items()}!=
            {k:Fraction(str(v)) for k,v in previous.five.RATES.items()}):
            raise ValueError('Same original signed arithmetic and own-rate units required')
        if background['original_units']!=previous.five.UNITS or not background['P0_not_merged_into_pressure_history']:
            raise ValueError('Same normalized five units and separate P0 required')
        with mp.workdps(self.c.dps+40):
            actual_P2=(self.c.exp(40)+11)*2;old=self.coordinates.logP_squared
            if max(ep(actual_P2)[0],ep(old)[0])>min(ep(actual_P2)[1],ep(old)[1]):
                raise ValueError('Upstream and original O2 Pstar definitions disagree')
            self.incoming_values={k:restore_common_source(v,self.coordinates) for k,v in
                self.incoming['actual_original_inlet_to_O2_inlet_correction_C0'].items()}
            self.incoming_Z={k:restore_common_source(v,self.coordinates) for k,v in
                self.incoming['actual_original_inlet_to_O2_inlet_correction_Z'].items()}
            self.P0=restore_common_source(background['original_separate_P0_over_Pstar_squared'],self.coordinates)
            self.P0_Z=restore_common_source(background['original_separate_P0_Z_over_Pstar_squared'],self.coordinates)
        if set(self.incoming_values)!=set(KEYS) or set(self.incoming_Z)!=set(KEYS):
            raise ValueError('All five actual correction value/Z functions required')
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def phase_level(self,count,N):
        if N!=self.incoming['candidate_N']:raise ValueError('Do not mix diagnostic N7 with N1024 upstream histories')
        level,saved=self.ordered.levels[count]
        origin=self.original.radius.evaluate(y=0,N=N);c=self.c
        rows=[]
        with mp.workdps(c.dps+40):
            for i,row in enumerate(saved['whole_source_cells']):
                covers=self.pressure.phase_boxes(c,origin,N,c.mpf(i)/count,c.mpf(i+1)/count)
                rows.append(dict(row,true_common_N_phase_boxes=[base.encoded(v) for v in covers]))
        return level,dict(saved,explicit_candidate_N=N,original_true_radius_phase_at_y0=origin,
            whole_source_cells=rows),origin

    def compose(self,record):
        if record['source_family']!=self.family or record['explicit_candidate_N']!=N:
            raise ValueError('Same-N source functional record required')
        current.require_identity(record,self.family,self.family['datum_enclosure_sha256'])
        c=self.c;coordinates=self.coordinates;sourceZ=interval(c,self.incoming['Z_box'])
        zl,zh=(base.point.pressure.exact_Z(v) for v in record['exact_Z_range'])
        with mp.workdps(c.dps+40):
            target=c.mpf((ep(c.mpf(int(zl.p))/int(zl.q))[0],ep(c.mpf(int(zh.p))/int(zh.q))[1]))
            if not ep(sourceZ)[0]<=ep(target)[0]<=ep(target)[1]<=ep(sourceZ)[1]:
                raise ValueError('Actual incoming cover must include entire target Z tile')
            original={k:coordinates.scalar(interval(c,record['source_defined_original_histories_at_y1'][k])) for k in KEYS}
            originalZ={k:coordinates.scalar(interval(c,record['source_defined_original_history_Z_at_y1'][k])) for k in KEYS}
            modulation={k:coordinates.scalar(interval(c,record['five_original_C0_integral_contributions'][k])) for k in KEYS}
            modulationZ={k:coordinates.scalar(interval(c,record['five_genuine_ordinary_Z_integral_contributions'][k])) for k in KEYS}
            correction={};correctionZ={};own={};ownZ={};budgets={}
            for key,rate in previous.five.RATES.items():
                decay=coordinates.decay(1,rate)
                correction[key]=decay*self.incoming_values[key]+modulation[key]
                correctionZ[key]=decay*self.incoming_Z[key]+modulationZ[key]
                own[key]=original[key]+correction[key];ownZ[key]=originalZ[key]+correctionZ[key]
                budgets[key]=dict(incoming_C0_log_absolute_upper=absolute_log_upper(self.incoming_values[key]),
                    incoming_Z_log_absolute_upper=absolute_log_upper(self.incoming_Z[key]),
                    O2_C0_log_absolute_upper=absolute_log_upper(modulation[key]),
                    O2_Z_log_absolute_upper=absolute_log_upper(modulationZ[key]),
                    logarithmic_envelope_comparison_not_physical_point_values=True)
            final=dict(source_family=self.family,candidate_N=N,exact_Z_range=record['exact_Z_range'],
                true_unit_O2_y_window=['0','1'],normalized_own_units=previous.five.UNITS,
                original_P0_datum_sha256=self.family['datum_enclosure_sha256'],
                actual_incoming_source_binding=dict(manifest=middle.NAME,manifest_sha256=sha(middle.NAME),
                    receipt=middle.RECEIPT,receipt_sha256=sha(middle.RECEIPT),source_record_key='whole_Z',
                    actual_upstream_Z_cover=self.incoming['Z_box'],source_family=self.family,candidate_N=N,
                    common_coordinates=coordinates.record(),same_exact_incoming_function_graph=self.exact_binding,
                    input_is_correction_not_original_plus_correction=True,
                    wider_whole_Z_cover_used_on_subtile_not_a_refined_function_value=True),
                actual_source_incoming_C0={k:v.record() for k,v in self.incoming_values.items()},
                actual_source_incoming_Z={k:v.record() for k,v in self.incoming_Z.items()},
                genuine_O2_window_C0={k:v.record() for k,v in modulation.items()},
                genuine_O2_window_Z={k:v.record() for k,v in modulationZ.items()},
                original_background_history_C0={k:v.record() for k,v in original.items()},
                original_background_history_Z={k:v.record() for k,v in originalZ.items()},
                propagated_actual_correction_C0={k:v.record() for k,v in correction.items()},
                propagated_actual_correction_Z={k:v.record() for k,v in correctionZ.items()},
                actual_original_plus_propagated_correction_C0={k:v.record() for k,v in own.items()},
                actual_original_plus_propagated_correction_Z={k:v.record() for k,v in ownZ.items()},
                separate_P0=self.P0.record(),separate_P0_Z=self.P0_Z.record(),
                actual_absolute_pressure=(self.P0+own['p']).record(),
                actual_absolute_pressure_Z=(self.P0_Z+ownZ['p']).record(),
                conservative_source_envelope_budget=budgets,
                no_incoming_zero_fallback_or_N7_integral_reused=True,
                common_factored_context_ledger_retained_no_enormous_scalar_cast=True,
                original_background_added_once_after_correction_transport=True,
                rate_zero_pressure_correction_retains_entire_actual_upstream_memory=True,
                full_Rc_route_or_functional_terminal_identity_solved=False,
                **dict.fromkeys(current.FLAGS,False))
        return final

    def integrate(self,count,*,Z_lower,Z_upper,N=1024,bits=24):
        saved=self.ordered.levels[count]
        level,new_phase,origin=self.phase_level(count,N)
        # The original y/J/mass/source cache is N-independent. Only a temporary
        # in-memory phase cover is replaced; no accepted producer file changes.
        self.ordered.levels[count]=(level,new_phase)
        try:report=self.parent.integrate(count,Z_lower=Z_lower,Z_upper=Z_upper,N=N,bits=bits)
        finally:self.ordered.levels[count]=saved
        report['fresh_same_N_true_radius_phase_origin']=origin
        report['phase_cover_regenerated_for_actual_upstream_N']=True
        report['accepted_N_independent_y_source_cache_phase_candidate']=saved[1]['explicit_candidate_N']
        report['actual_source_incoming_and_own_history_transport']=self.compose(report)
        return report


def save_integral(report,tag):
    count=report['ordered_source_cells'];rows=report.pop('whole_original_density_and_density_Z_source_records');archives=[]
    for left in range(0,count,512):
        right=min(count,left+512);name=PREFIX+'current_original_O2_source_incoming_common_N_'+tag+'_'+str(count)+'_cells_'+str(left)+'_'+str(right)+'.json.gz'
        chunk=dict(source_family=report['source_family'],candidate_N=N,exact_Z_range=report['exact_Z_range'],
            source_cell_count=count,source_cell_index_range=[left,right],records=rows[left:right],
            fresh_same_N_true_radius_phase_origin=report['fresh_same_N_true_radius_phase_origin'])
        raw=json.dumps(base.encoded(chunk),indent=2).encode()+b'\n';compressed=gzip.compress(raw,compresslevel=9,mtime=0)
        if len(compressed)>=100*1024*1024:raise ArithmeticError('Split complete evidence into smaller lossless chunks')
        (HERE/name).write_bytes(compressed)
        archives.append(dict(filename=name,compressed_bytes=len(compressed),uncompressed_bytes=len(raw),
            lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),source_cell_index_range=[left,right]))
    report['whole_original_density_and_density_Z_source_archives']=archives
    report['every_complete_source_cell_saved_losslessly']=True
    return report


def run():
    begin=time.monotonic();owner=OriginalO2SourceIncomingCommonN();reports=[]
    for tag,zl,zh in (('positive','.36','.38'),('negative','-.38','-.36')):
        for count in (256,2048):
            report=owner.integrate(count,Z_lower=zl,Z_upper=zh,N=N)
            reports.append(save_integral(report,tag))
            print('Actual N1024 source incoming and genuine O2 C1 integral archived:',tag,count,flush=True)
    result=dict(**{GATE:True},source_family=owner.family,candidate_N=N,
        original_complete_same_N_source_incoming_O2_integral_refinements=reports,
        exact_upstream_source_function_incoming_binding=owner.exact_binding,
        genuine_N1024_O2_phase_and_integrals_not_N7_relabelling=True,
        accepted_upstream_whole_Z_function_covers_restored_as_factored_covers=True,
        upstream_cover_widths_not_improved_by_downstream_O2_integration=True,
        all_chart_terminal_matching_or_global_N_or_true_recursion_admitted=False,
        **dict.fromkeys(current.FLAGS,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-begin,
        scope='Same candidate N1024 original upstream sc/2-to-O2-slope0 correction function covers and exact 13-cell graph roots feed fresh continuous strict-sign O2 y[0,1] genuine value/Z integrals. Factored original-plus-correction and separate P0/pressure memory. Upstream bounds remain conservative; no Z0 atlas/full Rc matching/controls/global N/stress/recursion/full NS.')
    (HERE/NAME).write_text(json.dumps(base.encoded(result),indent=2)+'\n',encoding='utf8');return result


if __name__=='__main__':run()
