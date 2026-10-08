"""Actual N1024 full-predicate density covers on five exact O2 intervals.

Restrict the measure of accepted source-density covers at original seams.
Each piece keeps its whole-bin source enclosure and receives a fresh true
Duhamel mass. Previously integrated whole-bin contributions are not copied.
"""
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_O2_full_predicate_Rc_continuation as preceding

native=preceding.native;history=preceding.history;packets=preceding.downstream.packets
HERE,PREFIX,sha=preceding.HERE,preceding.PREFIX,preceding.sha
ep=preceding.ep;encode=preceding.encode;KEYS=preceding.KEYS;N=1024
five=preceding.driver.integrals.five;iv=preceding.driver.integrals.interval
NAME=PREFIX+'current_original_O2_five_slope_partitions.json.gz'
RECEIPT=PREFIX+'current_original_O2_five_slope_partitions_check.json'
GATE='actual_N1024_full_predicate_O2_five_original_slope_interval_C0_Z_operators_installed'
EDGES=tuple(map(Fraction,('0','3/25','13/100','7/50','3/20','1')))
LABELS=('slope_0_12','slope_12_13','slope_13_14','slope_14_15','slope_15_1')


def rational(c,value):
    value=Fraction(value);return c.mpf(value.numerator)/value.denominator


def records(group):return {key:value.record() for key,value in group.items()}


def exact_pieces(left,right):
    left,right=Fraction(left),Fraction(right)
    if not 0<=left<right<=1:raise ValueError('Exact positive O2 interval inside[0,1] required')
    result=[]
    for index in range(64):
        lo=max(left,Fraction(index,64));hi=min(right,Fraction(index+1,64))
        if lo<hi:result.append((index,lo,hi))
    if not result or result[0][1]!=left or result[-1][2]!=right or any(a[2]!=b[1] for a,b in zip(result,result[1:])):
        raise ArithmeticError('Exact O2 piece partition must have no gap or overlap')
    return result


class OriginalO2FiveSlopePartitions:
    def __init__(self,bridge=None):
        self.source_owner=preceding.OriginalO2FullPredicateRcContinuation(bridge)
        owner=self.source_owner;self.ctx=c=owner.ctx;self.coordinates=owner.coordinates;self.family=owner.family
        self.hashes=dict(owner.hashes);self.driver=owner.saved
        self.direct=owner.manifest['actual_fresh_original_N1024_direct_integrals']
        if self.direct['explicit_candidate_N']!=N or self.direct['source_family']!=self.family or self.direct['ordered_source_cells']!=64:
            raise ValueError('Accepted actual original full-predicate 64-cell N1024 density archive required')
        if self.direct['exact_y_window']!=['0','1'] or self.direct['exact_Z_range']!=['-1','1']:
            raise ValueError('Whole original y and Z domain required')
        if self.direct['normalized_own_units']!=five.UNITS or self.direct['own_rates']!=five.RATES:
            raise ValueError('Original own-rate density units required')
        common=self.driver['actual_upstream_source_binding']['common_coordinates']
        if common['source_family']!=self.family:
            raise ValueError('Same actual thirteen-cell upstream family required')
        for actual,record in zip(self.coordinates.bases,common['common_log_bases'],strict=True):
            if ep(actual)!=ep(packets.interval(c,record)):raise ValueError('Actual upstream source factor changed')
        restore=lambda record:preceding.driver.upstream.restore_common_source(record,self.coordinates)
        self.incoming=dict(values={key:restore(row) for key,row in self.driver['actual_upstream_C0'].items()},
            Z_derivatives={key:restore(row) for key,row in self.driver['actual_upstream_Z'].items()})
        if any(set(group)!=set(KEYS) for group in self.incoming.values()):
            raise ValueError('All five actual upstream correction C0/Z rows required')
        inlet=self.direct['source_defined_original_inlet']
        convert=lambda row:owner.native_to_common(preceding.driver.integrals.restore_value(owner.atlas,row))
        self.background=dict(values={key:convert(row) for key,row in inlet['native_original_inlet_sources'].items()},
            Z_derivatives={key:convert(row) for key,row in inlet['native_original_inlet_Z_sources'].items()})
        self.Z_unit=convert(self.direct['changed_ordinary_Z_unit'])
        self.P0=convert(self.driver['separate_analytic_P0']);self.P0_Z=convert(self.driver['separate_analytic_P0_Z'])
        self.source_cells=self.direct['whole_source_density_phase_mass_records']
        if len(self.source_cells)!=64:raise ValueError('All64 accepted source-density cells required')
        self.densities=[]
        for index,cell in enumerate(self.source_cells):
            if tuple(map(Fraction,cell['exact_y_cell']))!=(Fraction(index,64),Fraction(index+1,64)):
                raise ValueError('Accepted density source partition changed')
            if not cell['branch_union_not_sum_or_duplicate_integral'] or not cell['whole_source_function_ranges_integrated_not_point_quadrature']:
                raise ValueError('Accepted whole-function branch union required')
            if {row['branch'] for row in cell['predicate_density_records']}!={'regular','positive','negative'}:
                raise ValueError('Axis and both signs of actual original predicate union required')
            groups=cell['normalized_four_density_unions']
            if len(groups)!=4 or any(set(group)!=set(KEYS) for group in groups):
                raise ValueError('Changed and original C0/Z density cover rows required')
            self.densities.append([{key:(self.Z_unit if j==1 else self.coordinates.scalar(1))*iv(c,value)
                for key,value in group.items()} for j,group in enumerate(groups)])
        self.origin=dict(self.direct['actual_original_radius_phase_endpoints'][0])
        if self.origin['source_family']!=self.family or self.origin['explicit_candidate_N']!=N or self.origin['original_y_exact']!='0':
            raise ValueError('Same original N1024 radius phase at slope y0 required')
        self.origin['true_original_phase_directed_boxes']=[
            dict(lower=mp.make_mpf(tuple(box['lower']['exact_mpf_tuple'])),
                upper=mp.make_mpf(tuple(box['upper']['exact_mpf_tuple'])))
            for box in self.origin['true_original_phase_directed_boxes']]
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.binding=dict(manifest=preceding.driver.NAME,manifest_sha256=sha(preceding.driver.NAME),
            receipt=preceding.driver.RECEIPT,receipt_sha256=sha(preceding.driver.RECEIPT),candidate_N=N,source_family=self.family,
            source_density_record='actual_fresh_original_N1024_direct_integrals.whole_source_density_phase_mass_records',
            original_source_cells=64,original_predicate_union=['regular','positive','negative'],
            exact_y_window=['0','1'],exact_Z_range=['-1','1'],original_common_coordinates=self.coordinates.record(),
            source_density_covers_restricted_by_measure_not_replaced_by_selected_points=True,
            source_covers_valid_on_every_contained_piece_without_replaying_source_queries=True,
            accepted_N160_ancestor_not_used_as_N1024_density=True,
            accepted_whole_bin_integrated_contributions_not_reused_as_piece_integrals=True,
            same_actual_thirteen_cell_incoming_correction_not_O2_exit_history=True,
            ordinary_Z_source_unit_retained_not_differentiated_as_cap=True)

    @native.inlet.source_precision
    def integrate(self,*,N=1024):
        if type(N) is not int or N!=1024:raise ValueError('Actual inherited covers require unchanged integer N1024')
        c=self.ctx;coords=self.coordinates;incoming=self.incoming;background=self.background
        cumulative=history.C1DuhamelOperator(coords);bg_operator=history.C1DuhamelOperator(coords)
        cells=[];livecells=[];piece_records=[];global_totals=[{key:coords.scalar(0) for key in KEYS} for unused in range(4)]
        with mp.workdps(c.dps+40):
            for label,a,b in zip(LABELS,EDGES[:-1],EDGES[1:],strict=True):
                totals=[{key:coords.scalar(0) for key in KEYS} for unused in range(4)];pieces=[]
                for index,left,right in exact_pieces(a,b):
                    width=rational(c,right-left);lo,hi=rational(c,left),rational(c,right)
                    local,final,decay=five.masses(c,width,lo,hi)
                    suffix={key:c.exp(-c.mpf(rate)*rational(c,b-right)) for key,rate in five.RATES.items()}
                    additions=[{key:value*local[key]*suffix[key] for key,value in group.items()} for group in self.densities[index]]
                    final_additions=[{key:value*final[key] for key,value in group.items()} for group in self.densities[index]]
                    for j in range(4):
                        for key in KEYS:
                            totals[j][key]+=additions[j][key];global_totals[j][key]+=final_additions[j][key]
                    row=dict(original_source_bin=index,exact_y_piece=[str(left),str(right)],exact_piece_width=str(right-left),
                        actual_phase_boxes_same_global_origin=five.pressure.phase_boxes(c,self.origin,N,lo,hi),
                        source_density_cover_domain=self.source_cells[index]['exact_y_cell'],original_partition_label=label,
                        positive_piece_local_endpoint_masses=local,piece_to_partition_endpoint_decay=suffix,
                        positive_piece_global_endpoint_masses=final,
                        four_piece_partition_endpoint_contributions=[records(group) for group in additions],
                        four_piece_global_endpoint_contributions=[records(group) for group in final_additions],
                        original_source_N1024_predicate_density_union_cover_index=index,
                        same_source_whole_bin_density_hull_valid_on_contained_piece=True,
                        actual_signed_whole_function_C0_Z_not_point_or_quadrature=True)
                    pieces.append(row);piece_records.append(row)
                width=rational(c,b-a);operator=history.C1DuhamelOperator(coords)
                operator.append(width,totals[0],totals[1],self.family)
                cumulative.append(width,totals[0],totals[1],self.family)
                original_operator=history.C1DuhamelOperator(coords);original_operator.append(width,totals[2],totals[3],self.family)
                bg_operator.append(width,totals[2],totals[3],self.family)
                correction=operator.apply(incoming['values'],incoming['Z_derivatives'],self.family)
                original=original_operator.apply(background['values'],background['Z_derivatives'],self.family)
                own={key:original['values'][key]+correction['values'][key] for key in KEYS}
                ownZ={key:original['Z_derivatives'][key]+correction['Z_derivatives'][key] for key in KEYS}
                geometry=self.source_owner.shell.transfer.geometry.cell('O2_slope',str(a),str(b))
                if ep(geometry['regular'])!=ep(width):raise ValueError('Actual original slope Jacobian must be exactly1')
                record=dict(label=label,chart='O2_slope',candidate_N=N,source_family=self.family,Z_box=c.mpf((-1,1)),
                    exact_y_partition=[str(a),str(b)],actual_true_cell_geometry=geometry['record'],source_bin_piece_records=pieces,
                    changed_C0_contributions=records(totals[0]),changed_Z_contributions=records(totals[1]),
                    original_C0_contributions=records(totals[2]),original_Z_contributions=records(totals[3]),
                    actual_cell_C1_operator=operator.record(),original_background_cell_C1_operator=original_operator.record(),
                    actual_inherited_correction_C0=records(incoming['values']),actual_inherited_correction_Z=records(incoming['Z_derivatives']),
                    actual_right_correction_C0=records(correction['values']),actual_right_correction_Z=records(correction['Z_derivatives']),
                    original_right_background_C0=records(original['values']),original_right_background_Z=records(original['Z_derivatives']),
                    actual_right_own_history_C0=records(own),actual_right_own_history_Z=records(ownZ),
                    original_separate_P0=self.P0.record(),original_separate_P0_Z=self.P0_Z.record(),
                    actual_absolute_pressure_C0=(self.P0+own['p']).record(),actual_absolute_pressure_Z=(self.P0_Z+ownZ['p']).record(),
                    exact_original_slope_cell_no_gap_or_duplicate_measure=True,
                    original_background_once_and_separate_P0=True,rate_zero_pressure_keeps_actual_upstream_memory=True,
                    all_original_phase_and_predicate_error_covers_retained=True,
                    restricted_density_hulls_are_function_covers_not_sharper_source_oracles=True)
                cells.append(record);livecells.append(dict(record=record,operator=operator,correction=correction,background=original))
                incoming=correction;background=original
                print('Actual full-predicate original O2 partition:',label,flush=True)
            if len(piece_records)!=68 or [cell['label'] for cell in cells]!=list(LABELS):
                raise ArithmeticError('Exactly68 disjoint pieces in the five original slope intervals required')
            report=dict(source_family=self.family,candidate_N=N,exact_Z_range=['-1','1'],exact_original_slope_edges=[str(edge) for edge in EDGES],
                accepted_actual_N1024_density_and_upstream_binding=self.binding,ordinary_Z_large_source_unit=self.Z_unit.record(),
                actual_original_five_slope_cell_records=cells,exact_source_bin_piece_count=len(piece_records),
                actual_five_slope_correction_C0_Z=[records(incoming['values']),records(incoming['Z_derivatives'])],
                actual_five_slope_background_C0_Z=[records(background['values']),records(background['Z_derivatives'])],
                actual_five_slope_cumulative_C1_operator=cumulative.record(),original_background_cumulative_C1_operator=bg_operator.record(),
                four_direct_piece_global_endpoint_totals=[records(group) for group in global_totals],
                accepted_full_spatial_O2_alternative_remains_separate_global_endpoint_only=True,
                whole_O2_spatial_IBP_bound_not_sliced_as_partial_cell_increment=True,
                actual_13_upstream_plus_five_O2_plus_six_Rc_records_not_yet_assembled=True,
                full24_original_C1_integral_range_transport_enclosed=False,actual_five_controls_installed=False,
                current_whole_N_selected=False,functional_terminal_identity_solved=False)
        return dict(report=report,cells=livecells,correction=incoming,background=background,cumulative=cumulative,global_totals=global_totals)


@native.inlet.source_precision
def run():
    began=time.monotonic();bridge,construction=native.inlet.native_bridge_owner()
    with native.inlet.CheckedSourceRuntime() as runtime:
        owner=OriginalO2FiveSlopePartitions(bridge);live=owner.integrate()
    report=dict(**{GATE:True},source_family=owner.family,candidate_N=N,actual_original_O2_five_slope_partitions=live['report'],
        original_bridge_construction=construction,original_source_runtime=runtime.record(),
        actual_five_original_O2_slope_partition_operators_installed=True,full24_original_C1_integral_range_transport_enclosed=False,
        actual_five_controls_installed=False,functional_terminal_identity_solved=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual N1024 full-predicate whole-Z signed density covers are restricted by exact measures at original five O2 slope boundaries. Fresh positive piece/partition masses and actual13-cell incoming correction; source backgrounds/ordinary-Z unit/P0 and phase origin retained. No double counting/replayed upstream or source-point substitution. Original24 record/control adapter, functional terminal/global N/recursion/full NS remain open.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(encode(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Actual full-predicate five original O2 slope operators generated',flush=True);return report


if __name__=='__main__':run()
