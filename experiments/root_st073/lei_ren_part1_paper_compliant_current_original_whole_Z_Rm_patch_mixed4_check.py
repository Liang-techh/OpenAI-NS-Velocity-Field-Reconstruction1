"""Audit full axial/radial coverage, original mixed algebra and actual joins."""
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_patch_mixed4 as current
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_whole_Z_switch_finite_N_check import exact_row
from lei_ren_part1_paper_compliant_current_original_Rm_patch_mixed4_cells_check import symbolic_sources, fixture


def encoded(value):return current.upstream.source.bridge._encode(current.serialized(value))


def source_contract(op,value,counts):
    p=op.op;f=p.flow
    assert value['whole_Z_original_leading_patch_function_provider_installed']
    assert not value['actual_patch_finite_N_density_oracle_installed']
    assert value['original_P0_normalized_axial5'] is p.P0
    assert value['raw_radial_Z_order4_only_no_selected_Z5']
    assert value['physical_prefactors_differentiated_before_grid']
    assert value['primitive_source_index_not_shifted']
    assert value['common_R_derivative_grid_is_x_grid_with_unit_Rm_power_minus_k']
    assert value['common_pressure_increment_grid_is_Mp_with_exact_same_separate_P0']
    Q=value['actual_Q_x_derivative_axial4']
    assert len(Q)==5 and all(len(row)==5 for row in Q)
    keys={str(k)+'_Z'+str(n) for k in range(5) for n in range(5-k)}
    for coordinate in ('x','y','R'):
        groups=(value['physical_velocity_'+coordinate+'_Z_mixed4'],
            {'P':value['physical_pressure_'+coordinate+'_Z_mixed4']},
            value['physical_five_primitive_'+coordinate+'_Z_mixed4'])
        for group in groups:
            for grid in group.values():
                assert {key[1:] for key in grid}==keys
                assert all(v.scale.bases is f.logs and v.ledger is f.ledger for v in grid.values())
                counts['physical_mixed4_source_rows']+=len(grid)
    assert len(value['raw_current_radius_y_derivative_axial_coefficients']['histories'])==5
    assert len(value['raw_current_radius_y_derivative_axial_coefficients']['velocity']['radial'][0])==5
    counts['independent_P0_object_joins']+=1


def run():
    began=time.monotonic();saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE] and saved['full_radial_mixed4_leading_patch_installed']
    assert saved['whole_Z_leading_Rm_implicit_controls_solved']
    assert not saved['current_same_N_finite_correction_after_Rm_installed']
    assert all(saved[key] is False for key in current.OPEN)
    for name,digest in saved['input_hashes'].items():assert current.sha(name)==digest,'Source changed: '+name
    counts=dict(whole_Z_cells=0,closed_radial_cells=0,fresh_actual_source_queries=0,
        physical_mixed4_source_rows=0,exact_directed_interval_rows=0,
        independent_P0_object_joins=0,exact_flat_support_point_rows=0,
        actual_Rm_Rh_background_join_rows=0,typed_rejections=0)
    with mp.workdps(540):
        # Independent nonzero, closed-integral fixtures check the original
        # primitive source index and all x/y/R derivative conversions.
        symbolic=symbolic_sources()
        independent=fixture(symbolic,'1.53','.2')
        assert independent['passed']
        independent['exact_closed_primitive_and_coordinate_identities']=symbolic['exact_closed_primitive_and_coordinate_identities']
    with mp.workdps(540),patch.object(current.mixed,'OriginalRmPatchMixed4Cells',forbidden), \
            patch.object(current.upstream.leading,'OriginalRmDefectPatchInverse',forbidden):
        owner=current.WholeZRmPatchMixed4();c=owner.c;right=mp.mpf(-1)
        assert saved['source_family']==owner.identity and saved['original_source_bindings']==encoded(owner.bindings)
        assert saved['exact_x_partition']==encoded(current.PARTITION)
        for ends,packet in zip(current.upstream.source.CELLS,saved['source_cells']):
            op=owner.owner(ends);p=op.op;f=op.flow
            assert current.ep(p.reference.Z)[0]==right;right=current.ep(p.reference.Z)[1]
            assert packet['exact_Z_cell']==list(ends) and packet['source_identity']==owner.identity
            assert p.P0 is p.reference.P0 and p.P0 is p.Rm['original_P0_normalized_axial5']
            for row,record in zip(p.P0,packet['exact_common_P0_axial5']):exact_row(f,row,record)
            xright=Fraction(1)
            for left,right_,record in zip(current.PARTITION,current.PARTITION[1:],packet['source_radial_cells']):
                assert Fraction(*left)==xright
                live=owner.query(ends,left,right_);same_source(encoded(live),record,counts)
                assert live['closed_x_cell_not_a_sample_interpolation'] and not live['geometry']['point']
                assert live['actual_partial_primitive_source_memory']['inlet_never_reset']
                source_contract(op,live,counts)
                xright=None if right_=='Rh' else Fraction(*right_)
                terminal=Fraction(*left)>=Fraction(71,40)
                assert live['exact_terminal_mean_and_radial_Q_refinement_from_same_map']==terminal
                if terminal:
                    assert all(current.ep(v)==(0,0) for rows in live['actual_gamma_ordinary_x_derivatives'] for v in rows)
                    assert all(v.zero for row in live['actual_Q_x_derivative_axial4'][1:] for v in row)
                counts['closed_radial_cells']+=1;counts['fresh_actual_source_queries']+=1
            assert xright is None
            for coordinate,key in (((1,1),'actual_Rm_function'),('Rh','actual_Rh_function')):
                live=owner.query(ends,coordinate);same_source(encoded(live),packet[key],counts)
                source_contract(op,live,counts);counts['fresh_actual_source_queries']+=1
                if coordinate==(1,1):
                    parent=p.evaluate((1,1))
                    mapping=dict(Utheta='Utheta',Uz='Uz',Ur='Ur')
                    for name,oldname in mapping.items():
                        for a,b in zip(live['physical_velocity_x_derivative_axial_coefficients'][name][0],parent['physical_velocity_axial_coefficients'][oldname]):
                            difference=a-b
                            assert difference.zero or current.ep(difference.coefficient)[0]<=0<=current.ep(difference.coefficient)[1]
                            counts['actual_Rm_Rh_background_join_rows']+=1
                else:
                    assert live['actual_partial_primitive_source_memory']['exact_terminal_identity_of_same_leading_map']
            exact_row(f,p.Rm_factor,packet['exact_Rm_radius'])
            exact_row(f,p.Rm_factor*c.exp(1),packet['exact_Rh_radius'])
            # Original beta is flat at every exact rational support edge.
            for edge in ((49,40),(51,40),(59,40),(61,40),(69,40),(71,40)):
                value=owner.query(ends,edge)
                assert all(current.ep(v)==(0,0) for rows in value['actual_gamma_ordinary_x_derivatives'] for v in rows)
                counts['exact_flat_support_point_rows']+=15
            counts['whole_Z_cells']+=1
            print('Whole-Z actual closed mixed4 patch atlas audit: '+str(ends),flush=True)
        assert right==1
        for call in (lambda:owner.owner(('0','0')),lambda:owner.query(('-1','-.5'),(0,1)),
            lambda:owner.query(('-1','-.5'),(2,1),(1,1)),lambda:owner.query(('-1','-.5'),(3,1))):
            try:call()
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Unadmitted axial/radial domain must reject')
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,
        full_closed_real_axial_and_Rm_Rh_radial_cover_checked=True,
        actual_current_source_control_and_mixed4_replays_checked=True,
        independent_nonzero_closed_primitive_x_y_R_mixed4_fixture=encoded(independent),
        original_flat_edges_same_P0_Rm_Rh_and_formal_units_checked=True,
        finite_N_correction_not_substituted_for_leading_source=True,
        full_radial_mixed4_leading_patch_installed=True,current_same_N_finite_correction_after_Rm_installed=False,
        actual_patch_finite_N_density_oracle_installed=False,replay_counts=counts,
        **dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Whole-Z actual leading Rm..Rh mixed4 atlas checks passed',flush=True);return result


if __name__=='__main__':run()
