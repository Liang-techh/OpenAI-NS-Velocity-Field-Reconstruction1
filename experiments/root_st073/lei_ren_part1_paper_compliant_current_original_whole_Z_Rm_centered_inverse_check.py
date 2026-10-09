"""Independent centering identities, common-box inverse and actual source audit."""
import gzip
import json
from pathlib import Path
import time
from unittest.mock import patch
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_centered_inverse as current
from lei_ren_part1_paper_compliant_current_original_whole_Z_bridge_source_check import same_source, forbidden
from lei_ren_part1_paper_compliant_current_original_whole_Z_switch_finite_N_check import exact_row
from lei_ren_part1_paper_compliant_current_original_Rm_defect_patch_inverse_check import exact_identities


def encoded(value):return current.source.bridge._encode(current.upstream.reference.serialized(current.leading.pack(value)))


def independent_centering_identities():
    z,E,mu,aa,H,Hi,Ki,V=sy.symbols('z E mu aa H Hi Ki V')
    actual_V=4*z+E;M=actual_V+mu;A=actual_V**2+aa
    assert sy.expand((M-4*z)-(E+mu))==0
    assert sy.expand((A-8*z*M+16*z*z)-(E*E+aa-8*z*mu))==0
    kernel=sy.symbols('kernel');actual_H=Hi+kernel;actual_K=Ki+actual_V*kernel
    assert sy.expand((actual_K-4*z*actual_H)-(E*actual_H+Ki-actual_V*Hi))==0
    # The macro driver is identical after changing only its affine seed.
    j,drive=sy.symbols('j drive');n=sy.symbols('n',integer=True,positive=True)
    rows=sy.Function('U');radial=sy.Sum(rows(n)*4**n,(n,1,24))
    assert sy.expand((4*z+j+radial+drive-4*z)-(j+radial+drive))==0
    return dict(passed=True,exact_actual_source_centering_identities=4,
        normalization_and_original_five_map_identities=exact_identities())


def run():
    began=time.monotonic();saved=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    assert saved[current.GATE] and saved['whole_Z_leading_Rm_implicit_controls_solved']
    assert saved['common_box_unique_smooth_leading_family'] and all(saved[key] is False for key in current.OPEN)
    assert not saved['full_radial_mixed4_leading_patch_installed']
    assert not saved['current_same_N_finite_correction_after_Rm_installed']
    for name,digest in saved['input_hashes'].items():assert current.sha(name)==digest,'Source changed: '+name
    identities=independent_centering_identities()
    counts=dict(whole_Z_cells=0,exact_directed_interval_rows=0,actual_factored_control_rows=0,
        normalized_defect_rows=0,common_unit_boxes=0,source_defined_pressure_joins=0,
        same_Jacobian_higher_orders=0,actual_patch_endpoint_queries=0,
        original_terminal_functional_identities=0,typed_rejections=0)
    margins={};defect_comparison={}
    with mp.workdps(540),patch.object(current.leading,'OriginalRmDefectPatchInverse',forbidden), \
            patch.object(current.upstream.endpoint,'OriginalLongReshapeEndpoint',forbidden), \
            patch.object(current.upstream.background,'OriginalReferenceRestoreFunctions',forbidden):
        owner=current.WholeZRmCenteredInverse();c=owner.c;right=mp.mpf(-1)
        assert saved['source_family']==owner.identity and saved['original_source_bindings']==encoded(owner.bindings)
        assert saved['candidate_N_not_used_for_leading_map']==owner.N
        for ends,cell in zip(current.source.CELLS,saved['source_cells']):
            op=owner.owner(ends);p=op.patch;f=p.flow
            assert current.ep(op.reference.Z)[0]==right;right=current.ep(op.reference.Z)[1]
            assert cell['exact_Z_cell']==list(ends) and cell['source_identity']==owner.identity
            same_source(encoded(op.source_centering_proof),cell['source_centering_proof'],counts)
            same_source(encoded(p.coefficients()),cell['coefficient_solution'],counts)
            inv=p.inverse['initial_C1_inverse']
            assert inv['certified'] and inv['uniform_implicit_C1_family_exists'] and inv['unique_in_initial_box']
            assert inv['self_map_strictly_inside'] and current.ep(inv['contraction_bound'])[1]<1
            assert current.ep(p.inverse['higher_inverse_contraction'])[1]<1
            assert all(current.ep(box)==(-1,1) for box in inv['initial_box'])
            assert all(current.ep(scale)==(1,1) for scale in inv['scales'])
            counts['common_unit_boxes']+=1
            assert [proof['order'] for proof in p.inverse['higher_order_proofs']]==[2,3,4,5]
            for proof in p.inverse['higher_order_proofs']:
                assert proof['inverse_contraction']._mpi_==p.inverse['higher_inverse_contraction']._mpi_
                counts['same_Jacobian_higher_orders']+=1
            # All cells use the SAME source equation and SAME unit box.
            # At shared endpoints uniqueness in that common box identifies
            # the local smooth solutions. Interval overlap is not the proof.
            assert p.P0 is op.reference.P0 and p.P0 is op.upstream.reference.P0
            assert p.P0 is p.Rm['original_P0_normalized_axial5'];counts['source_defined_pressure_joins']+=3
            assert p.angularUnit.scale.powers==(0,-1,0,0,0) and current.ep(p.angularUnit.coefficient)[0]>0
            assert all(v.scale.bases is f.logs and v.ledger is f.ledger for row in p.controls for v in row)
            assert p.defects[0] is p.centered['mean_error'] and p.defects[1] is p.centered['mixed_error']
            assert p.defects[2] is p.centered['angular_error']
            for i,row in enumerate(p.controls):
                for n,value in enumerate(row):
                    exact_row(f,value,cell['coefficient_solution']['actual_factored_coefficient_axial5'][i][n])
                    counts['actual_factored_control_rows']+=1
            for i,jet in enumerate(p.normalized_defects):
                for n,value in enumerate(jet.coefficients):
                    assert value._mpi_==current.read(c,cell['coefficient_solution']['normalized_five_defect_axial5'][i]['coefficients'][n])._mpi_
                    counts['normalized_defect_rows']+=1
            for row,record in zip(p.P0,cell['exact_common_P0_axial5']):exact_row(f,row,record)
            # Compare the tightened actual-centered source with the old
            # broad extension. An overlap is only a consistency check.
            broad=op.upstream.reference.postrestore((-6,1))['actual_centered_histories']
            broad_defects=[broad['mean_error'],broad['mixed_error'],broad['angular_error'],
                f.add(f.multiply(broad['axial_square'],p.invAm2),f.scale(broad['swirl_error'],-c.mpf('.5'))),
                f.scale(broad['pressure_error'],c.mpf('.5'))]
            old=[];new=[]
            for i,row in enumerate(broad_defects):
                divided=f.ordinary_cover(row[0].positive_divide(p.units[i],p.units[i].scale.evaluate()+c.ln(p.units[i].coefficient)))
                old.append(c.mpf(max(abs(v) for v in current.ep(divided))))
                new.append(c.mpf(max(abs(v) for v in current.ep(p.normalized_defects[i][0]))))
                for n in range(6):
                    difference=row[n]-p.defects[i][n]
                    assert difference.zero or current.ep(difference.coefficient)[0]<=0<=current.ep(difference.coefficient)[1]
            assert current.ep(new[0])[1]<current.ep(old[0])[0]
            assert current.ep(new[3])[1]<current.ep(old[3])[0]
            defect_comparison[str(ends)]=dict(old_normalized_C0_abs_upper=old,new_normalized_C0_abs_upper=new)
            for q,record in zip(((1,1),(2,1)),cell['original_partial_patch_endpoint_queries']):
                live=p.evaluate(q);same_source(encoded(live),record,counts)
                assert live['original_P0_normalized_axial5'] is p.P0
                assert live['same_original_P0_retained'] and live['source_caps_not_used_as_field_values']
                if q==(1,1):
                    assert all(v.zero for row in live['signed_bump_correction_sectors']['five_partial_primitive_changes'] for v in row)
                else:
                    assert live['exact_terminal_identity_of_same_unique_leading_map']
                    assert live['terminal_identity_not_an_inlet_reset_or_zero_containment_proof']
                    assert all(v.zero for row in live['actual_recovered_five_defects'] for v in row)
                    assert current.ep(live['normalized_primitives']['theta'][0].coefficient)[0]>0
                    counts['original_terminal_functional_identities']+=5
                counts['actual_patch_endpoint_queries']+=1
            margins[str(ends)]=dict(initial_contraction=inv['contraction_bound'],
                higher_Jacobian_contraction=p.inverse['higher_inverse_contraction'])
            counts['whole_Z_cells']+=1
            print('Whole-Z centered actual defects/common-box controls: '+str(ends),flush=True)
        assert right==1
        for call in (lambda:owner.owner(('0','0')),lambda:owner.owner(('-2','-1'))):
            try:call()
            except ValueError:counts['typed_rejections']+=1
            else:raise AssertionError('Unadmitted axial source cell must reject')
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.identity,
        genuine_source_centering_not_zero_target_substitution=True,
        core_coupled_rows_tail_and_all_bridge_switch_increments_retained=True,
        same_common_unit_box_proves_unique_whole_Z_smooth_leading_family=True,
        actual_leading_controls_and_same_Jacobian_axial5_checked=True,
        independent_source_centering_and_original_map_algebra=identities,
        actual_partial_patch_endpoint_functions_and_same_P0_checked=True,
        residual_zero_containment_not_used_as_functional_closure_proof=True,
        whole_Z_leading_Rm_implicit_controls_solved=True,contraction_margins=encoded(margins),
        actual_centered_defect_sharpening=encoded(defect_comparison),replay_counts=counts,
        full_radial_mixed4_leading_patch_installed=False,current_same_N_finite_correction_after_Rm_installed=False,
        **dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Whole-Z actual source-centered leading Rm controls checks passed',flush=True);return result


if __name__=='__main__':run()
