"""Actual common-source value/first-Z Picard inverse, without uniform claims."""
import json
import pickle
from types import SimpleNamespace
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_five_moment_reference_background import ReferenceDefectBackground
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_inverse_field import FiveBumpInverseField
from lei_ren_part1_paper_axial_correction import signed_log


def run():
    cache=Path(__file__).resolve().parents[2]/'work_paper_cache'/'reference_defect_Z03.pkl'
    with cache.open('rb') as stream:snapshot=pickle.load(stream)
    if snapshot['version']!=1 or snapshot['Z']!='.3' or snapshot['x']!='1.25':
        raise ValueError('Unsupported actual source snapshot')
    with mp.workdps(snapshot['precision']):
        data=snapshot['data'];baseline=snapshot['baseline']
        moment_map=FiveBumpMomentMap(precision=snapshot['precision'],order=96)
        def at_Z(z):
            if mp.mpf(z)!=mp.mpf('.3'):raise ValueError('Snapshot supports Z=.3 only')
            return data
        def at_x(x,z):
            at_Z(z)
            if mp.mpf(x)!=mp.mpf('1.25'):raise ValueError('Snapshot baseline supports x=1.25 only')
            return baseline
        reference=SimpleNamespace(Rm=snapshot['Rm'],delta=snapshot['delta'],
                                  defect_data=at_Z,evaluate_x=at_x)
        adapter=FiveBumpInverseField(reference,moment_map,steps=10)
        field=adapter.evaluate_x('1.25','.3')
        _,result=adapter.inverse_data('.3')
        assert field['P0']==baseline['P0'] and field['P0_Z']==baseline['P0_Z']
        assert field['stress'] is not None,field['stress_error']
        direct=moment_map.apply(result['h'],adapter.inverse_data('.3')[0])
        errors=[]
        for i in range(5):
            difference=direct[i]+data['d_dual'][i+1]-result['terminal_residual'][i]
            errors.extend(abs(v) for jet in (difference.value,difference.tangent) for v in jet.atoms.values())
        replay_error=max(errors,default=mp.mpf(0))
        assert replay_error<mp.mpf('1e-400'), 'Actual retained-jet residual replay mismatch'
        def encode_jet(jet):return {str(k):signed_log(v,70) for k,v in jet.atoms.items()}
        def encode_dual(dual):return dict(value=encode_jet(dual.value),first_Z=encode_jet(dual.tangent))
        history=[]
        for receipt in result['iterations']:
            history.append(dict(update=receipt['update'],
                increment=[encode_dual(v) for v in receipt['increment']],
                residual=[encode_dual(v) for v in receipt['residual']]))
            maximum=max(abs(v.value.evaluate()) for v in receipt['residual'])
            history[-1]['max_nominal_residual_value']=signed_log(maximum,70)
            history[-1]['max_nominal_residual_first_Z']=signed_log(
                max(abs(v.tangent.evaluate()) for v in receipt['residual']),70)
            print('update',receipt['update'],'max nominal centered residual',mp.nstr(maximum,12),flush=True)
        report=dict(Z='.3',x='1.25',precision=snapshot['precision'],pressure_order=9,width_order=2,
            inverse_method='incremental Picard on actual common-source five-bump map',
            nonlinear_updates=10,initial_linear_response=[encode_dual(v) for v in result['increments'][0]],
            iterations=history,terminal_residual=[encode_dual(v) for v in result['terminal_residual']],
            corrected_fields={k:encode_jet(field[k]) for k in ('F','Uz','P','P0','Ur')},
            P0_preserved=True,first_Z_from_same_source=True,source_part_count=69,
            callable_inverse_adapter_used=True,
            same_map_aggregate_replay_max_atom_error=signed_log(replay_error,70),
            inverse_increment_parts_separately_retained=True,
            terminal_residual_scope='direct polarization, not subtraction of physical moments',
            field_materialization_can_round_tiny_parts=True,
            scalar_nominal_residuals_not_relative_flat_source_bounds=True,
            uniform_Z_certified=False,infinite_response_error_enclosed=False,
            source_error_enclosed=False,quadrature_enclosed=False,
            functional_closure=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        return report


if __name__=='__main__':run()
