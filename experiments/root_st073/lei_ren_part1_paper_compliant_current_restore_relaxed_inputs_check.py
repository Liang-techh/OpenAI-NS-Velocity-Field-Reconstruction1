"""Focused actual source, full signed N/bw and whole restore-domain checks."""
import json
from pathlib import Path
import lei_ren_part1_paper_compliant_current_restore_relaxed_inputs as source


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed actual restoration source: '+name)
    owner=source.CurrentRestoreRelaxedInputs()
    if source.packets.encode(owner.theorem)!=data['original_current_full_signed_restoration_source_theorem']:
        raise ValueError('Complete signed restoration source theorem changed')
    if source.packets.encode(owner.proof)!=data['current_whole_Rz_Rm_relaxed_input_and_bounds']:
        raise ValueError('Actual source pressure/N/relaxed bound changed')
    if owner.conditions!=data['current_source_function_attachment_conditions']:raise ValueError('Actual source function conditions changed')
    queries={chart:owner.query(chart,(0,1) if chart=='axial_restore' else (-7,-6)) for chart in source.CHARTS}
    if source.packets.encode(queries)!=data['analytic_source_subbox_examples']:raise ValueError('Actual original source queries changed')
    bad=(('actual_patch',1,0),('axial_restore',-1,0),('axial_restore',2,0),
        ('restore_buffer',-5,0),('restore_buffer','inf',0),('axial_restore','.5',2))
    for chart,coordinate,Z in bad:
        try:owner.query(chart,coordinate,Z)
        except ValueError:pass
        else:raise ArithmeticError('Invalid original restore/buffer source query admitted')
    c=owner.ctx;read=lambda value:source.packets.interval(c,value)
    proof=owner.proof;join=owner.reshape.rows['reference_join_bounds']
    eps=read(owner.reshape.inner.rows['K1_ledger']['epsilon0'])
    if source.packets.recovery.endpoints(proof['original_signed_b_times_Pstar_absolute_upper'])[1]>=source.packets.recovery.endpoints(144*eps)[0]:
        raise ArithmeticError('Independent original b*Pstar normalization budget failed')
    if source.packets.recovery.endpoints(proof['original_signed_bw_absolute_upper'])[1]>source.packets.recovery.endpoints(read(join['axial_restoration_bw_upper']))[1]:
        raise ArithmeticError('Same current original signed bw budget enlarged')
    if any(data[k] or proof[k] for k in source.OPEN) or proof['active_patch_relaxed_input_certified']:
        raise ArithmeticError('Original restore gate promoted an active patch/whole modified target')
    r=dict(all_passed=True,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),source_family=owner.family,
        current_source_function_attachment_conditions=len(owner.conditions),
        exact_original_restore_and_terminal_AST_bindings=len(owner.theorem['original_source_AST_bindings']),
        exact_full_signed_N_bw_and_H0_source_identities=len(owner.theorem['exact_identities']),
        positive_actual_pressure_inherited_N_norm_gates=len(proof['actual_current_pressure_and_inherited_N_norm_gates']),
        positive_actual_relaxed_margins=len(proof['positive_directed_relaxed_margins']),
        original_b_Pstar_and_bw_normalization_budgets_checked=2,
        actual_source_queries_checked=len(queries),invalid_source_and_patch_queries_rejected=len(bad),
        source_radius_amplitude_and_width_not_materialized=True,source_ancestor_constructors_called=False,
        full_signed_p2_pressure_energy_and_meridional_terms_retained=True,
        p1_p2_whole_path_norm_bounds_certified=False,strict_completed_tensor_cone_new_regions_admitted=0,
        active_patch_relaxed_input_certified=False,scope=data['scope'],
        input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(r,indent=2)+'\n',encoding='utf8')
    print('Current actual full signed restoration/buffer relaxed input and source gates PASS',flush=True)
    return r


if __name__=='__main__':run()
