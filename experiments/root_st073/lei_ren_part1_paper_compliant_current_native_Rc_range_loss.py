"""Locate dominant fixed-N source and target range loss in the live route."""
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_native_cutoff_range_transport as current

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha;packets=current.packets;ep=current.ep
NAME=PREFIX+'current_native_Rc_range_loss.json'
GATE='current_original_fixed_N_Rc_dominant_range_loss_inventory_executed'


def upper(value):
    return None if value.zero else ep(value.record()['log_absolute_upper'])[1]


def run(owner,live):
    began=time.monotonic();checked=json.loads((HERE/current.RECEIPT).read_bytes())
    if not checked.get('all_passed') or not checked.get(current.GATE) or checked['source_family']!=owner.family:
        raise ValueError('Checked same-original-family complete fixed-N range route required')
    if len(live['cells'])!=24 or live['target'] is None or not live['record']['full24_original_C1_integral_range_transport_enclosed']:
        raise ValueError('Actual complete live original24-cell route required')
    hashes={**checked['input_hashes'],current.RECEIPT:sha(current.RECEIPT),Path(__file__).name:sha(Path(__file__).name)}
    for name,digest in hashes.items():
        if sha(name)!=digest:raise ValueError('Changed range-loss prerequisite: '+name)
    c=owner.ctx;coords=owner.coordinates;target=live['target'];suffix={key:coords.scalar(1) for key in current.RATES}
    rows=[];dominant={};primitive_dominant={};primitive_count=0
    for cell in reversed(live['cells']):
        values={key:cell['values'][key]*suffix[key] for key in current.RATES}
        jets={key:cell['Z_derivatives'][key]*suffix[key] for key in current.RATES}
        got=current.fixed_N_target_rows(values,jets,target['A'],target['AZ'],target['logA'],coords.scalar(target['mu']),target['logmu'])
        label=cell['record']['label'];per={}
        for group in ('values','Z_derivatives'):
            for key,value in got[group].items():
                cap=upper(value);name=key+('' if group=='values' else '_Z')
                per[name]=dict(exact_zero=value.zero,log_absolute_upper=None if cap is None else c.mpf(cap))
                if cap is not None and (name not in dominant or cap>dominant[name][0]):dominant[name]=(cap,label)
        primitives={}
        if cell['frame'] is not None:
            source=cell['frame'].record['actual_original_spatial_source']
            for branch in source['original_cutoff_branch_density_queries']:
                for phase in branch['original_conditional_spatial_density']['spatial_signed_density_branch_cells']:
                    for jet in phase['original_conditional_branch_first_jets']:
                        for key in ('A','A_Z','B_over_Pstar','B_Z_over_Pstar'):
                            value=jet['original_A_B_first_derivative_enclosures'][key]
                            if value['exact_zero']:continue
                            cap=ep(packets.interval(c,value['log_absolute_upper']))[1];primitive_count+=1
                            if key not in primitives or cap>ep(primitives[key])[1]:primitives[key]=c.mpf(cap)
                            if key not in primitive_dominant or cap>primitive_dominant[key][0]:primitive_dominant[key]=(cap,label)
        rebase={}
        if cell['frame'] is not None:
            for group in ('C0','Z'):
                for key,value in cell['frame'].values[group].items():
                    bound=coords.rebase(value,owner.family);before,after=upper(value),upper(bound)
                    rebase[group+'_'+key]=dict(exact_zero=value.zero,
                        native_log_absolute_upper=None if before is None else c.mpf(before),
                        common_log_absolute_upper=None if after is None else c.mpf(after),
                        upper_change=None if before is None else c.mpf(after-before))
        rows.append(dict(label=label,chart=cell['record']['chart'],normalized_final_target_contribution_ranges=per,
            original_phase_primitive_maximum_log_upper=primitives,source_to_common_basis_rebase=rebase,
            final_suffix_decay={key:value.record() for key,value in suffix.items()}))
        suffix={key:cell['factors'][key]['decay']*suffix[key] for key in current.RATES}
    Amax=primitive_dominant.get('A')
    result=dict(source_family=owner.family,**{GATE:True},candidate_N=live['record']['candidate_N'],
        original_continuous_cells=24,weighted_original_cell_ranges=list(reversed(rows)),
        dominant_normalized_target_contribution_cells={key:dict(label=label,log_absolute_upper=c.mpf(cap)) for key,(cap,label) in dominant.items()},
        dominant_original_primitive_range_cells={key:dict(label=label,log_absolute_upper=c.mpf(cap)) for key,(cap,label) in primitive_dominant.items()},
        original_phase_primitive_ranges_inspected=primitive_count,
        known_original_A_support_log_upper=c.ln(c.mpf(5)/4),
        original_A_enclosure_exceeds_known_support=Amax is not None and Amax[0]>ep(c.ln(c.mpf(5)/4))[1],
        source_enclosure_width_not_actual_NS_residual=True,
        rank_uses_range_upper_coordinates_not_field_values=True,
        per_cell_target_hulls_do_not_prove_joint_cancellation=True,
        next_production='Intersect original C0 primitives with proved source support while retaining tighter tiny formal factors and untouched derivative rows; then compare actual continuous route targets.',
        **dict.fromkeys(packets.OPEN,False),execution_seconds=time.monotonic()-began,input_hashes=hashes)
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Dominant actual original range cells:',{key:row['label'] for key,row in result['dominant_normalized_target_contribution_cells'].items()},flush=True)
    return result
