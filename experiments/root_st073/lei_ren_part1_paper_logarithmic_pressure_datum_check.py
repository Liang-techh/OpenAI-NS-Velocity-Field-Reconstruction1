"""Common-source overlap, global jet parity, and positive tail retention."""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_logarithmic_pressure_datum import LogarithmicPressureDatum,BETA2,BETA0
from lei_ren_part1_paper_Md11_pressure_datum import pressure_jets
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def run():
    name='lei_ren_part1_paper_large_Md_screen.json'
    old=json.loads((HERE/name).read_bytes());comparisons=[]
    with mp.workdps(240):
        for md in ('4','5','6'):
            new=LogarithmicPressureDatum(md);c=new.ctx
            trial=next(t for t in old['trials'] if t['Md']==md)
            prior=trial['fourteen_stage_trial_pressure_datum']
            jets=pressure_jets(c,c.mpf('.5'),8,prior)['physical_pressure_coefficients']
            scale=read_interval(c,prior['Pstar_squared'])
            fresh=new.normalized_jets(c.mpf('.5'),8)['normalized_pressure_coefficients']
            for n,(a,b) in enumerate(zip(jets,fresh)):
                al,au=endpoints(a/scale);bl,bu=endpoints(b)
                assert max(al,bl)<=min(au,bu)
            comparisons.append(dict(Md=md,common_source_coefficients_overlapping=9,
                overlap_is_consistency_check_not_containment_claim=True))
        datum=LogarithmicPressureDatum('40');c=datum.ctx
        assert len(datum.stages)==14 and len(BETA2)==6 and len(BETA0)==7
        assert endpoints(datum.tail_upper)[0]>0
        for name,row in datum.stages.items():
            assert endpoints(row['mass'])[0]>=0 and endpoints(row['mass'])[1]>0
        axis=datum.normalized_jets(c.mpf(0),127)['normalized_pressure_coefficients']
        assert all(endpoints(axis[n])==(mp.mpf(0),mp.mpf(0)) for n in range(1,128,2))
        for center in ('0','.5','1'):
            row=datum.normalized_jets(c.mpf(center),3)
            assert endpoints(row['normalized_pressure_coefficients'][0])[1]<0
        local=datum.normalized_jets(c.mpf(['.49','.51']),127)
        assert len(local['normalized_pressure_coefficients'])==128
        # Source identity is independent of the numerical enclosure precision.
        assert datum.source_sha==LogarithmicPressureDatum('40',precision=120).source_sha
        report=dict(common_source_checks=comparisons,stage_count=14,
            all_true_mass_boxes_nonnegative_with_positive_upper=True,
            strictly_positive_tail_retained=True,axis_odd_coefficients_exactly_zero=64,
            local_Taylor_order=127,global_real_center_values_negative=True,
            source_identity_independent_of_enclosure_precision=True,all_passed=True,
            local_C127_normalized_pressure=local,
            no_new_core_claimed=True,
            input_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in
                (Path(__file__).name,'lei_ren_part1_paper_logarithmic_pressure_datum.py',
                 'lei_ren_part1_paper_logarithmic_outer_parameters.py',
                 'lei_ren_part1_paper_Md11_pressure_datum.py',
                 'lei_ren_part1_paper_large_Md_screen.json')})
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(report)),indent=2)+'\n',encoding='utf-8')
        print('Logarithmic pressure common-source, positive-tail and C127 checks PASS',flush=True)
        return report


if __name__=='__main__':run()
