"""Actual common preheat finite-atom interval envelope; not quadrature certification."""
import json
from pathlib import Path
from lei_ren_part1_paper_continuous_preheat_pressure_check import source_profile,serialize
from lei_ren_part1_paper_continuous_preheat_pressure import ContinuousPreheatPressure
from lei_ren_part1_paper_preheat_interval_bounds import preheat_interval_bounds
import mpmath as mp


def run():
    print('building actual common preheat interval envelope',flush=True)
    source=source_profile();adapter=ContinuousPreheatPressure(source,quadrature_order=192)
    with mp.workdps(adapter.precision):
        receipt=adapter.taylor_components(2,center=0)
        bounds=preheat_interval_bounds(receipt,axial_radius='.8')
        assert len(bounds['atoms'])>=192
        assert any(v['region']=='post_Rv' and v['mass']>0 for v in bounds['atoms'].values())
        Path(__file__).with_suffix('.json').write_text(json.dumps(serialize(bounds),indent=2)+'\n',encoding='utf-8')
        print('actual finite preheat atom envelope saved',bounds['atom_count'],flush=True)


if __name__=='__main__':run()
