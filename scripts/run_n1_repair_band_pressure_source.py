"""Reproduce only the actual compact repair-band pressure calculation."""
from pathlib import Path
import sys
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments/root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_repair_band_pressure_source as current
import lei_ren_part1_paper_compliant_current_original_n1_repair_band_pressure_source_check as check


def main():
    with mp.workdps(540):
        owner,packets=current.run()
        receipt=check.run(owner,packets)
    public=current.CurrentOriginalN1RepairBandPressureSource()
    view=public.report(public.evaluate(('-1','-.5')))
    assert view[current.GATE] and receipt[current.GATE]
    print('ACCEPTED_REPAIR_BAND_PUBLIC_QUERY',mp.mp.dps,flush=True)


if __name__=='__main__':main()
