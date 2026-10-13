"""Compose saved checked lower and upper pressure models; no re-integration."""
from pathlib import Path
import sys
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments/root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_external_pressure_values as current
import lei_ren_part1_paper_compliant_current_original_n1_external_pressure_values_check as check


def main():
    owner,raw=current.run()
    receipt=check.run(owner,raw)
    with mp.workdps(15):
        public=current.CurrentOriginalN1ExternalPressureValues()
        view=public.report(public.evaluate('.371'))
        assert view[current.GATE] and receipt[current.GATE] and mp.mp.dps==15
    print('ACCEPTED_COMPLETE_EXTERNAL_PRESSURE_QUERY',flush=True)


if __name__=='__main__': main()
