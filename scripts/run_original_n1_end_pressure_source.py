"""Compute/check original end pressure without unchanged ancestor reruns."""
from pathlib import Path
import sys
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments/root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_end_pressure_source as current
import lei_ren_part1_paper_compliant_current_original_n1_end_pressure_source_check as check


def main():
    owner,samples,models=current.run()
    receipt=check.run(owner,samples,models)
    with mp.workdps(15):
        public=current.CurrentOriginalN1EndPressureSource(owner.coefficients)
        view=public.report(public.evaluate('.371'))
        assert view[current.GATE] and receipt[current.GATE] and mp.mp.dps==15
    print('ACCEPTED_ORIGINAL_N1_END_PRESSURE_PUBLIC',flush=True)


if __name__=='__main__':main()
