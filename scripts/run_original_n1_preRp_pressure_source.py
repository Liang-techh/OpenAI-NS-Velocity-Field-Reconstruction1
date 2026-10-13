"""Compute remaining original pressure and combine checked outer-pulse models."""
from pathlib import Path
import sys
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments/root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_source as current
import lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_source_check as check


def main():
    physical=current.physical
    modified=physical.partial.parent.CurrentOriginalN1ModifiedLeadingC3(physical.rg.CurrentC3TargetRanges())
    partial=physical.partial.CurrentOriginalN1PartialMeanQ(modified)
    atlas=physical.all_chart.CurrentOriginalN1AllChartPartialMeanQ(partial)
    fields=physical.CurrentOriginalN1AllChartPhysicalPressureFields(atlas)
    whole=current.pulse.CurrentOriginalN1PulsePressureSource()
    owner,raw=current.run(fields,whole)
    receipt=check.run(owner,raw)
    with mp.workdps(15):
        public=current.CurrentOriginalN1PreRpPressureSource(fields,whole)
        view=public.report(public.evaluate('.371'))
        assert view[current.GATE] and receipt[current.GATE] and mp.mp.dps==15
    print('ACCEPTED_ORIGINAL_RMINUS_TO_INFINITY_PRESSURE_PUBLIC',flush=True)


if __name__=='__main__':main()
