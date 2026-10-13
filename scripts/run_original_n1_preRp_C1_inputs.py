"""Restore accepted owners once and build/check new pre-Rp C1 inputs."""
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import run_original_n1_five_defect_values as restore
import lei_ren_part1_paper_compliant_current_original_n1_preRp_C1_inputs as current
import lei_ren_part1_paper_compliant_current_original_n1_preRp_C1_inputs_check as check


def main():
    cache={}
    fields=restore.build_checked(current.physical.CurrentOriginalN1AllChartPhysicalPressureFields,cache)
    whole=restore.build_checked(current.pulse.CurrentOriginalN1PulsePressureSource,cache)
    owner,packets=current.run(fields,whole)
    receipt=check.run(owner,packets)
    if not receipt['all_passed'] or not receipt[current.GATE]:
        raise RuntimeError('Original pre-Rp C1 inputs rejected')
    print('Checked source input records: 68 baseline unions and 4 compact cells.',flush=True)


if __name__=='__main__':main()
