"""Compose/check complete original pressure-moment first-Z derivatives."""
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_complete_pressure_C1 as current
import lei_ren_part1_paper_compliant_current_original_n1_complete_pressure_C1_check as check


def main():
    owner,_=current.run()
    receipt=check.run(owner)
    if not receipt['all_passed'] or not receipt[current.GATE]:
        raise RuntimeError('Complete original pressure moment C1 rejected')
    accepted=current.CurrentOriginalN1CompletePressureC1()
    accepted.derivative_enclosure('.371')
    print('Actual signed defect functions and n1 repair coefficients remain open.',flush=True)


if __name__=='__main__':main()
