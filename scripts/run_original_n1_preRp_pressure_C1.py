"""Consume checked Fourier C3 inputs in original pre-Rp pressure C1."""
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_C1 as current
import lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_C1_check as check


def main():
    owner,_=current.run()
    receipt=check.run(owner)
    if not receipt['all_passed'] or not receipt[current.UPPER_GATE]:
        raise RuntimeError('Original upper pressure C1 rejected')
    accepted=current.CurrentOriginalN1PreRpPressureC1()
    accepted.derivative_enclosure('.371')
    print('Checked original r_minus-to-infinity pressure derivatives.',flush=True)


if __name__=='__main__':main()
