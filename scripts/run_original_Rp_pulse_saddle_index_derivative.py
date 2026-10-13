"""Compute/check the actual finite-mu original saddle index derivative."""
from pathlib import Path
import sys
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'experiments/root_st073'))
import lei_ren_part1_paper_compliant_current_original_Rp_pulse_saddle_index_derivative as current
import lei_ren_part1_paper_compliant_current_original_Rp_pulse_saddle_index_derivative_check as check


def main():
    owner, view = current.run()
    receipt = check.run(owner, view)
    with mp.workdps(15):
        public = current.CurrentOriginalRpPulseSaddleIndexDerivative(owner.saddle)
        accepted = public.report(public.evaluate())
        assert accepted[current.GATE] and receipt[current.GATE] and mp.mp.dps == 15
    print('ACCEPTED_ORIGINAL_RP_SADDLE_INDEX_DERIVATIVE_PUBLIC', flush=True)


if __name__ == '__main__':
    main()
