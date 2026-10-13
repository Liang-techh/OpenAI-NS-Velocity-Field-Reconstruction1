"""Select/check refined original pulse coefficients without ancestor reruns."""
from pathlib import Path
import sys
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'experiments/root_st073'))
import lei_ren_part1_paper_compliant_current_original_Rp_pulse_coefficients_refined as current
import lei_ren_part1_paper_compliant_current_original_Rp_pulse_coefficients_refined_check as check


def main():
    owner, packets = current.run()
    receipt = check.run(owner, packets)
    with mp.workdps(15):
        public = current.CurrentOriginalRpPulseCoefficientsRefined(owner.selected, owner.basis)
        view = public.report(public.evaluate('.371'))
        assert view[current.GATE] and receipt[current.GATE]
        assert mp.mp.dps == 15
    print('ACCEPTED_ORIGINAL_PULSE_REFINED_COEFFICIENTS_PUBLIC', flush=True)


if __name__ == '__main__':
    main()
