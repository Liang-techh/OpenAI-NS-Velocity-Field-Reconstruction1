"""Refine and check original pulse basis/energy inputs without ancestor reruns."""
from pathlib import Path
import sys
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'experiments/root_st073'))
import lei_ren_part1_paper_compliant_current_original_Rp_pulse_basis_energy as current
import lei_ren_part1_paper_compliant_current_original_Rp_pulse_basis_energy_check as check


def main():
    owner, view = current.run()
    receipt = check.run(owner, view)
    with mp.workdps(15):
        public = current.CurrentOriginalRpPulseBasisEnergy()
        query = public.report(public.evaluate())
        assert query[current.GATE] and receipt[current.GATE]
        assert mp.mp.dps == 15
    print('ACCEPTED_ORIGINAL_PULSE_BASIS_ENERGY_PUBLIC', flush=True)


if __name__ == '__main__':
    main()
