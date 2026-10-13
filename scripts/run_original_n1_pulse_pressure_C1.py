"""Restore one checked pulse owner and build/check its C1 pressure adapter."""
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import run_original_n1_five_defect_values as restore
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_source as pulse
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C1 as current
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C1_check as check


def main():
    pulse_owner=restore.build_checked(pulse.CurrentOriginalN1PulsePressureSource)
    owner,_,models=current.run(pulse_owner)
    receipt=check.run(owner)
    if not receipt['all_passed'] or not receipt[current.LOCAL_GATE]:
        raise RuntimeError('Original pulse C1 pressure source or local accuracy rejected')
    accepted=current.CurrentOriginalN1PulsePressureC1Source(pulse_owner)
    accepted.report(accepted.evaluate('.371'))
    print('Checked Rp-to-Rv C0/C1 grid sizes:',len(models),len(accepted.c1_grid),flush=True)
    print('Full pressure derivative accuracy and n1 repaired fields remain open.',flush=True)


if __name__=='__main__':main()
