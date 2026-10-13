"""Restore checked original fields and produce/check finite endpoint C1 data."""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'experiments' / 'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1 as current
import lei_ren_part1_paper_compliant_current_original_n1_preRc_endpoint_C1_check as check


def main():
    owner, packets = current.run()
    receipt = check.run(owner)
    if not receipt['all_passed'] or not receipt[current.GATE]:
        raise RuntimeError('Original finite endpoint C1 source rejected')
    accepted = current.CurrentOriginalN1PreRcEndpointC1(owner.fields)
    accepted.report(accepted.evaluate(('0', '.5')))
    print('Checked finite endpoint C1 axial cells:', len(packets), flush=True)
    print('Full pressure derivatives and actual n1 repaired fields remain open.', flush=True)


if __name__ == '__main__':
    main()
