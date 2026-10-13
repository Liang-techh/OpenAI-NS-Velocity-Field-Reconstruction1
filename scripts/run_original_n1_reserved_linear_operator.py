"""Build/check the new original n1 linear operator, without ancestor runs."""
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator as current
import lei_ren_part1_paper_compliant_current_original_n1_reserved_linear_operator_check as check


def main():
    owner,raw=current.run()
    receipt=check.run(owner,raw)
    owner.assert_graph()
    if not receipt['all_passed']:raise RuntimeError('Original n1 linear repair operator rejected')
    print('Basis/operator only; full signed C1 defects and repaired fields remain open.',flush=True)


if __name__=='__main__':main()
