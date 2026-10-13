"""Load checked ancestors and compute only continuous axial prefix models."""
from pathlib import Path
import argparse
import gzip
import json
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments/root_st073'))
from run_n1_unchanged_pressure_prefix import build_checked_owner
import lei_ren_part1_paper_compliant_current_original_n1_unchanged_axial_pressure_models as current
import lei_ren_part1_paper_compliant_current_original_n1_unchanged_axial_pressure_models_check as check


@current.source_precision
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-existing',action='store_true',help='Check saved models without recomputing them')
    args=parser.parse_args()
    prefix=build_checked_owner(current.prior.CurrentOriginalN1UnchangedPressurePrefix)
    if args.check_existing:
        owner=current.CurrentOriginalN1UnchangedAxialPressureModels(prefix,False)
        data=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    else:
        owner,data=current.run(prefix)
    result=check.run(owner,data)
    public=current.CurrentOriginalN1UnchangedAxialPressureModels(prefix)
    for query in ('-1','0','.371','1',['.211','.371'],['-1','1']):
        assert public.evaluate(query)[current.GATE]
    assert result[current.GATE]
    print('ACCEPTED_UNCHANGED_AXIAL_PRESSURE_MODELS_PUBLIC',public.acceptance_loaded,flush=True)


if __name__=='__main__':main()
