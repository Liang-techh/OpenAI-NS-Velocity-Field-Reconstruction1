"""Restore saved true cutoff own-field data or execute a new original query."""
import argparse
import json
from pathlib import Path
import sys
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_cutoff_numeric_bridge as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evaluate',action='store_true')
    parser.add_argument('--chart',choices=('core','first'),default='first')
    parser.add_argument('--coordinate',default='1',help='Core rho or first fraction q; first q=0 is core rho=4')
    parser.add_argument('--Z',default='.371')
    parser.add_argument('--cells',type=int,default=4)
    parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    raw,receipt=current.restore_checked_archive()
    if not args.evaluate:
        print(json.dumps(dict(checked=True,source_family=raw['source_family'],
            actual_original_queries=receipt['query_results'],
            own_fields=['M1','mean1','K1','Q1','V1','P1'],five_signed_inner_cumulatives=True,
            original_phase_identity='s=q*sc/4; first q=0 is rho=4',
            pressure_m13_m15_exterior_still_needed=True,
            positive_order_five_moment_repair_closed=False,temporal_recursion=False,
            full_corrected_Navier_Stokes_solution=False),indent=2))
        return
    from run_original_n1_five_defect_values import build_checked
    with mp.workdps(540):
        cache={};numeric=build_checked(current.legacy.CurrentOriginalN1RadialCutoffExtension,cache)
        hierarchy=cache[current.inner.hierarchy.CurrentOriginalCoreHierarchySource]
        full=current.full.CurrentOriginalN1FullDefectLinearFunctions()
        inner=current.inner.CurrentOriginalN1InnerNumericBridge(full,numeric.regular,hierarchy)
        owner=current.CurrentOriginalN1CutoffNumericBridge(inner,numeric)
        view=owner.report(owner.evaluate(args.Z,args.coordinate,args.chart,args.cells))
        content=(json.dumps(current.op.snapshot(view),indent=2)+'\n').encode('utf8')
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_bytes(content)
            print('TRUE_CUTOFF_OWN_FIELDS_QUERY_SAVED',args.output.resolve(),flush=True)
        else:print(content.decode('utf8'))


if __name__=='__main__':main()
