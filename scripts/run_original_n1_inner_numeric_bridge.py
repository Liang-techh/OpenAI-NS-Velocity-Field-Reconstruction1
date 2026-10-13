"""Restore saved true n1 queries; optionally execute a new enclosing query."""
import argparse
import json
from pathlib import Path
import sys

import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_inner_numeric_bridge as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evaluate',action='store_true',help='Restore numerical/source owners once and run one actual query')
    parser.add_argument('--chart',choices=('core','first'),default='core')
    parser.add_argument('--coordinate',default='2',help='Core rho in [0,4], or first fraction q in [0,1]')
    parser.add_argument('--Z',default='.371')
    parser.add_argument('--output',type=Path,help='Save the full factored/log interval query as JSON')
    args=parser.parse_args()
    raw,receipt=current.restore_checked_archive()
    if not args.evaluate:
        print(json.dumps(dict(checked=True,source_family=raw['source_family'],
            saved_true_original_queries=receipt['query_results'],
            physical_conjugation=receipt['physical_conjugation_proof'],
            exact_first_coordinate='s=q*sc/4; q in [0,1]',
            field_scope='unrepaired common inner n1 only',
            original_numeric_n1_solution_oracle_executed=True,
            n1_repaired_own_M_V_P_installed=False,
            positive_order_five_moment_repair_closed=False,
            temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))
        return
    from run_original_n1_five_defect_values import build_checked
    with mp.workdps(540):
        cache={}
        numeric=build_checked(current.regular.CurrentOriginalN1RegularSeriesFunction,cache)
        hierarchy=cache[current.hierarchy.CurrentOriginalCoreHierarchySource]
        functions=current.full.CurrentOriginalN1FullDefectLinearFunctions()
        bridge=current.CurrentOriginalN1InnerNumericBridge(functions,numeric,hierarchy)
        view=bridge.report(bridge.evaluate(args.Z,args.coordinate,args.chart))
        content=json.dumps(current.regular.op.snapshot(view),indent=2)+'\n'
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True)
            args.output.write_bytes(content.encode('utf8'))
            print('TRUE_N1_QUERY_SAVED',args.output.resolve(),flush=True)
        else:
            print(content)


if __name__=='__main__':
    main()
