"""Restore true lower-pressure C1 enclosures, or execute an original query."""
import argparse
import json
from pathlib import Path
import sys
import mpmath as mp

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_lower_numeric_bridge as current


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evaluate',action='store_true')
    parser.add_argument('--Z',default='.371')
    parser.add_argument('--include-inner',action='store_true',help='Also compose a new actual Rin packet through rminus')
    parser.add_argument('--output',type=Path)
    args=parser.parse_args();raw,receipt=current.restore_checked_archive()
    if not args.evaluate:
        print(json.dumps(dict(checked=True,source_family=raw['source_family'],
            actual_original_queries=[dict(name=q['name'],Z_bounds=[q['Z']['lower'],q['Z']['upper']],
                true_signed_pressure_value_Z1_rows=q['true_signed_pressure_value_Z1_rows'],
                composed_with_actual_Rin_history=q['composed_with_actual_Rin_history'])
                for q in receipt['query_results']],
            signed_pressure_C1_functions=['m13','m15'],radial_scope='Rin=sc/4 to rminus=sc/2',
            complete_external_pressure_transfer=False,positive_order_five_moment_repair_closed=False,
            temporal_recursion=False,full_corrected_Navier_Stokes_solution=False),indent=2))
        return
    from run_original_n1_five_defect_values import build_checked
    with mp.workdps(540):
        cache={};numeric=build_checked(current.cut.legacy.CurrentOriginalN1RadialCutoffExtension,cache)
        hierarchy=cache[current.cut.inner.hierarchy.CurrentOriginalCoreHierarchySource]
        full=current.full.CurrentOriginalN1FullDefectLinearFunctions()
        inner=current.cut.inner.CurrentOriginalN1InnerNumericBridge(full,numeric.regular,hierarchy)
        cut=current.cut.CurrentOriginalN1CutoffNumericBridge(inner,numeric)
        prefix=current.prior.CurrentOriginalN1UnchangedPressurePrefix(numeric.regular.mean)
        owner=current.CurrentOriginalN1LowerNumericBridge(cut,prefix)
        Rin=cut.evaluate(args.Z,'1','first') if args.include_inner else None
        view=owner.report(owner.evaluate(args.Z,Rin_packet=Rin))
        content=(json.dumps(current.op.snapshot(view),indent=2)+'\n').encode('utf8')
        if args.output:
            args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_bytes(content)
            print('TRUE_LOWER_PRESSURE_C1_QUERY_SAVED',args.output.resolve(),flush=True)
        else:print(content.decode('utf8'))


if __name__=='__main__':main()
