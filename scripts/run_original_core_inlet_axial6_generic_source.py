"""Inspect the sealed core-inlet C6/generic C5 source without rebuilding."""
import json
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source as current


def main():
    raw=current.parent.linear.read(current.NAME);receipt=current.parent.linear.read(current.RECEIPT)
    pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(current.NAME,None)
    checker_name=Path(current.__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
    if (not receipt.get('all_passed') or not all(receipt.get(k) for k in current.GATES)
        or any(receipt.get(k) for k in current.OPEN) or pins!=raw['input_hashes']
        or archive!=current.sha(current.NAME) or checker!=current.sha(checker_name)
        or raw['source_family']!=receipt['source_family']):
        raise ValueError('Checked actual core-inlet source archive required')
    for name,digest in receipt['input_hashes'].items():
        if current.sha(name)!=digest:raise ValueError('Changed core-inlet source '+name)
    if not all(check['passed'] for check in receipt['checks'].values()):
        raise ValueError('Scoped source/unit/calculus/numerical checks required')
    views=raw['numeric_core_inlet_views']
    for view in views.values():
        if (not all(len(v)==7 for v in view['actual_complete_seven_Taylor_sources'].values())
            or len(view['actual_radial_Q_axial5'])!=6 or len(view['actual_full_signed_p2_axial5'])!=6
            or view['all_17_chart_source_installed']):
            raise ValueError('Genuine inlet-only seven-source/six-output rows required')
    print(json.dumps(dict(checked=True,radial_scope='rho=4 core inlet',
        actual_complete_C6_source_functions=True,actual_Q_inertial_p2_C5_functions=True,
        actual_sealed_numeric_views=list(views),P0_higher_derivatives_executable=True,
        all17_chart_axial6_source=False,compact_C5_targets_controls=False,
        complete_five_defect_C3=False,global_n1=False,temporal_recursion=False,
        full_corrected_Navier_Stokes_solution=False),indent=2))


if __name__=='__main__':main()
