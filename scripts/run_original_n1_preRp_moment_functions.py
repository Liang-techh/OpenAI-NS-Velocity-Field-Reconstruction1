"""Restore checked signed pre-Rp C1 functions without rebuilding ancestors."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'experiments' / 'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions as functions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixtures', action='store_true', help='Run independent physical-pressure and relocation fixtures')
    args = parser.parse_args()
    if args.fixtures:
        import lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions_check as check
        check.independent_pressure_fixture()
        check.relocation_fixture()
    owner = functions.SavedCurrentOriginalN1PreRpMomentFunctions()
    owner.assert_graph()
    print(json.dumps(dict(checked=owner.acceptance_loaded,
        radial_scope=owner.raw['radial_scope'], graph_nodes=len(owner.g.nodes),
        signed_C1_functions=functions.encode(owner.functions),
        complete_five_defects_installed=False, selected_n1_coefficients_installed=False), indent=2))


if __name__ == '__main__':
    main()
