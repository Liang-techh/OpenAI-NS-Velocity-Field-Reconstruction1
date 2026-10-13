"""Restore actual original core functions; optionally enclose their values."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'experiments'/'root_st073'))
import lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions as current


def run():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--evaluate',action='store_true',help='Attach the checked actual nonlinear oracle and return ranges, including full radial tails.')
    parser.add_argument('--rho',default='2')
    parser.add_argument('--Z',default='.371')
    parser.add_argument('--degree',type=int,default=6)
    parser.add_argument('--depth',type=int,default=3)
    args=parser.parse_args()
    owner=current.CurrentOriginalN1CoreFixedPointFunctions();owner.assert_graph()
    result=dict(checked=owner.acceptance_loaded,pulse_checked=owner.pulse.acceptance_loaded,
        graph_nodes=len(owner.g.nodes),core_functions=current.pulse.own.encode(owner.functions),
        exact_source_parameters=current.pulse.own.encode(owner.parameters),
        actual_operator_term_counts={key:len(rows) for key,rows in owner.terms.items()},
        radial_scope='0<=rho=Lambda*R<=4.1; -1<=Z<=1',
        anchored_G_is_unique_H_root_primitive=True,
        core_epsilon_is_inverse_Lambda=True,physical_Uz_has_no_extra_Pstar=True,
        source_functions_and_directed_enclosures_are_distinct=True,
        full_inner_cutoff_lower_collar_five_defect_adapter_installed=False,
        uniform_original_core_point_numerical_accuracy_claimed=False)
    if args.evaluate:
        samples={key:owner.component_enclosure(key,args.rho,args.Z,degree=args.degree,depth=args.depth)
            for key in ('Phi','Psi','mean_Psi')}
        owner.assert_graph()
        result['same_exact_source_parameter_identity']=owner.oracle_parameter_identity
        result['directed_enclosures']=current.common.encode(samples)
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':run()
