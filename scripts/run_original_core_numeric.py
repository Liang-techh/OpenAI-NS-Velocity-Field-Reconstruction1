"""Compute original normalized nonlinear-core approximations and true errors."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
import mpmath as mp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from openai_ns_reconstruction.original_core_numeric import OriginalCoreNumeric, _box, _bounds

DIRECTORY = ROOT/"experiments/root_st073"
NAME = "lei_ren_part1_paper_compliant_current_original_core_numeric_approximation.json"


def solve(precision,nodes,degree):
    print(f"Original numerical core: precision={precision}, quadrature panels={nodes}, degree={degree}",flush=True)
    solver = OriginalCoreNumeric(precision=precision,quadrature_nodes=nodes)
    c = solver.ctx
    early = _box(solver.iv,solver.pressure["stages"]["slope_transition_ref"]["mass"])
    mass_text = c.nstr(solver.mass_point,precision)
    mass_error = max(abs(v) for v in _bounds(c,early-solver.iv.mpf(mass_text)))
    results = []
    for Z in ("-0.3","0","0.3"):
        packet = solver.coefficients(Z,degree=degree)
        for rho in ("0","1","2","4"):
            results.append(solver.evaluate(packet,rho))
        print(f"Completed actual coupled prefix at Z={Z}",flush=True)
    return solver,dict(precision=precision,quadrature_panels=nodes,radial_degree=degree,
        nominal_early_pressure_atom_proxy=mass_text,
        early_proxy_uses_numerical_cutoff_primitive_interpolation=True,
        early_proxy_quadrature_error_independently_certified=False,
        early_proxy_error_comes_from_accepted_original_interval=True,
        early_pressure_atom_directed_absolute_error_upper=c.nstr(mass_error*(1+c.mpf('1e-30')),40),
        early_pressure_atom_error_exact_mpf_tuple=list(mass_error._mpf_),
        exact_original_j_fraction=str(solver.j_fraction),samples=results)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-report",action="store_true")
    args=parser.parse_args()
    started=time.monotonic()
    coarse,first=solve(50,96,12)
    fine,second=solve(70,192,18)
    # Compare serialized decimals above both working precisions. Otherwise
    # cancellation in a difference can exceed nearly saturated error bars.
    c=mp.mp.clone()
    c.dps=120
    comparisons=[]
    nonzero=0
    max_relative={field:c.mpf(0) for field in ('Phi','Psi','mean_Psi','Uz','Q','P0_over_Pstar_squared')}
    for old,new in zip(first['samples'],second['samples']):
        assert (old['Z'],old['rho'])==(new['Z'],new['rho'])
        for field,view in new['profiles'].items():
            point=c.mpf(view['approximation'])
            if point:
                nonzero+=1
                relative=c.mpf(view['relative_error_upper'])
                assert relative<c.mpf('1e-12'),(new['Z'],new['rho'],field,relative)
                max_relative[field]=max(max_relative[field],relative)
            oldview=old['profiles'][field]
            difference=abs(point-c.mpf(oldview['approximation']))
            joint=c.mpf(view['absolute_error_upper'])+c.mpf(oldview['absolute_error_upper'])
            assert difference<=joint,(new['Z'],new['rho'],field)
            comparisons.append(dict(Z=new['Z'],rho=new['rho'],field=field,
                difference=c.nstr(difference,35),combined_error_upper=c.nstr(joint,35),consistent=True))
    # Independent special-function evaluation checks the numeric angular
    # prefix against its analytic leading model and the actual error budget.
    bessel_checks=0
    for sample in second['samples']:
        z,r=c.mpf(sample['Z']),c.mpf(sample['rho'])
        H=(1-fine.delta)*z/2+(1-z*z)*(4*z+fine.j)
        chi=H*H/(H*H+(fine.j/500)**2)
        q=chi*r/2
        model=c.mpf(1) if not q else c.besselj(1,2*c.sqrt(q))/c.sqrt(q)
        delta=abs(c.mpf(sample['profiles']['Phi']['approximation'])-model)
        assert delta<=c.mpf(sample['profiles']['Phi']['absolute_error_upper'])
        bessel_checks+=1
    mass_difference=abs(c.mpf(first['nominal_early_pressure_atom_proxy'])-c.mpf(second['nominal_early_pressure_atom_proxy']))
    assert mass_difference<=c.mpf(first['early_pressure_atom_directed_absolute_error_upper'])+c.mpf(second['early_pressure_atom_directed_absolute_error_upper'])
    paths=[Path(__file__),ROOT/'src/openai_ns_reconstruction/original_core_numeric.py']
    report=dict(source_family=fine.family,levels=[first,second],checks=dict(passed=True,
        nonzero_normalized_profile_approximations=nonzero,finite_prefix_error_comparisons=len(comparisons),
        independent_Bessel_comparisons=bessel_checks,
        maximum_refined_relative_error_upper={key:c.nstr(value,30) for key,value in max_relative.items()},
        nominal_pressure_proxy_refinement_error_bars_consistent=True,
        no_interval_midpoint_selected=True,exact_rational_j_rederived=True,
        same_current_n1_source_family=True),
        comparisons=comparisons,input_hashes=fine.hashes,
        implementation_hashes={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
        limitations=dict(original_normalized_leading_core_numerical_approximation=True,
            point_proxy_swirl_source_zero_with_actual_nonzero_source_retained_in_error_branch=True,
            all14_pressure_atoms_retained_in_error_branch=True,
            exact_full_original_source_point_oracle=False,original_F0_evaluated=False,
            full_physical_velocity_evaluated=False,global_matching_completed=False,
            temporal_recursion=False,full_corrected_NS_solution=False),
        execution_seconds=round(time.monotonic()-started,3))
    if args.write_report:
        (DIRECTORY/NAME).write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8',newline='\n')
    print(json.dumps(dict(checks=report['checks'],execution_seconds=report['execution_seconds'],
        refined_nominal_pressure_atom=second['nominal_early_pressure_atom_proxy'],
        refined_pressure_atom_absolute_error_upper=second['early_pressure_atom_directed_absolute_error_upper']),indent=2),flush=True)


if __name__=='__main__':
    main()
