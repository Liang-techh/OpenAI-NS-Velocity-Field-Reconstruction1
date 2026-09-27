"""Evaluate fixed-mean wave growth using the saved higher-grid operators."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def decode(values):
    values = np.asarray(values,float)
    return values[...,0]+1j*values[...,1]


def run(candidate_path,output_path):
    candidate_raw = Path(candidate_path).read_bytes()
    candidate = json.loads(candidate_raw)
    operator_raw = (ROOT/'wave_growth_consistency.json').read_bytes()
    operator = json.loads(operator_raw)['operator_matrices']
    geometry = json.loads((ROOT/'wave_stress_growth_codesign.json').read_text())
    c = decode(candidate['selected']['coefficients_original'])
    mass,stiffness,strain = [decode(operator[k]) for k in ('mass','stiffness','strain')]
    quadratic = lambda matrix:float(np.vdot(c,matrix@c).real)
    energy = .25*quadratic(mass)
    production = -.5*quadratic(strain)
    dissipation = .5*geometry['viscosity']*quadratic(stiffness)
    result = dict(accepted=False,pde_validated=False,scale_recursion_established=False,
        candidate=Path(candidate_path).name,candidate_sha256=hashlib.sha256(candidate_raw).hexdigest(),
        operator_source='wave_growth_consistency.json',operator_sha256=hashlib.sha256(operator_raw).hexdigest(),
        quadrature=operator['label'],point_count=operator['point_count'],
        coefficients_original=candidate['selected']['coefficients_original'],
        physical_energy=energy,shear_production=production,viscous_dissipation=dissipation,
        growth_lambda=(production-dissipation)/(2*energy),
        scope='Fixed geometry and original mean velocity; only wave coefficients changed. Higher finite quadrature, not continuum convergence or wave integration.')
    Path(output_path).write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:result[k] for k in ('physical_energy','shear_production','viscous_dissipation','growth_lambda')}))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('candidate',type=Path)
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    run(args.candidate,args.output)
