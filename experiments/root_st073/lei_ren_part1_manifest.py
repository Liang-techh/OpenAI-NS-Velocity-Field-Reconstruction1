"""Declared Part I reconstruction seed; no assembled field is certified."""
import json
import math
from pathlib import Path
from lei_ren_part1 import parameter_map, SOURCE


def build_manifest(h=0.001, nu=0.01):
    parameters = parameter_map(h)
    if not parameters['in_part1_stated_delta_range'] or not math.isfinite(nu) or nu <= 0:
        raise ValueError('Require 0<h<1/400 and nu>0')
    delta = parameters['delta']
    return {
        'candidate_id': 'lr1-background-seed-001',
        'source': SOURCE, 'source_version': '2609.35406v1',
        'parameters': dict(parameters, nu=float(nu), critical_time=1.0),
        'sampling': {'tau': [2.0**(-k) for k in range(3, 9)],
                     'core_R': [0.0, 0.05], 'interior_Z': [-0.5, 0.5],
                     'scope': 'Fixed interior similarity sector; outer domain unresolved'},
        'fixed_profile_expected_tau_exponents': {
            'radial_length': 0.5,
            'axial_length': (1-delta)/2,
            'axial_to_radial_aspect': -delta/2,
            'angular_and_axial_velocity': -(1+delta)/2,
            'angular_rotation_rate': -1-delta/2,
            'moving_sector_energy_dominant_components': (1-3*delta)/2,
        },
        'exponent_scope': 'Dimensional predictions at fixed (R,Z), for a nonzero fixed leading profile. Not measured dynamics, global energy, or nonlinear recursion evidence.',
        'construction_contract': {
            'radial_velocity': 'Derive from axial profile by incompressibility; no independent free fit',
            'pressure': 'Compute P0(Z) from the assembled exterior and use the same value in the core',
            'forcing': 'No residual-defined force allowed; final smooth forcing prescription remains unresolved',
            'nontriviality': 'Reject an identically zero profile and amplitude-collapse improvements; freeze a declared profile normalization before optimization',
            'layers': ['regular core', 'transition and annular moment repair', 'outer and heat exterior'],
        },
        'acceptance_stages': [
            'geometry similarity and measured exponent fits',
            'self-similar background: smooth axis, divergence, finite energy, support and matching',
            'stress-resolved background: profile-derived admissible stress and measured remainder scale laws',
            'full oscillatory correction: complete momentum max and volume L2 <=1e-3 on declared independent holdout',
        ],
        'unresolved': ['coupled delta<=d_corr*mu restriction and profile compatibility',
                       'actual compatible outer/heat profile', 'P0 and renormalized infinite-tail moments',
                       'regular nonlinear core under that P0', 'stress cone and cutoff-aware remainder',
                       'global energy and forcing', 'actual oscillatory correction'],
        'background_assembled': False, 'pde_validated': False,
        'scale_recursion_established': False,
    }


if __name__ == '__main__':
    manifest = build_manifest()
    destination = Path(__file__).with_suffix('.json')
    destination.write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf-8')
    print(destination)
