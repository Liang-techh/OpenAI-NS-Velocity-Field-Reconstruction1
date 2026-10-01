"""Analytic radial-tail propagation into the fresh-family core inlet.

Bounds are relative to accepted construction data. This module encloses the
infinite core at a fixed scaled radius, not comparison ODE discretization.
Taylor coefficient k is derivative k divided by k!, including error jets.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_candidate_shared_inlet import normalized_inlet
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode


def absolute_jet_envelope(rows, radius):
    """Uniform derivative/factorial bounds on [0,radius] for finite rows."""
    ctx = rows[0].ctx
    values = []
    for k in range(rows[0].order+1):
        total = ctx.mpf(0)
        for n, row in enumerate(rows):
            lo, hi = endpoints(row[k])
            total += ctx.mpf(max(abs(lo),abs(hi)))*radius**n
        # Symmetry encloses the function and all derivative signs; products
        # of these jets bound the Leibniz convolution of derivatives.
        upper = endpoints(total)[1]
        values.append(ctx.mpf([-upper,upper]))
    return IntervalTaylor(ctx, values)


def error_jet(calc, name, scale=1):
    ctx = calc.ctx
    entries = {row['axial_order']:row[name] for row in calc.tail['rows']
               if row['scaled_radial_order'] == 0}
    values = []
    for k in range(4):
        bound = entries[k]*ctx.mpf(scale)/math.factorial(k)
        upper = endpoints(bound)[1]
        values.append(ctx.mpf([-upper,upper]))
    return IntervalTaylor(ctx, values)


def enclose(calc, radius='4'):
    c = calc.ctx
    with mp.workdps(calc.precision+40):
        s = c.mpf(radius)
        lo, hi = endpoints(s)
        if lo != hi or lo <= 0 or hi > mp.mpf('4.1'):
            raise ValueError('Fixed scaled radius in (0,4.1] required')
        finite = calc.core_state(s)
        ep = error_jet(calc, 'Phi_tail')
        eu = error_jet(calc, 'Psi_tail', calc.eps)
        bp = absolute_jet_envelope(calc.phi,s)
        bu = absolute_jet_envelope(calc.u,s)
        pp_error = bp*ep*2+ep*ep
        pu_error = bp*eu+bu*ep+ep*eu
        uu_error = bu*eu*2+eu*eu
        # Integrating a uniform axial-derivative envelope multiplies each
        # Taylor coefficient by the exact positive radial weight integral.
        full = dict(phi=finite['phi']+ep, U=finite['U']+eu,
            theta=finite['theta']+ep*s*s,
            z=finite['z']+eu*s,
            theta_z=finite['theta_z']+pu_error*s*s,
            p=finite['p']+pp_error*s,
            u_squared=finite['u_squared']+uu_error*s,
            weighted_phi_squared=finite['weighted_phi_squared']+pp_error*(s*s/2))
        pressure = calc.p0+calc.S*full['p']*calc.eps
        positive = endpoints(full['phi'][0])[0] > 0
        inlet = None
        if positive:
            zero = full['phi']*0
            # Only return D / I_z / pressure: these expressions do not use
            # radial endpoint derivatives. No stress/cone claim is made.
            transfer = normalized_inlet(calc.phi,calc.u,full,calc.S,calc.ell,calc.p0,
                calc.lam,s,calc.z,calc.delta,
                endpoints_override=dict(phi_exit=full['phi'],u_exit=full['U'],
                                        phi_s=zero,u_s=zero))
            inlet = dict(D=transfer['ratio'],I_z=transfer['iz']*c.sqrt(calc.eps),
                         pressure=transfer['pressure'])
        return dict(scaled_radius=radius, center_family=calc.center_family,
            completed_degree=calc.degree, state_sha256=calc.state_hash,
            analytic_tail_target_met=calc.tail['target_met'],
            full_phi_positive=positive, finite_core=finite,
            full_core_enclosures=full, physical_pressure=pressure,
            core_inlet=inlet, core_inlet_division_available=positive,
            Phi_error_jet=list(ep.coefficients), U_error_jet=list(eu.coefficients),
            analytic_radial_tail_propagated_into_core_moments=True,
            relative_to_accepted_construction_data=True,
            construction_parameter_errors_enclosed=False,
            comparison_ODE_error_enclosed=False, stress_cone_certified=False,
            temporal_recursion=False)


def run():
    calc = IntervalComparisonJets(4,diagnostic_partial=True)
    out = enclose(calc)
    here = Path(__file__).parent
    out['input_hashes'] = {name:hashlib.sha256((here/name).read_bytes()).hexdigest() for name in
        ('lei_ren_part1_paper_interval_core_inlet_bound.py',
         'lei_ren_part1_paper_interval_comparison_jets.py',
         'lei_ren_part1_paper_candidate_shared_inlet.py',
         'lei_ren_part1_paper_candidate_interval_core.py')}
    def pack(value):
        if isinstance(value,IntervalTaylor):
            return list(value.coefficients)
        if isinstance(value,dict):
            return {key:pack(item) for key,item in value.items()}
        return value
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(out)),indent=2)+'\n',encoding='utf-8')
    print('Core moment tail enclosure at degree',calc.degree,
          '; positive full-core inlet:',out['full_phi_positive'],flush=True)
    return out


if __name__ == '__main__':
    run()
