"""Analytic first-width collar coefficients over directed axial atom jets.

Outputs are coefficients divided by h_b, with actual Z derivatives in each
slot. No whole pressure/width field is materialized. Source errors must be
supplied in the inlet atoms by the caller; these formulas do not invent them.
"""
from lei_ren_part1_paper_interval_axial_second_jet import IntervalAxialSecondJet


def first_width_coefficients(inlet, *, s=2, switch_integral='.5'):
    """Return all eight normalized first-width coefficient jets.

J=integral_0^s chi0 is an enclosure supplied by the caller. For s=2,
J=1/2 exactly. Fractional s needs a separately directed switch primitive.
"""
    u=inlet['Uz']
    if not isinstance(u,IntervalAxialSecondJet):
        raise TypeError('A directed second-Z inlet is required')
    ctx=u.ctx
    def constant(value):
        return IntervalAxialSecondJet(ctx,value,pressure_order=u.pressure_order,
                                     width_order=u.width_order)
    # Keep a supplied directed interval instead of converting through a float.
    sv=constant(s);J=constant(switch_integral)
    return dict(g=-inlet['D']*J/2,
                u=-(inlet['R']/2).sqrt()*inlet['I_z']*J,
                theta=2*sv,mz=sv*u,mixed=2*sv*u,
                axial=sv*u*u,swirl=sv,p=sv)


def encode_coefficients(coefficients):
    from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
    return {state:{slot:{str(key):dict(nominal=encode(atom.nominal),
                                      perturbation=encode(atom.difference))
                        for key,atom in getattr(jet,slot).atoms.items()}
                   for slot in ('value','tangent','second')}
            for state,jet in coefficients.items()}
