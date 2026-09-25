"""Attach the implicit-jet extended core to a wider axial collar."""
from dataclasses import asdict
import json

import numpy as np

from axial_compact_join import AxiallyCompactField
from extended_axis_jet import ExtendedFullRadialField, ExtendedParameters
from full_radial import FullRadialField
from heat_exterior import physical
from heated_interior_join import HeatedInteriorJoin
from joined_field import JoinedField, independent_fd
from radial_continuation import ROOT
from wide_pressure_fit import PressureBubbleField


def load_extended_compact_candidate(eta_max=.97, eta_flat=.7,
                                    eta_outer=.94,
                                    swirl_bubble_amplitude=0.,
                                    outer_swirl_bubble_amplitude=0.):
    source = FullRadialField.load(
        ROOT/'NS_ST073_Full_Local_Recurrence/data/ST073-V-wide14.json')
    params = asdict(source.p)
    params['eta_max'] = eta_max
    inner = ExtendedFullRadialField(ExtendedParameters(**params))
    fit = json.loads((ROOT/'wide_pressure_fit.json').read_text())
    radial = JoinedField(inner=inner, join_X=fit['join_X'],
                         outer_ratio=fit['outer_ratio'],
                         swirl_bubble_amplitude=swirl_bubble_amplitude,
                         outer_swirl_bubble_amplitude=outer_swirl_bubble_amplitude)
    pressure = PressureBubbleField(radial, fit['pressure_amplitudes'])
    return AxiallyCompactField(pressure, eta_flat=eta_flat,
                               eta_outer=eta_outer)


def load_extended_heated_candidate(swirl_bubble_amplitude=0.,
                                   outer_swirl_bubble_amplitude=0.):
    return HeatedInteriorJoin(compact=load_extended_compact_candidate(
        swirl_bubble_amplitude=swirl_bubble_amplitude,
        outer_swirl_bubble_amplitude=outer_swirl_bubble_amplitude),
                              heat_age=2048.)


def run():
    old = HeatedInteriorJoin()
    new = load_extended_heated_candidate()
    k = 5.5
    tau = .5*2**(-k)
    X = np.full(6, .03)
    eta = np.array([.3, .48, .65, .75, .85, .9])
    points = new.compact.joined.inner.from_similarity(X, eta, tau)
    rows = []
    for i, (point, e) in enumerate(zip(points, eta)):
        R, d = independent_fd(new, point[None, :], tau,
                              .001*np.sqrt(new.nu*tau), .00025*tau)
        row = dict(eta=float(e), X=float(X[i]),
                   new_residual_norm=float(np.linalg.norm(R[0])),
                   new_divergence=float(d[0]))
        if e < old.compact.eta_outer:
            old_R, _ = independent_fd(old, point[None, :], tau,
                                      .001*np.sqrt(old.nu*tau), .00025*tau)
            row['old_residual_norm'] = float(np.linalg.norm(old_R[0]))
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = dict(k=k, tau=tau, rows=rows,
                  eta_flat=new.compact.eta_flat,
                  eta_outer=new.compact.eta_outer,
                  scope='Inner X=.03 sample comparison and pure-exterior analytic-form bound; radial bridge and full-volume gates not evaluated.',
                  accepted=False)
    tau_min = .5/64
    r_edge = float(new.attachment_radius(tau_min))
    edge_p = float(physical(np.array([[r_edge, 0., 0.]]), tau_min,
                            c=new.heat.amplitude, h=new.heat.h,
                            nu=new.nu)['pressure'][0])
    age = new.heat.heat_age+.5-tau_min
    report['pure_exterior_spatial_time_bound'] = float(
        abs(edge_p)*new.heat.heat_age/age*np.exp(-.5)
        /np.sqrt(new.nu*age))
    report['axial_energy_factor_at_tau_min'] = new.heat.axial_energy_factor(tau_min)
    report['outer_attachment_radius_at_tau_min'] = r_edge
    qmax = tau_min/(1-new.compact.eta_outer**2)
    report['compact_axial_height_at_tau_min'] = float(
        np.sqrt(new.nu)*qmax**(.5-new.heat.h)*new.compact.eta_outer)
    (ROOT/'extended_compact_join.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
