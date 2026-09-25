"""Necessary outer radius for positive-swirl repair of S and Cp moments.

The current two-moment patch ends at X=16. On a reserved pure-heat interval
[16,R], U_base=0. Any nonnegative E_new that restores both outgoing S and Cp
must satisfy a simple mass-versus-weighted-mass inequality. This is a lower
bound, not a construction or an impossibility result for other profiles.
"""
import json

import numpy as np
from scipy.optimize import brentq

from joined_field import JoinedField, ROOT
from heat_exterior import tail_moments


def run():
    field = JoinedField()
    audit = json.loads((ROOT/'remote_five_moment_audit.json').read_text())
    rows = []
    for source in audit['rows']:
        if source['variant'] != 'two_moment':
            continue
        eta = source['eta']
        delta = source['patch_end_moment_difference']
        cp_delta, s_delta = delta['Cp'], delta['S']
        cp16, s16, _ = tail_moments(16., eta, c=field.c,
                                    h=field.inner.h, n=96)
        cp16, s16 = float(cp16), float(s16)

        def available(R):
            cpR, sR, _ = tail_moments(R, eta, c=field.c,
                                       h=field.inner.h, n=96)
            cp_room = cp16-float(cpR)-cp_delta
            base_e2_mass = 2*(s16-float(sR))
            # Cp_new >= (1/(2R)) integral E_new^2, while restoring S
            # needs integral E_new^2 >= base mass + 2 delta S.
            mass_gap = 2*R*cp_room-base_e2_mass-2*s_delta
            return cp_room, base_e2_mass, mass_gap

        cp_only_min_radius = brentq(lambda R: available(R)[0], 16., 1e7)
        mass_min_radius = brentq(lambda R: available(R)[2],
                                 cp_only_min_radius*1.000001, 1e7)
        row = {'eta': eta, 'delta_S': s_delta, 'delta_Cp': cp_delta,
               'heat_Cp_tail_at_16': cp16,
               'heat_S_tail_at_16': s16,
               'Cp_only_min_R': cp_only_min_radius,
               'S_and_Cp_positivity_min_R': mass_min_radius,
               'at_R_64': dict(zip(
                   ('available_Cp_after_repair', 'base_E2_mass', 'mass_gap'),
                   available(64.))),
               'physical_radius_at_tau_0084': float(np.sqrt(
                   2*field.nu*(.0084/(1-eta**2))*mass_min_radius)),
               'physical_radius_at_tau_0128': float(np.sqrt(
                   2*field.nu*(.128/(1-eta**2))*mass_min_radius))}
        rows.append(row)
        print(json.dumps(row), flush=True)
    report = {'rows': rows,
              'derivation': 'On [16,R] the base has U=0 and E>=0. Restoring delta S>0 requires integral E_new^2 >= integral E_base^2 + 2 delta S (with any added U increasing this lower bound). Restoring delta Cp>0 requires integral E_new^2/(2X) = integral E_base^2/(2X) - delta Cp. Since X<=R and E_new^2>=0, the latter is at least integral E_new^2/(2R). Root of the resulting necessary inequality gives the listed lower radius.',
              'scope': 'Necessary positivity bound for repairing the current two-moment patch back to the unchanged heat exterior using corrections confined to [16,R]. It assumes the registered heat tail and unchanged inner patch. Other inner/exterior profiles or pressure data are outside this bound. Not a five-moment solution or full PDE acceptance.',
              'accepted': False}
    (ROOT/'five_moment_repair_bound.json').write_bytes(
        (json.dumps(report, indent=2)+'\n').encode())


if __name__ == '__main__':
    run()
