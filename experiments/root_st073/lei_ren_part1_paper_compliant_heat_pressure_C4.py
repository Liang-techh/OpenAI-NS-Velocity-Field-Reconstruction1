"""Conditional absolute-pressure candidate on the Gamma exterior, axial5/mixed4.

The admitted five-moment closure identifies P0+Mp with the negative
remaining pressure integral.  Use that equivalent expression before
interval enclosure, retaining the original forward history for comparison.
Source/history transfer remains conditional; default construction rejects it.
No pressure datum, velocity, or angular correction coefficient is changed.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import CompliantCollarGammaC4, product_rows
from lei_ren_part1_paper_compliant_heat_pressure_source_bridge import source_bridge
from lei_ren_part1_paper_compliant_pulse_physical_bounds import P
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def pressure_y_rows(K_rows, pressure_numerator, rate, scale, offset):
    """Ordinary y derivatives, each containing ordinary axial Taylor jets.

    P/scale = -exp(-rate*y)*A_p(y), and
    d_y(P/scale) = exp(-rate*y)*H(y)^2/2.
    A_p is the full infinite Gamma integral, not a finite inverse-radius
    series. Derivatives one through four follow from the defining pressure
    integral and the exact Leibniz rule.
    """
    c = pressure_numerator.ctx
    factor = scale*c.exp(-rate*offset)
    rows = [-pressure_numerator*factor]
    square = product_rows(K_rows, K_rows)
    for order in range(1, 5):
        row = K_rows[0]*0
        for j in range(order):
            row += square[j]*(math.comb(order-1, j)*(-rate)**(order-1-j))
        rows.append(row*(factor/2))
    return rows


class CompliantHeatPressureC4:
    """Callable exterior pressure view attached to the accepted C4 field."""
    def __init__(self, cells=64, *, allow_conditional=False):
        self.heat = CompliantCollarGammaC4(cells)
        self.ctx = c = self.heat.ctx
        self.family, self.source = self.heat.family, self.heat.source
        self.hashes = dict(self.heat.hashes)
        requirements = {
            'compliant_collar_Gamma_C4_check': ('all_passed', 'waiting_collar_and_collar_Gamma_joins_certified'),
            'compliant_absolute_moment_closure_check': ('all_passed', 'all_five_terminal_moment_identities_certified'),
            'compliant_absolute_moment_closure': ('full_global_pressure_terminal_identity_verified',
                                                  'original_forward_Mp_and_P0_retained')}
        for stem, gates in requirements.items():
            name = PREFIX+stem+'.json'
            record = json.loads((HERE/name).read_bytes())
            if any(not record.get(gate) for gate in gates):
                raise ValueError('Heat pressure prerequisite not admitted: '+name)
            if record['actual_five_defect_family_sha256'] != self.family:
                raise ValueError('Heat pressure family mismatch: '+name)
            source = record.get('implicit_source_sha256', record.get('source_sha256'))
            if source != self.source:
                raise ValueError('Heat pressure datum mismatch: '+name)
            for filename, digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/filename).read_bytes()).hexdigest() != digest:
                    raise ValueError('Heat pressure source changed: '+filename)
            self.hashes.update(record['input_hashes'])
            self.hashes[name] = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        self.bridge = source_bridge()
        if not self.bridge['complete_defining_function_history_bridge_verified'] and not allow_conditional:
            raise ValueError('Pressure source/history bridge incomplete; candidate use requires allow_conditional=True')
        bridge_name = PREFIX+'compliant_heat_pressure_source_bridge.py'
        self.hashes[bridge_name] = hashlib.sha256((HERE/bridge_name).read_bytes()).hexdigest()

    def exterior(self, Z, t):
        c = self.ctx
        Z, t = c.mpf(Z), c.mpf(t)
        if endpoints(t)[0] < 3:
            raise ValueError('Exact Gamma exterior offset t>=3 required')
        point = dict(self.heat.exterior(Z, t))
        local = point['full_local_Gamma_future_source']
        rows = pressure_y_rows(local['K_rows'], local['pressure_numerator'],
                               self.heat.prate, self.heat.pressure_scale, t)
        # Copy the two nested containers so the original provider packet is
        # still independently accessible with its original pressure history.
        mixed = dict(point['physical_mixed_derivatives_total_order_le4'])
        fields = dict(point['physical_velocity_and_pressure_y_derivative_Taylor'])
        point['original_forward_pressure_over_Pstar_squared_Taylor'] = point['pressure_over_Pstar_squared_Taylor']
        point['original_forward_pressure_mixed_bounds'] = mixed[P]
        mixed[P] = {'y'+str(k)+'_Z'+str(n): jet[n]*math.factorial(n)
                    for k, jet in enumerate(rows) for n in range(5-k)}
        fields[P] = [jet.truncate(4-k) for k, jet in enumerate(rows)]
        point.update(physical_mixed_derivatives_total_order_le4=mixed,
                     physical_velocity_and_pressure_y_derivative_Taylor=fields,
                     pressure_over_Pstar_squared_Taylor=rows[0],
                     pressure_y_derivative_axial5_Taylor=rows,
                     pressure_over_Utheta_squared_Taylor=-local['pressure_numerator']/(local['K_rows'][0]**2),
                     exact_pressure_definition='Candidate P=-integral_R^infinity Utheta^2/(2rho) drho; equality with retained P0+Mp requires the unresolved source/history bridge',
                     pressure_infinity_offset_exactly_zero_from_source_closure=False,
                     original_pressure_datum_and_forward_history_retained=True,
                     defining_source_and_pressure_scale_bridge_verified=False,
                     absolute_pressure_same_source_mixed4_available=False,
                     candidate_pressure_tail_mixed4_available=True,
                     source_history_transfer_conditional=True,
                     exact_Gamma_exterior_pressure_only=True,
                     heat_exterior_stress_identity_certified=False,
                     global_admissible_stress_lift_constructed=False,
                     temporal_recursion=False)
        return point

    def report(self):
        def pressure_summary(point):
            keys = ('Z', 'coordinate', 'pressure_over_Pstar_squared_Taylor',
                    'pressure_over_Utheta_squared_Taylor',
                    'original_forward_pressure_over_Pstar_squared_Taylor',
                    'original_forward_pressure_mixed_bounds',
                    'exact_positive_Ev0Pstar2_log_parts', 'exact_relative_velocity_log_parts',
                    'exact_pressure_definition', 'pressure_infinity_offset_exactly_zero_from_source_closure',
                    'original_pressure_datum_and_forward_history_retained',
                    'defining_source_and_pressure_scale_bridge_verified',
                    'candidate_pressure_tail_mixed4_available', 'source_history_transfer_conditional',
                    'absolute_pressure_same_source_mixed4_available', 'exact_Gamma_exterior_pressure_only',
                    'heat_exterior_stress_identity_certified', 'global_admissible_stress_lift_constructed',
                    'temporal_recursion')
            summary = {key: point[key] for key in keys}
            summary['physical_mixed_derivatives_total_order_le4'] = {P: point['physical_mixed_derivatives_total_order_le4'][P]}
            return summary
        with mp.workdps(270):
            return dict(actual_five_defect_family_sha256=self.family,
                        implicit_source_sha256=self.source,
                        scope='Conditional pressure-tail candidate: Gamma exterior t=log(R/Rtail)>=3, Z in[-1,1]; axial5 and logR/Z mixed4',
                        samples=[pressure_summary(self.exterior(z, t)) for z, t in (('0', '3'), ('.5', '3'), ('.5', '4'))],
                        whole_Z_inlet=pressure_summary(self.exterior([-1, 1], 3)),
                        whole_Z_unbounded_exterior=pressure_summary(self.exterior([-1, 1], [3, mp.inf])),
                        same_source_absolute_pressure_closure_used=True,
                        defining_source_bridge=self.bridge,
                        pressure_datum_or_velocity_changed=False,
                        exact_pressure_tail_not_a_fitted_correction=True,
                        absolute_pressure_same_source_mixed4_available=False,
                        candidate_pressure_tail_mixed4_available=True,
                        source_history_transfer_conditional=True,
                        heat_exterior_stress_identity_certified=False,
                        global_admissible_stress_lift_constructed=False,
                        temporal_recursion=False, input_hashes=self.hashes)


def run():
    result = CompliantHeatPressureC4(allow_conditional=True).report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)), indent=2)+'\n', encoding='utf8')
    print('Conditional Gamma pressure candidate generated; source/history transfer remains open', flush=True)
    return result


if __name__ == '__main__':
    run()
