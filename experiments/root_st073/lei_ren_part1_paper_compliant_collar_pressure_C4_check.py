"""Focused full-collar pressure history, integral and mixed-jet checks."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_collar_pressure_C4 import collar_pressure_bridge, collar_pressure_rows
from lei_ren_part1_paper_compliant_heat_pressure_source_bridge import source_bridge
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import CompliantCollarGammaC4
from lei_ren_part1_paper_compliant_pulse_physical_bounds import P
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'
NAME = PREFIX+'compliant_collar_pressure_C4.json'


def independent_full_collar_fixture():
    """Actual full Gamma/flat factors at moderate parameters, not NS data."""
    with mp.workdps(45):
        c = MPIntervalContext(); c.dps = 100
        a,S,eps,Z = map(mp.mpf, ('.15','.04','.002','.4'))
        rate = 1+2*a
        heat = CompliantCollarGammaC4.__new__(CompliantCollarGammaC4)
        heat.ctx=c; heat.a=c.mpf(a); heat.Scap=c.mpf(S); heat.S=c.mpf([0,S])
        heat.eps=c.mpf(eps); heat.k=1-heat.a; heat.delta=2*heat.a; heat.prate=1+heat.delta
        heat.cells=64; heat.shape_cache={}; heat.tail_cache={}; heat.gamma_cache={}

        def H(xi):
            if not xi:return mp.mpf(1)
            return xi**(-1-a)*mp.hyperu(1+a,2,1/xi)

        def sigma(v):
            if v<=0:return mp.mpf(0)
            if v>=1:return mp.mpf(1)
            return 1/(1+mp.exp(1/v**2-1/(1-v)**2))

        def K(v,z):
            phi = mp.exp(-4/(3-v)**2) if v<3 else mp.mpf(0)
            sig = sigma(v)
            return (1-sig)*(1-eps)+sig*H(2*(1-z*z)*S*mp.exp(-v))*(1-eps*phi)

        integral_checks = derivative_checks = 0
        for t in map(mp.mpf, ('.5','2','3')):
            future = heat.collar_tails(c.mpf(Z),c.mpf(t))['remaining_pressure_in_Rtail_units']
            shape = heat.shape(c.mpf(Z),c.mpf(t))
            rows = collar_pressure_rows(shape['K_rows'],future,c.mpf(rate),c.mpf(1),c.mpf(t))
            breaks = [t]+([mp.mpf(1)] if t<1 else [])+[mp.mpf(3)]
            collar = mp.quad(lambda v: mp.exp(-rate*v)*K(v,Z)**2/2,breaks) if t<3 else mp.mpf(0)
            xi3 = 2*(1-Z*Z)*S*mp.exp(-3)
            exterior = mp.exp(-3*rate)*mp.quad(lambda q:q**(rate-1)*H(xi3*q)**2,[0,'.25',1])/2
            lo,hi = endpoints(rows[0][0])
            if not lo<=-collar-exterior<=hi:
                raise ArithmeticError('Independent full collar/future pressure outside enclosure')
            integral_checks += 1
            prime = lambda y,z: mp.exp(-rate*y)*K(y,z)**2/2
            for order,n in ((1,0),(1,1),(1,2),(2,0),(2,1),(3,0),(4,0)):
                value = mp.diff(lambda z:mp.diff(lambda y:prime(y,z),t,order-1),Z,n)/math.factorial(n)
                lo,hi = endpoints(rows[order][n])
                if not lo<=value<=hi:
                    raise ArithmeticError('Independent full collar mixed derivative outside enclosure: '+str((t,order,n)))
                derivative_checks += 1
        return dict(full_sigma_phi_Gamma_future_integral_checks=integral_checks,
                    independent_mixed_derivative_checks=derivative_checks,
                    moderate_fixture_only=True, passed=True)


def run():
    record = json.loads((HERE/NAME).read_bytes())
    hashes = dict(record['input_hashes'])
    for filename,digest in hashes.items():
        if hashlib.sha256((HERE/filename).read_bytes()).hexdigest()!=digest:
            raise ValueError('Collar pressure dependency changed: '+filename)
    bridge = collar_pressure_bridge(source_bridge())
    if bridge!=record['collar_pressure_source_bridge']:
        raise ValueError('Collar pressure source bridge changed')
    c = MPIntervalContext(); c.dps = 270
    read = lambda value:read_interval(c,value)
    finite = overlaps = 0
    for point in record['samples']+[record['whole_Z_inlet'],record['whole_Z_terminal'],record['whole_Z_collar']]:
        for flag in ('original_pressure_datum_and_forward_history_retained',
                     'absolute_pressure_same_source_mixed4_available','complete_collar_pressure_history_transfer_verified',
                     'pressure_waiting_collar_and_collar_exterior_mixed4_joins_verified',
                     'full_sigma_phi_Gamma_future_integral_used'):
            if not point[flag]:raise ValueError('Missing collar pressure source/history gate: '+flag)
        if point['source_history_transfer_conditional'] or point['heat_exterior_stress_identity_certified'] or point['global_admissible_stress_lift_constructed'] or point['temporal_recursion']:
            raise ValueError('Collar pressure scope overclaimed')
        if endpoints(read(point['pressure_over_Utheta_squared_Taylor']['coefficients'][0]))[1]>=0:
            raise ArithmeticError('Full remaining pressure integral positivity lost')
        grid = point['physical_mixed_derivatives_total_order_le4'][P]
        if len(grid)!=15 or len(point['pressure_y_derivative_axial5_Taylor'])!=5:
            raise ValueError('Collar pressure mixed/axial derivative coverage incomplete')
        for key,value in grid.items():
            lo,hi=endpoints(read(value))
            if not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Nonfinite collar pressure derivative')
            oldlo,oldhi=endpoints(read(point['original_forward_pressure_mixed_bounds'][key]))
            if max(lo,oldlo)>min(hi,oldhi):raise ArithmeticError('Retained forward pressure diagnostic disjoint: '+key)
            finite+=1; overlaps+=1
    # Boundary values are the same full future integral by the bridge.
    # These independent enclosure comparisons are only diagnostics.
    oldheat = json.loads((HERE/(PREFIX+'compliant_heat_pressure_C4.json')).read_bytes())
    waiting = json.loads((HERE/(PREFIX+'compliant_steep_waiting_C4.json')).read_bytes())
    joins = 0
    for point,reference in ((record['whole_Z_inlet'],waiting['whole_Z']['waiting']['terminal']),
                            (record['whole_Z_terminal'],oldheat['whole_Z_inlet'])):
        for key,value in point['physical_mixed_derivatives_total_order_le4'][P].items():
            lo,hi=endpoints(read(value)); rl,rh=endpoints(read(reference['physical_mixed_derivatives_total_order_le4'][P][key]))
            if max(lo,rl)>min(hi,rh):raise ArithmeticError('Pressure interface diagnostic disjoint: '+key)
            joins+=1
    fixture = independent_full_collar_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],
                implicit_source_sha256=record['implicit_source_sha256'], all_passed=True,
                collar_pressure_source_bridge=bridge,
                complete_collar_absolute_pressure_history_bridge_verified=True,
                pressure_waiting_collar_and_collar_exterior_mixed4_joins_verified=True,
                absolute_pressure_same_source_mixed4_available=True,
                actual_finite_pressure_mixed_bounds_checked=finite,
                forward_pressure_overlap_diagnostics=overlaps, interface_overlap_diagnostics=joins,
                interval_overlap_used_as_functional_proof=False,
                independent_full_collar_fixture=fixture,
                global_admissible_stress_lift_constructed=False, whole_outer_cone_certified=False,
                physical_energy_integral_certified=False, temporal_recursion=False,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Full collar absolute pressure/history/flat joins and independent full-integral fixture PASS',flush=True)
    return result


if __name__=='__main__':
    run()
