"""Independent full-integral moment/stress checks of the exterior companion."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_heat_terminal_stress_identities import terminal_stress_identities
from lei_ren_part1_paper_compliant_heat_stress_equations import heat_stress_equations
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE = Path(__file__).parent
NAME = 'lei_ren_part1_paper_compliant_heat_stress_C4.json'


def independent_terminal_integral_fixture():
    """Numerically integrate all full future moments before taking stress.

    Power substitutions remove the integrable endpoint singularities.
    Adaptive precision preserves the true H deficit at tiny arguments;
    no inverse-radius series or stress-free target defines these moments.
    """
    with mp.workdps(55):
        a, R, c = map(mp.mpf, ('.25', '20', '1.3'))
        delta, k, b = 2*a, 1-a, (1-2*a)/2

        def H(x, order=0):
            if not x:
                return (-1)**order*mp.rf(a,order)*mp.rf(1+a,order)
            return ((-1)**order*mp.rf(a,order)*mp.rf(1+a,order)
                    *x**(-1-a-order)*mp.hyperu(1+a+order,2,1/x))

        def divided_deficit(x):
            if not x:
                return a*(1+a)
            lost = max(0, int(mp.ceil(-mp.log10(x))))
            with mp.workdps(55+lost):
                return (1-H(x))/x

        cache = {}
        rows = []
        for Z in map(mp.mpf, ('.4', '-.4', '1', '-1')):
            d, L = 1-Z*Z, 1-delta*Z*Z
            xi = 2*d/R
            if xi not in cache:
                h, hp = H(xi), H(xi,1)
                if xi:
                    angular_deficit = xi/a*mp.quad(
                        lambda w:divided_deficit(xi*w**(1/a)), [0,1])
                    angular_derivative = mp.quad(
                        lambda w:H(xi*w**(1/a),1), [0,1])/a
                    def future(rate):
                        mass = mp.quad(lambda w:H(xi*w**(1/rate))**2,[0,1])/rate
                        derivative = mp.quad(
                            lambda w:H(xi*w**(1/(rate+1)))*H(xi*w**(1/(rate+1)),1),
                            [0,1])/(rate+1)
                        return mass, derivative
                    BE, JE = future(delta)
                    BP, JP = future(1+delta)
                else:
                    angular_deficit, angular_derivative = mp.mpf(0), hp/a
                    BE, JE = 1/delta, hp/(delta+1)
                    BP, JP = 1/(1+delta), hp/(delta+2)
                cache[xi] = h,hp,angular_deficit,angular_derivative,BE,JE,BP,JP
            h,hp,J,I,BE,JE,BP,JP = cache[xi]
            Mtheta = mp.sqrt(2)*c*R**(1-a)*(1/k+J)
            Mtheta_Z = mp.sqrt(2)*c*R**(1-a)*(4*Z*I/R)
            Utheta = c*R**(-mp.mpf('.5')-a)*h
            Itheta = ((1-a)*Mtheta-b*Z*Mtheta_Z-R*mp.sqrt(2*R)*Utheta)/(2*L*R)
            Stheta = -mp.sqrt(2)*c*R**(-1-a)*((1+a)*h+xi*hp)
            Menergy = c*c*R**(-delta)*BE/2
            Menergy_Z = -4*c*c*Z*R**(-delta-1)*JE
            P = -c*c*R**(-1-delta)*BP/2
            P_Z = 4*c*c*Z*R**(-delta-2)*JP
            Iz = (2*delta*Z*Menergy-d*Menergy_Z
                  +R*(2*(1+delta)*Z*P-d*P_Z))/(L*mp.sqrt(2*R))
            theta_residual, axial_residual = Itheta+Stheta, Iz
            if max(abs(theta_residual),abs(axial_residual)) > mp.mpf('1e-35'):
                raise ArithmeticError('Independent full terminal moments do not cancel stress')
            kappa = 2+delta+2*xi*hp/h
            if not 2 <= kappa <= 2+delta:
                raise ArithmeticError('Canonical positive Gamma heat shear bound failed')
            rows.append(dict(Z=str(Z), theta_residual=str(theta_residual),
                             axial_residual=str(axial_residual), kappa=str(kappa)))
        return dict(full_future_moment_quadrature_checked=True,
                    full_hyperu_function_not_a_finite_series=True,
                    independent_original_3_16_3_17_stress_checks=len(rows),
                    symmetry_and_axial_endpoints_included=True,
                    actual_project_parameters_certified_by_source_bridge_not_this_fixture=True,
                    moderate_fixture_only=True, samples=rows, passed=True)


def run():
    from lei_ren_part1_paper_compliant_heat_terminal_history_bridge import terminal_history_bridge
    record = json.loads((HERE/NAME).read_bytes())
    hashes = dict(record['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest() != digest:
            raise ValueError('Heat stress dependency changed: '+name)
    bridge = terminal_history_bridge()
    theorem = terminal_stress_identities()
    equations = heat_stress_equations()
    if bridge != record['terminal_history_bridge'] or not bridge['complete_terminal_moment_history_bridge_verified']:
        raise ValueError('Actual terminal moment source/history bridge not admitted')
    for gate in ('renormalized_energy_constant_zero_from_selected_equation',
                 'actual_exterior_energy_full_Gamma_integral_and_half_normalization_verified',
                 'actual_zero_meridional_terminal_histories_verified',
                 'physical_terminal_energy_preserved_positive'):
        if not bridge[gate]:
            raise ValueError('Actual energy/meridional source proof missing: '+gate)
    if theorem != record['terminal_stress_theorem'] or equations != record['canonical_heat_stress_equations']:
        raise ValueError('Exact heat stress equations changed')
    c = MPIntervalContext(); c.dps=270
    read = lambda q:read_interval(c,q)
    pulse_name='lei_ren_part1_paper_compliant_pulse_mixed_C4.json'
    pulse=json.loads((HERE/pulse_name).read_bytes())['whole_Z_terminal']
    actual_terminal_primitive_zero_bounds=0
    for name in ('Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared'):
        for row in pulse['primitive_y_derivative_Taylor'][name]:
            for value in row['coefficients']:
                if endpoints(read(value))!=(mp.mpf(0),mp.mpf(0)):
                    raise ArithmeticError('Actual C4 terminal primitive jet is not zero: '+name)
                actual_terminal_primitive_zero_bounds+=1
    if endpoints(read(pulse['positive_terminal_future_energy_Taylor']['coefficients'][0]))[0]<=0:
        raise ArithmeticError('Actual positive pulse terminal energy was zeroed')
    hashes[pulse_name]=hashlib.sha256((HERE/pulse_name).read_bytes()).hexdigest()
    angular_overlap_diagnostics = 0
    zero_mixed_bounds = 0
    points = record['samples']+[record['whole_Z_inlet'],record['whole_Z_unbounded_exterior']]
    for point in points:
        for key in ('all_five_terminal_moments_same_source_verified',
                    'homogeneous_constants_eliminated_from_actual_terminal_moments',
                    'heat_exterior_stress_identity_certified',
                    'original_angular_energy_and_pressure_histories_retained',
                    'exterior_energy_full_Gamma_integral_and_half_normalization_verified',
                    'renormalized_energy_constant_zero_from_selected_equation',
                    'physical_terminal_energy_preserved_positive',
                    'zero_stress_follows_from_terminal_source_history_and_full_integral_identities'):
            if not point[key]:
                raise ValueError('Actual source-bound stress gate lost: '+key)
        for key in ('global_admissible_stress_lift_constructed','whole_outer_cone_certified',
                    'physical_energy_integral_certified','temporal_recursion',
                    'actual_compliant_field_physical_heat_region_certified'):
            if point[key]:
                raise ValueError('Exterior stress companion overclaims: '+key)
        for jet in point['exterior_meridional_moment_Taylor'].values():
            if any(endpoints(read(value))!=(mp.mpf(0),mp.mpf(0)) for value in jet['coefficients']):
                raise ArithmeticError('Inherited exterior meridional history lost exact zero')
        if endpoints(read(point['energy_Taylor']['coefficients'][0]))[0]<=0:
            raise ArithmeticError('Actual Gamma exterior energy is not positive')
        for new,old in zip(point['angular_Taylor']['coefficients'],
                           point['original_forward_angular_Taylor']['coefficients']):
            nl,nh = endpoints(read(new)); ol,oh=endpoints(read(old))
            if max(nl,ol)>min(nh,oh):
                raise ArithmeticError('Equivalent angular histories have disjoint enclosures')
            angular_overlap_diagnostics += 1
        for bound in point['original_retained_angular_tail_constant_defect']['coefficients']:
            lo,hi=endpoints(read(bound))
            if not lo<=0<=hi:
                raise ArithmeticError('Retained angular constant enclosure excludes exact zero')
        kappa = read(point['theta_inertial_stress_over_F_y_derivative_Taylor'][0]['coefficients'][0])
        if endpoints(kappa)[0] <= 0:
            raise ArithmeticError('Heat angular inertial stress is not positive')
        for label,grid in point['exterior_stress_mixed_derivatives_total_order_le4'].items():
            if len(grid)!=15:
                raise ValueError('Incomplete exterior stress mixed4 grid: '+label)
            for value in grid.values():
                if endpoints(read(value)) != (mp.mpf(0),mp.mpf(0)):
                    raise ArithmeticError('Exact source-bound heat stress is not zero')
                zero_mixed_bounds += 1
    fixture = independent_terminal_integral_fixture()
    hashes[NAME] = hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result = dict(
        actual_five_defect_family_sha256=record['actual_five_defect_family_sha256'],
        implicit_source_sha256=record['implicit_source_sha256'],
        actual_terminal_moment_history_bridge_recomputed=True,
        actual_terminal_meridional_primitive_zero_bounds=actual_terminal_primitive_zero_bounds,
        actual_exterior_energy_full_integral_and_half_normalization_verified=True,
        selected_energy_equation_residual_zero_verified=True,
        physical_terminal_energy_preserved_positive=True,
        full_terminal_integral_stress_identities_recomputed=True,
        independent_terminal_integral_fixture=fixture,
        actual_forward_angular_overlap_diagnostics=angular_overlap_diagnostics,
        exact_zero_exterior_stress_mixed_bounds=zero_mixed_bounds,
        all_five_terminal_moments_same_source_verified=True,
        homogeneous_constants_eliminated_from_actual_terminal_moments=True,
        heat_exterior_stress_identity_certified=True,
        actual_compliant_field_physical_heat_region_certified=False,
        main_dispatcher_and_collar_adoption_complete=False,
        global_admissible_stress_lift_constructed=False,
        whole_outer_cone_certified=False,
        physical_energy_integral_certified=False,
        temporal_recursion=False,
        all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Source-bound full Gamma terminal moments and exterior zero stress PASS',flush=True)
    return result


if __name__ == '__main__':
    run()
