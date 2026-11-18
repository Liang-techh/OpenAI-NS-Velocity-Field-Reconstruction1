"""Independent whole-axis bound for the zero-swirl finite inlet of ST073."""
from __future__ import annotations

import ast
from math import factorial, prod

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext

from .original_collar_whole_z import _dyadic_inv_p2, hull, read, require, sign


MAX_ORDER = 24
MASS_BOUND = 2**40
TUBE_GAP = 200  # The loss per radial step is 1/200.
SEED_FACTOR = 2**10
HIGHER_FACTOR = 2**14


def _definition(text, name):
    return next(n for n in ast.parse(text).body if isinstance(n,ast.FunctionDef) and n.name==name)


def _assigned(fn, name, expression):
    expected = ast.dump(ast.parse(expression, mode='eval').body)
    return any(isinstance(node,ast.Assign)
               and any(isinstance(t,ast.Name) and t.id==name for t in node.targets)
               and ast.dump(node.value)==expected for node in ast.walk(fn))


def _source_guards(inputs):
    older = inputs['finite']
    check = older['functional_difference']
    require(older['passed'] and check['ast_checks']['zero_nominal_swirl_seed']
            and check['ast_checks']['only_nominal_beta2_pressure']
            and check['ast_checks']['radius_four_derivative_inputs']
            and check['ast_checks']['same_mean_derivative_recipes'],
            'finite zero-swirl and radius-four provenance required')
    source = _definition(inputs['step_code'],'advance_scaled_one')
    initial = ast.unparse(_definition(inputs['initial_code'],'initial_rows'))
    q_jets = ast.unparse(_definition(inputs['pressure_code'],'q_jets'))
    tests = {
        'initial_rows': "'Uz': [_coerce_vector(ctx, fixed['U0_Z_taylor'])]" in initial
                        and "'P': [_coerce_vector(ctx, fixed['P0_Z_taylor'])]" in initial,
        'pressure_poles': 'a = 1 + center ** 2' in q_jets and 'b = 2 * center' in q_jets
                          and 'coefficients = [1 / a ** 2]' in q_jets
                          and 'coefficients.append(-(b * (n + 2) * current + (n + 3) * previous) / (a * (n + 1)))' in q_jets,
        'mean_derivative_W': _assigned(source,'wi',
            'scale(add(scale(mul(z,u[i]),1-delta_iv),mul(d,diff(u[i]))),-convert(1)/(i+1))')
            and 'wi = add(const(1), wi)' in ast.unparse(source),
        'axial_weight_H': _assigned(source,'hi','mul(d,u[i])')
            and 'hi = add(hi, scale(z, (1 - delta_iv) / 2))' in ast.unparse(source),
        'ux_products': _assigned(source,'rhs_u',
            'add(rhs_u,scale(mul(W[i],u[n-i]),eps*(n-i)),scale(mul(H[i],diff(u[n-i])),eps))')
            and _assigned(source,'uu','add(uu,mul(u[i],u[n-i]))'),
        'pressure_and_quadratic': _assigned(source,'rhs_u',
              'add(rhs_u,scale(add(u[n],scale(mul(z,uu),-2)),eps*(1+delta_iv)/2),'
              'mul(d,diff(pressure[n])),scale(mul(z,pressure[n]),-2*(1+delta_iv)))'),
        'only_angular_input_to_U_is_swirl': _assigned(source,'rhs_u',
            'add(rhs_u,scale(mul(z,mul(squared,aa_previous)),-2))'),
        'higher_pressure_exact_zero': _assigned(source,'p_next','scale(mul(squared,aa_current),convert(1)/(n+1))'),
        'axial_denominator': _assigned(source,'u_next','scale(mul(invL,rhs_u),convert(1)/(2*(n+1)**2))'),
        'ordinary_taylor_axial_derivative': '[(k + 1) * values[k + 1] if k + 1 < len(values) else zero for k in range(K)]'
                                            in ast.unparse(source),
    }
    require(all(tests.values()), 'finite axial recurrence changed: '+', '.join(k for k,v in tests.items() if not v))
    return tests


def _bounds(iv, q2, core):
    logs = q2['actual_correlated_q2_inlet']
    require(logs['source_log_basis_names'] == ['hb','F0_anchor','Pstar_squared','Lambda','sigma_a','sigma_2a'],
            'unrecognized physical log basis')
    hb_log, _, p2_log, lambda_log, _, _ = (read(iv,r) for r in logs['source_log_bases'])
    delta = read(iv,q2['actual_selected_leading_pressure_primitive']['delta_enclosure'])
    j = read(iv,core['required_j'])
    require(0 < delta.a and delta.b < iv.mpf(1)/8 and 0 < j.a and j.b < 1
            and p2_log.a > 0 and hb_log.b < -1024,
            'complex tube or logarithmic geometry not admitted')

    # U1/(epsilon Pstar^2) <= 521+352 |Mf| on the radius-1/4 capsule,
    # using |z|<=2, |d|<=5, |L|>=1/2, |1+z^2|>=1/2.
    require(521+352*MASS_BOUND <= SEED_FACTOR*(1+MASS_BOUND), 'seed inequality failed')
    base = SEED_FACTOR*(1+MASS_BOUND)
    require((iv.ln(base)+p2_log-lambda_log).b < 0 and iv.ln(HIGHER_FACTOR).b < lambda_log.a,
            'zero-swirl induction does not close')
    # With M0=10 and positive-index norm at most D, the upper bound on
    # |rhs/(epsilon D)| is 10(n+1)(1003 n +1204)+1. Divide by (n+1)^2.
    # (1003 n+1204)/(n+1)<=1204, so a single constant works at all n.
    require(10*1204+1 < HIGHER_FACTOR and 50-MAX_ORDER >= 25, 'induction or final Cauchy capsule failed')

    # Use a rigorous positive scale for epsilon and epsilon*Pstar^2.
    # These are interval factors, never signed field values.
    eps = _dyadic_inv_p2(iv,lambda_log/2).b
    ep2 = _dyadic_inv_p2(iv,(lambda_log-p2_log)/2).b
    invp2 = _dyadic_inv_p2(iv,read(iv,core['logPstar'])).b
    D = (base*ep2).b
    higher = (HIGHER_FACTOR*base*eps).b
    return dict(D=D,higher=higher,invp2=invp2,S=read(iv,logs['selected_phase']),
                hb_upper=iv.mpf(2)**(-1024),delta=delta)


def _poly_derivative_mass(i, start=2):
    """Derivative weight at rho=4 for finite orders start through 24."""
    return sum(prod(range(n-i+1,n+1))*4**(n-i) for n in range(max(i,start),MAX_ORDER+1))


def _inlet_viscous_radius(iv, p, curvature):
    # With |Z|<=1, |d|<=1 and |L|>=1/2, bound the derivative of the
    # original quotient separately from its numerator.
    a,az,azz = p
    b,bz = curvature
    n0 = az+4*a+8*b
    n1 = azz+2*az+4*(a+az)+8*(b+bz)
    return (32*n0).b, (32*(n1+n0)).b


def _axis_and_finite_inertia(iv, data):
    # On the terminal coefficient capsule, perturb U and its own mean
    # (and their radial derivatives) by the complete positive radial sector.
    # Three 1/32 losses cover Q, Q_Z, and the final ordinary derivative.
    D = data['D']
    u = (_poly_derivative_mass(0,start=1)*D).b
    ur = (_poly_derivative_mass(1,start=1)*D).b
    q0, q0z = 160, 32*160
    q, qz, qy = 332*u,32*332*u,4*332*ur
    b0 = 2*(4*q0+5*q0z)
    b = 2*(4*q+5*qz+4*qy)
    axis = 2*(q0+2*q0z)+q0**2//2+10*b0
    gap = 2*(q+2*qz+qy)+q*(q0+q/2)+(q0+q)*qy+u*b0+(10+u)*b
    return (8*data['invp2']*gap).b, axis


def _whole_remainder(iv, inputs, data):
    core, tube = inputs['core'], inputs['tube']
    p = inputs['pressure']['compliant_source']
    actual = inputs['whole']['directed_radii']
    require(core['contraction_proved'] and p['stages']['z_flatten']['beta']=='variable [0,2]'
            and len(inputs['whole']['cells'])==128
            and data['S'].a == data['S'].b and 0 < data['S'].a < 1,
            'wrong actual datum or radius-four source')
    norm = read(iv,core['scaled_map_size_upper'])/(1-read(iv,core['scaled_map_Lipschitz_upper']))
    h = read(iv,tube['Xh_parameter'])
    flat = read(iv,p['stages']['z_flatten']['mass']).b
    m0 = read(iv,p['fixed_beta0_mass_normalized']).b
    # On the radius-1/4 capsule the late pressure slope obeys
    # |Z*(gamma*M0+d*I1+gamma*I0)/L| <= 8*M0+352*flat_mass.
    late_slope = 8*m0+352*flat
    late = tuple((factorial(k)*8**k*late_slope).b for k in range(3))
    true_nonlinear = tuple((data['invp2']*norm*factorial(1+k)/
                      (20*(k+1)**2*h**k*(iv.mpf(4)/5)**(2+k))).b for k in range(3))
    finite_high = {i:tuple((factorial(k)*8**k*_poly_derivative_mass(i)*data['higher']).b
                           for k in range(3 if i == 1 else 2)) for i in (1,2)}
    finite_visc = _inlet_viscous_radius(iv,finite_high[1],finite_high[2])
    true_late_visc = _inlet_viscous_radius(iv,late,(0,0))
    true_nonlinear_visc = _inlet_viscous_radius(iv,true_nonlinear,(0,0))

    finite_inertia, axis_norm = _axis_and_finite_inertia(iv,data)
    y_max = 2*data['hb_upper']*data['S']
    radius_square_gap = 16*y_max/(1-y_max)
    axis_geometry = (data['invp2']*axis_norm*radius_square_gap/2).b
    radii = dict(late_pressure=true_late_visc, true_inlet_nonlinear=true_nonlinear_visc,
        finite_higher_viscous=finite_visc, finite_non_axis_inertia=(finite_inertia,(32*finite_inertia).b),
        shared_axis_radial_weight=(axis_geometry,(32*axis_geometry).b))
    radii.update({key:tuple(read(iv,r).b for r in actual[key])
                  for key in ('regular_viscous','inertial_nonlinear','cutoff_swirl')})
    return radii,finite_high


def evaluate(*, inputs, pins):
    guards = _source_guards(inputs)
    iv,mp = MPIntervalContext(),MPContext()
    iv.prec = 1800
    q2, core, finite, whole = (inputs[n] for n in ('q2','core','finite','whole'))
    datum = inputs['pressure']['compliant_source']
    family = q2['source_family']
    require(q2['checks']['passed'] and whole['passed'] and finite['passed']
            and all(datum[k] == core[k] == family[k] for k in ('implicit_source_sha256','datum_enclosure_sha256'))
            and whole['frozen_inputs']['experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_correlated_q2_inlet.json']
                == pins['q2']
            and finite['frozen_inputs']['experiments/root_st073/original_collar_regular_whole_z_v2_receipt.json.gz']
                == pins['whole'], 'whole-collar join or core family changed')
    data = _bounds(iv,q2,core)
    masses = _whole_remainder(iv,inputs,data)
    radii,higher_jets=masses
    combined = tuple(sum(row[k] for row in radii.values()).b for k in range(2))

    def encode(value):
        return dict(lower_exact_mpf_tuple=list(value._mpi_[0]),upper_exact_mpf_tuple=list(value._mpi_[1]),
                    rounded=[mp.nstr(mp.make_mpf(bits),17) for bits in value._mpi_],sign=sign(value))

    reference = -iv.mpf(finite['finite_reference_pressure_decimal'])
    require(0 < reference.a <= reference.b < MASS_BOUND,'saved reference mass outside theorem')
    full_reference_gap = []
    for row in finite['pressure_mass_sector_cells']:
        full_reference_gap.append(dict(axial_cell=row['axial_cell'], ordinary_C0_C1=[
            encode(read(iv,part)+hull(iv,-radius,radius))
            for part,radius in zip(row['positive_beta2_mass_gap_C0_C1'],combined,strict=True)]))
    require(len(full_reference_gap) == 128,'missing whole-axis pressure mass sector')

    return dict(passed=True, source_guards=guards, nominal_maximum_radial_order=MAX_ORDER,
                domain=['-1','1'], ordinary_Z_orders=[0,1],
                finite_mass_gate=f'abs(Mf)<={MASS_BOUND}',
                finite_mass_gate_for_exact_legacy_binary_scalar_proved=False,
                formula='actual_row4 - nominal_real_arithmetic_row4 = (M2-Mf)*viscous_mass_coefficient + E',
                normalized_unit='epsilon^2*hb*Pstar^2*S, S=sc/2',
                actual_and_finite_hydrodynamic_inlet_cancel_algebraically=True,
                higher_nominal_u_over_epsilon_p2_upper=encode(data['higher']),
                higher_inlet_ordinary_jet_radii={str(i):[encode(v) for v in row] for i,row in higher_jets.items()},
                whole_axis_remainder_C0_C1=[encode(v) for v in combined],
                remainder_components={key:[encode(v) for v in vals] for key,vals in radii.items()},
                explicit_reference_mass=encode(reference),
                full_gap_to_finite_nominal_model_at_explicit_reference_mass=full_reference_gap,
                explicit_reference_is_the_legacy_internal_scalar=False,
                archival_floating_point_roundoff_certified=False,
                inherited_error_packets_removed=False, forward212_terms_added_or_removed=0,
                independent_terminal_matching_closed=False, ancestor_calculations_replayed=False)
