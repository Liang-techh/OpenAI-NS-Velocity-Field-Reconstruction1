"""Analytic axial envelopes of the common positive preheat-pressure atoms.

Finite Gauss atom masses do not enclose the radial integral. Keep labeled
atoms, numerical mass-error budgets and analytic factor bounds separate.
"""
import mpmath as mp


def _nonnegative(value,name):
    v=mp.mpf(value)
    if not mp.isfinite(v) or v<0:raise ValueError(name+' must be finite and nonnegative')
    return v


def q_power_derivative_bounds(beta,axial_radius):
    """Upper bounds for derivatives 0..3 of (1+Z^2)^-beta on [-a,a]."""
    b=_nonnegative(beta,'beta');a=_nonnegative(axial_radius,'axial_radius')
    if a>=1:raise ValueError('Require compact axial_radius<1')
    return [mp.mpf(1),2*b*a,2*b+4*b*(b+1)*a*a,
            12*b*(b+1)*a+8*b*(b+1)*(b+2)*a**3]


def positive_stage_mass_upper(length,log_left_upper,log_right_upper,*,density_factor='.5'):
    """Conditional radial mass upper bound from schedule slope envelopes.

    Assume -3/2 <= d_y log(A/Pstar) <= 11/10 and supplied endpoint
    upper bounds for log(A/Pstar). The density is <= factor*(A/Pstar)^2.
    Prefix/flatten use factor=1/2; post-flatten U=A/2 uses factor=1/8.
    Endpoint bounds must be enclosed separately; this routine does not obtain
    them from quadrature-evaluated schedule endpoints.
    """
    L=_nonnegative(length,'length');factor=_nonnegative(density_factor,'density_factor')
    left=mp.mpf(log_left_upper);right=mp.mpf(log_right_upper)
    if not mp.isfinite(left) or not mp.isfinite(right):raise ValueError('Finite endpoint bounds required')
    upper=min(left+mp.mpf('1.1')*L,right+mp.mpf('1.5')*L)
    return dict(log_amplitude_upper=upper,mass_upper=L*factor*mp.exp(2*upper),
                endpoint_bounds_verified=False,rounding_enclosed=False)


def positive_stage_pressure_bounds(length,log_left_upper,log_right_upper,*,
                                   axial_radius,beta_upper=2,density_factor='.5'):
    """Conditional true-stage derivative bounds with beta(y) in [0,beta_upper].

    Unlike per-node mass errors, this bound applies to the entire positive
    radial integral even when beta varies between quadrature nodes. It is
    deliberately loose; it is not an estimate of the Gauss quadrature error.
    """
    b=_nonnegative(beta_upper,'beta_upper')
    if b>2:raise ValueError('Preheat beta_upper must be <=2')
    result=positive_stage_mass_upper(length,log_left_upper,log_right_upper,
                                     density_factor=density_factor)
    result['derivative_bounds']=[result['mass_upper']*f for f in q_power_derivative_bounds(b,axial_radius)]
    result['variable_beta_allowed']=True
    result['source_scope']='conditional positive radial integral, given endpoint and slope bounds'
    return result


def preheat_interval_bounds(receipt,*,axial_radius,mass_error_bounds=None):
    """Bound finite positive-atom pressure and optional per-atom mass errors.

    Consume taylor_components() raw MP data. Error keys are stage or stage:i
    for flattened quadrature atoms. All keys must be supplied when errors
    are provided, including exact-zero budgets. A mass error is assumed to
    multiply the same beta factor; a variable-beta radial quadrature error
    cannot in general be represented this way. Supply a separate universal
    derivative error envelope for that stage before claiming true pressure
    closure. This function never certifies quadrature or actual source data.
    """
    a=_nonnegative(axial_radius,'axial_radius')
    if a>=1:raise ValueError('Require axial_radius<1')
    rows={}
    for stage,component in receipt['components'].items():
        kind=component['kind']
        if kind=='q_power_atoms':
            atoms=[(stage+':'+str(i),atom['atom_at_Z0'],atom['beta'])
                   for i,atom in enumerate(component['atoms'])]
        elif kind=='q_power':atoms=[(stage,component['value_at_Z0'],component['beta'])]
        elif component.get('Z_independent',False):atoms=[(stage,component['value_at_Z0'],0)]
        else:raise ValueError('Unknown analytic axial factor for '+stage)
        for label,mass,beta in atoms:
            if label in rows:raise ValueError('Duplicate atom label')
            w=_nonnegative(mass,'mass');b=_nonnegative(beta,'beta')
            if b>2:raise ValueError('Preheat beta must be in [0,2]')
            factors=q_power_derivative_bounds(b,a)
            error=None if mass_error_bounds is None else _nonnegative(mass_error_bounds[label],'mass_error')
            rows[label]=dict(stage=stage,region=component.get('region'),beta=b,mass=w,
                derivative_bounds=[w*f for f in factors],
                mass_error=error,error_derivative_bounds=None if error is None else [error*f for f in factors])
    if mass_error_bounds is not None and set(mass_error_bounds)!=set(rows):
        raise ValueError('Mass error keys must exactly match retained atom labels')
    total=[mp.fsum(row['derivative_bounds'][i] for row in rows.values()) for i in range(4)]
    errors=None if mass_error_bounds is None else [mp.fsum(row['error_derivative_bounds'][i] for row in rows.values()) for i in range(4)]
    return dict(axial_radius=a,atoms=rows,normalized_derivative_bounds=total,
        normalized_mass_error_derivative_bounds=errors,
        normalized_C2_bound=total[0]+total[1]+total[2]/2,
        Pstar_squared=receipt['Pstar_squared'],atom_count=len(rows),
        scope='analytic factors of the declared finite positive-atom representation',
        supplied_mass_errors_verified=False,variable_beta_quadrature_error_enclosed=False,
        quadrature_error_enclosed=False,actual_source_norm_enclosed=False,
        rounding_enclosed=False,retained_atom_labels=True)
