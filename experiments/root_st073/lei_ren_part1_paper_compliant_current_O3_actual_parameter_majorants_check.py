"""Actual scalar enclosure, factored finite-N and physical density checks."""
import json
import math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_O3_actual_parameter_majorants as source


def independent_product_rule():
    x,z=s.symbols('local_logR local_Z')
    A=[[s.Symbol('a%d_%d'%(j,k)) for k in range(6)] for j in range(5)]
    B=[[s.Symbol('b%d_%d'%(j,k)) for k in range(6)] for j in range(5)]
    polynomial=lambda rows:sum(rows[j][k]*x**j*z**k/(s.factorial(j)*s.factorial(k))
        for j in range(5) for k in range(6))
    direct=s.Poly(s.expand(polynomial(A)*polynomial(B)),x,z)
    caps=source.product_caps(A,B)
    for j in range(5):
        for k in range(6):
            exact=direct.coeff_monomial(x**j*z**k)*s.factorial(j)*s.factorial(k)
            if s.expand(caps[j][k]-exact)!=0:
                raise ArithmeticError('Independent two-variable density product differs')
    return 30


def moderate_parameter_fixture():
    # Independent high-precision evaluation of the original symbolic
    # bounds lies inside the directed adapter with non-microscopic mu.
    c=MPIntervalContext();c.dps=75;mu=c.mpf('.001');N=37
    norms=source.normalized_source_norms(c,mu)
    actual=source.normalized_envelopes(c,mu,N,norms['C'],norms['Fhat'],norms['L1'])
    original=source.actual_source_norms(s.Rational(1,1000),1)
    reference=source.profile_majorants(s.Rational(1,1000),N,
        original['actual_cutoff_ordinary_logR_majorants'],
        original['original_theta_over_Pstar_logR_axial_majorants'],
        original['original_log_theta_logR_majorant'])
    count=0
    with mp.workdps(100):
        for key,rows in actual.items():
            target=reference[key.replace('_over_Ad','')]
            pairs=[]
            if isinstance(rows,list):
                for j,row in enumerate(rows):
                    pairs.extend(zip(row,target[j]) if isinstance(row,list) else [(row,target[j])])
            else:pairs=[(rows,target)]
            for row,exact in pairs:
                value=mp.mpf(str(s.N(exact,95)))
                lo,hi=source.numeric.transport.endpoints(row)
                if not lo<=value<=hi:
                    raise ArithmeticError('Independent symbolic source cap missed by directed recurrence: '+key)
                count+=1
    return count


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest: raise ValueError('Changed actual finite-N scalar/source: '+name)
    field=source.CurrentO3ActualParameterMajorants(require_checked=False)
    if data['actual_parameter_majorant_definition_sha256']!=field.definition:
        raise ValueError('Foreign scalar-bound finite-N definition')
    if source.encoded(field.theorem)!=data['exact_actual_parameter_and_factor_theorem']:
        raise ValueError('Actual parameter/factor theorem differs')
    for N in (1,10**12):
        if source.encoded(field.bounds(N))!=data['examples'][str(N)]:
            raise ValueError('Saved actual finite-N example differs')
    if data['actual_parameter_finite_N_modulation_and_density_bounds_available'] is not True:
        raise ValueError('Actual parameter query gate missing')
    if any(data.get(key) is not False for key in source.OPEN):
        raise ValueError('Finite-N source adapter admits unfinished repair/cone/recursion')
    fresh=field.bounds(37);physical_count=0
    for label,(a,p,r,h) in source.UNITS.items():
        row=fresh['physical_five_defect_density_per_dX_majorants'][label]
        if (row['Ad_power'],row['Pstar_power'],row['logR_power'],row['sqrt2_power'])!=(a,p,r,h):
            raise ValueError('Actual density factor ledger differs')
        for values in row['log_absolute_upper']:
            for value in values:
                if not all(mp.isfinite(v) for v in source.numeric.transport.endpoints(value)):
                    raise ValueError('Nonfinite physical density bound')
                physical_count+=1
    original=fresh['source_scalar_enclosures']
    if original['mu']._mpi_!=field.data['mu']._mpi_ or original['exact_J_at1']._mpi_!=field.ctx.mpf('.5')._mpi_:
        raise ValueError('Actual scalar uncertainty or exact original integral lost')
    # The tiny positive increment must survive; exp(H0)-1 would round to
    # zero at this precision. The integral bound H0*exp(H0) remains positive.
    increment=fresh['finite_N_modulation_majorants']['theta_increment_over_Pstar_over_Ad_majorants'][0][0]
    if source.numeric.transport.endpoints(increment)[0]<=0:
        raise ArithmeticError('Microscopic positive increment underflowed or cancelled')
    rejected=0
    for N in (0,-1,True,1.5,'37',mp.inf,None):
        try:field.bounds(N)
        except (ValueError,TypeError):rejected+=1
        else:raise ArithmeticError('Invalid finite N accepted')
    products=independent_product_rule();fixture=moderate_parameter_fixture()
    hashes=dict(data['input_hashes'])
    hashes[source.NAME]=source.sha(source.NAME);hashes[Path(__file__).name]=source.sha(Path(__file__).name)
    receipt=dict(all_passed=True,actual_parameter_majorant_definition_sha256=field.definition,
        source_family=field.data['source']['accepted']['source_family'],
        exact_actual_source_parameter_factor_identities=len(field.theorem['identities']),
        independent_two_variable_density_product_identities=products,
        independent_moderate_parameter_symbolic_enclosures=fixture,
        new_N37_actual_physical_density_derivative_bounds=physical_count,
        actual_mu_and_log_mu_defining_equations_enclosed=True,original_J_at1_exact=True,
        microscopic_positive_increment_retained_without_exp_subtraction=True,
        original_amplitude_and_giant_radius_exponentials_not_materialized=True,
        invalid_N_inputs_rejected=rejected,
        actual_parameter_finite_N_modulation_and_density_bounds_available=True,
        **{key:False for key in source.OPEN},input_hashes=hashes)
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Actual finite-N parameter bounds PASS:',len(field.theorem['identities']),'identities;',
        fixture,'independent scalar caps;',physical_count,'new physical density rows;',rejected,'guards',flush=True)
    return receipt


if __name__=='__main__':run()
