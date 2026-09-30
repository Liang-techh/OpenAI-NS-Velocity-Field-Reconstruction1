"""Arbitrary-exponent Section 7.5 affine and energy solves.

Caller must supply actual normalized incoming/pulse rows and the actual
energy target. No synthetic targets or freely chosen axial coefficients.
Quadrature weights have finite precision; receipts do not certify a cone.
"""
import json
from pathlib import Path
import mpmath as mp
import numpy as np
from lei_ren_part1_paper_outer_closure import _paper_raw_bump
from lei_ren_part1_paper_axial_pulse import pulse_K_p, pulse_value


# The unified profile replays the translated end bumps with order 192.  Use
# that same order for the coefficient solve so its exterior primitive is a
# replay of the solved equations, rather than a different quadrature problem.
CONSISTENT_QUADRATURE_ORDER = 192


def signed_log(value,precision=120):
    return {'sign':int(mp.sign(value)),
            'log_abs':mp.nstr(mp.log(abs(value)),precision) if value else None,
            'arbitrary_exponent_value':mp.nstr(value,precision)}


def from_signed_log(row):
    # A rounded log of size 1e152 loses 152 digits of relative information
    # on exponentiation. Preserve the serialized mantissa when available;
    # tiny nonzero exterior means must not be altered by a log roundtrip.
    if 'arbitrary_exponent_value' in row:
        value = mp.mpf(row['arbitrary_exponent_value'])
        if int(mp.sign(value)) != int(row['sign']):
            raise ValueError('Signed-log value and declared sign disagree')
        return value
    return mp.mpf(0) if not row['sign'] else row['sign']*mp.exp(mp.mpf(row['log_abs']))


def solve_actual_axial(mu, base_rows, pulse_rows, energy_target,*,precision=120,order=128):
    with mp.workdps(precision):
        mu=mp.mpf(str(mu)); target=mp.mpf(str(energy_target))
        if not 0<mu<=mp.mpf(1)/60 or target<=0:
            raise ValueError('Source mu and actual positive energy target required')
        requested_order = int(order)
        if requested_order < 16 or requested_order != order:
            raise ValueError('order must be an integer at least 16')
        effective_order = max(requested_order, CONSISTENT_QUADRATURE_ORDER)
        nodes,weights=np.polynomial.legendre.leggauss(effective_order)
        ell=mp.mpf('.15')
        raw=[mp.mpf(str(float(_paper_raw_bump(float(node))))) for node in nodes]
        normalization=sum(mp.mpf(str(float(w)))*v for w,v in zip(weights,raw))
        support=[ell*mp.mpf(str(float(n))) for n in nodes]
        measure=[ell*mp.mpf(str(float(w))) for w in weights]
        beta=[v/(ell*normalization) for v in raw]
        lambdas=[mp.mpf('.5')-i*mu for i in (1,2)]
        B=[sum(w*mp.exp(l*s)*v for w,s,v in zip(measure,support,beta)) for l in lambdas]
        matrix=[[mp.exp(-3*l)*b,mp.exp(-l)*b] for l,b in zip(lambdas,B)]
        determinant=-2*mp.exp(-2+6*mu)*mp.sinh(mu)*B[0]*B[1]
        base=[from_signed_log(row) for row in base_rows]
        pulse=[mp.exp(mp.mpf(row['log_normalized_pulse_integral'])) for row in pulse_rows]
        def affine(rhs):
            return [(-rhs[0]*matrix[1][1]+matrix[0][1]*rhs[1])/determinant,
                    (-matrix[0][0]*rhs[1]+rhs[0]*matrix[1][0])/determinant]
        u=affine(base); v=affine(pulse)
        gram=sum(w*mp.exp(-2*mu*s)*b*b for w,s,b in zip(measure,support,beta))
        K=[mp.exp(-26+factor*mu)*gram for factor in (6,2)]
        Kp=mp.mpf(str(pulse_K_p(quadrature_order=effective_order)))
        quadratic=Kp+mu*sum(k*b*b for k,b in zip(K,v))
        linear=2*mu*sum(k*a*b for k,a,b in zip(K,u,v))
        constant=mu*sum(k*a*a for k,a in zip(K,u))-target
        discriminant=linear*linear-4*quadratic*constant
        if constant>=0 or discriminant<=0:
            raise ArithmeticError('Actual axial energy target lacks positive branch')
        a=-2*constant/(linear+mp.sqrt(discriminant))
        if not mp.mpf('.9')<a<mp.mpf('1.2'):
            raise ArithmeticError('Actual energy root outside source interval')
        c=[x+a*y for x,y in zip(u,v)]
        # Replay per row normalized by its own RHS: a small absolute error
        # must not be mistaken for closure of a tiny moment.
        relative=[]
        for index in range(2):
            rhs=base[index]+a*pulse[index]
            defect=sum(matrix[index][j]*c[j] for j in range(2))+rhs
            relative.append(mp.nstr(abs(defect/rhs),30) if rhs else None)
        energy=quadratic*a*a+linear*a+constant
        return {'a_p':mp.nstr(a,precision),'c':[signed_log(x,precision) for x in c],
                'linear_matrix':[[mp.nstr(x,precision) for x in row] for row in matrix],
                'linear_rhs_inputs':{'base':base_rows,'pulse':pulse_rows},
                'input_mu':mp.nstr(mu,precision),
                'u':[signed_log(x,precision) for x in u],
                'v':[signed_log(x,precision) for x in v],
                'energy_target':mp.nstr(target,precision),
                'linear_relative_replay':relative,
                'energy_relative_replay':mp.nstr(abs(energy/target),30),
                # Keep the replay weights at solve precision.  A 30-digit
                # receipt is enough for display, but is not a faithful input
                # when a downstream energy replay asks for 1e-40 closure.
                'K_p':mp.nstr(Kp,precision),
                'K_bump':[mp.nstr(x,precision) for x in K],
                'K_p_source_precision':'Inherited from pulse_K_p float quadrature; '
                                       'the shared replay is exact for this stored input, '
                                       'not a claim of arbitrary-precision source integration.',
                'requested_quadrature_order':requested_order,
                'quadrature_order':effective_order,'precision':precision,
                'coefficient_and_profile_quadrature_consistent': True,
                'full_outer_closed':False,
                'scope':'Actual-input algebra solve; inherited quadrature/heat uncertainties, corrected pressure, derivatives/cone/core still open.'}


def axial_factor(receipt,*,xi=None,end_offset=None,precision=120):
    """U_z/E in pulse xi=mu*t or end-bump offset t-13/mu.

    Separate coordinates avoid losing unit-width end bumps at huge t.
    """
    if (xi is None)==(end_offset is None):
        raise ValueError('Supply exactly one of xi or end_offset')
    with mp.workdps(precision):
        if xi is not None:
            result=mp.mpf(receipt['a_p'])*mp.mpf(str(pulse_value(xi)))
        else:
            # The pulse is zero here, since its support ends at xi=11.
            nodes,weights=np.polynomial.legendre.leggauss(receipt['quadrature_order'])
            # Match the solve and the cumulative primitive's canonical
            # normalization. A binary64 sum here would change the actual
            # pointwise bump while leaving the moment equations unchanged.
            normalization=sum(mp.mpf(str(float(w)))*mp.mpf(str(float(_paper_raw_bump(float(n)))))
                for n,w in zip(nodes,weights))
            result=mp.mpf(0)
            offset=end_offset if isinstance(end_offset,mp.mpf) else mp.mpf(str(end_offset))
            for data,center in zip(receipt['c'],(-3.,-1.)):
                unit=(offset-mp.mpf(str(center)))/mp.mpf('.15')
                bump=mp.mpf(str(float(_paper_raw_bump(float(unit)))))/(mp.mpf('.15')*normalization)
                if bump:result+=from_signed_log(data)*bump
        return signed_log(result,precision)


def independent_end_bump_replay(mu,receipt,base_rows,pulse_rows,*,order=192,precision=160):
    """Reintegrate actual translated bumps with a different quadrature.

    Keep lambda and the nearly cancelling products at arbitrary precision.
    A binary64 replay would lose the O(mu) determinant entirely.
    """
    with mp.workdps(precision):
        mu=mp.mpf(str(mu)); a=mp.mpf(receipt['a_p'])
        c=[from_signed_log(row) for row in receipt['c']]
        nodes,weights=np.polynomial.legendre.leggauss(order)
        raw=[mp.mpf(str(float(_paper_raw_bump(float(x))))) for x in nodes]
        normalization=sum(mp.mpf(str(float(w)))*v for w,v in zip(weights,raw))
        errors=[]; absolute=[]
        for row_index in (1,2):
            lam=mp.mpf('.5')-row_index*mu; lhs=mp.mpf(0)
            for center,coefficient in zip((-3,-1),c):
                for node,weight,value in zip(nodes,weights,raw):
                    s=mp.mpf('.15')*mp.mpf(str(float(node)))
                    lhs+=mp.mpf(str(float(weight)))*mp.exp(lam*(center+s))\
                         *value/normalization*coefficient
            rhs=from_signed_log(base_rows[row_index-1])+a*mp.exp(
                mp.mpf(pulse_rows[row_index-1]['log_normalized_pulse_integral']))
            defect=lhs+rhs
            absolute.append(signed_log(defect,precision))
            errors.append(mp.nstr(abs(defect/rhs),30))
        return {'row_relative_differences':errors,
                'row_absolute_defects':absolute,
                'quadrature_order':order,
                'scope':'Independent translated-bump quadrature; incoming/pulse integration uncertainty remains inherited.'}


def run():
    from lei_ren_part1_paper_axial_energy_tail import build_default_tail
    from lei_ren_part1_paper_axial_incoming import incoming_axial_moments
    from lei_ren_part1_paper_axial_pulse_moments import normalized_pulse_integral
    tail=build_default_tail(precision=160,quadrature_order=128)
    schedule=tail.schedule
    pulse=[normalized_pulse_integral(schedule.mu,i,precision=160) for i in (1,2)]
    rows=[]
    for z in (0.,.5,-.5):
        incoming=incoming_axial_moments(schedule,z,order=64,precision=160)
        future=tail.evaluate(z,quadrature_order=128)
        base=[incoming['row_normalization'][key] for key in ('scaled_base_m1','scaled_base_m2')]
        with mp.workdps(160):
            prior=from_signed_log(incoming['E_prior_mu_Mztheta_over_RpEp2'])
            target=(1-mp.exp(-26))/4-prior+mp.mpf(future['energy_target_contribution_nominal'])
            receipt=solve_actual_axial(schedule.mu,base,pulse,mp.nstr(target,160),precision=160)
        for value in receipt['linear_relative_replay']:
            if value is not None and mp.mpf(value)>mp.mpf('1e-80'):
                raise ArithmeticError('Tiny source-row algebra replay failed')
        independent=independent_end_bump_replay(schedule.mu,receipt,base,pulse)
        if any(mp.mpf(v)>mp.mpf('1e-8') for v in independent['row_relative_differences']):
            raise ArithmeticError('Independent end-bump quadrature failed moment replay')
        receipt.update({'Z':z,'incoming_moments':incoming,'future_angular_energy':future,
                        'independent_end_bump_replay':independent,
                        'pulse_rows':pulse,
                        'sample_axial_factors':[axial_factor(receipt,xi=5.,precision=160),
                            axial_factor(receipt,end_offset=-3.,precision=160),
                            axial_factor(receipt,end_offset=-1.,precision=160)]})
        rows.append(receipt)
    report={'source':'https://arxiv.org/html/2609.35406v1','equations':'7.31,7.34',
            'schedule':schedule.metadata(),'actual_input_axial_corrections':rows,
            'scope':'Actual candidate-input pulse/affine/energy solve with inherited numerical and heat bounds. Not an exact global five-moment, pressure, cone or recursive PDE certificate.',
            'full_outer_closed':False}
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({'a_p':[row['a_p'][:24] for row in rows],
                      'linear_replay':[row['linear_relative_replay'] for row in rows],
                      'energy_replay':[row['energy_relative_replay'] for row in rows]}))
    return report


if __name__=='__main__':run()
