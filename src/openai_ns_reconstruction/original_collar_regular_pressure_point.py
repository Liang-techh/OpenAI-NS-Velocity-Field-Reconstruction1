"""Normalize the selected C1 regular-collar pressure packet in formal scales."""
from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext

BASIS = ['hb','F0_anchor','Pstar_squared','Lambda','sigma_a','sigma_2a']
SOURCES = {'direct_pressure':1, 'cut_swirl':1, 'principal':-1, 'remainder':-1}


def _read(iv, row):
    return iv.make_mpf(tuple(tuple(row[f'{side}_exact_mpf_tuple']) for side in ('lower','upper')))


def _hull(iv, lo, hi):
    return iv.make_mpf((lo._mpi_[0],hi._mpi_[1]))


def _record(mp, value):
    keys=('lower_exact_mpf_tuple','upper_exact_mpf_tuple')
    return {**{key:list(v) for key,v in zip(keys,value._mpi_,strict=True)},
            'display':[mp.nstr(mp.make_mpf(v),16) for v in value._mpi_]}


def _scale(iv, coefficient, log_scale, floor, ceiling):
    # Apply a positive factor to each endpoint in log space, so a very small
    # signed term keeps an outward interval even when its weight is formal.
    def endpoint(x):
        if x == 0:
            return iv.mpf(0)
        log_value=iv.ln(abs(x))+log_scale
        if log_value.b < floor.a:
            value=_hull(iv,iv.mpf(0),ceiling)
        else:
            if not -2048 < log_value.a <= log_value.b < 2048:
                raise ValueError('non-negligible coefficient log did not cancel into the target unit')
            value=iv.exp(log_value)
        return value if x>0 else -value

    return _hull(iv,endpoint(coefficient.a).a,endpoint(coefficient.b).b)


def enclose(q2_report):
    q=q2_report['actual_correlated_q2_inlet']
    packet=q['shared_q2_inlet']
    if (not q['passed'] or q['source_log_basis_names']!=BASIS or packet['Z']!='371/1000'
            or packet['q']!='2' or packet['components'][3]['source_weights']!=SOURCES):
        raise ValueError('actual selected q2 regular pressure sheet required')
    iv,mp=MPIntervalContext(),MPContext()
    iv.prec=1800
    bases=[_read(iv,v) for v in q['source_log_bases']]
    S=_read(iv,q['selected_phase'])
    if not (S.a == S.b and S.a > 0):
        raise ValueError('fixed selected q2 endpoint required')
    reference=(1,0,1,-2,0,0)
    scale_log_s=iv.ln(S)
    floor,ceiling=-1024*iv.ln(2),iv.mpf(2)**(-1024)
    def scale(row):
        powers=tuple(row['source_exponents'])
        if len(powers)==5:
            powers=(*powers,0)
        log_ratio=sum(((j-k)*x for j,k,x in zip(powers,reference,bases,strict=True) if j!=k),-scale_log_s)
        return _scale(iv,_read(iv,row['coefficient_interval']),log_ratio,floor,ceiling)

    values={}
    for order in ('0','1'):
        leaves=q['ordinary_Z'][order]['source_leaves']
        approximation=radius=iv.mpf(0)
        radii={}
        for name,sgn in SOURCES.items():
            leaf=leaves[name]
            model=sum((scale(row) for row in leaf['signed_approximation']),iv.mpf(0))
            bound=sum((abs(scale(row)).b for row in leaf['error_radius']),iv.mpf(0)).b
            approximation+=sgn*model
            radius+=bound
            radii[name]=_record(mp,bound)
        interval=approximation+_hull(iv,-radius,radius)
        values[order]=dict(interval=_record(mp,interval), signed_model_interval=_record(mp,approximation),
                           uncertainty_by_source=radii, uncertainty_total_upper=_record(mp,radius.b),
                           strict_sign='negative' if interval.b<0 else 'positive' if interval.a>0 else 'unresolved')
    return dict(Z='371/1000', q='2', units='epsilon^2*hb*Pstar^2*S', ordinary_Z=values,
                signed_weights=SOURCES, small_terms_retained_as_intervals=True,
                uncertainty_radii_are_function_values=False, legacy_ledger_equality=False)
