"""Directed rescaling of the intact selected cut-axial product packet."""
from __future__ import annotations

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext

from .original_collar_whole_z import hull, read, require, sign


def _ceiling(bits):
    negative,mantissa,exponent,_ = bits
    signed = -mantissa if negative else mantissa
    return signed << exponent if exponent >= 0 else -((-signed) >> -exponent)


def enclose(q2):
    """Enclose the actual selected row-three source in the old physical unit.

    The saved signed inlet-product model and each of its error sectors are
    scaled independently. Shifts to a dyadic working unit do not change sign.
    """
    packet = q2['actual_correlated_q2_inlet']
    require(q2['checks']['passed'] and packet['passed'] and packet['shared_q2_inlet']['Z']=='371/1000'
            and packet['shared_q2_inlet']['q']=='2'
            and packet['source_log_basis_names'][:4]==['hb','F0_anchor','Pstar_squared','Lambda'],
            'wrong selected cut axial packet')
    iv,mp = MPIntervalContext(),MPContext()
    iv.prec = 1800
    basis = [read(iv,v) for v in packet['source_log_bases'][:4]]
    a = read(iv,packet['selected_phase'])/2
    require(a.a==a.b and a.a>0,'expected exact sc/4')
    unit = (1,0,0,-2)
    ln2 = iv.ln(2)
    groups = {}
    for k in ('0','1'):
        leaf = packet['ordinary_Z'][k]['source_leaves']['cut_axial']
        groups[k] = {}
        for name in ('signed_approximation','error_radius'):
            rows = []
            for row in leaf[name]:
                powers = row['source_exponents']
                require(len(powers)==5 and powers[-1]==0,'unexpected flat-scale atom in cut axial sector')
                log_factor = sum(((p-q)*log_atom for p,q,log_atom in zip(powers[:-1],unit,basis,strict=True) if p!=q),
                                 -iv.ln(a))
                rows.append((read(iv,row['coefficient_interval']),log_factor))
            groups[k][name] = rows

    # The largest possible scaled magnitude sets the working exponent. Tiny
    # terms still get directed intervals even far below this working scale.
    largest = max((iv.ln(abs(c).b)+f).b for group in groups.values() for rows in group.values()
                  for c,f in rows if c != 0)
    binary = largest/ln2
    require(-10**20 < binary.a <= binary.b < 10**20, 'unresolved scale of the axial packet')
    shift = _ceiling(binary._mpi_[1])
    factor = iv.make_mpf(((0,1,shift,1),)*2)
    tiny = iv.mpf(2)**(-2048)

    def endpoint(x, log_factor):
        if x == 0:
            return iv.mpf(0)
        work = iv.ln(abs(x))+log_factor-shift*ln2
        if work.b < -2048*ln2.b:
            value = hull(iv,iv.mpf(0),tiny)
        else:
            require(-1500 < work.a <= work.b < 10,'unresolved finite axial mantissa')
            value = iv.exp(work)
        return value if x>0 else -value

    def scaled(c,f):
        return hull(iv,endpoint(c.a,f).a,endpoint(c.b,f).b)

    def encode(x):
        return dict(lower_exact_mpf_tuple=list(x._mpi_[0]),upper_exact_mpf_tuple=list(x._mpi_[1]),
                    rounded=[mp.nstr(mp.make_mpf(bits),17) for bits in x._mpi_],sign=sign(x))

    ordinary = {}
    for k, group in groups.items():
        model = sum((scaled(c,f) for c,f in group['signed_approximation']),iv.mpf(0))
        error = sum((abs(scaled(c,f)).b for c,f in group['error_radius']),iv.mpf(0)).b
        verified = (model+hull(iv,-error,error))*factor
        ordinary[k] = dict(signed_approximation=encode(model*factor),error_radius=encode(error*factor),
                           actual_interval=encode(verified))
    return dict(Z='371/1000',unit='epsilon^2*hb*(sc/4)',ordinary_Z=ordinary,
                existing_error_sectors_retained=True, error_radius_used_as_field_value=False,
                signed_model_and_actual_directed_interval_distinct=True,
                working_binary_exponent=shift,whole_axis_coverage=False)
