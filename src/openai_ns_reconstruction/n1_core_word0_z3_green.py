"""Signed whole-cell core Green integral for ordinary-Z3 forcing sectors.

The supplied cell forcings must be outward interval enclosures of the same
actual core sources. F0(Z) and Lambda are shared positive formal factors;
no source exponential is materialized and no cell coefficient is replaced
with its midpoint or absolute value.
"""
from __future__ import annotations


def integrate(c, cells):
    """Integrate at X=2, with eight t-cells and rho_source=4*t**2.

    ``cells`` are ordered mappings with ``t``, ``rho`` and ``g``. The latter
    maps rows 4,5,6 to sector lists containing ``F0_power``, ``Lambda_power``
    and directed ``coefficient``. Endpoint Z-independence and source provenance
    must be validated by the caller.
    """
    if len(cells) != 8:
        raise ValueError('exactly eight core t-cells required')
    fields = (('F1','5'),('Uz1','6'),('K1','6'),('P1','4'),('d_x_F1','5'),('d_x_Uz1','6'))
    output = {name: [] for name, _ in fields}
    for index, cell in enumerate(cells):
        t = cell['t']
        lower = c.mpf(index)/8
        upper = c.mpf(index+1)/8
        if not (t.a == lower.a and t.b == upper.b):
            raise ValueError('wrong directed t-cell endpoints')
        if not (t.a >= 0 and t.b <= 1):
            raise ValueError('core Green parameter outside [0,1]')
        expected_rho = 4*t*t
        if cell['rho']._mpi_ != expected_rho._mpi_:
            raise ValueError('cell rho must be the directed enclosure 4*t*t')
        if index == 0:
            # 0 <= t log(1/t) <= 1/e < 1/2 on [0,1].
            tlog = c.mpf([0, '0.5'])
        else:
            tlog = t*c.log(1/t)
        K = t*(1-t*t)  # X*t*(1-t²)/2 with X=2
        kernels = {'F1': K, 'Uz1': 2*tlog, 'K1': -K,
                   'P1': c.mpf(1), 'd_x_F1': t**3, 'd_x_Uz1': t}
        measure = c.mpf(2)/8
        for field, row in fields:
            for term in cell['g'][row]:
                output[field].append({
                    'cell': index, 'source_row': row,
                    'F0_power': term['F0_power'], 'Lambda_power': term['Lambda_power'],
                    'coefficient_interval': kernels[field]*measure*term['coefficient'],
                })
    return {'fields': output, 'coordinate': 'core endpoint rho=4, X=2',
        'ordinary_Z_order': 3, 'cell_count': 8,
        'integration': 'directed whole-cell inclusion, not point quadrature',
        'source_exponentials_materialized': False,
        'cell_sectors_remain_separate': True,
        'no_first_collar_or_allword_tail_added': True,
        'absolute_accuracy_or_sign_claim': False}
