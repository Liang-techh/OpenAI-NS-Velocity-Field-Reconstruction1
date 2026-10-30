"""Propagate separate raw word and tail bounds through the seam operator.

All input intervals are caller-owned and source-bound.  A positive radius is
never interpreted as a signed source value.  In particular, the finite-word
known signed sectors and the omitted all-word tail remain separate.
"""
from __future__ import annotations


def propagate(previous, *, word_zz_upper, tail_upper, tail_z_upper,
              tail_zz_upper, absolute_product_upper):
    """Return independent source-sector radii for rows 4, 5, and 6.

    ``absolute_product_upper(saved_coefficient_interval, positive_bound)``
    must return an outward-rounded, nonnegative product upper bound.  The
    function preserves each positive exponential's saved log record and never
    merges distinct factors.  The caller must validate common source, tube,
    units and all Cauchy premises before invoking this adapter.
    """
    matrices = previous['matrices']
    out = {}
    for row in ('4', '5', '6'):
        saved = previous['rows'][row]
        missing = saved['Z1_omitted_remainder']['sectors']
        word_errors = []
        for term in missing:
            if term['fully_evaluated'] is not False or term['cannot_omit'] is not True:
                raise ValueError('missing W1_ZZ must remain symbolic')
            word_errors.append({
                'column': term['column'],
                'source_exponents': term['source_exponents'],
                'combined_source_log_interval': term['combined_source_log_interval'],
                'radius': absolute_product_upper(term['coefficient_interval'], word_zz_upper),
                'origin': '-B1[row,column]*W1_ZZ[column]',
            })

        def terms(matrix, order, bound, origin):
            result = []
            for key, factors in matrices[matrix]['Z' + str(order)].items():
                rr, column = (int(x) for x in key.split(','))
                if rr != int(row):
                    continue
                for term in factors:
                    result.append({
                        'column': column,
                        'source_exponents': term['source_exponents'],
                        'combined_source_log_interval': term['combined_source_log_interval'],
                        'radius': absolute_product_upper(term['coefficient_interval'], bound),
                        'origin': origin,
                    })
            return result

        tail_z0 = (terms('C0', 0, tail_upper, 'C0*E')
                   + terms('B1', 0, tail_z_upper, 'B1*E_Z'))
        tail_z1 = (terms('C0', 1, tail_upper, 'C0_Z*E')
                   + terms('C0', 0, tail_z_upper, 'C0*E_Z')
                   + terms('B1', 1, tail_z_upper, 'B1_Z*E_Z')
                   + terms('B1', 0, tail_zz_upper, 'B1*E_ZZ'))
        out[row] = {
            'finite_Z1_known_signed_sectors': saved['Z1_known']['sectors'],
            'finite_Z1_missing_independent_radii': word_errors,
            'finite_Z1_enclosure_formula': 'known signed sectors + sum of symmetric missing-term radii',
            'all_word_tail_Z0_operator_radii': tail_z0,
            'all_word_tail_Z1_operator_radii': tail_z1,
        }
    return {
        'rows': out,
        'word_and_tail_budgets_separate': True,
        'radii_are_not_signed_values': True,
        'source_log_terms_not_collapsed': True,
        'all_word_tail_not_added_to_finite_word': True,
        'actual_m15_integral_evaluated': False,
        'all_Z_functional_coverage': False,
    }
