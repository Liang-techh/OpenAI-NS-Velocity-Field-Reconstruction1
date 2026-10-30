"""Signed known part of the second-Z first-word feedback source.

This leaf evaluates a supplied source-factored jet at one caller-selected
coordinate. It does not integrate the Green kernel, produce W1_ZZ, or replace
the missing W0_ZZZ with a bound. The caller must bind coordinates, sources,
physical units, derivative conventions and outward interval arithmetic.
"""
from __future__ import annotations


def assemble(B0, B1, W0):
    """Form the five available Leibniz terms, preserving each source sector.

    Matrix jets are ``B0[order][(row,column)]`` and ``B1`` likewise. Vector
    jets are ``W0[order][column]``. Every term contains source_exponents,
    interval-valued coefficient, and interval-valued combined_source_log.
    Inputs must be directed intervals in a common interval context.
    """
    specifications = (
        ('C0_ZZ_W0', B0, 2, 0, 1),
        ('2_C0_Z_W0_Z', B0, 1, 1, 2),
        ('C0_W0_ZZ', B0, 0, 2, 1),
        ('B1_ZZ_W0_Z', B1, 2, 1, 1),
        ('2_B1_Z_W0_ZZ', B1, 1, 2, 2),
    )
    output = {}
    for row in (4, 5, 6):
        known = []
        for name, matrix, matrix_order, word_order, multiplier in specifications:
            for (r, col), left_terms in matrix[matrix_order].items():
                if r != row:
                    continue
                if col not in W0[word_order]:
                    raise ValueError('missing actual W0 Z jet column')
                for left in left_terms:
                    for right in W0[word_order][col]:
                        exponent = [a+b for a, b in zip(left['source_exponents'], right['source_exponents'])]
                        if len(exponent) != 4:
                            raise ValueError('source factor dimension must be four')
                        known.append({
                            'origin': name, 'column': col,
                            'source_exponents': exponent,
                            'coefficient_interval': multiplier*left['coefficient_interval']*right['coefficient_interval'],
                            'combined_source_log_interval': left['combined_source_log_interval']+right['combined_source_log_interval'],
                        })
        missing = []
        for (r, col), factors in B1[0].items():
            if r == row:
                for term in factors:
                    missing.append({**term, 'column': col, 'symbolic_factor': f'W0_ZZZ[{col}]',
                                    'cannot_omit': True})
        output[str(row)] = {'known_signed_sector_terms': known, 'missing_B1_W0_ZZZ': missing,
            'known_formula': 'C0_ZZ W0 + 2 C0_Z W0_Z + C0 W0_ZZ + B1_ZZ W0_Z + 2 B1_Z W0_ZZ',
            'missing_formula': 'B1 W0_ZZZ'}
    return {'rows': output, 'location_scope': 'caller-specified point only',
        'known_terms_are_enclosures_not_point_values': True, 'individual_source_terms_preserved': True,
        'W0_ZZZ_computed': False, 'W1_ZZ_integrated': False,
        'whole_collar_coverage': False, 'all_word_tail_added': False}
