"""Join a core W0_Z3 enclosure to the inlet B1 feedback factors.

The caller supplies source-bound interval terms and verifies that the core
endpoint is the actual first-collar inlet. No logarithmic factors are
collapsed or materialized; W0's current-Z F0 and B1's original source factor
remain separately identified even when they refer to the same function.
"""
from __future__ import annotations


def join(previous_rows, core_fields, *, interval):
    fields={1:'F1',2:'Uz1',3:'K1',4:'P1',5:'d_x_F1',6:'d_x_Uz1'}
    out={}
    for row in ('4','5','6'):
        previous=previous_rows[row]; known=previous['known_signed_sector_terms']; completed=[]
        for b1 in previous['missing_B1_W0_ZZZ']:
            if b1['cannot_omit'] is not True:
                raise ValueError('missing feedback factor must remain identified')
            column=b1['column']; name=fields[column]
            word=core_fields[name]
            if not word:
                raise ValueError('missing actual W0_Z3 field sectors')
            for w0 in word:
                completed.append({
                    'origin':'B1_W0_ZZZ', 'column':column, 'word_field':name,
                    'word_cell':w0['cell'],
                    'B1_source_exponents':b1['source_exponents'],
                    'B1_combined_source_log_interval':b1['combined_source_log_interval'],
                    'W0_F0_current_Z_power':w0['F0_power'],
                    'W0_Lambda_power':w0['Lambda_power'],
                    'coefficient_interval':interval(b1['coefficient_interval'])*interval(w0['coefficient_interval']),
                })
        out[row]={'previous_known_signed_terms':known,'B1_W0_ZZZ_enclosed_terms':completed,
            'full_source_formula':'sum(previous_known_signed_terms) + sum(B1_W0_ZZZ_enclosed_terms)',
            'termwise_interval_enclosure_not_a_point_value':True}
    return {'rows':out,'F0_and_B1_source_factors_kept_separate':True,
        'same_source_correlation_not_assumed_from_interval_overlap':True,
        'existing_Cauchy_bound_not_removed':True,'existing_allword_tail_not_removed':True,
        'W1_Z2_Green_integral_computed':False,'m15_integral_computed':False,
        'scope':'one-sided rho=4 seam and previously selected narrow real Z interval'}
