"""Admit actual current angular full stress and regional completed tensor."""
import copy
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_current_angular_background_stress import (
    CurrentAngularBackgroundStress,VIEWS,GATES,OPEN,NAME,RECEIPT,HERE,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def finite(value):
    values=endpoints(value) if hasattr(value,'_mpi_') else (value,)
    if not all(mp.isfinite(x) for x in values):
        raise ValueError('Nonfinite actual current normalized function or physical source bound')


def check_view(field,view):
    if any(view[key] for key in OPEN) or not view['actual_full_stress_not_local_difference']:
        raise ValueError('Actual angular tensor scope differs')
    for name in ('A','E','P','K'):
        rows=view['current_actual_normalized_full_moment_rows'][name]
        if len(rows)!=5 or any(row.order!=5 for row in rows):raise ValueError('Full current mixed4 with axial5 moment rows required')
        for row in rows:
            for coefficient in row.coefficients:finite(coefficient)
    for value in view['normalized_energy_native_consistency'].coefficients:
        if not endpoints(value)[0]<=0<=endpoints(value)[1]:
            raise ValueError('Original full future energy consistency failed')
    for name,grid in view['current_actual_normalized_stress_mixed3'].items():
        if len(grid)!=10:raise ValueError('Original full stress mixed3 layout missing: '+name)
        for value in grid.values():finite(value)
    pressure=view['stable_actual_absolute_pressure_mixed4_factored']
    if len(pressure)!=15:raise ValueError('Stable same absolute-pressure mixed4 required')
    for value in pressure.values():finite(value)
    groups=[]
    for key,count in (('physical_cylindrical_stress_mixed3',10),('physical_stress_divergence_mixed2',6)):
        if set(view[key])!={'theta','axial'}:raise ValueError('Both prescribed stress components required')
        for rows in view[key].values():
            if len(rows)!=count:raise ValueError('Actual physical derivative layout differs')
            groups.extend(rows.values())
    for key in ('completed_theta_theta_stress_mixed2','physical_axial_viscosity_remainder_mixed2'):
        if len(view[key])!=6:raise ValueError('Actual completion/remainder mixed2 required')
        groups.extend(view[key].values())
    for key in ('completed_background_tensor_cartesian_components',
            'physical_completed_stress_divergence_cartesian','physical_remainder_cartesian',
            'physical_momentum_residual_decomposition_cartesian'):
        for rows in view[key].values():groups.extend(rows)
    if len(groups)!=65:raise ValueError('Actual angular tensor/divergence/remainder factored layout incomplete')
    nonzero=0
    for row in groups:
        finite(row['signed_coefficient']);finite(row['physical_lambda_exponent'])
        finite(row['physical_viscosity_exponent']);finite(row['radial_log_prefactor'])
        for value in row['actual_source_log_parts'].values():finite(value)
        if endpoints(row['physical_lambda_exponent'])[1]>0:
            raise ValueError('Whole-Z compact-time upper requires nonpositive lambda power')
        if row['exact_zero']:
            if row['log_absolute_upper'] is not None:raise ValueError('Exact-zero physical row has nonzero bound')
        else:
            finite(row['log_absolute_upper']);nonzero+=1
        if not row['positive_source_factors_not_materialized'] or not row['source_factors_combined_before_enclosure']:
            raise ValueError('Actual positive source log factors must be retained')
    if not nonzero:raise ValueError('Full actual tensor was replaced by the zero local support difference')
    return len(groups),nonzero


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentAngularBackgroundStress(require_checked=False)
    field.assert_graph()
    if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k!='current_actual_angular_tensor_views'}:
        raise ValueError('Actual angular moment/source/pressure manifest differs')
    if any(raw[key] for key in GATES+OPEN):raise ValueError('Producer claims its own acceptance')
    if set(raw['current_actual_angular_tensor_views'])!=set(VIEWS):raise ValueError('Entire angular domain and fresh physical parameters required')
    count=nonzero=0
    for name,args in VIEWS.items():
        view=field.angular(*args)
        if encode(pack(view))!=raw['current_actual_angular_tensor_views'][name]:
            raise ValueError('Actual angular tensor source replay differs: '+name)
        n,m=check_view(field,view);count+=n;nonzero+=m
    rejected=[]
    for label,kwargs in (('outside_Z',dict(Z=2)),('outside_angular',dict(offset=-5)),
            ('nonfinite_source',dict(offset='-inf')),('nonfinite_time',dict(log_tau='-inf')),
            ('nonpositive_viscosity',dict(viscosity=0)),('nonfinite_angle',dict(theta='inf'))):
        args=dict(Z='.2',offset='-2',log_tau='-1',theta=None,viscosity='1');args.update(kwargs)
        try:field.angular(**args)
        except ValueError:rejected.append(label);continue
        raise ValueError('Nonphysical angular request accepted: '+label)
    for attribute in ('outer','pressure'):
        clone=copy.copy(field);setattr(clone,attribute,object())
        try:clone.assert_graph()
        except ValueError:rejected.append('different_current_'+attribute+'_owner');continue
        raise ValueError('Foreign current source owner accepted')
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_actual_angular_tensor_view_count=len(VIEWS),
        current_actual_angular_physical_tensor_divergence_remainder_rows_checked=count,
        current_actual_angular_nonzero_physical_contributions_checked=nonzero,
        original_full_defect_AST_normalized_mixed4_identity_count=len(field.moment_proof['original_full_defect_AST_normalized_mixed4_identities']),
        original_full_stress_AST_baseline_homogeneity_identity_count=len(field.baselines['actual_original_stress_mixed3_baseline_cancellation_and_KR_homogeneity']),
        source_theorems={key:raw[key] for key in ('current_full_moment_normalization_AST_theorem',
            'current_full_moment_baseline_and_normalization_AST_theorem','current_angular_transport_source_theorem',
            'current_original_absolute_pressure_future_source_theorem','original_general_K_physical_tensor_remainder_theorem',
            'exact_current_angular_source_log_cancellation')},
        invalid_domains_and_different_source_owners_rejected=rejected,
        same_current_analytic_P0_Pin_nonzero_histories_complete_Gamma_future_retained=True,
        native_interval_consistency_is_diagnostic_not_source_identity=True,
        actual_pressure_remaining_representation_from_replayed_current_Cp_zero=True,
        actual_tensor_uses_full_moments_after_exact_original_unit_baseline_cancellation=True,
        source_domain=raw['domain'],
        scope='Actual current angular full A/E/P/K stress mixed3, symmetric completed physical tensor, divergence and axial-viscosity remainder mixed2 over the original s[-4,0], whole source Z and finite compact positive-time sectors. This is a regional source-bound decomposition; other chart tensors, joins, global cone/flatness/NS/energy, resolved points and n-dependent temporal recursion remain open.',
        input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)},
        all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current actual angular full stress and regional completed physical tensor/remainder PASS; global/temporal gates remain open',flush=True)
    return result


if __name__=='__main__':run()
