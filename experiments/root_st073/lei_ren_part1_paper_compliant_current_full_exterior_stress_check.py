"""Source-bound current full exterior stress admission, not global admission."""
import copy
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_full_exterior_stress import (
    CurrentFullExteriorStress,current_full_history_transfer,HERE,PREFIX,NAME,RECEIPT,
    GATES,OPEN,VIEWS,sha,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_current_postpulse_energy_history import energy_transport_source_proof
from lei_ren_part1_paper_compliant_heat_terminal_stress_identities import terminal_stress_identities
from lei_ren_part1_paper_compliant_heat_stress_equations import heat_stress_equations
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_Gamma_endpoint_binding,collar_moment_stress_identities
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def rejected_history_mutations(field):
    rejected=[]
    for label,change in (
        ('unchecked_energy_history',lambda h:setattr(h,'acceptance_loaded',False)),
        ('missing_energy_source_integral',lambda h:h.proof.update(selected_forward_cumulative_energy_equals_same_remaining_integral_by_FTC=False)),
        ('missing_zero_meridional_source',lambda h:h.proof.update(zero_meridional_histories_propagate_from_selected_terminal_by_FTC=False)),
        ('unchecked_pressure_source',lambda h:setattr(h.selected.pressure,'acceptance_loaded',False))):
        h=copy.copy(field.history);h.proof=dict(h.proof);h.selected=copy.copy(h.selected)
        h.selected.pressure=copy.copy(h.selected.pressure)
        change(h)
        try:current_full_history_transfer(h)
        except ValueError:rejected.append(label);continue
        raise ArithmeticError('Unproved full exterior history accepted: '+label)
    return dict(rejected=rejected,mutation_count=len(rejected),passed=True)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentFullExteriorStress(require_checked=False)
    if (raw['actual_five_defect_family_sha256'],raw['implicit_source_sha256'],raw['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
        raise ValueError('Full exterior source/family/datum differs')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Full exterior producer claims admission')
    # Recompute the actual current function bridge before invoking the
    # generic canonical theorem. Numeric zero boxes never supply this link.
    energy=energy_transport_source_proof(field.history)
    if encode(pack(energy))!=encode(pack(field.history.proof)):
        raise ValueError('Current energy/meridional source proof changed')
    proofs={
        'current_full_terminal_history_transfer':current_full_history_transfer(field.history),
        'canonical_full_Gamma_stress_theorem':terminal_stress_identities(),
        'canonical_heat_equations':heat_stress_equations(),
        'actual_native_pre_override_stress_AST':collar_Gamma_endpoint_binding(),
        'original_stress_units':collar_moment_stress_identities()}
    for key,value in proofs.items():
        if encode(pack(value))!=raw[key]:raise ValueError('Full exterior source/theorem changed: '+key)
    mutations=rejected_history_mutations(field)
    indices={'y%d_Z%d'%(j,n) for j in range(5) for n in range(5-j)}
    stress_rows=constant_rows=pressure_rows=0
    for name,args in VIEWS.items():
        point=field.exterior(*args)
        if encode(pack(point))!=raw['current_full_exterior_views'][name]:raise ValueError('Current full exterior packet changed: '+name)
        if any(point[k] for k in OPEN):raise ValueError('Full exterior packet exceeds scope')
        if endpoints(point['energy_Taylor'][0])[0]<=0:raise ArithmeticError('Terminal energy was zeroed')
        for key in ('original_native_Dtheta_enclosure','original_native_Cp_enclosure'):
            for v in point[key].coefficients:
                lo,hi=endpoints(v)
                if not lo<=0<=hi:raise ArithmeticError('Proved exact constant contradicts native enclosure: '+key)
        for key in ('actual_Dtheta_Taylor','actual_Cp_Taylor'):
            for v in point[key].coefficients:
                if endpoints(v)!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Terminal constant zero identity lost')
                constant_rows+=1
        for key in ('Mz','Mtheta_z'):
            if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in point['exterior_meridional_moment_Taylor'][key].coefficients):
                raise ArithmeticError('Actual zero meridional history lost')
        for rows in point['exterior_stress_mixed4'].values():
            if set(rows)!=indices:raise ValueError('Full exterior stress derivative omitted')
            for v in rows.values():
                if endpoints(v)!=(mp.mpf(0),mp.mpf(0)):raise ArithmeticError('Full exterior stress identity lost')
                stress_rows+=1
        rows=point['stable_current_absolute_pressure_mixed4']
        if set(rows)!=indices:raise ValueError('Current stable pressure mixed4 omitted')
        for v in rows.values():
            if not all(mp.isfinite(x) for x in endpoints(v)):raise ArithmeticError('Nonfinite exact pressure derivative')
            pressure_rows+=1
        # The raw native packet remains intact, including its original
        # pressure grids; the equivalent stable representation is separate.
        native=field.history.evaluate('heat_exterior',*args)['source_packet']
        for key,value in native.items():
            if encode(pack(point[key]))!=encode(pack(value)):
                raise ValueError('Full stress wrapper changed original native field: '+key)
        for native_key,stable_key in (('angular_Taylor','stable_current_angular_Taylor'),
            ('pressure_over_Pstar_squared_Taylor','stable_current_absolute_pressure_Taylor')):
            for a,b in zip(point[native_key].coefficients,point[stable_key].coefficients):
                al,ah=endpoints(a);bl,bh=endpoints(b)
                if max(al,bl)>min(ah,bh):raise ArithmeticError('Exact canonical transfer contradicts native diagnostic')
        if not point['all_current_five_terminal_histories_consumed'] or not point['stress_zero_is_source_identity_not_interval_overlap']:
            raise ValueError('Full exterior source history omitted')
    hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)}
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,**proofs,
        missing_history_mutations=mutations,full_exterior_stress_mixed4_rows=stress_rows,
        source_proved_zero_terminal_constant_C5_rows=constant_rows,
        stable_current_absolute_pressure_mixed4_rows=pressure_rows,
        original_native_field_and_datum_retained=True,
        current_positive_full_energy_half_retained=True,
        zero_stress_scope='source-bound original similarity tensor theta/axial exterior components; every Z in [-1,1], log(R/Rtail)>=3',
        input_hashes=hashes,all_passed=True,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current full exterior stress PASS: '+str(stress_rows)+' zero stress mixed4 rows; physical/global/time remain open',flush=True)
    return result


if __name__=='__main__':run()
