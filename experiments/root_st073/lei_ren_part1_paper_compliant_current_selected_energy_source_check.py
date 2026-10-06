"""Current future/selection ownership, implicit equations and pulse consumption."""
import copy
import json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_selected_energy_source import (
    CurrentSelectedEnergySource,original_algorithm_source_bindings,HERE,PREFIX,NAME,RECEIPT,
    GATES,OPEN,VIEWS,PULSE_VIEWS,sha,binding,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_axial_high_jets_check import independent_root_fixture
from lei_ren_part1_paper_compliant_fifth_axial_jets_check import branch_fixtures
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def source_proof(field):
    inherited=original_algorithm_source_bindings()
    specs={
        'replay_complete_future':{'out.params':'exact.repair.params','out.mu':'exact.repair.mu',
            'out.delta':'exact.repair.delta','out.repair':'exact.repair',
            'out.Ntail':'out.Nt*c.exp(-out.delta*out.angular.waiting)',
            'out.waiting_energy':'decay_integral(c,out.delta,out.angular.waiting)',
            'out.tail_multiplier':'out.Ntail*c.exp(-2*out.angular.waiting_logone)',
            'out.weighted_factor':'out.mu*out.pulse_energy_attenuation/2'},
        '__init__':{'self.future':'replay_complete_future(self.exact)',
            'self.amplitude.future':'self.future','self.energy4.base':'self.future',
            'self.energy4.angular':'self.exact.angular4','self.energy4.cache':'{}',
            'self.axial4.energy':'self.energy4','self.axial4.base':'self.amplitude','self.axial4.cache':'{}',
            'self.fifth.fourth':'self.axial4','self.fifth.angular4':'self.exact.angular4',
            'self.fifth.angular_cache':'{}','self.fifth.energy_cache':'{}','self.fifth.cache':'{}',
            'self.pulse.fifth':'self.fifth','self.pulse.selection':'_SelectedSource(self.pulse.high)',
            'self.pulse.data_cache':'{}','self.amplitude.whole':'self.amplitude.select([-1,1])'}}
    for method,rows in specs.items():
        for target,value in rows.items():
            binding('compliant_current_selected_energy_source',method,target,value)
            inherited[method+':'+target]=True
    # The positive selected quadratic is precisely the forward energy
    # condition, including both nonzero formal end tails and incoming energy.
    K,a,nu1,nu2,u1,u2,v1,v2,mu,ein,F=s.symbols('K a nu1 nu2 u1 u2 v1 v2 mu ein F',real=True)
    base=(1-s.exp(-26))/4
    target=mu*s.exp(-26)*F/2-mu*ein+base
    selected=K*a*a+nu1*(u1+v1*a)**2+nu2*(u2+v2*a)**2
    forward=s.exp(26)*(ein+(selected-base)/mu)
    if s.simplify(forward.subs(selected,target)-F/2)!=0:
        raise ArithmeticError('Selected energy equation does not yield full future/2')
    # Equivalently use the selected equation as one scalar residual.
    if s.simplify(forward-F/2-s.exp(26)*(selected-target)/mu)!=0:
        raise ArithmeticError('Forward/full-future residual identity changed')
    return dict(actual_algorithm_and_adapter_AST_bindings=inherited,
        checked_defining_parameter_identity_consumed=field.exact.parameter_bridge['passed'],
        same_current_positive_quadratic_gives_forward_full_future_half=True,
        nonzero_formal_end_energy_and_incoming_terms_retained=True,
        original_source_integrals_not_replaced_by_enclosure_values=True,passed=True)


def rejected_owner_mutations(field):
    """Reject old future/selected objects even when numeric boxes overlap."""
    rejected=[]
    def reject(label,change):
        out=copy.copy(field);change(out)
        try:out.assert_graph()
        except ValueError:rejected.append(label);return
        raise ArithmeticError('Stale owner mutation accepted: '+label)
    def old_future(out):
        out.energy4=copy.copy(out.energy4);out.energy4.base=out.original_sources[3].base
    def old_selection(out):
        out.axial4=copy.copy(out.axial4);out.axial4.base=out.original_sources[2]
    def old_cache(out):
        out.fifth=copy.copy(out.fifth);out.fifth.cache=out.original_sources[0].cache
    def old_pulse(out):
        out.pulse=copy.copy(out.pulse);out.pulse.fifth=out.original_sources[0]
    for name,change in (('old_C4_future',old_future),('old_selected_amplitude',old_selection),
                        ('old_selected_cache',old_cache),('old_pulse_fifth',old_pulse)):
        reject(name,change)
    return dict(rejected=rejected,mutation_count=len(rejected),passed=True)


def contains_zero(value,label):
    lo,hi=endpoints(value)
    if not lo<=0<=hi:raise ArithmeticError('Current exact equation enclosure differs: '+label)


@source_precision
def run(field=None):
    raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
    field=field if field is not None else CurrentSelectedEnergySource(require_checked=False)
    if (raw['actual_five_defect_family_sha256'],raw['implicit_source_sha256'],raw['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
        raise ValueError('Current selected source/family/datum differs')
    if any(raw[k] for k in GATES+OPEN):raise ValueError('Current selected producer claims admission')
    proof=source_proof(field);mutations=rejected_owner_mutations(field)
    c4_fixture=independent_root_fixture();c5_fixture=branch_fixtures()
    future_rows=selected_rows=equation_rows=pulse_rows=0
    for name,Z in VIEWS.items():
        packet=field.evaluate(Z)
        if encode(pack(packet))!=raw['current_selected_energy_views'][name]:raise ValueError('Current selected packet differs: '+name)
        if not all(packet['source_owner_graph'].values()) or any(packet[k] for k in OPEN):raise ValueError('Owner graph or scope differs')
        e4=packet['current_future_C4']['complete_future_energy_Taylor']
        e5=packet['current_future_C5']['complete_future_energy_Taylor']
        if e4.order!=4 or e5.order!=5 or endpoints(e5[0])[0]<=0:raise ArithmeticError('Complete current energy order/positivity lost')
        for n in range(5):
            if endpoints(e4[n])!=endpoints(e5[n]):raise ValueError('Current energy C5 changed current C4 prefix')
        future_rows+=6
        fourth=packet['current_selected_C4'];fifth=packet['current_selected_C5'];c=field.ctx
        ap=fifth['selected_ap_Taylor'];controls=fifth['selected_scaled_end_coefficient_Taylor']
        quad=fifth['quadratic_coefficients'];residual=ap*ap*quad['A2']+ap*quad['A1']+quad['A0']
        for j,(before,after) in enumerate(zip(
                [fourth['selected_ap_Taylor']]+fourth['selected_scaled_end_coefficient_Taylor'],[ap]+controls)):
            if before.order!=4 or after.order!=5:raise ValueError('Current selected coefficient order lost')
            for n in range(5):
                if endpoints(before[n])!=endpoints(after[n]):raise ValueError('Current selected C5 changed C4 prefix')
            for value in after.coefficients:
                if not all(mp.isfinite(v) for v in endpoints(value)):raise ArithmeticError('Current selected derivative not finite')
                selected_rows+=1
        if not mp.mpf('.9')<endpoints(ap[0])[0]<=endpoints(ap[0])[1]<mp.mpf('1.2'):
            raise ArithmeticError('Current positive amplitude lost paper bracket')
        if endpoints(fifth['positive_root_derivative_denominator'])[0]<=0:raise ArithmeticError('Positive root inverse lost')
        if not endpoints(controls[0][0])[1]<0<endpoints(controls[1][0])[0]:raise ArithmeticError('Current end coefficient signs lost')
        for n in range(6):contains_zero(residual[n],'selected quadratic order'+str(n));equation_rows+=1
        # Same exact normalized native two-row system; each coefficient uses
        # the unchanged fixed factor, plus the current selected amplitude.
        A=field.amplitude.pulse.basis['first_row'];D=field.amplitude.pulse.basis['exact_divided_difference_row']
        Pi=field.amplitude.pulse.rows['scaled_full_rows'];Q=fifth['actual_row_normalized_incoming_Taylor']
        row1=controls[0]*A[0]+controls[1]*A[1]+ap*(field.amplitude.mu*Pi[0])+Q[0]*field.amplitude.mu
        row2=controls[0]*D[0]+controls[1]*D[1]+ap*(Pi[1]-Pi[0])+(Q[1]-Q[0])
        for jet in (row1,row2):
            for n in range(6):contains_zero(jet[n],'selected linear row');equation_rows+=1
    for name,args in PULSE_VIEWS.items():
        point=field.pulse_point(*args)
        if encode(pack(point))!=raw['current_rebound_pulse_views'][name]:raise ValueError('Rebound pulse packet differs: '+name)
        if not point['pulse_all_mixed_derivatives_total_order_le4_available']:raise ValueError('Current pulse mixed4 unavailable')
        for rows in point['physical_mixed_derivatives_total_order_le4'].values():
            if len(rows)!=15:raise ValueError('Current pulse mixed derivative order lost')
            for value in rows.values():
                if not all(mp.isfinite(v) for v in endpoints(value)):raise ArithmeticError('Current pulse derivative not finite')
                pulse_rows+=1
        if name=='terminal':
            for key in ('Uz_over_Utheta','Mz_over_R_Utheta',
                        'Mtheta_z_over_sqrt2_R_3half_Utheta_squared','Ur_over_sqrt_R_over_2_Utheta'):
                if any(endpoints(v)!=(mp.mpf(0),mp.mpf(0)) for v in point[key].coefficients):
                    raise ArithmeticError('Current meridional terminal history not zero: '+key)
            future=field.fifth.future(args[1])['complete_future_energy_Taylor']/2
            energy=point['Mztheta_over_R_Utheta_squared']
            if any(endpoints(energy[n])!=endpoints(future[n]) for n in range(6)):
                raise ValueError('Current terminal energy not same full future/2')
            if endpoints(energy[0])[0]<=0:raise ArithmeticError('Physical terminal energy was zeroed')
    hashes=dict(raw['input_hashes']);hashes[NAME]=sha(NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_future_selected_source_proof=proof,
        owner_mutation_rejections=mutations,independent_positive_root_C4_fixture=c4_fixture,
        independent_implicit_fifth_fixture=c5_fixture,
        current_future_C5_rows=future_rows,current_selected_C5_rows=selected_rows,
        current_quadratic_and_linear_equation_rows=equation_rows,current_pulse_mixed4_rows=pulse_rows,
        current_zero_meridional_terminal_histories=True,current_positive_terminal_full_future_half=True,
        all_passed=True,input_hashes=hashes,**dict.fromkeys(GATES,True),**dict.fromkeys(OPEN,False))
    (HERE/RECEIPT).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Current future/selection source and restricted pulse mixed4 PASS: '+str(selected_rows)+' selected rows, '+str(pulse_rows)+' pulse rows',flush=True)
    return result


if __name__=='__main__':run()
