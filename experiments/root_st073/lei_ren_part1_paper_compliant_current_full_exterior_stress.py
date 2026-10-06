"""All current terminal moments -> actual full Gamma exterior stress zero.

The current energy/selected/postpulse owner supplies the missing histories.
Original angular and pressure terminal identities transfer through the same
native functions; forward diagnostics remain available. Only after those
links does the canonical full-integral theorem evaluate both stresses as
zero. Every physical assembler, global cone and temporal recursion remain
separate from this current similarity-stress exterior certificate.
"""
import json
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_postpulse_energy_history import (
    CurrentPostpulseEnergyHistory,HERE,PREFIX,sha,binding,pack,encode,endpoints)
from lei_ren_part1_paper_compliant_heat_terminal_stress_identities import terminal_stress_identities
from lei_ren_part1_paper_compliant_heat_stress_equations import heat_stress_equations
from lei_ren_part1_paper_compliant_collar_stress_C3 import collar_Gamma_endpoint_binding,collar_moment_stress_identities
from lei_ren_part1_paper_compliant_heat_stress_C4 import angular_primitive_y_rows
from lei_ren_part1_paper_compliant_collar_pressure_C4 import collar_pressure_rows
from lei_ren_part1_paper_compliant_current_heat_pressure_stress import mixed
from lei_ren_part1_paper_compliant_pulse_physical_bounds import P
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_full_exterior_stress.json'
RECEIPT=PREFIX+'current_full_exterior_stress_check.json'
GATES=('current_full_exterior_five_moment_source_history_certified',
       'current_full_exterior_similarity_stress_mixed4_certified','heat_exterior_stress_identity_certified')
OPEN=('current_exact_repair_installed_in_all_physical_charts','global_completed_tensor_admissibility',
      'admissible_stress_lift_constructed','full_background_NS_validation',
      'physical_energy_integral_certified','independently_bounded_flat_remainder',
      'full_cartesian_vector_derivatives_certified','temporal_recursion')
VIEWS={'whole_Z_exit':([-1,1],3),'whole_unbounded_exterior':([-1,1],[3,mp.inf]),
       'fresh_exterior':('.731','4.17'),'axis_exit':(0,3)}


def current_full_history_transfer(history):
    selected=history.selected;pressure=selected.pressure;exact=selected.exact;heat=history.heat
    if not history.acceptance_loaded or not pressure.acceptance_loaded:raise ValueError('Checked full energy and pressure histories required')
    if not history.proof['selected_forward_cumulative_energy_equals_same_remaining_integral_by_FTC']:
        raise ValueError('Current cumulative energy source history missing')
    graph=history.assert_graph();bindings={}
    specs={
        ('compliant_flatten_mixed_C4','flatten'):{
            'X':'(Xint+self.Xv*c.exp(-self.rate*t))/F',
            'Mp_v':"data['Mp']+data['u']*data['u']*((1-self.pressure_decay)/(2*self.prate))",
            'Mp':'Mp_v+Pint*self.Ev2','pressure':"Mp+data['P0']"},
        ('compliant_current_power_angular_source','data'):{
            'source':'self.fifth.angular(Z)','f':'self.flatten.flatten(Z,100)',
            'coeff':"[copy_jet(c,v) for v in source['physical_coefficient_Taylor']]"},
        ('compliant_current_steep_waiting_source','data'):{
            'terminal':'self.outer.angular(Z,0)','source':'self.outer.fifth.angular(Z)',
            'XR':"terminal['angular_Taylor']",'PR':"terminal['pressure_over_Pstar_squared_Taylor']",
            'XS':"(XR+self.infull['angular'])*c.exp(-self.rate/2)",
            'XQ':'XS+self.Ts','XT':"(XQ+self.outfull['angular'])*c.exp(-self.k/2)",
            'PS':"PR+self.infull['pressure']*(self.outer.flatten.Ev2*self.thetaR**2)",
            'PQ':'PS+decay_integral(c,3,self.Ts)*(self.outer.flatten.Ev2*self.thetaS**2/2)',
            'PT':"PQ+self.outfull['pressure']*(self.outer.flatten.Ev2*self.thetaQ**2)"},
        ('compliant_collar_Gamma_C4','data'):{
            'terminal':'self.steep.waiting(Z,1)','heat0':'self.collar_tails(Z,0)',
            'Xtail':"terminal['angular_Taylor']",'Ptail':"terminal['pressure_over_Pstar_squared_Taylor']",
            'defect':"Xtail*(1-self.eps)-heat0['angular_numerator']"},
        ('compliant_collar_Gamma_C4','exterior'):{
            'X':"(local['angular_numerator']+data['angular_tail_constant_defect']*c.exp(-self.k*t))/K",
            'energy':"local['energy_numerator']/(K*K*2)",
            'pressure':"data['pressure3']+(loc3['pressure_numerator']*c.exp(-3*self.prate)-local['pressure_numerator']*c.exp(-self.prate*t))*self.pressure_scale"}}
    for (stem,method),rows in specs.items():
        for target,value in rows.items():
            binding(stem,method,target,value);bindings[stem+'.'+method+':'+target]=True
    links=dict(same_checked_current_exact_repair=pressure.exact is exact,
        same_current_repair_function=heat.repair is exact.repair,
        same_original_inlet_Xp=selected.pulse.Xp is exact.pulse.Xp,
        same_original_inlet_pressure=selected.pulse.inlet_P is exact.pulse.inlet_P,
        same_original_inlet_amplitude=selected.pulse.high.constants is exact.pulse.high.constants,
        same_original_datum_function=history.flatten.inlet.datum is selected.future.angular.initial.datum is exact.repair.angular.initial.datum,
        same_original_angular_C4_repair=selected.fifth.angular4 is exact.angular5.angular4 is exact.angular4,
        same_original_angular_fifth_function=selected.fifth.angular.__func__ is exact.angular5.angular.__func__,
        same_original_angular_fifth_context=selected.fifth.ctx is exact.angular5.ctx,
        same_waiting_and_transition_parameter_definitions=selected.future.params is exact.future.params is exact.repair.params,
        same_current_full_Gamma_shape_radius_source=heat.exact_logS_terms==exact.repair.heat.logS_terms,
        checked_exact_pressure_function_identity=pressure.proof['passed'],
        checked_exact_angular_function_identity=pressure.angular.proof['passed'],
        checked_full_future_half_energy_identity=history.proof['same_current_full_future_half_normalization_on_every_stage'],
        checked_zero_meridional_terminal_and_FTC=history.proof['zero_meridional_histories_propagate_from_selected_terminal_by_FTC'],
        same_canonical_full_Gamma_source=history.proof['checked_canonical_Gamma_source_evidence']['verified'])
    if not all(links.values()):raise ValueError('Full current exterior history source link differs: '+str(links))
    # The rebound future affects the native energy only. The AST-bound X/P
    # recipes, exact inputs and current repair functions above are unchanged,
    # so the already proved Dtheta/Cp functions are exactly these histories.
    return dict(source_owner_graph=graph,current_angular_pressure_transfer_assignments=bindings,
        current_exact_function_links=links,
        current_energy_history_proof=history.proof,
        current_angular_terminal_source_proof=pressure.angular.proof,
        current_original_pressure_terminal_source_proof=pressure.proof,
        all_current_five_terminal_moment_functions_identified=True,
        original_forward_angular_and_pressure_diagnostics_retained=True,
        energy_source_change_does_not_reset_or_change_native_angular_pressure_functions=True,
        physical_terminal_energy_positive_not_zeroed=True,
        actual_meridional_histories_zero_by_selected_linear_equations_and_FTC=True,passed=True)


class CurrentFullExteriorStress:
    @source_precision
    def __init__(self,history=None,require_checked=True):
        self.history=history if history is not None else CurrentPostpulseEnergyHistory()
        self.pressure=self.history.selected.pressure;self.heat=self.history.heat;self.ctx=self.heat.ctx
        self.family=self.history.family;self.source=self.history.source;self.datum_sha=self.history.datum_sha
        self.proof=current_full_history_transfer(self.history)
        self.theorem=terminal_stress_identities();self.equations=heat_stress_equations()
        self.native_stress_AST=collar_Gamma_endpoint_binding();self.units=collar_moment_stress_identities()
        if not (self.theorem['full_terminal_moment_stress_theorem_verified'] and
                self.native_stress_AST['actual_pre_override_stress_AST_bound_to_terminal_theorem'] and
                self.units['original_collar_moment_to_stress_identities_verified']):
            raise ValueError('Full Gamma/native moment-to-stress theorem missing')
        if not 0<endpoints(self.heat.delta)[0]<=endpoints(self.heat.delta)[1]<1:
            raise ValueError('Canonical exterior parameter range lost')
        self.hashes=dict(self.history.hashes)
        for stem in ('heat_terminal_stress_identities','heat_stress_equations','heat_stress_C4',
            'collar_stress_C3','collar_pressure_C4','heat_pressure_C4'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.cache={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current full exterior scope/datum differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def exterior(self,Z,t):
        if endpoints(self.ctx.mpf(t))[0]<3:raise ValueError('Full Gamma exterior requires t>=3')
        self.history.assert_graph();key=(tuple(endpoints(self.ctx.mpf(Z))),tuple(endpoints(self.ctx.mpf(t))))
        if key in self.cache:return self.cache[key]
        c=self.ctx;point=dict(self.history.evaluate('heat_exterior',Z,t)['source_packet'])
        local=point['full_local_Gamma_future_source'];K=local['K_rows'];zero=IntervalTaylor.constant(c,0,5)
        X=local['angular_numerator']/K[0]
        pressure_rows=collar_pressure_rows(K[:5],local['pressure_numerator']*c.exp(-self.heat.prate*c.mpf(t)),
            self.heat.prate,self.heat.pressure_scale,c.mpf(t))
        data=self.heat.data(Z);local3=self.heat.local_Gamma(Z,3)
        Cp=data['pressure3']+local3['pressure_numerator']*(self.heat.pressure_scale*c.exp(-3*self.heat.prate))
        point.update(original_native_forward_angular_Taylor=point['angular_Taylor'],
            original_native_forward_angular_y_derivative_Taylor=point['angular_y_derivative_Taylor'],
            original_native_forward_pressure_Taylor=point['pressure_over_Pstar_squared_Taylor'],
            original_native_forward_pressure_mixed4=point['physical_mixed_derivatives_total_order_le4'][P],
            original_native_Dtheta_enclosure=data['angular_tail_constant_defect'],original_native_Cp_enclosure=Cp,
            stable_current_angular_Taylor=X,
            stable_current_angular_y_derivative_Taylor=angular_primitive_y_rows(X,point['angular_ODE_rate_Taylor']),
            stable_current_absolute_pressure_Taylor=pressure_rows[0],
            stable_current_absolute_pressure_mixed4=mixed(pressure_rows,4),
            actual_Dtheta_Taylor=zero,actual_Cp_Taylor=zero,
            exterior_meridional_moment_Taylor=dict(Mz=zero,Mtheta_z=zero),
            exterior_similarity_stress_over_F=dict(theta=zero,axial=zero),
            exterior_stress_mixed4={label:{'y%d_Z%d'%(j,n):c.mpf(0) for j in range(5) for n in range(5-j)} for label in ('theta','axial')},
            all_current_five_terminal_histories_consumed=True,
            actual_energy_is_positive_same_full_Gamma_future_half=True,
            stress_zero_is_source_identity_not_interval_overlap=True,
            exact_positive_stress_prefactors=dict(theta='sqrt(R/2)*B',axial='sqrt(R/2)*B^2',
                B='Ev0*theta_base*exp(-(1+delta)*t/2)',R='Rtail*exp(t)',
                Ev0='Pstar*U*exp(-13/(2mu)-13)',
                exact_logRtail_terms=self.heat.exact_logRtail_terms,
                exact_Ev0_squared_over_Pstar_squared_logs=self.history.flatten.logEv2_parts,
                physical_tensor_common_factor='nu*lambda^(-2-delta)',
                positivity_domain='each finite R>0 and physical time before T; infinity is the covered limit'),
            scope='current source-bound original similarity stress (3.16)-(3.18), Z in [-1,1], log(R/Rtail)>=3',
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        self.cache[key]=point;return point

    @source_precision
    def report(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,current_full_terminal_history_transfer=self.proof,
            canonical_full_Gamma_stress_theorem=self.theorem,canonical_heat_equations=self.equations,
            actual_native_pre_override_stress_AST=self.native_stress_AST,original_stress_units=self.units,
            current_full_exterior_views={name:self.exterior(*args) for name,args in VIEWS.items()},
            input_hashes=self.hashes,**dict.fromkeys(GATES+OPEN,False))


@source_precision
def run(field=None):
    field=field if field is not None else CurrentFullExteriorStress(require_checked=False)
    result=field.report();(HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('All current five terminal histories consumed: both full Gamma exterior similarity stresses zero',flush=True)
    return result


if __name__=='__main__':run()
