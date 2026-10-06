"""Actual current heat pressure and stress, retaining terminal constants.

Forward history is rewritten using the SAME full future integral. Neither
the pressure at infinity nor the angular integration constant is set to
zero. Positive stress factors remain source-defined, not cap values.
"""
import ast
import json
import math
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_core_first_interface import (
    CurrentCoreFirstInterface,HERE,PREFIX,sha,function,binding,pack)
from lei_ren_part1_paper_compliant_current_heat_source import current_source_bindings
from lei_ren_part1_paper_compliant_collar_pressure_C4 import collar_pressure_rows
from lei_ren_part1_paper_compliant_collar_stress_C3 import (
    axial_derivative,collar_defect_rows,collar_stress_rows)
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_heat_pressure_stress.json'
RECEIPT=PREFIX+'current_heat_pressure_stress_check.json'
GATES=('current_heat_forward_absolute_pressure_mixed4_companion_certified',
       'current_collar_actual_stress_mixed3_recovered',
       'current_exterior_actual_stress_mixed4_recovered')
SCOPES=('current_heat_terminal_constants_eliminated','heat_exterior_stress_identity_certified',
    'global_completed_tensor_admissibility','admissible_stress_lift_constructed',
    'full_background_NS_validation','physical_energy_integral_certified',
    'independently_bounded_flat_remainder','full_point_physical_field_evaluation',
    'full_cartesian_vector_derivatives_certified','temporal_recursion')
VIEWS={'whole_collar':('heat_collar',[0,3]),'collar_Gamma_join':('heat_collar',3),
       'whole_exterior':('heat_exterior',[3,mp.inf]),'fresh_exterior':('heat_exterior',4)}


def augmented_binding(stem,method,target,expression):
    rows=[node.value for node in ast.walk(function(stem,method)) if isinstance(node,ast.AugAssign)
          and isinstance(node.op,ast.Add) and ast.unparse(node.target)==target]
    wanted=ast.dump(ast.parse(expression,mode='eval').body)
    if len(rows)!=1 or ast.dump(rows[0])!=wanted:raise ValueError('Original pressure density changed')


def history_source_bindings():
    specs={
        'data':{'terminal':'self.steep.waiting(Z,1)','heat0':'self.collar_tails(Z,0)',
            'Xtail':"terminal['angular_Taylor']",'Ptail':"terminal['pressure_over_Pstar_squared_Taylor']",
            'defect':"Xtail*(1-self.eps)-heat0['angular_numerator']",
            'pressure3':'self.forward_pressure(Z,3,Ptail)'},
        'forward_pressure':{'K':"self.shape(Z,v,False)['K_rows'][0]"},
        'collar':{'X':"(tails['angular_numerator']+data['angular_tail_constant_defect']*c.exp(-self.k*t))/K",
            'energy':"tails['remaining_energy_in_Rtail_units']*c.exp(self.delta*t)/(K*K*2)",
            'pressure':"self.forward_pressure(Z,t,data['Ptail'])"},
        'exterior':{'X':"(local['angular_numerator']+data['angular_tail_constant_defect']*c.exp(-self.k*t))/K",
            'energy':"local['energy_numerator']/(K*K*2)",
            'pressure':"data['pressure3']+(loc3['pressure_numerator']*c.exp(-3*self.prate)-local['pressure_numerator']*c.exp(-self.prate*t))*self.pressure_scale"},
        'collar_tails':{'square':'D*(pre*2-D*(self.a*self.S))',
            'P':"one*(c.exp(-self.prate*t)/(2*self.prate)-self.eps*atoms['PW']+self.eps**2*atoms['PW2']/2)-pressure*(self.a*self.S)"},
        'shape':{'W':'product_rows(sr,fr)','pre':'[unity[j]-W[j]*self.eps for j in range(5)]',
            'K':'[pre[j]-D[j]*(self.a*self.S) for j in range(len(D))]'}}
    # W has an additional (1-sigma) contribution in the actual shape.
    specs['shape'].pop('W')
    for method,assignments in specs.items():
        for target,expression in assignments.items():binding('compliant_collar_Gamma_C4',method,target,expression)
    augmented_binding('compliant_collar_Gamma_C4','forward_pressure','integral',
        'K*K*(ds*c.exp(-self.prate*v)/2)')
    fn=function('compliant_collar_Gamma_C4','forward_pressure')
    returns=[node for node in fn.body if isinstance(node,ast.Return)]
    if len(returns)!=1 or ast.dump(returns[0].value)!=ast.dump(ast.parse('Ptail+integral*self.pressure_scale',mode='eval').body):
        raise ValueError('Actual forward absolute pressure changed')
    return dict(actual_current_waiting_terminal_and_both_constants_AST_bound=True,
        same_full_collar_pressure_density_and_Gamma_tail=True,
        actual_angular_and_absolute_pressure_histories_not_reset=True,
        actual_energy_full_future_half_normalization_retained=True,passed=True)


def mixed(rows,order):
    return {'y%d_Z%d'%(j,n):rows[j][n]*math.factorial(n)
            for j in range(order+1) for n in range(order+1-j)}


def constant_stress_rows(heat,Z,t,angular_defect,pressure_infinity,order):
    """Original moment-stress change due to two retained constants.

    Theta uses Qtheta=sqrt(R/2)*B. Pressure uses the distinct positive
    Qpressure=sqrt(R/2)*Pstar^2; division by an Ev2 box is unnecessary.
    Rows include the ordinary derivatives of these source factors.
    """
    c=heat.ctx;z=IntervalTaylor.variable(c,Z,5);d=1-z*z;L=1-z*z*heat.delta
    b=(1-heat.delta)/2
    theta=(angular_defect*heat.k-z*axial_derivative(angular_defect)*b)/L
    pressure=(z*pressure_infinity*(2*heat.prate)-d*axial_derivative(pressure_infinity))/L
    decay=c.exp(-heat.k*t)
    return dict(theta=[theta*decay*(-1)**j for j in range(order+1)],
        axial_pressure_constant=[pressure*c.mpf('.5')**j for j in range(order+1)])


class CurrentHeatPressureStress:
    @source_precision
    def __init__(self,joined=None,require_checked=True):
        self.joined=joined if joined is not None else CurrentCoreFirstInterface()
        if not self.joined.acceptance_loaded:raise ValueError('Checked current common core and four bridge joins required')
        self.owner=self.joined.common.source.bridge.owner
        self.source_owner=self.owner.source_assembly;self.heat=self.source_owner.heat
        self.ctx=self.heat.ctx;self.family=self.joined.family;self.source=self.joined.source_sha;self.datum_sha=self.joined.datum_sha
        self.graph=dict(same_current_heat_and_physical_owner=self.heat is self.owner.heat,
            checked_current_heat_source=self.source_owner.acceptance_loaded,
            common_core_is_same_current_physical_owner=self.joined.core is self.owner.core.original,
            same_current_source_datum=(self.family,self.source,self.datum_sha)==
                (self.source_owner.family,self.source_owner.source,self.source_owner.datum_sha),
            same_current_waiting=self.heat.steep is self.source_owner.before.steep,
            same_current_future=self.heat.future is self.heat.steep.future,
            same_selected_raw_heat=self.heat.exact_heat is self.heat.future.heat,
            same_original_native_context=self.ctx is self.heat.steep.ctx is self.heat.outer.ctx)
        if not all(self.graph.values()):raise ValueError('Heat companion must use the common current source graph')
        self.current_bindings=current_source_bindings(self.source_owner)
        self.history_bindings=history_source_bindings();self.hashes=dict(self.joined.hashes)
        for stem in ('current_heat_source','collar_Gamma_C4','collar_pressure_C4','heat_pressure_C4',
            'collar_stress_C3','heat_terminal_stress_identities','heat_stress_equations','exact_heat_component'):
            name=PREFIX+stem+'.py';self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.constants={};self.runtime={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in SCOPES):
                raise ValueError('Heat companion admission source or scope differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def terminal_constants(self,Z):
        c=self.ctx;z=c.mpf(Z);key=z._mpi_
        if key not in self.constants:
            data=self.heat.data(z);tail0=self.heat.collar_tails(z,0);gamma3=self.heat.local_Gamma(z,3)
            pressure_infinity=data['Ptail']+tail0['remaining_pressure_in_Rtail_units']*self.heat.pressure_scale
            from_trace3=data['pressure3']+gamma3['pressure_numerator']*(self.heat.pressure_scale*c.exp(-3*self.heat.prate))
            self.constants[key]=dict(angular_defect=data['angular_tail_constant_defect'],
                pressure_infinity=pressure_infinity,pressure_infinity_from_actual_trace3=from_trace3,
                actual_current_waiting_Ptail=data['Ptail'],actual_current_waiting_Xtail=data['Xtail'],
                full_collar_future_pressure=tail0['remaining_pressure_in_Rtail_units'],
                full_collar_future_angular=tail0['angular_numerator'],
                zero_constants_not_assumed=True)
        return self.constants[key]

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        if chart not in ('heat_collar','heat_exterior'):raise ValueError('Original heat chart required')
        c=self.ctx;z=c.mpf(Z);t=c.mpf(coordinate);lo,hi=endpoints(t)
        if endpoints(z)[0]<-1 or endpoints(z)[1]>1 or lo<0 or (chart=='heat_collar' and hi>3) or (chart=='heat_exterior' and lo<3):
            raise ValueError('Original whole-Z collar/exterior domain required')
        key=(chart,z._mpi_,t._mpi_)
        if key in self.runtime:return self.runtime[key]
        constants=self.terminal_constants(z);native=self.source_owner.evaluate(chart,z,t)['source_packet']
        if chart=='heat_collar':
            shape=self.heat.shape(z,t);tails=self.heat.collar_tails(z,t)
            remaining=tails['remaining_pressure_in_Rtail_units']
            defects=collar_defect_rows(self.heat,shape,tails,t)
            canonical=collar_stress_rows(self.heat,shape,defects,z,t);order=3
        else:
            shape=self.heat.local_Gamma(z,t)
            remaining=shape['pressure_numerator']*c.exp(-self.heat.prate*t);order=4
            # The projected full Gamma moments cancel independently of the
            # actual two constants. Those actual constants are retained below.
            zero=IntervalTaylor.constant(c,0,5)
            canonical=dict(theta=[zero]*5,axial=[zero]*5)
        pressure=collar_pressure_rows(shape['K_rows'],remaining,self.heat.prate,self.heat.pressure_scale,t)
        pressure[0]+=constants['pressure_infinity']
        extra=constant_stress_rows(self.heat,z,t,constants['angular_defect'],constants['pressure_infinity'],order)
        theta=[canonical['theta'][j]+extra['theta'][j] for j in range(order+1)]
        stress=dict(theta_Qtheta=mixed(theta,order),axial_Qz=mixed(canonical['axial'],order),
            axial_pressure_Qpressure=mixed(extra['axial_pressure_constant'],order))
        result=dict(chart=chart,Z=z,coordinate=t,
            actual_current_terminal_constants=constants,
            original_forward_pressure=native['pressure_over_Pstar_squared_Taylor'],
            original_forward_angular=native['angular_Taylor'],
            stable_same_forward_pressure_rows=pressure,
            stable_same_forward_pressure_mixed4=mixed(pressure,4),
            original_current_velocity_pressure_mixed4=native['physical_mixed_derivatives_total_order_le4'],
            actual_stress_factored_mixed_rows=stress,actual_stress_mixed_order=order,
            positive_source_stress_factors=dict(theta_Qtheta='sqrt(R/2)*B',axial_Qz='sqrt(R/2)*B^2',
                axial_pressure_Qpressure='sqrt(R/2)*Pstar^2',
                B='Ev0*theta_base*exp(-(1+delta)*t/2)',
                Ev0='Pstar*U*exp(-13/(2mu)-13)',R='Rtail*exp(t)',
                exact_logRtail_terms=self.heat.exact_heat.logradius_terms,
                exact_Ev0_squared_over_Pstar_squared_logs=self.heat.outer.flatten.logEv2_parts,
                physical_tensor_prefactor='nu*lambda^(-2-delta); transfer is a later global obligation'),
            actual_angular_constant_contribution_retained=True,
            actual_pressure_infinity_contribution_retained=True,
            exterior_theta_constant_stress_radial_power=-1,
            exterior_pressure_constant_stress_radial_power=c.mpf('.5'),
            full_collar_and_infinite_Gamma_integrals_used=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(SCOPES,False))
        self.runtime[key]=result;return result


@source_precision
def run(field=None):
    field=field if field is not None else CurrentHeatPressureStress(require_checked=False)
    views={name:field.evaluate(chart,[-1,1] if name!='fresh_exterior' else '.381',coordinate)
           for name,(chart,coordinate) in VIEWS.items()}
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,current_common_heat_source_graph=field.graph,
        current_heat_source_bindings=field.current_bindings,retained_history_bindings=field.history_bindings,
        actual_current_heat_companion_views=views,
        remaining_dependency='Eliminate actual current angular/pressure constants by original selected-source equations; extend collar shear to mixed4 before global stress/cone/flatness/energy and temporal recursion',
        **dict.fromkeys(GATES,False),**dict.fromkeys(SCOPES,False),input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(encode(pack(result)),indent=2)+'\n').encode('utf8'))
    print('Built actual current heat pressure and factored stress with both original constants retained',flush=True)
    return result


if __name__=='__main__':run()
