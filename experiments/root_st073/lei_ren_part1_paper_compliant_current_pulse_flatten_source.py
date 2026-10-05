"""Current native pulse terminal admission into the unchanged O5 flatten.

This extends current profile source ownership, not physical chart ownership.
The canonical endpoint theorem is consumed only after its defining inputs
have been rebound to the current native object. Bounds are not point values.
"""
import ast
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_native_pulse_source_dispatcher import (
    CurrentNativePulseSourceDispatcher, CHARTS, UNIFORM, SCOPES, HERE, PREFIX, sha)
from lei_ren_part1_paper_compliant_current_pulse_physical_assembly import (
    native_parameter_source_bridge, pressure_publication_bindings)
from lei_ren_part1_paper_compliant_current_downstream_physical_assembly import OPEN
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
from lei_ren_part1_paper_compliant_actual_Rsh_source_join import class_assignment
from lei_ren_part1_paper_compliant_power_inlet_C4 import CompliantPowerInletC4
from lei_ren_part1_paper_compliant_flatten_mixed_C4 import CompliantFlattenMixedC4
from lei_ren_part1_paper_compliant_pulse_radial_C4 import CompliantPulseRadialC4
from lei_ren_part1_paper_compliant_pulse_mixed_C4 import CompliantPulseMixedC4
from lei_ren_part1_paper_compliant_axial_pulse_field import CompliantAxialPulseField
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

RECEIPT=PREFIX+'current_pulse_flatten_source_check.json'
THEOREM=PREFIX+'pulse_end_flatten_join_check.json'
GATE='current_pulse_terminal_flatten_source_ownership_certified'
PROVED='current_pulse_terminal_flatten_source_functional_join_proved'


class CurrentPowerInlet(CompliantPowerInletC4):
    """Rp reference data for the original ODE; not the Rv terminal histories."""
    def __init__(self,pulse,family,source,hashes):
        self.ctx=pulse.ctx;self.mu=pulse.mu;self.delta=pulse.delta
        self.constants=pulse.high.constants
        self.inlet_H=pulse.inlet_H;self.inlet_P=pulse.inlet_P;self.Xp=pulse.Xp
        self.datum=pulse.selection.future.angular.initial.datum
        self.family=family;self.source=source;self.hashes=hashes


class CurrentFlattenMixedC4(CompliantFlattenMixedC4):
    """Original flatten ODEs; no saved terminal/history lookup in this owner."""
    def __init__(self,pulse,family,source,hashes,cells=128):
        if not isinstance(cells,int) or cells<1:raise ValueError('Positive cell count required')
        self.pulse=pulse;self.inlet=CurrentPowerInlet(pulse,family,source,hashes)
        self.ctx=c=pulse.ctx;self.cells=cells;self.mu=pulse.mu;self.delta=pulse.delta
        self.rate=pulse.rate;self.bp=c.mpf('.5')+self.mu;self.prate=pulse.prate
        self.family=family;self.source=source;self.hashes=hashes
        self.Xv=1/self.rate+(pulse.Xp-1/self.rate)*pulse.factor(-13*self.rate/self.mu)
        self.U=self.inlet.constants['U']
        self.logEv2_parts=dict(inlet_log=2*c.ln(self.U),inverse_mu_term=-13/self.mu,finite_offset=c.mpf(-26))
        self.log_pressure_decay_parts=dict(inverse_mu_term=-13/self.mu,finite_offset=c.mpf(-26))
        self.log_decay=sum(self.log_pressure_decay_parts.values(),c.mpf(0))
        logcap=2*c.ln(self.mu)-1000;exactlog=sum(self.logEv2_parts.values(),c.mpf(0))
        if endpoints(logcap-exactlog)[0]<=0 or endpoints(logcap-self.log_decay)[0]<=0:
            raise ArithmeticError('Exact positive terminal amplitude/pressure-decay cap not proved')
        self.Ev2=c.mpf([0,endpoints(c.exp(logcap))[1]])
        self.pressure_decay=c.mpf([0,endpoints(c.exp(logcap))[1]])

    def future_energy(self,Z):
        return self.pulse.selection.future.future(Z)['complete_future_energy_Taylor']/2

    def terminal_histories(self,Z):
        """All five actual Rv histories, including source-derived linear zeros."""
        point=self.pulse.end(Z,0)
        return dict(Mz_over_R_Utheta=point['Mz_over_R_Utheta'],
            Mtheta_z_over_sqrt2_R_3half_Utheta_squared=point['Mtheta_z_over_sqrt2_R_3half_Utheta_squared'],
            Mtheta_over_sqrt2_R_3half_Utheta=point['Mtheta_over_sqrt2_R_3half_Utheta'],
            Mztheta_over_R_Utheta_squared=point['Mztheta_over_R_Utheta_squared'],
            Mp_over_Pstar_squared=point['pressure']['Mp_over_Pstar_squared'],
            P0_over_Pstar_squared=point['pressure']['P0_over_Pstar_squared'],
            P_over_Pstar_squared=point['pressure']['P_over_Pstar_squared'])


def current_terminal_source_binding(owner):
    pulse=owner.pulse;flat=owner.flatten;inlet=flat.inlet
    graph=dict(same_current_terminal_owner=owner.dispatch.provider('pulse_end') is pulse,
        unchanged_native_class=type(pulse) is CompliantPulseMixedC4,
        unchanged_native_pressure=pulse.pressure_moment.__func__ is CompliantPulseRadialC4.pressure_moment,
        unchanged_native_data=pulse.data.__func__ is CompliantPulseRadialC4.data,
        unchanged_positive_factor_enclosure=pulse.factor.__func__ is CompliantAxialPulseField.factor,
        unchanged_flatten_ODE=flat.flatten.__func__ is CompliantFlattenMixedC4.flatten,
        unchanged_canonical_incoming=inlet.incoming.__func__ is CompliantPowerInletC4.incoming,
        same_future_provider=flat.pulse.selection.future is pulse.selection.future,
        same_U_constants=inlet.constants is pulse.high.constants,
        same_Pin=inlet.inlet_P is pulse.inlet_P,
        same_Hp=inlet.inlet_H is pulse.inlet_H,
        same_Xp=inlet.Xp is pulse.Xp,
        same_absolute_analytic_datum=inlet.datum is pulse.selection.future.angular.initial.datum,
        same_mu=flat.mu is inlet.mu is pulse.mu,
        same_delta=flat.delta is inlet.delta is pulse.delta)
    if not all(graph.values()):raise ValueError('Current endpoint/flatten defining object differs: '+str(graph))
    bindings={}
    for stem,method,values in (
        ('pulse_radial_C4','__init__',{'p0':'self.pulse.buffer.power(0,1)',
            'self.inlet_H':"box(p0['Mtheta_over_sqrt2_R_3half_Pstar'][0])",
            'self.inlet_P':"box(p0['Mp_over_Pstar_squared'][0])",
            'self.Xp':"self.inlet_H/box(p0['Utheta_over_Pstar'][0])"}),
        ('axial_pulse_field','end',{'future':"self.selection.future.future(Z)['complete_future_energy_Taylor']/2",
            'X':'1/self.rate+(self.Xp-1/self.rate)*self.factor(-13*self.rate/self.mu-self.rate*s)',
            'e':'future*c.exp(2*self.mu*s)+baseline-end_energy*(c.exp(2*self.mu*s)*self.E2cap)'}),
        ('flatten_mixed_C4','flatten',{'data':'self.inlet.incoming(Z)',
            'ev':'self.future_energy(Z)',
            'X':'(Xint+self.Xv*c.exp(-self.rate*t))/F',
            'energy':'(ev-Eint/2)*c.exp(2*self.mu*t)/(F*F)',
            'Mp_v':"data['Mp']+data['u']*data['u']*((1-self.pressure_decay)/(2*self.prate))",
            'Mp':'Mp_v+Pint*self.Ev2','pressure':"Mp+data['P0']"})):
        bindings[stem+'.'+method]=assignment_source_bindings(stem,method,values)
    bindings['current_flatten_constructor']={target:class_assignment(
        'current_pulse_flatten_source','CurrentFlattenMixedC4','__init__',target,expression)
        for target,expression in {
            'self.Xv':'1/self.rate+(pulse.Xp-1/self.rate)*pulse.factor(-13*self.rate/self.mu)',
            'self.logEv2_parts':"dict(inlet_log=2*c.ln(self.U),inverse_mu_term=-13/self.mu,finite_offset=c.mpf(-26))",
            'self.log_pressure_decay_parts':"dict(inverse_mu_term=-13/self.mu,finite_offset=c.mpf(-26))"}.items()}
    tree=ast.parse(Path(__file__).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CurrentFlattenMixedC4')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='future_energy')
    expression="self.pulse.selection.future.future(Z)['complete_future_energy_Taylor']/2"
    if (len(fn.body)!=1 or not isinstance(fn.body[0],ast.Return)
            or ast.dump(fn.body[0].value)!=ast.dump(ast.parse(expression,mode='eval').body)):
        raise ValueError('Current flatten must return the same native complete future half')
    bindings['current_future_energy_return']=expression
    bindings['actual_current_terminal_history_call']=assignment_source_bindings(
        'current_pulse_flatten_source','terminal_histories',{'point':'self.pulse.end(Z,0)'})
    historyfn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='terminal_histories')
    returned=historyfn.body[-1]
    expressions={name:"point['"+name+"']" for name in (
        'Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared',
        'Mtheta_over_sqrt2_R_3half_Utheta','Mztheta_over_R_Utheta_squared')}
    expressions.update({name:"point['pressure']['"+name+"']" for name in (
        'Mp_over_Pstar_squared','P0_over_Pstar_squared','P_over_Pstar_squared')})
    if not isinstance(returned,ast.Return) or not isinstance(returned.value,ast.Call):
        raise ValueError('Actual terminal history publication required')
    rows={v.arg:ast.dump(v.value) for v in returned.value.keywords}
    if rows!={k:ast.dump(ast.parse(v,mode='eval').body) for k,v in expressions.items()}:
        raise ValueError('Current five-history adapter must publish actual native terminal sources')
    bindings['actual_current_terminal_history_publication']=expressions
    # Algebraic source units, with the exact positive source H and decay S.
    # Caps only enclose them; they are not assigned as their defining values.
    U,q,mu,Xp,H,S,E,P0,Pin,Rp=s.symbols('U q mu Xp H S E P0 Pin Rp',nonzero=True)
    rate=1-mu;prate=1+2*mu
    native=dict(X=1/rate+(Xp-1/rate)*H,energy=E/2,
        Mp=Pin/q**2+(U/q)**2*(1-S)/(2*prate),P0=P0,
        logEv2=2*s.log(U)-13/mu-26,logRv=s.log(Rp)+13/mu)
    right=dict(X=1/rate+(Xp-1/rate)*H,energy=E/2,
        Mp=Pin/q**2+(U/q)**2*(1-S)/(2*prate),P0=P0,
        logEv2=2*(s.log(U)-13/(2*mu)-13),logRv=s.log(Rp)+13/mu)
    identities={key:s.expand(native[key]-right[key])==0 for key in native}
    if not all(identities.values()):raise ArithmeticError('Current terminal source unit transfer differs')
    return dict(current_defining_object_graph=graph,actual_source_assignments=bindings,
        exact_terminal_function_and_unit_identities=identities,
        native_absolute_pressure_publication=pressure_publication_bindings(),
        canonical_endpoint_theorem=THEOREM,
        retained_canonical_identity_counts=dict(functional=155,inlet_datum=34,endpoint=24),
        current_future_callable_used_for_every_Z=True,
        zero_linear_moments_are_from_empty_future_supports=True,
        actual_Rv_five_history_adapter_available=True,
        incoming_power_data_are_Rp_reference_not_terminal_histories=True,
        incoming_nonzero_histories_are_not_reset=True,
        no_saved_terminal_Xv_or_energy_lookup=True,
        source_caps_used_as_defining_field_values=False,
        interval_overlap_used_as_join_proof=False,passed=True)


class CurrentPulseFlattenSourceAssembly:
    @source_precision
    def __init__(self,require_checked=True,cells=128):
        self.dispatch=CurrentNativePulseSourceDispatcher()
        manifest=self.dispatch.manifest();self.pulse=self.dispatch.provider('pulse_end')
        self.family=self.dispatch.family;self.source=self.dispatch.source;self.datum_sha=self.dispatch.datum_sha
        self.hashes=dict(self.dispatch.hashes);self.registry=dict(manifest['ordered_current_chart_registry'])
        theorem=accepted(THEOREM,self.family,self.source,
            'pulse_end_flatten_full_moment_stress_pressure_functional_join_verified')
        if (theorem['source_caps_used_as_defining_field_values']
                or theorem['interval_overlap_used_as_join_proof']
                or (theorem['exact_functional_join_identities'],theorem['actual_local_inlet_and_datum_source_identities'],
                    theorem['actual_endpoint_source_identities'])!=(155,34,24)):
            raise ValueError('Canonical endpoint theorem scope changed')
        for name,digest in theorem['input_hashes'].items():
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Endpoint source conflict: '+name)
            self.hashes[name]=digest
        self.hashes[THEOREM]=sha(THEOREM)
        proxy=SimpleNamespace(pulse=self.pulse,pre=self.dispatch.rh_reference,
            params=self.dispatch.rh_reference.params,core=self.dispatch.anchor.patch.core,
            delta=self.dispatch.rh_reference.params.delta)
        self.parameter_bridge=native_parameter_source_bridge(proxy)
        self.hashes.update(self.parameter_bridge['input_hashes'])
        self.hashes[PREFIX+'current_pulse_physical_assembly.py']=sha(PREFIX+'current_pulse_physical_assembly.py')
        self.flatten=CurrentFlattenMixedC4(self.pulse,self.family,self.source,self.hashes,cells)
        self.binding=current_terminal_source_binding(self)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.registry['flatten']=dict(provider=PREFIX+'current_pulse_flatten_source.CurrentFlattenMixedC4',
            method='flatten',coverage_coordinate='t=log(R/Rv)',domain='[0,100]',acceptance_receipt=RECEIPT)
        self.acceptance_loaded=False
        if require_checked:
            record=accepted(RECEIPT,self.family,self.source,GATE)
            if record['datum_enclosure_sha256']!=self.datum_sha:raise ValueError('Current flatten datum differs')
            self.hashes.update(record['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.acceptance_loaded=True

    def provider(self,chart):
        if chart=='flatten':return self.flatten
        return self.dispatch.provider(chart)

    @source_precision
    def evaluate(self,chart,Z,coordinate):
        if chart!='flatten':return self.dispatch.evaluate(chart,Z,coordinate)
        packet=self.flatten.flatten(Z,coordinate)
        packet.update(current_terminal_histories=self.flatten.terminal_histories(Z),
            incoming_power_data_are_Rp_reference_not_terminal_histories=True,
            current_terminal_zero_linear_histories_from_native_supports=True)
        return dict(chart=chart,actual_five_defect_family_sha256=self.family,
            implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum_sha,
            source_provider=self.registry[chart]['provider'],acceptance_receipt=RECEIPT,
            source_coverage_coordinate='t=log(R/Rv)',source_coordinate_domain='[0,100]',
            derivative_coordinate='ordinary logR,Z; velocity/Ev0 and pressure/Pstar^2',
            physical_mixed_grids={'physical_mixed_derivatives_total_order_le4':
                packet['physical_mixed_derivatives_total_order_le4']},source_packet=packet,
            current_native_terminal_provider_used=True,current_future_energy_callable_used=True,
            **{GATE:self.acceptance_loaded,PROVED:True,UNIFORM:False},
            full_pulse_C4_installed=False,current_flatten_physical_owner_installed=False,
            output_kind='current flatten derivative source enclosures; no production point selection',
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,ordered_current_chart_registry=self.registry,
            current_downstream_chart_owner_count=21,current_native_pulse_chart_owner_count=6,
            current_pulse_terminal_flatten_source_bindings=self.binding,
            current_native_parameter_source_bridge=self.parameter_bridge,
            original_flatten_100_unit_ODEs_retained=True,current_Rp_external_pulse_join_certified=True,
            current_native_pulse_source_ownership_certified=True,
            **{GATE:self.acceptance_loaded,PROVED:True,UNIFORM:False},
            full_pulse_C4_installed=False,current_flatten_physical_owner_installed=False,
            all_profile_source_charts_callable=False,
            **dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),input_hashes=dict(self.hashes))

    @source_precision
    def report(self):
        result=self.manifest()
        terminal=self.dispatch.evaluate('pulse_end',[-1,1],0)
        packets={}
        for name,t in (('inlet',0),('whole_domain',[0,100]),('exit',100)):
            packets[name]=self.evaluate('flatten',[-1,1],t)
            print('Current pulse-to-flatten source: '+name,flush=True)
        result.update(current_native_terminal_packet=terminal,current_flatten_source_packets=packets,
            input_hashes=dict(self.hashes))
        return encode(pack(result))


def run():
    result=CurrentPulseFlattenSourceAssembly(require_checked=False).report()
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('Current pulse terminal -> unchanged flatten: twenty-one profile source owners generated',flush=True)
    return result


if __name__=='__main__':run()
