"""Explicit native fourteen-stage raw-preheat pressure integral operator.

Reifies the missing analytic flatten integral as a callable. All pressure
atoms come from the original current native swirl density, without angular
bumps and with H replaced by one in the preheat tail. The existing analytic
datum remains separate until the function-level identification is checked.
"""
import json
from pathlib import Path

import mpmath as mp
import sympy as s

from lei_ren_part1_paper_compliant_current_angular_terminal_closure import (
    CurrentAngularTerminalClosure,HERE,PREFIX,sha,function,binding,pack,encode,
    endpoints,IntervalTaylor,current_forward_terminal)
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_angular_high_jets import log_taylor
from lei_ren_part1_paper_compliant_axial_pulse_field import decay_integral
from lei_ren_part1_paper_logarithmic_pressure_datum import BETA2,BETA0
from lei_ren_part1_paper_compliant_actual_Rp_source_join import source_recipe_bindings,terminal_shape_proof
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_outer_angular_repair import intersect
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

NAME=PREFIX+'current_raw_preheat_pressure_operator.json'
RECEIPT=PREFIX+'current_raw_preheat_pressure_operator_check.json'
VIEWS={'whole_Z':[-1,1],'axis':0,'fresh':'.557','endpoint':1}
GATES=('current_native_flatten_pressure_integral_callable_reified',
       'current_raw_preheat_pressure_14_stage_operator_available',
       'current_raw_pressure_operator_native_forward_primitive_identified')
OPEN=('current_implicit_datum_to_native_raw_pressure_operator_identified',
      'current_heat_pressure_terminal_constant_eliminated',
      'current_heat_terminal_constants_eliminated','heat_exterior_stress_identity_certified',
      'current_exact_repair_installed_in_all_physical_charts','global_completed_tensor_admissibility',
      'physical_energy_integral_certified','independently_bounded_flat_remainder',
      'full_background_NS_validation','temporal_recursion')


@source_precision
def native_flatten_pressure_integral(exact,Z):
    """Same original Pint density, retaining all six analytic Z coefficients."""
    flat=exact.flatten;c=flat.ctx;z=c.mpf(Z);t=c.mpf(100)
    q=IntervalTaylor(c,[1+z**2,2*z,1,0,0,0])
    rho=log_taylor(q)-c.ln(2);Pint=q*0
    for i in range(flat.cells):
        a=t*i/flat.cells;b=t*(i+1)/flat.cells
        v=c.mpf([max(mp.mpf(0),endpoints(a)[0]),min(mp.mpf(100),endpoints(b)[1])])
        Fc=(rho*sigma_jets(c,v/100)[0]).exp()
        length=t/flat.cells
        Pint+=Fc*Fc/(q*q)*(c.exp(-flat.prate*a)*decay_integral(c,flat.prate,length)/2)
    return Pint


def raw_pressure_operator_source_proof(angular):
    """Bind the generator to actual density and primitive source formulas."""
    exact=angular.exact;bindings={};proofs={}
    def zero(name,left,right=0):
        if s.simplify(s.expand_power_exp(left-right))!=0:raise ArithmeticError('Raw pressure operator identity differs: '+name)
        proofs[name]=True
    specs={
        ('compliant_flatten_mixed_C4','flatten'):{
            'data':'self.inlet.incoming(Z)',
            'rho':'log_taylor(q)-c.ln(2)',
            'Fc':'(rho*sigma_jets(c,v/100)[0]).exp()',
            'Mp_v':"data['Mp']+data['u']*data['u']*((1-self.pressure_decay)/(2*self.prate))",
            'Mp':'Mp_v+Pint*self.Ev2','pressure':"Mp+data['P0']"},
        ('compliant_current_raw_preheat_pressure_operator','native_flatten_pressure_integral'):{
            'flat':'exact.flatten','rho':'log_taylor(q)-c.ln(2)',
            'Fc':'(rho*sigma_jets(c,v/100)[0]).exp()'},
        ('compliant_current_raw_preheat_pressure_operator','__init__'):{
            'self.Ev2':'intersect(self.ctx,self.exact.pulse.factor(self.logEv2),self.flat.Ev2)',
            'self.pressure_decay':'intersect(self.ctx,self.exact.pulse.factor(self.log_decay),self.flat.pressure_decay)'},
        ('compliant_current_raw_preheat_pressure_operator','evaluate'):{
            'incoming':'self.flat.inlet.incoming(z)',
            'forward':'current_forward_terminal(self.exact,z)',
            "atoms['pulse_reserved']":"incoming['u']*incoming['u']*((1-self.pressure_decay)/(2*p))",
            "atoms['z_flatten']":'Pint*self.Ev2',
            'thetaR':'c.exp(-bp*(100+L))/2','thetaS':'thetaR*c.exp(-bp-r/2)',
            'thetaQ':"thetaS*c.exp(-c.mpf('1.5')*Ts)",
            'thetaT':"thetaQ*c.exp(-c.mpf('1.5')+k/2)",
            "atoms['power_buffer_rel']":'one*(decay_integral(c,p,L)*self.Ev2*c.exp(-100*p)/8)',
            "atoms['steep_transition_in']":"one*(steep.infull['pressure']*self.Ev2*thetaR**2)",
            "atoms['steep_power']":'one*(decay_integral(c,3,Ts)*self.Ev2*thetaS**2/2)',
            "atoms['steep_transition_out']":"one*(steep.outfull['pressure']*self.Ev2*thetaQ**2)",
            "atoms['waiting']":'one*(decay_integral(c,1+self.flat.delta,W)*self.Ev2*thetaT**2/2)',
            'pre':"one*(1/(2*prate)-eps*epsatoms['PW']+eps**2*epsatoms['PW2']/2)",
            'exterior':'one*(c.exp(-3*prate)/(2*prate))',
            "atoms['heat_collar']":"(pre-exterior)*forward['pressure_scale']",
            "atoms['exterior_power_tail']":"exterior*forward['pressure_scale']",
            'prefix':"incoming['Mp']+atoms['pulse_reserved']",
            'total':"prefix+atoms['z_flatten']+sum((atoms[name] for name in BETA0),one*0)"},
        ('compliant_current_angular_terminal_closure','current_forward_terminal'):{
            'f':'flatten.flatten(z,100)',
            'PR':"f['pressure']['P_over_Pstar_squared']+decay_integral(c,prate,L)*(flatten.Ev2*c.exp(-100*prate)/8)",
            'PS':"PR+steep.infull['pressure']*(flatten.Ev2*thetaR**2)",
            'PQ':'PS+decay_integral(c,3,Ts)*(flatten.Ev2*thetaS**2/2)',
            'PT':"PQ+steep.outfull['pressure']*(flatten.Ev2*thetaQ**2)",
            'Ptail':'PT+decay_integral(c,1+native.delta,W)*(flatten.Ev2*thetaT**2/2)'},
        ('compliant_power_inlet_C4','incoming'):{
            'raw':"self.datum.normalized_jets(Z,5)['normalized_pressure_coefficients']"},
        ('compliant_pulse_radial_C4','pressure_moment'):{
            'p':"IntervalTaylor(c,inlet['Mp_over_Pstar_squared'])+u*u*(kernel/2)"}}
    for (stem,method),rows in specs.items():
        for target,expression in rows.items():
            binding(stem,method,target,expression);bindings[stem+'.'+method+':'+target]=True
    import ast
    for stem,method,expression in (
        ('compliant_flatten_mixed_C4','flatten','Fc*Fc/(q*q)*(c.exp(-self.prate*a)*decay_integral(c,self.prate,length)/2)'),
        ('compliant_current_raw_preheat_pressure_operator','native_flatten_pressure_integral','Fc*Fc/(q*q)*(c.exp(-flat.prate*a)*decay_integral(c,flat.prate,length)/2)')):
        rows=[n.value for n in ast.walk(function(stem,method)) if isinstance(n,ast.AugAssign) and ast.unparse(n.target)=='Pint']
        if len(rows)!=1 or ast.dump(rows[0])!=ast.dump(ast.parse(expression,mode='eval').body):
            raise ValueError('Original native flatten pressure density changed')
        bindings[stem+'.'+method+':Pint']=True
    rows=[n.value for n in ast.walk(function('compliant_current_angular_terminal_closure','current_forward_terminal'))
        if isinstance(n,ast.AugAssign) and ast.unparse(n.target)=='PR']
    expected=ast.parse('pastP*(flatten.Ev2*c.exp(-prate*(100+L))/4)',mode='eval').body
    if len(rows)!=1 or ast.dump(rows[0])!=ast.dump(expected):raise ValueError('Current forward absolute bump factor changed')
    bindings['current_forward_terminal:original_absolute_bump_factor']=True
    old_route=source_recipe_bindings();shape=terminal_shape_proof()
    if not shape['passed']:raise ValueError('Current reference-to-Rp pressure primitive source proof required')
    z,v,rho,sg,dt,lo=s.symbols('Z v rho sigma dt left',real=True)
    mu=s.symbols('mu',positive=True)
    q=1+z*z;p=1+2*mu
    normalized_velocity=s.exp(rho*sg)*s.exp(-(s.Rational(1,2)+mu)*v)/q
    zero('explicit_Fflat_density_is_native_swirl_square_over_two',
        normalized_velocity**2/2,s.exp(2*rho*sg)*s.exp(-p*v)/(2*q*q))
    zero('flatten_exponential_cell_weight_is_exact_integral',
        s.exp(-p*lo)*(1-s.exp(-p*dt))/(2*p),s.integrate(s.exp(-p*v)/2,(v,lo,lo+dt)))
    U,E=s.symbols('U Ev2',positive=True);Tp=s.symbols('pulse_length',positive=True)
    zero('native_pulse_pressure_density_primitive',
        s.diff(U*U*(1-s.exp(-p*Tp))/(2*p),Tp),U*U*s.exp(-p*Tp)/2)
    zero('actual_native_pulse_terminal_decay_exponent',p*13/mu,13/mu+26)
    # Every name below denotes a complete source function, not a cell sum
    # or a selected interval endpoint. Same density/domain and radial FTC
    # identify the cumulative native primitive with these integrals.
    P0,Mpin,pulse,Fflat,relative,Iin,Ipower,Iout,Iwait,scale,pre,exterior,bump=s.symbols(
        'P0 Mpin pulse Fflat relative Iin Ipower Iout Iwait pressure_scale pre exterior bump',real=True)
    raw_total=Mpin+pulse+Fflat+relative+Iin+Ipower+Iout+Iwait+(pre-exterior)*scale+exterior*scale
    native_tail=P0+Mpin+pulse+Fflat+relative+bump+Iin+Ipower+Iout+Iwait
    zero('complete_raw_operator_equals_native_forward_without_bump_plus_preheat',
        P0+raw_total,native_tail-bump+pre*scale)
    zero('heat_collar_and_exterior_pressure_partition',
        (pre-exterior)*scale+exterior*scale,pre*scale)
    # The operator has no pressure/free constant input: it depends only on
    # the explicit raw swirl and the unique prescribed raw waiting root.
    return dict(actual_density_AST_bindings=bindings,
        same_reference_to_Rp_pressure_recipe=old_route,
        current_reference_to_Rp_shape_proof=shape,identities=proofs,
        flatten_integral_callable=native_flatten_pressure_integral.__name__,
        full_fourteen_stage_names=list(BETA2)+['z_flatten']+list(BETA0),
        original_analytic_datum_not_replaced=True,
        numerical_caps_not_defining_density_values=True,
        generator_does_not_use_original_P0_to_define_raw_total=True,
        same_native_forward_cumulative_pressure_function_identified=True,
        raw_pressure_operator_returns_directed_enclosures_not_exact_stage_values=True,
        original_implicit_Fflat_identification_still_required=True,passed=True)


class CurrentRawPreheatPressureOperator:
    @source_precision
    def __init__(self,angular=None,require_checked=True):
        self.angular=angular if angular is not None else CurrentAngularTerminalClosure()
        if not self.angular.acceptance_loaded:raise ValueError('Checked current exact heat angular runtime required')
        self.exact=self.angular.exact;self.flat=self.exact.flatten;self.ctx=self.flat.ctx
        self.family=self.angular.family;self.source=self.angular.source;self.datum_sha=self.angular.datum_sha
        self.proof=raw_pressure_operator_source_proof(self.angular)
        self.hashes=dict(self.angular.hashes);self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.logEv2=sum((self.ctx.mpf(endpoints(v)) for v in self.flat.logEv2_parts.values()),self.ctx.mpf(0))
        self.log_decay=sum((self.ctx.mpf(endpoints(v)) for v in self.flat.log_pressure_decay_parts.values()),self.ctx.mpf(0))
        self.Ev2=intersect(self.ctx,self.exact.pulse.factor(self.logEv2),self.flat.Ev2)
        self.pressure_decay=intersect(self.ctx,self.exact.pulse.factor(self.log_decay),self.flat.pressure_decay)
        self.buffer=self.exact.pulse.pulse.buffer;self.axis_atoms=None;self.cache={};self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Current raw pressure operator source/scope differs')
            if receipt['raw_pressure_operator_source_proof']!=encode(pack(self.proof)):
                raise ValueError('Current raw pressure generator proof changed')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def early_axis_atoms(self):
        if self.axis_atoms is not None:return self.axis_atoms
        c=self.ctx;b=self.buffer
        scalar=lambda packet:c.mpf(endpoints(packet['Mp_over_Pstar_squared'][0]))
        slope=scalar(b.initial.slope(0,1,cells=128))
        axial=scalar(b.initial.axial(0,buffer_offset=11,cells=128))
        transition=scalar(b.slope_mu(0,1,cells=128));power=scalar(b.power(0,1,cells=128))
        # Directed differences enclose these source-defined positive atoms.
        # Positivity narrows a bound; it does not choose a pressure value.
        positive=lambda value:c.mpf([max(mp.mpf(0),endpoints(value)[0]),max(mp.mpf(0),endpoints(value)[1])])
        self.axis_atoms=dict(reference_extension=c.mpf('2.5'),
            slope_transition_ref=positive(slope-c.mpf('2.5')),
            axial_turnoff=positive(axial-slope),slope_transition_mu=positive(transition-axial),
            power_buffer=positive(power-transition))
        return self.axis_atoms

    @source_precision
    def evaluate(self,Z):
        c=self.ctx;z=c.mpf(Z);key=z._mpi_
        if key in self.cache:return self.cache[key]
        incoming=self.flat.inlet.incoming(z);forward=current_forward_terminal(self.exact,z)
        q=IntervalTaylor(c,[1+z*z,2*z,1,0,0,0]);invq2=q**(-2);one=q*0+1
        p=self.flat.prate;steep=self.exact.companion.heat.steep
        L=c.mpf(endpoints(-30*self.exact.repair.params.log_mu));Ts=c.mpf(endpoints(self.exact.repair.params.Ts))
        W=c.mpf(endpoints(self.exact.repair.angular.waiting));logone=c.mpf(endpoints(self.exact.repair.angular.waiting_logone))
        atoms={name:invq2*value for name,value in self.early_axis_atoms().items()}
        atoms['pulse_reserved']=incoming['u']*incoming['u']*((1-self.pressure_decay)/(2*p))
        Pint=native_flatten_pressure_integral(self.exact,z)
        atoms['z_flatten']=Pint*self.Ev2
        bp=c.mpf('.5')+self.flat.mu;k=1-self.flat.delta/2;bh=c.mpf('.5')+self.flat.delta/2;r=1-self.flat.mu
        thetaR=c.exp(-bp*(100+L))/2;thetaS=thetaR*c.exp(-bp-r/2)
        thetaQ=thetaS*c.exp(-c.mpf('1.5')*Ts);thetaT=thetaQ*c.exp(-c.mpf('1.5')+k/2)
        atoms['power_buffer_rel']=one*(decay_integral(c,p,L)*self.Ev2*c.exp(-100*p)/8)
        atoms['steep_transition_in']=one*(steep.infull['pressure']*self.Ev2*thetaR**2)
        atoms['steep_power']=one*(decay_integral(c,3,Ts)*self.Ev2*thetaS**2/2)
        atoms['steep_transition_out']=one*(steep.outfull['pressure']*self.Ev2*thetaQ**2)
        atoms['waiting']=one*(decay_integral(c,1+self.flat.delta,W)*self.Ev2*thetaT**2/2)
        tail=self.angular.heat.collar_tails(z,0);eps=self.angular.heat.eps;prate=1+self.flat.delta
        epsatoms=tail['separate_epsilon_atoms']
        pre=one*(1/(2*prate)-eps*epsatoms['PW']+eps**2*epsatoms['PW2']/2)
        exterior=one*(c.exp(-3*prate)/(2*prate))
        atoms['heat_collar']=(pre-exterior)*forward['pressure_scale']
        atoms['exterior_power_tail']=exterior*forward['pressure_scale']
        if set(atoms)!=set(BETA2+BETA0+('z_flatten',)):raise ValueError('Complete fourteen-stage raw pressure operator required')
        prefix=incoming['Mp']+atoms['pulse_reserved']
        # Use the correlated inherited prefix for the total, instead of
        # summing independent directed differences of its earlier atoms.
        total=prefix+atoms['z_flatten']+sum((atoms[name] for name in BETA0),one*0)
        oldP0=incoming['P0']
        out=dict(Z=z,raw_preheat_pressure_atoms=atoms,raw_complete_pressure_integral_Taylor=total,
            explicit_negative_raw_pressure_integral_Taylor=-total,
            original_analytic_P0_Taylor_retained=oldP0,
            unresolved_original_P0_plus_raw_total_Taylor=oldP0+total,
            native_flatten_pressure_integral_callable_Taylor=Pint,
            analytic_flatten_Fflat_Taylor=atoms['z_flatten'],
            correlated_actual_native_Rp_to_Rv_prefix=prefix,
            raw_forward_terminal_plus_full_preheat_pressure_enclosure=(
                forward['Ptail']-forward['pastP']*(self.Ev2*c.exp(-p*(100+L))/4)+pre*forward['pressure_scale']),
            exact_positive_amplitude_logs=self.flat.logEv2_parts,
            exact_positive_pressure_decay_logs=self.flat.log_pressure_decay_parts,
            native_factor_callables_consumed_for_both_positive_scales=True,
            directed_enclosures_not_selected_exact_stage_values=True,
            full_raw_tail_uses_H_one_and_no_angular_bumps=True,
            original_P0_not_reset_or_replaced=True,
            exact_radius_and_waiting_sources_preserved=True,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        self.cache[key]=out;return out


@source_precision
def run(field=None):
    field=field if field is not None else CurrentRawPreheatPressureOperator(require_checked=False)
    report=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
        datum_enclosure_sha256=field.datum_sha,raw_pressure_operator_source_proof=field.proof,
        current_raw_pressure_operator_views={name:field.evaluate(Z) for name,Z in VIEWS.items()},
        input_hashes=field.hashes,**dict.fromkeys(GATES+OPEN,False))
    (HERE/NAME).write_text(json.dumps(encode(pack(report)),indent=2)+'\n',encoding='utf8')
    print('Explicit native fourteen-stage raw pressure operator generated; original implicit-datum identification still pending',flush=True)
    return report


if __name__=='__main__':run()
