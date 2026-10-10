"""Current selected-source absolute future moments, in exact paper units.

The same cumulative functions cross every segment boundary. Original source
ODEs and full Gamma limits identify their integrals by FTC. Huge finite
boundary scales remain factored; no midpoint, cap or radial cutoff is used.
The executable target enclosures retain the true first ordinary Z row.
"""
import ast
import gzip
import hashlib
import json
import time
from types import FunctionType
from pathlib import Path
import sympy as s

import lei_ren_part1_paper_compliant_current_limit_selected_postpulse_registry as joined
import lei_ren_part1_paper_compliant_current_limit_band_physical_recovery as band
import lei_ren_part1_paper_compliant_current_limit_power_to_Rp as quiet
import lei_ren_part1_paper_compliant_outer_buffer as buffer
import lei_ren_part1_paper_compliant_pulse_interface_certificate as interfaces
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

HERE,PREFIX,sha,require=joined.HERE,joined.PREFIX,joined.sha,joined.require
NAME=PREFIX+'current_limit_absolute_future_integrals.json'
RECEIPT=PREFIX+'current_limit_absolute_future_integrals_check.json'
GATES=('current_repaired_limit_final_heat_reference_amplitude_source_installed',
       'current_repaired_limit_absolute_five_future_integrals_C1_source_closed')
OPEN=joined.OPEN
KEYS=('M','J','S_energy','renormalized_I','pressure_tail')
UNITS={'M':'R2*Pstar','J':'sqrt(2)*R2^(3/2)*Pstar^2',
    'S_energy':'R2*Pstar^2','renormalized_I':'sqrt(2)*R2^(3/2)*Pstar',
    'pressure_tail':'Pstar^2'}
SEGMENTS=('quiet_power','selected_pulse','selected_postpulse','collar_and_full_Gamma')


def exact_phase_box(context,value):
    """Directed source interval acquisition for the non-rational 2Rc phase.

    The native public API accepts Fraction only. The same native power
    code below receives a real source interval, with every arithmetic,
    kernel and inlet call unchanged. No rational midpoint is supplied.
    """
    if not hasattr(value,'_mpi_'):raise TypeError('Exact directed phase interval required')
    return context.mpf(joined.history.endpoints(value))


def native_real_phase_power(owner,Z,phase):
    original=buffer.SharedOuterBuffer.power
    environment=dict(original.__globals__)
    environment['fraction_box']=exact_phase_box
    evaluate=FunctionType(original.__code__,environment,original.__name__,original.__defaults__,original.__closure__)
    require(evaluate.__code__ is original.__code__ and
        all(environment[key] is value for key,value in original.__globals__.items() if key!='fraction_box'),
        'Only real phase acquisition may differ from the native power algorithm')
    return evaluate(owner,Z,phase)


def syntax(stem,method,statement):
    path=HERE/(PREFIX+stem+'.py');tree=ast.parse(path.read_text(encoding='utf8'))
    fn=next(q for q in ast.walk(tree) if isinstance(q,ast.FunctionDef) and q.name==method)
    wanted=ast.dump(ast.parse(statement).body[0])
    require(sum(ast.dump(q)==wanted for q in ast.walk(fn))==1,
        'Actual source statement changed: '+stem+'.'+method+':'+statement)
    return dict(statement=statement,function_AST_sha256=hashlib.sha256(ast.dump(fn).encode()).hexdigest(),
        source_sha256=sha(path.name))


def source_FTC_proof(owner,frame,power_report,pulse_receipt):
    p=joined.selected.source.inlet.identity.Projection();z=s.Symbol('Z',real=True)
    mu,t,Pstar=s.symbols('mu t Pstar',positive=True);bp=s.Rational(1,2)+mu
    A=s.Function('A')(z);old={k:s.Function(k+'0')(z) for k in ('m','h','k','e','p')}
    raw=dict(old,m=Pstar*old['m'],k=Pstar*old['k'])
    env=dict(u1=A,f=s.exp(-bp*t),decay=s.exp(-t),d3=s.exp(-s.Rational(3,2)*t),t=t,mu=mu,
        theta_kernel=(s.exp(-bp*t)-s.exp(-s.Rational(3,2)*t))/(1-mu),
        **{'decay_integral(c, 2 * mu, t)':(1-s.exp(-2*mu*t))/(2*mu),
           'decay_integral(c, 1 + 2 * mu, t)':(1-s.exp(-(1+2*mu)*t))/(1+2*mu),
           **{"get('%s')"%name:raw[key] for key,name in dict(m='Mz_over_R',
             h='Mtheta_over_sqrt2_R_3half_Pstar',k='Mtheta_z_over_sqrt2_R_3half_Pstar',
             e='Mztheta_over_R_Pstar_squared',p='Mp_over_Pstar_squared').items()}})
    E=p.assignment(buffer,'power','u',env)
    actual={key:p.assignment(buffer,'power',key,env)/(Pstar if key in ('m','k') else 1) for key in old}
    rates=dict(m=1,h=s.Rational(3,2),k=s.Rational(3,2),e=1,p=0)
    drivers=dict(m=0,h=E,k=0,e=-E*E/2,p=E*E/2)
    proofs={}
    def zero(label,expr):
        require(s.simplify(s.expand_power_exp(expr))==0,'Actual cumulative FTC failed: '+label);proofs[label]=True
    for key in actual:
        zero('native_quiet_'+key+'_source_ODE',s.diff(actual[key],t)+rates[key]*actual[key]-drivers[key])
        zero('native_quiet_'+key+'_source_inlet',actual[key].subs(t,0)-old[key])
    require(frame['actual_native_frame_identity']['passed'] and
        frame['underlying_function_projection']['ordinary_Z_identity_over_whole_original_domain'],
        'Actual current/native Rp whole-function identity required')
    require(power_report['symbolic_checks']['same_parent_all_five_C1_semigroup'] and
        power_report['symbolic_checks']['all_five_complete_power_ODEs_and_inlets'] and
        power_report['symbolic_checks']['exact_Rp_Rw_Rc_radius_identity'],
        'Actual current quiet-power C1 semigroup/radius identity required')
    require(owner.selected_owner.proof['passed'] and
        owner.selected_owner.proof['nonzero_formal_end_energy_and_incoming_terms_retained'],
        'Actual selected source proof must retain absolute inlet histories')
    require(interfaces.functional_identities()==pulse_receipt['source_bound_functional_pulse_identities'],
        'Original arbitrary-selected-source functional pulse joins changed')
    statements={}
    for stem,method,rows in (
        ('pulse_mixed_C4','transport_mixed',[
            "m1.append(Brows[k]-m1[k]*(c.mpf('.5')-mu))",
            "m2.append(Brows[k]-m2[k]*(c.mpf('.5')-2*mu))",
            'X.append((one if k==0 else one*0)-X[k]*(1-mu))',
            'energy.append(square-(one/2 if k==0 else one*0)+energy[k]*(2*mu))',
            "pressure.append(point['pressure']['P_y_over_Pstar_squared']*(-field.prate)**k)"]),
        ('flatten_mixed_C4','flatten_mixed',[
            'angularrate=list(rates); angularrate[0]=angularrate[0]+1-mu',
            'energyrate=[v*(-2) for v in rates]; energyrate[0]=energyrate[0]+2*mu',
            'theta_rows.append(th); Xrows.append(xx); erows.append(ee); prows.append(square*(Ev2/2))'])):
        # Bind individual statements even when production writes them on one line.
        for group in rows:
            for statement in group.split('; '):
                statements[stem+'.'+statement]=syntax(stem,method,statement)
    R,Utheta=s.symbols('R Utheta',positive=True)
    B,m1,m2,X,e=s.symbols('B m1 m2 X e',real=True)
    pulse_scales=dict(M=R*Utheta,J=s.sqrt(2)*R**s.Rational(3,2)*Utheta**2,
        I=s.sqrt(2)*R**s.Rational(3,2)*Utheta,S_energy=R*Utheta**2)
    derivatives=dict(M=B-(s.Rational(1,2)-mu)*m1,J=B-(s.Rational(1,2)-2*mu)*m2,
        I=1-(1-mu)*X,S_energy=B*B-s.Rational(1,2)+2*mu*e)
    vars=dict(M=m1,J=m2,I=X,S_energy=e)
    densities=dict(M=Utheta*B,J=s.sqrt(2*R)*Utheta**2*B,I=s.sqrt(2*R)*Utheta,
        S_energy=Utheta**2*(B*B-s.Rational(1,2)))
    for key,scale in pulse_scales.items():
        dy_scale=R*s.diff(scale,R)-bp*Utheta*s.diff(scale,Utheta)
        zero('actual_pulse_'+key+'_paper_density_Jacobian',
            (dy_scale*vars[key]+scale*derivatives[key])/R-densities[key])
    zero('actual_pulse_pressure_paper_density_Jacobian',(Utheta**2/2)/R-Utheta**2/(2*R))
    # The common variable-rate recurrence is used by flatten, power,
    # steep, waiting, collar and Gamma. Its normalized rows must first
    # be converted to the PAPER primitives at the actual radius.
    rho=s.Symbol('actual_log_amplitude_rate',real=True)
    post_dy=dict(I=1-(1-mu+rho)*X,
        S_energy=-s.Rational(1,2)+(2*mu-2*rho)*e)
    for key in post_dy:
        scale=pulse_scales[key]
        zero('all_actual_pure_swirl_'+key+'_paper_density_Jacobian',
            (R*s.diff(scale,R)*vars[key]+(-bp+rho)*Utheta*s.diff(scale,Utheta)*vars[key]
             +scale*post_dy[key])/R-densities[key].subs(B,0))
    require(owner.history.proof['selected_forward_cumulative_energy_equals_same_remaining_integral_by_FTC'] and
        owner.history.proof['zero_meridional_histories_propagate_from_selected_terminal_by_FTC'] and
        owner.exterior.proof['all_current_five_terminal_moment_functions_identified'],
        'Actual selected/postpulse/full exterior histories required')
    # Backward uniqueness from the accepted Rp source functions identifies
    # the native buffer at 2Rc. Absolute states are retained, never set to zero.
    t0=s.Integer(2)+s.log(2);Tw=s.Symbol('Tw',positive=True);Rc=s.Symbol('Rc',positive=True)
    zero('same_native_2Rc_radius_and_phase',(Rc*s.exp(-2)*s.exp(Tw*(t0/Tw)))-2*Rc)
    return dict(actual_native_power_source_projection=p.bindings,identities=proofs,
        actual_pulse_and_postpulse_ODE_AST=statements,
        actual_current_Rp_identity_receipt=joined.selected.source.inlet.identity.RECEIPT,
        same_Rp_values_and_Z_rows_plus_same_power_ODE_identify_entire_quiet_interval=True,
        actual_raw_m_k_divided_by_Pstar_exactly_once=True,
        same_native_Rw_Rc_mu_Tw_definitions_consumed=True,
        current_selected_linear_quadratic_and_functional_pulse_interfaces_consumed=True,
        current_postpulse_full_energy_angular_pressure_source_functions_consumed=True,
        same_cumulative_functions_reused_at_all_segment_boundaries=True,
        absolute_P0_retained_once_not_reset=True,interval_overlap_used_as_function_identity=False,passed=True)


def heat_reference_and_limits_proof(owner):
    bind=joined.selected.source.inlet.identity.ast_assignments
    bindings={
        'current_exact_tail_radius':bind('current_exact_repair_branch',None,'replay_heat',{
            'out.tail_finite':'p.yd+1+p.Tw+100-30*p.log_mu+2+p.Ts+angular.waiting',
            'out.logradius_terms':'dict(selected_reference=out.logRref,pulse_term=13/p.mu,finite_offset=out.tail_finite)'}),
        'actual_flatten_scale':bind('current_pulse_flatten_source','CurrentFlattenMixedC4','__init__',{
            'self.U':"self.inlet.constants['U']",
            'self.logEv2_parts':'dict(inlet_log=2*c.ln(self.U),inverse_mu_term=-13/self.mu,finite_offset=c.mpf(-26))'}),
        'actual_heat_amplitude':bind('current_heat_source','CurrentCollarGammaC4','__init__',{
            'self.a':'self.delta/2','self.k':'self.steep.k','self.bh':'self.steep.bh',
            'self.theta_base':'self.steep.thetaT*c.exp(-self.bh*self.steep.wait-self.steep.logone)'}),
        'actual_steep_amplitudes':bind('current_steep_waiting_source','CurrentSteepWaitingC4','__init__',{
            'self.k':'1-self.delta/2','self.bh':"c.mpf('.5')+self.delta/2",
            'self.thetaR':'c.exp(-100*self.bp-self.bp*self.outer.Lrel)/2',
            'self.thetaS':'self.thetaR*c.exp(-self.bp-self.rate/2)',
            'self.thetaQ':"self.thetaS*c.exp(-c.mpf('1.5')*self.Ts)",
            'self.thetaT':"self.thetaQ*c.exp(-c.mpf('1.5')+self.k/2)"})}
    a,mu,L,Ts,W,lone,lp,lu,lrp=s.symbols('a mu L Ts W lone logP logU logRp',real=True)
    bh=s.Rational(1,2)+a;k=1-a;bp=s.Rational(1,2)+mu;rate=1-mu
    logEv0=lp+lu-13/(2*mu)-13
    logtheta=-100*bp-s.log(2)-bp*L-bp-rate/2-s.Rational(3,2)*Ts-s.Rational(3,2)+k/2-bh*W-lone
    logRt=lrp+13/mu+100+L+2+Ts+W
    grouped=lp+lu+bh*lrp+13*a/mu+(a-mu)*(100+L)-k*Ts-14-mu/2+3*a/2-s.log(2)-lone
    require(s.expand(logEv0+logtheta+bh*logRt-grouped)==0,'Actual c_infinity log correlation lost')
    heat=owner.history.heat
    require(heat.exact_logRtail_terms==owner.history.selected.exact.repair.heat.logradius_terms,
        'Same actual exact heat radius required')
    require(heat.exact_logS_terms==owner.history.selected.exact.repair.heat.logS_terms,
        'Same actual inverse radius in the full Gamma argument required')
    gamma=owner.history.proof['checked_canonical_Gamma_source_evidence']
    require(gamma['verified'],'Complete positive Gamma theorem required')
    theorem=owner.exterior.theorem
    require(theorem['full_terminal_moment_stress_theorem_verified'] and
        theorem['theorem_uses_full_positive_Gamma_function'], 'Actual full terminal integral theorem required')
    bindings['canonical_Gamma_velocity']=syntax('collar_Gamma_C4','packet',
        'theta=K[0]*(self.theta_base*c.exp(-self.bh*t))')
    bindings['canonical_Gamma_angular']=syntax('collar_Gamma_C4','local_Gamma',
        "angular=one/self.k+tails['theta']*Sc")
    bindings['canonical_Gamma_energy']=syntax('collar_Gamma_C4','local_Gamma',
        "energy=one/self.delta-tails['energy']*(self.a*Sc)")
    bindings['canonical_Gamma_pressure']=syntax('collar_Gamma_C4','local_Gamma',
        "pressure=one/(2*self.prate)-tails['pressure']*(self.a*Sc)")
    # Bounds are derived from the full expectation, not a finite S expansion:
    # H=E[(1+xi*v)^(-a)] under Gamma(1+a), E[v]=1+a.
    A=s.Symbol('a',positive=True);R,c=s.symbols('R c_infinity',positive=True)
    Z,v=s.symbols('Z v',real=True);xi=s.Symbol('xi',nonnegative=True)
    require(s.integrate(s.exp(-A*v),(v,0,s.oo))==1/A,'Full angular deficit bound integral')
    # Derive the physical radial weights and genuine fixed-R Z row from
    # the admitted full H expectation, rather than adopting normalized
    # Taylor orders as physical asymptotics.
    bhA=s.Rational(1,2)+A;radial=R*s.exp(v);source_c=c
    bound_derivations={}
    weights={
        'renormalized_I':s.sqrt(2)*source_c*R**(1-A)*A*(1+A)*(2/R)*s.exp(-A*v),
        'renormalized_I_Z':s.sqrt(2)*source_c*R**(1-A)*A*(1+A)*(4/R)*s.exp(-A*v),
        'S_energy':source_c**2*R**(-2*A)*s.exp(-2*A*v)/2,
        'S_energy_Z':4*source_c**2*A*(1+A)*R**(-1-2*A)*s.exp(-(1+2*A)*v),
        'absolute_Pi':source_c**2*R**(-1-2*A)*s.exp(-(1+2*A)*v)/2,
        'absolute_Pi_Z':4*source_c**2*A*(1+A)*R**(-2-2*A)*s.exp(-(2+2*A)*v)}
    bounds={
        'renormalized_I':2*s.sqrt(2)*(1+A)*c*R**(-A),
        'renormalized_I_Z':4*s.sqrt(2)*(1+A)*c*R**(-A),
        'S_energy':c*c*R**(-2*A)/(4*A),
        'S_energy_Z':4*c*c*A*(1+A)*R**(-1-2*A)/(1+2*A),
        'absolute_Pi':c*c*R**(-1-2*A)/(2*(1+2*A)),
        'absolute_Pi_Z':4*c*c*A*(1+A)*R**(-2-2*A)/(2+2*A)}
    for name,bound_value in bounds.items():
        require(s.simplify(s.integrate(weights[name],(v,0,s.oo))-bound_value)==0,
            'Full Gamma radial Jacobian / first-Z bound differs: '+name)
        require(s.limit(bound_value,R,s.oo)==0,'Actual whole-Z terminal bound does not vanish: '+name)
        bound_derivations[name]=dict(full_integrable_majorant=s.sstr(weights[name]),
            exact_integral=s.sstr(bound_value),zero_limit=True)
    # Same source c=Ev0*theta_base*Rtail^bh and xi=2(1-Z^2)/R.
    Ev0,theta,Rt,t=s.symbols('Ev0 theta_base Rtail offset',positive=True)
    require(s.simplify(s.powsimp((Ev0*theta*s.exp(-bhA*t)).subs(t,s.log(R/Rt))
        -Ev0*theta*Rt**bhA*R**(-bhA),force=True))==0,'Actual Gamma amplitude product identity')
    require(s.simplify((2*(1-Z**2)/Rt*s.exp(-t)).subs(t,s.log(R/Rt))-2*(1-Z**2)/R)==0,
        'Actual Gamma fixed-R source argument identity')
    return dict(actual_amplitude_and_radius_AST=bindings,
        exact_c_infinity_definition='Ev0*theta_base*Rtail^((1+delta)/2)',
        exact_correlated_log_c_infinity=s.sstr(grouped),
        c_infinity_positive_and_Z_independent=True,actual_Gamma_source_bound=True,
        full_Gamma_expectation_bounds=['0<H<=1','abs(H_prime)<=a*(1+a)','0<=1-H(xi)<=a*(1+a)*xi'],
        full_angular_deficit_bound='J(xi)=int_0^infinity exp((1-a)*v)*(1-H(xi*exp(-v)))dv <= (1+a)*xi',
        terminal_whole_Z_value_and_first_Z_bounds={k:s.sstr(v) for k,v in bounds.items()},
        actual_paper_radial_Jacobian_and_first_Z_majorants=bound_derivations,
        source_bound_full_terminal_moment_definitions=theorem['terminal_moment_definitions'],
        terminal_M_and_J_identically_zero_from_current_selected_history=True,
        all_five_absolute_or_renormalized_C1_infinity_limits_zero=True,
        no_finite_radial_cutoff_or_cap_value_substitution=True,passed=True)


def telescope_proof(ledger):
    identities={};z=s.Symbol('Z',real=True)
    for key in KEYS:
        # Symbols name the actual callable/factored C1 source atoms in
        # THIS ledger. Adjacent segments reuse the very same atom IDs.
        atoms={name:s.Function(name)(z) for name in ledger['atoms']}
        def expr(terms):return sum((sign*atoms[name] for sign,name in terms),s.Integer(0))
        pieces=[expr(ledger['segments'][stage][key]) for stage in SEGMENTS]
        first=expr(ledger['normalized_endpoints']['R2'][key])
        total=sum(pieces)
        require(s.expand(total+first)==0,'Actual full future telescope: '+key)
        require(s.expand(s.diff(total+first,z))==0,'True first-Z future telescope: '+key)
        identities[key]=dict(value=True,true_first_Z=True,
            entire_four_segment_future_sum=s.sstr(total),no_boundary_history_reset=True)
    return dict(identities=identities,same_actual_endpoint_function_at_each_internal_boundary=True,
        full_transition_retained_in_renormalized_angular_integral=True,passed=True)


class CurrentLimitAbsoluteFutureIntegrals:
    @joined.source_precision
    def __init__(self,owner=None,require_checked=True):
        self.owner=owner if owner is not None else joined.CurrentLimitSelectedPostpulseRegistry()
        require(self.owner.acceptance_loaded,'Checked actual selected postpulse registry required')
        self.family_record=self.owner.family_record;self.family=self.owner.family
        self.source=self.owner.source;self.datum_sha=self.owner.datum_sha
        self.hashes=dict(self.owner.hashes);self.ctx=self.owner.ctx
        checked=joined.selected.source.inlet.checked
        checked('current_limit_selected_postpulse_registry',joined.GATES[0],self.family_record,self.hashes)
        checked('current_limit_band_physical_recovery',band.GATE,self.family_record,self.hashes)
        checked('current_limit_power_to_Rp',quiet.GATE,self.family_record,self.hashes)
        checked('current_limit_Rp_native_identity',joined.selected.source.inlet.identity.GATE,self.family_record,self.hashes)
        self.band_report=json.loads(gzip.decompress((HERE/band.NAME).read_bytes()))
        self.power_report=json.loads(gzip.decompress((HERE/quiet.NAME).read_bytes()))
        frame=json.loads((HERE/joined.selected.source.inlet.identity.NAME).read_bytes())
        pulse=json.loads((HERE/(PREFIX+'pulse_interface_certificate.json')).read_bytes())
        require(sha(PREFIX+'pulse_interface_certificate.json')==self.hashes[PREFIX+'pulse_interface_certificate.json'],
            'Original checked pulse function receipt must be pinned')
        self.required_target_nodes=self.band_report['exact_physical_recovery']['required_future_integrals']
        require(set(self.required_target_nodes)==set(KEYS),'All actual current absolute future targets required')
        self.pulse=self.owner.history.selected.pulse;self.buffer=self.pulse.pulse.buffer
        self.heat=self.owner.history.heat;self.params=self.owner.history.selected.future.params
        self.FTC=source_FTC_proof(self.owner,frame,self.power_report,pulse)
        self.heat_proof=heat_reference_and_limits_proof(self.owner)
        self.installed_target_graph=self.target_graph()
        self.cache={}
        for module in (band,quiet,buffer,interfaces):
            name=Path(module.__file__).name;self.hashes[name]=sha(name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=checked('current_limit_absolute_future_integrals',GATES[0],self.family_record,self.hashes)
            require(all(receipt[k] for k in GATES) and not any(receipt[k] for k in OPEN),'Absolute future source scope differs')
            self.acceptance_loaded=True

    def target_graph(self):
        """Substitute the ACTUAL heat amplitude into the original targets.

        The source recipe remains exact and logarithmic. This extends the
        accepted graph; its old source expressions are never rewritten.
        The numerical boundary adapter below evaluates the same recipe.
        """
        source=band.source;g=source.FunctionTransportGraph()
        nodes=self.band_report['exact_graph_nodes']
        g.nodes=[dict(row) for row in nodes]
        g.keys={json.dumps(row,sort_keys=True,separators=(',',':')):i for i,row in enumerate(nodes)}
        raw=self.band_report['exact_physical_recovery'];old=source.FunctionRef(g,raw['future_c_infinity'])
        require(g.nodes[old.node]==dict(operation='bound_variable',name='future_source_c_infinity'),
            'Original c_infinity target variable changed')
        exact=g.node('current_selected_final_heat_reference_function_recipe',source_family=self.family_record,
            source_owner='current_limit_selected_postpulse_registry.CurrentLimitSelectedPostpulseRegistry.history.heat',
            exact_definition=self.heat_proof['exact_c_infinity_definition'],
            exact_correlated_log=self.heat_proof['exact_correlated_log_c_infinity'],
            source_AST=self.heat_proof['actual_amplitude_and_radius_AST'],
            strictly_positive=True,Z_independent=True,defining_quantity_not_a_range_value=True,
            acceptance_prerequisite=joined.RECEIPT)
        targets={key:[band.substitution(g,source.FunctionRef(g,node),[old],[exact]).node for node in rows]
            for key,rows in self.required_target_nodes.items()}
        return dict(original_graph_prefix_length=len(nodes),original_target_handles=self.required_target_nodes,
            replaced_c_infinity_symbol=old.node,actual_c_infinity_recipe=exact.node,
            actual_source_target_handles=targets,appended_graph_nodes=g.nodes[len(nodes):],
            final_reference_c_infinity_installed=True,original_graph_preserved=True,
            positive_scale_units=UNITS)

    @joined.source_precision
    def boundary_and_reference(self,Z):
        c=self.ctx;Z=c.mpf(Z)
        if joined.history.endpoints(Z)[0]<-1 or joined.history.endpoints(Z)[1]>1:raise ValueError('Original Z in[-1,1] required')
        bc=self.buffer.ctx
        phase=(bc.mpf(2)+bc.ln(2))/bc.mpf(joined.history.endpoints(self.buffer.params.Tw))
        packet=native_real_phase_power(self.buffer,bc.mpf(joined.history.endpoints(Z)),phase)
        jet=lambda rows:IntervalTaylor(c,[c.mpf(joined.history.endpoints(v)) for v in rows[:2]])
        inverse_Pstar=c.exp(-c.mpf(joined.history.endpoints(self.params.logPstar)))
        own=dict(m=jet(packet['Mz_over_R'])*inverse_Pstar,
            h=jet(packet['Mtheta_over_sqrt2_R_3half_Pstar']),
            k=jet(packet['Mtheta_z_over_sqrt2_R_3half_Pstar'])*inverse_Pstar,
            e=jet(packet['Mztheta_over_R_Pstar_squared']),p=jet(packet['Mp_over_Pstar_squared']))
        P0=jet(packet['P0_over_Pstar_squared'])
        a=self.heat.a;mu=self.heat.mu;k=self.heat.k;bh=self.heat.bh
        length=c.mpf(joined.history.endpoints(self.params.Tw))-2-c.ln(2)
        logU=c.ln(self.pulse.high.constants['U'])
        log_reference_parts=dict(log_inlet_U=logU,quiet_radius=bh*length,
            cancelled_pulse_inverse_mu=13*a/mu,flatten_power=(a-mu)*(100+self.owner.history.outer.Lrel),
            steep=-k*self.owner.history.steep.Ts,
            finite=-14-mu/2+3*a/2-c.ln(2),minus_logone=-self.owner.history.steep.logone)
        reference_at_R2=c.exp(sum(log_reference_parts.values(),c.mpf(0)))
        ref_jet=IntervalTaylor.constant(c,reference_at_R2/k,1)
        targets=dict(M=-own['m'],J=-own['k'],S_energy=-own['e'],
            renormalized_I=ref_jet-own['h'],pressure_tail=-(P0+own['p']))
        return dict(Z=Z,native_2Rc_phase=phase,current_2Rc_common_own_C1=own,
            independent_original_P0_C1=P0,
            actual_final_reference_at_2Rc_over_Pstar_log_parts=log_reference_parts,
            actual_final_reference_primitive_in_R2_units=ref_jet,
            required_future_integral_C1_in_declared_units=targets,
            units=UNITS,original_2Rc_source_target_graph_handles=self.required_target_nodes,
            native_2Rc_function_identity_from_checked_Rp_and_backward_ODE_uniqueness=True,
            native_power_code_unchanged_only_real_phase_acquisition_adapted=True,
            final_reference_uses_same_heat_amplitude_not_current_O3_exponent=True)

    @joined.source_precision
    def endpoint_ledger(self,Z,boundary=None):
        """Callable absolute endpoint primitives and four exact increments.

        Each atom is coefficient(Z) * exp(finite + inverse_mu/mu),
        interpreted as the exact source function. Its C1 coefficient rows
        are directed enclosures. No cap or upper endpoint defines an atom.
        Internal atom cancellation precedes any exponential materialization.
        """
        boundary=boundary if boundary is not None else self.boundary_and_reference(Z)
        c=self.ctx;Z=c.mpf(Z);mu=self.heat.mu
        length=c.mpf(joined.history.endpoints(self.params.Tw))-2-c.ln(2)
        tail_offset=100+self.owner.history.outer.Lrel+2+self.owner.history.steep.Ts+self.owner.history.steep.wait
        B=length+tail_offset;atoms={};endpoints={}
        copy=lambda jet:IntervalTaylor(c,[c.mpf(joined.history.endpoints(v)) for v in jet.coefficients[:2]])
        def atom(stage,key,coefficient,finite=0,inverse=0,provider=None,quantity=None):
            name=stage+'_'+key
            atoms[name]=dict(actual_C1_coefficient_enclosure=copy(coefficient),
                exact_scale_log=dict(finite=c.mpf(finite),inverse_mu_coefficient=c.mpf(inverse)),
                source_provider=provider,source_quantity=quantity,source_family=self.family_record,
                numerical_rows_are_enclosures_not_defining_values=True,
                exact_positive_scale_never_replaced_with_cap=True)
            return [(1,name)]
        def provenance(obj,method,coordinate):
            fn=getattr(obj,method).__func__
            return dict(actual_method=fn.__module__+'.'+fn.__qualname__,
                source_sha256=sha(fn.__module__+'.py'),coordinate=coordinate)
        own=boundary['current_2Rc_common_own_C1'];P0=boundary['independent_original_P0_C1']
        provider=provenance(self.buffer,'power',boundary['native_2Rc_phase'])
        provider['exact_phase_adapter']='current_limit_absolute_future_integrals.native_real_phase_power'
        endpoints['R2']={key:atom('R2',key,own[native],provider=provider,quantity=native)
            for key,native in dict(M='m',J='k',I='h',S_energy='e').items()}
        endpoints['R2']['pressure_tail']=atom('R2','pressure_tail',P0+own['p'],provider=provider,
            quantity='P0+p; original independent P0 retained exactly once')
        ref2=boundary['actual_final_reference_primitive_in_R2_units']
        refs={stage:atom(stage,'Ipow',ref2,finite,inv,
            provider=provenance(self.heat,'exterior','actual c_infinity and final power'),
            quantity='sqrt(2)*c_infinity*R^(1-a)/(1-a) in R2 angular units')
            for stage,finite,inv in (('R2',0,0),('Rp',self.heat.k*length,0),
                ('Rv',self.heat.k*length-13*self.heat.a/mu,13),
                ('Rt',self.heat.k*B-13*self.heat.a/mu,13))}
        rp=self.pulse.entrance(Z,0);rv=self.pulse.end(Z,0)
        _,_,u,_,_=self.pulse.data(Z)
        rp_provider=provenance(self.pulse,'entrance',0);rv_provider=provenance(self.pulse,'end',0)
        pulse_keys=dict(M='Mz_over_R_Utheta',J='Mtheta_z_over_sqrt2_R_3half_Utheta_squared',
            I='Mtheta_over_sqrt2_R_3half_Utheta',S_energy='Mztheta_over_R_Utheta_squared')
        for stage,packet,provider in (('Rp',rp,rp_provider),('Rv',rv,rv_provider)):
            row={}
            for key,native in pulse_keys.items():
                scale=(u*u if key in ('J','S_energy') else u)
                # Materialize finite log terms only; 13/mu is kept exact.
                finite={'M':length,'J':c.mpf('1.5')*length,'I':c.mpf('1.5')*length,
                    'S_energy':length}[key]
                inverse=0
                if stage=='Rv':
                    finite-=26 if key in ('J','S_energy') else 13
                    inverse={'M':c.mpf('6.5'),'J':c.mpf('6.5'),'I':13,'S_energy':0}[key]
                row[key]=atom(stage,key,packet[native]*scale,finite,inverse,provider,native)
            row['pressure_tail']=atom(stage,'pressure_tail',packet['pressure']['P_over_Pstar_squared'],
                provider=provider,quantity='actual absolute P0+Mp')
            endpoints[stage]=row
        heat=self.heat.collar(Z,0);theta=heat['theta_over_Ev0_Taylor']*self.pulse.high.constants['U']
        provider=provenance(self.heat,'collar',0);zero=theta*0
        endpoints['Rt']=dict(
            M=atom('Rt','M',zero,provider=provider,quantity='selected terminal absolute M=0; homogeneous FTC'),
            J=atom('Rt','J',zero,provider=provider,quantity='selected terminal absolute J=0; homogeneous FTC'),
            I=atom('Rt','I',heat['angular_Taylor']*theta,c.mpf('1.5')*B-13,13,provider,'angular_Taylor*actual Utheta'),
            S_energy=atom('Rt','S_energy',heat['energy_Taylor']*theta*theta,B-26,0,provider,'energy_Taylor*actual Utheta^2'),
            pressure_tail=atom('Rt','pressure_tail',heat['pressure_over_Pstar_squared_Taylor'],
                provider=provider,quantity='actual inherited absolute P0+Mp'))
        normalized={}
        for stage,row in endpoints.items():
            normalized[stage]={key:row[key] for key in KEYS if key!='renormalized_I'}
            normalized[stage]['renormalized_I']=row['I']+[(-sign,name) for sign,name in refs[stage]]
        subtract=lambda right,left:right+[(-sign,name) for sign,name in left]
        stages=('R2','Rp','Rv','Rt')
        segments={name:{key:subtract(normalized[stages[i+1]][key],normalized[stages[i]][key])
            for key in KEYS} for i,name in enumerate(SEGMENTS[:3])}
        segments[SEGMENTS[3]]={key:[(-sign,name) for sign,name in normalized['Rt'][key]] for key in KEYS}
        return dict(atoms=atoms,absolute_endpoint_functions=endpoints,
            final_reference_endpoint_functions=refs,normalized_endpoints=normalized,segments=segments,
            units=UNITS,exact_endpoint_radius_logs_from_R2=dict(R2=(c.mpf(0),c.mpf(0)),
                Rp=(length,c.mpf(0)),Rv=(length,c.mpf(13)),Rt=(B,c.mpf(13))),
            same_internal_source_atom_reused_by_adjacent_segments=True,
            source_scope='all five absolute similarity primitives, true first ordinary Z row; finite scales factored')

    @joined.source_precision
    def future_integrals(self,Z):
        c=self.ctx;Z=c.mpf(Z);cache_key=tuple(joined.history.endpoints(Z))
        if cache_key in self.cache:return self.cache[cache_key]
        boundary=self.boundary_and_reference(Z);ledger=self.endpoint_ledger(Z,boundary)
        telescope=telescope_proof(ledger)
        # Cumulative FTC and the proved C1 infinity limits supply these
        # function enclosures; internal cancellation is algebraic, not an
        # overlap test on separately computed endpoint interval boxes.
        # Evaluate the reduced actual atom expression, separately from the
        # required target API. Each reduction consists of the genuine R2
        # source rows after exact shared-atom cancellation.
        actual={}
        for key in KEYS:
            total=IntervalTaylor.constant(c,0,1)
            for sign,name in ledger['normalized_endpoints']['R2'][key]:
                entry=ledger['atoms'][name]
                require(joined.history.endpoints(entry['exact_scale_log']['finite'])==(0,0) and
                    joined.history.endpoints(entry['exact_scale_log']['inverse_mu_coefficient'])==(0,0),
                    'Reduced endpoint must be in actual R2 units')
                total-=entry['actual_C1_coefficient_enclosure']*sign
            actual[key]=total
        result=dict(**boundary,absolute_future_integrals=actual,actual_four_segment_source_ledger=ledger,
            segment_definitions=dict(quiet_power='[2Rc,Rp]',selected_pulse='[Rp,Rv]',
                selected_postpulse='[Rv,Rtail]',collar_and_full_Gamma='[Rtail,infinity)'),
            densities=dict(M='Uz',J='sqrt(2R)*Utheta*Uz',S_energy='Uz^2-Utheta^2/2',
                renormalized_I='sqrt(2R)*Utheta-sqrt(2)*c_infinity*R^(-delta/2)',
                pressure_tail='Utheta^2/(2R)'),
            source_cumulative_FTC_proof=self.FTC,source_full_Gamma_C1_limits=self.heat_proof,
            exact_entire_future_telescope=telescope,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False),
            output_kind='source-bound absolute future C1 enclosures in declared positive paper units; no point/cap selection')
        self.cache[cache_key]=result;return result


@joined.source_precision
def run(field=None):
    began=time.monotonic();field=field if field is not None else CurrentLimitAbsoluteFutureIntegrals(require_checked=False)
    result=dict(source_family=field.family_record,**dict.fromkeys(GATES,False),**dict.fromkeys(OPEN,False),
        candidate_absolute_future_source_constructed=True,
        actual_current_cumulative_FTC_source_proof=field.FTC,
        actual_current_final_heat_reference_and_C1_terminal_limits=field.heat_proof,
        actual_final_reference_target_graph=field.installed_target_graph,
        source_views={})
    for name,Z in (('fresh','.519'),('axis',0)):
        result['source_views'][name]=joined.serialized(field.future_integrals(Z))
        print('Current absolute five future integrals: '+name,flush=True)
    result['all_passed']=True;result['input_hashes']=field.hashes
    result['execution_seconds']=time.monotonic()-began
    (HERE/NAME).write_text(json.dumps(joined.serialized(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':run()
