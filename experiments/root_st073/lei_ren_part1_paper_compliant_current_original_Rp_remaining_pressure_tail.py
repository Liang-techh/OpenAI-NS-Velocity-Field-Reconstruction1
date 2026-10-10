"""Original pulse pressure as a correlated negative remaining integral.

The original P0 and cumulative primitive are identified before any box
arithmetic. Positive late atoms stay exact functions with directed bounds;
they are never discarded or used to patch the pressure datum.
"""
import copy
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
from types import MethodType
import time

import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rp_centered_scale_arithmetic as centered
import lei_ren_part1_paper_compliant_current_pressure_terminal_closure as terminal_source
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision

refined, signed, correlated = centered.refined, centered.signed, centered.correlated
box, physical = refined.box, refined.physical
HERE, PREFIX, sha, ends = centered.HERE, centered.PREFIX, centered.sha, centered.ends
NAME=PREFIX+'current_original_Rp_remaining_pressure_tail.json.gz'
RECEIPT=PREFIX+'current_original_Rp_remaining_pressure_tail_check.json'
GATES=('current_original_Rp_exact_remaining_pulse_pressure_function_installed',
       'current_original_Rp_positive_late_pressure_atoms_bounded',
       'current_original_Rp_correlated_tail_pressure_physical_rows_installed')
OPEN=centered.OPEN
LATE=('z_flatten','power_buffer_rel','steep_transition_in','steep_power',
      'steep_transition_out','waiting','heat_collar','exterior_power_tail')


class OriginalRemainingPulsePressureFunction:
    """Pure exact witness with the original raw waiting root selected."""
    def __init__(self):
        self.witness=terminal_source.ExactOriginalPreheatPressureOperator()
        w=self.witness;p=w.partition;self.xi=s.Symbol('xi',real=True)
        rate=1+2*p['mu'];t=p['t'];q=p['q']
        U2=s.exp(2*w.native_log_up)
        self.pulse_partial=U2*(s.exp(-rate*self.xi/p['mu'])-
            s.exp(-rate*13/p['mu']))/(2*rate*q*q)
        self.pulse_prefix=U2*(1-s.exp(-rate*self.xi/p['mu']))/(2*rate*q*q)
        self.late_atoms={name:w.stage_integral(name,True) for name in LATE}
        self.remaining=self.pulse_partial+sum(self.late_atoms.values(),s.Integer(0))
        self.incoming=sum((w.stage_integral(name,True) for name in w.native_densities
                          if name not in LATE+('pulse_reserved',)),s.Integer(0))
        self.partial_cumulative=self.incoming+self.pulse_prefix
        complete=U2*(1-s.exp(-rate*13/p['mu']))/(2*rate*q*q)
        assert s.expand(self.pulse_partial+self.pulse_prefix-complete)==0
        assert s.simplify(s.diff(self.pulse_prefix,self.xi)*p['mu']-
            w.native_densities['pulse_reserved'].subs(t,self.xi/p['mu']))==0
        self.proof=dict(pulse_partition_exact=True,pulse_partial_FTC_exact=True,
            full_P0_identity_consumed_separately=True,
            identity='P0 + incoming_prefix + partial_pulse_integral = -remaining_pressure_tail',
            all_late_atoms=list(LATE),raw_waiting_root=str(w.raw_waiting_root),
            flatten_variable_q_exponent_retained=True,
            exact_pressure_units='Pstar^2; normalized witness is pressure/Pstar^2',
            domain='real Z in [-1,1], 0 <= xi < 13, mu > 0; original raw waiting root',
            no_interval_difference_defines_a_source_function=True)

    def remaining_pressure_after_pulse(self,Z,xi):
        Z,xi=Fraction(Z),Fraction(xi)
        if not -1<=Z<=1 or not 0<=xi<13:raise ValueError('Original pulse domain required')
        p=self.witness.partition
        return self.remaining.subs({p['z']:s.Rational(Z.numerator,Z.denominator),
            self.xi:s.Rational(xi.numerator,xi.denominator),p['W']:self.witness.raw_waiting_root})

    def snapshot(self):
        w=self.witness
        return (s.srepr(self.remaining),s.srepr(self.partial_cumulative),
                tuple((name,s.srepr(v)) for name,v in w.native_densities.items()),
                tuple((name,str(v)) for name,v in w.partition['domains'].items()),
                s.srepr(w.raw_waiting_root))


def tail_majorant_proof():
    """Analytic bounds on exact later integrals, independent of quadrature."""
    a,rho,q0=s.symbols('a rho q0',positive=True)
    j=s.symbols('j',integer=True,nonnegative=True)
    # Re(1+(Z+w)^2) >= 1-rho^2 on every real center |Z|<=1.
    # For 0<=sigma<=1, |q^(2sigma-2) 2^(-2sigma)| <= q0^-2.
    # Post-flatten logarithmic slope is <= -1/2. The raw collar factor
    # is in [1-eps,1]; its sole amplitude boost is 1/(1-eps).
    return dict(rho='1/4',complex_q_real_lower='15/16',
        flatten_density='exp(-p*t)*q^(2*sigma(t/100)-2)*2^(-2*sigma(t/100))/2',
        flatten_integral_upper='(15/16)^(-2)/(2*p)',
        late_Z_independent_integral_upper='exp(-100*p)/(8*(1-epsilon)^2)',
        whole_flatten_interval_not_truncated=True,
        all_seven_postflatten_stages_controlled_by_one_integrable_envelope=True,
        late_infinite_tail_not_dropped=True,
        Cauchy_ordinary_Taylor_upper='flatten_upper * (1/4)^(-j); j=1..5',
        postflatten_axial_derivatives_exactly_zero=True,
        source_slope_bounds=dict(power='-(1/2+mu)',steep_in='-(1/2+mu)-(1-mu)*sigma',
            steep_power='-3/2',steep_out='-3/2+(1-delta/2)*sigma',waiting='-(1+delta)/2'),
        raw_heat_collar_factor='1-epsilon*(1-sigma+sigma*phi)',
        source_parameters='0<mu<1, 0<delta<1, 0<epsilon<1/2; raw W>0')


def global_postflatten_envelope(witness):
    """Identify each exact density with one disjoint global-time envelope."""
    p=witness.partition;t=p['t'];mu=p['mu'];delta=p['delta'];eps=p['epsilon']
    L,Ts,W=p['L'],p['Ts'],p['W'];J=p['J'](t)
    rr,JJ,gap=s.symbols('r J gap',nonnegative=True)
    common=-2*mu*(L+1)-rr
    cones=dict(power_buffer_rel=-2*mu*t,
        steep_transition_in=-2*mu*(L+t)-2*rr*JJ,
        steep_power=common-2*t,
        steep_transition_out=common-2*Ts-2*gap-delta*JJ,
        waiting=common-2*Ts-1-delta/2-delta*t,
        heat_collar=common-2*Ts-1-delta/2-delta*(W+t),
        exterior_power_tail=common-2*Ts-1-delta/2-delta*(W+t))
    offsets=dict(power_buffer_rel=0,steep_transition_in=L,steep_power=L+1,
        steep_transition_out=L+1+Ts,waiting=L+Ts+2,
        heat_collar=L+Ts+2+W,exterior_power_tail=L+Ts+2+W)
    prefactor=s.exp(witness.native_log_ev2-(1+2*mu)*100)/8
    bracket=1-eps*(1-p['sigma'](t)+p['sigma'](t)*p['phi'](t))
    identities={};domains=[]
    for name,cone in cones.items():
        coefficients=s.Poly(s.expand(-cone),mu,L,Ts,W,t,rr,JJ,gap,delta).coeffs()
        if not all(v>=0 for v in coefficients):raise ArithmeticError('Nonpositive global density cone required')
        exponent=cone.subs({rr:1-mu,JJ:J,gap:t-J})
        factor=bracket**2/(1-eps)**2 if name=='heat_collar' else 1/(1-eps)**2 if name=='exterior_power_tail' else 1
        expected=prefactor*s.exp(-(offsets[name]+t)+exponent)*factor
        difference=s.expand_log(s.expand_power_exp((witness.native_densities[name]-expected).rewrite(s.exp)),force=True)
        if s.simplify(difference,doit=False)!=0:raise ArithmeticError('Original global density envelope differs: '+name)
        identities[name]=dict(exact_native_density_identity=True,global_offset=str(offsets[name]),
            nonpositive_exponent=str(cone),positive_cone_coefficient_test=True,
            prefactor=str(prefactor),raw_collar_factor_retained=name=='heat_collar')
        a,b=p['domains'][name]
        domains.append((s.sympify(offsets[name])+a,s.sympify(offsets[name])+b))
    if domains[0][0]!=0 or domains[-1][1]!=s.oo:raise ArithmeticError('Full global postflatten interval required')
    for left,right in zip(domains,domains[1:]):
        if s.simplify(left[1]-right[0])!=0:raise ArithmeticError('Disjoint consecutive original stages required')
    return dict(exact_stage_density_identities=identities,
        global_stage_domains=[list(map(str,row)) for row in domains],
        positive_lengths='L=-30log(mu)>0; Ts=4log(2/delta)>0; original raw W>0',
        primitive_gap='0<=J_sigma(t)<=t follows from 0<=sigma<=1',
        disjoint_global_domain='[0,infinity); collar [0,3] and exterior [3,infinity) share their tail origin',
        normalization='Uend^2 * exp(-100*p)/8; density normalization, not amplitude',
        envelope='Uend^2 * exp(-100*p)/(8*(1-epsilon)^2) * exp(-global_time)',
        exact_envelope_integral='Uend^2 * exp(-100*p)/(8*(1-epsilon)^2)',
        all_seven_native_stage_prefactors_and_origins_identified=True)


def original_cutoff_inequalities():
    """Bind the original formulas and derive their full-domain bounds."""
    source=terminal_source.bind_original_master_source()
    a,b,D=s.symbols('a b D',positive=True)
    ratio=a/(a+b)
    if s.ask(s.Q.positive(ratio)) is not True or s.ask(s.Q.positive(s.simplify(1-ratio))) is not True:
        raise ArithmeticError('Positive edge ratio bounds required')
    if s.ask(s.Q.nonpositive(-4/D**2)) is not True:
        raise ArithmeticError('Original flat phi exponent must be nonpositive')
    return dict(actual_original_formula_bindings=source,
        sigma_inside_support='a/(a+b), a>0, b>0; both sigma and 1-sigma positive',
        sigma_outside_support='original flat extensions 0 on x<=0 and 1 on x>=1',
        sigma_full_domain_bound='0<=sigma<=1',
        phi_full_collar_bound='phi=exp(-4/(3-t)^2) for t<3, otherwise0; 0<=phi<=1 by exp monotonicity',
        primitive_bound='J(t)=integral_0^t sigma; integrating 0<=sigma<=1 gives 0<=J(t)<=t for t>=0',
        source_ratio_and_phi_exponent_signs_checked=True,
        no_sampled_cutoff_values_used_as_uniform_bounds=True)


class CurrentOriginalRpRemainingPressureTail:
    @source_precision
    def __init__(self,before,require_checked=True):
        if type(before) is not centered.CurrentOriginalRpCenteredScaleArithmetic or not before.acceptance_loaded:
            raise ValueError('Accepted exact-centered current source required')
        self.before,self.refined=before,before.before
        self.amplitude,self.ctx,self.graph=before.amplitude,self.refined.ctx,before.graph
        self.family_record=before.family_record
        self.terminal=self.refined.raw.closed.pressure
        self.function=OriginalRemainingPulsePressureFunction()
        self._function_snapshot=self.function.snapshot()
        self._terminal_proof_snapshot=json.dumps(self.terminal.proof,sort_keys=True)
        self._parameter_sources=dict(mu=self.terminal.raw.flat.mu,
            delta=self.terminal.raw.flat.delta,
            epsilon=self.terminal.raw.flat.inlet.datum.parameters.epsilon,
            waiting=self.terminal.exact.repair.angular.waiting,
            log_mu=self.terminal.exact.repair.params.log_mu,Ts=self.terminal.exact.repair.params.Ts)
        self._parameter_snapshots={name:value._mpi_ for name,value in self._parameter_sources.items()}
        self.hashes=dict(before.hashes)
        for name in (centered.NAME,centered.RECEIPT,Path(__file__).name):self.hashes[name]=sha(name)
        self.acceptance_loaded=False;self._sources={}
        self.majorant_proof=tail_majorant_proof()
        self.majorant_proof['original_disjoint_global_time_envelope']=global_postflatten_envelope(self.function.witness)
        self.majorant_proof['original_cutoff_inequalities']=original_cutoff_inequalities()
        self.L_function=self.amplitude.radius.functions['Lrel']
        self.logmu_function=self.amplitude.radius.functions['logmu']
        self.expected_L=self.graph.mul(self.graph.constant(-30),self.logmu_function)
        bare=correlated.locator.CancelledGraphBounds(self.graph,self.ctx,{},())
        if bare.polynomial(self.graph.sub(self.L_function,self.expected_L)):
            raise ValueError('Original same-graph Lrel=-30logmu source identity required')
        Lbox=self.amplitude.reader().at(self.L_function)
        if ends(Lbox)[0]<=0:raise ValueError('Original Lrel must be positive')
        self._L_snapshot=Lbox._mpi_
        self.L_identity=self.graph.node('current_original_Rp_tail_exact_positive_Lrel',
            actual_Lrel_function=self.L_function.node,original_logmu_function=self.logmu_function.node,
            exact_minus_30_logmu_function=self.expected_L.node,
            same_graph_polynomial_identity=True,source_L_positive_from_logmu_negative=True,
            original_schedule_binding=self.terminal.proof['original_master_source']['original_axial_flatten_and_tail_AST_bindings']['PaperOuterSchedule.__init__:self.y_rel'])
        self._L_identity_snapshot=copy.deepcopy(self.graph.nodes[self.L_identity.node])
        self.proof=self.graph.node('current_original_Rp_exact_remaining_pressure_tail_identity',
            original_P0_binding=self.refined.before.before.binding,
            original_complete_native_integral=self.terminal.proof['exact_native_complete_pressure_integral'],
            original_raw_waiting_root=str(self.function.witness.raw_waiting_root),
            partial_pulse_integral=str(self.function.pulse_partial),
            exact_late_integrals={name:str(v) for name,v in self.function.late_atoms.items()},
            exact_source_function_identity=self.function.proof,
            directed_late_integral_majorant=self.majorant_proof,
            exact_positive_Lrel_identity=self.L_identity.node,
            existing_pressure_radial_FTC_rows_preserved=True,
            algebraic_rebase='Pstar^2 pressure = Pstar^2 pulse_F^2 times same-source tail coefficient',
            no_datum_or_forcing_replaced=True)
        self._proof_snapshot=copy.deepcopy(self.graph.nodes[self.proof.node])
        self.assert_graph()
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or not all(receipt[k] for k in GATES) or \
                    any(receipt[k] for k in OPEN) or receipt['source_family']!=self.family_record:
                raise ValueError('Current remaining pressure receipt/scope differs')
            box.pulse.radius.post.selected.inlet.add_hashes(self.hashes,receipt['input_hashes'])
            self.acceptance_loaded=True

    @source_precision
    def assert_graph(self):
        self.before.assert_graph()
        p=self.terminal;w=self.function.witness
        parameters=dict(mu=p.raw.flat.mu,delta=p.raw.flat.delta,
            epsilon=p.raw.flat.inlet.datum.parameters.epsilon,waiting=p.exact.repair.angular.waiting,
            log_mu=p.exact.repair.params.log_mu,Ts=p.exact.repair.params.Ts)
        checks=dict(accepted_current_terminal=p.acceptance_loaded and p.proof['passed'],
            same_original_complete_integral=str(w.integral(True))==p.proof['exact_native_complete_pressure_integral'],
            same_original_P0_function=str(w.datum())==p.proof['exact_original_P0_function'],
            same_current_selected_P0=p.raw.flat.inlet.datum is self.refined.selected.datum,
            same_current_family=p.family_record==self.family_record,
            exact_function_not_mutated=self.function.snapshot()==self._function_snapshot,
            complete_original_identity_retained=json.dumps(p.proof,sort_keys=True)==self._terminal_proof_snapshot,
            actual_complete_identity=p.proof['identities']['original_exact_P0_plus_native_exact_complete_integral'],
            actual_axial5_identity=all(p.proof['identities']['original_P0_plus_native_raw_integral_axial_'+str(j)] for j in range(6)),
            same_original_pressure_FTC=p.raw.proof['identities']['native_pulse_pressure_density_primitive'],
            same_bound_source_parameters=all(parameters[name]._mpi_==value for name,value in self._parameter_snapshots.items()),
            original_positive_parameter_domain=all(0<ends(parameters[name])[0]<=ends(parameters[name])[1]<1
                for name in ('mu','delta')) and 0<ends(parameters['epsilon'])[0]<=ends(parameters['epsilon'])[1]<self.ctx.mpf('0.5')
                and ends(parameters['waiting'])[0]>0 and ends(parameters['log_mu'])[1]<0 and ends(parameters['Ts'])[0]>0,
            same_original_mu_binding=p.raw.flat.mu._mpi_==self.ctx.mpf(ends(self.amplitude.reader().bindings[self.amplitude.mu.node]))._mpi_,
            same_proved_original_U0_box=self.refined.constants['U']._mpi_==self.ctx.mpf(ends(self.amplitude.U0_box))._mpi_,
            same_original_Lrel_functions=self.L_function==self.amplitude.radius.functions['Lrel']
                and self.logmu_function==self.amplitude.radius.functions['logmu'],
            unchanged_exact_Lrel_proof=self.graph.nodes[self.L_identity.node]==self._L_identity_snapshot,
            unchanged_positive_Lrel=self.amplitude.reader().at(self.L_function)._mpi_==self._L_snapshot
                and ends(self.amplitude.reader().at(self.L_function))[0]>0,
            unchanged_exact_tail_proof_node=self.graph.nodes[self.proof.node]==self._proof_snapshot)
        if not all(checks.values()):raise ValueError('Current original remaining pressure source differs: '+str(checks))
        return checks

    @source_precision
    def tail_rows(self,delivery):
        self.assert_graph()
        reader,token,request=self.amplitude.before._validate_delivery(delivery)
        if request.chart not in ('pulse_main','pulse_exit','pulse_gap'):
            raise ValueError('Current tail adapter owns main/exit/gap xi-coordinate charts')
        jac=box.pulse.radius.FunctionRef(self.graph,delivery['actual_source_view']['geometry']['native_to_log_radius_jacobian'])
        expected_jac=self.graph.quotient(self.graph.one,self.amplitude.mu,'original positive mu')
        if reader.polynomial(self.graph.sub(jac,expected_jac)):
            raise ValueError('Original xi-to-logR Jacobian must be exactly 1/mu')
        jacobian_proof=self.graph.node('current_original_Rp_tail_native_logR_FTC_identity',
            actual_native_jacobian=jac.node,exact_reciprocal_mu=expected_jac.node,
            proof='exact same-graph polynomial reciprocal normalization; mu*d_xi=d_logR',
            pressure_F2_definition='exp(-(1+2mu)*xi/mu)',
            existing_radial_product_rows_retained=True)
        c=self.ctx;mu=c.mpf(self.amplitude.reader(delivery).bindings[self.amplitude.mu.node])
        xi=c.mpf(ends(token.enclosure));Z=c.mpf(ends(reader.at(delivery['physical_inverse']['coordinate_functions']['Z'])))
        p=1+2*mu;gap=13-xi
        if ends(gap)[0]<=0:raise ValueError('Positive exact remaining pulse interval required')
        logeta=-p*gap/mu
        if ends(logeta)[1]>-signed.EXP_LOG_LIMIT:
            raise ValueError('This bounded tail adapter requires its proved large-gap ratio guard')
        eta=c.mpf([0,ends(c.exp(-signed.EXP_LOG_LIMIT))[1]])
        epsilon=c.mpf(ends(self.terminal.raw.flat.inlet.datum.parameters.epsilon))
        if not 0<ends(mu)[0]<=ends(mu)[1]<1 or not 0<ends(epsilon)[0]<=ends(epsilon)[1]<c.mpf('0.5'):
            raise ValueError('Original positive source parameters required')
        rho=c.mpf(1)/4;qlo=1-rho*rho
        flatten_upper=1/(2*p*qlo*qlo)
        post_upper=c.exp(-100*p)/(8*(1-epsilon)**2)
        late0=flatten_upper+post_upper
        late=IntervalTaylor(c,[c.mpf([0,ends(late0)[1]])]+[
            c.mpf([-ends(flatten_upper/rho**j)[1],ends(flatten_upper/rho**j)[1]]) for j in range(1,6)])
        q=IntervalTaylor(c,[1+Z*Z,2*Z,1,0,0,0]);invq2=q**(-2)
        U=self.refined.constants['U'];leading=invq2*(U*U/(2*p))
        coefficients=-leading+(leading-late*(U*U))*eta
        logF=self.refined.raw.scale(token.chart,token,(0,2,0))
        eta_function=self.graph.unary('exp',self.graph.quotient(
            self.graph.mul(self.graph.constant(-1),self.graph.add(self.graph.one,
                self.graph.mul(self.graph.constant(2),self.amplitude.mu)),
                self.graph.sub(self.graph.constant(13),token.function)),self.amplitude.mu,'original positive mu'))
        function=self.graph.node('current_original_Rp_pointed_remaining_preheat_pressure',
            original_tail_identity=self.proof.node,exact_Z_function=delivery['physical_inverse']['coordinate_functions']['Z'].node,
            exact_xi_function=token.function.node,exact_positive_ratio=eta_function.node,
            original_U0_function=self.amplitude.raw.U0.node,
            exact_pressure_function='P/Pstar^2=-U0^2*exp(-p*xi/mu)*(q^-2/(2p)*(1-eta)+eta*late_integral)',
            late_integral_function='all eight exact original later stage integrals divided by U0^2*exp(-13p/mu)',
            source_function_is_not_defined_by_interval_coefficients=True)
        return dict(coefficients=coefficients,leading_pressure_coefficient=leading,
            late_integral_Taylor_enclosure=late,flatten_majorant=flatten_upper,postflatten_majorant=post_upper,
            positive_ratio_enclosure=eta,exact_positive_ratio_function=eta_function.node,
            ratio_log_enclosure=logeta,exact_pressure_function=function.node,
            exact_native_logR_FTC_identity=jacobian_proof.node,
            retained_exact_pulse_F2_log_parts={k:v.node for k,v in logF.items()},
            ratio_upper_is_bound_not_value=True,all_late_positive_integrals_retained=True,
            majorant_proof=self.majorant_proof,original_P0_source=self.refined.before.before.source_binding(delivery))

    @source_precision
    def source(self,delivery):
        self.assert_graph();reader,token,request=self.amplitude.before._validate_delivery(delivery)
        record=self._sources.get(id(delivery))
        if record is not None:
            if record[0] is not delivery or record[3]!=self.amplitude.before._fingerprint(box.report(record[1])) or \
                    any(self.graph.nodes[node]!=value for node,value in record[4].items()):
                raise ValueError('Unchanged live remaining-tail view required')
            return record[1]
        original=self.refined.source(delivery);tail=self.tail_rows(delivery)
        view=dict(original);base=dict(original['original_factorized_values'])
        base['pressure']=self.refined.raw.factor('pressure',token.chart,token,(0,2,2),tail['coefficients'])
        view['original_factorized_values']=base
        for section,native in (('log_radius_mixed_rows',False),('native_coordinate_mixed_rows',True)):
            view[section]=dict(original[section]);rows=dict(original[section]['pressure'])
            for j in range(5):
                coefficient=IntervalTaylor.constant(self.ctx,tail['coefficients'][j]*math.factorial(j),0)
                rows[('n' if native else 'y')+'0_Z'+str(j)]=self.refined.transport.factor(
                    'pressure',token.chart,token,view['geometry'],(0,2,2),coefficient,0,j,native)
            view[section]['pressure']=rows
        view['source_equivalent_remaining_pressure_tail']=tail
        node_snapshots={node:copy.deepcopy(self.graph.nodes[node]) for node in
            (tail['exact_pressure_function'],tail['exact_positive_ratio_function'])}
        self._sources[id(delivery)]=(delivery,view,tail,self.amplitude.before._fingerprint(box.report(view)),node_snapshots)
        return view

    @source_precision
    def physical_rows(self,delivery):
        proxy=box._Proxy(self.refined);owner=self
        proxy.source=lambda observed:owner.source(observed)
        return refined.CurrentOriginalRpRefinedPulseCoefficients.physical_rows(proxy,delivery)

    @source_precision
    def evaluate(self,delivery,relative_width_target='1/100000000'):
        self.source(delivery)
        proxy=box._Proxy(self.before);sourceproxy=box._Proxy(self.refined)
        sourceproxy.physical_rows=self.physical_rows;proxy.before=sourceproxy
        result=centered.CurrentOriginalRpCenteredScaleArithmetic.evaluate(proxy,delivery,relative_width_target)
        result.update(source_equivalent_remaining_pressure_tail=self.source(delivery)['source_equivalent_remaining_pressure_tail'],
            exact_original_pressure_tail_identity=self.proof.node,source_coefficient_provider='accepted_refined_velocity_and_correlated_pressure_tail',
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))
        return result

    @source_precision
    def velocity_pressure(self,delivery,relative_width_target='1/100000000'):
        _,_,request=self.amplitude.before._validate_delivery(delivery)
        view=self.evaluate(delivery,relative_width_target);rows=view['physical_value_rows']['Cartesian_spatial_rows']
        return dict(values={name:rows[key]['x0_y0_z0'] for name,key in (('u','ux'),('v','uy'),('w','uz'),('p','p'))},
            physical_coordinates={key:value.node for key,value in request.forward_coordinates.items()},
            source_family=self.family_record,exact_original_pressure_tail_identity=self.proof.node,
            original_P0_source_binding=view['original_P0_source_binding'],
            full_certified_physical_accuracy=False,unrestricted_physical_point_API=False,
            **dict.fromkeys(GATES,self.acceptance_loaded),**dict.fromkeys(OPEN,False))


@source_precision
def run(before,deliveries,relative_width_target='1/1000'):
    began=time.monotonic();owner=CurrentOriginalRpRemainingPressureTail(before,require_checked=False)
    views={name:owner.evaluate(delivery,relative_width_target) for name,delivery in deliveries.items()}
    result=dict(source_family=owner.family_record,exact_remaining_pressure_identity=owner.function.proof,
        positive_tail_majorant=owner.majorant_proof,actual_pressure_tail_views=views,
        source_graph_assertions=owner.assert_graph(),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(GATES+OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(correlated.report(result),separators=(',',':'))+'\n').encode(),mtime=0))
    return owner,views
