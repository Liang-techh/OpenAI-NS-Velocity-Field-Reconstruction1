"""Actual core Banach map, implicit analytic pressure and common radial tails.

Retains exact source functions and nonnegative radial measures symbolically.
Directed boxes enclose their coefficients; no mass, amplitude or point field
is selected. This adapter adds provenance to unchanged current consumers.
"""
import ast
import copy
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_core_nonlinear_operator import (
    CurrentCoreNonlinearOperator,GATE as OPERATOR_GATE,NORM_BOUNDS_GATE,function,binding,sha)
from lei_ren_part1_paper_compliant_current_core_scaled_swirl_source import OPEN
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_logarithmic_pressure_datum import BETA2,BETA0
from lei_ren_part1_paper_analytic_radial_tail import tail_factor

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
NAME=PREFIX+'current_core_common_fixed_point.json'
RECEIPT=PREFIX+'current_core_common_fixed_point_check.json'
BANACH_GATE='current_core_actual_operator_Banach_fixed_point_admitted'
GATE='finite_rows_and_tails_bound_to_same_nonlinear_fixed_point'
TAIL_GATE='current_core_actual_operator_contraction_and_tail_certified'
STILL_OPEN=tuple(key for key in OPEN if key!=GATE)


def raw_flatten_pressure_function():
    """Exact local density and q-kernel of the original 100-unit flatten.

    A function-valued integral is the source. It is not a quadrature fit or
    the previously stored flatten upper bound. Huge inverse-mu exponents
    remain separated from the finite logarithmic terms.
    """
    t,z,mu,yd,Tw=s.symbols('offset Z mu yd Tw',real=True)
    x=t/100
    edge=lambda v:s.exp(-1/v**2)
    sigma=s.Piecewise((0,x<=0),(1,x>=1),(edge(x)/(edge(x)+edge(1-x)),True))
    beta=2-2*sigma
    yv=yd+1+Tw+13/mu
    # sigma(x)+sigma(1-x)=1 gives int_0^1 sigma=1/2 exactly.
    logA=-(s.Rational(1,2)+mu)*(yv+t)+(s.Rational(3,5)+mu)/2+mu*yd
    # Expand only the radial logarithm. Keeping sigma intact preserves its
    # defining functional identity rather than expanding Piecewise branches.
    logdensity=s.expand(2*logA)-2*sigma*s.log(2)-s.log(2)
    finite=s.expand(2*logA+13/mu)-2*sigma*s.log(2)-s.log(2)
    q=1+z*z
    M2,M0=s.symbols('true_M_beta2 true_M_beta0',nonnegative=True)
    flatten=s.Integral(s.exp(logdensity)*q**(-beta),(t,0,100))
    normalized=-(M2/q**2+M0+flatten)
    return dict(offset=t,Z=z,mu=mu,yd=yd,Tw=Tw,sigma=sigma,beta=beta,
                normalized_log_A=logA,log_density=logdensity,finite_log_density=finite,
                inverse_mu_log_density_coefficient=s.Integer(-13),
                normalized_pressure=normalized,flatten_integral=flatten,
                fixed_mass_symbols=(M2,M0))


def all_order_pressure_source_bindings():
    def class_binding(stem,class_name,target,expression):
        tree=ast.parse((HERE/('lei_ren_part1_paper_'+stem+'.py')).read_bytes())
        cls=next(node for node in tree.body if isinstance(node,ast.ClassDef) and node.name==class_name)
        fn=next(node for node in cls.body if isinstance(node,ast.FunctionDef) and node.name=='__init__')
        values=[node.value for node in ast.walk(fn) if isinstance(node,ast.Assign)
                and any(ast.unparse(t)==target for t in node.targets)]
        if len(values)!=1 or ast.dump(values[0])!=ast.dump(ast.parse(expression,mode='eval').body):
            raise ValueError('Original flatten coordinate/parameter definition changed: '+target)
    for target,expression in [('self.y_w','self.y_d+Decimal(1)'),('self.y_p','self.y_w+self.Tw'),
            ('self.y_v','self.y_p+Decimal(13)/self.mu'),('self.Tf','Decimal(100)'),
            ('self.y_f','_decimal_add_exact(self.y_v,self.Tf)'),('self.y_rel','self.y_f-Decimal(30)*self.log_mu')]:
        class_binding('outer','PaperOuterSchedule',target,expression)
    binding('candidate_pressure_function','q_jets','a','1+center**2')
    binding('candidate_pressure_function','q_jets','b','2*center')
    binding('candidate_pressure_function','q_jets','coefficients','[1/a**2]')
    fn=function('candidate_pressure_function','q_jets')
    append=[node for node in ast.walk(fn) if isinstance(node,ast.Call)
            and ast.unparse(node.func)=='coefficients.append']
    expected=ast.parse('-(b*(n+2)*current+(n+3)*previous)/(a*(n+1))',mode='eval').body
    if len(append)!=1 or ast.dump(append[0].args[0])!=ast.dump(expected):
        raise ValueError('All-order q inverse-square recurrence changed')
    # Bind the actual loop and all its parity/Cauchy branches, not only one
    # finite jet or one assignment from a different branch.
    fn=function('logarithmic_pressure_datum','normalized_jets')
    loops=[node for node in ast.walk(fn) if isinstance(node,ast.For)]
    loop=ast.parse("""for n,a in enumerate(q_jets(c,z,order)):
    if n==0:remainder=self.stages['z_flatten']['mass']
    elif n%2 and endpoints(z)==(mp.mpf(0),mp.mpf(0)):
        remainder=c.mpf(0)
    else:
        bound=self.flatten_complex_upper/self.rho**n
        remainder=c.mpf([endpoints(-bound)[0],endpoints(bound)[1]])
    rows.append(-(self.m2*a+(self.m0 if n==0 else 0)+remainder))
""").body[0]
    if len(loops)!=1 or ast.dump(loops[0])!=ast.dump(loop):
        raise ValueError('Original all-order pressure source/parity/Cauchy loop changed')
    # Original raw amplitude and variable-beta dependence.
    fn=function('outer','_log_A')
    returns=[node for node in ast.walk(fn) if isinstance(node,ast.Return)]
    expression="self.logPstar+Decimal('0.1')*y-Decimal('0.6')*self._J(y)-self.mu*self._J(y-self.y_d)-(Decimal(1)-self.mu)*self._J(y-self.y_rel)+(Decimal(1)-self.delta/Decimal(2))*self._J(y-self.y_rel-Decimal(1)-self.Ts)"
    if len(returns)!=1 or ast.dump(returns[0].value)!=ast.dump(ast.parse(expression,mode='eval').body):
        raise ValueError('Original radial log-amplitude changed')
    binding('outer','_evaluate_log_radius','flat_sigma','_sigma((y-self.y_v)/self.Tf)')
    binding('outer','_evaluate_log_radius','log_u',
            'log_a-Decimal(str(z_factor_log))+Decimal(str(flat_sigma*(z_factor_log-math.log(2.0))))',True)
    tree=ast.parse((HERE/(PREFIX+'pressure_source.py')).read_bytes())
    classes=[node for node in tree.body if isinstance(node,ast.ClassDef) and node.name=='CompliantPressureDatum']
    if len(classes)!=1:raise ValueError('Unique current pressure source class required')
    initializers=[node for node in classes[0].body if isinstance(node,ast.FunctionDef) and node.name=='__init__']
    if len(initializers)!=1:raise ValueError('Unique current pressure source initializer required')
    for target,expression in [('self.flatten_complex_upper','self.tail_upper/(1-self.rho**2)**2'),
                              ('self.rho',"c.mpf('.25')")]:
        values=[node.value for node in ast.walk(initializers[0]) if isinstance(node,ast.Assign)
                and any(ast.unparse(t)==target for t in node.targets)]
        if len(values)!=1 or ast.dump(values[0])!=ast.dump(ast.parse(expression,mode='eval').body):
            raise ValueError('Current pressure source all-order bound changed: '+target)
    kernel=raw_flatten_pressure_function()
    return dict(original_q_jets_all_order_recurrence_AST_bound=True,
                original_normalized_pressure_complete_loop_AST_bound=True,
                original_radial_amplitude_and_variable_beta_AST_bound=True,
                original_flatten_coordinates_and_length_AST_bound=True,
                exact_raw_flatten_function={key:str(kernel[key]) for key in
                    ('sigma','beta','log_density','finite_log_density','normalized_pressure')},
                fixed_true_masses_not_replaced_by_interval_values=True,
                flatten_true_integral_not_replaced_by_Cauchy_cap=True,
                raw_policy='original implicit preheat ansatz with H replaced by1; fixed source parameters and radial measures',
                infinite_radial_mass_integrals_remain_implicit_not_point_selected=True,passed=True)


@source_precision
def replay_production_majorants(core):
    """Replay the actual term function/calls at original producer precision."""
    major=core.records['core_transfer']
    c=MPIntervalContext();c.dps=major['precision']
    read=lambda key:read_interval(c,major[key])
    env=dict(ctx=c,eps=read('epsilon'),dt=c.mpf('1e-200'),
             product=read('product_constant'),J1=read('J1_norm_upper'),J2=read('J2_norm_upper'),
             axial=read('integrated_axial_derivative_constant'),radial=read('integrated_radial_derivative_constant'),
             fixed_multiplier=read('fixed_analytic_multiplier_factor_upper'),
             Bphi=read('Phi_ball_norm_upper'),Bpsi=read('Psi_ball_norm_upper'),
             Pcal=read('restored_pressure_coefficient_norm_upper'),Fsquare=read('F0_squared_Xh_norm_upper'),
             coefficients={key:read_interval(c,row) for key,row in major['fixed_coefficient_norm_upper'].items()},rows=[])
    fn=function('compliant_core_transfer','run')
    definitions=[node for node in ast.walk(fn) if isinstance(node,ast.FunctionDef) and node.name=='term']
    calls=[node for node in ast.walk(fn) if isinstance(node,ast.Expr) and isinstance(node.value,ast.Call)
           and isinstance(node.value.func,ast.Name) and node.value.func.id=='term']
    if len(definitions)!=1 or len(calls)!=20:
        raise ValueError('Complete actual term function and twenty calls required')
    tree=ast.Module(body=copy.deepcopy(definitions+calls),type_ignores=[])
    exec(compile(ast.fix_missing_locations(tree),'<actual-twenty-majorants>','exec'),env)
    rows=env['rows']
    for actual,stored in zip(rows,major['terms']):
        for key in ('component','term','phi_power','psi_power'):
            if actual[key]!=stored[key]:raise ValueError('Changed actual majorant row identity')
        for key in ('coefficient_upper','size_upper','Lipschitz_upper'):
            if actual[key]._mpi_!=read_interval(c,stored[key])._mpi_:
                raise ValueError('Actual term arithmetic does not replay stored majorant: '+actual['term']+' '+key)
    sums={}
    for component in ('theta','z'):
        for key,tag in (('size_upper','size'),('Lipschitz_upper','lip')):
            sums[component+'_'+tag]=sum((row[key] for row in rows if row['component']==component),c.mpf(0))
    R=read('angular_resolvent_norm_upper')
    size=(R*sums['theta_size']+sums['z_size'])/2
    lip=(R*sums['theta_lip']+sums['z_lip'])/2
    actual_m,actual_q=env['eps']*size,env['eps']*lip
    for key,value in [('nonlinear_map_size_upper',size),('nonlinear_map_Lipschitz_upper',lip),
                      ('scaled_map_size_upper',actual_m),('scaled_map_Lipschitz_upper',actual_q)]:
        if value._mpi_!=read(key)._mpi_:raise ValueError('Combined actual map bound changed: '+key)
    if endpoints(actual_m)[1]>mp.mpf('.5') or endpoints(actual_q)[1]>=mp.mpf('.5'):
        raise ValueError('Actual closed-ball/contraction gate fails')
    correction=actual_m/(1-actual_q)
    if endpoints(correction)[1]>1:
        raise ValueError('Conservative correction no longer inside unit ball')
    return dict(precision=c.dps,actual_twenty_rows=rows,angular_resolvent_norm_upper=R,
                size_bound=actual_m,Lipschitz_bound=actual_q,closed_ball_radius=c.mpf(1),
                correction_bound=correction,strict_contraction_margin=1-actual_q,
                all_twenty_rows_and_combined_bounds_bit_exact_replayed=True,
                actual_m_le_one_half_and_actual_q_lt_one_half=True,passed=True)


@source_precision
def same_source_models_and_domains(operator,replay):
    source=operator.source;core=source.core;c=core.ctx
    major=core.records['core_transfer'];tube=core.records['shared_analytic_tube'];linear=core.records['shared_linear_resolvent']
    read=lambda row,key:read_interval(c,row[key])
    same=[(core.epsilon,read(major,'epsilon')),(core.Lambda,read(major,'Lambda')),
          (core.logLambda,read(major,'logLambda')),(core.logP,read(major,'logPstar')),
          (core.j,read(tube,'j')),(core.h,read(tube,'Xh_parameter')),
          (source.radius,read(linear,'Cauchy_disk_radius'))]
    if any(a._mpi_!=b._mpi_ for a,b in same):
        raise ValueError('Actual source parameters differ from shared model/tube definitions')
    # Sigma is recomputed from the same j in the consumer's higher-precision
    # context. Bind the defining formulas and require directed containment
    # in the producer bound; differently rounded intervals are not identities.
    binding('compliant_core_physical_field','__init__','self.sigma','self.j/500')
    binding('shared_analytic_tube','run','sigma','j/500')
    sigma_lo,sigma_hi=endpoints(core.sigma);tube_lo,tube_hi=endpoints(read(tube,'sigma'))
    if core.sigma._mpi_!=(core.j/500)._mpi_ or sigma_lo<tube_lo or sigma_hi>tube_hi:
        raise ValueError('Same defining sigma is outside the declared tube bound')
    if (endpoints(core.delta)[0]<0 or endpoints(core.delta)[1]>endpoints(c.mpf('1e-200'))[1]
        or not tube['common_complex_axis_poles_excluded'] or source.radius._mpi_!=(read(tube,'complex_tube_radius')/2)._mpi_
        or core.datum.source_sha!=source.source or core.datum.datum_sha!=source.datum_sha
        or len(core.datum.stages)!=14):
        raise ValueError('Actual continuous pressure/tube/source scope differs')
    for name,row in core.datum.stages.items():
        if endpoints(row['mass'])[0]<0:raise ValueError('Negative radial pressure measure enclosure: '+name)
    # Current exact source lies under the baseline majorant by the selected
    # Cstar surplus. Its weighted physical F0² bound stays in logarithms.
    physical_log_bound=source.required_log+2*core.logLambda
    baseline_log=c.ln(read(major,'F0_squared_Xh_norm_upper'))
    if endpoints(baseline_log-physical_log_bound)[0]<=0:
        raise ValueError('Actual selected amplitude is outside the nonlinear-map norm bound')
    # Complex pressure bounds use the current exact radial measures. The
    # old/new envelope is a numerical upper bound, never a definition of P0.
    eta=read(tube,'complex_tube_radius');qmin=read(tube,'pressure_q_modulus_lower')
    datum=core.datum
    p=datum.parameters
    if (endpoints(p.mu)[0]<=0 or endpoints(p.mu)[1]>=1
        or endpoints(p.Tw)[0]<0 or endpoints(p.Ts)[0]<=0 or endpoints(p.yd)[0]<1):
        raise ValueError('Original flatten radial amplitude branch premises fail')
    source_premises=dict(zero_lt_mu_lt_one=True,Tw_nonnegative=True,Ts_positive=True,yd_ge_one=True,
        flatten_log_branch_proof=[
            'y_v=y_d+1+Tw+13/mu; for 0<=offset<=100 both y and y-y_d exceed1.',
            'y_rel=y_v+100-30log(mu), so y-y_rel=offset-100+30log(mu)<0.',
            'The final J argument y-y_rel-1-Ts is also negative.',
            'Thus J(y)=y-1+J(1), J(y-y_d)=y-y_d-1+J(1), and both late J terms vanish.'])
    pressure_bound=c.exp(2*core.logP)*((c.mpf(datum.m2)+c.mpf(datum.stages['z_flatten']['mass']))/qmin**2+c.mpf(datum.m0))
    derivative_bound=4*(1+eta)*pressure_bound/qmin
    if (endpoints(pressure_bound)[1]>endpoints(read(major,'fresh_pressure_norm_upper'))[1]
        or endpoints(derivative_bound)[1]>endpoints(read(major,'fresh_pressure_derivative_upper'))[1]):
        raise ValueError('Current all-order pressure source exceeds the admitted model envelope')
    # AST bindings to actual leading models and the correction consumed by
    # both finite/tail algorithms. No equality is inferred from box overlap.
    for target,expr in [('u','4*z+self.j'),('d','one-z*z'),('L','one-(z*z)*self.delta'),
        ('H','z*((1-self.delta)/2)+d*u'),('denominator','h2+self.sigma**2'),
        ('physicalP','pressure*c.exp(2*self.logP)'),('pz','IntervalTaylor(c,[(k+1)*physicalP[k+1] for k in range(6)])'),
        ('g','-((1-2*z*u)*u*((1+self.delta)/2)+4*H+d*pz-z*physicalP*(2*(1+self.delta)))'),
        ('slope','-g/(2*L)')]:
        binding('compliant_core_physical_field','axis_inputs',target,expr)
    binding('compliant_core_physical_field','normalized_jets','model',
            "(source['slope'][k]*math.factorial(k))*(rho if i==0 else 1 if i==1 else 0)")
    binding('compliant_core_physical_field','__init__','self.correction',"read(major,'scaled_map_size_upper')/(1-q)")
    # Replay at the consumer context as well, to avoid comparing differently
    # rounded producer/consumer intervals by overlap.
    consumer_correction=read(major,'scaled_map_size_upper')/(1-read(major,'scaled_map_Lipschitz_upper'))
    if consumer_correction._mpi_!=core.correction._mpi_:
        raise ValueError('Actual physical and finite-tail consumers do not use the newly bound map correction')
    return dict(same_current_original_source_graph=source.current_source_graph(),
                actual_parameter_model_and_tube_bindings=True,current_fourteen_nonnegative_measure_enclosures=True,
                same_sigma_definition_and_directed_producer_bound_containment=True,
                actual_flatten_parameter_branch_premises=source_premises,
                actual_pressure_complex_modulus_upper=pressure_bound,actual_pressure_complex_derivative_upper=derivative_bound,
                selected_physical_F0_squared_Xh_log_upper=physical_log_bound,
                selected_source_log_majorant_margin=baseline_log-physical_log_bound,
                leading_Psi_is_B_rho_over_two=True,
                leading_Phi_is_same_chi_Bessel_model=True,
                consumer_correction_rederived_from_actual_map=consumer_correction,
                old_new_envelope_is_only_an_upper_bound_not_a_source_definition=True,
                core_and_pressure_source_epsilons_distinct=True,passed=True)


def consumer_source_bindings():
    for target,expr in [('average','psi.radial_average()'),('pressure','physical_S*(phi*phi).radial_integral()'),
                       ('physical_S','scaled_S/(eps*eps)')]:
        binding('compliant_current_core_nonlinear_operator','source_terms',target,expr)
    for target,expr in [('theta',"sum(packet['terms']['theta'],packet['phi'].constant(0))"),
                       ('axial',"sum(packet['terms']['z'],packet['psi'].constant(0))"),
                       ('chi_phi',"packet['chi']*packet['phi']"),
                       ('drive',"packet['B'] if n==0 else packet['B'].constant(0)"),
                       ('squared',"packet['scaled_S']*(packet['phi']*packet['phi'])")]:
        binding('compliant_current_core_nonlinear_operator','next_rows',target,expr)
    returns=[node for node in ast.walk(function('compliant_current_core_nonlinear_operator','next_rows')) if isinstance(node,ast.Return)]
    expected="dict(A=[(eps*v-c)/(2*(n+1)*(n+2)) for v,c in zip(theta.rows[n],chi_phi.rows[n])],Uz=[eps*(b+eps*v)/(2*(n+1)**2) for b,v in zip(drive.rows[n],axial.rows[n])],P=[v/(n+1) for v in squared.rows[n]])"
    if len(returns)!=1 or ast.dump(returns[0].value)!=ast.dump(ast.parse(expected,mode='eval').body):
        raise ValueError('Actual arbitrary-order radial extraction changed')
    for target,expr in [('first','(rmax**(n-i)/(2**n*math.factorial(n-i)*math.factorial(n+1))*(k+1)*(n+1)**k*(1+M)**k)'),
                       ('ratio','rmax/2*(c.mpf(n+2)/(n+1))**k/((n+1-i)*(n+2))'),
                       ('tail','first/(1-ratio)*math.factorial(k)')]:
        binding('compliant_core_physical_field','model_phi_jets',target,expr)
    fn=function('compliant_core_coefficient_rebuild','rebuild')
    calls=[node for node in ast.walk(fn) if isinstance(node,ast.Call) and ast.unparse(node.func)=='advance_scaled_one']
    if len(calls)!=1 or ast.unparse(calls[0])!='advance_scaled_one(c, fixed, rows, n, c.mpf(Z), self.core.delta, self.core.epsilon)':
        raise ValueError('Finite production recurrence consumer changed')
    binding('compliant_core_coefficient_rebuild','profile','factor',
        "tail_factor(c,degree=N,radial_order=i,axial_order=k,radius=c.mpf(endpoints(r)[1]),h=self.core.h)['tail_per_Xh_norm']")
    for expr in ("self.model_tail_cache[key][gridkey(i,k)]+self.core.correction*factor",'self.core.epsilon*self.core.correction*factor'):
        binding('compliant_core_coefficient_rebuild','profile','tail',expr,True)
    binding('compliant_rooted_core_field','rooted_profile','tail','self.core.correction*factor')
    binding('compliant_core_coefficient_rebuild','values','ptail',
        'self.S_bound*self.phi_norm**2*r**(N+1)/(20**N*(1-r/20))')
    return dict(actual_all_order_map_and_integrals_AST_bound=True,
                original_factorial_model_tail_consumer_AST_bound=True,
                original_finite_recurrence_callable_AST_bound=True,
                original_model_plus_nonlinear_norm_tail_consumers_AST_bound=True,
                rooted_nonlinear_tail_uses_same_actual_map_correction=True,
                same_pressure_primitive_coefficient_tail_consumer_AST_bound=True,passed=True)


@source_precision
def true_pressure_measure_bindings(datum):
    """Bind current true measures, using the old receipt only on its unchanged early atom."""
    c=datum.ctx;p=datum.parameters
    expected2=('reference_extension','slope_transition_ref','axial_turnoff','slope_transition_mu','power_buffer','pulse_reserved')
    expected0=('power_buffer_rel','steep_transition_in','steep_power','steep_transition_out','waiting','heat_collar','exterior_power_tail')
    if BETA2!=expected2 or BETA0!=expected0 or set(datum.stages)!=set(expected2+expected0+('z_flatten',)):
        raise ValueError('True radial measure partition differs')
    old_name='lei_ren_part1_paper_coherent_uniform_fixed_beta_error.json'
    old=json.loads((HERE/old_name).read_bytes());_verify_hashes(old)
    refined_name='lei_ren_part1_paper_candidate_pressure_mass_refinement.json'
    refined=json.loads((HERE/refined_name).read_bytes());_verify_hashes(refined)
    if (not old['accepted_source_alignment_certified'] or not refined['pressure_source_unchanged']
        or old['accepted_schedule_sha256']!=refined['accepted_schedule_sha256']):
        raise ValueError('Independent early true mass alignment missing')
    tail=c.exp(c.mpf('.6')-p.yd)/(2*(1-p.epsilon)**2)
    early=read_interval(c,refined['stages']['slope_transition_ref']['refined_mass'])
    turnoff=c.exp(c.mpf('.6'))*(c.exp(-1)-c.exp(-p.yd))/2
    exact_masses={'reference_extension':c.mpf('2.5'),'slope_transition_ref':early,'axial_turnoff':turnoff}
    for name,row in datum.stages.items():
        wanted=exact_masses.get(name,c.mpf([0,endpoints(tail)[1]]))
        if row['mass']._mpi_!=wanted._mpi_:raise ValueError('Current true mass enclosure formula differs: '+name)
    if (datum.m2._mpi_!=sum((datum.stages[n]['mass'] for n in BETA2),c.mpf(0))._mpi_
        or datum.m0._mpi_!=sum((datum.stages[n]['mass'] for n in BETA0),c.mpf(0))._mpi_):
        raise ValueError('Common fixed-beta sums omit a true measure')
    if (endpoints(p.epsilon)[0]<=0 or endpoints(p.epsilon)[1]>=1 or endpoints(p.delta)[0]<=0):
        raise ValueError('Positive terminal measure parameters missing')
    return dict(all_fourteen_current_true_pressure_measures_bound=True,
        true_M2_defined_by_stages=list(BETA2),true_M0_defined_by_stages=list(BETA0),
        flatten_measure_is_raw_functional_integral=True,
        accepted_old_alignment_used_only_for_parameter_independent_early_atom=True,
        early_exact_density='exp(y/5-6J(y)/5)/2 on0<=y<=1; J is the same exact cutoff primitive',
        early_measure_receipt=refined_name,
        new_source_tail_proof=[
            'For y<0 the reference density is exp(y/5)/2; its integral equals5/2.',
            'On0<=y<=1 the normalized source is exp(y/10-3J(y)/5), independent of Md,mu,delta,epsilon and waiting.',
            'On1<=y<=yd it is exp(3/10-y/2), since exact J(1)=1/2; direct integration gives the turnoff mass.',
            'For y>=yd, split the remaining log amplitude into -mu*(J(y-yd)-Jlate) -(Jrel-Jlate) -delta*Jlate/2.',
            'Every bracket is nonnegative because J is nonnegative and increasing, yd<yrel<yrel+1+Ts.',
            'Thus A/Pstar<=exp(3/10-y/2); flatten adds the factor2^-sigma<=1.',
            'With H replaced by1, the terminal heat factor K lies in[1-epsilon,1]; its matching normalization costs at most1/(1-epsilon).',
            'Each remaining disjoint stage mass is at most integral_yd^infinity exp(3/5-y)/(2*(1-epsilon)^2) dy.',
            'This equals the current positive tail_upper. The old late source and old epsilon are never transferred.'],passed=True)


class CurrentCoreCommonFixedPoint:
    @source_precision
    def __init__(self,operator=None,require_checked=True):
        self.operator=operator if operator is not None else CurrentCoreNonlinearOperator()
        if not self.operator.acceptance_loaded:raise ValueError('Checked actual twenty-term operator required')
        self.source=self.operator.source;self.core=self.source.core;self.rebuild=self.source.rebuild;self.ctx=self.core.ctx
        self.family=self.source.family;self.source_sha=self.source.source;self.datum_sha=self.source.datum_sha
        self.hashes=dict(self.operator.hashes)
        parent_name=PREFIX+'current_core_nonlinear_operator_check.json'
        parent=accepted(parent_name,self.family,self.source_sha,OPERATOR_GATE);_verify_hashes(parent)
        if not parent[NORM_BOUNDS_GATE]:raise ValueError('Actual operator norm incidence proof required')
        self.parent=parent
        steep_name=PREFIX+'current_steep_waiting_source_check.json'
        steep=accepted(steep_name,self.family,self.source_sha,'current_steep_waiting_source_ownership_certified')
        _verify_hashes(steep)
        steep_raw_name=PREFIX+'current_steep_waiting_source.json'
        steep_raw=json.loads((HERE/steep_raw_name).read_bytes())
        kernels=steep_raw['current_steep_waiting_source_bindings']['exact_current_kernel_and_parameter_source_bindings']
        reflection={key:kernels['actual_shared_kernel_functional_identities'][key] for key in
            ('original_sigma_reflection_identity','original_J1_half_by_reflected_exact_integral')}
        if not all(reflection.values()) or not kernels['same_exact_original_logistic_sigma_function']:
            raise ValueError('Exact reflected integral J(1)=1/2 source proof missing')
        self.exact_reflection=reflection
        self.hashes[steep_name]=sha(steep_name);self.hashes[steep_raw_name]=sha(steep_raw_name)
        self.hashes[parent_name]=sha(parent_name)
        for name in [Path(__file__).name,'lei_ren_part1_paper_candidate_pressure_function.py',
                     'lei_ren_part1_paper_logarithmic_pressure_datum.py','lei_ren_part1_paper_outer.py',
                     'lei_ren_part1_paper_analytic_radial_tail.py','lei_ren_part1_paper_shared_commuting_resolvent.py',
                     'lei_ren_part1_paper_shared_commuting_resolvent.json',
                     '../../src/openai_ns_reconstruction/schedule_pressure.py']:
            self.hashes[name]=sha(name)
        self.pressure_bindings=all_order_pressure_source_bindings()
        self.replay=replay_production_majorants(self.core)
        self.model_bindings=same_source_models_and_domains(self.operator,self.replay)
        self.consumer_bindings=consumer_source_bindings()
        self.measure_bindings=true_pressure_measure_bindings(self.core.datum)
        for name in ('lei_ren_part1_paper_coherent_uniform_fixed_beta_error.json',
                     'lei_ren_part1_paper_candidate_pressure_mass_refinement.json'):
            self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source_sha,GATE);_verify_hashes(receipt)
            if (receipt['datum_enclosure_sha256']!=self.datum_sha or not receipt[BANACH_GATE]
                or not receipt[TAIL_GATE] or any(receipt[key] for key in STILL_OPEN)):
                raise ValueError('Common core fixed-point receipt source/scope differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT)
            self.acceptance_loaded=True

    @source_precision
    def evaluate(self,Z,rho='4.1',degree=6,depth=5):
        if degree<6 or depth<1:raise ValueError('Degree>=6 and depth>=1 required for decomposed model/nonlinear tail')
        packet=self.rebuild.rebuild(Z,degree=degree,depth=depth)
        selected=[]
        for i,k in ((0,0),(1,0),(0,1),(2,0),(1,1)):
            selected.append(self.rebuild.profile(packet,rho,radial_order=i,axial_order=k))
        return dict(Z=self.ctx.mpf(Z),rho=self.ctx.mpf(rho),current_original_finite_packet=packet,
                    common_solution_profile_packets=selected,
                    common_solution_pressure_profile_packets=[self.pressure_profile(packet,rho,i,k)
                        for i,k in ((0,0),(1,0),(0,1),(2,0),(1,1)) if k<=depth],
                    same_production_finite_coefficients_and_model_nonlinear_tails=True,
                    source_pressure_and_swirl_not_point_selected=True,
                    **dict.fromkeys((BANACH_GATE,GATE,TAIL_GATE),self.acceptance_loaded),
                    **dict.fromkeys(STILL_OPEN,False))

    @source_precision
    def pressure_profile(self,packet,rho,radial_order=0,axial_order=0):
        """Same pressure primitive, with both complete axial convolutions in its tail."""
        c=self.ctx;r=c.mpf(rho);i=radial_order;k=axial_order;N=packet['radial_degree']
        if (not isinstance(i,int) or not isinstance(k,int) or min(i,k)<0 or i>N or k>packet['axial_depth']
            or endpoints(r)[0]<0 or endpoints(r)[1]>endpoints(c.mpf('4.1'))[1]):
            raise ValueError('Pressure derivative outside original core/packet domain')
        polynomial=sum((packet['rows']['P'][n][k]*(math.factorial(n)//math.factorial(n-i))*math.factorial(k)*r**(n-i)
                        for n in range(i,N+1)),c.mpf(0))
        major=self.core.records['core_transfer']
        fixed=read_interval(c,major['fixed_analytic_multiplier_factor_upper'])
        norm=80*256*fixed*self.rebuild.S_bound*self.rebuild.phi_norm**2
        # n>N>=i makes every omitted rho^(n-i) term vanish on the axis.
        # The generic geometric helper intentionally requires positive radius.
        factor=(c.mpf(0) if endpoints(r)[1]==0 else
            tail_factor(c,degree=N,radial_order=i,axial_order=k,radius=c.mpf(endpoints(r)[1]),h=self.core.h)['tail_per_Xh_norm'])
        tail=norm*factor;upper=endpoints(tail)[1]
        return dict(radial_order=i,axial_order=k,finite_polynomial=polynomial,
                    mixed_pressure_increment_Xh_norm_upper=norm,infinite_tail_upper=tail,
                    pressure_scaled_enclosure=polynomial+c.mpf([-upper,upper]),
                    pressure_units='epsilon_core*physical_P',
                    same_S_Phi_product_primitive=True,source_values_not_selected=True)


@source_precision
def run(field=None):
    field=field if field is not None else CurrentCoreCommonFixedPoint(require_checked=False)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source_sha,
                datum_enclosure_sha256=field.datum_sha,all_order_pressure_source_bindings=field.pressure_bindings,
                accepted_exact_sigma_reflection_and_J1=field.exact_reflection,
                actual_operator_majorant_replay=field.replay,same_actual_model_source_domains=field.model_bindings,
                finite_and_tail_consumer_bindings=field.consumer_bindings,
                current_true_pressure_measure_bindings=field.measure_bindings,
                fresh_common_core_runtime=field.evaluate('.327'),
                mathematical_scope='Unique analytic fixed point for each fixed original implicit source within admitted directed data; coefficient/tail enclosures, not selected nonlinear point values',
                **dict.fromkeys((BANACH_GATE,GATE,TAIL_GATE),False),**dict.fromkeys(STILL_OPEN,False),
                input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(encode(result),indent=2)+'\n').encode('utf8'))
    print('Built actual map replay and common finite/tail adapter; separate source proof checker required',flush=True)
    return result


if __name__=='__main__':run()
