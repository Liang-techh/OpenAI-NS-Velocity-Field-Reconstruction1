"""Original small angular branch and complete future swirl energy C1.

All fixed weights are integrals and all RHSs are genuine same-graph
functions. The two angular coefficients are the admitted small quadratic
root, not centers of Newton boxes. Complete future energy includes the
original signed bumps, epsilon atoms and nonzero full Gamma deficit.
"""
import copy
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_Gamma_future_functions as Gamma
import lei_ren_part1_paper_compliant_outer_angular_repair as angular_source
import lei_ren_part1_paper_compliant_future_swirl_energy as energy_source
import lei_ren_part1_paper_compliant_current_heat_source as heat_bridge
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

parent, source, controls = Gamma.parent, Gamma.source, Gamma.controls
HERE, PREFIX, sha = Gamma.HERE, Gamma.PREFIX, Gamma.sha
NAME = PREFIX + 'current_original_n1_complete_future_functions.json.gz'
RECEIPT = PREFIX + 'current_original_n1_complete_future_functions_check.json'
GATE = 'current_original_n1_complete_future_energy_C1_function_installed'
ANGULAR = 'current_original_n1_actual_small_angular_C1_functions_installed'
PULSE = 'current_original_n1_actual_selected_pulse_coefficient_C1_functions_installed'
OPEN = Gamma.pulse.OPEN


def sqrt_pair(g, pair, certificate):
    if type(pair) is not source.C1Function:
        raise TypeError('A genuine C1 radicand required')
    g.ids([pair.value, pair.Z])
    root = g.unary('sqrt', pair.value)
    return source.C1Function(root, g.quotient(pair.Z, g.mul(g.constant(2), root), certificate))


def angular_weights(g, mu, normalization):
    """Exact A/B/D/E/F in the original normalized beta measure."""
    g.ids([mu, normalization])
    c = g.constant
    ell, t = c('3/20'), g.symbol('original_future_angular_beta_t')
    beta = g.node('compact_raw_beta', argument=t.node,
        definition='exp(-1/(1-t^2)) on abs(t)<1, zero otherwise',
        defining_module=PREFIX+'outer_pulse_map.py',
        defining_module_sha256=sha(PREFIX+'outer_pulse_map.py'))
    r, p = g.sub(g.one, mu), g.add(g.one, g.mul(c(2), mu))
    out = {}
    for name, rate, power in (('A',r,1),('B',g.neg(p),1),('D',g.neg(p),2),
                             ('E',g.mul(c(-2),mu),1),('F',g.mul(c(-2),mu),2)):
        numerator = controls.integral(g, g.mul(g.unary('exp', g.mul(rate, ell, t)),
            *([beta]*power)), g.nodes[t.node]['name'], c(-1), g.one,
            original_normalized_angular_bump_weight=name, endpoint_independent_of_Z=True)
        denominator = normalization if power == 1 else g.mul(ell, normalization, normalization)
        out[name] = g.quotient(numerator, denominator, 'positive original beta normalization and ell')
    return out


def small_quadratic_functions(g, b1, b2, p, q, kappa, scale):
    """p*x+y=b1, x+q*y+kappa*scale*(x^2+q*y^2)=b2."""
    if any(type(pair) is not source.C1Function for pair in (b1,b2)):
        raise TypeError('Two genuine angular C1 right sides required')
    g.ids([p,q,kappa,scale,b1.value,b1.Z,b2.value,b2.Z])
    a = Gamma.C1Algebra(g)
    fixed = a.fixed
    ks = g.mul(kappa,scale)
    A = g.mul(ks,g.add(g.one,g.mul(q,p,p)))
    B = a.sub(fixed(g.sub(g.one,g.mul(p,q))),a.scale(g.mul(g.constant(2),ks,q,p),b1))
    C = a.sub(a.add(a.scale(q,b1),a.scale(g.mul(ks,q),a.mul(b1,b1))),b2)
    disc = a.sub(a.mul(B,B),a.scale(g.mul(g.constant(4),A),C))
    root = sqrt_pair(g,disc,'positive admitted small-branch discriminant square root')
    # This exact rationalization never divides by the microscopic A.
    x = a.div(a.scale(-2,C),a.add(B,root),'positive B+sqrt(discriminant) on admitted small branch')
    y = a.sub(b1,a.scale(p,x))
    return dict(x=x,y=y,d1=a.scale(scale,x),d2=a.scale(scale,y)),dict(
        A=A,B=B,C=C,discriminant=disc,b1=b1,b2=b2,p=p,q=q,kappa=kappa,scale=scale)


def angular_functions(g, schedule, histories, heat_functions, weights):
    a, c = Gamma.C1Algebra(g), g.constant
    mu, delta, r, k = (schedule[key] for key in ('mu','delta','r','k'))
    logscale = g.mul(c(30),r,g.unary('log',mu))
    scale = g.unary('exp',logscale)
    z = g.symbol('Z')
    Xf0 = g.node('function_substitution', expression=histories['Xf'].value.node,
        variable=z.node,value=g.zero.node,Z_independent_substitution=True)
    pre = a.sub(a.fixed(Xf0),histories['Xf'])
    logtail_profile = g.add(g.mul(c('-3/2'),g.add(c(2),schedule['Ts'])),
        g.mul(c('1/2'),g.add(r,k)),
        g.neg(g.mul(g.add(c('1/2'),g.mul(c('1/2'),delta)),schedule['waiting'])))
    logtheta = g.add(g.mul(c('1/2'),g.add(r,k)),g.mul(k,schedule['waiting']),
        g.neg(schedule['log_one_minus_epsilon']))
    logpressure = g.add(g.mul(c(2),logtail_profile),
        g.mul(c(-2),schedule['log_one_minus_epsilon']),g.unary('log',g.mul(c('1/2'),delta)))
    heat_factor = lambda logmult: g.unary('exp',g.sub(g.sub(logmult,schedule['logRtail']),logscale))
    r_scaled = a.add(pre,a.scale(heat_factor(logtheta),heat_functions['theta']))
    s_scaled = a.scale(heat_factor(logpressure),heat_functions['pressure'])
    b1 = a.scale(g.quotient(g.unary('exp',r),weights['A'],'positive original angular A'),r_scaled)
    prate = g.add(g.one,g.mul(c(2),mu))
    b2 = a.scale(g.quotient(g.unary('exp',g.mul(c(-3),prate)),weights['B'],
        'positive original angular B'),s_scaled)
    p, q = g.unary('exp',g.mul(c(-2),r)),g.unary('exp',g.mul(c(-2),prate))
    kappa = g.quotient(weights['D'],g.mul(c(2),weights['B']),'positive 2B')
    selected, quadratic = small_quadratic_functions(g,b1,b2,p,q,kappa,scale)
    return selected,dict(quadratic=quadratic,r_scaled=r_scaled,s_scaled=s_scaled,
        preheat_difference=pre,logscale=logscale,log_theta_multiplier=logtheta,
        log_pressure_multiplier=logpressure,log_tail_profile_relative=logtail_profile)


def complete_energy_functions(g, schedule, heat_functions, angular, weights):
    a, c = Gamma.C1Algebra(g), g.constant
    fixed = a.fixed
    mu, delta, eps = (schedule[key] for key in ('mu','delta','epsilon'))
    z = source.C1Function(g.symbol('Z'),g.one)
    q = a.add(fixed(1),a.mul(z,z))
    t, v = g.symbol('original_future_energy_t'),g.symbol('original_future_energy_primitive_v')
    logratio = a.log(a.scale('1/2',q),'positive q/2')
    flatten = a.integral(a.scale(g.unary('exp',g.mul(c(-2),mu,t)),
        a.exp(a.scale(g.mul(c(2),Gamma.sigma(g,g.mul(c('1/100'),t))),logratio))),
        t,g.zero,c(100),original_complete_future_piece='flatten')
    Nf = a.scale(g.mul(c('1/4'),g.unary('exp',g.mul(c(-200),mu))),a.mul(q,q))
    Nrel = a.scale(g.mul(c('1/4'),g.unary('exp',g.mul(c(-2),mu,g.add(c(100),schedule['Lrel'])))),a.mul(q,q))
    decay = lambda rate,length: g.mul(length,g.unary('exprel',g.neg(g.mul(rate,length))))
    power = a.scale(decay(g.mul(c(2),mu),schedule['Lrel']),Nf)
    change = a.add(*(a.scale(g.unary('exp',g.mul(c(-2),mu,c(center))),
        a.add(a.scale(g.mul(c(2),weights['E']),dj),a.scale(weights['F'],a.mul(dj,dj))))
        for dj,center in zip((angular['d1'],angular['d2']),(-3,-1))))
    change = a.mul(Nrel,change)
    J = controls.integral(g,Gamma.sigma(g,v),g.nodes[v.node]['name'],g.zero,t,
        original_unit_steep_primitive=True,endpoint_independent_of_Z=True)
    Iin = controls.integral(g,g.unary('exp',g.add(g.mul(c(-2),mu,t),
        g.mul(c(-2),schedule['r'],J))),g.nodes[t.node]['name'],g.zero,g.one,
        original_complete_future_piece='steep_in')
    Iout = controls.integral(g,g.unary('exp',g.add(g.mul(c(-2),t),
        g.mul(c(2),schedule['k'],J))),g.nodes[t.node]['name'],g.zero,g.one,
        original_complete_future_piece='steep_out')
    Ns = g.unary('exp',g.sub(c(-1),mu))
    Nq = g.mul(Ns,g.unary('exp',g.mul(c(-2),schedule['Ts'])))
    Nt = g.mul(Nq,g.unary('exp',g.sub(c(-1),g.mul(c('1/2'),delta))))
    Ntail = g.mul(Nt,g.unary('exp',g.neg(g.mul(delta,schedule['waiting']))))
    tail_multiplier = g.mul(Ntail,g.unary('exp',g.mul(c(-2),schedule['log_one_minus_epsilon'])))
    atoms = schedule['preheat_atoms']['energy']
    baseline = g.quotient(g.one,delta,'strictly positive original delta')
    eps_atom = g.mul(c(-2),eps,atoms['W'])
    eps2_atom = g.mul(eps,eps,atoms['W_squared'])
    preheat = g.add(baseline,eps_atom,eps2_atom)
    postrel = a.scale(g.add(Iin,g.mul(Ns,decay(c(2),schedule['Ts'])),g.mul(Nq,Iout),
        g.mul(Nt,decay(delta,schedule['waiting'])),g.mul(tail_multiplier,preheat)),Nrel)
    heat = a.mul(Nrel,a.scale(g.mul(tail_multiplier,g.mul(c('1/2'),delta),schedule['S']),heat_functions['energy']))
    total = a.sub(a.add(flatten,power,change,postrel),heat)
    weighted = a.scale(g.mul(c('1/2'),mu,g.unary('exp',c(-26))),total)
    return dict(complete=total,weighted=weighted),dict(flatten=flatten,power=power,
        angular_change=change,steep_waiting_preheat_tail=postrel,Gamma_heat_difference=heat,
        Nf=Nf,Nrel=Nrel,Ns=Ns,Nq=Nq,Nt=Nt,Ntail=Ntail,tail_multiplier=tail_multiplier,
        preheat_baseline=baseline,epsilon_atom=eps_atom,epsilon_squared_atom=eps2_atom,
        steep_in_integral=Iin,steep_out_integral=Iout)


def checked_sources(family, hashes):
    out={}
    for stem,gate in (('outer_angular_repair','actual_implicit_functional_angular_pressure_repair_independently_checked'),
                      ('future_swirl_energy','whole_Z_positive_C1_complete_energy_checked')):
        name,receipt=PREFIX+stem+'.json',PREFIX+stem+'_check.json'
        raw,check=parent.linear.read(name),parent.linear.read(receipt)
        if not check.get('all_passed') or not check.get(gate):
            raise ValueError('Accepted original angular/energy source required')
        for key,value in family.items():
            if raw.get(key)!=value or (key in check and check[key]!=value):
                raise ValueError('Original angular/energy family differs')
        for record in (raw,check,dict(input_hashes={name:sha(name),receipt:sha(receipt)})):
            for path,digest in record['input_hashes'].items():
                if sha(path)!=digest or (path in hashes and hashes[path]!=digest):
                    raise ValueError('Changed angular/energy source '+path)
                hashes[path]=digest
        out[stem]=dict(raw=raw,receipt=receipt,gate=gate)
    return out


def small_branch_certificate(raw):
    """Uniform bound for the exact rationalized root on the admitted ball."""
    c=MPIntervalContext();c.dps=200
    def box(row):
        if 'exact_mpf_tuple' in row:
            q=tuple(row['exact_mpf_tuple']);return c.make_mpf((q,q))
        return c.make_mpf((tuple(row['lower_exact_mpf_tuple']),tuple(row['upper_exact_mpf_tuple'])))
    with mp.workdps(240):
        W={key:box(value) for key,value in raw['bump_weights'].items()}
        scale,R,ball=(box(raw[key]) for key in ('common_coefficient_scale','scaled_rhs_sup','scaled_coefficient_ball_radius'))
        d=-box(raw['rescaled_linear_determinant'])
        kappa=W['D']/(2*W['B'])
        b1=c.exp(1)*R/W['A'];b2=R/W['B']
        Amax=2*kappa*scale
        Bmin=d-2*kappa*scale*b1
        Cmax=b1+kappa*scale*b1*b1+b2
        disc=Bmin*Bmin-4*Amax*Cmax
        ep=Gamma.pulse.parent.linear.ep
        if ep(Bmin)[0]<=0 or ep(disc)[0]<=0:
            raise ValueError('Small angular root branch not proved')
        xmax=2*Cmax/(Bmin+c.sqrt(disc));ymax=b1+xmax
        derivative=Bmin-2*Amax*ball
        if ep(xmax-ball)[1]>=0 or ep(ymax-ball)[1]>=0 or ep(derivative)[0]<=0:
            raise ValueError('Rationalized root not in admitted original ball')
        return Gamma.pulse.coefficients.record(dict(passed=True,B_positive=Bmin,
            discriminant_positive=disc,x_absolute_upper=xmax,y_absolute_upper=ymax,
            admitted_ball=ball,positive_quadratic_derivative=derivative,
            exact_root_in_original_unique_small_branch=True,
            uniform_whole_Z_source_bounds_used=True,source_caps_are_evidence_not_functions=True))


def source_identity_proof(owner):
    """Identify exact new RHSs with the functions enclosed by prior sources."""
    frame=parent.linear.read(parent.quiet.FRAME_NAME)
    receipt=parent.linear.read(parent.quiet.FRAME_RECEIPT)
    native=receipt['independent_actual_six_row_native_inlet']
    if not (native['current_complete_Rp_qshapes_identified_from_actual_leading_functions']
            and native['independent_actual_native_incoming_Taylor_rows']==42):
        raise ValueError('Actual native incoming-function identity required')
    Xp=sy.sympify(frame['actual_analytic_native_inlet_functions']['X'])
    if any(symbol.name=='Z' for symbol in sy.cancel(Xp).free_symbols):
        raise ValueError('Actual incoming angular ratio must be Z-independent')
    eta,K=sy.symbols('eta K',nonnegative=True)
    Q=1+eta
    q,Xv,decay,k=sy.symbols('q Xv decay k',positive=True)
    # The endpoint flatten factor is Q/2. Its integrating factor gives
    # Xf(Z)=2Xv exp(-100r)/Q + integral exp(-r(100-t))*2^k Q^-k dt.
    difference=(2*Xv*decay+2**k)-(2*Xv*decay/Q+2**k*Q**(-k))
    expected=2*Xv*decay*eta/Q+2**k*(1-Q**(-k))
    if sy.simplify(difference-expected)!=0:
        raise ValueError('Actual flatten difference identity failed')
    primitive=-(1+sy.Symbol('s')*eta)**(-K)/eta
    s=sy.Symbol('s')
    if sy.simplify(sy.diff(primitive,s)-K*(1+s*eta)**(-1-K))!=0:
        raise ValueError('Positive divided flatten primitive identity failed')
    if sy.simplify(eta*(primitive.subs(s,1)-primitive.subs(s,0))-(1-Q**(-K)))!=0:
        raise ValueError('Original divided flatten RHS differs')
    if sy.simplify(sy.diff(1-Q**(-K),eta)-K*Q**(-1-K))!=0:
        raise ValueError('True first-Z flatten source differs')
    # Actual production assignments determine the heat conversions and
    # future composition; the audit does not compare numerical boxes.
    mu,a,Ts,tau,lone=sy.symbols('mu a Ts tau lone',real=True)
    r,k=1-mu,1-a
    tail=-sy.Rational(3,2)*(2+Ts)+(r+k)/2-(sy.Rational(1,2)+a)*tau
    env={'self.rate':r,'self.angular.restore_rate':k,'self.angular.waiting':tau,
        'self.angular.waiting_logone':lone,'self.logEtail_over_Erel':tail,'self.delta':2*a}
    for target,want in (('self.log_theta_multiplier',(r+k)/2+k*tau-lone),
                        ('self.log_pressure_multiplier',2*tail-2*lone+sy.log(a))):
        got=assignment('compliant_outer_angular_repair','__init__',target,env)
        if sy.simplify(got-want)!=0:
            raise ValueError('Original heat conversion differs: '+target)
    heat=heat_bridge.same_exact_Gamma_future_bindings()
    if not heat['passed'] or not all(heat['exact_full_Gamma_future_identities'].values()):
        raise ValueError('Original complete Gamma enclosure bindings required')
    c=MPIntervalContext();c.dps=80
    ep=parent.linear.ep
    if ep(-4*(c.exp(40)+11)-30+200*c.ln(10))[1]>=0:
        raise ValueError('Original strict delta minimum branch not proved')
    # The identified native callback's Xp bound is independent of Z. Prove
    # its TRUE pulse memory is contained; never materialize exp(-13r/mu).
    c.dps=200
    power=parent.linear.read(PREFIX+'power_inlet_C4.json')
    candidate=parent.linear.read(PREFIX+'outer_angular_candidate.json')
    Xp_box=read_interval(c,power['whole_Z_inlet']['canonical_incoming_Taylor']['X']['coefficients'][0])
    def parameter(i):
        row=owner.g.nodes[i];op=row['operation']
        if op=='exact_rational':return c.mpf(row['numerator'])/row['denominator']
        if op=='sum':return sum((parameter(j) for j in row['arguments']),c.mpf(0))
        if op=='product':
            value=c.mpf(1)
            for j in row['arguments']:value*=parameter(j)
            return value
        if op=='negative':return -parameter(row['argument'])
        if op=='analytic_unary':
            return (c.exp if row['name']=='exp' else c.ln)(parameter(row['argument']))
        raise ValueError('Original scalar parameter recipe required')
    def symbolic_parameter(i):
        row=owner.g.nodes[i];op=row['operation']
        if op=='exact_rational':return sy.Rational(row['numerator'],row['denominator'])
        if op=='sum':return sum(symbolic_parameter(j) for j in row['arguments'])
        if op=='product':return sy.prod(symbolic_parameter(j) for j in row['arguments'])
        if op=='negative':return -symbolic_parameter(row['argument'])
        if op=='analytic_unary':
            return (sy.exp if row['name']=='exp' else sy.log)(symbolic_parameter(row['argument']))
        raise ValueError('Exact original scalar parameter expression required')
    actual_mu=symbolic_parameter(owner.schedule['mu'].node)
    expected_mu=sy.exp(sy.log(sy.Rational(1,1000))-4*(sy.exp(40)+11))
    if sy.simplify(actual_mu/expected_mu-1)!=0:
        raise ValueError('Actual graph mu differs from original selected source recipe')
    mu_box=parameter(owner.schedule['mu'].node);r_box=1-mu_box;eq=1/r_box
    selected_mu=read_interval(c,power['selected_mu'])
    if max(ep(mu_box)[0],ep(selected_mu)[0])>min(ep(mu_box)[1],ep(selected_mu)[1]):
        raise ValueError('Selected mu enclosure contradicts its exact defining recipe')
    absdiff=max(mp.make_mpf((0,value[1],value[2],value[3])) for value in (Xp_box-eq)._mpi_)
    log_memory=c.ln(c.mpf(absdiff))-13*r_box/mu_box if absdiff else c.mpf(-1001)
    if ep(log_memory)[1]>=-1000:
        raise ValueError('Actual native pulse memory not in original cap')
    cap=c.exp(-1000)
    upper=cap._mpi_[1]
    Xv_box=eq+c.make_mpf(((1,upper[1],upper[2],upper[3]),upper))
    original_Xv=read_interval(c,candidate['angular_Xv_enclosure'])
    if not ep(original_Xv)[0]<=ep(Xv_box)[0]<=ep(Xv_box)[1]<=ep(original_Xv)[1]:
        raise ValueError('Actual native Xv enclosure not in archived candidate')
    flat,power_piece,bump,post,heat_piece=sy.symbols('flatten power bump post heat')
    total=assignment('compliant_fifth_axial_jets','future','total',dict(
        flatten=flat,power=power_piece,angular_change=bump,postrel=post,heatdifference=heat_piece))
    if sy.simplify(total-(flat+power_piece+bump+post-heat_piece))!=0:
        raise ValueError('Complete C5 future decomposition differs')
    Q,L,dj,E,F,center=sy.symbols('Q L dj E F center',real=True)
    env={'q':Q,'box(f.mu)':mu,'box(f.Lrel)':L,'dj':dj,
        "box(w['E'])":E,"box(w['F'])":F,'center':center}
    termwise={}
    for target,want in (('Nf',Q**2*sy.exp(-200*mu)/4),
                        ('Nrel',Q**2*sy.exp(-2*mu*(100+L))/4)):
        got=assignment('compliant_fifth_axial_jets','future',target,env)
        if sy.simplify(got-want)!=0:raise ValueError('Actual C5 '+target+' differs')
        termwise[target]=True
    got=assignment('compliant_fifth_axial_jets','future','angular_change',env,augmented=True)
    if sy.simplify(got-(2*E*dj+F*dj*dj)*sy.exp(-2*mu*center))!=0:
        raise ValueError('Actual C5 angular energy increment differs')
    termwise['signed_angular_increment']=True
    return dict(passed=True,actual_native_Xp_is_Z_independent=True,
        exact_native_incoming_source_projection_receipt=parent.quiet.FRAME_RECEIPT,
        true_pulse_memory_is_enclosed_by_original_source_candidate=True,
        flatten_value_and_first_Z_equal_original_positive_divided_integral=True,
        eta_zero_continuation_of_divided_integral_retained=True,
        complete_Gamma_value_and_first_Z_source_integrands_identified=True,
        actual_nonzero_S_same_waiting_and_radius_recipe=True,
        original_angular_heat_conversion_assignments_replayed=True,
        active_delta_minimum_branch_proved_for_actual_Md40=True,
        actual_native_Xv_inside_archived_candidate_interval=True,
        graph_mu_exact_recipe_matches_original_selected_parameter=True,
        actual_graph_mu_enclosure=Gamma.pulse.coefficients.record(mu_box),
        pulse_memory_containment_log_bound_strictly_below_minus_1000=True,
        actual_C5_future_has_same_five_signed_complete_pieces=True,
        actual_C5_termwise_source_identities=termwise,
        original_scaled_rhs_bound_applies_to_new_exact_functions=True,
        heat_source_assignment_evidence=heat)


def pulse_branch_certificate(owner, source_identity):
    """Fresh whole-Z interval proof for the actual same-function inputs."""
    if not source_identity['passed']:
        raise ValueError('Actual complete future identity required')
    c=MPIntervalContext();c.dps=200
    ep=parent.linear.ep
    with mp.workdps(240):
        read=lambda record:read_interval(c,record)
        basis=owner.pulse.basis_raw['view']
        W=basis['bump_basis'];mu=read(basis['original_mu'])
        graph_mu=read(source_identity['actual_graph_mu_enclosure'])
        if max(ep(mu)[0],ep(graph_mu)[0])>min(ep(mu)[1],ep(graph_mu)[1]):
            raise ValueError('Pulse basis mu contradicts exact original graph recipe')
        A=list(map(read,W['first_row']));D=list(map(read,W['exact_divided_difference_row']))
        det=read(W['divided_difference_determinant'])
        if ep(det)[1]>=0:raise ValueError('Original negative pulse determinant required')
        saddle=parent.linear.read(Gamma.pulse.basis.saddle.NAME)['view']
        Pi=list(map(read,saddle['refined_scaled_full_rows']))
        raw=parent.linear.read(PREFIX+'power_inlet_C4.json')
        native=raw['whole_Z_inlet']['canonical_incoming_Taylor']
        incoming=[read(native[key]['coefficients'][0]) for key in ('m1','m2')]
        energy=read(native['energy']['coefficients'][0])
        amplitude=parent.linear.read(PREFIX+'axial_amplitude_selection.json')
        factors=list(map(read,amplitude['actual_incoming_row_factor_positive_caps']))
        rows=[f*v for f,v in zip(factors,incoming)]
        inverse=lambda x,y:Gamma.pulse.coefficients.linear_inverse(A,D,det,x,y)
        u=inverse(-mu*rows[0],-(rows[1]-rows[0]))
        v=inverse(-mu*Pi[0],-(Pi[1]-Pi[0]))
        future=parent.linear.read(PREFIX+'future_swirl_energy.json')['whole_Z_C1']
        weighted=read(future['Section7_34_weighted_future_Taylor']['coefficients'][0])
        target=weighted-mu*energy+(1-c.exp(-26))/4
        cap=read(amplitude['positive_end_energy_cap'])
        upper=cap._mpi_[1]
        nu=c.make_mpf(((0,0,0,0),upper))
        K=read(basis['fixed_pulse_energy_Kpulse'])
        A2=K+sum((nu*j*j for j in v),c.mpf(0))
        A1=sum((2*nu*j*k for j,k in zip(u,v)),c.mpf(0))
        A0=-target+sum((nu*j*j for j in u),c.mpf(0))
        if ep(A2)[0]<=0 or ep(A0)[1]>=0:
            raise ValueError('Actual original pulse root conditions not proved')
        half=A1/(2*A2);disc=half*half-A0/A2
        if ep(disc)[0]<=0:raise ValueError('Actual completed-square discriminant not positive')
        ap=-half+c.sqrt(disc);denominator=2*A2*ap+A1
        low=A2*c.mpf('.9')**2+A1*c.mpf('.9')+A0
        high=A2*c.mpf('1.2')**2+A1*c.mpf('1.2')+A0
        C=[j+ap*k for j,k in zip(u,v)]
        if not mp.mpf('.9')<ep(ap)[0]<=ep(ap)[1]<mp.mpf('1.2'):
            raise ValueError('Actual pulse amplitude outside original bracket')
        if ep(low)[1]>=0 or ep(high)[0]<=0 or ep(denominator)[0]<=0:
            raise ValueError('Actual whole-Z pulse bracket/derivative denominator lost')
        if ep(C[0])[1]>=0 or ep(C[1])[0]<=0:
            raise ValueError('Actual pulse end coefficient signs not proved')
        return Gamma.pulse.coefficients.record(dict(passed=True,whole_Z_domain=[-1,1],
            fresh_actual_native_incoming_and_complete_future_bounds_used=True,
            A2=A2,A1=A1,A0=A0,discriminant=disc,ap=ap,C1=C[0],C2=C[1],
            positive_derivative_denominator=denominator,
            negative_lower_endpoint_residual=low,positive_upper_endpoint_residual=high,
            true_nonzero_nu_enclosed_not_selected_zero=True,
            nu_lower_zero_is_only_a_conservative_enclosure=True,
            actual_nu_strictly_positive_by_exact_product_mu_Kj_exp_2logE=True,
            exact_nonzero_nu_function_nodes=parent.encode(owner.pulse.functions['exact_nonzero_end_energy']),
            graph_and_basis_mu_same_original_recipe_enclosures_consistent=True,
            original_positive_branch_point_nine_to_one_point_two_proved=True,
            coefficient_first_Z_functions_from_exact_equations=True,
            derivative_accuracy_not_claimed=True))


class CurrentOriginalN1CompleteFutureFunctions:
    def __init__(self,Gamma_owner=None,require_checked=True):
        self.Gamma=Gamma_owner or Gamma.CurrentOriginalN1GammaFutureFunctions()
        if not self.Gamma.acceptance_loaded or not self.Gamma.pulse.acceptance_loaded:
            raise ValueError('Checked original Gamma and native pulse owners required')
        self.Gamma.assert_graph()
        self.g,self.family=self.Gamma.g,copy.deepcopy(self.Gamma.family)
        self.hashes=dict(self.Gamma.hashes)
        for module in (Gamma,Gamma.pulse):
            for path,digest in parent.linear.read(module.RECEIPT)['input_hashes'].items():
                if sha(path)!=digest or (path in self.hashes and self.hashes[path]!=digest):
                    raise ValueError('Changed checked ancestor source '+path)
                self.hashes[path]=digest
            self.hashes[module.NAME]=sha(module.NAME)
            self.hashes[module.RECEIPT]=sha(module.RECEIPT)
        records=checked_sources(self.family,self.hashes)
        self.source_receipts={key:dict(raw=PREFIX+key+'.json',receipt=value['receipt'],gate=value['gate'])
            for key,value in records.items()}
        self.source_identity=source_identity_proof(self.Gamma)
        self.branch_certificate=small_branch_certificate(records['outer_angular_repair']['raw'])
        self.weights=angular_weights(self.g,self.Gamma.schedule['mu'],self.Gamma.pulse.functions['normalization'])
        self.angular,self.angular_equations=angular_functions(self.g,self.Gamma.schedule,
            self.Gamma.histories,self.Gamma.functions,self.weights)
        self.functions,self.pieces=complete_energy_functions(self.g,self.Gamma.schedule,
            self.Gamma.functions,self.angular,self.weights)
        self.pulse_branch_certificate=pulse_branch_certificate(self.Gamma,self.source_identity)
        self.pulse_functions,self.pulse_equations=self.Gamma.pulse.original_coefficient_functions(self.functions['weighted'])
        for module in (Gamma,angular_source,energy_source,heat_bridge):
            self.hashes[Path(module.__file__).name]=sha(Path(module.__file__).name)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(self.g.nodes);self.snapshot=self._snapshot()
        self.acceptance_loaded=False
        if require_checked:
            receipt,saved=parent.linear.read(RECEIPT),parent.linear.read(NAME)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or not receipt.get(ANGULAR) or not receipt.get(PULSE)
                    or receipt['source_family']!=self.family or any(receipt.get(key) for key in OPEN)
                    or saved['exact_graph_nodes']!=self.g.nodes or saved['function_bindings']!=self.snapshot):
                raise ValueError('Checked original complete future C1 graph required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed complete future source '+name)
            self.acceptance_loaded=True

    def _snapshot(self):
        return parent.encode(dict(weights=self.weights,angular=self.angular,
            angular_equations=self.angular_equations,energy=self.functions,pieces=self.pieces,
            selected_pulse=self.pulse_functions,pulse_equations=self.pulse_equations,
            source_receipts=self.source_receipts,source_identity=self.source_identity,
            branch_certificate=self.branch_certificate,pulse_branch_certificate=self.pulse_branch_certificate,
            input_hashes=self.hashes))

    def assert_graph(self):
        self.Gamma.assert_graph()
        if self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot:
            raise ValueError('Original complete future function graph changed')


def run():
    began=time.monotonic();owner=CurrentOriginalN1CompleteFutureFunctions(require_checked=False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,function_bindings=owner.snapshot,
        actual_small_angular_root_not_a_box_center=True,exact_nonzero_heat_energy_retained=True,
        entire_Gamma_exterior_and_both_epsilon_atoms_retained=True,
        heat_Nrel_prefactor_occurs_exactly_once=True,pulse_branch_certificate_defined=True,
        **dict.fromkeys((GATE,ANGULAR,PULSE)+OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_COMPLETE_FUTURE_C1_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw


if __name__=='__main__':run()
