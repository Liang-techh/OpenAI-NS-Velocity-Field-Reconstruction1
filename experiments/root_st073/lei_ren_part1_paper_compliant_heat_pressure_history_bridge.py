"""Bind the C4 forward pressure to the original absolute history.

Every symbolic atom below names a production integral, not a fitted value.
The key identity is obtained by composing the actual forward assignments
with the original future decomposition and splitting the SAME collar
integral. Neither Ptail nor P3 is defined by the desired negative tail.
"""
import ast
import hashlib
from pathlib import Path

import sympy as s
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment
from lei_ren_part1_paper_compliant_power_inlet_C4_check import source_identities as inlet_identities
from lei_ren_part1_paper_compliant_collar_Gamma_C4_check import source_branch_expression

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def pressure_history_bridge():
    proofs, hashes = {}, {}

    def tree(stem):
        path = HERE/(PREFIX+stem+'.py')
        hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        return ast.parse(path.read_text(encoding='utf8'))

    def syntax(stem, method, target, expected, augmented=False):
        fn = next(n for n in ast.walk(tree(stem)) if isinstance(n, ast.FunctionDef) and n.name == method)
        nodes = [n.value for n in ast.walk(fn)
                 if (not augmented and isinstance(n, ast.Assign)
                     and any(ast.unparse(t) == target for t in n.targets))
                 or (augmented and isinstance(n, ast.AugAssign)
                     and ast.unparse(n.target) == target and isinstance(n.op, ast.Add))]
        exprs = expected if isinstance(expected, list) else [expected]
        if len(nodes) != len(exprs) or any(ast.dump(n) != ast.dump(ast.parse(e, mode='eval').body)
                                         for n, e in zip(nodes, exprs)):
            raise ValueError('Pressure history defining source changed: '+stem+':'+target)
        proofs['production_'+stem+'_'+method+'_'+target] = True

    def imported(stem, module, name, alias=None):
        if not any(isinstance(n, ast.ImportFrom) and n.module == PREFIX+module
                   and any(a.name == name and a.asname == alias for a in n.names)
                   for n in tree(stem).body):
            raise ValueError('Pressure history source import changed: '+stem+':'+name)
        proofs['import_'+stem+'_'+name] = True

    def actual(stem, method, target, env, augmented=False):
        tree(stem)
        return assignment(stem, method, target, env, augmented)

    def zero(name, expr):
        # These are exact, already identified integrals. Preserve each
        # identical integral as one atom while checking coefficient algebra;
        # do not ask Sympy to evaluate the Gamma/Piecewise integrals.
        expr = expr.xreplace({integral:s.Dummy('source_integral')
                             for integral in expr.atoms(s.Integral)})
        if s.simplify(s.expand_power_exp(expr)) != 0:
            raise ArithmeticError('Pressure history identity failed: '+name)
        proofs[name] = True

    # This lemma binds the O2/O3 expressions and proves Pin/q^2 and U/q
    # on the whole Z domain. The Z=.5 receipt encloses a proved constant;
    # it is not a sampled polynomial or an independently selected inlet.
    inherited_inlet = inlet_identities()
    if not inherited_inlet['actual_O2_O3_source_chain_has_exact_canonical_whole_Z_shapes']:
        raise ValueError('Canonical original inlet shapes not proved')
    tree('compliant_power_inlet_C4_check')
    for stem in ('compliant_outer_initial', 'compliant_outer_buffer'):
        tree(stem)
    syntax('compliant_power_inlet_C4', '__init__', 'sample', "buffer['samples'][-1]")
    syntax('compliant_power_inlet_C4', '__init__', 'q', "c.mpf('1.25')")
    syntax('compliant_power_inlet_C4', '__init__', 'self.inlet_P',
           "read_interval(c,sample['Mp_over_Pstar_squared'][0])*q*q")
    syntax('compliant_outer_buffer', 'report', 'samples',
           "[self.power('.5',t) for t in ('0','.5','1')]", augmented=True)
    syntax('compliant_power_inlet_C4', 'incoming', 'raw',
           "self.datum.normalized_jets(Z,5)['normalized_pressure_coefficients']")
    imported('compliant_power_inlet_C4', 'compliant_pressure_source', 'CompliantPressureDatum')
    imported('compliant_five_moment_repair', 'compliant_pressure_source',
             'CompliantPressureDatum', 'LogarithmicPressureDatum')
    syntax('compliant_power_inlet_C4', '__init__', 'self.datum', "CompliantPressureDatum('40',precision=160)")
    syntax('compliant_five_moment_repair', '__init__', 'self.datum', "LogarithmicPressureDatum('40',precision=160)")
    syntax('compliant_outer_initial', '__init__', 'self.repair', 'SharedFiveMomentRepair()')
    syntax('compliant_outer_initial', '__init__', 'self.datum', 'self.repair.datum')
    syntax('compliant_outer_buffer', '__init__', 'self.initial', 'SharedOuterInitial()')
    syntax('compliant_outer_angular_candidate', '__init__', 'self.initial', 'self.buffer.initial')
    syntax('compliant_outer_angular_repair', '__init__', 'self.angular', 'self.heat.angular')
    syntax('compliant_future_swirl_energy', '__init__', 'self.angular', 'self.repair.angular')
    syntax('compliant_axial_amplitude_selection', '__init__', 'self.future', 'CompliantFutureSwirlEnergy()')
    syntax('compliant_axial_pulse_field', '__init__', 'self.selection', 'CompliantAxialAmplitude()')
    syntax('compliant_axial_pulse_field', 'data', 'inlet', 'self.pulse.buffer.power(Z,1)')
    syntax('compliant_axial_pulse_field', 'pressure_moment', 'p0rows',
           "self.selection.future.angular.initial.datum.normalized_jets(c.mpf(Z),1)['normalized_pressure_coefficients']")

    mu, a, L, Ts, wait, lone, z = s.symbols('mu a Lrel Ts tau lone Z', real=True)
    bp, bh, pr, p, rate, k = s.Rational(1,2)+mu, s.Rational(1,2)+a, 1+2*mu, 1+2*a, 1-mu, 1-a
    q = 1+z*z
    U, Pin = s.symbols('original_O3_U original_O3_Pin', real=True)
    P0 = s.Function('original_analytic_P0')(z)
    Ev2 = U**2*s.exp(-13/mu-26)
    decay = s.exp(-13/mu-26)
    MpRv = actual('compliant_flatten_mixed_C4', 'flatten', 'Mp_v',
                  {"data['Mp']":Pin/q**2, "data['u']":U/q,
                   'self.pressure_decay':decay, 'self.prate':pr})
    pulseMp = actual('compliant_axial_pulse_field', 'pressure_moment', 'p',
                     {"IntervalTaylor(c, inlet['Mp_over_Pstar_squared'])":Pin/q**2,
                      'u':U/q, 'kernel':(1-decay)/pr})
    zero('same_original_MpRv_not_a_negative_tail_reset', MpRv-pulseMp)
    syntax('compliant_axial_pulse_field', 'end', 't', '13/self.mu+s')
    endfn = next(n for n in ast.walk(tree('compliant_axial_pulse_field'))
                 if isinstance(n, ast.FunctionDef) and n.name == 'end')
    calls = [n for n in ast.walk(endfn) if isinstance(n, ast.Call)
             and ast.unparse(n.func) == 'self.pressure_moment']
    expected = ast.parse('self.pressure_moment(Z,t,inlet,u,dict(inverse_mu_term=-13*self.prate/self.mu,finite_offset=-self.prate*s))', mode='eval').body
    if len(calls) != 1 or ast.dump(calls[0]) != ast.dump(expected):
        raise ValueError('Original pulse terminal pressure coordinate changed')
    zero('pulse_Rv_exact_pressure_decay', -13*pr/mu-(-13/mu-26))
    proofs['production_original_pulse_terminal_pressure_coordinate'] = True

    # Bind the actual cell integrands. Exact exponential cell masses in
    # C4 and directed whole integrands in the closure enclose one integral.
    syntax('compliant_flatten_mixed_C4', 'flatten', 'Fc', '(rho*sigma_jets(c,v/100)[0]).exp()')
    syntax('compliant_flatten_mixed_C4', 'flatten', 'Pint',
           'Fc*Fc/(q*q)*(c.exp(-self.prate*a)*decay_integral(c,self.prate,length)/2)', augmented=True)
    syntax('compliant_corrected_outer_field', 'flatten_tails', 'F',
           '(ratio*(2*sigma_enclosure(c,v/100))).exp()/(q*q)')
    syntax('compliant_corrected_outer_field', 'flatten_tails', 'pressure',
           'F*(ds*c.exp(-self.prate*v)/2)', augmented=True)
    sig, v = s.symbols('sigma v', real=True)
    rho = s.log(q/2)
    zero('same_original_flatten_pressure_density',
         s.exp(rho*sig)**2/q**2-s.exp(2*rho*sig)/q**2)
    syntax('compliant_flat_pulse_derivatives', '_sigma_left', 'odds', '1/(1-x)**2-1/x**2')
    syntax('compliant_axial_pulse_field', 'point', 'odds', '-1/y**2+1/(1-y)**2')
    x = s.symbols('x', real=True)
    zero('original_sigma_expression_shared',
         1/(1-x)**2-1/x**2-(-1/x**2+1/(1-x)**2))

    Iflat = s.Integral(s.exp(-pr*v+2*rho*s.Function('original_sigma')(v/100))/(2*q*q), (v,0,100))
    MpFlat = actual('compliant_flatten_mixed_C4', 'flatten', 'Mp',
                    {'Mp_v':MpRv,'Pint':Iflat,'self.Ev2':Ev2})
    Pflat = actual('compliant_flatten_mixed_C4', 'flatten', 'pressure',
                   {'Mp':MpFlat,"data['P0']":P0})
    syntax('compliant_power_angular_C4', 'data', 'f', 'self.flatten.flatten(Z,100)')
    syntax('compliant_steep_waiting_C4', 'data', 'terminal', 'self.outer.angular(Z,0)')
    syntax('compliant_steep_waiting_C4', 'data', 'PR', "terminal['pressure_over_Pstar_squared_Taylor']")
    syntax('compliant_collar_Gamma_C4', 'data', 'terminal', 'self.steep.waiting(Z,1)')
    syntax('compliant_collar_Gamma_C4', 'data', 'Ptail', "terminal['pressure_over_Pstar_squared_Taylor']")

    env = {'self.bp':bp,'self.bh':bh,'self.rate':rate,'self.k':k,
           'self.outer.Lrel':L,'self.Lrel':L,'self.Ts':Ts,'self.params.Ts':Ts,
           'self.mu':mu,'self.a':a,'self.prate':pr,'self.phrate':p,
           'self.angular.waiting':wait,'self.angular.waiting_logone':lone}
    weights = {}
    for name in ('self.Pf','self.Prel','self.Ps','self.Pq','self.Pt','self.Ptail','self.tailmult'):
        weights[name] = actual('compliant_corrected_outer_field','__init__',name,{**env,**weights})
    theta = {}
    for name in ('self.thetaR','self.thetaS','self.thetaQ','self.thetaT'):
        theta[name] = actual('compliant_steep_waiting_C4','__init__',name,{**env,**theta})
    theta_base = actual('compliant_collar_Gamma_C4','__init__','self.theta_base',
                        {'self.steep.thetaT':theta['self.thetaT'],'self.bh':bh,
                         'self.steep.wait':wait,'self.steep.logone':lone})
    for label, weight in (('R',weights['self.Prel']),('S',weights['self.Prel']*weights['self.Ps']),
                           ('Q',weights['self.Prel']*weights['self.Pq']),('T',weights['self.Prel']*weights['self.Pt'])):
        zero('same_actual_pressure_amplitude_'+label,theta['self.theta'+label]**2-weight)
    zero('same_actual_pressure_amplitude_heat',theta_base**2-
         weights['self.Prel']*weights['self.Ptail']*weights['self.tailmult'])
    D = lambda r,t:(1-s.exp(-r*t))/r
    # BP is the actual sum of the two original bump integrals; its origin
    # and signed formula are bound in both evaluators below.
    BP, Pintr, Poutr = s.symbols('actual_bump_pressure actual_entry_integral actual_exit_integral', real=True)
    syntax('compliant_power_angular_C4','angular','pastP',
           "(dj*past['B']+dj*dj*past['D']/2)*c.exp(-self.prate*center)",augmented=True)
    syntax('compliant_corrected_outer_field','correction','P',
           "(dj*w['B']+dj*dj*w['D']/2)*c.exp(-self.prate*center)",augmented=True)
    baseline = actual('compliant_power_angular_C4','angular','baseline',
                      {'decay_integral(c, self.prate, y)':D(pr,L),'self.flatten.Ev2':Ev2,'self.prate':pr})
    pangenv = {"data['flatten_exit_pressure']['P_over_Pstar_squared']":Pflat,
               'baseline':baseline,'pastP':BP,'self.flatten.Ev2':Ev2,'self.prate':pr,'self.Lrel':L}
    PR = actual('compliant_power_angular_C4','angular','pressure',pangenv)
    PR += actual('compliant_power_angular_C4','angular','pressure',pangenv,True)
    pscale = Ev2*theta_base**2
    ps = actual('compliant_steep_waiting_C4','data','PS',
                {'PR':PR,"self.infull['pressure']":Pintr,'self.outer.flatten.Ev2':Ev2,'self.thetaR':theta['self.thetaR']})
    pq = actual('compliant_steep_waiting_C4','data','PQ',
                {'PS':ps,'decay_integral(c, 3, self.Ts)':D(3,Ts),'self.outer.flatten.Ev2':Ev2,'self.thetaS':theta['self.thetaS']})
    pt = actual('compliant_steep_waiting_C4','data','PT',
                {'PQ':pq,"self.outfull['pressure']":Poutr,'self.outer.flatten.Ev2':Ev2,'self.thetaQ':theta['self.thetaQ']})
    Ptail = actual('compliant_steep_waiting_C4','waiting','pressure',
                   {"data['PT']":pt,'decay_integral(c, 1 + self.delta, t)':D(p,wait),
                    'self.outer.flatten.Ev2':Ev2,'self.thetaT':theta['self.thetaT']})
    # Forward and backward unit transition kernels use identical densities.
    for expression in ('ds*c.exp(-(1+2*mu)*v-2*rate*jc)/2','ds*c.exp(-3*v+2*rate*jc)/2'):
        for stem,method in (('compliant_steep_waiting_C4','transition_kernels'),('compliant_corrected_outer_field','transition_tails')):
            fn = next(n for n in ast.walk(tree(stem)) if isinstance(n,ast.FunctionDef) and n.name==method)
            nodes = [ast.dump(n.value) for n in ast.walk(fn) if isinstance(n,ast.AugAssign) and ast.unparse(n.target)=='pressure']
            if ast.dump(ast.parse(expression,mode='eval').body) not in nodes:
                raise ValueError('Original transition pressure density changed: '+stem)
    proofs['same_exact_entry_and_exit_pressure_integrands'] = True

    # The common collar bracket and infinite-tail decomposition are bound
    # to their production expressions by the previously admitted collar
    # source identities. Rebind the pressure expression here explicitly.
    S, That, Phat, Ehat, JW, EW, EW2, PW, PW2, eps = s.symbols('S ThetaHat PressureHat EnergyHat JW EW EW2 PW PW2 epsilon',real=True)
    tailenv = {'one':1,'self.k':k,'self.S':S,'self.a':a,'self.eps':eps,'self.delta':2*a,'self.prate':p,
               't':0,'angular':That,'energy':Ehat,'pressure':Phat,
               "atoms['JW']":JW,"atoms['EW']":EW,"atoms['EW2']":EW2,"atoms['PW']":PW,"atoms['PW2']":PW2}
    new_tail0 = actual('compliant_collar_Gamma_C4','collar_tails','P',tailenv)
    oldenv = {'theta_hat':That,'energy_hat':Ehat,'pressure_hat':Phat,'S':S,'self.a':a,'self.eps':eps,
              'self.eps ** 2':eps**2,'self.delta':2*a,'self.phrate':p,'self.k':k,'t':0,
              "atoms['JW']":JW,"atoms['EW']":EW,"atoms['EW2']":EW2,"atoms['PW']":PW,"atoms['PW2']":PW2}
    old_tail0 = source_branch_expression('heat_tails','P','pressure_hat *',oldenv)
    zero('same_production_complete_collar_pressure_tail',new_tail0-old_tail0)
    # Identify the actual K in both point and tail evaluators. At row zero
    # product_rows is ordinary multiplication, so these exact production
    # expressions have the scalar bracket below. Higher rows are Leibniz
    # derivatives of the same bracket, not new functions.
    for target,expression in (
            ('W','1-sig+sig*phi'),('C','1-self.eps*phi'),('pre','1-self.eps*W'),
            ('D',"IntervalTaylor(c,self.heat.deficit(Z,t)['deficit_scaled_Taylor'].coefficients)*(sig*C)")):
        syntax('compliant_corrected_outer_field','collar_shape',target,expression)
    for target,expression in (
            ('W','[unity[j]-sr[j]+v for j,v in enumerate(product_rows(sr,fr))]'),
            ('C','[unity[j]-fr[j]*self.eps for j in range(5)]'),
            ('pre','[unity[j]-W[j]*self.eps for j in range(5)]'),
            ('D',['product_rows(product_rows(sr,C),dr)','[sr[0]*C[0]*dr[0]]']),
            ('dr','gamma_deficit_mixed(c,Z,self.a,self.Scap,t,4 if high else 0)'),
            ('K','[pre[j]-D[j]*(self.a*self.S) for j in range(len(D))]')):
        syntax('compliant_collar_Gamma_C4','shape',target,expression)
    shfn = next(n for n in ast.walk(tree('compliant_corrected_outer_field'))
                if isinstance(n,ast.FunctionDef) and n.name=='collar_shape')
    oldKs = [n.value for n in ast.walk(shfn) if isinstance(n,ast.keyword) and n.arg=='K']
    if len(oldKs)!=1 or ast.dump(oldKs[0])!=ast.dump(ast.parse('self.jet(pre)-D*(self.a*S)',mode='eval').body):
        raise ValueError('Original heat collar bracket return changed')
    syntax('compliant_collar_Gamma_C4','collar_tails','shape','self.shape(Z,v,False)')
    syntax('compliant_collar_Gamma_C4','collar_tails','square','D*(pre*2-D*(self.a*self.S))')
    syntax('compliant_collar_Gamma_C4','collar_tails','pressure',
           'square*(ds*c.exp(-self.prate*v)/2)',augmented=True)
    syntax('compliant_corrected_outer_field','heat_tails','sh','self.collar_shape(Z,v)')
    syntax('compliant_corrected_outer_field','heat_tails','sq',"zero+2*sh['pre']-D*(self.a*S)")
    syntax('compliant_corrected_outer_field','heat_tails','pressure_hat',
           'D*sq*(ds*c.exp(-self.phrate*v)/2)',augmented=True)
    # The epsilon atoms integrate W and W^2 with the same pressure weight.
    syntax('compliant_collar_Gamma_C4','collar_tails',"atoms[label + 'W']",'ds*c.exp(-rate*v)*W',augmented=True)
    syntax('compliant_collar_Gamma_C4','collar_tails',"atoms[label + 'W2']",'ds*c.exp(-rate*v)*W**2',augmented=True)
    sig0,phi0,H0,D0 = s.symbols('sigma phi canonical_Gamma_H Dhat',real=True)
    W0 = 1-sig0+sig0*phi0
    pre0 = 1-eps*W0
    Dcollar = sig0*(1-eps*phi0)*D0
    K0 = pre0-a*S*Dcollar
    zero('actual_collar_K_is_original_full_Gamma_bracket',
         K0.subs(D0,(1-H0)/(a*S))-((1-sig0)*(1-eps)+sig0*H0*(1-eps*phi0)))
    zero('actual_collar_pressure_deficit_and_epsilon_atom_density',
         K0**2-(1-2*eps*W0+eps**2*W0**2-a*S*Dcollar*(2*pre0-a*S*Dcollar)))
    zero('actual_full_Gamma_pressure_deficit_density',
         (1-a*S*D0)**2-(1-a*S*D0*(2-a*S*D0)))
    zero('actual_flat_collar_K_equals_Gamma_on_entire_exterior',
         K0.subs({sig0:1,phi0:0,D0:(1-H0)/(a*S)})-H0)
    syntax('compliant_collar_Gamma_C4','forward_pressure','K',"self.shape(Z,v,False)['K_rows'][0]")
    syntax('compliant_collar_Gamma_C4','forward_pressure','integral',
           'K*K*(ds*c.exp(-self.prate*v)/2)',augmented=True)
    fn = next(n for n in ast.walk(tree('compliant_collar_Gamma_C4')) if isinstance(n,ast.FunctionDef) and n.name=='forward_pressure')
    returns = [n.value for n in ast.walk(fn) if isinstance(n,ast.Return)]
    if len(returns)!=1 or ast.dump(returns[0])!=ast.dump(ast.parse('Ptail+integral*self.pressure_scale',mode='eval').body):
        raise ValueError('C4 pressure forward return changed')
    proofs['production_C4_forward_integral_added_to_retained_Ptail'] = True
    syntax('compliant_collar_Gamma_C4','data','pressure3','self.forward_pressure(Z,3,Ptail)')
    # Use actual integrals for the common bracket, not free P3/A3 constants.
    # Here K is the just-bound source expression. H is the canonical
    # expectation whose evaluator/radius bridge is separately required by
    # the caller. The generic symbol is used only for algebra above.
    sigma = s.Piecewise((0,v<=0),(1,v>=1),
                        (1/(1+s.exp(1/v**2-1/(1-v)**2)),True))
    phi = s.Piecewise((0,v>=3),(s.exp(-4/(3-v)**2),True))
    H = s.Function('canonical_positive_Gamma_H')(2*(1-z*z)*S*s.exp(-v))
    collar = (1-sigma)*(1-eps)+sigma*H*(1-eps*phi)
    Ic = s.Integral(s.exp(-p*v)*collar**2/2,(v,0,3))
    Iex = s.Integral(s.exp(-p*v)*collar**2/2,(v,3,s.oo))
    full_tail0 = Ic+Iex
    # K is the exact Gamma bracket for v>=3 (original flat phi=0,
    # sigma=1); the pressure tail rescales by substitution v=3+w.
    proofs['complete_production_collar_tail_identified_by_pointwise_density_and_integral_linearity'] = True
    proofs['collar_integral_split_uses_same_K_on_0_3_and_3_infinity'] = True
    P3 = Ptail+pscale*Ic
    postenv = {"heat0['pressure']":full_tail0,'self.Ptail':weights['self.Ptail'],
               'self.tailmult':weights['self.tailmult'],"self.fullin['pressure']":Pintr,
               'self.Ps':weights['self.Ps'],'decay_integral(c, 3, self.params.Ts)':D(3,Ts),
               'self.Pq':weights['self.Pq'],"self.fullout['pressure']":Poutr,
               'self.Pt':weights['self.Pt'],'decay_integral(c, self.phrate, self.angular.waiting)':D(p,wait)}
    Ppost = actual('compliant_corrected_outer_field','data','Ppost',postenv)
    Pf = actual('compliant_corrected_outer_field','data','Pf',
                {'Ppost':Ppost,'self.Prel':weights['self.Prel'],"future['P']":BP,
                 'self.Pf':weights['self.Pf'],'decay_integral(c, self.prate, self.Lrel)':D(pr,L)})
    Prv = actual('compliant_corrected_outer_field','data',"item['Prv']",{'Pf':Pf,'fp':Iflat})
    zero('actual_forward_P3_plus_actual_remaining_integral_is_original_infinity_offset',
         P3+pscale*Iex-(P0+MpRv+Ev2*Prv))
    # Bind the original closure packet offset; it is the accumulated
    # original Mp/P0 plus the same Prv, not a newly chosen pressure datum.
    prefix = source_branch_expression('packet','prefix',"data['Prv'] - Pglobal",
                                      {"data['Prv']":Prv,'Pglobal':theta_base**2*Iex})
    Mp3 = actual('compliant_corrected_outer_field','packet','Mp',
                 {"IntervalTaylor(c, terminal['Mp_over_Pstar_squared'].coefficients)":MpRv,
                  'prefix':prefix,'self.Ev0Pstar2':Ev2})
    oldP3 = actual('compliant_corrected_outer_field','packet','Pforward',{'Mp':Mp3,'P0':P0})
    zero('actual_C4_P3_equals_original_retained_forward_packet',P3-oldP3)
    # Once the ORIGINAL infinity offset is zero, the same integral yields
    # equivalence on the whole exterior and therefore all axial derivatives.
    u = s.symbols('offset',real=True)
    Iu = s.Integral(s.exp(-p*v)*collar**2/2,(v,u,s.oo))
    from3 = Iex-Iu
    zero('entire_exterior_history_offset_is_the_original_constant',
         P3+pscale*from3+pscale*Iu-(P0+MpRv+Ev2*Prv))
    return dict(identities=proofs, canonical_inlet_source_lemma=inherited_inlet,
                complete_retained_pressure_history_bridge_verified=True,
                offset_identity='C4_P3+C*integral_3^infinity=original_P0+original_MpRv+Ev2*original_Prv',
                original_pressure_datum_reset=False, interval_overlap_used_as_proof=False,
                full_collar_integral_definition=str(full_tail0), input_hashes=hashes)
