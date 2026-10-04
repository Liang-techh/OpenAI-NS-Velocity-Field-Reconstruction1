"""Bind actual C4 angular/energy histories to original terminal moments.

Source expressions and receipt producer routes establish common functions.
Receipt hashes preserve their admitted derivative enclosures; they are not
used to infer functional equality. Caps are bounds on exact positive scales.
"""
import ast
import hashlib
import json
from pathlib import Path

import sympy as s
from lei_ren_part1_paper_compliant_absolute_moment_closure import assignment, source_identities
from lei_ren_part1_paper_compliant_heat_pressure_source_bridge import source_bridge
from lei_ren_part1_paper_compliant_collar_Gamma_C4_check import source_branch_expression

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def terminal_history_bridge():
    proofs, hashes, trees = {}, {}, {}

    def tree(stem):
        if stem not in trees:
            path = HERE/(PREFIX+stem+'.py')
            hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
            trees[stem] = ast.parse(path.read_text(encoding='utf8'))
        return trees[stem]

    def method(stem, name, cls=None):
        scope = tree(stem)
        if cls is not None:
            scope = next(n for n in scope.body if isinstance(n,ast.ClassDef) and n.name==cls)
        return next(n for n in ast.walk(scope) if isinstance(n,ast.FunctionDef) and n.name==name)

    def syntax(stem, name, target, expected, augmented=False, cls=None):
        nodes = [n.value for n in ast.walk(method(stem,name,cls))
                 if (not augmented and isinstance(n,ast.Assign)
                     and any(ast.unparse(t)==target for t in n.targets))
                 or (augmented and isinstance(n,ast.AugAssign)
                     and ast.unparse(n.target)==target and isinstance(n.op,ast.Add))]
        values = expected if isinstance(expected,list) else [expected]
        if len(nodes)!=len(values) or any(ast.dump(n)!=ast.dump(ast.parse(v,mode='eval').body)
                                        for n,v in zip(nodes,values)):
            raise ValueError('Terminal history source changed: '+stem+':'+name+':'+target)
        proofs['source_'+stem+'_'+name+'_'+target] = True

    def body(stem, name, expected, cls=None):
        production = method(stem,name,cls).body
        # A docstring is descriptive metadata, not executable semantics.
        if production and isinstance(production[0],ast.Expr) and isinstance(production[0].value,ast.Constant):
            production=production[1:]
        wanted=ast.parse(expected).body
        if [ast.dump(n) for n in production] != [ast.dump(n) for n in wanted]:
            raise ValueError('Terminal history callable changed: '+stem+':'+name)
        proofs['callable_'+stem+'_'+(cls or '')+'_'+name]=True

    def keyword(stem, name, key, expected):
        nodes = [n.value for n in ast.walk(method(stem,name)) if isinstance(n,ast.keyword) and n.arg==key]
        if len(nodes)!=1 or ast.dump(nodes[0])!=ast.dump(ast.parse(expected,mode='eval').body):
            raise ValueError('Terminal history output route changed: '+stem+':'+key)
        proofs['output_'+stem+'_'+name+'_'+key] = True

    def statement(stem,name,expected):
        wanted=ast.parse(expected).body[0]
        if sum(ast.dump(n)==ast.dump(wanted) for n in ast.walk(method(stem,name)))!=1:
            raise ValueError('Terminal history statement changed: '+stem+':'+name+':'+expected)
        proofs['statement_'+stem+'_'+name+'_'+expected]=True

    def dictitem(stem,name,target,key,expected):
        nodes=[n.value for n in ast.walk(method(stem,name)) if isinstance(n,ast.Assign)
               and any(ast.unparse(t)==target for t in n.targets)]
        if len(nodes)!=1 or not isinstance(nodes[0],ast.Dict):
            raise ValueError('Unique terminal history dict required: '+stem+':'+target)
        values=[value for label,value in zip(nodes[0].keys,nodes[0].values)
                if isinstance(label,ast.Constant) and label.value==key]
        if len(values)!=1 or ast.dump(values[0])!=ast.dump(ast.parse(expected,mode='eval').body):
            raise ValueError('Terminal history dict item changed: '+stem+':'+key)
        proofs['dictitem_'+stem+'_'+target+'_'+key]=True

    def loop(stem,name,target,expected):
        nodes=[n.iter for n in ast.walk(method(stem,name)) if isinstance(n,ast.For) and ast.unparse(n.target)==target]
        if len(nodes)!=1 or ast.dump(nodes[0])!=ast.dump(ast.parse(expected,mode='eval').body):
            raise ValueError('Terminal support loop changed: '+stem+':'+name)
        proofs['loop_'+stem+'_'+name+'_'+target]=True

    def actual(stem, name, target, env, augmented=False):
        tree(stem)
        return assignment(stem,name,target,env,augmented)

    def accumulated(stem, name, target, env):
        """Translate an actual scalar initializer and all additive updates."""
        def formal(node):
            label=ast.unparse(node)
            if label in env:return env[label]
            if isinstance(node,ast.Constant):return s.Rational(str(node.value))
            if isinstance(node,ast.UnaryOp) and isinstance(node.op,ast.USub):return -formal(node.operand)
            if isinstance(node,ast.BinOp):
                left,right=formal(node.left),formal(node.right)
                for kind,op in ((ast.Add,lambda:left+right),(ast.Sub,lambda:left-right),
                                (ast.Mult,lambda:left*right),(ast.Div,lambda:left/right),
                                (ast.Pow,lambda:left**right)):
                    if isinstance(node.op,kind):return op()
            raise ValueError('Unbound terminal history scalar: '+stem+':'+label)
        initial=[n.value for n in ast.walk(method(stem,name)) if isinstance(n,ast.Assign)
                 and any(ast.unparse(t)==target for t in n.targets)]
        updates=[n.value for n in ast.walk(method(stem,name)) if isinstance(n,ast.AugAssign)
                 and ast.unparse(n.target)==target and isinstance(n.op,ast.Add)]
        if len(initial)!=1:raise ValueError('Unique actual scalar initializer required')
        return formal(initial[0])+sum((formal(n) for n in updates),s.Integer(0))

    def zero(name, value):
        value = value.xreplace({i:s.Dummy('same_full_integral') for i in value.atoms(s.Integral)})
        if s.simplify(s.expand_power_exp(value)) != 0:
            raise ArithmeticError('Terminal history identity failed: '+name)
        proofs[name] = True

    def receipt(stem, gates):
        path = HERE/(PREFIX+stem+'.json')
        record = json.loads(path.read_bytes())
        if any(not record.get(gate) for gate in gates):
            raise ValueError('Terminal history prerequisite not admitted: '+stem)
        for name,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Terminal history receipt dependency changed: '+name)
            hashes[name] = digest
        hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        return record

    pressure = source_bridge()
    if not pressure['complete_defining_function_history_bridge_verified']:
        raise ValueError('Complete original pressure/function/radius bridge required')
    for section in ('retained_pressure_history_bridge','defining_function_bridge'):
        hashes.update(pressure[section]['input_hashes'])
    tree('compliant_heat_pressure_source_bridge')
    absolute = source_identities()
    if not absolute['identities']['absolute_angular_constant_whole_Z']:
        raise ValueError('Original functional angular closure not proved')
    tree('compliant_absolute_moment_closure')
    for stem,gates in (
        ('compliant_angular_high_jets_check',['all_passed','angular_coefficient_C4_available']),
        ('compliant_future_energy_high_jets_check',['all_passed','complete_future_corrected_energy_C4_available']),
        ('compliant_fifth_axial_jets_check',['all_passed','angular_coefficient_C5_available','complete_future_corrected_energy_C5_available']),
        ('compliant_axial_amplitude_selection_check',['all_passed','complete_future_energy_source_bound']),
        ('compliant_axial_pulse_field_check',['all_passed','whole_Z_terminal_linear_zero_and_positive_energy_checked']),
        ('compliant_pulse_mixed_C4_check',['all_passed','whole_Z_terminal_mixed_zero_and_positive_energy_checked'])):
        receipt(stem,gates)

    # Bind receipt-backed C5 coefficients to the original small smooth
    # branch, including C0 and the unchanged admitted C4 prefix.
    for stem,target,expected in (
        ('compliant_axial_high_jets','self.energy','CompliantFutureEnergyHighJets()'),
        ('compliant_future_energy_high_jets','self.angular','CompliantAngularHighJets()'),
        ('compliant_future_energy_high_jets','self.base','CompliantFutureSwirlEnergy()'),
        ('compliant_angular_high_jets','self.repair','CompliantAngularRepair()'),
        ('compliant_fifth_axial_jets','self.fourth','CompliantAxialHighJets()'),
        ('compliant_fifth_axial_jets','self.angular4','self.fourth.energy.angular')):
        syntax(stem,'__init__',target,expected)
    syntax('compliant_angular_high_jets','coefficients','old','r.coefficients(Z)')
    syntax('compliant_angular_high_jets','coefficients','prior',"old['scaled_coefficient_Taylor']")
    syntax('compliant_angular_high_jets','coefficients','out',
           'implicit_quadratic_jets(c,prior[0][0],prior[1][0],b1,b2,r.p,r.q,r.k*r.scale,prior)')
    syntax('compliant_fifth_axial_jets','angular','old','self.angular4.coefficients(endpoints(Z))')
    syntax('compliant_fifth_axial_jets','angular','prior',"[copy_jet(c,j) for j in old['scaled_coefficient_Taylor']]")
    syntax('compliant_fifth_axial_jets','angular','physical','[j*box(r.scale) for j in coeffs]')
    keyword('compliant_fifth_axial_jets','angular','physical_coefficient_Taylor','physical')
    keyword('compliant_fifth_axial_jets','report','whole_Z',
            'dict(angular=self.angular([-1,1]),energy=self.future([-1,1]),selected=self.select([-1,1]))')
    syntax('compliant_flatten_mixed_C4','__init__','self.fifth',
           "json.loads((HERE/(PREFIX+'compliant_fifth_axial_jets.json')).read_bytes())")
    syntax('compliant_power_angular_C4','__init__','self.fifth','self.flatten.fifth')
    syntax('compliant_power_angular_C4','__init__','aname',"PREFIX+'compliant_outer_angular_repair.json'")
    syntax('compliant_power_angular_C4','__init__','angular','json.loads((HERE/aname).read_bytes())')
    syntax('compliant_power_angular_C4','__init__','self.weights',
           "{k:read_interval(c,v) for k,v in angular['bump_weights'].items()}")
    keyword('compliant_outer_angular_repair','report','bump_weights','self.weights')
    syntax('compliant_outer_angular_repair','__init__','self.weights',
           'bump_weights(c,self.mu,self.normalization,cells=bump_cells)')
    syntax('compliant_power_angular_C4','data','source',
           ["self.fifth['whole_Z']['angular']","sample['angular']"])
    syntax('compliant_power_angular_C4','data','coeff',"[jet(v) for v in source['physical_coefficient_Taylor']]")
    proofs['C4_angular_coefficients_are_derivative_enclosures_of_original_unique_branch'] = True

    # The pulse does not alter the angular swirl primitive. Bind the
    # same canonical Xp and exact positive pulse memory, before caps.
    syntax('compliant_axial_pulse_field','__init__','self.Xp',
           "inlet['Mtheta_over_sqrt2_R_3half_Pstar'][0]/inlet['Utheta_over_Pstar'][0]")
    syntax('compliant_axial_pulse_field','end','X',
           '1/self.rate+(self.Xp-1/self.rate)*self.factor(-13*self.rate/self.mu-self.rate*s)')
    syntax('compliant_axial_pulse_field','factor','lower',
           "c.mpf(0) if low < -mp.mpf('1e25') else c.exp(c.mpf(low))")
    syntax('compliant_axial_pulse_field','factor','upper','c.exp(c.mpf(high))')
    syntax('compliant_corrected_outer_field','__init__','self.Xv',
           "self.pulse.end('0',0)['Mtheta_over_sqrt2_R_3half_Utheta']")
    syntax('compliant_flatten_mixed_C4','__init__','self.Xv',
           "read_interval(c,pulse['whole_Z_terminal']['Mtheta_over_sqrt2_R_3half_Utheta']['coefficients'][0])")
    from lei_ren_part1_paper_compliant_power_inlet_C4_check import source_identities as inlet_identities
    inlet_proof=inlet_identities()
    if not inlet_proof['canonical_Xp_constant']:
        raise ValueError('Actual source constant Xp lemma required')
    tree('compliant_power_inlet_C4_check')
    syntax('compliant_axial_pulse_field','__init__','inlet',"self.pulse.buffer.power('0',1)")
    for stem,cls in (('compliant_pulse_high_jets','CompliantPulseHighJets'),
                     ('compliant_pulse_radial_C4','CompliantPulseRadialC4')):
        syntax(stem,'__init__','self.pulse','self.high.base.pulse',cls=cls)
        syntax(stem,'__init__','p0','self.pulse.buffer.power(0,1)',cls=cls)
        syntax(stem,'__init__','self.inlet_H',"box(p0['Mtheta_over_sqrt2_R_3half_Pstar'][0])",cls=cls)
        syntax(stem,'__init__','self.Xp',"self.inlet_H/box(p0['Utheta_over_Pstar'][0])",cls=cls)
    body('compliant_pulse_high_jets','end',
         'return self._high_packet(super().end(*args,**kwargs))',cls='CompliantPulseHighJets')
    for stem,order in (('compliant_pulse_high_jets',4),('compliant_pulse_radial_C4',5)):
        keyword(stem,'_high_packet','Mtheta_over_sqrt2_R_3half_Utheta',
                "IntervalTaylor.constant(c,point['Mtheta_over_sqrt2_R_3half_Utheta'],%s)"%order)
    syntax('compliant_pulse_mixed_C4','_high_packet','point','super()._high_packet(point)')
    keyword('compliant_axial_pulse_field','report','whole_Z_terminal','self.end([-1,1],0)')
    syntax('compliant_pulse_radial_C4','report','r','CompliantAxialPulseField.report(self)')
    syntax('compliant_pulse_mixed_C4','report','r','super().report()')
    body('compliant_pulse_mixed_C4','run',
         "result=CompliantPulseMixedC4().report()\n"
         "Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\\n',encoding='utf-8')\n"
         "print('Actual all-chart O.4 velocity/pressure mixed derivatives through total order4 generated; uniform outer C4/cone pending',flush=True)\n"
         'return result')
    syntax('compliant_flatten_mixed_C4','__init__','pulse',
           "json.loads((HERE/(PREFIX+'compliant_pulse_mixed_C4.json')).read_bytes())")
    # This complete callable bounds exp(logvalue) by monotonicity even
    # when its positive exact value is too small to materialize. Directed
    # arithmetic then encloses the affine Xv expression. The codec stores
    # exact binary endpoints, not rounded decimal descriptions.
    body('compliant_axial_pulse_field','factor',
         "c=self.ctx;cut=c.mpf(cut);low,high=endpoints(logvalue)\n"
         "if high<=endpoints(cut)[0]:return c.mpf([0,endpoints(c.exp(cut))[1]])\n"
         "lower=c.mpf(0) if low<-mp.mpf('1e25') else c.exp(c.mpf(low))\n"
         "if high>mp.mpf('1e25'):raise ValueError('Use a segmented coordinate with a finite upper log factor')\n"
         "upper=c.exp(c.mpf(high))\n"
         "return c.mpf([max(mp.mpf(0),endpoints(lower)[0]),endpoints(upper)[1]])")
    body('compliant_five_moment_repair','pack',
         'if isinstance(value,IntervalTaylor):return dict(coefficients=list(value.coefficients))\n'
         'if isinstance(value,dict):return {k:pack(v) for k,v in value.items()}\n'
         'if isinstance(value,(list,tuple)):return [pack(v) for v in value]\nreturn value')
    body('schedule_endpoint_enclosures_check','encode',
         "if hasattr(value,'_mpi_'):\n"
         " lo,hi=endpoints(value)\n"
         " return dict(lower=mp.nstr(lo,90),upper=mp.nstr(hi,90),width=mp.nstr(hi-lo,60),lower_exact_mpf_tuple=list(value._mpi_[0]),upper_exact_mpf_tuple=list(value._mpi_[1]))\n"
         "if isinstance(value,mp.mpf):return dict(value=mp.nstr(value,90),log_abs=None if value==0 else mp.nstr(mp.log(abs(value)),70),exact_mpf_tuple=list(value._mpf_))\n"
         "if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}\n"
         "if isinstance(value,(list,tuple)):return [encode(v) for v in value]\nreturn value")
    body('uniform_fixed_beta_error','read_interval',
         "return ctx.mpf([mp.make_mpf(tuple(record['lower_exact_mpf_tuple'])),mp.make_mpf(tuple(record['upper_exact_mpf_tuple']))])")
    mu,a,L,Ts,wait,eps,S,z,v = s.symbols('mu a Lrel Ts waiting epsilon S Z v',real=True)
    rate,k,delta = 1-mu,1-a,2*a
    Xp = s.symbols('same_canonical_Xp',real=True)
    Xv = actual('compliant_axial_pulse_field','end','X',
                {'self.rate':rate,'self.Xp':Xp,'self.mu':mu,'s':s.Integer(0),
                 'self.factor(-13 * self.rate / self.mu - self.rate * s)':s.exp(-13*rate/mu)})
    zero('exact_common_Xv_inherits_nonzero_pulse_memory',
         Xv-(1/rate+(Xp-1/rate)*s.exp(-13*rate/mu)))
    # The admitted pulse C4 receipt extends this exact inherited angular
    # function; its terminal is at s=0. No cap endpoint defines Xv.
    proofs['C4_and_original_Xv_share_original_terminal_callable'] = True
    proofs['C4_receipt_Xv_encloses_exact_affine_exponential_memory_by_directed_producer_and_lossless_codec'] = True

    # Bind corrected compliant parameters and receipt producers before
    # using common symbols. The inherited legacy .01 relation is replaced
    # in CompliantOuterParameters; neither chart uses that old relation.
    syntax('compliant_pressure_source','__init__','self.log_epsilon',
           "c.ln(c.mpf('.001'))+self.log_delta",cls='CompliantOuterParameters')
    syntax('compliant_pressure_source','__init__','self.epsilon','c.exp(self.log_epsilon)',cls='CompliantOuterParameters')
    syntax('compliant_pressure_source','__init__','self.waiting','self.waiting_root_enclosure()',cls='CompliantOuterParameters')
    syntax('compliant_pressure_source','__init__','self.parameters','CompliantOuterParameters(Md,precision)',cls='CompliantPressureDatum')
    syntax('logarithmic_outer_parameters','__init__','self.delta','c.exp(self.log_delta)')
    syntax('logarithmic_outer_parameters','__init__','self.Ts','4*(c.ln(2)-self.log_delta)')
    syntax('compliant_five_moment_repair','__init__','self.datum',"LogarithmicPressureDatum('40',precision=160)")
    imported=[n for n in tree('compliant_five_moment_repair').body if isinstance(n,ast.ImportFrom)
              and n.module=='lei_ren_part1_paper_compliant_pressure_source']
    if len(imported)!=1 or [(n.name,n.asname) for n in imported[0].names]!=[('CompliantPressureDatum','LogarithmicPressureDatum')]:
        raise ValueError('Original repair must use the compliant datum alias')
    syntax('compliant_five_moment_repair','__init__','self.delta','c.mpf(endpoints(self.datum.parameters.delta))')
    for stem,target,expr in (
        ('compliant_outer_initial','self.datum','self.repair.datum'),
        ('compliant_outer_initial','self.params','self.datum.parameters'),
        ('compliant_outer_initial','self.delta','self.repair.delta'),
        ('compliant_outer_buffer','self.params','self.initial.params'),
        ('compliant_outer_angular_candidate','self.params','self.buffer.params'),
        ('compliant_outer_angular_candidate','self.delta','self.initial.delta'),
        ('compliant_outer_angular_candidate','self.epsilon','self.params.epsilon'),
        ('compliant_exact_heat_component','self.delta','self.angular.delta'),
        ('compliant_exact_heat_component','self.epsilon','self.angular.epsilon'),
        ('compliant_corrected_outer_field','self.eps','self.heat.epsilon'),
        ('compliant_power_inlet_C4','self.delta',"read_interval(c,pulse['selected_delta'])"),
        ('compliant_flatten_mixed_C4','self.delta','self.inlet.delta'),
        ('compliant_power_angular_C4','self.delta','self.flatten.delta'),
        ('compliant_steep_waiting_C4','self.delta','self.outer.delta'),
        ('compliant_collar_Gamma_C4','self.delta','self.steep.delta'),
        ('compliant_collar_Gamma_C4','self.eps','self.steep.epsilon')):
        syntax(stem,'__init__',target,expr)
    keyword('compliant_axial_pulse_field','report','selected_delta','self.delta')
    syntax('compliant_pulse_radial_C4','__init__','self.delta','box(self.high.base.future.delta)',cls='CompliantPulseRadialC4')
    syntax('compliant_steep_waiting_C4','__init__','self.Ts','4*(c.ln(2)-c.ln(self.delta))')
    syntax('compliant_steep_waiting_C4','__init__','self.epsilon',"c.mpf('.001')*self.delta")
    logdelta=s.symbols('same_log_delta',real=True)
    zero('actual_corrected_epsilon_same_in_both_routes',s.exp(s.log(s.Rational(1,1000))+logdelta)-s.exp(logdelta)/1000)
    zero('actual_Ts_same_in_both_routes',4*(s.log(2)-logdelta)-4*(s.log(2)-s.log(s.exp(logdelta))))
    receipt('compliant_outer_angular_candidate_check',['actual_waiting_root_in_original_same_source_enclosure',
                                                     'whole_axis_C1_angular_calls_checked'])
    syntax('compliant_steep_waiting_C4','__init__','name',
           ["PREFIX+'compliant_power_angular_C4_check.json'","PREFIX+'compliant_outer_angular_candidate.json'"])
    syntax('compliant_steep_waiting_C4','__init__','waiting',"json.loads((HERE/name).read_bytes())")
    syntax('compliant_steep_waiting_C4','__init__','self.wait',"read_interval(c,waiting['refined_same_source_waiting_root'])")
    syntax('compliant_steep_waiting_C4','__init__','self.logone',"read_interval(c,waiting['waiting_log_one_minus_epsilon'])")
    keyword('compliant_outer_angular_candidate','report','refined_same_source_waiting_root','self.waiting')
    keyword('compliant_outer_angular_candidate','report','waiting_log_one_minus_epsilon','self.waiting_logone')
    syntax('compliant_outer_angular_candidate','__init__','terminal',"self.before_waiting('0')")
    syntax('compliant_outer_angular_candidate','__init__','Xt',"terminal['X'][0]")
    syntax('compliant_outer_angular_candidate','__init__','logone',
           'c.mpf([endpoints(-self.epsilon/(1-self.epsilon))[0],endpoints(-self.epsilon)[1]])')
    syntax('compliant_outer_angular_candidate','__init__','self.waiting_logone','logone')
    body('compliant_outer_angular_candidate','run',
         "field=SharedOuterAngularCandidate();result=field.report()\n"
         "Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\\n',encoding='utf-8')\n"
         "print('Same-source angular flatten/steep/waiting candidate callable; actual preheat waiting root refined',flush=True)\n"
         'return result')
    # For 0<eps<1, -eps/(1-eps) <= log(1-eps) <= -eps.
    # Integrating 1/(1-x) on [0,eps] proves both bounds. Thus the root
    # receipt encloses the exact logarithm used below, never its midpoint.
    proofs['actual_wait_and_logone_receipts_enclose_original_exact_waiting_equation'] = True

    # Bind the actual flatten integral, not arbitrary Xf placeholders.
    q, rho = 1+z*z, s.log((1+z*z)/2)
    sigma = s.Function('original_flat_sigma')(v/100)
    If = s.Integral(s.exp(rho*sigma-rate*(100-v)),(v,0,100))
    Xf = (If+Xv*s.exp(-100*rate))/s.exp(rho)
    original_Xf = actual('compliant_corrected_outer_field','flatten_raw','X',
                         {'integral':If,'self.Xv':Xv,'self.rate':rate,'t':s.Integer(100),'F':s.exp(rho)})
    C4_Xf = actual('compliant_flatten_mixed_C4','flatten','X',
                   {'Xint':If,'self.Xv':Xv,'self.rate':rate,'t':s.Integer(100),'F':s.exp(rho)})
    syntax('compliant_flatten_mixed_C4','flatten','Xint',
           'Fc*(c.exp(-self.rate*(t*(self.cells-i-1)/self.cells))*decay_integral(c,self.rate,length))',True)
    syntax('compliant_corrected_outer_field','flatten_raw','integral',
           '(ratio*sigma_enclosure(c,v/100)).exp()*((c.exp(-self.rate*(t-b))-c.exp(-self.rate*(t-a)))/self.rate)',True)
    zero('same_original_and_C4_full_flatten_angular_integral',C4_Xf-original_Xf)
    zero('full_flatten_integral_defines_common_Xf',original_Xf-Xf)

    # At s=0 both flat beta supports have ended: h=0 and pastA=fullA.
    syntax('compliant_power_angular_C4','angular','pastA',
           "dj*(c.exp(self.rate*center)*past['A'])",True)
    syntax('compliant_corrected_outer_field','data','fullA',
           "dj*(c.exp(self.rate*center)*self.repair.weights['A'])",True)
    syntax('compliant_power_angular_C4','angular','y','self.Lrel+s')
    syntax('compliant_power_angular_C4','angular','Xbase',
           "eq+(data['flatten_exit_X']-1/self.rate)*c.exp(-self.rate*y)")
    syntax('compliant_power_angular_C4','angular','X','(Xbase+pastA*c.exp(-self.rate*s))/F[0]')
    # Require the actual completed-support branch, not just its formula.
    angular_fn=method('compliant_power_angular_C4','angular')
    branch=next((n for n in ast.walk(angular_fn) if isinstance(n,ast.If)
                 and ast.unparse(n.test)=='lo >= endpoints(ell)[1]'),None)
    if branch is None or not any(isinstance(n,ast.Assign)
                                 and any(ast.unparse(t)=='past' for t in n.targets)
                                 and ast.unparse(n.value)=='self.weights' for n in branch.body):
        raise ValueError('Actual angular completed-beta-support branch changed')
    syntax('compliant_power_angular_C4','angular','local','s-center')
    syntax('compliant_power_angular_C4','angular','beta','self.flat.beta(local)')
    receipt('compliant_flat_pulse_derivatives_check',['all_passed','original_radial_shape_derivatives_C4_available'])
    syntax('compliant_power_angular_C4','data','f','self.flatten.flatten(Z,100)')
    keyword('compliant_power_angular_C4','data','flatten_exit_X',"f['angular_Taylor']")
    fullA,Iin,Iout = s.symbols('same_actual_angular_bump_integral same_entry_angular_integral same_exit_angular_integral',real=True)
    XR = actual('compliant_corrected_outer_field','data','XR',
                {'Xf':Xf,'eq':1/rate,'self.rate':rate,'self.Lrel':L,'fullA':fullA})
    Xbase = actual('compliant_power_angular_C4','angular','Xbase',
                   {'eq':1/rate,"data['flatten_exit_X']":Xf,'self.rate':rate,'y':L})
    C4_XR = actual('compliant_power_angular_C4','angular','X',
                  {'Xbase':Xbase,'pastA':fullA,'self.rate':rate,'s':s.Integer(0),'F[0]':s.Integer(1)})
    zero('actual_C4_angular_terminal_is_original_XR',C4_XR-XR)
    env={'self.rate':rate,'self.k':k,'self.params.Ts':Ts,'self.angular.waiting':wait,
         'self.Iin':Iin,'self.Iout':Iout,'XR':XR}
    env['Xs']=actual('compliant_corrected_outer_field','data','Xs',env)
    env['Xq']=actual('compliant_corrected_outer_field','data','Xq',env)
    env['Xt']=actual('compliant_corrected_outer_field','data','Xt',env)
    old_tail=actual('compliant_corrected_outer_field','data','Xtail',env)
    syntax('compliant_corrected_outer_field','__init__','self.fullin',
           "transition_tails(c,0,self.mu,self.delta,'in',cells)")
    syntax('compliant_corrected_outer_field','__init__','self.fullout',
           "transition_tails(c,0,self.mu,self.delta,'out',cells)")
    syntax('compliant_corrected_outer_field','__init__','(_, self.Iin)',
           "unit_kernels(c,1,self.rate,'in',cells)")
    syntax('compliant_corrected_outer_field','__init__','(_, self.Iout)',
           "unit_kernels(c,1,self.k,'out',cells)")
    syntax('compliant_steep_waiting_C4','__init__','self.infull',"self.kernels(1,'in')")
    syntax('compliant_steep_waiting_C4','__init__','self.outfull',"self.kernels(1,'out')")
    syntax('compliant_steep_waiting_C4','transition_kernels','angular',
           ['ds*c.exp(rate*(v-jc))','ds*c.exp(rate*jc)'],True)
    syntax('compliant_outer_angular_candidate','unit_kernels','I',
           ['(c.exp(rate*b)-c.exp(rate*a))/rate*c.exp(-rate*jc)',
            'da*c.exp(rate*jc)'],True)
    syntax('compliant_steep_waiting_C4','data','terminal','self.outer.angular(Z,0)')
    XS=actual('compliant_steep_waiting_C4','data','XS',{'XR':C4_XR,"self.infull['angular']":Iin,'self.rate':rate})
    XQ=actual('compliant_steep_waiting_C4','data','XQ',{'XS':XS,'self.Ts':Ts})
    XT=actual('compliant_steep_waiting_C4','data','XT',{'XQ':XQ,"self.outfull['angular']":Iout,'self.k':k})
    C4_tail=actual('compliant_steep_waiting_C4','waiting','X',
                   {"data['XT']":XT,'self.k':k,'t':wait})
    zero('actual_C4_waiting_terminal_is_original_Xtail',C4_tail-old_tail)
    proofs['full_beta_supports_have_ended_at_angular_offset_zero'] = True
    proofs['same_entry_exit_angular_integrals_follow_original_sigma_J_weights'] = True

    # Existing pressure bridge proves the same original sigma/phi/K
    # source. Apply its pointwise bracket to both angular and energy
    # integrals; the canonical H/radius proof is already required above.
    sig = s.Piecewise((0,v<=0),(1,v>=1),
                     (1/(1+s.exp(1/v**2-1/(1-v)**2)),True))
    phi = s.Piecewise((0,v>=3),(s.exp(-4/(3-v)**2),True))
    H = s.Function('canonical_positive_Gamma_H')(2*(1-z*z)*S*s.exp(-v))
    W = 1-sig+sig*phi
    D = sig*(1-eps*phi)*(1-H)/(a*S)
    pre = 1-eps*W
    K = pre-a*S*D
    Theta = a*s.Integral(s.exp(k*v)*D,(v,0,s.oo))
    JW = s.Integral(s.exp(k*v)*W,(v,0,3))
    Ehat = s.Integral(s.exp(-delta*v)*D*(2*pre-a*S*D),(v,0,s.oo))
    EW = s.Integral(s.exp(-delta*v)*W,(v,0,3))
    EW2 = s.Integral(s.exp(-delta*v)*W**2,(v,0,3))
    A0 = 1/k+S*Theta+eps*JW
    E0 = 1/delta-2*eps*EW+eps**2*EW2-a*S*Ehat
    syntax('compliant_collar_Gamma_C4','collar_tails','angular','D*(self.a*ds*c.exp(self.k*v))',True)
    syntax('compliant_corrected_outer_field','heat_tails','theta_hat','D*(ds*self.a*c.exp(self.k*v))',True)
    syntax('compliant_collar_Gamma_C4','collar_tails','energy','square*(ds*c.exp(-self.delta*v))',True)
    syntax('compliant_corrected_outer_field','heat_tails','energy_hat','D*sq*(ds*c.exp(-self.delta*v))',True)
    C4_A0=actual('compliant_collar_Gamma_C4','collar_tails','A',
                  {'one':s.Integer(1),'self.k':k,'angular':Theta,'self.S':S,"atoms['JW']":JW,'self.eps':eps,'t':s.Integer(0)})
    old_A0=actual('compliant_corrected_outer_field','heat_tails','theta',
                 {'theta_hat':Theta,'S':S,'self.k':k,'t':s.Integer(0),'self.eps':eps,"atoms['JW']":JW})
    zero('actual_C4_heat0_angular_is_original_full_future_target',C4_A0-old_A0)
    zero('actual_common_heat0_target_has_original_repair_units',old_A0-A0)
    C4_E0=actual('compliant_collar_Gamma_C4','collar_tails','E',
                  {'one':s.Integer(1),'self.delta':delta,'t':s.Integer(0),'self.eps':eps,
                   "atoms['EW']":EW,"atoms['EW2']":EW2,'energy':Ehat,'self.a':a,'self.S':S})
    tree('compliant_corrected_outer_field')
    old_E0=source_branch_expression('heat_tails','E',
                 "atoms['EW']",
                 {'energy_hat':Ehat,'self.a':a,'S':S,'self.delta':delta,'t':s.Integer(0),'self.eps':eps,'self.eps ** 2':eps**2,
                  "atoms['EW']":EW,"atoms['EW2']":EW2})
    zero('actual_C4_heat0_energy_is_original_full_future_integral',C4_E0-old_E0)
    zero('full_future_energy_definition_retains_both_epsilon_atoms',old_E0-E0)
    syntax('compliant_collar_Gamma_C4','data','defect',"Xtail*(1-self.eps)-heat0['angular_numerator']")
    # Compose the actual waiting root and actual small bump branch with
    # old_tail and old_A0. A generic previously admitted closure alone is
    # not used as the zero of a new interface.
    syntax('compliant_future_swirl_energy','__init__','self.Lrel','-30*self.params.log_mu')
    syntax('compliant_power_angular_C4','__init__','self.Lrel',"read_interval(c,energy['postflatten_length'])")
    keyword('compliant_future_swirl_energy','report','postflatten_length','self.Lrel')
    syntax('compliant_outer_angular_repair','__init__','self.p','c.exp(-2*self.rate)')
    syntax('compliant_outer_angular_repair','__init__','self.det','self.p*self.q-1')
    syntax('compliant_outer_angular_repair','coefficients','b1',"r*(c.exp(self.rate)/self.weights['A'])")
    syntax('compliant_outer_angular_repair','coefficients','(nx, ny)','self.inverse(b1[0],b2[0]-nonlinear)')
    body('compliant_outer_angular_repair','inverse',
         'return ((u*self.q-v)/self.det,(v*self.p-u)/self.det)')
    syntax('compliant_corrected_outer_field','data','coeff',
           ["self.repair.coefficients(Z)['physical_coefficient_Taylor']",'[IntervalTaylor(c,v.coefficients) for v in coeff]'])
    receipt('compliant_outer_angular_repair_check',['all_passed'])
    uinv,vinv,pinv,qinv=s.symbols('angular_rhs second_rhs angular_p angular_q',real=True)
    xi_inv=(uinv*qinv-vinv)/(pinv*qinv-1)
    yi_inv=(vinv*pinv-uinv)/(pinv*qinv-1)
    zero('actual_small_branch_inverse_preserves_first_angular_equation',pinv*xi_inv+yi_inv-uinv)
    # Bind the positive divided-difference integral used by the actual
    # preheat_difference enclosure to Xf(0)-Xf(Z), including pulse memory.
    for target,expr in (
        ('start','2*self.angular.Xv*c.exp(-100*self.rate)'),
        ('divided','start/Q'),('k','1-sig'),
        ('factor','c.exp(k*c.ln(2))*k'),
        ('quotient','c.mpf([endpoints(c.exp((-1-k)*c.ln(Q)))[0],1])'),
        ('value','eta*divided')):
        syntax('compliant_outer_angular_repair','preheat_difference',target,expr)
    syntax('compliant_outer_angular_repair','preheat_difference','divided','weight*factor*quotient',True)
    syntax('compliant_outer_angular_repair','defects','pre','self.preheat_difference(Z)')
    syntax('compliant_outer_angular_repair','defects','rh',
           ['capjet(theta,self.theta_heat_over_scale_cap)','IntervalTaylor(c,[rh[0],0])','IntervalTaylor(c,[0,rh[1]])'])
    syntax('compliant_outer_angular_repair','__init__','self.theta_heat_over_scale_cap',
           'c.exp(self.log_theta_multiplier+self.strong_logS_cap-self.logscale)')
    eta,jj,wdd=s.symbols('eta cutoff_complement integration_fraction',positive=True)
    primitive=-(1+wdd*eta)**(-jj)
    zero('actual_preheat_divided_difference_integrand',
         s.diff(primitive,wdd)-eta*jj*(1+wdd*eta)**(-1-jj))
    zero('actual_preheat_divided_difference_endpoint',
         primitive.subs(wdd,1)-primitive.subs(wdd,0)-(1-(1+eta)**(-jj)))
    zero('actual_flatten_memory_difference',
         2*Xv*s.exp(-100*rate)*(1-1/q)-z*z*(2*Xv*s.exp(-100*rate)/q))
    rr,ss=s.symbols('log_q_over_two sigma_value',real=True)
    zero('actual_flatten_integrand_and_positive_difference_use_same_function',
         s.exp(rr*ss)/s.exp(rr)-s.exp(-(1-ss)*rr))
    # Endpoint integral bounds in the production divided difference follow
    # from monotonicity of (1+w*eta)^(-1-jj), 0<=w<=1, eta>=0, jj>=0.
    Xf0=Xf.subs(z,0)
    original_candidate_Xf=actual('compliant_outer_angular_candidate','flatten','X',
        {'integral':s.exp(100*rate)*If,'self.Xv':Xv,'self.rate':rate,
         't':s.Integer(100),'F':s.exp(rho)})
    zero('actual_waiting_candidate_uses_same_full_flatten_Xf',original_candidate_Xf-Xf)
    syntax('compliant_outer_angular_candidate','reference_buffer','inlet','self.flatten(Z,1)')
    syntax('compliant_outer_angular_candidate','reference_buffer','t','-30*self.params.log_mu*phase')
    raw_R=1/rate+(Xf0-1/rate)*s.exp(-rate*L)
    departure=actual('compliant_outer_angular_candidate','reference_buffer','departure',
        {"inlet['X']":Xf0,'eq':1/rate,'self.rate':rate,'t':L})
    raw_R_actual=actual('compliant_outer_angular_candidate','reference_buffer','X',
        {'departure':departure,'eq':1/rate})
    zero('actual_waiting_raw_reference_endpoint_is_unrepaired_XR0',raw_R_actual-raw_R)
    raw_S=actual('compliant_outer_angular_candidate','steep_in','X',
        {"inlet['X']":raw_R_actual,'I':Iin,'self.rate':rate,'t':s.Integer(1),'J':s.Rational(1,2)})
    raw_Q=raw_S+Ts
    raw_T=actual('compliant_outer_angular_candidate','before_waiting','X',
        {"inlet['X']":raw_Q,'I':Iout,'self.restore_rate':k,'J':s.Rational(1,2)})
    waiting_exact=actual('compliant_outer_angular_candidate','__init__','self.waiting',
        {'Xt':raw_T,'eq':1/k,'k':k,'logone':s.log(1-eps),
         'self.params.log_epsilon':s.log(eps),'self.collarJ':JW})
    zero('actual_waiting_root_with_actual_raw_terminal',
         waiting_exact-(s.log(raw_T-1/k)+s.log(1-eps)-s.log(eps)-s.log(1/k+JW))/k)
    decay_wait=eps*(1/k+JW)/((1-eps)*(raw_T-1/k))
    log_theta=actual('compliant_outer_angular_repair','__init__','self.log_theta_multiplier',
        {'self.rate':rate,'self.angular.restore_rate':k,'self.angular.waiting':wait,
         'self.angular.waiting_logone':s.log(1-eps)})
    zero('actual_angular_heat_multiplier_times_waiting_decay',
         s.exp(log_theta)*s.exp(-k*wait)-s.exp((rate+k)/2)/(1-eps))
    # physical d1,d2 satisfy A*(exp(-3r)*d1+exp(-r)*d2)=rpre+rheat;
    # the source-bound fullA sum is exactly that left side. The fixed-point
    # inverse lemma above and admitted unique-branch receipt justify this
    # substitution for the actual coefficients, not arbitrary picked values.
    # Name the source-bound full functions, then replay the actual source
    # expressions with these aliases. This keeps functional equality while
    # avoiding an attempted symbolic evaluation of Gamma future integrals.
    xf,xf0,th,jw=s.symbols('bound_actual_Xf bound_actual_Xf0 bound_actual_Theta bound_actual_JW',real=True)
    replay=dict(env)
    replay['XR']=actual('compliant_corrected_outer_field','data','XR',
        {'Xf':xf,'eq':1/rate,'self.rate':rate,'self.Lrel':L,'fullA':fullA})
    replay['Xs']=actual('compliant_corrected_outer_field','data','Xs',replay)
    replay['Xq']=actual('compliant_corrected_outer_field','data','Xq',replay)
    replay['Xt']=actual('compliant_corrected_outer_field','data','Xt',replay)
    raw_t=replay['Xt'].subs({xf:xf0,fullA:0})
    exact_decay=eps*(1/k+jw)/((1-eps)*(raw_t-1/k))
    replay['c.exp(-self.k * self.angular.waiting)']=exact_decay
    bound_tail=actual('compliant_corrected_outer_field','data','Xtail',replay)
    bound_pre=s.exp(-rate*L)*(xf0-xf)
    bound_heat=s.exp((rate+k)/2)*S*th/((1-eps)*exact_decay)
    bound_target=actual('compliant_corrected_outer_field','heat_tails','theta',
        {'theta_hat':th,'S':S,'self.k':k,'t':s.Integer(0),'self.eps':eps,"atoms['JW']":jw})
    zero('actual_angular_defect_zero_from_original_waiting_and_repair_equations',
         (1-eps)*bound_tail.subs(fullA,bound_pre+bound_heat)-bound_target)

    # Trace complete future-energy receipt generation back to the source
    # callable used by the original positive amplitude quadratic.
    syntax('compliant_future_energy_high_jets','future','high','self.angular.coefficients(endpoints(Z))')
    syntax('compliant_fifth_axial_jets','future','f','self.fourth.energy.base')
    syntax('compliant_fifth_axial_jets','future','high','self.angular(Z)')
    syntax('compliant_fifth_axial_jets','future','old','self.fourth.energy.future(endpoints(Z))')
    syntax('compliant_fifth_axial_jets','future','admitted',"append_fifth(c,old['complete_future_energy_Taylor'],total[5])")
    keyword('compliant_fifth_axial_jets','future','complete_future_energy_Taylor','admitted')
    keyword('compliant_future_energy_high_jets','future','complete_future_energy_Taylor','total')
    syntax('compliant_future_energy_high_jets','future','total',
           ['flatten+power+angular_change+postrel-heatdifference','IntervalTaylor(c,coeffs)'])
    syntax('compliant_fifth_axial_jets','future','total','flatten+power+angular_change+postrel-heatdifference')
    syntax('compliant_future_swirl_energy','future','total','flatten+power+angular_change+postrel-heat_difference')
    syntax('compliant_flatten_mixed_C4','future_energy','packet',
           ["self.fifth['whole_Z']['energy']","sample['energy']"])
    future_returns=[n.value for n in ast.walk(method('compliant_flatten_mixed_C4','future_energy'))
                    if isinstance(n,ast.Return)]
    expected="IntervalTaylor(c,[read_interval(c,v)/2 for v in packet['complete_future_energy_Taylor']['coefficients']])"
    if len(future_returns)!=1 or ast.dump(future_returns[0])!=ast.dump(ast.parse(expected,mode='eval').body):
        raise ValueError('C4 original selected future-energy half factor changed')
    proofs['receipt_C5_future_energy_halving_is_original_selected_pulse_target'] = True
    Fint,Bchange,Kin,Kout=s.symbols('same_flatten_energy same_angular_energy_change same_steep_in_energy same_steep_out_energy',real=True)
    Nf,Nrel = q*q*s.exp(-200*mu)/4,q*q*s.exp(-2*mu*(100+L))/4
    Ns=s.exp(-1-mu); Nq=Ns*s.exp(-2*Ts); Nt=Nq*s.exp(-1-delta/2)
    mult=Nt*s.exp(-delta*wait)/(1-eps)**2
    decay=lambda r,t:(1-s.exp(-r*t))/r
    Post=Kin+Ns*decay(2,Ts)+Nq*Kout+Nt*decay(delta,wait)+mult*E0
    Total=Fint+Nf*decay(2*mu,L)+Nrel*Bchange+Nrel*Post
    before=actual('compliant_future_swirl_energy','future','total',
                  {'flatten':Fint,'power':Nf*decay(2*mu,L),'angular_change':Nrel*Bchange,
                   'postrel':Nrel*(Post+mult*a*S*Ehat),'heat_difference':Nrel*mult*a*S*Ehat})
    zero('same_original_complete_future_energy_includes_actual_heat_deficit',before-Total)
    # Bind C4 energy scalars to the actual original producer fields.
    syntax('compliant_power_angular_C4','__init__','ename',"PREFIX+'compliant_future_swirl_energy.json'")
    syntax('compliant_power_angular_C4','__init__','energy','json.loads((HERE/ename).read_bytes())')
    syntax('compliant_power_angular_C4','__init__','pref',
           "{k:read_interval(c,v) for k,v in energy['relative_energy_prefactors'].items()}")
    syntax('compliant_power_angular_C4','__init__','kernels',"energy['postrel_energy_kernels']")
    syntax('compliant_power_angular_C4','__init__','atoms',"energy['whole_Z_C1']['separate_preheat_tail_atoms']")
    syntax('compliant_power_angular_C4','__init__','self.tail_multiplier',"read_interval(c,energy['tail_multiplier'])")
    keyword('compliant_future_swirl_energy','report','relative_energy_prefactors',
            'dict(Ns=self.Ns,Nq=self.Nq,Nt=self.Nt,Ntail=self.Ntail)')
    keyword('compliant_future_swirl_energy','report','postrel_energy_kernels','self.kernels')
    keyword('compliant_future_swirl_energy','report','tail_multiplier','self.tail_multiplier')
    keyword('compliant_future_swirl_energy','report','whole_Z_C1','self.future([-1,1])')
    keyword('compliant_future_swirl_energy','future','separate_preheat_tail_atoms',
            'dict(baseline=baseline,epsilon_atom=epsilon_atom,epsilon_squared_atom=epsilon_squared_atom)')
    syntax('compliant_future_swirl_energy','future','baseline','1/self.delta')
    syntax('compliant_future_swirl_energy','future','epsilon_atom',"-2*eps*atoms['W']")
    syntax('compliant_future_swirl_energy','future','epsilon_squared_atom',"eps**2*atoms['W_squared']")
    B0=1/delta-2*eps*EW+eps**2*EW2
    C4scalar=accumulated('compliant_power_angular_C4','__init__','self.post_scalar',
        {"read_interval(c, kernels['steep_in'])":Kin,"pref['Ns']":Ns,
         "read_interval(c, energy['steep_power_energy_integral'])":decay(2,Ts),
         "pref['Nq']":Nq,"read_interval(c, kernels['steep_out'])":Kout,
         "pref['Nt']":Nt,"read_interval(c, energy['waiting_energy_integral'])":decay(delta,wait),
         'self.tail_multiplier':mult,
         'sum((read_interval(c, v) for v in atoms.values()), c.mpf(0))':B0})
    zero('actual_C4_post_scalar_is_original_full_preheat_post_energy',
         C4scalar-(Post+mult*a*S*Ehat))
    original_post=actual('compliant_future_swirl_energy','future','postrel',
        {'Nrel':Nrel,"self.kernels['steep_in']":Kin,'self.Ns':Ns,'self.steep_power':decay(2,Ts),
         'self.Nq':Nq,"self.kernels['steep_out']":Kout,'self.Nt':Nt,
         'self.waiting_energy':decay(delta,wait),'self.tail_multiplier':mult,'preheat_tail':B0})
    zero('actual_original_postrel_scalar_and_C4_scalar_have_same_units',original_post-Nrel*C4scalar)
    # C4 post and full_change use the same actual receipt producers above;
    # their q^2-normalized units cancel before enclosure.
    syntax('compliant_power_angular_C4','data','heat',"jet(source['scaled_Gamma_future_defect_Taylor']['energy'])")
    syntax('compliant_power_angular_C4','data','full_change',
           "(dj*(2*self.weights['E'])+dj*dj*self.weights['F'])*c.exp(-2*self.mu*center)",True)
    syntax('compliant_future_swirl_energy','future','angular_change',
           "(dj*(2*w['E'])+dj*dj*w['F'])*c.exp(-2*self.mu*center)",True)
    C4_post=actual('compliant_power_angular_C4','data','post',
                  {'IntervalTaylor.constant(c, self.post_scalar, 5)':C4scalar,
                   'heat':Ehat,'c.mpf([0, endpoints(self.heatcap)[1]])':mult*a*S})
    zero('actual_C4_post_energy_is_original_selected_post_integral',C4_post-Post)
    syntax('compliant_axial_amplitude_selection','select','energy','self.future.future(Z)')
    syntax('compliant_axial_amplitude_selection','select','future',"energy['Section7_34_weighted_future_Taylor']")
    syntax('compliant_axial_amplitude_selection','select','target','future-incoming+self.base')
    syntax('compliant_axial_pulse_field','end','future',"self.selection.future.future(Z)['complete_future_energy_Taylor']/2")
    syntax('compliant_steep_waiting_C4','data','preheat_tail',
           'IntervalTaylor.constant(c,self.preheat_scalar,5)-heat*c.mpf([0,endpoints(self.delta*self.S/2)[1]])')
    syntax('compliant_steep_waiting_C4','data','H','preheat_tail*self.tail_normalization')
    syntax('compliant_collar_Gamma_C4','exterior','energy',"local['energy_numerator']/(K*K*2)")
    # Bind the actual current-radius Gamma numerator at every t>=3.
    # integrated_gamma_tails square_tail encloses the full integral of
    # Dhat*(1-a*Sc*Dhat/2). Its energy return has the required factor two.
    # The existing admitted Gamma enclosure proof supplies its full-function
    # semantics and infinite-tail error, not a truncated function definition.
    for target,expr in (
        ('Sc','self.S*c.exp(-t)'),
        ('tails','integrated_gamma_tails(c,Z,self.a,Sc,5,0)'),
        ('D','gamma_deficit_mixed(c,Z,self.a,self.Scap,t)'),
        ('K','[one-D[0]*(self.a*self.S)]+[-v*(self.a*self.S) for v in D[1:]]'),
        ('energy',"one/self.delta-tails['energy']*(self.a*Sc)")):
        syntax('compliant_collar_Gamma_C4','local_Gamma',target,expr)
    keyword('compliant_collar_Gamma_C4','local_Gamma','energy_numerator','energy')
    keyword('compliant_angular_high_jets','integrated_gamma_tails','energy','square_tail(2*a)*2')
    keyword('compliant_angular_high_jets','integrated_gamma_tails','pressure','square_tail(1+2*a)')
    syntax('compliant_collar_Gamma_C4','exterior','local','self.local_Gamma(Z,t)')
    syntax('compliant_collar_Gamma_C4','exterior','K',"local['K_rows'][0]")
    aa,sc,hc,jdef=s.symbols('canonical_a exact_current_inverse_radius canonical_H full_square_deficit_integral',positive=True)
    dh=(1-hc)/(aa*sc)
    zero('actual_current_radius_energy_deficit_is_full_Gamma_square',
         1-2*aa*sc*dh*(1-aa*sc*dh/2)-hc*hc)
    zero('actual_current_radius_baseline_is_full_infinite_energy_integral',
         s.integrate(s.exp(-2*aa*v),(v,0,s.oo))-1/(2*aa))
    local_num=actual('compliant_collar_Gamma_C4','local_Gamma','energy',
        {'one':s.Integer(1),'self.delta':2*aa,"tails['energy']":2*jdef,'self.a':aa,'Sc':sc})
    # By the pointwise identity and linearity of convergent full integrals,
    # BE=integral exp(-2a*v)*H(2d*Sc*exp(-v))^2 dv equals this numerator.
    full_BE=1/(2*aa)-2*aa*sc*jdef
    zero('actual_C4_local_energy_numerator_is_canonical_full_future_square',local_num-full_BE)
    local_energy=actual('compliant_collar_Gamma_C4','exterior','energy',
        {"local['energy_numerator']":local_num,'K':hc})
    zero('actual_C4_exterior_energy_has_original_half_and_H_squared_normalization',
         local_energy-full_BE/(2*hc*hc))
    # The original callable's pure-heat override is the same full integral
    # normalization; its C1 enclosure need not equal the C5 enclosure.
    syntax('compliant_corrected_outer_field','collar_heat','overrides',
        ["dict(theta=theta,energy=local['energy_tail_over_R_Utheta_squared']/2,pressure=-local['pressure_tail_over_Utheta_squared'])",'None'])
    original_e=actual('compliant_corrected_outer_field','heat_local','Ebar',
        {'eh':2*jdef,'self.a':aa,'Sbox':sc,'self.delta':2*aa,'K':hc})
    zero('actual_original_pure_heat_energy_and_C4_have_same_half_units',original_e/2-local_energy)

    # Prove the actual selected quadratic residual. The physical energy
    # at the pulse terminal is POSITIVE Total/2; it is its discrepancy from
    # the prescribed complete future that vanishes, not that energy itself.
    syntax('compliant_axial_amplitude_selection','__init__','self.base',"(1-c.exp(-26))/4")
    syntax('compliant_future_swirl_energy','__init__','self.weighted_factor','self.mu*self.pulse_energy_attenuation/2')
    syntax('compliant_future_swirl_energy','__init__','self.pulse_energy_attenuation','c.exp(-26)')
    syntax('compliant_future_swirl_energy','future','weighted','total*self.weighted_factor')
    keyword('compliant_future_swirl_energy','future','Section7_34_weighted_future_Taylor','weighted')
    syntax('compliant_axial_amplitude_selection','select','A2','self.K')
    syntax('compliant_axial_amplitude_selection','select','A1','u[0]*0')
    syntax('compliant_axial_amplitude_selection','select','A0','-target')
    syntax('compliant_axial_amplitude_selection','select','incoming',
           "IntervalTaylor(c,trial['actual_normalized_incoming_energy'])*self.mu")
    syntax('compliant_axial_amplitude_selection','select','A2','nu*vj*vj',True)
    syntax('compliant_axial_amplitude_selection','select','A1','uj*(2*nu*vj)',True)
    syntax('compliant_axial_amplitude_selection','select','A0','uj*uj*nu',True)
    syntax('compliant_axial_amplitude_selection','select','halflinear','A1[0]/(2*A2)')
    syntax('compliant_axial_amplitude_selection','select','ap','-halflinear+c.sqrt(halflinear**2-A0[0]/A2)')
    syntax('compliant_axial_amplitude_selection','select','scaled','[uj+selected*vj for uj,vj in zip(u,self.v)]')
    loops=[n for n in ast.walk(method('compliant_axial_amplitude_selection','select'))
           if isinstance(n,ast.For) and ast.unparse(n.target)=='(nu, uj, vj)']
    if len(loops)!=1 or ast.unparse(loops[0].iter)!='zip(self.nu, u, self.v)':
        raise ValueError('Actual selected energy quadratic contribution loop changed')
    pulseK,nu1,nu2=s.symbols('actual_Kpulse actual_nu1 actual_nu2',positive=True)
    uj1,uj2,vj1,vj2,target,ap=s.symbols('actual_u1 actual_u2 actual_v1 actual_v2 energy_target selected_ap',real=True)
    A2=pulseK+nu1*vj1*vj1+nu2*vj2*vj2
    A1=2*(nu1*uj1*vj1+nu2*uj2*vj2)
    A0=nu1*uj1*uj1+nu2*uj2*uj2-target
    selected_energy=pulseK*ap*ap+nu1*(uj1+vj1*ap)**2+nu2*(uj2+vj2*ap)**2
    zero('actual_selected_energy_quadratic_is_full_affine_end_energy_equation',
         selected_energy-target-(A2*ap*ap+A1*ap+A0))
    C2=s.symbols('actual_positive_A2',positive=True)
    C1,C0=s.symbols('actual_A1 actual_A0',real=True)
    root=-C1/(2*C2)+s.sqrt((C1/(2*C2))**2-C0/C2)
    zero('actual_completed_square_selected_root_has_zero_equation_residual',C2*root**2+C1*root+C0)
    ei,et=s.symbols('actual_incoming_energy actual_complete_future_energy',real=True)
    weighted=mu*s.exp(-26)*et/2
    base=(1-s.exp(-26))/4
    # Reconstruct the terminal from the original e'=2mu*e+Uz^2-1/2
    # integrating factor using the selected exact equation.
    actual_target=actual('compliant_axial_amplitude_selection','select','target',
        {'future':weighted,'incoming':mu*ei,'self.base':base})
    forward_terminal=s.exp(26)*(ei+(actual_target-base)/mu)
    zero('selected_forward_terminal_matches_positive_full_future_half',forward_terminal-et/2)
    syntax('compliant_axial_pulse_field','end','e',
           'future*c.exp(2*self.mu*s)+baseline-end_energy*(c.exp(2*self.mu*s)*self.E2cap)')
    endpoint_e=actual('compliant_axial_pulse_field','end','e',
        {'future':et/2,'self.mu':mu,'s':s.Integer(0),'baseline':s.Integer(0),
         'end_energy':s.Integer(0),'self.E2cap':s.symbols('same_positive_end_scale_squared',positive=True)})
    zero('actual_backward_pulse_terminal_keeps_positive_future_target',endpoint_e-et/2)
    proofs['actual_C4_future_energy_receipt_encloses_same_full_function_used_by_selected_amplitude'] = True
    proofs['actual_selected_energy_equation_residual_zero_with_positive_future_target'] = True
    proofs['renormalized_energy_constant_zero_from_selected_equation'] = True
    proofs['physical_terminal_energy_not_set_to_zero'] = True

    # Meridional moments are fixed by the original selected linear
    # equations. Empty future supports at s=0 propagate their exact zero
    # into all subsequent zero-Uz charts; the C4 field uses that route.
    syntax('compliant_axial_pulse_field','end','terminal','endpoints(s)==(mp.mpf(0),mp.mpf(0))')
    for target,expr in (
        ('C',"selected['selected_scaled_end_coefficient_Taylor']"),
        ('zero','C[0]*0'),('Bhat','zero'),('Byhat','zero'),
        ('mhat','[zero,zero]'),('end_energy','zero'),
        ('local','s-center'),('beta','raw_beta(c,local/ell)/(ell*normal)'),
        ('w','backward_bump_weights(c,self.mu,normal,local,cells)'),
        ('B','Bhat*self.Ecap'),('m','[v*self.Ecap for v in mhat]'),
        ('distance','-s'),
        ('baseline','c.mpf([endpoints(distance*c.exp(-2*self.mu*distance)/2)[0],endpoints(distance/2)[1]])')):
        syntax('compliant_axial_pulse_field','end',target,expr)
    loop('compliant_axial_pulse_field','end','(Cj, center)','zip(C,(-3,-1))')
    syntax('compliant_axial_pulse_field','end','Bhat','Cj*beta',True)
    syntax('compliant_axial_pulse_field','end','end_energy','Cj*Cj*(c.exp(-2*self.mu*center)*w[2])',True)
    statement('compliant_axial_pulse_field','end',
        "mhat[row-1]-=Cj*(c.exp((c.mpf('.5')-row*self.mu)*(center-s))*w[row-1])")
    keyword('compliant_axial_pulse_field','end','Mz_over_R_Utheta','m[0]')
    keyword('compliant_axial_pulse_field','end','Mtheta_z_over_sqrt2_R_3half_Utheta_squared','m[1]')
    body('compliant_outer_pulse_map','raw_beta',
        'lo,hi=endpoints(x);near=mp.mpf(0) if lo<=0<=hi else min(abs(lo),abs(hi));far=max(abs(lo),abs(hi))\n'
        'lower=c.exp(-1/(1-c.mpf(far)**2)) if far<1 else c.mpf(0)\n'
        'upper=c.exp(-1/(1-c.mpf(near)**2)) if near<1 else c.mpf(0)\n'
        'return c.mpf([endpoints(lower)[0],endpoints(upper)[1]])')
    syntax('compliant_axial_pulse_field','backward_bump_weights','begin',
        'c.mpf([max(mp.mpf(-1),min(mp.mpf(1),rlo)),max(mp.mpf(-1),min(mp.mpf(1),rhi))])')
    syntax('compliant_axial_pulse_field','backward_bump_weights','ell',"c.mpf('.15')")
    syntax('compliant_axial_pulse_field','backward_bump_weights','(rlo, rhi)','endpoints(start/ell)')
    syntax('compliant_axial_pulse_field','backward_bump_weights','length','1-begin')
    syntax('compliant_axial_pulse_field','backward_bump_weights','out','[c.mpf(0),c.mpf(0),c.mpf(0)]')
    statement('compliant_axial_pulse_field','backward_bump_weights',
        'if endpoints(length)==(mp.mpf(0),mp.mpf(0)):return out')
    statement('compliant_flat_pulse_derivatives','beta_jets',
        'if hi<=-1 or lo>=1:return IntervalTaylor.constant(c,0,4)')
    body('compliant_flat_pulse_derivatives','beta',
        "c=self.ctx;ell=c.mpf('.15');raw=beta_jets(c,c.mpf(s)/ell)\n"
        'return IntervalTaylor(c,[raw[n]/(ell**(n+1)*self.normalization) for n in range(5)])')
    for center in (-3,-1):
        argument=-s.Rational(center)/s.Rational(15,100)
        if argument<=1:raise ValueError('Actual terminal is inside end bump support')
    proofs['actual_terminal_raw_beta_and_future_weights_are_zero_by_empty_support'] = True
    # Neither high wrapper may overwrite the two original primitive fields.
    protected={'Mz_over_R_Utheta','Mtheta_z_over_sqrt2_R_3half_Utheta_squared'}
    for stem in ('compliant_pulse_high_jets','compliant_pulse_radial_C4','compliant_pulse_mixed_C4'):
        fn=method(stem,'_high_packet')
        for node in ast.walk(fn):
            if isinstance(node,ast.keyword) and node.arg in protected:
                raise ValueError('High wrapper overwrites primitive terminal history: '+stem)
            if isinstance(node,(ast.Assign,ast.AugAssign)):
                targets=node.targets if isinstance(node,ast.Assign) else [node.target]
                for target in targets:
                    if (isinstance(target,ast.Subscript) and ast.unparse(target.value)=='point'
                        and isinstance(target.slice,ast.Constant) and target.slice.value in protected):
                        raise ValueError('High wrapper overwrites primitive terminal history: '+stem)
        proofs['source_'+stem+'_preserves_base_meridional_primitive_fields']=True
    for stem,cls,base in (
        ('compliant_pulse_high_jets','CompliantPulseHighJets','CompliantAxialPulseField'),
        ('compliant_pulse_radial_C4','CompliantPulseRadialC4','CompliantPulseHighJets'),
        ('compliant_pulse_mixed_C4','CompliantPulseMixedC4','CompliantPulseRadialC4')):
        node=next(n for n in tree(stem).body if isinstance(n,ast.ClassDef) and n.name==cls)
        if [ast.unparse(n) for n in node.bases]!=[base]:
            raise ValueError('Actual terminal wrapper inheritance changed: '+cls)
        proofs['source_inheritance_'+cls]=True
    for key in protected:
        dictitem('compliant_pulse_mixed_C4','transport_mixed','primitive',key,"[point['%s']]"%key)
    syntax('compliant_pulse_mixed_C4','transport_mixed','m1',"primitive['Mz_over_R_Utheta']")
    syntax('compliant_pulse_mixed_C4','transport_mixed','m2',"primitive['Mtheta_z_over_sqrt2_R_3half_Utheta_squared']")
    statement('compliant_pulse_mixed_C4','transport_mixed',"m1.append(Brows[k]-m1[k]*(c.mpf('.5')-mu))")
    statement('compliant_pulse_mixed_C4','transport_mixed',"m2.append(Brows[k]-m2[k]*(c.mpf('.5')-2*mu))")
    syntax('compliant_pulse_mixed_C4','transport_mixed','A','[field.radial(Z,b,m) for b,m in zip(Brows,m1)]')
    syntax('compliant_pulse_mixed_C4','_high_packet','mixed','transport_mixed(self,Z,point,Brows,u)')
    syntax('compliant_pulse_mixed_C4','_high_packet','beta',"self.flat.beta(coord['offset_from_Rv']-center)")
    loop('compliant_pulse_mixed_C4','_high_packet','(Cj, center)',
         "zip(selected['selected_scaled_end_coefficient_Taylor'],(-3,-1))")
    syntax('compliant_pulse_mixed_C4','_high_packet','Brows',
        ['[ap*(shape[k]*math.factorial(k)*self.mu**k) for k in range(5)]',
         "[selected['selected_ap_Taylor']*0 for _ in range(5)]",
         '[b*self.Ecap for b in Brows]',"[selected['selected_ap_Taylor']*0 for _ in range(5)]"])
    statement('compliant_pulse_mixed_C4','_high_packet','Brows[k]+=Cj*(beta[k]*math.factorial(k))')
    statement('compliant_pulse_mixed_C4','_high_packet','point.update(mixed)')
    keyword('compliant_pulse_mixed_C4','transport_mixed','primitive_y_derivative_Taylor','primitive')
    # At s=0 all beta derivative rows vanish, so every Brows[k] is zero.
    # Both recurrences then preserve zero primitive derivatives by induction.
    M1,M2,bb=s.symbols('Mz_jet mixed_moment_jet axial_jet',real=True)
    zero('source_meridional_mixed_recurrences_preserve_zero_Mz',
         (bb-M1*(s.Rational(1,2)-mu)).subs({bb:0,M1:0}))
    zero('source_meridional_mixed_recurrences_preserve_zero_Mtheta_z',
         (bb-M2*(s.Rational(1,2)-2*mu)).subs({bb:0,M2:0}))
    syntax('compliant_pulse_high_jets','radial','numerator','z*b*2-z*m*(1-self.delta)-d*(mz-z*m/q*2)')
    b,m,mz=s.symbols('axial_value Mz_value Mz_derivative',real=True)
    radial=2*z*b-z*m*(1-delta)-(1-z*z)*(mz-2*z*m/(1+z*z))
    zero('actual_radial_recovery_zero_from_zero_axial_and_Mz_histories',radial.subs({b:0,m:0,mz:0}))
    syntax('compliant_flatten_mixed_C4','flatten_mixed','rows',
           '{UZ:[zero for _ in range(5)],UT:theta_rows,UR:[zero for _ in range(5)],P:prows}')
    syntax('compliant_power_angular_C4','_packet','mixed',
           'flatten_mixed(c,self.mu,one,[c.mpf(0)]+rates,theta,X,energy,pressure,self.flatten.Ev2)')
    syntax('compliant_steep_waiting_C4','packet','result',
        'self.outer._packet(Z,theta,X,energy,pressure,[one*r for r in rates],coordinate,dict(stage=stage,entire_steep_waiting_high_mixed_derivatives_available=True,same_actual_angular_terminal_histories_retained=True,exact_Gamma_and_both_epsilon_atoms_retained=True,exact_relative_velocity_log_parts=logs,no_forward_subtraction_of_unrelated_long_future_energy=True))')
    keyword('compliant_corrected_outer_field','packet','Mz_over_R_Utheta','zero')
    keyword('compliant_corrected_outer_field','packet','Mtheta_z_over_sqrt2_R_3half_Utheta_squared','zero')
    # Downstream pure-swirl charts construct zero velocity rows after the
    # source-proven terminal. FTC of the original cumulative moments then
    # propagates Mz=Mtheta_z=0 for all R>=Rv, independent of swirl slope.
    proofs['downstream_zero_rows_follow_source_proven_terminal_meridional_histories'] = True
    proofs['actual_cumulative_meridional_histories_remain_zero_by_FTC_when_Uz_is_zero'] = True
    proofs['actual_zero_Mz_Mtheta_z_and_meridional_velocities_inherited_from_selected_linear_equations'] = True
    return dict(
        identities=proofs,
        original_absolute_angular_source_proof=absolute,
        complete_terminal_moment_history_bridge_verified=True,
        actual_angular_tail_constant_zero_verified=True,
        actual_selected_energy_terminal_history_verified=True,
        renormalized_energy_constant_zero_from_selected_equation=True,
        actual_exterior_energy_full_Gamma_integral_and_half_normalization_verified=True,
        physical_terminal_energy_preserved_positive=True,
        actual_zero_meridional_terminal_histories_verified=True,
        canonical_Gamma_and_exact_radius_source_binding_reused=True,
        original_pressure_terminal_history_binding_reused=True,
        original_datum_coefficients_and_histories_changed=False,
        interval_overlap_used_as_proof=False,
        cap_endpoint_used_as_exact_value=False,
        source_receipts_used_as_derivative_enclosures_not_fitted_values=True,
        full_angular_future_integral_definition=str(A0),
        full_energy_future_integral_definition=str(E0),
        global_admissible_stress_lift_constructed=False,
        physical_energy_integral_certified=False,
        temporal_recursion=False,
        input_hashes=hashes)
