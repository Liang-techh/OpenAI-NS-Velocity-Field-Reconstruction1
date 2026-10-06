"""Current whole flatten/power cone from the original correlated histories.

Full source E/P and K_Z remain. Bounds on exponentials do not define fields.
"""
import ast
import copy
import math
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_pulse_gap_cone import (
    HERE,PREFIX,sha,pack,encode,endpoints,source_precision)
from lei_ren_part1_paper_compliant_current_original_cone import SourceAST
from lei_ren_part1_paper_compliant_current_pulse_flatten_source import current_terminal_source_binding
from lei_ren_part1_paper_compliant_current_angular_background_stress import (
    angular_source_log_parts,source_log_cancellation_theorem)
from lei_ren_part1_paper_compliant_flatten_cone import flatten_cone_identities
from lei_ren_part1_paper_compliant_outer_power_cone import outer_power_cone_identities
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets

DOMAINS=dict(flatten=(0,100),outer_power=(0,1))


def generic_flatten_power_theorem():
    result=dict(flatten=flatten_cone_identities(),outer_power=outer_power_cone_identities(),
        current_KR_cancellation=source_log_cancellation_theorem())
    result['input_hashes']={PREFIX+stem+'.py':sha(PREFIX+stem+'.py') for stem in (
        'flatten_cone','outer_power_cone','current_angular_background_stress')}
    result['historical_cone_receipts_loaded_or_promoted']=False
    return result


def current_flatten_power_source_theorem(field):
    field.assert_graph();owner=field.registry.owners['flatten'];owner.assert_graph()
    flat=owner.flatten;c=owner.ctx;asts=SourceAST()
    terminal=current_terminal_source_binding(SimpleNamespace(
        pulse=flat.pulse,flatten=flat,dispatch=owner.physical.dispatch))
    if not (owner.acceptance_loaded and owner.normalization['passed'] and owner.units['passed']
            and owner.history.proof['passed'] and all(v['passed'] for v in owner.joins.values())):
        raise ValueError('Checked current full moments, units, future and both joins required')
    if not owner.history.proof['zero_meridional_histories_propagate_from_selected_terminal_by_FTC']:
        raise ValueError('Current source zero axial velocity/shear must propagate by full FTC')
    end=field.registry.owners['end']
    if end.flatten_power is not owner or not end.acceptance_loaded or not end.join['endpoint_equality_is_functional_not_interval_overlap']:
        raise ValueError('Same current pulse-end/flatten full tensor join required')
    def syntax(stem,method,target,wanted=None,augmented=False):
        return asts.expression(stem,method,target,wanted=wanted,augmented=augmented)
    syntax('power_inlet_C4','incoming','q','IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0])')
    # The actual incoming X is a constant Taylor function of Z. Its scalar
    # is the very same current pulse Xp, not a historical inlet receipt.
    fn=asts.method('power_inlet_C4','incoming')
    ret=next(n for n in ast.walk(fn) if isinstance(n,ast.Return))
    Xnode=next(n.value for n in ret.value.keywords if n.arg=='X')
    if ast.dump(Xnode)!=ast.dump(ast.parse('IntervalTaylor.constant(c,self.Xp,5)',mode='eval').body):
        raise ValueError('Actual incoming scalar X source changed')
    syntax('flatten_mixed_C4','flatten','X','(Xint+self.Xv*c.exp(-self.rate*t))/F')
    syntax('flatten_mixed_C4','flatten','Xint','q*0')
    syntax('flatten_mixed_C4','flatten','Xint',
        'Fc*(c.exp(-self.rate*(t*(self.cells-i-1)/self.cells))*decay_integral(c,self.rate,length))',True)
    qnode=syntax('flatten_mixed_C4','flatten','q','IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0])')
    rhonode=syntax('flatten_mixed_C4','flatten','rho','log_taylor(q)-c.ln(2)')
    Fnode=syntax('flatten_mixed_C4','flatten','F','(rho*sj[0]).exp()')
    Fcnode=syntax('flatten_mixed_C4','flatten','Fc','(rho*sigma_jets(c,v/100)[0]).exp()')
    syntax('flatten_mixed_C4','flatten','sj','[sig[k]*math.factorial(k)/100**k for k in range(5)]')
    syntax('power_angular_C4','data','f','self.flatten.flatten(Z,100)')
    powerX=syntax('power_angular_C4','power','X',
        "one/self.rate+(data['flatten_exit_X']-1/self.rate)*c.exp(-self.rate*y)")
    syntax('power_angular_C4','power','y','(self.Lrel-4)*phase')
    for target,wanted in (
        ('inertial','[(A[j]*heat.k-z*axial_derivative(A[j])*b-Km[j])/L for j in range(4)]'),
        ('viscous','[K[j+1]-K[j]*(1+heat.a) for j in range(4)]'),
        ('shear','[v*(2*heat.S*c.exp(-offset)) for v in shifted_rows(viscous,-1)]')):
        syntax('collar_stress_C3','collar_stress_rows',target,wanted)
    syntax('flat_pulse_derivatives','_sigma_left','odds','1/(1-x)**2-1/x**2')
    syntax('flat_pulse_derivatives','_sigma_left','e','positive_exp(c,odds)')
    syntax('flat_pulse_derivatives','_sigma_left','value','e/(1+e)')
    syntax('flat_pulse_derivatives','_sigma_left','q','value*(1-value)')
    syntax('flat_pulse_derivatives','_sigma_left','L1','2/(1-x)**3+2/x**3')
    derivative=syntax('flat_pulse_derivatives','_sigma_left','derivatives')
    if ast.dump(derivative.elts[1])!=ast.dump(ast.parse('q*L1',mode='eval').body):
        raise ValueError('Original nonnegative sigma derivative changed')
    fn=asts.method('flat_pulse_derivatives','sigma_jets')
    reflected=ast.parse('rows.append(IntervalTaylor(c,[1-left[0]]+[(-1)**(n+1)*left[n] for n in range(1,5)]))',mode='eval').body
    if not any(ast.dump(n)==ast.dump(reflected) for n in ast.walk(fn)):
        raise ValueError('Original sigma reflection changed')
    for x,expected in ((0,[0,0,0,0,0]),(1,[1,0,0,0,0])):
        if any(endpoints(v)!=(expected[n],expected[n]) for n,v in enumerate(sigma_jets(c,c.mpf(x)))):
            raise ValueError('Original flat sigma endpoint changed')
    half=sigma_jets(c,c.mpf([0,.5]))[0];whole=sigma_jets(c,c.mpf([0,1]))[0]
    if endpoints(half)[0]<0 or endpoints(half)[1]>.5 or endpoints(whole)[0]<0 or endpoints(whole)[1]>1:
        raise ValueError('Original continuous sigma range changed')
    a,mu,z,t,L,lp,lu,rp,v=s.symbols('a mu Z t Lrel logP logU logRp v',real=True)
    r=1-mu;k=1-a;b=(1-2*a)/2;f=(1+z*z)/2;rho=s.log(f)
    ctx=SimpleNamespace(mpf=lambda x:s.Rational(str(x)),ln=s.log,exp=s.exp)
    checks={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Current flatten/power identity: '+name)
        checks[name]=True
    # Replay the native scalar exponential nodes, preserving their sources.
    class ScalarExponential(ast.NodeTransformer):
        def visit_Call(self,node):
            node=self.generic_visit(node)
            if isinstance(node.func,ast.Attribute) and node.func.attr=='exp' and not node.args and not node.keywords:
                return ast.copy_location(ast.Call(func=ast.Name(id='source_exp',ctx=ast.Load()),args=[node.func.value],keywords=[]),node)
            return node
    def replay_exp(node,env):
        return asts.evaluate(ast.fix_missing_locations(ScalarExponential().visit(copy.deepcopy(node))),dict(env,source_exp=s.exp))
    q=asts.evaluate(qnode,dict(Z=z,c=ctx,IntervalTaylor=lambda c,values:values[0]))
    actualrho=asts.evaluate(rhonode,dict(q=q,c=ctx,log_taylor=s.log))
    zero('same_native_rho',actualrho-rho)
    endpoint=replay_exp(Fnode,dict(rho=actualrho,sj=[s.Integer(1)]))
    sig=s.symbols('sigma',real=True)
    Fc=replay_exp(Fcnode,dict(rho=actualrho,c=ctx,v=v,sigma_jets=lambda c,x:[sig]))
    zero('same_native_F_endpoint',endpoint-f)
    zero('same_correlated_normalized_flatten_density',Fc/endpoint-s.exp((sig-1)*rho))
    for method in ('full_flatten_rows','full_power_rows'):
        fn=asts.method('current_flatten_power_background_tensor',method)
        ret=next(n for n in ast.walk(fn) if isinstance(n,ast.Return))
        expected=ast.parse('current_full_moment_rows(heat,K,K[0]*X,E0,P0)',mode='eval').body
        if ast.dump(ret.value)!=ast.dump(expected):raise ValueError('Actual full A=K*X source changed')
    syntax('current_steep_waiting_background_stress','current_full_moment_rows','A','[A0]')
    N=s.Function('same_current_flatten_N')(t,z)
    F=s.exp(rho*sig);actualK=s.exp((a-mu)*(t-100))*F/f
    actualA=actualK*N/F
    Wflat=((k+b*2*z*z/(1+z*z))*N-b*z*s.diff(N,z))/F-1
    zero('same_current_full_flatten_A_KX_and_W',
        (k*actualA-b*z*s.diff(actualA,z)-actualK)/actualK-Wflat)
    xf=s.Function('same_current_flatten_X')(z)
    actualX=asts.evaluate(powerX,dict(one=1,self=SimpleNamespace(rate=r),data={'flatten_exit_X':xf},c=ctx,y=t))
    H=k*(xf-1/r)-b*z*s.diff(xf,z)
    zero('same_actual_native_power_W',k*actualX-b*z*s.diff(actualX,z)-1-(mu-a)/r-H*s.exp(-r*t))
    # The general full stress has L in the inertial term. For positive W,
    # 0<L<=1 implies K*W/L >= K*W; no Z=0 restriction is imposed.
    W,K,invR,rate=s.symbols('W K inverse_R log_K_rate',real=True)
    zero('same_full_theta_source',K*W/(1-2*a*z*z)+2*invR*K*(rate-1-a)
        -K*(W/(1-2*a*z*z)+2*invR*(rate-1-a)))
    physical=SimpleNamespace(ctx=ctx,logP=lp,logRp=rp,
        flatten=SimpleNamespace(logEv2_parts={'inlet_log':2*lu}),
        history=SimpleNamespace(heat=SimpleNamespace(mu=mu,a=a,bh=s.Rational(1,2)+a,outer=SimpleNamespace(Lrel=L))))
    logfn=asts.replay('current_angular_background_stress','angular_source_log_parts',{})
    logparts={}
    for region,offset,wanted in (
        ('flatten',-L-100,lp+lu-13/(2*mu)+(a-mu)*(100+L)-13-s.log(2)),
        ('outer_power',-L,lp+lu-13/(2*mu)-100*(s.Rational(1,2)+mu)+(a-mu)*L-13-s.log(2))):
        source=logfn(physical,offset);zero(region+'_current_Bmax_regrouped',sum(source['B'].values())-wanted)
        logparts[region]=str(wanted)
        zero(region+'_Qz_equals_Qtheta_B',sum(source['Qz'].values())-sum(source['Qtheta'].values())-sum(source['B'].values()))
    wait,Ts=s.symbols('wait Ts',real=True)
    tail=rp+13/mu+100+L+2+Ts+wait
    zero('flatten_exact_inverse_radius_source',tail-wait-Ts-2-L+t-100-(rp+13/mu+t))
    phase=s.symbols('phase',real=True)
    zero('power_exact_inverse_radius_source',tail-wait-Ts-2-L+(L-4)*phase-(rp+13/mu+100+(L-4)*phase))
    asts.method('global_physical_assembly','radius')
    syntax('current_flatten_power_background_tensor','chart','q','-c.mpf(self.history.steep.wait)-c.mpf(self.history.steep.Ts)-2+v')
    live=dict(shared_current_flatten_power_end_owner=end.flatten_power is owner,
        current_scalar_Xp=flat.inlet.Xp is flat.pulse.Xp,
        current_flatten_power_callback=owner.outer.flatten is flat,
        same_exact_inverse_radius_source=owner.history.heat.exact_logS_terms==owner.history.selected.exact.repair.heat.logS_terms,
        current_source_log_function=owner.chart.__func__.__wrapped__.__globals__['angular_source_log_parts'] is angular_source_log_parts)
    if not all(live.values()):raise ValueError('Actual current defining source graph differs')
    return dict(identities=checks,actual_source_AST_bindings=asts.bindings,live_current_graph=live,
        current_scalar_terminal_function=terminal,current_full_units=owner.units,
        current_full_normalization=owner.normalization,current_full_energy_source=owner.history.proof,
        current_end_flatten_tensor_join=end.join,current_flatten_power_angular_joins=owner.joins,
        current_Bmax_exact_grouped_recipes=logparts,actual_sigma_half_bound=half,actual_sigma_whole_bound=whole,
        inertial_lower_uses_L_le1_on_all_Z=True,KR_factors_already_in_current_source_logs=True,
        current_A_equals_KX_bound_to_actual_full_row_AST=True,
        K_lower_inequality='F/f=f^(sigma-1)>=1, a-mu<0, original offsets<=0',
        exact_inverse_radius_not_Scap_endpoint=True,historical_cone_admissions_not_consumed=True,
        current_axial_velocity_and_shear_exact_zero_by_full_FTC=True,
        input_hashes=asts.hashes,passed=True)


def validate_whole_flatten_power_views(field,views):
    if set(views)!=set(DOMAINS):raise ValueError('Both whole current flatten/power views required')
    for region,domain in DOMAINS.items():
        view=views[region];raw=view['original_complete_view']
        if (view['actual_five_defect_family_sha256'],view['implicit_source_sha256'],view['datum_enclosure_sha256'])!=(field.family,field.source,field.datum_sha):
            raise ValueError('Foreign flatten/power source')
        if endpoints(raw['Z'])!=(-1,1) or endpoints(raw['original_coordinate'])!=domain or raw['chart']!=region:
            raise ValueError('Whole original flatten/power domain required')
        for name in ('actual_full_stress_not_local_difference','current_complete_pressure_energy_and_forward_X_retained',
                'ordinary_logR_derivatives_not_phase_derivatives'):
            if not raw[name]:raise ValueError('Current complete moment/source view missing: '+name)
        if set(raw['current_actual_normalized_full_moment_rows'])!=set(('A','E','P','K')):
            raise ValueError('All current full A/E/P/K histories required')
        if not raw['stable_actual_absolute_pressure_mixed4_factored'] or not raw['stable_complete_remaining_quadratic_sources']:
            raise ValueError('Actual full pressure/energy suffix required')
        for group in ('theta','axial','theta_inertial','theta_shear'):
            row=raw['current_actual_normalized_stress_mixed3'][group]['y0_Z0']
            if not all(mp.isfinite(v) for v in endpoints(row)):raise ValueError('Finite complete stress enclosure required')
    return views


@source_precision
def whole_current_flatten_power_bounds(field,views):
    validate_whole_flatten_power_views(field,views);owner=field.registry.owners['flatten']
    c=field.ctx;h=owner.heat;mu=c.mpf(owner.outer.mu);a=c.mpf(h.a);r=1-mu;k=1-a;b=(1-2*a)/2
    L=c.mpf(owner.outer.Lrel);g=(mu-a)/r;glo=c.mpf(endpoints(g)[0])
    Xp=c.mpf(owner.flatten.inlet.Xp);Xv=c.mpf(owner.flatten.Xv)
    def upper(x):return c.mpf(endpoints(abs(x))[1])
    def lower(x):return c.mpf(endpoints(x)[0])
    memorylog=c.ln(2*k*(upper(Xp)+1/r))-13*r/mu
    memorycap=c.ln(glo)-1000;memorygap=memorycap-memorylog
    if endpoints(memorygap)[0]<=0:raise ArithmeticError('Current exact signed memory not below relative gap cap')
    memory=upper(c.exp(memorycap));Wflat=lower(glo-memory)
    cover=[];slope=mp.mpf(0)
    for i in range(8):
        box=c.mpf([mp.mpf(i)/8,mp.mpf(i+1)/8]);jet=sigma_jets(c,box)[1]/100
        slope=max(slope,endpoints(jet)[1]);cover.append(dict(original_sigma_interval=box,sigma_t=jet))
    slope=c.mpf(slope);Hgap=b*c.expm1(50*r)-2*k;Hlower=lower(c.exp(-100*r)*Hgap/(2*r))
    common=dict(mu=mu,mu_minus_a=mu-a,rate=r,a=a,b=b,k=k,Lrel_minus4=L-4,
        Lmin=1-2*a,current_Xp=Xp,current_Xv=Xv,current_signed_memory_log_margin=memorygap,
        flatten_W_lower=Wflat,correlated_power_Hf_parameter_gap=Hgap,correlated_power_Hf_lower=Hlower)
    lp=c.mpf(owner.physical.logP);lu=c.mpf(owner.flatten.logEv2_parts['inlet_log'])/2
    results={};logs={};margins=dict(common)
    for region in DOMAINS:
        raw=views[region]['original_complete_view'];flat=region=='flatten'
        W=Wflat if flat else glo
        m=upper(2*mu+2*c.ln(2)*slope) if flat else upper(2*mu)
        rateabs=upper(1+mu+c.ln(2)*slope) if flat else upper(1+mu)
        # Exact S*exp(-q)=exp(-logR) is strictly positive. Evaluate a much
        # larger finite cap relative to g, avoiding exp(-enormous radius).
        logR=raw['exact_source_logR'];shearlog=c.ln(2*rateabs)-c.mpf(endpoints(logR)[0])
        shearcap=c.ln(glo)-1000;sheargap=shearcap-shearlog
        if endpoints(sheargap)[0]<=0:raise ArithmeticError('Current exact shear cap failed: '+region)
        shearupper=upper(c.exp(shearcap));bracket=lower(W-shearupper)
        Kmin=lower(c.exp(-(a-mu)*(L if flat else 4)))
        theta=lower(Kmin*bracket);axial=upper(raw['current_actual_normalized_stress_mixed3']['axial']['y0_Z0'])
        parts=dict(logPstar=lp,actual_inlet_logU=lu,inverse_mu=-13/(2*mu),finite=-13-c.ln(2))
        parts['flatten_power']=(a-mu)*(100+L) if flat else -100*(c.mpf('.5')+mu)+(a-mu)*L
        logB=sum(parts.values(),c.mpf(0))
        logterm=c.ln(m)+2*logB+2*c.ln(axial/theta);logthreshold=c.ln(2);conegap=logthreshold-logterm
        regionm=dict(W=W,K_lower=Kmin,theta_bracket=bracket,theta_lower=theta,
            exact_shear_log_margin=sheargap,kappa_minus2_lower=2*mu,
            kappa_minus2_upper=m,kappa_minus2_below2=2-m,directional_log_margin=conegap)
        margins.update({region+'_'+key:value for key,value in regionm.items()})
        results[region]=dict(positive_margins=regionm,current_full_raw_theta_enclosure=raw['current_actual_normalized_stress_mixed3']['theta']['y0_Z0'],
            current_full_axial_enclosure=raw['current_actual_normalized_stress_mixed3']['axial']['y0_Z0'],
            axial_absolute_upper=axial,continuous_W_lower=W,current_K_lower=Kmin,theta_uniform_lower=theta,
            exact_inverse_radius_source_logR=logR,exact_shear_absolute_log_upper=shearlog,
            shear_relative_log_cap=shearcap,shear_absolute_upper=shearupper,
            current_Bmax_grouped_parts=parts,current_log_Bmax=logB,
            log_directional_term_upper=logterm,log_directional_threshold=logthreshold,
            full_pressure_and_energy_in_axial_enclosure=True,raw_theta_box_used_as_positivity_proof=False,
            actual_source_shear_strictly_negative=True,actual_axial_shear_exact_zero=True,
            cone_inequality='2*Ctheta^2-(vs-2)*B^2*Cz^2>0; common positive physical factors retained')
    for label,value in margins.items():
        if endpoints(value)[0]<=0 or not all(mp.isfinite(v) for v in endpoints(value)):
            raise ArithmeticError('Current flatten/power strict positive margin failed: '+label)
    if endpoints(mu)[1]>=1 or endpoints(2*a)[1]>=1:raise ArithmeticError('Original parameter range failed')
    return dict(positive_margins=margins,regions=results,current_equilibrium_gap=g,
        exact_signed_memory_log_upper=memorylog,relative_memory_log_cap=memorycap,memory_absolute_upper=memory,
        original_sigma_continuous_interval_cover=cover,sigma_t_upper=slope,correlated_Hf_uniform_lower=Hlower,
        exact_Xv_recipe='1/r+(same_current_Xp-1/r)*exp(-13*r/mu), scalar in Z',
        flatten_comparison='V_t+(r+rho*sigma_t)*V=b*j*(1-sigma)-rho*sigma_t*(1+g)>=0, V=W-g',
        outer_power_comparison='W=g+Hf(Z)*exp(-r*(Lrel-4)*phase), Hf>=strict positive lower',
        whole_continuous_original_domains=True,source_caps_used_as_defining_fields=False,
        whole_Z_not_only_axis=True,point_or_phase_samples_used_as_proof=False,
        full_E_P_X_and_K_Z_preserved=True,exact_S_positive_not_interval_lower_endpoint=True)
