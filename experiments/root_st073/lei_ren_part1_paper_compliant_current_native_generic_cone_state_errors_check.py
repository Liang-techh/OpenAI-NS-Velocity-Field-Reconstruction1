"""Independent full-recovery identities and continuous cone-error checks."""
import ast
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_native_generic_cone_state_errors as producer
import lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery as recovery

HERE, PREFIX, sha = producer.HERE, producer.PREFIX, producer.sha
packets, ep = producer.packets, producer.ep


def exact_original_recovery():
    tree = ast.parse((HERE/(PREFIX+'current_generic_shear_moment_recovery.py')).read_text(encoding='utf8'))
    cls = next(q for q in tree.body if isinstance(q, ast.ClassDef) and q.name=='GenericMomentRecovery')
    fn = next(q for q in cls.body if isinstance(q, ast.FunctionDef) and q.name=='field')
    names = ('L','d','pressure','transport','Q','theta_linear','theta_quadratic','axial_linear','axial_quadratic')
    expressions = {q.targets[0].id:q.value for q in fn.body if isinstance(q,ast.Assign)
                   and len(q.targets)==1 and isinstance(q.targets[0],ast.Name) and q.targets[0].id in names}
    if set(expressions)!=set(names): raise ArithmeticError('Original recovery expressions required')
    z, de, p0 = s.symbols('Z delta original_P0',real=True)
    originals = {key:s.Function('original_'+key)(z) for key in recovery.RATES}
    defects = {key:s.Function('defect_'+key)(z) for key in recovery.RATES}
    E,V,dE,dV = [s.Function(key)(z) for key in ('E','V','delta_E','delta_V')]
    def evaluate(changed):
        env = dict(Z=z,de=de,E=E+dE if changed else E,V=V+dV if changed else V,
                   self=type('Datum',(),{'P0':p0})(),
                   axial_derivative=lambda f:s.diff(f,z))
        env.update({key:originals[key]+defects[key] if changed else originals[key] for key in recovery.RATES})
        for name in names:
            env[name] = eval(compile(ast.Expression(expressions[name]),'<independent original recovery AST>','eval'),{},env)
        return env
    old,new = evaluate(False),evaluate(True)
    got = producer.stress_increment(E,V,originals['m'],s.diff(originals['m'],z),dE,dV,
                                  defects,{key:s.diff(value,z) for key,value in defects.items()},z,de)
    keys = dict(transport='transport',pressure='pressure',radial='Q',theta_linear='theta_linear',
                theta_quadratic='theta_quadratic',axial_linear='axial_linear',axial_quadratic='axial_quadratic')
    for name, key in keys.items():
        if s.cancel(s.expand(new[key]-old[key]-got[name]))!=0:
            raise ArithmeticError('Independent original full recovery increment differs: '+name)
    return dict(passed=True,exact_full_original_recovery_AST_difference_identities=len(keys),
                original_P0_cancels_symbolically=True,nonzero_original_V_and_transport_retained=True,
                source_ast_assignments={key:ast.unparse(value) for key,value in expressions.items()})


def frozen_polynomial_and_quotient_identities():
    a,b,p1,p2,t,v,E,R,I,dI,A,N = s.symbols('a b p1 p2 t v E R I dI A N',nonzero=True)
    G = [a,a*a+b*b-2*a,p1*a-p2*b-a*a-b*b]
    G.append(2*a*G[2]**2-G[1]*(p1*b+p2*a)**2)
    aa = v/(1+t*t);bb = -aa*t;D=p1+p2*t-v;J=p2-p1*t
    targets = [aa,aa*(v-2),aa*D,aa**3*(2*D*D-(v-2)*J*J)]
    for expr,want in zip(G,targets):
        if s.cancel(expr.subs({a:aa,b:bb})-want)!=0:
            raise ArithmeticError('Active frozen G polynomial normalization differs')
    gradient_coefficients = []
    for expr in G:
        total = 0
        for var in (a,b,p1,p2):
            poly = s.Poly(s.diff(expr,var),a,b,p1,p2)
            if poly.total_degree()>5:raise ArithmeticError('Stability gradient degree exceeds H^5')
            total += sum(abs(c) for c in poly.coeffs())
        if total>208:raise ArithmeticError('Independent gradient coefficient sum exceeds 208')
        gradient_coefficients.append(int(total))
    oldp = R*I/E
    actualp = R*(I+dI)/(E*s.exp(A/N))
    difference = s.exp(-A/N)*(R*dI/E-oldp*(s.exp(A/N)-1))
    if s.simplify(actualp-oldp-difference)!=0:raise ArithmeticError('Full normalized stress quotient error differs')
    bL,Ay,By = s.symbols('bL Ay By')
    if s.simplify((s.exp(-A/N)*(bL+2*By/(N*E))-bL)
                  -((s.exp(-A/N)-1)*bL+s.exp(-A/N)*2*By/(N*E)))!=0:
        raise ArithmeticError('Axial shear exponential denominator differs')
    return dict(passed=True,frozen_G_normalization_identities=4,normalized_stress_and_axial_shear_identities=2,
                independent_gradient_coefficient_sums=gradient_coefficients,
                independent_gradient_bound='208*H^5 for H>=1; 208/512<1/2',
                stability_fraction_208_over512_less_than_half=True,H_positive_at_least2=True)


def independent_prefix_integrals():
    c = mp.mp.clone();c.dps=70;count=0
    for rate in (c.mpf(0),c.mpf(1),c.mpf('1.5')):
        for w in (c.mpf('.01'),c.mpf('2.5')):
            for incoming in (c.mpf('-.7'),c.mpf('.9')):
                density = lambda t:c.mpf('.2')+c.mpf('.8')*c.cos(3*t)
                cap = c.mpf(1)
                fullmass = w if rate==0 else -c.expm1(-rate*w)/rate
                for theta in (c.mpf(0),c.mpf('.2'),c.mpf(1)):
                    end=w*theta
                    actual=c.exp(-rate*end)*incoming+c.quad(lambda t:c.exp(-rate*(end-t))*density(t),[0,end])
                    if abs(actual)>abs(incoming)+cap*fullmass+c.mpf('1e-65'):
                        raise ArithmeticError('Independent continuous partial Duhamel bound failed')
                    count+=1
    return dict(passed=True,independent_signed_partial_Duhamel_integrals=count,
                distinct_rates_0_1_1p5_and_both_incoming_signs=True,interior_points_and_endpoints=True)


def manufactured_full_state_references():
    c=MPIntervalContext();c.dps=90
    C=lambda v:producer.LogUpper.constant(c,v)
    poly=lambda first,second=0:producer.Poly(c,{-1:C(abs(first)),-2:C(abs(second))})
    scalar=lambda row:c.mpf(ep(row)[0])
    mag=lambda row:max(abs(v) for v in ep(row))
    jt=lambda value,derivative=0:recovery.IntervalTaylor(c,[c.mpf(value),c.mpf(derivative)])
    checks=state_checks=0
    first=dict(m=c.mpf('-.3'),h=c.mpf('.2'),k=c.mpf('-.1'),e=c.mpf('.4'),p=c.mpf('-.2'))
    second=dict(m=c.mpf('.05'),h=c.mpf('-.02'),k=c.mpf('.04'),e=c.mpf('-.03'),p=c.mpf('.01'))
    firstZ={k:v*c.mpf('-.7') for k,v in first.items()};secondZ={k:v*c.mpf('.6') for k,v in second.items()}
    family=dict(actual_five_defect_family_sha256='manufactured',implicit_source_sha256='manufactured',datum_enclosure_sha256='manufactured')
    for n in (160,257,2048):
        for aval in ('-1.1','0','.9'):
            for zval in ('-.8','0','.7'):
                for de in ('.03','.49'):
                    E0,V0,m0,mZ0,B0,R0,S0 = map(c.mpf,('3','-1.7','.4','-.9','-.6','1.4','1.8'))
                    A0=c.mpf(aval);x=A0/n;Z=jt(zval,1)
                    E=jt(E0,'.2');V=jt(V0,'-.3')
                    original={k:jt('.2','.1') for k in recovery.RATES};original['m']=jt(m0,mZ0)
                    D0={k:first[k]/n+second[k]/(n*n) for k in recovery.RATES}
                    DZ0={k:firstZ[k]/n+secondZ[k]/(n*n) for k in recovery.RATES}
                    defect={k:jt(D0[k],DZ0[k]) for k in recovery.RATES}
                    theta_delta=E0*c.expm1(x);axial_delta=B0/n
                    old=recovery.GenericMomentRecovery(c,source_family=family,P0=jt('-.8','.17'),original=original)
                    new=recovery.GenericMomentRecovery(c,source_family=family,P0=old.P0,original=original,defect=defect)
                    oldfield=old.field(Z=Z,delta=de,E=E,V=V,E_y=jt('.4'),V_y=jt('-.2'))
                    newfield=new.field(Z=Z,delta=de,E=E+jt(theta_delta),V=V+jt(axial_delta),E_y=jt('.4'),V_y=jt('-.2'))
                    D={k:poly(first[k],second[k]) for k in recovery.RATES};DZ={k:poly(firstZ[k],secondZ[k]) for k in recovery.RATES}
                    dE=poly(E0*c.mpf('1.25')*c.exp(c.mpf(1)/128));dV=poly(B0)
                    bounds=producer.inertial_error_envelopes(c,C(E0),C(abs(V0)),C(abs(m0)+abs(mZ0)),dE,dV,D,DZ)
                    names=dict(pressure='absolute_pressure_over_S_squared',radial='Ur_over_S_sqrt_R_over_2',
                        theta_linear='inertial_theta_linear',theta_quadratic='inertial_theta_quadratic',
                        axial_linear='inertial_axial_linear',axial_quadratic='inertial_axial_quadratic')
                    logN=c.ln(n)
                    for key,name in names.items():
                        value=newfield[name][0]-oldfield[name][0]
                        bound=bounds[key].evaluate(logN)
                        cap=c.mpf(0) if bound.log is None else c.exp(bound.log)
                        if mag(value)>ep(cap)[0]:raise ArithmeticError('Independent full source error not enclosed: '+key)
                        checks+=1
                    ef=producer.LogUpper(c,c.mpf(1)/128);relative=poly(c.mpf('1.25')*c.exp(c.mpf(1)/128))
                    for axis in ('theta','axial'):
                        oldI=oldfield['inertial_'+axis+'_linear'][0]+S0*oldfield['inertial_'+axis+'_quadratic'][0]
                        newI=newfield['inertial_'+axis+'_linear'][0]+S0*newfield['inertial_'+axis+'_quadratic'][0]
                        p0=R0*oldI/E0;pN=R0*newI/(E0*c.exp(x));value=pN-p0
                        budget=((bounds[axis+'_linear']+bounds[axis+'_quadratic'].scale(C(S0))).scale(C(R0/E0))
                                +relative.scale(C(mag(p0)))).scale(ef).evaluate(logN)
                        if mag(value)>ep(c.exp(budget.log))[0]:raise ArithmeticError('Independent normalized p error not enclosed')
                        state_checks+=1
                    Ay,By,bL=map(c.mpf,('-.4','.3','-.8'))
                    for value,budget in ((-2*Ay/n,poly(2*abs(Ay))),
                        (c.exp(-x)*(bL+2*By/(n*E0))-bL,relative.scale(C(abs(bL)))+poly(2*abs(By)/E0).scale(ef))):
                        bound=budget.evaluate(logN)
                        if mag(value)>ep(c.exp(bound.log))[0]:raise ArithmeticError('Independent shear state error not enclosed')
                        state_checks+=1
    return dict(passed=True,independent_full_original_recovery_error_enclosures=checks,
                independent_normalized_stress_and_shear_error_enclosures=state_checks,
                N_values=[160,257,2048],manufactured_operator_references_only=True,
                original_field_point_values_not_evaluated=True,nonzero_V_P0_transport_and_all_history_Z_rows=True)


def original_live_route(owner,allNlive,got):
    if len(got['cells'])!=24 or got['branch_count']!=35:raise ArithmeticError('Whole original route/branch count differs')
    branch_count=state_count=prefix_count=0
    for original,row in zip(allNlive['cells'],got['cells']):
        if row['record']['label']!=original['record']['label']:raise ArithmeticError('Original cell order changed')
        if not original['geometry']['record']['width_and_endpoints_independent_of_Z']:
            raise ArithmeticError('Partial history Z bounds need Z-independent widths')
        D,DZ=producer.prefix_history_polynomials(owner.ctx,original)
        for key in producer.allN.RATES:
            for expected,actual in ((D[key],row['D'][key]),(DZ[key],row['DZ'][key])):
                if packets.encode(expected.record())!=packets.encode(actual.record()):
                    raise ArithmeticError('Continuous prefix cap provenance differs')
                prefix_count+=1
        for source,branch in zip(original['branches'],row['branches']):
            if source['conditional']['source']['packet'].source_family!=owner.family:raise ArithmeticError('Original branch family differs')
            for key,poly in branch['state'].items():
                if any(p>=0 for p in poly.terms):raise ArithmeticError('State error includes an N-growing/constant term')
                req=branch['requirements'][key];total=producer.LogUpper.add(owner.ctx,list(poly.terms.values()))
                if total.log is not None and ep(req['sufficient_log_N_lower']-total.log+got['margin']['log_rho_stability_positive'])[0]<0:
                    raise ArithmeticError('Directed negative-power threshold insufficient')
                state_count+=1
            branch_count+=1
    if branch_count!=35:raise ArithmeticError('Dropped original conditional branch')
    for threshold in (got['active_state_logN'],got['source_repair_conditions']['source_and_repair_sufficient_log_N_lower']):
        if ep(got['combined_logN']-threshold)[0]<0:raise ArithmeticError('Combined source/repair/active threshold insufficient')
    return dict(passed=True,original_continuous_cells=24,original_conditional_branches=branch_count,
                genuine_continuous_prefix_C0_Z_polynomials=prefix_count,
                normalized_negative_N_power_state_error_polynomials=state_count,
                original_P0_preservation_inherited_but_difference_cancellation_independently_checked=True,
                global_q_flat_or_repair_band_margin_not_claimed=True)


@producer.allN.paired.native.inlet.source_precision
def run(owner,allNlive,got):
    began=time.monotonic()
    result=dict(all_passed=True,**{producer.GATE:True},source_family=owner.family,
        independent_original_full_recovery=exact_original_recovery(),
        exact_frozen_cone_and_quotient=frozen_polynomial_and_quotient_identities(),
        continuous_prefix_reference=independent_prefix_integrals(),
        independent_state_references=manufactured_full_state_references(),
        current_original_live_source_provenance=original_live_route(owner,allNlive,got),
        read_only_review_model=dict(model='gpt-5.6-luna',reasoning_effort='max'),
        inactive_input_and_post_modulation_cone_margins_certified=False,
        correction_band_partial_moments_and_cone_errors_certified=False,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes={**owner.hashes,producer.NAME:sha(producer.NAME),Path(__file__).name:sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (HERE/producer.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Independent original full recovery, continuous prefixes and active cone state budget PASS',flush=True)
    return result
