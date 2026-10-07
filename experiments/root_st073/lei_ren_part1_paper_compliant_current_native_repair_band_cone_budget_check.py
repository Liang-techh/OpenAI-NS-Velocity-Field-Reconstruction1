"""Independent partial physical quadrature, terminal identities and bounds."""
import json
from pathlib import Path
import time
from types import SimpleNamespace
import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_native_repair_band_cone_budget as producer
import lei_ren_part1_paper_compliant_current_native_band_function_evaluator as evaluator

HERE, PREFIX, sha = producer.HERE, producer.PREFIX, producer.sha
packets, ep = producer.packets, producer.ep


def exact_density_and_shear_identities():
    x,mu,Z=s.symbols('x mu Z',positive=True)
    A,F,G=[s.Function(q)(Z) for q in ('A','F','G')]
    power=x**(-s.Rational(1,2)-mu)
    original=producer.packets.recovery.increment_densities(A*power,0,A*F,A*G)
    class SymbolicGraph:
        constant=staticmethod(s.Rational)
        c1add=staticmethod(lambda a,b:producer.source.C1Function(a.value+b.value,a.Z+b.Z))
        c1mul=staticmethod(lambda a,b:producer.source.C1Function(a.value*b.value,a.Z*b.value+a.value*b.Z))
        c1scale=staticmethod(lambda c,a:producer.source.C1Function(c*a.value,c*a.Z))
    pair=lambda v:producer.source.C1Function(v,s.diff(v,Z))
    actual=producer.density_pairs(SymbolicGraph(),pair(A*power),pair(A*F),pair(A*G))
    for key in original:
        if s.expand(actual[key].value-original[key])!=0 or s.simplify(actual[key].Z-s.diff(original[key],Z))!=0:
            raise ArithmeticError('Independent original density/Z mismatch: '+key)
    P,F0,Fy,Gy,a,p,dI,R=s.symbols('P F0 Fy Gy a p dI R',positive=True)
    old_ratio=-(s.Rational(1,2)+mu)
    new_ratio=(old_ratio*P+Fy)/(P+F0)
    if s.simplify(-2*(new_ratio-old_ratio)+2*(Fy+(s.Rational(1,2)+mu)*F0)/(P+F0))!=0:
        raise ArithmeticError('Band shear normalization differs')
    if s.simplify(R*(p*P/R+dI)/(P+F0)-p-(R*dI/P-p*F0/P)/(1+F0/P))!=0:
        raise ArithmeticError('Band full p normalization differs')
    return dict(passed=True,independent_original_C0_and_Z_density_identities=10,
                exact_band_shear_and_full_stress_quotient_identities=2,
                actual_amplitude_Z_product_rules_checked=True)


def manufactured_function_graph(owner,c):
    g=producer.source.FunctionTransportGraph();z=g.symbol('Z')
    pair=lambda value,jet:producer.source.C1Function(value,jet)
    linear=lambda a,b:pair(g.add(g.constant(a),g.mul(g.constant(b),z)),g.constant(b))
    A=linear('7/5','2/5');mu=g.node('original_source_parameter',name='mu')
    N=g.node('shared_positive_integer',name='N',lower=160);invN=g.quotient(g.one,N,'manufactured N positive')
    rows={key:linear(a,b) for key,a,b in (('M','1/50','3/1000'),('D=(J-M)/mu','3/100','-1/500'),
        ('I','9/1000','1/1000'),('S','1/100','-1/1000'),('Cp','1/200','1/2000'))}
    # Names are supplied by the authoritative original repair convention.
    if producer.controls.ROWS[1]!='D=(J-M)/mu':raise ArithmeticError('Original divided row convention differs')
    AA=g.c1mul(A,A)
    history={key:g.c1scale(invN,g.c1mul(amplitude,rows[row])) for key,row,amplitude in
        (('m','M',A),('h','I',A),('e','S',AA),('p','Cp',AA))}
    history['k']=g.c1scale(invN,g.c1mul(AA,g.c1add(rows['M'],g.c1scale(mu,rows[producer.controls.ROWS[1]]))))
    built=dict(graph=g,history=history,amplitude=A,parameters=dict(mu=mu),N=N,
        N_scaled_targets=rows,source_family=owner.family,source_graph_sha256='synthetic_reference')
    built=producer.exact_partial_band_graph(producer.controls.exact_control_graph(built,iterations=2))
    L=c.log(2);ell=L/40;centers=[L/5,L/2,4*L/5]
    class Oracle:
        mode='synthetic_reference';source_family=owner.family
        def parameter(self,name):
            if name=='mu':return c.mpf('.04')
            raise KeyError(name)
        def source(self,*a,**kw):raise ArithmeticError('No original source points in manufactured band reference')
        def integrate(self,fn,lo,hi):
            if lo==hi:return c.mpf(0)
            points=[lo,hi]
            if lo>=1:
                points += [c.exp(ci+q*ell) for ci in centers for q in (-1,0,1) if lo<c.exp(ci+q*ell)<hi]
            elif lo<0<hi:points.append(c.mpf(0))
            return c.quad(fn,sorted(set(points)))
    return built,Oracle(),centers,ell


def independent_partial_and_terminal_reference(owner):
    c=mp.mp.clone();c.dps=45;built,oracle,centers,ell=manufactured_function_graph(owner,c)
    raw=lambda t:c.exp(-1/(1-t*t)) if abs(t)<1 else c.mpf(0)
    J0=c.quad(raw,[-1,0,1]);g=lambda i,x:raw((c.log(x)-centers[i])/ell)/(ell*J0*x)
    near=lambda a,b:abs(a-b)<=c.mpf('1e-32')*max(1,abs(a),abs(b))
    checks=terminal=0
    for n in (160,257):
        for z in (c.mpf('-.6'),c.mpf('.7')):
            value=evaluator.BandFunctionEvaluator(built,oracle=oracle,Z=z,N=n,ctx=c)
            at=lambda root,x:value(root,dict(Z=z,repair_x=x))
            h=[at(q.value,1) for q in built['finite_picard_sequence'][-1]]
            hz=[at(q.Z,1) for q in built['finite_picard_sequence'][-1]]
            A=c.mpf('1.4')+c.mpf('.4')*z;AZ=c.mpf('.4');mu=c.mpf('.04')
            def physical(t):
                F=sum(h[i+2]*g(i,t) for i in range(3))/n
                G=(h[0]*g(0,t)+h[1]*g(2,t))/n
                FZ=sum(hz[i+2]*g(i,t) for i in range(3))/n
                GZ=(hz[0]*g(0,t)+hz[1]*g(2,t))/n
                E=A*t**(-c.mpf('.5')-mu);EZ=AZ*t**(-c.mpf('.5')-mu)
                dE,dV=A*F,A*G;dEZ,dVZ=AZ*F+A*FZ,AZ*G+A*GZ
                densities=producer.packets.recovery.increment_densities(E,c.mpf(0),dE,dV)
                jets=dict(m=dVZ,h=dEZ,k=(EZ+dEZ)*dV+(E+dE)*dVZ,
                    e=2*dV*dVZ-EZ*dE-E*dEZ-dE*dEZ,p=EZ*dE+E*dEZ+dE*dEZ)
                return densities,jets
            for x in (c.mpf(1),c.exp(centers[1]),c.mpf(2)):
                for key,rate in producer.allN.RATES.items():
                    rr=c.mpf(rate.numerator)/rate.denominator
                    for row,order in (('value',0),('Z',1)):
                        inc=at(getattr(built['history'][key],row),1)
                        independent=x**(-rr)*(inc+oracle.integrate(lambda t:t**(rr-1)*physical(t)[order][key],c.mpf(1),x))
                        observed=at(getattr(built['partial_band_histories'][key],row),x)
                        if not near(observed,independent):raise ArithmeticError('Partial exact function not equal to physical quadrature: '+key+'/'+row)
                        checks+=1
                        if x==2:
                            endpoint=at(getattr(built['finite_terminal_residual_identities'][key],row),x)
                            if not near(endpoint,independent):raise ArithmeticError('Finite terminal residual identity differs: '+key+'/'+row)
                            terminal+=1
    return dict(passed=True,independent_partial_physical_C0_Z_comparisons=checks,
        independent_finite_terminal_residual_identity_comparisons=terminal,
        manufactured_Z_values=['-.6','.7'],N_values=[160,257],endpoints_and_inside_bump_checked=True,
        amplitude_Z_and_nonzero_incoming_pressure_memory_checked=True,
        manufactured_reference_only=True,original_fixed_point_or_point_field_not_evaluated=True)


def independent_bump_and_state_bounds():
    c=MPIntervalContext();c.dps=65;scalar=mp.mp.clone();scalar.dps=50
    L=scalar.log(2);ell=L/40;raw=lambda t:scalar.exp(-1/(1-t*t)) if abs(t)<1 else scalar.mpf(0)
    J0=scalar.quad(raw,[-1,0,1]);caps=producer.bump_caps(c,dict(radius=c.mpf(str(ell)),raw_normalization=c.mpf(str(J0))))
    count=0
    for center in (L/5,L/2,4*L/5):
        for t in map(scalar.mpf,('-.99','-.8','-.5','0','.5','.8','.99')):
            x=scalar.exp(center+ell*t);beta=raw(t);beta_prime=-2*t*beta/(1-t*t)**2
            bump=beta/(ell*J0*x);dy=beta_prime/(ell*ell*J0*x)-bump
            if abs(bump)>ep(caps['value'])[0] or abs(dy)>ep(caps['log_radius_derivative'])[0]:
                raise ArithmeticError('Independent bump or log-radius derivative outside cap')
            count+=2
    # Check the quotient factor and both shear bounds for signed perturbations.
    C=lambda x:producer.LogUpper.constant(c,x)
    poly=lambda x:producer.Poly(c,{-1:C(x)})
    errors={key:poly(v) for key,v in (('theta_linear','2.1'),('theta_quadratic','.7'),
        ('axial_linear','1.3'),('axial_quadratic','.9'))}
    state=producer.band_normalized_state(c,errors,C('3.2'),C('1.7'),C('8'),C('5'),poly('3'),poly('7'),poly('2'))
    comparisons=0
    for n in (160,257):
        for sign in (-1,1):
            F=sign*scalar.mpf(2)/n;Fy=-sign*scalar.mpf(7)/n;Gy=sign*scalar.mpf(7)/n
            power=scalar.power(2,-scalar.mpf(2)/3);alpha=scalar.mpf('.66')
            values=dict(a=-2*(Fy+alpha*F)/(power+F),b=2*Gy/(power+F))
            rel=sign*scalar.mpf(3)/n
            values['p1']=(scalar.mpf('3.2')*(scalar.mpf('2.1')+scalar.mpf('1.7')*scalar.mpf('.7'))/n-8*rel)/(1+rel)
            values['p2']=(scalar.mpf('3.2')*(scalar.mpf('1.3')+scalar.mpf('1.7')*scalar.mpf('.9'))/n-5*rel)/(1+rel)
            for key,v in values.items():
                if abs(v)>ep(c.exp(state[key].evaluate(c.ln(n)).log))[0]:raise ArithmeticError('Independent band normalized state outside bound')
                comparisons+=1
    return dict(passed=True,independent_bump_value_and_log_derivative_comparisons=count,
                independent_signed_normalized_band_state_comparisons=comparisons,
                manufactured_reference_only=True)


def original_live_provenance(owner,built,live,cone_live,got):
    if got['original_query']['source']['packet'].source_family!=owner.family:raise ArithmeticError('Original band family differs')
    if built['partial_integration_variable']==built['band_variable']:raise ArithmeticError('Partial integral captures endpoint variable')
    c=owner.ctx;checks=0
    for key,rate in producer.allN.RATES.items():
        rr=c.mpf(rate.numerator)/rate.denominator
        expected=c.ln(2) if not rate else -c.expm1(-rr*c.ln(2))/rr
        if expected._mpi_!=got['masses'][key]._mpi_:raise ArithmeticError('Own true rate mass differs')
        for row in (got['D'][key],got['DZ'][key]):
            if any(p>=0 for p in row.terms):raise ArithmeticError('Partial moment has a nondecaying N power')
            checks+=1
    for req,state,margin in ((got['requirements'],got['state'],got['margin']),
        (got['quiet'][0]['requirements'],got['quiet'][0]['state'],got['quiet'][0]['margin'])):
        for key,poly in state.items():
            if any(p>=0 for p in poly.terms):raise ArithmeticError('State includes N growth')
            total=producer.LogUpper.add(c,list(poly.terms.values()))
            if total.log is not None and ep(req[key]['sufficient_log_N_lower']-total.log+margin['log_rho_stability_positive'])[0]<0:
                raise ArithmeticError('Directed local cone threshold insufficient')
    for threshold in (cone_live['combined_logN'],got['bandN'],got['positivityN'],got['quiet'][0]['logN']):
        if ep(got['combined_logN']-threshold)[0]<0:raise ArithmeticError('Combined source/repair/active/quiet/band threshold insufficient')
    if packets.encode(live['mu'])!=packets.encode(owner.target.repair_mu):raise ArithmeticError('Band uses different source mu')
    if len(built['finite_picard_sequence'][-1])!=5 or not live['conditions']['B_mu_Z_independent_and_C1_product_norm_submultiplicative']:
        raise ArithmeticError('All five value/Z rows need the common accepted C1-ball contract')
    rho=packets.interval(c,live['conditions']['formal_control_C1_ball_radius_log']['log_absolute_upper'])
    if rho._mpi_!=got['control_radius'].log._mpi_:raise ArithmeticError('Control-ball radius source differs')
    # Inactive quiet profile is identical; incoming histories retain the full route.
    old=next(row for row in live['cells'] if row['record']['label']=='power_to_Rc')
    if any(not row.zero for key in old['values'] for row in old['values'][key].values()):
        raise ArithmeticError('Original reserved quiet profile has nonzero modulation density')
    return dict(passed=True,partial_C0_Z_negative_power_history_polynomials=checks,
        normalized_quiet_and_band_state_polynomials=8,original_source_mu_and_amplitude_route_retained=True,
        original_quiet_no_new_density_with_incoming_memory=True,
        original_band_source_query_without_ancestor_rebuild=True,
        conditional_control_ball_not_instantiated_controls=True)


@producer.allN.paired.native.inlet.source_precision
def run(owner,built,live,cone_live,got):
    began=time.monotonic()
    report=json.loads((HERE/producer.NAME).read_bytes())
    if report['source_family']!=owner.family or not report[producer.GATE]:raise ArithmeticError('Same current band producer required')
    result=dict(all_passed=True,**{producer.GATE:True},source_family=owner.family,
        independent_exact_density_and_quotient=exact_density_and_shear_identities(),
        independent_partial_and_terminal=independent_partial_and_terminal_reference(owner),
        independent_bump_and_state=independent_bump_and_state_bounds(),
        actual_original_live_provenance=original_live_provenance(owner,built,live,cone_live,got),
        read_only_review_model=dict(model='gpt-5.6-luna',reasoning_effort='max'),
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        all_chart_q_flat_input_cone_margins_certified=False,
        **dict.fromkeys(packets.OPEN,False),input_hashes={**report['input_hashes'],producer.NAME:sha(producer.NAME),
            Path(evaluator.__file__).name:sha(Path(evaluator.__file__).name),Path(__file__).name:sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (HERE/producer.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Independent repair-band partial physical integrals, finite terminal identities and local budgets PASS',flush=True)
    return result


def recorded_original_provenance(report):
    """Check saved native output against its exact producer and input graph.

    This path avoids reconstructing live ancestors after native production.
    It compares directed budgets and exact graph ancestry, not point fields.
    """
    c=MPIntervalContext();c.dps=300
    original=json.loads((HERE/producer.allN.NAME).read_bytes())
    cone=json.loads((HERE/producer.state.NAME).read_bytes())
    for name in (producer.allN.NAME,producer.state.NAME,Path(producer.__file__).name):
        if sha(name)!=report['input_hashes'][name]:raise ArithmeticError('Native generation input changed: '+name)
    if report['source_family']!=original['source_family'] or report['source_family']!=cone['source_family']:
        raise ArithmeticError('Saved original source families differ')
    oldgraph=original['exact_function_graph_nodes'];newgraph=report['exact_function_graph_nodes']
    if newgraph[:len(oldgraph)]!=oldgraph:raise ArithmeticError('Band did not extend the exact accepted original all-N graph')
    for kind in ('exact_partial_band_history_roots','exact_finite_terminal_residual_identity_roots'):
        if set(report[kind])!=set(producer.allN.RATES):raise ArithmeticError('Saved band omitted original history row')
    for key,rate in producer.allN.RATES.items():
        rr=c.mpf(rate.numerator)/rate.denominator
        expected=c.ln(2) if not rate else -c.expm1(-rr*c.ln(2))/rr
        actual=packets.interval(c,report['original_own_rate_band_masses'][key])
        if ep(actual)[0]>ep(expected)[0] or ep(actual)[1]<ep(expected)[1]:
            raise ArithmeticError('Saved directed band mass does not enclose the higher-precision own rate')
    control=report['control_C1_ball_radius']
    if control!=original['actual_uniform_repair_C1_log_N_conditions']['formal_control_C1_ball_radius_log']:
        raise ArithmeticError('Saved control budget uses a different C1 ball')
    def verify_state(rows,requirements,margin):
        rho=packets.interval(c,margin['log_rho_stability_positive'])
        for key,row in rows.items():
            terms={term['N_power']:producer.LogUpper(c,packets.interval(c,term['coefficient']['log_absolute_upper'])) for term in row['terms']}
            if any(p>=0 for p in terms):raise ArithmeticError('Saved state has nondecaying N power')
            total=producer.LogUpper.add(c,list(terms.values()))
            threshold=packets.interval(c,requirements[key]['sufficient_log_N_lower'])
            if total.log is not None and ep(threshold-total.log+rho)[0]<0:
                raise ArithmeticError('Saved local directed frequency requirement insufficient')
    verify_state(report['band_normalized_state_error_polynomials'],report['band_state_frequency_requirements'],report['original_band_power_margin'])
    quiet=report['quiet_power_records']
    if len(quiet)!=1 or quiet[0]['original_cell_label']!='power_to_Rc' or not quiet[0]['incoming_five_history_and_pressure_memory_not_zeroed']:
        raise ArithmeticError('Saved quiet interval or history reservation differs')
    verify_state(quiet[0]['state_error_polynomials'],quiet[0]['requirements'],quiet[0]['margin'])
    combined=packets.interval(c,report['source_repair_active_quiet_band_sufficient_log_N_lower'])
    for record in (cone['original_source_repair_and_active_loop_log_N_lower'],report['band_sufficient_log_N_lower'],
                   report['band_positivity_sufficient_log_N_lower'],quiet[0]['sufficient_log_N_lower']):
        if ep(combined-packets.interval(c,record))[0]<0:raise ArithmeticError('Saved combined threshold insufficient')
    for name in ('band_partial_history_C0_error_polynomials','band_partial_history_Z_error_polynomials'):
        for key,row in report[name].items():
            if any(term['N_power']>=0 for term in row['terms']):raise ArithmeticError('Saved partial history has N growth')
    return dict(passed=True,partial_C0_Z_negative_power_history_polynomials=10,
        normalized_quiet_and_band_state_polynomials=8,exact_original_all_N_graph_prefix_nodes=len(oldgraph),
        saved_native_source_and_producer_hash_binding=True,same_original_control_C1_ball_preserved=True,
        own_rate_masses_and_directed_local_frequency_bounds_independently_checked=True,
        native_output_checked_from_saved_records=True,original_ancestor_or_live_owner_rebuild=False,
        original_field_point_evaluation_or_terminal_closure_claimed=False)


def run_recorded():
    began=time.monotonic();report=json.loads((HERE/producer.NAME).read_bytes())
    if not report[producer.GATE]:raise ArithmeticError('Current native band output required')
    owner=SimpleNamespace(family=report['source_family'])
    provenance=recorded_original_provenance(report)
    result=dict(all_passed=True,**{producer.GATE:True},source_family=owner.family,
        independent_exact_density_and_quotient=exact_density_and_shear_identities(),
        independent_partial_and_terminal=independent_partial_and_terminal_reference(owner),
        independent_bump_and_state=independent_bump_and_state_bounds(),
        actual_original_live_provenance=provenance,
        read_only_review_model=dict(model='gpt-5.6-luna',reasoning_effort='max'),
        checker_uses_saved_native_output_without_reconstructing_original_owners=True,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        all_chart_q_flat_input_cone_margins_certified=False,**dict.fromkeys(packets.OPEN,False),
        input_hashes={**report['input_hashes'],producer.NAME:sha(producer.NAME),
            Path(evaluator.__file__).name:sha(Path(evaluator.__file__).name),Path(__file__).name:sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (HERE/producer.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Independent repair-band physical/terminal references and saved original native budget PASS',flush=True)
    return result


if __name__=='__main__':run_recorded()
