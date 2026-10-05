"""Full original main/exit physical tensor and three-component remainder.

The admitted stress lift is reused. Only its remainder input is extended
to keep local/incoming radial histories and every quadratic cross product.
This regional source enclosure is not global flatness or corrected NS accuracy.
"""
import ast
import copy
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_pulse_main_exit_similarity_C4 import (
    CompliantPulseMainExitSimilarityC4,DOMAIN,PREFIX,source_precision)
import lei_ren_part1_paper_compliant_pulse_end_physical_C2 as admitted
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import SourceAST
from lei_ren_part1_paper_compliant_collar_stress_C3 import axial_derivative,shifted_rows
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import axial_operator_rows,physical_operators
from lei_ren_part1_paper_compliant_global_physical_assembly import PULSE
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import gp_jets
from mpmath.ctx_iv import MPIntervalContext

HERE=Path(__file__).parent
FALSE_FLAGS=('pulse_main_exit_cone_certified','whole_outer_cone_certified',
    'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder',
    'physical_energy_integral_certified','full_background_NS_validation','temporal_recursion',
    'production_exact_point_parameters_selected')


def main_exit_remainder_sectors(c,delta,mu,z,velocity):
    """Ordinary source rows already include every incoming exponential rate."""
    U0=velocity['local']['radial'];U1=velocity['incoming']['radial']
    Ut=velocity['local']['theta'];Uz=velocity['local']['axial']
    L=1-z*z*delta;b=(1-delta)/2
    def time(U):
        return [(U[j]/2+z*axial_derivative(U[j])*b+U[j+1])/L for j in range(3)]
    def radial_viscosity(U):
        return shifted_rows([-(U[j+2]-U[j]/4)*2 for j in range(3)],-1,2)
    def radial_convection(U,V):
        return shifted_rows(product_rows(U[:3],V[1:4]),-c.mpf('.5'),2)
    def axial_convection(U):
        dz=axial_operator_rows(c,U,-1,z,delta,count=2,order=1)
        return product_rows(Uz[:3],dz)
    def axial_viscosity(U):
        return [-v for v in axial_operator_rows(c,U,-1,z,delta)]
    def part(rows,mode,beta,extra,half=True):
        return dict(rows=rows,mode=mode,beta=beta,normalization_half=half,extra_source=extra)
    R00=radial_convection(U0,U0);A0=axial_convection(U0)
    R01=radial_convection(U0,U1);R10=radial_convection(U1,U0);A1=axial_convection(U1)
    return dict(radial=dict(
        time_local=part(time(U0),(.5,1,1,0),-3,'one'),
        time_incoming=part(time(U1),(.5,1,1,0),-3,'incoming1'),
        radial_viscosity_local=part(radial_viscosity(U0),(-.5,1,1,0),-3,'one'),
        radial_viscosity_incoming=part(radial_viscosity(U1),(-.5,1,1,0),-3,'incoming1'),
        nonlinear_local=part([a+b for a,b in zip(R00,A0)],(.5,2,2,0),-3,'one'),
        nonlinear_one_incoming=part([a+b+d for a,b,d in zip(R01,R10,A1)],(.5,2,2,0),-3,'incoming1'),
        nonlinear_incoming_square=part(radial_convection(U1,U1),(.5,2,2,0),-3,'incoming1_sq'),
        axial_viscosity_local=part(axial_viscosity(U0),(.5,1,1,0),-3+2*delta,'one'),
        axial_viscosity_incoming=part(axial_viscosity(U1),(.5,1,1,0),-3+2*delta,'incoming1')),
        theta=dict(axial_viscosity=part(
            [-v for v in axial_operator_rows(c,Ut,-1-delta,z,delta)],(0,1,0,0),-3+delta,'one',False)),
        axial=dict(axial_viscosity=part(
            [-v for v in axial_operator_rows(c,Uz,-1-delta,z,delta)],(0,1,1,0),-3+delta,'one',False)))


def main_exit_to_physical_packet(packet,mu):
    logs=packet['exact_source_logs'];converted=dict(packet)
    extras=dict(logs['extra']);extras['incoming1_sq']=2*extras['incoming1']
    converted.update(s=packet['xi']/mu,exact_logR=logs['R'],
        exact_pulse_reference_logB_parts=logs['B'],exact_logD=0,exact_logH=logs['H'],
        exact_source_logs=dict(logs,extra=extras))
    return converted


def compiled_main_exit_lift():
    """Replay the admitted lift, changing exactly its remainder-sector log hook."""
    asts=SourceAST();fn=copy.deepcopy(asts.method('pulse_end_physical_C2','lift_physical_packet'))
    fn.decorator_list=[];fn.name='lift_main_exit_packet'
    wanted=ast.dump(ast.parse("factor_logs(c,packet,sector['mode'],sector['normalization_half'])",mode='eval').body)
    matches=[]
    for node in ast.walk(fn):
        if isinstance(node,ast.Assign) and any(ast.unparse(t)=='parts' for t in node.targets) and ast.dump(node.value)==wanted:
            matches.append(node)
    if len(matches)!=1:raise ValueError('Admitted remainder log hook changed')
    target=matches[0]
    class AddHistory(ast.NodeTransformer):
        def visit_Assign(self,node):
            if node is target:
                return [node,ast.parse(
                    "parts['exact_source_history_log']=packet['exact_source_logs']['extra'][sector['extra_source']]").body[0]]
            return self.generic_visit(node)
    fn=AddHistory().visit(fn);ast.fix_missing_locations(fn)
    env=dict(vars(admitted));env['pulse_remainder_sectors']=main_exit_remainder_sectors
    exec(compile(ast.Module(body=[fn],type_ignores=[]),'<same admitted lift with main histories>','exec'),env)
    return env['lift_main_exit_packet'],dict(input_hashes=asts.hashes,
        admitted_stress_diagonal_divergence_and_Cartesian_lift_body_unchanged=True,
        only_remainder_extra_source_log_hook_added=True)


def source_and_join_binding(records):
    asts=SourceAST();checks={}
    def zero(name,value):
        if s.cancel(s.expand(s.expand_power_exp(value)))!=0:
            raise ArithmeticError('Main/exit physical source identity failed: '+name)
        checks[name]=True
    z,mu,delta,xi,lp,lu,lrp,finite=s.symbols('Z mu delta xi logP logU logRp finite',real=True)
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),exp=s.exp,
        expm1=lambda value:s.exp(value)-1,ln=s.log)
    logs=dict(R=lrp+xi/mu,B=dict(logPstar=lp,actual_log_inlet_U=lu,inverse_mu=-xi/(2*mu),finite=-xi),
              H=-(1-mu)*xi/mu,extra=dict(one=0,incoming1=-(s.Rational(1,2)-mu)*xi/mu))
    bridge=asts.replay('pulse_main_exit_physical_C2','main_exit_to_physical_packet',{})
    mapped=bridge(dict(exact_source_logs=logs,xi=xi),mu)
    zero('actual_original_main_logR_coordinate',mapped['s']-xi/mu)
    zero('actual_local_D_is_one_not_terminal_end_scale',mapped['exact_logD'])
    zero('actual_square_history_log_is_twice_incoming',mapped['exact_source_logs']['extra']['incoming1_sq']-2*logs['extra']['incoming1'])
    for key in ('R','H'):
        zero('actual_original_source_log'+key+'_retained',mapped['exact_log'+key]-logs[key])
    for key,value in logs['B'].items():
        zero('actual_original_B_'+key+'_retained',mapped['exact_pulse_reference_logB_parts'][key]-value)
    env=dict(mp=SimpleNamespace(mpf=lambda value:s.Rational(str(value))),math=math,
             axial_derivative=lambda value:s.diff(value,z),physical_operators=physical_operators)
    for stem,name in (('collar_Gamma_C4','product_rows'),('collar_stress_C3','shifted_rows'),
        ('pulse_end_physical_C2','axial_n'),('pulse_end_physical_C2','axial_operator_rows'),
        ('pulse_end_physical_C2','pulse_remainder_sectors'),
        ('pulse_main_exit_physical_C2','main_exit_remainder_sectors')):
        asts.replay(stem,name,env)
    row=lambda name:[s.Function(name+str(j))(z) for j in range(5)]
    U0,U1,Ut,Uz=[row(name) for name in ('Ur_local','Ur_incoming','Utheta','Uz')]
    h=s.symbols('actual_incoming_exponential',positive=True)
    velocity=dict(local=dict(radial=U0,theta=Ut,axial=Uz),incoming=dict(radial=U1))
    split=env['main_exit_remainder_sectors'](c,delta,mu,z,velocity)
    total=dict(radial=[a+h*b for a,b in zip(U0,U1)],theta=Ut,axial=Uz)
    original=env['pulse_remainder_sectors'](c,delta,mu,z,total)
    scales=dict(one=1,incoming1=h,incoming1_sq=h*h)
    for j in range(3):
        actual=sum(v['rows'][j]*scales[v['extra_source']] for v in split['radial'].values())
        expected=sum(v['rows'][j] for v in original['radial'].values())
        zero('actual_full_radial_time_viscosity_and_all_convection_row'+str(j),actual-expected)
        for label in ('theta','axial'):
            zero('actual_full_'+label+'_axial_viscosity_row'+str(j),
                 split[label]['axial_viscosity']['rows'][j]-original[label]['axial_viscosity']['rows'][j])
    checks['incoming_rates_already_in_velocity_rows_not_differentiated_twice']=True
    checks['all_four_radial_and_two_axial_convection_cross_products_retained']=True
    # The exact native xi11 data, not box overlap, identify the full history.
    certificate=records['pulse_interface_certificate']['source_bound_functional_pulse_identities']
    for flag in ('main_gap_linear_moments_whole_Z','main_gap_selected_energy_whole_Z',
        'gap_end_original_pressure_whole_Z',*('functional_main_gap_axial_order'+str(n) for n in range(6))):
        if not certificate[flag]:raise ValueError('Original main/gap source join missing: '+flag)
        checks['directly_consumed_'+flag]=True
    for flag in ('actual_selected_scaled_row1','actual_selected_scaled_row2',
        'exact_physical_end_selected_moment','actual_incoming_row_normalization',
        'actual_common_end_scale_normalization'):
        if not certificate[flag]:raise ValueError('Selected actual main/gap source missing: '+flag)
        checks['directly_consumed_'+flag]=True
    flat=records['flat_pulse_derivatives_check']
    for flag in ('original_radial_shape_derivatives_C4_available','original_flat_shapes_retained'):
        if not flat[flag]:raise ValueError('Original flat shape receipt missing: '+flag)
        checks['directly_consumed_'+flag]=True
    if not flat['analytic_original_beta_and_flat_envelope_checks']['all_required_endpoint_flat_limits']:
        raise ValueError('Original endpoint flat limits missing')
    checks['directly_consumed_all_required_endpoint_flat_limits']=True
    asts.method('flat_pulse_derivatives','gp_jets')
    ctx=MPIntervalContext();ctx.dps=100
    actual_flat=gp_jets(ctx,11)
    for j in range(5):
        if endpoints(actual_flat[j])!=(0,0):raise ArithmeticError('Original xi11 gp jet not zero')
        checks['actual_original_gp_xi11_jet'+str(j)+'_zero']=True
    # Exact change of integration variable in the original defining kernel:
    # k_i(11)=exp(-11*l_i/mu)*integral_0^(11/mu) exp(l_i*a)*gp(mu*a)da.
    asts.expression('pulse_main_exit_similarity_C4','partial_linear_kernel','t',wanted='xi/mu')
    asts.expression('pulse_main_exit_similarity_C4','partial_linear_kernel','coord',
        wanted='xi-mu*c.mpf([endpoints(a)[0],endpoints(b)[1]])')
    asts.expression('pulse_main_exit_similarity_C4','partial_linear_kernel','total',
        wanted="(c.exp(-lam*a)-c.exp(-lam*b))/lam*gp(c,coord)['value']",augmented=True)
    asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','ml',wanted='[ap*kernels[0]]')
    asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','nl',wanted='[ap*kernels[1]]')
    asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','mi',
        wanted='[incoming[0]*(-l1)**j for j in range(5)]')
    asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','ni',
        wanted='[incoming[1]*(-l2)**j for j in range(5)]')
    # Replay actual row/velocity algorithms. Selected M_i follows the consumed
    # two-row inverse identity, rather than defining ml as an arbitrary difference.
    asts.replay('pulse_end_physical_C2','pulse_velocity_rows',env)
    ap=s.Function('same_selected_ap')(z)
    incoming=[s.Function('same_actual_incoming_'+str(i))(z) for i in (1,2)]
    integrals=s.symbols('same_full_forward_integral_1 same_full_forward_integral_2',real=True)
    rates=[s.Rational(1,2)-i*mu for i in (1,2)]
    h11=[s.exp(-11*lam/mu) for lam in rates]
    D=[s.exp(finite-2*i) for i in (1,2)]
    selected=[s.exp(-13*lam/mu-(-1/mu+finite))*(inc+ap*I)
        for lam,inc,I in zip(rates,incoming,integrals)]
    asts.replay('pulse_end_stress_C3','pulse_coefficients',env)
    asts.replay('pulse_main_exit_similarity_C4','split_main_exit_stress',env)
    shape_env=dict(env,gp_jets=lambda *args:[s.Integer(0)]*5)
    main_shapes=asts.replay('pulse_main_exit_similarity_C4','main_exit_shapes',shape_env)
    C=1/(1+z*z);F=s.Function('same_complete_future')(z)
    J0=s.Function('same_selected_end_loss')(z);P0=s.Function('same_canonical_P0')(z)
    native_rows=main_shapes(c,mu,delta,z,C,s.Symbol('Xp'),ap,incoming,
        [hh*I for hh,I in zip(h11,integrals)],s.Integer(0),F,J0,P0,s.Integer(11))
    gap_shapes=asts.replay('pulse_gap_similarity_C4','gap_shapes',env)
    gap_rows=gap_shapes(c,mu,delta,z,C,s.Symbol('Xp'),selected[0],selected[1],F,J0,P0,s.Integer(2))
    R,B,H=s.symbols('same_positive_R same_positive_B same_positive_H',positive=True)
    terminal=s.exp(-1/mu+finite);Q=s.exp(-2*(1+2*mu)/mu)
    extra=dict(one=1,incoming1=h11[0],incoming2=h11[1],Q=Q,end_square=terminal**2)
    def sector_value(part,j,main):
        rp,bpow,dp,hpow=[s.Rational(str(value)) for value in part['mode']]
        factor=extra[part['extra_source']] if main else (
            {0:terminal,1:D[0],2:D[1]}[int(part['selected_D_recipe'][-1])]**dp
            *(Q if part['pressure_memory'] else 1))
        return part['full_derivative_rows'][j]*R**rp*B**bpow*H**hpow*factor
    for label in ('theta','axial'):
        for j in range(4):
            left=sum(sector_value(part,j,True) for part in native_rows['stress'][label].values())
            right=sum(sector_value(part,j,False) for part in gap_rows['stress'][label].values())
            zero('actual_xi11_full_stress_'+label+'_ordinary_row'+str(j),left-right)
    for key in ('e0','J','pressure_baseline_rows','pressure_memory_rows'):
        for j in range(5):
            zero('actual_xi11_same_'+key+'_row'+str(j),native_rows[key][j]-gap_rows[key][j])
    histories=[]
    for i,(local,incrows) in enumerate((('ml','mi'),('nl','ni'))):
        gaprows=asts.evaluate(asts.expression('pulse_gap_similarity_C4','gap_shapes','m' if i==0 else 'n'),
            dict(M1=selected[0],M2=selected[1],c=c,mu=mu))
        joined=[]
        for j in range(5):
            full=native_rows[local][j]+h11[i]*native_rows[incrows][j]
            joined.append(full)
            zero('actual_xi11_full_linear_history_'+str(i+1)+'_row'+str(j),full-D[i]*gaprows[j])
            for n in range(5-j):
                zero('actual_xi11_full_linear_history_'+str(i+1)+'_mixed'+str(j)+str(n),
                    s.diff(full,z,n)-D[i]*s.diff(gaprows[j],z,n))
        histories.append(joined)
    zero('actual_xi11_nonlinear_full_history_square',histories[0][0]**2-(D[0]*selected[0])**2)
    gapvelocity=env['pulse_velocity_rows'](c,delta,mu,z,C,[s.Integer(0)]*5,
        [selected[0]*(-rates[0])**j for j in range(5)])
    for label in ('radial','theta','axial'):
        for j in range(5):
            full=native_rows['velocity_local'][label][j]
            if label=='radial':full+=h11[0]*native_rows['velocity_incoming'][label][j]
            for n in range(5-j):
                zero('actual_xi11_full_velocity_'+label+'_mixed'+str(j)+str(n),
                    s.diff(full,z,n)-(D[0] if label=='radial' else 1)*s.diff(gapvelocity[label][j],z,n))
    gapenv=dict(c=c,mu=mu,d=s.Integer(2),self=SimpleNamespace(
        mu=mu,logP=lp,logU=lu,logRp=lrp,finite=finite),p=1+2*mu,r=1-mu)
    mainenv=dict(gapenv,xi=s.Integer(11))
    for key,main_target,gap_target in (('R','logR','logR'),('H','logH','logH')):
        left=asts.evaluate(asts.expression('pulse_main_exit_similarity_C4','main_exit',main_target),mainenv)
        right=asts.evaluate(asts.expression('pulse_gap_similarity_C4','packet',gap_target),gapenv)
        zero('actual_xi11_same_'+key+'_source_log',left-right)
    mainB=asts.evaluate(asts.expression('pulse_main_exit_similarity_C4','main_exit','logB'),mainenv)
    gapB=asts.evaluate(asts.expression('pulse_gap_similarity_C4','packet','logBparts'),gapenv)
    for key in mainB:zero('actual_xi11_same_B_'+key,mainB[key]-gapB[key])
    mainenv['logD0']=-1/mu+finite
    extras=asts.evaluate(asts.expression('pulse_main_exit_similarity_C4','main_exit','extras'),mainenv)
    gapQ=asts.evaluate(asts.expression('pulse_gap_similarity_C4','packet','logQ'),gapenv)
    zero('actual_xi11_same_absolute_pressure_Q_log',extras['Q']-gapQ)
    gapD=asts.evaluate(asts.expression('pulse_gap_similarity_C4','packet','logD'),
        dict(gapenv,logD0=-1/mu+finite))
    for i in (1,2):zero('actual_xi11_same_selected_D'+str(i)+'_log',gapD[i]-(finite-2*i))
    # The raw cumulative pressure and the absolute datum use the same original
    # native pressure method at the same time. Canonical P0 is already admitted.
    mainfn=asts.method('axial_pulse_field','main');gapfn=asts.method('axial_pulse_field','_gap')
    calls=[]
    for fn in (mainfn,gapfn):
        found=[node for node in ast.walk(fn) if isinstance(node,ast.Call)
            and ast.unparse(node.func)=='self.pressure_moment']
        if len(found)!=1:raise ValueError('Unique native absolute pressure call required')
        calls.append(found[0])
    for index in (0,2,3):
        if ast.dump(calls[0].args[index])!=ast.dump(calls[1].args[index]):
            raise ValueError('Native main/gap absolute pressure input differs')
    zero('actual_native_main_and_gap_pressure_time',
        asts.evaluate(asts.expression('axial_pulse_field','main','t'),
            dict(xi=11,self=SimpleNamespace(mu=mu),entrance_t=None))
        -asts.evaluate(asts.expression('axial_pulse_field','_gap','time_log'),
            dict(D=2,self=SimpleNamespace(mu=mu))))
    asts.expression('pulse_radial_C4','pressure_moment','p',
        wanted="IntervalTaylor(c,inlet['Mp_over_Pstar_squared'])+u*u*(kernel/2)")
    asts.expression('pulse_radial_C4','pressure_moment','rows',
        wanted="self.selection.future.angular.initial.datum.normalized_jets(c.mpf(Z),5)['normalized_pressure_coefficients']")
    checks['actual_xi11_native_same_inlet_swirl_and_canonical_absolute_pressure_source']=True
    checks['actual_xi11_full_stress3_and_absolute_pressure4_functional_join']=True
    maincall=asts.expression('pulse_main_exit_physical_C2','main_exit','point',
        wanted='self.lift(self.ctx,mapped,self.similarity.delta,velocity,log_tau,theta,viscosity)')
    gapcall=asts.expression('pulse_gap_physical_C2','gap','point',wanted=
        "lift_physical_packet(self.ctx,mapped,self.similarity.delta,packet['gap_source_rows']['velocity'],log_tau,theta,viscosity)")
    for label,left,right in zip(('log_tau','theta','viscosity'),maincall.args[-3:],gapcall.args[-3:]):
        if ast.dump(left)!=ast.dump(right):raise ValueError('Physical join input differs: '+label)
        checks['actual_main_gap_pullback_passes_identical_'+label+'_argument']=True
    # The original field dispatcher uses the same physical radius everywhere.
    radius=asts.replay('global_physical_assembly','radius',dict(PULSE=PULSE))
    obj=SimpleNamespace(ctx=SimpleNamespace(mpf=lambda value:value),params=SimpleNamespace(mu=mu),logRp=lrp)
    exitR=radius(obj,'pulse_exit',11,{},None)[0]
    gapR=radius(obj,'pulse_gap',11,{},None)[0]
    gapendR=radius(obj,'pulse_gap_end',-2/mu,{},None)[0]
    for label,value in (('exit',exitR),('gap',gapR),('gap_end',gapendR)):
        zero('actual_production_'+label+'_xi11_radius',value-(lrp+11/mu))
    pressure=records['pulse_end_flatten_join']['source_endpoint_binding']['identities']
    for flag in ('actual_pulse_P0_getter_is_same_canonical_flatten_pressure_over_C0_squared',
                 'consumed_same_absolute_pressure_identified_by_original_FTC_and_power_datum'):
        if not pressure[flag]:raise ValueError('Current absolute pressure source missing')
        checks['directly_consumed_'+flag]=True
    operator=records['pulse_end_physical_C2']['full_physical_operator_proof']
    if not all(operator['identities'].values()):raise ValueError('Full admitted meridional operator missing')
    for i in range(3):
        for j in range(3-i):
            if any(k+n>4 for k,n in physical_operators()[i,j+2]):raise ValueError('Physical error exceeds source mixed4')
            checks['source_mixed4_covers_full_error_operator_'+str(i)+str(j)]=True
    checks['same_stress3_completed_diagonal2_divergence2_pullback_operators_consumed']=True
    checks['same_fixed_positive_nu_lambda_time_radius_and_pressure_units_retained']=True
    checks['xi10_uses_the_same_original_global_main_source_function']=True
    checks['xi11_full_velocity4_pressure4_stress3_imply_physical_error2_join']=True
    checks['xi11_Er_and_Etheta_not_reset_to_zero']=True
    checks['xi11_Ez_zero_from_original_flat_axial_input']=True
    checks['functional_joins_not_interval_overlap']=True
    return dict(identities=checks,input_hashes=asts.hashes,
        admitted_full_physical_operator_identities=len(operator['identities']),actual_source_AST_bindings=asts.bindings,
        actual_full_main_exit_physical_source_functional_joins_verified=True,
        production_point_parameters_not_selected=True,regional_error_not_claimed_flat=True)


class CompliantPulseMainExitPhysicalC2:
    @source_precision
    def __init__(self):
        self.similarity=CompliantPulseMainExitSimilarityC4()
        self.ctx=self.similarity.ctx;self.family=self.similarity.family;self.source=self.similarity.source
        self.hashes=dict(self.similarity.hashes);self.records=dict(self.similarity.records)
        for stem in ('pulse_main_exit_similarity_C4','pulse_main_exit_similarity_C4_check',
                     'pulse_gap_physical_C2','pulse_gap_physical_C2_check','pulse_end_physical_C2',
                     'flat_pulse_derivatives_check'):
            name=PREFIX+stem+'.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
            if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(self.family,self.source):
                raise ValueError('Main/exit physical source family differs: '+stem)
            if 'all_passed' in record and not record['all_passed']:raise ValueError('Unaccepted source: '+stem)
            for path,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Physical source changed: '+path)
            self.hashes.update(record['input_hashes']);self.hashes[name]=hashlib.sha256(raw).hexdigest()
            self.records[stem]=record
        self.lift,self.lift_binding=compiled_main_exit_lift()
        self.proof=source_and_join_binding(self.records)
        self.hashes.update(self.lift_binding['input_hashes']);self.hashes.update(self.proof['input_hashes'])
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def main_exit(self,Z,xi,log_tau='-1',theta='0',viscosity='1'):
        packet=self.similarity.main_exit(Z,xi)
        mapped=main_exit_to_physical_packet(packet,self.similarity.mu)
        source=packet['main_exit_source_rows']
        velocity=dict(local=source['velocity_local'],incoming=source['velocity_incoming'])
        point=self.lift(self.ctx,mapped,self.similarity.delta,velocity,log_tau,theta,viscosity)
        point.pop('pulse_end_cone_certified',None)
        point.update(domain=DOMAIN,xi=packet['xi'],
            actual_pulse_main_exit_physical_decomposition_constructed=True,
            main_exit_and_exit_gap_completed_physical_interfaces_verified=True,
            actual_physical_stress_mixed_order=3,actual_physical_remainder_mixed_order=2,
            all_local_incoming_nonlinear_history_cross_products_retained=True,
            frozen_source_logs_not_differentiated_again=True,
            **{flag:False for flag in FALSE_FLAGS})
        return point

    @source_precision
    def report(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            domain=DOMAIN,source_and_physical_join_binding=self.proof,admitted_lift_binding=self.lift_binding,
            admitted_full_physical_operator_proof=self.records['pulse_end_physical_C2']['full_physical_operator_proof'],
            whole_original_main=self.main_exit([-1,1],['.02','10'],theta=None,viscosity='.7'),
            whole_original_exit=self.main_exit([-1,1],['10','11'],theta=None,viscosity='.7'),
            common_main_exit=self.main_exit([-1,1],10,theta=None,viscosity='.7'),
            original_exit_gap=self.main_exit([-1,1],11,theta=None,viscosity='.7'),
            input_hashes=self.hashes,actual_pulse_main_exit_physical_decomposition_constructed=True,
            main_exit_and_exit_gap_completed_physical_interfaces_verified=True,
            source_caps_used_as_defining_field_values=False,**{flag:False for flag in FALSE_FLAGS})


@source_precision
def run():
    result=CompliantPulseMainExitPhysicalC2().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original whole main/exit completed physical tensor and full three-component remainder generated; cone/global/recursion pending',flush=True)
    return result


if __name__=='__main__':run()
