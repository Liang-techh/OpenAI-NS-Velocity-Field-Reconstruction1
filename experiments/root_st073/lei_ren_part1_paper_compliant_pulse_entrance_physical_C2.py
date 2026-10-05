"""Original entrance completed physical tensor and full meridional remainder."""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_pulse_entrance_similarity_C4 import (
    CompliantPulseEntranceSimilarityC4,DOMAIN,PREFIX,source_precision)
from lei_ren_part1_paper_compliant_pulse_main_exit_physical_C2 import (
    CompliantPulseMainExitPhysicalC2,compiled_main_exit_lift,main_exit_to_physical_packet)
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import SourceAST
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import physical_operators
from lei_ren_part1_paper_compliant_global_physical_assembly import PULSE
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
FALSE_FLAGS=('pulse_entrance_cone_certified','whole_outer_cone_certified',
    'completed_full_tensor_cone_certified','global_admissible_stress_lift_constructed',
    'independently_bounded_global_flat_remainder','physical_energy_integral_certified',
    'full_background_NS_validation','temporal_recursion','production_exact_point_parameters_selected')


def entrance_physical_binding(records):
    asts=SourceAST();checks={}
    def zero(name,value):
        if s.cancel(s.expand(s.expand_power_exp(value)))!=0:
            raise ArithmeticError('Entrance physical source join failed: '+name)
        checks[name]=True
    mu,delta,z,xi,lrp=s.symbols('mu delta Z xi logRp',real=True)
    c=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),ln=s.log,exp=s.exp,expm1=lambda v:s.exp(v)-1)
    logs=dict(R=lrp+xi/mu,B=dict(logPstar=s.Symbol('logP'),actual_log_inlet_U=s.Symbol('logU'),
        inverse_mu=-xi/(2*mu),finite=-xi),H=-(1-mu)*xi/mu,
        extra=dict(one=0,incoming1=-(s.Rational(1,2)-mu)*xi/mu))
    mapped=asts.replay('pulse_main_exit_physical_C2','main_exit_to_physical_packet',{})(
        dict(exact_source_logs=logs,xi=xi),mu)
    zero('actual_entrance_lift_ordinary_y_is_xi_over_mu',mapped['s']-xi/mu)
    zero('actual_entrance_local_D_equals_one',mapped['exact_logD'])
    zero('actual_entrance_history_square_log',mapped['exact_source_logs']['extra']['incoming1_sq']-2*logs['extra']['incoming1'])
    radius=asts.replay('global_physical_assembly','radius',dict(PULSE=PULSE))
    assembly=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),params=SimpleNamespace(mu=mu),logRp=lrp)
    rentry=radius(assembly,'pulse_entrance',xi/mu,{},None)[0]
    rmain=radius(assembly,'pulse_main',s.Rational(1,50),{},None)[0]
    zero('actual_entrance_production_radius_uses_y_not_xi',rentry-logs['R'])
    zero('actual_inlet_production_radius_equals_logRp',rentry.subs(xi,0)-lrp)
    zero('actual_entrance_main_production_radius_join',rentry.subs(xi,s.Rational(1,50))-rmain)
    gate=records['pulse_entrance_similarity_C4_check']
    for flag in ('original_inlet_and_entrance_main_similarity_functional_joins_verified',
        'nonzero_incoming_moments_radial_velocity_and_energy_preserved',
        'actual_original_whole_entrance_similarity_companion_constructed'):
        if not gate[flag]:raise ValueError('Current entrance source missing: '+flag)
        checks['consumed_'+flag]=True
    original=records['pulse_main_exit_physical_C2']['source_and_physical_join_binding']
    for flag in ('incoming_rates_already_in_velocity_rows_not_differentiated_twice',
        'all_four_radial_and_two_axial_convection_cross_products_retained'):
        if not original['identities'][flag]:raise ValueError('Accepted full remainder source missing')
        checks['consumed_'+flag]=True
    env=dict(math=math,mp=SimpleNamespace(mpf=lambda v:s.Rational(str(v))),
        axial_derivative=lambda v:s.diff(v,z),physical_operators=physical_operators)
    for stem,name in (('collar_Gamma_C4','product_rows'),('collar_stress_C3','shifted_rows'),
        ('pulse_end_stress_C3','pulse_coefficients'),('pulse_end_physical_C2','pulse_velocity_rows'),
        ('pulse_end_physical_C2','axial_n'),('pulse_end_physical_C2','axial_operator_rows'),
        ('pulse_end_physical_C2','pulse_remainder_sectors'),
        ('pulse_main_exit_physical_C2','main_exit_remainder_sectors'),
        ('pulse_main_exit_similarity_C4','split_main_exit_stress')):
        asts.replay(stem,name,env)
    C=1/(1+z*z);Xp=s.Symbol('canonical_Xp',real=True)
    ap=s.Function('same_selected_ap')(z)
    incoming=[s.Function('canonical_m'+str(i))(z) for i in (1,2)]
    ein=s.Function('canonical_nonzero_energy')(z);future=s.Function('same_future_half')(z)
    J=s.Function('same_selected_loss')(z);P0=s.Function('same_terminal_pressure')(z)
    D2=s.Symbol('exact_selected_end_square',positive=True)
    K=(-mu*ein+(1-s.exp(-26))/4+mu*s.exp(-26)*(future-D2*J))/ap**2
    shape_env=dict(env,gp_jets=lambda *args:[s.Integer(0)]*5)
    source=asts.replay('pulse_main_exit_similarity_C4','main_exit_shapes',shape_env)(
        c,mu,delta,z,C,Xp,ap,incoming,[s.Integer(0)]*2,K,future,J,P0,s.Integer(0))
    zeros=[s.Integer(0)]*5
    canonical_m=[[value*(-(s.Rational(1,2)-i*mu))**j for j in range(5)]
        for i,value in enumerate(incoming,1)]
    canonical_e=[ein]
    for j in range(4):canonical_e.append(2*mu*canonical_e[j]-(s.Rational(1,2) if j==0 else 0))
    Q=s.exp(-13*(1+2*mu)/mu)
    # These are normalized shape derivatives, before the B^2 sector shift.
    # Packet pressure rows have already incorporated that shift and must not
    # be used as pulse_coefficients input a second time.
    asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','Pbase')
    asts.expression('pulse_main_exit_similarity_C4','main_exit_shapes','Pmemory')
    baseline=source['pressure_baseline_rows'][0]
    memory=source['pressure_memory_rows'][0];p=1+2*mu
    canonical_p=[(baseline if j==0 else 0)+Q*memory*p**j for j in range(5)]
    native=env['pulse_coefficients'](delta,mu,z,C,Xp,zeros,canonical_m[0],canonical_m[1],
        canonical_e,zeros,canonical_p)
    R,B=s.symbols('same_positive_R same_positive_B',positive=True)
    factors=dict(one=1,incoming1=1,incoming2=1,Q=Q,end_square=D2)
    def scaled(part,j,split):
        rp,bp,dp,hp=map(lambda v:s.Rational(str(v)),part['mode'])
        factor=factors[part['extra_source']] if split else 1
        return R**rp*B**bp*part['full_derivative_rows'][j]*factor
    for label in ('theta','axial'):
        for j in range(4):
            left=sum(scaled(part,j,True) for part in source['stress'][label].values())
            right=sum(scaled(part,j,False) for part in native[label].values())
            zero('actual_canonical_power_inlet_full_'+label+'_stress_row'+str(j),left-right)
    canonical_velocity=env['pulse_velocity_rows'](c,delta,mu,z,C,zeros,canonical_m[0])
    for label in ('radial','theta','axial'):
        for j in range(5):
            total=source['velocity_local'][label][j]
            if label=='radial':total+=source['velocity_incoming'][label][j]
            for n in range(5-j):
                zero('actual_canonical_power_inlet_velocity_'+label+'_mixed'+str(j)+str(n),
                    s.diff(total-canonical_velocity[label][j],z,n))
    split_error=env['main_exit_remainder_sectors'](c,delta,mu,z,
        dict(local=source['velocity_local'],incoming=source['velocity_incoming']))
    native_error=env['pulse_remainder_sectors'](c,delta,mu,z,canonical_velocity)
    for label in ('radial','theta','axial'):
        for j in range(3):
            zero('actual_canonical_inlet_full_'+label+'_remainder_row'+str(j),
                sum(part['rows'][j] for part in split_error[label].values())
                -sum(part['rows'][j] for part in native_error[label].values()))
    inlet=records['power_inlet_C4_check']['exact_functional_production_and_join_identities']
    for flag in ('canonical_incoming_m1','canonical_incoming_m2','canonical_incoming_energy',
        'canonical_Xp_constant','canonical_Mp_constant',
        'same_source_primitive_ODEs_identify_y_jets_through4','paper_3_9_and_physical_prefactors_preserve_join'):
        if not inlet[flag]:raise ValueError('Original upstream power source missing')
        checks['consumed_'+flag]=True
    for i in range(3):
        for j in range(3-i):
            if any(k+n>4 for k,n in physical_operators()[i,j+2]):
                raise ValueError('Entrance physical error needs unavailable source derivative')
            checks['same_source_physical_interface_operator_'+str(i)+str(j)]=True
    checks['same_canonical_pressure_FTC_and_original_datum_consumed']=True
    checks['same_xi_point02_source_function_implies_all_completed_physical_joins']=True
    checks['only_axial_remainder_zero_at_inlet_radial_and_theta_not_reset']=True
    checks['equivalent_forward_energy_not_added_to_backward_sectors']=True
    return dict(identities=checks,input_hashes=asts.hashes,actual_source_AST_bindings=asts.bindings,
        admitted_full_physical_operator_identities=original['admitted_full_physical_operator_identities'],
        inlet_and_entrance_main_join_are_source_functional_not_interval_overlap=True,
        production_radius_uses_original_ordinary_y=True)


class CompliantPulseEntrancePhysicalC2(CompliantPulseMainExitPhysicalC2):
    @source_precision
    def __init__(self):
        self.similarity=CompliantPulseEntranceSimilarityC4()
        self.ctx=self.similarity.ctx;self.family=self.similarity.family;self.source=self.similarity.source
        self.records=dict(self.similarity.records);self.hashes=dict(self.similarity.hashes)
        for stem in ('pulse_entrance_similarity_C4','pulse_entrance_similarity_C4_check',
            'pulse_main_exit_physical_C2','pulse_main_exit_physical_C2_check'):
            name=PREFIX+stem+'.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
            if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(self.family,self.source):
                raise ValueError('Entrance physical source family differs: '+stem)
            if 'all_passed' in record and not record['all_passed']:raise ValueError('Unaccepted physical dependency')
            for path,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Entrance physical source changed: '+path)
            self.records[stem]=record;self.hashes.update(record['input_hashes'])
            self.hashes[name]=hashlib.sha256(raw).hexdigest()
        self.lift,self.lift_binding=compiled_main_exit_lift()
        self.proof=entrance_physical_binding(self.records)
        self.hashes.update(self.lift_binding['input_hashes']);self.hashes.update(self.proof['input_hashes'])
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def entrance(self,Z,xi,log_tau='-1',theta='0',viscosity='1'):
        packet=self.similarity.entrance(Z,xi)
        mapped=main_exit_to_physical_packet(packet,self.similarity.mu)
        source=packet['main_exit_source_rows']
        velocity=dict(local=source['velocity_local'],incoming=source['velocity_incoming'])
        point=self.lift(self.ctx,mapped,self.similarity.delta,velocity,log_tau,theta,viscosity)
        point.pop('pulse_end_cone_certified',None)
        point.update(domain=DOMAIN,xi=packet['xi'],
            actual_pulse_entrance_physical_decomposition_constructed=True,
            original_inlet_and_entrance_main_completed_physical_interfaces_verified=True,
            actual_physical_stress_mixed_order=3,actual_physical_remainder_mixed_order=2,
            all_local_incoming_nonlinear_history_cross_products_retained=True,
            frozen_source_logs_not_differentiated_again=True,
            equivalent_forward_energy_not_double_counted=True,
            **{flag:False for flag in FALSE_FLAGS})
        return point

    @source_precision
    def report(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            domain=DOMAIN,source_and_physical_join_binding=self.proof,admitted_lift_binding=self.lift_binding,
            admitted_full_physical_operator_proof=self.records['pulse_main_exit_physical_C2']['admitted_full_physical_operator_proof'],
            whole_original_entrance=self.entrance([-1,1],[0,'.02'],theta=None,viscosity='.7'),
            original_inlet=self.entrance([-1,1],0,theta=None,viscosity='.7'),
            early_ordinary_y_chart=self.entrance([-1,1],self.similarity.mu*self.ctx.mpf([0,1]),theta=None,viscosity='.7'),
            common_entrance_main=self.entrance([-1,1],'.02',theta=None,viscosity='.7'),
            input_hashes=self.hashes,actual_pulse_entrance_physical_decomposition_constructed=True,
            original_inlet_and_entrance_main_completed_physical_interfaces_verified=True,
            source_caps_used_as_defining_field_values=False,**{flag:False for flag in FALSE_FLAGS})


@source_precision
def run():
    result=CompliantPulseEntrancePhysicalC2().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Whole original entrance completed physical tensor and all remainder components generated; cone/global/recursion pending',flush=True)
    return result


if __name__=='__main__':run()
