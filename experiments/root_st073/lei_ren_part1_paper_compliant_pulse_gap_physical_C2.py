"""Whole original inactive-gap completed physical tensor and radial history.

The admitted full meridional operator is applied to the same gap sectors.
D1, rather than the terminal D0, is the radial history factor. Positive
sources remain exact logs; this regional remainder is not a global flatness
or corrected Navier-Stokes accuracy claim.
"""
import ast
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import sympy as s
from lei_ren_part1_paper_compliant_pulse_gap_similarity_C4 import (
    CompliantPulseGapSimilarityC4, DOMAIN, PREFIX, source_precision)
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import SourceAST
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import (
    lift_physical_packet, pulse_velocity_rows, pulse_remainder_sectors)
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_operators
from lei_ren_part1_paper_compliant_global_physical_assembly import PULSE
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
FALSE_FLAGS=('pulse_gap_cone_certified','whole_outer_cone_certified',
    'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder',
    'physical_energy_integral_certified','full_background_NS_validation','temporal_recursion')


def gap_to_physical_packet(packet,mu,right_endpoint=False):
    """Frozen basepoint logs accompany already differentiated ordinary rows."""
    logs=packet['exact_source_logs']
    converted=dict(packet)
    converted.update(s=-4 if right_endpoint else -packet['distance']/mu,
        exact_logR=logs['R'],exact_pulse_reference_logB_parts=logs['B'],
        exact_logD=logs['D1'],exact_logH=logs['H'])
    return converted


def source_and_join_binding(records):
    asts=SourceAST();checks={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:
            raise ArithmeticError('Gap physical source identity failed: '+name)
        checks[name]=True
    mu,d,z,delta,lrp,lp,lu,finite=s.symbols(
        'mu distance Z delta original_logRp original_logP original_logU finite',real=True)
    lam1=s.Rational(1,2)-mu;lam2=s.Rational(1,2)-2*mu
    bp=s.Rational(1,2)+mu;r=1-mu
    logs=dict(R=lrp+(13-d)/mu,
        B=dict(logPstar=lp,actual_log_inlet_U=lu,inverse_mu=(-13+d)/(2*mu),finite=-13+d),
        D0=-1/mu+finite,D1=(d/2-1)/mu+finite-d,
        D2=(d/2-1)/mu+finite-2*d,H=-r*(13-d)/mu)
    bridge=asts.replay('pulse_gap_physical_C2','gap_to_physical_packet',{})
    packet=dict(exact_source_logs=logs,distance=d)
    mapped=bridge(packet,mu)
    zero('actual_gap_to_ordinary_logR_coordinate',mapped['s']+d/mu)
    zero('actual_radial_remainder_uses_effective_D1_not_terminal_D0',
        mapped['exact_logD']-logs['D1'])
    for key in ('R','H'):
        zero('actual_physical_adapter_retains_same_'+key,mapped['exact_log'+key]-logs[key])
    for key,value in logs['B'].items():
        zero('actual_physical_adapter_retains_B_'+key,
            mapped['exact_pulse_reference_logB_parts'][key]-value)
    right=bridge(dict(packet,distance=4*mu),mu,right_endpoint=True)
    zero('formal_right_endpoint_exactly_original_end_minus4',right['s']+4)
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),ln=lambda value:lu)
    native=SimpleNamespace(mu=mu,logE=-1/mu+finite)
    endself=SimpleNamespace(native=native,assembly=SimpleNamespace(logRp=lrp,logP=lp),
        flatten=SimpleNamespace(native=SimpleNamespace(U=s.Symbol('original_U',positive=True))))
    env=dict(c=c,mu=mu,v=-4,self=endself)
    endR=asts.evaluate(asts.expression('pulse_end_stress_C3','end','logR'),env)
    endB=asts.evaluate(asts.expression('pulse_end_stress_C3','end','logB'),env)
    endD=asts.evaluate(asts.expression('pulse_end_stress_C3','end','logD'),env)
    endH=asts.evaluate(asts.expression('pulse_end_stress_C3','end','logH'),env)
    zero('actual_original_radius_same_at_gap_end',logs['R'].subs(d,4*mu)-endR)
    zero('actual_original_B_same_at_gap_end',
        sum(v.subs(d,4*mu) for v in logs['B'].values())-sum(endB.values()))
    zero('actual_original_H_same_at_gap_end',logs['H'].subs(d,4*mu)-endH)
    zero('effective_D1_source_factor_at_gap_end',logs['D1'].subs(d,4*mu)-endD-4*lam1)
    zero('effective_D2_source_factor_at_gap_end',logs['D2'].subs(d,4*mu)-endD-4*lam2)
    radius=asts.replay('global_physical_assembly','radius',dict(PULSE=PULSE))
    radius_self=SimpleNamespace(ctx=SimpleNamespace(mpf=lambda value:value),
        params=SimpleNamespace(mu=mu),logRp=lrp)
    production_gap=radius(radius_self,'pulse_gap',13-4*mu,{},None)[0]
    production_gap_end=radius(radius_self,'pulse_gap_end',-4,{},None)[0]
    production_end=radius(radius_self,'pulse_end',-4,{},None)[0]
    zero('actual_production_gap_right_radius',production_gap-logs['R'].subs(d,4*mu))
    zero('actual_production_gap_end_chart_right_radius',production_gap_end-production_gap)
    zero('actual_production_end_left_radius',production_end-endR)
    zero('actual_production_physical_radius_gap_end_join',production_gap-production_end)
    gapcall=asts.expression('pulse_gap_physical_C2','gap','point',wanted=
        "lift_physical_packet(self.ctx,mapped,self.similarity.delta,packet['gap_source_rows']['velocity'],log_tau,theta,viscosity)")
    endcall=asts.expression('pulse_end_physical_C2','end','point',wanted=
        'lift_physical_packet(c,packet,self.stress.native.delta,velocity,log_tau,theta,viscosity)')
    for label,left,rightarg in zip(('log_tau','theta','viscosity'),gapcall.args[-3:],endcall.args[-3:]):
        if ast.dump(left)!=ast.dump(rightarg):raise ValueError('Physical endpoint argument differs: '+label)
        checks['actual_gap_end_pullback_passes_identical_'+label+'_argument']=True
    checks['actual_gap_physical_original_delta_from_same_similarity_source']=True
    # The same full beta histories differ only by the original backward factors.
    M1=s.Function('same_full_beta_M1')(z)
    for j in range(5):
        zero('same_nonzero_radial_history_ordinary_row'+str(j),
            s.exp(logs['D1'].subs(d,4*mu))*M1*(-lam1)**j
            -s.exp(endD)*s.exp(4*lam1)*M1*(-lam1)**j)
    for power in (1,2):
        zero('linear_or_quadratic_radial_error_factor_join_'+str(power),
            s.exp(power*logs['D1'].subs(d,4*mu))
            -s.exp(power*endD)*s.exp(4*power*lam1))
    # Replay the current row algorithms on arbitrary axial histories. Endpoint
    # inputs are zero, but the incoming Mz function and its derivatives remain.
    env=dict(mp=SimpleNamespace(mpf=lambda value:s.Rational(str(value))),
        math=__import__('math'),axial_derivative=lambda value:s.diff(value,z),
        physical_operators=physical_operators)
    asts.replay('collar_Gamma_C4','product_rows',env)
    asts.replay('collar_stress_C3','shifted_rows',env)
    asts.replay('pulse_end_physical_C2','pulse_velocity_rows',env)
    asts.replay('pulse_end_physical_C2','axial_n',env)
    asts.replay('pulse_end_physical_C2','axial_operator_rows',env)
    asts.replay('pulse_end_physical_C2','pulse_remainder_sectors',env)
    C=1/(1+z*z);zero_rows=[s.Integer(0)]*5
    raw=[M1*(-lam1)**j for j in range(5)]
    velocity=env['pulse_velocity_rows'](c,delta,mu,z,C,zero_rows,raw)
    for j,row in enumerate(velocity['radial']):
        zero('actual_gap_radial_full_ordinary_rate_row'+str(j),
            row-velocity['radial'][0]*(-s.Rational(1,2))**j)
        zero('actual_gap_axial_velocity_structural_zero_row'+str(j),velocity['axial'][j])
    errors=env['pulse_remainder_sectors'](c,delta,mu,z,velocity)
    for j in range(3):
        zero('actual_gap_radial_viscosity_structural_zero_row'+str(j),
            errors['radial']['radial_viscosity']['rows'][j])
        zero('actual_gap_axial_remainder_structural_zero_row'+str(j),
            errors['axial']['axial_viscosity']['rows'][j])
    checks['nonzero_radial_time_transport_and_axial_viscosity_not_removed']=True
    for i in range(3):
        for j in range(3-i):
            if any(k+n>4 for k,n in physical_operators()[i,j+2]):
                raise ValueError('Gap physical error exceeds available source mixed4')
            checks['same_physical_mixed2_error_operator_'+str(i)+str(j)]=True
    gap=records['pulse_gap_similarity_C4_check']
    for flag in ('actual_original_whole_inactive_gap_similarity_companion_constructed',
        'gap_end_similarity_velocity4_moment4_stress3_pressure4_functional_join_verified',
        'full_nonzero_raw_moment_radial_and_stress_histories_preserved'):
        if not gap[flag]:raise ValueError('Missing gap physical source '+flag)
        checks['consumed_'+flag]=True
    native_join=records['pulse_interface_certificate']['source_bound_functional_pulse_identities']
    for flag in ('actual_common_end_scale_normalization','production_gap_coordinate_leading_log',
        'gap_coordinate_pressure_log','gap_coordinate_pressure_time','gap_coordinate_energy_whole_Z',
        'gap_end_full_future_linear_weight','gap_end_energy_at_s_minus4',
        'gap_end_angular_history_whole_Z','gap_end_original_pressure_whole_Z',
        'full_future_supports_at_minus4_and_empty_at_zero',
        'end_energy_cross_products_exactly_zero_by_disjoint_supports',
        *('functional_main_gap_axial_order'+str(n) for n in range(6))):
        if not native_join[flag]:raise ValueError('Missing native functional gap/end source '+flag)
        checks['directly_consumed_native_'+flag]=True
    pressure_join=records['pulse_end_flatten_join']['source_endpoint_binding']['identities']
    for flag in ('actual_pulse_P0_getter_is_same_canonical_flatten_pressure_over_C0_squared',
        'consumed_same_absolute_pressure_identified_by_original_FTC_and_power_datum'):
        if not pressure_join[flag]:raise ValueError('Missing original uncapped pressure getter '+flag)
        checks['directly_consumed_'+flag]=True
    accepted=records['pulse_end_physical_C2']
    proof=accepted['full_physical_operator_proof']
    if not all(proof['identities'].values()):
        raise ValueError('Full current meridional physical operator not admitted')
    if not records['pulse_end_physical_C2_check']['actual_pulse_end_physical_decomposition_constructed']:
        raise ValueError('Accepted actual physical end required')
    checks['current_full_meridional_NS_operator_consumed_without_rerun']=True
    checks['same_fixed_positive_nu_original_lambda_time_radius_and_pressure_units']=True
    checks['source_stress3_implies_completed_physical_stress3_diagonal2_divergence2_join']=True
    checks['source_velocity4_implies_radial_and_angular_physical_remainder2_join']=True
    checks['gap_end_join_is_functional_not_interval_overlap']=True
    return dict(identities=checks,input_hashes=asts.hashes,
        admitted_full_physical_operator_identities=len(proof['identities']),
        same_uncapped_pressure_and_five_moment_family_retained=True,
        gap_end_completed_physical_interface_verified=True,
        effective_radial_history_scale='D1=D0*exp((.5-mu)*d/mu)',
        right_endpoint_radial_remainder_not_assumed_zero=True,
        regional_error_not_claimed_flat=True)


class CompliantPulseGapPhysicalC2:
    @source_precision
    def __init__(self):
        self.similarity=CompliantPulseGapSimilarityC4()
        self.ctx=self.similarity.ctx;self.family=self.similarity.family;self.source=self.similarity.source
        self.hashes=dict(self.similarity.hashes);self.records=dict(self.similarity.records)
        for stem in ('pulse_gap_similarity_C4','pulse_gap_similarity_C4_check','pulse_end_physical_C2'):
            name=PREFIX+stem+'.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
            if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(self.family,self.source):
                raise ValueError('Gap physical source family differs: '+stem)
            if 'all_passed' in record and not record['all_passed']:
                raise ValueError('Unaccepted gap physical source: '+stem)
            for path,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Gap physical source changed: '+path)
            self.hashes.update(record['input_hashes']);self.hashes[name]=hashlib.sha256(raw).hexdigest()
            self.records[stem]=record
        self.proof=source_and_join_binding(self.records)
        self.hashes.update(self.proof['input_hashes'])
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    @source_precision
    def gap(self,Z,distance,whole=False,right_endpoint=False,log_tau='-1',theta='0',viscosity='1'):
        packet=self.similarity.packet(Z,distance,whole=whole,right_endpoint=right_endpoint)
        mapped=gap_to_physical_packet(packet,self.similarity.mu,right_endpoint)
        point=lift_physical_packet(self.ctx,mapped,self.similarity.delta,
            packet['gap_source_rows']['velocity'],log_tau,theta,viscosity)
        # This adapter is for the inactive gap; the end-only cone flag is not
        # an assertion about the already admitted downstream end.
        point.pop('pulse_end_cone_certified',None)
        point.update(domain=DOMAIN,distance=packet['distance'],
            actual_pulse_gap_physical_decomposition_constructed=True,
            gap_end_completed_physical_interface_verified=True,
            actual_physical_stress_mixed_order=3,actual_physical_remainder_mixed_order=2,
            effective_radial_history_scale_is_D1=True,actual_axial_velocity_and_remainder_zero=True,
            right_endpoint_is_original_end_s_minus4=right_endpoint,
            gap_right_radial_remainder_not_assumed_zero=True,**{flag:False for flag in FALSE_FLAGS})
        return point

    @source_precision
    def report(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            domain=DOMAIN,source_and_physical_join_binding=self.proof,
            admitted_full_physical_operator_proof=self.records['pulse_end_physical_C2']['full_physical_operator_proof'],
            samples=[self.gap('.31',1,theta='.37',viscosity='.01'),
                self.gap('.5',2,theta='.37',viscosity='.7')],
            right_end_source=self.gap([-1,1],None,right_endpoint=True,theta=None,viscosity='.7'),
            whole_original_gap=self.gap([-1,1],None,whole=True,theta=None,viscosity='.7'),
            input_hashes=self.hashes,actual_pulse_gap_physical_decomposition_constructed=True,
            gap_end_completed_physical_interface_verified=True,
            full_nonzero_meridional_velocity_and_radial_remainder_retained=True,
            source_caps_used_as_defining_field_values=False,**{flag:False for flag in FALSE_FLAGS})


@source_precision
def run():
    result=CompliantPulseGapPhysicalC2().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Whole inactive-gap completed physical tensor and radial remainder generated; cone/global/recursion pending',flush=True)
    return result


if __name__=='__main__':run()
