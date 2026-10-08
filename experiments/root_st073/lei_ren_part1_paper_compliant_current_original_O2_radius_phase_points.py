"""Original O2 radius phase points, with a retained positive origin offset.

Exact integer periods are removed without expanding the selected logCstar.
Actual defining exp/log arithmetic gives the point approximation. Directed
arithmetic and the checked original hb*s_c/2 bound enclose the true phase.
No saved radius cap is a point value and no loop scales/global N are selected.
"""
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_factored_point_inputs as point
import lei_ren_part1_paper_compliant_current_native_spatial_phase as phase

HERE,PREFIX,sha=point.HERE,point.PREFIX,point.sha
NAME=PREFIX+'current_original_O2_radius_phase_points.json'
RECEIPT=PREFIX+'current_original_O2_radius_phase_points_check.json'
GATE='original_O2_true_radius_phase_points_with_directed_modular_and_positive_origin_offset_errors'


def interval_record(value):
    lo,hi=point.endpoints(value)
    return dict(lower=lo,upper=hi)


class OriginalO2RadiusPhasePoints:
    mode='original_O2_radius_phase_point_with_directed_error_budget'
    def __init__(self,dps=80):
        if type(dps) is not int or dps<40:raise ValueError('At least 40 phase digits required')
        self.dps=dps;self.frame=point.source.OriginalO2SourceParameterFrame(dps);self.family=self.frame.family
        definitions=self.frame.definitions
        if definitions['Md']!=40 or s.simplify(definitions['logPstar']-(s.exp(40)+11))!=0:
            raise ValueError('Phase exp/log recipe must equal the original parameter frame')
        if s.expand(self.frame.radius_log(0)-s.log(110)-10*(definitions['logCstar']+definitions['logPstar']))!=0:
            raise ValueError('Phase reference radius must equal the original parameter frame')
        self.hashes={};self.cache={}
        for stem,gate in (('current_original_O2_factored_point_inputs',point.GATE),
            ('current_native_spatial_phase',phase.GATE),
            ('current_generic_shear_loop_domain','current_original_generic_loop_two_sided_collars_and_repair_geometry_certified'),
            ('current_inner_exit_strict_collar','current_inner_exit_strict_collar_attached_to_current_source_graph_certified')):
            name=PREFIX+stem+'_check.json';receipt=json.loads((HERE/name).read_bytes())
            family=receipt.get('source_family') or {key:receipt[key] for key in self.family}
            if not receipt.get('all_passed') or not receipt.get(gate) or family!=self.family:
                raise ValueError('Same checked original O2/radius/phase family required: '+stem)
            for filename,digest in {**receipt['input_hashes'],name:sha(name)}.items():
                if filename in self.hashes and self.hashes[filename]!=digest:
                    raise ValueError('Original phase dependencies disagree: '+filename)
                if sha(filename)!=digest:raise ValueError('Original phase dependency changed: '+filename)
                self.hashes[filename]=digest
        self.domain=json.loads((HERE/phase.DOMAIN).read_bytes())['current_original_generic_loop_domain']
        self.collar=json.loads((HERE/phase.COLLAR).read_bytes())
        self.proof=self.collar['explicit_current_inner_exit_strict_collar']
        formal=self.domain['formal_radii']
        if formal['Ra']!='4*epsilon_core' or formal['r_minus']!='Ra*exp(hb*s_c/2)':
            raise ValueError('Original source phase origin changed')
        attachment=self.collar['exact_current_exit_source_attachment']
        if not attachment['passed'] or attachment['exact_width']!='hb=epsilon_b=cstar*K^-100; K is the original global physical norm sum':
            raise ValueError('Original positive hb source width required')
        self.sc_bits=self.proof['selected_first_phase_endpoint']
        if self.sc_bits['lower_exact_mpf_tuple']!=self.sc_bits['upper_exact_mpf_tuple']:
            raise ValueError('Original explicitly selected s_c, not a point chosen from its enclosure, required')
        self.sc_tuple=tuple(self.sc_bits['lower_exact_mpf_tuple'])
        if self.sc_tuple[0] or not self.sc_tuple[1]:raise ValueError('Positive original selected s_c required')
        self.loghb=self.proof['exact_source_width_log_enclosure']
        if self.domain['source_family']!=self.family:raise ValueError('Original r_minus family differs')
        sign,man,exponent,bits=self.frame.selected_logCstar_mpf_tuple
        if sign or exponent<0:raise ValueError('Integer original selected logCstar is required for this reduction')
        for name in (phase.DOMAIN,phase.COLLAR,Path(__file__).name):self.hashes[name]=sha(name)

    def evaluate(self,*,y,N):
        y=point.source.exact_coordinate(y)
        if type(N) is not int or N<1 or N.bit_length()>4096:
            raise ValueError('Explicit positive integer candidate N with at most 4096 bits required')
        key=(y,N)
        if key in self.cache:return self.cache[key]
        digits=(N.bit_length()*30103+99999)//100000
        c=mp.mp.clone();c.dps=self.dps+digits
        iv=MPIntervalContext();iv.dps=max(260,c.dps+40)
        with mp.workdps(iv.dps+40):
            chosen=c.make_mpf(self.frame.selected_logCstar_mpf_tuple)
            unused,integer_proof=phase.binary_mod_one(iv,iv.mpf(chosen),Fraction(10*N))
            if integer_proof['exact_fraction']!={'numerator':0,'denominator':1}:
                raise ArithmeticError('Selected logCstar integer period was not removed exactly')
            iy=iv.mpf(int(y.p))/int(y.q);cy=c.mpf(int(y.p))/int(y.q)
            main=c.mpf(N)*(14*(c.exp(40)+11)+c.ln(c.mpf(110)/4)+cy)
            main_box=iv.mpf(N)*(14*(iv.exp(40)+11)+iv.ln(iv.mpf(110)/4)+iy)
            p=mp.mp.clone();p.dps=iv.dps+30
            period=int(p.floor(p.mpf(point.endpoints(main_box)[0])))
            reduced=main_box-iv.mpf(period)
            # hb is the actual source width. Its accepted directed enclosure
            # is used as an error bound only, never as an evaluated point.
            hlo=p.make_mpf(tuple(self.loghb['lower_exact_mpf_tuple']))
            hhi=p.make_mpf(tuple(self.loghb['upper_exact_mpf_tuple']))
            sc=p.make_mpf(self.sc_tuple)
            offset_log=iv.ln(N)+iv.mpf([hlo,hhi])+iv.ln(iv.mpf(sc))-iv.ln(2)
            offset_log_upper=point.endpoints(offset_log)[1]
            buffer=iv.mpf(10)**-(self.dps+20)
            if offset_log_upper>=point.endpoints(iv.ln(buffer))[0]:
                raise ArithmeticError('Original positive microscopic origin offset exceeds phase budget')
            budget_upper=point.endpoints(buffer)[1]
            # This positive numerical error budget bounds N*hb*s_c/2. It
            # does not replace the original origin in the source expression.
            corrected=reduced-iv.mpf([0,budget_upper])
            periodic=phase.ordinary_mod_one(iv,corrected)
            approximate=main-c.floor(main)
            finite_error=point.endpoints(abs(main_box-iv.mpf(main)))[1]
        result=dict(source_family=self.family,mode=self.mode,original_y_exact=str(y),explicit_candidate_N=N,
            approximate_original_fractional_phase=approximate,
            true_original_phase_directed_boxes=[interval_record(box) for box in periodic['boxes']],
            phase_arithmetic_absolute_error_upper=finite_error,
            positive_original_origin_offset=dict(source='N*hb*s_c/2',strictly_positive=True,
                sign_in_original_affine_phase=-1,error_log_upper=offset_log_upper,
                representable_absolute_error_budget_upper=budget_upper,
                error_budget_dominates_offset_by_directed_log_comparison=True,
                actual_hb_sc_product_not_materialized_or_set_to_zero=True),
            periodic_projection_full_period=periodic['full_period'],
            error_contract='circular distance(true phase, approximate phase) <= arithmetic_error + exp(origin_offset_error_log_upper)',
            exact_source_phase='frac(N*(log(110/4)+14*logPstar+10*logCstar+1000+y-hb*s_c/2))',
            original_phase_origin=dict(Ra='4*epsilon_core',r_minus='Ra*exp(hb*s_c/2)',
                log_r_minus='logRa+hb*s_c/2',log_positive_log_radius_offset='loghb+logsc-log2',
                exact_width='hb=cstar*K^-100, same actual global physical norm sum',
                selected_sc_exact_mpf_tuple=list(self.sc_tuple)),
            exact_removed_integer_periods=dict(logCstar=integer_proof,constant_1000N=1000*N),
            same_original_radius_phase_bound=True,source_radius_caps_not_consumed=True,
            huge_selected_integer_and_native_radius_not_materialized=True,
            candidate_N_not_global_selection=True,conditioned_loop_inverse_or_primitives_installed=False,
            conditional_on_explicit_candidate_N_and_proved_origin_offset_budget=True,
            numerical_original_source_point_or_integral_oracle_installed=False)
        self.cache[key]=result;return result


def run():
    began=time.monotonic();owner=OriginalO2RadiusPhasePoints()
    samples=[owner.evaluate(y=y,N=N) for y,N in (('0',1),('.53',7),('1',257),('.53',(1<<200)+3))]
    report=dict(**{GATE:True},source_family=owner.family,mode=owner.mode,
        original_O2_radius_phase_queries=samples,
        source_phase_origin_bound_and_integer_periods_exactly_removed=True,
        actual_defining_exp_log_phase_points_with_directed_arithmetic_errors=True,
        positive_original_origin_offset_retained_and_rigorously_bounded=True,
        original_logPstar_logRref_recipe_checked_against_common_parameter_frame=True,
        phase_service_conditional_on_candidate_N_and_proved_offset_budget=True,
        conditioned_native_factor_and_loop_evaluation_installed=False,
        numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='True O2 frac(N log(R/r_minus)) point approximations and directed periodic covers at explicit candidate integer N. Original positive microscopic origin offset is retained with error. No conditioned loop inverse, A/B primitives, signed density integral, full oracle, global N or recursion is installed.')
    (HERE/NAME).write_text(json.dumps(point.source.inertial.profiles.loop.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Original O2 true radius phase points and directed positive-origin error budgets connected',flush=True)
    return report


if __name__=='__main__':run()
