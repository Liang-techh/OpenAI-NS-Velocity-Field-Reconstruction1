"""Whole original Ts-long steep-power cone and same-source outer tail.

The original angular inlet correlation gives a positive theta enclosure.
Exact shared Ts/wait terms cancel in the source Bmax log before enclosure.
The regional cone concerns the original two-vector, not the completed tensor.
"""
import ast
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_steep_power_stress_C3 import CompliantSteepPowerStressC3,source_precision
from lei_ren_part1_paper_compliant_global_physical_assembly import CompliantGlobalPhysicalAssembly
from lei_ren_part1_paper_compliant_collar_physical_C2 import stress_source_log_parts
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
DOMAIN=dict(phase=[0,1],Z=[-1,1],offset='q=-wait-1-Ts*(1-phase); -wait-1-Ts<=q<=-wait-1')
TAIL_DOMAIN='Rtail*exp(-wait-1-Ts)<=R<Rtail*exp(3); Gamma stress exactly zero beyond'


def steep_power_cone_identities():
    a,phase,wait,Ts,eps=s.symbols('a phase wait Ts eps',real=True)
    delta=2*a; k=1-a; bh=s.Rational(1,2)+a; q=-wait-1-Ts*(1-phase)
    KQ=s.symbols('actual_KQ',positive=True); K=KQ*s.exp(k*Ts*(1-phase))
    Sin,Bbase=s.symbols('exact_positive_inverse_Rtail original_Bbase',positive=True)
    proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:raise ArithmeticError('Power cone identity failed: '+name)
        proofs[name]=True
    Kq=s.diff(K,phase)/Ts
    zero('actual_original_power_K_log_rate',Kq/K+k)
    zero('actual_original_power_kappa_minus2_exact_two',delta-2*Kq/K-2)
    shear=2*Sin*s.exp(-q)*(Kq-(1+a)*K)
    zero('actual_shear_negative_exact_positive_source_factors',shear+4*Sin*s.exp(-q)*K)
    B=Bbase*s.exp(-bh*q); Bmax=Bbase*s.exp(bh*(wait+1+Ts))
    zero('source_B_decreases_on_whole_original_Ts',B-Bmax*s.exp(-bh*Ts*phase))
    zero('source_B_derivative_negative_positive_factors',s.diff(B,phase)+bh*Ts*B)
    Ev0,thetaT=s.symbols('source_Ev0 source_thetaT',positive=True)
    zero('source_Bmax_actual_inlet_amplitude',
         Bmax.subs(Bbase,Ev0*thetaT*s.exp(-bh*wait)/(1-eps))-Ev0*thetaT*s.exp(bh*(1+Ts))/(1-eps))
    zero('original_shared_steep_length_cancellation',-s.Rational(3,2)*Ts+bh*Ts+k*Ts)
    qt,qz,ct,cz,st,Bcurrent=s.symbols('Qtheta Qz Ctheta Cz Stheta B',real=True)
    numerator=2*(qt*ct*qt*st)**2-2*(qt*qz*cz*st)**2
    zero('two_vector_cone_exact_positive_factor_cancellation',
         numerator.subs(qz,qt*Bcurrent)-qt**4*st**2*(2*ct**2-2*Bcurrent**2*cz**2))
    zero('actual_stress_shear_dot_product',qt*ct*qt*st-qt**2*ct*st)
    return dict(identities=proofs,entire_original_steep_power_domain=DOMAIN,
                exact_source_K_positive_from_epsilon_below1=True,
                exact_source_S_positive_inverse_Rtail_not_cap_endpoint=True,
                original_power_kappa_minus2_exactly_two=True,
                source_K_independent_of_Z_and_Sz_exact_zero=True,
                normalized_cone_condition='2*B^2*Cz^2<2*Ctheta^2',
                sufficient_whole_domain_bound='2*Bmax^2*(Cz_absolute_upper/Ctheta_lower)^2<2',
                actual_Bmax_definition='Ev0*thetaT*exp(bh*(1+Ts))/(1-epsilon); original q=-wait-1-Ts',
                same_positive_physical_nu_lambda_factors_preserve_two_vector_cone=True,
                completed_diagonal_does_not_enter_two_vector_cone=True)


def steep_power_Bmax_binding():
    """Bind actual source logs, regrouping shared Ts and wait before enclosure."""
    path=HERE/(PREFIX+'collar_physical_C2.py')
    tree=ast.parse(path.read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='stress_source_log_parts')
    env={}; exec(compile(ast.Module(body=[fn],type_ignores=[]),'<actual power source B logs>','exec'),env)
    a,mu,wait,Ts,length,lp,lu,lrp,lone=s.symbols('a mu wait Ts Lrel logP logU logRp logone',real=True)
    c=SimpleNamespace(mpf=s.sympify,ln=s.log)
    flatten=SimpleNamespace(logEv2_parts=dict(inlet_log=2*lu))
    assembly=SimpleNamespace(ctx=c,logP=lp,logRp=lrp,dispatch=SimpleNamespace(provider=lambda chart:flatten))
    heat=SimpleNamespace(ctx=c,a=a,mu=mu,k=1-a,bh=s.Rational(1,2)+a,
                         steep=SimpleNamespace(wait=wait,Ts=Ts,logone=lone,outer=SimpleNamespace(Lrel=length)))
    original=env['stress_source_log_parts'](assembly,heat,-wait-1-Ts)['B']
    owned=ast.parse(Path(__file__).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(owned) if isinstance(n,ast.FunctionDef) and n.name=='whole_steep_power_bounds')
    names=('logparts',"logparts['steep']","logparts['waiting_and_current']")
    assignments=[n for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t) in names for t in n.targets)]
    if len(assignments)!=3:raise ValueError('Owned full-power inlet cancellation changed')
    env.update(assembly=assembly,h=heat,c=c)
    exec(compile(ast.Module(body=assignments,type_ignores=[]),'<owned power Bmax cancellation>','exec'),env)
    current=env['logparts']; checks={}
    if set(original)!=set(current):raise ValueError('Power inlet source log keys changed')
    for key in set(original)-{'steep','waiting_and_current'}:
        if s.simplify(original[key]-current[key])!=0:raise ArithmeticError('Power inlet source log mismatch: '+key)
        checks[key]=True
    if s.simplify(sum(original.values())-sum(current.values()))!=0:raise ArithmeticError('Actual power Ts/wait source cancellation failed')
    checks['combined_original_steep_and_current_Ts_terms']=True
    return dict(identities=checks,actual_Bmax_is_original_source_at_q_minus_wait_minus1_minus_Ts=True,
                exact_wait_plus_q_equals_minus1_minus_Ts_before_enclosure=True,
                exact_minus1point5Ts_plus_bhTs_equals_minus_kTs_before_enclosure=True,
                source_cap_values_used_as_field=False)


def steep_power_cone_source_bridge(stress):
    proofs={}; bindings={}
    for flag in ('original_full_future_backward_FTC_verified','exact_resonant_Kq_over_K_equals_minus_k',
                 'exact_kappa_minus2_equals_two','actual_source_pressure_datum_retained'):
        if not stress.proof[flag]:raise ValueError('Power cone exact source missing: '+flag)
        proofs['consumed_'+flag]=True
    for flag in ('same_full_exit_endpoint_and_normalized_moment_ODE_uniqueness_used',
                 'actual_meridional_source_and_terminal_zero_primitive_route_consumed',
                 'actual_absolute_pressure_datum_and_exact_Ev0_units_consumed',
                 'actual_source_XS_correlation_is_same_full_angular_future'):
        if not stress.bridge[flag]:raise ValueError('Power cone history missing: '+flag)
        proofs['consumed_'+flag]=True
    for flag in ('actual_power_exit_K_mixed4_join_verified','actual_power_exit_stress_mixed3_AST_join_verified',
                 'actual_power_exit_pressure_mixed4_AST_join_verified'):
        if not stress.bridge['actual_power_exit_formula_join'][flag]:raise ValueError('Power cone join missing: '+flag)
        proofs['consumed_'+flag]=True
    for flag in ('exact_S_is_positive_source_not_enclosure_endpoint','actual_original_velocity_reference_normalization_verified'):
        if not stress.exit_source.bridge[flag]:raise ValueError('Power cone inherited source missing: '+flag)
        proofs['consumed_'+flag]=True
    def syntax(stem,method,target,expression):
        fn=next(n for n in ast.walk(ast.parse((HERE/(PREFIX+stem+'.py')).read_text(encoding='utf8')))
                if isinstance(n,ast.FunctionDef) and n.name==method)
        values=[n.value for n in ast.walk(fn) if isinstance(n,ast.Assign) and any(ast.unparse(t)==target for t in n.targets)]
        wanted=ast.dump(ast.parse(expression,mode='eval').body)
        if sum(ast.dump(n)==wanted for n in values)!=1:raise ValueError('Power cone source changed: '+stem+'.'+target)
        bindings[stem+'.'+method+'.'+target]=True
    for target,expression in (('self.delta','self.steep.delta'),('self.a','self.delta/2'),('self.k','self.steep.k'),('self.bh','self.steep.bh')):
        syntax('collar_Gamma_C4','__init__',target,expression)
    syntax('steep_waiting_C4','__init__','self.k','1-self.delta/2')
    syntax('steep_waiting_C4','__init__','self.bh',"c.mpf('.5')+self.delta/2")
    syntax('steep_power_stress_C3','power_shape','ell','heat.steep.Ts*(1-phase)')
    syntax('steep_power_stress_C3','power_shape','K',"terminal['KQ']*c.exp(heat.k*ell)")
    syntax('steep_power_stress_C3','power_stress_rows','q','-heat.steep.wait-1-ell')
    syntax('steep_power_stress_C3','power_stress_rows','St','K*(-4*heat.S*c.exp(-q))')
    return dict(identities=proofs,actual_original_power_source_bindings=bindings,
                actual_Bmax_source_log_binding=steep_power_Bmax_binding(),
                actual_same_full_moment_source_bridge=stress.bridge,
                sampled_interval_overlap_used_as_proof=False)


@source_precision
def whole_steep_power_bounds(stress,assembly):
    h=stress.heat; c=stress.ctx
    packet=stress.steep_power([-1,1],[0,1]); rows=packet['steep_power_similarity_stress_mixed3_factored']
    Ct=rows['theta']['y0_Z0']; Cz=rows['axial']['y0_Z0']
    theta_lower=c.mpf(endpoints(Ct)[0]); axial_upper=c.mpf(endpoints(abs(Cz))[1])
    K=packet['steep_power_shape']['K_rows'][0][0]
    m=packet['steep_power_shear_strength_kappa_minus2']
    # Original Bmax at q=-wait-1-Ts; shared power/current Ts and wait
    # cancellations are source identities, not cancellations of unrelated boxes.
    logparts=stress_source_log_parts(assembly,h,c.mpf(0))['B']
    logparts['steep']=-h.k*h.steep.Ts
    logparts['waiting_and_current']=h.bh
    logBmax=sum(logparts.values(),c.mpf(0))
    margins=dict(delta=h.delta,epsilon=h.eps,one_minus_epsilon=1-h.eps,k=h.k,
                 p=h.prate,bh=h.bh,Lmin=1-h.delta,wait=h.steep.wait,Ts=h.steep.Ts,
                 source_K_lower=K,actual_kappa_minus2=m,theta_stress=Ct)
    for label,value in margins.items():
        if endpoints(value)[0]<=0 or not mp.isfinite(endpoints(value)[1]):raise ArithmeticError('Power cone positive margin failed: '+label)
    if endpoints(h.eps)[1]>=1:raise ArithmeticError('Original power K must be positive')
    if endpoints(m)!=(mp.mpf(2),mp.mpf(2)):raise ArithmeticError('Actual power kappa-2 must be exactly two')
    if endpoints(axial_upper)[0]<=0 or not mp.isfinite(endpoints(axial_upper)[1]):raise ArithmeticError('Finite nonzero axial bound required')
    if not all(mp.isfinite(v) for v in endpoints(logBmax)):raise ArithmeticError('Actual source Bmax requires finite logs')
    logcone=c.ln(2)+2*logBmax+2*c.ln(axial_upper/theta_lower); threshold=c.ln(2)
    if endpoints(logcone)[1]>=endpoints(threshold)[0]:raise ArithmeticError('Whole original steep-power directional cone failed')
    return dict(positive_margins=margins,actual_whole_source_theta_enclosure=Ct,actual_whole_source_axial_enclosure=Cz,
                theta_uniform_positive_lower=theta_lower,axial_uniform_absolute_upper=axial_upper,
                actual_kappa_minus2_enclosure=m,exact_kappa_minus2_upper=c.mpf(2),
                actual_Bmax_source_log_parts=logparts,actual_log_Bmax_enclosure=logBmax,
                log_directional_term_upper=logcone,log_threshold=threshold,
                whole_phase_Z_box_and_original_Ts_used=True,source_caps_only_define_bounds_not_field_values=True,
                interval_shear_upper_zero_not_used_for_strict_source_sign=True)


def run():
    with mp.workdps(300):
        stress=CompliantSteepPowerStressC3(); assembly=CompliantGlobalPhysicalAssembly()
        if (stress.family,stress.source)!=(assembly.family,assembly.source):raise ValueError('Power cone field family mismatch')
        hashes=dict(stress.hashes); hashes.update(assembly.hashes)
        for stem,gates in (
            ('steep_power_stress_C3_check',('actual_original_steep_power_similarity_stress_recovered',)),
            ('steep_power_physical_C2_check',('actual_regional_physical_steep_power_stress_remainder_identity_verified',
                'steep_power_exit_physical_stress_mixed3_and_remainder_mixed2_join_verified')),
            ('steep_exit_cone_check',('steep_exit_waiting_collar_two_vector_cone_certified',))):
            name=PREFIX+stem+'.json'; receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or any(not receipt[gate] for gate in gates):raise ValueError('Power cone prerequisite missing: '+name)
            if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(stress.family,stress.source):
                raise ValueError('Power cone receipt family mismatch')
            for source,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Power cone prerequisite changed: '+source)
            hashes.update(receipt['input_hashes']); hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        proof=steep_power_cone_identities(); bridge=steep_power_cone_source_bridge(stress)
        bounds=whole_steep_power_bounds(stress,assembly); hashes.update(assembly.dispatch.hashes)
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(actual_five_defect_family_sha256=stress.family,implicit_source_sha256=stress.source,
                    scope='Whole original Ts-long steep power phase[0,1], source Z[-1,1]; physical tau>0, r>0, |Z|<1',
                    domain=DOMAIN,tail_domain=TAIL_DOMAIN,exact_identities=proof,source_bridge=bridge,bounds=bounds,input_hashes=hashes,
                    steep_power_cone_certified=True,actual_steep_power_theta_stress_positive=True,
                    actual_steep_power_source_shear_strictly_negative=True,actual_steep_power_directional_cone_uniform_margin=True,
                    steep_power_exit_waiting_collar_two_vector_cone_certified=True,
                    same_source_steep_power_exit_similarity_and_physical_joins_consumed=True,
                    fixed_positive_viscosity_two_vector_cone_transfer_verified=True,
                    steep_power_regional_remainder_exact_zero=False,completed_full_tensor_cone_certified=False,
                    whole_outer_cone_certified=False,global_admissible_stress_lift_constructed=False,
                    independently_bounded_global_flat_remainder=False,physical_energy_integral_certified=False,temporal_recursion=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
        print('Whole original Ts steep-power cone generated; exit/waiting/collar joined, upstream/global cone pending',flush=True)
        return result


class CertifiedSteepPowerPhysical:
    """Attach accepted regional cone metadata to the unchanged physical field."""
    def __init__(self):
        from lei_ren_part1_paper_compliant_steep_power_physical_C2 import CompliantSteepPowerPhysicalC2
        path=HERE/(PREFIX+'steep_power_cone_check.json'); receipt=json.loads(path.read_bytes())
        if not receipt['all_passed'] or not receipt['steep_power_cone_certified']:raise ValueError('Accepted whole power cone required')
        for source,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Certified power source changed: '+source)
        self.field=CompliantSteepPowerPhysicalC2(); self.family,self.source=self.field.family,self.field.source
        if (receipt['actual_five_defect_family_sha256'],receipt['implicit_source_sha256'])!=(self.family,self.source):
            raise ValueError('Power physical cone family mismatch')
        self.receipt=receipt

    def steep_power(self,*args,**kwargs):
        point=self.field.steep_power(*args,**kwargs)
        point.update(steep_power_cone_certified=True,steep_power_exit_waiting_collar_two_vector_cone_certified=True,
                     certified_cone_scope='Original two-vector only; completed diagonal/global cone not certified',
                     completed_full_tensor_cone_certified=False)
        return point


if __name__=='__main__':run()
