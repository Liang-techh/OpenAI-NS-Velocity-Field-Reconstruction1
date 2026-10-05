"""Continuous original inactive-gap regional two-vector stress cone.

All four angular corrections and four axial histories use their actual
whole-gap source factors. Axial shear is structurally zero only here;
nonzero radial and cumulative histories remain. Global tensor admissibility
and temporal recursion are separate unfinished requirements.
"""
import ast
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_gap_similarity_C4 import (
    sources as gap_sources, DOMAIN, PREFIX, source_precision)
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import SourceAST
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

HERE=Path(__file__).parent
TAIL_DOMAIN='Rtail*exp(-2/mu-wait-Ts-102-Lrel)<=R<Rtail*exp(3); Gamma stress exactly zero beyond'
ADMISSIONS=('pulse_gap_cone_certified','actual_gap_theta_stress_positive',
    'actual_gap_source_shear_strictly_negative','actual_gap_zero_axial_shear_kappa_verified',
    'actual_gap_continuous_directional_cone_verified',
    'gap_end_flatten_power_angular_entry_exit_waiting_collar_two_vector_cone_certified',
    'same_source_gap_end_completed_physical_interface_consumed',
    'fixed_positive_viscosity_two_vector_cone_transfer_verified')
FALSE_FLAGS=('completed_full_tensor_cone_certified','whole_outer_cone_certified',
    'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder',
    'physical_energy_integral_certified','full_background_NS_validation','temporal_recursion')
# Relative to Qbase=sqrt(R/2)*B; tuple=(R,B,D,H,selected D,pressure Q).
FACTORS={
 'theta':{
    'signed_original_memory':(0,0,0,1,'D0',False),
    'meridional_Mz_history':(0,1,1,0,'D1',False),
    'mixed_angular_axial_history':(0,1,1,0,'D2',False),
    'radial_shear':(-1,0,0,0,'D0',False)},
 'axial':{
    'full_unperturbed_energy_and_pressure':(0,1,0,0,'D0',False),
    'same_absolute_pressure_memory':(0,1,0,0,'D0',True),
    'linear_axial_moment_history':(0,0,1,0,'D1',False),
    'selected_backward_energy_loss':(0,1,2,0,'D0',False)}}


def current_sources():
    records,hashes,family=gap_sources()
    for stem in ('pulse_gap_similarity_C4','pulse_gap_similarity_C4_check',
        'pulse_gap_physical_C2','pulse_gap_physical_C2_check'):
        name=PREFIX+stem+'.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
        if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=family:
            raise ValueError('Gap cone source family differs: '+stem)
        if 'all_passed' in record and not record['all_passed']:
            raise ValueError('Unaccepted gap cone source: '+stem)
        for path,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('Gap cone source changed: '+path)
        hashes.update(record['input_hashes']);hashes[name]=hashlib.sha256(raw).hexdigest()
        records[stem]=record
    if records['pulse_gap_similarity_C4']['domain']!=DOMAIN:
        raise ValueError('Whole original gap required')
    for flag in ('actual_pulse_gap_physical_decomposition_constructed',
        'gap_end_completed_physical_interface_verified',
        'actual_whole_gap_nonzero_radial_history_and_remainder_preserved',
        'actual_gap_axial_velocity_and_remainder_structural_zero'):
        if not records['pulse_gap_physical_C2_check'][flag]:
            raise ValueError('Gap physical source prerequisite missing: '+flag)
    return records,hashes,family


def relative_log_envelopes(c,mu,finite,lp,lu,lrp):
    """Exact endpoint substitution after grouping; never subtract log boxes."""
    return {
      'theta_signed_original_memory':-11*(1-mu)/mu,
      'theta_meridional_Mz_history':lp+lu+finite-c.mpf('5.5')/mu-13,
      'theta_mixed_angular_axial_history':lp+lu+finite-c.mpf('5.5')/mu-15,
      'theta_radial_shear':-lrp-11/mu,
      'axial_full_unperturbed_energy_and_pressure':lp+lu-c.mpf('5.5')/mu-11,
      'axial_same_absolute_pressure_memory':lp+lu-c.mpf('6.5')/mu-15-4*mu,
      'axial_linear_axial_moment_history':finite-2,
      'axial_selected_backward_energy_loss':lp+lu+2*finite-c.mpf('7.5')/mu-11}


def gap_cone_identities(records):
    asts=SourceAST();proofs={}
    mu,delta,z,d,lp,lu,lrp,finite=s.symbols('mu delta Z distance logP logU logRp finite',real=True)
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),ln=s.log)
    native=SimpleNamespace(mu=mu,finite=finite,logP=lp,logU=lu,logRp=lrp)
    env=dict(c=c,self=native,d=d,p=1+2*mu,r=1-mu)
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:
            raise ArithmeticError('Gap cone identity failed: '+name)
        proofs[name]=True
    logB=asts.evaluate(asts.expression('pulse_gap_similarity_C4','packet','logBparts'),env)
    logR=asts.evaluate(asts.expression('pulse_gap_similarity_C4','packet','logR'),env)
    logH=asts.evaluate(asts.expression('pulse_gap_similarity_C4','packet','logH'),env)
    env['logD0']=asts.evaluate(asts.expression('pulse_gap_similarity_C4','packet','logD0'),env)
    logD=asts.evaluate(asts.expression('pulse_gap_similarity_C4','packet','logD'),env)
    logQ=asts.evaluate(asts.expression('pulse_gap_similarity_C4','packet','logQ'),env)
    grouped={
      'theta_signed_original_memory':-(1-mu)*(13-d)/mu,
      'theta_meridional_Mz_history':lp+lu+finite+(-s.Rational(15,2)+d)/mu-13,
      'theta_mixed_angular_axial_history':lp+lu+finite+(-s.Rational(15,2)+d)/mu-13-d,
      'theta_radial_shear':-lrp-(13-d)/mu,
      'axial_full_unperturbed_energy_and_pressure':lp+lu+(-13+d)/(2*mu)-13+d,
      'axial_same_absolute_pressure_memory':lp+lu+(-13-d)/(2*mu)-13-d,
      'axial_linear_axial_moment_history':finite+(d/2-1)/mu-d,
      'axial_selected_backward_energy_loss':lp+lu+2*finite+(-17+d)/(2*mu)-13+d}
    envelopes=asts.replay('pulse_gap_cone','relative_log_envelopes',{})(c,mu,finite,lp,lu,lrp)
    rates={}
    for label,rows in FACTORS.items():
        for name,(rp,bp,dp,hp,which,pressure) in rows.items():
            key=label+'_'+name
            source=rp*logR+bp*sum(logB.values())+dp*logD[int(which[-1])]+hp*logH+(logQ if pressure else 0)
            zero('actual_'+key+'_grouped_source_log',source-grouped[key])
            rate=s.diff(grouped[key],d);rates[key]=rate
            endpoint=4*mu if pressure else s.Integer(2)
            zero('actual_'+key+'_monotone_endpoint_log_envelope',grouped[key].subs(d,endpoint)-envelopes[key])
            original=tuple(s.Rational(str(v))+w for v,w in zip((rp,bp,dp,hp),(s.Rational(1,2),1,0,0)))
            R,B,D,H=s.symbols('R B D H',positive=True)
            zero('actual_'+key+'_relative_positive_factor',
                R**original[0]*B**original[1]*D**original[2]*H**original[3]/s.sqrt(2)
                /(s.sqrt(R/2)*B)-R**rp*B**bp*D**dp*H**hp)
    asts.expression('pulse_gap_similarity_C4','gap_shapes','Bh',wanted='[zero for _ in range(5)]')
    physical=records['pulse_gap_physical_C2']['source_and_physical_join_binding']['identities']
    for j in range(5):
        flag='actual_gap_axial_velocity_structural_zero_row'+str(j)
        if not physical[flag]:raise ValueError('Actual gap zero axial source missing')
        proofs['directly_consumed_'+flag]=True
    C=1/(1+z*z);r=1-mu;b0=(1-delta)/2;L=1-delta*z*z;g=(mu-delta/2)/r
    equilibrium=asts.evaluate(asts.expression('pulse_end_stress_C3','pulse_coefficients','equilibrium'),
        dict(C=C,mu=mu,delta=delta,r=r,z=z,b=b0,L=L,axial_derivative=lambda value:s.diff(value,z)))
    zero('actual_whole_gap_positive_equilibrium_rewrite',
        equilibrium-(C*g+2*b0*z*z/(r*(1+z*z)**2))/L)
    path=HERE/'lei_ren_part1_paper_mp_stress.py';asts.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    original=next(n for n in ast.walk(ast.parse(path.read_text(encoding='utf8')))
        if isinstance(n,ast.FunctionDef) and n.name=='evaluate_mp_stress')
    R,B=s.symbols('R B',positive=True);F=B*C/s.sqrt(2*R);A=(2+2*mu)*F
    shearenv=dict(a=B*C,ay=-(s.Rational(1,2)+mu)*B*C,by=s.Integer(0),root=s.sqrt(2*R),R=R)
    shear={}
    for target in ('Stheta','Sz'):
        node=next(n.value for n in ast.walk(original) if isinstance(n,ast.Assign)
            and any(ast.unparse(t)==target for t in n.targets))
        if isinstance(node,ast.IfExp):node=node.body
        shear[target]=asts.evaluate(node,shearenv)
    zero('actual_negative_theta_shear',shear['Stheta']+A)
    zero('actual_zero_axial_shear_from_gap_source',shear['Sz'])
    zero('actual_gap_kappa_minus2_exact_without_rounded_subtraction',
        -(shear['Stheta']**2+shear['Sz']**2)/(F*shear['Stheta'])-2-2*mu)
    theta,axial,Q,aa=s.symbols('theta axial Q A',real=True)
    dot=-Q*aa*theta;cross=-Q*aa*axial
    zero('actual_gap_negative_dot_normalization',-dot-Q*aa*theta)
    zero('actual_gap_cross_direction_normalization',-cross-Q*aa*axial)
    zero('actual_gap_directional_cone_normalization',
        2*dot**2-2*mu*cross**2-(Q*aa)**2*(2*theta**2-2*mu*axial**2))
    nu,fac,ff,st,sz=s.symbols('nu lambda_factor F Stheta Sz',positive=True)
    zero('actual_fixed_positive_viscosity_physical_kappa',
        -((nu*fac*st)**2+(nu*fac*sz)**2)/(nu*(fac*ff)*(nu*fac*st))+(st**2+sz**2)/(ff*st))
    zero('actual_physical_directional_margin_positive_common_factor',
        2*((nu*fac)**2*dot)**2-2*mu*((nu*fac)**2*cross)**2-(nu*fac)**4*(2*dot**2-2*mu*cross**2))
    tree=ast.parse((HERE/(PREFIX+'global_physical_assembly.py')).read_text(encoding='utf8'))
    radiusenv={}
    for node in tree.body:
        if isinstance(node,ast.Assign):
            for target in node.targets:
                if isinstance(target,ast.Name) and target.id in ('MICRO','POST','PULSE'):
                    radiusenv[target.id]=ast.literal_eval(node.value)
    radius=asts.replay('global_physical_assembly','radius',radiusenv)
    wait,Ts,length=s.symbols('wait Ts Lrel',real=True)
    assembly=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRp=lrp,params=SimpleNamespace(mu=mu))
    heat=SimpleNamespace(steep=SimpleNamespace(wait=wait,Ts=Ts,outer=SimpleNamespace(Lrel=length)))
    left=radius(assembly,'pulse_gap',11,{},None)[0]
    right=radius(assembly,'pulse_gap',13-4*mu,{},None)[0]
    end=radius(assembly,'pulse_end',-4,{},None)[0]
    tail=radius(assembly,'heat_collar',0,{},heat)[0]
    zero('actual_production_gap_left_radius',left-lrp-11/mu)
    zero('actual_production_gap_right_equals_end_left_radius',right-end)
    zero('actual_production_gap_left_tail_offset',left-tail+2/mu+wait+Ts+102+length)
    for stem,flag in (('pulse_gap_similarity_C4_check','full_nonzero_raw_moment_radial_and_stress_histories_preserved'),
        ('pulse_gap_physical_C2_check','gap_end_completed_physical_interface_verified'),
        ('pulse_end_cone_check','pulse_end_flatten_power_angular_entry_exit_waiting_collar_two_vector_cone_certified'),
        ('fifth_axial_jets_check','actual_selected_ap_c1_c2_C5_available')):
        if not records[stem][flag]:raise ValueError('Cone composition source missing: '+flag)
        proofs['directly_consumed_'+flag]=True
    return dict(identities=proofs,input_hashes=asts.hashes,actual_source_AST_bindings=asts.bindings,
        exact_grouped_relative_log_recipes={key:str(value) for key,value in grouped.items()},
        source_log_derivative_rates={key:str(value) for key,value in rates.items()},
        monotonicity_scope='0<mu<1/4; actual coupled domain4mu<=d<=2; BQ max at4mu, others at2',
        equilibrium_lower_formula='theta_equilibrium>=g/2, g=(mu-delta/2)/(1-mu)>0',
        actual_kappa='kappa-2=2mu from actual zero Uz_y; radial history remains',
        physical_kappa='-|Sphys|^2/(nu*Fphys*Sphys_theta), identical fixed nu>0',
        regional_two_vector_cone_distinct_from_completed_tensor_cone=True)


def absolute_bound(c,value):return c.mpf(max(abs(v) for v in endpoints(value)))


@source_precision
def whole_gap_bounds(records):
    c=MPIntervalContext();c.dps=240
    selected=records['pulse_mixed_C4'];mu=read_interval(c,selected['selected_mu'])
    delta=read_interval(c,selected['selected_delta']);r=1-mu;b0=(1-delta)/2;g=(mu-delta/2)/r
    parameters=dict(mu=mu,delta=delta,quarter_minus_mu=c.mpf('.25')-mu,
        rate=r,one_minus_delta=1-delta,b0=b0,equilibrium_gap=g)
    def positive(name,value):
        lo,hi=endpoints(value)
        if lo<=0 or not mp.isfinite(hi):raise ArithmeticError('Whole gap positive margin failed: '+name)
    for name,value in parameters.items():positive(name,value)
    g_lower=c.mpf(endpoints(g)[0]);theta_lower=g_lower/4
    eq_upper=c.mpf(endpoints((g+2*b0/r)/(1-delta))[1]);theta_upper=eq_upper+g_lower/8
    epsilon=c.exp(-1000);axial_cap=theta_lower*epsilon
    source=records['pulse_gap_similarity_C4']['whole_original_gap']
    if not source['zero_axial_input_does_not_zero_radial_history']:
        raise ValueError('Actual gap radial history missing')
    eqlogs=source['full_meridional_stress_log_sectors']['theta']['equilibrium']['exact_source_log_parts']
    lp=read_interval(c,eqlogs['logPstar']);lu=read_interval(c,eqlogs['actual_log_inlet_U'])
    rows=records['outer_pulse_map']['pulse_rows']
    finite=-3*read_interval(c,rows['saddle_L'])+2*c.ln(read_interval(c,rows['saddle_u0']))-c.ln(6)/2-2*c.ln(mu)
    logC=read_interval(c,records['outer_initial']['selected_logCstar'])
    Tw=read_interval(c,records['fifth_axial_jets']['whole_Z']['selected']['incoming']['Z_independent_constant_definitions']['Tw'])
    lrp=c.ln(110)+10*(logC+lp)+lp+1+Tw
    envelopes=relative_log_envelopes(c,mu,finite,lp,lu,lrp);bounds={};margins={}
    sectors=source['full_meridional_stress_log_sectors']
    for label,names in FACTORS.items():
        if set(sectors[label])!=set(names)|({'equilibrium'} if label=='theta' else set()):
            raise ValueError('All original gap stress sectors required')
        for name,(rp,bp,dp,hp,which,pressure) in names.items():
            sector=sectors[label][name]
            expected=tuple(v+w for v,w in zip((rp,bp,dp,hp),(.5,1,0,0)))
            if tuple(sector['original_mode'])!=expected or sector['selected_D_recipe']!=which:
                raise ValueError('Original gap relative factor changed: '+name)
            if ('same_absolute_pressure_memory' in sector['exact_source_log_parts'])!=pressure:
                raise ValueError('Original pressure memory factor changed')
            key=label+'_'+name
            coefficient=read_interval(c,sector['full_stress_mixed3_coefficient_enclosures']['s0_Z0'])
            upper=absolute_bound(c,coefficient)
            if endpoints(upper)[1]<=0:raise ValueError('Expected nonzero gap history bound')
            actual_log=envelopes[key]+c.ln(upper)
            target=c.ln(g_lower/64) if label=='theta' else c.ln(theta_lower)-1000-c.ln(4)
            gap=target-actual_log;positive(key,gap);margins[key]=gap
            bounds[key]=dict(actual_signed_coefficient=coefficient,coefficient_absolute_upper=upper,
                relative_mode=(rp,bp,dp,hp),selected_D_recipe=which,pressure_memory_retained=pressure,
                exact_whole_gap_relative_log_envelope=envelopes[key],
                actual_relative_log_absolute_upper=actual_log,certified_log_cap=target,positive_log_gap=gap)
    algebraic=dict(theta_lower=theta_lower,theta_upper=theta_upper,
        equilibrium_lower_slack=g_lower/2-4*g_lower/64-theta_lower,
        equilibrium_upper_slack=g_lower/8-4*g_lower/64,
        negative_dot_lower=theta_lower,cross_absolute_upper=axial_cap,
        kappa_minus2=2*mu,kappa_minus2_below2=2-2*mu,
        directional_margin=2*theta_lower**2-2*mu*axial_cap**2)
    for name,value in algebraic.items():positive(name,value)
    return dict(positive_parameter_margins=parameters,positive_log_margins=margins,
        positive_algebraic_margins=algebraic,normalized_stress_sector_bounds=bounds,
        source_reduced_log_parameters=dict(logPstar=lp,actual_log_inlet_U=lu,finite=finite,logRp=lrp),
        actual_whole_gap_log_envelopes=envelopes,axial_uniform_absolute_upper=axial_cap,
        exact_axial_shear_ratio=c.mpf(0),normalized_cone='2*theta^2-2mu*axial^2>0',
        continuous_whole_original_gap_Z_domain_covered=True,
        exact_source_correlations_grouped_before_interval_enclosure=True,
        all_four_angular_and_four_axial_histories_retained=True,
        zero_axial_shear_proved_from_actual_gap_source=True,
        radial_velocity_and_five_moment_histories_not_zeroed=True,
        main_pulse_axial_shear_not_assumed_zero=True,
        phase_samples_used_as_proof=False,source_caps_used_as_defining_field_values=False,
        global_temporal_flatness_not_inferred=True)


@source_precision
def run():
    records,hashes,family=current_sources()
    proof=gap_cone_identities(records);hashes.update(proof['input_hashes'])
    bounds=whole_gap_bounds(records);hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result=dict(actual_five_defect_family_sha256=family[0],implicit_source_sha256=family[1],
        domain=DOMAIN,tail_domain=TAIL_DOMAIN,exact_source_cone_proof=proof,bounds=bounds,input_hashes=hashes,
        current_selected_C5_source_check_directly_consumed=True,
        regional_two_vector_cone_distinct_from_completed_tensor_cone=True,
        **{flag:True for flag in ADMISSIONS},**{flag:False for flag in FALSE_FLAGS})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Continuous whole original gap cone generated; main/entrance/global/recursion pending',flush=True)
    return result


class CertifiedPulseGapPhysical:
    """Add current regional admission to the same completed physical evaluator."""
    def __init__(self):
        from lei_ren_part1_paper_compliant_pulse_gap_physical_C2 import CompliantPulseGapPhysicalC2
        self.receipt=json.loads(Path(__file__).with_name(PREFIX+'pulse_gap_cone_check.json').read_bytes())
        if not self.receipt['all_passed'] or not self.receipt['pulse_gap_cone_certified']:
            raise ValueError('Accepted whole gap cone required')
        for path,digest in self.receipt['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('Certified gap physical source changed: '+path)
        self.physical=CompliantPulseGapPhysicalC2()
        if (self.physical.family,self.physical.source)!=(
            self.receipt['actual_five_defect_family_sha256'],self.receipt['implicit_source_sha256']):
            raise ValueError('Certified gap physical family differs')
    def gap(self,*args,**kwargs):
        result=self.physical.gap(*args,**kwargs)
        result.update(**{flag:True for flag in ADMISSIONS},**{flag:False for flag in FALSE_FLAGS})
        return result


if __name__=='__main__':run()
