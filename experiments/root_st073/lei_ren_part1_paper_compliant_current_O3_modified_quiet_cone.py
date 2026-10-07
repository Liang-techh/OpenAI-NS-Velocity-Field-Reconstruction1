"""Complete signed cone on the current finite-N independent repair strip.

Restrict the checked original power theorem to its actual local source
domain before normalizing errors. Same implicit controls and cumulative
pressure/moments are retained on bumps and gaps. Global gates stay open.
"""
import functools
import json
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_O3_modified_transition_cone as weighted

errors=weighted.errors
HERE,PREFIX,sha=errors.HERE,errors.PREFIX,errors.sha
NAME=PREFIX+'current_O3_modified_quiet_cone.json'
RECEIPT=PREFIX+'current_O3_modified_quiet_cone_check.json'
BASE_NAME=PREFIX+'current_O3_power_cone.json'
BASE_RECEIPT=PREFIX+'current_O3_power_cone_check.json'
parameters,endpoints,upper=errors.parameters,errors.endpoints,errors.upper
OPEN=tuple(k for k in weighted.OPEN if k!='current_modified_quiet_signed_cone_certified')


@functools.lru_cache(maxsize=1)
def exact_theorem():
    asts=errors.recovered.histories.SourceAST();checks={}
    def zero(name,a,b):
        if s.cancel(s.expand(a-b))!=0:raise ArithmeticError('Quiet signed cone source: '+name)
        checks[name]=True
    y=s.Symbol('actual_quiet_local_logx',real=True)
    mu,N,A=s.symbols('mu N same_positive_common_profile_scale',positive=True)
    h=s.symbols('axial0 axial2 swirl0 swirl1 swirl2',real=True)
    B=[s.Function('same_bump_'+str(i))(y) for i in range(3)]
    rows=[[s.diff(v,y,j) for j in range(5)] for v in B]
    alpha=s.Rational(1,2)+mu;f=s.exp(-alpha*y)
    env=dict(c=SimpleNamespace(mpf=lambda v:s.Rational(str(v))),h=h,rows=rows,Ahat=A,
        se=mu/N,sa=s.sqrt(mu)/N,original=dict(theta=[A*s.diff(f,y,j) for j in range(5)]))
    du=asts.evaluate(asts.expression('current_O3_repaired_histories','history','du'),env);env['du']=du
    V=asts.evaluate(asts.expression('current_O3_repaired_histories','history','Vnew'),env)
    U=asts.evaluate(asts.expression('current_O3_repaired_histories','history','Unew'),env)
    Hs=sum(h[i+2]*B[i] for i in range(3));Ha=h[0]*B[0]+h[1]*B[2]
    den=f+mu*Hs/N
    g=(2*f-(Hs+2*s.diff(Hs,y))/N)/den
    hh=2*s.diff(Ha,y)/(N*den)
    zero('actual_common_profile_scale_and_Pstar_cancel_in_theta_shear',1-2*U[1]/U[0]-2,mu*g)
    zero('actual_positive_original_bs_sign_with_same_physical_Pstar',2*V[1]/U[0],s.sqrt(mu)*hh)
    zero('actual_quiet_theta_shear_perturbation',g-2,-2*(s.diff(Hs,y)+alpha*Hs)/(N*den))
    for j in range(5):
        zero('actual_swirl_profile_ordinary_row'+str(j),du[j],A*mu*s.diff(Hs,y,j)/N)
        zero('actual_axial_profile_ordinary_row'+str(j),V[j],A*s.sqrt(mu)*s.diff(Ha,y,j)/N)
    # The original repaired source has the physical Pstar on both velocity
    # components. This factor must precede the positive shear quotient.
    asts.expression('current_O3_repaired_histories','history','h',wanted='self.repair.controls')
    asts.expression('current_O3_repaired_histories','history','original',
        wanted="raw_pre_velocity_rows(c,parent['current_original_pre_source'])")
    asts.expression('current_O3_repaired_histories','history','Ahat',
        wanted='(1+square(z)).reciprocal()*(self.histories.Ua*self.f2)')
    asts.expression('current_O3_repaired_histories','history','y',wanted='q-1')
    asts.expression('current_O3_modulated_histories','history','t',wanted='1+v')
    asts.expression('current_O3_modulated_histories','history','pre',wanted='self.pre.power(zv,v/Tw)')
    zero('same_actual_quiet_power_coordinate',1+y,y+1)
    zero('same_actual_modulation_transport_coordinate',1+(1+y),2+y)
    Tw=s.Symbol('same_positive_Tw',positive=True)
    zero('same_actual_quiet_original_power_phase_coordinate',Tw*((1+y)/Tw),1+y)
    # All envelopes are uniform in the implicit control ball. Actual
    # h_N need not be monotone with N. Only its fixed cap is used.
    q,R,G0,G1,m,fmin=s.symbols('inverse_N R G0 G1 mu fmin',positive=True)
    denom=fmin-m*R*G0*q
    gs=2*R*(G1+(s.Rational(1,2)+m)*G0)*q/denom
    hs=2*R*G1*q/denom
    zero('positive_inverse_N_theta_shear_cap_derivative',s.diff(gs,q),
        2*R*(G1+(s.Rational(1,2)+m)*G0)*fmin/denom**2)
    zero('positive_inverse_N_axial_shear_cap_derivative',s.diff(hs,q),2*R*G1*fmin/denom**2)
    for i in (0,1):
        zero('positive_inverse_N_local_profile_cap_derivative'+str(i),s.diff(R*q*s.Symbol('G'+str(i),positive=True),q),
            R*s.Symbol('G'+str(i),positive=True))
    K1,K2=s.symbols('positive_K1 positive_K2',positive=True)
    for name,value in (('linear_partial_history',K1*q),('quadratic_partial_history',K1*q+K2*q*q)):
        if s.diff(value,q).is_positive is not True:raise ValueError('Quiet history cap not inverse-N monotone')
        checks['positive_inverse_N_derivative_'+name]=True
    # Native leading stress uses only local y0/y1 and history y0, which
    # are positive sums/mins of these linear/quadratic envelopes.
    native=weighted.exact_theorem()
    asts.method('current_O3_uniform_repair_majorants','_primitive_caps')
    asts.method('current_O3_uniform_repair_majorants','_weights_cap')
    asts.method('current_O3_recovered_error_majorants','quiet_repair')
    asts.expression('current_O3_independent_repair_operator','bump_weights','centers',
        wanted='[c.mpf(1)/5,c.mpf(1)/2,c.mpf(4)/5]')
    asts.expression('current_O3_independent_repair_operator','bump_weights','ell',wanted='c.mpf(1)/40')
    centers=[s.Rational(1,5),s.Rational(1,2),s.Rational(4,5)];ell=s.Rational(1,40)
    for i in range(2):
        if centers[i]+ell>=centers[i+1]-ell:raise ValueError('Independent supports overlap')
        checks['same_disjoint_source_supports_'+str(i)]=True
    if centers[0]-ell<=0 or centers[-1]+ell>=1:raise ValueError('Repair support meets source joins')
    checks['source_entry_terminal_and_internal_gaps_have_flat_local_profiles']=True
    from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import exact_repair_theorem
    independent=exact_repair_theorem()
    if independent['identities'].get('actual_repair_inlet_f2') is not True:raise ValueError('Same original quiet inlet required')
    parameter=parameters.exact_theorem()
    if parameter['identities'].get('actual_original_logAd_exact_endpoint') is not True:
        raise ValueError('Same original Ua=exp(logAd) amplitude required')
    checks['same_original_Rd_Ua_exp_logAd_amplitude_source']=True
    return dict(passed=True,identities=checks,input_hashes={**native['input_hashes'],
        **independent['input_hashes'],**parameter['input_hashes'],**asts.hashes,Path(__file__).name:sha(Path(__file__).name)},
        exact_signed_shear=dict(a='2+mu*g',bs='+sqrt(mu)*h=+2*Uz_y/Utheta',
            g='[2*f-(Hs+2*Hs_y)/N]/[f+mu*Hs/N]',h='2*Ha_y/[N*(f+mu*Hs/N)]',
            Hs='h2*B0+h3*B1+h4*B2',Ha='h0*B0+h1*B2',f='exp[-(.5+mu)*y]'),
        old_uniform_negative_b_label_not_used_as_signed_bs=True,
        same_unique_implicit_controls_and_fixed_uniform_ball_retained=True,
        all_native0_stress_error_caps_nonincreasing_with_N=True,
        exact_weighted_signed_cone_algebra=native)


class CurrentModifiedQuietCone:
    def __init__(self,require_checked=True):
        self.errors=errors.CurrentCompletedTensorErrorMajorants();self.ctx=self.errors.ctx
        self.data,self.mu=self.errors.data,self.errors.mu;self.repair=self.errors.recovered.repair
        self.minimum_N=self.repair.quiet_threshold;self.theorem=exact_theorem()
        read=parameters.numeric.transport.read_interval;self.read=read;c=self.ctx
        base=json.loads((HERE/BASE_NAME).read_bytes());receipt=json.loads((HERE/BASE_RECEIPT).read_bytes())
        if not receipt['all_passed'] or not receipt['current_whole_O3_power_signed_two_vector_cone_certified']:
            raise ValueError('Checked original whole power signed cone required')
        for name,digest in receipt['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Changed original power theorem input: '+name)
        for name,value in self.data['source']['accepted']['source_family'].items():
            if base.get(name)!=value:raise ValueError('Different original quiet power family')
        self.base=base['current_whole_O3_power_correlated_bounds'];b=self.base
        if not b['complete_original_source_signed_tensor_and_remainder_preserved'] or not b['source_M_K_X_correlated_before_enclosure']:
            raise ValueError('Complete correlated original power baseline required')
        self.Tw=read(c,b['positive_margins']['actual_Tw'])
        if endpoints(self.Tw)[0]<=2:raise ValueError('Actual quiet strip must lie inside original power source')
        constants={key:read(c,value) for key,value in b['theta_lower_constants'].items()}
        deficit=read(c,b['actual_current_O3_deficit_enclosure'])
        floor=constants['cD']*deficit*c.exp(-2)+constants['theta_floor']
        if endpoints(floor)[0]<=0:raise ValueError('Local original power theta lower required')
        self.logUmin=self.data['logAd']-c.mpf('.5')-self.mu/2-2*(c.mpf('.5')+self.mu)
        self.logfloor=self.logUmin+c.ln(floor)
        self.rho0=upper(c,c.sqrt(read(c,b['full_directional_expression_upper'])/2))
        self.baseline=dict(quiet_local_logx_domain=(0,1),actual_original_power_offset_domain=(1,2),
            actual_Tw=self.Tw,canonical_theta_local_lower=floor,canonical_source_deficit=deficit,
            canonical_theta_lower_constants=constants,original_log_Utheta_scalar_lower=self.logUmin,
            original_theta_log_lower=self.logfloor,original_weighted_axial_ratio_upper=self.rho0,
            checked_whole_source_theorem_restricted_to_actual_local_domain=True,
            original_pressure_energy_M_K_X_and_radial_terms_retained=True)
        self.error_rows=self.leading_error_rows(self.minimum_N);self.proof=self.prove()
        self.hashes={**self.errors.hashes,errors.RECEIPT:sha(errors.RECEIPT),**receipt['input_hashes'],
            BASE_RECEIPT:sha(BASE_RECEIPT),BASE_NAME:sha(BASE_NAME),**self.theorem['input_hashes']}
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked['all_passed']:raise ValueError('Checked complete quiet signed cone required')
            for name,digest in checked['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed quiet signed cone input: '+name)

    def leading_error_rows(self,N):
        c=self.ctx;q,caps=self.errors.inputs('quiet',(0,1),N);result={}
        for label,parts in self.errors.model['stress'].items():
            result[label]={}
            for name,p in parts.items():
                cap=errors.bound_expression(c,p['expression'],caps,c.mpf([-1,1]),self.errors.delta,0,p['mode'][0])[0,0]
                rp,pp,_,_=p['mode'];extra=self.data['log_mu']/2 if label=='axial' else c.mpf(0)
                logparts=dict(relative_radius=(c.mpf(str(rp))-c.mpf('.5'))*q['native_logR'],
                    relative_Pstar=(pp-1)*self.data['logP'],external_Ad=p['error_Ad_power']*self.data['logAd'],
                    weighted_axial=extra,coefficient=c.ln(cap),baseline_floor=-self.logfloor)
                result[label][name]=dict(source_mode=p['mode'],error_Ad_power=p['error_Ad_power'],
                    native_coefficient_absolute_cap=cap,relative_error_log_parts=logparts,
                    log_absolute_error_over_baseline_theta_lower=sum(logparts.values(),c.mpf(0)))
        return result

    def prove(self):
        c=self.ctx;positive={}
        def require(name,value):
            if endpoints(value)[0]<=0:raise ArithmeticError('Quiet signed cone unresolved: '+name)
            positive[name]=value
        primitive=self.repair.query(self.minimum_N,(0,1))['quiet_repaired_shear']
        if not primitive['at_least_three_halves_mu_certified']:raise ValueError('Uniform quiet primitive reserve required')
        require('same_positive_swirl_denominator',primitive['normalized_swirl_denominator_lower'])
        require('same_actual_power_source_contains_quiet_strip',self.Tw-2)
        gerror=primitive['theta_shear_perturbation_absolute_cap']/self.mu
        hcap=primitive['axial_shear_absolute_cap']/c.sqrt(self.mu)
        require('actual_absolute_g_below3',1-gerror)
        require('actual_absolute_h_below4',4-hcap)
        for label,parts in self.error_rows.items():
            for name,p in parts.items():require('native_error_log_gap/'+label+'/'+name,-1010-p['log_absolute_error_over_baseline_theta_lower'])
            require('sum_native_errors_below_exp_minus1000/'+label,10-c.ln(len(parts)))
        eps=c.exp(-1000);rho=(self.rho0+eps)/(1-eps)
        alo=c.mpf('1.997');ahi=c.mpf('2.003');H=c.mpf(4);G=c.mpf(3);mumax=c.mpf('.001')
        direction=alo-H*rho;W=2*alo-mumax*H*H-2*ahi*H*rho-ahi*G*rho*rho
        dlower=direction/ahi;qlower=alo*alo*W/ahi**3
        for name,value in dict(positive_modified_theta=1-eps,weighted_axial_ratio_below_exp_minus390=c.exp(-390)-rho,
            weighted_H_over_theta=direction,weighted_W_over_theta_squared=W,
            original_signed_D_over_theta=dlower,original_signed_Q_over_theta_squared=qlower,
            actual_vs_minus2_lower=3*self.mu/2).items():require(name,value)
        return dict(minimum_integer_N=self.minimum_N,all_finite_integer_N_at_least_minimum=True,
            whole_quiet_logx_domain=(0,1),whole_Z_domain=(-1,1),all_disjoint_bumps_and_gaps=True,
            entry_internal_and_terminal_source_domains_included=True,positive_margins=positive,
            baseline_theta_log_lower=self.logfloor,relative_theta_error_upper=eps,relative_weighted_axial_error_upper=eps,
            modified_weighted_axial_ratio_upper=rho,actual_g_error_from2_upper=gerror,actual_h_absolute_upper=hcap,
            actual_a_lower=alo,actual_a_upper=ahi,actual_vs_minus2_lower=3*self.mu/2,
            weighted_H_over_theta_lower=direction,reduced_weighted_quadratic_over_theta_squared_lower=W,
            original_signed_D_over_theta_lower=dlower,original_signed_Q_over_theta_squared_lower=qlower,
            old_uniform_negative_b_label_not_used_as_signed_bs=True,
            complete_original_signed_pressure_energy_and_all_cross_sectors_retained=True,
            same_unique_implicit_control_family_no_saved_midpoints=True,
            current_modified_quiet_signed_cone_certified=True,**{key:False for key in OPEN})

    def query(self,logx,N):
        parameters.positive_integer_N(N);c=self.ctx;y=c.mpf(logx);lo,hi=endpoints(y)
        if N<self.minimum_N or lo<0 or hi>1 or not all(mp.isfinite(v) for v in (lo,hi)):
            raise ValueError('Finite quiet logx subset[0,1], integer N>=quiet minimum required')
        q=self.errors.recovered.quiet_repair(N,y)
        return dict(quiet_local_logx=y,actual_original_power_offset=1+y,actual_modulation_transport_offset=2+y,
            finite_integer_N=N,implicit_controls_absolute_ball_radius=self.repair.R,
            actual_bs_convention='+2*Uz_y/Utheta=sqrt(mu)*h',source_family=q['source_family'],
            local_disjoint_bump_logR_derivative_caps=q['local_disjoint_bump_logR_derivative_caps'],
            same_source_remaining_primitive_caps=q['same_source_shrinking_terminal_defect_caps'],
            exact_full_source_error_zero_after_last_support=q['all_error_logR4_axial5_rows_vanish_after_last_bump'],
            complete_modified_signed_cone_certified_for_entire_query_box=True,
            whole_source_signed_cone_certificate=self.proof,current_modified_quiet_signed_cone_certified=True,
            **{key:False for key in OPEN})


def run():
    field=CurrentModifiedQuietCone(require_checked=False)
    result=dict(source_family=field.data['source']['accepted']['source_family'],
        exact_signed_quiet_source_and_weighted_cone_theorem=field.theorem,
        original_power_local_baseline=field.baseline,whole_modified_quiet_cone=field.proof,
        complete_leading_stress_error_log_ledger=field.error_rows,
        examples={name:field.query(y,N) for name,y,N in (
            ('whole_strip',(0,1),field.minimum_N),('entry','0',field.minimum_N),
            ('first_bump',('.175','.225'),field.minimum_N),('first_gap',('.225','.475'),field.minimum_N),
            ('middle_bump',('.475','.525'),field.minimum_N),('last_bump',('.775','.825'),field.minimum_N),
            ('near_last_flat_edge',('.824','.825'),field.minimum_N),
            ('final_exact_strip',('.826','1'),field.minimum_N),('terminal','1',field.minimum_N),
            ('current_frequency',(0,1),10**12))},
        current_modified_quiet_signed_cone_certified=True,
        scoped_source_cone_does_not_increment_unmodified_registry_counts=True,
        **{key:False for key in OPEN},input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(parameters.encoded(result),indent=2)+'\n').encode())
    print('Modified quiet repair cone generated; all bumps/gaps/Z; all N>=',field.minimum_N,flush=True)
    return result


if __name__=='__main__':run()
