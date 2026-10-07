"""Signed modified O2 taper cone, with the degenerate flat edge explicit.

The original negative-offset buffer baseline is derived separately from
O3. Weighted full stress/error algebra is shared. This is not a strict
closed O2 or quiet/global cone certificate.
"""
import functools
import importlib
import json
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_O3_modified_transition_cone as o3
from lei_ren_part1_paper_compliant_outer_initial import stable_sigma

HERE,PREFIX,sha=o3.HERE,o3.PREFIX,o3.sha
NAME=PREFIX+'current_O2_modified_taper_cone.json'
RECEIPT=PREFIX+'current_O2_modified_taper_cone_check.json'
MINIMUM_N=o3.MINIMUM_N
OPEN=o3.OPEN
parameters,endpoints,upper=o3.parameters,o3.endpoints,o3.upper


@functools.lru_cache(maxsize=1)
def exact_theorem():
    asts=o3.errors.recovered.histories.SourceAST();checks={}
    def zero(name,a,b):
        if s.cancel(s.expand(s.expand_power_exp(a-b)))!=0:raise ArithmeticError('O2 taper identity: '+name)
        checks[name]=True
    t,z,delta=s.symbols('negative_offset Z delta',real=True)
    Ua,Ma,da=s.symbols('same_Ua same_Ma same_da',positive=True)
    EZa,EQa=s.symbols('same_EZa same_EQa',real=True)
    C=1/(1+z*z);L=1-delta*z*z;A=(1-delta/2)*C+(1-delta)*z*z*C*C
    B=((2*delta*z*z-1)*C+2*(1-z*z)*z*z*C*C)/L
    U=Ua*s.exp(-t/2);D=da*s.exp(-t);M=Ma*s.exp(-t)
    H=U*(1-D);K=U*(M-4*D);AZ=EZa/Ua**2;AQ=EQa/Ua**2-t/2
    CE=z*z*AZ+C*C*AQ
    zero('same_original_pure_buffer_U_ODE',s.diff(U,t),-U/2)
    zero('same_original_pure_buffer_m_ODE',s.diff(M,t),-M)
    zero('same_original_pure_buffer_h_ODE',s.diff(H,t),U-s.Rational(3,2)*H)
    zero('same_original_pure_buffer_k_ODE',s.diff(K,t),-s.Rational(3,2)*K)
    zero('same_original_pure_buffer_e_ODE',s.diff(U**2*CE,t),-U**2*C*C/2-U**2*CE)
    zero('same_actual_buffer_U_squared_over_deficit',U**2/D,Ua**2/da)
    H0=s.Symbol('same_future_H0',real=True);Hfuture=1+(H0-1)*s.exp(t)
    Pmemory=s.Function('same_absolute_pressure_memory')(z)
    Pabs=-U**2*C*C*Hfuture/2+Pmemory
    zero('same_backward_pure_buffer_future_H_ODE',s.diff(Hfuture,t),Hfuture-1)
    zero('same_absolute_pressure_FTC',s.diff(Pabs,t),U**2*C*C/2)
    zero('same_future_H_convex_combination',Hfuture,(1-s.exp(t))+s.exp(t)*H0)
    shifted=importlib.import_module('lei_ren_part1_paper_compliant_collar_stress_C3').shifted_rows
    product=importlib.import_module('lei_ren_part1_paper_compliant_collar_Gamma_C4').product_rows
    op=asts.replay('current_pre_pulse_stress_operator','raw_pre_stress_rows',dict(
        axial_derivative=lambda value:s.diff(value,z),product_rows=product,shifted_rows=shifted))
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)),exp=s.exp)
    rows=lambda expr:[s.diff(expr,t,j) for j in range(5)]
    histories=dict(m=rows(M*z),h=rows(H*C),k=rows(K*z*C),e=rows(U**2*CE))
    raw=op(c,delta,z,rows(U*C),[s.Integer(0)]*5,histories,rows(Pabs))
    theta=(A*(1-D)-C)/L+C*M+B*(M-4*D)
    zero('same_complete_original_O2_inertial_theta',
        sum(p['shape'][0] for name,p in raw['theta'].items() if name!='variable_radial_shear')/U,theta)
    zero('same_O2_correlated_positive_theta_split',theta,(A-C)/L+(-A/L-4*B)*D+(C+B)*M)
    zero('same_O2_radial_shear',raw['theta']['variable_radial_shear']['shape'][0]/U,-2*C)
    IE=(2*delta*z*CE-(1-z*z)*s.diff(CE,z))/L
    IP=(2*(1+delta)*z*(-C*C*Hfuture/2)-(1-z*z)*s.diff(-C*C*Hfuture/2,z))/L
    mem=(2*(1+delta)*z*Pmemory-(1-z*z)*s.diff(Pmemory,z))/L
    zero('same_full_O2_energy_and_absolute_pressure',
        (raw['axial']['retained_full_energy']['shape'][0]+raw['axial']['actual_absolute_pressure']['shape'][0])/U**2,
        IE+IP+mem/U**2)
    zero('same_full_O2_energy_future_pressure_odd_Z_factor',IE+IP,
        z*((2*delta*z*z-2*(1-z*z))*AZ+(2*delta*C*C+4*(1-z*z)*C**3)*AQ
            -Hfuture*((1+delta)*C*C+2*(1-z*z)*C**3))/L)
    for name in ('local_axial_transport','nonlinear_meridional_transport','retained_linear_axial_moment','axial_radial_shear'):
        for j,value in enumerate(raw['axial'][name]['shape']):zero('same_original_O2_'+name+'_zero'+str(j),value,0)
    # Replay the actual source assignments on the post-turnoff buffer.
    # The defining kernel is int_1^y exp(s-y) B(s)^j ds; B(y)=0 here,
    # hence each nonzero accumulated kernel satisfies K_y=-K.
    T,T0=s.symbols('actual_buffer_logR actual_buffer_inlet_logR',real=True)
    Uc,mc,hc,kc,ec,qc,pc,bm,bq,invP2=s.symbols('Uc mc hc kc ec qc pc bm bq invP2',real=True)
    env=dict(c=c,z=z,self=SimpleNamespace(invP2=invP2),t=T,d=s.exp(-T),root=s.exp(-T/2),d3=s.exp(-3*T/2),
        u1=Uc*C,old=dict(m=mc*z,h=hc*C,k=kc*z*C,e=ec*z*z+qc*C*C,p=pc),
        K=dict(B_mass=bm*s.exp(-T),B_squared_mass=bq*s.exp(-T)),square=lambda v:v*v)
    source_u=asts.evaluate(asts.expression('pre_pulse_mixed_C4','axial','u',wanted='u1*root'),env)
    source_hist=asts.evaluate(asts.expression('pre_pulse_mixed_C4','axial','hist'),env)
    for label,target in dict(m=-source_hist['m'],h=source_u-s.Rational(3,2)*source_hist['h'],
        k=-s.Rational(3,2)*source_hist['k'],e=-source_u**2/2-source_hist['e'],p=source_u**2/2).items():
        zero('same_replayed_actual_buffer_'+label+'_ODE',s.diff(source_hist[label],T),target)
    source_at={name:value.subs(T,T0) for name,value in source_hist.items()};uat=source_u.subs(T,T0)
    relative=dict(m=source_at['m']*s.exp(-t),h=source_at['h']*s.exp(-3*t/2)+uat*(s.exp(-t/2)-s.exp(-3*t/2)),
        k=source_at['k']*s.exp(-3*t/2),e=source_at['e']*s.exp(-t)-uat**2*t*s.exp(-t)/2,
        p=source_at['p']+uat**2*(1-s.exp(-t))/2)
    for label,value in relative.items():zero('same_replayed_actual_buffer_backward_'+label,source_hist[label].subs(T,T0+t),value)
    zero('same_replayed_actual_buffer_backward_U',source_u.subs(T,T0+t),uat*s.exp(-t/2))
    asts.expression('pre_pulse_mixed_C4','axial','B',wanted='[c.mpf(0)]*5')
    asts.expression('pre_pulse_mixed_C4','axial','parent',wanted='self.inlet(Z)')
    asts.expression('pre_pulse_mixed_C4','axial','y',wanted='c.exp(md)+selector')
    asts.expression('pre_pulse_mixed_C4','axial','t',wanted='y-1')
    asts.expression('pre_pulse_mixed_C4','axial','K',wanted='turnoff_kernels(c,y,self.params.Md,self.cells,self.window)')
    asts.expression('pre_pulse_mixed_C4','slope_mu','parent',wanted='self.axial(Z,buffer_offset=11)')
    asts.expression('current_O3_transition_background_tensor','chart','P',
        wanted="[copy_jet(c,raw['p'][0])+datum]+[copy_jet(c,row) for row in raw['p'][1:]]")
    asts.method('outer_initial','turnoff_kernels')
    # Pabs includes the original absolute datum. Matching its exact
    # t=0 value to the checked O3 future-integral decomposition and the
    # same FTC determines the constant memory and Hfuture, not Pabs/U².
    zero('same_source_absolute_pressure_backward_seam_transport',
        Pabs,Pabs.subs(t,0)+Ua**2*C*C*(1-s.exp(-t))/2)
    zero('same_source_absolute_pressure_seam',Pabs.subs(t,0),-Ua**2*C*C*H0/2+Pmemory)
    mu,Tw=s.symbols('same_mu same_Tw',positive=True);p=1+2*mu
    Pv=s.Function('same_actual_signed_Rv')(z)
    U1=asts.evaluate(asts.expression('current_O3_transition_direction_operator','actual_variable_transition_theorem','U1'),dict(s=s,Ua=Ua,mu=mu))
    memory=asts.evaluate(asts.expression('current_O3_transition_direction_operator','actual_variable_transition_theorem','Pmemory'),
        dict(s=s,U1=U1,Pv=Pv,C=C,p=p,mu=mu,Tw=Tw))
    zero('same_bound_original_absolute_pressure_memory',memory,
        Ua**2*s.exp(-1-mu)*(Pv+C*C/(2*p))*s.exp(-p*(13/mu+Tw)))
    at_seam=asts.evaluate(asts.expression('current_O3_transition_direction_operator','actual_variable_transition_theorem','Pabs'),
        dict(U=Ua,C=C,Hfuture=H0,Pmemory=memory))
    zero('same_bound_original_absolute_pressure_seam',at_seam,-Ua**2*C*C*H0/2+memory)
    datum=s.Function('same_original_axis_pressure_datum')(z)
    old_at_seam=dict(m=Ma*z,h=Ua*(1-da)*C,k=Ua*(Ma-4*da)*z*C,e=EZa*z*z+EQa*C*C,p=at_seam-datum)
    seam_env=dict(c=c,old=old_at_seam,u1=Ua*C,d=1,d3=1,K=dict(theta=0,energy=0,pressure=0),square=lambda v:v*v)
    transition_seam=asts.evaluate(asts.expression('pre_pulse_mixed_C4','slope_mu','hist'),seam_env)
    for name,value in old_at_seam.items():zero('same_actual_buffer_transition_history_seam/'+name,transition_seam[name],value)
    zero('same_actual_buffer_transition_absolute_pressure_seam',transition_seam['p']+datum,at_seam)
    parameter_theorem=parameters.exact_theorem()
    if parameter_theorem['identities'].get('actual_original_logAd_exact_endpoint') is not True:
        raise ValueError('Original canonical Ua=exp(logAd) source identity required')
    checks['same_source_Ua_exp_logAd_exact_endpoint']=True
    asts.expression('current_O3_finite_frequency_profiles','cutoff_rows','left',wanted='sigma_jets(c,t+2)')
    asts.expression('current_O3_finite_frequency_profiles','cutoff_rows','right',wanted='sigma_jets(c,4*t-1)')
    asts.expression('outer_initial','stable_sigma','derivative',
        wanted='value*(b/(a+b))*(2/x**3+2/(1-x)**3)')
    checks['actual_negative_buffer_chi_equals_monotone_sigma_offset_plus2']=True
    # Reuse checked exact nonnegative rational geometry identities, not
    # the O3 interval domain theorem itself.
    from lei_ren_part1_paper_compliant_current_O3_power_cone_operator import original_correlated_O2_O3_theorem
    geometry=original_correlated_O2_O3_theorem()
    for key in ('actual_nonnegative_M_coefficient','actual_positive_deficit_polynomial',
        'actual_deficit_shape_numerator_lower_positive_polynomial','actual_deficit_shape_denominator_below_four',
        'actual_AMGM_nonnegative_square'):
        if geometry['identities'].get(key) is not True:raise ValueError('Original O2 geometry proof omitted')
        checks['shared_exact_geometry/'+key]=True
    weighted=o3.exact_theorem()
    hashes={**geometry['input_hashes'],**weighted['input_hashes'],**parameter_theorem['input_hashes'],**asts.hashes,
        Path(__file__).name:sha(Path(__file__).name)}
    return dict(passed=True,identities=checks,input_hashes=hashes,
        negative_buffer_source_domain=(-2,0),same_actual_canonical_O2_inlet=True,
        canonical_history=dict(U='Ua*exp(-t/2)',D='da*exp(-t)',M='Ma*exp(-t)',
            h_over_U='1-D',k_over_U='M-4D',AZ='EZa/Ua²',AQ='EQa/Ua²-t/2'),
        future_pressure='H(t)=1+(H(0)-1)*exp(t); 1/(1+2mu)<=H<=1',
        exact_weighted_cone_and_all_N_monotone_error_theorem=weighted,
        strict_primitive_floor='vs-2 >= mu*chi(t)^2/4 >0 on (-2,0]',
        exact_left_flat_edge_has_vs_minus2_zero=True)


class CurrentModifiedO2TaperCone:
    def __init__(self,require_checked=True):
        self.o3=o3.CurrentModifiedO3TransitionCone();self.errors=self.o3.errors;self.ctx=self.o3.ctx
        self.data,self.mu,self.read=self.o3.data,self.o3.mu,self.o3.read
        self.theorem=exact_theorem();self.baseline=self.baseline_bounds()
        self.logfloor=self.baseline['theta_log_lower'];self.rho0=self.baseline['weighted_axial_ratio_upper']
        self.error_rows=self.leading_error_rows(MINIMUM_N);self.proof=self.prove()
        self.hashes={**self.o3.hashes,o3.NAME:sha(o3.NAME),o3.RECEIPT:sha(o3.RECEIPT),**self.theorem['input_hashes']}
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed']:raise ValueError('Checked actual O2 taper required')
            for name,digest in record['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed actual O2 taper input: '+name)

    def baseline_bounds(self):
        c=self.ctx;b=self.o3.base;read=self.read;d=self.data;positive={}
        def require(name,value):
            if endpoints(value)[0]<=0:raise ArithmeticError('O2 baseline bound unresolved: '+name)
            positive[name]=value
        da=read(c,b['actual_da']);cd=read(c,b['cD']);cx=read(c,b['cX']);delta=self.errors.delta
        H0bounds={key:read(c,value) for key,value in b['original_future_pressure_H_bounds'].items()}
        if H0bounds['upper']._mpi_!=c.mpf(1)._mpi_ or endpoints(H0bounds['lower'])[0]<=0:
            raise ValueError('Same checked original Hfuture(0) positive and <=1 required')
        require('source_bound_future_H0_positive',H0bounds['lower'])
        dcap=delta/(2*(1-delta));rlog=c.ln(2)-d['radius_logs']['logRd']+2-c.ln(cd*da)
        require('radial_relative_log_gap',-1000-rlog)
        rad=cd*da*c.exp(-1000);reserve=cd*da-dcap-rad;floor=cd*da+reserve
        require('positive_theta_reserve',reserve);require('canonical_theta_positive_lower',floor)
        ua=read(c,b['actual_Ua']);eqa=read(c,b['full_energy_EQa']);eza=read(c,b['full_energy_EZa'])
        aq=upper(c,abs(eqa/ua**2)+1);az=upper(c,abs(eza/ua**2))
        lz=upper(c,((2+2*delta)*az+(4+2*delta)*aq+3+delta)/(1-delta))
        elog=d['log_mu']+2*d['logP']+2*d['logAd']+2*c.ln(lz)-c.ln(2*cd*cx)-read(c,b['actual_da_log'])
        require('full_energy_relative_direction_log_gap',-1000-elog)
        # The same absolute memory is constant along the original buffer.
        # U(t)>=Ua, so the O3 memory bound with its U-minimum discarded is
        # conservative here; replace only its canonical theta floor.
        plog=read(c,b['full_pressure_memory_log_upper'])+c.ln(read(c,b['full_theta_uniform_lower']))-c.ln(floor)
        require('full_absolute_pressure_memory_log_gap',-1000-plog)
        direction=2*self.mu*(c.sqrt(c.exp(-1000)/(2*self.mu))+c.exp(-1000))**2
        return dict(offset_domain=(-2,0),positive_margins=positive,Ua=ua,Ma=read(c,b['actual_Ma']),da=da,
            cD=cd,cX=cx,delta_absolute_cap=dcap,radial_relative_log_upper=rlog,
            radial_error_upper=rad,canonical_theta_positive_reserve=reserve,canonical_theta_lower=floor,
            theta_lower_formula='Theta >= cX*Z²+cD*da*exp(-t)+positive_reserve',
            theta_log_lower=d['logAd']+c.ln(floor),full_energy_EQa=eqa,full_energy_EZa=eza,
            full_energy_AQ_upper=aq,full_energy_AZ_upper=az,energy_future_pressure_Z_factor_upper=lz,
            checked_original_future_H0_bounds=H0bounds,
            negative_buffer_future_H_bounds=H0bounds,
            original_Hfuture_is_future_integral_not_absolute_pressure_over_U_squared=True,
            full_energy_relative_direction_log_upper=elog,absolute_pressure_memory_relative_log_upper=plog,
            full_directional_expression_upper=direction,weighted_axial_ratio_upper=upper(c,c.sqrt(direction/2)),
            same_original_O2_pressure_energy_moments_and_all_stress_sectors=True,
            original_O3_positive_offset_baseline_not_extended_without_new_source_identity=True)

    def leading_error_rows(self,N):
        c=self.ctx;q,caps=self.errors.inputs('O2',(-2,0),N);result={}
        for label,parts in self.errors.model['stress'].items():
            result[label]={}
            for name,p in parts.items():
                cap=o3.errors.bound_expression(c,p['expression'],caps,c.mpf([-1,1]),self.errors.delta,0,p['mode'][0])[0,0]
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
            if endpoints(value)[0]<=0:raise ArithmeticError('O2 taper inequality unresolved: '+name)
            positive[name]=value
        for label,parts in self.error_rows.items():
            for name,p in parts.items():require('native_error_log_gap/'+label+'/'+name,-1010-p['log_absolute_error_over_baseline_theta_lower'])
            require('sum_native_errors_below_exp_minus1000/'+label,10-c.ln(len(parts)))
        eps=c.exp(-1000);rho=(self.rho0+eps)/(1-eps)
        alo=c.mpf('1.997');ahi=c.mpf('2.003');H=c.mpf(4);G=c.mpf(3);mumax=c.mpf('.001')
        direction=alo-H*rho;W=2*alo-mumax*H*H-2*ahi*H*rho-ahi*G*rho*rho
        dlower=direction/ahi;qlower=alo*alo*W/ahi**3
        for name,value in dict(positive_modified_theta=1-eps,weighted_axial_ratio_below_exp_minus390=c.exp(-390)-rho,
            weighted_H_over_theta=direction,weighted_W_over_theta_squared=W,
            original_signed_D_over_theta=dlower,original_signed_Q_over_theta_squared=qlower).items():require(name,value)
        return dict(minimum_integer_N=MINIMUM_N,all_finite_integer_N_at_least_minimum=True,
            strict_signed_cone_offset_domain='(-2,0]',whole_direction_and_alignment_offset_domain=(-2,0),
            whole_Z_domain=(-1,1),all_actual_phases=True,positive_margins=positive,
            baseline_theta_log_lower=self.logfloor,relative_theta_error_upper=eps,relative_weighted_axial_error_upper=eps,
            modified_weighted_axial_ratio_upper=rho,actual_a_lower=alo,actual_a_upper=ahi,
            actual_g_absolute_upper=G,actual_h_absolute_upper=H,
            actual_vs_minus2_source_lower='mu*chi(offset)^2/4; positive at every offset>-2',
            no_positive_constant_primitive_floor_over_whole_open_taper_claimed=True,
            weighted_H_over_theta_lower=direction,reduced_weighted_quadratic_over_theta_squared_lower=W,
            original_signed_D_over_theta_lower=dlower,original_signed_Q_over_theta_squared_lower=qlower,
            exact_flat_left_endpoint_has_a2_bs0_vs2=True,
            left_endpoint_strict_cone_not_claimed=True,original_pressure_energy_moments_and_all_cross_terms_retained=True,
            current_modified_O2_open_taper_signed_two_vector_cone_certified=True,**{key:False for key in OPEN})

    def query(self,offset,N):
        parameters.positive_integer_N(N);c=self.ctx;t=c.mpf(offset);lo,hi=endpoints(t)
        if N<MINIMUM_N or lo < -2 or hi>0 or not all(mp.isfinite(v) for v in (lo,hi)):
            raise ValueError('Finite O2 offset subset[-2,0], integer N>=22 required')
        p=self.o3.primitive.query(t,N);strict=lo>-2
        # Generic wide cutoff jets can enclose chi=0 even on a strictly
        # positive box. On this negative interval chi=sigma(t+2) is
        # monotone; its left endpoint yields a directed primitive floor.
        # The ordinary-jet helper deliberately uses a zero lower tail
        # cap below exp(-1000). Evaluate the equivalent defining positive
        # fraction here, rather than treating its numerical cap as chi.
        chi_left=stable_sigma(c,c.mpf(lo)+2)[0]
        floor=self.mu*chi_left**2/4
        return dict(offset=t,finite_integer_N=N,actual_phase=c.mpf(N)*t,
            signed_current_g_enclosure=p['finite_N_theta_shear_excess']/self.mu,
            signed_current_h_enclosure=p['finite_N_axial_shear']/c.sqrt(self.mu),
            actual_vs_minus2_source_lower=floor,
            generic_cutoff_jet_primitive_lower_enclosure=p['source_certified_joint_margin_lower_envelope'],
            monotone_negative_cutoff_left_endpoint_used_for_directed_floor=True,
            strict_primitive_positive_by_source_for_entire_query_box=strict,
            closed_box_contains_degenerate_left_endpoint=lo==-2,
            complete_modified_signed_cone_certified_for_entire_query_box=strict,
            whole_direction_and_alignment_certificate=self.proof,source_family=p['source_family'],
            **{key:False for key in OPEN})


def run():
    field=CurrentModifiedO2TaperCone(require_checked=False)
    result=dict(source_family=field.data['source']['accepted']['source_family'],
        exact_original_O2_baseline_and_weighted_source_theorem=field.theorem,
        original_whole_closed_O2_taper_baseline_bounds=field.baseline,
        whole_modified_O2_taper_cone=field.proof,complete_leading_stress_error_log_ledger=field.error_rows,
        examples={name:field.query(t,N) for name,t,N in (
            ('closed_direction_domain',(-2,0),22),('exact_degenerate_left_edge','-2',22),
            ('very_small_left_taper','-1.999',22),('whole_positive_taper',('-1.99','0'),22),
            ('old_seam','0',22),('current_frequency',('-1.999','0'),10**12))},
        current_modified_O2_open_taper_signed_two_vector_cone_certified=True,
        scoped_source_cone_does_not_increment_unmodified_registry_counts=True,
        **{key:False for key in OPEN},input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(parameters.encoded(result),indent=2)+'\n').encode())
    print('Modified O2 open taper cone generated; all phases/Z/offset(-2,0], all integer N>=22; flat edge explicit',flush=True)
    return result


if __name__=='__main__':run()
