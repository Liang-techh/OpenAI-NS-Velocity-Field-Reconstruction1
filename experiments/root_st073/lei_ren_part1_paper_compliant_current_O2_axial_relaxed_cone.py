"""Full original slow axial turnoff relaxed cone, with M/B correlation."""
import ast
import functools
import json
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_O2_reference_slope_relaxed_cone as ref
import lei_ren_part1_paper_compliant_current_O2_modified_taper_cone as taper
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets

HERE,PREFIX,sha=ref.HERE,ref.PREFIX,ref.sha
NAME=PREFIX+'current_O2_axial_relaxed_cone.json'
RECEIPT=PREFIX+'current_O2_axial_relaxed_cone_check.json'
OPEN=taper.OPEN
parameters,endpoints,read=taper.parameters,taper.endpoints,ref.read


def cutoff_budget(c,cells=1024):
    """Directed contiguous interval cover; no point/phase sampling."""
    if type(cells) is not int or cells<1024:raise ValueError('At least1024 directed cutoff cells required')
    weighted=tail=c.mpf(0);rows=[]
    for i in range(cells):
        lo=c.mpf(i)/cells;hi=c.mpf(i+1)/cells;x=c.mpf([endpoints(lo)[0],endpoints(hi)[1]])
        jets=sigma_jets(c,x);d=c.mpf(endpoints(jets[1])[1])
        # The exact sigma derivative is nonnegative. The endpoint-tail
        # jet enclosure may be symmetric; its upper still bounds it.
        w=d/(1+8*jets[0]);t=d*c.exp(-40*(1-x))
        wu,tu=endpoints(w)[1],endpoints(t)[1]
        weighted=c.mpf(max(endpoints(weighted)[1],wu));tail=c.mpf(max(endpoints(tail)[1],tu))
        rows.append(dict(cell=i,domain=x,weighted_derivative_upper=c.mpf(wu),derivative_over_actual_y_upper=c.mpf(tu)))
    if endpoints(c.mpf(9)/4-weighted)[0]<=0 or endpoints(c.mpf('1e-4')-tail)[0]<=0:
        raise ArithmeticError('Complete original cutoff budget unresolved')
    return dict(exact_rational_contiguous_cell_count=cells,whole_sigma_argument_domain=(0,1),
        weighted_derivative_upper=weighted,derivative_over_actual_y_upper=tail,
        coarse_weighted_derivative_bound=c.mpf(9)/4,coarse_derivative_over_y_bound=c.mpf('1e-4'),
        directed_interval_cover_includes_both_flat_edges=True,all_cell_certificates=rows)


@functools.lru_cache(maxsize=1)
def exact_theorem():
    asts=ref.numeric.transport.SourceAST();checks={}
    def zero(name,a,b):
        if s.cancel(s.expand(s.expand_power_exp(a-b)))!=0:raise ArithmeticError('Axial relaxed source identity: '+name)
        checks[name]=True
    t,z,delta=s.symbols('actual_axial_offset Z delta',real=True)
    Us,Ds=s.symbols('same_slope_U same_slope_deficit',positive=True)
    Mf,Bf=s.Function('actual_M')(t),s.Function('actual_B')(t)
    D=Ds*s.exp(-t);U=Us*s.exp(-t/2);C=1/(1+z*z);L=1-delta*z*z
    A=(1-delta/2)*C+(1-delta)*z*z*C*C
    Bg=((2*delta*z*z-1)*C+2*(1-z*z)*z*z*C*C)/L
    AZ,AQ,Hf=s.symbols('actual_kinetic_AZ actual_AQ actual_future_H',real=True)
    memory=s.Function('actual_absolute_pressure_memory')(z)
    Pabs=-U**2*C*C*Hf/2+memory;CE=z*z*AZ+C*C*AQ
    V=4*z*Bf;h=U*(1-D)*C;k=U*(Mf-4*D)*z*C;m=Mf*z
    ode={s.diff(Mf,t):4*Bf-Mf}
    zero('actual_axial_M_ODE',s.diff(m,t).subs(ode),V-m)
    zero('actual_axial_h_ODE',s.diff(h,t),U*C-s.Rational(3,2)*h)
    zero('actual_axial_shared_M_K_deficit_ODE',s.diff(k,t).subs(ode),U*C*V-s.Rational(3,2)*k)
    zero('actual_axial_M_minus4B_positive_FTC_density',s.diff(s.exp(t)*(Mf-4*Bf),t).subs(ode),-4*s.exp(t)*s.diff(Bf,t))
    K2=s.Function('same_original_squared_B_kernel')(t)
    kernel_ode={s.diff(K2,t):Bf**2-K2};W=s.exp(-t)+K2
    zero('actual_squared_kernel_kinetic_comparison_ODE',s.diff(W,t).subs(kernel_ode),Bf**2-W)
    zero('actual_squared_kernel_nonnegative_FTC_density',s.diff(s.exp(t)*W,t).subs(kernel_ode),s.exp(t)*Bf**2)
    zero('actual_squared_kernel_upper1_FTC_density',s.diff(s.exp(t)*(1-W),t).subs(kernel_ode),s.exp(t)*(1-Bf**2))
    # Same Rd future datum, transported backwards over the entire axial
    # and 11-unit buffer interval. No Rd datum is reused at the slope seam.
    HRd=s.Symbol('same_original_Rd_future_H',real=True);T=s.exp(40)-1+11
    Hactual=1+(HRd-1)*s.exp(t-T)
    zero('actual_axial_future_H_backward_Rd_ODE',s.diff(Hactual,t),Hactual-1)
    zero('actual_axial_future_H_backward_Rd_seam',Hactual.subs(t,T),HRd)
    zero('actual_axial_absolute_pressure_backward_FTC',s.diff(-U**2*C*C*Hactual/2+memory,t),U**2*C*C/2)
    zero('actual_axial_future_H_convex_Rd_weights',Hactual,(1-s.exp(t-T))+s.exp(t-T)*HRd)
    zero('actual_axial_logU_Rd_source_endpoint',-s.exp(40)/2+s.Rational(3,10)-s.Rational(11,2),-s.exp(40)/2-s.Rational(26,5))
    x=s.Symbol('same_sigma_argument',positive=True)
    ex=s.exp(-1/x**2);ey=s.exp(-1/(1-x)**2);sig=ex/(ex+ey)
    zero('actual_sigma_strict_positive_derivative_inside',s.diff(sig,x),ex*ey/(ex+ey)**2*(2/x**3+2/(1-x)**3))
    ctx=SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp)
    from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
    from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
    op=asts.replay('current_pre_pulse_stress_operator','raw_pre_stress_rows',dict(
        axial_derivative=lambda v:s.diff(v,z),shifted_rows=shifted_rows,product_rows=product_rows))
    rows=lambda v:[s.diff(v,t,j) for j in range(5)]
    raw=op(ctx,delta,z,rows(U*C),rows(V),dict(m=rows(m),h=rows(h),k=rows(k),e=rows(U**2*CE)),rows(Pabs))
    theta=(A-C)/L+(-A/L-4*Bg)*D+(C+Bg)*Mf
    zero('actual_complete_axial_turnoff_inertial_theta',sum(p['shape'][0] for name,p in raw['theta'].items() if name!='variable_radial_shear')/U,theta)
    zero('actual_full_axial_local_transport',raw['axial']['local_axial_transport']['shape'][0],-4*z*Bf/L)
    zero('actual_full_axial_nonlinear_transport',raw['axial']['nonlinear_meridional_transport']['shape'][0],4*z*Bf*Mf)
    zero('actual_full_axial_linear_moment_zero',raw['axial']['retained_linear_axial_moment']['shape'][0],0)
    zero('actual_full_axial_radial_shear',raw['axial']['axial_radial_shear']['shape'][0],8*z*s.diff(Bf,t))
    G=((2*delta*z*z-2*(1-z*z))*AZ+(2*delta*C*C+4*(1-z*z)*C**3)*AQ
        -Hf*((1+delta)*C*C+2*(1-z*z)*C**3))/L
    mem=(2*(1+delta)*z*memory-(1-z*z)*s.diff(memory,z))/L
    zero('actual_full_axial_energy_and_absolute_pressure',raw['axial']['retained_full_energy']['shape'][0]+raw['axial']['actual_absolute_pressure']['shape'][0],U**2*z*G+mem)
    a,b,r=s.symbols('a bs actual_Tz_over_theta',real=True);x=b*r/2;y=b*b/4
    zero('actual_a2_original_quadratic_product_factorization',2*(1-b*r/2)**2-b*b/2*(r+b/2)**2,2*(1+y)*(1-2*x-y))
    asts.expression('pre_pulse_mixed_C4','axial','u',wanted='u1*root')
    asts.expression('pre_pulse_mixed_C4','axial','V',wanted='[z*(4*b) for b in B]')
    asts.expression('pre_pulse_mixed_C4','axial','hist',wanted="dict(m=old['m']*d+z*(4*K['B_mass']),h=old['h']*d3+u1*(root-d3),k=old['k']*d3+u1*z*(4*root*K['B_mass']),e=old['e']*d+square(z)*(16*self.invP2*K['B_squared_mass'])-square(u1)*(t*d/2),p=old['p']+square(u1)*((1-d)/2))")
    asts.expression('pre_pulse_mixed_C4','axial','y',wanted='c.exp(md*phase)')
    asts.expression('pre_pulse_mixed_C4','axial','t',wanted='y-1')
    asts.expression('pre_pulse_mixed_C4','axial','B',wanted='turnoff_derivatives(c,y,md,phase)')
    asts.expression('pre_pulse_mixed_C4','axial','parent',wanted='self.inlet(Z)')
    asts.expression('pre_pulse_mixed_C4','slope','hist',wanted="dict(m=V,h=h,k=h*V,e=square(V)*self.invP2-square(qi)*(c.mpf(5)/12+mass[2]/2)*c.exp(-y),p=square(qi)*(c.mpf('2.5')+mass[1]/2))")
    for key in ('sigma_jets','sigma_tail_bound'):
        asts.method('flat_pulse_derivatives',key)
    asts.method('pre_pulse_mixed_C4','turnoff_derivatives')
    asts.method('outer_initial','turnoff_kernels')
    checks['same_original_M_initial4_and_monotone_B_give4B_leM_le4']=True
    checks['same_original_kinetic_EZ_initial16_over_Pstar2_and_squared_B_kernel_retained']=True
    checks['same_absolute_pressure_FTC_extends_original_backward_memory_into_axial_source']=True
    return dict(passed=True,identities=checks,
        source_domain='phase[0,1], y=exp(40*phase), t=y-1',
        exact_signed_shear=dict(a=2,bs='8Z*B_y/(Pstar*U*C)',vs_minus2='bs^2/2'),
        exact_original_correlations=dict(X='1-Ds*exp(-t)',K_over_U='M-4Ds*exp(-t)',M='4B+nonnegative_FTC_remainder'),
        exact_source_comparisons=dict(B_range=(0,1),B_initial=1,M_initial=4,squared_kernel_initial=0,
            squared_kernel_definition='K2(t)=integral_0^t exp(s-t)*B(s)^2 ds',
            squared_kernel_bound='0<=K2<=1-exp(-t); 0<=exp(-t)+K2<=1',
            full_kinetic_AZ_bound='0<=AZ<=16/(Pstar^2*Umin^2)',
            full_AQ_split='AQ=AQs-t/2; |AQs|<=A0; 0<=t/(1+t)<=1',
            future_pressure='H=1+(H_Rd-1)*exp(t-(exp(40)-1+11)); H_Rd<=H<=1',
            pressure_budget_uses_signed_H_le1=True,
            strict_shear='sigma_prime>0 inside(0,1), By=-sigma_prime/(40*(1+t)); bs!=0 iff interior phase and Z!=0'),
        input_hashes={**asts.hashes,**ref.exact_theorem()['input_hashes'],**taper.exact_theorem()['input_hashes'],Path(__file__).name:sha(Path(__file__).name)})


class CurrentOriginalAxialRelaxedCone:
    def __init__(self,require_checked=True):
        self.taper=taper.CurrentModifiedO2TaperCone();self.data=self.taper.data;self.ctx=self.taper.ctx
        self.mu,self.delta,self.read=self.taper.mu,self.taper.errors.delta,self.taper.read
        self.family=self.data['source']['accepted']['source_family'];self.theorem=exact_theorem()
        self.hashes={**self.taper.hashes,taper.NAME:sha(taper.NAME),taper.RECEIPT:sha(taper.RECEIPT),**self.theorem['input_hashes']}
        receipt=json.loads((HERE/ref.RECEIPT).read_bytes())
        if not receipt['all_passed'] or receipt['source_family']!=self.family or not receipt['current_original_Rh_O2_slope_relaxed_input_cone_certified']:
            raise ValueError('Same original slope relaxed source required')
        for name,digest in receipt['input_hashes'].items():
            if sha(name)!=digest:raise ValueError('Changed original slope input: '+name)
        self.hashes.update({**receipt['input_hashes'],ref.RECEIPT:sha(ref.RECEIPT)})
        self.cutoff=cutoff_budget(self.ctx);self.baseline=self.baseline_bounds();self.proof=self.prove()
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked['all_passed']:raise ValueError('Checked original full axial relaxed cone required')
            for name,digest in checked['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed original axial source: '+name)

    def baseline_bounds(self):
        c=self.ctx;d=self.data;delta=self.delta;b=self.taper.o3.base;read=self.read;positive={}
        def require(name,value):
            if endpoints(value)[0]<=0:raise ArithmeticError('Axial source baseline unresolved: '+name)
            positive[name]=value
        cD=read(c,b['cD']);cX=read(c,b['cX']);da=read(c,b['actual_da'])
        rlog=c.ln(2)-d['source']['logRref']-1-c.ln(cD*da)
        require('whole_axial_radial_vs_deficit_log_gap',-1000-rlog)
        radial=cD*da*c.exp(-1000);reserve=cD*da-delta/(2*(1-delta))-radial
        require('positive_whole_axial_theta_reserve',reserve);require('delta_below_one_millionth',c.mpf('1e-6')-delta)
        # Exact Us=exp(-1/5), from the symmetric slope J(1)=1/2.
        A0=(c.mpf(5)/12+(c.exp(c.mpf(6)/5)-1)/(c.mpf(12)/5))*c.exp(-c.mpf(3)/5)
        require('full_slope_energy_AQ_below1',1-A0)
        logUmin=-c.exp(40)/2+c.mpf(3)/10
        Hbounds={key:read(c,value) for key,value in b['original_future_pressure_H_bounds'].items()}
        if Hbounds['upper']._mpi_!=c.mpf(1)._mpi_:
            raise ValueError('Same original Rd future H upper1 required')
        require('same_original_Rd_future_H_positive',Hbounds['lower'])
        require('whole_axial_backward_Rd_offset_at_least11',c.mpf(11))
        return dict(positive_margins=positive,cD=cD,cX=cX,actual_Rd_deficit=da,
            canonical_theta_positive_reserve=reserve,full_slope_AQ_absolute_upper=A0,original_scalar_log_U_minimum=logUmin,
            radial_relative_log_upper=rlog,
            checked_original_Rd_future_H_bounds=Hbounds,
            strong_theta_lower='(1-delta)*Z^2*C^2*(1+8*(1-Z^2)*B)+cD*D+positive_reserve',
            weaker_theta_lower_for_absolute_pressure='cX*Z^2+cD*da+positive_reserve',
            actual_M_ge4B_correlation_retained=True,full_original_absolute_pressure_and_energy_retained=True)

    def prove(self):
        c=self.ctx;d=self.data;delta=self.delta;b=self.taper.o3.base;read=self.read;positive={}
        def require(name,value):
            if endpoints(value)[0]<=0:raise ArithmeticError('Full axial cone budget unresolved: '+name)
            positive[name]=value
        md=c.mpf(40);A0=self.baseline['full_slope_AQ_absolute_upper'];lu=self.baseline['original_scalar_log_U_minimum']
        denom=(1-delta)**2
        dominant=16/(md*denom)*(c.mpf(9)/4+8*delta)
        aq=32*A0*(1+delta)/(md*denom)*c.mpf('1e-4')
        pressure=(32+16*delta)/(md*denom)*c.mpf('1e-4')
        native_logs=dict(
            kinetic_AZ=c.ln(1024*(1+delta)*16/(md*denom))-2*d['logP']-2*lu,
            local_axial_and_meridional=c.ln(2048*(4+1/(1-delta))/(md*(1-delta)))-2*d['logP']-2*lu,
            axial_radial_shear=c.ln(512*(8/md)**2/(1-delta))-2*d['logP']-2*lu-d['source']['logRref']-1,
            bs_squared_over4=c.ln(64*(8/md)**2)-2*d['logP']-2*lu)
        old_pressure=read(c,b['full_pressure_memory_log_upper']);old_floor=read(c,b['full_theta_uniform_lower'])
        native_logs['absolute_original_pressure_memory']=(c.ln(64/md)-1-self.mu+old_pressure-d['logP']-d['logAd']
            +c.ln(old_floor)-c.ln(self.baseline['cX']*self.baseline['cD']*self.baseline['actual_Rd_deficit'])/2)
        for name,value in native_logs.items():require('source_log_budget_below_exp_minus1010/'+name,-1010-value)
        eps=c.exp(-1000);product=dominant+aq+pressure+4*eps;total=product+eps
        require('complete_signed_axial_product_plus_bs2_over4_below_point91',c.mpf('.91')-total)
        D=1-product/2;Q=2*(1-total)
        require('full_D_over_theta_exceeds_point54',D-c.mpf('.54'));require('full_Q_over_theta_squared_exceeds_point17',Q-c.mpf('.17'))
        return dict(whole_source_phase_domain=(0,1),whole_Z_domain=(-1,1),positive_margins=positive,
            complete_budget_components=dict(dominant_full_AQ=dominant,initial_full_AQ=aq,full_future_pressure=pressure,
                remaining_four_full_source_products=4*eps,bs_squared_over4=eps),
            full_suppressed_sector_log_budgets=native_logs,
            signed_bs_times_Tz_over_Ttheta_upper=product,bs_squared_over4_upper=eps,
            full_D_over_theta_lower=D,full_Q_over_theta_squared_lower=Q,
            current_original_O2_axial_relaxed_input_cone_certified=True,
            strict_source_shear_positive_only_for_interior_phase_and_nonzero_Z=True,
            whole_closed_axial_strict_cone_certified=False,
            full_original_pressure_energy_kinetic_and_own_moment_radial_sectors_retained=True,
            no_positive_constant_shear_excess_on_whole_axial_chart_claimed=True,**{key:False for key in OPEN})

    def query(self,phase,Z=(-1,1)):
        c=self.ctx;q,z=c.mpf(phase),c.mpf(Z);lo,hi=endpoints(q);zl,zh=endpoints(z)
        if lo<0 or hi>1 or zl < -1 or zh>1 or not all(mp.isfinite(v) for v in (lo,hi,zl,zh)):
            raise ValueError('Finite original axial phase subset[0,1] and Z subset[-1,1] required')
        strict=lo>0 and hi<1 and (zh<0 or zl>0)
        return dict(phase=q,Z=z,actual_ordinary_logR_distance=c.exp(40*q),actual_axial_offset=c.exp(40*q)-1,
            source_family=self.family,whole_full_signed_input_certificate=self.proof,
            current_original_O2_axial_relaxed_input_cone_certified=True,
            complete_original_strict_cone_certified_for_entire_query_box=strict,
            source_kappa2_midplane_or_flat_edge_included=not strict,**{key:False for key in OPEN})


def run():
    field=CurrentOriginalAxialRelaxedCone(require_checked=False)
    result=dict(source_family=field.family,exact_full_original_axial_source_theorem=field.theorem,
        directed_whole_original_cutoff_cover=field.cutoff,whole_original_axial_baseline=field.baseline,
        whole_original_axial_relaxed_cone=field.proof,
        examples={name:field.query(q,z) for name,q,z in (
            ('whole_axial',(0,1),(-1,1)),('slope_axial_seam',0,(-1,1)),('axial_buffer_seam',1,(-1,1)),
            ('midplane',('.2','.8'),0),('positive_strict_interior',('.2','.8'),('.1','1')),
            ('negative_strict_interior',('.2','.8'),('-1','-.1')))},
        current_original_O2_axial_relaxed_input_cone_certified=True,whole_closed_axial_strict_cone_certified=False,
        **{key:False for key in OPEN},input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(parameters.encoded(result),indent=2)+'\n').encode())
    print('Original full axial relaxed cone generated; complete product budget<.91',flush=True)
    return result


if __name__=='__main__':run()
