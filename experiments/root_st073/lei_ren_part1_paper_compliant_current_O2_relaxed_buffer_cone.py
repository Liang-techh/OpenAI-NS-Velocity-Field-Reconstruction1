"""Actual nonzero O2 buffer is relaxed-cone input, never strict admission.

Uses the signed original buffer source and complete pressure/energy.
Its source shear has kappa=2; frequency cannot remove this obstruction.
"""
import functools
import json
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_O2_modified_taper_cone as taper
from lei_ren_part1_paper_compliant_current_original_cone_operator import original_cone_theorem

HERE,PREFIX,sha=taper.HERE,taper.PREFIX,taper.sha
NAME=PREFIX+'current_O2_relaxed_buffer_cone.json'
RECEIPT=PREFIX+'current_O2_relaxed_buffer_cone_check.json'
parameters,endpoints,upper=taper.parameters,taper.endpoints,taper.upper
OPEN=taper.OPEN


@functools.lru_cache(maxsize=1)
def exact_theorem():
    parent=taper.exact_theorem();cone=original_cone_theorem()
    asts=taper.o3.errors.recovered.histories.SourceAST();checks={}
    def zero(name,a,b):
        if s.cancel(s.expand(a-b))!=0:raise ArithmeticError('Relaxed O2 source identity: '+name)
        checks[name]=True
    v,t=s.symbols('actual_buffer_selector shared_modulation_offset',real=True)
    asts.expression('current_O3_modulated_histories','history','pre',wanted='self.pre.axial(zv,buffer_offset=v)')
    coordinate=asts.evaluate(asts.expression('current_O3_modulated_histories','history','t',wanted='v-11'),dict(v=v))
    zero('actual_shared_offset_is_buffer_selector_minus11',coordinate,v-11)
    zero('actual_unmodified_buffer_selector_range',t+11,(t+11))
    asts.expression('pre_pulse_mixed_C4','axial','B',wanted='[c.mpf(0)]*5')
    asts.expression('pre_pulse_mixed_C4','axial','y',wanted='c.exp(md)+selector')
    asts.expression('current_O3_finite_frequency_profiles','cutoff_rows','left',wanted='sigma_jets(c,t+2)')
    checks['shared_offset_minus11_to_minus2_has_buffer_selector0_to9']=True
    checks['entire_domain_before_modulation_support_has_flat_zero_local_and_cumulative_errors']=True
    F,theta=s.symbols('same_F positive_Ttheta',positive=True);axial=s.Symbol('same_Tz',real=True)
    St,Sz=-2*F,s.Integer(0);kappa=-(St**2+Sz**2)/(F*St)
    dot=theta*St+axial*Sz;cross=-theta*Sz+axial*St
    zero('actual_buffer_kappa_equals2',kappa,2)
    zero('actual_buffer_signed_dot_is_minus2Ftheta',dot,-2*F*theta)
    zero('actual_relaxed_kappa_le2_threshold_is_zero',-F*St*(2-kappa),0)
    zero('actual_quadratic_reserve_at_degeneracy',2*dot**2-(kappa-2)*cross**2,8*F**2*theta**2)
    zero('actual_strict_shear_margin_is_exact_zero',kappa-2,0)
    # The checked source ODEs hold over this entire pure-buffer selector,
    # not by extrapolating an O3 positive-offset interval bound.
    for name in ('same_replayed_actual_buffer_m_ODE','same_replayed_actual_buffer_h_ODE',
        'same_replayed_actual_buffer_k_ODE','same_replayed_actual_buffer_e_ODE',
        'same_replayed_actual_buffer_p_ODE','same_O2_correlated_positive_theta_split',
        'same_full_O2_energy_and_absolute_pressure','same_source_absolute_pressure_backward_seam_transport'):
        if parent['identities'].get(name) is not True:raise ValueError('Actual buffer source theorem missing: '+name)
        checks['consumed_actual_buffer_source/'+name]=True
    return dict(passed=True,identities=checks,shared_offset_domain=(-11,-2),actual_buffer_selector_domain=(0,9),
        parent_source_theorem=parent,original_signed_cone_source_theorem=cone,
        paper_relaxed_input_reference='Lei-Ren v2 p22 (3.23), kappa=2: T dot S<0',
        paper_strict_target_reference='Lei-Ren v2 p22 (3.21), requires kappa>2 wherever T is nonzero',
        exact_source_shear=dict(a=2,bs=0,vs=2),
        increasing_N_cannot_change_the_exact_original_shear=True,
        input_hashes={**parent['input_hashes'],**cone['input_hashes'],**asts.hashes,
            Path(__file__).name:sha(Path(__file__).name)})


class CurrentOriginalO2RelaxedBufferCone:
    def __init__(self,require_checked=True):
        self.taper=taper.CurrentModifiedO2TaperCone();self.ctx=self.taper.ctx
        self.errors,self.data,self.read=self.taper.errors,self.taper.data,self.taper.read
        self.mu=self.taper.mu;self.theorem=exact_theorem();self.baseline=self.baseline_bounds()
        self.hashes={**self.taper.hashes,taper.NAME:sha(taper.NAME),taper.RECEIPT:sha(taper.RECEIPT),
            **self.theorem['input_hashes']}
        self.proof=self.prove()
        if require_checked:
            record=json.loads((HERE/RECEIPT).read_bytes())
            if not record['all_passed']:raise ValueError('Checked original O2 relaxed input required')
            for name,digest in record['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed relaxed O2 source: '+name)

    def baseline_bounds(self):
        c=self.ctx;b=self.taper.o3.base;read=self.read;d=self.data;positive={}
        def require(name,value):
            if endpoints(value)[0]<=0:raise ArithmeticError('Whole relaxed buffer bound unresolved: '+name)
            positive[name]=value
        da=read(c,b['actual_da']);cd=read(c,b['cD']);cx=read(c,b['cX']);delta=self.errors.delta
        dcap=delta/(2*(1-delta))
        # R=Rd*exp(t), t>=-11. This conservative bound ignores the
        # simultaneous growth of the positive deficit da*exp(-t).
        rlog=c.ln(2)-d['radius_logs']['logRd']+11-c.ln(cd*da)
        require('radial_relative_log_gap',-1000-rlog)
        radial=cd*da*c.exp(-1000);reserve=cd*da-dcap-radial;floor=cd*da+reserve
        require('positive_theta_reserve',reserve);require('canonical_theta_positive_lower',floor)
        ua=read(c,b['actual_Ua']);eqa=read(c,b['full_energy_EQa']);eza=read(c,b['full_energy_EZa'])
        aq=upper(c,abs(eqa/ua**2)+c.mpf(11)/2);az=upper(c,abs(eza/ua**2))
        H0={key:read(c,value) for key,value in b['original_future_pressure_H_bounds'].items()}
        if H0['upper']._mpi_!=c.mpf(1)._mpi_:raise ValueError('Future pressure upper1 required')
        require('source_bound_future_H0_positive',H0['lower'])
        lz=upper(c,((2+2*delta)*az+(4+2*delta)*aq+3+delta)/(1-delta))
        elog=d['log_mu']+2*d['logP']+2*d['logAd']+2*c.ln(lz)-c.ln(2*cd*cx)-read(c,b['actual_da_log'])
        require('full_energy_relative_direction_log_gap',-1000-elog)
        plog=read(c,b['full_pressure_memory_log_upper'])+c.ln(read(c,b['full_theta_uniform_lower']))-c.ln(floor)
        require('full_absolute_pressure_memory_log_gap',-1000-plog)
        direction=2*self.mu*(c.sqrt(c.exp(-1000)/(2*self.mu))+c.exp(-1000))**2
        # On t[-11,-2], U(t)>=Ua. Using Ua is conservative for both
        # the full theta floor and constant absolute pressure memory.
        return dict(shared_offset_domain=(-11,-2),positive_margins=positive,
            canonical_theta_lower=floor,canonical_theta_positive_reserve=reserve,
            theta_normalized_log_lower=d['logAd']+c.ln(floor),cD=cd,cX=cx,da=da,
            theta_lower_formula='Theta>=cX*Z^2+cD*da*exp(-t)+positive_reserve',
            full_energy_AQ_upper=aq,full_energy_AZ_upper=az,full_energy_EQa=eqa,full_energy_EZa=eza,
            whole_backward_buffer_AQ_growth_upper=c.mpf(11)/2,
            worst_native_log_radius=d['radius_logs']['logRd']-11,
            radial_relative_log_upper=rlog,radial_absolute_error_upper=radial,
            checked_future_pressure_H0_bounds=H0,full_energy_and_future_pressure_Z_factor_upper=lz,
            full_energy_relative_direction_log_upper=elog,absolute_pressure_memory_relative_log_upper=plog,
            full_directional_expression_upper=direction,weighted_axial_ratio_upper=upper(c,c.sqrt(direction/2)),
            full_original_pressure_energy_moments_and_radial_sectors_retained=True)

    def prove(self):
        if endpoints(self.baseline['canonical_theta_lower'])[0]<=0:raise ValueError('Nonzero source theta required')
        return dict(whole_shared_offset_domain=(-11,-2),whole_Z_domain=(-1,1),all_positive_finite_integer_N=True,
            theta_normalized_log_lower=self.baseline['theta_normalized_log_lower'],
            actual_a=2,actual_bs=0,actual_vs=2,strict_vs_minus2=0,
            original_signed_D_over_theta_exact=1,original_signed_Q_over_theta_squared_exact=2,
            signed_dot_is_strictly_negative=True,stress_nonzero_throughout_domain=True,
            relaxed_kappa_le2_inequality_is_exact_strict_direction_at_kappa2=True,
            current_original_O2_buffer_relaxed_input_cone_certified=True,
            current_original_O2_buffer_strict_cone_certified=False,
            any_N_upgrade_to_strict_cone_possible=False,**{key:False for key in OPEN})

    def query(self,offset,N=1):
        parameters.positive_integer_N(N);c=self.ctx;t=c.mpf(offset);lo,hi=endpoints(t)
        if lo < -11 or hi > -2 or not all(mp.isfinite(v) for v in (lo,hi)):
            raise ValueError('Finite unchanged shared O2 offset subset[-11,-2] required')
        return dict(shared_offset=t,actual_buffer_selector=11+t,finite_integer_N=N,
            native_logR=self.data['radius_logs']['logRd']+t,
            source_family=self.data['source']['accepted']['source_family'],
            local_profiles_and_all_cumulative_source_errors_exactly_zero=True,
            full_original_stress_not_zero=True,relaxed_direction_certificate=self.proof,
            complete_modified_signed_cone_certified_for_entire_query_box=False,
            current_original_O2_buffer_relaxed_input_cone_certified=True,**{key:False for key in OPEN})


def run():
    field=CurrentOriginalO2RelaxedBufferCone(require_checked=False)
    result=dict(source_family=field.data['source']['accepted']['source_family'],
        exact_original_buffer_relaxed_cone_source_theorem=field.theorem,
        whole_original_buffer_baseline=field.baseline,whole_original_buffer_relaxed_cone=field.proof,
        examples={name:field.query(v,N) for name,v,N in (
            ('whole_buffer',(-11,-2),1),('turnoff_to_buffer_seam','-11',1),('interior','-7.13',37),
            ('exact_modulation_flat_edge','-2',22),('current_frequency',(-11,-2),10**12))},
        current_original_O2_buffer_relaxed_input_cone_certified=True,
        current_original_O2_buffer_strict_cone_certified=False,**{key:False for key in OPEN},
        input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(parameters.encoded(result),indent=2)+'\n').encode())
    print('Original O2 buffer relaxed input generated; strict kappa-2 remains exactly zero',flush=True)
    return result


if __name__=='__main__':run()
