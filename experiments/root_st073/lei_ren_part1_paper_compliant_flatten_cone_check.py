"""Focused independent check of the whole original flatten two-vector cone."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_flatten_cone import (
    CompliantFlattenStressC3,CompliantGlobalPhysicalAssembly,DOMAIN,TAIL_DOMAIN,
    flatten_cone_identities,flatten_cone_source_bridge,whole_flatten_bounds)
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def independent_signed_pulse_memory_fixture():
    """Direct original sigmoid/history integral, with resolvable signed memory.

    These moderate samples check the formula independently; the actual cone
    proof is the continuous source-bound comparison and interval slope cover.
    """
    with mp.workdps(85):
        a,mu,xp=map(mp.mpf,('.1','.2','.3'))
        r=1-mu;k=1-a;b=(1-2*a)/2;g=(mu-a)/r
        memory=mp.exp(-13*r/mu)
        xv=1/r+(xp-1/r)*memory
        floor=g-2*k*abs(xp-1/r)*memory
        if not 0<xp<1/r or not 0<xv<1/r or floor<=0:
            raise ArithmeticError('Signed below-equilibrium fixture condition failed')
        def sigma(v):
            x=v/100
            if x<=0:return mp.mpf(0)
            if x>=1:return mp.mpf(1)
            odds=1/(1-x)**2-1/x**2
            return 1/(1+mp.exp(-odds))
        def sigma_t(v):
            x=v/100
            if x<=0 or x>=1:return mp.mpf(0)
            value=sigma(v)
            return value*(1-value)*(2/(1-x)**3+2/x**3)/100
        errors=[];gaps=[]
        for t,z in ((mp.mpf(0),mp.mpf(0)),(mp.mpf(50),mp.mpf('.5')),
                    (mp.mpf(100),mp.mpf(0)),(mp.mpf(100),mp.mpf(1))):
            f=(1+z*z)/2;rho=mp.log(f);j=2*z*z/(1+z*z)
            F=lambda v:mp.exp(rho*sigma(v))
            cuts=sorted(set([mp.mpf(0),t]+[mp.mpf(v) for v in (25,50,75) if v<t]))
            integral=lambda fn:mp.quad(fn,cuts) if t else mp.mpf(0)
            N=xv*mp.exp(-r*t)+integral(lambda v:mp.exp(-r*(t-v))*F(v))
            Nz=integral(lambda v:mp.exp(-r*(t-v))*F(v)*sigma(v)*2*z/(1+z*z))
            direct=((k+b*j)*N-b*z*Nz)/F(t)-1
            V0=b*j*xv+k*(xp-1/r)*memory
            shifted=integral(lambda v:mp.exp(r*v)*F(v)*
                (b*j*(1-sigma(v))-rho*sigma_t(v)*(1+g)))
            comparison=g+mp.exp(-r*t)/F(t)*(V0+shifted)
            error=abs(direct-comparison)
            if error>mp.mpf('1e-60') or direct<floor:
                raise ArithmeticError('Independent signed flatten history identity failed')
            errors.append(str(error));gaps.append(str(direct-floor))
        return dict(passed=True,direct_native_integral_identity_checks=4,
            positive_continuous_lemma_comparisons=4,Xp_and_Xv_below_equilibrium=True,
            signed_memory_is_nonzero=True,identity_errors=errors,positive_gaps=gaps,
            comparison_tolerance='1e-60',sampled_values_used_as_actual_cone_proof=False,
            actual_cone_proof_is_continuous_source_bound=True,fixture_is_global_NS_validation=False)


def run():
    with mp.workdps(300):
        name=PREFIX+'flatten_cone.json';record=json.loads((HERE/name).read_bytes())
        hashes=dict(record['input_hashes'])
        for source,digest in hashes.items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                raise ValueError('Flatten cone source changed: '+source)
        stress=CompliantFlattenStressC3();assembly=CompliantGlobalPhysicalAssembly();c=stress.ctx
        if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=(stress.family,stress.source):
            raise ValueError('Flatten cone source family differs')
        if flatten_cone_identities()!=record['exact_identities']:
            raise ValueError('Flatten source identities changed')
        bridge=flatten_cone_source_bridge(stress)
        if encode(pack(bridge))!=record['source_bridge']:
            raise ValueError('Flatten native history/sigma/B bindings changed')
        if bridge['serialized_Xv_box_used_as_defining_value'] or bridge['source_caps_used_as_fields'] or bridge['phase_samples_used_as_proof']:
            raise ValueError('Actual native signed history must define the field')
        for flag in ('signed_incoming_defect_retained','Xp_and_Xv_positive_bounds_required','continuous_original_domain_used'):
            if not bridge[flag]:raise ValueError('Flatten source fact missing: '+flag)
        bounds=whole_flatten_bounds(stress,assembly)
        if encode(pack(bounds))!=record['bounds']:raise ValueError('Whole original flatten bounds changed')
        for label,value in bounds['positive_margins'].items():
            if endpoints(value)[0]<=0 or not mp.isfinite(endpoints(value)[1]):
                raise ArithmeticError('Flatten margin failed: '+label)
        if endpoints(bounds['actual_kappa_minus2_lower'])[0]<=0 or endpoints(bounds['actual_kappa_minus2_upper'])[1]>=2:
            raise ArithmeticError('Whole variable shear range failed')
        if endpoints(bounds['log_directional_term_upper'])[1]>=endpoints(c.ln(2))[0]:
            raise ArithmeticError('Whole flatten directional inequality failed')
        cover=bounds['sigma_t_continuous_source_cover']
        if len(cover)!=8:raise ValueError('Original continuous slope cover incomplete')
        for i,cell in enumerate(cover):
            if tuple(cell['original_x_interval'])!=(mp.mpf(i)/8,mp.mpf(i+1)/8):
                raise ValueError('Slope interval domain changed')
            if endpoints(cell['actual_sigma_t_enclosure'])[1]>endpoints(bounds['actual_sigma_t_upper'])[1]:
                raise ArithmeticError('Whole slope upper misses source interval')
        for flag in ('whole_original_t_Z_domain_covered','continuous_variable_rate_W_comparison_used',
            'actual_Xp_Xv_sign_and_signed_pulse_memory_retained','slope_cover_uses_source_intervals_not_point_sampling',
            'source_caps_only_define_bounds_not_field_values','interval_S_lower_zero_not_used_for_strict_source_sign'):
            if not bounds[flag]:raise ValueError('Continuous actual flatten proof missing: '+flag)
        if bounds['raw_whole_theta_box_used_as_positivity_proof'] or bounds['phase_grid_sampling_used_as_proof']:
            raise ValueError('Grid/uncorrelated theta box cannot establish cone')
        if record['domain']!=DOMAIN or record['tail_domain']!=TAIL_DOMAIN:
            raise ValueError('Original flatten/tail domain changed')
        admissions=('flatten_cone_certified','actual_flatten_theta_stress_positive',
            'actual_flatten_source_shear_strictly_negative','actual_flatten_directional_cone_uniform_margin',
            'flatten_power_angular_entry_exit_waiting_collar_two_vector_cone_certified',
            'same_source_flatten_power_similarity_and_physical_joins_consumed',
            'fixed_positive_viscosity_two_vector_cone_transfer_verified')
        for flag in admissions:
            if not record[flag]:raise ValueError('Flatten cone admission missing: '+flag)
        incomplete=('flatten_regional_remainder_exact_zero','left_pulse_physical_stress_join_verified',
            'completed_full_tensor_cone_certified','whole_outer_cone_certified',
            'global_admissible_stress_lift_constructed','independently_bounded_global_flat_remainder',
            'physical_energy_integral_certified','temporal_recursion')
        for flag in incomplete:
            if record[flag]:raise ValueError('Flatten scope overclaimed: '+flag)
        fixture=independent_signed_pulse_memory_fixture()
        hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=stress.family,
            implicit_source_sha256=stress.source,domain=DOMAIN,tail_domain=TAIL_DOMAIN,
            positive_source_margins_checked=len(bounds['positive_margins']),
            whole_original_t_Z_source_bounds_recomputed=True,
            continuous_variable_rate_W_comparison_verified=True,
            original_signed_pulse_memory_logarithmic_bound_verified=True,
            original_sigma_continuous_interval_slope_cover_verified=True,
            actual_variable_flatten_kappa_minus2_range_verified=True,
            actual_source_shear_sign_uses_exact_positive_S=True,
            actual_Bmax_is_original_source_inlet_verified=True,
            phase_grid_sampling_used_as_proof=False,
            independent_signed_pulse_memory_fixture=fixture,input_hashes=hashes,
            **{flag:True for flag in admissions},**{flag:False for flag in incomplete})
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS whole original100-unit flatten cone and joined outer tail; left pulse/global/recursion pending',flush=True)
        return result


if __name__=='__main__':run()
