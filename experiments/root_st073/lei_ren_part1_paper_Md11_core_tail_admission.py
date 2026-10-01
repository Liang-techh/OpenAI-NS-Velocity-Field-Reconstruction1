"""Admit the new analytic pressure into unchanged universal core majorants.

The old fixed-point solution and finite coefficients are not transferred.
Only upper bounds for the same fixed-axis contraction map are reused, after
checking every pressure-dependent norm is dominated by its old majorant.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_analytic_radial_tail import tail_factor
from lei_ren_part1_paper_candidate_general_center_factory import PARAMETERS
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def run():
    names=('lei_ren_part1_paper_Md11_pressure_datum.json',
        'lei_ren_part1_paper_complex_pressure_bound.json',
        'lei_ren_part1_paper_analytic_radial_tail.json',
        'lei_ren_part1_paper_linear_resolvent_bound.json',
        'lei_ren_part1_paper_commuting_resolvent_bound.json',
        'lei_ren_part1_paper_nonlinear_candidate_Lambda120.json')
    loaded={};hashes={}
    def visit(name):
        if name in loaded:return loaded[name]
        raw=(HERE/name).read_bytes();hashes[name]=hashlib.sha256(raw).hexdigest()
        record=json.loads(raw) if name.endswith('.json') else {}
        loaded[name]=record
        for dep,digest in record.get('input_hashes',{}).items():
            visit(dep)
            if hashes[dep]!=digest:raise ValueError('analytic dependency changed:'+dep)
        return record
    new,old,tube,linear,commuting,candidate=[visit(n) for n in names]
    if not candidate['contraction_proved'] or not candidate['all_paper_8_50_terms_included']:
        raise ValueError('complete universal contraction majorant required')
    c=MPIntervalContext();c.dps=160
    with mp.workdps(220):
        get=lambda r,n:read_interval(c,r[n])
        eta=get(tube,'complex_tube_radius');h=get(tube,'Xh_parameter')
        qlower=get(tube,'pressure_q_modulus_lower')
        m2=get(new,'fixed_beta2_mass_upper_normalized');m0=get(new,'fixed_beta0_mass_upper_normalized')
        flat=get(new,'flatten_true_mass_upper_normalized');scale=get(new,'Pstar_squared')
        if len(new['stages'])!=14 or any(endpoints(m)[0]<0 for m in (m2,m0,flat)):
            raise ValueError('positive fourteen-stage analytic pressure required')
        if (not new['all_14_true_pressure_stages_included'] or new['pressure_units']!='physical_P'
            or new['finite_pressure_quadrature_approximation_used'] or new['pressure_parameter_order_truncated']
            or any(endpoints(get(row,'mass'))[0]<0 for row in new['stages'].values())):
            raise ValueError('complete untruncated true positive physical pressure required')
        if PARAMETERS!=dict(j='1e-14',Lambda='1e120',logC='5e151',delta='1e-200',sigma_denominator='500'):
            raise ValueError('fixed-axis data changed; no monotonicity transfer')
        if new['new_schedule_inputs']['logPstar'] not in ('14','14.0'):
            raise ValueError('this pressure admission assumes unchanged logPstar=14')
        sl,sh=endpoints(scale);el,eh=endpoints(c.exp(28))
        if not el<=sl<=sh<=eh:
            raise ValueError('physical pressure scale must be exp(28)')
        if endpoints(get(candidate,'Lambda'))!=endpoints(c.mpf('1e120')):
            raise ValueError('candidate Lambda changed')
        normalized=(m2+flat)/qlower**2+m0
        pressure=scale*normalized;derivative=4*(1+eta)*pressure/qlower
        a=1+eta/2;dt=c.mpf(PARAMETERS['delta']);j=c.mpf(PARAMETERS['j'])
        U=4*a+j;H=get(tube,'H0_modulus_upper');L=get(tube,'L_modulus_lower')
        g=(1+dt)/2*(1+2*a*U)*U+4*H+2*(1+dt)*a*pressure+(1+a*a)*derivative
        g_over_L=g/L;ratio=h/(eta/2)
        weight=c.mpf([1,max(mp.mpf(1),endpoints(4*ratio)[1])])
        gnorm=g_over_L*weight;psi=40*gnorm
        fresh=dict(normalized_complex_pressure_modulus_upper=normalized,
            physical_P0_complex_modulus_upper=pressure,
            P0_first_derivative_modulus_upper_on_half_capsule=derivative,
            g_over_L_modulus_upper_on_half_capsule=g_over_L,
            g_over_L_Xh_norm_upper=gnorm,Psi_model_Xh_norm_upper=psi)
        admission={}
        for key,value in fresh.items():
            # A directed upper endpoint is itself a conservative scalar bound.
            oldbound=c.mpf(endpoints(get(old,key))[1]);newbound=c.mpf(endpoints(value)[1])
            slack=oldbound-newbound
            if endpoints(slack)[0]<=0:raise ValueError('pressure majorant not dominated:'+key)
            admission[key]=dict(new_bound=value,reused_old_upper_bound=oldbound,strict_slack=slack)
        correction=get(candidate,'scaled_map_size_upper')
        lip=get(candidate,'scaled_map_Lipschitz_upper')
        if endpoints(correction)[1]>=1 or endpoints(lip)[1]>=1:
            raise ValueError('universal contraction/self-map gates not strict')
        # The old termwise map is polynomial with nonnegative coefficients in
        # Bphi,Bpsi; external pressure enters only the Psi center bound.
        # Same axis/tube/resolvents and dominated Bpsi retain both map gates.
        phi_norm=c.mpf([0,endpoints(get(linear,'Phi_model_Xh_norm_upper'))[1]])+correction
        psi_norm=c.mpf([0,endpoints(get(old,'Psi_model_Xh_norm_upper'))[1]])+correction
        rows=[];worst=c.mpf(0)
        for total in range(4):
            for i in range(total+1):
                k=total-i
                factor=tail_factor(c,degree=124,radial_order=i,axial_order=k,radius='4.1',h=h)['tail_per_Xh_norm']
                pb=phi_norm*factor;ub=psi_norm*factor
                worst=c.mpf([0,max(endpoints(worst)[1],endpoints(pb)[1],endpoints(ub)[1])])
                rows.append(dict(scaled_radial_order=i,axial_order=k,Phi_tail=pb,Psi_tail=ub))
        if endpoints(worst)[1]>mp.mpf('1e-12'):raise ValueError('new degree124 tail gate not met')
        output=dict(new_schedule_sha256=new['new_schedule_sha256'],parameters=PARAMETERS,
            pressure_admission=admission,all_pressure_norm_gates_strict=True,
            unchanged_axis_tube_and_resolvent_data=True,
            majorant_transfer='external pressure enters Eq.8.50 size/Lipschitz polynomials only through Bpsi; all powers and coefficients nonnegative',
            source_audit='nonlinear_map_majorant.py lines38-92; complex_pressure_bound.py lines34-81',
            scaled_map_size_upper=correction,scaled_map_Lipschitz_upper=lip,
            analytic_fixed_point_exists_for_new_fixed_datum=True,
            Xh_parameter=h,analytic_Phi_Xh_norm_upper=phi_norm,analytic_Psi_Xh_norm_upper=psi_norm,
            radial_degree=124,maximum_normalized_mixed_C3_tail=worst,tail_rows=rows,target_met=True,
            real_scaled_radial_domain=['0','4.1'],analytic_axis_domain=['-1','1'],
            old_fixed_point_solution_transferred=False,old_finite_coefficients_reused=False,
            original_parameter_errors_enclosed=False,unknown_outer_smallness_constants_verified=False,
            waiting_root_certified=False,full_NS_residual_bound=False,temporal_recursion=False)
        for name in (Path(__file__).name,'lei_ren_part1_paper_complex_pressure_bound.py',
            'lei_ren_part1_paper_nonlinear_map_majorant.py','lei_ren_part1_paper_candidate_core_tail_budget.py',
            'lei_ren_part1_paper_analytic_radial_tail.py','lei_ren_part1_paper_candidate_general_center_factory.py'):
            hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        output['input_hashes']=hashes
        Path(__file__).with_suffix('.json').write_text(json.dumps(encode(output),indent=2)+'\n',encoding='utf-8')
        print('Md1.1 pressure admitted to universal contraction; degree124 mixedC3 tail',mp.nstr(endpoints(worst)[1],16),flush=True)
        return output

if __name__=='__main__':run()
