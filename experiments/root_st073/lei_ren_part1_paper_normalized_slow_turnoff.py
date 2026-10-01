"""Uniform O.2 bounds from the exact normalized propagation equations.

This is a conditional parameter-family certificate, not a generated core or
fourteen-stage numerical pressure datum. The pressure envelope is for the
same complete PREHEAT backward integral, never a fitted pressure substitute.
No enormous physical amplitude or Decimal stage coordinate is constructed.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_outer_slope_field import evaluate_transition
from lei_ren_part1_paper_interval_long_reshape_field import sigma_value_derivative
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent


def upper(c,v):return c.mpf(endpoints(v)[1])


def lower(c,v):return c.mpf(endpoints(v)[0])


def sigma_derivative_upper(c,left,right):
    """Enclose a full closed phase cell, including the flat endpoints."""
    if left==0 or right==1:
        a=right if left==0 else 1-left
        if a>c.mpf('.25'):return c.mpf(32)
        # exp(-1/q^2)/q^3 increases on (0,.25]. Symmetry handles q=1.
        return upper(c,2*c.exp(-1/a**2+1/(1-a)**2)*(1/a**3+1/(1-a)**3))
    return upper(c,sigma_value_derivative(c,c.mpf([endpoints(left)[0],endpoints(right)[1]]))[1])


def certificate(Md,center=('.49','.51'),inlet_cells=512,phase_cells=512):
    c=MPIntervalContext();c.dps=120
    with mp.workdps(160):
        md=c.mpf(Md)
        if endpoints(md)[0]<2:raise ValueError('uniform amplitude bound requires Md>=2')
        z=IntervalTaylor(c,[c.mpf(list(center)),c.mpf(1)])
        if endpoints(z[0])[0]<=0 or endpoints(z[0])[1]>=1:
            raise ValueError('this certificate requires a positive compact axial slab inside (0,1)')
        zz=z[0];za=upper(c,zz);d=1-zz*zz
        delta=c.mpf([0,'1e-200']);epsilon=c.mpf([0,'1e-202']);L=1-delta*zz*zz
        # Only the exact unit angular shape is used to obtain Q at y=1.
        # I_theta and Q do not depend on pressure or on angular amplitude.
        A=(1+z*z).reciprocal()
        inlet=evaluate_transition(c,z,delta,c.mpf(1),A,z*0,1,inlet_cells)
        Q1=inlet['stress']['I_theta']*L/(inlet['F'][0]*inlet['R'])
        H0=(c.mpf(5)/12+inlet['dimensionless_increment_integrals'][2]/2)*c.exp(c.mpf('-.6'))
        # If log(Pstar)>=exp(Md)+11, then on 1<=y<=exp(Md),
        # 1/u^2 <= 4 exp(-exp(Md)-22.6) <= this Md>=2 bound.
        inv_u2=upper(c,4*c.exp(-c.exp(c.mpf(2))-c.mpf('22.6')))
        inflation=1/(1-epsilon)**2
        # Complete preheat angular profile is bounded by the continuation
        # u(y)*exp(-(s-y)/2)/(1-epsilon); its logarithmic Z derivative is
        # bounded by 2|Z|/(1+Z^2). Integrate both bounds to infinity.
        Pi_abs=upper(c,za*(1+delta+2*d/(1+zz*zz))*inflation)
        swirl_J1=upper(c,za*H0*(2*delta+4*d/(1+zz*zz)))
        axial_J1=upper(c,4*za*(5+(8+4*delta)*za*za)*inv_u2)
        J1_abs=upper(c,Pi_abs+swirl_J1+axial_J1)
        # Equation 6.24, |W|<=3, 0<=f<=4 and |f'|<=128/(Md*y).
        # On this positive slab Pi<=0 for the SAME preheat integral.
        # Thus |Z+Pi|<=max(Z_upper, Pi_abs-Z_lower), not Z_upper+Pi_abs.
        pressure_force_abs=c.mpf(max(endpoints(za)[1],endpoints(Pi_abs-lower(c,zz))[1]))
        Jprime_abs=upper(c,pressure_force_abs+za*inv_u2*(384/md+4+16*(1+delta)))
        J_over_y_abs=c.mpf(max(endpoints(J1_abs)[1],endpoints(Jprime_abs)[1]))
        forcing_lower=lower(c,(1-delta)*zz*zz/(1+zz*zz)-delta/2)
        Q_lower=c.mpf(min(endpoints(Q1)[0],endpoints(forcing_lower)[0]))
        x_lower=lower(c,Q_lower/upper(c,L))
        # Paper b_p=2f'Z/u is signed. Here h=-b_p, as in F18.
        # |h*N_z|=|2f'ZJ/L|; the tiny common u cancels exactly.
        product_abs=upper(c,256*za*J_over_y_abs/(md*lower(c,L)))
        b_squared_upper=upper(c,(256*za/md)**2*inv_u2)
        # Rref>=exp(10) gives R>=exp(11) over the entire cutoff and buffer.
        inverse_R_upper=upper(c,c.exp(-11))
        tau_lower=lower(c,x_lower-2*inverse_R_upper)
        finite_direction_lower=lower(c,2*x_lower-product_abs-(4+b_squared_upper)*inverse_R_upper)
        finite_bracket_lower=lower(c,(8-b_squared_upper**2/2)*x_lower
            -(8+2*b_squared_upper)*product_abs-(b_squared_upper+4)**2*inverse_R_upper)
        coarse_bracket_lower=finite_bracket_lower
        phase_rows=[]
        if phase_cells:
            if not isinstance(phase_cells,int) or phase_cells<4:raise ValueError('phase_cells>=4 required')
            for i in range(phase_cells):
                left=c.mpf(i)/phase_cells;right=c.mpf(i+1)/phase_cells
                ds=sigma_derivative_upper(c,left,right)
                q=c.mpf([endpoints(left)[0],endpoints(right)[1]])
                inv_y=c.exp(-md*q)
                jy=upper(c,Jprime_abs+(J1_abs-Jprime_abs)*inv_y)
                prod=upper(c,8*za*ds*jy/(md*lower(c,L)))
                b2=upper(c,(8*za*ds/md*inv_y)**2*inv_u2)
                # g(B) in Q'+Q=g(B) is affine increasing in B. Since
                # B decreases, the positive convolution gives
                # Q(y)>=min(Q1,g(B(y))) without computing exp(-exp(Md*q)).
                cutoff_lower=lower(c,1-sigma_value_derivative(c,right)[0])
                forcing_cell=lower(c,(1-delta)*zz*zz/(1+zz*zz)-delta/2
                    +cutoff_lower*(8*d*zz*zz/(1+zz*zz)+4*delta*zz*zz))
                Q_cell=c.mpf(min(endpoints(Q1)[0],endpoints(forcing_cell)[0]))
                x_cell=lower(c,Q_cell/upper(c,L))
                direction=lower(c,2*x_cell-prod-(4+b2)*inverse_R_upper)
                bracket=lower(c,(8-b2**2/2)*x_cell-(8+2*b2)*prod-(b2+4)**2*inverse_R_upper)
                phase_rows.append(dict(index=i,phase=q,sigma_derivative_upper=ds,
                    J_over_y_absolute_upper=jy,absolute_b_Nz_upper=prod,b_squared_upper=b2,
                    Q_lower=Q_cell,Ntheta_lower=x_cell,
                    finite_negative_dot_per_R_lower=direction,finite_strong_bracket_lower=bracket))
            finite_bracket_lower=c.mpf(min(endpoints(r['finite_strong_bracket_lower'])[0] for r in phase_rows))
            finite_direction_lower=c.mpf(min(endpoints(r['finite_negative_dot_per_R_lower'])[0] for r in phase_rows))
        passed=all(endpoints(v)[0]>0 for v in
            (Q_lower,tau_lower,finite_direction_lower,finite_bracket_lower))
        result=dict(Md=Md,axial_slab=list(center),inlet_cells=inlet_cells,
            Q_at_slope_endpoint=Q1,swirl_energy_inlet_H0=H0,
            inverse_u_squared_upper=inv_u2,pressure_Pi_absolute_upper=Pi_abs,
            J_inlet_absolute_upper=J1_abs,J_radial_derivative_absolute_upper=Jprime_abs,
            J_over_y_absolute_upper=J_over_y_abs,Q_uniform_lower=Q_lower,
            Ntheta_uniform_lower=x_lower,absolute_b_Nz_upper=product_abs,
            b_squared_upper=b_squared_upper,finite_tau_lower=tau_lower,
            finite_negative_dot_per_R_lower=finite_direction_lower,
            finite_strong_bracket_lower=finite_bracket_lower,
            coarse_finite_strong_bracket_lower=coarse_bracket_lower,
            full_phase_interval_cells=phase_rows,phase_cell_count=phase_cells,
            whole_cutoff_phase_conditionally_certified=passed,
            cutoff_phase_domain='q in [0,1], equivalently y in [1,exp(Md)]',
            zero_axial_buffer_domain='y in [exp(Md),exp(Md)+11]',
            buffer_inverse_u_squared_bound_used=False,
            retained_moment_zero_Uz_buffer_conditionally_certified=endpoints(tau_lower)[0]>0,
            whole_O2_slow_turnoff_conditionally_certified=passed and endpoints(tau_lower)[0]>0,
            endpoint_branch='b=0 uses weak kappa=2 inequality; interior b!=0 uses strong margin',
            cone_scope='relaxed condition (3.23), strong branch where kappa>2, weak branch at b=0',
            admissible_stress_lift_constructed=False,
            parameter_family=dict(logPstar='>= exp(Md)+11',delta='min(1e-200,exp(-4logPstar-30))',
                c_mu='.001',c_delta='.001',c_epsilon='.01',logRref='>=10',waiting_length='any nonnegative value'),
            pressure_scope='same complete reference-plus-outer preheat backward integral',
            pressure_envelope_inflation=inflation,
            actual_epsilon_relation='epsilon=c_epsilon*delta=.01*delta',
            preheat_pressure_sign_lemma_used='P<=0 and P_Z>=0 on positive Z from 0<=vartheta<=1',
            pressure_envelope_added_or_fitted=False,
            explicit_conditions=['complete angular ansatz of Section 6.1 with H replaced by 1',
                'inlet exact reference cumulative moments after functional repair',
                'global preheat pressure equals its complete backward integral',
                'parameter hierarchy and nonnegative waiting length'],
            new_fourteen_stage_datum_generated=False,new_core_generated=False,
            existing_Md11_core_reused=False,assembled_background_admissibility_certified=False,
            whole_axis_certified=False,unknown_paper_constants_verified=False,
            actual_heat_candidate_pressure_envelope_certified=False,temporal_recursion=False)
        print('Normalized uniform O.2 Md',Md,'certified conditionally',passed,
            'bracket lower',mp.nstr(endpoints(finite_bracket_lower)[0],10),flush=True)
        return result


def run():
    trials=[certificate(md) for md in ('16','24','32','40','48','64')]
    names=(Path(__file__).name,'lei_ren_part1_paper_interval_outer_slope_field.py',
        'lei_ren_part1_paper_interval_long_reshape_field.py','lei_ren_part1_paper_interval_repaired_reference_field.py',
        'lei_ren_part1_paper_interval_taylor.py','lei_ren_part1_paper_interval_comparison_enclosure.py')
    result=dict(trials=trials,source_equations='6.16-6.18 and 6.24-6.26',
        full_closed_phase_cells_not_point_sampling=True,
        input_hashes={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names})
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    return result


if __name__=='__main__':run()
