"""Original O2 phase-held Z primitives, with actual source factors/errors.

The implicit inverse is differentiated at fixed true radius phase. Positive
rho, s and inverse-h survive the large-u limit as source factors. Saved
accepted defining-quadrature coefficients are reused, never field caps.
"""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_O2_conditioned_primitives as base

HERE,PREFIX,sha=base.HERE,base.PREFIX,base.sha
NAME=PREFIX+'current_original_O2_conditioned_slow_Z.json'
RECEIPT=PREFIX+'current_original_O2_conditioned_slow_Z_check.json'
GATE='original_O2_true_radius_phase_held_Z_primitive_enclosures_connected'
ep=base.point.endpoints
square=base.current.square


def read_scalar(c,record):
    if isinstance(record,dict) and 'exact_mpf_tuple' in record:
        return c.make_mpf(tuple(record['exact_mpf_tuple']))
    return c.mpf(record)


def restore_point(c,saved):
    result=dict(saved);rows={}
    for key,pair in saved['inputs'].items():
        values=[]
        for row in pair:
            terms=[]
            for term in row['terms']:
                late=term['late_pressure_error']['log_upper']
                terms.append(base.point.FactoredPointTerm(tuple(term['original_R_Pstar_delta_L_powers']),
                    read_scalar(c,term['approximate_source_point_coefficient']),
                    read_scalar(c,term['directed_finite_coefficient_absolute_error_upper']),
                    None if late is None else read_scalar(c,late)))
            values.append(base.point.FactoredPointRow(tuple(terms),row['original_y'],row['original_Z']))
        rows[key]=tuple(values)
    result['inputs']=rows
    return result


def exact_derivative_identities():
    r,s,q,chi,psi,rz=sy.symbols('r s q chi psi rz',nonzero=True,real=True)
    T1=q*sy.sqrt(s)/r*(chi-psi)
    H=(2-3*s)*chi+s*psi+2*r*sy.sin(chi)
    T2=q*q/r**2*H
    derivative=lambda F:(sy.diff(F,r)-2*r*sy.diff(F,s)+2*sy.sin(chi)/s*sy.diff(F,chi))*rz
    want1=q*sy.sqrt(s)/r*(2*rz*sy.sin(chi)/s-rz/(r*s)*(chi-psi))
    want2=q*q/r**2*(-2*r*rz*(psi-3*chi)+2*rz*sy.sin(chi)
        +2*rz*sy.sin(chi)/s*(2-3*s+2*r*sy.cos(chi))-2*rz/r*H)
    assert sy.simplify((derivative(T1)-want1).subs(s,1-r*r))==0
    assert sy.simplify(derivative(T2)-want2)==0
    for k in range(3,9):
        G=r**k+s*(k-1)*r**(k-2)/2
        expected=r**(k-1)+(k-1)*(k-2)*s*r**(k-3)/2
        assert sy.expand((sy.diff(G,r)-2*r*sy.diff(G,s)-expected).subs(s,1-r*r))==0
    return dict(passed=True,original_transformed_T1_T2_total_r_derivatives=True,
        original_small_r_Fourier_derivative_coefficients=True,
        phase_held_implicit_identity='psi_Z=-T2_Z/(1+t²)',
        original_O2_a_Z_b_Z_q_Z_exact_zero=True,
        primitive_units='A_Z dimensionless; B_Z_slow is derivative of B_over_Pstar')


def slow_values(kernel,roots,coordinate,chart):
    """Enclose derivatives of the source function, not of an interval selector."""
    c=kernel.c;scalar=kernel.scalar;q=kernel.q;a=kernel.a
    if any(not roots[key][(0,1)].zero for key in ('a','b')) or not kernel.t0.zero:
        raise ValueError('This layer requires the exact original O2 a_Z=b_Z=t0=0 identities')
    x=c.mpf(coordinate);lo,hi=ep(x)
    p2Z=roots['p2'][(0,1)];EZ=roots['E'][(0,1)]
    uZ=(p2Z*q).positive_divide(kernel.dstar,kernel.dstar.scale.evaluate())
    if kernel.flat or lo==hi and lo in (0,mp.mpf('.5'),1):
        return dict(psi_Z=scalar(0),A_Z_slow=scalar(0),B_Z_slow=scalar(0)),dict(
            branch='exact_O2_symmetry_or_flat',u_Z=uZ.record(),p2_Z_retained=True,E_Z_term_retained=True)
    psi_fraction,chi_fraction=kernel.angles(x,chart);psi=2*c.pi*psi_fraction
    if kernel.geometry=='small_r_series':
        r=kernel.r;s=kernel.s;hinv=kernel.hinv
        if kernel.u.zero:
            T1=q*(2*c.sin(psi));T1Z=q*uZ*c.sin(2*psi)
            T2Z=square(q)*uZ*(4*(c.sin(psi)+c.sin(3*psi)/3))
            t=q*(2*c.cos(psi));branch='exact_midplane_nonzero_p2_Z'
        else:
            R=max(abs(v) for v in ep(r));M=48
            if R>mp.mpf('.25'):raise ArithmeticError('Regular Fourier branch range exceeded')
            W1=c.mpf(0);W1r=c.mpf(0);S2=c.mpf(0)
            for k in range(1,M+1):
                sine=c.sin(k*psi)/k
                W1+=r**(k-1)*sine
                if k>=2:W1r+=(k-1)*r**(k-2)*sine
                coeff=r**(k-1)
                if k>=3:coeff+=(k-1)*(k-2)*s*r**(k-3)/2
                S2+=coeff*sine
            if R:
                rc=c.mpf(R)
                tails=(rc**M/((M+1)*(1-rc)),rc**(M-1)/(1-rc),
                    rc**M/((M+1)*(1-rc))+rc**(M-2)/2*((M+1)/(1-rc)+rc/(1-rc)**2))
                W1,W1r,S2=(v+c.mpf((-ep(tail)[1],ep(tail)[1])) for v,tail in zip((W1,W1r,S2),tails))
            T1=q*hinv*(2*W1)
            T1Z=q*uZ*(2*(-r*s*W1+s*s*W1r))
            rZ=uZ*hinv*s;T2Z=square(q)*rZ*(4*S2)
            D=1-2*r*c.cos(psi)+r*r
            if ep(D)[0]<=0:raise ArithmeticError('Positive small-r denominator lost')
            t=q*hinv*(2*(c.cos(psi)-r)/D);branch='regular_small_r_Fourier'
    elif kernel.geometry=='signed_Mobius':
        r=kernel.r;rho=kernel.rho;s=kernel.s_source;hinv=kernel.hinv
        chi=2*c.pi*chi_fraction;K=uZ*hinv;rZ=K*s;SZ=-rZ*(2*r)
        difference=chi-psi
        if chart=='E':
            sinchi=scalar(c.sin(chi));chiZ=K*(2*c.sin(chi))
            if kernel.sign>0:
                NE=scalar(2*c.cos(chi/2)**2)-rho
                lam=scalar(4*c.cos(chi/2)**2)-rho*(2*c.cos(chi))-s*3
            else:
                NE=rho-scalar(2*c.sin(chi/2)**2)
                lam=scalar(4*c.sin(chi/2)**2)+rho*(2*c.cos(chi))-s*3
            t=(q*NE*2).positive_divide(hinv,hinv.scale.evaluate()+c.ln(c.mpf(ep(hinv.coefficient)[0])))
        else:
            half=c.sin(psi/2)**2 if kernel.sign>0 else c.cos(psi/2)**2
            D=square(rho)+scalar(4*abs(r)*half)
            lower=2*(rho.scale.evaluate()+c.ln(c.mpf(ep(rho.coefficient)[0])))
            sinchi=(s*c.sin(psi)).positive_divide(D,lower)
            chiZ=(rZ*(2*c.sin(psi))).positive_divide(D,lower)
            numerator=rho-scalar(2*half) if kernel.sign>0 else scalar(2*half)-rho
            t=(q*hinv*numerator*2).positive_divide(D,lower)
            lam=(s*(2*(1-r*c.cos(psi)))).positive_divide(D,lower)-s*3
        H=scalar(2*chi)+s*(psi-3*chi)+sinchi*(2*r)
        T1=q*hinv*(difference/r)
        T1Z=q*hinv*(1/r)*(chiZ-K*(difference/r))
        T2Z=square(q)*(1/r**2)*(SZ*(psi-3*chi)+rZ*sinchi*2+chiZ*lam-rZ*H*(2/r))
        branch='conditioned_signed_'+chart
    else:raise ValueError('Refine signed original source inputs before differentiating')
    psiZ=(-T2Z).positive_divide(scalar(1)+square(t),0)
    AZ=-a*psiZ*(1/(4*c.pi))
    BZ=-a*(EZ*T1+kernel.E*(T1Z+t*psiZ))*(1/(4*c.pi))
    return dict(psi_Z=psiZ,A_Z_slow=AZ,B_Z_slow=BZ),dict(branch=branch,
        u_Z=uZ.record(),T1_fixed_angle=T1.record(),T1_Z_fixed_angle=T1Z.record(),
        T2_Z_fixed_angle=T2Z.record(),direction=t.record(),p2_Z_retained=True,E_Z_term_retained=True,
        positive_implicit_denominator_lower=1,positive_rho_s_hinv_not_rounded_to_zero=True,
        differentiated_actual_source_functions_not_interval_midpoints=True)


class OriginalO2ConditionedSlowZ:
    mode='original_O2_true_radius_phase_held_Z_primitive_enclosures'
    def __init__(self):
        receipt=json.loads((HERE/base.RECEIPT).read_bytes())
        if not receipt.get('all_passed') or not receipt.get(base.GATE):raise ValueError('Accepted actual O2 C0 receipt required')
        self.owner=base.OriginalO2ConditionedPrimitives();self.family=self.owner.family
        if receipt['source_family']!=self.family:raise ValueError('Accepted O2 C0 family differs')
        self.hashes=dict(self.owner.hashes)
        for name,digest in {**receipt['input_hashes'],base.RECEIPT:sha(base.RECEIPT)}.items():
            if sha(name)!=digest:raise ValueError('Accepted source cache dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('C0/source dependencies disagree')
            self.hashes[name]=digest
        self.saved=json.loads((HERE/base.NAME).read_bytes())
        if self.saved['effective_numerical_precisions']!=self.owner.precisions:
            raise ValueError('Accepted cached quadrature precision and live error frame differ')
        for sample in self.saved['actual_original_O2_conditioned_point_queries']:
            saved=sample['actual_original_O2_factored_inputs']
            if saved['source_family']!=self.family or saved['mode']!=self.owner.inputs.mode:
                raise ValueError('Only actual defining-quadrature point/error caches are permitted')
            y=base.point.source.exact_coordinate(saved['original_y_exact']);Z=base.point.pressure.exact_Z(saved['original_Z_exact'])
            self.owner.inputs.point_cache[(y,Z)]=restore_point(self.owner.inputs.ctx,saved)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)

    def evaluate(self,*,y,Z,N,bits=80):
        c=self.owner.ctx;phase=self.owner.radius.evaluate(y=y,N=N)
        reused=(base.point.source.exact_coordinate(y),base.point.pressure.exact_Z(Z)) in self.owner.inputs.point_cache
        with mp.workdps(c.dps+40):
            query=self.owner.query(y=y,Z=Z);kernel=query['kernel'];rows=[]
            for box in phase['true_original_phase_directed_boxes']:
                inverse=kernel.evaluate(c.mpf([box['lower'],box['upper']]),bits=bits)
                if inverse['status']!='enclosed':raise ArithmeticError('Actual source inverse needs refinement')
                inverse.update(free_phase_parameter_not_spatial_phase=False,
                    original_common_N_and_radius_phase_bound=True,
                    actual_O2_defining_point_inputs_and_errors_consumed=True,
                    source_caps_or_midpoints_used_as_field_values=False)
                selected=inverse['selected_inverse']
                values,proof=slow_values(kernel,query['roots'],selected['coordinate_interval'],selected['chart'])
                rows.append(dict(C0=inverse,slow_Z={k:v.record() for k,v in values.items()},derivative_contract=proof))
        return dict(source_family=self.family,mode=self.mode,original_y_exact=str(base.point.source.exact_coordinate(y)),
            original_Z_exact=str(base.point.pressure.exact_Z(Z)),explicit_candidate_N=N,
            actual_original_radius_phase=phase,source_factor_basis=query['basis_contract'],
            effective_numerical_precisions=self.owner.precisions,actual_phase_held_Z_enclosures=rows,
            actual_original_O2_factored_inputs=base.point.record_query(query['point']),
            source_geometry=kernel.geometry_record(),numerical_arithmetic_ledger=query['ledger'],
            slow_Z_primitives_installed_on_this_O2_point=True,true_radius_phase_Z_exact_zero=True,
            cached_actual_defining_quadrature_bits_and_errors_reused=reused,source_caps_or_midpoints_used_as_field_values=False,
            actual_changed_five_moment_integral_evaluated=False,numerical_original_source_point_or_integral_oracle_installed=False,
            actual_five_controls_installed=False,current_whole_N_selected=False)


def run():
    began=time.monotonic();owner=OriginalO2ConditionedSlowZ();records=[]
    for y,Z,N in (('.53','-.37',7),('.53','0',7),('.53','.37',7),('1','.37',257)):
        row=owner.evaluate(y=y,Z=Z,N=N);records.append(row)
        print('True-radius O2 slow-Z:',y,Z,N,row['actual_phase_held_Z_enclosures'][0]['derivative_contract']['branch'],flush=True)
    # Choose one explicit diagnostic N whose genuine radius phase lies in
    # its source-signed narrow peak. This is coverage, not common-N admission.
    c=owner.owner.ctx
    with mp.workdps(c.dps+40):
        kernel=owner.owner.query(y='.53',Z='.37')['kernel']
        cap=ep(kernel.nq)[0]/2
        for N in range(1,65):
            phase=owner.owner.radius.evaluate(y='.53',N=N)
            boxes=phase['true_original_phase_directed_boxes']
            if len(boxes)!=1:continue
            lo,hi=boxes[0]['lower'],boxes[0]['upper']
            peak=(0<lo<=hi<cap or 1-cap<lo<=hi<1) if kernel.sign>0 else (mp.mpf('.5')-cap<lo<=hi<mp.mpf('.5')+cap)
            if not peak:continue
            row=owner.evaluate(y='.53',Z='.37',N=N);records.append(row)
            assert row['actual_phase_held_Z_enclosures'][0]['derivative_contract']['branch']=='conditioned_signed_E'
            print('True-radius O2 slow-Z peak:',N,'conditioned_signed_E',flush=True)
            break
        else:raise ArithmeticError('Explicit diagnostic N did not cover the native narrow peak')
    report=dict(**{GATE:True},source_family=owner.family,mode=owner.mode,
        actual_original_O2_slow_Z_point_queries=records,exact_original_derivative_identities=exact_derivative_identities(),
        slow_Z_primitives_installed_on_actual_O2_points=True,source_auxiliary_scales_and_full_point_Z_errors_preserved=True,
        no_original_defining_quadratures_or_ancestor_constructors_reexecuted=True,
        actual_changed_five_moment_integral_evaluated=False,numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(base.point.source.inertial.profiles.loop.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,
        scope='Local original O2 true-radius phase-held Z primitive enclosures. Cached accepted source quadrature coefficients/errors retain native factors; not densities, integrals, all-chart oracle, controls, common N or recursion.')
    (HERE/NAME).write_text(json.dumps(base.encoded(report),indent=2)+'\n',encoding='utf8')
    return report


if __name__=='__main__':run()
