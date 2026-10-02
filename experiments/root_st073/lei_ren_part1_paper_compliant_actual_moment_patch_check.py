"""Independent source-tail transport and implicit axial5 patch checks."""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_actual_moment_patch import (
    CompliantActualMomentPatch,implicit_axial_jets,IntervalTaylor)
from lei_ren_part1_paper_interval_five_bump_inverse import weighted_norm,mag
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_actual_moment_patch.json'


def symbolic_source_transport():
    gap,E,invA,I1,I16,I2=s.symbols('gap E invA I1 I16 I2',real=True)
    em,eh,ek,aa,eb,ep=s.symbols('em eh ek aa eb ep',real=True)
    rz=[E+s.exp(-gap)*(em-E),s.exp(-s.Rational(8,5)*gap)*eh,
        s.Rational(5,8)*E+s.exp(-s.Rational(8,5)*gap)*(ek-s.Rational(5,8)*E),
        E**2+s.exp(-gap)*(aa-E**2),s.exp(-s.Rational(6,5)*gap)*eb,s.exp(-gap/5)*ep]
    end=[s.exp(-1)*(rz[0]+E*I1),s.exp(-s.Rational(8,5))*rz[1],
        s.exp(-s.Rational(8,5))*(rz[2]+E*I16),s.exp(-1)*(rz[3]+E**2*I2),
        s.exp(-s.Rational(6,5))*rz[4],s.exp(-s.Rational(1,5))*rz[5]]
    rates=(1,s.Rational(8,5),s.Rational(8,5),1,s.Rational(6,5),s.Rational(1,5))
    rh=[q*s.exp(-2*k) for q,k in zip(end,rates)]
    direct=[s.exp(1)*rh[0],s.exp(s.Rational(8,5))*rh[2],s.exp(s.Rational(8,5))*rh[1],
        s.exp(1)*rh[3]*invA-s.exp(s.Rational(6,5))*rh[4]/2,s.exp(s.Rational(1,5))*rh[5]/2]
    alpha=s.exp(-gap-2)
    targets=[s.exp(-2)*(1+I1)*E+alpha*(em-E),
        s.exp(-s.Rational(16,5))*(s.Rational(5,8)+I16)*E+alpha**s.Rational(8,5)*(ek-s.Rational(5,8)*E),
        alpha**s.Rational(8,5)*eh,
        invA*(s.exp(-2)*(1+I2)*E**2+alpha*(aa-E**2))-alpha**s.Rational(6,5)*eb/2,
        alpha**s.Rational(1,5)*ep/2]
    if any(s.simplify(a-b)!=0 for a,b in zip(direct,targets)):
        raise ArithmeticError('Exact actual transport differs from old dominant-plus-tail functions')
    # The coefficient of h_n in the original quadratic map is precisely
    # the SAME pointwise Jacobian for all n; all amplitude derivatives
    # occur in the known part. This checks the algorithm independently.
    t=s.symbols('t');hc=s.symbols('h0:30');ac=s.symbols('a0:6')
    f,g,p,q,r=s.symbols('f g p q r',real=True)
    hs=[sum(hc[6*i+k]*t**k for k in range(6)) for i in range(5)]
    A=sum(ac[k]*t**k for k in range(6))
    Q=[s.Integer(0),f*hs[0]*hs[2]+g*hs[1]*hs[4],s.Integer(0),
       A*(p*hs[0]**2+q*hs[1]**2)-(r*hs[2]**2+f*hs[3]**2+g*hs[4]**2)/2,
       (p*hs[2]**2+q*hs[3]**2+r*hs[4]**2)/2]
    zero_symbols=s.symbols('c0:5');A0=s.symbols('A0')
    Q0=s.Matrix([0,f*zero_symbols[0]*zero_symbols[2]+g*zero_symbols[1]*zero_symbols[4],0,
        A0*(p*zero_symbols[0]**2+q*zero_symbols[1]**2)
            -(r*zero_symbols[2]**2+f*zero_symbols[3]**2+g*zero_symbols[4]**2)/2,
        (p*zero_symbols[2]**2+q*zero_symbols[3]**2+r*zero_symbols[4]**2)/2])
    J=Q0.jacobian(zero_symbols).subs({**{zero_symbols[i]:hc[6*i] for i in range(5)},A0:ac[0]})
    count=0
    for n in range(2,6):
        variables=[hc[6*i+n] for i in range(5)];sub={hc[6*i+k]:0 for i in range(5) for k in range(n,6)}
        for i,row in enumerate(Q):
            full=s.expand(row).coeff(t,n);known=s.expand(row.subs(sub)).coeff(t,n)
            if s.expand(full-known-sum(J[i,j]*variables[j] for j in range(5)))!=0:
                raise ArithmeticError('Higher implicit coefficient missing a nonlinear/amplitude term')
            count+=1
    return dict(exact_actual_dominant_tail_transport_identities=5,
        exact_higher_implicit_Jacobian_convolution_identities=count,passed=True)


def independent_coefficient_fixture(provider):
    with mp.workdps(90):
        c=MPIntervalContext();c.dps=115;tol=mp.mpf('1e-75');z=mp.mpf('.3')
        normalization=mp.quad(lambda r:mp.exp(-1/(1-r*r)) if abs(r)<1 else 0,[-1,0,1])
        radius=mp.mpf('.025');centers=[mp.mpf('1.25'),mp.mpf('1.5'),mp.mpf('1.75')]
        def beta(r):return mp.exp(-1/(1-r*r))/(radius*normalization) if abs(r)<1 else mp.mpf(0)
        def integral(i,p,m=1):
            return radius*mp.quad(lambda r:(centers[i]+radius*r)**p*beta(r)**m,[-1,0,1])
        # Independently integrated actual physical perturbation integrands.
        w06=[integral(i,mp.mpf('.6')) for i in (0,2)]
        w05=[integral(i,mp.mpf('.5')) for i in range(3)]
        w01=[integral(i,mp.mpf('.1')) for i in range(3)]
        wm9=[integral(i,mp.mpf('-.9')) for i in range(3)]
        fg=[integral(i,mp.mpf('.5'),2) for i in (0,2)]
        gg=[integral(i,0,2) for i in (0,2)]
        ff=[integral(i,0,2) for i in range(3)]
        ffx=[integral(i,-1,2) for i in range(3)]
        hfn=[lambda zz:mp.mpf('1e-7')*mp.sin(zz+mp.mpf('.8')),
             lambda zz:mp.mpf('1.5e-7')*mp.cos(zz/2),
             lambda zz:mp.mpf('1e-7')*mp.exp(zz/5),
             lambda zz:mp.mpf('2e-7')*(1+zz**3/10),
             lambda zz:mp.mpf('-.8e-7')*mp.sin(2*zz+1)]
        amplitude=lambda zz:(1+zz*zz)**2/mp.mpf('2.4')**2
        def dfn(zz):
            a,b,f,g,h=[fn(zz) for fn in hfn];ax=(a,b);ang=(f,g,h)
            return [-(a+b),-(w06[0]*a+w06[1]*b+fg[0]*a*f+fg[1]*b*h),
                -sum(w05[i]*ang[i] for i in range(3)),
                -amplitude(zz)*sum(gg[i]*ax[i]**2 for i in range(2))
                    +sum(w01[i]*ang[i]+ff[i]*ang[i]**2/2 for i in range(3)),
                -sum(wm9[i]*ang[i]+ffx[i]*ang[i]**2/2 for i in range(3))]
        def jet(fn):return IntervalTaylor(c,[c.mpf([v-tol,v+tol])/math.factorial(k)
            for k in range(6) for v in (mp.diff(fn,z,k),)])
        d=[jet(lambda zz,i=i:dfn(zz)[i]) for i in range(5)]
        W={n:([[c.mpf(endpoints(v)) for v in row] for row in values] if n=='L' else
             [c.mpf(endpoints(v)) for v in values]) for n,values in provider.W.items()}
        result=implicit_axial_jets(c,W,d,jet(amplitude),[c.mpf('1e-5')]*2+[c.mpf('2e-5')]*3)
        count=0
        for i,row in enumerate(result['controls']):
            for k in range(6):
                expected=mp.diff(hfn[i],z,k)/math.factorial(k);lo,hi=endpoints(row[k])
                if not lo-tol*1000<=expected<=hi+tol*1000:
                    raise ArithmeticError('Independent physical implicit coefficient derivative failed')
                count+=1
        return dict(independent_physical_bump_weight_integrals=21,
            independent_smooth_implicit_coefficient_Taylor_coefficients=count,
            fixture_inverse_contraction_upper=mp.nstr(endpoints(result['initial_C1_inverse']['contraction_bound'])[1],20),
            variable_inverse_amplitude_derivatives_retained=True,
            finite_fixture_only=True,actual_source_admission=False,passed=True)


def run():
    with mp.workdps(280):
        receipt=json.loads((HERE/NAME).read_bytes());provider=CompliantActualMomentPatch(require_checked=False);c=provider.ctx
        for name,digest in receipt['input_hashes'].items():
            if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Actual connected patch source changed: '+name)
        if not (receipt['actual_five_moment_patch_connected'] and receipt['actual_five_functional_terminal_identities_connected']
                and receipt['actual_five_defect_family_sha256']==provider.family and receipt['implicit_source_sha256']==provider.source):
            raise ValueError('Actual patch source/closure scope missing')
        inverse,data=provider.coefficients([-1,1]);tail_count=0;proof_count=0
        admitted=provider.reference.defect_admission
        prior=json.loads((HERE/'lei_ren_part1_paper_compliant_reference_restore_profiles_check.json').read_bytes())
        if not (prior['all_passed'] and prior['actual_source_binding']['same_original_cutoff_complement_identity']
                and admitted['directed_fixed_cutoff_integrals_not_midpoint_fit']
                and admitted['actual_defect_functions_defined_by_source_primitives']):
            raise ValueError('Original full restoration integral functions are not source bound')
        J=provider.reference.kernels(1)
        for name,part,rate in (('mean','linear_1',1),('mixed','linear_8_over_5','1.6'),('square','square_1',1)):
            original=read_interval(c,admitted['directed_restoration_integrals'][part])*c.exp(-c.mpf(rate))
            jl,jh=endpoints(J[name]);lo,hi=endpoints(original)
            if jl<lo or jh>hi:raise ArithmeticError('Reference endpoint refinement no longer uses identical admitted integral')
        expected_kernels=dict(d1=c.exp(-2)*(1+read_interval(c,admitted['directed_restoration_integrals']['linear_1'])),
            d2=c.exp(c.mpf('-3.2'))*(c.mpf('.625')+read_interval(c,admitted['directed_restoration_integrals']['linear_8_over_5'])),
            d4=c.exp(-2)*(1+read_interval(c,admitted['directed_restoration_integrals']['square_1'])))
        for name,expected in expected_kernels.items():
            lo,hi=endpoints(expected);kl,kh=endpoints(data['actual_dominant_kernels'][name])
            if max(lo,kl)>min(hi,kh):raise ArithmeticError('Admitted source dominant kernel formula inconsistent')
        for row,bound in zip(data['actual_tail_functions'],data['actual_tail_C1_norm_upper']):
            actual=mag(c,row[0])+mag(c,row[1])
            if endpoints(actual)!=endpoints(bound) or endpoints(actual)[1]>endpoints(provider.oldtail)[0]:
                raise ArithmeticError('Actual source tails exceed admitted correlated C1 cap')
            tail_count+=1
        for proof in data['retained_source_decay_proofs']:
            value=-proof['rate']*proof['source_gap_enclosure']+c.ln(proof['input_absolute_upper'])
            if endpoints(value)[1]>endpoints(provider.reference.logcap)[0]:raise ArithmeticError('Actual inherited source decay cap failed')
            proof_count+=1
        norm=weighted_norm(c,inverse['preconditioned_error'],provider.scales)
        if endpoints(norm)[1]>=1 or endpoints(norm)!=endpoints(inverse['higher_inverse_contraction']):
            raise ArithmeticError('Actual common higher-jet inverse not contractive')
        for proof in inverse['higher_order_proofs']:
            bnorm=max((mag(c,v)/provider.scales[i] for i,v in enumerate(proof['preconditioned_rhs'])),key=lambda v:endpoints(v)[1])
            if endpoints(bnorm/(1-norm))[1]>endpoints(proof['weighted_radius'])[1]:
                raise ArithmeticError('Higher implicit derivative radius insufficient')
        coefficient_count=0
        for live,stored in zip(inverse['controls'],receipt['actual_implicit_coefficient_inverse_axial5']['controls']):
            for a,b in zip(live.coefficients,stored['coefficients']):
                if endpoints(a)!=endpoints(read_interval(c,b)):raise ArithmeticError('Actual coefficient data changed')
                coefficient_count+=1
        profile_count=0
        packets=receipt['actual_patch_samples']+[receipt['whole_actual_patch'],receipt['actual_terminal_Rh']]
        groups=('Utheta_over_Pstar_axial5','Uz_axial5','Utheta_y_over_Pstar_axial5','Uz_y_axial5',
            'Ur_over_sqrt_R_over_2_axial4','Mz_over_R_axial5','Mtheta_over_sqrt2_Rm_1p5_Am_axial5',
            'Mtheta_z_over_sqrt2_Rm_1p5_Am_axial5','Mztheta_minus_8ZMz_plus_16Z2R_over_RmAm2_axial5',
            'Mp_over_Am2_axial5','P_over_Pstar2_axial5','original_P0_axial5')
        for packet in packets:
            for name in groups:
                row=packet[name]
                if len(row)!=(5 if name=='Ur_over_sqrt_R_over_2_axial4' else 6):raise ValueError('Actual patch axial order lost')
                for value in row:
                    lo,hi=endpoints(read_interval(c,value))
                    if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):raise ArithmeticError('Nonfinite patch field')
                    profile_count+=1
            if endpoints(read_interval(c,packet['Utheta_over_Pstar_axial5'][0]))[0]<=0:raise ArithmeticError('Actual patch positive swirl lost')
            if packet['radial_mixed4_certified'] or packet['full_cartesian_vector_derivatives_certified'] or packet['temporal_recursion']:
                raise ValueError('Unbuilt patch/full/temporal scope promoted')
            if packet['exact_terminal_closure_from_same_implicit_map']:
                for row in packet['actual_normalized_five_defects']:
                    if any(endpoints(read_interval(c,value))!=(mp.mpf(0),mp.mpf(0)) for value in row['coefficients']):
                        raise ArithmeticError('Functional terminal identities not exactly zero')
                if not packet['terminal_refinement_is_functional_identity_not_moment_reset']:
                    raise ValueError('Terminal closure proof scope lost')
        inlet=provider.evaluate(1,[-1,1]);old=provider.reference.terminal([-1,1],-6)
        if any(endpoints(a)!=endpoints(b) for a,b in zip(inlet['original_P0_axial5'],old['pressure_axis_axial5_coefficients'])):
            raise ArithmeticError('Actual inlet pressure datum changed')
        shape=old['actual_normalized_moment_shape_axial5_coefficients'];centered=old['actual_centered_moment_axial5_coefficients']
        convert=lambda row:IntervalTaylor(c,row)
        target=[convert(shape['mean']),convert(shape['theta']),convert(shape['theta_z']),
            convert(centered['axial_square'])*data['invAm2']-convert(shape['swirl'])/2,convert(shape['pressure'])/2]
        names=('Mz_over_R_axial5','Mtheta_over_sqrt2_Rm_1p5_Am_axial5','Mtheta_z_over_sqrt2_Rm_1p5_Am_axial5',
            'Mztheta_minus_8ZMz_plus_16Z2R_over_RmAm2_axial5','Mp_over_Am2_axial5')
        for row,name in zip(target,names):
            for a,b in zip(row.coefficients,inlet[name]):
                if max(endpoints(a)[0],endpoints(b)[0])>min(endpoints(a)[1],endpoints(b)[1]):
                    raise ArithmeticError('Actual unpatched Rm moment inconsistent with patch inlet')
        # This numerical consistency check complements the exact source
        # transport identities; overlap is not the source identity proof.
        for key in ('local_relaxed_cone_lower_margin','local_relaxed_bw_upper','local_relaxed_kappa_upper'):
            value=read_interval(c,receipt[key])
            if (endpoints(value)[0]<=0 if key.endswith('margin') else endpoints(value)[1]>=(mp.mpf('.1') if key.endswith('bw_upper') else 1)):
                raise ArithmeticError('Actual patch local cone theorem gate failed')
        result=dict(actual_five_defect_family_sha256=provider.family,implicit_source_sha256=provider.source,
            source_transport_and_implicit_structure=symbolic_source_transport(),
            original_restoration_integrals_directly_source_bound=True,
            independent_coefficient_fixture=independent_coefficient_fixture(provider),
            actual_source_tail_C1_caps_checked=tail_count,actual_source_factored_decay_caps_checked=proof_count,
            actual_implicit_coefficient_axial5_bounds_checked=coefficient_count,
            actual_velocity_pressure_moment_bounds_checked=profile_count,
            exact_actual_Rm_inlet_and_original_P0_bound=True,
            actual_five_functional_terminal_identities_connected=True,
            actual_five_moment_patch_connected=True,actual_patch_coefficient_axial5_enclosures_available=True,
            actual_patch_velocity_pressure_moment_axial5_enclosures_available=True,actual_patch_radial_recovery_axial4_available=True,
            actual_patch_local_relaxed_cone_theorem_bound=True,
            residual_zero_containment_not_used_as_closure_proof=True,
            radial_mixed4_certified=False,all_radial_endpoint_joins_certified=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,
            temporal_recursion=False,all_passed=True,
            input_hashes={**receipt['input_hashes'],NAME:hashlib.sha256((HERE/NAME).read_bytes()).hexdigest(),
                Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual source tails, implicit axial5 coefficients and terminal five identities PASS',flush=True)
    return result


if __name__=='__main__':run()
