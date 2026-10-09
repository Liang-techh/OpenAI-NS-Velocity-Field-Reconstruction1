"""Factored Section 11 fixed-phase Z derivatives for varying source shear.

The caller supplies enclosures of defining functions and ordinary Z jets,
not chart caps selected as point values. eta and d_star are the same frozen
positive source constants as in the C0 inverse. This is a local derivative
backend; it does not close the active-patch inherited source histories.
"""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_O2_conditioned_slow_Z as old

base=old.base;current=base.current;conditioned=base.conditioned;prior=base.prior
HERE,PREFIX,sha=old.HERE,old.PREFIX,old.sha;ep=old.ep;square=current.square
ZERO,ZROW=(0,0),(0,1)
NAME=PREFIX+'current_general_conditioned_slow_Z.json'
RECEIPT=PREFIX+'current_general_conditioned_slow_Z_check.json'
GATE='general_source_factored_phase_held_Z_primitives_with_varying_shear_installed'
OPEN=('actual_active_patch_genuine_coefficients_installed',
      'actual_complete_patch_integral_installed','actual_upstream_histories_closed',
      'continuous_Z_source_integral_contract_installed',
      'numerical_original_source_point_or_integral_oracle_installed',
      'actual_five_controls_installed','current_whole_N_selected',
      *base.point.source.inertial.profiles.loop.OPEN)


def exact_identities():
    r,s,q,qZ,rZ,chi,psi=sy.symbols('r s q qZ rZ chi psi',nonzero=True,real=True)
    J1=q*sy.sqrt(s)*(chi-psi)/r
    H=(2-3*s)*(chi-psi)+2*r*sy.sin(chi)
    J2osc=q*q*H/r**2
    chiZ=2*rZ*sy.sin(chi)/s;SZ=-2*r*rZ
    derivative=lambda F:qZ*sy.diff(F,q)+rZ*(sy.diff(F,r)-2*r*sy.diff(F,s)+2*sy.sin(chi)/s*sy.diff(F,chi))
    want1=qZ*sy.sqrt(s)*(chi-psi)/r+q*sy.sqrt(s)/r*(chiZ-rZ*(chi-psi)/(r*s))
    want2=2*q*qZ*H/r**2+q*q/r**2*(-3*SZ*(chi-psi)+2*rZ*sy.sin(chi)
        +chiZ*(2-3*s+2*r*sy.cos(chi))-2*rZ*H/r)
    assert sy.simplify((derivative(J1)-want1).subs(s,1-r*r))==0
    assert sy.simplify(derivative(J2osc)-want2)==0
    for k in range(3,9):
        G=r**k+s*(k-1)*r**(k-2)/2
        expected=r**(k-1)+(k-1)*(k-2)*s*r**(k-3)/2
        assert sy.expand((sy.diff(G,r)-2*r*sy.diff(G,s)-expected).subs(s,1-r*r))==0
    uZ=sy.symbols('uZ',real=True)
    W1=sum(r**(k-1)*sy.sin(k*psi)/k for k in range(1,4))
    V2=sum((r**k+(1-r*r)*(k-1)*r**(k-2)/2 if k>=2 else r)*sy.sin(k*psi)/k for k in range(1,4))
    smallJ1=2*q*sy.sqrt(1-r*r)*W1;smallJ2=4*q*q*V2
    smallD=lambda F:(qZ*sy.diff(F,q)+uZ*sy.diff(F,r)).subs(r,0)
    assert sy.simplify(smallJ1.subs(r,0)-2*q*sy.sin(psi))==0
    assert sy.simplify(smallJ2.subs(r,0)-q*q*sy.sin(2*psi))==0
    assert sy.simplify(smallD(smallJ1)-(2*qZ*sy.sin(psi)+q*uZ*sy.sin(2*psi)))==0
    assert sy.simplify(smallD(smallJ2)-(2*q*qZ*sy.sin(2*psi)+4*q*q*uZ*(sy.sin(psi)+sy.sin(3*psi)/3)))==0
    a,aZ,t0,t0Z,E,EZ,phi,psiZ,J,JZ,OZ=sy.symbols('a aZ t0 t0Z E EZ phi psiZ J JZ OZ',real=True)
    A=a*(phi-psi/(2*sy.pi))/2;AZ=aZ*A/a-a*psiZ/(4*sy.pi)
    nuZ=2*t0*t0Z+4*q*qZ
    original=2*sy.pi*phi*nuZ-(2*t0*t0Z*psi+2*t0Z*J+2*t0*JZ+4*q*qZ*psi+OZ)
    factored=4*sy.pi*A/a*nuZ-2*t0Z*J-2*t0*JZ-OZ
    assert sy.expand(original-factored)==0
    b=-a*t0;bZ=-aZ*t0-a*t0Z;t=t0+sy.Symbol('osc_t',real=True)
    T1=t0*psi+J;T1Z=t0Z*psi+JZ
    M=-a*T1/(2*sy.pi)-b*phi
    originalB=EZ*M/2+E*(-aZ*T1/(2*sy.pi)-a*(T1Z+t*psiZ)/(2*sy.pi)-bZ*phi)/2
    factoredB=EZ*(t0*A-a*J/(4*sy.pi))+E*(t0Z*A+t0*AZ-(aZ*J+a*(JZ+(t-t0)*psiZ))/(4*sy.pi))
    assert sy.expand(originalB-factoredB)==0
    return dict(passed=True,general_Mobius_J1_and_secular_free_J2_Z=True,
        general_fixed_phase_implicit_numerator=True,general_A_and_B_Z_product_rules=True,
        small_r_Fourier_derivative_coefficients=True,
        exact_r_zero_J1_J2_osc_and_both_Z_continuity_identities=True,
        implicit_identity='psi_Z=(2*pi*phi*nu_Z-T2_Z)/(1+t^2)',
        factored_numerator='4*pi*(A/a)*nu_Z-2*t0_Z*J1-2*t0*J1_Z-J2_osc_Z',
        original_eta_and_dstar_Z_independent=True)


def shear_source_Z(roots,eta_log,log_a_lower):
    """Differentiate the original q recipe, including its flat cutoff."""
    a=roots['a'][ZERO];aZ=roots['a'][ZROW];b=roots['b'][ZERO];bZ=roots['b'][ZROW]
    c=a.ctx;diva=lambda value:value.positive_divide(a,log_a_lower)
    t0=diva(-b);t0Z=diva(-bZ)-diva(t0*aZ)
    Delta=a+diva(square(b))-2
    kappaZ=aZ+diva(b*bZ*2)-diva(diva(square(b))*aZ)
    loop=current.q_enclosure(a,Delta,eta_log,log_a_lower)
    if loop['branch']=='requires_source_box_refinement':
        raise ArithmeticError('Split the source box at the original active/flat q boundary')
    q=loop['q']
    if loop['branch']=='flat':
        qZ=a.scalar(0);cutoff='exact_flat_source_q_and_q_Z'
    else:
        gamma=loop['positive_active_gamma_enclosure']
        root=current.nonnegative_sqrt(gamma.positive_divide(a*2,c.mpf(log_a_lower)+c.ln(2)))
        rootZ=-root*(kappaZ.positive_divide(gamma,eta_log)+diva(aZ))*c.mpf('.5')
        if Delta.zero or ep(Delta.coefficient)[1]<=0:
            sigma=c.mpf(1);sigmaZ=a.scalar(0);cutoff='exact_sigma_one_with_zero_source_derivative'
        else:
            eta=prior.ScaledEnclosure(prior.FormalScale(a.scale.bases,offset=eta_log),1,a.ledger)
            argument=1-conditioned.bounded_value(Delta.positive_divide(eta,eta_log))
            jets=prior.sigma_jets(c,argument);sigma=jets[0]
            sigmaZ=(-kappaZ).positive_divide(eta,eta_log)*jets[1]
            cutoff='original_directed_sigma_and_sigma_prime'
        if (root*sigma).record()!=q.record():
            raise ArithmeticError('Refine source cutoff; C0 q and differentiated root*sigma must be identical')
        qZ=root*sigmaZ+rootZ*sigma
    return dict(q=q,q_Z=qZ,t0=t0,t0_Z=t0Z,kappa_minus2=Delta,kappa_Z=kappaZ,
        loop=loop,cutoff_contract=cutoff)


def small_series(kernel,psi):
    c=kernel.c;r=kernel.r;s=kernel.s;R=max(abs(v) for v in ep(r));M=48
    if R>mp.mpf('.25'):raise ArithmeticError('Refine source before the small-r derivative series')
    W1=W1r=V2=S2=c.mpf(0)
    for k in range(1,M+1):
        sine=c.sin(k*psi)/k
        W1+=r**(k-1)*sine
        if k>=2:W1r+=(k-1)*r**(k-2)*sine
        G=r**k
        if k>=2:G+=s*(k-1)*r**(k-2)/2
        V2+=G*sine
        Gp=r**(k-1)
        if k>=3:Gp+=s*(k-1)*(k-2)*r**(k-3)/2
        S2+=Gp*sine
    if R:
        rc=c.mpf(R)
        tails=(rc**M/((M+1)*(1-rc)),rc**(M-1)/(1-rc),
            rc**(M+1)/((M+1)*(1-rc))+rc**(M-1)/(2*(1-rc)),
            rc**M/((M+1)*(1-rc))+rc**(M-2)/2*((M+1)/(1-rc)+rc/(1-rc)**2))
        W1,W1r,V2,S2=(v+c.mpf((-ep(tail)[1],ep(tail)[1]))
            for v,tail in zip((W1,W1r,V2,S2),tails))
    return W1,W1r,V2,S2


class GeneralConditionedSlowZ:
    """Enclosure backend, conditional on caller-owned source hypotheses.

    A common arithmetic frame is checked here; it does not establish a
    source family or validate a paper cone or an auxiliary-scale recipe.
    The caller must separately bind the defining source functions, their
    ordinary jets, frozen constants and valid positive a lower theorem.
    original_O2_records performs that binding for its genuine O2 inputs.
    """
    mode='general_factored_phase_held_Z_enclosure_backend'
    def __init__(self,roots,*,eta_log,dstar_log,log_a_lower):
        a=roots['a'][ZERO];self.c=c=a.ctx;self.roots={key:dict(row) for key,row in roots.items()}
        for key in ('a','b','p2','E'):
            for order in (ZERO,ZROW):
                value=self.roots[key][order]
                if not isinstance(value,prior.ScaledEnclosure) or value.ctx is not c or value.scale.bases is not a.scale.bases or value.ledger is not a.ledger:
                    raise ValueError('One source context, fixed factor basis and directed ledger required')
        for value in (eta_log,dstar_log,log_a_lower):
            if hasattr(value,'ctx') and value.ctx is not c:raise ValueError('Source auxiliary logarithms require the same context')
            if any(not mp.isfinite(x) for x in ep(c.mpf(value))):raise ValueError('Finite frozen source logarithms required')
        self.shear=shear_source_Z(self.roots,c.mpf(eta_log),c.mpf(log_a_lower))
        self.roots['t0']={ZERO:self.shear['t0'],ZROW:self.shear['t0_Z']}
        self.query=dict(q=self.shear['q'],roots=self.roots)
        self.kernel=conditioned.ConditionedPhase(self.query,c.mpf(dstar_log))
        if self.kernel.geometry=='requires_signed_source_refinement':raise ArithmeticError('Refine source sign/geometry before Z differentiation')
        self.log_a_lower=c.mpf(log_a_lower)

    def values(self,coordinate,chart='psi'):
        """Z rows at Phi(coordinate,chart); phi is held constant in Z.

        When coordinate is an issued inverse bracket, these enclose the Z
        derivatives at every fixed phase represented by that bracket. No
        derivative of a selected inverse interval endpoint is taken.
        """
        k=self.kernel;c=self.c;scalar=k.scalar;roots=self.roots;sh=self.shear
        if k.source is not self.query or k.roots is not roots:raise ValueError('Same live C0 inverse and defining source required')
        if chart not in (('psi','E') if k.geometry=='signed_Mobius' else ('psi',)):
            raise ValueError('Angle chart is not defined for this source geometry')
        x=c.mpf(coordinate);lo,hi=ep(x)
        if not 0<=lo<=hi<=1:raise ValueError('Closed-period coordinate box required')
        q=k.q;qZ=sh['q_Z'];a=k.a;aZ=roots['a'][ZROW];t0=k.t0;t0Z=sh['t0_Z']
        uZ=(roots['p2'][ZROW]*q+roots['p2'][ZERO]*qZ).positive_divide(k.dstar,k.dstar.scale.evaluate())
        proof=dict(q_Z=qZ.record(),t0_Z=t0Z.record(),kappa_Z=sh['kappa_Z'].record(),u_Z=uZ.record(),
            original_cutoff_derivative=sh['cutoff_contract'],a_Z_b_Z_q_Z_t0_Z_terms_retained=True,
            original_E_Z_product_term_retained=True,positive_implicit_denominator_lower=1,
            differentiated_defining_functions_not_interval_endpoints=True,
            source_owner_cone_scale_and_a_lower_hypotheses_require_caller_binding=True,
            C0_q_and_differentiated_original_root_sigma_identical=True,
            native_active_patch_source_histories_closed=False)
        if k.flat or lo==hi and lo in (0,mp.mpf('.5'),1):
            return dict(psi_Z=scalar(0),A_Z_slow=scalar(0),B_Z_slow=scalar(0)),dict(proof,
                branch='exact_general_flat_or_symmetry')
        psi_fraction,chi_fraction=k.angles(x,chart);psi=2*c.pi*psi_fraction
        if k.geometry=='small_r_series':
            r=k.r;s=k.s;hinv=k.hinv
            W1,W1r,V2,S2=small_series(k,psi)
            J1=q*hinv*(2*W1)
            J1Z=qZ*hinv*(2*W1)+q*uZ*(2*(-r*s*W1+s*s*W1r))
            rZ=uZ*hinv*s
            J2oscZ=q*qZ*(8*V2)+square(q)*rZ*(4*S2)
            D=1-2*r*c.cos(psi)+r*r
            if ep(D)[0]<=0:raise ArithmeticError('Positive small-r denominator lost')
            osc_t=q*hinv*(2*(c.cos(psi)-r)/D)
            branch='general_exact_midplane' if k.u.zero else 'general_small_r_Fourier'
        else:
            r=k.r;rho=k.rho;s=k.s_source;hinv=k.hinv
            chi=2*c.pi*chi_fraction;K=uZ*hinv;rZ=K*s;SZ=-rZ*(2*r)
            difference=chi-psi
            if chart=='E':
                sinchi=scalar(c.sin(chi));chiZ=K*(2*c.sin(chi))
                if k.sign>0:
                    NE=scalar(2*c.cos(chi/2)**2)-rho
                    lam=scalar(4*c.cos(chi/2)**2)-rho*(2*c.cos(chi))-s*3
                else:
                    NE=rho-scalar(2*c.sin(chi/2)**2)
                    lam=scalar(4*c.sin(chi/2)**2)+rho*(2*c.cos(chi))-s*3
                osc_t=(q*NE*2).positive_divide(hinv,hinv.scale.evaluate()+c.ln(c.mpf(ep(hinv.coefficient)[0])))
            else:
                half=c.sin(psi/2)**2 if k.sign>0 else c.cos(psi/2)**2
                D=square(rho)+scalar(4*abs(r)*half)
                lower=2*(rho.scale.evaluate()+c.ln(c.mpf(ep(rho.coefficient)[0])))
                sinchi=(s*c.sin(psi)).positive_divide(D,lower)
                chiZ=(rZ*(2*c.sin(psi))).positive_divide(D,lower)
                numerator=rho-scalar(2*half) if k.sign>0 else scalar(2*half)-rho
                osc_t=(q*hinv*numerator*2).positive_divide(D,lower)
                lam=(s*(2*(1-r*c.cos(psi)))).positive_divide(D,lower)-s*3
            Hosc=scalar(2*difference)-s*(3*difference)+sinchi*(2*r)
            J1=q*hinv*(difference/r)
            J1Z=qZ*hinv*(difference/r)+q*hinv*(1/r)*(chiZ-K*(difference/r))
            J2oscZ=q*qZ*(2/r**2)*Hosc+square(q)*(1/r**2)*(
                SZ*(-3*difference)+rZ*sinchi*2+chiZ*lam-rZ*Hosc*(2/r))
            branch='general_signed_'+chart
        C0=k.primitives(x,chart);A=C0['A'];diva=lambda value:value.positive_divide(a,self.log_a_lower)
        nuZ=t0*t0Z*2+q*qZ*4
        numerator=diva(A)*(4*c.pi)*nuZ-t0Z*J1*2-t0*J1Z*2-J2oscZ
        t=t0+osc_t;psiZ=numerator.positive_divide(scalar(1)+square(t),0)
        AZ=diva(aZ)*A-a*psiZ*(1/(4*c.pi))
        Mhalf=t0*A-a*J1*(1/(4*c.pi))
        BZ=roots['E'][ZROW]*Mhalf+k.E*(t0Z*A+t0*AZ-(aZ*J1+a*(J1Z+osc_t*psiZ))*(1/(4*c.pi)))
        proof.update(branch=branch,J1_fixed_angle=J1.record(),J1_Z_fixed_angle=J1Z.record(),
            J2_osc_Z_fixed_angle=J2oscZ.record(),nu_Z=nuZ.record(),direction=t.record(),
            secular_terms_cancelled_by_exact_source_identity=True,
            positive_rho_s_hinv_retained_as_source_factors=True)
        return dict(psi_Z=psiZ,A_Z_slow=AZ,B_Z_slow=BZ),proof

    def evaluate(self,phase,bits=80):
        with mp.workdps(self.c.dps+40):
            C0=self.kernel.evaluate(phase,bits=bits)
            if C0['status']!='enclosed':raise ArithmeticError('Source inverse needs refinement')
            selected=C0['selected_inverse']
            values,proof=self.values(selected['coordinate_interval'],selected['chart'])
        return dict(mode=self.mode,C0=C0,slow_Z={key:value.record() for key,value in values.items()},
            derivative_contract=proof,free_fixed_phase_parameter=True,
            source_owner_binding_not_established_by_generic_constructor=True,
            genuine_native_active_patch_source_installed=False)


def fixture(values,derivatives,*,dps=140,bases_values=None,source_powers=None):
    """Explicit diagnostic source functions only; never native parameters."""
    c=MPIntervalContext();c.dps=dps
    bases=tuple(c.mpf(v) for v in (bases_values or (0,0,0,0,0)))
    ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
        positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    roots={};source_powers=source_powers or {}
    for key in ('a','b','p2','E'):
        roots[key]={}
        for order,v in ((ZERO,values[key]),(ZROW,derivatives[key])):
            scale=prior.FormalScale(bases,source_powers.get((key,order),(0,0,0,0,0)))
            roots[key][order]=prior.ScaledEnclosure(scale,v,ledger)
    return c,roots


def original_O2_records():
    owner=old.OriginalO2ConditionedSlowZ();rows=[]
    for saved in owner.saved['actual_original_O2_conditioned_point_queries']:
        y=saved['original_y_exact'];Z=saved['original_Z_exact'];N=saved['explicit_candidate_N']
        query=owner.owner.query(y=y,Z=Z);c=owner.owner.ctx;logs=owner.owner.scales.logs
        new=GeneralConditionedSlowZ(query['roots'],eta_log=logs['eta'],dstar_log=logs['d_star'],log_a_lower=logs['a_min'])
        with mp.workdps(c.dps+40):
            if new.kernel.q.record()!=query['kernel'].q.record():raise ValueError('Same original O2 q source required')
            phase=owner.owner.radius.evaluate(y=y,N=N);items=[]
            for box in phase['true_original_phase_directed_boxes']:
                result=new.evaluate(c.mpf([box['lower'],box['upper']]))
                result['C0'].update(free_phase_parameter_not_spatial_phase=False,original_common_N_and_radius_phase_bound=True)
                result['free_fixed_phase_parameter']=False
                items.append(result)
        rows.append(dict(y=y,Z=Z,N=N,general_derivatives_at_genuine_original_O2_points=items,
            defining_point_coefficients_and_errors=saved['actual_original_O2_factored_inputs'],
            actual_original_radius_phase=phase,source_family=owner.family,
            original_source_basis_contract=query['basis_contract'],
            accepted_source_owner_and_auxiliary_logs_bound=True,
            native_caps_or_midpoints_selected_as_field_values=False))
    return owner,rows


def run():
    began=time.monotonic();owner,rows=original_O2_records()
    accepted=json.loads((HERE/old.RECEIPT).read_bytes())
    if not accepted.get('all_passed') or not accepted.get(old.GATE) or accepted['source_family']!=owner.family:
        raise ValueError('Same accepted original O2 reduced derivative receipt required')
    hashes=dict(owner.hashes)
    for name,digest in accepted['input_hashes'].items():
        if sha(name)!=digest or name in hashes and hashes[name]!=digest:
            raise ValueError('Original O2 derivative source dependencies disagree: '+name)
        hashes[name]=digest
    report=dict(**{GATE:True},source_family=owner.family,mode=GeneralConditionedSlowZ.mode,
        exact_general_derivative_identities=exact_identities(),
        genuine_original_O2_general_backend_queries=rows,
        accepted_original_O2_source_scale_contract=owner.owner.scales.record(),
        original_source_point_errors_and_auxiliary_factors_retained=True,
        no_ancestor_quadratures_or_producers_executed=True,
        active_patch_inherited_coefficients_still_require_defining_source_evaluation=True,
        generic_backend_requires_caller_owned_source_cone_and_scale_hypotheses=True,
        genuine_O2_path_binds_accepted_point_errors_family_logs_and_true_phase=True,
        **dict.fromkeys(OPEN,False),input_hashes={**hashes,Path(__file__).name:sha(Path(__file__).name),
            old.RECEIPT:sha(old.RECEIPT)},execution_seconds=time.monotonic()-began,
        scope='General factored fixed-phase ordinary-Z backend, connected to genuine accepted O2 source inputs. Varying-shear cases require independent focused checks. Active-patch source histories, full integrals and full reconstruction remain open.')
    (HERE/NAME).write_text(json.dumps(base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('General varying-shear factored Z backend connected to genuine O2 sources',flush=True)
    return report


if __name__=='__main__':run()
