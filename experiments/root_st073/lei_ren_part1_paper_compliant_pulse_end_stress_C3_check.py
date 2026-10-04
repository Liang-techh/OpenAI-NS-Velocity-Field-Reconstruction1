"""Independent full meridional pulse-stress and radial-source checks."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import (
    CompliantPulseEndStressC3,DOMAIN,original_formula_identities,pulse_coefficients)
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def independent_meridional_fixture():
    """Original bump integrals and independent (3.12)-(3.13) radial source.

    Finite moderate parameters test algebra/derivatives, not project cone,
    blow-up, global physical energy or final corrected NS residual.
    """
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=85
        mu,delta,xp,B0,D,R0=map(mp.mpf,('.2','.12','.3','.45','.07','2.8'))
        r=1-mu;bp=mp.mpf('.5')+mu;prate=1+2*mu;ell=mp.mpf('.15')
        normal=mp.quad(lambda x:mp.exp(-1/(1-x*x)),[-1,0,1])
        def beta(v):
            x=v/ell
            return mp.exp(-1/(1-x*x))/(ell*normal) if abs(x)<1 else mp.mpf(0)
        def C(z):return 1/(1+z*z)
        coeff=(lambda z:mp.mpf('.1')*(1+z),lambda z:-mp.mpf('.07')*(1-z*z))
        errors=[];count=0;nonzero_cross_terms=[]
        for sv,zv in ((mp.mpf(-3),mp.mpf('.3')),(mp.mpf(-1),mp.mpf('.5')),(mp.mpf(-2),mp.mpf('.3'))):
            R=R0*mp.exp(sv);B=B0*mp.exp(-bp*sv);H=mp.exp(-13*r/mu-r*sv)
            weights=[[],[],[]]
            for center in (-3,-1):
                lo=max(sv,center-ell);hi=center+ell
                cuts=[lo]+([mp.mpf(center)] if lo<center<hi else [])+[hi]
                if lo>=hi:
                    vals=[mp.mpf(0)]*3
                else:
                    vals=[mp.quad(lambda v,lam=lam:mp.exp(lam*(v-sv))*beta(v-center),cuts)
                          for lam in (mp.mpf('.5')-mu,mp.mpf('.5')-2*mu)]
                    vals.append(mp.quad(lambda v:mp.exp(2*mu*(sv-v))*beta(v-center)**2,cuts))
                for k,value in enumerate(vals):weights[k].append(value)
            shape=[sum((coeff[i](zv)*mp.diff(beta,sv-center,j) for i,center in enumerate((-3,-1))),mp.mpf(0))
                   for j in range(5)]
            def bh(z,j=0):
                return sum((coeff[i](z)*mp.diff(beta,sv-center,j) for i,center in enumerate((-3,-1))),mp.mpf(0))
            def m(z,k):return -sum(coeff[i](z)*weights[k][i] for i in range(2))
            def J(z):return sum(coeff[i](z)**2*weights[2][i] for i in range(2))
            def ev(z):return (mp.mpf('.8')+mp.mpf('.1')*z*z)*mp.exp(2*mu*sv)-mp.expm1(2*mu*sv)/(4*mu)
            def pressure(z):return -mp.mpf('.35')*C(z)**2*mp.exp(prate*sv)+C(z)**2*mp.expm1(prate*sv)/(2*prate)
            def jet(fn):return IntervalTaylor(c,[c.mpf(v) for v in mp.taylor(fn,zv,5)])
            cj=jet(C);zj=IntervalTaylor.variable(c,c.mpf(zv),5)
            Bh=[jet(lambda z,j=j:bh(z,j)) for j in range(5)]
            m1=[jet(lambda z:m(z,0))];m2=[jet(lambda z:m(z,1))]
            e0=[jet(ev)];loss=[jet(J)];P=[jet(pressure)]
            for j in range(4):
                m1.append(Bh[j]-m1[j]*c.mpf(mp.mpf('.5')-mu))
                m2.append(Bh[j]-m2[j]*c.mpf(mp.mpf('.5')-2*mu))
                square=sum((Bh[l]*Bh[j-l]*mp.binomial(j,l) for l in range(j+1)),cj*0)
                e0.append(e0[j]*c.mpf(2*mu)-(cj*0+mp.mpf('.5') if j==0 else cj*0))
                loss.append(loss[j]*c.mpf(2*mu)-square)
                P.append(P[j]*c.mpf(prate)+(cj*cj/2 if j==0 else cj*0))
            sectors=pulse_coefficients(c.mpf(delta),c.mpf(mu),zj,cj,c.mpf(xp),Bh,m1,m2,e0,loss,P)
            def direct(z,j,label):
                ut=B*C(z);uz=B*C(z)*D*bh(z)
                utz=B*mp.diff(C,z);uzz=B*D*(mp.diff(C,z)*bh(z)+C(z)*mp.diff(bh,z))
                uzy=B*C(z)*D*(bh(z,1)-bp*bh(z))
                X=1/r+(xp-1/r)*H
                moments_fn={
                    'theta':lambda zz:mp.sqrt(2)*R**mp.mpf('1.5')*B*C(zz)*X,
                    'z':lambda zz:R*B*C(zz)*D*m(zz,0),
                    'theta_z':lambda zz:mp.sqrt(2)*R**mp.mpf('1.5')*(B*C(zz))**2*D*m(zz,1),
                    'z_theta':lambda zz:R*(B*C(zz))**2*(ev(zz)-D**2*J(zz)),
                    'p':lambda zz:B*B*pressure(zz)}
                result=evaluate_mp_stress(mp.log(R),z,delta,Utheta=ut,Uz=uz,Utheta_y=-bp*ut,
                    Utheta_Z=utz,Uz_y=uzy,Uz_Z=uzz,moments={k:fn(z) for k,fn in moments_fn.items()},
                    moments_Z={k:mp.diff(fn,z) for k,fn in moments_fn.items()},
                    P=B*B*pressure(z),P_Z=B*B*mp.diff(pressure,z),precision=mp.mp.dps)
                if not j:return result['T_theta' if label=='theta' else 'T_z']
                if label=='theta':return result['N_theta']-result['I_theta']-(1+mu)*result['S_theta']
                shear_y=mp.sqrt(2/R)*B*C(z)*D*(bh(z,2)-2*bp*bh(z,1)+bp**2*bh(z))-.5*result['S_z']
                return result['N_z']-.5*result['I_z']+shear_y
            for label,parts in sectors.items():
                for j in range(2):
                    for n in range(3-j):
                        actual=c.mpf(0)
                        for part in parts.values():
                            rp,bpwr,dp,hp=part['mode']
                            factor=c.mpf(R)**rp*c.mpf(B)**bpwr*c.mpf(D)**dp*c.mpf(H)**hp/c.sqrt(2)
                            actual+=part['full_derivative_rows'][j][n]*mp.factorial(n)*factor
                        expected=mp.diff(lambda z:direct(z,j,label),zv,n)
                        lo,hi=endpoints(actual);miss=max(lo-expected,expected-hi,mp.mpf(0))
                        if miss>mp.mpf('1e-55'):raise ArithmeticError('Independent full meridional pulse stress failed: '+label+' '+str((sv,j,n,miss)))
                        count+=1;errors.append(str(miss))
            nonzero_cross_terms.append(dict(s=str(sv),Bhat=str(bh(zv)),Mz_hat=str(m(zv,0)),Mtheta_z_hat=str(m(zv,1)),quadratic_loss=str(J(zv))))
        if any(mp.mpf(v['Mz_hat'])==0 or mp.mpf(v['Mtheta_z_hat'])==0 for v in nonzero_cross_terms):
            raise ArithmeticError('Fixture must preserve nonzero meridional histories')
        return dict(passed=True,independent_full_stress_and_radial_source_mixed_checks=count,
            tolerance='1e-55',misses=errors,inside_both_supports_and_between_supports=nonzero_cross_terms,
            meridional_histories_nonzero=True,signed_original_memory_nonzero=True,
            fixture_is_actual_project_cone_or_NS_validation=False)


def run():
    with mp.workdps(300):
        name=PREFIX+'pulse_end_stress_C3.json';record=json.loads((HERE/name).read_bytes())
        for source,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                raise ValueError('Pulse end source changed: '+source)
        field=CompliantPulseEndStressC3();fresh=encode(pack(field.report()))
        if fresh!=record:raise ValueError('Current full meridional pulse stress report differs')
        if record['domain']!=DOMAIN:raise ValueError('Original pulse-end domain changed')
        if not record['reference_amplitude_bridge']['pulse_flatten_reference_amplitude_AST_bridge_verified']:
            raise ValueError('Same native B0/reference pressure scale must be AST-bound')
        if not record['source_bound_full_formula_proof']['all_meridional_cross_terms_retained']:
            raise ValueError('Full original stress formulas required')
        count=0
        for point in record['samples']+[record['whole_original_end']]:
            for parts in point['full_meridional_stress_log_sectors'].values():
                for part in parts.values():
                    for row in part['full_stress_mixed3_coefficient_enclosures'].values():
                        if not all(mp.isfinite(mp.mpf(row[key])) for key in ('lower','upper')):
                            raise ArithmeticError('Pulse stress coefficient not finite')
                        count+=1
        for flag in ('pulse_flatten_full_stress_join_verified','pulse_end_physical_decomposition_constructed',
            'pulse_end_cone_certified','global_admissible_stress_lift_constructed',
            'physical_energy_integral_certified','temporal_recursion','source_caps_used_as_defining_field_values'):
            if record[flag]:raise ValueError('Pulse-end stress scope overclaimed: '+flag)
        fixture=independent_meridional_fixture()
        hashes=dict(record['input_hashes']);hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        result=dict(all_passed=True,actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
            domain=DOMAIN,current_full_meridional_stress_source_report_recomputed=True,
            finite_stress_sector_mixed3_rows=count,actual_original_pulse_end_similarity_stress_recovered=True,
            pulse_flatten_reference_amplitude_AST_bridge_verified=True,original_absolute_pressure_datum_preserved=True,
            whole_original_end_domain_retained=True,all_meridional_cross_terms_retained=True,
            source_caps_used_as_defining_field_values=False,independent_meridional_fixture=fixture,
            pulse_flatten_full_stress_join_verified=False,pulse_end_physical_decomposition_constructed=False,
            pulse_end_cone_certified=False,global_admissible_stress_lift_constructed=False,
            physical_energy_integral_certified=False,temporal_recursion=False,input_hashes=hashes)
        Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
        print('PASS whole original pulse-end full meridional similarity stress; functional join/physical/cone pending',flush=True)
        return result


if __name__=='__main__':run()
