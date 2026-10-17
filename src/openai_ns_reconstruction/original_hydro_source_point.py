"""Same-generator hydro source from saved core jets and finite-width transport.

Only new bounded algebra is evaluated.  The formal epsilon sectors survive
before physical scales are combined; saved mass and ancestor computations
are never replayed.  This is not the complete pressure/hydro/swirl source.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_abs, mpf_cmp

from .original_leading_viscous_kernel import decode_pair, encode_pair
from .original_modal_source_point import _function_ast_guard
from .original_signed_C1_source import macro_nonlinear_enclosure, main_coefficient
from .original_source_jet_kernel import ordinary_product, ordinary_quotient

ROOT=Path(__file__).resolve().parents[2]
DIRECTORY=ROOT/'experiments/root_st073'
PARENT='lei_ren_part1_paper_compliant_current_original_signed_C1_source_numeric.json'
PARENT_SHA256='343d4b84e23374132467763850839f426efc3e80e97859bf6acf7db23f15723d'
FIRST='lei_ren_part1_paper_compliant_current_original_first_collar_jets.json'
BRIDGE='lei_ren_part1_paper_compliant_actual_bridge_integrals.json'
DOMAIN='lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain.json.gz'
FIXED='lei_ren_part1_paper_shared_bump_constants.json'
MASS='lei_ren_part1_paper_compliant_current_original_source_jet_mass_numeric.json'
MIXED='lei_ren_part1_paper_compliant_bridge_mixed_C4.json'
CORE='lei_ren_part1_paper_compliant_current_original_core_numeric_approximation.json'
MICRO='lei_ren_part1_paper_compliant_microswitch_mixed_C4.json'


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':'),ensure_ascii=False).encode()).hexdigest()


def load(path):
    if path.suffix=='.gz':
        with gzip.open(path,'rt',encoding='utf-8') as stream:return json.load(stream)
    return json.loads(path.read_bytes())


def exact_iv(ctx,iv,row):
    raw=row.get('source_enclosure_exact_mpf_tuples')
    if raw is None:raw=[row[k+'_exact_mpf_tuple'] for k in ('lower','upper')]
    return iv.mpf([ctx.make_mpf(tuple(v)) for v in raw])


def upper(ctx,value):
    rows=value._mpi_ if hasattr(value,'_mpi_') else [value._mpf_]*2
    a,b=(mpf_abs(v) for v in rows)
    return ctx.make_mpf(a if mpf_cmp(a,b)>=0 else b)


def ordinary_norm(ctx,iv,rows,weight):
    return sum((iv.mpf(upper(ctx,row))*weight**k/math.factorial(k)
                for k,row in enumerate(rows)),iv.mpf(0))


def width_errors(ctx,iv,coefficient,weight,logh,*,orders,source_proved,error_cap='1e-1000'):
    """Bound hb times a jet norm before admitting each ordinary error cap."""
    if source_proved is not True:raise ValueError('A same-source transport proof is required')
    threshold=iv.mpf(error_cap)
    if threshold.a<=0 or threshold.b>=1:raise ValueError('A positive error cap smaller than one is required')
    if coefficient.a<0 or weight.a<=0:raise ValueError('Positive norm bound and jet weight required')
    errors=[];proofs=[]
    for k in range(orders+1):
        bound=coefficient*math.factorial(k)/weight**k
        combined=logh+iv.log(bound)
        if mpf_cmp(combined._mpi_[1],iv.log(threshold)._mpi_[0])>=0:
            raise ValueError('Original width fails the requested transport error cap')
        errors.append(iv.mpf([-1,1])*threshold)
        proofs.append({'ordinary_Z_order':k,'hb_coefficient_bound_exact_mpf_tuples':[list(v) for v in bound._mpi_],
                       'combined_original_log_width_bound_exact_mpf_tuples':[list(v) for v in combined._mpi_],
                       'representable_absolute_error_cap':error_cap,'source_width_selected':False,'passed':True})
    return errors,proofs


def finite_width_comparison(ctx,iv,core,domain,weight,hcap,logh,*,source_proved,error_cap='1e-1000'):
    """Ordinary normalized psi/mpsi/apsi comparison jets at native phase 2.

    Positive IBP weights bound psi by the same core tube.  The own M/A
    ODEs give exact normalized equations m_y=psi-m, a_y=psi^2-a.
    Their positive exponential kernels apply in the weighted jet algebra.
    """
    if source_proved is not True:raise ValueError('Same-generator comparison/own-moment proofs required')
    if hcap.a<=0 or (4*iv.exp(2*hcap)).b>=iv.mpf('4.1').a:
        raise ValueError('The original smoothing must remain in the same core continuation tube')
    d=domain['joint_core_domain']
    rz=exact_iv(ctx,iv,d['Cauchy_Z_radius']);rr=iv.mpf(d['Cauchy_rho_radius'])
    major=exact_iv(ctx,iv,d['Xh_norm_majorants']['Psi'])*d['core_value_factor']
    # All required ordinary derivatives through Z3 are covered on the same
    # continuation tube; the controls and original hb are Z-independent.
    psi_cover=major*sum((weight**k/rz**k for k in range(4)),iv.mpf(0))
    m0=ordinary_norm(ctx,iv,core['Mpsi'],weight)
    a0=ordinary_norm(ctx,iv,core['Apsi'],weight)
    changes={'Psi':8*iv.exp(2*hcap)*psi_cover/rr,
             'Mpsi':2*(psi_cover+m0),'Apsi':2*(psi_cover**2+a0)}
    truth={};proofs={}
    for name,bound in changes.items():
        errors,proof=width_errors(ctx,iv,bound,weight,logh,orders=3,source_proved=source_proved,error_cap=error_cap)
        truth[name]=[v+e for v,e in zip(core[name],errors)];proofs[name]=proof
    _,rho_proof=width_errors(ctx,iv,8*iv.exp(2*hcap),weight,logh,orders=0,source_proved=source_proved,error_cap=error_cap)
    rho=iv.mpf(4)+iv.mpf([0,1])*iv.mpf(error_cap)
    return truth,rho,{'passed':True,'phase':'2','normalized_own_ODEs':['m_y=psi-m','a_y=psi^2-a'],
                      'positive_IBP_and_ODE_kernels':True,'fixed_Z_independent_width_and_controls':True,
                      'ordinary_error_proofs':proofs,'rho_error_proof':rho_proof,
                      'source_width_materialized_or_selected':False,'core_recomputed':False}


def four_product(left,right):
    if len(left)!=4 or len(right)!=4:raise ValueError('Ordinary Z0..3 rows required')
    return [sum((math.comb(k,j)*left[j]*right[k-j] for j in range(k+1)),left[0]*0)
            for k in range(4)]


class OriginalHydroSourcePoint:
    def __init__(self,*,directory=DIRECTORY):
        self.directory=Path(directory)
        if sha(self.directory/PARENT)!=PARENT_SHA256:raise ValueError('Pinned accepted pressure parent required')
        self.parent=load(self.directory/PARENT)
        if self.parent['checks']['passed'] is not True:raise ValueError('Passed parent required')
        self.input_hashes=dict(self.parent['input_hashes']);self.input_hashes[PARENT]=PARENT_SHA256
        self.implementation_hashes=dict(self.parent['implementation_hashes'])
        if any(name not in self.input_hashes for name in (FIRST,BRIDGE,DOMAIN,FIXED,MASS,MIXED,CORE,MICRO)):
            raise ValueError('Every hydrated source must have an accepted input pin')
        for base,group in ((self.directory,self.input_hashes),(ROOT,self.implementation_hashes)):
            for name,expected in group.items():
                if sha(base/name)!=expected:raise ValueError('Changed accepted pin: '+name)
        self.first=load(self.directory/FIRST);self.bridge=load(self.directory/BRIDGE)
        self.domain=load(self.directory/DOMAIN);self.fixed=load(self.directory/FIXED)
        self.mass=load(self.directory/MASS)
        self.mixed=load(self.directory/MIXED);self.core=load(self.directory/CORE)
        self.micro=load(self.directory/MICRO)
        fractions=[Fraction(v['exact_original_j_fraction']) for v in self.core['levels']]
        if not fractions or any(v!=fractions[0] for v in fractions):raise ValueError('Saved exact axis source changed')
        self.axis_j=fractions[0]
        if digest(self.mass['interior_mass'])!=self.parent['inherited_mass_packet_digest']:
            raise ValueError('Saved mass packet changed')
        self.ctx=mp.mp.clone();self.ctx.dps=180
        self.iv=MPIntervalContext();self.iv.dps=220
        self.weight=exact_iv(self.ctx,self.iv,self.bridge['axial_jet_weight'])
        self.hcap=exact_iv(self.ctx,self.iv,self.bridge['width_enclosure_cap'])
        self.logh=exact_iv(self.ctx,self.iv,self.bridge['source_log_hb_enclosure'])

    def source_guards(self):
        base=ROOT/'experiments/root_st073';src=ROOT/'src/openai_ns_reconstruction'
        specs=(
            ('hydro_modes',base/'lei_ren_part1_paper_compliant_macro_signed_integrals.py','_mode_rows',None,
             ["[base['drive_hydro'], one['drive_hydro']-base['drive_hydro'], _zero(base['drive_hydro'])]"]),
            ('mode_targets',base/'lei_ren_part1_paper_compliant_macro_signed_integrals.py','_target_moments',None,
             ["{MTH:phi,MZ:value,MTHZ:phi*value,MZT:{'axial':value*value,'swirl':phi*phi/2},MP:phi*phi}"]),
            ('hydro_direction',base/'lei_ren_part1_paper_compliant_inner_bridge_profiles.py','direction',None,
             ["W=1-(z*m)*(1-delta)-d*derivative(m)",
              "(-W*v+(m-z*derivative(m))*((1-delta)/2)+(z*a)*(2*delta)-d*derivative(a))/(2*L)"]),
            ('macro_origin',base/'lei_ren_part1_paper_compliant_actual_bridge_integrals.py','prepare','CompliantActualBridgeIntegrals',
             ["end=self.comparison.micro(Z,2,root)","modes=_mode_rows(self.direction_context,z,inputs,phi2,V2,moment2)",
              "R0=self.bridge.r*c.exp(c.mpf([0,endpoints(2*self.cap)[1]]))"]),
            ('macro_weights',src/'original_macro_pressure_rows.py','source_Uy','OriginalMacroPressureRows',
             ["modes=[prepared['R0']**j*R**(power-j) for j in range(3)]",
              "row=-prepared['quotient']*module.weighted_modes(prepared['modes']['drive_'+name],modes)"]),
            ('core_normalized_moments',src/'original_first_collar_jets.py','inlet_jets','OriginalFirstCollarJets',
             ["psi=[J.constant(c,0,order)]+[J(c,rows['Uz'][n][:order+1])/eps for n in range(1,degree+1)]",
              "linear(psi,0,4)","product(psi,psi,0,4)"]),
            ('fixed_axis_j',src/'original_core_numeric.py','__init__','OriginalCoreNumeric',
             ["j=Fraction(1,16000000*KN*CS*CA)","self.j_fraction=j"]),
            ('comparison_field_origin',base/'lei_ren_part1_paper_compliant_comparison_point_integrals.py','inlet','CompliantComparisonPointIntegrals',
             ["packet=self.atoms.field.build_root_rows(24,6) if root else self.atoms.rebuild.rebuild(Z,24,6)",
              "local=[self.source_profile_jet(packet,4,i,root) for i in range(3)]",
              "{name:IntervalTaylor(c,values) for name,values in atoms['actual_core_atom_axial_coefficients'].items()}"]),
            ('comparison_own_moments',base/'lei_ren_part1_paper_compliant_comparison_point_integrals.py','micro','CompliantComparisonPointIntegrals',
             ["rhs=dict(H=2*uniform_phi,M=uniform_V,K=2*uniform_phi*uniform_V,A=uniform_V*uniform_V,B=uniform_phi*uniform_phi,C=uniform_phi*uniform_phi)",
              "first=r0[name]-rate*initial"]),
            ('exact_comparison_source',base/'lei_ren_part1_paper_compliant_bridge_mixed_C4.py','formal_comparison_source','CompliantBridgeMixedC4',
             ["dict(initial='V_core(4,Z)',integrand='alpha(t)*D_t V_core(4exp(t),Z)',bounds=['0','y'])",
              "'1-sigma((t-hb)/hb)'", "'exp(-lambda*y)'", "'exp(-lambda*(y-t))'",
              "('M',1,'V_bar(t,Z)')", "('A',1,'V_bar(t,Z)^2')"]),
            ('same_core_Cauchy_majorant',base/'lei_ren_part1_paper_compliant_current_original_n1_inner_analytic_domain.py','domain','CurrentOriginalN1InnerAnalyticDomain',
             ["Bpsi=upper(c,s.rebuild.psi_norm)","dict(Phi=Bphi,Psi=Bpsi)","h/4",
              "[upper(c,(4*z+core.j if i==0 else 0)+core.epsilon*4*Bpsi*math.factorial(i)*2**i) for i in range(5)]"]),
        )
        guards={name:_function_ast_guard(path,function,fragments,class_name=cls)
                for name,path,function,cls,fragments in specs}
        for name,path,_,_,_ in specs:
            expected=(self.input_hashes.get(path.name) if path.parent==base else
                      self.implementation_hashes.get(path.relative_to(ROOT).as_posix()))
            if not expected or sha(path)!=expected:raise ValueError('Unbound source AST: '+name)
        family=self.parent['source_family']
        families=[self.first['source_family'],self.domain['source_family'],{k:self.bridge[k] for k in family}]
        checks={'same_family':all(v==family for v in families),
                'pressure_parent_specific_guards':self.parent['source_jet']['source_bindings']['passed'],
                'same_saved_mass':digest(self.mass['interior_mass'])==self.parent['inherited_mass_packet_digest'],
                'original_M_A_ODEs':self.bridge['exact_own_moment_ODEs']['M']['source']=='V' and
                                  self.bridge['exact_own_moment_ODEs']['A']['source']=='V^2',
                'core_seed_not_actual_first_collar_exit':guards['comparison_field_origin']['passed'],
                'domain_majorant_dominates_saved_same_core_Psi_norm':mpf_cmp(
                    tuple(self.domain['joint_core_domain']['Xh_norm_majorants']['Psi']['upper_exact_mpf_tuple']),
                    tuple(self.mixed['source_psi_Xh_norm']['upper_exact_mpf_tuple']))>=0,
                'mixed_same_core_family_and_implicit_source':all(self.mixed[k]==family[k] for k in
                    ('actual_five_defect_family_sha256','implicit_source_sha256')),
                'saved_exact_j_matches_pinned_integer_definition':self.axis_j==Fraction(
                    1,16000000*self.fixed['KN']*self.fixed['CS']*self.fixed['CA']),
                'saved_positive_Z_independent_IBP_weights':self.micro['comparison_moment_enclosure_source']['weights']==
                    'alpha>=0; -alpha_prime>=0; alpha(y)+integral_0^y[-alpha_prime]=1; Z independent',
                'saved_IBP_V_source':self.micro['comparison_moment_enclosure_source']['smoothing_V']==
                    'alpha(y)*V_core(y)+integral_0^y[-alpha_prime(t)]*V_core(t)dt',
                'saved_C5_Taylor_available':all(s['axial_Taylor_order']>=3 and 'k!' in s['units']
                                              for s in self.first['samples'])}
        return {'passed':all(r['passed'] for r in guards.values()) and all(checks.values()),'AST':guards,'checks':checks}

    def source_row(self,ztext,guards):
        from .original_hydro_source_kernel import (
            hydro_drive_ordinary_epsilon_sectors,
            hydro_moment_difference_ordinary_epsilon_sectors,
            evaluate_ordinary_sectors,
        )
        c,iv=self.ctx,self.iv
        sample=next(s for s in self.first['samples'] if s['Z']==ztext and s['phase']=='0')
        point={};truth={}
        for name in ('Psi','Mpsi','Apsi'):
            raw=sample['normalized_fields_and_histories'][name]['inlet'][:4]
            if len(raw)!=4:raise ValueError('Saved ordinary Z3 input is missing')
            pairs=[decode_pair(c,iv,row) for row in raw]
            point[name]=[p*math.factorial(k) for k,(p,t) in enumerate(pairs)]
            truth[name]=[t*math.factorial(k) for k,(p,t) in enumerate(pairs)]
        truth,rho,transport=finite_width_comparison(c,iv,truth,self.domain,self.weight,self.hcap,self.logh,
                                                  source_proved=guards['passed'])
        le,lei=decode_pair(c,iv,sample['source_scales']['log_epsilon'])
        denominator_j=self.axis_j.denominator
        parent_q=self.parent['source_jet']['comparison_q']['correlated_shared_phi_proof']
        qproof=parent_q[('-.3','0','.3').index(ztext)]
        if not qproof['passed'] or not qproof['shared_Phi0_cancellation_is_algebraic']:
            raise ValueError('The saved same-generator quotient proof is required')
        qtruth=[iv.mpf(1 if k==0 else 0)+iv.mpf([-1,1])*iv.mpf(qproof['orders'][k]['representable_absolute_error_cap'])
                for k in range(3)]
        qs=([c.mpf(1),c.mpf(0),c.mpf(0)],qtruth)
        evaluated=[];packets=[]
        for ctx,values,radius,logeps,q in ((c,point,c.mpf(4),le,qs[0]),(iv,truth,rho,lei,qs[1])):
            z=ctx.mpf(ztext);zero=ctx.mpf(0)
            u0=[4*z+ctx.mpf(1)/denominator_j,ctx.mpf(4),zero,zero]
            psi=values['Psi'];mpsi=values['Mpsi'];apsi=values['Apsi']
            u0sq=four_product(u0,u0);u0psi=four_product(u0,psi);psisq=four_product(psi,psi)
            v={0:u0[:3],1:psi[:3]}
            target_m={0:u0,1:psi}
            target_a={0:u0sq,1:[2*x for x in u0psi],2:psisq}
            kappa=ctx.exp(970)
            base=hydro_drive_ordinary_epsilon_sectors(z,None,v,target_m,target_a,delta_sectors={1:kappa})
            dm=[a-b for a,b in zip(mpsi,psi)]
            da1=[2*x for x in four_product(u0,dm)]
            da2=[a-b for a,b in zip(apsi,psisq)]
            change=hydro_moment_difference_ordinary_epsilon_sectors(
                z,None,v,{1:dm},{1:da1,2:da2},delta_sectors={1:kappa})
            # The exact moment-difference mode has no epsilon^0 sector.
            # This follows directly from dm/da, before interval subtraction.
            if 0 in change['numerator_sectors'] and any(x!=0 for x in change['numerator_sectors'][0]):
                raise ValueError('A spurious moment-difference baseline was introduced')
            numerator={}
            for degree,row in base['numerator_sectors'].items():numerator[degree]=[100*x for x in row]
            for degree,row in change['numerator_sectors'].items():
                previous=numerator.get(degree+1,[zero]*3)
                numerator[degree+1]=[a+radius*b for a,b in zip(previous,row)]
            eps=ctx.exp(logeps)
            denominator=evaluate_ordinary_sectors(base['denominator_sectors'],eps)
            contributions={}
            for degree,row in sorted(numerator.items()):
                signed=ordinary_product(q,ordinary_quotient([-v for v in row],denominator))
                contributions[degree]=[v*ctx.exp((degree-1)*logeps) for v in signed]
            result=[sum((row[k] for row in contributions.values()),zero) for k in range(3)]
            evaluated.append(result)
            packets.append((numerator,contributions,base['denominator_sectors'],base['numerator_sectors'],change['numerator_sectors']))
        sp,st=evaluated
        st,nonlinear=macro_nonlinear_enclosure(c,iv,st,self.weight,self.hcap,
                                              macro_norm_bound_proved=guards['passed'])
        np,cp,dp,bp,mp=packets[0];nt,ct,dt,bt,mt=packets[1]
        def encoded(left,right):
            keys=sorted(set(left)|set(right))
            if set(left)!=set(right):raise ValueError('Point/truth formal sectors differ')
            return [{'epsilon_degree':k,'ordinary_rows':[encode_pair(c,iv,c.mpf(p),iv.mpf(t)) for p,t in zip(left[k],right[k])]}
                    for k in keys]
        return {'Z':ztext,'ordinary_Z_orders':[0,1,2],
                'source_rows':[encode_pair(c,iv,p,t) for p,t in zip(sp,st)],
                'source_scope':'same_generator_hydro_sector_q1_factored_macro',
                'full_source_scope':False,'finite_width_comparison':transport,
                'epsilon_factored_numerator':encoded(np,nt),'epsilon_factored_denominator':encoded(dp,dt),
                'epsilon_factored_target_mode':encoded(bp,bt),'epsilon_factored_moment_difference_mode':encoded(mp,mt),
                'normalized_source_contributions':[dict(row,normalized_source_epsilon_power=row['epsilon_degree']-1)
                                                   for row in encoded(cp,ct)],
                'macro_nonlinear_remainder':nonlinear,
                'actual_source_identity':'Bhydro=-(Phi0/barPhi2)/epsilon*(100*Dtarget+R0*Ddifference)*exp(ell)',
                'approximation_identity':'q=1, ell=0, rho2=4 and saved numerical core coefficients; all finite-width changes in truth',
                'expell_factor_enclosed_in_truth':True,
                'mode_origin':'native2 comparison exit; R0=epsilon*rho2, rho2=4*exp(2hb)',
                'epsilon_delta_relation':'delta=epsilon*exp(970)',
                'exact_axis_j_rational':{'numerator':1,'denominator':denominator_j},
                'source_quotient_proof_inherited_digest':digest(qproof),
                'no_formal_sectors_discarded':True,'moment_difference_baseline_zero_by_shared_source':True,
                'approximation_is_exact_source_value':False,'source_width_selected':False}

    def integral_rows(self,rows):
        c,iv=self.ctx,self.iv
        masses={c.mpf(row['phase']):row for row in self.mass['interior_mass']['rows']}
        hydro={row['Z']:row for row in rows}
        pressure={row['Z']:row for row in self.parent['source_jet']['rows']}
        output=[]
        for inherited in self.parent['signed_source_integral']['rows']:
            z,phase=inherited['Z'],inherited['phase']
            sample=next(s for s in self.first['samples'] if s['Z']==z and s['phase']=='0')
            lp,lt=decode_pair(c,iv,sample['source_scales']['log_epsilon'])
            delta,di=c.exp(lp+970),iv.exp(lt+970)
            hp,ht=zip(*(decode_pair(c,iv,row) for row in hydro[z]['source_rows']))
            pp,pt=zip(*(decode_pair(c,iv,row) for row in pressure[z]['source_rows']))
            mass,mi=decode_pair(c,iv,masses[c.mpf(phase)])
            def pair(p,t):
                cp=main_coefficient(c.mpf(z),delta,p,mass,mass_Z_independent=True)
                ct=main_coefficient(iv.mpf(z),di,t,mi,mass_Z_independent=True)
                return [encode_pair(c,iv,a,b) for a,b in zip(cp,ct)]
            output.append({'Z':z,'phase':phase,'ordinary_Z_orders':[0,1],'raw_number':46,
                           'hydro_main_component_C0_C1':pair(hp,ht),
                           'pressure_plus_hydro_main_component_C0_C1':pair([a+b for a,b in zip(hp,pp)],
                                                                          [a+b for a,b in zip(ht,pt)]),
                           'physical_scale':inherited['physical_scale'],'mass_Z_independent':True,
                           'full_source_scope':False,'weighted_raw47_remainder_retained':True})
        return output

    def build(self):
        guards=self.source_guards()
        if guards['passed'] is not True:raise ValueError('Hydro source guards failed: '+str(guards))
        rows=[self.source_row(z,guards) for z in ('-.3','0','.3')]
        integrals=self.integral_rows(rows)
        return {'source_family':self.parent['source_family'],'parent_report':PARENT,'parent_report_sha256':PARENT_SHA256,
                'source_guards':guards,'hydro_source_rows':rows,'signed_source_integral_rows':integrals,
                'inherited_m15_packet':self.parent['inherited_m15_packet'],
                'inherited_m15_packet_digest':self.parent['inherited_m15_packet_digest'],
                'inherited_mass_packet_digest':self.parent['inherited_mass_packet_digest'],
                'inherited_pressure_source_digest':digest(self.parent['source_jet']),
                'inherited_pressure_C1_digest':digest(self.parent['signed_source_integral']),
                'scope':{'hydro_B_ZZ_numerical_approximation_and_enclosure':True,
                         'pressure_plus_hydro_C0_C1':True,'swirl_B_ZZ_recovered':False,
                         'full_source_scope':False,'complete_weighted_primitive':False},
                'contract':{'accepted_ancestor_replayed':False,'mass_recomputed':False,'source_width_selected':False,
                            'positive_cap_as_signed_value':False,'R100_modes_reused_as_R0_modes':False,
                            'raw213_reduced183_retained':True,'all_prior_errors_retained':True},
                'limits':dict(self.parent['limits'],hydro_pressure_is_full_source=False,full_goal_complete=False)}
