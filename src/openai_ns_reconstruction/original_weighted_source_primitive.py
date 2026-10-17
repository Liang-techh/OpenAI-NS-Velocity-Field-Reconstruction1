"""First-switch weighted source queries over pinned, saved original data.

The signed width coefficients are new bounded algebra.  The accepted mass,
core, macro sources and m15 ledger are hydrated without producer replay.
Neither a logarithmic width nor the positive swirl amplitude is materialized.
"""
from __future__ import annotations

import math
from fractions import Fraction

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_hydro_source_point import (
    ROOT, DIRECTORY, FIRST, BRIDGE, DOMAIN, MASS, CORE, sha, digest, load,
    exact_iv, upper, ordinary_norm, width_errors,
)
from .original_leading_viscous_kernel import decode_pair, encode_pair
from .original_modal_source_point import _function_ast_guard
from .original_signed_C1_source import main_coefficient
from .original_source_jet_kernel import ordinary_product, ordinary_quotient
from .original_weighted_angular_kernel import angular_mode_numerator_epsilon_sectors
from .original_weighted_pulse_moments import build_weighted_moments
from .original_weighted_remainder_kernel import signed_remainder_coefficients, enclose_weighted_remainder

PARENT='lei_ren_part1_paper_compliant_current_original_swirl_source_C1_numeric.json'
PARENT_SHA256='b870f9f8957295e9c22a96d077b144a564ffe97278ccec4fe1f13c26dbab14aa'
HYDRO='lei_ren_part1_paper_compliant_current_original_hydro_source_C1_numeric.json'
PRESSURE='lei_ren_part1_paper_compliant_current_original_signed_C1_source_numeric.json'
ANCHORED='lei_ren_part1_paper_compliant_anchored_axis_amplitude.json'
SWITCH='lei_ren_part1_paper_compliant_current_original_first_switch_numeric.json'


def serial(value):
    """Serialize proof arithmetic as exact endpoints, never a selected value."""
    if hasattr(value,'_mpi_'):
        return {'directed_exact_mpf_tuples':[list(v) for v in value._mpi_]}
    if hasattr(value,'_mpf_'):return {'exact_mpf_tuple':list(value._mpf_)}
    if isinstance(value,dict):return {str(k):serial(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serial(v) for v in value]
    return value


def gradient_point(ctx,z,delta,j):
    """Explicit same-source rational G' jet; no Bell-jet cancellation."""
    H=[-4*z**3-j*z*z+(9-delta)*z/2+j,
       -12*z*z-2*j*z+(9-delta)/2,-24*z-2*j]
    L=[1-delta*z*z,-2*delta*z,-2*delta]
    den=ordinary_product(H,H);den[0]+=(j/500)**2
    return ordinary_quotient(ordinary_product(L,H),den)


class OriginalWeightedSourcePrimitive:
    def __init__(self,*,directory=DIRECTORY):
        self.directory=directory
        if sha(directory/PARENT)!=PARENT_SHA256:raise ValueError('Pinned accepted all-sector parent required')
        self.parent=load(directory/PARENT)
        if self.parent['checks']['passed'] is not True:raise ValueError('Passed saved parent required')
        self.input_hashes=dict(self.parent['input_hashes']);self.input_hashes[PARENT]=PARENT_SHA256
        self.implementation_hashes=dict(self.parent['implementation_hashes'])
        for base,group in ((directory,self.input_hashes),(ROOT,self.implementation_hashes)):
            for name,expected in group.items():
                if sha(base/name)!=expected:raise ValueError('Changed accepted pin: '+name)
        required=(FIRST,BRIDGE,DOMAIN,MASS,CORE,HYDRO,PRESSURE,ANCHORED,SWITCH)
        if any(name not in self.input_hashes for name in required):raise ValueError('Missing saved source pin')
        self.first,self.bridge,self.domain,self.mass,self.core,self.hydro,self.pressure,self.anchored,self.switch=[load(directory/n) for n in required]
        if digest(self.mass['interior_mass'])!=self.parent['inherited_mass_packet_digest']:raise ValueError('Saved M0 changed')
        self.ctx=mp.mp.clone();self.ctx.dps=180
        self.iv=MPIntervalContext();self.iv.dps=220
        c,iv=self.ctx,self.iv
        self.weight=exact_iv(c,iv,self.bridge['axial_jet_weight'])
        self.hcap=exact_iv(c,iv,self.bridge['width_enclosure_cap'])
        self.logh=exact_iv(c,iv,self.bridge['source_log_hb_enclosure'])
        fractions={Fraction(v['exact_original_j_fraction']) for v in self.core['levels']}
        if len(fractions)!=1:raise ValueError('Original exact j is inconsistent')
        self.j=fractions.pop()

    def source_guards(self):
        base=ROOT/'experiments/root_st073'
        specs=(
            ('switch_geometry',base/'lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py',
             '__init__','CurrentOriginalPrescribedR100RmC6Functions',
             ["R=g.mul(c(100),g.unary('exp',g.mul(hb,s)))",
              "self.native_parameter('hbS','original switch_owner.flow.h; not merged with hbB by name')",
              "actual=first_switch_fields(a,known,initial,s,hb,R,F0,Pstar,barphi)"]),
            ('first_switch',base/'lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py',
             'first_switch_fields',None,
             ["angular=a.scale(known['D_over_R'],g.mul(g.constant('-1/2'),hb,hb,radius))",
              "logratio=a.integral(a.at(angular,phase,t),t,g.zero,phase,original_first_switch_full_current_radius_direction=True)",
              "chi_axial=g.mul(hb,g.sub(g.one,sigma(g,phase)))",
              "Uy=a.scale(a.mul(ratio,drive),g.neg(chi_axial))",
              "a.integral(a.scale(a.at(Uy,phase,v),hb),v,g.zero,phase,original_first_switch_axial_control=True,source_measure_applied_once=True)"]),
            ('full_drive',base/'lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions.py',
             'full_drive',None,
             ["a.scale(known['drive_hydro'],radius)","a.scale(known['drive_pressure'],g.mul(radius,Pstar,Pstar))",
              "a.mul(a.scale(a.mul(F0,F0),g.mul(radius,radius)),known['drive_swirl'])"]),
            ('angular_direction',base/'lei_ren_part1_paper_compliant_inner_bridge_profiles.py','direction',None,
             ["f=dress(phi,ratios)","h=dress(moments[MTH],ratios)","k=dress(moments[MTHZ],ratios)",
              "(-W+angular/(2*f))/L"]),
            ('angular_gauge',base/'lei_ren_part1_paper_compliant_inner_bridge_profiles.py','inputs','CompliantInnerBridgeProfiles',
             ["gradient=source['L']*H/den", "-multiplier*self.core.Lambda*gradient[j]*out[n-1-j]",
              "ratios(1)"]),
            ('angular_modes',base/'lei_ren_part1_paper_compliant_macro_signed_integrals.py','_mode_rows',None,
             ["[base['D_over_R'],one['D_over_R']-base['D_over_R'],two['D_over_R']-base['D_over_R']]"]),
            ('own_H_K',base/'lei_ren_part1_paper_compliant_comparison_point_integrals.py','micro','CompliantComparisonPointIntegrals',
             ["rhs=dict(H=2*uniform_phi,M=uniform_V,K=2*uniform_phi*uniform_V,A=uniform_V*uniform_V,B=uniform_phi*uniform_phi,C=uniform_phi*uniform_phi)",
              "first=r0[name]-rate*initial"]),
            ('normalized_K_source',ROOT/'src/openai_ns_reconstruction/original_first_collar_jets.py',
             '_refined_increments','OriginalFirstCollarJets',
             ["U0.base+a['Psi']*eps","U0.base*a['H']+a['Kpsi']*eps"]),
        )
        guards={name:_function_ast_guard(path,fn,fragments,class_name=cls) for name,path,fn,cls,fragments in specs}
        for name,path,_,_,_ in specs:
            expected=(self.input_hashes.get(path.name) if path.parent==base else
                      self.implementation_hashes.get(path.relative_to(ROOT).as_posix()))
            if expected!=sha(path):raise ValueError('Unbound actual source: '+name)
        width_identity=self.switch['samples'][0]['background_control_contract']['original_source_width_identity']
        width_gates=('shared_upstream_is_same_object','bridge_logh_is_actual_upstream_object',
                     'switch_width_is_exact_source_reboxing','bridge_hbB_switch_hbS_kept_distinct')
        checks={'accepted_parent_source_guards':self.parent['source_guards']['passed'],
                'same_saved_M0':digest(self.mass['interior_mass'])==self.parent['inherited_mass_packet_digest'],
                'same_actual_flat_cutoff':self.mass['interior_mass']['source']['formula']==
                    'gp(t)=exp(-1/t^2) for t>0; sigma=gp(t)/(gp(t)+gp(1-t))',
                'distinct_width_atoms_same_guarded_log_source':all(width_identity[k] for k in width_gates) and
                    width_identity['shared_K1_ledger_sha256']==self.input_hashes['lei_ren_part1_paper_compliant_K1_ledger.json'],
                'ordinary_gradient_source_bound':self.parent['source_guards']['AST']['anchored_gradient']['passed'],
                'three_actual_sectors_available':len(self.parent['swirl_source_rows'])==len(self.hydro['hydro_source_rows'])==3}
        return {'passed':all(v['passed'] for v in guards.values()) and all(checks.values()),'AST':guards,'checks':checks,
                'inherited_width_identity':width_identity,
                'inherited_cutoff_AST_digest':digest(self.mass['interior_mass']['source_pins']['original_sigma_ast']),
                'source_identity':{'R':'100*exp(hbS*s)',
                    'local_ell':'-hbS^2/2*integral_0^s sum_j A_j*exp((1-j)*hbS*t) dt',
                    'angular_control':'1','axial_control':'1-sigma(s)',
                    'weighted_kernel':'integral_0^s (1-sigma(t))*exp((q-j+1)*hbS*t+local_ell(t)) dt',
                    'hbB_and_hbS_distinct_source_atoms':True,'shared_checked_log_width_only':True,
                    'hbB_held_fixed_in_hbS_Taylor_derivatives':True}}

    def angular_row(self,ztext,guards):
        c,iv=self.ctx,self.iv
        sample=next(s for s in self.first['samples'] if s['Z']==ztext and s['phase']=='0')
        point={};truth={}
        for name in ('Phi','H','Psi','Mpsi','Kpsi'):
            rows=sample['normalized_fields_and_histories'][name]['inlet'][:4]
            if len(rows)!=4:raise ValueError('Missing ordinary Z3 seed')
            pairs=[decode_pair(c,iv,row) for row in rows]
            point[name]=[p*math.factorial(k) for k,(p,t) in enumerate(pairs)]
            truth[name]=[t*math.factorial(k) for k,(p,t) in enumerate(pairs)]
        index=('-.3','0','.3').index(ztext)
        saved_sw=self.parent['swirl_source_rows'][index]
        saved_h=self.hydro['hydro_source_rows'][index]
        inherited={}
        for name,packet in (('Phi',saved_sw),('Psi',saved_h),('Mpsi',saved_h)):
            proof=packet['finite_width_comparison']['ordinary_error_proofs'][name]
            if len(proof)!=4 or not all(r['passed'] for r in proof):raise ValueError('Saved prefix proof missing')
            truth[name]=[v+iv.mpf([-1,1])*iv.mpf(p['representable_absolute_error_cap']) for v,p in zip(truth[name],proof)]
            inherited[name]=digest(proof)
        # H_y+2H=2Phi and Kpsi_y+2Kpsi=2Phi*Psi.  Their positive
        # exponential kernels cover the same native prefix 0<=y<=2hbB.
        phi_cover=ordinary_norm(c,iv,truth['Phi'],self.weight)
        psi_cover=ordinary_norm(c,iv,truth['Psi'],self.weight)
        newproof={}
        for name,forcing in (('H',phi_cover),('Kpsi',phi_cover*psi_cover)):
            bound=4*(ordinary_norm(c,iv,truth[name],self.weight)+forcing)
            errors,proof=width_errors(c,iv,bound,self.weight,self.logh,orders=3,source_proved=guards['passed'])
            truth[name]=[a+b for a,b in zip(truth[name],errors)];newproof[name]=proof
        le,lei=decode_pair(c,iv,sample['source_scales']['log_epsilon'])
        gradient_truth=[exact_iv(c,iv,v)*math.factorial(k) for k,v in enumerate(
            self.anchored['anchored_amplitude_packets'][ztext]['G_gradient_taylor_coefficients'][:3])]
        jpoint=c.mpf(self.j.numerator)/self.j.denominator
        gradient=gradient_point(c,c.mpf(ztext),c.exp(le+970),jpoint)
        rho_proof=saved_h['finite_width_comparison']['rho_error_proof'][0]
        if not rho_proof['passed']:raise ValueError('Saved native radius proof required')
        rho=iv.mpf(4)+iv.mpf([0,1])*iv.mpf(rho_proof['representable_absolute_error_cap'])
        results=[];packets=[]
        for ctx,values,gp,logeps,radius in ((c,point,gradient,le,c.mpf(4)),(iv,truth,gradient_truth,lei,rho)):
            z=ctx.mpf(ztext);j=ctx.mpf(self.j.numerator)/self.j.denominator
            u0=[4*z+j,ctx.mpf(4),ctx.mpf(0),ctx.mpf(0)]
            p=angular_mode_numerator_epsilon_sectors(z,u0,values['Phi'],values['H'],values['Psi'],values['Mpsi'],values['Kpsi'],gp,ctx.exp(970))
            den=[sum((row[k]*ctx.exp(degree*logeps) for degree,row in p['denominator_sectors'].items()),ctx.mpf(0)) for k in range(3)]
            modes=[]
            for mode,sectors in enumerate(p['numerator_modes']):
                num=[sum((ctx.mpf(100)**(1-mode)*radius**mode*ctx.exp((degree+mode)*logeps)*row[k]
                          for degree,row in sectors.items()),ctx.mpf(0)) for k in range(3)]
                modes.append(ordinary_quotient(num,den))
            results.append(modes);packets.append(p)
        encoded=[{'mode_j':j,'radial_rate':1-j,'ordinary_rows':[encode_pair(c,iv,p,t) for p,t in zip(results[0][j],results[1][j])],
                  'Laurent_numerator_degrees':sorted(packets[0]['numerator_modes'][j])} for j in range(3)]
        return results[0],results[1],{'Z':ztext,'ordinary_Z_orders':[0,1,2],'modes':encoded,
            'definition':'A_j=100^(1-j)*rho2^j*sum_e epsilon^(j+e)*N_j,e/(2*Phi*L)',
            'source':'same native2 comparison; rho2=4exp(2hbB)',
            'native2_modes_propagated_once':'H100=Phi2+(R0/100)^2*(H2-Phi2); M100=V2+(R0/100)*(M2-V2)',
            'saved_prefix_proof_digests':inherited,'new_H_Kpsi_prefix_width_proofs':newproof,
            'H_Kpsi_own_ODEs':['H_y+2H=2Phi','Kpsi_y+2Kpsi=2Phi*Psi'],
            'gauge':'g=-Gprime/epsilon; signed rational Gprime jet; no Bell cancellation',
            'anchored_gradient_digest':digest(self.anchored['anchored_amplitude_packets'][ztext]['G_gradient_taylor_coefficients'][:3]),
            'no_Laurent_sector_discarded':True,'source_width_selected':False,'passed':True}

    def build(self):
        from .original_weighted_axial_modes import hydrate_axial_modes
        c,iv=self.ctx,self.iv
        guards=self.source_guards()
        if not guards['passed']:raise ValueError('Actual weighted source guards failed: '+str(guards))
        phases=[row['phase'] for row in self.mass['interior_mass']['rows']]
        weighted=build_weighted_moments(c,iv,phases)
        if not weighted['passed']:raise ValueError('New weighted moments failed')
        masses={c.mpf(row['phase']):row for row in self.mass['interior_mass']['rows']}
        moments={(c.mpf(row['phase']),row['power']):row for row in weighted['rows']}
        angular=[];axial=[];queries=[];byz={}
        for z in ('-.3','0','.3'):
            ap,at,ar=self.angular_row(z,guards);angular.append(ar)
            sample=next(s for s in self.first['samples'] if s['Z']==z and s['phase']=='0')
            hydrated=hydrate_axial_modes(c,iv,z,self.pressure,self.hydro,self.parent,self.bridge,self.weight,self.hcap, sample=sample)
            byz[z]=(ap,at,hydrated,sample)
            regressions=[]
            index=('-.3','0','.3').index(z)
            for label,modes,saved in (
                ('hydro',[m for m in hydrated['modes'] if m['sector']=='hydro'],self.hydro['hydro_source_rows'][index]['source_rows']),
                ('relative_swirl',hydrated['relative_swirl_modes'],self.parent['swirl_source_rows'][index]['relative_source_rows'])):
                for k,ref in enumerate(saved):
                    oldpoint,_=decode_pair(c,iv,ref)
                    newpoint=sum((m['point'][k] for m in modes),c.mpf(0))
                    scaled_difference=abs(newpoint-oldpoint)/max(c.mpf(1),abs(oldpoint))
                    if scaled_difference>c.mpf('1e-75'):raise ValueError('Saved numerical mode regression failed: '+label)
                    regressions.append({'sector':label,'ordinary_Z_order':k,
                        'scaled_approximation_difference':c.nstr(scaled_difference,30),'threshold':'1e-75','passed':True,
                        'numerical_regression_only_not_source_equality_proof':True})
            axial.append(dict(serial(hydrated['evidence']),numerical_mode_regressions=regressions,modes=[
                {'sector':m['sector'],'mode_j':m['mode_j'],'rate':m['rate'],
                 'ordinary_rows':m['serialized_rows'],'identity':m['identity']} for m in hydrated['modes']]))
        for old in self.parent['signed_source_integral_rows']:
            z,phase=old['Z'],old['phase'];ap,at,hydrated,sample=byz[z]
            mass,mi=decode_pair(c,iv,masses[c.mpf(phase)])
            w1,w1i=decode_pair(c,iv,moments[(c.mpf(phase),1)])
            w2,w2i=decode_pair(c,iv,moments[(c.mpf(phase),2)])
            pm=[{'ordinary_rows':m['point'],'rate':m['rate']} for m in hydrated['modes']]
            tm=[{'ordinary_rows':m['truth'],'rate':m['rate']} for m in hydrated['modes']]
            cp=signed_remainder_coefficients(pm,ap,mass,w1,w2,True)
            ct=signed_remainder_coefficients(tm,at,mi,w1i,w2i,True)
            rt,certificate=enclose_weighted_remainder(c,iv,coeffpackettruth=ct,axial_modestruth=tm,
                angular_modestruth=at,W2positiveIV=w2i,weight=self.weight,logh=self.logh,source_proved=True)
            le,lei=decode_pair(c,iv,sample['source_scales']['log_epsilon'])
            delta,di=c.exp(le+970),iv.exp(lei+970)
            pp=main_coefficient(c.mpf(z),delta,cp['leading'],c.mpf(1),mass_Z_independent=True)
            pt=main_coefficient(iv.mpf(z),di,rt,iv.mpf(1),mass_Z_independent=True)
            rel=hydrated['relative_swirl_modes']
            rp=signed_remainder_coefficients([{'ordinary_rows':m['point'],'rate':m['rate']} for m in rel],ap,mass,w1,w2,True)
            ri=signed_remainder_coefficients([{'ordinary_rows':m['truth'],'rate':m['rate']} for m in rel],at,mi,w1i,w2i,True)
            queries.append({'Z':z,'phase':phase,'ordinary_Z_orders':[0,1,2],
                'R_definition':'(sum_j B_j*K_(q-j+1)(s,hbS)-sum_j B_j*M0(s))/hbS',
                'signed_R_ordinary_rows':[encode_pair(c,iv,p,t) for p,t in zip(cp['leading'],rt)],
                'factored_R0_rows':[encode_pair(c,iv,p,t) for p,t in zip(cp['leading'],ct['leading'])],
                'factored_R1_rows':[encode_pair(c,iv,p,t) for p,t in zip(cp['next'],ct['next'])],
                'relative_swirl_R0_rows':[encode_pair(c,iv,p,t) for p,t in zip(rp['leading'],ri['leading'])],
                'relative_swirl_R1_rows':[encode_pair(c,iv,p,t) for p,t in zip(rp['next'],ri['next'])],
                'swirl_absolute_atom':hydrated['absolute_swirl_atom'],
                'raw47_signed_C0_C1':[encode_pair(c,iv,p,t) for p,t in zip(pp,pt)],
                'physical_scale':'R100*epsilon*hbS^3; R100=100',
                'raw46_main_inherited_digest':digest(old['all_sector_main_component_C0_C1']),
                'weighted_J_factored_identity':'J=sum_j B_j*M0+hbS*R; R-weighted integral=100*J',
                'finite_width_certificate':serial(certificate),
                'M0_hydrated_not_recomputed':True,'W1_W2_independent_of_M0':True,
                'old_raw47_error_budget_replaced':False,'complete_functional_closure':False})
        return {'source_family':self.parent['source_family'],'parent_report':PARENT,'parent_report_sha256':PARENT_SHA256,
            'source_guards':guards,'weighted_pulse_moments':weighted,'angular_source_rows':angular,
            'axial_mode_evidence':axial,'weighted_remainder_queries':queries,
            'inherited_m15_packet':self.parent['inherited_m15_packet'],
            'inherited_m15_packet_digest':self.parent['inherited_m15_packet_digest'],
            'inherited_mass_packet_digest':self.parent['inherited_mass_packet_digest'],
            'inherited_all_sector_source_digest':digest(self.parent['all_sector_source_rows']),
            'inherited_main_integrals_digest':digest(self.parent['signed_source_integral_rows']),
            'contract':{'accepted_ancestor_replayed':False,'M0_recomputed':False,
                'raw213_reduced183_and_all_errors_retained':True,'old_raw47_budget_replaced':False,
                'positive_cap_selected_as_value':False,'hb_or_absolute_swirl_amplitude_materialized':False,
                'local_width_expansion_is_temporal_n_recursion':False},
            'scope':{'weighted_source_enclosed_at_19_queries':True,'ordinary_R_through_Z2':True,
                'raw47_C_operator_through_ordinary_Z1':True,'new_signed_R0_R1_coefficients':True},
            'limits':dict(self.parent['limits'],full_weighted_functional_primitive=False,
                raw47_ledger_replacement=False,global_five_moment_closure=False,full_goal_complete=False)}
