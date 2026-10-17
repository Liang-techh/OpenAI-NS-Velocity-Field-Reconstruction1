"""Signed macro weighted source primitive over saved ordinary source data.

The native lower endpoint is a separate query, not the macro source port.
No macro endpoint error is replaced or cancelled. The accepted signed183
views stay exact. No accepted owner or numerical quadrature is replayed.
"""
from __future__ import annotations

import copy
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_hydro_source_point import ROOT,DIRECTORY,sha,digest,load,exact_iv
from .original_leading_viscous_kernel import decode_pair,encode_pair,endpoint_coefficient,mean_source_coefficient
from .original_modal_source_point import _function_ast_guard

PARENT='lei_ren_part1_paper_compliant_current_original_weighted_leaf_m15_numeric.json'
PARENT_SHA256='3076cf02cfd9a10ce12021eeafd3549db5e43f71659eed9fcdf9abe5476a33ac'
WEIGHTED='lei_ren_part1_paper_compliant_current_original_weighted_source_primitive_numeric.json'
MASS='lei_ren_part1_paper_compliant_current_original_source_jet_mass_numeric.json'
FIRST='lei_ren_part1_paper_compliant_current_original_first_collar_jets.json'
BRIDGE='lei_ren_part1_paper_compliant_actual_bridge_integrals.json'
PRESSURE='lei_ren_part1_paper_compliant_current_original_signed_C1_source_numeric.json'


class OriginalMacroWeightedM15:
    def __init__(self,*,directory=DIRECTORY):
        self.directory=Path(directory)
        if sha(self.directory/PARENT)!=PARENT_SHA256:raise ValueError('Accepted signed weighted views required')
        self.parent=load(self.directory/PARENT)
        if not self.parent['checks']['passed']:raise ValueError('Passed saved parent required')
        self.input_hashes=dict(self.parent['input_hashes']);self.input_hashes[PARENT]=PARENT_SHA256
        self.implementation_hashes=dict(self.parent['implementation_hashes'])
        for base,group in ((self.directory,self.input_hashes),(ROOT,self.implementation_hashes)):
            for name,expected in group.items():
                if sha(base/name)!=expected:raise ValueError('Changed accepted source pin: '+name)
        names=(WEIGHTED,MASS,FIRST,BRIDGE,self.parent['inherited_m15_packet']['report'],PRESSURE)
        if any(n not in self.input_hashes for n in names):raise ValueError('Missing same-source data pin')
        self.weighted,self.mass,self.first,self.bridge,self.ledger,self.pressure=[load(self.directory/n) for n in names]
        self.ctx=mp.mp.clone();self.ctx.dps=180
        self.iv=MPIntervalContext();self.iv.dps=220

    def source_guards(self):
        src=ROOT/'src/openai_ns_reconstruction';base=ROOT/'experiments/root_st073'
        specs=(
            ('proxy_mean',src/'original_macro_m15.py','proxy_modes',None,
             ["qa=-(a*(z*(1-delta))+derivative(a)*D)/A"]),
            ('actual_macro_endpoints',src/'original_macro_m15.py','leading','OriginalMacroM15',
             ["dy=[inputs['source_Uy']['truth'].truncate(2),lower_dy]",
              "lower_dy=(lower['actual']['Phi']/lower['comparison']['Phi'][0]*lower['comparison']['Psi'][1]*lower['rho']).truncate(2)"]),
            ('original_macro_leaf_scales',src/'original_macro_m15.py','evaluate','OriginalMacroM15',
             ["scale_atom=physical(1,1) if i==0 else o._summary(c.log(4)-2*o.core.logLambda,iv.log(4)-2*o.core.ll_iv)",
              "append('macro_leading_viscous_'+name+'_error',zero,leading['viscous_error_per_epsilon_hb'][i],[scale_atom,h],dict(binding,error_enclosure_only=True,endpoint=name,proxy_R_Qy_increment_zero_by_exact_exponential_identity=True,actual_Uy_and_nonzero_mean_endpoint_retained=True))"]),
            ('uniform_macro_angular',base/'lei_ren_part1_paper_compliant_actual_bridge_integrals.py','prepare','CompliantActualBridgeIntegrals',
             ["ell_box=self.scaled_jet_box(self.symmetric_norm_row(B),1,norm_only=True)",
              "dict(data=data,Z=z,inputs=inputs,micros=micros,end=end,modes=modes,R0=R0,L=L,quotient=data['phi0'].truncate(5)/phi2.truncate(5),angular_prefix_norm_per_hb=B,uniform_ell_box=ell_box)"]),
            ('uniform_angular_source',base/'lei_ren_part1_paper_compliant_actual_bridge_integrals.py','contributions','CompliantActualBridgeIntegrals',
             ["exp_bound=c.exp(6*self.cap)","angular_norm=p['angular_prefix_norm_per_hb']"]),
        )
        guards={name:_function_ast_guard(path,fn,fragments,class_name=cls) for name,path,fn,cls,fragments in specs}
        for name,path,_,_,_ in specs:
            expected=self.input_hashes.get(path.name) if path.parent==base else self.implementation_hashes.get(path.relative_to(ROOT).as_posix())
            if expected!=sha(path):raise ValueError('Unbound macro source AST: '+name)
        lower=self.mass['transported_lower_dy']['source_bindings']
        proof_formula='norm_order2(B)*6*cap*exp(6*cap); ordinary error k!/weight^k'
        proofs=[]
        for row in self.weighted['axial_mode_evidence']:
            proofs.extend(row['hydro_mode_proofs'].values())
            proofs.extend(v['macro_proof'] for v in row['swirl_mode_proofs'])
        proofs.extend(v['macro_nonlinear_remainder'] for v in self.pressure['source_jet']['rows'])
        proof_fields=('passed','macro_ell_source_bound_used','first_switch_ell_not_substituted','point_unchanged')
        checks={'same_family':self.parent['source_family']==self.weighted['source_family']==self.mass['source_family'],
                'inherited_macro1_first0_source_and_FTC':self.parent['source_guards']['passed'],
                'transported_lower_source':lower['passed'] and lower['modified_Qy_not_used'] and lower['graph']['source_handles_passed'],
                'rho2_ordinary_Z_derivatives_zero':all(v['rho']['ordinary_Z_derivatives']==['0','0'] and v['rho']['epsilon_and_hb_factored'] for v in self.mass['transported_lower_dy']['rows']),
                'ordinary_Z2_lower':all(v['ordinary_Z_orders']==[0,1,2] for v in self.mass['transported_lower_dy']['rows']),
                'new_native2_modes_have_uniform_macro_enclosure':all(v['passed'] and v['all_original_sectors_retained'] for v in self.weighted['axial_mode_evidence']),
                'uniform_prefix_proof_fields':len(proofs)==18 and all(v['formula']==proof_formula and all(v[k] is True for k in proof_fields) for v in proofs),
                'six_modes_rates_and_ordinary_Z2':all(v['ordinary_Z_orders']==[0,1,2] and [m['rate'] for m in v['modes']]==[2,2,1,3,2,1] and all(len(m['ordinary_rows'])==3 for m in v['modes']) for v in self.weighted['axial_mode_evidence']),
                'distinct_original_width_geometry':self.weighted['source_guards']['checks']['distinct_width_atoms_same_guarded_log_source']}
        return {'passed':all(v['passed'] for v in guards.values()) and all(checks.values()),'AST':guards,'checks':checks,
            'inherited_actual_FTC_and_seam_digest':digest(self.parent['source_guards']),
            'inherited_lower_source_digest':digest(lower),'uniform_prefix_proofs_digest':digest(proofs),'uniform_prefix_proof_count':len(proofs),
            'proxy_identity':'Rproxy*Qy_proxy=-4*epsilon^2*C_delta(a0); Rproxy=4epsilon*exp(Yq), a_proxy=a0*exp(-Yq)',
            'macro_FTC_identity':'Hmacro0-Hmacro1=Emacro0-Emacro1+integral_macro R*C_delta(Uy)dy',
            'native_lower_is_macro_source_port_proved':False,
            'native_lower_and_macro0_value_continuity_proved':False,
            'coupled_raw23_raw24_raw44_replacement_applied':False,
            'unresolved_seam':'Native lower radial dy and macro source_Uy are distinct ports. A separate value-continuity/jump identity is required before substituting the native endpoint into the macro FTC.',
            'uniform_enclosure_semantics':'C_j100*exp(ell_macro(y)) enclosed for every macro prefix; this is not a frozen endpoint exp(ell100) selection',
            'source_caps_or_interval_overlap_used_as_function_identity':False}

    def evaluate_sources(self,z,guards):
        from .original_macro_weighted_kernel import evaluate_macro_modes
        c,iv=self.ctx,self.iv
        source=next(r for r in self.weighted['axial_mode_evidence'] if r['Z']==z)
        lower=next(r for r in self.mass['transported_lower_dy']['rows'] if r['Z']==z)
        sample=next(r for r in self.first['samples'] if r['Z']==z and r['phase']=='0')
        le,lei=decode_pair(c,iv,sample['source_scales']['log_epsilon'])
        # The saved radius error only widens rho2; it never supplies hbB.
        cap=c.make_mpf(tuple(lower['rho']['width_cap']['upper_exact_mpf_tuple']))
        rho=(c.mpf(4),iv.mpf(4)+iv.mpf([0,cap]))
        lengths=(c.log(100)-le-c.log(rho[0]),iv.log(100)-lei-iv.log(rho[1]))
        modes=[{'rate':v['rate'],'pairs':[decode_pair(c,iv,r) for r in v['ordinary_rows']]} for v in source['modes']]
        inputs=[decode_pair(c,iv,r) for r in lower['result']]
        results=[]
        for index,ctx in enumerate((c,iv)):
            packet=evaluate_macro_modes(ctx,[{'rate':m['rate'],'ordinary_rows':[r[index] for r in m['pairs']]} for m in modes],
                lengths[index],source_proved=guards['passed'],L_Z_independent=True)
            de=ctx.exp((le if index==0 else lei)+970)
            response=packet['response'];C=mean_source_coefficient(ctx.mpf(z),de,response)
            E=[v*rho[index]/4 for v in endpoint_coefficient(ctx.mpf(z),de,[r[index] for r in inputs])]
            results.append({'weighted_response':response,'weighted_C':C,'lower_E':E,'kernels':packet['kernels']})
        def encoded(name):return [encode_pair(c,iv,p,t) for p,t in zip(results[0][name],results[1][name])]
        return {'Z':z,'ordinary_source_orders':[0,1,2],'ordinary_coefficient_orders':[0,1],
            'weighted_response_rows':encoded('weighted_response'),'weighted_C0_C1':encoded('weighted_C'),
            'native_lower_endpoint_C0_C1':encoded('lower_E'),
            'kernel_rows':[encode_pair(c,iv,p,t) for p,t in zip(results[0]['kernels'],results[1]['kernels'])],
            'mode_rates':[m['rate'] for m in modes],
            'physical_integral_scale':'R100*epsilon*hbB','physical_lower_scale':'4*epsilon^2*hbB',
            'length':encode_pair(c,iv,lengths[0],lengths[1]),'length_definition':'log(100/(epsilon*rho2)); rho2=4exp(2hbB)',
            'rho2':encode_pair(c,iv,rho[0],rho[1]),'rho_inside_dy_and_external_R0_have_distinct_roles':True,
            'lower_dy_source_digest':digest(lower),'modal_source_digest':digest(source),
            'native_endpoint_is_macro0_endpoint':False,'native_endpoint_used_in_ledger':False,
            'weighted_primitive_function':'sum_j integral_0^L Cj100(Z)*exp(rate_j*(u-L))*exp(ell_macro(u,Z))du',
            'point_approximation':'sum_j Cj100_point(Z)*K_rate_j(Lpoint); ell=0 and rho2=4 only in approximation',
            'uniform_macro_angular_enclosure_not_frozen_endpoint_value':True,
            'kernel_positive_weights_and_L_Z_independent':True,'old_mass_or_core_computed':False}

    def build(self):
        guards=self.source_guards()
        if not guards['passed']:raise ValueError('Macro primitive source guards failed: '+str(guards))
        sources=[self.evaluate_sources(z,guards) for z in ('-.3','0','.3')]
        # The missing native2/macro0 source identity is not supplied by the
        # interval enclosures. Keep all 183 accepted signed entries untouched.
        views=copy.deepcopy(self.parent['signed_m15_query_views'])
        seam=self.ledger['source_seam_proofs']['joins']['micro_macro']
        if not seam['passed'] or not seam['leading_Qy_common']['passed']:raise ValueError('Saved physical seam witness failed')
        partial_seam={'saved_physical_seam_passed':True,'saved_physical_seam_digest':digest(seam),
            'left_scope':copy.deepcopy(seam['left_scope']),'right_scope':copy.deepcopy(seam['right_scope']),
            'common_Qy_operator':copy.deepcopy(seam['leading_Qy_common']),
            'explicit_normalized_Psi_Mpsi_ledger_theorem_bound_here':False,
            'raw23_raw24_signed_endpoint_ledger_binding_bound_here':False,
            'source_values_selected':False,
            'next_source_AST_locations':['original_macro_numeric.py:_macro_increments,_background',
                'original_second_micro_numeric.py:_second_increments',
                'lei_ren_part1_paper_compliant_actual_bridge_integrals.py:contributions',
                'original_macro_pressure_rows.py:inputs,background'],
            'partial_source_observation':'Macroq0 prefix has first/second phases1 and zero macro integral; mean theta0=1,mass0=0. Formalize the normalized coefficient seam and signed endpoint ledger binding before contraction.'}

        return {'source_family':self.parent['source_family'],'parent_report':PARENT,'parent_report_sha256':PARENT_SHA256,
            'source_guards':guards,'macro_weighted_source_rows':sources,'signed_m15_query_views':views,
            'inherited_m15_packet':self.parent['inherited_m15_packet'],
            'inherited_m15_packet_digest':self.parent['inherited_m15_packet_digest'],
            'inherited_signed183_views_digest':digest(self.parent['signed_m15_query_views']),
            'source_primitive_only':True,'signed183_views_unchanged':True,
            'native_macro0_partial_seam_evidence':partial_seam,
            'endpoint_coupling':{'applied':False,'raw23_raw24_raw44_retained':True,
                'native_lower_equals_macro0_Uy_proved':False,
                'native_lower_macro0_value_and_jump_identity_proved':False,
                'required_identity':'If Psi and Mpsi continue exactly, Hnative-Hmacro0=Enative-Emacro0. This conditional identity has not been applied.',
                'no_interval_overlap_or_approximation_agreement_used_as_proof':True},
            'contract':{'accepted_ancestor_or_fixture_replayed':False,'mass_core_or_old_integrals_recomputed':False,
                'raw213_reduced183_archives_unchanged':True,'signed_raw46_raw47_unchanged':True,
                'six_inner_lower_reserved_once_target_and_N_unchanged':True,
                'positive_caps_midpoints_or_interval_overlap_used_as_values':False,
                'hbB_hbS_or_absolute_swirl_amplitude_materialized':False,
                'raw23_raw24_raw43_raw44_raw45_and_all_other_errors_retained':True},
            'limits':{'numeric_global_axial_function_primitive':False,'numeric_Rp_terminal_join':False,
                'global_five_moment_common_N_closure':False,'exact_heat_energy_stress_flat':False,
                'true_temporal_n_recursion':False,'mean_oscillatory_correction':False,
                'independent_corrected_Cartesian_NS':False,'full_goal_complete':False}}
