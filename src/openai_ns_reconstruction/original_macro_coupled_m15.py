"""Correlate actual macro leading endpoints before numerical enclosure.

Saved signed primitives are hydrated, never recomputed. Native and macro
radial derivatives remain distinct; their value-continuity seam removes
the common mean part of the endpoint difference as a function identity.
"""
from __future__ import annotations

import copy
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from .original_hydro_source_point import ROOT,DIRECTORY,sha,digest,load
from .original_leading_viscous_kernel import decode_pair,encode_pair,endpoint_coefficient
from .original_modal_source_point import _function_ast_guard

PARENT='lei_ren_part1_paper_compliant_current_original_macro_weighted_m15_numeric.json'
PARENT_SHA256='6835984804c4db15884b1832b4d95b16ccf1e394b2fd9145f612509164679cb7'
MASS='lei_ren_part1_paper_compliant_current_original_source_jet_mass_numeric.json'
MACRO='lei_ren_part1_paper_compliant_current_original_macro_pressure_rows.json'
LEADING='lei_ren_part1_paper_compliant_current_original_leading_viscous_source_numeric.json'


def materialize_coupled_view(original_terms,entries):
    """Expand the contracted view against immutable original reduced183."""
    if len(original_terms)!=183 or len(entries)!=182:raise ValueError('Expected original183 and coupled182 slots')
    if [v['parent_reduced_index'] for v in entries]!=[i for i in range(183) if i!=39]:raise ValueError('Wrong contracted endpoint')
    out=[]
    for row in entries:
        old=original_terms[row['parent_reduced_index']]
        if digest(old)!=row['parent_term_digest']:raise ValueError('Original source leaf changed')
        if row['operation']=='unchanged_parent_reference':out.append(copy.deepcopy(old))
        elif row['operation']=='same_source_signed_leaf_replacement':out.append(copy.deepcopy(row['term']))
        else:raise ValueError('Unknown source operation')
    return out


class OriginalMacroCoupledM15:
    def __init__(self,*,directory=DIRECTORY):
        self.directory=Path(directory)
        if sha(self.directory/PARENT)!=PARENT_SHA256:raise ValueError('Accepted macro primitive required')
        self.parent=load(self.directory/PARENT)
        if not self.parent['checks']['passed']:raise ValueError('Passed saved parent required')
        if not self.parent['source_guards']['checks']['rho2_ordinary_Z_derivatives_zero']:raise ValueError('Z-independent radius required')
        self.input_hashes=dict(self.parent['input_hashes']);self.input_hashes[PARENT]=PARENT_SHA256
        self.implementation_hashes=dict(self.parent['implementation_hashes'])
        for base,group in ((self.directory,self.input_hashes),(ROOT,self.implementation_hashes)):
            for name,expected in group.items():
                if sha(base/name)!=expected:raise ValueError('Changed accepted pin: '+name)
        names=(MASS,MACRO,self.parent['inherited_m15_packet']['report'],LEADING)
        if any(n not in self.input_hashes for n in names):raise ValueError('Missing source pin')
        self.mass,self.macro,self.ledger,self.ftcreport=[load(self.directory/n) for n in names]
        self.ctx=mp.mp.clone();self.ctx.dps=180
        self.iv=MPIntervalContext();self.iv.dps=220

    def ledger_binding(self,seam):
        """Bind the whole signed endpoint increment, never a value of a cap."""
        path=ROOT/'src/openai_ns_reconstruction/original_macro_m15.py'
        fragments=[
            "dy=[inputs['source_Uy']['truth'].truncate(2),lower_dy]",
            "ey=(source*(2*z)-((epsi-emean)*(z*(1-o.core.delta_iv))+o._derivative(epsi-emean)*D))/A",
            "gr=2*(1-t)*iv.exp(2*o.hcap*(1-t))",
            "viscous.append(self.symmetric(ey.truncate(1))*(1+o.hcap*gr)+self.symmetric((-ta*decay).truncate(1))*gr)",
        ]
        guards={'macro_endpoints':_function_ast_guard(path,'leading',fragments,class_name='OriginalMacroM15')}
        specs=[(path,'proxy_modes',None,["qa=-(a*(z*(1-delta))+derivative(a)*D)/A"]),
            (path,'kernels',None,["x=Y*q","E=ctx.exp(-x)"]),
            (path,'evaluate','OriginalMacroM15',[
                "scale_atom=physical(1,1) if i==0 else o._summary(c.log(4)-2*o.core.logLambda,iv.log(4)-2*o.core.ll_iv)",
                "append('macro_leading_viscous_'+name+'_error',zero,leading['viscous_error_per_epsilon_hb'][i],[scale_atom,h],dict(binding,error_enclosure_only=True,endpoint=name,proxy_R_Qy_increment_zero_by_exact_exponential_identity=True,actual_Uy_and_nonzero_mean_endpoint_retained=True))"]),
            (ROOT/'src/openai_ns_reconstruction/original_first_switch_m15.py','leading','OriginalFirstSwitchM15',[
                "low=macro['source_Uy']",
                "endpoints=[(coefficient(c,up['point'],o.core.delta)*(-1),coefficient(iv,up['truth'],o.core.delta_iv)*(-1)),(coefficient(c,low['point'],o.core.delta),coefficient(iv,low['truth'],o.core.delta_iv))]"]),
            (ROOT/'src/openai_ns_reconstruction/original_first_switch_m15.py','evaluate','OriginalFirstSwitchM15',[
                "append('first_switch_viscous_'+endpoint,p,t,[physical(1,1),hbS if k==0 else hbB],actual_Uy_source_endpoint=True,endpoint=endpoint)"])]
        for n,(source,fn,cls,parts) in enumerate(specs):
            if self.implementation_hashes.get(source.relative_to(ROOT).as_posix())!=sha(source):raise ValueError('Unpinned ledger source')
            guards[fn+'_'+str(n)]=_function_ast_guard(source,fn,parts,class_name=cls)
        if self.implementation_hashes.get(path.relative_to(ROOT).as_posix())!=sha(path):raise ValueError('Unpinned macro source')
        ftc=self.ftcreport['independent_leading_FTC']
        density=ROOT/ftc['physical_density_FTC_source']
        if self.input_hashes.get(density.name)!=sha(density):raise ValueError('Unpinned actual physical density')
        guards['physical_m15_density']=_function_ast_guard(density,'pressure_density',["a.scale(Py,g.neg(R))"])
        guards['physical_RQy_FTC']=_function_ast_guard(density,'pressure_FTC',["a.scale(a.add(Qy,Qyy),R)"])
        ftc_ok=(ftc['passed'] and ftc['physical_R_once'] and ftc['physical_Pstar_once'] and ftc['ordinary_Z1']
            and ftc['absolute_mean_retained_in_H'] and ftc['first1_H_not_zeroed']
            and ftc['mean_ODE_partition']=='Hlower-Hupper=Elower-Eupper+integral R*C(Uy) dy'
            and ftc['endpoint_operator']=='E=2Z*R*Uy/L'
            and ftc['source_operator']=='C(Uy)=[(1-delta)Z*Uy+(1-Z^2)*Uy_Z]/L'
            and all(v['modified_Qy_alias'] is False and v['ordinary_Z_orders']==[0,1] for v in ftc['endpoint_recipes'].values())
            and all(v['passed'] and len(v['mean_ODE_rows'])==3 and all(r['passed'] for r in v['mean_ODE_rows']) for v in ftc['independent_leading_charts'].values())
            and self.parent['source_guards']['passed']
            and all(v['ordinary_source_orders']==[0,1,2] and v['ordinary_coefficient_orders']==[0,1]
                and v['mode_rates']==[2,2,1,3,2,1] and v['kernel_positive_weights_and_L_Z_independent'] for v in self.parent['macro_weighted_source_rows']))
        return {'passed':all(v['passed'] for v in guards.values()) and seam['passed'] and ftc_ok,'AST':guards,
            'independent_leading_FTC_digest':digest(ftc),'independent_leading_FTC_passed':ftc_ok,
            'physical_density_FTC_source':ftc['physical_density_FTC_source'],
            'mean_ODE_partition':ftc['mean_ODE_partition'],'absolute_mean_retained_in_H':True,
            'saved_macro_prefix_source_proof_digest':digest(self.parent['source_guards']),
            'signed_endpoint_orientation':{'native_lower':1,'macro_upper':-1,'first_switch_lower_E':1},
            'actual_minus_proxy_category_reconstructed_from_original_physical_FTC':True,
            'old_upper_lower_enclosure_approximations_not_read_as_source_values':True,
            'old_raw_indices':[23,24,44],'old_reduced_indices':[19,20,39],
            'whole_increment':'(Hnative-Hmacro1)-(Hproxy0-Hproxy1)+Emacro1',
            'proxy_increment_zero':'Rproxy*Qy_proxy is independent of q',
            'value_seam':'Hnative-Hmacro0=Enative-Emacro0',
            'macro_FTC':'Hmacro0-Hmacro1=Emacro0-Emacro1+I_macro',
            'new_signed_function':'Enative+I_macro',
            'individual_positive_caps_used_as_signed_values':False,
            'actual_lower_radius_rho2_over4_retained':True,
            'native_macro_derivative_ports_equated':False,'actual_radius_ordinary_Z_derivatives_zero':True}

    def defect_rows(self,seam):
        from .original_macro_defect_kernel import evaluate_seam_defect
        c,iv=self.ctx,self.iv;rows=[]
        for saved in self.parent['macro_weighted_source_rows']:
            z=saved['Z']
            lower=next(v for v in self.mass['transported_lower_dy']['rows'] if v['Z']==z)
            macro=next(v for v in self.macro['samples'] if v['Z']==z and v['coordinate']['value']=='0')
            bp=[decode_pair(c,iv,v) for v in lower['result'][:2]]
            mpairs=[decode_pair(c,iv,v) for v in macro['leading_U_y_factored']['coefficient_ordinary_Z']]
            if len(mpairs)!=2:raise ValueError('Macro0 ordinary source Z1 required')
            rho=decode_pair(c,iv,saved['rho2']);outputs=[]
            # delta is the same original fixed parameter recovered from epsilon.
            first_scale=macro['leading_U_y_factored']['log_scale_terms'][0]
            le=decode_pair(c,iv,first_scale)
            for index,ctx in enumerate((c,iv)):
                delta=ctx.exp(le[index]+970)
                result=evaluate_seam_defect(ctx,ctx.mpf(z),delta,[v[index] for v in bp],[v[index] for v in mpairs],
                    values_source_proved=seam['passed'],same_radius_proved=seam['passed'])
                scale=rho[index]/4
                outputs.append({'defect':[v*scale for v in result['defect_C0_C1']],
                    'macro_E0':[v*scale for v in endpoint_coefficient(ctx.mpf(z),delta,[v[index] for v in mpairs])]})
            def encoded(name):return [encode_pair(c,iv,p,t) for p,t in zip(outputs[0][name],outputs[1][name])]
            rows.append({'Z':z,'ordinary_Z_orders':[0,1],
                'defect_C0_C1':encoded('defect'),'macro0_endpoint_C0_C1':encoded('macro_E0'),
                'physical_scale':'4*epsilon^2*hbB','external_actual_radius_factor':'rho2/4',
                'native_endpoint_C0_C1':copy.deepcopy(saved['native_lower_endpoint_C0_C1']),
                'correlated_defect_plus_macro0_endpoint':'Enative; same source algebra before interval arithmetic',
                'correlated_result_C0_C1':copy.deepcopy(saved['native_lower_endpoint_C0_C1']),
                'native_source_digest':digest(lower),'macro0_source_digest':digest(macro),
                'source_ports_distinct':True,'independent_interval_sum_used_as_identity':False,
                'defect_added_as_an_extra_ledger_term':False})
        return rows

    def build(self):
        from .original_native_macro_seam import build_native_macro_seam
        seam=build_native_macro_seam(self.parent,root=ROOT,directory=self.directory)
        if not seam['passed']:raise ValueError('Normalized native/macro seam failed')
        binding=self.ledger_binding(seam)
        if not binding['passed']:raise ValueError('Signed macro ledger binding failed')
        defects=self.defect_rows(seam);views=[]
        for before in self.parent['signed_m15_query_views']:
            old=next(v for v in self.ledger['samples'] if v['Z']==before['Z'])
            source=next(v for v in self.parent['macro_weighted_source_rows'] if v['Z']==before['Z'])
            entries=[];removed=[]
            for entry in before['entries']:
                index=entry['parent_reduced_index']
                if index==39:removed.append(dict(copy.deepcopy(entry),raw_index=44,active=False,provenance_only=True,no_cap_value_or_midpoint=True));continue
                new=copy.deepcopy(entry)
                if index in (19,20):
                    term=copy.deepcopy(old['correlated_reduced_m15_C0_C1_terms'][index])
                    term['name']='macro_first_leading_weighted_source_integral' if index==19 else 'macro_leading_native2_source_endpoint'
                    term['coefficient_C0_C1']=copy.deepcopy(source['weighted_C0_C1'] if index==19 else source['native_lower_endpoint_C0_C1'])
                    archived_binding=term['actual_function_binding']
                    term['actual_function_binding']=dict(error_enclosure_only=False,
                        source_function='I_macro' if index==19 else 'E_native_lower',
                        physical_scale_definition='100*epsilon*hbB' if index==19 else '4*epsilon^2*hbB; rho2/4 retained inside coefficient',
                        interval='whole original macro0..1 plus first-switch lower E source port',
                        actual_signed_coupled_macro_source=True,native_macro_value_seam_digest=digest(seam),
                        whole_signed_endpoint_binding_digest=digest(binding),old_raw_indices=[23,24,44],
                        archived_error_enclosure_binding_digest=digest(archived_binding),
                        signed_function_source_digest=digest(source),individual_cap_selected_as_value=False,
                        source_errors_and_actual_geometry_in_signed_primitive=True,
                        proxy_increment_cancelled_by_function_identity=True)
                    new.update(operation='same_source_signed_leaf_replacement',term=term,coupled_raw_indices=[23,24,44])
                entries.append(new)
            expanded=materialize_coupled_view(old['correlated_reduced_m15_C0_C1_terms'],entries)
            views.append({'Z':before['Z'],'parent_original_reduced183_digest':before['archived_reduced183_digest'],
                'raw213_digest':before['raw213_digest'],'parent_signed183_view_digest':before['expanded_signed_view_digest'],
                'entries':entries,'view_terms':182,'unchanged_original_references':178,
                'four_signed_replacement_slots':[19,20,41,42],'contracted_parent_slot':39,
                'contracted_source_endpoint_provenance':removed,'expanded_signed182_digest':digest(expanded),
                'retained_first_switch_signed_pair_digest':digest([v for v in entries if v['parent_reduced_index'] in (41,42)]),
                'independent_target':copy.deepcopy(before['independent_target']),
                'independent_target_digest':before['independent_target_digest'],
                'selected_N_evidence_digest':before['selected_N_evidence_digest'],
                'six_inner_lower_digest':before['six_inner_lower_digest'],
                'largest_unchanged_caps':copy.deepcopy(before['largest_unchanged_caps']),
                'signed_discrepancy_recipe':'sum(signed182 physical leaves)-unchanged independent negative complete-pulse m15 target',
                'global_terminal_closure':False})
        return {'source_family':self.parent['source_family'],'parent_report':PARENT,'parent_report_sha256':PARENT_SHA256,
            'native_macro_seam':seam,'signed_ledger_binding':binding,'native_macro_defect_rows':defects,
            'macro_weighted_source_rows':copy.deepcopy(self.parent['macro_weighted_source_rows']),
            'macro_weighted_source_rows_digest':digest(self.parent['macro_weighted_source_rows']),
            'signed_m15_query_views':views,'inherited_m15_packet':self.parent['inherited_m15_packet'],
            'inherited_m15_packet_digest':self.parent['inherited_m15_packet_digest'],
            'inherited_signed183_views_digest':digest(self.parent['signed_m15_query_views']),
            'contract':{'accepted_ancestor_or_fixture_replayed':False,'saved_macro_primitive_recomputed':False,
                'mass_core_M0_W1_W2_recomputed':False,'raw213_reduced183_archives_unchanged':True,
                'signed_raw46_raw47_unchanged':True,'six_inner_lower_reserved_once_target_and_N_unchanged':True,
                'positive_caps_midpoints_or_interval_overlap_used_as_values':False,
                'hbB_hbS_or_absolute_swirl_amplitude_materialized':False,
                'raw43_raw45_and_all_other_errors_retained':True},
            'limits':{'numeric_global_axial_function_primitive':False,'numeric_Rp_terminal_join':False,
                'global_five_moment_common_N_closure':False,'exact_heat_energy_stress_flat':False,
                'true_temporal_n_recursion':False,'mean_oscillatory_correction':False,
                'independent_corrected_Cartesian_NS':False,'full_goal_complete':False}}
