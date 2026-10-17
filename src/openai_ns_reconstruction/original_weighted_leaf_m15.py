"""Signed raw46/raw47 source leaves and a reference-preserving m15 view.

This adapter hydrates the accepted weighted queries and original ledger.
The replaced positive remainder is retained as provenance, never selected
as a value.  All other ledger rows and the independent target remain exact
references to their pinned parent bytes.
"""
from __future__ import annotations

import copy
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from mpmath.libmp import mpf_cmp

from .original_hydro_source_point import ROOT,DIRECTORY,sha,digest,load,exact_iv,upper
from .original_leading_viscous_kernel import decode_pair,encode_pair
from .original_modal_source_point import _function_ast_guard
from .original_weighted_source_primitive import serial

PARENT='lei_ren_part1_paper_compliant_current_original_weighted_source_primitive_numeric.json'
PARENT_SHA256='76e686626ca0310e56751f28a6e80194d77e5852ce2adf3bbf679ee1cf687f5b'
SOURCE='lei_ren_part1_paper_compliant_current_original_swirl_source_C1_numeric.json'
FTC='lei_ren_part1_paper_compliant_current_original_leading_viscous_source_numeric.json'
FIRST='lei_ren_part1_paper_compliant_current_original_first_collar_jets.json'
MASS='lei_ren_part1_paper_compliant_current_original_source_jet_mass_numeric.json'
NAMES=('first_switch_viscous_source_integral','first_switch_viscous_source_integral_error')


def materialize_view(parent_terms,entries):
    """Expand references without recalculation; changed rows stay explicit."""
    if len(parent_terms)!=183 or len(entries)!=183:raise ValueError('Original 183-slot view required')
    output=[]
    for index,entry in enumerate(entries):
        if entry['parent_reduced_index']!=index:raise ValueError('Ledger order changed')
        old=parent_terms[index]
        if digest(old)!=entry['parent_term_digest']:raise ValueError('Changed original ledger leaf')
        if entry['operation']=='unchanged_parent_reference':output.append(copy.deepcopy(old))
        elif entry['operation']=='same_source_signed_leaf_replacement':output.append(copy.deepcopy(entry['term']))
        else:raise ValueError('Unknown ledger operation')
    return output


class OriginalWeightedLeafM15:
    def __init__(self,*,directory=DIRECTORY):
        self.directory=Path(directory)
        if sha(self.directory/PARENT)!=PARENT_SHA256:raise ValueError('Pinned accepted weighted source report required')
        self.parent=load(self.directory/PARENT)
        if self.parent['checks']['passed'] is not True:raise ValueError('Passed weighted source report required')
        self.input_hashes=dict(self.parent['input_hashes']);self.input_hashes[PARENT]=PARENT_SHA256
        self.implementation_hashes=dict(self.parent['implementation_hashes'])
        for base,group in ((self.directory,self.input_hashes),(ROOT,self.implementation_hashes)):
            for name,expected in group.items():
                if sha(base/name)!=expected:raise ValueError('Changed accepted pin: '+name)
        self.ledger_name=self.parent['inherited_m15_packet']['report']
        for name in (SOURCE,FTC,FIRST,MASS,self.ledger_name):
            if name not in self.input_hashes:raise ValueError('Missing inherited source pin: '+name)
        self.source=load(self.directory/SOURCE);self.ftc=load(self.directory/FTC);self.ledger=load(self.directory/self.ledger_name)
        self.first=load(self.directory/FIRST);self.mass=load(self.directory/MASS)
        if self.parent['inherited_m15_packet']['report_sha256']!=self.input_hashes[self.ledger_name]:raise ValueError('Ledger pin changed')
        self.ctx=mp.mp.clone();self.ctx.dps=180
        self.iv=MPIntervalContext();self.iv.dps=220

    def source_guards(self):
        src=ROOT/'src/openai_ns_reconstruction';base=ROOT/'experiments/root_st073'
        specs=(
            ('actual_source_baseline',src/'original_first_switch_m15.py','leading','OriginalFirstSwitchM15',
             ["v=self.rows.inputs(Z,s)","macro=v['incoming']","low=macro['source_Uy']", "C(iv,low['truth'],o.core.delta_iv)*norms['mass']",
              "((B.truncate(2)*(zz*(1-delta))+o._derivative(B.truncate(2))*(1-zz*zz))/(1-zz*zz*delta)).truncate(1)",
              "coefficient(iv,low['truth'],o.core.delta_iv)"]),
            ('original_leaf_append',src/'original_first_switch_m15.py','evaluate','OriginalFirstSwitchM15',
             ["append('first_switch_viscous_source_integral',zero if empty else lead['vis_point'],zeroi if empty else lead['vis_truth'],[physical(1,1),hbS,hbS],actual_positive_complement_mass=True,incoming_mean_cancelled_by_exact_ODE_only=True)",
              "append('first_switch_viscous_source_integral_error',zero,zeroi if empty else lead['vis_remainder'],[physical(1,1),hbS,hbS,hbS],error_enclosure_only=True)",
              "Rpower*iv.log(100)-degree*o.core.ll_iv-mode*self.rows.micro.Y"]),
            ('actual_source_modal_weight',src/'original_first_switch_m15.py','source_norms','OriginalFirstSwitchM15',
             ["rate=power-j+radius_weight", "p['R0']**j*b.ctx.mpf(100)**(power-j)",
              "mass,_=b.pulse_masses(phase)","m=o._transfer(mass)",
              "variation=s*m*weighted_derivative+o.hcap*angular*iv.exp(ell)*m*weighted"]),
            ('shared_macro_source_port',src/'original_first_switch_pressure_rows.py','inputs','OriginalFirstSwitchPressureRows',
             ["incoming=self.macro_rows.inputs(Z,iv.mpf(1))"]),
            ('actual_B_modes',src/'original_macro_pressure_rows.py','source_Uy','OriginalMacroPressureRows',
             ["modes=[prepared['R0']**j*R**(power-j) for j in range(3)]",
              "row=-prepared['quotient']*module.weighted_modes(prepared['modes']['drive_'+name],modes)",
              "term=b.source_term(row,1,scale)",
              "o._scaled(v,term,divide_eps=True)","truth=total+o._norm_error_jet(nonlinear,order)"]),
            ('physical_density',base/'lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions.py',
             'pressure_density',None,["a.scale(Py,g.neg(R))"]),
            ('physical_FTC',base/'lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions.py',
             'pressure_FTC',None,["a.scale(a.add(Qy,Qyy),R)"]),
        )
        guards={name:_function_ast_guard(path,fn,fragments,class_name=cls) for name,path,fn,cls,fragments in specs}
        for name,path,_,_,_ in specs:
            expected=self.input_hashes.get(path.name) if path.parent==base else self.implementation_hashes.get(path.relative_to(ROOT).as_posix())
            if expected!=sha(path):raise ValueError('Unbound leaf source AST: '+name)
        inherited=self.ftc['independent_leading_FTC']
        checks={'weighted_source_bound':self.parent['source_guards']['passed'],
                'same_family':self.parent['source_family']==self.ledger['source_family']==self.ftc['source_family']==self.source['source_family'],
                'independent_leading_FTC':inherited['passed'] and inherited['ordinary_Z1'] and inherited['physical_R_once'] and inherited['physical_Pstar_once'],
                'mean_ODE_partition':inherited['mean_ODE_partition']=='Hlower-Hupper=Elower-Eupper+integral R*C(Uy) dy',
                'no_modified_Qy_alias':all(not v['modified_Qy_alias'] for v in inherited['endpoint_recipes'].values()),
                'same_actual_M0_function':self.parent['source_guards']['checks']['same_actual_flat_cutoff'],
                'same_raw213_reduced183_witness':self.parent['inherited_m15_packet']==self.source['inherited_m15_packet']}
        return {'passed':all(v['passed'] for v in guards.values()) and all(checks.values()),'AST':guards,'checks':checks,
            'independent_leading_FTC_digest':digest(inherited),
            'macro1_first0_leading_seam_digest':inherited['macro1_first0_leading_input_identity_sha256'],
            'identity':{'actual_viscous_FTC':'Hlower-Hupper=Elower-Eupper+integral R*C(Uy) dy',
                'source_primitive':'integral_0^s R*C(Uy)dy=R100*epsilon*hbS^2*C(J)',
                'signed_split':'C(J)=C(sum_j B_j)*M0+hbS*C(R)',
                'R_definition':'R=(sum_j B_j*K_(q-j+1)-sum_j B_j*M0)/hbS',
                'physical_main':'R100*epsilon*hbS^2','physical_remainder':'R100*epsilon*hbS^3',
                'M0_and_width_and_R_epsilon_axially_constant':True,
                'same_shared_lower_source':'macro.source_Uy is the first-switch lower source port',
                'shared_B_identity':'B_low=sum_j B_j as signed source functions; numerical truth boxes are enclosures and need not be equal',
                'B_j_source_recipe':'-(Phi0/barPhi2)*exp(ell_macro)*(sector_scale/epsilon)*R0^j*100^(q-j)*D_sector,j',
                'sector_scales':{'hydro':'1','pressure':'Pstar^2','swirl':'F0(Z)^2'},
                'actual_Uy_recipe':'epsilon*hbS*kappa(s)*exp(ell_local(s))*sum_j B_j*exp((q-j)*hbS*s)',
                'M0_source_recipe':'integral_0^s (1-sigma(t))dt from the original gp ratio',
                'source_recipes_not_registered_graph_nodes':True,
                'old_positive_raw47_bound_is_not_a_signed_function':True,
                'raw43_raw44_raw45_and_raw23_raw24_retained':True}}

    def log_cap(self,term,target,order,*,error_only=False):
        c,iv=self.ctx,self.iv
        unit=sum((exact_iv(c,iv,row) for row in target['log_scale_terms']),iv.mpf(0))
        scale=sum((exact_iv(c,iv,row) for row in term['log_scale_terms']),iv.mpf(0))-unit
        row=term['coefficient_C0_C1'][order]
        point,truth=decode_pair(c,iv,row)
        magnitude=upper(c,truth-iv.mpf(row['approximation']) if error_only else truth)
        if magnitude==0:return None
        coefficient_log=iv.log(iv.mpf(magnitude));result=scale+coefficient_log
        return {'upper_approximation':c.nstr(c.make_mpf(result._mpi_[1]),45),
                'upper_exact_mpf_tuple':list(result._mpi_[1]),'bound_only_not_a_function_value':True,
                'normalization':'independent complete pulse m15 target unit reference*A*Rv^2',
                'coefficient_log_upper_exact_mpf_tuple':list(coefficient_log._mpi_[1]),
                'coefficient_log_upper_approximation':c.nstr(c.make_mpf(coefficient_log._mpi_[1]),45),
                'coefficient_error_only':error_only}

    def build(self):
        from .original_weighted_leaf_kernel import build_factored_leaf_records
        guards=self.source_guards()
        if not guards['passed']:raise ValueError('Same-source weighted leaf binding failed: '+str(guards))
        main={(r['Z'],r['phase']):r for r in self.source['signed_source_integral_rows']}
        sources={r['Z']:r for r in self.source['all_sector_source_rows']}
        masses={self.ctx.mpf(r['phase']):r for r in self.mass['interior_mass']['rows']}
        queries=[]
        for row in self.parent['weighted_remainder_queries']:
            key=row['Z'],row['phase'];baseline=main[key]
            sample=next(r for r in self.first['samples'] if r['Z']==row['Z'] and r['phase']=='0')
            lp,lt=decode_pair(self.ctx,self.iv,sample['source_scales']['log_epsilon'])
            bs=[decode_pair(self.ctx,self.iv,r) for r in sources[row['Z']]['source_rows']]
            rs=[decode_pair(self.ctx,self.iv,r) for r in row['signed_R_ordinary_rows']]
            mass=decode_pair(self.ctx,self.iv,masses[self.ctx.mpf(row['phase'])])
            evaluated=[]
            for index,ctx in enumerate((self.ctx,self.iv)):
                de=ctx.exp((lp if index==0 else lt)+970)
                evaluated.append(build_factored_leaf_records([v[index] for v in bs],[v[index] for v in rs],
                    ctx.mpf(row['Z']),de,mass[index],R100=ctx.mpf(100),epsilon=ctx.exp(lp if index==0 else lt),
                    width_symbol='hbS',source_proved=guards['passed'],shared_source_proved=True,mass_Z_independent=True))
            leaves=[];regressions=[]
            for rawindex,label,reference in ((46,'raw46',baseline['all_sector_main_component_C0_C1']),(47,'raw47',row['raw47_signed_C0_C1'])):
                pp=evaluated[0][label]['unscaled_coefficient_rows'];tt=evaluated[1][label]['unscaled_coefficient_rows']
                encoded=[encode_pair(self.ctx,self.iv,p,t) for p,t in zip(pp,tt)]
                leaves.append({'raw_index':rawindex,'coefficient_C0_C1':encoded,'width_power':rawindex-44,
                    'source_ordinary_rows':sources[row['Z']]['source_rows'] if rawindex==46 else row['signed_R_ordinary_rows'],
                    'physical_atom':'R100*epsilon*hbS^'+str(rawindex-44)})
                for k,(p,ref) in enumerate(zip(pp,reference)):
                    rp,_=decode_pair(self.ctx,self.iv,ref)
                    err=abs(p-rp)/max(self.ctx.mpf(1),abs(rp))
                    if err>self.ctx.mpf('1e-150'):raise ValueError('Inherited signed operator regression failed')
                    regressions.append({'raw_index':rawindex,'ordinary_Z_order':k,'passed':True,
                        'scaled_approximation_difference':self.ctx.nstr(err,25),'threshold':'1e-150',
                        'numerical_regression_not_function_identity':True})
            pair={'passed':True,'leaves':leaves,'width_symbol':'hbS','width_materialized':False,
                  'source_rows_preserved_through_Z2':True,'same_source_identity_proved':True,
                  'operator_regressions':regressions,'combined_factored_identity':'R100*epsilon*hbS^2*(main+hbS*remainder)',
                  'independent_kernel_proof':serial(evaluated[1]['proof'])}
            queries.append({'Z':row['Z'],'phase':row['phase'],'signed_pair':pair,
                'weighted_query_digest':digest(row),'baseline_main_digest':digest(baseline['all_sector_main_component_C0_C1']),
                'same_source_FTC_identity_digest':digest(guards['identity']),
                'old_positive_raw47_selected_as_value':False,'query_is_global_closure':False})
        views=[]
        witnesses={r['Z']:r for r in self.parent['inherited_m15_packet']['samples']}
        for sample in self.ledger['samples']:
            z=sample['Z'];raw=sample['raw_m15_C0_C1_terms'];old=sample['correlated_reduced_m15_C0_C1_terms']
            witness=witnesses[z]
            if len(raw)!=213 or len(old)!=183 or digest(raw)!=witness['raw_sha256'] or digest(old)!=witness['reduced_sha256']:
                raise ValueError('Original cumulative ledger changed')
            if digest(raw[:6])!=witness['signed_inner_lower_six_sha256'] or digest(sample['independent_target'])!=witness['independent_target_sha256']:
                raise ValueError('Inner/lower or independent target changed')
            query=next(r for r in queries if r['Z']==z and self.ctx.mpf(r['phase'])==1)
            replacement_rows=query['signed_pair']['leaves']
            entries=[];provenance=[];local=[]
            for index,term in enumerate(old):
                entry={'parent_reduced_index':index,'parent_term_digest':digest(term),'name':term['name'],
                       'operation':'unchanged_parent_reference'}
                if term['name'] in NAMES:
                    k=NAMES.index(term['name']);rawindex=46+k
                    if term!=raw[rawindex]:raise ValueError('Unexpected raw/reduced leaf mapping')
                    new=copy.deepcopy(term);new['coefficient_C0_C1']=copy.deepcopy(replacement_rows[k]['coefficient_C0_C1'])
                    if k==1:new['name']='first_switch_viscous_source_integral_signed_remainder'
                    new['actual_function_binding']=dict(new['actual_function_binding'],
                        error_enclosure_only=False,actual_signed_weighted_source_leaf=True,
                        original_positive_bound_provenance_only=(k==1),same_source_identity_digest=digest(guards['identity']),
                        ordinary_Z_derivatives=[0,1],positive_amplitude_retained=True)
                    entry.update(operation='same_source_signed_leaf_replacement',term=new,raw_index=rawindex)
                    provenance.append({'raw_index':rawindex,'reduced_index':index,'parent_term_digest':digest(term),
                                       'original_term_retained_in_pinned_parent':True,
                                       'old_error_enclosure_selected_as_value':False})
                    local.append({'raw_index':rawindex,'name':term['name'],
                        'before_caps':[self.log_cap(term,sample['independent_target'],k) for k in (0,1)],
                        'after_caps':[self.log_cap(new,sample['independent_target'],k) for k in (0,1)],
                        'before_error_caps':[self.log_cap(term,sample['independent_target'],k,error_only=True) for k in (0,1)],
                        'after_error_caps':[self.log_cap(new,sample['independent_target'],k,error_only=True) for k in (0,1)]})
                entries.append(entry)
            expanded=materialize_view(old,entries)
            if len(provenance)!=2 or sum(t['name']=='reserved_n1_own_swirl_once' for t in expanded)!=1:raise ValueError('Required leaf/reserved multiplicity changed')
            unchanged=[i for i,t in enumerate(entries) if t['operation']=='unchanged_parent_reference']
            if len(unchanged)!=181 or any(expanded[i]!=old[i] for i in unchanged):raise ValueError('Foreign ledger leaf changed')
            views.append({'Z':z,'parent_report':self.ledger_name,'parent_reduced_terms_digest':digest(old),
                'raw213_digest':digest(raw),'archived_reduced183_digest':digest(old),'view_term_count':183,
                'unchanged_parent_references':181,'entries':entries,'replacement_provenance':provenance,
                'expanded_signed_view_digest':digest(expanded),'six_inner_lower_digest':digest(raw[:6]),
                'independent_target_digest':digest(sample['independent_target']),
                'independent_target':copy.deepcopy(sample['independent_target']),
                'selected_N_evidence_digest':digest(sample['selected_N_evidence']),
                'query_pair_digest':digest(query),'new_local_cap_comparisons':local,
                'signed_discrepancy_recipe':'sum(expanded183 signed physical leaves)-unchanged independent negative complete-pulse m15 target',
                'largest_unchanged_caps':copy.deepcopy(sample['largest_remaining_individual_log_caps']),
                'closure_or_small_residual_claim':False})
        return {'source_family':self.parent['source_family'],'parent_report':PARENT,'parent_report_sha256':PARENT_SHA256,
            'source_guards':guards,'signed_leaf_queries':queries,'signed_m15_query_views':views,
            'inherited_m15_packet':self.parent['inherited_m15_packet'],
            'inherited_m15_packet_digest':self.parent['inherited_m15_packet_digest'],
            'inherited_weighted_queries_digest':digest(self.parent['weighted_remainder_queries']),
            'inherited_weighted_moment_digest':digest(self.parent['weighted_pulse_moments']),
            'contract':{'accepted_ancestor_or_fixture_replayed':False,'M0_W1_W2_or_core_recomputed':False,
                'old_raw213_and_reduced183_archives_unchanged':True,'other181_leaves_unchanged':True,
                'six_inner_lower_reserved_once_target_and_selected_N_unchanged':True,
                'raw23_raw24_raw43_raw44_raw45_and_all_other_errors_retained':True,
                'positive_caps_or_midpoints_selected_as_function_values':False,
                'hb_or_amplitude_materialized':False},
            'scope':{'signed_raw46_raw47_pairs_at19queries':True,'signed_Rp_ledger_views_at3Z':True,
                'same_source_FTC_partition_bound':True,'raw47_signed_replacement_in_query_views':True},
            'limits':{'global_axial_function_numeric_primitive':False,'numeric_Rp_terminal_join':False,
                'global_five_moment_closure':False,'exact_heat_energy_stress_flat':False,
                'true_temporal_n_recursion':False,'mean_oscillatory_correction':False,
                'independent_corrected_Cartesian_NS':False,'full_goal_complete':False}}
