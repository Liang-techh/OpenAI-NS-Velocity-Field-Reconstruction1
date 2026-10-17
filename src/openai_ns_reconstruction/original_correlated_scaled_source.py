"""Same-history normalized C5 pressure differences with physical scales split.

The correction is recovered directly from signed V and mean-history ports.
Independent full/background boxes are never subtracted to manufacture it.
"""
from __future__ import annotations

from pathlib import Path
import hashlib
import json

from .original_finite_N_signed_source import (
    COMPACT,RM,RM_BASE,load,sha,digest,require,ast_guard)
from .original_signed_density_expression import SignedDensityExpressionGraph

PARENT='lei_ren_part1_paper_compliant_current_original_finite_N_signed_source.json'
PARENT_SHA256='3c309218eb3798ba077f2ae0c1c74be9e5f0a0dce023e38c3e1dee3283d309c7'


def build_core_pressure_identity(root,parent):
    """Bind the two *true* pressure functions; retain numerical enclosures.

    The q-coefficient ODE proof is for axial Taylor coefficients, not a
    temporal Navier--Stokes coefficient recursion.
    """
    from .original_inner_collar_source import _ast_guard
    import sympy as s
    root=Path(root);work=root/'experiments/root_st073'
    names={
        'join':'lei_ren_part1_paper_compliant_actual_Rp_source_join.json',
        'native':'lei_ren_part1_paper_compliant_current_limit_Rp_native_identity.py',
        'scaling':'lei_ren_part1_paper_compliant_current_core_recurrence_source_check.json',
        'pressure':'lei_ren_part1_paper_compliant_pressure_source.json',
        'q':'lei_ren_part1_paper_candidate_pressure_function.py',
        'canonical':'lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions.py',
    }
    pins={}
    for role,name in names.items():
        actual=sha(work/name);require(actual==parent['input_hashes'][name],'Changed core pressure source: '+name)
        pins[role]=dict(path='experiments/root_st073/'+name,sha256=actual)
    forward='src/openai_ns_reconstruction/original_core_numeric.py'
    require(sha(root/forward)==parent['implementation_hashes'][forward],'Changed forward core pressure seed')
    pins['forward']=dict(path=forward,sha256=sha(root/forward))
    proof=load(work/names['join'])['current_Rh_pressure_defining_function_proof']
    require(all(proof[k] is True for k in ('passed','same_fourteen_implicit_stage_sources',
        'common_flatten_function_retained_not_replaced_by_zero','exact_same_inherited_normalized_jets_callable',
        'output_enclosures_not_selected_pressure_values','callable_constructor_and_equations_not_hash_only_proof')),
        'Accepted true analytic native pressure identity missing')
    require(proof['exact_analytic_defining_function']=='P0/Pstar^2=-(M2*(1+Z^2)^-2+M0+F_flat(Z))',
        'Different analytic native pressure function')
    pressure=load(work/names['pressure'])['compliant_source']
    definition_hash=hashlib.sha256(json.dumps(pressure['implicit_source_definition'],sort_keys=True).encode()).hexdigest()
    require(definition_hash==pressure['implicit_source_sha256']==parent['source_family']['implicit_source_sha256'] and
        pressure['datum_enclosure_sha256']==parent['source_family']['datum_enclosure_sha256'] and
        pressure['stage_count_total']==14,'Pressure definition/family mismatch')
    guards={
        'q_initial_and_recurrence':ast_guard(root,'experiments/root_st073/'+names['q'],'q_jets',(
            'a = 1 + center**2','b = 2*center','coefficients = [1/a**2]','previous = ctx.mpf(0)',
            'coefficients.append(-(b*(n+2)*current+(n+3)*previous)/(a*(n+1)))','previous = current')),
        'forward_seed':_ast_guard(root/forward,root,'_seed',cls='OriginalCoreNumeric',text=(
            'qs = self.q_jets(c,z0,n)', 'M2 = _box(c,self.pressure["fixed_beta2_mass_normalized"])',
            'M0 = _box(c,self.pressure["fixed_beta0_mass_normalized"])',
            'flat = _box(c,self.pressure["stages"]["z_flatten"]["mass"])',
            'pressure.append(-(M2*q+(M0 if k == 0 else 0)+remainder))',
            'P0_Z_taylor=[c.exp(2*lp-ll)*v for v in pressure]')),
        'canonical_physical_pressure':_ast_guard(work/names['canonical'],root,'__init__',
            cls='CurrentOriginalN1CoreFixedPointFunctions',text=(
                "P0=a.scale(self.pulse.own.jets.native['P0'],g.unary('exp',g.mul(c(2),logP)))",)),
        'native_true_pressure':ast_guard(root,'experiments/root_st073/'+names['native'],'pressure_source_proof',(
            "P0=-(M2/(1+Z**2)**2+M0+flat)",
            'common_P0_kept_independent_of_cumulative_pressure=True')),
    }
    require(all(v['passed'] is True for v in guards.values()),'Cross-producer pressure AST binding failed')
    a,b,t,n,Cprev,C,Cnext=s.symbols('a b t n Cprev C Cnext')
    q=(a+b*t+t*t)**(-2)
    ode=s.simplify((a+b*t+t*t)*s.diff(q,t)+2*(b+2*t)*q)==0
    coefficient=a*(n+1)*Cnext+b*n*C+(n-1)*Cprev+2*b*C+4*Cprev
    recurrence=s.simplify(coefficient-(a*(n+1)*Cnext+b*(n+2)*C+(n+3)*Cprev))==0
    initial=s.simplify(q.subs(t,0)-a**(-2))==0
    require(ode and recurrence and initial,'All-order q source ODE identity failed')
    scaling=load(work/names['scaling'])['pressure_and_seed_bindings']
    required=('Lambda_times_core_epsilon_exactly_one',
        'P0_step_equals_core_epsilon_times_Pstar_squared_times_normalized_pressure',
        'no_second_epsilon_applied_to_ell','core_radial_epsilon_distinct_from_pressure_datum_parameter_epsilon')
    require(all(scaling[k] is True for k in required),'Core pressure scaling identity unavailable')
    return dict(passed=True,source_executed=False,source_pins=pins,AST=guards,
        true_normalized_pressure_function='p(Z)=-(M2/(1+Z^2)^2+M0+F_flat(Z))',
        true_forward_and_canonical_pressure_same_function=True,
        analytic_definition_sha256=definition_hash,all_fourteen_pressure_sources_retained=True,
        all_order_q_Taylor_identity=dict(passed=True,ODE='(a+b*t+t^2)q_t+2(b+2t)q=0',
            a='1+Z0^2 > 0 for real Z0',b='2Z0',initial='c0=a^-2; c[-1]=0',
            recurrence='a(n+1)c[n+1]+b(n+2)c[n]+(n+3)c[n-1]=0',
            local_analytic_uniqueness_implies_all_orders=True,ordinary_derivative_from_Taylor='q^(k)=k!*c[k]',
            original_q_jets_or_core_recurrence_executed=False),
        scaling=dict(canonical='physical_P0=Pstar^2*p',forward='P0_Z_taylor=epsilon_core*physical_P0_Taylor',
            epsilon_core='exp(-logLambda)',pressure_datum_epsilon_is_distinct=True,
            source_scaling_flags={k:scaling[k] for k in required}),
        forward_directed_rows_are_enclosures_of_same_true_function=True,
        forward_point_approximations_equal_true_function_claim=False,
        flat_remainder_or_mass_enclosures_selected_as_values=False,
        canonical_core_equals_forward_full_solution_proved=False,delta_six_zero_claim=False,
        temporal_coefficient_recursion_claim=False)


def build_correlated_scaled_source(root,parent):
    from .original_finite_N_density_kernel import ordinary_jet
    from .original_scaled_signed_density import correlated_scaled_density
    root=Path(root);work=root/'experiments/root_st073'
    require(sha(work/PARENT)==PARENT_SHA256,'Accepted finite-N source parent required')
    require(digest(load(work/PARENT))==digest(parent),'Supplied parent is not the accepted source')
    reports={}
    for name in (COMPACT,RM,RM_BASE):
        require(sha(work/name)==parent['input_hashes'][name],'Changed actual function source')
        reports[name]=load(work/name)
    compact=reports[COMPACT]['true_compact_C5_bindings'];rm=reports[RM]['Rm_patch_finite_N_C5_bindings']
    rm_base=reports[RM_BASE]['true_Rm_patch_C6_bindings']['chart']
    csrc='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions.py'
    rsrc='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_true_Rm_patch_finite_N_C5_functions.py'
    radial='experiments/root_st073/lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions.py'
    for rel in (csrc,rsrc,radial):require(sha(root/rel)==parent['input_hashes'][Path(rel).name],'Changed signed recovery source')
    guards={
        'compact_profiles':ast_guard(root,csrc,'profiles',('dV=a.mul(a.trim(A),G)',"radial.append(dict(F=F,G=G,original_E=E,delta_E=dE,delta_V=dV,E=a.add(E,dE),V=dV))")),
        'total_radial':ast_guard(root,radial,'total_radial_functions',("Vy=a.add(a.trim(bg['radial']['V_y']),a.scale(slow['B_over_Pstar'],invN),chart['loop']['B_over_Pstar_phi'])",)),
        'physical_Q':ast_guard(root,radial,'physical_pressure_fields',('my=a.sub(U,m)',"Q=div(numerator(U,m));Qy=div(numerator(Uy,my))")),
    }
    # The full V and m sums are checked as actual graph nodes, so their common
    # background cannot be independently re-boxed by this interface.
    charts={}
    for label,name,b,bg in (('compact',COMPACT,compact,compact['background']),('Rm_patch',RM,rm,rm_base)):
        raw=reports[name];nodes=raw['exact_graph_nodes'];g=SignedDensityExpressionGraph()
        ref=lambda n:g.reference(name,parent['input_hashes'][name],n,source_nodes=nodes)
        rows=lambda q:[ref(n) for n in q]
        pair=lambda q:[q.value.node,q.dz.node]
        if label=='compact':
            dV=b['profiles']['radial_y_rows'][0]['delta_V'];dVy=b['profiles']['radial_y_rows'][1]['delta_V']
            dm=b['history']['m'];full=b['modified_fields'];correction_V_y=rows(dVy[:2])
            require(full['V']==dV,'Compact V is not the signed source correction')
            radial_rule=dict(kind='exact_compact_delta_V_y_rows',rows=dVy[:2],additional_N_division=False)
        else:
            dV=b['chart']['finite_N']['delta_V'];dm=b['chart']['history']['m'];full=b['chart']['modified_fields']
            slow=b['radial']['slow_B_y'];fast=b['chart']['loop']['B_over_Pstar_phi'];N=ref(b['N'])
            correction_V_y=[ref(slow[k])/N+ref(fast[k]) for k in (0,1)]
            radial_rule=dict(kind='slow_B_y_over_N_plus_B_phi',slow=slow[:2],fast=fast[:2],N=b['N'],
                original_radial_N_division_applied_once=True,delta_V_not_divided_again=True)
        sum_proofs=[]
        for key,increment,order in (('V',dV,3),('m',dm,4)):
            for k in range(order):
                expected=[bg['fields'][key][k],increment[k]]
                node=nodes[full[key][k]]
                # The graph may suppress an exact zero background summand.
                if bg['fields'][key][k]==0:ok=full[key][k]==increment[k]
                else:ok=node.get('operation')=='sum' and sorted(node.get('arguments',()))==sorted(expected)
                require(ok,'Unproved same-history full='+key+'0+delta'+key+' at '+str(k))
                sum_proofs.append(dict(field=key,ordinary_Z_order=k,full=full[key][k],summands=expected,passed=True))
        z=ordinary_jet(ref(20),g.one,g.zero);de=ref(25220)
        result=correlated_scaled_density(rows(bg['fields']['V'][:3]),rows(bg['generic']['Q'][:3]),
            rows(dV[:3]),rows(dm[:4]),correction_V_y,z,de)
        # Pstar and R stay named exact source atoms. Numerical work can act on
        # bounded normalized coefficients without exponentiating these atoms.
        R=parent['finite_N_signed_source']['charts'][label]['physical_radius']
        charts[label]=dict(report=name,report_sha256=parent['input_hashes'][name],
            physical_scale_atoms=dict(Pstar=33572,radius=R,N=18),
            original_signed_ports=dict(delta_V=dV[:3],delta_m=dm[:4],baseline_V=bg['fields']['V'][:3],
                baseline_q=bg['generic']['Q'][:3],delta_V_y=radial_rule),
            same_history_sum_proofs=sum_proofs,
            normalized_local_roots={key:pair(result[key]) for key in ('delta_q','delta_q_y','v1','v2','F1','F2')},**g.packet(),
            physical_recombination=dict(delta_v='Pstar*v1+Pstar^2*v2',delta_F='Pstar*F1+Pstar^2*F2',
                delta_Q='Pstar*delta_q',delta_Q_y='Pstar*delta_q_y',
                delta_boundary='R^2*(Pstar*F1+Pstar^2*F2)-R*Pstar*delta_q_y'),
            ordinary_derivatives=True,full_minus_background_boxes_used=False,
            numerical_original_source_values_installed=False,caps_removed_or_narrowed=False,
            nonzero_incoming_pressure_retained=True)
    return dict(passed=True,source_executed=False,charts=charts,AST=guards,
        direct_signed_history_recovery=True,physical_scale_materialized=False,
        numerical_original_source_values_installed=False,finite_N_caps_active=True)
