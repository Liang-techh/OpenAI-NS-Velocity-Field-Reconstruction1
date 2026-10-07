"""Focused checks of centered phase and correlated inverse product rules."""
import json
from pathlib import Path
import time
import types
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_native_centered_phase_conditioning as current
import lei_ren_part1_paper_compliant_current_native_active_kappa_mixed_conditioning_check as preceding_checks
import lei_ren_part1_paper_compliant_current_native_periodic_mixed_conditioning_check as original_fixtures
from lei_ren_part1_paper_compliant_current_generic_shear_loop import flat_step

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
packets,ep=current.packets,current.ep;require,same=preceding_checks.require,preceding_checks.same


def centered_identity_checks():
    y,Z=sy.symbols('y Z');aa=sy.Function('a')(y,Z);al=sy.Function('alpha')(y,Z);q=al/sy.sqrt(aa)
    first=sy.diff(al,y)-al*sy.diff(aa,y)/(2*aa)
    mixed=sy.diff(al,y,Z)-(sy.diff(al,y)*sy.diff(aa,Z)+sy.diff(al,Z)*sy.diff(aa,y)+al*sy.diff(aa,y,Z))/(2*aa)+3*al*sy.diff(aa,y)*sy.diff(aa,Z)/(4*aa**2)
    require(sy.simplify(sy.sqrt(aa)*sy.diff(q,y)-first)==0,'Exact sqrt(a)*q_y identity failed')
    require(sy.simplify(sy.sqrt(aa)*sy.diff(q,y,Z)-mixed)==0,'Exact sqrt(a)*q_yZ identity failed')
    a,ay,aZ,ayZ,nu,nuy,nuZ,t,ty,tZ,tp,Cy,CZ,CyZ,ch=sy.symbols('a ay aZ ayZ nu nuy nuZ t ty tZ tp Cy CZ CyZ chi')
    F=1+t*t;v=a*nu;vy=ay*nu+a*nuy;vZ=aZ*nu+a*nuZ
    py=-4*nu*Cy/F;pZ=-4*nu*CZ/F
    Lyz=-4*(nuy*CZ+nuZ*Cy+nu*CyZ)
    pyZ=(Lyz-2*t*ty*pZ-2*t*tZ*py-2*t*tp*py*pZ)/F
    Aoriginal=ayZ*ch/2-(ay*pZ+aZ*py+a*pyZ)/(4*sy.pi)
    Acentered=ayZ*ch/2+(v*CyZ+vy*CZ+vZ*Cy)/(sy.pi*F)-2*v*t*(ty*CZ+tZ*Cy)/(sy.pi*F**2)+8*v*nu*t*tp*Cy*CZ/(sy.pi*F**3)
    require(sy.factor(Aoriginal-Acentered)==0,'Complete A_yZ inverse product cancellation failed')
    Minverse=-(ay*t*pZ+aZ*t*py+a*(ty*pZ+tZ*py+tp*py*pZ+t*pyZ))/(2*sy.pi)
    Mcentered=2*v*(1-t*t)*(ty*CZ+tZ*Cy)/(sy.pi*F**2)-8*v*nu*(1-t*t)*tp*Cy*CZ/(sy.pi*F**3)+2*t*(v*CyZ+vy*CZ+vZ*Cy)/(sy.pi*F)
    require(sy.factor(Minverse-Mcentered)==0,'Complete original M_yZ inverse product cancellation failed')
    a=sy.Function('a')(y,Z);b=sy.Function('b')(y,Z);U=sy.Function('qP')(y,Z);psi,phi=sy.symbols('psi phi')
    T1=-b*psi/a+2*U
    for variables in ((y,),(y,Z)):
        require(sy.simplify(-sy.diff(a*T1,*variables)/(2*sy.pi)-sy.diff(b,*variables)*phi
            -(-sy.diff(b,*variables)*(phi-psi/(2*sy.pi))-sy.diff(a*U,*variables)/sy.pi))==0,
            'Centered free-angle M product cancellation failed')
    return dict(passed=True,independent_exact_centered_product_identities=6,
        both_nu_first_cross_terms_in_mixed_L_retained=True,
        complete_original_M_inverse_correction_includes_b_sensitivity=True,
        fixed_free_angle_and_total_inverse_derivatives_distinguished=True)


def transition_v_references():
    p=mp.mp.clone();p.dps=80;count=0
    for eta in map(p.mpf,('.001','.1','.5')):
        def v(k):
            Delta=k-2
            if Delta<=0:return 2+2*eta
            if Delta>=eta:return k
            sigma=flat_step(p,1-Delta/eta)
            return k+sigma*sigma*(2*eta-Delta)
        for fraction in map(p.mpf,('.01','.05','.1','.25','.5','.75','.9','.95','.99')):
            k=2+eta*fraction
            require(2<=v(k)<=3,'Original active v range failed')
            require(abs(p.diff(v,k))<=129,'Original transition v first cap failed')
            require(abs(p.diff(v,k,2))<=11392/eta,'Original transition v second cap failed');count+=3
    return dict(passed=True,independent_active_v_range_and_derivative_comparisons=count,
        original_smooth_cutoff_used=True,body_constant_and_transition_branches_retained=True)


def fixture_references():
    namespace=dict(original_fixtures.primitive_references.__globals__);namespace['current']=current
    fixture=types.FunctionType(original_fixtures.primitive_references.__code__,namespace,'centered_original_scalar_primitive_references')
    return fixture()


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_Z_query_count']==2,'Completed original active-kappa stage required')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed active-kappa original source: '+name)
    require(not saved['actual_controls_or_terminal_closure_installed'] and not any(saved.get(k) for k in packets.OPEN),
        'Mixed bound improvement cannot admit controls/global field/recursion')
    identities=centered_identity_checks();vrefs=transition_v_references();fixtures=fixture_references()
    components=preceding_checks.leading_component_identity_check()
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        targets=current.current.NativeRcParameterTargets(current.current.preceding.NativeRcC1Histories(current.current.preceding.preceding.NativeO2C1Histories(current.current.preceding.preceding.make_middle_owner(bridge))))
        owner=current.NativeCenteredPhaseConditioning(current.previous.NativeSignedAveraging(targets))
        require(saved['exact_correlated_periodic_kernel_theorem']==packets.encode(owner.theorem),'New exact correlation theorem changed')
        direct=json.loads((HERE/current.current.NAME).read_bytes())['actual_original_Rc_parameter_target_records']
        rows=inherited=quiet=pressure=terms=weighted=targets_count=active=0;regions={};direct_comparisons={}
        for region,record in saved['actual_native_centered_phase_conditioning_records'].items():
            got=owner.route(packets.interval(owner.ctx,record['Z_box']))
            require(packets.encode(got['record'])==record,'Same original conditioned route changed: '+region)
            require(len(got['cells'])==24 and record['original_charts']==17,'Whole24-cell/17-chart route required')
            incoming={key:owner.coordinates.scalar(0) for key in current.RATES};incomingZ=dict(incoming)
            for cell in got['cells']:
                label=cell['record']['label'];chart=cell['record']['chart']
                if label!='initial_flat_collar':
                    proof=cell['record']['original_source']['native_periodic_mixed_proof']
                    if not proof['original_whole_support_flat']:
                        require(proof['native_active_kappa_correlation']['active_support_only']
                            and proof['joint_weighted_kernel_proof']['no_q_or_chi_division_used']
                            and proof['centered_original_phase_proof']['a_nu_equals_v_used_before_bounding']
                            and proof['centered_original_phase_proof']['complete_mixed_phase_and_inverse_cross_terms_retained'],
                            'Native correlation scope or flat regularity lost');active+=1
                    count=sum(len(group) for groups in cell['leading_terms'].values() for group in groups.values())
                    require(count==56,'Missing original slow product-rule terms');terms+=count
                for key in current.RATES:
                    decay=cell['factors'][key]['decay']
                    same(cell['incoming'][key],incoming[key],'Incoming C0 history reset')
                    same(cell['incoming_Z'][key],incomingZ[key],'Incoming Z history reset')
                    same(cell['cumulative'][key],decay*incoming[key]+cell['values'][key],'C0 route transport mismatch')
                    same(cell['cumulative_Z'][key],decay*incomingZ[key]+cell['Z_derivatives'][key],'Z route transport mismatch')
                    rows+=4;inherited+=2
                    if chart=='O3_power':
                        require(cell['values'][key].zero and cell['Z_derivatives'][key].zero,'Quiet modulation must be exact zero');quiet+=2
                        if key=='p':
                            same(cell['cumulative'][key],incoming[key],'Quiet C0 pressure memory lost')
                            same(cell['cumulative_Z'][key],incomingZ[key],'Quiet Z pressure memory lost');pressure+=2
                incoming,incomingZ=cell['cumulative'],cell['cumulative_Z']
            require(all(v['new_certificate_strictly_sharper'] for v in record['same_source_averaging_bound_comparison'].values()),
                'All five averaged target certificates should improve')
            require(len(record['first_bridge_weighted_slow_term_budgets'])==56
                and record['slow_term_budgets_are_algebraic_absolute_budget_contributions_not_separate_control_functions'],
                'Complete first-bridge weighted slow attribution required');weighted+=56
            comparison={}
            for key in current.repair.ROWS:
                require(got['targets']['values'][key][-1].zero and got['targets']['Z_derivatives'][key][-1].zero,'N^-2 target order lost');targets_count+=4
                new=current.previous.read_cap(owner.ctx,record['actual_uniform_N_scaled_repair_C1_caps']['transformed_N_scaled_target_C1_caps'][key])
                old=current.previous.read_cap(owner.ctx,direct[region]['actual_uniform_N_scaled_repair_C1_caps']['transformed_N_scaled_target_C1_caps'][key])
                comparison[key]=dict(new_averaging_strictly_sharper_than_direct_at_floor=(new.log is None or (old.log is not None and ep(new.log)[1]<ep(old.log)[0])),
                    direct_N_scaled_C1_cap=old.record(),new_averaged_N_scaled_C1_cap_at_floor=new.record(),bound_comparison_not_actual_error=True)
            direct_comparisons[region]=comparison
            regions[region]=dict(all_five_averaged_targets_improved=True,dominant_cells=record['dominant_target_bound_cells'],
                dominant_first_bridge_components=record['first_bridge_dominant_budget_components'],dominant_first_bridge_slow_terms=record['first_bridge_dominant_slow_terms'],
                whole_N_scaled_target_C1_cap=record['actual_uniform_N_scaled_repair_C1_caps']['whole_target_C1_cap'])
            print('Native centered-phase correlation checked:',region,flush=True)
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,native_Z_queries_checked=2,
        original_true_cells_checked=24,original_charts_checked=17,actual_C0_Z_transport_rows_checked=rows,
        actual_inherited_rows_checked=inherited,exact_quiet_rows_checked=quiet,preserved_pressure_memory_rows_checked=pressure,
        active_correlated_primitive_cells_checked=active,original_leading_slow_source_terms_checked=terms,
        first_bridge_weighted_slow_terms_checked=weighted,normalized_target_rows_checked=targets_count,
        independent_centered_phase_identity_checks=identities,independent_transition_v_derivative_references=vrefs,
        accepted_unchanged_weighted_kernel_receipt=current.preceding.RECEIPT,independent_exact_leading_components=components,
        independent_original_scalar_primitive_references=fixtures,exact_correlated_periodic_kernel_theorem=owner.theorem,
        same_source_direct_floor_comparisons=direct_comparisons,regions=regions,
        actual_controls_or_terminal_closure_installed=False,point_inverse_or_higher_jets_installed=False,
        one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes,scope=saved['scope'])
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Native centered-phase correlation PASS',rows,inherited,terms,weighted,targets_count,flush=True)
    return result


if __name__=='__main__':run()
