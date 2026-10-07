"""Focused source correlation, weighted kernel and native route checks."""
import json
from pathlib import Path
import time
import types
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_native_active_kappa_mixed_conditioning as current
import lei_ren_part1_paper_compliant_current_native_periodic_mixed_conditioning_check as preceding_checks

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
packets,ep=current.packets,current.ep
require,same=preceding_checks.require,preceding_checks.same


def weighted_kernel_references():
    p=mp.mp.clone();p.dps=70;count=0;angles=0
    def theta(u,psi):
        h=p.sqrt(1+u*u);r=u/h
        return 2*p.atan2((1+r)*p.sin(psi/2),(1-r)*p.cos(psi/2))
    def d(u,psi):
        h=p.sqrt(1+u*u);r=u/h
        return (p.cos(psi)-r)/(h*(1-2*r*p.cos(psi)+r*r))
    def points(end):
        return sorted(set([p.mpf(0),end]+[v for v in (p.mpf('.01'),p.mpf('.1'),p.pi-p.mpf('.1'),p.pi,p.pi+p.mpf('.1')) if 0<v<end]))
    for u in map(p.mpf,('0','.4','-.4','2','-2','20','-20')):
        h=p.sqrt(1+u*u)
        for psi in map(p.mpf,('.73','4.9')):
            du=lambda x:p.diff(lambda z:d(z,x),u)
            duu=lambda x:p.diff(lambda z:d(z,x),u,2)
            values=(p.quad(lambda x:d(u,x),points(psi)),p.quad(du,points(psi)),p.quad(duu,points(psi)),
                p.quad(lambda x:d(u,x)**2,points(psi)),p.quad(lambda x:2*d(u,x)*du(x),points(psi)),
                p.quad(lambda x:2*(du(x)**2+d(u,x)*duu(x)),points(psi)))
            caps=(4*p.pi/h,(80*p.pi/3)/h**2,160*p.pi/h**3,p.pi,60*p.pi/h,800*p.pi/h**2)
            for value,cap in zip(values,caps):
                require(abs(value)<=cap+p.mpf('1e-45'),'Independent weighted partial-angle kernel cap failed');count+=1
            th=theta(u,psi)
            require(abs(p.diff(lambda z:theta(z,psi),u)-2*p.sin(th)/h)<p.mpf('1e-45'),'Exact signed theta_u identity failed')
            require(abs(p.diff(lambda z:theta(z,psi),u,2)-(4*p.sin(th)*p.cos(th)/h**2-2*u*p.sin(th)/h**3))<p.mpf('1e-45'),
                'Exact signed theta_uu identity failed');angles+=2
    return dict(passed=True,independent_weighted_partial_angle_P_H_derivative_comparisons=count,
        independent_signed_theta_derivative_identities=angles,
        signed_u_and_both_r_sectors_and_partial_angles_tested=True,uniformity_comes_from_analytic_sector_proofs=True)


def leading_component_identity_check():
    c=MPIntervalContext();c.dps=100;y,Z=sy.symbols('y Z');f={name:sy.Function(name)(y,Z) for name in ('E','V','A','B')}
    jets={name:current.SlowJet(c,{k:current.LogUpper.constant(c,2+index+2*k[0]+k[1]) for k in current.ORDERS}) for index,name in enumerate(f)}
    groups=current.leading_component_caps(jets);density=dict(m=f['B'],h=f['E']*f['A'],k=f['E']*(f['V']*f['A']+f['B']),
        e=2*f['V']*f['B']-f['E']**2*f['A'],p=f['E']**2*f['A'])
    count=0
    for derivative,rows in groups.items():
        for key,terms in rows.items():
            expression=0
            for term in terms:
                product=sy.Integer(term['sign']*term['coefficient'])
                require(sum(v['source'] in ('A','B') for v in term['factors'])==1,'Every leading term must preserve original zero-mean primitive')
                for factor in term['factors']:
                    product*=sy.diff(f[factor['source']],y,factor['ordinary_y_order'],Z,factor['ordinary_Z_order'])
                expression+=product;count+=1
            expected=sy.diff(density[key],y) if derivative=='C0' else sy.diff(density[key],y,Z)
            require(sy.expand(expected-expression)==0,'Exact signed leading ordinary derivative term identity failed')
    require(count==56,'Expected complete first/mixed leading density product terms')
    return dict(passed=True,independent_exact_signed_slow_product_terms=count,
        all_ten_density_derivative_identites_checked=True,ordinary_derivatives_not_taylor_coefficients=True)


def fixture_references():
    # Reuse the accepted independent original scalar loop fixtures with this
    # new cap backend. No accepted module or runtime function is replaced.
    namespace=dict(preceding_checks.primitive_references.__globals__);namespace['current']=current
    fixture=types.FunctionType(preceding_checks.primitive_references.__code__,namespace,'active_kappa_scalar_primitive_references')
    return fixture()


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_Z_query_count']==2,'Completed original active-kappa stage required')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed active-kappa original source: '+name)
    require(not saved['actual_controls_or_terminal_closure_installed'] and not any(saved.get(k) for k in packets.OPEN),
        'Mixed bound improvement cannot admit controls/global field/recursion')
    kernels=weighted_kernel_references();components=leading_component_identity_check();fixtures=fixture_references()
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        targets=current.current.NativeRcParameterTargets(current.current.preceding.NativeRcC1Histories(current.current.preceding.preceding.NativeO2C1Histories(current.current.preceding.preceding.make_middle_owner(bridge))))
        owner=current.NativeActiveKappaMixedConditioning(current.previous.NativeSignedAveraging(targets))
        require(saved['exact_correlated_periodic_kernel_theorem']==packets.encode(owner.theorem),'New exact correlation theorem changed')
        direct=json.loads((HERE/current.current.NAME).read_bytes())['actual_original_Rc_parameter_target_records']
        rows=inherited=quiet=pressure=terms=weighted=targets_count=active=0;regions={};direct_comparisons={}
        for region,record in saved['actual_native_active_kappa_mixed_conditioning_records'].items():
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
                            and proof['joint_weighted_kernel_proof']['no_q_or_chi_division_used'],
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
            print('Native active-kappa and joint weighted kernels checked:',region,flush=True)
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,native_Z_queries_checked=2,
        original_true_cells_checked=24,original_charts_checked=17,actual_C0_Z_transport_rows_checked=rows,
        actual_inherited_rows_checked=inherited,exact_quiet_rows_checked=quiet,preserved_pressure_memory_rows_checked=pressure,
        active_correlated_primitive_cells_checked=active,original_leading_slow_source_terms_checked=terms,
        first_bridge_weighted_slow_terms_checked=weighted,normalized_target_rows_checked=targets_count,
        independent_weighted_kernel_references=kernels,independent_exact_leading_components=components,
        independent_original_scalar_primitive_references=fixtures,exact_correlated_periodic_kernel_theorem=owner.theorem,
        same_source_direct_floor_comparisons=direct_comparisons,regions=regions,
        actual_controls_or_terminal_closure_installed=False,point_inverse_or_higher_jets_installed=False,
        one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes,scope=saved['scope'])
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Native active-kappa and joint weighted kernels PASS',rows,inherited,terms,weighted,targets_count,flush=True)
    return result


if __name__=='__main__':run()
