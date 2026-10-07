"""Independent periodic kernel/implicit loop fixtures and native route checks."""
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_native_periodic_mixed_conditioning as current
import lei_ren_part1_paper_compliant_current_native_signed_averaging_check as preceding_checks
from lei_ren_part1_paper_compliant_current_generic_shear_loop import GenericLoopScales,GenericShearLoop

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
packets,ep=current.packets,current.ep;require,same=preceding_checks.require,preceding_checks.same


def kernel_references():
    p=mp.mp.clone();p.dps=70;count=0;point_count=0
    def d(u,psi):
        h=p.sqrt(1+u*u);r=u/h;D=1-2*r*p.cos(psi)+r*r
        return (p.cos(psi)-r)/(h*D)
    def pieces(end):
        return sorted(set([p.mpf(0),end]+[v for v in (p.mpf('.02'),p.mpf('.1'),p.pi-p.mpf('.1'),p.pi,p.pi+p.mpf('.1')) if 0<v<end]))
    for u in (p.mpf(0),p.mpf('.4'),p.mpf('-.4'),p.mpf(2),p.mpf(-2),p.mpf(20),p.mpf(-20)):
        for psi in (p.mpf('.73'),p.mpf('4.9')):
            du=lambda x:p.diff(lambda uu:d(uu,x),u)
            duu=lambda x:p.diff(lambda uu:d(uu,x),u,2)
            values=(p.quad(lambda x:d(u,x),pieces(psi)),p.quad(du,pieces(psi)),p.quad(duu,pieces(psi)),
                p.quad(lambda x:d(u,x)**2,pieces(psi)),p.quad(lambda x:2*d(u,x)*du(x),pieces(psi)),
                p.quad(lambda x:2*(du(x)**2+d(u,x)*duu(x)),pieces(psi)))
            for value,cap in zip(values,(4*p.pi,20*p.pi,100*p.pi,p.pi,50*p.pi,600*p.pi)):
                require(abs(value)<=cap+p.mpf('1e-45'),'Independent uniform partial-angle P/H derivative bound failed');count+=1
            require(abs(du(psi))<=6 and abs(duu(psi))<=50,'Independent point kernel derivative cap failed');point_count+=2
            # Independent direct theta/psi differentiation tests the coupled
            # normalized curvature cap, including signed u and shifted t0.
            for t0,q in ((p.mpf('-.8'),p.mpf('.5')),(p.mpf('2.5'),p.mpf('.03'))):
                t=t0+2*q*d(u,psi);tpsi=2*q*p.diff(lambda x:d(u,x),psi);chi=abs(u/q)
                cap=4*q*(1+chi*(1+abs(t0)))**p.mpf('1.5')
                require(abs(tpsi)/(1+t*t)**2<=cap and 2*abs(t*tpsi)/(1+t*t)**3<=cap,
                    'Independent correlated normalized inverse curvature failed');point_count+=2
    for u in (p.mpf('.4'),p.mpf(-2)):
        total=p.quad(lambda x:d(u,x)**2,pieces(2*p.pi))
        require(abs(total-p.pi)<p.mpf('1e-45'),'Original full-period H normalization lost');count+=1
    return dict(passed=True,independent_partial_angle_integrated_P_H_derivative_comparisons=count,
        independent_point_kernel_and_normalized_curvature_comparisons=point_count,
        signed_u_cases=['0','.4','-.4','2','-2','20','-20'],
        references_use_direct_original_kernel_integrals_not_saved_caps=True)


def primitive_references():
    iv=MPIntervalContext();iv.dps=110;p=mp.mp.clone();p.dps=85;count=0;branches=[]
    for kind,p20 in (('body_positive_r','.5'),('body_negative_r','-.5'),('body_r_crossing','0'),('cutoff_transition','.5'),('flat','.5')):
        scales=GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.1',
            t0_abs_max='1',p1_abs_max='10',p2_abs_max='3',dps=85);p=scales.ctx
        if kind=='cutoff_transition':a0=2+scales.eta/2;ay=scales.eta/50;az=scales.eta/80;ayz=scales.eta/100;b0=p.mpf(0);by=bz=byz=p.mpf(0)
        elif kind=='flat':a0=p.mpf('3.1');ay=p.mpf('.01');az=p.mpf('.02');ayz=p.mpf('.001');b0=p.mpf('.1');by=p.mpf('.002');bz=p.mpf('-.002');byz=p.mpf('.0001')
        else:a0=p.mpf('.8');ay=p.mpf('.01');az=p.mpf('.02');ayz=p.mpf('.001');b0=p.mpf('-.1');by=p.mpf('.002');bz=p.mpf('-.002');byz=p.mpf('.0001')
        cv=lambda x:iv.mpf(p.nstr(p.mpf(x),95));Y=iv.mpf(('-.001','.001'));Z=iv.mpf(('-.001','.001'))
        bases=tuple(iv.mpf(0) for _ in range(5));ledger={};scalar=lambda value:current.prior.ScaledEnclosure(current.prior.FormalScale(bases),iv.mpf(value),ledger)
        polynomial=lambda v,vy,vz,vyz:{current.ZERO:scalar(cv(v)+cv(vy)*Y+cv(vz)*Z+cv(vyz)*Y*Z),
            current.DY:scalar(cv(vy)+cv(vyz)*Z),current.DZ:scalar(cv(vz)+cv(vyz)*Y),current.DYZ:scalar(cv(vyz))}
        a=polynomial(a0,ay,az,ayz);b=polynomial(b0,by,bz,byz);p2=polynomial(p.mpf(p20),p.mpf('.01'),p.mpf('.02'),p.mpf('.001'))
        logE=cv('.03')*Z+(1-cv(a0)-cv(az)*Z)*Y/2-cv(ay)*Y*Y/4-cv(ayz)*Y*Y*Z/4
        E0=scalar(iv.exp(logE));EZ=E0*(cv('.03')-cv(az)*Y/2-cv(ayz)*Y*Y/4)
        E={current.ZERO:E0,current.DY:(scalar(1)-a[current.ZERO])*E0*iv.mpf('.5'),current.DZ:EZ,
            current.DYZ:(scalar(1)-a[current.ZERO])*EZ*iv.mpf('.5')-a[current.DZ]*E0*iv.mpf('.5')}
        zero=scalar(0);logamin=iv.ln(iv.mpf('.7'))
        def product(left,right):return {(j,k):sum((left[(i,ell)]*right[(j-i,k-ell)] for i in range(j+1) for ell in range(k+1)),zero) for j,k in current.ORDERS}
        def quotient(num,den):
            result={current.ZERO:num[current.ZERO].positive_divide(den[current.ZERO],logamin)}
            for k in (current.DY,current.DZ):result[k]=(num[k]-result[current.ZERO]*den[k]).positive_divide(den[current.ZERO],logamin)
            result[current.DYZ]=(num[current.DYZ]-result[current.DY]*den[current.DZ]-result[current.DZ]*den[current.DY]-result[current.ZERO]*den[current.DYZ]).positive_divide(den[current.ZERO],logamin)
            return result
        t0=quotient({k:-v for k,v in b.items()},a);ratio=quotient(product(b,b),a)
        Delta={k:a[k]+ratio[k]-(2 if k==current.ZERO else 0) for k in current.ORDERS}
        roots=dict(a=a,b=b,p2=p2,E=E,t0=t0,kappa_minus2=Delta)
        logeta=iv.ln(cv(scales.eta));logd=iv.ln(cv(scales.d_star))
        primitive=current.serial.whole_period_C1(roots,logeta,logamin,logd)
        bounded=current.primitive_caps(iv,roots,dict(log_actual_a_positive_lower=packets.encode(logamin)),logeta,logd,
            {1:current.LogUpper.constant(iv,32),2:current.LogUpper.constant(iv,1792)},primitive)
        def evaluate(y,z,phi):
            aa=a0+ay*y+az*z+ayz*y*z;bb=b0+by*y+bz*z+byz*y*z
            pp=p.mpf(p20)+p.mpf('.01')*y+p.mpf('.02')*z+p.mpf('.001')*y*z
            ee=p.exp(p.mpf('.03')*z+(1-a0-az*z)*y/2-ay*y*y/4-ayz*y*y*z/4)
            loop=GenericShearLoop(scales,a=aa,b=bb,p1=8,p2=pp,Utheta=ee)
            value=loop.evaluate(phi)
            return loop,value,ee
        h=p.mpf('1e-6')
        for phi in (p.mpf('.23'),p.mpf('.67')):
            loop,center,ee=evaluate(0,0,phi)
            if kind=='flat':require(loop.q==0 and ee>0,'Flat fixture must preserve nonzero underlying field')
            else:require(loop.q>0,'Body/transition independent active fixture lost')
            yp=evaluate(h,0,phi)[1];ym=evaluate(-h,0,phi)[1];zp=evaluate(0,h,phi)[1];zm=evaluate(0,-h,phi)[1]
            pp=evaluate(h,h,phi)[1];pm=evaluate(h,-h,phi)[1];mpv=evaluate(-h,h,phi)[1];mm=evaluate(-h,-h,phi)[1]
            for key in ('A','B'):
                values={current.ZERO:center[key],current.DY:(yp[key]-ym[key])/(2*h),current.DZ:(zp[key]-zm[key])/(2*h),
                    current.DYZ:(pp[key]-pm[key]-mpv[key]+mm[key])/(4*h*h)}
                for order,value in values.items():
                    cap=bounded[key][order]
                    require(abs(value)<p.mpf('1e-45') if cap.log is None else p.log(abs(value)+p.mpf('1e-75'))<=p.mpf(str(ep(cap.log)[1])),
                        'Independent fixed-phi primitive derivative outside new cap: '+kind+' '+key+' '+str(order));count+=1
            branches.append(kind)
    return dict(passed=True,independent_scalar_original_loop_C0_y_Z_yZ_comparisons=count,
        native_fields_not_defined_by_fixture_caps=True,body_both_signed_r_and_r_crossing_and_transition_flat_tested=True,
        exact_flat_primitive_jets_with_nonzero_underlying_field_preserved=True,branches=branches)


@current.native.inlet.source_precision
def run(bridge=None):
    began=time.monotonic();saved=json.loads((HERE/current.NAME).read_bytes())
    require(saved[current.GATE] and saved['native_Z_query_count']==2,'Native mixed conditioned route stage required')
    for name,digest in saved['input_hashes'].items():require(sha(name)==digest,'Changed original conditioned source: '+name)
    require(not saved['actual_controls_or_terminal_closure_installed'] and not any(saved.get(k) for k in packets.OPEN),
        'Periodic slow caps cannot admit controls/global field/recursion')
    kernels=kernel_references();fixtures=primitive_references();exact=preceding_checks.exact_mixed_source_theorem()
    if bridge is None:bridge,_=current.native.inlet.native_bridge_owner()
    with current.native.inlet.CheckedSourceRuntime():
        targets=current.current.NativeRcParameterTargets(current.current.preceding.NativeRcC1Histories(current.current.preceding.preceding.NativeO2C1Histories(current.current.preceding.preceding.make_middle_owner(bridge))))
        owner=current.NativePeriodicMixedConditioning(current.previous.NativeSignedAveraging(targets))
        require(saved['exact_correlated_periodic_kernel_theorem']==packets.encode(owner.theorem),'Correlated kernel/source theorem changed')
        direct=json.loads((HERE/current.current.NAME).read_bytes())['actual_original_Rc_parameter_target_records']
        rows=inherited=quiet=pressure=mixed_rows=budget_rows=target_rows=0;regions={};direct_comparison={}
        for name,record in saved['actual_native_periodic_mixed_conditioning_records'].items():
            got=owner.route(packets.interval(owner.ctx,record['Z_box']))
            require(packets.encode(got['record'])==record,'Conditioned native route changed: '+name)
            require(len(got['cells'])==24 and record['original_charts']==17,'Full original24-cell/17-chart route required')
            incoming={key:owner.coordinates.scalar(0) for key in current.RATES};incomingZ=dict(incoming)
            for cell in got['cells']:
                label=cell['record']['label'];chart=cell['record']['chart']
                if label!='initial_flat_collar':
                    proof=cell['record']['original_source']['native_normalized_active_velocity_correlations']
                    require(proof['original_velocity_outside_q_support_not_clipped_or_zeroed']
                        and proof['restrictions_used_only_in_modulation_densities_and_primitive_jets'],
                        'Active velocity bounds cannot change the background field outside q support')
                    mixed_rows+=8
                for key in current.RATES:
                    decay=cell['factors'][key]['decay']
                    same(cell['incoming'][key],incoming[key],'Native conditioning reset incoming C0 memory')
                    same(cell['incoming_Z'][key],incomingZ[key],'Native conditioning reset incoming Z memory')
                    same(cell['cumulative'][key],decay*incoming[key]+cell['values'][key],'Conditioned C0 route transport mismatch')
                    same(cell['cumulative_Z'][key],decay*incomingZ[key]+cell['Z_derivatives'][key],'Conditioned Z route transport mismatch')
                    rows+=4;inherited+=2
                    if chart=='O3_power':
                        require(cell['values'][key].zero and cell['Z_derivatives'][key].zero,'Conditioning introduced nonzero quiet modulation');quiet+=2
                        if key=='p':
                            same(cell['cumulative'][key],incoming[key],'Conditioning lost quiet C0 pressure memory')
                            same(cell['cumulative_Z'][key],incomingZ[key],'Conditioning lost quiet Z pressure memory');pressure+=2
                    if label!='initial_flat_collar':
                        require(cell['record']['per_density_budget'][key]['same_original_uniform_coefficient_cover_intersection'],
                            'Old/new cap intersection requires the same original source function')
                        budget_rows+=10
                incoming,incomingZ=cell['cumulative'],cell['cumulative_Z']
            require(all(v['new_certificate_strictly_sharper'] for v in record['same_source_averaging_bound_comparison'].values()),
                'Expected current conditioned averaging certificates to improve all five target rows')
            require(set(record['first_bridge_native_weighted_target_budget_components'])=={
                'right_endpoint','decayed_left_endpoint','slow_variable_integral','kernel_derivative_integral','quadratic_remainder_integral'},
                'Separate first-bridge endpoint/slow/kernel/quadratic budgets required')
            floor_comparison={}
            for key in current.repair.ROWS:
                require(got['targets']['values'][key][-1].zero and got['targets']['Z_derivatives'][key][-1].zero,'N^-2 target certificate lost')
                target_rows+=4
                new=current.previous.read_cap(owner.ctx,record['actual_uniform_N_scaled_repair_C1_caps']['transformed_N_scaled_target_C1_caps'][key])
                old=current.previous.read_cap(owner.ctx,direct[name]['actual_uniform_N_scaled_repair_C1_caps']['transformed_N_scaled_target_C1_caps'][key])
                floor_comparison[key]=dict(conditioned_average_strictly_sharper_than_prior_direct_uniform_cap=(new.log is None or (old.log is not None and ep(new.log)[1]<ep(old.log)[0])),
                    prior_direct_N_scaled_C1_cap=old.record(),conditioned_average_N_scaled_C1_cap_at_floor=new.record(),
                    bound_comparison_not_true_error_measurement=True)
            direct_comparison[name]=floor_comparison
            regions[name]=dict(all_five_target_bounds_improved_against_prior_averaging=True,
                dominant_target_cells=record['dominant_target_bound_cells'],first_bridge_dominant_components=record['first_bridge_dominant_budget_components'],
                new_uniform_N_scaled_target_C1_cap=record['actual_uniform_N_scaled_repair_C1_caps']['whole_target_C1_cap'])
            print('Native periodic mixed conditioning checked:',name,flush=True)
    hashes=dict(owner.service.hashes);hashes[current.NAME]=sha(current.NAME);hashes[Path(__file__).name]=sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=owner.family,**{current.GATE:True},native_Z_queries_checked=2,
        original_true_cells_checked=24,original_charts_checked=17,actual_C0_Z_contribution_and_transport_rows_checked=rows,
        actual_inherited_memory_rows_checked=inherited,exact_quiet_rows_checked=quiet,preserved_pressure_memory_rows_checked=pressure,
        native_primitive_four_slow_rows_checked=mixed_rows,separate_density_budget_rows_checked=budget_rows,
        normalized_target_order_rows_checked=target_rows,
        independent_original_kernel_references=kernels,independent_original_scalar_primitive_references=fixtures,
        independent_mixed_source_identity_check=exact,exact_correlated_periodic_kernel_theorem=owner.theorem,
        prior_direct_and_conditioned_floor_bound_comparisons=direct_comparison,regions=regions,
        actual_controls_or_terminal_closure_installed=False,point_inverse_or_higher_jets_installed=False,
        one_global_finite_N_or_control_field_admitted=False,global_inlet_to_Rc_histories_admitted=False,**dict.fromkeys(packets.OPEN,False),
        execution_seconds=time.monotonic()-began,input_hashes=hashes,scope=saved['scope'])
    (HERE/current.RECEIPT).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Native correlated periodic mixed conditioning PASS',rows,inherited,quiet,pressure,mixed_rows,budget_rows,target_rows,flush=True)
    return result


if __name__=='__main__':run()
