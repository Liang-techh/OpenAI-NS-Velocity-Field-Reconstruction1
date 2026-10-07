"""Current-source attachment, branch conditioning and independent loop fixtures."""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_shear_loop_jet_bounds as source
from lei_ren_part1_paper_compliant_current_generic_shear_loop import GenericLoopScales,GenericShearLoop

packets,bounds=source.packets,source.bounds


def independent_fixtures():
    """Only modest artificial sources are evaluated, never actual source logs."""
    c=mp.mp.clone();c.dps=90
    iv=MPIntervalContext();iv.dps=110;ep=packets.recovery.endpoints
    total=0;branches=[]
    for kind in ('active','cutoff_transition','flat'):
        scales=GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.1',
            t0_abs_max='1',p1_abs_max='10',p2_abs_max='3',dps=90)
        c=scales.ctx
        if kind=='active':a0=c.mpf('.8');ay=c.mpf('.02');az=c.mpf('.015');b0=c.mpf('-.1')
        elif kind=='cutoff_transition':a0=2+scales.eta/2;ay=scales.eta/50;az=scales.eta/80;b0=c.mpf(0)
        else:a0=c.mpf(3);ay=c.mpf('.02');az=c.mpf('.01');b0=c.mpf('.1')
        def fields(y,Z):
            a=a0+ay*y+az*Z+c.mpf('.000001')*y*Z+c.mpf('.000001')*y*y
            b=b0+c.mpf('.001')*y-c.mpf('.002')*Z
            p2=c.mpf('.5')+c.mpf('.01')*y+c.mpf('.02')*Z
            E=1+c.mpf('.02')*y+c.mpf('.03')*Z+c.mpf('.002')*y*y
            return a,b,p2,E
        def jet(values):
            return source.BoundJet(iv,{k:bounds.LogUpper.constant(iv,iv.mpf(str(v))) for k,v in zip(source.ORDERS,values)})
        # These direct polynomial caps hold on |y|,|Z|<=.01.
        aj=jet((a0+c.mpf('.01')*(abs(ay)+abs(az))+c.mpf('.000000001'),
            abs(ay)+c.mpf('.00000003'),abs(az)+c.mpf('.00000001'),'.000002','.000001'))
        bj=jet((abs(b0)+c.mpf('.00003'),'.001','.002',0,0))
        pj=jet(('.501','.01','.02',0,0));Ej=jet(('1.001','.02004','.03','.004',0))
        logamin=iv.ln(iv.mpf('.7'))
        t0j=bj*aj.inverse_positive(logamin)
        result=source.loop_bounds(iv,dict(a=aj,b=bj,p2=pj,t0=t0j,E=Ej),
            log_a_min=logamin,log_d=iv.ln(iv.mpf(str(scales.d_star))),
            log_eta=iv.ln(iv.mpf(str(scales.eta))),log_q_star=iv.ln(iv.mpf(str(scales.q_star))),
            sigma_caps={1:bounds.LogUpper.constant(iv,32),2:bounds.LogUpper.constant(iv,1792)})
        def evaluate(y,Z,phi):
            a,b,p2,E=fields(y,Z)
            loop=GenericShearLoop(scales,a=a,b=b,p1=8,p2=p2,Utheta=E)
            # Public scalar constructor/inverter is independent of log bounds.
            primitives=loop.evaluate(phi)
            return loop,dict(psi=loop.angle_at_phase(phi),A=primitives['A'],B=primitives['B'])
        h=c.mpf('.00001')
        for phi in (c.mpf('.23'),c.mpf('.67')):
            loop,center=evaluate(0,0,phi)
            actualbranch='flat' if loop.q==0 else 'active'
            if kind=='flat' and actualbranch!='flat':raise ArithmeticError('Independent flat branch fixture lost')
            if kind!='flat' and actualbranch!='active':raise ArithmeticError('Independent active branch fixture lost')
            branches.append(kind)
            yp=evaluate(h,0,phi)[1];ym=evaluate(-h,0,phi)[1]
            zp=evaluate(0,h,phi)[1];zm=evaluate(0,-h,phi)[1]
            pp=evaluate(h,h,phi)[1];pm=evaluate(h,-h,phi)[1]
            mpv=evaluate(-h,h,phi)[1];mm=evaluate(-h,-h,phi)[1]
            for key in ('psi','A','B'):
                actual={source.ZERO:center[key],(1,0):(yp[key]-ym[key])/(2*h),
                    (0,1):(zp[key]-zm[key])/(2*h),
                    (2,0):(yp[key]-2*center[key]+ym[key])/(h*h),
                    (1,1):(pp[key]-pm[key]-mpv[key]+mm[key])/(4*h*h)}
                for order,value in actual.items():
                    record=result['slow_phase_held_log_bounds'][key]['y%d_Z%d'%order]
                    cap=c.exp(ep(record['log_absolute_upper'])[1])
                    if abs(value)>cap+c.mpf('1e-12'):
                        raise ArithmeticError('Independent phase-held loop derivative exceeds cap: '+kind+' '+key+str(order))
                    total+=1
            lam=loop.a*(1+loop.direction(loop.angle_at_phase(phi))**2)/(2*c.pi*loop.v)
            lower=c.exp(ep(result['positive_conditioning']['log_lambda_positive_lower'])[0])
            if lam<lower*(1-c.mpf('1e-75')):raise ArithmeticError('Correlated phase lower exceeds independently evaluated lambda')
            total+=1
    return dict(comparisons=total,branches=sorted(set(branches)),
        finite_difference_step='1e-5, artificial sources only; absolute1e-12 allowance',
        actual_source_factors_or_scales_not_evaluated=True)


def run():
    manifest=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in manifest['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed current loop derivative source: '+name)
    provider=source.CurrentLoopJetBounds();c=provider.ctx;ep=packets.recovery.endpoints
    charts=manifest['current_actual_loop_jet_log_bounds_by_chart']
    if set(charts)!=set(provider.scales['whole_source_chart_inventory']) or len(charts)!=17:
        raise ArithmeticError('Whole current modification source cover incomplete')
    count=0
    for chart,record in charts.items():
        actual=provider.chart(chart)
        if packets.encode(actual)!=record:raise ValueError('Actual source loop log bound differs: '+chart)
        for rows in actual['slow_phase_held_log_bounds'].values():
            if set(rows)!={'y%d_Z%d'%k for k in source.ORDERS}:raise ArithmeticError('Unsupported loop derivative order exported')
            for bound in rows.values():
                if not bound['exact_zero'] and not all(mp.isfinite(v) for v in ep(bound['log_absolute_upper'])):
                    raise ArithmeticError('Nonfinite actual loop derivative log cap')
                count+=1
        for bound in actual['phase_and_first_slow_phase_log_bounds'].values():
            if not bound['exact_zero'] and not all(mp.isfinite(v) for v in ep(bound['log_absolute_upper'])):
                raise ArithmeticError('Nonfinite actual phase derivative log cap')
            count+=1
        for key in ('log_lambda_positive_lower','log_Poisson_D_positive_lower'):
            if not all(mp.isfinite(v) for v in ep(actual['positive_conditioning'][key])):
                raise ArithmeticError('True positive logarithmic conditioning lower missing')
        if not actual['q_flat_branch']['A_and_B_all_jets_exact_zero']:
            raise ArithmeticError('Flat edge must retain exact zero primitive jets')
    if any(manifest[k] for k in source.OPEN) or manifest['signed_current_point_loop_or_inverse_jets_installed']:
        raise ArithmeticError('Derivative log bounds promoted missing installed physical loop stages')
    if manifest['Z2_y2Z_or_full_mixed4_loop_bounds_certified']:
        raise ArithmeticError('Unsupported high derivative rows admitted')
    if packets.encode(provider.theorem)!=manifest['exact_inverse_primitive_theorem']:
        raise ValueError('Current phase-held recurrence theorem differs')
    # Reject eta above the actual Section11 support range.
    fixture=source.BoundJet.constant(c,1)
    try:source.loop_bounds(c,{k:fixture for k in ('a','b','p2','t0','E')},
        log_a_min=c.mpf(0),log_d=c.mpf(0),log_eta=c.mpf(0),log_q_star=c.mpf(0),sigma_caps=provider.sigma_caps)
    except ValueError:pass
    else:raise ArithmeticError('Invalid eta support cap admitted')
    independent=independent_fixtures()
    result=dict(all_passed=True,source_family=provider.family,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),
        phase_held_loop_primitive_derivative_bounds_certified=True,
        source_cover_charts=17,finite_current_source_derivative_log_caps=count,
        certified_phase_held_orders=[list(k) for k in source.ORDERS],
        original_fixed_step_derivative_orders=[1,2],
        exact_implicit_and_primitive_chain_rules=len(provider.theorem['exact_chain_rule_identities']),
        independent_scalar_loop_fixture=independent,invalid_eta_guard_passed=True,
        B_bounds_in_original_common_Pstar_unit=True,
        signed_current_point_loop_or_inverse_jets_installed=False,
        Z2_y2Z_or_full_mixed4_loop_bounds_certified=False,
        actual_eta_amplitude_radius_width_inverse_or_log_caps_not_exponentiated=True,
        source_graph_ancestor_constructors_called=False,
        input_hashes={**manifest['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Actual phase-held inverse/primitive log bounds PASS:17 charts; independent active/transition/flat fixtures PASS',flush=True)
    return result


if __name__=='__main__':run()
