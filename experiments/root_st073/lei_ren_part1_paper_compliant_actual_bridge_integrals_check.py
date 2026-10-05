"""Focused independent checks of finite-width signed bridge/own histories.

Checks source families, endpoint inheritance, exact finite-width radius
kernels, all logarithmic cap proofs, and direct original scalar integrals.
Moderate fixture parameters exercise the algebra; they are not selected
production parameters and do not establish global stress admissibility.
"""
import hashlib
import json
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode
from lei_ren_part1_paper_compliant_actual_bridge_integrals import (
    radius_kernels, volterra_primitive_kernels, own_moment_feedback, pulse_control_masses,
)

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_actual_bridge_integrals'


def sha(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def same_rows(left,right):
    return all(a.get('lower_exact_mpf_tuple')==b.get('lower_exact_mpf_tuple')
               and a.get('upper_exact_mpf_tuple')==b.get('upper_exact_mpf_tuple')
               for a,b in zip(left,right)) and len(left)==len(right)


def run():
    with mp.workdps(240):
        c=MPIntervalContext();c.dps=240
        r=json.loads((HERE/(NAME+'.json')).read_bytes())
        hashes=dict(r['input_hashes'])
        for name,digest in hashes.items():
            if sha(name)!=digest:raise ValueError('Current bridge source changed: '+name)
        hashes[NAME+'.json']=sha(NAME+'.json')
        hashes[NAME+'_check.py']=sha(NAME+'_check.py')
        comparison=json.loads((HERE/'lei_ren_part1_paper_compliant_comparison_point_integrals.json').read_bytes())
        for key in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256'):
            if r[key]!=comparison[key]:raise ValueError('Bridge changes the defining family: '+key)
        if not (r['full_finite_width_signed_bridge_integral_enclosures_available']
                and r['actual_own_six_moment_feedback_enclosures_available']
                and r['known_direction_not_replaced_by_actual_inertial_stress']):
            raise ValueError('Actual known-direction source/own-history gate missing')
        for key in ('actual_bridge_mixed4_feedback_installed','R100_R110_actual_feedback_composed',
                    'full_implicit_leading_inputs_recomputed','global_completed_tensor_admissibility',
                    'exact_production_point_parameters_selected','temporal_recursion'):
            if r[key]:raise ValueError('Unimplemented downstream claim: '+key)
        bindings=r['original_control_AST_bindings']
        if not all(row['verified'] for row in bindings.values()):raise ValueError('Original controls not bound')
        cap=read_interval(c,r['width_enclosure_cap']);logh=read_interval(c,r['source_log_hb_enclosure'])
        if endpoints(logh)[1]>=endpoints(c.ln(cap))[0]:raise ValueError('Cap used as source width')
        for proof in r['source_log_product_cap_proofs']:
            if not proof['passed'] or endpoints(read_interval(c,proof['product_log_upper']))[1]>endpoints(c.ln(cap))[0]:
                raise ValueError('Unproved source-factor enclosure')
        joins=0;core_atoms=0;finite_rows=0;signed_terms=0
        fields=('actual_log_F_over_inlet_axial5','actual_phi_axial5','actual_raw_V_axial5',
                'actual_radial_Q_axial4','pressure_axis_axial5',
                'actual_pressure_increment_axial5_divided_by_R_F0_squared')
        for label,packets in r['packets'].items():
            core=packets['core']
            for name,values in core['actual_own_six_moments_axial5'].items():
                if not same_rows(values,core['actual_core_atom_inlet'][name]):
                    raise ValueError('Original core atom reset: '+name)
                expected=comparison['comparison_point_packets'][label]['micro'][0]['comparison_own_six_moments'][name]['signed_hb_power_axial_coefficients'][0][:6]
                # Fresh atoms are independently recomputed by the producer.
                # Overlap with accepted old atom enclosures is required; exact
                # equality is not inferred from differently rounded receipts.
                for a,b in zip(values,expected):
                    al,ah=endpoints(read_interval(c,a));bl,bh=endpoints(read_interval(c,b))
                    if max(al,bl)>min(ah,bh):raise ValueError('Fresh core atom leaves accepted source')
                    core_atoms+=1
            for left,right in (('first_exit','second_inlet'),('second_exit','macro_inlet')):
                for key in fields:
                    if not same_rows(packets[left][key],packets[right][key]):
                        raise ValueError('Original source join changed: '+left+'/'+key)
                    joins+=1
                for name in core['actual_own_six_moments_axial5']:
                    if not same_rows(packets[left]['actual_own_six_moments_axial5'][name],
                                     packets[right]['actual_own_six_moments_axial5'][name]):
                        raise ValueError('Own moment history reset at join: '+name)
                    joins+=1
            end=packets['R100']
            if not end['coordinate_R100_endpoint'] or endpoints(read_interval(c,end['radius_enclosure_only']))!=(mp.mpf(100),mp.mpf(100)):
                raise ValueError('Finite-width macro endpoint is not exactly R100')
            for p in packets.values():
                if not (p['same_original_core_atoms_used'] and p['known_comparison_direction_preserved']
                        and p['actual_moment_enclosures_from_prescribed_FV']
                        and p['nonlinear_angular_prefix_remainder_controlled']):
                    raise ValueError('Original field/own-history source not retained')
                if p['comparison_moments_substituted_for_actual'] or p['cap_used_as_field_value'] or p['actual_point_moment_history_recovered']:
                    raise ValueError('Invalid actual-field substitution')
                rows=[p[key] for key in fields]
                rows+=list(p['actual_own_six_moments_axial5'].values())
                for row in rows:
                    for value in row:
                        lo,hi=endpoints(read_interval(c,value))
                        if not mp.isfinite(lo) or not mp.isfinite(hi) or lo>hi:
                            raise ValueError('Invalid finite-width bridge row')
                        finite_rows+=1
                for term in p['signed_angular_integral_terms']+sum(p['signed_axial_integral_terms'].values(),[]):
                    if term['width_power'] not in (1,2) or not term['width_not_materialized']:
                        raise ValueError('Finite-width signed factor lost')
                    if endpoints(read_interval(c,term['nonlinear_prefix_error_weighted_norm_per_next_hb']))[0]<0:
                        raise ValueError('Negative nonlinear error bound')
                    signed_terms+=1
        whole=r['whole_axis_R100']
        if endpoints(read_interval(c,whole['Z']))!=(mp.mpf(-1),mp.mpf(1)):
            raise ValueError('Whole axial source domain missing')
        if not whole['signed_macro_own_moment_Volterra_feedback_integrated']:
            raise ValueError('Whole-axis actual own-history feedback missing')
        if not whole['coordinate_R100_endpoint'] or whole['R100_functional_join_to_existing_switch_installed']\
                or whole['actual_point_moment_history_recovered']:
            raise ValueError('Whole-axis endpoint scope overclaimed')
        if endpoints(read_interval(c,whole['radius_enclosure_only']))!=(mp.mpf(100),mp.mpf(100)):
            raise ValueError('Whole-axis radius is not the fixed original R100')
        for row in [whole[key] for key in fields]+list(whole['actual_own_six_moments_axial5'].values()):
            for value in row:
                lo,hi=endpoints(read_interval(c,value))
                if not mp.isfinite(lo) or not mp.isfinite(hi) or lo>hi:
                    raise ValueError('Invalid whole-axis finite-width field')
                finite_rows+=1

        # Independent finite-width scalar fixtures. Original angular prefix is
        # integrated directly; no producer modal helper or cap is its value.
        comparisons=0;max_miss=mp.mpf(0);tol=mp.mpf('1e-55')
        def enclose(value,interval):
            nonlocal comparisons,max_miss
            lo,hi=endpoints(interval)
            miss=max(lo-value,value-hi,mp.mpf(0))
            max_miss=max(max_miss,miss);comparisons+=1
            if miss>tol:raise ValueError('Independent source integral misses enclosure')
        with mp.workdps(90):
            for text in ('.1','.3','.8'):
                phase=mp.mpf(text)
                masses=pulse_control_masses(c,c.mpf(text))
                def sigma(t):
                    if t==0:return mp.mpf(0)
                    if t==1:return mp.mpf(1)
                    return 1/(1+mp.exp(1/t**2-1/(1-t)**2))
                w=mp.quad(lambda t:1-sigma(t),[0,phase/2,phase])
                enclose(w,masses[0]);enclose(phase-w,masses[1])
            h=mp.mpf('.002');Ra=mp.mpf('.7');R0=Ra*mp.exp(2*h);R1=mp.mpf(100)
            length=mp.log(R1/R0)
            kernels={}
            for power in (1,2):
                actual=radius_kernels(c,c.mpf(mp.nstr(R0,85)),c.mpf(100),c.mpf(mp.nstr(length,85)),power)
                kernels[power]=[mp.quad(lambda t,j=j:(R0*mp.exp(t))**power*mp.exp(-j*t),[0,length/2,length]) for j in range(3)]
                for value,box in zip(kernels[power],actual):enclose(value,box)
                for rate in (1,2):
                    boxes=volterra_primitive_kernels(c,c.mpf(mp.nstr(R0,85)),c.mpf(100),c.mpf(mp.nstr(length,85)),power,rate)
                    for j,box in enumerate(boxes):
                        a=power-j
                        primitive=lambda t,a=a:R0**power*(mp.expm1(a*t)/a if a else t)
                        value=mp.quad(lambda t:mp.exp(-rate*(length-t))*primitive(t),[0,length/2,length])
                        enclose(value,box)
            D=[mp.mpf('.9'),mp.mpf('.1'),mp.mpf('.05')]
            start=mp.mpf('.7')
            B=start+sum(abs(a)*k for a,k in zip(D,kernels[1]))/2
            def ell(t):
                parts=[R0*mp.expm1(t),R0*t,R0*(-mp.expm1(-t))]
                return -h*(start+sum(a*k for a,k in zip(D,parts))/2)
            Q=mp.mpf('1.1')
            for power,drive in ((1,['-.3','.2','-.1']),(1,['.08','0','0']),(2,['.01','-.02','.015'])):
                drive=list(map(mp.mpf,drive))
                signed=-Q*sum(a*k for a,k in zip(drive,kernels[power]))
                error=h*Q*B*mp.exp(h*B)*sum(abs(a)*k for a,k in zip(drive,kernels[power]))
                original=-mp.quad(lambda t:Q*mp.exp(ell(t))*(R0*mp.exp(t))**power
                    *sum(a*mp.exp(-j*t) for j,a in enumerate(drive)),[0,length/2,length])
                enclose(original,c.mpf([mp.nstr(signed-error,85),mp.nstr(signed+error,85)]))
            # Own history positive kernels retain all nonlinear cross products.
            phi0=mp.mpf('.8');V0=mp.mpf('.4');y=mp.mpf('1.3')
            initial=dict(H=mp.mpf('.6'),M=mp.mpf('.3'),K=mp.mpf('.2'),
                         A=mp.mpf('.18'),B=mp.mpf('.25'),C=mp.mpf('.5'))
            constant=lambda v:IntervalTaylor.constant(c,v,0)
            own=own_moment_feedback(constant(phi0),constant(V0),{n:constant(v) for n,v in initial.items()},
                c.mpf(mp.nstr(mp.exp(-y),85)),constant(c.mpf(['-.02','.02'])),constant(c.mpf(['-.03','.03'])))
            for name,rate,source in (
                ('H',2,lambda p,v:2*p),('M',1,lambda p,v:v),('K',2,lambda p,v:2*p*v),
                ('A',1,lambda p,v:v*v),('B',2,lambda p,v:p*p),('C',1,lambda p,v:p*p)):
                value=initial[name]*mp.exp(-rate*y)+mp.quad(lambda t:mp.exp(-rate*(y-t))
                    *source(phi0+mp.mpf('.02')*mp.sin(t),V0+mp.mpf('.03')*mp.cos(t)),[0,y])
                enclose(value,own['actual'][name][0])
        result=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],
            implicit_source_sha256=r['implicit_source_sha256'],datum_enclosure_sha256=r['datum_enclosure_sha256'],
            original_source_control_AST_bindings=len(bindings),current_input_hashes_checked=len(hashes),
            logarithmic_source_product_proofs=len(r['source_log_product_cap_proofs']),
            fresh_actual_core_atom_coefficients_checked=core_atoms,functional_endpoint_joins_checked=joins,
            finite_actual_field_and_own_moment_coefficients_checked=finite_rows,
            signed_finite_width_source_terms_checked=signed_terms,
            independent_original_scalar_integral_comparisons=comparisons,
            maximum_positive_enclosure_miss=max_miss,tolerance=tol,
            full_finite_width_signed_bridge_integral_enclosures_available=True,
            actual_own_six_moment_feedback_enclosures_available=True,
            actual_bridge_mixed4_feedback_installed=False,R100_R110_actual_feedback_composed=False,
            full_implicit_leading_inputs_recomputed=False,global_completed_tensor_admissibility=False,
            exact_production_point_parameters_selected=False,temporal_recursion=False,
            whole_real_axial_domain_R100_checked=True,
            all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    print('PASS actual finite-width signed bridge and own-history focused check',flush=True)
    return result


if __name__=='__main__':
    run()
