"""Focused source inequality, flat endpoints and complete mixed density check."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_correlated_qy_density_jets as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

ep=current.ep;C0,Y,Z,YZ=current.ORDERS


def source_inequality_calculus():
    sigma,complement,eta,a,S=s.symbols('sigma complement eta a S',positive=True)
    q2=(s.Rational(3,5)*complement+eta)/a
    ay=s.Rational(12,5)*sigma*complement*S
    qy2=(1+eta)**2*ay**2/(4*a**4*q2)
    local=s.Rational(36,25)*(1+eta)**2*s.Rational(5,3)/a**3*sigma**2*complement*S**2
    assert s.cancel(qy2/local-3*complement/(3*complement+5*eta))==0
    assert s.Rational(36,25)*4*s.Rational(5,3)*s.Rational(5,4)**3==s.Rational(75,4)
    x,c=s.symbols('x c',positive=True)
    assert s.simplify(s.diff(s.log(x**-3*s.exp(-c/x**2)),x)-(2*c/x**3-3/x))==0
    for count in (4,64,128):
        assert s.Rational(1,count)**2<s.Rational(1,3)
    yy=s.symbols('y',positive=True);odds=yy**-2-(1-yy)**-2
    assert s.cancel(s.diff(odds,yy)+2*(yy**-3+(1-yy)**-3))==0
    return dict(passed=True,exact_source_qy_squared_ratio_and_global_constant_identities=2,
        exact_logistic_odds_derivative_identity=1,paired_flat_product_log_derivative_identity=1,
        proof='positive eta gives ratio<1; eta<=1 and a>=4/5 give coefficient sqrt(75/4); paired endpoints monotone for count>=4',
        actual_positive_eta_not_set_to_zero=True)


def independent_source_qy_fixtures(owner):
    p=mp.mp.clone();p.dps=180;c=owner.ctx;cases=0;checks=0
    assert ep(current.source.ordered.interval(c,owner.source.owner.scales.logs['eta']))[1]<=0
    points=(s.Rational(1,128),s.Rational(1,64),s.Rational(1,10),s.Rational(1,4),
        s.Rational(1,2),s.Rational(3,4),s.Rational(9,10),s.Rational(63,64),s.Rational(127,128))
    for point in points:
        y=p.mpf(int(point.p))/int(point.q);odds=1/y**2-1/(1-y)**2
        # Independent stable logistic pair retains the tiny complementary
        # factor even when adding it to one would round it away.
        if odds>=0:
            t=p.exp(-odds);sigma=t/(1+t);complement=1/(1+t)
        else:
            t=p.exp(odds);sigma=1/(1+t);complement=t/(1+t)
        S=y**-3+(1-y)**-3;ay=p.mpf(12)/5*sigma*complement*S;a=p.mpf(4)/5+p.mpf(6)/5*sigma
        index=min(63,int(s.floor(point*64)))
        logcap,unused=current.correlated_qy_log_cap(c,64,index)
        for eta_text in ('.5','.01','1e-30','1e-5000'):
            eta=p.mpf(eta_text);q=p.sqrt((p.mpf(3)/5*complement+eta)/a)
            actual=-(1+eta)*ay/(2*a*a*q)
            assert actual<=0 and p.log(abs(actual))<=p.mpf(ep(logcap)[1])+p.mpf('1e-120')
            cases+=1
        if point not in (s.Rational(1,64),s.Rational(63,64)):
            logsig,logcomp=current.original_log_sigma_complement(c,point)
            saved.contains(logsig,p.log(sigma),p.mpf('1e-120'))
            saved.contains(logcomp,p.log(complement),p.mpf('1e-120'));checks+=2
    return dict(passed=True,independent_original_positive_eta_qy_cases=cases,
        independent_log_sigma_and_complement_comparisons=checks,
        reflected_flat_cells_and_positive_eta_dominant_or_cutoff_dominant_regimes_included=True,
        finite_parameter_diagnostics_not_replacements_of_original_scale=True)


def independent_mixed_density_fixtures(owner):
    frame=owner.source.frame(64,0,branch='regular');a=frame.roots['q'].atlas;c=owner.ctx
    p=mp.mp.clone();p.dps=160;cases=0;comparisons=0
    for N in (160,257):
        for ztext in ('-.37','.37'):
            z0=p.mpf(ztext);y0=p.mpf('.41')
            for Atext in ('-.7','0','1e-50','.6'):
                rows=dict(E=('1.1','.2','-.3','.04'),V=('-.21','.07','.4','-.05'),
                    A=(Atext,'.17','-.13','.11'),B=('.53','-.09','.23','-.08'))
                jets={name:current.MixedJet(a,{order:a.scalar(c.mpf(value)) for order,value in zip(current.ORDERS,values,strict=True)})
                    for name,values in rows.items()}
                got=current.mixed_changed_densities(a,**jets,N=N)
                def function(name,y,z):
                    v,dy,dz,dyz=map(p.mpf,rows[name]);return v+dy*(y-y0)+dz*(z-z0)+dyz*(y-y0)*(z-z0)
                def densities(E,V):return dict(m=V,h=E,k=E*V,e=V*V-E*E/2,p=E*E/2)
                def changed(key,y,z):
                    E,V,A,B=(function(name,y,z) for name in ('E','V','A','B'))
                    return densities(E*p.exp(A/N),V+B/N)[key]-densities(E,V)[key]
                for key,jet in got.items():
                    targets={C0:changed(key,y0,z0),Y:p.diff(lambda y:changed(key,y,z0),y0),
                        Z:p.diff(lambda z:changed(key,y0,z),z0),
                        YZ:p.diff(lambda y:p.diff(lambda z:changed(key,y,z),z0),y0)}
                    for order,target in targets.items():
                        saved.contains(jet[order].finite_interval(),target,p.mpf('1e-120'));comparisons+=1
                cases+=1
    return dict(passed=True,independent_exponential_history_C0_y_Z_yZ_comparisons=comparisons,
        compatible_mixed_source_function_fixtures=cases,nonzero_mixed_input_jets_and_all_quadratic_cross_terms=True,
        zero_and_tiny_A_stable_expm1_included=True,ordinary_mixed_derivatives_not_derivatives_of_interval_caps=True)


def actual_source_phase_density_contract(owner,manifest):
    c=owner.ctx;rows=manifest['whole_original_O2_correlated_qy_mixed_density_cells']
    assert len(rows)==64 and manifest['explicit_candidate_N']==160
    previous=json.loads(gzip.decompress((current.HERE/current.phase.NAME).read_bytes()))
    assert previous['source_family']==owner.family
    identical=0;native=0;frames=0;rejected=0
    for index,row in enumerate(rows):
        assert set(row)==set(current.integrals.BRANCHES)
        for branch,item in row.items():
            frame=owner.source.frame(64,index,branch=branch);a=frame.roots['q'].atlas;record=item['source']
            assert record['source_family']==owner.family and record['source_level']==64 and record['source_index']==index
            assert record['correlated_qy_source_bound_installed'] and record['actual_q_y_function_reenclosed_not_replaced_or_selected']
            assert record['original_positive_q_C0_and_eta_unchanged'] and record['original_nu_y_exact_source_cancellation_retained']
            assert record['exact_y_cell']==[str(s.Rational(index,64)),str(s.Rational(index+1,64))]
            logcap,proof=current.correlated_qy_log_cap(c,64,index)
            assert record['correlated_qy_source_proof']['recipe']==proof['recipe']
            qy=current.integrals.restore_value(a,record['correlated_qy_native_source'])
            assert qy.scale.powers==(0,0,0,0,0) and ep(qy.coefficient)[1]==0 and ep(qy.coefficient)[0]<0
            assert frame.gamma_q.scale.powers==(0,0,0,-1,0)
            assert frame.roots['q'][C0] is frame.kernel.q and frame.roots['q'][Z].zero and frame.roots['q'][YZ].zero
            saved.contains(qy.scale.offset,logcap,mp.mpf('1e-200'))
            # Only source-y information changed. Original value/Z A/B must
            # agree exactly with the accepted whole true-graph archive.
            old=previous['whole_original_O2_predicate_mixed_cells'][index][branch]['actual_whole_phase_mixed']['complete_original_A_B_mixed']
            actual=item['actual_whole_phase_mixed']['complete_original_A_B_mixed']
            for name in ('A','B'):
                for order in (C0,Z):
                    oldrange=current.integrals.export_range(a,current.integrals.restore_value(a,old[name][str(order)]),derivative_unit=order==Z)
                    newrange=current.integrals.export_range(a,current.integrals.restore_value(a,actual[name][str(order)]),derivative_unit=order==Z)
                    assert ep(oldrange)==ep(newrange),(index,branch,name,order);identical+=1
            assert item['all_original_density_nonlinear_cross_terms_retained']
            assert item['slow_y_and_yZ_at_fixed_true_phase_not_full_rapid_spatial_derivative']
            for name,jet in item['actual_five_changed_density_C0_y_Z_yZ'].items():
                assert name in current.integrals.KEYS and set(jet)==set(map(str,current.ORDERS))
                for value in jet.values():
                    assert value['encloses_original_source_function'] and not value['point_value_selected']
                    restored=current.integrals.restore_value(a,value)
                    assert restored.scale.bases is a.bases and restored.ledger is a.ledger;native+=1
            assert owner.source.frame(64,index,branch=branch) is frame;frames+=1
    frame=owner.source.frame(64,0,branch='regular')
    for bad in (replace_frame(frame),):
        try:owner.source.describe(bad)
        except ValueError:rejected+=1
    for args in ((3,0),(64,-1),(True,0)):
        try:current.correlated_qy_log_cap(c,*args)
        except ValueError:rejected+=1
    a=frame.roots['q'].atlas
    for value in (-161,161):
        try:current.mixed_exponential_change(a,current.MixedJet.constant(a,a.scalar(value)),N=160)
        except ValueError:rejected+=1
    assert rejected==6
    return dict(passed=True,issued_original_predicate_frames=frames,native_actual_five_C0_y_Z_yZ_density_rows=native,
        unchanged_accepted_original_phase_A_B_C0_Z_exact_range_comparisons=identical,
        copied_owner_source_partition_and_unproved_exponential_guards=rejected,
        qy_q_power_zero_gamma_q_minus_one_with_original_positive_q=True,
        original_nu_y_pressure_family_and_transverse_source_rows_retained=True,
        native_slow_density_derivatives_not_sharp_averaging_or_terminal_norms=True)


def replace_frame(frame):
    from dataclasses import replace
    return replace(frame)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('sharp_O2_phase_averaging_installed','actual_all_route_incoming_histories_installed','functional_terminal_identity_solved',
        'current_whole_N_selected','all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed',
        *current.source.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.CorrelatedO2MixedDensity();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (('source_correlated_qy_and_paired_flat_product_calculus',source_inequality_calculus),
                ('independent_original_qy_positive_eta_diagnostics',lambda:independent_source_qy_fixtures(owner)),
                ('independent_complete_nonlinear_mixed_density',lambda:independent_mixed_density_fixtures(owner)),
                ('actual_issued_source_phase_mixed_density_contract',lambda:actual_source_phase_density_contract(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,**checks,
        global_qy_source_magnitude_upper=manifest['global_qy_source_magnitude_upper'],
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw)),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Original-source correlated finite q_y bound, paired flat endpoint proof, actual mixed phase propagation, unchanged C0/Z ranges and complete nonlinear five slow-y/yZ density rows. Not sharp phase averaging, terminal/global N, all-route inherited history or recursion.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('Original source-correlated q_y and complete mixed five density PASS',flush=True);return report


if __name__=='__main__':run()
