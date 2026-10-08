"""Source-scoped O2 calculus, continuous domains and factored q evidence."""
from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_mixed_source_jets as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

C0,Y,Z,YZ=current.C0,current.Y,current.Z,current.YZ
ep=current.ep


def defining_profile_calculus():
    y,v=s.symbols('y v',real=True);J=s.Function('J');sigma=s.Function('sigma')
    f=s.exp(y/10-s.Rational(3,5)*J(y));a=s.Rational(4,5)+s.Rational(6,5)*sigma(y)
    H=s.exp(-3*y/2)*(s.Rational(5,8)+s.Integral(s.exp(8*v/5-3*J(v)/5),(v,0,y)))
    D=s.exp(-y)*(s.Rational(5,12)+s.Integral(s.exp(6*v/5-6*J(v)/5),(v,0,y))/2)
    P=s.Rational(5,2)+s.Integral(s.exp(v/5-6*J(v)/5),(v,0,y))/2
    rates=((1-a)*f/2,f-3*H/2,f*f/2-D,f*f/2)
    for value,expected in zip((f,H,D,P),rates,strict=True):
        assert s.simplify((s.diff(value,y)-expected).subs(s.diff(J(y),y),sigma(y)))==0
    return dict(passed=True,independent_original_defining_integral_derivatives=4,
        fundamental_theorem_and_original_J_y_sigma_used=True,
        original_incoming_constants_retained=['5/8','5/12','5/2'])


def full_original_row_calculus(owner):
    z,f,H,D,P,p0,p0Z,p0ZZ=owner.original['inputs'];y=s.Symbol('y',real=True)
    symbols=(f,H,D,P);functions=tuple(s.Function('original_'+str(x))(y) for x in symbols)
    rates=((1-owner.a_symbol)*f/2,f-3*H/2,f*f/2-D,f*f/2)
    back=dict(zip(functions,symbols));back.update(dict(zip((s.diff(x,y) for x in functions),rates)))
    radial=axial=affine=correlated=0
    for key,terms in owner.original['rows'].items():
        derived=owner.derivative_templates[key]
        for (powers,expr),(dpowers,dexpr) in zip(terms,derived,strict=True):
            assert powers==dpowers
            independent=s.diff(s.exp(powers[0]*y)*expr.xreplace(dict(zip(symbols,functions))),y)/s.exp(powers[0]*y)
            assert s.cancel(independent.xreplace(back)-dexpr)==0;radial+=1
            for u in (p0,p0Z,p0ZZ):
                for v in (p0,p0Z,p0ZZ):assert s.diff(dexpr,u,v)==0
            affine+=1
    # Independently check that the original Z rows differentiate L(Z) exactly
    # once, including full pressure first/second jets. Radial profiles stay fixed.
    R,ps,delta=s.symbols('R Pstar delta',positive=True);pressure=s.Function('P0')(z)
    substitution={p0:pressure,p0Z:s.diff(pressure,z),p0ZZ:s.diff(pressure,z,2)}
    def complete(name,order):
        return sum(R**r*ps**p*delta**d*(1-delta*z*z)**ell*expr.subs(substitution)
            for (r,p,d,ell),expr in owner.original['rows'][(name,order)])
    for name in {key[0] for key in owner.original['rows']}:
        assert s.cancel(s.diff(complete(name,0),z)-complete(name,1))==0;axial+=1
    # Differentiate the exact alpha=P+W function relation, W_y=-P_y.
    W=owner.symbols[Y][-2]
    for key,terms in owner.rows[C0].items():
        dterms=owner.rows[Y][key]
        for (powers,expr,sens),(dpowers,dexpr,dsens) in zip(terms,dterms,strict=True):
            independent=powers[0]*expr+sum(s.diff(expr,u)*v for u,v in zip(symbols,rates,strict=True))-s.diff(expr,W)*f*f/2
            assert powers==dpowers and s.cancel(independent-dexpr)==0;correlated+=1
    return dict(passed=True,full_original_radial_coefficient_identities=radial,
        independent_original_Z_derivative_rows=axial,full_pressure_affine_error_rows=affine,
        exact_prefix_suffix_derivative_identities=correlated,
        original_radius_y_and_single_L_Z_derivative_retained=True)


def varying_parameter_calculus(owner):
    a,eta,ay=s.symbols('a eta a_y',positive=True)
    q=s.sqrt((1+eta)/a-s.Rational(1,2));qy=-(1+eta)*ay/(2*a*a*q)
    assert s.simplify(s.diff(q,a)*ay-qy)==0
    nu=1+2*q*q;assert s.simplify(nu-2*(1+eta)/a)==0
    assert s.simplify(s.diff(nu,a)*ay+2*(1+eta)*ay/a**2)==0
    assert s.simplify(4*q*qy+2*(1+eta)*ay/a**2)==0
    p=mp.mp.clone();p.dps=100;eta0=p.exp(-12);comparisons=0
    def sigma(y):
        left=p.exp(-1/y**2);right=p.exp(-1/(1-y)**2)
        return left/(left+right)
    def av(y):return p.mpf(4)/5+p.mpf(6)/5*sigma(y)
    def qv(y):return p.sqrt((1+eta0)/av(y)-p.mpf('.5'))
    for text in ('.23','.53','.97'):
        y=p.mpf(text);a0=av(y);ay0=p.diff(av,y)
        expected=-(1+eta0)*ay0/(2*a0*a0*qv(y))
        assert abs(p.diff(qv,y)-expected)<p.mpf('1e-85');comparisons+=1
        expectednu=-2*(1+eta0)*ay0/(a0*a0)
        assert abs(p.diff(lambda yy:1+2*qv(yy)**2,y)-expectednu)<p.mpf('1e-85');comparisons+=1
    return dict(passed=True,exact_q_and_phase_normalization_nu_identities=4,
        independent_finite_original_flat_profile_derivative_diagnostics=comparisons,
        finite_eta_diagnostic_not_selected_original_eta=True,
        nu_is_phase_normalization_not_physical_viscosity=True)


def continuous_source_contract(owner,manifest):
    c=owner.ctx;count=manifest['ordered_source_cells'];records=manifest['actual_whole_O2_mixed_source_cells']
    assert count==min(owner.parent.parent.levels) and len(records)==count
    terms=errors=roots_checked=0;last=None
    for index,row in enumerate(records):
        assert row['source_family']==owner.family and row['original_P0_datum_sha256']==owner.family['datum_enclosure_sha256']
        assert row['exact_y_cell']==[str(s.Rational(index,count)),str(s.Rational(index+1,count))]
        assert row['exact_Z_range']==['-1','1'] and row['source_index']==index and row['source_level']==count
        assert row['ordinary_Z_yZ_not_scaled_coordinate_derivatives']
        assert row['separate_source_hulls_not_exact_selected_compatible_field']
        assert not row['mixed_phase_inverse_or_primitives_installed'] and not row['actual_five_controls_installed']
        roots=row['actual_root_jets'];assert set(roots)=={'E','V','p1','p2','a','b','t0','q','nu'}
        for name,jets in roots.items():
            assert set(jets)=={str(x) for x in (C0,Y,Z,YZ)}
            for jet in jets.values():
                assert jet['encloses_original_source_function'] and not jet['point_value_selected']
                assert jet['formal_positive_scale']['radius_power']==0;roots_checked+=1
        q=roots['q'][str(C0)];assert q['sign']=='positive' and not q['exact_zero']
        assert q['formal_positive_scale']['source_exponents']==[0,0,0,1]
        assert ep(saved.interval(c,q['coefficient_interval']))==(1,1)
        logq=saved.interval(c,row['native_positive_logq']);assert all(mp.isfinite(x) for x in ep(logq)) and ep(logq)[1]<0
        assert roots['q'][str(Y)]['formal_positive_scale']['source_exponents']==[0,0,0,-1]
        assert ep(saved.interval(c,roots['a'][str(Y)]['coefficient_interval']))[0]>=0
        assert ep(saved.interval(c,roots['q'][str(Y)]['coefficient_interval']))[1]<=0
        assert ep(saved.interval(c,roots['nu'][str(Y)]['coefficient_interval']))[1]<=0
        for name in ('a','q','nu'):
            assert all(roots[name][str(key)]['exact_zero'] for key in (Z,YZ))
            assert not roots[name][str(Y)]['exact_zero']
        for name in ('b','t0'):
            assert all(v['exact_zero'] for v in roots[name].values())
        assert roots['V'][str(Y)]['exact_zero'] and roots['V'][str(YZ)]['exact_zero']
        assert not roots['V'][str(Z)]['exact_zero']
        for item in row['actual_original_mixed_source_terms']:
            assert item['ordinary_y_order'] in (0,1) and item['ordinary_Z_order'] in (0,1)
            for term in item['terms']:
                terms+=1
                for error in term['full_pressure_error_terms']:
                    assert error['full_original_error_source']['formal_positive_scale']['radius_power']==0;errors+=1
        assert row['pressure_prefix_suffix_function_correlation']['passed'];last=row
    assert terms==38*count and errors>0
    endpoint=current.positive.endpoint_logq(c,s.Rational(1),current.ordered.interval(c,owner.owner.scales.logs['eta']))
    assert ep(saved.interval(c,last['native_positive_logq']))[0]==ep(endpoint)[0]
    extras=manifest['actual_axis_and_signed_source_cells'];assert len(extras)==5
    for row in extras[:3]:
        assert row['exact_Z_range']==['0','0']
        jets=row['actual_root_jets']['p2']
        assert jets[str(C0)]['exact_zero'] and jets[str(Y)]['exact_zero']
        assert not jets[str(Z)]['exact_zero'] and not jets[str(YZ)]['exact_zero']
    # Only new interface guards are exercised; accepted ancestors are not rerun.
    frame=owner.source_frame(count,0,Z_lower=0,Z_upper=0)
    rejected=0
    for call in (lambda:owner.describe(replace(frame)),lambda:owner.source_frame(3,0,Z_lower=0,Z_upper=0),
        lambda:owner.source_frame(count,count,Z_lower=0,Z_upper=0),
        lambda:owner.source_frame(count,0,Z_lower=-2,Z_upper=0),
        lambda:owner.source_frame(count,0,Z_lower=1,Z_upper=-1)):
        try:call()
        except (ValueError,TypeError):rejected+=1
    assert rejected==5
    return dict(passed=True,continuous_whole_original_source_cells=count,
        source_coefficient_terms_retaining_original_factors=terms,late_pressure_error_terms=errors,
        original_mixed_root_rows_checked=roots_checked,axis_and_signed_source_frames=5,
        unsupported_domain_or_unissued_frame_rejections=rejected,
        positive_original_endpoint_q_preserved_formally=True,
        near_endpoint_q_y_enclosure_may_be_very_wide=True,
        no_native_huge_scale_materialization_or_constant_q_substitution=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('mixed_phase_inverse_or_primitives_installed','actual_changed_five_integrals_installed',
        'all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed','current_whole_N_selected',
        *current.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2MixedSources();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (('original_defining_profile_calculus',defining_profile_calculus),
            ('original_full_coefficient_and_pressure_calculus',lambda:full_original_row_calculus(owner)),
            ('actual_varying_q_nu_calculus',lambda:varying_parameter_calculus(owner)),
            ('whole_O2_mixed_source_contract',lambda:continuous_source_contract(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            uncompressed_bytes=len(raw),compressed_bytes=path.stat().st_size),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],
            current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Original O2 slope y[0,1]/Z[-1,1] source C0/y/Z/yZ plus varying a,q,phase-normalization nu. Full original pressure and positive endpoint q; source interface only. Mixed inverse, primitives, changed five integrals, all-route controls and full reconstruction remain open.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('Original whole O2 mixed sources and actual q/nu derivatives PASS',flush=True)
    return report


if __name__=='__main__':run()
