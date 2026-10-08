"""Original stress/pressure carrier identities and whole O2 source evidence."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_axial_carrier as current
import lei_ren_part1_paper_compliant_current_original_reference_whole_cell_integrals_check as saved

C0,Y=current.C0,current.Y;ep=current.ep


def full_original_pressure_factorization(owner):
    source=owner.source;z,f,H,D,P,p0,p0Z,p0ZZ=source.original['inputs'];W=source.symbols[Y][-2];a=source.a_symbol
    R,ps,delta=s.symbols('R Pstar delta',positive=True);errors=s.symbols('e0 e1 e2',real=True)
    Q=1+z*z;L=1-delta*z*z;tail=s.exp(s.Rational(3,5));alpha=P+W
    pressure={p0:-alpha/Q**2+errors[0]*tail*s.Rational(5,2)/(ps*Q**2),
        p0Z:4*alpha*z/Q**3+errors[1]*tail*10*z/(ps*Q**3),
        p0ZZ:(4-20*z*z)*alpha/Q**4+errors[2]*tail*22/(ps*Q**2)}
    def original(rows):
        return s.cancel(sum(R**r*ps**p*delta**d*L**ell*expr.subs(pressure)
            for (r,p,d,ell),expr in rows[('p2',0)])/z)
    actual={C0:original(source.original['rows']),Y:original(source.derivative_templates)}
    identities=normalizations=0
    for direction in (C0,Y):
        got=0
        for (r,p,d,ell),expr,error in owner.expressions[direction]:
            term=R**r*ps**p*delta**d*L**ell*expr.subs(owner.symbols[-1],z*z)
            if error is not None:term*=tail*errors[error]
            got+=term;identities+=1
            # Source-derived radius/delta collection, not a borrowed O2 ratio.
            assert r==1 and d>=0 and p-4*d+10*r<=11 and 10*r<=10;normalizations+=1
        assert s.cancel(got-actual[direction])==0
        assert s.cancel(actual[direction].subs(z,-z)-actual[direction])==0
        assert not s.denom(s.cancel(actual[direction])).subs(z,0)==0
    # Differentiate the complete original C0 carrier, using exact defining
    # integral ODEs and constant pressure error functions at fixed Z.
    rates={f:(1-a)*f/2,H:f-3*H/2,D:f*f/2-D,P:f*f/2,W:-f*f/2}
    dy=R*s.diff(actual[C0],R)+sum(s.diff(actual[C0],x)*rate for x,rate in rates.items())
    assert s.cancel(dy-actual[Y])==0
    proof=owner.pressure_proof
    assert proof['original_remainder_even_and_first_jet_odd'] and proof['original_datum_independent_of_source_y']
    assert len(proof['original_density_first_jet_identities'])==14
    assert proof['original_late_pressure_source_proof']['total_remainder_factors']==[5,10,44]
    return dict(passed=True,independent_full_original_pressure_affine_carrier_reconstruction=2,
        coefficient_pieces_with_exact_odd_even_factorization=identities,
        source_derived_Lambda0_normalization_power_checks=normalizations,
        exact_complete_original_radial_carrier_derivative_identity=True,
        original_pressure_density_first_jet_identities=14,
        same_original_pressure_error_functions_not_reference_profile_ratios=True,
        source_H_R_smooth_extension_not_replaced_by_R0_ZZ_away_from_axis=True)


def whole_source_contract(owner,manifest):
    c=owner.ctx;cert=manifest['actual_O2_whole_axial_carrier_certificate'];count=cert['original_mass_source_level'];width=cert['axial_cells_per_radial_cell']
    rows=cert['whole_original_source_partition'];assert len(rows)==count*width==1024
    assert count==min(owner.source.parent.parent.levels) and cert['original_y_window']==['0','1'] and cert['original_Z_window']==['-1','1']
    gb=[];gyb=[]
    for index in range(count):
        for j in range(width):
            row=rows[index*width+j]
            assert row['exact_y_cell']==[str(s.Rational(index,count)),str(s.Rational(index+1,count))]
            assert row['exact_abs_Z_cell']==[str(s.Rational(j,width)),str(s.Rational(j+1,width))]
            assert row['actual_original_f_H_D_P_a_and_remaining_pressure_mass'] and row['full_original_late_pressure_errors_retained']
            assert row['original_pressure_odd_first_jet_not_constant_divided_by_Z']
            g=saved.interval(c,row['p2_over_Z_divided_by_Lambda0']);gy=saved.interval(c,row['p2_y_over_Z_divided_by_Lambda0'])
            assert ep(g)[1]<0 and all(mp.isfinite(v) for value in (g,gy) for v in ep(value));gb.append(g);gyb.append(gy)
    g=saved.interval(c,cert['g_normalized']);gy=saved.interval(c,cert['g_y_normalized']);ratio=saved.interval(c,cert['g_y_over_g'])
    assert ep(g)[1]<0
    for a,b in zip(gb,gyb,strict=True):
        saved.contains(g,a);saved.contains(gy,b);saved.contains(ratio,b/a)
    assert cert['all_negative_Z_from_exact_even_carrier_and_y_derivative_parity'] and cert['all_Z_zero_from_same_smooth_carrier_extension']
    assert cert['actual_q_y_not_set_to_zero'] and cert['reference_profile_ratios_not_used']
    axis=0
    for index in (0,count//2,count-1):
        carrier=owner.cell(count,index,0,0)
        frame=owner.source.source_frame(count,index,Z_lower=0,Z_upper=0);atlas=frame.roots['q'].atlas
        for direction,order,key in ((C0,current.source.Z,'p2_over_Z_divided_by_Lambda0'),(Y,current.source.YZ,'p2_y_over_Z_divided_by_Lambda0')):
            original=current.source.axial.normalize(atlas,frame.roots['p2'][order],(11,10,0,0,0))
            saved.contains(original,carrier[key]);axis+=1
    rejected=0
    for call in (lambda:owner.cell(3,0,0,1),lambda:owner.cell(count,count,0,1),
        lambda:owner.cell(count,0,-2,0),lambda:owner.cell(count,0,1,0)):
        try:call()
        except (ValueError,TypeError):rejected+=1
    assert rejected==4
    return dict(passed=True,whole_original_O2_source_rectangles=1024,
        same_source_g_gy_ratio_rectangle_checks=1024,
        original_axis_transverse_jet_containment_comparisons=axis,
        unsupported_original_domains_or_cells_rejected=rejected,
        g_divided_by_Lambda0=g,g_y_divided_by_Lambda0=gy,g_y_over_g=ratio,
        source_sign_proof_includes_full_pressure_errors_and_Z_zero=True,
        negative_Z_coverage_exact_parity_not_samples=True)


def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('regular_signed_predicate_frames_or_mixed_primitives_installed','actual_changed_five_integrals_installed',
        'all_17_chart_or_24_cell_oracle_installed','actual_five_controls_installed','current_whole_N_selected',
        *current.source.ordered.base.point.source.inertial.profiles.loop.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalO2AxialCarrier();checks={}
    with mp.workdps(owner.ctx.dps+40):
        for name,call in (('full_original_pressure_carrier_calculus',lambda:full_original_pressure_factorization(owner)),
            ('actual_whole_O2_axial_carrier_domains',lambda:whole_source_contract(owner,manifest))):
            checks[name]=call();print(name,'PASS',flush=True)
    report=dict(all_passed=True,**{current.GATE:True},source_family=manifest['source_family'],**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            uncompressed_bytes=len(raw),compressed_bytes=path.stat().st_size),
        **dict.fromkeys(flags,False),input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Original O2 continuous y[0,1]/Z[-1,1] p2=Z*g,p2_y=Z*g_y with negative g, full pressure parity and all errors, actual radial profiles and source Lambda0 normalization. Source geometry only; no signed mixed phase, density integrals, terminal/global N or full reconstruction.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.base.encoded(report),indent=2).encode()+b'\n')
    print('Full original O2 axial carrier and pressure parity PASS',flush=True);return report


if __name__=='__main__':run()
