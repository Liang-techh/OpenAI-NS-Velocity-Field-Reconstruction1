"""Actual finite-N Rc joint targets, original amplitude and independent jets."""
import gzip
import hashlib
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_Rc_full_predicate_joint_targets as current

ep=current.ep;iv=current.packets.interval


def same(record,value):
    scale=record['formal_positive_scale'];c=value.ctx
    assert tuple(scale['source_exponents'])+(scale['radius_power'],)==value.scale.powers
    assert ep(iv(c,scale['additional_log_interval']))==ep(value.scale.offset)
    assert ep(iv(c,record['coefficient_interval']))==ep(value.coefficient)
    assert record['exact_zero']==value.zero and not record['point_value_selected'] and record['encloses_original_source_function']


def overlaps(a,b):
    assert a.ctx is b.ctx and a.scale.bases is b.scale.bases and a.ledger is b.ledger
    anchor=a.ctx.mpf(max(ep(a.scale.evaluate())[1],ep(b.scale.evaluate())[1]))
    left=a.coefficient*a.bounded_exp(a.scale.evaluate()-anchor)
    right=b.coefficient*b.bounded_exp(b.scale.evaluate()-anchor)
    return max(ep(left)[0],ep(right)[0])<=min(ep(left)[1],ep(right)[1])


def independent_target_algebra():
    z=s.Symbol('Z',real=True);mu=s.Symbol('mu',positive=True)
    A=s.Function('A')(z);m,h,k,e,p=[s.Function(key)(z) for key in ('m','h','k','e','p')]
    C=k-A*m;CZ=s.diff(k,z)-s.diff(A,z)*m-A*s.diff(m,z)
    assert s.simplify(s.diff(C,z)-CZ)==0
    rows=(m/A,C/(mu*A**2),h/A,e/A**2,p/A**2)
    expected=(s.diff(m,z)/A-m*s.diff(A,z)/A**2,
        (CZ-2*C*s.diff(A,z)/A)/(mu*A**2),s.diff(h,z)/A-h*s.diff(A,z)/A**2,
        s.diff(e,z)/A**2-2*e*s.diff(A,z)/A**3,s.diff(p,z)/A**2-2*p*s.diff(A,z)/A**3)
    for value,derivative in zip(rows,expected,strict=True):assert s.simplify(s.diff(value,z)-derivative)==0
    return dict(passed=True,independent_whole_function_joint_Z_and_five_quotient_identities=6,
        mu_Z_exact_zero_for_original_constant_source_parameter=True)


def polynomial_jet_fixtures(owner):
    # These finite polynomials test differentiation/normalization only. They
    # are explicitly separate from original source functions and target data.
    z=s.Symbol('Z');A=2+z*z;mu=s.Rational(1,7)
    H=dict(m=1+3*z+z*z,h=2-z+z**3,k=4+2*z-z*z,e=3+z**4,p=5-2*z+2*z*z)
    R=dict(M=H['m']/A,I=H['h']/A,S=H['e']/(A*A),Cp=H['p']/(A*A))
    R[current.parameters.repair.ROWS[1]]=(H['k']-A*H['m'])/(mu*A*A)
    c=owner.ctx;coords=owner.coordinates;scalar=coords.scalar;count=0
    for point in (s.Rational(-2,3),s.Rational(0),s.Rational(1,3),s.Rational(1)):
        parse=lambda value:c.mpf(str(s.N(value.subs(z,point),c.dps+40)))
        history={key:scalar(parse(value)) for key,value in H.items()}
        jets={key:scalar(parse(s.diff(value,z))) for key,value in H.items()}
        av=parse(A);az=parse(s.diff(A,z));mv=parse(mu)
        got=current.fixed.fixed_N_target_rows(history,jets,scalar(av),scalar(az),c.ln(av),scalar(mv),c.ln(mv))
        for field,derivative in (('values',False),('Z_derivatives',True)):
            for key,value in R.items():
                exact=parse(s.diff(value,z) if derivative else value)
                enclosed=got[field][key].finite_interval()
                assert max(ep(enclosed)[0],ep(exact)[0])<=min(ep(enclosed)[1],ep(exact)[1]);count+=1
    return dict(passed=True,independent_finite_polynomial_quotient_C0_Z_comparisons=count,
        fixtures_not_original_source_or_control_data=True)


def actual_source_inheritance_and_targets(owner,saved):
    assert saved['actual_accepted_Rc_correction_binding']==current.encode(owner.binding)
    ancestor=owner.manifest['actual_full_predicate_O2_correction_to_Rc_continuation'];inherited=0
    for row,field in zip(owner.corrections,current.FIELDS,strict=True):
        for index,group in enumerate(('values','Z_derivatives')):
            for key,value in row[group].items():
                same(ancestor[field][index][key],value);inherited+=1
                assert value.ctx is owner.ctx and value.scale.bases is owner.coordinates.bases and value.ledger is owner.coordinates.ledger
    live=owner.targets()
    assert current.encode(live['report'])==saved
    amplitude=live['amplitude'];A,AZ,mu,logA,logmu=[amplitude[key] for key in ('A','AZ','mu_source','logA','logmu')]
    assert ep(amplitude['mu'])==ep(owner.target_owner.repair_mu)
    assert amplitude['record']['original_Rc_power_offset']==2
    comparisons=0;caps=0
    for index,(correction,target) in enumerate(zip(owner.corrections,live['targets'],strict=True)):
        H,J=correction['values'],correction['Z_derivatives'];record=saved['actual_five_joint_target_alternatives'][index]
        C=H['k']-A*H['m'];CZ=J['k']-AZ*H['m']-A*J['m']
        same(record['actual_joint_numerator'],C);same(record['actual_joint_numerator_Z'],CZ)
        same(record['actual_joint_numerator'],target['joint_numerator']);same(record['actual_joint_numerator_Z'],target['joint_numerator_Z'])
        den=A*A;denmu=den*mu
        expected=dict(M=H['m'].positive_divide(A,logA),
            I=H['h'].positive_divide(A,logA),S=H['e'].positive_divide(den,2*logA),Cp=H['p'].positive_divide(den,2*logA))
        expected[current.parameters.repair.ROWS[1]]=C.positive_divide(denmu,2*logA+logmu)
        expectedZ={out:(J[key]*A-degree*H[key]*AZ).positive_divide(den if degree==1 else den*A,(degree+1)*logA)
            for out,key,degree in (('M','m',1),('I','h',1),('S','e',2),('Cp','p',2))}
        expectedZ[current.parameters.repair.ROWS[1]]=(CZ*A-2*C*AZ).positive_divide(A*A*A*mu,3*logA+logmu)
        for key in current.parameters.repair.ROWS:
            same(record['actual_target_C0'][key],target['values'][key]);same(record['actual_target_Z'][key],target['Z_derivatives'][key])
            assert overlaps(expected[key],target['values'][key]);assert overlaps(expectedZ[key],target['Z_derivatives'][key]);comparisons+=2
            cap=current.parameters.repair.LogUpper.add(owner.ctx,[current.parameters.magnitude(target['values'][key]*current.N),
                current.parameters.magnitude(target['Z_derivatives'][key]*current.N)])
            assert current.encode(cap.record())==record['actual_N_scaled_target_C1_upper_caps'][key];caps+=1
        assert not H['p'].zero and not J['p'].zero
    assert not amplitude['record']['original_separate_P0']['exact_zero']
    assert saved['exact_original_24_cell_control_range_precondition_met'] is False
    assert saved['full24_original_C1_integral_range_transport_enclosed'] is False
    return dict(passed=True,lossless_actual_Rc_correction_C0_Z_inheritance_rows=inherited,
        fresh_actual_original_Rc_amplitude_C0_Z_source_queries=1,actual_joint_numerator_C0_Z_comparisons=4,
        independently_rearranged_five_target_C0_Z_overlaps=comparisons,actual_N_scaled_target_C1_caps=caps,
        actual_nonzero_pressure_correction_and_separate_P0_verified=True,
        six_cell_and_upstream_O2_numerical_suites_not_replayed=True,original_five_slope_cell_precondition_stays_open=True)


def guards(owner):
    count=0
    for callback in (lambda:owner.targets(N=True),lambda:owner.targets(N=7),lambda:owner.targets(N=2048),
        lambda:owner.normalise(dict(owner.corrections[0]),owner.amplitude),
        lambda:owner.normalise(owner.corrections[0],dict(owner.amplitude))):
        try:callback()
        except ValueError:count+=1
        else:raise AssertionError('Unissued source object or mismatched candidate admitted')
    original=owner.corrections[0]['values']['m']
    foreign=current.downstream.history.CommonSourceCoordinates(owner.ctx,owner.coordinates.logP_squared,owner.family)
    try:
        owner.corrections[0]['values']['m']=foreign.scalar(1)
        try:owner.normalise(owner.corrections[0],owner.amplitude)
        except ValueError:count+=1
        else:raise AssertionError('Foreign source ledger admitted')
    finally:owner.corrections[0]['values']['m']=original
    return dict(passed=True,candidate_issued_source_object_and_ledger_guards=count)


@current.native.inlet.source_precision
def run():
    began=time.monotonic();path=current.HERE/current.NAME;raw=gzip.decompress(path.read_bytes());manifest=json.loads(raw)
    assert manifest[current.GATE] and manifest['candidate_N']==1024
    flags=('full24_original_C1_integral_range_transport_enclosed','actual_five_controls_installed',
        'functional_terminal_identity_solved','current_whole_N_selected',*current.packets.OPEN)
    assert all(manifest[flag] is False for flag in flags)
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    bridge,_=current.native.inlet.native_bridge_owner();checks={}
    with current.native.inlet.CheckedSourceRuntime():
        owner=current.OriginalRcFullPredicateJointTargets(bridge);saved=manifest['actual_full_predicate_Rc_joint_targets']
        with mp.workdps(owner.ctx.dps+40):
            for label,callback in (('independent_joint_and_ordinary_Z_quotient_algebra',independent_target_algebra),
                ('independent_finite_polynomial_target_jets',lambda:polynomial_jet_fixtures(owner)),
                ('actual_Rc_source_inheritance_amplitude_and_joint_targets',lambda:actual_source_inheritance_and_targets(owner,saved)),
                ('candidate_issued_source_and_ledger_guards',lambda:guards(owner))):
                checks[label]=callback();print(label+' PASS',flush=True)
    receipt=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,candidate_N=1024,**checks,
        compressed_producer_report=dict(filename=current.NAME,lossless_original_json_sha256=hashlib.sha256(raw).hexdigest(),
            compressed_bytes=path.stat().st_size,uncompressed_bytes=len(raw)),
        actual_original_Rc_joint_target_C0_Z_covers_installed=True,**dict.fromkeys(flags,False),
        input_hashes={**manifest['input_hashes'],current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Actual source-owned full-predicate N1024 whole-Z correction errors yield five finite joint targets with fresh original Rc amplitude/ordinary-Z and separate P0. Independent quotient identities/jets and source inheritance. Missing five original slope-cell control adapter, quantitative cancellation, fixed point, terminal/global N/recursion/full NS remain open.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.encode(receipt),indent=2).encode()+b'\n')
    print('Actual original full-predicate Rc five joint C0/Z targets PASS',flush=True);return receipt


if __name__=='__main__':run()
