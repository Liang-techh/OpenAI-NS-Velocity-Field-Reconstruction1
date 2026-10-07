"""Independent differentiation checks for the explicit scalar jet backend."""
import json
from pathlib import Path
import time
import sympy as s
import lei_ren_part1_paper_compliant_current_generic_loop_point_Z as point

loop=point.loop;HERE,NAME=point.HERE,point.RECEIPT


def symbolic():
    r,psi=s.symbols('r psi',real=True);den=1-2*r*s.cos(psi)+r*r
    w=(s.cos(psi)-r)/den
    assert s.simplify(s.diff(w,r)-(-den+2*(s.cos(psi)-r)**2)/den**2)==0
    u=s.symbols('u',real=True)
    assert s.simplify(s.diff(u/s.sqrt(1+u*u),u)-(1+u*u)**s.Rational(-3,2))==0
    a,b,az,bz=s.symbols('a b az bz',real=True)
    assert s.simplify(s.diff(-b/a,a)*az+s.diff(-b/a,b)*bz-(-bz/a+b*az/a**2))==0
    return dict(passed=True,independent_exact_chain_rule_identities=3)


def references():
    scales=loop.GenericLoopScales(a_min='1',margin_min='2',boundary_kappa_excess_min='.02',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='1',dps=60)
    c=scales.ctx;eta=scales.eta;step=c.mpf('1e-10');tolerance=c.mpf('1e-21')
    fixtures=[('positive-small-r','1.5','.2','.05'),('negative-r','1.5','-.2','-.5'),
        ('zero-r','1.5','.2','0'),('large-r','1.3','-.2','.9'),
        ('q-seam-Delta0',c.mpf(2),'0','.2'),
        ('q-transition',2+eta/2,'0','.2'),('q-flat',2+eta,'0','.2')]
    jets=dict(a_Z=c.mpf('.03'),b_Z=c.mpf('-.04'),p2_Z=c.mpf('.02'),E_Z=c.mpf('.1'))
    values=dict(p1=c.mpf(6),E=c.mpf(2))
    def at(a,b,p2,z):
        return loop.GenericShearLoop(scales,a=c.mpf(a)+z*jets['a_Z'],b=c.mpf(b)+z*jets['b_Z'],
            p2=c.mpf(p2)+z*jets['p2_Z'],p1=values['p1'],Utheta=values['E']+z*jets['E_Z'])
    def derivative(fn):
        return (-fn(2*step)+8*fn(step)-8*fn(-step)+fn(-2*step))/(12*step)
    errors=[];implicit=[];rows=[];direction_count=0;endpoint_count=0
    for label,a,b,p2 in fixtures:
        backend=point.GenericLoopPointZ(scales,a=a,b=b,p2=p2,**values,**jets)
        for phase in (c.mpf('.17'),c.mpf('.53')):
            got=backend.evaluate(phase)
            for key,jet in (('A',got.A_Z_slow),('B',got.B_Z_slow),('angle',got.angle_Z)):
                reference=derivative(lambda z:at(a,b,p2,z).evaluate(phase)[key])
                error=abs(reference-jet);assert error<tolerance,(label,phase,key,error)
                errors.append(error)
            source=at(a,b,p2,c.mpf(0));expected=source.evaluate(phase)
            assert abs(got.A-expected['A'])<c.mpf('1e-50')
            assert abs(got.B-expected['B'])<c.mpf('1e-50')
            assert got.approximate_only and not got.certified_original_point_oracle
            psi=got.angle
            total_phase_Z=derivative(lambda z:at(a,b,p2,z).phase_at_angle(psi))+expected['phase_angle_derivative']*got.angle_Z
            assert abs(total_phase_Z)<tolerance,(label,total_phase_Z)
            implicit.append(abs(total_phase_Z))
            rows.append(dict(fixture=label,phase=phase,A=got.A,B=got.B,A_Z=got.A_Z_slow,B_Z=got.B_Z_slow,
                approximate_manufactured_point_reference=True))
        for psi in (c.mpf('.31'),c.pi,c.mpf('5.8')):
            reference=derivative(lambda z:at(a,b,p2,z).direction(psi))
            assert abs(reference-backend.direction_Z(psi))<tolerance,(label,'direction')
            direction_count+=1
        for phase in (0,1,2,-1):
            q=backend.evaluate(phase)
            assert q.A==q.B==q.A_Z_slow==q.B_Z_slow==q.angle_Z==0
            endpoint_count+=1
    failed=False
    try:point.GenericLoopPointZ(scales,a='1.5',b='.2',p2='.05',**values,**{**jets,'E_Z':'nan'})
    except ValueError:failed=True
    assert failed
    return dict(passed=True,fixture_count=len(fixtures),primitive_and_inverse_Z_comparisons=len(errors),
        direction_Z_comparisons=direction_count,fixed_phase_implicit_identity_comparisons=len(implicit),
        periodic_endpoint_zero_rows=endpoint_count,maximum_absolute_Z_reference_error=max(errors),
        maximum_fixed_phase_implicit_identity_error=max(implicit),
        nonconstant_a_b_p2_E_signed_r_q0_q_transition_and_small_r_covered=True,
        explicit_normalization_B_inherits_E_units=True,
        numeric_error_or_native_source_certification_claimed=False,examples=rows)


def run():
    began=time.monotonic();record=json.loads((HERE/point.NAME).read_bytes())
    assert record[point.GATE]
    for name,digest in record['input_hashes'].items():assert point.sha(name)==digest,name
    hashes={**record['input_hashes'],point.NAME:point.sha(point.NAME),Path(__file__).name:point.sha(Path(__file__).name)}
    result=dict(all_passed=True,**{point.GATE:True},source_family=record['source_family'],
        independent_symbolic=symbolic(),independent_point_references=references(),
        **{key:record[key] for key in ('quadrature_or_roundoff_certified','original_native_parameter_point_values_selected',
            'numerical_original_source_point_or_integral_oracle_installed','actual_five_controls_installed',
            'certified_actual_fixed_point_tail_installed','actual_terminal_Z_function_closure_installed','current_whole_N_selected',*loop.OPEN)},
        input_hashes=hashes,execution_seconds=time.monotonic()-began,
        scope='Explicit-input original scalar formula checks only; no native source point value selection, quadrature certification, global N or actual control/field admission.')
    assert not any(result[key] for key in loop.OPEN)
    (HERE/NAME).write_text(json.dumps(loop.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Scalar original loop phase-held A/B/Z and inverse/direction references PASS',flush=True)
    return result


if __name__=='__main__':run()
