"""Compare live five-basis recovery with unchanged generic recovery math."""
import gzip
import json
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_Rm_generic_inputs as current
import lei_ren_part1_paper_compliant_current_generic_shear_source_packets as packets
import lei_ren_part1_paper_compliant_current_generic_shear_inputs as inputs
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow

fields,ep=current.fields,current.ep


def fixture(zvalue):
    """Independent existing generic packet on finite diagnostic source bases."""
    c=MPIntervalContext();c.dps=160
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    z=packets.IntervalTaylor.variable(c,c.mpf(zvalue),5);de=c.mpf('.2')
    zero=packets.IntervalTaylor.constant(c,0,5)
    source=dict(m=1+z+z*z,h=2-z+z**3*c.mpf('.2'),k=3+z*c.mpf('.4')-z*z*c.mpf('.2'),
        e=z*z*c.mpf('1.4')-c.mpf('.5'),p=-z*c.mpf('.2')+z**3+c.mpf('.3'))
    theta=z*c.mpf('.2')+z**3+c.mpf('.8');theta_y=z*z*c.mpf('.07')-c.mpf('.1')
    axial=-z*c.mpf('.4')+z*z*c.mpf('.1')+c.mpf('1.2');axial_y=z*z*c.mpf('.1')+c.mpf('.02')
    P0=z*c.mpf('.3')-z*z*c.mpf('.2')+c.mpf('.9')
    raw=dict(histories={key:[f.jet(row)]+[f.jet(zero)]*4 for key,row in source.items()},
        velocity=dict(theta=[f.jet(theta),f.jet(theta_y)]+[f.jet(zero)]*3,
            axial=[f.jet(axial),f.jet(axial_y)]+[f.jet(zero)]*3,
            radial=[f.jet(zero)[:5]]*5))
    op=SimpleNamespace(flow=f,c=c,zrows=f.jet(z),reference=SimpleNamespace(z=z,delta=de),
        P0=f.jet(P0),Pstar=f.factor((0,.5,0,0,0)),Rm_factor=f.scalar(7))
    packet=dict(raw_current_radius_y_derivative_axial_coefficients=raw,
        original_P0_normalized_axial5=op.P0,geometry=dict(point=True,exact_x=[5,4]))
    new=current.recover_inputs(op,packet)
    algebra=packets.FactoredAlgebra(c,(c.mpf(0),c.ln(9),c.mpf(0),c.mpf(0)),[])
    rows=lambda row:algebra.lift(row)
    native_m={key:(rows(row),)*5 for key,row in source.items()}
    native_v=dict(theta=(rows(theta),rows(theta_y),rows(zero),rows(zero),rows(zero)),
        axial=(rows(axial),rows(axial_y),rows(zero),rows(zero),rows(zero)),radial=(rows(zero),)*5)
    common_m={key:tuple(algebra.shift(row,packets.INVERSE_S) if key in ('m','k') else row for row in values)
              for key,values in native_m.items()}
    common_v={key:tuple(algebra.shift(row,packets.INVERSE_S) if key in ('axial','radial') else row for row in values)
              for key,values in native_v.items()}
    old=packets.CurrentSourcePacket('actual_patch',{},algebra,z,common_v,common_m,
        (rows(P0+source['p']),)+(rows(zero),)*4,rows(P0),native_v,native_m,{})
    expected=inputs.from_packet(old,de);recovered=old.recover_original(de);count=0
    def compare(actual,reference):
        nonlocal count
        reference=algebra.resolve(reference)
        assert len(actual)==reference.order+1
        for a,b in zip(actual,reference.coefficients):
            assert ep(f.ordinary_cover(a)-b)[0]<=0<=ep(f.ordinary_cover(a)-b)[1]
            count+=1
    for key in ('E','C','B','inertial_theta','inertial_axial','positive_denominator',
                'kappa','kappa_minus2','H0_minus2','D','J'):
        polynomial=getattr(expected,{'inertial_theta':'inertial_theta','inertial_axial':'inertial_axial',
            'positive_denominator':'denominator','kappa':'kappa_numerator','kappa_minus2':'kappa_excess_numerator',
            'H0_minus2':'stronger_numerator','D':'direction_numerator','J':'transverse_numerator'}.get(key,key))
        value=sum((row*(c.mpf(35)/4)**power for power,row in polynomial.terms.items()),algebra.lift(0))
        compare(new['actual_generic_source_numerators'][key],value)
    for name,key in (('theta_linear','inertial_theta_linear'),('theta_quadratic','inertial_theta_quadratic'),
                     ('axial_linear','inertial_axial_linear'),('axial_quadratic','inertial_axial_quadratic')):
        compare(new['full_signed_inertial_sectors_axial4'][name],recovered[key])
    compare(new['common_radial_Q_axial4'],recovered['Ur_over_S_sqrt_R_over_2'])
    compare(new['common_absolute_pressure_axial5'],old.absolute_pressure[0])
    assert new['common_original_P0_axial5'] is op.P0
    assert not algebra.final_rows and not algebra.proofs
    return dict(passed=True,fixture_Z=zvalue,unchanged_generic_recovery_and_all_signed_numerator_comparisons=count,
        finite_diagnostic_scales_only=True,native_source_factors_not_resolved=True)


def native(owner,report):
    counts=0;P0=0;scope=0
    for label in ('0','.5'):
        op=owner.upstream.owner(label).op
        for name,call in (('active_point',lambda:owner.evaluate(label,(5,4))),
                          ('terminal_cell',lambda:owner.cell(label,(71,40),'Rh'))):
            packet=call();record=report['frames'][label][name]
            assert current.base.encoded(fields.serialized(packet))==record
            assert packet['common_original_P0_axial5'] is op.P0;P0+=1
            assert packet['no_source_basis_projection_or_directed_offset_materialization']
            assert packet['full_inertial_linear_quadratic_pressure_meridional_sectors_retained']
            assert not packet['actual_positive_E_C_and_generic_cone_admission_certified']
            assert not packet['actual_phase_or_finite_N_density_integrals_installed'];scope+=1
            for key,rows in packet['actual_generic_source_numerators'].items():
                assert len(rows)==(6 if key in ('E','C','B','positive_denominator','kappa','kappa_minus2') else 5)
                assert all(v.scale.bases is op.flow.logs and v.ledger is op.flow.ledger for v in rows)
                counts+=len(rows)
            assert len(packet['common_radial_Q_axial4'])==5
    return dict(passed=True,actual_live_full_signed_generic_numerator_coefficients=counts,
        exact_separate_P0_object_checks=P0,native_source_scope_checks=scope,
        original_five_source_bases_and_directed_offsets_retained=True)


def run():
    began=time.monotonic();fixtures=[]
    with mp.workdps(180):
        for z in ('.3','-.4'):
            fixtures.append(fixture(z));print('Independent unchanged generic recovery fixture PASS',z,flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalRmGenericInputs(require_checked=False);evidence=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    assert report['original_full_recovery_source_bindings']==owner.bindings
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name,Path(packets.__file__).name,Path(inputs.__file__).name):
        fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        unchanged_generic_source_fixtures=fixtures,actual_live_source=evidence,
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual Rm full signed generic source inputs PASS',flush=True);return result


if __name__=='__main__':run()
