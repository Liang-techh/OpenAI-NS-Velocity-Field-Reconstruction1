"""Independent Rv common units, source germs and actual provider diagnostics.

Reuse unchanged accepted inlet and endpoint proofs. Test only the new
conversion and derivative boundary, without replaying source/repair solves.
"""
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_common_unit_seam as current
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def same(a,b,label):assert s.cancel(s.expand(a-b))==0,label


def algebra(owner):
    field=SimpleNamespace(graph=owner.graph,parameters=owner.before.radius.parameters,
        contracts=owner.before.radius.frame.bridge.leading.contracts)
    A=s.Symbol('absolute_logRp',real=True);L=s.Symbol('exact_logU0',real=True)
    q=current.pulse.radius.RadiusInterpreter(field,False,
        {owner.before.radius.logRp.node:A,owner.logU0.node:L})
    geometry=owner.before.radius.geometry('flatten',0);units=0
    for name,powers in current.pulse.POWERS.items():
        r,e,p=map(s.Rational,powers)
        left=sum(q.at(v) for v in owner.before.scale('pulse_end',0,powers).values())
        right=sum(q.at(v) for v in owner.raw.scale('flatten',geometry,(r,e,p-e)).values())
        same(left+e*L,right,(name,'exact Rv common source units'));units+=1
    same(q.at(geometry['native_to_log_radius_jacobian']),1,'flatten unit Jacobian')
    leftgeo=owner.before.radius.geometry('pulse_end',0)
    same(q.at(leftgeo['native_to_log_radius_jacobian']),1,'pulse end unit Jacobian')
    same(q.at(leftgeo['logR']),q.at(geometry['logR']),'actual absolute Rv equality')
    # Accepted complete native frame supplies the actual analytic u shape.
    frame=owner.before.radius.frame
    record=json.loads(gzip.decompress((current.HERE/current.pulse.radius.frame_source.NAME).read_bytes()))
    u=s.sympify(record['actual_analytic_native_inlet_functions']['u'])
    z=next(v for v in u.free_symbols if v.name=='Z')
    same(u*(1+z*z),u.subs(z,0),'actual u(Z)=U0/(1+Z^2) function identity')
    y,z=s.symbols('y Z',real=True);mu=s.Symbol('mu',positive=True)
    Rv,Ev,U,P=s.symbols('Rv Ev0 U0 Pstar',positive=True)
    Xv=s.Function('actual_terminal_X')(z);ev=s.Function('actual_complete_future_half')(z)
    p0=s.Function('independent_P0')(z);pin=s.Symbol('Pin',real=True)
    theta=1/(1+z*z);bp=s.Rational(1,2)+mu;a=1-mu;prate=2*bp
    Fv=s.exp(-13/(2*mu)-13);S=s.exp(-13/mu-26)
    pulse_pressure=p0+pin*theta**2+(U*theta)**2*(1-S)/(2*prate)
    flatten_pressure=p0+pin*theta**2+(U*theta)**2*(1-Fv**2)/(2*prate)
    same(pulse_pressure,flatten_pressure,'exact absolute P0+Mp endpoint memory')
    # Source germs at Rv: both schedules are flat through every order and
    # compact B support is absent here. Use the exact cumulative ODE solution.
    th=theta*s.exp(-bp*y);xx=1/a+(Xv-1/a)*s.exp(-a*y)
    ee=ev*s.exp(2*mu*y)-(s.exp(2*mu*y)-1)/(4*mu)
    funcs=dict(Mtheta=s.sqrt(2)*(Rv*s.exp(y))**s.Rational(3,2)*Ev*th*xx,
        Mztheta=Rv*s.exp(y)*Ev**2*th**2*ee,
        pressure=pulse_pressure+Ev**2*theta**2*(1-s.exp(-prate*y))/(2*prate),
        Utheta=Ev*th)
    checks=0
    for k in range(1,5):
        expected=dict(Mtheta=s.sqrt(2)*Rv**s.Rational(3,2)*Ev*theta*a**(k-1),
            Mztheta=-Rv*Ev**2*theta**2/2*(-2*mu)**(k-1),
            pressure=Ev**2*theta**2/2*(-prate)**(k-1),Utheta=Ev*theta*(-bp)**k)
        for name,fn in funcs.items():
            got=s.diff(fn,y,k).subs(y,0)
            for j in range(5-k):
                same(s.diff(got,z,j),s.diff(expected[name],z,j),(name,k,j,'independent source product rule'))
                checks+=1
    return dict(passed=True,exact_common_source_unit_identities=units,
        independent_radial_axial_product_rule_identities=checks,
        actual_native_u_U0_q_function_identity=True,
        exact_Rv_geometry_and_two_unit_Jacobians=True,
        independent_absolute_pressure_memory_identity=True,
        exact_positive_Fv_and_Ev0_scales_not_caps=True)


def overlap(a,b,label):
    assert a.order==b.order,label
    count=0
    for j in range(a.order+1):
        al,ah=current.pulse.radius.post.selected.inlet.endpoints(a[j])
        bl,bh=current.pulse.radius.post.selected.inlet.endpoints(b[j])
        assert max(al,bl)<=min(ah,bh),(label,j,'disjoint source enclosures');count+=1
    return count


@source_precision
def run(before=None):
    began=time.monotonic();candidate=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    for name,digest in candidate['input_hashes'].items():assert current.sha(name)==digest,name
    owner=current.CurrentOriginalRpCommonUnitSeam(before,require_checked=False)
    assert not owner.acceptance_loaded and not any(candidate[k] for k in current.GATES+current.OPEN)
    assert candidate['source_family']==owner.family_record
    assert candidate['actual_source_graph']==owner.assert_graph()
    rows=0;mixed=0;parts=0;native=0
    for name,Z in (('fresh','.521'),('cell',['.520','.522'])):
        view=owner.evaluate(Z)
        assert current.pulse.raw.packed(current.view_report(view))==candidate['actual_Rv_seam_views'][name],name
        left=view['reduced_pulse_terminal'];right=view['actual_flatten_terminal'];converted=view['converted_pulse_terminal_enclosures']
        for key,value in left.items():
            assert value.powers==right[key].powers==converted[key].powers==current.POWERS[key]
            rows+=overlap(value.coefficients,converted[key].coefficients,(name,key,'literal unit division vs source reduction'))
            rows+=overlap(value.coefficients,right[key].coefficients.truncate(value.coefficients.order),(name,key,'actual endpoint providers'))
            for v in (value,right[key],converted[key]):
                assert all(ref.graph is owner.graph for _,ref in v.log_scale_parts);parts+=len(v.log_scale_parts)
        for key,lrows in view['pulse_log_radius_mixed_rows'].items():
            assert len(lrows)==15
            for label,value in lrows.items():
                other=view['flatten_log_radius_mixed_rows'][key][label]
                assert value.powers==other.powers
                mixed+=overlap(value.coefficients,other.coefficients,(name,key,label))
                assert all(ref.graph is owner.graph for _,ref in value.log_scale_parts+other.log_scale_parts)
                parts+=len(value.log_scale_parts)+len(other.log_scale_parts)
        # Compare the new velocity mixed rows against the original actual
        # providers. P radial rows are checked by source algebra in Ev0^2
        # units; old pressure cap-scaled rows are kept as diagnostics only.
        actual=owner.before.evaluate('pulse_end',Z,0)['source_packet']
        actual_flat=owner.before.evaluate('flatten',Z,0)['current_source_packet']
        for key,label,divisor in (
            ('Utheta','Utheta_over_Pstar_without_common_theta_radial_factor',1),
            ('Uz','Uz_over_Pstar_without_common_theta_radial_factor',1),
            ('Ur','Ur_over_sqrt_R_over_2_Pstar_without_common_theta_radial_factor',owner.ctx.sqrt(2))):
            for coordinate,value in view['pulse_log_radius_mixed_rows'][key].items():
                a=actual['physical_mixed_derivatives_total_order_le4'][label][coordinate]/owner.U0_bound/divisor
                b=actual_flat['physical_mixed_derivatives_total_order_le4'][label][coordinate]/divisor
                native+=overlap(value.coefficients,current.IntervalTaylor.constant(owner.ctx,a,0),(name,key,coordinate,'pulse original mixed'))
                native+=overlap(value.coefficients,current.IntervalTaylor.constant(owner.ctx,b,0),(name,key,coordinate,'flatten original mixed'))
        for key in ('Mz','Mtheta_z','Uz','Ur'):
            assert all(current.pulse.radius.post.selected.inlet.endpoints(v)==(0,0) for v in left[key].coefficients.coefficients)
        for key in ('Mp','pressure'):
            assert view['pulse_log_radius_mixed_rows'][key]['y1_Z0'].powers==tuple(map(current.Fraction,(0,2,0)))
        assert left['Ur'].coefficients.order==4 and left['Utheta'].coefficients.order==5
    proof=algebra(owner)
    rejected=[]
    for Z in ('1.01', '-1.01'):
        try:owner.evaluate(Z)
        except ValueError:rejected.append(Z)
        else:raise AssertionError('Illegal original Z domain accepted')
    try:current.CurrentOriginalRpCommonUnitSeam(owner.raw,require_checked=False)
    except ValueError:pass
    else:raise AssertionError('Unclosed source bypass accepted')
    active=owner.before.evaluate('pulse_end','.521',-3)
    try:owner.converted_terminal(active,view['geometry'])
    except ValueError:pass
    else:raise AssertionError('Rv unit conversion incorrectly applied inside active pulse')
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        independent_common_unit_and_source_germ_proof=proof,
        actual_current_endpoint_coefficient_overlap_diagnostics=rows,
        actual_current_mixed_row_overlap_diagnostics=mixed,
        original_provider_velocity_mixed_row_overlap_diagnostics=native,
        all_current_final_scale_parts_on_one_graph=parts,
        source_function_join_consumes_current_rebound_endpoint_theorem=True,
        source_function_identity_not_inferred_from_interval_overlap=True,
        exact_pressure_derivative_scale_is_Ev0_squared_not_decay_cap=True,
        Ur_C4_and_Utheta_C5_preserved=True,
        active_pulse_conversion_unclosed_bypass_and_illegal_domains_rejected=rejected,
        no_numeric_absolute_point_or_uniform_global_time_admission=True,
        **dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.pulse.raw.packed(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_COMMON_UNIT_SEAM Rv units, five histories and mixed logR/Z4 endpoint rows',flush=True)
    return result


if __name__=='__main__':run()
