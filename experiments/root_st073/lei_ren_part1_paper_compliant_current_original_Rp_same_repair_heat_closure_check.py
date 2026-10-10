"""Changed-boundary checks for the existing repair and closed raw heat API.

The accepted upstream and Gamma/contraction/pressure proofs are reused.
Independent cumulative primal equations check signs and physical units.
"""
import ast
import copy
import json
from pathlib import Path
import time
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_same_repair_heat_closure as current
import lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery as generic
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def same(a,b,label):assert s.simplify(s.expand(a-b))==0,label


def independent_closed_cumulative_proof():
    t,v=s.symbols('t v',real=True);delta=s.Symbol('delta',positive=True)
    K=s.Function('actual_Gamma_collar_bracket')
    Rt,Ev0,Pstar,B=s.symbols('Rtail Ev0 Pstar theta_base',positive=True)
    bh=(1+delta)/2;k=1-delta/2;prate=1+delta
    # The angular defect integral converges; its constant 1/k is retained.
    # Energy and pressure use the complete future, including infinity.
    A=1/k+s.exp(-k*t)*s.Integral(s.exp(k*v)*(1-K(v)),(v,t,s.oo))
    E=s.Integral(s.exp(-delta*v)*K(v)**2,(v,t,s.oo))
    Q=s.Integral(s.exp(-prate*v)*K(v)**2/2,(v,t,s.oo))
    R=Rt*s.exp(t);theta=B*K(t)*s.exp(-bh*t);X=A/K(t)
    energy=E*s.exp(delta*t)/(2*K(t)**2)
    P0=s.Function('independent_analytic_P0')(s.Symbol('Z',real=True))*Pstar**2
    totalP=-Ev0**2*B**2*Q;Mp=totalP-P0
    primal=dict(Mz=s.Integer(0),Mtheta=s.sqrt(2)*R**s.Rational(3,2)*Ev0*theta*X,
        Mtheta_z=s.Integer(0),Mztheta=R*Ev0**2*theta**2*energy,Mp=Mp)
    expected=dict(Mz=0,Mtheta=s.sqrt(2)*R**s.Rational(3,2)*Ev0*theta,
        Mtheta_z=0,Mztheta=-R*Ev0**2*theta**2/2,Mp=Ev0**2*theta**2/2)
    for name,value in primal.items():same(s.diff(value,t),expected[name],name+' closed cumulative derivative')
    same(s.diff(A,t)+k*A,K(t),'closed angular ODE with retained 1/k')
    same(s.diff(totalP,t),expected['Mp'],'absolute pressure derivative')
    same(s.diff(P0,t),0,'independent P0 derivative')
    theta0=s.Symbol('theta',real=True)
    assert generic.history_densities(theta0,0)==dict(m=0,h=theta0,k=0,e=-theta0**2/2,p=theta0**2/2)
    return dict(passed=True,independent_complete_future_primal_integrals=True,
        actual_angular_defect_integral_has_retained_one_over_k=True,
        five_signed_cumulative_equations_and_pressure_units_proved=True,
        independent_P0_retained_without_axis_reset=True,
        heat_pressure_scale_is_Ev0_squared_theta_base_squared=True)


def identical(a,b,label):
    assert a.order==b.order==5,label
    for j in range(6):assert a[j]._mpi_==b[j]._mpi_,(label,j)
    return 6


def zero_rows(value,label):
    for j in range(6):assert current.raw_source.radius.post.selected.inlet.endpoints(value[j])==(0,0),(label,j)
    return 6


@source_precision
def run(before=None):
    began=time.monotonic();raw=json.loads((current.HERE/current.NAME).read_bytes())
    for name,digest in raw['input_hashes'].items():assert current.sha(name)==digest,name
    owner=current.CurrentOriginalRpSameRepairHeatClosure(before,require_checked=False)
    assert not owner.acceptance_loaded and not any(raw[k] for k in current.GATES+current.OPEN)
    assert raw['source_family']==owner.family_record
    assert raw['actual_same_selected_source_graph']==owner.assert_graph()
    assert raw['actual_exact_input_view_graph']==owner.exact.graph
    assert raw['actual_absolute_radius_and_amplitude_source_binding']==current.raw_source.packed(owner.absolute_binding)
    source=ast.parse(Path(current.__file__).read_text(encoding='utf8'))
    calls={ast.unparse(node.func).split('.')[-1] for node in ast.walk(source) if isinstance(node,ast.Call)}
    assert not calls & {'replay_heat','replay_repair','replay_future','CurrentExactRepairBranch','CurrentLimitHeatPressureBridge'}
    upstream_nodes=copy.deepcopy(owner.input_raw.graph.nodes)
    views={};rows=0;zeros=0;stresszeros=0
    for name,chart,coordinate in current.VIEWS:
        view=owner.evaluate(chart,raw['fresh_Z'],coordinate);views[name]=view
        assert current.raw_source.packed(current.view_report(view))==raw['actual_five_closed_heat_views'][name],name
        packet=view['source_packet'];h=view['closed_normalized_histories'];theta=h['theta']
        expected=dict(Mz=theta*0,Mtheta=theta*h['angular']*owner.ctx.sqrt(2),Mtheta_z=theta*0,
            Mztheta=theta**2*h['energy'],Mp=h['closed_pressure']-h['original_P0'],
            P0=h['original_P0'],pressure=h['closed_pressure'],Utheta=theta)
        for key,value in view['raw_histories'].items():
            rows+=identical(value.in_exact_units(current.raw_source.POWERS[key]),expected[key],(name,key))
        assert view['raw_histories']['P0'].coefficients is packet['original_analytic_P0_Taylor_retained']
        assert packet['original_axis_pressure_not_replaced_or_tail_patched']
        assert packet['pressure_is_original_forward_function_after_terminal_identity']
        assert packet['whole_Z_pressure_function_identity_consumed']
        for key in ('actual_Dtheta_Taylor','actual_Cp_Taylor'):zeros+=zero_rows(packet[key],(name,key))
        for key in ('Mz','Mtheta_z'):zeros+=zero_rows(view['raw_histories'][key].coefficients,(name,key))
        derivative=dict(Mz=theta*0,Mtheta=theta*owner.ctx.sqrt(2),Mtheta_z=theta*0,
            Mztheta=-theta**2/2,Mp=theta**2/2,P0=theta*0,pressure=theta**2/2)
        for key,value in derivative.items():
            rows+=identical(view['first_native_radial_derivatives'][key].coefficients,value,(name,'radial',key))
        assert view['first_native_radial_derivatives']['Mp'].powers==tuple(map(current.raw_source.Fraction,(0,2,0)))
        assert view['raw_histories']['pressure'].powers==tuple(map(current.raw_source.Fraction,(0,0,2)))
        if chart=='heat_exterior':
            for mixed in packet['actual_stress_factored_mixed4'].values():
                for interval in mixed.values():
                    assert current.raw_source.radius.post.selected.inlet.endpoints(interval)==(0,0)
                    stresszeros+=1
    assert owner.input_raw.graph.nodes==upstream_nodes,'Closure changed the upstream expression graph'
    assert owner.exact.repair is owner.input_raw.post.before.seed.exact.repair
    assert owner.exact.future is owner.input_raw.post.before.future
    # A numerical cap or a one-unit shift cannot masquerade as the source
    # heat radius. Check the changed boundary without solving any repair.
    shifted=copy.copy(owner.before);shifted.radius=copy.copy(owner.before.radius)
    shifted.graph=copy.deepcopy(owner.before.graph);shifted.radius.graph=shifted.graph
    shifted.radius.maps={key:dict(value) for key,value in owner.before.radius.maps.items()}
    old=shifted.radius.maps['heat_exterior']['logR']
    shifted.radius.maps['heat_exterior']['logR']=shifted.graph.add(
        current.raw_source.radius.FunctionRef(shifted.graph,old.node),shifted.graph.one)
    try:current.absolute_source_binding(shifted,owner.exact)
    except ValueError:pass
    else:raise AssertionError('Shifted absolute heat source accepted')
    wrong=copy.copy(owner.exact);wrong.future=owner.exact.prior.future
    try:wrong.assert_graph()
    except ValueError:pass
    else:raise AssertionError('A stale future branch accepted')
    for chart,coordinate in (('heat_collar',4),('heat_exterior',2),('waiting',0)):
        try:owner.evaluate(chart,raw['fresh_Z'],coordinate)
        except ValueError:pass
        else:raise AssertionError('Non-heat or invalid closure chart accepted')
    proof=independent_closed_cumulative_proof()
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        actual_absolute_source_binding=owner.absolute_binding,
        independent_closed_cumulative_equations=proof,
        actual_five_closed_heat_call_coefficient_checks=rows,
        source_proved_constant_and_linear_zero_coefficient_checks=zeros,
        pure_Gamma_exterior_zero_stress_mixed_rows=stresszeros,
        same_existing_repair_future_and_P0_preserved=True,
        no_repair_constructor_or_replay_calls=True,upstream_expression_graph_unmodified=True,
        shifted_absolute_heat_radius_stale_future_and_invalid_domains_rejected=True,
        unchanged_angular_pressure_Gamma_and_contraction_math_not_replayed=True,
        scope='Same-existing-repair heat source closure and factorized raw histories; unrestricted numeric points, global mixed/Cartesian field and full corrected NS remain open',
        **dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name),
            Path(generic.__file__).name:current.sha(Path(generic.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.raw_source.packed(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_SAME_REPAIR_HEAT_CLOSURE existing branch, exact absolute source binding, closed raw heat histories',flush=True)
    return result


if __name__=='__main__':run()
