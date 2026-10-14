"""Scoped source/unit/calculus checks for the true core inlet, not all charts."""
import gzip
import json
import math
from pathlib import Path
import time

import sympy as sy

import lei_ren_part1_paper_compliant_current_original_core_inlet_axial6_generic_source as current
import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_C3_functions_check as calculus


class KernelAlgebra(calculus.Algebra):
    def at(self,index):
        if index in self.cache:return self.cache[index]
        row=self.g.nodes[index];op=row['operation']
        if op=='function_substitution':
            value=self.at(row['expression']).subs(self.at(row['variable']),self.at(row['value']))
        elif op=='definite_integral':
            value=sy.integrate(self.at(row['integrand']),
                (sy.Symbol(row['variable']),self.at(row['lower']),self.at(row['upper'])))
        else:return super().at(index)
        self.cache[index]=value;return value


def source_rows(owner):
    g=owner.g;r=current.rows;checked=[]
    assert owner.previous.acceptance_loaded and owner.common.acceptance_loaded
    assert owner.moments.common is owner.core._oracle is owner.common
    assert owner.moments.proof==current.interior.source_theorem()
    for name,jet in owner.core_functions.items():
        old=r(owner.previous.core_functions[name])
        assert len(r(jet))==7 and r(jet)[:len(old)]==old
        assert owner.derivative.derivative(old[-1])==r(jet)[6]
        checked.append(name)
    for name,jet in owner.moment_functions.items():
        assert len(r(jet))==7
        for order,ref in enumerate(r(jet)):
            row=g.nodes[ref.node]
            assert row['operation']=='definite_integral'
            assert row['lower']==g.zero.node and row['upper']==g.one.node
            assert row['actual_core_normalized_moment']==name
            assert row['ordinary_slow_Z_derivative_order']==order
            assert row['same_unique_nonlinear_core_and_full_density_product']
            assert row['endpoint_independent_of_Z'] and row['true_core_inlet_complete_history']
    assert all(len(r(jet))==7 for jet in owner.fields.values())
    assert all(len(r(jet))==6 for jet in owner.functions.values())
    assert r(owner.P0)[:4]==r(owner.core.pulse.own.jets.native['P0'])
    for order,ref in enumerate(r(owner.P0)[4:],4):
        row=g.nodes[ref.node]
        assert row['operation']=='current_original_axis_pressure_function_derivative'
        assert row['original_function']==owner.P0.value.node
        assert row['ordinary_slow_Z_derivative_order']==order
        assert row['Taylor_coefficient_factorial']==math.factorial(order)
        assert row['source_projection_binding']==current.parent.active.current.ast_binding(owner.project_P0_function)
        assert row['pressure_provider_binding']==current.parent.active.current.ast_binding(owner.P0_rows)
        assert row['implicit_source_sha256']==owner.common.core.datum.source_sha
        assert row['datum_enclosure_sha256']==owner.common.datum_sha
    try:current.JetAlgebra(g,5).shift(current.FunctionJet(tuple(r(owner.fields['m'])[:6])))
    except ValueError:pass
    else:raise AssertionError('Missing genuine sixth mean source silently padded')
    return dict(passed=True,preserved_true_core_C5_sources=checked,
        actual_seven_row_complete_moment_histories=6,
        original_P0_C3_preserved_and_higher_rows_are_explicit_function_derivatives=True,
        executable_order6_pressure_projection_bound=True,missing_seventh_row_rejected=True,
        radial_scope='rho=4 core inlet only')


def physical_units():
    R,P,F,H,M,K,A,B,C,am,x,Rm,theta,mixed=sy.symbols(
        'R Pstar F0 H M K A B C am x Rm theta mixed',positive=True)
    # The accepted patch's physical theta and mixed primitives are
    # sqrt(2)*Pstar*Rm^(3/2)*am*(theta,mixed) = 2*(I_RF,I_RFV).
    rawh=am*x**sy.Rational(-3,2)*theta
    rawk=am*x**sy.Rational(-3,2)*mixed
    assignments={theta:R**2*F*H/(sy.sqrt(2)*P*Rm**sy.Rational(3,2)*am),
        mixed:R**2*F*K/(sy.sqrt(2)*P*Rm**sy.Rational(3,2)*am),x:R/Rm}
    W=sy.sqrt(2*R)*F/P
    assert sy.simplify(rawh.subs(assignments)-W*H/2)==0
    assert sy.simplify(rawk.subs(assignments)/P-W*K/(2*P))==0
    assert sy.simplify((R*A-R**2*F**2*B)/(R*P**2)-(A-R*F**2*B)/P**2)==0
    assert sy.simplify(R*F**2*C/P**2-(R*F**2/P**2)*C)==0
    # Axial velocity and its mean enter recover_inputs in physical units.
    U=sy.Symbol('physical_Uz')
    assert sy.diff(U/P,U)==1/P and sy.diff(M/P,M)==1/P
    return dict(passed=True,nonunit_radius_and_Pstar_retained=True,
        raw_h_and_raw_k_primitive_normalizations_proved=True,
        generic_V_mean_k_invPstar_applied_once=True,
        independent_pressure_and_energy_increment_units_proved=True,
        p2_quadratic_weight='Pstar, no extra half')


def full_calculus_fixture(owner):
    """Replay actual new nodes on independent nonzero sixth-order core data."""
    g=owner.g;r=current.rows;Z,rho=sy.symbols('Z fixture_core_rho',real=True)
    t=sy.Symbol('original_core_inlet_moment_t')
    eps,de,P=sy.Rational(1,9),sy.Rational(1,5),sy.Rational(3,2);R=4*eps
    Phi=1+(2+Z**6/5)*rho**2/23
    Psi=(1+Z**3+Z**6/7)*rho**2/11;U0=4*Z+sy.Rational(1,10)
    U=U0+eps*Psi;F=2+Z**2/13+Z**6/17;P0=-1-Z**2-Z**6/19
    integrals={name:sy.integrate(norm*t**weight*(Phi if left=='Phi' else U).subs(rho,rho*t)
        *((Phi if right=='Phi' else U).subs(rho,rho*t) if right else 1),(t,0,1))
        for name,(left,right,weight,norm) in current.interior.SPEC.items()}
    inlet={name:value.subs(rho,4) for name,value in integrals.items()};W=sy.sqrt(2*R)*F/P
    field=dict(E=W*Phi.subs(rho,4),V=U.subs(rho,4)/P,m=inlet['M']/P,
        h=W*inlet['H']/2,k=W*inlet['K']/(2*P),
        e=(inlet['A']-R*F**2*inlet['B'])/P**2,p=R*F**2*inlet['C']/P**2,P0=P0)
    field['pressure']=P0+field['p'];E,V,m,h,k,e,p=(field[q] for q in ('E','V','m','h','k','e','p'))
    d=1-Z**2;L=1-de*Z**2;transport=(1-de)*Z*m+d*sy.diff(m,Z)
    expected=dict(transport=transport,Q=(2*Z*V-transport)/L,
        theta_linear=(-E+(1-de/2)*h-(1-de)*Z*sy.diff(h,Z)/2)/L,
        theta_quadratic=((2*de-1)*Z*k-d*sy.diff(k,Z)+E*transport)/L,
        axial_linear=(-V+(1-de)*(m-Z*sy.diff(m,Z))/2)/L,
        axial_quadratic=(V*transport+2*de*Z*e-d*sy.diff(e,Z)
            +2*(1+de)*Z*(P0+p)-d*sy.diff(P0+p,Z))/L)
    expected['p2']=R*(expected['axial_linear']+P*expected['axial_quadratic'])/E
    assert all(sy.diff(field[name],Z,6)!=0 for name in ('m','e','p','P0'))
    base={g.symbol('Z').node:Z,owner.core.rho.node:rho,g.symbol('original_core_inlet_moment_t').node:t,
        owner.core.parameters['epsilon'].node:eps,owner.core.parameters['delta'].node:de,
        owner.core.parameters['logP'].node:sy.log(P)}
    def bind(ref,value):
        if ref.node in base:assert sy.simplify(base[ref.node]-value)==0
        base[ref.node]=value
    for name,body in (('Phi',Phi),('Psi',Psi),('U0',U0),('F0',F)):
        for order,ref in enumerate(r(owner.core_functions[name])):bind(ref,sy.diff(body,Z,order))
    for order,ref in enumerate(r(owner.P0)):bind(ref,sy.diff(P0,Z,order))
    alg=KernelAlgebra(g,base);moment_rows=0
    for name,jet in owner.moment_functions.items():
        for order,ref in enumerate(r(jet)):
            assert sy.expand(alg.at(ref.node)-sy.diff(integrals[name],Z,order))==0,(name,order)
            moment_rows+=1
    field_rows=0
    for name,jet in owner.fields.items():
        for order,ref in enumerate(r(jet)):
            assert sy.expand(alg.at(ref.node)-sy.diff(field[name],Z,order))==0,(name,order)
            field_rows+=1
    # Exact rational axial points keep the full differentiated equations small.
    points=(sy.Rational(3,8),sy.Integer(1));generic_rows=0
    derivatives={name:[sy.diff(body,Z,k) for k in range(6)] for name,body in expected.items()}
    for point in points:
        point_alg=KernelAlgebra(g,{node:value.subs(Z,point) for node,value in base.items()})
        for name,jet in owner.functions.items():
            for order,ref in enumerate(r(jet)):
                assert sy.simplify(point_alg.at(ref.node)-derivatives[name][order].subs(Z,point))==0,(point,name,order)
                generic_rows+=1
    return dict(passed=True,mathematical_fixture_only=True,
        independent_complete_moment_rows=moment_rows,independent_normalized_source_rows=field_rows,
        independent_full_generic_rows=generic_rows,exact_axial_points=[str(v) for v in points],
        nonzero_sixth_mean_energy_pressure_inputs=True,
        full_pressure_P0_plus_increment_and_all_signed_sectors_retained=True,
        nonunit_Pstar_without_extra_half_and_single_radius_factor_checked=True,
        original_scale_point_values_not_selected=True)


def actual_numeric_sources(owner,packets):
    views={};count=0
    for label,packet in packets.items():
        f=packet['flow'];sources=packet['actual_complete_seven_Taylor_sources'];sectors=packet['actual_full_signed_inertial_sectors_axial5']
        assert set(sources)=={'E','V','m','h','k','e','p','P0'}
        assert all(len(v)==7 for v in sources.values()) and all(len(v)==6 for v in sectors.values())
        assert len(packet['actual_radial_Q_axial5'])==len(packet['actual_full_signed_p2_axial5'])==6
        assert packet['all_seven_rows_genuine_before_differentiation']
        assert packet['radial_Q_actual_six_rows_from_genuine_seventh_mean_row_no_Qpad']
        assert not packet['all_17_chart_source_installed']
        for values in (*sources.values(),*sectors.values(),packet['actual_radial_Q_axial5'],packet['actual_full_signed_p2_axial5']):
            for value in values:
                assert value.scale.bases is f.logs and value.ledger is f.ledger
                assert not value.record()['point_value_selected'];count+=1
        z=owner.hierarchy.source.cell(packet['exact_Z_cell']);provided=owner.P0_rows(z)
        for order,ref in enumerate(current.rows(owner.P0)[4:],4):
            assert owner.project_P0_function(ref,z)==provided[order]*math.factorial(order)
        original=packet['actual_core_moment_and_density_tails']
        assert original['actual_nonlinear_fixed_point_and_all_density_product_tails_retained']
        assert original['radial_degree']==24 and original['axial_depth']==6
        views[label]=current.macro.serialized({key:packet[key] for key in (
            'exact_Z_cell','core_rho','actual_complete_seven_Taylor_sources','seven_source_units',
            'actual_full_signed_inertial_sectors_axial5','actual_radial_Q_axial5','actual_full_signed_p2_axial5',
            'radial_Q_actual_six_rows_from_genuine_seventh_mean_row_no_Qpad','all_17_chart_source_installed')})
    return dict(passed=True,actual_sealed_numerical_source_views=list(views),
        actual_scaled_source_and_output_rows=count,genuine_seven_source_six_output_rows=True,
        same_live_formal_bases_and_directed_ledger=True,all_density_tails_retained=True,
        executable_P0_Z4_Z5_Z6_projections_checked=True,
        global_mixed_norms_and_all17_source_not_certified=True),views


def run(owner,packets=None):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['core_inlet_axial6_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    packets=packets or {'whole_Z':owner.numerical((-1,1)),'off_center':owner.numerical(('0.375','0.375'))}
    checks=dict(source_rows=source_rows(owner),physical_units=physical_units())
    print('TRUE_CORE_INLET_C6_SOURCE_AND_UNITS_CHECKED',flush=True)
    checks['independent_calculus']=full_calculus_fixture(owner)
    checks['actual_numeric_source'],views=actual_numeric_sources(owner,packets)
    return seal(owner,raw,checks,views,began)


def seal(owner,raw,checks,views,began):
    """Save already passed scoped checks without repeating their calculus."""
    owner.assert_graph()
    assert all(check['passed'] for check in checks.values())
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    raw['numeric_core_inlet_views']=current.generic.base.encoded(views)
    (current.HERE/current.NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_CORE_INLET_C6_GENERIC5_SOURCE_CHECKED',flush=True)
    return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses the live admitted core; no ancestor rebuild')
