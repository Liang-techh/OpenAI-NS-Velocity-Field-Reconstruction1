"""Independent finite-function bounds and current source/domain attachment checks."""
import json
from pathlib import Path
from types import SimpleNamespace
import math
import sympy as s
import lei_ren_part1_paper_compliant_current_generic_shear_source_bounds as source


def independent_finite_quotient_fixture(c):
    y,Z=s.symbols('ordinary_logR Z');n=3+2*y+y*Z-Z**2;d=2+y*y+Z*Z/4
    actual=n/d;domains=(c.mpf('.2'),c.mpf('.3'))
    def bound(expr):
        poly=s.Poly(expr,y,Z)
        result=sum((c.mpf(str(abs(coeff)))*domains[0]**power[0]*domains[1]**power[1]
            for power,coeff in poly.terms()),c.mpf(0))
        return source.LogUpper.constant(c,result)
    numerator={key:bound(s.diff(n,y,key[0],Z,key[1])) for key in source.ORDERS}
    denominator={key:bound(s.diff(d,y,key[0],Z,key[1])) for key in source.ORDERS}
    table=source.quotient_table(c,numerator,denominator,c.ln(2))
    count=0
    for yy in (s.Rational(-1,5),s.Integer(0),s.Rational(1,5)):
        for zz in (s.Rational(-3,10),s.Integer(0),s.Rational(3,10)):
            for key,cover in table.items():
                value=abs(s.diff(actual,y,key[0],Z,key[1]).subs({y:yy,Z:zz}))
                if value:
                    exact=c.mpf(int(s.numer(value)))/int(s.denom(value))
                    if cover.log is None or source.packets.recovery.endpoints(c.ln(exact))[1]>source.packets.recovery.endpoints(cover.log)[0]:
                        raise ArithmeticError('Independent signed rational-function quotient bound failed')
                count+=1
    return count


def independent_modal_fixture(c):
    logs=(c.ln(2),2*c.ln(3),2*c.ln(5),2*c.ln(7))
    algebra=source.packets.FactoredAlgebra(c,logs,[])
    z=source.packets.IntervalTaylor.variable(c,c.mpf([-1,1]),5)
    row=algebra.shift(algebra.lift(2+z),(-1,.5,0,0))+algebra.shift(algebra.lift(2-z),(0,0,1,0))
    packet=SimpleNamespace(algebra=algebra,provenance={'logR_cover':c.mpf(['-.5','.5'])})
    count=0
    for k in (0,1,2):
        cover=source.modal_partial_bound(packet,row,k,radius_power=1)
        for Z in (-1,0,1):
            for y in ('-.5','0','.5'):
                actual=((c.mpf(3)/2*(2+Z)+25*(2-Z)) if k==0 else c.mpf(3)/2-25 if k==1 else c.mpf(0))*c.exp(c.mpf(y))
                if source.packets.recovery.endpoints(actual)!=(0,0):
                    if cover.log is None or source.packets.recovery.endpoints(c.ln(abs(actual)))[1]>source.packets.recovery.endpoints(cover.log)[0]:
                        raise ArithmeticError('Independent source-mode/radius normalization bound failed')
                elif k==2 and cover.log is not None:raise ArithmeticError('Exact-zero derivative lost')
                count+=1
    # A formal factor with an enormous log stays a finite log bound.
    extreme=source.packets.FactoredAlgebra(c,(c.mpf('1e800'),c.mpf('-1e800'),c.mpf(0),c.mpf(0)),[])
    row=extreme.shift(extreme.lift(z),(1,-1,0,0))
    packet=SimpleNamespace(algebra=extreme,provenance={'logR_cover':c.mpf('1e801')})
    cover=source.modal_partial_bound(packet,row,1,radius_power=1)
    if cover.log is None or not all(source.mp.isfinite(v) for v in source.packets.recovery.endpoints(cover.log)):
        raise ArithmeticError('Formal extreme source factor was resolved or overflowed')
    return count+1


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Current source norm dependency changed: '+name)
    owner=source.CurrentGenericSourceBounds();c=owner.ctx
    if source.packets.encode(owner.theorem)!=data['exact_quotient_and_radius_derivative_theorem']:
        raise ValueError('Current mixed quotient/radius theorem changed')
    if source.packets.encode(owner.lower)!=data['positive_noncore_source_denominator_theorem']:
        raise ValueError('Current positive original source lower changed')
    rebuilt={chart:owner.chart(chart) for chart in source.packets.CHARTS}
    if source.packets.encode(rebuilt)!=data['current_original_source_log_bound_charts']:
        raise ValueError('Current complete signed log norm chart changed')
    if len(data['admitted_quotient_source_charts'])!=15 or data['positive_noncore_source_denominator_theorem']['unresolved_source_charts']!=['core']:
        raise ArithmeticError('All15 noncore source quotients required, core remains separate')
    rows=0
    for chart,record in rebuilt.items():
        admitted=record['admitted_original_quotient_log_norms']
        if chart=='core':
            if admitted or record['quotient_derivatives_bound_on_checked_whole_chart']:
                raise ArithmeticError('Positive compact core cache admitted an axis or loop denominator')
            continue
        if set(admitted)!=set(('a','b','p1','p2','t0')):
            raise ArithmeticError('Full signed source quotient omitted')
        for table in admitted.values():
            if set(table)!=set('y%d_Z%d'%key for key in source.ORDERS):
                raise ArithmeticError('Required ordinary quotient derivative missing')
            rows+=len(table)
        if not all(record[key] for key in ('numerator_all_signed_inertial_pressure_energy_meridional_terms_retained',
            'original_radial_prefactor_inertial_shift_applied_once','original_width_inverse_and_source_factor_correlation_retained')):
            raise ArithmeticError('Current source term/derivative coordinate omitted')
    invalid=0
    for method,args in ((owner.chart,('unknown',)),
        (source.quotient_table,(c,{},{},c.mpf('inf'))),
        (source.LogUpper,(c,c.mpf('nan')))):
        try:method(*args)
        except ValueError:invalid+=1
        else:raise ArithmeticError('Invalid source/log norm input admitted')
    if any(data[key] or any(row[key] for row in rebuilt.values()) for key in source.OPEN):
        raise ArithmeticError('Original log norms admitted a modified/global/recursive field')
    for key in ('whole_upstream_source_derivative_norms_certified','whole_generic_scales_instantiated',
        'phase_held_loop_primitive_derivative_bounds_certified'):
        if data[key]:raise ArithmeticError('Incomplete loop domain/scales/primitive norms promoted')
    rational=independent_finite_quotient_fixture(c);modal=independent_modal_fixture(c)
    result=dict(all_passed=True,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),source_family=owner.family,
        exact_mixed_quotient_and_inertial_radius_identities=len(owner.theorem['exact_mixed_derivative_identities']),
        original_inner_positive_source_AST_bindings=len(owner.inner['original_actual_core_bridge_switch_AST_bindings']),
        original_O2_positive_source_AST_bindings=len(owner.outer['original_actual_O2_profile_AST_bindings']),
        raw_original_whole_cover_charts_checked=len(rebuilt),analytic_positive_denominator_noncore_charts_checked=15,
        full_signed_quotient_derivative_log_bounds_checked=rows,
        independent_exact_rational_function_derivative_envelopes=rational,
        independent_moderate_modal_radius_and_extreme_log_fixtures=modal,invalid_source_and_log_inputs_rejected=invalid,
        original_core_axis_or_loop_denominator_not_promoted=True,full_O2_buffer_offset_0_11_retained=True,
        whole_upstream_source_derivative_norms_certified=False,whole_generic_scales_instantiated=False,
        phase_held_loop_primitive_derivative_bounds_certified=False,source_ancestor_constructors_called=False,
        scope=data['scope'],input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),
            Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Current original15-chart quotient derivative log bounds and independent normalized fixtures PASS',flush=True)
    return result


if __name__=='__main__':run()
