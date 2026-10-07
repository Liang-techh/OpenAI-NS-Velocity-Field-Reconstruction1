"""Whole analytic attachment and independent finite full-primitive checks."""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_reshape_relaxed_inputs as source


def independent_power_primitive_fixture():
    """Artificial moderate full axis primitives, not current source values.

    Exact positive kernel integrals are evaluated independently of the
    current raw-stress and recovery programs. Both signs of Z and inlet
    mean errors are retained. This checks units/formulas, not admission.
    """
    P=source.packets;c=MPIntervalContext();c.dps=115;S=c.mpf(3);R0=c.mpf(110)
    count=0
    for zp in ('-.7','-.1','.1','.7'):
        z=P.IntervalTaylor.variable(c,zp,5);u0=1/(1+z*z);V=4*z+c.mpf('.0003')
        algebra=P.FactoredAlgebra(c,(c.mpf(0),2*c.ln(S),c.mpf(0),c.mpf(0)),[])
        lift=algebra.lift
        for y in ('0','.7','2'):
            theta=c.exp(c.mpf(y));R=R0*theta;u=u0*theta**c.mpf('.1')
            mt=u0*(c.sqrt(2)*R0**c.mpf('1.5'))*(theta**c.mpf('1.6')-1)/c.mpf('1.6')
            mtz=V*mt;mz=(4*z-c.mpf('.0001'))*R0+V*(R-R0)
            mzt=z*z*(R0*16)+V*V*(R-R0)-u0*u0*R0*(theta**c.mpf('1.2')-1)/c.mpf('2.4')
            mp0=u0*u0*(theta**c.mpf('.2')-1)/c.mpf('.4');P0=(3+z)/(S*S)
            histories=dict(m=mz/(R*S),h=mt/(c.sqrt(2)*R**c.mpf('1.5')*S),
                k=mtz/(c.sqrt(2)*R**c.mpf('1.5')*S*S),e=mzt/(R*S*S),p=mp0/(S*S))
            erows=[lift(u/S*c.mpf('.1')**j) for j in range(5)]
            vrows=[lift(V/S)]+[lift(0)]*4
            family=dict(zip(P.FAMILY_KEYS,('finite_fixture','finite_source','finite_datum')))
            packet=P.CurrentSourcePacket('finite_fixture',family,algebra,z,
                dict(theta=erows,axial=vrows,radial=[lift(0)]*5),
                {k:[lift(v)]+[lift(0)]*4 for k,v in histories.items()},
                [lift(P0+histories['p'])]+[lift(0)]*4,lift(P0),{}, {},
                dict(original_radius_source='artificial finite radius',cache_cover=True))
            actual=source.inputs.from_packet(packet,'.0005')
            def value(row):
                out=P.IntervalTaylor.constant(c,0,row.order)
                for powers,jet in row.terms.items():
                    out=out+jet*c.exp(sum((log*power for log,power in zip(algebra.logs,powers)),c.mpf(0)))
                return out
            def polynomial(poly):
                out=P.IntervalTaylor.constant(c,0,4)
                for power,row in poly.terms.items():out=out+value(row)*R**power
                return out
            quotient=actual.quotients()['p1'];p1=polynomial(quotient.numerator)/polynomial(quotient.denominator)
            dz=lambda row:P.recovery.axial_derivative(row)
            delta=c.mpf('.0005');d=1-z*z;L=1-z*z*delta
            angular=mt*(1-delta/2)-z*dz(mt)*(1-delta)/2-d*dz(mtz)+z*mtz*(2*delta-1)
            mean=mz/R;W=1-z*mean*(1-delta)-d*dz(mean)
            independent=(-W+angular/(u*(c.sqrt(2)*R**c.mpf('1.5'))))*R/L
            for j in range(min(p1.order,independent.order)+1):
                al,ah=P.recovery.endpoints(p1[j]);bl,bh=P.recovery.endpoints(independent[j])
                if ah<bl or bh<al:raise ArithmeticError('Independent full primitive angular input coefficient mismatch')
                count+=1
            if P.recovery.endpoints(p1[0])[0]<=3:raise ArithmeticError('Finite angular barrier fixture failed')
            if algebra.proofs or algebra.final_rows:raise ArithmeticError('Production source factors resolved')
    return count


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Current whole reshape source changed: '+name)
    owner=source.CurrentReshapeRelaxedInputs()
    if source.packets.encode(owner.theorem)!=data['original_current_angular_barrier_source_theorem']:
        raise ValueError('Current angular source/ODE identity changed')
    if source.packets.encode(owner.proof)!=data['current_whole_R110_Rz_relaxed_input_and_bounds']:
        raise ValueError('Current analytic barrier/parameter bound changed')
    if owner.conditions!=data['current_source_function_attachment_conditions']:raise ValueError('Current source attachment conditions changed')
    queries={chart:owner.query(chart,(0,1)) for chart in source.CHARTS}
    queries['reshape_zero']=owner.query('reshape',0);queries['reshape_right']=owner.query('reshape',1)
    if source.packets.encode(queries)!=data['analytic_source_subbox_examples']:raise ValueError('Actual source domain examples changed')
    if queries['reshape_zero']['angular_Q_above_half_on_entire_query'] or queries['reshape']['angular_Q_above_half_on_entire_query']:
        raise ArithmeticError('One-log-unit angular Q bound admitted at inlet')
    if not all(queries[name]['angular_Q_above_half_on_entire_query'] for name in ('inner_reference','reshape_right')):
        raise ArithmeticError('Actual post-one-unit angular source bound missing')
    bad=(('core',0,0),('reshape',-1,0),('reshape',2,0),('inner_reference','inf',0),('reshape','.5',2))
    for chart,coordinate,Z in bad:
        try:owner.query(chart,coordinate,Z)
        except ValueError:pass
        else:raise ArithmeticError('Invalid actual reshape source query admitted')
    if any(data[k] or owner.proof[k] for k in source.OPEN):raise ArithmeticError('Current original relaxed input promoted a whole modified target')
    r=dict(all_passed=True,**{source.GATE:True},**dict.fromkeys(source.OPEN,False),source_family=owner.family,
        current_source_function_attachment_conditions=len(owner.conditions),
        exact_source_AST_bindings=len(owner.theorem['original_source_AST_bindings']),
        exact_full_angular_source_ODE_and_control_identities=len(owner.theorem['exact_identities']),
        independent_full_power_primitive_fixture_coefficients=independent_power_primitive_fixture(),
        positive_current_barrier_margins=len(owner.proof['positive_directed_barrier_margins']),
        actual_source_queries_checked=len(queries),invalid_actual_source_queries_rejected=len(bad),
        angular_Q_not_radial_velocity_Q=True,angular_Q_half_bound_excludes_inlet=True,
        actual_source_radius_and_amplitude_not_materialized=True,source_ancestor_constructors_called=False,
        strict_completed_tensor_cone_new_regions_admitted=0,p1_p2_whole_path_norm_bounds_certified=False,
        axial_restoration_or_active_patch_relaxed_input_certified=False,
        scope=data['scope'],
        input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(r,indent=2)+'\n',encoding='utf8')
    print('Current actual whole reshape/reference analytic barrier and independent primitives PASS',flush=True)
    return r


if __name__=='__main__':run()
