"""Independent R100 pressure/velocity/moment normalization and provenance."""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_R100_endpoint as current
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

fields=current.fields;ep=current.ep


def finite(row):
    return row.ctx.mpf(0) if row.zero else row.coefficient*row.ctx.exp(row.scale.evaluate())


def contains(row,value):
    lo,hi=ep(finite(row));assert lo<=value<=hi,(mp.nstr(lo,20),mp.nstr(value,20),mp.nstr(hi,20))


def fixtures():
    x=sy.symbols('x');R=sy.Integer(100);F0=sy.Rational(2,7);P2=sy.Rational(5,4);delta=sy.Rational(1,13)
    g=x/7-x**2/11+x**3/13-x**4/17+x**5/19
    exponential=sy.series(sy.exp(g),x,0,6).removeO()
    exponential2=sy.series(sy.exp(2*g),x,0,6).removeO()
    polynomials={name:sy.Rational(j+3,11)+sum(sy.Rational((-1)**(j+n)*(j+1),37*(n+1))*x**n
            for n in range(1,6)) for j,name in enumerate(('phi','V','H','M','K','A','B','C'))}
    c=MPIntervalContext();c.dps=90;comparisons=velocity=pressure=moments=0;first_case=None
    def numeric(value):return mp.mpf(str(sy.N(value,125)))
    def jet(poly):return IntervalTaylor(c,[c.mpf(str(sy.N(sy.expand(poly).coeff(x,n),125))) for n in range(6)])
    with mp.workdps(130):
        for Z in (sy.Integer(0),sy.Rational(1,2),sy.Integer(1),sy.Integer(-1)):
            logRa=c.ln(c.mpf('.7'));flow=fields.MacroFlow(c,c.ln(c.mpf('.0002')),c.ln(c.mpf(5)/4),
                2*c.ln(c.mpf(2)/7),logRa,c.ln(100)-logRa,c.mpf('.1'))
            source_fields={name:flow.jet(jet(polynomials[name])) for name in ('phi','V')}
            histories={name:flow.jet(jet(polynomials[name])) for name in current.moments.RATES}
            p0=jet(1+2*x/9-x**2/10+x**3/12-x**4/15+x**5/18)
            kwargs=dict(flow=flow,Z=c.mpf(str(sy.N(Z,125))),delta=c.mpf(1)/13,p0=p0,
                amplitude_ratios=jet(exponential),amplitude_squared_ratios=jet(exponential2),
                actual_fields=source_fields,actual_histories=histories)
            result=current.recover_endpoint(**kwargs)
            z=Z+x;M=polynomials['M'];V=polynomials['V'];phi=polynomials['phi']
            Q=sy.series((2*z*V-z*M*(1-delta)-(1-z*z)*sy.diff(M,x))/(1-z*z*delta),x,0,5).removeO()
            expressions=dict(Ur=sy.sqrt(R/2)*Q,Utheta=sy.sqrt(2*R)*F0*phi*exponential,Uz=V)
            for name,expr in expressions.items():
                for n,row in enumerate(result['physical_velocity_axial_coefficients'][name]):
                    contains(row,numeric(sy.expand(expr).coeff(x,n)));velocity+=1
            p0expr=1+2*x/9-x**2/10+x**3/12-x**4/15+x**5/18
            inc=R*F0**2*polynomials['C']*exponential2
            for key,expr in (('physical_pressure_axis_axial5',P2*p0expr),
                             ('physical_pressure_radial_increment_axial5',inc),
                             ('physical_total_pressure_axial5',P2*p0expr+inc)):
                for n,row in enumerate(result[key]):contains(row,numeric(sy.expand(expr).coeff(x,n)));pressure+=1
            exprs=dict(Mtheta=F0*R**2*polynomials['H']*exponential,Mz=R*M,
                Mtheta_z=F0*R**2*polynomials['K']*exponential,
                Mztheta=R*polynomials['A']-R**2*F0**2*polynomials['B']*exponential2,
                Mp=inc)
            for name,expr in exprs.items():
                for n,row in enumerate(result['physical_cumulative_moment_axial5'][name]):
                    contains(row,numeric(sy.expand(expr).coeff(x,n)));moments+=1
            for n,row in enumerate(result['actual_Q_axial4_coefficients']):
                contains(row,numeric(Q.coeff(x,n)));comparisons+=1
            assert len(result['actual_Q_axial4_coefficients'])==5
            assert all(len(row)==6 for row in result['physical_cumulative_moment_axial5'].values())
            if first_case is None:first_case=kwargs
    return dict(passed=True,independent_radial_Q_Taylor_comparisons=comparisons,
        independent_physical_velocity_Taylor_comparisons=velocity,
        independent_separate_and_total_pressure_Taylor_comparisons=pressure,
        independent_physical_cumulative_moment_Taylor_comparisons=moments,
        midplane_interior_and_both_axis_boundaries_checked=True,
        nonzero_all_orders_and_amplitude_derivatives_used=True),first_case


def same_row(flow,saved,row):
    restored=current.restore_row(flow,saved)
    assert restored.scale.powers==row.scale.powers
    assert restored.scale.offset._mpi_==row.scale.offset._mpi_
    assert restored.coefficient._mpi_==row.coefficient._mpi_


def genuine():
    owner=current.OriginalR100Endpoint();saved=json.loads((current.HERE/current.NAME).read_bytes())
    imports=Qorders=pressurejoins=amplitudeguards=0
    for label in ('0','.5'):
        with mp.workdps(owner.c.dps+40):flow,proof,value=owner.owner(label)
        raw=owner.saved['packets'][label][-1]['evaluation']
        for name,row in value['normalized_actual_fields'].items():
            for old,v in zip(raw['same_original_macro_field_functions'][name],row):same_row(flow,old,v);imports+=1
        for name,row in value['normalized_actual_own_six_moments'].items():
            for old,v in zip(raw['actual_six_moment_functions'][name],row):same_row(flow,old,v);imports+=1
        assert value['radius'].scale.powers==(0,0,0,0,0) and ep(value['radius'].coefficient)==(100,100)
        assert len(value['actual_Q_axial4_coefficients'])==5;Qorders+=1
        for old,v in zip(proof['original_P0_coefficients'],value['pressure_axis_over_Pstar_squared_axial5']):
            assert v.scale.powers==(0,0,0,0,0) and v.coefficient._mpi_==old._mpi_;pressurejoins+=1
        # Actual G defines the log amplitude; no Gbar cap or value selector.
        original=owner.upstream.fields.records['anchored_axis_amplitude']['anchored_amplitude_packets'][label]
        want=2*fields.previous.read_interval(owner.c,original['logF0'])
        assert flow.logs[2]._mpi_==want._mpi_;amplitudeguards+=1
        encoded=saved['packets'][label]['function_evaluation']
        for key in ('actual_Q_axial4_coefficients','physical_pressure_axis_axial5',
                    'physical_pressure_radial_increment_axial5','F_actual_true_axial5_divided_by_F0'):
            for old,v in zip(encoded[key],value[key]):same_row(flow,old,v)
        assert saved['packets'][label]['source_family']==owner.family
        assert value['analytic_P0_kept_separate'] and value['amplitude_derivatives_dressed_once']
        assert value['physical_amplitudes_not_materialized']
    return dict(passed=True,exact_original_macro_endpoint_rows_reused=imports,
        exact_original_P0_coefficients_preserved=pressurejoins,
        actual_anchored_F0_squared_source_identities=amplitudeguards,
        Q_order4_only_guards=Qorders,whole_axis_micro_switch_R110_and_recursion_remain_open=True),owner


def guards(case,owner):
    count=0
    def reject(call):
        nonlocal count
        try:call()
        except (ValueError,TypeError,KeyError):count+=1
        else:raise AssertionError('Invalid endpoint request admitted')
    for key,val in (('Z',2),('delta',1),('delta',-1),('p0',IntervalTaylor.constant(case['flow'].c,1,4)),
                    ('amplitude_ratios',IntervalTaylor.constant(case['flow'].c,2,5))):
        reject(lambda key=key,val=val:current.recover_endpoint(**dict(case,**{key:val})))
    fields2=dict(case['actual_fields']);fields2['phi']=fields2['phi'][:5]
    reject(lambda:current.recover_endpoint(**dict(case,actual_fields=fields2)))
    histories=dict(case['actual_histories']);histories.pop('B')
    reject(lambda:current.recover_endpoint(**dict(case,actual_histories=histories)))
    reject(lambda:owner.owner('whole_Z'))
    record=owner.saved['packets']['0'][-1]['evaluation']['same_original_macro_field_functions']['phi'][0]
    for replacement in (dict(point_value_selected=True),dict(encloses_original_source_function=False),
                        dict(exact_zero=not record['exact_zero'])):
        reject(lambda replacement=replacement:current.restore_row(case['flow'],dict(record,**replacement)))
    other_c=MPIntervalContext();other_c.dps=90
    reject(lambda:current.recover_endpoint(**dict(case,p0=IntervalTaylor.constant(other_c,1,5))))
    saved_fraction=owner.saved['packets']['0'][-1]['evaluation']['geometry']['fraction']
    owner.owners.pop('0');owner.saved['packets']['0'][-1]['evaluation']['geometry']['fraction']=[1,2]
    reject(lambda:owner.owner('0'))
    owner.saved['packets']['0'][-1]['evaluation']['geometry']['fraction']=saved_fraction
    return dict(passed=True,invalid_geometry_ratio_pressure_history_or_selected_value_requests_rejected=count)


def run():
    began=time.monotonic();fixture,case=fixtures();print('Independent R100 physical normalization PASS',flush=True)
    native,owner=genuine();guard=guards(case,owner)
    report=json.loads((current.HERE/current.NAME).read_bytes())
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    fields.previous.bind(owner.hashes,current.NAME,current.sha(current.NAME))
    fields.previous.bind(owner.hashes,Path(__file__).name,current.sha(Path(__file__).name))
    assert report[current.GATE] and report['source_bindings']['passed']
    assert all(report[key] is False for key in fields.previous.OPEN)
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_normalization=fixture,genuine_source_binding=native,guards=guard,
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(fields.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Same-source factored R100 endpoint and original pressure PASS',flush=True);return result


if __name__=='__main__':run()
