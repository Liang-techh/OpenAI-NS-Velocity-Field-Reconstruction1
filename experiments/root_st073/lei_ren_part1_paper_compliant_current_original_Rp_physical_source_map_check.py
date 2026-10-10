"""Focused direct raw-row physical map checks; no source-owner replay."""
import copy
from dataclasses import replace
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
from types import SimpleNamespace
import sympy as s

import lei_ren_part1_paper_compliant_current_original_Rp_physical_source_map as current
from lei_ren_part1_paper_compliant_global_physical_assembly import (
    cartesian_source_row,time_source_row,angular_polynomial,CS,SN)
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


class PhysicalInterpreter(current.mixed.pulse.radius.RadiusInterpreter):
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings:return self.bindings[i]
        if i in self.memo:return self.memo[i]
        node=self.g.nodes[i]
        if node['operation']=='exact_operator_integer_power':
            assert node['exponent']<0
            assert node['nonzero_certificate']=='original operator denominator; 1-delta*Z^2 > 0'
            value=self.at(node['base'])**node['exponent']
        elif node['operation']=='analytic_unary' and node['name'] in ('sin','cos'):
            value=(s.sin if node['name']=='sin' else s.cos)(self.at(node['argument']))
        else:return super().at(i)
        self.memo[i]=value;return value


def interpretation(owner):
    delta=s.Symbol('same_actual_delta',positive=True)
    field=SimpleNamespace(graph=owner.graph,parameters=owner.radius.parameters,
        contracts=owner.radius.frame.bridge.leading.contracts)
    reader=PhysicalInterpreter(field,False,{
        owner.radius.logRp.node:s.Symbol('same_actual_logRp',real=True),
        owner.radius.functions['mu'].node:s.Symbol('same_actual_mu',positive=True),
        owner.delta_function.node:delta,
        owner.radius.functions['logP'].node:s.Symbol('same_actual_logPstar',real=True)})
    return reader,delta


def powers(row,radial=0):
    # Independent legacy linear pullback uses fixed bookkeeping bases.
    # The four actual source powers are tags, never numeric log values.
    return tuple(s.Rational(v.numerator,v.denominator) for v in row.powers)+(s.Integer(0),s.Rational(radial),s.Integer(0))


def exact_and_signed_rows(owner,sources,views):
    reader,delta=interpretation(owner);c=owner.ctx
    ends=current.mixed.pulse.radius.post.selected.inlet.endpoints
    exact_cache={};spatial=time_rows=terms=groups=pressure=0;coordinates=0
    def equal(left,right):
        key=(left,right)
        if key not in exact_cache:exact_cache[key]=s.simplify(left-right)==0
        assert exact_cache[key],key
    def grouped(items):
        result={}
        for key,value in items:
            key=tuple(Fraction(str(v)) for v in key)
            result[key]=result.get(key,c.mpf(0))+value
        return result
    for chart,view in views.items():
        source=sources[chart];owner.source_view(source,'.521');coords=view['coordinates']
        Z=s.Rational(521,1000);theta=s.Rational(7,10)
        ell=(-10-s.log(1-Z*Z))/2
        equal(reader.at(coords['loglambda']),ell)
        equal(reader.at(coords['log_r']),ell+(s.log(2)+reader.at(coords['logR']))/2)
        equal(reader.at(coords['x']),s.exp(reader.at(coords['log_r']))*s.cos(theta))
        equal(reader.at(coords['y']),s.exp(reader.at(coords['log_r']))*s.sin(theta))
        equal(reader.at(coords['z']),Z*s.exp((1-delta)*ell))
        equal(reader.at(coords['L']),1-delta*Z*Z)
        equal(reader.at(coords['tau']),s.exp(-10));equal(reader.at(coords['t']),1-s.exp(-10))
        equal(s.exp(2*ell)*(1-Z*Z),s.exp(-10))
        coordinates+=9
        zc=c.mpf('0.521');tc=c.mpf('0.7');cs,sn=c.cos(tc),c.sin(tc)
        grids={label:{row.derivative:[(powers(row),row.coefficients[0])]
            for row in source['log_radius_mixed_rows'][name].values()}
            for label,name in current.SOURCE.items()}
        amplitudes={label:{} for label in current.SOURCE}
        for component,rows in view['Cartesian_spatial_rows'].items():
            assert set(rows)=={'x%d_y%d_z%d'%index for index in current.INDICES}
            for key,row in rows.items():
                i,j,b=row.derivative;degree=i+j
                assert type(row) is current.PhysicalSourceRow and row.component==component
                reference=cartesian_source_row(c,grids,component,i,j,b,zc,owner.delta,cs,sn,amplitudes)
                expected_terms=[]
                for (label,a,q),angular in current.cartesian_templates()[component,i,j,b].items():
                    beta=-1 if label==current.UR else -2-2*delta if label==current.P else -1-delta
                    angle=angular.subs({CS:s.cos(theta),SN:s.sin(theta)})
                    for index,expression in current.physical_operators()[a,b].items():
                        expected_terms.append((label,index,angle*2**s.Rational(a-q,2)*
                            expression.subs({current.ZSYM:Z,current.DSYM:delta,current.BSYM:beta}),
                            beta-degree+b*(delta-1)))
                assert len(expected_terms)==len(row.terms)
                for term,(label,index,weight,gamma) in zip(row.terms,expected_terms):
                    base=source['log_radius_mixed_rows'][current.SOURCE[label]]['y%d_Z%d'%index]
                    assert term.source_row is base and term.source_label==label
                    assert term.radial_power==Fraction(-degree,2)
                    assert term.log_scale_parts[:-2] is not None and term.log_scale_parts[:-2]==base.log_scale_parts
                    equal(reader.at(term.operator_function),weight)
                    equal(reader.at(term.lambda_exponent),gamma)
                    equal(reader.at(term.log_scale_parts[-2][1]),-s.Rational(degree,2)*reader.at(coords['logR']))
                    equal(reader.at(term.log_scale_parts[-1][1]),gamma*ell)
                    product=base.coefficients[0]*term.operator_coefficient[0]
                    assert (product-term.signed_coefficient[0]).a<=0<=(product-term.signed_coefficient[0]).b
                    terms+=1
                for label,(legacy,gamma_box) in reference.items():
                    expected=grouped(legacy)
                    actual=grouped((powers(t.source_row,t.radial_power),t.signed_coefficient[0])
                        for t in row.terms if t.source_label==label)
                    assert set(expected)==set(actual)
                    for pp,value in actual.items():
                        difference=value-expected[pp]
                        assert ends(difference)[0]<=0<=ends(difference)[1],(chart,component,key,label)
                        groups+=1
                # Grouping keeps actual signs and includes zero-containing
                # terms. No absolute cap replaces a coefficient.
                assert sum(len(group['contributors']) for group in row.groups())==len(row.terms)
                spatial+=1
        for component,row in view['fixed_x_time_rows'].items():
            assert row.component==component and row.derivative==('t',)
            expected_terms=[];reference={}
            for (label,a,q),angular in current.cartesian_templates()[component,0,0,0].items():
                beta=-1 if label==current.UR else -2-2*delta if label==current.P else -1-delta
                angle=angular.subs({CS:s.cos(theta),SN:s.sin(theta)})
                L=1-delta*Z*Z
                for index,coefficient in (((0,0),-beta/(2*L)),((0,1),(1-delta)*Z/(2*L)),((1,0),1/L)):
                    expected_terms.append((label,index,angle*coefficient,beta-2))
                legacy,gamma=time_source_row(c,grids,label,zc,owner.delta,amplitudes)
                reference[label]=[(pp,value*angular_polynomial(c,angular,cs,sn)) for pp,value in legacy]
            assert len(expected_terms)==len(row.terms)
            for term,(label,index,weight,gamma) in zip(row.terms,expected_terms):
                assert term.source_row is source['log_radius_mixed_rows'][current.SOURCE[label]]['y%d_Z%d'%index]
                assert term.radial_power==0 and term.source_label==label
                equal(reader.at(term.operator_function),weight);equal(reader.at(term.lambda_exponent),gamma)
                equal(reader.at(term.log_scale_parts[-2][1]),0)
                equal(reader.at(term.log_scale_parts[-1][1]),gamma*ell);terms+=1
            for label,legacy in reference.items():
                expected=grouped(legacy)
                actual=grouped((powers(t.source_row),t.signed_coefficient[0]) for t in row.terms if t.source_label==label)
                assert set(expected)==set(actual)
                for pp,value in actual.items():
                    lo,hi=ends(value-expected[pp]);assert lo<=0<=hi;groups+=1
            time_rows+=1
        # This is a directed consistency check; the independent analytic
        # P0 and exact pressure definition are inherited source functions.
        grid=source['log_radius_mixed_rows']
        for n in range(5):
            label='y0_Z'+str(n)
            p,mp,p0=(grid[key][label] for key in ('pressure','Mp','P0'))
            assert p.powers==mp.powers==p0.powers and p.log_scale_parts==mp.log_scale_parts==p0.log_scale_parts
            lo,hi=ends(p.coefficients[0]-mp.coefficients[0]-p0.coefficients[0]);assert lo<=0<=hi
            pressure+=1
        for k in range(1,5):
            for n in range(5-k):
                label='y%d_Z%d'%(k,n);p,mp,p0=(grid[key][label] for key in ('pressure','Mp','P0'))
                assert p.powers==mp.powers and p.log_scale_parts==mp.log_scale_parts
                lo,hi=ends(p.coefficients[0]-mp.coefficients[0]);assert lo<=0<=hi
                assert ends(p0.coefficients[0])==(0,0);pressure+=1
    assert spatial==2100 and time_rows==60 and pressure==225 and coordinates==135
    return dict(passed=True,actual_current_charts=15,Cartesian_spatial_rows=spatial,fixed_x_time_rows=time_rows,
        exact_signed_source_operator_terms=terms,legacy_linear_pullback_signed_group_diagnostics=groups,
        exact_physical_coordinate_identities=coordinates,pressure_Mp_P0_directed_consistency_rows=pressure,
        all_current_full_source_units_and_lazy_scales_preserved=True,
        pressure_value_Pstar2_and_radial_theta2_scales_remain_separate=True,
        no_extra_Ur_sqrt2_or_native_Jacobian_or_factorial=True,
        accepted_original_implicit_fixture_reused_without_source_replay=True,
        finite_correction_N_not_the_transverse_derivative_order=True)


@source_precision
def run(before=None,observed_owner=None,observed_source_views=None,observed_views=None):
    began=time.monotonic();candidate=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    for name,digest in candidate['input_hashes'].items():assert current.sha(name)==digest,name
    owner=current.CurrentOriginalRpPhysicalSourceMap(before,require_checked=False)
    assert not any(candidate[k] for k in current.GATES+current.OPEN)
    assert candidate['source_family']==owner.family_record and candidate['actual_source_graph']==owner.assert_graph()
    assert candidate['original_linear_operator_definitions']==owner.definitions
    assert candidate['original_N_definition']==owner.original_N_definition
    assert candidate['accepted_independent_original_physical_fixture']==owner.canonical_fixture
    if observed_owner is None:
        sources={chart:owner.transport.evaluate(chart,'.521',coordinate) for _,chart,coordinate in current.mixed.pulse.VIEWS}
        views={chart:owner.map_source_view(source,'.521') for chart,source in sources.items()}
    else:
        assert type(observed_owner) is type(owner) and observed_owner.before is before
        assert observed_owner.transport is owner.transport and observed_owner.graph is owner.graph
        assert observed_owner.hashes==owner.hashes and all(observed_owner.assert_graph().values())
        sources,views=observed_source_views,observed_views
    assert set(sources)==set(views)==set(current.mixed.CHARTS)
    graph=candidate['exact_physical_expression_graph']
    assert graph==owner.graph.nodes[:len(graph)]
    for chart,view in views.items():
        assert current.mixed.pulse.raw.packed(current.report(view))==candidate['actual_current_fifteen_chart_physical_views'][chart]
        assert not any(view[k] for k in current.GATES+current.OPEN)
    exact=exact_and_signed_rows(owner,sources,views)
    rejected=[];source=sources['flatten']
    for name,action in (
        ('axial_boundary_not_finite_coordinate',lambda:owner.map_source_view(source,'1')),
        ('wrong_axial_point',lambda:owner.map_source_view(source,'.5')),
        ('float_native_coordinate',lambda:owner.coordinates(source['geometry'],Fraction(521,1000),-10.0,'7/10')),
        ('interval_log_time_is_not_an_exact_function',lambda:owner.coordinates(source['geometry'],Fraction(521,1000),owner.ctx.mpf(-10),'7/10'))):
        try:action()
        except (ValueError,TypeError):rejected.append(name)
        else:raise AssertionError('Wrong physical input accepted: '+name)
    foreign=current.mixed.pulse.radius.FunctionRef(copy.deepcopy(owner.graph),0)
    try:owner.coordinates(source['geometry'],Fraction(521,1000),foreign,'7/10')
    except ValueError:rejected.append('foreign_log_time_graph')
    else:raise AssertionError('Foreign time graph accepted')
    for name,change in (
        ('native_Jacobian_instead_of_ordinary_row',lambda row:replace(row,powers=row.powers[:-1]+(Fraction(1),))),
        ('wrong_pressure_radial_units',lambda row:replace(row,powers=(Fraction(0),Fraction(0),Fraction(2),Fraction(0))))):
        bad=dict(source);bad['log_radius_mixed_rows']=dict(source['log_radius_mixed_rows'])
        component='Ur' if name.startswith('native') else 'pressure'
        bad['log_radius_mixed_rows'][component]=dict(source['log_radius_mixed_rows'][component])
        bad['log_radius_mixed_rows'][component]['y1_Z0']=change(bad['log_radius_mixed_rows'][component]['y1_Z0'])
        try:owner.map_source_view(bad,'.521')
        except ValueError:rejected.append(name)
        else:raise AssertionError('Wrong physical source accepted: '+name)
    bad=copy.copy(owner);bad.delta_function=owner.graph.constant(0)
    try:bad.assert_graph()
    except ValueError:rejected.append('delta_cap_or_zero_as_function_definition')
    else:raise AssertionError('Changed defining delta accepted')
    # An exact time correlated with the same enormous source radius can
    # cancel before numeric arithmetic. No finite logR endpoint is needed.
    reader,delta=interpretation(owner);z=Fraction(521,1000)
    correlated=owner.coordinates(source['geometry'],z,owner.graph.neg(views['flatten']['coordinates']['logR']),'7/10')
    assert s.simplify(reader.at(correlated['log_r'])-(s.log(2)-s.log(1-s.Rational(521,1000)**2))/2)==0
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        original_linear_operator_definitions=owner.definitions,
        original_source_and_divergence_identities=owner.canonical_source_identities,
        accepted_independent_original_physical_fixture=owner.canonical_fixture,
        direct_actual_raw_row_physical_operator_transfer=exact,
        exact_same_graph_radius_time_cancellation_checked=True,
        actual_same_source_typed_observations_not_receipt_scalar_rows_used=True,
        original_N_definition=owner.original_N_definition,rejected_sources_and_domains=rejected,
        **dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.mixed.pulse.raw.packed(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_PHYSICAL_SOURCE_MAP 2100 Cartesian + 60 fixed-x time rows',flush=True)
    return result


if __name__=='__main__':run()
