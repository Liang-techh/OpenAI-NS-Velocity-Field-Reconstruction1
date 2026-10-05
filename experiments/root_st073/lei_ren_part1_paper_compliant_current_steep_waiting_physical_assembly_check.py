"""Focused current O7 physical source admission; retain prior operator fixtures."""
import ast
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_steep_waiting_physical_assembly import (
    CurrentSteepWaitingPhysicalAssembly, CHARTS, GATE, COMPOSED, SOURCE_GATE,
    UNIFORM, SCOPES, OPEN, HERE, PREFIX, sha)
from lei_ren_part1_paper_compliant_global_physical_assembly import (
    INDICES,COMPONENTS,UZ,UT,UR,P,cartesian_templates,POST)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_steep_waiting_physical_assembly.json'


def actual_dispatch_call_bindings():
    tree=ast.parse((HERE/(PREFIX+'current_steep_waiting_physical_assembly.py')).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='_CurrentSteepWaitingPhysicalDispatch')
    for method,expression in (('provider','self.current.provider(chart)'),
            ('evaluate','self.current.evaluate(chart,Z,coordinate)')):
        fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name==method)
        if len(fn.body)!=1 or not isinstance(fn.body[0],ast.Return) or (
                ast.dump(fn.body[0].value)!=ast.dump(ast.parse(expression,mode='eval').body)):
            raise ValueError('Physical dispatcher must call same checked current source: '+method)
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CurrentSteepWaitingPhysicalAssembly')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='evaluate')
    expression='BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=False)'
    if sum(ast.dump(n)==ast.dump(ast.parse(expression,mode='eval').body)
            for n in ast.walk(fn) if isinstance(n,ast.Call))!=1:
        raise ValueError('Original full physical evaluate call required')
    if not {'steep_entry','steep_power','steep_exit','waiting'}<=set(POST):raise ValueError('Original post fixed-unit route differs')
    return dict(current_provider_and_evaluate_same_source_AST_bound=True,
        actual_unchanged_full_physical_operator_call_AST_bound=True,
        four_current_O7_charts_use_original_post_fixed_unit_route=True,passed=True)


def independent_post_radius_and_fixed_unit_fixture():
    """Differentiate finite full fields; exercise four actual O7 routes."""
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=110;tol=mp.mpf('1e-65')
        mu=mp.mpf('.23');z=mp.mpf('.31');L=mp.mpf('18.4');T=mp.mpf('7.3');W=mp.mpf('5.7')
        field=object.__new__(CurrentSteepWaitingPhysicalAssembly)
        field.ctx=c;field.params=SimpleNamespace(mu=c.mpf(mu));field.logP=c.ln(2);field.logRp=c.ln(4)
        flat=SimpleNamespace(logEv2_parts=dict(inlet_log=2*c.ln(3),inverse_mu_term=-13/c.mpf(mu),finite_offset=c.mpf(-26)))
        steep=SimpleNamespace(outer=SimpleNamespace(Lrel=c.mpf(L)),Ts=c.mpf(T),wait=c.mpf(W))
        field.dispatch=SimpleNamespace(provider=lambda chart:flat)
        ev=2*3*mp.exp(-13/(2*mu)-13)
        functions={UT:lambda y,z:ev*mp.exp(-mp.mpf('.19')*y)/(1+z*z),
            UZ:lambda y,z:ev*(mp.sin(z)+mp.mpf('.02')*y*y*z**3),
            UR:lambda y,z:ev*(mp.cos(z)+mp.mpf('.03')*y*z*z),
            P:lambda y,z:4*(1+z**4+mp.sin(y)*z**3)}
        units={UT:ev,UZ:ev,UR:ev,P:mp.mpf(4)};counts={}
        for chart,value in (('steep_entry',mp.mpf('.37')),('steep_power',mp.mpf('.23')),
                ('steep_exit',mp.mpf('.63')),('waiting',mp.mpf('.57'))):
            offset=(100+L+value if chart=='steep_entry' else 100+L+1+T*value if chart=='steep_power'
                else 100+L+1+T+value if chart=='steep_exit' else 100+L+2+T+W*value)
            y=mp.log(4)+13/mu+offset;raw={};targets={}
            for label,fn in functions.items():
                raw[label]={}
                for k in range(5):
                    for n in range(5-k):
                        target=mp.diff(fn,(y,z),(k,n));targets[label,k,n]=target
                        v=target/units[label];raw[label]['y'+str(k)+'_Z'+str(n)]=c.mpf([v-tol,v+tol])
            packet={'physical_mixed_derivatives_total_order_le4':raw}
            logR,_=field.radius(chart,c.mpf(value),packet,steep)
            with mp.workdps(160):
                reference_offset=(100+L+value if chart=='steep_entry' else 100+L+1+T*value if chart=='steep_power'
                    else 100+L+1+T+value if chart=='steep_exit' else 100+L+2+T+W*value)
                reference=mp.log(4)+13/mu+reference_offset
            if not endpoints(logR)[0]<=reference<=endpoints(logR)[1]:raise ValueError('Independent current POST radius differs')
            grids,logs,amplitudes=field.normalized_sources(chart,c.mpf(z),c.mpf(value),packet,steep,logR)
            count=0
            for label in functions:
                unit=c.exp(sum((logs[i]*power for i,power in amplitudes[label].items()),c.mpf(0)))
                for (k,n),terms in grids[label].items():
                    if len(terms)!=1 or any(terms[0][0]):raise ValueError('Full ordinary POST derivative normalization differs')
                    lo,hi=endpoints(terms[0][1]*unit)
                    if not lo<=targets[label,k,n]<=hi:raise ValueError('Independent full POST fixed-unit restoration differs')
                    count+=1
            counts[chart]=count
        return dict(independent_full_field_mixed4_rows_by_current_POST_chart=counts,
            independent_full_field_mixed4_rows_checked=sum(counts.values()),
            four_original_O7_radius_branches_checked=True,nonzero_synthetic_all_velocity_components_checked=True,
            fixed_Ev0_velocity_and_absolute_Pstar_squared_pressure_checked=True,
            finite_fixture_only=True,actual_Md40_point_admission=False,passed=True)


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        field=CurrentSteepWaitingPhysicalAssembly(require_checked=False)
        if encode(pack(field.manifest()))!={k:v for k,v in raw.items() if k!='whole_current_steep_waiting_physical_maps'}:
            raise ValueError('Fresh current O7 physical graph or retained evidence differs')
        if (not all(raw['current_provider_graph_identity'].values())
                or not all(raw['original_physical_operators_retained'].values())
                or not raw['native_current_parameter_source_bridge']['passed']
                or not raw['current_post_fixed_unit_source_binding']['passed']
                or not raw['current_steep_waiting_radius_functional_binding']['passed']
                or not raw['current_steep_waiting_profile_source_binding']['passed']
                or raw[GATE] or raw[COMPOSED] or raw['current_steep_waiting_physical_owner_installed']
                or not raw[SOURCE_GATE]
                or tuple(raw['current_physical_chart_owners'])!=CHARTS
                or raw[UNIFORM] or raw['full_pulse_C4_installed'] or any(raw[k] for k in SCOPES+OPEN)):
            raise ValueError('Current O7 physical source scope differs')
        calls=actual_dispatch_call_bindings();c=MPIntervalContext();c.dps=240
        indices={'x'+str(i)+'_y'+str(j)+'_z'+str(k) for i,j,k in INDICES}
        groups=raw['whole_current_steep_waiting_physical_maps'];counts={};terms=0;zeros=0
        if set(groups)!={'steep_entry','steep_power','steep_exit','waiting'}:
            raise ValueError('Four current O7 physical source owners required')
        if any(set(views)!={'inlet','whole_domain','exit'} for views in groups.values()):
            raise ValueError('Original O7 inlet/domain/exit coverage required')
        packets={chart+'_'+view:packet for chart,views in groups.items() for view,packet in views.items()}
        expected={chart+'_'+view:chart for chart in groups for view in ('inlet','whole_domain','exit')}
        if set(packets)!=set(expected):raise ValueError('Whole steep/waiting physical coverage omitted')
        def check(row,label):
            nonlocal terms,zeros
            if bool(row['exact_zero'])!=(row['log_absolute_upper'] is None):raise ValueError('Source zero/bound differs')
            if row['exact_zero']:
                if row['terms']:raise ValueError('Exact zero has source terms')
                zeros+=1;return
            if label in (UR,UZ):raise ValueError('Current O7 meridional source was not exactly zero')
            lo,hi=endpoints(read_interval(c,row['log_absolute_upper']))
            if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):raise ValueError('Nonfinite physical log bound')
            if (not row['positive_source_exponentials_not_materialized']
                    or row['source_row_mode']!='provider_prebounded_mixed_rows'
                    or endpoints(read_interval(c,row['physical_lambda_exponent']))[1]>=0):
                raise ValueError('Original fixed-unit/lambda scope differs')
            for term in row['terms']:
                powers=[mp.make_mpf(tuple(v['exact_mpf_tuple'])) for v in term['source_log_exponents']]
                if len(powers)!=7 or powers[5]>0:raise ValueError('Original physical source bases/radial factors omitted')
                if label==P:
                    if powers[1]!=2 or powers[4]!=0:raise ValueError('Absolute pressure Pstar squared omitted')
                elif powers[4]!=1 or powers[1]!=0:raise ValueError('Fixed Ev0 velocity unit omitted')
                lo,hi=endpoints(read_interval(c,term['signed_coefficient']))
                if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):raise ValueError('Nonfinite signed coefficient')
                terms+=1
        for name,packet in packets.items():
            chart=expected[name]
            if (packet['chart']!=chart or packet[GATE] or packet[COMPOSED]
                    or packet['current_steep_waiting_physical_owner_installed'] or not packet[SOURCE_GATE]
                    or packet['current_source_owner']!=field.source_owners[chart]['provider']
                    or packet['current_source_acceptance_receipt']!=field.source_owners[chart]['acceptance_receipt']
                    or packet['datum_enclosure_sha256']!=field.datum_sha
                    or not packet['current_steep_provider_used_for_derivative_grid_and_radius']
                    or not packet['same_current_flatten_provider_used_for_fixed_unit']
                    or not packet['original_fixed_Ev0_velocity_and_Pstar_squared_pressure_units_restored']
                    or any(packet[k] for k in SCOPES+OPEN)):
                raise ValueError('Current physical O7 source owner/scope differs')
            spatial=packet['physical_spatial_cartesian_mixed4'];time=packet['first_fixed_x_physical_time_derivative'];count=0
            if set(spatial)!=indices or set(time)!=set(COMPONENTS):raise ValueError('Spatial4/time1 coverage omitted')
            for index,components in spatial.items():
                i,j,b=(int(v[1:]) for v in index.split('_'))
                if set(components)!=set(COMPONENTS):raise ValueError('Cartesian component omitted')
                for component,parts in components.items():
                    expected_labels={label for label,a,q in cartesian_templates()[component,i,j,b]}
                    if set(parts)!=expected_labels:raise ValueError('Moving-basis contribution omitted')
                    for label,row in parts.items():check(row,label);count+=1
            for component,parts in time.items():
                expected_labels={UR,UT} if component in ('ux','uy') else {UZ} if component=='uz' else {P}
                if set(parts)!=expected_labels:raise ValueError('Fixed-position time contribution omitted')
                for label,row in parts.items():check(row,label);count+=1
            if count!=216 or not packet['moving_cylindrical_basis_differentiated'] or not packet['normalization_not_differentiated_twice']:
                raise ValueError('Original physical operators/units omitted')
            counts[name]=count
        for chart in ('core','bridge_first','heat_collar','heat_exterior'):
            try:field.evaluate(chart,0,0)
            except ValueError:pass
            else:raise ValueError('Unadmitted physical owner accepted')
        # The same original post normalization/operator serves flatten and
        # all four O7 charts. Its source-hashed independent60-row fixture and
        # Cartesian/time fixture remain unchanged; do not rerun them.
        fixture=raw['reused_independent_post_fixed_unit_fixture']
        if (not fixture['passed'] or fixture['independently_differentiated_full_flatten_mixed4_rows']!=60
                or not raw['reused_independent_coordinate_fixture']['passed']):
            raise ValueError('Retained independent operator/unit evidence omitted')
        post_fixture=independent_post_radius_and_fixed_unit_fixture()
        old=raw['retained_twenty_three_chart_physical_evidence']
        result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
            datum_enclosure_sha256=field.datum_sha,current_physical_chart_owners_checked=27,
            current_steep_waiting_spatial_and_time_source_contributions_checked=counts,
            total_new_steep_waiting_source_contributions_checked=sum(counts.values()),
            twenty_seven_regular_source_contributions_checked=old['total_regular_source_contributions']
                +sum(counts[chart+'_whole_domain'] for chart in groups),
            retained_twenty_three_chart_physical_evidence=old,retained_signed_source_terms_checked=terms,
            exact_zero_source_contributions_checked=zeros,actual_current_dispatch_call_bindings=calls,
            independent_current_O7_radius_and_fixed_unit_fixture=post_fixture,
            reused_independent_post_fixed_unit_fixture=fixture,
            reused_independent_coordinate_fixture=raw['reused_independent_coordinate_fixture'],
            **{GATE:True,COMPOSED:True,SOURCE_GATE:True,UNIFORM:False},current_steep_waiting_physical_owner_installed=True,
            full_pulse_C4_installed=False,**dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),all_passed=True,
            input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS current27 physical owners: steep/waiting spatial4/time1, 2592 new source contributions; prior23 evidence retained',flush=True)
    return result


if __name__=='__main__':run()
