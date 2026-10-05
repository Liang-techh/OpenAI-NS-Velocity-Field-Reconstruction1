"""Current flatten fixed-unit, Cartesian spatial4/time1 focused acceptance."""
import ast
import json
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext

from lei_ren_part1_paper_compliant_current_flatten_physical_assembly import (
    CurrentFlattenPhysicalAssembly, CHARTS, GATE, COMPOSED,
    UNIFORM, SCOPES, OPEN, HERE, PREFIX, sha)
from lei_ren_part1_paper_compliant_global_physical_assembly import (
    INDICES,COMPONENTS,UZ,UT,UR,P,cartesian_templates)
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

NAME=PREFIX+'current_flatten_physical_assembly.json'


def actual_dispatch_call_bindings():
    tree=ast.parse((HERE/(PREFIX+'current_flatten_physical_assembly.py')).read_text(encoding='utf8'))
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='_CurrentFlattenPhysicalDispatch')
    for method,expression in (('provider','self.current.provider(chart)'),
            ('evaluate','self.current.evaluate(chart,Z,coordinate)')):
        fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name==method)
        if len(fn.body)!=1 or not isinstance(fn.body[0],ast.Return) or (
                ast.dump(fn.body[0].value)!=ast.dump(ast.parse(expression,mode='eval').body)):
            raise ValueError('Current physical dispatcher call must use the accepted same source: '+method)
    cls=next(n for n in tree.body if isinstance(n,ast.ClassDef) and n.name=='CurrentFlattenPhysicalAssembly')
    fn=next(n for n in cls.body if isinstance(n,ast.FunctionDef) and n.name=='evaluate')
    expression='BASE.evaluate(self,chart,Z,coordinate,log_tau=log_tau,theta=theta,axis=False)'
    if sum(ast.dump(n)==ast.dump(ast.parse(expression,mode='eval').body)
            for n in ast.walk(fn) if isinstance(n,ast.Call))!=1:
        raise ValueError('Unchanged original full physical evaluate call required')
    return dict(current_provider_and_evaluate_same_source_AST_bound=True,
        actual_unchanged_full_physical_operator_call_AST_bound=True,passed=True)


def independent_flatten_fixed_unit_fixture():
    """Differentiate finite full fields first, then restore frozen native units."""
    with mp.workdps(85):
        c=MPIntervalContext();c.dps=110;tol=mp.mpf('1e-65')
        mu=mp.mpf('.23');z=mp.mpf('.31');t=mp.mpf('2.7')
        field=object.__new__(CurrentFlattenPhysicalAssembly)
        field.ctx=c;field.params=SimpleNamespace(mu=c.mpf(mu));field.logP=c.ln(2);field.logRp=c.ln(4)
        flat=SimpleNamespace(logEv2_parts=dict(inlet_log=2*c.ln(3),inverse_mu_term=-13/c.mpf(mu),finite_offset=c.mpf(-26)))
        field.dispatch=SimpleNamespace(provider=lambda chart:flat)
        ev=2*3*mp.exp(-13/(2*mu)-13);y=mp.log(4)+13/mu+t
        functions={UT:lambda y,z:ev*mp.exp(-mp.mpf('.19')*y)*(1+z*z)**-1,
            UZ:lambda y,z:ev*(mp.sin(z)+mp.mpf('.02')*y*y*z**3),
            UR:lambda y,z:ev*(mp.cos(z)+mp.mpf('.03')*y*z*z),
            P:lambda y,z:4*(1+z**4+mp.sin(y)*z**3)}
        units={UT:ev,UZ:ev,UR:ev,P:mp.mpf(4)};targets={};raw={}
        for label,fn in functions.items():
            raw[label]={}
            for k in range(5):
                for n in range(5-k):
                    target=mp.diff(fn,(y,z),(k,n));targets[label,k,n]=target
                    v=target/units[label]
                    raw[label]['y'+str(k)+'_Z'+str(n)]=c.mpf([v-tol,v+tol])
        packet={'physical_mixed_derivatives_total_order_le4':raw}
        logR,_=field.radius('flatten',c.mpf(t),packet,flat)
        with mp.workdps(160):reference_y=mp.log(4)+13/mu+t
        if not endpoints(logR)[0]<=reference_y<=endpoints(logR)[1]:raise ArithmeticError('Independent Rv radius mismatch')
        grids,logs,amplitudes=field.normalized_sources('flatten',c.mpf(z),c.mpf(t),packet,flat,logR)
        count=0
        for label in functions:
            unit=c.exp(sum((logs[i]*power for i,power in amplitudes[label].items()),c.mpf(0)))
            for (k,n),terms in grids[label].items():
                if len(terms)!=1 or any(terms[0][0]):raise ValueError('Ordinary full derivative normalization changed')
                restored=terms[0][1]*unit
                if not endpoints(restored)[0]<=targets[label,k,n]<=endpoints(restored)[1]:
                    raise ArithmeticError('Full flatten derivative/fixed Ev0 or pressure unit mismatch')
                count+=1
        return dict(independently_differentiated_full_flatten_mixed4_rows=count,
            nonzero_synthetic_all_velocity_components_checked=True,
            fixed_full_Ev0_velocity_units_checked=True,absolute_pressure_Pstar_squared_units_checked=True,
            source_Rv_radius_checked=True,finite_fixture_only=True,passed=True)


def run():
    with mp.workdps(300):
        raw=json.loads((HERE/NAME).read_bytes());_verify_hashes(raw)
        field=CurrentFlattenPhysicalAssembly(require_checked=False)
        manifest=encode(pack(field.manifest()))
        if manifest!={k:v for k,v in raw.items() if k!='whole_current_flatten_physical_maps'}:
            raise ValueError('Fresh current flatten source graph, parameters or retained evidence differs')
        if (not all(raw['current_provider_graph_identity'].values())
                or not all(raw['original_physical_operators_retained'].values())
                or not raw['native_current_parameter_source_bridge']['passed']
                or not raw['current_flatten_fixed_unit_source_binding']['passed']
                or raw[GATE] or raw[COMPOSED] or raw['current_flatten_physical_owner_installed']
                or not raw['current_pulse_terminal_flatten_source_ownership_certified']
                or tuple(raw['current_physical_chart_owners'])!=CHARTS
                or raw[UNIFORM] or raw['full_pulse_C4_installed'] or any(raw[k] for k in SCOPES+OPEN)):
            raise ValueError('New physical flatten source scope differs')
        calls=actual_dispatch_call_bindings();c=MPIntervalContext();c.dps=240
        indices={'x'+str(i)+'_y'+str(j)+'_z'+str(k) for i,j,k in INDICES}
        packets=raw['whole_current_flatten_physical_maps'];counts={};terms=0;zeros=0
        if set(packets)!={'inlet','whole_domain','exit'}:raise ValueError('Original full flatten coverage omitted')
        def check(row,label):
            nonlocal terms,zeros
            if bool(row['exact_zero'])!=(row['log_absolute_upper'] is None):raise ValueError('Source zero/bound differs')
            if row['exact_zero']:
                if row['terms']:raise ValueError('Exact zero has source terms')
                zeros+=1;return
            if label in (UR,UZ):raise ValueError('Flatten meridional source was not exactly zero')
            lo,hi=endpoints(read_interval(c,row['log_absolute_upper']))
            if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):raise ValueError('Nonfinite physical log bound')
            if (not row['positive_source_exponentials_not_materialized']
                    or row['source_row_mode']!='provider_prebounded_mixed_rows'
                    or endpoints(read_interval(c,row['physical_lambda_exponent']))[1]>=0):
                raise ValueError('Original physical fixed-unit/lambda scope differs')
            for term in row['terms']:
                powers=[mp.make_mpf(tuple(v['exact_mpf_tuple'])) for v in term['source_log_exponents']]
                if len(powers)!=7 or powers[5]>0:raise ValueError('Full physical source bases/radial derivative factors omitted')
                if label==P:
                    if powers[1]!=2 or powers[4]!=0:raise ValueError('Absolute pressure Pstar squared unit omitted')
                elif powers[4]!=1 or powers[1]!=0:raise ValueError('Fixed Ev0 velocity unit omitted')
                lo,hi=endpoints(read_interval(c,term['signed_coefficient']))
                if not lo<=hi or not mp.isfinite(lo) or not mp.isfinite(hi):raise ValueError('Nonfinite signed source coefficient')
                terms+=1
        for name,packet in packets.items():
            if (packet['chart']!='flatten' or packet[GATE] or packet[COMPOSED]
                    or packet['current_flatten_physical_owner_installed']
                    or packet['current_source_owner']!=field.source_owners['flatten']['provider']
                    or packet['current_source_acceptance_receipt']!=field.source_owners['flatten']['acceptance_receipt']
                    or packet['datum_enclosure_sha256']!=field.datum_sha
                    or not packet['current_flatten_provider_used_for_derivative_grid_and_unit']
                    or not packet['original_fixed_Ev0_velocity_and_Pstar_squared_pressure_units_restored']
                    or any(packet[k] for k in SCOPES+OPEN)):
                raise ValueError('Current physical flatten owner or unfinished scope differs')
            spatial=packet['physical_spatial_cartesian_mixed4'];time=packet['first_fixed_x_physical_time_derivative'];count=0
            if set(spatial)!=indices or set(time)!=set(COMPONENTS):raise ValueError('Spatial4/time1 coverage omitted')
            for index,components in spatial.items():
                i,j,b=(int(v[1:]) for v in index.split('_'))
                if set(components)!=set(COMPONENTS):raise ValueError('Cartesian component omitted')
                for component,parts in components.items():
                    expected={label for label,a,q in cartesian_templates()[component,i,j,b]}
                    if set(parts)!=expected:raise ValueError('Moving-basis contribution omitted')
                    for label,row in parts.items():check(row,label);count+=1
            for component,parts in time.items():
                expected={UR,UT} if component in ('ux','uy') else {UZ} if component=='uz' else {P}
                if set(parts)!=expected:raise ValueError('Fixed-position time contribution omitted')
                for label,row in parts.items():check(row,label);count+=1
            if count!=216 or not packet['moving_cylindrical_basis_differentiated'] or not packet['normalization_not_differentiated_twice']:
                raise ValueError('Original full physical operators/units omitted')
            counts[name]=count
        for chart in ('core','bridge_first','outer_power','heat_exterior'):
            try:field.evaluate(chart,0,0)
            except ValueError:pass
            else:raise ValueError('Unadmitted physical owner accepted')
        fixture=independent_flatten_fixed_unit_fixture();old=raw['retained_twenty_chart_physical_evidence']
        result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source,
            datum_enclosure_sha256=field.datum_sha,current_physical_chart_owners_checked=21,
            current_flatten_spatial_and_time_source_contributions_checked=counts,
            total_new_flatten_source_contributions_checked=sum(counts.values()),
            twenty_one_regular_source_contributions_checked=old['total_regular_source_contributions']+counts['whole_domain'],
            retained_twenty_chart_physical_evidence=old,retained_signed_source_terms_checked=terms,
            exact_zero_source_contributions_checked=zeros,actual_current_dispatch_call_bindings=calls,
            independent_flatten_fixed_unit_fixture=fixture,
            **{GATE:True,COMPOSED:True,UNIFORM:False},current_flatten_physical_owner_installed=True,
            current_pulse_terminal_flatten_source_ownership_certified=True,
            full_pulse_C4_installed=False,**dict.fromkeys(SCOPES,False),**dict.fromkeys(OPEN,False),all_passed=True,
            input_hashes={**raw['input_hashes'],NAME:sha(NAME),Path(__file__).name:sha(Path(__file__).name)})
    Path(__file__).with_suffix('.json').write_bytes((json.dumps(result,indent=2)+'\n').encode('utf8'))
    print('PASS current21 physical owners: flatten spatial4/time1, current source routing and 60 full-field fixed-unit rows',flush=True)
    return result


if __name__=='__main__':run()
