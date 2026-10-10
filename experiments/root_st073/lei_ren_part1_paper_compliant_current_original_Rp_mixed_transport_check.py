"""General cumulative product rules and actual fifteen-chart mixed checks.

Only new source/scale transport is checked. Reuses accepted original phase,
inlet, selection, repair and closure; no upstream solves or legacy replay.
"""
import gzip
import json
import math
from pathlib import Path
import time
import sympy as s
from types import SimpleNamespace

import lei_ren_part1_paper_compliant_current_original_Rp_mixed_transport as current
from lei_ren_part1_paper_compliant_power_angular_C4 import quotient_log_rates
from lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery import history_densities
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision


def general_equations():
    y,v,z=s.symbols('y v Z',real=True);R0=s.Symbol('R0',positive=True)
    E=s.Function('true_Utheta');V=s.Function('true_Uz')
    seed=s.symbols('incoming_z incoming_theta incoming_thetaz incoming_ztheta incoming_p')
    dens=dict(Mz=V(v,z),Mtheta=s.sqrt(2)*E(v,z),
        Mtheta_z=s.sqrt(2)*E(v,z)*V(v,z),Mztheta=V(v,z)**2-E(v,z)**2/2,Mp=E(v,z)**2/2)
    radial=dict(Mz=1,Mtheta=s.Rational(3,2),Mtheta_z=s.Rational(3,2),Mztheta=1,Mp=0)
    checked=0
    for i,(name,density) in enumerate(dens.items()):
        r=radial[name];f=density.subs(v,y)
        primal=seed[i]+s.Integral((R0*s.exp(v))**r*density,(v,0,y))
        for k in range(1,5):
            wanted=(R0*s.exp(y))**r*sum(s.binomial(k-1,j)*r**(k-1-j)*s.diff(f,y,j) for j in range(k))
            got=s.diff(primal,y,k)
            for n in range(5-k):
                assert s.simplify(s.diff(got-wanted,z,n))==0,(name,k,n);checked+=1
    P0=s.Function('independent_analytic_P0')(z)
    Mp=seed[-1]+s.Integral(E(v,z)**2/2,(v,0,y))
    for k in range(1,5):
        for n in range(5-k):
            assert s.simplify(s.diff(Mp+P0,y,k,z,n)-s.diff(E(y,z)**2/2,y,k-1,z,n))==0
            checked+=1
    return dict(passed=True,independent_full_cumulative_mixed_equations=checked,
        nonzero_axial_and_signed_quadratic_terms_retained=True,
        arbitrary_incoming_histories_and_independent_P0_retained=True,
        radial_powers_differentiated_once_and_velocity_scale_derivatives_not_duplicated=True)


def overlap(a,b,label):
    al,ah=current.pulse.radius.post.selected.inlet.endpoints(a)
    bl,bh=current.pulse.radius.post.selected.inlet.endpoints(b)
    assert max(al,bl)<=min(ah,bh),(label,'disjoint directed source rows')
    return 1


def native_affine_maps(owner):
    radius=owner.owner.radius
    field=SimpleNamespace(graph=owner.graph,parameters=radius.parameters,
        contracts=radius.frame.bridge.leading.contracts)
    v=s.Symbol('actual_native_coordinate',real=True)
    q=current.pulse.radius.RadiusInterpreter(field,False,{radius.native.node:v,
        radius.logRp.node:s.Symbol('same_absolute_logRp',real=True),
        radius.functions['waiting'].node:s.Symbol('true_Z_independent_waiting',positive=True)})
    for chart,row in radius.maps.items():
        y=q.at(row['logR']);J=q.at(row['native_to_log_radius_jacobian'])
        assert s.simplify(s.diff(y,v)-J)==0,(chart,'true native Jacobian')
        assert s.diff(y,v,2)==0 and s.diff(J,q.z)==0,(chart,'affine and Z independent')
    return dict(passed=True,actual_affine_Z_independent_current_maps=len(radius.maps),
        true_native_Jacobian_function_checked=True)


def pulse_primitive_diagnostics(owner,view):
    source=view['original_forward_source_packet'];u=source['inlet_Utheta_over_Pstar_Taylor'];c=owner.ctx
    mu=owner.owner.selected.pulse.mu
    primitive=source['primitive_y_derivative_Taylor']
    specs=dict(Mz=('Mz_over_R_Utheta',c.mpf('1/2')-mu,u),
        Mtheta=('Mtheta_over_sqrt2_R_3half_Utheta',1-mu,u*c.sqrt(2)),
        Mtheta_z=('Mtheta_z_over_sqrt2_R_3half_Utheta_squared',c.mpf('1/2')-2*mu,u*u*c.sqrt(2)),
        Mztheta=('Mztheta_over_R_Utheta_squared',-2*mu,u*u))
    count=0
    for name,(key,rate,prefactor) in specs.items():
        for k in range(1,5):
            # The original provider derives normalized primitives from its
            # own source ODEs; now differentiate their remaining raw scale.
            expected=prefactor*sum((primitive[key][j]*(math.comb(k,j)*rate**(k-j)) for j in range(k+1)),primitive[key][0]*0)
            for n in range(5-k):
                count+=overlap(view['log_radius_mixed_rows'][name]['y'+str(k)+'_Z'+str(n)].coefficients[0],
                    expected[n]*math.factorial(n),(view['chart'],name,k,n,'original primitive product'))
    return count


def heat_lograte_diagnostics(owner,view,coordinate):
    c=owner.ctx;close=owner.owner.closed;heat=close.angular.heat;Z=view['Z']
    exact=current.pulse.radius.exact_coordinate(coordinate);t=c.mpf(exact.numerator)/exact.denominator
    shape=(current.pulse.closed.closure.angular.shape_radial5(heat,Z,t) if view['chart']=='heat_collar'
        else heat.local_Gamma(Z,t))
    rates=quotient_log_rates(shape['K_rows'][:5])
    # This independent log-rate ODE computes theta derivatives rather than
    # repeating the new adapter's direct K*exp(-bh*t) product formula.
    rates[0]-=heat.bh
    rows=[view['original_factorized_values']['Utheta'].coefficients];count=0
    for k in range(4):
        value=rows[0]*0
        for j in range(k+1):value+=rates[j]*rows[k-j]*math.comb(k,j)
        rows.append(value)
    for k,row in enumerate(rows):
        for n in range(5-k):
            count+=overlap(view['log_radius_mixed_rows']['Utheta']['y'+str(k)+'_Z'+str(n)].coefficients[0],
                row[n]*math.factorial(n),(view['chart'],k,n,'independent Gamma log-rate'))
    return count


@source_precision
def run(before=None,refresh_candidate=False,resume_owner=None,resume_views=None):
    began=time.monotonic();candidate=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    changed=[]
    for name,digest in candidate['input_hashes'].items():
        if current.sha(name)!=digest:changed.append(name)
    # A source-only adapter guard change can refresh the artifact in the
    # same fifteen actual calls used for checks. All defining dependencies
    # remain fixed, and every output must equal the recorded source values.
    if refresh_candidate:
        assert set(changed)<= {Path(current.__file__).name},changed
    else:assert not changed,changed
    owner=current.CurrentOriginalRpMixedTransport(before,require_checked=False)
    if resume_views is not None:
        # Resume a live same-source observation after a heat-only metadata
        # guard repair. Do not reconstruct arrays from a receipt or accept
        # changes to the mixed algebra, evaluator or scale factor methods.
        assert resume_owner.before is before and resume_owner.graph is owner.graph
        for method in ('histories','factor','evaluate'):
            a=getattr(type(resume_owner),method)
            b=getattr(type(owner),method)
            a=getattr(a,'__wrapped__',a).__code__;b=getattr(b,'__wrapped__',b).__code__
            assert (a.co_code,a.co_consts,a.co_names)==(b.co_code,b.co_consts,b.co_names),method
        assert len(resume_owner.call_trace)==len(resume_views)
        owner.call_trace=list(resume_owner.call_trace)
    assert not owner.acceptance_loaded and not any(candidate[k] for k in current.GATES+current.OPEN)
    assert candidate['source_family']==owner.family_record and candidate['actual_source_graph']==owner.assert_graph()
    independent=general_equations();affine=native_affine_maps(owner)
    rows=0;velocity=0;first=0;primitive=0;heat=0;parts=0;views={}
    for name,chart,coordinate in current.pulse.VIEWS:
        if resume_views is not None and name in resume_views:
            view=resume_views[name]
            assert view['chart']==chart
            packet=view['original_forward_source_packet']
            rawview={'source_packet':packet,'current_source_packet':packet}
            admitted,_=owner.velocity_rows(chart,candidate['fresh_Z'],coordinate,rawview)
            for key,values in admitted.items():
                for a,b in zip(values,view['velocity_radial_source_rows'][key]):
                    assert all(x._mpi_==y._mpi_ for x,y in zip(a.coefficients,b.coefficients))
        else:view=owner.evaluate(chart,candidate['fresh_Z'],coordinate)
        views[name]=view
        assert current.pulse.raw.packed(current.view_report(view))==candidate['actual_fifteen_chart_mixed_views'][name],name
        for key,mixed in view['log_radius_mixed_rows'].items():
            assert len(mixed)==15
            for label,value in mixed.items():
                k,n=value.derivative
                native=view['native_coordinate_mixed_rows'][key]['n'+str(k)+'_Z'+str(n)]
                assert value.powers[-1]==0 and native.powers[-1]==k
                assert value.powers[:-1]==native.powers[:-1] and value.source_units==native.source_units
                assert native.coefficients is value.coefficients,(chart,key,k,n,'native J must be exact scale')
                ldict=dict(value.log_scale_parts);ndict=dict(native.log_scale_parts)
                assert all(v.graph is owner.graph for v in tuple(ldict.values())+tuple(ndict.values()))
                assert all(ldict[t].node==ndict[t].node for t in ldict if t!='native_coordinate_Jacobian')
                J=current.pulse.radius.FunctionRef(owner.graph,view['geometry']['native_to_log_radius_jacobian'])
                wanted=owner.graph.mul(owner.graph.constant(k),owner.graph.unary('log',J)) if k else owner.graph.zero
                assert ndict['native_coordinate_Jacobian'].node==wanted.node
                parts+=len(ldict)+len(ndict);rows+=1
                if k==0:
                    expected=view['original_factorized_values'][key].coefficients[n]*math.factorial(n)
                    assert value.coefficients[0]._mpi_==expected._mpi_
                if key in ('Mp','pressure') and k:
                    expected=(0,2,2 if chart in current.pulse.PULSE else 0)
                    assert value.powers[:-1]==tuple(map(current.Fraction,expected))
                if key in view['velocity_radial_source_rows']:
                    expected=view['velocity_radial_source_rows'][key][k][n]*math.factorial(n)
                    velocity+=overlap(value.coefficients[0],expected,(chart,key,k,n,'actual velocity rows'))
        # First derivatives agree with the accepted previous density caller,
        # accounting for its directed J box without selecting a point J.
        packet=view['original_forward_source_packet']
        if chart in current.pulse.PULSE:
            E=packet['inlet_Utheta_over_Pstar_Taylor'];V=E*packet['Uz_over_Utheta']
        else:
            E=view['original_factorized_values']['Utheta'].coefficients;V=E*0
        density=history_densities(E,V)
        expected=dict(Mz=density['m'],Mtheta=density['h']*owner.ctx.sqrt(2),
            Mtheta_z=density['k']*owner.ctx.sqrt(2),Mztheta=density['e'],
            Mp=density['p'],P0=E*0,pressure=density['p'])
        J=owner.owner.radius.local_step(chart,coordinate,coordinate)['native_to_log_radius_jacobian_bound']
        for key,value in expected.items():
            for n in range(4):
                actual=view['log_radius_mixed_rows'][key]['y1_Z'+str(n)].coefficients[0]*J
                first+=overlap(actual,value[n]*math.factorial(n)*J,(chart,key,n,'original first density source'))
        if chart in current.pulse.PULSE:primitive+=pulse_primitive_diagnostics(owner,view)
        elif chart.startswith('heat_'):heat+=heat_lograte_diagnostics(owner,view,coordinate)
        print('Checked current raw mixed source:',name,flush=True)
    main=views['main'];gap=views['gap']
    assert current.pulse.radius.post.selected.inlet.endpoints(main['original_factorized_values']['Uz'].coefficients[0])[0]>0
    assert all(current.pulse.radius.post.selected.inlet.endpoints(v)==(0,0) for v in gap['original_factorized_values']['Uz'].coefficients.coefficients)
    assert any(current.pulse.radius.post.selected.inlet.endpoints(v)!=(0,0) for v in gap['original_factorized_values']['Mz'].coefficients.coefficients)
    rv=owner.before.evaluate(candidate['fresh_Z']);rvchecks=0
    for side,chart in (('pulse','pulse_end'),('flatten','flatten')):
        view=owner.evaluate(chart,candidate['fresh_Z'],0)
        target=rv[side+'_log_radius_mixed_rows']
        for key,grid in view['log_radius_mixed_rows'].items():
            for label,value in grid.items():
                coefficient=value.coefficients[0]
                if side=='pulse':coefficient/=owner.before.U0_bound**int(value.powers[1])
                rvchecks+=overlap(coefficient,target[key][label].coefficients[0],(side,key,label,'accepted Rv common seam'))
    rejected=[]
    for chart,coordinate in (('core',0),('pulse_end',1),('heat_collar',4)):
        try:owner.evaluate(chart,candidate['fresh_Z'],coordinate)
        except ValueError:rejected.append(chart)
        else:raise AssertionError('Illegal source chart/domain accepted')
    try:current.CurrentOriginalRpMixedTransport(owner.owner,require_checked=False)
    except ValueError:pass
    else:raise AssertionError('Unadmitted Rv/closed source bypass')
    # Unit metadata is a required function adapter contract even though
    # retained legacy labels contain Pstar. Exercise rejection without
    # mutating a live source provider or replaying quadrature.
    fake=dict(views['power']);packet=dict(fake['original_forward_source_packet']);packet['velocity_units']='U/Pstar'
    fake['current_source_packet']=packet
    try:owner.velocity_rows('outer_power',candidate['fresh_Z'],'1/2',fake)
    except ValueError:pass
    else:raise AssertionError('Wrong Pstar units silently treated as Ev0')
    if refresh_candidate:
        candidate['input_hashes']=owner.hashes
        candidate['actual_source_call_trace']=owner.call_trace[:15]
        candidate['execution_seconds']=time.monotonic()-began
        data=json.dumps(current.pulse.raw.packed(candidate),separators=(',',':'))+'\n'
        (current.HERE/current.NAME).write_bytes(gzip.compress(data.encode(),mtime=0))
    result=dict(all_passed=True,source_family=owner.family_record,**dict.fromkeys(current.GATES,True),
        independent_general_cumulative_mixed_equations=independent,
        independent_current_native_affine_Jacobian_maps=affine,
        actual_log_radius_mixed_rows_checked=rows,actual_native_rows_with_exact_Jacobian_power_checked=rows,
        actual_current_velocity_row_overlap_diagnostics=velocity,
        accepted_first_history_density_overlap_diagnostics=first,
        actual_pulse_primitive_radial_product_overlap_diagnostics=primitive,
        independent_closed_Gamma_lograte_overlap_diagnostics=heat,
        accepted_common_unit_Rv_mixed_seam_overlap_diagnostics=rvchecks,
        final_scale_parts_all_on_one_current_graph=parts,
        active_axial_terms_and_gap_history_memory_preserved=True,
        high_pressure_rows_do_not_differentiate_caps=True,
        native_Jacobian_bound_not_used_as_defining_multiplier=True,
        wrong_postpulse_Pstar_velocity_unit_metadata_rejected=True,
        all_fifteen_current_routes_and_same_closed_heat_consumed=True,
        uniform_global_physical_time_and_point_value_admission_remain_open=True,
        rejected_source_domains=rejected,**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.pulse.raw.packed(result),indent=2)+'\n',encoding='utf8',newline='\n')
    print('PASS_CURRENT_ORIGINAL_RP_MIXED_TRANSPORT fifteen actual charts, five histories and exact native mixed scales',flush=True)
    return result


if __name__=='__main__':run()
