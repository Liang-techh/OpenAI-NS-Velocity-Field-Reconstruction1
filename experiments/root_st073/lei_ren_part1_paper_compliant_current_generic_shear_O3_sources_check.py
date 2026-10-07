"""Original O3 row, common-unit recovery and reserved-domain focused checks."""
import copy
import gzip
import json
import math
from pathlib import Path
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_shear_O3_sources as source
from lei_ren_part1_paper_compliant_current_generic_shear_source_packets_check import overlap,modal_overlap

packets=source.packets


def original_saved_checks(service,chart,packet,field):
    view,_=service.saved_view(chart);original=view[source.SPEC[chart][2]];c=service.ctx
    physical=original['physical_velocity_pressure_y_derivative_Taylor'];count=0
    labels=dict(theta='Utheta_over_Pstar',axial='Uz',radial='Ur_over_current_sqrt_R_over_2')
    for name,rows in packet.native_velocity.items():
        for j,row in enumerate(rows):
            expected=packets.jet(c,physical[labels[name]][j])
            native=row.terms.get(packets.ZERO,packets.IntervalTaylor.constant(c,0,row.order))
            for k in range(expected.order+1):
                overlap(native[k],expected[k],chart+' raw physical '+name);count+=1
    for j,row in enumerate(packet.absolute_pressure):
        expected=packets.jet(c,physical['P_over_Pstar2'][j]);native=row.terms.get(packets.ZERO,packets.IntervalTaylor.constant(c,0,row.order))
        for k in range(expected.order+1):
            overlap(native[k],expected[k],chart+' original absolute pressure');count+=1
    for name,rows in packet.velocity.items():
        for j,row in enumerate(rows):
            native=packet.algebra.shift(row,(0,.5,0,0)) if name in ('axial','radial') else row
            if packets.encode(native)!=packets.encode(packet.native_velocity[name][j]):
                raise ArithmeticError('Original velocity common units changed')
            count+=modal_overlap(row,field['physical_velocity_pressure_ordinary_y_rows'][name][j],chart+' recovered '+name)
    for name,rows in packet.histories.items():
        for j,row in enumerate(rows):
            native=packet.algebra.shift(row,(0,.5,0,0)) if name in ('m','k') else row
            if packets.encode(native)!=packets.encode(packet.native_histories[name][j]):
                raise ArithmeticError('Original history common units changed')
            count+=modal_overlap(row,field['own_normalized_history_ordinary_y_rows'][name][j],chart+' recovered '+name)
    for j,row in enumerate(packet.absolute_pressure):
        count+=modal_overlap(row,field['physical_velocity_pressure_ordinary_y_rows']['pressure'][j],chart+' recovered pressure')
    if field['original_axis_pressure_over_S_squared'] is not packet.P0:
        raise ArithmeticError('Original pressure datum replaced')
    if packet.algebra.proofs or packet.algebra.final_rows:
        raise ArithmeticError('Production source factors resolved')
    if any(v!=(0,0) for row in packet.native_velocity['axial'] for jet in row.terms.values()
           for value in jet.coefficients for v in [packets.recovery.endpoints(value)]):
        raise ArithmeticError('Original O3 axial-zero identity lost')
    return count


def independent_raw_fixture():
    """Direct differentiation of a variable log-amplitude and physical radial factor."""
    c=MPIntervalContext();c.dps=90;y,z=s.symbols('y z',real=True)
    ell=-y/s.Integer(2)-y**3/s.Integer(70)
    u=(1+z*z/3)*s.exp(ell);V=(1+z/7)*(1+y/5)
    Q=(1-z*z/4)*(1+y/9+y*y/13)
    y0=s.Rational(1,5);z0=s.Rational(2,7)
    def value(expr):return c.mpf(str(s.N(expr.subs({y:y0,z:z0}),85)))
    def jet(expr):return packets.IntervalTaylor(c,[value(s.diff(expr,z,k))/math.factorial(k) for k in range(6)])
    history={key:[jet((i+1)*(1+z/11)*s.exp(-(i+1)*y/9)).truncate(5) for j in range(5)]
        for i,key in enumerate(packets.recovery.RATES)}
    original=dict(Utheta_over_Pstar_axial5_coefficients=jet(u).coefficients,
        log_Utheta_ordinary_y_derivatives=[jet(s.diff(ell,y,j)) for j in range(1,5)],
        Uz_ordinary_y_derivative_axial5=[jet(s.diff(V,y,j)) for j in range(5)],
        actual_Q_y_derivative_axial4=[jet(s.diff(Q,y,j)) for j in range(5)],
        actual_normalized_primitive_y_derivative_axial5=history,
        original_P0_axial5_coefficients=jet(2+z/17).coefficients)
    actual,_,_,_=source.raw_rows(c,original);count=0
    expected=dict(theta=u,axial=V,radial=s.exp(y/2)*Q)
    for name,rows in actual.items():
        for j,row in enumerate(rows):
            for k in range(5):
                expr=s.diff(expected[name],y,j,z,k)
                if name=='radial':expr*=s.exp(-y/2)
                reference=value(expr);left=row[k]*math.factorial(k)
                # Independent decimal conversion error is explicitly enclosed.
                error=c.mpf('1e-75')*(1+abs(reference))
                overlap(left,reference+ c.mpf([-1,1])*error,'independent raw fixture');count+=1
    return count


def guards(service):
    view,origin=service.saved_view('O3_power')
    calls=[lambda:service.saved('unknown')]
    foreign=copy.deepcopy(view);foreign['actual_upstream_physical_spatial4_time1_packet']['datum_enclosure_sha256']='foreign'
    calls.append(lambda:service.adapt('O3_power',foreign,origin))
    bad=copy.deepcopy(view);bad['coverage_coordinate']=2
    calls.append(lambda:service.adapt('O3_power',bad,origin))
    lost=copy.deepcopy(view);del lost[source.SPEC['O3_power'][2]]['actual_normalized_primitive_y_derivative_axial5']['k']
    calls.append(lambda:service.adapt('O3_power',lost,origin))
    for call in calls:
        try:call()
        except (ValueError,KeyError,TypeError):pass
        else:raise ArithmeticError('Invalid O3 source admitted')
    return len(calls)


def run():
    manifest=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in manifest['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed current original O3 dependency: '+name)
    service=source.CurrentO3Sources();views=json.loads(gzip.decompress((source.HERE/source.VIEWS).read_bytes()))
    count=0;normcount=0
    if tuple(views)!=source.CHARTS:raise ValueError('Original O3 chart inventory differs')
    for chart in source.CHARTS:
        packet=service.saved(chart);field=packet.recover_original(service.service.data['delta'])
        if packets.encode(dict(packet=packet.record(),original_recovered_field=field))!=views[chart]:
            raise ValueError('Original source row export differs: '+chart)
        count+=original_saved_checks(service,chart,packet,field)
        actual=service.log_bounds(packet,field)
        if packets.encode(actual)!=manifest['original_O3_quotient_log_norms'][chart]:
            raise ValueError('Current O3 original signed source bound differs')
        normcount+=sum(len(v) for v in actual['admitted_original_quotient_log_norms'].values())
        if any(v['exact_zero'] is not True for key in ('b','t0') for v in
            actual['admitted_original_quotient_log_norms'][key].values()):
            raise ArithmeticError('Source axial-zero quotients lost')
    reservation=service.right_reservation();endpoints=packets.recovery.endpoints
    if packets.encode(reservation)!=manifest['right_edge_and_new_repair_reservation']:
        raise ValueError('Reserved original power geometry differs')
    if endpoints(reservation['strict_kappa_minus2_source'])[0]<=0:
        raise ArithmeticError('Canonical positive shear excess lost')
    for value in reservation['positive_log_radius_gaps'].values():
        if endpoints(value)[0]<=0:raise ArithmeticError('Reserved radius order failed')
    if any(manifest.get(k) for k in source.OPEN) or reservation['complete_Section11_loop_domain_certified']:
        raise ArithmeticError('Local O3 result exceeds scope')
    result=dict(all_passed=True,source_family=service.family,
        **{source.GATE:True,source.RESERVE_GATE:True},**dict.fromkeys(source.OPEN,False),
        original_chart_count=2,ordinary_signed_quotient_derivative_bounds=normcount,
        original_and_recovered_source_coefficients_checked=count,
        independent_variable_log_amplitude_and_radial_prefactor_derivatives=independent_raw_fixture(),
        invalid_source_family_chart_domain_history_guards=guards(service),
        original_source_factors_compared_without_materialization=True,
        separate_pressure_datum_preserved=True,strict_excess_preserved_without_subtracting_rounded_two=True,
        positive_reserved_log_radius_gaps=3,source_graph_ancestor_constructors_called=False,
        complete_Section11_loop_domain_certified=False,
        input_hashes={**manifest['input_hashes'],source.NAME:source.sha(source.NAME),
            Path(__file__).name:source.sha(Path(__file__).name),
            packets.PREFIX+'current_generic_shear_source_packets_check.py':source.sha(packets.PREFIX+'current_generic_shear_source_packets_check.py')})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Original O3 extension PASS:',count,'source coefficients,',normcount,'quotient bounds; right repair reservation PASS',flush=True)
    return result


if __name__=='__main__':run()
