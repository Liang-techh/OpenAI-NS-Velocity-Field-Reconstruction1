"""Independent signed Duhamel kernels and actual whole-source envelope scope."""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_five_defect_bounds as source


def independent_kernels():
    """Signed oscillatory source, Z jets, ODE derivatives and quiet memory."""
    c=mp.mp.clone();c.dps=80;N=7;omega=2*c.pi*N
    y=c.mpf('.37');Z=c.mpf('-.23');comparisons=0;pressure_defect=None
    for j,(key,rate0) in enumerate(source.RATES.items()):
        r=c.mpf(rate0);A=c.mpf((-1)**j*(j+1));B=c.mpf(j-2)*c.mpf('.3');C=-c.mpf(j+1)*c.mpf('.2')
        poly=lambda z:z*z+c.mpf('.2')*z+c.mpf('.1')
        f=lambda t,z:(A*c.cos(omega*t)+B*t+C*poly(z))/N
        fZ=lambda t,z:C*(2*z+c.mpf('.2'))/N
        fy=lambda t,z:(-A*omega*c.sin(omega*t)+B)/N
        def exact(t,z):
            if rate0==0:
                Ic=c.sin(omega*t)/omega;Iy=t*t/2;I1=t
            else:
                Ic=(r*c.cos(omega*t)+omega*c.sin(omega*t)-r*c.exp(-r*t))/(r*r+omega*omega)
                I1=-c.expm1(-r*t)/r;Iy=t/r-I1/r
            return (A*Ic+B*Iy+C*poly(z)*I1)/N
        cutpoints=[c.mpf(0),*[c.mpf(k)/N for k in range(1,N) if c.mpf(k)/N<y],y]
        D=c.quad(lambda t:c.exp(-r*(y-t))*f(t,Z),cutpoints)
        DZ=c.quad(lambda t:c.exp(-r*(y-t))*fZ(t,Z),cutpoints)
        Dy=c.diff(lambda t:exact(t,Z),y);DyZ=c.diff(lambda z:c.diff(lambda t:exact(t,z),y),Z)
        Dyy=c.diff(lambda t:exact(t,Z),y,2)
        pairs=[(D,exact(y,Z)),(DZ,c.diff(lambda z:exact(y,z),Z)),
            (Dy,f(y,Z)-r*D),(DyZ,fZ(y,Z)-r*DZ),(Dyy,fy(y,Z)-r*Dy)]
        for actual,expected in pairs:
            if abs(actual-expected)>c.mpf('1e-65')*(1+abs(expected)):
                raise ArithmeticError('Independent signed source kernel/ODE mismatch: '+key)
            comparisons+=1
        # Triangle bound is independent of the producer: every source value
        # on[0,y] is covered; the positive kernel cannot enlarge its mass.
        cap=(abs(A)+abs(B)*y+abs(C)*abs(poly(Z)))/N
        mass=y if rate0==0 else 1/r
        if abs(D)>cap*mass:raise ArithmeticError('Independent1/r/length defect bound failed')
        comparisons+=1
        continued=lambda w:c.exp(-r*w)*D
        if abs(c.diff(continued,c.mpf(1))+r*continued(1))>c.mpf('1e-70'):
            raise ArithmeticError('Independent signed quiet transport mismatch')
        comparisons+=1
        if rate0==0:
            pressure_defect=D
            if D==0 or continued(1)!=D or c.diff(continued,c.mpf(1))!=0:
                raise ArithmeticError('Nonzero pressure memory was reset')
            comparisons+=1
    return dict(comparisons=comparisons,source='signed oscillatory plus affine-y/quadratic-Z; common N7',
        all_five_source_signs_and_nonzero_pressure_memory_exercised=True,
        independent_quadrature_vs_closed_primitives_allowance='1e-65 relative/absolute, modest kernels only',
        actual_source_scales_or_point_cover_values_materialized=False,
        nonzero_fixture_pressure_defect=c.nstr(pressure_defect,24))


def independent_union():
    c=MPIntervalContext();c.dps=80;U=source.LogUpper
    poly=lambda terms:source.packets.encode(source.source.FrequencyLogBound(c,{p:U.constant(c,v) for p,v in terms.items()}).record())
    row,witness=source.union_poly(c,{'left':poly({-1:2,-2:7}),'right':poly({-1:5,-2:3}),'quiet':poly({})})
    ep=source.packets.recovery.endpoints;comparisons=0
    if witness!={'-2':'left','-1':'right'}:raise ArithmeticError('Cover union maxima bound a sum rather than a union')
    for N in (1,7,83):
        cap=ep(c.exp(row.evaluate(c.ln(N)).log))[1]
        if cap<max(2/N+7/(N*N),5/N+3/(N*N)):
            raise ArithmeticError('Independent frequency union bound failed')
        comparisons+=1
    return dict(comparisons=comparisons,per_power_maximum_witnesses=witness,
        no_overlapping_cover_widths_added=True,fixture_only_bound_exponentials_evaluated=True)


def run():
    manifest=json.loads((source.HERE/source.NAME).read_bytes());provider=source.CurrentFiveDefectBounds()
    if manifest!=source.packets.encode(provider.compute()):raise ValueError('Actual whole-source envelopes differ from their defining source attachment')
    if not manifest[source.GATE] or manifest[source.FUNCTION_GATE] or any(manifest[k] for k in source.OPEN):
        raise ArithmeticError('Envelope promoted a source-function/seam/global stage')
    if manifest['actual_changed_five_moment_transport_integrated'] or manifest['actual_own_histories_radial_pressure_stress_point_evaluated']:
        raise ArithmeticError('No actual point integration may be admitted by a bound receipt')
    terms=0
    for key,orders in manifest['actual_whole_defect_frequency_majorants'].items():
        if set(orders)!=set(source.ORDERS):raise ArithmeticError('Cumulative defect derivative scope differs')
        for order in ('value','Z','y','yZ'):
            if any(t['N_power']>=0 for t in orders[order]['terms']):raise ArithmeticError('Cumulative ODE/source bound lost negative N powers')
            terms+=len(orders[order]['terms'])
        if not any(t['N_power']==0 for t in orders['yy']['terms']):raise ArithmeticError('Fast second-y source term omitted')
    pressure=manifest['actual_new_repair_inlet_Rc_pre_repair_defect_majorants']['p']
    oldpressure=manifest['actual_whole_defect_frequency_majorants']['p']
    if pressure['value']!=oldpressure['value'] or pressure['Z']!=oldpressure['Z'] or not pressure['pressure_memory_exactly_preserved']:
        raise ArithmeticError('Actual pressure quiet bound lost incoming memory')
    if not all(pressure[k]['exact_zero'] for k in ('y','yZ','yy')):raise ArithmeticError('Quiet pressure has spurious y derivatives')
    for key,contract in manifest['source_linked_defect_integral_contracts'].items():
        if len(contract['source_cover_references'])!=17 or contract['actual_source_owner_and_seam_identity_admitted']:
            raise ArithmeticError('Integral contract interpreted covers as source functions')
    kernels=independent_kernels();union=independent_union()
    hashes=dict(manifest['input_hashes']);hashes[source.NAME]=source.sha(source.NAME)
    hashes[Path(__file__).name]=source.sha(Path(__file__).name)
    result=dict(all_passed=True,source_family=provider.family,**{source.GATE:True,source.FUNCTION_GATE:False},
        **dict.fromkeys(source.OPEN,False),actual_whole_source_charts=17,
        cumulative_negative_N_terms_for_value_Z_y_yZ=terms,
        exact_Duhamel_rate_quiet_identities=len(provider.theorem['exact_identities']),
        independent_signed_kernel= kernels,independent_cover_union=union,
        overlapping_source_cover_widths_never_added=True,
        actual_log_radius_span_upper_bound_preserves_original_formal_left_offset=True,
        positive_rate_kernel_mass_vs_zero_rate_pressure_length_separated=True,
        actual_new_repair_inlet_Rc_before_repair_bound=True,
        pressure_quiet_memory_and_P0_preserved=True,
        derivative_scope_is_chart_interiors_and_one_sided_traces=True,
        source_function_replay_or_seam_identity_admitted=False,
        actual_changed_five_moment_transport_integrated=False,current_whole_N_selected=False,
        sufficient_mixed4_velocity_or_mixed3_stress_admitted=False,
        source_graph_ancestor_constructors_called=False,input_hashes=hashes)
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Whole five-defect Duhamel log bounds PASS:17 source cover union, signed kernels and quiet pressure memory',flush=True)
    return result


if __name__=='__main__':run()
