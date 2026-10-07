"""Independent Duhamel, no-reset, narrow-cell and cached-source checks."""
import json
import gzip
import math
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_generic_shear_moment_recovery as source


FAMILY=dict(actual_five_defect_family_sha256='fixture',implicit_source_sha256='fixture_source',datum_enclosure_sha256='fixture_P0')


def contains(interval,value,label):
    lo,hi=source.endpoints(interval)
    if not lo<=value<=hi:
        raise ArithmeticError(label+': independently integrated value outside source cover')


def transport_fixtures():
    c=MPIntervalContext();c.dps=130
    constant=lambda value:IntervalTaylor.constant(c,value,5)
    # This annular fixture is not a current source value or an axis datum.
    inlet={key:constant(value) for key,value in dict(m='.31',h='.52',k='.23',e='-.17',p='.41').items()}
    P0=constant('-1.9');state=source.GenericMomentRecovery(c,source_family=FAMILY,P0=P0,original=inlet)
    z=IntervalTaylor.variable(c,'.2',5)
    E=1+z*z/3;V=constant('.7')+z/4
    dE=constant('.03')+z/50;dV=constant('-.02')+z*z/40
    widths=('.3','1e-140');counts=dict(analytic_transport_coefficients=0,quiet_no_reset_coefficients=0,
                                    thin_width_nonzero_injection_cases=0,split_transport_coefficients=0,
                                    independent_physical_y_derivative_values=0)
    with mp.workdps(210):
        def scalar_jet(value):return [mp.make_mpf(v._mpi_[0]) for v in value.coefficients]
        # Independent scalar convolution, avoiding the interval operator.
        Ej,Vj,ej,uj=(scalar_jet(v) for v in (E,V,dE,dV))
        def product(a,b):return [sum(a[j]*b[n-j] for j in range(n+1)) for n in range(6)]
        def add(*rows):return [sum(row[n] for row in rows) for n in range(6)]
        def times(a,x):return [v*x for v in a]
        old=dict(m=Vj,h=Ej,k=product(Ej,Vj),e=add(product(Vj,Vj),times(product(Ej,Ej),-mp.mpf('.5'))),p=times(product(Ej,Ej),mp.mpf('.5')))
        inc=dict(m=uj,h=ej,k=add(product(Vj,ej),product(Ej,uj),product(ej,uj)),
                 e=add(times(product(Vj,uj),2),product(uj,uj),times(product(Ej,ej),-1),times(product(ej,ej),-mp.mpf('.5'))),
                 p=add(product(Ej,ej),times(product(ej,ej),mp.mpf('.5'))))
        for width in widths:
            next_state=state.advance(width=width,E=E,V=V,delta_E=dE,delta_V=dV,source_family=FAMILY)
            w=mp.mpf(width)
            for key,rate in source.RATES.items():
                r=mp.mpf(rate);decay=mp.exp(-r*w);mass=w if r==0 else -mp.expm1(-r*w)/r
                for j in range(6):
                    base=scalar_jet(inlet[key])[j]*decay+old[key][j]*mass
                    defect=inc[key][j]*mass
                    contains(next_state.original[key][j],base,'original '+key)
                    contains(next_state.defect[key][j],defect,'changed '+key)
                    contains(next_state.own()[key][j],base+defect,'own '+key)
                    counts['analytic_transport_coefficients']+=1
            if source.endpoints(next_state.defect['p'][0])[0]<=0:
                raise ArithmeticError('Actual positive pressure-square increment erased')
            if width=='1e-140':counts['thin_width_nonzero_injection_cases']+=1
        changed=state.advance(width='.3',E=E,V=V,delta_E=dE,delta_V=dV,source_family=FAMILY)
        quiet=changed.advance(width='1.2',E=E,V=V,delta_E=E*0,delta_V=V*0,source_family=FAMILY)
        for key,rate in source.RATES.items():
            r=mp.mpf(rate);first_mass=mp.mpf('.3') if not r else -mp.expm1(-r*mp.mpf('.3'))/r
            for j in range(6):
                target=inc[key][j]*first_mass*mp.exp(-r*mp.mpf('1.2'))
                contains(quiet.defect[key][j],target,'no reset '+key)
                counts['quiet_no_reset_coefficients']+=1
            if key!='p' and source.endpoints(quiet.defect[key][0])==(0,0):
                raise ArithmeticError('Quiet gap erased a nonzero incoming defect')
        if quiet.P0 is not P0:raise ArithmeticError('Original axis pressure object replaced')
        # Multiple source cells still enclose the same continuous history.
        split=state
        for _ in range(12):split=split.advance(width='.025',E=E,V=V,delta_E=dE,delta_V=dV,source_family=FAMILY)
        for key,rate in source.RATES.items():
            r=mp.mpf(rate);w=mp.mpf('.3');mass=w if not r else -mp.expm1(-r*w)/r
            for j in range(6):
                contains(split.defect[key][j],inc[key][j]*mass,'continuous split '+key)
                counts['split_transport_coefficients']+=1
        # Exact constant-in-y profile reference with the changed incoming
        # histories retained through a quiet gap. Recover physical radial
        # derivatives independently, including the one +1/2 shift.
        truth={}
        for key,rate in source.RATES.items():
            r=mp.mpf(rate);w=mp.mpf('1.5');first=mp.mpf('.3')
            mass=w if not r else -mp.expm1(-r*w)/r
            first_mass=first if not r else -mp.expm1(-r*first)/r
            value=add(times(scalar_jet(inlet[key]),mp.exp(-r*w)),times(old[key],mass),
                      times(inc[key],first_mass*mp.exp(-r*mp.mpf('1.2'))))
            rows=[value]
            for j in range(4):
                rows.append(add(old[key] if j==0 else [mp.mpf(0)]*6,times(rows[-1],-r)))
            truth[key]=rows
        zero=E*0
        recovered=quiet.field_rows(Z=z,delta='.001',E_rows=[E,zero,zero,zero,zero],V_rows=[V,zero,zero,zero,zero])
        zp=mp.mpf('.2');de=mp.mpf('.001');L=1-de*zp*zp;d=1-zp*zp
        Q=[(2*zp*(Vj[0] if j==0 else 0)-(1-de)*zp*truth['m'][j][0]-d*truth['m'][j][1])/L for j in range(5)]
        for j in range(5):
            target=sum(mp.binomial(j,k)*mp.mpf('.5')**(j-k)*Q[k] for k in range(j+1))
            contains(recovered['physical_velocity_pressure_ordinary_y_rows']['radial'][j][0],target,'physical radial derivative shift')
            ptarget=mp.mpf('-1.9')+truth['p'][0][0] if j==0 else truth['p'][j][0]
            contains(recovered['physical_velocity_pressure_ordinary_y_rows']['pressure'][j][0],ptarget,'pressure derivative convention')
            counts['independent_physical_y_derivative_values']+=2
    field=quiet.field(Z=z,delta='.001',E=E,V=V,E_y=E*0,V_y=V*0)
    if field['absolute_pressure_over_S_squared'].coefficients!=(quiet.P0+quiet.own()['p']).coefficients:
        raise ArithmeticError('Recovery lost the original absolute-pressure convention')
    if any(field[key] is not False for key in source.OPEN):raise ValueError('Recovery promoted missing global loop stages')
    return counts


def actual_cache_checks(cache):
    count=0
    for chart in ('Rh_reference','O2_slope','O2_axial','O2_buffer'):
        field=cache.recover(chart)
        pre=cache.records[chart+'_whole']['actual_upstream_original_pre_O2_source']
        physical=pre['physical_velocity_pressure_y_derivative_Taylor']
        # Compare coefficients against separately serialized original source
        # covers. Overlap alone is a range sanity check; the source program
        # identities supply the function/unit equality.
        pairs=(('theta','Utheta_over_Pstar',1),('axial','Uz',cache.invS),
               ('radial','Ur_over_current_sqrt_R_over_2',cache.invS),('pressure','P_over_Pstar2',1))
        for name,raw_name,factor in pairs:
            rows=field['physical_velocity_pressure_ordinary_y_rows'][name]
            for derivative,raw_row in enumerate(physical[raw_name]):
                actual=rows[derivative];expected=cache.jet(raw_row)*factor
                for j in range(min(actual.order,expected.order)+1):
                    al,ah=source.endpoints(actual[j]);bl,bh=source.endpoints(expected[j])
                    if ah<bl or bh<al:raise ArithmeticError('Converted current source cover differs: '+chart+' '+name)
                    count+=1
        for name,raw_rows in pre['actual_normalized_primitive_y_derivative_axial5'].items():
            factor=cache.invS if name in ('m','k') else 1
            for derivative,raw in enumerate(raw_rows):
                actual=field['own_normalized_history_ordinary_y_rows'][name][derivative];expected=cache.jet(raw)*factor
                for j in range(min(actual.order,expected.order)+1):
                    al,ah=source.endpoints(actual[j]);bl,bh=source.endpoints(expected[j])
                    if ah<bl or bh<al:raise ArithmeticError('Current history derivative units differ')
                    count+=1
        if field['source_family']!=cache.family or not field['covers_saved_original_functions_only']:
            raise ValueError('Foreign or point-relabeled current cache')
        if any(field[key] is not False for key in source.OPEN):raise ValueError('Original source promoted generic changed field')
    return count


def guards():
    c=MPIntervalContext();c.dps=90;one=IntervalTaylor.constant(c,1,5);zero=one*0
    original={key:zero for key in source.RATES}
    state=source.GenericMomentRecovery(c,source_family=FAMILY,P0=one,original=original)
    inputs=dict(width='.1',E=one,V=one,delta_E=zero,delta_V=zero,source_family=FAMILY)
    calls=[]
    for key,value in (('width',-1),('width','inf'),('source_family',{**FAMILY,'datum_enclosure_sha256':'foreign'}),
                      ('E',IntervalTaylor.constant(c,1,0))):
        kwargs={**inputs,key:value};calls.append(lambda kwargs=kwargs:state.advance(**kwargs))
    calls.append(lambda:source.GenericMomentRecovery(c,source_family=FAMILY,P0=one,original={'m':one}))
    calls.append(lambda:state.field(Z=IntervalTaylor.variable(c,2,5),delta='.01',E=one,V=one,E_y=zero,V_y=zero))
    calls.append(lambda:state.field(Z=IntervalTaylor.variable(c,0,5),delta=1,E=one,V=one,E_y=zero,V_y=zero))
    for call in calls:
        try:call()
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid own-history recovery input admitted')
    return len(calls)


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed own-history source: '+name)
    theorem=source.exact_theorem()
    if source.encode(theorem)!=data['exact_own_five_history_and_full_recovery_theorem']:
        raise ValueError('Own-history/full signed recovery theorem changed')
    cache=source.CurrentO2RecoveryCache()
    views=json.loads(gzip.decompress((source.HERE/source.VIEWS).read_bytes()))
    if data['current_original_O2_common_unit_cover_views']!=source.VIEWS or list(views)!=data['current_original_O2_saved_domains']:
        raise ValueError('Current common-unit source views changed')
    for chart in views:
        if source.encode(cache.recover(chart))!=views[chart]:
            raise ValueError('Current O2 source conversion changed: '+chart)
    result=dict(all_passed=True,**{source.GATE:True,source.CACHE_GATE:True},**dict.fromkeys(source.OPEN,False),
                source_family=cache.family,exact_source_and_recovery_identities=len(theorem['identities']),
                **transport_fixtures(),actual_current_original_cache_coefficients_checked=actual_cache_checks(cache),
                invalid_source_and_coordinate_guards=guards(),
                scope='Conditional source-cover Duhamel and exact recovery identities; current original O2 cover attachment only. No modified full-current function/moment/repair/N admission.',
                input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_text(json.dumps(result,indent=2)+'\n',encoding='utf8')
    print('Own five-history transport, quiet memory, thin width and actual O2 source recovery PASS',flush=True)
    return result


if __name__=='__main__':run()
