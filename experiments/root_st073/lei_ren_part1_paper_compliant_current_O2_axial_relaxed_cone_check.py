"""Focused full axial source sectors, cutoff cover and degenerate branches."""
import json
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_O2_axial_relaxed_cone as source
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_stress_rows
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_outer_initial import turnoff_kernels
from lei_ren_part1_paper_compliant_pre_pulse_mixed_C4 import turnoff_derivatives


def independent_a2_product_fixtures():
    c=MPIntervalContext();c.dps=100;count=0
    with mp.workdps(130):
        for bv in ('-.7','-.1','0','.1','.7'):
            for rv in ('-20','-.5','0','.5','4','20'):
                b,r=mp.mpf(bv),mp.mpf(rv);D=1-b*r/2
                direct=2*D*D-b*b/2*(r+b/2)**2
                factor=2*(1+b*b/4)*(1-b*r-b*b/4)
                if abs(direct-factor)>mp.mpf('1e-110'):raise ArithmeticError('Independent signed a2 factorization differs')
                result=source.ref.relaxed_cone_margins(c,2,bv,1,rv,vs_minus2=c.mpf(bv)**2/2)
                expected=D>0 and (direct>0 if b else True)
                if result['admitted']!=expected:raise ArithmeticError('Independent a2 relaxed/strict branch differs')
                if b==0 and result['strict_source_shear_for_entire_box']:raise ArithmeticError('Degenerate source branch marked strict')
                count+=1
    return count


def independent_full_axial_source_fixtures():
    """Moderate synthetic source data, all original sectors nontrivial.

    These check exact source units. Whole-current-source admission comes
    from the signed identities/log budget and directed cutoff cover.
    """
    c=MPIntervalContext();c.dps=100;count=0
    with mp.workdps(130):
        Us=c.exp(-c.mpf(1)/5);J,mass=source.ref.slope_masses(c,c.mpf(1),256)
        Xs=(c.mpf(5)/8+mass[0])*c.exp(-c.mpf(3)/2)/Us
        Ds=1-Xs;AQs=-(c.mpf(5)/12+mass[2]/2)*c.exp(-1)/Us**2
        Pstar=c.exp(50)
        for tv in ('0','1','10','40'):
            t=c.mpf(tv);y=1+t;phase=c.ln(y)/40
            jets=turnoff_derivatives(c,y,c.mpf(40),phase);B,By=jets[:2]
            K=turnoff_kernels(c,y,40,cells=256,window=800)
            U=Us*c.exp(-t/2);D=Ds*c.exp(-t);M=4*c.exp(-t)+4*K['B_mass']
            # Source ODE/FTC comparisons, separate from directed quadrature
            # overestimation: enclosures must overlap these analytic ranges.
            def overlaps(value,lo,hi,label):
                vl,vh=source.endpoints(value)
                if vh<lo or vl>hi:raise ArithmeticError('Original source comparison differs: '+label)
            overlaps(B,0,1,'B_range');overlaps(M,source.endpoints(4*B)[0],4,'M_ge4B_le4')
            overlaps(K['B_squared_mass'],0,source.endpoints(1-c.exp(-t))[1],'K2_le1_minus_exp')
            if tv=='0':
                if source.endpoints(M)!=(4,4) or source.endpoints(K['B_squared_mass'])!=(0,0):
                    raise ArithmeticError('Same source initial M4/K2zero lost')
            AZ=16/Pstar**2*(c.exp(-t)+K['B_squared_mass'])/U**2
            AQ=AQs-t/2
            # A deliberately negative conditional H tests the SIGNED
            # upper-H budget algebra; it is not the current Rd datum.
            Hf=1-c.mpf('.001')*c.exp(t-20)
            for zv in ('-1','-.373','0','.51','1'):
                constant=lambda v:IntervalTaylor(c,[v]+[c.mpf(0)]*5)
                z=IntervalTaylor(c,[c.mpf(zv),c.mpf(1)]+[c.mpf(0)]*4);zero=z*0
                delta=constant(c.mpf('1e-40'));C=(1+z*z).reciprocal();L=1-delta*z*z
                A=(1-delta/2)*C+(1-delta)*z*z*C*C
                Bg=((2*delta*z*z-1)*C+2*(1-z*z)*z*z*C*C)/L
                uu=constant(U)*C;V=z*(4*B);mem=(1+z*z)*c.mpf('1e-40')
                energy=(z*z*AZ+C*C*AQ)*U**2;pressure=C*C*(-U**2*Hf/2)+mem
                rows=lambda v:[v]+[zero]*4
                ur=[uu,-uu/2]+[zero]*3;vr=[V,z*(4*By)]+[zero]*3
                raw=raw_pre_stress_rows(c,delta,z,ur,vr,dict(m=rows(z*M),h=rows(uu*(1-D)),k=rows(z*uu*(M-4*D)),e=rows(energy)),rows(pressure))
                theta=((A-C)/L+(-A/L-4*Bg)*D+(C+Bg)*M)
                actual_theta=sum(p['shape'][0] for name,p in raw['theta'].items() if name!='variable_radial_shear')/U
                G=((2*delta*z*z-2*(1-z*z))*AZ+(2*delta*C*C+4*(1-z*z)*C**3)*AQ-((1+delta)*C*C+2*(1-z*z)*C**3)*Hf)/L
                from lei_ren_part1_paper_compliant_collar_stress_C3 import axial_derivative
                expected_axial=z*G*U**2+(2*(1+delta)*z*mem-(1-z*z)*axial_derivative(mem))/L
                actual_axial=raw['axial']['retained_full_energy']['shape'][0]+raw['axial']['actual_absolute_pressure']['shape'][0]
                for label,actual,expected in (('theta',actual_theta[0],theta[0]),('energy_pressure',actual_axial[0],expected_axial[0]),
                    ('local',raw['axial']['local_axial_transport']['shape'][0][0],(-V/L)[0]),
                    ('meridional',raw['axial']['nonlinear_meridional_transport']['shape'][0][0],(V*M)[0]),
                    ('radial_shear',raw['axial']['axial_radial_shear']['shape'][0][0],(z*(8*By))[0])):
                    al,ah=source.endpoints(actual);el,eh=source.endpoints(expected)
                    if ah<el or eh<al:raise ArithmeticError('Independent full axial sector differs: '+label)
                count+=1
    return count


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Changed full original axial source: '+name)
    field=source.CurrentOriginalAxialRelaxedCone(require_checked=False);c=field.ctx
    for key,value in (('exact_full_original_axial_source_theorem',field.theorem),
        ('directed_whole_original_cutoff_cover',field.cutoff),('whole_original_axial_baseline',field.baseline),('whole_original_axial_relaxed_cone',field.proof)):
        if source.parameters.encoded(value)!=data[key]:raise ValueError('Full original axial source budget changed: '+key)
    for name,value in {**field.baseline['positive_margins'],**field.proof['positive_margins']}.items():
        if source.endpoints(value)[0]<=0:raise ArithmeticError('Full signed axial budget missing: '+name)
    rows=field.cutoff['all_cell_certificates'];previous=c.mpf(0)
    for i,row in enumerate(rows):
        lo,hi=source.endpoints(row['domain'])
        if row['cell']!=i or lo!=source.endpoints(previous)[0]:raise ValueError('Directed cutoff cover has a gap')
        if source.endpoints(c.mpf(9)/4-row['weighted_derivative_upper'])[0]<=0 or source.endpoints(c.mpf('1e-4')-row['derivative_over_actual_y_upper'])[0]<=0:
            raise ArithmeticError('Source cutoff cell does not satisfy complete budget')
        previous=c.mpf(hi)
    if source.endpoints(previous)!=(1,1):raise ValueError('Directed cover omits cutoff terminal edge')
    specs=(('whole_axial',(0,1),(-1,1)),('slope_axial_seam',0,(-1,1)),('axial_buffer_seam',1,(-1,1)),
        ('midplane',('.2','.8'),0),('positive_strict_interior',('.2','.8'),('.1','1')),
        ('negative_strict_interior',('.2','.8'),('-1','-.1')))
    for name,q,z in specs:
        query=field.query(q,z)
        if source.parameters.encoded(query)!=data['examples'][name]:raise ValueError('Original axial source query changed: '+name)
        expected='strict_interior' in name
        if query['complete_original_strict_cone_certified_for_entire_query_box']!=expected:
            raise ValueError('Axial midplane or flat edge admitted strictly')
    invalid=((-1,0),(2,0),(('.99','1.01'),0),('.5',2),('.5',('-1.01','0')),(mp.inf,0),('.5',mp.inf))
    for q,z in invalid:
        try:field.query(q,z)
        except (ValueError,TypeError):pass
        else:raise ArithmeticError('Invalid axial source box admitted')
    for key in source.OPEN:
        if data.get(key) is not False:raise ValueError('Original axial relaxed cone promoted global gate')
    comparison=field.theorem['exact_source_comparisons']
    if comparison['B_range']!=(0,1) or comparison['M_initial']!=4 or comparison['squared_kernel_initial']!=0:
        raise ValueError('Original kernel/cutoff comparison source conditions missing')
    if not comparison['pressure_budget_uses_signed_H_le1']:
        raise ValueError('Signed pressure budget replaced by an absolute H bound')
    # Both Rd-to-axial transport weights are nonnegative everywhere:
    # t-T<=-11. The accepted Rd H range remains source-bound in baseline.
    if source.endpoints(field.baseline['checked_original_Rd_future_H_bounds']['lower'])[0]<=0:
        raise ValueError('Current Rd future datum not positive')
    algebra=independent_a2_product_fixtures();full=independent_full_axial_source_fixtures()
    result=dict(all_passed=True,source_family=data['source_family'],
        exact_full_axial_source_and_cone_identities=len(field.theorem['identities']),
        strict_full_axial_baseline_and_cone_budget_inequalities=len(field.baseline['positive_margins'])+len(field.proof['positive_margins']),
        directed_contiguous_cutoff_cells_checked=len(rows),independent_a2_signed_product_branch_fixtures=algebra,
        independent_full_axial_source_sector_fixtures=full,source_queries_recomputed=len(specs),invalid_queries_rejected=len(invalid),
        actual_M_ge4B_correlated_source_budget_retained=True,
        complete_original_pressure_energy_kinetic_local_and_radial_stress_retained=True,
        current_original_O2_axial_relaxed_input_cone_certified=True,whole_closed_axial_strict_cone_certified=False,
        **{key:False for key in source.OPEN},input_hashes={**data['input_hashes'],source.NAME:source.sha(source.NAME),Path(__file__).name:source.sha(Path(__file__).name)})
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(source.parameters.encoded(result),indent=2)+'\n').encode())
    print('Original full axial relaxed cone PASS:',len(rows),'directed cells;',full,'full source cases',flush=True)
    return result


if __name__=='__main__':run()
