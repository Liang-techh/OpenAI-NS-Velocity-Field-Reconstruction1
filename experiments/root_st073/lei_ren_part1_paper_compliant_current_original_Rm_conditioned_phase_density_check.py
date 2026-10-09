"""Independent scalar Z/phase derivatives, radius identities and live Rm cells."""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_Rm_conditioned_phase_density as current
import lei_ren_part1_paper_compliant_current_generic_shear_loop as original
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow

fields,ep=current.fields,current.ep


def scalar_fixtures():
    c=MPIntervalContext();c.dps=200;z=sy.Symbol('Z',real=True);R=sy.Rational
    scales=original.GenericLoopScales(a_min='.7',margin_min='1',boundary_kappa_excess_min='.019',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='2',dps=140)
    p=scales.ctx;eta=sy.Rational(p.nstr(scales.eta,145));dstar=c.ln(c.mpf(p.nstr(scales.d_star,145)))
    f=MacroFlow(c,c.mpf(-30),c.ln(9),c.ln(4),c.ln(5),c.ln(100)-c.ln(5),c.mpf('.01'))
    allowance=c.mpf(['-1e-27','1e-27']);tol=mp.mpf('1e-130');phi=p.mpf('.137');h=p.mpf('1e-8')
    evidence=[]
    for label,bb,pp in (('positive_Mobius',R(1,5),R(1,5)),('negative_Mobius',-R(1,5),-R(1,5)),
        ('small_r',R(1,5),R(1,50)),('zero_crossing',R(1,5),sy.Integer(0))):
        a=R(4,5)+R(7,100)*z;b=bb-R(1,20)*z;p2=pp+R(7,100)*z;E=R(13,10)+z/25
        Delta=a+b*b/a-2;q=sy.sqrt((2*eta-Delta)/(2*a))
        def pair(expr):
            rows={}
            for n in (0,1):
                value=sy.lambdify(z,sy.diff(expr,z,n),'mpmath')(mp.mpf(0))
                rows[(0,n)]=f.scalar(c.mpf([value-tol,value+tol]))
            return rows
        roots={name:pair(expr) for name,expr in dict(a=a,t0=-b/a,p2=p2,E=E).items()};qr=pair(q)
        source=dict(roots=roots,q=qr[(0,0)])
        got=current.Z_FIRST(source,qr,dstar,c.mpf('.137'))
        assert got['record']['status']=='enclosed' and set(got['values'])==set(current.OUTPUTS)
        af,bf,pf,ef=[sy.lambdify(z,expr,'mpmath') for expr in (a,b,p2,E)]
        def reference(zz,target,guess=None):
            loop=original.GenericShearLoop(scales,a=af(zz),b=bf(zz),p1=5,p2=pf(zz),Utheta=ef(zz))
            if target in (0,p.mpf('.5'),1):angle=2*p.pi*target
            else:
                start=2*p.pi*target if guess is None else guess
                try:angle=p.findroot(lambda v:loop.phase_at_angle(v)-target,(start,start+p.mpf('.001')),tol=p.mpf('1e-125'))
                except (ValueError,ArithmeticError):angle=loop.angle_at_phase(target)
            return loop.at_angle(angle),angle
        ref,guess=reference(p.mpf(0),phi);comparisons=0
        for key,rkey in (('A','A'),('B_over_Pstar','B')):
            lo,hi=ep(current.first.phase.bounded_value(got['values'][key])+allowance)
            assert lo<=ref[rkey]<=hi,(label,key);comparisons+=1
        for direction in ('Z','phi'):
            samples=[reference(shift*h if direction=='Z' else p.mpf(0),phi+(shift*h if direction=='phi' else 0),guess)[0] for shift in (-2,-1,1,2)]
            for key,rkey in (('A','A'),('B','B')):
                value=(samples[0][rkey]-8*samples[1][rkey]+8*samples[2][rkey]-samples[3][rkey])/(12*h)
                name=key+'_'+direction+('_over_Pstar' if key=='B' else '')
                lo,hi=ep(current.first.phase.bounded_value(got['values'][name])+allowance)
                assert lo<=value<=hi,(label,name);comparisons+=1
        N=257;native=current.densities.density_Z_kernels(roots['E'][(0,0)],roots['E'][(0,1)],f.scalar('.4'),f.scalar('.03'),got['values'],N)
        def densities_at(zz):
            primitive,_=reference(zz,phi,guess);ee=ef(zz);vv=p.mpf('.4')+p.mpf('.03')*zz
            dE=ee*p.expm1(primitive['A']/N);dV=primitive['B']/N;cross=ee*dE+dE*dE/2
            return dict(m=dV,h=dE,k=vv*dE+ee*dV+dE*dV,e=2*vv*dV+dV*dV-cross,p=cross)
        c0=densities_at(p.mpf(0));samples=[densities_at(shift*h) for shift in (-2,-1,1,2)]
        for key in ('m','h','k','e','p'):
            value=(samples[0][key]-8*samples[1][key]+8*samples[2][key]-samples[3][key])/(12*h)
            for row,want in ((native['kernels'][key],c0[key]),(native['Z_derivatives'][key],value)):
                lo,hi=ep(current.first.phase.bounded_value(row)+allowance);assert lo<=want<=hi,(label,key);comparisons+=1
        assert not any('_y' in key for key in got['values'])
        if label=='zero_crossing':assert not got['values']['A_Z'].zero
        evidence.append(dict(passed=True,case=label,geometry=got['record']['geometry'],
            independent_original_primitive_Z_phase_and_five_density_comparisons=comparisons,
            finite_difference='five point, h=1e-8, scalar precision140, allowance1e-27; diagnostic only'))
    return evidence


def radius_identities():
    P,C,B,sc,x=sy.symbols('P C hb sc x',positive=True)
    logR=sy.log(110)+10*(C+P)-6+sy.log(x);logminus=sy.log(4)-4*P-1000+B*sc/2
    expected=sy.log(sy.Rational(110,4))+14*P+10*C+994+sy.log(x)-B*sc/2
    assert sy.simplify(logR-logminus-expected)==0
    assert sy.simplify(sy.diff(expected,x)-1/x)==0
    rejects=0
    for N in (True,159,160.0,0,-1):
        try:current.candidate_N(N)
        except ValueError:rejects+=1
        else:raise AssertionError('Nonadmitted candidate N accepted')
    assert current.candidate_N(160)==160
    c=MPIntervalContext();c.dps=100
    seam=current.radius_phase.periodic_add(c,[c.mpf(['.99','1'])],c.mpf(['0','.02']))
    assert len(seam['boxes'])==2 and not seam['full_period']
    whole=current.radius_phase.periodic_add(c,[c.mpf('.3')],c.mpf([0,2]))
    assert whole['full_period'] and ep(whole['boxes'][0])==(0,1)
    return dict(passed=True,exact_original_log_radius_minus_inlet_identity=True,original_dx_over_x_Jacobian=True,
        invalid_N_rejected=rejects,periodic_seam_union_and_full_period_cover_checked=True)


def native(owner,report):
    boxes=primitive_rows=density_rows=source_checks=full_period=closed_cells=0
    for label in ('0','.5'):
        for name,call in (('active',lambda:owner.spatial_point(label,(5,4),257)),
            ('terminal',lambda:owner.spatial_point(label,(2,1),257)),('Rh',lambda:owner.spatial_point(label,'Rh',257)),
            ('active_cell',lambda:owner.spatial_cell(label,(124999999,100000000),(125000001,100000000),257)),
            ('terminal_cell',lambda:owner.spatial_cell(label,(71,40),'Rh',257))):
            data=call();assert current.base.encoded(fields.serialized(data))==report['frames'][label][name]
            op=owner.upstream.upstream.upstream.owner(label).op;f=op.flow;geometry=data['actual_original_Rm_radius_phase']
            assert geometry['exact_source_Rm_factor'] is op.Rm_factor
            assert geometry['exact_source_logRa'] is f.logs[3]
            assert geometry['actual_Rm_radius_phase_Z_independent'] and geometry['actual_positive_width_not_zeroed']
            assert all(geometry['actual_parameter_source_binding'].values());source_checks+=1
            if name.endswith('_cell'):
                assert data['actual_closed_radial_source_cell'] and not geometry['source_geometry']['point'];closed_cells+=1
                assert geometry['source_radial_measure']=='dy=dx/x' and ep(geometry['source_radial_logarithmic_cell_width'])[0]>0
                if name=='terminal_cell':
                    assert geometry['whole_period_cover'] and ep(geometry['phase_boxes'][0])==(0,1);full_period+=1
            assert data['actual_original_spatial_Z_density_interface_installed']
            assert not data['actual_spatial_y_derivative_or_radial_integral_installed']
            for cell in data['actual_source_bound_phase_density_cells']:
                assert cell['source_geometry']==geometry['source_geometry'] and cell['original_common_P0_axial5'] is op.P0
                assert cell['original_phase_Z_only_result']['status']=='enclosed'
                values=cell['actual_original_primitive_Z_values'];density=cell['actual_original_five_signed_density_C0_Z'];boxes+=1
                assert set(values)==set(current.OUTPUTS) and not any('_y' in key for key in values)
                rows=[*values.values(),*density['kernels'].values(),*density['Z_derivatives'].values(),*density['velocities'].values()]
                assert all(row.ctx is f.c and row.scale.bases is f.logs and row.ledger is f.ledger for row in rows)
                assert set(density['kernels'])==set(density['Z_derivatives'])=={'m','h','k','e','p'}
                primitive_rows+=len(values);density_rows+=len(density['kernels'])+len(density['Z_derivatives'])
                assert not cell['phase_is_independent_candidate_parameter']
                assert cell['source_q_Z_includes_all_b_and_b_Z_terms'] and cell['no_fabricated_y_derivative_exported']
                assert not cell['full_stress_cone_or_global_N_or_density_integrals_admitted']
                assert all(cell[key] is False for key in fields.previous.OPEN)
            assert all(data[key] is False for key in fields.previous.OPEN)
    return dict(passed=True,actual_source_bound_phase_density_boxes=boxes,actual_original_primitive_Z_phi_rows=primitive_rows,
        actual_five_signed_C0_Z_density_rows=density_rows,same_actual_radius_parameter_and_P0_binding_checks=source_checks,
        actual_closed_radial_source_cells=closed_cells,exact_terminal_whole_period_phase_covers=full_period,
        conditional_two_axial_frames_only=True,global_N_cone_integrals_and_recursion_not_admitted=True)


def run():
    began=time.monotonic()
    with mp.workdps(180):fixtures=scalar_fixtures();identities=radius_identities()
    print('Independent original phase/Z primitives and five densities PASS',flush=True)
    report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    owner=current.OriginalRmConditionedPhaseDensity(require_checked=False);evidence=native(owner,report)
    assert report[current.GATE] and report['source_family']==owner.family
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        independent_scalar_phase_Z_density_fixtures=fixtures,exact_radius_and_candidate_N_guards=identities,
        actual_live_source=evidence,original_Z_only_compiler_binding=owner.bindings,
        original_parameter_source_binding=owner.parameter_bindings,
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual Rm source-bound conditioned phase and five signed Z densities PASS',flush=True)
    return result


if __name__=='__main__':run()
