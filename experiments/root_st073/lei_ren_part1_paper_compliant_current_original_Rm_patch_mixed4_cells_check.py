"""Independent closed primitive derivatives, native cells, raw units and joins."""
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
import time
from types import SimpleNamespace
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_Rm_patch_mixed4_cells as current
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

fields,ep=current.fields,current.ep


def finite(row):return row.ctx.mpf(0) if row.zero else row.coefficient*row.ctx.exp(row.scale.evaluate())


def contains(row,value):
    lo,hi=ep(finite(row));tol=mp.power(10,-mp.mp.dps+30)*max(1,abs(value))
    assert lo<=value+tol and value-tol<=hi,(mp.nstr(lo,16),mp.nstr(value,16),mp.nstr(hi,16))


def symbolic_sources():
    x,z=sy.symbols('x z',positive=True);R=sy.Rational
    a=R(1,500)*(1+z+z**3/7);b=R(3,1000)*(1-z/3+z*z/5);am=sy.exp(-R(3,5))/(1+z*z)
    integ=lambda p:(x**(p+1)-1)/(p+1)
    G=4*z+R(1,500)*(1+z*z)+4*z*(x-1)+b*integ(3);m=G/x
    theta=R(5,8)+z/1000+integ(R(3,5))+a*integ(R(7,2))
    mixed=R(5,2)*z+R(1,1000)*(1+z*z)+4*z*(theta-R(5,8)-z/1000)+b*integ(R(18,5))+a*b*integ(R(13,2))
    Pstar=R(12,5);Rm=R(3);delta=R(1,2000)
    energy=-R(2,5)+z*z/100+b*b*integ(6)/(Pstar*am)**2-integ(R(1,5))/2-a*integ(R(31,10))-a*a*integ(6)/2
    pm=R(5,2)+z**3/1000+integ(-R(4,5))/2+a*integ(R(21,10))+a*a*integ(5)/2
    H=x**R(1,10)+a*x**3;g=b*x**3;V=4*z+g
    Q=(2*z*V-z*m*(1-delta)-(1-z*z)*sy.diff(m,z))/(1-delta*z*z)
    p0=R(3,100)+z/100+z**4/1000;Mp=(Pstar*am)**2*pm
    primitives=dict(Mz=Rm*G,Mtheta=sy.sqrt(2)*Rm**R(3,2)*Pstar*am*theta,
        Mtheta_z=sy.sqrt(2)*Rm**R(3,2)*Pstar*am*mixed,
        Mztheta=Rm*(Pstar*am)**2*energy+8*z*Rm*G-16*z*z*Rm*x,Mp=Mp)
    physical=dict(Utheta=Pstar*am*H,Uz=V,Ur=sy.sqrt(Rm*x/2)*Q,P=Pstar**2*p0+Mp,**primitives)
    raw=dict(m=m,h=am*theta/x**R(3,2),k=am*mixed/x**R(3,2),e=primitives['Mztheta']/(Rm*x*Pstar**2),
        p=am**2*pm,P=physical['P']/Pstar**2,theta=am*H,axial=V,radial=sy.sqrt(x)*Q)
    initial=dict(mass=m,theta=theta,mixed=mixed,energy=energy,pressure=pm)
    def dy(expr):return x*sy.diff(expr,x)
    exact=0
    for expr,rhs in ((G,V),(theta,sy.sqrt(x)*H),(mixed,sy.sqrt(x)*H*V),
        (energy,g*g/(Pstar*am)**2-H*H/2),(pm,H*H/(2*x))):
        assert sy.simplify(sy.diff(expr,x)-rhs)==0;exact+=1
    F=sy.Function('F')(x);value=F
    for k in range(1,5):
        value=dy(value);expected=sum(current.original.stirling_second(k,j)*x**j*sy.diff(F,x,j) for j in range(1,k+1))
        assert sy.simplify(value-expected)==0;exact+=1
    for k in range(5):
        value=F
        for _ in range(k):value=dy(value)+value/2
        other=sy.sqrt(x)*F
        for _ in range(k):other=dy(other)
        assert sy.simplify(value-other/sy.sqrt(x))==0;exact+=1
    de=sy.Symbol('delta',real=True)
    terminal_Q=(8*z*z-4*z*z*(1-de)-4*(1-z*z))/(1-de*z*z)
    assert sy.simplify(terminal_Q-4*((2+de)*z*z-1)/(1-de*z*z))==0;exact+=1
    targets={};raw_targets={}
    for name,expr in physical.items():
        ys=expr
        for k in range(5):
            if k:ys=dy(ys)
            for n in range(5-k):
                targets[name,'x',k,n]=sy.lambdify((x,z),sy.diff(expr,x,k,z,n),'mpmath')
                targets[name,'y',k,n]=sy.lambdify((x,z),sy.diff(ys,z,n),'mpmath')
                targets[name,'R',k,n]=sy.lambdify((x,z),sy.diff(expr,x,k,z,n)/Rm**k,'mpmath')
    for name,expr in raw.items():
        ys=expr
        for k in range(5):
            if k:ys=dy(ys)
            for n in range(5-k):
                v=sy.diff(ys,z,n)/(sy.sqrt(x) if name=='radial' else 1)
                raw_targets[name,k,n]=sy.lambdify((x,z),v,'mpmath')
    return dict(x=x,z=z,a=a,b=b,am=am,p0=p0,initial=initial,targets=targets,raw_targets=raw_targets,
        exact_closed_primitive_and_coordinate_identities=exact)


def fixture(source,xvalue,zvalue):
    c=MPIntervalContext();c.dps=180;xx=mp.mpf(xvalue);zz=mp.mpf(zvalue);tol=mp.mpf('1e-95')
    x,z=source['x'],source['z'];f=MacroFlow(c,c.ln(c.mpf('1e-5')),c.ln(c.mpf(144)/25),c.mpf(0),c.mpf(0),c.ln(100),c.mpf(1)/16)
    def jet(expr):
        return IntervalTaylor(c,[c.mpf([v-tol,v+tol])/math.factorial(n) for n in range(6)
            for v in (sy.lambdify((x,z),sy.diff(expr,z,n),'mpmath')(xx,zz),)])
    zjet=IntervalTaylor.variable(c,c.mpf(zz),5);am=jet(source['am']);zero=f.jet(zjet*0)
    controls=[f.jet(jet(source['b'])),zero,f.jet(jet(source['a'])),zero,zero]
    poly=(1+zjet*zjet)*(1+zjet*zjet)*c.exp(c.mpf('1.2'))
    op=SimpleNamespace(flow=f,c=c,controls=controls,zrows=f.jet(zjet),amrows=f.jet(am),P2=f.factor((0,1,0,0,0)),
        Pstar=f.factor((0,.5,0,0,0)),Rm_factor=f.scalar(3),P0=f.jet(jet(source['p0'])),
        invAm2=f.scale(f.jet(poly),f.factor((0,-1,0,0,0))),reference=SimpleNamespace(Z=c.mpf(zz),delta=c.mpf(1)/2000))
    initial={name:f.jet(jet(expr)) for name,expr in source['initial'].items()}
    gamma=[[c.mpf([v-tol,v+tol]) for v in (xx**3,3*xx*xx,6*xx,mp.mpf(6),mp.mpf(0))],[c.mpf(0)]*5,[c.mpf(0)]*5]
    packet=current.mixed_functions(op,c.mpf([xx-tol,xx+tol]),initial,gamma);count=rawcount=0
    for group,names in (('velocity',('Utheta','Uz','Ur')),('pressure',('P',)),('five_primitive',('Mz','Mtheta','Mtheta_z','Mztheta','Mp'))):
        for letter in ('x','y','R'):
            for name in names:
                grid=packet['physical_'+group+'_'+letter+'_Z_mixed4']
                if group!='pressure':grid=grid[name]
                for key,row in grid.items():
                    k,n=(int(v[1:]) for v in key.split('_'));contains(row,source['targets'][name,letter,k,n](xx,zz));count+=1
    raw=packet['raw_current_radius_y_derivative_axial_coefficients']
    for name,rows in (*raw['histories'].items(),('P',raw['absolute_pressure']),*raw['velocity'].items()):
        for k,row in enumerate(rows):
            for n in range(5-k):contains(row[n]*math.factorial(n),source['raw_targets'][name,k,n](xx,zz));rawcount+=1
    assert len(packet['actual_Q_x_derivative_axial4'])==5 and all(len(row)==5 for row in packet['actual_Q_x_derivative_axial4'])
    return dict(passed=True,x=xx,Z=zz,independent_closed_integral_physical_x_y_R_mixed4_comparisons=count,
        independent_raw_current_history_pressure_velocity_mixed4_comparisons=rawcount,
        nonzero_current_partial_primitives_not_Rm_reset=True,separate_pressure_datum_with_Z4_derivative=True,
        finite_polynomial_basis_fixture_not_native_beta_parameter_selection=True,
        numerical_reference_comparison_tolerance_only_not_interval_certificate=True)


def unit(op,name,letter,k):
    f=op.flow;c=op.c;Rm=op.Rm_factor
    powers=dict(Utheta=0,Uz=0,Ur=.5,P=0,Mz=1,Mtheta=1.5,Mtheta_z=1.5,Mztheta=1,Mp=0)
    amplitude=dict(Utheta=op.Pstar,Uz=f.scalar(1),Ur=f.scalar(1/c.sqrt(2)),P=op.P2,
        Mz=f.scalar(1),Mtheta=op.Pstar*c.sqrt(2),Mtheta_z=op.Pstar*c.sqrt(2),Mztheta=op.P2,Mp=op.P2)
    return amplitude[name]*f.radial_power(Rm,powers[name]-(k if letter=='R' else 0))


def normalized_cover(f,row):
    if not row.zero and ep(row.scale.evaluate())[1]>1000:
        raise ArithmeticError('Normalized source comparison attempted an enormous exponential')
    return f.ordinary_cover(row)


def same_row(a,b):
    return a.scale.powers==b.scale.powers and a.scale.offset._mpi_==b.scale.offset._mpi_ and a.coefficient._mpi_==b.coefficient._mpi_


def native_checks(report,owner):
    records=cellchecks=memory=P0=flat=joins=rawrows=unitchecks=0
    with mp.workdps(owner.c.dps+40):
        for label in ('0','.5'):
            op=owner.owner(label);f=op.flow;source=op.op
            units={}
            inlet=op.evaluate((1,1));parent=source.evaluate((1,1));Rh=op.evaluate('Rh')
            assert inlet['original_P0_normalized_axial5'] is source.P0 and Rh['original_P0_normalized_axial5'] is source.P0;P0+=2
            for name,pname in (('mass','mean'),('theta','theta'),('mixed','mixed'),('energy','centered_energy'),('pressure','pressure')):
                assert inlet['actual_partial_primitive_source_memory']['actual_point_partial_primitive_object_memory'][name] is parent['normalized_primitives'][pname];memory+=1
            for edge in (49,51,59,61,69,71):
                # Check the actual derivative interface, not unrelated costly
                # partial primitives at a directed near-support endpoint.
                x,_,_=op.coordinate((edge,40));q=Fraction(edge,40)
                rows=[current.gamma_rows(source,x,center,abs(q-center)==Fraction(1,40))
                      for center in (Fraction(5,4),Fraction(3,2),Fraction(7,4))]
                assert all(ep(v)==(0,0) for row in rows for v in row);flat+=15
            # Exact open zero-correction neighbourhood at Rm. Derivatives
            # and primitive own-rate transport follow the original source
            # equations; source closures, not overlaps, establish the join.
            assert all(v.zero for row in inlet['actual_g_x_derivative_axial5'] for v in row)
            assert Rh['actual_partial_primitive_source_memory']['exact_terminal_identity_of_same_leading_map'];joins+=2
            for packet in report['frames'][label]['points']+report['frames'][label]['cells']:
                v=packet['function_evaluation'];assert not packet['whole_axis_functions_installed']
                assert v['original_P0_only_at_radial_order0'] and v['raw_radial_Z_order4_only_no_selected_Z5']
                assert v['shared_Rm_offset_never_independently_subtracted_in_physical_grid']
                assert not v['actual_patch_finite_N_density_oracle_installed'] and not v['whole_axis_function_provider_or_finite_N_Rc_patch_installed']
                for group in ('velocity','five_primitive'):
                    for letter in ('x','y','R'):
                        for grid in v['physical_'+group+'_'+letter+'_Z_mixed4'].values():assert len(grid)==15;records+=15
                for letter in ('x','y','R'):assert len(v['physical_pressure_'+letter+'_Z_mixed4'])==15;records+=15
                raw=v['raw_current_radius_y_derivative_axial_coefficients'];assert all(len(row)==5 for row in raw['velocity']['radial']);rawrows+=25
                assert len(v['physical_pressure_axis_axial5'])==6 and len(v['original_P0_normalized_axial5'])==6;P0+=1
            choices=(((6,5),(13,10),(5,4)),((29,20),(31,20),(3,2)),((17,10),(9,5),(7,4)),((71,40),'Rh',(15,8)))
            for a,b,q in choices:
                cell=op.cell(a,b);point=op.evaluate(q)
                assert cell['original_P0_normalized_axial5'] is source.P0;P0+=1
                for group,names in (('velocity',('Utheta','Uz','Ur')),('pressure',('P',)),('five_primitive',('Mz','Mtheta','Mtheta_z','Mztheta','Mp'))):
                    for letter in ('x','y','R'):
                        for name in names:
                            cg=cell['physical_'+group+'_'+letter+'_Z_mixed4'];pg=point['physical_'+group+'_'+letter+'_Z_mixed4']
                            if group!='pressure':cg=cg[name];pg=pg[name]
                            for key,row in pg.items():
                                k,n=(int(v[1:]) for v in key.split('_'));ukey=(name,letter,k)
                                if ukey not in units:units[ukey]=unit(source,name,letter,k)
                                u=units[ukey]
                                coordinate='y' if letter=='y' else 'x';commonkey=coordinate+str(k)+'_Z'+str(n)
                                cn=cell['physical_common_unit_normalized_'+coordinate+'_Z_mixed4'][name][commonkey]
                                pn=point['physical_common_unit_normalized_'+coordinate+'_Z_mixed4'][name][commonkey]
                                # The original shared unit is attached ONCE.
                                # This is an algebraic prefactor check; no huge
                                # interval log is subtracted from itself.
                                assert same_row(cn*u,cg[key]) and same_row(pn*u,row),(label,name,key,'shared unit');unitchecks+=2
                                if name=='P':
                                    # Identical P0 translation preserves cell
                                    # inclusion; compare the finite increment.
                                    assert cell['original_P0_normalized_axial5'] is point['original_P0_normalized_axial5']
                                    cn=cell['physical_common_unit_normalized_'+coordinate+'_Z_mixed4']['Mp'][commonkey]
                                    pn=point['physical_common_unit_normalized_'+coordinate+'_Z_mixed4']['Mp'][commonkey]
                                cb=ep(normalized_cover(f,cn));pb=ep(normalized_cover(f,pn))
                                assert cb[0]<=pb[0]<=pb[1]<=cb[1],(label,a,b,name,key);cellchecks+=1
    return dict(passed=True,actual_factored_physical_x_y_R_mixed4_records=records,
        closed_radial_cell_contains_native_point_enclosure_comparisons=cellchecks,current_partial_primitive_object_joins=memory,
        original_P0_separate_object_and_record_checks=P0,exact_support_edge_beta_derivative_zero_checks=flat,
        exact_Rm_Rh_source_join_checks=joins,raw_radial_records_without_Z5=rawrows,
        canonical_normalization_applied_before_cell_comparison=True,
        algebraic_shared_positive_physical_unit_checks=unitchecks,
        huge_Rm_offset_not_independently_subtracted_or_materialized=True,
        pressure_cell_inclusion_uses_increment_and_exact_same_P0_translation=True)


def guards(owner):
    n=0
    def reject(fn):
        nonlocal n
        try:fn()
        except (ValueError,TypeError):n+=1
        else:raise AssertionError('Invalid mixed4/cell query accepted')
    op=owner.owner('0')
    for q in ((0,1),(3,1),(1,0),[1,1],'Rm'):reject(lambda q=q:op.evaluate(q))
    for a,b in (((1,1),(1,1)),((2,1),(1,1)),((0,1),(2,1)),((1,1),(3,1)),([1,1],(2,1)),('Rh',(2,1))):reject(lambda a=a,b=b:op.cell(a,b))
    for label in ('whole_Z','1','-.5'):reject(lambda label=label:owner.owner(label))
    bad=SimpleNamespace(**op.op.__dict__);bad.P0=list(op.op.P0);reject(lambda:current._MixedPatchOwner(bad))
    return dict(passed=True,invalid_coordinate_cell_source_or_P0_requests_rejected=n)


def run():
    began=time.monotonic();report=json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()));refs=[]
    with mp.workdps(140):
        source=symbolic_sources()
        for x,z in (('1.3','.3'),('1.7','.5')):
            refs.append(fixture(source,x,z));print('Independent closed primitive x/y/R mixed4 and raw units PASS',x,z,flush=True)
    owner=current.OriginalRmPatchMixed4Cells(require_checked=False);native=native_checks(report,owner);guard=guards(owner)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name):fields.previous.bind(owner.hashes,name,current.sha(name))
    assert report[current.GATE] and report['source_family']==owner.family and report['original_source_bindings']['passed']
    assert all(report[key] is False for key in fields.previous.OPEN)
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        exact_closed_primitive_Stirling_and_radial_shift_identities=source['exact_closed_primitive_and_coordinate_identities'],
        independent_closed_integral_and_raw_unit_fixtures=refs,genuine_actual_source_mixed4_cells=native,guards=guard,
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual leading Rm patch full mixed4, closed cells, raw units and source joins PASS',flush=True);return result


if __name__=='__main__':run()
