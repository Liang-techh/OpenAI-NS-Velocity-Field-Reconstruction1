"""Independent normalized-map algebra and direct finite bump/field fixtures."""
import copy
from fractions import Fraction
import json
import math
from pathlib import Path
import time
from types import SimpleNamespace
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_original_Rm_defect_patch_inverse as current
from lei_ren_part1_paper_compliant_current_original_bridge_macro_functions import MacroFlow
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

fields,ep=current.fields,current.ep


def exact_identities():
    t,a,I,z,P,x,Rm=sy.symbols('t a I z P x Rm',positive=True)
    h=sy.Matrix(sy.symbols('h0:5'));u=sy.Matrix(sy.symbols('u0:5'));d=sy.Matrix(sy.symbols('d0:5'))
    fg=sy.symbols('fg0:2');gg=sy.symbols('gg0:2');ff=sy.symbols('ff0:3');fx=sy.symbols('fx0:3')
    L=sy.Matrix([[1,1,0,0,0],[sy.Symbol('b0'),sy.Symbol('b1'),0,0,0],
        [0,0,*sy.symbols('c0:3')],[0,0,*sy.symbols('e0:3')],[0,0,*sy.symbols('p0:3')]])
    S=sy.diag(t,t,a,a,a)
    def nonlinear(v,iv,fgv,ffv,fxv):return sy.Matrix([0,fgv[0]*v[0]*v[2]+fgv[1]*v[1]*v[4],0,
        iv*(gg[0]*v[0]**2+gg[1]*v[1]**2)-sum(ffv[i]*v[i+2]**2 for i in range(3))/2,
        sum(fxv[i]*v[i+2]**2 for i in range(3))/2])
    Q=nonlinear(h,I,fg,ff,fx);Qn=nonlinear(u,I*t*t/a,[a*v for v in fg],[a*v for v in ff],[a*v for v in fx])
    count=0
    for v in L*S-S*L:assert sy.expand(v)==0;count+=1
    for v in S.inv()*Q.subs(dict(zip(h,S*u)))-Qn:assert sy.expand(v)==0;count+=1
    J=S.inv()*Q.jacobian(h).subs(dict(zip(h,S*u)))*S
    for v in J-Qn.jacobian(u):assert sy.expand(v)==0;count+=1
    assert sy.simplify((sy.exp(sy.Rational(6,5))*(1+z*z)**2/P**2)*t*t/(10**10*t*t/P**2)
        -sy.exp(sy.Rational(6,5))*(1+z*z)**2/10**10)==0;count+=1
    em,ek,eh,aa,eb,epres=sy.symbols('em ek eh aa eb epres')
    am=P*sy.exp(-sy.Rational(3,5))/(1+z*z);R=Rm*x
    # Derive incoming rows directly from physical reference primitives.
    H=sy.Rational(5,8)+eh;M=4*z+em;K=4*z*H+ek
    A=16*z*z+8*z*em+aa;B=sy.Rational(5,6)+eb;C=5+epres
    Mtheta=sy.sqrt(2)*Rm**sy.Rational(3,2)*am*H
    Mmixed=sy.sqrt(2)*Rm**sy.Rational(3,2)*am*K
    Mz=Rm*M;Mztheta=Rm*A-Rm*am**2*B/2;Mp=am**2*C/2
    derived=[Mz/Rm-4*z,(Mmixed-4*z*Mtheta)/(sy.sqrt(2)*Rm**sy.Rational(3,2)*am),
        Mtheta/(sy.sqrt(2)*Rm**sy.Rational(3,2)*am)-sy.Rational(5,8),
        (Mztheta-8*z*Mz+16*z*z*Rm)/(Rm*am**2)+sy.Rational(5,12),Mp/am**2-sy.Rational(5,2)]
    for got,want in zip(derived,(em,ek,eh,aa/am**2-eb/2,epres/2)):
        assert sy.simplify(got-want)==0;count+=1
    m,theta,mixed,energy,pm=sy.symbols('m theta mixed energy pm')
    physical=[R*m,sy.sqrt(2)*Rm**sy.Rational(3,2)*am*theta,
        sy.sqrt(2)*Rm**sy.Rational(3,2)*am*mixed,Rm*am**2*energy+8*z*R*m-16*z*z*R,am**2*pm]
    recovered=[physical[0]/R,physical[1]/(sy.sqrt(2)*Rm**sy.Rational(3,2)*am),
        physical[2]/(sy.sqrt(2)*Rm**sy.Rational(3,2)*am),
        (physical[3]-8*z*physical[0]+16*z*z*R)/(Rm*am**2),physical[4]/am**2]
    for got,want in zip(recovered,(m,theta,mixed,energy,pm)):assert sy.simplify(got-want)==0;count+=1
    return dict(passed=True,block_commutation_normalized_map_Jacobian_amplitude_defect_and_physical_unit_identities=count,
        exact_positive_formal_units_used_for_division=True,ordinary_weight_cover_may_include_zero_without_defining_scale=True)


def finite(row):
    return row.ctx.mpf(0) if row.zero else row.coefficient*row.ctx.exp(row.scale.evaluate())


def contains(row,value,reference_scale=1):
    # mp.diff/quad are finite numerical references, not directed intervals.
    # Reserve30 reference digits for their roundoff, including exact-zero
    # odd derivatives. This does not enlarge a production certificate.
    eps=mp.power(10,-mp.mp.dps+30)*max(1,abs(value),reference_scale)
    lo,hi=ep(finite(row));assert lo<=value+eps and value-eps<=hi,(mp.nstr(lo,15),mp.nstr(value,15),mp.nstr(hi,15))


def direct_bump_weights():
    radius=mp.mpf('.025');centers=tuple(mp.mpf(v) for v in ('1.25','1.5','1.75'))
    def b(s):return mp.mpf(0) if abs(s)>=1 else mp.exp(-1/(1-s*s))
    norm=mp.quad(b,[-1,0,1]);cache={}
    def weight(k,p,m=1,x=None):
        p=mp.mpf(p);upper=mp.mpf(1) if x is None else min(mp.mpf(1),max(mp.mpf(-1),(x-centers[k])/radius))
        key=(k,p,m,upper)
        if key not in cache:
            if upper==-1:out=mp.mpf(0)
            elif upper==1 and p==0 and m==1:out=mp.mpf(1)
            else:out=mp.quad(lambda s:(centers[k]+radius*s)**p*b(s)**m/(radius**(m-1)*norm**m),[-1,(upper-1)/2,upper])
            cache[key]=out
        return cache[key]
    def beta(k,x):
        off=x-centers[k];s=off/radius
        if abs(s)>=1:return mp.mpf(0),mp.mpf(0)
        value=b(s)/(radius*norm);return value,-2*off*value/(radius**2*(1-s*s)**2)
    return weight,beta,norm


def diagnostic(Z,accepted,weight,beta):
    c=MPIntervalContext();c.dps=160;Z=mp.mpf(Z)
    logRa=c.mpf(0);flow=MacroFlow(c,c.ln(c.mpf('1e-5')),c.ln(10**8),c.mpf(0),logRa,c.ln(100)-logRa,c.mpf(1)/16)
    z=IntervalTaylor.variable(c,c.mpf(Z),5);t=c.mpf(1)/8192;a=10**10*mp.mpf(ep(t)[0])**2/10**8
    W=current.weights(c,accepted.wraw)
    Wref={'L':[[mp.mpf(0) for _ in range(5)] for _ in range(5)],
        'fg':[weight(0,'.5',2),weight(2,'.5',2)],'gg':[weight(0,'0',2),weight(2,'0',2)],
        'ff':[weight(k,'0',2) for k in range(3)],'ff_over_x':[weight(k,'-1',2) for k in range(3)]}
    Wref['L'][0][:2]=[mp.mpf(1)]*2;Wref['L'][1][:2]=[weight(0,'.6'),weight(2,'.6')]
    for i in range(3):Wref['L'][2][i+2]=weight(i,'.5');Wref['L'][3][i+2]=-weight(i,'.1');Wref['L'][4][i+2]=weight(i,'-.9')
    for name,rows in Wref.items():
        for i,row in enumerate(rows):
            for j,v in enumerate(row if name=='L' else [row]):
                box=W[name][i][j] if name=='L' else W[name][i];lo,hi=ep(box);assert lo<=v<=hi
    funcs=(lambda v:mp.mpf('.001')*mp.cos(v),lambda v:mp.mpf('.002')*mp.sin(v),
        lambda v:-mp.mpf('.00003')*mp.exp(v),lambda v:mp.mpf('.00004')/(1+v*v),lambda v:-mp.mpf('.00002')*mp.cos(v))
    units=[mp.mpf(ep(t)[0])]*2+[a]*3
    def truth(v):return [fun(v)*units[i] for i,fun in enumerate(funcs)]
    def dfunc(v):
        h=truth(v);iv=mp.exp(mp.mpf('1.2'))*(1+v*v)**2/10**8
        Q=[0,Wref['fg'][0]*h[0]*h[2]+Wref['fg'][1]*h[1]*h[4],0,
            iv*sum(Wref['gg'][i]*h[i]**2 for i in range(2))-sum(Wref['ff'][i]*h[i+2]**2 for i in range(3))/2,
            sum(Wref['ff_over_x'][i]*h[i+2]**2 for i in range(3))/2]
        return [-sum(Wref['L'][i][j]*h[j] for j in range(5))-Q[i] for i in range(5)]
    def jet(fn):
        values=[mp.diff(fn,Z,n)/mp.factorial(n) for n in range(6)];eps=mp.mpf('1e-90')
        return IntervalTaylor(c,[c.mpf([v-eps,v+eps]) for v in values])
    d=[jet(lambda v,i=i:dfunc(v)[i]) for i in range(5)]
    poly=(1+z*z)*(1+z*z)*c.exp(c.mpf('1.2'));zero=[flow.scalar(0)]*6
    centered=dict(mean_error=flow.jet(d[0]),mixed_error=flow.jet(d[1]),angular_error=flow.jet(d[2]),
        axial_square=flow.scale(flow.multiply(flow.jet(d[3]),flow.jet(poly.reciprocal())),flow.factor((0,1,0,0,0))),
        swirl_error=zero,pressure_error=flow.scale(flow.jet(d[4]),2))
    P0=flow.jet(z*z*c.mpf('.13')+z*c.mpf('.07')+c.mpf('.11'))
    packet=dict(actual_centered_histories=centered,original_P0_normalized_axial5=P0,geometry=dict(exact_reference_offset=[-6,1]))
    ref=SimpleNamespace(flow=flow,z=z,Z=c.mpf(Z),delta=c.mpf('.07'),P0=P0,logC=c.ln(2)/10,postrestore=lambda offset:packet)
    bump=current.SharedFiveMomentRepair.__new__(current.SharedFiveMomentRepair);bump.ctx=c
    bump.normalization=current.restore_value(c,accepted.wraw['normalization'])
    bump.beta_sup=fields.previous.read_interval(c,accepted.fixed['bump_sup_bound']);bump.beta_deriv_sup=fields.previous.read_interval(c,accepted.fixed['bump_first_derivative_sup_bound'])
    bump.full_weights={(Fraction(r['center']),Fraction(r['power']),r['multiplicity']):current.restore_value(c,r['weight_interval']) for r in accepted.wraw['weight_records'].values()}
    op=current._ActualRmLeadingPatch(ref,t,W,bump);coeffs=fieldsrows=weightsrows=0
    for i,row in enumerate(op.inverse['controls']):
        for n,v in enumerate(row.coefficients):
            want=mp.diff(funcs[i],Z,n)/mp.factorial(n);lo,hi=ep(v);assert lo<=want<=hi;coeffs+=1
    for qt in ((1,1),(5,4),(3,2),(7,4),(2,1)):
        x=mp.mpf(qt[0])/qt[1];got=op.evaluate(qt)
        # These radial weights are Z-independent. mp.diff raises working
        # precision internally; it must not trigger fresh quadratures.
        partial_table={(k,p,m):weight(k,p,m,x) for k in range(3)
            for p,m in (('0',1),('.6',1),('.5',1),('.5',2),('0',2),('.1',1),('-.9',1),('-1',2))}
        partial=lambda k,p,m=1:partial_table[k,p,m]
        for key,v in got['actual_partial_weights'].items():
            import ast
            k,p,m=ast.literal_eval(key);want=partial(k,p,m);lo,hi=ep(v);assert lo<=want<=hi;weightsrows+=1
        def changes(v):
            h=truth(v);iv=mp.exp(mp.mpf('1.2'))*(1+v*v)**2/10**8
            return [sum(h[i]*partial(k,'0') for i,k in enumerate((0,2))),
                sum(h[i]*partial(k,'.6')+h[i]*h[k+2]*partial(k,'.5',2) for i,k in enumerate((0,2))),
                sum(h[k+2]*partial(k,'.5') for k in range(3)),
                iv*sum(h[i]**2*partial(k,'0',2) for i,k in enumerate((0,2)))
                -sum(h[k+2]*partial(k,'.1')+h[k+2]**2*partial(k,'0',2)/2 for k in range(3)),
                sum(h[k+2]*partial(k,'-.9')+h[k+2]**2*partial(k,'-1',2)/2 for k in range(3))]
        def rows(v):
            h=truth(v);dd=[dfunc(v)[i]+changes(v)[i] for i in range(5)] if x<mp.mpf(71)/40 else [mp.mpf(0)]*5
            amp=mp.exp(mp.mpf('-.6'))/(1+v*v);Rm=110*2*mp.mpf(10**4)**10*mp.exp(-6);R=Rm*x
            fc=sum(h[k+2]*beta(k,x)[0] for k in range(3));gc=sum(h[i]*beta(k,x)[0] for i,k in enumerate((0,2)))
            m=4*v+dd[0]/x;theta=mp.mpf(5)/8*x**mp.mpf('1.6')+dd[2];mixed=4*v*theta+dd[1]
            energy=dd[3]-mp.mpf(5)/12*x**mp.mpf('1.2');pm=dd[4]+mp.mpf(5)/2*x**mp.mpf('.2')
            # Row0 is linear: its explicit analytic derivative avoids
            # nested numerical differentiation inside every outer Z jet.
            mz=mp.mpf(4)
            if x<mp.mpf(71)/40:
                mz+=units[0]*(-mp.mpf('.001')*mp.sin(v)*(partial(0,'0')-1)
                    +mp.mpf('.002')*mp.cos(v)*(partial(2,'0')-1))/x
            V=4*v+gc;Q=(2*v*V-v*m*mp.mpf('.93')-(1-v*v)*mz)/(1-mp.mpf('.07')*v*v)
            Mp=10**8*amp*amp*pm;axis=10**8*(mp.mpf('.13')*v*v+mp.mpf('.07')*v+mp.mpf('.11'))
            return dict(velocity=dict(Ur=mp.sqrt(R/2)*Q,Utheta=10**4*amp*(x**mp.mpf('.1')+fc),Uz=V),
                first_y=dict(Utheta=10**4*amp*(x**mp.mpf('.1')/10+x*sum(h[k+2]*beta(k,x)[1] for k in range(3))),
                    Uz=x*sum(h[i]*beta(k,x)[1] for i,k in enumerate((0,2)))),
                primitive=dict(Mz=R*m,Mtheta=mp.sqrt(2)*Rm**mp.mpf('1.5')*10**4*amp*theta,
                    Mtheta_z=mp.sqrt(2)*Rm**mp.mpf('1.5')*10**4*amp*mixed,
                    Mztheta=Rm*10**8*amp*amp*energy+8*v*R*m-16*v*v*R,Mp=Mp),
                pressure=dict(axis=axis,increment=Mp,total=axis+Mp))
        for group,out in (('velocity',got['physical_velocity_axial_coefficients']),('first_y',got['physical_velocity_first_y_axial5']),('primitive',got['physical_cumulative_moment_axial5'])):
            for name,row in out.items():
                char=max(1,*(abs(rows(probe)[group][name]) for probe in (Z,Z-mp.mpf('.5'),Z+mp.mpf('.5'))))
                for n,rv in enumerate(row):contains(rv,mp.diff(lambda v:rows(v)[group][name],Z,n)/mp.factorial(n),char);fieldsrows+=1
        for name,key in (('axis','physical_pressure_axis_axial5'),('increment','physical_pressure_radial_increment_axial5'),('total','physical_total_pressure_axial5')):
            char=max(1,*(abs(rows(probe)['pressure'][name]) for probe in (Z,Z-mp.mpf('.5'),Z+mp.mpf('.5'))))
            for n,rv in enumerate(got[key]):contains(rv,mp.diff(lambda v:rows(v)['pressure'][name],Z,n)/mp.factorial(n),char);fieldsrows+=1
    return dict(passed=True,Z=Z,independent_manufactured_implicit_coefficient_Taylor_comparisons=coeffs,
        independent_direct_partial_bump_integral_comparisons=weightsrows,independent_physical_velocity_pressure_primitive_Taylor_comparisons=fieldsrows,
        finite_diagnostic_parameters_only_not_native_source_selection=True,
        numerical_reference_roundoff_budget='10^(-reference_dps+30)*max(1,abs(reference),finite_fixture_field_scale); only diagnostic comparison tolerance, not a source interval or rigorous reference-error proof',
        diagnostic_quadrature_is_not_rigorous_reference_certificate=True)


def genuine(report,owner):
    controls=defects=joins=P0=terminal=0;bounds={}
    with mp.workdps(owner.c.dps+40):
        for label in ('0','.5'):
            op=owner.owner(label);saved=report['frames'][label]['coefficient_solution'];inv=op.inverse['initial_C1_inverse']
            assert inv['certified'] and inv['self_map_strictly_inside'] and ep(inv['contraction_bound'])[1]<1
            bounds[label]=dict(initial_contraction=inv['contraction_bound'],higher_Jacobian_contraction=op.inverse['higher_inverse_contraction'])
            assert op.defects[0] is op.centered['mean_error'] and op.defects[1] is op.centered['mixed_error'] and op.defects[2] is op.centered['angular_error'];joins+=3
            assert ep(op.angularUnit.coefficient)[0]>0 and op.angularUnit.scale.powers==(0,-1,0,0,0)
            assert op.P0 is op.reference.P0 and op.P0 is op.Rm['original_P0_normalized_axial5'];P0+=2
            for i,row in enumerate(saved['actual_factored_coefficient_axial5']):
                for n,value in enumerate(row):
                    restored=current.endpoint.restore_row(op.flow,value);live=op.controls[i][n]
                    assert restored.coefficient._mpi_==live.coefficient._mpi_ and restored.scale.powers==live.scale.powers and restored.scale.offset._mpi_==live.scale.offset._mpi_;controls+=1
            for i,row in enumerate(saved['normalized_five_defect_axial5']):
                for n,value in enumerate(row['coefficients']):assert fields.previous.read_interval(owner.c,value)._mpi_==op.normalized_defects[i][n]._mpi_;defects+=1
            first=op.evaluate((1,1));assert all(v.zero for row in first['signed_bump_correction_sectors']['five_partial_primitive_changes'] for v in row)
            for i,row in enumerate(first['actual_recovered_five_defects']):
                assert all(v is op.defects[i][n] for n,v in enumerate(row));joins+=1
            for qt in ((71,40),(2,1)):
                packet=op.evaluate(qt);assert packet['exact_terminal_identity_of_same_unique_leading_map']
                assert all(v.zero for row in packet['actual_recovered_five_defects'] for v in row);terminal+=1
            for packet in report['frames'][label]['partial_patch']:
                assert not packet['whole_axis_functions_installed']
                value=packet['function_evaluation'];assert value['same_original_P0_retained'] and not value['whole_axis_leading_patch_or_finite_N_Rc_patch_certified']
                assert not value['full_radial_mixed4_or_finite_N_density_oracle_installed']
            assert not saved['finite_N_Rc_targets_or_controls_installed'] and not saved['whole_axis_implicit_family_certified']
    return dict(passed=True,actual_factored_control_Taylor_records=controls,actual_normalized_defect_Taylor_records=defects,
        exact_actual_Rm_inlet_row_joins=joins,original_P0_object_checks=P0,exact_full_weight_terminal_identity_checks=terminal,
        source_frame_contraction_bounds=bounds,no_native_source_value_selected=True,
        canonical_angular_unit_strictly_positive_and_outer_cover_never_divisor=True)


def guards(owner):
    n=0
    def reject(fn):
        nonlocal n
        try:fn()
        except (ValueError,TypeError):n+=1
        else:raise AssertionError('Invalid leading patch request accepted')
    op=owner.owner('0')
    for label in ('whole_Z','1','-.5'):reject(lambda label=label:owner.owner(label))
    for q in ((0,1),(3,1),(1,0),[1,1]):reject(lambda q=q:op.evaluate(q))
    for cells in (0,15,True):reject(lambda cells=cells:op.evaluate((1,1),cells))
    for t in (0,-1,op.c.mpf([1,2])):reject(lambda t=t:current._ActualRmLeadingPatch(op.reference,t,op.W,op.bump))
    bad={key:copy.deepcopy(value) for key,value in op.W.items()};bad['L'][0][2]=op.c.mpf(1)
    reject(lambda:current._ActualRmLeadingPatch(op.reference,op.t,bad,op.bump))
    reject(lambda:current.normalized_weights(op.c,op.W,op.c.mpf([-1,1])))
    return dict(passed=True,invalid_source_frame_coordinate_cells_scale_or_weight_requests_rejected=n)


def run():
    began=time.monotonic();report=json.loads((current.HERE/current.NAME).read_bytes());identities=exact_identities()
    owner=current.OriginalRmDefectPatchInverse(require_checked=False);refs=[]
    with mp.workdps(210):
        weight,beta,norm=direct_bump_weights()
        for Z in ('0','.5'):
            refs.append(diagnostic(Z,owner,weight,beta));print('Independent finite leading inverse and direct partial patch PASS',Z,flush=True)
    native=genuine(report,owner);guard=guards(owner)
    for name,digest in report['input_hashes'].items():fields.previous.bind(owner.hashes,name,digest)
    for name in (current.NAME,Path(__file__).name,'lei_ren_part1_paper_compliant_current_original_bridge_macro_functions.py'):
        fields.previous.bind(owner.hashes,name,current.sha(name))
    assert report[current.GATE] and report['source_family']==owner.family and report['original_source_bindings']['passed']
    assert all(report[key] is False for key in fields.previous.OPEN)
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,exact_normalization_defect_and_physical_unit_identities=identities,
        independent_finite_leading_inverse_and_partial_patch=refs,genuine_actual_Rm_source_and_implicit_solutions=native,
        guards=guard,**dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(current.base.encoded(result),indent=2)+'\n',encoding='utf8')
    print('Actual original Rm leading normalized inverse, partial primitives and physical units PASS',flush=True);return result


if __name__=='__main__':run()
