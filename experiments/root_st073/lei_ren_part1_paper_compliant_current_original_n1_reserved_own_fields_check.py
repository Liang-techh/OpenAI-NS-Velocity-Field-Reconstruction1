"""Source, fixed-forcing pressure, signed histories and n1 recovery checks."""
import json
from pathlib import Path
from types import SimpleNamespace
import time

import mpmath as mp
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_reserved_own_fields as current
import lei_ren_part1_paper_compliant_current_original_n1_full_defect_linear_functions_check as physical


def source_check(owner):
    owner.assert_graph();g=owner.g;a=current.parent.C1Algebra(g);f=owner.full
    assert owner.patch is f.repair_velocity
    assert owner.fields['M1']==owner.cumulative['m11'] and owner.fields['P1']==owner.cumulative['m13']
    for key in ('m11','m12','m14'):assert owner.incoming[key]==f.moments[key]
    for key in ('m13','m15'):
        assert owner.incoming[key]==a.add(f.segments['inner'][key],f.lower.functions[key],f.preRp.pieces['active'][key])
        assert owner.terminal[key]==a.add(owner.band_endpoint[key],f.preRp.pieces['quiet'][key],f.pulse.functions[key])
    complete,details=current.pressure_prefix_in_band(g,f,owner.s,g.constant(2))
    assert complete==f.preRp.pieces['compact']
    original=f.preRp.raw['defining_bindings']['compact']
    assert current.parent.encode(details['upper_fields'])==original['finite_endpoints'][1]
    assert current.parent.encode(details['integrals'])==original['physical_integrals']
    for key,pair in owner.patch_histories.items():
        for order,ref in enumerate((pair.value,pair.Z)):
            row=g.nodes[ref.node]
            assert row['operation']=='definite_integral' and row['lower']==g.zero.node
            assert row['upper']==owner.s.node and row['variable']==g.nodes[owner.patch_bodies['variable'].node]['name']
            assert row['integrand']==(owner.patch_bodies['bodies'][key].value,
                owner.patch_bodies['bodies'][key].Z)[order].node
            assert row['ordinary_slow_Z_derivative_order']==order
            assert row['history_not_reset_between_bumps'] and row['endpoint_independent_of_Z']
    assert type(owner.fields['Q1']) is current.source.FunctionRef
    assert type(owner.fields['V1']) is current.source.FunctionRef
    assert all(type(owner.fields[k]) is current.source.C1Function for k in ('mean1','M1','K1','P1'))
    assert f.reserved_background['nonlocal_mean_pressure_and_leading_N1p_not_zeroed']
    assert f.reserved_background['leading_swirl_equals_original_power_on_new_support']
    assert f.reserved_background['exact_new_and_old_supports']==current.full.linear.exact_supports()
    for q in ('0','3/10','1'):
        query=owner.at_band_fraction(q)
        assert query['C1_mean_M_K_P'] and query['Q_V_value_only']
        assert query['total_pressure_radial_derivative_not_installed']
    for q in ('-1/100','101/100'):
        try:owner.at_band_fraction(q)
        except ValueError:pass
        else:raise AssertionError('Out-of-band own-field query accepted')
    return dict(passed=True,actual_selected_patches_not_range_centers=True,
        ten_true_signed_history_integrals=True,original_axis_lower_active_pressure_datum_once=True,
        full_band_pressure_prefix_identical_to_original_four_C1_rows=True,
        all_sixteen_finite_endpoint_rows_and_twenty_integral_rows_identical=True,
        fixed_leading_N1p_retained_and_swirl_patch_added=True,
        C1_histories_and_C0_QV_only=True,exact_band_query_domain_enforced=True,
        original_leading_power_and_disjoint_support_source_guards_consumed=True,
        patch_pressure_radial_derivative_not_total=True)


def exact_density_check(owner):
    g=owner.g;f=owner.full;z=sy.Symbol('Z',real=True)
    Rc=sy.Symbol('Rc',positive=True);s=sy.Symbol('s',real=True);mu=sy.Symbol('mu',positive=True)
    U,T,A=(sy.Function(k)(z) for k in ('Uz1','Utheta1','Astar'))
    R=Rc*sy.exp(s);theta0=A*(2*R)**(-sy.Rational(1,2)-mu)
    bindings={owner.s.node:s,f.linear.background['logRc']:sy.log(Rc),
        f.linear.functions['mu'].node:mu,
        owner.patch['Uz1'].value.node:U,owner.patch['Uz1'].Z.node:sy.diff(U,z),
        owner.patch['Utheta1'].value.node:T,owner.patch['Utheta1'].Z.node:sy.diff(T,z),
        f.linear.background['Astar_ordinary_Z0_Z3'][0]:A,
        f.linear.background['Astar_ordinary_Z0_Z3'][1]:sy.diff(A,z)}
    expected=dict(m11=R*U,m12=R*sy.sqrt(2*R)*T,m13=theta0*T,
        m14=R*A*(2*R)**(-mu)*U,m15=-R*theta0*T)
    cache={};checked=[]
    for key,pair in owner.patch_bodies['per_ds'].items():
        for order,ref in enumerate((pair.value,pair.Z)):
            expression=physical.symbolic(g,ref,bindings,cache)
            difference=sy.expand_power_base(sy.powsimp(expression-sy.diff(expected[key],z,order),force=True),force=True)
            assert sy.simplify(difference)==0,(key,order,str(difference))
            checked.append(key+'_Z'+str(order))
    return dict(passed=True,actual_selected_density_graph_rows=checked,
        physical_dR_equals_R_ds=True,physical_F1_divided_by_sqrt_2R_once=True,
        m13_factor_two_and_m15_negative_weight=True,actual_Astar_Z_product_rule=True)


def exact_recovery_check(owner):
    g=owner.g;z=sy.Symbol('Z',real=True);de=sy.Symbol('delta',positive=True);R=sy.Symbol('R',positive=True)
    M,P,U=(sy.Function(k)(z) for k in ('M','P','U'))
    bindings={g.symbol('Z').node:z,owner.R.node:R,owner.full.cutoff.inner.parameters['delta'].node:de}
    for pair,expr in ((owner.cumulative['m11'],M),(owner.cumulative['m13'],P),(owner.patch['Uz1'],U)):
        bindings[pair.value.node]=expr;bindings[pair.Z.node]=sy.diff(expr,z)
    Q=(2*z*U-(1+de)*z*M/R-(1-z*z)*sy.diff(M,z)/R)/(1-de*z*z)
    expected=dict(mean1=M/R,K1=M/R-U,M1=M,P1=P,mean1_y=U-M/R,M1_y=R*U)
    cache={};checked=[]
    for key,expr in expected.items():
        pair=owner.fields[key]
        for order,ref in enumerate((pair.value,pair.Z)):
            assert sy.cancel(physical.symbolic(g,ref,bindings,cache)-sy.diff(expr,z,order))==0,(key,order)
            checked.append(key+'_Z'+str(order))
    for key,expr in (('Q1',Q),('V1',R*Q)):
        assert sy.cancel(physical.symbolic(g,owner.fields[key],bindings,cache)-expr)==0,key
        checked.append(key+'_value')
    return dict(passed=True,actual_own_recovery_graph_rows=checked,n1_mean_weight_one_plus_delta=True,
        Q_uses_full_mean_Z_not_only_K_Z=True,V_equals_R_Q_without_sqrt_factor=True,
        mean_y_plus_mean_equals_Uz=True,mass_y_equals_R_Uz_by_actual_FTC_integrand=True,
        coefficient_Z2_not_fabricated=True)


class HistoryInterpreter(physical.linear_check.Interpreter):
    """Fixture integration resolves compact supports instead of crossing them blindly."""
    def __init__(self,g,bindings,supports):
        super().__init__(g,bindings);self.supports=supports
    def at(self,node,local=None):
        if local is not None and node in local:return local[node]
        if node in self.bindings:return self.bindings[node]
        row=self.g.nodes[node];op=row['operation']
        if op=='function_substitution':
            value=self.at(row['value'],local)
            return self.at(row['expression'],{**(local or {}),row['variable']:value})
        if op=='definite_integral':
            var=next(i for i,r in enumerate(self.g.nodes) if r==dict(operation='bound_variable',name=row['variable']))
            lo,hi=self.at(row['lower'],local),self.at(row['upper'],local)
            if lo==hi:return mp.mpf(0)
            sign=1 if hi>lo else -1;lower,upper=sorted((lo,hi))
            breaks=sorted({lower,upper,*[v for v in self.supports if lower<v<upper]})
            return sign*mp.quad(lambda x:self.at(row['integrand'],{**(local or {}),var:x}),breaks)
        return super().at(node,local)


def finite_history_fixture():
    """Independent finite-mu partial integrals, nonzero incoming histories and gaps."""
    with mp.workdps(60):
        g=current.source.FunctionTransportGraph();c=g.constant;zr=g.symbol('Z');z=sy.Symbol('Z',real=True)
        Rc=mp.mpf('3.7');mu=mp.mpf('.08');delta=c('1/5')
        Av=sy.Rational(23,10)*(1+z*z)
        bind=lambda q:current.source.expression(g,sy.sympify(q),{'Z':zr})
        bg=dict(parameters={},logRc=g.unary('log',c('37/10')).node,
            log2Rc=g.unary('log',c('37/5')).node,
            Astar_ordinary_Z0_Z3=[bind(sy.diff(Av,z,k)).node for k in range(4)])
        # graph_weights needs the original mu FunctionRef, not a numeric matrix.
        muv=c('2/25');bg['parameters']['mu']=muv.node
        W=current.full.linear.graph_weights(g,bg)
        exprs=[sy.Rational(3,10)+z,-sy.Rational(1,5)+z*z,
            sy.Rational(2,5)-z,sy.Rational(3,20)+z*z,-sy.Rational(1,10)+z**3]
        pairs=[current.source.C1Function(bind(e),bind(sy.diff(e,z))) for e in exprs]
        patch=current.full.repair_velocity_functions(g,W,pairs,bg)
        fake=SimpleNamespace(repair_velocity=patch,linear=SimpleNamespace(background=bg,functions=W))
        histories,bodies=current.patch_histories(g,fake)
        a=current.parent.C1Algebra(g)
        incoming_M=sy.Rational(7,10)+sy.Rational(2,5)*z*z;incoming_P=sy.Rational(6,5)-z
        pair=lambda e:current.source.C1Function(bind(e),bind(sy.diff(e,z)))
        mass=a.add(pair(incoming_M),histories['m11']);pressure=a.add(pair(incoming_P),histories['m13'])
        fields=current.recovery_fields(g,patch,mass,pressure,delta,bodies['per_ds']['m13'])
        ell=mp.log(2)/400;centers=[j*mp.log(2)/40 for j in current.full.linear.CENTERS]
        beta=lambda t:mp.exp(-1/(1-t*t)) if abs(t)<1 else mp.mpf(0)
        integrate=lambda fn:mp.quad(fn,[-1,-.75,0,.75,1])
        normal=integrate(lambda t:beta(t)*mp.exp(ell*t/2))
        funcs=[sy.lambdify(z,e,'mpmath') for e in exprs];A=sy.lambdify(z,Av,'mpmath')
        def integral(j,k,s):
            cap=min(mp.mpf(1),max(mp.mpf(-1),(s-centers[j])/ell))
            if cap==-1:return mp.mpf(0)
            J=mp.quad(lambda t:beta(t)*mp.exp((k+1)*ell*t),[-1,*[v for v in (-.75,0,.75) if v<cap],cap])
            return mp.sqrt(2*Rc)*(2*Rc)**k*mp.exp((k+mp.mpf('.5'))*centers[j])*J/normal
        def direct(key,s,zz):
            cs=[fn(zz) for fn in funcs];aa=A(zz)
            if key=='m11':return sum(cs[i]*integral(j,0,s) for i,j in enumerate(current.full.linear.AXIAL))
            if key=='m14':return aa*sum(cs[i]*integral(j,-mu,s) for i,j in enumerate(current.full.linear.AXIAL))
            power=dict(m12=mp.mpf('.5'),m13=-mp.mpf('1.5')-mu,m15=-mp.mpf('.5')-mu)[key]
            factor=dict(m12=mp.mpf(1),m13=2*aa,m15=-aa)[key]
            return factor*sum(cs[i+2]*integral(j,power,s) for i,j in enumerate(current.full.linear.ANGULAR))
        points=[mp.mpf(0),centers[0],(centers[0]+ell+centers[1]-ell)/2,centers[2],mp.log(2)]
        supports=[v for center in centers for v in (center-ell,center,center+ell)]
        zz=mp.mpf('.371');worst=mp.mpf(0);rows=0
        for s in points:
            interp=HistoryInterpreter(g,{zr.node:zz,patch['log_R_over_Rc'].node:s,W['normal'].node:normal},supports)
            for key,p in histories.items():
                for order,ref in enumerate((p.value,p.Z)):
                    expected=mp.diff(lambda at:direct(key,s,at),zz,order)
                    error=abs(interp.at(ref.node)-expected);worst=max(worst,error);rows+=1
                    assert error<mp.mpf('1e-48'),(key,s,order,error)
            rr=Rc*mp.exp(s);M=incoming_M.subs(z,sy.Rational(str(zz)))
            M=mp.mpf(str(M))+direct('m11',s,zz)
            MZ=mp.mpf('.8')*zz+mp.diff(lambda at:direct('m11',s,at),zz)
            u=interp.at(patch['Uz1'].value.node)
            q=(2*zz*u-mp.mpf('1.2')*zz*M/rr-(1-zz*zz)*MZ/rr)/(1-mp.mpf('.2')*zz*zz)
            for key,expected in (('mean1',M/rr),('M1',M),('K1',M/rr-u),('Q1',q),('V1',rr*q)):
                ref=fields[key].value if type(fields[key]) is current.source.C1Function else fields[key]
                error=abs(interp.at(ref.node)-expected);worst=max(worst,error);rows+=1
                assert error<mp.mpf('1e-48'),(key,s,error)
        return dict(passed=True,mathematical_fixture_only=True,original_field_numerical_accuracy_not_claimed=True,
            finite_mu='.08',nonunit_Rc='3.7',original_physical_dR_and_signed_pressure_weights=True,
            independently_translated_partial_bump_integrals=True,ordinary_Z_product_rows=True,
            nonzero_incoming_mass_and_pressure=True,inactive_gap_keeps_cumulative_history=True,
            radial_points=len(points),checked_rows=rows,max_fixture_difference=str(worst))


def run(owner):
    began=time.monotonic();raw=current.parent.linear.read(current.NAME);owner.assert_graph()
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    assert raw['own_field_bindings']==owner.snapshot
    checks=dict(source=source_check(owner),densities=exact_density_check(owner),recovery=exact_recovery_check(owner))
    print('RESERVED_OWN_TRUE_HISTORIES_AND_RECOVERY_IDENTITIES_CHECKED',flush=True)
    checks['finite_histories']=finite_history_fixture()
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,**{current.GATE:True},
        **dict.fromkeys(current.OPEN,False),actual_field_numerical_accuracy_certified=False,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name)},execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_RESERVED_OWN_FIELDS_CHECKED',flush=True)
    return receipt


if __name__=='__main__':
    raise SystemExit('run(owner) reuses the retained true own-field graph')
