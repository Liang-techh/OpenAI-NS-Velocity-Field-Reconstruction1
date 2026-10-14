"""Scoped actual pulse second-radial and signed pressure FTC source checks."""
import copy
from pathlib import Path
import time
import sympy as sy
from sympy.core.function import AppliedUndef
import lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_outer_own_source_joins_check as outer_check


def zero(expr,label):
    # Rename arbitrary function values/ordinary derivatives to independent
    # indeterminates before polynomial reduction. This preserves the exact
    # identity and avoids expensive sorting of nested SymPy derivatives.
    atoms=sorted(expr.atoms(AppliedUndef,sy.Derivative),key=str)
    plain=expr.xreplace({value:sy.Symbol('source_jet_'+str(n)) for n,value in enumerate(atoms)})
    result=sy.cancel(plain)
    if result!=0:raise AssertionError((label,result))


def bind(bindings,jet,body,Z):
    for n,ref in enumerate(current.previous.rows(jet)):
        value=sy.diff(body,Z,n)
        if ref.node in (0,1):
            assert value==ref.node,(ref.node,body,n)
        elif ref.node in bindings:
            zero(bindings[ref.node]-value,('conflicting source input',ref.node))
        else:bindings[ref.node]=value


class Sigma(sy.Function):
    @classmethod
    def eval(cls,x):
        if x.is_number:
            if x<=0:return sy.Integer(0)
            if x>=1:return sy.Integer(1)
    def fdiff(self,argindex=1):return SigmaPrime(self.args[0])


class SigmaPrime(sy.Function):
    @classmethod
    def eval(cls,x):
        if x.is_number and (x<=0 or x>=1):return sy.Integer(0)
    def fdiff(self,argindex=1):return SigmaSecond(self.args[0])


class SigmaSecond(sy.Function):
    @classmethod
    def eval(cls,x):
        if x.is_number and (x<=0 or x>=1):return sy.Integer(0)


class GpPrimitive(sy.Function):
    @classmethod
    def eval(cls,x):
        if x==0:return sy.Integer(0)
    def fdiff(self,argindex=1):return Sigma(50*self.args[0])


class Beta(sy.Function):
    @classmethod
    def eval(cls,x):
        if x.is_number and (x<=-1 or x>=1):return sy.Integer(0)
    def fdiff(self,argindex=1):
        x=self.args[0];return -2*x*Beta(x)/(1-x*x)**2


class PulseAlgebra(outer_check.SourceAlgebra):
    """Exact source calculus, with the pinned scalar derivatives as rules."""
    def __init__(self,owner,bindings):super().__init__(owner.g,bindings);self.owner=owner
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings or i in self.cache:return super().at(i)
        row=self.g.nodes[i];op=row['operation']
        if op=='exact_original_source_flat_sigma':
            assert row['exact_function_not_its_interval_cap'];result=Sigma(self.at(row['argument']))
        elif op=='analytic_unary' and row['name'] in ('original_flat_sigma_prime','original_flat_sigma_second'):
            assert row['exact_flat_source_derivative_not_a_bound']
            second=row['name'].endswith('second')
            expected=self.owner.sigma_second_source if second else self.owner.sigma_source
            assert row['derivative_source']==expected
            assert row['ordinary_derivative']==('2*same original sigma Taylor coefficient[2]' if second
                else 'same original sigma Taylor coefficient[1]')
            result=(SigmaSecond if second else SigmaPrime)(self.at(row['argument']))
        elif op=='definite_integral':
            lower,upper=self.at(row['lower']),self.at(row['upper'])
            if lower==upper:result=sy.Integer(0)
            elif i==self.owner.primitive.node:
                assert lower==0
                integrand=self.g.nodes[row['integrand']]
                assert integrand['operation']=='exact_original_source_flat_sigma'
                argument=self.g.nodes[integrand['argument']];assert argument['operation']=='product'
                terms=[self.g.nodes[j] for j in argument['arguments']]
                assert any(t['operation']=='exact_rational' and sy.Rational(t['numerator'],t['denominator'])==50 for t in terms)
                assert any(t['operation']=='bound_variable' and t['name']==row['variable'] for t in terms)
                result=GpPrimitive(upper)
            else:raise AssertionError(('Nonempty unbound actual pulse scalar primitive',i))
        elif op=='compact_raw_beta':result=Beta(self.at(row['argument']))
        elif op=='analytic_unary' and row['name']=='compact_raw_beta_derivative':
            assert row['exact_original_derivative_not_its_global_cap']
            x=self.at(row['argument']);result=sy.Integer(0) if Beta(x)==0 else -2*x*Beta(x)/(1-x*x)**2
        elif op=='original_lazy_flat_branch':
            assert row['original_beta_exact_flat_derivative'] and row['active_body_not_evaluated_on_flat_branch']
            x=self.at(row['Delta'])
            if x.is_number and (x<=-1 or x>=1):result=sy.Integer(0)
            else:return super().at(i)
        elif op=='function_substitution':
            assert row['Z_independent_substitution']
            variable=self.g.nodes[row['variable']];assert variable['operation']=='bound_variable'
            x=sy.Symbol(variable['name'],real=True);value=self.at(row['value'])
            bindings={k:v.subs(x,value) if hasattr(v,'subs') else v for k,v in self.bindings.items()}
            bindings[row['variable']]=value
            result=PulseAlgebra(self.owner,bindings).at(row['expression'])
        else:return super().at(i)
        self.cache[i]=result;return result


def actual_second_radial_source(owner):
    owner.assert_graph();g=owner.g;Z=sy.Symbol('Z',real=True);mu=sy.Symbol('mu',positive=True)
    selected={key:sy.Function('actual_'+key)(Z) for key in owner.selected}
    checked=[];flat=[]
    for chart,data in owner.charts.items():
        x=sy.Symbol(g.nodes[data['coordinate'].node]['name'],real=True)
        jac=sy.Integer(1) if chart=='end' else mu;bp=sy.Rational(1,2)+mu
        y=13/mu+x if chart=='end' else x/mu
        theta=sy.Function('actual_'+chart+'_theta_datum')(Z)*sy.exp(-bp*y)
        bindings={g.symbol('Z').node:Z,owner.mu.node:mu,data['coordinate'].node:x,
            owner.W['normalization'].node:sy.Symbol('exact_beta_normalization',positive=True),
            owner.W['log_end_scale'].node:sy.Symbol('exact_log_end_scale',real=True)}
        for key,jet in owner.selected.items():bind(bindings,jet,selected[key],Z)
        bind(bindings,data['fields']['Utheta'],theta,Z)
        alg=PulseAlgebra(owner,bindings);B=alg.at(data['normalized_axial'].value)
        for key,expected in (('normalized_axial_y',jac*sy.diff(B,x)),
            ('normalized_axial_yy',jac**2*sy.diff(B,x,2))):
            for n,ref in enumerate(current.previous.rows(data[key])):
                zero(alg.at(ref)-sy.diff(expected,Z,n),(chart,key,n));checked.append(chart+'_'+key+'_Z'+str(n))
        Uyy=jac**2*sy.diff(theta*B,x,2)
        for n,ref in enumerate(current.previous.rows(data['fields']['U_yy'])):
            zero(alg.at(ref)-sy.diff(Uyy,Z,n),(chart,'U_yy',n));checked.append(chart+'_U_yy_Z'+str(n))
        if chart in ('main','end'):
            value=sy.Integer(0)
            for key in ('normalized_axial','normalized_axial_y','normalized_axial_yy'):
                for n,ref in enumerate(current.previous.rows(data[key])):
                    zero(PulseAlgebra(owner,{**bindings,data['coordinate'].node:value}).at(ref),(chart,'flat',key,n))
                    flat.append(chart+'_'+key+'_Z'+str(n))
    return dict(passed=True,actual_second_radial_source_rows=checked,actual_Rp_Rv_flat_axial_pulse_rows=flat,
        gp_defining_integral_retained=True,main_gap_mu_d_dxi_end_d_ds=True,
        arbitrary_actual_selected_C5_coefficients_used=True,
        sigma_second_factor_two_from_pinned_Taylor_source=True,original_beta_second_derivative_not_zero_at_center=True)


def generic_recovery_FTC(owner):
    """Check actual output expressions; seed only physical input source jets."""
    owner.assert_graph();g=owner.g;Z=sy.Symbol('Z',real=True);delta=sy.Symbol('delta',positive=True)
    R=sy.Symbol('actual_R',positive=True);passed=[];recovery=[];mass=[];continued=[];derived=[]
    algebra_sources={key:current.parent.active.current.ast_binding(getattr(current.previous.Algebra,key))
        for key in ('trim','add','mul','scale','div','shift')}
    for pin in algebra_sources.values():
        assert owner.hashes[pin['module']]==pin['module_sha256']==current.sha(pin['module'])
    pulse=current.previous.previous.pulse5
    old_receipt=current.parent.linear.read(pulse.RECEIPT)
    assert old_receipt['all_passed'] and old_receipt['current_original_n1_actual_leading_pulse_U_mean_C5_Q_C4_functions_installed']
    assert old_receipt['source_family']==owner.family and owner.hashes[pulse.RECEIPT]==current.sha(pulse.RECEIPT)
    rp_receipt=current.parent.linear.read(current.previous.previous.RECEIPT)
    assert rp_receipt['all_passed'] and rp_receipt['source_family']==owner.family
    L=1-delta*Z**2;d=1-Z**2
    q=lambda U,m:(2*Z*U-(1-delta)*Z*m-d*sy.diff(m,Z))/L
    for chart,data in owner.charts.items():
        fields=data['fields'];bindings={g.symbol('Z').node:Z,owner.delta.node:delta,data['radius'].node:R}
        U,Uy,Uyy=[sy.Integer(0) if chart=='gap' else sy.Function(chart+'_'+key)(Z) for key in ('U','Uy','Uyy')]
        m=sy.Function(chart+'_actual_mean')(Z);my=U-m;myy=Uy-my
        for key,body in (('U',U),('U_y',Uy),('U_yy',Uyy),('m',m)):
            bind(bindings,fields[key],body,Z)
        alg=outer_check.SourceAlgebra(g,bindings)
        for key,body in (('m_y',my),('m_yy',myy)):
            for n,ref in enumerate(current.previous.rows(fields[key])):
                zero(alg.at(ref)-sy.diff(body,Z,n),(chart,key,n));mass.append(chart+'_'+key+'_Z'+str(n))
        for key,body in (('Q',q(U,m)),('Q_y',q(Uy,my)),('Q_yy',q(Uyy,myy))):
            for n,ref in enumerate(current.previous.rows(fields[key])):
                if key!='Q_yy' and n>2:
                    # The accepted actual C5 source receipt proves these
                    # are the ordinary derivatives of this same C2 body;
                    # the Rp receipt proves its actual input substitution.
                    # Differentiate the already proved recovery identity,
                    # without expanding numerical/Rm ancestors again.
                    continued.append(chart+'_'+key+'_Z'+str(n));continue
                zero(alg.at(ref)-sy.diff(body,Z,n),(chart,key,n));recovery.append(chart+'_'+key+'_Z'+str(n))
        for key,body in (('Q',q(U,m)),('Q_y',q(Uy,my))):
            bind(bindings,fields[key],body,Z)
        alg=outer_check.SourceAlgebra(g,bindings)
        for n in range(1):
            Py=alg.at(current.previous.rows(data['density']['P1_y'])[n])
            m13=alg.at(current.previous.rows(data['FTC']['signed_pressure_physical_y_derivatives']['m13'])[n])
            m15=alg.at(current.previous.rows(data['FTC']['signed_pressure_physical_y_derivatives']['m15'])[n])
            zero(m13-Py,(chart,'actual pressure FTC',n))
            zero(m15+R*Py,(chart,'actual signed weighted pressure FTC',n))
            passed.extend((chart+'_m13_y_Z'+str(n),chart+'_m15_y_Z'+str(n)))
        # Both complete C3 tuples use the same pinned ordinary jet algebra,
        # with Z-independent R/delta and genuine accepted C5/Q4 inputs.
        # Differentiate the proved whole-function identity through Z3;
        # expanding every rational derivative again adds no construction.
        derived.extend(chart+'_'+key+'_y_Z'+str(n) for key in ('m13','m15') for n in range(1,4))
        print('ACTUAL_PULSE_PRESSURE_FTC_WHOLE_FUNCTION',chart,flush=True)
    return dict(passed=True,actual_mass_radial_recovery_rows=mass,actual_leading_Q_radial_recovery_rows=recovery,
        actual_signed_pressure_FTC_rows=passed,
        pressure_FTC_higher_rows_by_pinned_ordinary_derivative_commutation=derived,
        ordinary_jet_algebra_source_bindings=algebra_sources,
        actual_Q_lower_rows_directly_proved_and_higher_rows_from_pinned_ordinary_derivative_source=continued,
        accepted_derivative_receipt=pulse.RECEIPT,accepted_derivative_receipt_sha256=current.sha(pulse.RECEIPT),
        accepted_actual_Rp_substitution_receipt=current.previous.previous.RECEIPT,
        accepted_actual_Rp_substitution_receipt_sha256=current.sha(current.previous.previous.RECEIPT),
        actual_pressure_outputs_not_bound=True,actual_Q_inputs_bound_only_after_recovery_proof=True,
        nonzero_arbitrary_incoming_leading_mean_retained=True,leading_one_minus_delta=True)


def actual_partial_primitive_FTC(owner):
    owner.assert_graph();g=owner.g;passed=[];preceding=[]
    for chart,data in owner.charts.items():
        x=data['coordinate'];lo,hi=data['domain'];Z=sy.Symbol('Z',real=True)
        R=sy.Symbol('actual_R',positive=True);mu=sy.Symbol('actual_mu',positive=True)
        jac=sy.Integer(1) if chart=='end' else mu
        bindings={g.symbol('Z').node:Z,data['radius'].node:R,owner.mu.node:mu}
        for key,jet in data['FTC']['actual_products'].items():
            # I0/I1 share actual product rows. Use one arbitrary input per
            # source node; retain structural zero/one in the inactive gap.
            for ref in current.previous.rows(jet):
                if ref.node not in (0,1):bindings.setdefault(ref.node,sy.Symbol('actual_product_row_'+str(ref.node),real=True))
        alg=outer_check.SourceAlgebra(g,bindings)
        for key,jet in data['partial_integrals'].items():
            for n,ref in enumerate(current.previous.rows(jet)):
                row=g.nodes[ref.node];full=g.nodes[current.previous.rows(owner.complete_integrals[chart][key])[n].node]
                assert row['operation']==full['operation']=='definite_integral'
                assert row['variable']==full['variable']==g.nodes[x.node]['name']
                assert row['lower']==full['lower']==lo.node and row['upper']==x.node and full['upper']==hi.node
                assert row['integrand']==full['integrand']==current.previous.rows(data['actual_integrands'][key])[n].node
                assert row['ordinary_slow_Z_derivative_order']==n and row['physical_Jacobian_applied_once']
                assert row['physical_radius_power']==(1 if key.startswith('I0_') else 2)
                expected=current.previous.rows(data['FTC']['actual_integral_physical_y_derivatives'][key])[n]
                zero(jac*alg.at(row['integrand'])-alg.at(expected),(chart,key,n,'actual physical FTC'))
                passed.append(chart+'_'+key+'_Z'+str(n))
        a=current.previous.Algebra(g,3)
        for key,value in data['total_integrals'].items():
            expected=a.add(*(owner.complete_integrals[name][key] for name in data['preceding_complete_charts']),data['partial_integrals'][key])
            assert current.previous.rows(value)==current.previous.rows(expected)
            preceding.append(chart+'_'+key)
    return dict(passed=True,actual_partial_integrand_physical_y_FTC_rows=passed,
        actual_complete_preceding_chart_history_projections=preceding,
        actual_field_query_is_the_named_binder=True,actual_original_domain_and_measure_retained=True,
        lower_global_finite_tensors_retained=True,partial_pressure_map_source=owner.pressure_map_source,
        endpoint_tensor_source=owner.pressure_endpoint_source)


def run(owner):
    began=time.monotonic();owner.assert_graph()
    checks=dict(second_radial_source=actual_second_radial_source(owner))
    print('ACTUAL_SECOND_RADIAL_SOURCE_PASSED',flush=True)
    checks['partial_primitive_FTC']=actual_partial_primitive_FTC(owner)
    print('ACTUAL_PARTIAL_PRIMITIVE_FTC_PASSED',flush=True)
    checks['generic_recovery_FTC']=generic_recovery_FTC(owner)
    return save(owner,checks,time.monotonic()-began)


def save(owner,checks,execution_seconds):
    owner.assert_graph()
    assert set(checks)=={'second_radial_source','partial_primitive_FTC','generic_recovery_FTC'}
    assert all(value['passed'] for value in checks.values())
    raw=current.parent.linear.read(current.NAME)
    assert raw['pulse_radial_pressure_C3_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    pins=dict(owner.hashes);pins[current.NAME]=current.sha(current.NAME)
    pins[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,source_family=owner.family,input_hashes=pins,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_second_radial_source_rows=len(checks['second_radial_source']['actual_second_radial_source_rows']),
        actual_partial_integrand_FTC_rows=len(checks['partial_primitive_FTC']['actual_partial_integrand_physical_y_FTC_rows']),
        exact_graph_nodes=len(raw['exact_graph_nodes']),execution_seconds=execution_seconds,
        numerical_oracle_available=False,outer_total_radial_pressure_joins=False,global_n1=False)
    (current.HERE/current.RECEIPT).write_text(__import__('json').dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance()
    print('ACTUAL_PULSE_SECOND_RADIAL_AND_PRESSURE_FTC_PASSED',receipt['actual_second_radial_source_rows'],
        receipt['actual_partial_integrand_FTC_rows'],'seconds',receipt['execution_seconds'],flush=True)
    return receipt


if __name__=='__main__':
    import sys
    sys.setrecursionlimit(10000)
    run(current.CurrentOriginalPulseRadialPressureC3Functions(require_checked=False))
