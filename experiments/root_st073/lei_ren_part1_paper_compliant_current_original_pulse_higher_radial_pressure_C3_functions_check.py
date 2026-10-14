"""Actual flat scalar calculus, higher radial recovery and pressure source."""
from pathlib import Path
import json
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_pulse_higher_radial_pressure_C3_functions as current
import lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions_check as radial_check
import lei_ren_part1_paper_compliant_current_original_n1_outer_radial_pressure_source_joins_check as join_check
import lei_ren_part1_paper_compliant_current_original_n1_outer_own_source_joins_check as outer_check

zero=join_check.zero;bind=radial_check.bind


def sigma_formula(x,n):
    s=Sigma(x);q=s*(1-s);a=1-2*s
    L1=2/x**3+2/(1-x)**3
    if n==1:return q*L1
    L2=-6/x**4+6/(1-x)**4
    if n==2:return q*(L2+a*L1**2)
    L3=24/x**5+24/(1-x)**5
    if n==3:return q*(L3+3*a*L1*L2+(1-6*q)*L1**3)
    if n!=4:raise ValueError('Actual scalar derivative order1..4 required')
    L4=-120/x**6+120/(1-x)**6
    return q*(L4+a*(4*L1*L3+3*L2**2)+6*(1-6*q)*L1**2*L2+a*(1-12*q)*L1**4)


class Sigma(sy.Function):
    @classmethod
    def eval(cls,x):
        if x.is_number:
            if x<=0:return sy.Integer(0)
            if x>=1:return sy.Integer(1)
    def fdiff(self,argindex=1):return sigma_formula(self.args[0],1)


class GpPrimitive(sy.Function):
    @classmethod
    def eval(cls,x):
        if x==0:return sy.Integer(0)
    def fdiff(self,argindex=1):return Sigma(50*self.args[0])


class FormalSigma(radial_check.Sigma):
    def fdiff(self,argindex=1):return FormalSigmaPrime(self.args[0])
class FormalSigmaPrime(radial_check.SigmaPrime):
    def fdiff(self,argindex=1):return FormalSigmaSecond(self.args[0])
class FormalSigmaSecond(radial_check.SigmaSecond):
    def fdiff(self,argindex=1):return FormalSigmaThird(self.args[0])
class FormalSigmaThird(radial_check.SigmaPrime):
    def fdiff(self,argindex=1):return FormalSigmaFourth(self.args[0])
class FormalSigmaFourth(radial_check.SigmaPrime):pass
class FormalGpPrimitive(GpPrimitive):
    def fdiff(self,argindex=1):return FormalSigma(50*self.args[0])
class FormalBeta(radial_check.Beta):
    def fdiff(self,argindex=1):return FormalBetaPrime(self.args[0])
class FormalBetaPrime(radial_check.Beta):
    def fdiff(self,argindex=1):return FormalBetaSecond(self.args[0])
class FormalBetaSecond(radial_check.Beta):
    def fdiff(self,argindex=1):return FormalBetaThird(self.args[0])
class FormalBetaThird(radial_check.Beta):
    def fdiff(self,argindex=1):return FormalBetaFourth(self.args[0])
class FormalBetaFourth(radial_check.Beta):pass


class HigherAlgebra(radial_check.PulseAlgebra):
    def __init__(self,owner,bindings,flat_nodes=None,formal_sigma=False):
        super().__init__(owner.pulse,bindings);self.high_owner=owner
        self.flat_nodes=flat_nodes if flat_nodes is not None else set()
        self.formal_sigma=formal_sigma
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings or i in self.cache:return super().at(i)
        row=self.g.nodes[i];op=row['operation']
        if op=='exact_original_source_flat_sigma':
            assert row['exact_function_not_its_interval_cap']
            result=(FormalSigma if self.formal_sigma else Sigma)(self.at(row['argument']))
        elif op=='analytic_unary' and row['name'] in ('original_flat_sigma_prime','original_flat_sigma_second'):
            n=2 if row['name'].endswith('second') else 1
            expected=self.owner.sigma_second_source if n==2 else self.owner.sigma_source
            assert row['derivative_source']==expected
            x=self.at(row['argument'])
            if self.formal_sigma:result=(FormalSigmaPrime if n==1 else FormalSigmaSecond)(x)
            else:result=sy.Integer(0) if x.is_number and (x<=0 or x>=1) else sigma_formula(x,n)
        elif op=='analytic_unary' and row['name'].startswith(('original_flat_sigma_order','original_compact_beta_order')):
            n=row['derivative_order'];assert n in (3,4)
            sigma=row['name'].startswith('original_flat_sigma')
            expected=self.high_owner.sigma_source if sigma else self.high_owner.beta_source
            assert row['derivative_source']==expected and row['ordinary_from_Taylor_factor']==__import__('math').factorial(n)
            assert row['exact_function_not_a_numerical_bound'] and row['active_body_not_evaluated_on_flat_branch']
            assert row['active_interval']==([0,1] if sigma else [-1,1]) and row['flat_value']==0
            if not sigma:assert row['polynomial_coefficients']==list(current.flat.BETA_POLYNOMIALS[n])
            x=self.at(row['argument']);lo,hi=row['active_interval']
            if x.is_number and (x<=lo or x>=hi):
                self.flat_nodes.add(i);result=sy.Integer(0)
            elif self.formal_sigma:
                derivative=(FormalSigmaThird if n==3 else FormalSigmaFourth) if sigma else (FormalBetaThird if n==3 else FormalBetaFourth)
                result=derivative(x)
            else:result=self.at(row['interior_body'])
        elif op=='compact_raw_beta' and self.formal_sigma:
            result=FormalBeta(self.at(row['argument']))
        elif op=='definite_integral' and i==self.owner.primitive.node:
            assert self.at(row['lower'])==0
            assert row['variable']=='original_pulse_gp_primitive_v'
            result=(FormalGpPrimitive if self.formal_sigma else GpPrimitive)(self.at(row['upper']))
        elif op=='function_substitution':
            assert row['Z_independent_substitution']
            variable=self.g.nodes[row['variable']];assert variable['operation']=='bound_variable'
            x=sy.Symbol(variable['name'],real=True);value=self.at(row['value'])
            bindings={k:v.subs(x,value) if hasattr(v,'subs') else v for k,v in self.bindings.items()}
            bindings[row['variable']]=value
            result=HigherAlgebra(self.high_owner,bindings,self.flat_nodes,self.formal_sigma).at(row['expression'])
        else:return super().at(i)
        self.cache[i]=result;return result


def actual_scalar_and_axial_source(owner):
    owner.assert_graph();g=owner.g;Z=sy.Symbol('Z',real=True);mu=sy.Symbol('mu',positive=True)
    selected={key:sy.Function('actual_'+key)(Z) for key in owner.pulse.selected}
    pins=(owner.sigma_source,owner.beta_source,owner.beta_polynomial_source)
    for pin in pins:assert owner.hashes[pin['module']]==pin['module_sha256']==current.sha(pin['module'])
    scalar_rows=[];actual=[];derived=[];flat=[]
    # Compare actual interior bodies with independent repeated derivatives
    # of the same exact scalar definition. Never evaluate the inactive
    # interior quotient at a support boundary.
    scalar_ids=[i for i,row in enumerate(g.nodes[:len(owner.prefix)]) if row['operation']=='analytic_unary'
        and row['name'].startswith(('original_flat_sigma_order','original_compact_beta_order'))]
    for i in scalar_ids:
        row=g.nodes[i];x=sy.Symbol('actual_flat_argument',real=True)
        alg=HigherAlgebra(owner,{row['argument']:x})
        sigma=row['name'].startswith('original_flat_sigma');n=row['derivative_order']
        expected=sy.diff(Sigma(x) if sigma else radial_check.Beta(x),x,n)
        zero(alg.at(i)-expected,('actual flat scalar interior',i,n));scalar_rows.append(i)
        lo,hi=row['active_interval']
        for value in (lo-1,lo,hi,hi+1):
            guard=HigherAlgebra(owner,{row['argument']:sy.Integer(value)})
            assert guard.at(i)==0 and i in guard.flat_nodes
    for name,data in owner.charts.items():
        x=sy.Symbol(g.nodes[data['coordinate'].node]['name'],real=True)
        jac=1 if name=='end' else mu;y=13/mu+x if name=='end' else x/mu;bp=sy.Rational(1,2)+mu
        theta=sy.Function('actual_'+name+'_theta_datum')(Z)*sy.exp(-bp*y)
        bindings={g.symbol('Z').node:Z,owner.pulse.mu.node:mu,data['coordinate'].node:x,
            owner.pulse.W['normalization'].node:sy.Symbol('actual_beta_normalization',positive=True),
            owner.pulse.W['log_end_scale'].node:sy.Symbol('actual_log_end_scale',real=True)}
        for key,jet in owner.pulse.selected.items():bind(bindings,jet,selected[key],Z)
        bind(bindings,data['fields']['Utheta'],theta,Z)
        # The scalar derivative bodies and guards have just been proved.
        # Consume their exact derivative relations in the Leibniz rule,
        # avoiding a second expansion of the same sigma rational bodies.
        alg=HigherAlgebra(owner,bindings,formal_sigma=True);B=alg.at(owner.pulse.charts[name]['normalized_axial'].value)
        for n in (3,4):
            suffix='yyy' if n==3 else 'yyyy'
            key='normalized_axial_'+suffix;jet=data[key]
            zero(alg.at(jet.value)-jac**n*sy.diff(B,x,n),(name,key,'actual whole radial source'))
            actual.append(name+'_'+key+'_Z0');derived.extend(name+'_'+key+'_Z'+str(k) for k in range(1,6))
        # All B0..4 derivative sources are now admitted. Compare the
        # actual velocity binomial body using those source inputs once.
        BX=sy.Integer(0) if name=='gap' else sy.Function('actual_'+name+'_B')(x,Z)
        ub=dict(bindings)
        Bjets=[owner.pulse.charts[name][key] for key in ('normalized_axial','normalized_axial_y','normalized_axial_yy')]
        Bjets.extend(data[key] for key in ('normalized_axial_yyy','normalized_axial_yyyy'))
        for k,jet in enumerate(Bjets):bind(ub,jet,jac**k*sy.diff(BX,x,k),Z)
        ualg=HigherAlgebra(owner,ub)
        for n in (3,4):
            key='U_'+('yyy' if n==3 else 'yyyy')
            zero(ualg.at(data['fields'][key].value)-jac**n*sy.diff(theta*BX,x,n),(name,key,'actual binomial source body'))
            actual.append(name+'_'+key+'_Z0');derived.extend(name+'_'+key+'_Z'+str(k) for k in range(1,6))
        # All four pulse support boundaries: Rp/main0, main11,
        # end-4 and Rv/end0. The gap is identically inactive.
        for edge,value in zip(('left','right'),data['domain']):
            if name=='gap':continue
            env=dict(bindings);env[data['coordinate'].node]=outer_check.SourceAlgebra(g,bindings).at(value)
            alg0=HigherAlgebra(owner,env)
            for suffix in ('yyy','yyyy'):
                for key in ('normalized_axial_'+suffix,'U_'+suffix):
                    jet=data[key] if key.startswith('normalized') else data['fields'][key]
                    zero(alg0.at(jet.value),(name,edge,key,'actual flat source'))
                    flat.append(name+'_'+edge+'_'+key+'_Z0')
        print('ACTUAL_THIRD_FOURTH_AXIAL_SOURCE',name,flush=True)
    ordinary=current.parent.linear.read(current.radial.RECEIPT)['checks']['generic_recovery_FTC']['ordinary_jet_algebra_source_bindings']
    for pin in ordinary.values():assert owner.hashes[pin['module']]==pin['module_sha256']==current.sha(pin['module'])
    return dict(passed=True,actual_scalar_interior_body_rows=scalar_rows,actual_whole_third_fourth_axial_rows=actual,
        higher_axial_Z_rows_by_pinned_ordinary_derivative_commutation=derived,actual_support_boundary_flat_rows=flat,
        actual_inactive_body_not_evaluated=True,ordinary_jet_algebra_source_bindings=ordinary,
        actual_original_gp_integral_and_finite_mu_retained=True,
        formal_scalar_derivative_rules_consumed_only_after_actual_body_proof=True,
        velocity_binomial_inputs_bound_only_after_B_derivative_source_proof=True,
        actual_scalar_ordinary_Taylor_conversion_factors=[6,24],end_ell_powers_for_order3_4=[4,5])


def actual_higher_mass_Q_source(owner):
    owner.assert_graph();g=owner.g;Y=sy.Symbol('physical_y',real=True);Z=sy.Symbol('Z',real=True)
    delta=sy.Symbol('delta',positive=True);L=1-delta*Z**2;d=1-Z**2;actual=[];derived=[];mass=[]
    for name,data in owner.charts.items():
        # The gap has structural U=0 and a nonzero exponential mean.
        H=sy.Function('actual_'+name+'_mean_datum')(Z)*sy.exp(-Y) if name=='gap' else sy.Function('actual_'+name+'_mass')(Y,Z)
        U=H+sy.diff(H,Y);bindings={g.symbol('Z').node:Z,owner.pulse.delta.node:delta}
        for key,k in (('U',0),('U_y',1),('U_yy',2),('U_yyy',3),('U_yyyy',4)):
            bind(bindings,data['fields'][key],sy.diff(U,Y,k),Z)
        for key,k in (('m',0),('m_y',1),('m_yy',2)):
            bind(bindings,data['fields'][key],sy.diff(H,Y,k),Z)
        alg=outer_check.SourceAlgebra(g,bindings)
        for key,k in (('m_yyy',3),('m_yyyy',4)):
            for n,ref in enumerate(current.outer.rows(data['fields'][key])):
                zero(alg.at(ref)-sy.diff(H,Y,k,Z,n),(name,key,n,'actual mass ODE source'));mass.append(name+'_'+key+'_Z'+str(n))
        Q=(2*Z*U-(1-delta)*Z*H-d*sy.diff(H,Z))/L
        for key,k in (('Q',0),('Q_y',1),('Q_yy',2),('Q_yyy',3),('Q_yyyy',4)):
            zero(alg.at(data['fields'][key].value)-sy.diff(Q,Y,k),(name,key,'actual whole Q radial recovery'))
            actual.append(name+'_'+key+'_Z0');derived.extend(name+'_'+key+'_Z'+str(n) for n in range(1,5))
        print('ACTUAL_HIGHER_MASS_Q_RADIAL_SOURCE',name,flush=True)
    return dict(passed=True,actual_mass_third_fourth_ODE_source_rows=mass,actual_whole_Q_radial_source_rows=actual,
        higher_Q_Z_rows_by_pinned_ordinary_derivative_commutation=derived,
        actual_Q_outputs_not_bound=True,actual_U3_U4_inputs_consumed_from_proved_axial_source=True,
        gap_structural_zero_U_and_nonzero_exponential_mean_retained=True,leading_one_minus_delta=True)


def actual_pressure_radial_source(owner):
    owner.assert_graph();g=owner.g;Y=sy.Symbol('physical_y',real=True);Z=sy.Symbol('Z',real=True)
    mu=sy.Symbol('mu',positive=True);delta=sy.Symbol('delta',positive=True);Rp=sy.Symbol('actual_Rp',positive=True)
    actual=[];derived=[];geometry=[]
    for name,data in owner.charts.items():
        x=sy.Symbol(g.nodes[data['coordinate'].node]['name'],real=True);jac=1 if name=='end' else mu
        rb={owner.pulse.mu.node:mu,owner.inlet.geometry['Rp'].node:Rp,data['coordinate'].node:x}
        realR=outer_check.SourceAlgebra(g,rb).at(data['radius'])
        for k in (1,2):
            zero(jac**k*sy.diff(realR,x,k)-realR,(name,'actual R radial derivative',k));geometry.append(name+'_R_y'+str(k))
        Q=sy.Function('actual_'+name+'_Q_datum')(Z)*sy.exp(-Y) if name=='gap' else sy.Function('actual_'+name+'_Q')(Y,Z)
        U=sy.Integer(0) if name=='gap' else sy.Function('actual_'+name+'_U')(Y,Z)
        R=sy.Symbol('actual_'+name+'_radial_origin',positive=True)*sy.exp(Y)
        bindings={g.symbol('Z').node:Z,owner.pulse.delta.node:delta,data['radius'].node:R}
        for key,k in (('U',0),('U_y',1),('U_yy',2)):
            bind(bindings,data['fields'][key],sy.diff(U,Y,k),Z)
        for key,k in (('Q',0),('Q_y',1),('Q_yy',2),('Q_yyy',3),('Q_yyyy',4)):
            bind(bindings,data['fields'][key],sy.diff(Q,Y,k),Z)
        alg=outer_check.SourceAlgebra(g,bindings);old=owner.pulse.charts[name]['density']
        S=alg.at(old['S'].value);Py=alg.at(old['P1_y'].value);Ny=alg.at(old['N1p'].value)
        expect=dict(S_y=sy.diff(S,Y),S_yy=sy.diff(S,Y,2),P1_yy=sy.diff(Py,Y),P1_yyy=sy.diff(Py,Y,2),
            N1p_y=sy.diff(Ny,Y),N1p_yy=sy.diff(Ny,Y,2),
            m15_yy=sy.diff(-R*Py,Y),m15_yyy=sy.diff(-R*Py,Y,2))
        for key,jet in data['higher_density'].items():
            zero(alg.at(jet.value)-expect[key],(name,key,'actual whole pressure radial source'))
            actual.append(name+'_'+key+'_Z0');derived.extend(name+'_'+key+'_Z'+str(n) for n in range(1,4))
        print('ACTUAL_PRESSURE_SECOND_THIRD_RADIAL_SOURCE',name,flush=True)
    return dict(passed=True,actual_radius_y_and_yy_source_identities=geometry,actual_whole_higher_pressure_source_rows=actual,
        higher_pressure_Z_rows_by_pinned_ordinary_derivative_commutation=derived,
        actual_Pyy_Pyyy_outputs_not_bound=True,actual_Q_inputs_bound_only_after_radial_recovery_proof=True,
        same_retained_own_pressure_FTC_source_consumed=True,all_R_derivatives_and_signed_terms_retained=True,
        R_y_equals_R_without_extra_mu=True)


def admitted_receipt(owner,module):
    """Consume an accepted ancestor theorem with its exact source pins."""
    rec=current.parent.linear.read(module.RECEIPT)
    assert rec['all_passed'] and rec['source_family']==owner.family
    assert owner.hashes[module.RECEIPT]==current.sha(module.RECEIPT)
    for name,digest in rec['input_hashes'].items():
        assert owner.hashes[name]==digest==current.sha(name),name
    return rec


def ordinary_source_contract(owner):
    rec=admitted_receipt(owner,current.radial)
    pins=rec['checks']['generic_recovery_FTC']['ordinary_jet_algebra_source_bindings']
    assert set(pins)=={'trim','add','mul','scale','div','shift'}
    for key,pin in pins.items():
        assert pin==current.parent.active.current.ast_binding(getattr(current.outer.Algebra,key))
        assert owner.hashes[pin['module']]==pin['module_sha256']==current.sha(pin['module'])
    for data in owner.charts.values():
        for prefix,count in (('U',6),('m',6),('Q',5)):
            for suffix in ('','_y','_yy','_yyy','_yyyy'):
                assert len(current.outer.rows(data['fields'][prefix+suffix]))==count
        assert all(len(current.outer.rows(jet))==4 for jet in data['higher_density'].values())
    return pins


def actual_quiet_higher_source(owner):
    owner.assert_graph();g=owner.g;q=owner.quiet['source'];mod=current.rp.quiet
    ordinary=ordinary_source_contract(owner)
    quiet_rec=admitted_receipt(owner,mod);compact_rec=admitted_receipt(owner,mod.previous)
    assert quiet_rec['checks']['native_mass_identity']['actual_compact_derived_K_equals_original_Rp_Pstar_kappa']
    assert quiet_rec['checks']['native_mass_identity']['actual_Rh_source_identity_and_same_original_native_mass_FunctionRef_used']
    assert compact_rec['checks']['actual_source']['original_compact_log_transport_measure_and_relative_only_terminal_theorem']
    # Only endpoint values of the correction are replaced by the admitted
    # terminal theorem. No correction in the compact interior is zeroed.
    _,c,q0,endpoint_inputs,normalizers=outer_check.quiet_start_context(owner.joins.outer)
    assert q0==q and q['incoming']==c['terminal']
    Z=sy.Symbol('Z',real=True);s=sy.Symbol(g.nodes[q['s']]['name'],real=True)
    mu=sy.Symbol('mu',positive=True);P=sy.Symbol('actual_Pstar',positive=True)
    R2=sy.Symbol('actual_R2',positive=True);kappa=sy.Symbol('actual_normalized_mass_coefficient',nonzero=True,real=True)
    length=sy.Symbol('actual_quiet_length',positive=True);delta=sy.Symbol('delta',positive=True)
    terminal_body={key:(kappa*Z if key=='m' else sy.Integer(0) if key=='V'
        else sy.Function('actual_terminal_'+key)(Z)) for key in ('E','V','m','h','k','e','p','P0')}
    # Every substitution below is restricted to x=2. The leading mass
    # family is the independently proved original linear-Z native source.
    inputs=dict(endpoint_inputs)
    for i in c['modified_fields']['m']:inputs.pop(i,None)
    for key,ids in c['background']['fields'].items():
        if key in terminal_body:bind(inputs,owner.jet(ids),terminal_body[key],Z)
    for key,certificate in c['terminal_certificates'].items():
        row=g.nodes[certificate]
        assert row['operation']=='proved_original_compact_relative_C5_terminal_identity'
        assert row['actual_endpoint_rows']==c['terminal']['correction'][key]
        assert row['equivalent_residual_rows']==c['terminal_rhs'][key]
        assert row['zero_rows']==[g.zero.node]*6 and row['relative_correction_only']
        assert row['absolute_leading_histories_and_pressure_not_zeroed']
        bind(inputs,owner.jet(c['history'][key]),sy.Integer(0),Z)
    alg=outer_check.EndpointSourceAlgebra(g,inputs);terminal=[];flat=[]
    for key in ('E','V','m','h','k','e','p','P0'):
        for n,(left,right) in enumerate(zip(q['incoming']['modified_fields'][key],q['leading_incoming'][key])):
            zero(alg.at(left)-alg.at(right),('actual modified-to-leading terminal',key,n))
            terminal.append(key+'_Z'+str(n))
    # Independently check the actual compact profile correction at the
    # terminal, rather than inferring it from a pressure or mass output.
    env=dict(inputs)
    x=g.nodes[q['incoming']['modified_fields']['m'][0]]['variable'];env[x]=sy.Integer(2)
    flat_alg=outer_check.EndpointSourceAlgebra(g,env)
    for key in ('delta_E','delta_V'):
        for n,i in enumerate(c['profiles']['radial_y_rows'][0][key]):
            zero(flat_alg.at(i),('actual compact flat profile correction',key,n));flat.append(key+'_Z'+str(n))
    bindings={g.symbol('Z').node:Z,q['s']:s,q['mu']:mu,q['Pstar']:P,q['R2']:R2,
        q['length']:length,q['delta']:delta}
    for key,body in terminal_body.items():
        bind(bindings,owner.jet(q['incoming']['modified_fields'][key]),body,Z)
        bind(bindings,owner.jet(q['leading_incoming'][key]),body,Z)
    alg=outer_check.SourceAlgebra(g,bindings);flow=[];physical=[]
    R=R2*sy.exp(s);K=R2*P*kappa
    zero(alg.at(q['K'])-K,'actual original quiet K recipe')
    zero(alg.at(q['R'])-R,'actual quiet radius recipe')
    for j in range(5):
        for key,packet in (('V','profiles_y0_to_y4'),('m','histories_y0_to_y4')):
            actual=q['flow'][packet][j][key];full_leading=q['leading_flow'][packet][j][key]
            assert len(actual)==6 and len(full_leading)==7
            leading=full_leading[:6]
            body=sy.Integer(0) if key=='V' else (-1)**j*kappa*Z*sy.exp(-s)
            for n,(left,right) in enumerate(zip(actual,leading)):
                zero(alg.at(left)-alg.at(right),('actual modified-leading quiet flow',j,key,n))
                zero(alg.at(left)-sy.diff(body,Z,n),('actual quiet flow source',j,key,n))
                flow.append('y'+str(j)+'_'+key+'_Z'+str(n))
    for prefix,base in (('U',sy.Integer(0)),('m',K*Z/R),('Q',-K/R)):
        for j,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy')):
            key=prefix+suffix;body=sy.diff(base,s,j)
            for n,ref in enumerate(current.outer.rows(owner.quiet['fields'][key])):
                zero(alg.at(ref)-sy.diff(body,Z,n),('actual quiet physical higher source',key,n))
                physical.append(key+'_Z'+str(n))
    expected=dict(S_y=K**2/R**2,S_yy=-2*K**2/R**2,
        P1_yy=-K**2/(4*R),P1_yyy=K**2/(4*R),
        N1p_y=-K**2/(2*R**2),N1p_yy=K**2/R**2,m15_yy=sy.Integer(0),m15_yyy=sy.Integer(0))
    pressure=[]
    for key,jet in owner.quiet['higher_density'].items():
        for n,ref in enumerate(current.outer.rows(jet)):
            zero(alg.at(ref)-sy.diff(expected[key],Z,n),('actual quiet higher pressure',key,n))
            pressure.append(key+'_Z'+str(n))
    # Keep the complete incoming pressure as an arbitrary, nonzero function.
    # Check its query independence before binding it in this source recipe.
    incoming={key:sy.Function('actual_incoming_'+key)(Z) for key in ('m13','m15')}
    for key,body in incoming.items():
        for i in q['incoming_pressure'][key]:
            assert not outer_check.free_query_dependency(owner.joins.outer,i,q['s'])
        bind(bindings,owner.jet(q['incoming_pressure'][key]),body,Z)
    alg=outer_check.SourceAlgebra(g,bindings);closed=[]
    increments=dict(m13=K**2*(1/R2-1/R)/4,m15=-K**2*s/4)
    for key,ids in q['closed_pressure'].items():
        for n,i in enumerate(ids):
            zero(alg.at(i)-sy.diff(incoming[key]+increments[key],Z,n),('actual full quiet pressure',key,n))
            closed.append(key+'_Z'+str(n))
    P1=alg.at(q['closed_pressure']['m13'][0])
    for n,key in ((2,'P1_yy'),(3,'P1_yyy')):
        zero(sy.diff(P1,s,n)-alg.at(owner.quiet['higher_density'][key].value),('actual retained full pressure derivative',key))
    return dict(passed=True,actual_compact_terminal_modified_leading_rows=terminal,
        actual_compact_flat_profile_correction_rows=flat,actual_quiet_modified_leading_flow_rows=flow,
        actual_quiet_physical_y0_to_y4_source_rows=physical,actual_quiet_higher_pressure_C3_rows=pressure,
        actual_retained_full_quiet_pressure_C3_rows=closed,actual_original_K_recipe_and_Pstar_once=True,
        normalized_terminal_mass_bound_only_after_original_linear_Z_theorem=True,
        actual_modified_and_leading_terminal_identity_before_flow_input_binding=True,
        compact_relative_zero_used_only_at_terminal=True,
        incoming_pressure_query_independence_before_arbitrary_function_binding=True,
        pressure_and_Q_outputs_not_bound=True,ordinary_source_derivative_factors_not_reapplied=True,
        admitted_original_native_mass_receipt=mod.RECEIPT,
        admitted_relative_terminal_receipt=mod.previous.RECEIPT,
        ordinary_jet_algebra_source_bindings=ordinary,
        exact_modified_C5_and_leading_C6_source_lengths_checked=True,
        fixed_actual_compact_normalizer_nodes=normalizers),dict(bindings=bindings,Z=Z,s=s,mu=mu,P=P,R2=R2,
            kappa=kappa,length=length,delta=delta,K=K)


def actual_internal_higher_pressure_joins(owner):
    owner.assert_graph();g=owner.g;rec=admitted_receipt(owner,current.previous)
    ordinary=ordinary_source_contract(owner)
    lower=rec['checks']['leading_internal_sources'];assert lower['passed']
    assert len(rec['checks']['own_second_radial_source']['actual_own_second_radial_source_rows'])==90
    Z=sy.Symbol('Z',real=True);mu=sy.Symbol('mu',positive=True);delta=sy.Symbol('delta',positive=True)
    pressure=[];flat=[];recovery=[];flat_derived=[];recovery_derived=[]
    definitions=(('main_gap','main',sy.Integer(11),'gap',sy.Integer(11)),
        ('gap_end','gap',13-4*mu,'end',sy.Integer(-4)))
    for seam,left_name,left_x,right_name,right_x in definitions:
        H=sy.Function('actual_'+seam+'_physical_mean')(Z);R=sy.Symbol('actual_'+seam+'_radius',positive=True)
        qbody=(-(1-delta)*Z*H-(1-Z**2)*sy.diff(H,Z))/(1-delta*Z**2)
        algs=[]
        for side,name,x in (('left',left_name,left_x),('right',right_name,right_x)):
            data=owner.charts[name];bindings={g.symbol('Z').node:Z,owner.pulse.mu.node:mu,
                owner.pulse.delta.node:delta,data['coordinate'].node:x,data['radius'].node:R,
                owner.pulse.W['normalization'].node:sy.Symbol('actual_beta_normalization',positive=True),
                owner.pulse.W['log_end_scale'].node:sy.Symbol('actual_log_end_scale',real=True)}
            for key,jet in owner.pulse.selected.items():bind(bindings,jet,sy.Function('actual_selected_'+key)(Z),Z)
            bind(bindings,data['fields']['Utheta'],sy.Function('actual_'+seam+'_theta')(Z),Z)
            flat_alg=HigherAlgebra(owner,bindings)
            for key in ('U_yyy','U_yyyy'):
                zero(flat_alg.at(data['fields'][key].value),(seam,side,key,'actual flat high axial source'))
                flat.append(seam+'_'+side+'_'+key+'_Z0')
                flat_derived.extend(seam+'_'+side+'_'+key+'_Z'+str(n) for n in range(1,6))
            # Consume only the already proved actual low-order fields,
            # radius identity, and the high-order flat source just checked.
            for j,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy')):
                bind(bindings,data['fields']['U'+suffix],sy.Integer(0),Z)
                if j<3:bind(bindings,data['fields']['m'+suffix],(-1)**j*H,Z)
            alg=outer_check.SourceAlgebra(g,bindings)
            for key,body in (('m_yyy',-H),('m_yyyy',H),('Q_yyy',-qbody),('Q_yyyy',qbody)):
                zero(alg.at(data['fields'][key].value)-body,(seam,side,key,'actual high mass Q recovery'))
                recovery.append(seam+'_'+side+'_'+key+'_Z0')
                recovery_derived.extend(seam+'_'+side+'_'+key+'_Z'+str(n)
                    for n in range(1,len(current.outer.rows(data['fields'][key]))))
            abstract_Q=sy.Function('actual_'+seam+'_recovered_Q')(Z)
            for j,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy')):
                bind(bindings,data['fields']['Q'+suffix],(-1)**j*abstract_Q,Z)
            algs.append(outer_check.SourceAlgebra(g,bindings))
        for key in ('P1_yy','P1_yyy','N1p_y','N1p_yy','m15_yy','m15_yyy'):
            left=owner.charts[left_name]['higher_density'][key];right=owner.charts[right_name]['higher_density'][key]
            for n,(lr,rr) in enumerate(zip(current.outer.rows(left),current.outer.rows(right))):
                zero(algs[0].at(lr)-algs[1].at(rr),(seam,key,n,'actual higher pressure source join'))
                pressure.append(seam+'_'+key+'_Z'+str(n))
        print('ACTUAL_INTERNAL_HIGHER_PRESSURE_JOIN',seam,flush=True)
    return dict(passed=True,actual_internal_higher_pressure_C3_rows=pressure,
        actual_high_axial_flat_whole_source_rows=flat,actual_high_mass_Q_whole_source_rows=recovery,
        higher_flat_Z_rows_by_pinned_ordinary_derivative_commutation=flat_derived,
        higher_recovery_Z_rows_by_pinned_ordinary_derivative_commutation=recovery_derived,
        ordinary_jet_algebra_source_bindings=ordinary,
        exact_input_row_lengths_and_ordinary_algebra_consumed_before_jet_bindings=True,
        actual_low_field_and_common_radius_lemmas=current.previous.RECEIPT,
        nonzero_mean_and_Z_dependent_Q_retained=True,actual_pressure_outputs_not_bound=True,
        high_Q_bound_only_after_actual_mass_recovery_source_proof=True)


def actual_outer_higher_pressure_joins(owner,quiet_proof,context):
    owner.assert_graph();assert quiet_proof['passed'];g=owner.g
    ordinary=ordinary_source_contract(owner)
    rec=admitted_receipt(owner,current.previous.previous)
    assert rec['checks']['Rp_endpoint_recovery']['passed'] and rec['checks']['Rv_endpoint_recovery']['passed']
    q=owner.quiet['source'];Z=context['Z'];mu=context['mu'];delta=context['delta'];P=context['P']
    K=context['K'];R=context['R2']*sy.exp(context['length']);M=K*Z/R
    pressure=[];recovery=[];flat=[];flat_derived=[];recovery_derived=[]
    # The quiet whole source was proved from actual normalized input rows.
    # Keep its saved endpoint substitution rather than resetting pressure.
    quiet_bindings=dict(context['bindings']);quiet_alg=outer_check.SourceAlgebra(g,quiet_bindings)
    data=owner.charts['main'];bindings={g.symbol('Z').node:Z,owner.pulse.mu.node:mu,owner.pulse.delta.node:delta,
        data['coordinate'].node:sy.Integer(0),data['radius'].node:R,
        owner.pulse.W['normalization'].node:sy.Symbol('actual_beta_normalization',positive=True),
        owner.pulse.W['log_end_scale'].node:sy.Symbol('actual_log_end_scale',real=True)}
    for key,jet in owner.pulse.selected.items():bind(bindings,jet,sy.Function('actual_selected_'+key)(Z),Z)
    bind(bindings,data['fields']['Utheta'],sy.Function('actual_Rp_theta')(Z),Z)
    alg=HigherAlgebra(owner,bindings)
    for key in ('U_yyy','U_yyyy'):
        zero(alg.at(data['fields'][key].value),('Rp actual flat high axial source',key));flat.append('Rp_'+key+'_Z0')
        flat_derived.extend('Rp_'+key+'_Z'+str(n) for n in range(1,6))
    for j,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy')):
        bind(bindings,data['fields']['U'+suffix],sy.Integer(0),Z)
        if j<3:bind(bindings,data['fields']['m'+suffix],(-1)**j*M,Z)
    alg=outer_check.SourceAlgebra(g,bindings)
    for key,body in (('m_yyy',-M),('m_yyyy',M),('Q_yyy',K/R),('Q_yyyy',-K/R)):
        zero(alg.at(data['fields'][key].value)-body,('Rp actual high mass Q recovery',key));recovery.append('Rp_'+key+'_Z0')
        recovery_derived.extend('Rp_'+key+'_Z'+str(n) for n in range(1,len(current.outer.rows(data['fields'][key]))))
    for j,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy')):
        bind(bindings,data['fields']['Q'+suffix],(-1)**(j+1)*K/R,Z)
    pulse_alg=outer_check.SourceAlgebra(g,bindings)
    for key,left in owner.outer_pressure_traces['Rp']['left'].items():
        right=owner.outer_pressure_traces['Rp']['right'][key]
        for n,(lr,rr) in enumerate(zip(current.outer.rows(left),current.outer.rows(right))):
            zero(quiet_alg.at(lr)-pulse_alg.at(rr),('actual Rp higher pressure source join',key,n))
            pressure.append('Rp_'+key+'_Z'+str(n))
    # The Rv side retains the accepted nonzero pressure value and its
    # query independence; only its radial derivatives are zero.
    data=owner.charts['end'];bindings={g.symbol('Z').node:Z,owner.pulse.mu.node:mu,owner.pulse.delta.node:delta,
        data['coordinate'].node:sy.Integer(0),data['radius'].node:sy.Symbol('actual_Rv',positive=True),
        owner.pulse.W['normalization'].node:sy.Symbol('actual_beta_normalization',positive=True),
        owner.pulse.W['log_end_scale'].node:sy.Symbol('actual_log_end_scale',real=True)}
    for key,jet in owner.pulse.selected.items():bind(bindings,jet,sy.Function('actual_selected_'+key)(Z),Z)
    bind(bindings,data['fields']['Utheta'],sy.Function('actual_Rv_theta')(Z),Z)
    alg=HigherAlgebra(owner,bindings)
    for key in ('U_yyy','U_yyyy'):
        zero(alg.at(data['fields'][key].value),('Rv actual flat high axial source',key));flat.append('Rv_'+key+'_Z0')
        flat_derived.extend('Rv_'+key+'_Z'+str(n) for n in range(1,6))
    for j,suffix in enumerate(('', '_y','_yy','_yyy','_yyyy')):
        bind(bindings,data['fields']['U'+suffix],sy.Integer(0),Z)
        if j<3:bind(bindings,data['fields']['m'+suffix],sy.Integer(0),Z)
    alg=outer_check.SourceAlgebra(g,bindings)
    for key in ('m_yyy','m_yyyy','Q_yyy','Q_yyyy'):
        zero(alg.at(data['fields'][key].value),('Rv actual high mass Q recovery',key));recovery.append('Rv_'+key+'_Z0')
        recovery_derived.extend('Rv_'+key+'_Z'+str(n) for n in range(1,len(current.outer.rows(data['fields'][key]))))
    for suffix in ('','_y','_yy','_yyy','_yyyy'):bind(bindings,data['fields']['Q'+suffix],sy.Integer(0),Z)
    alg=outer_check.SourceAlgebra(g,bindings)
    for key,left in owner.outer_pressure_traces['Rv']['left'].items():
        right=owner.outer_pressure_traces['Rv']['right'][key]
        for n,(lr,rr) in enumerate(zip(current.outer.rows(left),current.outer.rows(right))):
            zero(alg.at(lr)-alg.at(rr),('actual Rv higher pressure source join',key,n));pressure.append('Rv_'+key+'_Z'+str(n))
    assert owner.retained_tail_pressure is owner.joins.tail_pressure
    for ref in current.outer.rows(owner.retained_tail_pressure):
        assert not outer_check.free_query_dependency(owner.joins.outer,ref.node,owner.joins.tail_coordinate.node)
    return dict(passed=True,actual_Rp_Rv_higher_pressure_C3_rows=pressure,
        actual_outer_high_axial_flat_whole_source_rows=flat,actual_outer_high_mass_Q_whole_source_rows=recovery,
        higher_outer_flat_Z_rows_by_pinned_ordinary_derivative_commutation=flat_derived,
        higher_outer_recovery_Z_rows_by_pinned_ordinary_derivative_commutation=recovery_derived,
        ordinary_jet_algebra_source_bindings=ordinary,
        exact_input_row_lengths_and_ordinary_algebra_consumed_before_jet_bindings=True,
        actual_quiet_input_source_identity_consumed=True,
        actual_native_nonzero_linear_mass_and_shared_radius_lemmas=current.previous.previous.RECEIPT,
        actual_tail_pressure_value_retained_and_radial_query_independence_checked=True,
        pressure_outputs_not_bound=True,high_Q_bound_only_after_actual_mass_recovery=True)


CHECKS={'actual_scalar_and_axial_source','actual_higher_mass_Q_source','actual_pressure_radial_source',
    'actual_quiet_higher_source','actual_internal_higher_pressure_joins','actual_outer_higher_pressure_joins'}


def save(owner,checks,execution_seconds):
    owner.assert_graph();assert set(checks)==CHECKS and all(proof['passed'] for proof in checks.values())
    raw=current.parent.linear.read(current.NAME)
    assert raw['pulse_higher_radial_pressure_C3_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    pins={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)}
    rec=dict(all_passed=True,source_family=owner.family,input_hashes=pins,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_internal_higher_pressure_C3_rows=len(checks['actual_internal_higher_pressure_joins']['actual_internal_higher_pressure_C3_rows']),
        actual_Rp_Rv_higher_pressure_C3_rows=len(checks['actual_outer_higher_pressure_joins']['actual_Rp_Rv_higher_pressure_C3_rows']),
        exact_graph_nodes=len(raw['exact_graph_nodes']),numerical_oracle_available=False,
        execution_seconds=execution_seconds)
    assert rec['actual_internal_higher_pressure_C3_rows']==48 and rec['actual_Rp_Rv_higher_pressure_C3_rows']==48
    (current.HERE/current.RECEIPT).write_text(json.dumps(rec,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ACTUAL_PULSE_HIGHER_RADIAL_PRESSURE_AND_SOURCE_JOINS_CHECKED',flush=True)
    return rec


def run(owner):
    began=time.monotonic();checks={}
    for fn in (actual_scalar_and_axial_source,actual_higher_mass_Q_source,actual_pressure_radial_source):
        checks[fn.__name__]=fn(owner)
    checks['actual_quiet_higher_source'],context=actual_quiet_higher_source(owner)
    checks['actual_internal_higher_pressure_joins']=actual_internal_higher_pressure_joins(owner)
    checks['actual_outer_higher_pressure_joins']=actual_outer_higher_pressure_joins(owner,checks['actual_quiet_higher_source'],context)
    return save(owner,checks,time.monotonic()-began)


if __name__=='__main__':
    import sys
    sys.setrecursionlimit(max(10000,sys.getrecursionlimit()))
    run(current.CurrentOriginalPulseHigherRadialPressureC3Functions(require_checked=False))
