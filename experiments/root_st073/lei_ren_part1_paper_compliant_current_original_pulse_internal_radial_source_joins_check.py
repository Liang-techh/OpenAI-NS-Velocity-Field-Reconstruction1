"""Actual fixed-kernel transport and internal physical-y source joins."""
from pathlib import Path
import time
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_pulse_internal_radial_source_joins as current
import lei_ren_part1_paper_compliant_current_original_pulse_radial_pressure_C3_functions_check as radial_check
import lei_ren_part1_paper_compliant_current_original_n1_outer_radial_pressure_source_joins_check as join_check
import lei_ren_part1_paper_compliant_current_original_n1_outer_own_source_joins_check as outer_check

zero=join_check.zero;bind=radial_check.bind


class KernelAlgebra(radial_check.PulseAlgebra):
    def at(self,i):
        i=i.node if hasattr(i,'node') else i
        if i in self.bindings or i in self.cache:return super().at(i)
        row=self.g.nodes[i]
        if (row['operation']=='definite_integral' and row.get('ordinary_slow_Z_derivative_order',0)>0
            and row['variable'] in ('original_own_main_history_v','original_own_end_history_v','original_pulse_gp_primitive_v')):
            assert row['endpoint_independent_of_Z'] and row['derivative_of_same_defining_integral']
            assert row['actual_higher_source_body_not_interval_endpoint']
            # The genuine scalar geometry has no Z dependence. Expanded
            # higher source rows can retain unsimplified integral-of-zero
            # nodes; prove their actual integrand zero before using that fact.
            zero(self.at(row['integrand']),('actual scalar kernel derivative integrand',i))
            result=sy.Integer(0)
        elif row['operation']=='analytic_unary' and row['name']=='exprel':
            x=self.at(row['argument']);result=sy.Integer(1) if x==0 else (sy.exp(x)-1)/x
        elif row['operation']=='function_substitution':
            assert row['Z_independent_substitution']
            variable=self.g.nodes[row['variable']];assert variable['operation']=='bound_variable'
            x=sy.Symbol(variable['name'],real=True);value=self.at(row['value'])
            bindings={k:v.subs(x,value) if hasattr(v,'subs') else v for k,v in self.bindings.items()}
            bindings[row['variable']]=value
            result=KernelAlgebra(self.owner,bindings).at(row['expression'])
        else:return super().at(i)
        self.cache[i]=result;return result


def fixed_kernel_source(owner):
    owner.assert_graph();g=owner.g;W=owner.W
    mu=sy.Symbol('mu',positive=True);Z=sy.Symbol('Z',real=True);p=sy.Symbol('actual_logpref',real=True)
    norm=sy.Symbol('actual_beta_normalization',positive=True);t=sy.Symbol('original_pulse_beta_coordinate',real=True)
    ell=sy.Rational(3,20);lams=[sy.Rational(1,2)-i*mu for i in (1,2)]
    js=sy.symbols('actual_J1 actual_J2',positive=True);pis=sy.symbols('actual_Pi1 actual_Pi2',real=True)
    bindings={g.symbol('Z').node:Z,owner.pulse.mu.node:mu,W['normalization'].node:norm,
        W['log_common_prefactor'].node:p}
    alg=KernelAlgebra(owner.pulse,bindings)
    zero(alg.at(W['ell'])-ell,'actual ell')
    zero(alg.at(W['log_end_scale'])-(p-sy.log(mu)),'actual log end scale')
    for i,ref in enumerate(W['lambda_rows']):zero(alg.at(ref)-lams[i],('actual lambda',i))
    # Pin the existing exact arbitrary-function transport theorem and its
    # complete source dependencies; no finite fixture supplies a value.
    name=current.PREFIX+'pulse_interface_certificate.json'
    rec=current.parent.linear.read(name)
    assert owner.hashes[name]==current.sha(name) and rec['all_passed']
    assert rec['actual_five_defect_family_sha256']==owner.family['actual_five_defect_family_sha256']
    assert rec['implicit_source_sha256']==owner.family['implicit_source_sha256']
    needed=('actual_first_row_exponential_weight','actual_second_row_divided_difference',
        'actual_common_end_scale_normalization','actual_incoming_row_normalization',
        'main_gap_linear_moments_whole_Z','gap_end_full_future_linear_weight',
        'full_future_supports_at_minus4_and_empty_at_zero')
    assert all(rec['source_bound_functional_pulse_identities'][key] for key in needed)
    for file,digest in rec['input_hashes'].items():assert owner.hashes[file]==digest==current.sha(file)
    # Compare the actual first/divided beta bodies against their exact
    # defining recipes, then consume linearity/translation of that measure.
    kernel=current.parent.linear.read(current.radial.kernels.RECEIPT)
    assert kernel['all_passed'] and kernel['source_family']==owner.family
    assert kernel['checks']['direct_basis']['negative_determinant_identity_for_every_positive_mu']
    beta=radial_check.Beta(t);j_nodes=[];checks=[]
    for column,ref in enumerate(W['first_row']):
        row=g.nodes[ref.node];assert row['operation']=='product'
        quotient=[i for i in row['arguments'] if g.nodes[i]['operation']=='positive_quotient']
        assert len(quotient)==1
        q=g.nodes[quotient[0]];integral=g.nodes[q['numerator']]
        assert q['denominator']==W['normalization'].node
        assert integral['operation']=='definite_integral' and integral['variable']=='original_pulse_beta_coordinate'
        assert integral['lower']==g.constant(-1).node and integral['upper']==g.one.node
        zero(alg.at(integral['integrand'])-beta*sy.exp(lams[0]*ell*t),('actual J1 body',column))
        j_nodes.append(q['numerator']);bindings[q['numerator']]=norm*js[0]
    # The divided body is pointwise (row2-row1)/mu in the same full beta
    # measure. It defines J2 and not a separately rounded difference.
    for column,ref in enumerate(W['divided_row']):
        q=g.nodes[ref.node];assert q['operation']=='positive_quotient' and q['denominator']==W['normalization'].node
        integral=g.nodes[q['numerator']];assert integral['operation']=='definite_integral'
        assert integral['variable']=='original_pulse_beta_coordinate'
        assert integral['lower']==g.constant(-1).node and integral['upper']==g.one.node
        center=sy.Integer(-3 if column==0 else -1);s=center+ell*t
        exact=beta*(sy.exp(lams[1]*s)-sy.exp(lams[0]*s))/mu
        zero(KernelAlgebra(owner.pulse,bindings).at(integral['integrand'])-exact,('actual divided beta body',column))
        bindings[q['numerator']]=norm*(sy.exp(lams[1]*center)*js[1]-sy.exp(lams[0]*center)*js[0])/mu
        checks.append('actual_full_beta_column_'+str(column))
    alg=KernelAlgebra(owner.pulse,bindings)
    A=[alg.at(ref) for ref in W['first_row']];D=[alg.at(ref) for ref in W['divided_row']]
    # J2 occurs explicitly in the positive determinant; identify its
    # actual defining full beta integral, independently of the D result.
    closure=owner.joins.outer.operand_closure(W['positive_determinant_absolute'].node)
    for i in closure:
        row=g.nodes[i]
        if row['operation']=='definite_integral' and row['variable']=='original_pulse_beta_coordinate' and i not in bindings:
            assert row['lower']==g.constant(-1).node and row['upper']==g.one.node
            zero(alg.at(row['integrand'])-beta*sy.exp(lams[1]*ell*t),'actual J2 determinant body')
            bindings[i]=norm*js[1]
    alg=KernelAlgebra(owner.pulse,bindings)
    zero(alg.at(W['positive_determinant_absolute'])-(D[0]*A[1]-A[0]*D[1]),'actual positive divided determinant')
    xi=sy.Symbol('original_pulse_xi',real=True)
    gp=radial_check.GpPrimitive(xi)*radial_check.Sigma(11-xi)
    kernel_checks=[]
    for n,ref in enumerate(W['scaled_main_rows']):
        q=g.nodes[ref.node];assert q['operation']=='positive_quotient' and q['denominator']==owner.pulse.mu.node
        row=g.nodes[q['numerator']];assert row['operation']=='definite_integral'
        assert row['variable']=='original_pulse_xi' and row['lower']==g.zero.node and row['upper']==g.constant(11).node
        assert row['endpoint_independent_of_Z']
        zero(alg.at(row['integrand'])-sy.exp(lams[n]*(xi-13)/mu-p)*gp,('actual scaled main moment body',n))
        bindings[ref.node]=pis[n];kernel_checks.append('actual_scaled_main_row_'+str(n))
    # Only the specific source history kernels are recognized. A generic
    # nonempty definite integral still fails closed in KernelAlgebra.
    vmain=sy.Symbol('original_own_main_history_v',real=True);vend=sy.Symbol('original_own_end_history_v',real=True)
    x=sy.Symbol('original_own_pulse_main_xi',real=True);s=sy.Symbol('original_own_pulse_end_s',real=True)
    alg=KernelAlgebra(owner.pulse,bindings)
    for n,packets in enumerate(owner.history_kernels['main']):
        row=g.nodes[packets[0].node]
        assert row['actual_complete_own_pulse_history'] and row['endpoint_independent_of_Z']
        zero(alg.at(row['lower']),'actual main history lower')
        zero(alg.at(row['upper'])-x/mu,'actual main history upper')
        expected=sy.exp(-lams[n]*(x/mu-vmain))*radial_check.GpPrimitive(mu*vmain)*radial_check.Sigma(11-mu*vmain)
        zero(alg.at(row['integrand'])-expected,('actual forward main kernel',n))
        # xi=mu*v, dv=dxi/mu, fixed full support [0,11].
        converted=expected.subs(vmain,xi/mu)/mu
        scaled=sy.exp(lams[n]*(xi-13)/mu-p)*gp/mu
        zero(converted-sy.exp(lams[n]*(13-x)/mu+p)*scaled,('actual main fixed integral transport',n))
        kernel_checks.append('actual_forward_main_transport_'+str(n))
    end_columns={}
    for n,packets in enumerate(owner.history_kernels['end']):
        for ref in packets:
            row=g.nodes[ref.node]
            assert row['actual_complete_own_pulse_history'] and row['endpoint_independent_of_Z']
            zero(alg.at(row['lower'])-s,'actual backward end lower');zero(alg.at(row['upper']),'actual backward end upper')
            body=alg.at(row['integrand']);matches=[]
            for column,center in enumerate((-3,-1)):
                expected=sy.exp(lams[n]*(vend-s))*radial_check.Beta((vend-center)/ell)/(ell*norm)
                if sy.simplify(body-expected)==0:matches.append(column)
            assert len(matches)==1;column=matches[0];end_columns[ref.node]=(n,column)
            assert -4<sy.Rational(-3 if column==0 else -1)-ell and sy.Rational(-3 if column==0 else -1)+ell<0
            kernel_checks.append('actual_backward_end_full_support_'+str(n)+'_'+str(column))
    # Abstraction is now a proved scalar input lemma: retain the actual
    # determinant relation while avoiding repeated exponential expansion.
    aa=sy.symbols('actual_A0 actual_A1',positive=True);dd=sy.symbols('actual_D0 actual_D1',positive=True)
    for ref,value in zip(W['first_row'],aa):bindings[ref.node]=value
    for ref,value in zip(W['divided_row'],dd):bindings[ref.node]=value
    bindings[W['positive_determinant_absolute'].node]=dd[0]*aa[1]-aa[0]*dd[1]
    return dict(passed=True,actual_full_beta_recipe_columns=checks,actual_history_source_transport=kernel_checks,
        arbitrary_function_interface_receipt=name,arbitrary_function_interface_receipt_sha256=current.sha(name),
        actual_history_producer=owner.source_history_producer,actual_selector=owner.source_selector,
        actual_positive_finite_mu_divided_determinant=True,actual_log_end_scale_equals_logpref_minus_logmu=True,
        actual_main_change_of_variable_and_once_Jacobian=True,
        actual_end_full_support_at_minus4=True,numerical_fixture_not_used_as_source=True,
        scalar_input_abstraction_only_after_actual_recipe_and_determinant_proof=True),dict(
        bindings=bindings,Z=Z,mu=mu,p=p,A=aa,D=dd,pis=pis,lams=lams,end_columns=end_columns)


def selected_linear_source(owner,context):
    g=owner.g;bindings=dict(context['bindings']);Z=context['Z'];mu=context['mu'];p=context['p']
    a0=sy.Symbol('actual_native_swirl_amplitude',positive=True)
    h1,h2=sy.symbols('actual_native_m1_coefficient actual_native_m2_coefficient',real=True,nonzero=True)
    bodies=dict(u=a0/(1+Z**2),m1=h1*Z*(1+Z**2),m2=h2*Z*(1+Z**2))
    for key,body in bodies.items():bind(bindings,owner.inlet.current_frame[key],body,Z)
    ap=sy.Function('actual_selected_ap')(Z);bind(bindings,owner.pulse.selected['ap'],ap,Z)
    A,D=context['A'],context['D'];det=D[0]*A[1]-A[0]*D[1]
    alg=KernelAlgebra(owner.pulse,bindings)
    incoming=[bodies[key]*alg.at(f) for key,f in zip(('m1','m2'),owner.W['incoming_factors'])]
    rhs0=-mu*(incoming[0]+ap*context['pis'][0])
    rhs1=-(incoming[1]-incoming[0]+ap*(context['pis'][1]-context['pis'][0]))
    C=[(A[1]*rhs1-D[1]*rhs0)/det,(D[0]*rhs0-A[0]*rhs1)/det]
    checked=[]
    for key,body in zip(('C1','C2'),C):
        for n,ref in enumerate(current.outer.rows(owner.pulse.selected[key])):
            zero(alg.at(ref)-sy.diff(body,Z,n),('actual selected linear source',key,n));checked.append(key+'_Z'+str(n))
        bind(bindings,owner.pulse.selected[key],body,Z)
    zero(sum(A[j]*C[j] for j in range(2))-rhs0,'actual first selected row')
    zero(sum(D[j]*C[j] for j in range(2))-rhs1,'actual divided selected row')
    return dict(passed=True,actual_selected_C1_C2_C5_source_rows=checked,
        actual_first_and_divided_linear_equations=True,
        actual_ap_bound_as_admitted_input_not_numeric_midpoint=True,
        exact_current_native_input_projection_receipt=current.rp.RECEIPT,
        actual_end_scale_once=True),dict(context,bindings=bindings,native_bodies=bodies,C=C,ap=ap)


def leading_internal_sources(owner,context):
    owner.assert_graph();g=owner.g;Z=context['Z'];mu=context['mu'];p=context['p'];A=context['A'];D=context['D']
    P=sy.Symbol('actual_Pstar',positive=True);Rp=sy.Symbol('actual_Rp',positive=True);delta=sy.Symbol('delta',positive=True)
    base=dict(context['bindings']);base.update({owner.inlet.geometry['Pstar'].node:P,
        owner.inlet.geometry['Rp'].node:Rp,owner.pulse.delta.node:delta})
    rows=[];flat=[];pressure=[];radii=[];endpoint_lemmas=[];mass_ode=[];recovery_rows=[]
    field_derived=[];mass_derived=[];recovery_derived=[]
    derivative_rec=current.parent.linear.read(current.radial.RECEIPT)
    assert derivative_rec['all_passed'] and derivative_rec['source_family']==owner.family
    ordinary_source=derivative_rec['checks']['generic_recovery_FTC']['ordinary_jet_algebra_source_bindings']
    for pin in ordinary_source.values():assert owner.hashes[pin['module']]==pin['module_sha256']==current.sha(pin['module'])
    weights=[A,[A[j]+mu*D[j] for j in range(2)]]
    definitions=(('main_gap','main',sy.Integer(11),'gap',sy.Integer(11)),
        ('gap_end','gap',13-4*mu,'end',sy.Integer(-4)))
    for seam,left_name,left_x,right_name,right_x in definitions:
        algs=[];packets=[]
        for name,value in ((left_name,left_x),(right_name,right_x)):
            data=owner.pulse.charts[name];bindings=dict(base);bindings[data['coordinate'].node]=value
            if name=='main':
                for n,kernels in enumerate(owner.history_kernels['main']):
                    bindings[kernels[0].node]=sy.exp(context['lams'][n]*(13-value)/mu+p)*context['pis'][n]
            elif name=='end':
                for packets0 in owner.history_kernels['end']:
                    for ref in packets0:
                        n,column=context['end_columns'][ref.node]
                        bindings[ref.node]=sy.exp(-context['lams'][n]*value)*weights[n][column]
            # The first physical mean is the same source history times
            # the actual swirl. No mean or Q output has been bound yet.
            alg=KernelAlgebra(owner.pulse,bindings);algs.append(alg);packets.append(data)
        left,right=algs
        for key in ('Utheta','U','U_y','U_yy','m','m_y','m_yy'):
            for n,(lref,rref) in enumerate(zip(current.outer.rows(packets[0]['fields'][key]),current.outer.rows(packets[1]['fields'][key]))):
                if n>0:
                    # Both are genuine ordinary rows of the same just-proved
                    # C0 function, built from admitted C5 inputs and the
                    # pinned ordinary product/quotient source algebra.
                    field_derived.append(seam+'_'+key+'_Z'+str(n));continue
                zero(left.at(lref)-right.at(rref),(seam,'actual leading endpoint',key,n))
                rows.append(seam+'_'+key+'_Z'+str(n))
                if key.startswith('U') and key!='Utheta':
                    zero(left.at(lref),(seam,'actual flat axial',key,n));flat.append(seam+'_'+key+'_Z'+str(n))
        print('ACTUAL_INTERNAL_WHOLE_FIELD_IDENTITIES',seam,flush=True)
        Rleft=left.at(packets[0]['radius']);Rright=right.at(packets[1]['radius'])
        zero(Rleft-Rright,(seam,'actual physical radius'));radii.append(seam)
        for side,alg,data in (('left',left,packets[0]),('right',right,packets[1])):
            for n,(m,my,myy) in enumerate(zip(*(current.outer.rows(data['fields'][key]) for key in ('m','m_y','m_yy')))):
                if n>0:
                    mass_derived.extend((seam+'_'+side+'_m_y_plus_m_Z'+str(n),seam+'_'+side+'_m_yy_minus_m_Z'+str(n)));continue
                zero(alg.at(my)+alg.at(m),(seam,side,n,'actual flat mass first ODE'))
                zero(alg.at(myy)-alg.at(m),(seam,side,n,'actual flat mass second ODE'))
                mass_ode.extend((seam+'_'+side+'_m_y_plus_m_Z'+str(n),seam+'_'+side+'_m_yy_minus_m_Z'+str(n)))
        # Bind only the actual mean/U rows just compared as functions.
        # This retains the nonlinear Q recovery while avoiding expanding
        # the same selected-kernel sum inside every rational derivative.
        H=sy.Function('actual_'+seam+'_physical_mean')(Z)
        for alg,data in zip(algs,packets):
            for key,body in (('U',sy.Integer(0)),('U_y',sy.Integer(0)),('U_yy',sy.Integer(0)),
                ('m',H),('m_y',-H),('m_yy',H)):
                bind(alg.bindings,data['fields'][key],body,Z)
            alg.cache.clear()
        for key in ('Q','Q_y','Q_yy'):
            for n,(lref,rref) in enumerate(zip(current.outer.rows(packets[0]['fields'][key]),current.outer.rows(packets[1]['fields'][key]))):
                if n>0:
                    recovery_derived.append(seam+'_'+key+'_Z'+str(n));continue
                zero(left.at(lref)-right.at(rref),(seam,'actual leading Q recovery',key,n));recovery_rows.append(seam+'_'+key+'_Z'+str(n))
        # Q's exact endpoint recovery is now a proved source lemma.
        qbody=(-(1-delta)*Z*H-(1-Z**2)*sy.diff(H,Z))/(1-delta*Z**2)
        abstract_Q=sy.Function('actual_'+seam+'_recovered_Q')(Z)
        for alg,data in zip(algs,packets):
            for key,body in (('Q',qbody),('Q_y',-qbody),('Q_yy',qbody)):
                zero(alg.at(data['fields'][key].value)-body,(seam,key,'actual flat-gap whole-function recovery lemma'))
            for key in ('Q','Q_y','Q_yy'):
                bind(alg.bindings,data['fields'][key],-abstract_Q if key=='Q_y' else abstract_Q,Z)
            alg.cache.clear()
        for key in ('P1_y','N1p','m15_y'):
            for n,(lref,rref) in enumerate(zip(current.outer.rows(packets[0]['density'][key]),current.outer.rows(packets[1]['density'][key]))):
                zero(left.at(lref)-right.at(rref),(seam,'actual radial pressure source',key,n));pressure.append(seam+'_'+key+'_Z'+str(n))
        endpoint_lemmas.append(dict(seam=seam,actual_nonzero_gap_mean_retained=True,
            actual_complete_forward_and_backward_kernel_sources_consumed=True,
            mean_inputs_bound_only_after_actual_source_comparison=True,Q_inputs_bound_only_after_recovery_proof=True))
    return dict(passed=True,actual_leading_field_and_Q_internal_rows=rows,actual_axial_flat_rows=flat,
        higher_leading_field_rows_by_pinned_ordinary_derivative_commutation=field_derived,
        actual_flat_mass_ODE_endpoint_rows=mass_ode,
        higher_flat_mass_ODE_rows_by_pinned_ordinary_derivative_commutation=mass_derived,
        actual_Q_Qy_Qyy_internal_rows_by_proved_input_source_composition=recovery_rows,
        higher_Q_Qy_Qyy_rows_by_pinned_ordinary_derivative_commutation=recovery_derived,
        ordinary_jet_algebra_source_bindings=ordinary_source,
        accepted_actual_radial_source_receipt=current.radial.RECEIPT,
        accepted_actual_radial_source_receipt_sha256=current.sha(current.radial.RECEIPT),
        actual_radial_pressure_and_weighted_density_C3_internal_rows=pressure,
        actual_internal_radius_identities=radii,actual_endpoint_source_lemmas=endpoint_lemmas,
        original_distinct_main_gap_and_gap_end_coordinates=True,
        leading_recovery_one_minus_delta=True,actual_pressure_outputs_not_bound=True,
        pressure_rows_by_proved_input_source_and_pinned_map_composition=True,
        no_second_full_ancestor_pressure_expansion_claimed=True)


def own_second_radial_source(owner):
    owner.assert_graph();g=owner.g;Z=sy.Symbol('Z',real=True);mu=sy.Symbol('mu',positive=True)
    delta=sy.Symbol('delta',positive=True);Rp=sy.Symbol('actual_Rp',positive=True);M=sy.Function('actual_retained_own_mass')(Z)
    checked=[]
    for name,data in owner.pulse.charts.items():
        x=sy.Symbol(g.nodes[data['coordinate'].node]['name'],real=True);jac=1 if name=='end' else mu
        bindings={g.symbol('Z').node:Z,owner.pulse.mu.node:mu,owner.pulse.delta.node:delta,
            owner.inlet.geometry['Rp'].node:Rp,data['coordinate'].node:x}
        fields=owner.joins.current_own[name]['fields']
        assert current.outer.rows(fields['M1'])==current.outer.rows(owner.joins.Rp_cumulative['m11'])
        for ref in current.outer.rows(fields['M1']):
            assert not outer_check.free_query_dependency(owner.joins.outer,ref.node,data['coordinate'].node)
        bind(bindings,fields['M1'],M,Z)
        alg=outer_check.SourceAlgebra(g,bindings)
        for key,original in (('Uz1_yy','Uz1'),('Utheta1_yy','Utheta1'),('F1_yy','F1'),
            ('M1_yy','M1'),('mean1_yy','mean1'),('K1_yy','K1'),('Q1_yy','Q1'),('V1_yy','V1')):
            for n,(ref,bodyref) in enumerate(zip(current.outer.rows(owner.own_second_radial[name][key]),current.outer.rows(fields[original]))):
                zero(alg.at(ref)-jac**2*sy.diff(alg.at(bodyref),x,2),(name,key,n,'actual own second radial source'))
                checked.append(name+'_'+key+'_Z'+str(n))
    return dict(passed=True,actual_own_second_radial_source_rows=checked,
        actual_saved_mass_equals_same_query_independent_Rp_m11_before_binding=True,
        retained_nonzero_mass_arbitrary_function=True,original_query_Jacobian=True,
        own_one_plus_delta_recovery_retained=True,
        no_pressure_second_radial_claim=True,no_extra_Q1_Z3_claim=True)


def internal_primitive_composition(owner,leading_proof):
    owner.assert_graph();assert leading_proof['passed'];g=owner.g
    a=current.outer.Algebra(g,3);Z=sy.Symbol('Z',real=True);delta=sy.Symbol('delta',positive=True)
    integral_rows=[];tensor_rows=[];pressure_rows=[]
    for seam,left_name,right_name in (('main_gap','main','gap'),('gap_end','gap','end')):
        left=owner.pulse.charts[left_name];right=owner.pulse.charts[right_name]
        assert right['preceding_complete_charts']==left['preceding_complete_charts']+[left_name]
        bindings={g.symbol('Z').node:Z,owner.pulse.delta.node:delta}
        # Actual complete preceding histories remain common source inputs.
        for name,packet in owner.pulse.complete_integrals.items():
            for key,jet in packet.items():
                for ref in current.outer.rows(jet):
                    row=g.nodes[ref.node]
                    bindings.setdefault(ref.node,sy.Integer(0) if row['integrand']==0
                        else sy.Symbol('actual_complete_integral_'+str(ref.node),real=True))
        for data,is_left in ((left,True),(right,False)):
            for key,jet in data['partial_integrals'].items():
                expected=a.add(*(owner.pulse.complete_integrals[name][key] for name in data['preceding_complete_charts']),jet)
                assert current.outer.rows(expected)==current.outer.rows(data['total_integrals'][key])
                for n,ref in enumerate(current.outer.rows(jet)):
                    row=g.nodes[ref.node];full=current.outer.rows(owner.pulse.complete_integrals[left_name][key])[n]
                    assert row['operation']=='definite_integral' and row['variable']==g.nodes[data['coordinate'].node]['name']
                    assert row['lower']==data['domain'][0].node and row['upper']==data['coordinate'].node
                    assert row['endpoint_independent_of_Z'] and row['ordinary_slow_Z_derivative_order']==n
                    if is_left:
                        complete=g.nodes[full.node]
                        assert complete['variable']==row['variable'] and complete['lower']==row['lower']
                        assert complete['upper']==data['domain'][1].node and complete['integrand']==row['integrand']
                        assert complete['ordinary_slow_Z_derivative_order']==n
                    bindings[ref.node]=bindings[full.node] if is_left else sy.Integer(0)
                    integral_rows.append(seam+'_'+('left_full' if is_left else 'right_empty')+'_'+key+'_Z'+str(n))
        Q=sy.Function('actual_'+seam+'_recovered_Q')(Z);R=sy.Symbol('actual_'+seam+'_matched_radius',positive=True)
        for data in (left,right):
            bindings[data['radius'].node]=R
            for key,body in (('U',sy.Integer(0)),('Q',Q),('Q_y',-Q)):
                bind(bindings,data['fields'][key],body,Z)
        for jet in owner.pulse.global_lower_tensors.values():
            for ref in current.outer.rows(jet):
                if ref.node not in (0,1):bindings.setdefault(ref.node,sy.Symbol('retained_global_lower_'+str(ref.node),real=True))
        alg=outer_check.SourceAlgebra(g,bindings)
        for key in left['upper_tensors']:
            for n,(lref,rref) in enumerate(zip(current.outer.rows(left['upper_tensors'][key]),current.outer.rows(right['upper_tensors'][key]))):
                zero(alg.at(lref)-alg.at(rref),(seam,key,n,'actual upper source tensor'));tensor_rows.append(seam+'_'+key+'_Z'+str(n))
        for key in ('m13','m15'):
            for n,(lref,rref) in enumerate(zip(current.outer.rows(left['signed_pressure_increment'][key]),current.outer.rows(right['signed_pressure_increment'][key]))):
                zero(alg.at(lref)-alg.at(rref),(seam,key,n,'actual retained pressure primitive'));pressure_rows.append(seam+'_'+key+'_Z'+str(n))
            for data,name in ((left,left_name),(right,right_name)):
                expected=a.add(owner.joins.Rp_cumulative[key],data['signed_pressure_increment'][key])
                assert current.outer.rows(expected)==current.outer.rows(owner.joins.current_own[name]['cumulative'][key])
    return dict(passed=True,actual_left_complete_right_empty_integral_endpoint_rows=integral_rows,
        actual_common_upper_tensor_C3_endpoint_rows=tensor_rows,
        actual_retained_m13_m15_C3_internal_pressure_primitive_rows=pressure_rows,
        actual_full_preceding_histories_and_finite_global_lower_tensors_retained=True,
        actual_common_Rp_cumulative_anchor_not_reset=True,
        integral_endpoints_identified_by_exact_same_body_and_named_binder=True,
        actual_pressure_outputs_not_bound=True,
        actual_own_P1_value_and_first_radial_derivative_internal_C3_joins=True,
        actual_first_radial_derivative_uses_accepted_pressure_FTC=True)


def run(owner):
    began=time.monotonic();owner.assert_graph()
    kernel,context=fixed_kernel_source(owner)
    print('ACTUAL_FIXED_KERNEL_TRANSPORT_PASSED',flush=True)
    selected,context=selected_linear_source(owner,context)
    print('ACTUAL_SELECTED_LINEAR_C5_SOURCE_PASSED',flush=True)
    leading=leading_internal_sources(owner,context)
    print('ACTUAL_INTERNAL_RADIAL_PRESSURE_SOURCE_JOINS_PASSED',flush=True)
    checks=dict(fixed_kernel_source=kernel,selected_linear_source=selected,leading_internal_sources=leading,
        internal_primitive_composition=internal_primitive_composition(owner,leading),
        own_second_radial_source=own_second_radial_source(owner))
    return save(owner,checks,time.monotonic()-began)


def save(owner,checks,execution_seconds):
    owner.assert_graph()
    assert set(checks)=={'fixed_kernel_source','selected_linear_source','leading_internal_sources',
        'internal_primitive_composition','own_second_radial_source'}
    assert all(proof['passed'] for proof in checks.values())
    raw=current.parent.linear.read(current.NAME)
    assert raw['pulse_internal_radial_source_join_bindings']==owner.snapshot
    assert raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    pins=dict(owner.hashes);pins[current.NAME]=current.sha(current.NAME);pins[Path(__file__).name]=current.sha(Path(__file__).name)
    receipt=dict(all_passed=True,source_family=owner.family,input_hashes=pins,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        actual_internal_radial_pressure_and_density_rows=len(checks['leading_internal_sources']['actual_radial_pressure_and_weighted_density_C3_internal_rows']),
        actual_internal_pressure_primitive_rows=len(checks['internal_primitive_composition']['actual_retained_m13_m15_C3_internal_pressure_primitive_rows']),
        actual_own_second_radial_rows=len(checks['own_second_radial_source']['actual_own_second_radial_source_rows']),
        exact_graph_nodes=len(raw['exact_graph_nodes']),execution_seconds=execution_seconds,
        numerical_oracle_available=False,combined_outer_seam_accepted=False,global_n1=False)
    (current.HERE/current.RECEIPT).write_text(__import__('json').dumps(receipt,indent=2)+'\n',encoding='utf-8',newline='\n')
    owner.load_acceptance()
    print('ACTUAL_INTERNAL_RADIAL_SOURCE_JOINS_ACCEPTED',receipt['exact_graph_nodes'],
        'seconds',receipt['execution_seconds'],flush=True)
    return receipt


if __name__=='__main__':
    import sys
    sys.setrecursionlimit(10000)
    run(current.CurrentOriginalPulseInternalRadialSourceJoins(require_checked=False))
