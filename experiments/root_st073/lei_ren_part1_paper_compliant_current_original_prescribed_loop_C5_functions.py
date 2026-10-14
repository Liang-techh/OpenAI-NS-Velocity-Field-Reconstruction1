"""True prescribed-source Section 11 C5 loops and finite-N transport to Rm.

Ordinary axial rows are genuine function derivatives. The phase inverse is
differentiated by its defining implicit equation. Original N/eta/d_star,
the microscopic collar, and every incoming correction history are retained.
Quantitative high-order admission and the remaining seven charts stay open.
"""
import copy
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
import time

import lei_ren_part1_paper_compliant_current_original_prescribed_R100_Rm_C6_functions as previous
import lei_ren_part1_paper_compliant_current_generic_loop_function_sources as original

source,parent,controls=previous.source,previous.parent,previous.controls
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_prescribed_loop_C5_functions.json.gz'
RECEIPT=PREFIX+'current_original_prescribed_loop_C5_functions_check.json'
GATES=('current_original_true_prescribed_ten_chart_signed_loop_C5_functions_installed',
       'current_original_true_prescribed_ten_chart_finite_N_density_C5_functions_installed',
       'current_original_true_prescribed_rminus_to_Rm_correction_history_C5_functions_installed')
OPEN=tuple(dict.fromkeys(previous.OPEN+('current_original_prescribed_loop_C5_numeric_containment_installed',
    'current_original_all_17_chart_finite_N_C5_functions_installed',
    'current_original_actual_modified_bridge_C6_functions_installed',
    'current_original_actual_modified_R100_Rm_C6_functions_installed')))
FunctionJet,rows,encode=previous.FunctionJet,previous.rows,previous.encode
CHARTS=('first_micro','second_micro','frozen_macro','first_switch','second_switch','post_power',
        'long_reshape','reference','restoration','postrestore')
RATES=dict(m=Fraction(1),h=Fraction(3,2),k=Fraction(3,2),e=Fraction(1),p=Fraction(0))


class LoopAlgebra(previous.Algebra):
    def sqrt(self,jet):
        q=rows(self.trim(jet));g=self.g;out=[g.unary('positive_sqrt',q[0])]
        for n in range(1,self.order+1):
            rest=g.add(*(g.mul(g.constant(math.comb(n,j)),out[j],out[n-j]) for j in range(1,n)))
            out.append(g.quotient(g.sub(q[n],rest),g.mul(g.constant(2),out[0]),
                'original active positive square root; lazy flat extension outside active domain'))
        return FunctionJet(tuple(out))

    def lazy(self,jet,predicate,Delta,eta,coordinate=None,collar=None):
        g=self.g
        return FunctionJet(tuple(g.node('original_lazy_flat_branch',active_body=q.node,
            Delta=Delta.node,eta=eta.node,flat_value=g.zero.node,flat_predicate=predicate.node,
            original_collar_coordinate=coordinate.node if collar is not None else None,
            original_collar_upper=collar.node if collar is not None else None,
            active_body_not_evaluated_on_flat_branch=True,ordinary_slow_Z_derivative_order=n,
            same_original_flat_extension_through_order5=True) for n,q in enumerate(rows(self.trim(jet)))))

    def sigma(self,x):
        """Original value plus logistic ODE derivatives, with exact flat ends."""
        g=self.g;one=self.fixed(1)
        value=previous.sigma(g,rows(x)[0])
        predicate=g.node('logical_or',arguments=[
            g.node('exact_real_comparison',operator='<=',left=rows(x)[0].node,right=g.zero.node).node,
            g.node('exact_real_comparison',operator='>=',left=rows(x)[0].node,right=g.one.node).node])
        inverse=self.div(one,x,'interior sigma argument x>0')
        complement=self.div(one,self.sub(one,x),'interior sigma argument x<1')
        odds=self.sub(self.mul(complement,complement),self.mul(inverse,inverse))
        out=[value]
        for n in range(self.order):
            # d sigma/dZ = sigma*(1-sigma)*d(log odds)/dZ.
            a=LoopAlgebra(g,n)
            p=FunctionJet(tuple(out));v=a.mul(p,a.sub(a.fixed(1),p))
            body=g.add(*(g.mul(g.constant(math.comb(n,j)),rows(v)[j],rows(odds)[n-j+1]) for j in range(n+1)))
            out.append(g.node('original_lazy_flat_branch',active_body=body.node,
                Delta=rows(x)[0].node,eta=g.one.node,flat_value=g.zero.node,
                flat_predicate=predicate.node,active_body_not_evaluated_on_flat_branch=True,
                ordinary_slow_Z_derivative_order=n+1,original_sigma_logistic_ODE=True))
        return FunctionJet(tuple(out))


class MixedTaylor:
    """Normalized Z/angle coefficients of total degree <= order."""
    def __init__(self,g,order=5):
        self.g,self.order=g,order
        self.indices=tuple((k,j) for k in range(order+1) for j in range(order+1-k))
    def constant(self,ref):
        if type(ref) is not source.FunctionRef:ref=self.g.constant(ref)
        self.g.ids([ref]);return {(0,0):ref}
    def axial(self,jet):
        if len(rows(jet))<self.order+1:raise ValueError('Genuine axial rows required')
        return {(k,0):self.g.mul(self.g.constant(Fraction(1,math.factorial(k))),rows(jet)[k]) for k in range(self.order+1)}
    def add(self,*values):
        return {ij:self.g.add(*(v.get(ij,self.g.zero) for v in values)) for ij in self.indices}
    def scale(self,value,ref):
        if type(ref) is not source.FunctionRef:ref=self.g.constant(ref)
        return {ij:self.g.mul(ref,v) for ij,v in value.items()}
    def mul(self,left,right):
        g=self.g
        return {(k,j):g.add(*(g.mul(left.get((a,b),g.zero),right.get((k-a,j-b),g.zero))
            for a in range(k+1) for b in range(j+1))) for k,j in self.indices}
    def div(self,left,right):
        g=self.g;out={}
        for degree in range(self.order+1):
            for k in range(degree+1):
                j=degree-k
                rest=g.add(*(g.mul(right.get((a,b),g.zero),out.get((k-a,j-b),g.zero))
                    for a in range(k+1) for b in range(j+1) if a+b>0))
                out[k,j]=g.quotient(g.sub(left.get((k,j),g.zero),rest),right[0,0],
                    '1-2*r*cos(psi)+r^2=(cos(psi)-r)^2+sin(psi)^2>0 since |r|<1')
        return out
    def trig(self,name,psi):
        g=self.g;cos=g.unary('cos',psi);sin=g.unary('sin',psi)
        cycle=(cos,g.neg(sin),g.neg(cos),sin) if name=='cos' else (sin,cos,g.neg(sin),g.neg(cos))
        return {(0,j):g.mul(g.constant(Fraction(1,math.factorial(j))),cycle[j%4]) for j in range(self.order+1)}


def at_inverse(g,value,inverse,angle):
    return g.node('substitute_original_inverse_angle',body=value.node,inverse_angle=inverse.node,
        angle_variable=g.nodes[angle.node]['name'])


def series_powers(g,delta,order):
    """Coefficients of powers of a series with zero constant term."""
    powers=[[g.one]+[g.zero]*order]
    for j in range(1,order+1):
        before=powers[-1]
        powers.append([g.add(*(g.mul(before[k],delta[n-k]) for k in range(n+1))) for n in range(order+1)])
    return powers


def implicit_inverse(g,partials,inverse,angle,order=5):
    """Solve each normalized coefficient of Phi(Z+h,psi(Z+h))=phi."""
    at={ij:at_inverse(g,value,inverse,angle) for ij,value in partials.items()}
    lam=at[0,1];delta=[g.zero]*(order+1)
    for n in range(1,order+1):
        powers=series_powers(g,delta,n)
        known=g.add(*(g.mul(at[k,j],powers[j][n-k])
            for k in range(n+1) for j in range(n-k+1) if (k,j)!=(0,0)))
        delta[n]=g.neg(g.quotient(known,lam,'original phase inverse lambda=K*(1+t^2)>0'))
    ordinary=FunctionJet((inverse,)+tuple(g.mul(g.constant(math.factorial(n)),delta[n]) for n in range(1,order+1)))
    return ordinary,delta,at


def compose_inverse(g,partials,inverse,angle,delta,order=5):
    at={ij:at_inverse(g,value,inverse,angle) for ij,value in partials.items()}
    powers=series_powers(g,delta,order)
    return FunctionJet(tuple(g.mul(g.constant(math.factorial(n)),
        g.add(*(g.mul(at[k,j],powers[j][n-k]) for k in range(n+1) for j in range(n-k+1))))
        for n in range(order+1)))


def loop_functions(g,roots,phase,eta,dstar,chart,coordinate,sc):
    """Unchanged signed loop bodies with all implicit C5 cross terms."""
    a=LoopAlgebra(g,5);one=a.fixed(1);c=g.constant
    shear=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)
    collar=g.mul(c('3/4'),sc) if chart=='first_micro' else None
    predicate=shear
    if collar is not None:
        predicate=g.node('logical_or',arguments=[shear.node,
            g.node('exact_real_comparison',operator='<=',left=coordinate.node,right=collar.node).node])
    lazy=lambda jet:a.lazy(jet,predicate,roots['Delta'].value,eta,coordinate,collar)
    gamma=a.sub(a.fixed(g.mul(c(2),eta)),roots['Delta'])
    cut=a.sub(one,a.scale(roots['Delta'],g.quotient(g.one,eta,'same positive canonical eta')))
    q=lazy(a.mul(a.sigma(cut),a.sqrt(a.div(gamma,a.scale(roots['a'],2),'same original active a>0'))))
    u=a.scale(a.mul(roots['p2'],q),g.quotient(g.one,dstar,'same positive canonical d_star'))
    hinv=a.div(one,a.sqrt(a.add(one,a.mul(u,u))),'positive sqrt(1+u^2)')
    r=a.mul(u,hinv);alpha=a.scale(a.mul(q,hinv),2)
    angle=g.symbol('true_loop_C5_angle_'+chart);m=MixedTaylor(g,5)
    rr=m.axial(r);cc=m.trig('cos',angle)
    numerator=m.add(cc,m.scale(rr,-1))
    den=m.add(m.constant(1),m.scale(m.mul(rr,cc),-2),m.mul(rr,rr))
    t=m.add(m.axial(roots['t0']),m.mul(m.axial(alpha),m.div(numerator,den)))
    t2=m.mul(t,t)
    ratio=a.add(one,a.mul(roots['t0'],roots['t0']),a.scale(a.mul(q,q),2))
    pi=g.node('mathematical_pi');twopi=g.mul(c(2),pi)
    K=a.div(one,a.scale(ratio,twopi),'same positive 2*pi*(1+t0^2+2*q^2)')
    integrate=lambda ref:controls.integral(g,ref,g.nodes[angle.node]['name'],g.zero,angle,
        measure='original angle dpsi with Z and native radius held fixed',
        true_loop_C5_original_integrand=True,angle_endpoints_independent_of_Z=True)
    T1={(k,0):integrate(t[k,0]) for k in range(6)}
    T2=FunctionJet(tuple(g.mul(c(math.factorial(k)),integrate(t2[k,0])) for k in range(6)))
    P=a.add(a.fixed(angle),T2);Phi=a.mul(K,P)
    lam=m.mul(m.axial(K),m.add(m.constant(1),t2))
    partials={(k,0):g.mul(c(Fraction(1,math.factorial(k))),rows(Phi)[k]) for k in range(6)}
    for k,j in m.indices:
        if j:
            partials[k,j]=g.mul(c(Fraction(1,j)),lam[k,j-1])
            T1[k,j]=g.mul(c(Fraction(1,j)),t[k,j-1])
    inverse=g.node('original_monotone_phase_inverse',phase_function=Phi.value.node,
        angle_variable=g.nodes[angle.node]['name'],target_phase=phase.node,
        angle_lower=g.zero.node,angle_upper=twopi.node,phase_endpoints=[0,1],
        derivative_positive_certificate='K*(1+t^2)>0; same original a/q-flat construction',
        flat_inverse='2*pi*phase',Delta=roots['Delta'].value.node,eta=eta.node,
        branch_inverse_is_original_not_support_cap=True,
        true_prescribed_source_C5_implicit_inverse=True)
    psi,delta,at=implicit_inverse(g,partials,inverse,angle)
    T1hat=compose_inverse(g,T1,inverse,angle,delta)
    that=compose_inverse(g,t,inverse,angle,delta)
    lamhat=compose_inverse(g,lam,inverse,angle,delta)
    chi=a.sub(a.fixed(phase),a.scale(psi,g.quotient(g.one,twopi,'positive 2*pi')))
    M=a.sub(a.neg(a.scale(a.mul(roots['a'],T1hat),g.quotient(g.one,twopi,'positive 2*pi'))),
        a.scale(roots['b'],phase))
    A=lazy(a.scale(a.mul(roots['a'],chi),'1/2'))
    B=lazy(a.scale(a.mul(roots['E'],M),'1/2'))
    invlambda=a.div(one,lamhat,'original lambda at true inverse positive')
    Aphi=lazy(a.scale(a.mul(roots['a'],a.sub(one,a.scale(invlambda,g.quotient(g.one,twopi,'positive 2*pi')))),'1/2'))
    Bphi=lazy(a.scale(a.mul(roots['E'],a.sub(a.neg(a.scale(a.mul(roots['a'],a.mul(that,invlambda)),
        g.quotient(g.one,twopi,'positive 2*pi'))),roots['b'])),'1/2'))
    return dict(q=q,inverse=psi,A=A,B_over_Pstar=B,A_phi=Aphi,B_over_Pstar_phi=Bphi,
        angle=angle,Phi=Phi,phase_lambda_at_inverse=lamhat,T1_at_inverse=T1hat,
        direction_at_inverse=that,phase_partial_coefficients=partials,
        inverse_normalized_increment=delta,flat_predicate=predicate,original_collar_upper=collar,
        positive_denominator_exact_recipe='D=1-2*r*cos(psi)+r^2; use correlated original Poisson evaluation for numerical enclosures',
        no_quantitative_C5_inverse_conditioning_claimed=True)


def generic_roots(g,chart,Pstar):
    a=LoopAlgebra(g,5);field=chart['fields'];radial=chart['radial'];generic=chart['generic']
    E,V=(a.trim(field[key]) for key in ('E','V'))
    C=a.sub(E,a.scale(a.trim(radial['E_y']),2));B=a.scale(a.trim(radial['V_y']),2)
    divide=lambda jet:a.div(jet,E,'same strictly positive actual prescribed swirl E')
    aa,bb=divide(C),divide(B)
    p1=divide(a.scale(a.add(generic['theta_linear'],a.scale(generic['theta_quadratic'],Pstar)),chart['radius']))
    # Preserve the already proved full radius/Pstar p2 function by identity.
    p2=generic['p2'];t0=a.neg(a.div(B,C,'same original a>0 and C=a*E'))
    kap=a.add(aa,a.div(a.mul(bb,bb),aa,'same original positive a'))
    return dict(E=E,V=V,a=aa,b=bb,p1=p1,p2=p2,t0=t0,Delta=a.sub(kap,a.fixed(2)),kappa=kap)


def finite_N_functions(g,roots,loop,N):
    a=LoopAlgebra(g,5);invN=g.quotient(g.one,N,'same exact positive integer N on every source chart')
    x=a.scale(loop['A'],invN);ex=a.exp(x)
    # Value uses expm1; all higher ordinary rows use the full exp recurrence.
    exm1=FunctionJet((g.unary('expm1',x.value),)+rows(ex)[1:])
    r=a.mul(roots['E'],exm1);s=a.scale(loop['B_over_Pstar'],invN)
    E,V=roots['E'],roots['V']
    inc=dict(m=s,h=r,k=a.add(a.mul(V,r),a.mul(E,s),a.mul(r,s)),
        e=a.sub(a.add(a.scale(a.mul(V,s),2),a.mul(s,s)),a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2'))),
        p=a.add(a.mul(E,r),a.scale(a.mul(r,r),'1/2')))
    return dict(E_N=a.mul(E,ex),V_N=a.add(V,s),delta_E=r,delta_V=s,A_over_N=x,
        five_signed_density_increments=inc)


def transport(g,densities,incoming,coordinate,lower,jacobian,label):
    a=LoopAlgebra(g,5);t=g.symbol('true_loop_C5_transport_'+label)
    width=g.mul(jacobian,g.sub(coordinate,lower));out={};parts={};memory={}
    for name,rate in RATES.items():
        decay=g.unary('exp',g.mul(g.constant(-rate),width)) if rate else g.one
        kernel=g.unary('exp',g.mul(g.constant(-rate),jacobian,g.sub(coordinate,t))) if rate else g.one
        body=a.scale(a.at(densities[name],coordinate,t),g.mul(jacobian,kernel))
        contribution=a.integral(body,t,lower,coordinate,actual_finite_N_signed_correction=True,
            original_rate=str(rate),original_chart=label,physical_Jacobian_applied_once=True,
            zero_correction_inlet_only_at_original_rminus=True,predecessor_memory_never_reset=True)
        out[name]=a.add(a.scale(incoming[name],decay),contribution)
        parts[name]=contribution;memory[name]=decay
    return out,parts,memory


class CurrentOriginalPrescribedLoopC5Functions:
    def __init__(self,background,require_checked=True):
        if type(background) is not previous.CurrentOriginalPrescribedR100RmC6Functions or not background.acceptance_loaded:
            raise TypeError('Checked true prescribed background C6 through Rm required')
        background.assert_graph();self.background=background;self.g=g=background.g
        self.family=copy.deepcopy(background.family);self.hashes=dict(background.hashes)
        for name,digest in parent.linear.read(previous.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for name in (previous.NAME,previous.RECEIPT,Path(original.__file__).name):self._pin(name,sha(name))
        self.bindings=dict(original_loop=parent.active.current.ast_binding(previous.connection.original_loop),
            original_generic_roots=previous.previous.inlet.generic.source_bindings(),
            original_density_graph=parent.active.current.ast_binding(original.build))
        self.N=self.unique('exact_positive_integer_expression',None)
        self.eta=self.unique('current_original_source_parameter','eta')
        self.dstar=self.unique('current_original_source_parameter','d_star')
        self.sc=self.unique('current_original_source_parameter','sc')
        if not g.nodes[self.N.node].get('one_common_N_for_all_charts_and_control_denominators'):
            raise ValueError('Same canonical common repair-only N source required')
        self.Pstar=background.Pstar;self.P0=background.P0
        self.phase_origin=g.mul(background.bridge.hb,self.sc,g.constant('1/2'))
        bridge=background.bridge;self.source_charts={
            'first_micro':bridge.micro,'second_micro':bridge.micro,'frozen_macro':bridge.macro,
            'first_switch':background.charts['first_switch'],'second_switch':background.charts['second_switch'],
            'post_power':background.charts['post_switch'],'long_reshape':background.charts['long_reshape'],
            'reference':background.charts['reference'],'restoration':background.charts['restoration'],
            'postrestore':background.charts['post_restore']}
        self.charts={};self.endpoints={};a=LoopAlgebra(g,5)
        incoming={key:a.fixed(0) for key in RATES};self.zero_inlet=incoming
        for name in CHARTS:
            chart=self.source_charts[name];coordinate=chart['phase']
            lower=g.mul(self.sc,g.constant('1/2')) if name=='first_micro' else g.one if name=='second_micro' else g.zero
            upper=g.constant(2) if name=='second_micro' else g.one
            offset=g.sub(g.unary('log',g.quotient(chart['radius'],bridge.Ra,'same positive actual R/Ra')),self.phase_origin)
            phase=g.unary('fractional_part',g.mul(self.N,offset))
            roots=generic_roots(g,chart,self.Pstar)
            loop=loop_functions(g,roots,phase,self.eta,self.dstar,name,coordinate,self.sc)
            finite=finite_N_functions(g,roots,loop,self.N)
            history,parts,memory=transport(g,finite['five_signed_density_increments'],incoming,
                coordinate,lower,chart['physical_measure'],name)
            fields={key:a.trim(jet) for key,jet in chart['fields'].items()}
            fields.update(E=finite['E_N'],V=finite['V_N'])
            for key in RATES:fields[key]=a.add(a.trim(chart['fields'][key]),history[key])
            fields['pressure']=a.add(fields['P0'],fields['p'])
            self.charts[name]=dict(source_chart=chart,coordinate=coordinate,domain=[lower,upper],
                radius=chart['radius'],physical_measure=chart['physical_measure'],global_offset=offset,
                actual_phase=phase,roots=roots,loop=loop,finite_N=finite,incoming=incoming,
                history=history,contributions=parts,memory=memory,modified_fields=fields,
                same_background_pressure_datum=True,total_Z_phase_derivative_exact_zero=True)
            self.endpoints[name]={key:a.at(jet,coordinate,upper) for key,jet in history.items()}
            incoming=self.endpoints[name]
        self.Rm_correction=incoming
        self.Rm_fields={key:a.at(jet,self.charts['postrestore']['coordinate'],g.one)
            for key,jet in self.charts['postrestore']['modified_fields'].items()}
        self.definition=dict(true_prescribed_background_source_not_legacy_recipe_derivative_padding=True,
            all_ten_signed_root_q_inverse_A_B_and_phase_primitive_rows_through_Z5=True,
            inverse_rows_from_exact_implicit_normalized_series_not_interval_endpoints=True,
            original_first_micro_start_sc_half_and_collar_three_sc_quarters=True,
            one_original_common_N_eta_dstar_and_global_log_R_rminus_phase=True,
            all_five_signed_finite_N_density_cross_terms_and_expm1_retained=True,
            complete_original_rminus_to_Rm_Duhamel_correction_histories_C5=True,
            correction_zero_inlet_only_at_rminus_background_nonzero_memories_unchanged=True,
            original_rate_zero_pressure_memory_retained=True,
            new_functions_preserve_original_fixed_leading_source_family=True,
            original_repair_only_N_not_global_high_derivative_frequency_admission=True,
            numeric_C5_containment_and_mixed_y_bounds_not_certified=True,
            remaining_actual_Rm_patch_and_six_outer_charts_not_installed=True,
            all17_compact_controls_global_n1_recursion_and_corrected_NS_open=True,
            shared_graph_unreferenced_suffix_allowed_after_exact_pinned_prefix=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.background,self.g);self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
            checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not receipt.get('all_passed') or not all(receipt.get(k) for k in GATES) or any(receipt.get(k) for k in OPEN)
                or receipt['source_family']!=self.family or pins!=self.hashes or raw['input_hashes']!=self.hashes
                or archive!=sha(NAME) or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['prescribed_loop_C5_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked true prescribed-loop C5 source required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed loop C5 source '+name)
            self.acceptance_loaded=True
    def unique(self,operation,name):
        found=[source.FunctionRef(self.g,i) for i,row in enumerate(self.g.nodes)
            if row['operation']==operation and (name is None or row.get('name')==name)]
        if len(found)!=1:raise ValueError('Unique original source parameter required '+str(name))
        if name is not None:
            row=self.g.nodes[found[0].node]
            if not row['Z_independent'] or row['source_family']!=self.family:
                raise ValueError('Same original Z-independent loop parameter required')
        return found[0]
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed loop source '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        # Tuple mixed-partial keys are explicitly encoded, not string-eval'd.
        charts={name:{key:value for key,value in chart.items() if key!='source_chart'} for name,chart in self.charts.items()}
        charts=copy.copy(charts)
        for name,chart in list(charts.items()):
            chart=dict(chart);loop=dict(chart['loop'])
            loop['phase_partial_coefficients']=[dict(Z=k,angle=j,node=value.node)
                for (k,j),value in loop['phase_partial_coefficients'].items()]
            chart['loop']=loop;charts[name]=chart
        return encode(dict(charts=charts,endpoints=self.endpoints,Rm_correction=self.Rm_correction,
            Rm_fields=self.Rm_fields,N=self.N,eta=self.eta,dstar=self.dstar,sc=self.sc,
            phase_origin=self.phase_origin,definition=self.definition,bindings=self.bindings))
    def assert_graph(self):
        self.background.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot
            or self.hashes!=self._hashes or self._owners!=(self.background,self.g)):
            raise ValueError('True loop C5 source, inverse or transported memory changed')
        expected=self.zero_inlet
        for name in CHARTS:
            chart=self.charts[name]
            if chart['source_chart'] is not self.source_charts[name] or chart['incoming'] is not expected:
                raise ValueError('True background or correction predecessor identity changed')
            expected=self.endpoints[name]


def run(background):
    began=time.monotonic();owner=CurrentOriginalPrescribedLoopC5Functions(background,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,
        prescribed_loop_C5_bindings=owner.snapshot,input_hashes=owner.hashes,
        **dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_PRESCRIBED_TEN_CHART_LOOP_AND_TRANSPORT_C5_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw


if __name__=='__main__':raise SystemExit('run(background) reuses the checked true prescribed Rm owner')
