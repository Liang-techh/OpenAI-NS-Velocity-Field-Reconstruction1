"""True core-to-collar histories and signed physical n1 pressure functions.

All integrals have actual bodies. Original binary sc and formal hb remain
exact and nonzero; neither the core comparison nor tiny caps define a field.
"""
import copy
import gzip
import json
from pathlib import Path
import time
from types import SimpleNamespace

import lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions as core
import lei_ren_part1_paper_compliant_current_original_n1_unchanged_pressure_prefix as prior
import lei_ren_part1_paper_compliant_current_original_n1_common_inner_selection as selection
import lei_ren_part1_paper_compliant_current_original_micro_finite_N as modification

parent,source,controls,c3=core.parent,core.source,core.controls,core.c3
own=core.pulse.own
HERE,PREFIX,sha=core.HERE,core.PREFIX,core.sha
NAME=PREFIX+'current_original_n1_lower_collar_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_lower_collar_functions_check.json'
GATE='current_original_n1_Rin_to_rminus_true_lower_collar_C1_functions_installed'
OPEN=core.OPEN


def at(g,jet,variable,value):
    return c3.C3Function(*(g.node('function_substitution',expression=q.node,
        variable=variable.node,value=value.node,Z_independent_substitution=True)
        for q in c3.rows(jet)))


def integral(g,jet,variable,lower,upper,**metadata):
    return c3.C3Function(*(controls.integral(g,q,g.nodes[variable.node]['name'],lower,upper,
        endpoint_independent_of_Z=True,ordinary_slow_Z_derivative_order=k,
        actual_complete_collar_history=True,numerical_cap_not_defining_integrand=True,
        **metadata) for k,q in enumerate(c3.rows(jet))))


def exp3(g,q):
    v=g.unary('exp',q.value);c=g.constant
    return c3.C3Function(v,g.mul(v,q.Z),g.mul(v,g.add(q.ZZ,g.mul(q.Z,q.Z))),
        g.mul(v,g.add(q.ZZZ,g.mul(c(3),q.Z,q.ZZ),g.mul(q.Z,q.Z,q.Z))))


def exact_binary_scalar(g,bits):
    """man*2**exponent without allocating an astronomically large integer."""
    sign,man,exponent,_=bits;c=g.constant
    return g.mul(c(-man if sign else man),g.unary('exp',
        g.mul(c(exponent),g.unary('log',c(2)))))


def true_core_rows(owner):
    g=owner.g;rows={}
    for quantity in ('Phi','Psi','mean_Psi'):
        template=copy.deepcopy(g.nodes[owner.components[quantity].value.node])
        operation=template.pop('operation');jets=[]
        for i in range(3):
            jets.append(c3.C3Function(*(g.node(operation,**{
                **template,'radial_order':i,'ordinary_slow_Z_derivative_order':k}) for k in range(4))))
        rows[quantity]=jets
    return rows


def collar_fields(g,rho_variable,core_rows,core_functions,parameters,sigma=None,sigma_prime=None):
    """Exact prescribed-shear continuation, including own nonlocal histories."""
    a=c3.C3Algebra(g);c=g.constant;s=g.symbol('original_lower_collar_s')
    hb,eps=parameters['hb'],parameters['epsilon'];de=parameters['delta']
    sigma=sigma or (lambda x:own.jets.complete.Gamma.sigma(g,x))
    sigma_prime=sigma_prime or (lambda x:g.node('analytic_unary',
        name='original_flat_sigma_prime',argument=x.node,
        derivative_source=own.c2.current.ast_binding(own.c2.flat_source.sigma_jets),
        ordinary_derivative='same original sigma Taylor coefficient[1]',
        exact_flat_source_derivative_not_a_bound=True))
    rho=lambda x:g.mul(c(4),g.unary('exp',g.mul(hb,x)))
    chi=lambda x:g.add(g.one,g.mul(g.sub(hb,g.one),sigma(x)))
    rs=rho(s);R=g.mul(eps,rs);jacobian=g.mul(hb,rs)
    at_core=lambda q,x:at(g,q,rho_variable,rho(x))
    fc=[at_core(q,s) for q in core_rows['Phi']]
    pc=[at_core(q,s) for q in core_rows['Psi']]
    A=a.div(fc[1],fc[0]);v=g.symbol('original_lower_collar_ratio_v')
    ratio_body=a.scale(a.div(at_core(core_rows['Phi'][1],v),at_core(core_rows['Phi'][0],v)),
        g.mul(g.sub(chi(v),g.one),hb,rho(v)))
    logratio=integral(g,ratio_body,v,g.zero,s,original_ratio_ODE=True,
        rho4_inlet_ratio_exact_one=True,reference_core_not_selected_as_collar=True)
    ratio=exp3(g,logratio)
    X=a.fixed(chi(s));Xm1=a.fixed(g.sub(chi(s),g.one))
    Xr=a.fixed(g.quotient(g.mul(g.sub(hb,g.one),sigma_prime(s)),
        g.mul(hb,rs),'positive exact hb and rho'))
    phi=[a.mul(ratio,fc[0]),a.mul(ratio,a.mul(X,fc[1])),
        a.mul(ratio,a.add(a.mul(Xr,fc[1]),a.mul(X,a.mul(Xm1,a.mul(A,fc[1]))),a.mul(X,fc[2])))]
    t=g.symbol('original_lower_collar_Psi_t')
    psi_body=a.scale(a.mul(at(g,ratio,s,t),at_core(core_rows['Psi'][1],t)),g.mul(chi(t),hb,rho(t)))
    psi0=a.add(at(g,core_rows['Psi'][0],rho_variable,c(4)),
        integral(g,psi_body,t,g.zero,s,actual_Psi_radial_transport=True))
    psi=[psi0,a.mul(ratio,a.mul(X,pc[1])),
        a.mul(ratio,a.add(a.mul(Xr,pc[1]),a.mul(X,a.mul(Xm1,a.mul(A,pc[1]))),a.mul(X,pc[2])))]
    w=g.symbol('original_lower_collar_mean_w')
    mean0=a.div(a.add(a.scale(at(g,core_rows['mean_Psi'][0],rho_variable,c(4)),4),
        integral(g,a.scale(at(g,psi0,s,w),g.mul(hb,rho(w))),w,g.zero,s,
            actual_own_mean_FTC=True,nonzero_rho4_incoming_mean_retained=True)),a.fixed(rs))
    mean1=a.div(a.add(psi0,a.neg(mean0)),a.fixed(rs))
    mean2=a.div(a.add(psi[1],a.neg(a.scale(mean1,2))),a.fixed(rs))
    mean=[mean0,mean1,mean2]
    U=[a.add(core_functions['U0'],a.scale(psi0,eps)),*[a.scale(q,eps) for q in psi[1:]]]
    m=[a.add(core_functions['U0'],a.scale(mean0,eps)),*[a.scale(q,eps) for q in mean[1:]]]
    Q=[own.radial_Q(g,de,u,mm) for u,mm in zip(U,m)]
    b=own.c2.C2Algebra(g);Qy=b.scale(Q[1],rs)
    Qyy=b.add(Qy,b.scale(Q[2],g.mul(rs,rs)))
    u=g.symbol('original_lower_collar_core_pressure_u')
    fcu=at(g,core_rows['Phi'][0],rho_variable,u)
    core_pressure=integral(g,a.mul(fcu,fcu),u,g.zero,c(4),same_true_core_Phi_squared=True)
    p=g.symbol('original_lower_collar_pressure_p');fp=at(g,phi[0],s,p)
    collar_pressure=integral(g,a.scale(a.mul(fp,fp),g.mul(hb,rho(p))),p,g.zero,s,
        actual_same_collar_Phi_squared=True,core_pressure_history_not_reset=True)
    F0=core_functions['F0'];F=[a.mul(F0,q) for q in phi]
    pressure=a.add(core_functions['P0'],a.scale(a.mul(a.mul(F0,F0),
        a.add(core_pressure,collar_pressure)),eps))
    return dict(query=s,domain=dict(lower=parameters['sc_quarter'],upper=parameters['sc_half']),
        radius=R,rho=rs,physical_measure=g.mul(eps,jacobian),
        weighted_physical_measure=g.mul(eps,eps,rs,jacobian),
        log_radius=g.add(g.neg(parameters['logLambda']),g.unary('log',rs)),
        ratio=ratio,logratio=logratio,chi=chi(s),chi_rho=Xr.value,
        Phi_rho_rows=phi,Psi_rho_rows=psi,mean_Psi_rho_rows=mean,F_rho_rows=F,
        U_rho_rows=U,m_rho_rows=m,Q_rho_rows=Q,U=U[0],U_y=a.scale(U[1],rs),
        m=m[0],m_y=a.add(U[0],a.neg(m[0])),Q=Q[0],Q_y=Qy,Q_yy=Qyy,
        Utheta=a.scale(F[0],g.unary('sqrt',g.mul(c(2),R))),
        Ur=b.scale(Q[0],g.unary('sqrt',g.mul(c('1/2'),R))),
        absolute_swirl_base_pressure=pressure,
        no_extra_Pstar_in_physical_U_m_Q=True,all_histories_start_at_true_core_rho4=True)


def direct_pressure_functions(g,fields,parameters):
    """Original differential N1p with both radial derivatives and true C1."""
    a=parent.C1Algebra(g);c=g.constant
    c1=lambda q:source.C1Function(q.value,q.Z)
    z=source.C1Function(g.symbol('Z'),g.one);de=a.fixed(parameters['delta'])
    d=a.sub(a.fixed(1),a.mul(z,z));L=a.sub(a.fixed(1),a.mul(de,z,z))
    h=a.sub(a.fixed(1),de)
    U,Q,Qy,Qyy=[c1(fields[key]) for key in ('U','Q','Q_y','Q_yy')]
    Qz=source.C1Function(fields['Q'].Z,fields['Q'].ZZ)
    S=a.add(a.div(a.add(Q,a.scale('1/2',a.mul(h,z,Qz)),Qy),L),
        a.mul(Q,a.add(a.scale('1/2',Q),Qy)),
        a.div(a.mul(U,a.add(a.scale(-2,a.mul(z,Q)),a.mul(d,Qz),a.scale(-2,a.mul(z,Qy)))),L))
    N=a.add(a.scale('-1/2',S),a.scale(g.quotient(g.one,fields['radius'],'positive physical collar R'),a.add(Qyy,Qy)))
    s=fields['query'];domain=fields['domain'];functions={};bodies={}
    for key,measure,sign in (('m13',fields['physical_measure'],1),('m15',fields['weighted_physical_measure'],-1)):
        body=a.scale(g.mul(c(sign),measure),N);bodies[key]=body
        functions[key]=source.C1Function(*(controls.integral(g,q,g.nodes[s.node]['name'],domain['lower'],domain['upper'],
            original_direct_N1p=True,ordinary_slow_Z_derivative_order=k,
            endpoint_independent_of_Z=True,all_finite_viscous_derivatives_retained=True,
            original_m15_minus_applied_once=True,signed_integrand_not_a_cap=True)
            for k,q in enumerate((body.value,body.Z))))
    return functions,dict(N1p=N,S=S,signed_physical_integrands=bodies,
        direct_differential_form_not_endpoint_subtraction=True)


class CurrentOriginalN1LowerCollarFunctions:
    def __init__(self,core_owner=None,require_checked=True):
        self.core=core_owner or core.CurrentOriginalN1CoreFixedPointFunctions()
        if not self.core.acceptance_loaded:raise ValueError('Checked true nonlinear core required')
        self.core.assert_graph();self.g=g=self.core.g;self.family=copy.deepcopy(self.core.family)
        self.hashes=dict(self.core.hashes)
        for name,digest in parent.linear.read(core.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for name in (core.NAME,core.RECEIPT):self._pin(name,sha(name))
        for stem,gate in (('current_inner_exit_strict_collar','current_inner_exit_strict_collar_attached_to_current_source_graph_certified'),
            ('current_generic_shear_loop_domain','current_original_generic_loop_two_sided_collars_and_repair_geometry_certified'),
            ('current_original_n1_common_inner_selection',selection.GATE),
            ('current_original_n1_unchanged_pressure_prefix',prior.GATE)):
            module_name=PREFIX+stem+'.py';receipt_name=PREFIX+stem+'_check.json'
            receipt=parent.linear.read(receipt_name)
            identity=receipt.get('source_family') or {key:receipt[key] for key in self.family}
            if not receipt.get('all_passed') or not receipt.get(gate) or identity!=self.family:
                raise ValueError('Checked same-source collar prerequisite required '+stem)
            for name,digest in receipt['input_hashes'].items():self._pin(name,digest)
            for name in (module_name,receipt_name):self._pin(name,sha(name))
        exit_name=PREFIX+'current_inner_exit_strict_collar.json';exit_record=parent.linear.read(exit_name);self._pin(exit_name,sha(exit_name))
        bits=exit_record['explicit_current_inner_exit_strict_collar']['selected_first_phase_endpoint']
        if bits['lower_exact_mpf_tuple']!=bits['upper_exact_mpf_tuple']:raise ValueError('Original exact selected sc required')
        self.sc_bits=bits['lower_exact_mpf_tuple']
        hb=next(source.FunctionRef(g,i) for i,n in enumerate(g.nodes) if n['operation']=='current_original_source_parameter' and n['name']=='hbB')
        if g.nodes[hb.node]['source_family']!=self.family:raise ValueError('Original bridge width differs')
        attachment=exit_record['exact_current_exit_source_attachment']
        if not attachment['passed'] or attachment['exact_width']!='hb=epsilon_b=cstar*K^-100; K is the original global physical norm sum':
            raise ValueError('Actual collar width definition required')
        sc=exact_binary_scalar(g,self.sc_bits);c=g.constant
        self.parameters={**self.core.parameters,'hb':hb,'sc':sc,
            'sc_quarter':g.mul(c('1/4'),sc),'sc_half':g.mul(c('1/2'),sc)}
        ctx=core.common.MPIntervalContext();ctx.dps=500;sc_value=selection.read(ctx,bits)
        selector=SimpleNamespace(c=ctx,macro=SimpleNamespace(sc=sc_value))
        support=modification.OriginalMicroFiniteN.modification_support(selector,'first_micro',0,sc_value/2)
        if not support['increments_exact_zero']:raise ValueError('Later source changes must vanish throughout this collar')
        self.geometry_certificate=dict(sc_exact_mpf_tuple=self.sc_bits,
            original_width_node=hb.node,original_width_definition=attachment['exact_width'],
            sc_and_hb_Z_independent=True,physical_Ra_is_4_over_Lambda_not_4=True,
            physical_Rin_is_Ra_exp_hb_sc_quarter=True,physical_rminus_is_Ra_exp_hb_sc_half=True,
            original_modification_support=core.common.encode(support),
            nonzero_original_sc_not_a_selected_cover_center=True)
        self.core_radial_rows=true_core_rows(self.core)
        self.fields=collar_fields(g,self.core.rho,self.core_radial_rows,self.core.functions,self.parameters)
        self.functions,self.pressure_bindings=direct_pressure_functions(g,self.fields,self.parameters)
        for module in (modification,selection,prior):self._pin(Path(module.__file__).name,sha(Path(module.__file__).name))
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(k) for k in OPEN)
                or receipt['source_family']!=self.family or raw['exact_graph_nodes']!=g.nodes
                or raw['lower_collar_bindings']!=self.snapshot):raise ValueError('Checked true lower-collar functions required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed lower-collar source '+name)
            self.acceptance_loaded=True

    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed collar dependency '+name)
        self.hashes[name]=digest

    def _snapshot(self):return own.encode(dict(fields=self.fields,functions=self.functions,
        core_radial_rows=self.core_radial_rows,parameters=self.parameters,
        pressure_bindings=self.pressure_bindings,geometry_certificate=self.geometry_certificate))

    def assert_graph(self):
        self.core.assert_graph()
        if self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes:
            raise ValueError('True lower-collar graph changed')


def run():
    began=time.monotonic();owner=CurrentOriginalN1LowerCollarFunctions(require_checked=False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,lower_collar_bindings=owner.snapshot,
        actual_integral_function_bodies=True,numerical_collar_point_accuracy_claimed=False,
        **dict.fromkeys((GATE,)+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_TRUE_LOWER_COLLAR_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw


if __name__=='__main__':run()
