"""Actual quiet-flow source, C5 calculus, native mass and pressure identities."""
import json
from pathlib import Path
import time
from types import SimpleNamespace
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_true_quiet_Rp_C5_pressure_C3_functions as current
import lei_ren_part1_paper_compliant_current_original_true_compact_C5_functions_check as prior
import lei_ren_part1_paper_compliant_current_original_outer_leading_C3_identity_check as leading_check

lexical=prior.lexical


class Algebra(prior.Algebra):
    def at(self,i):
        if i in self.cache:return self.cache[i]
        row=self.g.nodes[i]
        if row['operation']=='analytic_unary' and row['name']=='exprel':
            x=self.at(row['argument']);self.cache[i]=sy.Integer(1) if x==0 else (sy.exp(x)-1)/x
            return self.cache[i]
        return super().at(i)


def actual_source_geometry(owner):
    g=owner.g;free={i:set(v) for i,v in owner.compact.free_source_seeds.items()};r=current.rows;native=g.nodes[owner.s.node]['name']
    for i,row in enumerate(g.nodes):
        if row['operation']=='compact_raw_beta':free[i]=lexical.free_variables(g,row['argument'],free)
    assert owner.compact.acceptance_loaded and owner.lower.acceptance_loaded
    assert owner.incoming is owner.compact.terminal and owner.incoming_integrals is owner.compact.terminal['integrals']
    assert owner.relative_incoming is owner.compact.terminal['correction'] and owner.initial is owner.compact.initial
    assert owner.Tw is owner.compact.outer.Tw
    for name in ('N','mu','Pstar','P0','delta'):assert getattr(owner,name) is getattr(owner.compact,name)
    count=0
    for data,order in ((owner.fields.values(),5),(owner.radial.values(),5),(owner.generic.values(),4),
        (owner.relative.values(),5),(owner.pressure.values(),3),(owner.Rp_frame.values(),5),
        (owner.lower_pressure.values(),3),(owner.partitioned_Rin_Rp.values(),3)):
        for q in data:
            assert len(r(q))==order+1
            for ref in r(q):assert lexical.free_variables(g,ref,free)<={'Z',native},ref.node
            count+=order+1
    assert lexical.free_variables(g,owner.K,free)==set()
    for q in owner.quiet_closed.values():assert r(q)[1:]==(g.zero,)*3
    for flow in (owner.flow,owner.leading_flow):
        for group in (flow['profiles_y0_to_y4'],flow['histories_y0_to_y4']):
            assert len(group)==5 and all(len(r(q))==flow['order']+1 for row in group for q in row.values())
    for key,q in owner.pressure_local.items():
        for n,ref in enumerate(r(q)):
            row=g.nodes[ref.node]
            assert row['lower']==g.zero.node and row['upper']==owner.s.node
            assert row['physical_radius_power']==(1 if key=='S0' else 2)
            assert row['physical_Jacobian_applied_once'] and row['ordinary_slow_Z_derivative_order']==n
    Rc=sy.Symbol('Rc',positive=True);Tw=sy.Symbol('Tw',positive=True);S=sy.Symbol('s')
    alg=Algebra(g,{owner.compact.Rc.node:Rc,owner.Tw.node:Tw,owner.s.node:S})
    assert sy.simplify(alg.at(owner.R2.node)-2*Rc)==0
    assert sy.simplify(alg.at(owner.length.node)-(Tw-2-sy.log(2)))==0
    assert sy.simplify(alg.at(owner.Rp.node)-Rc*sy.exp(Tw-2))==0
    assert owner.compact.outer.mu_proof['positive_mu_and_power_2mu_ge_eta_and_Rp_reservation']
    assert owner.definition['inner_0_Rin_C3_and_lower_active_field_seam_not_proved']
    return dict(passed=True,actual_quiet_frame_radial_pressure_partition_rows=count,
        true_compact_exit_and_all_original_nonzero_histories_preserved=True,
        exact_2Rc_s_Tw_minus2_minus_log2_Rp_geometry=True,
        original_constant_K_source_Z_independent=True,
        lower_component_separate_and_inner_C3_native_pulse_seam_norms_not_claimed=True)


def actual_flow_calculus(owner):
    g=owner.g;Z=sy.Symbol('Z');S=sy.Symbol('quiet_s',real=True);point=sy.Rational(3,8);mu=sy.Rational(1,13)
    A=2+Z+Z**6/19;H={k:3+i+Z**2/(23+i)+Z**6/(31+i) for i,k in enumerate(current.RATES)}
    alpha=sy.Rational(1,2)+mu;E=A*sy.exp(-alpha*S)
    decay=lambda r:sy.exp(-sy.Rational(r)*S)
    mass=lambda k:(1-sy.exp(-k*S))/k
    hist=dict(m=H['m']*decay(1),h=H['h']*decay(sy.Rational(3,2))+A*decay(sy.Rational(3,2))*(sy.exp((1-mu)*S)-1)/(1-mu),
        k=H['k']*decay(sy.Rational(3,2)),e=H['e']*decay(1)-A*A*decay(1)*mass(2*mu)/2,
        p=H['p']+A*A*mass(1+2*mu)/2)
    count=0
    for flow,initial in ((owner.flow,owner.incoming['modified_fields']),(owner.leading_flow,owner.leading_incoming)):
        bindings={owner.s.node:S,g.symbol('Z').node:point,owner.mu.node:mu}
        for key,body in dict(E=A,**H).items():
            for n,ref in enumerate(current.rows(initial[key])):bindings[ref.node]=sy.diff(body,Z,n).subs(Z,point)
        alg=Algebra(g,bindings)
        for j,(profiles,history) in enumerate(zip(flow['profiles_y0_to_y4'],flow['histories_y0_to_y4'])):
            for key,q in profiles.items():
                body=E if key=='E' else sy.Integer(0)
                for n,ref in enumerate(current.rows(q)):
                    assert sy.simplify(alg.at(ref.node)-sy.diff(body,S,j,Z,n).subs(Z,point))==0,('profile',flow['order'],j,key,n);count+=1
            for key,q in history.items():
                for n,ref in enumerate(current.rows(q)):
                    assert sy.simplify(alg.at(ref.node)-sy.diff(hist[key],S,j,Z,n).subs(Z,point))==0,('history',flow['order'],j,key,n);count+=1
    # Execute the actual legacy defining helper, with the same functions.
    fg=current.source.FunctionTransportGraph();a3=current.original_flow.target.C3Function
    make=lambda body:a3(*current.rows(prior.jet(fg,body,3)))
    fs=fg.symbol('fixture_s');old=current.original_flow.power_flow(fg,make(A),{k:make(v) for k,v in H.items()},fg.constant(mu),fs)
    oldalg=Algebra(fg,{fs.node:S,fg.symbol('Z').node:point})
    for j,row in enumerate(old['complete_histories_y0_y1_y2_y3_y4_C3']):
        for k,q in row.items():
            for n,ref in enumerate(current.rows(q)):assert sy.simplify(oldalg.at(ref.node)-sy.diff(hist[k],S,j,Z,n).subs(Z,point))==0
    return dict(passed=True,actual_C5_C6_profile_history_y0_to_y4_calculus_rows=count,
        original_defining_C3_helper_prefix_rows=100,finite_mu_exprel_and_all_incoming_histories=True,
        higher_fixture_source_rows_not_original_numeric_magnitudes=True)


def actual_inlet(owner):
    g=owner.g;Z=sy.Symbol('Z');point=sy.Rational(3,8);bindings={g.symbol('Z').node:point,owner.s.node:0}
    bodies={k:2+i+Z**2/(17+i)+Z**5/(23+i) for i,k in enumerate(owner.incoming['modified_fields'])}
    for k,q in owner.incoming['modified_fields'].items():
        for n,ref in enumerate(current.rows(q)):bindings[ref.node]=sy.diff(bodies[k],Z,n).subs(Z,point)
    for n,ref in enumerate(current.rows(owner.P0)):bindings[ref.node]=sy.diff(bodies['P0'],Z,n).subs(Z,point)
    # The original compact bumps are flat at2Rc: its incoming V is exactly0.
    for ref in current.rows(owner.incoming['modified_fields']['V']):bindings[ref.node]=0
    bodies['V']=sy.Integer(0);alg=Algebra(g,bindings);count=0
    for k,q in owner.inlet['fields'].items():
        body=bodies['P0']+bodies['p'] if k=='pressure' else bodies[k]
        for n,ref in enumerate(current.rows(q)):
            assert sy.simplify(alg.at(ref.node)-sy.diff(body,Z,n).subs(Z,point))==0,(k,n);count+=1
    for k,q in owner.relative_incoming.items():
        for n,ref in enumerate(current.rows(q)):bindings[ref.node]=sy.Rational(7+n,31+n)
    alg=Algebra(g,bindings)
    for k,q in owner.inlet['relative'].items():
        for n,ref in enumerate(current.rows(q)):assert alg.at(ref.node)==sy.Rational(7+n,31+n);count+=1
    return dict(passed=True,actual_complete_and_relative_2Rc_inlet_rows=count,
        compact_flat_velocity_seam_and_nonzero_pressure_history_not_reset=True,
        relative_rows_tested_nonzero_to_reject_zero_inlet_replacement=True,mathematical_fixture_only=True)


def native_mass_identity(owner):
    """Actual source-node K equality, with the original turnoff integral."""
    module=current.previous.previous.original_leading
    lead=current.parent.linear.read(module.NAME);frame=current.parent.linear.read(current.original_quiet.FRAME_NAME)
    quiet=current.parent.linear.read(current.original_flow.NAME)
    og=prior.graph(frame['exact_graph_nodes']);params=quiet['actual_C3_repaired_exit_power_and_pulse_functions']['parameters']
    field=SimpleNamespace(graph=og,contracts=module.source_contracts(),parameters={k:current.source.FunctionRef(og,v) for k,v in params.items()})
    q=leading_check.Interpreter(field);Tw=sy.Symbol('Tw',positive=True)
    ff=frame['actual_Rp_native_frame_functions'];q.bindings[ff['Tw']]=Tw
    axial=lead['actual_original_leading_C3_source_functions']['charts']['O2_axial'];kernel=axial['kernels']['original_B_mass']
    oldfn=q.kernel_functions.setdefault(kernel,sy.Function('original_turnoff_mass'))
    native_mean=q.at(ff['actual_identified_Rp_histories_C3']['m'][0]);native_kappa=sy.cancel(native_mean/q.z)
    bg=owner.compact.outer.background['O2_axial'];local=bg['history_contributions']['m'].value
    X=sy.Symbol(owner.g.nodes[bg['phase'].node]['name'],real=True);T=sy.Symbol('turnoff_t',real=True)
    yc=sy.Symbol(og.nodes[axial['kernels']['original_y_endpoint']]['name'],real=True)
    oldrow=og.nodes[kernel];u=sy.Symbol(oldrow['variable'],real=True)
    class MassAlgebra(Algebra):
        def at(self,i):
            if i in self.cache:return self.cache[i]
            row=self.g.nodes[i]
            if row['operation']=='exact_original_source_flat_sigma':
                self.cache[i]=leading_check.Sigma(self.at(row['argument']));return self.cache[i]
            return super().at(i)
    B=4*q.z*sy.exp(-q.logP);nb={owner.g.symbol('Z').node:q.z,owner.Pstar.node:sy.exp(q.logP),owner.mu.node:q.mu}
    Rh=owner.compact.outer.leading_Rh['fields']
    for key in ('m','V'):
        for n,ref in enumerate(current.rows(Rh[key])):nb[ref.node]=sy.diff(B,q.z,n)
    newrow=owner.g.nodes[local.node];tref=owner.g.symbol(newrow['variable'])
    body_bind=dict(nb);body_bind[tref.node]=T;body_bind[bg['phase'].node]=X
    actual_body=MassAlgebra(owner.g,body_bind).at(newrow['integrand'])/B
    old_body=q.at(oldrow['integrand']).subs({u:sy.exp(40*T),yc:sy.exp(40*X)})*40*sy.exp(40*T)
    assert sy.simplify(actual_body-old_body)==0
    assert MassAlgebra(owner.g,body_bind).at(newrow['lower'])==0 and MassAlgebra(owner.g,body_bind).at(newrow['upper'])==X
    assert q.at(oldrow['lower'])==1 and q.at(oldrow['upper'])==yc
    # Use only the just-identified exact kernel, not a selected range value.
    nb[local.node]=B*oldfn(sy.exp(40*X));nb[owner.Tw.node]=Tw;nb[bg['phase'].node]=X
    Rc=sy.Symbol('Rc',positive=True);nb[owner.compact.Rc.node]=Rc
    actual_K=MassAlgebra(owner.g,nb).at(owner.K.node)
    expected_K=Rc*sy.exp(Tw-2)*sy.exp(q.logP)*native_kappa
    assert leading_check.zero(actual_K-expected_K)
    return dict(passed=True,actual_compact_derived_K_equals_original_Rp_Pstar_kappa=True,
        exact_O2_turnoff_kernel_change_y_equals_exp40x_and_Jacobian=True,
        actual_Rh_source_identity_and_same_original_native_mass_FunctionRef_used=True,
        common_Rp_Rc_ratio_and_finite_mu_Tw_recipe_bound=True,
        no_numeric_kernel_value_or_range_endpoint_used=True)


def actual_physical_pressure(owner):
    g=owner.g;Z=sy.Symbol('Z');S=sy.Symbol('quiet_s',real=True);point=sy.Rational(3,8)
    P=sy.Rational(3,2);R2=sy.Rational(7,3);kap=sy.Rational(5,11);K=kap*R2;mean=kap*Z*sy.exp(-S)
    bindings={g.symbol('Z').node:point,owner.s.node:S,owner.Pstar.node:P,owner.delta.node:sy.Rational(1,5),owner.R2.node:R2,owner.K.node:K}
    bind=lambda q,body:[bindings.__setitem__(ref.node,sy.diff(body,Z,n).subs(Z,point)) for n,ref in enumerate(current.rows(q))]
    bind(owner.fields['m'],mean/P);bind(owner.fields['V'],sy.Integer(0));bind(owner.radial['V_y'],sy.Integer(0));bind(owner.radial['V_yy'],sy.Integer(0))
    bind(owner.radial['m_y'],-mean/P)
    I0=13+Z**3/17;I1=19+Z**3/23;B0=29+Z**3/31;B1=37+Z**3/41
    bind(owner.incoming_integrals['S0'],I0);bind(owner.incoming_integrals['S1'],I1)
    bind(owner.initial['Q_plus_Q_y'],B0);bind(owner.initial['R_Q_y'],B1)
    initial=dict(m13=-I0/2-B0,m15=I1/2-K+B1)
    for k,q in owner.incoming_pressure.items():bind(q,initial[k])
    alg=Algebra(g,bindings);wantQ=dict(Q=-kap*sy.exp(-S),Q_y=kap*sy.exp(-S),Q_yy=-kap*sy.exp(-S))
    for k,body in wantQ.items():
        for n,ref in enumerate(current.rows(owner.physical[k])):assert sy.simplify(alg.at(ref.node)-sy.diff(body,Z,n))==0,(k,n)
    for n,ref in enumerate(current.rows(owner.N1p)):assert sy.simplify(alg.at(ref.node)-(kap**2*sy.exp(-2*S)/4 if n==0 else 0))==0
    local=dict(m13=K*K*(1/R2-sy.exp(-S)/R2)/4,m15=-K*K*S/4)
    for k,q in owner.pressure.items():
        for n,ref in enumerate(current.rows(q)):
            want=sy.diff(initial[k]+local[k],Z,n).subs(Z,point)
            assert sy.simplify(alg.at(ref.node)-want)==0,(k,n,'actual general')
            assert sy.simplify(alg.at(current.rows(owner.closed_pressure[k])[n].node)-want)==0,(k,n,'actual closed')
    for k,q in owner.lower_pressure.items():bind(q,43+Z**3/47)
    for k,q in owner.Rp_endpoint['pressure'].items():bind(q,53+Z**3/59)
    alg=Algebra(g,bindings)
    for k,q in owner.partitioned_Rin_Rp.items():
        for n,ref in enumerate(current.rows(q)):assert sy.simplify(alg.at(ref.node)-sy.diff(96+Z**3/47+Z**3/59,Z,n).subs(Z,point))==0
    return dict(passed=True,actual_Q_Qy_Qyy_C4_and_direct_N1p_C3_identity=True,
        actual_general_pressure_equals_closed_original_quiet_primitive_C3_rows=8,
        finite_nonzero_initial_pressure_and_original_rminus_boundary_retained=True,
        lower_plus_preRp_partition_C3_rows=8,inner_prefix_not_added_or_claimed=True,
        physical_R_and_Rweighted_signs_checked=True,mathematical_fixture_only=True)


def run(owner):
    began=time.monotonic();owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['true_quiet_Rp_C5_C3_bindings']==owner.snapshot and raw['exact_graph_nodes']==owner.g.nodes[:len(raw['exact_graph_nodes'])]
    checks={}
    for name,fn in (('actual_source_geometry',actual_source_geometry),('actual_flow',actual_flow_calculus),
        ('actual_inlet',actual_inlet),('native_mass_identity',native_mass_identity),('actual_physical_pressure',actual_physical_pressure)):
        checks[name]=fn(owner);print('TRUE_QUIET_'+name.upper()+'_CHECKED',flush=True)
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **dict.fromkeys(current.GATES,True),**dict.fromkeys(current.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('TRUE_QUIET_RP_C5_AND_PRE_RP_PRESSURE_C3_CHECKED',flush=True);return receipt


if __name__=='__main__':raise SystemExit('run(owner) reuses accepted actual compact/lower sources')
