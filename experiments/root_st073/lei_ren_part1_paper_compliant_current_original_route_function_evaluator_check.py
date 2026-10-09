"""Fresh issued graph dispatch and exact affine-boundary preservation checks."""
from dataclasses import replace
import json
from pathlib import Path
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_route_function_evaluator as current
import lei_ren_part1_paper_compliant_current_native_Rc_functional_controls as controls

ep=current.ep


def interval(c,row):
    return c.mpf((mp.mp.make_mpf(tuple(row['lower_exact_mpf_tuple'])),
                  mp.mp.make_mpf(tuple(row['upper_exact_mpf_tuple']))))


def contains(a,b):
    lo,hi=ep(a);left,right=ep(b);assert lo<=left<=right<=hi,(lo,left,right,hi)


def numerical(owner,manifest):
    g=owner.built['graph'];c=owner.c;integrals=operators=decays=points=conditional=0
    frames=[]
    for saved,history in zip(manifest['actual_original_integral_expression_frames'],manifest['actual_original_affine_history_operator_frames']):
        N=saved['explicit_candidate_N'];Z=saved['original_Z_exact']
        assert sy.Rational(Z)==sy.Rational(37,100)
        assert sy.Rational(history['original_Z_exact'])==sy.Rational(Z) and history['explicit_candidate_N']==N
        frame=owner.history_frame(Z=Z,N=N);frames.append(frame)
        with mp.workdps(c.dps+40):
            values={}
            for node in owner.integrals:
                value=owner.evaluate(current.transport.FunctionRef(g,node),Z=Z,N=N)
                assert value.ctx is c and not hasattr(value,'scale')
                assert ep(value)==ep(interval(c,saved['actual_issued_integral_values'][str(node)]))
                assert ep(value)==ep(owner.integrate(g.nodes[node],Z=Z,N=N))
                values[node]=value;integrals+=1
            high=MPIntervalContext();high.dps=c.dps+100
            prior={(key,jet):high.mpf(0) for key in current.transport.RATES for jet in ('C0','Z')}
            saved_ops={(row['route_index'],row['key'],row['ordinary_order']):row for row in history['actual_six_route_forcing_and_homogeneous_intervals']}
            scalar=mp.mp.clone();scalar.dps=c.dps+100
            for index,cell in enumerate(owner.route):
                left,_=current.reference.symbolic_node(owner.provider,cell['lower'])
                right,_=current.reference.symbolic_node(owner.provider,cell['upper'])
                width=right-left;total=sy.Integer(5) if index==0 else sy.Integer(5)+right
                for key,rate in current.transport.RATES.items():
                    r=high.mpf(rate.numerator)/rate.denominator
                    step=high.exp(-r*high.mpf(int(width.p))/int(width.q))
                    true_decay=scalar.exp(-scalar.mpf(rate.numerator)/rate.denominator*scalar.mpf(int(total.p))/int(total.q))
                    for jet,label in (('C0','value'),('Z','Z')):
                        role=index,key,jet;row=saved_ops[role]
                        assert ep(frame.forcing[role])==ep(interval(c,row['additive_original_integral_contribution']))
                        assert ep(frame.decay[role])==ep(interval(c,row['original_upstream_multiplier']))
                        assert ep(frame.decay[role])[0]<=true_decay<=ep(frame.decay[role])[1];decays+=1
                        prior[key,jet]=step*prior[key,jet]+high.mpf(ep(values[cell['contributions'][key][label]]))
                        contains(frame.forcing[role],prior[key,jet]);operators+=1
            # Explicit hypothetical nonzero boundary intervals exercise the
            # affine action only; they are never identified as original data.
            incoming={role:c.mpf((str((i+1)/1000),str((i+2)/1000))) for i,role in enumerate(owner.boundaries.values())}
            action=owner.apply_assumed_boundary(frame,incoming)
            for role,value in action.items():
                assert ep(value)==ep(frame.decay[role]*incoming[role[1:]]+frame.forcing[role]);conditional+=1
            assert frame.record['unknown_incoming_remains_explicit_not_zeroed']
            assert frame.record['separate_original_P0_and_P0_Z_unchanged']
            assert frame.record['additive_operator_component_not_actual_full_history']
    for row in manifest['actual_original_source_expression_points']:
        node=row['node'];source=g.nodes[node];variable=g.nodes[source['coordinate']]['name']
        Z=manifest['actual_original_integral_expression_frames'][0]['original_Z_exact']
        value=owner.evaluate(current.transport.FunctionRef(g,node),Z=Z,N=160,variables={variable:row['coordinate']})
        assert ep(value)==ep(interval(c,row['value']))
        assert ep(value)==ep(owner.source(source,coordinate=row['coordinate'],Z=Z,N=160));points+=1
    assert len(owner.identities)==60 and all(row['exact_original_affine_history_identity_passed'] for row in owner.identities)
    return dict(passed=True,fresh_current_N_affine_frames=len(frames),actual_issued_graph_integral_expression_results=integrals,
        actual_source_expression_results=points,exact_original_affine_outgoing_identities=60,
        independent_higher_precision_segment_recurrences_enclosed_by_graph_forcing=operators,
        independent_higher_precision_homogeneous_decay_checks=decays,
        conditional_nonzero_boundary_action_checks=conditional,
        conditional_probe_is_explicit_hypothesis_not_original_upstream_data=True,
        true_source_phase_derived_by_owner_and_not_naive_graph_fractional_part=True,
        parent_producers_and_quadrature_suites_not_reexecuted=True)


def parameters(owner):
    c=owner.c;phase=owner.provider.phase;b=phase.binder;checked=0
    with mp.workdps(c.dps+40):
        for name in owner.built['parameters']:
            value=owner.parameter(name);assert value.ctx is c and not hasattr(value,'scale')
            lo,hi=ep(value)
            if name in ('logP','Tw','Md','sc'):
                original=phase.analytic(phase.symbols[name],c);a,z=ep(original);assert lo<=a<=z<=hi
            elif name in ('logC','T'):
                a,z=ep(c.mpf(b.fixed[name]));assert lo<=a<=z<=hi
            else:
                assert lo==0 and hi>0
                log=(c.mpf(b.loghB if name=='hbB' else b.loghS) if name in ('hbB','hbS')
                    else c.ln(c.mpf(1)/1000)-4*phase.analytic(phase.symbols['logP'],c))
                assert ep(log)[1]<ep(c.ln(c.mpf(hi)))[0]
            checked+=1
    return dict(passed=True,original_source_recipe_parameter_enclosures=checked,
        analytic_logP_Tw_not_selected_as_values=True,selected_dyadic_logC_T_sc_retained=True,
        microscopic_widths_and_mu_have_nonzero_directed_upper_tails_not_zero_values=True,
        zero_lower_tail_does_not_supply_strict_positive_quotient_certificate=True)


def guards(owner):
    g=owner.built['graph'];frame=owner.history_frame(Z='.37',N=160)
    root=current.transport.FunctionRef(g,next(iter(owner.integrals)));source=owner.provider.rows['O2_slope','density_m_C0']
    before=current.leaves.digest_rows(g.nodes);rejections=0
    def reject(call):
        nonlocal rejections
        try:call()
        except (ValueError,TypeError,ArithmeticError,NotImplementedError,current.transport.SourceOracleRequired):rejections+=1
        else:raise AssertionError('Unsupported or invalid original evaluator request accepted')
    saved=json.loads((current.HERE/current.transport.NAME).read_bytes())
    unsupported=next(row for row in g.nodes if row['operation']=='original_function_graph' and row['chart'] not in current.leaves.SUPPORTED)
    for call in (
        lambda:owner.evaluate(current.transport.FunctionRef(current.transport.FunctionTransportGraph(),0),Z='.37',N=160),
        lambda:owner.evaluate(root,Z=0,N=160),lambda:owner.evaluate(root,Z=2,N=160),
        lambda:owner.evaluate(root,Z='.37',N=159),lambda:owner.evaluate(root,Z='.37',N=True),
        lambda:owner.source(source,coordinate='.537',Z='.37',N=160,phase=owner.c.mpf('.1')),
        lambda:owner.source(unsupported,coordinate='.5',Z='.37',N=160),
        lambda:owner.integrate(dict(g.nodes[root.node]),Z='.37',N=160),
        lambda:owner.apply_assumed_boundary(replace(frame),{}),lambda:owner.apply_assumed_boundary(frame,{}),
        lambda:owner.ordinary({'value':0}),lambda:owner.ordinary(mp.mp.mpf(1)),
        lambda:owner.parameter('unknown'),
        lambda:owner.evaluate(current.transport.FunctionRef(g,source['phase']),Z='.37',N=160),
        lambda:current.transport.evaluate(owner.built,root,oracle=owner,Z='.37',N=160,ctx=owner.c),
        lambda:controls.FunctionEvaluator(owner.built,oracle=owner,Z='.37',N=160,ctx=owner.c)):
        reject(call)
    for cell in owner.route:
        for pair in cell['outgoing'].values():
            for node in pair.values():reject(lambda node=node:owner.evaluate(current.transport.FunctionRef(g,node),Z='.37',N=160))
    for pair in saved['exact_Rc_target_roots'].values():
        for node in pair.values():reject(lambda node=node:owner.evaluate(current.transport.FunctionRef(g,node),Z='.37',N=160))
    row=g.nodes[root.node];old=row['upper'];row['upper']=row['lower']
    try:reject(lambda:owner.evaluate(root,Z='.37',N=160))
    finally:row['upper']=old
    assert before==current.leaves.digest_rows(g.nodes)==owner.provider.graph_digest
    return dict(passed=True,invalid_unsupported_and_full_oracle_rejections=rejections,
        all60_actual_history_roots_and_all10_terminal_target_roots_reject_missing_upstream_sources=True,
        cached_results_do_not_bypass_graph_integrity=True,original_graph_unchanged=True,
        generic_full_function_and_control_evaluators_still_reject_partial_mode=True)


@current.precise.phase.native.inlet.source_precision
def run():
    began=time.monotonic();manifest=json.loads((current.HERE/current.NAME).read_bytes())
    assert manifest[current.GATE]
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    flags=('actual_original_upstream_history_installed','full_scalar_control_evaluator_compatibility_installed',
        'full_17_chart_source_or_24_cell_integral_oracle_installed','actual_five_controls_installed',
        'actual_terminal_Z_function_closure_installed','current_whole_N_selected',*current.precise.phase.packets.OPEN)
    assert all(manifest[key] is False for key in flags)
    owner=current.OriginalRouteFunctionEvaluator();assert owner.family==manifest['source_family']
    assert owner.source_graph_sha256==manifest['source_graph_sha256']
    checks=dict(genuine_graph_evaluation_and_affine_operator=numerical(owner,manifest),
        original_parameter_recipe_and_tail_enclosures=parameters(owner),
        issued_graph_upstream_boundary_and_partial_scope_guards=guards(owner))
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        source_graph_sha256=owner.source_graph_sha256,**checks,**dict.fromkeys(flags,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Fresh genuine issued six-route source/integral evaluation and original affine history operators, explicit unknown upstream functions and conditional action only. No full upstream history, control compatibility, full-Z terminal closure, global N, actual recursion or corrected NS.')
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.precise.encode(result),indent=2).encode()+b'\n')
    print('Original node-aware routes and explicit upstream affine history PASS',flush=True)
    return result


if __name__=='__main__':run()
