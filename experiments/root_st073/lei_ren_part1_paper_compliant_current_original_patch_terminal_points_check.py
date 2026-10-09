"""Terminal patch source, native graph, finite scalar and Rh seam checks."""
import ast
import copy
from dataclasses import replace
import inspect
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_patch_terminal_points as current
import lei_ren_part1_paper_compliant_current_original_point_source_leaves_check as prior

ep=current.ep;leaves=current.leaves;base=current.base;slow=current.slow


def terminal_inertial_identity(owner):
    """Bind terminal p1/p2 directly to full five-primitive stress formulas."""
    t=owner.owner.inputs.frame.owner.template
    x=s.Symbol('terminal_patch_x',positive=True);z=t['z'];R=t['R'];ps=t['Pstar'];delta=t['delta']
    f=s.exp(-s.Rational(3,5))*x**s.Rational(1,10);C=1/(1+z*z)
    E=C*f;u=ps*E;V=4*z;L=1-delta*z*z;d=1-z*z
    h=s.Rational(5,8)*E;k=V*h;m=V
    energy=V*V/ps**2-s.Rational(5,12)*E*E;pressure=s.Rational(5,2)*E*E
    mt=s.sqrt(2)*R**s.Rational(3,2)*ps*h;mtz=s.sqrt(2)*R**s.Rational(3,2)*ps*k
    ee=ps**2*energy;P=ps**2*(t['P0']+pressure)
    W=1-(1-delta)*z*m-d*s.diff(m,z)
    Q=-W+((1-delta/2)*mt-(1-delta)*z*s.diff(mt,z)/2-d*s.diff(mtz,z)+(2*delta-1)*z*mtz)/(s.sqrt(2)*R**s.Rational(3,2)*u)
    N=-W*V+(1-delta)*(m-z*s.diff(m,z))/2+2*delta*z*ee-d*s.diff(ee,z)+2*(1+delta)*z*P-d*s.diff(P,z)
    expected={'p1':R*Q/L,'p2':R*N/(L*u)}
    substitutions={t['f']:f,t['H']:s.Rational(5,8)*f,t['D']:s.Rational(5,12)*f*f,t['P']:s.Rational(5,2)*f*f}
    identities={}
    for key in ('p1','p2'):
        for order in (0,1):
            expression=t[key if not order else key+'_Z'].xreplace(substitutions)
            wanted=s.diff(expected[key],z,order)
            assert s.cancel(expression-wanted)==0,(key,order)
            identities[key+'_Z'+str(order)]=True
    # These are the exact template object and compiler used by the projected
    # point inputs, rather than an independently selectable coefficient box.
    templates=owner.owner.inputs.templates
    assert all(templates['source_template_reconstruction_identities'].values())
    for filename,digest in t['original_program_AST_hashes'].items():assert owner.hashes[filename]==digest==current.sha(filename)
    return dict(passed=True,source_bound_full_primitive_p1_p2_and_ordinary_Z_identities=identities,
        exact_current_R_and_full_Mztheta_pressure_P0_retained=True,
        actual_projected_compiler_uses_same_source_template=True,
        original_stress_AST_hashes=t['original_program_AST_hashes'])


def independent_density_reference(owner,manifest):
    # Reuse the independently implemented finite scalar loop/Jacobian test.
    # Only its sample selection changes; reverse AST equality keeps its math.
    original=ast.parse(inspect.getsource(prior.finite_density_reference))
    changed=copy.deepcopy(original)
    selections=[n for n in ast.walk(changed) if isinstance(n,ast.ListComp) and isinstance(n.elt,ast.Name) and n.elt.id=='row']
    if len(selections)!=1 or len(selections[0].generators[0].ifs)!=1:raise ValueError('Independent finite test selection changed')
    old=selections[0].generators[0].ifs[0]
    selections[0].generators[0].ifs[0]=ast.parse("row['chart']=='actual_patch'",mode='eval').body
    reverse=copy.deepcopy(changed)
    next(n for n in ast.walk(reverse) if isinstance(n,ast.ListComp) and isinstance(n.elt,ast.Name) and n.elt.id=='row').generators[0].ifs[0]=old
    assert ast.dump(reverse)==ast.dump(original)
    env=dict(vars(prior));exec(compile(ast.fix_missing_locations(changed),'<same independent finite density test; terminal selection>','exec'),env)
    result=env['finite_density_reference'](owner,manifest)
    assert result['new_reference_and_arbitrary_O2_scalar_density_cases']==2*len(manifest['actual_original_point_source_frames'])
    return dict(passed=True,terminal_patch_finite_scalar_cases=result['new_reference_and_arbitrary_O2_scalar_density_cases'],
        independent_changed_minus_original_and_Z_Jacobian_comparisons=result['independent_changed_minus_original_and_Z_Jacobian_comparisons'],
        maximum_reference_outside_enclosure_discrepancy=result['maximum_reference_outside_enclosure_discrepancy'],
        only_test_sample_selection_projected_not_scalar_or_density_math=True,
        finite_diagnostic_units_only=result['finite_diagnostic_units_only'],
        finite_probe_not_selected_as_native_source_phase_or_parameters=True)


def native_points_and_seam(owner,manifest):
    p=mp.mp.clone();p.dps=120;frames=[];closed=0;dispatch=0;ordinary=0;formal=0
    with mp.workdps(350):
        for sample in manifest['actual_original_point_source_frames']:
            x=current.terminal_coordinate(s.E if sample['original_coordinate_exact']=='E' else sample['original_coordinate_exact'])
            frame=owner.frame(coordinate=current.phase_coordinate(x),Z=sample['original_Z_exact'],N=sample['explicit_candidate_N']);frames.append(frame)
            rows=slow.restore_point(p,sample['actual_defining_point_inputs_and_errors'])['inputs']
            z=s.Rational(sample['original_Z_exact']);z=p.mpf(int(z.p))/int(z.q)
            f=p.exp(current.patch_offset(p,x)/10);C=1/(1+z*z)
            for key,order,wanted in (('E',0,C*f),('E',1,-2*z*C*C*f),('a',0,p.mpf(4)/5),('a',1,0),('b',0,0),('b',1,0)):
                got=sum((v.coefficient for v in rows[key][order].terms),p.mpf(0))
                assert abs(got-wanted)<p.mpf('1e-45');closed+=1
            assert frame.query['roots']['b'][0,0].zero and frame.query['roots']['b'][0,1].zero
            assert not frame.query['roots']['E'][0,0].zero
            if frame.Z==0:
                assert frame.query['roots']['p2'][0,0].zero and not frame.query['roots']['p2'][0,1].zero
                assert not frame.values['Z']['p'].zero
            for key in leaves.transport.RATES:
                for order in ('C0','Z'):
                    row=owner.rows['actual_patch','density_'+key+'_'+order]
                    got=owner.dispatch(row,frame);assert got is frame.values[order][key]
                    assert got.ledger is frame.query['ledger'] and got.ctx is frame.ctx;dispatch+=1
                    try:
                        value=owner.source(row,coordinate=current.phase_coordinate(x),Z=str(frame.Z),N=frame.N)
                        assert hasattr(value,'_mpi_') and not hasattr(value,'scale')
                        assert all(mp.isfinite(v) for v in ep(value));ordinary+=1
                    except ArithmeticError:formal+=1
            for piece in frame.pieces:
                target=piece['inverse']['phase'];image=piece['inverse']['selected_inverse']['phase_image']
                assert ep(image)[0]<=ep(target)[0]<=ep(target)[1]<=ep(image)[1]
                assert piece['derivative']['p2_Z_retained'] and piece['derivative']['E_Z_term_retained']
            assert frame.record['original_pressure_coefficient_and_late_flatten_errors_retained']
        left=next(f for f in frames if f.coordinate==s.E)
        right=owner.provider.frame(chart='Rh_reference',coordinate=-5,Z=str(left.Z),N=left.N)
        assert all(ep(a)==ep(b) for a,b in zip(left.query['kernel'].q.scale.bases,right.query['kernel'].q.scale.bases,strict=True))
        seam=0
        for order in ('C0','Z'):
            for key in leaves.transport.RATES:
                a=left.values[order][key];b=prior.scalar_same_basis(right.values[order][key],a.scale.bases,a.ledger)
                lo,hi=ep((a-b).coefficient);assert lo<=0<=hi;seam+=1
        x=owner.phase.x
        assert s.simplify(owner.phase.maps['actual_patch'].subs(x,s.E)-owner.phase.maps['Rh_reference'].subs(x,-5))==0
    assert dispatch==manifest['original_role_bound_point_dispatches']==50
    assert ordinary==manifest['ordinary_interval_callbacks'] and formal==len(manifest['unmaterializable_source_values_retained_factored'])
    return dict(passed=True,true_closed_source_coefficients_checked=closed,actual_original_root_dispatches=dispatch,
        actual_directed_ordinary_interval_callbacks=ordinary,unmaterializable_source_values_retained_factored=formal,
        native_Rh_seam_C0_Z_overlap_diagnostics=seam,exact_original_phase_expression_at_Rh_equal=True,
        overlap_not_used_as_function_identity=True,terminal_background_closure_not_incoming_correction_reset=True),frames


def guards(owner,frames):
    row=owner.rows['actual_patch','density_m_C0'];frame=frames[0];g=owner.built['graph'];before=leaves.digest_rows(g.nodes)
    calls=[lambda:current.terminal_coordinate('1.7'),lambda:current.terminal_coordinate('2.72'),
        lambda:current.terminal_coordinate({'original_patch_log_offset':'.5'}),
        lambda:current.terminal_coordinate({'original_patch_log_offset':'1.01'}),
        lambda:current.terminal_coordinate('nan'),lambda:owner.frame(coordinate=2,Z='.37',N=159),
        lambda:owner.frame(coordinate=2,Z='1.01',N=257),lambda:owner.frame(coordinate=2,Z='.37',N=True),
        lambda:owner.frame(coordinate=2,Z='.37',N=257,bits=1),
        lambda:owner.frame(chart='Rh_reference',coordinate=2,Z='.37',N=257),
        lambda:owner.source(row,coordinate=2,Z='.37',N=257,phase='.3'),
        lambda:owner.source_factored(dict(row),coordinate=2,Z='.37',N=257),
        lambda:owner.dispatch(row,replace(frame,N=frame.N+1)),lambda:owner.parameter('P0'),
        lambda:owner.integrate(lambda x:x,1,2)]
    rejected=0
    for call in calls:
        try:call()
        except (ValueError,TypeError,NotImplementedError):rejected+=1
    assert rejected==len(calls)
    old=row['source_node'];row['source_node']=old+1
    try:
        try:owner.require_row(row)
        except ValueError:rejected+=1
        else:raise AssertionError('Mutated source graph admitted')
    finally:row['source_node']=old
    assert leaves.digest_rows(g.nodes)==before
    return dict(passed=True,invalid_domain_N_Z_precision_phase_row_frame_and_full_oracle_rejections=rejected,
        exact_original_graph_unchanged=True,active_support_domain_rejected=True)


@current.precise.phase.native.inlet.source_precision
def run():
    began=time.monotonic();manifest=json.loads((current.HERE/current.NAME).read_bytes())
    assert manifest[current.GATE] and manifest['mode']=='original_patch_terminal_point_partial'
    for name,digest in manifest['input_hashes'].items():assert current.sha(name)==digest,name
    owner=current.OriginalPatchTerminalPointLeaves()
    assert owner.family==manifest['source_family'] and owner.graph_digest==manifest['exact_original_transport_graph_sha256']
    assert owner.recipe==manifest['terminal_source_recipe']
    native,frames=native_points_and_seam(owner,manifest)
    flags=('actual_patch_active_support_source_installed','actual_patch_terminal_integral_installed',
        'full_17_chart_source_or_integral_oracle_installed','actual_five_controls_installed',
        'actual_terminal_Z_function_closure_installed','current_whole_N_selected',*current.precise.phase.packets.OPEN)
    assert all(manifest[key] is False for key in flags)
    report=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        terminal_source_to_full_inertial_compiler_identities=terminal_inertial_identity(owner),
        terminal_native_source_and_Rh_seam_checks=native,
        independent_finite_scalar_density_and_Z_checks=independent_density_reference(owner,manifest),
        source_domain_and_partial_oracle_guards=guards(owner,frames),**dict.fromkeys(flags,False),
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),Path(__file__).name:current.sha(Path(__file__).name),
            Path(prior.__file__).name:current.sha(Path(prior.__file__).name)},execution_seconds=time.monotonic()-began,
        scope=manifest['scope'])
    (current.HERE/current.RECEIPT).write_bytes(json.dumps(current.precise.encode(report),indent=2).encode()+b'\n')
    print('Original actual-patch terminal: source, density/Z, native Rh seam and guards PASS',flush=True)
    return report


if __name__=='__main__':run()
