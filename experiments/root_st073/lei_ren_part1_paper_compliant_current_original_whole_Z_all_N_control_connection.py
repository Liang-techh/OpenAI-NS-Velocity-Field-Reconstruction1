"""Current 17-chart integral functions and a repair-only exact frequency.

Ranges certify the defining functions; they never become function values.
The double-dyadic integer is exact but deliberately not materialized. This
connects the current source recipes to an exact C1 control limit, not a
numerically evaluated field or a globally admissible reconstruction.
"""
from dataclasses import dataclass
from fractions import Fraction
import ast
import gzip
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
import time
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions as current
import lei_ren_part1_paper_compliant_current_native_Rc_all_N_function_controls as algebra

controls,source,repair=algebra.controls,algebra.source,algebra.repair
packets=repair.packets
HERE,PREFIX,sha,ep=current.HERE,current.PREFIX,current.sha,current.ep
NAME=PREFIX+'current_original_whole_Z_all_N_control_connection.json.gz'
RECEIPT=PREFIX+'current_original_whole_Z_all_N_control_connection_check.json'
GATE='current_original_17_chart_all_N_integral_graph_repair_integer_and_C1_limit_connected'
CHARTS=('first_micro','second_micro','frozen_macro','first_switch','second_switch','post_power',
    'long_reshape','reference','restoration','postrestore','actual_patch',*current.CHARTS)
RATES={key:Fraction(value) for key,value in current.RATES.items()}


def require(test,message):
    if not test:raise ValueError(message)


def ast_binding(fn):
    text=textwrap.dedent(inspect.getsource(fn))
    name=Path(inspect.getfile(fn)).name
    return dict(module=name,module_sha256=sha(name),function=fn.__qualname__,
        AST_sha256=hashlib.sha256(ast.dump(ast.parse(text)).encode()).hexdigest())


@dataclass(frozen=True)
class DoubleDyadicInteger:
    """Strict finite integer expression N=2^(2^J), not a float/string/int alias."""
    J: int

    def __post_init__(self):
        if type(self.J) is not int or self.J<12:raise ValueError('Exact integer J>=12 required')

    def log_interval(self,c):
        lo,hi=ep(c.ln(2))
        # ldexp changes a binary exponent; it never constructs 2^J or N.
        return c.mpf((mp.ldexp(lo,self.J),mp.ldexp(hi,self.J)))

    def record(self):
        return dict(type='DoubleDyadicInteger',J=self.J,definition='N=2^(2^J)',
            exact_positive_integer=True,Z_independent=True,
            outer_integer_and_inner_power_not_materialized=True)

    def prove(self,c,N0,required):
        require(type(N0) is int and 1<=N0 and N0.bit_length()<=4096,'Admitted fixed prefix N0 required')
        require(ep(self.log_interval(c))[0]>=ep(required)[1],'Directed log(N) fails repair threshold')
        return dict(**self.record(),log_N=self.log_interval(c),required_log_N_upper=required,
            N_ge_N0=True,N0_bit_length=N0.bit_length(),
            prefix_proof='J>=12 => 2^J>=4096>=bit_length(N0); N=2^(2^J)>N0',
            log_N_lower_ge_required_upper=True,
            proof='For binary upper L<2^E choose J>=E+2; log(2)>1/2 gives log(N)>2^(E+1)>L',
            legacy_Python_int_guards_not_bypassed=True,
            global_higher_y_stress_join_frequency_conditions_not_admitted=True)


def choose_integer(c,N0,required):
    upper=ep(required)[1]
    require(upper>=0,'Nonnegative directed log-frequency threshold required')
    sign,man,exponent,bits=upper._mpf_
    require(not sign,'Positive binary threshold required')
    E=0 if not man else exponent+bits
    return DoubleDyadicInteger(max(12,E+2))


def exact_geometry():
    """Actual current local coordinates, independently of the older 24-cell route."""
    P,C,T,B,S,sc,x=sy.symbols('logP logC T hbB hbS sc x',real=True)
    Y=sy.log(25)+1000+4*P
    L=Y+sy.log(sy.Rational(11,10))-B*sc/2
    Rm=L+10*(C+P)-6
    offsets=dict(first_micro=B*(x-sc/2),second_micro=B*(x-sc/2),
        frozen_macro=Y*x+B*(2*(1-x)-sc/2),
        first_switch=Y+S*x-B*sc/2,second_switch=Y+S*(1+x)-B*sc/2,
        post_power=Y+sy.log(sy.Rational(11,10))*x+2*S*(1-x)-B*sc/2,
        long_reshape=L+T*x,reference=L+T+(10*(C+P)-T-8)*x,
        restoration=L+10*(C+P)-8+x,postrestore=L+10*(C+P)-7+x,
        actual_patch=Rm+sy.log(x),Rh_reference=Rm+6+x,O2_slope=Rm+6+x,
        O2_axial=Rm+6+sy.exp(40*x),O2_buffer=Rm+6+sy.exp(40)+x,
        O3_transition=Rm+17+sy.exp(40)+x,O3_power=Rm+18+sy.exp(40)+x)
    domains={key:(sy.Integer(0),sy.Integer(1)) for key in CHARTS}
    domains.update(first_micro=(sc/2,sy.Integer(1)),second_micro=(sy.Integer(1),sy.Integer(2)),
        actual_patch=(sy.Integer(1),sy.E),Rh_reference=(sy.Integer(-5),sy.Integer(0)),
        O2_buffer=(sy.Integer(0),sy.Integer(11)),O3_power=(sy.Integer(0),sy.Integer(2)))
    return offsets,domains,x


def native_partitions():
    result={chart:partition for chart,partition,_ in current.core.current.WholeZSharpBridgeFunctions.partitions(None)}
    result.update({chart:current.previous.previous.PARTITION for chart in CHARTS[3:6]})
    result.update({chart:current.previous.LONG_PARTITION for chart in CHARTS[6:10]})
    result['actual_patch']=current.previous.PATCH_PARTITION
    result.update(current.PARTITIONS)
    return result


def symbolic_coordinate(value):
    if value==('sc_half',):return sy.Symbol('sc',real=True)/2
    if value=='Rh':return sy.E
    return sy.Rational(*value)


def recipe_manifest():
    """Bind the exact source recipes and quantity paths, never their range outputs."""
    core=current.core;switch=current.previous.previous;long=current.previous
    providers={}
    groups=((CHARTS[:3],core.current.WholeZSharpBridgeFunctions.source_query),
        (CHARTS[3:6],switch.WholeZAllNSwitchFunctions.leading_packet),
        (CHARTS[6:11],long.WholeZAllNLongPatchFunctions.leading_packet),
        (CHARTS[11:],current.WholeZAllNOuterRcFunctions.leading_packet))
    for charts,fn in groups:
        binding=ast_binding(fn)
        for chart in charts:providers[chart]=dict(binding,
            defining_provider=fn.__qualname__,
            native_chart=chart,
            quantity_paths=dict(E='original_generic_source.common_velocity_E_axial5[Z_order]',
                V='original_generic_source.common_velocity_V_axial5[Z_order]',
                a='original_roots.a[(0,Z_order)]',t0='original_roots.t0[(0,Z_order)]',
                p2='original_roots.p2[(0,Z_order)]',
                b='original_full_source_quotients.actual_b_axial5[Z_order]',
                Delta='original_full_source_quotients.actual_Delta_axial5[Z_order]'),
            leading_source_fixed_input_N0=True,range_endpoints_not_field_coefficients=True,
            original_independent_P0_and_parent_leading_recipes_retained=True)
    providers['O3_power']['quiet_source_proof']=dict(
        source_frontend_binding=current.FRONTEND_BINDINGS['o3'],
        admission=ast_binding(current.o3.source_mu_admission),
        checked_receipt=current.RECEIPT,checked_receipt_sha256=sha(current.RECEIPT),
        source_identity='2*mu>=eta implies q=q_Z=0 throughout current power offset [0,2]',
        original_primitive_A_B_zero_but_incoming_histories_retained=True)
    for chart in ('reference','restoration','postrestore','Rh_reference'):
        providers[chart]['quantity_paths']['b']='original_full_source_quotients.actual_nonzero_b_axial5[Z_order]'
    providers['long_reshape']['quantity_paths'].update(
        b='exact zero from original_full_source_quotients.exact_b_and_t0_zero',
        Delta='original_full_source_quotients.Delta_axial5[Z_order]')
    return providers


def current_range_containment(saved):
    """Bind the current four-stage range proof to these same defining functions.

    The prior checkers already prove each original inverse/C1 density range.
    Here identical incoming/outgoing records, every original partition and
    source receipt connect that induction to the newly assembled exact graph.
    """
    modules=(current.core,current.previous.previous,current.previous,current)
    stage_reports=[];bindings=[]
    for module in modules:
        data=json.loads(gzip.decompress((HERE/module.NAME).read_bytes()))
        receipt=json.loads((HERE/module.RECEIPT).read_bytes())
        require(data[module.GATE] and receipt[module.GATE] and receipt['all_passed'],'Current source stage unaccepted')
        require(data['source_family']==saved['source_family'] and data['fixed_leading_input_N0']==saved['fixed_leading_input_N0'],
            'Changed leading source family or fixed prefix')
        key=next(k for k in data if k.startswith('actual_four_Z_all_N_'))
        stage_reports.append({tuple(row['exact_Z_cell']):row for row in data[key]})
        bindings.append(dict(report=module.NAME,report_sha256=sha(module.NAME),
            receipt=module.RECEIPT,receipt_sha256=sha(module.RECEIPT)))
    count=seams=0;partitions=native_partitions();manifest=recipe_manifest()
    for ends in saved['exact_Z_partition']:
        before=None;charts=[];cell_count=0
        for data in stage_reports:
            row=data[tuple(ends)]
            windows=next(v for k,v in row.items() if k.startswith('actual_all_N_') and k.endswith('_windows'))
            for window in windows:
                if before is None:
                    require(row['source_owned_zero_inlet_proved_for_every_N'] and
                        window['actual_N_scaled_incoming_C0_Z']==row['actual_N_scaled_zero_inlet_C0_Z'],
                        'Original exact zero correction inlet changed')
                if before is not None:
                    require(before==window['actual_N_scaled_incoming_C0_Z'],'Actual normalized predecessor history changed')
                    seams+=1
                before=window['actual_N_scaled_outgoing_C0_Z'];charts.append(window['actual_chart'])
                cells=window['actual_all_N_source_cells'];cell_count+=len(cells)
                expected=current.serialized(list(zip(partitions[window['actual_chart']],partitions[window['actual_chart']][1:])))
                require([[q['exact_left'],q['exact_right']] for q in cells]==expected,'Actual exact native partition changed')
                for cell in cells:
                    require(cell['physical_Jacobian_applied_once'],'Original Jacobian multiplicity changed')
                    require(cell['all_N_phase_cover_not_fixed_N_phase_reuse'],'Phase cover not all-N')
        require(tuple(charts)==CHARTS and cell_count==57,'Actual current source partition changed')
        require(before==row['actual_uniform_N_scaled_Rc_correction_C0_Z'],'Actual Rc correction range changed')
        count+=cell_count
    return dict(original_current_stage_bindings=bindings,current_cells_bound=count,current_window_seams_bound=seams,
        exact_native_partitions=current.serialized(partitions),source_recipe_manifest=manifest,
        exact_density_identity=current.core.density_binding(),
        exact_graph_range_induction=[
            'Same original source recipe quantities, original q cutoff/collar and monotone inverse give A/B C0/Z functions.',
            'Full-phase original inverse ranges contain the exact fractional_part(N*offset) for every N>=N0.',
            'Normalized density equals coefficient[-1]+coefficient[-2]/N, including exprel and all cross terms.',
            'Original positive exp(-rate*(right-offset))*Jacobian integral is bounded by existing cell mass/suffix ranges.',
            'Zero inlet and the identical checked predecessor records inductively contain every exact chart history.',
            'Same terminal positive amplitude/mu and joint quotient rules contain the exact normalized Rc targets.'],
        original_range_receipts_not_used_as_point_values=True,
        new_exact_graph_and_existing_uniform_range_use_same_defining_functions=True)


def original_loop(g,roots,phase,eta,dstar,chart,coordinate,sc):
    """Exact original phase inverse and C1 primitives, with lazy flat branch.

    The formulas are the Z specialization of generic_loop_function_sources.build.
    All roots below are original function recipes, not range magnitudes.
    """
    one,two,half=g.one,g.constant(2),g.constant('1/2')
    mul,add,sub,neg=g.mul,g.add,g.sub,g.neg
    div=lambda a,b:g.quotient(a,b,'same original active positive denominator')
    fun=g.unary
    collar=None
    if chart=='first_micro':
        collar=g.mul(g.constant('3/4'),sc)
    shear_flat=g.node('exact_real_comparison',operator='>=',left=roots['Delta'].value.node,right=eta.node)
    predicate=shear_flat
    if collar is not None:
        collar_flat=g.node('exact_real_comparison',operator='<=',left=coordinate.node,right=collar.node)
        predicate=g.node('logical_or',arguments=[shear_flat.node,collar_flat.node])
    flat=lambda v:g.node('original_lazy_flat_branch',active_body=v.node,Delta=roots['Delta'].value.node,
        eta=eta.node,flat_value=g.zero.node,flat_predicate=predicate.node,
        original_collar_coordinate=coordinate.node if collar is not None else None,
        original_collar_upper=collar.node if collar is not None else None,
        original_collar_extra_flat_condition='first_micro coordinate<=3*sc/4 on domain>=sc/2',
        original_collar_definition_module=PREFIX+'current_original_whole_Z_R100_finite_N.py',
        original_collar_proof_module=Path(current.core.__file__).name,
        original_collar_receipt=current.core.RECEIPT,original_collar_receipt_sha256=sha(current.core.RECEIPT),
        active_body_not_evaluated_on_flat_branch=True)
    a,b,E,p2,t0,Delta=[roots[key] for key in ('a','b','E','p2','t0','Delta')]
    gamma=sub(mul(two,eta),Delta.value);cut=sub(one,div(Delta.value,eta))
    sigma=fun('original_flat_sigma',cut)
    root=fun('positive_sqrt',div(gamma,mul(two,a.value)))
    q=flat(mul(sigma,root))
    sigZ=neg(mul(fun('original_flat_sigma_prime',cut),div(Delta.Z,eta)))
    LZ=neg(mul(half,add(div(Delta.Z,gamma),div(a.Z,a.value))))
    qZ=flat(mul(root,add(sigZ,mul(sigma,LZ))))
    u=div(mul(p2.value,q),dstar);uZ=div(add(mul(p2.Z,q),mul(p2.value,qZ)),dstar)
    hinv=div(one,fun('positive_sqrt',add(one,mul(u,u))))
    r=mul(u,hinv);rZ=mul(uZ,hinv,hinv,hinv)
    alpha=mul(two,q,hinv)
    alphaZ=mul(two,add(mul(qZ,hinv),neg(mul(q,u,uZ,hinv,hinv,hinv))))
    psi_name='loop_angle_'+chart;psi=g.symbol(psi_name)
    nn=sub(fun('cos',psi),r)
    den=add(one,neg(mul(two,r,fun('cos',psi))),mul(r,r))
    w=div(nn,den);wr=div(add(neg(den),mul(two,nn,nn)),mul(den,den))
    t=add(t0.value,mul(alpha,w));tZ=add(t0.Z,mul(alphaZ,w),mul(alpha,wr,rZ))
    pi=g.node('mathematical_pi');twopi=mul(two,pi)
    ratio=add(one,mul(t0.value,t0.value),mul(two,q,q))
    K=div(one,mul(twopi,ratio))
    KZ=neg(div(mul(K,add(mul(two,t0.value,t0.Z),mul(g.constant(4),q,qZ))),ratio))
    integral=lambda v:controls.integral(g,v,psi_name,g.zero,psi,
        measure='original angle dpsi at fixed slow source parameters')
    T1,T2=integral(t),integral(mul(t,t))
    P=add(psi,T2);Phi=mul(K,P)
    inverse=g.node('original_monotone_phase_inverse',phase_function=Phi.node,
        angle_variable=psi_name,target_phase=phase.node,angle_lower=g.zero.node,angle_upper=twopi.node,
        derivative_positive_certificate='K*(1+t^2)>0; a>0 and original q-flat extension',
        phase_endpoints=[0,1],flat_inverse='2*pi*phase',
        Delta=Delta.value.node,eta=eta.node,branch_inverse_is_original_not_support_cap=True)
    at=lambda v:g.node('substitute_original_inverse_angle',body=v.node,inverse_angle=inverse.node,
        angle_variable=psi_name)
    lam=at(mul(K,add(one,mul(t,t))))
    PhiZ=add(mul(KZ,P),mul(K,integral(mul(two,t,tZ))))
    psiZ=neg(div(at(PhiZ),lam))
    chi=sub(phase,div(inverse,twopi));T1hat=at(T1);that=at(t)
    M=sub(neg(div(mul(a.value,T1hat),twopi)),mul(b.value,phase))
    T1hatZ=add(at(integral(tZ)),mul(that,psiZ))
    MZ=sub(neg(div(add(mul(a.Z,T1hat),mul(a.value,T1hatZ)),twopi)),mul(b.Z,phase))
    A=source.C1Function(flat(mul(half,a.value,chi)),
        flat(mul(half,sub(mul(a.Z,chi),div(mul(a.value,psiZ),twopi)))))
    B=source.C1Function(flat(mul(half,E.value,M)),
        flat(mul(half,add(mul(E.Z,M),mul(E.value,MZ)))))
    return A,B,dict(inverse=inverse.node,A=[A.value.node,A.Z.node],B=[B.value.node,B.Z.node],
        q=[q.node,qZ.node],Phi=Phi.node,phase_Z_exact_zero=True)


def parameter_identity(owner):
    providers=(owner.owner.source,owner.owner.switch_source,owner.owner.long_source,owner.owner.patch_source)
    first=providers[0]
    for row in providers:
        require(row.identity==owner.identity and row.N==owner.N0,'Canonical parameter source family changed')
        require(row.eta_log._mpi_==first.eta_log._mpi_ and row.dstar_log._mpi_==first.dstar_log._mpi_,
            'Original chart cutoff/inverse parameters are not identical source tuples')
    return packets.encode(dict(source_family=owner.identity,fixed_leading_input_N0=owner.N0,
        original_canonical_eta_log=first.eta_log,original_canonical_dstar_log=first.dstar_log,
        source_provider_types=[type(row).__name__ for row in providers],
        exact_original_parameter_tuples_equal=True,parameter_ranges_not_source_field_values=True))


def build_graph(family,N0,selection,manifest,parameter_proof):
    g=source.FunctionTransportGraph()
    param=lambda name,definition,**binding:g.node('current_original_source_parameter',name=name,
        definition=definition,source_family=family,fixed_leading_input_N0=N0,
        quantity_is_exact_original_parameter_not_range_endpoint=True,Z_independent=True,
        current_source_receipt=current.RECEIPT,current_source_receipt_sha256=sha(current.RECEIPT),**binding)
    P=g.add(g.unary('exp',g.constant(40)),g.constant(11))
    parameters=dict(logP=P,
        logC=param('logC','same selected_logCstar and CurrentLongRadiusPhase exact singleton tuple'),
        T=param('T','same original T=400*A singleton; all-N long window_length exact tuple guards'),
        hbB=param('hbB','original bridge_owner.flow.h; positive source ledger'),
        hbS=param('hbS','original switch_owner.flow.h; not merged with hbB by name'),
        sc=param('sc','same original selected_first_phase_endpoint exact singleton'))
    mu=g.unary('exp',g.sub(g.unary('log',g.constant('1/1000')),g.mul(g.constant(4),P)))
    parameters['mu']=mu
    N=g.node('exact_positive_integer_expression',**selection.record(),
        one_common_N_for_all_charts_and_control_denominators=True)
    invN=g.quotient(g.one,N,'DoubleDyadicInteger exact positive integer')
    Z=g.symbol('Z')
    offsets,domains,x=exact_geometry()
    require(parameter_proof['source_family']==family and parameter_proof['fixed_leading_input_N0']==N0,
        'Same original canonical parameter proof required')
    eta=param('eta','same canonical original source eta_log; all current owners share identical source tuples',
        canonical_parameter_proof=parameter_proof,source_attribute='eta_log')
    dstar=param('d_star','same canonical original source dstar_log; all current owners share identical source tuples',
        canonical_parameter_proof=parameter_proof,source_attribute='dstar_log')
    H={key:source.C1Function(g.zero,g.zero) for key in RATES};windows=[];leaf_nodes=[]
    partitions=native_partitions()
    for chart in CHARTS:
        variable='native_'+chart;t=g.symbol(variable)
        offset=source.expression(g,offsets[chart],parameters,t)
        lower,upper=[source.expression(g,v,parameters,t) for v in domains[chart]]
        left=source.expression(g,offsets[chart].subs(x,domains[chart][0]),parameters)
        right=source.expression(g,offsets[chart].subs(x,domains[chart][1]),parameters)
        width=g.sub(right,left)
        jac=source.expression(g,sy.diff(offsets[chart],x),parameters,t)
        phase=g.unary('fractional_part',g.mul(N,offset))
        def leaf(quantity,order):
            handle=g.node('current_original_leading_function_recipe',source_family=family,
                native_chart=chart,recipe=manifest[chart],quantity=quantity,Z_order=order,
                coordinate=t.node,Z_variable=Z.node,phase_independent_leading_source=True,
                defining_quantity_not_a_range_value=True)
            leaf_nodes.append(handle.node);return handle
        raw={key:source.C1Function(leaf(key,k),leaf(key,1)) for key,k in
            ((name,0) for name in ('E','V','a','b','p2','t0','Delta'))}
        if chart=='O3_power':
            density={key:source.C1Function(g.zero,g.zero) for key in RATES};loop=None
        else:
            A,B,loop=original_loop(g,raw,phase,eta,dstar,chart,t,parameters['sc'])
            coefficients,U=algebra.coefficient_pairs(g,raw['E'],raw['V'],A,B,N)
            density={key:g.c1add(coefficients[-1][key],g.c1scale(invN,coefficients[-2][key])) for key in RATES}
        before=H;after={};contributions={};memories={};pieces=[]
        for a,b in zip(partitions[chart],partitions[chart][1:]):
            sa,sb=symbolic_coordinate(a),symbolic_coordinate(b)
            pa,pb=[source.expression(g,v,parameters) for v in (sa,sb)]
            pieces.append(dict(exact_left=current.serialized(a),exact_right=current.serialized(b),
                lower=pa.node,upper=pb.node,
                left_offset=source.expression(g,offsets[chart].subs(x,sa),parameters).node,
                right_offset=source.expression(g,offsets[chart].subs(x,sb),parameters).node))
        for key,rate in RATES.items():
            kernel=g.unary('exp',g.neg(g.mul(g.constant(rate),g.sub(right,offset))))
            memory=g.unary('exp',g.neg(g.mul(g.constant(rate),width))) if rate else g.one
            pair=density[key]
            integrate=lambda h:g.zero if h==g.zero else g.add(*(controls.integral(g,g.mul(kernel,h,jac),
                variable,source.FunctionRef(g,p['lower']),source.FunctionRef(g,p['upper']),
                measure='original dlogR: Jacobian exactly once',exact_native_left=p['exact_left'],
                exact_native_right=p['exact_right'],actual_current_chart=chart,
                N_scaled_density=True,Z_derivative_order_held_phase=True) for p in pieces))
            contribution=source.C1Function(integrate(pair.value),integrate(pair.Z))
            after[key]=g.c1add(g.c1scale(memory,before[key]),contribution)
            contributions[key]=contribution;memories[key]=memory
        H=after
        windows.append(dict(chart=chart,domain=[lower.node,upper.node],offset=offset.node,
            left_offset=left.node,right_offset=right.node,width=width.node,Jacobian=jac.node,
            phase=phase.node,source_roots=raw,loop=loop,density=density,actual_current_source_cells=pieces,
            incoming=before,contributions=contributions,memory=memories,outgoing=H,
            quiet_source_proof=manifest[chart].get('quiet_source_proof'),
            quiet_local_source_zero=chart=='O3_power',predecessor_memory_never_reset=True))
    # Terminal amplitude is the actual same power source at coordinate 2.
    terminal_coordinate=g.constant(2)
    terminal=lambda k:g.node('current_original_leading_function_recipe',source_family=family,
        native_chart='O3_power',recipe=manifest['O3_power'],quantity='E',Z_order=k,
        coordinate=terminal_coordinate.node,Z_variable=Z.node,
        defining_quantity_not_a_range_value=True,terminal_positive_A_rc=True,
        original_positive_amplitude_receipt=current.RECEIPT,
        original_positive_amplitude_receipt_sha256=sha(current.RECEIPT))
    amplitude=source.C1Function(terminal(0),terminal(1))
    targets,joint=algebra.normalize_history(g,H,amplitude,mu)
    built=dict(graph=g,N=N,parameters=parameters,amplitude=amplitude,N_scaled_targets=targets,
        history={key:g.c1scale(invN,pair) for key,pair in H.items()},normalized_history=H,
        windows=windows,joint_numerator=joint,source_family=family,source_leaf_nodes=leaf_nodes)
    return controls.exact_control_graph(built,iterations=3)


def encode_graph(value):
    if isinstance(value,source.FunctionRef):return value.node
    if isinstance(value,source.C1Function):return [value.value.node,value.Z.node]
    if isinstance(value,dict):return {key:encode_graph(v) for key,v in value.items() if key!='graph'}
    if isinstance(value,(list,tuple)):return [encode_graph(v) for v in value]
    return value


def load_inputs():
    saved=json.loads(gzip.decompress((HERE/current.NAME).read_bytes()))
    receipt=json.loads((HERE/current.RECEIPT).read_bytes())
    require(saved[current.GATE] and receipt[current.GATE] and receipt['all_passed'],'Checked current Rc family required')
    require(saved['source_family']==receipt['source_family'],'Same current source family required')
    hashes=dict(receipt['input_hashes'])
    hashes[current.NAME]=sha(current.NAME);hashes[current.RECEIPT]=sha(current.RECEIPT)
    modules=(current,algebra,controls,source,repair,source.sources)
    for module in modules:hashes[Path(module.__file__).name]=sha(Path(module.__file__).name)
    hashes[Path(__file__).name]=sha(Path(__file__).name)
    for name,digest in hashes.items():require(sha(name)==digest,'Changed current prerequisite: '+name)
    return saved,hashes


def frequency_connection(saved,c):
    rows=saved['actual_four_Z_all_N_outer_Rc_transports'];N0=saved['fixed_leading_input_N0']
    logmu=packets.interval(c,rows[0]['actual_positive_mu_log'])
    for row in rows:
        require(packets.interval(c,row['actual_positive_mu_log'])._mpi_==logmu._mpi_,'Same source mu required')
        require(row['target_C1_caps']['uniform_for_all_integer_N_ge_N0'],'Uniform original target function bounds required')
    defining=c.ln(c.mpf('1/1000'))-4*(c.exp(40)+11)
    require(defining._mpi_==logmu._mpi_,'Exact original source mu definition required')
    mu=c.exp(logmu);W=repair.fresh_weights(c,mu,cells=512);B=repair.fresh_linear_inverse(c,mu,W)
    bounds={}
    for key in repair.ROWS:
        caps=[row['target_C1_caps']['transformed_N_scaled_target_C1_caps'][key] for row in rows]
        logs=[packets.interval(c,q['log_absolute_upper']) for q in caps if not q['exact_zero']]
        local=repair.LogUpper(c,None if not logs else c.mpf(max(ep(q)[1] for q in logs)))
        # Global sup|f| and sup|f_Z| may occur in different Z cells.
        # Twice the max of the local sum bounds the original global C1 norm.
        bounds[key]=local*repair.LogUpper.constant(c,2)
    nonzero=[q.log for q in bounds.values() if q.log is not None]
    D=repair.LogUpper(c,None if not nonzero else c.mpf(max(ep(q)[1] for q in nonzero)))
    target=dict(transformed_N_scaled_target_C1_caps={k:v.record() for k,v in bounds.items()},whole_target_C1_cap=D.record())
    if D.log is None:
        # No log(0); zero actual target has unique zero control in a small ball.
        raise ValueError('Current nonzero family expected; zero-target contract requires a separate zero graph')
    cond=repair.contraction_log_conditions(c,target,B,W,c.mpf(ep(logmu)[0]),c.ln(N0))
    required=c.mpf(ep(cond['repair_sufficient_common_log_N_lower'])[1])
    chosen=choose_integer(c,N0,required);selection=chosen.prove(c,N0,required)
    read=lambda key:packets.interval(c,cond[key]['log_absolute_upper'])
    CA,CQ,rho,Dlog=[read(key) for key in ('exact_integral_matrix_inverse_log_cap',
        'general_transformed_quadratic_C1_log_cap','formal_control_C1_ball_radius_log','actual_target_C1_log_cap')]
    logN=chosen.log_interval(c)
    logL=c.mpf(ep(c.ln(2)+CA+CQ+rho-logN)[1])
    require(ep(logL+c.ln(2))[1]<=0,'Original C1 map fails half contraction')
    # The original sufficient threshold and radius definitions prove image/positivity.
    require(ep(logN-required)[0]>=0,'Repair image/positivity conditions fail')
    limit=dict(exact_definition='h*=lim_C1 h_n; h0=0; h_next=-B(mu)^-1*(d_exact(N,Z)+Q(mu,h)/N)',
        target='17-chart current exact source integral functions normalized directly as N*r',
        global_Z_domain=[-1,1],whole_Z_target_cap=target,
        uniform_control_ball_log=rho,uniform_contraction_log_upper=logL,contraction_at_most='1/2',
        exact_C1_limit_exists_unique_for_this_repair_only_integer=True,
        limit_residual_identity='B*h*+d_exact+Q(h*)/N=0 in C1(Z)',
        implicit_Z_identity='(B+DQ(h*)/N)*h*_Z=-d_exact_Z',
        tail='norm_C1(h*-h_n)<=rho*2^-n',
        preconditioned_residual_tail='norm_C1(B^-1*(B*h_n+d+Q(h_n)/N))<=3*rho*2^(-n-1)',
        tail_radius_log=rho,tail_factors_must_stay_separate_from_astronomical_log=True,
        numerical_quadrature_or_phase_errors_not_covered=True,
        global_higher_derivative_source_stress_and_heat_conditions_not_proved=True)
    return chosen,dict(actual_source_log_mu=logmu,actual_source_mu=mu,
        original_whole_Z_target_caps=target,
        global_C1_cap_aggregation='2*max_cell(sup_cell|d|+sup_cell|d_Z|); extrema may lie in different cells',
        fresh_exact_integral_weights=W,
        fresh_linear_inverse=B,repair_conditions=cond,selected_integer=selection,C1_limit=limit)


def factored_tail(limit,depth):
    if type(depth) is not int or not 0<=depth<=4096:raise ValueError('Exact finite depth in[0,4096] required')
    return dict(depth=depth,radius_log=limit['tail_radius_log'],
        relative_control_tail=dict(numerator=1,denominator=2**depth),
        relative_preconditioned_residual_tail=dict(numerator=3,denominator=2**(depth+1)),
        exact_mathematical_C1_tail_only=True)


def run():
    began=time.monotonic();saved,hashes=load_inputs()
    c=MPIntervalContext();c.dps=500
    with mp.workdps(540):
        containment=current_range_containment(saved)
        chosen,frequency=frequency_connection(saved,c)
        manifest=recipe_manifest()
        parameter_proof=parameter_identity(current.WholeZAllNOuterRcFunctions())
        for row in manifest.values():hashes[row['module']]=row['module_sha256']
        graph=build_graph(saved['source_family'],saved['fixed_leading_input_N0'],chosen,manifest,parameter_proof)
        result=dict(**{GATE:True},source_family=saved['source_family'],
            fixed_leading_input_N0=saved['fixed_leading_input_N0'],exact_Z_partition=saved['exact_Z_partition'],
            exact_current_charts=CHARTS,current_range_cells_per_Z=57,current_range_cells_whole_Z=228,
            exact_source_recipe_manifest=manifest,original_generic_loop_binding=ast_binding(source.sources.build),
            exact_current_source_uniform_range_containment=containment,
            canonical_original_source_parameter_proof=parameter_proof,
            exact_function_graph_nodes=graph['graph'].nodes,exact_current_integral_and_control_graph=encode_graph(graph),
            actual_whole_Z_frequency_connection=frequency,
            exact_C1_tail_examples=[factored_tail(frequency['C1_limit'],n) for n in (0,3,32)],
            exact_function_graph_not_numerical_point_oracle=True,
            actual_leading_source_ranges_not_used_as_function_values=True,
            actual_repair_only_integer_expression_selected=True,
            actual_global_frequency_admitted=False,actual_numeric_controls_evaluated=False,
            actual_terminal_Z_function_closure_installed=False,actual_power_repair_band_source_extended=False,
            physical_original_exterior_five_targets_closed=False,**dict.fromkeys(current.OPEN,False),
            input_hashes=hashes,execution_seconds=time.monotonic()-began,
            scope='Current 17-chart exact source recipe/phase-inverse/integral graph, direct N-scaled relative '
                'targets, exact repair-only integer descriptor and source-bound C1 control limit. '
                'No numeric source oracle/controls, actual repair-band join, global frequency or full reconstruction admission.')
        (HERE/NAME).write_bytes(gzip.compress((json.dumps(packets.encode(result),separators=(',',':'))+'\n').encode(),mtime=0))
    print('CURRENT_17_CHART_CONTROL_CONNECTION',len(graph['graph'].nodes),'nodes; exact J',chosen.J,flush=True)
    return result


if __name__=='__main__':run()
