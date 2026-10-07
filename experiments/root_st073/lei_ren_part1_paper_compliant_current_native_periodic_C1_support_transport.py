"""Original whole-period Z derivative supports before nonlinear density.

These derivative bounds use the original implicit phase and parameter
identities. They are not derivatives of C0 caps and do not define a field.
"""
import ast
import inspect
import json
from pathlib import Path
import textwrap
import time
import lei_ren_part1_paper_compliant_current_native_collected_q2_transport as accepted
import lei_ren_part1_paper_compliant_current_native_serial_C1_cell as serial

cutoff=accepted.cutoff;cover=accepted.cover;density=accepted.density;native=accepted.native;packets=accepted.packets
HERE,PREFIX,sha=accepted.HERE,accepted.PREFIX,accepted.sha;ep=accepted.ep;ZERO=accepted.ZERO;DZ=(0,1)
NAME=PREFIX+'current_native_periodic_C1_support_transport.json'
RECEIPT=PREFIX+'current_native_periodic_C1_support_transport_check.json'
GATE='current_original_whole_period_Z_support_before_fixed_N_transport_executed'


def periodic_C1_support(roots,qrows,q2rows,log_a_lower,dstar_log):
    """The accepted signed/small Poisson estimates on one original cutoff branch."""
    c=roots['a'][ZERO].ctx;scalar=roots['a'][ZERO].scalar
    theorem=serial.periodic_parameter_theorem()
    if qrows[ZERO].zero and qrows[DZ].zero:
        if any(not value.zero for value in q2rows.values()):raise ValueError('Original flat q squared consistency required')
        zero=scalar(0)
        return dict(values=dict(A_Z=zero,B_Z_over_Pstar=zero),record=dict(
            original_periodic_parameter_C1_theorem=theorem,original_q_and_q_Z_flat_proof=True,
            original_entire_fractional_phase_cover=c.mpf((0,1))))
    # The active original source has a<=kappa<2+eta<=5/2 and
    # |b|<=v/2<=3/2. Restrict only their C0 outer ranges; a_Z,b_Z stay original.
    a,ar=accepted.accepted.signed_C0_support_range(roots['a'][ZERO],scalar(c.mpf((0,'2.5'))))
    b,br=accepted.accepted.signed_C0_support_range(roots['b'][ZERO],scalar(c.mpf(('-1.5','1.5'))))
    t0=(-b).positive_divide(a,log_a_lower)
    t0_Z=(-roots['b'][DZ]-t0*roots['a'][DZ]).positive_divide(a,log_a_lower)
    absolute=serial.absolute_upper;minimum=serial.minimum_upper
    T=minimum(t0,roots['t0'][ZERO]);TZ=minimum(t0_Z,roots['t0'][DZ])
    q,q_Z=qrows[ZERO],qrows[DZ]
    dstar=accepted.prior.ScaledEnclosure(accepted.prior.FormalScale(a.scale.bases,offset=dstar_log),1,a.ledger)
    u_Z=(roots['p2'][DZ]*q+roots['p2'][ZERO]*q_Z).positive_divide(dstar,dstar_log)
    Q,QZ,AU,AZ,UZ,R,RZ=[absolute(value) for value in
        (q,q_Z,a,roots['a'][DZ],u_Z,q2rows[ZERO],q2rows[DZ])]
    pi=c.pi
    A_cap=minimum(scalar(c.mpf('1.25')),AU*(T*Q*4+R*101))
    B_cap=minimum(scalar(c.mpf('1.5')),T*A_cap+AU*Q*2)
    # The accepted |W1|, |W1_u| and |(sW2)_u| bounds give the
    # original fixed-angle T2_Z numerator. Collect q*q_Z=R_Z/2
    # before taking absolute ranges. Linear q_Z remains genuinely present.
    firstZ=(TZ*Q+T*QZ+T*Q*UZ)*(4*pi)+T*Q*UZ*(16*pi)
    nuZ=T*TZ*2+RZ*2
    T2Z=T*TZ*(4*pi)+firstZ*4+RZ*(400*pi)+R*UZ*(4000*pi)
    L=absolute(nuZ*(2*pi)+T2Z)
    AZ_cap=absolute(AZ*c.mpf('.5')+AU*L*(1/(4*pi)))
    E=roots['E'][ZERO];EZ=absolute(roots['E'][DZ]);EU=absolute(E)
    BZ_cap=absolute(EZ*B_cap+EU*(TZ*A_cap+T*AZ_cap+AZ*Q*2+AU*QZ*2
        +AU*Q*UZ*10+AU*(T+1)*L*(1/(4*pi))))
    values=dict(A_Z=serial.symmetric_bound(AZ_cap),B_Z_over_Pstar=serial.symmetric_bound(BZ_cap))
    record=dict(original_periodic_parameter_C1_theorem=theorem,
        original_entire_fractional_phase_cover=c.mpf((0,1)),original_active_a_b_C0_support=dict(a=ar,b=br),
        original_A_and_B_over_E_C0_supports=dict(A=A_cap.record(),B_over_E=B_cap.record()),
        original_fixed_angle_T2_Z_upper=T2Z.record(),original_nu_Z_upper=nuZ.record(),
        original_implicit_phase_Z_numerator_upper=L.record(),
        original_whole_period_Z_derivative_supports={key:value.record() for key,value in values.items()},
        q_times_q_Z_collected_as_direct_q_squared_Z_over2=True,
        genuine_linear_q_Z_and_u_Z_terms_retained=True,
        original_a_Z_b_Z_p2_Z_E_Z_rows_not_clipped_or_erased=True,
        derivative_bounds_follow_original_implicit_phase_not_C0_cap_differentiation=True,
        signed_small_Poisson_parameter_bounds_uniform=True,
        support_ranges_not_defining_field_values=True)
    return dict(values=values,record=record)


def compile_periodic_query(source_query,kernel,log_a_lower,support_rows):
    cache={}
    def supported_first_jets(source,qrows,q2rows,dstar,phi,branch):
        got=accepted.compile_collected_first_jets(source,qrows,q2rows,dstar,phi,branch)
        if got['values'] is None:return got
        key=id(source)
        if key not in cache:cache[key]=periodic_C1_support(source['roots'],qrows,q2rows,log_a_lower,dstar)
        support=cache[key];original=got['values'];changed=dict(original);decisions={}
        for name in ('A_Z','B_Z_over_Pstar'):
            changed[name],decisions[name]=accepted.accepted.signed_C0_support_range(original[name],support['values'][name])
        if changed['A'] is not original['A'] or changed['B_over_Pstar'] is not original['B_over_Pstar']:
            raise ValueError('Derivative support must retain original C0 primitive objects')
        trace=dict(original_whole_period_derivative_support=support['record'],ordinary_Z_range_decisions=decisions,
            original_C0_A_B_objects_retained=True,original_linear_q_Z_and_all_source_cross_rows_retained=True,
            support_is_outer_derivative_range_not_new_source_value=True)
        support_rows.append(trace);got['values']=changed
        got['record'].update(original_first_jet_rows_before_periodic_Z_support={key:value.record() for key,value in original.items()},
            original_A_B_first_derivative_enclosures={key:value.record() for key,value in changed.items()},
            original_whole_period_Z_support_before_density=trace)
        return got
    tree=ast.parse(textwrap.dedent(inspect.getsource(cover.NativeSignedUDensityCover.spatial_query)))
    tree.body[0].name='periodic_supported_original_spatial_query'
    tree,counts=density.replace_expressions(tree,[
        ("self.first_owner.owner.query(chart,Z,geometry['raw']['coordinate'])","conditional_source_query(chart,Z,geometry['raw']['coordinate'])"),
        ("branch_first_jets(source,qsource['rows'],self.first_owner.dstar,phi,branch)","supported_first_jets(source,qsource['rows'],qsource['q2_rows'],self.first_owner.dstar,phi,branch)"),
        ("density.density_Z_kernels(E,E_Z,V,V_Z,got['values'],N)","supported_density_kernels(E,E_Z,V,V_Z,got['values'],N)")])
    scope=dict(vars(cover));scope.update(conditional_source_query=source_query,supported_density_kernels=kernel,
        supported_first_jets=supported_first_jets)
    exec(compile(tree,'<original-spatial-density-with-whole-period-Z-support>','exec'),scope)
    return scope['periodic_supported_original_spatial_query'],counts


class NativePeriodicC1Density(accepted.NativeCollectedQ2Density):
    def __init__(self,owner):
        super().__init__(owner);self.periodic_support_rows=[]
        self.hashes={**self.hashes,Path(serial.__file__).name:sha(Path(serial.__file__).name),
            Path(__file__).name:sha(Path(__file__).name)};self.service.bind_hashes(self.hashes)

    @native.inlet.source_precision
    def spatial_query(self,chart,Z,coordinate,N):
        root_owner=self.qcover.owner.owner.owner
        positive=root_owner.decode(root_owner.inventory[chart]['actual_positive_denominator_theorem'])
        tree=ast.parse(textwrap.dedent(inspect.getsource(density.NativeCutoffDensityCover.spatial_query)))
        tree.body[0].name='periodic_supported_cutoff_spatial_query'
        tree,counts=density.replace_expressions(tree,[
            ('compile_conditional_query(source_query,self.kernel)','compile_original_periodic_query(source_query,self.kernel)')])
        scope=dict(vars(density));scope['compile_original_periodic_query']=lambda source,kernel:compile_periodic_query(
            source,kernel,positive['log_actual_a_positive_lower'],self.periodic_support_rows)
        exec(compile(tree,'<original-cutoff-query-with-whole-period-Z-support>','exec'),scope)
        start=len(self.support_rows);pstart=len(self.periodic_support_rows)
        got=scope['periodic_supported_cutoff_spatial_query'](self,chart,Z,coordinate,N)
        rows=self.support_rows[start:];got['record'].update(
            original_C0_A_B_source_support_ranges_before_nonlinear_density=rows,
            original_whole_period_Z_support_ranges_before_nonlinear_density=self.periodic_support_rows[pstart:],
            derivative_support_not_C0_cap_differentiation=True,original_spatial_query_AST_replacements=counts)
        got['record']['original_A_over_N_formal_factor_and_A_Z_preserved']=all(row['A']['original_tighter_formal_range_retained'] for row in rows)
        return got


class NativePeriodicC1Oracle(accepted.NativeCollectedQ2Oracle):
    def __init__(self,role_owner,built=None):
        super().__init__(role_owner,built);self.owner=NativePeriodicC1Density(self.original_density_owner)
        self.hashes={**self.hashes,**self.owner.hashes};self.service.bind_hashes(self.hashes)

    @native.inlet.source_precision
    def density_frame(self,*,chart,Z,coordinate,N):
        frame=super().density_frame(chart=chart,Z=Z,coordinate=coordinate,N=N)
        frame.record['phase_solver_backend']='original branch-local q squared first jets with original whole-period ordinary Z support'
        return frame


class NativePeriodicC1Transport(accepted.NativeCollectedQ2Transport):
    def __init__(self,role_owner):
        super().__init__(role_owner)
        checked=json.loads((HERE/accepted.RECEIPT).read_bytes())
        if not checked.get('all_passed') or not checked.get(accepted.GATE) or checked['source_family']!=self.family:
            raise ValueError('Checked original collected q squared baseline required')
        self.oracle=NativePeriodicC1Oracle(role_owner,self.built)
        self.hashes={**self.hashes,**self.oracle.hashes,**checked['input_hashes'],accepted.RECEIPT:sha(accepted.RECEIPT),
            Path(__file__).name:sha(Path(__file__).name)};self.service.bind_hashes(self.hashes)


@native.inlet.source_precision
def run(role_owner,*,return_live=False):
    began=time.monotonic();owner=NativePeriodicC1Transport(role_owner);got=owner.route()
    if got['history'] is None or len(got['cells'])!=24:raise ArithmeticError('All24 original continuous cells required')
    baseline=json.loads((HERE/accepted.NAME).read_bytes());comparison=accepted.accepted.compare_targets(owner.ctx,got['record'],baseline)
    result=dict(**got['record'],**{GATE:True},original_whole_period_parameter_theorem=serial.periodic_parameter_theorem(),
        comparison_with_checked_collected_q_squared_target_ranges=comparison,
        strict_target_absolute_upper_reductions=sum(row['strict_absolute_upper_reduction'] for row in comparison.values()),
        derivative_bounds_not_obtained_by_differentiating_C0_caps=True,
        useful_repair_contraction_or_actual_controls_established=False,
        execution_seconds=time.monotonic()-began,input_hashes=owner.hashes,
        scope='Original branch-local whole-period ordinary Z primitive support derived from accepted implicit phase/Poisson parameter identities and direct q squared jets, before nonlinear density and complete24-cell fixed-N original transport. No actual controls/global N/closure/recursion admission.')
    (HERE/NAME).write_text(json.dumps(packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original whole-period Z supports transported through24 cells;strict target upper reductions:',result['strict_target_absolute_upper_reductions'],flush=True)
    return (result,owner,got) if return_live else result
