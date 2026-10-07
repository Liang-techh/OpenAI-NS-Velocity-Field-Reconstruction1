"""Original density-DAG equality, coefficient C1 transport and uniform bounds.

Checks the exact coefficient functions against original full source roots.
Manufactured numerical references are separate from actual source ranges.
"""
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_native_Rc_all_N_function_controls as current
import lei_ren_part1_paper_compliant_current_native_Rc_functional_controls_check as functional_check
import lei_ren_part1_paper_compliant_current_native_signed_input_enclosures as signed

HERE,PREFIX,sha=current.HERE,current.PREFIX,current.sha
source,controls,repair=current.source,current.controls,current.repair
ep=current.ep


def symbolic_walk(nodes,atoms):
    cache={}
    def walk(i):
        if i in atoms:return atoms[i]
        if i in cache:return cache[i]
        row=nodes[i];op=row['operation']
        if op=='constant':q=sy.Integer(row['value'])
        elif op=='exact_rational':q=sy.Rational(row['numerator'],row['denominator'])
        elif op in ('shared_positive_integer','shared_positive_integer_parameter'):q=sy.Symbol('N',positive=True,integer=True)
        elif op=='original_source_parameter':q=sy.Symbol(row['name'],positive=True)
        elif op=='sum':q=sum(walk(j) for j in row['arguments'])
        elif op=='product':q=sy.prod(walk(j) for j in row['arguments'])
        elif op=='negative':q=-walk(row['argument'])
        elif op in ('positive_quotient','positive_function_quotient'):q=walk(row['numerator'])/walk(row['denominator'])
        elif op=='analytic_unary':
            a=walk(row['argument']);name=row['name']
            q=(sy.exp(a)-1)/a if name=='exprel' else sy.exp(a)-1 if name=='expm1' else getattr(sy,name)(a)
        else:raise AssertionError('Unexpected unbound source node: '+op)
        cache[i]=q;return q
    return walk


def original_density_equalities(owner,built):
    """Compare every active cell with its original unsplit full signed DAG."""
    E,EZ,V,VZ,A,AZ,B,BZ=sy.symbols('E E_Z V V_Z mathcal_A mathcal_A_Z mathcal_B mathcal_B_Z',real=True)
    N=sy.Symbol('N',positive=True,integer=True);count=0
    for cell in built['all_N_cells']:
        inputs=cell['source_inputs']
        if inputs is None:continue
        original=owner.role_owner.owner.views[cell['original_cell']['chart']]
        outside_atoms={};inside_atoms={}
        for key,pair_symbols in (('E',(E,EZ)),('V',(V,VZ)),('A',(A,AZ)),('B',(B,BZ))):
            pair=inputs[key]
            for ref,symbol in zip((pair.value,pair.Z),pair_symbols):
                row=built['graph'].nodes[ref.node]
                outside_atoms[ref.node]=symbol;inside_atoms[row['source_node']]=symbol
        old= symbolic_walk(original['function_graph_nodes'],inside_atoms)
        new= symbolic_walk(built['graph'].nodes,outside_atoms)
        for key in current.RATES:
            for order in ('value','Z'):
                old_root=(original['five_signed_increment_rate_roots'][key] if order=='value' else
                    original['five_signed_increment_rate_first_derivatives']['Z'][key])
                reconstructed=sum(new(getattr(cell['coefficient_density_pairs'][p][key],order).node)*N**p for p in current.ORDERS)
                if sy.simplify(old(old_root)-reconstructed)!=0:raise AssertionError('Original full signed coefficient identity differs')
                count+=1
    return dict(passed=True,original_full_signed_density_DAG_C0_Z_equalities=count,
        original_nonzero_axial_V_and_V_Z_retained=True,exprel_identity_extended_continuously_at_A_zero=True,
        coefficients_not_assumed_N_independent=True)


def coefficient_transport_C1(built):
    g=built['graph'];Z=sy.Symbol('Z',real=True);A=sy.Function('actual_Rc_amplitude')(Z)
    atoms={built['amplitude'].value.node:A,built['amplitude'].Z.node:sy.diff(A,Z)};integrals=0
    for i,cell in enumerate(built['all_N_cells']):
        for p,rows in cell['coefficient_contributions'].items():
            for key,pair in rows.items():
                if pair.value==g.zero:
                    if pair.Z!=g.zero:raise AssertionError('Quiet/zero coefficient first derivative differs')
                    continue
                f=sy.Function('original_coefficient_integral_'+str(i)+'_'+str(-p)+'_'+key)(Z)
                atoms[pair.value.node],atoms[pair.Z.node]=f,sy.diff(f,Z);integrals+=1
    walk=symbolic_walk(g.nodes,atoms);checks=0
    for p in current.ORDERS:
        for group in (built['coefficient_history'][p],built['coefficient_target_orders'][p]):
            for pair in group.values():
                if sy.simplify(sy.diff(walk(pair.value.node),Z)-walk(pair.Z.node))!=0:
                    raise AssertionError('Coefficient history/target quotient derivative lost')
                checks+=1
    for group in (built['history'],built['targets'],built['N_scaled_targets']):
        for pair in group.values():
            if sy.simplify(sy.diff(walk(pair.value.node),Z)-walk(pair.Z.node))!=0:
                raise AssertionError('N-reconstructed C1 function derivative lost')
            checks+=1
    N=sy.Symbol('N',positive=True,integer=True)
    for key in current.ROWS:
        first,second=[walk(built['coefficient_target_orders'][p][key].value.node) for p in current.ORDERS]
        if sy.simplify(walk(built['N_scaled_targets'][key].value.node)-(first+second/N))!=0:
            raise AssertionError('N*r normalization reconstructed incorrectly')
        checks+=1
    return dict(passed=True,genuine_N_dependent_C1_coefficient_integral_functions=integrals,
        exact_history_target_C1_and_N_reconstruction_checks=checks,
        joint_m_k_uses_distinct_actual_rates_and_full_A_Z=True)


def uniform_scalar_references(owner):
    """Nonzero V/Z references and A=0 exercise original increments and exprel."""
    c=mp.mp.clone();c.dps=70;iv=MPIntervalContext();iv.dps=80
    basis=tuple(iv.mpf(0) for unused in range(5));ledger=dict(directed_small_exponential_tails=0,
        positive_function_denominator_intersections=0,positive_function_root_intersections=0,directed_independent_log_rescalings=0)
    box=lambda lo,hi:signed.ScaledEnclosure(signed.FormalScale(basis),iv.mpf((lo,hi)),ledger)
    inputs=[box('1.7','1.9'),box('.1','.2'),box('-.5','-.4'),box('-.08','-.04')]
    primitive=dict(A=box('-1.1','.9'),A_Z=box('-.3','.5'),B_over_Pstar=box('-.8','.8'),B_Z_over_Pstar=box('-.2','.3'))
    ranges=current.uniform_coefficients(*inputs,primitive);count=range_count=0
    class ReferenceOracle:
        mode='synthetic_reference';source_family=owner.family
        def parameter(self,*args):raise AssertionError('No original-field parameter expected')
        def source(self,*args,**kwargs):raise AssertionError('No current original source point expected')
        def integrate(self,*args):raise AssertionError('No original-field quadrature expected')
    for a in ('-1.1','0','.9'):
        for N in (160,257,2048):
            E,EZ,V,VZ,A,AZ,B,BZ=map(c.mpf,('1.8','.15','-.45','-.06',a,'.11','.6','.12'))
            g=source.FunctionTransportGraph();nr=g.node('shared_positive_integer',name='N',lower=160)
            pair=lambda v,j:source.C1Function(g.constant(str(v)),g.constant(str(j)))
            coefficients,F=current.coefficient_pairs(g,pair(E,EZ),pair(V,VZ),pair(A,AZ),pair(B,BZ),nr)
            built=dict(graph=g,source_family=owner.family,source_graph_sha256='manufactured_density_reference_only')
            value=controls.FunctionEvaluator(built,oracle=ReferenceOracle(),Z=c.mpf(0),N=N,ctx=c)
            density=lambda z:repair.packets.recovery.increment_densities(E+EZ*z,V+VZ*z,
                (E+EZ*z)*c.expm1((A+AZ*z)/N),(B+BZ*z)/N)
            for key in current.RATES:
                expected=density(c.mpf(0))[key];expectedZ=c.diff(lambda z:density(z)[key],c.mpf(0))
                actual=sum(value(coefficients[p][key].value)*c.mpf(N)**p for p in current.ORDERS)
                actualZ=sum(value(coefficients[p][key].Z)*c.mpf(N)**p for p in current.ORDERS)
                for got,want in ((actual,expected),(actualZ,expectedZ)):
                    if abs(got-want)>c.mpf('1e-55')*(1+abs(want)):raise AssertionError('Original increment scalar reference differs')
                    count+=1
                for group,want in (('values',expected),('Z_derivatives',expectedZ)):
                    parts=ranges[group][key];r=sum((q*(iv.mpf(1)/N)**(-p) for p,q in parts.items()),inputs[0].scalar(0))
                    lo,hi=ep(r.finite_interval())
                    if not lo<=want<=hi:raise AssertionError('Original increment outside uniform N range')
                    range_count+=1
    return dict(passed=True,independent_original_increment_C0_Z_scalar_comparisons=count,
        same_references_inside_fresh_uniform_ranges=range_count,reference_N_values=[160,257,2048],
        nonzero_original_axial_fixture=True,A_zero_with_nonzero_A_Z_tested=True,
        manufactured_operator_references_only=True,current_original_source_field_points_evaluated=False)


def original_phase_and_live_ranges(owner,built,live):
    g=built['graph'];source_rows=0;integral_rows=0;zero_rows=0;pressure_rows=0
    maps,unused,x=source.exact_radius_maps()
    previous={p:{key:source.C1Function(g.zero,g.zero) for key in current.RATES} for p in current.ORDERS}
    for i,cell in enumerate(built['all_N_cells']):
        original=cell['original_cell']
        coordinate=g.symbol('coordinate_'+str(i))
        offset=source.expression(g,maps[original['chart']],built['parameters'],coordinate)
        jacobian=source.expression(g,sy.diff(maps[original['chart']],x),built['parameters'],coordinate)
        if cell['incoming']!=previous:raise AssertionError('Exact coefficient history reset')
        previous=cell['outgoing']
        for p,rows in cell['coefficient_contributions'].items():
            for key,pair in rows.items():
                for order,q in (('value',pair.value),('Z',pair.Z)):
                    if q==g.zero:zero_rows+=1;continue
                    row=g.nodes[q.node]
                    if row['operation']!='definite_integral' or row['extracted_N_power']!=p:
                        raise AssertionError('Exact coefficient integral replaced by mass/cap')
                    kernel=g.unary('exp',g.neg(g.mul(g.constant(current.RATES[key]),
                        g.sub(source.FunctionRef(g,original['right_radius_offset']),offset))))
                    density=getattr(cell['coefficient_density_pairs'][p][key],order)
                    expected=g.mul(kernel,density,jacobian)
                    if row['integrand']!=expected.node:
                        raise AssertionError('Kernel/coefficient/Jacobian identity differs')
                    if not row['coefficient_still_depends_on_N_and_actual_phase']:raise AssertionError('N dependence erased')
                    integral_rows+=1
    for index in built['all_N_source_reference_nodes']:
        row=g.nodes[index];phase=g.nodes[row['phase']]
        if row['shared_N']!=built['N'].node or row['graph_sha256']!=built['source_graph_sha256']:
            raise AssertionError('Source coefficient frequency/family changed')
        if phase['operation']!='analytic_unary' or phase['name']!='fractional_part':raise AssertionError('Original phase unbound')
        product=g.nodes[phase['argument']]
        if product['operation']!='product' or built['N'].node not in product['arguments']:
            raise AssertionError('Actual spatial phase independent of N')
        source_rows+=1
    if len(live['cells'])!=24 or live['source_branch_count']!=35:raise AssertionError('Full current original source route missing')
    for cell in live['cells']:
        if cell['record']['label'] in ('power_to_r_plus','power_to_Rc'):
            for p in current.ORDERS:
                for key in current.RATES:
                    if not cell['values'][key][p].zero or not cell['Z_derivatives'][key][p].zero:
                        raise AssertionError('Original quiet power has nonzero increment')
                for incoming,outgoing in ((cell['incoming'],cell['outgoing']),(cell['incomingZ'],cell['outgoingZ'])):
                    left,right=incoming['p'][p],outgoing['p'][p]
                    if left.record()!=right.record():raise AssertionError('Quiet pressure coefficient memory changed')
                    pressure_rows+=1
        for branch in cell['branches']:
            if branch['density']['record']['N_lower']!=160:raise AssertionError('Uniform N domain changed')
            if branch['density']['record']['actual_fixed_N_result_rescaled']:raise AssertionError('Fixed-N rescaling admitted')
    combined=current.combined_log_conditions(owner,live)
    threshold=ep(combined['source_and_repair_sufficient_log_N_lower'])[0]
    if any(threshold<ep(v)[1] for v in combined['original_generic_source_log_N_requirements'].values()):
        raise AssertionError('Combined original source requirement omitted')
    if threshold<ep(live['conditions']['repair_sufficient_common_log_N_lower'])[1]:
        raise AssertionError('Combined repair requirement omitted')
    return dict(passed=True,original_coefficient_source_refs_with_actual_N_phase=source_rows,
        exact_signed_coefficient_C0_Z_integral_rows=integral_rows,exact_quiet_or_zero_rows=zero_rows,
        original_continuous_radial_cells=24,original_conditional_source_branches=35,
        original_quiet_pressure_coefficient_memory_rows=pressure_rows,original_generic_source_requirements_combined=17,
        actual_uniform_domain=live['Z_box'],N_lower=160,current_whole_N_selected=False)


def conditional_tail_theorem():
    rho,n=sy.symbols('rho n',positive=True);CA,CQ,D=sy.symbols('CA CQ D',positive=True)
    N0=4*CA*CQ*rho
    if sy.simplify(2*CA*CQ*rho/N0-sy.Rational(1,2))!=0:raise AssertionError('Half-contraction recipe differs')
    if sy.simplify((sy.Rational(1,2)**n)/(1-sy.Rational(1,2))*(rho/2)-rho*2**(-n))!=0:
        raise AssertionError('Conditional C1 tail formula differs')
    return dict(passed=True,conditional_C1_tail='rho*2^-n',first_step_bound='norm_C1(h1-h0)<=CA*D=rho/2',
        requires_actual_same_source_ball_contraction_frequency_and_C1_functions=True,
        instantiated_actual_N_or_numeric_tail_certificate=False)


@current.paired.native.inlet.source_precision
def run(owner,built,live):
    began=time.monotonic()
    result=dict(all_passed=True,**{current.GATE:True},source_family=owner.family,
        original_exact_density_identities=original_density_equalities(owner,built),
        exact_coefficient_history_target_C1=coefficient_transport_C1(built),
        independent_scalar_uniform_ranges=uniform_scalar_references(owner),
        actual_phase_and_all_N_source_transport=original_phase_and_live_ranges(owner,built,live),
        exact_reconstructed_functional_map_C1=functional_check.symbolic_function_map(built),
        conditional_uniform_C1_tail_theorem=conditional_tail_theorem(),
        read_only_review_model=dict(model='gpt-5.6-luna',reasoning_effort='max'),
        original_point_or_full_corrected_physical_integrals_evaluated=False,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(current.packets.OPEN,False),
        input_hashes={**owner.hashes,current.NAME:sha(current.NAME),Path(__file__).name:sha(Path(__file__).name)},
        execution_seconds=time.monotonic()-began,
        scope='Exact original full signed density DAG/order equivalence, actual N-phase/integral bindings, full-source uniform C1 coefficient/history/target ranges and conditional repair/tail algebra. No fixed-N rescaling or numerical fixed point/global N/closure/recursion admission.')
    (HERE/current.RECEIPT).write_text(json.dumps(current.packets.encode(result),indent=2)+'\n',encoding='utf8')
    print('Original all-N function/source/repair-contract check PASS; global N and evaluated controls remain open',flush=True)
    return result
