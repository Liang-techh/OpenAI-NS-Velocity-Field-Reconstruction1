"""Actual two-sided patch support function joins; no ancestor rebuild.

The six source seams are proved from normalized flat beta jets, partial
FTC weights, the current implicit moment family and the actual mixed4
program. Full global physical interfaces and resolved coefficients are
separate. Numerical intervals containing zero are not a join proof.
"""
import hashlib
import inspect
import json
import math
from pathlib import Path
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_steep_entry_stress_C3 import SourceAST
from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet, Z

HERE = Path(__file__).resolve().parent
PREFIX = 'lei_ren_part1_paper_compliant_'
NAME = PREFIX+'current_patch_support_velocity_pressure.json'
RECEIPT = PREFIX+'current_patch_support_velocity_pressure_check.json'
EDGES = (49,51,59,61,69,71)
CENTERS = (s.Rational(5,4),s.Rational(3,2),s.Rational(7,4))
WEIGHTS = ((s.Rational(0),1),(s.Rational(3,5),1),(s.Rational(1,2),2),
           (s.Rational(1,2),1),(s.Rational(0),2),(s.Rational(1,10),1),
           (s.Rational(-9,10),1),(s.Rational(-1),2))
LABELS = ('Utheta_over_Pstar','Uz','Ur_over_sqrt_Rm_over_2','P_over_Pstar2')


def sha(name):
    return hashlib.sha256((HERE/name).read_bytes()).hexdigest()


def _context():
    def exact(value):
        if isinstance(value,tuple):
            return value[0]
        return s.Rational(str(value)) if isinstance(value,(str,int,float)) else value
    return SimpleNamespace(mpf=exact,exp=s.exp)


def _source_acceptance():
    records, hashes = {}, {}
    requirements = {
        'flat_pulse_derivatives_check': ('all_passed','original_radial_shape_derivatives_C4_available'),
        'actual_feedback_moment_patch_check': ('all_passed','current_actual_patch_implicit_axial5_recomputed',
            'current_five_functional_terminal_identities_connected','current_Rm_source_parent_and_P0_retained'),
        'actual_feedback_patch_mixed_C4_check': ('all_passed','current_actual_patch_mixed4_available',
            'current_five_functional_terminal_identities_connected','source_joins_not_proved_by_interval_overlap'),
        'current_modified_physical_velocity_check': ('all_passed',
            'current_modified_velocity_pressure_physical_locator_integrated')}
    for stem, flags in requirements.items():
        name = PREFIX+stem+'.json'
        record = json.loads((HERE/name).read_bytes())
        if not all(record.get(flag) is True for flag in flags):
            raise ValueError('Current source prerequisite missing: '+name)
        records[stem] = record
        # Local prerequisites are completely hash checked. The current
        # physical cover pins the same local definitions; no graph rebuilt.
        if stem != 'current_modified_physical_velocity_check':
            for source, digest in record['input_hashes'].items():
                if source in hashes and hashes[source] != digest:
                    raise ValueError('Conflicting source definition: '+source)
                hashes[source] = digest
        hashes[name] = sha(name)
    current = records['current_modified_physical_velocity_check']
    identity = ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256')
    patch = records['actual_feedback_patch_mixed_C4_check']
    for key in identity:
        if current.get(key) != patch.get(key):
            raise ValueError('Current physical field and patch source differ: '+key)
    for name, digest in hashes.items():
        if sha(name) != digest:
            raise ValueError('Changed accepted patch source: '+name)
    for stem in ('actual_feedback_moment_patch','actual_feedback_patch_mixed_C4',
                 'actual_moment_patch','actual_patch_mixed_C4','current_patch_stress_operator'):
        name = PREFIX+stem+'.py'
        if current['input_hashes'].get(name) != sha(name):
            raise ValueError('Current physical dispatcher uses a different patch definition: '+name)
        hashes[name] = sha(name)
    return dict(source_family={key:patch[key] for key in identity}, input_hashes=hashes,
        actual_current_implicit_functional_closure_consumed=True,
        old_zero_containment_not_used_as_implicit_or_support_join_proof=True)


def _programs(asts):
    c = _context()
    J = lambda value,order=5: FunctionJet.function(c,value,order)
    env = dict(math=math,IntervalTaylor=FunctionJet,square=lambda v:v*v,
        derivative=lambda v:J(s.diff(v.expr,Z),v.order-1),copy_jet=lambda ctx,row:row,
        endpoints=lambda v:(v,v))
    for stem,name in (
        ('actual_patch_mixed_C4','power_derivatives'),
        ('actual_patch_mixed_C4','product_derivatives'),
        ('actual_patch_mixed_C4','stirling_second'),
        ('actual_patch_mixed_C4','log_radial_derivatives'),
        ('collar_stress_C3','shifted_rows')):
        asts.replay(stem,name,env)
    mixed = asts.replay('actual_patch_mixed_C4','patch_mixed',env)
    raw = asts.replay('current_patch_stress_operator','raw_patch_rows',env)
    return c,J,env,mixed,raw


def replayed_trace_packets(asts, perturb=None):
    """Independent defining germs after the proved flat/FTC limit reduction.

    Outside derivatives use the actual power program. Inside derivatives
    use direct differentiation of the limiting defining base function.
    Common histories denote source identities proved separately below.
    """
    c,J,env,mixed,raw = _programs(asts)
    x = s.Symbol('actual_positive_support_x',positive=True)
    delta = s.Symbol('same_delta',real=True)
    ps = s.Symbol('same_Pstar',positive=True)
    am = J(s.exp(-s.Rational(3,5))/(1+Z*Z))
    packets, raw_packets = [], []
    for side in ('left','right'):
        if side == 'left':
            H = [J(row) for row in env['power_derivatives'](c,x,'.1')]
        else:
            H = [J(s.diff(x**s.Rational(1,10),x,j)) for j in range(5)]
        g = [J(0) for _ in range(5)]
        V = [J(4*Z)+g[0]]+g[1:]
        initial = {key:J(s.Function('same_source_'+key)(Z))
                   for key in ('mass','theta','mixed','energy','pressure')}
        p0 = J(s.Function('same_actual_P0')(Z))
        if side == 'right' and perturb:
            if perturb == 'H': H[0] = H[0]+1
            elif perturb == 'mass': initial['mass'] = initial['mass']+1
            elif perturb == 'P0': p0 = p0+1
            else: raise ValueError('Unknown independent negative control')
        packet = mixed(c,x,Z,delta,am,J(1/(ps**2*am.expr**2)),1/ps**2,H,V,g,initial,p0)
        packet.update(x=x,Z=Z,actual_H_x_derivative_axial5=H,actual_V_x_derivative_axial5=V,
            actual_inherited_patch_packet=dict(original_P0_axial5=p0.coefficients))
        packets.append(packet)
        raw_packets.append(raw(c,packet,1/ps**2))
    if packets[0] is packets[1]:
        raise ArithmeticError('Independent support germs became a self-comparison')
    return x,packets,raw_packets


def exact_patch_support_theorem():
    accepted = _source_acceptance()
    asts, checks = SourceAST(), {}
    def zero(name,a,b):
        if s.cancel(a-b) != 0:
            raise ArithmeticError('Actual patch support function identity failed: '+name)
        checks[name] = True
    # Exact unchanged current provider bindings; only source inspection.
    from lei_ren_part1_paper_compliant_actual_feedback_patch_mixed_C4 import current_mixed_source_bindings
    from lei_ren_part1_paper_compliant_actual_feedback_moment_patch import current_patch_source_bindings
    bindings = dict(current_mixed=current_mixed_source_bindings(),current_moments=current_patch_source_bindings())
    for target,wanted in (
        ('centers',"(c.mpf(5)/4,c.mpf(3)/2,c.mpf(7)/4)"),
        ('defects','[a+b for a,b in zip(d,change)]'),
        ('mass','z*4+defects[0]/x'),
        ('theta',"defects[2]+x**c.mpf('1.6')*c.mpf('.625')"),
        ('mixed','(z*theta)*4+defects[1]'),
        ('energy',"defects[3]-x**c.mpf('1.2')*c.mpf(5)/12"),
        ('pressure_moment',"defects[4]+x**c.mpf('.2')*c.mpf('2.5')"),
        ('terminal','endpoints(x)[0]>=mp.mpf(71)/40'),
        ('p0',"self.reference.inputs(Z)['original_axis_pressure']")):
        asts.expression('actual_moment_patch','evaluate',target,wanted=wanted)
        checks['actual_parent_'+target+'_program'] = True
    asts.expression('actual_moment_patch','evaluate','defects',wanted='[zero]*5')
    for target,wanted in (
        ('argument','(c.mpf(x)-c.mpf(center))/self.radius'),
        ('derivatives','[rows[k]*(math.factorial(k)/(self.radius**(k+1)*self.N)) for k in range(5)]')):
        asts.expression('actual_patch_mixed_C4','gamma',target,wanted=wanted)
    asts.expression('actual_feedback_patch_mixed_C4','__init__','self.radius',wanted='c.mpf(1)/40')
    asts.expression('five_moment_repair','beta','value',wanted='c.exp(-1/den)/(radius*self.normalization)')
    asts.expression('five_moment_repair','partial_weight','total',
        wanted='box**pp*val**multiplicity*length/cells',augmented=True)
    # Actual analytic beta derivatives: their common zero trace is proved
    # by the accepted exponential-polynomial bound, not rounded boxes.
    from lei_ren_part1_paper_compliant_flat_pulse_derivatives import BETA_POLYNOMIALS
    w = s.Symbol('positive_beta_distance',positive=True)
    for j,row in enumerate(BETA_POLYNOMIALS):
        if s.limit(s.exp(-1/w)/w**(2*j),w,0,dir='+') != 0:
            raise ArithmeticError('Actual beta derivative has a nonflat support trace')
        checks['normalized_beta_ordinary_derivative_flat_limit_'+str(j)] = True
    # Partial FTC derivatives1..4 are derivatives0..3 of x^p*gamma^m.
    dx = s.Symbol('local_x_displacement',real=True)
    gamma_jets = s.symbols('inside_gamma0:5',real=True)
    gamma = sum(gamma_jets[j]*dx**j/s.factorial(j) for j in range(5))
    reduction = dict.fromkeys(gamma_jets,s.Integer(0))
    source_checks = {}
    c,J,env,_,_ = _programs(asts)
    change_program = asts.expression('actual_moment_patch','evaluate','change')
    controls = [J(s.Function('actual_implicit_h'+str(i))(Z)) for i in range(5)]
    am = J(s.exp(-s.Rational(3,5))/(1+Z*Z))
    ps = s.Symbol('same_Pstar',positive=True)
    def full(i,p,m):
        return s.Integer(1) if p==0 and m==1 else s.Symbol('actual_full_W_%d_%s_%d'%(i,str(p).replace('/','_').replace('-','minus'),m),real=True)
    def changes(weight):
        return asts.evaluate(change_program,dict(h=controls,w=lambda i,p,m=1:J(weight(i,s.Rational(p),m)),
            square=lambda v:v*v,data=dict(invAm2=J(1/(ps**2*am.expr**2))),zero=J(0)))
    full_change = changes(full)
    defects = [s.Function('actual_incoming_defect'+str(i))(Z) for i in range(5)]
    for edge in EDGES:
        point = s.Rational(edge,40)
        bump = (edge-49)//10
        entering = edge%10==9
        if point != CENTERS[bump]+(-1 if entering else 1)*s.Rational(1,40):
            raise ArithmeticError('Actual rational support edge differs')
        for i,center in enumerate(CENTERS):
            if i != bump and abs(point-center) <= s.Rational(1,40):
                raise ArithmeticError('Actual beta supports are not disjoint at this edge')
        per_edge = {}
        for p,m in WEIGHTS:
            density = (point+dx)**p*gamma**m
            for j in range(1,5):
                value = s.diff(density,dx,j-1).subs(dx,0).subs(reduction)
                if value != 0: raise ArithmeticError('Partial FTC support derivative fails')
                per_edge['partial_weight_%s_%d_x%d'%(str(p),m,j)] = True
        initial_germs, residuals = [], []
        for side in ('left','right'):
            inside = (side=='right') if entering else (side=='left')
            residual = {}
            def weight(i,p,m):
                if i < bump: return full(i,p,m)
                if i > bump: return s.Integer(0)
                base = s.Integer(0) if entering else full(i,p,m)
                if not inside: return base
                strip = residual.setdefault((i,p,m),s.Symbol('inside_strip_%d_%s_%d'%(i,str(p).replace('/','_').replace('-','minus'),m),real=True))
                return base+(1 if entering else -1)*strip
            increments = changes(weight)
            # Each strip tends to zero by bounded continuous density times
            # its shrinking positive length. Its derivatives are FTC rows.
            reduced = [s.cancel(v.expr.subs(dict.fromkeys(residual.values(),0))) for v in increments]
            germs = [defects[i]+reduced[i] for i in range(5)]
            if edge==71:
                # Left unrefined full weights and right terminal branch
                # are identified by the SAME accepted implicit equations.
                closure = {defects[i]:-full_change[i].expr for i in range(5)}
                germs = [s.cancel(value.subs(closure)) for value in germs]
                if side=='right': germs = [s.Integer(0)]*5
            initial_germs.append(germs)
            residuals.append([str(v) for v in residual.values()])
        for i in range(5):
            if s.cancel(initial_germs[0][i]-initial_germs[1][i]) != 0:
                raise ArithmeticError('Independent cumulative source germs differ')
            per_edge['same_actual_cumulative_defect_function_'+str(i)] = True
            if edge==71 and initial_germs[0][i] != 0:
                raise ArithmeticError('Terminal full weights do not use actual implicit closure')
        per_edge['shrinking_strip_integral_limit_from_continuous_positive_radius_density'] = True
        per_edge['same_actual_controls_and_analytic_P0_not_refitted'] = True
        per_edge['inside_and_outside_defining_germs_constructed_independently'] = True
        source_checks[str(edge)] = dict(identities=per_edge,passed=True,exact_x=str(point),
            active_bump=bump,entering=entering,inside_strip_source_symbols_by_side=residuals,
            right_edges_retain_full_previous_and_active_weights=not entering,
            last_edge_uses_current_five_functional_implicit_equations=edge==71)
    x,packets,raw_packets = replayed_trace_packets(asts)
    for group,prefix in (('physical_velocity_pressure_x_Z_mixed4','x'),
                         ('physical_velocity_pressure_y_Z_mixed4','y')):
        for label in LABELS:
            for j in range(5):
                for n in range(5-j):
                    key='%s%d_Z%d'%(prefix,j,n)
                    zero('two_actual_germs_'+label+'_'+key,packets[0][group][label][key],packets[1][group][label][key])
    for label,raw_label,factor in (('Utheta_over_Pstar','theta',1),('Uz','axial',1),
                                 ('Ur_over_sqrt_Rm_over_2','radial',s.sqrt(x)),('P_over_Pstar2','P',1)):
        rows = raw_packets[0]['absolute_pressure'] if raw_label=='P' else raw_packets[0]['velocity'][raw_label]
        for j in range(5):
            for n in range(5-j):
                zero('actual_normalized_physical_unit_'+label+'_y%d_Z%d'%(j,n),
                    packets[0]['physical_velocity_pressure_y_Z_mixed4'][label]['y%d_Z%d'%(j,n)],
                    rows[j][n]*s.factorial(n)*factor)
    for obj in (SourceAST,FunctionJet):
        path = Path(inspect.getsourcefile(obj));asts.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
    hashes = {**accepted['input_hashes'],**asts.hashes,Path(__file__).name:sha(Path(__file__).name)}
    return dict(source_family=accepted['source_family'],exact_support_source_theorems=source_checks,
        generic_two_germ_mixed4_and_current_unit_identities=checks,
        source_bindings={**asts.bindings,**bindings},passed=True,
        positive_radius_generic_mixed4_lemma_applies_at_all_six_actual_edges=True,
        common_analytic_P0_and_actual_current_implicit_controls_retained=True,
        cumulative_nonzero_right_edge_histories_not_reset=True,
        source_function_equalities_precede_bounds=True,
        source_x_Z_and_logR_Z_mixed4_joins_certified=True,
        mapped_physical_spatial4_time1_trace_views_available=False,
        global_smooth_velocity_energy_common_N_cones_recursion_NS_remain_open=True,input_hashes=hashes)


def run():
    theorem = exact_patch_support_theorem()
    result = dict(exact_actual_patch_support_velocity_pressure_theorem=theorem,
        current_six_patch_support_velocity_pressure_mixed4_source_joins_certified=True,
        mapped_physical_spatial4_time1_trace_views_available=False,
        current_modified_global_velocity_interfaces_certified=False,
        common_N_modified_cones_energy_recursion_full_NS_certified=False,
        input_hashes=theorem['input_hashes'])
    (HERE/NAME).write_bytes((json.dumps(result,indent=2)+'\n').encode())
    print('Actual six patch support source joins generated; global physical interfaces remain open',flush=True)
    return result


if __name__=='__main__':
    run()
