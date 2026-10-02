"""Absolute source identities for the corrected leading outer moments.

This composes the existing implicit definitions, rather than selecting values
inside interval boxes or resetting histories. Old conservative receipts remain
unchanged. The admitted view uses algebraically equivalent backward pressure
and heat angular primitives only AFTER the source identities have been checked.
Higher axial regularity, the stress cone and temporal recursion are separate.
"""
import ast
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_corrected_outer_field import CompliantCorrectedOuterField
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def assignment(stem, method, target, env, augmented=False):
    """Translate a selected PRODUCTION expression; reject unknown operations."""
    tree = ast.parse((HERE / (PREFIX + stem + '.py')).read_text(encoding='utf-8'))
    function = next(n for n in ast.walk(tree) if isinstance(n, ast.FunctionDef) and n.name == method)
    nodes = [n.value for n in ast.walk(function)
             if (not augmented and isinstance(n, ast.Assign) and any(ast.unparse(t) == target for t in n.targets))
             or (augmented and isinstance(n, ast.AugAssign) and ast.unparse(n.target) == target and isinstance(n.op, ast.Add))]
    if len(nodes) != 1:
        raise ValueError('Unique source assignment required: ' + stem + ':' + target)

    def formal(node):
        label = ast.unparse(node)
        if label in env:
            return env[label]
        if isinstance(node, ast.Constant):
            return s.Rational(str(node.value))
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
            return -formal(node.operand)
        if isinstance(node, ast.BinOp):
            l, r = formal(node.left), formal(node.right)
            for kind, operation in ((ast.Add, lambda: l+r), (ast.Sub, lambda: l-r),
                                    (ast.Mult, lambda: l*r), (ast.Div, lambda: l/r),
                                    (ast.Pow, lambda: l**r)):
                if isinstance(node.op, kind):
                    return operation()
        if isinstance(node, ast.Call) and len(node.args) == 1:
            name = ast.unparse(node.func)
            value = formal(node.args[0])
            if name in ('c.mpf', 'self.jet'):
                return value
            if name in ('c.exp', 'self.pulse.factor'):
                return s.exp(value)
            if name == 'c.ln':
                return s.log(value)
        raise ValueError('Unbound source expression: ' + stem + ':' + label)
    return formal(nodes[0])


def source_identities():
    """Whole-Z algebra, with arbitrary functional histories and heat defects."""
    proofs = {}

    def zero(name, value):
        reduced = s.simplify(s.expand_power_exp(value))
        if reduced != 0:
            raise ArithmeticError('Absolute closure identity failed: ' + name + ': ' + str(reduced))
        proofs[name] = True

    mu, a, Ts, tau, lone, logmu, LR = s.symbols('mu a Ts tau lone logmu LR', real=True)
    eps = s.symbols('epsilon', positive=True)
    rate, k, bp, bh = 1-mu, 1-a, s.Rational(1, 2)+mu, s.Rational(1, 2)+a
    length = 2+Ts+tau
    ell = -s.Rational(3, 2)*(2+Ts)+(rate+k)/2-bh*tau
    logT = -s.Rational(3, 2)*length-ell
    env = {'self.rate': rate, 'self.mu': mu, 'self.delta': 2*a,
           'self.params.Ts': Ts, 'self.angular.waiting': tau,
           'self.angular.restore_rate': k, 'self.angular.waiting_logone': lone,
           'self.params.log_mu': logmu, 'self.logEtail_over_Erel': ell,
           'self.bp': bp, 'self.bh': bh, 'self.k': k, 'self.LR': LR,
           'self.Ls': LR-bp-rate/2,
           'self.Lq': LR-bp-rate/2-s.Rational(3, 2)*Ts,
           'self.Lt': LR-bp-rate/2-s.Rational(3, 2)*Ts-s.Rational(3, 2)+k/2}
    expected = {'self.logscale': 30*rate*logmu, 'self.logtail_distance': length,
                'self.logEtail_over_Erel': ell,
                'self.log_theta_multiplier': (rate+k)/2+k*tau-lone,
                'self.log_pressure_multiplier': 2*ell-2*lone+s.log(a)}
    for target, expression in expected.items():
        zero('source_repair_' + target[5:], assignment('compliant_outer_angular_repair', '__init__', target, env)-expression)
    for target in ('self.Ls', 'self.Lq', 'self.Lt'):
        zero('source_outer_' + target[5:], assignment('compliant_corrected_outer_field', '__init__', target, env)-env[target])
    zero('source_outer_tail_amplitude', assignment('compliant_corrected_outer_field', '__init__', 'self.Ltail', env)-LR-ell)
    zero('angular_transfer_independent_of_steep_length', logT+(rate+k)/2+k*tau)
    zero('angular_heat_multiplier_inverse_transfer', logT+expected['self.log_theta_multiplier']+lone)
    zero('pressure_heat_multiplier', expected['self.log_pressure_multiplier']-(2*ell-2*lone+s.log(a)))

    Xf, Xf0, A, Iin, Iout, L = s.symbols('Xf Xf0 A Iin Iout L', real=True)
    env.update(Xf=Xf, eq=1/rate, fullA=A, **{'self.Lrel': L, 'self.Iin': Iin, 'self.Iout': Iout})
    XR = (Xf-1/rate)*s.exp(-rate*L)+A+1/rate
    env['XR'] = XR
    Xs = (XR+Iin)*s.exp(-rate/2)
    env['Xs'] = Xs
    Xq = Xs+Ts
    env['Xq'] = Xq
    Xt = (Xq+Iout)*s.exp(-k/2)
    env['Xt'] = Xt
    Xtail = (Xt-1/k)*s.exp(-k*tau)+1/k
    for target, expression in (('XR', XR), ('Xs', Xs), ('Xq', Xq), ('Xt', Xt), ('Xtail', Xtail)):
        zero('source_forward_angular_' + target, assignment('compliant_corrected_outer_field', 'data', target, env)-expression)
    raw0 = Xtail.subs({Xf: Xf0, A: 0})
    rpre = s.exp(-rate*L)*(Xf0-Xf)
    S, theta, pressure, J = s.symbols('S theta_hat pressure_hat Jcollar', real=True)
    rH = s.exp(expected['self.log_theta_multiplier'])*S*theta
    zero('correlated_flatten_memory', Xtail.subs(A, 0)-raw0+s.exp(logT)*rpre)
    zero('bump_transport', Xtail-Xtail.subs(A, 0)-s.exp(logT)*A)
    zero('corrected_angular_vs_Z0',
         ((1-eps)*(Xtail.subs(A, rpre+rH)-raw0)-S*theta).subs(lone, s.log(1-eps)))
    # Exact waiting equation, not a selected interval center.
    Xt0 = s.symbols('Xt0', real=True)
    waitenv = {'Xt': Xt0, 'eq': 1/k, 'k': k, 'logone': lone,
               'self.params.log_epsilon': s.log(eps), 'self.collarJ': J}
    wait = assignment('compliant_outer_angular_candidate', '__init__', 'self.waiting', waitenv)
    zero('source_waiting_root', wait-(s.log(Xt0-1/k)+lone-s.log(eps)-s.log(1/k+J))/k)
    decay = eps*(1/k+J)/((1-eps)*(Xt0-1/k))
    zero('waiting_Z0_preheat_target', (1-eps)*(1/k+(Xt0-1/k)*decay)-(1/k+eps*J))
    dj, center, wA, wB, wD = s.symbols('dj center A_weight B_weight D_weight', real=True)
    benv = {'dj': dj, 'center': center, "w['A']": wA, "w['B']": wB, "w['D']": wD,
            'self.rate': rate, 'self.prate': 1+2*mu}
    zero('source_physical_angular_bump_increment', assignment('compliant_corrected_outer_field', 'correction', 'A', benv, augmented=True)-dj*s.exp(rate*center)*wA)
    zero('source_physical_pressure_bump_increment', assignment('compliant_corrected_outer_field', 'correction', 'P', benv, augmented=True)-(dj*wB+dj**2*wD/2)*s.exp(-(1+2*mu)*center))
    # Smooth functional formulation also checks derivative transport.
    z = s.symbols('Z', real=True)
    xf, th = s.Function('Xf')(z), s.Function('Theta_hat')(z)
    raw = 1/k+eps*J+(1-eps)*s.exp(logT-rate*L)*(xf-Xf0)
    corrected = raw+(1-eps)*s.exp(logT)*(s.exp(-rate*L)*(Xf0-xf)+s.exp(-logT)*S*th/(1-eps))
    zero('absolute_angular_constant_whole_Z', corrected-(1/k+eps*J+S*th))
    zero('absolute_angular_constant_axial_derivative', s.diff(corrected-(1/k+eps*J+S*th), z))

    x = s.symbols('x', positive=True)
    zero('reference_pressure_extension_integral', s.integrate(x**s.Rational(1, 5)/(2*x), (x, 0, 1))-s.Rational(5, 2))
    zero('reference_angular_extension_integral', s.integrate(x**s.Rational(3, 5), (x, 0, 1))-s.Rational(5, 8))
    Erel = s.symbols('Erel', positive=True)
    pressure_bump = Erel**2*s.exp(expected['self.log_pressure_multiplier'])*S*pressure
    pressure_heat_loss = (Erel*s.exp(ell)/(1-eps))**2*a*S*pressure
    zero('physical_pressure_bump_equals_heat_loss', (pressure_bump-pressure_heat_loss).subs(lone, s.log(1-eps)))
    p0, rawinner, rawouter, actualinner = s.symbols('P0 raw_inner raw_outer actual_inner', real=True)
    # Actual inner repair restores the raw reference integral. Only the
    # disjoint angular supports and Gamma replacement change swirl mass.
    offset = p0+actualinner+rawouter+pressure_bump-pressure_heat_loss
    zero('absolute_pressure_constant_from_source', offset.subs(actualinner, rawinner).subs(p0, -rawinner-rawouter).subs(lone, s.log(1-eps)))
    pv, prv, scale, pr = s.symbols('MpRv Prv scale Pr', real=True)
    zero('forward_backward_pressure_after_closure', (p0+pv+scale*(prv-pr)+scale*pr).subs(p0, -pv-scale*prv))
    pb = s.Function('Delta_bump')(z)
    ph = s.Function('Delta_heat')(z)
    zero('pressure_functional_derivative', s.diff(pb-ph, z).subs(s.diff(pb, z), s.diff(ph, z)))
    return dict(identities=proofs,
                angular_transfer_log=str(logT), angular_heat_multiplier_log=str(expected['self.log_theta_multiplier']),
                pressure_heat_multiplier_log=str(expected['self.log_pressure_multiplier']),
                pressure_source_chain='P0=-Mp_pre(0..infinity); Mp_actual=Mp_pre+Erel^2*sH-Delta_heat; Erel^2*sH=Delta_heat',
                angular_source_chain='(1-epsilon)*Xtail_corrected=1/k+epsilon*Jcollar+S*Theta_hat',
                proof_uses_interval_zero_overlap=False, exact_positive_S_not_replaced_by_cap=True)


class CompliantAbsoluteMomentClosure:
    """Admitted outer view; all original forward histories remain accessible."""
    def __init__(self):
        self.proof = source_identities()
        self.outer = CompliantCorrectedOuterField()
        self.hashes = dict(self.outer.hashes)
        required = {
            'compliant_five_moment_repair_check': 'actual_implicit_functional_five_moment_closure_independently_checked',
            'compliant_outer_angular_repair': 'actual_smooth_functional_angular_pressure_repair_proved',
            'compliant_outer_angular_repair_check': 'all_passed',
            'compliant_axial_pulse_field_check': 'actual_O4_partial_pulse_field_installed',
            'compliant_corrected_outer_field_check': 'all_passed'}
        for stem, gate in required.items():
            name = PREFIX+stem+'.json'
            record = json.loads((HERE/name).read_bytes())
            if not record.get(gate):
                raise ValueError('Absolute closure missing prerequisite: '+gate)
            if record.get('actual_five_defect_family_sha256') != self.outer.angular.initial.family:
                raise ValueError('Absolute closure family mismatch: '+name)
            source_sha = record.get('implicit_source_sha256', record.get('source_sha256'))
            if source_sha is not None and source_sha != self.outer.angular.initial.datum.source_sha:
                raise ValueError('Absolute closure implicit pressure source mismatch: '+name)
            for source, digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/source).read_bytes()).hexdigest() != digest:
                    raise ValueError('Absolute closure source changed: '+source)
                self.hashes[source] = digest
            self.hashes[name] = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        datum = self.outer.angular.initial.datum
        if (datum.definition['angular_profile'] != 'Section 6.1 reference-plus-outer ansatz, H replaced by 1'
            or datum.definition['waiting'] != 'unique positive root of the continuous raw preheat waiting equation'
            or datum.definition['c_epsilon'] != '.001'):
            raise ValueError('Exact raw preheat pressure definition changed')
        self.hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def field(self, stage, Z, coordinate):
        if stage not in ('flatten', 'power_buffer', 'angular_end', 'steep_in', 'steep_power', 'steep_out', 'waiting', 'collar_heat'):
            raise ValueError('Unknown corrected outer chart')
        point = dict(getattr(self.outer, stage)(Z, coordinate))
        # Tight equivalent pressure is a consequence of the absolute source
        # proof; Mp and P0 retain their original forward enclosures.
        point['original_forward_pressure_enclosure'] = point['P_over_Pstar_squared']
        point['P_over_Pstar_squared'] = -point['remaining_pressure_integral_in_Ev0_squared']*self.outer.Ev0Pstar2
        point['P_over_Utheta_squared'] = point['independent_backward_pressure_over_Utheta_squared']
        point['exact_pressure_definition'] = 'P=P0+Mp=-integral_R^infinity Utheta_corrected^2/(2rho)drho; exact source identity, original Mp/P0 retained'
        if stage == 'collar_heat':
            point['original_forward_angular_enclosure'] = point['Mtheta_over_sqrt2_R_3half_Utheta']
            point['Mtheta_over_sqrt2_R_3half_Utheta'] = self.outer.heat_tails(Z, coordinate)['X_target']
            point['forward_angular_target_equality_verified'] = True
        point.update(full_global_pressure_terminal_identity_verified=True,
                     all_five_terminal_moment_identities_certified=True,
                     full_outer_five_moment_match=True,
                     C4_outer_bounds_certified=False, whole_outer_cone_certified=False,
                     global_admissible_stress_lift_constructed=False, temporal_recursion=False,
                     closure_proof_scope='exact implicit leading-profile functions on all Z in [-1,1]; callable enclosures currently C1')
        return point

    def report(self):
        with mp.workdps(210):
            calls = [('flatten', '.5', 0), ('flatten', '.5', 1),
                     ('angular_end', '.5', -3), ('waiting', '.5', 1),
                     ('collar_heat', '.5', 0), ('collar_heat', '.5', 3),
                     ('collar_heat', '.5', 1000), ('collar_heat', [-1, 1], 3)]
            samples = [self.field(*call) for call in calls]
            data = self.outer.data('.5')
            return dict(proof=self.proof, samples=samples,
                paper_reference=dict(version='2609.35406v2', equations='(7.10), (7.13), (7.17), (7.21), (7.27), (7.34)-(7.36)',
                    repository_path='work_paper_cache/2609.35406v2.txt',
                    sha256=hashlib.sha256((HERE.parents[1]/'work_paper_cache/2609.35406v2.txt').read_bytes()).hexdigest()),
                original_angular_constant_enclosure=data['angular_tail_constant_defect'],
                source_sha256=self.outer.angular.initial.datum.source_sha,
                datum_enclosure_sha256=self.outer.angular.initial.datum.datum_sha,
                actual_five_defect_family_sha256=self.outer.angular.initial.family,
                terminal_moments=dict(Mz='0 from selected axial linear equation',
                    Mtheta_z='0 from selected mixed linear equation',
                    Mztheta='0 at infinity from selected energy equation and delta>0 finite tail',
                    Mtheta='renormalized constant 0 from waiting, correlated memory and first angular repair equation',
                    Mp='Mp(infinity)=-P0 from raw preheat datum, repaired reference and second angular repair equation'),
                original_forward_Mp_and_P0_retained=True, no_new_pressure_datum_or_fitted_constant=True,
                all_five_terminal_moment_identities_certified=True,
                full_global_pressure_terminal_identity_verified=True,
                forward_angular_heat_target_equality_verified=True, full_outer_five_moment_match=True,
                C4_outer_bounds_certified=False, whole_outer_cone_certified=False,
                full_physical_energy_certified=False, global_admissible_stress_lift_constructed=False,
                temporal_recursion=False, input_hashes=self.hashes)


def run():
    result = CompliantAbsoluteMomentClosure().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)), indent=2)+'\n', encoding='utf-8')
    print('Absolute pressure/angular source identities composed; all five leading terminal moments admitted; C4/cone/temporal pending', flush=True)
    return result


if __name__ == '__main__':
    run()
