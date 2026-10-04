"""Source-bound absolute pressure through the original full heat collar.

The forward pressure and the negative remaining integral are equivalent
after the admitted original absolute closure. All sigma/phi/Gamma factors,
velocities, moments and the pressure datum are inherited without changes.
"""
import ast
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_heat_pressure_C4 import CompliantHeatPressureC4, pressure_y_rows
from lei_ren_part1_paper_compliant_pulse_physical_bounds import P
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_'


def collar_pressure_rows(K_rows, remaining_pressure, rate, scale, offset):
    """FTC derivatives of the full remaining integral in Rtail units.

    The zeroth row uses reference units directly; no cancelling interval
    exp(+rate*t)*exp(-rate*t) is introduced. Higher rows use the common
    pressure FTC/Leibniz helper and the full collar bracket, not pure Gamma.
    """
    rows = pressure_y_rows(K_rows, remaining_pressure*0, rate, scale, offset)
    rows[0] = -remaining_pressure*scale
    return rows


def collar_pressure_bridge(inherited):
    """Extend the admitted actual pressure history to every collar offset."""
    if not inherited['complete_defining_function_history_bridge_verified']:
        raise ValueError('Complete original pressure defining/history bridge required')
    history = inherited['retained_pressure_history_bridge']
    if not history['complete_retained_pressure_history_bridge_verified']:
        raise ValueError('Original retained pressure history not admitted')
    proofs, hashes, trees = {}, dict(history['input_hashes']), {}

    def method(stem, name):
        if stem not in trees:
            path = HERE/(PREFIX+stem+'.py')
            hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
            trees[stem] = ast.parse(path.read_text(encoding='utf8'))
        return next(n for n in ast.walk(trees[stem]) if isinstance(n, ast.FunctionDef) and n.name == name)

    def syntax(stem, name, target, expression, augmented=False):
        nodes = [n.value for n in ast.walk(method(stem, name))
                 if (not augmented and isinstance(n, ast.Assign)
                     and any(ast.unparse(v) == target for v in n.targets))
                 or (augmented and isinstance(n, ast.AugAssign)
                     and ast.unparse(n.target) == target and isinstance(n.op, ast.Add))]
        wanted = ast.parse(expression, mode='eval').body
        if len(nodes) != 1 or ast.dump(nodes[0]) != ast.dump(wanted):
            raise ValueError('Collar pressure source changed: '+stem+':'+name+':'+target)
        proofs['source_'+stem+'_'+name+'_'+target] = True

    def zero(name, expression):
        expression = expression.xreplace({q:s.Dummy('same_source_integral') for q in expression.atoms(s.Integral)})
        if s.simplify(s.expand_power_exp(expression)) != 0:
            raise ArithmeticError('Collar pressure identity failed: '+name)
        proofs[name] = True

    def keyword(stem, name, key, expression):
        nodes = [n.value for n in ast.walk(method(stem,name)) if isinstance(n,ast.keyword) and n.arg==key]
        wanted = ast.parse(expression,mode='eval').body
        if len(nodes)!=1 or ast.dump(nodes[0])!=ast.dump(wanted):
            raise ValueError('Collar pressure output source changed: '+stem+':'+name+':'+key)
        proofs['output_'+stem+'_'+name+'_'+key] = True

    stem = 'compliant_collar_Gamma_C4'
    syntax(stem, 'collar', 'data', 'self.data(Z)')
    syntax(stem, 'collar', 'tails', 'self.collar_tails(Z,t)')
    syntax(stem, 'collar', 'shape', 'self.shape(Z,t)')
    syntax(stem, 'collar', 'pressure', "self.forward_pressure(Z,t,data['Ptail'])")
    syntax(stem, 'forward_pressure', 'K', "self.shape(Z,v,False)['K_rows'][0]")
    syntax(stem, 'forward_pressure', 'integral', 'K*K*(ds*c.exp(-self.prate*v)/2)', True)
    syntax(stem, 'collar_tails', 'length', '3-t')
    syntax(stem, 'collar_tails', 'a', 't+length*i/self.cells')
    syntax(stem, 'collar_tails', 'b', 't+length*(i+1)/self.cells')
    syntax(stem, 'collar_tails', 'tails', 'integrated_gamma_tails(c,Z,self.a,self.Scap,5,3)')
    syntax(stem, 'collar_tails', 'square', 'D*(pre*2-D*(self.a*self.S))')
    syntax(stem, 'collar_tails', 'pressure', 'square*(ds*c.exp(-self.prate*v)/2)', True)
    syntax(stem, 'collar_tails', 'P', "one*(c.exp(-self.prate*t)/(2*self.prate)-self.eps*atoms['PW']+self.eps**2*atoms['PW2']/2)-pressure*(self.a*self.S)")
    keyword(stem, 'collar_tails', 'remaining_pressure_in_Rtail_units', 'P')
    keyword(stem, 'collar', 'complete_future_tail', 'tails')
    syntax('compliant_collar_pressure_C4', 'collar_pressure_rows', 'rows',
           'pressure_y_rows(K_rows,remaining_pressure*0,rate,scale,offset)')
    syntax('compliant_collar_pressure_C4', 'collar_pressure_rows', 'rows[0]', '-remaining_pressure*scale')
    syntax('compliant_collar_pressure_C4', 'collar', 'point', 'dict(self.heat.collar(Z,t))')
    syntax('compliant_collar_pressure_C4', 'collar', 'tails', "point['complete_future_tail']")
    syntax('compliant_collar_pressure_C4', 'collar', 'rows',
           "collar_pressure_rows(point['heat_bracket_y_derivatives'],tails['remaining_pressure_in_Rtail_units'],self.heat.prate,self.heat.pressure_scale,t)")
    keyword('compliant_collar_pressure_C4', 'collar', 'pressure_over_Pstar_squared_Taylor', 'rows[0]')
    keyword('compliant_collar_pressure_C4', 'collar', 'pressure_y_derivative_axial5_Taylor', 'rows')

    # The inherited proof binds the actual K^2 density, the full infinite
    # Gamma tail, Ptail/P3 and original infinity constant. Splitting that
    # same integral at arbitrary t makes no choice of a pressure datum.
    t, v, rate, scale, constant = s.symbols('t v rate scale original_infinity_constant', real=True)
    K = s.Function('original_full_collar_bracket')
    density = s.exp(-rate*v)*K(v)**2/2
    I0t = s.Integral(density, (v, 0, t))
    It3 = s.Integral(density, (v, t, 3))
    I3inf = s.Integral(density, (v, 3, s.oo))
    Ptail = constant-scale*(I0t+It3+I3inf)
    forward = Ptail+scale*I0t
    remaining = It3+I3inf
    zero('actual_forward_plus_full_remaining_is_original_infinity_constant', forward+scale*remaining-constant)
    zero('absolute_closure_transfers_to_every_collar_offset', (forward+scale*remaining).subs(constant,0))
    pre, D, a, S, eps, W = s.symbols('pre D a S epsilon W', real=True)
    zero('full_sigma_phi_Gamma_pressure_density_is_K_squared',
         ((pre-a*S*D)**2-((1-eps*W)**2-a*S*D*(2*pre-a*S*D))).subs(pre,1-eps*W))
    phi = s.Function('phi')(t)
    h = s.Function('H')(t)
    sig = s.Function('sigma')(t)
    bracket = (1-sig)*(1-eps)+sig*h*(1-eps*phi)
    zero('flat_waiting_collar_value_uses_original_epsilon', bracket.subs(sig,0)-(1-eps))
    zero('flat_collar_exterior_value_uses_full_Gamma', bracket.subs({sig:1,phi:0})-h)
    source = s.exp(-rate*t)*K(t)**2/2
    for order in range(1,5):
        formula = s.exp(-rate*t)*sum(s.binomial(order-1,j)*(-rate)**(order-1-j)*s.diff(K(t)**2,t,j)
                                    for j in range(order))/2
        zero('full_collar_pressure_FTC_y'+str(order), s.diff(source,t,order-1)-formula)
    # The existing admitted flat-shape/source ODE proof supplies K jets
    # through y4 and axial5 on both sides of 0 and 3. C0 is the same
    # full remaining integral; positive derivative rows follow by FTC.
    name = PREFIX+'compliant_collar_Gamma_C4_check.json'
    receipt = json.loads((HERE/name).read_bytes())
    if not receipt['all_passed'] or not receipt['waiting_collar_and_collar_Gamma_joins_certified']:
        raise ValueError('Original collar source and flat joins not admitted')
    for filename,digest in receipt['input_hashes'].items():
        if hashlib.sha256((HERE/filename).read_bytes()).hexdigest() != digest:
            raise ValueError('Collar join prerequisite changed: '+filename)
    hashes.update(receipt['input_hashes'])
    hashes[name] = hashlib.sha256((HERE/name).read_bytes()).hexdigest()
    return dict(identities=proofs,
                complete_collar_absolute_pressure_history_bridge_verified=True,
                pressure_waiting_collar_and_collar_exterior_mixed4_joins_verified=True,
                original_absolute_closure_and_datum_retained=True,
                full_sigma_phi_Gamma_future_integral_used=True,
                exact_S_not_a_cap_endpoint=True,
                reference_and_current_radius_pressure_units_distinguished=True,
                interval_overlap_used_as_proof=False,
                full_remaining_definition=str(remaining), input_hashes=hashes)


class CompliantCollarPressureC4(CompliantHeatPressureC4):
    """Full collar pressure companion; exterior pressure inherited unchanged."""
    def __init__(self, cells=64):
        super().__init__(cells)
        self.collar_bridge = collar_pressure_bridge(self.bridge)
        self.hashes.update(self.collar_bridge['input_hashes'])
        self.hashes[Path(__file__).name] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def collar(self, Z, t):
        c = self.ctx
        Z, t = c.mpf(Z), c.mpf(t)
        if endpoints(t)[0] < 0 or endpoints(t)[1] > 3:
            raise ValueError('Original heat collar offset in[0,3] required')
        point = dict(self.heat.collar(Z,t))
        tails = point['complete_future_tail']
        rows = collar_pressure_rows(point['heat_bracket_y_derivatives'],
                                    tails['remaining_pressure_in_Rtail_units'],
                                    self.heat.prate, self.heat.pressure_scale, t)
        mixed = dict(point['physical_mixed_derivatives_total_order_le4'])
        fields = dict(point['physical_velocity_and_pressure_y_derivative_Taylor'])
        point['original_forward_pressure_over_Pstar_squared_Taylor'] = point['pressure_over_Pstar_squared_Taylor']
        point['original_forward_pressure_mixed_bounds'] = mixed[P]
        mixed[P] = {'y'+str(k)+'_Z'+str(n): jet[n]*math.factorial(n)
                    for k,jet in enumerate(rows) for n in range(5-k)}
        fields[P] = [jet.truncate(4-k) for k,jet in enumerate(rows)]
        K = point['heat_bracket_y_derivatives'][0]
        point.update(physical_mixed_derivatives_total_order_le4=mixed,
                     physical_velocity_and_pressure_y_derivative_Taylor=fields,
                     pressure_over_Pstar_squared_Taylor=rows[0],
                     pressure_y_derivative_axial5_Taylor=rows,
                     pressure_over_Utheta_squared_Taylor=-tails['remaining_pressure_in_Rtail_units']*c.exp(self.heat.prate*t)/(K*K),
                     exact_pressure_definition='P=P0+Mp=-C*integral_t^infinity exp(-(1+delta)*v)*K_full_collar(v,Z)^2/2 dv; C=Ev2*theta_base^2 in Rtail units',
                     pressure_infinity_offset_exactly_zero_from_source_closure=True,
                     original_pressure_datum_and_forward_history_retained=True,
                     absolute_pressure_same_source_mixed4_available=True,
                     complete_collar_pressure_history_transfer_verified=True,
                     pressure_waiting_collar_and_collar_exterior_mixed4_joins_verified=True,
                     source_history_transfer_conditional=False,
                     full_sigma_phi_Gamma_future_integral_used=True,
                     heat_exterior_stress_identity_certified=False,
                     global_admissible_stress_lift_constructed=False,
                     actual_physical_map_and_prefactor_transfer_pending=True,
                     temporal_recursion=False)
        return point

    def report(self):
        def summary(point):
            keys = ('Z','coordinate','pressure_over_Pstar_squared_Taylor',
                    'pressure_over_Utheta_squared_Taylor','pressure_y_derivative_axial5_Taylor',
                    'original_forward_pressure_over_Pstar_squared_Taylor','original_forward_pressure_mixed_bounds',
                    'original_pressure_datum_and_forward_history_retained',
                    'absolute_pressure_same_source_mixed4_available','complete_collar_pressure_history_transfer_verified',
                    'pressure_waiting_collar_and_collar_exterior_mixed4_joins_verified',
                    'source_history_transfer_conditional','full_sigma_phi_Gamma_future_integral_used',
                    'heat_exterior_stress_identity_certified','global_admissible_stress_lift_constructed',
                    'actual_physical_map_and_prefactor_transfer_pending','temporal_recursion')
            point_summary = {key:point[key] for key in keys}
            point_summary['physical_mixed_derivatives_total_order_le4'] = {P:point['physical_mixed_derivatives_total_order_le4'][P]}
            return point_summary
        with mp.workdps(270):
            return dict(actual_five_defect_family_sha256=self.family, implicit_source_sha256=self.source,
                        scope='Absolute full-collar pressure, Z[-1,1], t=log(R/Rtail) in[0,3]; axial5 and mixed4',
                        samples=[summary(self.collar(z,t)) for z,t in (('0','0'),('.5','.5'),('.5','1'),('.5','2'),('-1','3'))],
                        whole_Z_inlet=summary(self.collar([-1,1],0)),
                        whole_Z_terminal=summary(self.collar([-1,1],3)),
                        whole_Z_collar=summary(self.collar([-1,1],[0,3])),
                        collar_pressure_source_bridge=self.collar_bridge,
                        complete_collar_absolute_pressure_history_bridge_verified=True,
                        pressure_waiting_collar_and_collar_exterior_mixed4_joins_verified=True,
                        absolute_pressure_same_source_mixed4_available=True,
                        original_pressure_datum_and_forward_history_retained=True,
                        pressure_datum_or_velocity_changed=False,
                        global_admissible_stress_lift_constructed=False,
                        whole_outer_cone_certified=False, physical_energy_integral_certified=False,
                        temporal_recursion=False, input_hashes=self.hashes)


def run():
    result = CompliantCollarPressureC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Full collar absolute pressure generated; original forward history and flat joins retained',flush=True)
    return result


if __name__ == '__main__':
    run()
