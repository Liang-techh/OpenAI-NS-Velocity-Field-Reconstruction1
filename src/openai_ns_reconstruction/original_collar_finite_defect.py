"""Source-level finite-pressure defect for the fourth m15 contribution."""
from __future__ import annotations

import ast

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext
import sympy as sp

from .original_collar_whole_z import hull, read, require, sign


CORE = 'src/openai_ns_reconstruction/original_core_numeric.py'
INNER = 'src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py'
JETS = 'src/openai_ns_reconstruction/original_first_collar_jets.py'
TRANSPORT = 'src/openai_ns_reconstruction/original_first_collar_transport.py'


def _add(*pairs):
    return tuple(map(sum, zip(*pairs, strict=True)))


def _scale(a, factor):
    return tuple(value*factor for value in a)


def _sub(a, b):
    return tuple(x-y for x,y in zip(a,b,strict=True))


def _mul(a, b):
    return a[0]*b[0], a[1]*b[0]+a[0]*b[1]


def _div(a, b):
    value = a[0]/b[0]
    return value, (a[1]-value*b[1])/b[0]


def _viscous(T, K, z, d, L, delta):
    return _div(_add(_scale(_mul(z,T[:2]),1-delta), _mul(d,T[1:]), _scale(_mul(z,K),-2)), L)


def _inertial_difference(z, d, L, delta, true, finite):
    # Each tuple stores ordinary jets for U, Q, Q_Z, and rho Q_rho.
    u, q, qz, qr = true
    uf, qf, qfz, qfr = finite
    du, dq, dqz, dqr = tuple(_sub(a,b) for a,b in zip(true,finite,strict=True))

    def transport(q, qz, qr):
        return _div(_add(_scale(_mul(z,_add(q,qr)),-2), _mul(d,qz)), L)

    time = _div(_add(dq,_scale(_mul(z,dqz),(1-delta)/2),dqr),L)
    square = tuple(v/2 for v in _mul(dq,_add(q,qf)))
    radial = _add(_mul(q,dqr),_mul(dq,qfr))
    advection = _add(_mul(du,transport(qf,qfz,qfr)),_mul(u,transport(dq,dqz,dqr)))
    return _add(time,square,radial,advection)


def _inertial(z, d, L, delta, fields):
    u,q,qz,qr = fields
    time = _div(_add(q,_scale(_mul(z,qz),(1-delta)/2),qr),L)
    transport = _div(_add(_scale(_mul(z,_add(q,qr)),-2),_mul(d,qz)),L)
    return _add(time,_mul(q,_add(tuple(v/2 for v in q),qr)),_mul(u,transport))


def difference_density(*, Z, delta, rho_square, rho_square_gap, chi, inv_p2,
                       true_inlet, finite_inlet, actual_T_K, actual_flow, finite_flow,
                       cutoff_product_over_p2):
    """Return signed regular viscous, inertial and cutoff gap C1 jets.

    Averaging these densities over s in [0,S] gives the gap in
    epsilon^2*hb*Pstar^2*S. rho_square_gap is the exact rho^2-16, an
    independent input so a correlated expm1 representation can be retained.
    The inlet pairs are (Psi_rho/Pstar^2 through Z2, Psi_rhorho/Pstar^2
    through Z1); actual_T_K uses the same Pstar^2 normalization. Flow jets
    use the physical U and Q. The cutoff product includes eta and is extended
    by zero outside its true support.
    """
    z, d, L = (Z,1), (1-Z**2,-2*Z), (1-delta*Z**2,-2*delta*Z)

    def initial(jets):
        p, p_rho = jets
        return _scale(p,16), _add(_scale(p[:2],32),_scale(p_rho,64))

    T0, K0 = initial(true_inlet)
    Tf, Kf = initial(finite_inlet)
    T, K = actual_T_K
    radial = _viscous(_sub(T,T0),_sub(K,K0),z,d,L,delta)
    gate = _viscous(T0,K0,z,d,L,delta)
    inlet = _viscous(_sub(T0,Tf),_sub(K0,Kf),z,d,L,delta)
    viscous = _add(_scale(radial,chi),_scale(gate,chi-1),inlet)
    inertial = _add(_scale(_inertial(z,d,L,delta,finite_flow),rho_square_gap),
                    _scale(_inertial_difference(z,d,L,delta,actual_flow,finite_flow),rho_square))
    return dict(viscous=viscous, inertial=_scale(inertial,inv_p2/2),
                cutoff=_scale(cutoff_product_over_p2,-2*rho_square))


def _mass_response(Z, delta):
    """Coefficient of the fixed beta=2 mass in the inlet viscous source."""
    X = Z**2
    x, d, L = 1+X, 1-X, 1-delta*X
    A = 3+delta+(delta-1)*X
    G = x**3*L
    G1 = 3*x**2*L-delta*x**3
    G2 = 6*x*L-6*delta*x**2
    B = A/G
    B1 = (delta-1-B*G1)/G
    B2 = (-2*B1*G1-B*G2)/G
    C = d*(B+2*X*B1)-(3+delta)*X*B
    C1 = -(B+2*X*B1)+d*(3*B1+2*X*B2)-(3+delta)*(B+X*B1)
    return 16*C/L, 32*Z*(C1+delta*C/L)/L


def _function(source, owner, name):
    tree = ast.parse(source)
    scope = next(node.body for node in tree.body if isinstance(node,ast.ClassDef) and node.name==owner) if owner else tree.body
    return next(node for node in scope if isinstance(node,ast.FunctionDef) and node.name==name)


def _assigned(fn, variable, text):
    expected = ast.dump(ast.parse(text,mode='eval').body)
    return any(isinstance(node,ast.Assign) and any(isinstance(target,ast.Name) and target.id==variable for target in node.targets)
               and ast.dump(node.value)==expected for node in ast.walk(fn))


def _source_checks(inputs, pins):
    finite = inputs['inner']
    jets = inputs['jets']
    q2 = inputs['q2']
    whole = inputs['whole']
    canonical = inputs['canonical']
    family = q2['source_family']
    require(finite['checks']['passed'] and jets['checks']['passed'] and q2['checks']['passed'] and whole['passed']
            and canonical['passed'] and family == finite['source_family'] == jets['source_family']
            and all(finite['implementation_hashes'][name] == pins[role]
                    for role,name in (('core_code',CORE),('inner_code',INNER),('jets_code',JETS),('transport_code',TRANSPORT)))
            and jets['implementation_hashes'][CORE] == pins['core_code']
            and jets['implementation_hashes'][JETS] == pins['jets_code']
            and whole['frozen_inputs']['experiments/root_st073/original_collar_regular_pressure_v2_receipt.json.gz']
            == pins['canonical'], 'accepted finite/canonical sources changed')
    owner = 'OriginalN1InnerM15Numeric'
    seed = _function(inputs['core_code'],'OriginalCoreNumeric','_seed')
    prep = _function(inputs['inner_code'],owner,'_prepare')
    cell = _function(inputs['inner_code'],owner,'source_cell')
    integ = _function(inputs['inner_code'],owner,'regular_integrals')
    frame = _function(inputs['inner_code'],None,'leading_regular_frame')
    stencil = _function(inputs['inner_code'],None,'normalized_pressure_source')
    jets_fn = _function(inputs['jets_code'],'OriginalFirstCollarJets','inlet_jets')
    trans = _function(inputs['transport_code'],'OriginalFirstCollarTransport','_inputs')
    nominal_branch = next(node.orelse for node in seed.body if isinstance(node,ast.If)
                          and isinstance(node.test,ast.Name) and node.test.id=='directed')
    zero_branch = ast.FunctionDef(name='_nominal',args=seed.args,body=nominal_branch,decorator_list=[])
    checks = {
        'zero_nominal_swirl_seed': _assigned(zero_branch,'S','[c.mpf(0)]*(n+1)'),
        'only_nominal_beta2_pressure': _assigned(zero_branch,'pressure','[-self.M2_point*q for q in qs]'),
        'same_finite_zero_width_frame': _assigned(prep,'vp',
            "leading_regular_frame(inp['values'],inp['radial'],z,J.constant(c,4,5),J.constant(c,2,5),de,e,P2,J.constant(c,0,5),c.mpf(0),self.o.core.j)"),
        'nominal_seed_coupling': any(isinstance(n,ast.Assign)
            and ast.unparse(n.value)=='self.o.core._seed(Z, 24, 5, directed=False)' for n in prep.body)
            and any(isinstance(n,ast.keyword) and n.arg=='point_S'
                    and ast.unparse(n.value)=="J(c, fixed['S_Z_taylor'][:6])" for n in ast.walk(prep)),
        'nominal_cut_swirl_multiplies_seed': _assigned(cell,'couplingp',
            "(a['point_S']*a['point_frame']['Phi'][0].rows[0]*fp*(2/(self.o.core.epsilon**3*c.exp(2*self.o.core.logP)))).truncate(1)"),
        'nominal_pressure_after_zero_swirl': _assigned(cell,'pressurep',"(-(a['point_Np0'].truncate(1)+couplingp)*16).truncate(1)"),
        'same_nominal_pressure_operator': _assigned(prep,'psource','normalized_pressure_source(vp).rows[0]')
            and 'return -flow / (2 * P2) + (qr * 2 + qrr * r) / (e * P2)' in ast.unparse(stencil),
        'radius_four_derivative_inputs': "radial[name][i] +=" in ast.unparse(jets_fn)
            and 'J, r = (self.Jet, c.mpf(4))' in ast.unparse(jets_fn)
            and _assigned(trans,'inlet','o.inlet_jets(Z,degree=24,order=order)'),
        'same_mean_derivative_recipes': _assigned(frame,'pr',"ratio*chi*comparison['Psi'][1]")
            and _assigned(frame,'mr',"(values['Psi']-values['Mpsi'])/rho")
            and _assigned(frame,'mean',"[Laurent({0:values['Mpsi']}),Laurent({0:mr}),Laurent({0:(pr-2*mr)/rho})]"),
        'fixed_mesh_sum': _assigned(integ,'ends',
            "[iv.mpf(0)]+[iv.mpf(v) for v in ('.5','.75','1','2') if Fraction(v)<Fraction(q)]+[iv.mpf(q)]")
            and _assigned(integ,'width','(hi-lo)/self.cells')
            and _assigned(integ,'wp','sum(_bounds(c,width))/2')
            and 'point[j] += p[j] * wp' in ast.unparse(integ),
    }
    failed = [key for key,value in checks.items() if not value]
    require(not failed, 'finite point recipe changed: '+', '.join(failed))
    ports = canonical['canonical_function_bridge']['ports']
    require(ports['row4_regular_plus_cut'] == [109,111] and ports['cut_swirl'] == [105,107]
            and not ports['legacy_forward_finite_response_equals_all_word'], 'missing physical row-four sources')
    return dict(ast_checks=checks, real_row_four=ports['row4_regular_plus_cut'], real_cut_swirl=ports['cut_swirl'],
                source_functions_numeric_owner_called=False, nominal_mesh_mass_at_q2=2,
                normalized_finite_pressure='-16*N1p_finite(4,Z)/Pstar^2',
                ordinary_Z_orders=[0,1], cutoff_swirl_in_nominal_arithmetic_is_exact_zero=True)


def _algebra_checks():
    z, delta = sp.symbols('Z delta')
    chi,r2,kappa = sp.Rational(3,5),sp.Rational(81,4),sp.Rational(2,7)
    L, d = 1-delta*z**2, 1-z**2
    pf, p2f, mr = (sp.Function(name)(z) for name in ('p_f','p2_f','mean_r'))
    mrr = (pf-2*mr)/4
    qrf = (2*z*pf-(1-delta)*z*mr-d*sp.diff(mr,z))/L
    qrrf = (2*z*p2f-(1-delta)*z*mrr-d*sp.diff(mrr,z))/L
    finite_visc_rhs = 16*(d*sp.diff(pf,z)-(3+delta)*z*pf-8*z*p2f)/L
    mean_error = sp.cancel(-16*(2*qrf+4*qrrf)-finite_visc_rhs)
    require(mean_error == 0 and sp.diff(mean_error,z) == 0,'finite mean elimination changed')

    def jet(x, order=1):
        return (x,sp.diff(x,z),sp.diff(x,z,2)) if order==2 else (x,sp.diff(x,z))

    # Independent polynomial sources for both sides, including nonzero
    # curvature, nonzero radial Q, and nonzero actual cutoff swirl.
    T, K, p, pr = z**2+z+2, z**2-z/2+1, z**3-z+2, z**2+2
    pf, p2f = z**2+z+1, z**3-3*z+2
    q, qf = z**2+2*z+1, z**3-z+2
    actual = (z**3+1,q,sp.diff(q,z),z+4)
    finite = (z+1,qf,sp.diff(qf,z),2*z-3)
    product = z**2+3*z-1
    vals = difference_density(Z=z, delta=delta, rho_square=r2, rho_square_gap=r2-16,
        chi=chi, inv_p2=kappa, true_inlet=(jet(p,2),jet(pr)), finite_inlet=(jet(pf,2),jet(p2f)),
        actual_T_K=(jet(T,2),jet(K)), actual_flow=tuple(map(jet,actual)),
        finite_flow=tuple(map(jet,finite)), cutoff_product_over_p2=jet(product))

    def independent_flow(u,q,qz,qr):
        return (q+(1-delta)*z*qz/2+qr)/L+q*(q/2+qr)+u*(-2*z*(q+qr)+d*qz)/L

    expected = {
        'viscous': chi*((1-delta)*z*T+d*sp.diff(T,z)-2*z*K)/L
                    -16*(d*sp.diff(pf,z)-(3+delta)*z*pf-8*z*p2f)/L,
        'inertial': kappa*(r2*independent_flow(*actual)-16*independent_flow(*finite))/2,
        'cutoff': -2*r2*product,
    }
    count = 2
    for key,pair in vals.items():
        for order, actual in enumerate(pair):
            require(sp.cancel(actual-sp.diff(expected[key],z,order)) == 0, 'signed density algebra changed: '+key)
            count+=1

    pressure_slope = z*(2*d/(1+z**2)**3+(1+delta)/(1+z**2)**2)/L
    mass = 16*(d*sp.diff(pressure_slope,z)-(3+delta)*z*pressure_slope)/L
    for order,row in enumerate(_mass_response(z,delta)):
        require(sp.cancel(row-sp.diff(mass,z,order)) == 0, 'ordinary Z order in mass coefficient changed')
        count+=1
    return count


def evaluate(*, inputs, pins):
    binding = _source_checks(inputs,pins)
    symbolic_checks = _algebra_checks()
    iv, mp = MPIntervalContext(), MPContext()
    iv.prec = 1800
    pressure = inputs['pressure']['compliant_source']
    require(all(pressure[k] == inputs['q2']['source_family'][k] for k in ('implicit_source_sha256','datum_enclosure_sha256'))
            and inputs['whole']['frozen_inputs']['experiments/root_st073/lei_ren_part1_paper_compliant_pressure_source.json']
                == pins['pressure'], 'pressure seed source changed')
    sample, = (row for row in inputs['jets']['samples'] if (row['Z'],row['phase']) == ('0','0'))
    finite_reference_text = sample['pressure']['P0_over_Pstar_squared'][0]['approximation']
    reference = -iv.mpf(finite_reference_text)
    M2 = read(iv,pressure['fixed_beta2_mass_normalized'])
    delta = read(iv,inputs['q2']['actual_selected_leading_pressure_primitive']['delta_enclosure'])
    defect = M2-reference
    require(reference.a > 0 and defect.a > 0 and 0 < delta.a <= delta.b < 1, 'signed pressure mass gap unresolved')

    def encode(value):
        endpoints = value._mpi_
        return dict(lower_exact_mpf_tuple=list(endpoints[0]),upper_exact_mpf_tuple=list(endpoints[1]),
                    rounded=[mp.nstr(mp.make_mpf(bits),17) for bits in endpoints],sign=sign(value))

    def mass_gap(z):
        return _scale(_mass_response(z,delta),defect)

    cells = []
    for i in range(-64,64):
        z = hull(iv,iv.mpf(i)/64,iv.mpf(i+1)/64)
        cells.append(dict(axial_cell=encode(z),positive_beta2_mass_gap_C0_C1=[encode(v) for v in mass_gap(z)]))
    queries = {}
    for z, old in inputs['whole']['legacy_fourth_row_diagnostics'].items():
        response = mass_gap(iv.mpf(z))
        require(old['ledger_count'] == 212 and len(old['ordinary_Z']) == 2,
                'fourth source was not read from an intact 212-row ledger')
        diagnostic = []
        for k, row in enumerate(old['ordinary_Z']):
            old_gap = read(iv,row['canonical_minus_saved_finite'])
            diagnostic.append(dict(ordinary_Z_order=k, full_archived_finite_difference=encode(old_gap),
                signed_mass_sector=encode(response[k]), remaining_point_difference=encode(old_gap-response[k])))
        queries[z] = dict(ledger_sha256=old['ledger_sha256'],ordinary_Z=diagnostic)
    return dict(passed=True,normalized_unit='epsilon^2*hb*Pstar^2*S, S=sc/2',
                functional_difference=binding,exact_algebra_checks=symbolic_checks,
                beta2_reference_is_archived_point_pressure_proxy=True,
                proxy_is_identified_with_finite_producer_coefficient=False,
                finite_reference_pressure_decimal=finite_reference_text,actual_M2_minus_reference_M2=encode(defect),
                pressure_mass_sector_cells=cells, legacy_point_gap_minus_mass_sector=queries,
                complete_defect_enclosure_over_whole_Z=False,
                roundoff_in_forward_values_certified_here=False, inherited_error_packets_removed=False,
                six_original_contributions_preserved=True, forward212_ledger_mutated=False,
                reserved_own_swirl_reapplied=False, independent_terminal_matching_closed=False,
                ancestor_calculations_replayed=False)
