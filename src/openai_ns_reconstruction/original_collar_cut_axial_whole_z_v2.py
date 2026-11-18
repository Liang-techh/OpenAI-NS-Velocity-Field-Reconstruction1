"""Whole-axis C1 enclosure for the actual third m15 contribution."""
from __future__ import annotations

import ast
import hashlib
import json

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext
import sympy as sp

from .original_collar_whole_z import hull, read, require, sign


PREFIX = 'lei_ren_part1_paper_compliant_current_original_n1_'
UNIT = 'epsilon^2*hb*(sc/4)'


def _axis_response(Z,delta,j):
    """C=-2*T_-2*T_(-1-delta)(4Z+j), with ordinary derivatives."""
    L,L1,L2 = 1-delta*Z**2,-2*delta*Z,-2*delta
    d = 1-Z**2
    A = 4-(8+4*delta)*Z**2-(1+delta)*j*Z
    A1,A2 = -2*(8+4*delta)*Z-(1+delta)*j,-2*(8+4*delta)
    T = A/L
    T1 = (A1-T*L1)/L
    T2 = (A2-2*T1*L1-T*L2)/L
    numerator = -2*Z*T+d*T1
    n1 = -2*T-4*Z*T1+d*T2
    C = -2*numerator/L
    C1 = (-2*n1-C*L1)/L
    return C,C1


def _algebra():
    Z,delta,j = sp.symbols('Z delta j')
    D,L = 1-Z**2,1-delta*Z**2
    a = (-(1+delta)*Z*(4*Z+j)+4*D)/L
    direct = -2*(-2*Z*a+D*sp.diff(a,Z))/L
    require(all(sp.cancel(value-sp.diff(direct,Z,k)) == 0
                for k,value in enumerate(_axis_response(Z,delta,j))),
            'axis response or ordinary derivative sign changed')
    X,t = sp.symbols('X t',positive=True)
    ku = t*sp.log(X/t)
    primitive = t*t*sp.log(X/t)/2+t*t/4
    require(sp.simplify(sp.diff(primitive,t)-ku) == 0, 'row-six Green primitive derivative changed')
    require(sp.limit(sp.expand_log(primitive,force=True),t,0,dir='+') == 0
            and sp.simplify(primitive.subs(t,X)-X**2/4) == 0, 'row-six Green endpoint mass changed')
    h,q,r,a,eps,U,C = sp.symbols('hb q rho a epsilon U C')
    physical_axis = 2*eps*h*r*U*(eps*C*r/4)*a
    require(sp.cancel(physical_axis/(eps*eps*h*a) - r*r*U*C/2) == 0,
            'axial term used a wrong physical measure')
    return dict(actual_axis_response_and_ordinary_derivative=True,green_primitive_and_endpoints=True,
                physical_axis_term_in_legacy_unit=True,zero_source_radius_for_axis_model_not_assumed=True)


def _source_guards(inputs,pins):
    q2,domain,natural,coupled,cut = (inputs[k] for k in ('q2','domain','natural','coupled','cut'))
    archive_check = inputs['natural_check']
    family = q2['source_family']
    require(q2['checks']['passed'] and cut['passed'] and archive_check['all_passed']
            and all(bundle['source_family']==family for bundle in (domain,natural,coupled,archive_check))
            and natural['input_hashes'][PREFIX+'inner_analytic_domain.json.gz']==pins['domain']
            and archive_check['input_hashes'][PREFIX+'natural_component_enclosure.json.gz']==pins['natural']
            and archive_check['input_hashes'][PREFIX+'natural_component_enclosure.py']==pins['natural_code'],
            'natural axial source detached from the actual leading domain')
    canonical = q2['actual_first_raw_n1_function']
    history = q2['actual_n1_cutoff_history_identity']
    prefix = 'experiments/root_st073/'+PREFIX
    require(all(canonical['source_pins'][prefix+name]['expected_sha256']==pins[role]
                for role,name in (('raw_code','inner_volterra_functions.py'),('coupled','coupled_first_step.json.gz'))),
            'the raw first word is not pinned by the canonical source')
    require(canonical['passed'] and canonical['all_word_binding']['raw_all_word_convergence_passed']
            and canonical['actual_function_handles']['first']['Uz1'][:2]==cut['raw_all_word_Uz1_C0_C1']
            and history['passed'] and history['cutoff_identity']['ordinary_Z_orders']==[0,1,2,3]
            and cut['leading_U_C0_C1']==q2['actual_native_collar_transport']['canonical_function_handles']['U_C0_C3'][:2]
            and cut['legacy_unit']==UNIT and cut['bounds_of_support_in_s']==['0','3*sc/16']
            and q2['actual_cutoff_complete_support_measure']['exact_unweighted_q_mass']=='1/2+(1/4)*(1/2)=5/8',
            'actual axial source, kernel or cutoff mass changed')

    def source_function(code,name):
        text = inputs[code]
        return ast.unparse(next(value for value in ast.parse(text).body
                                if isinstance(value,ast.FunctionDef) and value.name==name))

    raw = source_function('raw_code','raw_system')
    green = source_function('raw_code','green_apply')
    natural_def = source_function('natural_code','natural_majorants')
    require(
        "g6 = g.sub(g.mul(c(2), eps, Nz), g.quotient(g.mul(c(4), zv, r, eps, eps, Np.value), Lv, 'positive original L'))" in raw
        and "Nz = g.neg(axial_second(g, U, r, de, g.sub(c(-1), de)))" in raw
        and "N1z=derivative.c3(Nz)" in raw
        and "zip(('5', '6', '6', '4', '5', '6'), kernels)" in green
        and "log = g.mul(X, t, g.unary('log', g.quotient(g.one, t, 'positive t inside unit integral')))" in green
        and "log = g.mul(c(2), t, g.add(g.mul(c('1/2'), hb, s), g.unary('log', g.quotient(g.one, t, 'positive t inside unit integral'))))" in green
        and "log = g.mul(c('1/2'), x, distance)" in green and "measure = g.mul(c('1/2'), hb, x)" in green
        and "g6rest = eps(axial_source(rows) * 2) + eps(Np * (4 * z * r / L))" in natural_def
        and "ic = {k: v * 2 for k, v in sc.items()}" in natural_def
        and "ib = {k: algebra.width(v) * (a / 4) for k, v in sb.items()}" in natural_def,
        'the accepted axial source L1 cover no longer matches raw g6/epsilon')
    proof = natural['domain']['uniform_proof']
    cover = natural['domain']['axial_nonaxis_integral_cover']
    require(proof['remaining_axial_Psi_and_pressure_sources_not_set_zero']
            and proof['original_own_mean_and_Q_epsilon_sectors_used']
            and proof['inverse_hb_cancelled_before_norm'] and cover['all_covers_have_directed_log_proof'],
            'no certified nonaxis source on the common phase interval')


def _positive_exp_upper(iv,value):
    # Use only an exact dyadic upper scale. A log of an exponentially small
    # factor can itself be unrepresentably large as a binary integer.
    log2 = value/iv.ln(2)
    if log2.b < -10**20:
        return iv.make_mpf(((0,1,-10**20,1),)*2)
    require(log2.b < 10**20,'normalized log factor unresolved')
    negative,man,exp,_ = log2._mpi_[1]
    integer = -man if negative else man
    upper = integer << exp if exp >= 0 else -((-integer) >> -exp)
    return iv.make_mpf(((0,1,upper,1),)*2)


def _uniform_data(iv,inputs):
    domain = inputs['domain']
    natural_domain = inputs['natural']['domain']
    coupled = inputs['coupled']['samples']['actual_first_keep']['analytic_preconditioned_tail']
    q2 = inputs['q2']
    packet = q2['actual_correlated_q2_inlet']
    require(packet['source_log_basis_names'][:4]==['hb','F0_anchor','Pstar_squared','Lambda'],
            'unexpected raw axial scales')
    log_hb,_,_,log_lambda = (read(iv,v) for v in packet['source_log_bases'][:4])
    a = read(iv,packet['selected_phase'])/2
    eta = read(iv,domain['common_tubes']['outer_Z_radius'])
    rho_cap = read(iv,domain['actual_inner_geometry']['real_rho_upper_bound']).b
    j = read(iv,inputs['core']['required_j'])
    delta = read(iv,q2['actual_selected_leading_pressure_primitive']['delta_enclosure'])
    require(a.a==a.b and 0<a.a<iv.mpf(1)/4 and log_hb.b<0 and log_lambda.a>0
            and 0<eta.a and eta.b<iv.mpf(1)/4
            and domain['joint_core_domain']['Z_real_centers']==['-1','1']
            and domain['joint_core_domain']['Cauchy_rho_radius']=='.5'
            and domain['actual_inner_geometry']['actual_bridge_not_replaced_by_comparison']
            and 4<rho_cap<5 and 0<delta.a<=delta.b<iv.mpf(1)/8 and 0<j.a<=j.b<1,
            'whole-axis tube, cutoff geometry or actual parameters changed')
    proof = natural_domain['uniform_proof']
    require(read(iv,proof['original_outer_eta'])._mpi_==eta._mpi_
            and read(iv,proof['source_M0_tube'])._mpi_==(eta/8)._mpi_
            and read(iv,coupled['solution_Z_radius'])._mpi_==(eta/8)._mpi_
            and coupled['retained_preconditioned_indices']==[0,1]
            and coupled['all_k_derivative_words_allowed']
            and coupled['arbitrary_k_term_formula']=='Mbase*[a*(1+a)*(Cp+C1*k/loss)]^k/k!',
            'nonaxis or raw word norm lost its axial tube')
    gap = eta.a/16
    sqrtr = iv.sqrt(rho_cap).b
    rest = read(iv,natural_domain['axial_nonaxis_integral_cover']['modulus_upper']).b
    require(rest>0,'nonaxis word-zero remainder must retain pressure and axial sources')
    # The raw axial Green kernel t log(X/t) is bounded by X, so a bound
    # on the normalized row-six source L1 gives a direct word-zero bound.
    w0_bound = (sqrtr*rest).b
    mbase, radius, oneJ, Cp, C1, loss, tail, beta = (read(iv,coupled[k]).b for k in (
        'initial_base_common_tube_bound','kinematic_J_norm_upper','I_plus_J_norm_upper',
        'off_kinematic_B0_norm_upper','B1_norm_upper','total_available_Z_loss',
        'omitted_sum_bound','geometric_tail_ratio_upper'))
    # For denominators use their lower exact endpoints. The same bound
    # applies uniformly through the support on the real radial interval.
    loss = read(iv,coupled['total_available_Z_loss']).a
    require(0<loss and 0<beta<1 and 0<tail and 0<Cp and 0<C1 and 0<mbase,'invalid all-word bounds')
    w1_raw = (mbase*radius*oneJ*(Cp+C1/loss)).b
    w1_bound = _positive_exp_upper(iv,iv.ln(w1_raw)+log_lambda)
    tail_bound = _positive_exp_upper(iv,iv.ln(tail)+log_lambda)

    # Reuse the core degree-independent norm and the leading radial ODE:
    # |Psi|<=4Bpsi at the inlet, |Psi_rho|<=8Bpsi*|Phi/Phi_core|
    # throughout the actual real radial collar. Both bounds hold uniformly
    # on the admitted axial complex disks.
    Bpsi = read(iv,domain['joint_core_domain']['Xh_norm_majorants']['Psi']).b
    Phi = read(iv,domain['actual_collar_radial_bounds']['Phi'][0]['0']['coefficient']).b
    floor = read(iv,domain['denominator_bounds']['Phi_core_modulus_lower']).a
    require(0<Bpsi and 0<Phi and 0<floor,'zero core or collar sup norm')
    leading = _positive_exp_upper(iv,iv.ln(Bpsi*(4+4*rho_cap*Phi/floor))-log_lambda)
    hb = _positive_exp_upper(iv,log_hb)
    width = (3*hb*a/2).b
    require(0<width<iv.mpf('0.5'),'cutoff collar radius ratio unresolved')
    # Integral of rho^2*eta/2 in q uses rho^2=16 exp(2 hb a q),
    # q<=3/4 and the exact unweighted mass 5/8.
    axis_mass = hull(iv,iv.mpf(5),(5/(1-width)).b)
    return dict(delta=delta,j=j,axis_mass=axis_mass,gap=gap,rho_cap=rho_cap,
                radius_a=a,log_hb=log_hb,log_lambda=log_lambda,
                field_radii={'leading_nonaxis':(leading,(leading/gap).b),
                    'word0_nonaxis':(w0_bound,(w0_bound/gap).b),
                    'word1':(w1_bound,(w1_bound/gap).b),
                    'omitted_after_word1':(tail_bound,(tail_bound/gap).b)})


def _integrate_cell(iv,data,Z):
    delta,j = data['delta'],data['j']
    u0 = 4*Z+j
    C,C1 = _axis_response(Z,delta,j)
    radii = data['field_radii']
    dU,dUZ = radii['leading_nonaxis']
    rho = data['rho_cap']
    radial_weight = 5*rho/4

    # The two signed factors in the axis term share the same geometry;
    # interval multiplication by its positive weighted integral is exact.
    axis = (data['axis_mass']*u0*C, data['axis_mass']*(4*C+u0*C1))
    u = abs(u0).b+dU
    uz = 4+dUZ
    leading = (radial_weight*dU*rho*abs(C).b/4,
               radial_weight*rho*(dUZ*abs(C).b+dU*abs(C1).b)/4)
    errors = {'leading_nonaxis_times_axis_word0':tuple(v.b for v in leading)}
    for name in ('word0_nonaxis','word1','omitted_after_word1'):
        e0,e1 = radii[name]
        errors[name] = ( (radial_weight*u*e0).b, (radial_weight*(uz*e0+u*e1)).b)
    totals = [sum(v[k] for v in errors.values()).b for k in range(2)]
    return axis,errors,tuple(v+hull(iv,-r,r) for v,r in zip(axis,totals,strict=True))


def _archived_comparison(iv,inputs,data,encode,samples):
    quiet,inner,whole = (inputs[k] for k in ('quiet','inner','row4_whole'))
    family = inputs['q2']['source_family']
    require(quiet['checks']['passed'] and inner['checks']['passed'] and whole['passed']
            and quiet['source_family']==inner['source_family']==family,'forward archive family changed')
    result = {}
    exact_log_factors = (-2*data['log_lambda'],data['log_hb'],iv.ln(data['radius_a']))
    for old_Z,Z in (('-.3','-0.3'),('0','0'),('.3','0.3')):
        terminal, = (row for row in quiet['samples']
                     if (row['Z'],row['chart'],row['native_coordinate'])==(old_Z,'quiet','1'))
        ledger = terminal['m15_C0_C1_terms']
        digest = hashlib.sha256(json.dumps(ledger,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        parents = [row for row in inner['samples'] if (row['Z'],row['q'],row['cells_per_piece'])==(old_Z,'2',16)
                   and ledger[:6]==row['m15_C0_C1_terms']]
        require(len(ledger)==212 and parents and digest==whole['legacy_fourth_row_diagnostics'][Z]['ledger_sha256'],
                'old source is no longer the third contribution of the saved full ledger')
        parent = parents[0]
        require(tuple(parent['exact_phase_tuple'])==data['radius_a']._mpi_[0],
                'old axial contribution uses a different chart width')
        row = ledger[2]
        binding = row['actual_function_binding']
        require(row['name']=='first_lower_axial_width' and len(row['log_scale_terms'])==3
                and row['log_scale_terms'][1].get('source_log_atom_only') is True
                and binding['density']=='2*rho*U_lead*u_cut' and binding['measure']=='hb*(sc/4)*dq'
                and binding['own_cutoff']=='1-sigma(4q-2)',
                'wrong archived axial density or cutoff')
        log_overlap = []
        for item,actual in zip(row['log_scale_terms'],exact_log_factors,strict=True):
            archived = iv.make_mpf(tuple(tuple(bits) for bits in item['source_enclosure_exact_mpf_tuples']))
            log_overlap.append(bool(actual.a<=archived.b and archived.a<=actual.b))
        require(all(log_overlap),'incompatible recorded scale enclosures')
        values = []
        for k,(item,v) in enumerate(zip(row['coefficient_C0_C1'],samples[Z],strict=True)):
            enclosure = iv.make_mpf(tuple(tuple(bits) for bits in item['source_enclosure_exact_mpf_tuples']))
            model = iv.mpf(item['approximation'])
            values.append(dict(order=k,finite_approximation=encode(model),archived_enclosure=encode(enclosure),
                actual_minus_stored_finite=encode(v-model),archived_enclosure_contains_actual=bool(
                    enclosure.a<=v.a and v.b<=enclosure.b),overlap=bool(enclosure.a<=v.b and v.a<=enclosure.b)))
        result[Z] = dict(ledger_count=len(ledger),ledger_sha256=digest,raw_third_row=row,
                         same_physical_unit=UNIT,paired_scale_interval_overlap=log_overlap,ordinary_Z=values)
    return result


def evaluate(*,inputs,pins):
    _source_guards(inputs,pins)
    proof = _algebra()
    iv,mp = MPIntervalContext(),MPContext()
    iv.prec = 1800
    data = _uniform_data(iv,inputs)

    def encode(value):
        return dict(lower_exact_mpf_tuple=list(value._mpi_[0]),upper_exact_mpf_tuple=list(value._mpi_[1]),
                    rounded=[mp.nstr(mp.make_mpf(bits),17) for bits in value._mpi_],sign=sign(value))

    def at(z):
        axis,parts,actual = _integrate_cell(iv,data,z)
        return actual,dict(axial_cell=encode(z),signed_axis_response_C0_C1=[encode(v) for v in axis],
            error_radii={k:[encode(v) for v in pair] for k,pair in parts.items()},
            actual_axial_response_C0_C1=[encode(v) for v in actual])

    cells = [at(hull(iv,iv.mpf(k)/128,iv.mpf(k+1)/128))[1] for k in range(-128,128)]
    saved = {}
    samples = {}
    for label in ('-0.3','0','0.3','0.371'):
        samples[label],saved[label] = at(iv.mpf(label))
    diagnostics = _archived_comparison(iv,inputs,data,encode,samples)
    previous = inputs['cut']['selected_axial_enclosure']
    require(previous['Z']=='371/1000' and previous['unit']==UNIT,'wrong saved axial point')
    selected = []
    for k,new in enumerate(samples['0.371']):
        old = read(iv,previous['ordinary_Z'][str(k)]['actual_interval'])
        selected.append(dict(order=k,previous_contains_new=bool(old.a<=new.a and new.b<=old.b),
                             intervals_overlap=bool(old.a<=new.b and new.a<=old.b)))
    require(all(v['intervals_overlap'] for v in selected)
            and all(v.a>0 for v in samples['0.371']) and all(v.a>0 for v in samples['0']),
            'actual axial model conflicts with saved sign or phase')
    return dict(passed=True, same_physical_unit=UNIT, real_Z_domain=['-1','1'], ordinary_Z_orders=[0,1],
                exact_algebra=proof, exact_axis_formula='C=-2*T_-2*T_(-1-delta)(4*Z+j); U1_axis/epsilon=C*rho/4',
                raw_word_zero_nonaxis_is_independently_bounded=True,
                raw_word1_bounded_separately=True, omitted_raw_words_start_at=2,
                field_error_radii_C0_C1={name:[encode(v) for v in pair] for name,pair in data['field_radii'].items()},
                axis_geometry_factor=encode(data['axis_mass']), cells=cells, point_queries=saved,
                selected_point_crosscheck=selected, legacy_forward_third_row_comparison=diagnostics,
                finite_legacy_coefficients_equal_canonical_functions=False,
                archived_forward_terms_mutated=False, inherited_error_packets_removed=False,
                independent_terminal_matching_closed=False, ancestor_calculations_replayed=False)
