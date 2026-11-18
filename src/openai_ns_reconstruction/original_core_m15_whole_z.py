"""Bound the first two actual m15 contributions on the full real axial interval."""
from __future__ import annotations

import ast
import hashlib
import json
from math import factorial

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext
import sympy as sp

from .original_collar_whole_z import _dyadic_inv_p2, hull, read, require, sign
from .original_collar_cut_axial_whole_z_v2 import _positive_exp_upper
from .original_signed_density_expression import SignedDensityExpressionGraph


PREFIX = 'lei_ren_part1_paper_compliant_current_original_n1_'
CUTOFF = PREFIX+'inner_cutoff_functions.json.gz'
UNIT = ('epsilon^2','epsilon^2*Pstar^2')


def _prod(a,b):
    return a[0]*b[0],a[0]*b[1]+a[1]*b[0]


def _quot(a,b):
    value = a[0]/b[0]
    return value,(a[1]-value*b[1])/b[0]


def _sum(*terms):
    return tuple(map(sum,zip(*terms,strict=True)))


def _scale(pair,x):
    return tuple(x*v for v in pair)


def _axis_inertia(Z,delta,j):
    u,u1 = 4*Z+j,4
    L,L1,L2 = 1-delta*Z**2,-2*delta*Z,-2*delta
    d,d1 = 1-Z**2,-2*Z
    N = (1+delta)*Z*u-4*d
    N1 = (1+delta)*(u+4*Z)+8*Z
    N2 = 8*(2+delta)
    Q = N/L
    Q1 = (N1-Q*L1)/L
    Q2 = (N2-2*Q1*L1-Q*L2)/L
    Qjet,djet,Ljet,zjet = (Q,Q1),(d,d1),(L,L1),(Z,1)
    B = _quot(_sum(_scale(_prod(zjet,Qjet),-2),_prod(djet,(Q1,Q2))),Ljet)
    time = _quot(_sum(Qjet,_scale(_prod(zjet,(Q1,Q2)),(1-delta)/2)),Ljet)
    return _sum(time,tuple(x/2 for x in _prod(Qjet,Qjet)),
                _prod((u,u1),B))


def _algebra():
    r,z,eps,delta,j = sp.symbols('rho Z epsilon delta j')
    b = sp.Function('b')(r,z)
    S = sp.Function('S')(r,z)
    co = -eps**2*r*(-S/2+2*sp.diff(b,r)+r*sp.diff(b,r,2))
    ftc = sp.diff(-eps**2*r**2*sp.diff(b,r),r)+eps**2*r*S/2
    require(sp.expand(co-ftc)==0 and sp.expand(sp.diff(co-ftc,z))==0,'core pressure FTC changed')
    p,m = sp.symbols('p mean_rho')
    d,L = 1-z**2,1-delta*z**2
    dp,md,dd = sp.symbols('p_Z mean_rho_Z mean_Z_defect')
    b0 = (2*z*p-(1-delta)*z*m-d*md)/L
    model = ((3+delta)*z*p-d*dp)/(2*L)-((1-delta)*z*(m-p/2)+d*dd)/L
    require(sp.cancel(b0.subs(md,dp/2+dd)-model)==0,'own mean endpoint difference changed')
    t = sp.symbols('t',real=True)
    require(sp.integrate(4*t*(1-t),(t,0,1))==sp.Rational(2,3),'wrong mean correction weight')
    require(sp.integrate(2*r/4,(r,0,4))==4 and sp.integrate(r/2,(r,0,4))==4,
            'wrong core axial or inertia volume weight')
    u = 4*z+j
    Q = ((1+delta)*z*u-4*d)/L
    q1 = sp.diff(Q,z)
    expectation = (Q+(1-delta)*z*q1/2)/L+Q**2/2+u*(-2*z*Q+d*q1)/L
    require(all(sp.cancel(actual-sp.diff(expectation,z,k))==0
                for k,actual in enumerate(_axis_inertia(z,delta,j))),'axis inertia derivative changed')
    return dict(regular_pressure_FTC_C0_C1=True,own_mean_endpoint_reduction=True,
                actual_weight_2_over_3=True,core_axial_and_inertia_weights=True,
                axis_inertia_C0_C1=True)


def _source_graph(inputs,pins):
    lower,cut,q2 = (inputs[k] for k in ('lower','cutoff','q2'))
    old,regular = inputs['row3'],inputs['regular']
    history = q2['actual_n1_cutoff_history_identity']
    certified = q2['actual_core_pressure_m15_primitive']
    raw = q2['actual_first_raw_n1_function']['actual_function_handles']
    family = q2['source_family']
    require(lower['source_family']==cut['source_family']==family and history['passed']
            and certified['actual_core_pressure_m15_functional_identity_bound']
            and certified['actual_canonical_core_endpoint_handle_bound']
            and regular['passed'] and old['passed']
            and q2['input_hashes'][CUTOFF]==pins['cutoff']
            and q2['input_hashes'][PREFIX+'lower_collar_functions.py']==pins['lower_code'],
            'wrong core, lower or cutoff sources')
    original = cut['exact_graph_nodes']
    require(original[:len(lower['exact_graph_nodes'])]==lower['exact_graph_nodes'],
            'the collar does not extend the exact core graph')
    low = lower['lower_collar_bindings']
    p,f,radial = low['parameters'],low['fields'],low['core_radial_rows']
    text = inputs['lower_code']
    tree = ast.parse(text)
    definition = ast.unparse(next(node for node in tree.body if isinstance(node,ast.FunctionDef)
                                 and node.name=='collar_fields'))
    require("psi0 = a.add(at(g, core_rows['Psi'][0], rho_variable, c(4))" in definition
            and "psi = [psi0, a.mul(ratio, a.mul(X, pc[1]))" in definition
            and "a.scale(at(g, core_rows['mean_Psi'][0], rho_variable, c(4)), 4)" in definition
            and "mean1 = a.div(a.add(psi0, a.neg(mean0)), a.fixed(rs))" in definition
            and "Q = [own.radial_Q(g, de, u, mm) for u, mm in zip(U, m)]" in definition
            and "rho4_inlet_ratio_exact_one=True" in definition
            and low['geometry_certificate']['sc_and_hb_Z_independent']
            and f['all_histories_start_at_true_core_rho4']
            and certified['canonical_actual_core_endpoint_handle']['requested']==[53125,53126]
            and history['graph_binding']['history_nodes']['rho4_port']==[53125,53126],
            'core and lower limits no longer reuse the true inlet')
    rawU = raw['core']['Uz1'][:2]
    require(all(original[node]['operation']=='current_original_unique_regular_n1_component'
                and original[node]['quantity']=='Uz1' and original[node]['chart']=='core'
                and original[node]['coordinate']==raw['core_coordinate']
                and original[node]['ordinary_slow_Z_derivative_order']==k
                for k,node in enumerate(rawU)), 'wrong core all-word axial functions')
    g = SignedDensityExpressionGraph()
    def ref(node):
        return g.reference(CUTOFF,pins['cutoff'],node,source_nodes=original)
    def at(function,variable,point):
        return g.insert(dict(operation='function_substitution',expression=function.node,
                             variable=variable.node,value=point.node))
    def node_pair(pair):
        return [v.node for v in pair]
    eps,delta,j,Z = (ref(i) for i in (p['epsilon'],p['delta'],p['j'],20))
    r = ref(raw['core_coordinate'])
    four,d,L = g.constant(4),1-Z**2,1-delta*Z**2
    p0,p1 = (at(ref(node),r,four) for node in radial['Psi'][1][:2])
    m0,m1,m2 = (at(ref(node),r,four) for node in radial['mean_Psi'][1][:3])
    # The core derivative is the collar derivative at s=0 by the exact
    # shared inlet and own-mean identity. One E0 serves both orientations.
    b0 = (2*Z*p0-(1-delta)*Z*m0-d*m1)/L
    b1 = (2*(p0+Z*p1)-(1-delta)*(m0+Z*m1)+2*Z*m1-d*m2+2*delta*Z*b0)/L
    E0 = (16*eps**2*b0,16*eps**2*b1)
    s,S = ref(f['query']),ref(p['sc_half'])
    rho_S = at(ref(f['rho']),s,S)
    ES = tuple(eps*rho_S**2*at(ref(node),s,S) for node in f['Q_rho_rows'][1][:2])
    core_viscous,collar_viscous,paired = tuple(-e for e in E0),tuple(a-b for a,b in zip(E0,ES,strict=True)),tuple(-e for e in ES)

    U,UZ = (4*Z+j+eps*ref(radial['Psi'][0][0]),4+eps*ref(radial['Psi'][0][1]))
    u,uZ = map(ref,rawU)
    var = g.insert(dict(operation='bound_variable',name='m15_core_rho'))
    core_axial = tuple(g.insert(dict(operation='definite_integral',
        integrand=at(body,r,var).node,lower=g.zero.node,upper=four.node,variable=var.node,
        ordinary_Z_order=k,endpoint_independent_of_Z=True)) for k,body in enumerate((2*eps*U*u,2*eps*(UZ*u+U*uZ))))
    core_all = tuple(map(ref,(53125,53126)))
    core_pressure = tuple(x-y for x,y in zip(core_all,core_axial,strict=True))
    return dict(**g.packet(),core_axial_C0_C1=node_pair(core_axial),core_pressure_C0_C1=node_pair(core_pressure),
                E0_shared_C0_C1=node_pair(E0),ES_C0_C1=node_pair(ES),
                core_viscous_signed_C0_C1=node_pair(core_viscous),
                collar_viscous_signed_C0_C1=node_pair(collar_viscous),
                paired_viscous_signed_C0_C1=node_pair(paired),
                exact_shared_source_identity='core:-E0, collar:+E0-ES, same E0 at rho=4; sum=-ES',
                remaining_collar_source_positions=[3,4,5,6],
                linked_regular_collar_receipt=pins['regular'],
                linked_row3_receipt=pins['row3'])


def _data(iv,inputs):
    q2, core, tube = (inputs[k] for k in ('q2','core','tube'))
    packet = q2['actual_correlated_q2_inlet']
    signed = q2['actual_signed_C_inner_m15']
    row3 = inputs['row3']
    family = q2['source_family']
    require(packet['passed'] and signed['passed'] and signed['source_family']==family
            and signed['actual_axis_Psi_mean_Qpsi_zero_receipt_and_source_checked']
            and signed['axis_endpoint_rho_squared_b_rho_limit_zero']
            and signed['core_mean_not_reset'] and signed['core_source_handles']==[53125,53126]
            and signed['raw_core_bounds_proof']['exact_core_source_and_solution_Z_tubes_match']
            and signed['raw_core_bounds_proof']['whole_core_radial_source_disk_checked']
            and signed['core_bounds_proof']['core_radial_interval']==[0,4]
            and signed['core_bounds_proof']['no_inverse_rho_bound_used']
            and core['contraction_proved'] and core['analytic_core_family_sha256']==tube['analytic_core_family_sha256']
            and all(core[k]==family[k] for k in ('implicit_source_sha256','datum_enclosure_sha256')),
            'same actual whole-core pressure and swirl bounds required')
    require(row3['raw_word_zero_nonaxis_is_independently_bounded']
            and row3['raw_word1_bounded_separately'] and row3['omitted_raw_words_start_at']==2
            and packet['source_log_basis_names'][:4]==['hb','F0_anchor','Pstar_squared','Lambda'],
            'core response lost separately admitted raw terms')
    loghb,_,logp2,loglam = (read(iv,v) for v in packet['source_log_bases'][:4])
    invp2 = _dyadic_inv_p2(iv,read(iv,core['logPstar']))
    delta = read(iv,q2['actual_selected_leading_pressure_primitive']['delta_enclosure'])
    j = read(iv,core['required_j'])
    S = read(iv,packet['selected_phase'])
    require(S.a==S.b and 0<S.a<iv.mpf(1)/4 and 0<delta.a<=delta.b<iv.mpf(1)/8
            and 0<j.a<=j.b<1 and loghb.b<0, 'core parameters or cutoff width changed')
    hb = hull(iv,iv.mpf(0),_positive_exp_upper(iv,loghb))
    sigma = hull(iv,iv.mpf(0),read(iv,packet['sigma_S_atom']['positive_cap']).b)
    require(sigma.b>0 and packet['sigma_S_atom']['sigma_is_strictly_positive'],
            'missing positive nonmaterialized flat atom')
    norm = read(iv,core['scaled_map_size_upper'])/(1-read(iv,core['scaled_map_Lipschitz_upper']))
    h = read(iv,tube['Xh_parameter'])
    # The linear core model has zero rho^2 derivative. Drop (n+1)^-2
    # from the accepted X_h embedding and integrate the difference between
    # Psi_rho(4t) and Psi_rho(4) against t dt to bound the own-mean error.
    mean_defect = []
    for k in range(3):
        bound = invp2*norm*factorial(2+k)/(600*(k+1)**2*h**k*(iv.mpf(4)/5)**(3+k))
        mean_defect.append(bound.b)
    magnitudes = {}
    for key in ('core_inertial_error','core_swirl_coupling_unsigned_error'):
        values=[]
        for order in ('0','1'):
            radius=iv.mpf(0)
            for sector in signed['ordinary_Z'][order][key]:
                hp,fp,pp,lp,sig = sector['source_exponents']
                coeff = read(iv,sector['coefficient_interval'])
                require(hp==fp==sig==0 and -4<=lp<=-2 and 0<=pp<=2
                        and coeff.a>=0 and coeff.b>0 and sector['coefficient_interval']['directed'],
                        'unexpected core remainder power or signed remainder value')
                # These exact polynomial bounds are uniform in Z, even
                # though the old report displayed a selected signed model.
                radius += _positive_exp_upper(iv,iv.ln(coeff.b)+(pp-1)*logp2+(lp+2)*loglam)
            values.append(radius.b)
        magnitudes[key]=tuple(values)
    radii={name:tuple(read(iv,r).b for r in pair) for name,pair in row3['field_error_radii_C0_C1'].items()}
    axis_mass=read(iv,row3['axis_geometry_factor'])
    require(0<axis_mass.a<=axis_mass.b, 'the exact positive third-row axis weight is required')
    return dict(invp2=invp2,delta=delta,j=j,S=S,hb=hb,sigma=sigma,
                cut_axis_mass=axis_mass,
                loglam=loglam,logp2=logp2, mean_defect=tuple(mean_defect),
                axial_radii=radii,pressure_radii=magnitudes)


def _core_cells(iv,data,Z,inlet,archived_axis):
    delta,j = data['delta'],data['j']
    U = 4*Z+j
    axis = _scale(archived_axis,4/data['cut_axis_mass'])
    # On the radius 1/4 capsule |T_(-1-delta)u0|<=112. A 1/8
    # Cauchy loss bounds its derivative by 896, so |C|<=19712 on
    # radius 1/8 and |C_Z|<=8*19712 on the real axis.
    C,C1 = 19712,8*19712
    source = data['axial_radii']
    eu,eu1 = source['leading_nonaxis']
    axial_radii = {'leading_nonaxis_times_axis_word0':( (4*eu*C).b,
                                      (4*(eu1*C+eu*C1)).b)}
    Umax,U1max = abs(U).b+eu,4+eu1
    for key in ('word0_nonaxis','word1','omitted_after_word1'):
        e,e1 = source[key]
        axial_radii[key] = ( (8*Umax*e).b, (8*(Umax*e1+U1max*e)).b)
    axial = tuple(a+hull(iv,-r,r) for a,r in zip(axis,(sum(x[i] for x in axial_radii.values()).b
                                                      for i in range(2)),strict=True))
    dm,dm1,dm2 = tuple(hull(iv,-v,v) for v in data['mean_defect'])
    p,p1,p2 = inlet
    d,d1 = 1-Z**2,-2*Z
    L,L1 = 1-delta*Z**2,-2*delta*Z
    # Keep the common p/2 part of the own mean correlated analytically.
    numerator = ((3+delta)*Z*p-d*p1)/2-(1-delta)*Z*dm-d*dm1
    n1 = ((3+delta)*(p+Z*p1)-d1*p1-d*p2)/2-(1-delta)*(dm+Z*dm1)-d1*dm1-d*dm2
    b,b1 = _quot((numerator,n1),(L,L1))
    viscous = -16*b,-16*b1
    axis_inertial = _scale(_axis_inertia(Z,delta,j),4*data['invp2'])
    radii = data['pressure_radii']
    pressure = _sum(viscous,axis_inertial,tuple(hull(iv,-t,t) for t in
                     (sum(pair[i] for pair in radii.values()).b for i in range(2))))
    return dict(axis=axial,pressure=pressure,axis_model=axis,axis_errors=axial_radii,
                viscous=viscous,axis_inertia=axis_inertial,E0_normalized=(-viscous[0],-viscous[1]))


def _archived(iv,inputs,data,queries,encode):
    raw_report,inner,previous = (inputs[k] for k in ('quiet','inner','row3'))
    family = inputs['q2']['source_family']
    require(raw_report['checks']['passed'] and inner['checks']['passed']
            and raw_report['source_family']==inner['source_family']==family,
            'wrong forward family')
    diagnostics={}
    row_names = ('core_axial_coupling','core_pressure_integral_and_finite_boundary')
    for src_z,z in (('-.3','-0.3'),('0','0'),('.3','0.3')):
        term, = (r for r in raw_report['samples'] if (r['Z'],r['chart'],r['native_coordinate'])==(src_z,'quiet','1'))
        ledger = term['m15_C0_C1_terms']
        dig = hashlib.sha256(json.dumps(ledger,sort_keys=True,separators=(',',':')).encode()).hexdigest()
        parent = [r for r in inner['samples'] if (r['Z'],r['q'],r['cells_per_piece'])==(src_z,'2',16)
                  and ledger[:6]==r['m15_C0_C1_terms']]
        require(len(ledger)==212 and parent
                and dig==previous['legacy_forward_third_row_comparison'][z]['ledger_sha256'],
                'old first pair lost its q2 source and exact ledger')
        comparisons=[]
        for k,row in enumerate(ledger[:2]):
            require(row['name']==row_names[k] and len(row['log_scale_terms'])==k+1,
                    'wrong old core moment or pressure scale')
            atoms = (-2*data['loglam'],data['logp2'])[:k+1]
            for saved,actual in zip(row['log_scale_terms'],atoms,strict=True):
                old = iv.make_mpf(tuple(tuple(x) for x in saved['source_enclosure_exact_mpf_tuples']))
                require(old.a<=actual.b and actual.a<=old.b, 'old physical unit mismatch')
            result=[]
            for item,actual in zip(row['coefficient_C0_C1'],queries[z][('axis','pressure')[k]],strict=True):
                old = iv.make_mpf(tuple(tuple(x) for x in item['source_enclosure_exact_mpf_tuples']))
                finite = iv.mpf(item['approximation'])
                result.append(dict(actual_minus_stored_finite=encode(actual-finite),
                    archived_enclosure_contains_actual=bool(old.a<=actual.a and actual.b<=old.b),
                    intersects=bool(old.a<=actual.b and actual.a<=old.b)))
            comparisons.append(dict(unit=UNIT[k],ordinary_Z=result,raw_undivided_row=row))
        diagnostics[z]=dict(ledger_count=len(ledger),ledger_sha256=dig,core_rows=comparisons)
    return diagnostics


def _common_six(iv,data,Z,fields,pressure_cell,third,fifth_remainder):
    invp2,a = data['invp2'],data['S']/2
    hb,sigma,delta = data['hb'],data['sigma'],data['delta']
    p0,p1 = (read(iv,r) for r in pressure_cell['inlet_ordinary_C0_C1_C2'][:2])
    principal = _scale(_quot((32*Z*p0,32*(p0+Z*p1)),(1-delta*Z**2,-2*delta*Z)),sigma)
    row4 = tuple(read(iv,r) for r in pressure_cell['row_4_ordinary_C0_C1'])
    return ( _scale(fields['axis'],invp2),fields['pressure'],_scale(third,hb*a*invp2),
             _scale(row4,hb*data['S']),principal,_scale(fifth_remainder,hb*sigma))


def evaluate(*,inputs,pins):
    source = _source_graph(inputs,pins)
    proof = _algebra()
    iv,mp = MPIntervalContext(),MPContext()
    iv.prec = 1800
    data = _data(iv,inputs)
    row3,row4,row6 = (inputs[n] for n in ('row3','row4','row6'))
    require(all(report['passed'] for report in (row3,row4,row6))
            and len(row3['cells'])==256 and len(row4['cells'])==128
            and len(row6['cells'])==64 and row6['covers_every_real_Z']
            and row6['units']['primitive']=='epsilon^2*hb*Pstar^2*sigma(S)',
            'same actual axial cells and pressure-remainder units required')

    def encode(value):
        return dict(lower_exact_mpf_tuple=list(value._mpi_[0]),upper_exact_mpf_tuple=list(value._mpi_[1]),
                    rounded=[mp.nstr(mp.make_mpf(bits),17) for bits in value._mpi_],sign=sign(value))

    def entry(Z,pressure_cell,third,axis_base,remainder):
        true_inlet = tuple(read(iv,r) for r in pressure_cell['inlet_ordinary_C0_C1_C2'])
        fields = _core_cells(iv,data,Z,true_inlet,axis_base)
        common = _common_six(iv,data,Z,fields,pressure_cell,third,remainder)
        packed = {key:[encode(v) for v in fields[key]] for key in
                  ('axis','pressure','axis_model','viscous','axis_inertia','E0_normalized')}
        packed['axis_errors'] = {key:[encode(v) for v in pair] for key,pair in fields['axis_errors'].items()}
        packed['axial_cell'] = encode(Z)
        packed['actual_q2_six_in_common_pressure_unit'] = [[encode(v) for v in row] for row in common]
        packed['actual_q2_total_common_pressure_unit'] = [encode(v) for v in _sum(*common)]
        return fields,packed

    cells=[]
    for i,pressure_cell in enumerate(row4['cells']):
        Z=read(iv,pressure_cell['axial_cell'])
        left,right = row3['cells'][2*i:2*i+2]
        c0,c1 = (read(iv,r['axial_cell']) for r in (left,right))
        rc = row6['cells'][i//2]
        cover = read(iv,rc['axial_cell'])
        require(c0.a==Z.a and c0.b==c1.a and c1.b==Z.b
                and cover.a<=Z.a and Z.b<=cover.b,'misaligned closed axial cells')
        def merged(key):
            pair=[]
            for a,b in zip(left[key],right[key],strict=True):
                av,bv = read(iv,a),read(iv,b)
                pair.append(hull(iv,min(av.a,bv.a),max(av.b,bv.b)))
            return pair
        rem=tuple(read(iv,r) for r in rc['ordinary_C0_C1']['primitive'])
        cells.append(entry(Z,pressure_cell,merged('actual_axial_response_C0_C1'),
                           merged('signed_axis_response_C0_C1'),rem)[1])

    points,queries = {},{}
    for Z in ('-0.3','0','0.3','0.371'):
        third=tuple(read(iv,r) for r in row3['point_queries'][Z]['actual_axial_response_C0_C1'])
        axis_base=tuple(read(iv,r) for r in row3['point_queries'][Z]['signed_axis_response_C0_C1'])
        rem=tuple(read(iv,r) for r in row6['exact_samples'][Z]['ordinary_C0_C1']['primitive'])
        queries[Z],points[Z]=entry(iv.mpf(Z),row4['samples'][Z],third,axis_base,rem)
    actual = _archived(iv,inputs,data,queries,encode)
    selected=[]
    for k in range(2):
        previous = inputs['q2']['actual_signed_C_inner_m15']['ordinary_Z'][str(k)]
        previous = read(iv,previous['actual_core_m15']['normalized']['normalized_integral_enclosure'])
        axis,pressure = queries['0.371']['axis'][k],queries['0.371']['pressure'][k]
        now = axis*data['invp2']+pressure
        selected.append(dict(order=k,directed_enclosure=encode(now),
                             overlaps_earlier_directed=bool(now.a<=previous.b and previous.a<=now.b)))
    require(all(r['overlaps_earlier_directed'] for r in selected),'mistimed earlier core point or source')
    return dict(passed=True, actual_graph=source, exact_algebra=proof,
        real_Z_domain=['-1','1'], ordinary_Z_orders=[0,1], separate_core_units=list(UNIT),
        numeric_cells=cells, point_queries=points, selected_core_crosscheck=selected,
        archived_first_two_row_comparisons=actual,
        core_pressure_volume_radii={key:[encode(v) for v in pair] for key,pair in data['pressure_radii'].items()},
        own_mean_rho_minus_half_Psi_rho_ordinary_C0_C1_C2_radius=[encode(v) for v in data['mean_defect']],
        signed_six_same_source_core_and_lower_cover=True,
        true_positive_flat_atom_used_as_directed_enclosure_only_in_numeric_rows=True,
        cumulative_rows_appended_or_removed=0, inherited_error_packets_removed=False,
        legacy_finite_approximation_is_exact_function=False, independent_terminal_matching_closed=False,
        ancestor_calculations_replayed=False)
