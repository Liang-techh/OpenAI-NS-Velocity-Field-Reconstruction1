"""Keep the common beta=2 pressure mass together in the core endpoint."""
from math import factorial

from mpmath.ctx_iv import MPIntervalContext
from mpmath.ctx_mp import MPContext
import sympy as sp

from .original_collar_whole_z import _dyadic_inv_p2, hull, read, require, sign


UNIT = 'epsilon^2*Pstar^2'


def _coefficient(Z,delta):
    x = Z**2
    u,d,L = 1+x,1-x,1-delta*x
    A = 3+delta+(delta-1)*x
    G = u**3*L
    Gx,Gxx = 3*u**2*L-delta*u**3,6*u*L-6*delta*u**2
    B = A/G
    Bx = (delta-1-B*Gx)/G
    Bxx = (-2*Bx*Gx-B*Gxx)/G
    w = d*(B+2*x*Bx)-(3+delta)*x*B
    wx = -(B+2*x*Bx)+d*(3*Bx+2*x*Bxx)-(3+delta)*(B+x*Bx)
    return 8*w/L,16*Z*(wx+delta*w/L)/L


def _pair_proof():
    z,delta = sp.symbols('Z delta')
    d,L,x = 1-z**2,1-delta*z**2,1+z**2
    p = z*(2*d/x**3+(1+delta)/x**2)/L
    exact = 8*(d*sp.diff(p,z)-(3+delta)*z*p)/L
    require(all(sp.cancel(c-sp.diff(exact,z,k))==0
                for k,c in enumerate(_coefficient(z,delta))),
            'wrong common pressure mass derivative')
    return True


def evaluate(*,inputs,pins):
    old,nominal,q2,core = (inputs[k] for k in ('actual','nominal','q2','core'))
    p = inputs['pressure']['compliant_source']
    family=q2['source_family']
    require(old['passed'] and old['signed_six_same_source_core_and_lower_cover'] and nominal['passed']
            and old['frozen_inputs']['experiments/root_st073/lei_ren_part1_paper_compliant_current_original_n1_correlated_q2_inlet.json']==pins['q2']
            and nominal['frozen_inputs']['experiments/root_st073/lei_ren_part1_paper_compliant_pressure_source.json']==pins['pressure']
            and all(p[k]==core[k]==family[k] for k in ('implicit_source_sha256','datum_enclosure_sha256')),
            'the two pressure decompositions do not have the same actual mass')
    _pair_proof()

    iv,mp = MPIntervalContext(),MPContext()
    iv.prec=1800
    m2 = read(iv,p['fixed_beta2_mass_normalized'])
    delta = read(iv,q2['actual_selected_leading_pressure_primitive']['delta_enclosure'])
    j = read(iv,core['required_j'])
    invp2 = _dyadic_inv_p2(iv,read(iv,core['logPstar']))
    require(0<m2.a and 0<delta.a<=delta.b<iv.mpf(1)/8 and 0<j.a<=j.b<1,
            'invalid pressure or core parameters')
    curves = nominal['remainder_components']
    # The accepted late and nonlinear inlet viscosities are twice the core
    # endpoint operator. Bound the hydro p-slope on a radius-1/4 capsule by
    # 521 under |Z|<=2, |d|<=5, |L|>=1/2, |U_axis|<=9.
    late,nonlinear = (tuple(read(iv,v).b/2 for v in curves[name])
                      for name in ('late_pressure','true_inlet_nonlinear'))
    a,b,c = (521*factorial(k)*8**k*invp2.b for k in range(3))
    hydro0 = b+4*a
    hydro1 = c+2*b+4*(a+b)
    hydro = 16*hydro0,16*(hydro1+hydro0)
    e0,e1,e2 = (read(iv,v).b for v in old['own_mean_rho_minus_half_Psi_rho_ordinary_C0_C1_C2_radius'])
    mean = 32*(e0+e1),32*(2*e0+4*e1+e2)
    packed = {name:tuple(read(iv,v).b for v in pair)
              for name,pair in old['core_pressure_volume_radii'].items()}
    radii=dict(late=late,nonlinear=nonlinear,hydro=hydro,own_mean=mean,**packed)
    budget=tuple(sum(pair[k] for pair in radii.values()).b for k in range(2))

    def encode(v):
        return dict(lower_exact_mpf_tuple=list(v._mpi_[0]),upper_exact_mpf_tuple=list(v._mpi_[1]),
                    rounded=[mp.nstr(mp.make_mpf(bits),17) for bits in v._mpi_],sign=sign(v))

    def enclose(row,Z):
        result=[]
        new_rows=[]
        coefficients = _coefficient(Z,delta)
        for k,(coef,radius) in enumerate(zip(coefficients,budget,strict=True)):
            paired=m2*coef+read(iv,row['axis_inertia'][k])+hull(iv,-radius,radius)
            other=read(iv,row['pressure'][k])
            lo,hi=max(paired.a,other.a),min(paired.b,other.b)
            require(lo<=hi,'two enclosures of the same actual pressure are disjoint')
            result.append(hull(iv,lo,hi))
            new_rows.append(paired)
        return tuple(result),tuple(new_rows)

    cells=[]
    for row in old['numeric_cells']:
        Z = read(iv,row['axial_cell'])
        pressure,independent = enclose(row,Z)
        sums=[]
        for k,value in enumerate(pressure):
            others = sum((read(iv,c[k]) for i,c in enumerate(row['actual_q2_six_in_common_pressure_unit']) if i!=1),iv.mpf(0))
            sums.append(encode(value+others))
        cells.append(dict(axial_cell=row['axial_cell'],pressure_C0_C1=[encode(v) for v in pressure],
                          independent_mass_pair_C0_C1=[encode(v) for v in independent],
                          actual_q2_total_common_pressure_unit=sums))

    points={}
    for z,parent in old['archived_first_two_row_comparisons'].items():
        bound,independent = enclose(old['point_queries'][z],iv.mpf(z))
        legacy=parent['core_rows'][1]
        require(parent['ledger_count']==212 and legacy['unit']==UNIT,'wrong core pressure archive locator')
        records=[]
        for k,(v,new,item) in enumerate(zip(bound,independent,legacy['raw_undivided_row']['coefficient_C0_C1'],strict=True)):
            archive=iv.make_mpf(tuple(tuple(bits) for bits in item['source_enclosure_exact_mpf_tuples']))
            finite=iv.mpf(item['approximation'])
            records.append(dict(order=k,actual_minus_stored_finite=encode(v-finite),
                archived_contains_independent_mass_pair=bool(archive.a<=new.a and new.b<=archive.b),
                archived_contains_refined=bool(archive.a<=v.a and v.b<=archive.b),
                earlier_archive_contains_unpaired=legacy['ordinary_Z'][k]['archived_enclosure_contains_actual']))
        points[z]=dict(ledger_sha256=parent['ledger_sha256'],ordinary_Z=records)

    return dict(passed=True,unit=UNIT,exact_mass_coefficient_and_derivative=True,
                one_shared_actual_M2_used_per_ordinary_order=True,whole_axis_cells=cells,
                pressure_error_radii={k:[encode(v.b) for v in row] for k,row in radii.items()},
                legacy_pressure_diagnostics=points, prior_intersection_is_of_proven_same_function=True,
                inherited_error_packets_removed=False, legacy_binary_roundoff_bound=False,
                independent_terminal_matching_closed=False, ancestor_calculations_replayed=False)
