"""Independent Fourier/implicit-derivative checks and source consumption."""
import copy
import gzip
import json
import math
from pathlib import Path
import time

import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_generic_shear_loop as loop
import lei_ren_part1_paper_compliant_current_original_n1_preRp_pressure_C1 as current


def symbolic_derivatives():
    z = sy.Symbol('z', real=True)
    x = sy.Symbol('x')
    f = {(i,j):sy.Symbol('f%d%d'%(i,j)) for i in range(4) for j in range(4-i)}
    p1, p2, p3 = sy.symbols('p1 p2 p3')
    p = p1*z+p2*z*z/2+p3*z**3/6
    # An arbitrary multivariate Taylor jet proves the general C3 chain rule
    # without relying on SymPy's names for Subs/Derivative dummy variables.
    F = sum(value*x**i*z**j/(sy.factorial(i)*sy.factorial(j)) for (i,j),value in f.items())
    actual = sy.diff(F.subs(x,p),z,3).subs(z,0)
    expected = (f[0,3]+3*f[1,2]*p1+3*f[2,1]*p1*p1+f[3,0]*p1**3
                +3*f[1,1]*p2+3*f[2,0]*p1*p2+f[1,0]*p3)
    if sy.simplify(actual-expected) != 0:
        raise AssertionError('Independent total-third implicit/composition chain rule failed')
    x, y = sy.symbols('x y', real=True)
    cov = 1/sy.cosh(x-y)
    norms = [sy.simplify(sy.diff(cov, x, k, y, k).subs({x:0, y:0})) for k in range(4)]
    if norms != [1, 1, 5, 61]:
        raise AssertionError('Parseval hyperbolic covariance derivative norms failed')
    u = sy.Symbol('u', real=True)
    h = 1/sy.sqrt(1+u*u)
    r = u*h
    if sy.simplify(sy.diff(r, u)-h**3) != 0 or sy.simplify(sy.diff(h, u)+r*h*h) != 0:
        raise AssertionError('Genuine u derivatives of Fourier coefficients failed')
    if not (float(sy.sqrt(5)+1) < 4 and float(sy.sqrt(61)+3*sy.sqrt(5)+2) < 17):
        raise AssertionError('Third u-chain L2 majorants failed')
    return dict(implicit_and_composed_total_Z3=True, covariance_norms=[1,1,5,61],
        u_chain_L2_constants=[1,1,4,17],
        point_bound_proof=dict(chi_u='2*(2+4)=12', chi_uu='8h+40h+32h<=80',
            chi_psiu='8/h^2+32/h^2=40/h^2', chi_psipsi='32/h^5',
            geometric_source='1/(1-abs(r))=(1+abs(r))/h^2<=2/h^2'),
        independent_symbolic_differentiation=True)


def _fixture_B(z, base_u, phase):
    """Synthetic nonunit input, evaluated by the original analytic loop."""
    obj = object.__new__(loop.GenericShearLoop)
    obj.ctx = mp.mp
    obj.a = mp.mpf('1.3')+mp.mpf('.07')*z+mp.mpf('.02')*z*z
    obj.t0 = mp.mpf('.2')+mp.mpf('.3')*z-mp.mpf('.04')*z*z+mp.mpf('.01')*z**3
    obj.q = mp.mpf('.4')+mp.mpf('.08')*z+mp.mpf('.03')*z*z-mp.mpf('.01')*z**3
    obj.u = mp.mpf(base_u)+mp.mpf('.2')*z-mp.mpf('.1')*z*z+mp.mpf('.03')*z**3
    obj.h = mp.sqrt(1+obj.u**2)  # The loop engine names sqrt(1+u^2) h.
    obj.r = obj.u/obj.h
    obj.one_minus_abs_r = 1/(obj.h*(obj.h+abs(obj.u)))
    obj.one_minus_r2 = 1/obj.h**2
    obj.b = -obj.a*obj.t0
    obj.Utheta = mp.mpf('1.7')+mp.mpf('.11')*z+mp.mpf('.02')*z*z
    nu = 1+obj.t0**2+2*obj.q**2
    def equation(psi):
        return psi+obj.integrals(psi)[1]-2*mp.pi*phase*nu
    psi = mp.findroot(equation, (2*mp.pi*phase*mp.mpf('.9'), 2*mp.pi*phase*mp.mpf('1.1')),
                      tol=mp.eps*16)
    if not 0 < psi < 2*mp.pi:
        raise AssertionError('Fixture inverse left its phase period')
    B = obj.Utheta/2*(-obj.a*obj.integrals(psi)[0]/(2*mp.pi)-obj.b*phase)
    return B


def finite_loop_fixtures():
    results = []
    with mp.workdps(65):
        c = MPIntervalContext()
        c.dps = 100
        bd = current.rg.Bounds(c)
        for base_u in ('-3', '0', '3'):
            a = lambda z: mp.mpf('1.3')+mp.mpf('.07')*z+mp.mpf('.02')*z*z
            t0 = lambda z: mp.mpf('.2')+mp.mpf('.3')*z-mp.mpf('.04')*z*z+mp.mpf('.01')*z**3
            q = lambda z: mp.mpf('.4')+mp.mpf('.08')*z+mp.mpf('.03')*z*z-mp.mpf('.01')*z**3
            u = lambda z: mp.mpf(base_u)+mp.mpf('.2')*z-mp.mpf('.1')*z*z+mp.mpf('.03')*z**3
            E = lambda z: mp.mpf('1.7')+mp.mpf('.11')*z+mp.mpf('.02')*z*z
            functions = dict(a=a, b=lambda z:-a(z)*t0(z), E=E, t0=t0,
                             p2=lambda z:mp.mpf('.7')*u(z)/q(z))
            def records(fn):
                return [current.record(bd.constant(c.mpf(mp.nstr(abs(mp.diff(fn, 0, k)), 80)))) for k in range(4)]
            cell = dict(A=[current.record(bd.constant(c.exp(1000))) for _ in range(4)],
                B=[current.record(bd.constant(c.exp(1000))) for _ in range(4)], proof=dict(
                full_original_phase_box=[0,1], live_source_phase_Z_exact_zero=True, exact_same_source_P0=True,
                branch='synthetic_active_fixture', actual_source_ordinary_Z3_rows={k:records(v) for k,v in functions.items()},
                q_bounds=records(q), actual_dstar_log=current.record(c.log(c.mpf('.7'))),
                actual_C3_function_handles={'synthetic_fixture_only':True}))
            bound = current.sharpen_primitive(c, cell)['candidate_B'][3]
            for phase in (mp.mpf('.21'), mp.mpf('.77')):
                derivative = mp.diff(lambda z:_fixture_B(z, base_u, phase), 0, 3)
                if mp.log(abs(derivative)) > current.ep(bound.log)[0]:
                    raise AssertionError('Original analytic B third derivative exceeds new bound')
                results.append(dict(u0=base_u, phase=mp.nstr(phase, 8),
                    third_derivative=mp.nstr(derivative, 30), synthetic_fixture_only=True))
    return results


def finite_pressure_fixture():
    """Differentiate actual reduced integrals with a nonzero finite delta."""
    with mp.workdps(70):
        c = MPIntervalContext()
        c.dps = 100
        bd = current.rg.Bounds(c)
        delta, Rc, rmin, center = map(mp.mpf,('.2','3.7','.23','.37'))
        U = lambda z:mp.mpf('.4')+mp.mpf('.13')*z+mp.mpf('.07')*z*z+mp.mpf('.03')*z**3
        q = lambda z:((1+delta)*z*U(z)-(1-z*z)*mp.diff(U,z))/(1-delta*z*z)
        def reduced(z):
            L,d,a = 1-delta*z*z,1-z*z,1-delta
            u,uz,qq,qz = U(z),mp.diff(U,z),q(z),mp.diff(q,z)
            w0,w1 = Rc-rmin,(Rc*Rc-rmin*rmin)/2
            s0 = w0*(qq/L+a*z*qz/(2*L)+qq*qq/2+d*u*qz/L
                -2*z*z*u*u/(L*L)-2*delta*z*z*u*u/(L*L)+2*z*d*u*uz/(L*L))
            s1 = w1*(2*qq/L-qq/L+a*z*qz/(2*L)+qq*qq-qq*qq/2+d*u*qz/L
                -4*z*z*u*u/(L*L)+2*a*z*z*u*u/(L*L)+2*z*d*u*uz/(L*L))
            return -s0/2,s1/2
        rows = lambda fn,n:[bd.constant(c.mpf(mp.nstr(abs(mp.diff(fn,center,k)),90))) for k in range(n)]
        view = dict(U=rows(U,4),q=rows(q,3))
        endpoint_caps = dict(q_plus_q_y_Z=rows(q,2)[1],q_y_Z=bd.zero)
        end = dict(lower=dict(caps=endpoint_caps),upper=dict(caps=endpoint_caps))
        bounds = current.active_pressure_C1(c,('.37','.37'),[view],end,c.mpf('.2'))
        actual = [mp.diff(lambda z:reduced(z)[i],center) for i in range(2)]
        for key,power,value in zip(('m13','m15'),(1,2),actual):
            if abs(value)/Rc**power > c.exp(bounds[key].log).a:
                raise AssertionError('Actual finite reduced pressure derivative exceeds C1 bound')
        return dict(Rc='3.7',rmin='.23',delta='.2',Z='.37',
            differentiated_actual_integrals=True,nonconstant_U_q_and_coefficients=True,
            actual_m13_Z=mp.nstr(actual[0],30),actual_m15_Z=mp.nstr(actual[1],30))


@current.source_precision
def run(owner=None):
    began = time.monotonic()
    owner = owner or current.CurrentOriginalN1PreRpPressureC1(False)
    raw = json.loads(gzip.decompress((current.HERE/current.NAME).read_bytes()))
    if raw['source_family'] != owner.family or any(raw.get(k) for k in current.OPEN):
        raise AssertionError('Original source family or open gate changed')
    for name, digest in raw['input_hashes'].items():
        if current.sha(name) != digest:
            raise AssertionError('Changed original C1 source '+name)
    symbolic = symbolic_derivatives()
    fixtures = finite_loop_fixtures()
    pressure_fixture = finite_pressure_fixture()
    rebuilt = current.record(current.build_data(owner.c, **owner.parents))
    if rebuilt != raw['source_data']:
        raise AssertionError('Actual source transport/pressure C1 archive differs')
    entries = rebuilt['native_primitive_bounds']
    charts = rebuilt['upgraded_chart_views']
    if len(entries) != 228 or len(charts) != 68:
        raise AssertionError('Original phase/radial/axial coverage missing')
    for row in entries:
        if not (row['predecessor_never_reset'] and row['rate_one_kernel_upper']==1
                and row['actual_whole_chart_rate_zero_width_used_as_subcell_upper']):
            raise AssertionError('Original rate-one history dropped')
        v = current.logcap(owner.c, row['new_physical_local_B_ZZZ_over_N'])
        if v.log is not None and current.ep(v.log)[1] >= -1000:
            raise AssertionError('Actual same-N physical third increment not controlled')
    grid = owner.upper_C1_grid()
    maximum = {}
    if len(grid) != 4096:
        raise AssertionError('Complete accepted derivative grid missing')
    for key in ('m13','m15'):
        maximum[key] = max(current.ep(v['fields'][key]['total_error'])[1] for v in grid)
        if maximum[key] >= current.ep(owner.c.mpf('.001'))[0]:
            raise AssertionError('Actual rminus-to-infinity derivative accuracy failed')
        parent = owner.pulse_data['c1_grid']
        for new, old in zip(grid, parent):
            if new['fields'][key]['center_slope'] != old['fields'][key]['center_slope']:
                raise AssertionError('A cap replaced the signed derivative source')
            added = new['fields'][key]['preRp_nonzero_source_error']
            if current.ep(added)[0] <= 0:
                raise AssertionError('A bounded pre-Rp term was assigned zero')
    before = mp.mp.dps
    with mp.workdps(15):
        view = owner.derivative_enclosure('.371')
        if mp.mp.dps != 15 or view['radial_scope'] != 'original r_minus to infinity':
            raise AssertionError('Public default-precision derivative query failed')
    if mp.mp.dps != before:
        raise AssertionError('Public query changed caller precision')
    pins = dict(owner.hashes)
    for name in (current.NAME, Path(current.__file__).name, Path(__file__).name):
        pins[name] = current.sha(name)
    receipt = dict(all_passed=True, source_family=owner.family,
        checks=dict(symbolic=symbolic, original_loop_finite_third_fixtures=fixtures,
            independent_nonunit_pressure_C1=pressure_fixture,
            actual_source_replay_equal=True, native_cells=228, complete_charts=68,
            no_history_reset=True, common_Pstar_and_selected_N_once=True,
            original_q_ZZ_and_pressure_C1_consumption=True, derivative_grid_count=4096,
            unchanged_signed_derivative_sources=True, nonzero_preRp_allowance=True,
            default_precision_restored=True), maximum_normalized_derivative_error=maximum,
        input_hashes=pins, execution_seconds=time.monotonic()-began,
        **{current.GATE:True,current.UPPER_GATE:True}, **dict.fromkeys(current.OPEN,False))
    (current.HERE/current.RECEIPT).write_bytes((json.dumps(current.record(receipt),separators=(',',':'))+'\n').encode())
    print('PASS_PRE_RP_PRESSURE_C1',len(charts),len(entries),len(grid),flush=True)
    return receipt


if __name__ == '__main__':
    run()
