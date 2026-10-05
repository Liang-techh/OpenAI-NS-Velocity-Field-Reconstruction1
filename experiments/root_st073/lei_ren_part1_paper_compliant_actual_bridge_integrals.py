"""Finite-width original prescribed-shear bridge and its own moment histories.

The known direction is the comparison direction, as in (4.34)--(4.36).
Actual moments are integrated from actual F/V; they never replace that known
direction. Signed source integrals retain hb, Pstar^2 and F0^2 in exact log
factors. The nonlinear angular prefix has a separate hb^2 error bound.
No cap, interval endpoint or parameter midpoint defines a production field.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp

from lei_ren_part1_paper_compliant_comparison_point_integrals import (
    CompliantComparisonPointIntegrals, WidthPolynomial, jet_norm, upper,
)
import lei_ren_part1_paper_compliant_comparison_point_integrals as comparison_module
import lei_ren_part1_paper_compliant_core_integral_atoms as atom_module
from lei_ren_part1_paper_compliant_core_integral_atoms import CompliantCoreIntegralAtoms
from lei_ren_part1_paper_compliant_core_physical_field import intersection
from lei_ren_part1_paper_compliant_inner_bridge_profiles import direction, symmetric
from lei_ren_part1_paper_compliant_macro_signed_integrals import _mode_rows
from lei_ren_part1_paper_compliant_frozen_comparison_field import (
    MTH, MZ, MTHZ, MZT, MP, derivative, dress, square,
)
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_candidate_exact_amplitude import _encode

HERE = Path(__file__).parent
PREFIX = 'lei_ren_part1_paper_compliant_'
RATES = dict(H=2, M=1, K=2, A=1, B=2, C=1)


def real_chi_model_tail(c,chi,*args):
    """Keep the exact real H^2/(H^2+sigma^2) correlation before tail bounds.

    Naive division loses the shared H^2 and can give an upper value >1 on
    real Z intervals. The original ratio is in[0,1], since sigma^2>0.
    Only coefficient0 is intersected; every derivative row is unchanged.
    """
    correlated=IntervalTaylor(c,[intersection(c,chi[0],c.mpf([0,1]))]
                              +list(chi.coefficients[1:]))
    return atom_module.bessel_radial_tail_coefficients(c,correlated,*args)


def real_model_replay(module,method):
    """Replay the unchanged admitted method with a local correlated bound.

    No module/global monkeypatch or source-file mutation. All original
    arithmetic is compiled unchanged; only its final model-tail callback
    supplies the independently exact real-ratio intersection above.
    """
    tree=ast.parse(Path(module.__file__).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name==method)
    env=dict(module.__dict__)
    env['bessel_radial_tail_coefficients']=real_chi_model_tail
    exec(compile(ast.Module(body=[fn],type_ignores=[]),module.__file__,'exec'),env)
    return env[method]


class BridgeCoreAtoms(CompliantCoreIntegralAtoms):
    model_tail_coefficients=real_model_replay(atom_module,'model_tail_coefficients')


class BridgeComparison(CompliantComparisonPointIntegrals):
    source_profile_jet=real_model_replay(comparison_module,'source_profile_jet')

    def __init__(self):
        super().__init__()
        self.atoms=BridgeCoreAtoms();self.core=self.atoms.core;self.ctx=self.atoms.ctx
        # Whole-axis axial interval rows are much larger than point rows.
        # A stronger bound on the same exact hb controls their exponential
        # remainders; no source parameter is changed or selected.
        self.weight=self.atoms.axial_weight;self.cap=self.ctx.mpf('1e-50000')
        if endpoints(self.bridge.logh)[1]>=endpoints(self.ctx.ln(self.cap))[0]:
            raise ValueError('Exact source width does not admit whole-axis remainder cap')
        self.hashes.update(self.atoms.hashes)
        if self.bridge.source!=self.core.source or self.bridge.family!=self.core.family:
            raise ValueError('Correlated real source must preserve the defining family')
        self.cache={};self.weight_cache={}


def sha(name):
    return hashlib.sha256((HERE / name).read_bytes()).hexdigest()


def pulse_control_masses(c,phase,cells=64):
    """Directed first-chart masses, including non-dyadic decimal phases.

    Integrate to the upper phase endpoint, then use 0<=1-sigma<=1 to
    enclose every requested phase. Quadrature nodes are not field parameters.
    """
    s=c.mpf(phase);lo,hi=endpoints(s)
    if lo<0 or hi>1:raise ValueError('First-chart phase in[0,1] required')
    if hi==0:return c.mpf(0),c.mpf(0)
    if lo==hi==1:return c.mpf('.5'),c.mpf('.5')
    step=c.mpf(hi)/cells;total=c.mpf(0);error=c.mpf(0)
    for n in range(cells):
        a=step*n;b=step*(n+1);mid=(a+b)/2
        values=[1-sigma_jets(c,t)[0] for t in (a,mid,b)]
        cover=sigma_jets(c,c.mpf([endpoints(a)[0],min(mp.mpf(1),endpoints(b)[1])]))
        total+=step*(values[0]+4*values[1]+values[2])/6
        error+=step**5*upper(c,abs(cover[4])*math.factorial(4))/2880
    complement=total+symmetric(c,error)-c.mpf([0,hi-lo])
    positive=c.mpf([0,hi])
    complement=intersection(c,complement,positive)
    sigma=intersection(c,s-complement,positive)
    return complement,sigma


def named_moments(values):
    return {MTH: values['H'], MZ: values['M'], MTHZ: values['K'],
            MZT: dict(axial=values['A'], swirl=values['B']), MP: values['C']}


def value_enclosure(poly):
    """Interval extension of the admitted hb^2 + controlled hb^3 source.

    The interval [0,cap] bounds the same exact hb in every polynomial. It is
    not a selected width. The retained rows are also exported separately.
    """
    c = poly.ctx
    hbox = c.mpf([0, endpoints(poly.cap)[1]])
    result = poly.rows[0] + poly.rows[1]*hbox + poly.rows[2]*hbox**2
    return result + IntervalTaylor(c, [
        symmetric(c, poly.remainder*hbox**3/poly.weight**k)
        for k in range(result.order+1)])


def radius_kernels(c, radius0, radius1, length, power):
    """Exact int R^power exp(-j Delta)dDelta, j=0,1,2.

    R1=R0 exp(length) is simplified before arithmetic. This avoids an
    enormous intermediate exp(length), and keeps the fixed R100 endpoint.
    """
    if power == 1:
        return [radius1-radius0, radius0*length,
                radius0-radius0*radius0/radius1]
    if power == 2:
        return [(radius1*radius1-radius0*radius0)/2,
                radius0*radius1-radius0*radius0, radius0*radius0*length]
    raise ValueError('Original bridge uses only radius powers 1 and 2')


def weighted_modes(rows, kernels):
    return sum((row*kernel for row, kernel in zip(rows, kernels)), rows[0]*0)


def weighted_mode_norm(rows, kernels, weight):
    c = rows[0].ctx
    return sum((jet_norm(row, weight)*upper(c, abs(kernel))
                for row, kernel in zip(rows, kernels)), c.mpf(0))


def volterra_primitive_kernels(c,radius0,radius1,length,power,rate):
    """int exp(-rate*(L-t))*int_0^t R^p exp(-j*s)ds dt.

    The inner primitive, including its resonant t factor, is integrated
    before interval enclosure. R0^p exp((p-j)L)=R0^j R1^(p-j).
    """
    decay=c.exp(-rate*length);weight=(1-decay)/rate;rows=[]
    for j in range(3):
        a=power-j
        if a==0:
            value=radius0**power*(length/rate-weight/rate)
        else:
            integral=(radius0**power*length*decay if a+rate==0 else
                      (radius0**j*radius1**a-radius0**power*decay)/(a+rate))
            value=(integral-radius0**power*weight)/a
        rows.append(value)
    return rows


def own_moment_feedback(phi0, V0, initial, theta, delta_phi, delta_V):
    """Own six Volterra histories, with no comparison-history substitution.

    Delta fields enclose every actual prefix on [0,y]. Positive Volterra
    kernels preserve these signed coefficient bounds. Mixed and quadratic
    changes are formed before enclosure; initial atom uncertainty cancels
    from the feedback rather than being subtracted twice.
    """
    targets = dict(H=phi0, M=V0, K=phi0*V0, A=square(V0),
                   B=square(phi0)/2, C=square(phi0))
    changes = dict(H=2*delta_phi, M=delta_V,
        K=2*(phi0*delta_V+V0*delta_phi+delta_phi*delta_V),
        A=2*V0*delta_V+square(delta_V),
        B=2*phi0*delta_phi+square(delta_phi),
        C=2*phi0*delta_phi+square(delta_phi))
    baseline = {}; feedback = {}; actual = {}
    for name, rate in RATES.items():
        decay = theta**rate
        baseline[name] = initial[name]*decay+targets[name]*(1-decay)
        feedback[name] = changes[name]*((1-decay)/rate)
        actual[name] = baseline[name]+feedback[name]
    return dict(actual=actual, baseline=baseline, feedback=feedback,
                source_changes=changes)


def bridge_control_bindings():
    """Bind the known direction and distinct original microscopic controls."""
    path = HERE / (PREFIX+'bridge_mixed_C4.py')
    tree = ast.parse(path.read_text(encoding='utf8'))
    fn = next(n for n in ast.walk(tree)
              if isinstance(n, ast.FunctionDef) and n.name == 'evaluate')
    wanted = {
        'first_chi': ('chi', "[algebra.lift(1-sigma[0])+h*sigma[0]]+[-algebra.lift(sigma[k])+h*sigma[k] for k in range(1,4)]"),
        'second_and_macro_chi': ('chi', '[h]+[h*0]*3'),
        'comparison_direction': ('(directions, momentrows)', "comparison_directions(algebra,Z,self.core.delta,barphi,barV,comparison['moments'],inp['p0'],inp['F0_ratios'],inp['F0_squared_ratios'],scale)"),
    }
    result = {}
    for name, (target, expression) in wanted.items():
        expr = ast.dump(ast.parse(expression, mode='eval').body)
        assignments = [n for n in ast.walk(fn) if isinstance(n, ast.Assign)
                       and any(ast.unparse(t) == target for t in n.targets)]
        count = sum(ast.dump(n.value) == expr for n in assignments)
        expected = 2 if name == 'second_and_macro_chi' else 1
        if count != expected:
            raise ValueError('Original bridge control changed: '+name)
        result[name] = dict(verified=True, matching_assignments=count,
                            expression=expression)
    actual=ast.parse((HERE/(PREFIX+'inner_bridge_profiles.py')).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(actual) if isinstance(n,ast.FunctionDef) and n.name=='actual')
    source=next(kw.value for n in ast.walk(fn) if isinstance(n,ast.Call)
                for kw in n.keywords if kw.arg=='source_integral_definitions')
    expected=dict(F='F=f*exp(-.5*integral_0^y chi*Dbar dt)',
        V='V=v-integral_0^y chi*(phi_actual/phi_bar)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl) dt')
    for name,text in expected.items():
        entries=[kw.value for kw in source.keywords if kw.arg==name]
        if len(entries)!=1 or not isinstance(entries[0],ast.Constant) or entries[0].value!=text:
            raise ValueError('Original actual source integral changed: '+name)
        result['actual_'+name+'_source_integral']=dict(verified=True,expression=text)
    return result


class CompliantActualBridgeIntegrals:
    def __init__(self):
        self.comparison = BridgeComparison()
        self.bridge = self.comparison.bridge.bridge
        self.core = self.comparison.core
        self.ctx = c = self.comparison.ctx
        self.direction_context = SimpleNamespace(ctx=c,delta=self.core.delta)
        self.weight = self.comparison.weight
        self.cap = self.comparison.cap
        self.logh = self.bridge.logh
        self.Y = c.ln(100/self.bridge.r)
        self.pressure_log = 2*self.core.logP
        # Gbar is a bound, not G(Z). Preserve the full exact F0(Z)^2 range.
        self.swirl_log = c.mpf([endpoints(-2*self.core.logC-2*self.core.Lambda*self.core.Gbar)[0],
                               endpoints(-2*self.core.logC)[1]])
        self.hashes = dict(self.comparison.hashes)
        for name in (PREFIX+'comparison_point_integrals.json',
                     PREFIX+'comparison_point_integrals_check.json',
                     PREFIX+'macro_signed_integrals.py', Path(__file__).name):
            self.hashes[name] = sha(name)
        check = json.loads((HERE/(PREFIX+'comparison_point_integrals_check.json')).read_bytes())
        if not check.get('all_passed'):
            raise ValueError('Accepted comparison finite-width remainder required')
        for name, digest in check['input_hashes'].items():
            if sha(name) != digest:
                raise ValueError('Comparison source changed: '+name)
        self.bindings = bridge_control_bindings()
        self.proofs = []
        self.cache = {}
        self.pulse_cache = {}

    def scaled_jet_box(self, row, power, extra_log=0, norm_only=False):
        """Log proof before enclosing an unmaterializable positive factor."""
        c = self.ctx
        bound = jet_norm(row, self.weight)
        if endpoints(bound)[1] == 0:
            return row*0
        log_bound = power*self.logh+c.mpf(extra_log)+c.ln(bound)
        if endpoints(log_bound)[1] > endpoints(c.ln(self.cap))[0]:
            raise ArithmeticError('Source-weighted bridge term exceeds admitted cap')
        self.proofs.append(dict(width_power=power, extra_source_log=c.mpf(extra_log),
            normalized_weighted_jet_norm_upper=bound, product_log_upper=log_bound,
            enclosure_norm_cap=self.cap, passed=True))
        values = []
        for k, value in enumerate(row.coefficients):
            cap = endpoints(self.cap/self.weight**k)[1]
            lo, hi = endpoints(value)
            values.append(c.mpf([-cap if lo < 0 or norm_only else 0,
                                cap if hi > 0 or norm_only else 0]))
        return IntervalTaylor(c, values)

    def symmetric_norm_row(self, bound, order=5):
        c = self.ctx
        return IntervalTaylor(c, [symmetric(c, bound/self.weight**k)
                                 for k in range(order+1)])

    def comparison_micro_cover(self, data, a, b):
        """Uniform finite-width comparison source on the entire phase cell."""
        c = self.ctx
        s = c.mpf([a, b]); smax = c.mpf(b)
        if b <= 1:
            A0 = s; A1 = s*s/2; J = s*s/2
        else:
            # alpha, u*alpha and int_0^s A0 are nonnegative. These bounds
            # cover the full first and second charts, including their join.
            A0 = c.mpf([a if a <= 1 else 1, b])
            A1 = c.mpf([a*a/2 if a <= 1 else '.5', b*b/2])
            J = s*A0-A1
        phi, V = self.comparison.fields(data, s, A0, A1)
        p, v, L1, V1 = data['phi0'], data['V0'], data['L1'], data['V1']
        rhs = dict(H=2*phi, M=V, K=2*phi*V, A=V*V, B=phi*phi, C=phi*phi)
        r0 = dict(H=2*p, M=v, K=2*p*v, A=v*v, B=p*p, C=p*p)
        r1 = dict(H=2*p*L1, M=V1, K=2*(p*L1*v+p*V1),
                  A=2*v*V1, B=2*p*p*L1, C=2*p*p*L1)
        moments = {}
        for name, initial in data['moments'].items():
            rate = RATES[name]
            first = r0[name]-rate*initial
            second = r1[name]*J-rate*first*s*s/2
            bound = (jet_norm(initial,self.weight)*(rate*smax)**3/6
                +rate**2*smax**3*jet_norm(r0[name],self.weight)/6
                +rate*smax*smax*jet_norm(rhs[name].rows[1],self.weight)/2
                +smax*(jet_norm(rhs[name].rows[2],self.weight)
                       +self.cap*rhs[name].remainder))*c.exp(rate*self.cap*smax)
            moments[name] = WidthPolynomial([initial,first*s,second],
                                           bound,self.weight,self.cap)
        return dict(phi=phi, V=V, moments=moments)

    def prepare(self, Z, root=False):
        c = self.ctx
        key = ('root' if root else c.mpf(Z)._mpi_)
        if key in self.cache:
            return self.cache[key]
        data = self.comparison.inlet(Z, root)
        z = data['Z']
        inputs = {name:IntervalTaylor(c,value.coefficients) if isinstance(value,IntervalTaylor) else value
                  for name,value in self.bridge.inputs(z).items()}
        covers = [self.comparison_micro_cover(data,0,1),
                  self.comparison_micro_cover(data,1,2)]
        micros = []
        for cover in covers:
            phi = value_enclosure(cover['phi']); V = value_enclosure(cover['V'])
            moments = named_moments({n:value_enclosure(p) for n,p in cover['moments'].items()})
            dirs = direction(c,z,self.core.delta,phi,V,moments,inputs['p0'],
                             inputs['F0_ratios'],inputs['F0_squared_ratios'])
            micros.append(dict(direction=dirs, quotient=data['phi0'].truncate(5)/phi.truncate(5),
                               comparison=cover))
        end = self.comparison.micro(Z,2,root)
        phi2 = value_enclosure(end['phi']); V2 = value_enclosure(end['V'])
        moment2 = named_moments({n:value_enclosure(p) for n,p in end['moments'].items()})
        modes = _mode_rows(self.direction_context,z,inputs,phi2,V2,moment2)
        R0 = self.bridge.r*c.exp(c.mpf([0,endpoints(2*self.cap)[1]]))
        L = self.Y-c.mpf([0,endpoints(2*self.cap)[1]])
        kernels = radius_kernels(c,R0,c.mpf(100),L,1)
        # Uniform in every actual prefix, not only the final signed integral.
        B = (jet_norm(micros[0]['direction']['D_over_R'],self.weight)
             *upper(c,R0)*(1+self.cap)
             +jet_norm(micros[1]['direction']['D_over_R'],self.weight)
             *upper(c,R0)*self.cap
             +weighted_mode_norm(modes['D_over_R'],kernels,self.weight))/2
        ell_box = self.scaled_jet_box(self.symmetric_norm_row(B),1,norm_only=True)
        # Dbar>0 is an inherited same-family theorem, independent of boxes.
        ell_box = IntervalTaylor(c,[c.mpf([-endpoints(abs(ell_box[0]))[1],0])]
                                +list(ell_box.coefficients[1:]))
        prepared = dict(data=data,Z=z,inputs=inputs,micros=micros,end=end,
            modes=modes,R0=R0,L=L,quotient=data['phi0'].truncate(5)/phi2.truncate(5),
            angular_prefix_norm_per_hb=B,uniform_ell_box=ell_box)
        uniform_delta_V = data['phi0'].truncate(5)*0
        for name, powerR, scale in (('hydro',1,c.mpf(0)),
                                    ('pressure',1,self.pressure_log),
                                    ('swirl',2,self.swirl_log)):
            # Integral of absolute modal norms bounds every prefix. A signed
            # endpoint integral cannot bound earlier prefixes after cancellation.
            n = (jet_norm(micros[0]['quotient']*micros[0]['direction']['drive_'+name],self.weight)
                 *upper(c,R0**powerR)*(1+self.cap)
                 +jet_norm(micros[1]['quotient']*micros[1]['direction']['drive_'+name],self.weight)
                 *upper(c,R0**powerR)*self.cap
                 +jet_norm(prepared['quotient'],self.weight)
                 *weighted_mode_norm(modes['drive_'+name],
                     radius_kernels(c,R0,c.mpf(100),L,powerR),self.weight))
            uniform_delta_V += self.scaled_jet_box(self.symmetric_norm_row(n*c.exp(6*self.cap)),
                                                  1,scale,norm_only=True)
        prepared['uniform_delta_V'] = uniform_delta_V
        self.cache[key] = prepared
        return prepared

    def pulse_masses(self, phase):
        c = self.ctx
        lo, hi = endpoints(c.mpf(phase))
        if lo < 0 or hi > 1:
            raise ValueError('First-chart phase in[0,1] required')
        if hi == 0:
            return c.mpf(0),c.mpf(0)
        if lo == hi == 1:
            return c.mpf('.5'),c.mpf('.5')
        key = c.mpf(phase)._mpi_
        if key not in self.pulse_cache:
            self.pulse_cache[key] = pulse_control_masses(c,phase)
        return self.pulse_cache[key]

    def source_term(self, row, power, extra_log=0, error_norm=0):
        extra=self.ctx.mpf(extra_log)
        source=('Pstar^2' if extra._mpi_==self.pressure_log._mpi_ else
                'F0(Z)^2=Cstar^-2*exp(-2Lambda*G(Z))' if extra._mpi_==self.swirl_log._mpi_ else '1')
        return dict(normalized_axial_coefficients=list(row.coefficients),
            width_power=power, additional_positive_scale_log_enclosure=self.ctx.mpf(extra_log),
            positive_scale_log_enclosure=power*self.logh+self.ctx.mpf(extra_log),
            additional_positive_scale_source=source,
            nonlinear_prefix_error_weighted_norm_per_next_hb=self.ctx.mpf(error_norm),
            nonlinear_prefix_error_width_power=power+1,
            source_width='hb=cstar*K^-100',width_not_materialized=True)

    def contributions(self, prepared, chart, coordinate):
        c = self.ctx; p = prepared; zero = p['data']['phi0'].truncate(5)*0
        angular = []; axial = {n:[] for n in ('hydro','pressure','swirl')}
        s = c.mpf(coordinate)
        angular_norm = p['angular_prefix_norm_per_hb']
        exp_bound = c.exp(6*self.cap)
        scales = dict(hydro=c.mpf(0),pressure=self.pressure_log,swirl=self.swirl_log)
        powerR = dict(hydro=1,pressure=1,swirl=2)
        first_phase = s if chart == 'first' else c.mpf(1)
        wc, ws = self.pulse_masses(first_phase)
        Rmicro = self.bridge.r*c.exp(c.mpf([0,endpoints(2*self.cap)[1]]))
        for width_power, mass in ((1,wc),(2,ws)):
            d = p['micros'][0]['direction']
            angular.append(self.source_term(d['D_over_R']*(-Rmicro*mass/2),width_power))
            for name in axial:
                row = -(p['micros'][0]['quotient']*d['drive_'+name])*Rmicro**powerR[name]*mass
                error = jet_norm(p['micros'][0]['quotient']*d['drive_'+name],self.weight)\
                    *upper(c,Rmicro**powerR[name]*mass)*angular_norm*exp_bound
                axial[name].append(self.source_term(row,width_power,scales[name],error))
        second_phase = s-1 if chart == 'second' else c.mpf(1 if chart == 'macro' else 0)
        d = p['micros'][1]['direction']
        angular.append(self.source_term(d['D_over_R']*(-Rmicro*second_phase/2),2))
        for name in axial:
            row = -(p['micros'][1]['quotient']*d['drive_'+name])*Rmicro**powerR[name]*second_phase
            error = jet_norm(p['micros'][1]['quotient']*d['drive_'+name],self.weight)\
                *upper(c,Rmicro**powerR[name]*second_phase)*angular_norm*exp_bound
            axial[name].append(self.source_term(row,2,scales[name],error))
        if chart == 'macro':
            length = p['L']*s
            R1 = c.mpf(100) if endpoints(s)==(mp.mpf(1),mp.mpf(1)) else self.bridge.r*c.exp(self.Y*s+c.mpf([0,endpoints(2*self.cap*(1-s))[1]]))
            kernels1 = [c.mpf(0)]*3 if endpoints(s)==(mp.mpf(0),mp.mpf(0)) else radius_kernels(c,p['R0'],R1,length,1)
            angular.append(self.source_term(weighted_modes(p['modes']['D_over_R'],kernels1)*(-c.mpf('.5')),1))
            for name in axial:
                kernels = [c.mpf(0)]*3 if endpoints(s)==(mp.mpf(0),mp.mpf(0)) else radius_kernels(c,p['R0'],R1,length,powerR[name])
                row = -p['quotient']*weighted_modes(p['modes']['drive_'+name],kernels)
                error = jet_norm(p['quotient'],self.weight)*angular_norm*exp_bound\
                    *weighted_mode_norm(p['modes']['drive_'+name],kernels,self.weight)
                axial[name].append(self.source_term(row,1,scales[name],error))
        return angular,axial

    def field_changes(self, prepared, angular, axial):
        c = self.ctx
        ell = prepared['data']['phi0'].truncate(5)*0
        for term in angular:
            row = IntervalTaylor(c,term['normalized_axial_coefficients'])
            ell += self.scaled_jet_box(row,term['width_power'])
        ell=IntervalTaylor(c,[c.mpf([endpoints(ell[0])[0],min(mp.mpf(0),endpoints(ell[0])[1])])]
                            +list(ell.coefficients[1:]))
        # Uniform bound controls exp(ell)-1 in the truncated axial jet ring.
        delta_phi_norm = jet_norm(prepared['data']['phi0'].truncate(5),self.weight)\
            *prepared['angular_prefix_norm_per_hb']*c.exp(6*self.cap)
        delta_phi = self.scaled_jet_box(self.symmetric_norm_row(delta_phi_norm),1,norm_only=True)
        delta_phi = IntervalTaylor(c,[c.mpf([-endpoints(abs(delta_phi[0]))[1],0])]
                                  +list(delta_phi.coefficients[1:]))
        delta_V = ell*0
        for terms in axial.values():
            for term in terms:
                row = IntervalTaylor(c,term['normalized_axial_coefficients'])
                delta_V += self.scaled_jet_box(row,term['width_power'],term['additional_positive_scale_log_enclosure'])
                error = self.symmetric_norm_row(term['nonlinear_prefix_error_weighted_norm_per_next_hb'])
                delta_V += self.scaled_jet_box(error,term['width_power']+1,
                                             term['additional_positive_scale_log_enclosure'],norm_only=True)
        return ell,delta_phi,delta_V

    def macro_moments(self,p,q,inlet,delta_phi,delta_V):
        """Actual own history with signed finite-width Volterra feedback.

        Propagate the actual micro endpoint. Integrate all signed linear
        angular/axial prefix modes with their true positive Volterra kernel.
        Exponential and quadratic feedback have separate bounded remainder.
        """
        c=self.ctx;length=p['L']*q
        R1=c.mpf(100) if endpoints(q)==(mp.mpf(1),mp.mpf(1)) else self.bridge.r*c.exp(self.Y*q+c.mpf([0,endpoints(2*self.cap*(1-q))[1]]))
        phi0=p['data']['phi0'].truncate(5);V0=p['data']['V0'].truncate(5)
        start={n:IntervalTaylor(c,[c.mpf(v) for v in values]) for n,values in inlet.items()}
        angular0,axial0=self.contributions(p,'second',c.mpf(2))
        targets=dict(H=phi0,M=V0,K=phi0*V0,A=square(V0),B=square(phi0)/2,C=square(phi0))
        zero=phi0*0;baseline={};feedback={};actual={};ledger={};quadratic={}
        B=p['angular_prefix_norm_per_hb'];exp_bound=c.exp(6*self.cap)
        scales=dict(hydro=c.mpf(0),pressure=self.pressure_log,swirl=self.swirl_log)
        powers=dict(hydro=1,pressure=1,swirl=2)
        for name,rate in RATES.items():
            decay=c.exp(-rate*length);W=(1-decay)/rate
            baseline[name]=start[name]*decay+targets[name]*(1-decay)
            Fterms=[];Vterms=[]
            for term in angular0:
                row=phi0*IntervalTaylor(c,term['normalized_axial_coefficients'])*W
                Fterms.append(self.source_term(row,term['width_power']))
            kernels=volterra_primitive_kernels(c,p['R0'],R1,length,1,rate)
            Fterms.append(self.source_term(phi0*weighted_modes(p['modes']['D_over_R'],kernels)*(-c.mpf('.5')),1))
            # exp(ell)-1-ell: true full-prefix norm <= hb^2 B^2 exp(|ell|)/2.
            Ferr=jet_norm(phi0,self.weight)*B*B*exp_bound*upper(c,abs(W))/2
            Fterms.append(self.source_term(self.symmetric_norm_row(Ferr),2))
            for part,scale in scales.items():
                for term in axial0[part]:
                    Vterms.append(self.source_term(IntervalTaylor(c,term['normalized_axial_coefficients'])*W,
                        term['width_power'],scale,
                        term['nonlinear_prefix_error_weighted_norm_per_next_hb']*upper(c,abs(W))))
                kernels=volterra_primitive_kernels(c,p['R0'],R1,length,powers[part],rate)
                row=-p['quotient']*weighted_modes(p['modes']['drive_'+part],kernels)
                err=jet_norm(p['quotient'],self.weight)*B*exp_bound\
                    *weighted_mode_norm(p['modes']['drive_'+part],kernels,self.weight)
                Vterms.append(self.source_term(row,1,scale,err))
            terms=[]
            fm={'H':c.mpf(2),'M':c.mpf(0),'K':2*V0,'A':c.mpf(0),'B':2*phi0,'C':2*phi0}[name]
            vm={'H':c.mpf(0),'M':c.mpf(1),'K':2*phi0,'A':2*V0,'B':c.mpf(0),'C':c.mpf(0)}[name]
            for multiplier,rows in ((fm,Fterms),(vm,Vterms)):
                for term in rows:
                    signed=IntervalTaylor(c,term['normalized_axial_coefficients'])*multiplier
                    error=term['nonlinear_prefix_error_weighted_norm_per_next_hb']\
                        *(jet_norm(multiplier,self.weight) if isinstance(multiplier,IntervalTaylor) else abs(multiplier))
                    terms.append(self.source_term(signed,term['width_power'],term['additional_positive_scale_log_enclosure'],error))
            feedback[name]=zero
            for term in terms:
                feedback[name]+=self.scaled_jet_box(IntervalTaylor(c,term['normalized_axial_coefficients']),
                    term['width_power'],term['additional_positive_scale_log_enclosure'])
                feedback[name]+=self.scaled_jet_box(self.symmetric_norm_row(term['nonlinear_prefix_error_weighted_norm_per_next_hb']),
                    term['width_power']+1,term['additional_positive_scale_log_enclosure'],norm_only=True)
            quadratic[name]={'H':zero,'M':zero,'K':2*delta_phi*delta_V,
                             'A':square(delta_V),'B':square(delta_phi),'C':square(delta_phi)}[name]*W
            feedback[name]+=quadratic[name]
            actual[name]=baseline[name]+feedback[name]
            ledger[name]=terms
        return dict(actual=actual,baseline=baseline,feedback=feedback,
                    signed_linear_volterra_terms=ledger,quadratic_remainder=quadratic)

    def packet(self, Z, coordinate, chart='macro', root=False):
        c = self.ctx; q = c.mpf(coordinate); lo,hi=endpoints(q)
        if chart not in ('first','second','macro') or (chart=='second' and (lo<1 or hi>2))\
                or (chart!='second' and (lo<0 or hi>1)):
            raise ValueError('Coordinate enclosure in the original chart required')
        p = self.prepare(Z,root)
        angular,axial = self.contributions(p,chart,q)
        ell,delta_phi,delta_V = self.field_changes(p,angular,axial)
        if chart=='first' and hi==0:
            ell=ell*0;delta_phi=delta_phi*0;delta_V=delta_V*0
        if chart=='macro':
            theta = self.bridge.r/100 if lo==hi==1 else c.exp(-self.Y*q-c.mpf([0,endpoints(2*self.cap*(1-q))[1]]))
            radius = c.mpf(100) if lo==hi==1 else self.bridge.r/theta
            formal_y='2hb+q*(log(100/Ra)-2hb)'
        else:
            theta=c.exp(-c.mpf([0,endpoints(self.cap*q)[1]]));radius=self.bridge.r/theta
            formal_y='hb*s'
        phi0=p['data']['phi0'].truncate(5);V0=p['data']['V0'].truncate(5)
        initial={n:j.truncate(5) for n,j in p['data']['moments'].items()}
        history_delta_V=p['uniform_delta_V'] if not (chart=='first' and hi==0) else delta_V*0
        own=own_moment_feedback(phi0,V0,initial,theta,delta_phi,history_delta_V)
        if chart=='macro' and hi!=0:
            inlet=self.packet(Z,2,'second',root)['actual_own_six_moments_axial5']
            own=self.macro_moments(p,q,inlet,delta_phi,history_delta_V)
        exp_ell=ell.exp()
        exp_minus_one=IntervalTaylor(c,[c.expm1(ell[0])]+list(exp_ell.coefficients[1:]))
        point_delta_phi=phi0*exp_minus_one
        phi=phi0+point_delta_phi;V=V0+delta_V
        z=IntervalTaylor.variable(c,p['Z'],5);d=1-square(z);L=1-square(z)*self.core.delta
        Q=(2*z*V-z*own['actual']['M']*(1-self.core.delta)-d*derivative(own['actual']['M']))/L
        pressure=dress(own['actual']['C'],p['inputs']['F0_squared_ratios'])
        return dict(Z=p['Z'],chart=chart,coordinate=q,exact_log_radius_over_Ra=formal_y,
            radius_enclosure_only=radius,theta_enclosure_only=theta,
            signed_angular_integral_terms=angular,signed_axial_integral_terms=axial,
            angular_prefix_weighted_norm_per_hb=p['angular_prefix_norm_per_hb'],
            actual_log_F_over_inlet_axial5=list(ell.coefficients),
            actual_phi_axial5=list(phi.coefficients),actual_raw_V_axial5=list(V.coefficients),
            actual_delta_phi_axial5=list(point_delta_phi.coefficients),actual_delta_V_axial5=list(delta_V.coefficients),
            uniform_prefix_delta_phi_axial5=list(delta_phi.coefficients),
            uniform_prefix_delta_V_axial5=list(history_delta_V.coefficients),
            actual_own_six_moments_axial5={n:list(j.coefficients) for n,j in own['actual'].items()},
            zero_width_baseline_own_six_moments_axial5={n:list(j.coefficients) for n,j in own['baseline'].items()},
            finite_width_feedback_own_six_moments_axial5={n:list(j.coefficients) for n,j in own['feedback'].items()},
            signed_linear_own_moment_volterra_terms=own.get('signed_linear_volterra_terms'),
            nonlinear_quadratic_own_moment_feedback_axial5={n:list(j.coefficients) for n,j in own.get('quadratic_remainder',{}).items()},
            actual_radial_Q_axial4=list(Q.coefficients),
            pressure_axis_axial5=list(p['inputs']['p0'].truncate(5).coefficients),
            actual_pressure_increment_axial5_divided_by_R_F0_squared=list(pressure.coefficients),
            actual_core_atom_inlet={n:list(j.coefficients) for n,j in initial.items()},
            finite_width_comparison_exit={n:j.packet() for n,j in p['end']['moments'].items()},
            same_original_core_atoms_used=True,known_comparison_direction_preserved=True,
            actual_moment_enclosures_from_prescribed_FV=True,
            actual_point_moment_history_recovered=False,
            signed_macro_own_moment_Volterra_feedback_integrated=chart=='macro' and hi!=0,
            exact_core_inlet=chart=='first' and hi==0,
            coordinate_R100_endpoint=chart=='macro' and lo==hi==1,
            R100_functional_join_to_existing_switch_installed=False,
            comparison_moments_substituted_for_actual=False,
            signed_source_integrals_include_both_micro_charts_and_macro=chart=='macro',
            nonlinear_angular_prefix_remainder_controlled=True,
            source_width_not_materialized=True,cap_used_as_field_value=False,
            exact_production_point_parameters_selected=False,
            actual_bridge_mixed4_feedback_installed=False,
            full_implicit_leading_inputs_recomputed=False,temporal_recursion=False)

    def report(self):
        points={}
        for name,Z,root in (('0',0,False),('.5','.5',False),('exact_shared_root',0,True)):
            points[name]=dict(core=self.packet(Z,0,'first',root),
                first_exit=self.packet(Z,1,'first',root),
                second_inlet=self.packet(Z,1,'second',root),
                second_exit=self.packet(Z,2,'second',root),
                macro_inlet=self.packet(Z,0,'macro',root),
                R100=self.packet(Z,1,'macro',root))
            print('Actual finite-width signed bridge and own six histories: '+name,flush=True)
        return dict(actual_five_defect_family_sha256=self.core.family,
            implicit_source_sha256=self.core.source,datum_enclosure_sha256=self.core.datum.datum_sha,
            original_bridge_source='Lei-Ren v2 (4.34)-(4.36), (9.23): comparison q is known; actual moments recover I and stress',
            exact_source_width='hb=epsilon_b=cstar*K^-100',source_log_hb_enclosure=self.logh,
            axial_jet_weight=self.weight,width_enclosure_cap=self.cap,
            exact_macro_source=dict(R0='Ra*exp(2hb)',L='log(100/Ra)-2hb',
                comparison_modes=[0,-1,-2],actual_angular='ell=-.5*integral_0^y chi*Dbar dt',
                actual_axial='V=Va-integral_0^y chi*(phi_a/barphi)*exp(ell)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl)dt'),
            exact_own_moment_ODEs={n:dict(rate=r,source={'H':'2phi','M':'V','K':'2phi*V','A':'V^2','B':'phi^2','C':'phi^2'}[n],
                initial='actual_core_integral_atom_'+n) for n,r in RATES.items()},
            full_finite_width_signed_bridge_integral_enclosures_available=True,
            actual_own_six_moment_feedback_enclosures_available=True,
            signed_macro_own_moment_Volterra_feedback_integrated=True,
            actual_point_moment_history_recovered=False,
            known_direction_not_replaced_by_actual_inertial_stress=True,
            original_control_AST_bindings=self.bindings,packets=points,
            whole_axis_R100=self.packet([-1,1],1,'macro'),
            source_axial_domain=[-1,1],
            exact_swirl_scale_source='F0(Z)^2=Cstar^-2*exp(-2Lambda*G(Z)); Gbar is a bound only',
            swirl_scale_log_enclosure=self.swirl_log,
            correlated_real_chi0_bound_applied_before_model_tail=True,
            real_model_tail_source_replay=dict(model_tail_coefficients_source=sha(Path(atom_module.__file__).name),
                                              source_profile_jet_source=sha(Path(comparison_module.__file__).name),
                                              only_final_model_tail_callback_changed=True,
                                              real_identity='0<=H(Z)^2/(H(Z)^2+sigma^2)<=1'),
            source_log_product_cap_proofs=self.proofs,
            actual_bridge_mixed4_feedback_installed=False,
            R100_R110_actual_feedback_composed=False,
            full_implicit_leading_inputs_recomputed=False,
            global_completed_tensor_admissibility=False,
            exact_production_point_parameters_selected=False,temporal_recursion=False,
            input_hashes=self.hashes)


def run():
    with mp.workdps(400):
        result=CompliantActualBridgeIntegrals().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(_encode(result),indent=2)+'\n',encoding='utf8')
    return result


if __name__=='__main__':
    run()
