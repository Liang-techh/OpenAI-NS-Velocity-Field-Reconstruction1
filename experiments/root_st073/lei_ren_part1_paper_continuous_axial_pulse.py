"""Continuous source pulse and MP-centered full-row integral atoms.

Full integral definitions include all support. Numerical centered evaluations
retain an explicit positive omitted-piece bound and uncertified quadrature
error. Startup partial rows use an endpoint-scaled integration-by-parts
primitive with the same ``value_jet`` source and a positive tail bound. These
atoms are not yet installed in the global mean/energy solve.
"""
import json
from pathlib import Path
import mpmath as mp


class ContinuousAxialPulse:
    def __init__(self,*,precision=160):self.precision=precision

    @staticmethod
    def sigma_pair(x):
        if x<=0:return mp.mpf(0),mp.mpf(0)
        if x>=1:return mp.mpf(1),mp.mpf(0)
        phase=1/x**2-1/(1-x)**2
        if phase>=0:
            tiny=mp.exp(-phase);value=tiny/(1+tiny)
        else:
            tiny=mp.exp(phase);value=1/(1+tiny)
        product=tiny/(1+tiny)**2
        derivative=product*(2/x**3+2/(1-x)**3)
        return value,derivative

    def value_jet(self,xi):
        with mp.workdps(self.precision):
            xi=mp.mpf(xi)
            if xi<=0 or xi>=11:return dict(value=mp.mpf(0),derivative=mp.mpf(0))
            if xi>=mp.mpf('.02'):
                primitive=xi-mp.mpf('.01');primitive_derivative=mp.mpf(1)
            else:
                q=50*xi
                if q<=mp.mpf('.5'):
                    # On the flat left endpoint, sigma(u) <= exp(-phase(q))
                    # for 0 <= u <= q.  If even the resulting integral upper
                    # bound is below the active MP resolution, adaptive
                    # quadrature must not manufacture a tiny nonzero pulse.
                    phase=1/q**2-1/(1-q)**2
                    log_upper=mp.log(q)-phase-mp.log(50)
                    if log_upper < -mp.mpf(max(8,self.precision-8))*mp.log(10):
                        raise ArithmeticError(
                            'Startup pulse is below the active MP precision; '
                            'a smaller-xi log evaluator is required'
                        )
                    primitive=mp.quad(lambda u:self.sigma_pair(u)[0],[0,q])/50
                else:
                    # Symmetric integral, preserving the positive correction.
                    primitive=xi-mp.mpf('.01')+mp.quad(
                        lambda u:self.sigma_pair(u)[0],[0,1-q])/50
                primitive_derivative=self.sigma_pair(q)[0]
                if primitive==0:
                    raise ArithmeticError('Nonzero startup pulse unresolved; a bounded endpoint integral is required')
            cutoff,cutoff_derivative=self.sigma_pair(11-xi)
            return dict(value=primitive*cutoff,
                derivative=primitive_derivative*cutoff-primitive*cutoff_derivative)

    def full_row(self,mu,row,*,band=48):
        """Evaluate (1/mu) int[0,11] exp(-lambda*(13-v)/mu) gp(v) dv."""
        with mp.workdps(self.precision):
            mu=mp.mpf(mu);band=mp.mpf(band)
            if row not in (1,2) or not 0<mu<=mp.mpf('1e-6'):
                raise ValueError('Centered source rows require small positive mu')
            if band<=0:
                raise ValueError('Centered source band must be positive')
            lam=mp.mpf('.5')-row*mu;k=lam/mu
            u0=mp.root(2/k,3);L=1/u0**2;w=1/mp.sqrt(6*L)
            if band*w>=mp.mpf('.25'):raise ValueError('Band exceeds centered small-mu domain')
            def centered(x):
                v=1+x*w;u=u0*v
                difference=x*x*(2*v+1)/(6*v*v)
                smooth=1/(1-u)**2
                return (mp.mpf('10.99')-u)*mp.exp(smooth-difference)/(1+mp.exp(-1/u**2+smooth))
            total=mp.quad(centered,[-band,-8,-4,0,4,8,band])
            logvalue=-2*k-3*L+2*mp.log(u0)-mp.log(6)/2-mp.log(mu)+mp.log(total)
            phases=[band**2*(2*(1+sign*band*w)+1)/(6*(1+sign*band*w)**2) for sign in (-1,1)]
            bounds=[mp.log(mp.mpf('5.5'))+4-mp.log(mu)-2*k-3*L-min(phases),
                mp.log(mp.mpf('5.5'))-mp.log(mu)-mp.mpf('2.5')*k,
                mp.log(11/lam)-3*k]
            largest=max(bounds);logbound=largest+mp.log(sum(mp.exp(v-largest) for v in bounds))
            n=lambda x:mp.nstr(x,self.precision)
            return dict(row=row,mu=n(mu),lambda_value=n(lam),band=n(band),
                log_normalized_pulse_integral=n(logvalue),
                log_relative_omitted_positive_bound=n(logbound-logvalue),
                centered_integral=n(total),log_omitted_piece_bounds=[n(v) for v in bounds],
                continuous_definition='Integral of the same value_jet pulse over full support [0,11]',
                quadrature_enclosure_certified=False)

    def partial_row(self,mu,row,xi,*,band=48):
        """Same weighted primitive on startup, plateau, and cutoff.

        A finite centered band retains a positive omitted contribution bound;
        it is never replaced by zero. Extremely small endpoint pulses can
        remain below the active MP precision and raise explicitly.
        """
        with mp.workdps(self.precision):
            mu=mp.mpf(mu);xi=mp.mpf(xi);band=mp.mpf(band)
            if row not in (1,2) or not 0<mu<=mp.mpf('1e-6'):
                raise ValueError('Centered source rows require small positive mu')
            if band<=0:
                raise ValueError('Centered source band must be positive')
            if xi<=0:
                return dict(value=mp.mpf(0),derivative=mp.mpf(0),exact_support_zero=True)
            if xi>=11:
                result=self.full_row(mu,row,band=band)
                result['derivative']=mp.mpf(0)
                return result
            lam=mp.mpf('.5')-row*mu;k=lam/mu
            if xi<=mp.mpf('.02'):
                # Endpoint scaling resolves the startup contribution without
                # subtracting the exponentially small primitive from a bulk
                # formula.  With s=k*(xi-v),
                #
                # J(xi)=exp(-k*(13-xi))/(mu*k)
                #        * integral_0^(k*xi) exp(-s) gp(xi-s/k) ds.
                #
                # The source value is the shared value_jet, including its
                # nested MP startup primitive.  A finite band truncates the
                # positive sigma correction after integration by parts; the
                # omitted correction therefore has negative sign.
                source_at_endpoint=self.value_jet(xi)['value']
                if source_at_endpoint==0:
                    raise ArithmeticError(
                        'Startup endpoint pulse is below the active MP precision; '
                        'a smaller-xi log evaluator is required'
                    )
                span_total=k*xi
                span=min(span_total,band)
                if span<=0:
                    raise ArithmeticError(
                        'Startup endpoint span collapsed at positive xi; '
                        'a smaller-xi log evaluator is required'
                    )
                # Integrating the raw startup value_jet at every outer node
                # would nest a second adaptive quadrature.  Integrate by
                # parts once instead.  If G(v)=gp(v), G'=sigma(50v), then
                #
                # (1/mu) int exp(-k*(xi-v))G(v)dv
                #   = G(xi)/lambda - mu/lambda^2 * int exp(-s)
                #       sigma(50*(xi-s/k)) ds.
                # The remaining endpoint integral is single-level and uses
                # the same sigma_pair definition as value_jet.
                def sigma_scaled(s):
                    v=xi-s/k
                    if v<=0:return mp.mpf(0)
                    return mp.exp(-s)*self.sigma_pair(50*v)[0]
                knots=[mp.mpf(0)]+[point for point in (4,8,16,32)
                    if mp.mpf(point)<span]+[span]
                sigma_total=mp.quad(sigma_scaled,knots)
                lam=mp.mpf('.5')-row*mu
                endpoint_homogeneous=(source_at_endpoint/lam-
                    mu*sigma_total/(lam*lam))
                if endpoint_homogeneous<=0:
                    raise ArithmeticError(
                        'Startup endpoint integration lost a positive weighted primitive; '
                        'a higher precision or smaller-xi evaluator is required'
                    )
                logvalue=-k*(13-xi)+mp.log(endpoint_homogeneous)
                # The retained integration-by-parts correction is a
                # subtraction.  When span < k*xi, the missing sigma tail
                # therefore gives an absolute error with negative sign,
                # rather than a positive omitted source-mass contribution.
                if span_total>span:
                    sigma_at_endpoint=self.sigma_pair(50*xi)[0]
                    if sigma_at_endpoint==0:
                        raise ArithmeticError(
                            'Startup sigma endpoint is below the active MP precision; '
                            'a smaller-xi log evaluator is required'
                        )
                    logbound=(-k*(13-xi)+mp.log(mu*sigma_at_endpoint/(lam*lam))-span)
                    log_relative_bound=logbound-logvalue
                    omitted_sign=-1
                else:
                    log_relative_bound=None
                    omitted_sign=0
                log_source=-k*(13-xi)+mp.log(source_at_endpoint)-mp.log(mu)
                return dict(log_normalized_pulse_integral=mp.nstr(logvalue,self.precision),
                    log_derivative=mp.nstr(log_source,self.precision),
                    log_relative_omitted_positive_bound=None,
                    log_relative_omitted_absolute_bound=(
                        mp.nstr(log_relative_bound,self.precision)
                        if log_relative_bound is not None else None),
                    omitted_correction_sign=omitted_sign,
                    xi=mp.nstr(xi,self.precision),region='startup_endpoint_scaled',
                    startup_span=mp.nstr(span,self.precision),
                    startup_total_span=mp.nstr(span_total,self.precision),
                    quadrature_enclosure_certified=False,
                    continuous_definition=(
                        'Integral of value_jet over [0,xi] using s=k*(xi-v) '
                        'endpoint scaling with a signed omitted IBP correction bound'),
                    startup_omitted_bound='For s>=startup_span, sigma(50v)<=sigma(50xi) '
                        'gives a negative omitted IBP correction bounded by '
                        'exp(-k*(13-xi))*mu*sigma(50xi)*exp(-startup_span)/lambda^2; '
                        'quadrature error is not enclosed')
            if mp.mpf('.02')<xi<=10:
                # Exact antiderivative on the plateau; the positive startup
                # atom is bounded explicitly, not claimed identically zero.
                delta=xi-mp.mpf('.02')
                lower_factor=mp.mpf('.01')/k-1/k**2
                bracket=delta/k+lower_factor*(-mp.expm1(-k*delta))
                logvalue=-k*(13-xi)-mp.log(mu)+mp.log(bracket)
                logbound=mp.log(mp.mpf('.02')/(mu*k))-k*mp.mpf('12.98')
                return dict(log_normalized_pulse_integral=mp.nstr(logvalue,self.precision),
                    log_derivative=mp.nstr(-k*(13-xi)+mp.log(xi-mp.mpf('.01'))-mp.log(mu),self.precision),
                    log_relative_omitted_positive_bound=mp.nstr(logbound-logvalue,self.precision),
                    xi=mp.nstr(xi,self.precision),region='plateau',
                    quadrature_enclosure_certified=False,
                    continuous_definition='Integral of value_jet over [0,xi]; analytic plateau plus bounded positive startup')
            u0=mp.root(2/k,3);L=1/u0**2;w=1/mp.sqrt(6*L)
            if band*w>=mp.mpf('.25'):
                raise ValueError('Band exceeds centered small-mu domain')
            lower=((11-xi)/u0-1)/w
            if lower<=-band:
                result=self.full_row(mu,row,band=band)
                result.update(region='after_saddle_window',
                    xi=mp.nstr(xi,self.precision),
                    log_derivative=mp.nstr(-k*(13-xi)+mp.log(self.value_jet(xi)['value'])-mp.log(mu),self.precision),
                    continuous_definition='Integral through xi; full central window plus bounded positive remainder')
                return result
            if lower>=band:
                return self._cutoff_prefix(mu,k,xi,band)
            def centered(x):
                v=1+x*w;u=u0*v
                difference=x*x*(2*v+1)/(6*v*v)
                smooth=1/(1-u)**2
                return (mp.mpf('10.99')-u)*mp.exp(smooth-difference)/(1+mp.exp(-1/u**2+smooth))
            knots=[lower]+[x for x in (-8,-4,0,4,8) if lower<x<band]+[band]
            total=mp.quad(centered,knots)
            logvalue=-2*k-3*L+2*mp.log(u0)-mp.log(6)/2-mp.log(mu)+mp.log(total)
            full=self.full_row(mu,row,band=band)
            bounds=[mp.mpf(x) for x in full['log_omitted_piece_bounds']]
            largest=max(bounds)
            logbound=largest+mp.log(sum(mp.exp(x-largest) for x in bounds))
            logderivative=-k*(13-xi)+mp.log(self.value_jet(xi)['value'])-mp.log(mu)
            return dict(log_normalized_pulse_integral=mp.nstr(logvalue,self.precision),
                log_derivative=mp.nstr(logderivative,self.precision),
                log_relative_omitted_positive_bound=mp.nstr(logbound-logvalue,self.precision),
                xi=mp.nstr(xi,self.precision),centered_lower=mp.nstr(lower,self.precision),
                quadrature_enclosure_certified=False,
                continuous_definition='Integral of value_jet over [0,xi] with the full-row weight')

    def _cutoff_prefix(self,mu,k,xi,band):
        """Resolve an increasing weighted cutoff integral at its upper endpoint."""
        u=11-xi;extent=max(mp.mpf(64),band)
        n=lambda x:mp.nstr(x,self.precision)
        if u>=mp.mpf('.5'):
            # In this half of the cutoff gp is O(1); exponential weight sets
            # the endpoint width. A global gp<=11 gives the omitted bound.
            total=mp.quad(lambda s:mp.exp(-s)*self.value_jet(xi-s/k)['value'],
                          [0,4,16,extent])
            logvalue=-k*(13-xi)-mp.log(mu*k)+mp.log(total)
            bounds=[-k*(13-xi)-mp.log(mu*k)+mp.log(11)-extent]
        else:
            slope=k-2/u**3
            if slope<=0:raise ArithmeticError('Endpoint branch crossed saddle')
            if u+extent/slope>=mp.mpf('.5'):
                raise ArithmeticError('Endpoint window exceeds flat cutoff half')
            def scaled(s):
                du=s/slope;v=u+du
                # Stable reciprocal-square difference, avoiding large close
                # exponent subtraction for the actual extreme-mu candidate.
                phase=-k*du+du*(2*u+du)/(u*u*v*v)
                smooth=1/(1-v)**2
                return (mp.mpf('10.99')-v)*mp.exp(phase+smooth)/(1+mp.exp(-1/v**2+smooth))
            total=mp.quad(scaled,[0,4,16,extent])
            common=-2*k-k*u-1/u**2-mp.log(mu)
            logvalue=common-mp.log(slope)+mp.log(total)
            # ku+1/u^2 is convex with derivative >=slope for v>=u.
            # On v<=.5, the logistic smooth exponent is <=4.
            bounds=[common+mp.log(11)+4-mp.log(slope)-extent,
                    mp.log(mp.mpf('5.5')/mu)-mp.mpf('2.5')*k,
                    mp.log(11/(mu*k))-3*k]
        largest=max(bounds)
        logbound=largest+mp.log(sum(mp.exp(x-largest) for x in bounds))
        return dict(log_normalized_pulse_integral=n(logvalue),
            log_derivative=n(-k*(13-xi)+mp.log(self.value_jet(xi)['value'])-mp.log(mu)),
            log_relative_omitted_positive_bound=n(logbound-logvalue),
            log_omitted_piece_bounds=[n(x) for x in bounds],xi=n(xi),
            region='cutoff_endpoint',quadrature_enclosure_certified=False,
            continuous_definition='Integral through xi of shared pulse; endpoint quadrature plus positive omitted-piece bound')


def run():
    source=json.loads(Path(__file__).with_name('lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    pulse=ContinuousAxialPulse(precision=160);higher=ContinuousAxialPulse(precision=200)
    rows=[]
    with mp.workdps(200):
        mu=mp.mpf(source['axial']['input_mu'])
        for index in (1,2):
            row=pulse.full_row(mu,index);fine=higher.full_row(mu,index)
            old=source['axial']['linear_rhs_inputs']['pulse'][index-1]
            row['log_precision_refinement']=mp.nstr(mp.mpf(row['log_normalized_pulse_integral'])-mp.mpf(fine['log_normalized_pulse_integral']),60)
            row['log_change_from_float_centered_input']=mp.nstr(mp.mpf(row['log_normalized_pulse_integral'])-mp.mpf(old['log_normalized_pulse_integral']),60)
            rows.append(row)
        x=mp.mpf('10.5');jet=pulse.value_jet(x);h=mp.mpf('1e-6')
        derivative=(pulse.value_jet(x-2*h)['value']-8*pulse.value_jet(x-h)['value']+
            8*pulse.value_jet(x+h)['value']-pulse.value_jet(x+2*h)['value'])/(12*h)
        diagnostic_mu=mp.mpf('1e-12')
        center=11-mp.root(2/((mp.mpf('.5')-diagnostic_mu)/diagnostic_mu),3)
        central=pulse.partial_row(diagnostic_mu,1,center,band=16)
        partial_errors=[]
        for step in (mp.mpf('1e-12'),mp.mpf('5e-13')):
            left=pulse.partial_row(diagnostic_mu,1,center-step,band=16)
            right=pulse.partial_row(diagnostic_mu,1,center+step,band=16)
            scale=mp.mpf(central['log_derivative'])
            ratio=(mp.exp(mp.mpf(right['log_normalized_pulse_integral'])-scale)-
                mp.exp(mp.mpf(left['log_normalized_pulse_integral'])-scale))/(2*step)
            partial_errors.append(mp.nstr(abs(ratio-1),60))
        plateau=pulse.partial_row(diagnostic_mu,1,mp.mpf(5))
        plateau_errors=[]
        for step in (mp.mpf('1e-17'),mp.mpf('5e-18')):
            left=pulse.partial_row(diagnostic_mu,1,5-step)
            right=pulse.partial_row(diagnostic_mu,1,5+step)
            scale=mp.mpf(plateau['log_derivative'])
            ratio=(mp.exp(mp.mpf(right['log_normalized_pulse_integral'])-scale)-
                mp.exp(mp.mpf(left['log_normalized_pulse_integral'])-scale))/(2*step)
            plateau_errors.append(mp.nstr(abs(ratio-1),60))
        endpoint_x=11-2*mp.root(2/((mp.mpf('.5')-diagnostic_mu)/diagnostic_mu),3)
        endpoint=pulse.partial_row(diagnostic_mu,1,endpoint_x,band=16)
        endpoint_errors=[]
        for step in (mp.mpf('1e-17'),mp.mpf('5e-18')):
            left=pulse.partial_row(diagnostic_mu,1,endpoint_x-step,band=16)
            right=pulse.partial_row(diagnostic_mu,1,endpoint_x+step,band=16)
            scale=mp.mpf(endpoint['log_derivative'])
            ratio=(mp.exp(mp.mpf(right['log_normalized_pulse_integral'])-scale)-
                mp.exp(mp.mpf(left['log_normalized_pulse_integral'])-scale))/(2*step)
            endpoint_errors.append(mp.nstr(abs(ratio-1),60))
        candidate_k=(mp.mpf('.5')-mu)/mu
        candidate_cutoff=[pulse.partial_row(mu,1,x) for x in
            (mp.mpf('10.25'),mp.mpf('10.75'),11-2*mp.root(2/candidate_k,3),
             11-mp.root(2/candidate_k,3)/2)]
        # Startup samples are evaluated at one moderate and the actual shared
        # candidate mu.  Keep the receipts compact because the logs can carry
        # 10^28-scale exponents.  The finite-band tail bound is positive and
        # explicit; it is not a claim of uniform startup accuracy.
        startup_receipts={}
        startup_derivative_diagnostics={}
        startup_points=(mp.mpf('.005'),mp.mpf('.01'),mp.mpf('.015'),mp.mpf('.02'))
        startup_cases=(('moderate_mu_1e-6',mp.mpf('1e-6')),('actual_candidate_mu',mu))
        for case_name,startup_mu in startup_cases:
            by_row={}
            for index in (1,2):
                compact=[]
                for point in startup_points:
                    atom=pulse.partial_row(startup_mu,index,point,band=16)
                    compact.append(dict(xi=atom['xi'],region=atom['region'],
                        log_normalized_pulse_integral=atom['log_normalized_pulse_integral'],
                        log_derivative=atom['log_derivative'],
                        log_relative_omitted_positive_bound=atom['log_relative_omitted_positive_bound'],
                        log_relative_omitted_absolute_bound=atom['log_relative_omitted_absolute_bound'],
                        omitted_correction_sign=atom['omitted_correction_sign'],
                        startup_span=atom['startup_span'],
                        startup_total_span=atom['startup_total_span'],
                        quadrature_enclosure_certified=atom['quadrature_enclosure_certified']))
                by_row[str(index)]=compact
                probe=mp.mpf('.01');central_start=pulse.partial_row(startup_mu,index,probe,band=16)
                scale=mp.mpf(central_start['log_derivative'])
                errors=[]
                startup_k=(mp.mpf('.5')-startup_mu)/startup_mu
                # Resolve the endpoint exponential rather than taking a
                # fixed xi step.  At the actual candidate, a 1e-6 step would
                # move the log by O(1e22) and overflow the scaled diagnostic.
                step_base=min(mp.mpf('1e-6'),mp.mpf('.05')/startup_k)
                for step in (step_base,step_base/2):
                    left=pulse.partial_row(startup_mu,index,probe-step,band=16)
                    right=pulse.partial_row(startup_mu,index,probe+step,band=16)
                    ratio=(mp.exp(mp.mpf(right['log_normalized_pulse_integral'])-scale)-
                        mp.exp(mp.mpf(left['log_normalized_pulse_integral'])-scale))/(2*step)
                    errors.append(mp.nstr(abs(ratio-1),60))
                startup_derivative_diagnostics.setdefault(case_name,{})[str(index)]=dict(
                    xi=mp.nstr(probe,60),steps=[mp.nstr(x,60) for x in (step_base,step_base/2)],
                    relative_errors=errors,
                    comparison='Centered derivative of the finite-band log primitive versus its direct source-term log_derivative; signed omitted IBP correction and quadrature effects remain',
                    log_relative_omitted_absolute_bound=central_start['log_relative_omitted_absolute_bound'],
                    omitted_correction_sign=central_start['omitted_correction_sign'])
            startup_receipts[case_name]=by_row
        report=dict(rows=rows,value_jet_replay_relative_error=mp.nstr(abs(derivative/jet['derivative']-1),60),
            partial_derivative_diagnostic_mu=str(diagnostic_mu),
            partial_derivative_relative_errors=partial_errors,
            plateau_primitive_derivative_relative_errors=plateau_errors,
            plateau_primitive_receipt=plateau,
            cutoff_endpoint_derivative_relative_errors=endpoint_errors,
            cutoff_endpoint_primitive_receipt=endpoint,
            actual_candidate_cutoff_partial_receipts=candidate_cutoff,
            shared_continuous_pulse_definition=True,
            startup_partial_receipts=startup_receipts,
            startup_derivative_diagnostics=startup_derivative_diagnostics,
            partial_weighted_primitives_implemented='Startup endpoint-scaled, plateau, full cutoff via endpoint/saddle/after-window branches, and exact support endpoints',
            installed_in_global_profile=False,finite_energy_certified=False,
            scope='MP pointwise/full-row and startup, plateau, and cutoff partial primitives with positive omitted-piece bounds where applicable and signed startup correction bounds; ordinary startup samples resolve at the tested precision, while tiny endpoint pulses, quadrature enclosure, global installation, and finite-energy claims remain open. Continuous energy atom exists in separate provider.')
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(dict(log_refinement=[r['log_precision_refinement'] for r in rows],
            log_input_change=[r['log_change_from_float_centered_input'] for r in rows],
            log_omitted_bound=[r['log_relative_omitted_positive_bound'] for r in rows])),flush=True)
        return report


if __name__=='__main__':run()
