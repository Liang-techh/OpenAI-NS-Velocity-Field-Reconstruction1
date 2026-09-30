"""Continuous source pulse and MP-centered full-row integral atoms.

Full integral definitions include all support. Numerical centered evaluations
retain an explicit positive omitted-piece bound and uncertified quadrature
error. These atoms are not yet installed in the global mean/energy solve.
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
        """Same weighted primitive, resolved inside the saddle window.

        The positive contribution before the window is bounded, not erased.
        Outside that window an endpoint-specific evaluator is still required.
        """
        with mp.workdps(self.precision):
            mu=mp.mpf(mu);xi=mp.mpf(xi);band=mp.mpf(band)
            if row not in (1,2) or not 0<mu<=mp.mpf('1e-6'):
                raise ValueError('Centered source rows require small positive mu')
            if xi<=0:
                return dict(value=mp.mpf(0),derivative=mp.mpf(0),exact_support_zero=True)
            if xi>=11:
                result=self.full_row(mu,row,band=band)
                result['derivative']=mp.mpf(0)
                return result
            lam=mp.mpf('.5')-row*mu;k=lam/mu
            u0=mp.root(2/k,3);L=1/u0**2;w=1/mp.sqrt(6*L)
            if band*w>=mp.mpf('.25'):
                raise ValueError('Band exceeds centered small-mu domain')
            lower=((11-xi)/u0-1)/w
            if not -band<lower<band:
                raise ValueError('Partial primitive outside resolved saddle window')
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
        report=dict(rows=rows,value_jet_replay_relative_error=mp.nstr(abs(derivative/jet['derivative']-1),60),
            partial_derivative_diagnostic_mu=str(diagnostic_mu),
            partial_derivative_relative_errors=partial_errors,
            shared_continuous_pulse_definition=True,
            partial_weighted_primitives_implemented='Resolved saddle window and exact support endpoints only',
            installed_in_global_profile=False,finite_energy_certified=False,
            scope='MP pointwise/full-row and saddle-window partial primitives with positive omitted-piece bounds; off-window partial evaluation, energy atom and global installation remain open.')
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(dict(log_refinement=[r['log_precision_refinement'] for r in rows],
            log_input_change=[r['log_change_from_float_centered_input'] for r in rows],
            log_omitted_bound=[r['log_relative_omitted_positive_bound'] for r in rows])),flush=True)
        return report


if __name__=='__main__':run()
