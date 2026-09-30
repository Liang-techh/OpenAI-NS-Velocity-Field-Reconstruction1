"""Continuous end-bump/pulse coefficient solve; incoming data remain inherited.

This solves numerical continuous integral atoms, not a certified exact identity.
No coefficients are installed into the global velocity/mean provider here.
"""
import json
from pathlib import Path
from functools import lru_cache
import mpmath as mp
from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log
from lei_ren_part1_paper_continuous_axial_basis import ContinuousAxialBump
from lei_ren_part1_paper_continuous_axial_pulse import ContinuousAxialPulse


class ContinuousEndCorrection:
    """Point values and accumulated row contributions of the solved bumps."""
    def __init__(self,receipt):
        self.precision=receipt['algebra_precision']
        self.basis=ContinuousAxialBump(precision=self.precision)
        with mp.workdps(self.precision):
            self.mu=mp.mpf(receipt['input_mu'])
            self.c=[from_signed_log(x) for x in receipt['c']]

    def value_jet(self,offset):
        with mp.workdps(self.precision):
            jets=[self.basis.values(mp.mpf(offset)-center) for center in (-3,-1)]
            return {key:sum(c*jet[key] for c,jet in zip(self.c,jets))
                    for key in ('beta','beta_s','beta_ss')}

    def weighted_primitive(self,row,offset):
        with mp.workdps(self.precision):
            if row not in (1,2):raise ValueError('Expected row 1 or 2')
            lam=mp.mpf('.5')-row*self.mu
            return sum(c*mp.exp(lam*center)*self.basis.primitive(lam,mp.mpf(offset)-center)
                       for c,center in zip(self.c,(-3,-1)))

    def weighted_tail(self,row,offset,*,coefficients=None):
        """Integrate the actual end correction from offset to infinity."""
        with mp.workdps(self.precision):
            if row not in (1,2):raise ValueError('Expected row 1 or 2')
            lam=mp.mpf('.5')-row*self.mu
            c=self.c if coefficients is None else coefficients
            return sum(value*mp.exp(lam*center)*self.basis.tail(lam,mp.mpf(offset)-center)
                       for value,center in zip(c,(-3,-1)))


class ContinuousAxialCorrection(ContinuousEndCorrection):
    """One pulse/end provider for values and normalized cumulative row means.

    Incoming rows are supplied explicitly; their uncertainty is not enclosed.
    This is a component adapter, not the physical exterior installation.
    """
    def __init__(self,receipt,base_rows):
        super().__init__(receipt)
        self.pulse=ContinuousAxialPulse(precision=self.precision)
        with mp.workdps(self.precision):
            self.a=mp.mpf(receipt['a_p'])
            self.base=[from_signed_log(x) for x in base_rows]

    def pulse_value_jet(self,xi):
        with mp.workdps(self.precision):
            return {key:self.a*value for key,value in self.pulse.value_jet(xi).items()}

    @lru_cache(maxsize=8)
    def terminal_balance(self,row,evaluation_precision=None):
        """Cache the evaluated full-atom residual; never substitute zero."""
        with mp.workdps(max(self.precision,evaluation_precision or self.precision)):
            if row not in (1,2):raise ValueError('Expected row 1 or 2')
            atom=self.pulse.full_row(self.mu,row)
            integral=mp.exp(mp.mpf(atom['log_normalized_pulse_integral']))
            end=self.weighted_primitive(row,0)
            return dict(value=self.base[row-1]+self.a*integral+end,
                        pulse_integral=integral,pulse_atom=atom,
                        full_end_primitive=end,materialized_residual_retained=True,
                        exact_functional_identity='M c + b + a p = 0 for exact continuous atoms',
                        functional_identity_certified_for_materialized_coefficients=False)

    def cumulative_row(self,row,*,xi=None,end_offset=None):
        if (xi is None)==(end_offset is None):
            raise ValueError('Supply pulse xi or end-bump offset')
        with mp.workdps(self.precision):
            if row not in (1,2):raise ValueError('Expected row 1 or 2')
            atom=self.pulse.partial_row(self.mu,row,xi) if xi is not None else self.pulse.full_row(self.mu,row)
            integral=(mp.exp(mp.mpf(atom['log_normalized_pulse_integral']))
                      if 'log_normalized_pulse_integral' in atom else mp.mpf(0))
            end=self.weighted_primitive(row,end_offset) if end_offset is not None else mp.mpf(0)
            if end_offset is not None and mp.mpf(end_offset)>=mp.mpf('-3.15'):
                # Pulse support is complete throughout the end-bump region.
                # N(s)=rho-int_s^infinity end, where rho is the evaluated
                # full balance. Keeping rho makes this a representation of
                # the same materialized field, not an imposed terminal mask.
                mean=self.terminal_balance(row)['value']-self.weighted_tail(row,end_offset)
                method='retained_terminal_balance_minus_direct_end_tail'
            else:
                mean=self.base[row-1]+self.a*integral+end
                method='forward_continuous_primitive'
            relative_bound=atom.get('log_relative_omitted_absolute_bound',
                                    atom.get('log_relative_omitted_positive_bound'))
            correction_sign=atom.get('omitted_correction_sign',1)
            logbound=(mp.log(self.a)+mp.mpf(atom['log_normalized_pulse_integral'])+
                      mp.mpf(relative_bound) if relative_bound is not None else None)
            return dict(nominal=signed_log(mean,self.precision),method=method,
                log_pulse_omitted_absolute_bound=mp.nstr(logbound,self.precision) if logbound is not None else None,
                pulse_omitted_correction_sign=correction_sign,
                log_pulse_omitted_positive_bound=(mp.nstr(logbound,self.precision)
                    if logbound is not None and correction_sign==1 else None),
                incoming_uncertainty_enclosed=False,quadrature_enclosure_certified=False,
                terminal_mean_forced_zero=False)


def solve_continuous(source, *, precision=200, energy_precision=100):
    from lei_ren_part1_paper_continuous_pulse_energy import continuous_pulse_energy
    with mp.workdps(precision):
        old=source['axial'];mu=mp.mpf(old['input_mu'])
        target=mp.mpf(old['energy_target'])
        basis=ContinuousAxialBump(precision=precision)
        pulse_provider=ContinuousAxialPulse(precision=precision)
        matrix,det=basis.matrix(mu)
        base=[from_signed_log(x) for x in old['linear_rhs_inputs']['base']]
        pulse_rows=[pulse_provider.full_row(mu,i) for i in (1,2)]
        pulses=[mp.exp(mp.mpf(x['log_normalized_pulse_integral'])) for x in pulse_rows]
        def affine(rhs):
            return [(-rhs[0]*matrix[1][1]+matrix[0][1]*rhs[1])/det,
                    (-matrix[0][0]*rhs[1]+rhs[0]*matrix[1][0])/det]
        u=affine(base);v=affine(pulses)
        gram=basis.energy_gram(mu)
        K=[mp.exp(-26+i*mu)*gram for i in (6,2)]
        Kp=continuous_pulse_energy(precision=energy_precision)
        quadratic=Kp+mu*sum(k*b*b for k,b in zip(K,v))
        linear=2*mu*sum(k*a*b for k,a,b in zip(K,u,v))
        constant=mu*sum(k*a*a for k,a in zip(K,u))-target
        discriminant=linear*linear-4*quadratic*constant
        if constant>=0 or discriminant<=0:
            raise ArithmeticError('Actual inherited target has no positive source branch')
        a=-2*constant/(linear+mp.sqrt(discriminant))
        if not mp.mpf('.9')<a<mp.mpf('1.2'):
            raise ArithmeticError('Continuous source energy root outside admissible source interval')
        c=[x+a*y for x,y in zip(u,v)]
        row_errors=[];old_continuous_errors=[]
        old_c=[from_signed_log(x) for x in old['c']];old_a=mp.mpf(old['a_p'])
        for i in range(2):
            rhs=base[i]+a*pulses[i]
            row_errors.append(mp.nstr(abs((sum(matrix[i][j]*c[j] for j in range(2))+rhs)/rhs),40))
            old_rhs=base[i]+old_a*pulses[i]
            old_continuous_errors.append(mp.nstr(abs((sum(matrix[i][j]*old_c[j] for j in range(2))+old_rhs)/old_rhs),40))
        energy=quadratic*a*a+linear*a+constant
        return dict(input_mu=mp.nstr(mu,precision),a_p=mp.nstr(a,precision),
            c=[signed_log(x,precision) for x in c],
            coefficient_relative_changes=[mp.nstr(x/y-1,60) for x,y in zip(c,old_c)],
            amplitude_relative_change=mp.nstr(a/old_a-1,60),
            linear_matrix=[[mp.nstr(x,precision) for x in row] for row in matrix],
            pulse_rows=pulse_rows,K_p=mp.nstr(Kp,energy_precision),
            K_bump=[mp.nstr(x,precision) for x in K],
            linear_relative_replay=row_errors,
            old_coefficients_continuous_row_relative_defects=old_continuous_errors,
            energy_relative_replay=mp.nstr(abs(energy/target),40),
            energy_atom_precision=energy_precision,algebra_precision=precision,
            energy_atom_quadrature_order=max(48,min(128,energy_precision*4//5)),
            energy_atom_precision_is_not_certified_accuracy=True,
            inherited_base_rows=True,inherited_energy_target=True,
            quadrature_enclosure_certified=False,installed_in_global_profile=False,
            global_mean_closed=False,finite_energy_certified=False)


def run():
    source=json.loads(Path(__file__).with_name('lei_ren_part1_paper_seeded_shared_candidate.json').read_text())
    report=solve_continuous(source)
    correction=ContinuousEndCorrection(report)
    full_component=ContinuousAxialCorrection(report,source['axial']['linear_rhs_inputs']['base'])
    report['continuous_component_cumulative_rows_after_support']=[
        full_component.cumulative_row(row,end_offset=0) for row in (1,2)]
    with mp.workdps(report['algebra_precision']):
        offset=mp.mpf('-3.0375');jet=correction.value_jet(offset)
        lam=mp.mpf('.5')-correction.mu
        derivative=mp.exp(lam*offset)*jet['beta']
        errors=[]
        for h in (mp.mpf('1e-5'),mp.mpf('5e-6')):
            values={i:correction.weighted_primitive(1,offset+i*h) for i in (-2,-1,1,2)}
            numerical=(values[-2]-8*values[-1]+8*values[1]-values[2])/(12*h)
            errors.append(mp.nstr(abs(numerical/derivative-1),40))
        report['end_component_primitive_derivative_errors']=errors
        report['end_component_provider_available']=True
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('a_p','linear_relative_replay',
        'old_coefficients_continuous_row_relative_defects','energy_relative_replay',
        'amplitude_relative_change','coefficient_relative_changes')}),flush=True)
    return report


if __name__=='__main__':run()
