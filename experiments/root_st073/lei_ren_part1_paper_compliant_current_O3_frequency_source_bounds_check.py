"""Independent original sigma/cutoff/base derivatives versus whole caps."""
import json
from pathlib import Path
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_O3_frequency_source_bounds as source


def run():
    data=json.loads((source.HERE/source.NAME).read_bytes())
    for name,digest in data['input_hashes'].items():
        if source.sha(name)!=digest:raise ValueError('Actual whole-source norm program changed: '+name)
    theorem=source.exact_source_bound_theorem()
    if source._encode(theorem)!=data['exact_actual_modulation_source_norm_theorem']:
        raise ValueError('Whole-support source theorem differs')
    mu=s.Symbol('same_actual_mu',positive=True,finite=True)
    N=s.Symbol('finite_integer_N',positive=True,integer=True,finite=True)
    Ad=s.Symbol('same_actual_amplitude_at_Rd',positive=True,finite=True)
    if source._encode(source.bounded_modulation(mu,N,Ad))!=data['actual_source_norm_and_modulation_bounds']:
        raise ValueError('Actual published source norms or N envelopes differ')
    if data.get('whole_support_original_cutoff_and_base_source_norm_formulas_available') is not True:
        raise ValueError('Whole-support source formula scope missing')
    for flag in ('actual_numeric_graph_parameter_enclosures_substituted',
                 'repair_completed_tensor_common_N_modified_cones_certified'):
        if data.get(flag) is not False:
            raise ValueError('Factored source formulas exceed completed scope')
    tested_sigma,tested_cutoff,tested_base=0,0,0
    with mp.workdps(85):
        def sigma(x):
            if x<=0:return mp.mpf(0)
            if x>=1:return mp.mpf(1)
            if x>mp.mpf('.5'):return 1-sigma(1-x)
            odds=1/(1-x)**2-1/x**2
            e=mp.exp(odds)
            return e/(1+e)
        def cutoff(t):return sigma(t+2)*(1-sigma(4*t-1))
        def real(expr):return mp.mpf(str(s.N(expr,82)))
        caps=source.sigma_source_caps()
        C=source.actual_source_norms(s.Rational(1,37),1)['actual_cutoff_ordinary_logR_majorants']
        for value in ('.001','.1','.25','.5','.75','.99'):
            x=mp.mpf(value)
            for j in range(5):
                if abs(mp.diff(sigma,x,j))>real(caps[j])*(1+mp.mpf('1e-65')):
                    raise ArithmeticError('Independent original sigma exceeds whole source cap')
                tested_sigma+=1
        for value in ('-2','-1.731','-1.127','-1','0','.25','.337','.499','.5'):
            t=mp.mpf(value)
            for j in range(5):
                if abs(mp.diff(cutoff,t,j))>real(C[j])*(1+mp.mpf('1e-65')):
                    raise ArithmeticError('Independent actual cutoff exceeds whole source cap')
                tested_cutoff+=1
        t,z,m=s.symbols('t z mu',real=True)
        J=s.Function('actual_J')(t)
        f=s.exp(-t/2-m*J)
        frows=[s.diff(f,t,j) for j in range(5)]
        qrows=[s.diff(1/(1+z*z),z,k) for k in range(6)]
        for muv in (s.Rational(1,100),s.Integer(1),s.Integer(7)):
            norms=source.actual_source_norms(muv,1)
            for tv in (s.Rational(-2),s.Rational(-371,1000),s.Rational(0),s.Rational(127,1000),s.Rational(1,2)):
                tm=real(tv)
                jm=mp.mpf(0) if tm<=0 else mp.quad(sigma,[0,tm])
                derivative_map={J:s.Float(str(jm),82)}
                for j in range(1,5):
                    derivative_map[s.diff(J,t,j)]=s.Float(str(mp.diff(sigma,tm,j-1) if tm>0 else mp.mpf(0)),82)
                values=[abs(real(row.subs(derivative_map,simultaneous=True).subs({t:tv,m:muv}).doit())) for row in frows]
                for zv in (s.Rational(-1),s.Rational(237,1000),s.Rational(1)):
                    for j in range(5):
                        for k in range(6):
                            actual=values[j]*abs(real(qrows[k].subs(z,zv)))
                            cap=real(norms['original_theta_over_Pstar_logR_axial_majorants'][j][k])
                            if actual>cap*(1+mp.mpf('1e-65')):
                                raise ArithmeticError('Independent actual base derivative exceeds whole norm')
                            tested_base+=1
    rejected=0
    for m,a in ((0,1),(-1,1),(1,0),(1,s.oo)):
        try:source.actual_source_norms(m,a)
        except ValueError:rejected+=1
        else:raise ArithmeticError('Invalid same-source parameter accepted')
    for m,a in ((N,Ad),(mu,N)):
        try:source.bounded_modulation(m,N,a)
        except ValueError:rejected+=1
        else:raise ArithmeticError('N-dependent original slow source accepted')
    hashes=dict(data['input_hashes'])
    for name in (source.NAME,Path(__file__).name):hashes[name]=source.sha(name)
    receipt=dict(all_passed=True,source_family=theorem['source_family'],
        exact_source_norm_identities=len(theorem['identities']),
        independent_original_sigma_derivative_comparisons=tested_sigma,
        independent_actual_cutoff_derivative_comparisons=tested_cutoff,
        independent_base_source_mixed_derivative_comparisons=tested_base,
        actual_amplitude_kept_factored_checks_normalized_by_amplitude=True,
        invalid_or_N_dependent_original_parameter_inputs_rejected=rejected,
        whole_support_original_cutoff_and_base_source_norm_formulas_available=True,
        numeric_actual_graph_parameter_enclosures_substituted=False,
        repair_completed_tensor_common_N_modified_cones_certified=False,input_hashes=hashes)
    (source.HERE/source.RECEIPT).write_bytes((json.dumps(receipt,indent=2)+'\n').encode())
    print('Actual whole-support source norms PASS:',len(theorem['identities']),'identities;',
        tested_sigma+tested_cutoff+tested_base,'independent derivative comparisons',flush=True)
    return receipt


if __name__=='__main__':run()
