"""Paper (9.16)/(9.4)/(9.5) shared symbolic bridge parameter.

K is the full physical C3 norm sum, not a guessed finite value. Its full
norm bounds and supplied-input tests remain explicit admission conditions.
A directed outer bound follows solely from K>=Cstar and cstar<=1/2.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
K_TERMS=['1e6','1/ell_c','1/Ra','Cstar','Pstar','A',
    'norm_C3(Fcore)','norm_C3(1/Fcore)','norm_C3_Z(P0)',
    'sum_five_norm_C3(core_moments)','norm_C3(Df)','norm_C3(Ef)','norm_C3(1/Df)']
class SharedExitParameters:
    def __init__(self,c):
        self.ctx=c
        name='lei_ren_part1_paper_logarithmic_core_majorant.json'
        self.major=json.loads((HERE/name).read_bytes())
        for n,d in self.major['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=d:raise ValueError('majorant source changed: '+n)
        kp_name='lei_ren_part1_paper_shared_pressure_Kp.json'
        self.kp=json.loads((HERE/kp_name).read_bytes())
        for n,d in self.kp['input_hashes'].items():
            if hashlib.sha256((HERE/n).read_bytes()).hexdigest()!=d:raise ValueError('pressure constant dependency changed: '+n)
        if self.kp['implicit_source_sha256']!=self.major['implicit_source_sha256'] or self.kp['datum_enclosure_sha256']!=self.major['datum_enclosure_sha256']:
            raise ValueError('pressure constant source mismatch')
        if not self.kp['same_source_preheat_pressure_Kp_certified'] or self.kp['pressure_units']!='P0/Pstar^2':
            raise ValueError('pressure constant units/admission missing')
        with mp.workdps(c.dps+40):
            self.logLambda=read_interval(c,self.major['logLambda'])
            self.Gbar=read_interval(c,self.major['Gupper_in_logC_definition'])
            if endpoints(self.Gbar)[0]<0:raise ValueError('logC floor requires nonnegative Gbar')
            if self.major['logC_definition']!='Lambda*Gupper + 2*logLambda + 1000':
                raise ValueError('symbolic Cstar definition changed')
            self.logC_lower=2*self.logLambda+1000
            self.log_h_upper=-100*self.logC_lower-c.ln(2)
            self.h_upper=c.exp(self.log_h_upper)
            self.h_box=c.mpf([0,endpoints(self.h_upper)[1]])
            if endpoints(self.h_upper)[1]<=0:raise ValueError('positive parameter envelope lost')
            self.definition=dict(K_terms=K_TERMS,
                A='10+norm_C3(log(Cstar Fcore))+norm_C3(Uzcore)',
                pressure_Kp=self.kp['pressure_Kp'],
                epsilon0='1/(1e6*(1+KN))',KN=self.kp['KN'],gamma='.01',
                eta_tol='min(epsilon0/4,e_star/100)',
                cstar='0.5*min(epsilon0*gamma^2/(1e6*K1),eta_tol/(12*K1),1)',
                h_b='cstar*K^-100',epsilon_b='h_b (the identical scalar)',
                K1='positive absolute moment constant >=1',
                core_norm_domain='0<=R/Ra<=exp(ell_c), full real Z domain',
                frozen_norm_domain='0<=log(R/Ra)<=log(110/Ra), full real Z domain',
                source_sha=self.major['implicit_source_sha256'],
                source_datum_sha=self.major['datum_enclosure_sha256'],
                source_logC_definition=self.major['logC_definition'])
            self.sha=hashlib.sha256(json.dumps(self.definition,sort_keys=True).encode()).hexdigest()
        self.input_hashes={**self.major['input_hashes'],**self.kp['input_hashes'],kp_name:hashlib.sha256((HERE/kp_name).read_bytes()).hexdigest(),name:hashlib.sha256((HERE/name).read_bytes()).hexdigest(),
            Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    def report(self):
        return dict(parameter_family_sha256=self.sha,definition=self.definition,input_hashes=self.input_hashes,
            logC_lower=self.logC_lower,log_h_upper=self.log_h_upper,h_upper=self.h_upper,
            shared_h_epsilon_box=self.h_box,
            h_equals_epsilon_by_symbolic_identity=True,cstar_restrictions_satisfied_by_symbolic_min=True,
            h_upper_proof='K>=Cstar; cstar<=1/2; logCstar>=2logLambda+1000',
            no_K_numeric_value_invented=True,actual_positive_h_not_materialized=True,
            h_box_zero_lower_is_enclosure_not_zero_selection=True,
            conditional_h_positive_requires_full_norm_K_finite=True,
            full_physical_C3_K_norms_certified=False,full_Df_Ef_input_tests_certified=False,
            preheat_pressure_Kp_certified=True,pressure_Kp=self.kp['pressure_Kp'],
            KN=self.kp['KN'],epsilon0=self.kp['epsilon0'],
            K1_e_star_numeric_bounds_certified=False,source_j_eta_tol_relation_verified=False,
            full_Section9_parameter_admission=False,
            scope='conditional shared-parameter family; core/connection calculations stay local in Z')
def run():
    from mpmath.ctx_iv import MPIntervalContext
    c=MPIntervalContext();c.dps=160;p=SharedExitParameters(c);r=p.report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(r),indent=2)+'\n',encoding='utf-8')
    print('Shared h_b=epsilon_b=cstar K^-100: symbolic identity restored; full K admission still pending')
    return r
if __name__=='__main__':run()
