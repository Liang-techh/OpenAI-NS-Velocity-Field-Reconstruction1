"""Finite-N bounds with actual scalars and a factored original Rd amplitude.

No ancestor graph is reconstructed. This bounds the supported modulation
and its five defect densities before the independent moment repair. It
does not choose a sufficient common N or admit the completed stress cone.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_patch_physical_numeric_bounds as numeric
from lei_ren_part1_paper_compliant_current_O3_frequency_source_bounds import actual_source_norms
from lei_ren_part1_paper_compliant_current_O3_frequency_majorants import profile_majorants
from lei_ren_part1_paper_compliant_current_physical_tensor_locator_operator import (
    original_native_radius_recipes, SYMBOLS, V)

HERE, PREFIX, sha = numeric.HERE, numeric.PREFIX, numeric.sha
NAME = PREFIX+'current_O3_actual_parameter_majorants.json'
RECEIPT = PREFIX+'current_O3_actual_parameter_majorants_check.json'
NORMS_RECEIPT = PREFIX+'current_O3_frequency_source_bounds_check.json'
PRE_NAME = PREFIX+'pre_pulse_mixed_C4.json'
UNITS = {'M':(1,1,'0',0), 'I':(1,1,'.5',1),
         'J':(2,2,'.5',1), 'S':(2,2,'0',0), 'Cp':(2,2,'-1',0)}
OPEN = ('independent_repair_bounds_available', 'completed_signed_tensor_error_bounds_available',
        'common_N_modified_cones_certified', 'actual_coefficient_recursion_certified')


def encoded(value):
    return numeric.encode(numeric.pack(value))


def contains(outer, inner):
    a,b = numeric.transport.endpoints(outer); u,v = numeric.transport.endpoints(inner)
    return a <= u <= v <= b


def positive_integer_N(N):
    if type(N) is not int or N < 1:
        raise ValueError('One positive finite Python integer N required')
    return N


def maximum(c, a, b):
    if isinstance(a, s.Basic) or isinstance(b, s.Basic):
        return s.Max(a,b)
    al,ah = numeric.transport.endpoints(a); bl,bh = numeric.transport.endpoints(b)
    return c.mpf([max(al,bl),max(ah,bh)])


def normalized_source_norms(c, mu):
    B = [c.mpf(1),c.mpf(32),c.mpf(1792),
         256*c.exp(-c.mpf(1)/2)*(3/c.sqrt(2))**9,3104*c.exp(-c.mpf(2))*6**6]
    C = [B[0]]+[4**j*B[j] for j in range(1,5)]
    a = (1+mu)/2; b,d,e = (mu*B[j] for j in range(1,4))
    P = [c.mpf(1),a,a*a+b,a**3+3*a*b+d,a**4+6*a*a*b+3*b*b+4*a*d+e]
    T = [c.exp(c.mpf(1))]+[maximum(c,c.exp(c.mpf(1))/2**j,P[j]) for j in range(1,5)]
    return dict(B=B,C=C,T=T,Fhat=[[T[j]*math.factorial(k) for k in range(6)] for j in range(5)],L1=a)


def normalized_envelopes(c, mu, N, C, F, L):
    """Ordinary derivative caps; every F-linear result has Ad outside."""
    square = [sum(math.comb(j,k)*C[k]*C[j-k] for k in range(j+1)) for j in range(5)]
    H = [mu/(8*c.pi*N)*sum(math.comb(j,k)*square[j-k]*(4*c.pi*N)**k
        for k in range(j+1)) for j in range(5)]
    bell = [c.mpf(1)]
    for j in range(1,5):
        bell.append(sum(math.comb(j-1,k)*H[j-k]*bell[k] for k in range(j)))
    exp_cap = c.exp(H[0]); G = [exp_cap*b for b in bell]; dG = [exp_cap*H[0]]+G[1:]
    theta = [[sum(math.comb(j,i)*F[j-i][k]*G[i] for i in range(j+1)) for k in range(6)] for j in range(5)]
    inc = [[sum(math.comb(j,i)*F[j-i][k]*dG[i] for i in range(j+1)) for k in range(6)] for j in range(5)]
    axial = [[c.sqrt(mu)/(2*c.pi*N)*sum(math.comb(j,i)*math.comb(j-i,h)*F[j-i-h][k]*C[h]*(2*c.pi*N)**i
        for i in range(j+1) for h in range(j-i+1)) for k in range(6)] for j in range(5)]
    return dict(exponent_A_over_N_derivative_majorants=H,exponential_derivative_majorants=G,
        exponential_increment_majorants=dG,modified_theta_over_Pstar_over_Ad_majorants=theta,
        theta_increment_over_Pstar_over_Ad_majorants=inc,modified_axial_over_Pstar_over_Ad_majorants=axial,
        theta_shear_error_vs_same_periodic_loop_upper=mu*C[0]*C[1]/(2*c.pi*N),
        axial_shear_error_vs_same_periodic_loop_upper=exp_cap/N*(2*c.sqrt(mu)*C[0]*(mu*C[0]**2/(8*c.pi))
            +c.sqrt(mu)/c.pi*(C[0]*L+C[1])))


def product_caps(A, B):
    return [[sum(math.comb(j,i)*math.comb(k,h)*A[i][h]*B[j-i][k-h]
        for i in range(j+1) for h in range(k+1)) for k in range(6)] for j in range(5)]


def density_caps(F, bounds):
    u = bounds['modified_axial_over_Pstar_over_Ad_majorants']
    e = bounds['theta_increment_over_Pstar_over_Ad_majorants']
    t = bounds['modified_theta_over_Pstar_over_Ad_majorants']
    uu,Fe,ee = product_caps(u,u),product_caps(F,e),product_caps(e,e)
    positive_energy = [[Fe[j][k]+ee[j][k]/2 for k in range(6)] for j in range(5)]
    return dict(M=u,I=e,J=product_caps(u,t),
        S=[[uu[j][k]+positive_energy[j][k] for k in range(6)] for j in range(5)],Cp=positive_energy)


def inputs():
    source = numeric.inputs(); c = source['ctx']; hashes = dict(source['hashes'])
    receipt = json.loads((HERE/NORMS_RECEIPT).read_bytes())
    if not receipt['all_passed'] or not receipt['whole_support_original_cutoff_and_base_source_norm_formulas_available']:
        raise ValueError('Checked whole-support source norm formulas required')
    if receipt['source_family'] != source['accepted']['source_family']:
        raise ValueError('Finite-N norms belong to a different actual source family')
    for name,digest in receipt['input_hashes'].items():
        if name in hashes and hashes[name] != digest:
            raise ValueError('Conflicting norm/current scalar dependency: '+name)
        if sha(name) != digest: raise ValueError('Changed actual source norm input: '+name)
        hashes[name] = digest
    hashes[NORMS_RECEIPT] = sha(NORMS_RECEIPT)
    if hashes.get(PRE_NAME) != sha(PRE_NAME): raise ValueError('Current original endpoint packet not bound')
    pressure = json.loads((HERE/numeric.PRESSURE_NAME).read_bytes())['compliant_source']
    if pressure['implicit_source_definition']['c_mu'] != '.001': raise ValueError('Current defining c_mu changed')
    logmu = numeric.transport.read_interval(c,pressure['parameter_bounds']['log_mu'])
    mu = numeric.transport.read_interval(c,source['accepted']['records'][numeric.transport.COVER_NAME]
        ['actual_current_source_cover_hypotheses']['source_mu'])
    derived_logmu = c.ln(c.mpf('.001'))-4*source['logP']
    if not contains(logmu,derived_logmu) or not contains(mu,c.exp(derived_logmu)):
        raise ValueError('Actual mu/log_mu enclosures miss their original defining equation')
    ml,mh = numeric.transport.endpoints(mu)
    if ml <= 0 or mh > numeric.transport.endpoints(c.mpf('.001'))[1] or not all(mp.isfinite(v) for v in (ml,mh)):
        raise ValueError('Actual finite 0<mu<=.001 required')
    j1 = numeric.transport.read_interval(c,json.loads((HERE/PRE_NAME).read_bytes())['join_packets']['slope_exit']['original_J'])
    if j1._mpi_ != c.mpf('.5')._mpi_: raise ValueError('Exact original J(1)=1/2 required')
    logAd = -c.exp(c.mpf(40))/2-c.mpf(26)/5
    logC = numeric.transport.read_interval(c,json.loads((HERE/numeric.NORM_NAME).read_bytes())['selected_logCstar'])
    radius = dict(logRm=source['logRref']-6,logRd=source['logRref']+source['logP'],
        logRw=source['logRref']+source['logP']+1,logRp=source['logRref']+source['logP']+1-60*logmu)
    return dict(source=source,ctx=c,mu=mu,log_mu=logmu,exact_J_at1=j1,logAd=logAd,
        logC=logC,logP=source['logP'],radius_logs=radius,hashes=hashes)


def exact_theorem():
    c = SimpleNamespace(mpf=lambda v:s.Rational(str(v)),exp=s.exp,sqrt=s.sqrt,pi=s.pi)
    mu = s.Symbol('same_actual_mu',positive=True,finite=True)
    N = s.Symbol('one_finite_N',positive=True,integer=True,finite=True)
    Ad = s.Symbol('same_original_Rd_amplitude',positive=True,finite=True)
    norms = normalized_source_norms(c,mu); expected = actual_source_norms(mu,1); checks = {}
    def zero(name,a,b):
        difference=s.expand(a-b)
        if difference != 0 and s.cancel(difference) != 0:
            raise ArithmeticError('Actual scalar/factor identity failed: '+name)
        checks[name] = True
    for label,key in (('B','original_sigma_ordinary_majorants'),('C','actual_cutoff_ordinary_logR_majorants'),('T','original_base_time_factor_majorants')):
        for j,value in enumerate(norms[label]): zero('same_actual_'+label+str(j),value,expected[key][j])
    for j in range(5):
        for k in range(6): zero('same_actual_Fhat_%d_%d'%(j,k),norms['Fhat'][j][k],expected['original_theta_over_Pstar_logR_axial_majorants'][j][k])
    zero('same_actual_L1',norms['L1'],expected['original_log_theta_logR_majorant'])
    # Replay the actual interval recurrence with symbolic positive source
    # variables, compared against the already proved majorant program.
    C = s.symbols('C0:5',positive=True,finite=True)
    F = [[s.Symbol('F%d_%d'%(j,k),positive=True,finite=True) for k in range(6)] for j in range(5)]
    L = s.Symbol('L1',positive=True,finite=True)
    got = normalized_envelopes(c,mu,N,C,F,L)
    want = profile_majorants(mu,N,C,[[Ad*v for v in row] for row in F],L)
    for key,rows in got.items():
        other = key.replace('_over_Ad',''); target = want[other]
        divisor = Ad if '_over_Ad' in key else 1
        if isinstance(rows,list):
            for j,row in enumerate(rows):
                if isinstance(row,list):
                    for k,value in enumerate(row): zero('finite_N_'+key+'_%d_%d'%(j,k),value,target[j][k]/divisor)
                else: zero('finite_N_'+key+'_'+str(j),row,target[j]/divisor)
        else: zero('finite_N_'+key,rows,target)
    asts = numeric.transport.SourceAST()
    asts.expression('pre_pulse_mixed_C4','slope','(J, mass)',wanted='slope_masses(c,y,self.cells)')
    asts.expression('pre_pulse_mixed_C4','inlet','self.cache[key]',wanted='self.slope(Z,1)')
    asts.method('pre_pulse_mixed_C4','slope_masses')
    endpoint_name='lei_ren_part1_paper_interval_outer_slope_field.py'
    node=next(n for n in ast.walk(ast.parse((HERE/endpoint_name).read_text())) if isinstance(n,ast.FunctionDef) and n.name=='transition_integrals')
    endpoint=[n for n in ast.walk(node) if isinstance(n,ast.If) and ast.unparse(n.test)=='y == 1']
    if len(endpoint)!=1 or ast.unparse(endpoint[0].body[0]) != "j = c.mpf('.5')":
        raise ValueError('Exact symmetric original J endpoint source changed')
    x=s.Symbol('original_sigma_argument',real=True)
    odds=1/(1-x)**2-1/x**2
    zero('original_sigma_odds_reflection',odds.subs(x,1-x),-odds)
    w=s.Symbol('positive_sigma_odds_exponential',positive=True)
    zero('original_sigma_complement_symmetry',(1/w)/(1+1/w),1-w/(1+w))
    checks['actual_source_integral_symmetric_J_at1_assignment_bound']=True
    asts.hashes[endpoint_name]=sha(endpoint_name)
    density_node=asts.expression('current_O3_finite_frequency_profiles','profile','densities',wanted="{'M':dict(Pstar_power=1,logR_power=0,sqrt2_power=0,coefficient=u),'I':dict(Pstar_power=1,logR_power='.5',sqrt2_power=1,coefficient=e),'J':dict(Pstar_power=2,logR_power='.5',sqrt2_power=1,coefficient=u*(E+e)),'S':dict(Pstar_power=2,logR_power=0,sqrt2_power=0,coefficient=u*u-E*e-e*e/2),'Cp':dict(Pstar_power=2,logR_power=-1,sqrt2_power=0,coefficient=E*e+e*e/2)}")
    E,u,e=s.symbols('original_E own_axial increment',real=True)
    actual=asts.evaluate(density_node,dict(E=Ad*E,u=Ad*u,e=Ad*e))
    for label,(a,p,r,h) in UNITS.items():
        row=actual[label]
        if (row['Pstar_power'],s.Rational(str(row['logR_power'])),row['sqrt2_power'])!=(p,s.Rational(r),h):
            raise ValueError('Actual moment density physical units changed')
        zero('moment_density_Ad_degree_'+label,s.diff(row['coefficient'],Ad)*Ad,a*row['coefficient'])
    recipes,proof = original_native_radius_recipes(); q=SYMBOLS
    zero('same_absolute_Rd_radius_recipe',recipes['O3_slope_mu'],q['log_Rref']+q['log_P']+V)
    zero('same_absolute_Rw_power_recipe',recipes['O3_power'],q['log_Rref']+q['log_P']+1+q['Tw']*V)
    M,J = s.symbols('Md original_J_at1',real=True)
    zero('actual_original_logAd_exact_endpoint',
        (s.Rational(1,10)-s.Rational(3,5)*J-(s.exp(M)+10)/2).subs(J,s.Rational(1,2)),
        -s.exp(M)/2-s.Rational(26,5))
    for label,(a,p,r,h) in UNITS.items():
        r=s.Rational(r);lc,t=s.symbols('log_C same_logR_offset',real=True);lp=s.exp(M)+11
        zero('physical_factor_combined_before_enclosure_'+label,
            a*(-s.exp(M)/2-s.Rational(26,5))+p*lp+r*(s.log(110)+10*lc+11*lp+t)+s.Rational(h,2)*s.log(2),
            (-s.Rational(a,2)+p+11*r)*s.exp(M)-s.Rational(26*a,5)+11*p+121*r+10*r*lc+r*s.log(110)+r*t+s.Rational(h,2)*s.log(2))
    return dict(passed=True,identities=checks,source_bindings=asts.bindings,
        input_hashes={**asts.hashes,**proof['input_hashes'],Path(__file__).name:sha(Path(__file__).name)},
        original_Ad_kept_separate_from_N_corrected_profiles=True,
        full_signed_tensor_and_independent_repair_not_certified=True)


class CurrentO3ActualParameterMajorants:
    def __init__(self,require_checked=True):
        self.data=inputs();self.ctx=self.data['ctx'];self.theorem=exact_theorem()
        self.hashes={**self.data['hashes'],**self.theorem['input_hashes']}
        self.definition=hashlib.sha256(json.dumps(self.hashes,sort_keys=True).encode()).hexdigest()
        self.norms=normalized_source_norms(self.ctx,self.data['mu'])
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt['all_passed'] or receipt['actual_parameter_majorant_definition_sha256']!=self.definition:
                raise ValueError('Checked same-source actual parameter majorants required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest: raise ValueError('Changed checked actual parameter envelope: '+name)

    def bounds(self,N):
        positive_integer_N(N);c=self.ctx;d=self.data;n=self.norms
        bounds=normalized_envelopes(c,d['mu'],N,n['C'],n['Fhat'],n['L1'])
        densities=density_caps(n['Fhat'],bounds); physical={}
        t=c.mpf([-2,c.mpf(1)/2])
        for label,(a,p,r,h) in UNITS.items():
            r=c.mpf(r)
            # Combine correlated giant source logs algebraically first.
            factor=(-c.mpf(a)/2+p+11*r)*c.exp(40)-c.mpf(26*a)/5+11*p+121*r+10*r*d['logC']+r*c.ln(110)+r*t+c.mpf(h)/2*c.ln(2)
            rows=[[sum(math.comb(j,i)*abs(r)**(j-i)*densities[label][i][k] for i in range(j+1))
                for k in range(6)] for j in range(5)]
            physical[label]=dict(Ad_power=a,Pstar_power=p,logR_power=str(UNITS[label][2]),sqrt2_power=h,
                factored_log_scale=factor,ordinary_logR_axial_majorants_after_R_power_derivatives=rows,
                log_absolute_upper=[[c.mpf(numeric.transport.endpoints(factor+c.ln(value))[1]) for value in row] for row in rows])
        return dict(finite_integer_N=N,source_family=d['source']['accepted']['source_family'],
            actual_parameter_majorant_definition_sha256=self.definition,
            whole_support_logR_offset=['-2','1/2'],whole_axial_interval=['-1','1'],
            source_scalar_enclosures=dict(mu=d['mu'],log_mu=d['log_mu'],Md='40',logPstar=d['logP'],
                exact_J_at1=d['exact_J_at1'],original_log_amplitude_at_Rd=d['logAd'],radius_logs=d['radius_logs']),
            source_amplitude_factor_log=d['logAd'],shear_error_caps_use_normalized_loop_units=True,
            original_source_norms_normalized_by_Ad=n,finite_N_modulation_majorants=bounds,
            five_defect_density_majorants_normalized_by_Ad_power=densities,
            physical_five_defect_density_per_dX_majorants=physical,
            actual_mu_and_original_log_amplitude_enclosures_substituted=True,
            original_amplitude_and_giant_radius_exponentials_not_materialized=True,
            density_majorants_are_absolute_triangle_bounds_not_signed_errors=True,
            N_is_query_parameter_not_selected_common_admissible_frequency=True,
            **{key:False for key in OPEN})


def run():
    field=CurrentO3ActualParameterMajorants(require_checked=False)
    result=dict(actual_parameter_majorant_definition_sha256=field.definition,
        exact_actual_parameter_and_factor_theorem=field.theorem,
        actual_parameter_finite_N_modulation_and_density_bounds_available=True,
        examples={str(N):field.bounds(N) for N in (1,10**12)},
        example_N_not_a_sufficient_common_N_certificate=True,**{key:False for key in OPEN},input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(encoded(result),indent=2)+'\n').encode())
    print('Actual finite-N source bounds generated: mu interval; original logAd; five physical density factors',flush=True)
    return result


if __name__=='__main__':run()
