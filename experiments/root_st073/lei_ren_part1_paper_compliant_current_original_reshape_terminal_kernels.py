"""Original full finite long-reshape endpoint kernels, with nonzero errors.

The exact endpoint integral is 1/k plus a bounded remainder, not a reset to
1/k. Source-owned B ordinary jets and the original frozen T are retained.
This completes three defining kernels, not the unresolved bridge/switch
integrals or inlet histories which feed their affine moment transport.
"""
import ast
import copy
import json
import math
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
from mpmath.ctx_iv import MPIntervalContext
import lei_ren_part1_paper_compliant_current_general_conditioned_slow_Z as general
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_inner_bridge_profiles import logarithm,square
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval

base=general.base;prior=general.prior;conditioned=general.conditioned
HERE,PREFIX,sha=general.HERE,general.PREFIX,general.sha;ep=general.ep
NAME=PREFIX+'current_original_reshape_terminal_kernels.json'
RECEIPT=PREFIX+'current_original_reshape_terminal_kernels_check.json'
GATE='actual_original_full_finite_terminal_reshape_kernel_functions_with_factored_errors_installed'
SOURCE=PREFIX+'actual_long_reshape_mixed_C4.json'
SOURCE_CHECK=PREFIX+'actual_long_reshape_mixed_C4_check.json'
PARAMS=PREFIX+'physical_norm_family.json'
CORE=PREFIX+'core_transfer.json'
STEP='lei_ren_part1_paper_shared_fixed_step_bound.json'
KINDS=dict(theta=('1.6',1,'1.55'),pressure=('.2',2,'.1'),swirl=('1.2',2,'1.1'))
OPEN=('actual_bridge_and_first_switch_integral_functions_installed',
      'actual_R110_source_histories_closed','actual_active_patch_genuine_coefficients_installed',
      'actual_complete_patch_integral_installed','numerical_original_source_point_or_integral_oracle_installed',
      'actual_five_controls_installed','current_whole_N_selected',
      *base.point.source.inertial.profiles.loop.OPEN)


def bind(hashes,name,digest):
    if sha(name)!=digest or name in hashes and hashes[name]!=digest:
        raise ValueError('Original reshape source dependency differs: '+name)
    hashes[name]=digest


def source_bindings():
    from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
    bindings=assignment_source_bindings('long_reshape_profiles','evaluate',{
        'y':'self.T*phase','inp':'self.inputs(Z)',
        'logu':'B*(1-sig)-logq+(y/10-self.core.logC-self.core.logP)',
        'logq':'logarithm(1+square(z))',
        'kernels':"{name:backward_kernel(c,B,self.T,y,m,low,high,self.proofs) for name,m,low,high in (('theta',1,'1.55','1.65'),('pressure',2,'.1','.3'),('swirl',2,'1.1','1.3'))}"})
    bindings.update(assignment_source_bindings('long_reshape_profiles','evaluate',{
        'decays':"{name:normalized_amplitude_decay(c,B*(sig*m)-y*k,low,high,y,self.proofs) for name,m,k,low,high in (('theta',1,'1.6','1.55','1.65'),('pressure',2,'.2','.1','.3'),('swirl',2,'1.2','1.1','1.3'))}"}))
    bindings.update(assignment_source_bindings('actual_long_reshape_mixed_C4','__init__',{
        'self.A':"read_interval(c,self.core.records['physical_norm_family']['A_upper'])",'self.T':'400*self.A'}))
    bindings.update(assignment_source_bindings('long_reshape_profiles','inputs',{
        'inlet':"self.switch.inlet(Z)",'original':"jet(inlet['actual_R110_log_shape_axial5_coefficients'])"}))
    path=HERE/(PREFIX+'long_reshape_profiles.py');tree=ast.parse(path.read_text(encoding='utf8'))
    fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='backward_kernel')
    doc=ast.get_docstring(fn)
    if 'exp(-k*t + m*B*(sigma(y/T)-sigma((y-t)/T)))' not in doc:
        raise ValueError('Original full finite defining kernel changed')
    evaluate=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='evaluate')
    outputs=[kw.value for node in ast.walk(evaluate) if isinstance(node,ast.Call)
        for kw in node.keywords if kw.arg=='log_Utheta_over_Pstar_axial5_coefficients']
    if len(outputs)!=1 or ast.dump(outputs[0])!=ast.dump(ast.parse('list(logu.coefficients)',mode='eval').body):
        raise ValueError('Original source phase0 cached logu rows no longer expose the defining function')
    bindings.update(assignment_source_bindings('long_reshape_mixed_C4','evaluate',{
        'parent':'self.reshape.evaluate(Z,phase=phase)','T':'self.reshape.T'}))
    B,q,C,P=sy.symbols('B logq logC logP');phase0logu=B-q-C-P
    if sy.expand(phase0logu+q+C+P-B)!=0:raise ArithmeticError('Original phase0 B recovery identity failed')
    return dict(passed=True,original_source_assignment_bindings=bindings,
        endpoint='y=T; sigma(1)=1; 1-sigma(1-t/T)=sigma(t/T)',
        source_B='same actual R110 normalization function, recovered from original phase0 logu',
        exact_phase0_B_recovery_identity='(B-log(1+Z^2)-logCstar-logPstar)+log(1+Z^2)+logCstar+logPstar=B',
        same_defining_B_function_identity_not_interval_overlap=True,
        ordinary_Z_jet_coefficients_at_every_point_not_a_single_center_Taylor_model=True)


def exact_identities():
    x=sy.symbols('x',positive=True);odds=1/(1-x)**2-1/x**2
    assert sy.simplify(odds.subs(x,1-x)+odds)==0
    L,r=sy.symbols('L r',positive=True)
    moments=[]
    for p in range(6):
        F=sy.exp(-r*L)*sum(sy.binomial(p,j)*L**(p-j)*sy.factorial(j)/r**(j+1) for j in range(p+1))
        assert sy.simplify(sy.diff(F,L)+L**p*sy.exp(-r*L))==0
        moments.append(str(F))
    z,d=sy.symbols('z d');b=sy.symbols('b1:6');rows=[sy.Integer(1)]
    for n in range(1,6):rows.append(sy.expand(sum(j*b[j-1]*d*rows[n-j] for j in range(1,n+1))/n))
    series=sy.series(sy.exp(d*sum(b[j-1]*z**j for j in range(1,6))),z,0,6).removeO()
    assert all(sy.expand(series.coeff(z,n)-rows[n])==0 for n in range(6))
    return dict(passed=True,original_sigma_odds_reflection_identity=True,
        complete_tail_moments_through5= moments,
        ordinary_Taylor_exponential_Bell_coefficients_through5=True,
        full_finite_decomposition='K=1/k+body_difference+tail_[L,T]-exp(-k*L)/k',
        sigma_body_upper='sigma(t/T)<=exp(4-(T/L)^2) for 0<=t<=L<=T/2',
        sigma_global_upper='sigma(t/T)<=8*t/T by sigma(0)=0 and original derivative bound8')


class TerminalKernelEnclosures:
    """Local source-jet enclosure algorithm; caller supplies valid B/T."""
    def __init__(self,c,B,T):
        if not isinstance(B,IntervalTaylor) or B.ctx is not c or B.order!=5:
            raise ValueError('Same-context ordinary source B Taylor jets0..5 required')
        self.c=c;self.B=B;self.T=c.mpf(T)
        if ep(self.T)[0]<=0 or any(not mp.isfinite(x) for x in ep(self.T)):
            raise ValueError('Original frozen T must be finite and positive')
        if any(not mp.isfinite(x) for row in B.coefficients for x in ep(row)):
            raise ValueError('Finite ordinary source B coefficients required')
        self.ledger=dict(directed_small_exponential_tails=0,positive_function_denominator_intersections=0,
            positive_function_root_intersections=0,directed_independent_log_rescalings=0)
        self.bases=(c.ln(self.T),c.mpf(0),c.mpf(0),c.mpf(0),c.mpf(0))
        self.scalar=lambda value:prior.ScaledEnclosure(prior.FormalScale(self.bases),value,self.ledger)
        self.Tsource=prior.ScaledEnclosure(prior.FormalScale(self.bases,(1,0,0,0,0)),1,self.ledger)
        self.bounds=[self.scalar(c.mpf(max(abs(v) for v in ep(row)))) for row in B.coefficients]

    def exp_source(self,logvalue):
        return prior.ScaledEnclosure(prior.FormalScale(self.bases,offset=logvalue),1,self.ledger)

    def bell_bounds(self,m):
        polynomials=[{0:self.scalar(1)}]
        for n in range(1,6):
            row={}
            for j in range(1,n+1):
                factor=self.bounds[j]*(self.c.mpf(m*j)/n)
                for p,previous in polynomials[n-j].items():
                    row[p+1]=row.get(p+1,self.scalar(0))+factor*previous
            polynomials.append(row)
        return polynomials

    def powers(self,value,count=5):
        rows=[self.scalar(1)]
        for _ in range(count):rows.append(rows[-1]*value)
        return rows

    def evaluate(self,kind,*,L=4096):
        if kind not in KINDS:raise ValueError('An original theta/pressure/swirl kernel is required')
        if type(L) is not int or not 1<=L<=10**6:raise ValueError('Explicit finite integration split1..10^6 required')
        c=self.c;k0,m,r0=KINDS[kind];k=c.mpf(k0);r=c.mpf(r0);ell=c.mpf(L)
        if ep(self.T)[0]<2*L:raise ValueError('Body split must satisfy L<=T/2 throughout source box')
        beta0=max(abs(v) for v in ep(self.B[0]))
        ratio=8*m*c.mpf(beta0)/self.T
        if ep(ratio)[1]>ep(k-r)[0]:raise ValueError('Original full integrand positive decay rate is not source-certified')
        logeps=4-(self.T/ell)**2;epsilon=self.exp_source(logeps)
        smallness=self.bounds[0]*epsilon*m
        if not smallness.zero and ep(smallness.scale.evaluate()+c.ln(c.mpf(ep(smallness.coefficient)[1])))[1]>0:
            raise ValueError('Body exponential perturbation must be at most1')
        polys=self.bell_bounds(m);epspowers=self.powers(epsilon)
        invT=self.scalar(8).positive_divide(self.Tsource,self.bases[0]);invTpowers=self.powers(invT)
        tail_exp=self.exp_source(-r*ell);baseline_tail=self.exp_source(-k*ell)*(1/k)
        body=[];tails=[];errors=[];coefficients=[]
        for n in range(6):
            if n==0:
                b=self.bounds[0]*epsilon*(c.exp(1)*m/k)
                tail=tail_exp*(1/r)+baseline_tail
            else:
                b=self.scalar(0);tail=self.scalar(0)
                for p,value in polys[n].items():
                    if p<1:raise ArithmeticError('Positive-order Bell polynomial must vanish at delta=0')
                    b+=value*epspowers[p]*(c.exp(1)/k)
                    mass=sum((c.mpf(math.comb(p,j))*ell**(p-j)*math.factorial(j)/r**(j+1) for j in range(p+1)),c.mpf(0))
                    tail+=value*invTpowers[p]*tail_exp*mass
            error=b+tail
            symmetric=prior.ScaledEnclosure(error.scale,error.coefficient*c.mpf([-1,1]),self.ledger)
            coefficients.append(symmetric+(1/k if n==0 else 0))
            body.append(b);tails.append(tail);errors.append(error)
        decay=self.exp_source(-k*self.T+self.B[0]*m)
        decay_bounds=[]
        for n,row in enumerate(polys):
            bound=sum(row.values(),self.scalar(0))*decay
            decay_bounds.append(bound if n==0 else prior.ScaledEnclosure(bound.scale,bound.coefficient*c.mpf([-1,1]),self.ledger))
        return dict(kind=kind,coefficients=coefficients,absolute_errors=errors,body_errors=body,tail_errors=tails,
            positive_incoming_decay_Z_bounds=decay_bounds,
            proof=dict(full_original_finite_T_integral_enclosed=True,source_T_not_shortened=True,
                finite_body_split=L,body_log_epsilon=logeps,body_smallness=smallness.record(),
                source_decay_guard='8*m*sup|B0|/T <= k-rate_min',source_decay_ratio=ratio,
                k=k,m=m,rate_min=r,original_frozen_T_Z_exact_zero=True,
                ordinary_Z_Taylor_coefficients_not_derivatives=True,
                original_nonzero_B_Z_jet_error_sectors_retained=True,
                complete_original_tail_not_zeroed=True,incoming_decay_not_reset_to_zero=True))


def record(row):
    return {key:([v.record() for v in value] if key in ('coefficients','absolute_errors','body_errors','tail_errors','positive_incoming_decay_Z_bounds') else value)
        for key,value in row.items()}


class OriginalLongReshapeTerminalKernels:
    mode='original_full_finite_terminal_reshape_kernel_functions_partial'
    def __init__(self,dps=500):
        if type(dps) is not int or dps<500:raise ValueError('At least500 digits preserve original source cache bits')
        self.c=c=MPIntervalContext();c.dps=dps;self.hashes={}
        checked=json.loads((HERE/SOURCE_CHECK).read_bytes())
        if not checked.get('all_passed') or not checked.get('current_actual_long_reshape_mixed4_available'):
            raise ValueError('Accepted current actual long-reshape source required')
        for name,digest in checked['input_hashes'].items():bind(self.hashes,name,digest)
        bind(self.hashes,SOURCE_CHECK,sha(SOURCE_CHECK))
        source=json.loads((HERE/SOURCE).read_bytes());bind(self.hashes,SOURCE,sha(SOURCE))
        self.family=checked['actual_five_defect_family_sha256'];self.source=checked['implicit_source_sha256'];self.datum=checked['datum_enclosure_sha256']
        for key,want in (('actual_five_defect_family_sha256',self.family),('implicit_source_sha256',self.source),('datum_enclosure_sha256',self.datum)):
            if source[key]!=want:raise ValueError('Current same-family reshape/datum required')
        if source['shared_exact_axial_source']['formal_integrals_numerically_reconstructed']:
            raise ValueError('This successor preserves unresolved original upstream integrals')
        self.source_binding=source_bindings();self.old_source=source
        self.inlet=source['actual_R110_inlet']['actual_inherited_axial5_packet']
        if self.inlet['source_full_kernel']!='integral_0^y exp(-k*t+m*B*(sigma(y/T)-sigma((y-t)/T)))dt; (k,m)=(1.6,1),(.2,2),(1.2,2)':
            raise ValueError('Original cache full finite kernel definition differs')
        if ep(read_interval(c,self.inlet['phase_enclosure']))!=(0,0):raise ValueError('Original actual R110 phase0 cache required')
        self.Z=read_interval(c,self.inlet['Z'])
        if ep(self.Z)!=(-1,1):raise ValueError('Accepted original whole-Z source range required')
        params=json.loads((HERE/PARAMS).read_bytes());core=json.loads((HERE/CORE).read_bytes())
        for name in (PARAMS,CORE,STEP):bind(self.hashes,name,sha(name))
        A=read_interval(c,params['A_upper']);T=read_interval(c,source['actual_Rsh_exit']['original_T'])
        expected=A*400
        if not ep(T)[0]<=ep(expected)[0]<=ep(expected)[1]<=ep(T)[1]:
            raise ValueError('Original cached T must contain the directed same-source400*Abar')
        if source['actual_Rsh_exit']['source_T']!='400*selected_Abar; all y derivatives keep original inverse-T factors':
            raise ValueError('Original frozen T source binding changed')
        for packet in (source['whole_reshape'],source['actual_R110_inlet'],source['actual_Rsh_exit'],*source['interior_packets']):
            if read_interval(c,packet['original_T'])._mpi_!=T._mpi_:
                raise ValueError('The same original frozen T is required at every Z/radius source packet')
        self.source_binding.update(original_cached_T_contains_directed_same_source400_Abar=True,
            fixed_T_constructor_AST_bound=True,all_original_source_packets_share_identical_T=True)
        logu=IntervalTaylor(c,[read_interval(c,row) for row in self.inlet['log_Utheta_over_Pstar_axial5_coefficients']])
        z=IntervalTaylor.variable(c,self.Z,5)
        B=logu+logarithm(1+square(z))+read_interval(c,params['selected_logCstar'])+read_interval(c,core['logPstar'])
        rows=[]
        for n,value in enumerate(B.coefficients):
            if n<=2:
                bound=2*A/math.factorial(n);lo,hi=ep(value);upper=ep(bound)[1]
                if max(lo,-upper)>min(hi,upper):raise ValueError('Recovered actual B conflicts with original same-source C2 theorem')
                value=c.mpf((max(lo,-upper),min(hi,upper)))
            rows.append(value)
        self.B=IntervalTaylor(c,rows);self.T=T;self.backend=TerminalKernelEnclosures(c,self.B,T)
        step=json.loads((HERE/STEP).read_bytes())
        if not step.get('global_derivative_bound_certified') or step.get('global_derivative_upper')!=8:
            raise ValueError('Original fixed sigma derivative8 certificate required')
        for module in (general,base,prior):bind(self.hashes,Path(module.__file__).name,sha(Path(module.__file__).name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))

    def evaluate(self,L=4096):
        with mp.workdps(self.c.dps+40):
            kernels={kind:self.backend.evaluate(kind,L=L) for kind in KINDS}
        return dict(source_family=self.family,implicit_source_sha256=self.source,datum_enclosure_sha256=self.datum,
            mode=self.mode,original_Z_domain=self.Z,original_frozen_T=self.T,
            same_actual_R110_B_ordinary_Z_coefficients=list(self.B.coefficients),
            source_function_binding=self.source_binding,
            actual_full_finite_endpoint_kernels={key:record(value) for key,value in kernels.items()},
            source_owned_affine_moment_transport=dict(
                theta='K_theta + D_theta * unknown_original_theta_R110',
                theta_z='same_V110 * K_theta + D_theta * unknown_original_theta_z_R110',
                pressure='K_pressure + D_pressure * unknown_original_pressure_R110',
                swirl='K_swirl + D_swirl * unknown_original_swirl_R110',
                mean='same_V110 + exp(-T)*(unknown_original_mean_R110-same_V110)',
                axial='same_V110^2 + exp(-T)*(unknown_original_axial_R110-same_V110^2)'),
            separate_P0_retained=True,original_shared_V110_E_function_namespace=self.old_source['current_actual_source_namespace'],
            upstream_inlet_history_functions_still_unresolved=True,
            original_bridge_switch_formal_reconstruction_flag_preserved_false=True,
            no_original_ancestor_constructors_or_producers_executed=True)


def run():
    began=time.monotonic();owner=OriginalLongReshapeTerminalKernels()
    report=dict(**{GATE:True},source_family=owner.family,actual_original_terminal_kernel_evaluation=owner.evaluate(),
        exact_source_kernel_identities=exact_identities(),**dict.fromkeys(OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Three original full finite long-reshape terminal kernels and ordinary-Z0..5 enclosures on accepted whole axis. Nonzero factored remainders and inherited-history multipliers remain. Unresolved actual bridge/switch and R110 inlet functions are not promoted to numerical source closure.')
    (HERE/NAME).write_text(json.dumps(base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Original full finite reshape endpoint kernels computed with factored source errors',flush=True)
    return report


if __name__=='__main__':run()
