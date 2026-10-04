"""Original inactive gap: five moments, absolute pressure and full stress.

The complete gap 11<=xi<=13-4mu is kept. Zero axial input does not
remove radial velocity or incoming histories. Exact positive factors stay
as source logs; coefficients are enclosures, never defining field values.
"""
import ast
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_pulse_end_cone import current_sources
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import (
    SourceAST,source_precision,pulse_coefficients,backward_bump_weights)
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import pulse_velocity_rows,ordinary_grid
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
DOMAIN=dict(xi='[11,13-4mu]',distance='d=13-xi in[4mu,2]',end_s='[-2/mu,-4]',
    gap_main_xi='[11,12]',gap_end_s='[-1/mu,-4]',Z=[-1,1],
    ordinary_derivative='d_s=d_logR=-mu*d_d')
FALSE_FLAGS=('pulse_gap_physical_decomposition_constructed','pulse_gap_cone_certified',
    'whole_outer_cone_certified','global_admissible_stress_lift_constructed',
    'independently_bounded_global_flat_remainder','physical_energy_integral_certified',
    'full_background_NS_validation','temporal_recursion')


def sources():
    records,hashes,family=current_sources()
    for stem in ('pulse_end_cone_check','fifth_axial_jets','outer_pulse_map',
        'pulse_interface_certificate','power_inlet_C4','power_inlet_C4_check',
        'pulse_end_flatten_join','outer_initial','outer_buffer'):
        name=PREFIX+stem+'.json';raw=(HERE/name).read_bytes();record=json.loads(raw)
        if (record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])!=family:
            raise ValueError('Inactive gap source family differs: '+stem)
        if 'all_passed' in record and not record['all_passed']:
            raise ValueError('Unaccepted inactive gap source: '+stem)
        for source,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:
                raise ValueError('Inactive gap source changed: '+source)
        hashes.update(record['input_hashes']);hashes[name]=hashlib.sha256(raw).hexdigest();records[stem]=record
    certificate=records['pulse_interface_certificate']
    for flag in ('exact_uncapped_selected_sources_used',
        'exact_functional_main_gap_and_gap_end_identities_certified',
        'functional_axial_derivatives_through5_identified'):
        if not certificate[flag]:raise ValueError('Actual inactive gap source functions missing: '+flag)
    join=records['pulse_end_flatten_join']['source_endpoint_binding']['identities']
    for flag in ('actual_pulse_P0_getter_is_same_canonical_flatten_pressure_over_C0_squared',
        'consumed_same_absolute_pressure_identified_by_original_FTC_and_power_datum'):
        if not join[flag]:raise ValueError('Actual uncapped absolute pressure source missing: '+flag)
    inlet=records['power_inlet_C4_check']['exact_functional_production_and_join_identities']
    for flag in ('canonical_Xp_constant','canonical_Mp_constant',
        'single_sample_only_encloses_proved_Z_independent_Hp_Pin_constants'):
        if not inlet[flag]:raise ValueError('Canonical original inlet function not bound: '+flag)
    return records,hashes,family


def decode_jet(c,value):
    return IntervalTaylor(c,[read_interval(c,row) for row in value['coefficients']])


def gap_shapes(c,mu,delta,z,C,Xp,M1,M2,future,J0,P0,distance):
    """Full ordinary rows in the common uncapped source reference."""
    zero=C*0;p=1+2*mu;K=c.exp(-2*distance)
    Bh=[zero for _ in range(5)]
    m=[M1*(-(c.mpf('.5')-mu))**j for j in range(5)]
    n=[M2*(-(c.mpf('.5')-2*mu))**j for j in range(5)]
    e0=[future*K-c.expm1(-2*distance)/(4*mu)]
    e0 += [(future-1/(4*mu))*K*(2*mu)**j for j in range(1,5)]
    loss=[J0*K*(2*mu)**j for j in range(5)]
    baseline=-C*C/(2*p);Ptilde=P0-baseline
    regularP=[baseline]+[zero]*4
    memoryP=[Ptilde*p**j for j in range(5)]
    allzero=[zero]*5
    regular=pulse_coefficients(delta,mu,z,C,Xp,Bh,allzero,allzero,e0,loss,regularP)
    first=pulse_coefficients(delta,mu,z,C,Xp,Bh,m,allzero,allzero,allzero,allzero)
    second=pulse_coefficients(delta,mu,z,C,Xp,Bh,allzero,n,allzero,allzero,allzero)
    pressure=pulse_coefficients(delta,mu,z,C,Xp,Bh,allzero,allzero,allzero,allzero,memoryP)
    selected={
        'theta':{'equilibrium':(regular['theta']['equilibrium'],'D0',False),
            'signed_original_memory':(regular['theta']['signed_original_memory'],'D0',False),
            'meridional_Mz_history':(first['theta']['meridional_transport'],'D1',False),
            'mixed_angular_axial_history':(second['theta']['meridional_transport'],'D2',False),
            'radial_shear':(regular['theta']['radial_shear'],'D0',False)},
        'axial':{'full_unperturbed_energy_and_pressure':(regular['axial']['full_energy_and_pressure'],'D0',False),
            'same_absolute_pressure_memory':(pressure['axial']['full_energy_and_pressure'],'D0',True),
            'linear_axial_moment_history':(first['axial']['linear_axial_moment'],'D1',False),
            'selected_backward_energy_loss':(regular['axial']['selected_backward_energy_loss'],'D0',False)}}
    stress={label:{name:dict(part,selected_D_recipe=which,pressure_memory=pmem)
        for name,(part,which,pmem) in rows.items()} for label,rows in selected.items()}
    velocity=pulse_velocity_rows(c,delta,mu,z,C,Bh,m)
    return dict(stress=stress,velocity=velocity,m=m,n=n,e0=e0,J=loss,
        pressure_baseline_rows=[baseline*(-p)**j for j in range(5)],
        pressure_memory_rows=[Ptilde]+[zero]*4,
        pressure_memory_coefficient=Ptilde,K=K)


def source_proof(records):
    asts=SourceAST();proofs={}
    def zero(name,value):
        if s.simplify(s.expand_power_exp(value))!=0:
            raise ArithmeticError('Inactive gap source identity failed: '+name)
        proofs[name]=True
    mu,d,z,L,u0=s.symbols('mu d Z saddle_L saddle_u0',positive=True)
    r=1-mu;p=1+2*mu;C=1/(1+z*z)
    finite=-3*L+2*s.log(u0)-s.log(6)/2-2*s.log(mu)
    c=SimpleNamespace(mpf=lambda x:s.Rational(str(x)),ln=s.log,exp=s.exp,expm1=lambda v:s.exp(v)-1)
    logpref=asts.evaluate(asts.expression('outer_pulse_map','pulse_rows','logpref'),
        dict(k0=1/(2*mu),L=L,u0=u0,mu=mu,c=c))
    logscale=asts.evaluate(asts.expression('outer_pulse_map','__init__','self.logscale'),
        dict(self=SimpleNamespace(rows={'common_logpref':logpref},mu=mu),c=c))
    zero('actual_logE_equals_reduced_minus_inverse_mu_plus_finite',logscale-(-1/mu+finite))
    for name in ('gap','gap_from_end'):
        node=asts.expression('axial_pulse_field',name,'finite')
        actual=asts.evaluate(node,dict(L=L,u0=u0,c=c,self=SimpleNamespace(mu=mu)))
        zero('actual_'+name+'_same_finite_source_recipe',actual-finite)
    leading=asts.evaluate(asts.expression('axial_pulse_field','gap','leading'),
        dict(D=d,finite=finite,self=SimpleNamespace(mu=mu)))
    end_leading=asts.evaluate(asts.expression('axial_pulse_field','gap_from_end','leading'),
        dict(s=-d/mu,finite=finite,self=SimpleNamespace(mu=mu)))
    zero('actual_both_gap_chart_reduced_leading_sources_identical',end_leading-leading)
    for i in (1,2):
        native=leading-i*d
        zero('actual_reduced_D'+str(i)+'_moment_scale',native-(logscale+(s.Rational(1,2)-i*mu)*d/mu))
    asts.expression('axial_pulse_field','_gap','future',
        wanted="self.selection.future.future(Z)['complete_future_energy_Taylor']/2")
    asts.expression('axial_pulse_field','_gap','full_weights',
        wanted='backward_bump_weights(c,self.mu,normal,-4)')
    weightfn=asts.method('axial_pulse_field','backward_bump_weights')
    if ast.literal_eval(weightfn.args.defaults[-1])!=256:
        raise ValueError('Native full-beta weight quadrature default changed')
    proofs['same_native_gap_weight_cells_256_required']=True
    asts.expression('axial_pulse_field','_gap','zero',wanted='C[0]*0')
    gapfn=asts.method('axial_pulse_field','_gap')
    expected=ast.dump(ast.parse('Cj*(c.exp(lam*center)*w[row-1])',mode='eval').body)
    terms=[node for node in ast.walk(gapfn) if isinstance(node,ast.AugAssign)
        and isinstance(node.op,ast.Sub) and ast.unparse(node.target)=='total'
        and ast.dump(node.value)==expected]
    if len(terms)!=1:raise ValueError('Actual signed full-beta gap moment changed')
    asts.bindings['axial_pulse_field._gap.actual_negative_full_beta_sum']=True
    asts.expression('axial_pulse_field','end','end_energy',
        wanted='Cj*Cj*(c.exp(-2*self.mu*center)*w[2])',augmented=True)
    asts.expression('axial_pulse_field','__init__','self.full_end_energy_weights',
        wanted="[c.exp(-2*self.mu*center)*self.pulse.basis['energy_gram'] for center in (-3,-1)]")
    asts.expression('pulse_end_stress_C3','end','rawP0',
        wanted="packet['flatten_future_defect_zeroth_Taylor']['pressure_defect_rows']")
    asts.expression('pulse_end_stress_C3','end','P0',
        wanted='-(IntervalTaylor.constant(c,1,5)/(2*self.heat.prate)+rawP)/(scale*scale)')
    asts.expression('pulse_end_stress_C3','end','P',
        wanted='[P0*c.exp(prate*v)+C*C*(c.expm1(prate*v)/(2*prate))]')
    F=s.Function('same_complete_future')(z);J0=s.Function('same_full_beta_squared_mass')(z)
    P0=s.Function('same_canonical_absolute_pressure_at_Rv')(z)
    Xp=s.symbols('same_original_Xp',real=True)
    M1=s.Function('same_full_linear_beta_weight_1')(z);M2=s.Function('same_full_linear_beta_weight_2')(z)
    K=s.exp(-2*d);Q=s.exp(-p*d/mu);dy=lambda value:-mu*s.diff(value,d)
    e0=K*F+(1-K)/(4*mu);J=K*J0
    P=-C*C/(2*p)+Q*(P0+C*C/(2*p));H=s.exp(-r*(13-d)/mu)
    X=1/r+(Xp-1/r)*H
    for i,M in ((1,M1),(2,M2)):
        lam=s.Rational(1,2)-i*mu;actual=M*s.exp(lam*d/mu)
        for j in range(5):zero('actual_gap_m'+str(i)+'_ordinary_row'+str(j),
            (-mu)**j*s.diff(actual,d,j)-actual*(-lam)**j)
        zero('actual_gap_end_m'+str(i)+'_complete_beta_source_at_minus4',
            actual.subs(d,4*mu)-M*s.exp(4*lam))
    zero('actual_unperturbed_energy_ODE',dy(e0)-2*mu*e0+s.Rational(1,2))
    zero('actual_full_energy_loss_ODE',dy(J)-2*mu*J)
    zero('actual_absolute_pressure_ODE',dy(P)-p*P-C*C/2)
    zero('actual_signed_angular_history_ODE',dy(X)-(1-r*X))
    zero('actual_gap_end_unperturbed_energy_source_at_minus4',
        e0.subs(d,4*mu)-(s.exp(-8*mu)*F+(1-s.exp(-8*mu))/(4*mu)))
    zero('actual_gap_end_loss_source_at_minus4',J.subs(d,4*mu)-s.exp(-8*mu)*J0)
    zero('actual_gap_end_pressure_source_at_minus4',
        P.subs(d,4*mu)-(P0*s.exp(-4*p)+C*C*(s.exp(-4*p)-1)/(2*p)))
    mainD=asts.evaluate(asts.expression('axial_pulse_field','gap','D'),dict(xi=s.Integer(12)))
    endD=asts.evaluate(asts.expression('axial_pulse_field','gap_from_end','D'),
        dict(s=-1/mu,self=SimpleNamespace(mu=mu)))
    zero('actual_gap_main_distance_at_xi12',mainD-1)
    zero('actual_gap_end_distance_at_s_minus_inverse_mu',endD-mainD)
    B0,Rv,D0=s.symbols('actual_Bv actual_Rv actual_end_scale',positive=True)
    B=B0*s.exp(p*d/(2*mu));R=Rv*s.exp(-d/mu)
    D1=D0*s.exp((s.Rational(1,2)-mu)*d/mu)
    D2=D0*s.exp((s.Rational(1,2)-2*mu)*d/mu)
    zero('actual_raw_Mz_history_constant',dy(R*B*D1))
    zero('actual_raw_Mtheta_z_history_constant',dy(R**s.Rational(3,2)*B**2*D2))
    zero('actual_raw_angular_history_constant',dy(R**s.Rational(3,2)*B*H))
    zero('actual_absolute_pressure_memory_BsquaredQ_constant',dy(B**2*Q))
    # Replay the actual new splitter and the admitted full stress formula
    # with arbitrary axial functions. This identifies every ordinary row,
    # rather than using overlap of independently enclosed histories.
    ns=dict(math=math,mp=SimpleNamespace(mpf=s.Rational),
        axial_derivative=lambda v:s.diff(v,z),product_rows=product_rows,shifted_rows=shifted_rows)
    asts.replay('pulse_end_stress_C3','pulse_coefficients',ns)
    asts.replay('pulse_end_physical_C2','pulse_velocity_rows',ns)
    asts.replay('pulse_gap_similarity_C4','gap_shapes',ns)
    delta=s.symbols('delta',real=True)
    split=ns['gap_shapes'](c,mu,delta,z,C,Xp,M1,M2,F,J0,P0,d)
    zeros=[s.Integer(0)]*5
    native=ns['pulse_coefficients'](delta,mu,z,C,Xp,zeros,
        [M1*s.exp((s.Rational(1,2)-mu)*d/mu)*(-(s.Rational(1,2)-mu))**j for j in range(5)],
        [M2*s.exp((s.Rational(1,2)-2*mu)*d/mu)*(-(s.Rational(1,2)-2*mu))**j for j in range(5)],
        split['e0'],split['J'],[(-mu)**j*s.diff(P,d,j) for j in range(5)])
    scales={'D0':D0,'D1':D1,'D2':D2}
    def value(part,j,which):
        rp,bp,dp,hp=map(s.Rational,part['mode'])
        return R**rp*B**bp*which**dp*H**hp/s.sqrt(2)*part['full_derivative_rows'][j]
    for label in ('theta','axial'):
        for j in range(4):
            reconstructed=sum(value(part,j,scales[part['selected_D_recipe']])*(Q if part['pressure_memory'] else 1)
                for part in split['stress'][label].values())
            actual=sum(value(part,j,D0) for part in native[label].values())
            zero('actual_full_'+label+'_stress_row'+str(j)+'_split_source_identity',reconstructed-actual)
    logC,logP,Tw=s.symbols('same_logC same_logP same_Tw',real=True)
    unit=SimpleNamespace(logC=logC,logP=logP,params=SimpleNamespace(Tw=Tw))
    unit.logRref=asts.evaluate(asts.expression('global_physical_assembly','__init__','self.logRref'),dict(self=unit,c=c))
    originalRp=asts.evaluate(asts.expression('global_physical_assembly','__init__','self.logRp'),dict(self=unit))
    zero('actual_direct_production_logRp_recipe',originalRp-(s.log(110)+10*(logC+logP)+logP+1+Tw))
    # Actual cumulative pressure has its original inlet and the same datum.
    Pin,U,datum=s.symbols('same_Pin same_U same_analytic_datum',real=True)
    t=(13-d)/mu;decay=s.exp(-p*t)
    mp_raw=C*C*(Pin+U*U*(1-decay)/(2*p))
    node=asts.expression('pulse_radial_C4','pressure_moment','p',
        wanted="IntervalTaylor(c,inlet['Mp_over_Pstar_squared'])+u*u*(kernel/2)")
    actual_mp=asts.evaluate(node,dict(IntervalTaylor=lambda ctx,value:value,c=c,
        inlet={'Mp_over_Pstar_squared':C*C*Pin},u=C*U,kernel=(1-decay)/p))
    zero('actual_cumulative_Mp_retains_original_Pin_and_decay',actual_mp-mp_raw)
    zero('actual_cumulative_Mp_FTC',dy(mp_raw)-C*C*U*U*decay/2)
    proofs['same_uncapped_pressure_datum_route_consumed']=True
    proofs['same_functional_main_gap_and_gap_end_source_certificate_consumed']=True
    return dict(identities=proofs,input_hashes=asts.hashes,actual_source_AST_bindings=asts.bindings,
        arbitrary_axial_functions_used_not_chosen_enclosure_values=True,
        original_full_beta_supports_inside_minus4_to0=True,
        ordinary_rows_through4_imply_same_axial_mixed_rows=True,
        exact_positive_exponentials_kept_as_defining_sources=True,
        caps_and_baseline_intervals_not_used_as_exact_boundary_values=True)


class CompliantPulseGapSimilarityC4:
    @source_precision
    def __init__(self,cells=256):
        if cells!=256:raise ValueError('The same native full-beta gap source requires cells=256')
        self.records,self.hashes,family=sources();self.family,self.source=family
        self.ctx=c=MPIntervalContext();c.dps=240
        selected=self.records['pulse_mixed_C4']
        self.mu=read_interval(c,selected['selected_mu']);self.delta=read_interval(c,selected['selected_delta'])
        self.normal=read_interval(c,selected['bump_normalization'])
        fifth=self.records['fifth_axial_jets']['whole_Z']
        self.controls=[decode_jet(c,value) for value in fifth['selected']['selected_scaled_end_coefficient_Taylor']]
        self.future=decode_jet(c,fifth['energy']['complete_future_energy_Taylor'])/2
        self.weights=backward_bump_weights(c,self.mu,self.normal,-4,cells)
        basis=self.records['outer_pulse_map']['bump_basis']
        gram=read_interval(c,basis['energy_gram'])
        zero=self.controls[0]*0;self.M=[zero,zero];self.J0=zero
        for Cj,center in zip(self.controls,(-3,-1)):
            for i in (1,2):
                self.M[i-1]-=Cj*(c.exp((c.mpf('.5')-i*self.mu)*center)*self.weights[i-1])
            self.J0+=Cj*Cj*(c.exp(-2*self.mu*center)*gram)
        self.U=read_interval(c,fifth['selected']['incoming']['Z_independent_constant_definitions']['U'])
        sample=self.records['outer_buffer']['samples'][-1]
        if not sample['pulse_inlet'] or endpoints(read_interval(c,sample['Z']))!=(mp.mpf('.5'),mp.mpf('.5')):
            raise ValueError('Original canonical inlet sample changed')
        canonical_q=c.mpf('1.25')
        self.Hp=read_interval(c,sample['Mtheta_over_sqrt2_R_3half_Pstar'][0])*canonical_q
        self.Xp=self.Hp/self.U
        self.Pin=read_interval(c,sample['Mp_over_Pstar_squared'][0])*canonical_q**2
        whole=self.records['pulse_end_stress_C3']['whole_original_end']
        # The whole end includes s=0. This is a conservative bound for the
        # SAME defining P0 getter, never an exact value chosen from that box.
        self.P0=decode_jet(c,whole['signed_pressure_relative_pure_swirl_reference_rows'][0])
        self.rows=self.records['outer_pulse_map']['pulse_rows']
        L=read_interval(c,self.rows['saddle_L']);u0=read_interval(c,self.rows['saddle_u0'])
        self.finite=-3*L+2*c.ln(u0)-c.ln(6)/2-2*c.ln(self.mu)
        source_eq=whole['full_meridional_stress_log_sectors']['theta']['equilibrium']['exact_source_log_parts']
        self.logP=read_interval(c,source_eq['logPstar'])
        self.logU=read_interval(c,source_eq['actual_log_inlet_U'])
        logC=read_interval(c,self.records['outer_initial']['selected_logCstar'])
        Tw=read_interval(c,fifth['selected']['incoming']['Z_independent_constant_definitions']['Tw'])
        self.logRp=c.ln(110)+10*(logC+self.logP)+self.logP+1+Tw
        self.proof=source_proof(self.records);self.hashes.update(self.proof['input_hashes'])
        self.cells=cells

    @source_precision
    def packet(self,Z,distance,whole=False,right_endpoint=False):
        c=self.ctx;z0=c.mpf(Z)
        if endpoints(z0)[0]<-1 or endpoints(z0)[1]>1:raise ValueError('Original Z[-1,1] required')
        d=(c.mpf([endpoints(4*self.mu)[0],2]) if whole else
            4*self.mu if right_endpoint else c.mpf(distance))
        if not whole and not right_endpoint and (
            endpoints(d)[0]<endpoints(4*self.mu)[1] or endpoints(d)[1]>2):
            raise ValueError('Original d=13-xi in[4mu,2] required')
        z=IntervalTaylor.variable(c,z0,5)
        C=IntervalTaylor(c,[1+z0**2,2*z0,1,0,0,0]).reciprocal()
        rows=gap_shapes(c,self.mu,self.delta,z,C,self.Xp,self.M[0],self.M[1],
            self.future,self.J0,self.P0,d)
        # Keep every defining scale. Common inverse-mu terms are grouped
        # symbolically before numerical enclosure.
        p=1+2*self.mu;r=1-self.mu
        logBparts=dict(logPstar=self.logP,actual_log_inlet_U=self.logU,
            inverse_mu=(-13+d)/(2*self.mu),finite=-13+d)
        logR=self.logRp+(13-d)/self.mu
        logH=-r*(13-d)/self.mu
        logD0=-1/self.mu+self.finite
        logD={0:logD0,1:(d/2-1)/self.mu+self.finite-d,
            2:(d/2-1)/self.mu+self.finite-2*d}
        logQ=-p*d/self.mu
        stress={}
        for label,sectors in rows['stress'].items():
            stress[label]={}
            for name,part in sectors.items():
                rp,bp,dp,hp=part['mode'];which=int(part['selected_D_recipe'][-1])
                parts=dict(source_logR=rp*logR,**{key:bp*value for key,value in logBparts.items()},
                    exact_selected_moment_scale=dp*logD[which],signed_original_memory_log=hp*logH,
                    normalization=-c.ln(2)/2)
                if part['pressure_memory']:parts['same_absolute_pressure_memory']=logQ
                stress[label][name]=dict(original_mode=part['mode'],selected_D_recipe=part['selected_D_recipe'],
                    exact_source_log_parts=parts,full_stress_mixed3_coefficient_enclosures={
                        key.replace('y','s',1):value for key,value in ordinary_grid(part['full_derivative_rows'],3).items()})
        pressure={
            'radial_swirl_particular':dict(exact_source_log_parts={key:2*value for key,value in logBparts.items()},
                full_pressure_mixed4_coefficient_enclosures=ordinary_grid(rows['pressure_baseline_rows'],4)),
            'same_absolute_datum_memory':dict(exact_source_log_parts={
                'logPstar':2*self.logP,'actual_log_inlet_U':2*self.logU,
                'inverse_mu':-13/self.mu,'finite':c.mpf(-26)},
                full_pressure_mixed4_coefficient_enclosures=ordinary_grid(rows['pressure_memory_rows'],4))}
        # The unchanged cumulative pressure primitive, separate from the datum.
        zero=C*0;Pinpart=C*C*self.Pin;tail=-C*C/(2*p)
        mp_rows={'same_original_inlet':dict(exact_source_log_parts={'logPstar':2*self.logP},
            full_moment_mixed4_coefficient_enclosures=ordinary_grid([Pinpart]+[zero]*4,4)),
            'same_original_swirl_limit':dict(exact_source_log_parts={
                'logPstar':2*self.logP,'actual_log_inlet_U':2*self.logU},
                full_moment_mixed4_coefficient_enclosures=ordinary_grid([-tail]+[zero]*4,4)),
            'retained_forward_swirl_decay':dict(exact_source_log_parts={
                'logPstar':2*self.logP,'actual_log_inlet_U':2*self.logU,
                'exact_forward_time_decay':-p*(13-d)/self.mu},
                full_moment_mixed4_coefficient_enclosures=ordinary_grid([tail*(-p)**j for j in range(5)],4))}
        velocity={}
        velocity_modes={'radial':(.5,1,1,0),'theta':(0,1,0,0),'axial':(0,1,1,0)}
        for label,ordinary in rows['velocity'].items():
            rp,bp,dp,hp=velocity_modes[label]
            logs=dict(source_logR=rp*logR,**{key:bp*value for key,value in logBparts.items()},
                exact_selected_moment_scale=dp*logD[1],
                normalization=-c.ln(2)/2 if label=='radial' else c.mpf(0))
            velocity[label]=dict(exact_source_log_parts=logs,
                full_velocity_mixed4_coefficient_enclosures=ordinary_grid(ordinary,4))
        # Raw cumulative linear/energy moments have exact cancellations in d.
        constant_m1=dict(logRp=self.logRp,logPstar=self.logP,actual_log_inlet_U=self.logU,
            inverse_mu=c.mpf('5.5')/self.mu,finite=self.finite-13)
        constant_m2=dict(logRp=c.mpf('1.5')*self.logRp,logPstar=2*self.logP,
            actual_log_inlet_U=2*self.logU,inverse_mu=c.mpf('5.5')/self.mu,
            finite=self.finite-26,normalization=c.ln(2)/2)
        energy_base=dict(logRp=self.logRp,logPstar=2*self.logP,actual_log_inlet_U=2*self.logU,finite=c.mpf(-26))
        angular_logs=dict(source_logR=c.mpf('1.5')*logR,**logBparts,normalization=c.ln(2)/2)
        angular_memory=dict(logRp=c.mpf('1.5')*self.logRp,logPstar=self.logP,actual_log_inlet_U=self.logU,normalization=c.ln(2)/2)
        baseline_energy=C*C*(c.expm1(2*d)/(4*self.mu))
        energy_derivatives=[baseline_energy]+[
            C*C*(-c.exp(2*d)/2)*(-2*self.mu)**(j-1) for j in range(1,5)]
        def moment(logs,ordinary):return dict(exact_source_log_parts=logs,
            full_moment_mixed4_coefficient_enclosures=ordinary_grid(ordinary,4))
        moments=dict(theta=dict(equilibrium=moment(angular_logs,[C/r*r**j for j in range(5)]),
                signed_history=moment(angular_memory,[C*(self.Xp-1/r)]+[zero]*4)),
            z=dict(full_nonzero_history=moment(constant_m1,[C*self.M[0]]+[zero]*4)),
            theta_z=dict(full_nonzero_history=moment(constant_m2,[C*C*self.M[1]]+[zero]*4)),
            z_theta=dict(complete_future=moment(energy_base,[C*C*self.future]+[zero]*4),
                gap_radial_integral=moment(energy_base,energy_derivatives),
                selected_end_loss=moment(dict(energy_base,exact_selected_end_square_scale=2*logD0),
                    [-C*C*self.J0]+[zero]*4)),
            p=mp_rows)
        return dict(Z=z0,distance=d,domain=DOMAIN,
            full_meridional_stress_log_sectors=stress,full_absolute_pressure_log_sectors=pressure,
            full_velocity_log_sectors=velocity,five_raw_cumulative_moment_log_sectors=moments,
            exact_source_logs=dict(R=logR,B=logBparts,H=logH,D0=logD0,D1=logD[1],D2=logD[2],Q=logQ),
            gap_source_rows=rows,
            all_five_cumulative_moments_and_same_absolute_datum_retained=True,
            zero_axial_input_does_not_zero_radial_history=True,
            original_analytic_P0_source_getter_retained=True,
            P0_whole_end_box_only_bounds_the_exact_terminal_function=True,
            original_G1_G2_positive_factors_not_materialized=True,
            numerical_caps_used_as_field_values=False,
            **{flag:False for flag in FALSE_FLAGS})

    @source_precision
    def report(self):
        result=dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            domain=DOMAIN,source_function_proof=self.proof,
            whole_original_gap=self.packet([-1,1],None,whole=True),
            right_end_source_packet=self.packet([-1,1],None,right_endpoint=True),
            same_chart_switch_packet=self.packet([-1,1],1),
            full_beta_linear_weights=self.weights,ordinary_C5_controls=self.controls,
            current_complete_future=self.future,current_full_beta_loss=self.J0,
            same_absolute_P0_function_bound=self.P0,input_hashes=self.hashes,
            actual_original_whole_inactive_gap_similarity_companion_constructed=True,
            gap_end_similarity_velocity4_moment4_stress3_pressure4_functional_join_verified=True,
            main_gap_and_gap_end_original_coordinate_interfaces_consumed=True,
            source_caps_used_as_defining_field_values=False,
            **{flag:False for flag in FALSE_FLAGS})
        result['input_hashes'][Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        return result


@source_precision
def run():
    result=CompliantPulseGapSimilarityC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Whole original inactive-gap moments/velocity/stress/pressure generated; physical/cone/global/recursion pending',flush=True)
    return result


if __name__=='__main__':run()
