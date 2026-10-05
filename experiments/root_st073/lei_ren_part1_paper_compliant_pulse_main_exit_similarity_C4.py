"""Original main/exit pulse: uncapped source selection and full similarity data.

All rows are bounds for the defining source functions. The interval selected
root, terminal datum and quadrature bounds are never chosen as point values.
The ordinary variable is y=log(R/Rp), xi=mu*y; both original charts remain.
"""
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_pulse_gap_similarity_C4 import (
    CompliantPulseGapSimilarityC4,decode_jet)
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import (
    SourceAST,source_precision,pulse_coefficients)
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import pulse_velocity_rows,ordinary_grid
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import gp_jets
from lei_ren_part1_paper_compliant_axial_pulse_field import gp,decay_integral
from lei_ren_part1_paper_compliant_axial_high_jets import positive_quadratic_jets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
DOMAIN=dict(main_xi=['.02','10'],exit_xi=['10','11'],Z=[-1,1],
    ordinary_derivative='d_y=d_logR=mu*d_xi',angular_normalization='Mtheta/(sqrt(2)*R^(3/2)*Utheta)')
FALSE_FLAGS=('pulse_main_exit_physical_decomposition_constructed','pulse_main_exit_cone_certified',
    'whole_outer_cone_certified','global_admissible_stress_lift_constructed',
    'independently_bounded_global_flat_remainder','physical_energy_integral_certified',
    'full_background_NS_validation','temporal_recursion','production_exact_point_parameters_selected')


def uncapped_selected_pulse_source(c,mu,K,incoming,energy_in,future_weighted,
                                  A,D,full_rows,end_weights,logG):
    """Same positive quadratic, evaluated with actual exp sources, never caps.

    Scalar inputs give a point source. Consistent Taylor inputs give all
    supplied axial orders by the exact implicit recurrence. For the extreme
    production parameters the callable remains a defining source; the
    separate directed bounds are used without materializing exp(logG).
    """
    det=A[0]*D[1]-A[1]*D[0]
    bounds=lambda value:endpoints(value) if hasattr(value,'_mpi_') else (value,value)
    def inverse(r1,r2):
        return [(r1*D[1]-r2*A[1])/det,(r2*A[0]-r1*D[0])/det]
    factors=[c.exp(-13*(c.mpf('.5')-i*mu)/mu-logG) for i in (1,2)]
    Q=[value*factor for value,factor in zip(incoming,factors)]
    u=inverse(Q[0]*(-mu),-(Q[1]-Q[0]))
    v=inverse(-mu*full_rows[0],-(full_rows[1]-full_rows[0]))
    logE=logG-c.ln(mu)
    nu=[mu*weight*c.exp(2*logE) for weight in end_weights]
    target=energy_in*(-mu)+future_weighted+(1-c.exp(-26))/4
    A2=K;A1=u[0]*0;A0=-target
    for weight,uj,vj in zip(nu,u,v):
        A2+=weight*vj*vj
        A1+=uj*(2*weight*vj)
        A0+=uj*uj*weight
    if hasattr(A1,'order'):
        disc=A1[0]**2-4*A2*A0[0]
        if bounds(A2)[0]<=0 or bounds(A0[0])[1]>=0 or bounds(disc)[0]<=0:
            raise ArithmeticError('Uncapped selected source positive branch lost')
        half=A1[0]/(2*A2)
        root=-half+c.sqrt(half*half-A0[0]/A2)
        ap,denominator=positive_quadratic_jets(c,root,A2,A1,A0)
    else:
        disc=A1**2-4*A2*A0
        if bounds(A2)[0]<=0 or bounds(A0)[1]>=0 or bounds(disc)[0]<=0:
            raise ArithmeticError('Uncapped selected source positive branch lost')
        half=A1/(2*A2)
        ap=-half+c.sqrt(half*half-A0/A2)
        denominator=2*A2*ap+A1
    if bounds(denominator)[0]<=0:
        raise ArithmeticError('Uncapped selected source derivative inverse lost')
    controls=[uj+ap*vj for uj,vj in zip(u,v)]
    return dict(ap=ap,controls=controls,logE=logE,row_factors=factors,nu=nu,
                affine_incoming=u,slopes=v,A2=A2,A1=A1,A0=A0,denominator=denominator)


def partial_linear_kernel(c,mu,xi,lam,cells=64,window=1000):
    """Directed full original integral_0^(xi/mu) exp(-lam*v) gp(xi-mu*v)dv."""
    t=xi/mu;length=min(endpoints(t)[0],mp.mpf(window))
    if length<0 or endpoints(lam)[0]<=0:raise ValueError('Positive original kernel required')
    total=c.mpf(0)
    for i in range(cells):
        a=c.mpf(length)*i/cells;b=c.mpf(length)*(i+1)/cells
        coord=xi-mu*c.mpf([endpoints(a)[0],endpoints(b)[1]])
        total+=(c.exp(-lam*a)-c.exp(-lam*b))/lam*gp(c,coord)['value']
    tail=c.mpf(0)
    if endpoints(t)[1]>length:
        tail=c.mpf([0,endpoints(11*c.exp(-lam*length)/lam)[1]])
    return dict(enclosure=total+tail,positive_omitted_tail_bound=tail,window=length,
                exact_definition='integral_0^(xi/mu) exp(-lambda_i*v)*gp(xi-mu*v) dv',
                omitted_history_not_deleted=True)


def partial_future_energy(c,xi,cells=64):
    """Same integral_xi^11 exp(-2a) gp(a)^2 da; affine cell endpoints retain xi."""
    length=11-xi
    if endpoints(length)[0]<0:raise ValueError('Original xi<=11 required')
    if endpoints(length)==(0,0):return c.mpf(0)
    total=c.mpf(0)
    for i in range(cells):
        a=xi*((c.mpf(cells)-i)/cells)+c.mpf(11)*i/cells
        b=xi*((c.mpf(cells)-i-1)/cells)+c.mpf(11)*(i+1)/cells
        value=gp(c,c.mpf([endpoints(a)[0],endpoints(b)[1]]))['value']
        total+=c.exp(-2*a)*decay_integral(c,2,length/cells)*value**2
    return total


def split_main_exit_stress(delta,mu,z,C,Xp,Bh,ml,mi,nl,ni,e0,J,Pbase,Pmemory):
    """Separate actual source histories; all meridional cross terms remain."""
    zero=C*0;zeros=[zero]*5
    base=pulse_coefficients(delta,mu,z,C,Xp,Bh,ml,nl,e0,zeros,Pbase)
    first=pulse_coefficients(delta,mu,z,C,Xp,Bh,mi,zeros,zeros,zeros,zeros)
    second=pulse_coefficients(delta,mu,z,C,Xp,zeros,zeros,ni,zeros,zeros,zeros)
    pressure=pulse_coefficients(delta,mu,z,C,Xp,zeros,zeros,zeros,zeros,zeros,Pmemory)
    loss=pulse_coefficients(delta,mu,z,C,Xp,zeros,zeros,zeros,zeros,J,zeros)
    result={label:{name:dict(part,extra_source='one') for name,part in parts.items()
                   if name!='selected_backward_energy_loss'} for label,parts in base.items()}
    result['theta']['incoming_Mz_transport']=dict(first['theta']['meridional_transport'],extra_source='incoming1')
    result['theta']['incoming_Mtheta_z_transport']=dict(second['theta']['meridional_transport'],extra_source='incoming2')
    for name in ('linear_axial_moment','nonlinear_meridional_transport'):
        result['axial']['incoming_Mz_'+name]=dict(first['axial'][name],extra_source='incoming1')
    result['axial']['same_absolute_pressure_memory']=dict(pressure['axial']['full_energy_and_pressure'],extra_source='Q')
    result['axial']['selected_end_energy_loss']=dict(loss['axial']['selected_backward_energy_loss'],extra_source='end_square')
    return result


def main_exit_shapes(c,mu,delta,z,C,Xp,ap,incoming,kernels,future_energy,future,J0,P0,xi):
    shape=gp_jets(c,xi)
    Bh=[ap*(shape[j]*math.factorial(j)*mu**j) for j in range(5)]
    ml=[ap*kernels[0]];nl=[ap*kernels[1]]
    l1=c.mpf('.5')-mu;l2=c.mpf('.5')-2*mu
    mi=[incoming[0]*(-l1)**j for j in range(5)]
    ni=[incoming[1]*(-l2)**j for j in range(5)]
    for j in range(4):
        ml.append(Bh[j]-ml[j]*l1)
        nl.append(Bh[j]-nl[j]*l2)
    distance=13-xi;K=c.exp(-2*distance);p=1+2*mu
    e0=[future*K-c.expm1(-2*distance)/(4*mu)-ap*ap*(c.exp(2*xi)*future_energy/mu)]
    square=product_rows(Bh,Bh)
    for j in range(4):e0.append(square[j]-C*0-(c.mpf('.5') if j==0 else 0)+e0[j]*(2*mu))
    loss=[J0*K*(2*mu)**j for j in range(5)]
    baseline=-C*C/(2*p);Ptilde=P0-baseline;zero=C*0
    Pbase=[baseline]+[zero]*4
    Pmemory=[Ptilde*p**j for j in range(5)]
    stress=split_main_exit_stress(delta,mu,z,C,Xp,Bh,ml,mi,nl,ni,e0,loss,Pbase,Pmemory)
    local_velocity=pulse_velocity_rows(c,delta,mu,z,C,Bh,ml)
    incoming_velocity=pulse_velocity_rows(c,delta,mu,z,C,[zero]*5,mi)
    return dict(Bh=Bh,ml=ml,mi=mi,nl=nl,ni=ni,e0=e0,J=loss,stress=stress,
                velocity_local=local_velocity,velocity_incoming=incoming_velocity,
                pressure_baseline_rows=[baseline*(-p)**j for j in range(5)],
                pressure_memory_rows=[Ptilde]+[zero]*4,
                original_gp_xi_derivative_Taylor=shape,actual_full_axial_shear_rows=shifted_rows(Bh,-(c.mpf('.5')+mu)))


def source_proof(records):
    asts=SourceAST();proofs={}
    def zero(name,value):
        if s.cancel(s.expand(s.expand_power_exp(value)))!=0:
            raise ArithmeticError('Main/exit source identity failed: '+name)
        proofs[name]=True
    asts.expression('axial_pulse_field','main','ap',wanted="selected['selected_ap_Taylor']")
    asts.expression('axial_pulse_field','main','B',wanted="ap*shape['value']")
    asts.expression('axial_pulse_field','main','By',wanted="ap*(self.mu*shape['derivative'])")
    asts.expression('axial_pulse_field','main','incoming',wanted='rows[row-1]*self.factor(-lam*t,cut)')
    asts.expression('axial_pulse_field','main','kernel',
        wanted="(c.exp(-lam*a)-c.exp(-lam*b))/lam*gp(c,xi-self.mu*v)['value']",augmented=True)
    asts.expression('pulse_mixed_C4','_high_packet','Brows',
        wanted='[ap*(shape[k]*math.factorial(k)*self.mu**k) for k in range(5)]')
    asts.expression('pulse_radial_C4','pressure_moment','p',
        wanted="IntervalTaylor(c,inlet['Mp_over_Pstar_squared'])+u*u*(kernel/2)")
    asts.expression('pulse_radial_C4','pressure_moment','rows',
        wanted="self.selection.future.angular.initial.datum.normalized_jets(c.mpf(Z),5)['normalized_pressure_coefficients']")
    for target in ('emu','end_energy'):
        asts.expression('axial_pulse_field','main',target)
    asts.expression('axial_pulse_field','main','e',count=2)
    for target in ('A2','A1','A0','logE','factors','nu','target'):
        asts.expression('pulse_main_exit_similarity_C4','uncapped_selected_pulse_source',target)
    mu,t,z,xi=s.symbols('mu t Z xi',positive=True);lam=s.symbols('lambda',positive=True)
    ap,ein,F,J,Kgp,IF=s.symbols('actual_ap incoming_energy complete_future end_loss Kgp remaining_gp_squared',real=True)
    g=s.Function('actual_gp')
    linear=s.exp(-lam*t)*s.Integral(s.exp(lam*s.Symbol('v'))*g(mu*s.Symbol('v')),(s.Symbol('v'),0,t))
    zero('actual_forward_linear_kernel_ODE',s.diff(linear,t)+lam*linear-g(mu*t))
    decay=s.exp(-2*(13-xi))
    energy=F*decay+(1-decay)/(4*mu)-ap**2*s.exp(2*xi)*IF/mu
    actual_dIF=-s.exp(-2*xi)*g(xi)**2
    zero('actual_selected_backward_main_energy_ODE',
         mu*(s.diff(energy,xi)+s.diff(energy,IF)*actual_dIF)-2*mu*energy-ap**2*g(xi)**2+s.Rational(1,2))
    forward=s.exp(2*xi)*(ein+ap**2*(Kgp-IF)/mu-(1-s.exp(-2*xi))/(4*mu))
    selected=(1-s.exp(-26))/4-mu*ein+mu*s.exp(-26)*(F-J)
    zero('actual_forward_backward_energy_same_selected_source',
         s.expand(mu*(forward-(energy-decay*J))).subs(ap**2*Kgp,selected))
    h1,h2,Q,D2=s.symbols('incoming_factor_1 incoming_factor_2 same_pressure_memory selected_end_square',real=True)
    delta,Xp=s.symbols('delta Xp',real=True);C=1/(1+z*z)
    row=lambda name:[s.Function(name+str(j))(z) for j in range(5)]
    Bh,ml,mi,nl,ni,e0,loss,Pm=[row(name) for name in ('Bh','ml','mi','nl','ni','e','J','Pmemory')]
    zeros=[s.Integer(0)]*5;Pb=[-C*C/(2*(1+2*mu))]+zeros[1:]
    env=dict(mp=SimpleNamespace(mpf=s.Rational),axial_derivative=lambda value:s.diff(value,z),
             product_rows=product_rows,shifted_rows=shifted_rows)
    asts.replay('pulse_end_stress_C3','pulse_coefficients',env)
    asts.replay('pulse_main_exit_similarity_C4','split_main_exit_stress',env)
    split=env['split_main_exit_stress'](delta,mu,z,C,Xp,Bh,ml,mi,nl,ni,e0,loss,Pb,Pm)
    native=env['pulse_coefficients'](delta,mu,z,C,Xp,Bh,[a+h1*b for a,b in zip(ml,mi)],
        [a+h2*b for a,b in zip(nl,ni)],e0,[D2*v for v in loss],[a+Q*b for a,b in zip(Pb,Pm)])
    factors=dict(one=1,incoming1=h1,incoming2=h2,Q=Q,end_square=D2)
    # D is one for the main pulse. The incoming exponential and pressure
    # rates have ALREADY been differentiated inside the source rows.
    R,B,H=s.symbols('R B H',positive=True)
    def scaled(part,j):
        rp,bp,dp,hp=map(s.Rational,part['mode'])
        return R**rp*B**bp*H**hp/s.sqrt(2)*part['full_derivative_rows'][j]
    for label in ('theta','axial'):
        for j in range(4):
            zero('actual_full_'+label+'_ordinary_stress_row'+str(j),
                 sum(scaled(part,j)*factors[part['extra_source']] for part in split[label].values())
                 -sum(scaled(part,j) for part in native[label].values()))
    l1=s.Rational(1,2)-mu;l2=s.Rational(1,2)-2*mu;p=1+2*mu;r=1-mu
    zero('actual_raw_Mz_incoming_radial_rate_zero',1-(s.Rational(1,2)+mu)-l1)
    zero('actual_raw_Mtheta_z_incoming_radial_rate_zero',s.Rational(3,2)-2*(s.Rational(1,2)+mu)-l2)
    zero('actual_raw_angular_memory_radial_rate_zero',s.Rational(3,2)-(s.Rational(1,2)+mu)-r)
    zero('actual_absolute_pressure_memory_radial_rate_zero',-2*(s.Rational(1,2)+mu)+p)
    zero('actual_selected_energy_loss_raw_radial_rate_zero',1-2*(s.Rational(1,2)+mu)+2*mu)
    zero('actual_main_exit_coordinate_xi10_shared',(xi/mu).subs(xi,10)-10/mu)
    zero('actual_xi11_gap_distance2',13-11-2)
    certificate=records['pulse_interface_certificate']['source_bound_functional_pulse_identities']
    for name in ('main_gap_linear_moments_whole_Z','main_gap_selected_energy_whole_Z',
                 'gap_end_original_pressure_whole_Z'):
        if not certificate[name]:raise ValueError('Actual selected terminal join unavailable: '+name)
        proofs['consumed_'+name]=True
    for n in range(6):
        name='functional_main_gap_axial_order'+str(n)
        if not certificate[name]:raise ValueError('Actual C5 terminal function unavailable')
        proofs['consumed_'+name]=True
    if not records['fifth_axial_jets']['actual_selected_ap_c1_c2_C5_available']:
        raise ValueError('Actual implicit C5 bounds unavailable')
    proofs['current_C5_selected_bounds_only_enclose_uncapped_source']=True
    proofs['ordinary_rows_use_mu_scaled_gp_derivatives_not_xi_derivatives']=True
    proofs['same_canonical_absolute_pressure_source_getter_retained']=True
    finite=s.symbols('same_finite_saddle_part',real=True)
    logG=-1/mu+finite+s.log(mu)
    for i in (1,2):
        zero('actual_uncapped_row'+str(i)+'_reduced_source_log',
             -13*(s.Rational(1,2)-i*mu)/mu-logG-(-s.Rational(11,2)/mu+13*i-finite-s.log(mu)))
    zero('actual_uncapped_end_scale_log',logG-s.log(mu)-(-1/mu+finite))
    Kj=s.symbols('exact_original_Rp_energy_weight',positive=True)
    zero('actual_uncapped_nu_original_Rp_Rv_units',mu*Kj*s.exp(2*(logG-s.log(mu)))-Kj*s.exp(2*logG)/mu)
    asts.expression('outer_pulse_map','__init__','self.logscale',
        wanted="self.rows['common_logpref']-c.ln(self.mu)")
    asts.expression('outer_pulse_map','correction_basis','det',wanted='A[0]*D[1]-A[1]*D[0]')
    asts.expression('outer_pulse_map','correction_basis','x',wanted='d-s')
    asts.expression('outer_pulse_map','correction_basis','base',wanted='c.exp(-x/2+mu*x)')
    asts.expression('outer_pulse_map','correction_basis','A[n]',wanted='mass*base',augmented=True)
    asts.expression('outer_pulse_map','correction_basis','D[n]',wanted='mass*base*quotient',augmented=True)
    asts.expression('outer_pulse_map','correction_basis','gram',
        wanted='(c.mpf(2)/cells)*beta**2*c.exp(-2*mu*s)/(ell*normalization**2)',augmented=True)
    asts.expression('outer_pulse_map','pulse_rows','core[row - 1]',
        wanted='(c.mpf(2)*band/panels)*base*c.exp(row*(2+u))',augmented=True)
    asts.expression('axial_amplitude_selection','__init__','logexact',
        wanted='self.future.params.log_mu+c.ln(Kj)+2*self.log_end_scale')
    if not records['outer_pulse_map']['bump_basis']['identity'].startswith('row2=row1+mu*D'):
        raise ValueError('Exact divided-difference integral basis missing')
    proofs['same_original_full_gp_beta_and_Gram_integral_definitions_retained']=True
    center,v=s.symbols('center original_beta_coordinate',real=True)
    x=-center-v
    zero('actual_first_row_original_beta_translation',s.exp(-x/2+mu*x)-s.exp((s.Rational(1,2)-mu)*(center+v)))
    zero('actual_divided_difference_original_beta_translation',
         s.exp(-x/2+mu*x)*(s.exp(mu*x)-1)/mu
         -(s.exp((s.Rational(1,2)-2*mu)*(center+v))-s.exp((s.Rational(1,2)-mu)*(center+v)))/mu)
    return dict(identities=proofs,input_hashes=asts.hashes,actual_source_AST_bindings=asts.bindings,
        arbitrary_axial_functions_used_not_chosen_enclosure_values=True,
        uncapped_positive_quadratic_source_callable_installed=True,
        production_exact_point_evaluation=False,source_caps_used_as_defining_field_values=False)


class CompliantPulseMainExitSimilarityC4(CompliantPulseGapSimilarityC4):
    @source_precision
    def __init__(self,cells=64):
        super().__init__(256)
        name=PREFIX+'pulse_gap_similarity_C4_check.json';raw=(HERE/name).read_bytes();gate=json.loads(raw)
        if not gate['all_passed'] or not gate['actual_original_whole_inactive_gap_similarity_companion_constructed']:
            raise ValueError('Current admitted gap source companion required')
        if (gate['actual_five_defect_family_sha256'],gate['implicit_source_sha256'])!=(self.family,self.source):
            raise ValueError('Main/exit and admitted gap source family differ')
        for path,digest in gate['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:
                raise ValueError('Admitted gap source changed: '+path)
        self.hashes.update(gate['input_hashes']);self.hashes[name]=hashlib.sha256(raw).hexdigest()
        self.cells=cells
        fifth=self.records['fifth_axial_jets']['whole_Z']['selected']
        self.ap=decode_jet(self.ctx,fifth['selected_ap_Taylor'])
        self.incoming=[decode_jet(self.ctx,v) for v in fifth['incoming']['moment_Taylor']]
        self.incoming_energy=decode_jet(self.ctx,fifth['incoming']['energy_Taylor'])
        self.mainproof=source_proof(self.records);self.hashes.update(self.mainproof['input_hashes'])

    @source_precision
    def main_exit(self,Z,xi):
        c=self.ctx;z0=c.mpf(Z);xi=c.mpf(xi)
        if endpoints(z0)[0]<-1 or endpoints(z0)[1]>1:raise ValueError('Original Z[-1,1] required')
        if endpoints(xi)[0]<endpoints(c.mpf('.02'))[0] or endpoints(xi)[1]>11:
            raise ValueError('Original main/exit xi[.02,11] required')
        z=IntervalTaylor.variable(c,z0,5)
        C=IntervalTaylor(c,[1+z0**2,2*z0,1,0,0,0]).reciprocal()
        kernel=[partial_linear_kernel(c,self.mu,xi,c.mpf('.5')-i*self.mu,self.cells) for i in (1,2)]
        remaining=partial_future_energy(c,xi,self.cells)
        rows=main_exit_shapes(c,self.mu,self.delta,z,C,self.Xp,self.ap,self.incoming,
                             [v['enclosure'] for v in kernel],remaining,self.future,self.J0,self.P0,xi)
        bp=c.mpf('.5')+self.mu;p=1+2*self.mu;r=1-self.mu
        logR=self.logRp+xi/self.mu
        logB=dict(logPstar=self.logP,actual_log_inlet_U=self.logU,inverse_mu=-xi/(2*self.mu),finite=-xi)
        logH=-r*xi/self.mu;logD0=-1/self.mu+self.finite
        extras=dict(one=c.mpf(0),incoming1=-(c.mpf('.5')-self.mu)*xi/self.mu,
                    incoming2=-(c.mpf('.5')-2*self.mu)*xi/self.mu,
                    Q=-p*(13-xi)/self.mu,end_square=2*logD0)
        def logs(rp,bpow,hpow=0,extra='one',norm=0):
            return dict(source_logR=rp*logR,**{key:bpow*value for key,value in logB.items()},
                        signed_original_memory_log=hpow*logH,exact_source_history_log=extras[extra],normalization=norm)
        stress={}
        for label,sectors in rows['stress'].items():
            stress[label]={}
            for name,part in sectors.items():
                rp,bpow,dp,hpow=part['mode']
                stress[label][name]=dict(original_mode=part['mode'],exact_extra_source=part['extra_source'],
                    exact_source_log_parts=logs(rp,bpow,hpow,part['extra_source'],-c.ln(2)/2),
                    full_stress_mixed3_coefficient_enclosures=ordinary_grid(part['full_derivative_rows'],3))
        def velocity(rp,bpow,ordinary,extra='one',norm=0):
            return dict(exact_source_log_parts=logs(rp,bpow,0,extra,norm),
                        full_velocity_mixed4_coefficient_enclosures=ordinary_grid(ordinary,4))
        velocities=dict(radial=dict(local=velocity(c.mpf('.5'),1,rows['velocity_local']['radial'],norm=-c.ln(2)/2),
                incoming_history=velocity(c.mpf('.5'),1,rows['velocity_incoming']['radial'],'incoming1',-c.ln(2)/2)),
            theta=dict(actual=velocity(0,1,rows['velocity_local']['theta'])),
            axial=dict(actual=velocity(0,1,rows['velocity_local']['axial'])))
        pressure=dict(radial_swirl_particular=dict(exact_source_log_parts=logs(0,2),
                full_pressure_mixed4_coefficient_enclosures=ordinary_grid(rows['pressure_baseline_rows'],4)),
            same_absolute_datum_memory=dict(exact_source_log_parts=dict(logPstar=2*self.logP,
                actual_log_inlet_U=2*self.logU,inverse_mu=-13/self.mu,finite=c.mpf(-26)),
                full_pressure_mixed4_coefficient_enclosures=ordinary_grid(rows['pressure_memory_rows'],4)))
        zero=C*0
        def moment(parts,ordinary):
            return dict(exact_source_log_parts=parts,
                        full_moment_mixed4_coefficient_enclosures=ordinary_grid(ordinary,4))
        raw_in1=dict(logRp=self.logRp,logPstar=self.logP,actual_log_inlet_U=self.logU)
        raw_in2=dict(logRp=c.mpf('1.5')*self.logRp,logPstar=2*self.logP,
                     actual_log_inlet_U=2*self.logU,normalization=c.ln(2)/2)
        angular_memory=dict(logRp=c.mpf('1.5')*self.logRp,logPstar=self.logP,
                            actual_log_inlet_U=self.logU,normalization=c.ln(2)/2)
        energy_constant=dict(logRp=self.logRp,logPstar=2*self.logP,actual_log_inlet_U=2*self.logU,
                             finite=c.mpf(-26),exact_selected_end_square_scale=2*logD0)
        pinpart=C*C*self.Pin;swirllimit=C*C/(2*p)
        moments=dict(theta=dict(equilibrium=moment(logs(c.mpf('1.5'),1,norm=c.ln(2)/2),
                [C/r*r**j for j in range(5)]),
                signed_original_history=moment(angular_memory,[C*(self.Xp-1/r)]+[zero]*4)),
            z=dict(local_pulse_integral=moment(logs(1,1),
                shifted_rows([C*v for v in rows['ml']],c.mpf('.5')-self.mu,4)),
                original_incoming=moment(raw_in1,[C*self.incoming[0]]+[zero]*4)),
            theta_z=dict(local_pulse_integral=moment(logs(c.mpf('1.5'),2,norm=c.ln(2)/2),
                shifted_rows([C*C*v for v in rows['nl']],c.mpf('.5')-2*self.mu,4)),
                original_incoming=moment(raw_in2,[C*C*self.incoming[1]]+[zero]*4)),
            z_theta=dict(full_unperturbed_energy=moment(logs(1,2),
                shifted_rows([C*C*v for v in rows['e0']],-2*self.mu,4)),
                selected_end_energy_loss=moment(energy_constant,[-C*C*self.J0]+[zero]*4)),
            p=dict(original_inlet=moment(dict(logPstar=2*self.logP),[pinpart]+[zero]*4),
                swirl_limit=moment(dict(logPstar=2*self.logP,actual_log_inlet_U=2*self.logU),[swirllimit]+[zero]*4),
                original_swirl_decay=moment(dict(logPstar=2*self.logP,actual_log_inlet_U=2*self.logU,
                    exact_forward_time_decay=-p*xi/self.mu),[-swirllimit*(-p)**j for j in range(5)])))
        return dict(Z=z0,xi=xi,domain=DOMAIN,full_meridional_stress_log_sectors=stress,
            full_velocity_log_sectors=velocities,full_absolute_pressure_log_sectors=pressure,
            five_raw_cumulative_moment_log_sectors=moments,main_exit_source_rows=rows,
            original_partial_linear_kernel_bounds=kernel,original_partial_future_energy_bound=remaining,
            exact_source_logs=dict(R=logR,B=logB,H=logH,D0=logD0,extra=extras),
            selected_source_is_uncapped_positive_quadratic=True,
            coefficients_are_only_enclosures_of_defining_sources=True,
            original_full_gp_and_axial_shear_retained=True,
            original_forward_linear_histories_not_replaced_by_end_only_histories=True,
            same_absolute_pressure_getter_and_raw_Mp_kept_separate=True,
            **{flag:False for flag in FALSE_FLAGS})

    @source_precision
    def report(self):
        c=self.ctx
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,domain=DOMAIN,
            whole_original_main=self.main_exit([-1,1],['.02','10']),
            whole_original_exit=self.main_exit([-1,1],['10','11']),
            original_main_left=self.main_exit([-1,1],'.02'),
            common_main_exit=self.main_exit([-1,1],10),
            original_exit_gap=self.main_exit([-1,1],11),
            source_function_proof=self.mainproof,current_selected_ap_C5_enclosure=self.ap,
            current_incoming_C5_enclosures=self.incoming,current_incoming_energy_C5_enclosure=self.incoming_energy,
            exact_selected_source_log_definitions=dict(
                common_logG=-1/self.mu+self.finite+c.ln(self.mu),
                end_scale=-1/self.mu+self.finite,
                incoming_rows=[-c.mpf('5.5')/self.mu+13*i-self.finite-c.ln(self.mu) for i in (1,2)],
                end_energy='nu_j=mu*Kj*exp(2logE); Kj is original Rp-weight exp(-26)*exp(-2mu*center)*Gram'),
            exact_selected_source_integral_definitions=dict(
                Kpulse='integral_0^11 exp(-2a)*gp(a)^2 da',
                full_rows='Pi_i=exp(-13lambda_i/mu-logG)*integral_0^(11/mu) exp(lambda_i*v)*gp(mu*v) dv',
                first_row='A_j=integral exp(lambda_1*(center_j+v))*beta(v) dv',
                divided_difference='D_j=integral exp(lambda_1*(center_j+v))*(exp(-mu*(center_j+v))-1)/mu*beta(v) dv',
                Gram='integral exp(-2mu*v)*beta(v)^2 dv; beta has original width .15 and exact normalization',
                incoming='same exact C5 inlet functions, not chosen coefficient box values',
                future='mu*exp(-26)*complete_future_energy/2 from original corrected flatten/angular/heat sources',
                pressure='same original analytic datum getter, transported by exact FTC from canonical Rv pressure'),
            same_canonical_terminal_pressure_enclosure=self.P0,
            actual_original_whole_main_exit_similarity_companion_constructed=True,
            main_exit_and_exit_gap_similarity_source_functional_joins_consumed=True,
            source_caps_used_as_defining_field_values=False,input_hashes=dict(self.hashes,
                **{Path(__file__).name:hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}),
            **{flag:False for flag in FALSE_FLAGS})


@source_precision
def run():
    result=CompliantPulseMainExitSimilarityC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original main/exit full moments/velocity/pressure/stress source companion generated; physical/cone pending',flush=True)
    return result


if __name__=='__main__':run()
