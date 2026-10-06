"""Source-bound finite-frequency O2/O3 profile candidate; not a completed field.

The actual profiles cross the original zero-shear seam. Their local logR/Z
jets and five moment increment densities are callable. Integrated moments,
common-axis pressure, a uniform admissible N and independent repair remain
uninstalled, so no modified radial velocity or full NS field is exported.
"""
import json
import math
import operator
import ast
from pathlib import Path
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_transition_direction import (
    CurrentO3TransitionDirection,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,_verify_hashes)
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import raw_pre_velocity_rows
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_current_O3_transition_direction_operator import SourceAST,upper

NAME=PREFIX+'current_O3_finite_frequency_profiles.json'
RECEIPT=PREFIX+'current_O3_finite_frequency_profiles_check.json'
OPEN=('current_O3_transition_finite_N_modified_profiles_installed',
      'current_O3_transition_modified_five_moment_repair_certified',
      'current_whole_O3_transition_signed_two_vector_cone_certified',
      'current_O2_buffer_modified_whole_cone_certified',
      'current_modified_radial_velocity_pressure_and_complete_tensor_installed')

def cutoff_rows(c,t):
    """Ordinary logR derivatives, including both exact flat edges."""
    left=sigma_jets(c,t+2);right=sigma_jets(c,4*t-1)
    a=[left[j]*math.factorial(j) for j in range(5)]
    b=[1-right[0]]+[-right[j]*4**j*math.factorial(j) for j in range(1,5)]
    return product_rows(a,b)

def trig_rows(c,phase,N,harmonic,cosine=False):
    angle=harmonic*c.pi*phase;rate=harmonic*c.pi*N
    return [rate**j*(c.cos(angle+j*c.pi/2) if cosine else c.sin(angle+j*c.pi/2)) for j in range(5)]

def exact_profile_theorem():
    v,mu,N,phi=s.symbols('offset mu integer_N phase',real=True)
    E=s.Function('actual_E0')(v);chi=s.Function('same_chi')(v);a0=1-2*s.diff(E,v)/E
    A=-mu*chi**2*s.sin(4*s.pi*phi)/(8*s.pi)
    B=-E*s.sqrt(mu)*chi*s.cos(2*s.pi*phi)/(2*s.pi)
    da=mu*chi**2*s.cos(4*s.pi*phi);bL=2*s.sqrt(mu)*chi*s.sin(2*s.pi*phi)
    checks={}
    def zero(name,a,b):
        if s.simplify(s.expand(a-b))!=0:raise ArithmeticError('Finite-frequency identity failed: '+name)
        checks[name]=True
    zero('supported_A_periodic_primitive',s.diff(A,phi),-da/2)
    zero('supported_B_periodic_primitive',s.diff(B,phi),E*bL/2)
    zero('supported_A_zero_mean',s.integrate(A,(phi,0,1)),0)
    zero('supported_B_zero_mean',s.integrate(B,(phi,0,1)),0)
    AP=A.subs(phi,N*v);BP=B.subs(phi,N*v)
    EN=E*s.exp(AP/N);UN=BP/N
    aN=a0+da.subs(phi,N*v)-2*s.diff(A,v).subs(phi,N*v)/N
    bN=s.exp(-AP/N)*(bL.subs(phi,N*v)+2*s.diff(B,v).subs(phi,N*v)/(N*E))
    zero('actual_finite_N_theta_shear_chain_rule',1-2*s.diff(EN,v)/EN,aN)
    zero('actual_finite_N_axial_shear_chain_rule',2*s.diff(UN,v)/EN,bN)
    zero('actual_slow_B_derivative_retains_E0_derivative',s.diff(B,v)/E,
         -s.sqrt(mu)*(chi*s.diff(E,v)/E+s.diff(chi,v))*s.cos(2*s.pi*phi)/(2*s.pi))
    oldE,de,du,x=s.symbols('E0 actual_swirl_increment actual_axial_increment X',real=True)
    newE=oldE+de;newU=du
    moments={'M':(newU,du),'I':(s.sqrt(2*x)*(newE-oldE),s.sqrt(2*x)*de),
      'J':(s.sqrt(2*x)*newU*newE,s.sqrt(2*x)*du*(oldE+de)),
      'S':(newU**2-newE**2/2+oldE**2/2,du**2-oldE*de-de**2/2),
      'Cp':((newE**2-oldE**2)/(2*x),(oldE*de+de**2/2)/x)}
    for name,(a,b) in moments.items():zero('actual_five_moment_increment_'+name,a,b)
    asts=SourceAST()
    asts.expression('pre_pulse_mixed_C4','axial','B',wanted='[c.mpf(0)]*5')
    packet_calls=[node for node in ast.walk(asts.method('pre_pulse_mixed_C4','axial'))
        if isinstance(node,ast.Call) and ast.unparse(node.func)=='self.packet']
    wanted=ast.parse("[zero-c.mpf('.5')]+[zero]*3",mode='eval').body
    if len(packet_calls)!=1 or ast.dump(packet_calls[0].args[5])!=ast.dump(wanted):
        raise ValueError('Actual O2 buffer constant log-shear source changed')
    checks['actual_O2_buffer_original_log_shear_minus_half']=True
    asts.expression('pre_pulse_mixed_C4','slope_mu','logU',
      wanted="[zero-c.mpf('.5')-mu*sig[0]]+[zero-mu*sig[k]*math.factorial(k) for k in range(1,4)]")
    # Bind the phase coordinate to the actual physical radius program,
    # rather than treating a chart selector as a logR derivative.
    radius=asts.replay('global_physical_assembly','radius',dict(PULSE=()))
    rr,md,w=s.symbols('logRref Md buffer_offset',real=True)
    rd=rr+s.exp(md)+11
    owner=SimpleNamespace(ctx=SimpleNamespace(mpf=s.sympify),logRref=rr,logP=s.exp(md)+11)
    left=radius(owner,'O2_buffer',w,dict(actual_y=s.exp(md)+w,exact_radius_source='original'),None)[0]
    right=radius(owner,'O3_slope_mu',v,{},None)[0]
    zero('actual_O2_modulation_coordinate_is_logR_minus_logRd',left-rd,w-11)
    zero('actual_O3_modulation_coordinate_is_logR_minus_logRd',right-rd,v)
    zero('actual_shared_phase_and_all_derivatives_at_old_seam',left.subs(w,11),right.subs(v,0))
    sy=s.Symbol('phase_sin_squared',real=True);cc=s.Symbol('same_cutoff',real=True)
    zero('supported_loop_shear_uniform_lower_expression',
       2*mu*s.Function('same_sigma')(v)+mu*cc**2*(1-2*sy)+4*mu*cc**2*sy/(2+3*mu),
       2*mu*s.Function('same_sigma')(v)+mu*cc**2-6*mu**2*cc**2*sy/(2+3*mu))
    return dict(identities=checks,original_eq4_38_phase='N*logR, translated by the constant -N*logRd',
      actual_shared_coordinate='O2_buffer_offset-11 = log(R/Rd); O3_offset = log(R/Rd)',
      cutoff='sigma(t+2)*(1-sigma(4*t-1))',compact_support=[-2,'1/2'],plateau=[-1,'1/4'],
      exact_flat_original_profile_edges=['O2_buffer9','O3_offset1/2'],
      original_internal_seam_now_modulated=True,
      normalized_theta='E_N/Pstar=(E0/Pstar)*exp(A/N)',
      normalized_axial='U_N/Pstar=-(E0/Pstar)*sqrt(mu)*chi*cos(2*pi*N*t)/(2*pi*N)',
      positive_profiles='E_N>0 for every finite integer N>=1, by E0>0 and the exponential',
      cumulative_five_moments_and_pressure_not_copied_from_original=True,
      uniform_N_and_pressure_moment_repair_and_cone_not_certified=True,
      input_hashes=asts.hashes,passed=True)

class CurrentO3FiniteFrequencyProfiles:
    @source_precision
    def __init__(self,direction):
        if type(direction) is not CurrentO3TransitionDirection or not direction.acceptance_loaded:
            raise ValueError('Same checked current O3 variable direction source required')
        direction.assert_graph();self.direction=direction;self.registry=direction.registry;self.ctx=direction.ctx
        self.family=direction.family;self.source=direction.source;self.datum_sha=direction.datum_sha
        self.pre=self.registry.owners['o3'].physical.pre
        if self.registry.owners['o2'].physical.pre is not self.pre:raise ValueError('Same O2/O3 pre source required')
        self.mu=self.ctx.mpf(self.registry.owners['o3'].pulse.mu);self.theorem=exact_profile_theorem()
        self.hashes={**direction.hashes,**self.theorem['input_hashes'],
            Path(__file__).name:sha(Path(__file__).name)}
        # The translated phase is globally smooth across the original seam.
        c=self.ctx;q=sigma_jets(c,c.mpf(1)/4)[0]
        floor=2*self.mu*q
        if endpoints(q)[0]<=0 or endpoints(self.mu/2-floor)[0]<=0:raise ArithmeticError('Supported loop floor unresolved')
        self.supported_loop_vs_excess_lower=floor

    @source_precision
    def profile(self,region,Z,coordinate,N):
        """Finite-N candidate coefficients with exact common Pstar factored."""
        self.direction.assert_graph()
        if isinstance(N,bool):raise ValueError('Finite positive integer frequency required')
        try:N=operator.index(N)
        except TypeError as exc:raise ValueError('Finite positive integer frequency required') from exc
        if N<1:raise ValueError('Finite positive integer frequency required')
        c=self.ctx;v=c.mpf(coordinate);z=c.mpf(Z);lo,hi=endpoints(v)
        if endpoints(z)[0]<-1 or endpoints(z)[1]>1:raise ValueError('Current Z[-1,1] required')
        if region=='O2_buffer' and 0<=lo<=hi<=11:
            source=self.pre.axial(z,buffer_offset=v);t=v-11
        elif region=='O3_slope_mu' and 0<=lo<=hi<=1:
            source=self.pre.slope_mu(z,v);t=v
        else:raise ValueError('Original O2 buffer or O3 transition chart required')
        old=raw_pre_velocity_rows(c,source)
        if any(any(endpoints(value)!=(0,0) for value in row.coefficients) for row in old['axial']):
            raise ValueError('This supported modulation requires the exact original zero axial source')
        chi=cutoff_rows(c,t);n=c.mpf(N);phase=n*t;mu=self.mu
        sin4=trig_rows(c,phase,n,4);cos2=trig_rows(c,phase,n,2,cosine=True)
        square=product_rows(chi,chi)
        Arows=[-mu*r/(8*c.pi) for r in product_rows(square,sin4)]
        A=IntervalTaylor(c,[r/math.factorial(j) for j,r in enumerate(Arows)])
        G=(A/n).exp()
        grows=[G[j]*math.factorial(j) for j in range(5)]
        # exp(A/N)-1 is kept as a signed integral. Subtraction from 1
        # would erase its source when the finite N exceeds precision.
        x=A[0]/n;mag=upper(c,abs(x))
        dgrows=[x*c.exp(c.mpf([-endpoints(mag)[1],endpoints(mag)[1]]))]+grows[1:]
        Enew=product_rows(old['theta'],grows)
        einc=product_rows(old['theta'],dgrows)
        Vnew=[r*(-c.sqrt(mu)/(2*c.pi*n)) for r in product_rows(old['theta'],product_rows(chi,cos2))]
        sig=sigma_jets(c,v)[0] if region=='O3_slope_mu' else c.mpf(0)
        da=mu*chi[0]**2*c.cos(4*c.pi*phase)
        slowA=-mu*chi[0]*chi[1]*c.sin(4*c.pi*phase)/(4*c.pi)
        logEy=source['log_Utheta_ordinary_y_derivatives'][0][0]
        beta=2*c.sqrt(mu)*chi[0]*c.sin(2*c.pi*phase)
        slowB_over_E=-c.sqrt(mu)*(chi[0]*logEy+chi[1])*c.cos(2*c.pi*phase)/(2*c.pi)
        aexcess=2*mu*sig+da-2*slowA/n
        b=c.exp(-A[0]/n)*(beta+2*slowB_over_E/n)
        u=Vnew[0];e=einc[0];E=old['theta'][0]
        densities={
          'M':dict(Pstar_power=1,logR_power=0,sqrt2_power=0,coefficient=u),
          'I':dict(Pstar_power=1,logR_power='.5',sqrt2_power=1,coefficient=e),
          'J':dict(Pstar_power=2,logR_power='.5',sqrt2_power=1,coefficient=u*(E+e)),
          'S':dict(Pstar_power=2,logR_power=0,sqrt2_power=0,coefficient=u*u-E*e-e*e/2),
          'Cp':dict(Pstar_power=2,logR_power=-1,sqrt2_power=0,coefficient=E*e+e*e/2)}
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
          datum_enclosure_sha256=self.datum_sha,region=region,coordinate=v,shared_logR_offset=t,
          finite_integer_N=N,actual_phase=n*t,phase_is_actual_translated_N_logR=True,
          spatial_cutoff_ordinary_logR_rows=chi,
          original_E0_over_Pstar_rows=old['theta'],
          modified_theta_over_Pstar_ordinary_logR_axial5=Enew,
          modified_axial_over_Pstar_ordinary_logR_axial5=Vnew,
          swirl_increment_over_Pstar_ordinary_logR_axial5=einc,
          exact_log_source_Pstar=c.mpf(self.pre.params.logPstar),
          finite_N_shear=dict(a_minus2=aexcess,b=b,slow_DXA=slowA,slow_DXB_over_E0=slowB_over_E),
          actual_five_moment_increment_densities_per_dX=densities,
          original_axial_is_exact_zero=True,
          original_profile_jets_unchanged_outside_support=endpoints(t)[1]<=-2 or endpoints(t)[0]>=.5,
          modified_profiles_available_as_local_candidate=True,
          moments_pressure_radial_recovery_and_completed_tensor_not_installed=True,
          **dict.fromkeys(OPEN,False))

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
          datum_enclosure_sha256=self.datum_sha,exact_finite_frequency_profile_theorem=self.theorem,
          supported_periodic_loop_whole_O3_vs_excess_lower=self.supported_loop_vs_excess_lower,
          supported_loop_floor_proof='v<=1/4: chi=1 and vs-2>=mu/2; v>=1/4: sigma(v)>=sigma(1/4), vs-2>=2mu*sigma(1/4).',
          supported_loop_uses_same_current_inviscid_direction_bound=self.direction.bounds['source_correlated_periodic_shear_loop'],
          source_bound_seam_crossing_local_finite_N_profiles_available=True,
          phase_samples_not_used_to_certify_any_cone=True,
          original_five_moments_not_reused_for_modified_profiles=True,
          current_strict_nonzero_whole_regions_including_inherited=15,
          remaining_registry_regions_without_current_whole_cone=17,
          scope='Local finite-N E/U profiles with actual phase and ordinary logR/Z jets on both sides of O2/O3 seam. Five exact increment density sources factored by current Pstar/R. Whole O3 periodic target loop stays strict but actual finite-N pN, integrated moments, shared pressure, independent five-bump repair, modified radial velocity, O2 taper and a common admitted finite N remain open. Original counts do not increase.',
          input_hashes=self.hashes,**dict.fromkeys(OPEN,False))

@source_precision
def run(field):
    result=field.manifest()
    # A finite working frequency exercises the construction, not admission.
    result['local_candidate_examples']=[field.profile(region,'.7',v,10**12) for region,v in (
       ('O2_buffer',9),('O2_buffer',10),('O2_buffer',11),
       ('O3_slope_mu',0),('O3_slope_mu','.1337000000001337'),('O3_slope_mu','.4'),('O3_slope_mu','.5'),('O3_slope_mu',1))]
    (HERE/NAME).write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual finite-N O2/O3 profile candidate constructed; uniform N, moments/pressure and cone remain open',flush=True)
    return result
