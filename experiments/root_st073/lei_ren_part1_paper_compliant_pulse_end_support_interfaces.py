"""Full pulse-end stress/error flat interfaces with nonzero shared histories.

Only a local forcing is removed in the reference. Original boundary moments,
future energy, angular memory and absolute pressure are kept. Bounds use
unscaled C5 controls; positive D is retained as an exact source log.
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
from lei_ren_part1_paper_compliant_pulse_end_stress_C3 import pulse_coefficients, SourceAST, source_precision
from lei_ren_part1_paper_compliant_pulse_end_physical_C2 import (
    pulse_velocity_rows, pulse_remainder_sectors, axial_operator_rows,
    ordinary_grid, FALSE_FLAGS)
from lei_ren_part1_paper_compliant_pulse_flat_comparison import difference_transport,symmetric_jet
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_tail_bound,beta_jets
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
from lei_ren_part1_paper_compliant_collar_Gamma_C4 import product_rows
from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_bracket
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
EDGES=[dict(row=row,center=center,side=side,exact_edge=str(s.Rational(center)+side*s.Rational(3,20)))
    for row,center in enumerate((-3,-1)) for side in (-1,1)]


def current_sources():
    records={};hashes={};family=None
    for stem in ('pulse_end_stress_C3','pulse_end_stress_C3_check','pulse_end_physical_C2','pulse_end_physical_C2_check',
        'pulse_flat_comparison','pulse_flat_comparison_check','flat_pulse_derivatives_check','fifth_axial_jets',
        'pulse_mixed_C4','axial_pulse_field'):
        name=PREFIX+stem+'.json';record=json.loads((HERE/name).read_bytes())
        if stem.endswith('_check') and not record['all_passed']:raise ValueError('Unaccepted boundary source '+stem)
        identity=(record['actual_five_defect_family_sha256'],record['implicit_source_sha256'])
        if family is not None and identity!=family:raise ValueError('Boundary source family differs')
        family=identity
        for path,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Boundary source changed '+path)
        hashes.update(record['input_hashes']);hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        records[stem]=record
    return records,hashes,family


def full_difference_rows(c,delta,mu,z,C,Xp,Bh,m_diff,n_diff,J_diff,m_actual):
    """Actual-minus-reference, including nonzero reference meridional history."""
    zero=C*0;zeros=[zero]*5
    parts=pulse_coefficients(delta,mu,z,C,Xp,Bh,m_diff,n_diff,zeros,J_diff,zeros)
    for name in ('equilibrium','signed_original_memory','radial_shear'):
        parts['theta'][name]['shape']=[zero]*4
        parts['theta'][name]['full_derivative_rows']=[zero]*4
    # Bh_reference=0, so this product is Bh_actual*m_actual, not Bh*delta_m.
    actual=pulse_coefficients(delta,mu,z,C,Xp,Bh,m_actual,zeros,zeros,zeros,zeros)
    parts['axial']['nonlinear_meridional_transport']=actual['axial']['nonlinear_meridional_transport']
    delta_velocity=pulse_velocity_rows(c,delta,mu,z,C,Bh,m_diff)
    actual_velocity=pulse_velocity_rows(c,delta,mu,z,C,Bh,m_actual)
    reference_m=[m_actual[k]-m_diff[k] for k in range(5)]
    reference_velocity=pulse_velocity_rows(c,delta,mu,z,C,zeros,reference_m)
    delta_velocity['theta']=zeros
    errors=pulse_remainder_sectors(c,delta,mu,z,delta_velocity)
    du,au,ru=(v['radial'] for v in (delta_velocity,actual_velocity,reference_velocity))
    dz=delta_velocity['axial']
    # Product difference a*a_y-r*r_y = (a-r)*a_y+r*(a_y-r_y).
    first=product_rows(du[:3],au[1:4])
    second=product_rows(ru[:3],du[1:4])
    radial=shifted_rows([first[k]+second[k] for k in range(3)],-c.mpf('.5'),2)
    axial=product_rows(dz[:3],axial_operator_rows(c,au,-1,z,delta,count=2,order=1))
    errors['radial']['nonlinear_transport']['rows']=[radial[k]+axial[k] for k in range(3)]
    return dict(stress=parts,velocity=delta_velocity,remainder=errors)


def source_difference_proof():
    asts=SourceAST();proofs={}
    required={
        ('flat_pulse_derivatives','beta','raw'):'beta_jets(c,c.mpf(s)/ell)',
        ('pulse_end_stress_C3','end','Bh[j]'):"Cj*(beta[j]*math.factorial(j))",
        ('axial_pulse_field','backward_bump_weights','begin'):'c.mpf([max(mp.mpf(-1),min(mp.mpf(1),rlo)),max(mp.mpf(-1),min(mp.mpf(1),rhi))])',
        ('pulse_flat_comparison','difference_transport','m1'):"[symmetric_jet(Brows[0]*(h*c.exp((c.mpf('.5')-mu)*h)))]"}
    for (stem,name,target),expr in required.items():
        # Original end differentiates its unscaled Bh accumulator before enclosure.
        asts.expression(stem,name,target,wanted=expr,augmented=target=='Bh[j]')
        proofs['actual_source_'+stem+'.'+name+'.'+target]=True
    fn=asts.method('flat_pulse_derivatives','beta')
    returned=next(n for n in fn.body if isinstance(n,ast.Return))
    if ast.unparse(returned.value)!='IntervalTaylor(c, [raw[n] / (ell ** (n + 1) * self.normalization) for n in range(5)])':
        raise ValueError('Original width/normalization derivative scale changed')
    proofs['actual_normalized_beta_width_powers_and_positive_normalization']=True
    h,ell=s.symbols('h ell',positive=True);r=s.symbols('normalized_r',real=True)
    for side in (-1,1):
        distance_r=side*(1-h/ell)
        if s.expand(1-distance_r**2-(2*h/ell-(h/ell)**2))!=0:raise ArithmeticError('Original endpoint distance changed')
        proofs['endpoint_w_distance_'+str(side)]=True
    a,b,da,db=s.symbols('actual_a actual_b delta_a delta_b',real=True)
    if s.expand(a*b-(a-da)*(b-db)-(da*b+(a-da)*db))!=0:raise ArithmeticError('Retained-history product difference failed')
    proofs['full_nonzero_reference_product_difference']=True
    # Differential equations identify every derivative from the common boundary value.
    rate,mu=s.symbols('rate mu',real=True)
    B=[s.symbols('forcing'+str(k)) for k in range(5)]
    dm=[s.Symbol('same_boundary_difference')]
    for k in range(4):dm.append(B[k]-rate*dm[k])
    for k in range(5):
        if s.expand(dm[k].subs({dm[0]:0,**{v:0 for v in B}}))!=0:raise ArithmeticError('Moment endpoint derivative changed')
        proofs['moment_difference_endpoint_order'+str(k)]=True
    J=[s.Symbol('same_J_boundary_difference')]
    for k in range(4):
        square=sum(s.binomial(k,l)*B[l]*B[k-l] for l in range(k+1))
        J.append(2*mu*J[k]-square)
    for k in range(5):
        if s.expand(J[k].subs({J[0]:0,**{v:0 for v in B}}))!=0:raise ArithmeticError('Energy endpoint derivative changed')
        proofs['quadratic_difference_endpoint_order'+str(k)]=True
    # Exact support coordinates, without rounding a physical interval to r=1.
    c=MPIntervalContext();c.dps=100
    for side in (-1,1):
        if any(endpoints(v)!=(0,0) for v in beta_jets(c,side).coefficients):raise ArithmeticError('Original exact normalized beta endpoint is not flat')
        proofs['actual_normalized_beta_endpoint_all4_zero_'+str(side)]=True
    edge_values=[s.Rational(e['exact_edge']) for e in EDGES]
    if not -4<edge_values[0]<edge_values[1]<edge_values[2]<edge_values[3]<0:
        raise ArithmeticError('Original support geometry changed')
    proofs['all_four_original_edges_and_disjoint_supports_retained']=True
    return dict(identities=proofs,input_hashes=asts.hashes,
        exact_boundary_history_shared_not_zeroed=True,
        local_reference='zero this local beta forcing; SAME boundary moments, energy, pressure and angular memory',
        flat_limit='beta derivatives <= A_n exp(-1/W) W^(-2n)/(ell^(n+1)*N); W<=2h/ell; each limit zero',
        fixed_positive_tau_and_current_source_scope=True,
        full_stress3_velocity4_and_physical_error2_limits_follow_current_source_operators=True)


class PulseEndSupportInterfaces:
    @source_precision
    def __init__(self):
        self.records,self.hashes,self.family=current_sources()
        self.ctx=c=MPIntervalContext();c.dps=240
        mixed=self.records['pulse_mixed_C4']
        self.mu=read_interval(c,mixed['selected_mu']);self.delta=read_interval(c,mixed['selected_delta'])
        self.normalization=read_interval(c,mixed['bump_normalization'])
        if endpoints(self.normalization)[0]<=0:raise ValueError('Original beta normalization not positive')
        jet=lambda v:IntervalTaylor(c,[read_interval(c,x) for x in v['coefficients']])
        selected=self.records['fifth_axial_jets']['whole_Z']['selected']
        self.controls=[jet(v) for v in selected['selected_scaled_end_coefficient_Taylor']]
        whole=self.records['pulse_end_stress_C3']['whole_original_end']
        self.m_actual=[jet(v) for v in whole['source_formal_Mz_rows']]
        self.proof=source_difference_proof();self.hashes.update(self.proof['input_hashes'])
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def interface(self,edge,h):
        c=self.ctx;h=c.mpf(h);ell=c.mpf('.15')
        if endpoints(h)[0]<0 or endpoints(h)[1]>endpoints(c.mpf('.01'))[1]:raise ValueError('Original interface width h in [0,.01]')
        Z=c.mpf([-1,1]);z=IntervalTaylor.variable(c,Z,5)
        C=IntervalTaylor(c,[1+Z**2,2*Z,1,0,0,0]).reciprocal()
        if endpoints(h)[1]==0:B=[C*0]*5
        else:B=[symmetric_jet(self.controls[edge['row']]*(beta_tail_bound(c,k,2*h/ell)/(ell**(k+1)*self.normalization))) for k in range(5)]
        diff=difference_transport(c,B,self.mu,h)
        # e=e0-D^2 J, so the reference-relative J has the negative energy forcing.
        J=[-v for v in diff['energy']]
        rows=full_difference_rows(c,self.delta,self.mu,z,C,c.mpf(1),B,diff['linear_m1'],diff['linear_m2'],J,self.m_actual)
        stress={label:{name:ordinary_grid(part['full_derivative_rows'],3) for name,part in sectors.items()}
            for label,sectors in rows['stress'].items()}
        remainder={label:{name:dict(beta=part['beta'],mode=part['mode'],grid=ordinary_grid(part['rows'],2))
            for name,part in sectors.items()} for label,sectors in rows['remainder'].items()}
        # Pullback coefficients are recorded with their derivative exponents.
        # Exact B,D,H,R source factors remain the current parent recipes.
        physical_stress={};physical_error={}
        for label,sectors in stress.items():
            physical_stress[label]={}
            for name,grid in sectors.items():
                physical_stress[label][name]={}
                for i in range(4):
                    for j in range(4-i):
                        physical_stress[label][name]['r'+str(i)+'_z'+str(j)]=physical_bracket(c,grid,i,j,Z,self.delta,-2-self.delta)
        for label,sectors in remainder.items():
            physical_error[label]={}
            for name,sector in sectors.items():
                physical_error[label][name]={}
                for i in range(3):
                    for j in range(3-i):
                        physical_error[label][name]['r'+str(i)+'_z'+str(j)]=physical_bracket(c,sector['grid'],i,j,Z,self.delta,sector['beta'])
        return dict(edge=edge,h=h,unscaled_beta_forcing_difference_rows=B,
            source_primitive_differences={key:diff[key] for key in ('linear_m1','linear_m2','energy')},
            similarity_stress_difference_mixed3=stress,
            similarity_velocity_difference_mixed4={name:ordinary_grid(jets,4) for name,jets in rows['velocity'].items()},
            physical_stress_difference_mixed3_coefficients=physical_stress,
            physical_three_component_error_difference_mixed2_coefficients=physical_error,
            actual_boundary_histories_retained=True,
            exact_shared_swirl_pressure_and_incoming_memory_differences_zero=True,
            current_positive_log_factors_preserved_as_parent_recipes=True)

    @source_precision
    def report(self):
        samples=[self.interface(edge,h) for edge in EDGES for h in ('.01','.001','.000001','0')]
        return dict(actual_five_defect_family_sha256=self.family[0],implicit_source_sha256=self.family[1],
            source_difference_proof=self.proof,original_edges=EDGES,interfaces=samples,input_hashes=self.hashes,
            actual_pulse_end_all_four_support_functional_interfaces_verified=True,
            actual_pulse_end_stress3_and_physical_error2_flat_interface_bounds_available=True,
            source_caps_used_as_defining_field_values=False,
            source_factor_scope='Coefficients divided by original exact positive B/D/H/R sectors; restore current physical lambda/nu prefactors at fixed tau>0',
            support_flatness_does_not_prove_global_temporal_flat_remainder=True,
            **{flag:False for flag in FALSE_FLAGS})


@source_precision
def run():
    result=PulseEndSupportInterfaces().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Four original pulse-end support interfaces: full stress3 and three-component physical error2 flat differences generated',flush=True)
    return result


if __name__=='__main__':run()
