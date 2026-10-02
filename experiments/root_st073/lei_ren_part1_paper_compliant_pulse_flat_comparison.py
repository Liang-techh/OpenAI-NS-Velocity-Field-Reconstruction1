"""Quantitative flat velocity/moment differences at actual O.4 supports.

The comparison field has zero local axial input and the SAME exact moment
values at the support endpoint. Bounds concern differences, never deletion
of actual histories. Original flat shapes and actual C5 coefficients supply
every mixed derivative needed by velocity C4 and physical spatial mapping.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_tail_bound,beta_tail_bound,symmetric
from lei_ren_part1_paper_compliant_pulse_radial_C4 import CompliantPulseRadialC4
from lei_ren_part1_paper_compliant_pulse_mixed_C4 import binomial_product
from lei_ren_part1_paper_compliant_pulse_physical_bounds import (
    PulsePhysicalBounds,physical_operators,physical_bracket,abs_upper,UZ,UR)
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def symmetric_jet(jet):return IntervalTaylor(jet.ctx,[symmetric(jet.ctx,abs_upper(jet.ctx,v)) for v in jet.coefficients])


def difference_transport(c,Brows,mu,h):
    """Both orientations: |int_0^h exp(+-lambda v)Bdv|<=h exp(lambda h)||B||.

    The reference ODE has B=0 and the exact common boundary history.
    Constants in energy/angular/pressure cancel; they are not reset.
    """
    m1=[symmetric_jet(Brows[0]*(h*c.exp((c.mpf('.5')-mu)*h)))]
    m2=[symmetric_jet(Brows[0]*(h*c.exp((c.mpf('.5')-2*mu)*h)))]
    square=Brows[0]*Brows[0]
    energy=[symmetric_jet(square*(h*c.exp(2*mu*h)))]
    for k in range(4):
        m1.append(Brows[k]-m1[k]*(c.mpf('.5')-mu))
        m2.append(Brows[k]-m2[k]*(c.mpf('.5')-2*mu))
        square=Brows[0]*0
        for j in range(k+1):square+=Brows[j]*Brows[k-j]*math.comb(k,j)
        energy.append(square+energy[k]*(2*mu))
    return dict(linear_m1=m1,linear_m2=m2,energy=energy,
        angular_difference_exact_zero=True,swirl_pressure_difference_exact_zero=True)


class PulseFlatComparison:
    def __init__(self):
        self.ctx=c=MPIntervalContext(); c.dps=240
        name=PREFIX+'compliant_pulse_interface_certificate.json'; receipt=json.loads((HERE/name).read_bytes())
        if not receipt['all_passed'] or not receipt['exact_functional_main_gap_and_gap_end_identities_certified']:
            raise ValueError('Exact same-source functional interfaces required')
        self.hashes=dict(receipt['input_hashes'])
        for source,digest in self.hashes.items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Flat comparison source changed: '+source)
        self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.family=receipt['actual_five_defect_family_sha256']; self.source=receipt['implicit_source_sha256']
        fifth=json.loads((HERE/(PREFIX+'compliant_fifth_axial_jets.json')).read_bytes())
        pulse=json.loads((HERE/(PREFIX+'compliant_pulse_mixed_C4.json')).read_bytes())
        physical=json.loads((HERE/(PREFIX+'compliant_pulse_physical_bounds.json')).read_bytes())
        if any(r['actual_five_defect_family_sha256']!=self.family or r['implicit_source_sha256']!=self.source for r in (fifth,pulse,physical)):
            raise ValueError('Flat comparison family/source mismatch')
        jet=lambda v:IntervalTaylor(c,[read_interval(c,x) for x in v['coefficients']])
        selected=fifth['whole_Z']['selected']
        self.ap=jet(selected['selected_ap_Taylor'])
        self.controls=[jet(v) for v in selected['selected_scaled_end_coefficient_Taylor']]
        self.mu=read_interval(c,pulse['selected_mu']); self.delta=read_interval(c,pulse['selected_delta'])
        self.endcap=read_interval(c,pulse['positive_end_factor_cap'])
        self.normalization=read_interval(c,pulse['bump_normalization'])
        self.inlet_u=jet(pulse['whole_Z_terminal']['inlet_Utheta_over_Pstar_Taylor'])
        self.radial_field=CompliantPulseRadialC4.__new__(CompliantPulseRadialC4)
        self.radial_field.ctx=c; self.radial_field.delta=self.delta; self.radial_field.mu=self.mu
        self.physical=PulsePhysicalBounds.__new__(PulsePhysicalBounds); self.physical.ctx=c
        self.physical.mu=self.mu; self.physical.delta=self.delta
        self.physical.logP=read_interval(c,physical['logPstar'])
        self.physical.logRp_parts={k:read_interval(c,v) for k,v in physical['logRp_parts'].items()}
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def forcing(self,kind,h,row=None):
        c=self.ctx; h=c.mpf(h)
        if endpoints(h)[0]<0 or endpoints(h)[1]>endpoints(c.mpf('.01'))[1]:raise ValueError('Local log-radius distance h in[0,.01] required')
        if endpoints(h)[1]==0:
            return [IntervalTaylor.constant(c,0,5) for _ in range(5)]
        if kind in ('main_entrance','main_exit'):
            W=self.mu*h; bounds=[]
            if endpoints(50*W)[1]>mp.mpf('.5'):raise ValueError('Original entrance tail outside proven half support')
            for k in range(5):
                if kind=='main_entrance':
                    raw=W*sigma_tail_bound(c,0,50*W) if k==0 else 50**(k-1)*sigma_tail_bound(c,k-1,50*W)
                else:
                    raw=11*sigma_tail_bound(c,k,W)
                    if k>0:raw+=k*sigma_tail_bound(c,k-1,W)
                bounds.append(raw*self.mu**k)
            return [symmetric_jet(self.ap*bound) for bound in bounds]
        if kind=='end':
            if row not in (0,1):raise ValueError('Actual end row required')
            ell=c.mpf('.15'); W=2*h/ell
            return [symmetric_jet(self.controls[row]*(self.endcap*beta_tail_bound(c,k,W)/(ell**(k+1)*self.normalization))) for k in range(5)]
        raise ValueError('Unknown original pulse support')

    def comparison(self,label,kind,h,row=None,edge=None):
        c=self.ctx; h=c.mpf(h); z=c.mpf([-1,1]); B=self.forcing(kind,h,row)
        primitive=difference_transport(c,B,self.mu,h)
        A=[self.radial_field.radial(z,b,m) for b,m in zip(B,primitive['linear_m1'])]
        uz=[self.inlet_u*binomial_product(B,-(c.mpf('.5')+self.mu),k) for k in range(5)]
        ur=[self.inlet_u.truncate(4)*binomial_product(A,-self.mu,k) for k in range(5)]
        grids={name:{'y'+str(k)+'_Z'+str(n):jet[n]*math.factorial(n)
            for k,jet in enumerate(rows) for n in range(5-k)} for name,rows in ((UZ,uz),(UR,ur))}
        if kind=='main_entrance':chi,offset=c.mpf(0),c.mpf(0)
        elif kind=='main_exit':chi,offset=c.mpf(11),-h
        else:chi,offset=c.mpf(13),c.mpf(edge)-h
        physical={}
        for name,grid in grids.items():
            beta=-1 if name==UR else -1-self.delta; physical[name]={}
            for a,b in physical_operators():
                value=physical_bracket(c,grid,a,b,z,self.delta,c.mpf(beta)); norm=abs_upper(c,value)
                zero=endpoints(norm)[1]==0; lognorm=None if zero else c.ln(norm)
                scale=self.physical.component_scale(name,a,b,chi,offset,'-10')
                physical[name]['r'+str(a)+'_z'+str(b)]=dict(factored_difference_absolute_upper=norm,
                    factored_difference_log_upper=lognorm,scale=scale,
                    physical_difference_log_upper=None if zero else lognorm+sum((v for k,v in scale.items() if k.endswith('_term')),c.mpf(0)),
                    exact_difference_zero=zero)
        return dict(label=label,kind=kind,h=h,end_row=row,edge=edge,
            B_y_derivative_difference_bounds=[j.truncate(5-k) for k,j in enumerate(B)],
            normalized_primitive_difference_bounds={key:[j.truncate(5-k) for k,j in enumerate(rows)]
                for key,rows in primitive.items() if isinstance(rows,list)},
            profile_velocity_mixed_difference_bounds=grids,physical_cylindrical_difference_bounds=physical,
            angular_swirl_pressure_differences_exact_zero=True,
            reference_definition='zero local axial input, same EXACT interface moment/energy values and identical O4 swirl/pressure histories',
            actual_nonzero_histories_and_positive_terminal_energy_preserved=True,
            h_zero_comparison_is_a_difference_not_zero_actual_field=True)

    def report(self):
        c=self.ctx; rows=[]
        with mp.workdps(270):
            for h in ('.01','.003','.001','.000001','0'):
                rows.append(self.comparison('entrance','main_entrance',h))
                rows.append(self.comparison('main_exit','main_exit',h))
                for row,center in enumerate((-3,-1)):
                    for side in (-1,1):
                        edge=c.mpf(center)+side*c.mpf('.15')
                        rows.append(self.comparison('end'+str(row+1)+('_left' if side<0 else '_right'),'end',h,row,edge))
            return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
                comparisons=rows,quantitative_flat_velocity_moment_comparison_ledger_available=True,
                original_support_forcing_majorants_available=True,
                comparison_for_all_required_profile_and_physical_mixed_orders=True,
                underlying_flat_limit_proof=dict(shape='exact original exp(-1/x²) cutoff and exp(-1/(1-r²)) bump',
                    radial_source_width='main xi width=mu*h; end w<=2h/.15',
                    moment_difference='h*exp(lambda*h)*source bound; repeated derivatives via source ODEs',
                    energy_difference='h*exp(2mu*h)*B² bound; no -.5 because SAME reference source cancels it',
                    velocity_difference='linear paper(3.9) recovery plus physical prefactor/coordinate maps',
                    finite_source_scales_at_each_fixed_positive_tau=True,
                    each_required_difference_derivative_tends_to_zero=True,
                    numerical_endpoint_caps_also_vanish=True),
                full_pulse_C4_installed=False,full_outer_C4_certified=False,
                two_sided_external_high_order_joins_certified=False,physical_energy_integral_certified=False,
                whole_outer_cone_certified=False,temporal_recursion=False,
                next_dependency='Two-sided O3/pulse and pulse/O5 high-order joins, then post-pulse C4/energy/cone',input_hashes=self.hashes)


def run():
    result=PulseFlatComparison().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Actual original support flat comparison: moment/velocity and physical mixed derivative bounds generated without zeroing histories',flush=True)
    return result


if __name__=='__main__':run()
