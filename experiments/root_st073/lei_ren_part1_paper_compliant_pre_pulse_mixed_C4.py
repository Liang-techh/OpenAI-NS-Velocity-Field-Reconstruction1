"""Actual repaired Rh -> O.1/O.2/O.3 -> Rp, physical logR/Z mixed4.

Reference terminal identities come from the ACTUAL inner five-bump patch.
All later primitive histories are transported, including after Uz vanishes.
Absolute Rref and exp(logCstar) are never materialized. Positive amplitudes
in these outer charts have finite logarithms and are evaluated directly;
the microscopic bridge's hb is not used as an outer coordinate width.
"""
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_outer_buffer import SharedOuterBuffer, transition_kernels, decay_integral
from lei_ren_part1_paper_compliant_outer_initial import turnoff_kernels
from lei_ren_part1_paper_interval_outer_slope_field import transition_integrals
from lei_ren_part1_paper_compliant_inner_bridge_profiles import square, derivative, logarithm
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import exponential_derivatives
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


def rate_rows(rows, rate, k):
    return sum((rows[j]*(math.comb(k,j)*rate**(k-j)) for j in range(k+1)),rows[0]*0)


def product_rows(a,b,k):
    return sum((a[j]*b[k-j]*math.comb(k,j) for j in range(k+1)),a[0]*0)


def physical_mixed(c,Z,delta,u,logU,V,history,p0,invP2):
    """Differentiate physical quantities, with fixed current-basepoint units.

    history = m=Mz/R, h=Mtheta/(sqrt2 R^1.5 Pstar), k=analogous
    mixed moment, e=Mztheta/(R Pstar^2), p=Mp/Pstar^2. The ordinary
    y derivative equations are exact; Taylor coefficients below concern Z.
    """
    z=IntervalTaylor.variable(c,c.mpf(Z),5);zero=z*0
    U=[u*row for row in exponential_derivatives(logU)]
    rows={name:[value] for name,value in history.items()}
    for j in range(4):
        rows['m'].append(V[j]-rows['m'][j])
        rows['h'].append(U[j]-rows['h'][j]*c.mpf('1.5'))
        rows['k'].append(product_rows(U,V,j)-rows['k'][j]*c.mpf('1.5'))
        rows['e'].append(product_rows(V,V,j)*invP2-product_rows(U,U,j)/2-rows['e'][j])
        rows['p'].append(product_rows(U,U,j)/2)
    Q=[(2*z*V[j]-(z*rows['m'][j])*(1-delta)-(1-square(z))*derivative(rows['m'][j]))/(1-square(z)*delta)
       for j in range(5)]
    physical=dict(Utheta_over_Pstar=U,Uz=V,
        Ur_over_current_sqrt_R_over_2=[rate_rows(Q,c.mpf('.5'),j) for j in range(5)],
        P_over_Pstar2=[p0+rows['p'][0]]+rows['p'][1:])
    primitives=dict(Mz_over_current_R=[rate_rows(rows['m'],1,j) for j in range(5)],
        Mtheta_over_current_sqrt2_R_1p5_Pstar=[rate_rows(rows['h'],c.mpf('1.5'),j) for j in range(5)],
        Mtheta_z_over_current_sqrt2_R_1p5_Pstar=[rate_rows(rows['k'],c.mpf('1.5'),j) for j in range(5)],
        Mztheta_over_current_R_Pstar2=[rate_rows(rows['e'],1,j) for j in range(5)],
        Mp_over_Pstar2=rows['p'])
    grid=lambda values:{'y'+str(j)+'_Z'+str(n):row[n]*math.factorial(n)
        for j,row in enumerate(values) for n in range(5-j)}
    return dict(physical_velocity_pressure_y_Z_mixed4={name:grid(values) for name,values in physical.items()},
        physical_five_primitive_y_Z_mixed4={name:grid(values) for name,values in primitives.items()},
        physical_velocity_pressure_y_derivative_Taylor={name:[row.truncate(4-j) for j,row in enumerate(values)] for name,values in physical.items()},
        actual_normalized_primitive_y_derivative_axial5=rows,actual_Q_y_derivative_axial4=Q,
        radial_prefactors_differentiated_before_mixed_grid=True)


def slope_masses(c,y,cells):
    """Positive monotone integral hull for every point in an interval y."""
    lo,hi=endpoints(c.mpf(y))
    if lo<0 or hi>1:raise ValueError('O2 slope y in[0,1] required')
    def point(value):
        sign,mantissa,exponent,_=value._mpf_
        exact=Fraction((-1 if sign else 1)*mantissa)*(Fraction(2)**exponent)
        return transition_integrals(c,exact,cells)
    left=point(lo);right=left if lo==hi else point(hi)
    hull=lambda a,b:c.mpf([endpoints(a)[0],endpoints(b)[1]])
    return hull(left[0],right[0]),[hull(a,b) for a,b in zip(left[1],right[1])]


def turnoff_derivatives(c,y,Md,phase):
    """D_y^k sigma(1-log(y)/Md), not D_phase^k, through four."""
    sig=sigma_jets(c,1-phase);ordinary=[sig[k]*math.factorial(k) for k in range(5)]
    # Signed Stirling numbers: (y D_y)^j and ordinary D_y^k conversion.
    signed=((1,),(0,1),(0,-1,1),(0,2,-3,1),(0,-6,11,-6,1))
    return [ordinary[0]]+[sum((ordinary[j]*((-1/c.mpf(Md))**j)*signed[k][j] for j in range(1,k+1)),c.mpf(0))/y**k for k in range(1,5)]


class CompliantPrePulseMixedC4:
    def __init__(self,cells=256,window=800):
        self.buffer=SharedOuterBuffer();self.initial=self.buffer.initial;self.ctx=c=self.initial.ctx
        self.family=self.initial.family;self.source=self.initial.datum.source_sha
        self.datum=self.initial.datum;self.delta=self.initial.delta;self.invP2=self.initial.invP2
        self.params=self.initial.params;self.cells=cells;self.window=window
        self.hashes=dict(self.buffer.hashes);self.cache={}
        for part in ('actual_moment_patch_check','actual_patch_mixed_C4_check','outer_buffer_check'):
            name=PREFIX+part+'.json';record=json.loads((HERE/name).read_bytes())
            if record.get('actual_five_defect_family_sha256')!=self.family:
                raise ValueError('Actual pre-pulse family mismatch: '+part)
            if record.get('implicit_source_sha256',self.source)!=self.source:
                raise ValueError('Actual pre-pulse source mismatch: '+part)
            if part=='actual_moment_patch_check' and not record.get('actual_five_functional_terminal_identities_connected'):
                raise ValueError('Actual incoming histories must have functional terminal closure')
            if part=='actual_patch_mixed_C4_check' and not (record.get('all_passed') and record.get('patch_Rm_and_Rh_functional_profile_joins_certified')):
                raise ValueError('Accepted actual Rh functional mixed4 join required')
            if part=='outer_buffer_check' and not (record.get('whole_axis_C1_buffer_API_checked') and record.get('retained_axial_history_at_pulse_inlet_checked')):
                raise ValueError('Accepted original whole-axis buffer source required')
            for path,digest in record['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Pre-pulse source changed: '+path)
            self.hashes.update(record['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        self.source_chain=dict(actual_terminal='actual_moment_patch: uniform unique implicit solution of the actual transported five defects',
            Rh='Rh=e*Rm=e^-5*Rref; beyond all original bump supports',
            reference='u/Pstar=exp(y/10)/(1+Z^2), V=4Z; y=log(R/Rref)',
            terminal_moments='m=4Z,h=5u/8,k=4Zh,e=16Z^2/Pstar^2-5u^2/12,p=5u^2/2',
            pressure='SAME analytic P0(Z) at all radii; no new pressure datum',
            O2='original Eq4.6 slope integral J(y), followed by sigma(1-log(y)/Md) turnoff and 11-unit buffer',
            O3='original slope-mu transition, then Tw=-60log(mu) pure-power buffer',
            histories='all five integral ODEs; no zero reset when Uz turns off')

    def coordinates(self,Z):
        c=self.ctx;Z=c.mpf(Z);lo,hi=endpoints(Z)
        if lo<-1 or hi>1:raise ValueError('Z in[-1,1] required')
        z=IntervalTaylor.variable(c,Z,5);qi=(1+square(z)).reciprocal()
        return z,qi

    def packet(self,Z,chart,coordinate,u,logu,logU,V,history,extra=None):
        c=self.ctx;raw=self.datum.normalized_jets(endpoints(c.mpf(Z)),5)['normalized_pressure_coefficients']
        p0=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in raw])
        z,qi=self.coordinates(Z);logshape=-logarithm(1+square(z))
        logjet=[logu]+list(logshape.coefficients[1:])
        data=physical_mixed(c,Z,self.delta,u,logU,V,history,p0,self.invP2)
        result=dict(Z=c.mpf(Z),chart=chart,coverage_coordinate=coordinate,**data,
            Utheta_over_Pstar_axial5_coefficients=list(u.coefficients),
            log_Utheta_over_Pstar_base_source=logu,
            log_Utheta_over_Pstar_axial5_coefficients=logjet,
            log_Utheta_ordinary_y_derivatives=logU,Uz_ordinary_y_derivative_axial5=V,
            original_P0_axial5_coefficients=list(p0.coefficients),
            exact_formal_prefactors=dict(Utheta='Pstar',Uz='1',Ur='sqrt(current R/2)',pressure='Pstar^2',
                Mz='current R',Mtheta='sqrt2*current R^1.5*Pstar',Mtheta_z='sqrt2*current R^1.5*Pstar',Mztheta='current R*Pstar^2',Mp='Pstar^2'),
            derivative_coordinate='y=logR and paper Z; phase/fraction only selects coverage',
            actual_five_histories_and_analytic_pressure_retained=True,actual_Rh_terminal_identities_bound=True,
            all_pre_pulse_mixed4_available=True,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False)
        if extra:result.update(extra)
        return result

    def reference(self,Z,offset):
        c=self.ctx;y=c.mpf(offset);lo,hi=endpoints(y)
        if lo<-5 or hi>0:raise ValueError('Rh-to-Rref offset in[-5,0] required')
        z,qi=self.coordinates(Z);u=qi*c.exp(y/10);zero=z*0;V=z*4
        h=u*c.mpf('.625');hist=dict(m=V,h=h,k=h*V,e=square(V)*self.invP2-square(u)*c.mpf(5)/12,p=square(u)*c.mpf('2.5'))
        return self.packet(Z,'Rh_to_Rref',y,u,y/10-c.ln(1+square(z)[0]),[zero+c.mpf('.1')]+[zero]*3,
            [V]+[zero]*4,hist,dict(exact_radius_source='R=Rref*exp(offset), Rh=e^-5 Rref'))

    def slope(self,Z,y):
        c=self.ctx;y=c.mpf(y);z,qi=self.coordinates(Z);zero=z*0
        J,mass=slope_masses(c,y,self.cells);factor=c.exp(y/10-c.mpf('.6')*J);u=qi*factor;V=z*4
        sig=sigma_jets(c,y);logU=[zero+c.mpf('.1')-sig[0]*c.mpf('.6')]+[zero-sig[k]*math.factorial(k)*c.mpf('.6') for k in range(1,4)]
        h=qi*(c.mpf('.625')+mass[0])*c.exp(-c.mpf('1.5')*y)
        hist=dict(m=V,h=h,k=h*V,e=square(V)*self.invP2-square(qi)*(c.mpf(5)/12+mass[2]/2)*c.exp(-y),
            p=square(qi)*(c.mpf('2.5')+mass[1]/2))
        return self.packet(Z,'O2_slope',y,u,y/10-c.mpf('.6')*J-c.ln(1+square(z)[0]),logU,[V]+[zero]*4,hist,
            dict(original_J=J,original_dimensionless_masses=mass,exact_radius_source='R=Rref*exp(y)'))

    def inlet(self,Z):
        key=tuple(c for c in self.ctx.mpf(Z)._mpi_)
        if key not in self.cache:self.cache[key]=self.slope(Z,1)
        return self.cache[key]

    def axial(self,Z,phase=None,buffer_offset=None):
        c=self.ctx;z,qi=self.coordinates(Z);zero=z*0;md=c.mpf(self.params.Md)
        if (phase is None)==(buffer_offset is None):raise ValueError('Select phase or buffer_offset')
        if phase is not None:
            phase=c.mpf(phase)
            if endpoints(phase)[0]<0 or endpoints(phase)[1]>1:raise ValueError('Axial turnoff phase in[0,1] required')
            y=c.exp(md*phase);B=turnoff_derivatives(c,y,md,phase);chart='O2_axial_turnoff';selector=phase
        else:
            selector=c.mpf(buffer_offset)
            if endpoints(selector)[0]<0 or endpoints(selector)[1]>11:raise ValueError('O2 buffer offset in[0,11] required')
            y=c.exp(md)+selector;B=[c.mpf(0)]*5;chart='O2_11_unit_buffer'
        parent=self.inlet(Z);u1=IntervalTaylor(c,parent['Utheta_over_Pstar_axial5_coefficients'])
        old={key:rows[0] for key,rows in parent['actual_normalized_primitive_y_derivative_axial5'].items()}
        t=y-1;d=c.exp(-t);root=c.exp(-t/2);d3=c.exp(-c.mpf('1.5')*t)
        u=u1*root;V=[z*(4*b) for b in B];K=turnoff_kernels(c,y,self.params.Md,self.cells,self.window)
        hist=dict(m=old['m']*d+z*(4*K['B_mass']),h=old['h']*d3+u1*(root-d3),
            k=old['k']*d3+u1*z*(4*root*K['B_mass']),
            e=old['e']*d+square(z)*(16*self.invP2*K['B_squared_mass'])-square(u1)*(t*d/2),
            p=old['p']+square(u1)*((1-d)/2))
        logu=parent['log_Utheta_over_Pstar_base_source']-t/2
        return self.packet(Z,chart,selector,u,logu,[zero-c.mpf('.5')]+[zero]*3,V,hist,
            dict(actual_y=y,original_turnoff_kernel_enclosures=K,original_cutoff_ordinary_y_derivatives=B,
                exact_radius_source='R=Rref*exp(y); y=exp(Md*phase) or exp(Md)+offset',
                source_B='sigma(1-log(y)/Md)',retained_far_tail_not_reset=True))

    def slope_mu(self,Z,offset):
        c=self.ctx;t=c.mpf(offset);z,qi=self.coordinates(Z);zero=z*0;mu=c.mpf(self.params.mu)
        if endpoints(t)[0]<0 or endpoints(t)[1]>1:raise ValueError('O3 transition offset in[0,1] required')
        parent=self.axial(Z,buffer_offset=11);u1=IntervalTaylor(c,parent['Utheta_over_Pstar_axial5_coefficients'])
        old={key:rows[0] for key,rows in parent['actual_normalized_primitive_y_derivative_axial5'].items()}
        K=transition_kernels(c,t,mu,self.cells);sig=sigma_jets(c,t)
        factor=c.exp(-t/2-mu*K['J']);d=c.exp(-t);d3=c.exp(-c.mpf('1.5')*t);u=u1*factor
        logU=[zero-c.mpf('.5')-mu*sig[0]]+[zero-mu*sig[k]*math.factorial(k) for k in range(1,4)]
        hist=dict(m=old['m']*d,h=old['h']*d3+u1*(d3*K['theta']),k=old['k']*d3,
            e=old['e']*d-square(u1)*(d*K['energy']/2),p=old['p']+square(u1)*(K['pressure']/2))
        return self.packet(Z,'O3_slope_mu',t,u,parent['log_Utheta_over_Pstar_base_source']-t/2-mu*K['J'],logU,[zero]*5,hist,
            dict(original_transition_kernels=K,selected_positive_mu=mu,exact_radius_source='R=Rd*exp(offset)',
                exact_logU_source='logu_d-t/2-mu*J(t); mu remains positive even below arithmetic addition precision'))

    def power(self,Z,phase):
        c=self.ctx;phase=c.mpf(phase);z,qi=self.coordinates(Z);zero=z*0;mu=c.mpf(self.params.mu)
        if endpoints(phase)[0]<0 or endpoints(phase)[1]>1:raise ValueError('O3 power phase in[0,1] required')
        parent=self.slope_mu(Z,1);u1=IntervalTaylor(c,parent['Utheta_over_Pstar_axial5_coefficients'])
        old={key:rows[0] for key,rows in parent['actual_normalized_primitive_y_derivative_axial5'].items()}
        t=self.params.Tw*phase;f=c.exp((-c.mpf('.5')-mu)*t);d=c.exp(-t);d3=c.exp(-c.mpf('1.5')*t)
        theta=(f-d3)/(1-mu);theta=c.mpf([max(mp.mpf(0),endpoints(theta)[0]),endpoints(theta)[1]])
        u=u1*f;hist=dict(m=old['m']*d,h=old['h']*d3+u1*theta,k=old['k']*d3,
            e=old['e']*d-square(u1)*(d*decay_integral(c,2*mu,t)/2),
            p=old['p']+square(u1)*(decay_integral(c,1+2*mu,t)/2))
        return self.packet(Z,'O3_power_to_Rp',phase,u,parent['log_Utheta_over_Pstar_base_source']+(-c.mpf('.5')-mu)*t,
            [zero-c.mpf('.5')-mu]+[zero]*3,[zero]*5,hist,
            dict(local_offset=t,total_log_length=self.params.Tw,exact_radius_source='R=Rw*exp(Tw*phase); Rp=Rw*exp(Tw)',
                pulse_inlet=endpoints(phase)==(mp.mpf(1),mp.mpf(1)),selected_positive_mu=mu,
                exact_logU_source='logu_w+(-1/2-mu)*Tw*phase; derivative coordinate is logR'))

    def evaluate(self,chart,Z,value):
        routes={'reference':self.reference,'slope':self.slope,'axial':self.axial,
            'buffer':lambda z,t:self.axial(z,buffer_offset=t),'slope_mu':self.slope_mu,'power':self.power}
        if chart not in routes:raise ValueError('Unknown pre-pulse source chart: '+str(chart))
        return routes[chart](Z,value)

    def report(self):
        charts={'reference':[-5,0],'slope':[0,1],'axial':[0,1],'buffer':[0,11],'slope_mu':[0,1],'power':[0,1]}
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum.datum_sha,actual_source_chain=self.source_chain,
            whole_charts={name:self.evaluate(name,[-1,1],domain) for name,domain in charts.items()},
            join_packets={name:self.evaluate(chart,[-1,1],value) for name,chart,value in (
                ('Rh','reference',-5),('Rref_left','reference',0),('Rref_right','slope',0),
                ('slope_exit','slope',1),('turnoff_inlet','axial',0),('turnoff_exit','axial',1),
                ('buffer_inlet','buffer',0),('Rd_left','buffer',11),('Rd_right','slope_mu',0),
                ('Rw_left','slope_mu',1),('Rw_right','power',0),('Rp','power',1))},
            all_pre_pulse_mixed4_available=True,actual_Rh_terminal_identities_bound=True,
            pre_pulse_functional_mixed4_joins_certified=False,Rp_accepted_pulse_source_join_certified=False,
            point_core_coefficients_reconstructed=False,full_inner_dispatcher_installed=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,
            temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(280):result=CompliantPrePulseMixedC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Actual Rh-to-Rp full six-chart physical mixed4 generated; all five histories/P0 retained',flush=True)
    return result


if __name__=='__main__':run()
