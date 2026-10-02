"""Original microscopic switches: physical phase/Z4 and formal logR scales.

hb and hb^-1 are never materialized. All D_s bounds enclose the original
s=log(R/100)/hb source, with exact D_y^k=hb^-k D_s^k conversion ledgers.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_inner_switch_profiles import (
    CompliantInnerSwitchProfiles,IntervalTaylor,MTH,MTHZ,MZ,MZT,MP,power_transport,as_initial)
from lei_ren_part1_paper_compliant_inner_bridge_profiles import square,derivative,dress
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import exponential_derivatives,scaled_positive_source
from lei_ren_part1_paper_compliant_long_reshape_mixed_C4 import reshape_mixed
from lei_ren_part1_paper_compliant_inner_bridge_profiles import logarithm
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'


class FactoredJet:
    """Finite source polynomial; hb/Pstar/F0/u-base stay in log factors.

    Exponents multiply (loghb, 2logPstar, 2logF0base, 2logu-base).
    All bases are fixed at the Z basepoint. True Z derivatives belong to
    the coefficient jets. No cap is taken while building Bell/Leibniz rows.
    """
    def __init__(self,algebra,terms,order=5):
        self.algebra=algebra;self.ctx=algebra.ctx;self.order=order
        self.terms={key:row.truncate(order) for key,row in terms.items()
            if any(endpoints(value)!=(0,0) for value in row.coefficients)}
    def coerce(self,value):return self.algebra.lift(value,self.order)
    def __add__(self,value):
        other=self.coerce(value);order=min(self.order,other.order);terms=dict(self.terms)
        for key,row in other.terms.items():terms[key]=terms[key]+row if key in terms else row
        return FactoredJet(self.algebra,terms,order)
    __radd__=__add__
    def __neg__(self):return FactoredJet(self.algebra,{key:-row for key,row in self.terms.items()},self.order)
    def __sub__(self,value):return self+-self.coerce(value)
    def __rsub__(self,value):return self.coerce(value)+-self
    def __mul__(self,value):
        other=self.coerce(value);terms={};order=min(self.order,other.order)
        for a,left in self.terms.items():
            for b,right in other.terms.items():
                key=tuple(x+y for x,y in zip(a,b));row=left*right
                terms[key]=terms[key]+row if key in terms else row
        return FactoredJet(self.algebra,terms,order)
    __rmul__=__mul__
    def reciprocal(self):
        if set(self.terms)!={(0,0,0,0)}:raise ValueError('Only unscaled coefficient denominators may be inverted')
        return self.algebra.lift(self.terms[(0,0,0,0)].reciprocal())
    def __truediv__(self,value):return self*self.coerce(value).reciprocal()
    def __pow__(self,n):
        if n<0:raise ValueError('Negative factored powers must use explicit source shifts')
        out=self.coerce(1)
        for _ in range(n):out=out*self
        return out
    def Zderivative(self):
        return FactoredJet(self.algebra,{key:derivative(row) for key,row in self.terms.items()},self.order-1)


class FactoredAlgebra:
    def __init__(self,c,logs,proofs):self.ctx=c;self.logs=tuple(logs);self.proofs=proofs;self.final_rows=[]
    def lift(self,value,order=5):
        if isinstance(value,FactoredJet):
            if value.algebra is not self:raise ValueError('Source log bases differ')
            return value
        row=value if isinstance(value,IntervalTaylor) else IntervalTaylor.constant(self.ctx,value,order)
        return FactoredJet(self,{(0,0,0,0):row},row.order)
    def shift(self,value,powers):
        value=self.lift(value)
        return FactoredJet(self,{tuple(a+b for a,b in zip(key,powers)):row for key,row in value.terms.items()},value.order)
    def width(self,value,power=1):return self.shift(value,(power,0,0,0))
    def resolve(self,value,label=None):
        value=self.lift(value);out=IntervalTaylor.constant(self.ctx,0,value.order);ledger=[]
        for key,row in sorted(value.terms.items()):
            logbase=sum((self.logs[k]*power for k,power in enumerate(key) if power),self.ctx.mpf(0))
            out=out+scaled_positive_source(self.ctx,logbase,row,self.proofs)
            ledger.append(dict(source_exponents=list(key),combined_positive_source_log=logbase,
                final_coefficient_jet=list(row.coefficients),all_Bell_Leibniz_and_axial_factors_already_combined=True))
        if label is not None:self.final_rows.append(dict(physical_row=label,terms=ledger))
        return out
    def resolve_coefficient(self,value,n,label):
        value=self.lift(value);out=self.ctx.mpf(0);ledger=[]
        for key,row in sorted(value.terms.items()):
            logbase=sum((self.logs[k]*power for k,power in enumerate(key) if power),self.ctx.mpf(0))
            coefficient=IntervalTaylor.constant(self.ctx,row[n],0)
            out+=scaled_positive_source(self.ctx,logbase,coefficient,self.proofs)[0]
            ledger.append(dict(source_exponents=list(key),combined_positive_source_log=logbase,
                final_ordinary_coefficient=row[n],all_Bell_Leibniz_and_axial_factors_already_combined=True))
        self.final_rows.append(dict(physical_row=label,terms=ledger))
        return out
    def tree(self,value):
        if isinstance(value,FactoredJet):return self.resolve(value)
        if isinstance(value,dict):return {key:self.tree(row) for key,row in value.items()}
        if isinstance(value,list):return [self.tree(row) for row in value]
        return value


class RadialTaylor:
    """Three y Taylor orders with directed axial Taylor coefficients."""
    def __init__(self,rows):self.rows=list(rows);self.ctx=rows[0].ctx;self.order=len(rows)-1
    def coerce(self,other):
        if isinstance(other,RadialTaylor):return other
        value=other if isinstance(other,IntervalTaylor) else self.rows[0]*0+other
        return RadialTaylor([value]+[value*0]*self.order)
    def __add__(self,other):
        other=self.coerce(other);return RadialTaylor([a+b for a,b in zip(self.rows,other.rows)])
    __radd__=__add__
    def __neg__(self):return RadialTaylor([-row for row in self.rows])
    def __sub__(self,other):return self+-self.coerce(other)
    def __rsub__(self,other):return self.coerce(other)+-self
    def __mul__(self,other):
        other=self.coerce(other)
        return RadialTaylor([sum((self.rows[j]*other.rows[k-j] for j in range(k+1)),self.rows[0]*0) for k in range(self.order+1)])
    __rmul__=__mul__
    def reciprocal(self):
        out=[self.rows[0].reciprocal()]
        for k in range(1,self.order+1):out.append(-out[0]*sum((self.rows[j]*out[k-j] for j in range(1,k+1)),out[0]*0))
        return RadialTaylor(out)
    def __truediv__(self,other):return self*self.coerce(other).reciprocal()
    def Zderivative(self):return RadialTaylor([derivative(row) for row in self.rows])
    def ordinary(self):return [row*math.factorial(k) for k,row in enumerate(self.rows)]


def comparison_radial_directions(c,Z,delta,phi,v,moments,p0,ratios,ratios2):
    """Enclose exact comparison moment ODEs after ORIGINAL smoothing ends.

    barF and barV are constant radially here; their own inherited moments
    are not reset to unsmoothed frozen moments or used as actual moments.
    """
    lift=lambda row:RadialTaylor([row]+[row*0]*3)
    z=lift(IntervalTaylor.variable(c,c.mpf(Z),6));d=1-z*z;L=1-z*z*delta
    update={MTH:(phi*2,2),MTHZ:(phi*v*2,2),MZ:(v,1),
        'A':(square(v),1),'B':(square(phi),2),MP:(square(phi),1)}
    initial={MTH:moments[MTH],MTHZ:moments[MTHZ],MZ:moments[MZ],
        'A':moments[MZT]['axial'],'B':moments[MZT]['swirl'],MP:moments[MP]}
    rows={}
    for name,value in initial.items():
        rhs,rate=update[name];ordinary=[value]+[(rhs-value*rate)*((-rate)**(k-1)) for k in range(1,4)]
        rows[name]=RadialTaylor([row/math.factorial(k) for k,row in enumerate(ordinary)])
    true=lambda radial,factors:RadialTaylor([dress(row,factors) for row in radial.rows])
    f=true(lift(phi),ratios);h=true(rows[MTH],ratios);k=true(rows[MTHZ],ratios)
    b=true(rows['B'],ratios2);p=true(rows[MP],ratios2);m=rows[MZ];a=rows['A'];vv=lift(v);pressure=lift(p0)
    W=1-(z*m)*(1-delta)-d*m.Zderivative()
    angular=h*(1-delta/2)-(z*h.Zderivative())*((1-delta)/2)-d*k.Zderivative()+(z*k)*(2*delta-1)
    result=dict(D_over_R=(-W+angular/(f*2))/L,
        hydro=(-W*vv+(m-z*m.Zderivative())*((1-delta)/2)+(z*a)*(2*delta)-d*a.Zderivative())/(L*2),
        pressure=((z*pressure)*(2*(1+delta))-d*pressure.Zderivative())/(L*2),
        swirl=(-(z*b)*(2*delta)+d*b.Zderivative()+(z*p)*(2*(1+delta))-d*p.Zderivative())/(L*2))
    return {name:radial.ordinary() for name,radial in result.items()}


def rate_rows(rows,rate,k):
    return sum((rows[j]*(math.comb(k,j)*rate**(k-j)) for j in range(k+1)),rows[0]*0)


def product_rows(a,b,k):
    return sum((a[j]*b[k-j]*math.comb(k,j) for j in range(k+1)),a[0]*0)


def switch_controls(c,branch,cutoff,Dbar_y,drive_y,quotient,scale_width):
    """Ordinary phase derivatives of a, logF and raw V through four."""
    if branch not in ('first','second'):raise ValueError('Original first/second branch required')
    sig=[cutoff[k]*math.factorial(k) for k in range(5)]
    complement=[1-sig[0]]+[-row for row in sig[1:]]
    tinyD=[scale_width(row,k+1) for k,row in enumerate(Dbar_y)]
    if branch=='first':a=tinyD
    else:a=[sum((tinyD[j]*(math.comb(k,j)*complement[k-j]) for j in range(k+1)),tinyD[0]*0)+sig[k]*c.mpf('.8') for k in range(4)]
    logF=[scale_width(row,1)*(-c.mpf('.5')) for row in a]
    ratio_rows=[quotient*row for row in exponential_derivatives(logF)]
    if branch=='second':Vphase=[quotient*0]*4
    else:
        Vphase=[]
        # hb^2 * D_s^j drive = hb^(j+2) * D_y^j drive;
        # physical Pstar^2/F0^2 scales have already been factored in logs.
        driver=[drive_y(j,j+2) for j in range(4)]
        for k in range(4):
            value=quotient*0
            for i in range(k+1):
                for j in range(k-i+1):
                    ell=k-i-j;weight=math.factorial(k)//(math.factorial(i)*math.factorial(j)*math.factorial(ell))
                    value+=ratio_rows[j]*driver[ell]*(complement[i]*weight)
            Vphase.append(-value)
    one=a[0]*0+1
    logU=[scale_width(one-a[0],1)/2]+[scale_width(row,1)*(-c.mpf('.5')) for row in a[1:]]
    return dict(a_phase_derivatives=a,logF_phase_derivatives=logF,
        logUtheta_phase_derivatives=logU,Uz_positive_phase_derivatives=Vphase)


def phase_physical(c,Z,delta,h,amp,logu0,logU,V,shapes,p0,invP2,scale_width,proofs):
    """Physical primitive equations in exact source phase; no hb inverse."""
    algebra=h.algebra if isinstance(h,FactoredJet) else None
    lift=algebra.lift if algebra else lambda value:value
    z=lift(IntervalTaylor.variable(c,c.mpf(Z),5));zero=z*0
    amp=lift(amp);V=[lift(row) for row in V];shapes={name:lift(row) for name,row in shapes.items()};p0=lift(p0)
    U=exponential_derivatives(logU)
    shifted=lambda rate,m:exponential_derivatives([logU[0]*m+h*rate]+[row*m for row in logU[1:]])
    theta_rhs=shifted(c.mpf('1.5'),1);swirl_rhs=shifted(1,2);pressure_rhs=shifted(0,2)
    mean=[shapes['mean']]
    for k in range(4):mean.append(scale_width(V[k]-mean[k],1))
    dz=lambda row:row.Zderivative() if isinstance(row,FactoredJet) else derivative(row)
    Q=[(2*z*V[k]-(z*mean[k])*(1-delta)-(1-z*z)*dz(mean[k]))/(1-z*z*delta) for k in range(5)]
    amp2=amp*amp
    scaled=(lambda row:algebra.shift(amp2*row,(0,0,0,1))) if algebra else lambda row:scaled_positive_source(c,2*logu0,amp2*row,proofs)
    pressure=[scaled(shapes['pressure'])/2]+[scaled(scale_width(pressure_rhs[k-1],1))/2 for k in range(1,5)]
    primitives=dict(Mtheta_over_current_sqrt2_R_1p5_Utheta=[amp*shapes['theta']]+[amp*scale_width(theta_rhs[k-1],1) for k in range(1,5)],
        Mtheta_z_over_current_sqrt2_R_1p5_Utheta=[amp*shapes['theta_z']]+[amp*scale_width(product_rows(theta_rhs,V,k-1),1) for k in range(1,5)],
        Mz_over_current_R=[shapes['mean']]+[scale_width(rate_rows(V,h,k-1),1) for k in range(1,5)],
        Mztheta_over_current_R_Pstar2=[shapes['axial']*invP2-scaled(shapes['swirl'])/2]
            +[scale_width(rate_rows([product_rows(V,V,n) for n in range(4)],h,k-1),1)*invP2
              -scaled(scale_width(swirl_rhs[k-1],1))/2 for k in range(1,5)],Mp_over_Pstar2=pressure)
    physical=dict(Utheta_over_current_Utheta=[amp*row for row in U],Uz=V,
        Ur_over_current_sqrt_R_over_2=[rate_rows(Q,h/2,k) for k in range(5)],
        P_over_Pstar2=[p0+pressure[0]]+pressure[1:])
    def grid(name,rows):
        result={}
        for k,row in enumerate(rows):
            for n in range(5-k):
                key='s'+str(k)+'_Z'+str(n)
                result[key]=(algebra.resolve_coefficient(row*math.factorial(n),n,name+'/'+key)
                    if algebra else row[n]*math.factorial(n))
        return result
    return dict(physical_velocity_pressure_phase_Z_mixed4={name:grid(name,rows) for name,rows in physical.items()},
        physical_five_primitive_phase_Z_mixed4={name:grid(name,rows) for name,rows in primitives.items()},
        actual_Q_phase_derivative_axial4=algebra.tree(Q) if algebra else Q,
        radial_prefactors_differentiated_before_mixed_grid=True,
        source_width_and_amplitude_products_capped_only_after_final_derivatives=bool(algebra),
        final_factored_physical_row_ledgers=algebra.final_rows if algebra else [])


class CompliantMicroswitchMixedC4:
    def __init__(self):
        self.switch=CompliantInnerSwitchProfiles();self.core=self.switch.core;self.ctx=c=self.switch.ctx
        self.family=self.switch.family;self.source=self.switch.source;self.hashes=dict(self.switch.hashes);self.proofs=[]
        for part in ('inner_switch_profiles_check','long_reshape_mixed_C4_check'):
            name=PREFIX+part+'.json';receipt=json.loads((HERE/name).read_bytes())
            if not receipt['all_passed'] or receipt['actual_five_defect_family_sha256']!=self.family or receipt['implicit_source_sha256']!=self.source:
                raise ValueError('Same actual switch/reshape prerequisites required')
            for path,digest in receipt['input_hashes'].items():
                if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Microswitch source changed: '+path)
            self.hashes.update(receipt['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.logh=self.switch.logh;self.invP2=c.exp(-2*self.core.logP)
        self.logF0=c.mpf([endpoints(-self.core.logC-self.core.Lambda*self.core.Gbar)[0],endpoints(-self.core.logC)[1]])
        self.h=scaled_positive_source(c,self.logh,IntervalTaylor.constant(c,1,5),self.proofs)[0]
        self.distance=c.ln(100/self.switch.bridge.r)
        if endpoints(self.distance)[0]<=endpoints(2*self.switch.cap)[1]:raise ValueError('Comparison smoothing must end before R100')
        self.continuation_distance=c.ln(c.mpf('4.1')/4)
        if endpoints(self.continuation_distance)[0]<=endpoints(2*self.switch.cap)[1]:raise ValueError('Comparison smoothing leaves continuation box')
        self.shared_axial_source=json.loads((HERE/(PREFIX+'long_reshape_mixed_C4.json')).read_bytes())['shared_exact_axial_source']
        self.hashes[PREFIX+'long_reshape_mixed_C4.json']=hashlib.sha256((HERE/(PREFIX+'long_reshape_mixed_C4.json')).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def width(self,row,power=1,extra_log=0):
        return scaled_positive_source(self.ctx,self.logh*power+extra_log,row,self.proofs)

    def evaluate(self,Z,phase,branch):
        c=self.ctx;phase=c.mpf(phase)
        lo,hi=endpoints(phase);start,end=(0,1) if branch=='first' else (1,2)
        if lo<start or hi>end:raise ValueError('Original branch phase interval required')
        parent=self.switch.phase(Z,phase);inp=self.switch.bridge.inputs(Z)
        R=c.mpf(100)*c.exp(self.h*phase);comparison=self.switch.bridge.comparison(c.mpf(Z),self.switch.bridge.r/R)
        directions=comparison_radial_directions(c,Z,self.core.delta,comparison['phi'],comparison['v'],
            comparison['moments'],inp['p0'],inp['F0_ratios'],inp['F0_squared_ratios'])
        Dbar=[rate_rows(directions['D_over_R'],c.mpf(1),k)*R for k in range(4)]
        jet=lambda row:IntervalTaylor(c,row)
        phi=jet(parent['F_actual_over_F0_axial5_coefficients']);v=jet(parent['Uz_actual_axial5_coefficients'])
        quotient=phi/comparison['phi'].truncate(5)
        logu0=c.ln(2*R)/2+self.logF0+c.ln(phi[0])-self.core.logP
        algebra=FactoredAlgebra(c,(self.logh,2*self.core.logP,2*self.logF0,2*logu0),self.proofs)
        def drive_y(k,width_power):
            raw=(algebra.lift(rate_rows(directions['hydro'],c.mpf(1),k)*R)
                +algebra.shift(rate_rows(directions['pressure'],c.mpf(1),k)*R,(0,1,0,0))
                +algebra.shift(rate_rows(directions['swirl'],c.mpf(2),k)*R**2,(0,0,1,0)))
            return algebra.width(raw,width_power)
        controls=switch_controls(c,branch,sigma_jets(c,phase-start),Dbar,drive_y,algebra.lift(quotient),algebra.width)
        V=[v]+controls['Uz_positive_phase_derivatives']
        dressed=dress(phi,inp['F0_ratios']);amp=dressed/phi[0]
        amp=IntervalTaylor(c,[c.mpf(1)]+list(amp.coefficients[1:]))
        raw=parent['actual_moment_shape_axial5_coefficients']
        shapes=dict(theta=jet(raw[MTH])/(phi*2),theta_z=jet(raw[MTHZ])/(phi*2),mean=jet(raw[MZ]),
            axial=jet(raw[MZT]['axial']),swirl=jet(raw[MZT]['swirl'])/square(phi),pressure=jet(raw[MP])/square(phi))
        packet=phase_physical(c,Z,self.core.delta,algebra.width(1),amp,logu0,controls['logUtheta_phase_derivatives'],V,shapes,
            jet(parent['pressure_axis_axial5_coefficients']),algebra.shift(1,(0,-1,0,0)),algebra.width,self.proofs)
        ledger={}
        final_rows={row['physical_row']:row for row in packet['final_factored_physical_row_ledgers']}
        for group in ('physical_velocity_pressure_phase_Z_mixed4','physical_five_primitive_phase_Z_mixed4'):
            ledger[group]={}
            for name,grid in packet[group].items():
                bounds={}
                for key,value in grid.items():
                    k=int(key.split('_')[0][1:]);size=max(abs(q) for q in endpoints(value))
                    term_logs=[]
                    for term in final_rows[name+'/'+key]['terms']:
                        coefficient=term['final_ordinary_coefficient'];magnitude=max(abs(q) for q in endpoints(coefficient))
                        if not magnitude:continue
                        exponents=list(term['source_exponents']);exponents[0]-=k
                        combined=sum((algebra.logs[j]*power for j,power in enumerate(exponents) if power),c.mpf(0))
                        term_logs.append(combined+c.ln(c.mpf(magnitude)))
                    # Triangle bound in logs; hb powers cancel BEFORE any cap.
                    # max(log|term|)+log(number of nonzero terms) avoids even
                    # materializing huge positive logR derivative magnitudes.
                    logupper=(c.mpf(max(endpoints(row)[1] for row in term_logs))+c.ln(len(term_logs))) if term_logs else None
                    bounds[key.replace('s','y',1)]=dict(exact_zero=not term_logs,log_absolute_upper=logupper,
                        log_absolute_upper_from_capped_phase=None if size==0 else c.ln(c.mpf(size))-k*self.logh,
                        formal_width_powers_cancelled_before_log_bound=True,nonzero_factored_terms=len(term_logs),
                        exact_source_conversion='D_y^k=hb^-k*D_s^k',source_phase_order=k)
                ledger[group][name]=bounds
        return dict(Z=c.mpf(Z),phase=phase,branch=branch,**packet,controls=algebra.tree(controls),
            comparison_direction_ordinary_y_axial5=directions,actual_parent_axial5_packet=parent,
            physical_logR_Z_mixed4_log_bound_ledger=ledger,exact_width_source='hb=cstar*K^-100; s=log(R/100)/hb',
            exact_positive_width_log=self.logh,width_enclosure_is_not_source=self.h,
            exact_source_R='100*exp(hb*s)',R_enclosure_only=R,exact_positive_swirl_source_log=logu0,
            formal_radius_tree=dict(op='exp',log=dict(op='sum',args=['log(100)',dict(op='mul',args=['hb','s'])])),
            factored_source_log_bases=list(algebra.logs),
            actual_axial_function_source=(dict(op='sum',args=self.shared_axial_source['V100']['args']+[
                dict(formal_integral='I_first_switch',shared_source=self.shared_axial_source['shared_source_namespace'],
                    bounds=['0','s'],source_integrand=self.shared_axial_source['formal_signed_integrals']['I_first_switch'])])
                if branch=='first' else self.shared_axial_source['V110']),
            first_switch_derivative_source=dict(
                integrand='-hb^2*(1-sigma(s))*(phi_actual/barphi)*(R*hydro+R*Pstar^2*pressure+R^2*F0^2*swirl)',
                kth_positive_order='V_s^(k+1)=D_s^k integrand; full three-factor multinomial Leibniz rule; k=0..3',
                D_s_on_comparison_drive='hb*D_y; R=100exp(hb*s)',
                quotient_derivatives='(phi_actual/barphi)*Bell(logF_s,...); barphi constant radially here'),
            comparison_is_constant_Fbar_Vbar_after_original_smoothing=True,
            comparison_own_moments_retained=True,actual_moments_not_replaced_by_comparison=True,
            physical_derivative_normalization='fixed current-basepoint factors; phase and Z independent',
            inverse_hb_not_materialized=True,phase_derivatives_do_not_differentiate_width_caps=True,
            switch_phase_mixed4_available=True,switch_logR_mixed4_log_ledger_available=True,
            bridge_switch_inlet_join_certified=False,complete_postswitch_power_installed=False,
            full_inner_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False)

    def report(self):
        packets=dict(whole_first=self.evaluate([-1,1],[0,1],'first'),whole_second=self.evaluate([-1,1],[1,2],'second'),
            actual_R100_inlet=self.evaluate([-1,1],0,'first'),first_side_phase1=self.evaluate([-1,1],1,'first'),
            second_side_phase1=self.evaluate([-1,1],1,'second'),actual_R2_exit=self.evaluate([-1,1],2,'second'))
        samples=[self.evaluate(z,s,branch) for z,s,branch in (('0','.5','first'),('.5','1.5','second'))]
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,**packets,interior_packets=samples,
            whole_postswitch_power=self.postpower([-1,1],[0,1]),actual_R2_power_inlet=self.postpower([-1,1],0),
            actual_R110_power_exit=self.postpower([-1,1],1),
            comparison_smoothing_log_distance_to_R100=self.distance,source_positive_width_log=self.logh,
            comparison_smoothing_continuation_log_distance=self.continuation_distance,
            comparison_moment_enclosure_source=dict(
                smoothing_log_phi='alpha(y)*logPhi_core(y)+integral_0^y[-alpha_prime(t)]*logPhi_core(t)dt',
                smoothing_V='alpha(y)*V_core(y)+integral_0^y[-alpha_prime(t)]*V_core(t)dt',
                weights='alpha>=0; -alpha_prime>=0; alpha(y)+integral_0^y[-alpha_prime]=1; Z independent',
                radial_history='positive kernels integrate covered smooth histories, not endpoint constant values',
                normalized_weights=dict(H='1-theta^2',mean='1-theta',K='1-theta^2',A='1-theta',B='(1-theta^2)/2',C='1-theta'),
                axial_Taylor_orders=list(range(7)),exact_moment_point_values_reconstructed=False),
            shared_exact_axial_source=self.shared_axial_source,
            factored_width_and_amplitude_cap_proofs=self.proofs,formal_source_radius_not_rounded=True,
            comparison_own_moments_retained=True,actual_source_moment_history_retained=True,
            switch_phase_mixed4_available=True,switch_logR_mixed4_log_ledger_available=True,
            bridge_switch_inlet_join_certified=False,complete_postswitch_power_installed=True,
            full_inner_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)

    def postpower(self,Z,fraction):
        """Complete actual R2..110 with exact formal positive R2 retained."""
        c=self.ctx;fraction=c.mpf(fraction)
        if endpoints(fraction)[0]<0 or endpoints(fraction)[1]>1:raise ValueError('Postswitch original fraction in[0,1] required')
        inlet=self.switch.phase(Z,2);inp=self.switch.bridge.inputs(Z)
        jet=lambda row:IntervalTaylor(c,row)
        phi2=jet(inlet['F_actual_over_F0_axial5_coefficients']);V=jet(inlet['Uz_actual_axial5_coefficients'])
        raw=inlet['actual_moment_shape_axial5_coefficients']
        moments2={name:({key:jet(row) for key,row in value.items()} if name==MZT else jet(value)) for name,value in raw.items()}
        length=c.ln(c.mpf(110)/100)-2*self.h
        if endpoints(length)[0]<=0:raise ValueError('Original R2 must precede110')
        zeta=length*fraction;theta=c.exp(-zeta)
        phi=phi2*theta**c.mpf('.4')
        moments=power_transport(c,theta,as_initial(moments2),phi2,V)
        R=c.mpf(100)*c.exp(2*self.h+zeta)
        if endpoints(fraction)==(mp.mpf(1),mp.mpf(1)):R=c.mpf(110)
        else:R=c.mpf([max(mp.mpf(100),endpoints(R)[0]),min(mp.mpf(110),endpoints(R)[1])])
        dressed=dress(phi,inp['F0_ratios']);ell=logarithm(dressed)
        logu0=c.ln(2*R)/2+self.logF0+c.ln(phi[0])-self.core.logP
        logu=IntervalTaylor(c,[logu0]+list(ell.coefficients[1:]))
        constant=lambda value:IntervalTaylor.constant(c,value,5)
        shapes=dict(theta=moments[MTH]/(phi*2),theta_z=moments[MTHZ]/(phi*2),mean=moments[MZ],
            axial=moments[MZT]['axial'],swirl=moments[MZT]['swirl']/square(phi),pressure=moments[MP]/square(phi))
        packet=reshape_mixed(c,Z,self.core.delta,logu,[constant('.1')]+[constant(0)]*3,V,
            shapes,jet(inlet['pressure_axis_axial5_coefficients']),self.invP2,self.proofs)
        return dict(Z=c.mpf(Z),source_fraction=fraction,**packet,actual_R2_parent_axial5_packet=inlet,
            actual_postswitch_phi_axial5=list(phi.coefficients),actual_postswitch_moment_shapes={name:(
                {key:list(row.coefficients) for key,row in value.items()} if isinstance(value,dict) else list(value.coefficients)) for name,value in moments.items()},
            actual_postswitch_V_axial5=list(V.coefficients),original_P0_axial5=inlet['pressure_axis_axial5_coefficients'],
            source_positive_R2='100*exp(2*hb)',source_zeta='log(R/R2)=fraction*(log(110/100)-2hb)',
            source_theta='exp(-zeta)=R2/R; no rounded R2',zeta_enclosure_only=zeta,R_enclosure_only=R,
            theta_enclosure_only=theta,width_enclosure_is_not_source=self.h,
            formal_log_radius_tree=dict(op='sum',args=['log(100)','2hb',dict(op='mul',args=['fraction','log(110/100)-2hb'])]),
            actual_axial_function_source=self.shared_axial_source['V110'],
            source_angular_identity='log(F/F100)=-.4*log(R/100)+.6*hb-.5*hb^2*JD',
            exact_positive_width_log=self.logh,inverse_hb_not_materialized=True,
            source_caps_not_radius_or_velocity_values=True,actual_R2_histories_retained=True,
            complete_postswitch_power_installed=True,full_inner_interfaces_certified=False,
            full_cartesian_vector_derivatives_certified=False,admissible_stress_lift_constructed=False,temporal_recursion=False)


def run():
    with mp.workdps(280):result=CompliantMicroswitchMixedC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original microswitch phase/Z4 and exact formal logR derivative ledgers generated',flush=True)
    return result


if __name__=='__main__':run()
