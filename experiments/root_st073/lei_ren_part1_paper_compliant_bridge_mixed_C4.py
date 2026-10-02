"""Original prescribed-shear bridge: both smoothing charts and frozen tail.

Mixed core derivatives missing from the old total5 packet are bounded by
the SAME admitted analytic Xh norm. No upstream source or field is changed.
Original widths/radii/integrals stay formal; caps only enclose final rows.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_microswitch_mixed_C4 import (
    CompliantMicroswitchMixedC4,FactoredAlgebra,FactoredJet,IntervalTaylor,
    exponential_derivatives,phase_physical,rate_rows,product_rows,MTH,MZ,MTHZ,MZT,MP)
from lei_ren_part1_paper_compliant_inner_bridge_profiles import square,derivative,dress
from lei_ren_part1_paper_compliant_core_physical_field import gridkey
from lei_ren_part1_paper_compliant_core_uniform_bounds import embedding
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import sigma_jets
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
S2={1:(0,1),2:(0,1,1),3:(0,1,3,1)}


class AxialSixAlgebra(FactoredAlgebra):
    def lift(self,value,order=6):return super().lift(value,order)


class NestedTaylor:
    """Finite radial Taylor coefficients which may themselves be factored."""
    def __init__(self,rows):self.rows=list(rows);self.ctx=rows[0].ctx;self.order=len(rows)-1
    def coerce(self,value):
        if isinstance(value,NestedTaylor):return value
        row=self.rows[0]*0+value
        return NestedTaylor([row]+[row*0]*self.order)
    def __add__(self,value):
        other=self.coerce(value);return NestedTaylor([a+b for a,b in zip(self.rows,other.rows)])
    __radd__=__add__
    def __neg__(self):return NestedTaylor([-row for row in self.rows])
    def __sub__(self,value):return self+-self.coerce(value)
    def __rsub__(self,value):return self.coerce(value)+-self
    def __mul__(self,value):
        other=self.coerce(value)
        return NestedTaylor([sum((self.rows[j]*other.rows[k-j] for j in range(k+1)),self.rows[0]*0)
            for k in range(self.order+1)])
    __rmul__=__mul__
    def reciprocal(self):
        out=[self.rows[0].reciprocal()]
        for k in range(1,self.order+1):out.append(-out[0]*sum((self.rows[j]*out[k-j] for j in range(1,k+1)),out[0]*0))
        return NestedTaylor(out)
    def __truediv__(self,value):return self*self.coerce(value).reciprocal()
    def Zderivative(self):return NestedTaylor([row.Zderivative() if isinstance(row,FactoredJet) else derivative(row) for row in self.rows])
    def ordinary(self):return [row*math.factorial(k) for k,row in enumerate(self.rows)]


def radial_core_log(phi_rows,rho):
    """(rho*D_rho)^k logPhi, k1..3, true axial6; no radial sampling."""
    phi=NestedTaylor([row/math.factorial(k) for k,row in enumerate(phi_rows)])
    # Differentiate the outer series, then divide as a series through2.
    top=NestedTaylor([(k+1)*phi.rows[k+1] for k in range(3)])
    bottom=NestedTaylor(phi.rows[:3]);quotient=(top/bottom).ordinary()
    return [sum((quotient[i-1]*(S2[k][i]*rho**i) for i in range(1,k+1)),quotient[0]*0) for k in range(1,4)]


def smoothed_comparison_rows(algebra,phi,v,core_log_y,core_v_y,alpha):
    """Original alpha interpolation; not a frozen radial field in smoothing."""
    logcore=[algebra.width(row,k) for k,row in enumerate(core_log_y,1)]
    vcore=[algebra.width(row,k) for k,row in enumerate(core_v_y,1)]
    logbar=[sum((logcore[k-j]*alpha[j]*math.comb(k,j) for j in range(k+1)),logcore[0]*0) for k in range(3)]
    vbar=[algebra.lift(v)]+[sum((vcore[k-j]*alpha[j]*math.comb(k,j) for j in range(k+1)),vcore[0]*0) for k in range(3)]
    bell=exponential_derivatives(logbar+[logbar[0]*0])
    return [algebra.lift(phi)*row for row in bell[:4]],vbar,logbar


def comparison_moment_rows(algebra,coordinate_scale,phi,V,moments):
    """Ordinary coordinate derivatives of own comparison moment shapes."""
    lift=algebra.lift;scale=lift(coordinate_scale)
    initials={MTH:moments[MTH],MZ:moments[MZ],MTHZ:moments[MTHZ],
        'A':moments[MZT]['axial'],'B':moments[MZT]['swirl'],MP:moments[MP]}
    rhs={MTH:[row*2 for row in phi],MZ:V,MTHZ:[product_rows(phi,V,k)*2 for k in range(4)],
        'A':[product_rows(V,V,k) for k in range(4)],'B':[product_rows(phi,phi,k) for k in range(4)],
        MP:[product_rows(phi,phi,k) for k in range(4)]}
    rates={MTH:2,MZ:1,MTHZ:2,'A':1,'B':2,MP:1};out={}
    for name,value in initials.items():
        rows=[lift(value)]
        for k in range(3):rows.append(scale*(rhs[name][k]-rows[k]*rates[name]))
        out[name]=rows
    return out


def comparison_directions(algebra,Z,delta,phi,V,moments,p0,ratios,ratios2,coordinate_scale):
    """General (9.13) with varying comparison Fbar/Vbar and own histories."""
    rows=comparison_moment_rows(algebra,coordinate_scale,phi,V,moments)
    nested=lambda ordinary:NestedTaylor([row/math.factorial(k) for k,row in enumerate(ordinary)])
    lift=lambda value:nested([algebra.lift(value)]+[algebra.lift(value)*0]*3)
    true=lambda ordinary,factors:nested([row*IntervalTaylor(algebra.ctx,
        [factors[k]/math.factorial(k) for k in range(7)]) for row in ordinary])
    z=lift(IntervalTaylor.variable(algebra.ctx,algebra.ctx.mpf(Z),6));d=1-z*z;L=1-z*z*delta
    f=true(phi,ratios);h=true(rows[MTH],ratios);k=true(rows[MTHZ],ratios)
    b=true(rows['B'],ratios2);p=true(rows[MP],ratios2);m=nested(rows[MZ]);a=nested(rows['A']);vv=nested(V);pressure=lift(p0)
    W=1-(z*m)*(1-delta)-d*m.Zderivative()
    angular=h*(1-delta/2)-(z*h.Zderivative())*((1-delta)/2)-d*k.Zderivative()+(z*k)*(2*delta-1)
    values=dict(D_over_R=(-W+angular/(f*2))/L,
        hydro=(-W*vv+(m-z*m.Zderivative())*((1-delta)/2)+(z*a)*(2*delta)-d*a.Zderivative())/(L*2),
        pressure=((z*pressure)*(2*(1+delta))-d*pressure.Zderivative())/(L*2),
        swirl=(-(z*b)*(2*delta)+d*b.Zderivative()+(z*p)*(2*(1+delta))-d*p.Zderivative())/(L*2))
    return {name:row.ordinary() for name,row in values.items()},rows


def bridge_controls(algebra,coordinate_scale,chi,Dbar,drive,quotient,bar_log):
    """Differentiate ORIGINAL chi*Dbar and chi*quotient*axial drive."""
    scale=algebra.lift(coordinate_scale)
    logF=[-scale*product_rows(chi,Dbar,k)/2 for k in range(4)]
    relative_log=[logF[k]-bar_log[k] for k in range(3)]+[logF[3]*0]
    quotient_rows=[algebra.lift(quotient)*row for row in exponential_derivatives(relative_log)]
    V=[]
    for k in range(4):
        value=quotient_rows[0]*0
        for i in range(k+1):
            for j in range(k-i+1):
                ell=k-i-j;weight=math.factorial(k)//(math.factorial(i)*math.factorial(j)*math.factorial(ell))
                value+=chi[i]*quotient_rows[j]*drive[ell]*weight
        V.append(-scale*value)
    return dict(logF_coordinate_derivatives=logF,
        logUtheta_coordinate_derivatives=[logF[0]+scale/2]+logF[1:],
        Uz_positive_order_coordinate_derivatives=V,quotient_coordinate_derivatives=quotient_rows[:4])


def export_logR_ledger(c,packet,algebra,width_conversion):
    rows={row['physical_row']:row for row in packet['final_factored_physical_row_ledgers']};ledger={}
    for group in ('physical_velocity_pressure_phase_Z_mixed4','physical_five_primitive_phase_Z_mixed4'):
        ledger[group]={}
        for name,grid in packet[group].items():
            ledger[group][name]={}
            for key,value in grid.items():
                k=int(key.split('_')[0][1:]);logs=[]
                for term in rows[name+'/'+key]['terms']:
                    size=max(abs(v) for v in endpoints(term['final_ordinary_coefficient']))
                    if not size:continue
                    powers=list(term['source_exponents']);powers[0]-=k if width_conversion else 0
                    logs.append(sum((algebra.logs[i]*p for i,p in enumerate(powers) if p),c.mpf(0))+c.ln(c.mpf(size)))
                ledger[group][name][key.replace('s','y',1)]=dict(exact_zero=not logs,nonzero_factored_terms=len(logs),
                    log_absolute_upper=None if not logs else c.mpf(max(endpoints(row)[1] for row in logs))+c.ln(len(logs)),
                    source_coordinate_order=k,exact_source_conversion='D_y^k=hb^-k*D_s^k' if width_conversion else 'D_y^k=D_coordinate^k',
                    formal_width_powers_cancelled_before_log_bound=True)
    return ledger


class CompliantBridgeMixedC4:
    def __init__(self):
        self.micro=CompliantMicroswitchMixedC4();self.bridge=self.micro.switch.bridge;self.core=self.bridge.core
        self.ctx=c=self.bridge.ctx;self.family=self.bridge.family;self.source=self.bridge.source
        self.hashes=dict(self.micro.hashes);self.proofs=[];self.core_bounds=[];self.cache={}
        name=PREFIX+'microswitch_mixed_C4_check.json';record=json.loads((HERE/name).read_bytes())
        if not record['all_passed'] or record['actual_five_defect_family_sha256']!=self.family or record['implicit_source_sha256']!=self.source:
            raise ValueError('Accepted SAME-source microswitch mixed4 required')
        for path,digest in record['input_hashes'].items():
            if hashlib.sha256((HERE/path).read_bytes()).hexdigest()!=digest:raise ValueError('Bridge mixed source changed: '+path)
        self.hashes.update(record['input_hashes']);self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.logh=self.micro.logh;self.h=self.micro.h;self.logF0=self.micro.logF0;self.r=self.bridge.r
        linear=self.core.records['shared_linear_resolvent'];major=self.core.records['core_transfer']
        self.phi_norm=read_interval(c,linear['Phi_model_Xh_norm_upper'])+self.core.correction
        self.psi_norm=read_interval(c,major['fresh_Psi_model_Xh_norm_upper'])+self.core.correction
        self.shared=self.micro.shared_axial_source
        self.core_norm_bindings={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in (
            'lei_ren_part1_paper_shared_linear_resolvent.json',PREFIX+'core_transfer.json',
            PREFIX+'core_physical_field.json','lei_ren_part1_paper_shared_analytic_tube.json')}
        self.hashes.update(self.core_norm_bindings)
        self.comparison_source=self.formal_comparison_source()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def formal_comparison_source(self):
        """Exact global functions, separately from the broad numeric covers."""
        fields=dict(logphi=dict(initial='log(Phi_core(4,Z))',
            integrand='alpha(t)*D_t log(Phi_core(4exp(t),Z))',bounds=['0','y']),
            V=dict(initial='V_core(4,Z)',integrand='alpha(t)*D_t V_core(4exp(t),Z)',bounds=['0','y']))
        initial=dict(H='(1/8)*integral_0^4 rho*Phi_core(rho,Z)drho',
            M='(1/4)*integral_0^4 V_core(rho,Z)drho',K='(1/8)*integral_0^4 rho*Phi_core(rho,Z)*V_core(rho,Z)drho',
            A='(1/4)*integral_0^4 V_core(rho,Z)^2 drho',B='(1/16)*integral_0^4 rho*Phi_core(rho,Z)^2 drho',
            C='(1/4)*integral_0^4 Phi_core(rho,Z)^2 drho')
        moments={name:dict(initial=initial[name],initial_weight='exp(-lambda*y)',
            kernel='exp(-lambda*(y-t))',rate=rate,integrand=integrand,bounds=['0','y'])
            for name,rate,integrand in [('H',2,'2*phi_bar(t,Z)'),('M',1,'V_bar(t,Z)'),
                ('K',2,'2*phi_bar(t,Z)*V_bar(t,Z)'),('A',1,'V_bar(t,Z)^2'),
                ('B',2,'phi_bar(t,Z)^2'),('C',1,'phi_bar(t,Z)^2')]}
        recipe=dict(core_analytic_family=self.core.records['core_transfer']['analytic_core_family_sha256'],
            implicit_source_sha256=self.source,core_norm_source_bindings=self.core_norm_bindings,
            alpha='1-sigma((t-hb)/hb)',core_radius='rho=4exp(t)',fields=fields,moments=moments)
        namespace=hashlib.sha256(json.dumps(recipe,sort_keys=True).encode('utf8')).hexdigest()
        return dict(**recipe,source_namespace=namespace,
            moment_ODEs=['H_y+2H=2phi_bar','M_y+M=V_bar','K_y+2K=2phi_bar*V_bar',
                'A_y+A=V_bar^2','B_y+2B=phi_bar^2','C_y+C=phi_bar^2'],
            numeric_covers_are_not_exact_point_values=True)

    def core_rectangular(self,rho,Z):
        """Ordinary rho0..3/Z0..6, old total5 plus admitted full-norm bounds."""
        c=self.ctx;rho=c.mpf(rho);Z=c.mpf(Z);key=(rho._mpi_,Z._mpi_)
        if key in self.cache:return self.cache[key]
        packet=self.core.normalized_jets(rho,Z);rmax=c.mpf(endpoints(rho)[1]);result={name:[] for name in ('Phi','Uz')}
        for i in range(4):
            for name in result:
                coefficients=[]
                for n in range(7):
                    if i+n<=5:value=packet[name][gridkey(i,n)]
                    else:
                        norm=self.phi_norm if name=='Phi' else self.psi_norm*self.core.epsilon
                        bound=norm*embedding(c,self.core.h,rmax,i,n);upper=endpoints(bound)[1]
                        value=c.mpf([-upper,upper])
                        self.core_bounds.append(dict(component=name,rho_order=i,axial_order=n,rho_max=rmax,
                            actual_analytic_Xh_norm=norm,embedding_factor=embedding(c,self.core.h,rmax,i,n),
                            ordinary_derivative_absolute_upper=bound,same_unique_analytic_core=True))
                    coefficients.append(value/math.factorial(n))
                result[name].append(IntervalTaylor(c,coefficients))
        self.cache[key]=result;return result

    def source_graph(self,coordinate):
        prefix=dict(formal_integral='I_bridge',shared_source=self.shared['shared_source_namespace'],bounds=['0',coordinate],
            source_integrand=self.shared['formal_signed_integrals']['I_bridge'])
        return dict(op='sum',args=self.shared['V100']['args'][:-1]+[prefix])

    def evaluate(self,Z,value,chart):
        c=self.ctx;Z=c.mpf(Z);value=c.mpf(value);lo,hi=endpoints(value)
        if chart not in ('first','second','macro'):raise ValueError('Original bridge chart required')
        limits=(0,1) if chart in ('first','macro') else (1,2)
        if lo<limits[0] or hi>limits[1]:raise ValueError('Original bridge coordinate interval required')
        microscopic=chart!='macro'
        if microscopic:
            y=self.h*value;R=self.r*c.exp(y);theta=c.exp(-y)
            if lo==hi==0:R=self.r;theta=c.mpf(1)
            formal_y='hb*s';formal_radius=dict(op='mul',args=['Ra',dict(op='exp',log='hb*s')])
        else:
            total=c.ln(100/self.r);length=total-2*self.h
            if endpoints(length)[0]<=0:raise ValueError('Frozen bridge must precede R100')
            y=2*self.h+length*value
            if lo==hi==1:y=total;R=c.mpf(100);theta=self.r/100
            else:
                R=self.r*c.exp(y);theta=c.exp(-y)
                R=c.mpf([max(endpoints(self.r)[0],endpoints(R)[0]),min(mp.mpf(100),endpoints(R)[1])])
                theta=c.mpf([max(endpoints(self.r/100)[0],endpoints(theta)[0]),min(mp.mpf(1),endpoints(theta)[1])])
            formal_y='2hb+fraction*(log(100/Ra)-2hb)'
            formal_radius=dict(op='mul',args=['Ra',dict(op='exp',log=formal_y)])
        theta=c.mpf([max(endpoints(self.r/100)[0],endpoints(theta)[0]),min(mp.mpf(1),endpoints(theta)[1])])
        actual=self.bridge.actual(Z,theta);inp=self.bridge.inputs(Z);comparison=self.bridge.comparison(Z,theta)
        jet=lambda row:IntervalTaylor(c,row)
        phi=jet(actual['F_actual_over_F0_axial5_coefficients']);V=jet(actual['Uz_actual_axial5_coefficients'])
        logu0=c.ln(2*R)/2+self.logF0+c.ln(phi[0])-self.core.logP
        algebra=AxialSixAlgebra(c,[self.logh,2*self.core.logP,2*self.logF0,2*logu0],self.proofs)
        h=algebra.width(1);scale=h if microscopic else algebra.lift(1)
        if microscopic:
            rho=c.mpf(4)*c.exp(y)
            if endpoints(rho)[1]>endpoints(c.mpf('4.1'))[1]:raise ValueError('Original smoothing leaves analytic core continuation')
            core=self.core_rectangular(rho,Z);logcore=radial_core_log(core['Phi'],rho)
            vcore=[sum((core['Uz'][i]*(S2[k][i]*rho**i) for i in range(1,k+1)),core['Uz'][0]*0) for k in range(1,4)]
            if chart=='first':
                alpha=[algebra.lift(1)]+[algebra.lift(0)]*3
                sigma=[row*math.factorial(k) for k,row in enumerate(sigma_jets(c,value))]
                chi=[algebra.lift(1-sigma[0])+h*sigma[0]]+[-algebra.lift(sigma[k])+h*sigma[k] for k in range(1,4)]
            else:
                sigma=[row*math.factorial(k) for k,row in enumerate(sigma_jets(c,value-1))]
                alpha=[algebra.lift(1-sigma[0])]+[-algebra.lift(sigma[k]) for k in range(1,4)]
                chi=[h]+[h*0]*3
            barphi,barV,barlog=smoothed_comparison_rows(algebra,comparison['phi'],comparison['v'],logcore,vcore,alpha)
        else:
            rho=None;barphi=[algebra.lift(comparison['phi'])]+[algebra.lift(0)]*3
            barV=[algebra.lift(comparison['v'])]+[algebra.lift(0)]*3;barlog=[algebra.lift(0)]*3
            alpha=[algebra.lift(0)]*4;chi=[h]+[h*0]*3
        directions,momentrows=comparison_directions(algebra,Z,self.core.delta,barphi,barV,comparison['moments'],inp['p0'],
            inp['F0_ratios'],inp['F0_squared_ratios'],scale)
        Dbar=[rate_rows(directions['D_over_R'],scale,k)*R for k in range(4)]
        drive=[rate_rows(directions['hydro'],scale,k)*R
            +algebra.shift(rate_rows(directions['pressure'],scale,k)*R,(0,1,0,0))
            +algebra.shift(rate_rows(directions['swirl'],scale*2,k)*R**2,(0,0,1,0)) for k in range(4)]
        controls=bridge_controls(algebra,scale,chi,Dbar,drive,phi/comparison['phi'].truncate(5),barlog)
        vrows=[V]+controls['Uz_positive_order_coordinate_derivatives']
        dressed=dress(phi,inp['F0_ratios']);amp=dressed/phi[0]
        amp=IntervalTaylor(c,[c.mpf(1)]+list(amp.coefficients[1:]))
        raw=actual['actual_moment_shape_axial5_coefficients']
        shapes=dict(theta=jet(raw[MTH])/(phi*2),theta_z=jet(raw[MTHZ])/(phi*2),mean=jet(raw[MZ]),
            axial=jet(raw[MZT]['axial']),swirl=jet(raw[MZT]['swirl'])/square(phi),pressure=jet(raw[MP])/square(phi))
        coordinate_power=(algebra.width if microscopic else lambda row,power=1:algebra.lift(row))
        packet=phase_physical(c,Z,self.core.delta,scale,amp,logu0,controls['logUtheta_coordinate_derivatives'],vrows,shapes,
            jet(actual['pressure_axis_axial5_coefficients']),algebra.shift(1,(0,-1,0,0)),coordinate_power,self.proofs)
        logledger=export_logR_ledger(c,packet,algebra,microscopic)
        if not microscopic:
            for group in ('physical_velocity_pressure_phase_Z_mixed4','physical_five_primitive_phase_Z_mixed4'):
                target=group.replace('phase','y')
                packet[target]={name:{key.replace('s','y',1):row for key,row in grid.items()} for name,grid in packet.pop(group).items()}
                logledger[target]=logledger.pop(group)
            for row in packet['final_factored_physical_row_ledgers']:row['physical_row']=row['physical_row'].replace('/s','/y')
            packet['actual_Q_y_derivative_axial4']=packet.pop('actual_Q_phase_derivative_axial4')
        return dict(Z=Z,coordinate=value,chart=chart,**packet,
            physical_derivative_coordinate='s=log(R/Ra)/hb' if microscopic else 'y=logR; fraction labels coverage only',
            controls=algebra.tree(controls),comparison_directions=algebra.tree(directions),
            comparison_logphi_positive_coordinate_orders=algebra.tree(barlog),
            comparison_own_moment_coordinate_rows=algebra.tree(momentrows),
            actual_parent_axial5_packet=actual,comparison_parent_axial6_packet=dict(phi=list(comparison['phi'].coefficients),
                V=list(comparison['v'].coefficients)),
            exact_comparison_function_source=dict(namespace=self.comparison_source['source_namespace'],argument_y=formal_y,
                fields=['logphi','V'],own_moments=['H','M','K','A','B','C']),
            physical_logR_Z_mixed4_log_bound_ledger=logledger,factored_source_log_bases=list(algebra.logs),
            source_width_log=self.logh,width_enclosure_is_not_source=self.h,
            source_radius_tree=formal_radius,source_logR_over_Ra=formal_y,
            y_enclosure_only=y,R_enclosure_only=R,theta_enclosure_only=theta,core_rho_enclosure_only=rho,
            source_chi='1-(1-hb)*sigma(y/hb)',source_alpha='1-sigma((y-hb)/hb)',
            actual_axial_function_source=self.source_graph(formal_y),
            angular_source='F=F_exit*exp[-.5*integral_0^y chi(t)*Dbar(t,Z)dt]',
            actual_moments_not_replaced_by_comparison=True,comparison_radially_frozen=not microscopic,
            all_source_width_and_amplitude_factors_retained_until_final_rows=True,
            bridge_mixed4_available=True,core_bridge_and_R100_local_functional_joins_certified=True,
            full_inner_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False)

    def report(self):
        packets=dict(whole_first=self.evaluate([-1,1],[0,1],'first'),whole_second=self.evaluate([-1,1],[1,2],'second'),
            core_exit=self.evaluate([-1,1],0,'first'),first_side_phase1=self.evaluate([-1,1],1,'first'),
            second_side_phase1=self.evaluate([-1,1],1,'second'),smoothing_exit=self.evaluate([-1,1],2,'second'),
            whole_macro=self.evaluate([-1,1],[0,1],'macro'),macro_inlet=self.evaluate([-1,1],0,'macro'),
            R100_exit=self.evaluate([-1,1],1,'macro'))
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,**packets,
            interior_packets=[self.evaluate(z,v,chart) for z,v,chart in (('0','.5','first'),('.5','1.5','second'),('0','.5','macro'))],
            same_core_rectangular_extension_proofs=self.core_bounds,source_phi_Xh_norm=self.phi_norm,source_psi_Xh_norm=self.psi_norm,
            factored_positive_width_amplitude_cap_proofs=self.proofs,shared_exact_axial_source=self.shared,
            exact_global_comparison_source=self.comparison_source,full_core_Xh_norm_source_bindings=self.core_norm_bindings,
            Psi_full_norm_uses_fresh_transfer_not_uncertified_linear_Psi_norm=True,
            smoothing_continuation_log_distance=self.micro.continuation_distance,
            original_core_stress_free_source='same admitted analytic fixed point solving (3.20)/(8.2); on alpha=1 the comparison is the exact core',
            core_bridge_join_source='chi(0)=1 with all positive derivatives zero; comparison=core; same actual exit moments/P0; stress-free core ODEs give identical jets',
            phase1_phase2_join_source='original flat sigma endpoint jets; same own comparison and actual moments/P0; alpha/chi share all one-sided jets',
            R100_join_source='late chi=hb; frozen comparison and same actual R100 histories; first switch D_s=hb*D_y and V_s=hb*V_y at its inlet',
            bridge_mixed4_available=True,core_bridge_and_R100_local_functional_joins_certified=True,
            full_inner_interfaces_certified=False,full_cartesian_vector_derivatives_certified=False,
            admissible_stress_lift_constructed=False,temporal_recursion=False,input_hashes=self.hashes)


def run():
    with mp.workdps(280):result=CompliantBridgeMixedC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf8')
    print('Original varying-comparison bridge three-chart mixed4 and factored logR ledgers generated',flush=True)
    return result


if __name__=='__main__':run()
