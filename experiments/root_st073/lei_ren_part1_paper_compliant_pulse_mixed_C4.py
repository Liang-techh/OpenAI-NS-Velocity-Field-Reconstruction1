"""Same-source O.4 velocities and pressure at every mixed order <= four.

y=log(R) is the paper similarity radial coordinate. Physical velocity
prefactors are differentiated before axial derivatives; this is not a
cylindrical-r or physical-z derivative certificate. Global outer C4, cone,
physical energy and temporal coefficient recursion remain separate tasks.
"""
import hashlib
import json
import math
from pathlib import Path

import mpmath as mp
from lei_ren_part1_paper_compliant_pulse_radial_C4 import CompliantPulseRadialC4
from lei_ren_part1_paper_compliant_pulse_high_jets import derivative
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import FlatPulseDerivatives, gp_jets
from lei_ren_part1_paper_compliant_five_moment_repair import pack
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_'


def binomial_product(rows,rate,k):
    result=rows[0]*0
    for j in range(k+1):result+=rows[j]*(math.comb(k,j)*rate**(k-j))
    return result


def transport_mixed(field,Z,point,Brows,u):
    """Repeated source ODE differentiation and true physical product rules."""
    c=field.ctx; mu=field.mu; bp=c.mpf('.5')+mu
    if len(Brows)!=5 or any(v.order!=5 for v in Brows):raise ValueError('True C5 axial source and y orders0..4 required')
    primitive={
        'Mz_over_R_Utheta':[point['Mz_over_R_Utheta']],
        'Mtheta_z_over_sqrt2_R_3half_Utheta_squared':[point['Mtheta_z_over_sqrt2_R_3half_Utheta_squared']],
        'Mtheta_over_sqrt2_R_3half_Utheta':[IntervalTaylor.constant(c,point['Mtheta_over_sqrt2_R_3half_Utheta'],5)
            if not isinstance(point['Mtheta_over_sqrt2_R_3half_Utheta'],IntervalTaylor)
            else point['Mtheta_over_sqrt2_R_3half_Utheta']],
        'Mztheta_over_R_Utheta_squared':[point['Mztheta_over_R_Utheta_squared']],
        'P_over_Pstar_squared':[point['pressure']['P_over_Pstar_squared']]}
    m1=primitive['Mz_over_R_Utheta']; m2=primitive['Mtheta_z_over_sqrt2_R_3half_Utheta_squared']
    X=primitive['Mtheta_over_sqrt2_R_3half_Utheta']; energy=primitive['Mztheta_over_R_Utheta_squared']
    pressure=primitive['P_over_Pstar_squared']; one=IntervalTaylor.constant(c,1,5)
    for k in range(4):
        m1.append(Brows[k]-m1[k]*(c.mpf('.5')-mu))
        m2.append(Brows[k]-m2[k]*(c.mpf('.5')-2*mu))
        X.append((one if k==0 else one*0)-X[k]*(1-mu))
        square=Brows[0]*0
        for j in range(k+1):square+=Brows[j]*Brows[k-j]*math.comb(k,j)
        energy.append(square-(one/2 if k==0 else one*0)+energy[k]*(2*mu))
        pressure.append(point['pressure']['P_y_over_Pstar_squared']*(-field.prate)**k)
    A=[field.radial(Z,b,m) for b,m in zip(Brows,m1)]
    physical={
        'Uz_over_Pstar_without_common_theta_radial_factor':[u*binomial_product(Brows,-bp,k) for k in range(5)],
        'Utheta_over_Pstar_without_common_theta_radial_factor':[u*(-bp)**k for k in range(5)],
        'Ur_over_sqrt_R_over_2_Pstar_without_common_theta_radial_factor':[
            u.truncate(4)*binomial_product(A,-mu,k) for k in range(5)],
        'P_over_Pstar_squared':pressure}
    # Store each derivative as an ordinary axial Taylor jet through 4-k,
    # and also expose all 15 actual mixed derivatives, not Taylor coefficients.
    grid={}
    for label,rows in physical.items():
        grid[label]={}
        for k,row in enumerate(rows):
            for n in range(5-k):grid[label]['y'+str(k)+'_Z'+str(n)]=row[n]*math.factorial(n)
    return dict(B_y_derivative_Taylor=Brows,primitive_y_derivative_Taylor=primitive,
        physical_velocity_and_pressure_y_derivative_Taylor={label:[row.truncate(4-k) for k,row in enumerate(rows)]
            for label,rows in physical.items()},physical_mixed_derivatives_total_order_le4=grid,
        normalized_radial_recovery_y_derivative_Taylor=A,
        physical_radial_y_derivative_Taylor=[binomial_product(A,-mu,k) for k in range(5)])


class CompliantPulseMixedC4(CompliantPulseRadialC4):
    def __init__(self):
        super().__init__(); self.flat=FlatPulseDerivatives(self.ctx)
        name=PREFIX+'compliant_flat_pulse_derivatives_check.json'; r=json.loads((HERE/name).read_bytes())
        if not r['all_passed'] or not r['original_radial_shape_derivatives_C4_available'] or r['actual_five_defect_family_sha256']!=self.flat.family:
            raise ValueError('Same-family admitted flat derivatives required')
        for source,digest in r['input_hashes'].items():
            if hashlib.sha256((HERE/source).read_bytes()).hexdigest()!=digest:raise ValueError('Mixed C4 source changed: '+source)
            self.hashes[source]=digest
        self.hashes[name]=hashlib.sha256((HERE/name).read_bytes()).hexdigest()
        self.hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()

    def _high_packet(self,point):
        point=super()._high_packet(point); c=self.ctx; Z=point['Z']
        selected,_,u,_,_=self.data(Z); coord=point['coordinate']; kind=coord['kind']
        if kind=='main':
            shape=gp_jets(c,coord['xi']); ap=selected['selected_ap_Taylor']
            Brows=[ap*(shape[k]*math.factorial(k)*self.mu**k) for k in range(5)]
        elif kind=='end':
            Brows=[selected['selected_ap_Taylor']*0 for _ in range(5)]
            for Cj,center in zip(selected['selected_scaled_end_coefficient_Taylor'],(-3,-1)):
                beta=self.flat.beta(coord['offset_from_Rv']-center)
                for k in range(5):Brows[k]+=Cj*(beta[k]*math.factorial(k))
            Brows=[b*self.Ecap for b in Brows]
        elif kind in ('gap_main','gap_end'):Brows=[selected['selected_ap_Taylor']*0 for _ in range(5)]
        else:raise ValueError('Unknown actual O.4 chart: '+kind)
        mixed=transport_mixed(self,Z,point,Brows,u)
        radial=mixed['physical_radial_y_derivative_Taylor']
        thetaZ=derivative(u)/u.truncate(4)
        point.update(mixed)
        point.update(Uz_over_Utheta=Brows[0],Uz_y_over_Utheta=binomial_product(Brows,-(c.mpf('.5')+self.mu),1),
            Ur_over_sqrt_R_over_2_Utheta=radial[0],Ur_y_over_sqrt_R_over_2_Utheta=radial[1],
            Ur_Z_over_sqrt_R_over_2_Utheta=derivative(radial[0])+radial[0].truncate(3)*thetaZ.truncate(3),
            Ur_yZ_over_sqrt_R_over_2_Utheta=derivative(radial[1])+radial[1].truncate(3)*thetaZ.truncate(3),
            factored_Uz_over_Pstar_Taylor=u*Brows[0],
            factored_Ur_over_sqrt_R_over_2_Pstar_Taylor=u.truncate(4)*radial[0],
            pulse_all_mixed_derivatives_total_order_le4_available=True,
            original_flat_pulse_support_derivatives_installed=True,
            mixed_coordinate_scope='y=log similarity R and paper axial Z; physical velocity components with differentiated radial prefactors',
            cartesian_spatial_derivatives_certified=False,
            full_pulse_C4_installed=False,full_outer_C4_certified=False,
            uniform_pulse_C4_chart_interface_certificate_available=False,
            next_dependency='Uniform all-chart mixed bounds and quantitative source-identified interfaces; post-pulse mixed C4 and stress cone')
        return point

    def report(self):
        with mp.workdps(210):
            r=super().report()
            # The lower endpoint of the reciprocal interval yields a
            # legal real s box. A small main-chart overlap covers its
            # possible displacement from the formal s=-1/mu boundary.
            gap_end_start=-endpoints(1/self.mu)[0]
            c=self.ctx; overlap_xi=13+self.mu*c.mpf(gap_end_start)
            if endpoints(overlap_xi)[1]>mp.mpf('12.0001'):
                raise ArithmeticError('Whole gap chart coverage not proved')
            r.update(whole_Z_entrance_box=self.main([-1,1],[0,'.02'],cells=64),
                whole_Z_main_box=self.main([-1,1],['.02','10'],cells=64),
                whole_Z_exit_box=self.main([-1,1],[10,11],cells=64),
                whole_Z_gap_main_box=self.gap([-1,1],[11,'12.0001']),
                whole_Z_gap_end_box=self.gap_from_end([-1,1],[gap_end_start,-4]),
                whole_Z_end_box=self.end([-1,1],[-4,0],cells=64),
                whole_gap_box_overlap_xi=overlap_xi,
                pulse_all_mixed_derivatives_total_order_le4_available=True,
                original_flat_pulse_support_derivatives_installed=True,
                factored_all_chart_mixed_derivative_boxes_available=True,
                all_chart_box_scope='finite factored interval bounds; physical norm ledger and functional interface proof still pending',
                full_pulse_C4_installed=False,full_outer_C4_certified=False,
                uniform_pulse_C4_chart_interface_certificate_available=False,
                next_dependency='Uniform all-chart mixed bounds/flat interfaces, then post-pulse C4 and admissible stress; genuine temporal recursion remains pending')
            return r


def run():
    result=CompliantPulseMixedC4().report()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Actual all-chart O.4 velocity/pressure mixed derivatives through total order4 generated; uniform outer C4/cone pending',flush=True)
    return result


if __name__=='__main__':run()
