"""Same-source directed switched comparison and scaled inlet transfer.

All amplitude/Lambda factors cancel algebraically before transfer. No F0
point amplitude is required. Radial profile comparison is not temporal
coefficient recursion or a repaired reference inlet.
"""
import hashlib
import json
import math
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_logarithmic_core_exit import LogarithmicCoreExit
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_candidate_shared_inlet import diff
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode

HERE=Path(__file__).parent

def scaled_transfer(c,phi,u,m,ell_scaled,S_scaled,p0_scaled,eps,s,z,dt,initial_phi):
    """Same normalized_inlet equations with e ell, e^2 S, e P0 inputs.

    I_z is physical I_z from the frozen comparison convention. The
    base driver B (chi=1) cancels its sqrt(e) factors exactly. Physical
    exit consumers MUST multiply both base drivers by their cutoff chi.
    """
    d=1-z*z;L=1-z*z*dt
    if endpoints(L[0])[0]<=0 or endpoints(phi[0])[0]<=0:raise ValueError('inlet denominator crosses zero')
    pressure=p0_scaled+S_scaled*m['p']
    mzts=m['u_squared']*eps-S_scaled*m['weighted_phi_squared']
    V=z*m['z']*(1-dt)+d*diff(m['z'])-s
    baseC=m['theta']*(1-dt/2)-z*diff(m['theta'])*((1-dt)/2)
    baseC-=d*diff(m['theta_z']);baseC+=z*m['theta_z']*(2*dt-1)
    Cscaled=baseC*eps-ell_scaled*(z*m['theta']*((1-dt)/2)+d*m['theta_z'])
    itheta=phi*V*eps/L+Cscaled/(2*L*s)
    D=itheta/phi
    Ns=(V*u+(m['z']-z*diff(m['z']))*((1-dt)/2))*eps
    Ns+=z*mzts*(2*dt)-d*diff(mzts)
    Ns+=(z*pressure*(2*(1+dt))-d*diff(pressure))*s
    B=-(initial_phi/phi)*Ns/(2*L)
    Iz=Ns/(L*c.sqrt(2*s)*c.sqrt(eps))
    Ur=(z*u*(2*s)-z*m['z']*(1-dt)-d*diff(m['z']))*c.sqrt(eps)/(L*c.sqrt(2*s))
    return dict(D=D,I_z=Iz,pressure_scaled=pressure,scaled_axial_numerator=Ns,
        unmodulated_driver_A=-D/2,unmodulated_driver_B=B,physical_Ur=Ur,
        angular_transfer=itheta,combined_fifth_moment_scaled=mzts)

class LogarithmicComparison:
    def __init__(self):
        self.core=LogarithmicCoreExit();self.ctx=c=self.core.ctx;self.precision=c.dps
        self.identity=self.core.state['identity'];self.center_family=self.identity['center_Z']
        with mp.workdps(c.dps+40):
            self.eps=read_interval(c,self.identity['epsilon']);self.delta=read_interval(c,self.identity['delta'])
            self.z=IntervalTaylor.variable(c,c.mpf(self.center_family),2)
            jet=lambda row:IntervalTaylor(c,row[:3])
            self.ell_scaled=jet(self.core.fixed['ell_Z_taylor'])
            self.S_scaled=jet(self.core.fixed['S_Z_taylor'])
            self.p0_scaled=jet(self.core.fixed['P0_Z_taylor'])
            self.initial_phi=self.field(c.mpf(4),'Phi')
        names=(Path(__file__).name,'lei_ren_part1_paper_interval_comparison_enclosure.py',
            'lei_ren_part1_paper_interval_exit_continuation_enclosure.py')
        self.input_hashes={**self.core.input_hashes,**{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}}

    def field(self,r,name,radial_order=0):
        c=self.ctx
        return IntervalTaylor(c,[self.core.mixed(r,radial_order,k)['full_'+name]/math.factorial(k)
            for k in range(3)])

    def core_state(self,r):
        result=self.core.core_moments(r)
        return dict(phi=self.field(r,'Phi'),U=self.field(r,'Uz'),
            **{n:IntervalTaylor(self.ctx,row) for n,row in result['full_moment_coefficients'].items()})

    def transfer(self,state,s):
        return scaled_transfer(self.ctx,state['phi'],state['U'],state,self.ell_scaled,
            self.S_scaled,self.p0_scaled,self.eps,s,self.z,self.delta,self.initial_phi)

    def integrate(self,cells=32):
        if not isinstance(cells,int) or cells<1:raise ValueError('positive cell count required')
        c=self.ctx
        with mp.workdps(c.dps+40):
            hb=c.mpf('.005');step=hb/cells
            start=self.core_state(4*c.exp(hb));zero=start['phi']*0
            G,J=zero,zero
            moments={n:v for n,v in start.items() if n not in ('phi','U')}
            saved=[];minimum=None
            for n in range(cells):
                left=hb*(1+c.mpf(n)/cells);right=hb*(1+c.mpf(n+1)/cells)
                y=c.mpf([endpoints(left)[0],endpoints(right)[1]]);s=4*c.exp(y)
                source=self.field(s,'Phi')
                lower=endpoints(source[0])[0]
                if lower<=0:raise ValueError('core source denominator crosses zero')
                minimum=lower if minimum is None else min(minimum,lower)
                alpha=alpha_box(c,y,hb)
                g=self.field(s,'Phi',1)*s*alpha/source
                j=self.field(s,'Uz',1)*s*alpha
                partial=c.mpf([0,endpoints(step)[1]])
                f=start['phi']*(G+g*partial).exp();u=start['U']+J+j*partial
                rhs=dict(theta=f*(2*s*s),z=u*s,theta_z=f*u*(2*s*s),
                    p=f*f*s,u_squared=u*u*s,weighted_phi_squared=f*f*s*s)
                box_moments={name:moments[name]+rhs[name]*partial for name in moments}
                cell_state=dict(phi=f,U=u,**box_moments)
                transfer=self.transfer(cell_state,s)
                saved.append(dict(y_interval=y,scaled_radius=s,source_phi=source,
                    alpha=alpha,log_slope=g,axial_slope=j,state=cell_state,transfer=transfer))
                moments={name:value+rhs[name]*step for name,value in moments.items()}
                G,J=G+g*step,J+j*step
            end=dict(phi=start['phi']*G.exp(),U=start['U']+J,**moments)
            if endpoints(end['phi'][0])[0]<=0:raise ValueError('comparison endpoint not positive')
            s=4*c.exp(2*hb);transfer=self.transfer(end,s)
            return dict(identity=self.identity,implicit_source_sha256=self.identity['implicit_source_sha256'],
                datum_enclosure_sha256=self.identity['datum_enclosure_sha256'],
                input_hashes=self.input_hashes,cells=cells,y_interval=['.005','.01'],
                core_branch_y_interval=['0','.005'],core_branch_definition='same admitted analytic core',
                endpoint_scaled_radius=s,endpoint_state=end,endpoint_transfer=transfer,
                unmodulated_driver_chi='1',initial_phi=self.initial_phi,
                minimum_source_phi_lower=c.mpf(minimum),Phi_log_increment=G,U_increment=J,
                cell_enclosures=saved,retained_axial_field_order=2,
                retained_axial_driver_order=transfer['D'].order,
                cell_integral_discretization_enclosed=True,core_analytic_tail_propagated=True,
                method='directed cell ranges and exact integrating-factor identity',
                amplitude_and_Lambda_factors_cancelled_before_transfer=True,
                physical_F0_not_materialized=True,new_source_only=True,
                comparison_hb_candidate='.005',paper_hb_cstar_K_inverse100_relation_verified=False,
                full_Section9_parameter_admission=False,
                terminal_five_moment_repair_completed=False,physical_exit_ODE_solved=False,
                admissible_stress_lift_constructed=False,whole_axis=False,temporal_recursion=False)

def run():
    calc=LogarithmicComparison();result=calc.integrate()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Md40 controlled comparison: 32 closed cells; fields C2, drivers C'+str(result['retained_axial_driver_order']))
    return result

if __name__=='__main__':run()
