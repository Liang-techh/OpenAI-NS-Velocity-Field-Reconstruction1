"""Coupled-parameter comparison in phase xi=y/h_b, not a fake width.

The same implicit positive h is used in radii, integration steps and the
physical exit plateau. Numeric boxes uniformly enclose its allowed family.
"""
import hashlib
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_logarithmic_comparison import LogarithmicComparison
from lei_ren_part1_paper_shared_exit_parameters import SharedExitParameters
from lei_ren_part1_paper_interval_comparison_enclosure import alpha_box
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
HERE=Path(__file__).parent
class SharedComparison(LogarithmicComparison):
    def __init__(self):
        super().__init__();self.parameters=SharedExitParameters(self.ctx)
        self.h=self.parameters.h_box
        self.input_hashes.update(self.parameters.input_hashes)
        self.input_hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    def integrate(self,cells=32):
        if not isinstance(cells,int) or cells<1:raise ValueError('positive cell count')
        c=self.ctx
        with mp.workdps(c.dps+40):
            step=self.h/cells;start=self.core_state(4*c.exp(self.h));zero=start['phi']*0
            G,J=zero,zero;moments={n:v for n,v in start.items() if n not in ('phi','U')}
            saved=[];minimum=None
            for n in range(cells):
                left=1+c.mpf(n)/cells;right=1+c.mpf(n+1)/cells
                xi=c.mpf([endpoints(left)[0],endpoints(right)[1]])
                y=self.h*xi;s=4*c.exp(y);source=self.field(s,'Phi')
                lower=endpoints(source[0])[0]
                if lower<=0:raise ValueError('shared comparison phi denominator')
                minimum=lower if minimum is None else min(minimum,lower)
                alpha=alpha_box(c,xi,c.mpf(1))
                g=self.field(s,'Phi',1)*s*alpha/source
                j=self.field(s,'Uz',1)*s*alpha
                partial=c.mpf([0,endpoints(step)[1]])
                f=start['phi']*(G+g*partial).exp();u=start['U']+J+j*partial
                rhs=dict(theta=f*(2*s*s),z=u*s,theta_z=f*u*(2*s*s),
                    p=f*f*s,u_squared=u*u*s,weighted_phi_squared=f*f*s*s)
                state=dict(phi=f,U=u,**{n:moments[n]+rhs[n]*partial for n in moments})
                saved.append(dict(phase_xi=xi,y_interval=y,step=step,scaled_radius=s,
                    alpha=alpha,state=state,transfer=self.transfer(state,s),source_phi=source))
                moments={n:v+rhs[n]*step for n,v in moments.items()}
                G,J=G+g*step,J+j*step
            end=dict(phi=start['phi']*G.exp(),U=start['U']+J,**moments)
            s=4*c.exp(2*self.h)
            return dict(identity=self.identity,input_hashes=self.input_hashes,
                parameter_family_sha256=self.parameters.sha,parameters=self.parameters.report(),
                phase_domain=['1','2'],core_phase_domain=['0','1'],cells=cells,
                endpoint_scaled_radius=s,endpoint_state=end,endpoint_transfer=self.transfer(end,s),
                initial_phi=self.initial_phi,cell_enclosures=saved,minimum_source_phi_lower=c.mpf(minimum),
                Phi_log_increment=G,U_increment=J,
                h_b_equals_exit_epsilon_by_definition=True,normalized_phase_cutoff_exact=True,
                positive_h_family_enclosed_without_dividing_by_zero=True,
                analytic_core_tail_propagated=True,comparison_discretization_enclosed=True,
                retained_axial_field_order=2,retained_axial_driver_order=1,
                full_Section9_parameter_admission=False,full_K_norm_certificate_missing=True,
                whole_axis=False,terminal_five_moment_repair_completed=False,temporal_recursion=False)
def run():
    calc=SharedComparison();r=calc.integrate()
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(r)),indent=2)+'\n',encoding='utf-8')
    print('Shared-parameter comparison: 32 phase cells; identical h_b and exit epsilon')
    return r
if __name__=='__main__':run()
