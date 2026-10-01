"""Fresh O.1 continuation and first O.2 slope transition (Eq. 4.6).

The terminal implicit moment identities, not midpoint repair controls, supply
the reference branch. Directed rectangles enclose J and all new radial masses.
Only y=log(R/Rref)<=1 is implemented here; axial turnoff follows this stage.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_interval_repaired_reference_field import validate_receipt
from lei_ren_part1_paper_interval_long_reshape_field import sigma_value_derivative, physical_packet, RM_LOG_EXACT
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_reference_endpoint_targets import accepted_parameters, ACCEPTED_SHA

HERE=Path(__file__).parent


def fraction_box(c,q):
    q=Fraction(q)
    return c.mpf(q.numerator)/q.denominator


def transition_integrals(c,y,cells=256):
    """Enclose J(y) and three dimensionless masses, without quadrature fits.

    sigma is monotone. Endpoint rectangles bound its primitive on each cell.
    On a cell s in [a,b], J(s) is between J_lower(a), J_upper(b).
    Thus exp(k*s-m*.6*J(s)) has a directed enclosing interval.
    """
    y=Fraction(y)
    if not 0<=y<=1 or not isinstance(cells,int) or cells<1:
        raise ValueError('y in [0,1] and positive integer cells required')
    j=c.mpf(0); masses=[c.mpf(0) for _ in range(3)]
    rates=(c.mpf('1.6'),c.mpf('.2'),c.mpf('1.2'))
    powers=(1,2,2)
    dy=fraction_box(c,y/cells)
    for i in range(cells):
        a=fraction_box(c,y*i/cells);b=fraction_box(c,y*(i+1)/cells)
        sa=sigma_value_derivative(c,a)[0];sb=sigma_value_derivative(c,b)[0]
        sj=c.mpf([endpoints(sa)[0],endpoints(sb)[1]])
        next_j=j+dy*sj
        jcell=c.mpf([endpoints(j)[0],endpoints(next_j)[1]])
        scell=c.mpf([endpoints(a)[0],endpoints(b)[1]])
        for k,(rate,power) in enumerate(zip(rates,powers)):
            masses[k]+=dy*c.exp(rate*scell-c.mpf('.6')*power*jcell)
        j=next_j
    # Symmetry sigma(1-s)=1-sigma(s) gives the exact endpoint J(1)=1/2.
    if y==1:j=c.mpf('.5')
    return j,masses


def reference_moments(c,z,R,u):
    theta=u*(c.mpf(5)/8*c.sqrt(2)*R**c.mpf('1.5'))
    return dict(z=z*(4*R),theta=theta,theta_z=z*theta*4,
        z_theta=z*z*(16*R)-u*u*(c.mpf(5)/12*R),p=u*u*c.mpf('2.5'))


def evaluate_transition(c,z,delta,Rref,A,P0,y,cells=256):
    """C1 axial-family field for first slope transition, with inherited masses."""
    q=Fraction(y);yy=fraction_box(c,q)
    j,mass=transition_integrals(c,q,cells)
    sig=sigma_value_derivative(c,yy)[0]
    factor=c.exp(yy/10-c.mpf('.6')*j)
    u=A*factor;slope=c.mpf('.1')-c.mpf('.6')*sig
    uy=u*slope;R=Rref*c.exp(yy)
    inlet=reference_moments(c,z,Rref,A)
    dtheta=A*(c.sqrt(2)*Rref**c.mpf('1.5')*mass[0])
    dp=A*A*(mass[1]/2)
    denergy=A*A*(Rref*mass[2]/2)
    # Entire axial contribution remains analytic: Uz=4Z on this stage.
    moments=dict(z=z*(4*R),theta=inlet['theta']+dtheta,
        theta_z=z*(inlet['theta']+dtheta)*4,
        z_theta=z*z*(16*R)-A*A*(c.mpf(5)/12*Rref)-denergy,
        p=inlet['p']+dp)
    normalized_shear=c.mpf(-2) if q==1 else (c.mpf('-.8') if q==0 else 2*slope-1)
    packet=physical_packet(c,z,delta,R,u,uy,z*4,moments,P0,normalized_shear)
    packet.update(y=str(q),stage='O.2 angular slope transition',J=j,
        dimensionless_increment_integrals=mass,inherited_reference_moments=inlet,
        slope=slope,directed_integral_cells=cells,
        all_five_moments_transported=True,P0_preserved=True,
        axial_turnoff_installed=False,whole_outer_cone_certified=False,
        exact_heat_exterior_installed=False,temporal_recursion=False)
    return packet


class IntervalOuterSlopeField:
    def __init__(self):
        self.accepted_parameters=accepted_parameters()
        self.calc=IntervalComparisonJets(4);c=self.calc.ctx
        self.inverse_name='lei_ren_part1_paper_interval_five_bump_inverse.json'
        self.defects_name='lei_ren_part1_paper_interval_functional_defects.json'
        for name in (self.inverse_name,self.defects_name):
            raw=json.loads((HERE/name).read_bytes());validate_receipt(raw,name)
            if raw['state_sha256']!=self.calc.state_hash:raise ValueError('core mismatch')
            if name==self.inverse_name:
                if not raw.get('certified') or not raw.get('uniform_implicit_C1_family_exists'):
                    raise ValueError('uniform terminal implicit identities required')
        with mp.workdps(self.calc.precision+60):
            self.z=self.calc.z.truncate(1);self.P0=self.calc.p0.truncate(1)
            self.A=(1+self.z*self.z).reciprocal()*c.exp(14)
            self.Rref=110*c.exp(fraction_box(c,RM_LOG_EXACT+6))

    def evaluate_reference_x(self,x):
        """Exact repaired-terminal reference branch, x=R/Rm from 2 to exp(6)."""
        with mp.workdps(self.calc.precision+60):
            c=self.calc.ctx;q=Fraction(x);xx=fraction_box(c,q)
            if q<2 or endpoints(c.ln(xx))[1]>6:raise ValueError('2<=x<=exp(6) required')
            R=self.Rref*xx*c.exp(-6);u=self.A*c.exp((c.ln(xx)-6)/10)
            m=reference_moments(c,self.z,R,u)
            packet=physical_packet(c,self.z,self.calc.delta,R,u,u/10,self.z*4,m,self.P0,c.mpf('-.8'))
            packet.update(x=str(q),stage='O.1 exact reference continuation',
                terminal_identity_source='uniform implicit five-moment inverse, not interval midpoint',
                all_five_moments_transported=True,P0_preserved=True,temporal_recursion=False)
            return packet

    def evaluate_y(self,y,cells=256):
        with mp.workdps(self.calc.precision+60):
            return evaluate_transition(self.calc.ctx,self.z,self.calc.delta,self.Rref,self.A,self.P0,y,cells)


def run():
    field=IntervalOuterSlopeField()
    samples=[field.evaluate_reference_x(x) for x in (2,10,100)]
    samples += [field.evaluate_y(y) for y in ('0','.5','1')]
    result=dict(state_sha256=field.calc.state_hash,center_family=field.calc.center_family,
        source_equation='Lei-Ren Part I Eq. (4.6); O.1 exact reference continuation',
        source_text_lines='lei_ren_part1.txt 3925-4055, 5941-6177, 17720-18026',
        fresh_C1_outer_slope_field_installed=True,all_five_cumulative_moments_retained=True,
        analytic_axis_pressure_retained=True,samples=samples,
        accepted_pressure_schedule_sha256=ACCEPTED_SHA,
        accepted_Md=field.accepted_parameters['Md'],
        accepted_Md_satisfies_paper_greater_than_one=Fraction(field.accepted_parameters['Md'])>1,
        Md_not_used_in_this_initial_slope_segment=True,
        original_parameter_errors_enclosed=False,whole_axis_field=False,
        whole_outer_admissibility_certified=False,axial_turnoff_installed=False,
        heat_exterior_matched=False,temporal_recursion=False)
    names=(Path(__file__).name,field.inverse_name,field.defects_name,
        'lei_ren_part1_paper_interval_long_reshape_field.py',
        'lei_ren_part1_paper_interval_repaired_reference_field.py',
        'lei_ren_part1_paper_reference_endpoint_targets.py',
        'lei_ren_part1_paper_coherent_pressure_source_alignment.json')
    result['input_hashes']={n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Fresh reference and angular slope field installed;',[(r.get('x',r.get('y')),r['cone']['status']) for r in samples],flush=True)
    return result

if __name__=='__main__':run()
