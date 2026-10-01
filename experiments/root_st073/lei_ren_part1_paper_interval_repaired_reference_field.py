"""Directed C1 five-bump reference-annulus field on the fresh axial family.

The controls enclose the unique implicit inverse family, not midpoint controls.
Cumulative moments and pressure use the same source defects and partial map.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_jet,_pack
from lei_ren_part1_paper_interval_exit_stress_enclosure import cone
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
HERE=Path(__file__).parent


def stress_values(c,z,delta,R,u,uy,uz,uzy,m,P):
    zv=z[0];d=1-zv*zv;L=1-delta*zv*zv;root=c.sqrt(2*R)
    transport=-R+(1-delta)*zv*m['z'][0]+d*m['z'][1]
    Itheta=u[0]*transport/(L*root)+((1-delta/2)*m['theta'][0]
        -(1-delta)*zv*m['theta'][1]/2-d*m['theta_z'][1]
        +(2*delta-1)*zv*m['theta_z'][0])/(2*L*R)
    Iz=(transport*uz[0]+(1-delta)*(m['z'][0]-zv*m['z'][1])/2
        +2*delta*zv*m['z_theta'][0]-d*m['z_theta'][1]
        +R*(2*(1+delta)*zv*P[0]-d*P[1]))/(L*root)
    F=u[0]/root
    st=(2*uy[0]-u[0])/u[0];sz=2*uzy[0]/u[0]
    tt=Itheta/F+st;tz=Iz/F+sz
    return dict(stress=dict(I_theta=Itheta,I_z=Iz,S_theta=F*st,S_z=F*sz,
        T_theta=F*tt,T_z=F*tz),normalized_stress=dict(S_theta_over_F=st,
        S_z_over_F=sz,T_theta_over_F=tt,T_z_over_F=tz),cone=cone(c,st,sz,tt,tz))


def evaluate_reference(c,z,delta,Rm,Am,P0,h,d,x,moment_map):
    """Point radial x, complete enclosed axial family; cumulative C1 moments."""
    q=Fraction(x);xx=c.mpf(q.numerator)/q.denominator
    if not Fraction(1)<=q<=Fraction(2):raise ValueError('reference repair x must lie in [1,2]')
    rows=moment_map.apply(h,(Am*Am).reciprocal(),q)
    bumps=moment_map.bump_values(q)
    beta,beta_x=bumps['values'],bumps['derivatives']
    zero=h[0]*0
    f=sum((h[i+2]*beta[i] for i in range(3)),zero)
    fx=sum((h[i+2]*beta_x[i] for i in range(3)),zero)
    g=h[0]*beta[0]+h[1]*beta[2];gx=h[0]*beta_x[0]+h[1]*beta_x[2]
    power=xx**c.mpf('.1');R=Rm*xx;root=c.sqrt(2*R)
    u=Am*(f+power);uy=Am*(fx*xx+power*c.mpf('.1'))
    uz=z*4+g;uzy=gx*xx
    F=u/root
    reference=dict(z=z*(4*R),theta=Am*(c.mpf(5)/8*c.sqrt(2)*R**c.mpf('1.5')*power))
    reference['theta_z']=z*reference['theta']*4
    reference['z_theta']=z*z*(16*R)-Am*Am*(c.mpf(5)/12*R*power*power)
    reference['p']=Am*Am*(c.mpf('2.5')*power*power)
    centered=[d[i]+rows[i] for i in range(5)]
    scale=Am*(c.sqrt(2)*Rm**c.mpf('1.5'))
    increments=dict(z=centered[0]*Rm,theta=centered[2]*scale,
        theta_z=(centered[1]+z*centered[2]*4)*scale,
        z_theta=centered[3]*Am*Am*Rm+z*centered[0]*(8*Rm),p=centered[4]*Am*Am)
    moments={key:reference[key]+increments[key] for key in reference}
    P=P0+moments['p']
    zv=z[0];L=1-delta*zv*zv;mz=moments['z']
    Ur=(2*zv*R*uz[0]-(1-delta)*zv*mz[0]-(1-zv*zv)*mz[1])/(L*root)
    Ur_R=((1+delta)*zv*uz[0]+2*zv*uzy[0]-(1-zv*zv)*uz[1])/(L*root)-Ur/(2*R)
    radial_moment_rhs=dict(z=uz,theta=u*root,theta_z=u*uz*root,
        z_theta=uz*uz-u*u/2,p=u*u/(2*R))
    stress=stress_values(c,z,delta,R,u,uy,uz,uzy,moments,P)
    return dict(x=str(q),R=R,logR=c.ln(R),F=F,Utheta=u,Uz=uz,P=P,P0=P0,
        Utheta_y=uy,Uz_y=uzy,F_y=(uy-u/2)/root,F_R=(uy-u/2)/(root*R),
        Ur_value=Ur,Ur_R_value=Ur_R,physical_moments=moments,reference_moments=reference,
        centered_terminal_rows=centered,partial_correction_rows=rows,
        correction_fields=dict(f=f,f_x=fx,g=g,g_x=gx),moment_radial_rhs=radial_moment_rhs,
        pressure_from_same_axis_datum_and_pressure_moment=True,
        moment_radial_identities_used_for_divergence_recovery=True,
        Ur_Z_available=False,retained_axial_order=1,**stress)


def validate_receipt(raw,label):
    for name,digest in raw.get('input_hashes',{}).items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:
            raise ValueError(label+' dependency changed:'+name)


class IntervalRepairedReferenceField:
    def __init__(self):
        from lei_ren_part1_paper_interval_partial_five_bump_map import IntervalPartialFiveBumpMap
        self.calc=IntervalComparisonJets(4);c=self.calc.ctx
        self.inverse_name='lei_ren_part1_paper_interval_five_bump_inverse.json'
        self.defects_name='lei_ren_part1_paper_interval_functional_defects.json'
        inverse=json.loads((HERE/self.inverse_name).read_text())
        defects=json.loads((HERE/self.defects_name).read_text())
        validate_receipt(inverse,'inverse');validate_receipt(defects,'defects')
        if not inverse.get('certified') or not inverse.get('uniform_implicit_C1_family_exists'):
            raise ValueError('uniform implicit inverse certificate required')
        if any(raw['state_sha256']!=self.calc.state_hash for raw in (inverse,defects)):
            raise ValueError('source core mismatch')
        self.inverse=inverse;self.defects=defects
        with mp.workdps(self.calc.precision+60):
            self.h=[restore_jet(c,v,order=1) for v in inverse['controls']]
            self.d=[restore_jet(c,defects['rows'][str(i)],order=1) for i in range(1,6)]
            self.z=self.calc.z.truncate(1);self.P0=self.calc.p0.truncate(1)
            self.Am=(1+self.z*self.z).reciprocal()*c.exp(c.mpf('13.4'))
            self.Rm=110*c.exp(10*(c.mpf('5e151')+14)-6)
            self.map=IntervalPartialFiveBumpMap(c)
    def evaluate_x(self,x):
        with mp.workdps(self.calc.precision+60):
            data=evaluate_reference(self.calc.ctx,self.z,self.calc.delta,self.Rm,
                self.Am,self.P0,self.h,self.d,x,self.map)
            data.update(center_family=self.calc.center_family,implicit_inverse_family_enclosed=True,
                controls_projected_to_midpoints=False,point_radial_coordinate=True,
                whole_axis_field=False,Cartesian_divergence_independently_validated=False,
                whole_reference_annulus_cone_certified=False,temporal_recursion=False)
            if Fraction(x)==2:
                data['terminal_moment_identities_implied_by_uniform_inverse']=True
                data['terminal_zero_containment_diagnostic']=all(endpoints(v[k])[0]<=0<=endpoints(v[k])[1]
                    for v in data['centered_terminal_rows'] for k in (0,1))
            return data


def run():
    field=IntervalRepairedReferenceField()
    points=('1','1.25','1.5','1.75','2')
    rows=[field.evaluate_x(x) for x in points]
    result=dict(center_family=field.calc.center_family,state_sha256=field.calc.state_hash,
        callable_reference_annulus_field_installed=True,controlled_partial_moment_integrals=True,
        P0_preserved=True,implicit_C1_controls_used=True,samples=rows,
        terminal_identity_scope='local axial family, reference-annulus five-bump map at x=2',
        terminal_moment_identities_implied_by_uniform_inverse=True,
        whole_axis_or_whole_profile_matching_complete=False,heat_exterior_matched=False,
        original_construction_parameter_errors_enclosed=False,temporal_recursion=False)
    names=(Path(__file__).name,field.inverse_name,field.defects_name,
        'lei_ren_part1_paper_interval_partial_five_bump_map.py',
        'lei_ren_part1_paper_bump_integral_enclosures_check.json',
        'lei_ren_part1_paper_mp_stress.py')
    result['input_hashes']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(_pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Repaired reference field installed;',[(r['x'],r['cone']['status']) for r in rows],flush=True)
    print('Terminal zero containment diagnostic:',rows[-1]['terminal_zero_containment_diagnostic'],flush=True)
    return result

if __name__=='__main__':run()
