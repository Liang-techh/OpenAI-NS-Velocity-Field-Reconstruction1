"""Independent physical-coordinate vector derivatives and original scale checks."""
import hashlib
import json
from pathlib import Path

import mpmath as mp
import sympy as s
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_compliant_cartesian_field import (
    CompliantCartesianField,cartesian_templates,cartesian_brackets,transverse_step,
    physical_time_bracket,INDICES,CS,SN,COMPONENTS,UR,UT,UZ,P)
from lei_ren_part1_paper_compliant_pulse_physical_bounds import PulsePhysicalBounds
from lei_ren_part1_paper_uniform_fixed_beta_error import read_interval
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints

HERE=Path(__file__).parent
NAME='lei_ren_part1_paper_compliant_cartesian_field.json'


def template_proof():
    """Derive basis derivatives in Cartesian coordinates, then check commutation."""
    x,y=s.symbols('x y',positive=True); r=s.sqrt(x*x+y*y)
    direct={('x',CS):s.diff(x/r,x),('x',SN):s.diff(y/r,x),
            ('y',CS):s.diff(x/r,y),('y',SN):s.diff(y/r,y)}
    expected={('x',CS):y*y/r**3,('x',SN):-x*y/r**3,
              ('y',CS):-x*y/r**3,('y',SN):x*x/r**3}
    for key,value in direct.items():
        if s.simplify(value-expected[key])!=0:raise ArithmeticError('Independent Cartesian basis derivative failed')
    template=cartesian_templates(); count=4
    for (component,i,j,b),row in template.items():
        if any(a+q!=i+j or a+b>4 for label,a,q in row):raise ArithmeticError('Cartesian degree incorrect')
        if i+j+b<=2:
            xy=transverse_step(transverse_step(row,'x'),'y')
            yx=transverse_step(transverse_step(row,'y'),'x')
            for key in set(xy)|set(yx):
                if s.expand(xy.get(key,0)-yx.get(key,0))!=0:raise ArithmeticError('Cartesian transverse derivatives do not commute')
            count+=1
    # Pure swirl still has nonzero Cartesian basis terms with zero radial velocity.
    dxuy=template[('uy',1,0,0)]
    if s.expand(dxuy[(UT,0,1)]-SN**2)!=0:raise ArithmeticError('Moving swirl basis term omitted')
    count+=1
    divergence={}
    for row in (template[('ux',1,0,0)],template[('uy',0,1,0)],template[('uz',0,0,1)]):
        for key,value in row.items():divergence[key]=divergence.get(key,0)+value
    # uz's z derivative is tracked by b, not the transverse row's radial a.
    del divergence[(UZ,0,0)]
    circle=s.groebner([CS**2+SN**2-1],CS,SN)
    expected_div={(UR,1,0):s.Integer(1),(UR,0,1):s.Integer(1)}
    for key in set(divergence)|set(expected_div):
        difference=s.Poly(s.expand(divergence.get(key,0)-expected_div.get(key,0)),CS,SN)
        if circle.reduce(difference.as_expr())[1]!=0:raise ArithmeticError('Cartesian/cylindrical divergence identity failed')
    count+=1
    return dict(independent_basis_and_commutator_checks=count,passed=True)


def implicit_cartesian_fixture():
    """Differentiate a Cartesian field after an independent positive lambda root."""
    with mp.workdps(85):
        c=MPIntervalContext(); c.dps=105; tol=mp.mpf('1e-60')
        delta=mp.mpf('.03'); tau=mp.mpf('.7'); z0=mp.mpf('.4'); y0=mp.mpf('.2'); theta=mp.mpf('.7')
        lam=mp.sqrt(tau/(1-z0*z0)); R=mp.exp(y0); rr=lam*mp.sqrt(2*R)
        xx=rr*mp.cos(theta); yy=rr*mp.sin(theta); zz=lam**(1-delta)*z0
        profiles={UR:lambda y,z:mp.exp(mp.mpf('.11')*y)*(1+mp.mpf('.3')*z+mp.mpf('.4')*z*z)+mp.exp(-mp.mpf('.09')*y)*mp.sin(z),
            UT:lambda y,z:mp.exp(-mp.mpf('.21')*y)*(1-mp.mpf('.2')*z+mp.mpf('.15')*z**3)+mp.exp(mp.mpf('.13')*y)*mp.cos(z),
            UZ:lambda y,z:mp.exp(mp.mpf('.07')*y)*(z+mp.mpf('.1')*mp.sin(z)),
            P:lambda y,z:mp.exp(-mp.mpf('.23')*y)*(1+z**4)+mp.exp(mp.mpf('.04')*y)*mp.sin(z)}
        beta={UR:mp.mpf(-1),UT:-1-delta,UZ:-1-delta,P:-2-2*delta}; grids={}
        for label,f in profiles.items():
            grids[label]={'y'+str(k)+'_Z'+str(n):c.mpf([v-tol,v+tol])
                for k in range(5) for n in range(5-k) for v in (mp.diff(f,(y0,z0),(k,n)),)}
        roots={}
        def physical(component,x,y,z,remaining_time=tau):
            # mp.diff raises its working precision for tiny finite differences.
            # Reusing a lower-precision root destroys high derivative fixtures.
            key=(z,remaining_time,mp.mp.prec)
            if key not in roots:
                roots[key]=mp.findroot(lambda l:l*l-l**(2*delta)*z*z-remaining_time,lam,tol=mp.eps*16,verify=True)
            ll=roots[key]; rad=mp.sqrt(x*x+y*y); Z=z/ll**(1-delta); Y=mp.log(rad*rad/(2*ll*ll))
            ur=ll**beta[UR]*profiles[UR](Y,Z); ut=ll**beta[UT]*profiles[UT](Y,Z)
            return {'ux':lambda:ur*x/rad-ut*y/rad,'uy':lambda:ur*y/rad+ut*x/rad,
                    'uz':lambda:ll**beta[UZ]*profiles[UZ](Y,Z),
                    'p':lambda:ll**beta[P]*profiles[P](Y,Z)}[component]()
        count=0; cache={}
        for i,j,b in INDICES:
            brackets=cartesian_brackets(c,grids,i,j,b,c.mpf(z0),c.mpf(delta),c.mpf(mp.cos(theta)),c.mpf(mp.sin(theta)),cache)
            for component,parts in brackets.items():
                mapped=c.mpf(0)
                for label,value in parts.items():
                    exponent=c.mpf(beta[label])-i-j+b*(c.mpf(delta)-1)
                    mapped+=value*c.exp(exponent*c.ln(c.mpf(lam))-(i+j)*c.ln(c.mpf(R))/2)
                target=mp.diff(lambda x,y,z:physical(component,x,y,z),(xx,yy,zz),(i,j,b))
                lo,hi=endpoints(mapped)
                if not lo<=target<=hi:raise ArithmeticError('Independent Cartesian vector derivative failed: '+str((component,i,j,b)))
                count+=1
        times={label:physical_time_bracket(c,grids[label],c.mpf(z0),c.mpf(delta),c.mpf(beta[label]))
               *c.exp((c.mpf(beta[label])-2)*c.ln(c.mpf(lam))) for label in profiles}
        mapped_times={'ux':times[UR]*c.mpf(mp.cos(theta))-times[UT]*c.mpf(mp.sin(theta)),
                      'uy':times[UR]*c.mpf(mp.sin(theta))+times[UT]*c.mpf(mp.cos(theta)),
                      'uz':times[UZ],'p':times[P]}
        for component,mapped in mapped_times.items():
            target=-mp.diff(lambda remaining:physical(component,xx,yy,zz,remaining),tau)
            if not endpoints(mapped)[0]<=target<=endpoints(mapped)[1]:raise ArithmeticError('Fixed-x physical time derivative failed')
        return dict(independent_cartesian_vector_derivatives_checked=count,
                    independent_fixed_x_time_derivatives_checked=4,
                    finite_parameter_fixture_only=True,actual_source_admission=False,passed=True)


def scale_fixture():
    with mp.workdps(80):
        c=MPIntervalContext(); c.dps=100; f=CompliantCartesianField.__new__(CompliantCartesianField)
        p=PulsePhysicalBounds.__new__(PulsePhysicalBounds); p.ctx=f.ctx=c
        p.mu=f.mu=c.mpf('.04'); p.delta=f.delta=c.mpf('.03'); p.logP=c.mpf('2.3')
        p.logRp_parts=dict(logCstar=c.mpf('3.7'),logPstar=c.mpf('23'),finite_outer_offset=c.mpf('4.1'))
        f.pulse=p; f.Evparts=dict(inlet_log=c.mpf('.34'),inverse_mu_term=-13/p.mu,finite_offset=c.mpf(-26))
        logRp=sum(p.logRp_parts.values(),c.mpf(0)); count=0
        for N,b in ((a,b) for a in range(5) for b in range(5-a)):
            for label in (UR,UT,UZ,P):
                for normalization,chi,offset in (('pulse','.02','0'),('pulse','13','-4'),('pulse','0','-1'),('Ev0','13','5.5')):
                    row=f.scale(label,N,b,normalization,c.mpf(chi),c.mpf(offset))
                    total=sum(row['log_prefactor_parts'].values(),c.mpf(0))
                    t=c.mpf(chi)/p.mu+c.mpf(offset)
                    if normalization=='pulse':
                        if label==UR:expected=p.logP+(1-N)*logRp/2-c.ln(2)/2-(p.mu+c.mpf(N)/2)*t
                        elif label==P:expected=2*p.logP-N*(logRp+t)/2
                        else:expected=p.logP-N*logRp/2-(c.mpf('.5')+p.mu+c.mpf(N)/2)*t
                    else:expected=(2*p.logP if label==P else p.logP+sum(f.Evparts.values(),c.mpf(0))/2)-N*(logRp+t)/2
                    if max(endpoints(total)[0],endpoints(expected)[0])>min(endpoints(total)[1],endpoints(expected)[1]):
                        raise ArithmeticError('Independent Cartesian scale normalization failed')
                    count+=1
        # First time derivatives add exactly lambda^-2, never a stage-time factor.
        for label in (UR,UT,UZ,P):
            row=f.scale(label,0,0,'Ev0',c.mpf(13),c.mpf('5.5'),time=True)
            base=f.scale(label,0,0,'Ev0',c.mpf(13),c.mpf('5.5'))
            if endpoints(row['physical_lambda_exponent']-base['physical_lambda_exponent'])[0]!=-2:
                # Directed interval subtraction need only enclose the exact equality.
                if not endpoints(row['physical_lambda_exponent']-base['physical_lambda_exponent'])[0]<=-2<=endpoints(row['physical_lambda_exponent']-base['physical_lambda_exponent'])[1]:
                    raise ArithmeticError('Physical time exponent incorrect')
        return dict(independent_original_normalization_scale_checks=count,passed=True)


def run():
    r=json.loads((HERE/NAME).read_bytes()); hashes=dict(r['input_hashes'])
    for name,digest in hashes.items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('Cartesian map source changed: '+name)
    c=MPIntervalContext(); c.dps=260; spatial=time=0
    if r['chart_count']!=16 or r['spatial_multiindex_count']!=35:raise ValueError('Cartesian chart/multiindex coverage missing')
    for chart,row in r['charts'].items():
        if len(row['spatial_multiindices'])!=35:raise ValueError('Missing Cartesian multiindices')
        for values in row['spatial_multiindices'].values():
            if set(values)!=set(COMPONENTS):raise ValueError('Missing Cartesian component')
            for contributions in values.values():
                for value in contributions.values():
                    norm=read_interval(c,value['absolute_upper'])
                    if any(not mp.isfinite(v) or v<0 for v in endpoints(norm)):raise ArithmeticError('Nonfinite Cartesian bracket norm')
                    if value['exactly_zero']!=(endpoints(norm)[1]==0):raise ArithmeticError('False Cartesian zero')
                    scale=row['shared_prefactor_bounds'][value['scale_key']]
                    if endpoints(read_interval(c,scale['physical_lambda_exponent']))[1]>=0:raise ArithmeticError('Lambda bound direction changed')
                    for part in scale['log_prefactor_parts'].values():
                        if any(not mp.isfinite(v) for v in endpoints(read_interval(c,part))):raise ArithmeticError('Lost formal finite Cartesian scale')
                    spatial+=1
        if chart not in ('local_O3_power','pulse_entrance','pulse_main','pulse_exit','pulse_gap_main','pulse_gap_end','pulse_end'):
            for values in row['spatial_multiindices'].values():
                for contributions in values.values():
                    for label,value in contributions.items():
                        if label in (UR,UZ) and not value['exactly_zero']:raise ArithmeticError('Inherited postpulse zero history lost')
        for value in row['first_physical_time_derivative_cylindrical'].values():
            if any(not mp.isfinite(v) for v in endpoints(read_interval(c,value['absolute_upper']))):raise ArithmeticError('Physical time bracket not finite')
            time+=1
    for key in ('full_cartesian_vector_derivatives_certified','core_axis_interfaces_certified','whole_background_installed',
                'physical_energy_integral_certified','whole_outer_cone_certified','temporal_recursion'):
        if r[key]:raise ValueError('Cartesian chart map scope overclaimed: '+key)
    proof=template_proof(); fixture=implicit_cartesian_fixture(); scales=scale_fixture()
    hashes[NAME]=hashlib.sha256((HERE/NAME).read_bytes()).hexdigest()
    hashes[Path(__file__).name]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    out=dict(actual_five_defect_family_sha256=r['actual_five_defect_family_sha256'],implicit_source_sha256=r['implicit_source_sha256'],
        template_checks=proof,independent_cartesian_fixture=fixture,original_scale_fixture=scales,
        actual_cartesian_contribution_bounds_checked=spatial,actual_fixed_x_time_brackets_checked=time,
        accepted_outer_chart_cartesian_spatial_C4_mapped=True,accepted_outer_chart_first_physical_time_derivative_mapped=True,
        full_cartesian_vector_derivatives_certified=False,core_axis_interfaces_certified=False,
        whole_background_installed=False,physical_energy_integral_certified=False,whole_outer_cone_certified=False,temporal_recursion=False,
        all_passed=True,input_hashes=hashes)
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print('Cartesian vector map: independent spatial4/fixed-x time derivatives, moving basis and original radius/amplitude scales PASS',flush=True)
    return out


if __name__=='__main__':run()
