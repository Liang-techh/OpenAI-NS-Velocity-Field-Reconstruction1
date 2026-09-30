"""Resolved reference join and similarity-derivative fixture."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_dual import AxialDual
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_field import evaluate_correction,join


def run():
    with mp.workdps(90):
        moment_map=FiveBumpMomentMap(precision=90,order=96)
        dual=lambda a,b=0:AxialDual(a,b,pressure_order=0,width_order=0)
        z=dual(mp.mpf('.3'),1);am=dual(5,mp.mpf('-.7'));rm=mp.mpf(16);dt=mp.mpf('.001')
        h=[dual(mp.mpf(a),mp.mpf(b)) for a,b in zip(['.00013','-.00007','.00011','-.00019','.00005'],['.00002','.00003','-.00004','.00001','.00006'])]
        val=lambda a:a.value.component(0,0) if isinstance(a,AxialDual) else a.component(0,0)
        keys=['theta','z','theta_z','z_theta','p'];records=[]
        for x in [mp.mpf('1.25'),mp.mpf('2.2')]:
            R=rm*x;root=mp.sqrt(2*R)
            theta=mp.sqrt(2)*rm**mp.mpf('1.5')*am*(x**mp.mpf('1.6')-1)/mp.mpf('1.6')
            moments=dict(theta=theta,z=4*z*rm*(x-1),theta_z=4*z*theta,
              z_theta=16*z*z*rm*(x-1)-rm*am*am*(x**mp.mpf('1.2')-1)/mp.mpf('2.4'),
              p=am*am*(x**mp.mpf('.2')-1)/mp.mpf('.4'))
            p0=dual(mp.mpf('-.2'),mp.mpf('.04'));p=p0+moments['p']
            raw=dict(axial=16*z*z*rm*(x-1),swirl=rm*am*am*(x**mp.mpf('1.2')-1)/mp.mpf('2.4'))
            base=dict(R=dual(R).value,Z='.3',logR=mp.log(R),delta=dt,
              P=p.value,P_Z=p.tangent,P0=p0.value,P0_Z=p0.tangent,
              moments={k:v.value for k,v in moments.items()},moments_Z={k:v.tangent for k,v in moments.items()},
              raw_quadratic_integrals={k:v.value for k,v in raw.items()},raw_quadratic_integrals_Z={k:v.tangent for k,v in raw.items()},
              moment_parts={'reference':{k:v.value for k,v in moments.items()}},moment_parts_Z={'reference':{k:v.tangent for k,v in moments.items()}},
              Utheta=(am*x**mp.mpf('.1')).value,Uz=(4*z).value,F=(am*x**mp.mpf('.1')/root).value)
            correction=evaluate_correction(moment_map,h,am,rm,z,x,delta=dt)
            result=join(base,correction,delta=dt)
            assert result['stress'] is not None
            assert result['P0']==p0.value and result['P0_Z']==p0.tangent
            assert abs(val(result['P'])-val(base['P']+correction['increments']['p']))<mp.mpf('1e-75')
            assert all(abs(val(result['moments'][k])-val(base['moments'][k]+correction['increments'][k]))<mp.mpf('1e-75') for k in keys)
            # Independent differentiation of the explicit angular field in log R.
            scalar_h=[val(v) for v in h]
            def u_at_y(y):
                coord=mp.exp(y)
                f=sum(scalar_h[j+2]*moment_map.beta(coord,moment_map.centers[j]) for j in range(3))
                return val(am)*(coord**mp.mpf('.1')+f)
            expected_y=mp.diff(u_at_y,mp.log(x))
            y_error=abs((val(result['Utheta_y'])-expected_y)/expected_y)
            assert y_error<mp.mpf('1e-65')
            expected_a=1-2*expected_y/val(result['Utheta'])
            assert abs(val(result['a'])-expected_a)<mp.mpf('1e-65')
            base_ur=(2*val(z)*R*val(4*z)-(1-dt)*val(z)*val(moments['z'])-(1-val(z)**2)*moments['z'].tangent.component(0,0))/((1-dt*val(z)**2)*root)
            ur_error=abs(val(result['Ur'])-base_ur-val(correction['delta_Ur']))
            assert ur_error<mp.mpf('1e-65')
            assert result['delta_Ur_Z'] is None
            corrected_raw=result['raw_quadratic_integrals']
            raw_identity=val(corrected_raw['axial'])-val(corrected_raw['swirl'])-val(result['moments']['z_theta'])
            assert abs(raw_identity)<mp.mpf('1e-75')
            assert result['raw_quadratic_correction_applied']
            assert all(abs(val(sum(rows[k] for rows in result['moment_parts'].values()))-val(result['moments'][k]))<mp.mpf('1e-75') for k in keys)
            assert all(abs(val(v))<mp.mpf('1e-75') for v in result['baseline_compatibility'].values())
            records.append(dict(x=str(x),logR_derivative_relative_error=mp.nstr(y_error,30),joined_radial_identity_absolute_error=mp.nstr(ur_error,30),P0_preserved=True))
        try:evaluate_correction(moment_map,h,am,rm,z,1,Rm_Z=1)
        except ValueError:fixed_radius_guard=True
        else:raise AssertionError('nonconstant Rm accepted')
        report=dict(records=records,fixed_radius_guard=fixed_radius_guard,raw_swirl_convention='integral(u_theta^2/2) dR; Mztheta=raw_axial-raw_swirl',actual_source_join=False,functional_closure=False,cone_certified=False,quadrature_enclosed=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2),flush=True)
        return report

if __name__=='__main__':run()
