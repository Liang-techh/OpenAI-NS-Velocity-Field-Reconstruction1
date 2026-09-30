"""Independent scalar field/moment replay of the termwise correction adapter."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_dual import AxialDual
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_five_bump_response import build_response
from lei_ren_part1_paper_five_bump_response_field import evaluate_response_correction


def run():
    with mp.workdps(90):
        moment_map=FiveBumpMomentMap(precision=90,order=96)
        rm=mp.mpf(16);z=mp.mpf('.3');dt=mp.mpf('.001')
        am=mp.mpf(5);amz=mp.mpf('-.7');degree=2
        amplitude=AxialDual(am,amz,pressure_order=0,width_order=0)
        response=build_response(moment_map,amplitude,degree=degree)
        d=[mp.mpf(v)*mp.mpf('1e-9') for v in [1,-2,3,-4,5]]
        dz=[mp.mpf(v)*mp.mpf('1e-10') for v in [2,1,-1,3,-2]]
        defects=[AxialDual(a,b,pressure_order=0,width_order=0) for a,b in zip(d,dz)]
        def scalar(v):
            if isinstance(v,AxialDual):v=v.value
            return v.component(0,0) if hasattr(v,'component') else v
        coeff=response.evaluate_terms(defects)
        h=[sum(rows[i] for rows in coeff.values()) for i in range(5)]
        ell=mp.mpf(1)/40;centers=[mp.mpf(5)/4,mp.mpf(3)/2,mp.mpf(7)/4]
        norm=ell*mp.quad(lambda s:mp.exp(-1/(1-s*s)),[-1,0,1])
        def beta(x,c):
            t=(x-c)/ell
            return mp.exp(-1/(1-t*t))/norm if abs(t)<1 else mp.mpf(0)
        def direct(coeff,a,x):
            end=min(x,mp.mpf(2));cuts=sorted({mp.mpf(1),end,*[c for s in centers for c in (s-ell,s,s+ell) if 1<c<end]})
            def field(y):
                b=[beta(y,c) for c in centers]
                f=sum(coeff[j+2]*b[j] for j in range(3));g=coeff[0]*b[0]+coeff[1]*b[2]
                return a*y**mp.mpf('.1'),a*f,g
            def densities(y,i):
                u,du,g=field(y)
                if i==0:return mp.sqrt(2*rm*y)*du*rm
                if i==1:return g*rm
                if i==2:return mp.sqrt(2*rm*y)*(4*z*du+(u+du)*g)*rm
                if i==3:return (8*z*g+g*g-(2*u*du+du*du)/2)*rm
                return (2*u*du+du*du)/(2*y)
            return [mp.quad(lambda y:densities(y,i),cuts) for i in range(5)]
        records=[]
        for x in [mp.mpf('1.25'),mp.mpf('2.2')]:
            result=evaluate_response_correction(response,defects,Rm=rm,Z=z,x=x,delta=dt)
            actual=[sum(rows[i] for rows in result['moment_terms'].values()) for i in range(5)]
            expected=direct([scalar(v) for v in h],am,x)
            relative=[abs((scalar(a)-b)/b) for a,b in zip(actual,expected)]
            assert max(relative)<mp.mpf('1e-12')
            field=[sum(rows[i] for rows in result['velocity_terms'].values()) for i in range(3)]
            step=mp.mpf('1e-10')
            # Independent Z differentiation of scalar moment integrals.
            def shifted(k):
                rr=build_response(moment_map,am+k*step*amz,degree=degree)
                cc=rr.evaluate_terms([a+k*step*b for a,b in zip(d,dz)])
                hh=[sum(rows[i] for rows in cc.values()) for i in range(5)]
                return direct(hh,am+k*step*amz,x)[1]
            samples={k:shifted(k) for k in (-2,-1,1,2)}
            mz_z=(-samples[2]+8*samples[1]-8*samples[-1]+samples[-2])/(12*step)
            R=rm*x
            expected_ur=(2*z*R*scalar(field[1])-(1-dt)*z*expected[1]-(1-z*z)*mz_z)/((1-dt*z*z)*mp.sqrt(2*R))
            actual_ur=sum(scalar(v) for v in result['radial_velocity_terms'].values())
            ur_relative=abs((actual_ur-expected_ur)/expected_ur)
            assert ur_relative<mp.mpf('1e-12')
            if x>2:assert all(scalar(v)==0 for v in field)
            assert any(sum(key)>degree and any(scalar(v)!=0 for v in rows) for key,rows in result['moment_terms'].items())
            records.append(dict(x=str(x),moment_relative_errors=[mp.nstr(v,30) for v in relative],radial_velocity_relative_error=mp.nstr(ur_relative,30),moment_monomial_count=len(result['moment_terms'])))
        report=dict(records=records,velocity_degree=degree,quadratic_moment_degree=2*degree,
          higher_quadratic_terms_nonzero=True,independent_replay='scalar normalized velocity-density integration and independent Z stencil',
          quadrature_enclosed=False,actual_source_field_installed=False,functional_closure=False,radial_first_Z_available=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2),flush=True)
        return report

if __name__=='__main__':run()
