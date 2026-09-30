"""Independent field-density replay of the finite Section10 five-bump map."""
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_five_bump_map import FiveBumpMomentMap
from lei_ren_part1_paper_axial_dual import AxialDual


def run():
    with mp.workdps(90):
        moment_map = FiveBumpMomentMap(precision=90, order=96)
        ell=mp.mpf(1)/40
        centers=[mp.mpf(5)/4,mp.mpf(3)/2,mp.mpf(7)/4]
        normalization=ell*mp.quad(lambda s:mp.exp(-1/(1-s*s)),[-1,0,1])
        def bump(x,center):
            s=(x-center)/ell
            return mp.exp(-1/(1-s*s))/normalization if abs(s)<1 else mp.mpf(0)
        h=[mp.mpf(s) for s in ['.00013','-.00007','.00011','-.00019','.00005']]
        hz=[mp.mpf(s) for s in ['.00002','.00003','-.00004','.00001','.00006']]
        z=mp.mpf('.3');am=mp.exp(mp.mpf('1.4'))/(1+z*z)
        amz=-2*z*am/(1+z*z)
        cuts=[1]+[v for c in centers for v in (c-ell,c,c+ell)]+[2]
        def direct(coeff,amplitude,left=mp.mpf(1),right=mp.mpf(2)):
            def fields(x):
                b=[bump(x,c) for c in centers]
                f=sum(coeff[j+2]*b[j] for j in range(3))
                g=coeff[0]*b[0]+coeff[1]*b[2]
                u=x**mp.mpf('.1')+f
                return f,g,u
            def density(x,row):
                f,g,u=fields(x)
                if row==0:return g
                if row==1:return mp.sqrt(x)*u*g
                if row==2:return mp.sqrt(x)*f
                if row==3:return (g/amplitude)**2-(u*u-x**mp.mpf('.2'))/2
                return (u*u-x**mp.mpf('.2'))/(2*x)
            interval_cuts=sorted({left,right,*[c for c in cuts if left<c<right]})
            return [mp.quad(lambda x:density(x,row),interval_cuts) for row in range(5)]
        expected=direct(h,am)
        actual=moment_map.apply(h,am)
        def scalar(v):
            if hasattr(v,'value'):v=v.value
            if hasattr(v,'evaluate'):v=v.evaluate(pressure=1,width=1)
            return v
        relative=[abs((scalar(a)-b)/b) for a,b in zip(actual,expected)]
        hd=[AxialDual(a,b,pressure_order=0,width_order=0) for a,b in zip(h,hz)]
        ad=AxialDual(am,amz,pressure_order=0,width_order=0)
        out=moment_map.apply(hd,ad)
        step=mp.mpf('1e-15')
        # Independent scalar four-point stencil differentiates the field densities.
        samples={k:direct([a+k*step*b for a,b in zip(h,hz)],am+k*step*amz) for k in (-2,-1,1,2)}
        expected_z=[(-samples[2][i]+8*samples[1][i]-8*samples[-1][i]+samples[-2][i])/(12*step) for i in range(5)]
        tangent_relative=[abs((scalar(a.tangent)-b)/b) for a,b in zip(out,expected_z)]
        rhs=[mp.mpf(s) for s in ['.001','-.002','.003','-.004','.005']]
        coeff=moment_map.linear_inverse(rhs)
        # Linear density replay, independently remove its quadratic field contributions.
        replay=direct(coeff,am)
        q=moment_map.quadratic(coeff,am)
        linear_relative=[abs((scalar(v)-scalar(n)-r)/r) for v,n,r in zip(replay,q,rhs)]
        quadratic=moment_map.quadratic(h,am)
        bilinear=moment_map.bilinear(h,h,am)
        assert all(a==b for a,b in zip(quadratic,bilinear))
        cross=moment_map.bilinear(h,hz,am)
        reverse=moment_map.bilinear(hz,h,am)
        assert max(abs(a-b) for a,b in zip(cross,reverse))<mp.mpf('1e-80')
        jac=moment_map.jacobian(h,am)
        fixed_amp=moment_map.apply(hd,AxialDual(am,0,pressure_order=0,width_order=0))
        jac_error=max(abs(sum(jac[i][j]*hz[j] for j in range(5))-scalar(fixed_amp[i].tangent)) for i in range(5))
        assert jac_error<mp.mpf('1e-70')
        split=mp.mpf('1.4')
        partial=moment_map.partial(h,am,1,split)
        expected_partial=direct(h,am,mp.mpf(1),split)
        partial_relative=[abs((a-b)/b) for a,b in zip(partial,expected_partial)]
        assert max(partial_relative)<mp.mpf('1e-12')
        empty=moment_map.partial(h,am,1,1)
        assert all(v==0 for v in empty)
        assert max(relative)<mp.mpf('1e-12')
        assert max(tangent_relative)<mp.mpf('1e-12')
        assert max(linear_relative)<mp.mpf('1e-12')
        report=dict(row_relative_errors=[mp.nstr(x,30) for x in relative],
          first_Z_relative_errors=[mp.nstr(x,30) for x in tangent_relative],
          partial_relative_errors=[mp.nstr(x,30) for x in partial_relative],jacobian_direction_absolute_error=mp.nstr(jac_error,30),
          bilinear_self_and_symmetry_pass=True,empty_partial_zero=True,
          inverse_relative_errors=[mp.nstr(x,30) for x in linear_relative],
          independent_replay='normalized bump field densities integrated with mp.quad',
          quadrature_enclosed=False,actual_defect_solve=False,functional_closure=False,cone_certified=False)
        Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n')
        print(json.dumps(report,indent=2),flush=True)
        return report

if __name__=='__main__':run()
