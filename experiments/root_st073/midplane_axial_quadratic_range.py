"""Resolved signed axial-moment quadratic range in the current finite basis."""
import json
import numpy as np
from adaptive_bridge_moment_fit import moment_slices, outer_moments
from adaptive_bridge_recursive_defect import build_fields
from midplane_resolved_feasibility import (ZeroBackground, RADIAL_BREAKS,
    _integrated_moment_coefficients, _evaluate)
from radial_continuation import ROOT
from separated_moment_modes import SeparatedMomentModes, RADIAL_WINDOWS_THREE


def run():
    inner, fields = build_fields()
    base = fields['two_sided_cone']
    a = np.array(json.loads((ROOT/'midplane_axial_cone_all_knots_repair.json').read_text())['amplitudes'])
    def field(v, background=base):
        return SeparatedMomentModes(background, v, windows=RADIAL_WINDOWS_THREE, knots=(11.,15.,19.))
    current = field(a)
    rows=[]
    for block,k in enumerate((11.,15.,19.)):
        indices=np.arange(12*block,12*(block+1))
        units=[field(np.eye(36)[i],ZeroBackground(base)) for i in indices]
        data=moment_slices(inner,base,current,orders=(k,),unit_fields=units,
            unit_fields_are_deltas=True,radial_breaks=RADIAL_BREAKS,n=24)[0]
        coeff=_integrated_moment_coefficients(data)
        probe=np.random.default_rng(73).normal(size=12)*.01
        np.testing.assert_allclose(_evaluate(coeff,probe),outer_moments(data,probe),rtol=1e-10,atol=1e-7)
        c=.5*(coeff[0][1]-coeff[0][3])
        b=.5*(coeff[1][1]-coeff[1][3])
        Q=.5*(coeff[2][1]-coeff[2][3])
        # Axisymmetric axial momentum has no swirl contribution. Inspect the
        # omitted coefficients explicitly before restricting to poloidal modes.
        swirl=max(np.max(abs(b[:6])),np.max(abs(Q[:6,:])),np.max(abs(Q[:,:6])))
        B=b[6:]; H=Q[6:,6:]
        eigenvalues=np.linalg.eigvalsh(H)
        delta=np.zeros(12)
        lower=None
        if eigenvalues[0]>0:
            delta[6:]=-.5*np.linalg.solve(H,B)
            lower=float(c+.5*B@delta[6:])
        candidate=a.copy();candidate[indices]+=delta
        replay=[]
        for n in (48,96):
            f=field(candidate)
            d=moment_slices(inner,base,f,orders=(k,),unit_fields=[f],radial_breaks=RADIAL_BREAKS,n=n)[0]
            values=outer_moments(d,np.zeros(1))
            replay.append(dict(order=n,moments=values.tolist(),signed_axial=float(.5*(values[1]-values[3]))))
        row=dict(k=k,constant=float(c),linear=B.tolist(),quadratic=H.tolist(),
            poloidal_eigenvalues=eigenvalues.tolist(),swirl_coefficient_max=float(swirl),
            positive_definite=bool(eigenvalues[0]>0),quadratic_minimum=lower,
            minimizing_local_delta=delta.tolist(),direct_replay=replay)
        row.update(range_diagnostics(c,B,H))
        rows.append(row)
        print(json.dumps(row),flush=True)
        (ROOT/'midplane_axial_quadratic_range.json').write_text(json.dumps(dict(scales=rows,
            scope='Numerical finite-basis signed axial-moment range with fixed pressure and fixed axial/time ansatz; not a continuum impossibility proof.',
            accepted=False,scale_recursion_established=False),indent=2)+'\n')

def range_diagnostics(c,b,H):
    bound=40.0
    diagonal_min=[]
    for bi,qi in zip(b,np.diag(H)):
        choices=[-bound,bound]
        if qi>0:
            choices.append(float(np.clip(-bi/(2*qi),-bound,bound)))
        diagonal_min.append(min(qi*x*x+bi*x for x in choices))
    off=H-np.diag(np.diag(H))
    lower=float(c+sum(diagonal_min)-bound**2*np.sum(abs(off)))
    eig,vec=np.linalg.eigh(H)
    ray=None
    if eig[0]<0:
        direction=vec[:,0]
        roots=np.roots([eig[0],b@direction,c])
        ray=dict(direction=direction.tolist(),roots=roots.tolist(),
                 smallest_abs_root=float(min(abs(roots))))
    return dict(poloidal_box_bound=bound,
        signed_axial_box_lower_bound=lower,negative_ray=ray,
        bound_scope='Floating-point quadratic model for poloidal deltas only, swirl fixed; not interval-certified or a global PDE obstruction.')

if __name__=='__main__':
    run()
