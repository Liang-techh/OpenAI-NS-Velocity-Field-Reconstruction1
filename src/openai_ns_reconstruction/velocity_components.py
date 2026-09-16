"""User-facing Cartesian velocity components for the frozen coupled candidate.

Coordinates/time are dimensionless model variables, not calibrated SI units.
The available time interval is [.25,.75]. Visual correspondence is not verified.
"""
from functools import lru_cache
import hashlib
import json
from pathlib import Path
import numpy as np
from .constrained_coupled import CoupledCandidate

DEFAULT_CANDIDATE=Path(__file__).resolve().parent/'data'/'velocity_candidate.json'


class VelocityField:
    def __init__(self,candidate_path=DEFAULT_CANDIDATE):
        self.path=Path(candidate_path).resolve()
        self.sha256=hashlib.sha256(self.path.read_bytes()).hexdigest()
        self.candidate=CoupledCandidate.load(self.path)

    def components(self,x,y,z,t):
        """Return (u,v,w), broadcasting coordinates/time like NumPy.

        Scalars return Python floats; array inputs return matching broadcast arrays.
        The exact compact zero exterior is handled before candidate evaluation.
        """
        x,y,z,t=np.broadcast_arrays(*[np.asarray(a,dtype=float) for a in (x,y,z,t)])
        if not all(np.all(np.isfinite(a)) for a in (x,y,z,t)):raise ValueError('coordinates and time must be finite')
        if np.any((t<.25)|(t>.75)):raise ValueError('available model time is [.25,.75]; no extrapolation')
        shape=x.shape;points=np.column_stack((x.ravel(),y.ravel(),z.ravel()));times=t.ravel()
        result=np.zeros((len(points),3))
        indices=np.flatnonzero((np.hypot(points[:,0],points[:,1])<2)&(np.abs(points[:,2])<2))
        for start in range(0,len(indices),4096):
            ids=indices[start:start+4096];result[ids]=self.candidate.velocity(points[ids],times[ids])
        if not np.all(np.isfinite(result)):raise FloatingPointError('candidate returned nonfinite velocities')
        values=tuple(result[:,i].reshape(shape) for i in range(3))
        return tuple(float(a) for a in values) if not shape else values

    def at_points(self,points,time):
        """Return (...,3) vectors for (...,3) Cartesian positions."""
        p=np.asarray(points,dtype=float)
        if p.ndim<1 or p.shape[-1]!=3:raise ValueError('points require shape (...,3)')
        return np.stack(self.components(*np.moveaxis(p,-1,0),time),axis=-1)

    def grid(self,x,y,z,times):
        """Return grid vectors in (time,x,y,z,component) order."""
        axes=[np.asarray(a,dtype=float) for a in (times,x,y,z)]
        if any(a.ndim!=1 or not a.size for a in axes):raise ValueError('grid axes and times must be nonempty1D arrays')
        tt,xx,yy,zz=np.meshgrid(*axes,indexing='ij')
        return np.stack(self.components(xx,yy,zz,tt),axis=-1)

    def metadata(self):
        return {'family':'coupled_velocity_v1','candidate_sha256':self.sha256,
                'components':['u','v','w'],'coordinates':'right-handed Cartesian x,y,z',
                'units':'dimensionless; no physical unit or OpenAI scene calibration claimed',
                'time_interval':[.25,.75],'compact_support':'x²+y²<4 and |z|<2; exact zero exterior',
                'visual_correspondence':'not verified','pde_acceptance':'failed; see coupled_joint validation',
                'grid_layout':['time','x','y','z','component']}


@lru_cache(maxsize=1)
def default_field():return VelocityField()

def velocity(x,y,z,t):return default_field().components(x,y,z,t)
def u(x,y,z,t):return velocity(x,y,z,t)[0]
def v(x,y,z,t):return velocity(x,y,z,t)[1]
def w(x,y,z,t):return velocity(x,y,z,t)[2]


def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate',type=Path,default=DEFAULT_CANDIDATE)
    parser.add_argument('--point',nargs=4,type=float,metavar=('X','Y','Z','T'))
    parser.add_argument('--export',type=Path,help='write compressed grid.npz and metadata.json here')
    parser.add_argument('--grid-size',type=int,default=17)
    args=parser.parse_args();field=VelocityField(args.candidate)
    if args.point is not None:print(json.dumps(dict(zip(('u','v','w'),field.components(*args.point))),indent=2))
    if args.export:
        if not 2<=args.grid_size<=129:parser.error('grid-size must be between2 and129')
        axes=np.linspace(-2,2,args.grid_size);times=np.array([.25,.5,.75]);values=field.grid(axes,axes,axes,times)
        args.export.mkdir(parents=True,exist_ok=True)
        np.savez_compressed(args.export/'grid.npz',x=axes,y=axes,z=axes,times=times,u=values[...,0],v=values[...,1],w=values[...,2])
        meta=field.metadata();meta['grid_shape']=list(values.shape);meta['grid_npz_sha256']=hashlib.sha256((args.export/'grid.npz').read_bytes()).hexdigest()
        (args.export/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
        print(str(args.export/'grid.npz'))
    if args.point is None and args.export is None:print(json.dumps(field.metadata(),indent=2))

if __name__=='__main__':main()
