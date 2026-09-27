"""Unified global diagnostic field; no NS trajectory or recursion acceptance."""
import hashlib
import json
from pathlib import Path
import numpy as np
from global_axial_extension import build_candidate
from global_collar_tangent import CollarCorrection, BOXES
from enriched_shape_replay import build_field
from endpoint_acceleration_projection import AccelerationCorrectionField, _unpack_acceleration_control

ROOT=Path(__file__).resolve().parent


def load():
    paths={name:ROOT/file for name,file in dict(
        balanced='balanced_refined_tangent.json',acceleration='endpoint_acceleration_projection.json',
        collar='global_collar_tangent.json').items()}
    raw={name:path.read_bytes() for name,path in paths.items()}
    reports={name:json.loads(value) for name,value in raw.items()}
    if any(r['status']!='completed' for r in reports.values()):
        raise ValueError('Only frozen completed source reports may be assembled')
    if reports['acceleration']['sources']['balanced']['sha256']!=hashlib.sha256(raw['balanced']).hexdigest():
        raise ValueError('Acceleration parent does not match balanced source')
    global_base,localized,_,_,snapshot,_=build_candidate(paths['balanced'])
    tau0=snapshot['inputs']['mean']['tau']
    if reports['collar']['tau']!=tau0 or reports['acceleration']['inputs']['tau0']!=tau0:
        raise ValueError('Reference times differ')
    g=snapshot['inputs']['wave'];c,w=np.asarray(g['center']),np.asarray(g['widths'])
    if not all(b[1]<c[0]-w[0] or b[0]>c[0]+w[0] or b[3]<c[1]-w[1] or b[2]>c[1]+w[1] for b in BOXES):
        raise ValueError('Collar correction intersects the inner wave patch')
    exterior=CollarCorrection(global_base,reports['collar']['control'],tau0)
    blocks=_unpack_acceleration_control(reports['acceleration']['ridge_fit']['coefficients'])
    # The unpacker returns mode -> (acceleration, pressure slope).
    acceleration={mode:values[0] for mode,values in blocks.items()}
    pressure={mode:values[1] for mode,values in blocks.items()}
    full=AccelerationCorrectionField(exterior,g['center'],g['widths'],g['carrier'],acceleration,pressure,tau0)
    return full,global_base,localized,snapshot,reports,{name:hashlib.sha256(value).hexdigest() for name,value in raw.items()}


def run():
    full,base,localized,snapshot,reports,hashes=load()
    tau=snapshot['inputs']['mean']['tau'];g=snapshot['inputs']['wave']
    local,_=build_field(reports['balanced'])
    blocks=_unpack_acceleration_control(reports['acceleration']['ridge_fit']['coefficients'])
    inner=AccelerationCorrectionField(local,g['center'],g['widths'],g['carrier'],
        {m:v[0] for m,v in blocks.items()},{m:v[1] for m,v in blocks.items()},tau)
    oldglobal,*_=build_candidate()
    outer=CollarCorrection(oldglobal,reports['collar']['control'],tau)
    c,w=np.asarray(g['center']),np.asarray(g['widths'])
    inner_points=np.array([[c[0],0,c[1]],[c[0],0,c[1]+.5*w[1]],[0,c[0],c[1]]])
    outer_points=np.array([[(b[0]+b[1])/2,0,(b[2]+b[3])/2] for b in BOXES])
    axes=np.array([[0.,0.,0.],[0.,0.,.0007],[0.,0.,.002]])
    report=dict(status='running',accepted=False,pde_validated=False,scale_recursion_established=False,
        scope='Assembly and support checks only. Constituent momentum metrics cover different spatial domains; no whole-domain or interval acceptance.',
        source_hashes=hashes,rows=[])
    for dk in (0.,1e-6):
        t=tau*2**(-dk)
        iu,ip=full.fields(inner_points,t);ru,rp=inner.fields(inner_points,t)
        ou,op=full.fields(outer_points,t);eu,ep=outer.fields(outer_points,t)
        au,ap=full.fields(axes,t)
        if not np.isfinite(au).all() or not np.isfinite(ap).all():
            raise ValueError('Nonfinite field on the symmetry axis')
        row=dict(delta_k=dk,tau=t,inner_velocity_error=float(np.max(abs(iu-ru))),
                 inner_pressure_error=float(np.max(abs(ip-rp))),outer_velocity_error=float(np.max(abs(ou-eu))),
                 outer_pressure_error=float(np.max(abs(op-ep))),axis_velocity=au.tolist(),axis_pressure=ap.tolist())
        report['rows'].append(row)
        print(json.dumps(row),flush=True)
    report['status']='completed'
    (ROOT/'global_acceleration_candidate.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    run()
