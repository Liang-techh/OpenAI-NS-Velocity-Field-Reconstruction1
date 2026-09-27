"""Check actual residual-stress geometry near the growing wave center."""
import json
import numpy as np
from fourier_shear_feasibility import load_saved_field
from midplane_outer_residual_source import outer_cones
from joined_field import coordinates
from radial_continuation import ROOT


def run():
    field,repair=load_saved_field()
    seed=json.loads((ROOT/'broad_shear_growth.json').read_text(encoding='utf-8'))
    r,z=seed['center']
    k=repair['k'];tau=.5*2.**-k
    co=coordinates(r/np.sqrt(field.nu),z/np.sqrt(field.nu),tau,field.inner.h)
    eta=float(co['eta'])
    ri=np.sqrt(2*field.nu*float(co['q'])*field.join_X)
    y=(r/ri-1)/(field.ratio-1)
    locations=[(eta,y),(eta-.03,y),(eta+.03,y),(eta,y-.07),(eta,y+.07)]
    rows=outer_cones(field,k,order=64,radial_breaks=repair['radial_breaks'],locations=locations)
    report=dict(accepted=False,source='fourier_shear_feasibility.json',k=k,
                center_similarity=dict(eta=eta,y=y),rows=rows,
                pass_count=sum(row['cone_pass'] for row in rows),
                scope='Five full-residual stress-cone samples near the growing wave center; distinct from the 22 outer controller nodes. No support-wide certificate.')
    (ROOT/'broad_shear_wave_cone.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(report),flush=True)
    return report


if __name__=='__main__':run()
