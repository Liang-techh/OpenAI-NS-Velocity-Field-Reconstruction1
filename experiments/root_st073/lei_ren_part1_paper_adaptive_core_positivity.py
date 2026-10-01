"""Resumable whole-axis Bernstein positivity enclosure for the finite core.

Accepted cells and pending cells partition [-1,1]. A run limit leaves the
proof incomplete, not false; every saved endpoint is an exact MP tuple.
"""
import argparse,json,time,hashlib
from pathlib import Path
import mpmath as mp
from mpmath.ctx_iv import MPIntervalContext
from lei_ren_part1_paper_uniform_axis_jets import uniform_axis_jets
from lei_ren_part1_paper_amplitude_factored_core import factored_core_coefficients
from lei_ren_part1_paper_factored_core_positivity import bernstein_enclosure,squared_axis_rows
from lei_ren_part1_paper_global_finite_core_bounds import accepted_pressure_axis_rows
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints


def pack(value):return list(value._mpf_)
def unpack(value):return mp.make_mpf(tuple(value))


def initial_cells():
    j=mp.mpf('1e-14');sigma=j/500
    root=mp.findroot(lambda z:(1-mp.mpf('1e-200'))*z/2+(1-z*z)*(4*z+j),(-j,mp.mpf(0)))
    edges={mp.mpf(-1),mp.mpf(1),root};distance=sigma
    while distance<3:
        for v in (root-distance,root+distance):
            if -1<v<1:edges.add(v)
        distance*=2
    edges=sorted(edges)
    return [dict(left=pack(l),right=pack(r),depth=0) for l,r in zip(edges,edges[1:])]


def coverage(state):
    cells=state['accepted']+state['pending']+state['unresolved'];cells=sorted(cells,key=lambda c:unpack(c['left']))
    edge=mp.mpf(-1)
    for cell in cells:
        if unpack(cell['left'])!=edge:raise AssertionError('Axial partition gap or overlap')
        edge=unpack(cell['right'])
    if edge!=1:raise AssertionError('Incomplete partition')
    length=sum((unpack(c['right'])-unpack(c['left']) for c in state['accepted']),mp.mpf(0))
    return length/2


def run(seconds=180,max_depth=24):
    base=Path(__file__).parent;path=Path(__file__).with_suffix('.json');precision=260
    files=[Path(__file__),base/'lei_ren_part1_paper_amplitude_factored_core.py',
        base/'lei_ren_part1_paper_factored_core_positivity.py',base/'lei_ren_part1_paper_uniform_axis_jets.py',
        base/'lei_ren_part1_paper_global_finite_core_bounds.py',base/'lei_ren_part1_paper_global_pressure_high_derivatives.json']
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    ctx=MPIntervalContext();ctx.dps=precision
    with mp.workdps(precision+40):
        if path.exists():
            state=json.loads(path.read_text())
            if state['input_hashes']!=hashes:raise ValueError('Saved proof inputs changed; preserve receipt and explicitly rebuild')
        else:
            state=dict(input_hashes=hashes,precision=precision,radial_degree=18,axis_domain=['-1','1'],
                normalized_radial_domain=['0','1'],accepted=[],pending=initial_cells(),unresolved=[],evaluations=0,
                infinite_radial_remainder_enclosed=False,full_K_certified=False,temporal_recursion=False)
        start=time.monotonic();degree=18;length=degree+4;lam=ctx.mpf('1e36');r=4/lam
        p0,alignment=accepted_pressure_axis_rows(ctx,ctx.mpf([-1,1]),length,files[-1])
        p0=[v.value for v in p0];state['accepted_schedule_sha256']=alignment['accepted_schedule']['sha256']
        def save():
            fraction=coverage(state)
            state['accepted_axial_length_fraction']=mp.nstr(fraction,30)
            state['whole_axis_finite_core_positivity_certified']=not state['pending'] and not state['unresolved']
            state['minimum_positive_lower_exact']=min((c['lower'] for c in state['accepted']),key=unpack) if state['accepted'] else None
            path.write_text(json.dumps(state,indent=2)+'\n')
        while state['pending'] and time.monotonic()-start<seconds:
            cell=state['pending'].pop();left=unpack(cell['left']);right=unpack(cell['right']);z=ctx.mpf([left,right])
            axis=uniform_axis_jets(ctx,radius=1,j='1e-14',Lambda=lam,logC='5e151',delta='1e-200',length=length,axial_interval=z)
            gradient=axis['gradient_coefficients']
            core=factored_core_coefficients(z,'1e-200',ell_Z_taylor=[-lam*v for v in gradient],
                S_Z_taylor=squared_axis_rows(ctx,gradient,lam,axis['F0_interval'],length),U0_Z_taylor=axis['U0'],
                P0_Z_taylor=p0,radial_degree=degree,precision=precision,scalar_converter=ctx.mpf)
            bound=bernstein_enclosure(ctx,[row[0]*r**n for n,row in enumerate(core['A'])]);lo,hi=endpoints(bound['range'])
            state['evaluations']+=1
            if lo>0:
                cell.update(lower=pack(lo),upper=pack(hi));state['accepted'].append(cell)
            elif cell['depth']>=max_depth:
                cell.update(last_lower=pack(lo),reason='inconclusive interval at refinement limit');state['unresolved'].append(cell)
            else:
                mid=(left+right)/2;depth=cell['depth']+1
                state['pending'].extend([dict(left=pack(left),right=pack(mid),depth=depth),dict(left=pack(mid),right=pack(right),depth=depth)])
            if state['evaluations']%10==0:
                save();print('evaluations',state['evaluations'],'accepted',len(state['accepted']),
                    'pending',len(state['pending']),'covered length fraction',state['accepted_axial_length_fraction'],flush=True)
        save();print('saved accepted',len(state['accepted']),'pending',len(state['pending']),
                     'unresolved',len(state['unresolved']),'whole-axis certificate',state['whole_axis_finite_core_positivity_certified'],flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--seconds',type=float,default=180);parser.add_argument('--max-depth',type=int,default=24)
    args=parser.parse_args();run(args.seconds,args.max_depth)
