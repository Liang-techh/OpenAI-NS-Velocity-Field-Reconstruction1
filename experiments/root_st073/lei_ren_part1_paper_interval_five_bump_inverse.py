"""Uniform directed contraction for the five-bump map on an axial family.

Certifies an implicit C1 coefficient family for the enclosed true map/data.
No midpoint data is substituted in the map; midpoints only choose a fixed
preconditioner. No repaired physical field is installed by this module.
"""
import hashlib
import json
from fractions import Fraction
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value,restore_jet
from lei_ren_part1_paper_interval_comparison_jets import IntervalComparisonJets
from lei_ren_part1_paper_interval_exit_continuation_enclosure import _pack as pack
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
HERE=Path(__file__).parent


def mag(c,v):
    lo,hi=endpoints(v)
    return c.mpf(max(abs(lo),abs(hi)))


def dot(row,v):
    return sum((a*b for a,b in zip(row,v)),row[0]*0)


def matmul(A,B):
    return [[dot(row,[B[k][j] for k in range(len(B))]) for j in range(len(B[0]))] for row in A]


def weights(c,raw):
    table={(Fraction(row['center']),Fraction(row['power']),row['multiplicity']):restore_value(c,row['weight_interval'])
        for row in raw['weight_records'].values()}
    centers=(Fraction(5,4),Fraction(3,2),Fraction(7,4))
    def w(i,p,k=1):return table[centers[i],Fraction(p),k]
    L=[[c.mpf(0) for j in range(5)] for i in range(5)]
    L[0][0]=L[0][1]=c.mpf(1)
    for j,i in enumerate((0,2)):L[1][j]=w(i,'.6')
    for i in range(3):
        L[2][i+2]=w(i,'.5');L[3][i+2]=-w(i,'.1');L[4][i+2]=w(i,'-.9')
    return dict(L=L,fg=[w(0,'.5',2),w(2,'.5',2)],gg=[w(0,'0',2),w(2,'0',2)],
        ff=[w(i,'0',2) for i in range(3)],ff_over_x=[w(i,'-1',2) for i in range(3)])


def nonlinear(c,W,h,invAm2):
    c1,c2,x1,x2,x3=h;cs=(c1,c2);xs=(x1,x2,x3)
    return [c.mpf(0),W['fg'][0]*c1*x1+W['fg'][1]*c2*x3,c.mpf(0),
        invAm2*sum((W['gg'][i]*cs[i]*cs[i] for i in range(2)),c.mpf(0))
        -sum((W['ff'][i]*xs[i]*xs[i] for i in range(3)),c.mpf(0))/2,
        sum((W['ff_over_x'][i]*xs[i]*xs[i] for i in range(3)),c.mpf(0))/2]


def nonlinear_jacobian(c,W,h,invAm2):
    J=[[c.mpf(0) for j in range(5)] for i in range(5)]
    J[1][0]=W['fg'][0]*h[2];J[1][1]=W['fg'][1]*h[4]
    J[1][2]=W['fg'][0]*h[0];J[1][4]=W['fg'][1]*h[1]
    for j in range(2):J[3][j]=invAm2*W['gg'][j]*h[j]*2
    for i in range(3):
        J[3][i+2]=-W['ff'][i]*h[i+2]
        J[4][i+2]=W['ff_over_x'][i]*h[i+2]
    return J


def identity_error(c,R,J):
    product=matmul(R,J)
    return [[c.mpf(int(i==j))-product[i][j] for j in range(5)] for i in range(5)]


def weighted_norm(c,A,scales):
    return max((sum((mag(c,A[i][j])*scales[j]/scales[i] for j in range(5)),c.mpf(0))
        for i in range(5)),key=lambda v:endpoints(v)[1])


def intersection(c,a,b):
    al,ah=endpoints(a);bl,bh=endpoints(b)
    if max(al,bl)>min(ah,bh):raise ValueError('fixed-point interval intersection empty')
    return c.mpf([max(al,bl),min(ah,bh)])


def certify(c,W,d,invAm2,scales=None,tightening=12):
    if scales is None:scales=[c.mpf('2e-14')]*2+[c.mpf('1e-32')]*3
    L=W['L']
    midpoint=mp.matrix([[(endpoints(v)[0]+endpoints(v)[1])/2 for v in row] for row in L])
    inv=mp.inverse(midpoint)
    R=[[c.mpf(inv[i,j]) for j in range(5)] for i in range(5)]
    A=identity_error(c,R,L)
    linear_error=weighted_norm(c,A,scales)
    if endpoints(linear_error)[1]>=1:raise ValueError('preconditioner nonsingularity not certified')
    X=[c.mpf([-endpoints(s)[1],endpoints(s)[1]]) for s in scales]
    Jq=nonlinear_jacobian(c,W,X,invAm2[0])
    J=[[L[i][j]+Jq[i][j] for j in range(5)] for i in range(5)]
    K=identity_error(c,R,J)
    contraction=weighted_norm(c,K,scales)
    def image(box):
        q=nonlinear(c,W,box,invAm2[0])
        forcing=[d[i][0]+q[i] for i in range(5)]
        return [dot(A[i],box)-dot(R[i],forcing) for i in range(5)]
    initial_image=image(X)
    included=all(endpoints(X[i])[0]<endpoints(initial_image[i])[0] and
        endpoints(initial_image[i])[1]<endpoints(X[i])[1] for i in range(5))
    passed=included and endpoints(contraction)[1]<1
    if not passed:
        return dict(certified=False,self_map_strictly_inside=included,contraction_bound=contraction,
            linear_preconditioner_error=linear_error,initial_box=X,initial_image=initial_image)
    h=X
    for _ in range(tightening):h=[intersection(c,a,b) for a,b in zip(h,image(h))]
    # h_Z solves J h_Z = -d_Z - partial_Z Q. Only row4 depends on Am.
    qz=[c.mpf(0)]*5
    qz[3]=invAm2[1]*sum((W['gg'][i]*h[i]*h[i] for i in range(2)),c.mpf(0))
    bz=[-dot(row,[d[i][1]+qz[i] for i in range(5)]) for row in R]
    Jq=nonlinear_jacobian(c,W,h,invAm2[0])
    K=identity_error(c,R,[[L[i][j]+Jq[i][j] for j in range(5)] for i in range(5)])
    q=weighted_norm(c,K,scales)
    if endpoints(q)[1]>=1:raise ValueError('derivative inverse contraction not certified')
    bnorm=max((mag(c,bz[i])/scales[i] for i in range(5)),key=lambda v:endpoints(v)[1])
    radius=bnorm/(1-q)
    hz=[c.mpf([-endpoints(radius*si)[1],endpoints(radius*si)[1]]) for si in scales]
    for _ in range(tightening):hz=[intersection(c,hz[i],bz[i]+dot(K[i],hz)) for i in range(5)]
    residual=[dot(L[i],h)+nonlinear(c,W,h,invAm2[0])[i]+d[i][0] for i in range(5)]
    differentiated=[dot(Jrow,hz)+d[i][1]+qz[i] for i,Jrow in enumerate(
        [[L[i][j]+Jq[i][j] for j in range(5)] for i in range(5)])]
    return dict(certified=True,uniform_implicit_C1_family_exists=True,unique_in_initial_box=True,
        self_map_strictly_inside=True,contraction_bound=contraction,derivative_contraction_bound=q,
        linear_preconditioner_error=linear_error,scales=scales,initial_box=X,initial_image=initial_image,
        controls=[IntervalTaylor(c,[h[i],hz[i]]) for i in range(5)],
        enclosure_residual=residual,derivative_enclosure_residual=differentiated,
        residual_zero_containment_only_diagnostic=True,closure_proof='strict self-map plus uniform contraction and implicit differentiation',
        preconditioner_only_uses_midpoints=True,actual_map_or_defects_projected_to_midpoints=False,
        field_installed=False,terminal_physical_five_moment_closure=False)


def run():
    calc=IntervalComparisonJets(4);c=calc.ctx
    dn='lei_ren_part1_paper_interval_functional_defects.json'
    wn='lei_ren_part1_paper_bump_integral_enclosures_check.json'
    raw=json.loads((HERE/dn).read_text());wraw=json.loads((HERE/wn).read_text())
    if raw['state_sha256']!=calc.state_hash:raise ValueError('functional defect core mismatch')
    for name,digest in raw['input_hashes'].items():
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('defect dependency changed:'+name)
    for name,digest in (('lei_ren_part1_paper_bump_integral_enclosures.py',wraw['module_sha256']),
        ('lei_ren_part1_paper_bump_integral_enclosures_check.py',wraw['check_source_sha256'])):
        if hashlib.sha256((HERE/name).read_bytes()).hexdigest()!=digest:raise ValueError('weight dependency changed:'+name)
    with mp.workdps(calc.precision+60):
        d=[restore_jet(c,raw['rows'][str(i)],order=1) for i in range(1,6)]
        z=calc.z.truncate(1);Am=(1+z*z).reciprocal()*c.exp(c.mpf('13.4'))
        invAm2=(Am*Am).reciprocal()
        result=certify(c,weights(c,wraw),d,invAm2)
        result.update(center_family=calc.center_family,state_sha256=calc.state_hash,
            coefficient_order=['c1','c2','xi1','xi2','xi3'],fixed_accepted_parameters=True,
            whole_axis_inverse=False,uniform_full_field_cone=False,temporal_recursion=False)
    result['input_hashes']={name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in
        (Path(__file__).name,dn,wn,'lei_ren_part1_paper_interval_functional_defects.py',
        'lei_ren_part1_paper_five_bump_map.py')}
    Path(__file__).with_suffix('.json').write_text(json.dumps(encode(pack(result)),indent=2)+'\n',encoding='utf-8')
    print('Uniform five-bump inverse certificate:',result['certified'],
        'contraction bound',[mp.nstr(v,12) for v in endpoints(result['contraction_bound'])],flush=True)
    return result

if __name__=='__main__':run()
