"""Finite rational lower sums independently checked against kernel-tail majorants."""
def run(leaf,c):
    from fractions import Fraction as Q
    from math import factorial
    def end(t):
        sign,mantissa,exponent,_=t; v=Q(int(mantissa))*Q(2)**int(exponent);return -v if sign else v
    out=leaf.build(c,a=c.mpf(2),chi_z1=c.mpf([-2,2]),chi_z2=c.mpf([-3,3])); sums=[[Q(0) for _ in range(3)] for _ in range(2)];r=Q(2)
    for k in range(2,42):
        t=r**k/(factorial(k)*factorial(k+1)); s1=k*2; s2=k*3+k*(k-1)*4
        for row,m in ((0,Q(1)),(1,Q(k))):
            sums[row][0]+=m*t;sums[row][1]+=m*t*s1;sums[row][2]+=m*t*s2
    checks={}
    for row,field in enumerate(('F1','d_x_F1')):
        for z in range(3): checks[field+'_Z'+str(z)]=end(out[field][z]._mpi_[1])>=sums[row][z]
    checks['tail_only']=out['source_or_feedback_included'] is False and out['bounds_are_not_field_values'] is True
    return {'passed':all(checks.values()),'check_count':len(checks),'checks':checks,'synthetic_only':True,'reference':'exact rational positive sums through j=41'}
