"""Necessary bound for nonincreasing positive swirl at frozen endpoints.
The functional inequality is analytic. Profile/tail values remain floating
numerical evaluations, not interval-certified input enclosures.
"""
from common import *
from annulus import Annulus
from heat_exterior import profile
from scipy.optimize import brentq

def bound(xc,xb,Fc,Fb,I):
    if not(0<xc<xb and Fc>Fb>=0):raise ValueError('Ordered positive endpoints required')
    lo=Fb*(xb*xb-xc*xc);hi=Fc*(xb*xb-xc*xc)
    if not(lo<=I<=hi):return {'I_feasible':False,'I_range':[lo,hi]}
    xstar=np.sqrt(xc*xc+(I-lo)/(Fc-Fb))
    cap=Fc*Fc*(xstar-xc)+Fb*Fb*(xb-xstar)
    lam=(Fc+Fb)/(2*xstar)
    return dict(I_feasible=True,xstar=float(xstar),Cp_upper=float(cap),dual_lambda=float(lam),I_range=[lo,hi])

def row(e,c=.25,n=96):
    a=Annulus(c=c,n=24);_,_,inner,target,_=a.setup_eta(e)
    fc=float(a.core.dval('F',a.Xc,e));fb=float(profile(a.Xb,e,c,n=n)['F']);diff=target-inner
    b=bound(a.Xc,a.Xb,fc,fb,float(diff[1]))
    fx=float(a.core.dval('F',a.Xc,e,nx=1));ux=float(a.core.dval('U',a.Xc,e,nx=1))
    aa=-2*a.Xc*fx/fc;bb=2*a.Xc*ux/(np.sqrt(2*a.Xc)*fc);v=aa+bb*bb/aa
    return dict(eta=float(e),heat_amplitude=c,Fc=fc,Fb=fb,I_annulus=float(diff[1]),Cp_required=float(diff[4]),Cp_excess=float(diff[4]-b.get('Cp_upper',np.nan)),inner_shear_a=aa,inner_shear_b=bb,inner_v=v,**b)

def main():
    ee=np.linspace(-.5,.5,101);rows=[row(e) for e in ee]
    threshold=[]
    for e in [0,.25,.5,-.5]:
        rt=brentq(lambda c:row(e,c)['Cp_excess'],.15,.5,xtol=1e-13)
        threshold.append(dict(eta=e,min_amplitude_for_necessary_bound=rt))
    result=dict(rows=rows,thresholds=threshold,positive_gap_count=sum(x['Cp_excess']>0 for x in rows),max_gap=max(x['Cp_excess'] for x in rows),scope='Analytic necessary inequality plus floating input values; not an interval certificate, a universal NS obstruction or sufficiency of a passing bound.')
    save('evidence/necessary_bound.json',result)
    print('eta0',row(0));print('thresholds',threshold);print('violations',result['positive_gap_count'],len(rows))
if __name__=='__main__':main()
