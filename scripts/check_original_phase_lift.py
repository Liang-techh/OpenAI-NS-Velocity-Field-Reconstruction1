"""Independent moderate-parameter fixtures for the smooth periodic loop."""
import json
import mpmath as mp


def run():
    c=mp.mp.clone();c.dps=40;pi=c.pi
    rows=[]
    for q in (c.mpf('0'),c.mpf('.6'),c.mpf('-.4')):
        a,E,b,p2,d=map(c.mpf,('.8','1.3','.25','.7','.9'))
        t0=-b/a;u=p2*q/d;r=u/c.sqrt(1+u*u);alpha=2*q/c.sqrt(1+u*u)
        K=1/(2*pi*(1+t0*t0+2*q*q))
        t=lambda x:t0+alpha*(c.cos(x)-r)/(1-2*r*c.cos(x)+r*r)
        I1=c.quad(t,[0,pi,2*pi]);I2=c.quad(lambda x:t(x)**2,[0,pi,2*pi])
        assert abs(I1-2*pi*t0)<c.mpf('1e-34')
        assert abs(I2-2*pi*(t0*t0+2*q*q))<c.mpf('1e-34')
        phase=lambda x:K*(x+c.quad(lambda z:t(z)**2,[0,x]))
        phi=c.mpf('.31')
        psi=c.findroot(lambda x:phase(x)-phi,(pi/2,pi),solver='secant')
        assert 0<psi<2*pi and abs(phase(psi)-phi)<c.mpf('1e-34')
        T=c.quad(t,[0,psi]);Tplus=c.quad(t,[0,pi,2*pi,2*pi+psi])
        A=lambda f,s:a*(f-s/(2*pi))/2
        B=lambda f,v:E*(-a*v/(2*pi)-b*f)/2
        assert abs(A(phi+1,psi+2*pi)-A(phi,psi))<c.mpf('1e-34')
        assert abs(B(phi+1,Tplus)-B(phi,T))<c.mpf('1e-34')
        wrong_B=lambda f,v:E*(-a*v/(2*pi)-(b+c.mpf('.1'))*f)/2
        assert abs(wrong_B(phi+1,Tplus)-wrong_B(phi,T))>c.mpf('.01')
        assert abs(2*K*(2*pi+I2)-1)>c.mpf('.9')
        if not q:assert abs(A(phi,psi))+abs(B(phi,T))<c.mpf('1e-34')
        rows.append(dict(q=str(q),period_and_shift=True,arbitrary_phase=str(phi),
            wrong_root_relation_rejected=True,wrong_K_rejected=True))
    return dict(passed=True,synthetic_parameters_only=True,actual_phase_value_selected=False,
        rows=rows,q_zero_included=True,negative_q_included=True,arbitrary_fast_phase=True)


if __name__=='__main__':print(json.dumps(run(),indent=2))
