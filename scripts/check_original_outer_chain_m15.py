"""Independent local-germ, varying-J and signed outer-chain fixture.

Trial periodic germs test the coordinate calculus. They are not evaluations
of the original Poisson field. The actual source graph proves its own germs.
"""
import mpmath as mp


def check_fixture():
    c=mp.mp.clone();c.dps=85;N=17;P=c.mpf('2.7');delta=c.mpf('.001');Z=c.mpf('.27')
    b=lambda z:1+z+z**3
    U0=lambda z:4*z
    initial=lambda z:c.mpf('.37')+z+z*z
    leading_initial=lambda z:c.mpf('.11')+z*z
    def U(y,z,phase):return U0(z)+P*b(z)*c.cos(N*y+phase)/N
    def mean(y,z,phase):
        particular=lambda q:P*b(z)/N*(c.cos(N*q+phase)+N*c.sin(N*q+phase))/(1+N*N)
        return U0(z)+particular(y)+(initial(z)-U0(z)-particular(0))*c.exp(-y)
    def lead_mean(y,z):return U0(z)+(leading_initial(z)-U0(z))*c.exp(-y)
    def qy(y,z,phase,leading=False):
        u=U0(z) if leading else U(y,z,phase)
        m=lead_mean(y,z) if leading else mean(y,z,phase)
        uy=c.mpf(0) if leading else c.diff(lambda q:U(q,z,phase),y)
        uz=c.diff(U0,z) if leading else c.diff(lambda q:U(y,q,phase),z)
        mz=c.diff(lambda q:lead_mean(y,q) if leading else mean(y,q,phase),z)
        return (2*z*uy-(1-delta)*z*(u-m)-(1-z*z)*(uz-mz))/(1-delta*z*z)
    identity=lambda x:x
    log_exit=lambda x:c.log(x)-1
    exponential_inlet=lambda x:c.exp(40*x)-1
    exponential_exit=lambda x:c.exp(40*x)-c.exp(40)
    # Each pair is an actual chart type at its common endpoint. Constant-J
    # translations are immaterial to these local derivative checks.
    seams=[('patch_reference',log_exit,c.e,lambda x:x+5,c.mpf(-5)),
           ('reference_slope',identity,c.mpf(0),identity,c.mpf(0)),
           ('slope_axial',lambda x:x-1,c.mpf(1),exponential_inlet,c.mpf(0)),
           ('axial_buffer',exponential_exit,c.mpf(1),identity,c.mpf(0)),
           ('buffer_transition',lambda x:x-11,c.mpf(11),identity,c.mpf(0)),
           ('transition_power',lambda x:x-1,c.mpf(1),identity,c.mpf(0))]
    differences=[];rows=[]
    def native_y_derivative(f,y,x):
        return c.diff(f,x)/c.diff(y,x)
    def native_y_second(f,y,x):
        return c.diff(lambda q:c.diff(f,q)/c.diff(y,q),x)/c.diff(y,x)
    for i,(name,yl,xl,yr,xr) in enumerate(seams):
        phase=c.mpf('.413')+c.mpf('.19')*i;row=[]
        for k in (0,1):
            left=c.diff(lambda z:native_y_derivative(lambda x:U(yl(x),z,phase),yl,xl),Z,k)
            right=c.diff(lambda z:native_y_derivative(lambda x:U(yr(x),z,phase),yr,xr),Z,k)
            physical=c.diff(lambda z:c.diff(lambda y:U(y,z,phase),0),Z,k)
            row.extend((abs(left-right),abs(left-physical)))
        for y,x in ((yl,xl),(yr,xr)):
            row.append(abs(native_y_second(lambda q:U(y(q),Z,phase),y,x)-c.diff(lambda q:U(q,Z,phase),0,2)))
        differences.extend(row);rows.append(dict(interface=name,maximum_difference=c.nstr(max(row),25)))
    phase=c.mpf('.413')
    for k in (0,1):
        differences.append(abs(c.diff(lambda z:c.diff(lambda y:mean(y,z,phase),0),Z,k)-c.diff(lambda z:U(0,z,phase)-initial(z),Z,k)))
    # Exact power-zone source: U=Uy=0, both means keep nonzero inlet memory.
    pm=lambda y,z:initial(z)*c.exp(-y)
    pl=lambda y,z:leading_initial(z)*c.exp(-y)
    def power_D(y,z):
        d=pm(y,z)-pl(y,z);dz=c.diff(lambda q:pm(y,q)-pl(y,q),z)
        return c.mpf('7.3')*c.exp(y)*((1-delta)*z*d+(1-z*z)*dz)/(1-delta*z*z)
    for k in (0,1):
        differences.append(abs(c.diff(lambda z:power_D(0,z)-power_D(2,z),Z,k)))
    # Arbitrary shared endpoint functions are legitimate inputs to the
    # additive telescoping identity. None is selected as an original value.
    trace=lambda i,z:c.mpf('7.3')*c.exp(c.mpf('.13')*i)*(qy(0,z,phase+c.mpf('.19')*i)-qy(0,z,phase+c.mpf('.19')*i,True))
    chain=[]
    for k in (0,1):
        signed=sum(c.diff(lambda z:trace(i,z)-trace(i+1,z),Z,k) for i in range(7))
        external=c.diff(lambda z:trace(0,z)-trace(7,z),Z,k)
        error=abs(signed-external);differences.append(error);chain.append(c.nstr(error,25))
    correct=c.diff(lambda y:U(y,Z,phase),0)
    wrong_fast=-P*b(Z)*c.sin(phase)/N
    correct_second=c.diff(lambda y:U(y,Z,phase),0,2)
    wrong_second=c.diff(lambda x:U(exponential_inlet(x),Z,phase),0,2)/(40*40)
    external=trace(0,Z)-trace(7,Z)
    controls=dict(extra_N_on_fast_term_rejected=abs(correct-wrong_fast)>c.mpf('.01'),
        frozen_J_second_derivative_rejected=abs(correct_second-wrong_second)>c.mpf('.01'),
        phase_reset_rejected=abs(U(0,Z,phase)-U(0,Z,0))>c.mpf('.01'),
        nonzero_incoming_mean=abs(initial(Z))>c.mpf('.01'),
        power_mean_not_zero=abs(pm(1,Z))>c.mpf('.01'),
        power_D_not_zero=abs(power_D(0,Z))>c.mpf('.01'),
        missing_physical_R_rejected=abs(power_D(0,Z)-power_D(0,Z)/c.mpf('7.3'))>c.mpf('.01'),
        wrong_telescoping_sign_rejected=abs(2*external)>c.mpf('.01'),
        inlet_retained=abs(trace(0,Z))>c.mpf('.01'),
        outlet_retained=abs(trace(7,Z))>c.mpf('.01'),
        leading_subtraction_retained=abs(qy(0,Z,phase,True))>c.mpf('.01'),
        genuine_ordinary_Z1=abs(c.diff(lambda z:power_D(0,z),Z))>c.mpf('.01'))
    maximum=max(differences)
    return dict(passed=maximum<c.mpf('1e-38') and all(controls.values()),finite_N=N,
        synthetic_local_germs_only=True,original_source_values_evaluated=False,
        original_Poisson_inverse_evaluated=False,NS_residual=False,
        interfaces=rows,interface_count=6,power_transport_with_nonzero_mean=True,
        maximum_supporting_difference=c.nstr(maximum,28),signed_chain_C0_Z1_differences=chain,
        controls=controls,control_count=len(controls))


if __name__=='__main__':
    import json
    result=check_fixture();print(json.dumps(result,indent=2));raise SystemExit(0 if result['passed'] else 1)
