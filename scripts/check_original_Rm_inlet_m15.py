"""Independent synthetic varying-J and signed-trace bookkeeping fixture."""
import mpmath as mp


def check_fixture():
    c=mp.mp.clone();c.dps=65;N=17;phase=c.mpf('.413');z=c.mpf('.27')
    P=c.mpf('2.7');R0=c.mpf('7.3');gap=c.mpf('3.2')
    b=lambda Z:1+Z+Z**3
    rate=lambda Z:c.mpf('.31')+Z**2
    V0=lambda Z:4*Z/P
    B=lambda y,Z,phi:b(Z)*c.exp(rate(Z)*y)*c.cos(phi)
    V=lambda y,Z:V0(Z)+B(y,Z,N*y+phase)/N
    m0=lambda Z:c.mpf('.37')+Z+Z**2
    def primitive(y,Z):
        k=1+rate(Z)
        return c.exp(k*y)*(k*c.cos(N*y+phase)+N*c.sin(N*y+phase))/(k*k+N*N)
    def mean(y,Z):
        return V0(Z)+(m0(Z)-V0(Z))*c.exp(-y)+b(Z)/N*c.exp(-y)*(primitive(y,Z)-primitive(0,Z))
    def Qy(y,Z):
        U=P*V(y,Z);M=P*mean(y,Z)
        Uy=P*c.diff(lambda s:V(s,Z),y)
        Uz=c.diff(lambda t:P*V(y,t),Z);Mz=c.diff(lambda t:P*mean(y,t),Z)
        d=c.mpf('.001')
        return (2*Z*Uy-(1-d)*Z*(U-M)-(1-Z*Z)*(Uz-Mz))/(1-d*Z*Z)
    def leading_Qy(y,Z):
        leading_mean=lambda Y,T:V0(T)+(c.mpf('.11')+T*T-V0(T))*c.exp(-Y)
        U=P*V0(Z);M=P*leading_mean(y,Z)
        Uz=P*c.diff(V0,Z);Mz=P*c.diff(lambda T:leading_mean(y,T),Z)
        d=c.mpf('.001')
        return (-(1-d)*Z*(U-M)-(1-Z*Z)*(Uz-Mz))/(1-d*Z*Z)
    ypost=lambda q:gap*(q-1)
    yrm=lambda x:c.log(x)
    postU=lambda q,Z:P*V(ypost(q),Z)
    rmU=lambda x,Z:P*V(yrm(x),Z)
    differences=[]
    for k in (0,1):
        post=c.diff(lambda Z:c.diff(lambda q:postU(q,Z),1)/gap,z,k)
        right=c.diff(lambda Z:c.diff(lambda x:rmU(x,Z),1),z,k)
        independent=c.diff(lambda Z:P*c.diff(lambda y:V(y,Z),0),z,k)
        differences.extend((abs(post-right),abs(right-independent)))
    native2=c.diff(lambda x:x*c.diff(lambda xx:rmU(xx,z),x),1)
    physical2=P*c.diff(lambda y:V(y,z),0,2)
    differences.append(abs(native2-physical2))
    for k in (0,1):
        differences.append(abs(c.diff(lambda Z:c.diff(lambda y:mean(y,Z),0),z,k)-c.diff(lambda Z:V(0,Z)-m0(Z),z,k)))
    D=lambda y,Z:R0*c.exp(y)*(Qy(y,Z)-leading_Qy(y,Z))
    # Actual graph values are not used here. Five synthetic source pieces
    # retain two nonzero exterior traces and four exact internal seams.
    endpoints=[c.mpf(v) for v in ('-.7','-.4','-.2','-.1','0','.3')]
    chain=[]
    for k in (0,1):
        signed=sum(c.diff(lambda Z:D(a,Z)-D(bb,Z),z,k) for a,bb in zip(endpoints[:-1],endpoints[1:]))
        external=c.diff(lambda Z:D(endpoints[0],Z)-D(endpoints[-1],Z),z,k)
        differences.append(abs(signed-external));chain.append(c.nstr(abs(signed-external),20))
    uy=P*c.diff(lambda y:V(y,z),0)
    slow=P*rate(z)*B(0,z,phase)/N
    fast=-P*b(z)*c.sin(phase)
    differences.append(abs(uy-slow-fast))
    wrong_J=c.diff(lambda x:rmU(x,z),1,2)
    wrong_fast=slow+fast/N
    correct_external=D(endpoints[0],z)-D(endpoints[-1],z)
    controls=dict(
        variable_J_derivative_required=abs(wrong_J-physical2)>c.mpf('.01'),
        no_extra_N_on_fast_term=abs(wrong_fast-uy)>c.mpf('.01'),
        incoming_memory_not_zero=abs(m0(z))>c.mpf('.01'),
        lower_minus_upper_sign=abs(-correct_external-correct_external)>c.mpf('.01'),
        physical_R_required=abs(correct_external-(correct_external/R0))>c.mpf('.01'),
        leading_subtraction_retained=abs(D(endpoints[0],z)-R0*c.exp(endpoints[0])*Qy(endpoints[0],z))>c.mpf('.01'),
        phase_not_reset=abs(B(0,z,phase)-B(0,z,0))>c.mpf('.01'),
        ordinary_mixed_Z_derivative_nonzero=abs(c.diff(lambda Z:rate(Z)*b(Z),z))>c.mpf('.01'),
        external_inlet_retained=abs(D(endpoints[0],z))>c.mpf('.01'),
        external_outlet_retained=abs(D(endpoints[-1],z))>c.mpf('.01'))
    maximum=max(differences)
    return dict(passed=maximum<c.mpf('1e-40') and all(controls.values()),
        synthetic_parameters_only=True,original_source_values_evaluated=False,finite_N=N,
        actual_Poisson_inverse_evaluated=False,variable_J='1/x',physical_y='log(x)',
        maximum_supporting_difference=c.nstr(maximum,25),chain_C0_Z1_differences=chain,
        controls=controls,control_count=len(controls),physical_second_derivative_keeps_J_x=True,
        nonzero_incoming_mean=True,ordinary_Z1_retained=True,not_an_NS_residual=True)


if __name__=='__main__':
    import json
    result=check_fixture();print(json.dumps(result,indent=2));raise SystemExit(0 if result['passed'] else 1)
