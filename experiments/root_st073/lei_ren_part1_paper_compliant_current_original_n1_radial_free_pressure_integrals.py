"""Exact physical leading-pressure integrals without internal radial jets.

Use the leading own mean identity m_y=U-m and q radial recovery. All
moments and boundary differences are in physical R units, at fixed Z.
Finite endpoint q_R terms remain; support-zero is not inferred here.
"""


def reduced_S_integrals(z,delta,moments,boundaries):
    """Return integral S dR and integral R*S dR as ordinary Z jets."""
    a=1-delta;d=1-z*z;L=1-delta*z*z
    I0=moments['unweighted'];I1=moments['weighted'];B=boundaries
    S0=(B['R_q']/L+a*z*I0['q_Z']/(2*L)+B['R_q_squared']/2
        +d*I0['U_q_Z']/L-2*z*z*B['R_U_squared']/(L*L)
        -2*delta*z*z*I0['U_squared']/(L*L)+2*z*d*I0['U_U_Z']/(L*L))
    S1=(B['R_squared_q']/L-I1['q']/L+a*z*I1['q_Z']/(2*L)
        +B['R_squared_q_squared']/2-I1['q_squared']/2+d*I1['U_q_Z']/L
        -2*z*z*B['R_squared_U_squared']/(L*L)
        +2*a*z*z*I1['U_squared']/(L*L)+2*z*d*I1['U_U_Z']/(L*L))
    return dict(S=S0,R_S=S1)


def fixed_pressure_integrals(z,delta,moments,boundaries):
    reduced=reduced_S_integrals(z,delta,moments,boundaries)
    return dict(N1p=-reduced['S']/2+boundaries['q_plus_R_q_R'],
        R_N1p=-reduced['R_S']/2+boundaries['R_squared_q_R'])


def original_identity_requirements():
    return dict(own_mean='m=Mz0/R; m_y=U-m; y=logR',
        radial_recovery='L*q=2Z*U-(1-delta)*Z*m-(1-Z^2)*m_Z; q=V0/R',
        fixed_forcing='N1p=-S/2+2*q_R+R*q_RR',
        original_physical_measure='I0(f)=integral f dR; I1(f)=integral R*f dR',
        boundary='B(f)=f(upper)-f(lower); all finite upper terms retained',
        value_inputs='U,U_Z,q,q_Z and true own mean recovery; no internal radial/phase derivative',
        Z1_inputs='U through Z2, own mean through Z3, q through Z2; endpoints fixed in Z',
        full_pressure_tail_support_or_five_moment_closure_implied=False)
