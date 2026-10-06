"""Actual five-moment patch in current-radius raw paper stress units."""
import math
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_current_O2_background_tensor import SourceAST
from lei_ren_part1_paper_compliant_current_O3_transition_background_tensor import copy_jet,IntervalTaylor,endpoints
from lei_ren_part1_paper_compliant_actual_patch_mixed_C4 import (
    power_derivatives,product_derivatives,log_radial_derivatives,stirling_second)
from lei_ren_part1_paper_compliant_collar_stress_C3 import shifted_rows
from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z
from lei_ren_part1_paper_compliant_inner_bridge_profiles import square


def raw_patch_rows(c,packet,invP2):
    """Convert full ordinary x/axial5 rows before physical differentiation."""
    x=c.mpf(packet['x']);z=IntervalTaylor.variable(c,c.mpf(packet['Z']),5)
    am=(1+square(z)).reciprocal()*c.exp(c.mpf('-.6'))
    source=packet['actual_normalized_primitive_x_derivative_axial5']
    raw={key:[copy_jet(c,row) for row in source[name]] for key,name in
        (('m','mean'),('h','theta'),('k','mixed'),('e','energy'),('p','pressure'))}
    H=[copy_jet(c,row) for row in packet['actual_H_x_derivative_axial5']]
    V=[copy_jet(c,row) for row in packet['actual_V_x_derivative_axial5']]
    G=[raw['m'][0]*x]+V[:4]
    inverse=power_derivatives(c,x,-1);angular=power_derivatives(c,x,'-1.5')
    xrows=[x,c.mpf(1),c.mpf(0),c.mpf(0),c.mpf(0)]
    energy=[(z*G[j]*8-square(z)*(16*xrows[j]))*invP2+square(am)*raw['e'][j] for j in range(5)]
    histories=dict(m=raw['m'],
        h=[am*product_derivatives(raw['h'],angular,j) for j in range(5)],
        k=[am*product_derivatives(raw['k'],angular,j) for j in range(5)],
        e=[product_derivatives(energy,inverse,j) for j in range(5)],
        p=[square(am)*row for row in raw['p']])
    p0=IntervalTaylor(c,[c.mpf(endpoints(value)) for value in
        packet['actual_inherited_patch_packet']['original_P0_axial5']])
    P=[p0+histories['p'][0]]+histories['p'][1:]
    Q=[copy_jet(c,row) for row in packet['actual_Q_x_derivative_axial4']]
    convert=lambda rows:log_radial_derivatives(c,x,rows)
    return dict(histories={key:convert(rows) for key,rows in histories.items()},
        absolute_pressure=convert(P),velocity=dict(theta=convert([am*row for row in H]),
            axial=convert(V),radial=shifted_rows(convert(Q),c.mpf('.5'),4)))


def actual_patch_normalization_theorem():
    """Replay the original patch program and this adapter in defining jets."""
    x=s.Symbol('x',positive=True);z=Z;delta=s.Symbol('delta',real=True)
    Ps=s.Symbol('Pstar',positive=True)
    c=SimpleNamespace(mpf=lambda v:(v[0] if isinstance(v,tuple) else
        s.Rational(str(v)) if isinstance(v,(str,int,float)) else v),exp=s.exp)
    J=lambda expr,order=5:FunctionJet.function(c,expr,order)
    am=J(s.exp(-s.Rational(3,5))/(1+z*z));zero=J(0)
    H=[J(s.Function('actual_H_'+str(j))(z)) for j in range(5)]
    g=[J(s.Function('actual_g_'+str(j))(z)) for j in range(5)]
    V=[J(4*z)+g[0]]+g[1:]
    initial={key:J(s.Function('actual_'+key)(z)) for key in ('mass','theta','mixed','energy','pressure')}
    p0=J(s.Function('same_actual_P0')(z))
    env=dict(math=math,IntervalTaylor=FunctionJet,square=lambda v:v*v,
        derivative=lambda v:J(s.diff(v.expr,z),v.order-1),copy_jet=lambda ctx,value:value,
        endpoints=lambda v:(v,v),shifted_rows=shifted_rows)
    asts=SourceAST()
    for name in ('power_derivatives','product_derivatives','stirling_second','log_radial_derivatives'):
        asts.replay('actual_patch_mixed_C4',name,env)
    original=asts.replay('actual_patch_mixed_C4','patch_mixed',env)
    adapter=asts.replay('current_patch_stress_operator','raw_patch_rows',env)
    packet=original(c,x,z,delta,am,J(1/(Ps**2*am.expr**2)),1/Ps**2,H,V,g,initial,p0)
    packet.update(x=x,Z=z,actual_H_x_derivative_axial5=H,actual_V_x_derivative_axial5=V,
        actual_inherited_patch_packet=dict(original_P0_axial5=p0.coefficients))
    result=adapter(c,packet,1/Ps**2);checks={}
    velocity=result['velocity'];raw=result['histories'];P=result['absolute_pressure']
    for label,rows,factor in (('Utheta_over_Pstar',velocity['theta'],1),('Uz',velocity['axial'],1),
            ('Ur_over_sqrt_Rm_over_2',velocity['radial'],s.sqrt(x)),('P_over_Pstar2',P,1)):
        for j in range(5):
            for n in range(5-j):
                left=packet['physical_velocity_pressure_y_Z_mixed4'][label]['y%d_Z%d'%(j,n)]
                right=rows[j][n]*math.factorial(n)*factor
                if s.cancel(left-right)!=0:raise ArithmeticError('Actual patch velocity/pressure units differ: '+label+str((j,n)))
                checks[label+'_y%d_Z%d'%(j,n)]=True
    mapping=(('Mz_over_Rm','m',1),('Mtheta_over_sqrt2_Rm_1p5_Pstar','h',s.Rational(3,2)),
        ('Mtheta_z_over_sqrt2_Rm_1p5_Pstar','k',s.Rational(3,2)),('Mztheta_over_Rm_Pstar2','e',1),('Mp_over_Pstar2','p',0))
    for label,key,power in mapping:
        rows=shifted_rows(raw[key],power,4)
        grid=packet['physical_five_primitive_x_Z_mixed4'][label]
        for j in range(5):
            for n in range(5-j):
                left=(grid['x0_Z%d'%n] if j==0 else sum(stirling_second(j,k)*x**k*grid['x%d_Z%d'%(k,n)] for k in range(1,j+1)))
                right=x**power*rows[j][n]*math.factorial(n)
                if s.cancel(left-right)!=0:raise ArithmeticError('Actual patch primitive current-radius units differ: '+label+str((j,n)))
                checks[label+'_y%d_Z%d'%(j,n)]=True
    u,v=velocity['theta'][0].expr,velocity['axial'][0].expr
    rhs=dict(m=v-raw['m'][0].expr,h=u-s.Rational(3,2)*raw['h'][0].expr,
        k=u*v-s.Rational(3,2)*raw['k'][0].expr,e=v*v/Ps**2-u*u/2-raw['e'][0].expr,p=u*u/2)
    for key,value in rhs.items():
        if s.cancel(raw[key][1].expr-value)!=0:raise ArithmeticError('Actual patch five source ODE differs: '+key)
        checks['raw_five_history_ODE_'+key]=True
    if s.cancel(P[1].expr-u*u/2)!=0:raise ArithmeticError('Absolute pressure source ODE differs')
    checks['absolute_pressure_y_ODE']=True
    return dict(identities=checks,original_patch_mixed_and_actual_unit_adapter_AST_replayed=True,
        arbitrary_current_H_g_five_histories_and_shared_P0_defining_jets_used=True,
        current_radius='R=Rm*x; Rm=Rref*exp(-6); D_y=x*D_x',
        radial_velocity='sqrt(Rm/2)*sqrt(x)*Q=sqrt(R/2)*Q; radial full rows=(D_y+1/2)^j Q',
        full_energy_cross_terms_and_Pstar_inverse_retained=True,
        absolute_P0_added_once_before_physical_differentiation=True,
        same_source_functions_not_interval_overlap=True,input_hashes=asts.hashes,passed=True)
