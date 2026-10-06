"""Signed original paper columns/remainder of the actual repaired source.

Raw Uz and m/k carry additional Pstar sectors. Only exact radius powers
are shifted; all local/history derivatives and nonlinear products remain.
"""
import math
from types import SimpleNamespace
import sympy as s
from lei_ren_part1_paper_compliant_current_O3_independent_repair_operator import SourceAST
from lei_ren_part1_paper_compliant_current_pre_pulse_stress_operator import (
    axial_derivative,shifted_rows,product_rows,axial_operator_rows,raw_pre_formula_theorem)

def sectors(parts,expected):
    powers=[part['Pstar_power'] for part in parts]
    if powers!=list(expected):raise ValueError('Actual modified Pstar source sectors changed')
    return {part['Pstar_power']:part['ordinary_logR_rows'] for part in parts}

def modified_pre_stress_rows(c,delta,z,view):
    velocity=view['modified_cylindrical_velocity_source_log_sectors']
    history=view['modified_five_histories_in_original_normalized_units_Pstar_sectors']
    u=sectors(velocity['theta'],(1,))[1];V=sectors(velocity['axial'],(1,))[1]
    m=sectors(history['m'],(0,1));k=sectors(history['k'],(0,1))
    h=sectors(history['h'],(0,))[0];e=sectors(history['e'],(0,))[0]
    P=view['modified_absolute_pressure_over_Pstar2_ordinary_logR_rows']
    d=1-z*z;L=1-z*z*delta;b=(1-delta)/2;kk=1-delta/2
    transport={power:[z*rows[j]*(1-delta)+d*axial_derivative(rows[j]) for j in range(5)]
      for power,rows in m.items()}
    theta={};axial={}
    def add(parts,name,rate,power,shape):
        parts[name+'_P'+str(power)]=dict(mode=(rate,power,0,0),shape=shape,
          full_derivative_rows=shifted_rows(shape,c.mpf(str(rate))))
    add(theta,'local_transport',.5,1,[-u[j]/L for j in range(4)])
    add(theta,'retained_angular_moment',.5,1,[(h[j]*kk-z*b*axial_derivative(h[j]))/L for j in range(4)])
    for power,rows in k.items():
        add(theta,'retained_mixed_moment',.5,1+power,[(z*(2*delta-1)*rows[j]-d*axial_derivative(rows[j]))/L for j in range(4)])
    for power,rows in transport.items():
        product=product_rows(u,rows)
        add(theta,'meridional_transport',.5,1+power,[product[j]/L for j in range(4)])
    add(theta,'variable_radial_shear',-.5,1,[2*u[j+1]-u[j] for j in range(4)])
    add(axial,'local_axial_transport',.5,1,[-V[j]/L for j in range(4)])
    for power,rows in transport.items():
        product=product_rows(V,rows)
        add(axial,'nonlinear_meridional_transport',.5,1+power,[product[j]/L for j in range(4)])
    for power,rows in m.items():
        add(axial,'retained_linear_axial_moment',.5,power,[(rows[j]-z*axial_derivative(rows[j]))*(1-delta)/(2*L) for j in range(4)])
    add(axial,'retained_full_energy',.5,2,[(z*(2*delta)*e[j]-d*axial_derivative(e[j]))/L for j in range(4)])
    add(axial,'actual_absolute_pressure',.5,2,[(z*(2*(1+delta))*P[j]-d*axial_derivative(P[j]))/L for j in range(4)])
    add(axial,'axial_radial_shear',-.5,1,[2*V[j+1] for j in range(4)])
    return dict(theta=theta,axial=axial)

def modified_pre_remainder_sectors(c,delta,z,view):
    velocity=view['modified_cylindrical_velocity_source_log_sectors']
    Ur=sectors(velocity['radial'],(0,1));Ut=sectors(velocity['theta'],(1,));Uz=sectors(velocity['axial'],(1,))
    L=1-z*z*delta;out=dict(radial={},theta={},axial={})
    def add(label,name,rate,power,beta,half,rows):
        out[label][name+'_P'+str(power)]=dict(mode=(rate,power,0,0),normalization_half=half,beta=beta,rows=rows)
    dzUr={power:axial_operator_rows(c,rows,-1,z,delta,count=2,order=1) for power,rows in Ur.items()}
    for power,rows in Ur.items():
        add('radial','time',.5,power,-3,True,[(rows[j]/2+z*axial_derivative(rows[j])*((1-delta)/2)+rows[j+1])/L for j in range(3)])
        add('radial','radial_viscosity',-.5,power,-3,True,shifted_rows([-(rows[j+2]-rows[j]/4)*2 for j in range(3)],-1,2))
        add('radial','axial_viscosity',.5,power,-3+2*delta,True,[-value for value in axial_operator_rows(c,rows,-1,z,delta)])
    nonlinear={power:[Ur[0][0]*0 for j in range(3)] for power in (0,1,2)}
    for pa,left in Ur.items():
        for pb,right in Ur.items():
            radial=shifted_rows(product_rows(left[:3],right[1:4]),-c.mpf('.5'),2)
            nonlinear[pa+pb]=[a+b for a,b in zip(nonlinear[pa+pb],radial)]
    for pa,left in Uz.items():
        for pb,right in dzUr.items():
            axial=product_rows(left[:3],right)
            nonlinear[pa+pb]=[a+b for a,b in zip(nonlinear[pa+pb],axial)]
    for power,rows in nonlinear.items():add('radial','nonlinear_transport',.5,power,-3,True,rows)
    for label,parts in (('theta',Ut),('axial',Uz)):
        for power,rows in parts.items():
            add(label,'axial_viscosity',0,power,-3+delta,False,[-value for value in axial_operator_rows(c,rows,-1-delta,z,delta)])
    return out

def exact_modified_pre_stress_theorem():
    """Replay actual sector operators against unchanged full raw programs."""
    paper=raw_pre_formula_theorem()
    if not paper['passed'] or not paper['original_full_paper_stress_AST_replayed']:
        raise ValueError('Original full variable-source paper stress identity is required')
    asts=SourceAST();checks={};P,delta=s.symbols('Pstar delta',real=True);z=s.Symbol('Z',real=True)
    c=SimpleNamespace(mpf=lambda value:s.Rational(str(value)))
    def zero(name,left,right):
        if s.cancel(left-right)!=0:raise ArithmeticError('Modified paper sector identity: '+name)
        checks[name]=True
    asts.expression('current_modified_pre_stress','__init__','self.source',wanted='source if source is not None else CurrentO3RepairedHistories()')
    asts.expression('current_modified_pre_stress','stress','view',wanted='self.source.history(region,Z,coordinate)')
    asts.expression('current_modified_pre_stress','stress','columns',wanted='modified_pre_stress_rows(c,self.delta,z,view)')
    asts.expression('current_modified_pre_stress','stress','remainder',wanted='modified_pre_remainder_sectors(c,self.delta,z,view)')
    checks['actual_checked_repaired_field_is_same_stress_and_remainder_source']=True
    row=lambda name:[s.Function(name+str(j))(z) for j in range(5)]
    u,V,m0,m1,h,k0,k1,e,p=(row(name) for name in ('u','Vhat','m0','m1','h','k0','k1','e','p'))
    part=lambda power,rows:dict(Pstar_power=power,ordinary_logR_rows=rows)
    view=dict(modified_cylindrical_velocity_source_log_sectors=dict(theta=[part(1,u)],axial=[part(1,V)]),
      modified_five_histories_in_original_normalized_units_Pstar_sectors=dict(m=[part(0,m0),part(1,m1)],
        h=[part(0,h)],k=[part(0,k0),part(1,k1)],e=[part(0,e)],p=[part(0,p)]),
      modified_absolute_pressure_over_Pstar2_ordinary_logR_rows=p)
    env=dict(math=math,axial_derivative=lambda value:s.diff(value,z),product_rows=product_rows,shifted_rows=shifted_rows)
    asts.replay('current_modified_pre_stress_operator','sectors',env)
    left=asts.replay('current_modified_pre_stress_operator','modified_pre_stress_rows',env)(c,delta,z,view)
    raw=dict(m=[a+P*b for a,b in zip(m0,m1)],h=h,k=[a+P*b for a,b in zip(k0,k1)],e=e,p=p)
    right=asts.replay('current_pre_pulse_stress_operator','raw_pre_stress_rows',env)(c,delta,z,u,[P*v for v in V],raw,p)
    for label,parts in right.items():
        for name,original in parts.items():
            chosen=[part for key,part in left[label].items() if key.rsplit('_P',1)[0]==name]
            if not chosen:raise ValueError('Missing complete original stress term: '+name)
            for j in range(4):
                value=sum(P**part['mode'][1]*part['full_derivative_rows'][j] for part in chosen)
                zero('actual_signed_stress_'+label+'_'+name+'_y'+str(j),value,P**original['mode'][1]*original['full_derivative_rows'][j])
    # Symbolic ordinary source projections preserve the very same axial
    # derivative operators, with arbitrary new and incoming radial rows.
    from lei_ren_part1_paper_compliant_actual_Rp_source_join import FunctionJet,Z
    from lei_ren_part1_paper_compliant_pulse_physical_bounds import physical_operators
    class Projection(FunctionJet):
        def __pow__(self,n):return self.function(self.ctx,self.expr**n,self.order)
    jet=lambda value,order=5:Projection.function(c,value,order)
    derivative=lambda value,n=1:jet(s.diff(value.expr,Z,n),value.order-n)
    velocity=dict(radial={0:[jet(v) for v in row('radial0')],1:[jet(v) for v in row('radial1')]},
      theta={1:[jet(v) for v in u]},axial={1:[jet(v) for v in V]})
    vview=dict(modified_cylindrical_velocity_source_log_sectors={key:[part(power,rows) for power,rows in parts.items()] for key,parts in velocity.items()})
    renv=dict(math=math,axial_derivative=derivative,axial_n=derivative,physical_operators=physical_operators,
      product_rows=product_rows,shifted_rows=shifted_rows)
    asts.replay('pulse_end_physical_C2','axial_operator_rows',renv)
    asts.replay('current_modified_pre_stress_operator','sectors',renv)
    rleft=asts.replay('current_modified_pre_stress_operator','modified_pre_remainder_sectors',renv)(c,delta,jet(Z),vview)
    total_velocity={key:[sum(P**power*rows[j] for power,rows in parts.items()) for j in range(5)] for key,parts in velocity.items()}
    # raw_pre expects normalized theta and raw physical axial/radial.
    total_velocity['theta']=velocity['theta'][1]
    rright=asts.replay('current_pre_pulse_stress_operator','raw_pre_remainder_sectors',renv)(c,delta,jet(Z),total_velocity)
    original_remainder=asts.replay('pulse_end_physical_C2','pulse_remainder_sectors',renv)(c,delta,None,jet(Z),total_velocity)
    for label,parts in rright.items():
        for name,original in parts.items():
            chosen=[part for key,part in rleft[label].items() if key.rsplit('_P',1)[0]==name]
            for part in chosen:
                if (part['mode'][0],part['beta'],part['normalization_half'])!=(original['mode'][0],original['beta'],original['normalization_half']):
                    raise ValueError('Modified remainder physical radius/beta/half mode differs')
            for j in range(3):
                value=sum(P**part['mode'][1]*part['rows'][j].expr for part in chosen)
                zero('actual_signed_remainder_'+label+'_'+name+'_y'+str(j),value,P**original['mode'][1]*original['rows'][j].expr)
                zero('actual_original_full_remainder_coefficient_'+label+'_'+name+'_y'+str(j),
                  original['rows'][j].expr,original_remainder[label][name]['rows'][j].expr)
    return dict(identities=checks,passed=True,input_hashes={**paper['input_hashes'],**asts.hashes},
      consumed_arbitrary_variable_original_full_paper_stress_theorem=paper,
      theta_stress_sector_count=len(left['theta']),axial_stress_sector_count=len(left['axial']),
      remainder_sector_counts={key:len(parts) for key,parts in rleft.items()},
      full_actual_raw_pre_paper_programs_replayed=True,
      physical_Uz_is_Pstar_times_Vhat_not_Vhat=True,
      all_radial_Pstar0_Pstar1_square_and_cross_terms_retained=True,
      ordinary_stress3_remainder2_and_axial_source_derivative_operators_preserved=True,
      completed_diagonal_tensor_energy_and_cone_not_admitted=True)
