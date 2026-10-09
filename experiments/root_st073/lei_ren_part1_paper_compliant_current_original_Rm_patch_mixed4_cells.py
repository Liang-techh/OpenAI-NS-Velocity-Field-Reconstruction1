"""Actual leading Rm patch: full mixed4, raw current units and radial cells.

Accepted native0,.5 coefficient/source jets feed the original primitive
RHSs. Closed rational x cells enclose the same functions across flat bump
edges. Whole-Z and finite-N Rc/source density closure remain separate.
"""
from fractions import Fraction
import gzip
import json
import math
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rm_defect_patch_inverse as previous
import lei_ren_part1_paper_compliant_actual_patch_mixed_C4 as original
from lei_ren_part1_paper_compliant_flat_pulse_derivatives import beta_jets
from lei_ren_part1_paper_compliant_core_physical_field import intersection

fields,base,ep=previous.fields,previous.base,previous.ep
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_Rm_patch_mixed4_cells.json.gz'
RECEIPT=PREFIX+'current_original_Rm_patch_mixed4_cells_check.json'
GATE='original_actual_two_frame_Rm_patch_mixed4_and_radial_cells_installed'
FLAT=PREFIX+'flat_pulse_derivatives_check.json'


def source_bindings():
    from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
    bound=assignment_source_bindings('actual_patch_mixed_C4','patch_mixed',{
        'G':"[initial['mass']*x]+V[:4]",
        'mass':"[sum((G[j]*inverses[k-j]*math.comb(k,j) for j in range(k+1)),zero) for k in range(5)]",
        'HV':"[product_derivatives(H,V,k) for k in range(4)]",
        'H2':"[product_derivatives(H,H,k) for k in range(4)]",
        'g2':"[product_derivatives(g,g,k) for k in range(4)]",
        'Q':"[(2*z*V[k]-(z*mass[k])*(1-delta)-d*derivative(mass[k]))/L for k in range(5)]",
        'P':"[p0+square(am)*pressure[0]]+[square(am)*row for row in pressure[1:]]",
        'physical_primitives':"dict(Mz_over_Rm=G,Mtheta_over_sqrt2_Rm_1p5_Pstar=[am*row for row in theta],Mtheta_z_over_sqrt2_Rm_1p5_Pstar=[am*row for row in mixed],Mztheta_over_Rm_Pstar2=[((z*G[k])*8-square(z)*(16*xder[k]))*invP2+square(am)*energy[k] for k in range(5)],Mp_over_Pstar2=[square(am)*row for row in pressure])"})
    # The append statements fix the k+1 derivative index of each primitive.
    import ast
    tree=ast.parse(Path(original.__file__).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='patch_mixed')
    appends=("theta.append(sum((H[j]*(sqrtx[k-j]*math.comb(k,j)) for j in range(k+1)),zero))",
        "mixed.append(sum((HV[j]*(sqrtx[k-j]*math.comb(k,j)) for j in range(k+1)),zero))",
        "energy.append(g2[k]*invAm2-H2[k]/2)",
        "pressure.append(sum((H2[j]*(inverses[k-j]*math.comb(k,j)/2) for j in range(k+1)),zero))")
    for text in appends:
        want=ast.dump(ast.parse(text,mode='eval').body)
        if sum(ast.dump(n)==want for n in ast.walk(fn) if isinstance(n,ast.Call))!=1:
            raise ValueError('Original patch primitive derivative RHS changed')
    return dict(passed=True,original_assignments=bound,original_primitive_append_RHS_count=len(appends),
        primitive_derivative_k_plus1_uses_source_derivative_k=True,
        log_coordinate='Dy=x*Dx; Stirling second-kind weights; DR^k=Rm^-k Dx^k',
        raw_current_units='h,k=am*x^-1.5*theta,mixed; e=physical_Mztheta/(R*Pstar2); p=Mp/Pstar2; radial=(Dy+.5)^k Q')


def product_derivatives(f,A,B,k):
    return f.add(*[f.scale(f.multiply(A[j],B[k-j]),math.comb(k,j)) for j in range(k+1)])


def scalar_product(f,rows,scalar_derivatives,k):
    return f.add(*[f.scale(rows[j],scalar_derivatives[k-j]*math.comb(k,j)) for j in range(k+1)])


def log_rows(f,x,rows):
    return [rows[0]]+[[sum((rows[j][n]*(original.stirling_second(k,j)*x**j)
        for j in range(1,k+1)),f.scalar(0)) for n in range(len(rows[0]))] for k in range(1,len(rows))]


def gamma_rows(op,x,center,exact_edge=False):
    c=op.c;r=c.mpf(1)/40;center=c.mpf(center.numerator)/center.denominator
    if exact_edge:return [c.mpf(0)]*5
    beta=beta_jets(c,(x-center)/r)
    rows=[beta[k]*(math.factorial(k)/(r**(k+1)*op.bump.normalization)) for k in range(5)]
    old=op.bump.beta(x-center)
    rows[0]=intersection(c,rows[0],old[0]);rows[1]=intersection(c,rows[1],old[1])
    return rows


def partial_initial(op,x,terminal):
    """Closed-cell version of the SAME accepted leading partial primitives."""
    c=op.c;f=op.flow;h=op.controls;zero=[f.scalar(0)]*6
    centers=(Fraction(5,4),Fraction(3,2),Fraction(7,4));weights={}
    def w(k,p,m=1):
        key=(k,p,m)
        if key not in weights:weights[key]=op.bump.partial_weight(x,centers[k],p,m,64)
        return weights[key]
    changes=[f.add(*[f.scale(h[i],w(k,'0')) for i,k in enumerate((0,2))]),
        f.add(*[f.add(f.scale(h[i],w(k,'.6')),f.scale(f.multiply(h[i],h[k+2]),w(k,'.5',2))) for i,k in enumerate((0,2))]),
        f.add(*[f.scale(h[k+2],w(k,'.5')) for k in range(3)]),
        f.add(*[f.scale(f.multiply(f.multiply(h[i],h[i]),op.invAm2),w(k,'0',2)) for i,k in enumerate((0,2))],
            *[f.scale(h[k+2],-w(k,'.1')) for k in range(3)],
            *[f.scale(f.multiply(h[k+2],h[k+2]),-w(k,'0',2)/2) for k in range(3)]),
        f.add(*[f.add(f.scale(h[k+2],w(k,'-.9')),f.scale(f.multiply(h[k+2],h[k+2]),w(k,'-1',2)/2)) for k in range(3)])]
    diagnostic=[f.add(a,b) for a,b in zip(op.defects,changes)]
    if terminal:
        if any(not ep(v.coefficient)[0]<=0<=ep(v.coefficient)[1] for row in diagnostic for v in row):
            raise ArithmeticError('Complete-support cell implicit-map diagnostic excludes zero')
        defects=[zero]*5
    else:defects=diagnostic
    unit=lambda v:[f.scalar(v)]+zero[1:]
    theta=f.add(defects[2],unit(c.mpf(5)/8*x**c.mpf('1.6')))
    initial=dict(mass=f.add(f.scale(op.zrows,4),f.scale(defects[0],1/x)),theta=theta,
        mixed=f.add(f.scale(f.multiply(op.zrows,theta),4),defects[1]),
        energy=f.add(defects[3],unit(-c.mpf(5)/12*x**c.mpf('1.2'))),
        pressure=f.add(defects[4],unit(c.mpf(5)/2*x**c.mpf('.2'))))
    return initial,dict(original_actual_Rm_defects=op.defects,partial_primitive_changes=changes,
        partial_weights={str(k):v for k,v in weights.items()},unrefined_defect_diagnostic=diagnostic,
        exact_terminal_identity_of_same_leading_map=terminal,inlet_never_reset=True)


def mixed_functions(op,x,initial,gamma,*,terminal=False):
    """Original radial primitive RHSs lifted to the canonical formal algebra."""
    f=op.flow;c=op.c;zero=[f.scalar(0)]*6;unit=lambda v:[f.scalar(v)]+zero[1:]
    powers=original.power_derivatives(c,x,'.1');H=[];g=[]
    for k in range(5):
        H.append(f.add(unit(powers[k]),*[f.scale(op.controls[i+2],gamma[i][k]) for i in range(3)]))
        g.append(f.add(*[f.scale(op.controls[i],gamma[j][k]) for i,j in enumerate((0,2))]))
    if ep(H[0][0].coefficient)[0]<=0:raise ArithmeticError('Positive source swirl lost in patch cell')
    V=[f.add(f.scale(op.zrows,4),g[0])]+g[1:]
    inverse=original.power_derivatives(c,x,-1);sqrtx=original.power_derivatives(c,x,'.5')
    G=[f.scale(initial['mass'],x)]+V[:4]
    mass=[scalar_product(f,G,inverse,k) for k in range(5)]
    if terminal:
        if any(ep(v)!=(0,0) for rows in gamma for v in rows):
            raise ValueError('Exact complete-support branch requires vanishing beta')
        # Functional full-weight identities of the same unique leading map:
        # G=4Z*x, m=4Z, V=4Z on the entire exact terminal window.
        G=[f.scale(op.zrows,4*x),f.scale(op.zrows,4)]+[zero]*3
        mass=[f.scale(op.zrows,4)]+[zero]*4
    theta=[initial['theta']];mixed=[initial['mixed']];energy=[initial['energy']];pressure=[initial['pressure']]
    HV=[product_derivatives(f,H,V,k) for k in range(4)];H2=[product_derivatives(f,H,H,k) for k in range(4)];g2=[product_derivatives(f,g,g,k) for k in range(4)]
    for k in range(4):
        theta.append(scalar_product(f,H,sqrtx,k));mixed.append(scalar_product(f,HV,sqrtx,k))
        energy.append(f.add(f.multiply(g2[k],op.invAm2),f.scale(H2[k],-c.mpf('.5'))))
        pressure.append(f.scale(scalar_product(f,H2,inverse,k),c.mpf('.5')))
    Q=[previous.previous.previous.radial_Q(f,op.reference.Z,op.reference.delta,V[k],mass[k]) for k in range(5)]
    if terminal:
        z=op.reference.z;delta=op.reference.delta
        Q=[f.jet(((z*z*(2+delta)-1)/(1-z*z*delta))*4)[:5]]+[zero[:5]]*4
    # Q has Z0..4 only. Pad internally for products, never export Z5.
    Qpad=[row+[f.scalar(0)] for row in Q]
    Ur=[scalar_product(f,Qpad,sqrtx,k)[:5] for k in range(5)]
    am2=f.multiply(op.amrows,op.amrows);P2=op.P2;Rm=op.Rm_factor
    axis=f.scale(op.P0,P2);inc=[f.scale(f.multiply(am2,row),P2) for row in pressure]
    P=[f.add(axis,inc[0])]+inc[1:]
    xder=[x,c.mpf(1)]+[c.mpf(0)]*3
    centered=[f.scale(f.multiply(am2,row),Rm*P2) for row in energy]
    angular=original.power_derivatives(c,x,'-1.5')
    raw_energy=[f.add(f.scale(f.add(f.scale(f.multiply(op.zrows,G[k]),8),
        f.scale(f.multiply(op.zrows,op.zrows),-16*xder[k])),f.factor((0,-1,0,0,0))),f.multiply(am2,energy[k])) for k in range(5)]
    # Rm is one common positive source constant. Its gigantic directed log
    # offset must not be independently subtracted from itself in sums or
    # coordinate conversion. Convert normalized source rows first and attach
    # Rm^(a-k) exactly once, after the final ordinary derivative grid.
    normalized=dict(Utheta=[f.multiply(op.amrows,row) for row in H],Uz=V,Ur=Ur,
        P=[f.add(op.P0,f.multiply(am2,pressure[0]))]+[f.multiply(am2,row) for row in pressure[1:]],
        Mz=G,Mtheta=[f.multiply(op.amrows,row) for row in theta],
        Mtheta_z=[f.multiply(op.amrows,row) for row in mixed],Mztheta=raw_energy,
        Mp=[f.multiply(am2,row) for row in pressure])
    radius_exponents=dict(Utheta=0,Uz=0,Ur=.5,P=0,Mz=1,Mtheta=1.5,Mtheta_z=1.5,Mztheta=1,Mp=0)
    amplitude=dict(Utheta=op.Pstar,Uz=f.scalar(1),Ur=f.scalar(1/c.sqrt(2)),P=P2,
        Mz=f.scalar(1),Mtheta=op.Pstar*c.sqrt(2),Mtheta_z=op.Pstar*c.sqrt(2),Mztheta=P2,Mp=P2)
    def physical_unit(name,k=0):return amplitude[name]*f.radial_power(Rm,radius_exponents[name]-k)
    physical_primitives={name:[f.scale(row,physical_unit(name)) for row in normalized[name]]
                         for name in ('Mz','Mtheta','Mtheta_z','Mztheta','Mp')}
    physical_velocity={name:[f.scale(row,physical_unit(name)) for row in normalized[name]]
                       for name in ('Utheta','Uz','Ur')}
    raw_hist_x=dict(m=mass,h=[f.multiply(op.amrows,scalar_product(f,theta,angular,k)) for k in range(5)],
        k=[f.multiply(op.amrows,scalar_product(f,mixed,angular,k)) for k in range(5)],
        e=[scalar_product(f,raw_energy,inverse,k) for k in range(5)],p=[f.multiply(am2,row) for row in pressure])
    raw_hist={name:log_rows(f,x,rows) for name,rows in raw_hist_x.items()}
    Qy=log_rows(f,x,Qpad)
    raw_radial=[f.add(*[f.scale(Qy[j],math.comb(k,j)*c.mpf('.5')**(k-j)) for j in range(k+1)])[:5] for k in range(5)]
    raw_pressure=[f.add(op.P0,raw_hist['p'][0])]+raw_hist['p'][1:]
    grid=lambda rows,letter:{letter+str(k)+'_Z'+str(n):row[n]*math.factorial(n) for k,row in enumerate(rows) for n in range(5-k)}
    common_x={name:grid(rows,'x') for name,rows in normalized.items()}
    common_y={name:grid(log_rows(f,x,rows),'y') for name,rows in normalized.items()}
    def physical_grid(name,letter):
        source=common_y[name] if letter=='y' else common_x[name]
        return {letter+key[1:]:row*physical_unit(name,int(key.split('_')[0][1:]) if letter=='R' else 0)
                for key,row in source.items()}
    return dict(actual_H_x_derivative_axial5=H,actual_V_x_derivative_axial5=V,actual_g_x_derivative_axial5=g,
        actual_normalized_primitive_x_derivative_axial5=dict(mean=mass,theta=theta,mixed=mixed,energy=energy,pressure=pressure),
        actual_Q_x_derivative_axial4=Q,
        physical_velocity_x_derivative_axial_coefficients=physical_velocity,
        physical_pressure_x_derivative_axial5=P,physical_five_primitive_x_derivative_axial5=physical_primitives,
        physical_velocity_x_Z_mixed4={name:physical_grid(name,'x') for name in physical_velocity},
        physical_velocity_y_Z_mixed4={name:physical_grid(name,'y') for name in physical_velocity},
        physical_velocity_R_Z_mixed4={name:physical_grid(name,'R') for name in physical_velocity},
        physical_pressure_x_Z_mixed4=physical_grid('P','x'),physical_pressure_y_Z_mixed4=physical_grid('P','y'),physical_pressure_R_Z_mixed4=physical_grid('P','R'),
        physical_five_primitive_x_Z_mixed4={name:physical_grid(name,'x') for name in physical_primitives},
        physical_five_primitive_y_Z_mixed4={name:physical_grid(name,'y') for name in physical_primitives},
        physical_five_primitive_R_Z_mixed4={name:physical_grid(name,'R') for name in physical_primitives},
        physical_common_unit_normalized_x_Z_mixed4=common_x,physical_common_unit_normalized_y_Z_mixed4=common_y,
        shared_positive_physical_unit_radius_exponents=radius_exponents,
        shared_positive_physical_unit_amplitudes=amplitude,
        common_R_derivative_grid_is_x_grid_with_unit_Rm_power_minus_k=True,
        common_pressure_increment_grid_is_Mp_with_exact_same_separate_P0=True,
        shared_Rm_offset_never_independently_subtracted_in_physical_grid=True,
        exact_terminal_mean_and_radial_Q_refinement_from_same_map=terminal,
        physical_centered_Mztheta_x_derivative_axial5=centered,
        original_P0_normalized_axial5=op.P0,physical_pressure_axis_axial5=axis,physical_pressure_increment_x_derivative_axial5=inc,
        raw_current_radius_y_derivative_axial_coefficients=dict(histories=raw_hist,absolute_pressure=raw_pressure,
            velocity=dict(theta=log_rows(f,x,[f.multiply(op.amrows,row) for row in H]),axial=log_rows(f,x,V),radial=raw_radial)),
        raw_radial_Z_order4_only_no_selected_Z5=True,physical_prefactors_differentiated_before_grid=True,
        original_P0_only_at_radial_order0=True,primitive_source_index_not_shifted=True)


class _MixedPatchOwner:
    def __init__(self,op):
        if op.P0 is not op.reference.P0 or op.P0 is not op.Rm['original_P0_normalized_axial5']:
            raise ValueError('Same original P0 object required on every point/cell path')
        self.op=op;self.flow=op.flow;self.c=op.c;self.cache={}

    def coordinate(self,coordinate):
        if coordinate=='Rh':return self.c.exp(1),True,dict(exact_x_source='exp(1)',point=True)
        q=previous.previous.fraction(coordinate)
        x=self.c.mpf(q.numerator)/q.denominator
        if q<1 or ep(x)[1]>ep(self.c.exp(1))[0]:raise ValueError('Exact patch x in[1,e] required')
        return x,q>=Fraction(71,40),dict(exact_x=[q.numerator,q.denominator],point=True)

    def evaluate(self,coordinate):
        x,terminal,geometry=self.coordinate(coordinate);key=('point',str(coordinate))
        if key in self.cache:return self.cache[key]
        if coordinate=='Rh':initial,memory=partial_initial(self.op,x,terminal)
        else:
            parent=self.op.evaluate(coordinate);p=parent['normalized_primitives']
            if parent['original_P0_normalized_axial5'] is not self.op.P0:
                raise ValueError('Same inherited point P0 object required')
            initial={key:p[name] for key,name in (('mass','mean'),('theta','theta'),('mixed','mixed'),('energy','centered_energy'),('pressure','pressure'))}
            memory=dict(actual_point_partial_primitive_object_memory=initial,
                original_actual_Rm_defects=self.op.defects,exact_terminal_identity_of_same_leading_map=terminal,
                original_P0_is_same_object=parent['original_P0_normalized_axial5'] is self.op.P0)
        return self.build(x,initial,memory,geometry,key,coordinate)

    def cell(self,left,right):
        lo=previous.previous.fraction(left);hi=None if right=='Rh' else previous.previous.fraction(right)
        if hi is not None and lo>=hi:raise ValueError('Strictly ordered exact patch cell required')
        a,_,_=self.coordinate(left);b,terminal,_=self.coordinate(right)
        x=self.c.mpf([ep(a)[0],ep(b)[1]]);terminal=lo>=Fraction(71,40)
        key=('cell',lo,'Rh' if hi is None else hi)
        if key in self.cache:return self.cache[key]
        initial,memory=partial_initial(self.op,x,terminal)
        geometry=dict(exact_left=[lo.numerator,lo.denominator],point=False)
        if hi is None:geometry['exact_right_source']='exp(1)'
        else:geometry['exact_right']=[hi.numerator,hi.denominator]
        return self.build(x,initial,memory,geometry,key,None)

    def build(self,x,initial,memory,geometry,key,coordinate):
        centers=(Fraction(5,4),Fraction(3,2),Fraction(7,4));q=None if coordinate in (None,'Rh') else Fraction(*coordinate)
        terminal=memory['exact_terminal_identity_of_same_leading_map']
        gamma=([[self.c.mpf(0)]*5 for _ in centers] if terminal else
               [gamma_rows(self.op,x,center,q is not None and abs(q-center)==Fraction(1,40)) for center in centers])
        value=mixed_functions(self.op,x,initial,gamma,terminal=terminal)
        value.update(geometry=geometry,actual_gamma_ordinary_x_derivatives=gamma,actual_partial_primitive_source_memory=memory,
            source_frames_conditional_on_same_accepted_Rm_inlet=True,closed_x_cell_not_a_sample_interpolation=not geometry['point'],
            original_flat_support_edges_and_Rm_Rh_leading_joins_retained=True,
            actual_leading_patch_mixed4_installed=True,actual_patch_finite_N_density_oracle_installed=False,
            whole_axis_function_provider_or_finite_N_Rc_patch_installed=False)
        self.cache[key]=value;return value


class OriginalRmPatchMixed4Cells:
    mode='genuine_original_actual_Rm_leading_patch_mixed4_and_closed_radial_cells'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalRmDefectPatchInverse(dps);self.c=self.upstream.c;self.family=self.upstream.family
        self.hashes=dict(self.upstream.hashes);self.owners={};bind=fields.previous.bind
        flat=json.loads((HERE/FLAT).read_bytes())
        if not all(flat.get(gate) for gate in ('all_passed','original_radial_shape_derivatives_C4_available','quantitative_flat_support_majorants_available','support_crossing_derivatives_available')):
            raise ValueError('Accepted original beta flat derivative and crossing envelopes required')
        if flat['actual_five_defect_family_sha256']!=self.family or flat['implicit_source_sha256']!=self.upstream.repair['implicit_source_sha256']:
            raise ValueError('Same original beta source family required')
        for path,digest in flat['input_hashes'].items():bind(self.hashes,path,digest)
        for name in (FLAT,Path(original.__file__).name,'lei_ren_part1_paper_compliant_flat_pulse_derivatives.py',
            'lei_ren_part1_paper_compliant_current_patch_stress_operator.py',Path(__file__).name):bind(self.hashes,name,sha(name))
        self.bindings=source_bindings()
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked.get('all_passed') or not checked.get(GATE) or checked['source_family']!=self.family:
                raise ValueError('Accepted actual patch mixed4/cell receipt required')
            for path,digest in checked['input_hashes'].items():bind(self.hashes,path,digest)
            bind(self.hashes,RECEIPT,sha(RECEIPT))

    def owner(self,label):
        if label not in self.owners:self.owners[label]=_MixedPatchOwner(self.upstream.owner(label))
        return self.owners[label]

    def evaluate(self,label,coordinate):
        with mp.workdps(self.c.dps+40):value=self.owner(label).evaluate(coordinate)
        return dict(mode=self.mode,source_family=self.family,source_frame=label,function_evaluation=fields.serialized(value),
            whole_axis_functions_installed=False,no_original_ancestor_producers_or_full_checkers_executed=True)

    def cell(self,label,left,right):
        with mp.workdps(self.c.dps+40):value=self.owner(label).cell(left,right)
        return dict(mode=self.mode,source_family=self.family,source_frame=label,function_evaluation=fields.serialized(value),
            whole_axis_functions_installed=False,no_original_ancestor_producers_or_full_checkers_executed=True)


def run():
    began=time.monotonic();owner=OriginalRmPatchMixed4Cells(require_checked=False)
    points=((1,1),(5,4),(51,40),(3,2),(7,4),(71,40),(2,1),'Rh')
    cells=(((6,5),(13,10)),((29,20),(31,20)),((17,10),(9,5)),((71,40),'Rh'))
    report=dict(**{GATE:True},source_family=owner.family,original_source_bindings=owner.bindings,
        frames={label:dict(points=[owner.evaluate(label,q) for q in points],cells=[owner.cell(label,a,b) for a,b in cells]) for label in ('0','.5')},
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual native0,.5 leading patch full physical x/y/R mixed4, raw current-radius history/velocity rows and closed radial source cells. Same incoming partial primitives and P0. Whole-Z, finite-N Rc/density integrator, all24/global N and real n-recursion remain open.')
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),indent=2)+'\n').encode(),mtime=0))
    print('Actual leading Rm patch mixed4, raw current-radius rows and closed radial cells generated',flush=True);return report


if __name__=='__main__':run()
