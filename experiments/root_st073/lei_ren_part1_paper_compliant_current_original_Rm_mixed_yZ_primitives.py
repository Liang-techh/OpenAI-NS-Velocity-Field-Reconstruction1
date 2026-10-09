"""Genuine actual Rm yZ roots and original all-signed-u mixed primitives.

Mixed derivatives come from ordinary raw source rows and exact quotient
calculus. Fourier absolute and L2 identities enclose original fixed-phi
inverse derivatives for all finite signed u. No O2 restriction, selected
inverse value, bound derivative or actual spatial mixed chain is used.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rm_weighted_averaging as previous

first=previous.previous;phase=first.phase
fields,base,ep=first.fields,first.base,first.ep
HERE,PREFIX,sha=first.HERE,first.PREFIX,first.sha
NAME=PREFIX+'current_original_Rm_mixed_yZ_primitives.json.gz'
RECEIPT=PREFIX+'current_original_Rm_mixed_yZ_primitives_check.json'
GATE='original_actual_Rm_genuine_yZ_roots_and_all_u_fixed_phi_mixed_primitives_installed'
C0,Y,Z,YZ=(0,0),(1,0),(0,1),(1,1)
ORDERS=(C0,Y,Z,YZ)


def serialized(value):
    if isinstance(value,dict):return {('y%d_Z%d'%key if isinstance(key,tuple) else key):serialized(row)
        for key,row in value.items()}
    if isinstance(value,(tuple,list)):return [serialized(row) for row in value]
    return fields.serialized(value)


class Mixed:
    """Ordinary source 1 by 1 jets, with all four product terms."""
    def __init__(self,rows):
        if set(rows)!=set(ORDERS):raise ValueError('Genuine C0,y,Z,yZ source rows required')
        self.rows=rows;phase.previous.same_source(type('Algebra',(),dict(c=rows[C0].ctx,
            logs=rows[C0].scale.bases,ledger=rows[C0].ledger))(),list(rows.values()))
    def __getitem__(self,key):return self.rows[key]
    def coerce(self,value):
        if isinstance(value,Mixed):self[C0].coerce(value[C0]);return value
        value=self[C0].coerce(value)
        return Mixed({key:value if key==C0 else value.scalar(0) for key in ORDERS})
    def __neg__(self):return Mixed({key:-row for key,row in self.rows.items()})
    def __add__(self,other):
        other=self.coerce(other);return Mixed({key:self[key]+other[key] for key in ORDERS})
    __radd__=__add__
    def __sub__(self,other):return self+-self.coerce(other)
    def __rsub__(self,other):return self.coerce(other)+-self
    def __mul__(self,other):
        other=self.coerce(other);zero=self[C0].scalar(0)
        return Mixed({(j,k):sum((self[(i,l)]*other[(j-i,k-l)] for i in range(j+1)
            for l in range(k+1)),zero) for j,k in ORDERS})
    __rmul__=__mul__


def source_mixed_frame(op,patch,generic,quotients):
    """Differentiate first source identities with genuine raw yZ rows."""
    source,qr,record=first.source_first_frame(op,patch,generic,quotients)
    dual=first.lift_generic_y(op,patch);f=op.flow;c=op.c
    proofs=quotients['source_positive_theorems'];n=generic['actual_generic_source_numerators']
    def divide(value,denominator,proof):
        if proof['source_row'] is not denominator:raise ValueError('Same live source positive denominator required')
        return value.positive_divide(denominator,proof['source_log_lower'])
    E,EZ=n['E'][:2];C,CZ=n['C'][:2]
    Ey=source['roots']['E'][Y];EyZ=dual['common_velocity_E_axial5'][1].dy
    Cy,CyZ=(dual['actual_generic_source_numerators']['C'][k].dy for k in (0,1))
    By,ByZ=(dual['actual_generic_source_numerators']['B'][k].dy for k in (0,1))
    shear=quotients['shear_quotients_axial5'];a,aZ=shear['a'][:2];b,bZ=shear['b'][:2];t0,t0Z=shear['t0'][:2]
    radial=record['original_correlated_shear_y'];ay,by,t0y=(radial[k] for k in ('a_y','b_y','t0_y'))
    x=phase.previous.source_coordinate(c,patch['geometry']);gamma=patch['actual_gamma_ordinary_x_derivatives']
    FyZ=sum((op.controls[j+2][1]*(x*(c.mpf(9)/10*gamma[j][1]+x*gamma[j][2])) for j in range(3)),f.scalar(0))
    correlated=quotients['actual_correlated_shear_source'];FZ=correlated['angular_correction_numerator'][1]
    F,H,Hy=radial['F'],radial['H'],radial['H_y']
    HZ=patch['actual_H_x_derivative_axial5'][0][1];HyZ=patch['actual_H_x_derivative_axial5'][1][1]*x
    Hproof=correlated['positive_H_source_theorem'];J=divide(Hy,H,Hproof)
    JZ=divide(HyZ-J*HZ,H,Hproof)
    ayZ=divide(FyZ-FZ*J-F*JZ+ay*HZ*c.mpf('.5'),H,Hproof)*(-2)
    byZ=divide(ByZ-bZ*Ey-b*EyZ-by*EZ,E,proofs['E'])
    t0yZ=divide(-ByZ-t0Z*Cy-t0*CyZ-t0y*CZ,C,proofs['C'])
    g=divide(phase.first.current.square(b),a,proofs['a'])
    gZ=divide(b*bZ*2-g*aZ,a,proofs['a'])
    Delta_y=radial['Delta_y']
    Delta_yZ=ayZ+divide(bZ*by*2+b*byZ*2-gZ*ay-g*ayZ-(Delta_y-ay)*aZ,a,proofs['a'])
    loop=quotients['original_shear_q'];q2=loop['q_squared_axial_coefficients'];q=qr[C0]
    q2y=radial['q_squared_y'];qy=qr[Y];qZ=qr[Z]
    twice_a=a*2;twice_q=q*2
    q2yZ=(-Delta_yZ-q2[1]*ay*2-q2[0]*ayZ*2-q2y*aZ*2).positive_divide(
        twice_a,phase.previous.positive_source(f,twice_a,'2a')['source_log_lower'])
    qyZ=(q2yZ-qy*qZ*2).positive_divide(twice_q,phase.previous.positive_source(f,twice_q,'2q')['source_log_lower'])
    sectors=dual['full_signed_inertial_sectors_axial4']
    IyZ=sectors['axial_linear'][1].dy+sectors['axial_quadratic'][1].dy*op.Pstar
    bar=quotients['full_signed_inertial_quotients_before_shared_R_axial4']['p2']
    bary=record['original_generic_recovery_y_tangents']['p2_before_R_y']
    baryZ=divide(IyZ-bar[1]*Ey-bar[0]*EyZ-bary*EZ,E,proofs['E'])
    R=quotients['exact_same_shared_positive_radius_factor'];p2yZ=R*(baryZ+bar[1])
    derivatives=dict(a=ayZ,t0=t0yZ,E=EyZ,p2=p2yZ)
    roots={name:dict(rows) for name,rows in source['roots'].items()}
    for name in roots:roots[name][YZ]=derivatives[name]
    qrows=dict(qr);qrows[YZ]=qyZ
    V=generic['common_velocity_V_axial5'];Vy=dual['common_velocity_V_axial5']
    velocity={C0:V[0],Y:Vy[0].dy,Z:V[1],YZ:Vy[1].dy}
    phase.previous.same_source(f,[row for rows in roots.values() for row in rows.values()]+list(qrows.values())+list(velocity.values()))
    return dict(q=q,roots=roots),qrows,velocity,dict(first_y_Z_source=record,
        genuine_mixed_raw_source_rows=dict(E_yZ=EyZ,C_yZ=CyZ,B_yZ=ByZ,I_axial_before_R_yZ=IyZ),
        original_correlated_shear_yZ=dict(F_yZ=FyZ,H_Z=HZ,H_yZ=HyZ,a_yZ=ayZ,b_yZ=byZ,t0_yZ=t0yZ,
            Delta_yZ=Delta_yZ,q_squared_yZ=q2yZ,q_yZ=qyZ),
        original_pre_R_p2_yZ=baryZ,original_full_p2_yZ=p2yZ,original_radius_y=R,original_radius_Z_exact_zero=True,
        exact_radius_mixed_recipe='p2_yZ=R*(p2bar_yZ+p2bar_Z), R_y=R and R_Z=0',
        raw_yZ_rows_not_relabelled_first_derivatives=True,three_raw_axial_coefficients_cover_exported_mixed_yZ=True,
        active_q_sigma_one_and_fixed_eta_dstar_required=True,
        genuine_yy_ZZ_or_total_spatial_yZ_installed=False,**dict.fromkeys(fields.previous.OPEN,False))


def all_u_mixed_bounds(f,source,qrows,dstar_log,phi):
    """Original implicit mixed calculus with all-finite-u Fourier bounds."""
    c=f.c;phibox=c.mpf(phi)
    if not 0<=ep(phibox)[0]<=ep(phibox)[1]<=1:raise ValueError('Original closed phi source interval required')
    if source['q'] is not qrows[C0] or any(set(rows)!=set(ORDERS) for rows in source['roots'].values()) or set(qrows)!=set(ORDERS):
        raise ValueError('Same genuine mixed original roots/q required')
    phase.previous.same_source(f,[row for rows in source['roots'].values() for row in rows.values()]+list(qrows.values()))
    roots={name:Mixed(rows) for name,rows in source['roots'].items()};q=Mixed(qrows)
    dstar=f.factor((0,0,0,0,0),dstar_log);u=roots['p2']*q
    u=Mixed({key:row.positive_divide(dstar,dstar_log) for key,row in u.rows.items()})
    absolute=lambda row:previous.magnitude(f,row)
    absjet=lambda jet:{key:absolute(row) for key,row in jet.rows.items()}
    a,t0,E=roots['a'],roots['t0'],roots['E'];qa,ua,ta,aa,ea=(absjet(jet) for jet in (q,u,t0,a,E))
    one=f.scalar(1);pi=c.pi
    nu=1+phase.first.current.square(t0[C0])+phase.first.current.square(q[C0])*2
    nuYZ=t0[Y]*t0[Z]*2+t0[C0]*t0[YZ]*2+q[Y]*q[Z]*4+q[C0]*q[YZ]*4
    G={key:qa[key]+qa[C0]*ua[key] for key in (Y,Z)}
    D={key:nu*ta[key]*2+nu*G[key]*(2*c.sqrt(2))+absolute(t0[C0]*t0[key])*2+
        absolute(q[C0]*q[key])*4 for key in (Y,Z)}
    psi={key:D[key]*(2*pi) for key in (Y,Z)}
    # Absolute Fourier coefficients: ||g||1<=1+2|u|,
    # ||g_u||1<=6, ||g_uu||1<=40, ||(k+1)g||1<=4h^3.
    peak=one+ua[C0]*2
    osc={key:qa[key]*peak*2+qa[C0]*ua[key]*12 for key in (Y,Z)}
    oscYZ=qa[YZ]*peak*2+(qa[Y]*ua[Z]+qa[Z]*ua[Y]+qa[C0]*ua[YZ])*12+qa[C0]*ua[Y]*ua[Z]*80
    tpoint={key:ta[key]+osc[key] for key in (Y,Z)}
    tpsi=qa[C0]*(one+ua[C0])*(one+ua[C0])*(one+ua[C0])*8
    L2={key:ta[key]*c.sqrt(2*pi)+G[key]*(2*c.sqrt(pi)) for key in (Y,Z)}
    oscillatoryL2YZ=(qa[YZ]+qa[Y]*ua[Z]+qa[Z]*ua[Y]+qa[C0]*ua[YZ]+qa[C0]*ua[Y]*ua[Z]*40)*(2*c.sqrt(pi))
    L2YZ=ta[YZ]*c.sqrt(2*pi)+oscillatoryL2YZ
    FYZ=absolute(nuYZ)*(2*pi)+(L2[Y]*L2[Z]+nu*L2YZ*c.sqrt(2*pi))*2
    psiYZ=FYZ+tpoint[Y]*psi[Z]+tpoint[Z]*psi[Y]+tpsi*psi[Y]*psi[Z]
    AYZ=aa[YZ]*c.mpf('.5')+(aa[Y]*psi[Z]+aa[Z]*psi[Y]+aa[C0]*psiYZ)*(c.mpf(1)/(4*pi))
    J={C0:qa[C0]*(c.mpf(1)/c.sqrt(2))}
    for key in (Y,Z):J[key]=G[key]*(c.mpf(1)/c.sqrt(2))+(one*c.mpf('.25')+ta[C0]*c.mpf('.5'))*D[key]
    J[YZ]=oscillatoryL2YZ*(c.mpf(1)/(2*c.sqrt(2*pi)))+(
        osc[Y]*psi[Z]+osc[Z]*psi[Y]+tpsi*psi[Y]*psi[Z]+qa[C0]*peak*psiYZ*2)*(c.mpf(1)/(4*pi))
    original=first.all_u_first_bounds(f,source,qrows,dstar_log,phibox)
    Ap={C0:absolute(original['values']['A']),Y:absolute(original['values']['A_y']),
        Z:absolute(original['values']['A_Z']),YZ:AYZ}
    M=Mixed(ta)*Mixed(Ap)+Mixed(aa)*Mixed(J)
    BYZ=(Mixed(ea)*M)[YZ]
    values=dict(original['values'],A_yZ=previous.symmetric(f,AYZ),B_yZ_over_Pstar=previous.symmetric(f,BYZ))
    if ep(phibox)[0]==ep(phibox)[1] and ep(phibox)[0] in (0,mp.mpf('.5'),1):
        values['A_yZ']=f.scalar(0);values['B_yZ_over_Pstar']=f.scalar(0)
    return dict(values=values,record=dict(status='original_all_signed_u_fixed_phi_mixed_derivatives_enclosed',
        original_genuine_mixed_u_rows=u.rows,original_genuine_mixed_q_rows=qrows,
        Fourier_coefficient_L1_bounds=dict(g='1+2|u|',g_u=6,g_uu=40,angle_weighted_g='4*(1+|u|)^3'),
        original_F_yZ_cap=FYZ,original_inverse_angle_y_Z_caps=psi,original_inverse_angle_yZ_cap=psiYZ,
        original_J_mixed_caps=J,original_A_yZ_cap=AYZ,original_B_yZ_over_Pstar_cap=BYZ,
        implicit_equation='F=psi+T2-2pi*phi*nu=0, F_psi=1+t^2>=1',
        original_inverse_yZ_recipe='-(F_yZ+2t*t_y*psi_Z+2t*t_Z*psi_y+2t*t_psi*psi_y*psi_Z)/(1+t^2)',
        endpoint_ratio='2|t|/(1+t^2)<=1',original_J_full_mixed_endpoint_terms_retained=True,
        original_B_full_E_t0_A_a_J_mixed_product_retained=True,nonzero_t0_and_b_allowed=True,
        every_finite_signed_u_and_zero_included=True,derivatives_enclose_original_functions_not_cap_derivatives=True,
        narrow_numerical_inverse_yZ_or_actual_spatial_mixed_chain_not_claimed=True,
        **dict.fromkeys(fields.previous.OPEN,False)))


def leading_mixed(E,V,values):
    A=Mixed({C0:values['A'],Y:values['A_y'],Z:values['A_Z'],YZ:values['A_yZ']})
    B=Mixed({C0:values['B_over_Pstar'],Y:values['B_y_over_Pstar'],Z:values['B_Z_over_Pstar'],YZ:values['B_yZ_over_Pstar']})
    E,V=Mixed(E),Mixed(V);EA=E*A;EEA=E*EA
    return {key:row.rows for key,row in dict(m=B,h=EA,k=V*EA+E*B,e=V*B*2-EEA,p=EEA).items()}


class OriginalRmMixedYZPrimitives:
    mode='actual_Rm_genuine_yZ_source_roots_and_general_original_all_u_fixed_phi_mixed_A_B_and_leading_densities'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalRmWeightedAveraging(dps);self.c=self.upstream.c;self.family=self.upstream.family
        self.hashes=dict(self.upstream.hashes);fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted actual Rm genuine mixed primitive receipt required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def query(self,label,left,right=None,phi=(0,1)):
        p=self.upstream.upstream.phase;mixed=p.upstream.upstream.upstream.owner(label);op=mixed.op;f=op.flow;c=op.c
        with mp.workdps(c.dps+40):
            patch=mixed.evaluate(left) if right is None else mixed.cell(left,right)
            generic=p.upstream.upstream.evaluate(label,left) if right is None else p.upstream.upstream.cell(label,left,right)
            quotient=phase.previous.recover_quotients(op,generic,p.upstream.eta_log,p.upstream.dstar_log,patch=patch)
            source,qr,V,record=source_mixed_frame(op,patch,generic,quotient)
            phibox=c.mpf(phi);got=all_u_mixed_bounds(f,source,qr,p.upstream.dstar_log,phibox)
            leading=leading_mixed(source['roots']['E'],V,got['values'])
        return dict(source_family=self.family,source_frame=label,source_geometry=generic['source_geometry'],
            exact_common_P0_axial5=op.P0,exact_source_Rm_factor=op.Rm_factor,original_fixed_phi_box=phibox,
            original_genuine_mixed_source_roots=source['roots'],original_genuine_mixed_q_rows=qr,
            original_genuine_mixed_V_rows=V,original_mixed_source_record=record,
            original_fixed_phi_A_B_mixed_values=got['values'],original_mixed_primitive_proof=got['record'],
            original_five_signed_leading_density_C0_y_Z_yZ=leading,
            genuine_original_yZ_source_and_primitive_enclosures_installed=True,
            actual_spatial_mixed_yZ_or_Z_averaged_integrals_installed=False,
            **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalRmMixedYZPrimitives(require_checked=False);frames={}
    for label in ('0','.5'):
        frames[label]=dict(active=owner.query(label,(5,4)),terminal=owner.query(label,(2,1)),Rh=owner.query(label,'Rh'),
            whole_active=owner.query(label,(1,1),(71,40)),whole_terminal=owner.query(label,(71,40),'Rh'))
        print('Actual Rm genuine mixed yZ roots and original primitive bounds',label,flush=True)
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},frames=serialized(frames),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()
