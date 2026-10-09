"""Full signed generic shear/inertial inputs from the actual leading Rm owner.

The five original formal bases and every directed offset remain in the
source algebra. No legacy patch owner, modal projection, quotient cover,
finite-N pulse, or selected source value is used.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rm_patch_mixed4_cells as previous
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

fields,base,ep=previous.fields,previous.base,previous.ep
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_Rm_generic_inputs.json.gz'
RECEIPT=PREFIX+'current_original_Rm_generic_inputs_check.json'
GATE='original_actual_two_frame_Rm_full_signed_generic_inputs_installed'


def source_bindings():
    from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
    bound=assignment_source_bindings('current_generic_shear_moment_recovery','field',{
        'L':'1-Z*Z*de','d':'1-Z*Z','pressure':'self.P0+p',
        'transport':'m*(Z*(1-de))+axial_derivative(m)*d',
        'Q':'(V*(2*Z)-transport)/L',
        'theta_linear':'(-E+h*(1-de/2)-axial_derivative(h)*(Z*((1-de)/2)))/L',
        'theta_quadratic':'(k*(Z*(2*de-1))-axial_derivative(k)*d+E*transport)/L',
        'axial_linear':'(-V+(m-axial_derivative(m)*Z)*((1-de)/2))/L',
        'axial_quadratic':'(V*transport+e*(Z*(2*de))-axial_derivative(e)*d+pressure*(Z*(2*(1+de)))-axial_derivative(pressure)*d)/L'})
    bound.update(assignment_source_bindings('current_generic_shear_inputs','from_packet',{
        'It':"with_R(field['inertial_theta_linear']+algebra.shift(field['inertial_theta_quadratic'],(0,.5,0,0)))",
        'Iz':"with_R(field['inertial_axial_linear']+algebra.shift(field['inertial_axial_quadratic'],(0,.5,0,0)))",
        'den':'C*E','kap':'C*C+B*B','excess':'kap-2*den',
        'Hnum':'C*It-B*Iz','stronger':'Hnum-2*den','D':'Hnum-kap','J':'C*Iz+B*It'}))
    return dict(passed=True,unchanged_original_full_recovery_assignments=bound,
        no_signed_sector_or_pressure_term_removed=True,
        source_jet_derivative_consumes_one_axial_order=True,
        canonical_five_bases_and_directed_offsets_preserved=True)


def recover_inputs(op,packet):
    """Same signed source equations with order-aware canonical row arithmetic."""
    f=op.flow;c=op.c
    def add(*rows):
        count=min(map(len,rows))
        return [sum((row[n] for row in rows),f.scalar(0)) for n in range(count)]
    def mul(a,b):
        return [sum((a[j]*b[n-j] for j in range(n+1)),f.scalar(0))
                for n in range(min(len(a),len(b)))]
    scale=lambda rows,value:[v*value for v in rows]
    dz=lambda rows:[rows[n+1]*(n+1) for n in range(len(rows)-1)]
    z=op.zrows;zjet=op.reference.z;de=op.reference.delta
    d=f.jet(1-zjet*zjet);L=1-zjet*zjet*de
    if ep(L[0])[0]<=0:raise ValueError('Actual positive original meridional denominator required')
    invL=f.jet(L.reciprocal());divideL=lambda rows:mul(rows,invL)
    raw=packet['raw_current_radius_y_derivative_axial_coefficients']
    invS=f.factor((0,-.5,0,0,0));S=op.Pstar
    E=raw['velocity']['theta'][0];Ey=raw['velocity']['theta'][1]
    V=scale(raw['velocity']['axial'][0],invS);Vy=scale(raw['velocity']['axial'][1],invS)
    histories={key:scale(rows[0],invS) if key in ('m','k') else rows[0]
               for key,rows in raw['histories'].items()}
    m,h,k,e,p=(histories[key] for key in ('m','h','k','e','p'))
    if packet['original_P0_normalized_axial5'] is not op.P0:
        raise ValueError('Exact same live P0 object required')
    pressure=add(op.P0,p)
    transport=add(scale(mul(m,z),1-de),mul(dz(m),d))
    Q=divideL(add(scale(mul(V,z),2),scale(transport,-1)))
    itl=divideL(add(scale(E,-1),scale(h,1-de/2),scale(mul(dz(h),z),-(1-de)/2)))
    itq=divideL(add(scale(mul(k,z),2*de-1),scale(mul(dz(k),d),-1),mul(E,transport)))
    izl=divideL(add(scale(V,-1),scale(add(m,scale(mul(dz(m),z),-1)),(1-de)/2)))
    izq=divideL(add(mul(V,transport),scale(mul(e,z),2*de),scale(mul(dz(e),d),-1),
        scale(mul(pressure,z),2*(1+de)),scale(mul(dz(pressure),d),-1)))
    C=add(E,scale(Ey,-2));B=scale(Vy,2)
    geometry=packet['geometry'];point=geometry['point']
    if point:
        x=c.exp(1) if geometry.get('exact_x_source')=='exp(1)' else c.mpf(geometry['exact_x'][0])/geometry['exact_x'][1]
    else:
        a=c.mpf(geometry['exact_left'][0])/geometry['exact_left'][1]
        b=c.exp(1) if geometry.get('exact_right_source')=='exp(1)' else c.mpf(geometry['exact_right'][0])/geometry['exact_right'][1]
        x=c.mpf([ep(a)[0],ep(b)[1]])
    R=op.Rm_factor*x
    It=scale(add(itl,scale(itq,S)),R);Iz=scale(add(izl,scale(izq,S)),R)
    den=mul(C,E);kap=add(mul(C,C),mul(B,B));excess=add(kap,scale(den,-2))
    H=add(mul(C,It),scale(mul(B,Iz),-1))
    D=add(H,scale(kap,-1));J=add(mul(C,Iz),mul(B,It))
    result=dict(common_velocity_E_axial5=E,common_velocity_V_axial5=V,
        common_own_five_histories_axial5=histories,common_original_P0_axial5=op.P0,
        common_absolute_pressure_axial5=pressure,common_radial_Q_axial4=Q,
        actual_generic_source_numerators=dict(E=E,C=C,B=B,inertial_theta=It,inertial_axial=Iz,
            positive_denominator=den,kappa=kap,kappa_minus2=excess,
            H0_minus2=add(H,scale(den,-2)),D=D,J=J),
        full_signed_inertial_sectors_axial4=dict(theta_linear=itl,theta_quadratic=itq,
            axial_linear=izl,axial_quadratic=izq),
        quotient_recipes=dict(a='C/E',b='B/E',t0='-B/C',p1='inertial_theta/E',
            p2='inertial_axial/E',kappa='kappa/(C*E)',D='D/(C*E)',J='J/(C*E)'),
        full_quadratic_recipe='(2*D^2*(C*E)-(kappa-2*C*E)*J^2)/(C*E)^3',
        original_physical_radius=R,source_geometry=geometry,
        source_log_basis_names=('log_hb','2_log_Pstar','log_F0_squared','log_Ra','zero'),
        original_fixed_source_log_bases=f.logs,source_ledger_is_same_object=True,
        actual_patch_source_owner=type(op).__name__,original_P0_is_same_live_object=True,
        derivative_orders=dict(profile_and_histories=5,radial_and_full_inertial=4),
        full_inertial_linear_quadratic_pressure_meridional_sectors_retained=True,
        no_source_basis_projection_or_directed_offset_materialization=True,
        signed_source_expressions_not_quotient_bounds=True,
        actual_positive_E_C_and_generic_cone_admission_certified=False,
        source_frame_conditional_on_same_actual_Rm_inlet=True,
        actual_phase_or_finite_N_density_integrals_installed=False,
        **dict.fromkeys(fields.previous.OPEN,False))
    for rows in (*histories.values(),E,V,pressure,op.P0,Q,itl,itq,izl,izq,*result['actual_generic_source_numerators'].values()):
        if any(v.scale.bases is not f.logs or v.ledger is not f.ledger for v in rows):
            raise ValueError('One live five-basis source algebra and ledger required')
    return result


class OriginalRmGenericInputs:
    mode='actual_leading_Rm_owner_full_signed_generic_shear_and_inertial_source_inputs'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalRmPatchMixed4Cells(dps);self.c=self.upstream.c
        self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes);self.bindings=source_bindings()
        bind=fields.previous.bind
        for name in (Path(__file__).name,PREFIX+'current_generic_shear_moment_recovery.py',PREFIX+'current_generic_shear_inputs.py'):
            bind(self.hashes,name,sha(name))
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted full signed actual Rm generic-input receipt required')
            for name,digest in receipt['input_hashes'].items():bind(self.hashes,name,digest)
            bind(self.hashes,RECEIPT,sha(RECEIPT))
    def evaluate(self,label,coordinate):
        owner=self.upstream.owner(label)
        with mp.workdps(self.c.dps+40):return recover_inputs(owner.op,owner.evaluate(coordinate))
    def cell(self,label,left,right):
        owner=self.upstream.owner(label)
        with mp.workdps(self.c.dps+40):return recover_inputs(owner.op,owner.cell(left,right))


def run():
    began=time.monotonic();owner=OriginalRmGenericInputs(require_checked=False);frames={}
    for label in ('0','.5'):
        frames[label]=dict(active_point=owner.evaluate(label,(5,4)),terminal_cell=owner.cell(label,(71,40),'Rh'))
        print('Actual Rm full signed generic inputs',label,flush=True)
    result=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},
        original_full_recovery_source_bindings=owner.bindings,frames=fields.serialized(frames),
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(result),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return result


if __name__=='__main__':run()
