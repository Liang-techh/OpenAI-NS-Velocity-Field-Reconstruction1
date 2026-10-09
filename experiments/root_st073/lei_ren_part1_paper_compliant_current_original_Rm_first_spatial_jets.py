"""Genuine actual Rm first y/Z source jets and the full spatial N chain.

Raw mixed4 rows differentiate the original generic recovery, while the
correlated shear identity supplies a_y and the active q branch supplies
q_y. No y row is inferred from a Z row. All-u bounds cover sign-crossing
source cells; numerical original inverse jets are used where executable.
"""
import ast
import copy
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rm_all_u_density_integrals as previous

phase=previous.phase;generic_module=phase.previous.previous
fields,base,ep=previous.fields,previous.base,previous.ep
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_Rm_first_spatial_jets.json.gz'
RECEIPT=PREFIX+'current_original_Rm_first_spatial_jets_check.json'
GATE='original_actual_Rm_genuine_first_y_Z_phase_and_five_density_spatial_chains_installed'


class FirstY:
    """An ordinary y tangent in the same canonical source algebra."""
    def __init__(self,value,derivative=None):
        self.v=value;self.dy=value.scalar(0) if derivative is None else derivative
        phase.previous.same_source(type('Algebra',(),dict(c=value.ctx,logs=value.scale.bases,ledger=value.ledger))(),[value,self.dy])
    @property
    def scale(self):return self.v.scale
    @property
    def ledger(self):return self.v.ledger
    def coerce(self,other):
        if isinstance(other,FirstY):self.v.coerce(other.v);return other
        return FirstY(self.v.coerce(other))
    def __neg__(self):return FirstY(-self.v,-self.dy)
    def __add__(self,other):
        other=self.coerce(other);return FirstY(self.v+other.v,self.dy+other.dy)
    __radd__=__add__
    def __sub__(self,other):return self+-self.coerce(other)
    def __rsub__(self,other):return self.coerce(other)+-self
    def __mul__(self,other):
        other=self.coerce(other);return FirstY(self.v*other.v,self.dy*other.v+self.v*other.dy)
    __rmul__=__mul__


class FirstYFlow:
    def __init__(self,flow):self.flow=flow
    def __getattr__(self,name):return getattr(self.flow,name)
    def scalar(self,value):return FirstY(self.flow.scalar(value))
    def factor(self,*args,**kwargs):return FirstY(self.flow.factor(*args,**kwargs))
    def jet(self,value):return [FirstY(row) for row in self.flow.jet(value)[:3]]


class FirstYOp:
    def __init__(self,op):self.original=op;self.flow=FirstYFlow(op.flow)
    def __getattr__(self,name):return getattr(self.original,name)


def lift_generic_y(op,packet):
    """Differentiate the unchanged recovery; raw rows already are Dy rows.

    Three axial Taylor coefficients suffice for the exported first y and
    first Z inputs. The one final physical R factor is handled explicitly
    below when forming p2_y, rather than duplicated in the raw sources.
    """
    raw=packet['raw_current_radius_y_derivative_axial_coefficients']
    def lift(rows):
        if len(rows)<2 or min(len(rows[0]),len(rows[1]))<3:
            raise ValueError('Genuine source y and axial rows required')
        return [FirstY(value,derivative) for value,derivative in zip(rows[0][:3],rows[1][:3])]
    velocity={}
    for name in ('theta','axial'):
        rows=raw['velocity'][name]
        if len(rows)<3:raise ValueError('Genuine second y velocity row required')
        velocity[name]=[lift(rows),lift(rows[1:])]
    lifted=dict(packet,raw_current_radius_y_derivative_axial_coefficients=dict(
        histories={key:[lift(rows)] for key,rows in raw['histories'].items()},velocity=velocity))
    return generic_module.recover_inputs(FirstYOp(op),lifted)


def source_first_frame(op,patch,generic,quotients):
    """Actual C0 and independent y/Z tangents of all original loop roots."""
    f=op.flow;c=op.c
    if any(packet['common_original_P0_axial5'] is not op.P0 for packet in (generic,quotients)):
        raise ValueError('Same live actual P0 source required')
    if patch['original_P0_normalized_axial5'] is not op.P0 or patch['geometry']!=generic['source_geometry'] or generic['source_geometry']!=quotients['source_geometry']:
        raise ValueError('Same actual radial source geometry and pressure owner required')
    dual=lift_generic_y(op,patch)
    for name in ('E','C','B'):
        a=dual['actual_generic_source_numerators'][name][0].v
        b=generic['actual_generic_source_numerators'][name][0]
        if base.encoded(fields.serialized(a))!=base.encoded(fields.serialized(b)):
            raise ValueError('Lifted source value differs from unchanged actual recovery: '+name)
    n=generic['actual_generic_source_numerators'];proofs=quotients['source_positive_theorems']
    E=n['E'][0];C=n['C'][0];B=n['B'][0]
    Ey=dual['common_velocity_E_axial5'][0].dy
    Cy=dual['actual_generic_source_numerators']['C'][0].dy
    By=dual['actual_generic_source_numerators']['B'][0].dy
    shear=quotients['shear_quotients_axial5'];a,b,t0=(shear[key][0] for key in ('a','b','t0'))
    def divide(value,denominator,proof):
        if proof['source_row'] is not denominator:raise ValueError('Exact same positive source denominator required')
        return value.positive_divide(denominator,proof['source_log_lower'])
    # a=.8-2F/H, F=sum control*(x*gamma_x-gamma/10).
    x=phase.previous.source_coordinate(c,patch['geometry'])
    gamma=patch['actual_gamma_ordinary_x_derivatives']
    H=patch['actual_H_x_derivative_axial5'][0][0]
    Hy=patch['actual_H_x_derivative_axial5'][1][0]*x
    correlated=quotients['actual_correlated_shear_source']
    F=correlated['angular_correction_numerator'][0]
    Fy=sum((op.controls[j+2][0]*(x*(c.mpf(9)/10*gamma[j][1]+x*gamma[j][2]))
            for j in range(3)),f.scalar(0))
    ay=divide(Fy-F*divide(Hy,H,correlated['positive_H_source_theorem']),H,
              correlated['positive_H_source_theorem'])*(-2)
    by=divide(By-b*Ey,E,proofs['E'])
    t0y=divide(-By-t0*Cy,C,proofs['C'])
    b2_over_a=divide(phase.first.current.square(b),a,proofs['a'])
    Delta_y=ay+divide(b*by*2-b2_over_a*ay,a,proofs['a'])
    loop=quotients['original_shear_q'];qjet=loop['q_axial_coefficients']
    if qjet is None or ep(shear['kappa_minus2'][0].coefficient)[1]>=0:
        raise ValueError('Actual sigma=1 source theorem required for genuine q_y')
    q=qjet[0];q2=loop['q_squared_axial_coefficients'][0]
    twice_a=a*2;twice_q=q*2
    q2y=(-Delta_y-q2*ay*2).positive_divide(twice_a,phase.previous.positive_source(f,twice_a,'2a')['source_log_lower'])
    qy=q2y.positive_divide(twice_q,phase.previous.positive_source(f,twice_q,'2q')['source_log_lower'])
    sectors=dual['full_signed_inertial_sectors_axial4']
    I_y=sectors['axial_linear'][0].dy+sectors['axial_quadratic'][0].dy*op.Pstar
    p2bar=quotients['full_signed_inertial_quotients_before_shared_R_axial4']['p2'][0]
    p2bar_y=divide(I_y-p2bar*Ey,E,proofs['E'])
    R=quotients['exact_same_shared_positive_radius_factor']
    p2y=R*(p2bar_y+p2bar)
    roots={name:{(0,0):rows[0],(0,1):rows[1],(1,0):derivative}
           for name,rows,derivative in (
               ('a',shear['a'],ay),('t0',shear['t0'],t0y),
               ('E',generic['common_velocity_E_axial5'],Ey),
               ('p2',quotients['actual_factored_p1_p2_axial4']['p2'],p2y))}
    qr={(0,0):q,(0,1):qjet[1],(1,0):qy}
    phase.previous.same_source(f,[row for rows in roots.values() for row in rows.values()]+list(qr.values()))
    return dict(q=q,roots=roots),qr,dict(source_geometry=generic['source_geometry'],
        original_common_P0_axial5=op.P0,original_P0_y_exactly_zero=True,
        original_generic_recovery_y_tangents=dict(E_y=Ey,C_y=Cy,B_y=By,
            inertial_axial_before_R_y=I_y,p2_before_R_y=p2bar_y),
        original_correlated_shear_y=dict(F=F,F_y=Fy,H=H,H_y=Hy,a_y=ay,b_y=by,t0_y=t0y,Delta_y=Delta_y,q_squared_y=q2y,q_y=qy),
        original_physical_radius=R,original_physical_radius_y=R,
        original_full_p2_y=p2y,exact_R_y_recipe='p2_y=R*(p2bar_y+p2bar)',
        genuine_y_from_mixed4_raw_rows_not_Z_surrogate=True,
        log_rows_already_include_x_derivatives_no_second_conversion=True,
        q_y_from_original_active_branch_with_fixed_eta=True,
        axial_export_order='C0 and ordinary first Z/y only; generic lift uses Z0..2',
        **dict.fromkeys(fields.previous.OPEN,False))


def compile_all_u_y_bounds():
    tree=ast.parse(Path(previous.__file__).read_text(encoding='utf8'))
    fn=copy.deepcopy(next(n for n in tree.body if getattr(n,'name',None)=='all_u_primitive_bounds'))
    indices=[];labels=[]
    class Direction(ast.NodeTransformer):
        def visit_Tuple(self,node):
            self.generic_visit(node)
            if all(isinstance(v,ast.Constant) for v in node.elts) and [v.value for v in node.elts]==[0,1]:
                node.elts=[ast.Constant(1),ast.Constant(0)];indices.append('genuine y source index')
            return node
        def visit_keyword(self,node):
            self.generic_visit(node)
            if node.arg and '_Z' in node.arg:
                labels.append(node.arg);node.arg=node.arg.replace('_Z','_y')
            return node
        def visit_Constant(self,node):
            if isinstance(node.value,str):node.value=node.value.replace('C0/Z','C0/y')
            return node
    fn=Direction().visit(fn)
    if len(indices)!=3 or len(labels)!=7:raise ValueError('Accepted all-u first-direction template sites changed')
    env=dict(vars(previous))
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
                 '<original all-u proof; genuine y direction>','exec'),env)
    return env['all_u_primitive_bounds'],dict(only_genuine_derivative_index_and_axis_labels_changed=True,
        original_all_u_mathematical_operations_unchanged=True,index_edits=indices,axis_labels=labels)


ALL_U_Y,Y_BINDING=compile_all_u_y_bounds()


def all_u_first_bounds(f,source,qr,dstar_log,phi):
    projected=lambda index:(dict(q=source['q'],roots={name:{(0,0):rows[(0,0)],index:rows[index]}
        for name,rows in source['roots'].items()}),{(0,0):qr[(0,0)],index:qr[index]})
    zsource,zqr=projected((0,1));ysource,yqr=projected((1,0))
    z=previous.all_u_primitive_bounds(f,zsource,zqr,dstar_log,phi)
    y=ALL_U_Y(f,ysource,yqr,dstar_log,phi)
    values=dict(z['values'],A_y=y['values']['A_y'],B_y_over_Pstar=y['values']['B_y_over_Pstar'])
    return dict(values=values,record=dict(status='enclosed_by_original_all_u_first_y_Z_theorem',
        original_Z_majorants=z['record'],original_y_majorants=y['record'],
        genuine_directional_source_projection=True,all_u_y_template_binding=Y_BINDING))


def density_derivative(E,Ed,V,Vd,primitives,N,direction):
    """The original five signed product rules for a genuine first direction."""
    c=E.ctx;factor=c.mpf(1)/phase.candidate_N(N)
    A=primitives['A']*factor;Ad=primitives['A_'+direction]*factor
    increment=phase.densities.density.factored_expm1(A)
    exponential=c.exp(phase.first.phase.bounded_value(A))
    dE=E*increment;dEd=Ed*increment+E*exponential*Ad
    dV=primitives['B_over_Pstar']*factor;dVd=primitives['B_'+direction+'_over_Pstar']*factor
    kernels=phase.densities.density.signed_density_kernels(E,V,dE,dV)
    cross=Ed*dE+E*dEd;half_square=dE*dEd
    derivatives=dict(m=dVd,h=dEd,
        k=Vd*dE+V*dEd+Ed*dV+E*dVd+dEd*dV+dE*dVd,
        e=Vd*dV*2+V*dVd*2+dV*dVd*2-cross-half_square,p=cross+half_square)
    return dict(kernels=kernels,derivatives=derivatives,original_direction=direction,
                original_velocity_increment_derivatives=dict(deltaE_d=dEd,deltaV_d=dVd))


class OriginalRmFirstSpatialJets:
    mode='actual_Rm_genuine_first_y_Z_phi_source_primitives_and_density_fast_spatial_chain'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalRmAllUDensityIntegrals(dps);self.phase=self.upstream.phase
        self.c=self.upstream.c;self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes)
        fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted actual Rm genuine first spatial jet receipt required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def query(self,label,left,right=None,N=257):
        N=phase.candidate_N(N);op=self.phase.upstream.upstream.upstream.owner(label).op;f=op.flow
        mixed=self.phase.upstream.upstream.upstream.owner(label)
        with mp.workdps(self.c.dps+40):
            if right is None:
                patch=mixed.evaluate(left);generic=self.phase.upstream.upstream.evaluate(label,left)
            else:
                patch=mixed.cell(left,right);generic=self.phase.upstream.upstream.cell(label,left,right)
            quotients=phase.previous.recover_quotients(op,generic,self.phase.upstream.eta_log,
                self.phase.upstream.dstar_log,patch=patch)
            source,qr,source_record=source_first_frame(op,patch,generic,quotients)
            if label not in self.phase.mappers:self.phase.mappers[label]=phase.RmRadiusPhase(op,self.family,self.phase.parameter_family)
            mapper=self.phase.mappers[label];radius=mapper.point(left,N) if right is None else mapper.cell(left,right,N)
            cells=[]
            Ey=source['roots']['E'][(1,0)];V=generic['common_velocity_V_axial5']
            Vy=lift_generic_y(op,patch)['common_velocity_V_axial5'][0].dy
            for phi in radius['phase_boxes']:
                got=phase.first.conditioned_first_jets(source,qr,self.phase.upstream.dstar_log,phi)
                route='original full conditioned first jets'
                if got['values'] is None:
                    got=all_u_first_bounds(f,source,qr,self.phase.upstream.dstar_log,phi)
                    route='original all-u genuine y/Z enclosures'
                values=got['values']
                spatial=dict(values,A_y=values['A_y']+values['A_phi']*N,
                    B_y_over_Pstar=values['B_y_over_Pstar']+values['B_phi_over_Pstar']*N)
                y_density=density_derivative(source['roots']['E'][(0,0)],Ey,V[0],Vy,spatial,N,'y')
                Z_density=density_derivative(source['roots']['E'][(0,0)],source['roots']['E'][(0,1)],
                    V[0],V[1],values,N,'Z')
                slow_density=density_derivative(source['roots']['E'][(0,0)],Ey,V[0],Vy,values,N,'y')
                cells.append(dict(source_family=self.family,source_geometry=generic['source_geometry'],candidate_N=N,
                    source_first_y_Z=source_record,original_primitive_route=route,original_first_primitive_record=got['record'],
                    original_fixed_phi_first_primitive_values=values,actual_spatial_first_primitive_values=spatial,
                    original_five_signed_density_C0=y_density['kernels'],
                    original_five_signed_density_fixed_phi_y=slow_density['derivatives'],
                    actual_five_signed_density_y=y_density['derivatives'],actual_five_signed_density_Z=Z_density['derivatives'],
                    exact_common_P0_axial5=op.P0,exact_same_shared_radius_factor=quotients['exact_same_shared_positive_radius_factor'],
                    actual_phase_y_exactly_N=True,actual_phase_Z_exactly_zero=True,
                    genuine_first_y_Z_only_not_mixed_second_jets=True,
                    **dict.fromkeys(fields.previous.OPEN,False)))
            return dict(source_family=self.family,source_frame=label,candidate_N=N,
                actual_original_Rm_radius_phase=radius,actual_first_spatial_source_cells=cells,
                genuine_first_y_Z_source_and_density_chain_installed=bool(cells),
                genuine_mixed_second_derivatives_or_sharp_integrals_installed=False,
                **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalRmFirstSpatialJets(require_checked=False);frames={}
    for label in ('0','.5'):
        frames[label]=dict(active=owner.query(label,(5,4)),terminal=owner.query(label,(2,1)),
            Rh=owner.query(label,'Rh'),whole_active=owner.query(label,(1,1),(71,40)),
            whole_terminal=owner.query(label,(71,40),'Rh'))
        print('Actual Rm genuine first y/Z and spatial density chain',label,flush=True)
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_all_u_y_template_binding=Y_BINDING,frames=fields.serialized(frames),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),
                                      compresslevel=6,mtime=0));return report


if __name__=='__main__':run()
