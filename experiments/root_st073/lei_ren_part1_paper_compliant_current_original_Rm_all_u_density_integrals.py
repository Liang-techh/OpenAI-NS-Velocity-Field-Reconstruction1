"""Original all-signed-u primitive/Z bounds complete the actual Rm atlas.

The original inverse identity and Poisson L2 identities give analytic
enclosures through u=0, preserving nonzero u_Z. These are bounds on the
original inverse-defined functions, not a replacement inverse or field.
Real incoming histories, sharp moment closure and whole-Z remain open.
"""
import ast
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_Rm_whole_density_integrals as previous

phase=previous.previous.previous
fields,base,ep=previous.fields,previous.base,previous.ep
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_Rm_all_u_density_integrals.json.gz'
RECEIPT=PREFIX+'current_original_Rm_all_u_density_integrals_check.json'
GATE='original_actual_Rm_all_signed_u_Z_bounds_and_complete_local_five_density_integrals_installed'
RATES=previous.RATES


def defining_source_theorem():
    """Bind the original graph and prove the inverse/Poisson algebra."""
    tree=ast.parse(Path(phase.first.__file__).read_text(encoding='utf8'))
    fn=next(n for n in tree.body if getattr(n,'name',None)=='conditioned_first_jets')
    expressions=(
        'nu=1+t0.square()+q.square()*2',
        'F=((1+t0.square())*rad+t0*q*hinv*W1*4+q.square()*s*W2*4).divide(nu,0)*(1/(2*c.pi))',
        'numerator=t0*q*hinv*W1*4+q.square()*(s*W2*4-psi*(4*c.pi))',
        'A=(a*numerator).divide(nu,0)*(1/(4*c.pi))',
        'B=E*(t0*A-a*q*hinv*W1*(1/(2*c.pi)))',
        'direction=loop.t0+(loop.q*loop.hinv*(loop.scalar(cosine)-r.v)*2).positive_divide(D,logD)',
        'Fx=(1+current.square(direction)).positive_divide(nu.v,0)')
    for text in expressions:
        wanted=ast.dump(ast.parse(text).body[0])
        count=2 if text.startswith('Fx=') else 1
        if sum(ast.dump(n)==wanted for n in ast.walk(fn) if isinstance(n,ast.Assign))!=count:
            raise ValueError('Original phase/primitive defining assignment changed: '+text)
    x,T,a,M=sy.symbols('x T a M',real=True);nu=1+M
    F=(2*sy.pi*x+T)/(2*sy.pi*nu)
    A=a*(T-2*sy.pi*x*M)/(4*sy.pi*nu)
    if sy.simplify(A-a*(F-x)/2)!=0:raise ArithmeticError('Original inverse A identity failed')
    r=sy.symbols('r',real=True);s=1-r*r
    S0=1/s;S1=r/s**2;S2=(1+r*r)/s**3
    dr_norm=r*r/s*S0-2*r*S1+s*S2
    if sy.simplify(dr_norm-1/s**2)!=0 or sy.simplify(dr_norm*s**3-s)!=0:
        raise ArithmeticError('Original normalized Poisson derivative identity failed')
    z=sy.symbols('z',real=True);h=sy.sqrt(1+z*z)
    if sy.simplify(sy.diff(z/h,z)-(1+z*z)**sy.Rational(-3,2))!=0:
        raise ArithmeticError('Original r_u identity failed')
    return dict(original_defining_phase_primitive_and_direction_AST_assignments=expressions,
        exact_inverse_identity='A=a/2*(phi-psi_fraction)',
        original_direction='t=t0+2*q*sqrt(s)*(cos(psi)-r)/(1-2*r*cos(psi)+r^2)',
        original_parameters='r=u/sqrt(1+u^2); s=1/(1+u^2); nu=1+t0^2+2*q^2',
        normalized_Poisson_coefficients='sqrt(1-r^2)*r^k, k>=0',
        exact_squared_coefficient_norm=1,exact_squared_u_derivative_coefficient_norm='s<=1',
        original_full_period_integral_t_squared='2*pi*(t0^2+2*q^2)',
        original_full_period_integral_t_u_squared='4*pi*q^2*s',
        original_full_period_integral_t_q_squared='4*pi',
        positive_phase_derivative='Phi_psi_fraction=(1+t^2)/nu>0',
        all_finite_signed_u_including_zero_covered=True,
        source_radius_huge_exponential_or_inverse_bracket_not_materialized=True)


THEOREM=defining_source_theorem()


def absolute(row):
    return phase.first.phase.absolute(row)


def symmetric_majorant(row):
    lo,hi=ep(row.coefficient)
    if lo<0 or not mp.isfinite(hi):raise ValueError('Finite nonnegative formal majorant required')
    return phase.first.prior.ScaledEnclosure(row.scale,row.ctx.mpf([-hi,hi]),row.ledger)


def all_u_primitive_bounds(f,source,qrows,dstar_log,phi):
    """Enclose original C0, fixed-phi Z, and phi derivatives without sign cuts.

    Parseval and Cauchy bound the original t integrals. At the unique
    inverse, |psi_fraction_Z| <= nu*|Phi_Z|. The endpoint factor in J_Z
    uses |(t-t0)/(1+t^2)|<=1/2+|t0|, avoiding a large-u pointwise peak.
    """
    c=f.c;phi=c.mpf(phi)
    if ep(phi)[0]<0 or ep(phi)[1]>1:raise ValueError('Actual closed phase box in [0,1] required')
    roots=source['roots']
    if set(roots)!= {'a','t0','E','p2'} or set(qrows)!={(0,0),(0,1)}:
        raise ValueError('Original genuine C0/Z roots and q rows required')
    rows=[row for values in roots.values() for row in values.values()]+list(qrows.values())
    phase.previous.same_source(f,rows)
    if source['q'] is not qrows[(0,0)]:raise ValueError('Same live q C0 source required')
    a,t0,E,p2=(roots[name][(0,0)] for name in ('a','t0','E','p2'))
    az,tz,Ez,pz=(roots[name][(0,1)] for name in ('a','t0','E','p2'))
    q,qz=qrows[(0,0)],qrows[(0,1)]
    dstar=f.factor((0,0,0,0,0),dstar_log)
    u=(p2*q).positive_divide(dstar,dstar_log)
    uz=(pz*q+p2*qz).positive_divide(dstar,dstar_log)
    aa,at,aq=(absolute(row) for row in (a,t0,q))
    aaz,atz,aqz,auz,aE,aEz=(absolute(row) for row in (az,tz,qz,uz,E,Ez))
    nu=f.scalar(1)+phase.first.current.square(t0)+phase.first.current.square(q)*2
    # sqrt(t0^2+2q^2)<=nu; all terms are positive source majorants.
    nu_Z=at*atz*2+aq*aqz*4
    source_q_change=aqz+aq*auz
    psi_Z=nu*atz*2+nu*source_q_change*(2*c.sqrt(2))+nu_Z
    half=c.mpf(1)/2;inverse_sqrt2=1/c.sqrt(2)
    A=aa*half;A_Z=aaz*half+aa*psi_Z*half
    J=aq*inverse_sqrt2
    J_Z=source_q_change*inverse_sqrt2+(f.scalar(c.mpf(1)/4)+at*half)*psi_Z
    B_over_E=at*A+aa*J
    B=aE*B_over_E
    B_Z=aEz*B_over_E+aE*(atz*A+at*A_Z+aaz*J+aa*J_Z)
    A_phi=aa*(f.scalar(1)+nu)*half
    J_phi=nu*(f.scalar(c.mpf(1)/4)+at*half)
    B_phi=aE*(at*A_phi+aa*J_phi)
    majorants=dict(A=A,B_over_Pstar=B,A_Z=A_Z,B_Z_over_Pstar=B_Z,
                   A_phi=A_phi,B_phi_over_Pstar=B_phi)
    values={key:symmetric_majorant(row) for key,row in majorants.items()}
    phase.previous.same_source(f,[*values.values(),u,uz,nu,psi_Z,J,J_Z])
    return dict(values=values,record=dict(status='enclosed_by_original_all_u_integral_theorem',
        geometry='all finite signed u, including source intervals spanning zero',actual_phase_box=phi,
        original_u=u.record(),original_u_Z=uz.record(),original_nu=nu.record(),
        original_genuine_Z_rows_retained=True,original_positive_phase_inverse_exists_on_each_source_point=True,
        numerical_inverse_bracket_or_selected_field_value_not_claimed=True,
        full_phase_range_majorants_not_phase_samples=True,
        implicit_psi_fraction_Z_majorant=psi_Z.record(),original_J_C0_majorant=J.record(),
        original_J_Z_majorant=J_Z.record(),original_primitive_majorants={key:row.record() for key,row in majorants.items()},
        original_A_B_C0_Z_phi_enclosures={key:row.record() for key,row in values.items()},
        original_inverse_identity_and_Poisson_L2_theorem=THEOREM,
        **dict.fromkeys(fields.previous.OPEN,False)))


class OriginalRmAllUDensityIntegrals:
    mode='actual_Rm_original_all_u_C0_Z_phi_majorants_complete_local_density_integral_atlas'
    def __init__(self,dps=500,require_checked=True):
        self.atlas=previous.OriginalRmWholeDensityIntegrals(dps);self.phase=self.atlas.phase
        self.c=self.atlas.c;self.family=self.atlas.family;self.hashes=dict(self.atlas.hashes);self.cache={}
        fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted actual original all-u Rm density-integral receipt required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def close_source_cell(self,label,unknown,N):
        op=self.phase.upstream.upstream.upstream.owner(label).op;f=op.flow
        geometry=unknown['source_geometry'];left=tuple(geometry['exact_left']);right=tuple(geometry['exact_right'])
        if unknown['source_family']!=self.family or unknown['exact_common_P0_axial5'] is not op.P0:
            raise ValueError('Same actual unresolved Rm source/P0 owner required')
        generic=self.phase.upstream.upstream.cell(label,left,right)
        quotients=self.phase.upstream.cell(label,left,right)
        source,qrows=phase.source_frame(op,generic,quotients)
        spatial=unknown['actual_source_bound_unresolved_phase_density'];radius=spatial['actual_original_Rm_radius_phase']
        if radius['candidate_N']!=N or radius['source_geometry']!=geometry or radius['exact_source_Rm_factor'] is not op.Rm_factor:
            raise ValueError('Same actual source geometry, phase and candidate N required')
        cells=[]
        for phi in radius['phase_boxes']:
            got=all_u_primitive_bounds(f,source,qrows,self.phase.upstream.dstar_log,phi)
            values=got['values'];E,V=generic['common_velocity_E_axial5'],generic['common_velocity_V_axial5']
            density=phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],values,N)
            cells.append(dict(source_geometry=geometry,source_family=self.family,candidate_N=N,
                original_phase_Z_only_result=got['record'],actual_original_primitive_Z_values=values,
                actual_original_five_signed_density_C0_Z=density,original_common_P0_axial5=op.P0,
                exact_same_shared_radius_factor=quotients['exact_same_shared_positive_radius_factor'],
                source_q_Z_includes_all_b_and_b_Z_terms=True,no_fabricated_y_derivative_exported=True,
                original_source_V_already_in_common_Pstar_units=True,source_radius_Z_independent=True,
                phase_is_independent_candidate_parameter=False,original_all_u_analytic_enclosures=True,
                **dict.fromkeys(fields.previous.OPEN,False)))
        closed=dict(source_family=self.family,actual_original_Rm_radius_phase=radius,
            actual_source_bound_phase_density_cells=cells,actual_closed_radial_source_cell=True,
            actual_original_spatial_Z_density_interface_installed=bool(cells),
            actual_spatial_y_derivative_or_sharp_integral_installed=False,
            **dict.fromkeys(fields.previous.OPEN,False))
        integral=previous.integrate_source_cell(f,closed,left,right,previous.JOIN,
            source_family=self.family,P0=op.P0,Rm_factor=op.Rm_factor)
        return dict(original_unresolved_source_geometry=geometry,actual_all_u_closed_density_source=closed,
                    actual_all_u_cell_to_join_integral=integral)

    def contribution(self,label,N=257):
        N=phase.candidate_N(N);key=(label,N)
        if key in self.cache:return self.cache[key]
        with mp.workdps(self.c.dps+40):
            parent=self.atlas.contribution(label,N=N);op=self.phase.upstream.upstream.upstream.owner(label).op;f=op.flow
            closed=[self.close_source_cell(label,cell,N) for cell in parent['actual_unresolved_active_cells']]
            cells=[]
            for entry in parent['complete_active_source_atlas']:
                if entry['density_integral_enclosed']:
                    cell=parent['actual_active_cells_to_join'][entry['resolved_cell_index']]
                else:cell=closed[entry['unresolved_cell_index']]['actual_all_u_cell_to_join_integral']
                if cell['source_geometry']!=entry['source_geometry']:
                    raise ValueError('Actual gap-free same-geometry source replacement required')
                cells.append(cell)
            active={name:[sum((cell['actual_cell_to_target_integral_C0_Z'][name][n] for cell in cells),f.scalar(0))
                          for n in range(2)] for name in RATES}
            whole=previous.previous.affine_transport(f,active,parent['actual_terminal_local_defect_integral_C0_Z'],
                                                    parent['actual_terminal_logarithmic_interval_width'])
            packet=dict(source_family=self.family,source_frame=label,candidate_N=N,
                exact_actual_active_partition=parent['exact_actual_active_partition'],actual_complete_active_cells_to_join=cells,
                actual_all_u_closed_edge_sources=closed,previous_unknown_integral_count=len(closed),
                actual_active_local_defect_integral_C0_Z=active,
                actual_terminal_local_defect_integral_C0_Z=parent['actual_terminal_local_defect_integral_C0_Z'],
                actual_whole_patch_local_defect_integral_C0_Z=whole,
                actual_terminal_source_partition=parent['actual_terminal_source_partition'],
                actual_terminal_logarithmic_interval_width=parent['actual_terminal_logarithmic_interval_width'],
                exact_common_P0_axial5=op.P0,full_original_logarithmic_interval_width=self.c.mpf(1),
                actual_leading_Rm_memory=parent['actual_leading_Rm_memory'],
                actual_leading_join_memory=parent['actual_leading_join_memory'],actual_leading_Rh_memory=parent['actual_leading_Rh_memory'],
                complete_local_whole_patch_C0_Z_integral_bounds_installed=True,unknown_local_source_integral_count=0,
                finite_N_Rm_incoming_defect_is_unsupplied_affine_argument=True,
                exact_incoming_recipe='deltaH_j(Rh)=exp(-lambda_j)*deltaH_j(Rm)+local_integral_j([1,e])',
                all_u_bounds_are_original_function_enclosures_not_a_modified_field=True,
                integral_bound_completeness_not_sharp_defect_or_five_moment_closure=True,
                full_Z_prefix_repair_cone_global_N_and_n_recursion_not_claimed=True,
                **dict.fromkeys(fields.previous.OPEN,False))
            self.cache[key]=packet;return packet

    def transport_supplied_incoming(self,label,incoming,N=257):
        if incoming is None:raise ValueError('Real finite-N correction at Rm cannot be silently reset')
        packet=self.contribution(label,N);op=self.phase.upstream.upstream.upstream.owner(label).op
        with mp.workdps(self.c.dps+40):
            rows=previous.previous.affine_transport(op.flow,incoming,packet['actual_whole_patch_local_defect_integral_C0_Z'],
                                                  packet['full_original_logarithmic_interval_width'])
        return dict(source_family=self.family,source_frame=label,exact_common_P0_axial5=op.P0,
            conditional_whole_patch_transported_defect_C0_Z=rows,
            supplied_Rm_incoming_enclosure_not_complete_upstream_prefix_proof=True,
            **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalRmAllUDensityIntegrals(require_checked=False);frames={}
    for label in ('0','.5'):
        frames[label]=owner.contribution(label)
        print('Original all-u Rm complete local integral bounds',label,
              len(frames[label]['actual_all_u_closed_edge_sources']),'edge cells enclosed',flush=True)
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_all_u_theorem=THEOREM,frames=fields.serialized(frames),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),
                                      compresslevel=6,mtime=0));return report


if __name__=='__main__':run()
