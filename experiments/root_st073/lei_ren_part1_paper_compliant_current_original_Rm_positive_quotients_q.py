"""Actual leading Rm positive source quotients and original shear q.

Whole radial source cells at the admitted axial frames retain all five
bases. Inertial quotients are differentiated before their one shared R
factor. Original eta/dstar are parameter definitions, not transferred cone
or phase admission theorems for the new owner.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rm_generic_inputs as previous
import lei_ren_part1_paper_compliant_current_native_correlated_shear_q as original

fields,base,ep=previous.fields,previous.base,previous.ep
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_Rm_positive_quotients_q.json.gz'
RECEIPT=PREFIX+'current_original_Rm_positive_quotients_q_check.json'
GATE='original_actual_two_frame_whole_radial_Rm_positive_quotients_and_q_installed'
UNIFORM=PREFIX+'current_generic_shear_uniform_inputs.json'
UNIFORM_CHECK=PREFIX+'current_generic_shear_uniform_inputs_check.json'


def add(f,*rows):
    return [sum((row[n] for row in rows),f.scalar(0)) for n in range(min(map(len,rows)))]


def multiply(f,a,b):
    return [sum((a[j]*b[n-j] for j in range(n+1)),f.scalar(0)) for n in range(min(len(a),len(b)))]


def scale(rows,value):return [v*value for v in rows]


def same_source(f,rows):
    if any(row.ctx is not f.c or row.scale.bases is not f.logs or row.ledger is not f.ledger for row in rows):
        raise ValueError('Same live context, canonical source bases and ledger required')


def source_coordinate(c,geometry):
    if geometry['point']:
        return c.exp(1) if geometry.get('exact_x_source')=='exp(1)' else c.mpf(geometry['exact_x'][0])/geometry['exact_x'][1]
    left=c.mpf(geometry['exact_left'][0])/geometry['exact_left'][1]
    right=c.exp(1) if geometry.get('exact_right_source')=='exp(1)' else c.mpf(geometry['exact_right'][0])/geometry['exact_right'][1]
    return c.mpf([ep(left)[0],ep(right)[1]])


def zero_mode_sqrt(value):
    if value.scale.powers!=(0,0,0,0,0):
        raise ValueError('Square root needs verified source rebase for nontrivial factor powers')
    return original.nonnegative_sqrt(value)


def positive_source(f,row,role):
    """A complete source enclosure, not an unsigned sector or selected value."""
    same_source(f,[row])
    lo,hi=ep(row.coefficient)
    if lo<=0 or not mp.isfinite(hi):raise ValueError('Strict positive whole-source enclosure required: '+role)
    lower=ep(row.scale.evaluate()+f.c.ln(f.c.mpf(lo)))[0]
    upper=ep(row.scale.evaluate()+f.c.ln(f.c.mpf(hi)))[1]
    return dict(role=role,source_row=row,source_log_lower=f.c.mpf(lower),source_log_upper=f.c.mpf(upper),
        strict_positive_complete_source_enclosure=True,source_value_not_selected=True)


def quotient(f,numerator,denominator,proof):
    if not denominator or proof['source_row'] is not denominator[0]:
        raise ValueError('Positive theorem must bind the exact same denominator source')
    count=min(len(numerator),len(denominator));result=[]
    for n in range(count):
        correction=sum((denominator[j]*result[n-j] for j in range(1,n+1)),f.scalar(0))
        result.append((numerator[n]-correction).positive_divide(denominator[0],proof['source_log_lower']))
    return result


def correlated_shear_a(op,patch):
    """Cancel am analytically and differentiate only the true bump change.

    E=am*(x^.1+f), so a=.8-2*(x*f_x-.1*f)/(x^.1+f).
    Forming f/f_x from actual controls avoids subtracting the reference
    profile from an interval enclosure of itself.
    """
    f=op.flow;c=op.c;x=source_coordinate(c,patch['geometry'])
    gamma=patch['actual_gamma_ordinary_x_derivatives']
    correction=add(f,*[scale(op.controls[j+2],gamma[j][0]) for j in range(3)])
    derivative=add(f,*[scale(op.controls[j+2],gamma[j][1]) for j in range(3)])
    numerator=add(f,scale(derivative,x),scale(correction,-c.mpf(1)/10))
    H=patch['actual_H_x_derivative_axial5'][0]
    proof=positive_source(f,H[0],'actual_H_without_am')
    change=scale(quotient(f,numerator,H,proof),-2)
    a=add(f,[f.scalar(c.mpf(4)/5)]+[f.scalar(0)]*5,change)
    return a,dict(original_exact_identity='a=4/5-2*(x*f_x-f/10)/(x^.1+f)',
        common_axial_amplitude_canceled_analytically=True,
        angular_correction_numerator=numerator,actual_shear_correction_axial5=change,
        exact_reference_shear=f.c.mpf(4)/5,positive_H_source_theorem=proof,
        reference_profile_not_subtracted_from_its_own_enclosure=True)


def shear_q_jets(f,a,Delta,eta_log,proof_a):
    """Original lazy C0 cutoff, with analytic jets on certified sigma=1/flat."""
    loop=original.q_enclosure(a[0],Delta[0],eta_log,proof_a['source_log_lower'])
    count=min(len(a),len(Delta));zero=[f.scalar(0)]*count
    if loop['branch']=='flat':
        return dict(original_C0_lazy_cutoff=loop,q_axial_coefficients=zero,q_squared_axial_coefficients=zero,
            analytic_axial_jets_available=True,jet_branch='exact original flat cutoff',q_positive_theorem=None)
    if ep(Delta[0].coefficient)[1]>0:
        return dict(original_C0_lazy_cutoff=loop,q_axial_coefficients=None,q_squared_axial_coefficients=None,
            analytic_axial_jets_available=False,jet_branch='cutoff source cell requires refinement',q_positive_theorem=None)
    eta=f.factor((0,0,0,0,0),eta_log)
    gamma=add(f,[eta*2]+zero[1:],scale(Delta,-1))
    twice_a=scale(a,2);denproof=positive_source(f,twice_a[0],'2a')
    q2=quotient(f,gamma,twice_a,denproof)
    positive_source(f,q2[0],'q_squared')
    q0=zero_mode_sqrt(q2[0]);q=[q0]
    twoq=[q0*2];qproof=positive_source(f,twoq[0],'2q')
    for n in range(1,count):
        correction=sum((q[j]*q[n-j] for j in range(1,n)),f.scalar(0))
        q.append((q2[n]-correction).positive_divide(twoq[0],qproof['source_log_lower']))
    return dict(original_C0_lazy_cutoff=loop,q_axial_coefficients=q,q_squared_axial_coefficients=q2,
        analytic_axial_jets_available=True,jet_branch='Delta<=0: original sigma identically1',
        q_positive_theorem=positive_source(f,q0,'q'),
        q0_square_root_has_no_nontrivial_source_powers_to_collapse=True)


def recover_quotients(op,source,eta_log,dstar_log,*,patch=None):
    f=op.flow;c=op.c;n=source['actual_generic_source_numerators']
    if c is not f.c or source['original_fixed_source_log_bases'] is not f.logs or not source['source_ledger_is_same_object']:
        raise ValueError('Same actual Rm source context, basis and ledger binding required')
    if source['actual_patch_source_owner']!=type(op).__name__:
        raise ValueError('Actual patch source owner binding required')
    if source['common_original_P0_axial5'] is not op.P0:
        raise ValueError('Exact same actual Rm P0 source required')
    for rows in (*n.values(),*source['full_signed_inertial_sectors_axial4'].values(),op.P0):same_source(f,rows)
    R=source['original_physical_radius'];radius_proof=positive_source(f,R,'actual_Rm_times_x')
    same_source(f,[op.Rm_factor])
    expected_R=op.Rm_factor*source_coordinate(c,source['source_geometry'])
    if R.scale.powers!=expected_R.scale.powers or R.scale.offset._mpi_!=expected_R.scale.offset._mpi_ or R.coefficient._mpi_!=expected_R.coefficient._mpi_:
        raise ValueError('Radius must be exactly the same factored actual Rm*x source')
    proofs={key:positive_source(f,n[key][0],key) for key in ('E','C','positive_denominator')}
    if patch is None:
        a=quotient(f,n['C'],n['E'],proofs['E']);shear=dict(common_axial_amplitude_canceled_analytically=False)
    else:
        if source['source_geometry']!=patch['geometry']:
            raise ValueError('Generic source and correlated patch must have the same exact geometry')
        if patch['original_P0_normalized_axial5'] is not op.P0:
            raise ValueError('Correlated H must belong to the same live P0/Rm owner')
        for rows in (*patch['actual_H_x_derivative_axial5'],*op.controls):same_source(f,rows)
        a,shear=correlated_shear_a(op,patch)
    b=quotient(f,n['B'],n['E'],proofs['E'])
    proofs['a']=positive_source(f,a[0],'a')
    t0=quotient(f,scale(n['B'],-1),n['C'],proofs['C'])
    # Square the same C0 source, never two independent signed boxes.
    b2=multiply(f,b,b);b2[0]=original.square(b[0])
    kappa=add(f,a,quotient(f,b2,a,proofs['a']))
    Delta=add(f,kappa,[f.scalar(-2)]+[f.scalar(0)]*(len(kappa)-1))
    sectors=source['full_signed_inertial_sectors_axial4']
    before_R=dict(p1=add(f,sectors['theta_linear'],scale(sectors['theta_quadratic'],op.Pstar)),
        p2=add(f,sectors['axial_linear'],scale(sectors['axial_quadratic'],op.Pstar)))
    inertial={key:quotient(f,row,n['E'],proofs['E']) for key,row in before_R.items()}
    full_inertial={key:scale(row,R) for key,row in inertial.items()}
    H0_over_R=add(f,inertial['p1'],multiply(f,inertial['p2'],t0))
    J_over_R=add(f,inertial['p2'],scale(multiply(f,inertial['p1'],t0),-1))
    loop=shear_q_jets(f,a,Delta,eta_log,proofs['a'])
    q=loop['original_C0_lazy_cutoff']['q']
    dstar=f.factor((0,0,0,0,0),dstar_log)
    dproof=positive_source(f,dstar,'original_selected_dstar_parameter')
    u_over_R=(inertial['p2'][0]*q).positive_divide(dstar,dproof['source_log_lower'])
    u_axial=None
    if loop['analytic_axial_jets_available']:
        drows=[dstar]+[f.scalar(0)]*(len(inertial['p2'])-1)
        u_axial=quotient(f,multiply(f,inertial['p2'],loop['q_axial_coefficients']),drows,dproof)
    result=dict(source_geometry=source['source_geometry'],original_P0_is_same_live_object=True,
        common_original_P0_axial5=op.P0,source_positive_theorems=proofs,actual_correlated_shear_source=shear,
        shear_quotients_axial5=dict(a=a,b=b,t0=t0,kappa=kappa,kappa_minus2=Delta),
        full_signed_inertial_quotients_before_shared_R_axial4=inertial,
        actual_factored_p1_p2_axial4=full_inertial,
        H0_over_shared_R_axial4=H0_over_R,J_over_shared_R_axial4=J_over_R,
        exact_signed_radius_recipes=dict(H0='R*H0_over_R',D='R*H0_over_R-kappa',J='R*J_over_R'),
        original_shear_q=loop,original_eta_log_parameter=eta_log,original_dstar_log_parameter=dstar_log,
        selected_dstar_is_parameter_not_new_owner_margin_proof=True,
        original_signed_u_over_shared_R_C0=u_over_R,original_signed_u_C0=u_over_R*R,
        original_signed_u_before_shared_R_axial4=u_axial,
        exact_same_shared_positive_radius_factor=R,
        exact_same_shared_radius_positive_theorem=radius_proof,
        source_context_basis_ledger_geometry_and_P0_bound_to_same_actual_owner=True,
        inertial_quotients_differentiated_before_attaching_one_radius_factor=True,
        five_source_bases_and_all_directed_offsets_retained=True,
        source_enclosure_endpoints_used_only_as_proven_bounds_not_function_values=True,
        original_selected_parameters_do_not_transfer_old_owner_margin_theorems=True,
        actual_positive_E_C_on_requested_source_cell_certified=True,
        source_frame_conditional_on_same_accepted_Rm_inlet=True,
        whole_axis_provider_or_new_dstar_cone_margin_certified=False,
        actual_phase_inverse_or_finite_N_density_integrals_installed=False,
        **dict.fromkeys(fields.previous.OPEN,False))
    return result


class OriginalRmPositiveQuotientsQ:
    mode='actual_Rm_source_positive_quotients_and_original_lazy_q_with_shared_radius'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalRmGenericInputs(dps);self.c=self.upstream.c
        self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes)
        repair=self.upstream.upstream.upstream.repair;bind=fields.previous.bind
        rows=[]
        for name in (UNIFORM,UNIFORM_CHECK,original.RECEIPT):
            data=json.loads((HERE/name).read_bytes())
            if name!=UNIFORM and not data.get('all_passed'):raise ValueError('Accepted original scalar/q definitions required')
            if any(data['source_family'][key]!=repair[key] for key in data['source_family']):
                raise ValueError('Same actual family, implicit source and pressure datum required')
            for path,digest in data['input_hashes'].items():bind(self.hashes,path,digest)
            bind(self.hashes,name,sha(name));rows.append(data)
        if not rows[1].get('current_original_whole_generic_input_margins_and_log_scales_certified') or not rows[2].get(original.GATE):
            raise ValueError('Accepted original logarithmic scalar recipe and lazy q implementation required')
        scales=rows[0]['current_actual_logarithmic_loop_scales']
        self.eta_log=fields.previous.read_interval(self.c,scales['selected_positive_eta_log'])
        self.dstar_log=fields.previous.read_interval(self.c,scales['logarithmic_selected_positive_lower_constants']['d_star'])
        for value in (self.eta_log,self.dstar_log):
            if ep(value)[0]!=ep(value)[1]:raise ValueError('Exact original selected log parameter required')
        for name in (Path(__file__).name,Path(original.__file__).name):bind(self.hashes,name,sha(name))
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted actual Rm positive quotient/q receipt required')
            for name,digest in receipt['input_hashes'].items():bind(self.hashes,name,digest)
            bind(self.hashes,RECEIPT,sha(RECEIPT))
    def evaluate(self,label,coordinate):
        owner=self.upstream.upstream.owner(label);op=owner.op
        with mp.workdps(self.c.dps+40):
            return recover_quotients(op,self.upstream.evaluate(label,coordinate),self.eta_log,self.dstar_log,patch=owner.evaluate(coordinate))
    def cell(self,label,left,right):
        owner=self.upstream.upstream.owner(label);op=owner.op
        with mp.workdps(self.c.dps+40):
            return recover_quotients(op,self.upstream.cell(label,left,right),self.eta_log,self.dstar_log,patch=owner.cell(left,right))


def run():
    began=time.monotonic();owner=OriginalRmPositiveQuotientsQ(require_checked=False);frames={}
    for label in ('0','.5'):
        frames[label]=dict(whole_patch=owner.cell(label,(1,1),'Rh'),
            first_support=owner.cell(label,(6,5),(13,10)),
            terminal=owner.cell(label,(71,40),'Rh'),active_point=owner.evaluate(label,(5,4)))
        print('Actual Rm positive source quotients and original q',label,flush=True)
    result=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},
        frames=fields.serialized(frames),**dict.fromkeys(fields.previous.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(result),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return result


if __name__=='__main__':run()
