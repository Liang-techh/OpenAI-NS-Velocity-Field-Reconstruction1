"""Current whole-Z Rh function join, with a correction-preserving continuation.

The leading equality follows from the current common-box implicit map and
full bump weights. Correction histories remain separate Duhamel functions;
their directed ranges are never used as the functions' defining values.
"""
import ast
from fractions import Fraction
import hashlib
import inspect
import json
from pathlib import Path
import textwrap
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_whole_Z_all_N_outer_Rc_functions as outer
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_patch_mixed4 as patch
import lei_ren_part1_paper_compliant_current_original_Rm_patch_mixed4_cells as mixed
import lei_ren_part1_paper_compliant_current_original_Rh_reference_finite_N as reference
import lei_ren_part1_paper_interval_five_bump_inverse as inverse

HERE,PREFIX,sha,bind=outer.HERE,outer.PREFIX,outer.sha,outer.bind
CELLS,RATES,OPEN=outer.CELLS,outer.RATES,outer.OPEN
Pair,add,scale=outer.Pair,outer.add,outer.scale
NAME=PREFIX+'current_original_whole_Z_Rh_functional_join.json'
RECEIPT=PREFIX+'current_original_whole_Z_Rh_functional_join_check.json'
GATES=('current_whole_Z_repaired_leading_Rh_reference_function_identity_checked',
       'current_whole_Z_Rh_C1_correction_preserving_reference_continuation_installed')


def serialized(value):
    return outer.core.current.encode(outer.serialized(value))


def ast_node(fn):
    return ast.parse(textwrap.dedent(inspect.getsource(fn))).body[0]


def require_statements(fn,texts):
    node=ast_node(fn);all_nodes=list(ast.walk(node))
    for text in texts:
        expected=ast.dump(ast.parse(text).body[0])
        if sum(ast.dump(n)==expected for n in all_nodes)!=1:
            raise ValueError('Actual defining equation changed: '+fn.__qualname__+' '+text)
    return dict(module=Path(inspect.getfile(fn)).name,
        module_sha256=sha(Path(inspect.getfile(fn)).name),
        function=fn.__qualname__,AST_sha256=hashlib.sha256(ast.dump(node).encode()).hexdigest(),
        required_defining_statements=texts)


def source_bindings():
    """Assert live formulas, not just a family hash or an old join receipt."""
    reference_binding=reference.source_bindings()
    # This legacy metadata is a different formal rebase. Bind the actual
    # current whole-Z phase below; do not transplant the old phase recipe.
    reference_binding.pop('exact_reference_phase')
    changes="""changes=[f.add(*[f.scale(h[i],w(k,'0')) for i,k in enumerate((0,2))]),
        f.add(*[f.add(f.scale(h[i],w(k,'.6')),f.scale(f.multiply(h[i],h[k+2]),w(k,'.5',2))) for i,k in enumerate((0,2))]),
        f.add(*[f.scale(h[k+2],w(k,'.5')) for k in range(3)]),
        f.add(*[f.scale(f.multiply(f.multiply(h[i],h[i]),op.invAm2),w(k,'0',2)) for i,k in enumerate((0,2))],
            *[f.scale(h[k+2],-w(k,'.1')) for k in range(3)],
            *[f.scale(f.multiply(h[k+2],h[k+2]),-w(k,'0',2)/2) for k in range(3)]),
        f.add(*[f.add(f.scale(h[k+2],w(k,'-.9')),f.scale(f.multiply(h[k+2],h[k+2]),w(k,'-1',2)/2)) for k in range(3)])]"""
    return dict(partial_primitives=require_statements(mixed.partial_initial,[changes,
            'diagnostic=[f.add(a,b) for a,b in zip(op.defects,changes)]',
            "theta=f.add(defects[2],unit(c.mpf(5)/8*x**c.mpf('1.6')))",
            "initial=dict(mass=f.add(f.scale(op.zrows,4),f.scale(defects[0],1/x)),theta=theta,mixed=f.add(f.scale(f.multiply(op.zrows,theta),4),defects[1]),energy=f.add(defects[3],unit(-c.mpf(5)/12*x**c.mpf('1.2'))),pressure=f.add(defects[4],unit(c.mpf(5)/2*x**c.mpf('.2'))))"]),
        current_inverse=require_statements(mixed.previous._ActualRmLeadingPatch.__init__,[
            "self.inverse=implicit_axial_jets(c,self.Wn,self.normalized_defects,self.normalized_invAm2,[c.mpf(1)]*5)",
            "self.controls=[f.scale(f.jet(row),self.units[i]) for i,row in enumerate(self.inverse['controls'])]",
            "self.invAm2=f.scale(f.jet(self.poly),f.factor((0,-1,0,0,0)))",
            "self.angularUnit=f.factor((0,-1,0,0,0))*(10**10*self.t*self.t)",
            "self.normalized_invAm2=self.poly/(10**10)"]),
        transformed_weights=require_statements(mixed.previous.normalized_weights,[
            "return {name:rows if name in ('L','gg') else [v*angular_cover for v in rows] for name,rows in W.items()}"]),
        fixed_matrix_from_same_full_weights=require_statements(inverse.weights,[
            'L[0][0]=L[0][1]=c.mpf(1)',
            "for j,i in enumerate((0,2)):L[1][j]=w(i,'.6')",
            "return dict(L=L,fg=[w(0,'.5',2),w(2,'.5',2)],gg=[w(0,'0',2),w(2,'0',2)],ff=[w(i,'0',2) for i in range(3)],ff_over_x=[w(i,'-1',2) for i in range(3)])"]),
        actual_quadratic_map=require_statements(inverse.nonlinear,[
            "return [c.mpf(0),W['fg'][0]*c1*x1+W['fg'][1]*c2*x3,c.mpf(0),invAm2*sum((W['gg'][i]*cs[i]*cs[i] for i in range(2)),c.mpf(0))-sum((W['ff'][i]*xs[i]*xs[i] for i in range(3)),c.mpf(0))/2,sum((W['ff_over_x'][i]*xs[i]*xs[i] for i in range(3)),c.mpf(0))/2]"]),
        current_patch_derivatives=require_statements(mixed.mixed_functions,[
            'G=[f.scale(initial[\'mass\'],x)]+V[:4]',
            "raw_hist_x=dict(m=mass,h=[f.multiply(op.amrows,scalar_product(f,theta,angular,k)) for k in range(5)],k=[f.multiply(op.amrows,scalar_product(f,mixed,angular,k)) for k in range(5)],e=[scalar_product(f,raw_energy,inverse,k) for k in range(5)],p=[f.multiply(am2,row) for row in pressure])",
            "P=[f.add(axis,inc[0])]+inc[1:]",
            "Q=[previous.previous.previous.radial_Q(f,op.reference.Z,op.reference.delta,V[k],mass[k]) for k in range(5)]"]),
        reference_functions=reference_binding,
        reference_background=require_statements(reference.background_cell,[
            'logE=-reference.background.logarithm(1+ref.z*ref.z)+x/10',
            'radius=op.Rm_factor*c.exp(x+6)',
            "shapes=dict(theta=unit(c.mpf(5)/8),theta_z=f.scale(ref.zrows,c.mpf(5)/2),mean=ref.Vref,axial=f.scale(f.multiply(ref.zrows,ref.zrows),16),swirl=unit(c.mpf(5)/6),pressure=unit(5))"]),
        full_weights=require_statements(mixed.previous.SharedFiveMomentRepair.partial_weight,[
            'if power==0 and multiplicity==1:return c.mpf(1)',
            'return self.full_weights[center,power,multiplicity]']),
        current_all_N_reference_frontend=outer.FRONTEND_BINDINGS['rh'],
        current_all_N_density=outer.core.density_binding(),
        current_all_N_patch_transport=require_statements(outer.previous.WholeZAllNLongPatchFunctions.transport,[
            "incoming=self.normalized_R110_incoming(ends)",
            "contribution=scale(query['N_scaled_density_C0_Z'][name],weight['original_mass']*weight['suffix_decay'])",
            'outgoing={name:add(scale(incoming[name],memory[name]),local[name]) for name in RATES}']),
        current_all_N_outer_transport=require_statements(outer.WholeZAllNOuterRcFunctions.transport,[
            'incoming=self.normalized_Rh_incoming(ends)',
            'initial=dict(incoming)',
            'outgoing={name:add(scale(incoming[name],memory[name]),local[name]) for name in RATES}']))


def exact_join_proof():
    """Independent full-weight map algebra and common physical functions."""
    h=s.Matrix(s.symbols('c1 c2 xi1 xi2 xi3'));iv=s.Symbol('invAm2',positive=True)
    w=lambda k,p,m=1:s.Symbol('w_'+str(k)+'_'+str(p).replace('-','n').replace('.','p')+'_'+str(m))
    L=s.zeros(5);L[0,0]=L[0,1]=1
    for i,k in enumerate((0,2)):L[1,i]=w(k,'.6')
    for k in range(3):L[2,k+2]=w(k,'.5');L[3,k+2]=-w(k,'.1');L[4,k+2]=w(k,'-.9')
    Q=s.Matrix([0,sum(h[i]*h[k+2]*w(k,'.5',2) for i,k in enumerate((0,2))),0,
        iv*sum(h[i]**2*w(k,'0',2) for i,k in enumerate((0,2)))-sum(h[k+2]**2*w(k,'0',2) for k in range(3))/2,
        sum(h[k+2]**2*w(k,'-1',2) for k in range(3))/2])
    partial=s.Matrix([h[0]+h[1],sum(h[i]*w(k,'.6')+h[i]*h[k+2]*w(k,'.5',2) for i,k in enumerate((0,2))),
        sum(h[k+2]*w(k,'.5') for k in range(3)),
        sum(iv*h[i]**2*w(k,'0',2) for i,k in enumerate((0,2)))-sum(h[k+2]*w(k,'.1')+h[k+2]**2*w(k,'0',2)/2 for k in range(3)),
        sum(h[k+2]*w(k,'-.9')+h[k+2]**2*w(k,'-1',2)/2 for k in range(3))])
    assert all(s.expand(v)==0 for v in partial-L*h-Q)
    t,a=s.symbols('t angularUnit',positive=True);S=s.diag(t,t,a,a,a)
    assert all(s.expand(v)==0 for v in L*S-S*L)
    u=s.Matrix(s.symbols('u0:5'))
    # gg and ff share physical bump integrals but transform differently:
    # gg belongs to the axial square, ff to the angular square.
    Qn=s.Matrix([0,a*sum(u[i]*u[k+2]*w(k,'.5',2) for i,k in enumerate((0,2))),0,
        iv*t*t/a*sum(u[i]**2*w(k,'0',2) for i,k in enumerate((0,2)))-a*sum(u[k+2]**2*w(k,'0',2) for k in range(3))/2,
        a*sum(u[k+2]**2*w(k,'-1',2) for k in range(3))/2])
    assert all(s.expand(v)==0 for v in S.inv()*Q.subs(dict(zip(h,S*u)))-Qn)
    y,z=s.symbols('y Z',real=True);P=s.Symbol('Pstar',positive=True);delta=s.Symbol('delta',real=True)
    P0=s.Function('P0')(z);x=s.exp(y);am=s.exp(-s.Rational(3,5))/(1+z*z)
    E=am*s.exp(y/10);V=4*z;theta=s.Rational(5,8)*s.exp(8*y/5)
    mixed0=V*theta;energy=-s.Rational(5,12)*s.exp(6*y/5);pressure=s.Rational(5,2)*s.exp(y/5)
    radial=(2*z*V-(1-delta)*z*V-(1-z*z)*s.diff(V,z))/(1-delta*z*z)
    left=dict(E=E,V=V,m=V,h=am*theta*s.exp(-3*y/2),k=am*mixed0*s.exp(-3*y/2),
        e=V*V/P**2+am*am*energy/x,p=am*am*pressure,P0=P0,Q=radial,
        Utheta=E,Uz=V,Ur=s.exp(y/2)*radial,absolute_pressure=P0+am*am*pressure,
        Mz=x*V,Mtheta=am*theta,Mtheta_z=am*mixed0,
        Mztheta=x*V*V/P**2+am*am*energy,Mp=am*am*pressure)
    er=s.exp((y-6)/10)/(1+z*z);hr=s.Rational(5,8)*er
    right=dict(E=er,V=V,m=V,h=hr,k=hr*V,e=V*V/P**2-s.Rational(5,12)*er*er,
        p=s.Rational(5,2)*er*er,P0=P0,Q=radial,Utheta=er,Uz=V,Ur=s.exp(y/2)*radial,
        absolute_pressure=P0+s.Rational(5,2)*er*er,Mz=x*V,
        Mtheta=s.exp(3*y/2)*hr,Mtheta_z=s.exp(3*y/2)*hr*V,
        Mztheta=x*(V*V/P**2-s.Rational(5,12)*er*er),Mp=s.Rational(5,2)*er*er)
    assert all(s.simplify(s.expand_power_exp(left[k]-right[k]))==0 for k in left)
    return dict(passed=True,full_support_primitive_changes_equal_current_Lh_plus_Q=True,
        exact_positive_block_normalization_commutes=True,original_defect_equation='d+Lh+Q(h,Am^-2)=0',
        full_weight_identity_uses_unique_current_smooth_leading_family=True,
        exact_common_functions=list(left),physical_mixed4_functions=['Utheta','Uz','Ur','absolute_pressure','Mz','Mtheta','Mtheta_z','Mztheta','Mp'],
        implied_ordinary_y_Z_mixed4_rows=135,ordinary_Z_derivatives_of_equal_functions_retained=True,
        coordinate_identity='x=exp(y), reference_offset=y-6, Rh:y=1, offset=-5',
        radius_identity='Rm*exp(1)=Rref*exp(-5), Rref=Rm*exp(6)',
        leading_terminal_neighborhood='x>71/40, with all actual bumps zero',
        common_physical_units_attached_after_differentiation=True,
        independent_P0_is_not_part_of_five_defect_map=True,
        interval_overlap_is_not_the_functional_identity_proof=True)


class WholeZRhFunctionalJoin:
    """One current all-N graph, pure leading seam and C1 correction memory."""
    def __init__(self,dps=500,owner=None,require_checked=True):
        self.outer=owner if owner is not None else outer.WholeZAllNOuterRcFunctions(dps)
        self.c=self.outer.c;self.identity=self.outer.identity;self.N0=self.outer.N0
        self.hashes=dict(self.outer.hashes);self.cache={};self.bindings=source_bindings();self.proof=exact_join_proof()
        self.leading=self.outer.owner.patch_source.source
        if type(self.leading) is not patch.WholeZRmPatchMixed4:raise ValueError('Current whole-Z leading provider required')
        for name in (patch.RECEIPT,patch.upstream.RECEIPT):
            raw=json.loads((HERE/name).read_bytes())
            if not raw.get('all_passed') or raw['source_family']!=self.identity:raise ValueError('Current leading admission required')
            for path,digest in raw['input_hashes'].items():bind(self.hashes,path,digest)
            bind(self.hashes,name,sha(name))
        admission=json.loads((HERE/patch.upstream.RECEIPT).read_bytes())
        if not admission['same_common_unit_box_proves_unique_whole_Z_smooth_leading_family']:
            raise ValueError('A current whole-axis uniqueness theorem is required')
        for module in (mixed,mixed.previous,inverse,reference,patch,outer):
            name=Path(module.__file__).name;bind(self.hashes,name,sha(name))
        bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        self.acceptance_loaded=False
        if require_checked:
            raw=json.loads((HERE/RECEIPT).read_bytes())
            if not raw.get('all_passed') or not all(raw.get(k) for k in GATES) or raw['source_family']!=self.identity:
                raise ValueError('Dedicated current Rh functional join receipt required')
            for path,digest in raw['input_hashes'].items():bind(self.hashes,path,digest)
            self.acceptance_loaded=True

    def seam(self,ends):
        key=tuple(ends)
        if key not in CELLS:raise ValueError('Original whole-Z cell required')
        if key in self.cache:return self.cache[key]
        op=self.outer.owner.owner(key);patch_op=self.leading.owner(key).op
        if patch_op is not op:raise ValueError('One live current implicit leading owner required')
        packet=self.leading.query(key,'Rh');background=self.outer.leading_packet(key,'Rh_reference',(-5,1))
        if packet['original_P0_normalized_axial5'] is not op.P0 or background['exact_common_P0_axial5'] is not op.P0:
            raise ValueError('One independent P0 object must pass both branches')
        terminal=packet['actual_partial_primitive_source_memory']
        if not terminal['exact_terminal_identity_of_same_leading_map']:raise ValueError('Actual complete-support source required')
        if packet['geometry'].get('exact_x_source')!='exp(1)' or outer.ep(self.c.exp(1))[0]<=mp.mpf(71)/40:
            raise ValueError('This function identity is only for the actual complete-support Rh seam')
        inverse_record=op.inverse['initial_C1_inverse']
        if not inverse_record['certified'] or not inverse_record['unique_in_initial_box']:
            raise ValueError('Actual unique current implicit coefficients required')
        incoming=self.outer.normalized_Rh_incoming(key)
        result=dict(source_identity=self.identity,exact_Z_cell=list(key),fixed_leading_input_N0=self.N0,
            current_unique_implicit_owner_shared_by_both_branches=True,exact_common_P0_axial5=op.P0,
            leading_terminal_source_function_proof=self.proof,
            exact_Rm_radius=op.Rm_factor,exact_Rh_radius=op.Rm_factor*self.c.exp(1),
            actual_leading_patch_raw_histories=packet['raw_current_radius_y_derivative_axial_coefficients'],
            actual_reference_histories=background['original_generic_source']['common_own_five_histories_axial5'],
            genuine_N_scaled_correction_incoming_ranges={k:outer.record(v) for k,v in incoming.items()},
            incoming_defining_source=dict(module=Path(outer.previous.__file__).name,
                method='WholeZAllNLongPatchFunctions.transport',terminal='actual_uniform_N_scaled_Rh_correction_C0_Z',
                original_core_zero_inlet_and_bridge_R110_predecessor_chain_retained=True,
                source_report=outer.previous.NAME,source_receipt=outer.previous.RECEIPT),
            current_global_phase='frac(N*(logRm+y-logRa-hb*s_c/2))',
            reference_global_phase='frac(N*(logRm+6+offset-logRa-hb*s_c/2))',
            phase_same_at_Rh_and_not_restarted=True,phase_Z_exact_zero=True,
            leading_identity_does_not_zero_correction_incoming=True,
            incoming_range_is_not_a_selected_function_value=True,
            source_function_identity_from_original_full_weight_unique_leading_map=True,
            shared_original_P0=True,overlap_only_consistency_not_functional_identity_proof=False)
        self.cache[key]=result;return result

    def continuation(self,ends,offset=(-5,1)):
        """Actual reference source recipe and directed C1 history ranges.

        For each j, H_j(y)=exp(-lambda_j*y) H_j(Rh)+
        integral_0^y exp(-lambda_j*(y-s)) D_j(s) ds. H=N*correction.
        The same inlet handle is used at y=0, including its ordinary Z row.
        """
        q=reference.fraction(offset);key=tuple(ends);join=self.seam(key)
        op=self.outer.owner.owner(key);f=op.flow;c=self.c
        y=c.mpf((q+5).numerator)/(q+5).denominator
        incoming=self.outer.normalized_Rh_incoming(key)
        background=self.outer.leading_packet(key,'Rh_reference',offset)['original_generic_source']['common_own_five_histories_axial5']
        density=None if q==-5 else self.outer.query(key,'Rh_reference',(-5,1),offset)['N_scaled_density_C0_Z']
        corrections={};complete={};recipes={};memories={};local={}
        if self.outer.owner.patch_source.dstar_log._mpi_!=self.outer.owner.source.dstar_log._mpi_:
            raise ValueError('The actual patch/reference inverse source parameter must agree')
        eps=outer.core.epsilon_cover(c,self.N0)
        for name,rate in RATES.items():
            if q==-5:
                memory=f.scalar(1);driver=Pair(f.scalar(0),f.scalar(0));value=incoming[name]
            else:
                weight=self.outer.weights(key,'Rh_reference',(-5,1),offset,rate)
                memory=weight['cell_decay'];driver=scale(density[name],weight['original_mass'])
                value=add(scale(incoming[name],memory),driver)
            memories[name]=memory;local[name]=outer.record(driver);corrections[name]=outer.record(value)
            complete[name]=outer.record(add(Pair(*background[name][:2]),scale(value,eps)))
            inlet_handle=dict(operation='actual_current_N_scaled_Rh_correction_function',name=name,
                source_family=self.identity,exact_Z_cell=list(key),source=join['incoming_defining_source'])
            recipes[name]=dict(operation='own_rate_Duhamel_continuation',inlet_handle=inlet_handle,
                exact_log_distance=[(q+5).numerator,(q+5).denominator],own_rate=str(rate),
                kernel='exp(-lambda*(y-s))',physical_measure='ds=d(logR)',
                source_density=dict(module=Path(outer.core.__file__).name,method='normalized_drivers',
                    actual_phase='frac(N*(logRh+s-logRa-hb*s_c/2))',
                    original_current_all_N_primitive_provider='WholeZAllNOuterRcFunctions.primitive_cell'),
                at_Rh_exactly_same_inlet_handle=q==-5,
                exact_Z_derivative_uses_same_inlet_Z_and_density_Z=True)
        return dict(source_identity=self.identity,exact_Z_cell=list(key),exact_reference_offset=list(offset),
            fixed_leading_input_N0=self.N0,epsilon=eps,exact_common_P0_axial5=op.P0,
            defining_correction_function_graph=recipes,N_scaled_correction_range_C0_Z=corrections,
            N_scaled_local_integral_range_C0_Z=local,own_rate_incoming_memory=memories,
            complete_background_plus_epsilon_correction_range_C0_Z=complete,
            complete_state_definition='canonical_leading_background+epsilon*H, epsilon=1/N',
            ranges_are_enclosures_not_defining_function_values=True,
            pressure_rate_zero_memory_is_one=True,original_P0_independent_and_included_once=True,
            full_finite_N_mixed4_or_global_field_installed=False,**dict.fromkeys(OPEN,False))


def run(owner=None):
    began=time.monotonic()
    with mp.workdps(540):
        field=WholeZRhFunctionalJoin(owner=owner,require_checked=False)
        cells=[]
        for ends in CELLS:
            cells.append(dict(join=field.seam(ends),reference_inlet=field.continuation(ends),
                reference_first_unit=field.continuation(ends,(-4,1))))
            print('Current functional Rh join and retained correction: '+str(ends),flush=True)
        report=dict(all_passed=True,candidate_current_Rh_functional_join_constructed=True,
            source_family=field.identity,fixed_leading_input_N0=field.N0,actual_source_bindings=field.bindings,
            exact_current_function_proof=field.proof,source_cells=cells,
            **dict.fromkeys(GATES+OPEN,False),input_hashes=field.hashes,execution_seconds=time.monotonic()-began)
        (HERE/NAME).write_text(json.dumps(serialized(report),indent=2)+'\n',encoding='utf8')
    return field


if __name__=='__main__':run()
