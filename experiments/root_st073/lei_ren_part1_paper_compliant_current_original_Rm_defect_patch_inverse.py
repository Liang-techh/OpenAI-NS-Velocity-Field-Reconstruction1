"""Actual two-frame Rm defects, normalized leading inverse and partial patch.

Coordinate normalization occurs in the common formal source basis BEFORE
ordinary interval materialization. No angular floor chooses a physical
coefficient. This is the leading Lh+Q+d map, not the finite-N Rc repair.
"""
import ast
from fractions import Fraction
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_reference_restore_functions as previous
import lei_ren_part1_paper_compliant_actual_moment_patch as original
from lei_ren_part1_paper_compliant_actual_moment_patch import implicit_axial_jets
from lei_ren_part1_paper_compliant_five_moment_repair import SharedFiveMomentRepair,pack
from lei_ren_part1_paper_interval_five_bump_inverse import weights
from lei_ren_part1_paper_interval_exit_continuation_enclosure import restore_value
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

fields,endpoint,base,ep=previous.fields,previous.endpoint,previous.base,previous.ep
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_Rm_defect_patch_inverse.json'
RECEIPT=PREFIX+'current_original_Rm_defect_patch_inverse_check.json'
GATE='original_actual_two_frame_Rm_leading_defect_inverse_and_partial_patch_installed'
ORDER=('c1','c2','xi1','xi2','xi3')


def source_bindings():
    tree=ast.parse(Path(previous.original.__file__).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='defects')
    text="[centered['mean_error']*c.exp(x),centered['mixed_error']*c.exp(c.mpf('1.6')*x),centered['angular_error']*c.exp(c.mpf('1.6')*x),(centered['axial_square']*invAm2)*c.exp(x)-centered['swirl_error']*(c.exp(c.mpf('1.2')*x)/2),centered['pressure_error']*(c.exp(c.mpf('.2')*x)/2)]"
    want=ast.dump(ast.parse(text,mode='eval').body)
    if sum(ast.dump(n.value)==want for n in ast.walk(fn) if isinstance(n,ast.Assign))!=1:
        raise ValueError('Original five normalized source defects changed')
    from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
    bound=assignment_source_bindings('actual_moment_patch','evaluate',{
        'change':"[sum((h[i]*w(k,'0') for i,k in enumerate((0,2))),zero),sum((h[i]*w(k,'.6')+h[i]*h[k+2]*w(k,'.5',2) for i,k in enumerate((0,2))),zero),sum((h[k+2]*w(k,'.5') for k in range(3)),zero),sum((square(h[i])*data['invAm2']*w(k,'0',2) for i,k in enumerate((0,2))),zero)-sum((h[k+2]*w(k,'.1')+square(h[k+2])*(w(k,'0',2)/2) for k in range(3)),zero),sum((h[k+2]*w(k,'-.9')+square(h[k+2])*(w(k,'-1',2)/2) for k in range(3)),zero)]",
        'mass':"z*4+defects[0]/x",'theta':"defects[2]+x**c.mpf('1.6')*c.mpf('.625')",
        'mixed':"(z*theta)*4+defects[1]",'energy':"defects[3]-x**c.mpf('1.2')*c.mpf(5)/12",
        'pressure_moment':"defects[4]+x**c.mpf('.2')*c.mpf('2.5')"})
    return dict(passed=True,original_defects_assignment=True,original_partial_recovery_assignments=bound,
        leading_map='Lh+Q(h,Am^-2)+d=0',finite_N_Rc_map_not_substituted=True,
        exact_units='S=diag(t,t,angularUnit,angularUnit,angularUnit); angularUnit=10^10*t^2*canonical_Pstar^-2',
        physical_units='Mtheta/Mtheta_z use fixed sqrt2*Rm^1.5*Am; centered Mztheta uses Rm*Am^2; Mp uses Am^2')


def normalized_weights(c,W,angular_cover):
    """Outer coefficient cover of S^-1 Q(Su), never a division unit.

    Its lower bound may be zero for an astronomically small positive source.
    Division uses the strict-positive FORMAL angularUnit, not this cover.
    Uniform contraction for the outer box includes the actual positive map.
    """
    if ep(angular_cover)[0]<0 or any(not mp.isfinite(v) for v in ep(angular_cover)):
        raise ValueError('Finite nonnegative directed cover of the positive angular unit required')
    return {name:rows if name in ('L','gg') else [v*angular_cover for v in rows]
            for name,rows in W.items()}


class _ActualRmLeadingPatch:
    def __init__(self,reference,t,W,bump):
        self.reference=reference;self.flow=f=reference.flow;self.c=c=f.c
        self.t=c.mpf(t);self.W=W;self.bump=bump
        if ep(self.t)[0]<=0 or ep(self.t)[0]!=ep(self.t)[1] or not mp.isfinite(ep(self.t)[1]):
            raise ValueError('Exact fixed positive scalar axial box coordinate required')
        if bump.ctx is not c:raise ValueError('Same-context original bump provider required')
        for i in range(5):
            for j in range(5):
                if (i<2)!=(j<2) and W['L'][i][j]!=0:
                    raise ValueError('Original commuting 2+3 block matrix required')
        self.Rm=reference.postrestore((-6,1));s=self.centered=self.Rm['actual_centered_histories']
        if self.Rm['geometry']['exact_reference_offset']!=[-6,1] or set(s)!=set(previous.RATES):
            raise ValueError('Actual Rm offset and original six centered row names required')
        self.P0=self.Rm['original_P0_normalized_axial5']
        if self.P0 is not reference.P0:raise ValueError('Same original analytic P0 object required')
        if any(len(row)!=6 or any(v.scale.bases is not f.logs or v.ledger is not f.ledger for v in row)
               for row in (*s.values(),self.P0)):
            raise ValueError('Same Rm basis/ledger and ordinary Z0..5 rows required')
        self.z=reference.z;self.zrows=f.jet(self.z)
        self.poly=(1+self.z*self.z)*(1+self.z*self.z)*c.exp(c.mpf('1.2'))
        self.invAm2=f.scale(f.jet(self.poly),f.factor((0,-1,0,0,0)))
        self.angularUnit=f.factor((0,-1,0,0,0))*(10**10*self.t*self.t)
        self.units=[f.scalar(self.t)]*2+[self.angularUnit]*3
        self.defects=[s['mean_error'],s['mixed_error'],s['angular_error'],
            f.add(f.multiply(s['axial_square'],self.invAm2),f.scale(s['swirl_error'],-c.mpf('.5'))),
            f.scale(s['pressure_error'],c.mpf('.5'))]
        self.normalized_defects=[IntervalTaylor(c,[f.ordinary_cover(v.positive_divide(self.units[i],
            self.units[i].scale.evaluate()+c.ln(self.units[i].coefficient))) for v in row])
            for i,row in enumerate(self.defects)]
        # ONLY the transformed weights need this directed nonzero cover.
        # Defects were divided in the exact source basis, cancelling Pstar^-2.
        self.angular_cover=f.ordinary_cover(self.angularUnit)
        self.Wn=normalized_weights(c,W,self.angular_cover)
        self.normalized_invAm2=self.poly/(10**10)
        self.inverse=implicit_axial_jets(c,self.Wn,self.normalized_defects,
            self.normalized_invAm2,[c.mpf(1)]*5)
        self.controls=[f.scale(f.jet(row),self.units[i]) for i,row in enumerate(self.inverse['controls'])]
        self.am=(1+self.z*self.z).reciprocal()*c.exp(c.mpf('-.6'))
        self.amrows=f.jet(self.am);self.Pstar=f.factor((0,.5,0,0,0));self.P2=f.factor((0,1,0,0,0))
        self.Rm_factor=f.factor((0,5,0,0,0),10*reference.logC+c.ln(110)-6)
        self.cache={}

    def coefficients(self):
        return dict(actual_Rm_centered_histories=self.centered,actual_five_defects=self.defects,
            normalized_five_defect_axial5=self.normalized_defects,canonical_invAm2=self.invAm2,
            original_fixed_axial_coordinate=self.t,canonical_positive_angular_unit=self.angularUnit,
            exact_row_and_control_units=self.units,directed_angular_weight_cover=self.angular_cover,
            normalized_invAm2=self.normalized_invAm2,original_weight_matrix=self.W,
            normalized_weight_enclosures=self.Wn,normalized_implicit_inverse_axial5=self.inverse,
            actual_factored_coefficient_axial5=self.controls,
            coefficients_are_enclosures_of_same_unique_implicit_functions=True,
            common_basis_division_before_materialization=True,canonical_Pstar_source_cancels_exactly=True,
            weight_cover_only_encloses_positive_source_and_is_never_used_as_divisor=True,
            scales_Z_independent_and_no_extra_product_terms=True,P0_excluded_from_five_defects=True,
            residual_zero_containment_only_diagnostic=True,whole_axis_implicit_family_certified=False,
            finite_N_Rc_targets_or_controls_installed=False)

    def evaluate(self,coordinate,partial_cells=64):
        q=previous.fraction(coordinate);c=self.c;f=self.flow
        if not 1<=q or q>Fraction(3):raise ValueError('Original exact rational patch coordinate required')
        x=c.mpf(q.numerator)/q.denominator
        if ep(x)[1]>ep(c.exp(1))[0]:raise ValueError('Original patch 1<=R/Rm<=e required')
        if type(partial_cells) is not int or partial_cells<16:raise ValueError('At least16 directed partial bump cells required')
        key=(q,partial_cells)
        if key in self.cache:return self.cache[key]
        zero=[f.scalar(0)]*6;h=self.controls
        unit=lambda v:[f.scalar(v)]+zero[1:]
        centers=(Fraction(5,4),Fraction(3,2),Fraction(7,4));beta=[]
        for cc in centers:beta.append(self.bump.beta(x-c.mpf(cc.numerator)/cc.denominator))
        fc=f.add(*[f.scale(h[k+2],beta[k][0]) for k in range(3)])
        fx=f.add(*[f.scale(h[k+2],beta[k][1]) for k in range(3)])
        gc=f.add(*[f.scale(h[i],beta[k][0]) for i,k in enumerate((0,2))])
        gx=f.add(*[f.scale(h[i],beta[k][1]) for i,k in enumerate((0,2))])
        ws={}
        def w(k,p,m=1):
            key=(k,p,m)
            if key not in ws:ws[key]=self.bump.partial_weight(x,centers[k],p,m,partial_cells)
            return ws[key]
        change=[f.add(*[f.scale(h[i],w(k,'0')) for i,k in enumerate((0,2))]),
            f.add(*[f.add(f.scale(h[i],w(k,'.6')),f.scale(f.multiply(h[i],h[k+2]),w(k,'.5',2))) for i,k in enumerate((0,2))]),
            f.add(*[f.scale(h[k+2],w(k,'.5')) for k in range(3)]),
            f.add(*[f.scale(f.multiply(f.multiply(h[i],h[i]),self.invAm2),w(k,'0',2)) for i,k in enumerate((0,2))],
                *[f.scale(h[k+2],-w(k,'.1')) for k in range(3)],
                *[f.scale(f.multiply(h[k+2],h[k+2]),-w(k,'0',2)/2) for k in range(3)]),
            f.add(*[f.add(f.scale(h[k+2],w(k,'-.9')),f.scale(f.multiply(h[k+2],h[k+2]),w(k,'-1',2)/2)) for k in range(3)])]
        diagnostic=[f.add(a,b) for a,b in zip(self.defects,change)];terminal=q>=Fraction(71,40)
        if terminal:
            for row in diagnostic:
                for v in row:
                    lo,hi=ep(v.coefficient)
                    if not lo<=0<=hi:raise ArithmeticError('Same full-weight implicit-map diagnostic excludes zero')
        # The exact implicit solution closes the same full-weight map.
        # Keep its unrefined enclosure diagnostic and every incoming row.
        defects=[zero]*5 if terminal else diagnostic
        mass=f.add(f.scale(self.zrows,4),f.scale(defects[0],1/x))
        theta=f.add(defects[2],unit(c.mpf(5)/8*x**c.mpf('1.6')))
        mixed=f.add(f.scale(f.multiply(self.zrows,theta),4),defects[1])
        energy=f.add(defects[3],unit(-c.mpf(5)/12*x**c.mpf('1.2')))
        pressure_moment=f.add(defects[4],unit(c.mpf(5)/2*x**c.mpf('.2')))
        H=f.add(unit(x**c.mpf('.1')),fc);V=f.add(f.scale(self.zrows,4),gc)
        if ep(H[0].coefficient)[0]<=0:raise ArithmeticError('Positive original patch swirl lost')
        Q=previous.previous.radial_Q(f,self.reference.Z,self.reference.delta,V,mass)
        Utheta=f.scale(f.multiply(self.amrows,H),self.Pstar)
        Utheta_y=f.scale(f.multiply(self.amrows,f.add(unit(x**c.mpf('.1')/10),f.scale(fx,x))),self.Pstar)
        R=self.Rm_factor*x;rootR=f.radial_power(R,.5)
        Am=f.scale(self.amrows,self.Pstar);Am2=f.scale(f.multiply(self.amrows,self.amrows),self.P2)
        Mz=f.scale(mass,R)
        centered_energy=f.scale(f.multiply(Am2,energy),self.Rm_factor)
        Mztheta=f.add(centered_energy,f.scale(f.multiply(self.zrows,Mz),8),
            f.scale(f.multiply(self.zrows,self.zrows),-16*R))
        Mp=f.multiply(Am2,pressure_moment);axis=f.scale(self.P0,self.P2);P=f.add(axis,Mp)
        primitive=dict(Mz=Mz,Mtheta=f.scale(f.multiply(Am,theta),self.Rm_factor*f.radial_power(self.Rm_factor,.5)*c.sqrt(2)),
            Mtheta_z=f.scale(f.multiply(Am,mixed),self.Rm_factor*f.radial_power(self.Rm_factor,.5)*c.sqrt(2)),Mztheta=Mztheta,Mp=Mp)
        physical_corrections=dict(Utheta=f.scale(f.multiply(self.amrows,fc),self.Pstar),Uz=gc,
            Mp=f.multiply(Am2,defects[4]),centered_Mztheta=f.scale(f.multiply(Am2,defects[3]),self.Rm_factor))
        value=dict(exact_x=[q.numerator,q.denominator],actual_Rm_source_memory=self.centered,
            signed_bump_correction_sectors=dict(f=fc,f_x=fx,g=gc,g_x=gx,five_partial_primitive_changes=change,
                physical_corrections=physical_corrections),actual_partial_weights={str(k):v for k,v in ws.items()},
            actual_five_defects_before_exact_terminal_identity=diagnostic,actual_recovered_five_defects=defects,
            normalized_primitives=dict(mean=mass,theta=theta,mixed=mixed,centered_energy=energy,pressure=pressure_moment),
            physical_velocity_axial_coefficients=dict(Ur=[v*rootR*(1/c.sqrt(2)) for v in Q],Utheta=Utheta,Uz=V),
            physical_velocity_first_y_axial5=dict(Utheta=Utheta_y,Uz=f.scale(gx,x)),
            physical_cumulative_moment_axial5=primitive,physical_centered_Mztheta_axial5=centered_energy,
            original_P0_normalized_axial5=self.P0,physical_pressure_axis_axial5=axis,
            physical_pressure_radial_increment_axial5=Mp,physical_total_pressure_axial5=P,
            formal_positive_geometry=dict(R=R,Rm=self.Rm_factor,Am=Am,sqrtR=rootR),
            exact_terminal_identity_of_same_unique_leading_map=terminal,
            terminal_identity_not_an_inlet_reset_or_zero_containment_proof=True,
            physical_prefactors_and_amplitude_Z_products_applied_before_export=True,
            same_original_P0_retained=True,source_caps_not_used_as_field_values=True,
            full_radial_mixed4_or_finite_N_density_oracle_installed=False,
            whole_axis_leading_patch_or_finite_N_Rc_patch_certified=False)
        self.cache[key]=value;return value


class OriginalRmDefectPatchInverse:
    mode='genuine_original_actual_Rm_leading_five_defect_inverse_and_partial_patch'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalReferenceRestoreFunctions(dps);self.c=c=self.upstream.c
        self.family=self.upstream.family;self.hashes=dict(self.upstream.hashes);self.owners={}
        bind=fields.previous.bind
        for name in (previous.NAME,previous.RECEIPT):
            row=json.loads((HERE/name).read_bytes())
            if not row.get(previous.GATE) or name==previous.RECEIPT and not row.get('all_passed'):
                raise ValueError('Accepted actual Rm inlet functions required')
            if row['source_family']!=self.family:raise ValueError('Same actual Rm family required')
            for path,digest in row['input_hashes'].items():bind(self.hashes,path,digest)
            bind(self.hashes,name,sha(name))
        self.repair=json.loads((HERE/(PREFIX+'five_moment_repair.json')).read_bytes())
        checked=json.loads((HERE/(PREFIX+'five_moment_repair_check.json')).read_bytes())
        admitted=json.loads((HERE/previous.ADMISSION).read_bytes())
        for key in ('actual_five_defect_family_sha256','implicit_source_sha256','datum_enclosure_sha256'):
            if self.repair[key]!=admitted[key]:raise ValueError('Same original leading map family/source/datum required')
        for gate in ('actual_implicit_functional_five_moment_closure_independently_checked',
            'corrected_partial_moments_and_same_pressure_independently_checked','residual_zero_containment_not_used_as_closure_proof'):
            if not checked.get(gate):raise ValueError('Original leading repair receipt required: '+gate)
        if checked['actual_five_defect_family_sha256']!=self.family:raise ValueError('Same original repair family required')
        for row in (self.repair,checked):
            for path,digest in row['input_hashes'].items():bind(self.hashes,path,digest)
        self.t=fields.previous.read_interval(c,self.repair['axial_box_scale'])
        self.admitted_logP=fields.previous.read_interval(c,admitted['logPstar'])
        wn='lei_ren_part1_paper_bump_integral_enclosures_check.json'
        self.wraw=json.loads((HERE/wn).read_bytes());self.W=weights(c,self.wraw)
        if not self.wraw.get('directed_integrals_certified'):raise ValueError('Original directed beta integrals required')
        self.fixed=json.loads((HERE/'lei_ren_part1_paper_shared_bump_constants.json').read_bytes())
        if any(ep(v)!=ep(fields.previous.read_interval(c,self.fixed['fixed_matrix'][i][j])) for i,row in enumerate(self.W['L']) for j,v in enumerate(row)):
            raise ValueError('Original fixed leading bump matrix required')
        for name,key in (('lei_ren_part1_paper_bump_integral_enclosures.py','module_sha256'),
            ('lei_ren_part1_paper_bump_integral_enclosures_check.py','check_source_sha256')):bind(self.hashes,name,self.wraw[key])
        for name in (PREFIX+'five_moment_repair.json',PREFIX+'five_moment_repair_check.json',wn,
            'lei_ren_part1_paper_shared_bump_constants.json',Path(original.__file__).name,
            'lei_ren_part1_paper_compliant_five_moment_repair.py','lei_ren_part1_paper_interval_five_bump_inverse.py',Path(__file__).name):bind(self.hashes,name,sha(name))
        # Hydrate ONLY the stateless beta/partial-weight methods. The old
        # pressure/core/repair constructors and old producers are not run.
        self.bump=SharedFiveMomentRepair.__new__(SharedFiveMomentRepair);self.bump.ctx=c
        self.bump.normalization=restore_value(c,self.wraw['normalization'])
        self.bump.beta_sup=fields.previous.read_interval(c,self.fixed['bump_sup_bound'])
        self.bump.beta_deriv_sup=fields.previous.read_interval(c,self.fixed['bump_first_derivative_sup_bound'])
        self.bump.full_weights={(Fraction(r['center']),Fraction(r['power']),r['multiplicity']):restore_value(c,r['weight_interval']) for r in self.wraw['weight_records'].values()}
        self.bindings=source_bindings()
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked.get('all_passed') or not checked.get(GATE) or checked['source_family']!=self.family:
                raise ValueError('Accepted actual two-frame leading inverse receipt required')
            for path,digest in checked['input_hashes'].items():bind(self.hashes,path,digest)
            bind(self.hashes,RECEIPT,sha(RECEIPT))

    def owner(self,label):
        if label in self.owners:return self.owners[label]
        op=self.upstream.owner(label)
        if not ep(self.admitted_logP)[0]<=ep(op.logP)[0]<=ep(op.logP)[1]<=ep(self.admitted_logP)[1]:
            raise ValueError('Canonical same-source Pstar lies outside original leading admission')
        with mp.workdps(self.c.dps+40):
            self.owners[label]=_ActualRmLeadingPatch(op,self.t,self.W,self.bump)
        return self.owners[label]

    def evaluate(self,label,coordinate):
        with mp.workdps(self.c.dps+40):value=self.owner(label).evaluate(coordinate)
        return dict(mode=self.mode,source_family=self.family,source_frame=label,
            function_evaluation=fields.serialized(value),whole_axis_functions_installed=False,
            no_original_ancestor_producers_or_full_checkers_executed=True)


def run():
    began=time.monotonic();owner=OriginalRmDefectPatchInverse(require_checked=False)
    with mp.workdps(owner.c.dps+40):
        frames={label:dict(coefficient_solution=fields.serialized(pack(owner.owner(label).coefficients())),
            partial_patch=[owner.evaluate(label,x) for x in ((1,1),(5,4),(3,2),(7,4),(71,40),(2,1))]) for label in ('0','.5')}
    report=dict(**{GATE:True},source_family=owner.family,original_source_bindings=owner.bindings,frames=frames,
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual native source frames0,.5: Rm five defects, canonical anisotropic leading inverse Z0..5 and partial bump velocities/pressure/primitives. Whole-axis source providers, full radial mixed4, finite-N Rc/all24/controls/global N and n-recursion stay open.')
    (HERE/NAME).write_text(json.dumps(base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Actual Rm leading defects, normalized axial5 inverse and partial patch generated',flush=True);return report


if __name__=='__main__':run()
