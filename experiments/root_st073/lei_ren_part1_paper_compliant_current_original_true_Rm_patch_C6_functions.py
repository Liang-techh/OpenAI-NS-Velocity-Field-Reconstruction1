"""True prescribed Rm inlet, original leading implicit C6 map and patch.

Exact finite bump integrals and the original unique common-box branch define
the control values. Higher rows solve the same full Jacobian through Z6.
No coefficient enclosure or its midpoint defines a field value. Section11
correction transport remains a separate layer with a nonzero incoming.
"""
import copy
from fractions import Fraction
import gzip
import itertools
import json
from pathlib import Path
import time

import mpmath as mp

import lei_ren_part1_paper_compliant_current_original_modified_Q4_pressure_C3_functions as previous
import lei_ren_part1_paper_compliant_current_original_Rm_defect_patch_inverse as leading
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_centered_inverse as admitted
import lei_ren_part1_paper_compliant_current_original_whole_Z_Rm_patch_mixed4 as old_patch
import lei_ren_part1_paper_compliant_flat_pulse_derivatives as shapes

source,parent,controls=previous.source,previous.parent,previous.controls
normalization=previous.previous.previous.previous
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_true_Rm_patch_C6_functions.json.gz'
RECEIPT=PREFIX+'current_original_true_Rm_patch_C6_functions_check.json'
GATES=('current_original_true_Rm_leading_implicit_control_C6_functions_installed',
       'current_original_true_Rm_patch_background_C6_generic_C5_functions_installed')
OPEN=tuple(dict.fromkeys(previous.OPEN+('current_original_true_Rm_patch_C6_numeric_containment_installed',
    'current_original_true_Rm_patch_finite_N_C5_transport_installed')))
FunctionJet,rows,encode=previous.FunctionJet,previous.rows,previous.encode
Algebra=previous.Algebra
CENTERS=(Fraction(5,4),Fraction(3,2),Fraction(7,4))
RADIUS=Fraction(1,40)


def power(g,x,p):return g.unary('exp',g.mul(g.constant(p),g.unary('log',x)))


def raw_beta(g,x):
    return g.node('compact_raw_beta',argument=x.node,
        exact_definition='exp(-1/(1-x^2)) for |x|<1; exact0 elsewhere',
        original_source=parent.active.current.ast_binding(leading.SharedFiveMomentRepair.beta),
        original_radius_not_absorbed_into_raw_shape=True)


def beta_first_two(g,x,beta):
    c=g.constant;d=g.sub(g.one,g.mul(x,x));d2=g.mul(d,d);d3=g.mul(d2,d);d4=g.mul(d2,d2)
    first=g.mul(beta,g.quotient(g.mul(c(-2),x),d2,'original beta interior 1-x^2>0'))
    second=g.mul(beta,g.add(g.quotient(c(-2),d2,'original beta interior'),
        g.quotient(g.mul(c(-8),x,x),d3,'original beta interior'),
        g.quotient(g.mul(c(4),x,x),d4,'original beta interior')))
    predicate=g.node('logical_or',arguments=[
        g.node('exact_real_comparison',operator='<=',left=x.node,right=c(-1).node).node,
        g.node('exact_real_comparison',operator='>=',left=x.node,right=g.one.node).node])
    lazy=lambda body:g.node('original_lazy_flat_branch',active_body=body.node,Delta=x.node,eta=g.one.node,
        flat_value=g.zero.node,flat_predicate=predicate.node,active_body_not_evaluated_on_flat_branch=True,
        original_beta_exact_flat_derivative=True)
    return lazy(first),lazy(second)


class ExactBumpWeights:
    def __init__(self,g):
        self.g=g;self.free_seeds={};self.cache={};self.partial={}
        s=g.symbol('true_Rm_beta_normalization_s');beta=raw_beta(g,s)
        self.free_seeds[beta.node]={g.nodes[s.node]['name']}
        self.normalization=controls.integral(g,beta,g.nodes[s.node]['name'],g.constant(-1),g.one,
            exact_original_positive_beta_normalization=True,range_or_midpoint_not_defining_integral=True)
        self.free_seeds[self.normalization.node]=set()
    def gamma(self,x,i):
        g=self.g;argument=g.mul(g.constant(40),g.sub(x,g.constant(CENTERS[i])))
        beta=raw_beta(g,argument);d1,d2=beta_first_two(g,argument,beta)
        # Issued coordinates and quadrature variables are Z-independent.
        self.free_seeds[beta.node]={g.nodes[x.node]['name']}
        scale=g.quotient(g.one,self.normalization,'positive exact original beta normalization')
        return (g.mul(g.constant(40),scale,beta),g.mul(g.constant(1600),scale,d1),
                g.mul(g.constant(64000),scale,d2))
    def weight(self,i,p,m=1,upper=None):
        g=self.g;key=(i,Fraction(p),m,upper.node if upper is not None else None)
        cache=self.partial if upper is not None else self.cache
        if key in cache:return cache[key]
        x=g.symbol('true_Rm_weight_%d_%s_%d%s'%(i,str(Fraction(p)).replace('/','_'),m,'_partial' if upper is not None else ''))
        gamma=self.gamma(x,i)[0]
        body=g.mul(power(g,x,p),*([gamma]*m))
        lower=g.one if upper is not None else g.constant(CENTERS[i]-RADIUS)
        right=upper if upper is not None else g.constant(CENTERS[i]+RADIUS)
        result=controls.integral(g,body,g.nodes[x.node]['name'],lower,right,
            exact_original_Rm_bump_weight=True,original_center=str(CENTERS[i]),original_radius=str(RADIUS),
            exact_power=str(Fraction(p)),multiplicity=m,full_original_support=upper is None,
            ordinary_measure='dx',Z_independent_weight=True,range_or_midpoint_not_defining_integral=True)
        if upper is None:self.free_seeds[result.node]=set()
        cache[key]=result;return result
    def table(self,upper=None):
        g=self.g;w=lambda i,p,m=1:self.weight(i,p,m,upper)
        # Integral gamma dx=1 exactly by the normalization substitution.
        mass=[g.one,g.one] if upper is None else [w(0,0),w(2,0)]
        return dict(L=[[mass[0],mass[1],g.zero,g.zero,g.zero],
            [w(0,Fraction(3,5)),w(2,Fraction(3,5)),g.zero,g.zero,g.zero],
            [g.zero,g.zero,*[w(i,Fraction(1,2)) for i in range(3)]],
            [g.zero,g.zero,*[g.neg(w(i,Fraction(1,10))) for i in range(3)]],
            [g.zero,g.zero,*[w(i,Fraction(-9,10)) for i in range(3)]]],
            fg=[w(i,Fraction(1,2),2) for i in (0,2)],gg=[w(i,0,2) for i in (0,2)],
            ff=[w(i,0,2) for i in range(3)],ff_over_x=[w(i,-1,2) for i in range(3)])


def nonlinear(a,W,h,alpha):
    zero=a.fixed(0)
    cross=a.add(a.scale(a.mul(h[0],h[2]),W['fg'][0]),a.scale(a.mul(h[1],h[4]),W['fg'][1]))
    axial=a.add(a.scale(a.mul(h[0],h[0]),W['gg'][0]),a.scale(a.mul(h[1],h[1]),W['gg'][1]))
    swirl=a.add(*(a.scale(a.mul(h[i+2],h[i+2]),W['ff'][i]) for i in range(3)))
    pressure=a.add(*(a.scale(a.mul(h[i+2],h[i+2]),W['ff_over_x'][i]) for i in range(3)))
    return [zero,cross,zero,a.sub(a.mul(alpha,axial),a.scale(swirl,'1/2')),a.scale(pressure,'1/2')]


def jacobian(g,W,h,alpha):
    J=[list(row) for row in W['L']];c=g.constant;add=g.add;mul=g.mul
    J[1][0]=add(J[1][0],mul(W['fg'][0],h[2]));J[1][1]=add(J[1][1],mul(W['fg'][1],h[4]))
    J[1][2]=mul(W['fg'][0],h[0]);J[1][4]=mul(W['fg'][1],h[1])
    for i,j in ((0,0),(1,1)):J[3][i]=mul(c(2),alpha,W['gg'][j],h[i])
    for i in range(3):
        J[3][i+2]=add(J[3][i+2],g.neg(mul(W['ff'][i],h[i+2])))
        J[4][i+2]=add(J[4][i+2],mul(W['ff_over_x'][i],h[i+2]))
    return J


def determinant(g,matrix):
    n=len(matrix)
    return g.add(*(g.mul(g.constant(-1 if sum(p[i]>p[j] for i in range(n) for j in range(i+1,n))%2 else 1),
        *(matrix[i][p[i]] for i in range(n))) for p in itertools.permutations(range(n))))


def inverse_matrix(g,J):
    det=determinant(g,J);den=g.mul(det,det);out=[]
    for i in range(5):
        row=[]
        for j in range(5):
            minor=[[J[r][c] for c in range(5) if c!=i] for r in range(5) if r!=j]
            adj=g.mul(g.constant((-1)**(i+j)),determinant(g,minor))
            row.append(g.quotient(g.mul(adj,det),den,
                'same original common-box Jacobian invertible; determinant squared strictly positive'))
        out.append(row)
    return out,det


def implicit_rows(g,W,defects,alpha,values,order=6):
    """Ordinary coefficient recurrence of the SAME full nonlinear map."""
    J=jacobian(g,W,values,alpha.value);inverse,det=inverse_matrix(g,J)
    out=[[value] for value in values];proof=[]
    for n in range(1,order+1):
        a=Algebra(g,n)
        # Only the unknown order n is temporarily zero in the equation.
        known=[FunctionJet(tuple(q)+(g.zero,)) for q in out]
        Q=nonlinear(a,W,known,a.trim(alpha))
        rhs=[g.add(rows(defects[i])[n],rows(Q[i])[n]) for i in range(5)]
        solved=[g.neg(g.add(*(g.mul(inverse[i][j],rhs[j]) for j in range(5)))) for i in range(5)]
        for i in range(5):out[i].append(solved[i])
        proof.append(dict(order=n,known_rhs=rhs,solved_ordinary_row=solved,
            only_unknown_coefficient_temporarily_zero=True,same_full_Jacobian_used=True))
    return [FunctionJet(tuple(q)) for q in out],dict(Jacobian=J,inverse=inverse,determinant=det,recurrence=proof)


def inlet_defects(g,fields,am,Pstar):
    a=Algebra(g,6);z=FunctionJet((g.symbol('Z'),g.one)+(g.zero,)*5)
    inv=g.quotient(g.one,Pstar,'same positive Pstar');inv2=g.mul(inv,inv)
    mass=a.scale(fields['m'],Pstar);theta=a.div(fields['h'],am,'same positive Rm normalized swirl amplitude')
    mixed=a.scale(a.div(fields['k'],am,'same positive Rm swirl'),Pstar)
    centered=a.add(fields['e'],a.scale(a.mul(z,fields['m']),g.mul(g.constant(-8),inv)),
        a.scale(a.mul(z,z),g.mul(g.constant(16),inv2)))
    energy=a.div(centered,a.mul(am,am),'same positive Rm normalized swirl squared')
    pressure=a.div(fields['p'],a.mul(am,am),'same positive Rm normalized swirl squared')
    return [a.sub(mass,a.scale(z,4)),a.sub(mixed,a.scale(a.mul(z,theta),4)),
        a.sub(theta,a.fixed('5/8')),a.add(energy,a.fixed('5/12')),a.sub(pressure,a.fixed('5/2'))]


def map_conjugacy(g,W,units,alpha,defects):
    """Exact diagonal conjugacy to the admitted original normalized map."""
    a=Algebra(g,0);b=Algebra(g,6);angular=units[2]
    u=[g.symbol('true_Rm_branch_proof_u%d'%i) for i in range(5)]
    h=[g.mul(units[i],u[i]) for i in range(5)]
    Wn={key:value if key in ('L','gg') else [g.mul(q,angular) for q in value] for key,value in W.items()}
    alpha_n=b.scale(alpha,g.quotient(g.mul(units[0],units[0]),angular,'positive original angular normalization unit'))
    dn=[b.scale(q,g.quotient(g.one,units[i],'positive original fixed normalization unit')) for i,q in enumerate(defects)]
    rawQ=nonlinear(a,W,[a.fixed(q) for q in h],a.trim(alpha))
    normQ=nonlinear(a,Wn,[a.fixed(q) for q in u],a.trim(alpha_n))
    rawF=[g.add(defects[i].value,rawQ[i].value,*(g.mul(W['L'][i][j],h[j]) for j in range(5))) for i in range(5)]
    normF=[g.add(dn[i].value,normQ[i].value,*(g.mul(Wn['L'][i][j],u[j]) for j in range(5))) for i in range(5)]
    rawJ=jacobian(g,W,h,alpha.value);normJ=jacobian(g,Wn,u,alpha_n.value)
    return dict(variables=u,raw_controls=h,normalized_weights=Wn,normalized_alpha=alpha_n,
        normalized_defects=dn,raw_map=rawF,normalized_map=normF,
        raw_map_divided_by_units=[g.quotient(q,units[i],'positive original fixed normalization unit') for i,q in enumerate(rawF)],
        conjugated_raw_Jacobian=[[g.quotient(g.mul(rawJ[i][j],units[j]),units[i],
            'positive original fixed normalization unit') for j in range(5)] for i in range(5)],
        normalized_Jacobian=normJ,exact_identity='G(u)=S^-1 F(Su); J_G=S^-1 J_F S',
        original_block_diagonal_L_commutes_with_S=True,
        normalized_alpha_identity='exp(6/5)*(1+Z^2)^2/10^10')


def patch_fields(g,x,am,Pstar,P0,defects,h,alpha,Wpartial,gamma):
    a=Algebra(g,6);z=FunctionJet((g.symbol('Z'),g.one)+(g.zero,)*5);inv=g.quotient(g.one,Pstar,'same positive Pstar')
    changes=nonlinear(a,Wpartial,h,alpha)
    residual=[a.add(defects[i],changes[i],*(a.scale(h[j],Wpartial['L'][i][j]) for j in range(5))) for i in range(5)]
    H=a.add(a.fixed(power(g,x,Fraction(1,10))),*(a.scale(h[i+2],gamma[i][0]) for i in range(3)))
    axial=a.add(*(a.scale(h[i],gamma[j][0]) for i,j in ((0,0),(1,2))))
    Hx=a.add(a.fixed(g.mul(g.constant('1/10'),power(g,x,Fraction(-9,10)))),
        *(a.scale(h[i+2],gamma[i][1]) for i in range(3)))
    axial_x=a.add(*(a.scale(h[i],gamma[j][1]) for i,j in ((0,0),(1,2))))
    mass=a.add(a.scale(z,4),a.scale(residual[0],g.quotient(g.one,x,'positive original x=R/Rm')))
    theta=a.add(residual[2],a.fixed(g.mul(g.constant('5/8'),power(g,x,Fraction(8,5)))))
    mixed=a.add(a.scale(a.mul(z,theta),4),residual[1])
    energy=a.sub(residual[3],a.fixed(g.mul(g.constant('5/12'),power(g,x,Fraction(6,5)))))
    pressure=a.add(residual[4],a.fixed(g.mul(g.constant('5/2'),power(g,x,Fraction(1,5)))))
    am2=a.mul(am,am);angular=power(g,x,Fraction(-3,2));invx=g.quotient(g.one,x,'positive original patch x')
    rawenergy=a.add(a.scale(a.mul(z,mass),g.mul(g.constant(8),inv,inv)),
        a.scale(a.mul(z,z),g.mul(g.constant(-16),inv,inv)),a.scale(a.mul(am2,energy),invx))
    fields=dict(E=a.mul(am,H),V=a.scale(a.add(a.scale(z,4),axial),inv),m=a.scale(mass,inv),
        h=a.scale(a.mul(am,theta),angular),k=a.scale(a.mul(am,mixed),g.mul(angular,inv)),
        e=rawenergy,p=a.mul(am2,pressure),P0=a.trim(P0))
    fields['pressure']=a.add(fields['P0'],fields['p'])
    radial=dict(E_y=a.scale(a.mul(am,Hx),x),V_y=a.scale(axial_x,g.mul(x,inv)))
    return fields,radial,dict(actual_partial_defects=residual,partial_changes=changes,H=H,axial_correction=axial,
        physical_mass=mass,theta=theta,mixed=mixed,centered_energy=energy,pressure_moment=pressure)


class CurrentOriginalTrueRmPatchC6Functions:
    def __init__(self,modified_owner,require_checked=True):
        if type(modified_owner) is not previous.CurrentOriginalModifiedQ4PressureC3Functions or not modified_owner.acceptance_loaded:
            raise TypeError('Checked true prescribed/modified source through Rm required')
        modified_owner.assert_graph();self.modified=modified_owner;self.background=modified_owner.background;self.g=g=modified_owner.g
        self.family=copy.deepcopy(modified_owner.family);self.hashes=dict(modified_owner.hashes)
        for name,digest in parent.linear.read(previous.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for name in (previous.NAME,previous.RECEIPT,admitted.NAME,admitted.RECEIPT,old_patch.RECEIPT,Path(admitted.__file__).name,
            Path(leading.__file__).name,Path(leading.original.__file__).name,Path(shapes.__file__).name):self._pin(name,sha(name))
        admission=parent.linear.read(admitted.RECEIPT);patch=parent.linear.read(old_patch.RECEIPT)
        if (not admission['all_passed'] or not admission[admitted.GATE] or admission['source_family']!=self.family
            or not admission['same_common_unit_box_proves_unique_whole_Z_smooth_leading_family']
            or not admission['actual_leading_controls_and_same_Jacobian_axial5_checked']
            or not patch['all_passed'] or patch['source_family']!=self.family):
            raise ValueError('Same original whole-Z unique leading branch and shape admission required')
        for name,digest in admission['input_hashes'].items():self._pin(name,digest)
        for name,digest in patch['input_hashes'].items():self._pin(name,digest)
        branch_source=parent.linear.read(admitted.NAME)
        self.fixed_N0=g.nodes[self.background.bridge.hb.node]['fixed_leading_input_N0']
        if (branch_source['source_family']!=self.family
            or branch_source['candidate_N_not_used_for_leading_map']!=self.fixed_N0
            or branch_source['exact_Z_domain']!=['-1','1']):
            raise ValueError('Same exact original leading N0 and full-Z source recipes required')
        c=modified_owner.loop.background.bridge.inlet.common.ctx
        for row in admission['contraction_margins'].values():
            for key in ('initial_contraction','higher_Jacobian_contraction'):
                q=c.mpf(mp.make_mpf(tuple(row[key]['upper_exact_mpf_tuple'])))
                if not 0<=q<1:raise ValueError('Original base branch/Jacobian contraction lost')
        self.bindings=dict(original_leading_map=leading.source_bindings(),source_centering=admitted.source_bindings(),
            original_nonlinear=parent.active.current.ast_binding(leading.original.nonlinear_series),
            original_implicit_recurrence=parent.active.current.ast_binding(leading.implicit_axial_jets),
            true_normalized_fields=parent.active.current.ast_binding(normalization.generic_fields),
            true_Rm_source=parent.active.current.ast_binding(type(self.background).__init__),
            original_centered_packet=parent.active.current.ast_binding(leading.previous.ActualReferenceRestoreFunctions.packet),
            original_postrestore=parent.active.current.ast_binding(leading.previous.ActualReferenceRestoreFunctions.postrestore),
            original_leading_branch=parent.active.current.ast_binding(leading._ActualRmLeadingPatch.__init__))
        if branch_source['original_source_bindings']!=self.bindings['source_centering']:
            raise ValueError('Original admitted centered source bindings must equal the actual source recipes')
        self.Pstar=self.background.Pstar;self.P0=self.background.P0;self.inlet=self.background.endpoints['Rm']
        a=Algebra(g,6);z=FunctionJet((g.symbol('Z'),g.one)+(g.zero,)*5)
        self.am=a.div(a.fixed(g.unary('exp',g.constant('-3/5'))),a.add(a.fixed(1),a.mul(z,z)),
            'positive canonical normalized Rm amplitude exp(-3/5)/(1+Z^2)')
        self.alpha=a.scale(a.div(a.fixed(1),a.mul(self.am,self.am),'same positive Rm normalized amplitude squared'),
            g.unary('exp',g.mul(g.constant(-2),self.background.core.parameters['logP'])))
        self.defects=inlet_defects(g,self.inlet['fields'],self.am,self.Pstar)
        self.bump=ExactBumpWeights(g);self.W=self.bump.table()
        weights_file='lei_ren_part1_paper_bump_integral_enclosures_check.json'
        weights_proof=parent.linear.read(weights_file);self._pin(weights_file,sha(weights_file))
        if not weights_proof['directed_integrals_certified']:raise ValueError('Exact original bump integrals lack containment proof')
        original_records={(Fraction(row['center']),Fraction(row['power']),row['multiplicity']):row
            for row in weights_proof['weight_records'].values()}
        weight_bindings=[]
        for (i,p,m,upper),ref in self.bump.cache.items():
            if upper is not None or (CENTERS[i],p,m) not in original_records:
                raise ValueError('Exact new weight must match an original certified integral recipe')
            row=g.nodes[ref.node]
            if (row['original_center']!=str(CENTERS[i]) or row['original_radius']!='1/40'
                or row['exact_power']!=str(p) or row['multiplicity']!=m or not row['full_original_support']):
                raise ValueError('Changed original exact bump integral recipe')
            weight_bindings.append(dict(node=ref.node,center=str(CENTERS[i]),power=str(p),multiplicity=m,
                exact_integral_contained_by_original_directed_record=True,
                original_directed_weight_interval=copy.deepcopy(original_records[(CENTERS[i],p,m)]['weight_interval']),
                containment_proof='Exact integral body/support/normalization equals this pinned original directed integral recipe; its certified directed interval therefore contains the same exact integral.'))
        self.branch_transfer=dict(
            original_whole_Z_unique_branch_receipt=admitted.RECEIPT,
            original_exact_weight_recipe_receipt=weights_file,exact_weight_recipe_bindings=weight_bindings,
            actual_leading_N0_equal_by_exact_integer=True,actual_source_Z_domain=['-1','1'],
            exact_defects_from_original_physical_primitives=True,
            raw_defect_recipes=['M-4Z','Theta_mixed-4Z*Theta','Theta-5/8',
                '(Mztheta-8Z*Mz+16Z^2*Rm)/(Rm*Am^2)+5/12','Mp/Am^2-5/2'],
            original_current_profile_recipes_bound_by_checked_C6_parent=True,
            original_raw_beta_definition_and_finite_support_normalization_unchanged=True,
            original_four_cell_C0_and_Jacobian_contractions_strictly_lt1=True,
            exact_source_correspondence=dict(
                normalized_raw_histories=['M','H/(2*phi)','K/(2*phi)','A','B/phi^2','C/phi^2'],
                actual_Rm_amplitude_identity='sqrt(2*Rm)*F0*phi/Pstar=exp(-3/5)/(1+Z^2)',
                centered_defects=['M-4Z','K/(2*phi)-4Z*H/(2*phi)','H/(2*phi)-5/8',
                    'Am^-2*(A-8Z*M+16Z^2)-(B/phi^2-5/6)/2','(C/phi^2-5)/2'],
                predecessor_theorem='Checked unique prescribed core/bridge/switch/finite-long/reference/restoration source with identical initial data and fixed N0; exact moment ODE normalization and original finite memories identify the two representations.',
                acceptance_requires_actual_node_normalization_and_exact_weight_recipe_checks=True,
                source_bindings=self.bindings),
            exact_Jacobian_invertibility_proof='Same exact L,Q,defects,amplitude,source N0 and normalized common box: original directed bound ||I-RJ||<1 applies to exact integral weights and exact defining source functions, hence the Neumann inverse exists and det(J)^2>0.',
            sixth_function_derivative_is_implicit_theorem_extension_not_new_numeric_norm=True)
        repair=parent.linear.read(PREFIX+'five_moment_repair.json');t=repair['axial_box_scale']
        if t['lower_exact_mpf_tuple']!=t['upper_exact_mpf_tuple']:raise ValueError('Exact fixed branch box scale required')
        sign,man,exponent,_=t['lower_exact_mpf_tuple'];fixed=Fraction((-1 if sign else 1)*man)*Fraction(2)**exponent
        self.t=g.constant(fixed);angular=g.mul(g.constant(10**10),self.t,self.t,
            g.unary('exp',g.mul(g.constant(-2),self.background.core.parameters['logP'])))
        self.units=[self.t,self.t,angular,angular,angular]
        self.conjugacy=map_conjugacy(g,self.W,self.units,self.alpha,self.defects)
        self.branch_transfer['raw_to_original_normalized_conjugacy']=encode(self.conjugacy)
        self.system=g.node('original_Rm_unique_leading_implicit_system',source_family=self.family,
            defining_map='L*h+Q(h,invAm2)+true_Rm_defects=0',exact_weight_nodes=encode(self.W),
            exact_defect_value_nodes=[q.value.node for q in self.defects],invAm2=self.alpha.value.node,
            exact_Z_independent_row_and_control_units=[q.node for q in self.units],
            original_common_box_branch_receipt=admitted.RECEIPT,branch_receipt_sha256=sha(admitted.RECEIPT),
            original_whole_Z_base_unique_branch_and_Jacobian_admitted=True,
            original_leading_fixed_point_not_Rc_finite_N_control=True,
            fixed_leading_input_N0=self.fixed_N0,exact_source_and_weight_branch_transfer=self.branch_transfer,
            source_projection_binding=parent.active.current.ast_binding(admitted.WholeZRmCenteredInverse.owner),
            source_bindings=self.bindings,numerical_enclosures_not_values=True)
        values=[g.node('original_Rm_unique_leading_implicit_control',system=self.system.node,component=i,
            Z_variable=g.symbol('Z').node,source_family=self.family,
            exact_value_of_same_unique_original_common_box_solution=True) for i in range(5)]
        self.control_functions,self.implicit=implicit_rows(g,self.W,self.defects,self.alpha,values)
        b=Algebra(g,0);base=[b.fixed(v) for v in values];Q=nonlinear(b,self.W,base,b.trim(self.alpha))
        self.base_residual=[g.add(self.defects[i].value,Q[i].value,
            *(g.mul(self.W['L'][i][j],values[j]) for j in range(5))) for i in range(5)]
        self.base_identity=dict(system=self.system,residual=self.base_residual,
            defining_exact_unique_solution_identity_not_numeric_zero_containment=True,
            same_actual_L_Q_defect_and_amplitude_source=True,source_transfer=self.branch_transfer)
        self.x=x=g.symbol('true_Rm_patch_x');self.gamma=[self.bump.gamma(x,i) for i in range(3)]
        self.partial_weights=self.bump.table(x)
        fields,radial,primitives=patch_fields(g,x,self.am,self.Pstar,self.P0,self.defects,
            self.control_functions,self.alpha,self.partial_weights,self.gamma)
        R=g.mul(self.inlet['radius'],x);jac=g.quotient(g.one,x,'positive original actual patch x')
        generic=previous.generic_functions(g,fields,self.background.delta,self.Pstar,R,order=5)
        self.chart=dict(phase=x,domain=[g.one,g.unary('exp',g.one)],radius=R,physical_measure=jac,
            fields=fields,radial=radial,generic=generic,primitives=primitives)
        self.endpoints={name:self.at(value) for name,value in (('Rm',g.one),('Rh',g.unary('exp',g.one)))}
        self.free_source_seeds=dict(self.bump.free_seeds)
        for ref in values:self.free_source_seeds[ref.node]={'Z'}
        self.free_source_seeds[self.system.node]={'Z'}
        self.definition=dict(true_prescribed_Rm_C6_input_not_old_axial5_enclosure_adapter=True,
            exact_original_beta_weights_normalization_centers_radius_and_support=True,
            original_unique_common_box_branch_same_full_Jacobian_through_Z6=True,
            all_six_implicit_derivatives_include_invAm2_Z_products=True,
            original_five_centered_defects_from_true_m_h_k_e_p_without_P0=True,
            leading_control_not_finite_N_correction_and_no_history_reset=True,
            actual_E_V_m_h_k_e_p_P0_pressure_and_Ey_Vy_C6=True,full_signed_generic_C5=True,
            original_x1_inlet_and_Rh_exp1_exit_and_log_Jacobian_1_over_x=True,
            original_branch_zero_residual_is_function_identity_not_zero_containment=True,
            control_C6_numeric_enclosure_and_mixed_norm_admission_not_certified=True,
            Section11_nonzero_incoming_correction_transport_not_installed=True,
            six_outer_charts_full_preRp_n1_heat_stress_recursion_and_NS_open=True,
            shared_graph_unreferenced_suffix_allowed_after_exact_pinned_prefix=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.modified,self.background,self.g,self.bump);self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None);checker=pins.pop(Path(__file__).stem+'_check.py',None)
            if (not receipt.get('all_passed') or not all(receipt.get(k) for k in GATES) or any(receipt.get(k) for k in OPEN)
                or not receipt.get('same_exact_original_unique_branch_theorem_bound')
                or receipt['source_family']!=self.family or pins!=self.hashes or raw['input_hashes']!=self.hashes
                or archive!=sha(NAME) or checker!=sha(Path(__file__).stem+'_check.py')
                or raw['true_Rm_patch_C6_bindings']!=self.snapshot or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked true Rm leading C6 control/patch source required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed true Rm source '+name)
            self.acceptance_loaded=True
    def at(self,value):
        a=Algebra(self.g,6);b=Algebra(self.g,5);chart=self.chart
        return dict(fields={k:a.at(jet,self.x,value) for k,jet in chart['fields'].items()},
            radial={k:a.at(jet,self.x,value) for k,jet in chart['radial'].items()},
            generic={k:b.at(jet,self.x,value) for k,jet in chart['generic'].items()},
            radius=self.g.node('function_substitution',expression=chart['radius'].node,
                variable=self.x.node,value=value.node,Z_independent_substitution=True))
    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed true Rm input '+name)
        self.hashes[name]=digest
    def _snapshot(self):
        return encode(dict(inlet=self.inlet,am=self.am,alpha=self.alpha,defects=self.defects,W=self.W,
            partial_weights=self.partial_weights,units=self.units,system=self.system,controls=self.control_functions,
            implicit=self.implicit,conjugacy=self.conjugacy,base_identity=self.base_identity,branch_transfer=self.branch_transfer,chart=self.chart,endpoints=self.endpoints,
            exact_beta_normalization=self.bump.normalization,
            free_source_seeds=[dict(node=i,variables=sorted(names)) for i,names in self.free_source_seeds.items()],
            definition=self.definition,bindings=self.bindings))
    def assert_graph(self):
        self.modified.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot
            or self.hashes!=self._hashes or self._owners!=(self.modified,self.background,self.g,self.bump)
            or self.inlet is not self.background.endpoints['Rm']):
            raise ValueError('True Rm inlet, exact implicit branch or patch functions changed')


def run(modified_owner):
    began=time.monotonic();owner=CurrentOriginalTrueRmPatchC6Functions(modified_owner,False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,true_Rm_patch_C6_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_RM_LEADING_CONTROL_AND_PATCH_C6_DEFINED',len(owner.g.nodes),flush=True);return owner,raw


if __name__=='__main__':raise SystemExit('run(modified_owner) reuses the admitted true Rm endpoint')
