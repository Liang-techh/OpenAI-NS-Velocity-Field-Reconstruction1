"""Source-owned signed Rm defects and unknown signed finite-N perturbations.

Center the raw moments analytically before interval arithmetic. In particular,
never subtract order-one rounded histories to recover an exponentially small
defect. This exports the inlet; it does not certify the finite-N active patch.
"""
import hashlib
import json
import math

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_reference_numeric import OriginalReferenceNumeric


def centered_frame(ctx, J, Z, j, epsilon, invP, am, six, base, angular,
                   increments=None, inverse_exp=None):
    """Taylor-Z2 centered identities, including the full signed perturbation.

    angular contains H/Phi-5/4, B/Phi^2-5/6 and C/Phi^2-5 separately.
    inverse_exp encloses exp(-a)-1 and exp(-2a)-1, where a=A/N.
    """
    z=J.variable(ctx,ctx.mpf(Z),2);zero=J.constant(ctx,0,2)
    inc={n:zero for n in ('m','h','k','e','p')} if increments is None else increments
    e1,e2=(zero,zero) if inverse_exp is None else inverse_exp
    mass=base['gamma1']*j+six['Mpsi']*epsilon
    mixed=(base['Kalpha_over_Phi']*j+six['Kpsi_over_Phi']*epsilon)/2
    square=base['gamma2']*(j*j)+base['Dpsi']*(2*j*epsilon)+six['Apsi']*(epsilon*epsilon)
    energy=square*(invP*invP)/(am*am)
    theta=angular['H']/2+ctx.mpf(5)/8
    dm=inc['m']/invP
    di=(inc['k']/invP-z*inc['h']*4)/am
    dh=inc['h']/am
    de=(inc['e']-z*inc['m']*(8*invP))/(am*am)
    dp=inc['p']/(am*am)
    # Expand every factor 1+ei. Forming it first would round an unknown
    # exponentially small perturbation into a much larger order-one ulp.
    return [mass+dm,
        mixed+mixed*e1+di+di*e1,
        angular['H']/2+theta*e1+dh+dh*e1,
        energy-angular['B']/2+energy*e2-angular['B']*e2/2-e2*(ctx.mpf(5)/12)+de+de*e2,
        angular['C']/2+(angular['C']+5)*e2/2+dp+dp*e2]


class OriginalRmDefectNumeric:
    """Arbitrary interior-Z Rm inlet in original source and Pstar units."""
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        self.reference=OriginalReferenceNumeric(precision=precision,
            amplitude_precision=amplitude_precision,directory=directory)
        f=self.reference;self.ctx,self.iv,self.J=f.ctx,f.iv,f.J
        self.family=f.family;self.hashes=dict(f.hashes);self.code_hashes=dict(f.code_hashes)
        path=f.leading.core.directory/(PREFIX+'compliant_current_original_reference_velocity.json')
        old=json.loads(path.read_text(encoding='utf-8'))
        if not old['checks']['passed'] or old['source_family']!=self.family:
            raise ValueError('Accepted actual reference/restoration source required')
        for name,digest in old['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed accepted reference implementation: '+name)
            self.code_hashes[name]=digest
        for name,digest in old['input_hashes'].items():
            if self.hashes.get(name)!=digest:
                if hashlib.sha256((f.leading.core.directory/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed accepted reference input: '+name)
                self.hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        for stem in ('compliant_current_original_prescribed_bridge_C6_functions',
                     'compliant_current_original_true_Rm_patch_C6_functions',
                     'compliant_current_generic_shear_moment_recovery'):
            path=f.leading.core.directory/(PREFIX+stem+'.py')
            digest=hashlib.sha256(path.read_bytes()).hexdigest()
            if self.hashes.get(path.name,digest)!=digest:raise ValueError('Changed Rm frame source')
            self.hashes[path.name]=digest
        self.logP,self.logPp=f.logP,f.logPp
        # This is an error scale, never an input field or its selected value.
        self.tail_log=-64*self.logP;self.tail=self.iv.exp(self.tail_log)
        self._cache={};self._datum={}
        name='src/openai_ns_reconstruction/original_Rm_defect_numeric.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _small(self,log,proof,label):
        if log is None:return self.iv.mpf(0)
        if _bounds(self.ctx,log)[1]>_bounds(self.ctx,self.tail_log)[0]:
            raise ArithmeticError('Original source does not admit the Rm error scale: '+label)
        proof.append(dict(quantity=label,source_log_upper=self.reference.owner._source_summary(log),
            bound_log=self.reference.owner._source_summary(self.tail_log),
            bound_is_not_signed_field_value=True))
        return self.tail

    def _angular(self,Z,proof):
        f=self.reference;c,iv,J=self.ctx,self.iv,self.J;inp=f.inputs(Z);outp={};outt={}
        # Normalized angular forcing is constant throughout these three
        # charts. Propagate the centered incoming memory across L+1+1.
        distance=f.length+2
        for short,name,target,rate in (('H','H_over_Phi','1.25','1.6'),
                ('B','B_over_Phi2',iv.mpf(5)/6,'1.2'),('C','C_over_Phi2',5,'.2')):
            initial=inp['six'][name][1]-J.constant(iv,iv.mpf(target),2)
            errors=[]
            for k in range(3):
                mag=f.leading._abs(initial[k])
                log=None if _bounds(c,mag)[1]==0 else iv.ln(mag)-iv.mpf(rate)*distance
                errors.append(f.leading._symmetric(self._small(log,proof,'centered_'+short+'_Z'+str(k))))
            outp[short]=J.constant(c,0,2);outt[short]=J(iv,errors)
        return outp,outt

    def _datum_pair(self,Z):
        if Z not in self._datum:
            core=self.reference.leading.core;c,iv,J=self.ctx,self.iv,self.J;out=[]
            for directed,ctx in ((False,c),(True,iv)):
                fixed,_,_,_=core._seed(Z,4,1,directed=directed)
                # Core ODE pressure = epsilon*physical P0. Its fixed rows
                # have unit exp(2logP-logLambda); generic datum is P0/Pstar².
                lp=core.lp_iv if directed else core.logP
                ll=core.ll_iv if directed else core.logLambda
                factor=ctx.exp(ctx.mpf(ll)-2*ctx.mpf(lp))
                out.append(J(ctx,[ctx.mpf(v)*factor for v in fixed['P0_Z_taylor'][:6]]))
            self._datum[Z]=out
        return self._datum[Z]

    def evaluate(self,Z):
        if Z in self._cache:return self._cache[Z]
        f=self.reference;c,iv,J=self.ctx,self.iv,self.J
        # This call validates Z and retains the entire original correction prefix.
        inlet=f.postrestore(Z,'1');proof=[]
        six,base,_=f._post_six(Z,c.mpf(1),iv.mpf(1),proof)
        ap,at=self._angular(Z,proof)
        inc={}
        for name,rows in inlet['original_five_correction_history_bounds'].items():
            errors=[f.leading._symmetric(self._small(f.owner._source_log(row),proof,
                'raw_'+name+'_Z'+str(k))/math.factorial(k)) for k,row in enumerate(rows[:3])]
            inc[name]=J(iv,errors)
        # Recover original A/N logs, not the serialized1e-1000 numeric cap.
        logs={p['quantity']:p['original_log_upper'] for p in inlet['small_error_proofs']
              if p.get('quantity','').startswith('reference_A_over_N_Z')}
        ab=[]
        for k in range(3):
            row=logs.get('reference_A_over_N_Z'+str(k))
            ab.append(iv.mpf(0) if row is None else self._small(f.owner._source_log(row),proof,
                'local_A_over_N_Z'+str(k))/math.factorial(k))
        if _bounds(c,ab[0])[1]>=c.mpf('.25'):raise ArithmeticError('Small exponential perturbation required')
        def exponential_increment(multiplier):
            a0,a1,a2=[x*multiplier for x in ab]
            # |exp(+-a0)|<2; ordinary derivatives converted to Taylor once.
            return J(iv,[f.leading._symmetric(2*a0),f.leading._symmetric(2*a1),
                         f.leading._symmetric(2*a2+a1*a1)])
        e1,e2=exponential_increment(1),exponential_increment(2)
        pp,tt={},{}
        for directed,ctx in ((False,c),(True,iv)):
            z=J.variable(ctx,ctx.mpf(Z),2)
            am=(1+z*z).reciprocal()*ctx.exp(ctx.mpf('-.6'))
            j=ctx.mpf(f.owner.amplitude.j_iv if directed else f.owner.amplitude.j)
            eps=ctx.mpf(f.leading.core.epsilon_iv if directed else f.leading.core.epsilon)
            invP=ctx.exp(-ctx.mpf(self.logP if directed else self.logPp))
            ss={n:pair[directed] for n,pair in six.items()}
            bb={n:(pair[directed] if n in ('Kalpha_over_Phi','Dpsi') else J.constant(ctx,pair[directed],2))
                for n,pair in base.items()}
            angular=at if directed else ap
            original=centered_frame(ctx,J,Z,j,eps,invP,am,ss,bb,angular)
            canonical=centered_frame(ctx,J,Z,j,eps,invP,am,ss,bb,angular,
                increments=inc if directed else None)
            modified=centered_frame(ctx,J,Z,j,eps,invP,am,ss,bb,angular,
                increments=inc if directed else None,inverse_exp=(e1,e2) if directed else None)
            # Raw generic S=Pstar moment units, with correlations preserved
            # separately in the signed centered frame above.
            zero=J.constant(ctx,0,2);dd=inc if directed else dict.fromkeys(inc,zero)
            mass_change=bb['gamma1']*j+ss['Mpsi']*eps
            M=z*4+mass_change
            centered=bb['gamma2']*(j*j)+bb['Dpsi']*(2*j*eps)+ss['Apsi']*(eps*eps)
            Kover2=(z*4)*(angular['H']/2+ctx.mpf(5)/8)+(bb['Kalpha_over_Phi']*j+ss['Kpsi_over_Phi']*eps)/2
            raw=dict(m=M*invP+dd['m'],h=am*(angular['H']/2+ctx.mpf(5)/8)+dd['h'],
                k=am*Kover2*invP+dd['k'],
                e=((z*z)*16+z*mass_change*8+centered)*(invP*invP)
                    -(am*am)*(angular['B']+ctx.mpf(5)/6)/2+dd['e'],
                p=(am*am)*(angular['C']+5)/2+dd['p'])
            data=dict(original=original,canonical=canonical,modified=modified,raw=raw,am=am,
                invAm2=(am*am).reciprocal()*(invP*invP))
            (tt if directed else pp).update(data)
        p0p,p0t=self._datum_pair(Z)
        packet=dict(Z=Z,chart='unpatched_Rm',source_family=self.family,
            selected_integer=inlet['selected_integer'],
            coordinate=inlet['coordinate'],
            signed_unmodified_defects=[f._rows(p,t) for p,t in zip(pp['original'],tt['original'])],
            actual_finite_N_defect_enclosures=[f._rows(p,t) for p,t in zip(pp['modified'],tt['modified'])],
            actual_finite_N_canonical_frame_defect_enclosures=[f._rows(p,t) for p,t in zip(pp['canonical'],tt['canonical'])],
            generic_own_moment_enclosures={n:f._rows(pp['raw'][n],tt['raw'][n]) for n in inc},
            normalized_Rm_amplitude=f._rows(pp['am'],tt['am']),
            inverse_physical_amplitude_squared=f._rows(pp['invAm2'],tt['invAm2']),
            actual_inverse_physical_amplitude_squared=f._rows(pp['invAm2'],tt['invAm2']+tt['invAm2']*e2),
            amplitude_inverse_exponential_increments=dict(exp_minus_a_minus1=f._rows(J.constant(c,0,2),e1),
                exp_minus_2a_minus1=f._rows(J.constant(c,0,2),e2)),
            finite_N_signed_unknown_history_errors={n:f._rows(J.constant(c,0,2),v) for n,v in inc.items()},
            independent_P0_over_Pstar_squared=f._rows(p0p,p0t),
            own_pressure_over_Pstar_squared=f._rows(p0p.truncate(2)+pp['raw']['p'],p0t.truncate(2)+tt['raw']['p']),
            centered_angular_memories={n:f._rows(ap[n],at[n]) for n in ap},
            original_five_correction_history_bounds=inlet['original_five_correction_history_bounds'],
            positive_order_n1_own_history=inlet['positive_order_n1_own_history'],
            source_contract=dict(generic_unit='S=Pstar',history_units=dict(m='Mz/(R*Pstar)',
                h='Mtheta/(sqrt(2)*R^1.5*Pstar)',k='Mtheta_z/(sqrt(2)*R^1.5*Pstar²)',
                e='Mztheta/(R*Pstar²)',p='Mp/Pstar²'),
                am='exp(-3/5)/(1+Z²); unmodified normalized amplitude',
                defects=['M-4Z','K/(2Phi)-4Z*H/(2Phi)','H/(2Phi)-5/8',
                    '(A-8ZM+16Z²)/(Pstar²*am²)-(B/Phi²-5/6)/2','(C/Phi²-5)/2'],
                modified_normalization='am_N=am*exp(A/N); inverse exponential perturbations retained',
                analytic_centering_before_materialization=True,
                datum='core ODE fixed P0 rows times exp(logLambda-2logP) = physical P0/Pstar²',
                P0_added_once_and_excluded_from_defects=True,
                finite_N_correction_centers_not_resolved=True,unknown_signed_errors_not_field_values=True,
                original_N_and_all_incoming_histories_preserved=True),
            small_error_proofs=f.previous._serializable_proof(proof),
            limits=dict(signed_defect_Z_order=2,P0_Z_order=5,raw_generic_histories_Z_order=2,
                finite_N_active_patch_uniform_inverse=False,patch_controls_or_fields=False,
                Z3_or_higher_defects=False,global_N_admission=False,global_field=False,
                genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._cache[Z]=packet
        return packet
