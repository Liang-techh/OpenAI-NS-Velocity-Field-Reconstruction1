"""Original n1 flat velocity cutoff with own mean/Q and pressure histories.

Coordinates are q=s/(original sc/4) on the selected inner chart. The raw
solution is used only inside its admitted domain. Beyond velocity support,
zero cut velocities retain their own primitive and known leading pressure.
"""
import hashlib
import importlib
import json

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _box, _bounds
from .original_first_collar_axial2 import OriginalFirstCollarAxial2
from .original_first_collar_transport import WidthChange


def stable_sigma_point(c,t):
    """Numerical center only; directed sigma_jets supplies the true value."""
    if t<=0:return c.mpf(0)
    if t>=1:return c.mpf(1)
    if t<c.mpf('.001'):return c.mpf(0)
    if t>c.mpf('.999'):return c.mpf(1)
    w=1/t**2-1/(1-t)**2
    if w>=0:
        e=c.exp(-w);return e/(1+e)
    e=c.exp(w);return 1/(1+e)


class OriginalN1CutoffNumeric(OriginalFirstCollarAxial2):
    def __init__(self,*,precision=85,directory=None,control_cells=128):
        super().__init__(precision=precision,directory=directory)
        if type(control_cells) is not int or not 32<=control_cells<=256:
            raise ValueError('Original cutoff control quadrature cells32..256 required')
        self.control_cells=control_cells
        receipts=(('compliant_current_original_n1_lower_collar_functions',
                   'current_original_n1_Rin_to_rminus_true_lower_collar_C1_functions_installed'),
                  ('compliant_current_original_n1_inner_cutoff_functions',
                   'current_original_n1_true_inner_cutoff_and_five_endpoint_functions_installed'),
                  ('compliant_current_original_n1_cutoff_numeric_bridge',
                   'current_original_true_cutoff_n1_own_fields_and_five_inner_moment_enclosures_installed'),
                  ('compliant_current_original_n1_common_inner_selection',
                   'current_original_n1_common_inner_geometry_below_known_modification_starts_selected'))
        for stem,gate in receipts:
            path=self.leading.core.directory/(PREFIX+stem+'_check.json')
            receipt=json.loads(path.read_text(encoding='utf-8'))
            if not receipt['all_passed'] or not receipt[gate] or receipt['source_family']!=self.family:
                raise ValueError('Accepted same-family original cutoff/selection receipt required')
            for name,digest in receipt['input_hashes'].items():
                if hashlib.sha256((self.leading.core.directory/name).read_bytes()).hexdigest()!=digest:
                    raise ValueError('Changed actual cutoff dependency: '+name)
                if self.hashes.get(name,digest)!=digest:raise ValueError('Conflicting cutoff source')
                self.hashes[name]=digest
            self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        path=self.leading.core.directory/(PREFIX+'compliant_current_inner_exit_strict_collar.json')
        record=json.loads(path.read_text(encoding='utf-8'))
        if not record['exact_current_exit_source_attachment']['passed']:
            raise ValueError('Original exact sc attachment required')
        bits=record['explicit_current_inner_exit_strict_collar']['selected_first_phase_endpoint']
        if bits['lower_exact_mpf_tuple']!=bits['upper_exact_mpf_tuple']:
            raise ValueError('Original selected sc must be an exact singleton')
        self.sc_tuple=bits['lower_exact_mpf_tuple']
        t=list(self.sc_tuple);t[2]-=2;self.phase_tuple=t
        self.phase_exact=mp.make_mpf(tuple(t))
        self.phase_iv=self.iv.mpf([self.phase_exact,self.phase_exact])
        self.phase_upper=_bounds(self.ctx,self.phase_iv)[1]
        if not 0<self.phase_exact<self.ctx.mpf('.0625'):
            raise ValueError('Original sc/4 must stay inside admitted raw inner domain')
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        relative='src/openai_ns_reconstruction/original_n1_cutoff_numeric.py'
        self.code_hashes[relative]=hashlib.sha256((ROOT/relative).read_bytes()).hexdigest()
        self._raw_cache={};self._mass_cache={}
        self.control_masses=importlib.import_module(PREFIX+'compliant_actual_bridge_integrals').pulse_control_masses

    def _jump(self,v,s):
        o,c,iv,J=self.leading,self.ctx,self.iv,self.Jet
        point,truth={},{}
        for directed,label in ((False,'point'),(True,'truth')):
            ctx=iv if directed else c;z=J.variable(ctx,ctx.mpf(v['Z']),2)
            de=o.core.delta_iv if directed else o.core.delta;P2=iv.exp(2*o.core.lp_iv) if directed else c.exp(2*o.core.logP)
            L=1-z*z*de;radial=v['inlet'][label]['radial']
            ph=radial['Phi'][1].truncate(2);ps=radial['Psi'][1].truncate(2)/P2
            rows=dict(p=z*ps*(-8)/L,fx=z*z*ph*32/(L*L),ux=z*z*ps*64/(L*L))
            (truth if directed else point).update(rows)
        # Finite-radius/ratio variations of the complete singular measure.
        C=lambda b,err=0:WidthChange(o,b,err)
        z=C(J.variable(iv,iv.mpf(v['Z']),2));de=v['de'];L=1-z*z*de
        rho=C(J.constant(iv,4,2),4*iv.mpf(s)*iv.exp(o.hcap*iv.mpf(s)))
        x=C(J.constant(iv,2,2),iv.mpf(s)*iv.exp(o.hcap*iv.mpf(s)/2))
        one=C(J.constant(iv,1,2),iv.mpf(1))
        a=v['inlet']['truth'];phi=C(a['values']['Phi'].truncate(2),v['delta']['Phi'])
        comp=C(a['radial']['Phi'][0].truncate(2),self._core_comparison_changes('Phi',0,s,2))
        ratio=phi/comp
        ph=C(a['radial']['Phi'][1].truncate(2),self._core_comparison_changes('Phi',1,s,2))
        ps=C(a['radial']['Psi'][1].truncate(2)/v['P2'],self._core_comparison_changes('Psi',1,s,2)/v['P2'])
        varying=dict(p=one*rho*z*ratio*ps*(-2)/L,
                     fx=one*rho*x*z*z*ratio*ph*4/(L*L),
                     ux=one*rho*x*z*z*ratio*ps*8/(L*L))
        sigma=stable_sigma_point(c,c.mpf(s))
        with mp.workdps(280):sigbox=self.sigma_jets(iv,iv.mpf(s))[0]
        return point,truth,{name:q.error for name,q in varying.items()},sigma,sigbox

    def _core_jets(self,Z):
        if Z not in self._core:self._core[Z]=self.n1.evaluate(Z,'4',degree=20,kernel_degree=28)
        packet=self._core[Z];point,truth={},{}
        for name,rows in packet['normalized_original_n1'].items():
            point[name]=self.Jet(self.ctx,[self.ctx.mpf(q['approximation']) for q in rows])
            truth[name]=self.Jet(self.iv,[self.iv.mpf([mp.make_mpf(tuple(t)) for t in q['original_source_enclosure_exact_mpf_tuples']]) for q in rows])
        return point,truth

    def _raw_uniform(self,Z):
        if Z in self._raw_cache:return self._raw_cache[Z]
        o,iv=self.leading,self.iv;s=self.phase_upper
        v=self._inputs(Z,s,order=5);v['Z']=Z
        forcing=self._forcing(v);linear=self._regular_operator_bound(v)
        _,jt,variation,_,_=self._jump(v,s)
        with mp.workdps(280):sig=self.sigma_jets(iv,iv.mpf([0,s]))[0]
        cp,ct=self._core_jets(Z);zero=self.Jet.constant(iv,0,2)
        scales=dict(f=iv.mpf(1),u=iv.mpf(1),k=iv.mpf(1),p=iv.mpf(1),fx=v['e']**2,ux=v['e']*v['P2'])
        errors={};truth={}
        for name in self.fields:
            source=forcing[name].rows.get(0,zero) if name in forcing else zero
            regular=o._width_error(iv.mpf(s)*iv.mpf('1.5')*(o._norm(linear[name])+o._norm(source)))
            singular=(o._norm(jt[name])*o._abs(sig)+o._width_error(o._abs(sig)*variation[name]))*scales[name] if name in jt else iv.mpf(0)
            errors[name]=regular+singular
            truth[name]=ct[name]+o._norm_error_jet(errors[name],2)
        self._raw_cache[Z]=(v,cp,ct,truth,errors)
        return self._raw_cache[Z]

    def _cutoff(self,q):
        c,iv=self.ctx,self.iv
        if _bounds(c,iv.mpf(q))[1]<=c.mpf('.5'):return c.mpf(1),iv.mpf(1),iv.mpf(0)
        if _bounds(c,iv.mpf(q))[0]>=c.mpf('.75'):return c.mpf(0),iv.mpf(0),iv.mpf(0)
        arg=4*iv.mpf(q)-2
        with mp.workdps(280):sig=self.sigma_jets(iv,arg)
        return 1-stable_sigma_point(c,4*c.mpf(q)-2),1-sig[0],-4*sig[1]

    def _mass(self,q):
        c,iv=self.ctx,self.iv;value=c.mpf(q)
        if _bounds(c,iv.mpf(q))[1]<=c.mpf('.5'):return value,iv.mpf(q)
        if _bounds(c,iv.mpf(q))[0]>=c.mpf('.75'):return c.mpf('0.625'),iv.mpf('0.625')
        if q not in self._mass_cache:
            arg=4*value-2
            point=c.mpf('.5')+c.quad(lambda t:1-stable_sigma_point(c,t),[0,arg/2,arg])/4
            with mp.workdps(280):box,_=self.control_masses(iv,4*iv.mpf(q)-2,cells=self.control_cells)
            self._mass_cache[q]=point,iv.mpf('.5')+box/4
        return self._mass_cache[q]

    def evaluate(self,Z,q):
        """Actual cut n1 on q in[0,1], with phase=q*original(sc/4)."""
        return self._evaluate(Z,q,None)

    def lower_collar(self,Z,fraction):
        """Actual Rin-to-rminus own history, fraction=s/(original sc/4) in[1,2]."""
        return self._evaluate(Z,None,fraction)

    def _evaluate(self,Z,q,s):
        c,iv,o=self.ctx,self.iv,self.leading
        coordinate=q if q is not None else s
        if (type(Z) is not str or type(coordinate) is not str or not c.isfinite(c.mpf(Z))
                or not c.isfinite(c.mpf(coordinate)) or not -1<c.mpf(Z)<1
                or not (1 if q is None else 0)<=c.mpf(coordinate)<=(2 if q is None else 1)):
            raise ValueError('Interior decimal Z and inner q[0,1] or lower fraction[1,2] required')
        exterior=q is None
        lo,hi=_bounds(c,iv.mpf(coordinate))
        if lo<(1 if exterior else 0) or hi>(2 if exterior else 1):
            raise ValueError('Directed selected fraction lies outside admitted chart')
        raw,cp,ct,rt,errors=self._raw_uniform(Z)
        current_phase=self.phase_iv*iv.mpf(coordinate)
        bound_phase=_bounds(c,current_phase)[1]
        if exterior:chi,chii,chider=c.mpf(0),iv.mpf(0),iv.mpf(0);qmass='1'
        else:chi,chii,chider=self._cutoff(q);qmass=q
        mass,massi=self._mass(qmass)
        v=self._inputs(Z,bound_phase,order=5);v['Z']=Z
        rho=v['rho'].truncate(2);x=v['x'].truncate(2)
        # Same actual cutoff history. Its primitive is exactly constant once q>=3/4.
        rho_change=o._width_error(4*self.phase_iv*iv.exp(o.hcap*self.phase_iv))
        integrand_error=o._norm(raw['rho'].truncate(2))*errors['u']+rho_change*o._norm(ct['u'])
        hp=cp['u']*(4*mass)
        ht=ct['u']*(4*massi)+o._norm_error_jet(o._abs(massi)*integrand_error,2)
        history_error=o._width_error(self.phase_iv*o._norm(ht))
        primitive_p=cp['mean_u']*4
        primitive_t=ct['mean_u']*4+o._norm_error_jet(history_error,2)
        meanp=primitive_p/4;meant=primitive_t/rho
        fp,ft=cp['f']*chi,rt['f']*chii
        up,ut=cp['u']*chi,rt['u']*chii
        kp,kt=meanp-up,meant-ut
        def recover(u,m,ctx,de):
            z=self.Jet.variable(ctx,ctx.mpf(Z),2)
            return ((2*z*u-z*(1+de)*m-(1-z*z)*o._derivative(m))/(1-z*z*de)).truncate(1)
        qp,qt=recover(up,meanp,c,o.core.delta),recover(ut,meant,iv,o.core.delta_iv)
        # Direct pressure history: original leading forcing once, plus cut F coupling.
        forcing=self._forcing(v);jp,jt,variation,sig,sigbox=self._jump(v,bound_phase)
        known_regular=o._width_error(current_phase*iv.mpf('1.5')*o._norm(forcing['p'].rows[0]))
        known_change=o._width_error(o._abs(sigbox)*variation['p'])
        coupling=(raw['S']*raw['Phi'][0].rows[0]*rt['f']*(2/(raw['e']**3*raw['P2']))).truncate(2)
        coupling_norm=o._abs(massi)*o._norm(raw['rho'].truncate(2)*coupling)
        coupling_error=o._width_error(self.phase_iv*coupling_norm)
        pp=cp['p']+jp['p']*sig
        pt=ct['p']+jt['p']*sigbox+o._norm_error_jet(known_regular+known_change+coupling_error,2)
        out={name:o._jet_summary(p,t) for name,p,t in
             (('f',fp,ft),('u',up,ut),('mean_u',meanp,meant),('k',kp,kt),('q',qp,qt),('p',pp,pt))}
        # First radial derivative sectors: hb^-1*phase_in^-1 is never materialized.
        def radial_rows(directed):
            ctx=iv if directed else c;U=ut if directed else up;M=meant if directed else meanp
            z=self.Jet.variable(ctx,ctx.mpf(Z),2);de=o.core.delta_iv if directed else o.core.delta
            r=rho if directed else ctx.mpf(4);xx=x if directed else ctx.mpf(2)
            state=rt if directed else cp;ch=chii if directed else chi
            b=chider if directed else c.mpf(0)
            if not directed and q is not None and c.mpf('.5')<c.mpf(q)<c.mpf('.75'):
                arg=4*c.mpf(q)-2;sigpoint=stable_sigma_point(c,arg)
                b=-4*sigpoint*(1-sigpoint)*(2/arg**3+2/(1-arg)**3)
            fr0=state['fx']*ch/(2*xx);ur0=state['ux']*ch/(2*xx)
            frs=state['f']*b/r;urs=state['u']*b/r
            mr=(U-M)/r
            qr0=((ur0*(2*z)-mr*(z*(1+de))-o._derivative(mr)*(1-z*z))/(1-z*z*de)).truncate(1)
            qrs=(urs*(2*z)/(1-z*z*de)).truncate(1)
            return dict(f_rho={(0,0):fr0,(-1,-1):frs},u_rho={(0,0):ur0,(-1,-1):urs},
                mean_u_rho={(0,0):mr},q_rho={(0,0):qr0,(-1,-1):qrs})
        pr,tr=radial_rows(False),radial_rows(True)
        radial={name:{str(a)+','+str(b):o._jet_summary(pr[name][(a,b)],jet) for (a,b),jet in rows.items()} for name,rows in tr.items()}
        return dict(Z=Z,coordinate=coordinate,chart='selected_Rin_to_rminus_zero_velocity_own_history' if exterior else 'selected_original_inner_cutoff',
            source_family=self.family,normalized_original_n1=out,
            cutoff=dict(value=o._summary(chi,chii),keep_fraction='.5',zero_fraction='.75',
                definition='1-sigma(4q-2)',only_f_and_u_multiplied=True),
            geometry=dict(exact_sc_mpf_tuple=self.sc_tuple,exact_phase_in_mpf_tuple=self.phase_tuple,
                actual_phase_definition='fraction*original_sc/4',
                rho='4exp(original_hb*actual_phase)',original_positive_hb_not_selected=True),
            normalized_own_primitive=dict(definition='rho*mean_u; physical M1=epsilon^2 times this primitive',
                core=o._jet_summary(primitive_p,ct['mean_u']*4),
                increment_per_original_hb_phase_in=o._jet_summary(hp,ht),
                actual_primitive=o._jet_summary(primitive_p,primitive_t),
                exact_same_terminal_history_used=exterior or c.mpf(q)>=c.mpf('.75')),
            cutoff_control_mass=o._summary(mass,massi),native_first_radial_width_sectors=radial,
            width_sector_definition='key a,b denotes original_hb^a*original_phase_in^b times coefficient',
            pressure_history=dict(core_pressure_retained=True,leading_N1p_retained_once=True,
                cut_F_coupling_integrated_only_over_actual_support=True,
                known_regular_error=o._summary(c.mpf(0),o._symmetric(known_regular)),
                known_singular_variation_error=o._summary(c.mpf(0),o._symmetric(known_change)),
                cut_F_coupling_error=o._summary(c.mpf(0),o._symmetric(coupling_error))),
            raw_uniform_own_domain=dict(phase_upper='original_sc/4',bound_errors={n:o._summary(c.mpf(0),o._symmetric(e)) for n,e in errors.items()},
                only_uses_existing_admitted_raw_half_phase=True),
            coefficient_convention='axial coefficient k=D_Z^k/k!; radial derivatives ordinary',
            limits=dict(full_actual_selected_cutoff_numeric=True,raw_uncut_domain_extended=False,
                unrepaired_own_history_continuation_only=True,complete_five_moment_repair=False,
                exterior_m13_m15_pressure_moment_tails_closed=False,physical_velocity_assembled=False,
                complete_Cartesian_C2_join=False,global_field=False,temporal_n_ge_2_recursion=False,full_corrected_NS=False))
