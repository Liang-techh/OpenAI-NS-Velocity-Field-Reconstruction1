"""Signed core-to-cutoff/lower m15 with original width and flat atoms.

The singular leading pressure source is integrated as d-sigma, after exact
inverse-hb cancellation. Regular pressure and cut-velocity sources retain
their physical hb measure. Cumulative core pressure is never subtracted.
"""
import copy
from fractions import Fraction
import hashlib

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_collar_transport import Laurent, WidthChange
from .original_n1_cutoff_numeric import OriginalN1CutoffNumeric, stable_sigma_point
from .original_n1_core_m15_numeric import OriginalN1CoreM15Numeric


def leading_regular_frame(values, comparison, z, rho, x, de, eps, P2, sigma, hb, axis_j):
    """Actual prescribed-shear radial rows, with the separate hb^-1 row removed.

    Only the original sigma-prime/hb term is omitted here. Its weighted
    pressure contribution is evaluated by the d-sigma primitive below.
    """
    chi=1-sigma+sigma*hb
    A=comparison['Phi'][1]/comparison['Phi'][0]
    Ar=comparison['Phi'][2]/comparison['Phi'][0]-A*A
    ratio=values['Phi']/comparison['Phi'][0]
    phi=[Laurent({0:values['Phi']}),Laurent({0:values['Phi']*chi*A}),
         Laurent({0:values['Phi']*(chi*chi*A*A+chi*Ar)})]
    pr=ratio*chi*comparison['Psi'][1]
    psi=[Laurent({0:values['Psi']}),Laurent({0:pr}),
         Laurent({0:ratio*(chi*(chi-1)*A*comparison['Psi'][1]+chi*comparison['Psi'][2])})]
    mr=(values['Psi']-values['Mpsi'])/rho
    mean=[Laurent({0:values['Mpsi']}),Laurent({0:mr}),Laurent({0:(pr-2*mr)/rho})]
    U0=4*z+axis_j
    U=[psi[0]*eps+U0,psi[1]*eps,psi[2]*eps]
    M=[mean[0]*eps+U0,mean[1]*eps,mean[2]*eps]
    d,L=1-z*z,1-z*z*de
    Q=[(u*(2*z)-m*(z*(1-de))-m.dz()*d)/L for u,m in zip(U,M)]
    return dict(Phi=phi,U=U,M=M,Q=Q,z=z,rho=rho,x=x,de=de,e=eps,P2=P2)


def normalized_pressure_source(frame):
    """Original first-collar N1p/Pstar^2, as a finite hb Laurent jet."""
    z,r,de,e,P2=(frame[k] for k in ('z','rho','de','e','P2'))
    d,L=1-z*z,1-z*z*de
    q,qr,qrr=frame['Q'];qrho=qr*r
    flow=(q+q.dz()*(z*((1-de)/2))+qrho)/L+q*(q/2+qrho)
    flow+=frame['U'][0]*(q*(z*(-2))+q.dz()*d-qrho*(z*2))/L
    return -flow/(2*P2)+(qr*2+qrr*r)/(e*P2)


class OriginalN1InnerM15Numeric:
    def __init__(self, *, cells=32, cutoff=None, core=None, **kwargs):
        if type(cells) is not int or not 8<=cells<=256:
            raise ValueError('Use8..256 directed cells per fixed cutoff piece')
        if cutoff is not None and kwargs:raise ValueError('Supply a retained cutoff owner or constructor options')
        self.cut=cutoff if cutoff is not None else OriginalN1CutoffNumeric(**kwargs)
        self.o=self.cut.leading
        self.ctx,self.iv,self.J=self.cut.ctx,self.cut.iv,self.cut.Jet
        self.audit=mp.mp.clone();self.audit.dps=self.iv.dps+30
        self.family,self.directory=self.cut.family,self.o.core.directory
        self.core=core if core is not None else OriginalN1CoreM15Numeric(solver=self.cut.n1)
        if self.core.family!=self.family:raise ValueError('Same original signed core and cutoff source required')
        self.cells=cells
        self.hashes,self.code_hashes=dict(self.cut.hashes),dict(self.cut.code_hashes)
        for name,digest in self.core.hashes.items():
            if self.hashes.get(name,digest)!=digest:raise ValueError('Core and cutoff m15 source family changed')
            self.hashes[name]=digest
        self.code_hashes.update(self.core.code_hashes)
        name='src/openai_ns_reconstruction/original_n1_inner_m15_numeric.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
        self._prepared,self._integrals,self._packets,self._source_cells={},{},{},{}

    def rows(self,point,truth):
        return self.o._jet_summary(point.truncate(1),truth.truncate(1))

    def atom(self,point,truth):return self.o._summary(point,truth)

    def _prepare(self,Z):
        if Z in self._prepared:return self._prepared[Z]
        c,iv,J=self.ctx,self.iv,self.J
        raw,cp,ct,rt,errors=self.cut._raw_uniform(Z)
        tmax=_bounds(c,self.cut.phase_iv*2)[1]
        v=self.cut._inputs(Z,tmax,order=5);v['Z']=Z
        # Point fields approximate the same original zero-width inlet.
        inp=v['inlet']['point'];z=J.variable(c,c.mpf(Z),5)
        de,e,P2=self.o.core.delta,self.o.core.epsilon,c.exp(2*self.o.core.logP)
        vp=leading_regular_frame(inp['values'],inp['radial'],z,J.constant(c,4,5),
            J.constant(c,2,5),de,e,P2,J.constant(c,0,5),c.mpf(0),self.o.core.j)
        psource=normalized_pressure_source(vp).rows[0]
        fixed,_,_,_=self.o.core._seed(Z,24,5,directed=False)
        C=lambda base,error=0:WidthChange(self.o,base,error)
        zi=C(J.variable(iv,iv.mpf(Z),2));L=1-zi*zi*v['de']
        rr=C(J.constant(iv,4,2),4*iv.mpf(tmax)*iv.exp(self.o.hcap*iv.mpf(tmax)))
        one=C(J.constant(iv,1,2),iv.mpf(1))
        inlet=v['inlet']['truth']
        phi=C(inlet['values']['Phi'].truncate(2),v['delta']['Phi'])
        comp=C(inlet['radial']['Phi'][0].truncate(2),self.cut._core_comparison_changes('Phi',0,tmax,2))
        psir=C(inlet['radial']['Psi'][1].truncate(2),self.cut._core_comparison_changes('Psi',1,tmax,2))
        singular=2*rr*rr*zi*one*(phi/comp)*psir/L
        basep=inp['radial']['Psi'][1].truncate(2)*z.truncate(2)*32/(1-z.truncate(2)**2*de)
        result=dict(v=v,cp=cp,rt=rt,point_frame=vp,point_Np0=psource,
            point_S=J(c,fixed['S_Z_taylor'][:6]),singular_point=basep,
            singular_truth=singular.base,singular_width_error=singular.error,
            raw_domain='n1 raw fields used only where q<3/4 and t<phase_in')
        self._prepared[Z]=result
        return result

    def source_cell(self,Z,left,right):
        """Signed regular pressure/axial density over an exact q cell."""
        c,iv,J=self.ctx,self.iv,self.J
        key=Z,left._mpi_,right._mpi_
        if key in self._source_cells:return self._source_cells[key]
        a=self._prepare(Z);v=a['v']
        qbox=iv.mpf([_bounds(self.audit,left)[0],_bounds(self.audit,right)[1]])
        mid=sum(_bounds(c,qbox))/2
        t=self.cut.phase_iv*qbox
        with mp.workdps(280):sig=self.cut.sigma_jets(iv,t)[0]
        vt=leading_regular_frame(v['actual'],v['comparison'],v['z'],v['rho'],v['x'],
            v['de'],v['e'],v['P2'],J.constant(iv,sig,5),iv.mpf([0,_bounds(c,self.o._width_error(iv.mpf(1)))[1]]),self.o.core.j_iv)
        Np=normalized_pressure_source(vt).rows[0].truncate(1)
        with mp.workdps(280):chi=1-self.cut.sigma_jets(iv,4*qbox-2)[0]
        chip=1-stable_sigma_point(c,4*mid-2)
        f=a['rt']['f']*chi;u=a['rt']['u']*chi
        fp=a['cp']['f']*chip;up=a['cp']['u']*chip
        coupling=(v['S']*vt['Phi'][0].rows[0]*f*(2/(v['e']**3*v['P2']))).truncate(1)
        couplingp=(a['point_S']*a['point_frame']['Phi'][0].rows[0]*fp*(2/(self.o.core.epsilon**3*c.exp(2*self.o.core.logP)))).truncate(1)
        axial=(vt['rho']*vt['U'][0].rows[0]*u*2).truncate(1)
        axialp=(a['point_frame']['U'][0].rows[0]*up*8).truncate(1)
        pressure=(-(Np+coupling)*vt['rho']*vt['rho']).truncate(1)
        pressurep=(-(a['point_Np0'].truncate(1)+couplingp)*16).truncate(1)
        out=(axialp,pressurep),(axial,pressure)
        self._source_cells[key]=out
        return out

    def regular_integrals(self,Z,q):
        key=Z,str(q)
        if key in self._integrals:return self._integrals[key]
        c,iv,J=self.ctx,self.iv,self.J
        point=[J.constant(c,0,1) for _ in range(2)]
        truth=[J.constant(iv,0,1) for _ in range(2)]
        ends=[iv.mpf(0)]+[iv.mpf(v) for v in ('.5','.75','1','2') if Fraction(v)<Fraction(q)]+[iv.mpf(q)]
        for lo,hi in zip(ends,ends[1:]):
            if _bounds(self.audit,hi)[1]<=_bounds(self.audit,lo)[0]:continue
            width=(hi-lo)/self.cells
            for k in range(self.cells):
                left,right=lo+width*k,lo+width*(k+1)
                p,t=self.source_cell(Z,left,right)
                wp=sum(_bounds(c,width))/2
                for j in range(2):point[j]+=p[j]*wp;truth[j]+=t[j]*width
        self._integrals[key]=point,truth
        return point,truth

    def sigma_atom(self,q):
        """Exact original flat value as exp(odds)*kappa; never exp(odds)."""
        c,iv=self.ctx,self.iv
        t=self.cut.phase_iv*iv.mpf(q)
        odds=1/(1-t)**2-1/t**2
        cap=t*iv.exp(-1000)
        if _bounds(c,odds)[1]>_bounds(c,iv.log(t)-1000)[0]:
            raise ArithmeticError('Original flat atom is outside the certified vanishing cap')
        kappa=1/(1+iv.mpf([0,_bounds(c,cap)[1]]))
        tp=c.make_mpf(tuple(self.cut.phase_tuple))*c.mpf(q)
        return self.atom(1/(1-tp)**2-1/tp**2,odds),kappa

    def evaluate(self,Z,q):
        c,iv,J=self.ctx,self.iv,self.J
        if (type(Z) is not str or not c.isfinite(c.mpf(Z)) or not -1<c.mpf(Z)<1
                or type(q) is not str or not c.isfinite(c.mpf(q)) or not 0<=Fraction(q)<=2):
            raise ValueError('Original decimal Z(-1,1) and first/lower fraction q[0,2] required')
        key=Z,q
        if key in self._packets:return self._packets[key]
        inlet=self.core.evaluate(Z,'4');terms=copy.deepcopy(inlet['m15_core_C0_C1_terms'])
        pp,tt=self.regular_integrals(Z,q)
        e=self.o.core
        ep=self.atom(-2*e.logLambda,-2*e.ll_iv)
        P=self.atom(2*e.logP,2*e.lp_iv)
        h=self.atom(sum(_bounds(c,self.o.logh))/2,self.o.logh)
        h['source_log_atom_only']=True
        phase=self.atom(c.log(c.make_mpf(tuple(self.cut.phase_tuple))),iv.log(self.cut.phase_iv))
        def append(name,p,t,atoms,binding):
            terms.append(dict(name=name,coefficient_C0_C1=self.rows(p,t),
                log_scale_terms=copy.deepcopy(atoms),actual_function_binding=binding))
        append('first_lower_axial_width',pp[0],tt[0],[ep,h,phase],
            dict(density='2*rho*U_lead*u_cut',measure='hb*(sc/4)*dq',
                 full_leading_axis_and_Psi=True,own_cutoff='1-sigma(4q-2)'))
        append('first_lower_regular_pressure_width',pp[1],tt[1],[ep,P,h,phase],
            dict(density='-rho^2*(Np0+2S*Phi*f_cut/(epsilon^3*Pstar^2))',
                 leading_N1p_after_velocity_support_retained=True,constant_core_pressure_cancels_exactly=True))
        if Fraction(q)==0:
            zero=J.constant(c,0,1);zeroi=J.constant(iv,0,1)
            append('first_lower_singular_pressure',zero,zeroi,[ep],dict(empty_dsigma_prefix=True))
            append('first_lower_singular_width_error',zero,zeroi,[ep,h],dict(empty_dsigma_prefix=True))
        else:
            a=self._prepare(Z);atom,kappa=self.sigma_atom(q)
            append('first_lower_singular_pressure',a['singular_point'],a['singular_truth']*kappa,[ep,atom],
                dict(density_after_inverse_hb_cancellation='A(s,h)*sigma_prime(s) ds',
                     A0='32Z*Psi_core_rho(4)/L',actual_phase='s=q*(sc/4)',
                     original_nonzero_flat_amplitude_retained_as_log_atom=True))
            err=self.o._norm_error_jet(a['singular_width_error'],2)*kappa
            append('first_lower_singular_width_error',J.constant(c,0,2),err,[ep,h,atom],
                dict(proof='|A-A0|_norm<=hb*E, integrate against positive d_sigma',
                     accumulated_pressure_history_not_reset=True))
        out=dict(Z=Z,q=q,chart='selected_first_cutoff' if Fraction(q)<=1 else 'selected_lower_collar',
            source_family=self.family,m15_C0_C1_terms=terms,core_inlet=inlet,cells_per_piece=self.cells,
            exact_sc_tuple=self.cut.sc_tuple,exact_phase_tuple=self.cut.phase_tuple,
            original_geometry=dict(rho='4exp(hb*(sc/4)*q)',R='epsilon*rho',measure='d_rho=hb*(sc/4)*rho*dq'),
            contract=dict(actual_signed_core_inlet_retained_once=True,
                original_leading_sigma_phase_not_own_cutoff=True,
                leading_N1p_and_cutF_pressure_source_retained=True,
                hb_inverse_cancelled_before_enclosure=True,
                singular_pressure_integrated_against_dsigma=True,
                constant_pressure_cancelled_without_endpoint_subtraction=True,
                regular_cell_ranges_include_field_and_quadrature_errors=True,
                original_nonzero_flat_sigma_log_atom_retained=True,
                actual_first_lower_source_and_native_phase_preserved=True,
                uncut_n1_field_not_extended_beyond_admitted_domain=True,
                ordinary_Z0_Z1_and_product_rules_retained=True),
            limits=dict(signed_core_first_lower_m15_with_C1_error=True,
                active_compact_signed_pressure=False,independent_preRp_signed_numeric_seam=False,
                global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
                stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out
        return out
