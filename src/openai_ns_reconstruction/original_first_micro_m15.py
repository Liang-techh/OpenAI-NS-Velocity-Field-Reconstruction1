"""Signed first-micro m15 pressure integral with original finite-N errors.

The leading singular source is integrated against positive d-sigma. The
actual modified-minus-leading correction uses the reduced S1 identity,
including quadratic cross errors and the finite R*Q_y endpoint.
"""
import copy
from fractions import Fraction
import hashlib
import json
import math

import mpmath as mp

from .original_core_numeric import ROOT, PREFIX, _bounds
from .original_first_collar_transport import WidthChange
from .original_first_micro_pressure_rows import OriginalFirstMicroPressureRows
from .original_n1_cutoff_numeric import stable_sigma_point
from .original_n1_inner_m15_numeric import leading_regular_frame,normalized_pressure_source,OriginalN1InnerM15Numeric


class OriginalFirstMicroM15:
    def __init__(self,*,cells=8,rows=None,**kwargs):
        if type(cells) is not int or not 4<=cells<=128:raise ValueError('Use4..128 directed cells')
        if rows is not None and kwargs:raise ValueError('Supply retained source rows or constructor options')
        self.rows=rows if rows is not None else OriginalFirstMicroPressureRows(**kwargs)
        self.t,self.o=self.rows.t,self.rows.o
        self.ctx,self.iv,self.J,self.audit=self.rows.ctx,self.rows.iv,self.rows.J,self.rows.audit
        self.directory,self.family=self.rows.directory,self.rows.family
        self.hashes,self.code_hashes=dict(self.rows.hashes),dict(self.rows.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_n1_inner_m15_numeric.json')
        data=json.loads(path.read_bytes())
        if (data['source_family']!=self.family or not data['checks']['passed']
                or not data['limits']['signed_core_first_lower_m15_C0_C1']):
            raise ValueError('Accepted original signed cumulative rminus source required')
        for root,group,target in ((self.directory,data['input_hashes'],self.hashes),
                                 (ROOT,data['implementation_hashes'],self.code_hashes)):
            for name,h in group.items():
                if target.get(name,h)!=h or hashlib.sha256((root/name).read_bytes()).hexdigest()!=h:
                    raise ValueError('Changed or conflicting original m15 source: '+name)
                target[name]=h
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self.inlets={Fraction(p['Z']):p for p in data['samples'] if Fraction(p['q'])==2 and p['cells_per_piece']==16}
        self._inner_engine=None;self.cells=cells
        self.s0=self.t.phase_iv*2
        # A local first-cell source bound must also cover every volume prefix.
        for frame in self.rows.micro.ranges:
            chart=next(v for v in frame['actual_17_chart_C3_transports'] if v['chart']=='first_micro')
            if chart['cells'][0]['exact_left']!=['sc_half'] or chart['cells'][0]['exact_right']!=[1,4]:
                raise ValueError('Original first source cell must cover [sc/2,1/4]')
        self._prepared,self._integrals,self._packets={},{},{}
        name='src/openai_ns_reconstruction/original_first_micro_m15.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def inlet(self,Z):
        key=Fraction(Z)
        if key not in self.inlets:
            if self._inner_engine is None:self._inner_engine=OriginalN1InnerM15Numeric(cells=16,cutoff=self.t)
            self.inlets[key]=self._inner_engine.evaluate(Z,'2')
            self.hashes.update(self._inner_engine.hashes);self.code_hashes.update(self._inner_engine.code_hashes)
        p=self.inlets[key]
        if (p['source_family']!=self.family or Fraction(p['Z'])!=key or Fraction(p['q'])!=2
                or p['exact_sc_tuple']!=self.t.sc_tuple or p['exact_phase_tuple']!=self.t.phase_tuple
                or len(p['m15_C0_C1_terms'])!=6):
            raise ValueError('Complete same-Z six-sector signed rminus inlet required')
        return p

    def prepare(self,Z,s):
        key=Z,s._mpi_
        if key in self._prepared:return self._prepared[key]
        c,iv,J=self.ctx,self.iv,self.J;upper=_bounds(self.audit,s)[1]
        v=self.t._inputs(Z,upper,order=5);inp=v['inlet']['point']
        C=lambda b,error=0:WidthChange(self.o,b,error)
        z=C(J.variable(iv,iv.mpf(Z),2));L=1-z*z*v['de']
        r=C(J.constant(iv,4,2),4*iv.mpf(upper)*iv.exp(self.o.hcap*iv.mpf(upper)))
        one=C(J.constant(iv,1,2),iv.mpf(1));inc=v['inlet']['truth']
        phi=C(inc['values']['Phi'].truncate(2),v['delta']['Phi'])
        comp=C(inc['radial']['Phi'][0].truncate(2),self.t._core_comparison_changes('Phi',0,upper,2))
        psir=C(inc['radial']['Psi'][1].truncate(2),self.t._core_comparison_changes('Psi',1,upper,2))
        singular=2*r*r*z*one*(phi/comp)*psir/L
        zp=J.variable(c,c.mpf(Z),2)
        ap=inp['radial']['Psi'][1].truncate(2)*zp*32/(1-zp*zp*self.o.core.delta)
        out=dict(v=v,inlet_point=inp,singular_point=ap,singular_truth=singular.base,width_error=singular.error)
        self._prepared[key]=out
        return out

    def frame(self,p,Z,sbox,*,point):
        c,iv,J=self.ctx,self.iv,self.J;v=p['v']
        if point:
            s=sum(_bounds(c,sbox))/2;inp=p['inlet_point']
            return leading_regular_frame(inp['values'],inp['radial'],J.variable(c,c.mpf(Z),5),
                J.constant(c,4,5),J.constant(c,2,5),self.o.core.delta,self.o.core.epsilon,
                c.exp(2*self.o.core.logP),J.constant(c,stable_sigma_point(c,s),5),c.mpf(0),self.o.core.j)
        with mp.workdps(280):sigma=self.t.sigma_jets(iv,sbox)[0]
        return leading_regular_frame(v['actual'],v['comparison'],v['z'],v['rho'],v['x'],
            v['de'],v['e'],v['P2'],J.constant(iv,sigma,5),
            iv.mpf([0,_bounds(self.audit,self.o._width_error(iv.mpf(1)))[1]]),self.o.core.j_iv)

    def regular_integral(self,Z,s):
        key=Z,s._mpi_
        if key in self._integrals:return self._integrals[key]
        c,iv,J=self.ctx,self.iv,self.J;p=self.prepare(Z,s)
        value=J.constant(c,0,1);truth=J.constant(iv,0,1)
        width=(s-self.s0)/self.cells
        for k in range(self.cells):
            left,right=self.s0+width*k,self.s0+width*(k+1)
            box=iv.mpf([_bounds(self.audit,left)[0],_bounds(self.audit,right)[1]])
            fp,ft=self.frame(p,Z,box,point=True),self.frame(p,Z,box,point=False)
            dp=-(normalized_pressure_source(fp).rows[0]*fp['rho']*fp['rho']).truncate(1)
            dt=-(normalized_pressure_source(ft).rows[0]*ft['rho']*ft['rho']).truncate(1)
            value+=dp*(sum(_bounds(c,width))/2);truth+=dt*width
        self._integrals[key]=value,truth
        return value,truth

    def sigma_atom(self,s):
        """sigma=exp(odds)*kappa, including unmaterializable lower amplitude."""
        c,iv=self.ctx,self.iv
        odds=1/(1-s)**2-1/s**2
        if _bounds(self.audit,odds)[1]<=_bounds(self.audit,iv.log(s)-1000)[0]:
            cap=s*iv.exp(-1000)
            kappa=1/(1+iv.mpf([0,_bounds(self.audit,cap)[1]]));point=c.mpf(1)
        else:
            kappa=1/(1+iv.exp(odds));point=sum(_bounds(c,kappa))/2
        atom=self.o._summary(sum(_bounds(c,odds))/2,odds)
        atom['original_nonzero_sigma_log_atom']=True
        return atom,point,kappa

    def magnitude(self,jet,order=2):
        out=[]
        for k in range(order+1):
            cap=max(abs(v) for v in _bounds(self.audit,jet[k]*math.factorial(k)))
            out.append(None if not cap else self.iv.log(self.iv.mpf(cap)))
        return out

    def product(self,a,b):
        return [self.rows.sum_logs([None if x is None or y is None else x+y
            for x,y in ((a[j],b[k-j]) for j in range(k+1))]) for k in (0,1)]

    def delta_product(self,a,b,ea,eb):
        pieces=(self.product(a,eb),self.product(b,ea),self.product(ea,eb))
        return [self.rows.sum_logs(v[k] for v in pieces) for k in (0,1)]

    def finite_N_errors(self,Z,s,packet):
        """S1/2 boundary and volume errors, and separate viscous endpoint."""
        iv,J=self.iv,self.J;_,base,rho=self.rows.background(Z,s)
        bounds={name:self.magnitude(base[name]) for name in ('U','Q')}
        errors={name:[self.rows.read_log(v) for v in packet['physical_pressure_rows'][name]['finite_N_absolute_error_log_bounds']]
                for name in ('U','Q','Q_y')}
        u,q=bounds['U'],bounds['Q'];eu,eq=errors['U'],errors['Q']
        density=dict(q=eq[:2],qZ=eq[1:],q2=self.delta_product(q,q,eq,eq),
            U2=self.delta_product(u,u,eu,eu),UqZ=self.delta_product(u,q[1:],eu,eq[1:]),
            UUZ=self.delta_product(u,u[1:],eu,eu[1:]))
        z=J.variable(iv,iv.mpf(Z),2);one=J.constant(iv,1,2)
        de=self.o.core.delta_iv;d=1-z*z;L=1-z*z*de
        bc=dict(q=one/(2*L),q2=one/4,U2=-z*z/(L*L))
        vc=dict(q=-one/(2*L),qZ=z*(1-de)/(4*L),q2=-one/4,
            UqZ=d/(2*L),U2=z*z*(1-de)/(L*L),UUZ=z*d/(L*L))
        lr=self.magnitude(rho,0)[0]
        boundary=[None if v is None else v+2*lr for v in self.rows.combine(bc,density,1)]
        volume=[None if v is None else v+2*lr for v in self.rows.combine(vc,density,1)]
        viscous=[None if v is None else v+lr for v in errors['Q_y']]
        return dict(endpoint=boundary,volume=volume,viscous=viscous)

    def evaluate(self,Z,coordinate,*,selected_sc_multiple=False):
        packet=self.rows.evaluate(Z,coordinate,selected_sc_multiple=selected_sc_multiple)
        s,flat,recipe=self.rows.micro._coordinate(coordinate,selected_sc_multiple)
        key=Z,coordinate,selected_sc_multiple
        if key in self._packets:return self._packets[key]
        c,iv,J=self.ctx,self.iv,self.J;inlet=self.inlet(Z)
        terms=copy.deepcopy(inlet['m15_C0_C1_terms']);empty=s._mpi_==self.s0._mpi_
        epsilon=self.o._summary(-self.o.core.logLambda,-self.o.core.ll_iv)
        eps2=self.o._summary(-2*self.o.core.logLambda,-2*self.o.core.ll_iv)
        P2=self.o._summary(2*self.o.core.logP,2*self.o.core.lp_iv)
        h=self.o._summary(sum(_bounds(c,self.o.logh))/2,self.o.logh);h['source_log_atom_only']=True
        zero,zeroi=J.constant(c,0,1),J.constant(iv,0,1)
        def append(name,p,t,atoms,binding):
            terms.append(dict(name=name,coefficient_C0_C1=self.o._jet_summary(p.truncate(1),t.truncate(1)),
                log_scale_terms=copy.deepcopy(atoms),actual_function_binding=binding))
        binding=dict(interval='actual first_micro [sc/2,S]',measure='dy=hb ds; R=epsilon*rho',
            own_cut_velocities_exact_zero=True,positive_order_pressure_history_not_reset=True)
        if empty:
            leading=[('first_micro_regular_pressure',[eps2,P2,h]),('first_micro_singular_upper',[eps2]),
                     ('first_micro_singular_lower',[eps2]),('first_micro_singular_width_error',[eps2,h])]
            for name,atoms in leading:append(name,zero,zeroi,atoms,dict(binding,exact_empty_increment=True))
            error={name:[None,None] for name in ('endpoint','volume','viscous')}
        else:
            prepared=self.prepare(Z,s);pp,tt=self.regular_integral(Z,s)
            append('first_micro_regular_pressure',pp,tt,[eps2,P2,h],dict(binding,density='-rho^2*Np_regular ds'))
            for name,end,sign in (('first_micro_singular_upper',s,1),('first_micro_singular_lower',self.s0,-1)):
                atom,kp,kt=self.sigma_atom(end)
                append(name,prepared['singular_point']*(kp*sign),prepared['singular_truth']*(kt*sign),[eps2,atom],
                    dict(binding,primitive='A0*(sigma(S)-sigma(sc/2))',endpoint_sign=sign,
                         entire_one_minus_h_in_width_error=True,pressure_Pstar_squared_cancelled=True))
            atom,_,kt=self.sigma_atom(s)
            err=self.o._norm_error_jet(prepared['width_error'],2)*kt
            append('first_micro_singular_width_error',zero,err,[eps2,h,atom],
                dict(binding,proof='|A-A0|_norm<=hb*E; positive dsigma mass <=sigma(S)',
                     leading_inverse_hb_cancelled_before_evaluation=True))
            error=self.finite_N_errors(Z,s,packet)
        # At the defining rminus lower endpoint both correction value and
        # all needed derivative rows vanish by the source flat-collar theorem.
        lower=self.rows.evaluate(Z,'.5',selected_sc_multiple=True)
        if not all(v['exact_zero'] for name in ('U','Q','Q_y')
                for v in lower['physical_pressure_rows'][name]['finite_N_absolute_error_log_bounds']):
            raise ValueError('Actual flat lower finite-N pressure endpoint must be zero')
        for name,atoms in (('endpoint',[eps2]),('volume',[eps2,h]),('viscous',[epsilon])):
            logs=error[name];bound=self.rows.sum_logs(logs)
            if bound is None:truth=zeroi
            else:
                # One common norm bound covers independently differentiated
                # value/Z1 errors; it is not differentiated as a cap function.
                upper=iv.mpf(_bounds(self.audit,bound)[1])
                atom=self.o._summary(_bounds(c,upper)[1],upper)
                atom.update(positive_error_bound_atom_only=True,source_value_selected=False)
                atoms=atoms+[atom]
                if name=='volume':
                    length=iv.log(s-self.s0)
                    atoms=atoms+[self.o._summary(sum(_bounds(c,length))/2,length)]
                truth=J(iv,[iv.mpf([-1,1]),iv.mpf([-1,1])])
            append('first_micro_finite_N_'+name+'_error',zero,truth,atoms,
                dict(binding,error_enclosure_only=True,ordinary_C0_C1_source_product_bounds=True,
                     quadratic_delta_delta_cross_terms_retained=True,lower_endpoint_exact_zero_checked=True,
                     S1_half_and_viscous_minus_retained=True,cap_not_selected_as_signed_function=True))
        out=dict(Z=Z,coordinate=recipe,source_family=self.family,m15_C0_C1_terms=terms,
            incoming_q2_m15_terms=copy.deepcopy(inlet['m15_C0_C1_terms']),cells_per_piece=self.cells,
            original_sc_tuple=self.t.sc_tuple,actual_selected_integer=copy.deepcopy(packet['selected_integer']),
            finite_N_error_logs_before_physical_measures={k:self.rows.log_rows(v) for k,v in error.items()},
            exact_empty_increment=empty,
            contract=dict(signed_core_cutoff_lower_inlet_consumed_once=True,actual_first_micro_pressure_integrated=True,
                original_nonzero_lower_sigma_endpoint_retained=True,leading_singular_dsigma_positive_measure=True,
                inverse_hb_cancelled_before_enclosure=True,finite_N_S1_endpoint_and_volume_errors_retained=True,
                finite_R_Qy_and_R_QyZ_viscous_errors_retained=True,all_quadratic_cross_errors_retained=True,
                actual_selected_N_and_phase_preserved=True,native_hb_ds_and_RdR_applied_once=True,
                ordinary_C0_C1_errors_not_cap_derivatives=True,unknown_signed_corrections_not_point_values=True),
            limits=dict(signed_first_micro_leading_m15_with_complete_factored_finite_N_errors=True,
                resolved_finite_N_signed_centers=False,useful_total_relative_accuracy=False,
                remaining_active_compact_signed_pressure=False,independent_preRp_signed_seam=False,
                global_five_moment_numeric_join=False,common_global_N=False,exact_heat_exterior=False,
                stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out
        return out
