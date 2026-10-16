"""Actual Rm patch pressure jets with genuine third mean/error rows."""
import copy
import hashlib
import json
import math
from fractions import Fraction

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_Rm_control_Z3 import OriginalRmControlZ3
from .original_Rm_finite_N_numeric import OriginalRmFiniteNNumeric,RATES
from .original_first_micro_pressure_rows import OriginalFirstMicroPressureRows,pressure_coefficients
from .original_reference_pressure_rows import OriginalReferencePressureRows


class OriginalRmPressureRows(OriginalFirstMicroPressureRows):
    def __init__(self, *, controls=None, **kwargs):
        if controls is not None and kwargs:raise ValueError('Supply genuine controls or constructor options')
        self.controls=controls if controls is not None else OriginalRmControlZ3(**kwargs)
        p=self.controls.patch;self.patch=p;self.reference=p.reference;self.owner=self.reference.owner
        self.inlet_rows=OriginalReferencePressureRows(reference=self.reference).next_field().next_field()
        for name in ('ctx','iv','J','audit','directory','family','o','t','micro'):
            setattr(self,name,getattr(self.inlet_rows,name))
        self.hashes,self.code_hashes=dict(self.controls.hashes),dict(self.controls.code_hashes)
        self.hashes.update(self.inlet_rows.hashes);self.code_hashes.update(self.inlet_rows.code_hashes)
        view=object.__new__(OriginalRmFiniteNNumeric)
        view.ctx,view.iv,view.J,view.owner,view.leading=self.ctx,self.iv,self.J,self.owner,p
        view.family=self.family
        self.cover_owner=view
        path=self.directory/(PREFIX+'compliant_current_generic_shear_loop_jet_bounds.json')
        self.loop=json.loads(path.read_bytes())['current_actual_loop_jet_log_bounds_by_chart']['actual_patch']
        if self.loop['source_family']!=self.family or self.loop['primitive_units']['B']!='original B/Pstar; E=Utheta/Pstar':
            raise ValueError('Original actual_patch slow/fast pressure derivative units required')
        A=self.read_log(self.loop['slow_phase_held_log_bounds']['A']['y0_Z0'])
        if A is not None and _bounds(self.audit,self.owner.logN)[0]<_bounds(self.audit,A)[1]:
            raise ValueError('Original N must admit patch A/N')
        path=self.directory/(PREFIX+'compliant_current_original_Rm_control_Z3_numeric.json')
        raw=path.read_bytes();data=json.loads(raw)
        if not data['checks']['passed'] or data['source_family']!=self.family:
            raise ValueError('Checked genuine Rm third controls required')
        for base,group,target in ((self.directory,data['input_hashes'],self.hashes),(ROOT,data['implementation_hashes'],self.code_hashes)):
            for name,digest in group.items():
                if hashlib.sha256((base/name).read_bytes()).hexdigest()!=digest or target.get(name,digest)!=digest:
                    raise ValueError('Changed Rm pressure prerequisite: '+name)
                target[name]=digest
        self.hashes[path.name]=hashlib.sha256(raw).hexdigest()
        self._backgrounds={};self._packets={};self._histories={};self._inlet_histories={};self._uniforms={}
        name='src/openai_ns_reconstruction/original_Rm_pressure_rows.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def geometry(self,x):
        c,iv,ref=self.ctx,self.iv,self.reference
        xp,xt,yp,yt=self.cover_owner._coordinate(x)
        terms=[(c.log(110),iv.log(110)),(10*ref.logCp,10*ref.logC),
               (10*ref.logPp,10*ref.logP),(c.mpf(-6),iv.mpf(-6)),(yp,yt)]
        return dict(x_point=xp,x=xt,y_point=yp,y=yt,atomic_log_radius_terms=terms,
            logR_point=sum((a for a,b in terms),c.mpf(0)),logR=sum((b for a,b in terms),iv.mpf(0)))

    def _fields(self,Z,ctx,xx,directed,terminal,uniform=False):
        p,J=self.patch,self.J;data=self.controls.coefficients(Z)
        h=data['truth' if directed else 'point']
        defects=data['source']['truth' if directed else 'point']['defects']
        z=J.variable(ctx,ctx.mpf(Z),3);zero=J.constant(ctx,0,3)
        if uniform:
            def whole(value):
                lo,hi=_bounds(self.audit,value);return ctx.mpf([min(0,lo),max(0,hi)])
            W={name:[[whole(v) for v in row] for row in values] if name=='L' else [whole(v) for v in values]
                for name,values in p.W.items()}
        else:W=p.W if directed and terminal else p._W(xx,directed,terminal)
        from lei_ren_part1_paper_compliant_actual_moment_patch import nonlinear_series
        invP=ctx.exp(-ctx.mpf(p.defect.logP if directed else p.defect.logPp))
        am=(1+z*z).reciprocal()*ctx.exp(ctx.mpf('-.6'));alpha=(am*am).reciprocal()*(invP*invP)
        Q=nonlinear_series(ctx,W,h,alpha)
        unrefined=[defects[i]+Q[i]+sum((h[j]*W['L'][i][j] for j in range(5)),zero) for i in range(5)]
        residual=[zero]*5 if terminal and not uniform else unrefined
        beta=[]
        for center in p.centers:
            cc=ctx.mpf(center.numerator)/center.denominator
            if uniform:
                beta.append((ctx.mpf([0,_bounds(self.audit,p.bump.beta_sup)[1]]),
                    ctx.mpf([-_bounds(self.audit,p.bump.beta_deriv_sup)[1],_bounds(self.audit,p.bump.beta_deriv_sup)[1]])))
            elif terminal:
                # Exact x>=71/40 is beyond every compact bump support,
                # including the right endpoint itself by flatness.
                beta.append((ctx.mpf(0),ctx.mpf(0)))
            elif directed:beta.append(p.bump.beta(xx-cc))
            else:
                u=(xx-cc)/ctx.mpf('.025')
                value=ctx.mpf(0) if abs(u)>=1 else ctx.exp(-1/(1-u*u))/(ctx.mpf('.025')*p.normp)
                derivative=ctx.mpf(0) if not value else -2*u*value/(ctx.mpf('.025')*(1-u*u)**2)
                beta.append((value,derivative))
        g=sum((h[i]*beta[k][0] for i,k in enumerate((0,2))),zero)
        gy=sum((h[i]*beta[k][1]*xx for i,k in enumerate((0,2))),zero)
        f=sum((h[i+2]*beta[i][0] for i in range(3)),zero)
        fy=sum((h[i+2]*beta[i][1]*xx for i in range(3)),zero)
        U=4*z+g;Uy=gy;M=4*z+residual[0]/xx
        de=self.o.core.delta_iv if directed else self.o.core.delta;A,D=1-z*z*de,1-z*z
        derivative=lambda row:J(ctx,[(k+1)*row[k+1] for k in range(row.order)])
        Q=((U*(2*z)-M*(z*(1-de))-derivative(M)*D)/A).truncate(2)
        Qy=((Uy*(2*z)-(U-M)*(z*(1-de))-derivative(U-M)*D)/A).truncate(2)
        if terminal and not uniform:
            # The unique full-weight branch has U=M=4Z and Uy=0.
            # Apply that identity before interval subtraction loses correlation.
            Qy=zero.truncate(2)
        return dict(U=U,m=M,Q=Q,U_y=Uy,Q_y=Qy,E=am*(f+xx**ctx.mpf('.1')),
            E_y=am*(fy+xx**ctx.mpf('.1')*ctx.mpf('.1')),partial_defects=residual,
            unrefined_partial_defects=unrefined)

    def background(self,Z,x):
        key=Z,x
        if key not in self._backgrounds:
            geo=self.geometry(x);terminal=x=='Rh' or Fraction(x)>=Fraction(71,40)
            self._backgrounds[key]=(self._fields(Z,self.ctx,geo['x_point'],False,terminal),
                self._fields(Z,self.iv,geo['x'],True,terminal),geo)
        return self._backgrounds[key]

    def uniform_background(self,Z,x):
        key=Z,x
        if key not in self._uniforms:
            geo=self.geometry(x);xx=self.iv.mpf([1,_bounds(self.audit,geo['x'])[1]])
            self._uniforms[key]=self._fields(Z,self.iv,xx,True,False,uniform=True)
        return self._uniforms[key]

    def incoming_history_logs(self,Z):
        if Z not in self._inlet_histories:
            local,prefix,_=self.reference._covers(Z,'postrestore',self.iv.mpf(1))
            _,_,hist=self.owner._correction_logs(local,prefix,False)
            if any(len(v)<4 for v in hist.values()):raise ValueError('All four genuine incoming C3 rows required')
            self._inlet_histories[Z]={name:list(v[:4]) for name,v in hist.items()}
        return self._inlet_histories[Z]

    def history_logs(self,Z,x):
        key=Z,x
        if key in self._histories:return self._histories[key]
        iv=self.iv;geo=self.geometry(x);frames,groups,local=self.cover_owner._cover(Z,geo['x'])
        hist={name:list(v) for name,v in self.incoming_history_logs(Z).items()};query=None if x=='Rh' else Fraction(x);steps=[]
        for group in groups:
            left=Fraction(*group['exact_left']);right=None if group['exact_right']=='Rh' else Fraction(*group['exact_right'])
            if query is not None and query<=left:continue
            partial=query is not None and (right is None or query<right)
            end=geo['x'] if partial else group['right'];width=iv.log(end/group['left'])
            lo,hi=_bounds(self.audit,width)
            if hi<=0:continue
            width=iv.mpf([max(0,lo),hi]);weights={}
            for name,rate in RATES.items():
                r=iv.mpf(rate);mass=width if r==0 else -iv.expm1(-r*width)/r
                updated=[]
                for k in range(4):
                    source=self.cover_owner._log_cover([cell['actual_density_ranges']['N_scaled'][name][k] for cell in group['cells']])
                    driver=None if source is None else source-self.owner.logN+iv.log(mass)
                    memory=None if hist[name][k] is None else hist[name][k]-r*width
                    updated.append(self.sum_logs((memory,driver)))
                hist[name]=updated;weights[name]=dict(rate=rate,decay_log=-r*width,positive_mass=mass)
            steps.append(dict(exact_left=group['exact_left'],exact_right=group['exact_right'],
                actual_log_width=width,overlapping_Z_covers=len(group['cells']),weights=weights))
            if partial:break
        out=dict(history=hist,frames=frames,groups=groups,local=local,steps=steps)
        self._histories[key]=out;return out

    def evaluate(self,Z,x):
        key=Z,x
        if key in self._packets:return self._packets[key]
        if type(Z) is not str or not self.ctx.isfinite(self.ctx.mpf(Z)) or not -1<self.ctx.mpf(Z)<1:
            raise ValueError('Interior decimal axial Z required')
        p,t,geo=self.background(Z,x);h=self.history_logs(Z,x);lp=self.o.core.lp_iv
        B=[self.cover_owner._log_cover([cell['actual_source_and_primitive_range']['B'][k] for cell in h['local']]) for k in range(4)]
        errors=dict(U=[None if v is None else v+lp-self.owner.logN for v in B],
            m=[None if v is None else v+lp for v in h['history']['m']])
        errors['U_Z'],errors['m_Z']=errors['U'][1:],errors['m'][1:]
        slow=self.loop['slow_phase_held_log_bounds']['B'];fast=self.loop['phase_and_first_slow_phase_log_bounds']
        slow_y=[self.read_log(slow[name]) for name in ('y1_Z0','y1_Z1')]
        phase_y=[self.read_log(fast[name]) for name in ('B_phi','B_Z_phi')]
        slow_y=[None if v is None else v+lp-self.owner.logN for v in slow_y]
        phase_y=[None if v is None else v+lp for v in phase_y]
        errors['U_y']=[self.sum_logs(values) for values in zip(slow_y,phase_y)]
        coeff=pressure_coefficients(self.J,self.iv,Z,self.o.core.delta_iv,2)
        errors['Q']=self.combine(coeff['Q'],errors,2);errors['Q_y']=self.combine(coeff['Q_y'],errors,1)
        fields={name:dict(signed_background_ordinary_Z=[self.o._summary(p[name][k]*math.factorial(k),t[name][k]*math.factorial(k)) for k in range(order+1)],
            finite_N_absolute_error_log_bounds=self.log_rows(errors[name][:order+1]),finite_N_log_bound_is_not_selected_signed_value=True)
            for name,order in (('U',2),('m',3),('Q',2),('U_y',1),('Q_y',1))}
        out=dict(Z=Z,native_x=x,chart='actual_patch',source_family=self.family,physical_pressure_rows=fields,
            selected_integer=copy.deepcopy(self.owner.integer),
            R_Q_y_endpoint=dict(signed_background_ordinary_Z=fields['Q_y']['signed_background_ordinary_Z'],
                log_physical_factor=self.o._summary(geo['logR_point'],geo['logR']),
                finite_N_coefficient_error_log_bounds=self.log_rows(errors['Q_y']),finite_N_has_no_extra_epsilon=True),
            geometry=dict(atomic_log_radius_terms=[self.o._summary(a,b) for a,b in geo['atomic_log_radius_terms']],
                native_x=self.o._summary(geo['x_point'],geo['x']),log_x=self.o._summary(geo['y_point'],geo['y']),
                log_physical_radius=self.o._summary(geo['logR_point'],geo['logR']),
                physical_measure='dy=dx/x once',physical_radius_materialized=False,exact_Rh_endpoint=x=='Rh'),
            leading_partial_defects=[self.reference._rows(a,b) for a,b in zip(p['partial_defects'],t['partial_defects'])],
            finite_N_history_error_logs={name:self.log_rows(v) for name,v in h['history'].items()},
            finite_N_incoming_history_error_logs={name:self.log_rows(v) for name,v in self.incoming_history_logs(Z).items()},
            physical_transport_steps=self.reference.previous._serializable_proof(h['steps']),
            finite_N_axial_radial_error_components=dict(slow_B_y_over_N=self.log_rows(slow_y),fast_B_phi=self.log_rows(phase_y),
                formula='deltaUy=Pstar*(By/N+Bphi); mixed Pstar*(ByZ/N+BZphi)',extra_hb_factor=False,extra_x_factor=False),
            source_cover=dict(exact_Z_cells=h['frames'],unique_native_cells=len(h['groups']),local_source_cells=len(h['local']),
                generic_derivative_chart='actual_patch',native_domain=['1','e'],original_phase='frac(N*(logR-logRa-hb*sc/2))',
                all_upstream_reference_histories_retained=True,overlapping_Z_cells_do_not_duplicate_width=True),
            contract=dict(genuine_original_controls_Z3=True,physical_defects_not_normalized_substitute=True,
                actual_partial_bump_weights=True,mean_Z3_and_QZZ_not_zero_padded=True,
                actual_Uy_and_Qy_Z1=True,all_four_finite_N_history_rows=True,both_incoming_and_local_errors=True,
                incoming_one_over_N_not_reapplied=True,each_new_driver_divided_by_N_once=True,
                original_Pstar_once=True,slow_fast_phase_terms_separate=True,fast_Bphi_not_divided_by_N=True,
                native_dx_over_x_once=True,ordinary_Z_factorial_once=True,
                leading_terminal_identity_does_not_reset_finite_N_history=True),
            limits=dict(actual_patch_UQ_Z2_mean_Z3_UyQy_Z1_with_complete_absolute_errors=True,
                active_patch_signed_m15_integrated=False,resolved_finite_N_signed_centers=False,
                useful_total_relative_accuracy=False,global_five_moment_numeric_join=False,common_global_N=False,
                exact_heat_exterior=False,stress_cone=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=out;return out
