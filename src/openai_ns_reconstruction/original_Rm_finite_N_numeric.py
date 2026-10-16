"""Original finite-N Rm patch values with directed signed-error transport.

Leading controls remain fixed. Magnitude caps bound unknown signed source
increments; they are never selected as field values or phase samples.
"""
import hashlib
import json
import math
from fractions import Fraction

import mpmath as mp

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_Rm_patch_numeric import OriginalRmPatchNumeric,encode


RATES=dict(m='1',h='1.5',k='1.5',e='1',p='0')


def transport_step(ctx,incoming,driver,width,rate):
    """One actual-correction Duhamel step; driver already includes1/N."""
    r=ctx.mpf(rate)
    decay=ctx.exp(-r*width)
    mass=width if r==0 else -ctx.expm1(-r*width)/r
    return incoming*decay+driver*mass,decay,mass


def canonical_defect_increment(ctx,J,Z,x,am,P,invP,history):
    """Change of the original canonical partial defects in raw moment units."""
    z=J.variable(ctx,ctx.mpf(Z),2);dm=history['m']
    # M=Mz/R=P*m, whereas d0=x*(M-4Z). Only the latter gains x.
    return [dm*(P*x),(history['k']*P-z*history['h']*4)*(x**ctx.mpf('1.5'))/am,
        history['h']*(x**ctx.mpf('1.5'))/am,
        (history['e']-z*dm*(8*invP))*x/(am*am),history['p']/(am*am)]


class OriginalRmFiniteNNumeric:
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None):
        self.leading=OriginalRmPatchNumeric(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        p=self.leading;self.ctx,self.iv,self.J=p.ctx,p.iv,p.J;self.family=p.family
        self.owner=p.reference.owner;self.hashes=dict(p.hashes);self.code_hashes=dict(p.code_hashes)
        path=p.directory/(PREFIX+'compliant_current_original_Rm_patch_numeric.json')
        saved=json.loads(path.read_text(encoding='utf-8'))
        if not saved['checks']['passed'] or saved['source_family']!=self.family:
            raise ValueError('Accepted actual original leading patch required')
        for name,digest in saved['implementation_hashes'].items():
            if hashlib.sha256((ROOT/name).read_bytes()).hexdigest()!=digest:
                raise ValueError('Changed original leading implementation: '+name)
            self.code_hashes[name]=digest
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self._leading={(v['Z'],v['native_x']):v for v in saved['samples']};self._packets={}
        self.tail_log=-64*p.defect.logP;self.tail=self.iv.exp(self.tail_log)
        for stem in ('compliant_current_original_C3_target_ranges',
                     'compliant_current_original_C3_source_targets',
                     'compliant_current_original_true_Rm_patch_finite_N_C5_functions'):
            path=p.directory/(PREFIX+stem+'.py');digest=hashlib.sha256(path.read_bytes()).hexdigest()
            if self.hashes.get(path.name,digest)!=digest:raise ValueError('Changed finite-N source: '+stem)
            self.hashes[path.name]=digest
        name='src/openai_ns_reconstruction/original_Rm_finite_N_numeric.py'
        self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _jet(self,rows,directed):return self.leading._jet(rows,directed)
    def _rows(self,p,t):return self.leading.reference._rows(p,t)
    def _small(self,log,proof,label):
        if log is None:return self.iv.mpf(0)
        if _bounds(self.ctx,log)[1]>_bounds(self.ctx,self.tail_log)[0]:
            raise ArithmeticError('Original N does not admit patch numerical error scale: '+label)
        proof.append(dict(quantity=label,source_log_upper=self.owner._source_summary(log),
            admitted_bound_log=self.owner._source_summary(self.tail_log),cap_is_not_signed_field_value=True))
        return self.tail

    def _coordinate(self,x):
        c,iv=self.ctx,self.iv
        if type(x) is not str:raise ValueError('Native x decimal string or Rh required')
        if x=='Rh':return c.exp(1),iv.exp(1),c.mpf(1),iv.mpf(1)
        xp,xt=c.mpf(x),iv.mpf(x)
        if not c.isfinite(xp) or Fraction(x)<1 or _bounds(c,xt)[1]>_bounds(c,iv.exp(1))[0]:
            raise ValueError('Certified original1<=x<e decimal or exact Rh required')
        return xp,xt,c.ln(xp),iv.ln(xt)

    def _endpoint(self,value):
        return self.iv.exp(1) if value=='Rh' else self.iv.mpf(value[0])/value[1]

    def _cover(self,Z,xt):
        c,iv=self.ctx,self.iv;zl,zh=_bounds(c,iv.mpf(Z));xl,xh=_bounds(c,xt)
        frames=[];groups={};local=[]
        for frame in self.owner.ranges:
            left,right=map(c.mpf,frame['exact_Z_cell'])
            if right<zl or left>zh:continue
            if frame['source_identity']!=self.family or not frame['uniform_for_all_integer_N_ge_N0'] or frame['epsilon_domain']!='[0,1/N0]':
                raise ValueError('Same fixed leading family and original all-N bound required')
            frames.append(frame['exact_Z_cell'])
            chart=next(q for q in frame['actual_17_chart_C3_transports'] if q['chart']=='actual_patch')
            for cell in chart['cells']:
                a,b=self._endpoint(cell['exact_left']),self._endpoint(cell['exact_right'])
                al,ah=_bounds(c,a);bl,bh=_bounds(c,b)
                if not bh<xl and not al>xh:local.append(cell)
                key=json.dumps([cell['exact_left'],cell['exact_right']])
                groups.setdefault(key,dict(left=a,right=b,exact_left=cell['exact_left'],exact_right=cell['exact_right'],cells=[]))['cells'].append(cell)
        if not frames or not local:raise ValueError('Complete original active-patch source cover required')
        ordered=sorted(groups.values(),key=lambda g:_bounds(c,g['left'])[0])
        if ordered[0]['exact_left']!=[1,1] or ordered[-1]['exact_right']!='Rh' or any(
                a['exact_right']!=b['exact_left'] for a,b in zip(ordered,ordered[1:])):
            raise ValueError('Contiguous original native partition required')
        return frames,ordered,local

    def _log_cover(self,rows):
        # Every member bounds the same source on an overlapping cover.
        # Select a positive majorant, never a signed source value.
        logs=[self.owner._source_log(v) for v in rows];logs=[v for v in logs if v is not None]
        return None if not logs else self.iv.mpf(max(_bounds(self.ctx,v)[1] for v in logs))

    def evaluate(self,Z,x):
        key=(Z,x)
        if key in self._packets:return self._packets[key]
        c,iv,J=self.ctx,self.iv,self.J;xp,xt,yp,yt=self._coordinate(x)
        if key not in self._leading:self._leading[key]=self.leading.evaluate(Z,x)
        base=self._leading[key]
        if (Z,'1') not in self._leading:self._leading[Z,'1']=self.leading.evaluate(Z,'1')
        inlet=self._leading[Z,'1'];frames,groups,local=self._cover(Z,xt);proof=[]
        query=None if x=='Rh' else Fraction(x)
        incoming={n:self._jet(v,True) for n,v in inlet['finite_N_incoming_history_errors'].items()}
        history=dict(incoming);steps=[]
        # Merge coincident Z covers into one native cell. They must not
        # duplicate physical widths or decay the incoming memory twice.
        for group in groups:
            al,ah=_bounds(c,group['left']);bl,bh=_bounds(c,group['right'])
            left=Fraction(*group['exact_left'])
            right=None if group['exact_right']=='Rh' else Fraction(*group['exact_right'])
            if query is not None and query<=left:continue
            partial=query is not None and (right is None or query<right)
            end=xt if partial else group['right']
            width=iv.ln(end/group['left']);lo,hi=_bounds(c,width)
            if hi<=0:continue
            width=iv.mpf([max(c.mpf(0),lo),hi])
            weights={}
            for name,rate in RATES.items():
                rows=[]
                for k in range(3):
                    log=self._log_cover([cell['actual_density_ranges']['N_scaled'][name][k] for cell in group['cells']])
                    log=None if log is None else log-self.owner.logN
                    bound=self._small(log,proof,'density_'+name+'_Z'+str(k))/math.factorial(k)
                    rows.append(self.leading.reference.leading._symmetric(bound))
                history[name],decay,mass=transport_step(iv,history[name],J(iv,rows),width,rate)
                weights[name]=dict(decay=decay,mass=mass,rate=rate)
            steps.append(dict(exact_left=group['exact_left'],exact_right=group['exact_right'],
                actual_end=self.leading.reference._summary(min(xp,c.mpf(_bounds(c,group['right'])[1])),end),
                physical_log_width=width,own_weights=weights,overlapping_Z_covers=len(group['cells']),
                physical_measure='dx/x once',incoming_already_contains_one_over_N=True,new_driver_divided_by_N_once=True))
            if partial:break
        ab=[];bb=[]
        for k in range(3):
            for label,target in (('A',ab),('B',bb)):
                log=self._log_cover([cell['actual_source_and_primitive_range'][label][k] for cell in local])
                log=None if log is None else log-self.owner.logN
                target.append(self._small(log,proof,'local_'+label+'_over_N_Z'+str(k))/math.factorial(k))
        if _bounds(c,ab[0])[1]>=c.mpf('.25'):raise ArithmeticError('Original small exponent cover required')
        abs_exp=[2*ab[0],2*ab[1],2*ab[2]+ab[1]*ab[1]]
        exp_error=J(iv,[self.leading.reference.leading._symmetric(v) for v in abs_exp])
        deltaV=J(iv,[self.leading.reference.leading._symmetric(v) for v in bb])
        original={n:self._jet(v,True) for n,v in base['generic_leading_fields'].items()}
        point={n:self._jet(v,False) for n,v in base['generic_leading_fields'].items()}
        deltaE=original['E']*exp_error
        modified={n:v+history[n] for n,v in original.items() if n in RATES}
        modified.update(E=original['E']+deltaE,V=original['V']+deltaV,
            pressure=original['pressure']+history['p'])
        z=J.variable(iv,iv.mpf(Z),2);P=iv.exp(self.leading.defect.logP);invP=iv.exp(-self.leading.defect.logP)
        delta=iv.mpf(self.leading.reference.leading.core.delta_iv)
        dm=history['m'];dmZ=J(iv,[dm[1],2*dm[2]])
        dQ=((z*deltaV*2-z*dm*(1-delta)-(1-z*z)*dmZ)/(1-z*z*delta)*P).truncate(1)
        sim_point={n:self._jet(base['physical_similarity'][n],False) for n in
            ('Uz','Utheta_over_Pstar','Ur_over_sqrt_R_over2','Mz_over_R')}
        sim_truth={n:self._jet(base['physical_similarity'][n],True) for n in sim_point}
        sim_truth['Uz']+=deltaV*P;sim_truth['Utheta_over_Pstar']+=deltaE
        sim_truth['Ur_over_sqrt_R_over2']+=dQ;sim_truth['Mz_over_R']+=dm*P
        am=(1+z*z).reciprocal()*iv.exp(iv.mpf('-.6'))
        changes=canonical_defect_increment(iv,J,Z,xt,am,P,invP,history)
        defects=[]
        for i,rows in enumerate(base['actual_partial_defects']):
            defects.append(self._rows(self._jet(rows,False),self._jet(rows,True)+changes[i]))
        zero=J.constant(c,0,2)
        packet=dict(Z=Z,chart='actual_finite_N_Rm_patch',native_x=x,source_family=self.family,
            selected_integer=base['selected_integer'],leading_controls_Z2=base['controls_Z2'],
            generic_modified_fields={n:self._rows(point[n],v) for n,v in modified.items()},
            physical_similarity={n:self._rows(sim_point[n],sim_truth[n]) for n in sim_point},
            actual_canonical_partial_defect_enclosures=defects,
            finite_N_history_errors={n:self._rows(zero,v) for n,v in history.items()},
            finite_N_local_errors=dict(E=self._rows(zero,deltaE),V=self._rows(zero,deltaV)),
            finite_N_nonzero_incoming_errors=inlet['finite_N_incoming_history_errors'],
            physical_transport_steps=encode(steps),
            independent_P0_over_Pstar_squared=base['independent_P0_over_Pstar_squared'],
            positive_order_n1_own_history=base['positive_order_n1_own_history'],
            coordinate=dict(log_x=self.leading.reference._summary(yp,yt),Rm_reference_offset=-6,
                global_y_log_terms=[self.leading.reference._summary(10*self.leading.reference.logCp,10*self.leading.reference.logC),
                    self.leading.reference._summary(10*self.leading.defect.logPp,10*self.leading.defect.logP),
                    self.leading.reference._summary(c.mpf(-6),iv.mpf(-6)),self.leading.reference._summary(yp,yt)]),
            source_cover=dict(exact_Z_cells=frames,unique_native_cells=len(groups),local_source_cells=len(local),
                original_all_N_ge_N0=True,original_phase='frac(N*(logR-logRa-hb*sc/2))',
                full_phase_enclosure=True,phase_not_selected=True),
            small_error_proofs=proof,
            contract=dict(leading_controls_unchanged=True,primitive_B_is_B_over_Pstar=True,
                local_source='E_N=E*exp(A/N); V_N=V+B_over_Pstar/N',
                actual_N_scaled_density_contains_all_quadratic_terms=True,
                finite_incoming_already_contains_one_over_N=True,new_density_divided_by_N_once=True,
                pressure_memory_rate_zero=True,physical_Jacobian_dx_over_x_applied_once=True,
                signed_integral_centers_unresolved=True,leading_point_proxies_have_full_directed_error=True,
                magnitude_caps_not_signed_values=True,original_P0_unchanged_and_added_once=True,
                leading_terminal_zero_does_not_reset_finite_N_moments=True),
            limits=dict(local_and_own_history_Z_order=2,radial_velocity_Z_order=1,
                actual_finite_N_patch_values_enclosed=True,resolved_signed_finite_N_centers=False,
                total_native_y_or_fast_phase_derivatives=False,physical_Cartesian_assembly=False,
                modified_terminal_five_moment_closure=False,new_whole_axis_numeric_cover=False,
                global_N_admission=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=packet
        return packet
