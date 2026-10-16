"""Original Rc-to-2Rc compact background, with correlated terminal bounds.

Finite-N signed centers are unresolved. Zero proxies carry source-owned
directed errors; neither a magnitude endpoint nor a finite Picard iterate
defines the unique original compact control function.
"""
from fractions import Fraction
import gzip
import hashlib
import json
import math

import mpmath as mp

from .original_core_numeric import ROOT,PREFIX,_bounds
from .original_outer_o3_numeric import OriginalOuterO3Numeric
from .original_Rm_finite_N_numeric import RATES
from .original_Rm_patch_numeric import encode


def supported_mass(iv,s,L,rate,direction,supports):
    """Own-rate mass over true bump supports in physical log radius du."""
    r=iv.mpf(rate);c=mp.mp.clone();c.dps=iv.dps+10
    lo,hi=_bounds(c,s);ll,lh=_bounds(c,L);total=iv.mpf(0);parts=[]
    for a,b in supports:
        al,ah=_bounds(c,a);bl,bh=_bounds(c,b)
        if direction=='forward':
            if hi<=al:continue
            left=iv.mpf([max(c.mpf(0),al),max(c.mpf(0),ah)])
            right=iv.mpf([min(lo,bl),min(hi,bh)])
            if _bounds(c,right)[1]<=_bounds(c,left)[0]:continue
            width=right-left;wl,wh=_bounds(c,width);width=iv.mpf([max(0,wl),max(0,wh)])
            mass=width if r==0 else -iv.expm1(-r*width)/r
            mass*=iv.exp(-r*(s-right))
        elif direction=='backward':
            if lo>=bh:continue
            left=iv.mpf([max(lo,al),max(hi,ah)])
            right=b
            width=right-left;wl,wh=_bounds(c,width)
            if wh<=0:continue
            width=iv.mpf([max(0,wl),max(0,wh)])
            mass=width if r==0 else iv.expm1(r*width)/r
            mass*=iv.exp(r*(left-s))
        else:raise ValueError('Forward or backward own-rate mass required')
        # A mass is nonnegative even if directed subtraction straddles0.
        ml,mh=_bounds(c,mass);mass=iv.mpf([max(0,ml),max(0,mh)])
        total+=mass;parts.append(dict(support_left=a,support_right=b,clipped_log_width=width,mass=mass))
    return total,parts


class OriginalCompactNumeric:
    def __init__(self,*,precision=85,amplitude_precision=110,directory=None,outer=None):
        self.outer=outer or OriginalOuterO3Numeric(precision=precision,amplitude_precision=amplitude_precision,directory=directory)
        o=self.outer;self.patch=o.patch;self.ctx,self.iv,self.J=o.ctx,o.iv,o.J
        self.family=o.family;self.owner=o.owner;self.directory=o.directory
        self.hashes=dict(o.hashes);self.code_hashes=dict(o.code_hashes)
        path=self.directory/(PREFIX+'compliant_current_original_outer_o3_numeric.json')
        accepted=json.loads(path.read_text(encoding='utf-8'))
        if not accepted['checks']['passed'] or accepted['source_family']!=self.family:raise ValueError('Accepted actual original Rc numerical inlet required')
        for n,d in accepted['implementation_hashes'].items():
            if hashlib.sha256((ROOT/n).read_bytes()).hexdigest()!=d:raise ValueError('Changed original O3 implementation: '+n)
            self.code_hashes[n]=d
        self.hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        self._inlets={p['Z']:p for p in accepted['samples'] if p['native_O3_chart']=='O3_power' and p['native_O3']=='2'}
        self._buffer_inlets=self.outer._inlets;self._packets={}
        specifications=(('current_original_C3_limit_controls',('current_original_same_Banach_limit_C3_controls_installed','current_original_same_repair_integer_quantitative_C3_Picard_tail_installed')),
            ('current_original_C3_repair_band',('current_original_same_limit_repair_band_C3_velocity_histories_installed','current_original_same_limit_repair_band_C3_relative_terminal_identities_installed','current_original_same_limit_repair_band_C3_directed_magnitude_ranges_installed')),
            ('current_original_true_compact_C5_functions',('current_original_true_same_branch_compact_controls_C5_functions_installed','current_original_true_Rc_to_2Rc_relative_five_C5_zero_identities_installed')))
        raw={}
        for stem,gates in specifications:
            receipt_path=self.directory/(PREFIX+'compliant_'+stem+'_check.json')
            receipt=json.loads(receipt_path.read_text(encoding='utf-8'))
            if not receipt['all_passed'] or receipt['source_family']!=self.family or not all(receipt[k] for k in gates):raise ValueError('Accepted same-original compact source required: '+stem)
            for n,d in receipt['input_hashes'].items():
                if hashlib.sha256((self.directory/n).read_bytes()).hexdigest()!=d or self.hashes.get(n,d)!=d:raise ValueError('Changed/conflicting original compact input: '+n)
                self.hashes[n]=d
            self.hashes[receipt_path.name]=hashlib.sha256(receipt_path.read_bytes()).hexdigest()
            archive=self.directory/(PREFIX+'compliant_'+stem+'.json.gz')
            self.hashes[archive.name]=hashlib.sha256(archive.read_bytes()).hexdigest()
            raw[stem]=json.loads(gzip.decompress(archive.read_bytes()))
            if raw[stem]['source_family']!=self.family:raise ValueError('Same original compact family required')
        self.bounds=raw['current_original_C3_limit_controls']['actual_quantitative_C3_bounds']
        self.ranges=raw['current_original_C3_repair_band']['actual_four_Z_C3_band_ranges']
        compact=raw['current_original_true_compact_C5_functions'];self.binding=compact['true_compact_C5_bindings'];nodes=compact['exact_graph_nodes']
        selected=self.bounds['actual_same_repair_integer'];bandN=raw['current_original_C3_repair_band']['actual_selected_repair_integer']
        inletN=next(iter(self._inlets.values()))['selected_integer']
        if not selected['J']==bandN['J']==nodes[self.binding['N']]['J']==inletN['recipe']['J']:raise ValueError('Unchanged exact original selected N required')
        self.certificates={}
        for n,node in self.binding['terminal_certificates'].items():
            certificate=nodes[node]
            if (certificate['operation']!='proved_original_compact_relative_C5_terminal_identity'
                or not certificate['relative_correction_only'] or not certificate['absolute_leading_histories_and_pressure_not_zeroed']
                or certificate['source_system']!=self.binding['system']):raise ValueError('Original signed compact terminal identity required')
            self.certificates[n]=dict(source_archive=PREFIX+'compliant_current_original_true_compact_C5_functions.json.gz',node=node,
                operation=certificate['operation'],source_system=certificate['source_system'],ordinary_axial_orders=list(range(6)),
                relative_only=True,zero_is_source_function_identity_not_rounding=True)
        self.L=self.iv.log(2);self.centers=[self.L*self.iv.mpf(v) for v in ('.2','.5','.8')];ell=self.L/40
        self.supports=[(v-ell,v+ell) for v in self.centers]
        self.control_source=dict(selected_integer=selected,exact_control_function_nodes=self.binding['h'],
            exact_normalized_target_nodes=self.binding['targets'],joint_numerator_nodes=self.binding['joint'],
            original_unique_branch_system=self.binding['system'],branch_transfer=self.binding['branch_transfer'],
            control_bounds=self.bounds,point_centers_resolved=False,cap_endpoints_never_used_as_signed_values=True)
        name='src/openai_ns_reconstruction/original_compact_numeric.py';self.code_hashes[name]=hashlib.sha256((ROOT/name).read_bytes()).hexdigest()

    def _rows(self,p,t):return self.outer._rows(p,t)
    def _jet(self,rows,directed):return self.outer._jet(rows,directed)

    def _cap(self,records,proof,label,k):
        log=self.patch._log_cover(records)
        # These archived repair bounds already contain1/N. In particular,
        # never multiply a relaxed incoming cap by the unmaterialized N.
        cap=self.patch._small(log,proof,label)/math.factorial(k)
        return self.patch.leading.reference.leading._symmetric(cap)

    def _intersect(self,a,b):
        c,iv,J=self.ctx,self.iv,self.J;out=[]
        for u,v in zip(a.coefficients,b.coefficients):
            lo,hi=_bounds(c,u);ll,hh=_bounds(c,v)
            if max(lo,ll)>min(hi,hh):raise ArithmeticError('Original compact forward/backward bounds disagree')
            out.append(iv.mpf([max(lo,ll),min(hi,hh)]))
        return J(iv,out)

    def _scalar(self,inlet,s_point,s_directed):
        scalar=self.outer.scalar;out={}
        for label,ctx,mu,s in (('point',self.ctx,self.outer.mu,s_point),('directed',self.iv,self.outer.mui,s_directed)):
            restore=lambda row:self.patch.leading.reference.previous._restore(row) if label=='directed' else self.ctx.mpf(row['approximation'])
            old={n:restore(v) for n,v in inlet['normalized_mu_profiles'].items()};q=ctx.mpf('2.5')+s
            result=dict(old,J=q,U1=ctx.expm1(-mu*q)/mu,U2=ctx.expm1(-2*mu*q)/mu)
            for name,alpha,c in (('C1',1,1),('C0',0,2),('Cm',-1,2)):
                # (.5+u+2)^n; each original scalar moment has a directed
                #100-term exp-series remainder. Here q<=2.5+log2<3.2.
                moments=[sum((math.comb(n,k)*ctx.mpf(2)**(n-k)*scalar._moment(s,alpha,k,ctx) for k in range(n+1)),ctx.mpf(0)) for n in range(1,5)]
                local=-c*moments[0]+c*c*mu*moments[1]/2-c**3*mu**2*moments[2]/6
                if label=='directed':
                    rem=c**4*mu**3*moments[3]/24;local+=ctx.mpf([0,_bounds(self.ctx,rem)[1]])
                result[name]=old[name]+ctx.exp(3*alpha)*local
            out[label]=result
        return out

    def evaluate(self,Z,x):
        c,iv,J=self.ctx,self.iv,self.J
        if type(x) is not str or not c.isfinite(c.mpf(x)) or not 1<=Fraction(x)<=2:raise ValueError('Original compact native1<=x=R/Rc<=2 required')
        key=Z,x
        if key in self._packets:return self._packets[key]
        if Z not in self._inlets:self._inlets[Z]=self.outer.evaluate(Z,'O3_power','2')
        if Z not in self._buffer_inlets:self._buffer_inlets[Z]=self.outer.buffer.evaluate(Z,'11')
        inlet=self._inlets[Z];sp,st=c.log(c.mpf(x)),iv.log(iv.mpf(x))
        zl,zh=_bounds(c,iv.mpf(Z));frames=[p for p in self.ranges if not (c.mpf(p['exact_Z_cell'][1])<zl or c.mpf(p['exact_Z_cell'][0])>zh)]
        if not frames:raise ValueError('Original compact axial source cover required')
        scalar=self._scalar(inlet,sp,st)
        point,sectorsp,Qp=self.outer._leading(Z,c.mpf(2)+sp,'O3_power',self._buffer_inlets[Z],scalar,c)
        truth,sectorst,Qt=self.outer._leading(Z,iv.mpf(2)+st,'O3_power',self._buffer_inlets[Z],scalar,iv)
        proof=[];forward={};backward={};history={};masses={};zero=J.constant(c,0,2)
        for n,rate in RATES.items():
            density=J(iv,[self._cap([p['actual_signed_density_magnitude_C3_bounds'][n][k] for p in frames],proof,'compact_density_'+n+'_Z'+str(k),k) for k in range(3)])
            initial=self._jet(inlet['finite_N_history_errors'][n],True)
            incoming=J(iv,[self._cap([p['actual_incoming_correction_C3_bounds'][n][k] for p in frames],proof,'compact_incoming_'+n+'_Z'+str(k),k) for k in range(3)])
            initial=self._intersect(initial,incoming)
            supports=[self.supports[i] for i in ((0,2) if n in ('m','k') else (0,1,2))]
            fm,fp=supported_mass(iv,st,self.L,rate,'forward',supports)
            bm,bp=supported_mass(iv,st,self.L,rate,'backward',supports)
            forward[n]=initial*iv.exp(-iv.mpf(rate)*st)+density*fm
            backward[n]=density*bm
            history[n]=self._intersect(forward[n],backward[n])
            if x=='1':history[n]=self._intersect(self._jet(inlet['finite_N_history_errors'][n],True),history[n])
            masses[n]=dict(forward=fm,backward=bm,forward_parts=fp,backward_parts=bp,rate=rate,
                incoming_already_contains_one_over_N=True,density_already_contains_one_over_N=True,physical_measure='du=dlog(x), once')
        local=[]
        sl,sh=_bounds(c,st)
        for i,(a,b) in enumerate(self.supports):
            if not (sh<_bounds(c,a)[0] or sl>_bounds(c,b)[1]):local.append(i)
        corrections={}
        for name,label in (('E','delta_E'),('V','delta_V')):
            active=bool(local) if name=='E' else bool(set(local)&{0,2})
            corrections[name]=J(iv,[self._cap([p['actual_selected_N_profile_y0_y1_y2_C3_bounds'][0][label][k] for p in frames],proof,'compact_local_'+name+'_Z'+str(k),k) if active else iv.mpf(0) for k in range(3)])
        if x=='1':corrections={n:self._jet(inlet['finite_N_local_errors'][n],True) for n in ('E','V')}
        p0p=self._jet(inlet['independent_P0_over_Pstar_squared'],False).truncate(2);p0t=self._jet(inlet['independent_P0_over_Pstar_squared'],True).truncate(2)
        point['pressure']=point['p']+p0p;truth['pressure']=truth['p']+p0t
        modified={n:truth[n]+history[n] for n in RATES};modified.update(E=truth['E']+corrections['E'],V=truth['V']+corrections['V'],pressure=truth['pressure']+history['p'])
        Pp,Pt=c.exp(self.patch.leading.defect.logPp),iv.exp(self.patch.leading.defect.logP)
        z=J.variable(iv,iv.mpf(Z),2);dm=history['m'];dmz=J(iv,[dm[1],2*dm[2]]);delta=iv.mpf(self.patch.leading.reference.leading.core.delta_iv)
        dQ=((z*corrections['V']*2-z*dm*(1-delta)-(1-z*z)*dmz)/(1-z*z*delta)*Pt).truncate(1)
        simpoint=dict(Uz=point['V']*Pp,Utheta_over_Pstar=point['E'],Mz_over_R=point['m']*Pp,Ur_over_sqrt_R_over2=Qp)
        simtruth=dict(Uz=modified['V']*Pt,Utheta_over_Pstar=modified['E'],Mz_over_R=modified['m']*Pt,Ur_over_sqrt_R_over2=Qt+dQ)
        ref=self.patch.leading.reference
        packet=dict(Z=Z,native_compact_x=x,chart='actual_original_compact_Rc_2Rc',source_family=self.family,selected_integer=inlet['selected_integer'],
            generic_leading_fields={n:self._rows(point[n],truth[n]) for n in point},generic_modified_fields={n:self._rows(point[n],modified[n]) for n in modified},
            physical_similarity={n:self._rows(simpoint[n],simtruth[n]) for n in simpoint},leading_small_sectors={n:self._rows(sectorsp[n],sectorst[n]) for n in sectorsp},
            normalized_mu_profiles={n:ref._summary(scalar['point'][n],scalar['directed'][n]) for n in scalar['point']},original_mu=inlet['original_mu'],
            finite_N_history_errors={n:self._rows(zero,v) for n,v in history.items()},finite_N_local_errors={n:self._rows(zero,v) for n,v in corrections.items()},
            incoming_actual_finite_N_history=inlet['finite_N_history_errors'],forward_history_enclosures={n:self._rows(zero,v) for n,v in forward.items()},
            backward_terminal_identity_enclosures={n:self._rows(zero,v) for n,v in backward.items()},physical_supported_masses=encode(masses),
            original_control_source=self.control_source,relative_terminal_source_identities=self.certificates,
            normalized_target_magnitude_ranges=[dict(exact_Z_cell=p['exact_Z_cell'],targets=p['actual_N_scaled_five_target_C3_ranges'],joint=p['actual_joint_angular_numerator_C3_range']) for p in self.owner.ranges if p['exact_Z_cell'] in [f['exact_Z_cell'] for f in frames]],
            independent_P0_over_Pstar_squared=inlet['independent_P0_over_Pstar_squared'],incoming_n1_history_not_extended=inlet['positive_order_n1_own_history'],
            coordinate=dict(global_y_log_terms=[ref._summary(10*ref.logCp,10*ref.logC),ref._summary(10*self.patch.leading.defect.logPp,10*self.patch.leading.defect.logP),
                ref._summary(c.exp(40),iv.exp(40)),ref._summary(c.mpf(14),iv.mpf(14)),ref._summary(sp,st)]),
            compact_source_cover=dict(exact_Z_cells=[p['exact_Z_cell'] for p in frames],active_bumps=local,log_bump_supports=encode(self.supports),same_actual_original_N=True),
            small_error_proofs=proof,contract=dict(unchanged_original_unique_control_branch=True,zero_is_only_point_proxy_for_unknown_signed_correction=True,
                incoming_and_density_not_divided_by_N_again=True,relaxed_inlet_not_multiplied_by_N=True,joint_target_not_split=True,
                true_forward_and_backward_own_rate_integrals=True,backward_zero_from_original_signed_function_identity=True,
                only_relative_five_histories_zero_at_terminal=True,P0_added_once=True,original_mu_retained=True),
            limits=dict(physical_Cartesian_assembly=False,resolved_signed_finite_N_centers=False,original_control_signed_values_resolved=False,
                compact_n1_numeric_fields=False,total_radial_or_fast_phase_derivatives=False,reduced_m13_m15_integrals_continued=False,
                relative_pressure_cancellation_accuracy=False,admissible_stress_cone_margin=False,beyond_2Rc=False,
                global_N_admission=False,global_field=False,genuine_temporal_n_ge_2=False,full_corrected_NS=False))
        self._packets[key]=packet;return packet
