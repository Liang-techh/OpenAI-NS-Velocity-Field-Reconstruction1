"""Issued full O2 source predicates with separate original q and |u| logs.

Each frame encloses the same original functions restricted by its named
predicate. The outer physical rectangle is not asserted entirely regular
or signed. Together abs(u)<=1/4 and u>=3/16/u<=-3/16 cover every source.
"""
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
from types import FunctionType,MappingProxyType
import time
import mpmath as mp
import sympy as s
import lei_ren_part1_paper_compliant_current_original_O2_regular_mixed_phase as regular
import lei_ren_part1_paper_compliant_current_original_O2_axial_carrier as carrier
import lei_ren_part1_paper_compliant_current_original_O2_positive_q_curvature as curve

source=regular.source;base,prior,ep=source.base,source.prior,source.ep
MixedJet=source.MixedJet;C0,Y,Z,YZ=source.C0,source.Y,source.Z,source.YZ
HERE,PREFIX,sha=source.HERE,source.PREFIX,source.sha
NAME=PREFIX+'current_original_O2_predicate_source_frames.json.gz'
RECEIPT=PREFIX+'current_original_O2_predicate_source_frames_check.json'
GATE='original_O2_full_regular_signed_issued_predicate_source_frames_bound'


class O2QUAtlas(source.O2MixedAtlas):
    """Last slot is now defined logabsu; no radius factor remains there."""
    def __init__(self,frame,*,lower,upper,logq,logu):
        super().__init__(frame,lower=lower,upper=upper,logq=logq)
        self.bases=self.bases[:4]+(self.ctx.mpf(ep(logu)),)
        if ep(self.bases[4])[0]<-2:raise ValueError('Bounded lower original logabsu or exact zero required')

    def add(self,left,right):
        left.scale.pair(right.scale)
        if left.scale.bases is not self.bases or left.ledger is not self.ledger or right.ledger is not self.ledger:
            raise ValueError('One original q/u basis and ledger required')
        if left.zero:return right
        if right.zero:return left
        lp,rp=left.scale.powers,right.scale.powers
        powers=(max(lp[0],rp[0]),max(lp[1],rp[1]),min(lp[2],rp[2]),min(lp[3],rp[3]),max(lp[4],rp[4]))
        anchor=prior.FormalScale(self.bases,powers)
        coeffs=[value.coefficient*value.bounded_exp((value.scale-anchor).evaluate()) for value in (left,right)]
        self.ledger['collected_O2_P_C_L_q_u_additions']=self.ledger.get('collected_O2_P_C_L_q_u_additions',0)+1
        return prior.ScaledEnclosure(anchor,coeffs[0]+coeffs[1],self.ledger)

    def record(self):
        return dict(source_family=self.family,exact_outer_Z_bounds=[str(v) for v in self.bounds],
            source_basis_order=['logPstar','selected_logCstar','logL','original_logq','defined_original_logabsu_or_zero'],
            defining_basis=self.bases,original_Q=self.Q,original_L=self.L,
            source_radius_collected_into_P_C_before_native_logs=True,
            legacy_scale_record_radius_power_field_is_logabsu_power_on_this_atlas=True,
            actual_radius_power_is_exact_zero_and_not_in_fifth_slot=True,
            last_slot_not_compatible_with_previous_zero_or_radius_slot=True,
            original_q_and_u_are_function_logs_not_selected_values=True)


def rebase_original_root(a,old,*,raw_frame,source_owner):
    source_owner.describe(raw_frame)
    raw_atlas=raw_frame.roots['q'].atlas
    if old.scale.bases is not raw_atlas.bases or old.ledger is not raw_atlas.ledger or old.scale.powers[4]!=0:
        raise ValueError('Issued original zero-last-slot O2 source value required')
    if raw_frame.family!=a.family:raise ValueError('Original source families differ')
    # Only the q log may be predicate-restricted; P/C/L stay the same source.
    for index in range(3):
        if ep(a.bases[index])!=ep(raw_atlas.bases[index]):raise ValueError('Original P/C/L source basis changed')
    if not ep(raw_atlas.bases[3])[0]<=ep(a.bases[3])[0]<=ep(a.bases[3])[1]<=ep(raw_atlas.bases[3])[1]:
        raise ValueError('Predicate q range must be a subset of the original positive source cover')
    return prior.ScaledEnclosure(prior.FormalScale(a.bases,old.scale.powers,a.copy_interval(old.scale.offset)),
        a.copy_interval(old.coefficient),a.ledger)


_old_init=source.axial.PREDICATE_PHASE_INIT
_closed_init=FunctionType(_old_init.__code__,dict(_old_init.__globals__,bounded_value=source.positive.bounded),
    _old_init.__name__,_old_init.__defaults__,_old_init.__closure__)
_closed_init.__kwdefaults__=_old_init.__kwdefaults__


class O2PredicatePhase(source.positive.PositiveLogQPhase):
    def __init__(self,query,dstar_log):
        if query.get('regular_predicate',False)==query.get('signed_predicate',False):
            raise ValueError('Exactly one proved original source geometry predicate required')
        q=query['q'];u=query['original_u_source'];c=q.ctx
        q.scale.pair(u.scale)
        if q.ledger is not u.ledger or q.zero or ep(q.coefficient)[0]<=0:
            raise ValueError('Same-source strictly positive original q required')
        if query['regular_predicate']:
            if max(abs(v) for v in ep(source.positive.bounded(u)))>ep(c.mpf(1)/4)[1]:
                raise ValueError('Original closed regular-u predicate not proved')
        else:
            lo,hi=ep(u.coefficient)
            if not lo*hi>0:raise ValueError('Strict original u sign required')
            lower=u.scale.evaluate()+c.ln(c.mpf(min(abs(lo),abs(hi))))
            if ep(lower)[0]<ep(c.ln(c.mpf(3)/16))[0]:raise ValueError('Original abs(u)>=3/16 predicate required')
            p2=query['roots']['p2'][C0]
            if not ep(p2.coefficient)[0]*ep(p2.coefficient)[1]>0 or (ep(p2.coefficient)[0]>0)!=(lo>0):
                raise ValueError('Original p2/u strict signs must agree')
        _closed_init(self,query,dstar_log)
        if self.geometry=='signed_Mobius':
            bound=3/c.sqrt(265);ar=source.base.conditioned.clipped(c,abs(self.r),ep(bound)[0],1)
            self.r=ar*self.sign


@dataclass(frozen=True)
class O2PredicateFrame:
    owner:object
    family:object
    count:int
    index:int
    branch:str
    roots:object
    u:object
    gamma_g:object
    gamma_q:object
    kernel:object
    record:dict


class OriginalO2PredicateSources:
    def __init__(self):
        self.source=source.OriginalO2MixedSources();self.owner=self.source.owner;self.ctx=self.source.ctx;self.family=self.source.family
        self.hashes=dict(self.source.hashes)
        for module in (source,carrier,curve,regular):
            receipt=json.loads((HERE/module.RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(module.GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted same-family source, carrier, curvature and regular phase required')
            for name,digest in {**receipt['input_hashes'],module.RECEIPT:sha(module.RECEIPT)}.items():
                if sha(name)!=digest:raise ValueError('Original O2 predicate dependency changed: '+name)
                if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original O2 predicate closures disagree')
                self.hashes[name]=digest
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        cert=json.loads(gzip.decompress((HERE/carrier.NAME).read_bytes()))['actual_O2_whole_axial_carrier_certificate']
        self.certificate=cert;self.g=source.ordered.interval(self.ctx,cert['g_normalized'])
        self.gamma=source.ordered.interval(self.ctx,cert['g_y_over_g'])
        if cert['sign']!=-1 or ep(self.g)[1]>=0:raise ValueError('Strict full original O2 carrier sign required')
        self.frames={};self.cache={}

    def frame(self,count,index,*,branch):
        if branch not in ('regular','positive','negative'):raise ValueError('Named actual original u predicate required')
        key=(count,index,branch)
        if key in self.cache:return self.cache[key]
        bounds=dict(regular=(-1,1),positive=(-1,0),negative=(0,1))[branch]
        raw=self.source.source_frame(count,index,Z_lower=bounds[0],Z_upper=bounds[1]);rawa=raw.roots['q'].atlas;c=self.ctx
        with mp.workdps(c.dps+40):
            logq=c.mpf(ep(rawa.bases[3]));original_logq=c.mpf(ep(logq));logd=rawa.copy_interval(self.owner.scales.logs['d_star'])
            logLambda=11*rawa.bases[0]+10*rawa.bases[1];g_upper=c.mpf(-ep(self.g)[0])
            logg=c.ln(g_upper);guard=c.ln(c.mpf(3)/16)
            q_predicate_lower=guard+logd-logLambda-logg
            if branch!='regular':
                required=max(ep(logq)[0],ep(q_predicate_lower)[0])
                if required>ep(logq)[1]:raise ValueError('Original signed source predicate has empty q enclosure')
                logq=c.mpf((required,ep(logq)[1]))
                upper=logLambda+logg+c.mpf(ep(logq)[1])-logd
                logu=c.mpf((ep(guard)[0],ep(upper)[1]))
            else:logu=c.mpf(0)
            a=O2QUAtlas(self.owner.inputs.frame,lower=bounds[0],upper=bounds[1],logq=logq,logu=logu)
            roots={name:MixedJet(a,{order:rebase_original_root(a,value,raw_frame=raw,source_owner=self.source)
                for order,value in row.rows.items()}) for name,row in raw.roots.items()}
            q=roots['q'][C0];d=prior.ScaledEnclosure(prior.FormalScale(a.bases,offset=a.copy_interval(logd)),1,a.ledger)
            if branch=='regular':u0=a.scalar(c.mpf((-c.mpf(1)/4,c.mpf(1)/4)))
            else:u0=prior.ScaledEnclosure(prior.FormalScale(a.bases,(0,0,0,0,1)),1 if branch=='positive' else -1,a.ledger)
            # These are source identities on a named restricted domain, never
            # derivatives of its bounding interval or selected source values.
            gamma_g=a.scalar(a.copy_interval(self.gamma))
            q_lower=q.scale.evaluate()+c.ln(c.mpf(ep(q.coefficient)[0]))
            gamma_q=roots['q'][Y].positive_divide(q,q_lower)
            p2=(d*u0).positive_divide(q,q_lower)
            old=roots['p2']
            roots['p2']=MixedJet(a,{C0:p2,Y:p2*gamma_g,Z:old[Z],YZ:old[YZ]})
            u=MixedJet(a,{C0:u0,Y:u0*a.add(gamma_g,gamma_q),
                Z:(roots['p2'][Z]*q).positive_divide(d,logd),
                YZ:a.add(roots['p2'][YZ]*q,roots['p2'][Z]*roots['q'][Y]).positive_divide(d,logd)})
            query=dict(q=q,roots={name:{C0:row[C0],Z:row[Z]} for name,row in roots.items()},
                original_u_source=u0,regular_predicate=branch=='regular',signed_predicate=branch!='regular')
            kernel=O2PredicatePhase(query,logd)
            expected='small_r_series' if branch=='regular' else 'signed_Mobius'
            if kernel.geometry!=expected:raise ArithmeticError('Original predicate backend differs')
            kernel.nu=roots['nu'][C0]
            kernel.n0=kernel.normalized(kernel.scalar(1),1,True)
            kernel.nt=kernel.normalized(kernel.scalar(0),1,True)
            kernel.nq=kernel.normalized(base.current.square(q),c.mpf('.5'),True)
            kernel.ntq=kernel.normalized(kernel.scalar(0),1/(2*c.sqrt(2)),False)
            record=dict(source_family=self.family,source_level=count,source_index=index,exact_y_cell=[str(raw.left),str(raw.right)],
                exact_outer_Z_bounds=[str(v) for v in bounds],branch=branch,
                actual_original_source_predicate='abs(u)<=1/4' if branch=='regular' else ('u>=3/16' if branch=='positive' else 'u<=-3/16'),
                original_u_definition='u=p2*q/dstar',atlas=a.record(),actual_root_jets={name:row.record() for name,row in roots.items()},
                actual_u_mixed=u.record(),same_source_gamma_g=gamma_g.record(),same_source_gamma_q=gamma_q.record(),
                exact_source_identities=['p2=dstar*u/q','p2_y=p2*(g_y/g)','u_y/u=g_y/g+q_y/q where u!=0',
                    'u_Z=p2_Z*q/dstar','u_yZ=(p2_yZ*q+p2_Z*q_y)/dstar'],
                p2_C0_y_reenclosed_from_actual_source_identities=True,original_p2_Z_yZ_template_and_pressure_error_rows_retained=True,
                original_positive_logq_before_predicate=original_logq,predicate_logq=logq,
                signed_predicate_q_lower_implication=q_predicate_lower,
                q_lower_implication_uses_actual_p2_upper='|p2|<=Lambda0*max|g/Lambda0|, |Z|<=1; |u|>=3/16 implies q>=(3/16)*dstar/(Lambda0*max|g/Lambda0|)',
                original_pressure_family_and_incoming_histories_unchanged=True,
                conditional_domain_not_claimed_entire_outer_rectangle=True,
                predicate_not_selected_field_or_derivative_of_cap=True,
                q_stays_original_positive_function_with_actual_y_derivative=True,
                source_function_derivatives_not_compatible_selected_interval_fields=True,
                same_bases_context_and_ledger_for_all_rebuilt_roots=True,
                exact_radius_collected_and_last_slot_semantics_changed=True,
                kernel_geometry=kernel.geometry_record(),actual_mixed_primitives_installed=False)
            result=O2PredicateFrame(self,self.family,count,index,branch,MappingProxyType(roots),u,gamma_g,gamma_q,kernel,record)
            self.frames[id(result)]=result;self.cache[key]=result;return result

    def describe(self,frame):
        if type(frame) is not O2PredicateFrame or frame.owner is not self or self.frames.get(id(frame)) is not frame:
            raise ValueError('Issued same-owner original O2 predicate frame required')
        return frame.record


def run():
    began=time.monotonic();owner=OriginalO2PredicateSources();count=min(owner.source.parent.parent.levels)
    rows=[{branch:owner.describe(owner.frame(count,index,branch=branch)) for branch in ('regular','positive','negative')} for index in range(count)]
    report=dict(**{GATE:True},source_family=owner.family,whole_original_O2_predicate_source_cells=rows,
        exact_outer_original_domain=dict(y=['0','1'],Z=['-1','1']),
        coverage='abs(u)<=1/4 OR u>=3/16 OR u<=-3/16; overlap3/16<=abs(u)<=1/4',
        all_real_original_u_including_zero_covered_by_predicates=True,
        actual_positive_q_and_q_y_nu_y_retained=True,actual_mixed_primitives_installed=False,
        actual_changed_five_integrals_installed=False,all_17_chart_or_24_cell_oracle_installed=False,
        actual_five_controls_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(source.ordered.base.point.source.inertial.profiles.loop.OPEN,False),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Issued actual original O2 full regular/positive/negative source predicates with separate logq/logabsu atlas, actual mixed roots, carrier-correlated y rows, original transverse rows and positive q. Predicates cover full original y/Z domain; source interface only, mixed primitives/integrals and full reconstruction remain open.')
    (HERE/NAME).write_bytes(gzip.compress(json.dumps(base.encoded(report),indent=2).encode()+b'\n',compresslevel=9,mtime=0))
    print('Full original O2 issued source frames:',count,'radial cells by3 predicates',flush=True);return report


if __name__=='__main__':run()
