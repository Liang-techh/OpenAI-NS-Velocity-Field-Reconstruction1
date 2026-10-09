"""Actual Rm source-bound phase primitives and endpoint-retaining C0 IBP.

The original leading vector has exact zero phase mean. Its bounded phase
primitive, true slow-y derivative and nonlinear remainder give N^-2
weighted covers. Ordinary Z keeps the accepted independent direct bound;
no mixed yZ or sharper Z averaging is invented.
"""
from fractions import Fraction
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rm_first_spatial_jets as previous
import lei_ren_part1_paper_compliant_current_original_Rm_terminal_density_integrals as terminal
import lei_ren_part1_paper_compliant_current_native_signed_averaging as identities

fields,base,ep=previous.fields,previous.base,previous.ep
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_Rm_weighted_averaging.json.gz'
RECEIPT=PREFIX+'current_original_Rm_weighted_averaging_check.json'
GATE='original_actual_Rm_source_phase_primitives_and_endpoint_retaining_C0_weighted_averaging_installed'
RATES=terminal.RATES
PARTITION=((1,1),(49,40),(51,40),(59,40),(61,40),(69,40),(71,40),(2,1),(9,4),(5,2),'Rh')


def magnitude(f,row):
    previous.phase.previous.same_source(f,[row])
    if row.zero:return f.scalar(0)
    c=f.c;scale=type(row.scale)(f.logs,row.scale.powers,c.mpf(ep(row.scale.offset)[1]))
    return type(row)(scale,max(abs(v) for v in ep(row.coefficient)),f.ledger)


def symmetric(f,bound):return magnitude(f,bound)*f.c.mpf([-1,1])


def tighter(f,left,right):
    """Choose a valid bound, never an endpoint as a source function value."""
    left,right=magnitude(f,left),magnitude(f,right)
    if left.zero:return left,'left_exact_zero'
    if right.zero:return right,'right_exact_zero'
    ratio=(left.scale-right.scale).evaluate()+f.c.ln(left.coefficient)-f.c.ln(right.coefficient)
    if ep(ratio)[0]>0:return right,'right_strictly_smaller'
    return left,'left_valid_bound'


def read_saved_source(f,row):
    scale=row['formal_positive_scale'];read=lambda value:fields.previous.read_interval(f.c,value)
    model=f.scalar(0)
    return type(model)(type(model.scale)(f.logs,tuple(scale['source_exponents'])+(scale['radius_power'],),
        read(scale['additional_log_interval'])),read(row['coefficient_interval']),f.ledger)


def leading_pairs(E,Ey,V,Vy,A,Ay,B,By):
    pair=lambda v,d:(v,d)
    def add(*pairs):return tuple(sum((v[n] for v in pairs),E.scalar(0)) for n in (0,1))
    def mul(l,r):return l[0]*r[0],l[1]*r[0]+l[0]*r[1]
    scale=lambda p,v:tuple(row*v for row in p)
    ee,vv,aa,bb=(pair(v,d) for v,d in ((E,Ey),(V,Vy),(A,Ay),(B,By)))
    EA=mul(ee,aa);EEA=mul(ee,EA)
    return dict(m=bb,h=EA,k=add(mul(vv,EA),mul(ee,bb)),
        e=add(scale(mul(vv,bb),2),scale(EEA,-1)),p=EEA)


def phase_distance(c,boxes):
    """sup min(phi,1-phi), for every actual endpoint phase enclosure."""
    if not boxes:raise ValueError('Actual endpoint phase boxes required')
    bounds=[]
    for box in boxes:
        lo,hi=ep(c.mpf(box))
        if lo<0 or hi>1:raise ValueError('Original closed fractional phase required')
        bounds.append(c.mpf('.5') if lo<=mp.mpf('.5')<=hi else c.mpf(max(min(lo,1-lo),min(hi,1-hi))))
    return c.mpf(max(ep(v)[1] for v in bounds))


def phase_primitive_enclosure(f,source,phi):
    """Enclose the original integral-defined G and its true slow-y derivative.

    Both endpoints vanish by the original phase reflection identity. The
    two complementary integrals provide the distance-to-endpoint bound;
    these intervals are not selected values of G or derivatives of caps.
    """
    if not source.get('original_leading_zero_phase_mean') or not source.get('phase_primitive_endpoints_exactly_zero'):
        raise ValueError('Original same-source zero-mean phase primitive required')
    pairs=source['original_leading_signed_density_pairs']
    if set(pairs)!=set(RATES):raise ValueError('All five genuine original leading pairs required')
    distance=phase_distance(f.c,[phi]);values={};slow={}
    for key,pair in pairs.items():
        previous.phase.previous.same_source(f,list(pair))
        values[key]=symmetric(f,magnitude(f,pair[0])*distance)
        slow[key]=symmetric(f,magnitude(f,pair[1])*distance)
    return dict(original_phase_primitive_C0=values,original_phase_primitive_slow_y=slow,
        original_closed_phi=f.c.mpf(phi),distance_to_original_period_endpoint_cap=distance,
        same_original_integral_definition=source['source_phase_primitive_definition'],
        genuine_y_derivative_definition=source['source_phase_primitive_y_definition'],
        source_primitive_function_enclosures_not_selected_values=True,
        genuine_fixed_phi_y_not_spatial_N_chain=True,**dict.fromkeys(fields.previous.OPEN,False))


def weighted_cap(f,leading,slow,remainder,width,suffix,rate,N,left_distance,right_distance):
    """Exact own-rate IBP bound; all true endpoint terms are retained."""
    c=f.c;N=previous.phase.candidate_N(N);lam=c.mpf(rate.numerator)/rate.denominator
    width,suffix=c.mpf(width),c.mpf(suffix)
    if ep(width)[0]<=0 or ep(suffix)[0]<0:raise ValueError('Ordered positive cell and downstream suffix required')
    for distance in (left_distance,right_distance):
        if ep(c.mpf(distance))[0]<0 or ep(c.mpf(distance))[1]>mp.mpf('.5'):
            raise ValueError('Original phase primitive distance bound required')
    mass=terminal.local.positive_kernel_mass(c,width,rate)
    decay=c.mpf(1) if not rate else c.exp(-lam*width)
    downstream=c.mpf(1) if not rate else c.exp(-lam*suffix)
    leading,slow,remainder=(magnitude(f,v) for v in (leading,slow,remainder))
    G=leading*c.mpf('.5');Gy=slow*c.mpf('.5')
    endpoints=leading*(c.mpf(right_distance)+decay*c.mpf(left_distance))
    slow_integral=(Gy+G*lam)*mass;rest=remainder*mass
    cap=(endpoints+slow_integral+rest)*(downstream/c.mpf(N*N))
    return cap,dict(positive_own_rate_mass=mass,incoming_cell_decay=decay,downstream_suffix_decay=downstream,
        endpoint_primitive_cap=endpoints,slow_and_weight_derivative_integral_cap=slow_integral,
        nonlinear_second_remainder_integral_cap=rest,phase_primitive_cap=G,slow_phase_primitive_y_cap=Gy,
        actual_endpoint_distance_caps=dict(left=left_distance,right=right_distance),
        exact_bound_recipe='suffix*(endpoint+mass*(G_y_slow+lambda*G+remainder))/N^2',
        original_log_radius_measure_once=True,all_cell_endpoints_retained=True,
        nonlinear_mean_not_zeroed=True,internal_endpoint_cancellation_not_claimed=True)


class OriginalRmWeightedAveraging:
    mode='actual_Rm_original_phase_primitive_and_own_rate_C0_averaging_with_direct_Z_preserved'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalRmFirstSpatialJets(dps);self.c=self.upstream.c;self.family=self.upstream.family
        self.hashes=dict(self.upstream.hashes);self.theorem=identities.exact_theorem();self.cache={}
        for name in (Path(__file__).name,Path(identities.__file__).name):fields.previous.bind(self.hashes,name,sha(name))
        self.baseline=json.loads(gzip.decompress((HERE/previous.previous.NAME).read_bytes()))
        if self.baseline['source_family']!=self.family or not self.baseline[previous.previous.GATE]:
            raise ValueError('Same accepted complete local source integral baseline required')
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted actual Rm weighted averaging receipt required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def source_cell(self,label,left,right,N):
        key=(label,left,right,N)
        if key in self.cache:return self.cache[key]
        p=self.upstream.phase;mixed=p.upstream.upstream.upstream.owner(label);op=mixed.op;f=op.flow;c=op.c
        patch=mixed.cell(left,right);generic=p.upstream.upstream.cell(label,left,right)
        quotients=previous.phase.previous.recover_quotients(op,generic,p.upstream.eta_log,p.upstream.dstar_log,patch=patch)
        source,qr,record=previous.source_first_frame(op,patch,generic,quotients)
        whole=previous.all_u_first_bounds(f,source,qr,p.upstream.dstar_log,c.mpf([0,1]))
        got=previous.phase.first.conditioned_first_jets(source,qr,p.upstream.dstar_log,c.mpf([0,1]))
        values=whole['values'] if got['values'] is None else got['values']
        caps={name:tighter(f,magnitude(f,value),magnitude(f,whole['values'][name]))[0] for name,value in values.items()}
        A,B,Ay,By=(symmetric(f,caps[name]) for name in ('A','B_over_Pstar','A_y','B_y_over_Pstar'))
        E=source['roots']['E'][(0,0)];Ey=source['roots']['E'][(1,0)];V=generic['common_velocity_V_axial5'][0]
        Vy=previous.lift_generic_y(op,patch)['common_velocity_V_axial5'][0].dy
        leading=leading_pairs(E,Ey,V,Vy,A,Ay,B,By)
        acap=caps['A'];argument=previous.phase.first.phase.bounded_value(acap*(c.mpf(1)/N))
        if ep(argument)[0]<0 or ep(argument)[1]>1:raise ValueError('Original bounded A/N source required')
        exponential=c.exp(argument);e0,v0=magnitude(f,E),magnitude(f,V);b0=caps['B_over_Pstar']
        Q=e0*previous.phase.first.current.square(acap)*exponential*c.mpf('.5')
        F=e0*acap*exponential
        remainder=dict(m=f.scalar(0),h=Q,k=v0*Q+F*b0,
            e=e0*Q+previous.phase.first.current.square(b0)+previous.phase.first.current.square(F)*c.mpf('.5'),
            p=e0*Q+previous.phase.first.current.square(F)*c.mpf('.5'))
        density=previous.phase.densities.density_Z_kernels(E,source['roots']['E'][(0,1)],V,
            generic['common_velocity_V_axial5'][1],dict(values,A=A,B_over_Pstar=B),N)
        if label not in p.mappers:p.mappers[label]=previous.phase.RmRadiusPhase(op,self.family,p.parameter_family)
        mapper=p.mappers[label];endpoints={which:mapper.point(x,N) for which,x in (('left',left),('right',right))}
        result=dict(source_family=self.family,source_geometry=generic['source_geometry'],candidate_N=N,
            exact_common_P0_axial5=op.P0,exact_source_Rm_factor=op.Rm_factor,genuine_slow_y_source=record,
            original_whole_phase_primitive_caps=caps,original_leading_signed_density_pairs=leading,
            original_nonlinear_second_remainder_caps=remainder,original_direct_density_C0=density['kernels'],
            actual_endpoint_phases=endpoints,
            source_phase_primitive_definition='G_j(y,Z,phi)=integral_0^phi f1_j(y,Z,s) ds of the same original inverse functions',
            source_phase_primitive_y_definition='G_j_y_slow=integral_0^phi f1_j_y_slow(y,Z,s) ds',
            phase_primitive_endpoints_exactly_zero=True,primitive_caps_are_bounds_not_function_values=True,
            original_leading_zero_phase_mean=True,original_full_density_mean_not_zeroed=True,
            original_slow_y_not_total_spatial_N_chain=True,
            original_all_u_y_source_theorem=whole['record'],original_conditioned_whole_phase_record=got['record'],
            **dict.fromkeys(fields.previous.OPEN,False))
        result['actual_endpoint_phase_primitive_enclosures']={key:[phase_primitive_enclosure(f,result,box)
            for box in endpoint['phase_boxes']] for key,endpoint in endpoints.items()}
        self.cache[key]=result;return result

    def contribution(self,label,N=257):
        N=previous.phase.candidate_N(N);op=self.upstream.phase.upstream.upstream.upstream.owner(label).op;f=op.flow;c=op.c
        cells=[];total={key:f.scalar(0) for key in RATES};direct={key:f.scalar(0) for key in RATES}
        with mp.workdps(c.dps+40):
            xs=[c.exp(1) if x=='Rh' else c.mpf(x[0])/x[1] for x in PARTITION]
            for left,right,xl,xr in zip(PARTITION,PARTITION[1:],xs,xs[1:]):
                source=self.source_cell(label,left,right,N);width=c.ln(xr/xl);suffix=c.mpf(0) if right=='Rh' else 1-c.ln(xr)
                distances={k:phase_distance(c,v['phase_boxes']) for k,v in source['actual_endpoint_phases'].items()}
                rows={};details={}
                for key,rate in RATES.items():
                    leading=source['original_leading_signed_density_pairs'][key]
                    cap,proof=weighted_cap(f,*leading,source['original_nonlinear_second_remainder_caps'][key],
                        width,suffix,rate,N,distances['left'],distances['right'])
                    mass=proof['positive_own_rate_mass'];decay=proof['downstream_suffix_decay']
                    simple=magnitude(f,source['original_direct_density_C0'][key])*mass*decay
                    used,route=tighter(f,cap,simple)
                    rows[key]=used;total[key]+=used;direct[key]+=simple
                    details[key]=dict(original_IBP=proof,weighted_Nminus2_cap=cap,original_direct_cap=simple,
                        used_valid_cap=used,used_route=route)
                cells.append(dict(source=source,source_log_width=width,source_suffix_to_Rh=suffix,
                    actual_cell_to_Rh_C0_caps=rows,original_own_rate_averaging=details))
            baseline=self.baseline['frames'][label]['actual_whole_patch_local_defect_integral_C0_Z']
            output={};comparisons={}
            for key in RATES:
                old=read_saved_source(f,baseline[key][0]);oldcap=magnitude(f,old)
                used,route=tighter(f,total[key],oldcap)
                output[key]=[symmetric(f,used),read_saved_source(f,baseline[key][1])]
                reduction=None if oldcap.zero or used.zero else (oldcap.scale-used.scale).evaluate()+c.ln(oldcap.coefficient)-c.ln(used.coefficient)
                comparisons[key]=dict(accepted_direct_whole_C0_cap=oldcap,new_cellwise_valid_C0_cap=total[key],
                    used_C0_cap=used,used_route=route,log_absolute_upper_reduction=reduction,
                    strict_absolute_upper_reduction=(not oldcap.zero and used.zero) or reduction is not None and ep(reduction)[0]>0,
                    ordinary_Z_accepted_direct_bound_unchanged=True)
        return dict(source_family=self.family,source_frame=label,candidate_N=N,exact_partition=PARTITION,
            actual_source_weighted_cells=cells,actual_averaged_whole_patch_local_integral_C0_Z=output,
            actual_direct_new_partition_C0_caps=direct,actual_C0_comparisons=comparisons,
            strict_C0_reductions=sum(row['strict_absolute_upper_reduction'] for row in comparisons.values()),
            endpoint_retaining_C0_Nminus2_averaging_installed=True,ordinary_Z_Nminus2_averaging_installed=False,
            genuine_mixed_yZ_derivative_still_required_for_Z_averaging=True,
            real_finite_N_Rm_incoming_correction_is_still_unsupplied=True,
            exact_common_P0_axial5=op.P0,source_integral_caps_not_five_moment_closure=True,
            **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalRmWeightedAveraging(require_checked=False);frames={}
    for label in ('0','.5'):
        frames[label]=owner.contribution(label)
        print('Actual Rm weighted C0 averaging',label,'strict reductions',frames[label]['strict_C0_reductions'],flush=True)
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_exact_leading_remainder_and_phase_primitive_theorem=owner.theorem,frames=fields.serialized(frames),
        input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()
