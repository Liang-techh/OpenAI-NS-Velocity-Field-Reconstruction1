"""Actual original Rm axial phase primitives and endpoint-retaining IBP.

Genuine fixed-phi Z/yZ rows differentiate the original zero-mean leading
functions. Every nonlinear second remainder is differentiated. Existing
C0 bounds are preserved and an accepted direct Z cover is retained when
the conservative mixed bound does not improve it. No inlet is supplied.
"""
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_Rm_mixed_yZ_primitives as previous

averaging=previous.previous
first,phase=previous.first,previous.phase
fields,base,ep=previous.fields,previous.base,previous.ep
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_Rm_weighted_Z_averaging.json.gz'
RECEIPT=PREFIX+'current_original_Rm_weighted_Z_averaging_check.json'
GATE='original_actual_Rm_genuine_G_Z_G_yZ_and_endpoint_retaining_weighted_Z_averaging_installed'
RATES,PARTITION=averaging.RATES,averaging.PARTITION
C0,Y,Z,YZ=previous.C0,previous.Y,previous.Z,previous.YZ


def require_phase_Z_zero(source):
    if source.get('original_radius_phase_Z_exact_zero') is not True:
        raise ValueError('Original radius phase must have exactly zero Z derivative')
    if not source.get('original_leading_zero_phase_mean') or not source.get('phase_primitive_endpoints_exactly_zero'):
        raise ValueError('Same original zero-mean phase primitive required')


def phase_primitive_Z_enclosure(f,source,phi):
    """Original G_Z/G_yZ integral enclosures, never derivatives of caps."""
    require_phase_Z_zero(source)
    rows=source['original_five_signed_leading_density_C0_y_Z_yZ']
    if set(rows)!=set(RATES) or any(set(row)!=set(previous.ORDERS) for row in rows.values()):
        raise ValueError('All five genuine original Z/yZ leading rows required')
    distance=averaging.phase_distance(f.c,[phi]);values={};mixed={}
    for key,row in rows.items():
        phase.previous.same_source(f,list(row.values()))
        values[key]=averaging.symmetric(f,averaging.magnitude(f,row[Z])*distance)
        mixed[key]=averaging.symmetric(f,averaging.magnitude(f,row[YZ])*distance)
    return dict(original_phase_primitive_Z=values,original_phase_primitive_slow_yZ=mixed,
        original_closed_phi=f.c.mpf(phi),distance_to_original_period_endpoint_cap=distance,
        original_G_Z_definition='integral_0^phi f1_Z(y,Z,s) ds',
        original_G_yZ_definition='integral_0^phi f1_yZ_fixed_phi(y,Z,s) ds',
        genuine_original_mixed_functions_not_derivatives_of_selected_bounds=True,
        phase_Z_exact_zero_required=True,**dict.fromkeys(fields.previous.OPEN,False))


def remainder_Z_caps(f,E,V,values,N):
    """Differentiate all exact second remainders; normalize B only once."""
    N=phase.candidate_N(N);c=f.c;absolute=lambda row:averaging.magnitude(f,row)
    e,ez,v,vz=(absolute(rows[key]) for rows,key in ((E,C0),(E,Z),(V,C0),(V,Z)))
    a,az,b,bz=(absolute(values[key]) for key in ('A','A_Z','B_over_Pstar','B_Z_over_Pstar'))
    square=phase.first.current.square
    argument=phase.first.phase.bounded_value(a*(c.mpf(1)/N))
    if ep(argument)[0]<0 or ep(argument)[1]>1:raise ValueError('Original bounded |A|/N required')
    exponential=c.exp(argument)
    Q=e*square(a)*exponential*c.mpf('.5')
    QZ=(ez*square(a)+e*a*az*2)*exponential*c.mpf('.5')+e*square(a)*az*exponential*(c.mpf(1)/(6*N))
    F=e*a*exponential
    # Exact F=N*E*expm1(A/N), hence F_Z=N*EZ*expm1(A/N)+E*exp(A/N)*AZ.
    # This retains the exprel' contribution through exprel+x*exprel'=exp.
    FZ=(ez*a+e*az)*exponential
    caps=dict(m=f.scalar(0),h=QZ,k=vz*Q+v*QZ+FZ*b+F*bz,
        e=ez*Q+e*QZ+b*bz*2+F*FZ,p=ez*Q+e*QZ+F*FZ)
    return caps,dict(original_R_E_cap=Q,original_R_E_Z_cap=QZ,original_F_N_cap=F,original_F_N_Z_cap=FZ,
        original_R2_integral_mass=c.mpf('.5'),original_R2_prime_integral_mass=c.mpf(1)/6,
        exact_F_N_Z_recipe='N*E_Z*expm1(A/N)+E*exp(A/N)*A_Z',
        full_signed_remainder_derivatives_retained=True,nonlinear_mean_derivative_not_zeroed=True,
        bounded_exponential_only=True)


class OriginalRmWeightedZAveraging:
    mode='actual_Rm_genuine_original_Z_yZ_phase_primitives_and_own_rate_weighted_Z_IBP_with_direct_fallback'
    def __init__(self,dps=500,require_checked=True):
        self.upstream=previous.OriginalRmMixedYZPrimitives(dps);self.c=self.upstream.c;self.family=self.upstream.family
        self.hashes=dict(self.upstream.hashes);self.cache={}
        fields.previous.bind(self.hashes,Path(__file__).name,sha(Path(__file__).name))
        self.baseline=json.loads(gzip.decompress((HERE/averaging.NAME).read_bytes()))
        if not self.baseline[averaging.GATE] or self.baseline['source_family']!=self.family:
            raise ValueError('Accepted same-source C0/Z local baseline required')
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted actual Rm weighted Z receipt required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def source_cell(self,label,left,right,N):
        N=phase.candidate_N(N);key=(label,left,right,N)
        if key in self.cache:return self.cache[key]
        p=self.upstream.upstream.upstream.phase
        op=p.upstream.upstream.upstream.owner(label).op;f=op.flow;c=op.c
        packet=self.upstream.query(label,left,right)
        if packet['source_family']!=self.family or packet['exact_common_P0_axial5'] is not op.P0 or packet['exact_source_Rm_factor'] is not op.Rm_factor:
            raise ValueError('Same actual mixed source owner, P0 and Rm required')
        if packet['original_mixed_source_record']['original_radius_Z_exact_zero'] is not True:
            raise ValueError('Original shared radius must be Z independent')
        with mp.workdps(c.dps+40):
            roots=packet['original_genuine_mixed_source_roots'];qr=packet['original_genuine_mixed_q_rows']
            velocity=packet['original_genuine_mixed_V_rows'];source=dict(q=qr[C0],roots=roots)
            got=phase.first.conditioned_first_jets(source,qr,p.upstream.dstar_log,c.mpf([0,1]))
            values=dict(packet['original_fixed_phi_A_B_mixed_values'])
            if got['values'] is not None:
                for name,row in got['values'].items():
                    cap,_=averaging.tighter(f,row,values[name]);values[name]=averaging.symmetric(f,cap)
            leading=previous.leading_mixed(roots['E'],velocity,values)
            remainder,proof=remainder_Z_caps(f,roots['E'],velocity,values,N)
            density=phase.densities.density_Z_kernels(roots['E'][C0],roots['E'][Z],velocity[C0],velocity[Z],values,N)
            if label not in p.mappers:p.mappers[label]=phase.RmRadiusPhase(op,self.family,p.parameter_family)
            endpoints={which:p.mappers[label].point(x,N) for which,x in (('left',left),('right',right))}
            for row in endpoints.values():
                if row.get('actual_Rm_radius_phase_Z_independent') is not True:
                    raise ValueError('Z-dependent radius phase is unsupported')
                if row['source_family']!=self.family or row['candidate_N']!=N or row['exact_source_Rm_factor'] is not op.Rm_factor:
                    raise ValueError('Same actual endpoint radius owner and N required')
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,
                source_geometry=packet['source_geometry'],exact_common_P0_axial5=op.P0,exact_source_Rm_factor=op.Rm_factor,
                genuine_original_mixed_source=packet,original_fixed_phi_A_B_mixed_values=values,
                original_conditioned_first_record=got['record'],original_five_signed_leading_density_C0_y_Z_yZ=leading,
                original_nonlinear_second_remainder_Z_caps=remainder,original_remainder_Z_proof=proof,
                original_direct_density_Z=density['Z_derivatives'],actual_endpoint_phases=endpoints,
                original_radius_phase_Z_exact_zero=True,original_leading_zero_phase_mean=True,
                phase_primitive_endpoints_exactly_zero=True,original_full_density_mean_not_zeroed=True,
                original_slow_yZ_not_total_spatial_yZ=True,**dict.fromkeys(fields.previous.OPEN,False))
            result['actual_endpoint_phase_primitive_Z_enclosures']={which:[phase_primitive_Z_enclosure(f,result,box)
                for box in row['phase_boxes']] for which,row in endpoints.items()}
        self.cache[key]=result;return result

    def contribution(self,label,N=257):
        N=phase.candidate_N(N)
        if N!=self.baseline['candidate_N']:
            raise ValueError('Comparison requires the same candidate N as the accepted local baseline')
        op=self.upstream.upstream.upstream.phase.upstream.upstream.upstream.owner(label).op
        f,c=op.flow,op.c;cells=[];total={key:f.scalar(0) for key in RATES};direct={key:f.scalar(0) for key in RATES}
        with mp.workdps(c.dps+40):
            xs=[c.exp(1) if x=='Rh' else c.mpf(x[0])/x[1] for x in PARTITION]
            for left,right,xl,xr in zip(PARTITION,PARTITION[1:],xs,xs[1:]):
                source=self.source_cell(label,left,right,N);require_phase_Z_zero(source)
                width=c.ln(xr/xl);suffix=c.mpf(0) if right=='Rh' else 1-c.ln(xr)
                distances={which:averaging.phase_distance(c,row['phase_boxes']) for which,row in source['actual_endpoint_phases'].items()}
                rows={};details={}
                for key,rate in RATES.items():
                    lead=source['original_five_signed_leading_density_C0_y_Z_yZ'][key]
                    cap,proof=averaging.weighted_cap(f,lead[Z],lead[YZ],source['original_nonlinear_second_remainder_Z_caps'][key],
                        width,suffix,rate,N,distances['left'],distances['right'])
                    simple=averaging.magnitude(f,source['original_direct_density_Z'][key])*proof['positive_own_rate_mass']*proof['downstream_suffix_decay']
                    used,route=averaging.tighter(f,cap,simple);rows[key]=used;total[key]+=used;direct[key]+=simple
                    details[key]=dict(original_Z_IBP=proof,weighted_Nminus2_Z_cap=cap,original_direct_Z_cap=simple,
                        used_valid_Z_cap=used,used_route=route)
                cells.append(dict(source=source,source_log_width=width,source_suffix_to_Rh=suffix,
                    actual_cell_to_Rh_Z_caps=rows,original_own_rate_Z_averaging=details))
            baseline=self.baseline['frames'][label]['actual_averaged_whole_patch_local_integral_C0_Z']
            output={};comparisons={}
            for key in RATES:
                old=averaging.read_saved_source(f,baseline[key][1]);oldcap=averaging.magnitude(f,old)
                used,route=averaging.tighter(f,total[key],oldcap)
                output[key]=[averaging.read_saved_source(f,baseline[key][0]),averaging.symmetric(f,used)]
                reduction=None if oldcap.zero or used.zero else (oldcap.scale-used.scale).evaluate()+c.ln(oldcap.coefficient)-c.ln(used.coefficient)
                comparisons[key]=dict(accepted_direct_whole_Z_cap=oldcap,new_cellwise_valid_Z_cap=total[key],
                    used_Z_cap=used,used_route=route,log_absolute_upper_reduction=reduction,
                    strict_absolute_upper_reduction=(not oldcap.zero and used.zero) or reduction is not None and ep(reduction)[0]>0,
                    accepted_C0_bound_unchanged=True)
        return dict(source_family=self.family,source_frame=label,candidate_N=N,exact_partition=PARTITION,
            actual_source_weighted_Z_cells=cells,actual_averaged_whole_patch_local_integral_C0_Z=output,
            actual_direct_new_partition_Z_caps=direct,actual_Z_comparisons=comparisons,
            strict_Z_reductions=sum(row['strict_absolute_upper_reduction'] for row in comparisons.values()),
            endpoint_retaining_genuine_Z_Nminus2_averaging_installed=True,
            actual_spatial_mixed_yZ_chain_installed=False,real_finite_N_Rm_incoming_correction_is_still_unsupplied=True,
            exact_common_P0_axial5=op.P0,source_integral_caps_not_five_moment_closure=True,
            **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalRmWeightedZAveraging(require_checked=False);frames={}
    for label in ('0','.5'):
        frames[label]=owner.contribution(label)
        print('Actual original Rm weighted Z averaging',label,'strict reductions',frames[label]['strict_Z_reductions'],flush=True)
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        frames=previous.serialized(frames),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        **dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()
