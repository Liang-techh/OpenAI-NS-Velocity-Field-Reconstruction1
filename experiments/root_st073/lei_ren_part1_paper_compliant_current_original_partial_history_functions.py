"""Genuine radial-variable reference/O2 affine histories and endpoint jets.

Exact rational subintervals reuse live entire-cell function/phase covers,
with newly computed masses to the requested endpoint. Reference factors
stay canonical until export; O2 hulls are already normalized. Unknown
upstream functions and the separate original pressure datum remain open.
"""
from dataclasses import dataclass
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import sympy as sy
import lei_ren_part1_paper_compliant_current_original_route_function_evaluator as current

O2=current.O2;reference=current.reference;leaves=current.leaves;precise=current.precise
transport=current.transport;HERE,PREFIX,sha,ep=current.HERE,current.PREFIX,current.sha,current.ep
NAME=PREFIX+'current_original_partial_history_functions.json.gz'
RECEIPT=PREFIX+'current_original_partial_history_functions_check.json'
GATE='actual_original_radial_variable_reference_O2_affine_histories_and_endpoint_jets_installed'


@dataclass(frozen=True)
class OriginalPartialHistoryFrame:
    left: object
    right: object
    Z: object
    N: int
    forcing: dict
    decay: dict
    record: dict


class OriginalPartialHistoryFunctions:
    mode='original_radial_history_partial'
    def __init__(self,evaluator=None):
        self.current=current.OriginalRouteFunctionEvaluator() if evaluator is None else evaluator
        if type(self.current) is not current.OriginalRouteFunctionEvaluator:
            raise TypeError('Genuine issued original route evaluator required')
        self.family=self.source_family=self.current.family;self.c=self.current.c
        self.source_graph_sha256=self.current.source_graph_sha256
        self.hashes=dict(self.current.hashes);self.seeds={};self.cache={};self.issued={};self.jet_cache={}
        accepted=json.loads((HERE/current.RECEIPT).read_bytes())
        if not accepted.get('all_passed') or not accepted.get(current.GATE) or accepted['source_family']!=self.family:
            raise ValueError('Accepted same-source original route evaluator required')
        self.bind_hashes({**accepted['input_hashes'],current.RECEIPT:sha(current.RECEIPT)})
        self.kernel_identities=self.partial_kernel_identities()
        self.bind_hashes({Path(__file__).name:sha(Path(__file__).name)})

    def bind_hashes(self,closure):
        for name,digest in closure.items():
            if sha(name)!=digest:raise ValueError('Original partial history dependency changed: '+name)
            if name in self.hashes and self.hashes[name]!=digest:raise ValueError('Original partial history closures disagree: '+name)
            self.hashes[name]=digest

    def coordinate(self,value):
        q=precise.rational(value)
        if not -5<=q<=1:raise ValueError('Exact native reference/O2 coordinate in[-5,1] required')
        return q

    def partial_kernel_identities(self):
        g=self.current.built['graph'];Y=sy.Symbol('requested_native_right',real=True);records=[]
        for node,owner in self.current.integrals.items():
            row=g.nodes[node];contract=owner.require_integral(row)
            if owner is self.current.reference:
                source=self.current.provider.rows['Rh_reference',contract['source_role']];B=sy.Integer(0)
            else:
                source=owner.sources[node];B=owner.bounds[contract['route_index']][1]
            expression,F=reference.symbolic_node(self.current.provider,row['integrand'],source)
            x=sy.Symbol(row['variable'],real=True);rate=transport.RATES[contract['key']]
            r=sy.Rational(rate.numerator,rate.denominator)
            if sy.simplify(sy.exp(r*(B-Y))*expression-sy.exp(-r*(Y-x))*F)!=0:
                raise ValueError('Original density/kernel cannot be localized to requested endpoint')
            records.append(dict(original_integral_node=node,original_source_node=contract['source_node'],
                key=contract['key'],ordinary_order=contract['ordinary_order'],original_kernel_endpoint=str(B),
                exact_requested_endpoint_kernel_identity_passed=True,original_native_Jacobian_applied_once=1,
                exact_partial_definition='integral_a^b exp(-rate*(b-x))*same_issued_original_density(x,Z,N) dx',
                endpoints_and_kernel_Z_independent=True))
        return records

    def seed(self,z,N):
        key=z,N
        if key not in self.seeds:
            live=self.current.history_frame(Z=str(z),N=N);s=self.current.settings
            ref=self.current.reference.frame(Z=str(z),N=N,count=s['reference_count'],bits=s['reference_bits'])
            slope=self.current.O2.frame(Z=str(z),N=N,count=s['O2_count'],bits=s['O2_bits'])
            if (ref.Z,ref.N)!=(z,N) or (slope.Z,slope.N)!=(z,N) or (live.Z,live.N)!=(z,N):
                raise ValueError('One original Z and candidate N required')
            if ref.owner.family!=self.family or slope.record['source_family']!=self.family:
                raise ValueError('Same original live source family required')
            if (ref.record['source_graph_sha256']!=self.source_graph_sha256 or
                slope.record['source_graph_sha256']!=self.source_graph_sha256 or
                self.current.reference.provider is not self.current.provider or
                self.current.O2.provider is not self.current.provider or
                self.current.reference.issued.get(id(ref)) is not ref or
                self.current.O2.issued.get(id(slope)) is not slope):
                raise ValueError('Issued live original source graph/parameter owners required')
            self.bind_hashes(self.current.reference.hashes);self.bind_hashes(self.current.O2.hashes)
            self.seeds[key]=ref,slope,live
        return self.seeds[key]

    def mass(self,left,right,endpoint,rate,ctx=None):
        """Stable positive mass even for exact subcells below decimal precision."""
        c=self.c if ctx is None else ctx
        width=right-left;distance=endpoint-right
        if not width>0 or distance<0:raise ValueError('Positive subcell below requested endpoint required')
        w=c.mpf(int(width.p))/int(width.q);d=c.mpf(int(distance.p))/int(distance.q)
        r=c.mpf(rate.numerator)/rate.denominator
        value=w if not rate else c.exp(-r*d)*(-c.expm1(-r*w))/r
        if ep(value)[0]<=0:raise ArithmeticError('Stable original positive partial mass lost')
        return value

    def subinterval(self,left,right,*,Z,N):
        a=self.coordinate(left);b=self.coordinate(right)
        if a>b:raise ValueError('Ordered native interval required')
        z=self.current.request(Z,N);key=a,b,z,N
        if key in self.cache:return self.cache[key]
        c=self.c;forcing={(k,j):c.mpf(0) for k in transport.RATES for j in ('C0','Z')};decay={};records=[]
        with mp.workdps(c.dps+40):
            width=b-a;w=c.mpf(int(width.p))/int(width.q)
            for density_key,rate in transport.RATES.items():
                decay[density_key]=c.exp(-(c.mpf(rate.numerator)/rate.denominator)*w)
            if a<b:
                ref,slope,live=self.seed(z,N);ro=ref.owner;atlas=ro.atlas
                total={order:{k:{j:atlas.scalar(0) for j in ('C0','Z')} for k in transport.RATES}
                    for order in reference.whole.points.exact.ORDERS}
                for item in ref.cells:
                    parent_left,parent_right=map(sy.Rational,item['source']['exact_reference_cell'])
                    left_part=max(a,parent_left);right_part=min(b,parent_right)
                    if not left_part<right_part:continue
                    # The same live owner retains the defining source program
                    # and its cached whole-cell coefficient/phase cover.
                    query=ro.cell(parent_left,parent_right,N=N,bits=ref.bits)
                    if query['record']['source_family']!=self.family or query['record']['candidate_N']!=N:
                        raise ValueError('Original live reference source/phase cell required')
                    masses={}
                    for k,rate in transport.RATES.items():
                        mass=self.mass(left_part,right_part,b,rate,ro.ctx);masses[k]=mass
                        for order in total:
                            for j in ('C0','Z'):
                                value=query['coefficients'][order][k][j]
                                if value.ctx is not ro.ctx or value.scale.bases is not atlas.bases or value.ledger is not atlas.ledger:
                                    raise ValueError('One canonical reference factor atlas required')
                                total[order][k][j]=atlas.add(total[order][k][j],value*mass)
                    records.append(dict(chart='Rh_reference',exact_parent_cell=[str(parent_left),str(parent_right)],
                        exact_restricted_cell=[str(left_part),str(right_part)],requested_kernel_endpoint=str(b),
                        exact_positive_own_rate_masses=masses,original_live_source_record=query['record'],
                        entire_source_and_true_phase_cover_restricted_conservatively=True,
                        original_exact_N_dependent_coefficients_and_pressure_errors_retained=True,
                        canonical_reference_factors_retained_until_combination=True))
                for k in transport.RATES:
                    for j in ('C0','Z'):
                        value=atlas.add(total[-1][k][j]*(ro.ctx.mpf(1)/N),total[-2][k][j]*(ro.ctx.mpf(1)/N**2))
                        forcing[k,j]=self.current.ordinary(leaves.base.conditioned.bounded_value(value))
                for item in slope.record['genuine_source_cell_phase_integral_records']:
                    parent_left,parent_right=map(sy.Rational,item['exact_y_cell'])
                    left_part=max(a,parent_left);right_part=min(b,parent_right)
                    if not left_part<right_part:continue
                    phase=item['actual_original_phase_origin_and_increment']
                    if (item['source']['source_family']!=self.family or phase['explicit_candidate_N']!=N or
                        phase['source_graph_sha256']!=self.source_graph_sha256 or
                        tuple(map(sy.Rational,item['source']['exact_Z_range']))!=(z,z)):
                        raise ValueError('Original live O2 source/phase cell required')
                    masses={}
                    for k,rate in transport.RATES.items():
                        mass=self.mass(left_part,right_part,b,rate);masses[k]=mass
                        for j in ('C0','Z'):
                            density=item['five_signed_C0_Z_density_hulls'][j][k]
                            if density.ctx is not c or hasattr(density,'scale'):
                                raise TypeError('Live normalized O2 density range required')
                            forcing[k,j]+=density*mass
                    records.append(dict(chart='O2_slope',original_route_index=item['route_index'],
                        exact_parent_cell=[str(parent_left),str(parent_right)],exact_restricted_cell=[str(left_part),str(right_part)],
                        requested_kernel_endpoint=str(b),exact_positive_own_rate_masses=masses,
                        original_live_source_record=item['source'],
                        actual_original_phase_cover=item['actual_original_phase_origin_and_increment'],
                        whole_parent_signed_C0_Z_density_hulls=item['five_signed_C0_Z_density_hulls'],
                        entire_source_and_true_phase_cover_restricted_conservatively=True,
                        original_quad_pressure_and_correlated_Z_errors_retained=True,
                        old_route_endpoint_masses_and_preintegrated_values_not_reused=True))
        record=dict(source_family=self.family,source_graph_sha256=self.source_graph_sha256,
            exact_native_interval=[str(a),str(b)],original_Z_exact=str(z),explicit_candidate_N=N,
            directed_partial_C0_Z_forcing={k:{j:forcing[k,j] for j in ('C0','Z')} for k in transport.RATES},
            directed_homogeneous_multipliers=decay,actual_partial_source_cell_records=records,
            exact_requested_endpoint_kernel_bindings=self.kernel_identities,
            unknown_incoming_functions_retained_at_exact_left_endpoint=True,
            prefix_original_incoming_nodes=({str(node):dict(key=k,ordinary_order=j) for node,(k,j) in self.current.boundaries.items()} if a==-5 else None),
            exact_zero_width_integral=a==b,separate_original_P0_and_P0_Z_unchanged=True,
            reference_exact_N_dependent_order_minus1_and_minus2_functions_retained=True,
            O2_full_original_DAG_density_hulls_at_actual_fixed_N_not_all_N_coefficients=True,
            fixed_nonzero_Z_radial_query_not_whole_Z_or_full_history=True,
            ordinary_Z_kernel_and_endpoints_independent_no_Leibniz_boundary_terms=True,
            live_source_owners_used_not_saved_report_hydration=True,
            full_scalar_control_evaluator_compatibility_installed=False,actual_original_upstream_history_installed=False)
        frame=OriginalPartialHistoryFrame(a,b,z,N,forcing,decay,record)
        self.cache[key]=frame;self.issued[id(frame)]=frame;self.current.unchanged();return frame

    def prefix(self,Y,*,Z,N):return self.subinterval(-5,Y,Z=Z,N=N)

    def endpoint_jet(self,frame,*,chart=None):
        self.current.unchanged()
        if type(frame) is not OriginalPartialHistoryFrame or self.issued.get(id(frame)) is not frame:
            raise ValueError('Issued original partial history frame required')
        if frame.right==0 and chart is None:
            raise ValueError('Specify original one-sided chart at the reference/O2 seam; no jet join admitted')
        chosen=('Rh_reference' if frame.right<0 else 'O2_slope') if chart is None else chart
        if chosen not in leaves.SUPPORTED or (chosen=='Rh_reference' and not -5<=frame.right<=0) or (chosen=='O2_slope' and not 0<=frame.right<=1):
            raise ValueError('Original endpoint chart/domain required')
        key=id(frame),chosen
        if key in self.jet_cache:return self.jet_cache[key]
        c=self.c;source={};forcing={};coefficient={}
        with mp.workdps(c.dps+40):
            for k,rate in transport.RATES.items():
                r=c.mpf(rate.numerator)/rate.denominator
                coefficient[k]=-r*frame.decay[k]
                for j in ('C0','Z'):
                    row=self.current.provider.rows[chosen,'density_'+k+'_'+j]
                    source[k,j]=self.current.source(row,coordinate=str(frame.right),Z=str(frame.Z),N=frame.N)
                    forcing[k,j]=source[k,j]-r*frame.forcing[k,j]
        result=dict(source_family=self.family,exact_native_endpoint=str(frame.right),endpoint_chart=chosen,
            original_Z_exact=str(frame.Z),explicit_candidate_N=frame.N,
            actual_point_C0_Z_density={k:{j:source[k,j] for j in ('C0','Z')} for k in transport.RATES},
            directed_endpoint_y_derivative_forcing={k:{j:forcing[k,j] for j in ('C0','Z')} for k in transport.RATES},
            directed_endpoint_y_derivative_incoming_multipliers=coefficient,
            defining_equation='H_y=original_density-rate*H; H_yZ=original_density_Z-rate*H_Z',
            fundamental_theorem_of_calculus_on_original_native_coordinate=True,
            actual_source_phase_derived_at_endpoint_not_free_angle=True,
            conditional_on_unknown_left_history=True,higher_velocity_stress_join_not_admitted=True)
        self.jet_cache[key]=result;return result

    def apply_assumed_boundary(self,frame,incoming):
        self.current.unchanged()
        if type(frame) is not OriginalPartialHistoryFrame or self.issued.get(id(frame)) is not frame:
            raise ValueError('Issued original partial history frame required')
        if set(incoming)!=set(frame.forcing) or any(not hasattr(v,'_mpi_') or v.ctx is not self.c or hasattr(v,'scale') for v in incoming.values()):
            raise TypeError('All ten directed same-context assumed boundary enclosures required')
        with mp.workdps(self.c.dps+40):
            return {(k,j):frame.decay[k]*incoming[k,j]+value for (k,j),value in frame.forcing.items()}


@precise.phase.native.inlet.source_precision
def run(*,return_live=False):
    began=time.monotonic();owner=OriginalPartialHistoryFunctions();frames=[];jets=[]
    for N in (160,257):
        for Y in ('-5','-2.337','0','.12','.13','.14','.15','.731','1'):
            frame=owner.prefix(Y,Z='.37',N=N);frames.append(frame)
            if Y in ('-2.337','.731'):jets.append(owner.endpoint_jet(frame))
        print('Original radial-variable prefix histories and endpoint jets:',N,flush=True)
    pairs=(('-5','-2.337'),('-2.337','0'),('0','.137'),('.137','.731'),('.12','.15'),('.731','1'),('.731','.731'),('.12',str(sy.Rational(3,25)+sy.Rational(1,10**1000))))
    for a,b in pairs:frames.append(owner.subinterval(a,b,Z='.37',N=160))
    result=dict(**{GATE:True},source_family=owner.family,source_graph_sha256=owner.source_graph_sha256,mode=owner.mode,
        actual_radial_variable_history_frames=[frame.record for frame in frames],actual_endpoint_y_Z_jet_records=jets,
        exact_partial_kernel_identities=owner.kernel_identities,
        arbitrary_exact_radial_query_and_subintervals_installed=True,
        actual_original_upstream_history_installed=False,full_scalar_control_evaluator_compatibility_installed=False,
        full_17_chart_source_or_24_cell_integral_oracle_installed=False,actual_five_controls_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(precise.phase.packets.OPEN,False),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,
        scope='Original reference/O2 affine history functions of arbitrary exact radial coordinate at fixed nonzero Z/current N, stable exact subinterval masses and true endpoint y/Z jets. Unknown upstream functions and separate P0 retained. No actual full history, continuous-Z terminal closure, global N, recursion or corrected NS.')
    raw=json.dumps(precise.encode(result),indent=2).encode()+b'\n'
    (HERE/NAME).write_bytes(gzip.compress(raw,compresslevel=9,mtime=0))
    return (result,owner,frames) if return_live else result


if __name__=='__main__':run()
