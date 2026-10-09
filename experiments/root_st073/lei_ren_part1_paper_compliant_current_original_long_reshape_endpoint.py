"""Actual R110 histories transported through the original finite long reshape.

The accepted finite terminal kernel backend is reused with the new anchored
B function. Incoming histories, decays and separate analytic P0 remain.
Native source frames0,.5 are conditional on admitted upstream inlets.
"""
import json
import math
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_second_switch_R110 as second
import lei_ren_part1_paper_compliant_current_original_reshape_terminal_kernels as kernels
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_inner_bridge_profiles import logarithm

fields,endpoint=second.fields,second.endpoint
prior,base,ep=second.prior,second.base,second.ep
HERE,PREFIX,sha=second.HERE,second.PREFIX,second.sha
NAME=PREFIX+'current_original_long_reshape_endpoint.json'
RECEIPT=PREFIX+'current_original_long_reshape_endpoint_check.json'
GATE='original_actual_two_frame_long_reshape_terminal_history_functions_installed'
NAMES=('theta','theta_z','pressure','swirl','mean','axial')


def source_bindings():
    from lei_ren_part1_paper_compliant_actual_switch_mixed_C4 import assignment_source_bindings
    bindings=assignment_source_bindings('long_reshape_profiles','inputs',{
        'normalized':"dict(theta=moments[MTH]/(phi*2),theta_z=moments[MTHZ]/(phi*2),pressure=moments[MP]/square(phi),swirl=moments[MZT]['swirl']/square(phi),mean=moments[MZ],axial=moments[MZT]['axial'])"})
    bindings.update(assignment_source_bindings('long_reshape_profiles','evaluate',{
        'mean':"inp['moments']['mean']*theta+inp['v']*(1-theta)",
        'axial':"inp['moments']['axial']*theta+square(inp['v'])*(1-theta)",
        'moment_shapes':"dict(theta=inherited['theta']+kernels['theta'],theta_z=inherited['theta_z']+inp['v']*kernels['theta'],pressure=inherited['pressure']+kernels['pressure'],swirl=inherited['swirl']+kernels['swirl'],mean=mean,axial=axial)"}))
    bindings.update(assignment_source_bindings('current_original_second_switch_R110','inlet',{
        'B':'G*(-Lambda)+logarithm(phi)+logarithm(1+z*z)+c.ln(2*R)/2'}))
    return dict(passed=True,actual_inlet_and_history_assignments=bindings,
        original_kernel_and_frozen_T_bindings=kernels.source_bindings(),
        phase_zero_function_identity='B-logq-logCstar=-Lambda*G+log(phi110)+log(220)/2=log(Utheta110)',
        identity_uses_anchored_F0_definition_not_interval_overlap=True,
        terminal_log_Utheta='T/10-logCstar-log(1+Z^2); B*(1-sigma(1)) vanishes exactly',
        original_V110_and_P0_source_namespace_preserved=True)


def reciprocal(flow,rows):
    """Factored ordinary Taylor reciprocal, with one positive source base."""
    if len(rows)!=6 or any(v.scale.bases is not flow.logs or v.ledger is not flow.ledger for v in rows):
        raise ValueError('Same-source ordinary denominator rows0..5 required')
    cover=flow.ordinary_cover(rows[0])
    if ep(cover)[0]<=0:raise ValueError('Actual positive angular base required; no selected denominator')
    loglower=flow.c.ln(flow.c.mpf(ep(cover)[0]))
    div=lambda row:row.positive_divide(rows[0],loglower)
    result=[div(flow.scalar(1))]
    for n in range(1,6):result.append(-div(sum((rows[j]*result[n-j] for j in range(1,n+1)),flow.scalar(0))))
    return result


def radial_Q(flow,Z,delta,V,M):
    c=flow.c;z=IntervalTaylor.variable(c,Z,5);d=1-z*z;L=1-z*z*delta
    if ep(L[0])[0]<=0:raise ValueError('Positive original radial denominator required')
    Mz=[M[n+1]*(n+1) for n in range(5)]+[flow.scalar(0)]
    return flow.multiply(flow.add(flow.scale(flow.multiply(flow.jet(z),V),2),
        flow.scale(flow.multiply(flow.jet(z),M),-(1-delta)),
        flow.scale(flow.multiply(flow.jet(d),Mz),-1)),flow.jet(L.reciprocal()))[:5]


class ActualLongReshapeEndpoint:
    def __init__(self,flow,Z,delta,phi,V,histories,B,T,logC,p0,*,L=4096):
        self.flow=f=flow;self.c=c=f.c;self.Z=c.mpf(Z);self.delta=c.mpf(delta)
        if ep(self.Z)[0]<-1 or ep(self.Z)[1]>1 or ep(self.delta)[0]<0 or ep(self.delta)[1]>=1:
            raise ValueError('Original axial coordinate and radial parameter required')
        if set(histories)!=set(second.moments.RATES):raise ValueError('All actual R110 six histories required')
        for row in (phi,V,*histories.values()):
            if len(row)!=6 or any(v.scale.bases is not f.logs or v.ledger is not f.ledger for v in row):
                raise ValueError('Actual R110 source basis/ledger and ordinary rows0..5 required')
        for jet in (B,p0):
            if not isinstance(jet,IntervalTaylor) or jet.ctx is not c or jet.order!=5:
                raise ValueError('Same-context B/P0 ordinary source jets0..5 required')
        self.phi=phi;self.V=V;self.histories=histories;self.B=B;self.p0=p0
        self.T=c.mpf(T);self.logC=c.mpf(logC);self.L=L
        self.backend=kernels.TerminalKernelEnclosures(c,B,self.T)
        invphi=reciprocal(f,phi);invphi2=f.multiply(invphi,invphi)
        self.inlets=dict(theta=f.scale(f.multiply(histories['H'],invphi),c.mpf('.5')),
            theta_z=f.scale(f.multiply(histories['K'],invphi),c.mpf('.5')),
            pressure=f.multiply(histories['C'],invphi2),swirl=f.multiply(histories['B'],invphi2),
            mean=histories['M'],axial=histories['A'])

    def import_kernel(self,row):
        # Directed change of arithmetic coordinates, not a selected source.
        # log(T) powers stay as their complete directed linear form in offset.
        f=self.flow
        return prior.ScaledEnclosure(prior.FormalScale(f.logs,offset=row.scale.evaluate()),row.coefficient,f.ledger)

    def decay_rows(self,k,m):
        c=self.c;f=self.flow
        bell=IntervalTaylor(c,[c.mpf(0)]+[v*m for v in self.B.coefficients[1:]]).exp()
        factor=f.factor((0,0,0,0,0),-c.mpf(k)*self.T+m*self.B[0])
        return f.scale(f.jet(bell),factor)

    def evaluate(self):
        f=self.flow;c=self.c
        raw={name:self.backend.evaluate(name,L=self.L) for name in kernels.KINDS}
        K={name:[self.import_kernel(v) for v in row['coefficients']] for name,row in raw.items()}
        decays={name:self.decay_rows(k,m) for name,(k,m,_) in kernels.KINDS.items()}
        inherited={name:f.multiply(self.inlets[name],decays['theta' if name=='theta_z' else name])
            for name in ('theta','theta_z','pressure','swirl')}
        dmean=f.factor((0,0,0,0,0),-self.T);V2=f.multiply(self.V,self.V)
        incoming_mean=f.scale(f.add(self.inlets['mean'],f.scale(self.V,-1)),dmean)
        incoming_axial=f.scale(f.add(self.inlets['axial'],f.scale(V2,-1)),dmean)
        shapes=dict(theta=f.add(inherited['theta'],K['theta']),
            theta_z=f.add(inherited['theta_z'],f.multiply(self.V,K['theta'])),
            pressure=f.add(inherited['pressure'],K['pressure']),swirl=f.add(inherited['swirl'],K['swirl']),
            mean=f.add(self.V,incoming_mean),axial=f.add(V2,incoming_axial))
        z=IntervalTaylor.variable(c,self.Z,5);logq=logarithm(1+z*z)
        logU=IntervalTaylor(c,[-logq[0]+self.T/10-self.logC]+[-v for v in logq.coefficients[1:]])
        ratio=IntervalTaylor(c,[c.mpf(0)]+list(logU.coefficients[1:])).exp()
        ratio2=IntervalTaylor(c,[c.mpf(0)]+[2*v for v in logU.coefficients[1:]]).exp()
        R=f.factor((0,0,0,0,0),self.T)*110
        U=f.factor((0,0,0,0,0),logU[0]);U2=U*U;P2=f.factor((0,1,0,0,0))
        Q=radial_Q(f,self.Z,self.delta,self.V,shapes['mean'])
        dressed={name:f.multiply(row,f.jet(ratio if name in ('theta','theta_z') else ratio2))
            for name,row in shapes.items() if name in ('theta','theta_z','pressure','swirl')}
        rroot=f.factor((0,0,0,0,0),self.T/2)*c.sqrt(110)
        axis=f.scale(f.jet(self.p0),P2);increment=f.scale(dressed['pressure'],U2*c.mpf('.5'))
        moment=dict(Mtheta=f.scale(dressed['theta'],R*rroot*U*c.sqrt(2)),
            Mtheta_z=f.scale(dressed['theta_z'],R*rroot*U*c.sqrt(2)),Mz=f.scale(shapes['mean'],R),
            Mztheta=f.add(f.scale(shapes['axial'],R),f.scale(dressed['swirl'],-R*U2*c.mpf('.5'))),Mp=increment)
        source=dict(theta=[f.scalar(1)]+[f.scalar(0)]*5,theta_z=self.V,
            pressure=[f.scalar(1)]+[f.scalar(0)]*5,swirl=[f.scalar(1)]+[f.scalar(0)]*5,
            mean=self.V,axial=V2)
        rates=dict(theta=c.mpf('1.6'),theta_z=c.mpf('1.6'),pressure=c.mpf('.2'),swirl=c.mpf('1.2'),mean=1,axial=1)
        ODE={name:f.add(source[name],f.scale(row,-rates[name])) for name,row in shapes.items()}
        return dict(actual_R110_normalized_inlet_shapes=self.inlets,
            actual_terminal_normalized_six_history_shapes=shapes,
            actual_terminal_log_radius_ODE_rows=ODE,
            actual_full_finite_kernels={name:kernels.record(row) for name,row in raw.items()},
            original_nonzero_incoming_decay_rows=decays,actual_inherited_terminal_history_contributions=inherited,
            original_mean_axial_decay=dmean,actual_mean_incoming_memory=incoming_mean,actual_axial_incoming_memory=incoming_axial,
            actual_Q_axial4_coefficients=Q,original_P0_normalized_axial5=f.jet(self.p0),
            terminal_log_Utheta_axial5=list(logU.coefficients),
            terminal_log_Utheta_over_Pstar_axial5=[logU[0]-f.logs[1]/2]+list(logU.coefficients[1:]),
            terminal_Utheta_derivative_ratios=list(ratio.coefficients),terminal_Utheta_squared_derivative_ratios=list(ratio2.coefficients),
            physical_velocity_axial_coefficients=dict(Ur=[v*rroot*(1/c.sqrt(2)) for v in Q],Utheta=f.scale(f.jet(ratio),U),Uz=self.V),
            physical_pressure_axis_axial5=axis,physical_pressure_radial_increment_axial5=increment,
            physical_total_pressure_axial5=f.add(axis,increment),physical_cumulative_moment_axial5=moment,
            source_geometry=dict(log_R_over110=self.T,original_fixed_T=self.T,Rsh=R,frozen_T_Z_exact_zero=True),
            actual_V110_preserved_exactly=True,separate_original_P0_preserved=True,
            original_nonzero_inlet_memory_and_body_tail_errors_retained=True,
            original_kernel_coordinate_change_is_directed_enclosure_only=True,
            amplitude_derivative_dressing_applied_once=True,physical_positive_factors_not_materialized=True,
            Q_ordinary_orders=list(range(5)),other_ordinary_orders=list(range(6)))


class OriginalLongReshapeEndpoint:
    mode='genuine_original_two_frame_actual_finite_long_reshape_terminal_histories'
    def __init__(self,dps=500):
        self.upstream=second.OriginalSecondSwitchR110(dps);self.c=c=self.upstream.c
        self.hashes=dict(self.upstream.hashes);self.family=self.upstream.family;self.saved=None;self.owners={}
        for module in (second,kernels):
            for name in (module.NAME,module.RECEIPT):
                row=json.loads((HERE/name).read_bytes())
                if not row.get(module.GATE) or name==module.RECEIPT and not row.get('all_passed'):
                    raise ValueError('Accepted actual R110 and complete finite kernel receipts required')
                if row['source_family']!=self.family:raise ValueError('Same R110/reshape source family required')
                for path,digest in row['input_hashes'].items():fields.previous.bind(self.hashes,path,digest)
                fields.previous.bind(self.hashes,name,sha(name))
                if name==second.NAME:self.saved=row
        params=json.loads((HERE/kernels.PARAMS).read_bytes());self.params=params
        self.A=fields.previous.read_interval(c,params['A_upper']);self.T=400*self.A
        self.logC=fields.previous.read_interval(c,params['selected_logCstar'])
        join_name=PREFIX+'reference_join_bounds.json';join=json.loads((HERE/join_name).read_bytes())
        if fields.previous.read_interval(c,join['B_C2_upper'])._mpi_!=(2*self.A)._mpi_:
            raise ValueError('Same original B norm theorem and selected A required')
        accepted=json.loads((HERE/kernels.NAME).read_bytes())['actual_original_terminal_kernel_evaluation']
        stored=fields.previous.read_interval(c,accepted['original_frozen_T'])
        if not ep(stored)[0]<=ep(self.T)[0]<=ep(self.T)[1]<=ep(stored)[1]:
            raise ValueError('Original fixed T must contain the same source400*A expression')
        for name in (kernels.PARAMS,join_name,Path(__file__).name):fields.previous.bind(self.hashes,name,sha(name))
        self.bindings=source_bindings()

    def owner(self,label):
        if label in self.owners:return self.owners[label]
        source=self.saved['post_power_packets'].get(label)
        if not source:raise ValueError('Admitted actual R110 frames0,.5 only')
        packet=source[-1];row=packet['post_power_function_evaluation']
        if packet['source_family']!=self.family or packet['source_frame']!=label or row['geometry']['exact_fixed_radius']!=[110,1]:
            raise ValueError('Same actual source R110 inlet required')
        first_owner=self.upstream.upstream.owner(label);flow=first_owner.flow;c=self.c
        _,proof,_=self.upstream.upstream.upstream.owner(label)
        restore=lambda values:[endpoint.restore_row(flow,v) for v in values]
        actual={name:restore(values) for name,values in row['actual_fields'].items()}
        histories={name:restore(values) for name,values in row['actual_six_histories'].items()}
        B=IntervalTaylor(c,[fields.previous.read_interval(c,v) for v in packet['original_anchored_log_shape_B_axial5']])
        # These actual native source enclosures already satisfy the original
        # theorem. No cap/intersection supplies a B field value.
        for n in range(3):
            upper=ep(2*self.A/math.factorial(n))[1]
            if max(abs(v) for v in ep(B[n]))>upper:raise ValueError('Actual B C2 norm not source-certified')
        p0=IntervalTaylor(c,proof['original_P0_coefficients'][:6])
        records=self.upstream.upstream.upstream.upstream.fields.records
        delta=c.exp(fields.previous.read_interval(c,records['pressure_source']['compliant_source']['parameter_bounds']['log_delta']))
        self.owners[label]=ActualLongReshapeEndpoint(flow,proof['Z'],delta,actual['phi'],actual['V'],histories,B,self.T,self.logC,p0)
        return self.owners[label]

    def evaluate(self,label):
        with mp.workdps(self.c.dps+40):op=self.owner(label);value=op.evaluate()
        return dict(mode=self.mode,source_frame=label,source_family=self.family,
            actual_R110_B_log_axial5=list(op.B.coefficients),function_evaluation=fields.serialized(value),
            actual_B_C2_source_norm_verified_without_cap=True,
            same_original_selected_A_and_fixed_T=True,whole_axis_functions_installed=False,
            original_upstream_micro_function_providers_complete=False,no_original_ancestor_producers_or_full_checkers_executed=True)


def run():
    began=time.monotonic();owner=OriginalLongReshapeEndpoint()
    report=dict(**{GATE:True},source_family=owner.family,original_source_bindings=owner.bindings,
        original_fixed_T=owner.T,terminal_packets={label:owner.evaluate(label) for label in ('0','.5')},
        **dict.fromkeys(fields.previous.OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began,
        scope='Actual source R110 six histories and anchored B feed complete finite original long reshape at native0,.5. Nonzero inherited memory, kernel/body/tail/source errors, separate analytic P0 and physical factors remain. Conditional upstream inlets; whole-Z, reference/restoration/global closure and actual n-recursion stay open.')
    (HERE/NAME).write_text(json.dumps(base.encoded(report),indent=2)+'\n',encoding='utf8')
    print('Actual original R110 histories transported through complete finite long reshape',flush=True);return report


if __name__=='__main__':run()
