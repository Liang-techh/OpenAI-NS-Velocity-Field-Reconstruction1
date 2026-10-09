"""True frozen-macro whole-cell source and finite-N local correction.

R0=Ra exp(2hb), L=log(100/Ra)-2hb are the defining physical geometry.
The macro correction inlet remains unknown; background inlets are not it.
"""
import ast
import gzip
import json
from pathlib import Path
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_switch_finite_N as downstream

long=downstream.long;fields,base,ep=downstream.fields,downstream.base,downstream.ep
parameters,primitives,phase,bounds=downstream.parameters,downstream.primitives,downstream.phase,downstream.bounds
macro=downstream.first.endpoint.moments
HERE,PREFIX,sha=downstream.HERE,downstream.PREFIX,downstream.sha
NAME=PREFIX+'current_original_macro_finite_N.json.gz'
RECEIPT=PREFIX+'current_original_macro_finite_N_check.json'
GATE='current_original_frozen_macro_whole_source_and_finite_N_local_driver_installed'
RATES,PARTITION,C0,Z=downstream.RATES,downstream.PARTITION,downstream.C0,downstream.Z


def source_bindings():
    tree=ast.parse(Path(macro.__file__).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='evaluate'
        and any(isinstance(q,ast.Assign) and any(ast.unparse(t)=='source_polys' for t in q.targets) for q in ast.walk(n)))
    equations=("source_polys=dict(H=phi.scale(2),M=V,K=(phi*V).scale(2),A=V*V,B=phi*phi,C=phi*phi)",
        "actual=f.add(inherited,signed,f.error_rows(integral_error))",
        "derivatives[name]=f.add(source_end[name],f.scale(actual,-rate))")
    for text in equations:
        wanted=ast.dump(ast.parse(text).body[0])
        if sum(ast.dump(n)==wanted for n in ast.walk(fn) if isinstance(n,ast.Assign))!=1:
            raise ValueError('Current original cumulative macro recipe changed: '+text)
    field_tree=ast.parse(Path(macro.fields.__file__).read_text(encoding='utf8'))
    flow=next(n for n in field_tree.body if isinstance(n,ast.ClassDef) and n.name=='MacroFlow')
    field_fn=next(n for n in flow.body if isinstance(n,ast.FunctionDef) and n.name=='evaluate')
    field_equations=("ell=self.add(self.ell_in,*[self.scale(row,self.h*I[1,l]*(-c.mpf('.5'))) for l,row in enumerate(self.d)])",
        "product=self.multiply(self.q,self.drive[part][j])", "parts[part]=self.scale(self.add(*sums),-scale)",
        "scale=self.factor((1,int(part=='pressure'),int(part=='swirl'),0,0))")
    for text in field_equations:
        wanted=ast.dump(ast.parse(text).body[0])
        if sum(ast.dump(n)==wanted for n in ast.walk(field_fn) if isinstance(n,ast.Assign))!=1:
            raise ValueError('Original coupled macro phi/V source changed: '+text)
    original=ast.parse((HERE/(PREFIX+'actual_bridge_integrals.py')).read_text(encoding='utf8'))
    geometry={'R0':'Ra*exp(2hb)','L':'log(100/Ra)-2hb'}
    call=next(n for n in ast.walk(original) if isinstance(n,ast.keyword) and n.arg=='exact_macro_source')
    for key,text in geometry.items():
        got=[n.value.value for n in call.value.keywords if n.arg==key and isinstance(n.value,ast.Constant)]
        if got!=[text]:raise ValueError('Original exact macro geometry changed: '+key)
    phase_origin={PREFIX+'current_generic_loop_function_sources.py':('phi','fractional_part(N*log(R/r_minus))'),
        PREFIX+'current_generic_shear_loop_domain.py':('r_minus','Ra*exp(hb*s_c/2)')}
    for name,(key,text) in phase_origin.items():
        origin=ast.parse((HERE/name).read_text(encoding='utf8'))
        got=[n.value.value for n in ast.walk(origin) if isinstance(n,ast.keyword) and n.arg==key and isinstance(n.value,ast.Constant)]
        if text not in got:raise ValueError('Original phase origin changed: '+name)
    return dict(actual_macro_history_assignments=equations,original_coupled_macro_field_assignments=field_equations,
        original_defining_macro_geometry=geometry,
        original_phase_and_positive_interior_origin_bindings=phase_origin,
        current_macro_source_controls=macro.fields.bridge.bridge_control_bindings(),
        unchanged_switch_recovery_and_general_q_bindings=downstream.source_bindings(),
        source_equation='phi_y=-hb*Dbar*phi/2; V_y=-hb*(phi/barphi)*drive(R,Z)',
        actual_geometry_not_selected_hb_cap=True,independent_source_geometry_and_ODE_receipt_required=True)


class WholeMacroFlow(macro.fields.MacroFlow):
    """Same live algebra, original functions, exact geometry on whole cells."""
    def __init__(self,flow,series):
        self.__dict__.update(flow.__dict__);self.original_flow=flow;self.series=series
        if series.flow is not flow:raise ValueError('One current macro/switch source flow required')
        self.length=self.scalar(self.Y)-self.h*2
        parameters.positive_source(self,self.length,'actual_positive_R0_R100_log_length')

    def coordinate(self,cell):
        if not isinstance(cell,tuple) or len(cell)!=2:raise ValueError('Two exact rational macro cell endpoints required')
        l,r=map(long.fraction,cell)
        if r<l:raise ValueError('Ordered macro source cell required')
        c=self.c;t=downstream.downstream.scalar_hull(c,c.mpf(l.numerator)/l.denominator,c.mpf(r.numerator)/r.denominator)
        return t,l==r==0,l==r==1

    def radius(self,t,root=False):
        p=.5 if root else 1
        return self.factor((0,0,0,p,0),p*self.Y*t)*self.series.scalar_series(2*p,1-t)[0]

    def geometry(self,cell):
        t,empty,terminal=self.coordinate(cell)
        R0=self.factor((0,0,0,1,0))*self.series.scalar_series(2,self.c.mpf(1))[0]
        R=self.scalar(100) if terminal else R0 if empty else self.radius(t)
        # S enters only ordinary polynomial weights; this outward source
        # cover is not a selected hb. Physical radius/measure retain hb.
        S=self.ordinary_cover(self.length)*t
        return S,R0,R,empty


def macro_cell(op,series,left,right):
    f,c=op.flow,op.flow.c;view=WholeMacroFlow(f,series);cell=(left,right)
    t,empty,terminal=view.coordinate(cell)
    evaluated=macro.ActualMacroMoments(view,op.inlet).evaluate(cell)
    _,R0,R,_=view.geometry(cell);rootR=f.scalar(10) if terminal else view.radius(t,root=True)
    if empty:rootR=f.factor((0,0,0,.5,0))*series.scalar_series(1,c.mpf(1))[0]
    phi,V=(evaluated['same_original_macro_field_functions'][name] for name in ('phi','V'))
    modes={power:[view.radial_power(R0,j)*view.radial_power(R,power-j) for j in range(3)] for power in (1,2)}
    D=f.add(*[f.scale(row,modes[1][j]) for j,row in enumerate(f.d)])
    dp=parameters.positive_source(f,D[0],'actual_frozen_macro_Dbar')
    a=f.scale(D,f.h)
    core=fields.IntervalTaylor(c,[f.ordinary_cover(row) for row in f.phi0])
    if ep(core[0])[0]<=0:raise ValueError('Current original core phi must be positive')
    inv_barphi=f.multiply(f.q,f.jet(core.reciprocal()));quotient=f.multiply(phi,inv_barphi)
    drives={part:f.add(*[f.scale(row,modes[2 if part=='swirl' else 1][j])
                        for j,row in enumerate(f.drive[part])]) for part in fields.PARTS}
    Vy=f.add(*[f.scale(f.multiply(quotient,row),-f.factor((1,int(part=='pressure'),int(part=='swirl'),0,0)))
               for part,row in drives.items()])
    return dict(chart='frozen_macro',phase=t,fields=dict(phi=phi,V=V),
        histories=evaluated['actual_six_moment_functions'],phi_y=f.scale(f.multiply(a,phi),-c.mpf('.5')),
        V_y=Vy,correlated_a_axial5=a,actual_Dbar_source=D,actual_positive_Dbar_proof=dp,
        original_current_drive_modes=drives,original_current_phi_over_barphi=quotient,
        radius=R,root_radius=rootR,window_length=view.length,
        physical_log_radius_over100=f.scalar(0) if terminal else f.scalar(-f.Y*(1-t))+f.h*(2*(1-t)),
        exact_log_radius_over_R0=view.length*t,actual_micro_exit_radius=R0,
        whole_source_geometry=evaluated['geometry'],original_cumulative_macro_evidence=evaluated['complete_integral_evidence'],
        source_interval_functions_not_endpoint_hulls=True,true_physical_measure='dy=(Y-2hb)*d_rho',
        actual_micro_background_inlet_preserved=True,real_finite_N_micro_exit_correction_supplied=False)


def weights(series,left,right,rate):
    f,c=series.flow,series.c;l,r=map(long.fraction,(left,right))
    if r<=l:raise ValueError('Positive actual macro cell width required')
    d=c.mpf((r-l).numerator)/(r-l).denominator;s=c.mpf((1-r).numerator)/(1-r).denominator
    L=f.scalar(f.Y)-f.h*2;width=L*d
    if not rate:return width,f.scalar(1),f.scalar(1)
    lam=c.mpf(rate.numerator)/rate.denominator
    decay=f.factor((0,0,0,0,0),-lam*f.Y*d)*series.scalar_series(2*lam,d)[0]
    tail=f.factor((0,0,0,0,0),-lam*f.Y*s)*series.scalar_series(2*lam,s)[0]
    mass=long.bounded_mass(f,f.ordinary_cover(width),lam,decay)
    return f.scalar(mass),decay,tail


def actual_phase_record(f,source,sc):
    lo,hi=ep(sc)
    if lo!=hi or not 0<lo<1:raise ValueError('Actual selected positive singleton s_c required')
    offset=f.h*(sc/2)
    return dict(actual_phase_definition='frac(N*((Y-2hb)*rho+hb*(2-s_c/2)))',
        actual_phase_Z_exact_zero=True,actual_selected_s_c=sc,original_positive_left_log_offset=offset,
        actual_unwrapped_log_radius_phase=source['exact_log_radius_over_R0']+f.h*2-offset)


class OriginalMacroFiniteN:
    mode='current_original_frozen_macro_source_and_complete_local_micro_exit_Rm_affine_driver'
    def __init__(self,dps=500,require_checked=True):
        self.downstream=downstream.OriginalSwitchFiniteN(dps);self.c=self.downstream.c;self.family=self.downstream.family
        self.moments=self.downstream.r100.upstream;self.hashes=dict(self.downstream.hashes);self.cache={}
        collar_name=PREFIX+'current_inner_exit_strict_collar.json';collar=json.loads((HERE/collar_name).read_bytes())
        origin=self.moments.fields
        for key,want in (('actual_five_defect_family_sha256',self.family),('implicit_source_sha256',origin.source),
                         ('datum_enclosure_sha256',origin.datum)):
            if collar[key]!=want:raise ValueError('Same current leading/collar phase-origin source required')
        self.collar=collar['explicit_current_inner_exit_strict_collar']
        self.sc=fields.previous.read_interval(self.c,self.collar['selected_first_phase_endpoint'])
        if ep(self.sc)[0]!=ep(self.sc)[1] or not 0<ep(self.sc)[0]<1:
            raise ValueError('Same selected strict positive singleton s_c required')
        checked_name=PREFIX+'current_inner_exit_strict_collar_check.json';checked=json.loads((HERE/checked_name).read_bytes())
        if not checked.get('all_passed') or checked['actual_five_defect_family_sha256']!=self.family:
            raise ValueError('Accepted same-source positive interior phase-origin receipt required')
        for name in (collar_name,checked_name,PREFIX+'current_generic_loop_function_sources.py',PREFIX+'current_generic_shear_loop_domain.py'):
            fields.previous.bind(self.hashes,name,sha(name))
        self.saved_downstream=json.loads(gzip.decompress((HERE/downstream.NAME).read_bytes()))
        if not self.saved_downstream[downstream.GATE] or self.saved_downstream['source_family']!=self.family:
            raise ValueError('Accepted same current R100/Rm local driver required')
        self.bindings=source_bindings()
        for name in (Path(__file__).name,Path(macro.__file__).name):fields.previous.bind(self.hashes,name,sha(name))
        if require_checked:
            receipt=json.loads((HERE/RECEIPT).read_bytes())
            if not receipt.get('all_passed') or not receipt.get(GATE) or receipt['source_family']!=self.family:
                raise ValueError('Accepted actual original macro finite-N receipt required')
            for name,digest in receipt['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def query(self,label,left,right,N=257):
        N=phase.candidate_N(N);key=(label,left,right,N)
        if key in self.cache:return self.cache[key]
        op,proof=self.moments.owner(label);series=self.downstream.first.owner(label)
        ref=self.downstream.downstream.reference.owner(label);f,c=ref.flow,ref.c
        _,_,inlet=self.downstream.r100.owner(label)
        if op.flow is not f or series.flow is not f:raise ValueError('One current live macro/endpoint source algebra required')
        if f.logs[0]._mpi_!=fields.previous.read_interval(c,self.collar['exact_source_width_log_enclosure'])._mpi_:
            raise ValueError('Exact current bridge/collar source-width tuple required')
        if downstream.canonical_expression(inlet['pressure_axis_over_Pstar_squared_axial5'])!=downstream.canonical_expression(ref.P0):
            raise ValueError('Same original live analytic P0 required')
        with mp.workdps(c.dps+40):
            source=macro_cell(op,series,left,right);proxy,raw,recovered,a=downstream.recover_source(ref,inlet,source)
            recovered.pop('source_frame_conditional_on_same_current_R100_background')
            recovered['source_frame_conditional_on_same_actual_macro_micro_inlet']=True
            roots,qr,proof=downstream.general_quotients(proxy,recovered,a,self.downstream.downstream.parameters.eta_log)
            got=primitives.all_u_primitive_bounds(f,roots,qr,self.downstream.downstream.parameters.dstar_log,c.mpf([0,1]))
            got=downstream.bounded_exponent_cover(f,got,N)
            if qr[C0].zero and qr[Z].zero:
                got['values']={key:f.scalar(0) for key in got['values']}
                got['record']['original_exact_flat_inverse_A_B_and_Z_zero']=True
            E,V=recovered['common_velocity_E_axial5'],recovered['common_velocity_V_axial5']
            density=phase.densities.density_Z_kernels(E[0],E[1],V[0],V[1],got['values'],N)
            result=dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=ref.P0,
                exact_cell=dict(left=left,right=right),actual_macro_background_source=source,original_raw_source=raw,
                original_generic_source=recovered,original_full_source_quotients=proof,
                original_roots=roots['roots'],original_q_C0_Z=qr,original_primitive_values=got['values'],
                original_primitive_proof=got['record'],original_signed_five_density_C0_Z=density,
                **actual_phase_record(f,source,self.sc),
                actual_phase_full_period_cover=c.mpf([0,1]),real_finite_N_micro_exit_boundary_supplied=False,
                **dict.fromkeys(fields.previous.OPEN,False))
        self.cache[key]=result;return result

    def contribution(self,label,N=257):
        N=phase.candidate_N(N)
        if N!=self.saved_downstream['candidate_N']:raise ValueError('Same current N required; no earlier-N transplant')
        series=self.downstream.first.owner(label);f,c=series.flow,series.c
        ref=self.downstream.downstream.reference.owner(label);total={name:[f.scalar(0),f.scalar(0)] for name in RATES};cells=[]
        with mp.workdps(c.dps+40):
            for left,right in zip(PARTITION,PARTITION[1:]):
                source=self.query(label,left,right,N);rows={};weight_rows={}
                for name,rate in RATES.items():
                    mass,decay,tail=weights(series,left,right,rate)
                    parameters.positive_source(f,mass,'actual_macro_cell_positive_Duhamel_mass')
                    density=source['original_signed_five_density_C0_Z']
                    pair=[bounds.symmetric(f,bounds.magnitude(f,density[part][name])*mass*tail)
                          for part in ('kernels','Z_derivatives')]
                    for n,row in enumerate(pair):total[name][n]+=row
                    rows[name]=pair;weight_rows[name]=dict(positive_full_mass=mass,incoming_decay=decay,
                        downstream_suffix_decay=tail,own_rate=str(rate),true_physical_measure_once=True)
                cells.append(dict(source=source,signed_cell_driver_C0_Z=rows,own_rate_weights=weight_rows))
                print('Current original macro source cell',label,left,right,flush=True)
            memory={name:weights(series,(0,1),(1,1),rate)[1] for name,rate in RATES.items()}
            accepted=self.saved_downstream['frames'][label];downstream.guard_saved_source(self.family,label,N,ref.P0,accepted)
            restore=lambda row:downstream.first.endpoint.restore_row(f,row)
            suffix={name:restore(row) for name,row in accepted['retained_R100_Rm_incoming_memory'].items()}
            child={name:[restore(row) for row in rows] for name,rows in accepted['actual_R100_Rm_local_driver_C0_Z'].items()}
            complete={name:[total[name][n]*suffix[name]+child[name][n] for n in range(2)] for name in RATES}
            complete_memory={name:memory[name]*suffix[name] for name in RATES}
        return dict(source_family=self.family,source_frame=label,candidate_N=N,exact_common_P0_axial5=ref.P0,
            actual_macro_source_cells=cells,actual_micro_exit_R100_local_driver_C0_Z=total,
            retained_micro_exit_R100_incoming_memory=memory,actual_micro_exit_Rm_local_driver_C0_Z=complete,
            retained_micro_exit_Rm_incoming_memory=complete_memory,
            original_micro_exit_radius='Ra*exp(2hb)',original_macro_length=f.scalar(f.Y)-f.h*2,
            exact_affine_boundary_formula='deltaH(Rm)=memory(R0,Rm)*deltaH(R0)+all_local_drivers(R0,Rm)',
            genuine_finite_N_micro_exit_correction_still_unsupplied=True,actual_finite_N_R100_correction_supplied=False,
            actual_finite_N_Rm_incoming_correction_supplied=False,accepted_downstream_rehydrated_in_same_live_source_algebra=True,
            **dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalMacroFiniteN(require_checked=False)
    frames={label:owner.contribution(label) for label in ('0','.5')}
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},candidate_N=257,
        original_source_bindings=owner.bindings,frames=downstream.downstream.serialized(frames),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()
