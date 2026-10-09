"""Hydrate-only original Ra microscopic prescribed-shear source functions.

Checked comparison width coefficients determine the known direction. Both
microscopic controls act in both field equations. Actual histories retain
the real core atoms; no ancestor constructor or producer is executed.
"""
import ast
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time
import mpmath as mp
import lei_ren_part1_paper_compliant_current_original_macro_finite_N as downstream

fields,parameters,ep=downstream.fields,downstream.parameters,downstream.ep
bridge=downstream.macro.fields.bridge
comparison=bridge.comparison_module
HERE,PREFIX,sha=downstream.HERE,downstream.PREFIX,downstream.sha
NAME=PREFIX+'current_original_micro_functions.json.gz'
RECEIPT=PREFIX+'current_original_micro_functions_check.json'
GATE='current_original_Ra_micro_hydrated_whole_source_functions_installed'
CHARTS=('first_micro','second_micro')


def serialized(value):
    if isinstance(value,comparison.WidthPolynomial):return serialized(value.packet())
    if isinstance(value,fields.IntervalTaylor):return serialized(list(value.coefficients))
    if isinstance(value,dict):return {key:serialized(row) for key,row in value.items()}
    if isinstance(value,(list,tuple)):return [serialized(row) for row in value]
    return fields.serialized(value)


def source_bindings():
    tree=ast.parse(Path(comparison.__file__).read_text(encoding='utf8'))
    fn=next(n for n in ast.walk(tree) if isinstance(n,ast.FunctionDef) and n.name=='fields')
    expressions=("a=data['L1']*A0", "b=data['L2']*A1",
        "Rg=data['L3_norm']*s**3/6",
        "Rphi=jet_norm(p,self.weight)*(Rg*c.exp(e)+na*nb+self.cap*nb*nb/2+(na+self.cap*nb)**3*c.exp(self.cap*na+self.cap**2*nb)/6)")
    for text in expressions:
        target=ast.parse(text).body[0]
        # a,b,Rg share an original semicolon-separated source line, but each
        # remains an independent assignment in the parsed source tree.
        if not any(ast.dump(n)==ast.dump(target) for n in ast.walk(fn)):
            raise ValueError('Defining comparison source changed: '+text)
    return dict(original_micro_controls=bridge.bridge_control_bindings(),
        recovered_source_coefficients='p, p*L1, p*(L2+L1^2)/2; V0,V1,V2/2 at comparison s=1',
        recovered_L3_bound='L3_norm <= (6/8)*Rphi_upper(s=2)/positive_lower(p[0])',
        recovered_V3_bound='V3_norm <= (6/8)*RV_upper(s=2)',
        original_comparison_remainder_assignments=expressions,
        actual_first_equations='phi_y=-((1-sigma)+hb*sigma)*Dbar*phi/2; V_y=-((1-sigma)+hb*sigma)*(phi/barphi)*drive',
        actual_second_equations='phi_y=-hb*Dbar*phi/2; V_y=-hb*(phi/barphi)*drive',
        actual_geometry='R=Ra*exp(hb*s), dy=hb*ds',
        recovered_norms_are_upper_bounds_not_source_values=True,
        no_later_R100_switch_controls_reused=True)


class HydratedComparison:
    """Only pure source arithmetic, sharing the current original flow."""
    def __init__(self,flow,records,label,pressure):
        self.flow=flow;self.ctx=c=flow.c;self.weight=flow.weight
        entry=records['comparison_point_integrals']['comparison_point_packets'][label]['micro']
        self.points={str(ep(fields.previous.read_interval(c,row['coordinate']))[0]):row for row in entry}
        zero,one,two=(self.points[str(mp.mpf(n))] for n in (0,1,2))
        self.cap=fields.previous.read_interval(c,one['width_remainder_cap'])
        if ep(flow.logs[0])[1]>=ep(c.ln(self.cap))[0]:raise ValueError('Same true hb must be bounded by the comparison error cap')
        self.Z=fields.previous.read_interval(c,zero['Z'])
        for point in entry:
            if (fields.previous.read_interval(c,point['Z'])._mpi_!=self.Z._mpi_
                or fields.previous.read_interval(c,point['axial_jet_weight'])._mpi_!=self.weight._mpi_
                or fields.previous.read_interval(c,point['original_positive_width_log_enclosure'])._mpi_!=flow.logs[0]._mpi_):
                raise ValueError('Same comparison frame, weight and source width required')
        jet=lambda rows:fields.IntervalTaylor(c,[fields.previous.read_interval(c,v) for v in rows])
        decode=lambda packet:[jet(row) for row in packet['signed_hb_power_axial_coefficients']]
        p,V=decode(one['comparison_phi']),decode(one['comparison_raw_V'])
        if any(v.order!=6 for v in p+V):raise ValueError('Original comparison Z0..6 coefficients required')
        lower=ep(p[0][0])[0]
        if lower<=0:raise ValueError('Strictly positive original core phi coefficient required')
        L1=p[1]/p[0];L2=2*p[2]/p[0]-L1*L1
        Rphi=fields.previous.read_interval(c,two['comparison_phi']['omitted_hb_cubed_weighted_axial_jet_norm_upper'])
        RV=fields.previous.read_interval(c,two['comparison_raw_V']['omitted_hb_cubed_weighted_axial_jet_norm_upper'])
        self.data=dict(Z=self.Z,phi0=p[0],V0=V[0],L1=L1,L2=L2,V1=V[1],V2=V[2]*2,
            L3_norm=c.mpf(ep(Rphi*c.mpf(3)/4/c.mpf(lower))[1]),
            V3_norm=c.mpf(ep(RV*c.mpf(3)/4)[1]),
            moments={name:decode(packet)[0] for name,packet in zero['comparison_own_six_moments'].items()})
        amplitude=records['anchored_axis_amplitude']['anchored_amplitude_packets'][label]
        if fields.previous.read_interval(c,amplitude['Z'])._mpi_!=self.Z._mpi_:raise ValueError('Same actual anchored amplitude frame required')
        core=records['core_transfer'];Lambda=fields.previous.read_interval(c,core['Lambda'])
        gradient=[fields.previous.read_interval(c,v) for v in amplitude['G_gradient_taylor_coefficients']]
        def ratios(multiplier):
            values=[c.mpf(1)]
            for n in range(1,7):values.append(sum((-multiplier*Lambda*gradient[j]*values[n-1-j] for j in range(n)),c.mpf(0))/n)
            return [v*math.factorial(n) for n,v in enumerate(values)]
        self.inputs=dict(p0=fields.IntervalTaylor(c,[c.mpf(ep(v)) for v in pressure.normalized_jets(ep(self.Z),6)['normalized_pressure_coefficients']]),
            F0_ratios=ratios(1),F0_squared_ratios=ratios(2))
        self.delta=c.exp(fields.previous.read_interval(c,records['pressure_source']['compliant_source']['parameter_bounds']['log_delta']))
        self.proof=dict(recovered_comparison_phase_one_width_coefficients=True,
            recovered_third_derivative_norms_are_upper_bounds=True,positive_phi_lower=lower,
            source_L3_norm_upper=self.data['L3_norm'],source_V3_norm_upper=self.data['V3_norm'],
            source_remainder_formula_recomputed_on_whole_cell=True,
            endpoint_remainder_not_used_as_whole_cell_remainder=True,
            comparison_source_namespace=zero['comparison_source_namespace'],
            source_width_not_selected_or_materialized=True)
        self.comparison=self;self.cache={}

    def fields(self,data,s,A0,A1):
        return comparison.CompliantComparisonPointIntegrals.fields(self,data,s,A0,A1)

    def cover(self,left,right):
        c=self.ctx;l,r=ep(c.mpf(left))[0],ep(c.mpf(right))[1]
        if not 0<=l<=r<=2:raise ValueError('Ordered original comparison micro cell in[0,2] required')
        key=(l,r)
        if key in self.cache:return self.cache[key]
        cover=bridge.CompliantActualBridgeIntegrals.comparison_micro_cover(self,self.data,l,r)
        phi,V=(bridge.value_enclosure(cover[key]) for key in ('phi','V'))
        own=bridge.named_moments({name:bridge.value_enclosure(row) for name,row in cover['moments'].items()})
        directions=bridge.direction(c,self.Z,self.delta,phi,V,own,self.inputs['p0'],self.inputs['F0_ratios'],self.inputs['F0_squared_ratios'])
        self.cache[key]=dict(comparison=cover,phi=phi,V=V,directions=directions,
            quotient=self.data['phi0'].truncate(5)/phi.truncate(5))
        return self.cache[key]


class MicroFunctions:
    """The two coupled original micro windows on one live current flow."""
    def __init__(self,flow,series,decoder,actual_core):
        self.flow=flow;self.c=c=flow.c;self.series=series;self.decoder=decoder
        if series.flow is not flow or decoder.flow is not flow:raise ValueError('One source basis and ledger required')
        jet=lambda rows:fields.IntervalTaylor(c,[fields.previous.read_interval(c,v) for v in rows])
        self.phi0=flow.jet(jet(actual_core['actual_phi_axial5']))
        self.V0=flow.jet(jet(actual_core['actual_raw_V_axial5']))
        self.initial={name:flow.jet(jet(row)) for name,row in actual_core['actual_own_six_moments_axial5'].items()}
        self.known={key:decoder.cover(*ends) for key,ends in zip(CHARTS,((0,1),(1,2)))}
        self.Ra=flow.factor((0,0,0,1,0));self.Rcover=self.Ra*series.scalar_series(1,c.mpf([0,2]))[0]
        self.parts=fields.PARTS;self.cache={}

    def exp_rows(self,ell):
        f,c=self.flow,self.c;bound=f.norm(ell);upper=f.small_upper(bound)
        if ep(upper)[1]>c.mpf('.5'):raise ValueError('Actual micro full log jet norm must be <=1/2')
        quadratic=f.scale(f.multiply(ell,ell),c.mpf('.5'))
        error=bound*bound*bound*(c.exp(upper)/6)
        exp=f.add([f.scalar(1)]+[f.scalar(0)]*5,ell,quadratic,f.error_rows(error))
        return exp,bound,error

    def prefix(self,chart,s):
        f,c=self.flow,self.c;s=c.mpf(s);lo,hi=ep(s)
        if chart not in CHARTS or (chart=='first_micro' and not 0<=lo<=hi<=1) or (chart=='second_micro' and not 1<=lo<=hi<=2):
            raise ValueError('Original first or second normalized coordinate required')
        if chart=='first_micro':wc,ws=bridge.pulse_control_masses(c,s)
        else:wc,ws=c.mpf('.5'),c.mpf('.5')
        terms=[('first_micro',f.h*wc),('first_micro',f.h*f.h*ws)]
        if chart=='second_micro':terms.append(('second_micro',f.h*f.h*(s-1)))
        angular=[];nominal={part:[] for part in self.parts};proof=[]
        for window,mass in terms:
            known=self.known[window];dirs=known['directions'];q=f.jet(known['quotient'])
            d=f.jet(dirs['D_over_R'].truncate(5))
            angular.append(f.scale(d,-self.Rcover*mass*c.mpf('.5')))
            for part in self.parts:
                power=2 if part=='swirl' else 1
                scale=f.factor((0,int(part=='pressure'),int(part=='swirl'),0,0))
                drive=f.jet(dirs['drive_'+part].truncate(5))
                radius=self.Rcover*self.Rcover if power==2 else self.Rcover
                nominal[part].append(f.scale(f.multiply(q,drive),-radius*mass*scale))
            proof.append(dict(window=window,true_positive_control_mass=mass,
                original_full_window_comparison=known['comparison'],source_measure_once=True))
        ell=f.add(*angular);exp,bound,error=self.exp_rows(ell)
        phi=f.multiply(self.phi0,exp)
        # Every partial prefix is controlled by absolute integral masses,
        # independently of any cancellation in the final signed sum.
        absolute_log=sum((f.norm(row) for row in angular),f.scalar(0))
        log_upper=f.small_upper(absolute_log);delta_exp=absolute_log*c.exp(log_upper)
        deltaV={};errors={}
        for part in self.parts:
            rows=nominal[part];row=f.add(*rows)
            errors[part]=sum((f.norm(item) for item in rows),f.scalar(0))*delta_exp
            deltaV[part]=f.add(row,f.error_rows(errors[part]))
        V=f.add(self.V0,*deltaV.values())
        return dict(phi=phi,V=V,ell=ell,deltaV=deltaV,full_exponential_error=error,
            full_axial_nonlinear_errors=errors,absolute_full_prefix_log_norm=absolute_log,source_control_terms=proof)

    def histories(self,chart,s):
        f,c=self.flow,self.c;s=c.mpf(s);lo,hi=ep(s)
        if lo==hi==0:return self.initial,dict(exact_original_core_atoms=True)
        full=self.prefix('first_micro',c.mpf([0,hi])) if hi<=1 else self.prefix('second_micro',c.mpf([1,hi]))
        if hi>1:
            first_full=self.prefix('first_micro',c.mpf([0,1]))
            # Covers both first and second actual prefixes before applying
            # positive Volterra weights; no endpoint replacement is made.
            def hull(a,b):
                difference=b-a
                return a+downstream.bounds.symmetric(f,downstream.bounds.magnitude(f,difference))
            full={name:[hull(a,b) for a,b in zip(first_full[name],full[name])] for name in ('phi','V')}
        phi,V=full['phi'],full['V']
        rhs=dict(H=f.scale(phi,2),M=V,K=f.scale(f.multiply(phi,V),2),A=f.multiply(V,V),B=f.multiply(phi,phi),C=f.multiply(phi,phi))
        result={};weights={}
        for name,rate in bridge.RATES.items():
            decay=self.series.scalar_series(-rate,s)[0]
            average,error=self.series.scalar_series(-rate,s,average=True)
            mass=f.h*s*average
            result[name]=f.add(f.scale(self.initial[name],decay),f.scale(rhs[name],mass))
            weights[name]=dict(true_incoming_decay=decay,true_positive_Volterra_mass=mass,complete_mass_error=error,
                original_core_memory_retained=True,source_rhs_covers_every_actual_prefix=True)
        return result,weights

    def evaluate(self,chart,left,right):
        f,c=self.flow,self.c;l,r=c.mpf(left),c.mpf(right)
        if ep(l)[0]!=ep(l)[1] or ep(r)[0]!=ep(r)[1]:raise ValueError('Exact source cell endpoints required')
        low,high=ep(l)[0],ep(r)[1]
        if (chart not in CHARTS or low>high or
            (chart=='first_micro' and not 0<=low<=high<=1) or
            (chart=='second_micro' and not 1<=low<=high<=2)):
            raise ValueError('Ordered source cell in the defining micro chart required')
        s=c.mpf([ep(l)[0],ep(r)[1]]);key=(chart,s._mpi_)
        if key in self.cache:return self.cache[key]
        at=self.prefix(chart,s);own,weights=self.histories(chart,s)
        local=self.decoder.cover(ep(l)[0],ep(r)[1]);dirs=local['directions']
        R=self.Ra*self.series.scalar_series(1,s)[0]
        rootR=f.factor((0,0,0,.5,0))*self.series.scalar_series(c.mpf('.5'),s)[0]
        D=f.scale(f.jet(dirs['D_over_R'].truncate(5)),R)
        positive=parameters.positive_source(f,D[0],'actual_micro_comparison_Dbar')
        if chart=='first_micro':
            sig=downstream.long.intersection(c,bridge.sigma_jets(c,s)[0],c.mpf([0,1]))
            chi=f.scalar(1-sig)+f.h*sig
        else:sig=c.mpf(0);chi=f.h
        a=f.scale(D,chi)
        if ep(f.logs[0])[1]>=0:raise ValueError('Original micro hb<1 required for correlated control bound')
        lower=positive['source_log_lower']+f.logs[0]
        if ep(a[0].coefficient)[0]<=0:a[0]=a[0].positive_intersection(lower)
        barphi=local['phi'].truncate(5)
        quotient=f.multiply(at['phi'],f.jet(barphi.reciprocal()))
        drives={part:f.scale(f.jet(dirs['drive_'+part].truncate(5)),R*R if part=='swirl' else R) for part in self.parts}
        Vy=f.add(*[f.scale(f.multiply(quotient,row),-chi*f.factor((0,int(part=='pressure'),int(part=='swirl'),0,0))) for part,row in drives.items()])
        result=dict(chart=chart,phase=s,fields=dict(phi=at['phi'],V=at['V']),histories=own,
            phi_y=f.scale(f.multiply(a,at['phi']),-c.mpf('.5')),V_y=Vy,correlated_a_axial5=a,
            actual_Dbar_source=D,actual_positive_Dbar_proof=positive,original_sigma=sig,original_coupled_chi=chi,
            radius=R,root_radius=rootR,window_length=f.h,physical_log_radius_over100=f.scalar(-f.Y)+f.h*s,
            exact_log_radius_over_Ra=f.h*s,original_prefix_field_evidence=at,original_core_Volterra_history_evidence=weights,
            original_comparison_source_evidence=local['comparison'],original_current_phi_over_barphi=quotient,
            original_current_drive_sources=drives,actual_core_background_inlet_preserved=True,
            true_physical_measure='dy=hb*ds',source_interval_functions_not_endpoint_hulls=True,
            real_finite_N_micro_exit_correction_supplied=False)
        self.cache[key]=result;return result


class OriginalMicroFunctions:
    mode='current_original_Ra_micro_whole_source_from_checked_width_coefficients'
    def __init__(self,dps=500,require_checked=True):
        self.downstream=downstream.OriginalMacroFiniteN(dps);self.c=self.downstream.c
        self.family=self.downstream.family;self.hashes=dict(self.downstream.hashes);self.cache={}
        self.bindings=source_bindings();self.origin=self.downstream.moments.fields
        for module in (comparison,bridge,Path(__file__)):
            name=module.name if isinstance(module,Path) else Path(module.__file__).name
            fields.previous.bind(self.hashes,name,sha(name))
        if require_checked:
            checked=json.loads((HERE/RECEIPT).read_bytes())
            if not checked.get('all_passed') or checked['source_family']!=self.family:raise ValueError('Accepted same-source micro function receipt required')
            for name,digest in checked['input_hashes'].items():fields.previous.bind(self.hashes,name,digest)
            fields.previous.bind(self.hashes,RECEIPT,sha(RECEIPT))

    def owner(self,label):
        if label not in ('0','.5'):raise ValueError('Current admitted axial frames0,.5 required')
        if label in self.cache:return self.cache[label]
        op,_=self.downstream.moments.owner(label);f=op.flow;series=self.downstream.downstream.first.owner(label)
        decoder=HydratedComparison(f,self.origin.records,label,self.origin.pressure)
        actual=self.origin.records['actual_bridge_integrals']['packets'][label]['core']
        self.cache[label]=MicroFunctions(f,series,decoder,actual);return self.cache[label]

    def query(self,label,chart,left,right):
        owner=self.owner(label)
        with mp.workdps(self.c.dps+40):value=owner.evaluate(chart,left,right)
        return dict(source_family=self.family,source_frame=label,original_micro_source=value,
            hydrate_only_comparison_decoder_proof=owner.decoder.proof,actual_source_coordinate_Z=owner.decoder.Z,
            actual_selected_positive_s_c=self.downstream.sc,
            true_phase_origin_log_offset=owner.flow.h*(self.downstream.sc/2),
            no_ancestor_constructors_or_producers_executed=True,**dict.fromkeys(fields.previous.OPEN,False))


def run():
    began=time.monotonic();owner=OriginalMicroFunctions(require_checked=False)
    frames={}
    for label in ('0','.5'):
        frames[label]=[]
        for chart,start in (('first_micro',0),('second_micro',1)):
            for l,r in ((0,.25),(.25,.5),(.5,.75),(.75,1)):
                frames[label].append(owner.query(label,chart,start+l,start+r))
                print('Current original micro source',label,chart,l,r,flush=True)
    report=dict(mode=owner.mode,source_family=owner.family,**{GATE:True},original_source_bindings=owner.bindings,
        frames=serialized(frames),input_hashes=owner.hashes,
        execution_seconds=time.monotonic()-began,**dict.fromkeys(fields.previous.OPEN,False))
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(downstream.base.encoded(report),separators=(',',':'))+'\n').encode(),compresslevel=6,mtime=0))
    return report


if __name__=='__main__':run()
