"""Genuine core-inlet C6 data and full signed generic Q/inertial/p2 C5.

The inlet is the true rho=4 core endpoint. Every moment retains the admitted
nonlinear solution and its complete integral history. This does not install
an order-six source on the seventeen downstream active charts.
"""
import copy
from dataclasses import dataclass
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace
import time

import mpmath as mp
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
import lei_ren_part1_paper_compliant_current_original_n1_analytic_core_lower_pressure_C3_functions as previous
import lei_ren_part1_paper_compliant_current_core_interior_moments as interior
import lei_ren_part1_paper_compliant_current_original_bridge_macro_functions as macro
import lei_ren_part1_paper_compliant_current_original_Rm_generic_inputs as generic

source,parent,controls=previous.source,previous.parent,previous.controls
HERE,PREFIX,sha=previous.HERE,previous.PREFIX,previous.sha
NAME=PREFIX+'current_original_core_inlet_axial6_generic_source.json.gz'
RECEIPT=PREFIX+'current_original_core_inlet_axial6_generic_source_check.json'
GATES=('current_original_true_core_inlet_C6_generic_source_functions_installed',
    'current_original_true_core_inlet_Q_inertial_p2_C5_functions_installed',
    'current_original_true_core_inlet_axial6_numeric_source_installed')
OPEN=tuple(dict.fromkeys(previous.OPEN+('current_original_all_17_chart_axial6_generic_source_installed',
    'current_original_compact_C5_target_and_control_functions_installed')))


@dataclass(frozen=True)
class FunctionJet:
    rows:tuple
    @property
    def value(self):return self.rows[0]
    @property
    def order(self):return len(self.rows)-1


def rows(jet):return jet.rows if type(jet) is FunctionJet else previous.higher.rows(jet)


def encode(value):
    if type(value) is FunctionJet:return [q.node for q in value.rows]
    if isinstance(value,dict):return {k:encode(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [encode(v) for v in value]
    return previous.higher.encode(value)


class JetAlgebra:
    """Ordinary source-function derivatives, with exact Leibniz factors."""
    def __init__(self,g,order):self.g,self.order=g,order
    def fixed(self,value):
        ref=value if type(value) is source.FunctionRef else self.g.constant(value)
        self.g.ids([ref]);return FunctionJet((ref,)+(self.g.zero,)*self.order)
    def trim(self,jet):
        if len(rows(jet))<self.order+1:raise ValueError('Required genuine derivative source rows absent')
        return FunctionJet(tuple(rows(jet)[:self.order+1]))
    def add(self,*jets):return FunctionJet(tuple(self.g.add(*(rows(j)[k] for j in jets)) for k in range(self.order+1)))
    def neg(self,jet):return FunctionJet(tuple(self.g.neg(q) for q in rows(self.trim(jet))))
    def sub(self,left,right):return self.add(left,self.neg(right))
    def mul(self,left,right):
        l,r=rows(self.trim(left)),rows(self.trim(right));g=self.g
        return FunctionJet(tuple(g.add(*(g.mul(g.constant(math.comb(k,j)),l[j],r[k-j])
            for j in range(k+1))) for k in range(self.order+1)))
    def scale(self,jet,value):return self.mul(jet,self.fixed(value))
    def div(self,left,right,certificate):
        l,r=rows(self.trim(left)),rows(self.trim(right));out=[];g=self.g
        for k in range(self.order+1):
            rest=g.add(*(g.mul(g.constant(math.comb(k,j)),r[j],out[k-j]) for j in range(1,k+1)))
            out.append(g.quotient(g.sub(l[k],rest),r[0],certificate))
        return FunctionJet(tuple(out))
    def shift(self,jet):
        if len(rows(jet))<self.order+2:raise ValueError('Axial differentiation requires the genuine next source row')
        return FunctionJet(tuple(rows(jet)[1:self.order+2]))
    def at(self,jet,variable,value):
        return FunctionJet(tuple(self.g.node('function_substitution',expression=q.node,
            variable=variable.node,value=value.node,Z_independent_substitution=True) for q in rows(self.trim(jet))))
    def integral(self,jet,variable,lower,upper,**attrs):
        g=self.g
        return FunctionJet(tuple(controls.integral(g,q,g.nodes[variable.node]['name'],lower,upper,
            endpoint_independent_of_Z=True,ordinary_slow_Z_derivative_order=k,
            true_core_inlet_complete_history=True,numerical_cap_not_defining_integrand=True,**attrs)
            for k,q in enumerate(rows(self.trim(jet)))))


def extend6(derivative,jet):
    out=list(rows(jet))
    while len(out)<7:out.append(derivative.derivative(out[-1]))
    if len(out)!=7:raise ValueError('A genuine order-six source is required')
    return FunctionJet(tuple(out))


def generic_functions(g,fields,delta,Pstar,radius):
    """The unchanged full generic equations, now with genuine sixth inputs."""
    a=JetAlgebra(g,5);z=FunctionJet((g.symbol('Z'),g.one)+(g.zero,)*4)
    one=a.fixed(1);de=a.fixed(delta);d=a.sub(one,a.mul(z,z));L=a.sub(one,a.mul(de,a.mul(z,z)))
    E,V,m,h,k,e,p,P0=(a.trim(fields[name]) for name in ('E','V','m','h','k','e','p','P0'))
    mz,hz,kz,ez,pz=(a.shift(fields[name]) for name in ('m','h','k','e','pressure'))
    pressure=a.trim(fields['pressure']);half=lambda j:a.scale(j,'1/2')
    transport=a.add(a.mul(a.sub(one,de),a.mul(z,m)),a.mul(d,mz))
    div=lambda j:a.div(j,L,'same positive L=1-delta*Z^2 on original real axis')
    Q=div(a.sub(a.scale(a.mul(z,V),2),transport))
    theta_linear=div(a.sub(a.add(a.neg(E),a.mul(a.sub(one,half(de)),h)),
        half(a.mul(a.sub(one,de),a.mul(z,hz)))))
    theta_quadratic=div(a.add(a.mul(a.sub(a.scale(de,2),one),a.mul(z,k)),a.neg(a.mul(d,kz)),a.mul(E,transport)))
    axial_linear=div(a.add(a.neg(V),half(a.mul(a.sub(one,de),a.sub(m,a.mul(z,mz))))))
    axial_quadratic=div(a.add(a.mul(V,transport),a.scale(a.mul(de,a.mul(z,e)),2),a.neg(a.mul(d,ez)),
        a.scale(a.mul(a.add(one,de),a.mul(z,pressure)),2),a.neg(a.mul(d,pz))))
    p2=a.scale(a.div(a.add(axial_linear,a.scale(axial_quadratic,Pstar)),E,
        'same actual positive core inlet swirl E'),radius)
    return dict(Q=Q,transport=transport,theta_linear=theta_linear,theta_quadratic=theta_quadratic,
        axial_linear=axial_linear,axial_quadratic=axial_quadratic,p2=p2)


class CoreSixMoments(interior.CurrentCoreInteriorMoments):
    """Reuse exact core moment/tail production without an unused tensor owner."""
    def __init__(self,common):
        if type(common) is not interior.CurrentCoreCommonFixedPoint or not common.acceptance_loaded:
            raise TypeError('The same checked actual nonlinear common core is required')
        self.common=common;self.core=common.core;self.rebuild=common.rebuild;self.ctx=common.ctx
        self.family,self.source,self.datum_sha=common.family,common.source_sha,common.datum_sha
        self._owners=(common,self.core,self.rebuild,self.ctx);self.proof=interior.source_theorem()
        self.cache={};self.rows_cache={};self.tail_cache={};self.assert_graph()
    def assert_graph(self):
        s=self.common.source
        if (self._owners!=(self.common,self.core,self.rebuild,self.ctx) or not self.common.acceptance_loaded
            or not s.acceptance_loaded or self.rebuild is not s.rebuild or self.core is not s.core
            or self.ctx is not self.core.ctx or not all(s.current_source_graph().values())):
            raise ValueError('Same admitted nonlinear core/moment production graph required')


class CurrentOriginalCoreInletAxial6GenericSource:
    def __init__(self,core_lower_owner,hierarchy_owner,require_checked=True):
        if type(core_lower_owner) is not previous.CurrentOriginalN1AnalyticCoreLowerPressureC3Functions or not core_lower_owner.acceptance_loaded:
            raise TypeError('Checked actual analytic core/lower C5 owner required')
        core_lower_owner.assert_graph();hierarchy_owner.assert_graph()
        self.previous=core_lower_owner;self.core=core_lower_owner.core;self.hierarchy=hierarchy_owner
        self.family=copy.deepcopy(core_lower_owner.family);self.g=g=core_lower_owner.g
        if hierarchy_owner.family!=self.family:raise ValueError('Same actual core/pressure source family required')
        self.hashes=dict(core_lower_owner.hashes)
        for module in (previous,interior):
            receipt=parent.linear.read(module.RECEIPT)
            if not receipt['all_passed']:raise ValueError('Checked true source dependency required')
            for name,digest in receipt['input_hashes'].items():self._pin(name,digest)
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self._pin(name,sha(name))
        for name,digest in hierarchy_owner.hashes.items():self._pin(name,digest)
        for module in (macro,generic):self._pin(Path(module.__file__).name,sha(Path(module.__file__).name))
        self.common=self.core._oracle
        if self.common is None:raise ValueError('Retained checked nonlinear core oracle required')
        self.parameter_identity=self.core.parameter_identity(self.common)
        self.moments=CoreSixMoments(self.common);self.numeric_cache={}
        self.derivative=derivative=previous.HigherCoreZDerivative(core_lower_owner.lower)
        for jet in core_lower_owner.core_functions.values():derivative.seed_rows(jet)
        for jets in core_lower_owner.core_radial_rows.values():
            for jet in jets:derivative.seed_rows(jet)
        self.core_functions={name:extend6(derivative,jet) for name,jet in core_lower_owner.core_functions.items()}
        self.P0=self.P0_functions();a=JetAlgebra(g,6);c=g.constant
        self.Phi=self.core_functions['Phi'];self.F0=self.core_functions['F0'];rho=self.core.rho
        self.U=a.add(self.core_functions['U0'],a.scale(self.core_functions['Psi'],self.core.parameters['epsilon']))
        t=g.symbol('original_core_inlet_moment_t');sample=g.mul(rho,t)
        phi=a.at(self.Phi,rho,sample);u=a.at(self.U,rho,sample)
        densities=dict(H=a.scale(phi,g.mul(c(2),t)),M=u,K=a.scale(a.mul(phi,u),g.mul(c(2),t)),
            A=a.mul(u,u),B=a.scale(a.mul(phi,phi),t),C=a.mul(phi,phi))
        self.moment_functions={key:a.integral(body,t,g.zero,g.one,actual_core_normalized_moment=key,
            same_unique_nonlinear_core_and_full_density_product=True) for key,body in densities.items()}
        self.inlet_moments={key:a.at(jet,rho,c(4)) for key,jet in self.moment_functions.items()}
        self.Pstar=g.unary('exp',self.core.parameters['logP']);self.radius=g.mul(c(4),self.core.parameters['epsilon'])
        F=a.at(self.F0,rho,c(4));P=a.at(self.Phi,rho,c(4));U=a.at(self.U,rho,c(4))
        invS=g.unary('exp',g.neg(self.core.parameters['logP']));invS2=g.unary('exp',g.mul(c(-2),self.core.parameters['logP']))
        common_swirl=a.scale(F,g.mul(g.unary('sqrt',g.mul(c(2),self.radius)),invS))
        FF=a.scale(a.mul(F,F),g.mul(self.radius,invS2));J=self.inlet_moments
        self.fields=dict(E=a.mul(common_swirl,P),V=a.scale(U,invS),m=a.scale(J['M'],invS),
            h=a.scale(a.mul(common_swirl,J['H']),'1/2'),k=a.scale(a.mul(common_swirl,J['K']),g.mul(c('1/2'),invS)),
            e=a.sub(a.scale(J['A'],invS2),a.mul(FF,J['B'])),p=a.mul(FF,J['C']),P0=self.P0)
        self.fields['pressure']=a.add(self.P0,self.fields['p'])
        self.functions=generic_functions(g,self.fields,self.core.parameters['delta'],self.Pstar,self.radius)
        self.definition=dict(radial_scope='true core inlet rho=4; R=4/Lambda',
            genuine_seven_source_rows_before_axial_differentiation=True,
            same_actual_P0_order6_provider_and_independent_pressure=True,
            full_signed_linear_quadratic_meridional_and_pressure_sectors_retained=True,
            generic_V_mean_and_k_invPstar_applied_once=True,
            p2_quadratic_weight='Pstar, with no extra 1/2',one_positive_physical_radius_factor=True,
            no_Qpad_row_or_enclosing_center_used=True,
            source_function_bodies_separate_from_numeric_enclosures=True,
            all_17_chart_axial6_source_installed=False,compact_C5_targets_controls_installed=False,
            complete_five_defect_C3_and_global_n1=False,
            shared_graph_unreferenced_suffix_allowed_after_exact_pinned_prefix=True)
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self._owners=(self.previous,self.core,self.hierarchy,self.common,self.moments,self.g)
        self.acceptance_loaded=False
        if require_checked:
            raw,receipt=parent.linear.read(NAME),parent.linear.read(RECEIPT)
            pins=dict(receipt.get('input_hashes',{}));archive=pins.pop(NAME,None)
            checker_name=Path(__file__).stem+'_check.py';checker=pins.pop(checker_name,None)
            if (not receipt.get('all_passed') or not all(receipt.get(k) for k in GATES)
                or any(receipt.get(k) for k in OPEN) or receipt['source_family']!=self.family
                or pins!=self.hashes or raw['input_hashes']!=self.hashes or archive!=sha(NAME)
                or checker!=sha(checker_name) or raw['core_inlet_axial6_bindings']!=self.snapshot
                or raw['exact_graph_nodes']!=g.nodes[:len(raw['exact_graph_nodes'])]):
                raise ValueError('Checked actual core-inlet axial6/generic5 source required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed inlet high-order source '+name)
            self.acceptance_loaded=True

    def P0_functions(self):
        old=self.core.pulse.own.jets.native['P0'];g=self.g
        row=g.nodes[old.value.node]
        if not row['coordinate_independent'] or not row['exact_source_op_P0_function_not_enclosure_endpoint']:
            raise ValueError('Same original independent normalized P0 required')
        result=list(previous.higher.rows(old))
        provider=self.hierarchy.source.admission.pressure;datum=self.common.core.datum
        if (provider.source_sha!=datum.source_sha or provider.datum_sha!=datum.datum_sha
            or provider.input_hashes!=datum.input_hashes
            or provider.normalized_jets.__func__ is not datum.normalized_jets.__func__):
            raise ValueError('Identical original fourteen-atom pressure definition required')
        binding=parent.active.current.ast_binding(self.project_P0_function)
        for k in (4,5,6):
            result.append(g.node('current_original_axis_pressure_function_derivative',
                source_family=self.family,original_function=old.value.node,
                Z_variable=g.symbol('Z').node,Z_order=k,ordinary_slow_Z_derivative_order=k,
                Taylor_coefficient_factorial=math.factorial(k),source_projection_binding=binding,
                pressure_provider_binding=parent.active.current.ast_binding(self.P0_rows),
                defining_equation='d_Z^Z_order original_function(Z)',
                source_projection=f'factorial({k})*P0_rows(Z).coefficients[{k}]',
                implicit_source_sha256=provider.source_sha,datum_enclosure_sha256=provider.datum_sha,
                coordinate_independent=True,defining_quantity_not_a_range_value=True,
                true_order6_pressure_provider=True,original_value_handle_index=old.value.node,
                same_actual_fourteen_atoms_and_analytic_flatten_function=True))
        return FunctionJet(tuple(result))

    def P0_rows(self,Z):
        self.hierarchy.assert_graph()
        data=self.hierarchy.source.inputs(Z)['p0']
        if data.order!=6 or len(data.coefficients)!=7:raise ValueError('Actual seven-row pressure provider required')
        return data

    def project_P0_function(self,ref,Z):
        """Enclose the declared ordinary derivative of the original P0 function."""
        self.assert_graph();self.g.ids([ref]);row=self.g.nodes[ref.node]
        if ref not in self.P0.rows[4:] or row['operation']!='current_original_axis_pressure_function_derivative':
            raise ValueError('Exact registered higher pressure derivative required')
        order=row['ordinary_slow_Z_derivative_order']
        if row['original_function']!=self.P0.value.node or row['Taylor_coefficient_factorial']!=math.factorial(order):
            raise ValueError('Original P0 function and ordinary derivative factor required')
        return self.P0_rows(Z).coefficients[order]*math.factorial(order)

    @interior.source_precision
    def numerical(self,ends):
        if not isinstance(ends,(tuple,list)) or len(ends)!=2:raise ValueError('Exact ordered Z endpoints required')
        self.assert_graph();c=self.common.ctx;z=self.hierarchy.source.cell(ends);z=c.mpf(macro.ep(z))
        key=z._mpi_
        if key in self.numeric_cache:return self.numeric_cache[key]
        packet=self.moments.evaluate(z,c.mpf(4));raw,density=self.moments.rows(z)
        amplitude=self.common.source.amplitude.evaluate(z)
        logRa=c.ln(4)-self.common.core.logLambda
        f=macro.MacroFlow(c,c.mpf(macro.ep(self.hierarchy.source.logh)),2*self.common.core.logP,
            2*amplitude['logF0'],logRa,c.ln(100)-logRa,c.mpf(macro.ep(self.hierarchy.source.weight)))
        lift=lambda values:[f.scalar(c.mpf(macro.ep(value))) for value in values]
        mul=lambda l,r:[sum((l[j]*r[k-j] for j in range(k+1)),f.scalar(0)) for k in range(min(len(l),len(r)))]
        scale=lambda values,factor:[value*factor for value in values]
        def profile(name,i=0):
            out=[]
            for k in range(7):
                tail=self.common.rebuild.profile(raw,c.mpf(4),radial_order=i,axial_order=k)['infinite_radial_tail_bounds'][name]
                finite=sum((row[k]*(math.factorial(n)//math.factorial(n-i))*c.mpf(4)**(n-i)
                    for n,row in enumerate(raw['rows']['A' if name=='Phi' else 'Uz']) if n>=i),c.mpf(0))
                out.append(finite+interior.symmetric(c,tail)/math.factorial(k))
            if name=='Phi' and i==0:
                out[0]=interior.intersection(c,out[0],c.mpf([macro.ep(self.common.core.phi_floor)[0],macro.ep(self.common.core.phi_ceiling)[1]]))
            return lift(out)
        phi,U=profile('Phi'),profile('Uz');phir,Ur=profile('Phi',1),profile('Uz',1)
        phi_y=[a*c.mpf('1/2')+4*b for a,b in zip(phi,phir)];U_y=scale(Ur,4)
        J={name:lift([packet['ordinary_rho4_axial6_moment_grids'][name][interior.key(0,k)]/math.factorial(k)
            for k in range(7)]) for name in interior.SPEC}
        F=lift([value/math.factorial(k) for k,value in enumerate(packet['original_F0_relative_axial6_ordinary_derivatives'])])
        F2=lift([value/math.factorial(k) for k,value in enumerate(packet['original_F0_squared_relative_axial6_ordinary_derivatives'])])
        swirl=f.factor((0,-.5,.5,.5,0))*c.sqrt(2);invS=f.factor((0,-.5,0,0,0))
        RFF=f.factor((0,-1,1,1,0));invS2=f.factor((0,-1,0,0,0))
        E=scale(mul(F,phi),swirl);Ey=scale(mul(F,phi_y),swirl)
        M=J['M'];h=scale(mul(F,J['H']),swirl*c.mpf('1/2'));raw_k=scale(mul(F,J['K']),swirl*c.mpf('1/2'))
        energy=[a*invS2-b*RFF for a,b in zip(J['A'],mul(F2,J['B']))]
        p=scale(mul(F2,J['C']),RFF)
        P0=lift(packet['original_analytic_axis_pressure_axial6_Taylor_coefficients'])
        if len(P0)!=7:raise ValueError('Actual independent pressure sixth row required')
        op=SimpleNamespace(flow=f,c=c,zrows=f.jet(IntervalTaylor.variable(c,z,5)),
            reference=SimpleNamespace(z=IntervalTaylor.variable(c,z,5),delta=self.common.core.delta),
            Pstar=f.factor((0,.5,0,0,0)),P0=P0,Rm_factor=f.factor((0,0,0,1,0)))
        actual=dict(raw_current_radius_y_derivative_axial_coefficients=dict(
            velocity=dict(theta=[E,Ey],axial=[U,U_y]),
            histories={name:[values] for name,values in dict(m=M,h=h,k=raw_k,e=energy,p=p).items()}),
            original_P0_normalized_axial5=P0,geometry=dict(point=True,exact_x=[1,1]))
        recovered=generic.recover_inputs(op,actual)
        sectors=recovered['full_signed_inertial_sectors_axial4'];axial=[l+q*op.Pstar
            for l,q in zip(sectors['axial_linear'],sectors['axial_quadratic'])]
        denominator=recovered['common_velocity_E_axial5'];lower=macro.ep(denominator[0].coefficient)[0]
        if lower<=0:raise ValueError('Positive actual core swirl coefficient required')
        loglower=denominator[0].scale.evaluate()+c.ln(c.mpf(lower));quotient=[]
        for n in range(6):
            rest=sum((denominator[j]*quotient[n-j] for j in range(1,n+1)),f.scalar(0))
            quotient.append((axial[n]-rest).positive_divide(denominator[0],loglower))
        # Keep the radius out of the quotient recurrence and apply it once.
        p2=[value*op.Rm_factor for value in quotient]
        sources=dict(E=E,V=scale(U,invS),m=scale(M,invS),h=h,
            k=scale(raw_k,invS),e=energy,p=p,P0=P0)
        if not all(len(values)==7 for values in sources.values()) or any(len(values)!=6 for values in sectors.values()):
            raise ValueError('Sixth input/fifth output source rows were truncated')
        result=dict(source_family=self.family,exact_Z_cell=list(ends),core_rho=4,flow=f,
            actual_complete_seven_Taylor_sources=sources,
            seven_source_units='post-recovery normalized E,V,m,h,k,e,p,P0; invPstar applied once',
            actual_full_signed_inertial_sectors_axial5=sectors,
            actual_radial_Q_axial5=recovered['common_radial_Q_axial4'],actual_full_signed_p2_axial5=p2,
            radial_Q_actual_six_rows_from_genuine_seventh_mean_row_no_Qpad=True,
            actual_core_moment_and_density_tails=packet,
            existing_generic_production_assignments_executed=True,all_seven_rows_genuine_before_differentiation=True,
            source_ranges_not_defining_function_values=True,all_17_chart_source_installed=False)
        self.numeric_cache[key]=result;return result

    def _pin(self,name,digest):
        if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Changed inlet source input '+name)
        self.hashes[name]=digest
    def _snapshot(self):return encode(dict(core_functions=self.core_functions,P0=self.P0,
        moments=self.moment_functions,inlet_moments=self.inlet_moments,fields=self.fields,
        functions=self.functions,Pstar=self.Pstar,radius=self.radius,definition=self.definition))
    def assert_graph(self):
        self.previous.assert_graph();self.hierarchy.assert_graph();self.moments.assert_graph()
        if (self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes
            or (self.previous,self.core,self.hierarchy,self.common,self.moments,self.g)!=self._owners):
            raise ValueError('Actual inlet C6/generic5 source changed')


def run(core_lower_owner,hierarchy_owner):
    began=time.monotonic();owner=CurrentOriginalCoreInletAxial6GenericSource(core_lower_owner,hierarchy_owner,False)
    data=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,core_inlet_axial6_bindings=owner.snapshot,
        input_hashes=owner.hashes,**dict.fromkeys(GATES+OPEN,False),execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(data,separators=(',',':'))+'\n').encode(),mtime=0))
    print('TRUE_CORE_INLET_AXIAL6_GENERIC5_SOURCE_DEFINED',len(owner.g.nodes),flush=True)
    return owner,data


if __name__=='__main__':raise SystemExit('run(core_lower_owner,hierarchy_owner) reuses the live source owners')
