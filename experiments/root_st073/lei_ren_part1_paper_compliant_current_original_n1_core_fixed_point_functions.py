"""Actual admitted core functions, twenty-term map and enclosing oracle.

Component nodes denote the unique analytic solution of the admitted actual
nonlinear equations. A lazy same-family oracle returns finite nonlinear
coefficient enclosures with full tails; no midpoint or model value is selected.
"""
import copy
import ast
import gzip
import hashlib
import inspect
import json
import math
from pathlib import Path
import time
import scipy.integrate

import lei_ren_part1_paper_compliant_current_original_n1_pulse_pressure_functions as pulse
import lei_ren_part1_paper_compliant_current_core_common_fixed_point as common
import lei_ren_part1_paper_compliant_current_core_nonlinear_operator as operator
import lei_ren_part1_paper_compliant_current_original_core_hierarchy_source as hierarchy

parent,source,controls=pulse.parent,pulse.source,pulse.controls
c3=pulse.own.c3
HERE,PREFIX,sha=pulse.HERE,pulse.PREFIX,pulse.sha
NAME=PREFIX+'current_original_n1_core_fixed_point_functions.json.gz'
RECEIPT=PREFIX+'current_original_n1_core_fixed_point_functions_check.json'
GATE='current_original_n1_true_core_fixed_point_functions_and_enclosure_bridge_installed'
OPEN=pulse.OPEN


def map_terms(g,values):
    """Compile all original ten angular and ten axial nonlinear terms."""
    original=operator.production_expressions()
    return {key:[source.expression(g,expr,values) for expr in original[label]]
        for key,label in (('theta','theta_terms'),('z','z_terms'))}


def component_polynomial(packet,rho,epsilon,quantity,radial_order,axial_order,convert):
    """Exact finite ordinary rows; strip the correlated axis before scaling."""
    i,k=radial_order,axial_order
    if quantity not in ('Phi','Psi','mean_Psi'):raise ValueError('Unknown actual core component')
    rows=packet['rows']['A' if quantity=='Phi' else 'Uz'];out=convert(0)
    # Uz row zero is exactly U0. Its subtraction must be structural: two
    # equal interval boxes subtracted numerically would lose correlation.
    lower=i if quantity=='Phi' else max(1,i)
    for n in range(lower,packet['radial_degree']+1):
        coefficient=rows[n][k]
        if quantity!='Phi':coefficient=coefficient/epsilon
        if quantity=='mean_Psi':coefficient=coefficient/(n+1)
        out+=coefficient*(math.factorial(n)//math.factorial(n-i))*math.factorial(k)*rho**(n-i)
    return out


class CurrentOriginalN1CoreFixedPointFunctions:
    def __init__(self,pulse_owner=None,oracle=None,require_checked=True):
        self.pulse=pulse_owner or pulse.CurrentOriginalN1PulsePressureFunctions()
        if not self.pulse.acceptance_loaded:raise ValueError('Checked complete original pulse functions required')
        self.pulse.assert_graph();self.g=g=self.pulse.g;self.family=copy.deepcopy(self.pulse.family)
        self.hashes=dict(self.pulse.hashes)
        self._oracle=None;self._packets={}
        for name,digest in parent.linear.read(pulse.RECEIPT)['input_hashes'].items():self._pin(name,digest)
        for name in (pulse.NAME,pulse.RECEIPT):self._pin(name,sha(name))
        for module,gate in ((common,common.BANACH_GATE),(operator,operator.NORM_BOUNDS_GATE),
                            (hierarchy,hierarchy.GATE)):
            receipt=parent.linear.read(module.RECEIPT)
            if not receipt.get('all_passed') or not receipt.get(gate):raise ValueError('Checked true core dependency required')
            identity=receipt.get('source_family',receipt)
            if any(identity.get(key)!=value for key,value in self.family.items()):raise ValueError('Core source family differs')
            for name,digest in receipt['input_hashes'].items():self._pin(name,digest)
            for name in (module.NAME,module.RECEIPT,Path(module.__file__).name):self._pin(name,sha(name))
        self.rho=g.symbol('original_n1_true_core_rho');z=g.symbol('Z');c=g.constant;a=c3.C3Algebra(g)
        ast=parent.active.current.ast_binding
        self.parameters={}
        for key,path in (('j','common.core.j'),('logC','common.core.logC'),('anchor','common.source.amplitude.roots')):
            self.parameters[key]=g.node('current_original_source_parameter',name='true_core_'+key,
                Z_independent=True,source_family=self.family,exact_fixed_source_parameter=True,
                source_path=path,source_provider=ast(self.attach_oracle),
                parameter_enclosure_not_selected_value=True,
                anchor_definition='unique real H0 root in [-j,0]' if key=='anchor' else None)
        Gamma=self.pulse.own.jets.complete.Gamma
        logP=source.FunctionRef(g,Gamma.pulse.preRp.operator.background['parameters']['logP'])
        logLambda=g.add(g.mul(c(4),logP),c(1000))
        Lambda,epsilon=g.unary('exp',logLambda),g.unary('exp',g.neg(logLambda))
        self.parameters.update(logP=logP,logLambda=logLambda,Lambda=Lambda,epsilon=epsilon,delta=Gamma.schedule['delta'])
        de=a.fixed(self.parameters['delta']);one=a.fixed(1);Z=c3.C3Function(z,g.one,g.zero,g.zero)
        U0=a.add(a.scale(Z,4),a.fixed(self.parameters['j']));d=a.add(one,a.neg(a.mul(Z,Z)))
        L=a.add(one,a.neg(a.mul(de,a.mul(Z,Z))))
        H=a.add(a.scale(a.mul(a.add(one,a.neg(de)),Z),'1/2'),a.mul(d,U0))
        sigma=g.quotient(self.parameters['j'],c(500),'positive fixed scalar 500')
        gradient=a.div(a.mul(L,H),a.add(a.mul(H,H),a.fixed(g.mul(sigma,sigma))))
        chi=a.div(a.mul(H,H),a.add(a.mul(H,H),a.fixed(g.mul(sigma,sigma))))
        v=g.symbol('original_true_core_gauge_v')
        gv=g.node('function_substitution',expression=gradient.value.node,variable=z.node,value=v.node,
            Z_independent_substitution=True)
        G=controls.integral(g,gv,g.nodes[v.node]['name'],self.parameters['anchor'],z,
            original_H_root_anchor=True,anchor_is_exact_root_not_disk_center=True,
            actual_anchored_gauge_integral=True)
        Gjet=c3.C3Function(G,gradient.value,gradient.Z,gradient.ZZ)
        logF=a.add(a.fixed(g.neg(self.parameters['logC'])),a.neg(a.scale(Gjet,Lambda)))
        diff=pulse.own.jets.OrdinaryZDerivative(g);diff.seed_jet(logF)
        F0value=g.unary('exp',logF.value)
        F0=diff.lift(source.C1Function(F0value,g.mul(F0value,logF.Z)))
        P0=a.scale(self.pulse.own.jets.native['P0'],g.unary('exp',g.mul(c(2),logP)))
        def component(quantity,i=0,k=0):
            return g.node('current_original_unique_analytic_core_component',quantity=quantity,
                radial_order=i,ordinary_slow_Z_derivative_order=k,rho=self.rho.node,Z_variable=z.node,
                source_family=self.family,actual_twenty_term_operator=True,
                unique_analytic_Banach_solution=True,axis_conditions='Phi(0,Z)=1; Psi(0,Z)=0',
                component_definition=('mean_Psi=integral_0^1 Psi(rho*t,Z) dt' if quantity=='mean_Psi'
                    else 'Phi=Phi_model+T_theta; Psi=B*rho/2+T_z'),
                source_provider=ast(self.component_enclosure),
                radial_recurrence_provider=ast(operator.next_rows),
                actual_source_provider=ast(operator.source_terms),
                actual_green_inverse_provider=ast(operator.green_inverse),
                admitted_solution_not_a_model_or_enclosure_center=True)
        self.components={key:c3.C3Function(*(component(key,0,k) for k in range(4)))
            for key in ('Phi','Psi','mean_Psi')}
        Phi,Psi,M=[self.components[key] for key in ('Phi','Psi','mean_Psi')]
        t=g.symbol('original_true_core_pressure_v')
        phi_v=g.node('function_substitution',expression=Phi.value.node,variable=self.rho.node,value=t.node,
            Z_independent_substitution=True)
        pi=controls.integral(g,g.mul(phi_v,phi_v),g.nodes[t.node]['name'],g.zero,self.rho,
            actual_same_Phi_squared_primitive=True,axis_pressure_separate=True)
        # The source system retains the exact primitive derivative, including
        # both F0 derivatives. It never asks the enclosing oracle to divide
        # its scaled pressure range by a possibly zero-lower F0^2 range.
        phiz_v=g.node('function_substitution',expression=Phi.Z.node,variable=self.rho.node,value=t.node,
            Z_independent_substitution=True)
        piZ=controls.integral(g,g.mul(c(2),phi_v,phiz_v),g.nodes[t.node]['name'],g.zero,self.rho)
        S=a.mul(F0,F0);Pcal=g.mul(S.value,pi);PcalZ=g.add(g.mul(S.Z,pi),g.mul(S.value,piZ))
        values=dict(s=self.rho,z=z,delta=de.value,L=L.value,d=d.value,epsilon=epsilon,
            g=gradient.value,S=S.value,U0=U0.value,U0_z=U0.Z,P0=P0.value,P0_z=P0.Z,
            H0=H.value,W0=g.sub(g.one,g.add(g.mul(g.sub(g.one,de.value),z,U0.value),g.mul(d.value,U0.Z))),
            Phi=Phi.value,Phi_s=component('Phi',1,0),Phi_z=Phi.Z,
            Psi=Psi.value,Psi_s=component('Psi',1,0),Psi_z=Psi.Z,
            M=M.value,M_z=M.Z,Pcal=Pcal,Pcal_z=PcalZ)
        self.terms=map_terms(g,values)
        U=a.add(U0,a.scale(Psi,epsilon));mean=a.add(U0,a.scale(M,epsilon))
        self.functions=dict(Phi=Phi,Psi=Psi,mean_Psi=M,U0=U0,gradient=gradient,chi=chi,
            G=Gjet,logF0=logF,F0=F0,P0=P0,U=U,m=mean,
            Q=pulse.own.radial_Q(g,de.value,U,mean))
        self.definition=dict(actual_source_terms=self.terms,primitive_Phi_squared=pi,
            primitive_Phi_squared_Z=piZ,physical_Pcal=Pcal,physical_Pcal_Z=PcalZ,
            source_variables=values,angular_equation='2(rho Phi_rr+2 Phi_r)+chi Phi=epsilon Etheta',
            B=source.expression(g,operator.production_expressions()['B'],values),
            axial_equation='2(rho Psi_rr+Psi_r)=B+epsilon Ez',
            actual_twenty_term_source_and_Green_inverse_retained=True,
            G_anchored_at_unique_H_root=True,Lambda_epsilon_exactly_one=True,
            same_implicit_native_P0=True,physical_Uz_has_no_extra_Pstar=True,
            exact_Phi_and_Psi_distinct_from_numeric_enclosures=True)
        self._U0=U0
        self.hashes[Path(__file__).name]=sha(Path(__file__).name)
        self.prefix=copy.deepcopy(g.nodes);self.snapshot=self._snapshot();self._hashes=dict(self.hashes)
        self.acceptance_loaded=False
        if oracle is not None:self.attach_oracle(oracle)
        if require_checked:
            receipt,raw=parent.linear.read(RECEIPT),parent.linear.read(NAME)
            if (not receipt.get('all_passed') or not receipt.get(GATE) or any(receipt.get(key) for key in OPEN)
                    or receipt['source_family']!=self.family or raw['exact_graph_nodes']!=g.nodes
                    or raw['core_fixed_point_bindings']!=self.snapshot):raise ValueError('Checked core fixed-point bridge required')
            for name,digest in receipt['input_hashes'].items():
                if sha(name)!=digest:raise ValueError('Changed true core bridge '+name)
            self.acceptance_loaded=True

    def _pin(self,name,digest):
        if sha(name)!=digest or (name in self.hashes and self.hashes[name]!=digest):raise ValueError('Changed core dependency '+name)
        self.hashes[name]=digest

    def _snapshot(self):return pulse.own.encode(dict(functions=self.functions,parameters=self.parameters,
        definition=self.definition,components=self.components))

    def assert_graph(self):
        self.pulse.assert_graph()
        if self.g.nodes[:len(self.prefix)]!=self.prefix or self._snapshot()!=self.snapshot or self.hashes!=self._hashes:
            raise ValueError('True original core function graph changed')

    def attach_oracle(self,oracle=None):
        oracle=oracle or common.CurrentCoreCommonFixedPoint()
        if type(oracle) is not common.CurrentCoreCommonFixedPoint or not oracle.acceptance_loaded:
            raise ValueError('Actual checked nonlinear core oracle required')
        identity=dict(actual_five_defect_family_sha256=oracle.family,
            implicit_source_sha256=oracle.source_sha,datum_enclosure_sha256=oracle.datum_sha)
        if identity!=self.family:raise ValueError('Core oracle and pulse source families differ')
        for name,digest in oracle.hashes.items():
            if sha(name)!=digest or self.hashes.get(name,digest)!=digest:raise ValueError('Core oracle dependency differs')
        self.oracle_parameter_identity=self.parameter_identity(oracle)
        self._oracle=oracle
        return oracle

    @common.source_precision
    def parameter_identity(self,oracle):
        """Exact graph recipes and direct source-datum identity, not box overlap."""
        g=self.g;c=g.constant;p=self.parameters;ctx=oracle.ctx;core=oracle.core
        expected_logP=g.add(g.unary('exp',c(40)),c(11))
        expected_delta=g.unary('exp',g.neg(g.add(g.mul(c(4),p['logP']),c(30))))
        assert p['logP']==expected_logP and p['delta']==expected_delta
        assert p['logLambda']==g.add(g.mul(c(4),p['logP']),c(1000))
        assert p['epsilon']==g.unary('exp',g.neg(p['logLambda']))
        exactLP=ctx.exp(40)+11;exactLL=4*exactLP+1000
        contains=lambda box,value:(common.endpoints(box)[0]<=common.endpoints(value)[0]
            and common.endpoints(value)[1]<=common.endpoints(box)[1])
        checks={key:contains(getattr(core,key),value) for key,value in
            (('logP',exactLP),('logLambda',exactLL),('delta',ctx.exp(-4*exactLP-30)),
             ('Lambda',ctx.exp(exactLL)),('epsilon',ctx.exp(-exactLL)))}
        if not all(checks.values()):raise ValueError('Oracle does not enclose the exact pulse/core parameter recipes '+str(checks))
        amplitude=oracle.source.amplitude
        read=lambda row,key:common.read_interval(ctx,row[key])
        scalar_checks=dict(
            j_matches_admitted_source=core.j._mpi_==read(core.records['core_transfer'],'required_j')._mpi_,
            logC_matches_original_selected_source=core.logC._mpi_==read(core.records['physical_norm_family'],'selected_logCstar')._mpi_,
            sigma_is_same_j_over_500=core.sigma._mpi_==(core.j/500)._mpi_,
            amplitude_uses_same_core=amplitude.core is core,
            amplitude_uses_same_nonlinear_rebuild=amplitude.rebuild is oracle.rebuild)
        scalar_checks.update({key+'_matches_anchored_amplitude':
            amplitude.params[key]._mpi_==getattr(core,key)._mpi_ for key in ('j','delta','sigma')})
        for key,path in (('j','common.core.j'),('logC','common.core.logC'),('anchor','common.source.amplitude.roots')):
            node=g.nodes[p[key].node]
            scalar_checks[key+'_leaf_uses_this_exact_source']=(
                node['operation']=='current_original_source_parameter'
                and node['source_path']==path and node['source_family']==self.family
                and node['exact_fixed_source_parameter'] and node['parameter_enclosure_not_selected_value'])
        roots=amplitude.roots;lo,hi=common.endpoints(roots['anchor_interval'])
        jhi=common.endpoints(core.j)[1]
        scalar_checks['unique_H_root_in_original_bracket']=(
            roots['unique_real_anchor_certified'] and roots['uniform_interval_parameter_root_certificates']
            and roots['nominal_parameters_used_only_for_provisional_centers']
            and -jhi<=lo<=hi<=0 and common.endpoints(roots['anchor_Hprime_lower'])[0]>0)
        if not all(scalar_checks.values()):raise ValueError('Actual core scalar/anchor identity differs '+str(scalar_checks))
        binding=parent.linear.read(hierarchy.native.frame.native.NAME)['current_P0_source_binding']
        actual_AST=hashlib.sha256(ast.dump(ast.parse(inspect.getsource(
            core.datum.normalized_jets.__func__).strip())).encode()).hexdigest()
        flags=('same_fourteen_continuous_atoms_and_analytic_flatten_function',
            'current_order6_first_six_and_native_order5_same_defining_rows',
            'common_P0_kept_independent_of_cumulative_pressure')
        if (not binding['passed'] or not all(binding[key] for key in flags)
                or actual_AST!=binding['inherited_normalized_jets_AST_sha256']
                or len(core.datum.stages)!=14 or core.datum.datum_sha!=self.family['datum_enclosure_sha256']
                or not contains(core.datum.parameters.logPstar,exactLP)):
            raise ValueError('Pulse native P0 and actual core implicit datum differ')
        return dict(passed=True,exact_graph_logP_delta_logLambda_epsilon_recipes_verified=True,
            actual_core_encloses_same_exact_recipes=checks,
            actual_source_scalar_and_unique_anchor_identity=scalar_checks,
            selected_j_logC_and_anchor_are_from_this_actual_core=True,
            same_fourteen_atom_P0_callable_AST_and_datum_hash=True,
            interval_overlap_not_used_as_a_function_identity=True)

    @common.source_precision
    def component_enclosure(self,quantity,rho,Z,radial_order=0,axial_order=0,degree=6,depth=3):
        """Ordinary derivatives of the actual solution, as directed ranges."""
        if not self._oracle:self.attach_oracle()
        oracle=self._oracle;c=oracle.ctx;r=c.mpf(rho);z=c.mpf(Z);i,k=radial_order,axial_order
        if quantity not in ('Phi','Psi','mean_Psi') or min(i,k)<0 or i+k>5 or k>depth:
            raise ValueError('Requested actual core derivative unavailable')
        if common.endpoints(z)[0]<-1 or common.endpoints(z)[1]>1 or degree<6:
            raise ValueError('Actual admitted core domain/degree required')
        key=(z._mpi_,degree,depth)
        if key not in self._packets:self._packets[key]=oracle.rebuild.rebuild(z,degree=degree,depth=depth)
        packet=self._packets[key];profile=oracle.rebuild.profile(packet,r,radial_order=i,axial_order=k)
        value=component_polynomial(packet,r,oracle.core.epsilon,quantity,i,k,c.mpf)
        tail=profile['infinite_radial_tail_bounds']['Phi' if quantity=='Phi' else 'Uz']
        if quantity!='Phi':tail=tail/oracle.core.epsilon
        bound=common.endpoints(tail)[1]
        return dict(quantity=quantity,radial_order=i,axial_order=k,finite_polynomial_enclosure=value,
            infinite_radial_tail_upper=tail,enclosure=value+c.mpf([-bound,bound]),
            all_rows_are_ordinary_derivatives=True,actual_twenty_term_nonlinear_solution_enclosed=True,
            numerical_center_not_selected=True,physical_Uz_has_no_extra_Pstar=True,
            source_family=self.family,radial_degree=degree,axial_depth=depth)


def run():
    began=time.monotonic();owner=CurrentOriginalN1CoreFixedPointFunctions(require_checked=False)
    raw=dict(source_family=owner.family,exact_graph_nodes=owner.g.nodes,core_fixed_point_bindings=owner.snapshot,
        actual_operator_term_counts={key:len(values) for key,values in owner.terms.items()},
        unique_actual_analytic_solution_function_definitions=True,component_oracle_is_lazy=True,
        **dict.fromkeys((GATE,)+OPEN,False),input_hashes=owner.hashes,execution_seconds=time.monotonic()-began)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(raw,separators=(',',':'))+'\n').encode(),mtime=0))
    print('ORIGINAL_N1_TRUE_CORE_FUNCTIONS_DEFINED',len(owner.g.nodes),flush=True)
    return owner,raw


if __name__=='__main__':run()
