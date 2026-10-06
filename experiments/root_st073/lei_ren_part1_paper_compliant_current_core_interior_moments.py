"""Same fixed-point six interior core moments, nonsingular axis and Euler jets."""
import ast
import gzip
import json
import math
from pathlib import Path
from types import SimpleNamespace

import mpmath as mp
import sympy as s
from lei_ren_part1_paper_compliant_current_bridge_background_tensor import (
    CurrentBridgeBackgroundTensor,HERE,PREFIX,sha,pack,encode,endpoints,source_precision,accepted,_verify_hashes)
from lei_ren_part1_paper_compliant_current_core_scaled_swirl_source import CurrentCoreScaledSwirlSource
from lei_ren_part1_paper_compliant_current_core_nonlinear_operator import CurrentCoreNonlinearOperator
from lei_ren_part1_paper_compliant_current_core_common_fixed_point import CurrentCoreCommonFixedPoint
from lei_ren_part1_paper_compliant_current_core_first_interface import defining_boundary_bindings,binding,function
from lei_ren_part1_paper_compliant_core_integral_atoms import finite_atom_coefficients
from lei_ren_part1_paper_compliant_core_physical_field import symmetric,intersection
from lei_ren_part1_paper_interval_taylor import IntervalTaylor

NAME=PREFIX+'current_core_interior_moments.json.gz'
RECEIPT=PREFIX+'current_core_interior_moments_check.json'
GATES=('current_core_six_interior_moment_functions_and_tails_available',
       'current_core_interior_moment_rho4_axial6_and_Euler4_available',
       'current_core_interior_normalized_axis_and_inlet_source_identities_certified')
OPEN=('current_core_full_background_tensor_available','current_core_bridge_completed_tensor_join_certified',
      'core_axis_tensor_remainder_limits_certified','global_completed_tensor_admissibility',
      'admissible_stress_lift_constructed','full_background_NS_validation','physical_energy_integral_certified',
      'full_point_physical_field_evaluation','actual_point_moment_history_recovered',
      'independently_bounded_flat_remainder','temporal_recursion')
SPEC={'H':('Phi',None,1,2),'M':('Uz',None,0,1),'K':('Phi','Uz',1,2),
      'A':('Uz','Uz',0,1),'B':('Phi','Phi',1,1),'C':('Phi','Phi',0,1)}
INDICES=[(i,k) for i in range(5) for k in range(7)]
VIEWS={'whole':([-1,1],[0,4]),'axis':([-1,1],0),'inlet':([-1,1],4),
       'fresh':('.359','.731'),'fresh_axis':('-.317',0),'fresh_inlet':('.537',4)}


def key(i,k):return 'rho'+str(i)+'_Z'+str(k)


def stirling(q,j):
    if q==0:return int(j==0)
    if j==0 or j>q:return 0
    return j*stirling(q-1,j)+stirling(q-1,j-1)


def convolve(c,left,right,order=6):
    return [[sum((left[n][ell]*right[m][k-ell]
             for n in range(len(left)) for m in range(len(right)) if n+m==total
             for ell in range(k+1)),c.mpf(0))
             for k in range(order+1)]
             for total in range(len(left)+len(right)-1)]


def integrated_polynomial(c,rows,rho,i,k,weight,normalization):
    return normalization*sum((row[k]*(math.factorial(n)//math.factorial(n-i))
            *rho**(n-i)/(n+weight+1) for n,row in enumerate(rows) if n>=i),c.mpf(0))


def mixed_product_tail(c,left_norm,right_norm,left_tail,right_tail,i,k):
    # Z entries are Taylor coefficients; rho entries are ordinary derivatives.
    return sum((math.comb(i,a)*(left_norm[a,ell]*right_tail[i-a,k-ell]
            +left_tail[a,ell]*right_norm[i-a,k-ell]+left_tail[a,ell]*right_tail[i-a,k-ell])
            for a in range(i+1) for ell in range(k+1)),c.mpf(0))


def relative_ratios(c,ell,multiplier=1):
    # ell[k] is the Taylor coefficient of (log F0)' at the same basepoint.
    rows=[c.mpf(1)]
    for n in range(1,7):
        rows.append(sum((multiplier*ell[j]*rows[n-1-j] for j in range(n)),c.mpf(0))/n)
    return [value*math.factorial(k) for k,value in enumerate(rows)]


def source_theorem():
    rho,t,z=s.symbols('rho t Z',real=True)
    phi=s.Function('Phi');uz=s.Function('Uz')
    integrals={name:s.Integral(norm*t**w*(phi(rho*t,z) if left=='Phi' else uz(rho*t,z))
            *((phi(rho*t,z) if right=='Phi' else uz(rho*t,z)) if right else 1),(t,0,1))
            for name,(left,right,w,norm) in SPEC.items()}
    # The transformed [0,1] integrals equal the six original physical primitives.
    rr,eps,F0=s.symbols('R epsilon F0',positive=True)
    raw=dict(I_V=eps*rho*s.Symbol('M'),I_V2=eps*rho*s.Symbol('A'),
             I_RF=eps**2*F0*rho**2*s.Symbol('H')/2,
             I_RFV=eps**2*F0*rho**2*s.Symbol('K')/2,
             I_RF2=eps**2*F0**2*rho**2*s.Symbol('B'),
             I_F2=eps*F0**2*rho*s.Symbol('C'))
    # Arbitrary finite coefficient arrays test the unchanged original endpoint program.
    P=[[s.Symbol('p_'+str(n)+'_'+str(k)) for k in range(3)] for n in range(4)]
    U=[[s.Symbol('u_'+str(n)+'_'+str(k)) for k in range(3)] for n in range(4)]
    class Exact:
        @staticmethod
        def mpf(value):return s.Rational(value)
    c=Exact();original=finite_atom_coefficients(c,P,U,2)
    density={'Phi':P,'Uz':U,'PhiUz':convolve(c,P,U,2),'PhiPhi':convolve(c,P,P,2),'UzUz':convolve(c,U,U,2)}
    endpoint={}
    for name,(left,right,w,norm) in SPEC.items():
        rows=density[left+right] if right else density[left]
        endpoint[name]=[s.expand(integrated_polynomial(c,rows,s.Integer(4),0,k,w,norm)-original[name][k])==0 for k in range(3)]
    # Exact differential operator identity on an arbitrary smooth function.
    f=s.Function('f')(rho);ordinary={}
    current=f
    for q in range(5):
        ordinary[str(q)]=s.simplify(current-sum(stirling(q,j)*rho**j*s.diff(f,rho,j) for j in range(q+1)))==0
        current=s.expand(rho*s.diff(current,rho))
    if not all(all(v) for v in endpoint.values()) or not all(ordinary.values()):raise ValueError('Core endpoint/Euler source identity failed')
    original_boundary=defining_boundary_bindings()
    binding('compliant_core_physical_field','profiles','coefficients',
        "(2*source['z']/source['L'],source['z']*(-(1-self.delta))/source['L'],-source['d']/source['L'])")
    binding('compliant_core_physical_field','profiles','grids[Q][index]',
        "sum((math.comb(k,j)*math.factorial(j)*(coefficients[0][j]*vel[gridkey(i,k-j)]+coefficients[1][j]*mean[gridkey(i,k-j)]+coefficients[2][j]*mean[gridkey(i,k-j+1)]) for j in range(k+1)),c.mpf(0))")
    binding('compliant_current_core_interior_moments','evaluate','coefficients',
        '(2*zj/L,-zj*(1-self.core.delta)/L,-d/L)')
    fn=function('compliant_current_core_interior_moments','evaluate')
    calls=[node for node in ast.walk(fn) if isinstance(node,ast.Call) and ast.unparse(node.func)=='values.append']
    expected=ast.parse("math.comb(k,b)*math.factorial(b)*(coefficients[0][b]*uv+coefficients[1][b]*grids['M'][key(i,k-b)]+coefficients[2][b]*grids['M'][key(i,k-b+1)])",mode='eval').body
    if len(calls)!=1 or ast.dump(calls[0].args[0])!=ast.dump(expected):raise ValueError('Original coefficientwise core Q source product changed')
    # Replay the original coefficients on arbitrary exact Z/delta before bounds.
    dt=s.Symbol('delta',real=True);LL=1-dt*z*z;dd=1-z*z
    old=function('compliant_core_physical_field','profiles')
    assignment=next(node for node in ast.walk(old) if isinstance(node,ast.Assign) and any(ast.unparse(v)=='coefficients' for v in node.targets))
    original_Q=eval(compile(ast.Expression(assignment.value),'<original core Q coefficients>','eval'),
        dict(source=dict(z=z,L=LL,d=dd),self=SimpleNamespace(delta=dt)))
    current_Q=(2*z/LL,-z*(1-dt)/LL,-dd/LL)
    Q_bindings=[s.cancel(a-b)==0 for a,b in zip(original_Q,current_Q)]
    if not all(Q_bindings):raise ValueError('Original core Q source coefficients differ')
    return dict(normalized_six_integral_definitions={k:str(v) for k,v in integrals.items()},
        original_physical_integral_recovery={k:str(v) for k,v in raw.items()},
        original_rho4_finite_atom_identities=endpoint,ordinary_logR_Euler4_operator_identities=ordinary,
        differentiated_normalized_integral='d_rho^i d_Z^k J = normalization*int_0^1 t^(i+weight)*d_rho^i d_Z^k density(rho*t,Z) dt',
        integral_tail_weight='normalization/(i+weight+1) times uniform mixed density tail; Z factorial included once',
        raw_axis_valuations=dict(I_V=1,I_V2=1,I_F2=1,I_RF=2,I_RFV=2,I_RF2=2),
        axis_normalized_values=dict(H='1',M='4Z+j',K='4Z+j',A='(4Z+j)^2',B='1/2',C='1'),
        pressure_identity='P=Pstar^2*P0+epsilon_core*F0base^2*dress(rho*C); P0 separate from C',
        original_core_first_boundary_bindings=original_boundary,
        original_core_Q_production_AST_and_coefficient_identities=Q_bindings,
        same_original_Q_mean_sixth_derivative_source_retained=True,
        source_functions_not_chosen_from_enclosures=True,passed=True)


class CurrentCoreInteriorMoments:
    @source_precision
    def __init__(self,bridge_tensor=None,common=None,require_checked=True):
        self.bridge_tensor=bridge_tensor if bridge_tensor is not None else CurrentBridgeBackgroundTensor()
        if not self.bridge_tensor.acceptance_loaded:raise ValueError('Checked current complete bridge tensor required')
        self.bridge=self.bridge_tensor.bridge
        if common is None:
            scaled=CurrentCoreScaledSwirlSource(bridge=self.bridge)
            common=CurrentCoreCommonFixedPoint(operator=CurrentCoreNonlinearOperator(source=scaled))
        self.common=common;self.core=common.core;self.rebuild=common.rebuild;self.ctx=self.core.ctx
        self.family=common.family;self.source=common.source_sha;self.datum_sha=common.datum_sha
        self.assert_graph();self.proof=source_theorem();self.cache={};self.rows_cache={};self.tail_cache={}
        first_name=PREFIX+'current_core_first_interface_check.json'
        first=accepted(first_name,self.family,self.source,'current_core_bridge_functional_mixed4_join_certified');_verify_hashes(first)
        if first['datum_enclosure_sha256']!=self.datum_sha or not first['all_current_bridge_functional_interfaces_certified']:
            raise ValueError('Same core/first source interface receipt required')
        self.proof.update(consumed_same_fixed_point_model_bindings=common.model_bindings,
            consumed_same_fixed_point_finite_tail_bindings=common.consumer_bindings,
            consumed_core_first_functional_receipt=first_name,
            core_first_source_receipt_does_not_admit_completed_tensor=True)
        self.hashes={**self.bridge_tensor.hashes,**common.hashes}
        self.hashes.update(first['input_hashes']);self.hashes[first_name]=sha(first_name)
        for name in (Path(__file__).name,PREFIX+'current_bridge_background_tensor_check.json',
                PREFIX+'core_integral_atoms.py',PREFIX+'core_physical_field.py',
                PREFIX+'core_coefficient_rebuild.py',PREFIX+'current_core_first_interface.py'):
            self.hashes[name]=sha(name)
        self.acceptance_loaded=False
        if require_checked:
            receipt=accepted(RECEIPT,self.family,self.source,GATES[0]);_verify_hashes(receipt)
            if receipt['datum_enclosure_sha256']!=self.datum_sha or not all(receipt[k] for k in GATES) or any(receipt[k] for k in OPEN):
                raise ValueError('Core interior moment receipt source/scope differs')
            self.hashes.update(receipt['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    def assert_graph(self):
        if not (self.common.acceptance_loaded and self.bridge_tensor.acceptance_loaded
                and self.common.source.bridge is self.bridge
                and self.rebuild is self.bridge.upstream.comparison.atoms.rebuild
                and self.core is self.bridge.upstream.core is self.bridge.core.original
                and self.ctx is self.bridge.ctx
                and (self.family,self.source,self.datum_sha)==(self.bridge.family,self.bridge.source,self.bridge.datum_sha)):
            raise ValueError('Same admitted nonlinear core/bridge/moment graph required')

    def rows(self,Z):
        z=self.ctx.mpf(Z)
        if not all(mp.isfinite(v) for v in endpoints(z)) or endpoints(z)[0]<-1 or endpoints(z)[1]>1:raise ValueError('Finite Z in[-1,1] required')
        if z._mpi_ not in self.rows_cache:
            packet=self.rebuild.rebuild(z,degree=24,depth=6);P=[row[:7] for row in packet['rows']['A'][:25]];U=[row[:7] for row in packet['rows']['Uz'][:25]]
            density={'Phi':P,'Uz':U,'PhiUz':convolve(self.ctx,P,U),'PhiPhi':convolve(self.ctx,P,P),'UzUz':convolve(self.ctx,U,U)}
            self.rows_cache[z._mpi_]=(packet,density)
        return self.rows_cache[z._mpi_]

    def norms_and_tails(self,packet,density,rmax):
        c=self.ctx;cache=(packet['Z']._mpi_,rmax._mpi_)
        if cache in self.tail_cache:return self.tail_cache[cache]
        norms={label:{} for label in ('Phi','Uz')};tails={label:{} for label in norms}
        for i,k in INDICES:
            original=self.rebuild.profile(packet,c.mpf([0,endpoints(rmax)[1]]),radial_order=i,axial_order=k)
            for label in norms:
                norms[label][i,k]=sum((abs(row[k])*(math.factorial(n)//math.factorial(n-i))*rmax**(n-i)
                    for n,row in enumerate(density[label]) if n>=i),c.mpf(0))
                tails[label][i,k]=original['infinite_radial_tail_bounds'][label]/math.factorial(k)
        self.tail_cache[cache]=(norms,tails);return norms,tails

    @source_precision
    def evaluate(self,Z,rho):
        self.assert_graph();c=self.ctx;r=c.mpf(rho)
        if not all(mp.isfinite(v) for v in endpoints(r)) or endpoints(r)[0]<0 or endpoints(r)[1]>4:raise ValueError('Finite prescribed core rho in[0,4] required')
        packet,density=self.rows(Z);z=packet['Z'];cache=(z._mpi_,r._mpi_)
        if cache in self.cache:return self.cache[cache]
        rmax=c.mpf(endpoints(r)[1]);norms,tails=self.norms_and_tails(packet,density,rmax)
        finite={name:{} for name in SPEC};errors={name:{} for name in SPEC};grids={name:{} for name in SPEC}
        for name,(left,right,w,norm) in SPEC.items():
            rows=density[left+right] if right else density[left]
            for i,k in INDICES:
                label=key(i,k);value=integrated_polynomial(c,rows,r,i,k,w,norm)
                error=(mixed_product_tail(c,norms[left],norms[right],tails[left],tails[right],i,k)
                        if right else tails[left][i,k])*norm/(i+w+1)
                finite[name][label]=value;errors[name][label]=error
                grids[name][label]=(value+symmetric(c,error))*math.factorial(k)
        for name,floor,ceiling in (('H',self.core.phi_floor,self.core.phi_ceiling),
                ('B',self.core.phi_floor**2/2,self.core.phi_ceiling**2/2),
                ('C',self.core.phi_floor**2,self.core.phi_ceiling**2)):
            grids[name][key(0,0)]=intersection(c,grids[name][key(0,0)],c.mpf([endpoints(floor)[0],endpoints(ceiling)[1]]))
        logR={name:{'y'+str(q)+'_Z'+str(k):sum((stirling(q,j)*r**j*grids[name][key(j,k)] for j in range(q+1)),c.mpf(0))
                for q in range(5) for k in range(7)} for name in SPEC}
        ell=[self.core.Lambda*v for v in packet['fixed']['ell_Z_taylor'][:6]]
        ratios=relative_ratios(c,ell);ratios2=relative_ratios(c,ell,2)
        datum=self.core.datum.normalized_jets(z,6)['normalized_pressure_coefficients']
        increment={}
        for i,k in INDICES:
            increment[key(i,k)]=sum((math.comb(k,j)*ratios2[j]*(r*grids['C'][key(i,k-j)]
                +(i*grids['C'][key(i-1,k-j)] if i else 0)) for j in range(k+1)),c.mpf(0))
        # Q uses a sixth axial derivative of the radial mean, retained above.
        zj=IntervalTaylor.variable(c,z,6);L=1-zj*zj*self.core.delta;d=1-zj*zj
        coefficients=(2*zj/L,-zj*(1-self.core.delta)/L,-d/L);Q={}
        for i in range(5):
            for k in range(6):
                # The full coefficientwise product uses all lower axial rows.
                values=[]
                for b in range(k+1):
                    uv=sum((packet['rows']['Uz'][n][k-b]*(math.factorial(n)//math.factorial(n-i))*r**(n-i)
                        for n in range(i,25)),c.mpf(0))
                    uv=(uv+symmetric(c,tails['Uz'][i,k-b]))*math.factorial(k-b)
                    values.append(math.comb(k,b)*math.factorial(b)*(coefficients[0][b]*uv
                        +coefficients[1][b]*grids['M'][key(i,k-b)]+coefficients[2][b]*grids['M'][key(i,k-b+1)]))
                Q[key(i,k)]=sum(values,c.mpf(0))
        result=dict(Z=z,rho=r,radial_degree=24,axial_depth=6,
            exact_normalized_integral_source=self.proof['normalized_six_integral_definitions'],
            finite_integrated_axial_Taylor_coefficient_grids=finite,
            directed_integrated_axial_Taylor_tail_bounds=errors,
            ordinary_rho4_axial6_moment_grids=grids,ordinary_logR4_axial6_moment_grids=logR,
            original_density_finite_absolute_coefficient_bounds={label:{key(*ik):v for ik,v in values.items()} for label,values in norms.items()},
            original_density_directed_radial_axial_Taylor_tail_bounds={label:{key(*ik):v for ik,v in values.items()} for label,values in tails.items()},
            original_analytic_axis_pressure_axial6_Taylor_coefficients=list(datum),
            original_F0_relative_axial6_ordinary_derivatives=ratios,original_F0_squared_relative_axial6_ordinary_derivatives=ratios2,
            normalized_pressure_increment_rho4_axial6=increment,normalized_radial_Q_rho4_axial5=Q,
            original_pressure_axis_and_increment_retained_separately=True,
            raw_physical_integrals_recovered_with_original_epsilon_and_F0=self.proof['original_physical_integral_recovery'],
            original_epsilon_positive_log=self.ctx.ln(self.core.epsilon),
            exact_F0_source='exp(-selected_logCstar-Lambda*G(Z))',
            actual_nonlinear_fixed_point_and_all_density_product_tails_retained=True,
            original_mean_axial6_retained_before_Q_axial5=True,
            normalized_axis_uses_no_division_by_rho=True,
            physical_axis_tensor_limits_not_admitted=True,source_enclosures_not_resolved_point_values=True,
            finite_polynomial_is_not_complete_solution=True,**dict.fromkeys(OPEN,False))
        self.cache[cache]=result;return result

    def manifest(self):
        return dict(actual_five_defect_family_sha256=self.family,implicit_source_sha256=self.source,
            datum_enclosure_sha256=self.datum_sha,source_theorem=self.proof,
            original_core_common_fixed_point_receipt=PREFIX+'current_core_common_fixed_point_check.json',
            existing_completed_tensor_inventory=dict(regions=32,adjacent=31,internal=10),
            current_velocity_pressure_atlas_inventory=dict(adjacent=14,internal=8),
            **dict.fromkeys(GATES+OPEN,False),input_hashes=self.hashes)


@source_precision
def run(field=None):
    field=field if field is not None else CurrentCoreInteriorMoments(require_checked=False)
    result=field.manifest();result['current_core_interior_views']={}
    for name,args in VIEWS.items():
        result['current_core_interior_views'][name]=field.evaluate(*args)
        print('Build same fixed-point core interior six moments: '+name,flush=True)
    (HERE/NAME).write_bytes(gzip.compress((json.dumps(encode(pack(result)),indent=2)+'\n').encode(),mtime=0))
    return result


if __name__=='__main__':run()
