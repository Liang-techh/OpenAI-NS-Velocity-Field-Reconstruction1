"""Actual twenty-term equations, correlated coefficient units and core oracle."""
import json
from pathlib import Path
import time
import scipy.integrate
import sympy as sy

import lei_ren_part1_paper_compliant_current_original_n1_core_fixed_point_functions as current
import lei_ren_part1_paper_compliant_current_original_n1_preRp_moment_functions_check as arithmetic


def equation_fixture():
    g=current.source.FunctionTransportGraph()
    names='s z delta L d epsilon g S U0 U0_z P0 P0_z H0 W0 Phi Phi_s Phi_z Psi Psi_s Psi_z M M_z Pcal Pcal_z'.split()
    symbols={name:sy.Symbol(name) for name in names};values={name:g.symbol(name) for name in names}
    terms=current.map_terms(g,values);v=symbols
    r,z,de,L,d,eps,gg,S,U,Uz,P,Pz,H0,W0,phi,phir,phiz,psi,psir,psiz,M,Mz,C,Cz=[v[name] for name in names]
    W=W0-eps*((1-de)*z*M+d*Mz);H=H0+eps*d*psi;u=U+eps*psi
    angular=eps*(W*(r*phir+phi)+H*phiz+de*phi/2-de*z*phi*u)/L-gg*H*phi/L
    axial=(W*eps*r*psir+H*(Uz+eps*psiz)+(1+de)*(u-2*z*u*u)/2
        +d*(Pz+eps*Cz)-2*(1+de)*z*(P+eps*C)-2*z*eps*r*S*phi*phi)/L
    B=(H0*Uz+(1+de)*(U-2*z*U*U)/2+d*Pz-2*(1+de)*z*P)/L
    totals={key:sum(arithmetic.symbolic(g,ref) for ref in rows) for key,rows in terms.items()}
    assert len(terms['theta'])==len(terms['z'])==10
    assert sy.expand(angular+gg*H0/L*phi-eps*totals['theta'])==0
    assert sy.expand(axial-B-eps*totals['z'])==0
    return dict(passed=True,independent_original_angular_and_axial_equations_recovered=True,
        all_twenty_terms_retained=True,centrifugal_pressure_and_its_Z_derivative_retained=True)


def coefficient_fixture():
    r,z=sy.symbols('rho Z');center=sy.Rational(3,7);point=sy.Rational(4,5);epsilon=sy.Rational(1,9)
    Phi=1+r*(1+z**2)+r**2*(z+z**3)/5+r**3*(1-z)/7
    Psi=r*(2+z**3)+r**2*(1+z**2)/3+r**3*(z-2)/5
    U0=4*z+sy.Rational(7,100);Uz=U0+epsilon*Psi
    mean=sy.integrate(Psi,(r,0,r))/r
    packet=dict(radial_degree=3,rows={})
    for key,expr in (('A',Phi),('Uz',Uz)):
        coefficients=sy.Poly(expr,r)
        packet['rows'][key]=[[sy.diff(coefficients.nth(n),z,k).subs(z,center)/sy.factorial(k)
            for k in range(4)] for n in range(4)]
    for key,expr in (('Phi',Phi),('Psi',Psi),('mean_Psi',mean)):
        for i,k in ((0,0),(1,0),(0,1),(2,1),(0,3)):
            got=current.component_polynomial(packet,point,epsilon,key,i,k,sy.sympify)
            want=sy.diff(expr,r,i,z,k).subs({r:point,z:center})
            assert sy.cancel(got-want)==0,(key,i,k)
    assert current.component_polynomial(packet,0,epsilon,'Psi',0,0,sy.sympify)==0
    assert current.component_polynomial(packet,0,epsilon,'mean_Psi',0,0,sy.sympify)==0
    return dict(passed=True,independent_nonunit_epsilon_finite_polynomial_fixture=True,
        coefficient_to_ordinary_mixed_derivative_conversion_once=True,
        affine_U0_removed_structurally_before_epsilon_division=True,
        exact_radial_average_factor_checked=True,missing_nonlinear_rows_not_zero_padded=True)


def primitive_fixture(owner):
    g=owner.g;r=sy.Symbol(g.nodes[owner.rho.node]['name']);z=sy.Symbol('Z')
    phi=1+r*z+r**2*(1+z**2);F=sy.exp(-z**2/10)
    Phi=owner.components['Phi'];S=owner.definition['source_variables']['S']
    # This finite mathematical substitution checks the actual graph's
    # product/primitive derivative; it is not an original core solution.
    bindings={Phi.value.node:phi,Phi.Z.node:sy.diff(phi,z),S.node:F**2}
    sZ=g.nodes[owner.definition['physical_Pcal_Z'].node]['arguments'][0]
    product=g.nodes[sZ]
    if product['operation']=='product':
        for ref in product['arguments']:
            if ref!=owner.definition['primitive_Phi_squared'].node:bindings[ref]=sy.diff(F**2,z)
    value=arithmetic.symbolic(g,owner.definition['physical_Pcal'],bindings=bindings)
    derivative=arithmetic.symbolic(g,owner.definition['physical_Pcal_Z'],bindings=bindings)
    want=F**2*sy.integrate(phi**2,(r,0,r))
    assert sy.simplify(value-want)==0
    assert sy.simplify(derivative-sy.diff(want,z))==0
    return dict(passed=True,mathematical_fixture_only=True,
        actual_same_Phi_squared_primitive_and_both_F0_derivatives_checked=True)


def actual_oracle_fixture(owner):
    print('ATTACHING_CHECKED_ACTUAL_CORE_ORACLE',flush=True)
    oracle=owner.attach_oracle();c=oracle.ctx
    assert owner.oracle_parameter_identity['passed']
    result={}
    for quantity in ('Phi','Psi','mean_Psi'):
        axis=owner.component_enclosure(quantity,'0','.371')
        want=c.mpf(1 if quantity=='Phi' else 0)
        assert axis['enclosure']._mpi_==want._mpi_
        for k in (0,1,3):
            sample=owner.component_enclosure(quantity,'2','.371',axial_order=k)
            assert sample['all_rows_are_ordinary_derivatives'] and sample['numerical_center_not_selected']
            lo,hi=current.common.endpoints(sample['enclosure']);assert lo<=hi
            assert current.common.endpoints(sample['infinite_radial_tail_upper'])[0]>=0
            if quantity=='Phi' and k==0:assert lo>0 and hi<1
            result[quantity+'_Z'+str(k)]=dict(enclosure=current.common.encode(sample['enclosure']),
                tail=current.common.encode(sample['infinite_radial_tail_upper']))
    assert len(owner._packets)==1
    # Same global gauge anchor: root centers guide disks only; they are
    # never substituted for the exact mathematical source root.
    roots=oracle.source.amplitude.roots
    assert roots['unique_real_anchor_certified'] and roots['G_nonnegative_on_real_domain']
    assert roots['nominal_parameters_used_only_for_provisional_centers']
    return dict(passed=True,actual_source_evaluation=True,Z='.371',rho='2',radial_degree=6,axial_depth=3,
        same_checked_nonlinear_finite_packet_shared=True,all_three_actual_components_Z0_Z1_Z3_enclosed=True,
        exact_axis_normalizations_retained=True,positive_nontrivial_Phi_enclosed=True,
        true_H_root_anchor_retained=True,finite_polynomial_plus_full_tail_not_selected_as_point=True,
        exact_parameter_and_native_P0_identity=owner.oracle_parameter_identity,
        samples=result)


def run(owner=None):
    began=time.monotonic();owner=owner or current.CurrentOriginalN1CoreFixedPointFunctions(require_checked=False)
    owner.assert_graph();raw=current.parent.linear.read(current.NAME)
    assert raw['exact_graph_nodes']==owner.g.nodes and raw['core_fixed_point_bindings']==owner.snapshot
    assert raw['actual_operator_term_counts']=={'theta':10,'z':10}
    assert not any(raw[key] for key in current.OPEN)
    checks=dict(original_equations=equation_fixture(),ordinary_coefficients=coefficient_fixture(),
        primitive=primitive_fixture(owner),actual_oracle=actual_oracle_fixture(owner))
    owner.assert_graph()
    for name,digest in owner.hashes.items():assert current.sha(name)==digest,name
    receipt=dict(all_passed=True,source_family=owner.family,checks=checks,
        **{current.GATE:True},**dict.fromkeys(current.OPEN,False),
        full_inner_cutoff_lower_collar_five_defect_adapter_installed=False,
        uniform_original_core_point_numerical_accuracy_claimed=False,
        source_functions_and_directed_enclosures_are_distinct=True,
        input_hashes={**owner.hashes,current.NAME:current.sha(current.NAME),
            Path(__file__).name:current.sha(Path(__file__).name),
            Path(arithmetic.__file__).name:current.sha(Path(arithmetic.__file__).name)},
        execution_seconds=time.monotonic()-began)
    (current.HERE/current.RECEIPT).write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8',newline='\n')
    print('ORIGINAL_N1_TRUE_CORE_FUNCTIONS_AND_ORACLE_CHECKED',flush=True)
    return receipt


if __name__=='__main__':run()
