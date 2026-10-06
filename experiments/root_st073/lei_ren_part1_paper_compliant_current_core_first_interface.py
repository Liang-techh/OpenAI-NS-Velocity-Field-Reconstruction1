"""Common analytic core -> original first prescribed-shear chart.

The actual atom source is retained. Pressure caps and covering profiles are
only enclosures of that source, not definitions of inlet moments or values.
An independent checker admits the stress-free integral identities and join.
"""
import ast
import copy
import json
import math
from pathlib import Path

import mpmath as mp

from lei_ren_part1_paper_compliant_current_core_common_fixed_point import (
    CurrentCoreCommonFixedPoint,GATE as COMMON_GATE,function,binding,sha)
from lei_ren_part1_paper_compliant_actual_reference_restore_mixed_C4 import accepted
from lei_ren_part1_paper_compliant_macro_signed_integrals import _verify_hashes
from lei_ren_part1_paper_compliant_switch_signed_integrals import source_precision
from lei_ren_part1_paper_schedule_endpoint_enclosures import endpoints
from lei_ren_part1_paper_schedule_endpoint_enclosures_check import encode
from lei_ren_part1_paper_compliant_frozen_comparison_field import dress
from lei_ren_part1_paper_interval_taylor import IntervalTaylor
from lei_ren_part1_paper_compliant_five_moment_repair import pack
import lei_ren_part1_paper_compliant_comparison_point_integrals as comparison_module

HERE=Path(__file__).parent
PREFIX='lei_ren_part1_paper_compliant_'
NAME=PREFIX+'current_core_first_interface.json'
RECEIPT=PREFIX+'current_core_first_interface_check.json'
GATE='current_core_bridge_functional_mixed4_join_certified'
ALL_GATE='all_current_bridge_functional_interfaces_certified'
SCOPES=('actual_point_moment_history_recovered','full_point_physical_field_evaluation',
    'full_inner_interfaces_certified','full_cartesian_vector_derivatives_certified',
    'global_completed_tensor_admissibility','admissible_stress_lift_constructed',
    'full_background_NS_validation','physical_energy_integral_certified',
    'independently_bounded_flat_remainder','temporal_recursion')
PRIOR_JOINS=('current_first_second_functional_mixed4_join_certified',
    'current_second_macro_functional_mixed4_join_certified',
    'current_R100_bridge_switch_functional_mixed4_join_certified')


def return_binding(stem,method,expression):
    fn=function(stem,method)
    rows=[node for node in fn.body if isinstance(node,ast.Return)]
    if len(rows)!=1 or ast.dump(rows[0].value)!=ast.dump(ast.parse(expression,mode='eval').body):
        raise ValueError('Complete production return formula changed: '+stem+'.'+method)


def defining_boundary_bindings():
    # All six atoms are actual common radial integrals. Comparison covers
    # from inner_bridge_profiles.inputs are deliberately not their source.
    return_binding('compliant_core_integral_atoms','finite_atom_coefficients',
        'dict(H=linear(phi,1,8),M=linear(uz,0,4),K=product(phi,uz,1,8),A=product(uz,uz,0,4),B=product(phi,phi,1,16),C=product(phi,phi,0,4))')
    binding('compliant_comparison_point_integrals','inlet','packet',
        'self.atoms.field.build_root_rows(24,6) if root else self.atoms.rebuild.rebuild(Z,24,6)')
    binding('compliant_comparison_point_integrals','inlet','atoms','self.atoms.atoms_from_packet(packet,6,shared_root=root)')
    binding('compliant_actual_bridge_integrals','prepare','data','self.comparison.inlet(Z,root)')
    binding('compliant_actual_bridge_integrals','packet','initial',"{n:j.truncate(5) for n,j in p['data']['moments'].items()}")
    binding('compliant_actual_bridge_integrals','packet','phi0',"p['data']['phi0'].truncate(5)")
    binding('compliant_actual_bridge_integrals','packet','pressure',"dress(own['actual']['C'],p['inputs']['F0_squared_ratios'])")
    binding('compliant_actual_bridge_integrals','packet','Q',
        "(2*z*V-z*own['actual']['M']*(1-self.core.delta)-d*derivative(own['actual']['M']))/L")
    binding('compliant_actual_bridge_integrals','packet','history_delta_V',
        "p['uniform_delta_V'] if not (chart=='first' and hi==0) else delta_V*0")
    fn=function('compliant_actual_bridge_integrals','packet')
    zero=ast.parse("if chart=='first' and hi==0:\n ell=ell*0;delta_phi=delta_phi*0;delta_V=delta_V*0").body[0]
    if not any(ast.dump(node)==ast.dump(zero) for node in ast.walk(fn)):
        raise ValueError('Exact first-inlet source increments no longer vanish')
    # Bind the core covering integral as an enclosure of the same primitive.
    # Its positive-radial-order derivatives are the differentiated square.
    binding('compliant_core_physical_field','profiles','covering','self.normalized_jets(c.mpf([0,endpoints(rho)[1]]),z)')
    binding('compliant_core_physical_field','profiles','K[gridkey(i,k)]',
        'rho*square_cover[gridkey(0,k)] if i==0 else square[gridkey(i-1,k)]')
    binding('compliant_core_physical_field','profiles','grids[PI][index]',
        'sum((math.comb(k,j)*source[\'relative_F0_squared_derivatives\'][j]*K[gridkey(i,k-j)] for j in range(k+1)),c.mpf(0))')
    binding('compliant_frozen_comparison_field','dress','c','jet.ctx')
    return_binding('compliant_frozen_comparison_field','dress',
        'jet*IntervalTaylor(c,[ratios[k]/math.factorial(k) for k in range(jet.order+1)])')
    return dict(all_six_true_atom_normalizations_AST_bound=True,
        original_shared_packet_comparison_inlet_AST_bound=True,
        actual_first_initial_values_and_moments_are_atoms_not_covers=True,
        exact_first_zero_source_increment_branch_bound=True,
        same_actual_pressure_C_and_radial_Q_recovery_AST_bound=True,
        original_covering_integral_is_only_an_enclosure=True,
        true_amplitude_product_derivatives_retained=True,
        pressure_identity='PI_core(4,Z)=4*dress(C); P_phys=Pstar^2*PD+epsilon_core*F0^2*PI_core=Pstar^2*PD+Ra*F0^2*dress(C)',
        axis_pressure_is_not_part_of_the_four_C_integral=True,passed=True)


def shared_packet_inlet(comparison,packet):
    """Replay original inlet arithmetic, replacing only its packet acquisition.

    A current common-core degree24/depth6 packet avoids a duplicate radial
    rebuild. No endpoint, coefficient, tail or moment arithmetic is changed.
    """
    if packet['radial_degree']!=24 or packet['axial_depth']!=6:
        raise ValueError('Unchanged original inlet degree24/depth6 required')
    original=function('compliant_comparison_point_integrals','inlet')
    fn=copy.deepcopy(original)
    rows=[node for node in ast.walk(fn) if isinstance(node,ast.Assign)
          and any(ast.unparse(t)=='packet' for t in node.targets)]
    if len(rows)!=1:raise ValueError('Unique original inlet packet acquisition required')
    old=copy.deepcopy(rows[0].value);rows[0].value=ast.Name(id='_common_packet',ctx=ast.Load())
    restored=copy.deepcopy(fn)
    next(node for node in ast.walk(restored) if isinstance(node,ast.Assign)
         and any(ast.unparse(t)=='packet' for t in node.targets)).value=old
    if ast.dump(restored)!=ast.dump(original):raise ValueError('Inlet replay changed derivative/atom/tail arithmetic')
    env=dict(comparison_module.__dict__);env['_common_packet']=packet
    exec(compile(ast.fix_missing_locations(ast.Module(body=[fn],type_ignores=[])),
        '<common packet original atom inlet>','exec'),env)
    return env['inlet'](comparison,packet['Z'],False)


class CurrentCoreFirstInterface:
    @source_precision
    def __init__(self,common=None,require_checked=True):
        self.common=common if common is not None else CurrentCoreCommonFixedPoint()
        if not self.common.acceptance_loaded:raise ValueError('Checked common nonlinear core required')
        self.source=self.common.source;self.bridge=self.source.bridge
        self.upstream=self.bridge.upstream;self.comparison=self.upstream.comparison;self.atoms=self.source.atoms
        self.core=self.source.core;self.ctx=self.core.ctx
        self.family=self.source.family;self.source_sha=self.source.source;self.datum_sha=self.source.datum_sha
        self.graph=dict(common_rebuild_is_actual_atom_rebuild=self.common.rebuild is self.atoms.rebuild is self.comparison.atoms.rebuild,
            common_original_core_is_actual_bridge_core=self.core is self.upstream.core is self.bridge.core.original,
            actual_comparison_uses_same_atoms=self.atoms is self.comparison.atoms,
            same_original_amplitude=self.source.amplitude is self.atoms.field.peak.amplitude,
            common_records_are_actual_current_records=self.core.records is self.bridge.core.records,
            common_source_domain_and_context=self.ctx is self.comparison.ctx is self.bridge.ctx,
            current_family_source_datum=(self.family,self.source_sha,self.datum_sha)==(self.bridge.family,self.bridge.source,self.bridge.datum_sha))
        if not all(self.graph.values()) or not self.bridge.acceptance_loaded:
            raise ValueError('Core/first must use one actual checked current source graph')
        self.hashes=dict(self.common.hashes)
        prior_name=PREFIX+'current_bridge_functional_joins_check.json'
        self.prior=accepted(prior_name,self.family,self.source_sha,PRIOR_JOINS[0]);_verify_hashes(self.prior)
        if self.prior['datum_enclosure_sha256']!=self.datum_sha or not all(self.prior[key] for key in PRIOR_JOINS):
            raise ValueError('Three same-source bridge joins required')
        self.hashes[prior_name]=sha(prior_name)
        for name in (Path(__file__).name,'lei_ren_part1_paper_core_recursion.py',
                     PREFIX+'core_integral_atoms.py',PREFIX+'comparison_point_integrals.py',
                     PREFIX+'bridge_mixed_C4.py',PREFIX+'inner_bridge_profiles.py',
                     PREFIX+'frozen_comparison_field.py',PREFIX+'microswitch_mixed_C4.py',
                     PREFIX+'K1_ledger.json'):
            self.hashes[name]=sha(name)
        self.bindings=defining_boundary_bindings();self.runtime={};self.acceptance_loaded=False
        if require_checked:
            report=accepted(RECEIPT,self.family,self.source_sha,GATE);_verify_hashes(report)
            if report['datum_enclosure_sha256']!=self.datum_sha or not report[ALL_GATE] or any(report[k] for k in SCOPES):
                raise ValueError('Fourth bridge interface receipt datum/scope differs')
            self.hashes.update(report['input_hashes']);self.hashes[RECEIPT]=sha(RECEIPT);self.acceptance_loaded=True

    @source_precision
    def boundary(self,Z='.359'):
        c=self.ctx;z=c.mpf(Z);key=z._mpi_
        if key not in self.runtime:
            packet=self.common.rebuild.rebuild(z,degree=24,depth=6)
            inlet=shared_packet_inlet(self.comparison,packet)
            first=self.bridge.evaluate(z,0,'first')
            raw=first['actual_parent_axial5_packet']
            p0=self.upstream.prepare(z)['inputs']['p0']
            ratios=self.upstream.prepare(z)['inputs']['F0_squared_ratios']
            C=inlet['moments']['C'];dressed=dress(C,ratios)
            pressure=dict(PD_normalized_axis_pressure=list(p0.coefficients),
                PI_core_normalized_increment=[4*v for v in dressed.coefficients],
                dressed_first_C=list(dressed.coefficients),
                pressure_units=dict(PD='P0/Pstar^2',PI_core='(P-P0)/(epsilon_core*F0base^2)',
                    first_C='(P-P0)/(Ra*F0base^2)'),
                prefactors_kept_in_logs=dict(axis=2*self.core.logP,increment=self.core.ctx.ln(self.core.epsilon),
                    F0_squared='2*(-selected_logCstar-Lambda*G(Z))'),axis_pressure_retained_separately=True)
            ledger=first['physical_logR_Z_mixed4_log_bound_ledger']
            groups={group:{name:len(grid) for name,grid in fields.items()} for group,fields in ledger.items()}
            self.runtime[key]=dict(Z=z,original_core_packet_degree=24,original_core_axial_depth=6,
                actual_core_atom_inlet=inlet['original_core_atoms'],common_pressure_boundary=pressure,
                actual_first_trace_axial5={name:raw[name] for name in (
                    'F_actual_over_F0_axial5_coefficients','Uz_actual_axial5_coefficients',
                    'actual_moment_shape_axial5_coefficients','actual_Q_axial4_coefficients',
                    'pressure_axis_axial5_coefficients','pressure_increment_true_axial5_divided_by_R_F0_squared')},
                original_first_phase0_controls=first['controls'],
                first_mixed4_component_row_counts=groups,
                source_radius_tree=first['source_radius_tree'],source_logR_over_Ra=first['source_logR_over_Ra'],
                same_packet_actual_comparison_acquisition=True,
                original_inlet_arithmetic_changed=False,
                **dict.fromkeys((GATE,ALL_GATE),self.acceptance_loaded),**dict.fromkeys(SCOPES,False))
        return self.runtime[key]


@source_precision
def run(field=None):
    field=field if field is not None else CurrentCoreFirstInterface(require_checked=False)
    result=dict(actual_five_defect_family_sha256=field.family,implicit_source_sha256=field.source_sha,
        datum_enclosure_sha256=field.datum_sha,current_common_core_bridge_graph=field.graph,
        original_boundary_source_bindings=field.bindings,
        fresh_original_first_boundary=field.boundary(),
        common_nonlinear_fixed_point_receipt=PREFIX+'current_core_common_fixed_point_check.json',
        retained_three_join_receipt=PREFIX+'current_bridge_functional_joins_check.json',
        **dict.fromkeys((GATE,ALL_GATE),False),**dict.fromkeys(SCOPES,False),input_hashes=field.hashes)
    (HERE/NAME).write_bytes((json.dumps(encode(pack(result)),indent=2)+'\n').encode('utf8'))
    print('Built common core -> original first inlet and pressure units; independent ODE/interface checker required',flush=True)
    return result


if __name__=='__main__':run()
