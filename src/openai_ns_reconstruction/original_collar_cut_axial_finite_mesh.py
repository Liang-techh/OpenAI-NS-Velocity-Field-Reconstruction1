"""Exact-real reduction of the archived third-row midpoint recipe."""
from __future__ import annotations

import ast
from fractions import Fraction

import sympy as sp


def _require(condition,text):
    if not condition:
        raise ValueError(text)


def _function(code,name):
    return ast.unparse(next(value for value in ast.walk(ast.parse(code))
                            if isinstance(value,ast.FunctionDef) and value.name==name))


def evaluate(*,inputs,pins):
    q2,whole,earlier = (inputs[k] for k in ('q2','whole','earlier_cut'))
    canonical = q2['actual_first_raw_n1_function']
    _require(q2['checks']['passed'] and canonical['passed'] and whole['passed'] and earlier['passed']
             and canonical['scale_conjugacy']['normalization']['u']=='Uz1/epsilon'
             and canonical['scale_conjugacy']['scale_matrix']['positional_map'][1]['inverse']=='1/epsilon'
             and q2['implementation_hashes'][pins['forward_path']]==pins['forward_code']
             and q2['implementation_hashes'][pins['cutoff_path']]==pins['cutoff_code'],
             'actual physical or finite normalization provenance changed')
    cell = _function(inputs['forward_code'],'source_cell')
    mesh = _function(inputs['forward_code'],'regular_integrals')
    seed = _function(inputs['forward_code'],'_prepare')
    cutoff = _function(inputs['cutoff_code'],'stable_sigma_point')
    uniform = _function(inputs['cutoff_code'],'_raw_uniform')
    inlet = _function(inputs['cutoff_code'],'_core_jets')
    _require("up = a['cp']['u'] * chip" in cell
             and "axialp = (a['point_frame']['U'][0].rows[0] * up * 8).truncate(1)" in cell
             and "chip = 1 - stable_sigma_point(c, 4 * mid - 2)" in cell
             and "mid = sum(_bounds(c, qbox)) / 2" in cell
             and "J.constant(c, 4, 5)" in seed
             and "cp, ct = self._core_jets(Z)" in uniform
             and "self.n1.evaluate(Z, '4', degree=20, kernel_degree=28)" in inlet
             and "'.5', '.75', '1', '2'" in mesh
             and "width = (hi - lo) / self.cells" in mesh and "point[j] += p[j] * wp" in mesh
             and "w = 1 / t ** 2 - 1 / (1 - t) ** 2" in cutoff
             and "return e / (1 + e)" in cutoff and "return 1 / (1 + e)" in cutoff,
             'old finite axial source or midpoint recipe changed')

    # The transition midpoint argument t=(2k+1)/32 occurs with its
    # reflection; no midpoint enters either numerical clipping branch.
    arguments = [Fraction(2*k+1,32) for k in range(16)]
    _require(min(arguments)>Fraction(1,1000) and max(arguments)<Fraction(999,1000)
             and all(a+b==1 for a,b in zip(arguments,reversed(arguments),strict=True)),
             'discrete cutoff reflection changed')
    E,t = sp.symbols('E t',positive=True)
    odds = 1/t**2-1/(1-t)**2
    _require(sp.cancel(odds+odds.subs(t,1-t))==0
             and sp.cancel(E/(1+E)+1/(1+E)-1)==0,'reflection is not exact')
    weight = Fraction(1,2)+Fraction(1,64)*8
    _require(weight==Fraction(5,8), 'old cutoff mass changed')

    eps,h,a,r,U,u = sp.symbols('epsilon hb a rho U Uz1_cut')
    physical = 2*eps*h*a*r*U*u
    _require(sp.cancel(physical-eps**2*h*a*(2*r*U*u/eps))==0,'row-three units changed')
    legacy = {}
    for Z,row in whole['legacy_forward_third_row_comparison'].items():
        _require(row['ledger_count']==212 and row['same_physical_unit']==whole['same_physical_unit'],
                 'legacy row left its actual ledger')
        legacy[Z] = row['ledger_sha256']
    return dict(passed=True, physical_to_legacy_scale_recipe='u_cut=Uz1cut/epsilon; ds=(sc/4)*dq',
                supersedes_earlier_reversed_scale_label=earlier['same_unit_recipe'],
                scaling_applied_in_whole_Z_calculation='Uz1/epsilon',
                mesh_cells_per_piece=16, transition_nodes=[str(a) for a in arguments],
                nominal_mass='5/8', nominal_finite_third_row='5*U_finite(4,Z)*u_finite(4,Z)',
                nominal_finite_ordinary_Z1='5*(U_finite_Z*u_finite+U_finite*u_finite_Z)',
                true_minus_nominal_parts=['actual_minus_5*Uaxis*C','5*(Uaxis-U_finite)*C',
                                          '5*U_finite*(C-u_finite)'],
                actual_or_binary_finite_inputs_equated=False,
                archival_floating_point_roundoff_certified=False,
                all_three_legacy_ledger_sha256=legacy, inherited_error_packets_removed=False,
                independent_terminal_matching_closed=False, ancestor_calculations_replayed=False)
