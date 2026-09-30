"""Exact dependency graph of the ACTUAL two quadrature moment equations.

Matrix entries are exact rationals of serialized canonical quadrature data.
Incoming/pulse rows stay named shared atoms: their exponents may have 28
digits, so converting their full values to enormous integers is forbidden.
This proves algebraic closure of the stored quadrature model, not continuous
integral closure of the pointwise velocity. No current tail is overwritten.
"""
from fractions import Fraction
import json
from pathlib import Path
import mpmath as mp
from lei_ren_part1_paper_axial_correction import from_signed_log, signed_log


class LinearAtoms:
    def __init__(self,terms=None):
        self.terms={k:Fraction(v) for k,v in (terms or {}).items() if v}

    def scale(self,value):
        return LinearAtoms({k:v*value for k,v in self.terms.items()})

    def __add__(self,other):
        result=dict(self.terms)
        for k,v in other.terms.items():result[k]=result.get(k,Fraction(0))+v
        return LinearAtoms(result)

    def evaluate(self,atoms):
        return sum((mp.mpf(v.numerator)/v.denominator)*atoms[k]
            for k,v in self.terms.items()) if self.terms else mp.mpf(0)

    def serial(self):
        return {k:dict(numerator=str(v.numerator),denominator=str(v.denominator))
            for k,v in self.terms.items()}

    def derivative(self):
        """Derivative for differentiable atoms and constant rational weights."""
        return LinearAtoms({k+'_Z':v for k,v in self.terms.items()})


class ExactQuadratureClosure:
    def __init__(self,axial):
        self.receipt=axial
        self.matrix=[[Fraction(x) for x in row] for row in axial['linear_matrix']]
        m=self.matrix
        determinant=m[0][0]*m[1][1]-m[0][1]*m[1][0]
        if not determinant:raise ArithmeticError('Actual stored matrix is singular')
        # Treat a_p(Z)*p_i as shared product atoms, so this is a functional
        # identity rather than freezing the energy root at a single Z.
        self.rhs=[LinearAtoms({'b1':1,'ap_p1':1}),LinearAtoms({'b2':1,'ap_p2':1})]
        rhs=self.rhs
        self.coefficients=[
            (rhs[0].scale(-m[1][1])+rhs[1].scale(m[0][1])).scale(1/determinant),
            (rhs[0].scale(m[1][0])+rhs[1].scale(-m[0][0])).scale(1/determinant)]
        self.rows=self.verify_matrix(self.matrix)
        if any(row.terms for row in self.rows):
            raise ArithmeticError('Exact quadrature row identity failed')
        self.rows_Z=[self.coefficients[0].derivative().scale(m[i][0])+
            self.coefficients[1].derivative().scale(m[i][1])+rhs[i].derivative()
            for i in range(2)]
        if any(row.terms for row in self.rows_Z):
            raise ArithmeticError('Differentiated shared-atom identity failed')

    def verify_matrix(self,matrix):
        return [self.coefficients[0].scale(matrix[i][0])+
            self.coefficients[1].scale(matrix[i][1])+self.rhs[i] for i in range(2)]

    def numeric_atoms(self):
        inputs=self.receipt['linear_rhs_inputs']
        a=mp.mpf(self.receipt['a_p'])
        return dict(b1=from_signed_log(inputs['base'][0]),
            b2=from_signed_log(inputs['base'][1]),
            ap_p1=a*mp.exp(mp.mpf(inputs['pulse'][0]['log_normalized_pulse_integral'])),
            ap_p2=a*mp.exp(mp.mpf(inputs['pulse'][1]['log_normalized_pulse_integral'])))

    def audit(self,precision):
        with mp.workdps(precision):
            atoms=self.numeric_atoms()
            c=[v.evaluate(atoms) for v in self.coefficients]
            old=[from_signed_log(row) for row in self.receipt['c']]
            drift=[mp.nstr(abs(x/y-1),40) if y else None for x,y in zip(c,old)]
            # Replay the rounded evaluated coefficients separately. Do not
            # confuse an exact dependency identity with this arithmetic sum.
            matrix=[[mp.mpf(x) for x in row] for row in self.receipt['linear_matrix']]
            rounded=[]
            for i in range(2):
                rhs=atoms['b'+str(i+1)]+atoms['ap_p'+str(i+1)]
                defect=sum(matrix[i][j]*c[j] for j in range(2))+rhs
                rounded.append(mp.nstr(abs(defect/rhs),40) if rhs else None)
            changed=[list(row) for row in self.matrix]
            changed[0][0]+=Fraction('1e-100')
            mismatch=self.verify_matrix(changed)
            if not mismatch[0].terms:
                raise ArithmeticError('Changed independent basis integral was silently zeroed')
            a=mp.mpf(self.receipt['a_p']); mu=mp.mpf(self.receipt['input_mu'])
            target=mp.mpf(self.receipt['energy_target'])
            energy=mp.mpf(self.receipt['K_p'])*a*a+mu*sum(
                mp.mpf(k)*value*value for k,value in zip(self.receipt['K_bump'],c))-target
            return dict(exact_row_terms=[row.serial() for row in self.rows],
                exact_differentiated_row_terms=[row.serial() for row in self.rows_Z],
                shared_atom_table=self.receipt['linear_rhs_inputs'],
                coefficient_expressions=[v.serial() for v in self.coefficients],
                numerical_coefficient_relative_drift=drift,
                rounded_coefficient_row_replay=rounded,
                original_row_replay=self.receipt['linear_relative_replay'],
                exact_stored_quadrature_row_identity=True,
                conditional_Z_derivative_row_identity=True,
                Z_derivative_identity_assumptions='Same Z-independent matrix and differentiable shared b_i(Z), a_p(Z)*p_i atoms; numerical/continuous Z derivatives are not certified.',
                independent_basis_mismatch_detected=True,
                independent_basis_mismatch_row_terms=mismatch[0].serial(),
                rounded_graph_energy_relative_replay=mp.nstr(abs(energy/target),40),
                exact_algebraic_energy_root_certified=False,
                continuous_integral_atoms_certified=False,
                installed_in_velocity_or_mean=False,
                physical_nonzero_tail_was_overwritten=False,
                finite_energy_certified=False,
                scope='Exact shared-atom algebra on measured quadrature inputs; actual continuous basis integrals and incoming primitive provenance remain open.')


def run():
    from lei_ren_part1_paper_joined_outer import build_joined_field
    from lei_ren_part1_paper_seeded_axial_inputs import solve_joined_seed
    print('collecting actual shared candidate matrix and inputs',flush=True)
    field=build_joined_field()
    solved=solve_joined_seed(field,'.3')
    graph=ExactQuadratureClosure(solved['axial'])
    report=graph.audit(field.precision)
    report['seed_origin']=solved['seed_origin']
    report['precision']=field.precision
    Path(__file__).with_suffix('.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    Path(__file__).with_name('lei_ren_part1_paper_seeded_shared_candidate.json').write_text(
        json.dumps(solved,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:report[k] for k in ('exact_row_terms','numerical_coefficient_relative_drift',
        'rounded_coefficient_row_replay','exact_stored_quadrature_row_identity','installed_in_velocity_or_mean')}),flush=True)
    return report


if __name__=='__main__':run()
