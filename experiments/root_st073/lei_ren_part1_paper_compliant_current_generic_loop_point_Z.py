"""Approximate Section 11 A/B and phase-held slow-Z point backend.

Inputs are explicit point functions and their jets, never saved range caps.
The existing scalar loop supplies phase inversion. Analytic chain rules and
angle quadrature supply the missing slow-Z primitive rows. This backend
does not select native parameters, certify quadrature, or install a global
original source oracle.
"""
from dataclasses import dataclass
import json
from pathlib import Path
import time
import lei_ren_part1_paper_compliant_current_generic_shear_loop as loop

HERE,PREFIX,sha=loop.HERE,loop.PREFIX,loop.sha
NAME=PREFIX+'current_generic_loop_point_Z.json'
RECEIPT=PREFIX+'current_generic_loop_point_Z_check.json'
GATE='original_Section11_scalar_phase_held_A_B_Z_backend_implemented'


@dataclass(frozen=True)
class PrimitivePointZ:
    phase: object
    angle: object
    angle_Z: object
    A: object
    B: object
    A_Z_slow: object
    B_Z_slow: object
    approximate_only: bool=True
    certified_original_point_oracle: bool=False


class GenericLoopPointZ:
    """Fixed eta/d_star; derivatives hold the fast phase and y constant.

    Native source dispatch must separately supply true a,b,p1,p2,E and
    ordinary Z derivatives. These values are neither whole-chart caps nor
    derivatives of a chart selector. p1 enters the accepted cone guard but
    not the defining primitive formulas. B inherits E's units: supply
    E=Utheta/Pstar to obtain the graph's B_over_Pstar row.
    """
    def __init__(self,scales,*,a,b,p1,p2,E,a_Z,b_Z,p2_Z,E_Z):
        self.base=loop.GenericShearLoop(scales,a=a,b=b,p1=p1,p2=p2,Utheta=E)
        self.ctx=c=scales.ctx
        self.a_Z,self.b_Z,self.p2_Z,self.E_Z=(c.mpf(v) for v in (a_Z,b_Z,p2_Z,E_Z))
        if not all(c.isfinite(v) for v in (self.a_Z,self.b_Z,self.p2_Z,self.E_Z)):
            raise ValueError('Finite true slow-Z point derivatives required')
        v=self.base;eta=scales.eta
        self.t0_Z=-self.b_Z/v.a+v.b*self.a_Z/v.a**2
        self.kappa_Z=self.a_Z+2*v.b*self.b_Z/v.a-v.b**2*self.a_Z/v.a**2
        self.flat=v.kappa>=2+eta
        if not self.flat and v.q==0:
            raise ArithmeticError('Active q rounded to zero; increase input precision')
        if self.flat:
            self.q_Z=c.mpf(0)
        else:
            x=(2+eta-v.kappa)/eta
            step=loop.flat_step(c,x)
            if 0<x<1:
                odds=1/(1-x)**2-1/x**2;e=c.exp(-abs(odds))
                slope=e/(1+e)**2*(2/x**3+2/(1-x)**3)
            else:slope=c.mpf(0)
            root=c.sqrt((2+2*eta-v.kappa)/(2*v.a))
            root_Z=-root*(self.kappa_Z/(2+2*eta-v.kappa)+self.a_Z/v.a)/2
            self.q_Z=-slope*self.kappa_Z/eta*root+step*root_Z
        self.u_Z=(self.p2_Z*v.q+v.p2*self.q_Z)/scales.d_star
        self.h_Z=v.u*self.u_Z/v.h
        self.r_Z=self.u_Z/v.h**3
        self.C=v.q/v.h
        self.C_Z=self.q_Z/v.h-v.q*self.h_Z/v.h**2
        self.nu=1+v.t0**2+2*v.q**2
        self.nu_Z=2*v.t0*self.t0_Z+4*v.q*self.q_Z

    def direction_Z(self,psi):
        c=self.ctx;v=self.base;psi=c.mpf(psi)
        if not c.isfinite(psi) or not 0<=psi<=2*c.pi:
            raise ValueError('Angle must lie in [0,2pi]')
        if self.flat:return self.t0_Z
        cosine=c.cos(psi);den=v._denominator(psi)
        # Stable correlated denominator; no subtraction of rounded units.
        numerator=(v.one_minus_abs_r-2*c.sin(psi/2)**2 if v.r>=0 else
                   2*c.cos(psi/2)**2-v.one_minus_abs_r)
        w=numerator/den
        w_Z=self.r_Z*(-den-2*numerator*(v.r-cosine))/den**2
        return self.t0_Z+2*self.C_Z*w+2*self.C*w_Z

    def angle_integral_Z(self,psi):
        c=self.ctx;v=self.base;psi=c.mpf(psi)
        if not c.isfinite(psi) or not 0<=psi<=2*c.pi:
            raise ValueError('Angle must lie in [0,2pi]')
        if psi==0:return c.mpf(0),c.mpf(0)
        if psi==2*c.pi:return 2*c.pi*self.t0_Z,2*c.pi*self.nu_Z
        if self.flat:return psi*self.t0_Z,2*psi*v.t0*self.t0_Z
        points={c.mpf(0),psi}
        if 0<c.pi<psi:points.add(c.pi)
        # Resolve the Poisson peak near 0 (r>0) or pi (r<0).
        width=v.one_minus_abs_r/c.sqrt(abs(v.r)) if v.r else c.mpf(1)
        if width<c.mpf('.25'):
            center=c.mpf(0) if v.r>0 else c.pi
            for factor in (1,4,16,64):
                for sign in (-1,1):
                    point=center+sign*factor*width
                    if 0<point<psi:points.add(point)
        points=sorted(points)
        T1_Z=c.quad(self.direction_Z,points)
        T2_Z=c.quad(lambda s:2*v.direction(s)*self.direction_Z(s),points)
        return T1_Z,T2_Z

    def evaluate(self,phase):
        c=self.ctx;v=self.base;phase=c.mpf(phase)
        if not c.isfinite(phase):raise ValueError('Finite fast phase required')
        original_phase=phase;phase-=c.floor(phase)
        if phase==0 and original_phase>0:phase=c.mpf(1)
        psi=v.angle_at_phase(phase)
        if self.flat or phase in (0,1):
            zero=c.mpf(0)
            return PrimitivePointZ(phase,psi,zero,zero,zero,zero,zero)
        T1,T2=v.integrals(psi);T1_Z,T2_Z=self.angle_integral_Z(psi)
        actual_phase=(psi+T2)/(2*c.pi*self.nu)
        phase_Z=T2_Z/(2*c.pi*self.nu)-actual_phase*self.nu_Z/self.nu
        phase_psi=(1+v.direction(psi)**2)/(2*c.pi*self.nu)
        psi_Z=-phase_Z/phase_psi
        chi=phase-psi/(2*c.pi);M=-v.a*T1/(2*c.pi)-v.b*phase
        A=v.a*chi/2;B=v.Utheta*M/2
        A_Z=self.a_Z*chi/2-v.a*psi_Z/(4*c.pi)
        M_Z=-self.a_Z*T1/(2*c.pi)-v.a*(T1_Z+v.direction(psi)*psi_Z)/(2*c.pi)-self.b_Z*phase
        B_Z=self.E_Z*M/2+v.Utheta*M_Z/2
        return PrimitivePointZ(phase,psi,psi_Z,A,B,A_Z,B_Z)


def run():
    began=time.monotonic();checked=json.loads((HERE/loop.RECEIPT).read_bytes())
    if not checked['all_passed'] or not checked[loop.GATE]:raise ValueError('Accepted original scalar loop receipt required')
    scales=loop.GenericLoopScales(a_min='1',margin_min='2',boundary_kappa_excess_min='.02',
        t0_abs_max='1',p1_abs_max='8',p2_abs_max='1',dps=60)
    backend=GenericLoopPointZ(scales,a='1.5',b='.2',p1='6',p2='.05',E='2',
        a_Z='.03',b_Z='-.04',p2_Z='.02',E_Z='.1')
    q=backend.evaluate('.37')
    record=dict(**{GATE:True},source_family=checked['current_source_family'],
        input_hashes={Path(__file__).name:sha(Path(__file__).name),Path(loop.__file__).name:sha(Path(loop.__file__).name),loop.RECEIPT:sha(loop.RECEIPT)},
        exact_original_loop_value_and_phase_inverse_reused=True,
        analytic_q_and_direction_slow_Z_chain_rules=True,
        phase_inverse_Z_from_original_implicit_equation=True,
        original_eta_and_d_star_Z_independent=True,
        point_primitive_A_B_and_slow_Z_backend_implemented=True,
        separate_original_point_input_provider_required=True,
        quadrature_or_roundoff_certified=False,
        original_native_parameter_point_values_selected=False,
        numerical_original_source_point_or_integral_oracle_installed=False,
        actual_five_controls_installed=False,certified_actual_fixed_point_tail_installed=False,
        actual_terminal_Z_function_closure_installed=False,current_whole_N_selected=False,
        **dict.fromkeys(loop.OPEN,False),
        example_explicit_point_reference=dict(phase=q.phase,A=q.A,B=q.B,A_Z_slow=q.A_Z_slow,B_Z_slow=q.B_Z_slow,
            manufactured_inputs_only=True,approximate_only=True),
        execution_seconds=time.monotonic()-began,
        scope='Original Section11 recipe point backend for A/B and slow-Z rows with explicit true scalar inputs; approximate angle quadrature and implicit inversion. Native source/parameter dispatch and numerical certification remain open.')
    (HERE/NAME).write_text(json.dumps(loop.encoded(record),indent=2)+'\n',encoding='utf8')
    print('Original loop point A/B and phase-held slow-Z backend implemented',flush=True)
    return record


if __name__=='__main__':run()
