"""Automatic axial tangents and moment-derived radial velocity at the exit.

The inherited prescribed ODE is differentiated, including its drivers and
initial moments. This is the derivative of the finite RK construction, not
an enclosure of the exact ODE or a global Cartesian divergence certificate.
"""
from functools import lru_cache
import mpmath as mp
from lei_ren_part1_paper_pressure_width_bridge import PressureWidthExitBridge
from lei_ren_part1_paper_pressure_width_jet import PressureWidthJet
from lei_ren_part1_paper_mp_stress import evaluate_mp_stress


class AxialPressureWidthExitBridge(PressureWidthExitBridge):
    def jet(self, value=0):
        return self.comparison.jet(value)

    @lru_cache(maxsize=64)
    def evaluate(self, s, Z):
        with mp.workdps(self.precision):
            raw = super().evaluate(s, Z)
            z = mp.mpf(str(Z)); delta = self.comparison.delta
            result = dict(raw)
            for name in ('y', 'R', 'F', 'Uz', 'P', 'log_F_over_Fa',
                         'F_R_over_F', 'Uz_R', 'chi'):
                result[name] = raw[name].value
                result[name + '_Z'] = raw[name].tangent
            result['moments'] = {k: v.value for k, v in raw['moments'].items()}
            result['moments_Z'] = {k: v.tangent for k, v in raw['moments'].items()}
            R = result['R']; root = (2*R).sqrt()
            F = result['F']; FZ = result['F_Z']
            gy = R*result['F_R_over_F']; uy = R*result['Uz_R']
            result['Utheta'] = root*F
            result['Utheta_Z'] = root*FZ
            result['Utheta_y'] = root*F*(gy + mp.mpf('.5'))
            result['Uz_y'] = uy
            convert = lambda value: PressureWidthJet(value,
                pressure_order=self.pressure_order, width_order=self.width_order)
            stress = evaluate_mp_stress(mp.log(self.r), z, delta,
                Utheta=result['Utheta'], Uz=result['Uz'],
                Utheta_y=result['Utheta_y'], Utheta_Z=result['Utheta_Z'],
                Uz_y=uy, Uz_Z=result['Uz_Z'], moments=result['moments'],
                moments_Z=result['moments_Z'], P=result['P'], P_Z=result['P_Z'],
                precision=self.precision, radius_override=R, scalar_converter=convert,
                shear_theta=2*F*gy, shear_z=root*uy/R, include_components=True)
            result['Ur'] = stress['U_r']; result['stress'] = stress
            L = 1-delta*z*z; d = 1-z*z
            result['Ur_R'] = ((1+delta)*z*result['Uz'] + 2*z*uy
                             - d*result['Uz_Z'])/(L*root) - result['Ur']/(2*R)
            # This identity uses moment derivatives prescribed by the ODE.
            # It does not differentiate the finite RK interpolant in R.
            result['ODE_divergence_numerator'] = (
                L*(root*result['Ur_R'] + result['Ur']/root)
                + d*result['Uz_Z'] - (1+delta)*z*result['Uz'] - 2*z*uy)
            result.update(Z_moment_derivatives_available=True,
                radial_velocity_recovered_from_moments=True,
                independent_Cartesian_divergence_verified=False)
            return result

    def metadata(self):
        result = super().metadata()
        result.update(Z_tangents_installed=True,
            derivative_of_finite_RK_construction=True,
            radial_velocity_recovered_from_moments=True,
            independent_Cartesian_divergence_verified=False,
            axial_tangent_remainder_enclosed=False)
        return result

    def physical_chart(self, s, Z, logq, *, nu='.01', phi=0):
        """Cartesian component-valued field on the same local exit chart."""
        with mp.workdps(self.precision):
            z = mp.mpf(str(Z)); q = mp.exp(mp.mpf(str(logq)))
            viscosity = mp.mpf(str(nu)); angle = mp.mpf(str(phi))
            if viscosity <= 0:
                raise ValueError('Positive viscosity required')
            field = self.evaluate(s, Z); delta = self.comparison.delta
            radius = (2*viscosity*q*field['R']).sqrt()
            axial = mp.sqrt(viscosity)*q**((1-delta)/2)*z
            scale = mp.sqrt(viscosity)*q**(-(1+delta)/2)
            ur = mp.sqrt(viscosity/q)*field['Ur']
            ut = scale*field['Utheta']; uz = scale*field['Uz']
            return dict(xyz=(radius*mp.cos(angle), radius*mp.sin(angle), axial),
                uvw=(ur*mp.cos(angle)-ut*mp.sin(angle),
                     ur*mp.sin(angle)+ut*mp.cos(angle), uz),
                pressure=viscosity*q**(-1-delta)*field['P'], tau=q*(1-z*z),
                scope='Local finite exit chart; global matching and error bounds remain open.')
