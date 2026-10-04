# Exact heat equations and conditional exterior stress closure

Update: the actual terminal-history bridge and stress companion satisfy the theorem's exterior inputs. The inward collar pressure/main routing are also complete. The separate CompliantHeatPhysicalC4 now admits the actual regional physical momentum identity with correct viscosity units; read docs/HEAT_PHYSICAL_NS_2026_10_03.md. This theorem module remains generic. Collar/global stress and recursion remain open.

The canonical Gamma heat equations and terminal-normalization theorem are implemented in experiments/root_st073/lei_ren_part1_paper_compliant_heat_stress_equations.py. The source is Lei-Ren Part I v2, equations (3.12), (3.13), (3.18) and (5.1)-(5.11), PDF pages 20-21 and 58-60.

For a=delta/2 and the full positive Gamma expectation H:

```
xi^2*H'' + [1+2*(1+a)*xi]*H' + a*(1+a)*H = 0
ur=uz=0
utheta=A*r^(-1-delta)*H(4*nu*tau/r^2)
p=-integral_r^infinity utheta(rho,tau)^2/rho drho
```

Gamma integration by parts proves the ODE, including vanishing endpoint terms. It implies the angular heat PDE and radial pressure balance. In the original similarity stress equations, the angular source cancels against the viscous shear derivative and the axial source vanishes. The remaining stress forms are:

```
Ttheta=Ctheta(Z)/R
Tz=Cz(Z)/sqrt(R)
```

Matching the heat velocity alone does not determine Ctheta or Cz.

## Terminal-normalization theorem

With the original exterior zero-meridional moments and exact renormalized angular, energy and pressure terminal targets, H<=1 and |H'|<=a*(1+a) yield:

```
R*Ttheta = O(c*R^(-a)) -> 0
sqrt(R)*Tz = O(c^2*R^(-2*a)) -> 0
```

For 0<a<1/2, |Z|<=1 and R>=1, these limits force both homogeneous constants to vanish. The theorem and explicit bounds are proved conditionally on those exact terminal inputs.

Actual pressure-history transfer is admitted by the separate pressure bridges. Actual angular, selected-energy and meridional terminal-history transfer is now admitted by the terminal-history bridge. The exterior companion reports homogeneous_constants_eliminated_from_actual_terminal_moments=True and heat_exterior_stress_identity_certified=True. actual_compliant_field_physical_heat_region_certified remains False. Global stress, cone, flat remainder and recursion are also uncompleted.

## Evidence and next task

```
python experiments/root_st073/lei_ren_part1_paper_compliant_heat_stress_equations.py
```

Source identities and 12 independent full positive Gamma quadrature ODE checks pass, covering a=.003,.15,.35 and xi=0,.002,.08,2. The finite fixture supplements exact identities; it does not replace terminal-history proof.

Next extend equivalent pressure inward through the collar, adopt the companions in the main dispatcher, transfer the original physical map/prefactors and continue collar/all-region stress assembly.
