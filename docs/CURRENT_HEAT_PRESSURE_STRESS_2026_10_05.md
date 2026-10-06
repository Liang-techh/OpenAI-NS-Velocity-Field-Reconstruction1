# Current heat pressure and actual stress with retained constants

Implementation and scoped receipts: commit [fc383401](https://github.com/Liang-techh/OpenAI-NS-Velocity-Field-Reconstruction1/commit/fc383401).

F57B now has a current-source pressure/stress companion. It acquires the checked common core and all four bridge joins, then uses that graph's waiting, collar and full Gamma exterior. The producer and independent checker pass; checked loading and a fresh current-source coordinate are checked separately. This is partial F57B completion: current exterior stress is recovered, but has not been proved zero.

Run the focused stage:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentheatstress
```

Files: `experiments/root_st073/lei_ren_part1_paper_compliant_current_heat_pressure_stress.py`, `_check.py`, and their matching JSON artifacts. The controller uses one companion object for production and checking, reusing accepted upstream receipts by hash.

## The actual two terminal constants

Let t=log(R/Rtail), a=delta/2, k=1-a, b=(1-delta)/2, p=1+delta, d=1-Z^2 and L=1-delta*Z^2. All are the original current parameters. The actual angular memory is

```text
Dtheta(Z) = (1-epsilon) Xtail_current(Z) - A_future_current(0,Z).
```

The current waiting source supplies Xtail and Ptail. The original full collar pressure primitive is

```text
A_p(t,Z) = integral_t^infinity exp(-p*v) K(v,Z)^2/2 dv
Cp(Z) = Ptail(Z) + pressure_scale*A_p(0,Z)
      = pressure3(Z) + pressure_scale*exp(-3p)*Gamma_pressure_numerator(3,Z)
P_current(t,Z)/Pstar^2 = Cp(Z) - pressure_scale*A_p(t,Z).
```

Here pressure_scale=Ev0^2*theta_base^2/Pstar^2. Equality of the two definitions of Cp follows from additivity of the same source integral, including the infinite Gamma tail. It does not follow from interval overlap. Both actual enclosures are retained. Neither Cp nor Dtheta is set to zero or selected from an interval midpoint.

The current checker explicitly binds both actual Cp callables, the collar's pressure-defect density and its split at t=3. It consumes the checked current owner's exact Gamma pressure-rescaling theorem, proves equality of the full remaining density with the forward K-squared density, and identifies the two constants as functions of Z. Different directed boxes are retained; they need not be equal as intervals.

The exact inverse radius is S=1/Rtail and the full Gamma argument is xi=2*(1-Z^2)*S*exp(-t). Source radius/amplitude definitions and positive factors remain in logarithmic form. The numerical S enclosure and amplitude cap boxes bound these exact source quantities; cap endpoints do not replace them.

## Stress recovery in original units

The full projected future moments give the original collar stresses from (3.16)-(3.18). Their canonical Gamma contribution cancels in the exterior. The actual current constants contribute additional terms that must be retained:

```text
Qtheta = sqrt(R/2)*B, B=Ev0*theta_base*exp(-(1/2+a)*t)
Qz = sqrt(R/2)*B^2
Qpressure = sqrt(R/2)*Pstar^2

Delta Ttheta/Qtheta = exp(-k*t)*(k*Dtheta-b*Z*Dtheta_Z)/L
Delta Tz/Qpressure = (2*p*Z*Cp-d*Cp_Z)/L.
```

The angular extra moment is radially constant. Including physical factors, its stress has radial power -1. The pressure-offset stress has radial power +1/2. Thus a finite factored coefficient bound is not a uniform bound on the full pressure-offset stress over R to infinity. This term is an explicit remaining obstruction to a current stress-free exterior and must be closed through the original pressure/source equations. We do not infer that a nonzero interval enclosure proves the exact constant nonzero.

The separate Qpressure factor avoids dividing by the Ev2 enclosure, which contains zero. Ordinary mixed derivatives include the derivatives of the stated physical factors. This recovery is a leading similarity-stress calculation; transfer into a complete physical tensor/NS residual remains a later obligation.

## Evidence and precise admission

Eight current graph checks pass. Four views cover the full collar, its Gamma join, the entire [3,infinity) exterior and a fresh Z=.381 exterior source. They contain 60 absolute-pressure mixed4 rows and 150 factored stress rows. Collar stress is mixed3: its viscous shear consumes K derivatives through radial order four. Exterior actual stress is mixed4. The checker independently differentiates a nonzero-Cp/nonzero-Dtheta fixture: 15 full future pressure rows and 30 direct physical stress rows pass. It also proves 30 physical-factor identities for the actual constants, three exact radius/Gamma identities, the full collar moment-stress equations and canonical Gamma cancellation identities.

Admitted gates are:

- `current_heat_forward_absolute_pressure_mixed4_companion_certified`
- `current_collar_actual_stress_mixed3_recovered`
- `current_exterior_actual_stress_mixed4_recovered`

Current terminal-constant elimination and `heat_exterior_stress_identity_certified` remain false. Global tensor/cone admissibility, stress lift, flat remainder, required-domain physical energy, complete nonlinear physical point evaluation, full Cartesian residual and temporal recursion remain false. The unlocalized whole-space Gamma field still has infinite kinetic energy; this receipt does not change the required-domain energy obligation.

## Next work: connect current repair inputs and extend collar derivatives

The current native pulse supplies Xp and the flatten endpoint Xv by a production callable. The repair path stores a legacy covering enclosure of Xv. Existing datum/hash compatibility is not yet a proof of equality of those functions. Furthermore, the fifth angular coefficients and complete future energy acquire separate repair constructors. Their source functions and weights must be identified, rather than equating their boxes.

The minimum current angular bridge must bind:

```text
Xv_current = Xv_repair
scale*r_pre_repair(Z) = Xf_current(0)-Xf_current(Z)
scale*r_heat_repair(Z) = S_exact*Theta_current(Z)
W_A_future_repair = W_A_angular4_repair
```

Then replay the actual XR -> XS -> XQ -> XT -> Xtail propagation and prove

```text
(1-epsilon)*Xtail_current = 1/k + epsilon*J_current + S_exact*Theta_current.
```

If the initial source functions differ, rebuild the current repair from the actual native inlet before attempting the identity. Do not overwrite Dtheta afterward. Separately connect current Ptail and the full collar/Gamma pressure future to the original common axis datum and transported pressure primitive; only that equation may eliminate Cp. Extend original shape/Gamma derivatives to K_y5 before admitting collar stress mixed4, including the flat phi crossing. The detailed queue is at the top of AGENT_TASKS.md.
