# Original pulse-end four support interfaces — 2026-10-04

All four original beta edges s=-3±.15 and s=-1±.15 now have explicit whole-Z flat difference bounds for full stress mixed3, velocity mixed4 and physical three-component error mixed2. The local reference removes only the beta input near that edge; it retains the exact common boundary moments, quadratic history, future energy, angular memory and analytic absolute pressure.

Run only this new layer:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage pulseendinterfaces
```

## Difference construction

Original beta derivatives through order four satisfy

`|partial_s^n beta(s-center)| <= A_n exp(-1/W) W^(-2n) / (ell^(n+1) N)`,

where ell=.15, N is the actual positive normalization and W<=2h/ell for distance h to an edge. The existing admitted analytic majorants and vanishing numerical caps are consumed; original beta functions and selected C5 axial coefficients remain unchanged.

For linear moment differences, the integrating-factor bound is h exp(lambda_i h) times the local forcing bound. Repeated ordinary-logR derivatives use delta_m' = delta_B - lambda_i delta_m. For the selected backward loss J, delta_J' = 2mu delta_J - delta(B^2). Constants in the unperturbed energy/pressure sources cancel only in the difference.

The current stress and physical operator require the full history in nonlinear terms. In particular:

- the axial nonlinear stress difference is B_actual times m_actual, since B_reference=0;
- the radial quadratic transport difference is delta_Ur times partial_y(Ur_actual) plus Ur_reference times partial_y(delta_Ur);
- the axial radial-error transport difference retains delta_Uz times the full axial operator on Ur_actual.

Replacing these products by delta_B delta_m or delta_Ur partial_y(delta_Ur) would discard generally nonzero reference histories. The implementation uses the full current end-chart C5 moment bounds. The original positive D factor is kept as its parent exact log recipe, rather than replacing it with the old cap.

Shared swirl, pressure, incoming angular memory and angular shear differences are exactly zero. This says the reference and actual functions agree in those components; neither actual stress nor actual moment histories are set to zero.

## Evidence

Focused producer/checker PASS:

- 415 current source hashes and 21 source/interface identities;
- 6,416 finite difference entries, 1,604 exact endpoint zeros and 4,812 monotone envelope comparisons;
- 1,448 independent original-beta quadrature/derivative comparisons: 800 full-stress, 360 velocity and 288 three-component-error rows;
- all four edges and both interior orientations, with eight nonzero reference radial histories;
- fixture comparison tolerance 1e-55; maximum positive enclosure miss zero.

The actual-source continuity argument is the analytic flat envelope plus source ODEs and current smooth physical pullbacks. Sample distances demonstrate the implementation of these bounds; they are not the proof of the whole flat limit. The accepted full-meridional Cartesian operator is consumed unchanged and is not rerun.

Stress mixed3 transfers to completed diagonal and divergence mixed2 through the current exact product operators. Velocity mixed4 transfers to physical error mixed2, including radial material and full viscosity terms. Normalized exact support coordinates identify endpoint beta zeros without treating a rounded numerical physical endpoint as an exact normalized r=±1.

## Scope

The bounds hold for the current source family, whole Z[-1,1], each original local support neighborhood and fixed positive tau. Positive B,D,H,R,lambda,nu factors are restored from current parent recipes. This closes the internal end-chart support interfaces.

These bounds do not provide an order-five radial Taylor remainder, uniform flatness as tau tends to zero, a whole-pulse cone, global completed-tensor admissibility, required-domain energy, actual n-dependent recursion or oscillatory cancellation.

The next task is continuous whole-end admissibility with the nonzero axial shear retained. In particular, the paper defines kappa=a+b^2/a, with a=2+2mu and b=2 partial_y(Uz)/Utheta. Inside a pulse, kappa-2 is 2mu+b^2/a; the pure-swirl simplification kappa-2=2mu cannot be used there.
