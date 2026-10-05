# Actual finite-width bridge histories composed through R110

Date: 2026-10-05. This supersedes the open R100-to-R110 companion transfer in ACTUAL_BRIDGE_INTEGRALS_2026_10_04.md. It completes a source-functional enclosure interface, not selection of production point histories or the full corrected NS solution.

Run the bounded stage:

```powershell
python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage actualbridgeswitch
```

The producer/checker are lei_ren_part1_paper_compliant_actual_bridge_switch{,_check}.py and their matching JSON receipts. CompliantActualBridgeSwitch.packet(Z) accepts real Z intervals in [-1,1]. Reports include the complete interval and Z=0,.5.

## What is now connected

ActualR100BridgeAdapter consumes upstream.packet(Z,1,"macro"). Exact source geometry simplifies 2hb+(log(100/Ra)-2hb) to log(100/Ra). Its named interface is R=100; the adapter rejects other actual-field radii.

The actual phi, raw axial V, six H/M/K/A/B/C own histories, canonical P0 and amplitude derivative ratios are passed to the unchanged original inner-switch inputs/phase/post/inlet methods. H/M/K/A/B/C map to Mtheta, Mz, Mtheta_z, axial quadratic, swirl quadratic and pressure shapes. C is a pressure primitive, never raw axial V. Original dress(C,F0_squared_ratios) and canonical P0 are retained.

The known comparison field remains the prescribed shear source. Its finite-width macro endpoint is consumed with exact_endpoint_at_R100, then frozen comparison histories continue with theta=100/R: H/K/B relax at rate2 and M/A/C at rate1, with targets phi,V,phi*V,V^2,phi^2/2,phi^2. Seven AST bindings guard the endpoint call, geometry flag, target/transport expressions and radius bounds. Actual own moments are not used to redefine Dbar/Ebar.

Original microscopic controls and hb=cstar*K^-100 remain unchanged. The original R2=100*exp(2hb) and exact a=.8,b=0 postpower carry all actual histories through R110. The output supplies axial5 velocity, six moment shapes, canonical pressure and actual R110 log shape, plus axial4 radial Q. The inherited signed-switch callable now uses the same new actual inlet and finite-width comparison data.

## Signed-integral corrections

Older macro, first-switch and complete-switch receipts now retain the full F0(Z)^2 logarithmic range [-2logC-2Lambda*Gbar,-2logC]. Gbar is a bound, not G(Z). Factored log packets explicitly describe enclosures of exact amplitudes. Complete-switch packets also export the separate amplitude log before adding the enormous width log.

Comparison inputs remain axial6 until the direction differentiates moments. Complete signed-switch directions therefore really contain six axial5 coefficients. The previous early truncate(5) lost the highest direction row. The checker now requires exact row counts and checks the separate full swirl-amplitude range.

Only affected signed receipts and the new actual-bridge dependency receipts were refreshed. Existing physical/cone/mixed4 source modules were not regenerated or relabeled.

## Evidence and limits

Focused composition PASS: 307 current hashes; 7 unchanged callable identities; 7 comparison-continuation AST bindings; 180 actual R100 velocity/moment/pressure coefficients transferred; 72 independent prescribed-direction coefficients; 720 finite actual endpoint coefficients. Six independent positive Volterra integrals check postpower moments against direct quadrature with tolerance 1e-70. All three older signed checks and the finite-width bridge check also pass.

Read-only GPT-5.6 Luna/max review accepted geometry, moment units, comparison rates and derivative orders. Its comparison-provenance concern was resolved by the seven explicit source bindings and endpoint guard.

The new companion sets R100_R110_actual_feedback_composed=True. Existing mixed4 providers are not yet replaced by this companion. Exact production point histories, downstream leading five-defect inputs/Jacobian/remainders, global completed-tensor admissibility, global temporal-flat/physical-volume/required-domain energy, actual n-dependent recursion and oscillatory correction remain open.

## Next work

1. Install the same source adapter into the existing switch mixed4 provider, including required rectangular derivatives and physical microscopic width factors. Preserve the original controls and endpoint function identities.
2. Transfer this actual R110 inlet through long reshape and restoration. Use its actual pressure and all six histories; distinguish raw Uz from V_pressure=4C.
3. Recompute the shared leading five-defect input, implicit Jacobian and nonlinear omitted-term bounds. Keep companion composition distinct from this acceptance gate.
4. Proceed to independent global tensor/flatness/energy gates, then the actual n=1 and n>=2 recovery equations and oscillatory correction. Coordinate scaling remains insufficient for recursion.
