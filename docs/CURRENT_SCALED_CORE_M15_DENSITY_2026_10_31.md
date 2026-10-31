# Scale-preserving core m15 density — 2026-10-31

The goal remains **ACTIVE / INCOMPLETE**. The new `src/openai_ns_reconstruction/original_core_scaled_m15_density_v1.py` is a local coefficient-space algebra for the signed core density, not an integral or a completed physical field. It does not append to the 212 inherited m15 terms, remove any of their errors, or establish an endpoint identity. P1 and m15 remain distinct.

## Implemented interface

At fixed rho, write `Delta = U_lead - U0` and let the supplied `F1` and `Uz1` already include the cutoff eta. Given ordinary Z0/Z1 input rows and `L=F0_Z/F0`, the leaf emits separate signed terms for

\[
 d_{\rm ax}=2\epsilon(U_0+\Delta)U_{z1},\qquad
 d_{\rm sw}=-2\epsilon^2F_0\rho\Phi F_1.
\]

It uses the ordinary product rule, in particular `F_lead_Z=F0*(Phi_Z+L*Phi)`. The Z1 input row of each supplied field is the derivative of the **full** input field, already expressed in its own formal powers; any derivative of an F0 factor inside a supplied response is the caller's responsibility. The leaf requires caller premises `rho_Z=epsilon_Z=eta_Z=0` and eta-included responses. These are caller claims, not independently certified geometry or source identities. Moving-endpoint derivatives are outside this interface.

Epsilon and F0 remain symbolic; explicit input powers are combined with the outer powers, including a synthetic check both with and without a factored F0 in `F1`. Each finite input and separately declared additive tail is expanded into labeled signed products; no tail is assumed absent without an explicit empty list, no source error is discarded, and multi-cell terms retain their cell IDs without summation. An original, **not Pstar-normalized**, N1p Z0/Z1 packet can be supplied explicitly as a separate `-epsilon^2*rho*N1p` pressure channel with per-term formal powers. A formal epsilon power `-1` can represent an already factored Lambda inside N1p without evaluating an inverse. The pressure channel is excluded by default and is never integrated or appended to the existing pressure/endpoint ledger.

This two-power epsilon/F0 ledger is not yet a verified conversion to the older four-source-power W1 basis; other source-scale metadata is retained as provenance, not silently folded into these powers. Each coefficient must cover the entire requested cell if used as an enclosure. The leaf neither certifies that coverage nor validates the backend's mathematical range semantics, binder identity, live source owner, or matching source/datum. It permits externally supplied certificates but never promotes their claims into its own proof. No midpoint or point-value claim is accepted.

## Check and evidence

The sole frozen entry ran 29 synthetic checks. A fresh rational-interval fixture enclosed independent synthetic coefficient witnesses `(axial Z0,Z1)=(8,28/15)` and `(swirl Z0,Z1)=(-5,-1/2)`, preserved tail products and source/scale labels, and exercised the optional separate pressure channel. It rejected missing derivative or tail declarations, eta-free responses, nonzero fixed-geometry derivatives, and missing pressure scale powers. These checks validate this algebra on synthetic inputs, not the original physical source or full-cell quadrature. The preceding historical sample at `(rho,Z)=(2,.371)` was retained only as provenance and was not used as a cell enclosure.

Frozen source SHA-256 `def2815851ef648733366134e5e93569fa7265a3e85dfd0f8c33c2d63de270b0`; entry SHA-256 `6aa147b829b453111bbe73a4f320b0e5152bec69573aacc507ee260291be4c41`; complete result SHA-256 `93a68a939535942858ee0a8aba2842205ea7fa937be21a48bf01ee0ef334deca`. The complete receipt is `experiments/root_st073/scaled_core_m15_density_v1_receipt.json.gz`, raw SHA-256 `eb0b3f0551ff2b860cfe26057741f53da41d4d925b018954435d9ece286e8a48`. Authenticate raw gzip bytes before decompression. Scope flags describe the new runner, not a process-wide trace.

## Next step

Construct a new, frozen source-backed adapter that supplies same-source full-rho-cell Phi, Delta, F1, Uz1 and N1p ordinary Z0/Z1 enclosures with all tails and the correct existing source powers. Link the output to the guarded integrator only after full-cell coverage, actual provider and binder semantics are established; maintain fixed-rho and moving-endpoint derivatives separately. The saved single-point oracle cannot supply this evidence. Integrate the result only after verifying which signed terms are already in the six inner/lower contributions and 212-item ledger, preserving all inherited errors and avoiding duplicate pressure or own-swirl terms. Functional five-moment terminal closure remains unproved.
