> Successor: [actual original-N first-switch values](CURRENT_ORIGINAL_FIRST_SWITCH_NUMERIC_2026_10_15.md) now continue R100 through R1 with unequal original controls and retained histories. This macro document remains the accepted prerequisite; next work starts at second_switch.

# Actual original-N macro velocity, radius and own mean — 2026-10-15

`OriginalMacroNumeric.evaluate(Z,q,...)` continues actual modified leading+n1
physical velocity enclosures from the accepted second-micro endpoint R0 to
the exact original R100 endpoint. The native fraction is q in[0,1].

```python
from openai_ns_reconstruction.original_macro_numeric import OriginalMacroNumeric
field = OriginalMacroNumeric()
inlet = field.evaluate('.3', '0', log_tau='-2', azimuth='.25')
relaxation = field.evaluate('.3', '1e-18', log_tau='-2', azimuth='.25')
exit = field.evaluate('.3', '1', log_tau='-2', azimuth='.25')
```

## Original geometry and actual source

The defining source geometry remains

```
Y = log(100/Ra) = log(25)+logLambda
Lmacro = Y-2hb
R0 = Ra*exp(2hb)
R(q) = Ra*exp(2hb+Lmacro*q)
rho(q) = 4*exp(Y*q+2hb*(1-q))
R(1) = 100
```

Unlike the microscopic charts, the ordinary Y*q geometry must be retained.
The physical reader keeps its full logarithm in radius, radial/swirl
velocities, Cartesian u/v and streamfunction weights. Only the original
positive hb contribution is represented by a proved small error. No hb
bound becomes the source width or radius. Neither huge original F0 nor
the double-dyadic frequency is materialized.

Only the comparison profiles are frozen. The actual Phi/U source still
evolves under the original hb control and known comparison direction. The
numerical core proxies are enclosed using the COMPLETE defining bridge
integrals: first micro, second micro, and macro. They are not restarted at
R0, and no complete micro contribution is added twice. All hydro, original
pressure and swirl drive parts remain. The macro radius kernels already
include their integration in physical log radius; one hb source factor
is applied, without a second Lmacro or hb Jacobian.

## True own mean and retained corrections

With U=U0+epsilon*Psi and M=U0+epsilon*Mpsi, the normalized history satisfies

```
dMpsi/dy = Psi-Mpsi
theta = exp(-(Y-2hb)*q)
Mpsi(q) = theta*Mpsi_micro2
          + (1-theta)*Psi_core
          + integral_0^(Lmacro*q) exp(-(Lmacro*q-y))*(Psi(y)-Psi_core) dy
```

The incoming unmodified Mpsi uses the accepted actual finite-width second
micro history, not a comparison mean. The last integral uses an absolute
modal norm uniform over EVERY macro prefix. A signed endpoint integral
cannot bound earlier prefixes after cancellation. Epsilon normalization
and the original hb log are applied before the small error cap, avoiding
loss of a tiny physical source in enormous normalized units.

Phi/Psi have Taylor Z0..3, Mpsi Z0..2 and Qpsi Z0..1. The common Q recovery
uses this same actual mean. The incoming first-chart selected-N mean Z3
cap remains unresolved, so no Mpsi Z3 or Qpsi ZZ promotion occurs.

The canonical original selected integer remains
`N=2^(2^1358356628656378313)`. Signed Section11 functions retain their
original phase `frac(N*(2hb+(Y-2hb)*q-hb*sc/2))`. Full-phase caps enclose
the numerical proxies; no loop angle is chosen or resolved.

All five selected-N correction histories include both complete micro
windows and the traversed macro prefix. Micro density uses hb*ds; macro
density uses (Y-2hb)*dq, once each. Positive Duhamel kernels are bounded
above by1; Lmacro<=Y is used only in a conservative error bound. Original
rate-zero pressure memory remains in the history packet. These are
magnitude bounds of the defining signed histories, not evaluated signed
integral centers or a completed pressure reconstruction.

The n1 cutoff velocities f/u remain exactly zero here. Their terminal
primitive persists, and mean1=primitive/current_rho uses the true large
macro rho, rather than the micro rho4 proxy. q1 retains weight1+delta.

## Focused run

Run `python scripts/run_original_macro_numeric.py --write-report`, then
`python scripts/run_original_macro_numeric.py --inspect-report` from the
repository with the existing runtime dependencies. The saved report is
`experiments/root_st073/lei_ren_part1_paper_compliant_current_original_macro_numeric.json`.

The focused driver checks modal kernels and resonances against independent
quadrature; nonzero-inlet Volterra mean relaxation and axial derivatives;
actual geometry and the distinct macro Jacobian; the accepted second-micro
join; persistent correction/n1 histories; and cross-precision coefficients
and log geometry. Small q=1e-18 resolves the macroscopic mean relaxation
near its inlet; q=.5 and q=1 cover the large radius expansion.

The accepted run exits0 in190.676s with13 original packets,60 independent
modal/Volterra/mean/geometry/measure checks,204 selected-N small-error proofs,
63 second-micro/macro native and physical coefficient overlaps,84 incoming
history/mean-relaxation/terminal-memory checks,15 current-rho and endpoint
checks, and26 precision60/85 coefficient/log comparisons. Seven unsupported
requests are rejected. Saved-report inspection exits0 in2.288s. The report
pins1751 source inputs and32 implementation files.

At Z=.3, the numerical Mpsi/Psi ratio changes from0.5 at q=0 to about
0.8049868424165964 at q=1e-18, and is1 within the reported error at q=1.
This resolves actual own-mean relaxation across the macro inlet; it is not
an imposed zero or a frozen comparison mean. The enormous normalized Psi
units retain their original epsilon physical factor.

## Remaining work

Tasks14/15 stay IN PROGRESS; the full goal stays ACTIVE / INCOMPLETE.
Next continue first_switch from its ACTUAL R100 Phi/U, own histories and
all five selected-N correction memories. Use the existing true prescribed
downstream/high-radial sources and original sigma control. Do not substitute
the earlier candidate N=257, frozen comparison data, or zero histories.

Also open: total radial/fast-phase derivatives and physical vorticity beyond
r_minus; tighter incoming mean Z3; later chart evaluation and complete
signed five-defect/selected moment repair/own pressure; full Cartesian C2
joins; global pressure and common N admission; exact heat exterior,
finite energy and radial tail; admissible stress and flat remainder;
genuine temporal n>=2 recursion and smooth sum; mean/two-family oscillatory
cancellation; and independent fixed-forcing corrected Cartesian NS tests.
Native axial Taylor derivatives and a radial chart extension do not prove
temporal coefficient recursion or a full Navier–Stokes blowup solution.
