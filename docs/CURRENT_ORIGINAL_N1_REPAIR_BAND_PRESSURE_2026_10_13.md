# Actual repair-band pressure moments, 2026-10-13

The original compact band Rc to 2Rc now has evaluated signed pressure contributions over all four original axial cells [-1,-0.5], [-0.5,0], [0,0.5] and [0.5,1]. This adds a required radial segment to the n1 pressure defects. The full reconstruction goal remains ACTIVE / INCOMPLETE.

## Completed

- [x] Bind the accepted current C3 repair-band velocity, complete own mean and finite-N corrections to one common physical Pstar scale. Reuse the accepted quiet coefficient k=K/Rc from its genuine original source kernel.
- [x] Restore the absolute original Rc from the accepted native frame. The phase offset alone is not an absolute radius. The source identity is log Rc = original phase offset + log r_minus, equivalently log(110)+10*selected_logCstar+11*exp(40)+124.
- [x] Apply the correct x=R/Rc measures: dR=Rc dx and R dR=R*Rc dx. No extra R is inserted in the unweighted or weighted integrals.
- [x] Retain the finite boundary terms from q_R. After normalizing m13 by Rc and m15 by Rc squared, their error bounds keep the required factor 1/Rc.
- [x] Use U U_Z in the radial-free identity. Keep actual tiny nonzero delta in directed arithmetic, including an underflow regression at delta=1e-1000. Explicit delta<1/2 gives the inverse-L bound used by the conservative boundary constants.
- [x] Preserve the actual finite-N velocity and own-mean corrections as nonzero magnitude errors. The compact baseline is q=-k/x with zero baseline U. Its exact pressure contributions are m13/Rc=k^2/8 and m15/Rc^2=-k^2*log(2)/4; these do not replace the corrected source.
- [x] Evaluate both corrected contributions on all four full axial cells. m13 is positive and m15 negative; each source-normalized coefficient enclosure has relative width about 2.8031565691e-4, below 0.001.
- [x] Protect public construction, evaluation and reports with the same 540-digit global precision context. Default-global-precision construction and evaluation pass. Directed interval endpoints are reboxed from their exact binary values rather than rounded decimal strings.
- [x] Check the nonlinear physical pressure identity, nonunit radius/Jacobian/endpoint scaling, derivative selection, tiny delta and signed error expansion in every axial cell. The public evaluated query passes.

The quartet uses the `lei_ren_part1_paper_compliant_current_original_n1_repair_band_pressure_source` stem under `experiments/root_st073`. Its raw candidate gate remains false; the accepted receipt enables `current_original_n1_repair_band_pressure_source_installed`. The independent checker and actual run passed. A standalone CLI assembles the same checked constructors and new run/check calls: `python scripts/run_n1_repair_band_pressure_source.py`. The CLI wrapper was compiled; the worker ran the equivalent production/check/public sequence directly.

The numeric certificate concerns values of these two segment integrals. Signed defining function handles retain their ordinary Z1 structure, but precise pressure derivative bounds and higher derivative accuracy are not established by this value certificate. Absolute radii and N are not materialized. The logarithmic coefficient enclosures retain the original Rc and Rc-squared physical scales for later composition.

## Next construction tasks

- [ ] Refine the actual outer pulse saddle integrals before claiming precise outer pressure values. A current source query at Z=0.371 showed about 5.862% width in scaled end coefficients, so reusing that enclosure would limit the next pressure calculation.
- [ ] Preserve correlation between the original two end coefficients when consuming refined input intervals. Their actual selected branch and moment/energy equations must remain unchanged; no representative coefficient may replace the source function.
- [ ] Evaluate remaining active pre-repair chart contributions and selected Rp-to-Rv pulse contributions. Reuse accepted lower/O3/compact-band/quiet pieces exactly once.
- [ ] Compose complete m13/m15 defects with correct physical factors and all finite endpoints, together with m11/m12/m14.
- [ ] Bind positive I2/Astar and normalized bump matrices, perform n1's original 2+3 repair and recover its own M/V/P.
- [ ] Verify functional five-moment closure, axis and derivative joins, and pressure compatibility before genuine n2/temporal recursion.
- [ ] Complete global admissible stress, oscillatory stress cancellation, flat remainder and independently checked corrected NS residuals.

Complete pressure tails, n1 second/third repairs, outer support installation and global physical pressure defect gates in this adapter remain false. This result does not certify the later temporal/stress/corrected-NS stages.
