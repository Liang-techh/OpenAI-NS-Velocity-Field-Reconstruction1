# ST066 — evolving log-scale profiles, not an accepted NS recursion

Task #1235. Completed bounded local continuation from actual ST065-C. No default, main, publication, other-agent or schedule changes. The original physical 1e-3 momentum target remains UNMET.

## Implemented representation

Static 9x12 F/G/P spatial profiles are extended with four Chebyshev terms in k/3-1: 1296 coefficients instead of324. The full physical time derivative includes dprofile/dk divided by tau*log(2), in addition to all coordinate/amplitude chain-rule terms. Initial-energy normalization is differentiated across all log-time coefficients and applied only once at k0. No per-frame renormalization. Prescribed q/h scales remain unchanged and are NOT dynamically recovered.

Original nu=.01,T=1,tau0=.5,k[0,6],h=.005,compact velocity AND pressure,independently prescribed force and absolute scientific thresholds retained. Full-domain audits include residual=-force outside shrinking velocity support. No exact annulus/exterior matching is implemented.

## Actual fits and independent phases

S is a static refit control,1000iterations/27.55s. D has four log terms,1600iterations/45.34s. B warms D with tighter energy retention and a quartic residual proxy,1400iterations/69.56s. All hit iteration limits; no convergence claim. Three endpoints frozen at2026-09-22T16:49:23.478255+00:00 before independent seeds9226691/9226692.

After phase-A diagnostics, W was separately registered: warm B, a35%terminal core-rotation target and axial sign/derivative penalties,1600iterations/84.53s. It was frozen at2026-09-22T16:55:38.407768+00:00 before NEW independent seeds9226695/9226696. This is sequential redesign with fresh holdouts, not untouched phase-A evidence. Tube B positivity is NOT fully achieved; its soft training guard retains a negative minimum. Targets are autonomous, not source-paper values.

Phase A:4fields x12scales=48 audits. Phase B:3fields x12scales=36 additional audits. Each has7boundaries plus5independent interior log-times,40/64/96full-domain quadrature,scale-resolved grid peaks and fresh core probes. Independent spatial/time finite-difference ladders cover parent/B in phase A and all3fields in phase B; S/D do not inherit those checks. Not the old ST063 4096point protocol or continuum certification.

## Physical residual tradeoff

|k|ST065-C volume L2|ST066-B volume L2|ST066-W volume L2|
|---|---:|---:|---:|
|0|4.01276822|3.48882874|3.73329586|
|3|19.26371018|17.97191750|21.15653766|
|6|92.84198336|80.64718981|105.07109295|

B L2 drops13.06%/6.71%/13.13%, but grid peaks change +.233%/-2.050%/+4.922%. W late-scale L2 worsens9.83%/13.17%; k6 peak1309.53031 versus parent1168.91297. Neither field is promoted. Conditioned training losses are not physical residuals. Absolute.001 gates fail badly.

## Actual winding and geometry

W characteristic diffusion-time turns at k0..6:
1.15106358,1.22455458,1.25486816,1.26232916,1.26758979,1.28972973,1.34930308 (+17.22%). This is N_nu,not actual completed particle rotations.

W re-seeded scale-aligned ensembles have mean interval turns .17768733,.18678672,.19544736,.20406709,.21304599,.22278635 (+25.38%). These are DIFFERENT parcels each interval. The SAME12initial particles reach only.38709808mean signed turns versus parent.32933506 (+17.54%). All12leave the shrinking observation core by k4; their per-interval means become smaller, not larger. W total absolute mean angular travel.38927767 exceeds signedmean,so some individual reversals remain. A special3particlemidplane control reaches1.33037turns and staysinside,notrepresentativeofwholevolume. Physical/similarity ODE and tolerance checks agree; no manufactured line-count evidence.

W fixed-window aspect at k0/k6:2.04457928/2.16579160 (+5.93%). Parent2.11625014/2.16071705: newinitialaspect is lower,andfinalonlyslightlyhigher. Do not claim5.93%betterthanparentatalltimes. W normalizedcoreprofiledrift reaches19.410%,so measuredcollapsenolongeriszero. A fitted varyingprofile plusimposedscalesisnotspontaneousrecursion.

Fresh finitecoreflow/pressure directionfractions1 atall12checkedscales,parent/B/W. Wmeanradialentrainment~.963-.967;limitedscope. Wsame-time effective-volume ratiosrange.990614..1.018811,andenergy ratios.991489..1.009493,withinauxiliarybounds onthesechecks,notalltime. Bminimumvolume.9898787 andenergy.9899208 slightlyviolate.99lowerguards; these failures retained. Noamplitudeorper-frameunitenergyshortcut.

## Visualization and actual tests

Standalone MATLAB start_evolving defaults toW; call start_evolving('ST066-B') or ('ST065-C') for controls. One log-scale slider, physical/isotropic/similarity views, actual same-particle paths and complete current-profile evaluation. It no longer reuses one static profile for every k. Allfive model datasets included. UI and evaluator have NOT been executed in nativeMATLAB. PythonMATroundtrip,coefficient/referencechecks were executed. PNGs are numericalPythonrenders,notMATLABscreenshots.

17focusedtests passed4.65s. Clean copied directory:17passed4.43s,compileall0; actual Wk6fullintegral/core/gridreplay completed1.265s andscientificexit1; entirereport equals original exactly. Resumeverified/skippedin.664s,exit1,nooptimizer. All84audits raninprimarydirectory; cleancheckdoesnotrerunall84. Originalparent728filemanifestverified. No cloudCI,Lean,fullhistoricalsuite ornativeMATLABrun.

Rawarrays saved every25fitsteps; no full L-BFGS internal-state restart claim. Originalsource-at-fit andmetadata-onlypostfreezeconfigcorrection are retained withdiff; numericalevaluator/arrays/thresholdsunchanged.

## Persistence boundary

Complete code, immutable profiles, fits, source snapshots, all reports, trajectory data, MATLAB data and figures are in the conversation archive NS_ST066_Evolving_Profiles.zip. The attempted numerical-source upload was blocked by the tool. This branch currently records this documentation, NOT a fully runnable source/data import. No PR/merge/fullsourcepublication is claimed.

WprofileSHA256:2ea8ed2eb1db076faf0daa19bd85a4de1873dbd46d0756a60b72f8182e9e77a8.
BprofileSHA256:7b9a7f6679dd70c35686951b81348bc654833a0c1ebc540b3ba2fbceb5e2c01e.
ParentCprofileSHA256:8b4bf31eae7fb48455ab0dc196688b670521f5fa959be7f31d52dd4b27cba54c.

From completeofflinebundle:
```bash
python -m pip install -r requirements.txt
python verify_delivery.py
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python -m pytest -q -W error tests
python replay.py --candidate ST066-W --k 6 --out outputs/W_k6.json
```
Finalcommand intentionallyexits1;--resumevalidatesexistingreports,nooptimizer.

All PDE/dynamicrecursion/sourceidentity/blowupflagsfalse. W is awinding-orientednegative-PDEcontrol;B is alower-volume-residualcontrol. Neitherestablishestheuserscale-recursioncapability.
