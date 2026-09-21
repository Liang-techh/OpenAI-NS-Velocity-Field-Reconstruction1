# ST063 — visible axial-core redistribution, with explicit tradeoffs

Task #939. Completed bounded local computation from ST061-P. No main/default/publication/schedule change. The original full NS1e-3 target remains UNMET.

## Actual geometry change

Two geometry controls were fitted inside the EXISTING2594stored-parameter compact axisymmetric family. No physical coordinate stretching and no new source-image numeric target. All1728velocity and864pressure variables were available; the original restricted force stayed fixed, with originalnu=.01,time[.25,.75],compact smoothu ANDp,rawbounds,E0=1 and original scientific gates unchanged.

An anisotropic angular-velocity target B used independent radial/axial scales, localized toward the axis and extending along z. It was fitted with existing basis columns, not introduced as a new Gaussian evaluated in the final field. The autonomous local swirl anchor was expanded from about20%to80%, while radial/axial anchors remained. The geometry-control screening explicitly allowed parent-relative sampled max/L2 up to3%higher; this is NOT a change to the scientific .001 thresholds.

Original G1/G2 L-BFGS fits reached650/900iterations (81.13/114.40seconds) and failed at least one structure screen. Their raw endpoints were not promoted. Separate preregistered336-direction poloidal/pressure restorations held raw swirl coefficients fixed and produced G1R/G2R in1/2restorationsteps. Energy normalization can still scale physical swirl. No optimizer convergence/global optimality claim.

All selected arrays frozen **2026-09-21T02:31:57.340710+00:00**, before either new seed or independent morphology sample. No holdout-driven retuning.

## Six complete original independent validations

Each4096Cartesianpoints,sixoriginaltimes,complete original spatial/time/energy-quadrature ladders. Report full-vector max and fixed-time spatialvolumeL2=sqrt(64mean|R|²),worsttime at h=.005,timeh=.0025.

|Seed|Field|Full sampled max|Spatial volume L2|
|---|---|---:|---:|
|9216391|ST061-P|.026279255060335002|.034035600656015505|
|9216391|ST063-G1R|.024546533709074312|.033767328775136926|
|9216391|ST063-G2R|.02438535329193401|.03425749818472021|
|9216392|ST061-P|.02566467314396125|.03358854497849373|
|9216392|ST063-G1R|.026861392722101148|.033883206585681475|
|9216392|ST063-G2R|.02250355909296729|.034206087967209475|

G2Rmax improves7.207%/12.317%, but L2 **WORSENS0.652%/1.839%**. G1Rsecondmax worsens4.663% and secondL2 worsens.877%. All6original scientificCLIs actuallyexit1 for momentum_max andmomentum_L2; otheroriginalsampledgatespass. No universal-best/defaultpromotion.

Independent81x133spacegrid/11times gives parentmax.02624635 andG2R.02385574. Each sampledpeak checked with originalCartesianFD at3steps; fine differences are numerical derivative checks, not PDE error bounds. No continuum maximum or new global location optimizer was proved.

## Fixed-window morphology and its limits

At64x96cylindricalquadrature in r<=.35,|z|<=.60, omega_z²-weighted centered axial variance divided by SINGLE-transverse variance r²/2 gives:

|Time|Parent aspect|G2R aspect|Gain|
|---|---:|---:|---:|
|.25|1.48755304|1.59680197|7.344%|
|.33713|1.55434662|1.65255503|6.318%|
|.5|1.65049839|1.73146643|4.906%|
|.64137|1.71165446|1.78030748|4.011%|
|.75|1.75122740|1.81256943|3.503%|

Two quadrature levels and a larger observation cylinder are retained. These are window-dependent moments, not intrinsic core boundaries or OpenAI-image similarity scores. The hoped-for10–20%aspect change is NOT achieved. At t=.5 the transverse RMS shrinks only1.20%,axialRMS grows3.64%.

The most visible change is the ACTUAL axial rotation profile: on-axis B at z=0 changes .09172868->.08589759 (center weakens), at z=.4 .04687709->.08446240, at z=.6 .02876447->.08350277. A longer moderate-rotation plateau, not stronger rotation everywhere.

At t=.5 and fixed absolute omega_z>=.15, the near-axis sampled band r=0,.025,.05 is connected across z[-.9,.9] for G2R versus parentcentralspan.488. Its reported1.8span is OBSERVATION-LIMITED, not the physical core length. At thresholds.25 and.35 neither field includes the axis midplane. General connected components can route around weak axis gaps; therefore component span alone is not treated as proof of a strong axial spine. Three-dimensional streamlines still look broadly similar and the middle disk remains.

## Structure and effective volume

G2Rpassesall5directionchecks on800newcorepoints×17times.1025newmidplaneprobes min signed shearretention.99698925; positivebiasminimumweakens .0001481113->.0001303017. Finalexpandedscaledprofiledrift improves .18657752->.18417910. Allfinite probes only. G1R has a radial-pressure sign miss (worstfraction.99875),which is not hidden.

Independent64x96energy participationvolume changes forG2R at.25/.5/.75:+.0518%,+.0373%,+.0173%; enstrophyvolumes:-.1239%,-.5281%,-.6369%. InitialE0=.9999999999996014,initialvelocityL2change1.5408%. Original support unchanged. No amplitude/volume-collapse shortcut, but concentration is not exactly unchanged.

## Visualization and execution boundary

`visualization/matlab/ns_compare_core.m` is a new native-MATLAB comparison UI: identical200seeds,physicalaxes,camera,time,absoluteisovalues and colorlimits; parentleft/childright; continuous originaltimepolynomial; bottomtimebar and actualon-axisprofile. Defaultomega_zlevel.25; predeclaredlevels.15/.35 are accessible. No cosmeticaxisstretch. Slice/whole-support modes are included. Streamlines are instantaneous,notparticles.

**This new MATLAB UI has NOT been executed natively.** The included st063_models.mat holds allthreeactualmodelcoefficientarrays; PythonMATroundtrip and numerical reference tests passed. MAT SHA2561f9402c4205b0548e789e64fb0c754d3f1c78830e3720b36167731c83d2f87c4. PNGpreviews are executed Pythonrenders,NOT MATLABscreenshots. The lightweight MATLAB ZIP includes evaluator,data,viewer and reference-test script; noPython/networkneededfornormalviewing.

Main14focusedtests passed12.59s,warnings-as-errors. Cleancompletepackage14passed12.70s,then G2R's full originalfirstseedvalidation AND independentgeometryaudit actually reran in58.77s. Scienceexit1 retained; entirevalidationreportdict matched exactly. No optimizer duringreplay. No cloudCI/nativeMATLAB/Lean/fullhistoricaltests.

Failedpermission-limitedplotstartup and initialunsimplifiedsymbolic-expression test are retained. Afterfixes, plotting andunchangedmathematicaltests passed; nofield/thresholdchanges. Fitcheckpoints saveactualarrays every10steps but do NOT serialize exact L-BFGS internal restart state.

## Persistent delivery

This branch currently stores the executed geometryaudit, newMATLABviewer and thisresultrecord. The COMPLETE actualoptimizer/repair code, dependencyruntime, rawparent/G1R/G2R, allregistrations/histories/discardedendpoints, sixfullreports, morphology/structure/volume/denseaudits, previews andMATdata are delivered in the conversationarchive. Do NOT claim allarraydata or dependencies are already committed. NoPR/merge/defaultchange.

Original raw hashes:
- Parent: d56e01b25dc0e6077687fe0b72c4608f4c332ae042230fb7e5a796c5ff189429
- G1R: 5482fc02f1b85d9f49a9d14c41cc4b5f25d191ad9fa0964eaf6773e9dc58b64c
- G2R: b9505453786af0ae0b0fd6e6d784d14c6892827b72cf2a70ace2f9f87acd13a0

Use COMPLETE offlinebundle: MATLAB `start_here`; Python `python verify_delivery.py`, `python -m pytest -q -W error`, then `python experiments/root_st063/replay_st063.py --id ST063-G2R --out outputs/G2R --seed 9216391 --validate --geometry`. Finalscienceexit1expected.

All PDE/sourceidentity/paperexact/blowupflagsremainfalse.
