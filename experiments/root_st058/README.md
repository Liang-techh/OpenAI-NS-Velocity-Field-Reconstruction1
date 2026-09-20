# ST058 — boundary-peak and volume-residual continuation

Task #867. Bounded numerical work from the complete frozen ST057-KS field. No default, visualization, publication, other-agent or schedule changes. Full NS 1e-3 remains UNMET.

## Six completed original independent validations

All chosen coefficients frozen at **2026-09-20T15:58:51.161236+00:00**, BEFORE BOTH fresh seeds. No subsequent fitting or selection from holdouts. Each seed4096Cartesianpoints, sixoriginaltimes, full separate original spatial/time/energy-quadrature ladders. Fullvector maximum and fixed-time spatialvolumeL2=sqrt(64mean|R|^2), worsttime at h=.005/timeh=.0025:

|Seed|Field|Full sampled max|Spatial volume L2|
|---|---|---:|---:|
|9205891|ST057-KS|.03980868206267804|.043225993603371426|
|9205891|ST058-S|.03435744643509721|.04276018060872198|
|9205891|ST058-IP|.028201833609731157|.04318486089140872|
|9205892|ST057-KS|.03895738235161998|.0436033413759781|
|9205892|ST058-S|.034071210267925665|.04253040118564458|
|9205892|ST058-IP|.028483452563298792|.04295015559880649|

S maximum improves13.69%/12.54%, L2 improves1.078%/2.461%. IP maximum improves29.16%/26.89%, L2 improves only0.095%/1.498%. S has lowerL2 on bothseeds, IP lowermaxima. Both improve both reported measures relative to parent on these samples, NOT at everypoint/time. All six actual scienceCLIs exit1 for momentum_max and momentum_L2; other original gates pass on these samples. No continuum bound or global-best claim.

## What changed in the solver

Same2594stored-parameter compact axisymmetric family. 308reduced coupled poloidal/swirl/pressure directions, not2592independent variables this round. Initial spatial distribution remains free, exact original E0 normalization and derivative retained, original prescribed force fixed. A numerical nullspace preserves raw velocity at only the historical single coreprobe at25times, not the full initialfield.

New dense trainingpool5562spacepoints x27times=150174, with improved sampling of both axial collars. Activepeak constraints feed local Gauss-Newton half-space QP steps (SciPy dual L-BFGS-B), followed by full nonlinear screening and backtracking. I reserves3%slack in linearized per-timeMSE and harmonic targets to allow useful nonlinear steps without crossing the final feasibility screens. No claim of a converged SOCP/continuous minimax or exact pressure elimination.

ActualS:10outeriterations,92.55s,harmonicweight2. I:warmS increment,12iterations,146.83s,harmonicweight8. IP:6fullpressure-onlysteps,864pressurevariables,8.40s,velocity/force exactlyI. I remains TRAINING-only intermediate, not separately fullvalidated. All registrations beforeexecution. Budgets, failed proposals and no-convergence status retained.

Originalnu=.01,t[.25,.75],smoothcompactu ANDp r<2/|z|<2,restricted solenoidalforce,rawbounds,E0=1,original max/L2 .001 gates unchanged. Autonomous finite screens preserve sampledcore/pressure,shear>=.995,bias>=.80,limiteddrift and enstrophymoments. Harmonic screening explicitly allows ratio1.0005: not strict monotonicity.

## Matched sampling diagnosis and independent peak refinement

SAME parent,old17times: oldspacepool sees.02948571, denser spacepool sees.03866781. Changing only time17to27 gives.03009668. The previous spatial collar sampling was missing peaks; this ablation is a design diagnostic, not independent acceptance.

NEW independent71x141spacegrid at19times: parentmax.03995965,S.03476979,IP.03060256. OriginalCartesianFD at every gridpeak agrees within1.65e-7. IP's largest gridpeak occurs atINITIALt=.25,r=.28562,z=1.93648 and exceeds its trainingpoolpeak.02825733. Training still misses some peaks; no continuous upperbound.

IPsecondholdoutpeak refines with h.005/.0025/.00125 to.02848345/.02849408/.02849477. The gain survives derivative refinement. Old full reports remain unchanged.

## Structure and actual volume checks

NEW800corepoints x17times: all five flow/pressure directionfractions1 for bothchildren. NEW1025midplaneprobes: signed shearretentionS.99905505/IP.99711039,above.995. Upwardbias remains positive, but minimumuz decreases parent.000163323 toS.000162528/IP.000151994. Finalexpandedprofiledrift slightly WORSENS .19234851 to.19244747/.19270092. Finitecore scopeR[.035,.215],|Z|[.025,.215],notwhole-domain/sourcecorrespondence.

IPinitialrelativevelocityL2change2.3203%; independent E0=.999999999999486. Bestuniformpoloidal/swirlscales leave nonuniformL2change.03279560: not justamplitude rescaling.

Declaredsupportcylindervolume remains16*pi. Separateparticipationvolumes use (integral density)^2/integral density^2. At t=.75,energydensityeffectivevolume15.15968 toIP15.26434(+.690%),enstrophyeffectivevolume13.04776to13.03951(-.063%). Across3times changes remain+.69%to+.81% and-.063%to-.330%,respectively. These measure concentration,not supportvolume/source-image similarity. No volume-collapse shortcut.

## Compatibility and remaining residual budget

96-order degree2-8 combined harmonic L2lower-bound ESTIMATES at t=.25/.5/.75:
parent .00570149/.00447788/.00554224;
IP .00549608/.00419555/.00533516.
Improved but still>.001. Exactweak/Gramidentities plusfloatingvelocityquadrature,NOTintervalcertifiedlowerbounds. Individualdegree2defects worsen. Sinitialdegree8slightlyworsens.00570297,within explicitauxiliaryslack. Newrandomtime and48/72/96quadraturechecks retained.

Disjoint residual-energy slabs at t=.75: IP |z|>1.9squaredresidualintegral falls about44%,but adjacent1.75<|z|<=1.9 worsensabout11%. Roughly85%of IPsquaredresidual remains in|z|<=1.75. This explains why a largepeakgain does not imply equallylargeL2gain. Separatecoarse/finequadrature,notcontinuumcertification.

## Executed tests and persistence boundary

18new+inheritedfocusedtests passed20.35s,warnings-as-errors. Fullobjective/constraint andpressure directionalcalibrations,coupledresidualderivatives,energy/support,immutablearrays andselectedsymbolicidentities checked. Sixfullvalidations and3completeauditbundles finished. One initialIfuture-import syntaxerror stoppedbeforefitting; failedtext/log andfixedexecutedsource retained. OptionalClarabelinstallation failedbeforeexecution; onlyinstalledSciPy used. No nativeMATLAB/cloudCI/Lean/fullhistoricalsuite wasrun.

This branch contains the original-validator wrapper (commit8b58669d9e8bb8c5ccfdab560a853daaa1c21262) and this finalresult record. Completeexecutedsource,rawfields,maps,histories,failures,sixreports andallaudits are in the conversationarchive/reviewablepatch. Do not infer that fulloptimizer orcandidatearrays have already been committed. No PR/merge/defaultpromotion isclaimed.

RawSHA256:
parentKS af06787a530f5152df01c7a361033bc939d5764498e2c2d04c9cb0fcb6a10237;
S e2c4001dfa3916cf3646891b4dbd93b35bac1872c63d5211f3bdeafc90850b5e;
I 40f336a0d322b7f0eb85574a8279216a59ebe32b93f644a7c496a0c6fe1aad96;
IP 5ea5d117d7330f34b0b32f7dd67828c72efd0e03ac007bf9ca3115ce9f21123f.

From complete delivery root:
```bash
python -m pip install -r requirements.txt
python verify_delivery.py
python -m pytest -q -W error experiments/root_st058/test_continuation.py experiments/root_st057/test_coupled.py
python experiments/root_st058/replay_st058.py --id ST058-IP --out outputs/IP --seed 9205891 --validate --structure --harmonic
```
Last command currently exits1 for failedmomentumgates. All scientificidentity/blowupflags remainfalse.