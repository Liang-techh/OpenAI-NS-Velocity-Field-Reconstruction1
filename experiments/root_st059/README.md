# ST059 — bulk-volume continuation from ST058-IP

Task #878. Original full NS1e-3 target remains UNMET. No default, publication, visualization, other-agent or schedule changes.

## Six actual full independent validations

B/E were frozen at2026-09-20T17:10:04.594773+00:00 before BOTH fresh seeds9205991/9205992. No later tuning or holdout-based coefficient selection. Each report uses4096Cartesianpoints, sixoriginaltimes and the FULL separate original spatial/time/energy-quadrature ladders. Norms below: fullvector Euclidean maximum and fixed-time spatialvolumeL2=sqrt(64mean|R|^2), worsttime at h=.005/timeh=.0025.

|Seed|Field|Full max|Spatial volume L2|
|---|---|---:|---:|
|9205991|ST058-IP|.029220349693429643|.04186240935141885|
|9205991|ST059-B|.029405093098592783|.041257954927351306|
|9205991|ST059-E|.02853430337661678|.03953542431106542|
|9205992|ST058-IP|.030474690716127684|.04194905366763824|
|9205992|ST059-B|.030537956834563448|.04155289982876589|
|9205992|ST059-E|.03047211027864723|.039743299779638995|

E L2 improves5.559%/5.258%; maxima improve2.348%/.00847%. Second peak is essentially unchanged, not a breakthrough. B improvesL2 1.444%/.944%, but WORSENS bothmaxima .632%/.208%. All6originalscienceCLIsexit1 for momentum_max andmomentum_L2; otheroriginalsampledgatespass. No continuum upperbound or uniformpointwiseimprovement.

NEW81x153axis-inclusivegrid at21times gives maxima parent.02996784339,B.03005434954,E.02996737426, all atinitialt=.25 near r=.28127,z=1.93855. Everygridpeak independently differentiated; worstFDvectorerror over63checks1.506e-7. Randompeaks can be larger. Refining the same randompeakchecks to h=.00125 gives parent .02922960099/.03051220349, E .02854346658/.03050966009. Tiny secondpeakadvantage survives refinement, not a supremumcertificate.

## Actual method and bounded execution

Same2594stored-coefficientcompactaxisymmetricfamily. Originalnu=.01,time[.25,.75],compactuANDp r<2/|z|<2,rawbounds,E0=1,restricted independently prescribedsolenoidalforce and0.001gates UNCHANGED. Forcefixed; initial spatialvelocityfree with exactoriginalnormalization and derivative.

B308reducedjointpoloidal/swirl/pressure directions; E578 after adding interior spatial directions. A numerical nullspace retains rawvelocity only at the historical singlecoreprobe at25times, not the fullinitialfield. This is not a solve with2592independentvariables.

Use .3timequadraturemean+.7smoothmax volume weights and separately screen squaredresidualintegrals in3absolute-z regions: bulk|z|<=1.75, adjacent1.75<|z|<=1.9, outer1.9<|z|<=2. Training aggregates bothzsigns at7times; independentaudits separate them. Regional ratio<=1.0000001 under actual nonlinear screening, denser262548pointpeakpool, originalpressure/core/shear/bias/profile/enstrophy/harmonic safeguards. Gauss-Newton finitehalfspaceQP directions are screened on the complete nonlinear field with backtracking; no continuousminimax/exactpressureelimination claim.

E uses192point spatialJacobianchunks, accumulating globalweakmoments before squaring. Dense/streamed assembly and directional derivatives calibrated. ActualB10steps69.31s; E16steps123.45s,15accepted afteraninitialrejection. Budgets differ, so not a pure matched-budget enrichmentablation. No optimizerconvergence claim. Legacy stats.loss is an inherited debugging objective; actualselection uses the registered weightedvolume/peak score andfeasibility.

## Independent regional result and counterevidence

Mapped40x32/64x48 quadrature per6signedslabs at5NEWinteriortimes+endpoints. At t=.75,E squaredresidualintegrals fall11.828%inbulk,11.412%inadjacentcollar,6.788%inoutercollar. These are squaredintegral reductions, NOT equal L2reductions. Bothbulkhalves improve atall7checkedtimes.

At othernewtimes, E two-sidedouterlayerregressesup to.6727%, one signedouterslab .8395%, adjacenttwo-sidedlayer .1118%. B's worstregionalregressions1.188%aggregated/1.274%signedslab. Training regionalnon-worsening did not generalize perfectly; noall-region/all-time claim.

## Structure and actual concentration

NEW800corepoints/17times: all5velocity/pressuredirection fractions1 for bothchildren. Core scope scaledR[.035,.215],|Z|[.025,.215],notwholefield. NEW1025midplaneprobes:Eminimumsignedshearretention.9996547828. Positivebias remains, but minimumuz falls .0001696702to.0001664573. Same-grid finalscaleddrift parent19.02462%toE18.84778%, a modest improvement,not source-matchpercentage; Bslightlyworsens19.03818%.

EinitialrelativevelocityL2change1.0648%, independent64x96E0=.9999999999994867. NonuniformL2changeafterbestglobalcomponentrescaling.01483659, so not onlyamplitudescaling. Prescribedcylinder16pi unchanged. Across3times energy participationvolume falls.497%-.714%; enstrophy participationvolume falls1.266%-1.635%. Smallconcentrationchanges, not exactvolumepreservation or an observedcollapse shortcut.

Degree2-8 harmonicL2lowerboundESTIMATES at t=.25/.5/.75 change .00549608/.00419555/.00533516 toE .00539326/.00399395/.00512368. Still>.001. Exactpolynomial/weakidentities plusfloatingvelocityquadrature, NOTintervalcertifiedlowerbounds.

## Pressure-only control and remaining dynamics

Separatelyregistered864pressurevariable EP fromE found NOacceptable nonzero step afterbacktracking,3.96s. Its coefficientvectorisEXACTLYE; it is not a thirdcandidate and doesnotproveoptimality.

IndependentweightedQR over108existing spatialpressuregradientdirections at each of3times dropscoefficient/core/timecouplingrestrictions ONLYforadiagnostic. WithEvelocity/forcefixed, the remainingL2 is approximately .037039/.029447/.036471 at t=.25/.5/.75. Removed fractions13.48%/16.31%/15.16%refer toSQUAREDresidual,notL2. Mosterrorisnotremoved bythosepressuredirections. No newcandidate, fullinfinite-dimensionalprojectioncertificate, intervalbound or family-wideimpossibilityresult.

## Actual evidence and remote scope

26new+inheritedfocusedlocaltests passed22.80s,warnings-as-errors; clean self-containeddelivery26passed22.24s. Sixfullvalidations and3completeauditbundlesfinished. PackagedEreplayactuallyexecuted fullfirstseedCartesianvalidation+newstructure+harmonicaudits in71.26s, scientificexit1, allspatial-laddermax/L2andgates reproducedEXACTLY(maxdifference0). Nooptimizerduringreplay. No cloudCI/nativeMATLAB/Lean/fullhistoricalsuiterun.

Remote branch currently contains the originalvalidatorwrapper (69bfb730f637c9794ec64c25f3a9462fac0e8563), executed mappedregion/pressureaudit (3ed7d03a4a8c7f540a1457a46278135eff5a1066) and thisresultsrecord. Completecoupledoptimizer/dependencies/rawarrays/maps/registrations/histories/failures/sixreports/allaudits are delivered in the conversationZIP and reviewablepatch, not allcommittedhere. The diagnosticimports therefore require that complete package/dependency integration. No PR/merge/defaultpromotion claim.

RawSHA256:
parentIP5ea5d117d7330f34b0b32f7dd67828c72efd0e03ac007bf9ca3115ce9f21123f;
Bfe2157923d7bda6b85261c4efbe31afbee214a718a149e9fa5d24e2ed69c7bea;
E422afc55fa7027ab5c76e268c4f383c41aaed51c6ff6ec4b7f22bd0aeb04efd0.

From complete package root:
```bash
python -m pip install -r requirements.txt
python verify_delivery.py
python experiments/root_st059/replay_st059.py --id ST059-E --out outputs/E --seed 9205991 --validate --structure --harmonic
```
Lastcommandcurrentlyexits1. Add --regions --pressure-budget for supplementarydiagnostics. Allscientificidentityandblowupflagsremainfalse.
