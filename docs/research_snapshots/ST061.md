# ST061 — alternative constrained methods with volume retention

Task #900. Actual bounded continuation from complete ST060-Q. Original full NS1e-3 target remains UNMET. No default, visualization, publication, other-agent task or schedule changes.

## Six completed independent original validations

All selected arrays frozen2026-09-20T20:46:14.692395+00:00 BEFORE fresh9206291/9206292. No later tuning or parameter selection. Each4096Cartesianpoints,sixoriginaltimes,FULL original separate spatial/time/energy-quadrature ladders. Fullvector maximum and fixed-time spatialvolumeL2=sqrt(64mean|R|^2),worsttime at h=.005/timeh=.0025:

|Seed|Field|Sampled maximum|Spatial volume L2|
|---|---|---:|---:|
|9206291|ST060-Q|.027559315226219006|.03374187106068435|
|9206291|ST061-D|.026788020110924435|.03333465903335943|
|9206291|ST061-P|.02580782220423295|.03347432641256566|
|9206292|ST060-Q|.027992393733693038|.03457784543286563|
|9206292|ST061-D|.027204469154934973|.03426968736612646|
|9206292|ST061-P|.02620795221578385|.034440441141806846|

Pmax improves6.355%/6.375%,L2 only.793%/.397%. Dmax improves2.799%/2.815%,L2 1.207%/.891%. D is lower-L2,P lower-peak. All6scienceCLIs actuallyexit1 for momentum_max/momentum_L2; otheroriginalsampledgatespass. No uniform/all-time/continuum acceptance claim. C is retained as TRAINING-only methodcontrol,not independentlyvalidated thisround.

## New implemented method and research alternatives

Same2594stored-parameter compactaxisymmetricfamily,737existingcoupleddirections; initialspatialvelocityfree with exactdifferentiatedE0normalization; prescribedrestrictedforcefixed. Originalnu=.01,t[.25,.75],compactu ANDp,r<2/|z|<2,rawbounds andoriginalscientificgates unchanged. Inheritedcore/pressure/shear/profile/morphology/harmonic/regionguardsretained. Energy/enstrophy participationvolume ratios[.99,1.02] remainfiniteauxiliaryconstraints,notnewPDEthresholds.

422800trainingpoints plus180declaredcollarsentinels informedbypriorpeaks,notfreshholdouts. C:oldquadraticsupportingcuts,6steps/2accepted,125.81s. D:exact quadraticenergies of LINEARIZED residuals in a96D adaptiveHessian-metricsubspace; installedSciPySLSQP withlazylinearconstraint exchange thenfullnonlinearscreening,6steps/4accepted,108.20s. Neitherconverged. D's trainingpeakisbetter,C's trainingL2slightlybetter. Onefinitebudgetcase,notuniversalspeedup.

Pwarm-startsDstep6 beforeholdouts. Parent-normalized peak+L2sum mustdecrease; eachstepmayraiseonenormbyatmost.5%,but allparent-basednonlinearsafeguards remain.10extrasteps/9accepted,210.09s. Finaltrainingmax.02625890/L2.03416480. Additionalwork/warmstart prevents pureone-factorattribution. Not a formalSQP-filter implementation orconvergenceproof.

Research read primaryClarabel docs andGoulart-Chen2024, SciPyleast_squares docs, andEspanol-Jeronimo2024 variableprojection. Clarabel installation failedtwice beforeexecution; noClarabel/CVXPY/OSQP/SCSsolveclaimed. Actual full864pressure unconstrainedLSelimination decreasesweighteddiscretesquaredresidualenergy14.44% andworsttrainingL2.03436501->.03202885; normal-equationrelativeerror3.39e-16. Butalltestedcoreaxialpressuredirectionsflipwrong (correctfraction1->0),so proposalREJECTED. Velocity/forceunchanged,Rthetaunchanged. Notconstrainedvariableprojection orcontinuumlowerbound. Completeauthorreferences andexperimentboundary areintheoffline METHOD_REVIEW.md.

## Independent peaks, structure and effective volume

New91x171spacegrid/20times:parentmax.02866041308 ->P.02654506941. Allgridpeaks independentlyCartesianFDchecked,largestvectorerror1.36e-7. Same11starts forboundedlocalpeaksearches find parent.02924590046 att=.25,r.28058,z1.93583;P.02654659407 att=.75,r.75241,z1.92762 (~9.23%lower). Parent8/11/P10/11localsearchesreportsuccess; stopsretained. FineFDvectorerrors<=8.65e-9. Neithergridnorlocalsearchisaglobalupperbound.

Fresh800corepoints/17times:all5flow/pressuredirections1forbothchildren. Fresh1025midplaneprobes:Pminimumsignedshearretention.99725730. PositivebiasminimumWEAKENS .000192734->.000173438. FinalexpandedprofiledriftslightlyWORSENS .18601103->.18630760. ScopeR[.035,.215],|Z|[.025,.215],notwhole-domain/sourceimagecorrespondence.

Independent64x96quadrature at.25/.5/.75:Penergy participationvolume changes[-.4774%,-.3269%,-.1739%],enstrophy[-.7459%,-.8123%,-.8874%]. BothD/Pstaywithin[.99,1.02]on3checks,notalltime. PinitialrelativevelocityL2redistribution1.5425%,E0=.9999999999995216; supportcylinder16pi unchanged. Smallconcentrationchanges,notliteralunchangedvolume oramplitudecollapse.

All25signedslab/time finequadrature residual-energyintegrals decreaseforP. Att=.75 bulk decreases.785%,adjacentcollars5.799%,outer6.313%;squaredintegrals,NOTL2percentages orpointwiseproof. Harmonicdegree2-8necessarylower-boundESTIMATES .00466022/.00305253/.00496685 ->.00456628/.00299216/.00482061 at.25/.5/.75;still>.001,notintervalcertified. Complete6-typeextraaudits onlyparent/P;Dhasoriginalfullvalidator pluscore/volumeaudits.

## Actual checkpoints, tests and evidence boundaries

Dreallystoppedaftergeneration2,thennewprocessresumedtogeneration6. C/D0-6,P0-10states containrawcandidate,delta,map,parent,sourcebinding,damping/history. Frozenfitsrejectfurtherfitting. C'soriginalsource-at-execution preservedbecausealaterscalingfix changedonlyitsUNUSEDdirectsolvermodule. Cutnumericsnotchangedorrelabeled.

Primary27focusedtests passed17.99s,warnings-as-errors. Independentcleanpackage27passed16.83s;compileandall-generationverificationexit0. CleanpackageP actuallyreplayedfullfirstseedoriginalvalidationin65.41s andreturnedscienceexit1; entirereportdict/space,time,energy laddersmatchEXACTLY. --resume took.55s,verified/skippedreport,exit1,nooptimizer. Supplementaryaudits actuallyraninprimarydirectory andcopiedbyte-identically;notreruninpackagesmoke.

InitialtoySLSQPstatusfailurefixedbyvariable/objectivescaling withoutlooseningtests. Auditsinitiallyfailedbeforeexecutionduemissinginheritedvolume_geometryimportpath;failedsource/logsretained,pathfixed,sixcompletedreportsverified/skipped,thenallunstartedauditscompleted. No candidate/operator/seed/thresholdchanges. No cloudCI, nativeMATLAB, Lean orfullhistoricalsuite executed.

Threeactualoptimizers committedat2deed372ad5339e41222fecfdd2189a9b541b2cc;returnedGitblobidentitiesmatchlocalbytes. Completeexecutedsource/dependencies,rawQ/C/D/P,allmaps/generations,negativeprojectioncontrol,sixfullreports andallaudits areintheconversationZIP. Do notclaimfullarrays/maps/dependenciesalreadyonGitHub. No PR/merge/defaultpromotion.

From COMPLETE offline delivery:
```bash
python -m pip install -r requirements.txt
python verify_delivery.py
python experiments/root_st061/session_status.py
python -m pytest -q -W error
python experiments/root_st061/replay_st061.py --id ST061-P --out outputs/P --seed 9206291 --validate --structure
```
Scientificexit1expected. Add--resume forverifiedcompletedreports;nooptimizer. Dstructure-replayscopeiscore/volumeonly.

RawSHA256:
D `728bd92b3cf321c46726c953ec72342cba12464a77cc0e2c8b1cd157d54b5b37`;
P `d56e01b25dc0e6077687fe0b72c4608f4c332ae042230fb7e5a796c5ff189429`.

Allpde/sourcecorrespondence/paperexact/blowupflagsremainfalse.