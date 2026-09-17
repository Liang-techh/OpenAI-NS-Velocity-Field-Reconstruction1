# Latest residual channel diagnosis

Read-only computation on continued_pressure using development seed 914027, 4096 points and validator step 0.005. These samples have informed development and are not final blind acceptance.

| Time | Maximum momentum residual | Azimuthal residual L2 |
| --- | ---: | ---: |
| .25 | 1.00275 | 2.329 |
| .3125 | .85313 | 2.004 |
| .4375 | .97319 | 2.124 |
| .5625 | 1.10107 | 2.304 |
| .6875 | 1.54433 | 2.430 |
| .75 | 1.94296 | 2.513 |

Azimuthal residual contributes 61–87% of squared momentum residual L2 over these times, including 82.1% at .75. Refinement from .01 to .005 changes the .75 maximum from 1.94577 to 1.94296. Axisymmetric pressure has zero azimuthal gradient; pressure-only optimization cannot repair this dominant channel.

Next CR003/CR005 implementation: introduce localized temporal swirl degrees of freedom (time slabs or cubic B-splines), preserving the divergence-free swirl form, initial field, core constraints, coefficient bounds and prescribed global torque. The current outer torque wrapper assumes a quadratic time polynomial for the base moment; replace that assumption with directly evaluated moment integration before using a new temporal basis. Compare against continued_pressure, retain fixed force and unchanged acceptance thresholds. Review external basis-capacity branch before implementing duplicate functionality; its evidence uses the older optimized_v4 candidate and is not current-candidate validation.
