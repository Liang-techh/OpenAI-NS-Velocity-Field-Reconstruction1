# Current native pulse chart ownership — 2026-10-05

CurrentNativePulseSourceDispatcher consumes the checked current Rp source join and exposes all six original native pulse charts through one CompliantPulseMixedC4 object. Together with the fourteen accepted background charts, the current downstream registry has twenty source owners. The fifteen-owner Rp adapter and fourteen-chart physical adapter retain their earlier accepted scopes.

## Callable fields and domains

| Chart | Native method | Source coordinate | Original domain |
| --- | --- | --- | --- |
| pulse_entrance | entrance | t=log(R/Rp) | [0,.02/mu] |
| pulse_main | main | xi=mu*log(R/Rp) | [.02,10] |
| pulse_exit | main | xi=mu*log(R/Rp) | [10,11] |
| pulse_gap | gap | xi=mu*log(R/Rp) | [11,12] |
| pulse_gap_end | gap_from_end | s=log(R/Rv) | [-1/mu,-4] |
| pulse_end | end | s=log(R/Rv) | [-4,0] |

The dispatcher enforces these subdomains explicitly. The underlying main and gap methods accept broader ranges, so their own guards alone do not define the advertised chart domains.

Each chart returns four velocity/pressure component grids with all fifteen ordinary logR/Z mixed derivatives through total order four, plus the original primitive histories and formal amplitude metadata. The same live inlet, selected coefficients, pressure datum, source parameters, caches and pulse object are retained.

## Full interval coverage

The generated report contains 360 whole-chart derivative rows and sixty supplemental overlap rows, all on Z∈[-1,1]. Entrance coverage uses the upper enclosure of .02/mu.

For gap_end, the exact reciprocal endpoint cannot be passed as an independent outward interval to the original native guard without extending outside that guard. A legal real box starts at -lower(1/mu). Its source coordinate image xi=13+mu*s is enclosed below 12.0001. A same-object main-coordinate gap box [12,12.0001] covers the small remaining boundary interval. This establishes coverage of the original exact domain union while leaving the defining endpoint unchanged. The supplemental packet is coverage for the existing gap owner, not a twenty-first source owner.

## Acceptance and limits

The focused checker verifies current hashes, the twenty-owner registry, unchanged native bound methods and object identity, each whole coordinate box, finite mixed grids, gap overlap and narrower route guards. It consumes the checked external Rp receipt and replays the existing pulse_interface_certificate source identities for selected moment/energy ODEs, coordinate changes and axial derivatives through five.

The existing functional certificate is narrower than a uniform native two-sided mixed-four interface certificate. Therefore uniform_pulse_C4_chart_interface_certificate_available, full_pulse_C4_installed and quantitative_flat_velocity_interface_bound_ledger_available stay false. Having all six native charts callable does not certify global C4, cone admissibility, energy, nonlinear production points or temporal recursion.

Run `python experiments/root_st073/lei_ren_part1_paper_compliant_reconstruction.py --stage currentpulsechain`. The default CurrentNativePulseSourceDispatcher loads its own acceptance before certifying ownership. Use `evaluate(chart,Z,coordinate)`; `provider(chart)` returns the same native object for every pulse chart.

Next, explicitly extend the current physical Cartesian/time adapter through these owners and close quantitative native interface bounds as required. Shared full leading-source/remainder admission, production fields, completed global stress/flatness/required-domain energy, true n-dependent coefficient recursion and oscillatory corrected dynamics remain open.
