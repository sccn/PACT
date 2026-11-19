# PACT local validation for PR #6

## Environment and reproduction

The checks ran locally on Windows 11 build 26200 with MATLAB R2024b Update 3
(24.2.0.2806996), EEGLAB 2024.0, FIRfilt 2.8, MATLAB Signal Processing Toolbox
and Statistics and Machine Learning Toolbox. The polar plots and historical
comparison use the bundled CircStat toolbox. The tests use synthetic data.

From the PACT root in PowerShell, substitute the path to your EEGLAB installation.

```powershell
matlab -singleCompThread -batch "addpath('C:/path/to/eeglab2024.0'); eeglab nogui; addpath(pwd,'code/tests'); pac_setPath; test_pact_math; test_pact_reference; test_pact_integration; test_pact_menus; validate_pact_science; validate_pact_compatibility"
```

To reproduce the numerical comparison, extract the historical function into a
temporary folder under a different name. Then add that folder to the MATLAB path
and call `validate_pact_compatibility(@pact_legacy_fixture)`.

```powershell
git show 7d39d9c1fdb03b4c72635e2360f0c6dc7e834c8b:pac_man.m |
    ForEach-Object { $_ -replace '^function EEG = pac_man', 'function EEG = pact_legacy_fixture' } |
    Set-Content (Join-Path $env:TEMP 'pact_legacy_fixture.m')
```

## Local results

| Check | Result |
| --- | --- |
| `test_pact_math` | Pass. Checks parsing, p-value adjustment and reproducible surrogate draws. |
| `test_pact_reference` | Pass. Independent Holm and BH calculations include ties and extreme p-values. Both FFT FIR helpers preserve the impulse center and symmetry for low-pass, high-pass and band-pass modes. Tests check stopband and passband gains away from transition edges. |
| `test_pact_integration` | Pass. Checks raw samples and event metadata, original-coordinate HAS indices, deterministic output, legacy MI and circular summary definitions, boundary exclusion, mixed numeric and text markers, and an injected event-window positive control with finite significant p-values. The 100% HAS case retains every eligible sample. |
| `test_pact_menus` | Pass. Runs the 11 registered callbacks with scripted inputs and renders the MATLAB figures. Checks all seven cancel branches and cancellation of Compute PAC's follow-up statistics dialog with empty responses. |
| `validate_pact_compatibility` | Pass. Local one-hour fixture and historical comparison detailed below. |

A separate batch smoke run opened and dismissed the seven real EEGLAB
modal dialogs and verified that the datasets were preserved. The committed
menu test scripts the dialog responses for repeatability. The menu harness
covers Compute PAC, raw/HAS plotting, Copy event markers,
statistics setup, MI plotting, bar and polar angular histograms, phase-sorted
amplitude, event-related PAC and both scans. MATLAB draws the figures with
visibility disabled. Windows display appearance, interactive zoom and every
combination of scan options were outside the test. The HFO-frequency scan
completed with FIRfilt warnings about filter order in the smallest bands.
Its effective frequency responses require separate validation for a given grid.

## Scientific simulations

Each scenario uses 800 independently generated 90-second continuous
realizations sampled at 100 Hz. The 3 Hz oscillator has a stochastic phase
trajectory. A 28 Hz carrier has an independent, autocorrelated log-amplitude
envelope with AR(1) parameter 0.96. The modulation multiplier is
`1 + strength*cos(LFO phase)`. The nonstationary null applies a step gain
change halfway through and adds eight random broadband impulses. The
same realizations are reused across HAS rates and modulation strengths.
There are 199 surrogate draws per realization. The table reports uncorrected
rejections at alpha 0.05 over all 800 generated realizations. Missing p-values
are counted in the denominator.

| Synthetic scenario | 1% HAS | 2% HAS | 5% HAS | 20% HAS | 100% HAS |
| --- | ---: | ---: | ---: | ---: | ---: |
| Independent autocorrelated envelope | 41 | 39 | 45 | 34 | 9 |
| Injected coupling, strength 0.12 | 338 | 552 | 759 | 800 | 800 |
| Injected coupling, strength 0.50 | 800 | 800 | 800 | 800 | 800 |
| No injected modulation, gain step and impulses | 44 | 48 | 45 | 23 | 7 |

Three p-values were missing at 1% HAS. One came from the 0.12-strength
alternative and two came from the artifact null. The exact conditional-null
simulation fixed 20 unequal block magnitudes and drew independent uniform block
orientations. It rejected 52 of 1,000 draws. The continuous nulls show
conservative rejection at 100% HAS (1.13% and 0.88%), with the artifact null
also conservative at 20% HAS (2.88%). The observed 1–5% HAS null rates span
4.88–6.00%. Sampling uncertainty and dependence across the HAS-rate columns
limit comparisons across rates.

The revised partition distributes a trailing remainder across full blocks and
requires at least four occupied blocks. It returns missing inference when a
selected segment cannot hold a full block. Five cycles is a design heuristic.
The randomization assumes independence and uniform orientation of the block
resultants conditional on their magnitudes. Those assumptions can fail with
persistent rhythms, phase-locked events, artifacts or long dependence.
The test measures cross-block phase consistency of weighted resultants.
Evidence of consistency alone cannot establish amplitude modulation.
For example, the phase-locked marker control used a constant-amplitude
28 Hz carrier, a 4 Hz LFO, and ±40 ms windows centered at the same LFO phase.
It produced pooled MI 0.8044 and an uncorrected p-value of 0.005 at 100% HAS.
The event plot reports pooled MI without presenting a significance claim.

## Historical comparison and performance

A seeded two-channel, 12,000-sample continuous fixture without boundaries at
5% HAS matched the historical unweighted vector length, mean phase and
circular standard deviation and amplitude-weighted MI to within 1e-12.
The one-hour synthetic recording with 8 channels and 360,000 samples per
channel completed in 1.507 seconds with 199 surrogate draws and
`-singleCompThread`.
Timing is specific to this machine and configuration. The benchmark does not
establish maximum supported recording size.

The coordinate masks, per-segment Hilbert transform, guards that cover the
FIR half-order, all-sample selection at 100% HAS, event-window HAS selection,
phase-bin summaries and surrogate inference intentionally differ from
historical behavior. Legacy phase and amplitude distribution-test fields
remain present with NaN values. The Kaiser design uses its calculated
window and group delay. Testing across supported
MATLAB and EEGLAB versions remains outstanding.

## Output dimensions

Let C be the channel count, T the original sample count, B the number
of phase bins, W the accepted event count and L the odd event window length
in samples.

| Field | Size and meaning |
| --- | --- |
| `analyticEEG` | C x T. HFO envelope times the unit phasor of the LFO phase. Use the masks before using edge or boundary samples. Very short segments contain zeros. |
| `validSampleMask`, `selectedSampleMask` | 1 x T. The first excludes guards. The second gives eligible samples before HAS, intersected with the accepted event windows in event mode. |
| `hfoIndex` | C x 1 cell array of one-based original EEG sample indices. Below 100%, the amplitude must exceed the relevant percentile threshold. Ties and pooled thresholds can change the per-channel fraction. |
| `mi`, `vectLength`, `vectLengthVar`, `angleMean`, `angleStd` | C x 1. Absolute mean weighted vector, unweighted resultant length, 1 minus that length, circular mean phase in radians in [0, 2 pi), and CircStat's angular deviation. MI retains the envelope's amplitude units. |
| `moduInd95CI`, `moduInd99CI` | C x 1. Historical names for the 95th and 99th surrogate null percentiles. |
| `phaseBinValues` | C x (B+1). Quantiles of selected phase on [-pi, pi]. |
| `phaseProbability`, `phaseSortedAmp`, `phaseSortedAmpSe` | C x 1 cell arrays of equal-width phase-bin counts, mean amplitudes and descriptive standard errors ordered from 0 to 2 pi. Empty or unobserved bins lack amplitude estimates. The standard errors do not account for temporal dependence. |
| `histBinMax` | C x 1. Maximum equal-width bin count. |
| `uncorrected`, `Bonferroni`, `BonfHolm`, `FDR` | Structs of C x 1 p-value vectors. Finite values define each correction family. BH-FDR assumes independence or suitable positive dependence. No adjustment across scanned bands or HAS rates is applied. |
| `currentMIpval` and other `current*` fields | Active C x 1 p-value vectors plus display labels for the selected methods. |
| `winEdgesAfterRejection` | W x 2, inclusive first and last indices in the original recording. |
| `windowMeanAmp`, `windowStdAmp`, `windowMeanPhase`, `windowStdPhase` | C x L, computed over accepted event clips. Amplitude uses ordinary mean and standard deviation. Phase uses circular mean and circular standard deviation in radians. Overlapping clips can share samples. HAS uses their union. |

Keep this PR in draft pending scientific review of the block null for intended
datasets, especially event-locked and strongly nonstationary designs.
