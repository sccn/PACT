# PACT local validation

Branch `claude/pact-validation-2026-11`, stacked on PR #6
(`refactor/pact-scientific-safety-2026-11` at `6ec7ba0`). Synthetic data only;
no clinical data were used for any result in this file.

## Environment

Windows 11 Pro build 26200; Intel i9-10900K, 96 GB RAM. MATLAB R2024b Update 3
(24.2.0.2806996) with Signal Processing, Statistics and Machine Learning,
Optimization, Image Processing, Parallel Computing and Wavelet toolboxes.
EEGLAB 2024.0 with the bundled firfilt plugin (`pop_eegfiltnew`,
`findboundaries`, `pop_firws`). CircStat: bundled in `external/CircStat`
(upstream commit `302acc1`, see `external/README.md`).

## How to run

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File code\tests\run_pact_tests.ps1 `
    -Eeglab C:/path/to/eeglab2024.0 -IncludeProcessRegression -IncludeScience
```

Every MATLAB run goes through `code/tests/gui/run_matlab_with_timeout.ps1`,
which terminates only the MATLAB process tree it started and exits 124 on
timeout. Suites:

| Suite | Contents | Real dialogs |
| --- | --- | --- |
| mocked | `test_pact_math`, `test_pact_circstat`, `test_pact_reference`, `test_pact_integration`, `test_pact_dialog_guard`, `test_pact_menus` | No. `inputgui` is stubbed. |
| gui | `test_pact_dialog_timeout`, `smoke_pact_real_dialogs` (in `code/tests/gui`) | Yes, with no mocks. |
| process | `test_process_timeout.ps1` (~2 min) | Reproduces the original hang. |
| science | `validate_pact_science` (~6 min) | No. |

Historical comparison: extract master's `pac_man` as `pact_legacy_fixture`
(see the end of this file) and run `validate_pact_compatibility(@pact_legacy_fixture)`.

## Results reproduced on this branch (19 Nov 2026)

| Check | Result |
| --- | --- |
| `test_pact_math` | PASS. Holm/BH/Bonferroni, numeric parsing, and `pac_blockSurrogate` against an independent reference. Checks invariance to a global phase offset, exact zero dMI for constant amplitude at 100% HAS with a skewed phase distribution, short runs, and input errors. |
| `test_pact_circstat` | PASS. Bundled CircStat functions used by PACT match closed forms. In a one-off check, old 2012a and upstream outputs were bit-identical on 2000 random inputs. |
| `test_pact_reference` | PASS. Independent Holm/BH with ties; FIR impulse centering, symmetry and gains. |
| `test_pact_integration` | PASS. Data and events unchanged; original-coordinate HAS; boundaries; mixed numeric and text events; event windows; `dmi`/`numBlocks`; `windowStdPhase` equals `circ_std`. Two statistical regressions: a non-sinusoidal null (40 channels, ≤ 6 rejections allowed) and an event-locked constant-envelope control (p > 0.05 required). The PR #6 code at `6ec7ba0` fails both, with 40/40 rejections and p = 0.005. |
| `test_pact_dialog_guard` | PASS. Static and recursive. A fixture of the old `inputdlg` callback is flagged. |
| `test_pact_menus` (mocked) | PASS. All 11 callbacks and 7 cancel branches. |
| `test_pact_dialog_timeout` (real) | PASS. On a live `inputdlg`, the original exact-tag search finds no `'cancel'` button (`'Cancel'` is present). The revised harness cancels it in 0.5 s. A dialog with no Cancel button is removed and reported as a timeout after 2.4 s (deadline 2 s). |
| `smoke_pact_real_dialogs` (real) | PASS. All 7 real PACT dialogs canceled in about 1.3 s each; EEG and ALLEEG unchanged. |
| `test_process_timeout.ps1` | PASS. The original harness logic against a real `inputdlg` hung as in November. The wrapper terminated it after 91 s (limit 90 s) and left no process. |
| `validate_pact_compatibility` | PASS. On a boundary-free fixture, `mi`, `vectLength`, `angleMean` and `angleStd` match master `7d39d9c` within 9e-16. |
| Performance (2000 draws, 6 boundaries, 2% HAS) | 64 ch × 10 min @ 1 kHz: 9.4 s; 32 ch × 1 h @ 500 Hz: 25.2 s; 128 ch × 20 min @ 500 Hz: 16.1 s. MATLAB memory ≤ 3.6 GB. |

The block-rotation table in the previous version of this file (commit
`6ec7ba0`) was reproduced exactly. For example, the null rejections were
41/39/45/34/9 of 800, matching the earlier report. Its conclusion is
superseded by the findings below.

## Why the MI test changed

Further simulations showed that the PR's block-rotation null is invalid when
the LFO phase distribution is not uniform. With a non-sinusoidal 3 Hz LFO
(harmonics; 2–10 Hz phase band) and HFO amplitude independent of phase, it
rejected 20.4 / 90.6 / 100 / 100% of null channels at HAS 1/5/20/100%
(480 channels per cell; nominal 5%). It was also conservative at 100% HAS
(0.8–1.2%). Partial-cycle sums of the mean amplitude telescope across
adjacent blocks, which makes them negatively correlated.

The historical PACT null permutes single phase samples, which ignores
autocorrelation. On the same simulations, it rejected 4–56% of null channels
(table below).

**Revised method** (`pac_blockSurrogate`): block-debiased MI. In each block,
the block's mean candidate-sample phasor m_j is subtracted, so each HAS sample
contributes a_t·(exp(iφ_t) − m_j). The statistic is dMI = |Σ| / N_HAS. The null
applies random sign flips to the block contributions, which assumes
independent, centrally symmetric blocks but not isotropy. Before adoption,
four variants were compared on six scenarios. Local debiasing fixed both the
non-sinusoidal bias and the 100% conservativeness, and kept power: 39/96/100/100%
at k=0.12 versus 40/96/100/100% for rotation. The sign-flip null was chosen
over rotation because it needs the weaker assumption.

## Scientific simulation (`validate_pact_science`, 800 channels per cell)

Synthetic data: 90 s at 100 Hz. The 3 Hz LFO has a random-walk phase. The
28 Hz carrier has an autocorrelated log-normal envelope; coupling multiplies
the envelope by (1 + k·cos φ). There were 199 sign-flip draws per channel, and
the rejection rate uses uncorrected p ≤ 0.05 among channels with a finite p.
"Missing" counts channels with fewer than four occupied blocks, whose p-value
is NaN. Holm FWER is the share of 8-channel recordings with any Holm-adjusted
p ≤ 0.05. LEGACY is the historical permutation null on the same HAS samples.

| Scenario | HAS 1% | 2% | 5% | 20% | 100% | Legacy null (1–100%) |
| --- | --- | --- | --- | --- | --- | --- |
| N1 independent, AR envelope τ≈0.25 s | 5.1 | 5.6 | 5.9 | 4.4 | 4.4 | 56.4 / 53.9 / 50.1 / 39.4 / 30.4 |
| N2 long bursts τ≈1 s | 1.7 (214 missing) | 3.3 (23) | 4.8 | 4.6 | 4.4 | 39.0 / 30.4 / 22.2 / 11.6 / 4.2 |
| N3 non-sinusoidal LFO, 2–10 Hz | 4.8 | 5.9 | 5.0 | 4.6 | 5.2 | 24.0 / 20.1 / 17.6 / 15.4 / 12.0 |
| N4 pink 1/f noise | 5.1 | 4.6 | 5.4 | 5.1 | 4.5 | 36.1 / 37.4 / 37.1 / 34.8 / 28.5 |
| N5 EMG bursts, gain drift and step, 4 boundaries | 3.2 (40) | 4.2 | 4.6 | 5.4 | 4.1 | 48.4 / 46.0 / 40.4 / 29.9 / 16.5 |
| N6 event-locked mode, unrelated events | 1.9 (281) | 3.9 (27) | 4.9 | 6.0 | 5.4 | n/a |
| C1 8 broadband impulses per channel | 4.8 | 5.5 | 5.1 | 4.6 | 4.2 | n/a |

All 35 null cells have 95% Clopper–Pearson intervals that include 5% or lie
below it (largest: 5.9% [4.3, 7.7]). Holm FWER was 0–7% per 100 recordings.

Power (% of channels rejected):

| Coupling | HAS 1% | 2% | 5% | 20% | 100% |
| --- | --- | --- | --- | --- | --- |
| k = 0.05 | 10.4 | 15.9 | 28.2 | 68.4 | 95.1 |
| k = 0.10 | 29.6 | 54.2 | 86.9 | 100 | 100 |
| k = 0.20 | 81.9 | 98.5 | 100 | 100 | 100 |
| k = 0.30 | 98.1 | 100 | 100 | 100 | 100 |
| k = 0.50 | 99.6 | 100 | 100 | 100 | 100 |
| event-locked, k = 0.30 inside windows | 19.2 (179 missing) | 76.7 | 99.5 | 100 | 100 |
| event-locked, k = 0.50 inside windows | 31.6 (133) | 92.6 | 100 | 100 | 100 |

In this simulation, the synthetic coupling is strongest at 100% HAS. Low HAS
rates select few samples and lose power. This reflects the synthetic
generative model, not a recommendation of HAS rate.

## The 18 Nov 2026 menu-test hang

* Two Codex MATLAB runs (started 23:11 and 23:31) were idle, not computing:
  about 70 s of CPU in 30–50 min, 0 s in a 10-s sample, 1.4 GB, no windows.
  Their logs (`%TEMP%\pact-dialogs.log` and `pact-delivery-final.log`) were
  still held open.
* The 23:11 run's last line was `CANCEL PASS Set up statistics`. Its next
  callback, *Plot event-related PAC*, then used MATLAB's `inputdlg`. The timer
  searched for Tag `'cancel'`, but `inputdlg` uses `'Cancel'`, and `findall`
  is case-sensitive. So `uiwait` waited forever on an invisible figure. This
  was a test-harness defect, not a PACT computation.
* The 23:31 run stopped at the start of the same timer-driven real-dialog
  phase. A sibling run 30 s earlier crashed with an access violation
  (`iqm.dll`) in that phase. Its exact blocking point was not established.
* Both stale runs were terminated with authorization. Unrelated MATLAB
  sessions were not touched.
* Prevention, in three independent layers: the static recursive guard,
  per-dialog deadlines in the real-GUI harness, and the process-level wall
  clock. The regression tests above demonstrate each layer.

## Unresolved limitations

* The sign-flip test is exploratory. It assumes independent, centrally
  symmetric block contributions. Event-locked waveforms that produce both
  phase concentration and amplitude change, broadband transients filtered
  into both bands, and dependence longer than a block can violate this.
  None of these p-values is evidence of physiological PAC on its own.
* The earlier event-locked control (constant 28 Hz envelope, ±40 ms windows
  at the same 4 Hz phase, 100% HAS) gave p = 0.005 under the rotation null.
  Under the revised test it gave 0/20 rejections (min p = 0.185) while raw
  MI was 0.82, because block-local debiasing removes the shared window
  phase. Debiasing does not cover evoked responses with genuine amplitude
  change at a fixed phase; those still test positive and need trial-level
  designs.
* With few occupied blocks (low HAS, short recordings, long bursts, short
  event windows), p-values are NaN or conservative.
* Correction is across channels only, not across scanned frequencies, HAS
  rates, windows or conditions.
* `pac_scanHfoByPhase` reports "Filter order too low" for narrow amplitude
  bands. Its frequency responses need validation for the chosen grid. The
  legacy `pac_scanLfoPhaseFreq` refuses data with boundaries.
* `analyticEEG` contains zeros in segments shorter than the filter guard;
  use `validSampleMask`.
* Only MATLAB R2024b with EEGLAB 2024.0 was tested.

## Output fields

C = channels, T = samples, B = phase bins, W = accepted events, L = event
window length (odd).

| Field | Size and meaning |
| --- | --- |
| `analyticEEG` | C × T. HFO envelope × unit phasor of LFO phase, in original coordinates. Excluded samples remain; use the masks. |
| `validSampleMask`, `selectedSampleMask` | 1 × T. Valid excludes filter guards at edges and boundaries (≥ 1 s or half the FIR order). Selected is valid ∩ event windows in event mode. |
| `hfoIndex` | C × 1 cell of one-based original sample indices above the HAS threshold (all eligible samples at 100%). |
| `mi` | C × 1. Historical MI, descriptive. Same as master on boundary-free data. |
| `dmi` | C × 1. Block-debiased MI, the statistic that the p-values test. |
| `numBlocks`, `blockLength` | C × 1 occupied blocks; scalar block length in samples. |
| `vectLength`, `vectLengthVar`, `angleMean`, `angleStd` | C × 1. Resultant length, 1 − length, mean phase in [0, 2π), angular deviation √(2(1−R)) (CircStat `circ_std`). |
| `moduInd95CI`, `moduInd99CI` | C × 1. 95th/99th percentiles of the sign-flip null of `dmi`. These are historical names, not confidence intervals. |
| `phaseBinValues`, `phaseProbability`, `phaseSortedAmp`, `phaseSortedAmpSe`, `histBinMax` | Descriptive phase-bin summaries ordered 0 to 2π. The SEs ignore autocorrelation. |
| `uncorrected`, `Bonferroni`, `BonfHolm`, `FDR` | Structs of C × 1 p-values. `MIpval` refers to `dmi`. Legacy distribution-test fields are NaN. |
| `winEdgesAfterRejection` | W × 2 inclusive window edges, centered on events. |
| `windowMeanAmp`, `windowStdAmp`, `windowMeanPhase`, `windowStdPhase` | C × L. Event-window summaries. Phase uses the circular mean and angular deviation, as in master. |
| `surrogateMethod`, `inferenceNotice`, `rngSeed` | Method description, caveats, and the fixed seed (0). The caller's RNG is restored. |

### Historical fixture

```powershell
git show 7d39d9c1fdb03b4c72635e2360f0c6dc7e834c8b:pac_man.m |
    ForEach-Object { $_ -replace '^function EEG = pac_man', 'function EEG = pact_legacy_fixture' } |
    Set-Content (Join-Path $env:TEMP 'pact_legacy_fixture.m')
```
