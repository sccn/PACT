# PACT validation (PR #6)

Round 2 of local validation, 19 November 2026, on branch
`claude/pact-validation-round2`, stacked on `07931e4` (round 1). Synthetic
data only, except a read-only technical check on local PREVeNT pilot files.
No result in this file comes from an external service or from data leaving
the machine. **Tested commit: see the PR comment for the exact SHA.**

Every result below was **executed locally** unless it is marked
*documented, not re-run*.

## 1. Environment (exact versions)

| Component | Version |
| --- | --- |
| OS / CPU | Windows 11 Pro 10.0.26200; Intel i9-10900K; 96 GB RAM |
| MATLAB | R2024b Update 3 (24.2.0.2806996) |
| Toolboxes | Signal Processing 24.2, Statistics and Machine Learning 24.2 (also installed: Optimization, Image Processing, Parallel Computing, Wavelet) |
| EEGLAB | 2024.0 |
| FIRfilt | 2.8 (bundled with EEGLAB 2024.0: `pop_eegfiltnew`, `pop_firws`, `firfilt`, `findboundaries`) |
| CircStat | `circstat/circstat-matlab` @ `302acc161d33bdb5c432243b946b178d88c361b8` (bundled in `external/CircStat`) |

## 2. How to run

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File code\tests\run_pact_tests.ps1 `
    -Eeglab C:/path/to/eeglab2024.0 -IncludeProcessRegression -IncludeScience
```

Each suite runs in its own MATLAB process under
`code/tests/gui/run_matlab_with_timeout.ps1`. That wrapper kills only the
process tree it started (`taskkill /T`) and exits 124 on timeout, so no
orphan MATLAB processes remain. The long grids run separately:

```matlab
validate_pact_null('csv','dev.csv')                % development grid, about 45 min
validate_pact_null('holdout','csv','holdout.csv')  % holdout grid, about 10 min
validate_pact_null('centering')                    % centering variants, about 1.5 min
validate_pact_events                               % event mode, about 10 min
validate_pact_filters('plots','images/validation_2026-11')
```

All scripts use fixed seeds (`rng(seed,'twister')`) and print the parameter
values they use.

**Timing note.** firfilt calls a bare `tic` internally, which resets
MATLAB's global timer. Any `tic ... toc` around PACT therefore under-reports
elapsed time. Round 1's `pac_man` timings (9.4 / 25.2 / 16.1 s) were about
20% too low. All timings below use timer handles (`t = tic; ... toc(t)`).

## 3. Pass/fail matrix (round 2, executed)

| Suite / script | Result |
| --- | --- |
| `test_pact_math` (exact sign-flip enumeration, p floor 2^(1−K), degenerate zero, blockCycles, Monte Carlo path, event shuffle vs enumerated reference) | PASS |
| `test_pact_circstat` | PASS |
| `test_pact_reference` (Holm/BH with ties; FIR impulse, symmetry, gains) | PASS |
| `test_pact_integration` (immutability, coordinates, boundaries, events, 3 statistical regressions, under-powered warning) | PASS |
| `test_pact_filters` (filt_fir1fft spec, scan transition widths, pop_eegfiltnew zero phase) | PASS |
| `test_pact_scans` (independent recomputation of both scans, immutability, dimensions, 100% HAS, mixed boundaries, unsupported configurations) | PASS |
| `test_pact_dialog_guard` (static, recursive) | PASS |
| `test_pact_menus` (mocked callbacks) | PASS |
| `test_pact_dialog_timeout` (real `inputdlg`) | PASS |
| `smoke_pact_real_dialogs` (7 real PACT dialogs) | PASS |
| `test_process_timeout.ps1` (original hang reproduced, killed at the limit) | PASS (re-run this round: killed at 91 s with a 90 s limit, no processes left) |
| `validate_pact_compatibility` vs master `7d39d9c` | PASS (mi, vectLength, angleMean, angleStd ≤ 9e-16) |
| `validate_pact_null` development grid | Calibrated within sampling error. Pooled 4.3% (191k tests); flagged excesses not reproduced with fresh seeds (§4.4) |
| `validate_pact_null` holdout grid (+ N02 replication) | Pooled 5.0%; the 1600-channel replication shows a slight liberal excess of about +0.5 pp (5.4–5.7%) with no identified mechanism (§4.4) |
| `validate_pact_null('centering')` | Candidate centering is the only calibrated *and* powerful variant (§4.3) |
| `validate_pact_events` | Calibrated; evoked-only responses not detected (§4.6) |
| `validate_pact_filters` | Measured; 3 scan filter defects found and fixed (§5) |

## 4. Statistical validation of the MI tests

### 4.1 Estimands and nulls

**`mi` (descriptive, historical).** For channel c with HAS sample set H,
`mi` = |N⁻¹ Σ_{t∈H} a_t e^{iφ_t}|, where a is the HFO envelope, φ the LFO
phase and N = |H|. It is unchanged from master (§6). No p-value refers to
it, except in event mode.

**`dmi` (tested in continuous modes).** Valid samples are split into blocks
B_1..B_J within boundary segments. Each block spans at least
max(0.5 s, `blockCycles`/f_low) (default 5 cycles of the lowest LFO band
edge); a segment shorter than that forms one block. Let C_j be the block's
candidate samples (valid, and inside event windows in event mode), and
m_j = |C_j|⁻¹ Σ_{t∈C_j} e^{iφ_t}. Then

  Z_j = Σ_{t∈H∩B_j} a_t (e^{iφ_t} − m_j),  dmi = |Σ_j Z_j| / N.

The **null** is that, over the K occupied blocks, the Z_j are independent
and centrally symmetric (Z_j and −Z_j are equally likely). Under it,
|Σ_j s_j Z_j| has the same distribution for every sign vector s. The
p-value is the fraction of the 2^(K−1) distinct sign patterns (the first
sign fixed) with |Σ s_j Z_j| ≥ |Σ Z_j|. It is exact for K ≤ 16, and Monte
Carlo, (1 + b)/(nSurro + 1), otherwise. The smallest attainable p is
2^(1−K) (`minAttainableP`), so K = 4 or 5 can never reach 0.05;
`pac_man` warns (`PACT:Underpowered`).

The null concerns **statistical dependence** between HFO amplitude and LFO
phase. It does not distinguish physiological PAC from waveform-shape or
common-input dependence (§4.5).

**Why centering by candidate samples.** If amplitude is independent of
phase and the HAS propensity is stationary within a block,
E[Z_j | φ] = c_j Σ_{t∈C_j}(e^{iφ_t} − m_j) = 0 exactly, *whatever the phase
distribution*. That removes the phase-clustering bias of non-sinusoidal
LFOs. It also removes the partial-cycle remainder of the mean amplitude,
which made the PR's block-rotation null conservative at 100% HAS.

**hfoTopRatio = 100 with constant amplitude.** Then H = C, so
Z_j = a Σ_{C_j}(e^{iφ} − m_j) = 0 identically for any phase distribution,
and dmi = 0. Before this round, rounding residue (~1e-16) was compared with
an equally tiny null, so p was arbitrary and could fall below 0.05. It is
now reported as dmi = 0, p = 1 (bug fix; regression in `test_pact_math`).

**Event mode (`hfoPool` 4) tests `mi`** against a trial-shuffle null
(`pac_eventShuffle`). For non-overlapping, equal-length windows w = 1..W,
the null pairs window w's amplitudes and HAS positions with window π(w)'s
phases at the same latency, over random permutations π (exact for W ≤ 8,
otherwise Monte Carlo). It tests event-to-event covariation of amplitude
with phase. Anything identical across events survives the shuffle and is
not reported as coupling.

### 4.2 Independent implementation and agreement

`validate_pact_null` contains its own generator and a reference
implementation of dmi and the sign-flip test (`refTest`). The reference uses
explicit loops and its own `RandStream`, and calls no PACT function. Over
the development grid (batch 1, all cells), production vs reference:
statistic max |Δ| = 7.1e-15; exact p mismatches 0 of 23,274; Monte Carlo |Δp| median 0.020 (n = 48,433). The Monte Carlo differences match independent random draws
with 199 permutations (SE ≈ 0.035 near p = 0.5). The grid path (pac_man at
100% HAS, then pac_blockSurrogate per HAS × block length) was checked
against end-to-end `pac_man` (identical HAS sets and dmi ≤ 1e-10).

### 4.3 Candidate vs HAS-only phase means (`validate_pact_null('centering')`)

Reference implementation; 200 channels per row; 60 s; 5-cycle blocks; % of
channels rejected at α = 0.05.

| Scenario | HAS | Candidate mean (production) | HAS-only mean | No centering |
| --- | --- | --- | --- | --- |
| N01 null | 1 / 5 / 20 / 100% | 5.5 / 4.0 / 3.5 / 2.5 | 5.0 / 5.5 / 5.0 / 2.5 | 6.0 / 4.0 / 2.5 / 1.0 |
| N02 null, non-sinusoidal | 1 / 5 / 20 / 100% | 5.5 / 5.0 / 6.0 / 5.5 | 7.0 / 4.5 / 4.0 / 5.5 | **23 / 77 / 100 / 100** |
| P01 k = 0.10 | 1 / 5 / 20 / 100% | 20 / 40.5 / 75 / 90 | **5.0 / 8.5 / 16** / 90 | 19.5 / 39 / 74 / 86.5 |
| P02 k = 0.30 | 1 / 5 / 20 / 100% | 85.5 / 100 / 100 / 100 | **6.5 / 28 / 81.5** / 100 | 85.5 / 100 / 100 / 100 |

Without centering, the test is invalid for non-sinusoidal LFOs. HAS-only
centering subtracts the HAS phase clustering, which is the signal, so power
collapses at low HAS. At 100% HAS the two centering variants coincide.

### 4.4 Null calibration grid (`validate_pact_null`)

**Simulator** (fixed seeds, base 1000; parameters printed by the script):
- fs 200 Hz; LFO 3 Hz with a random-walk phase (diffusion 0.25 rad²/s);
  phase band 2–4 Hz (2–10 Hz for N02).
- HFO carrier: 40–60 Hz band-limited noise, unit variance. Envelope:
  exp(σ·AR(1)) with τ 0.25 s and σ 0.5 (N07: σ 1; N05: 1/f log-envelope).
- Background: 0.3 × pink noise. Amplitude band 35–70 Hz.
- 100 channels × 4 batches = **400 independent channels per cell**.
- Grid: HAS 1/2/5/20/100% × blockCycles 3/5/10/20 × 60 s and 300 s;
  199 permutations. Holm FWER is computed over families of 10 channels.

Cells of one scenario × length share the same channels, so they are
correlated, and pooled rates are not independent replicates.

Development grid: 480 null cells; pooled 4.30% (8212/190879 tests, cells correlated); per-cell range 0.0–8.5% (median 4.8%); 7 cells with CI entirely above 5%, 76 entirely below; Holm FWER per cell 0.0–22.5% (median 2.5%); 1121 missing (0.6%). By block length: 3 cyc 4.81%, 5 cyc 4.84%, 10 cyc 4.67%, 20 cyc 2.87%. Few-block regime (K ≤ 8): 31866 channel-tests, 1.0% rejected, 23236 with minimum attainable p > 0.05.

Pooled by scenario × length (20 cells share 400 channels, so the effective n is about 400):

| Scenario / length | Type I (pooled) | Tests | Missing (excluded) | Holm FWER |
| --- | --- | --- | --- | --- |
| N01 / 300 s | 4.94% | 8000 | 0 (0.0%) | 6.4% |
| N01 / 60 s | 4.29% | 7995 | 5 (0.1%) | 3.5% |
| N02 / 300 s | 4.16% | 8000 | 0 (0.0%) | 3.4% |
| N02 / 60 s | 3.80% | 7993 | 7 (0.1%) | 3.0% |
| N03 / 300 s | 5.21% | 8000 | 0 (0.0%) | 4.5% |
| N03 / 60 s | 3.94% | 7990 | 10 (0.1%) | 3.8% |
| N04 / 300 s | 5.44% | 8000 | 0 (0.0%) | 6.1% |
| N04 / 60 s | 3.01% | 8000 | 0 (0.0%) | 2.6% |
| N05 / 300 s | 5.92% | 8000 | 0 (0.0%) | 3.1% |
| N05 / 60 s | 3.69% | 7961 | 39 (0.5%) | 3.2% |
| N06 / 300 s | 5.69% | 8000 | 0 (0.0%) | 3.9% |
| N06 / 60 s | 2.74% | 7997 | 3 (0.0%) | 1.5% |
| N07 / 300 s | 6.74% | 8000 | 0 (0.0%) | 8.5% |
| N07 / 60 s | 2.35% | 7968 | 32 (0.4%) | 2.5% |
| N08 / 300 s | 5.65% | 8000 | 0 (0.0%) | 7.2% |
| N08 / 60 s | 3.80% | 7533 | 467 (5.8%) | 2.9% |
| N09 / 300 s | 4.61% | 8000 | 0 (0.0%) | 3.6% |
| N09 / 60 s | 3.74% | 7988 | 12 (0.1%) | 2.1% |
| N10 / 300 s | 5.09% | 8000 | 0 (0.0%) | 4.0% |
| N10 / 60 s | 3.15% | 7990 | 10 (0.1%) | 2.4% |
| N11 / 300 s | 4.28% | 8000 | 0 (0.0%) | 4.6% |
| N11 / 60 s | 3.26% | 7998 | 2 (0.0%) | 3.2% |
| N12 / 300 s | 4.39% | 8000 | 0 (0.0%) | 2.6% |
| N12 / 60 s | 3.24% | 7466 | 534 (6.7%) | 3.3% |

Pooled by block length:

| Block length | Type I (pooled) | Tests | Missing (excluded) | Holm FWER |
| --- | --- | --- | --- | --- |
|  3 cycles | 4.81% | 47962 | 38 (0.1%) | 4.7% |
|  5 cycles | 4.84% | 47929 | 71 (0.1%) | 3.9% |
| 10 cycles | 4.67% | 47800 | 200 (0.4%) | 4.2% |
| 20 cycles | 2.87% | 47188 | 812 (1.7%) | 2.6% |

Pooled by HAS rate:

| HAS | Type I (pooled) | Tests | Missing (excluded) | Holm FWER |
| --- | --- | --- | --- | --- |
|   1% | 4.06% | 37459 | 941 (2.5%) | 3.2% |
|   2% | 4.26% | 38226 | 174 (0.5%) | 3.9% |
|   5% | 4.40% | 38394 | 6 (0.0%) | 4.5% |
|  20% | 4.39% | 38400 | 0 (0.0%) | 4.0% |
| 100% | 4.40% | 38400 | 0 (0.0%) | 3.5% |

**N01, all cells** — % rejected (α = 0.05), 400 channels per cell; bold = 95% CI entirely above 5%; missing = channels with < 4 occupied blocks or < 36 HAS samples (excluded from the denominator).

| Scenario | Length | HAS | 3 cyc | 5 cyc | 10 cyc | 20 cyc | missing | Holm FWER (5 cyc) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| N01 | 60 s | 1% | 5.2 | 5.0 | 7.0 | 0.0 | 5 | 7.5 |
| N01 | 60 s | 2% | 4.8 | 4.5 | 6.0 | 0.0 | 0 | 7.5 |
| N01 | 60 s | 5% | 4.8 | 6.0 | 6.0 | 0.0 | 0 | 7.5 |
| N01 | 60 s | 20% | 5.8 | 6.2 | 7.0 | 0.0 | 0 | 5.0 |
| N01 | 60 s | 100% | 5.5 | 6.8 | 5.2 | 0.0 | 0 | 0.0 |
| N01 | 300 s | 1% | 5.2 | 5.5 | 5.0 | 5.8 | 0 | 7.5 |
| N01 | 300 s | 2% | 5.2 | 5.2 | 4.8 | 4.8 | 0 | 7.5 |
| N01 | 300 s | 5% | 4.5 | 3.8 | 4.8 | 4.8 | 0 | 2.5 |
| N01 | 300 s | 20% | 4.8 | 4.5 | 4.8 | 6.0 | 0 | 5.0 |
| N01 | 300 s | 100% | 4.8 | 4.2 | 4.8 | 5.8 | 0 | 2.5 |

Holm FWER is estimated from 40 families per cell (95% CI roughly ±7 points at 5%). The maximum, 22.5% (9/40), comes from N07 300 s, the dataset whose excess did not replicate; the replication's per-cell FWER was 0–8.8% (median 3.8%).

Cells whose 95% CI lies entirely above 5%: N04 300 s HAS 5% 3 cyc (30/400), N05 300 s HAS 20% 3 cyc (30/400), N07 300 s HAS 1% 20 cyc (34/400), N07 300 s HAS 2% 20 cyc (30/400), N07 300 s HAS 20% 20 cyc (31/400), N07 300 s HAS 100% 5 cyc (30/400), N07 300 s HAS 100% 20 cyc (33/400). Under exact calibration, about 2.5% of 480 cells (12) would do so by chance. **Follow-up with fresh seeds:** N07 300 s re-run end-to-end in production (800 channels per cell): 40 null cells; pooled 4.13% (1319/31945 tests, cells correlated); per-cell range 0.0–6.2% (median 4.6%); 0 cells with CI entirely above 5%, 6 entirely below; Holm FWER per cell 0.0–8.8% (median 3.8%); 55 missing (0.2%). By block length: 3 cyc 4.76%, 5 cyc 4.85%, 10 cyc 4.68%, 20 cyc 2.21%. Few-block regime (K ≤ 8): 5391 channel-tests, 0.9% rejected, 4103 with minimum attainable p > 0.05. The N05 300 s excess (5.9%) was not reproduced in an independent probe (3.5–5.0%, 400 channels). Burst-straddling of block edges was ruled out: HAS weight within 50 ms of an edge was 4.3% vs 4.2% expected, and edges at amplitude minima changed nothing.

The 60 s recordings are conservative (pooled 3.4%) because short recordings and long blocks leave few blocks. With 20 cycles, 60 s gives about 5 blocks, so the minimum p (≥ 0.0625) is above α; `pac_man` now warns (`PACT:Underpowered`).

**Adjacent-block dependence.** The mean cos(angle(Z_j), angle(Z_{j+1})) at
100% HAS and 5 cycles was between −0.04 and +0.02 for every null and
control scenario (quick run, 100 channels each). There is no evidence that
consecutive block contributions are dependent at the default block length.

**Holdout grid** (never used during development). fs 500 Hz, LFO 6 Hz,
band 4–8 Hz (N02: 4–16 Hz), carrier 90–130 Hz, band 80–150 Hz, 120 s,
blockCycles 5 and 10, seed base 9000, 400 channels per cell:

Holdout grid: 60 null cells; pooled 5.03% (1206/24000 tests, cells correlated); per-cell range 3.2–8.0% (median 4.8%); 2 cells with CI entirely above 5%, 0 entirely below; Holm FWER per cell 0.0–15.0% (median 5.0%); 0 missing (0.0%). By block length: 5 cyc 5.11%, 10 cyc 4.94%. Few-block regime (K ≤ 8): 1 channel-tests, 0.0% rejected, 0 with minimum attainable p > 0.05.

Pooled by scenario:

| Scenario / length | Type I (pooled) | Tests | Missing (excluded) | Holm FWER |
| --- | --- | --- | --- | --- |
| N01 / 120 s | 4.03% | 4000 | 0 (0.0%) | 5.2% |
| N02 / 120 s | 6.45% | 4000 | 0 (0.0%) | 7.5% |
| N05 / 120 s | 4.83% | 4000 | 0 (0.0%) | 5.0% |
| N06 / 120 s | 5.50% | 4000 | 0 (0.0%) | 6.0% |
| N09 / 120 s | 4.55% | 4000 | 0 (0.0%) | 3.2% |
| N12 / 120 s | 4.80% | 4000 | 0 (0.0%) | 4.2% |

Holdout replication (N01, N02; 1600 channels per cell, fresh seeds): 20 null cells; pooled 5.55% (1776/32000 tests, cells correlated); per-cell range 5.1–6.2% (median 5.5%); 3 cells with CI entirely above 5%, 0 entirely below; Holm FWER per cell 3.8–10.6% (median 5.6%); 0 missing (0.0%). By block length: 5 cyc 5.53%, 10 cyc 5.57%. Few-block regime (K ≤ 8): 0 channel-tests, 0.0% rejected, 0 with minimum attainable p > 0.05.

Replication pooled by scenario × block length:

| Scenario / blocks | Type I (pooled) | Tests | Missing (excluded) | Holm FWER |
| --- | --- | --- | --- | --- |
| N01 / 10 cyc | 5.54% | 8000 | 0 (0.0%) | 5.4% |
| N01 / 5 cyc | 5.35% | 8000 | 0 (0.0%) | 5.4% |
| N02 / 10 cyc | 5.60% | 8000 | 0 (0.0%) | 6.5% |
| N02 / 5 cyc | 5.71% | 8000 | 0 (0.0%) | 7.0% |

In this configuration the test is slightly liberal, by about +0.4–0.7 percentage points (cells 5.1–6.2%; 3 of 20 cells have a CI just above 5%). A moment diagnostic found no mechanism. The ratio E|ΣZ_j|² / E Σ|Z_j|², which equals 1 for independent, mean-zero blocks, was 0.96–1.01 (± 0.035). Mean block skewness was 0.00 ± 0.02. Removing envelope leakage or background noise did not change the rate (all 4.5–5.6% at 800 channels).

### 4.5 Negative controls: dependence that is not physiological PAC

Here the HFO amplitude genuinely depends on LFO phase, through waveform shape
or shared input. Rejections show what the test cannot exclude. They are not
type-I errors.

**Dependence controls (development grid)** — % rejected, 400 channels per cell.

| Scenario | Length | HAS | 3 cyc | 5 cyc | 10 cyc | 20 cyc |
| --- | --- | --- | --- | --- | --- | --- |
| C01 | 60 s | 1% | 6.5 | 5.8 | 6.0 | 0.0 |
| C01 | 60 s | 2% | 5.5 | 5.8 | 6.0 | 0.0 |
| C01 | 60 s | 5% | 7.2 | 5.5 | 7.5 | 0.0 |
| C01 | 60 s | 20% | 8.2 | 7.5 | 7.2 | 0.0 |
| C01 | 60 s | 100% | 12.8 | 11.8 | 11.8 | 0.0 |
| C01 | 300 s | 1% | 5.0 | 6.2 | 4.8 | 5.2 |
| C01 | 300 s | 2% | 5.2 | 5.2 | 4.2 | 5.0 |
| C01 | 300 s | 5% | 5.5 | 5.2 | 5.2 | 5.8 |
| C01 | 300 s | 20% | 13.8 | 14.0 | 14.8 | 13.5 |
| C01 | 300 s | 100% | 45.0 | 47.0 | 47.5 | 47.0 |
| C02 | 60 s | 1% | 4.5 | 3.8 | 4.5 | 0.0 |
| C02 | 60 s | 2% | 5.2 | 4.8 | 5.2 | 0.0 |
| C02 | 60 s | 5% | 5.8 | 6.8 | 5.2 | 0.0 |
| C02 | 60 s | 20% | 3.8 | 5.0 | 4.8 | 0.0 |
| C02 | 60 s | 100% | 5.0 | 5.0 | 5.5 | 0.0 |
| C02 | 300 s | 1% | 7.0 | 7.5 | 7.8 | 7.2 |
| C02 | 300 s | 2% | 5.5 | 6.5 | 5.5 | 7.0 |
| C02 | 300 s | 5% | 7.0 | 7.2 | 7.0 | 6.0 |
| C02 | 300 s | 20% | 6.8 | 5.5 | 6.5 | 5.8 |
| C02 | 300 s | 100% | 5.0 | 7.0 | 6.2 | 6.2 |
| C03 | 60 s | 1% | 7.0 | 5.5 | 5.0 | 0.0 |
| C03 | 60 s | 2% | 7.0 | 4.8 | 6.2 | 0.0 |
| C03 | 60 s | 5% | 6.2 | 5.2 | 6.8 | 0.0 |
| C03 | 60 s | 20% | 8.0 | 9.8 | 8.8 | 0.0 |
| C03 | 60 s | 100% | 9.5 | 10.8 | 11.8 | 0.0 |
| C03 | 300 s | 1% | 4.8 | 5.2 | 5.5 | 6.2 |
| C03 | 300 s | 2% | 6.5 | 6.5 | 6.5 | 6.2 |
| C03 | 300 s | 5% | 8.2 | 7.0 | 7.5 | 6.5 |
| C03 | 300 s | 20% | 23.8 | 22.0 | 22.5 | 23.5 |
| C03 | 300 s | 100% | 29.8 | 29.0 | 30.5 | 30.5 |

C01 (sawtooth) is increasingly detected with HAS and length (up to 47% at 100% HAS, 300 s). C02 (spike-and-wave at random times) is barely detected. In continuous mode, C03 (an identical evoked ERP plus HFO burst every 2.5 s) is flagged in up to 30% of channels at 20–100% HAS. Use event mode for event-locked designs (§4.6).

**Power (development grid)** — % rejected, 400 channels per cell.

| Scenario | Length | HAS | 3 cyc | 5 cyc | 10 cyc | 20 cyc |
| --- | --- | --- | --- | --- | --- | --- |
| P01 | 60 s | 1% | 17.0 | 17.8 | 13.8 | 0.0 |
| P01 | 60 s | 2% | 25.0 | 24.2 | 20.8 | 0.0 |
| P01 | 60 s | 5% | 44.0 | 41.5 | 38.8 | 0.0 |
| P01 | 60 s | 20% | 77.2 | 75.5 | 73.2 | 0.0 |
| P01 | 60 s | 100% | 90.0 | 89.0 | 87.2 | 0.0 |
| P01 | 300 s | 1% | 68.5 | 67.5 | 69.0 | 67.5 |
| P01 | 300 s | 2% | 86.5 | 87.5 | 87.5 | 86.8 |
| P01 | 300 s | 5% | 97.5 | 98.2 | 99.0 | 98.2 |
| P01 | 300 s | 20% | 100.0 | 100.0 | 100.0 | 100.0 |
| P01 | 300 s | 100% | 100.0 | 100.0 | 100.0 | 100.0 |
| P02 | 60 s | 1% | 87.5 | 86.2 | 78.5 | 0.0 |
| P02 | 60 s | 2% | 98.2 | 97.8 | 96.2 | 0.0 |
| P02 | 60 s | 5% | 100.0 | 100.0 | 100.0 | 0.0 |
| P02 | 60 s | 20% | 100.0 | 100.0 | 100.0 | 0.0 |
| P02 | 60 s | 100% | 100.0 | 100.0 | 100.0 | 0.0 |
| P02 | 300 s | 1% | 100.0 | 100.0 | 100.0 | 100.0 |
| P02 | 300 s | 2% | 100.0 | 100.0 | 100.0 | 100.0 |
| P02 | 300 s | 5% | 100.0 | 100.0 | 100.0 | 100.0 |
| P02 | 300 s | 20% | 100.0 | 100.0 | 100.0 | 100.0 |
| P02 | 300 s | 100% | 100.0 | 100.0 | 100.0 | 100.0 |
| P03 | 60 s | 1% | 3.2 | 2.8 | 1.2 | 0.0 |
| P03 | 60 s | 2% | 2.8 | 1.8 | 0.2 | 0.0 |
| P03 | 60 s | 5% | 1.5 | 0.5 | 0.0 | 0.0 |
| P03 | 60 s | 20% | 0.2 | 0.0 | 0.0 | 0.0 |
| P03 | 60 s | 100% | 0.0 | 0.0 | 0.0 | 0.0 |
| P03 | 300 s | 1% | 6.0 | 5.8 | 3.0 | 1.5 |
| P03 | 300 s | 2% | 4.0 | 3.8 | 2.0 | 0.5 |
| P03 | 300 s | 5% | 3.8 | 1.5 | 0.5 | 0.0 |
| P03 | 300 s | 20% | 1.0 | 0.0 | 0.0 | 0.0 |
| P03 | 300 s | 100% | 0.0 | 0.0 | 0.0 | 0.0 |

P03 (k = 0.3 with the preferred phase rotating once per minute) is detected in 0–6% of channels, at or below the null rate, because block contributions cancel. The test requires a consistent preferred phase. With 20-cycle blocks on 60 s data, no channel can reach p ≤ 0.05 at all (0.0 columns; see `minAttainableP`).

### 4.6 Event mode (`validate_pact_events`, 400 channels per cell)

fs 200 Hz; events every 2.5 s ± 0.3 s for 120 s (about 47 events); ±500 ms
windows; HAS 5/20/100%; 199 permutations. E07 uses 6 events (exact,
720 permutations).

| Scenario | Expected | HAS 5% | HAS 20% | HAS 100% |
| --- | --- | --- | --- | --- |
| E01 events unrelated to LFO | ≈ 5% | 3.8 [2.1, 6.1] | 3.2 [1.7, 5.5] | 4.0 [2.3, 6.4] |
| E02 events at LFO peaks, flat envelope | ≈ 5% | 5.0 [3.1, 7.6] | 4.8 [2.9, 7.3] | 4.5 [2.7, 7.0] |
| E03 identical evoked ERP + 50 Hz burst | ≈ 5% | 4.5 [2.7, 7.0] | 5.2 [3.3, 7.9] | 5.5 [3.5, 8.2] |
| E07 E01 with 6 events (exact) | ≈ 5% | 4.8 [2.9, 7.3] | 5.2 [3.3, 7.9] | 4.0 [2.3, 6.4] |
| E04 PAC inside windows, k = 0.3 | power | 100 | 100 | 100 |
| E05 PAC inside windows, k = 0.6 | power | 100 | 100 | 100 |
| E06 PAC with phase reset at events, k = 0.6 | power | 99.0 | 100 | 100 |
| **E03c** E03 data, *continuous* block test | contrast | 5.8 | **12.0** | **13.8** |

The block test (continuous mode) flags identical evoked responses as
dependence. The trial shuffle does not, which is why event mode needs its
own null. The earlier event-locked control (constant envelope, ±40 ms
windows at the same phase), which gave p = 0.005 under the round-0
rotation null, now gives p > 0.05 (regression in `test_pact_integration`).
E06 shows that phase reset alone does not hide coupling. Phase still varies
across events within the window, and the shuffle uses that variability.
Only coupling whose phase is identical across events at every latency is
undetectable; it cannot be distinguished from an evoked response.

## 5. Filters and frequency scans (`validate_pact_filters`, `test_pact_scans`)

Impulse responses were measured through the production code, analysed with
`freqz` and checked with centre-frequency sinusoids. Plots:
`images/validation_2026-11/filters_*.png` (dotted lines show legacy orders).

| Path | Design | Findings |
| --- | --- | --- |
| `pac_man` (`pop_eegfiltnew`) | Hamming windowed sinc, zero-phase | Impulse response equals the design exactly (0 error) and is symmetric. −6 dB edges equal the reported cutoffs (e.g. PREVeNT 3–4 Hz → 2.00–5.00 Hz; 35–70 Hz → 30.66–74.35 Hz at 200 Hz). Passband ripple ≤ 0.04 dB; stopband −47 to −52 dB; sinusoid gain 0.999–1.003 with 0.00° lag. The guard max(1 s, order/2) covers the transient (up to 3.3 s for 0.5–1 Hz at 1 kHz). |
| `pac_scanHfoByPhase` (`pop_firws`), legacy | One transition width per family | **Defect.** On the GUI default grid (0.5–8 Hz, 10 log bins) the transition widths exceeded the low bins' widths. The passbands collapsed: centre gain −3.43 / −1.64 / −0.44 dB (sinusoid gain 0.67 / 0.83 / 0.95) for the three lowest phase bins and down to −0.73 dB for amplitude bins, with "Filter order too low" warnings. The guard was a fixed 1 s against up to 4 s of transient. |
| `pac_scanHfoByPhase`, revised | Each band's transition capped at the band width (`pac_firwsOrder`) | Centre gain −0.06 to +0.04 dB; −6 dB edges on the requested bin edges (±0.005 Hz); no warnings. Guard = half the longest filter (up to 9.1 s on the default grid at 500 Hz). Hilbert transforms run per boundary segment. |
| `pac_scanLfoPhaseFreq` (`filt_fir1fft`), legacy | Kaiser, deviation 0.1 | **Defect.** `kaiserord` chose β = 0 (a rectangular window). Measured passband-edge loss was 2.1 dB against a stated ±0.9 dB, with a 23 dB stopband. **No edge guard.** Numeric −99 boundaries were not refused. |
| `pac_scanLfoPhaseFreq`, revised | Kaiser, deviation 0.01, verified with `freqz` and lengthened until met | Ripple ≤ 0.09 dB; stopband ≤ −40.8 dB; zero phase (alignment ≤ 6e-17). Half the longest filter is excluded at each edge. Filters are about 3× longer for the lowest bands (12.6 s of support for 0.5–0.61 Hz at 500 Hz). −6 dB bands extend about 10% beyond each requested edge (20% transitions), so adjacent rows overlap. |

**Scan correctness** (`test_pact_scans`). Cells of both scans were
recomputed independently: firfilt or reflection padding arithmetic by
direct convolution, segment-wise Hilbert, guard, HAS and the statistics.
They agree with production to ≤ 1e-5 relative (single-precision storage)
and ≤ 1e-9 respectively. The coupled channel's maximum falls in the cell
containing 3 Hz × 45 Hz. Inputs are unchanged. Output dimensions are
C × (bins−1) × (bins−1). 100% HAS uses every valid sample. Mixed
numeric/text boundaries are handled. Too-short data and band edges above
Nyquist/1.2 raise `PACT:ScanUnsupported`.

**What the legacy scan reported.** On a 120 s fixture coupled at 3 Hz,
master's HAS × LFO scan gave MI ≈ 1.4 in the 4.16–6 Hz row (no 3 Hz
content), because the 20 dB filter leaked the 3 Hz rhythm into it. The
revised scan gives 0.04–0.24 there and keeps ≈ 1.4–1.6 in the 2.88–4.16 Hz
row.

**Runtime** (timer handles). 8 ch × 120 s @ 200 Hz: HFO × phase 10 × 10
grid 0.7 s; HAS × LFO 15 × 15 grid 3.0 s. `pac_man` (2% HAS, 2000
permutations, 6 boundaries): 64 ch × 10 min @ 1 kHz 11.4 s; 32 ch × 1 h
@ 500 Hz 27.3 s; 128 ch × 20 min @ 500 Hz 19.8 s; 21 ch × 5 min @ 200 Hz
0.3 s; ≤ 3.6 GB.

## 6. Old vs new numerical results (120 s fixture, 4 channels)

| Output | master `7d39d9c` → round 1 `07931e4` → round 2 |
| --- | --- |
| `mi`, `vectLength`, `angleMean`, `angleStd` | Identical (≤ 7e-18) across all three. |
| Continuous p-values | master permutation null [0, 0.030, 0.80, 0.33] → [0.005, 0.28, 0.985, 0.765] in round 1 and round 2. Round 2 is identical to round 1 with more than 16 occupied blocks; with ≤ 16 blocks it is now exact. |
| `dmi` | Identical between round 1 and round 2. |
| Event-mode p-values | Round 1 block test [0.005, 0.54, 0.12, 0.135] → round 2 trial shuffle [0.005, 0.685, 0.15, 0.075] (intentional change of null). |
| `pac_scanHfoByPhase` | Up to 2.2% relative change on this fixture, mainly in the lowest bins (revised filters and guard). |
| `pac_scanLfoPhaseFreq` | Up to 80% relative change: spurious cross-band values removed (§5). |

## 7. Bugs found and fixed in round 2

1. `pac_blockSurrogate`: when dmi is identically zero (constant amplitude at 100% HAS), rounding residue was compared with a null of equally tiny values, so p was arbitrary and could fall below 0.05. Now dmi = 0 and p = 1.
2. With few occupied blocks the Monte Carlo null was used instead of the exact null, and an untestable configuration (min p > α) gave silent 0% rejection; now exact, and `PACT:Underpowered` warns.
3. Event mode used the block test, which flags evoked responses (12–14% rejections on evoked-only data); it now uses a validated trial-shuffle null.
4. `pac_scanHfoByPhase`: collapsed passbands (up to −3.4 dB at the centre), a 1 s guard shorter than the filter transients, Hilbert transforms across boundaries, and amplitude envelopes recomputed for every phase band.
5. `filt_fir1fft`: violated its documented ripple spec (rectangular window), leaking out-of-band rhythms into neighbouring scan rows.
6. `pac_scanLfoPhaseFreq`: no edge guard; numeric −99 boundaries not refused; band edges above Nyquist/1.2 failed with an obscure error.
7. Benchmarks: firfilt's internal `tic` made `tic/toc` timings about 20% too low.

## 8. Remaining assumptions and limitations

* All p-values are exploratory. The continuous test assumes independent,
  centrally symmetric block contributions; the event test assumes
  exchangeable, non-overlapping events. In the holdout configuration (500 Hz, 6 Hz LFO, 120 s), large replications show a slight liberal excess of about +0.5 percentage points (5.4–5.7% at nominal 5%) with no identified mechanism. Treat p-values near 0.05 accordingly.
* Rejections mean statistical dependence, not physiological PAC.
  Sharp-edged waveforms and common-input transients produce genuine
  dependence (§4.5). Spectral or waveform-shape diagnostics are still the
  user's responsibility.
* Consistent preferred phase is required for power. A preferred phase that
  rotates over time (P03) cancels across blocks (§4.5).
* Correction is across channels only, not across scanned frequencies, HAS
  rates, block lengths, windows or conditions.
* With few blocks (short recordings, long `blockCycles`, sparse HAS) or few
  events, p-values are coarse or NaN; check `minAttainableP`.
* Scans are descriptive: no p-values, and adjacent bands overlap. Narrow
  low-frequency bins need long filters, and the guard removes their
  transient support from both ends and around boundaries.
* Only MATLAB R2024b + EEGLAB 2024.0 (FIRfilt 2.8) were tested. MATLAB
  R2022b is installed on this machine but was not tested.

## 9. Round 1 (documented, superseded where noted)

Round 1 (`07931e4`) replaced the block-rotation null after it failed with
non-sinusoidal LFOs (20/91/100/100% at HAS 1/5/20/100%). It moved
CircStat to `external/` and diagnosed the 18 Nov menu-test hang. That hang
was MATLAB `inputdlg`'s `'Cancel'` tag vs the harness's `'cancel'` search,
leaving an indefinite `uiwait`. It was prevented by a static guard,
real-dialog deadlines and a process-level timeout. Round 1's 800-channel
grid (N1–N6, C1; 1.7–5.9% Type I) remains valid for the continuous test.
Round 2's exact null changes p only when ≤ 16 blocks are occupied.

## 10. Output fields

C = channels, T = samples, B = phase bins, W = accepted events, L = event
window length (odd).

| Field | Size and meaning |
| --- | --- |
| `analyticEEG` | C × T. HFO envelope × unit phasor of LFO phase, in original coordinates. Excluded samples remain; use the masks. |
| `validSampleMask`, `selectedSampleMask` | 1 × T. Valid excludes filter guards at edges and boundaries (≥ 1 s or half the FIR order). Selected is valid ∩ event windows in event mode. |
| `hfoIndex` | C × 1 cell of one-based original sample indices above the HAS threshold (all eligible samples at 100%). |
| `mi` | C × 1. Historical MI, descriptive. Tested only in event mode. |
| `dmi` | C × 1. Block-debiased MI; tested in continuous modes. |
| `testedStatistic` | 'dmi' (continuous modes) or 'mi' (event mode): the statistic that the p-values and null percentiles refer to. |
| `numBlocks`, `blockLength`, `blockCycles` | C × 1 occupied blocks; block length in samples; block length in LFO cycles. |
| `numEvents` | Number of event windows used by the trial shuffle (event mode). |
| `exactNull`, `minAttainableP` | C × 1. Whether the null was enumerated, and the smallest p the null allows. |
| `vectLength`, `vectLengthVar`, `angleMean`, `angleStd` | C × 1. Resultant length, 1 − length, mean phase in [0, 2π), angular deviation √(2(1−R)) (CircStat `circ_std`). |
| `moduInd95CI`, `moduInd99CI` | C × 1. 95th/99th percentiles of the null of the tested statistic. These are historical names, not confidence intervals. |
| `phaseBinValues`, `phaseProbability`, `phaseSortedAmp`, `phaseSortedAmpSe`, `histBinMax` | Descriptive phase-bin summaries ordered 0 to 2π. The SEs ignore autocorrelation. |
| `uncorrected`, `Bonferroni`, `BonfHolm`, `FDR` | Structs of C × 1 p-values. `MIpval` refers to `testedStatistic`. Legacy distribution-test fields are NaN. |
| `winEdgesAfterRejection` | W × 2 inclusive window edges, centred on events. |
| `windowMeanAmp`, `windowStdAmp`, `windowMeanPhase`, `windowStdPhase` | C × L. Event-window summaries (circular mean and angular deviation for phase). |
| `surrogateMethod`, `inferenceNotice`, `rngSeed` | Method description, caveats, and the fixed seed (0). The caller's RNG is restored. |
| `pacScan.*` (scans) | Adds `phaseFilterOrder`, `ampFilterOrder` / `filterOrders`, `guardSamples`, `validSampleMask`, `hasRates`, `dataDimensions`. |

### Historical fixture

```powershell
git show 7d39d9c1fdb03b4c72635e2360f0c6dc7e834c8b:pac_man.m |
    ForEach-Object { $_ -replace '^function EEG = pac_man', 'function EEG = pact_legacy_fixture' } |
    Set-Content (Join-Path $env:TEMP 'pact_legacy_fixture.m')
```
