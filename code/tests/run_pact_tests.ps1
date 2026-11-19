<#
.SYNOPSIS
  Run the PACT test suites, each in its own time-limited MATLAB process.

.DESCRIPTION
  1. Mocked suite (no real dialogs): math, CircStat, references, integration,
     filter designs, independent scan recomputation, static modal-dialog
     guard, and mocked menu callbacks.
  2. Real-GUI suite: dialog-timeout regression and real-dialog smoke test.
  3. Process-level regression (optional, about 2 minutes): the original
     Cancel-tag hang is reproduced and must be terminated by the wrapper.
  4. Scientific validation (optional, about 1 hour): event-mode validation,
     filter-response audit, and the null grid in quick mode. The full null
     grid (validate_pact_null, validate_pact_null('holdout')) takes several
     hours; run it separately.
  Every MATLAB run goes through gui\run_matlab_with_timeout.ps1.

.EXAMPLE
  powershell -NoProfile -ExecutionPolicy Bypass -File code\tests\run_pact_tests.ps1 `
     -Eeglab C:/path/to/eeglab2024.0 -IncludeProcessRegression -IncludeScience

  Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.
#>
param(
    [Parameter(Mandatory = $true)][string]$Eeglab,
    [switch]$IncludeProcessRegression,
    [switch]$IncludeScience,
    [int]$TimeoutSeconds = 900
)
$tests = Split-Path -Parent $MyInvocation.MyCommand.Path
$root = (Resolve-Path (Join-Path $tests '..\..')).Path -replace '\\', '/'
$eeg = $Eeglab -replace '\\', '/'
$wrapper = Join-Path $tests 'gui\run_matlab_with_timeout.ps1'
$setup = "addpath('$eeg'); evalc('eeglab nogui'); addpath('$root'); pac_setPath; disp(version); disp(eeg_getversion);"
$suites = [ordered]@{
    mocked = "$setup cd('$root/code/tests'); test_pact_math; test_pact_circstat; test_pact_reference; test_pact_integration; test_pact_filters; test_pact_scans; test_pact_dialog_guard; test_pact_menus;"
    gui    = "$setup test_pact_dialog_timeout; smoke_pact_real_dialogs;"
}
if ($IncludeScience) {
    $suites['science'] = "$setup cd('$root/code/tests'); validate_pact_events; validate_pact_filters; validate_pact_null('quick');"
}

$failed = @()
foreach ($name in $suites.Keys) {
    $log = Join-Path $env:TEMP "pact-tests-$name.log"
    $limit = $TimeoutSeconds
    if ($name -eq 'science') { $limit = [Math]::Max($TimeoutSeconds, 10800) }
    & powershell -NoProfile -ExecutionPolicy Bypass -File $wrapper -Command $suites[$name] `
        -TimeoutSeconds $limit -LogPath $log | Write-Output
    $code = $LASTEXITCODE
    Select-String -Path $log -CaseSensitive -Pattern 'PASS|FAIL|RESULT|BENCH|LEGACY' | ForEach-Object { '  ' + $_.Line }
    if ($code -ne 0) { $failed += "$name (exit $code; see $log)" }
}
if ($IncludeProcessRegression) {
    & powershell -NoProfile -ExecutionPolicy Bypass -File (Join-Path $tests 'gui\test_process_timeout.ps1') | Write-Output
    if ($LASTEXITCODE -ne 0) { $failed += 'process timeout regression' }
}
if ($failed) { Write-Output ('FAILED: ' + ($failed -join '; ')); exit 1 }
Write-Output 'ALL PACT TEST SUITES PASSED'
exit 0
