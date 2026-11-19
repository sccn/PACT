<#
.SYNOPSIS
  Regression: the original Cancel-tag hang can no longer wait indefinitely.

.DESCRIPTION
  Runs repro_original_cancel_tag_hang (the original harness logic against a
  real inputdlg, which never returns on its own) under
  run_matlab_with_timeout.ps1 and checks that:
    * the repro actually reached the blocking dialog,
    * the wrapper reported a timeout (exit code 124) close to the limit,
    * no process from that run remains afterwards.
  Exit code 0 means PASS.

.EXAMPLE
  powershell -NoProfile -ExecutionPolicy Bypass -File code\tests\gui\test_process_timeout.ps1

  Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.
#>
param([int]$TimeoutSeconds = 90)

$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$log = Join-Path $env:TEMP ("pact-timeout-regression-{0}.log" -f [guid]::NewGuid())
$command = "addpath('{0}'); repro_original_cancel_tag_hang" -f ($here -replace '\\', '/')
$started = Get-Date
& powershell -NoProfile -ExecutionPolicy Bypass -File (Join-Path $here 'run_matlab_with_timeout.ps1') `
    -Command $command -TimeoutSeconds $TimeoutSeconds -LogPath $log | Write-Output
$code = $LASTEXITCODE
$elapsed = ((Get-Date) - $started).TotalSeconds

$failures = @()
if ($code -ne 124) { $failures += "expected timeout exit code 124, got $code" }
if ($elapsed -lt $TimeoutSeconds) { $failures += ("ended after {0:N0} s, before the limit" -f $elapsed) }
if ($elapsed -gt $TimeoutSeconds + 60) { $failures += ("took {0:N0} s, far beyond the limit" -f $elapsed) }
if (-not (Select-String -Path $log -Pattern 'REPRO waiting in inputdlg' -Quiet)) {
    $failures += 'repro never reached the blocking dialog'
}
if (Select-String -Path $log -Pattern 'REPRO unexpectedly returned' -Quiet) {
    $failures += 'repro returned; it no longer reproduces the hang'
}
Start-Sleep -Seconds 2
$left = Get-CimInstance Win32_Process | Where-Object { $_.CommandLine -like '*repro_original_cancel_tag_hang*' }
if ($left) { $failures += ('processes left running: ' + (($left | ForEach-Object ProcessId) -join ', ')) }

if ($failures) {
    Write-Output ('test_process_timeout: FAIL - ' + ($failures -join '; '))
    exit 1
}
Write-Output ("test_process_timeout: PASS (original hang reproduced and terminated after {0:N0} s)" -f $elapsed)
exit 0
