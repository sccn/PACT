<#
.SYNOPSIS
  Run a MATLAB -batch command with an independent wall-clock limit.

.DESCRIPTION
  Starts MATLAB in its own process, waits up to TimeoutSeconds, and on expiry
  terminates the whole process tree it started (matlab.exe launcher,
  MATLAB.exe and any child workers) with taskkill /T /F. Only that tree is
  touched; other MATLAB sessions are left alone. This guard does not rely on
  MATLAB timers or callbacks, so it also bounds hangs that MATLAB cannot
  recover from (modal waits, deadlocks, crashes that leave a process alive).

  Exit codes: the MATLAB exit code, or 124 on timeout (the coreutils
  convention), or 2 for invalid arguments.

.EXAMPLE
  powershell -NoProfile -ExecutionPolicy Bypass -File code\tests\gui\run_matlab_with_timeout.ps1 `
    -TimeoutSeconds 600 -LogPath $env:TEMP\pact-gui.log `
    -Command "addpath('C:/path/to/eeglab'); eeglab nogui; addpath(pwd); pac_setPath; smoke_pact_real_dialogs"

  Copyright (C) 2026 Makoto Miyakoshi. GPL-2.0-or-later.
#>
param(
    [Parameter(Mandatory = $true)][string]$Command,
    [int]$TimeoutSeconds = 600,
    [string]$WorkingDirectory = (Get-Location).Path,
    [string]$LogPath = (Join-Path $env:TEMP ("pact-matlab-{0}.log" -f [guid]::NewGuid())),
    [string]$Matlab = ''
)

if ($Command.Contains('"')) {
    Write-Error 'Use single quotes inside the MATLAB command; double quotes cannot be passed safely.'
    exit 2
}
if ($TimeoutSeconds -lt 1) { Write-Error 'TimeoutSeconds must be positive.'; exit 2 }
if (-not $Matlab) {
    $found = Get-Command matlab -ErrorAction SilentlyContinue
    if (-not $found) { Write-Error 'matlab not found on PATH; pass -Matlab.'; exit 2 }
    $Matlab = $found.Source
}

$errorLog = "$LogPath.stderr"
$started = Get-Date
$process = Start-Process -FilePath $Matlab -ArgumentList ('-batch "{0}"' -f $Command) `
    -WorkingDirectory $WorkingDirectory -NoNewWindow -PassThru `
    -RedirectStandardOutput $LogPath -RedirectStandardError $errorLog
$null = $process.Handle   # cache the handle so ExitCode is available later

if ($process.WaitForExit($TimeoutSeconds * 1000)) {
    $process.WaitForExit()
    $code = $process.ExitCode
    Write-Output ("MATLAB finished in {0:N1} s with exit code {1}. Log: {2}" -f `
        ((Get-Date) - $started).TotalSeconds, $code, $LogPath)
    exit $code
}

# Timeout: kill only the tree rooted at the process we started.
& taskkill.exe /T /F /PID $process.Id 2>&1 | Out-Null
$process.WaitForExit(30000) | Out-Null
Write-Output ("TIMEOUT after {0} s: terminated MATLAB process tree rooted at PID {1}. Log: {2}" -f `
    $TimeoutSeconds, $process.Id, $LogPath)
exit 124
