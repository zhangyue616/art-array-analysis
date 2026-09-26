param(
    [string]$PythonExe = 'python',
    [int]$Threads = 8,
    [string]$Aria2 = '',
    [string]$Fastp = '',
    [string]$Bowtie2 = '',
    [string]$Bowtie2Index = '',
    [switch]$StopAfterFirstLibrary
)

$ErrorActionPreference = 'Stop'

$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot '../..')).Path
$statePath = Join-Path $projectRoot 'data/source_metadata/p6_rnaseq/orchestrator_state.json'
$scriptPath = Join-Path $projectRoot 'scripts/p6_rnaseq/orchestrate_complete_runs.py'
$stdoutPath = Join-Path $projectRoot 'data/source_metadata/p6_rnaseq/orchestrator_process.stdout.log'
$stderrPath = Join-Path $projectRoot 'data/source_metadata/p6_rnaseq/orchestrator_process.stderr.log'

if (Test-Path $statePath) {
    $state = Get-Content -Raw -LiteralPath $statePath | ConvertFrom-Json
    if ($state.stage -eq 'all_runs_complete') {
        Write-Output 'All 12 runs are already marked complete; no orchestrator started.'
        exit 0
    }
    if ($state.orchestrator_pid) {
        $activeProcess = Get-Process -Id ([int]$state.orchestrator_pid) -ErrorAction SilentlyContinue
        if ($activeProcess -and $activeProcess.ProcessName -eq 'python') {
            Write-Output ("Orchestrator is already active as PID {0}; no duplicate started." -f $activeProcess.Id)
            exit 0
        }
    }
}

$orchestratorArguments = @($scriptPath, '--threads', [string]$Threads)
if ($Aria2) { $orchestratorArguments += @('--aria2', $Aria2) }
if ($Fastp) { $orchestratorArguments += @('--fastp', $Fastp) }
if ($Bowtie2) { $orchestratorArguments += @('--bowtie2', $Bowtie2) }
if ($Bowtie2Index) { $orchestratorArguments += @('--bowtie2-index', $Bowtie2Index) }
if ($StopAfterFirstLibrary) { $orchestratorArguments += '--stop-after-first-library' }

$newProcess = Start-Process `
    -FilePath $PythonExe `
    -ArgumentList $orchestratorArguments `
    -WorkingDirectory $projectRoot `
    -WindowStyle Hidden `
    -RedirectStandardOutput $stdoutPath `
    -RedirectStandardError $stderrPath `
    -PassThru

Write-Output ("Started resumable P6 RNA-seq orchestrator as PID {0}." -f $newProcess.Id)
