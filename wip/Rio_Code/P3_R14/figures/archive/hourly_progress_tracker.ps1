$ErrorActionPreference = "Continue"

$sourceRepo = "D:/LRHWork/Model_RH/PyBaMM-ECdrag2"
$portRepo = "D:/LRHWork/Model_RH/PyBaMM-port-work"
$logFile = "D:/LRHWork/Model_RH/PyBaMM-ECdrag2/wip/Rio_Code/P3_R14/hourly_progress.log"

function Write-Section {
    param([string]$title)
    Add-Content -Path $logFile -Value ""
    Add-Content -Path $logFile -Value ("==== " + $title + " ====")
}

while ($true) {
    $ts = Get-Date -Format "yyyy-MM-dd HH:mm:ss"
    Add-Content -Path $logFile -Value ""
    Add-Content -Path $logFile -Value ("[Tick] " + $ts)

    try {
        Write-Section "Port repo status"
        $portStatus = git -C $portRepo status --short --branch 2>&1
        $portStatus | ForEach-Object { Add-Content -Path $logFile -Value $_ }

        Write-Section "Port repo recent commits"
        $portLog = git -C $portRepo log --oneline -n 3 2>&1
        $portLog | ForEach-Object { Add-Content -Path $logFile -Value $_ }
    } catch {
        Add-Content -Path $logFile -Value ("[Error][Port repo] " + $_.Exception.Message)
    }

    try {
        Write-Section "Source repo status"
        $srcStatus = git -C $sourceRepo status --short --branch 2>&1
        $srcStatus | ForEach-Object { Add-Content -Path $logFile -Value $_ }
    } catch {
        Add-Content -Path $logFile -Value ("[Error][Source repo] " + $_.Exception.Message)
    }

    Add-Content -Path $logFile -Value "----------------------------------------"
    Start-Sleep -Seconds 3600
}

