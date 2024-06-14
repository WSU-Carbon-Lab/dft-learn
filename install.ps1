# INSTALL SCRIPT FOR THE CLUSTERING ALGORITHM
# -------------------------------------------

# Ensure that there is a WaveMetrics folder in the users documents folder
$folder = [Environment]::GetFolderPath("MyDocuments") + "\WaveMetrics"
if (-not (Test-Path $folder)) {
    # Path does not exist so exit
    Write-Host "The WaveMetrics folder does not exist in the users documents folder. Please create the folder and try again."
    exit
}

# If there are multiple versions of igor installed then the user will be prompted to select the version to install to
$igorVersions = Get-ChildItem "C:\Program Files\WaveMetrics\Igor Pro*" | Where-Object { $_.PSIsContainer } | Select-Object Name
if ($igorVersions.Count -gt 1) {
    Write-Host "Multiple versions of Igor Pro have been detected. Please select the version to install to:"
    $igorVersions | ForEach-Object { Write-Host $_.Name }
    $igorVersion = Read-Host "Enter the version to install to"
    $igorVersion = $igorVersion -replace "Igor Pro ", ""
    $igorVersion = $igorVersion -replace "Folder", ""
    $igorVersion = $igorVersion -replace "\.", ""
    $igorVersion = $igorVersion -replace " ", ""
} else {
    $igorVersion = $igorVersions.Name -replace "Igor Pro ", ""
    $igorVersion = $igorVersion -replace "Folder", ""
    $igorVersion = $igorVersion -replace "\.", ""
    $igorVersion = $igorVersion -replace " ", ""
}


# Link / Copy the clusteringPanel to the Igor Procs Folder
# Link / Copy the cluseringAlgorithm to the User Procs Folder
$igorProcsFolder = "$folder\Igor Pro $igorVersion User Files\Igor Procedures"
$userProcsFolder = "$folder\Igor Pro $igorVersion User Files\User Procedures"
$currentDirectory = Get-Location
# Link the clusteringPanel to the Igor Procs Folder
New-Item -ItemType SymbolicLink -Path "$igorProcsFolder\clusteringPanel.ipf" -Value "$currentDirectory\clusteringPanel v1.ipf" -Force
# Link the clusteringAlgorithm to the Igor Procs Folder
New-Item -ItemType SymbolicLink -Path "$userProcsFolder\DFT_Clustering" -Value "$currentDirectory\DFT_Clustering" -Force
# Link the Element Library to the User Procs Folder
New-Item -ItemType SymbolicLink -Path "$userProcsFolder\Element Library" -Value "$currentDirectory\Element Library" -Force
