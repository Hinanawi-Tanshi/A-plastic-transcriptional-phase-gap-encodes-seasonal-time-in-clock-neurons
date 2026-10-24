% Fig2C: run this file from MATLAB.
% Original analysis: Bootstrap/fig 2c/fig2c_flypaired_bagged_MET.m; fig2c_intensity_sixpanel_analysis.m
% MATLAB R2023a; Statistics and Machine Learning Toolbox.
% All paths are relative to this file; input files are read only.
% Output files are written under results/ with the effective configuration.
% Run/F5 uses the script path; F9 uses the active saved Editor file.
figureFile = string(mfilename('fullpath'));
if strlength(figureFile) == 0 && usejava('desktop')
    figureFile = string(matlab.desktop.editor.getActiveFilename);
end
[figureDir, figureName] = fileparts(figureFile);
if ~strcmp(figureName, 'Fig2C')
    error('bootstrap:ScriptLocation', ...
        'Cannot locate Fig2C.m. Run the full saved file (Run/F5 or run(fullPath)).');
end
scriptDir = fullfile(fileparts(figureDir), 'script');
if ~all(isfile(fullfile(figureDir, ["data_binary.xlsx", "data_intensity.xlsx"]))) || ~isfolder(scriptDir)
    error('bootstrap:FigureLayout', ...
        'Expected data_binary.xlsx, data_intensity.xlsx in %s and a sibling script folder.', figureDir);
end
addpath(scriptDir);
% Recording tab names: YYMMDD FLYn, matched across the two workbooks.
% Keep binary tab order (batch membership) and time/value column order.

cfgKDE = struct();
cfgKDE.outputDir = fullfile(figureDir, 'results');
cfgKDE.outputStem = "Fig2C_KDE";
cfgKDE.showFigures = true;
% Analysis settings: these affect the estimator or statistical inference.
cfgKDE.nInner = 1000;
cfgKDE.nOuter = 2000;
cfgKDE.rngSeed = 1;
cfgKDE.delta0 = 1; % hours
cfgKDE.alpha = 0.05;
cfgKDE.inputFile = fullfile(figureDir, 'data_binary.xlsx');
cfgKDE.intensityFile = fullfile(figureDir, 'data_intensity.xlsx');
% Consecutive binary tabs belong to the two original experimental batches.
% Intensity tabs are paired by date/FLY key, independently of their order.
cfgKDE.batchSheetCounts = [38, 15];
cfgKDE.batchNames = ["Batch1"; "Batch2"];
% Nominal start CT is inferred by rounding the first valid timestamp (hours).
% Preserve CT24 as CT24; do not wrap the recording time modulo 24.
% Fly/workbook is the biological unit. Four M cells are pooled within fly.
cfgKDE.recordingWindowHours = 2;
cfgKDE.representativeTimeOffset = 1;
cfgKDE.kdeWindowMinutes = 20;
cfgKDE.mainFigureBinWidth = 1;
cfgKDE.diagnosticBinWidth = 0.5;

cfgIntensity = struct();
cfgIntensity.outputDir = fullfile(figureDir, 'results');
cfgIntensity.outputStem = "Fig2C_Intensity";
cfgIntensity.showFigures = true;
% Analysis settings: these affect the estimator or statistical inference.
cfgIntensity.nInner = 1000;
cfgIntensity.nOuter = 2000;
cfgIntensity.rngSeed = 11;
cfgIntensity.delta0 = 1; % hours
cfgIntensity.alpha = 0.05;
cfgIntensity.inputFile = cfgKDE.inputFile;
cfgIntensity.batchNames = cfgKDE.batchNames;
cfgIntensity.intensityFile = cfgKDE.intensityFile;
cfgIntensity.batchSheetCounts = cfgKDE.batchSheetCounts;
cfgIntensity.nCTBootstrap = 5000;
cfgIntensity.recordingWindowHours = 2;
cfgIntensity.representativeTimeOffset = 1;
cfgIntensity.timestampToleranceHours = 1e-6;
cfgIntensity.panelCRequestedCT = [11, 13, 15, 17, 19];
cfgIntensity.deltaAmplitude = 0.20;
% Display settings; both analyses report the M-minus-E peak difference.
cfgIntensity.outerDifferenceBinWidth = 0.5;
cfgIntensity.outerPeakBinWidth = 0.5;
cfgIntensity.panelAPeakBinWidth = 1;
cfgIntensity.drawPanelA = true;
cfgIntensity.drawPanelC = true;
cfgIntensity.drawSixPanel = true;
cfgIntensity.drawDiagnostics = true;
cfgIntensity.panelCUseBoxplot = false;
cfgIntensity.M_color = [1, 128/255, 0];
cfgIntensity.E_color = [0, 0, 1];
cfgIntensity.diffColor = [0.45, 0.15, 0.55];

results = struct();
results.KDE = runFlyPairedKDEAnalysis(cfgKDE);
results.intensity = runFlyPairedIntensityAnalysis(cfgIntensity);
