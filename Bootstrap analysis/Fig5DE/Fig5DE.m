% Fig5DE: run this file from MATLAB.
% Original analysis: Bootstrap/fig 5de/fig5de_paired_bagged_MET.m
% MATLAB R2023a; Statistics and Machine Learning Toolbox.
% All paths are relative to this file; input files are read only.
% Output files are written under results/ with the effective configuration.
% Run/F5 uses the script path; F9 uses the active saved Editor file.
figureFile = string(mfilename('fullpath'));
if strlength(figureFile) == 0 && usejava('desktop')
    figureFile = string(matlab.desktop.editor.getActiveFilename);
end
[figureDir, figureName] = fileparts(figureFile);
if ~strcmp(figureName, 'Fig5DE')
    error('bootstrap:ScriptLocation', ...
        'Cannot locate Fig5DE.m. Run the full saved file (Run/F5 or run(fullPath)).');
end
scriptDir = fullfile(fileparts(figureDir), 'script');
if ~isfile(fullfile(figureDir, 'data.xlsx')) || ~isfolder(scriptDir)
    error('bootstrap:FigureLayout', ...
        'Expected data.xlsx in %s and a sibling script folder.', figureDir);
end
addpath(scriptDir);
% Sheet names are editable; keep sheet order and data-column positions.

cfg = struct();
cfg.outputDir = fullfile(figureDir, 'results');
cfg.outputStem = "Fig5DE";
cfg.showFigures = true;
% Analysis settings: these affect the estimator or statistical inference.
cfg.nInner = 1000;
cfg.nOuter = 2000;
cfg.rngSeed = 1;
cfg.delta0 = 1; % hours
cfg.alpha = 0.05;
cfg.dataFolder = figureDir;
cfg.inputFile = fullfile(figureDir, 'data.xlsx');
cfg.controlPath = cfg.inputFile;
cfg.experimentPaths = repmat(cfg.inputFile, 2, 1);
cfg.controlMSheet = struct('sheet',1,'range',"A1:U19");
cfg.controlESheet = struct('sheet',1,'range',"W1:AQ19");
cfg.experimentMSheets = [struct('sheet',2,'range',"A1:U13"); struct('sheet',3,'range',"A1:U13")];
cfg.experimentESheets = [struct('sheet',2,'range',"W1:AQ13"); struct('sheet',3,'range',"W1:AQ13")];
cfg.experimentTags = ["MPDFR"; "EPDFR"];
cfg.experimentLabels = ["MPDFR (experiment)"; "EPDFR (experiment)"];
cfg.controlTag = "LD1212_LD4";
cfg.controlLabel = "LD1212_LD4 (control)";
cfg.M_sheet = "sLNv";
cfg.E_sheet = "5th sLNv";
cfg.multiplicityAdjustment = "None; 4 raw MET p-values are reported";
cfg.M_color = [1, 128/255, 0];
cfg.E_color = [0, 0, 1];
cfg.controlColor = [128, 128, 128] / 255;
cfg.mainFigureBinWidth = 1;
cfg.plot.axesSizeCm = [1.1, 1.2];
% Each cell-type contrast is experiment peak minus control peak.
results = runPeakShiftComparison(cfg);
