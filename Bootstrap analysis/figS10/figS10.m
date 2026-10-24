% figS10: run this file from MATLAB.
% Original analysis: Bootstrap/fig s10/figs10_paired_bagged_MET.m
% MATLAB R2023a; Statistics and Machine Learning Toolbox.
% All paths are relative to this file; input files are read only.
% Output files are written under results/ with the effective configuration.
% Run/F5 uses the script path; F9 uses the active saved Editor file.
figureFile = string(mfilename('fullpath'));
if strlength(figureFile) == 0 && usejava('desktop')
    figureFile = string(matlab.desktop.editor.getActiveFilename);
end
[figureDir, figureName] = fileparts(figureFile);
if ~strcmp(figureName, 'figS10')
    error('bootstrap:ScriptLocation', ...
        'Cannot locate figS10.m. Run the full saved file (Run/F5 or run(fullPath)).');
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
cfg.outputStem = "figS10";
cfg.showFigures = true;
% Analysis settings: these affect the estimator or statistical inference.
cfg.nInner = 1000;
cfg.nOuter = 2000;
cfg.rngSeed = 1;
cfg.delta0 = 1; % hours
cfg.alpha = 0.05;
cfg.showFigures = false;
cfg.dataFolder = figureDir;
cfg.inputFile = fullfile(figureDir, 'data.xlsx');
cfg.datasetNames = ["DL4"; "DL6"; "DL9"; "DL10"; "DL11"; "DL12"; "DL13"; "DL14"; "DL15"; "DL18"; "DL20"];
cfg.MSheets = [struct('sheet',1,'range',"A2:U14"); struct('sheet',2,'range',"A2:U17"); struct('sheet',3,'range',"A2:U15"); struct('sheet',4,'range',"A2:U16"); struct('sheet',5,'range',"A2:U15"); struct('sheet',6,'range',"A2:U19"); struct('sheet',7,'range',"A2:U15"); struct('sheet',8,'range',"A2:U15"); struct('sheet',9,'range',"A2:U15"); struct('sheet',10,'range',"A2:U17"); struct('sheet',11,'range',"A2:U14")];
cfg.ESheets = [struct('sheet',1,'range',"W2:AQ14"); struct('sheet',2,'range',"W2:AQ17"); struct('sheet',3,'range',"W2:AQ15"); struct('sheet',4,'range',"W2:AQ16"); struct('sheet',5,'range',"W2:AQ15"); struct('sheet',6,'range',"W2:AQ19"); struct('sheet',7,'range',"W2:AQ15"); struct('sheet',8,'range',"W2:AQ15"); struct('sheet',9,'range',"W2:AQ15"); struct('sheet',10,'range',"W2:AQ17"); struct('sheet',11,'range',"W2:AQ14")];
cfg.filePattern = 'DL*.xlsx';
cfg.rngSeedBase = 1;
cfg.meanSEMCondition = "ZT";
cfg.meanSEMUnit = "hemibrains";
cfg.mainFigureBinWidth = 1;
cfg.diagnosticDiffBinWidth = 0.5;
cfg.mainXLimits = [6, 22];
% Each file: paired within a time row; rows are independently resampled.
results = runPairedBatchAnalysis(cfg);
