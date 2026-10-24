% figS5E: run this file from MATLAB.
% Original analysis: Bootstrap/fig s5d_per0/fig_s5d_per0_main.m
% MATLAB R2023a; Statistics and Machine Learning Toolbox.
% All paths are relative to this file; input files are read only.
% Output files are written under results/ with the effective configuration.
% Run/F5 uses the script path; F9 uses the active saved Editor file.
figureFile = string(mfilename('fullpath'));
if strlength(figureFile) == 0 && usejava('desktop')
    figureFile = string(matlab.desktop.editor.getActiveFilename);
end
[figureDir, figureName] = fileparts(figureFile);
if ~strcmp(figureName, 'figS5E')
    error('bootstrap:ScriptLocation', ...
        'Cannot locate figS5E.m. Run the full saved file (Run/F5 or run(fullPath)).');
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
cfg.outputStem = "figS5E";
cfg.showFigures = true;
cfg.perFile = fullfile(figureDir, 'data.xlsx');
cfg.wtFile = cfg.perFile;
cfg.sheetMap.per0 = struct('m',struct('sheet',1,'range',"A2:D280"),'e',struct('sheet',1,'range',"F2:I280"));
cfg.sheetMap.wt = struct('m',struct('sheet',2,'range',"A2:D380"),'e',struct('sheet',2,'range',"F2:I379"));
cfg.figureLabel = "figS5E";
cfg.modelName = "periodic_spline";
cfg.period = 24;
cfg.nSplineBasis = 5;
cfg.harmonicOrder = 2;
cfg.gridStep = 0.1;
cfg.nNullBootstrap = 2000;
cfg.nAmplitudeBootstrap = 2000;
cfg.nullMode = "pooled_bootstrap";
cfg.randomSeed = 20260728;
cfg.alpha = 0.05;
cfg.progressEvery = 250;
cfg.verbose = true;
cfg.maxDenominator = 12;
cfg.denominatorTolerance = 1e-8;
% Preserve the original data-relative margin and record how it is defined.
cfg.equivalenceFractionOfWT = 0.20;
cfg.equivalenceAbsolute = struct('m', nan, 'e', nan);
cfg.plotRawPoints = true;
cfg.plotAmplitudeScatter = true;

results = runClockBootstrapAnalysis(cfg);
