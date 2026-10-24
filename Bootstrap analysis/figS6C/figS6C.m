% figS6C: run this file from MATLAB.
% Original analysis: Bootstrap/fig s6c/figs6c_paired_bagged_MET.m
% MATLAB R2023a; Statistics and Machine Learning Toolbox.
% All paths are relative to this file; input files are read only.
% Output files are written under results/ with the effective configuration.
% Run/F5 uses the script path; F9 uses the active saved Editor file.
figureFile = string(mfilename('fullpath'));
if strlength(figureFile) == 0 && usejava('desktop')
    figureFile = string(matlab.desktop.editor.getActiveFilename);
end
[figureDir, figureName] = fileparts(figureFile);
if ~strcmp(figureName, 'figS6C')
    error('bootstrap:ScriptLocation', ...
        'Cannot locate figS6C.m. Run the full saved file (Run/F5 or run(fullPath)).');
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
cfg.outputStem = "figS6C";
cfg.showFigures = true;
% Analysis settings: these affect the estimator or statistical inference.
cfg.nInner = 1000;
cfg.nOuter = 2000;
cfg.rngSeed = 1;
cfg.delta0 = 0.5; % hours
cfg.alpha = 0.05;
cfg.dataSheet = "data";
cfg.timeColumn = 2;
cfg.MColumn = 4;
cfg.EColumn = 5;
% The reported single-group difference is M minus E, in hours.
% Plot settings below do not change the bootstrap estimator.
cfg.mainFigureBinWidth = 1;
cfg.plot.figureSizeCm = [5, 4];
cfg.plot.axesSizeCm = [2.6, 2];
cfg.plot.fontSize = 7;
cfg.plot.lineWidth = 0.8;
cfg.plot.autoYLimits = false;
cfg.plot.xLimits = [16, 24];
cfg.plot.xTicks = 16:2:24;
cfg.plot.minorTicks = true;

compartments = ["wholecell", "nuclear"];
dataSheets = [1, 2];
results = struct();
for k = 1:numel(compartments)
    cfg.inputFile = fullfile(figureDir, 'data.xlsx');
    cfg.dataSheet = dataSheets(k);
    cfg.outputStem = "figS6C_" + compartments(k);
    results.(compartments(k)) = runPairedGapAnalysis(cfg);
end
