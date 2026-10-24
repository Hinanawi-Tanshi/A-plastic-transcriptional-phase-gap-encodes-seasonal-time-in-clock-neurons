% figS5C: run this file from MATLAB.
% Original analysis: Bootstrap/fig s5d/figs5d_paired_bagged_MET.m
% MATLAB R2023a; Statistics and Machine Learning Toolbox.
% All paths are relative to this file; input files are read only.
% Output files are written under results/ with the effective configuration.
% Run/F5 uses the script path; F9 uses the active saved Editor file.
figureFile = string(mfilename('fullpath'));
if strlength(figureFile) == 0 && usejava('desktop')
    figureFile = string(matlab.desktop.editor.getActiveFilename);
end
[figureDir, figureName] = fileparts(figureFile);
if ~strcmp(figureName, 'figS5C')
    error('bootstrap:ScriptLocation', ...
        'Cannot locate figS5C.m. Run the full saved file (Run/F5 or run(fullPath)).');
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
cfg.outputStem = "figS5C";
cfg.showFigures = true;
% Analysis settings: these affect the estimator or statistical inference.
cfg.nInner = 1000;
cfg.nOuter = 2000;
cfg.rngSeed = 1;
cfg.delta0 = 1; % hours
cfg.alpha = 0.05;
cfg.inputFile = fullfile(figureDir, 'data.xlsx');
cfg.dataSheet = 1;
cfg.timeColumn = 2;
cfg.MColumn = 4;
cfg.EColumn = 5;
% The reported single-group difference is M minus E, in hours.
% Plot settings below do not change the bootstrap estimator.
cfg.mainFigureBinWidth = 1;
cfg.plot.figureSizeCm = [8, 6];
cfg.plot.axesSizeCm = [3.85, 2.5];
cfg.plot.fontSize = 8;
cfg.plot.lineWidth = 0.8;
cfg.plot.autoYLimits = false;
cfg.plot.xLimits = [10, 20];
cfg.plot.xTicks = 10:2:20;
cfg.plot.minorTicks = false;

results = runPairedGapAnalysis(cfg);
