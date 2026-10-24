% Fig5H: two independent dbt phase-gap contrasts.
% Run/F5 uses the script path; F9 uses the active saved Editor file.
figureFile = string(mfilename('fullpath'));
if strlength(figureFile) == 0 && usejava('desktop')
    figureFile = string(matlab.desktop.editor.getActiveFilename);
end
[figureDir, figureName] = fileparts(figureFile);
if ~strcmp(figureName, 'Fig5H')
    error('bootstrap:ScriptLocation', ...
        'Cannot locate Fig5H.m. Run the full saved file (Run/F5 or run(fullPath)).');
end
scriptDir = fullfile(fileparts(figureDir), 'script');
if ~isfile(fullfile(figureDir, 'data.xlsx')) || ~isfolder(scriptDir)
    error('bootstrap:FigureLayout', ...
        'Expected data.xlsx in %s and a sibling script folder.', figureDir);
end
addpath(scriptDir);
% Sheet names are editable; keep sheet order and data-column positions.
results = struct();

cfg = struct();
cfg.outputDir = fullfile(figureDir, 'results');
cfg.outputStem = "Fig5H_dbtS";
cfg.showFigures = true;
% Analysis settings: these affect the estimator or statistical inference.
cfg.nInner = 1000;
cfg.nOuter = 2000;
cfg.rngSeed = 1;
cfg.delta0 = 0.5; % hours
cfg.alpha = 0.05;
cfg.inputFile = fullfile(figureDir, 'data.xlsx');
cfg.controlPath = cfg.inputFile;
cfg.analysisMode = "Fig5H_dbtS";
cfg.conditionLabel = "dbtS (condition)";
cfg.controlLabel = "L12D12 (control)";
cfg.conditionColor = [0, 0.5019607843137255, 1];
cfg.controlColor = [128, 128, 128] / 255;
cfg.dataSheet = 1;
cfg.controlSheet = 2;
cfg.timeColumn = 2;
cfg.MColumn = 4;
cfg.EColumn = 5;
cfg.mainFigureBinWidth = 1;
% Contrast: (E-M in condition) minus (E-M in independent control).
results.dbtS = runGapChangeComparison(cfg);

cfg = struct();
cfg.outputDir = fullfile(figureDir, 'results');
cfg.outputStem = "Fig5H_dbtL";
cfg.showFigures = true;
% Analysis settings: these affect the estimator or statistical inference.
cfg.nInner = 1000;
cfg.nOuter = 2000;
cfg.rngSeed = 1;
cfg.delta0 = 0.5; % hours
cfg.alpha = 0.05;
cfg.inputFile = fullfile(figureDir, 'data.xlsx');
cfg.controlPath = cfg.inputFile;
cfg.analysisMode = "Fig5H_dbtL";
cfg.conditionLabel = "dbtL (condition)";
cfg.controlLabel = "L9D15 (control)";
cfg.conditionColor = [1, 0, 0];
cfg.controlColor = [128, 128, 128] / 255;
cfg.dataSheet = 3;
cfg.controlSheet = 4;
cfg.timeColumn = 2;
cfg.MColumn = 4;
cfg.EColumn = 5;
cfg.mainFigureBinWidth = 1;
% Contrast: (E-M in condition) minus (E-M in independent control).
results.dbtL = runGapChangeComparison(cfg);
