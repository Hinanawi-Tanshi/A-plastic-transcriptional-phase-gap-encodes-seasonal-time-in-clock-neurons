% Fig1H: run the saved file with Run/F5 or run(fullPath).
% MATLAB R2023a. Inputs are read only; outputs go into results/.
% F9 also works when this saved file is the active MATLAB Editor document.
figureFile = string(mfilename('fullpath'));
if strlength(figureFile) == 0 && usejava('desktop')
    figureFile = string(matlab.desktop.editor.getActiveFilename);
end
[figureDir, figureName] = fileparts(figureFile);
if ~strcmp(figureName, 'Fig1H')
    error('kde:ScriptLocation', ...
        'Cannot locate Fig1H.m. Run the full saved file (Run/F5 or run(fullPath)).');
end
scriptDir = fullfile(fileparts(figureDir), 'script');
if ~isfile(fullfile(figureDir, 'data_binary.xlsx')) || ~isfolder(scriptDir)
    error('kde:FigureLayout', 'Expected data_binary.xlsx beside this file and a sibling script folder.');
end
addpath(scriptDir);

cfg = struct();
cfg.inputFile = fullfile(figureDir, 'data_binary.xlsx');
cfg.traceLayout = "fly_cells";
% Each recording keeps its native time/value pairs. Retain tab/column order.
cfg.cellType = "lLNv";
% Nominal ZT is rounded from each recording's first valid timestamp.
cfg.outputDir = fullfile(figureDir, 'results');
cfg.outputStem = "Fig1H";
cfg.showFigures = true;
cfg.exportPNG = true;
% Method7/7b only. KDE is recomputed from time/signal inputs.
% Cell plot: one dot per cell, flies sorted within condition, pooled-cell SEM.
% Fly plot and Kruskal-Wallis: one value per animal, SEM across animals.
% Requires Statistics and Machine Learning Toolbox for kruskalwallis.
cfg.conditions = ["LD4 ZT4-6","LD4 ZT8-10","LD4 ZT12-14"];
cfg.kdeWindowMinutes = 20;
% Raw timestamps include duplicates; processing preserves the producer's rules.
cfg.rngSeed = 42;
cfg.flyJitterWidth = 0.15;
cfg.cellJitterWidth = 0.12;
cfg.figureSizeCm = [5,5];
cfg.axesSizeCm = [3,2.3];
result = runKDEMeanDotplot(cfg);
