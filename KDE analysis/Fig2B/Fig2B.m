% Fig2B: run the saved file with Run/F5 or run(fullPath).
% MATLAB R2023a. Inputs are read only; outputs go into results/.
% F9 also works when this saved file is the active MATLAB Editor document.
figureFile = string(mfilename('fullpath'));
if strlength(figureFile) == 0 && usejava('desktop')
    figureFile = string(matlab.desktop.editor.getActiveFilename);
end
[figureDir, figureName] = fileparts(figureFile);
if ~strcmp(figureName, 'Fig2B')
    error('kde:ScriptLocation', ...
        'Cannot locate Fig2B.m. Run the full saved file (Run/F5 or run(fullPath)).');
end
scriptDir = fullfile(fileparts(figureDir), 'script');
if ~isfile(fullfile(figureDir, 'data_binary.xlsx')) || ~isfolder(scriptDir)
    error('kde:FigureLayout', 'Expected data_binary.xlsx beside this file and a sibling script folder.');
end
addpath(scriptDir);

cfg = struct();
cfg.inputFile = fullfile(figureDir, 'data_binary.xlsx');
cfg.traceLayout = "native_single";
% Each recording keeps its native time/value pairs. Retain tab/column order.
cfg.outputDir = fullfile(figureDir, 'results');
cfg.outputStem = "Fig2B";
cfg.showFigures = true;
cfg.exportPNG = true;
% Native convolution, with the source timestamp order preserved.
% Input QA is exported; two source traces contain repeated/backwards times.
cfg.kdeWindowMinutes = 20;
cfg.figureSizeCm = [5,12];
result = runNativeKDEAnalysis(cfg);
