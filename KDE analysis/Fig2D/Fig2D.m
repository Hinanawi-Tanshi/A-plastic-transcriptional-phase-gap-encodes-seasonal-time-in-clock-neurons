% Fig2D: run the saved file with Run/F5 or run(fullPath).
% MATLAB R2023a. Inputs are read only; outputs go into results/.
% F9 also works when this saved file is the active MATLAB Editor document.
figureFile = string(mfilename('fullpath'));
if strlength(figureFile) == 0 && usejava('desktop')
    figureFile = string(matlab.desktop.editor.getActiveFilename);
end
[figureDir, figureName] = fileparts(figureFile);
if ~strcmp(figureName, 'Fig2D')
    error('kde:ScriptLocation', ...
        'Cannot locate Fig2D.m. Run the full saved file (Run/F5 or run(fullPath)).');
end
scriptDir = fullfile(fileparts(figureDir), 'script');
if ~isfile(fullfile(figureDir, 'data_binary.xlsx')) || ~isfolder(scriptDir)
    error('kde:FigureLayout', 'Expected data_binary.xlsx beside this file and a sibling script folder.');
end
addpath(scriptDir);

cfg = struct();
cfg.inputFile = fullfile(figureDir, 'data_binary.xlsx');
cfg.traceLayout = "continuous_cells";
% Each recording keeps its native time/value pairs. Retain tab/column order.
cfg.legacySheetCount = 8; % first 8 anonymous flies retain the original shared-time processing
cfg.outputDir = fullfile(figureDir, 'results');
cfg.outputStem = "Fig2D";
cfg.showFigures = true;
cfg.exportPNG = true;
% Both input formats; native absolute timestamps include recordings >24h.
cfg.dataSources = ["old","new"];
cfg.maxTimeHours = inf; % source's optional absolute-time cutoff; inf keeps full data
cfg.kdeWindowMinutes = 180;
cfg.upsampleMinutes = 3; % PCHIP after native KDE
% The source's summary selects absolute t<24h, pools cells, and zero-fills empty bins.
% Both mean and SEM receive this additional smoothing for display.
cfg.summarySmoothingSigmaMinutes = 30;
cfg.plotIndividualFlies = false; % full traces are exported regardless
cfg.MColor = [1,128/255,0];
cfg.EColor = [0,0,1];
cfg.figureSizeCm = [4,3];
cfg.axesSizeCm = [2.4,1.6];
result = runContinuousKDEAnalysis(cfg);
