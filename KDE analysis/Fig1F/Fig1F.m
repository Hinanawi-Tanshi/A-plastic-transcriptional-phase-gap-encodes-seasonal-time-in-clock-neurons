% Fig1F: run the saved file with Run/F5 or run(fullPath).
% MATLAB R2023a. Inputs are read only; outputs go into results/.
% F9 also works when this saved file is the active MATLAB Editor document.
figureFile = string(mfilename('fullpath'));
if strlength(figureFile) == 0 && usejava('desktop')
    figureFile = string(matlab.desktop.editor.getActiveFilename);
end
[figureDir, figureName] = fileparts(figureFile);
if ~strcmp(figureName, 'Fig1F')
    error('kde:ScriptLocation', ...
        'Cannot locate Fig1F.m. Run the full saved file (Run/F5 or run(fullPath)).');
end
scriptDir = fullfile(fileparts(figureDir), 'script');
if ~isfile(fullfile(figureDir, 'data_binary.xlsx')) || ~isfolder(scriptDir)
    error('kde:FigureLayout', 'Expected data_binary.xlsx beside this file and a sibling script folder.');
end
addpath(scriptDir);

cfg = struct();
cfg.inputFile = fullfile(figureDir, 'data_binary.xlsx');
cfg.traceLayout = "aligned_pairs";
% Each recording keeps its native time/value pairs. Retain tab/column order.
cfg.outputDir = fullfile(figureDir, 'results');
cfg.outputStem = "Fig1F";
cfg.showFigures = true;
cfg.exportPNG = true;
% Each pair is one recorded trace; animal identities are unspecified.
cfg.conditions = ["ZT4","ZT8","ZT12","ZT24"];
cfg.kdeWindowMinutes = 1;
cfg.targetStepSeconds = 1;
cfg.alignedDurationMinutes = 10;
% Figure and axes dimensions affect display only.
cfg.figureSizeCm = [3,3];
cfg.axesSizeCm = [0.9,1.6];
result = runAlignedKDEAnalysis(cfg);
