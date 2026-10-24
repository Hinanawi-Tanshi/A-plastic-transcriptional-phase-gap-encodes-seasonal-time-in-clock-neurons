% Fig3CFI: run the saved file with Run/F5 or run(fullPath).
% MATLAB R2023a. Inputs are read only; outputs go into results/.
% F9 also works when this saved file is the active MATLAB Editor document.
figureFile = string(mfilename('fullpath'));
if strlength(figureFile) == 0 && usejava('desktop')
    figureFile = string(matlab.desktop.editor.getActiveFilename);
end
[figureDir, figureName] = fileparts(figureFile);
if ~strcmp(figureName, 'Fig3CFI')
    error('kde:ScriptLocation', ...
        'Cannot locate Fig3CFI.m. Run the full saved file (Run/F5 or run(fullPath)).');
end
scriptDir = fullfile(fileparts(figureDir), 'script');
if ~isfile(fullfile(figureDir, 'data.xlsx')) || ~isfolder(scriptDir)
    error('kde:FigureLayout', 'Expected data.xlsx beside this file and a sibling script folder.');
end
addpath(scriptDir);

cfg = struct();
cfg.inputFile = fullfile(figureDir, 'data.xlsx');
cfg.outputDir = fullfile(figureDir, 'results');
cfg.outputStem = "Fig3CFI";
cfg.showFigures = true;
cfg.exportPNG = true;
% Current Figure 3 LP/Fig 3C/fig3c.m analysis (ZT15/18/21).
% Each signal column is one animal. Fractional activity is preserved.
% Nadaraya-Watson regression uses local weight normalization, without truncation.
cfg.ZTList = [15,18,21];
cfg.panelNames = ["Fig3C","Fig3F","Fig3I"];
cfg.kernelWindowMinutes = 30; % sigma = window/6 = 5 min
cfg.upsampleMinutes = 1;
cfg.lpPrePoints = 3;
cfg.lpPostPoints = 2;
% Curve axis follows E experiment, as in the source. Other axes must agree
% within recorded rounding precision; raw per-file axes remain in statistics.
cfg.timeAxisToleranceHours = 1e-5;
cfg.EColor = [192/255,0,0];
cfg.MColor = [192/255,0,0];
cfg.controlColor = [0.6,0.6,0.6];
cfg.exportIndividualPanels = true;
cfg.figureSizeCm = [6,10];
cfg.axesSizeCm = [4.5,1.8];
result = runLPKDEAnalysis(cfg);
