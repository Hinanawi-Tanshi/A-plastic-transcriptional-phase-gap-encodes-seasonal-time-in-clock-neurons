function value = kdeFigureVisibility(cfg)
% KDEFIGUREVISIBILITY Allow the same pipeline to run interactively or headlessly.
value = 'on';
if ~cfg.showFigures, value = 'off'; end
end
