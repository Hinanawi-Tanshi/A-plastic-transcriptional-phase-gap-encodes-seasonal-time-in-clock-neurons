function cleanup = bootstrapFigureContext(cfg)
% Create an output directory and restore figure visibility after the call.
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
previousVisibility = get(groot, 'DefaultFigureVisible');
if cfg.showFigures
    set(groot, 'DefaultFigureVisible', 'on');
else
    set(groot, 'DefaultFigureVisible', 'off');
end
cleanup = onCleanup(@() set(groot, 'DefaultFigureVisible', previousVisibility));
end
