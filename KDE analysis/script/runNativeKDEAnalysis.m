function result = runNativeKDEAnalysis(cfg)
% RUNNATIVEKDEANALYSIS Fig2B: native-order 20-min binary-activity smoothing.
% Source: YJJ code/code/fig2b.m. Cached KDE columns are excluded from inputs.
% Only NaN time rows are removed. Signals and time order retain source rules.
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
traces = readKDETraces(cfg); result = struct();
qa = cell(numel(traces),1); tables = cell(numel(traces),1);
for i = 1:numel(traces)
    tr = traces(i); qa{i} = kdeTraceQA(tr,tr.time,tr.signal);
    valid = ~isnan(tr.time); t = tr.time(valid); signal = tr.signal(valid);
    probability = compute_KDE(t*60,signal,cfg.kdeWindowMinutes);
    result.traces(i).meta = tr.meta;
    result.traces(i).time_h = t; result.traces(i).signal = signal;
    result.traces(i).KDE = probability;
    tables{i} = table(repmat(tr.meta.TraceID,numel(t),1), ...
        repmat(tr.meta.Condition,numel(t),1),repmat(tr.meta.CellType,numel(t),1), ...
        t,signal,probability,'VariableNames', ...
        {'TraceID','Condition','CellType','Time_h','Signal','KDE'});
end
result.qa = vertcat(qa{:});
writetable(vertcat(tables{:}),fullfile(cfg.outputDir,cfg.outputStem + "_traces.csv"));
writetable(result.qa,fullfile(cfg.outputDir,cfg.outputStem + "_input_QA.csv"));
warnKDEInputQA(result.qa,cfg);
for type = ["M","E"]
    idx = find(arrayfun(@(x) x.meta.CellType == type,result.traces));
    fig = figure('Color','w','Units','centimeters','Position',[2 2 cfg.figureSizeCm], ...
        'Visible',kdeFigureVisibility(cfg),'Name',char(type + " cell"));
    for j = 1:numel(idx)
        tr = result.traces(idx(j));
        ax = subplot(numel(idx),1,j,'Parent',fig);
        plot(ax,tr.time_h,tr.KDE,'k-','LineWidth',1);
        ylim(ax,[-.1 1]); ylabel(ax,tr.meta.Condition,'Interpreter','none');
        set(ax,'FontName','Arial','FontSize',7,'TickDir','out','Box','off');
    end
    xlabel(ax,'Time (h)');
    exportKDEFigure(fig,fullfile(cfg.outputDir,cfg.outputStem + "_" + type + "cell"),cfg);
end
result = saveKDEResults(result,cfg);
end
