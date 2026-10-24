function result = runAlignedKDEAnalysis(cfg)
% RUNALIGNEDKDEANALYSIS Fig1F: 1-min KDE at 1-s resolution, midpoint alignment.
% Source: YJJ code/code/batch.m. Every time/signal pair has its own time axis.
% No animal identities are inferred. Both time and signal NaNs are excluded.
% compute_KDE target-time mode is essential; do not replace with convolution.
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
traces = readKDETraces(cfg);
tRelative = (0:cfg.alignedDurationMinutes*60)'/60;
result = struct(); result.traces = traces;
qa = cell(numel(traces),1); fullTables = cell(numel(traces),1);
for g = 1:numel(cfg.conditions)
    idx = find(arrayfun(@(x) x.meta.Condition == cfg.conditions(g), traces));
    assert(~isempty(idx), 'kde:MissingCondition', 'No traces for %s.', cfg.conditions(g));
    aligned = NaN(numel(tRelative), numel(idx));
    for j = 1:numel(idx)
        k = idx(j); tr = traces(k);
        qa{k} = kdeTraceQA(tr, tr.time, tr.signal);
        valid = ~isnan(tr.time) & ~isnan(tr.signal);
        t = tr.time(valid); signal = tr.signal(valid);
        assert(numel(t) >= 2, 'kde:ShortTrace', '%s needs at least two samples.', tr.meta.TraceID);
        tt = (min(t):cfg.targetStepSeconds/3600:max(t))';
        fullKDE = compute_KDE(t, signal, cfg.kdeWindowMinutes/60, tt);
        midpoint = (min(t)+max(t))/2;
        target = midpoint + (tRelative-cfg.alignedDurationMinutes/2)/60;
        aligned(:,j) = compute_KDE(t, signal, cfg.kdeWindowMinutes/60, target);
        result.traces(k).fullTime_h = tt;
        result.traces(k).fullKDE = fullKDE;
        result.traces(k).alignedKDE = aligned(:,j);
        fullTables{k} = table(repmat(tr.meta.TraceID,numel(tt),1), ...
            repmat(tr.meta.Condition,numel(tt),1),tt,fullKDE, ...
            'VariableNames',{'TraceID','Condition','Time_h','KDE'});
    end
    mn = mean(aligned,2,'omitnan'); sd = std(aligned,0,2,'omitnan');
    n = sum(~isnan(aligned),2); sem = sd./sqrt(n);
    result.groups(g).condition = cfg.conditions(g);
    result.groups(g).traceIndices = idx;
    result.groups(g).relativeTime_min = tRelative;
    result.groups(g).aligned = aligned;
    result.groups(g).mean = mn; result.groups(g).std = sd;
    result.groups(g).sem = sem; result.groups(g).n = n;
    names = arrayfun(@(x) char(x.meta.TraceID), traces(idx), 'UniformOutput', false);
    out = [table(tRelative,'VariableNames',{'RelativeTime_min'}), ...
        array2table(aligned,'VariableNames',names), ...
        table(mn,sd,sem,n,'VariableNames',{'Mean','Std','SEM','NTraces'})];
    stem = fullfile(cfg.outputDir,cfg.outputStem + "_" + cfg.conditions(g));
    writetable(out,stem + "_aligned.csv");
    fig = figure('Color','w','Units','centimeters','Position',[2 2 cfg.figureSizeCm], ...
        'Visible',kdeFigureVisibility(cfg),'Name',char(cfg.conditions(g)));
    ax = axes(fig,'Units','centimeters','Position',[0.7 0.7 cfg.axesSizeCm]);
    hold(ax,'on');
    upper = min(1,mn+sem); lower = max(0,mn-sem);
    fill(ax,[tRelative;flipud(tRelative)],[upper;flipud(lower)], ...
        [.3 .3 .3],'EdgeColor','none','FaceAlpha',.4);
    plot(ax,tRelative,aligned,'Color',[.8 .8 .8],'LineWidth',.5);
    plot(ax,tRelative,mn,'k-','LineWidth',1);
    ylim(ax,[-.1 1]); xlim(ax,[0 cfg.alignedDurationMinutes]);
    set(ax,'LineWidth',1,'TickDir','out','FontName','Arial','FontSize',7,'Color','none');
    exportKDEFigure(fig,stem,cfg);
end
result.qa = vertcat(qa{:});
writetable(result.qa,fullfile(cfg.outputDir,cfg.outputStem + "_input_QA.csv"));
warnKDEInputQA(result.qa,cfg);
writetable(vertcat(fullTables{:}),fullfile(cfg.outputDir,cfg.outputStem + "_full_traces.csv"));
result = saveKDEResults(result,cfg);
end
