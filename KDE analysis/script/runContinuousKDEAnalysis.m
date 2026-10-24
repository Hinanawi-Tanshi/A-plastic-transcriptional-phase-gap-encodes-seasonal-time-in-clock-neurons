function result = runContinuousKDEAnalysis(cfg)
% RUNCONTINUOUSKDEANALYSIS Fig2D: full native recordings and 0-24h summary.
% Source: YJJ code/code/fig2d.m, current data_source='both'.
% Legacy recordings share one time axis; newer recordings retain each cell's
% native axis in separate time/value column pairs within the recording tab.
% New M signals are nearest-neighbor aligned to E (or the first M) timestamps.
% Time is absolute hours. Full traces retain DD2; only the summary uses t<24h.
% Summary SEM is across cells. Empty bins are zero, as in the source.
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
traces = readKDETraces(cfg);
formats = cfg.dataSources;
keep = arrayfun(@(tr) ismember(tr.meta.SourceFormat,formats),traces);
traces = traces(keep);
assert(~isempty(traces),'kde:NoData','No traces for the selected data sources.');
flyIDs = unique(arrayfun(@(tr) tr.meta.FlyIndex,traces),'stable');
result = struct(); qa = cell(numel(traces),1);
for i = 1:numel(traces), qa{i} = kdeTraceQA(traces(i),traces(i).time,traces(i).signal); end
result.qa = vertcat(qa{:});
longTables = cell(numel(flyIDs),1);
for s = 1:numel(flyIDs)
    idx = find(arrayfun(@(tr) tr.meta.FlyIndex == flyIDs(s),traces));
    cells = traces(idx); metadata = vertcat(cells.meta);
    [~,order] = sort(metadata.CellIndex); cells = cells(order); metadata = metadata(order,:);
    nM = sum(metadata.CellType == "M"); hasE = any(metadata.CellType == "E");
    assert(all(metadata.CellType(1:nM)=="M") && ...
        (~hasE || metadata.CellType(end)=="E"), ...
        'kde:CellOrder','M columns must precede a single E column.');
    if all(metadata.SourceFormat == "old")
        t = cells(1).time; data = NaN(numel(t),numel(cells));
        for j = 1:numel(cells)
            assert(isequaln(t,cells(j).time),'kde:OldTimeAxis','Old-format cell axes must agree.');
            data(:,j) = cells(j).signal;
        end
    else
        reference = find(metadata.IsReference == 1);
        assert(isscalar(reference),'kde:ReferenceAxis','Exactly one reference trace is required per new-format fly.');
        ref = cells(reference); good = isfinite(ref.time); t = ref.time(good);
        data = NaN(numel(t),numel(cells));
        for j = 1:numel(cells)
            tr = cells(j); good = isfinite(tr.time);
            if hasE && j == reference
                data(:,j) = tr.signal(good);
            else
                data(:,j) = interp1(tr.time(good),tr.signal(good),t,'nearest',NaN);
            end
        end
    end
    % Keep the source's optional absolute-time cutoff semantics.
    if isfinite(cfg.maxTimeHours) && max(t)-min(t) > cfg.maxTimeHours
        selected = t <= cfg.maxTimeHours; t = t(selected); data = data(selected,:);
    end
    result.flies(s).rawTime_h = t;
    result.flies(s).alignedRawSignal = data;
    data(isnan(data)) = 0;
    [tMinutes,order] = sort(t*60); data = data(order,:);
    nativeKDE = zeros(size(data));
    for j = 1:size(data,2)
        nativeKDE(:,j) = compute_KDE(tMinutes,data(:,j),cfg.kdeWindowMinutes);
    end
    start = ceil(tMinutes(1)/cfg.upsampleMinutes)*cfg.upsampleMinutes;
    finish = floor(tMinutes(end)/cfg.upsampleMinutes)*cfg.upsampleMinutes;
    target = (start:cfg.upsampleMinutes:finish)';
    upsampled = zeros(numel(target),size(data,2));
    for j = 1:size(data,2)
        upsampled(:,j) = interp1(tMinutes,nativeKDE(:,j),target,'pchip');
    end
    result.flies(s).metadata = metadata;
    result.flies(s).flyIndex = flyIDs(s);
    result.flies(s).nMCells = nM; result.flies(s).hasE = hasE;
    result.flies(s).nativeTime_min = tMinutes;
    result.flies(s).nativeKDE = nativeKDE;
    result.flies(s).time_min = target;
    result.flies(s).KDE = upsampled;
    nT = numel(target); nC = size(data,2);
    assert(height(metadata)==nC,'kde:CellCount','Cell metadata and signal columns differ.');
    traceLabels = repmat(metadata.TraceID',nT,1);
    typeLabels = repmat(metadata.CellType',nT,1);
    longTables{s} = table(traceLabels(:), ...
        repmat(flyIDs(s),nT*nC,1),typeLabels(:), ...
        repmat(target/60,nC,1),upsampled(:), ...
        'VariableNames',{'TraceID','FlyIndex','CellType','Time_h','KDE'});
    if cfg.plotIndividualFlies
        fig = figure('Color','w','Units','centimeters','Position',[2 2 9 5], ...
            'Visible',kdeFigureVisibility(cfg),'Name',sprintf('Fly%d',flyIDs(s)));
        ax = axes(fig); hold(ax,'on');
        if nM>0, plot(ax,target/60,upsampled(:,1:nM),'Color',cfg.MColor); end
        if hasE, plot(ax,target/60,upsampled(:,end),'Color',cfg.EColor); end
        xlabel(ax,'Absolute time (h)'); ylabel(ax,'KDE');
        set(ax,'FontName','Arial','FontSize',7,'TickDir','out','Box','off');
        exportKDEFigure(fig,fullfile(cfg.outputDir,cfg.outputStem + "_Fly" + flyIDs(s)),cfg);
    end
end

% Retain original absolute <24h selection, cell pooling and SEM denominator.
bins = (0:cfg.upsampleMinutes:1440-cfg.upsampleMinutes)';
M = cell(numel(bins),1); E = cell(numel(bins),1);
for s = 1:numel(result.flies)
    fly = result.flies(s); selected = fly.time_min < 1440;
    times = mod(fly.time_min(selected),1440); values = fly.KDE(selected,:);
    for i = 1:numel(times)
        b = find(bins == round(times(i)));
        if isempty(b), continue; end
        if fly.nMCells>0, M{b} = [M{b};values(i,1:fly.nMCells)']; end %#ok<AGROW>
        if fly.hasE, E{b} = [E{b};values(i,end)]; end %#ok<AGROW>
    end
end
mnM = zeros(numel(bins),1); semM = mnM; mnE = mnM; semE = mnM;
for b = 1:numel(bins)
    mnM(b) = mean(M{b},'omitnan'); semM(b) = std(M{b},'omitnan')/sqrt(numel(M{b}));
    mnE(b) = mean(E{b},'omitnan'); semE(b) = std(E{b},'omitnan')/sqrt(numel(E{b}));
end
mnM(isnan(mnM))=0; semM(isnan(semM))=0;
mnE(isnan(mnE))=0; semE(isnan(semE))=0;
% Additional display smoothing acts on BOTH mean and SEM (source sigma=30min).
if cfg.summarySmoothingSigmaMinutes > 0
    sigma = cfg.summarySmoothingSigmaMinutes/cfg.upsampleMinutes;
    half = ceil(3*sigma); kernelTime = (-half:half)';
    kernel = exp(-kernelTime.^2/(2*sigma^2)); kernel = kernel/sum(kernel);
else
    kernel = 1;
end
summary = table(bins/60,mnM,semM,mnE,semE,cellfun(@numel,M),cellfun(@numel,E), ...
    conv(mnM,kernel,'same'),conv(semM,kernel,'same'), ...
    conv(mnE,kernel,'same'),conv(semE,kernel,'same'), ...
    'VariableNames',{'Time_h','M_Mean','M_SEM','E_Mean','E_SEM','N_MCells','N_ECells', ...
    'M_DisplayMean','M_DisplaySEM','E_DisplayMean','E_DisplaySEM'});
result.summary24h = summary;
result.summarySmoothingKernel = kernel;
writetable(vertcat(longTables{:}),fullfile(cfg.outputDir,cfg.outputStem + "_full_traces.csv"));
writetable(summary,fullfile(cfg.outputDir,cfg.outputStem + "_summary_0_24h.csv"));
writetable(result.qa,fullfile(cfg.outputDir,cfg.outputStem + "_input_QA.csv"));
warnKDEInputQA(result.qa,cfg);
fig = figure('Color','w','Units','centimeters','Position',[2 2 cfg.figureSizeCm], ...
    'Visible',kdeFigureVisibility(cfg));
ax = axes(fig,'Units','centimeters','Position',[.7 .7 cfg.axesSizeCm]); hold(ax,'on');
for type = ["M","E"]
    color = cfg.MColor; if type == "E", color = cfg.EColor; end
    mn = summary.(type + "_DisplayMean"); sem = summary.(type + "_DisplaySEM");
    fill(ax,[bins/60;flipud(bins/60)],[mn-sem;flipud(mn+sem)], ...
        color,'FaceAlpha',.2,'EdgeColor','none');
    plot(ax,bins/60,mn,'Color',color,'LineWidth',1);
end
set(ax,'LineWidth',.5,'TickDir','out','FontName','Arial','FontSize',7,'Color','none','Box','off');
ylim(ax,[-.1 1]); xlim(ax,[0 24]); xticks(ax,0:6:24); yticks(ax,0:.5:1);
exportKDEFigure(fig,fullfile(cfg.outputDir,cfg.outputStem + "_summary_0_24h"),cfg);
result = saveKDEResults(result,cfg);
end
