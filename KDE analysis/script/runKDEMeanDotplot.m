function result = runKDEMeanDotplot(cfg)
% RUNKDEMEANDOTPLOT Fig1H, Method7 and Method7b of Task3_Suppl_Analysis.m.
% Recompute the producer's 20-min KDE from raw input, then preserve Method7.
% The pooled-cell plot is descriptive; the Kruskal-Wallis test uses animals.
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
oldRNG = rng; rngCleanup = onCleanup(@() rng(oldRNG)); %#ok<NASGU>
traces = readKDETraces(cfg);
conditions = cellstr(cfg.conditions); n_conditions = numel(conditions);
cond_varNames = strrep(conditions,'LD4 ','');
flyData = cell(n_conditions,1); qa = cell(numel(traces),1);
cellMeta = cell(0,6); traceTables = cell(0,1);
for c = 1:n_conditions
    idx = find(arrayfun(@(tr) tr.meta.Condition == cfg.conditions(c),traces));
    ids = unique(arrayfun(@(tr) tr.meta.FlyIndex,traces(idx)),'stable');
    assert(isequal(ids(:),(1:numel(ids))'),'kde:FlyOrder','Fly indices must be consecutive within condition.');
    for f_idx = 1:numel(ids)
        ii = idx(arrayfun(@(tr) tr.meta.FlyIndex == ids(f_idx),traces(idx)));
        [~,order] = sort(arrayfun(@(tr) tr.meta.CellIndex,traces(ii))); ii = ii(order);
        fly = struct(); fly.cells = cell(numel(ii),1); fly.n_cells = numel(ii);
        fly.fileName = char(traces(ii(1)).meta.SourceSheet);
        fly.flyName = sprintf('Fly%d_%s',f_idx,fly.fileName); fly.condition = conditions{c};
        for s_idx = 1:numel(ii)
            k = ii(s_idx); tr = traces(k); qa{k} = kdeTraceQA(tr,tr.time,tr.signal);
            valid = ~isnan(tr.time); t = tr.time(valid); sig = tr.signal(valid);
            kde = compute_KDE(t*60,sig,cfg.kdeWindowMinutes);
            fly.cells{s_idx} = struct('time',t,'signal',sig,'name',char(tr.meta.CellName), ...
                'KDE',kde,'time_min',t*60);
            cellMeta(end+1,:) = {tr.meta.TraceID,c,f_idx,s_idx,tr.meta.CellName,tr.meta.SourceFile}; %#ok<AGROW>
            traceTables{end+1,1} = table(repmat(tr.meta.TraceID,numel(t),1),t,kde, ...
                'VariableNames',{'TraceID','Time_h','KDE'}); %#ok<AGROW>
        end
        fly.time_hours = fly.cells{1}.time;
        flyData{c}{f_idx} = fly;
    end
end
fly_plot_x = cell(n_conditions,1);
cell_plot_x = NaN(size(cellMeta,1),1);
%% ========== 方法7: 果蝇KDE均值 (带jitter dotplot) ==========
% 先求每个细胞KDE的时间均值，再跨细胞平均，得到每只果蝇一个值。

fprintf('\n========== 方法7: 果蝇KDE均值 dotplot ==========\n');

% 收集各条件下每只果蝇的KDE均值
fly_kde_mean = cell(n_conditions, 1);
fly_kde_sem  = cell(n_conditions, 1);  % 果蝇内细胞间SEM

for c = 1:n_conditions
    n_flies_c = length(flyData{c});
    fly_kde_mean{c} = zeros(n_flies_c, 1);
    fly_kde_sem{c}  = zeros(n_flies_c, 1);
    
    for f_idx = 1:n_flies_c
        fly = flyData{c}{f_idx};
        n_cells = fly.n_cells;
        
        % 每细胞KDE的时间均值
        cell_kde_means = zeros(n_cells, 1);
        for s_idx = 1:n_cells
            cell_kde_means(s_idx) = nanmean(fly.cells{s_idx}.KDE);
        end
        
        fly_kde_mean{c}(f_idx) = mean(cell_kde_means);  % 果蝇均值
        fly_kde_sem{c}(f_idx)  = std(cell_kde_means) / sqrt(n_cells);  % 果蝇内SEM
        
        fprintf('  [%s] Fly%d: KDE均值=%.5f ± %.5f (SEM, n_cells=%d)\n', ...
            conditions{c}, f_idx, fly_kde_mean{c}(f_idx), fly_kde_sem{c}(f_idx), n_cells);
    end
    
    fprintf('  --- 条件汇总: mean=%.5f, std=%.5f (n_flies=%d) ---\n', ...
        mean(fly_kde_mean{c}), std(fly_kde_mean{c}), n_flies_c);
end

% ----- 带jitter的dotplot -----
flyFig = figure('Name', 'Fly KDE Mean Dotplot', 'Position', [100, 100, 450, 400], 'Visible', kdeFigureVisibility(cfg));
hold on;

colors_cond = lines(n_conditions);
jitter_width = cfg.flyJitterWidth;
rng(cfg.rngSeed);  % 固定jitter随机种子，保证可重复

for c = 1:n_conditions
    vals = fly_kde_mean{c};
    n = length(vals);
    
    % 生成jittered x坐标
    x_jitter = c + jitter_width * (randn(n, 1));
    fly_plot_x{c} = x_jitter;
    
    % 散点 (半透明)
    scatter(x_jitter, vals, 70, colors_cond(c,:), 'filled', ...
        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
    
    % 均值 ± SEM 条形
    mn = mean(vals);
    sem = std(vals) / sqrt(n);
    
    % 均值横线
    plot([c - 0.25, c + 0.25], [mn, mn], '-', 'Color', colors_cond(c,:), ...
        'LineWidth', 2.5);
    % SEM 误差线
    plot([c, c], [mn - sem, mn + sem], '-', 'Color', colors_cond(c,:), ...
        'LineWidth', 1.5);
    
    fprintf('  [%s] n=%d, mean=%.5f, SEM=%.5f\n', ...
        conditions{c}, n, mn, sem);
end

% Kruskal-Wallis 非参数检验（条件间差异）
if n_conditions >= 2
    all_vals_kw = [];
    all_groups_kw = [];
    for c = 1:n_conditions
        all_vals_kw = [all_vals_kw; fly_kde_mean{c}];
        all_groups_kw = [all_groups_kw; c * ones(length(fly_kde_mean{c}), 1)];
    end
    [p_kw, tbl_kw] = kruskalwallis(all_vals_kw, all_groups_kw, 'off');
    fprintf('  Kruskal-Wallis: chi2=%.3f, p=%.4f\n', tbl_kw{2,5}, p_kw);
    
    % 在图上标注p值
    y_max = max(cellfun(@max, fly_kde_mean));
    y_min = min(cellfun(@min, fly_kde_mean));
    y_range = y_max - y_min;
    if p_kw < 0.001
        p_str = 'p < 0.001';
    elseif p_kw < 0.01
        p_str = sprintf('p = %.3f', p_kw);
    else
        p_str = sprintf('p = %.3f', p_kw);
    end
    text(mean(1:n_conditions), y_max + 0.08 * y_range, p_str, ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    % Reserve headroom so the p-value stays inside the axes below the title.
    ylim([min(0,y_min), y_max + 0.2 * max(y_range,eps)]);
end

set(gca, 'XTick', 1:n_conditions, 'XTickLabel', cond_varNames);
xtickangle(20);
ylabel('Mean KDE Probability');
title('Per-fly mean KDE');
grid on;
box on;

exportKDEFigure(flyFig, fullfile(cfg.outputDir,cfg.outputStem + "_fly_dotplot"), cfg);

%% 方法7b: 散点图 (每个细胞一个点, 果蝇按均值排序, 每果蝇一列)
fig = figure('Name', 'Cell KDE Mean Dotplot', ...
    'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, cfg.figureSizeCm], 'Visible', kdeFigureVisibility(cfg));
hold on;

% 收集每个细胞的KDE均值，以及所属果蝇信息
cell_vals_all = [];     % KDE均值
cell_group_all = [];    % 条件编号
cell_fly_all = [];      % 果蝇编号 (同条件内)

for c = 1:n_conditions
    fly_counter = 0;
    for f_idx = 1:length(flyData{c})
        fly = flyData{c}{f_idx};
        fly_counter = fly_counter + 1;
        for s_idx = 1:fly.n_cells
            cell_val = nanmean(fly.cells{s_idx}.KDE);
            cell_vals_all(end+1) = cell_val;
            cell_group_all(end+1) = c;
            cell_fly_all(end+1) = fly_counter;
        end
    end
end

cell_vals_all = cell_vals_all(:);
cell_group_all = cell_group_all(:);
cell_fly_all = cell_fly_all(:);

% 果蝇颜色调色板 (跨条件可重复)
max_flies_per_cond = max(cellfun(@(x) length(x), flyData));
fly_colors = lines(max_flies_per_cond);
rng(cfg.rngSeed);

% 计算每个条件内每只果蝇的均值，并按均值排序
cond_fly_means = cell(n_conditions, 1);   % {c}(fly_idx) = mean KDE
cond_fly_order = cell(n_conditions, 1);   % 排序后的果蝇原始编号
cond_fly_pos   = cell(n_conditions, 1);   % 排序后每果蝇的x位置

x_offset = 0;
cond_centers = zeros(n_conditions, 1);   % 每组中心x位置
cond_x_ranges = zeros(n_conditions, 2);  % 每组x范围 [start, end]

for c = 1:n_conditions
    idx_cond = cell_group_all == c;
    flies_c = unique(cell_fly_all(idx_cond));
    n_flies_c = length(flies_c);
    
    % 计算每只果蝇的均值
    fly_means_c = zeros(n_flies_c, 1);
    for fi = 1:n_flies_c
        idx_f = idx_cond & (cell_fly_all == flies_c(fi));
        fly_means_c(fi) = mean(cell_vals_all(idx_f));
    end
    
    % 按均值升序排列
    [~, sort_idx] = sort(fly_means_c);
    cond_fly_order{c} = flies_c(sort_idx);          % 排序后的果蝇编号
    cond_fly_means{c} = fly_means_c(sort_idx);      % 排序后的均值
    
    % 分配x位置
    x_start = x_offset + 1;
    x_positions = x_start : (x_start + n_flies_c - 1);
    cond_fly_pos{c} = x_positions;
    
    cond_x_ranges(c, :) = [x_start - 0.5, x_start + n_flies_c - 1 + 1.5];
    cond_centers(c) = mean([x_start, x_start + n_flies_c - 1]);
    
    x_offset = x_start + n_flies_c + 1;  % 下一组起始 (组间留1间距 + 1列mean±SEM)
    
    % 绘制每只果蝇的细胞散点
    cell_jitter = cfg.cellJitterWidth;  % 果蝇列内细胞水平jitter
    for fi = 1:n_flies_c
        fly_id = cond_fly_order{c}(fi);
        idx_fly = idx_cond & (cell_fly_all == fly_id);
        cell_vals_fly = cell_vals_all(idx_fly);
        n_cells_fly = sum(idx_fly);
        
        x_fly = cond_fly_pos{c}(fi);
        x_jitter = x_fly + cell_jitter * randn(n_cells_fly, 1);
        cell_plot_x(idx_fly) = x_jitter;
        
        scatter(x_jitter, cell_vals_fly, 8, fly_colors(fi, :), 'filled', ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'k', 'LineWidth', 0.2);
    end
    
    % 绘制该条件的均值 ± SEM (放在每组最后一列右侧)
    vals_c = cell_vals_all(idx_cond);
    mn_c = mean(vals_c);
    sem_c = std(vals_c) / sqrt( length (vals_c) );
    
    x_right = cond_fly_pos{c}(end) + 1.0;  % 最后一列右侧+1, 独占一列
    bar_hw = 0.25;  % 均值横线半宽
    
    % 均值横线 (0.5pt)
    plot([x_right - bar_hw, x_right + bar_hw], [mn_c, mn_c], 'k-', 'LineWidth', 0.5);
    % SEM竖线 (0.5pt)
    plot([x_right, x_right], [mn_c - sem_c, mn_c + sem_c], 'k-', 'LineWidth', 0.5);
end

% X轴标签: 条件名放在每组中央
set(gca, 'XTick', cond_centers, 'XTickLabel', cond_varNames);
xtickangle(20);
xlim([-0.5, x_offset]);

% 组间分隔线
for c = 1:n_conditions-1
    sep_x = (cond_x_ranges(c, 2) + cond_x_ranges(c+1, 1)) / 2;
    xline(sep_x, ':', 'Color', [0.2, 0.2, 0.2], 'LineWidth', 0.8);
end

ax = gca;
ax.Units = 'centimeters';
ax.Position(3:4) = cfg.axesSizeCm;
ax.LineWidth = 0.8;
ax.TickDir = 'out';
ax.TickLength = [0.02, 0.02];
ax.XTickLabelRotation = 0;
ax.XRuler.TickLabelGapOffset = -1;

ax.FontName = 'Arial';
ax.FontSize = 7;
% ax.FontWeight = "bold";
ax.Color = 'none';
box off;
hold off;

%% 保存
set(findall(fig, '-property', 'FontUnits'), ...
    'FontUnits', 'points');

set(findall(fig, '-property', 'FontSize'), ...
    'FontSize', 7);

set(findall(fig, '-property', 'FontSizeMode'), ...
    'FontSizeMode', 'manual');

set(findall(fig, '-property', 'FontName'), ...
    'FontName', 'Arial');

% 让导出纸张尺寸与图窗物理尺寸完全一致
fig.PaperUnits = 'centimeters';
fig.PaperSize = cfg.figureSizeCm;
fig.PaperPosition = [0, 0, cfg.figureSizeCm];
fig.PaperPositionMode = 'manual';
fig.InvertHardcopy = 'off';

drawnow;

exportKDEFigure(fig, fullfile(cfg.outputDir,cfg.outputStem + "_cell_dotplot"), cfg);

cellTable = cell2table(cellMeta,'VariableNames', ...
    {'TraceID','ConditionIndex','FlyIndex','CellIndex','CellName','SourceFile'});
conditionNames = cfg.conditions(:);
cellTable.Condition = conditionNames(cellTable.ConditionIndex);
cellTable.MeanKDE = cell_vals_all; cellTable.PlotX = cell_plot_x;
flyRows = cell(0,7); summaryRows = cell(n_conditions,7);
for c = 1:n_conditions
    vals = fly_kde_mean{c}; cellVals = cell_vals_all(cell_group_all==c);
    for f_idx = 1:numel(vals)
        flyRows(end+1,:) = {conditions{c},f_idx,flyData{c}{f_idx}.n_cells, ...
            vals(f_idx),fly_kde_sem{c}(f_idx),fly_plot_x{c}(f_idx), ...
            find(cond_fly_order{c}==f_idx)}; %#ok<AGROW>
    end
    summaryRows(c,:) = {conditions{c},numel(vals),numel(cellVals), ...
        mean(vals),std(vals)/sqrt(numel(vals)),mean(cellVals),std(cellVals)/sqrt(numel(cellVals))};
end
flyTable = cell2table(flyRows,'VariableNames', ...
    {'Condition','FlyIndex','NCells','MeanKDE','WithinFlyCellSEM','PlotX','CellPlotOrder'});
summaryTable = cell2table(summaryRows,'VariableNames', ...
    {'Condition','NFlies','NCells','FlyMean','FlySEM','PooledCellMean','PooledCellSEM'});
result = struct(); result.flyData = flyData; result.cells = cellTable;
result.flies = flyTable; result.conditions = summaryTable;
result.flyKDEMean = fly_kde_mean; result.flyKDESEM = fly_kde_sem;
result.conditionFlyOrder = cond_fly_order; result.conditionFlyPosition = cond_fly_pos;
result.conditionCenters = cond_centers; result.conditionXRanges = cond_x_ranges;
result.flyPlotX = fly_plot_x; result.cellPlotX = cell_plot_x;
result.qa = vertcat(qa{:});
if n_conditions >= 2
    result.kruskalWallisP = p_kw; result.kruskalWallisTable = tbl_kw;
    writecell(tbl_kw,fullfile(cfg.outputDir,cfg.outputStem + "_KruskalWallis_animals.csv"));
end
writetable(cellTable,fullfile(cfg.outputDir,cfg.outputStem + "_cell_means.csv"));
writetable(flyTable,fullfile(cfg.outputDir,cfg.outputStem + "_fly_means.csv"));
writetable(summaryTable,fullfile(cfg.outputDir,cfg.outputStem + "_condition_summary.csv"));
writetable(vertcat(traceTables{:}),fullfile(cfg.outputDir,cfg.outputStem + "_KDE_traces.csv"));
writetable(result.qa,fullfile(cfg.outputDir,cfg.outputStem + "_input_QA.csv"));
warnKDEInputQA(result.qa,cfg);
result = saveKDEResults(result,cfg);
end
