function result = runLPKDEAnalysis(cfg)
% RUNLPKDEANALYSIS Current Fig3CFI kernel regression and LP statistics.
% Source: Figure 3 LP/Fig 3C/fig3c.m (current modified version).
% Each column is one animal. Signal NaNs become zero for curves, while rows
% missing in either cell type are excluded from the kernel and the curve gap.
% Statistics use RAW data: complete pre3/post2 endpoints, exact permutation,
% Welch test and Holm correction across the six ZT x cell-type comparisons.
% NW uses local normalization and does not binarize fractional activities.
if ~isfolder(cfg.outputDir), mkdir(cfg.outputDir); end
qaRows = cell(0,5);
ZT_list = cfg.ZTList;
n_ZT = length(ZT_list);

window_width = cfg.kernelWindowMinutes;          % KDE窗口 30分钟
sigma_width = window_width / 6;
lp_pre_points = cfg.lpPrePoints;          % LP前使用3个原始时间点进行统计
lp_post_points = cfg.lpPostPoints;         % LP后使用最初2个原始时间点进行统计

color_E = cfg.EColor;
color_M = cfg.MColor;
color_control_shade = cfg.controlColor;  % 灰色 - Control

%% 2. 批量读取与 KDE 处理
all_results = repmat(struct('zt', []), 1, n_ZT);

for z = 1:n_ZT
    zt = ZT_list(z);
    fprintf('处理 ZT%d ...\n', zt);
    
    rawdata_E = readLPMatrix(cfg.inputFile, sprintf('EcellZT%d',zt));
    rawdata_E_ctrl = readLPMatrix(cfg.inputFile, sprintf('EcellZT%d_control',zt));
    rawdata_M = readLPMatrix(cfg.inputFile, sprintf('McellZT%d',zt));
    rawdata_M_ctrl = readLPMatrix(cfg.inputFile, sprintf('McellZT%d_control',zt));
    inputMatrices = {rawdata_E,rawdata_E_ctrl,rawdata_M,rawdata_M_ctrl};
    names = ["Ecell","Ecell_control","Mcell","Mcell_control"];
    for q = 1:4
        a = inputMatrices{q};
        assert(size(a,1)==size(rawdata_E,1),'kde:LPTimeAxis','ZT%d input row counts differ.',zt);
        delta = max(abs(a(:,1)-rawdata_E(:,1)));
        assert(all(isfinite(a(:,1))) && delta<=cfg.timeAxisToleranceHours, ...
            'kde:LPTimeAxis','ZT%d %s time axis differs from E experiment by %g hours.',zt,names(q),delta);
        qaRows(end+1,:) = {zt,names(q),size(a,1),size(a,2)-1,delta}; %#ok<AGROW>
    end

    % 时间轴 (假设所有文件时间一致，以Ecell为准)
    t_hours = rawdata_E(:, 1);
    t_min = t_hours * 60;
    
    % Cell 数据 (第2列至末尾)
    E_data = rawdata_E(:, 2:end);
    M_data = rawdata_M(:, 2:end);
    
    % Control 数据 (第2列至末尾)
    E_control = rawdata_E_ctrl(:, 2:end);
    M_control = rawdata_M_ctrl(:, 2:end);

    % 保留未平滑、未填补的原始数据，供LP前后统计使用
    E_data_raw = E_data;
    E_control_raw = E_control;
    M_data_raw = M_data;
    M_control_raw = M_control;
    
    num_E_cell = size(E_data, 2);
    num_E_control = size(E_control, 2);
    num_M_cell = size(M_data, 2);
    num_M_control = size(M_control, 2);
    
    fprintf('  Ecell: %d animals, %d controls | Mcell: %d animals, %d controls\n', ...
        num_E_cell, num_E_control, num_M_cell, num_M_control);
    
    % -------- 缺失值处理 --------
    % 判断 Cell 和 Control 均无数据的时间段，用于后续断开绘图
    missing_idx_E = all(isnan(E_data), 2) & all(isnan(E_control), 2);
    missing_idx_M = all(isnan(M_data), 2) & all(isnan(M_control), 2);
    missing_idx = missing_idx_E | missing_idx_M;
    lp_gap_idx = missing_idx_E & missing_idx_M;
    
    % NaN转0用于KDE计算
    E_data(isnan(E_data)) = 0;
    E_control(isnan(E_control)) = 0;
    M_data(isnan(M_data)) = 0;
    M_control(isnan(M_control)) = 0;
    
    % -------- Nadaraya-Watson 核回归 (1分钟上采样) --------
    t_upsample = (t_min(1) : cfg.upsampleMinutes : t_min(end))';
    
    
    valid_mask = ~missing_idx;
    t_min_valid = t_min(valid_mask);
    E_valid = E_data(valid_mask, :);
    E_control_valid = E_control(valid_mask, :);
    M_valid = M_data(valid_mask, :);
    M_control_valid = M_control(valid_mask, :);
    
    E_kde_upsample = compute_NW_activity(t_min_valid,E_valid,sigma_width,t_upsample);
    E_control_kde_upsample = compute_NW_activity(t_min_valid,E_control_valid,sigma_width,t_upsample);
    M_kde_upsample = compute_NW_activity(t_min_valid,M_valid,sigma_width,t_upsample);
    M_control_kde_upsample = compute_NW_activity(t_min_valid,M_control_valid,sigma_width,t_upsample);

    % -------- 断开缺失时间段 --------
    if ~isempty(missing_idx) && sum(missing_idx) > 0
        idx_missing = find(missing_idx);
        first_missing = idx_missing(1);
        last_missing = idx_missing(end);
        
        if first_missing > 1
            t_start_drop = t_min(first_missing - 1);
        else
            t_start_drop = t_min(first_missing);
        end
        
        if last_missing < length(t_min)
            t_end_drop = t_min(last_missing + 1);
        else
            t_end_drop = t_min(last_missing);
        end
        
        drop_idx = (t_upsample > t_start_drop) & (t_upsample < t_end_drop);
        E_kde_upsample(drop_idx, :) = NaN;
        E_control_kde_upsample(drop_idx, :) = NaN;
        M_kde_upsample(drop_idx, :) = NaN;
        M_control_kde_upsample(drop_idx, :) = NaN;
    end
    
    % -------- 计算均值和SEM --------
    E_mean = mean(E_kde_upsample, 2, 'omitnan');
    E_sem = std(E_kde_upsample, 0, 2, 'omitnan') / sqrt(num_E_cell);
    E_control_mean = mean(E_control_kde_upsample, 2, 'omitnan');
    E_control_sem = std(E_control_kde_upsample, 0, 2, 'omitnan') / sqrt(num_E_control);
    
    M_mean = mean(M_kde_upsample, 2, 'omitnan');
    M_sem = std(M_kde_upsample, 0, 2, 'omitnan') / sqrt(num_M_cell);
    M_control_mean = mean(M_control_kde_upsample, 2, 'omitnan');
    M_control_sem = std(M_control_kde_upsample, 0, 2, 'omitnan') / sqrt(num_M_control);
    
    % -------- 存储结果 --------
    all_results(z).zt = zt;
    all_results(z).t_plot = t_upsample / 60;
    all_results(z).E_mean = E_mean;
    all_results(z).E_sem = E_sem;
    all_results(z).E_control_mean = E_control_mean;
    all_results(z).E_control_sem = E_control_sem;
    all_results(z).M_mean = M_mean;
    all_results(z).M_sem = M_sem;
    all_results(z).M_control_mean = M_control_mean;
    all_results(z).M_control_sem = M_control_sem;
    all_results(z).t_E_hours = rawdata_E(:, 1);
    all_results(z).t_E_control_hours = rawdata_E_ctrl(:, 1);
    all_results(z).t_M_hours = rawdata_M(:, 1);
    all_results(z).t_M_control_hours = rawdata_M_ctrl(:, 1);
    all_results(z).E_raw = E_data_raw;
    all_results(z).E_control_raw = E_control_raw;
    all_results(z).M_raw = M_data_raw;
    all_results(z).M_control_raw = M_control_raw;
    all_results(z).lp_gap_idx = lp_gap_idx;
    all_results(z).E_kde = E_kde_upsample;
    all_results(z).E_control_kde = E_control_kde_upsample;
    all_results(z).M_kde = M_kde_upsample;
    all_results(z).M_control_kde = M_control_kde_upsample;
end

%% 3. 绘图: Ecell (ZT15, ZT18, ZT21 合并子图)
figure('Color', 'w', 'Visible', kdeFigureVisibility(cfg));
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, cfg.figureSizeCm]);

subplot_h = cfg.axesSizeCm(2);    % cm
subplot_w = cfg.axesSizeCm(1);    % cm
left_margin = 0.8;
bottom_margin = 0.4;
gap = 1;

for z = 1:n_ZT
    if isempty(all_results(z).zt)
        continue;
    end
    
    res = all_results(z);
    t_plot = res.t_plot;
    
    bottom_pos = bottom_margin + (n_ZT - z) * (subplot_h + gap);
    ax = axes('Units', 'centimeters', ...
              'Position', [left_margin, bottom_pos, subplot_w, subplot_h]);
    hold(ax, 'on');
    
    plot_trace_with_shade(ax, t_plot, res.E_mean, res.E_sem, color_E, color_E);
    plot_trace_with_shade(ax, t_plot, res.E_control_mean, res.E_control_sem, color_control_shade, color_control_shade);
    
    hold(ax, 'off');
    xlim(ax, [t_plot(1)-0.1, min(round(t_plot(1)+3), t_plot(end)+0.1)]);
    xticks(ax, ceil(t_plot(1)-0.1):1:floor(min(round(t_plot(1)+3), t_plot(end)+0.1)));
    ylim(ax, [0 1]);
    ylabel(ax, sprintf('ZT%d', res.zt), 'FontName', 'Arial', 'FontSize', 7);
    set(ax, 'LineWidth', 1, 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 7, 'Color', 'none');
end
% 用 annotation 代替 tiledlayout title
annotation('textbox', [0, 0.92, 1, 0.06], 'String', 'Ecell', ...
    'FontName', 'Arial', 'FontSize', 9, 'HorizontalAlignment', 'center', ...
    'EdgeColor', 'none');
exportKDEFigure(gcf,fullfile(cfg.outputDir,cfg.outputStem + "_Ecell_allZT"),cfg);

%% 4. 绘图: Mcell (ZT15, ZT18, ZT21 合并子图)
figure('Color', 'w', 'Visible', kdeFigureVisibility(cfg));
set(gcf, 'Units', 'centimeters', 'Position', [5, 5, cfg.figureSizeCm]);

for z = 1:n_ZT
    if isempty(all_results(z).zt)
        continue;
    end
    
    res = all_results(z);
    t_plot = res.t_plot;
    
    bottom_pos = bottom_margin + (n_ZT - z) * (subplot_h + gap);
    ax = axes('Units', 'centimeters', ...
              'Position', [left_margin, bottom_pos, subplot_w, subplot_h]);
    hold(ax, 'on');
    
    plot_trace_with_shade(ax, t_plot, res.M_mean, res.M_sem, color_M, color_M);
    plot_trace_with_shade(ax, t_plot, res.M_control_mean, res.M_control_sem, color_control_shade, color_control_shade);
    
    hold(ax, 'off');
    xlim(ax, [t_plot(1)-0.1, min(round(t_plot(1)+3), t_plot(end)+0.1)]);
    xticks(ax, ceil(t_plot(1)-0.1):1:floor(min(round(t_plot(1)+3), t_plot(end)+0.1)));
    ylim(ax, [0 1]);
    ylabel(ax, sprintf('ZT%d', res.zt), 'FontName', 'Arial', 'FontSize', 7);
    set(ax, 'LineWidth', 1, 'TickDir', 'out', 'FontName', 'Arial', 'FontSize', 7, 'Color', 'none');
end
annotation('textbox', [0, 0.92, 1, 0.06], 'String', 'Mcell', ...
    'FontName', 'Arial', 'FontSize', 9, 'HorizontalAlignment', 'center', ...
    'EdgeColor', 'none');
exportKDEFigure(gcf,fullfile(cfg.outputDir,cfg.outputStem + "_Mcell_allZT"),cfg);

summary_blocks = cell(n_ZT * 2, 1);
animal_blocks = cell(n_ZT * 4, 1);
raw_blocks = cell(n_ZT * 4, 1);
summary_count = 0;
animal_count = 0;

for z = 1:n_ZT
    if isempty(all_results(z).zt)
        continue;
    end

    res = all_results(z);
    [pre_idx, post_idx, gap_start, gap_end] = select_lp_window( ...
        res.lp_gap_idx, lp_pre_points, lp_post_points);

    if isempty(pre_idx)
        warning('ZT%d无法找到具有%d个pre和%d个post时间点的共同LP缺口，跳过统计。', ...
            res.zt, lp_pre_points, lp_post_points);
        continue;
    end

    cell_types = ["Ecell", "Mcell"];
    experiment_data = {res.E_raw, res.M_raw};
    control_data = {res.E_control_raw, res.M_control_raw};
    experiment_times = {res.t_E_hours, res.t_M_hours};
    control_times = {res.t_E_control_hours, res.t_M_control_hours};

    for c = 1:numel(cell_types)
        cell_type = cell_types(c);
        [experiment_animals, experiment_raw] = build_animal_tables( ...
            res.zt, cell_type, "Experiment", experiment_times{c}, ...
            experiment_data{c}, pre_idx, post_idx);
        [control_animals, control_raw] = build_animal_tables( ...
            res.zt, cell_type, "Control", control_times{c}, ...
            control_data{c}, pre_idx, post_idx);

        animal_count = animal_count + 1;
        animal_blocks{animal_count} = experiment_animals;
        raw_blocks{animal_count} = experiment_raw;
        animal_count = animal_count + 1;
        animal_blocks{animal_count} = control_animals;
        raw_blocks{animal_count} = control_raw;

        experiment_valid = experiment_animals.Included;
        control_valid = control_animals.Included;
        delta_experiment = experiment_animals.Delta(experiment_valid);
        delta_control = control_animals.Delta(control_valid);

        [permutation_p, num_permutations, observed_difference] = ...
            exact_permutation_test(delta_experiment, delta_control);
        [welch_t, welch_df, welch_p, ci_lower, ci_upper, hedges_g] = ...
            welch_test_with_effect(delta_experiment, delta_control, 0.05);

        pre_times_text = string(strjoin( ...
            compose('%.5f', res.t_E_hours(pre_idx)), ', '));
        post_times_text = string(strjoin( ...
            compose('%.5f', res.t_E_hours(post_idx)), ', '));

        summary_count = summary_count + 1;
        summary_blocks{summary_count} = table( ...
            res.zt, cell_type, res.t_E_hours(gap_start), res.t_E_hours(gap_end), ...
            pre_times_text, post_times_text, ...
            sum(experiment_valid), sum(control_valid), ...
            mean(experiment_animals.PreMean(experiment_valid), 'omitnan'), ...
            mean(experiment_animals.PostMean(experiment_valid), 'omitnan'), ...
            mean(delta_experiment, 'omitnan'), ...
            mean(control_animals.PreMean(control_valid), 'omitnan'), ...
            mean(control_animals.PostMean(control_valid), 'omitnan'), ...
            mean(delta_control, 'omitnan'), observed_difference, ...
            ci_lower, ci_upper, permutation_p, NaN, num_permutations, ...
            welch_t, welch_df, welch_p, NaN, hedges_g, false, ...
            'VariableNames', { ...
            'ZT', 'CellType', 'GapStart_h', 'GapEnd_h', ...
            'PreTimes_h', 'PostTimes_h', 'N_Experiment', 'N_Control', ...
            'ExperimentPreMean', 'ExperimentPostMean', 'ExperimentDeltaMean', ...
            'ControlPreMean', 'ControlPostMean', 'ControlDeltaMean', ...
            'DifferenceInDelta', 'CI95Lower', 'CI95Upper', ...
            'PermutationP', 'PermutationP_Holm', 'NumPermutations', ...
            'WelchT', 'WelchDF', 'WelchP', 'WelchP_Holm', ...
            'HedgesG', 'PrimarySignificant_Holm05'});
    end
end

if summary_count == 0
    warning('未生成LP前后统计结果。');
else
    summary_table = vertcat(summary_blocks{1:summary_count});
    animal_stats_table = vertcat(animal_blocks{1:animal_count});
    selected_raw_table = vertcat(raw_blocks{1:animal_count});

    summary_table.PermutationP_Holm = holm_adjust(summary_table.PermutationP);
    summary_table.WelchP_Holm = holm_adjust(summary_table.WelchP);
    summary_table.PrimarySignificant_Holm05 = ...
        summary_table.PermutationP_Holm < 0.05;

    writetable(summary_table,fullfile(cfg.outputDir,cfg.outputStem + "_LP_statistics.csv"));
    writetable(animal_stats_table,fullfile(cfg.outputDir,cfg.outputStem + "_LP_animals.csv"));
    writetable(selected_raw_table,fullfile(cfg.outputDir,cfg.outputStem + "_LP_selected_raw.csv"));
end

result = struct(); result.groups = all_results;
result.inputQA = cell2table(qaRows,'VariableNames', ...
    {'ZT','Input','NSamples','NAnimals','MaxTimeDifference_h'});
writetable(result.inputQA,fullfile(cfg.outputDir,cfg.outputStem + "_input_QA.csv"));
if summary_count > 0
    result.statistics = summary_table; result.animalStatistics = animal_stats_table;
    result.selectedRaw = selected_raw_table;
end
for z = 1:n_ZT
    res = all_results(z);
    for type = ["E","M"]
        n = size(res.(type + "_kde"),2); nc = size(res.(type + "_control_kde"),2);
        values = [res.(type + "_kde"),res.(type + "_control_kde")];
        names = [cellstr(compose('Experiment_A%d',1:n)),cellstr(compose('Control_A%d',1:nc))];
        curves = [table(res.t_plot,'VariableNames',{'Time_h'}),array2table(values,'VariableNames',names), ...
            table(res.(type + "_mean"),res.(type + "_sem"), ...
            res.(type + "_control_mean"),res.(type + "_control_sem"), ...
            'VariableNames',{'ExperimentMean','ExperimentSEM','ControlMean','ControlSEM'})];
        stem = fullfile(cfg.outputDir,cfg.panelNames(z) + "_" + type + "cell");
        writetable(curves,stem + "_curves.csv");
        if cfg.exportIndividualPanels
            fig = figure('Color','w','Units','centimeters','Position',[2 2 6 3], ...
                'Visible',kdeFigureVisibility(cfg));
            ax = axes(fig,'Units','centimeters','Position',[.8 .7 cfg.axesSizeCm]); hold(ax,'on');
            color = cfg.EColor; if type=="M", color = cfg.MColor; end
            plot_trace_with_shade(ax,res.t_plot,res.(type+"_mean"),res.(type+"_sem"),color,color);
            plot_trace_with_shade(ax,res.t_plot,res.(type+"_control_mean"),res.(type+"_control_sem"),cfg.controlColor,cfg.controlColor);
            last = min(round(res.t_plot(1)+3),res.t_plot(end)+.1);
            xlim(ax,[res.t_plot(1)-.1,last]); xticks(ax,ceil(res.t_plot(1)-.1):floor(last));
            ylim(ax,[0 1]); ylabel(ax,sprintf('ZT%d',res.zt));
            set(ax,'LineWidth',1,'TickDir','out','FontName','Arial','FontSize',7,'Color','none');
            exportKDEFigure(fig,stem,cfg);
        end
    end
end
result = saveKDEResults(result,cfg);
end

function raw = readLPMatrix(file,sheet)
raw = readmatrix(file,'Sheet',sheet,'NumHeaderLines',1);
assert(size(raw,2)>=2,'kde:LPInput','%s needs time and at least one animal column.',sheet);
end

function plot_trace_with_shade(ax, t_plot, y_mean, y_sem, color_line, color_shade)
    valid_idx = ~isnan(y_mean);
    diff_valid = diff([0; valid_idx; 0]);
    starts = find(diff_valid == 1);
    ends = find(diff_valid == -1) - 1;
    for k = 1:length(starts)
        idx = starts(k):ends(k);
        t_f = [t_plot(idx); flipud(t_plot(idx))];
        y_f = [y_mean(idx) + y_sem(idx); flipud(y_mean(idx) - y_sem(idx))];
        fill(ax, t_f, y_f, color_shade, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
    end
    plot(ax, t_plot, y_mean, 'Color', color_line, 'LineWidth', 1);
end

function [pre_idx, post_idx, gap_start, gap_end] = ...
        select_lp_window(gap_mask, n_pre_points, n_post_points)
    pre_idx = [];
    post_idx = [];
    gap_start = [];
    gap_end = [];

    gap_mask = logical(gap_mask(:));
    gap_edges = diff([false; gap_mask; false]);
    gap_starts = find(gap_edges == 1);
    gap_ends = find(gap_edges == -1) - 1;
    if isempty(gap_starts)
        return;
    end

    [~, longest_gap] = max(gap_ends - gap_starts + 1);
    gap_start = gap_starts(longest_gap);
    gap_end = gap_ends(longest_gap);

    row_index = (1:numel(gap_mask))';
    valid_before = find(~gap_mask & row_index < gap_start);
    valid_after = find(~gap_mask & row_index > gap_end);
    if numel(valid_before) < n_pre_points || numel(valid_after) < n_post_points
        pre_idx = [];
        post_idx = [];
        return;
    end

    pre_idx = valid_before(end - n_pre_points + 1:end);
    post_idx = valid_after(1:n_post_points);
end

function [animal_table, raw_table] = build_animal_tables( ...
        zt, cell_type, group_name, t_hours, data, pre_idx, post_idx)
    n_animals = size(data, 2);
    pre_values = data(pre_idx, :);
    post_values = data(post_idx, :);
    n_pre = sum(~isnan(pre_values), 1)';
    n_post = sum(~isnan(post_values), 1)';
    included = n_pre == numel(pre_idx) & n_post == numel(post_idx);

    pre_mean = mean(pre_values, 1, 'omitnan')';
    post_mean = mean(post_values, 1, 'omitnan')';
    delta = post_mean - pre_mean;
    delta(~included) = NaN;

    animal_index = (1:n_animals)';
    animal_id = compose('%s_%s_ZT%d_A%02d', ...
        cell_type, group_name, zt, animal_index);
    animal_table = table( ...
        repmat(zt, n_animals, 1), repmat(cell_type, n_animals, 1), ...
        repmat(group_name, n_animals, 1), animal_index, animal_id, ...
        pre_mean, post_mean, delta, n_pre, n_post, included, ...
        'VariableNames', {'ZT', 'CellType', 'Group', 'AnimalIndex', ...
        'AnimalID', 'PreMean', 'PostMean', 'Delta', 'N_Pre', 'N_Post', 'Included'});

    selected_times = [t_hours(pre_idx); t_hours(post_idx)];
    selected_values = [pre_values; post_values];
    phase = [repmat("Pre", numel(pre_idx), 1); ...
        repmat("Post", numel(post_idx), 1)];
    point_number = [(1:numel(pre_idx))'; (1:numel(post_idx))'];
    n_selected = numel(selected_times);

    raw_table = table( ...
        repmat(zt, n_animals * n_selected, 1), ...
        repmat(cell_type, n_animals * n_selected, 1), ...
        repmat(group_name, n_animals * n_selected, 1), ...
        repelem(animal_index, n_selected), repelem(animal_id, n_selected), ...
        repmat(phase, n_animals, 1), repmat(point_number, n_animals, 1), ...
        repmat(selected_times, n_animals, 1), selected_values(:), ...
        repelem(included, n_selected), ...
        'VariableNames', {'ZT', 'CellType', 'Group', 'AnimalIndex', ...
        'AnimalID', 'Phase', 'PointNumber', 'Time_h', 'Value', 'UsedInEndpoint'});
end

function [p_value, n_permutations, observed_difference] = ...
        exact_permutation_test(experiment_delta, control_delta)
    experiment_delta = experiment_delta(~isnan(experiment_delta));
    control_delta = control_delta(~isnan(control_delta));
    n_experiment = numel(experiment_delta);
    n_control = numel(control_delta);

    if n_experiment == 0 || n_control == 0
        p_value = NaN;
        n_permutations = 0;
        observed_difference = NaN;
        return;
    end

    combined = [experiment_delta(:); control_delta(:)];
    observed_difference = mean(experiment_delta) - mean(control_delta);
    combinations = nchoosek(1:numel(combined), n_experiment);
    n_permutations = size(combinations, 1);
    experiment_sums = sum(combined(combinations), 2);
    permutation_differences = experiment_sums / n_experiment - ...
        (sum(combined) - experiment_sums) / n_control;

    tolerance = 10 * eps(max(1, abs(observed_difference)));
    p_value = sum(abs(permutation_differences) >= ...
        abs(observed_difference) - tolerance) / n_permutations;
end

function [t_stat, df, p_value, ci_lower, ci_upper, hedges_g] = ...
        welch_test_with_effect(experiment_delta, control_delta, alpha)
    experiment_delta = experiment_delta(~isnan(experiment_delta));
    control_delta = control_delta(~isnan(control_delta));
    n_experiment = numel(experiment_delta);
    n_control = numel(control_delta);

    t_stat = NaN;
    df = NaN;
    p_value = NaN;
    ci_lower = NaN;
    ci_upper = NaN;
    hedges_g = NaN;
    if n_experiment < 2 || n_control < 2
        return;
    end

    mean_difference = mean(experiment_delta) - mean(control_delta);
    variance_experiment = var(experiment_delta, 0);
    variance_control = var(control_delta, 0);
    standard_error_squared = variance_experiment / n_experiment + ...
        variance_control / n_control;

    if standard_error_squared == 0
        if mean_difference == 0
            t_stat = 0;
            p_value = 1;
        else
            t_stat = sign(mean_difference) * Inf;
            p_value = 0;
        end
        df = Inf;
        ci_lower = mean_difference;
        ci_upper = mean_difference;
    else
        t_stat = mean_difference / sqrt(standard_error_squared);
        df = standard_error_squared^2 / ( ...
            (variance_experiment / n_experiment)^2 / (n_experiment - 1) + ...
            (variance_control / n_control)^2 / (n_control - 1));
        p_value = betainc(df / (df + t_stat^2), df / 2, 0.5);
        beta_quantile = betaincinv(alpha, df / 2, 0.5);
        t_critical = sqrt(df * (1 - beta_quantile) / beta_quantile);
        margin = t_critical * sqrt(standard_error_squared);
        ci_lower = mean_difference - margin;
        ci_upper = mean_difference + margin;
    end

    pooled_variance = ((n_experiment - 1) * variance_experiment + ...
        (n_control - 1) * variance_control) / ...
        (n_experiment + n_control - 2);
    if pooled_variance == 0
        if mean_difference == 0
            hedges_g = 0;
        else
            hedges_g = sign(mean_difference) * Inf;
        end
    else
        cohens_d = mean_difference / sqrt(pooled_variance);
        correction = 1 - 3 / (4 * (n_experiment + n_control) - 9);
        hedges_g = correction * cohens_d;
    end
end

function adjusted_p = holm_adjust(raw_p)
    adjusted_p = NaN(size(raw_p));
    valid_index = find(~isnan(raw_p));
    if isempty(valid_index)
        return;
    end

    [sorted_p, sort_order] = sort(raw_p(valid_index));
    n_tests = numel(sorted_p);
    multipliers = (n_tests:-1:1)';
    sorted_adjusted = min(1, cummax(sorted_p .* multipliers));
    adjusted_p(valid_index(sort_order)) = sorted_adjusted;
end
