function publicResult = runPeakShiftComparison(cfg)
% Shared publication pipeline; settings are supplied by the figure entry script.
figureCleanup = bootstrapFigureContext(cfg); %#ok<NASGU>
scriptFolder = cfg.outputDir;
dataFolder = cfg.dataFolder;
smokeTestMode = false;
controlPath = cfg.controlPath;
experimentPaths = cfg.experimentPaths;
experimentTags = cfg.experimentTags;
experimentLabels = cfg.experimentLabels;
controlTag = cfg.controlTag;
controlLabel = cfg.controlLabel;
M_color = cfg.M_color;
E_color = cfg.E_color;
controlColor = cfg.controlColor;
M_sheet = cfg.M_sheet;
E_sheet = cfg.E_sheet;
nInner = cfg.nInner;
nOuter = cfg.nOuter;
delta0 = cfg.delta0;
alpha = cfg.alpha;
rngSeed = cfg.rngSeed;
mainFigureBinWidth = cfg.mainFigureBinWidth;
multiplicityAdjustment = cfg.multiplicityAdjustment;
outputExcel = fullfile(scriptFolder, cfg.outputStem + "_paired_bagged_results.xlsx");
outputMAT = fullfile(scriptFolder, cfg.outputStem + "_paired_bagged_results.mat");

if ~isfile(controlPath)
    error('fig5de:MissingControlFile', ...
        'Control workbook not found: %s', controlPath);
end
for g = 1:numel(experimentPaths)
    if ~isfile(experimentPaths(g))
        error('fig5de:MissingExperimentFile', ...
            'Experiment workbook not found: %s', experimentPaths(g));
    end
end

controlData = readPairedWideData(controlPath, cfg.controlMSheet, cfg.controlESheet, ...
    controlLabel);
% Preallocate from the identically shaped control struct so MATLAB accepts
% subsequent scalar-struct assignments without dissimilar-structure errors.
experimentData = repmat(controlData, numel(experimentPaths), 1);
for g = 1:numel(experimentPaths)
    experimentData(g) = readPairedWideData(experimentPaths(g), ...
        cfg.experimentMSheets(g), cfg.experimentESheets(g), experimentLabels(g));
end

fprintf('\n============================================================\n');
fprintf('Paired input-data summary: %s\n', controlLabel);
fprintf('============================================================\n');
printPairedDataSummary(controlData);

for g = 1:numel(experimentPaths)
    fprintf('\n============================================================\n');
    fprintf('Paired input-data summary: %s\n', experimentLabels(g));
    fprintf('============================================================\n');
    printPairedDataSummary(experimentData(g));

    if ~isequal(experimentData(g).time, controlData.time)
        warning('fig5de:DifferentTimeGrids', ...
            ['%s and control use different sampled-time grids. Each ' ...
             'bagged estimator remains on its own observed grid; no ' ...
             'interpolation or parametric curve fitting is performed.'], ...
            experimentTags(g));
    end
end

%% ========================================================================
% PART 2. Paired inner/outer bagging and M/E peak-shift METs
% ========================================================================
% One joint paired resampling run per biological group produces separate M
% and E bagged-peak sampling distributions. Thus, there are two inferential
% contrasts per experiment, while the genuine M/E sample pairing is retained.

fprintf('\nRunning control paired bagging: nInner=%d, nOuter=%d...\n', ...
    nInner, nOuter);
controlResult = runPairedGroupBagging(controlData, nInner, nOuter, ...
    rngSeed + 1, controlLabel);

experimentResults = repmat(controlResult, numel(experimentPaths), 1);

analysisStart = tic;
for g = 1:numel(experimentPaths)
    fprintf('\nRunning %s paired bagging: nInner=%d, nOuter=%d...\n', ...
        experimentTags(g), nInner, nOuter);
    experimentResults(g) = runPairedGroupBagging( ...
        experimentData(g), nInner, nOuter, ...
        rngSeed + 2000 + 100 * g, experimentLabels(g));
end

% Preallocate contrast struct arrays from the first fully formed contrast.
firstContrastM = computePeakShiftContrast(experimentResults(1), ...
    controlResult, "M", experimentTags(1), delta0, alpha);
firstContrastE = computePeakShiftContrast(experimentResults(1), ...
    controlResult, "E", experimentTags(1), delta0, alpha);
contrastM = repmat(firstContrastM, numel(experimentPaths), 1);
contrastE = repmat(firstContrastE, numel(experimentPaths), 1);

for g = 1:numel(experimentPaths)
    if g > 1
        currentContrastM = computePeakShiftContrast( ...
            experimentResults(g), controlResult, "M", ...
            experimentTags(g), delta0, alpha);
        currentContrastE = computePeakShiftContrast( ...
            experimentResults(g), controlResult, "E", ...
            experimentTags(g), delta0, alpha);
        contrastM(g) = currentContrastM;
        contrastE(g) = currentContrastE;
    end

    reportContrast(contrastM(g), experimentLabels(g), controlLabel);
    reportContrast(contrastE(g), experimentLabels(g), controlLabel);
end
totalElapsedSeconds = toc(analysisStart) + ...
    controlResult.outerElapsedSeconds;

% Warn without changing the requested Wald MET when an outer contrast is
% conspicuously non-normal.
for g = 1:numel(experimentPaths)
    if contrastM(g).stronglyNonNormal
        warning('fig5de:StronglyNonNormalMContrast', ...
            ['%s M outer-bootstrap contrast is strongly non-normal; ' ...
             'consider a studentized-bootstrap / bootstrap-t MET as a ' ...
             'sensitivity analysis. Wald MET remains the primary test.'], ...
            experimentTags(g));
    end
    if contrastE(g).stronglyNonNormal
        warning('fig5de:StronglyNonNormalEContrast', ...
            ['%s E outer-bootstrap contrast is strongly non-normal; ' ...
             'consider a studentized-bootstrap / bootstrap-t MET as a ' ...
             'sensitivity analysis. Wald MET remains the primary test.'], ...
            experimentTags(g));
    end
end

fprintf('\nTotal bagging elapsed time: %.1f seconds\n', ...
    totalElapsedSeconds);

%% ========================================================================
% PART 2B. Structured MATLAB, Excel, and MAT outputs
% ========================================================================
results = struct();
results.control = controlResult;
for g = 1:numel(experimentPaths)
    results.(experimentTags(g)) = experimentResults(g);
end

contrasts = struct();
for g = 1:numel(experimentPaths)
    contrasts.(experimentTags(g)).M = contrastM(g);
    contrasts.(experimentTags(g)).E = contrastE(g);
end

% Peak summary: one row per group and cell type.
peakRows = 2 * (1 + numel(experimentPaths));
peakGroup = strings(peakRows, 1);
peakRole = strings(peakRows, 1);
peakCellType = strings(peakRows, 1);
hardPeak = NaN(peakRows, 1);
baggedPeak = NaN(peakRows, 1);
outerSD = NaN(peakRows, 1);
CI_low = NaN(peakRows, 1);
CI_high = NaN(peakRows, 1);
peakConfidence = NaN(peakRows, 1);
peakEntropy = NaN(peakRows, 1);
peakStability = NaN(peakRows, 1);
innerHardPeakSD = NaN(peakRows, 1);
innerMonteCarloSE = NaN(peakRows, 1);

row = 0;
[row, peakGroup, peakRole, peakCellType, hardPeak, baggedPeak, ...
    outerSD, CI_low, CI_high, peakConfidence, peakEntropy, ...
    peakStability, innerHardPeakSD, innerMonteCarloSE] = ...
    addPeakRows(row, peakGroup, peakRole, peakCellType, hardPeak, ...
    baggedPeak, outerSD, CI_low, CI_high, peakConfidence, ...
    peakEntropy, peakStability, innerHardPeakSD, innerMonteCarloSE, ...
    controlTag, "Control", controlData, controlResult);

for g = 1:numel(experimentPaths)
    [row, peakGroup, peakRole, peakCellType, hardPeak, baggedPeak, ...
        outerSD, CI_low, CI_high, peakConfidence, peakEntropy, ...
        peakStability, innerHardPeakSD, innerMonteCarloSE] = ...
        addPeakRows(row, peakGroup, peakRole, peakCellType, hardPeak, ...
        baggedPeak, outerSD, CI_low, CI_high, peakConfidence, ...
        peakEntropy, peakStability, innerHardPeakSD, ...
        innerMonteCarloSE, experimentTags(g), "Experiment", ...
        experimentData(g), experimentResults(g));
end

peakSummaryTable = table(peakGroup, peakRole, peakCellType, hardPeak, ...
    baggedPeak, outerSD, CI_low, CI_high, peakConfidence, peakEntropy, ...
    peakStability, innerHardPeakSD, innerMonteCarloSE, ...
    'VariableNames', {'Group', 'Role', 'CellType', 'HardPeak_h', ...
    'BaggedPeak_h', 'OuterSD_h', 'CI95_low_h', 'CI95_high_h', ...
    'PeakConfidence', 'PeakEntropy', 'PeakStability', ...
    'InnerHardPeakSD_h', 'InnerMonteCarloSE_h'});

% Contrast summary: one row for M and one for E in every experiment.
nContrasts = 2 * numel(experimentPaths);
experimentColumn = strings(nContrasts, 1);
cellTypeColumn = strings(nContrasts, 1);
experimentPeak = NaN(nContrasts, 1);
controlPeak = NaN(nContrasts, 1);
peakShift = NaN(nContrasts, 1);
peakShiftSE = NaN(nContrasts, 1);
shiftCI_low = NaN(nContrasts, 1);
shiftCI_high = NaN(nContrasts, 1);
delta0Column = NaN(nContrasts, 1);
thresholdExcess = NaN(nContrasts, 1);
standardizedEffect = NaN(nContrasts, 1);
standardizedThresholdExcess = NaN(nContrasts, 1);
MET_z = NaN(nContrasts, 1);
MET_p = NaN(nContrasts, 1);
MET_z_near = NaN(nContrasts, 1);
MET_z_far = NaN(nContrasts, 1);
boundaryZ_lower = NaN(nContrasts, 1);
boundaryZ_upper = NaN(nContrasts, 1);
bootstrapSupport = NaN(nContrasts, 1);
outerSkewness = NaN(nContrasts, 1);
outerKurtosis = NaN(nContrasts, 1);
JB_p = NaN(nContrasts, 1);
stronglyNonNormal = false(nContrasts, 1);
decision = strings(nContrasts, 1);

row = 0;
for g = 1:numel(experimentPaths)
    for cellType = ["M", "E"]
        row = row + 1;
        if cellType == "M"
            c = contrastM(g);
        else
            c = contrastE(g);
        end
        experimentColumn(row) = experimentTags(g);
        cellTypeColumn(row) = cellType;
        experimentPeak(row) = c.experimentPeak;
        controlPeak(row) = c.controlPeak;
        peakShift(row) = c.diff;
        peakShiftSE(row) = c.SE;
        shiftCI_low(row) = c.CI_low;
        shiftCI_high(row) = c.CI_high;
        delta0Column(row) = c.delta0;
        thresholdExcess(row) = c.thresholdExcess;
        standardizedEffect(row) = c.standardizedEffect;
        standardizedThresholdExcess(row) = ...
            c.standardizedThresholdExcess;
        MET_z(row) = c.MET_z;
        MET_p(row) = c.MET_p;
        MET_z_near(row) = c.MET_z_near;
        MET_z_far(row) = c.MET_z_far;
        boundaryZ_lower(row) = c.boundaryZ_lower;
        boundaryZ_upper(row) = c.boundaryZ_upper;
        bootstrapSupport(row) = c.bootstrapSupport;
        outerSkewness(row) = c.outerSkewness;
        outerKurtosis(row) = c.outerKurtosis;
        JB_p(row) = c.JB_p;
        stronglyNonNormal(row) = c.stronglyNonNormal;
        decision(row) = c.MET_decision;
    end
end

contrastSummaryTable = table(experimentColumn, cellTypeColumn, ...
    experimentPeak, controlPeak, peakShift, peakShiftSE, ...
    shiftCI_low, shiftCI_high, delta0Column, thresholdExcess, ...
    standardizedEffect, standardizedThresholdExcess, MET_z, MET_p, ...
    MET_z_near, MET_z_far, boundaryZ_lower, boundaryZ_upper, ...
    bootstrapSupport, outerSkewness, outerKurtosis, JB_p, ...
    stronglyNonNormal, decision, 'VariableNames', { ...
    'Experiment', 'CellType', 'ExperimentPeak_h', 'ControlPeak_h', ...
    'PeakShift_exp_minus_ctrl_h', 'PeakShift_SE_h', 'CI95_low_h', ...
    'CI95_high_h', 'Delta0_h', 'ThresholdExcess_h', ...
    'StandardizedPeakShift', 'StandardizedThresholdExcess', ...
    'MET_Z', 'MET_p', 'MET_Z_near', 'MET_Z_far', ...
    'BoundaryZ_lower', 'BoundaryZ_upper', ...
    'BootstrapSupport_descriptive', 'OuterSkewness', ...
    'OuterKurtosis', 'JB_p', 'StrongNonNormalityFlag', ...
    'MET_decision'});

% One wide outer table keeps every replicate used for reproducibility.
outerBootstrapTable = table((1:nOuter).', ...
    controlResult.outerTbag_M, controlResult.outerTbag_E, ...
    'VariableNames', {'OuterReplicate', 'ControlTbag_M_h', ...
    'ControlTbag_E_h'});
for g = 1:numel(experimentPaths)
    tag = char(experimentTags(g));
    outerBootstrapTable.([tag '_Tbag_M_h']) = ...
        experimentResults(g).outerTbag_M;
    outerBootstrapTable.([tag '_Tbag_E_h']) = ...
        experimentResults(g).outerTbag_E;
    outerBootstrapTable.([tag '_minusCtrl_M_h']) = ...
        contrastM(g).outerDiff;
    outerBootstrapTable.([tag '_minusCtrl_E_h']) = ...
        contrastE(g).outerDiff;
end

configurationTable = table(nInner, nOuter, delta0, alpha, rngSeed, ...
    mainFigureBinWidth, multiplicityAdjustment, true, false, ...
    smokeTestMode, ...
    string(controlPath), string(dataFolder), ...
    'VariableNames', {'N_inner', 'N_outer', 'Delta0_h', 'Alpha', ...
    'RNG_seed', 'MainFigureBinWidth_h', ...
    'MultiplicityAdjustment', ...
    'PairedResamplingWithinGroup', 'PairingBetweenGroups', ...
    'SmokeTestMode', 'ControlDataFile', 'DataFolder'});

methodsNotes = { ...
    'Item', 'Definition'; ...
    'Hard peak estimator', 'Sampled time attaining the maximum of a bootstrap mean curve'; ...
    'Bagged peak estimator', 'Mean of inner-bootstrap discrete hard-argmax peak times'; ...
    'Within-group sampling', 'At every time and bootstrap layer, M and E use the same resampled sample-column indices'; ...
    'Between-group sampling', 'Experiment and control are independent and are resampled independently'; ...
    'M contrast', 'Tbag_M(experiment) - Tbag_M(control)'; ...
    'E contrast', 'Tbag_E(experiment) - Tbag_E(control)'; ...
    'Positive contrast', 'Experiment peaks later than control'; ...
    'Negative contrast', 'Experiment peaks earlier than control'; ...
    'Outer bootstrap', 'Complete inner bagging is repeated inside every outer nonparametric bootstrap dataset'; ...
    'Sampling SE', 'SD of the relevant outer-bootstrap contrast; never divided by sqrt(nOuter)'; ...
    '95% CI', 'Percentile interval of the relevant outer-bootstrap contrast'; ...
    'MET null', 'H0: |peak shift| <= delta0'; ...
    'MET alternative', 'H1: |peak shift| > delta0'; ...
    'MET p-value', 'Wald normal approximation at the least-favorable / nearest interval-null boundary'; ...
    'Bootstrap support', 'Descriptive proportion only; not the reported frequentist p-value'; ...
    'Standardized effect', 'Peak shift divided by its outer-bootstrap SE; not classical Cohen''s d'; ...
    'Time grids', 'Every group retains its observed sampling grid; no interpolation or parametric peak model is used'; ...
    'Histogram width', '1-h display bins for outer bagged peaks; not raw experimental time resolution'; ...
    'Figure colors', 'M=orange [255 128 0]; E=blue [0 0 255]; control alpha=0.3; experiment alpha=0.8'; ...
    'Multiplicity', char(multiplicityAdjustment); ...
    'Monte Carlo counts', 'nInner and nOuter affect numerical precision, not biological sample size'};

% Write result sheets. Unrelated files are untouched.
writetable(contrastSummaryTable, outputExcel, ...
    'Sheet', 'ContrastSummary', 'WriteMode', 'overwritesheet');
writetable(peakSummaryTable, outputExcel, ...
    'Sheet', 'PeakSummary', 'WriteMode', 'overwritesheet');
writetable(outerBootstrapTable, outputExcel, ...
    'Sheet', 'OuterBootstrap', 'WriteMode', 'overwritesheet');
writetable(configurationTable, outputExcel, ...
    'Sheet', 'Configuration', 'WriteMode', 'overwritesheet');
writecell(methodsNotes, outputExcel, ...
    'Sheet', 'MethodsNotes', 'WriteMode', 'overwritesheet');

writeGroupExcelSheets(outputExcel, controlTag, controlData, controlResult);
for g = 1:numel(experimentPaths)
    writeGroupExcelSheets(outputExcel, experimentTags(g), ...
        experimentData(g), experimentResults(g));
end

save(outputMAT, 'results', 'contrasts', 'controlData', ...
    'experimentData', 'nInner', 'nOuter', 'delta0', 'alpha', ...
    'rngSeed', 'mainFigureBinWidth');

fprintf('\nExcel results saved to:\n  %s\n', outputExcel);
fprintf('MAT results saved to:\n  %s\n', outputMAT);

%% ========================================================================
% PART 3. Four main figures: M and E peak distributions for MPDFR/EPDFR
% ========================================================================
% Each plot uses the outer-bootstrap distribution of the bagged estimator.
% Control is pale; the experiment uses the same saturated cell-type color.
% These are peak-time sampling distributions, not inner hard-peak-location
% probabilities.

mainFigurePaths = strings(numel(experimentPaths), 2);
for g = 1:numel(experimentPaths)
        outputM = fullfile(scriptFolder, cfg.outputStem + "_" + ...
            experimentTags(g) + "_M_paired_bagged.svg");
        outputE = fullfile(scriptFolder, cfg.outputStem + "_" + ...
            experimentTags(g) + "_E_paired_bagged.svg");


    plotOuterPeakComparison(controlResult.outerTbag_M, ...
        experimentResults(g).outerTbag_M, controlLabel, ...
        experimentLabels(g), controlColor, M_color, ...
        mainFigureBinWidth, "M", outputM, cfg.plot);
    plotOuterPeakComparison(controlResult.outerTbag_E, ...
        experimentResults(g).outerTbag_E, controlLabel, ...
        experimentLabels(g), controlColor, E_color, ...
        mainFigureBinWidth, "E", outputE, cfg.plot);

    mainFigurePaths(g, :) = [outputM, outputE];
    fprintf('%s M figure saved to:\n  %s\n', experimentTags(g), outputM);
    fprintf('%s E figure saved to:\n  %s\n', experimentTags(g), outputE);
end

figurePathTable = table(experimentTags, mainFigurePaths(:, 1), ...
    mainFigurePaths(:, 2), 'VariableNames', ...
    {'Experiment', 'M_figure', 'E_figure'});
writetable(figurePathTable, outputExcel, ...
    'Sheet', 'FigureFiles', 'WriteMode', 'overwritesheet');

%% ========================================================================
% PART 4. Diagnostic figures: contrast distributions and normal Q-Q plots
% ========================================================================
for g = 1:numel(experimentPaths)
        diagnosticPath = fullfile(scriptFolder, cfg.outputStem + "_" + ...
            experimentTags(g) + "_paired_bagged_diagnostics.svg");


    plotContrastDiagnostics(contrastM(g), contrastE(g), ...
        experimentTags(g), M_color, E_color, delta0, ...
        nInner, diagnosticPath);
    fprintf('%s diagnostic figure saved to:\n  %s\n', ...
        experimentTags(g), diagnosticPath);
end

%% ========================================================================
% Local helper functions
% ========================================================================

publicResult = struct('control', controlResult, 'experiments', experimentResults, 'contrastM', contrastM, 'contrastE', contrastE, 'summary', contrastSummaryTable);
publicResult.configuration = cfg;
if exist('outputMAT', 'var') && isfile(outputMAT)
    save(outputMAT, 'cfg', 'publicResult', '-append');
else
    save(fullfile(cfg.outputDir, cfg.outputStem + "_results.mat"), 'cfg', 'publicResult');
end
end

function data = readPairedWideData(fileName, M_sheet, E_sheet, groupLabel)
% Read paired M/E wide-format sheets.
% Row 1 contains sample IDs; column 1 contains sampled times; numeric values
% occupy rows 2:end and columns 2:end. M/E columns are matched by sample ID.

    M_cells = bootstrapReadCells(fileName,M_sheet);
    E_cells = bootstrapReadCells(fileName,E_sheet);

    if size(M_cells, 1) < 2 || size(M_cells, 2) < 2 || ...
            size(E_cells, 1) < 2 || size(E_cells, 2) < 2
        error('fig5de:EmptyWideSheet', ...
            '%s contains an empty or malformed M/E sheet.', groupLabel);
    end

    M_time = parseTimeCells(M_cells(2:end, 1));
    E_time = parseTimeCells(E_cells(2:end, 1));
    M_values = cellBlockToDouble(M_cells(2:end, 2:end));
    E_values = cellBlockToDouble(E_cells(2:end, 2:end));

    M_headers = normalizeHeaders(M_cells(1, 2:end));
    E_headers = normalizeHeaders(E_cells(1, 2:end));

    if numel(M_time) ~= numel(E_time) || ...
            any(~sameNumericWithNaN(M_time, E_time))
        error('fig5de:TimeGridMismatchWithinFile', ...
            '%s M and E sheets do not have identical time rows.', groupLabel);
    end
    if size(M_values, 2) ~= size(E_values, 2) || ...
            numel(M_headers) ~= numel(E_headers) || ...
            any(M_headers ~= E_headers)
        error('fig5de:SampleHeaderMismatch', ...
            '%s M and E sample columns are not identically aligned.', ...
            groupLabel);
    end

    finiteTime = isfinite(M_time) & isfinite(E_time);
    if ~all(finiteTime)
        warning('fig5de:InvalidTimeRows', ...
            '%s contains nonnumeric time rows; these rows are excluded.', ...
            groupLabel);
    end
    time = M_time(finiteTime);
    M_values = M_values(finiteTime, :);
    E_values = E_values(finiteTime, :);

    if numel(unique(time)) ~= numel(time)
        error('fig5de:DuplicateTimeRows', ...
            '%s contains duplicate sampled-time rows.', groupLabel);
    end
    [time, order] = sort(time);
    M_values = M_values(order, :);
    E_values = E_values(order, :);

    nTime = numel(time);
    M = cell(1, nTime);
    E = cell(1, nTime);
    nPairs = NaN(nTime, 1);
    mean_M = NaN(nTime, 1);
    sd_M = NaN(nTime, 1);
    mean_E = NaN(nTime, 1);
    sd_E = NaN(nTime, 1);
    pairedCorrelation = NaN(nTime, 1);

    for j = 1:nTime
        completePair = isfinite(M_values(j, :)) & ...
            isfinite(E_values(j, :));
        incompletePair = xor(isfinite(M_values(j, :)), ...
            isfinite(E_values(j, :)));
        if any(incompletePair)
            warning('fig5de:IncompletePair', ...
                ['%s time %.4g contains incomplete M/E sample columns. ' ...
                 'These columns are excluded at this time so bootstrap ' ...
                 'observations remain paired.'], groupLabel, time(j));
        end

        Mj = M_values(j, completePair).';
        Ej = E_values(j, completePair).';
        if isempty(Mj)
            error('fig5de:NoPairsAtTime', ...
                '%s has no complete M/E pairs at time %.4g.', ...
                groupLabel, time(j));
        end

        M{j} = Mj;
        E{j} = Ej;
        nPairs(j) = numel(Mj);
        mean_M(j) = mean(Mj);
        sd_M(j) = std(Mj);
        mean_E(j) = mean(Ej);
        sd_E(j) = std(Ej);
        if numel(Mj) > 2 && std(Mj) > 0 && std(Ej) > 0
            correlationMatrix = corrcoef(Mj, Ej);
            pairedCorrelation(j) = correlationMatrix(1, 2);
        end
    end

    [~, hardIndex_M] = max(mean_M);
    [~, hardIndex_E] = max(mean_E);

    data = struct();
    data.label = groupLabel;
    data.fileName = string(fileName);
    data.time = time(:);
    data.sampleHeaders = M_headers(:).';
    data.M = M;
    data.E = E;
    data.nPairs = nPairs;
    data.mean_M = mean_M;
    data.mean_E = mean_E;
    data.hardPeak_M = time(hardIndex_M);
    data.hardPeak_E = time(hardIndex_E);
    data.timepointTable = table(time(:), nPairs, mean_M, sd_M, ...
        mean_E, sd_E, pairedCorrelation, 'VariableNames', ...
        {'Time_h', 'N_pairs', 'Mean_M', 'SD_M', 'Mean_E', 'SD_E', ...
         'PairedCorrelation_M_E'});
end

function time = parseTimeCells(timeCells)
    time = NaN(numel(timeCells), 1);
    for j = 1:numel(timeCells)
        value = timeCells{j};
        if isnumeric(value) && isscalar(value) && isfinite(value)
            time(j) = double(value);
        else
            valueString = string(value);
            if ~ismissing(valueString)
                token = regexp(char(valueString), ...
                    '[-+]?\d*\.?\d+', 'match', 'once');
                if ~isempty(token)
                    time(j) = str2double(token);
                end
            end
        end
    end
end

function values = cellBlockToDouble(cellBlock)
    values = NaN(size(cellBlock));
    for k = 1:numel(cellBlock)
        item = cellBlock{k};
        if isnumeric(item) && isscalar(item)
            values(k) = double(item);
        elseif islogical(item) && isscalar(item)
            values(k) = double(item);
        else
            itemString = string(item);
            if ~ismissing(itemString)
                converted = str2double(itemString);
                if isfinite(converted)
                    values(k) = converted;
                end
            end
        end
    end
end

function headers = normalizeHeaders(headerCells)
    headers = strings(size(headerCells));
    for k = 1:numel(headerCells)
        item = headerCells{k};
        if isnumeric(item) && isscalar(item)
            headers(k) = compose('%.15g', double(item));
        else
            headers(k) = strtrim(string(item));
        end
    end
end

function equalMask = sameNumericWithNaN(a, b)
    equalMask = (a == b) | (isnan(a) & isnan(b));
end

function printPairedDataSummary(data)
    fprintf('%-8s %-8s %-12s %-12s %-12s\n', ...
        'Time', 'nPair', 'Mean M', 'Mean E', 'Pair r');
    for j = 1:numel(data.time)
        fprintf('%-8.4g %-8d %-12.4f %-12.4f %-12.4f\n', ...
            data.time(j), data.nPairs(j), data.mean_M(j), ...
            data.mean_E(j), ...
            data.timepointTable.PairedCorrelation_M_E(j));
    end
    fprintf('Observed hard peaks: M=%.4g h, E=%.4g h\n', ...
        data.hardPeak_M, data.hardPeak_E);
end

function result = runPairedGroupBagging(data, nInner, nOuter, ...
    seed, groupLabel)
% Original-data inner bagging followed by the outer nonparametric bootstrap.

    rng(seed, 'twister');
    [Tbag_M, Tbag_E, innerPeakTimes_M, innerPeakTimes_E, P_M, P_E] = ...
        computePairedBaggedPeaks(data.M, data.E, data.time, nInner);

    outerTbag_M = NaN(nOuter, 1);
    outerTbag_E = NaN(nOuter, 1);
    rng(seed + 1, 'twister');
    progressStep = max(1, floor(nOuter / 20));
    outerStart = tic;

    for b = 1:nOuter
        [outerData_M, outerData_E] = bootstrapPairedDataset( ...
            data.M, data.E);
        [outerTbag_M(b), outerTbag_E(b)] = ...
            computePairedBaggedPeaks(outerData_M, outerData_E, ...
            data.time, nInner);

        if mod(b, progressStep) == 0 || b == nOuter
            fprintf('  %s outer bootstrap: %d/%d complete\n', ...
                groupLabel, b, nOuter);
        end
    end
    outerElapsedSeconds = toc(outerStart);

    valid = isfinite(outerTbag_M) & isfinite(outerTbag_E);
    if sum(valid) < 2
        error('fig5de:InsufficientOuterReplicates', ...
            '%s has fewer than two valid outer replicates.', groupLabel);
    end

    [stability_M, entropy_M] = computePeakStability(P_M);
    [stability_E, entropy_E] = computePeakStability(P_E);

    result = struct();
    result.label = groupLabel;
    result.Tbag_M = Tbag_M;
    result.Tbag_E = Tbag_E;
    result.innerPeakTimes_M = innerPeakTimes_M;
    result.innerPeakTimes_E = innerPeakTimes_E;
    result.peakProbability_M = P_M;
    result.peakProbability_E = P_E;
    result.peakConfidence_M = max(P_M);
    result.peakConfidence_E = max(P_E);
    result.peakEntropy_M = entropy_M;
    result.peakEntropy_E = entropy_E;
    result.peakStability_M = stability_M;
    result.peakStability_E = stability_E;
    result.innerHardPeakSD_M = std(innerPeakTimes_M, 'omitnan');
    result.innerHardPeakSD_E = std(innerPeakTimes_E, 'omitnan');
    result.innerMonteCarloSE_M = result.innerHardPeakSD_M / ...
        sqrt(sum(isfinite(innerPeakTimes_M)));
    result.innerMonteCarloSE_E = result.innerHardPeakSD_E / ...
        sqrt(sum(isfinite(innerPeakTimes_E)));
    result.outerTbag_M = outerTbag_M;
    result.outerTbag_E = outerTbag_E;
    result.validOuter = valid;
    result.outerSD_M = std(outerTbag_M(valid));
    result.outerSD_E = std(outerTbag_E(valid));
    result.CI95_M = prctile(outerTbag_M(valid), [2.5 97.5]);
    result.CI95_E = prctile(outerTbag_E(valid), [2.5 97.5]);
    result.outerElapsedSeconds = outerElapsedSeconds;
end

function contrast = computePeakShiftContrast(experimentResult, ...
    controlResult, cellType, experimentTag, delta0, alpha)
% Independent-group outer contrast for one cell type.

    if cellType == "M"
        experimentPeak = experimentResult.Tbag_M;
        controlPeak = controlResult.Tbag_M;
        experimentOuter = experimentResult.outerTbag_M;
        controlOuter = controlResult.outerTbag_M;
    elseif cellType == "E"
        experimentPeak = experimentResult.Tbag_E;
        controlPeak = controlResult.Tbag_E;
        experimentOuter = experimentResult.outerTbag_E;
        controlOuter = controlResult.outerTbag_E;
    else
        error('fig5de:UnknownCellType', 'cellType must be M or E.');
    end

    outerDiff = experimentOuter - controlOuter;
    valid = isfinite(outerDiff);
    outerDiffValid = outerDiff(valid);
    diffEstimate = experimentPeak - controlPeak;
    bootstrapSE = std(outerDiffValid);
    CI95 = prctile(outerDiffValid, [2.5 97.5]);

    [MET_z, MET_p, z_near, z_far, boundaryZ_lower, ...
        boundaryZ_upper] = computeWaldMET(diffEstimate, ...
        bootstrapSE, delta0, "survival");

    outerSkewness = skewness(outerDiffValid, 0);
    outerKurtosis = kurtosis(outerDiffValid, 0);
    if isfinite(outerKurtosis) && outerKurtosis > 0
        bimodalityCoefficient = ...
            (outerSkewness ^ 2 + 1) / outerKurtosis;
    else
        bimodalityCoefficient = NaN;
    end

    if numel(outerDiffValid) >= 8 && bootstrapSE > 0
        [JB_reject, JB_p] = jbtest(outerDiffValid, alpha);
    else
        JB_reject = NaN;
        JB_p = NaN;
    end

    stronglyNonNormal = abs(outerSkewness) > 1 || ...
        outerKurtosis > 5 || ...
        (isfinite(bimodalityCoefficient) && ...
         bimodalityCoefficient > 5 / 9);

    if isfinite(MET_p) && MET_p < alpha
        decision = "Reject H0: absolute peak shift exceeds delta0";
    else
        decision = ...
            "Do not reject H0: insufficient evidence beyond delta0";
    end

    contrast = struct();
    contrast.experiment = experimentTag;
    contrast.cellType = cellType;
    contrast.definition = "Tbag_" + cellType + ...
        "(experiment) - Tbag_" + cellType + "(control)";
    contrast.experimentPeak = experimentPeak;
    contrast.controlPeak = controlPeak;
    contrast.diff = diffEstimate;
    contrast.SE = bootstrapSE;
    contrast.CI95 = CI95;
    contrast.CI_low = CI95(1);
    contrast.CI_high = CI95(2);
    contrast.delta0 = delta0;
    contrast.thresholdExcess = abs(diffEstimate) - delta0;
    if isfinite(bootstrapSE) && bootstrapSE > 0
        contrast.standardizedEffect = diffEstimate / bootstrapSE;
        contrast.standardizedThresholdExcess = ...
            contrast.thresholdExcess / bootstrapSE;
    else
        contrast.standardizedEffect = NaN;
        contrast.standardizedThresholdExcess = NaN;
    end
    contrast.MET_z = MET_z;
    contrast.MET_p = MET_p;
    contrast.MET_z_near = z_near;
    contrast.MET_z_far = z_far;
    contrast.boundaryZ_lower = boundaryZ_lower;
    contrast.boundaryZ_upper = boundaryZ_upper;
    contrast.bootstrapSupport = ...
        mean(abs(outerDiffValid) > delta0);
    contrast.outerDiff = outerDiff;
    contrast.validOuter = valid;
    contrast.outerSkewness = outerSkewness;
    contrast.outerKurtosis = outerKurtosis;
    contrast.bimodalityCoefficient = bimodalityCoefficient;
    contrast.JB_reject = JB_reject;
    contrast.JB_p = JB_p;
    contrast.stronglyNonNormal = stronglyNonNormal;
    contrast.MET_decision = decision;
end

function reportContrast(c, experimentLabel, controlLabel)
    fprintf('\n====================================================================\n');
    fprintf('%s %s peak versus %s\n', ...
        experimentLabel, c.cellType, controlLabel);
    fprintf('====================================================================\n');
    fprintf('Experiment bagged peak: %.4f h\n', c.experimentPeak);
    fprintf('Control bagged peak:    %.4f h\n', c.controlPeak);
    fprintf('Peak shift (experiment-control): %.4f +/- %.4f h\n', ...
        c.diff, c.SE);
    fprintf('95%% outer-bootstrap CI: [%.4f, %.4f] h\n', ...
        c.CI_low, c.CI_high);
    fprintf('H0: |peak shift| <= %.4f h\n', c.delta0);
    fprintf('H1: |peak shift| >  %.4f h\n', c.delta0);
    fprintf('MET Z: %.4f; MET p-value: %.6g\n', c.MET_z, c.MET_p);
    fprintf(['Bootstrap support beyond threshold: %.4f ' ...
        '(descriptive only)\n'], c.bootstrapSupport);
    fprintf('Decision: %s\n', c.MET_decision);
end

function [row, peakGroup, peakRole, peakCellType, hardPeak, ...
    baggedPeak, outerSD, CI_low, CI_high, peakConfidence, ...
    peakEntropy, peakStability, innerHardPeakSD, innerMonteCarloSE] = ...
    addPeakRows(row, peakGroup, peakRole, peakCellType, hardPeak, ...
    baggedPeak, outerSD, CI_low, CI_high, peakConfidence, ...
    peakEntropy, peakStability, innerHardPeakSD, innerMonteCarloSE, ...
    groupTag, role, data, result)

    row = row + 1;
    peakGroup(row) = groupTag;
    peakRole(row) = role;
    peakCellType(row) = "M";
    hardPeak(row) = data.hardPeak_M;
    baggedPeak(row) = result.Tbag_M;
    outerSD(row) = result.outerSD_M;
    CI_low(row) = result.CI95_M(1);
    CI_high(row) = result.CI95_M(2);
    peakConfidence(row) = result.peakConfidence_M;
    peakEntropy(row) = result.peakEntropy_M;
    peakStability(row) = result.peakStability_M;
    innerHardPeakSD(row) = result.innerHardPeakSD_M;
    innerMonteCarloSE(row) = result.innerMonteCarloSE_M;

    row = row + 1;
    peakGroup(row) = groupTag;
    peakRole(row) = role;
    peakCellType(row) = "E";
    hardPeak(row) = data.hardPeak_E;
    baggedPeak(row) = result.Tbag_E;
    outerSD(row) = result.outerSD_E;
    CI_low(row) = result.CI95_E(1);
    CI_high(row) = result.CI95_E(2);
    peakConfidence(row) = result.peakConfidence_E;
    peakEntropy(row) = result.peakEntropy_E;
    peakStability(row) = result.peakStability_E;
    innerHardPeakSD(row) = result.innerHardPeakSD_E;
    innerMonteCarloSE(row) = result.innerMonteCarloSE_E;
end

function writeGroupExcelSheets(outputExcel, groupTag, data, result)
    tag = char(groupTag);
    peakProbabilityTable = table(data.time(:), ...
        result.peakProbability_M(:), result.peakProbability_E(:), ...
        'VariableNames', {'Time_h', 'PeakProbability_M', ...
        'PeakProbability_E'});
    writetable(peakProbabilityTable, outputExcel, ...
        'Sheet', [tag '_PeakProb'], 'WriteMode', 'overwritesheet');
    writetable(data.timepointTable, outputExcel, ...
        'Sheet', [tag '_Timepoints'], 'WriteMode', 'overwritesheet');
end

function plotOuterPeakComparison(controlOuter, experimentOuter, ...
    controlLabel, experimentLabel, controlColor, experimentColor, ...
    binWidth, cellType, outputPath, plotConfig)

    controlOuter = controlOuter(isfinite(controlOuter));
    experimentOuter = experimentOuter(isfinite(experimentOuter));

    % Bins are centered at integer multiples of binWidth. This display
    % choice does not claim that raw observations were sampled at binWidth.
    xLimits = [10, 20];
    binEdges = (xLimits(1) - binWidth / 2) : binWidth : ...
        (xLimits(2) + binWidth / 2);

    [controlProb, ~] = histcounts(controlOuter, binEdges, ...
        'Normalization', 'probability');
    [experimentProb, ~] = histcounts(experimentOuter, binEdges, ...
        'Normalization', 'probability');
    maxProbability = max([controlProb(:); experimentProb(:)]);

    fig = figure('Name', sprintf('%s bagged peak comparison', cellType), ...
        'NumberTitle', 'off', 'Color', 'w', 'Units', 'centimeters', ...
        'Position', [2, 2, 3, 4]);
    hold on;
    % histogram(controlOuter, binEdges, 'FaceColor', controlColor, ...
    %     'EdgeColor', 'none', 'FaceAlpha', 0.8, ...
    %     'DisplayName', controlLabel, 'Normalization', 'probability');
    histogram(controlOuter, binEdges, 'FaceColor', experimentColor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
        'DisplayName', controlLabel, 'Normalization', 'probability');
    histogram(experimentOuter, binEdges, 'FaceColor', experimentColor, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8, ...
        'DisplayName', experimentLabel, 'Normalization', 'probability');
    hold off;

    ax = gca;
    ax.Units = 'centimeters';
    ax.Position(3:4) = plotConfig.axesSizeCm;
    ax.LineWidth = 0.8;
    ax.TickDir = 'out';
    ax.TickLength = [0.02, 0.02];
    ax.XTickLabelRotation = 0;
    ax.XRuler.TickLabelGapOffset = -1;
    ax.FontName = 'Arial';
    ax.FontSize = 7;
    ax.Color = 'none';

    % if maxProbability > 0.78
    %     ylim([0, 1]);
    %     yticks(0:0.5:1);
    % else
    %     ylim([0, 0.8]);
    %     yticks(0:0.4:0.8);
    % end
    ylim([0, 1]);
    yticks(0:0.5:1);
    xlim(xLimits);
    xticks(10:5:20);
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = setdiff(0:1:24, 10:5:20);
    box off;

    set(findall(fig, '-property', 'FontUnits'), ...
        'FontUnits', 'points');
    set(findall(fig, '-property', 'FontSize'), 'FontSize', 7);
    set(findall(fig, '-property', 'FontSizeMode'), ...
        'FontSizeMode', 'manual');
    set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');

    fig.PaperUnits = 'centimeters';
    fig.PaperSize = [3, 4];
    fig.PaperPosition = [0, 0, 3, 4];
    fig.PaperPositionMode = 'manual';
    fig.InvertHardcopy = 'off';
    drawnow;
    print(fig, outputPath, '-dsvg', '-vector');
end

function plotContrastDiagnostics(contrastM, contrastE, ...
    experimentTag, M_color, E_color, delta0, nInner, outputPath)

    fig = figure('Name', experimentTag + " peak-shift diagnostics", ...
        'NumberTitle', 'off', 'Color', 'w', 'Units', 'centimeters', ...
        'Position', [3, 3, 18, 14]);
    layout = tiledlayout(fig, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile(layout, 1);
    plotOneContrastDistribution(contrastM, M_color, delta0);
    title(sprintf('%s M shift', experimentTag));

    nexttile(layout, 2);
    valuesM = contrastM.outerDiff(contrastM.validOuter);
    qqplot(valuesM);
    xlabel('Standard-normal quantiles');
    ylabel('Outer M-shift quantiles (h)');
    title(sprintf('M Q-Q plot (JB p=%.3g)', contrastM.JB_p));
    box off;

    nexttile(layout, 3);
    plotOneContrastDistribution(contrastE, E_color, delta0);
    title(sprintf('%s E shift', experimentTag));

    nexttile(layout, 4);
    valuesE = contrastE.outerDiff(contrastE.validOuter);
    qqplot(valuesE);
    xlabel('Standard-normal quantiles');
    ylabel('Outer E-shift quantiles (h)');
    title(sprintf('E Q-Q plot (JB p=%.3g)', contrastE.JB_p));
    box off;

    sgtitle(layout, sprintf( ...
        'Outer-bootstrap peak-shift diagnostics (nInner=%d)', nInner));
    drawnow;
    print(fig, outputPath, '-dsvg', '-vector');
end

function plotOneContrastDistribution(c, experimentColor, delta0)
    values = c.outerDiff(c.validOuter);
    histogram(values, 30, 'Normalization', 'pdf', ...
        'FaceColor', [0.45 0.45 0.45], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.75, 'DisplayName', 'Outer contrast');
    hold on;
    if c.SE > 0 && range(values) > 0
        xGrid = linspace(min(values), max(values), 300);
        plot(xGrid, normpdf(xGrid, mean(values), c.SE), ...
            'Color', experimentColor, 'LineWidth', 1.2, ...
            'DisplayName', 'Normal fit');
    end
    xline(-delta0, '--k', '-delta0', 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
    xline(delta0, '--k', '+delta0', 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
    xline(c.diff, '-', 'Observed shift', 'Color', [0.75 0 0], ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');
    xlabel('Experiment - control peak shift (h)');
    ylabel('Density');
    legend('Location', 'best');
    box off;
end
