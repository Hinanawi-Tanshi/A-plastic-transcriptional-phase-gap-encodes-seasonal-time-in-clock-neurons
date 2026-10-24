function publicResult = runPairedBatchAnalysis(cfg)
% Shared publication pipeline; settings are supplied by the figure entry script.
figureCleanup = bootstrapFigureContext(cfg); %#ok<NASGU>
folderPath = cfg.outputDir;
outputFolder = cfg.outputDir;
dataFolder = cfg.dataFolder;
filePattern = cfg.filePattern;
nInner = cfg.nInner;
nOuter = cfg.nOuter;
delta0 = cfg.delta0;
alpha = cfg.alpha;
rngSeedBase = cfg.rngSeedBase;
meanSEMCondition = cfg.meanSEMCondition;
meanSEMUnit = cfg.meanSEMUnit;
mainFigureBinWidth = cfg.mainFigureBinWidth;
diagnosticDiffBinWidth = cfg.diagnosticDiffBinWidth;
mainXLimits = cfg.mainXLimits;
if cfg.showFigures, figureVisibility = 'on'; else, figureVisibility = 'off'; end
outputSuffix = "";
outputExcel = fullfile(outputFolder, cfg.outputStem + "_results.xlsx");
outputMeanSEMExcel = fullfile(outputFolder, cfg.outputStem + "_mean_sem.xlsx");

% The mean+/-SEM workbook is rebuilt from scratch on every run so a prior
% run with a different file subset cannot leave stale condition sheets.
if isfile(outputMeanSEMExcel)
    delete(outputMeanSEMExcel);
end

datasetNames = string(cfg.datasetNames(:));
files = repmat(struct('name','','folder',char(dataFolder)), numel(datasetNames), 1);
fileNumbers = zeros(numel(datasetNames),1);
for k = 1:numel(datasetNames)
    files(k).name = char(datasetNames(k) + ".xlsx");
    fileNumbers(k) = str2double(extractAfter(datasetNames(k), 'DL'));
end
% Names are specified in the original natural numeric order to preserve seeds.
assert(issorted(fileNumbers), 'Dataset names must be in natural numeric order.');
nFiles = numel(files);

fprintf('====================================================================\n');
fprintf('Fig. S10 paired batch bagged-bootstrap + Wald MET\n');
fprintf('Files: %d | nInner: %d | nOuter: %d | delta0: %.3f h\n', ...
    nFiles, nInner, nOuter, delta0);
fprintf('Input: %s\n', dataFolder);
fprintf('====================================================================\n');

% Per-file tables are collected in cells and concatenated once, avoiding
% repeated growth of large outer-bootstrap tables.
summaryRows = repmat(struct(), nFiles, 1);
outerTables = cell(nFiles, 1);
hardBootstrapTables = cell(nFiles, 1);
peakProbabilityTables = cell(nFiles, 1);
timepointTables = cell(nFiles, 1);
diagnosticRows = repmat(struct(), nFiles, 1);
batchResults = repmat(struct(), nFiles, 1);

batchStart = tic;

%% ========================================================================
% PART 2. Fully paired inner/outer bagging bootstrap and Wald-type MET
% ========================================================================
for f = 1:nFiles
    inputFile = cfg.inputFile;
    datasetName = datasetNames(f);
    datasetName = string(datasetName);
    fileSeed = rngSeedBase + 10000 * f;

    fprintf('\n[%d/%d] %s\n', f, nFiles, files(f).name);
    fprintf('------------------------------------------------------------\n');

    [time, M_perTime, E_perTime, timepointTable, dataInfo] = ...
        readPairedWideWorkbook(inputFile,cfg.MSheets(f),cfg.ESheets(f));
    timepointTable = addvars(timepointTable, ...
        repmat(datasetName, height(timepointTable), 1), ...
        'Before', 1, 'NewVariableNames', 'Dataset');
    timepointTables{f} = timepointTable;

    % Separate plotting-data export requested for downstream figure work.
    % n is the number of complete M/E pairs retained at each time point.
    % SEM describes the raw sample mean: SD/sqrt(n). It is unrelated to the
    % bootstrap SE of a hard or bagged peak estimator.
    meanSEMConditionColumn = repmat(meanSEMCondition, ...
        height(timepointTable), 1);
    meanSEMUnitColumn = repmat(meanSEMUnit, height(timepointTable), 1);
    meanSEM_M = table(meanSEMConditionColumn, timepointTable.Time_h, ...
        timepointTable.N_pairs, meanSEMUnitColumn, ...
        timepointTable.Mean_sLNv, ...
        timepointTable.SD_sLNv ./ sqrt(timepointTable.N_pairs), ...
        'VariableNames', {'condition', 'time', 'n', 'unit', 'mean', 'sem'});
    meanSEM_E = table(meanSEMConditionColumn, timepointTable.Time_h, ...
        timepointTable.N_pairs, meanSEMUnitColumn, ...
        timepointTable.Mean_5th_sLNv, ...
        timepointTable.SD_5th_sLNv ./ sqrt(timepointTable.N_pairs), ...
        'VariableNames', {'condition', 'time', 'n', 'unit', 'mean', 'sem'});

    % Source files are named DLxx, while the requested output convention is
    % LDxx_M / LDxx_E. All current inputs contain a numeric condition code.
    if isfinite(fileNumbers(f))
        meanSEMSheetBase = "LD" + string(fileNumbers(f));
    else
        meanSEMSheetBase = "LD_" + datasetName;
    end
    writetable(meanSEM_M, outputMeanSEMExcel, ...
        'Sheet', meanSEMSheetBase + "_M", ...
        'WriteMode', 'overwritesheet');
    writetable(meanSEM_E, outputMeanSEMExcel, ...
        'Sheet', meanSEMSheetBase + "_E", ...
        'WriteMode', 'overwritesheet');

    meanCurve_M = cellfun(@mean, M_perTime).';
    meanCurve_E = cellfun(@mean, E_perTime).';
    [~, hardIndex_M] = max(meanCurve_M);
    [~, hardIndex_E] = max(meanCurve_E);
    hardPeak_M = time(hardIndex_M);
    hardPeak_E = time(hardIndex_E);
    hardDiff = hardPeak_M - hardPeak_E;

    % Inner paired bootstrap: every replicate is a discrete hard argmax.
    % Their average is the bootstrap-smoothed / bagged estimator.
    rng(fileSeed + 1, 'twister');
    [Tbag_M, Tbag_E, innerPeakTimes_M, innerPeakTimes_E, P_M, P_E] = ...
        computePairedBaggedPeaks(M_perTime, E_perTime, time, nInner);
    baggedDiff = Tbag_M - Tbag_E;

    % Hard-argmax sensitivity analysis using the same paired single-layer
    % bootstrap replicates. The SD of bootstrap differences is already the
    % bootstrap SE; nInner is never used as an experimental sample size.
    bootHardDiff = innerPeakTimes_M - innerPeakTimes_E;
    validHard = isfinite(bootHardDiff);
    hardSE = std(bootHardDiff(validHard));
    hardCI95 = prctile(bootHardDiff(validHard), [2.5 97.5]);
    hardBootstrapSupport = mean(abs(bootHardDiff(validHard)) > delta0);
    [hardMET_z, hardMET_p, ~, ~, hardBoundaryZ_lower, ...
        hardBoundaryZ_upper] = computeWaldMET(hardDiff, hardSE, delta0);

    statsHard = struct();
    statsHard.peak_M = hardPeak_M;
    statsHard.peak_E = hardPeak_E;
    statsHard.diff = hardDiff;
    statsHard.SE = hardSE;
    statsHard.CI95 = hardCI95;
    statsHard.bootstrapSupport = hardBootstrapSupport;
    statsHard.MET_z = hardMET_z;
    statsHard.MET_p = hardMET_p;
    statsHard.boundaryZ_lower = hardBoundaryZ_lower;
    statsHard.boundaryZ_upper = hardBoundaryZ_upper;

    % Outer nonparametric bootstrap: one paired outer dataset is created,
    % then the complete inner-bagging estimator is recomputed on it.
    outerTbag_M = NaN(nOuter, 1);
    outerTbag_E = NaN(nOuter, 1);
    rng(fileSeed + 2, 'twister');
    outerStart = tic;
    progressStep = max(1, floor(nOuter / 10));
    for b = 1:nOuter
        [outerData_M, outerData_E] = ...
            bootstrapPairedDataset(M_perTime, E_perTime);
        [outerTbag_M(b), outerTbag_E(b)] = ...
            computePairedBaggedPeaks(outerData_M, outerData_E, ...
            time, nInner);

        if mod(b, progressStep) == 0 || b == nOuter
            fprintf('  outer bootstrap %d/%d\n', b, nOuter);
        end
    end
    outerElapsedSeconds = toc(outerStart);

    outerBaggedDiff = outerTbag_M - outerTbag_E;
    validOuter = isfinite(outerTbag_M) & isfinite(outerTbag_E) & ...
        isfinite(outerBaggedDiff);
    if sum(validOuter) < 2
        error('figs10:TooFewValidOuterReplicates', ...
            '%s produced fewer than two valid outer replicates.', ...
            files(f).name);
    end

    outerSD_M = std(outerTbag_M(validOuter));
    outerSD_E = std(outerTbag_E(validOuter));
    baggedSE = std(outerBaggedDiff(validOuter));
    CI95_M = prctile(outerTbag_M(validOuter), [2.5 97.5]);
    CI95_E = prctile(outerTbag_E(validOuter), [2.5 97.5]);
    CI95_diff = prctile(outerBaggedDiff(validOuter), [2.5 97.5]);

    [MET_z, MET_p, MET_z_near, MET_z_far, ...
        boundaryZ_lower, boundaryZ_upper] = ...
        computeWaldMET(baggedDiff, baggedSE, delta0);

    thresholdExcess = abs(baggedDiff) - delta0;
    if isfinite(baggedSE) && baggedSE > 0
        standardizedEffect = baggedDiff / baggedSE;
        standardizedThresholdExcess = thresholdExcess / baggedSE;
    else
        standardizedEffect = NaN;
        standardizedThresholdExcess = NaN;
    end
    bootstrapSupport = ...
        mean(abs(outerBaggedDiff(validOuter)) > delta0);

    [peakStability_M, peakEntropy_M] = computePeakStability(P_M);
    [peakStability_E, peakEntropy_E] = computePeakStability(P_E);
    peakConfidence_M = max(P_M, [], 'omitnan');
    peakConfidence_E = max(P_E, [], 'omitnan');

    % Monte Carlo diagnostics only. These sqrt terms describe simulation
    % integration precision and are never used for biological inference.
    innerMCSE_M = std(innerPeakTimes_M, 'omitnan') / ...
        sqrt(sum(isfinite(innerPeakTimes_M)));
    innerMCSE_E = std(innerPeakTimes_E, 'omitnan') / ...
        sqrt(sum(isfinite(innerPeakTimes_E)));
    outerSEMonteCarloErrorApprox = baggedSE / ...
        sqrt(2 * max(sum(validOuter) - 1, 1));

    if sum(validOuter) > 1
        outerCovarianceMatrix = cov(outerTbag_M(validOuter), ...
            outerTbag_E(validOuter));
        outerCorrelationMatrix = corrcoef(outerTbag_M(validOuter), ...
            outerTbag_E(validOuter));
        outerCovariance_M_E = outerCovarianceMatrix(1, 2);
        outerCorrelation_M_E = outerCorrelationMatrix(1, 2);
    else
        outerCovariance_M_E = NaN;
        outerCorrelation_M_E = NaN;
    end

    outerDiffValid = outerBaggedDiff(validOuter);
    outerSkewness = skewness(outerDiffValid, 0);
    outerKurtosis = kurtosis(outerDiffValid, 0); % normal = 3
    if isfinite(outerKurtosis) && outerKurtosis > 0
        outerBimodalityCoefficient = ...
            (outerSkewness ^ 2 + 1) / outerKurtosis;
    else
        outerBimodalityCoefficient = NaN;
    end

    if numel(outerDiffValid) >= 8 && baggedSE > 0
        [JB_reject, JB_p] = jbtest(outerDiffValid, alpha);
    else
        JB_reject = NaN;
        JB_p = NaN;
    end
    stronglyNonNormal = abs(outerSkewness) > 1 || ...
        outerBimodalityCoefficient > 0.70 || outerKurtosis > 6;
    if stronglyNonNormal
        warning('figs10:OuterBootstrapNonNormal', [ ...
            '%s: outer bootstrap distribution is strongly non-normal; ' ...
            'consider a studentized-bootstrap / bootstrap-t MET as a ' ...
            'sensitivity analysis.'], files(f).name);
    end

    if isfinite(MET_p) && MET_p < alpha
        MET_decision = "Reject H0: |Delta| exceeds delta0";
    else
        MET_decision = ...
            "Do not reject H0: minimum-effect evidence is insufficient";
    end

    statsBagged = struct();
    statsBagged.peak_M = Tbag_M;
    statsBagged.peak_E = Tbag_E;
    statsBagged.diff = baggedDiff;
    statsBagged.SE = baggedSE;
    statsBagged.CI95 = CI95_diff;
    statsBagged.CI_low = CI95_diff(1);
    statsBagged.CI_high = CI95_diff(2);
    statsBagged.MET_z = MET_z;
    statsBagged.MET_p = MET_p;
    statsBagged.MET_z_near = MET_z_near;
    statsBagged.MET_z_far = MET_z_far;
    statsBagged.boundaryZ_lower = boundaryZ_lower;
    statsBagged.boundaryZ_upper = boundaryZ_upper;
    statsBagged.standardizedEffect = standardizedEffect;
    statsBagged.thresholdExcess = thresholdExcess;
    statsBagged.standardizedThresholdExcess = ...
        standardizedThresholdExcess;
    statsBagged.peakProbability_M = P_M;
    statsBagged.peakProbability_E = P_E;
    statsBagged.peakConfidence_M = peakConfidence_M;
    statsBagged.peakConfidence_E = peakConfidence_E;
    statsBagged.peakStability_M = peakStability_M;
    statsBagged.peakStability_E = peakStability_E;
    statsBagged.bootstrapSupport = bootstrapSupport;
    statsBagged.nInner = nInner;
    statsBagged.nOuter = nOuter;

    % Store full distributions in the workspace result structure.
    batchResults(f).file = files(f).name;
    batchResults(f).time = time;
    batchResults(f).statsHard = statsHard;
    batchResults(f).statsBagged = statsBagged;
    batchResults(f).innerPeakTimes_M = innerPeakTimes_M;
    batchResults(f).innerPeakTimes_E = innerPeakTimes_E;
    batchResults(f).outerTbag_M = outerTbag_M;
    batchResults(f).outerTbag_E = outerTbag_E;
    batchResults(f).outerBaggedDiff = outerBaggedDiff;

    % One scalar summary row per workbook.
    summaryRows(f).Dataset = datasetName;
    summaryRows(f).SourceFile = string(files(f).name);
    summaryRows(f).DatasetNumber = fileNumbers(f);
    summaryRows(f).N_timepoints = numel(time);
    summaryRows(f).N_pairs_min = dataInfo.nPairsMin;
    summaryRows(f).N_pairs_max = dataInfo.nPairsMax;
    summaryRows(f).HardPeak_sLNv_h = hardPeak_M;
    summaryRows(f).HardPeak_5th_sLNv_h = hardPeak_E;
    summaryRows(f).HardDifference_M_minus_E_h = hardDiff;
    summaryRows(f).HardBootstrapSE_h = hardSE;
    summaryRows(f).HardCI95_low_h = hardCI95(1);
    summaryRows(f).HardCI95_high_h = hardCI95(2);
    summaryRows(f).HardMET_Z = hardMET_z;
    summaryRows(f).HardMET_p = hardMET_p;
    summaryRows(f).BaggedPeak_sLNv_h = Tbag_M;
    summaryRows(f).BaggedPeak_sLNv_outerSD_h = outerSD_M;
    summaryRows(f).BaggedPeak_sLNv_CI95_low_h = CI95_M(1);
    summaryRows(f).BaggedPeak_sLNv_CI95_high_h = CI95_M(2);
    summaryRows(f).BaggedPeak_5th_sLNv_h = Tbag_E;
    summaryRows(f).BaggedPeak_5th_sLNv_outerSD_h = outerSD_E;
    summaryRows(f).BaggedPeak_5th_sLNv_CI95_low_h = CI95_E(1);
    summaryRows(f).BaggedPeak_5th_sLNv_CI95_high_h = CI95_E(2);
    summaryRows(f).BaggedDifference_M_minus_E_h = baggedDiff;
    summaryRows(f).PhaseGap_E_minus_M_h = -baggedDiff;
    summaryRows(f).BaggedDifference_bootstrapSE_h = baggedSE;
    summaryRows(f).BaggedDifference_CI95_low_h = CI95_diff(1);
    summaryRows(f).BaggedDifference_CI95_high_h = CI95_diff(2);
    summaryRows(f).Delta0_h = delta0;
    summaryRows(f).ThresholdExcess_h = thresholdExcess;
    summaryRows(f).StandardizedPeakShift_not_Cohen_d = ...
        standardizedEffect;
    summaryRows(f).StandardizedThresholdExcess = ...
        standardizedThresholdExcess;
    summaryRows(f).BoundaryZ_lower = boundaryZ_lower;
    summaryRows(f).BoundaryZ_upper = boundaryZ_upper;
    summaryRows(f).MET_Z = MET_z;
    summaryRows(f).MET_p = MET_p;
    summaryRows(f).BootstrapSupport_descriptive = bootstrapSupport;
    summaryRows(f).PeakConfidence_sLNv = peakConfidence_M;
    summaryRows(f).PeakConfidence_5th_sLNv = peakConfidence_E;
    summaryRows(f).PeakStability_sLNv = peakStability_M;
    summaryRows(f).PeakStability_5th_sLNv = peakStability_E;
    summaryRows(f).OuterCorrelation_M_E = outerCorrelation_M_E;
    summaryRows(f).OuterSkewness_difference = outerSkewness;
    summaryRows(f).OuterKurtosis_difference = outerKurtosis;
    summaryRows(f).JB_p = JB_p;
    summaryRows(f).StrongNonNormalityFlag = stronglyNonNormal;
    summaryRows(f).N_inner = nInner;
    summaryRows(f).N_outer = nOuter;
    summaryRows(f).N_outer_valid = sum(validOuter);
    summaryRows(f).RNG_seed_file = fileSeed;
    summaryRows(f).MET_decision = MET_decision;

    diagnosticRows(f).Dataset = datasetName;
    diagnosticRows(f).OuterSkewness = outerSkewness;
    diagnosticRows(f).OuterKurtosis_normal3 = outerKurtosis;
    diagnosticRows(f).OuterBimodalityCoefficient = ...
        outerBimodalityCoefficient;
    diagnosticRows(f).JarqueBera_p = JB_p;
    diagnosticRows(f).JarqueBera_reject = JB_reject;
    diagnosticRows(f).StrongNonNormalityFlag = stronglyNonNormal;
    diagnosticRows(f).OuterCovariance_M_E = outerCovariance_M_E;
    diagnosticRows(f).OuterCorrelation_M_E = outerCorrelation_M_E;
    diagnosticRows(f).InnerMCSE_sLNv_h = innerMCSE_M;
    diagnosticRows(f).InnerMCSE_5th_sLNv_h = innerMCSE_E;
    diagnosticRows(f).ApproxMCerror_differenceSE_h = ...
        outerSEMonteCarloErrorApprox;
    diagnosticRows(f).OuterElapsedSeconds = outerElapsedSeconds;

    outerTables{f} = table(repmat(datasetName, nOuter, 1), ...
        (1:nOuter).', outerTbag_M, outerTbag_E, outerBaggedDiff, ...
        validOuter, 'VariableNames', {'Dataset', 'OuterReplicate', ...
        'OuterTbag_sLNv_h', 'OuterTbag_5th_sLNv_h', ...
        'OuterDifference_M_minus_E_h', 'ValidOuter'});

    hardBootstrapTables{f} = table(repmat(datasetName, nInner, 1), ...
        (1:nInner).', innerPeakTimes_M, innerPeakTimes_E, bootHardDiff, ...
        'VariableNames', {'Dataset', 'InnerReplicate', ...
        'HardPeak_sLNv_h', 'HardPeak_5th_sLNv_h', ...
        'HardDifference_M_minus_E_h'});

    peakProbabilityTables{f} = table( ...
        repmat(datasetName, numel(time), 1), time(:), P_M(:), P_E(:), ...
        'VariableNames', {'Dataset', 'Time_h', ...
        'PeakProbability_sLNv', 'PeakProbability_5th_sLNv'});

    fprintf('Bagged peaks: sLNv %.4f +/- %.4f h; 5th sLNv %.4f +/- %.4f h\n', ...
        Tbag_M, outerSD_M, Tbag_E, outerSD_E);
    fprintf('Difference M-E: %.4f +/- %.4f h; 95%% CI [%.4f, %.4f]\n', ...
        baggedDiff, baggedSE, CI95_diff(1), CI95_diff(2));
    fprintf('Wald MET: Z = %.4f, p = %.6g | %s\n', ...
        MET_z, MET_p, MET_decision);

    %% --------------------------------------------------------------------
    % PART 3. Per-file main figure in the original figs10.m visual style
    % ---------------------------------------------------------------------
    outputMainSVG = fullfile(outputFolder, ...
        datasetName + "_paired_bagged_MET" + outputSuffix + ".svg");
    createMainOuterFigure(outerTbag_M(validOuter), ...
        outerTbag_E(validOuter), mainFigureBinWidth, mainXLimits, ...
        figureVisibility, outputMainSVG, datasetName);

    %% --------------------------------------------------------------------
    % PART 4. Per-file outer-bootstrap diagnostics
    % ---------------------------------------------------------------------
    outputDiagnosticSVG = fullfile(outputFolder, ...
        datasetName + "_paired_bagged_MET_diagnostics" + ...
        outputSuffix + ".svg");
    createDiagnosticFigure(outerTbag_M(validOuter), ...
        outerTbag_E(validOuter), outerDiffValid, Tbag_M, Tbag_E, ...
        baggedDiff, baggedSE, delta0, JB_p, outerCorrelation_M_E, ...
        diagnosticDiffBinWidth, figureVisibility, ...
        outputDiagnosticSVG, datasetName);

    fprintf('Saved main figure: %s\n', outputMainSVG);
    fprintf('Saved diagnostics: %s\n', outputDiagnosticSVG);
end

%% ========================================================================
% PART 5. Combined Excel output for all workbooks
% ========================================================================
summaryTable = struct2table(summaryRows);
diagnosticTable = struct2table(diagnosticRows);
outerBootstrapTable = vertcat(outerTables{:});
hardBootstrapTable = vertcat(hardBootstrapTables{:});
peakProbabilityTable = vertcat(peakProbabilityTables{:});
timepointDataTable = vertcat(timepointTables{:});

configurationTable = table(string(dataFolder), string(outputFolder), ...
    string(filePattern), nFiles, nInner, nOuter, delta0, alpha, rngSeedBase, ...
    mainFigureBinWidth, diagnosticDiffBinWidth, true, ...
    'VariableNames', {'InputFolder', 'OutputFolder', 'FilePattern', 'N_files', ...
    'N_inner', 'N_outer', 'Delta0_h', 'Alpha', 'RNG_seed_base', ...
    'MainFigureBinWidth_h', 'DiagnosticDiffBinWidth_h', ...
    'PairedResampling'});

methodsNotes = { ...
    'Item', 'Definition'; ...
    'Sheet mapping', 'Sheet 1 = sLNv (M); Sheet 2 = 5th sLNv (E)'; ...
    'Pairing', 'The same row/column position across sheets is one biological pair; one common column index is used for M and E at both bootstrap layers'; ...
    'Time rows', 'Each sampled time row is resampled independently; columns are within-row sample indices, not longitudinal identities across time'; ...
    'Hard estimator', 'Sampled time attaining the maximum of the mean curve; ties use MATLAB max first-index behavior'; ...
    'Bagged estimator', 'Mean of inner-bootstrap discrete hard-argmax peak times'; ...
    'Uncertainty', 'Outer nonparametric bootstrap SD/SE and percentile confidence interval'; ...
    'Hard SE', 'SD of paired single-layer hard-peak bootstrap differences; never divided by sqrt(nInner)'; ...
    'Bagged SE', 'SD of paired outer-bootstrap bagged-peak differences; never divided by sqrt(nOuter)'; ...
    'MET null', 'H0: |Delta| <= delta0'; ...
    'MET alternative', 'H1: |Delta| > delta0'; ...
    'MET p-value', 'Wald normal approximation at the least-favorable interval-null boundary'; ...
    'Standardized effect', 'Bagged difference divided by bootstrap SE; standardized peak-shift statistic, not classical Cohen d'; ...
    'Bootstrap support', 'Descriptive only; not the reported frequentist p-value'; ...
    'Main figures', 'Outer-bootstrap distributions of bagged sLNv and 5th-sLNv peaks; display bin width is recorded in Configuration and is not raw temporal sampling resolution'; ...
    'Monte Carlo counts', 'nInner and nOuter control numerical precision, not biological sample size or power by themselves'};

% Overwrite the named analysis sheets on rerun; source workbooks are never
% modified. The smoke-test suffix prevents a test run from replacing the
% formal result workbook.
writetable(summaryTable, outputExcel, 'Sheet', 'Summary', ...
    'WriteMode', 'overwritesheet');
writetable(timepointDataTable, outputExcel, 'Sheet', 'TimepointData', ...
    'WriteMode', 'overwritesheet');
writetable(peakProbabilityTable, outputExcel, ...
    'Sheet', 'PeakProbabilities', 'WriteMode', 'overwritesheet');
writetable(hardBootstrapTable, outputExcel, 'Sheet', 'HardBootstrap', ...
    'WriteMode', 'overwritesheet');
writetable(outerBootstrapTable, outputExcel, 'Sheet', 'OuterBootstrap', ...
    'WriteMode', 'overwritesheet');
writetable(diagnosticTable, outputExcel, 'Sheet', 'Diagnostics', ...
    'WriteMode', 'overwritesheet');
writetable(configurationTable, outputExcel, 'Sheet', 'Configuration', ...
    'WriteMode', 'overwritesheet');
writecell(methodsNotes, outputExcel, 'Sheet', 'MethodsNotes', ...
    'WriteMode', 'overwritesheet');

batchElapsedSeconds = toc(batchStart);
fprintf('\n====================================================================\n');
fprintf('Completed %d files in %.1f seconds.\n', nFiles, batchElapsedSeconds);
fprintf('Combined Excel results: %s\n', outputExcel);
fprintf('Mean +/- SEM plotting data: %s\n', outputMeanSEMExcel);
fprintf('====================================================================\n');

%% ========================================================================
% Local helper functions
% ========================================================================

publicResult = struct('datasets', batchResults, 'summary', summaryTable, 'outerBootstrap', outerBootstrapTable, 'timepointData', timepointDataTable);
publicResult.configuration = cfg;
if exist('outputMAT', 'var') && isfile(outputMAT)
    save(outputMAT, 'cfg', 'publicResult', '-append');
else
    save(fullfile(cfg.outputDir, cfg.outputStem + "_results.mat"), 'cfg', 'publicResult');
end
end

function [time, Mdata, Edata, timepointTable, dataInfo] = ...
    readPairedWideWorkbook(filePath, M_sheet, E_sheet)
% Read the no-header, wide-format paired workbook and retain complete pairs.

    rawM = bootstrapReadMatrix(filePath,M_sheet);
    rawE = bootstrapReadMatrix(filePath,E_sheet);

    if size(rawM, 1) ~= size(rawE, 1) || ...
            size(rawM, 2) ~= size(rawE, 2)
        error('figs10:SheetSizeMismatch', ...
            'The two sheets in %s have different dimensions.', filePath);
    end
    if size(rawM, 2) < 2
        error('figs10:NoSampleColumns', ...
            '%s does not contain sample columns after the time column.', ...
            filePath);
    end

    timeM = rawM(:, 1);
    timeE = rawE(:, 1);
    if any(~isfinite(timeM)) || any(~isfinite(timeE))
        error('figs10:InvalidTimeColumn', [ ...
            '%s contains a nonnumeric/empty time value. These files are ' ...
            'expected to have no header row.'], filePath);
    end
    if any(timeM ~= timeE)
        error('figs10:TimeMismatch', ...
            'The two sheets in %s do not have identical time columns.', ...
            filePath);
    end
    if numel(unique(timeM)) ~= numel(timeM)
        error('figs10:DuplicateTime', ...
            '%s contains duplicated time rows.', filePath);
    end

    [time, sortIndex] = sort(timeM);
    rawM = rawM(sortIndex, :);
    rawE = rawE(sortIndex, :);

    nTime = numel(time);
    Mdata = cell(1, nTime);
    Edata = cell(1, nTime);
    nPairs = NaN(nTime, 1);
    mean_M = NaN(nTime, 1);
    sd_M = NaN(nTime, 1);
    mean_E = NaN(nTime, 1);
    sd_E = NaN(nTime, 1);
    pairedCorrelation = NaN(nTime, 1);

    for j = 1:nTime
        Mrow = rawM(j, 2:end);
        Erow = rawE(j, 2:end);
        finiteM = isfinite(Mrow);
        finiteE = isfinite(Erow);
        incompletePair = xor(finiteM, finiteE);
        if any(incompletePair)
            warning('figs10:IncompletePair', [ ...
                '%s, time %.4g: %d positions are present in only one ' ...
                'sheet and are excluded to preserve pairing.'], ...
                filePath, time(j), sum(incompletePair));
        end

        completePair = finiteM & finiteE;
        Mj = Mrow(completePair).';
        Ej = Erow(completePair).';
        if isempty(Mj)
            error('figs10:NoPairsAtTime', ...
                '%s has no complete M/E pairs at time %.4g.', ...
                filePath, time(j));
        end

        Mdata{j} = Mj;
        Edata{j} = Ej;
        nPairs(j) = numel(Mj);
        mean_M(j) = mean(Mj);
        sd_M(j) = std(Mj);
        mean_E(j) = mean(Ej);
        sd_E(j) = std(Ej);
        if numel(Mj) > 2 && std(Mj) > 0 && std(Ej) > 0
            corrMatrix = corrcoef(Mj, Ej);
            pairedCorrelation(j) = corrMatrix(1, 2);
        end
    end

    timepointTable = table(time(:), nPairs, mean_M, sd_M, mean_E, ...
        sd_E, pairedCorrelation, 'VariableNames', {'Time_h', 'N_pairs', ...
        'Mean_sLNv', 'SD_sLNv', 'Mean_5th_sLNv', 'SD_5th_sLNv', ...
        'PairedCorrelation'});
    dataInfo.nPairsMin = min(nPairs);
    dataInfo.nPairsMax = max(nPairs);
end

function edges = centeredHistogramEdges(values, binWidth)
% Construct edges whose bin centers are integer multiples of binWidth.
    values = values(isfinite(values));
    if isempty(values)
        edges = [-binWidth / 2, binWidth / 2];
        return;
    end
    lowestCenter = floor(min(values) / binWidth) * binWidth;
    highestCenter = ceil(max(values) / binWidth) * binWidth;
    edges = (lowestCenter - binWidth / 2) : binWidth : ...
        (highestCenter + binWidth / 2);
    if numel(edges) < 2
        edges = lowestCenter + [-binWidth / 2, binWidth / 2];
    end
end

function createMainOuterFigure(outerPeak_M, outerPeak_E, binWidth, ...
    xLimits, visibility, outputPath, datasetName)
% Original figs10.m style, now showing outer-bootstrap Tbag distributions.

    peakEdges = centeredHistogramEdges([outerPeak_M; outerPeak_E], binWidth);
    M_color = [255/255, 128/255, 0];
    E_color = [0, 0, 1];

    fig = figure('Name', char(datasetName + ...
        " outer-bootstrap bagged peaks"), 'NumberTitle', 'off', ...
        'Units', 'centimeters', 'Position', [2, 2, 6, 4], ...
        'Color', 'w', 'Visible', visibility);
    hold on;
    histogram(outerPeak_E, peakEdges, 'FaceColor', E_color, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8, ...
        'DisplayName', '5th sLNv', 'Normalization', 'probability');
    histogram(outerPeak_M, peakEdges, 'FaceColor', M_color, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8, ...
        'DisplayName', 'sLNv', 'Normalization', 'probability');
    hold off;

    ax = gca;
    ax.Units = 'centimeters';
    ax.Position(3:4) = [4.25, 2.5];
    ax.LineWidth = 1;
    ax.TickDir = 'out';
    ax.TickLength = [0.02, 0.02];
    ax.XTickLabelRotation = 0;
    ax.XRuler.TickLabelGapOffset = -1;
    ax.FontName = 'Arial';
    ax.FontSize = 7.5;
    ax.Color = 'none';
    box off;
    ylim([0, 1]);
    xlim(xLimits);
    xticks(xLimits(1):2:xLimits(2));
    yticks(0:0.5:1);

    set(findall(fig, '-property', 'FontUnits'), 'FontUnits', 'points');
    set(findall(fig, '-property', 'FontSize'), 'FontSize', 7.5);
    set(findall(fig, '-property', 'FontSizeMode'), ...
        'FontSizeMode', 'manual');
    set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');

    fig.PaperUnits = 'centimeters';
    fig.PaperSize = [6, 4];
    fig.PaperPosition = [0, 0, 6, 4];
    fig.PaperPositionMode = 'manual';
    fig.InvertHardcopy = 'off';
    drawnow;
    print(fig, char(outputPath), '-dsvg', '-vector');
    close(fig);
end

function createDiagnosticFigure(outerTbag_M, outerTbag_E, outerDiff, ...
    Tbag_M, Tbag_E, baggedDiff, baggedSE, delta0, JB_p, ...
    outerCorrelation, diffBinWidth, visibility, outputPath, datasetName)

    M_color = [255/255, 128/255, 0];
    E_color = [0, 0, 1];
    diagnosticFig = figure('Name', char(datasetName + ...
        " paired bagged-peak diagnostics"), 'NumberTitle', 'off', ...
        'Color', 'w', 'Units', 'centimeters', ...
        'Position', [3 3 18 14], 'Visible', visibility);
    layout = tiledlayout(diagnosticFig, 2, 2, ...
        'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile(layout, 1);
    diffEdges = centeredHistogramEdges(outerDiff, diffBinWidth);
    histogram(outerDiff, diffEdges, 'Normalization', 'pdf', ...
        'FaceColor', [0.45 0.45 0.45], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.75, 'DisplayName', 'Outer bootstrap');
    hold on;
    xMin = min(outerDiff);
    xMax = max(outerDiff);
    if xMin == xMax
        xMin = xMin - diffBinWidth;
        xMax = xMax + diffBinWidth;
    end
    xGrid = linspace(xMin, xMax, 300);
    if baggedSE > 0
        plot(xGrid, normpdf(xGrid, mean(outerDiff), baggedSE), ...
            'Color', M_color, 'LineWidth', 1.2, ...
            'DisplayName', 'Normal fit');
    end
    xline(-delta0, '--k', '-delta0', 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
    xline(delta0, '--k', '+delta0', 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
    xline(baggedDiff, '-', 'Observed diff', 'Color', [0.75 0 0], ...
        'LineWidth', 1.2, 'HandleVisibility', 'off');
    xlabel('Bagged peak difference M-E (h)');
    ylabel('Density');
    title('Outer-bootstrap distribution');
    legend('Location', 'best');
    box off;

    nexttile(layout, 2);
    if numel(outerDiff) >= 2
        qqplot(outerDiff);
    end
    xlabel('Standard-normal quantiles');
    ylabel('Outer-difference quantiles (h)');
    title(sprintf('Normal Q-Q plot (JB p=%.3g)', JB_p));
    box off;

    nexttile(layout, 3);
    scatter(outerTbag_M, outerTbag_E, 10, [0.25 0.25 0.25], ...
        'filled', 'MarkerFaceAlpha', 0.25);
    xline(Tbag_M, '--', 'sLNv estimate', 'Color', M_color, ...
        'LineWidth', 0.9);
    yline(Tbag_E, '--', '5th estimate', 'Color', E_color, ...
        'LineWidth', 0.9);
    xlabel('Outer Tbag sLNv (h)');
    ylabel('Outer Tbag 5th sLNv (h)');
    title(sprintf('Paired outer estimates (r=%.3f)', outerCorrelation));
    box off;

    nexttile(layout, 4);
    nValid = numel(outerDiff);
    runningMean = cumsum(outerDiff) ./ (1:nValid).';
    runningSD = NaN(nValid, 1);
    for k = 2:nValid
        runningSD(k) = std(outerDiff(1:k));
    end
    yyaxis left;
    plot(1:nValid, runningMean, 'Color', [0.2 0.2 0.2], ...
        'LineWidth', 1);
    yline(baggedDiff, '--', 'Final estimate', ...
        'Color', [0.2 0.2 0.2], 'HandleVisibility', 'off');
    ylabel('Cumulative mean difference (h)');
    yyaxis right;
    plot(1:nValid, runningSD, 'Color', [0.75 0 0], 'LineWidth', 1);
    yline(baggedSE, '--', 'Final SE', 'Color', [0.75 0 0], ...
        'HandleVisibility', 'off');
    ylabel('Cumulative bootstrap SE (h)');
    xlabel('Valid outer replicates');
    title('Outer Monte Carlo convergence');
    box off;

    sgtitle(layout, char(datasetName + " paired bagged-bootstrap diagnostics"));
    set(findall(diagnosticFig, '-property', 'FontName'), ...
        'FontName', 'Arial');
    drawnow;
    print(diagnosticFig, char(outputPath), '-dsvg', '-vector');
    close(diagnosticFig);
end
