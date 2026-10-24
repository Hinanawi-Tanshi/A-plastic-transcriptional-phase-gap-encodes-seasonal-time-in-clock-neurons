function publicResult = runFlyPairedIntensityAnalysis(cfg)
% Shared publication pipeline; settings are supplied by the figure entry script.
figureCleanup = bootstrapFigureContext(cfg); %#ok<NASGU>
bootstrapRawStore('configure',cfg.inputFile,cfg);
scriptFolder = cfg.outputDir;
nInner = cfg.nInner;
nOuter = cfg.nOuter;
nCTBootstrap = cfg.nCTBootstrap;
delta0 = cfg.delta0;
alpha = cfg.alpha;
rngSeed = cfg.rngSeed;
recordingWindowHours = cfg.recordingWindowHours;
representativeTimeOffset = cfg.representativeTimeOffset;
timestampToleranceHours = cfg.timestampToleranceHours;
outerDifferenceBinWidth = cfg.outerDifferenceBinWidth;
outerPeakBinWidth = cfg.outerPeakBinWidth;
drawPanelA = cfg.drawPanelA;
drawPanelC = cfg.drawPanelC;
drawSixPanel = cfg.drawSixPanel;
drawDiagnostics = cfg.drawDiagnostics;
panelCUseBoxplot = cfg.panelCUseBoxplot;
panelAPeakBinWidth = cfg.panelAPeakBinWidth;
panelCRequestedCT = cfg.panelCRequestedCT;
deltaAmplitude = cfg.deltaAmplitude;
M_color = cfg.M_color;
E_color = cfg.E_color;
diffColor = cfg.diffColor;
parentFolders = normalizeIntensitySourceFolders(fullfile(cfg.inputFile, cfg.batchNames));
parentFolder = parentFolders(1);
outputStem = char(cfg.outputStem);
outputExcel = fullfile(scriptFolder, [outputStem '_results.xlsx']);
outputMAT = fullfile(scriptFolder, [outputStem '_results.mat']);
outputMainSVG = fullfile(scriptFolder, [outputStem '.svg']);
outputDiagnosticSVG = fullfile(scriptFolder, [outputStem '_diagnostics.svg']);
outputPanelASVG = fullfile(scriptFolder, [outputStem '_panelA_bagged_peak.svg']);
outputPanelCSVG = fullfile(scriptFolder, [outputStem '_panelC_active_amplitude.svg']);

fprintf('Reading paired signal/intensity workbooks from %d source folders...\n', ...
    numel(parentFolders));
[cellMetricsTable, flyMetricsTable, importAuditTable, sourceSummaryTable] = ...
    loadPairedIntensityData(parentFolders, recordingWindowHours, ...
    timestampToleranceHours);
disp(sourceSummaryTable(:, {'SourceIndex', 'SourceName', ...
    'N_flyRecords', 'N_validPairedFlies', 'N_cellSheets', ...
    'N_ignoredXLSX', 'N_CTfoldersWithoutPairs'}));

validPairRows = flyMetricsTable.ValidPrimaryPair & ...
    isfinite(flyMetricsTable.AnalysisTime_midpoint_h);
if ~any(validPairRows)
    error('fig2cIntensity:NoValidPairs', ...
        'No complete fly-level M/E pairs were recovered.');
end

ctList = unique(flyMetricsTable.AnalysisTime_midpoint_h(validPairRows));
ctList = ctList(:).';
nCT = numel(ctList);
if nCT < 2
    error('fig2cIntensity:TooFewCT', ...
        'At least two valid midpoint CT groups are required.');
end

% Primary outcome: time-normalized gated output multiplied by 2 h. This is
% a 2-h-equivalent event-gated AUC. The normalization prevents shorter raw
% recordings from being interpreted as lower transcriptional output.
[outputM_byCT, outputE_byCT] = groupPairedMetricByCT( ...
    flyMetricsTable, ctList, 'GatedAUC_2hEquivalent_M', ...
    'GatedAUC_2hEquivalent_E', validPairRows);
[dutyM_byCT, dutyE_byCT] = groupPairedMetricByCT( ...
    flyMetricsTable, ctList, 'DutyCycle_M', 'DutyCycle_E', ...
    validPairRows);
[activeM_byCT, activeE_byCT] = groupPairedMetricByCT( ...
    flyMetricsTable, ctList, 'ActiveAmplitude_M', ...
    'ActiveAmplitude_E', validPairRows);
[contrastM_byCT, contrastE_byCT] = groupPairedMetricByCT( ...
    flyMetricsTable, ctList, 'ActiveContrast_M', ...
    'ActiveContrast_E', validPairRows);
[durationM_byCT, durationE_byCT] = groupPairedMetricByCT( ...
    flyMetricsTable, ctList, 'MeanObservedDuration_h_M', ...
    'MeanObservedDuration_h_E', validPairRows);

ctSummaryTable = buildCTSummaryTable(flyMetricsTable, ctList, ...
    validPairRows, recordingWindowHours);

fprintf('\n========== fig 2c paired intensity summary ==========\n');
fprintf('Matched fly pairs: %d; cell sheets: %d; midpoint CT groups: %d\n', ...
    sum(validPairRows), height(cellMetricsTable), nCT);
fprintf(['Primary outcome: event-gated, time-normalized, 2-h-equivalent ' ...
    'AUC (DeltaF/F*h)\n']);
fprintf('%-8s %-7s %-13s %-13s %-11s %-11s\n', ...
    'Mid-CT', 'nFly', 'M gated AUC', 'E gated AUC', 'M duty', 'E duty');
for j = 1:nCT
    fprintf('%-8.4g %-7d %-13.5f %-13.5f %-11.5f %-11.5f\n', ...
        ctList(j), numel(outputM_byCT{j}), ...
        mean(outputM_byCT{j}, 'omitnan'), ...
        mean(outputE_byCT{j}, 'omitnan'), ...
        mean(dutyM_byCT{j}, 'omitnan'), ...
        mean(dutyE_byCT{j}, 'omitnan'));
    if numel(outputM_byCT{j}) < 2
        warning('fig2cIntensity:SingleFlyCT', ...
            ['Midpoint CT %.4g contains one paired fly. Bootstrap ' ...
             'resampling cannot estimate within-CT biological variation.'], ...
            ctList(j));
    end
end

%% ========================================================================
% PART 2. Paired bagging, outer bootstrap, and Wald-type MET
% ========================================================================
% The primary peak estimator is applied to the fly-level 2-h-equivalent
% gated-output curve. Multiplying the time-normalized output by 2 h changes
% only its units, not its peak time.

[hardPeak_M, hardPeak_E] = computeHardPeaks( ...
    outputM_byCT, outputE_byCT, ctList);
hardDiff = hardPeak_M - hardPeak_E;

rng(rngSeed + 1, 'twister');
fprintf('\nComputing observed paired bagged peaks (nInner=%d)...\n', ...
    nInner);
[Tbag_M, Tbag_E, innerPeakTimes_M, innerPeakTimes_E, P_M, P_E] = ...
    computePairedBaggedPeaks(outputM_byCT, outputE_byCT, ctList, nInner, "finite");
baggedDiff = Tbag_M - Tbag_E;

outerHardPeak_M = NaN(nOuter, 1);
outerHardPeak_E = NaN(nOuter, 1);
outerTbag_M = NaN(nOuter, 1);
outerTbag_E = NaN(nOuter, 1);
rng(rngSeed + 2, 'twister');
outerStart = tic;
progressStep = max(1, floor(nOuter / 20));
fprintf('Running paired outer nonparametric bootstrap (nOuter=%d)...\n', ...
    nOuter);
for b = 1:nOuter
    [outerData_M, outerData_E] = ...
        bootstrapPairedDataset(outputM_byCT, outputE_byCT);
    [outerHardPeak_M(b), outerHardPeak_E(b)] = ...
        computeHardPeaks(outerData_M, outerData_E, ctList);
    [outerTbag_M(b), outerTbag_E(b)] = ...
        computePairedBaggedPeaks(outerData_M, outerData_E, ...
        ctList, nInner, "finite");
    if mod(b, progressStep) == 0 || b == nOuter
        fprintf('  outer bootstrap: %d/%d complete\n', b, nOuter);
    end
end
outerElapsedSeconds = toc(outerStart);

outerHardDiff = outerHardPeak_M - outerHardPeak_E;
outerBaggedDiff = outerTbag_M - outerTbag_E;
validOuterHard = isfinite(outerHardPeak_M) & ...
    isfinite(outerHardPeak_E) & isfinite(outerHardDiff);
validOuterBagged = isfinite(outerTbag_M) & ...
    isfinite(outerTbag_E) & isfinite(outerBaggedDiff);
if sum(validOuterBagged) < 2
    error('fig2cIntensity:TooFewOuter', ...
        'Fewer than two valid outer bootstrap replicates were obtained.');
end

% Bootstrap SE is the SD of bootstrap estimates. Never divide by
% sqrt(nOuter): nOuter controls Monte Carlo precision, not biological n.
hardSE = std(outerHardDiff(validOuterHard), 'omitnan');
hardCI95 = prctile(outerHardDiff(validOuterHard), [2.5 97.5]);
hardSD_M = std(outerHardPeak_M(validOuterHard), 'omitnan');
hardSD_E = std(outerHardPeak_E(validOuterHard), 'omitnan');
hardCI95_M = prctile(outerHardPeak_M(validOuterHard), [2.5 97.5]);
hardCI95_E = prctile(outerHardPeak_E(validOuterHard), [2.5 97.5]);
hardBootstrapSupport = mean(abs(outerHardDiff(validOuterHard)) > delta0);
[hardMET_z, hardMET_p, hardMET_z_near, hardMET_z_far, ...
    hardBoundaryZ_lower, hardBoundaryZ_upper] = ...
    computeWaldMET(hardDiff, hardSE, delta0);

baggedSE = std(outerBaggedDiff(validOuterBagged), 'omitnan');
outerSD_M = std(outerTbag_M(validOuterBagged), 'omitnan');
outerSD_E = std(outerTbag_E(validOuterBagged), 'omitnan');
baggedCI95 = prctile(outerBaggedDiff(validOuterBagged), [2.5 97.5]);
baggedCI95_M = prctile(outerTbag_M(validOuterBagged), [2.5 97.5]);
baggedCI95_E = prctile(outerTbag_E(validOuterBagged), [2.5 97.5]);
[MET_z, MET_p, MET_z_near, MET_z_far, ...
    boundaryZ_lower, boundaryZ_upper] = ...
    computeWaldMET(baggedDiff, baggedSE, delta0);
baggedBootstrapSupport = ...
    mean(abs(outerBaggedDiff(validOuterBagged)) > delta0);
thresholdExcess = abs(baggedDiff) - delta0;
standardizedEffect = baggedDiff / baggedSE; % not classical Cohen's d
standardizedThresholdExcess = thresholdExcess / baggedSE;

[peakStability_M, peakEntropy_M] = computePeakStability(P_M);
[peakStability_E, peakEntropy_E] = computePeakStability(P_E);
peakConfidence_M = max(P_M);
peakConfidence_E = max(P_E);

% These divisions by sqrt(number of Monte Carlo replicates) quantify only
% simulation noise, not experimental sampling uncertainty or statistical
% power. They are never used as biological SEs or in the MET.
innerMonteCarloSE_M = std(innerPeakTimes_M, 'omitnan') / ...
    sqrt(sum(isfinite(innerPeakTimes_M)));
innerMonteCarloSE_E = std(innerPeakTimes_E, 'omitnan') / ...
    sqrt(sum(isfinite(innerPeakTimes_E)));
outerSEMonteCarloErrorApprox = baggedSE / ...
    sqrt(2 * max(sum(validOuterBagged) - 1, 1));

outerDiffValid = outerBaggedDiff(validOuterBagged);
outerCovarianceMatrix = cov(outerTbag_M(validOuterBagged), ...
    outerTbag_E(validOuterBagged));
outerCovariance_M_E = outerCovarianceMatrix(1, 2);
outerCorrelationMatrix = corrcoef(outerTbag_M(validOuterBagged), ...
    outerTbag_E(validOuterBagged));
outerCorrelation_M_E = outerCorrelationMatrix(1, 2);
outerSkewness = skewness(outerDiffValid, 0);
outerKurtosis = kurtosis(outerDiffValid, 0);
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
    outerBimodalityCoefficient > 0.65 || outerKurtosis > 5;
if stronglyNonNormal
    warning('fig2cIntensity:OuterNonNormal', ...
        ['Outer bootstrap distribution is strongly non-normal; consider ' ...
         'a studentized-bootstrap / bootstrap-t MET as a sensitivity ' ...
         'analysis. The Wald-type MET remains unchanged.']);
end

if numel(ctList) >= 2 && ...
        ((P_M(1) > 0.05 && P_M(end) > 0.05) || ...
         (P_E(1) > 0.05 && P_E(end) > 0.05))
    warning('fig2cIntensity:CircularBoundaryMass', ...
        ['Peak probability occupies both CT1 and CT25. Inspect a circular ' ...
         'phase sensitivity analysis; no silent method change was made.']);
end

statsHard = struct();
statsHard.peak_M = hardPeak_M;
statsHard.peak_E = hardPeak_E;
statsHard.diff = hardDiff;
statsHard.SE = hardSE;
statsHard.CI95 = hardCI95;
statsHard.SD_M = hardSD_M;
statsHard.SD_E = hardSD_E;
statsHard.CI95_M = hardCI95_M;
statsHard.CI95_E = hardCI95_E;
statsHard.bootstrapSupport = hardBootstrapSupport;
statsHard.MET_z = hardMET_z;
statsHard.MET_p = hardMET_p;
statsHard.MET_z_near = hardMET_z_near;
statsHard.MET_z_far = hardMET_z_far;
statsHard.boundaryZ_lower = hardBoundaryZ_lower;
statsHard.boundaryZ_upper = hardBoundaryZ_upper;

statsBagged = struct();
statsBagged.peak_M = Tbag_M;
statsBagged.peak_E = Tbag_E;
statsBagged.diff = baggedDiff;
statsBagged.SE = baggedSE;
statsBagged.CI95 = baggedCI95;
statsBagged.CI_low = baggedCI95(1);
statsBagged.CI_high = baggedCI95(2);
statsBagged.SD_M = outerSD_M;
statsBagged.SD_E = outerSD_E;
statsBagged.CI95_M = baggedCI95_M;
statsBagged.CI95_E = baggedCI95_E;
statsBagged.MET_z = MET_z;
statsBagged.MET_p = MET_p;
statsBagged.MET_z_near = MET_z_near;
statsBagged.MET_z_far = MET_z_far;
statsBagged.boundaryZ_lower = boundaryZ_lower;
statsBagged.boundaryZ_upper = boundaryZ_upper;
statsBagged.t_lower_like = boundaryZ_lower;
statsBagged.t_upper_like = boundaryZ_upper;
statsBagged.standardizedEffect = standardizedEffect;
statsBagged.thresholdExcess = thresholdExcess;
statsBagged.standardizedThresholdExcess = ...
    standardizedThresholdExcess;
statsBagged.bootstrapSupport = baggedBootstrapSupport;
statsBagged.peakProbability_M = P_M;
statsBagged.peakProbability_E = P_E;
statsBagged.peakConfidence_M = peakConfidence_M;
statsBagged.peakConfidence_E = peakConfidence_E;
statsBagged.peakEntropy_M = peakEntropy_M;
statsBagged.peakEntropy_E = peakEntropy_E;
statsBagged.peakStability_M = peakStability_M;
statsBagged.peakStability_E = peakStability_E;
statsBagged.nInner = nInner;
statsBagged.nOuter = nOuter;
statsBagged.delta0 = delta0;
statsBagged.alpha = alpha;
statsBagged.pairedSamplingUnit = 'fly';
statsBagged.primaryOutcome = ...
    '2-h-equivalent event-gated AUC after observed-duration normalization';

%% ========================================================================
% Per-CT paired differences and tabular output
% ========================================================================
rng(rngSeed + 3, 'twister');
ctDifferenceTable = computeCTDifferenceTable(ctList, outputM_byCT, ...
    outputE_byCT, nCTBootstrap);

% Standalone Panel A tables. CT-level n/SD/SEM use independent paired
% flies. Outer SD is the bootstrap sampling SE of the peak estimator.
% Outer SD/sqrt(Nouter) is reported only as MonteCarloSEM and is never used
% for experimental inference or as the MET denominator.
panelACTStats = table(ctList(:), ctSummaryTable.N_pairedFlies, ...
    ctSummaryTable.Mean_GatedAUC_2hEquivalent_M, ...
    ctSummaryTable.SD_GatedAUC_2hEquivalent_M, ...
    ctSummaryTable.Nvalid_GatedAUC_2hEquivalent_M, ...
    ctSummaryTable.SEM_GatedAUC_2hEquivalent_M, ...
    ctSummaryTable.Mean_GatedAUC_2hEquivalent_E, ...
    ctSummaryTable.SD_GatedAUC_2hEquivalent_E, ...
    ctSummaryTable.Nvalid_GatedAUC_2hEquivalent_E, ...
    ctSummaryTable.SEM_GatedAUC_2hEquivalent_E, ...
    ctDifferenceTable.MeanDifference_MminusE, ...
    ctDifferenceTable.SD_pairedDifference, ...
    ctDifferenceTable.SEM_pairedDifference, ...
    ctDifferenceTable.BootstrapCI95_low, ...
    ctDifferenceTable.BootstrapCI95_high, ...
    'VariableNames', {'CT_midpoint_h', 'N_pairedFlies', ...
    'M_mean_DeltaFoverF_h', 'M_SD_DeltaFoverF_h', 'M_n', ...
    'M_SEM_DeltaFoverF_h', 'E_mean_DeltaFoverF_h', ...
    'E_SD_DeltaFoverF_h', 'E_n', 'E_SEM_DeltaFoverF_h', ...
    'MeanPairedDifference_MminusE', 'SD_pairedDifference', ...
    'SEM_pairedDifference', 'PairedBootstrapCI95_low', ...
    'PairedBootstrapCI95_high'});

validOuterM = isfinite(outerTbag_M);
validOuterE = isfinite(outerTbag_E);
nValidOuterM = sum(validOuterM);
nValidOuterE = sum(validOuterE);
panelABootstrapStats = table(["M"; "E"], [Tbag_M; Tbag_E], ...
    [mean(outerTbag_M(validOuterM), 'omitnan'); ...
     mean(outerTbag_E(validOuterE), 'omitnan')], ...
    [outerSD_M; outerSD_E], [nValidOuterM; nValidOuterE], ...
    [outerSD_M / sqrt(nValidOuterM); ...
     outerSD_E / sqrt(nValidOuterE)], ...
    [baggedCI95_M(1); baggedCI95_E(1)], ...
    [baggedCI95_M(2); baggedCI95_E(2)], ...
    [sum(validPairRows); sum(validPairRows)], ...
    [peakConfidence_M; peakConfidence_E], ...
    [peakStability_M; peakStability_E], ...
    'VariableNames', {'Condition', 'ObservedTbag_h', ...
    'OuterMeanTbag_h', 'OuterSD_asBootstrapSE_h', 'N_validOuter', ...
    'OuterMonteCarloSEM_h_notInference', 'OuterCI95_low_h', ...
    'OuterCI95_high_h', 'N_biologicalFlyPairs_total', ...
    'PeakConfidence', 'PeakStability'});

panelAMETTable = table(baggedDiff, baggedSE, baggedCI95(1), ...
    baggedCI95(2), delta0, alpha, thresholdExcess, ...
    standardizedEffect, standardizedThresholdExcess, ...
    boundaryZ_lower, boundaryZ_upper, MET_z_near, MET_z_far, ...
    MET_z, MET_p, baggedBootstrapSupport, sum(validOuterBagged), ...
    'VariableNames', {'Difference_MminusE_h', ...
    'OuterBootstrapSE_h', 'CI95_low_h', 'CI95_high_h', ...
    'MinimumEffectThreshold_delta0_h', 'Alpha', ...
    'ThresholdExcess_h', 'StandardizedPeakShift_notCohenD', ...
    'StandardizedThresholdExcess', 'BoundaryZ_lower', ...
    'BoundaryZ_upper', 'MET_Z_near', 'MET_Z_far', 'MET_Z', ...
    'MET_p', 'BootstrapSupport_descriptive_notPvalue', ...
    'N_validOuter'});

% Standalone Panel C uses the prespecified common active window. The main
% scientific question is whether active-state amplitude changes across CT
% within M and within E. M-versus-E paired tests remain secondary outputs.
panelCIndex = find(ismember(ctList, panelCRequestedCT));
panelCCT = ctList(panelCIndex);
missingPanelCCT = setdiff(panelCRequestedCT, panelCCT);
if ~isempty(missingPanelCCT)
    warning('fig2cIntensity:MissingPanelCCT', ...
        'Requested Panel C CT groups are absent: %s', ...
        strjoin(string(missingPanelCCT), ', '));
end
if numel(panelCCT) < 2
    error('fig2cIntensity:TooFewPanelCCT', ...
        'Panel C requires at least two requested CT groups.');
end
panelCActiveM_byCT = activeM_byCT(panelCIndex);
panelCActiveE_byCT = activeE_byCT(panelCIndex);
panelCContrastM_byCT = contrastM_byCT(panelCIndex);
panelCContrastE_byCT = contrastE_byCT(panelCIndex);

rng(rngSeed + 4, 'twister');
[panelCCTStats, panelCPairedTests] = computePanelCStatistics( ...
    panelCCT, panelCActiveM_byCT, panelCActiveE_byCT, ...
    nCTBootstrap, deltaAmplitude, alpha);
rng(rngSeed + 5, 'twister');
[panelCContrastStats, panelCContrastTests] = ...
    computePanelCStatistics(panelCCT, panelCContrastM_byCT, ...
    panelCContrastE_byCT, nCTBootstrap, deltaAmplitude, alpha);
[panelCLMEAnova, panelCLMECoefficients, panelCLMEInfo] = ...
    fitPanelCAmplitudeLME(flyMetricsTable, validPairRows, panelCCT);
rng(rngSeed + 6, 'twister');
panelCModulationTable = computePanelCModulation( ...
    panelCCT, panelCActiveM_byCT, panelCActiveE_byCT, ...
    nCTBootstrap, deltaAmplitude);
[panelCWithinCTOmnibus, panelCWithinCTPairwiseTOST] = ...
    computeWithinConditionCTTests(panelCCT, panelCActiveM_byCT, ...
    panelCActiveE_byCT, deltaAmplitude, alpha);

conditionTable = table(["M"; "E"], ...
    [hardPeak_M; hardPeak_E], [hardSD_M; hardSD_E], ...
    [hardCI95_M(1); hardCI95_E(1)], ...
    [hardCI95_M(2); hardCI95_E(2)], ...
    [Tbag_M; Tbag_E], [outerSD_M; outerSD_E], ...
    [baggedCI95_M(1); baggedCI95_E(1)], ...
    [baggedCI95_M(2); baggedCI95_E(2)], ...
    [peakConfidence_M; peakConfidence_E], ...
    [peakStability_M; peakStability_E], ...
    'VariableNames', {'Condition', 'HardPeak_midpointCT_h', ...
    'HardPeak_outerSD_h', 'HardPeak_CI95_low_h', ...
    'HardPeak_CI95_high_h', 'BaggedPeak_midpointCT_h', ...
    'BaggedPeak_outerSD_h', 'BaggedPeak_CI95_low_h', ...
    'BaggedPeak_CI95_high_h', 'PeakConfidence', 'PeakStability'});

summaryTable = table(hardDiff, hardSE, hardCI95(1), hardCI95(2), ...
    hardMET_z, hardMET_p, hardBootstrapSupport, ...
    baggedDiff, baggedSE, baggedCI95(1), baggedCI95(2), ...
    thresholdExcess, standardizedEffect, standardizedThresholdExcess, ...
    boundaryZ_lower, boundaryZ_upper, MET_z, MET_p, ...
    baggedBootstrapSupport, delta0, alpha, ...
    'VariableNames', {'HardDifference_MminusE_h', 'Hard_SE_h', ...
    'Hard_CI95_low_h', 'Hard_CI95_high_h', 'Hard_MET_Z', ...
    'Hard_MET_p', 'HardBootstrapSupport_descriptive', ...
    'BaggedDifference_MminusE_h', 'Bagged_SE_h', ...
    'Bagged_CI95_low_h', 'Bagged_CI95_high_h', ...
    'ThresholdExcess_h', 'StandardizedPeakShift_notCohenD', ...
    'StandardizedThresholdExcess', 'BoundaryZ_lower', ...
    'BoundaryZ_upper', 'MET_Z', 'MET_p', ...
    'BaggedBootstrapSupport_descriptive', 'Delta0_h', 'Alpha'});

peakProbabilityTable = table(ctList(:), P_M(:), P_E(:), ...
    'VariableNames', {'AnalysisTime_midpoint_h', ...
    'InnerPeakProbability_M', 'InnerPeakProbability_E'});

outerBootstrapTable = table((1:nOuter).', outerHardPeak_M, ...
    outerHardPeak_E, outerHardDiff, outerTbag_M, outerTbag_E, ...
    outerBaggedDiff, 'VariableNames', {'OuterReplicate', ...
    'OuterHardPeak_M_h', 'OuterHardPeak_E_h', ...
    'OuterHardDifference_MminusE_h', 'OuterTbag_M_h', ...
    'OuterTbag_E_h', 'OuterBaggedDifference_MminusE_h'});

diagnosticTable = table(outerSkewness, outerKurtosis, ...
    outerBimodalityCoefficient, JB_p, double(JB_reject), ...
    double(stronglyNonNormal), outerCovariance_M_E, ...
    outerCorrelation_M_E, innerMonteCarloSE_M, ...
    innerMonteCarloSE_E, outerSEMonteCarloErrorApprox, ...
    sum(validOuterBagged), outerElapsedSeconds, ...
    'VariableNames', {'OuterSkewness', 'OuterKurtosis_normal3', ...
    'OuterBimodalityCoefficient', 'JarqueBera_p', ...
    'JarqueBeraReject', 'StrongNonNormalityFlag', ...
    'OuterCovariance_M_E', 'OuterCorrelation_M_E', ...
    'InnerMonteCarloSE_M_h', 'InnerMonteCarloSE_E_h', ...
    'ApproxMonteCarloErrorOfBaggedSE_h', 'N_validOuter', ...
    'OuterElapsedSeconds'});

configurationTable = table(strjoin(parentFolders, newline), rngSeed, nInner, nOuter, ...
    nCTBootstrap, delta0, alpha, recordingWindowHours, ...
    representativeTimeOffset, timestampToleranceHours, ...
    outerDifferenceBinWidth, outerPeakBinWidth, panelAPeakBinWidth, ...
    string(strjoin(string(panelCRequestedCT), ',')), deltaAmplitude, ...
    drawPanelA, drawPanelC, drawSixPanel, drawDiagnostics, true, "fly", ...
    "equal-weight mean of four M-cell time-normalized outputs", ...
    "equal-weight mean of 5th/E-cell time-normalized outputs", ...
    'VariableNames', {'SourceFolder', 'RNG_seed', 'N_inner', 'N_outer', ...
    'N_CT_difference_bootstrap', 'Delta0_h', 'Alpha', ...
    'RecordingWindow_h', 'RepresentativeTimeOffset_h', ...
    'TimestampTolerance_h', 'OuterDifferenceBinWidth_h', ...
    'OuterPeakBinWidth_h', 'PanelA_peakBinWidth_h', ...
    'PanelC_requestedCT', 'PanelC_deltaAmplitude_DeltaFoverF', ...
    'DrawPanelA', 'DrawPanelC', 'DrawSixPanel', 'DrawDiagnostics', ...
    'PairedResampling', ...
    'BiologicalSamplingUnit', 'M_FlySummary', 'E_FlySummary'});

methodsNotes = { ...
    'Item', 'Definition'; ...
    'Raw inputs', 'Date/FLY tab-name pairing between binary and continuous-intensity workbooks; nominal CT inferred from the first valid timestamp'; ...
    'Multiple sources', 'Independent fly records are pooled by CT across parentFolders; SourceIndex, SourceFolder and OriginalFlyID retain provenance'; ...
    'Input coverage', 'InputSources reports per-source fly/cell counts and unrecognized XLSX counts; CT folders without recognized pairs issue a warning'; ...
    'Duplicate protection', 'Repeated canonical source paths or repeated recording-ID+CT recording keys cause an error before inference'; ...
    'Sheet pairing', 'Exact sheet-name match within each paired workbook'; ...
    'Sheet classification', 'Sheet containing 5th = E; every other sheet = M'; ...
    'Time window', 'Subtract first signal timestamp and retain relative times from 0 through 2 h'; ...
    'Analysis CT', 'Folder CTxx records CTxx to CTxx+2 and is represented by midpoint CTxx+1'; ...
    'Integration', 'Irregular-sampling trapezoidal weights computed from aligned timestamps'; ...
    'Duty cycle', 'Event-weighted duration divided by observed duration'; ...
    'Active amplitude', 'Event-gated intensity integral divided by event-weighted duration; undefined when no event occurs'; ...
    'Observed gated AUC', 'Integral of event(t) times DeltaF/F(t) over the observed window'; ...
    'Primary gated output', 'Observed gated AUC divided by observed duration, then multiplied by 2 h; a 2-h-equivalent AUC'; ...
    'Reason for duration normalization', 'Raw recordings differ in observed duration, so unnormalized AUC would confound output with recording length'; ...
    'Exact decomposition', 'Time-normalized gated output = duty cycle times active-state amplitude'; ...
    'M fly value', 'Equal-weight mean of valid non-5th cell outputs; conditional amplitude is reconstructed to preserve output = duty times amplitude'; ...
    'E fly value', 'Equal-weight mean of valid 5th-cell outputs; fly remains the biological unit'; ...
    'Hard peak', 'Sampled midpoint CT attaining the maximum CT-level mean primary-output curve'; ...
    'Bagged peak', 'Mean of inner-bootstrap discrete hard peak times'; ...
    'Outer bootstrap', 'Paired fly resampling within each CT followed by complete re-bagging'; ...
    'Bootstrap SE', 'SD of outer estimates; never divided by sqrt(nOuter)'; ...
    'MET null', 'H0: |Delta| <= delta0'; ...
    'MET alternative', 'H1: |Delta| > delta0'; ...
    'MET p-value', 'Wald normal approximation at the least-favorable/nearest interval-null boundary'; ...
    'Boundary Z', 'Standardized distances to null boundaries, not Student t statistics'; ...
    'Standardized effect', 'Peak difference divided by outer-bootstrap SE; not classical Cohen d'; ...
    'Bootstrap support', 'Descriptive only; not the reported frequentist p-value'; ...
    'Monte Carlo counts', 'nInner/nOuter improve numerical precision, not biological n or power'; ...
    'Panel A distribution', 'Outer-bootstrap bagged peak estimates; E is drawn before M with 1-h display bins'; ...
    'Panel A Monte Carlo SEM', 'Outer SD divided by sqrt(valid outer replicates describes simulation precision only and is never an experimental SE or MET denominator'; ...
    'Panel C window', 'Prespecified common active window CT11,13,15,17,19; not selected from statistical significance'; ...
    'Panel C biological unit', 'One complete paired M/E value per fly within CT'; ...
    'Panel C main CI', 'Percentile 95% CI from paired fly bootstrap within each CT'; ...
    'Panel C primary question', 'Whether active-state amplitude changes across CT11-19 within M and within E; this addresses temporal amplitude modulation'; ...
    'Panel C omnibus tests', 'Separate one-way categorical-CT ANOVA for M and E, with Kruskal-Wallis sensitivity and Brown-Forsythe variance check'; ...
    'Panel C no-change interpretation', 'ANOVA or Kruskal-Wallis p>0.05 means no detected CT effect and does not by itself prove amplitude stability'; ...
    'Panel C within-condition equivalence', 'All pairwise CT contrasts are evaluated by independent-sample Welch TOST using the prespecified deltaAmplitude margin and Holm adjustment within condition'; ...
    'Panel C global all-pairs equivalence', 'Intersection-union rule: all ten component TOST tests must pass at alpha; Holm-adjusted results are additionally provided for simultaneous pair-specific claims'; ...
    'Panel C pairwise sampling', 'Different CT groups contain different flies, so CT-to-CT comparisons are independent-sample rather than paired tests'; ...
    'Panel C equivalence margin', 'deltaAmplitude in DeltaF/F; user-adjustable and requires independent biological or technical justification'; ...
    'Panel C M-E TOST secondary', 'Paired M-E equivalence at a fixed CT is retained as a secondary analysis and is not the primary test of temporal amplitude modulation'; ...
    'Panel C TOST CI', 'Two-sided 90% t CI corresponding to alpha=0.05 TOST; distinct from the peak-timing MET interval null'; ...
    'Panel C multiplicity', 'Holm adjustment is reported across the five CT-specific tests for paired t, sign-rank, and TOST p-values'; ...
    'Panel C mixed model', 'Amplitude ~ Condition * categorical(CT) + (1|FlyID); p>0.05 alone is not evidence of equivalence'; ...
    'Panel C contrast sensitivity', 'Active amplitude minus inactive amplitude is summarized separately without replacing the prespecified raw active-amplitude outcome'; ...
    'Panel C range diagnostic', 'Bootstrap CI for max-minus-min CT mean and modulation index is descriptive; upper CI below deltaAmplitude is a bounded-modulation diagnostic, not a standalone equivalence test'};

% Remove only this script's own prior workbook so obsolete sheets cannot
% persist between runs. No raw or legacy analysis file is touched.
if isfile(outputExcel)
    delete(outputExcel);
end
writetable(summaryTable, outputExcel, 'Sheet', 'Summary');
writetable(conditionTable, outputExcel, 'Sheet', 'PeakConditions');
writetable(ctSummaryTable, outputExcel, 'Sheet', 'CTCurves');
writetable(ctDifferenceTable, outputExcel, 'Sheet', 'CTPairedDifferences');
writetable(flyMetricsTable, outputExcel, 'Sheet', 'FlyLevelData');
writetable(cellMetricsTable, outputExcel, 'Sheet', 'CellLevelData');
writetable(importAuditTable, outputExcel, 'Sheet', 'ImportAudit');
writetable(sourceSummaryTable, outputExcel, 'Sheet', 'InputSources');
writetable(peakProbabilityTable, outputExcel, 'Sheet', 'PeakProbabilities');
writetable(outerBootstrapTable, outputExcel, 'Sheet', 'OuterBootstrap');
writetable(diagnosticTable, outputExcel, 'Sheet', 'Diagnostics');
writetable(configurationTable, outputExcel, 'Sheet', 'Configuration');
writecell(methodsNotes, outputExcel, 'Sheet', 'MethodsNotes');
writetable(panelACTStats, outputExcel, 'Sheet', 'PanelA_CTStats');
writetable(panelABootstrapStats, outputExcel, 'Sheet', ...
    'PanelA_BootstrapStats');
writetable(panelAMETTable, outputExcel, 'Sheet', 'PanelA_MET');
writetable(panelCCTStats, outputExcel, 'Sheet', 'PanelC_CTStats');
writetable(panelCPairedTests, outputExcel, 'Sheet', ...
    'PanelC_PairedTests');
writetable(panelCContrastStats, outputExcel, 'Sheet', ...
    'PanelC_ContrastStats');
writetable(panelCContrastTests, outputExcel, 'Sheet', ...
    'PanelC_ContrastTests');
writetable(panelCLMEAnova, outputExcel, 'Sheet', 'PanelC_LME_ANOVA');
writetable(panelCLMECoefficients, outputExcel, 'Sheet', ...
    'PanelC_LME_Coeff');
writetable(panelCLMEInfo, outputExcel, 'Sheet', 'PanelC_LME_Info');
writetable(panelCModulationTable, outputExcel, 'Sheet', ...
    'PanelC_Modulation');
writetable(panelCWithinCTOmnibus, outputExcel, 'Sheet', ...
    'PanelC_WithinCT_Omnibus');
writetable(panelCWithinCTPairwiseTOST, outputExcel, 'Sheet', ...
    'PanelC_WithinCT_PairwiseTOST');

save(outputMAT, 'statsHard', 'statsBagged', 'cellMetricsTable', ...
    'flyMetricsTable', 'importAuditTable', 'sourceSummaryTable', ...
    'parentFolders', 'ctSummaryTable', ...
    'ctDifferenceTable', 'peakProbabilityTable', ...
    'outerBootstrapTable', 'diagnosticTable', 'configurationTable', ...
    'panelACTStats', 'panelABootstrapStats', 'panelAMETTable', ...
    'panelCCTStats', 'panelCPairedTests', 'panelCContrastStats', ...
    'panelCContrastTests', 'panelCLMEAnova', ...
    'panelCLMECoefficients', 'panelCLMEInfo', ...
    'panelCModulationTable', 'panelCWithinCTOmnibus', ...
    'panelCWithinCTPairwiseTOST', ...
    'innerPeakTimes_M', 'innerPeakTimes_E', 'outerHardPeak_M', ...
    'outerHardPeak_E', 'outerTbag_M', 'outerTbag_E', ...
    'outerBaggedDiff', 'ctList', 'panelCCT', 'deltaAmplitude');

fprintf('\n====================================================================\n');
fprintf('Intensity peak timing comparison (paired fly bootstrap)\n');
fprintf('Primary curve: 2-h-equivalent event-gated AUC\n');
fprintf('Minimum-effect threshold: delta0 = %.3f h\n', delta0);
fprintf('H0: |Delta| <= %.3f h\n', delta0);
fprintf('H1: |Delta| >  %.3f h\n', delta0);
fprintf('====================================================================\n');
fprintf('%-34s %15s %15s\n', 'Statistic', 'Hard argmax', 'Bagged peak');
fprintf('%-34s %15.3f %15.3f\n', 'Peak M (h)', hardPeak_M, Tbag_M);
fprintf('%-34s %15.3f %15.3f\n', 'Peak E (h)', hardPeak_E, Tbag_E);
fprintf('%-34s %15.3f %15.3f\n', 'Peak difference M-E (h)', ...
    hardDiff, baggedDiff);
fprintf('%-34s %15.3f %15.3f\n', 'Bootstrap SE (h)', hardSE, baggedSE);
fprintf('%-34s [%6.3f,%6.3f] [%6.3f,%6.3f]\n', ...
    '95% percentile CI (h)', hardCI95(1), hardCI95(2), ...
    baggedCI95(1), baggedCI95(2));
fprintf('%-34s %15.3f %15.3f\n', 'MET Z', hardMET_z, MET_z);
fprintf('%-34s %15.5g %15.5g\n', 'MET p-value', hardMET_p, MET_p);
fprintf('%-34s %15.3f %15.3f\n', ...
    'Bootstrap support (descriptive)', hardBootstrapSupport, ...
    baggedBootstrapSupport);
fprintf('%-34s %15s %15.3f\n', 'Peak confidence M', 'N/A', ...
    peakConfidence_M);
fprintf('%-34s %15s %15.3f\n', 'Peak confidence E', 'N/A', ...
    peakConfidence_E);
fprintf('%-34s %15s %15.3f\n', 'Peak stability M', 'N/A', ...
    peakStability_M);
fprintf('%-34s %15s %15.3f\n', 'Peak stability E', 'N/A', ...
    peakStability_E);
fprintf('====================================================================\n');
fprintf('Excel results: %s\n', outputExcel);
fprintf('MAT results:   %s\n', outputMAT);
fprintf(['Bootstrap support is descriptive only; the reported frequentist ' ...
    'p-value is MET_p.\n']);

fprintf('\nPanel C: within-condition active-amplitude modulation, CT11-19\n');
fprintf(['Omnibus p>0.05 means no detected CT effect; evidence of stability ' ...
    'requires the prespecified equivalence analysis.\n']);
for c = 1:height(panelCWithinCTOmnibus)
    fprintf([ ...
        '%s: ANOVA F(%g,%g)=%.3f, p=%.4g; ' ...
        'Welch F(%g,%.2f)=%.3f, p=%.4g; KW p=%.4g; ' ...
        'range=%.4f; all-pairs equivalent IUT=%d, Holm=%d\n'], ...
        char(panelCWithinCTOmnibus.Condition(c)), ...
        panelCWithinCTOmnibus.ANOVA_DF1(c), ...
        panelCWithinCTOmnibus.ANOVA_DF2(c), ...
        panelCWithinCTOmnibus.ANOVA_F(c), ...
        panelCWithinCTOmnibus.ANOVA_p(c), ...
        panelCWithinCTOmnibus.WelchANOVA_DF1(c), ...
        panelCWithinCTOmnibus.WelchANOVA_DF2(c), ...
        panelCWithinCTOmnibus.WelchANOVA_F(c), ...
        panelCWithinCTOmnibus.WelchANOVA_p(c), ...
        panelCWithinCTOmnibus.KruskalWallis_p(c), ...
        panelCWithinCTOmnibus.ObservedMaxMinusMinMean(c), ...
        panelCWithinCTOmnibus.AllPairwiseEquivalent_IUT(c), ...
        panelCWithinCTOmnibus.AllPairwiseEquivalent_Holm(c));
end

%% ========================================================================
% PART 3A. Standalone Panel A: outer-bootstrap bagged peak distribution
% ========================================================================
if drawPanelA
    panelAOuterM = outerTbag_M(validOuterBagged);
    panelAOuterE = outerTbag_E(validOuterBagged);
    panelAEdges = makeCenteredEdges([panelAOuterM; panelAOuterE], ...
        panelAPeakBinWidth);

    panelAFig = figure('Name', ...
        'Panel A: outer-bootstrap bagged peak distribution', ...
        'NumberTitle', 'off', 'Color', 'w', 'Units', 'centimeters', ...
        'Position', [2, 2, 5, 4]);
    hold on;
    % Match the formal fig 2c bootstrap rendering order: E first, M last.
    histogram(panelAOuterE, panelAEdges, 'FaceColor', E_color, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8, ...
        'DisplayName', 'Ecell', 'Normalization', 'probability');
    histogram(panelAOuterM, panelAEdges, 'FaceColor', M_color, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.8, ...
        'DisplayName', 'Mcell', 'Normalization', 'probability');

    panelAAx = gca;
    panelAAx.Units = 'centimeters';
    panelAAx.Position(3:4) = [3.5, 2.6];
    panelAAx.LineWidth = 0.8;
    panelAAx.TickDir = 'out';
    panelAAx.TickLength = [0.02, 0.02];
    panelAAx.XTickLabelRotation = 0;
    panelAAx.XRuler.TickLabelGapOffset = -1;
    panelAAx.FontName = 'Arial';
    panelAAx.FontSize = 7;
    panelAAx.Color = 'none';
    ylim([0, 1.0]);
    xlim([10, 20]);
    xticks(10:5:20);
    yticks(0:0.5:1.0);
    panelAAx.XMinorTick = 'on';
    panelAAx.XAxis.MinorTickValues = setdiff(10:1:20, 10:5:20);
    grid off;
    box off;
    hold off;

    applyFigureTypography(panelAFig);
    panelAFig.PaperUnits = 'centimeters';
    panelAFig.PaperSize = [5, 4];
    panelAFig.PaperPosition = [0, 0, 5, 4];
    panelAFig.PaperPositionMode = 'manual';
    panelAFig.InvertHardcopy = 'off';
    drawnow;
    print(panelAFig, outputPanelASVG, '-dsvg', '-vector');
    fprintf('Standalone Panel A saved to: %s\n', outputPanelASVG);
end

%% ========================================================================
% PART 3C. Standalone Panel C: M/E amplitude panels sharing the CT axis
% ========================================================================
if drawPanelC
    panelCFig = plotPanelCAmplitudeByCondition(panelCCT, ...
        panelCActiveM_byCT, panelCActiveE_byCT, panelCCTStats, ...
        M_color, E_color, outputPanelCSVG, panelCUseBoxplot);
end

%% ========================================================================
% PART 3. Six-panel intensity figure in the original fig 2c visual style
% ========================================================================
if drawSixPanel

mainFig = figure('Name', 'fig 2c intensity six-panel analysis', ...
    'NumberTitle', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 18, 12]);
mainLayout = tiledlayout(mainFig, 2, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: primary total-output endpoint.
ax1 = nexttile(mainLayout, 1);
plotMetricByCT(ax1, ctList, outputM_byCT, outputE_byCT, ...
    M_color, E_color, '2-h-equivalent gated AUC (DeltaF/F*h)', ...
    'Event-gated transcriptional output');
addPanelLabel(ax1, 'A');

% Panel 2: event occupancy/duration component.
ax2 = nexttile(mainLayout, 2);
plotMetricByCT(ax2, ctList, dutyM_byCT, dutyE_byCT, ...
    M_color, E_color, 'Duty cycle', 'Active-time fraction');
ylim(ax2, [0, max(0.05, 1.08 * maxFiniteY(ax2))]);
addPanelLabel(ax2, 'B');

% Panel 3: conditional intensity component. Zero-event flies remain NaN.
ax3 = nexttile(mainLayout, 3);
plotMetricByCT(ax3, ctList, activeM_byCT, activeE_byCT, ...
    M_color, E_color, 'Active-state DeltaF/F', ...
    'Conditional amplitude (event = 1)');
addPanelLabel(ax3, 'C');

% Panel 4: paired fly-level M-E difference at every CT.
ax4 = nexttile(mainLayout, 4);
plotPairedDifferenceByCT(ax4, ctList, outputM_byCT, outputE_byCT, ...
    ctDifferenceTable, diffColor);
addPanelLabel(ax4, 'D');

% Panel 5: inner hard-peak probabilities whose expectations define Tbag.
ax5 = nexttile(mainLayout, 5);
barGroup = bar(ax5, ctList, [P_E(:), P_M(:)], 'grouped', ...
    'EdgeColor', 'none');
barGroup(1).FaceColor = E_color;
barGroup(2).FaceColor = M_color;
hold(ax5, 'on');
xline(ax5, Tbag_E, '--', 'E Tbag', 'Color', E_color, ...
    'LineWidth', 0.8, 'LabelVerticalAlignment', 'bottom');
xline(ax5, Tbag_M, '--', 'M Tbag', 'Color', M_color, ...
    'LineWidth', 0.8, 'LabelVerticalAlignment', 'top');
xlabel(ax5, 'Recording-window midpoint CT (h)');
ylabel(ax5, 'P_{inner}(hard peak = CT)');
title(ax5, 'Inner peak-location probability');
ylim(ax5, [0 1]);
legend(ax5, {'Ecell', 'Mcell'}, 'Location', 'best');
styleAxes(ax5, ctList);
addPanelLabel(ax5, 'E');

% Panel 6: continuous outer-bagged difference and MET boundaries.
ax6 = nexttile(mainLayout, 6);
plotOuterDifference(ax6, outerDiffValid, outerDifferenceBinWidth, ...
    baggedDiff, delta0, MET_z, MET_p);
addPanelLabel(ax6, 'F');

applyFigureTypography(mainFig);
mainFig.PaperUnits = 'centimeters';
mainFig.PaperSize = [18, 12];
mainFig.PaperPosition = [0, 0, 18, 12];
mainFig.PaperPositionMode = 'manual';
mainFig.InvertHardcopy = 'off';
drawnow;
print(mainFig, outputMainSVG, '-dsvg', '-vector');
fprintf('Six-panel figure saved to: %s\n', outputMainSVG);
end

%% ========================================================================
% PART 4. Six-panel diagnostics
% ========================================================================
if drawDiagnostics
diagnosticFig = figure('Name', 'fig 2c intensity diagnostics', ...
    'NumberTitle', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [3, 3, 18, 12]);
diagnosticLayout = tiledlayout(diagnosticFig, 2, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% Diagnostic 1: actual coverage, motivating duration normalization.
dax1 = nexttile(diagnosticLayout, 1);
plotMetricByCT(dax1, ctList, durationM_byCT, durationE_byCT, ...
    M_color, E_color, 'Observed duration (h)', ...
    'Actual recording coverage');
hold(dax1, 'on');
yline(dax1, recordingWindowHours, '--k', '2 h target', ...
    'LineWidth', 0.8);
addPanelLabel(dax1, 'A');

% Diagnostic 2: active/inactive intensity separation at fly level.
dax2 = nexttile(diagnosticLayout, 2);
intensityGroups = {flyMetricsTable.InactiveAmplitude_M(validPairRows), ...
    flyMetricsTable.ActiveAmplitude_M(validPairRows), ...
    flyMetricsTable.InactiveAmplitude_E(validPairRows), ...
    flyMetricsTable.ActiveAmplitude_E(validPairRows)};
plotGroupedValues(dax2, intensityGroups, ...
    {'M inactive', 'M active', 'E inactive', 'E active'}, ...
    {M_color * 0.55, M_color, E_color * 0.55, E_color});
ylabel(dax2, 'Time-weighted DeltaF/F');
title(dax2, 'Active versus inactive intensity');
addPanelLabel(dax2, 'B');

% Diagnostic 3: mechanistic duty-amplitude decomposition. These raw
% correlations are descriptive and can be CT-confounded.
dax3 = nexttile(diagnosticLayout, 3);
dutyM_all = flyMetricsTable.DutyCycle_M(validPairRows);
ampM_all = flyMetricsTable.ActiveAmplitude_M(validPairRows);
dutyE_all = flyMetricsTable.DutyCycle_E(validPairRows);
ampE_all = flyMetricsTable.ActiveAmplitude_E(validPairRows);
scatter(dax3, dutyE_all, ampE_all, 18, E_color, 'filled', ...
    'MarkerFaceAlpha', 0.55, 'DisplayName', 'Ecell');
hold(dax3, 'on');
scatter(dax3, dutyM_all, ampM_all, 18, M_color, 'filled', ...
    'MarkerFaceAlpha', 0.55, 'DisplayName', 'Mcell');
rhoM = safeSpearman(dutyM_all, ampM_all);
rhoE = safeSpearman(dutyE_all, ampE_all);
xlabel(dax3, 'Duty cycle');
ylabel(dax3, 'Active-state DeltaF/F');
title(dax3, sprintf('Duty-amplitude relation (rho_M=%.2f, rho_E=%.2f)', ...
    rhoM, rhoE));
legend(dax3, 'Location', 'best');
styleAxes(dax3, []);
addPanelLabel(dax3, 'C');

% Diagnostic 4: normal approximation used by the Wald-type MET.
dax4 = nexttile(diagnosticLayout, 4);
qqplot(dax4, outerDiffValid);
xlabel(dax4, 'Standard-normal quantiles');
ylabel(dax4, 'Outer bagged-difference quantiles (h)');
title(dax4, sprintf('Outer-difference Q-Q plot (JB p=%.3g)', JB_p));
styleAxes(dax4, []);
addPanelLabel(dax4, 'D');

% Diagnostic 5: pairing-induced covariance of outer peak estimates.
dax5 = nexttile(diagnosticLayout, 5);
scatter(dax5, outerTbag_M(validOuterBagged), ...
    outerTbag_E(validOuterBagged), 10, [0.25 0.25 0.25], ...
    'filled', 'MarkerFaceAlpha', 0.25);
hold(dax5, 'on');
xline(dax5, Tbag_M, '--', 'M estimate', 'Color', M_color);
yline(dax5, Tbag_E, '--', 'E estimate', 'Color', E_color);
xlabel(dax5, 'Outer Tbag M (h)');
ylabel(dax5, 'Outer Tbag E (h)');
title(dax5, sprintf('Paired outer estimates (r=%.3f)', ...
    outerCorrelation_M_E));
styleAxes(dax5, []);
addPanelLabel(dax5, 'E');

% Diagnostic 6: cumulative outer Monte Carlo behavior.
dax6 = nexttile(diagnosticLayout, 6);
nValidOuter = numel(outerDiffValid);
runningMean = cumsum(outerDiffValid) ./ (1:nValidOuter).';
runningSD = NaN(nValidOuter, 1);
for k = 2:nValidOuter
    runningSD(k) = std(outerDiffValid(1:k), 'omitnan');
end
yyaxis(dax6, 'left');
plot(dax6, 1:nValidOuter, runningMean, 'Color', [0.2 0.2 0.2], ...
    'LineWidth', 0.9);
yline(dax6, baggedDiff, '--', 'Observed diff', 'Color', [0.75 0 0]);
ylabel(dax6, 'Cumulative mean difference (h)');
yyaxis(dax6, 'right');
plot(dax6, 1:nValidOuter, runningSD, 'Color', M_color, ...
    'LineWidth', 0.9);
yline(dax6, baggedSE, '--', 'Final SE', 'Color', M_color);
ylabel(dax6, 'Cumulative bootstrap SE (h)');
xlabel(dax6, 'Number of outer replicates');
title(dax6, 'Outer Monte Carlo convergence');
styleAxes(dax6, []);
addPanelLabel(dax6, 'F');

applyFigureTypography(diagnosticFig);
diagnosticFig.PaperUnits = 'centimeters';
diagnosticFig.PaperSize = [18, 12];
diagnosticFig.PaperPosition = [0, 0, 18, 12];
diagnosticFig.PaperPositionMode = 'manual';
diagnosticFig.InvertHardcopy = 'off';
drawnow;
print(diagnosticFig, outputDiagnosticSVG, '-dsvg', '-vector');
fprintf('Diagnostic figure saved to: %s\n', outputDiagnosticSVG);
end

%% ========================================================================
% Local helper functions
% ========================================================================

publicResult = struct('statistics', statsBagged, 'outerBootstrap', outerBootstrapTable, 'flyData', flyMetricsTable, 'cellData', cellMetricsTable, 'ctSummary', ctSummaryTable, 'panelCStatistics', panelCCTStats, 'panelCPairedTests', panelCPairedTests);
publicResult.configuration = cfg;
if exist('outputMAT', 'var') && isfile(outputMAT)
    save(outputMAT, 'cfg', 'publicResult', '-append');
else
    save(fullfile(cfg.outputDir, cfg.outputStem + "_results.mat"), 'cfg', 'publicResult');
end
end

function folders = normalizeIntensitySourceFolders(folders)
% Accept one char path, a string array, or a cell array of paths. Canonical
% paths prevent the same root being added twice under relative/absolute names.
    if ischar(folders)
        folders = string({folders});
    else
        folders = string(folders);
    end
    folders = strtrim(folders(:));
    if isempty(folders) || any(ismissing(folders) | strlength(folders) == 0)
        error('fig2cIntensity:EmptySourceFolder', ...
            'parentFolders must contain at least one nonempty folder path.');
    end
    for s = 1:numel(folders)
        pathObject = javaObject('java.io.File', char(folders(s)));
        folders(s) = string(pathObject.getCanonicalPath());
    end
    if numel(unique(lower(folders))) ~= numel(folders)
        error('fig2cIntensity:DuplicateSourceFolder', ...
            'The same source folder occurs more than once in parentFolders.');
    end
end

function tagged = tagIntensitySource(tagged, sourceIndex, sourceFolder)
% Keep the original workbook-derived identifier alongside a source-qualified
% FlyID. All measurements for an M/E fly pair remain in the same table row.
    tagged.SourceIndex = repmat(sourceIndex, height(tagged), 1);
    tagged.SourceFolder = repmat(sourceFolder, height(tagged), 1);
    tagged.OriginalFlyID = tagged.FlyID;
    tagged.FlyID = "source" + string(sourceIndex) + "::" + tagged.FlyID;
end

function [cellTable, flyTable, auditTable, sourceTable] = ...
        loadPairedIntensityData(parentFolders, recordingWindowHours, ...
        timestampToleranceHours)
% Reuse the original single-source preprocessing independently for each
% batch, then pool fly records within CT. The sampling unit stays the fly;
% sources are provenance labels, not a new bootstrap sampling level.
    parentFolders = normalizeIntensitySourceFolders(parentFolders);
    nSources = numel(parentFolders);
    cellParts = cell(nSources, 1);
    flyParts = cell(nSources, 1);
    auditParts = cell(nSources, 1);
    sourceNames = strings(nSources, 1);
    nFly = zeros(nSources, 1);
    nValid = zeros(nSources, 1);
    nCells = zeros(nSources, 1);
    nCT = zeros(nSources, 1);
    nIgnored = zeros(nSources, 1);
    nEmptyCT = zeros(nSources, 1);
    for s = 1:nSources
        sourceFolder = parentFolders(s);
        if ~bootstrapRawIsFolder(sourceFolder)
            error('fig2cIntensity:MissingDataFolder', ...
                'Source folder %d not found: %s', s, sourceFolder);
        end
        [~, sourceName] = fileparts(sourceFolder);
        sourceNames(s) = string(sourceName);
        fprintf('  Source %d/%d: %s\n', s, nSources, sourceFolder);
        [cells, flies, audit, discovery] = loadSinglePairedIntensityData( ...
            sourceFolder, recordingWindowHours, timestampToleranceHours);
        nFly(s) = height(flies);
        nValid(s) = sum(flies.ValidPrimaryPair);
        nCells(s) = height(cells);
        nCT(s) = numel(unique(flies.AnalysisTime_midpoint_h));
        nIgnored(s) = discovery.NIgnoredXLSX;
        nEmptyCT(s) = discovery.NCTfoldersWithoutPairs;
        cellParts{s} = tagIntensitySource(cells, s, sourceFolder);
        flyParts{s} = tagIntensitySource(flies, s, sourceFolder);
        auditParts{s} = tagIntensitySource(audit, s, sourceFolder);
    end
    cellTable = vertcat(cellParts{:});
    flyTable = vertcat(flyParts{:});
    auditTable = vertcat(auditParts{:});

    % Matching filenames in different roots can be copies of the same fly.
    % Stop rather than count one experimental recording twice or choose a
    % version silently. Distinct flies need distinct recording-ID identifiers.
    recordKeys = lower(flyTable.OriginalFlyID) + "|CT" + ...
        compose('%.12g', flyTable.RecordingStartCT_h);
    [uniqueKeys, ~, keyGroups] = unique(recordKeys);
    keyCounts = accumarray(keyGroups, 1);
    repeated = find(keyCounts > 1, 1);
    if ~isempty(repeated)
        duplicateRows = keyGroups == repeated;
        error('fig2cIntensity:DuplicateRecording', ...
            ['Recording %s occurs more than once. Resolve repeated ' ...
             'recordings before pooling these files:\n%s'], ...
            uniqueKeys(repeated), ...
            strjoin(flyTable.SignalFile(duplicateRows), newline));
    end
    cellTable = sortrows(cellTable, {'AnalysisTime_midpoint_h', ...
        'FlyID', 'Condition', 'SheetName'});
    flyTable = sortrows(flyTable, {'AnalysisTime_midpoint_h', 'FlyID'});
    auditTable = sortrows(auditTable, {'AnalysisTime_midpoint_h', 'FlyID'});
    sourceTable = table((1:nSources).', sourceNames, parentFolders, ...
        nFly, nValid, nCells, nCT, nIgnored, nEmptyCT, ...
        'VariableNames', {'SourceIndex', 'SourceName', 'SourceFolder', ...
        'N_flyRecords', 'N_validPairedFlies', 'N_cellSheets', ...
        'N_loadedCTgroups', 'N_ignoredXLSX', 'N_CTfoldersWithoutPairs'});
end

function [cellTable, flyTable, auditTable, discovery] = loadSinglePairedIntensityData( ...
        parentFolder, recordingWindowHours, timestampToleranceHours)
% Discover exact signal/intensity workbook pairs, align sheets/timestamps,
% calculate time-weighted cell metrics, and aggregate cell-time within fly.

    subfolders = bootstrapRawFolders(parentFolder);
    subfolders = subfolders([subfolders.isdir] & ...
        ~ismember({subfolders.name}, {'.', '..'}));

    pairTemplate = struct('StartCT', NaN, 'EndCT', NaN, ...
        'MidCT', NaN, 'FlyKey', "", 'SignalPath', "", ...
        'IntensityPath', "", 'SignalName', "", 'IntensityName', "", ...
        'Sheets', strings(0, 1));
    allPairs = repmat(pairTemplate, 0, 1);
    totalSheets = 0;
    discovery = struct('NIgnoredXLSX', 0, 'NCTfoldersWithoutPairs', 0);

    for k = 1:numel(subfolders)
        folderName = subfolders(k).name;
        ctToken = regexp(folderName, 'CT(\d+\.?\d*)', ...
            'tokens', 'once', 'ignorecase');
        if isempty(ctToken)
            continue;
        end
        startCT = str2double(ctToken{1});
        if ~isfinite(startCT)
            continue;
        end
        folderPath = fullfile(parentFolder, folderName);
        [folderPairs, nIgnored] = pairWorkbooksExactly(folderPath);
        discovery.NIgnoredXLSX = discovery.NIgnoredXLSX + nIgnored;
        if isempty(folderPairs)
            discovery.NCTfoldersWithoutPairs = ...
                discovery.NCTfoldersWithoutPairs + 1;
            if nIgnored > 0
                warning('fig2cIntensity:UnrecognizedCTWorkbooks', ...
                    ['No recognized signal/intensity workbook pairs in %s; ' ...
                     '%d XLSX files do not match the input naming rules. ' ...
                     'Expected: RecordingNNNN_binary.xlsx and ' ...
                     'RecordingNNNN_intensity.xlsx.'], folderPath, nIgnored);
            end
        end
        for f = 1:numel(folderPairs)
            signalSheets = string(bootstrapRawSheets(folderPairs(f).SignalPath));
            intensitySheets = string(bootstrapRawSheets(folderPairs(f).IntensityPath));
            if numel(signalSheets) ~= numel(intensitySheets) || ...
                    ~isempty(setxor(lower(signalSheets), ...
                    lower(intensitySheets)))
                error('fig2cIntensity:SheetSetMismatch', ...
                    ['Signal/intensity sheet sets differ for fly %s.\n' ...
                     'Signal: %s\nIntensity: %s'], folderPairs(f).FlyKey, ...
                    folderPairs(f).SignalPath, ...
                    folderPairs(f).IntensityPath);
            end

            pairRecord = pairTemplate;
            pairRecord.StartCT = startCT;
            pairRecord.EndCT = startCT + recordingWindowHours;
            pairRecord.MidCT = startCT + recordingWindowHours / 2;
            pairRecord.FlyKey = folderPairs(f).FlyKey;
            pairRecord.SignalPath = folderPairs(f).SignalPath;
            pairRecord.IntensityPath = folderPairs(f).IntensityPath;
            pairRecord.SignalName = folderPairs(f).SignalName;
            pairRecord.IntensityName = folderPairs(f).IntensityName;
            pairRecord.Sheets = signalSheets(:);
            allPairs(end + 1, 1) = pairRecord; %#ok<AGROW>
            totalSheets = totalSheets + numel(signalSheets);
        end
    end

    nFly = numel(allPairs);
    if nFly == 0 || totalSheets == 0
        error('fig2cIntensity:NoWorkbookPairs', ...
            'No exact signal/intensity workbook pairs were found in %s.', ...
            parentFolder);
    end

    cellTemplate = makeEmptyCellRecord();
    flyTemplate = makeEmptyFlyRecord();
    auditTemplate = struct('RecordingStartCT_h', NaN, ...
        'AnalysisTime_midpoint_h', NaN, 'FlyID', "", ...
        'SignalFile', "", 'IntensityFile', "", 'N_sheets', 0, ...
        'AllSheetsMatched', false, 'N_unmatchedTimestamps', 0);
    cellRecords = repmat(cellTemplate, totalSheets, 1);
    flyRecords = repmat(flyTemplate, nFly, 1);
    auditRecords = repmat(auditTemplate, nFly, 1);

    cellIndex = 0;
    for f = 1:nFly
        pair = allPairs(f);
        nSheets = numel(pair.Sheets);
        metricsThisFly = repmat(makeEmptyCellMetric(), nSheets, 1);
        isEThisFly = false(nSheets, 1);
        unmatchedThisFly = 0;

        intensitySheets = string(bootstrapRawSheets(pair.IntensityPath));
        for s = 1:nSheets
            sheetName = pair.Sheets(s);
            intensityIndex = find(strcmpi(intensitySheets, sheetName), 1);
            if isempty(intensityIndex)
                error('fig2cIntensity:MissingIntensitySheet', ...
                    'No intensity sheet %s for fly %s.', ...
                    sheetName, pair.FlyKey);
            end
            intensitySheet = intensitySheets(intensityIndex);
            [metrics, unmatchedCount] = readAlignedCellMetrics( ...
                pair.SignalPath, pair.IntensityPath, sheetName, ...
                intensitySheet, recordingWindowHours, ...
                timestampToleranceHours);
            metricsThisFly(s) = metrics;
            unmatchedThisFly = unmatchedThisFly + unmatchedCount;
            isE = contains(sheetName, '5th', 'IgnoreCase', true);
            isEThisFly(s) = isE;

            cellIndex = cellIndex + 1;
            record = cellTemplate;
            record.RecordingStartCT_h = pair.StartCT;
            record.RecordingEndCT_h = pair.EndCT;
            record.AnalysisTime_midpoint_h = pair.MidCT;
            record.FlyID = pair.FlyKey;
            record.SignalFile = pair.SignalPath;
            record.IntensityFile = pair.IntensityPath;
            record.SheetName = sheetName;
            if isE
                record.Condition = "E";
            else
                record.Condition = "M";
            end
            metricNames = fieldnames(metrics);
            for q = 1:numel(metricNames)
                record.(metricNames{q}) = metrics.(metricNames{q});
            end
            cellRecords(cellIndex) = record;
        end

        Mmetrics = metricsThisFly(~isEThisFly);
        Emetrics = metricsThisFly(isEThisFly);
        Magg = aggregateCellMetrics(Mmetrics, recordingWindowHours);
        Eagg = aggregateCellMetrics(Emetrics, recordingWindowHours);

        flyRecord = flyTemplate;
        flyRecord.RecordingStartCT_h = pair.StartCT;
        flyRecord.RecordingEndCT_h = pair.EndCT;
        flyRecord.AnalysisTime_midpoint_h = pair.MidCT;
        flyRecord.FlyID = pair.FlyKey;
        flyRecord.SignalFile = pair.SignalPath;
        flyRecord.IntensityFile = pair.IntensityPath;
        flyRecord.N_M_cells = numel(Mmetrics);
        flyRecord.N_E_cells = numel(Emetrics);
        flyRecord.N_active_M_cells = Magg.NActiveCells;
        flyRecord.N_active_E_cells = Eagg.NActiveCells;
        flyRecord = assignAggregateToFly(flyRecord, Magg, 'M');
        flyRecord = assignAggregateToFly(flyRecord, Eagg, 'E');
        flyRecord.ValidPrimaryPair = Magg.ValidPrimary && Eagg.ValidPrimary;
        flyRecords(f) = flyRecord;

        if numel(Mmetrics) ~= 4 || numel(Emetrics) ~= 1
            warning('fig2cIntensity:UnexpectedCellCount', ...
                '%s has %d M sheets and %d E sheets.', ...
                pair.SignalPath, numel(Mmetrics), numel(Emetrics));
        end
        if unmatchedThisFly > 0
            error('fig2cIntensity:UnmatchedTimestamps', ...
                ['%d signal timestamps in the analysis window had no exact ' ...
                 'intensity match for fly %s.'], ...
                unmatchedThisFly, pair.FlyKey);
        end

        auditRecord = auditTemplate;
        auditRecord.RecordingStartCT_h = pair.StartCT;
        auditRecord.AnalysisTime_midpoint_h = pair.MidCT;
        auditRecord.FlyID = pair.FlyKey;
        auditRecord.SignalFile = pair.SignalPath;
        auditRecord.IntensityFile = pair.IntensityPath;
        auditRecord.N_sheets = nSheets;
        auditRecord.AllSheetsMatched = true;
        auditRecord.N_unmatchedTimestamps = unmatchedThisFly;
        auditRecords(f) = auditRecord;
    end

    if cellIndex ~= totalSheets
        error('fig2cIntensity:InternalCountMismatch', ...
            'Internal sheet preallocation count mismatch.');
    end

    cellTable = struct2table(cellRecords);
    flyTable = struct2table(flyRecords);
    auditTable = struct2table(auditRecords);
    cellTable = sortrows(cellTable, {'AnalysisTime_midpoint_h', ...
        'FlyID', 'Condition', 'SheetName'});
    flyTable = sortrows(flyTable, ...
        {'AnalysisTime_midpoint_h', 'FlyID'});
    auditTable = sortrows(auditTable, ...
        {'AnalysisTime_midpoint_h', 'FlyID'});
end

function [pairs, nIgnored] = pairWorkbooksExactly(folderPath)
% Exact normalized recording-ID key matching; no contains/prefix matching.
    nIgnored = 0;
    allFiles = bootstrapRawFiles(folderPath);
    if isempty(allFiles)
        pairs = repmat(struct('FlyKey', "", 'SignalPath', "", ...
            'IntensityPath', "", 'SignalName', "", ...
            'IntensityName', ""), 0, 1);
        return;
    end
    names = string({allFiles.name});
    keep = ~startsWith(names, '~$');
    allFiles = allFiles(keep);
    names = string({allFiles.name});
    signalMask = ~cellfun(@isempty, regexpi(cellstr(names), ...
        '_binary\.xlsx$'));
    intensityMask = ~cellfun(@isempty, regexpi(cellstr(names), ...
        '_intensity\.xlsx$'));
    nIgnored = sum(~signalMask & ~intensityMask);
    signalFiles = allFiles(signalMask);
    intensityFiles = allFiles(intensityMask);

    signalMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    intensityMap = containers.Map('KeyType', 'char', 'ValueType', 'double');
    signalKeys = strings(numel(signalFiles), 1);
    intensityKeys = strings(numel(intensityFiles), 1);
    for i = 1:numel(signalFiles)
        signalKeys(i) = normalizedFlyKey(signalFiles(i).name, 'signal');
        keyChar = char(signalKeys(i));
        if isKey(signalMap, keyChar)
            error('fig2cIntensity:DuplicateSignalKey', ...
                'Duplicate signal key %s in %s.', signalKeys(i), folderPath);
        end
        signalMap(keyChar) = i;
    end
    for i = 1:numel(intensityFiles)
        intensityKeys(i) = normalizedFlyKey( ...
            intensityFiles(i).name, 'intensity');
        keyChar = char(intensityKeys(i));
        if isKey(intensityMap, keyChar)
            error('fig2cIntensity:DuplicateIntensityKey', ...
                'Duplicate intensity key %s in %s.', ...
                intensityKeys(i), folderPath);
        end
        intensityMap(keyChar) = i;
    end

    missingIntensity = setdiff(signalKeys, intensityKeys);
    orphanIntensity = setdiff(intensityKeys, signalKeys);
    if ~isempty(missingIntensity) || ~isempty(orphanIntensity)
        error('fig2cIntensity:WorkbookPairMismatch', ...
            ['Exact signal/intensity key mismatch in %s. Missing intensity: ' ...
             '%s. Orphan intensity: %s.'], folderPath, ...
            strjoin(missingIntensity, ', '), strjoin(orphanIntensity, ', '));
    end

    sortedKeys = sort(signalKeys);
    pairTemplate = struct('FlyKey', "", 'SignalPath', "", ...
        'IntensityPath', "", 'SignalName', "", 'IntensityName', "");
    pairs = repmat(pairTemplate, numel(sortedKeys), 1);
    for i = 1:numel(sortedKeys)
        keyChar = char(sortedKeys(i));
        sf = signalFiles(signalMap(keyChar));
        inf = intensityFiles(intensityMap(keyChar));
        pairs(i).FlyKey = sortedKeys(i);
        pairs(i).SignalPath = string(fullfile(folderPath, sf.name));
        pairs(i).IntensityPath = string(fullfile(folderPath, inf.name));
        pairs(i).SignalName = string(sf.name);
        pairs(i).IntensityName = string(inf.name);
    end
end

function key = normalizedFlyKey(fileName, fileType)
    [~, baseName, ~] = fileparts(fileName);
    if strcmpi(fileType, 'signal')
        baseName = regexprep(baseName, ...
            '_binary$', '', 'ignorecase');
    else
        baseName = regexprep(baseName, ...
            '_intensity$', '', 'ignorecase');
    end
    baseName = regexprep(strtrim(baseName), '\s+', ' ');
    key = lower(string(baseName));
end

function [metrics, unmatchedCount] = readAlignedCellMetrics( ...
        signalPath, intensityPath, signalSheet, intensitySheet, ...
        recordingWindowHours, timestampToleranceHours)
    try
        signalTable = bootstrapRawReadSheet(signalPath, signalSheet);
        intensityTable = bootstrapRawReadSheet(intensityPath, intensitySheet);
    catch ME
        error('fig2cIntensity:SheetReadFailed', ...
            'Cannot read paired sheets for %s: %s', signalPath, ME.message);
    end
    if width(signalTable) < 2 || width(intensityTable) < 2
        error('fig2cIntensity:TooFewColumns', ...
            'Both paired sheets must contain time and value columns.');
    end
    signalTime = signalTable{:, 1};
    eventValues = signalTable{:, 2};
    intensityTime = intensityTable{:, 1};
    intensityValues = intensityTable{:, 2};
    if ~isnumeric(signalTime) || ~isnumeric(eventValues) || ...
            ~isnumeric(intensityTime) || ~isnumeric(intensityValues)
        error('fig2cIntensity:NonNumericColumns', ...
            'Time, signal, and intensity columns must be numeric.');
    end
    signalTime = signalTime(:);
    eventValues = eventValues(:);
    intensityTime = intensityTime(:);
    intensityValues = intensityValues(:);
    metrics = makeEmptyCellMetric();
    metrics.NRawSignalRows = numel(signalTime);
    metrics.NRawIntensityRows = numel(intensityTime);
    unmatchedCount = 0;
    if isempty(signalTime)
        return;
    end

    firstValidTime = find(isfinite(signalTime), 1, 'first');
    if isempty(firstValidTime)
        return;
    end
    relativeTime = signalTime - signalTime(firstValidTime);
    keep = isfinite(relativeTime) & relativeTime >= 0 & ...
        relativeTime <= recordingWindowHours + timestampToleranceHours;
    signalTimeWindow = signalTime(keep);
    relativeTimeWindow = relativeTime(keep);
    eventWindow = eventValues(keep);
    metrics.NRowsInWindow = numel(relativeTimeWindow);
    if isempty(relativeTimeWindow)
        return;
    end

    [matched, intensityLocation] = ismembertol(signalTimeWindow, ...
        intensityTime, timestampToleranceHours, 'DataScale', 1);
    unmatchedCount = sum(~matched);
    if unmatchedCount > 0
        return;
    end
    intensityWindow = intensityValues(intensityLocation);

    [relativeTimeWindow, order] = sort(relativeTimeWindow);
    eventWindow = eventWindow(order);
    intensityWindow = intensityWindow(order);
    if any(diff(relativeTimeWindow) <= 0)
        error('fig2cIntensity:NonIncreasingTime', ...
            '%s [%s] contains duplicate or non-increasing timestamps.', ...
            signalPath, signalSheet);
    end
    metrics = calculateCellMetrics(relativeTimeWindow, eventWindow, ...
        intensityWindow, metrics);
end

function metrics = calculateCellMetrics(time, eventValues, intensity, metrics)
    if numel(time) < 2
        return;
    end
    validEvent = isfinite(eventValues);
    if ~all(validEvent)
        time = time(validEvent);
        intensity = intensity(validEvent);
        eventValues = eventValues(validEvent);
    end
    if numel(time) < 2
        return;
    end

    if any(abs(eventValues - round(eventValues)) > 1e-12 | ...
            (eventValues ~= 0 & eventValues ~= 1))
        warning('fig2cIntensity:NonBinarySignal', ...
            'Signal values were not exactly 0/1 and were binarized as >0.');
    end
    event = eventValues > 0;
    if any(event & ~isfinite(intensity))
        warning('fig2cIntensity:MissingActiveIntensity', ...
            ['At least one active frame has missing intensity; primary ' ...
             'metrics for this cell are set to NaN.']);
        return;
    end

    weights = trapezoidWeights(time);
    duration = sum(weights);
    if ~isfinite(duration) || duration <= 0
        return;
    end
    activeWeights = weights .* double(event);
    activeDuration = sum(activeWeights);

    gatedIntensity = zeros(size(intensity));
    gatedIntensity(event) = intensity(event);
    gatedAUC = sum(weights .* gatedIntensity, 'omitnan');
    if activeDuration > 0
        activeAmplitude = gatedAUC / activeDuration;
    else
        activeAmplitude = NaN;
    end

    finiteInactive = ~event & isfinite(intensity);
    inactiveWeightObserved = sum(weights(finiteInactive));
    if inactiveWeightObserved > 0
        inactiveAUCObserved = sum(weights(finiteInactive) .* ...
            intensity(finiteInactive));
        inactiveAmplitude = inactiveAUCObserved / inactiveWeightObserved;
    else
        inactiveAmplitude = NaN;
    end

    finiteIntensity = isfinite(intensity);
    finiteIntensityDuration = sum(weights(finiteIntensity));
    if finiteIntensityDuration > 0
        rawIntensityAUC = sum(weights(finiteIntensity) .* ...
            intensity(finiteIntensity));
        rawMeanIntensity = rawIntensityAUC / finiteIntensityDuration;
    else
        rawIntensityAUC = NaN;
        rawMeanIntensity = NaN;
    end

    dutyCycle = activeDuration / duration;
    gatedOutputRate = gatedAUC / duration;
    if isfinite(activeAmplitude)
        decompositionError = gatedOutputRate - dutyCycle * activeAmplitude;
    else
        decompositionError = NaN;
    end
    activeContrast = activeAmplitude - inactiveAmplitude;
    baselineCorrectedGatedRate = dutyCycle * activeContrast;

    metrics.NRowsUsed = numel(time);
    metrics.NActiveFrames = sum(event);
    metrics.NInactiveFrames = sum(~event);
    metrics.NBouts = sum(diff([false; event]) == 1);
    metrics.ObservedDuration_h = duration;
    metrics.MedianSamplingInterval_min = median(diff(time)) * 60;
    metrics.EventWeightedDuration_h = activeDuration;
    metrics.InactiveWeightedDuration_h = inactiveWeightObserved;
    metrics.DutyCycle = dutyCycle;
    metrics.ObservedGatedAUC = gatedAUC;
    metrics.GatedOutputRate = gatedOutputRate;
    metrics.ActiveAmplitude = activeAmplitude;
    metrics.InactiveAmplitude = inactiveAmplitude;
    metrics.ActiveContrast = activeContrast;
    metrics.BaselineCorrectedGatedRate = baselineCorrectedGatedRate;
    metrics.RawIntensityAUC = rawIntensityAUC;
    metrics.RawMeanIntensity = rawMeanIntensity;
    metrics.DecompositionError = decompositionError;
    metrics.ValidPrimary = isfinite(gatedOutputRate);
end

function weights = trapezoidWeights(time)
% Dot product weights equivalent to trapz(time, value).
    time = time(:);
    n = numel(time);
    weights = zeros(n, 1);
    if n < 2
        return;
    end
    weights(1) = (time(2) - time(1)) / 2;
    weights(n) = (time(n) - time(n - 1)) / 2;
    if n > 2
        weights(2:n-1) = (time(3:n) - time(1:n-2)) / 2;
    end
end

function aggregate = aggregateCellMetrics(metrics, recordingWindowHours)
    aggregate = makeEmptyAggregate();
    if isempty(metrics)
        return;
    end
    valid = [metrics.ValidPrimary].';
    metrics = metrics(valid);
    if isempty(metrics)
        return;
    end

    duration = [metrics.ObservedDuration_h].';
    eventDuration = [metrics.EventWeightedDuration_h].';
    gatedAUC = [metrics.ObservedGatedAUC].';
    gatedRate = [metrics.GatedOutputRate].';
    dutyCycle = [metrics.DutyCycle].';
    if ~any(isfinite(duration)) || ~any(isfinite(gatedRate))
        return;
    end

    aggregate.MeanObservedDuration_h = mean(duration, 'omitnan');
    aggregate.ObservedGatedAUC_perCell = mean(gatedAUC, 'omitnan');
    % Equal cell weighting reproduces the original within-fly rule: first
    % summarize each M sheet, then average the four M-cell summaries.
    aggregate.GatedOutputRate = mean(gatedRate, 'omitnan');
    aggregate.GatedAUC_2hEquivalent = ...
        aggregate.GatedOutputRate * recordingWindowHours;
    aggregate.DutyCycle = mean(dutyCycle, 'omitnan');
    if aggregate.DutyCycle > 0
        aggregate.ActiveAmplitude = ...
            aggregate.GatedOutputRate / aggregate.DutyCycle;
    end

    inactiveAmplitude = [metrics.InactiveAmplitude].';
    inactiveFraction = 1 - dutyCycle;
    validInactive = isfinite(inactiveAmplitude) & ...
        isfinite(inactiveFraction) & inactiveFraction > 0;
    if any(validInactive)
        aggregate.InactiveAmplitude = ...
            sum(inactiveAmplitude(validInactive) .* ...
            inactiveFraction(validInactive)) / ...
            sum(inactiveFraction(validInactive));
    end
    aggregate.ActiveContrast = ...
        aggregate.ActiveAmplitude - aggregate.InactiveAmplitude;
    aggregate.BaselineCorrectedGatedRate = ...
        aggregate.DutyCycle * aggregate.ActiveContrast;

    rawMean = [metrics.RawMeanIntensity].';
    validRaw = isfinite(rawMean);
    if any(validRaw)
        aggregate.RawMeanIntensity = mean(rawMean(validRaw));
    end
    aggregate.NActiveCells = sum(eventDuration > 0);
    aggregate.WithinCellSD_GatedOutputRate = ...
        std([metrics.GatedOutputRate].', 'omitnan');
    aggregate.WithinCellSD_DutyCycle = ...
        std([metrics.DutyCycle].', 'omitnan');
    aggregate.WithinCellSD_ActiveAmplitude = ...
        std([metrics.ActiveAmplitude].', 'omitnan');
    aggregate.DecompositionError = aggregate.GatedOutputRate - ...
        aggregate.DutyCycle * aggregate.ActiveAmplitude;
    aggregate.ValidPrimary = isfinite(aggregate.GatedOutputRate);
end

function flyRecord = assignAggregateToFly(flyRecord, aggregate, suffix)
    fieldMap = {'MeanObservedDuration_h', 'ObservedGatedAUC_perCell', ...
        'GatedOutputRate', 'GatedAUC_2hEquivalent', 'DutyCycle', ...
        'ActiveAmplitude', 'InactiveAmplitude', 'ActiveContrast', ...
        'BaselineCorrectedGatedRate', 'RawMeanIntensity', ...
        'WithinCellSD_GatedOutputRate', 'WithinCellSD_DutyCycle', ...
        'WithinCellSD_ActiveAmplitude', 'DecompositionError'};
    for i = 1:numel(fieldMap)
        destination = [fieldMap{i} '_' suffix];
        flyRecord.(destination) = aggregate.(fieldMap{i});
    end
end

function metric = makeEmptyCellMetric()
    metric = struct('NRawSignalRows', 0, 'NRawIntensityRows', 0, ...
        'NRowsInWindow', 0, 'NRowsUsed', 0, 'NActiveFrames', 0, ...
        'NInactiveFrames', 0, 'NBouts', 0, 'ObservedDuration_h', NaN, ...
        'MedianSamplingInterval_min', NaN, ...
        'EventWeightedDuration_h', NaN, ...
        'InactiveWeightedDuration_h', NaN, 'DutyCycle', NaN, ...
        'ObservedGatedAUC', NaN, 'GatedOutputRate', NaN, ...
        'ActiveAmplitude', NaN, 'InactiveAmplitude', NaN, ...
        'ActiveContrast', NaN, 'BaselineCorrectedGatedRate', NaN, ...
        'RawIntensityAUC', NaN, 'RawMeanIntensity', NaN, ...
        'DecompositionError', NaN, 'ValidPrimary', false);
end

function record = makeEmptyCellRecord()
    record = struct('RecordingStartCT_h', NaN, 'RecordingEndCT_h', NaN, ...
        'AnalysisTime_midpoint_h', NaN, 'FlyID', "", ...
        'SignalFile', "", 'IntensityFile', "", 'SheetName', "", ...
        'Condition', "", 'NRawSignalRows', 0, ...
        'NRawIntensityRows', 0, 'NRowsInWindow', 0, 'NRowsUsed', 0, ...
        'NActiveFrames', 0, 'NInactiveFrames', 0, 'NBouts', 0, ...
        'ObservedDuration_h', NaN, 'MedianSamplingInterval_min', NaN, ...
        'EventWeightedDuration_h', NaN, ...
        'InactiveWeightedDuration_h', NaN, 'DutyCycle', NaN, ...
        'ObservedGatedAUC', NaN, 'GatedOutputRate', NaN, ...
        'ActiveAmplitude', NaN, 'InactiveAmplitude', NaN, ...
        'ActiveContrast', NaN, 'BaselineCorrectedGatedRate', NaN, ...
        'RawIntensityAUC', NaN, 'RawMeanIntensity', NaN, ...
        'DecompositionError', NaN, 'ValidPrimary', false);
end

function aggregate = makeEmptyAggregate()
    aggregate = struct('MeanObservedDuration_h', NaN, ...
        'ObservedGatedAUC_perCell', NaN, 'GatedOutputRate', NaN, ...
        'GatedAUC_2hEquivalent', NaN, 'DutyCycle', NaN, ...
        'ActiveAmplitude', NaN, 'InactiveAmplitude', NaN, ...
        'ActiveContrast', NaN, 'BaselineCorrectedGatedRate', NaN, ...
        'RawMeanIntensity', NaN, 'NActiveCells', 0, ...
        'WithinCellSD_GatedOutputRate', NaN, ...
        'WithinCellSD_DutyCycle', NaN, ...
        'WithinCellSD_ActiveAmplitude', NaN, ...
        'DecompositionError', NaN, 'ValidPrimary', false);
end

function record = makeEmptyFlyRecord()
    record = struct('RecordingStartCT_h', NaN, 'RecordingEndCT_h', NaN, ...
        'AnalysisTime_midpoint_h', NaN, 'FlyID', "", ...
        'SignalFile', "", 'IntensityFile', "", ...
        'N_M_cells', 0, 'N_E_cells', 0, 'N_active_M_cells', 0, ...
        'N_active_E_cells', 0, ...
        'MeanObservedDuration_h_M', NaN, ...
        'ObservedGatedAUC_perCell_M', NaN, ...
        'GatedOutputRate_M', NaN, 'GatedAUC_2hEquivalent_M', NaN, ...
        'DutyCycle_M', NaN, 'ActiveAmplitude_M', NaN, ...
        'InactiveAmplitude_M', NaN, 'ActiveContrast_M', NaN, ...
        'BaselineCorrectedGatedRate_M', NaN, ...
        'RawMeanIntensity_M', NaN, ...
        'WithinCellSD_GatedOutputRate_M', NaN, ...
        'WithinCellSD_DutyCycle_M', NaN, ...
        'WithinCellSD_ActiveAmplitude_M', NaN, ...
        'DecompositionError_M', NaN, ...
        'MeanObservedDuration_h_E', NaN, ...
        'ObservedGatedAUC_perCell_E', NaN, ...
        'GatedOutputRate_E', NaN, 'GatedAUC_2hEquivalent_E', NaN, ...
        'DutyCycle_E', NaN, 'ActiveAmplitude_E', NaN, ...
        'InactiveAmplitude_E', NaN, 'ActiveContrast_E', NaN, ...
        'BaselineCorrectedGatedRate_E', NaN, ...
        'RawMeanIntensity_E', NaN, ...
        'WithinCellSD_GatedOutputRate_E', NaN, ...
        'WithinCellSD_DutyCycle_E', NaN, ...
        'WithinCellSD_ActiveAmplitude_E', NaN, ...
        'DecompositionError_E', NaN, 'ValidPrimaryPair', false);
end

function [M_byCT, E_byCT] = groupPairedMetricByCT( ...
        flyTable, ctList, Mvariable, Evariable, validRows)
    M_byCT = cell(1, numel(ctList));
    E_byCT = cell(1, numel(ctList));
    for j = 1:numel(ctList)
        rows = validRows & ...
            flyTable.AnalysisTime_midpoint_h == ctList(j);
        M_byCT{j} = flyTable.(Mvariable)(rows);
        E_byCT{j} = flyTable.(Evariable)(rows);
        if numel(M_byCT{j}) ~= numel(E_byCT{j})
            error('fig2cIntensity:PairLengthMismatch', ...
                'M/E fly counts differ at midpoint CT %.4g.', ctList(j));
        end
    end
end

function ctTable = buildCTSummaryTable(flyTable, ctList, validRows, ...
        recordingWindowHours)
    nCT = numel(ctList);
    nPairs = zeros(nCT, 1);
    ctTable = table(ctList(:) - recordingWindowHours / 2, ...
        ctList(:) + recordingWindowHours / 2, ctList(:), nPairs, ...
        'VariableNames', {'RecordingStartCT_h', 'RecordingEndCT_h', ...
        'AnalysisTime_midpoint_h', 'N_pairedFlies'});
    metrics = {'GatedAUC_2hEquivalent', 'GatedOutputRate', ...
        'ObservedGatedAUC_perCell', 'DutyCycle', 'ActiveAmplitude', ...
        'InactiveAmplitude', 'ActiveContrast', ...
        'BaselineCorrectedGatedRate', 'RawMeanIntensity', ...
        'MeanObservedDuration_h'};
    for j = 1:nCT
        rows = validRows & flyTable.AnalysisTime_midpoint_h == ctList(j);
        ctTable.N_pairedFlies(j) = sum(rows);
    end
    for m = 1:numel(metrics)
        base = metrics{m};
        for condition = {'M', 'E'}
            c = condition{1};
            variable = [base '_' c];
            meanValues = NaN(nCT, 1);
            sdValues = NaN(nCT, 1);
            semValues = NaN(nCT, 1);
            nValues = zeros(nCT, 1);
            for j = 1:nCT
                rows = validRows & ...
                    flyTable.AnalysisTime_midpoint_h == ctList(j);
                values = flyTable.(variable)(rows);
                finiteValues = values(isfinite(values));
                nValues(j) = numel(finiteValues);
                if ~isempty(finiteValues)
                    meanValues(j) = mean(finiteValues);
                    sdValues(j) = std(finiteValues);
                    semValues(j) = sdValues(j) / sqrt(nValues(j));
                end
            end
            ctTable.(['Mean_' base '_' c]) = meanValues;
            ctTable.(['SD_' base '_' c]) = sdValues;
            ctTable.(['SEM_' base '_' c]) = semValues;
            ctTable.(['Nvalid_' base '_' c]) = nValues;
        end
    end
end

function differenceTable = computeCTDifferenceTable( ...
        ctList, M_byCT, E_byCT, nBootstrap)
    nCT = numel(ctList);
    nPairs = zeros(nCT, 1);
    meanDifference = NaN(nCT, 1);
    sdDifference = NaN(nCT, 1);
    semDifference = NaN(nCT, 1);
    ciLow = NaN(nCT, 1);
    ciHigh = NaN(nCT, 1);
    for j = 1:nCT
        Mj = M_byCT{j}(:);
        Ej = E_byCT{j}(:);
        valid = isfinite(Mj) & isfinite(Ej);
        differences = Mj(valid) - Ej(valid);
        nPairs(j) = numel(differences);
        if isempty(differences)
            continue;
        end
        meanDifference(j) = mean(differences);
        sdDifference(j) = std(differences);
        semDifference(j) = sdDifference(j) / sqrt(nPairs(j));
        bootstrapIndex = randi(nPairs(j), nPairs(j), nBootstrap);
        bootstrapMeans = mean(differences(bootstrapIndex), 1);
        ci = prctile(bootstrapMeans, [2.5 97.5]);
        ciLow(j) = ci(1);
        ciHigh(j) = ci(2);
    end
    differenceTable = table(ctList(:), nPairs, meanDifference, ...
        sdDifference, semDifference, ciLow, ciHigh, ...
        'VariableNames', {'AnalysisTime_midpoint_h', 'N_pairedFlies', ...
        'MeanDifference_MminusE', 'SD_pairedDifference', ...
        'SEM_pairedDifference', 'BootstrapCI95_low', ...
        'BootstrapCI95_high'});
end

function [summaryTable, testsTable] = computePanelCStatistics( ...
        ctList, M_byCT, E_byCT, nBootstrap, equivalenceMargin, alpha)
% Complete-pair CT summaries for Panel C. Biological SEM uses sqrt(nFly).
% Percentile CIs resample paired flies with one common M/E index per CT.
% TOST is an amplitude equivalence analysis and is distinct from the
% interval-null minimum-effect test used for peak timing.
    nCT = numel(ctList);
    nPairs = zeros(nCT, 1);
    meanM = NaN(nCT, 1);
    sdM = NaN(nCT, 1);
    semM = NaN(nCT, 1);
    ciMlow = NaN(nCT, 1);
    ciMhigh = NaN(nCT, 1);
    meanE = NaN(nCT, 1);
    sdE = NaN(nCT, 1);
    semE = NaN(nCT, 1);
    ciElow = NaN(nCT, 1);
    ciEhigh = NaN(nCT, 1);
    meanDifference = NaN(nCT, 1);
    sdDifference = NaN(nCT, 1);
    semDifference = NaN(nCT, 1);
    ciDifferenceLow = NaN(nCT, 1);
    ciDifferenceHigh = NaN(nCT, 1);
    degreesFreedom = NaN(nCT, 1);
    pairedT = NaN(nCT, 1);
    pairedTP = NaN(nCT, 1);
    signrankP = NaN(nCT, 1);
    tostTLower = NaN(nCT, 1);
    tostPLower = NaN(nCT, 1);
    tostTUpper = NaN(nCT, 1);
    tostPUpper = NaN(nCT, 1);
    tostP = NaN(nCT, 1);
    ci90Low = NaN(nCT, 1);
    ci90High = NaN(nCT, 1);

    for j = 1:nCT
        Mj = M_byCT{j}(:);
        Ej = E_byCT{j}(:);
        if numel(Mj) ~= numel(Ej)
            error('fig2cIntensity:PanelCPairLengthMismatch', ...
                'Panel C M/E lengths differ at midpoint CT %.4g.', ...
                ctList(j));
        end
        valid = isfinite(Mj) & isfinite(Ej);
        Mj = Mj(valid);
        Ej = Ej(valid);
        differences = Mj - Ej;
        n = numel(differences);
        nPairs(j) = n;
        if n == 0
            continue;
        end

        meanM(j) = mean(Mj);
        sdM(j) = std(Mj);
        semM(j) = sdM(j) / sqrt(n);
        meanE(j) = mean(Ej);
        sdE(j) = std(Ej);
        semE(j) = sdE(j) / sqrt(n);
        meanDifference(j) = mean(differences);
        sdDifference(j) = std(differences);
        semDifference(j) = sdDifference(j) / sqrt(n);

        bootstrapIndex = randi(n, n, nBootstrap);
        bootstrapMeanM = mean(Mj(bootstrapIndex), 1);
        bootstrapMeanE = mean(Ej(bootstrapIndex), 1);
        bootstrapMeanDifference = bootstrapMeanM - bootstrapMeanE;
        ci = prctile(bootstrapMeanM, [2.5, 97.5]);
        ciMlow(j) = ci(1);
        ciMhigh(j) = ci(2);
        ci = prctile(bootstrapMeanE, [2.5, 97.5]);
        ciElow(j) = ci(1);
        ciEhigh(j) = ci(2);
        ci = prctile(bootstrapMeanDifference, [2.5, 97.5]);
        ciDifferenceLow(j) = ci(1);
        ciDifferenceHigh(j) = ci(2);

        if n >= 2 && isfinite(semDifference(j)) && ...
                semDifference(j) > 0
            degreesFreedom(j) = n - 1;
            pairedT(j) = meanDifference(j) / semDifference(j);
            pairedTP(j) = 2 * tcdf(-abs(pairedT(j)), ...
                degreesFreedom(j));
            signrankP(j) = safeSignrank(differences);

            % Paired two one-sided tests (TOST) for amplitude equivalence:
            % H0: Delta <= -margin OR Delta >= +margin
            % H1: -margin < Delta < +margin
            tostTLower(j) = (meanDifference(j) + equivalenceMargin) / ...
                semDifference(j);
            tostPLower(j) = 1 - tcdf(tostTLower(j), ...
                degreesFreedom(j));
            tostTUpper(j) = (meanDifference(j) - equivalenceMargin) / ...
                semDifference(j);
            tostPUpper(j) = tcdf(tostTUpper(j), degreesFreedom(j));
            tostP(j) = max(tostPLower(j), tostPUpper(j));
            criticalT = tinv(1 - alpha, degreesFreedom(j));
            ci90Low(j) = meanDifference(j) - ...
                criticalT * semDifference(j);
            ci90High(j) = meanDifference(j) + ...
                criticalT * semDifference(j);
        elseif n >= 2
            signrankP(j) = safeSignrank(differences);
        end
    end

    pairedTPHolm = holmAdjustPValues(pairedTP);
    signrankPHolm = holmAdjustPValues(signrankP);
    tostPHolm = holmAdjustPValues(tostP);
    equivalentUnadjusted = tostP < alpha & ...
        ci90Low > -equivalenceMargin & ci90High < equivalenceMargin;
    equivalentHolm = tostPHolm < alpha;

    summaryTable = table(ctList(:), nPairs, meanM, sdM, semM, ...
        ciMlow, ciMhigh, meanE, sdE, semE, ciElow, ciEhigh, ...
        meanDifference, sdDifference, semDifference, ...
        ciDifferenceLow, ciDifferenceHigh, ...
        'VariableNames', {'CT_midpoint_h', 'N_pairedFlies', ...
        'M_mean', 'M_SD', 'M_SEM', 'M_BootstrapCI95_low', ...
        'M_BootstrapCI95_high', 'E_mean', 'E_SD', 'E_SEM', ...
        'E_BootstrapCI95_low', 'E_BootstrapCI95_high', ...
        'Difference_mean_MminusE', 'Difference_SD', ...
        'Difference_SEM', 'Difference_BootstrapCI95_low', ...
        'Difference_BootstrapCI95_high'});

    testsTable = table(ctList(:), nPairs, meanDifference, ...
        sdDifference, semDifference, degreesFreedom, pairedT, ...
        pairedTP, pairedTPHolm, signrankP, signrankPHolm, ...
        repmat(equivalenceMargin, nCT, 1), tostTLower, tostPLower, ...
        tostTUpper, tostPUpper, tostP, tostPHolm, ci90Low, ci90High, ...
        equivalentUnadjusted, equivalentHolm, ...
        'VariableNames', {'CT_midpoint_h', 'N_pairedFlies', ...
        'MeanDifference_MminusE', 'SD_difference', 'SEM_difference', ...
        'DF', 'PairedT_t', 'PairedT_p', 'PairedT_p_Holm', ...
        'Signrank_p_sensitivity', 'Signrank_p_Holm', ...
        'EquivalenceMargin_DeltaFoverF', 'TOST_t_lower', ...
        'TOST_p_lower', 'TOST_t_upper', 'TOST_p_upper', 'TOST_p', ...
        'TOST_p_Holm', 'CI90_low', 'CI90_high', ...
        'Equivalent_unadjusted', 'Equivalent_Holm'});
end

function pValue = safeSignrank(differences)
% Wilcoxon signed-rank is a median-oriented sensitivity analysis; it does
% not replace the paired mean-difference TOST.
    try
        if all(differences == 0)
            pValue = 1;
        else
            pValue = signrank(differences, 0);
        end
    catch
        pValue = NaN;
    end
end

function adjustedP = holmAdjustPValues(rawP)
% Strong family-wise-error control across the finite CT-specific p-values.
    adjustedP = NaN(size(rawP));
    finiteIndex = find(isfinite(rawP));
    if isempty(finiteIndex)
        return;
    end
    [sortedP, order] = sort(rawP(finiteIndex));
    m = numel(sortedP);
    scaledP = (m - (1:m).' + 1) .* sortedP(:);
    adjustedSorted = min(1, cummax(scaledP));
    adjustedP(finiteIndex(order)) = adjustedSorted;
end

function [omnibusTable, pairwiseTable] = ...
        computeWithinConditionCTTests(ctList, M_byCT, E_byCT, ...
        equivalenceMargin, alpha)
% Primary Panel C inference: temporal amplitude modulation within M and
% within E. CT groups contain different flies, so CT-to-CT contrasts are
% independent-sample Welch tests, not paired tests.
%
% Omnibus p-values answer whether a CT effect is detected. They cannot
% prove stability when nonsignificant. Pairwise Welch TOST instead asks
% whether every CT contrast lies inside a prespecified meaningful margin.
    conditions = ["M"; "E"];
    conditionData = {M_byCT, E_byCT};
    nCondition = numel(conditions);
    nCT = numel(ctList);
    nTotal = zeros(nCondition, 1);
    observedRange = NaN(nCondition, 1);
    anovaF = NaN(nCondition, 1);
    anovaDF1 = NaN(nCondition, 1);
    anovaDF2 = NaN(nCondition, 1);
    anovaP = NaN(nCondition, 1);
    etaSquared = NaN(nCondition, 1);
    omegaSquared = NaN(nCondition, 1);
    welchF = NaN(nCondition, 1);
    welchDF1 = NaN(nCondition, 1);
    welchDF2 = NaN(nCondition, 1);
    welchP = NaN(nCondition, 1);
    kruskalWallisP = NaN(nCondition, 1);
    brownForsytheP = NaN(nCondition, 1);

    pairIndex = nchoosek(1:nCT, 2);
    nPairwisePerCondition = size(pairIndex, 1);
    nPairwiseRows = nCondition * nPairwisePerCondition;
    pairCondition = strings(nPairwiseRows, 1);
    CT1 = NaN(nPairwiseRows, 1);
    CT2 = NaN(nPairwiseRows, 1);
    n1 = zeros(nPairwiseRows, 1);
    n2 = zeros(nPairwiseRows, 1);
    mean1 = NaN(nPairwiseRows, 1);
    mean2 = NaN(nPairwiseRows, 1);
    difference = NaN(nPairwiseRows, 1);
    welchSE = NaN(nPairwiseRows, 1);
    welchDF = NaN(nPairwiseRows, 1);
    welchT = NaN(nPairwiseRows, 1);
    differenceP = NaN(nPairwiseRows, 1);
    ci90Low = NaN(nPairwiseRows, 1);
    ci90High = NaN(nPairwiseRows, 1);
    tostTLower = NaN(nPairwiseRows, 1);
    tostPLower = NaN(nPairwiseRows, 1);
    tostTUpper = NaN(nPairwiseRows, 1);
    tostPUpper = NaN(nPairwiseRows, 1);
    tostP = NaN(nPairwiseRows, 1);

    row = 0;
    for c = 1:nCondition
        dataByCT = conditionData{c};
        cleaned = cell(1, nCT);
        groupMeans = NaN(1, nCT);
        groupVariances = NaN(1, nCT);
        groupCounts = zeros(1, nCT);
        allValues = [];
        allGroups = [];
        for j = 1:nCT
            values = dataByCT{j}(:);
            values = values(isfinite(values));
            cleaned{j} = values;
            groupCounts(j) = numel(values);
            if ~isempty(values)
                groupMeans(j) = mean(values);
                groupVariances(j) = var(values);
                allValues = [allValues; values]; %#ok<AGROW>
                allGroups = [allGroups; repmat(j, numel(values), 1)]; %#ok<AGROW>
            end
        end
        nTotal(c) = numel(allValues);
        validGroups = groupCounts > 0 & isfinite(groupMeans);
        k = sum(validGroups);
        if k >= 2
            observedRange(c) = max(groupMeans(validGroups)) - ...
                min(groupMeans(validGroups));
        end

        % Conventional equal-variance one-way ANOVA, calculated directly
        % so all sums of squares and effect sizes remain auditable.
        if k >= 2 && nTotal(c) > k
            grandMean = mean(allValues);
            ssBetween = sum(groupCounts(validGroups) .* ...
                (groupMeans(validGroups) - grandMean).^2);
            ssWithin = 0;
            for j = find(validGroups)
                ssWithin = ssWithin + sum((cleaned{j} - ...
                    groupMeans(j)).^2);
            end
            ssTotal = ssBetween + ssWithin;
            anovaDF1(c) = k - 1;
            anovaDF2(c) = nTotal(c) - k;
            msBetween = ssBetween / anovaDF1(c);
            msWithin = ssWithin / anovaDF2(c);
            if msWithin > 0
                anovaF(c) = msBetween / msWithin;
                anovaP(c) = 1 - fcdf(anovaF(c), ...
                    anovaDF1(c), anovaDF2(c));
            elseif msBetween == 0
                anovaF(c) = 0;
                anovaP(c) = 1;
            else
                anovaF(c) = Inf;
                anovaP(c) = 0;
            end
            if ssTotal > 0
                etaSquared(c) = ssBetween / ssTotal;
                omegaSquared(c) = max(0, ...
                    (ssBetween - anovaDF1(c) * msWithin) / ...
                    (ssTotal + msWithin));
            end
        end

        % Welch's one-way ANOVA is robust to unequal CT-group variances.
        welchValid = groupCounts >= 2 & isfinite(groupVariances) & ...
            groupVariances > 0;
        kWelch = sum(welchValid);
        if kWelch >= 2
            meansW = groupMeans(welchValid);
            variancesW = groupVariances(welchValid);
            countsW = groupCounts(welchValid);
            weights = countsW ./ variancesW;
            weightSum = sum(weights);
            weightedMean = sum(weights .* meansW) / weightSum;
            correctionTerms = ((1 - weights / weightSum).^2) ./ ...
                (countsW - 1);
            correctionSum = sum(correctionTerms);
            welchDF1(c) = kWelch - 1;
            welchDF2(c) = (kWelch^2 - 1) / (3 * correctionSum);
            numerator = sum(weights .* (meansW - weightedMean).^2) / ...
                welchDF1(c);
            denominator = 1 + 2 * (kWelch - 2) / ...
                (kWelch^2 - 1) * correctionSum;
            welchF(c) = numerator / denominator;
            welchP(c) = 1 - fcdf(welchF(c), ...
                welchDF1(c), welchDF2(c));
        end

        try
            kruskalWallisP(c) = kruskalwallis( ...
                allValues, allGroups, 'off');
        catch
            kruskalWallisP(c) = NaN;
        end
        try
            brownForsytheP(c) = vartestn(allValues, allGroups, ...
                'TestType', 'BrownForsythe', 'Display', 'off');
        catch
            brownForsytheP(c) = NaN;
        end

        for q = 1:nPairwisePerCondition
            row = row + 1;
            j1 = pairIndex(q, 1);
            j2 = pairIndex(q, 2);
            values1 = cleaned{j1};
            values2 = cleaned{j2};
            pairCondition(row) = conditions(c);
            CT1(row) = ctList(j1);
            CT2(row) = ctList(j2);
            n1(row) = numel(values1);
            n2(row) = numel(values2);
            if n1(row) < 2 || n2(row) < 2
                continue;
            end
            mean1(row) = mean(values1);
            mean2(row) = mean(values2);
            difference(row) = mean1(row) - mean2(row);
            variance1 = var(values1);
            variance2 = var(values2);
            varianceTerm1 = variance1 / n1(row);
            varianceTerm2 = variance2 / n2(row);
            welchSE(row) = sqrt(varianceTerm1 + varianceTerm2);
            denominatorDF = varianceTerm1^2 / (n1(row) - 1) + ...
                varianceTerm2^2 / (n2(row) - 1);
            if welchSE(row) <= 0 || denominatorDF <= 0
                continue;
            end
            welchDF(row) = (varianceTerm1 + varianceTerm2)^2 / ...
                denominatorDF;
            welchT(row) = difference(row) / welchSE(row);
            differenceP(row) = 2 * tcdf(-abs(welchT(row)), ...
                welchDF(row));
            criticalT = tinv(1 - alpha, welchDF(row));
            ci90Low(row) = difference(row) - criticalT * welchSE(row);
            ci90High(row) = difference(row) + criticalT * welchSE(row);
            tostTLower(row) = (difference(row) + equivalenceMargin) / ...
                welchSE(row);
            tostPLower(row) = 1 - tcdf(tostTLower(row), welchDF(row));
            tostTUpper(row) = (difference(row) - equivalenceMargin) / ...
                welchSE(row);
            tostPUpper(row) = tcdf(tostTUpper(row), welchDF(row));
            tostP(row) = max(tostPLower(row), tostPUpper(row));
        end
    end

    differencePHolm = NaN(size(differenceP));
    tostPHolm = NaN(size(tostP));
    equivalentUnadjusted = false(size(tostP));
    equivalentHolm = false(size(tostP));
    nEquivalentUnadjusted = zeros(nCondition, 1);
    nEquivalentHolm = zeros(nCondition, 1);
    allEquivalentIUT = false(nCondition, 1);
    allEquivalentHolm = false(nCondition, 1);
    for c = 1:nCondition
        rows = pairCondition == conditions(c);
        differencePHolm(rows) = holmAdjustPValues(differenceP(rows));
        tostPHolm(rows) = holmAdjustPValues(tostP(rows));
        equivalentUnadjusted(rows) = tostP(rows) < alpha & ...
            ci90Low(rows) > -equivalenceMargin & ...
            ci90High(rows) < equivalenceMargin;
        equivalentHolm(rows) = tostPHolm(rows) < alpha;
        nEquivalentUnadjusted(c) = sum(equivalentUnadjusted(rows));
        nEquivalentHolm(c) = sum(equivalentHolm(rows));
        % Global all-pairs equivalence is an intersection-union test: every
        % component TOST must reject at alpha. Holm is also shown for the
        % stronger goal of simultaneous pair-specific claims.
        allEquivalentIUT(c) = all(equivalentUnadjusted(rows));
        allEquivalentHolm(c) = all(equivalentHolm(rows));
    end

    omnibusTable = table(conditions, repmat(nCT, nCondition, 1), ...
        nTotal, observedRange, anovaF, anovaDF1, anovaDF2, anovaP, ...
        etaSquared, omegaSquared, welchF, welchDF1, welchDF2, welchP, ...
        repmat(nCT - 1, nCondition, 1), kruskalWallisP, ...
        brownForsytheP, repmat(equivalenceMargin, nCondition, 1), ...
        repmat(nPairwisePerCondition, nCondition, 1), ...
        nEquivalentUnadjusted, nEquivalentHolm, allEquivalentIUT, ...
        allEquivalentHolm, ...
        'VariableNames', {'Condition', 'N_CT', 'N_total', ...
        'ObservedMaxMinusMinMean', 'ANOVA_F', 'ANOVA_DF1', ...
        'ANOVA_DF2', 'ANOVA_p', 'ANOVA_etaSquared', ...
        'ANOVA_omegaSquared', 'WelchANOVA_F', 'WelchANOVA_DF1', ...
        'WelchANOVA_DF2', 'WelchANOVA_p', 'KruskalWallis_DF', ...
        'KruskalWallis_p', 'BrownForsythe_p', ...
        'EquivalenceMargin_DeltaFoverF', 'N_pairwiseComparisons', ...
        'N_pairwiseEquivalent_unadjusted', ...
        'N_pairwiseEquivalent_Holm', 'AllPairwiseEquivalent_IUT', ...
        'AllPairwiseEquivalent_Holm'});

    pairwiseTable = table(pairCondition, CT1, CT2, n1, n2, mean1, ...
        mean2, difference, welchSE, welchDF, welchT, differenceP, ...
        differencePHolm, repmat(equivalenceMargin, nPairwiseRows, 1), ...
        ci90Low, ci90High, tostTLower, tostPLower, tostTUpper, ...
        tostPUpper, tostP, tostPHolm, equivalentUnadjusted, ...
        equivalentHolm, ...
        'VariableNames', {'Condition', 'CT1_h', 'CT2_h', 'N_CT1', ...
        'N_CT2', 'Mean_CT1', 'Mean_CT2', ...
        'Difference_CT1minusCT2', 'Welch_SE', 'Welch_DF', ...
        'Welch_t', 'Difference_p', 'Difference_p_Holm', ...
        'EquivalenceMargin_DeltaFoverF', 'CI90_low', 'CI90_high', ...
        'TOST_t_lower', 'TOST_p_lower', 'TOST_t_upper', ...
        'TOST_p_upper', 'TOST_p', 'TOST_p_Holm', ...
        'Equivalent_unadjusted_IUTcomponent', ...
        'Equivalent_Holm_simultaneous'});
end

function [anovaTable, coefficientTable, infoTable] = ...
        fitPanelCAmplitudeLME(flyTable, validRows, panelCCT)
% CT is between flies and Condition is paired within fly. The random
% intercept preserves the within-fly M/E dependence.
    selected = validRows & ...
        ismember(flyTable.AnalysisTime_midpoint_h, panelCCT) & ...
        isfinite(flyTable.ActiveAmplitude_M) & ...
        isfinite(flyTable.ActiveAmplitude_E);
    rowIndex = find(selected);
    nPair = numel(rowIndex);
    modelFormula = ...
        "Amplitude ~ Condition * CT + (1|FlyID)";
    if nPair < 3
        statusMessage = "Too few complete pairs for mixed model";
        anovaTable = table(statusMessage, 'VariableNames', {'Status'});
        coefficientTable = anovaTable;
        infoTable = table("Failed", statusMessage, modelFormula, ...
            "Satterthwaite", 2 * nPair, nPair, numel(panelCCT), ...
            'VariableNames', {'Status', 'Message', 'ModelFormula', ...
            'DFMethod', 'N_observations', 'N_pairedFlies', 'N_CT'});
        return;
    end

    amplitudeM = flyTable.ActiveAmplitude_M(selected);
    amplitudeE = flyTable.ActiveAmplitude_E(selected);
    ctValues = flyTable.AnalysisTime_midpoint_h(selected);
    pairID = "row" + string(rowIndex);
    amplitude = [amplitudeM; amplitudeE];
    condition = categorical([repmat("M", nPair, 1); ...
        repmat("E", nPair, 1)], ["E", "M"]);
    CT = categorical([string(ctValues); string(ctValues)]);
    FlyID = categorical([pairID; pairID]);
    modelData = table(amplitude, condition, CT, FlyID, ...
        'VariableNames', {'Amplitude', 'Condition', 'CT', 'FlyID'});

    try
        lme = fitlme(modelData, char(modelFormula), 'FitMethod', 'REML');
        anovaResult = anova(lme, 'DFMethod', 'Satterthwaite');
        if istable(anovaResult)
            anovaTable = anovaResult;
        else
            % R2023a returns classreg.regr.lmeutils.titleddataset.
            anovaTable = dataset2table(anovaResult);
        end
        if ~isempty(anovaTable.Properties.RowNames)
            anovaTable = addvars(anovaTable, ...
                string(anovaTable.Properties.RowNames), 'Before', 1, ...
                'NewVariableNames', 'Term');
            anovaTable.Properties.RowNames = {};
        end
        coefficientResult = lme.Coefficients;
        if istable(coefficientResult)
            coefficientTable = coefficientResult;
        else
            coefficientTable = dataset2table(coefficientResult);
        end
        if ~isempty(coefficientTable.Properties.RowNames)
            coefficientTable = addvars(coefficientTable, ...
                string(coefficientTable.Properties.RowNames), 'Before', 1, ...
                'NewVariableNames', 'Coefficient');
            coefficientTable.Properties.RowNames = {};
        end
        infoTable = table("OK", "", modelFormula, "Satterthwaite", ...
            height(modelData), nPair, numel(categories(CT)), ...
            'VariableNames', {'Status', 'Message', 'ModelFormula', ...
            'DFMethod', 'N_observations', 'N_pairedFlies', 'N_CT'});
    catch modelError
        warning('fig2cIntensity:PanelCLMEFailed', ...
            'Panel C mixed model could not be fitted: %s', ...
            modelError.message);
        statusMessage = string(modelError.message);
        anovaTable = table(statusMessage, 'VariableNames', {'Status'});
        coefficientTable = anovaTable;
        infoTable = table("Failed", statusMessage, modelFormula, ...
            "Satterthwaite", height(modelData), nPair, ...
            numel(categories(CT)), 'VariableNames', ...
            {'Status', 'Message', 'ModelFormula', 'DFMethod', ...
            'N_observations', 'N_pairedFlies', 'N_CT'});
    end
end

function modulationTable = computePanelCModulation( ...
        ctList, M_byCT, E_byCT, nBootstrap, amplitudeMargin)
% Descriptive CT-modulation diagnostic. Each bootstrap replicate resamples
% paired flies within every CT stratum. Range is max(CT mean)-min(CT mean).
    nCT = numel(ctList);
    observedMeansM = NaN(1, nCT);
    observedMeansE = NaN(1, nCT);
    bootstrapMeansM = NaN(nBootstrap, nCT);
    bootstrapMeansE = NaN(nBootstrap, nCT);
    nPairs = zeros(1, nCT);
    for j = 1:nCT
        Mj = M_byCT{j}(:);
        Ej = E_byCT{j}(:);
        valid = isfinite(Mj) & isfinite(Ej);
        Mj = Mj(valid);
        Ej = Ej(valid);
        n = numel(Mj);
        nPairs(j) = n;
        if n == 0
            continue;
        end
        observedMeansM(j) = mean(Mj);
        observedMeansE(j) = mean(Ej);
        bootstrapIndex = randi(n, n, nBootstrap);
        bootstrapMeansM(:, j) = mean(Mj(bootstrapIndex), 1).';
        bootstrapMeansE(:, j) = mean(Ej(bootstrapIndex), 1).';
    end

    conditions = ["M"; "E"];
    observedMeans = {observedMeansM, observedMeansE};
    bootstrapMeans = {bootstrapMeansM, bootstrapMeansE};
    observedRange = NaN(2, 1);
    rangeCILow = NaN(2, 1);
    rangeCIHigh = NaN(2, 1);
    observedModulationIndex = NaN(2, 1);
    modulationIndexCILow = NaN(2, 1);
    modulationIndexCIHigh = NaN(2, 1);
    CT_atMinimumMean = NaN(2, 1);
    CT_atMaximumMean = NaN(2, 1);
    bootstrapSupportRangeBelowMargin = NaN(2, 1);
    for c = 1:2
        meansObserved = observedMeans{c};
        meansBootstrap = bootstrapMeans{c};
        if all(isfinite(meansObserved))
            [minimumMean, minimumIndex] = min(meansObserved);
            [maximumMean, maximumIndex] = max(meansObserved);
            observedRange(c) = maximumMean - minimumMean;
            CT_atMinimumMean(c) = ctList(minimumIndex);
            CT_atMaximumMean(c) = ctList(maximumIndex);
            denominator = abs(mean(meansObserved));
            if denominator > 0
                observedModulationIndex(c) = observedRange(c) / denominator;
            end
        end
        validBootstrap = all(isfinite(meansBootstrap), 2);
        if ~any(validBootstrap)
            continue;
        end
        validMeans = meansBootstrap(validBootstrap, :);
        bootstrapRange = max(validMeans, [], 2) - min(validMeans, [], 2);
        bootstrapDenominator = abs(mean(validMeans, 2));
        bootstrapModulationIndex = bootstrapRange ./ bootstrapDenominator;
        bootstrapModulationIndex(bootstrapDenominator <= 0) = NaN;
        ci = prctile(bootstrapRange, [2.5, 97.5]);
        rangeCILow(c) = ci(1);
        rangeCIHigh(c) = ci(2);
        ci = prctile(bootstrapModulationIndex(isfinite( ...
            bootstrapModulationIndex)), [2.5, 97.5]);
        if ~isempty(ci)
            modulationIndexCILow(c) = ci(1);
            modulationIndexCIHigh(c) = ci(2);
        end
        bootstrapSupportRangeBelowMargin(c) = ...
            mean(bootstrapRange < amplitudeMargin);
    end

    modulationTable = table(conditions, repmat(nCT, 2, 1), ...
        repmat(sum(nPairs), 2, 1), CT_atMinimumMean, CT_atMaximumMean, ...
        observedRange, rangeCILow, rangeCIHigh, ...
        observedModulationIndex, modulationIndexCILow, ...
        modulationIndexCIHigh, repmat(amplitudeMargin, 2, 1), ...
        rangeCIHigh < amplitudeMargin, ...
        bootstrapSupportRangeBelowMargin, ...
        'VariableNames', {'Condition', 'N_CT', ...
        'N_pairedFlyRecords_acrossCT', 'CT_atMinimumMean', ...
        'CT_atMaximumMean', 'ObservedMaxMinusMinMean', ...
        'BootstrapCI95_low', 'BootstrapCI95_high', ...
        'ObservedModulationIndex', 'ModulationIndexCI95_low', ...
        'ModulationIndexCI95_high', 'AmplitudeMargin_DeltaFoverF', ...
        'RangeUpperCI_belowMargin', ...
        'BootstrapSupport_rangeBelowMargin_descriptive'});
end

function plotMetricByCT(ax, ctList, M_byCT, E_byCT, Mcolor, Ecolor, ...
        yLabelText, titleText)
    [meanM, semM] = summarizeCellArray(M_byCT);
    [meanE, semE] = summarizeCellArray(E_byCT);
    hold(ax, 'on');
    plotJitteredCTValues(ax, ctList, E_byCT, Ecolor, 0.12);
    plotJitteredCTValues(ax, ctList, M_byCT, Mcolor, -0.12);
    errorbar(ax, ctList, meanE, semE, '-o', 'Color', Ecolor, ...
        'MarkerFaceColor', Ecolor, 'MarkerSize', 3, 'LineWidth', 1.0, ...
        'DisplayName', 'Ecell mean +/- SEM');
    errorbar(ax, ctList, meanM, semM, '-o', 'Color', Mcolor, ...
        'MarkerFaceColor', Mcolor, 'MarkerSize', 3, 'LineWidth', 1.0, ...
        'DisplayName', 'Mcell mean +/- SEM');
    xlabel(ax, 'Recording-window midpoint CT (h)');
    ylabel(ax, yLabelText);
    title(ax, titleText);
    legend(ax, 'Location', 'best');
    styleAxes(ax, ctList);
end

function plotJitteredCTValues(ax, ctList, valuesByCT, color, offset)
    for j = 1:numel(ctList)
        values = valuesByCT{j};
        values = values(isfinite(values));
        n = numel(values);
        if n == 0
            continue;
        end
        if n == 1
            jitter = 0;
        else
            jitter = linspace(-0.055, 0.055, n).';
        end
        scatter(ax, ctList(j) + offset + jitter, values(:), 10, color, ...
            'filled', 'MarkerFaceAlpha', 0.35, ...
            'HandleVisibility', 'off');
    end
end

function [meanValues, semValues] = summarizeCellArray(valuesByCT)
    meanValues = NaN(1, numel(valuesByCT));
    semValues = NaN(1, numel(valuesByCT));
    for j = 1:numel(valuesByCT)
        values = valuesByCT{j};
        values = values(isfinite(values));
        if isempty(values)
            continue;
        end
        meanValues(j) = mean(values);
        semValues(j) = std(values) / sqrt(numel(values));
    end
end

function plotPairedDifferenceByCT(ax, ctList, M_byCT, E_byCT, ...
        differenceTable, color)
    hold(ax, 'on');
    for j = 1:numel(ctList)
        Mj = M_byCT{j}(:);
        Ej = E_byCT{j}(:);
        valid = isfinite(Mj) & isfinite(Ej);
        differences = Mj(valid) - Ej(valid);
        n = numel(differences);
        if n == 0
            continue;
        end
        if n == 1
            jitter = 0;
        else
            jitter = linspace(-0.08, 0.08, n).';
        end
        scatter(ax, ctList(j) + jitter, differences, 12, color, ...
            'filled', 'MarkerFaceAlpha', 0.40, ...
            'HandleVisibility', 'off');
    end
    lowerError = differenceTable.MeanDifference_MminusE - ...
        differenceTable.BootstrapCI95_low;
    upperError = differenceTable.BootstrapCI95_high - ...
        differenceTable.MeanDifference_MminusE;
    errorbar(ax, ctList, differenceTable.MeanDifference_MminusE, ...
        lowerError, upperError, '-o', 'Color', color, ...
        'MarkerFaceColor', color, 'MarkerSize', 3, 'LineWidth', 1.0, ...
        'DisplayName', 'Mean paired difference, 95% bootstrap CI');
    yline(ax, 0, '--k', 'LineWidth', 0.8, ...
        'HandleVisibility', 'off');
    xlabel(ax, 'Recording-window midpoint CT (h)');
    ylabel(ax, 'Gated AUC difference M-E (DeltaF/F*h)');
    title(ax, 'Paired fly-level output difference');
    legend(ax, 'Location', 'best');
    styleAxes(ax, ctList);
end

function plotOuterDifference(ax, values, binWidth, observedDiff, ...
        delta0, MET_z, MET_p)
    edges = makeCenteredEdges(values, binWidth);
    histogram(ax, values, edges, 'Normalization', 'probability', ...
        'FaceColor', [0.45 0.45 0.45], 'EdgeColor', 'none', ...
        'FaceAlpha', 0.8);
    hold(ax, 'on');
    if numel(unique(values)) > 1
        [density, xDensity] = ksdensity(values);
        plot(ax, xDensity, density * binWidth, 'Color', [0.15 0.15 0.15], ...
            'LineWidth', 1.0, 'DisplayName', 'KDE x bin width');
    end
    xline(ax, -delta0, '--k', '-delta0', 'LineWidth', 0.9);
    xline(ax, delta0, '--k', '+delta0', 'LineWidth', 0.9);
    xline(ax, observedDiff, '-', 'Observed diff', 'Color', [0.75 0 0], ...
        'LineWidth', 1.2);
    xlabel(ax, 'Outer bagged difference M-E (h)');
    ylabel(ax, 'Bootstrap proportion per bin');
    title(ax, sprintf('Outer difference: MET Z=%.2f, p=%.3g', ...
        MET_z, MET_p));
    styleAxes(ax, []);
end

function plotGroupedValues(ax, groups, labels, colors)
    hold(ax, 'on');
    for g = 1:numel(groups)
        values = groups{g};
        values = values(isfinite(values));
        n = numel(values);
        if n == 0
            continue;
        end
        if n == 1
            jitter = 0;
        else
            jitter = linspace(-0.12, 0.12, n).';
        end
        scatter(ax, g + jitter, values(:), 12, colors{g}, 'filled', ...
            'MarkerFaceAlpha', 0.40, 'HandleVisibility', 'off');
        groupMean = mean(values);
        groupSEM = std(values) / sqrt(n);
        errorbar(ax, g, groupMean, groupSEM, 'o', 'Color', colors{g}, ...
            'MarkerFaceColor', colors{g}, 'MarkerSize', 4, ...
            'LineWidth', 1.1, 'HandleVisibility', 'off');
    end
    xlim(ax, [0.5, numel(groups) + 0.5]);
    xticks(ax, 1:numel(groups));
    xticklabels(ax, labels);
    xtickangle(ax, 20);
    styleAxes(ax, []);
end

function rho = safeSpearman(x, y)
    valid = isfinite(x) & isfinite(y);
    if sum(valid) < 3 || numel(unique(x(valid))) < 2 || ...
            numel(unique(y(valid))) < 2
        rho = NaN;
        return;
    end
    rho = corr(x(valid), y(valid), 'Type', 'Spearman');
end

function panelCFig = plotPanelCAmplitudeByCondition(panelCCT, ...
        M_byCT, E_byCT, panelCCTStats, Mcolor, Ecolor, outputSVG, useBoxplot)
% Separate condition panels match the within-M / within-E omnibus question.
% Complete-pair fly selection and all saved numerical results are retained.
% CT is vertical (earliest at bottom); amplitude is horizontal.
% Default display: fly points with mean and paired-bootstrap 95% CI.
% Optional boxplot display: median and Q1-Q3, with whiskers extending to the
% most extreme values within 1.5 IQR of the box. All flies, including
% outliers, are overlaid once as jittered points in either display mode.
    panelCFig = figure('Name', ...
        'Panel C: active-state amplitude within M and within E', ...
        'NumberTitle', 'off', 'Color', 'w', 'Units', 'centimeters', ...
        'Position', [2, 2, 4.9, 4]);
    panelCAxes = gobjects(1, 2);
    panelCAxes(1) = axes(panelCFig, 'Units', 'centimeters', ...
        'Position', [0.82, 0.70, 1.66, 2.60]);
    panelCAxes(2) = axes(panelCFig, 'Units', 'centimeters', ...
        'Position', [2.90, 0.70, 1.66, 2.60]);

    completeM_byCT = cell(size(M_byCT));
    completeE_byCT = cell(size(E_byCT));
    for j = 1:numel(panelCCT)
        Mj = M_byCT{j}(:);
        Ej = E_byCT{j}(:);
        valid = isfinite(Mj) & isfinite(Ej);
        completeM_byCT{j} = Mj(valid);
        completeE_byCT{j} = Ej(valid);
    end
    % Retain the prior common amplitude limits for visual comparability.
    allAmplitude = [vertcat(completeM_byCT{:}); ...
        vertcat(completeE_byCT{:}); ...
        panelCCTStats.M_BootstrapCI95_low; ...
        panelCCTStats.M_BootstrapCI95_high; ...
        panelCCTStats.E_BootstrapCI95_low; ...
        panelCCTStats.E_BootstrapCI95_high];
    allAmplitude = allAmplitude(isfinite(allAmplitude));
    if isempty(allAmplitude)
        amplitudeLimits = [0, 1];
    else
        dataMinimum = min(allAmplitude);
        dataMaximum = max(allAmplitude);
        dataRange = max(dataMaximum - dataMinimum, 0.1);
        amplitudeLimits = [max(0, dataMinimum - 0.12 * dataRange), ...
            dataMaximum + 0.12 * dataRange];
    end

    conditionData = {completeM_byCT, completeE_byCT};
    conditionNames = {'M', 'E'};
    conditionColors = {Mcolor, Ecolor};
    meanValues = {panelCCTStats.M_mean, panelCCTStats.E_mean};
    lowerCI = {panelCCTStats.M_BootstrapCI95_low, ...
        panelCCTStats.E_BootstrapCI95_low};
    upperCI = {panelCCTStats.M_BootstrapCI95_high, ...
        panelCCTStats.E_BootstrapCI95_high};
    for c = 1:2
        ax = panelCAxes(c);
        hold(ax, 'on');
        if useBoxplot
            valuesByCT = conditionData{c};
            boxValues = vertcat(valuesByCT{:});
            boxCT = repelem(panelCCT(:), cellfun(@numel, valuesByCT(:)));
            if ~isempty(boxValues)
                % Native MATLAB boxchart uses the 1.5-IQR whisker convention.
                % Suppress its own outlier symbols: scatter shows all flies.
                boxchart(ax, boxCT, boxValues, 'Orientation', 'horizontal', ...
                    'BoxFaceColor', conditionColors{c}, ...
                    'BoxEdgeColor', conditionColors{c}, ...
                    'BoxMedianLineColor', conditionColors{c}, ...
                    'WhiskerLineColor', conditionColors{c}, ...
                    'BoxFaceAlpha', 0.15, 'BoxWidth', 0.60, ...
                    'LineWidth', 0.8, 'MarkerStyle', 'none', ...
                    'DisplayName', 'Median, IQR and 1.5-IQR whiskers');
            end
        end
        for j = 1:numel(panelCCT)
            values = conditionData{c}{j};
            nFly = numel(values);
            if nFly == 0
                continue;
            elseif nFly == 1
                ctJitter = 0;
            else
                ctJitter = linspace(-0.16, 0.16, nFly).';
            end
            scatter(ax, values, panelCCT(j) + ctJitter, ...
                11, conditionColors{c}, 'filled', ...
                'MarkerFaceAlpha', 0.50, 'MarkerEdgeColor', 'none', ...
                'HandleVisibility', 'off');
        end
        if ~useBoxplot
            errorbar(ax, meanValues{c}, panelCCT(:), ...
                meanValues{c} - lowerCI{c}, upperCI{c} - meanValues{c}, ...
                'horizontal', 'o', 'Color', conditionColors{c}, ...
                'MarkerFaceColor', conditionColors{c}, 'MarkerSize', 3.4, ...
                'CapSize', 3, 'LineWidth', 0.9, ...
                'DisplayName', 'Mean and 95% bootstrap CI');
        end
        xlim(ax, amplitudeLimits);
        ylim(ax, [min(panelCCT) - 1, max(panelCCT) + 1]);
        yticks(ax, panelCCT);
        ax.YDir = 'normal';
        ax.YMinorTick = 'on';
        ax.YAxis.MinorTickValues = setdiff( ...
            (min(panelCCT) - 1):1:(max(panelCCT) + 1), panelCCT);
        xlabel(ax, 'Active \DeltaF/F');
        title(ax, conditionNames{c}, 'Color', conditionColors{c}, ...
            'FontWeight', 'normal');
        styleAxes(ax, []);
    end
    ylabel(panelCAxes(1), 'CT (h)');
    panelCAxes(2).YAxis.Visible = 'off';
    % Link only CT. Each condition has its own amplitude axis; equal
    % initial limits keep magnitude comparisons visually consistent.
    linkaxes(panelCAxes, 'y');

    applyFigureTypography(panelCFig);
    panelCFig.PaperUnits = 'centimeters';
    panelCFig.PaperSize = [4.9, 4];
    panelCFig.PaperPosition = [0, 0, 4.9, 4];
    panelCFig.PaperPositionMode = 'manual';
    panelCFig.InvertHardcopy = 'off';
    drawnow;
    print(panelCFig, outputSVG, '-dsvg', '-vector');
    fprintf('Standalone Panel C saved to: %s\n', outputSVG);
end

function styleAxes(ax, ctList)
    ax.LineWidth = 0.8;
    ax.TickDir = 'out';
    ax.TickLength = [0.02, 0.02];
    ax.FontName = 'Arial';
    ax.FontSize = 7;
    ax.Color = 'none';
    box(ax, 'off');
    grid(ax, 'off');
    if ~isempty(ctList)
        xlim(ax, [min(ctList) - 0.8, max(ctList) + 0.8]);
        xticks(ax, ctList);
    end
end

function applyFigureTypography(fig)
    set(findall(fig, '-property', 'FontUnits'), 'FontUnits', 'points');
    set(findall(fig, '-property', 'FontSize'), 'FontSize', 7);
    set(findall(fig, '-property', 'FontSizeMode'), ...
        'FontSizeMode', 'manual');
    set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');
end

function addPanelLabel(ax, labelText)
    text(ax, -0.13, 1.06, labelText, 'Units', 'normalized', ...
        'FontName', 'Arial', 'FontSize', 8, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', ...
        'Clipping', 'off');
end

function maxY = maxFiniteY(ax)
    limits = ylim(ax);
    maxY = limits(2);
    if ~isfinite(maxY) || maxY <= 0
        maxY = 1;
    end
end
