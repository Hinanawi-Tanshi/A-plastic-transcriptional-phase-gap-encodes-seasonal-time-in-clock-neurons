function publicResult = runFlyPairedKDEAnalysis(cfg)
% Shared publication pipeline; settings are supplied by the figure entry script.
figureCleanup = bootstrapFigureContext(cfg); %#ok<NASGU>
bootstrapRawStore('configure',cfg.inputFile,cfg);
scriptFolder = cfg.outputDir;
inputRootFolders = fullfile(cfg.inputFile, cfg.batchNames);
nInner = cfg.nInner;
nOuter = cfg.nOuter;
delta0 = cfg.delta0;
alpha = cfg.alpha;
rngSeed = cfg.rngSeed;
recordingWindowHours = cfg.recordingWindowHours;
representativeTimeOffset = cfg.representativeTimeOffset;
kdeWindowMinutes = cfg.kdeWindowMinutes;
mainFigureBinWidth = cfg.mainFigureBinWidth;
diagnosticBinWidth = cfg.diagnosticBinWidth;
outputExcel = fullfile(scriptFolder, cfg.outputStem + "_results.xlsx");
outputMainSVG = fullfile(scriptFolder, cfg.outputStem + ".svg");
outputDiagnosticSVG = fullfile(scriptFolder, cfg.outputStem + "_diagnostics.svg");

if isempty(inputRootFolders)
    error('fig2cBagged:NoInputRoots', ...
        'At least one dataset root must be listed in inputRootFolders.');
end
if numel(unique(lower(inputRootFolders))) ~= numel(inputRootFolders)
    error('fig2cBagged:DuplicateInputRoot', ...
        'inputRootFolders contains a duplicate directory.');
end
missingRoot = ~bootstrapRawIsFolder(inputRootFolders);
if any(missingRoot)
    error('fig2cBagged:MissingDataFolder', ...
        'Data folder not found: %s', ...
        strjoin(inputRootFolders(missingRoot), ' | '));
end

computeKDEFile = string(which('compute_KDE'));

% Direct raw-XLSX import reproduces the statistical preprocessing used by
% fig2c.m. No Mcell_KDE.mat, Ecell_KDE.mat, or perSampleStats.mat is read.
% Virtual CT groups are inferred from the first timestamp, rounded to an
% hour. A group CTxx represents [CTxx, CTxx+2]. Because the
% analyzed quantity is the two-hour average KDE, its representative time is
% the interval midpoint CTxx+1, not the interval start CTxx.
[ctList, M_byCT, E_byCT, flyMetricsTable, cellMetricsTable, ...
    sourceRootsTable] = loadPairedFlyKDEMeansFromRaw(inputRootFolders, ...
    recordingWindowHours, kdeWindowMinutes);
nCT = numel(ctList);
if nCT < 2
    error('fig2cBagged:TooFewCT', ...
        'At least two valid CT groups are required.');
end

nPairs = cellfun(@numel, M_byCT).';
meanCurve_M = cellfun(@(v) mean(v, 'omitnan'), M_byCT).';
meanCurve_E = cellfun(@(v) mean(v, 'omitnan'), E_byCT).';
sdCurve_M = cellfun(@(v) std(v, 'omitnan'), M_byCT).';
sdCurve_E = cellfun(@(v) std(v, 'omitnan'), E_byCT).';
semCurve_M = sdCurve_M ./ sqrt(nPairs);
semCurve_E = sdCurve_E ./ sqrt(nPairs);

fprintf('========== fig 2c paired fly-level KDE summary ==========\n');
fprintf(['Raw input: %d roots, %d workbooks, %d sheets; ' ...
    'KDE window = %.4g min\n'], height(sourceRootsTable), ...
    height(flyMetricsTable), height(cellMetricsTable), kdeWindowMinutes);
for r = 1:height(sourceRootsTable)
    fprintf('  root %d: %s (%d CT folders, %d paired flies)\n', ...
        sourceRootsTable.SourceRootIndex(r), ...
        sourceRootsTable.SourceRoot(r), ...
        sourceRootsTable.N_CTFolders(r), ...
        sourceRootsTable.N_validPairedFlies(r));
end
fprintf('%-8s %-8s %-8s %-8s %-13s %-13s\n', ...
    'Start', 'End', 'Midpoint', 'nFly', 'Mean M KDE', 'Mean E KDE');
for j = 1:nCT
    fprintf('%-8.4g %-8.4g %-8.4g %-8d %-13.5f %-13.5f\n', ...
        ctList(j) - representativeTimeOffset, ...
        ctList(j) + representativeTimeOffset, ctList(j), nPairs(j), ...
        meanCurve_M(j), meanCurve_E(j));
    if nPairs(j) < 2
        warning('fig2cBagged:SingleFlyCT', ...
            ['CT %.4g contains one paired fly. Bootstrap resampling ' ...
             'cannot estimate within-CT biological variability there.'], ...
            ctList(j));
    end
end

recordingStartCT = ctList(:) - representativeTimeOffset;
recordingEndCT = ctList(:) + representativeTimeOffset;
ctSummaryTable = table(recordingStartCT, recordingEndCT, ctList(:), ...
    nPairs, meanCurve_M, sdCurve_M, ...
    semCurve_M, meanCurve_E, sdCurve_E, semCurve_E, ...
    'VariableNames', {'RecordingStartCT_h', 'RecordingEndCT_h', ...
    'AnalysisTime_midpoint_h', 'N_pairedFlies', 'MeanKDE_M', 'SD_M', ...
    'SEM_M', 'MeanKDE_E', 'SD_E', 'SEM_E'});

%% ========================================================================
% PART 2. Hard argmax, paired bagging, outer bootstrap, and Wald-type MET
% ========================================================================
[hardPeak_M, hardPeak_E] = computeHardPeaks(M_byCT, E_byCT, ctList);
hardDiff = hardPeak_M - hardPeak_E;

rng(rngSeed + 1, 'twister');
fprintf('\nComputing observed paired bagged peaks (nInner=%d)...\n', ...
    nInner);
[Tbag_M, Tbag_E, innerPeakTimes_M, innerPeakTimes_E, P_M, P_E] = ...
    computePairedBaggedPeaks(M_byCT, E_byCT, ctList, nInner, "finite");
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
        bootstrapPairedDataset(M_byCT, E_byCT);
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
    error('fig2cBagged:TooFewOuter', ...
        'Fewer than two valid outer replicates were obtained.');
end

% Bootstrap SEs are replicate SDs. Never divide by sqrt(nOuter).
hardSE = std(outerHardDiff(validOuterHard), 'omitnan');
hardCI95 = prctile(outerHardDiff(validOuterHard), [2.5 97.5]);
hardSD_M = std(outerHardPeak_M(validOuterHard), 'omitnan');
hardSD_E = std(outerHardPeak_E(validOuterHard), 'omitnan');
hardCI95_M = prctile(outerHardPeak_M(validOuterHard), [2.5 97.5]);
hardCI95_E = prctile(outerHardPeak_E(validOuterHard), [2.5 97.5]);
hardBootstrapSupport = ...
    mean(abs(outerHardDiff(validOuterHard)) > delta0, 'omitnan');
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
    mean(abs(outerBaggedDiff(validOuterBagged)) > delta0, 'omitnan');
thresholdExcess = abs(baggedDiff) - delta0;
standardizedEffect = baggedDiff / baggedSE; % not classical Cohen d
standardizedThresholdExcess = thresholdExcess / baggedSE;

[peakStability_M, peakEntropy_M] = computePeakStability(P_M);
[peakStability_E, peakEntropy_E] = computePeakStability(P_E);
peakConfidence_M = max(P_M);
peakConfidence_E = max(P_E);

% Monte Carlo diagnostics only; not experimental sampling SEs.
innerMCSE_M = std(innerPeakTimes_M, 'omitnan') / ...
    sqrt(sum(isfinite(innerPeakTimes_M)));
innerMCSE_E = std(innerPeakTimes_E, 'omitnan') / ...
    sqrt(sum(isfinite(innerPeakTimes_E)));
outerSEMonteCarloErrorApprox = baggedSE / ...
    sqrt(2 * max(sum(validOuterBagged) - 1, 1));

outerCovarianceMatrix = cov(outerTbag_M(validOuterBagged), ...
    outerTbag_E(validOuterBagged));
outerCovariance_M_E = outerCovarianceMatrix(1, 2);
outerCorrelationMatrix = corrcoef(outerTbag_M(validOuterBagged), ...
    outerTbag_E(validOuterBagged));
outerCorrelation_M_E = outerCorrelationMatrix(1, 2);
outerSkewness = skewness(outerBaggedDiff(validOuterBagged), 0);
outerKurtosis = kurtosis(outerBaggedDiff(validOuterBagged), 0);
if isfinite(outerKurtosis) && outerKurtosis > 0
    outerBimodalityCoefficient = ...
        (outerSkewness ^ 2 + 1) / outerKurtosis;
else
    outerBimodalityCoefficient = NaN;
end
if sum(validOuterBagged) >= 8 && baggedSE > 0
    [JB_reject, JB_p] = ...
        jbtest(outerBaggedDiff(validOuterBagged), alpha);
else
    JB_reject = NaN;
    JB_p = NaN;
end
stronglyNonNormal = abs(outerSkewness) > 1 || ...
    outerBimodalityCoefficient > 0.65 || outerKurtosis > 5;
if stronglyNonNormal
    warning('fig2cBagged:OuterNonNormal', ...
        ['Outer bootstrap distribution is strongly non-normal; consider ' ...
         'a studentized-bootstrap / bootstrap-t MET as a sensitivity ' ...
         'analysis. The Wald MET remains the primary test.']);
end

% The CT0-2 and CT24-26 windows have midpoint labels CT1 and CT25, which
% represent the same circadian phase modulo 24 h. Preserve the original
% linear argmax definition, but flag a split across both boundaries.
if numel(ctList) >= 2 && ...
        ((P_M(1) > 0.05 && P_M(end) > 0.05) || ...
         (P_E(1) > 0.05 && P_E(end) > 0.05))
    warning('fig2cBagged:CircularBoundaryMass', ...
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
statsBagged.recordingWindowHours = recordingWindowHours;
statsBagged.representativeTimeOffset = representativeTimeOffset;
statsBagged.kdeWindowMinutes = kdeWindowMinutes;
statsBagged.dataSource = 'raw XLSX workbooks';
statsBagged.sourceFolders = cellstr(inputRootFolders);
statsBagged.nSourceFolders = numel(inputRootFolders);

%% ========================================================================
% Excel output and command-window comparison
% ========================================================================
conditionTable = table(["M"; "E"], ...
    [hardPeak_M; hardPeak_E], [hardSD_M; hardSD_E], ...
    [hardCI95_M(1); hardCI95_E(1)], ...
    [hardCI95_M(2); hardCI95_E(2)], ...
    [Tbag_M; Tbag_E], [outerSD_M; outerSD_E], ...
    [baggedCI95_M(1); baggedCI95_E(1)], ...
    [baggedCI95_M(2); baggedCI95_E(2)], ...
    [peakConfidence_M; peakConfidence_E], ...
    [peakStability_M; peakStability_E], ...
    'VariableNames', {'Condition', 'HardPeak_midpointCT_h', 'HardPeak_SD_h', ...
    'HardPeak_CI95_low_h', 'HardPeak_CI95_high_h', ...
    'BaggedPeak_midpointCT_h', 'BaggedPeak_outerSD_h', ...
    'BaggedPeak_CI95_low_h', 'BaggedPeak_CI95_high_h', ...
    'PeakConfidence', 'PeakStability'});

summaryTable = table(hardDiff, hardSE, hardCI95(1), hardCI95(2), ...
    hardMET_z, hardMET_p, hardBootstrapSupport, ...
    baggedDiff, baggedSE, baggedCI95(1), baggedCI95(2), ...
    thresholdExcess, standardizedEffect, ...
    standardizedThresholdExcess, boundaryZ_lower, boundaryZ_upper, ...
    MET_z, MET_p, baggedBootstrapSupport, delta0, alpha, ...
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
    'VariableNames', {'AnalysisTime_midpoint_h', 'InnerPeakProbability_M', ...
    'InnerPeakProbability_E'});

outerBootstrapTable = table((1:nOuter).', outerHardPeak_M, ...
    outerHardPeak_E, outerHardDiff, outerTbag_M, outerTbag_E, ...
    outerBaggedDiff, ...
    'VariableNames', {'OuterReplicate', 'OuterHardPeak_M_h', ...
    'OuterHardPeak_E_h', 'OuterHardDifference_h', ...
    'OuterTbag_M_h', 'OuterTbag_E_h', ...
    'OuterBaggedDifference_MminusE_h'});

diagnosticTable = table(outerSkewness, outerKurtosis, ...
    outerBimodalityCoefficient, JB_p, double(JB_reject), ...
    double(stronglyNonNormal), outerCovariance_M_E, ...
    outerCorrelation_M_E, innerMCSE_M, innerMCSE_E, ...
    outerSEMonteCarloErrorApprox, sum(validOuterBagged), ...
    outerElapsedSeconds, ...
    'VariableNames', {'OuterSkewness', 'OuterKurtosis_normal3', ...
    'OuterBimodalityCoefficient', 'JarqueBera_p', ...
    'JarqueBeraReject', 'StrongNonNormalityFlag', ...
    'OuterCovariance_M_E', 'OuterCorrelation_M_E', ...
    'InnerMonteCarloSE_M_h', 'InnerMonteCarloSE_E_h', ...
    'ApproxMonteCarloErrorOfBaggedSE_h', 'N_validOuter', ...
    'OuterElapsedSeconds'});

sourceFolderSummary = strjoin(inputRootFolders, " | ");
configurationTable = table(numel(inputRootFolders), ...
    sourceFolderSummary, "raw XLSX", string(computeKDEFile), ...
    kdeWindowMinutes, rngSeed, nInner, ...
    nOuter, delta0, alpha, recordingWindowHours, ...
    representativeTimeOffset, mainFigureBinWidth, true, ...
    "fly", "mean of non-5th sheets", "mean of 5th sheets", ...
    'VariableNames', {'N_sourceFolders', 'SourceFolders', ...
    'DataSourceMode', 'ComputeKDEFile', 'KDEWindow_min', ...
    'RNG_seed', 'N_inner', ...
    'N_outer', 'Delta0_h', 'Alpha', 'RecordingWindow_h', ...
    'RepresentativeTimeOffset_h', 'MainFigureBinWidth_h', ...
    'PairedResampling', 'BiologicalSamplingUnit', ...
    'M_FlySummary', 'E_FlySummary'});

methodsNotes = { ...
    'Item', 'Definition'; ...
    'Source preprocessing', 'Direct import of raw signal-start XLSX workbooks; no analyzed MAT files are read'; ...
    'Multiple source roots', 'All directories listed in inputRootFolders are scanned independently'; ...
    'Cross-root CT pooling', 'Repeated CTxx folders across roots are pooled by analysis midpoint CTxx+1 before bootstrap'; ...
    'Source traceability', 'SourceRootIndex and SourceRoot remain in fly- and cell-level output tables'; ...
    'Raw file inclusion', 'All binary-recording tabs; continuous intensity is stored separately'; ...
    'Sheet classification', 'Sheet name containing 5th = E; every other sheet = M, matching fig2c.m'; ...
    'Time preprocessing', 'Subtract first sheet time and retain relative times from 0 through 2 h'; ...
    'KDE preprocessing', 'Original-sample-time compute_KDE with a 20-min Gaussian window; NaN/event binarization follows compute_KDE.m'; ...
    'Unused raw inputs', 'Intensity files and 1-min display resampling are not inputs to the peak-time estimator'; ...
    'Biological bootstrap unit', 'Fly/workbook'; ...
    'Pairing', 'Same fly index for M and E within CT at both bootstrap layers'; ...
    'M fly value', 'Mean of valid within-fly KDE means from all non-5th sheets'; ...
    'E fly value', 'Mean of valid within-fly KDE means from all 5th sheets'; ...
    'Recording window', 'Nominal CT is rounded from the first valid timestamp; each recording spans two hours'; ...
    'Representative time', 'Two-hour mean KDE is assigned to its interval midpoint CTxx+1'; ...
    'Hard peak', 'Sampled midpoint CT attaining maximum of CT-level mean KDE curve'; ...
    'Bagged peak', 'Mean of inner-bootstrap discrete hard midpoint CT values'; ...
    'Uncertainty', 'Outer nonparametric bootstrap SD and percentile CI'; ...
    'MET null', 'H0: |Delta| <= delta0'; ...
    'MET alternative', 'H1: |Delta| > delta0'; ...
    'MET p-value', 'Wald normal approximation at least-favorable interval-null boundary'; ...
    'Boundary Z values', 'Standardized boundary distances, not Student t statistics'; ...
    'Standardized effect', 'Peak difference / outer-bootstrap SE; not classical Cohen d'; ...
    'Bootstrap support', 'Descriptive only, not the reported frequentist p-value'; ...
    'Main figure', 'Outer-bootstrap distributions of Tbag_M and Tbag_E'; ...
    'Monte Carlo counts', 'nInner/nOuter affect numerical precision, not biological n'};

writetable(summaryTable, outputExcel, 'Sheet', 'Summary', ...
    'WriteMode', 'overwritesheet');
writetable(conditionTable, outputExcel, 'Sheet', 'Conditions', ...
    'WriteMode', 'overwritesheet');
writetable(ctSummaryTable, outputExcel, 'Sheet', 'CTCurves', ...
    'WriteMode', 'overwritesheet');
writetable(sourceRootsTable, outputExcel, 'Sheet', 'SourceRoots', ...
    'WriteMode', 'overwritesheet');
writetable(cellMetricsTable, outputExcel, 'Sheet', 'CellLevelData', ...
    'WriteMode', 'overwritesheet');
writetable(flyMetricsTable, outputExcel, 'Sheet', 'FlyLevelData', ...
    'WriteMode', 'overwritesheet');
writetable(peakProbabilityTable, outputExcel, ...
    'Sheet', 'PeakProbabilities', 'WriteMode', 'overwritesheet');
writetable(outerBootstrapTable, outputExcel, ...
    'Sheet', 'OuterBootstrap', 'WriteMode', 'overwritesheet');
writetable(diagnosticTable, outputExcel, 'Sheet', 'Diagnostics', ...
    'WriteMode', 'overwritesheet');
writetable(configurationTable, outputExcel, 'Sheet', 'Configuration', ...
    'WriteMode', 'overwritesheet');
writecell(methodsNotes, outputExcel, 'Sheet', 'MethodsNotes', ...
    'WriteMode', 'overwritesheet');

fprintf('\n====================================================================\n');
fprintf('fig 2c peak timing comparison (paired fly bootstrap)\n');
fprintf(['Two-hour mean KDE values are assigned to recording-window ' ...
    'midpoints (folder CTxx -> CTxx+1).\n']);
fprintf('Minimum-effect threshold: delta0 = %.3f h\n', delta0);
fprintf('H0: |Delta| <= %.3f h\n', delta0);
fprintf('H1: |Delta| >  %.3f h\n', delta0);
fprintf('====================================================================\n');
fprintf('%-34s %15s %15s\n', 'Statistic', 'Hard argmax', ...
    'Bagged peak');
fprintf('%-34s %15.3f %15.3f\n', 'Peak M (h)', hardPeak_M, Tbag_M);
fprintf('%-34s %15.3f %15.3f\n', 'Peak E (h)', hardPeak_E, Tbag_E);
fprintf('%-34s %15.3f %15.3f\n', 'Peak difference M-E (h)', ...
    hardDiff, baggedDiff);
fprintf('%-34s %15.3f %15.3f\n', 'Bootstrap SE (h)', ...
    hardSE, baggedSE);
fprintf('%-34s [%6.3f,%6.3f] [%6.3f,%6.3f]\n', ...
    '95% percentile CI (h)', hardCI95(1), hardCI95(2), ...
    baggedCI95(1), baggedCI95(2));
fprintf('%-34s %15.3f %15.3f\n', 'MET Z', hardMET_z, MET_z);
fprintf('%-34s %15.5g %15.5g\n', 'MET p-value', ...
    hardMET_p, MET_p);
fprintf('%-34s %15.3f %15.3f\n', ...
    'Bootstrap support (descriptive)', hardBootstrapSupport, ...
    baggedBootstrapSupport);
fprintf('%-34s %15s %15.3f\n', 'Peak confidence M', ...
    'N/A', peakConfidence_M);
fprintf('%-34s %15s %15.3f\n', 'Peak confidence E', ...
    'N/A', peakConfidence_E);
fprintf('%-34s %15s %15.3f\n', 'Peak stability M', ...
    'N/A', peakStability_M);
fprintf('%-34s %15s %15.3f\n', 'Peak stability E', ...
    'N/A', peakStability_E);
fprintf('====================================================================\n');
fprintf('Excel results saved to:\n  %s\n', outputExcel);
fprintf(['Bootstrap support is descriptive only. The reported frequentist ' ...
    'p-value is the Wald-type MET p-value.\n']);
fprintf(['Increasing nInner/nOuter improves Monte Carlo precision, not ' ...
    'biological sample size or power by itself.\n']);

%% ========================================================================
% PART 3. Main figure in the original fig 2c style
% ========================================================================
% The legacy figure plotted a one-level hard-peak distribution. This main
% figure instead plots outer-bootstrap sampling distributions of Tbag.
M_color = [255/255, 128/255, 0];
E_color = [0, 0, 1];
outerPeak_M = outerTbag_M(validOuterBagged);
outerPeak_E = outerTbag_E(validOuterBagged);
mainPeakEdges = makeCenteredEdges([outerPeak_M; outerPeak_E], ...
    mainFigureBinWidth);

fig = figure('Name', ...
    'Outer-bootstrap bagged peak midpoint-CT distribution', ...
    'NumberTitle', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [2, 2, 3, 4]);
hold on;
histogram(outerPeak_E, mainPeakEdges, 'FaceColor', E_color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.8, ...
    'DisplayName', 'Ecell', 'Normalization', 'probability');
histogram(outerPeak_M, mainPeakEdges, 'FaceColor', M_color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.8, ...
    'DisplayName', 'Mcell', 'Normalization', 'probability');

ax = gca;
ax.Units = 'centimeters';
ax.Position(3:4) = [1.8, 2.6];
ax.LineWidth = 0.8;
ax.TickDir = 'out';
ax.TickLength = [0.02, 0.02];
ax.XTickLabelRotation = 0;
ax.XRuler.TickLabelGapOffset = -1;
ax.FontName = 'Arial';
ax.FontSize = 7;
ax.Color = 'none';
ylim([0, 1.0]);
xlim([10, 20]);
xticks(10:5:20);
yticks(0:0.5:1.0);
ax.XMinorTick = 'on';
major_ticks = 10:5:20;
minor_ticks = setdiff(10:1:20, major_ticks);
ax.XAxis.MinorTickValues = minor_ticks;
grid off;
box off;
hold off;

set(findall(fig, '-property', 'FontUnits'), 'FontUnits', 'points');
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
print(fig, outputMainSVG, '-dsvg', '-vector');
fprintf('Main figure saved to:\n  %s\n', outputMainSVG);

%% ========================================================================
% PART 4. Numerical diagnostics and diagnostic figure
% ========================================================================
diagnosticFig = figure('Name', 'fig 2c paired bagged diagnostics', ...
    'NumberTitle', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [3, 3, 18, 12]);
diagnosticLayout = tiledlayout(diagnosticFig, 2, 3, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: observed midpoint-CT KDE curves used by the peak estimator.
nexttile(diagnosticLayout, 1);
errorbar(ctList, meanCurve_E, semCurve_E, '-o', 'Color', E_color, ...
    'MarkerFaceColor', E_color, 'MarkerSize', 3, 'LineWidth', 1.0, ...
    'DisplayName', 'E mean +/- SEM');
hold on;
errorbar(ctList, meanCurve_M, semCurve_M, '-o', 'Color', M_color, ...
    'MarkerFaceColor', M_color, 'MarkerSize', 3, 'LineWidth', 1.0, ...
    'DisplayName', 'M mean +/- SEM');
xline(Tbag_M, '--', 'M Tbag', 'Color', M_color, ...
    'HandleVisibility', 'off');
xline(Tbag_E, '--', 'E Tbag', 'Color', E_color, ...
    'HandleVisibility', 'off');
xlabel('Recording-window midpoint CT (h)');
ylabel('Fly-level mean KDE');
title('Observed midpoint-CT curves');
legend('Location', 'best');
box off;

% Panel 2: inner hard-peak probabilities used to construct Tbag.
nexttile(diagnosticLayout, 2);
barGroup = bar(ctList, [P_E(:), P_M(:)], 'grouped', ...
    'EdgeColor', 'none');
barGroup(1).FaceColor = E_color;
barGroup(2).FaceColor = M_color;
xlabel('Recording-window midpoint CT (h)');
ylabel('P_{inner}(hard peak = midpoint CT)');
title('Inner peak-location probabilities');
ylim([0 1]);
legend({'Ecell', 'Mcell'}, 'Location', 'best');
box off;

% Panel 3: outer bagged-difference distribution and null boundaries.
nexttile(diagnosticLayout, 3);
outerDiffValid = outerBaggedDiff(validOuterBagged);
diffEdges = makeCenteredEdges(outerDiffValid, diagnosticBinWidth);
histogram(outerDiffValid, diffEdges, 'Normalization', 'probability', ...
    'FaceColor', [0.45 0.45 0.45], 'EdgeColor', 'none', ...
    'FaceAlpha', 0.8);
hold on;
xline(-delta0, '--k', '-delta0', 'LineWidth', 0.9);
xline(delta0, '--k', '+delta0', 'LineWidth', 0.9);
xline(baggedDiff, '-', 'Observed diff', 'Color', [0.75 0 0], ...
    'LineWidth', 1.2);
xlabel('Outer Tbag difference M-E (h)');
ylabel('Bootstrap proportion');
title('Outer-bootstrap difference');
box off;

% Panel 4: Q-Q plot for the Wald normal approximation.
nexttile(diagnosticLayout, 4);
qqplot(outerDiffValid);
xlabel('Standard-normal quantiles');
ylabel('Outer-difference quantiles (h)');
title(sprintf('Normal Q-Q plot (JB p=%.3g)', JB_p));
box off;

% Panel 5: paired outer estimates and covariance.
nexttile(diagnosticLayout, 5);
scatter(outerTbag_M(validOuterBagged), ...
    outerTbag_E(validOuterBagged), 10, [0.25 0.25 0.25], ...
    'filled', 'MarkerFaceAlpha', 0.25);
hold on;
xline(Tbag_M, '--', 'M estimate', 'Color', M_color);
yline(Tbag_E, '--', 'E estimate', 'Color', E_color);
xlabel('Outer Tbag M (h)');
ylabel('Outer Tbag E (h)');
title(sprintf('Paired outer estimates (r=%.3f)', ...
    outerCorrelation_M_E));
box off;

% Panel 6: cumulative outer Monte Carlo convergence.
nexttile(diagnosticLayout, 6);
nValidOuter = numel(outerDiffValid);
runningMean = cumsum(outerDiffValid) ./ (1:nValidOuter).';
runningSD = NaN(nValidOuter, 1);
for k = 2:nValidOuter
    runningSD(k) = std(outerDiffValid(1:k), 'omitnan');
end
yyaxis left;
plot(1:nValidOuter, runningMean, 'Color', [0.2 0.2 0.2], ...
    'LineWidth', 0.9);
yline(baggedDiff, '--', 'Observed diff', 'Color', [0.75 0 0]);
ylabel('Cumulative mean difference (h)');
yyaxis right;
plot(1:nValidOuter, runningSD, 'Color', M_color, 'LineWidth', 0.9);
yline(baggedSE, '--', 'Final SE', 'Color', M_color);
ylabel('Cumulative bootstrap SE (h)');
xlabel('Number of outer replicates');
title('Outer Monte Carlo convergence');
box off;

set(findall(diagnosticFig, '-property', 'FontName'), ...
    'FontName', 'Arial');
set(findall(diagnosticFig, '-property', 'FontSize'), 'FontSize', 7);
diagnosticFig.PaperUnits = 'centimeters';
diagnosticFig.PaperSize = [18, 12];
diagnosticFig.PaperPosition = [0, 0, 18, 12];
diagnosticFig.PaperPositionMode = 'manual';
diagnosticFig.InvertHardcopy = 'off';
drawnow;
print(diagnosticFig, outputDiagnosticSVG, '-dsvg', '-vector');
fprintf('Diagnostic figure saved to:\n  %s\n', outputDiagnosticSVG);

%% ========================================================================
% Local helper functions
% ========================================================================

publicResult = struct('statistics', statsBagged, 'outerBootstrap', outerBootstrapTable, 'flyData', flyMetricsTable, 'cellData', cellMetricsTable, 'ctSummary', ctSummaryTable);
publicResult.configuration = cfg;
if exist('outputMAT', 'var') && isfile(outputMAT)
    save(outputMAT, 'cfg', 'publicResult', '-append');
else
    save(fullfile(cfg.outputDir, cfg.outputStem + "_results.mat"), 'cfg', 'publicResult');
end
end

function [ctList, M_byCT, E_byCT, flyTable, cellTable, ...
    sourceRootsTable] = loadPairedFlyKDEMeansFromRaw(parentFolders, ...
    recordingWindowHours, kdeWindowMinutes)
% Reproduce the raw-data path used by fig2c.m without reading analyzed MAT
% files. Each XLSX workbook is one fly. All non-5th sheets are M cells and
% all sheets containing "5th" are E cells. A folder CTxx is a recording
% window starting at CTxx; the analysis time is its interval midpoint.
% Identical CT values found under different parentFolders are pooled.

    if ~isscalar(recordingWindowHours) || ...
            ~isfinite(recordingWindowHours) || recordingWindowHours <= 0
        error('fig2cBagged:InvalidRecordingWindow', ...
            'recordingWindowHours must be a positive finite scalar.');
    end
    if ~isscalar(kdeWindowMinutes) || ...
            ~isfinite(kdeWindowMinutes) || kdeWindowMinutes <= 0
        error('fig2cBagged:InvalidKDEWindow', ...
            'kdeWindowMinutes must be a positive finite scalar.');
    end

    parentFolders = string(parentFolders(:));
    parentFolders = strtrim(parentFolders);
    parentFolders = parentFolders(strlength(parentFolders) > 0);
    if isempty(parentFolders)
        error('fig2cBagged:NoInputRoots', ...
            'At least one dataset root is required.');
    end
    if numel(unique(lower(parentFolders))) ~= numel(parentFolders)
        error('fig2cBagged:DuplicateInputRoot', ...
            'The same dataset root was supplied more than once.');
    end
    missingRoot = ~bootstrapRawIsFolder(parentFolders);
    if any(missingRoot)
        error('fig2cBagged:MissingDataFolder', ...
            'Data folder not found: %s', ...
            strjoin(parentFolders(missingRoot), ' | '));
    end

    nSourceRoots = numel(parentFolders);
    subfoldersByRoot = cell(nSourceRoots, 1);
    nFolderCandidates = 0;
    for r = 1:nSourceRoots
        subfolders = bootstrapRawFolders(parentFolders(r));
        subfolders = subfolders([subfolders.isdir] & ...
            ~ismember({subfolders.name}, {'.', '..'}));
        subfoldersByRoot{r} = subfolders;
        nFolderCandidates = nFolderCandidates + numel(subfolders);
    end

    folderTemplate = struct('sourceRootIndex', NaN, 'sourceRoot', '', ...
        'recordingStartCT', NaN, ...
        'recordingEndCT', NaN, 'analysisTime', NaN, ...
        'folderPath', '', 'files', []);
    folderRecords = repmat(folderTemplate, nFolderCandidates, 1);
    nValidFolders = 0;
    nFlyTotal = 0;
    nSheetTotal = 0;

    % First pass discovers the complete raw workbook/sheet set so all
    % output columns can be preallocated before reading numerical data.
    for r = 1:nSourceRoots
        parentFolder = parentFolders(r);
        subfolders = subfoldersByRoot{r};
        for k = 1:numel(subfolders)
            folderName = subfolders(k).name;
            ctToken = regexp(folderName, 'CT(\d+\.?\d*)', ...
                'tokens', 'once');
            if isempty(ctToken)
                continue;
            end
            recordingStartCT = str2double(ctToken{1});
            if ~isfinite(recordingStartCT)
                continue;
            end
            folderPath = fullfile(parentFolder, folderName);
            files = listRawDataWorkbooks(folderPath);
            if isempty(files)
                warning('fig2cBagged:NoRawWorkbook', ...
                    'No raw signal-start XLSX workbook found in %s.', ...
                    folderPath);
                continue;
            end

            fileTemplate = struct('name', '', 'path', '', 'sheets', {{}});
            fileRecords = repmat(fileTemplate, numel(files), 1);
            for f = 1:numel(files)
                sourcePath = fullfile(folderPath, files(f).name);
                try
                    sheets = bootstrapRawSheets(sourcePath);
                catch ME
                    error('fig2cBagged:SheetDiscoveryFailed', ...
                        'Cannot list sheets in %s: %s', ...
                        sourcePath, ME.message);
                end
                fileRecords(f).name = files(f).name;
                fileRecords(f).path = sourcePath;
                fileRecords(f).sheets = sheets;
                nSheetTotal = nSheetTotal + numel(sheets);
            end

            nValidFolders = nValidFolders + 1;
            folderRecords(nValidFolders).sourceRootIndex = r;
            folderRecords(nValidFolders).sourceRoot = char(parentFolder);
            folderRecords(nValidFolders).recordingStartCT = ...
                recordingStartCT;
            folderRecords(nValidFolders).recordingEndCT = ...
                recordingStartCT + recordingWindowHours;
            folderRecords(nValidFolders).analysisTime = ...
                recordingStartCT + recordingWindowHours / 2;
            folderRecords(nValidFolders).folderPath = char(folderPath);
            folderRecords(nValidFolders).files = fileRecords;
            nFlyTotal = nFlyTotal + numel(fileRecords);
        end
    end
    folderRecords = folderRecords(1:nValidFolders);
    if nFlyTotal == 0 || nSheetTotal == 0
        error('fig2cBagged:NoRawData', ...
            'No usable raw XLSX workbooks were found under any input root.');
    end

    % Fly-level columns.
    flySourceRootIndex = zeros(nFlyTotal, 1);
    flySourceRoot = strings(nFlyTotal, 1);
    flyStartCT = NaN(nFlyTotal, 1);
    flyEndCT = NaN(nFlyTotal, 1);
    flyAnalysisTime = NaN(nFlyTotal, 1);
    flyIDColumn = strings(nFlyTotal, 1);
    flySourceFile = strings(nFlyTotal, 1);
    nMCellColumn = zeros(nFlyTotal, 1);
    nECellColumn = zeros(nFlyTotal, 1);
    MMetricColumn = NaN(nFlyTotal, 1);
    EMetricColumn = NaN(nFlyTotal, 1);
    validPairColumn = false(nFlyTotal, 1);

    % Cell/sheet-level audit columns.
    cellSourceRootIndex = zeros(nSheetTotal, 1);
    cellSourceRoot = strings(nSheetTotal, 1);
    cellStartCT = NaN(nSheetTotal, 1);
    cellEndCT = NaN(nSheetTotal, 1);
    cellAnalysisTime = NaN(nSheetTotal, 1);
    cellFlyID = strings(nSheetTotal, 1);
    cellSourceFile = strings(nSheetTotal, 1);
    cellSheetName = strings(nSheetTotal, 1);
    cellCondition = strings(nSheetTotal, 1);
    cellRawRows = zeros(nSheetTotal, 1);
    cellWindowRows = zeros(nSheetTotal, 1);
    cellTotalEvents = zeros(nSheetTotal, 1);
    cellMeanKDE = NaN(nSheetTotal, 1);
    cellValidKDE = false(nSheetTotal, 1);

    flyIndex = 0;
    cellIndex = 0;
    for k = 1:numel(folderRecords)
        folderRecord = folderRecords(k);
        for f = 1:numel(folderRecord.files)
            fileRecord = folderRecord.files(f);
            sourcePath = fileRecord.path;
            sheets = fileRecord.sheets;
            [~, flyID, ~] = fileparts(fileRecord.name);

            mCellMetrics = NaN(numel(sheets), 1);
            eCellMetrics = NaN(numel(sheets), 1);
            nM = 0;
            nE = 0;
            for s = 1:numel(sheets)
                sheetName = sheets{s};
                [metric, nRawRows, nWindowRows, totalEvents] = ...
                    readRawSheetKDE(sourcePath, sheetName, ...
                    recordingWindowHours, kdeWindowMinutes);
                isE = contains(sheetName, '5th', 'IgnoreCase', true);
                if isE
                    nE = nE + 1;
                    eCellMetrics(nE) = metric;
                    condition = "E";
                else
                    nM = nM + 1;
                    mCellMetrics(nM) = metric;
                    condition = "M";
                end

                cellIndex = cellIndex + 1;
                cellSourceRootIndex(cellIndex) = ...
                    folderRecord.sourceRootIndex;
                cellSourceRoot(cellIndex) = string(folderRecord.sourceRoot);
                cellStartCT(cellIndex) = folderRecord.recordingStartCT;
                cellEndCT(cellIndex) = folderRecord.recordingEndCT;
                cellAnalysisTime(cellIndex) = folderRecord.analysisTime;
                cellFlyID(cellIndex) = string(flyID);
                cellSourceFile(cellIndex) = string(sourcePath);
                cellSheetName(cellIndex) = string(sheetName);
                cellCondition(cellIndex) = condition;
                cellRawRows(cellIndex) = nRawRows;
                cellWindowRows(cellIndex) = nWindowRows;
                cellTotalEvents(cellIndex) = totalEvents;
                cellMeanKDE(cellIndex) = metric;
                cellValidKDE(cellIndex) = isfinite(metric);
            end

            mFlyMetric = mean(mCellMetrics(1:nM), 'omitnan');
            eFlyMetric = mean(eCellMetrics(1:nE), 'omitnan');
            validPair = isfinite(mFlyMetric) && isfinite(eFlyMetric);

            flyIndex = flyIndex + 1;
            flySourceRootIndex(flyIndex) = folderRecord.sourceRootIndex;
            flySourceRoot(flyIndex) = string(folderRecord.sourceRoot);
            flyStartCT(flyIndex) = folderRecord.recordingStartCT;
            flyEndCT(flyIndex) = folderRecord.recordingEndCT;
            flyAnalysisTime(flyIndex) = folderRecord.analysisTime;
            flyIDColumn(flyIndex) = string(flyID);
            flySourceFile(flyIndex) = string(sourcePath);
            nMCellColumn(flyIndex) = nM;
            nECellColumn(flyIndex) = nE;
            MMetricColumn(flyIndex) = mFlyMetric;
            EMetricColumn(flyIndex) = eFlyMetric;
            validPairColumn(flyIndex) = validPair;

            if nM ~= 4 || nE ~= 1
                warning('fig2cBagged:UnexpectedCellCount', ...
                    ['%s has %d M and %d E sheets. Available valid cells ' ...
                     'are aggregated within fly; inspect CellLevelData.'], ...
                    sourcePath, nM, nE);
            end
            if ~validPair
                warning('fig2cBagged:IncompleteFlyPair', ...
                    'Excluding incomplete fly pair: %s', sourcePath);
            end
        end
    end
    if flyIndex ~= nFlyTotal || cellIndex ~= nSheetTotal
        error('fig2cBagged:RawImportCountMismatch', ...
            'Internal raw-import count mismatch.');
    end

    flyTable = table(flySourceRootIndex, flySourceRoot, ...
        flyStartCT, flyEndCT, flyAnalysisTime, flyIDColumn, ...
        flySourceFile, nMCellColumn, nECellColumn, MMetricColumn, ...
        EMetricColumn, validPairColumn, 'VariableNames', ...
        {'SourceRootIndex', 'SourceRoot', 'RecordingStartCT_h', ...
        'RecordingEndCT_h', 'AnalysisTime_midpoint_h', 'FlyID', ...
        'SourceFile', 'N_M_cells', 'N_E_cells', 'FlyMeanKDE_M', ...
        'FlyMeanKDE_E', 'ValidPair'});
    flyTable = sortrows(flyTable, ...
        {'AnalysisTime_midpoint_h', 'SourceRootIndex', 'FlyID'});

    cellTable = table(cellSourceRootIndex, cellSourceRoot, ...
        cellStartCT, cellEndCT, cellAnalysisTime, cellFlyID, ...
        cellSourceFile, cellSheetName, cellCondition, cellRawRows, ...
        cellWindowRows, cellTotalEvents, cellMeanKDE, cellValidKDE, ...
        'VariableNames', {'SourceRootIndex', 'SourceRoot', ...
        'RecordingStartCT_h', 'RecordingEndCT_h', ...
        'AnalysisTime_midpoint_h', 'FlyID', 'SourceFile', ...
        'SheetName', 'Condition', 'N_rawRows', ...
        'N_rowsInWindow', 'TotalEventsInWindow', 'MeanKDE', 'ValidKDE'});
    cellTable = sortrows(cellTable, {'AnalysisTime_midpoint_h', ...
        'SourceRootIndex', 'FlyID', 'Condition', 'SheetName'});

    sourceRootIndex = (1:nSourceRoots).';
    sourceRoot = parentFolders;
    nCTFolders = zeros(nSourceRoots, 1);
    nWorkbooks = zeros(nSourceRoots, 1);
    nValidPairedFlies = zeros(nSourceRoots, 1);
    nSheets = zeros(nSourceRoots, 1);
    for r = 1:nSourceRoots
        nCTFolders(r) = sum([folderRecords.sourceRootIndex] == r);
        nWorkbooks(r) = sum(flyTable.SourceRootIndex == r);
        nValidPairedFlies(r) = sum(flyTable.SourceRootIndex == r & ...
            flyTable.ValidPair);
        nSheets(r) = sum(cellTable.SourceRootIndex == r);
    end
    sourceRootsTable = table(sourceRootIndex, sourceRoot, nCTFolders, ...
        nWorkbooks, nValidPairedFlies, nSheets, ...
        'VariableNames', {'SourceRootIndex', 'SourceRoot', ...
        'N_CTFolders', 'N_workbooks', 'N_validPairedFlies', 'N_sheets'});

    validRows = flyTable.ValidPair & ...
        isfinite(flyTable.AnalysisTime_midpoint_h);
    if ~any(validRows)
        error('fig2cBagged:NoValidPairs', ...
            'No complete M/E fly pairs were recovered from raw workbooks.');
    end
    ctList = unique(flyTable.AnalysisTime_midpoint_h(validRows));
    ctList = ctList(:).';
    M_byCT = cell(1, numel(ctList));
    E_byCT = cell(1, numel(ctList));
    for j = 1:numel(ctList)
        useRows = validRows & ...
            flyTable.AnalysisTime_midpoint_h == ctList(j);
        M_byCT{j} = flyTable.FlyMeanKDE_M(useRows);
        E_byCT{j} = flyTable.FlyMeanKDE_E(useRows);
        if numel(M_byCT{j}) ~= numel(E_byCT{j})
            error('fig2cBagged:PairLengthMismatch', ...
                'M/E fly counts differ at midpoint CT %.4g.', ctList(j));
        end
    end
end

function files = listRawDataWorkbooks(folderPath)
% Match the original fig2c.m inclusion rule while ignoring Excel lock files.
    files = bootstrapRawFiles(folderPath);
    if isempty(files)
        return;
    end
    names = string({files.name});
    keep = ~startsWith(names, '~$') & ...
        ~contains(names, 'intensity', 'IgnoreCase', true) & ...
        ~contains(names, '_KDE', 'IgnoreCase', true);
    files = files(keep);
    if isempty(files)
        return;
    end
    [~, order] = sort(lower(string({files.name})));
    files = files(order);
end

function [metric, nRawRows, nWindowRows, totalEvents] = ...
    readRawSheetKDE(sourcePath, sheetName, recordingWindowHours, ...
    kdeWindowMinutes)
% Port of the statistical path in fig2c.m: read columns 1-2, align the
% sheet to its first time, retain 0-2 h, compute KDE at original sample
% times, then average that KDE. compute_KDE.m retains the original event
% binarization and NaN behavior.
    try
        sheetData = bootstrapRawReadSheet(sourcePath, sheetName);
    catch ME
        error('fig2cBagged:RawSheetReadFailed', ...
            'Cannot read %s [%s]: %s', sourcePath, sheetName, ME.message);
    end
    if width(sheetData) < 2
        error('fig2cBagged:TooFewRawColumns', ...
            '%s [%s] must contain time and event columns.', ...
            sourcePath, sheetName);
    end
    timeRaw = sheetData{:, 1};
    eventsRaw = sheetData{:, 2};
    if ~isnumeric(timeRaw) || ~isnumeric(eventsRaw)
        error('fig2cBagged:NonNumericRawColumns', ...
            '%s [%s] has nonnumeric time/event columns.', ...
            sourcePath, sheetName);
    end
    timeRaw = timeRaw(:);
    eventsRaw = eventsRaw(:);
    if numel(timeRaw) ~= numel(eventsRaw)
        error('fig2cBagged:RawLengthMismatch', ...
            '%s [%s] has mismatched time/event lengths.', ...
            sourcePath, sheetName);
    end

    nRawRows = numel(timeRaw);
    nWindowRows = 0;
    totalEvents = 0;
    metric = NaN;
    if isempty(timeRaw)
        return;
    end
    timeRelative = timeRaw - timeRaw(1);
    keep = timeRelative >= 0 & ...
        timeRelative <= recordingWindowHours;
    timeWindow = timeRelative(keep);
    eventsWindow = eventsRaw(keep);
    nWindowRows = numel(timeWindow);
    totalEvents = sum(eventsWindow > 0);
    if isempty(timeWindow)
        return;
    end

    kdeValues = compute_KDE(timeWindow * 60, eventsWindow, ...
        kdeWindowMinutes);
    metric = mean(kdeValues(:), 'omitnan');
end
