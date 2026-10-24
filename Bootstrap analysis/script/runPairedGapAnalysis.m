function publicResult = runPairedGapAnalysis(cfg)
% Shared publication pipeline; settings are supplied by the figure entry script.
figureCleanup = bootstrapFigureContext(cfg); %#ok<NASGU>
filePath = cfg.inputFile;
folderPath = cfg.outputDir;
outputExcel = fullfile(folderPath, cfg.outputStem + "_results.xlsx");
outputMainSVG = fullfile(folderPath, cfg.outputStem + ".svg");
outputDiagnosticSVG = fullfile(folderPath, cfg.outputStem + "_diagnostics.svg");
nInner = cfg.nInner;
nOuter = cfg.nOuter;
delta0 = cfg.delta0;
alpha = cfg.alpha;
rngSeed = cfg.rngSeed;
mainFigureBinWidth = cfg.mainFigureBinWidth;

rawdata = readmatrix(filePath, 'Sheet', cfg.dataSheet);
zt = rawdata(:, cfg.timeColumn);
Mcell = rawdata(:, cfg.MColumn);
Ecell = rawdata(:, cfg.EColumn);

% Exclude nonnumeric/header rows from the time list. Within each time point,
% keep complete M/E pairs so a common row index is always meaningful.
ztList = unique(zt(~isnan(zt)));
nZT = numel(ztList);
M_perSample = cell(1, nZT);
E_perSample = cell(1, nZT);

nPairs = NaN(nZT, 1);
mean_M = NaN(nZT, 1);
sd_M = NaN(nZT, 1);
mean_E = NaN(nZT, 1);
sd_E = NaN(nZT, 1);
pairedCorrelation = NaN(nZT, 1);

fprintf('========== Paired data summary by time point ==========\n');
fprintf('%-8s %-8s %-12s %-12s %-12s\n', ...
    'Time', 'nPair', 'Mean M', 'Mean E', 'Pair r');

for j = 1:nZT
    rowAtTime = (zt == ztList(j));
    Mj = Mcell(rowAtTime);
    Ej = Ecell(rowAtTime);
    completePair = ~isnan(Mj) & ~isnan(Ej);

    if any(xor(isnan(Mj), isnan(Ej)))
        warning('bootstrap:IncompletePair', ...
            ['Time %.4g contains incomplete M/E rows. Those rows are ' ...
             'excluded so all bootstrap observations remain paired.'], ...
            ztList(j));
    end

    Mj = Mj(completePair);
    Ej = Ej(completePair);
    Mj = Mj(:);
    Ej = Ej(:);

    if isempty(Mj)
        error('bootstrap:NoPairsAtTime', ...
            'No complete M/E pairs are available at time %.4g.', ztList(j));
    end

    M_perSample{j} = Mj;
    E_perSample{j} = Ej;
    nPairs(j) = numel(Mj);
    mean_M(j) = mean(Mj);
    mean_E(j) = mean(Ej);
    sd_M(j) = std(Mj);
    sd_E(j) = std(Ej);

    if numel(Mj) > 2 && std(Mj) > 0 && std(Ej) > 0
        pairCorrMatrix = corrcoef(Mj, Ej);
        pairedCorrelation(j) = pairCorrMatrix(1, 2);
    end

    fprintf('%-8.4g %-8d %-12.3f %-12.3f %-12.3f\n', ...
        ztList(j), nPairs(j), mean_M(j), mean_E(j), ...
        pairedCorrelation(j));
end

timepointTable = table(ztList(:), nPairs, mean_M, sd_M, mean_E, sd_E, ...
    pairedCorrelation, 'VariableNames', { ...
    'Time_h', 'N_pairs', 'Mean_M', 'SD_M', 'Mean_E', 'SD_E', ...
    'PairedCorrelation_M_E'});

% Observed discrete hard-argmax peaks are retained as a useful reference.
[~, hardIndex_M] = max(mean_M);
[~, hardIndex_E] = max(mean_E);
hardPeak_M = ztList(hardIndex_M);
hardPeak_E = ztList(hardIndex_E);
hardDiff = hardPeak_M - hardPeak_E;

%% ========================================================================
% PART 2. Fully paired inner/outer bagging bootstrap and Excel output
% ========================================================================
% Inner bootstrap:
%   - At each time point, resample paired rows using one common index matrix.
%   - Compute paired bootstrap mean curves for M and E.
%   - Retain the hard argmax in every replicate.
%   - Average the discrete peak times to obtain Tbag_M and Tbag_E.

fprintf('\nComputing paired bootstrap-smoothed / bagged peaks (nInner=%d)...\n', ...
    nInner);
rng(rngSeed + 1, 'twister');
[Tbag_M, Tbag_E, innerPeakTimes_M, innerPeakTimes_E, P_M, P_E] = ...
    computePairedBaggedPeaks(M_perSample, E_perSample, ztList, nInner);
baggedDiff = Tbag_M - Tbag_E;

% Outer nonparametric bootstrap:
%   - Resample paired rows once to form each outer dataset.
%   - Re-run the complete paired inner-bagging procedure on that dataset.
%   - The across-outer SD is the sampling SE of the bagged estimator.
outerTbag_M = NaN(nOuter, 1);
outerTbag_E = NaN(nOuter, 1);
rng(rngSeed + 2, 'twister');
outerStart = tic;
progressStep = max(1, floor(nOuter / 20));

fprintf('Running paired outer nonparametric bootstrap (nOuter=%d)...\n', ...
    nOuter);
for b = 1:nOuter
    [outerData_M, outerData_E] = ...
        bootstrapPairedDataset(M_perSample, E_perSample);
    [outerTbag_M(b), outerTbag_E(b)] = ...
        computePairedBaggedPeaks(outerData_M, outerData_E, ztList, nInner);

    if mod(b, progressStep) == 0 || b == nOuter
        fprintf('  outer bootstrap: %d/%d complete\n', b, nOuter);
    end
end
outerElapsedSeconds = toc(outerStart);

outerBaggedDiff = outerTbag_M - outerTbag_E;
validOuter = isfinite(outerTbag_M) & isfinite(outerTbag_E) & ...
    isfinite(outerBaggedDiff);

outerSD_M = std(outerTbag_M(validOuter));
outerSD_E = std(outerTbag_E(validOuter));
baggedSE = std(outerBaggedDiff(validOuter));

CI95_M = prctile(outerTbag_M(validOuter), [2.5 97.5]);
CI95_E = prctile(outerTbag_E(validOuter), [2.5 97.5]);
CI95_diff = prctile(outerBaggedDiff(validOuter), [2.5 97.5]);

[MET_z, MET_p, MET_z_near, MET_z_far, ...
    boundaryZ_lower, boundaryZ_upper] = ...
    computeWaldMET(baggedDiff, baggedSE, delta0);

bootstrapSupport = mean(abs(outerBaggedDiff(validOuter)) > delta0);
thresholdExcess = abs(baggedDiff) - delta0;
standardizedEffect = baggedDiff / baggedSE;
standardizedThresholdExcess = thresholdExcess / baggedSE;

[peakStability_M, peakEntropy_M] = computePeakStability(P_M);
[peakStability_E, peakEntropy_E] = computePeakStability(P_E);
peakConfidence_M = max(P_M);
peakConfidence_E = max(P_E);

% Numerical-integration diagnostics only. These quantities use sqrt of a
% Monte Carlo replicate count solely to describe simulation precision; they
% are not biological sampling SEs and are not used in MET inference.
innerMCSE_M = std(innerPeakTimes_M, 'omitnan') / ...
    sqrt(sum(isfinite(innerPeakTimes_M)));
innerMCSE_E = std(innerPeakTimes_E, 'omitnan') / ...
    sqrt(sum(isfinite(innerPeakTimes_E)));
outerSEMonteCarloErrorApprox = baggedSE / ...
    sqrt(2 * max(sum(validOuter) - 1, 1));

% Joint outer-bootstrap covariance is the information lost by an
% independent-M/E bootstrap. It enters Var(M-E)=Var(M)+Var(E)-2Cov(M,E).
outerCovarianceMatrix = cov(outerTbag_M(validOuter), ...
    outerTbag_E(validOuter));
outerCorrelationMatrix = corrcoef(outerTbag_M(validOuter), ...
    outerTbag_E(validOuter));
outerCovariance_M_E = outerCovarianceMatrix(1, 2);
outerCorrelation_M_E = outerCorrelationMatrix(1, 2);

% Distribution-shape diagnostics for the Wald normal approximation.
outerSkewness = skewness(outerBaggedDiff(validOuter), 0);
outerKurtosis = kurtosis(outerBaggedDiff(validOuter), 0); % normal = 3
outerBimodalityCoefficient = ...
    (outerSkewness ^ 2 + 1) / outerKurtosis;

if sum(validOuter) >= 8 && baggedSE > 0
    [JB_reject, JB_p] = jbtest(outerBaggedDiff(validOuter), alpha);
else
    JB_reject = NaN;
    JB_p = NaN;
end

% Flag pronounced skewness, multimodality, or heavy tails. This is a
% sensitivity warning only; the requested Wald MET remains the main test.
stronglyNonNormal = abs(outerSkewness) > 1 || ...
    outerBimodalityCoefficient > 0.70 || outerKurtosis > 6;
if stronglyNonNormal
    warning('bootstrap:OuterBootstrapNonNormal', [ ...
        'Outer bootstrap distribution is strongly non-normal; consider a ' ...
        'studentized-bootstrap / bootstrap-t MET as a sensitivity analysis.']);
end

if MET_p < alpha
    MET_decision = "Reject H0: absolute peak shift exceeds delta0";
else
    MET_decision = "Do not reject H0: minimum-effect evidence is insufficient";
end

% Store the complete result in a workspace structure.
statsBaggedPaired = struct();
statsBaggedPaired.estimator = ...
    'Bootstrap-smoothed / bagged paired discrete-argmax peak estimator';
statsBaggedPaired.peak_M = Tbag_M;
statsBaggedPaired.peak_E = Tbag_E;
statsBaggedPaired.diff = baggedDiff;
statsBaggedPaired.outerSD_M = outerSD_M;
statsBaggedPaired.outerSD_E = outerSD_E;
statsBaggedPaired.SE = baggedSE;
statsBaggedPaired.CI95_M = CI95_M;
statsBaggedPaired.CI95_E = CI95_E;
statsBaggedPaired.CI95_diff = CI95_diff;
statsBaggedPaired.MET_z = MET_z;
statsBaggedPaired.MET_p = MET_p;
statsBaggedPaired.MET_z_near = MET_z_near;
statsBaggedPaired.MET_z_far = MET_z_far;
statsBaggedPaired.boundaryZ_lower = boundaryZ_lower;
statsBaggedPaired.boundaryZ_upper = boundaryZ_upper;
statsBaggedPaired.delta0 = delta0;
statsBaggedPaired.alpha = alpha;
statsBaggedPaired.thresholdExcess = thresholdExcess;
statsBaggedPaired.standardizedEffect = standardizedEffect;
statsBaggedPaired.standardizedThresholdExcess = ...
    standardizedThresholdExcess;
statsBaggedPaired.bootstrapSupport = bootstrapSupport;
statsBaggedPaired.peakProbability_M = P_M;
statsBaggedPaired.peakProbability_E = P_E;
statsBaggedPaired.peakConfidence_M = peakConfidence_M;
statsBaggedPaired.peakConfidence_E = peakConfidence_E;
statsBaggedPaired.peakStability_M = peakStability_M;
statsBaggedPaired.peakStability_E = peakStability_E;
statsBaggedPaired.outerCovariance_M_E = outerCovariance_M_E;
statsBaggedPaired.outerCorrelation_M_E = outerCorrelation_M_E;
statsBaggedPaired.outerSkewness = outerSkewness;
statsBaggedPaired.outerKurtosis = outerKurtosis;
statsBaggedPaired.outerBimodalityCoefficient = ...
    outerBimodalityCoefficient;
statsBaggedPaired.JB_p = JB_p;
statsBaggedPaired.JB_reject = JB_reject;
statsBaggedPaired.nInner = nInner;
statsBaggedPaired.nOuter = nOuter;
statsBaggedPaired.rngSeed = rngSeed;
statsBaggedPaired.outerElapsedSeconds = outerElapsedSeconds;

fprintf('\n====================================================================\n');
fprintf('Fully paired bootstrap-smoothed / bagged peak results\n');
fprintf('====================================================================\n');
fprintf('M peak:          %.4f +/- %.4f h, 95%% CI [%.4f, %.4f]\n', ...
    Tbag_M, outerSD_M, CI95_M(1), CI95_M(2));
fprintf('E peak:          %.4f +/- %.4f h, 95%% CI [%.4f, %.4f]\n', ...
    Tbag_E, outerSD_E, CI95_E(1), CI95_E(2));
fprintf('M-E difference:  %.4f +/- %.4f h, 95%% CI [%.4f, %.4f]\n', ...
    baggedDiff, baggedSE, CI95_diff(1), CI95_diff(2));
fprintf('Minimum-effect threshold delta0: %.4f h\n', delta0);
fprintf('MET Z: %.4f\n', MET_z);
fprintf('MET p-value: %.6g\n', MET_p);
fprintf('Bootstrap support |diff| > delta0: %.4f (descriptive only)\n', ...
    bootstrapSupport);
fprintf('Outer paired correlation of Tbag_M and Tbag_E: %.4f\n', ...
    outerCorrelation_M_E);
fprintf('Decision at alpha=%.3f: %s\n', alpha, MET_decision);
fprintf('Outer bootstrap elapsed time: %.1f seconds\n', outerElapsedSeconds);
fprintf('====================================================================\n');

% Compact, human-readable report table. For M and E, ``Uncertainty`` is
% the SD across paired outer-bootstrap bagged peak estimates; for M-E it
% is the paired outer-bootstrap SE of the difference. NaN means that the
% column is not applicable to that row.
summaryStatistic = [ ...
    "Bagged peak M"; "Bagged peak E"; "Peak difference M-E"; ...
    "Minimum-effect threshold"; "Threshold excess"; ...
    "Standardized peak-shift statistic"; ...
    "Standardized threshold excess"; "MET Z"; "MET p-value"; ...
    "Nearest-boundary Z"; "Far-boundary Z"; ...
    "Boundary Z lower"; "Boundary Z upper"; ...
    "Bootstrap support"; "Peak confidence M"; "Peak confidence E"; ...
    "Peak stability M"; "Peak stability E"; ...
    "Outer paired correlation"; "nInner"; "nOuter"; "MET decision"];
summaryEstimate = [Tbag_M; Tbag_E; baggedDiff; delta0; ...
    thresholdExcess; standardizedEffect; standardizedThresholdExcess; ...
    MET_z; MET_p; MET_z_near; MET_z_far; boundaryZ_lower; ...
    boundaryZ_upper; bootstrapSupport; peakConfidence_M; ...
    peakConfidence_E; peakStability_M; peakStability_E; ...
    outerCorrelation_M_E; nInner; nOuter; NaN];
summaryUncertainty = [outerSD_M; outerSD_E; baggedSE; ...
    NaN(19, 1)];
summaryCI_low = [CI95_M(1); CI95_E(1); CI95_diff(1); NaN(19, 1)];
summaryCI_high = [CI95_M(2); CI95_E(2); CI95_diff(2); NaN(19, 1)];
summaryUnit = ["h"; "h"; "h"; "h"; "h"; repmat("", 17, 1)];
summaryNote = [ ...
    "Estimate +/- outer-bootstrap SD"; ...
    "Estimate +/- outer-bootstrap SD"; ...
    "Estimate +/- paired outer-bootstrap SE"; ...
    "delta0; defines H0: |Delta| <= delta0"; ...
    "|M-E| - delta0"; ...
    "Not classical Cohen's d"; ...
    "Equals MET Z when the observed shift exceeds delta0"; ...
    "Wald distance from the nearest interval-null boundary"; ...
    "Reported frequentist p-value; Wald normal approximation"; ...
    "Standardized distance to nearest null boundary"; ...
    "Standardized distance to far null boundary"; ...
    "Not a Student t statistic"; "Not a Student t statistic"; ...
    "Descriptive only; not a p-value"; ...
    "Maximum inner peak-location probability"; ...
    "Maximum inner peak-location probability"; ...
    "1=stable, 0=maximally diffuse"; ...
    "1=stable, 0=maximally diffuse"; ...
    "Correlation preserves paired outer covariance"; ...
    "Monte Carlo precision, not biological sample size"; ...
    "Monte Carlo precision, not biological sample size"; ...
    string(MET_decision)];
summaryTable = table(summaryStatistic, summaryEstimate, ...
    summaryUncertainty, summaryCI_low, summaryCI_high, summaryUnit, ...
    summaryNote, 'VariableNames', {'Statistic', 'Estimate', ...
    'Uncertainty_SD_or_SE', 'CI95_low', 'CI95_high', 'Unit', 'Note'});

% Full machine-readable one-row table retained for downstream scripts.
fullSummaryTable = table(Tbag_M, outerSD_M, CI95_M(1), CI95_M(2), ...
    Tbag_E, outerSD_E, CI95_E(1), CI95_E(2), ...
    baggedDiff, baggedSE, CI95_diff(1), CI95_diff(2), ...
    delta0, thresholdExcess, standardizedEffect, ...
    standardizedThresholdExcess, MET_z, MET_p, ...
    MET_z_near, MET_z_far, boundaryZ_lower, boundaryZ_upper, ...
    bootstrapSupport, peakConfidence_M, peakConfidence_E, ...
    peakStability_M, peakStability_E, outerCorrelation_M_E, ...
    nInner, nOuter, string(MET_decision), ...
    'VariableNames', { ...
    'BaggedPeak_M_h', 'OuterSD_M_h', 'CI95_M_low_h', 'CI95_M_high_h', ...
    'BaggedPeak_E_h', 'OuterSD_E_h', 'CI95_E_low_h', 'CI95_E_high_h', ...
    'Difference_M_minus_E_h', 'Difference_SE_h', ...
    'Difference_CI95_low_h', 'Difference_CI95_high_h', ...
    'Delta0_h', 'ThresholdExcess_h', 'StandardizedPeakShift', ...
    'StandardizedThresholdExcess', 'MET_Z', 'MET_p', ...
    'MET_Z_near', 'MET_Z_far', 'BoundaryZ_lower', 'BoundaryZ_upper', ...
    'BootstrapSupport_descriptive', 'PeakConfidence_M', ...
    'PeakConfidence_E', 'PeakStability_M', 'PeakStability_E', ...
    'OuterCorrelation_M_E', 'N_inner', 'N_outer', 'MET_decision'});

peakProbabilityTable = table(ztList(:), P_M(:), P_E(:), ...
    'VariableNames', {'Time_h', 'PeakProbability_M', 'PeakProbability_E'});

outerBootstrapTable = table((1:nOuter).', outerTbag_M, outerTbag_E, ...
    outerBaggedDiff, 'VariableNames', { ...
    'OuterReplicate', 'OuterTbag_M_h', 'OuterTbag_E_h', ...
    'OuterDifference_M_minus_E_h'});

conditionDiagnosticsTable = table(["M"; "E"], ...
    [hardPeak_M; hardPeak_E], [Tbag_M; Tbag_E], ...
    [outerSD_M; outerSD_E], [CI95_M(1); CI95_E(1)], ...
    [CI95_M(2); CI95_E(2)], ...
    [std(innerPeakTimes_M, 'omitnan'); std(innerPeakTimes_E, 'omitnan')], ...
    [innerMCSE_M; innerMCSE_E], ...
    [peakConfidence_M; peakConfidence_E], ...
    [peakEntropy_M; peakEntropy_E], ...
    [peakStability_M; peakStability_E], ...
    'VariableNames', {'Condition', 'ObservedHardPeak_h', 'BaggedPeak_h', ...
    'OuterBootstrapSD_h', 'CI95_low_h', 'CI95_high_h', ...
    'InnerHardPeakSD_h', 'InnerMonteCarloSE_h', 'PeakConfidence', ...
    'PeakEntropy', 'PeakStability'});

diagnosticNames = [ ...
    "Outer skewness"; ...
    "Outer kurtosis (normal=3)"; ...
    "Outer bimodality coefficient"; ...
    "Jarque-Bera p-value"; ...
    "Jarque-Bera reject at alpha"; ...
    "Strong non-normality heuristic flag"; ...
    "Outer covariance Tbag_M,Tbag_E"; ...
    "Outer correlation Tbag_M,Tbag_E"; ...
    "Approx MC error of difference SE"; ...
    "Number of valid outer replicates"; ...
    "Outer elapsed seconds"];
diagnosticValues = [outerSkewness; outerKurtosis; ...
    outerBimodalityCoefficient; JB_p; double(JB_reject); ...
    double(stronglyNonNormal); outerCovariance_M_E; ...
    outerCorrelation_M_E; outerSEMonteCarloErrorApprox; ...
    sum(validOuter); outerElapsedSeconds];
diagnosticTable = table(diagnosticNames, diagnosticValues, ...
    'VariableNames', {'Diagnostic', 'Value'});

configurationTable = table(string(filePath), string(folderPath), rngSeed, ...
    nInner, nOuter, delta0, alpha, mainFigureBinWidth, true, ...
    'VariableNames', {'SourceDataFile', 'OutputFolder', 'RNG_seed', ...
    'N_inner', 'N_outer', 'Delta0_h', 'Alpha', ...
    'MainFigureBinWidth_h', 'PairedResampling'});

methodsNotes = { ...
    'Item', 'Definition'; ...
    'Estimator', 'Bootstrap-smoothed / bagged discrete hard-argmax peak'; ...
    'Pairing', 'M and E from each biological sample are resampled with one common row index at every bootstrap layer'; ...
    'Peak', 'Sampled time attaining the maximum of the bootstrap mean curve'; ...
    'Bagged peak', 'Mean of inner-bootstrap discrete hard peak times'; ...
    'Uncertainty', 'SD and percentile CI across outer-bootstrap bagged estimates'; ...
    'MET null', 'H0: |Delta| <= delta0'; ...
    'MET alternative', 'H1: |Delta| > delta0'; ...
    'MET p-value', 'Wald normal approximation at the least-favorable interval-null boundary'; ...
    'Bootstrap support', 'Descriptive only; not the reported frequentist p-value'; ...
    'Inner Monte Carlo SE', 'Numerical-integration precision only; not biological sampling uncertainty'; ...
    'Main figure', 'Outer-bootstrap sampling distributions of Tbag_M and Tbag_E; display bin width is recorded in Configuration and is not raw temporal sampling resolution'};

% Each sheet is overwritten on rerun; unrelated files are untouched.
writetable(summaryTable, outputExcel, 'Sheet', 'Summary', ...
    'WriteMode', 'overwritesheet');
writetable(fullSummaryTable, outputExcel, 'Sheet', 'FullSummary', ...
    'WriteMode', 'overwritesheet');
writetable(conditionDiagnosticsTable, outputExcel, ...
    'Sheet', 'ConditionDiagnostics', 'WriteMode', 'overwritesheet');
writetable(peakProbabilityTable, outputExcel, ...
    'Sheet', 'PeakProbabilities', 'WriteMode', 'overwritesheet');
writetable(outerBootstrapTable, outputExcel, ...
    'Sheet', 'OuterBootstrap', 'WriteMode', 'overwritesheet');
writetable(timepointTable, outputExcel, ...
    'Sheet', 'TimepointData', 'WriteMode', 'overwritesheet');
writetable(diagnosticTable, outputExcel, ...
    'Sheet', 'Diagnostics', 'WriteMode', 'overwritesheet');
writetable(configurationTable, outputExcel, ...
    'Sheet', 'Configuration', 'WriteMode', 'overwritesheet');
writecell(methodsNotes, outputExcel, 'Sheet', 'MethodsNotes', ...
    'WriteMode', 'overwritesheet');

fprintf('Excel results saved to:\n  %s\n', outputExcel);

%% ========================================================================
% PART 3. Main figure: use the per-figure style configuration and save method
% ========================================================================
% The main figure now shows sampling uncertainty of the bagged estimator:
% each value is a complete inner-bagged peak calculated in one paired outer
% bootstrap dataset. This is distinct from innerPeakTimes_M/E, which are
% discrete hard-argmax locations used to construct one Tbag estimate.
outerPeak_M = outerTbag_M(validOuter);
outerPeak_E = outerTbag_E(validOuter);

% Display-bin width comes from the figure configuration. It does not
% change the observed time grid or the statistical estimator.
outerPeakAll = [outerPeak_M; outerPeak_E];
lowestBinCenter = floor(min(outerPeakAll) / mainFigureBinWidth) * ...
    mainFigureBinWidth;
highestBinCenter = ceil(max(outerPeakAll) / mainFigureBinWidth) * ...
    mainFigureBinWidth;
mainPeakEdges = (lowestBinCenter - mainFigureBinWidth / 2) : ...
    mainFigureBinWidth : ...
    (highestBinCenter + mainFigureBinWidth / 2);

% Ensure at least one complete bin even for a numerically degenerate outer
% distribution in which all replicates give the same Tbag value.
if numel(mainPeakEdges) < 2
    mainPeakEdges = lowestBinCenter + ...
        [-mainFigureBinWidth / 2, mainFigureBinWidth / 2];
end

[outerBinProb_M, ~] = histcounts(outerPeak_M, mainPeakEdges, ...
    'Normalization', 'probability');
[outerBinProb_E, ~] = histcounts(outerPeak_E, mainPeakEdges, ...
    'Normalization', 'probability');
maxOuterBinProbability = max([outerBinProb_M(:); outerBinProb_E(:)]);

M_color = [255/255, 128/255, 0];
E_color = [0, 0, 1];

fig = figure('Name', 'Paired outer-bootstrap bagged peak distribution', ...
    'NumberTitle', 'off', 'Units', 'centimeters', ...
    'Position', [2, 2, cfg.plot.figureSizeCm], 'Color', 'w');
hold on;

histogram(outerPeak_E, mainPeakEdges, 'FaceColor', E_color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.8, 'DisplayName', 'Ecell', ...
    'Normalization', 'probability');
histogram(outerPeak_M, mainPeakEdges, 'FaceColor', M_color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.8, 'DisplayName', 'Mcell', ...
    'Normalization', 'probability');

ax = gca;
ax.Units = 'centimeters';
ax.Position(3:4) = cfg.plot.axesSizeCm;
ax.LineWidth = cfg.plot.lineWidth;
ax.TickDir = 'out';
ax.TickLength = [0.02, 0.02];
ax.XTickLabelRotation = 0;
ax.XRuler.TickLabelGapOffset = -1;
ax.FontName = 'Arial';
ax.FontSize = cfg.plot.fontSize;
ax.Color = 'none';

% Preserve the original range whenever possible, but avoid clipping a bin
% probability above 0.8 when the outer distribution is highly concentrated.
if ~cfg.plot.autoYLimits || maxOuterBinProbability > 0.78
    ylim([0, 1]);
    yticks(0:0.5:1);
else
    ylim([0, 0.8]);
    yticks(0:0.4:0.8);
end
xlim(cfg.plot.xLimits);
xticks(cfg.plot.xTicks);

if cfg.plot.minorTicks
    ax.XMinorTick = 'on';
    ax.XAxis.MinorTickValues = setdiff(10:1:20, cfg.plot.xTicks);
end

box off;
hold off;

set(findall(fig, '-property', 'FontUnits'), 'FontUnits', 'points');
set(findall(fig, '-property', 'FontSize'), 'FontSize', cfg.plot.fontSize);
set(findall(fig, '-property', 'FontSizeMode'), 'FontSizeMode', 'manual');
set(findall(fig, '-property', 'FontName'), 'FontName', 'Arial');

fig.PaperUnits = 'centimeters';
fig.PaperSize = cfg.plot.figureSizeCm;
fig.PaperPosition = [0, 0, cfg.plot.figureSizeCm];
fig.PaperPositionMode = 'manual';
fig.InvertHardcopy = 'off';
drawnow;
print(fig, outputMainSVG, '-dsvg', '-vector');
fprintf('Main figure saved to:\n  %s\n', outputMainSVG);

%% ========================================================================
% PART 4. Numerical diagnostics and diagnostic figure
% ========================================================================
diagnosticFig = figure('Name', 'Paired bagged-peak diagnostics', ...
    'NumberTitle', 'off', 'Color', 'w', 'Units', 'centimeters', ...
    'Position', [3 3 18 14]);
diagnosticLayout = tiledlayout(diagnosticFig, 2, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% Panel 1: outer difference distribution, threshold, estimate, and normal fit.
nexttile(diagnosticLayout, 1);
histogram(outerBaggedDiff(validOuter), 30, 'Normalization', 'pdf', ...
    'FaceColor', [0.45 0.45 0.45], 'EdgeColor', 'none', ...
    'FaceAlpha', 0.75, 'DisplayName', 'Outer bootstrap');
hold on;
xGrid = linspace(min(outerBaggedDiff(validOuter)), ...
    max(outerBaggedDiff(validOuter)), 300);
if baggedSE > 0
    plot(xGrid, normpdf(xGrid, mean(outerBaggedDiff(validOuter)), ...
        baggedSE), 'Color', M_color, 'LineWidth', 1.2, ...
        'DisplayName', 'Normal fit');
end
xline(-delta0, '--k', '-delta0', 'LineWidth', 0.9, ...
    'HandleVisibility', 'off');
xline(delta0, '--k', '+delta0', 'LineWidth', 0.9, ...
    'HandleVisibility', 'off');
xline(baggedDiff, '-', 'Observed diff', 'Color', [0.75 0 0], ...
    'LineWidth', 1.2, 'HandleVisibility', 'off');
hold off;
xlabel('Bagged peak difference M-E (h)');
ylabel('Density');
title('Outer-bootstrap distribution');
legend('Location', 'best');
box off;

% Panel 2: Q-Q plot for the Wald normal approximation.
nexttile(diagnosticLayout, 2);
qqplot(outerBaggedDiff(validOuter));
xlabel('Standard-normal quantiles');
ylabel('Outer-difference quantiles (h)');
title(sprintf('Normal Q-Q plot (JB p=%.3g)', JB_p));
box off;

% Panel 3: paired outer bagged estimates and their covariance.
nexttile(diagnosticLayout, 3);
scatter(outerTbag_M(validOuter), outerTbag_E(validOuter), 10, ...
    [0.25 0.25 0.25], 'filled', 'MarkerFaceAlpha', 0.25);
xline(Tbag_M, '--', 'M estimate', 'Color', M_color, ...
    'LineWidth', 0.9);
yline(Tbag_E, '--', 'E estimate', 'Color', E_color, ...
    'LineWidth', 0.9);
xlabel('Outer Tbag M (h)');
ylabel('Outer Tbag E (h)');
title(sprintf('Paired outer estimates (r=%.3f)', ...
    outerCorrelation_M_E));
box off;

% Panel 4: cumulative Monte Carlo convergence across outer replicates.
nexttile(diagnosticLayout, 4);
outerDiffValid = outerBaggedDiff(validOuter);
nValidOuter = numel(outerDiffValid);
runningMean = cumsum(outerDiffValid) ./ (1:nValidOuter).';
runningSD = NaN(nValidOuter, 1);
for k = 2:nValidOuter
    runningSD(k) = std(outerDiffValid(1:k));
end

yyaxis left;
plot(1:nValidOuter, runningMean, 'Color', [0.2 0.2 0.2], ...
    'LineWidth', 0.9);
yline(baggedDiff, '--', 'Observed diff', 'Color', [0.75 0 0], ...
    'LineWidth', 0.8);
ylabel('Cumulative mean difference (h)');

yyaxis right;
plot(1:nValidOuter, runningSD, 'Color', M_color, 'LineWidth', 0.9);
yline(baggedSE, '--', 'Final SE', 'Color', M_color, ...
    'LineWidth', 0.8);
ylabel('Cumulative bootstrap SE (h)');
xlabel('Number of outer replicates');
title('Outer Monte Carlo convergence');
box off;

set(findall(diagnosticFig, '-property', 'FontName'), 'FontName', 'Arial');
set(findall(diagnosticFig, '-property', 'FontSize'), 'FontSize', 7);
diagnosticFig.PaperUnits = 'centimeters';
diagnosticFig.PaperSize = [18, 14];
diagnosticFig.PaperPosition = [0, 0, 18, 14];
diagnosticFig.PaperPositionMode = 'manual';
diagnosticFig.InvertHardcopy = 'off';
drawnow;
print(diagnosticFig, outputDiagnosticSVG, '-dsvg', '-vector');
fprintf('Diagnostic figure saved to:\n  %s\n', outputDiagnosticSVG);

%% ========================================================================
% Local helper functions
% ========================================================================

publicResult = struct('statistics', statsBaggedPaired, 'summary', fullSummaryTable, 'outerBootstrap', outerBootstrapTable, 'timepointData', timepointTable);
publicResult.configuration = cfg;
if exist('outputMAT', 'var') && isfile(outputMAT)
    save(outputMAT, 'cfg', 'publicResult', '-append');
else
    save(fullfile(cfg.outputDir, cfg.outputStem + "_results.mat"), 'cfg', 'publicResult');
end
end
