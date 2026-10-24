function publicResult = runClockBootstrapAnalysis(cfg)
% Shared publication pipeline; settings are supplied by the figure entry script.
figureCleanup = bootstrapFigureContext(cfg); %#ok<NASGU>
perFile = cfg.perFile;
wtFile = cfg.wtFile;

if ~isfile(perFile) || ~isfile(wtFile)
    error("ClockAnalysis:MissingInput", ...
        "Expected per0_long.xlsx and Fig2F LD3DD1_long.xlsx beside this script.");
end
if ~isfolder(cfg.outputDir)
    mkdir(cfg.outputDir);
end

rng(cfg.randomSeed, "twister");

%% ------------------------------ Analysis -------------------------------
results = struct();
results.per0.m = analyze_clock_dataset(perFile, "m", "per0", cfg);
results.per0.e = analyze_clock_dataset(perFile, "e", "per0", cfg);
results.wt.m = analyze_clock_dataset(wtFile, "m", "WT", cfg);
results.wt.e = analyze_clock_dataset(wtFile, "e", "WT", cfg);

signalNames = ["m", "e"];
for signalIndex = 1:numel(signalNames)
    signalName = signalNames(signalIndex);
    fieldName = char(signalName);
    absoluteMargin = cfg.equivalenceAbsolute.(fieldName);
    if isfinite(absoluteMargin)
        delta = absoluteMargin;
        deltaSource = "absolute";
    else
        delta = cfg.equivalenceFractionOfWT * results.wt.(fieldName).fit.amplitude;
        deltaSource = "fraction_of_WT";
    end

    results.per0.(fieldName).equivalenceDelta = delta;
    results.per0.(fieldName).equivalenceDeltaSource = deltaSource;
    results.per0.(fieldName).rhythmicDetected = ...
        results.per0.(fieldName).nullBootstrap.pValue < cfg.alpha;
    results.per0.(fieldName).equivalenceSupported = ...
        results.per0.(fieldName).amplitudeBootstrap.amplitudeUpper95 < delta;

    if results.per0.(fieldName).rhythmicDetected && ...
            ~results.per0.(fieldName).equivalenceSupported
        classification = "rhythmic";
    elseif ~results.per0.(fieldName).rhythmicDetected && ...
            results.per0.(fieldName).equivalenceSupported
        classification = "practically_arrhythmic";
    elseif results.per0.(fieldName).rhythmicDetected && ...
            results.per0.(fieldName).equivalenceSupported
        classification = "detectable_but_below_margin";
    else
        classification = "inconclusive";
    end
    results.per0.(fieldName).classification = classification;

    results.wt.(fieldName).equivalenceDelta = delta;
    results.wt.(fieldName).rhythmicDetected = ...
        results.wt.(fieldName).nullBootstrap.pValue < cfg.alpha;
    results.wt.(fieldName).equivalenceSupported = false;
    if results.wt.(fieldName).rhythmicDetected
        results.wt.(fieldName).classification = "rhythmic";
    else
        results.wt.(fieldName).classification = "rhythmicity_not_detected";
    end
end

%% --------------------------- Summary table -----------------------------
genotypeOrder = ["per0", "wt"];
rowCount = numel(genotypeOrder) * numel(signalNames);
Genotype = strings(rowCount, 1);
Signal = strings(rowCount, 1);
Model = strings(rowCount, 1);
N = zeros(rowCount, 1);
NTimePoints = zeros(rowCount, 1);
Denominator = zeros(rowCount, 1);
NullBootstrapMode = strings(rowCount, 1);
NullBootstrapRequested = zeros(rowCount, 1);
NullBootstrapAttempts = zeros(rowCount, 1);
NullBootstrapValid = zeros(rowCount, 1);
NullStatisticMean = zeros(rowCount, 1);
NullStatisticSD = zeros(rowCount, 1);
NullStatisticCI95Low = zeros(rowCount, 1);
NullStatisticCI95High = zeros(rowCount, 1);
ObservedStatistic = zeros(rowCount, 1);
RhythmicityP = zeros(rowCount, 1);
FittedAmplitude = zeros(rowCount, 1);
AmplitudeBootstrapRequested = zeros(rowCount, 1);
AmplitudeBootstrapAttempts = zeros(rowCount, 1);
AmplitudeBootstrapValid = zeros(rowCount, 1);
AmplitudeBootstrapMean = zeros(rowCount, 1);
AmplitudeBootstrapSD = zeros(rowCount, 1);
AmplitudeBootstrapMonteCarloSEM = zeros(rowCount, 1);
AmplitudeMedian = zeros(rowCount, 1);
AmplitudeCI95Low = zeros(rowCount, 1);
AmplitudeCI95High = zeros(rowCount, 1);
AmplitudeUpper95 = zeros(rowCount, 1);
EquivalenceDelta = zeros(rowCount, 1);
RhythmicDetected = false(rowCount, 1);
EquivalenceSupported = false(rowCount, 1);
Classification = strings(rowCount, 1);

row = 0;
for genotypeIndex = 1:numel(genotypeOrder)
    genotypeField = char(genotypeOrder(genotypeIndex));
    for signalIndex = 1:numel(signalNames)
        signalField = char(signalNames(signalIndex));
        row = row + 1;
        current = results.(genotypeField).(signalField);
        Genotype(row) = current.genotype;
        Signal(row) = upper(current.signal);
        Model(row) = current.fit.modelName;
        N(row) = current.data.n;
        NTimePoints(row) = current.data.nTimePoints;
        Denominator(row) = current.data.denominator;
        NullBootstrapMode(row) = current.nullBootstrap.mode;
        NullBootstrapRequested(row) = cfg.nNullBootstrap;
        NullBootstrapAttempts(row) = current.nullBootstrap.nAttempts;
        NullBootstrapValid(row) = current.nullBootstrap.nValid;
        NullStatisticMean(row) = mean(current.nullBootstrap.statistics);
        NullStatisticSD(row) = std(current.nullBootstrap.statistics, 0);
        NullStatisticCI95Low(row) = current.nullBootstrap.nullCI95(1);
        NullStatisticCI95High(row) = current.nullBootstrap.nullCI95(2);
        ObservedStatistic(row) = current.fit.observedStatistic;
        RhythmicityP(row) = current.nullBootstrap.pValue;
        FittedAmplitude(row) = current.fit.amplitude;
        AmplitudeBootstrapRequested(row) = cfg.nAmplitudeBootstrap;
        AmplitudeBootstrapAttempts(row) = current.amplitudeBootstrap.nAttempts;
        AmplitudeBootstrapValid(row) = current.amplitudeBootstrap.nValid;
        amplitudeValues = current.amplitudeBootstrap.amplitude;
        AmplitudeBootstrapMean(row) = mean(amplitudeValues);
        AmplitudeBootstrapSD(row) = std(amplitudeValues, 0);
        AmplitudeBootstrapMonteCarloSEM(row) = ...
            AmplitudeBootstrapSD(row) / sqrt(AmplitudeBootstrapValid(row));
        AmplitudeMedian(row) = current.amplitudeBootstrap.amplitudeMedian;
        AmplitudeCI95Low(row) = current.amplitudeBootstrap.amplitudeCI95(1);
        AmplitudeCI95High(row) = current.amplitudeBootstrap.amplitudeCI95(2);
        AmplitudeUpper95(row) = current.amplitudeBootstrap.amplitudeUpper95;
        EquivalenceDelta(row) = current.equivalenceDelta;
        RhythmicDetected(row) = current.rhythmicDetected;
        EquivalenceSupported(row) = current.equivalenceSupported;
        Classification(row) = current.classification;
    end
end

summaryTable = table(Genotype, Signal, Model, N, NTimePoints, Denominator, ...
    NullBootstrapMode, NullBootstrapRequested, NullBootstrapAttempts, ...
    NullBootstrapValid, NullStatisticMean, NullStatisticSD, ...
    NullStatisticCI95Low, NullStatisticCI95High, ObservedStatistic, ...
    RhythmicityP, FittedAmplitude, AmplitudeBootstrapRequested, ...
    AmplitudeBootstrapAttempts, AmplitudeBootstrapValid, ...
    AmplitudeBootstrapMean, AmplitudeBootstrapSD, ...
    AmplitudeBootstrapMonteCarloSEM, AmplitudeMedian, ...
    AmplitudeCI95Low, AmplitudeCI95High, AmplitudeUpper95, ...
    EquivalenceDelta, RhythmicDetected, EquivalenceSupported, Classification);

% Descriptive statistics for the original observations at every ZT.
% N is the number of original rows at that ZT, not the number of bootstrap
% replicates. SEM is calculated as SD/sqrt(N).
descriptiveParts = cell(rowCount, 1);
row = 0;
for genotypeIndex = 1:numel(genotypeOrder)
    genotypeField = char(genotypeOrder(genotypeIndex));
    for signalIndex = 1:numel(signalNames)
        signalField = char(signalNames(signalIndex));
        row = row + 1;
        current = results.(genotypeField).(signalField);
        currentSummary = summarize_clock_timecourse(current.data);
        nSummaryRows = height(currentSummary);
        currentSummary = addvars(currentSummary, ...
            repmat(current.genotype, nSummaryRows, 1), ...
            repmat(upper(current.signal), nSummaryRows, 1), ...
            'Before', 1, 'NewVariableNames', {'Genotype', 'Signal'});
        descriptiveParts{row} = currentSummary;
    end
end
descriptiveStatisticsTable = vertcat(descriptiveParts{:});

modelTag = regexprep(char(cfg.modelName), "[^A-Za-z0-9]+", "_");
summaryPath = fullfile(cfg.outputDir, ...
    sprintf("%s_summary_%s.csv", cfg.outputStem, modelTag));
writetable(summaryTable, summaryPath);
descriptiveStatisticsPath = fullfile(cfg.outputDir, ...
    sprintf("%s_descriptive_statistics_%s.csv", cfg.outputStem, modelTag));
writetable(descriptiveStatisticsTable, descriptiveStatisticsPath);
resultPath = fullfile(cfg.outputDir, ...
    sprintf("%s_results_%s.mat", cfg.outputStem, modelTag));
save(resultPath, "cfg", "results", "summaryTable", ...
    "descriptiveStatisticsTable");

%% ----------------------------- Main figure -----------------------------
colors = struct();
colors.wt = [0.10, 0.35, 0.75];
colors.per0 = [0.90, 0.35, 0.10];

mainFigure = figure("Name", cfg.figureLabel + " rhythmicity analysis", ...
    "Color", "white", "Units", "centimeters", ...
    "Position", [2, 2, 18.5, 11.5]);
layout = tiledlayout(mainFigure, 2, 3, ...
    "TileSpacing", "compact", "Padding", "compact");

for signalIndex = 1:numel(signalNames)
    signalField = char(signalNames(signalIndex));

    timeAxis = nexttile(layout, (signalIndex - 1) * 3 + 1);
    plot_clock_timecourse_panel(timeAxis, results.wt.(signalField), ...
        results.per0.(signalField), cfg, colors, signalIndex == 1);

    nullAxis = nexttile(layout, (signalIndex - 1) * 3 + 2);
    plot_null_bootstrap_panel(nullAxis, results.per0.(signalField), colors.per0);

    amplitudeAxis = nexttile(layout, (signalIndex - 1) * 3 + 3);
    plot_amplitude_bootstrap_panel(amplitudeAxis, results.wt.(signalField), ...
        results.per0.(signalField), results.per0.(signalField).equivalenceDelta, colors);
end

title(layout, sprintf("Flexible-shape rhythmicity analysis (%s; B_{null}=%d, B_{amp}=%d)", ...
    strrep(cfg.modelName, "_", " "), cfg.nNullBootstrap, cfg.nAmplitudeBootstrap), ...
    "FontName", "Arial", "FontSize", 9, "FontWeight", "normal");

local_export_figure(mainFigure, cfg.outputDir, ...
    sprintf("%s_%s", cfg.outputStem, modelTag));

%% -------------------- All null distributions figure -------------------
nullFigure = figure("Name", "All rhythmicity null distributions", ...
    "Color", "white", "Units", "centimeters", ...
    "Position", [2, 2, 12.5, 10.5]);
nullLayout = tiledlayout(nullFigure, 2, 2, ...
    "TileSpacing", "compact", "Padding", "compact");

plot_null_bootstrap_panel(nexttile(nullLayout, 1), results.per0.m, colors.per0);
plot_null_bootstrap_panel(nexttile(nullLayout, 2), results.per0.e, colors.per0);
plot_null_bootstrap_panel(nexttile(nullLayout, 3), results.wt.m, colors.wt);
plot_null_bootstrap_panel(nexttile(nullLayout, 4), results.wt.e, colors.wt);
title(nullLayout, "Periodic-model improvement under the no-rhythm null", ...
    "FontName", "Arial", "FontSize", 9, "FontWeight", "normal");

local_export_figure(nullFigure, cfg.outputDir, ...
    sprintf("%s_null_all_%s", cfg.outputStem, modelTag));

%% ------------------------------ Report ---------------------------------
disp(summaryTable);
disp(descriptiveStatisticsTable);
fprintf("\nSaved analysis to:\n  %s\n  %s\n  %s\n", ...
    resultPath, summaryPath, descriptiveStatisticsPath);

%% ---------------- Standalone amplitude-bootstrap boxplot ---------------
% Do not run a t-test on the bootstrap replicates: they are estimates
% generated from the same biological data, not independent observations.
% Compare genotype within the same cell type:
%   WT M versus per0 M, and WT E versus per0 E.
amplitudeComparisonM = bootstrap_amplitude_difference_test( ...
    results.wt.m.amplitudeBootstrap.amplitude, ...
    results.per0.m.amplitudeBootstrap.amplitude);
amplitudeComparisonE = bootstrap_amplitude_difference_test( ...
    results.wt.e.amplitudeBootstrap.amplitude, ...
    results.per0.e.amplitudeBootstrap.amplitude);

rawAmplitudeP = [amplitudeComparisonM.pValue; amplitudeComparisonE.pValue];
adjustedAmplitudeP = holm_adjust_pvalues(rawAmplitudeP);

results.amplitudeComparisons = struct();
results.amplitudeComparisons.m = amplitudeComparisonM;
results.amplitudeComparisons.e = amplitudeComparisonE;
results.amplitudeComparisons.rawP = rawAmplitudeP;
results.amplitudeComparisons.holmAdjustedP = adjustedAmplitudeP;

Comparison = ["WT M vs per0 M"; "WT E vs per0 E"];
RawBootstrapP = rawAmplitudeP;
HolmAdjustedP = adjustedAmplitudeP;
MedianDifference = [amplitudeComparisonM.medianDifference; ...
    amplitudeComparisonE.medianDifference];
DifferenceCI95Low = [amplitudeComparisonM.differenceCI95(1); ...
    amplitudeComparisonE.differenceCI95(1)];
DifferenceCI95High = [amplitudeComparisonM.differenceCI95(2); ...
    amplitudeComparisonE.differenceCI95(2)];
amplitudeComparisonTable = table(Comparison, RawBootstrapP, HolmAdjustedP, ...
    MedianDifference, DifferenceCI95Low, DifferenceCI95High);

amplitudeComparisonPath = fullfile(cfg.outputDir, ...
    sprintf("%s_amplitude_comparison_%s.csv", cfg.outputStem, modelTag));
writetable(amplitudeComparisonTable, amplitudeComparisonPath);

amplitudeFigure = figure("Name", "Amplitude bootstrap boxplot", ...
    "Color", "white", "Units", "centimeters", ...
    "Position", [2, 2, 7.2, 6.3]);
amplitudeAxis = axes(amplitudeFigure);
amplitudeAxis.Units = 'centimeters';
amplitudeAxis.Position(3:4) = [3.8, 2.5];
amplitudeAxis.FontName = "Arial";
amplitudeAxis.FontSize = 8;

% ax.Units = 'centimeters';
% ax.Position(3:4) = [1.4, 1.8];
% 
% 
% ax.LineWidth = 0.8;
% ax.TickDir = 'out';
% ax.TickLength = [0.02, 0.02];
% ax.XTickLabelRotation = 0;
% ax.XRuler.TickLabelGapOffset = -1;
% 
% ax.FontName = 'Arial';
% ax.FontSize = 7;
plot_amplitude_boxplot_panel(amplitudeAxis, results, adjustedAmplitudeP, ...
    colors, cfg.plotAmplitudeScatter);
local_export_figure(amplitudeFigure, cfg.outputDir, ...
    sprintf("%s_amplitude_boxplot_%s", cfg.outputStem, modelTag));

% Overwrite the MAT result once so it also contains the comparison results.
save(resultPath, "cfg", "results", "summaryTable", ...
    "descriptiveStatisticsTable", "amplitudeComparisonTable");
disp(amplitudeComparisonTable);

%% Local export utility (plotting/output only; analysis is in separate files)

publicResult = results;
publicResult.configuration = cfg;
if exist('outputMAT', 'var') && isfile(outputMAT)
    save(outputMAT, 'cfg', 'publicResult', '-append');
else
    save(fullfile(cfg.outputDir, cfg.outputStem + "_results.mat"), 'cfg', 'publicResult');
end
end

function local_export_figure(fig, outputDir, baseName)
savefig(fig, fullfile(outputDir, baseName + ".fig"));
exportgraphics(fig, fullfile(outputDir, baseName + ".png"), "Resolution", 600);
exportgraphics(fig, fullfile(outputDir, baseName + ".pdf"), "ContentType", "vector");
print(fig, fullfile(outputDir, baseName + ".svg"), "-dsvg");
end
