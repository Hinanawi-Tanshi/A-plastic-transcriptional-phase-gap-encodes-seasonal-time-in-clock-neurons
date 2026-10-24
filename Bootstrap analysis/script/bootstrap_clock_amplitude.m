function bootstrap = bootstrap_clock_amplitude(data, cfg, observedFit)
%BOOTSTRAP_CLOCK_AMPLITUDE Stratified biological-sample bootstrap by ZT.

arguments
    data struct
    cfg struct
    observedFit struct
end

targetCount = cfg.nAmplitudeBootstrap;
amplitude = nan(targetCount, 1);
rmsAmplitude = nan(targetCount, 1);
predictions = nan(targetCount, numel(observedFit.timeGrid));
accepted = 0;
attempts = 0;
maxAttempts = max(targetCount * 3, targetCount + 100);
timeLevels = unique(data.time, "sorted");

warningState = warning;
cleanup = onCleanup(@() warning(warningState)); %#ok<NASGU>
warning("off", "all");

while accepted < targetCount && attempts < maxAttempts
    attempts = attempts + 1;
    selectedRows = [];
    for timeIndex = 1:numel(timeLevels)
        candidates = find(data.time == timeLevels(timeIndex));
        draw = candidates(randi(numel(candidates), numel(candidates), 1));
        selectedRows = [selectedRows; draw]; %#ok<AGROW>
    end

    sampled = local_subset_data(data, selectedRows);
    try
        fitted = fit_clock_model(sampled, cfg);
        if isfinite(fitted.amplitude) && numel(fitted.predicted) == size(predictions, 2)
            accepted = accepted + 1;
            amplitude(accepted) = fitted.amplitude;
            rmsAmplitude(accepted) = fitted.rmsAmplitude;
            predictions(accepted, :) = fitted.predicted(:)';
        end
    catch
        % Retry failed/separated fits.
    end

    if cfg.verbose && mod(accepted, cfg.progressEvery) == 0 && accepted > 0
        fprintf("  Amplitude bootstrap: %d/%d valid fits\n", accepted, targetCount);
    end
end

if accepted < targetCount
    error("ClockAnalysis:BootstrapFailure", ...
        "Only %d/%d valid amplitude-bootstrap fits were obtained.", accepted, targetCount);
end

amplitude = amplitude(1:accepted);
rmsAmplitude = rmsAmplitude(1:accepted);
predictions = predictions(1:accepted, :);

bootstrap = struct();
bootstrap.amplitude = amplitude;
bootstrap.rmsAmplitude = rmsAmplitude;
bootstrap.amplitudeMedian = median(amplitude);
bootstrap.amplitudeCI95 = prctile(amplitude, [2.5, 97.5]);
bootstrap.amplitudeUpper95 = prctile(amplitude, 95);
bootstrap.rmsMedian = median(rmsAmplitude);
bootstrap.rmsCI95 = prctile(rmsAmplitude, [2.5, 97.5]);
bootstrap.timeGrid = observedFit.timeGrid;
bootstrap.curveMedian = prctile(predictions, 50, 1)';
bootstrap.curveCI95Low = prctile(predictions, 2.5, 1)';
bootstrap.curveCI95High = prctile(predictions, 97.5, 1)';
bootstrap.nValid = accepted;
bootstrap.nAttempts = attempts;
end

function subset = local_subset_data(data, selectedRows)
subset = data;
fieldsToSubset = ["ztLabel", "sampleID", "time", "phase", "signal", ...
    "trials", "successes"];
for fieldIndex = 1:numel(fieldsToSubset)
    fieldName = fieldsToSubset(fieldIndex);
    subset.(fieldName) = data.(fieldName)(selectedRows, :);
end
subset.n = numel(selectedRows);
subset.uniqueTimes = unique(subset.time, "sorted");
subset.nTimePoints = numel(subset.uniqueTimes);
subset.table = data.table(selectedRows, :);
end
