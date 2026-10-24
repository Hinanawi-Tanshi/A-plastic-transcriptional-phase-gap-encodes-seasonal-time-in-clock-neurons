function fit = fit_clock_model(data, cfg)
%FIT_CLOCK_MODEL Fit flat and selected periodic GLMs and compare deviances.

arguments
    data struct
    cfg struct
end

[X, state] = clock_model_design(data.time, cfg);
designFlat = ones(data.n, 1);
designPeriodic = [ones(data.n, 1), X];

if data.family == "binomial"
    response = [data.successes, data.trials];
    [betaFlat, devFlat] = glmfit(designFlat, response, "binomial", ...
        "Link", "logit", "Constant", "off");
    [betaPeriodic, devPeriodic] = glmfit(designPeriodic, response, "binomial", ...
        "Link", "logit", "Constant", "off");
    inverseLink = "logit";
else
    [betaFlat, devFlat] = glmfit(designFlat, data.signal, "normal", ...
        "Link", "identity", "Constant", "off");
    [betaPeriodic, devPeriodic] = glmfit(designPeriodic, data.signal, "normal", ...
        "Link", "identity", "Constant", "off");
    inverseLink = "identity";
end

if lower(string(cfg.modelName)) == "categorical"
    grid = state.levels;
    if any(abs(grid) < 1e-10)
        grid = [grid; cfg.period];
    end
else
    grid = (0:cfg.gridStep:cfg.period)';
    if grid(end) < cfg.period
        grid = [grid; cfg.period];
    end
end

[XGrid, ~] = clock_model_design(grid, cfg, state);
predicted = glmval(betaPeriodic, [ones(numel(grid), 1), XGrid], ...
    inverseLink, "Constant", "off");
flatPredicted = glmval(betaFlat, ones(numel(grid), 1), ...
    inverseLink, "Constant", "off");

predicted = min(max(predicted, 0), 1);
flatPredicted = min(max(flatPredicted, 0), 1);

observedStatistic = max(0, devFlat - devPeriodic);
amplitude = max(predicted) - min(predicted);
meanLevel = trapz(grid, predicted) / cfg.period;
rmsAmplitude = sqrt(trapz(grid, (predicted - meanLevel).^2) / cfg.period);

fit = struct();
fit.modelName = string(cfg.modelName);
fit.modelState = state;
fit.betaFlat = betaFlat;
fit.betaPeriodic = betaPeriodic;
fit.devFlat = devFlat;
fit.devPeriodic = devPeriodic;
fit.observedStatistic = observedStatistic;
fit.timeGrid = grid;
fit.predicted = predicted;
fit.flatPredicted = flatPredicted;
fit.amplitude = amplitude;
fit.rmsAmplitude = rmsAmplitude;
fit.nPeriodicParameters = size(designPeriodic, 2);
fit.nFlatParameters = size(designFlat, 2);
fit.designRank = rank(designPeriodic);
end
