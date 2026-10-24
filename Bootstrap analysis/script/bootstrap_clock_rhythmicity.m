function bootstrap = bootstrap_clock_rhythmicity(data, cfg, observedFit)
%BOOTSTRAP_CLOCK_RHYTHMICITY Calibrate periodic-vs-flat deviance under H0.
%
% cfg.nullMode:
%   "pooled_bootstrap" - resample responses from the pooled time-independent
%                        empirical distribution (default and robust)
%   "permutation"      - permute responses across fixed time labels
%   "parametric"       - simulate from the fitted flat GLM

arguments
    data struct
    cfg struct
    observedFit struct
end

targetCount = cfg.nNullBootstrap;
statistics = nan(targetCount, 1);
accepted = 0;
attempts = 0;
maxAttempts = max(targetCount * 3, targetCount + 100);

warningState = warning;
cleanup = onCleanup(@() warning(warningState)); %#ok<NASGU>
warning("off", "all");

while accepted < targetCount && attempts < maxAttempts
    attempts = attempts + 1;
    simulated = data;
    simulated.signal = local_null_response(data, cfg, observedFit);
    if data.family == "binomial"
        simulated.successes = round(simulated.signal .* simulated.trials);
    end

    try
        fitted = fit_clock_model(simulated, cfg);
        if isfinite(fitted.observedStatistic)
            accepted = accepted + 1;
            statistics(accepted) = fitted.observedStatistic;
        end
    catch
        % Retry until the requested number of valid bootstrap fits is reached.
    end

    if cfg.verbose && mod(accepted, cfg.progressEvery) == 0 && accepted > 0
        fprintf("  Null bootstrap: %d/%d valid fits\n", accepted, targetCount);
    end
end

statistics = statistics(1:accepted);
if accepted < targetCount
    error("ClockAnalysis:BootstrapFailure", ...
        "Only %d/%d valid null-bootstrap fits were obtained.", accepted, targetCount);
end

pValue = (1 + sum(statistics >= observedFit.observedStatistic)) / (accepted + 1);

bootstrap = struct();
bootstrap.mode = string(cfg.nullMode);
bootstrap.statistics = statistics;
bootstrap.observedStatistic = observedFit.observedStatistic;
bootstrap.pValue = pValue;
bootstrap.nullCI95 = prctile(statistics, [2.5, 97.5]);
bootstrap.nValid = accepted;
bootstrap.nAttempts = attempts;
end

function response = local_null_response(data, cfg, observedFit)
mode = lower(string(cfg.nullMode));
n = data.n;

switch mode
    case "pooled_bootstrap"
        if data.family == "binomial"
            response = nan(n, 1);
            trialLevels = unique(data.trials);
            for trialIndex = 1:numel(trialLevels)
                target = find(data.trials == trialLevels(trialIndex));
                pool = target;
                draw = pool(randi(numel(pool), numel(target), 1));
                response(target) = data.signal(draw);
            end
        else
            response = data.signal(randi(n, n, 1));
        end

    case "permutation"
        if data.family == "binomial"
            response = nan(n, 1);
            trialLevels = unique(data.trials);
            for trialIndex = 1:numel(trialLevels)
                target = find(data.trials == trialLevels(trialIndex));
                response(target) = data.signal(target(randperm(numel(target))));
            end
        else
            response = data.signal(randperm(n));
        end

    case "parametric"
        if data.family == "binomial"
            probability = observedFit.flatPredicted(1);
            response = zeros(n, 1);
            for row = 1:n
                successes = sum(rand(data.trials(row), 1) < probability);
                response(row) = successes / data.trials(row);
            end
        else
            meanLevel = observedFit.flatPredicted(1);
            residualScale = std(data.signal - meanLevel, 0);
            response = meanLevel + residualScale .* randn(n, 1);
            response = min(max(response, 0), 1);
        end

    otherwise
        error("ClockAnalysis:UnknownNullMode", ...
            "Unknown nullMode '%s'.", cfg.nullMode);
end
end
