function [X, state] = clock_model_design(time, cfg, state)
%CLOCK_MODEL_DESIGN Construct the selected periodic-model design matrix.
%
% Supported cfg.modelName values:
%   "periodic_spline" - periodic cubic B-spline basis
%   "fourier"         - K-harmonic Fourier basis
%   "cosinor"         - one-harmonic Fourier basis
%   "categorical"     - one coefficient per observed phase (minus reference)

arguments
    time (:,1) double
    cfg struct
    state struct = struct()
end

phase = mod(time, cfg.period);
modelName = lower(string(cfg.modelName));
training = isempty(fieldnames(state));

switch modelName
    case "periodic_spline"
        if training
            state.modelName = modelName;
            state.nBasis = cfg.nSplineBasis;
        end
        basis = local_periodic_cubic_basis(phase, cfg.period, state.nBasis);
        % The full periodic basis sums to one. Dropping its last column
        % avoids collinearity with the intercept while preserving the span.
        X = basis(:, 1:end-1);

    case "fourier"
        if training
            state.modelName = modelName;
            state.harmonicOrder = cfg.harmonicOrder;
        end
        X = local_fourier_basis(phase, cfg.period, state.harmonicOrder);

    case "cosinor"
        if training
            state.modelName = modelName;
            state.harmonicOrder = 1;
        end
        X = local_fourier_basis(phase, cfg.period, 1);

    case "categorical"
        if training
            state.modelName = modelName;
            state.levels = unique(phase, "sorted");
            if numel(state.levels) < 2
                error("ClockAnalysis:TooFewTimes", ...
                    "The categorical model requires at least two distinct phases.");
            end
        end
        X = zeros(numel(phase), numel(state.levels) - 1);
        for levelIndex = 2:numel(state.levels)
            X(:, levelIndex - 1) = abs(phase - state.levels(levelIndex)) < 1e-10;
        end

    otherwise
        error("ClockAnalysis:UnknownModel", ...
            "Unknown modelName '%s'. Use periodic_spline, fourier, cosinor, or categorical.", ...
            cfg.modelName);
end

state.nPredictors = size(X, 2);
end

function X = local_fourier_basis(phase, period, harmonicOrder)
X = zeros(numel(phase), 2 * harmonicOrder);
omega = 2 * pi / period;
for harmonic = 1:harmonicOrder
    X(:, 2 * harmonic - 1) = cos(harmonic * omega .* phase);
    X(:, 2 * harmonic) = sin(harmonic * omega .* phase);
end
end

function basis = local_periodic_cubic_basis(phase, period, nBasis)
if nBasis < 4
    error("ClockAnalysis:SplineBasis", "nSplineBasis must be at least 4.");
end

u = mod(phase, period) ./ period .* nBasis;
basis = zeros(numel(u), nBasis);

for basisIndex = 0:(nBasis - 1)
    wrappedValue = zeros(size(u));
    for wrap = -2:2
        distance = u - basisIndex + wrap * nBasis;
        wrappedValue = wrappedValue + local_centered_cubic_bspline(distance);
    end
    basis(:, basisIndex + 1) = wrappedValue;
end

rowSum = sum(basis, 2);
if any(rowSum <= 0)
    error("ClockAnalysis:SplineBasis", "Periodic spline basis normalization failed.");
end
basis = basis ./ rowSum;
end

function value = local_centered_cubic_bspline(x)
absoluteX = abs(x);
value = zeros(size(x));

insideOne = absoluteX < 1;
value(insideOne) = 2 / 3 - absoluteX(insideOne).^2 + ...
    0.5 .* absoluteX(insideOne).^3;

betweenOneAndTwo = absoluteX >= 1 & absoluteX < 2;
value(betweenOneAndTwo) = (2 - absoluteX(betweenOneAndTwo)).^3 ./ 6;
end
