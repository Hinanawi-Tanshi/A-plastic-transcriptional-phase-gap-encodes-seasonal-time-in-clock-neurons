function data = read_clock_long_sheet(filePath, sheetName, cfg)
%READ_CLOCK_LONG_SHEET Read and validate a four-column long-format clock sheet.
%
% Expected columns:
%   1. ZT label (for example "ZT4")
%   2. sample identifier
%   3. signal value
%   4. numeric ZT time

arguments
    filePath (1,1) string
    sheetName
    cfg struct
end

if ~isfile(filePath)
    error("ClockAnalysis:MissingInput", "Input file not found: %s", filePath);
end

[raw,sheetName] = bootstrapReadCells(filePath,sheetName);
if size(raw, 2) < 4
    error("ClockAnalysis:BadShape", ...
        "Sheet %s in %s must contain at least four columns.", sheetName, filePath);
end

signal = local_numeric_column(raw(:, 3));
time = local_numeric_column(raw(:, 4));
valid = isfinite(signal) & isfinite(time);

if ~all(valid)
    dropped = sum(~valid);
    warning("ClockAnalysis:DroppedRows", ...
        "Dropped %d non-data/header rows from %s [%s].", dropped, filePath, sheetName);
end

signal = signal(valid);
time = time(valid);
ztLabel = string(raw(valid, 1));
sampleID = string(raw(valid, 2));

if isempty(signal)
    error("ClockAnalysis:EmptySheet", "No numeric observations found in %s [%s].", filePath, sheetName);
end
if any(signal < 0 | signal > 1)
    error("ClockAnalysis:SignalRange", ...
        "Signal values in %s [%s] must lie between 0 and 1.", filePath, sheetName);
end
if any(time < 0 | time > cfg.period)
    error("ClockAnalysis:TimeRange", ...
        "ZT values in %s [%s] must lie between 0 and %.4g.", filePath, sheetName, cfg.period);
end

denominator = local_infer_denominator(signal, cfg.maxDenominator, cfg.denominatorTolerance);
if isfinite(denominator)
    family = "binomial";
    trials = repmat(denominator, numel(signal), 1);
    successes = round(signal .* trials);
    reconstructed = successes ./ trials;
    if max(abs(reconstructed - signal)) > cfg.denominatorTolerance
        error("ClockAnalysis:DenominatorInference", ...
            "Signal values cannot be represented by the inferred denominator.");
    end
else
    family = "normal";
    trials = nan(size(signal));
    successes = nan(size(signal));
end

labelTime = nan(size(time));
for row = 1:numel(ztLabel)
    token = regexp(ztLabel(row), "[-+]?\d*\.?\d+", "match", "once");
    if ~isempty(token)
        labelTime(row) = str2double(token);
    end
end
labelMismatch = isfinite(labelTime) & abs(labelTime - time) > cfg.denominatorTolerance;
if any(labelMismatch)
    warning("ClockAnalysis:ZTLabelMismatch", ...
        "%d rows in %s [%s] have inconsistent text and numeric ZT values.", ...
        sum(labelMismatch), filePath, sheetName);
end

key = string(time) + "|" + sampleID;
if numel(unique(key)) < numel(key)
    warning("ClockAnalysis:DuplicateKey", ...
        "Duplicate time/sample keys were detected in %s [%s].", filePath, sheetName);
end

data = struct();
data.filePath = filePath;
data.sheetName = sheetName;
data.ztLabel = ztLabel;
data.sampleID = sampleID;
data.time = time;
data.phase = mod(time, cfg.period);
data.signal = signal;
data.family = family;
data.denominator = denominator;
data.trials = trials;
data.successes = successes;
data.n = numel(signal);
data.uniqueTimes = unique(time, "sorted");
data.nTimePoints = numel(data.uniqueTimes);
data.table = table(ztLabel, sampleID, signal, time, ...
    'VariableNames', {'ZTLabel', 'SampleID', 'Signal', 'ZT'});
end

function values = local_numeric_column(cells)
values = nan(size(cells));
for row = 1:numel(cells)
    value = cells{row};
    if isnumeric(value) && isscalar(value)
        values(row) = double(value);
    elseif islogical(value) && isscalar(value)
        values(row) = double(value);
    elseif ischar(value) || isstring(value)
        parsed = str2double(string(value));
        if isfinite(parsed)
            values(row) = parsed;
        end
    end
end
end

function denominator = local_infer_denominator(signal, maxDenominator, tolerance)
denominator = nan;
for candidate = 1:maxDenominator
    if max(abs(signal .* candidate - round(signal .* candidate))) <= tolerance
        denominator = candidate;
        return;
    end
end
end
