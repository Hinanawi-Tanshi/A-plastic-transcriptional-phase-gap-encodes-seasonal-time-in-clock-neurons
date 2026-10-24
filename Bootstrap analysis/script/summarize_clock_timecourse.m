function summary = summarize_clock_timecourse(data)
%SUMMARIZE_CLOCK_TIMECOURSE Compute descriptive statistics by ZT.

time = unique(data.time, "sorted");
n = zeros(size(time));
meanSignal = zeros(size(time));
standardDeviation = zeros(size(time));
standardError = zeros(size(time));
medianSignal = zeros(size(time));
minimumSignal = zeros(size(time));
maximumSignal = zeros(size(time));

for timeIndex = 1:numel(time)
    values = data.signal(data.time == time(timeIndex));
    n(timeIndex) = numel(values);
    meanSignal(timeIndex) = mean(values);
    standardDeviation(timeIndex) = std(values, 0);
    standardError(timeIndex) = standardDeviation(timeIndex) / sqrt(n(timeIndex));
    medianSignal(timeIndex) = median(values);
    minimumSignal(timeIndex) = min(values);
    maximumSignal(timeIndex) = max(values);
end

lower = max(0, meanSignal - 1.96 .* standardError);
upper = min(1, meanSignal + 1.96 .* standardError);

summary = table(time, n, meanSignal, standardDeviation, standardError, ...
    medianSignal, minimumSignal, maximumSignal, lower, upper, ...
    'VariableNames', {'ZT', 'N', 'Mean', 'SD', 'SEM', 'Median', ...
    'Minimum', 'Maximum', 'CI95Low', 'CI95High'});
end
