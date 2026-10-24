function peakProb = computePeakProbability(peakTimes, time)
    peakProb = NaN(1, numel(time));
    validPeak = isfinite(peakTimes);
    if ~any(validPeak)
        return;
    end

    for j = 1:numel(time)
        peakProb(j) = mean(peakTimes(validPeak) == time(j));
    end
end
