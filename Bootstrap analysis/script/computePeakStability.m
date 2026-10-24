function [stability, entropyValue] = computePeakStability(peakProb)
% 1 = a point mass at one sampled time; 0 = uniform over all sampled times.
    K = numel(peakProb);
    if K <= 1
        entropyValue = 0;
        stability = 1;
        return;
    end

    p = peakProb(isfinite(peakProb) & peakProb > 0);
    if isempty(p)
        entropyValue = NaN;
        stability = NaN;
        return;
    end

    entropyValue = -sum(p .* log(p));
    stability = 1 - entropyValue / log(K);
end
