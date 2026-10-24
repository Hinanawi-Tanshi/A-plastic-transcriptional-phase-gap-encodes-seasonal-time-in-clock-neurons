function adjusted = holm_adjust_pvalues(pValues)
%HOLM_ADJUST_PVALUES Holm step-down family-wise error correction.

arguments
    pValues (:,1) double
end

if any(~isfinite(pValues) | pValues < 0 | pValues > 1)
    error("ClockAnalysis:BadPValue", "p-values must be finite and lie in [0, 1].");
end

[sortedP, order] = sort(pValues);
nTests = numel(pValues);
sortedAdjusted = zeros(size(sortedP));

for index = 1:nTests
    sortedAdjusted(index) = (nTests - index + 1) * sortedP(index);
end
sortedAdjusted = cummax(sortedAdjusted);
sortedAdjusted = min(sortedAdjusted, 1);

adjusted = zeros(size(pValues));
adjusted(order) = sortedAdjusted;
end
