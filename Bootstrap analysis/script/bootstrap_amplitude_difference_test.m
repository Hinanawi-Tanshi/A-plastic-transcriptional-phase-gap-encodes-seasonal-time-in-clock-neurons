function comparison = bootstrap_amplitude_difference_test(referenceAmplitude, testAmplitude)
%BOOTSTRAP_AMPLITUDE_DIFFERENCE_TEST Compare two independent bootstrap distributions.
%
% This is deliberately not a t-test on bootstrap replicates. The two
% independent bootstrap distributions are convolved to obtain the
% distribution of reference - test amplitude. The two-sided tail
% probability quantifies whether zero is compatible with that difference.

arguments
    referenceAmplitude (:,1) double
    testAmplitude (:,1) double
end

referenceAmplitude = referenceAmplitude(isfinite(referenceAmplitude));
testAmplitude = testAmplitude(isfinite(testAmplitude));
if isempty(referenceAmplitude) || isempty(testAmplitude)
    error("ClockAnalysis:EmptyBootstrap", ...
        "Amplitude comparison requires two non-empty bootstrap distributions.");
end

differences = referenceAmplitude - testAmplitude';
differences = differences(:);
nDifferences = numel(differences);

probabilityNonPositive = (1 + sum(differences <= 0)) / (nDifferences + 1);
probabilityNonNegative = (1 + sum(differences >= 0)) / (nDifferences + 1);
pValue = min(1, 2 * min(probabilityNonPositive, probabilityNonNegative));

comparison = struct();
comparison.direction = "reference_minus_test";
comparison.medianDifference = median(differences);
comparison.differenceCI95 = prctile(differences, [2.5, 97.5]);
comparison.pValue = pValue;
comparison.nPairwiseDifferences = nDifferences;
end
