function [Tbag_M, Tbag_E, peakTimes_M, peakTimes_E, P_M, P_E] = ...
    computePairedBaggedPeaks(Mdata, Edata, time, nInner, peakMode)
% Same row-index matrix is used for M and E at every time point and every
% inner replicate. Each replicate remains a hard discrete argmax.

    if nargin < 5, peakMode = "standard"; end
    nTime = numel(time);
    bootstrapMeans_M = NaN(nInner, nTime);
    bootstrapMeans_E = NaN(nInner, nTime);

    for j = 1:nTime
        Mj = Mdata{j};
        Ej = Edata{j};
        Mj = Mj(:);
        Ej = Ej(:);

        if numel(Mj) ~= numel(Ej)
            error('bootstrap:PairLengthMismatch', ...
                'M and E pair counts differ at time index %d.', j);
        end

        nPair = numel(Mj);
        if nPair > 0
            sampleIndex = randi(nPair, nPair, nInner);
            bootstrapMeans_M(:, j) = mean(Mj(sampleIndex), 1).';
            bootstrapMeans_E(:, j) = mean(Ej(sampleIndex), 1).';
        end
    end

    if peakMode == "finite"
        [peakValues_M, peakIndex_M] = max(bootstrapMeans_M, [], 2, 'omitnan');
        [peakValues_E, peakIndex_E] = max(bootstrapMeans_E, [], 2, 'omitnan');
    else
        [peakValues_M, peakIndex_M] = max(bootstrapMeans_M, [], 2);
        [peakValues_E, peakIndex_E] = max(bootstrapMeans_E, [], 2);
    end
    peakTimes_M = time(peakIndex_M);
    peakTimes_E = time(peakIndex_E);
    peakTimes_M = peakTimes_M(:);
    peakTimes_E = peakTimes_E(:);
    if peakMode == "finite"
        peakTimes_M(~isfinite(peakValues_M)) = NaN;
        peakTimes_E(~isfinite(peakValues_E)) = NaN;
    else
        peakTimes_M(isnan(peakValues_M)) = NaN;
        peakTimes_E(isnan(peakValues_E)) = NaN;
    end

    Tbag_M = mean(peakTimes_M, 'omitnan');
    Tbag_E = mean(peakTimes_E, 'omitnan');
    if nargout >= 5
        P_M = computePeakProbability(peakTimes_M, time);
        P_E = computePeakProbability(peakTimes_E, time);
    end
end
