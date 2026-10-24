function [peak_M, peak_E] = computeHardPeaks(Mdata, Edata, time)
    mean_M = cellfun(@(v) mean(v, 'omitnan'), Mdata);
    mean_E = cellfun(@(v) mean(v, 'omitnan'), Edata);
    [max_M, index_M] = max(mean_M, [], 'omitnan');
    [max_E, index_E] = max(mean_E, [], 'omitnan');
    if isempty(index_M) || ~isfinite(max_M)
        peak_M = NaN;
    else
        peak_M = time(index_M);
    end
    if isempty(index_E) || ~isfinite(max_E)
        peak_E = NaN;
    else
        peak_E = time(index_E);
    end
end
