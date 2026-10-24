function edges = makeCenteredEdges(values, binWidth)
% Histogram bins centered on integer multiples of binWidth.
    values = values(isfinite(values));
    if isempty(values) || ~isfinite(binWidth) || binWidth <= 0
        error('bootstrap:InvalidHistogramInput', ...
            'Histogram values and bin width must be valid.');
    end
    lowestCenter = floor(min(values) / binWidth) * binWidth;
    highestCenter = ceil(max(values) / binWidth) * binWidth;
    edges = (lowestCenter - binWidth / 2) : binWidth : ...
        (highestCenter + binWidth / 2);
    if numel(edges) < 2
        edges = lowestCenter + [-binWidth / 2, binWidth / 2];
    end
end
