function [bootstrapData_M, bootstrapData_E] = ...
    bootstrapPairedDataset(Mdata, Edata)
% Form one outer dataset using one common within-time row index for M and E.

    nTime = numel(Mdata);
    bootstrapData_M = cell(size(Mdata));
    bootstrapData_E = cell(size(Edata));

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
        if nPair == 0
            bootstrapData_M{j} = Mj;
            bootstrapData_E{j} = Ej;
        else
            sampleIndex = randi(nPair, nPair, 1);
            bootstrapData_M{j} = Mj(sampleIndex);
            bootstrapData_E{j} = Ej(sampleIndex);
        end
    end
end
