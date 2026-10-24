function [raw,label] = bootstrapReadCells(filePath,selector)
% Read a whole sheet, or one original table's rectangular block in a sheet.
% Numeric sheet positions keep manual tab renaming independent of analysis.
if isstruct(selector)
    validateattributes(selector.sheet,{'numeric'},{'scalar','integer','positive'});
    raw = readcell(filePath,'Sheet',selector.sheet,'Range',selector.range);
    label = "Sheet" + string(selector.sheet) + ":" + string(selector.range);
else
    raw = readcell(filePath,'Sheet',selector);
    label = string(selector);
end
end
