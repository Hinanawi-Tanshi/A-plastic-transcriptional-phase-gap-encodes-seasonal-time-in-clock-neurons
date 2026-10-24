function raw = bootstrapReadMatrix(filePath,selector)
% Read a numeric source block without depending on the user-editable tab name.
if isstruct(selector)
    validateattributes(selector.sheet,{'numeric'},{'scalar','integer','positive'});
    raw = readmatrix(filePath,'Sheet',selector.sheet,'Range',selector.range);
else
    raw = readmatrix(filePath,'Sheet',selector);
end
end
