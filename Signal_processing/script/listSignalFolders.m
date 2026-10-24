function [folders,indices,files]=listSignalFolders(rootpath,selection)
%LISTSIGNALFOLDERS Find actual per-cell result folders in numeric ID order.
% selection: vector of IDs, or a SigInfo struct with optional plotCells.
if nargin<2, selection=[]; end
if isstruct(selection)
    if isfield(selection,'plotCells'), selection=selection.plotCells;
    else, selection=[]; end
end
assert(isfolder(rootpath),'MS2:MissingSignals','Missing signal directory: %s',rootpath);
folders=dir(rootpath); folders=folders([folders.isdir]);
tokens=regexp({folders.name},'(?:^|_)cell_(\d+)$','tokens','once');
keep=~cellfun(@isempty,tokens); folders=folders(keep); tokens=tokens(keep);
indices=cellfun(@(token) str2double(token{1}),tokens);
assert(~isempty(indices),'MS2:MissingSignals','No cell result folders in %s.',rootpath);
assert(numel(unique(indices))==numel(indices),'MS2:DuplicateCell','Duplicate cell result folders.');
[indices,order]=sort(indices); folders=folders(order);
if ~isempty(selection)
    assert(all(ismember(selection,indices)),'MS2:MissingCell','A requested cell has no result folder.');
    keep=ismember(indices,selection); folders=folders(keep); indices=indices(keep);
end
files=cell(1,numel(indices));
for k=1:numel(indices)
    files{k}=fullfile(folders(k).folder,folders(k).name,sprintf('cell_%02d_tracked_signal.mat',indices(k)));
    assert(isfile(files{k}),'MS2:MissingCellData','Missing result file: %s',files{k});
end
end
