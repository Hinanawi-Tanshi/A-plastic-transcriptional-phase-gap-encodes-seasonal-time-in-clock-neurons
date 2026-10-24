function renamefolder(cellidfile,signalFolder)
%RENAMEFOLDER Rename once; reuse directories that already have the right name.
if nargin<2, signalFolder='tracked_signal'; end
data=load(cellidfile,'cellids'); cellids=data.cellids;
[folders,indices]=listSignalFolders(fullfile(pwd,signalFolder));
assert(size(cellids,1)==numel(indices),'MS2:CellIdentity','Cell_ID.mat must describe every cell.');
targets=cell(size(indices));
for k=1:numel(indices)
    i=indices(k);
    assert(strcmp(cellids{i,1},sprintf('cell_%02d',i)), ...
        'MS2:CellIdentity','Cell_ID.mat does not match tracked IDs.');
    label=char(cellids{i,2});
    assert(~isempty(strtrim(label)) && isempty(regexp(label,'[<>:"/\\|?*]','once')), ...
        'MS2:CellIdentity','Cell type must be a valid, nonempty directory-name component.');
    targets{k}=[label '_' cellids{i,1}];
    if ~strcmp(folders(k).name,targets{k})
        assert(strcmp(folders(k).name,cellids{i,1}),'MS2:AlreadyRenamed', ...
            'Cell %d already has a different type prefix. Inspect the mapping.',i);
        assert(~isfolder(fullfile(folders(k).folder,targets{k})), ...
            'MS2:ExistingCellFolder','Target cell directory already exists.');
    end
end
for k=1:numel(indices)
    if ~strcmp(folders(k).name,targets{k})
        movefile(fullfile(folders(k).folder,folders(k).name),fullfile(folders(k).folder,targets{k}));
    end
end
end
