function maxcellnum=TrackedDataSave(SigInfo)
%TRACKEDDATASAVE Apply tracked masks to images matched by time and Z index.
% Accepts e.g. Den_t001_z001.tif, DenRaw_t0001_z0001.tif and
% track_results_t0001_z0001.tif. Requires a complete, unique T-by-Z grid.
rawFiles=indexPlanes(fullfile(SigInfo.p,SigInfo.rawfolder),SigInfo);
maskFiles=indexPlanes(fullfile(SigInfo.p,SigInfo.trackfolder),SigInfo);
target=fullfile(SigInfo.p,SigInfo.targetfolder);
expectedSize=[SigInfo.xvalue SigInfo.yvalue];
celllabel=zeros([expectedSize SigInfo.zvalue],'uint16');
frameIds=cell(1,SigInfo.tvalue);
% Validate masks and image dimensions before creating per-cell outputs.
for k=1:numel(maskFiles)
    mask=imread(maskFiles{k});
    assert(ismatrix(mask) && isequal(size(mask),expectedSize), ...
        'MS2:PlaneSize','Mask size differs from the configured size: %s',maskFiles{k});
    assert(isnumeric(mask) && all(isfinite(mask(:))) && all(mask(:)>=0) && ...
        all(mask(:)==fix(mask(:))) && max(mask(:))<=65535, ...
        'MS2:MaskLabels','Masks must contain nonnegative integer IDs, not RGB colors.');
    rawInfo=imfinfo(rawFiles{k});
    assert(numel(rawInfo)==1 && isequal([rawInfo.Height rawInfo.Width],expectedSize), ...
        'MS2:PlaneSize','Expected one plane with configured size: %s',rawFiles{k});
    t=floor((k-1)/SigInfo.zvalue)+1; z=mod(k-1,SigInfo.zvalue)+1;
    ids=double(unique(mask(mask>0)))';
    frameIds{t}=unique([frameIds{t}(:);ids(:)])';
    if t==1, celllabel(:,:,z)=uint16(mask); end
end
ids=frameIds{1};
assert(~isempty(ids) && isequal(ids,1:max(ids)),'MS2:CellLabels', ...
    'The first volume must use consecutive cell IDs 1:N (0 is background).');
maxcellnum=max(ids);
for t=1:SigInfo.tvalue
    assert(isequal(frameIds{t},ids),'MS2:MissingCell', ...
        'Cell IDs at time %d differ from the first volume. Check tracking first.',t);
end
if ~isfolder(target), mkdir(target); end
celllabelmax=max(celllabel,[],3);
f=figure('Position',[50 50 SigInfo.yvalue SigInfo.xvalue]);
imagesc(celllabelmax); axis image; hold on
for cellIndex=ids
    [r,c]=find(celllabelmax==cellIndex);
    text(mean(c),mean(r),num2str(cellIndex),'Color','white');
end
saveas(f,fullfile(SigInfo.p,'cell_labels.tif'));
for k=1:numel(rawFiles)
    raw=imread(rawFiles{k});
    assert(ismatrix(raw),'MS2:PlaneSize','Raw TIFFs must be grayscale planes.');
    mask=imread(maskFiles{k});
    t=floor((k-1)/SigInfo.zvalue)+1; z=mod(k-1,SigInfo.zvalue)+1;
    cellIds=double(unique(mask(mask>0)))';
    for cellIndex=cellIds
        cellImage=uint16(mask==cellIndex).*raw; % historical arithmetic
        file=fullfile(target,sprintf('cell%d_t%04d_z%04d.tif',cellIndex,t,z));
        imwrite(cellImage,file,'Compression','none');
    end
end
end

function files=indexPlanes(folder,info)
assert(isfolder(folder),'MS2:MissingImages','Missing image directory: %s',folder);
entries=dir(folder); entries=entries(~[entries.isdir]);
isTiff=~cellfun(@isempty,regexpi({entries.name},'\.tiff?$','once'));
entries=entries(isTiff);
files=cell(1,info.tvalue*info.zvalue);
for k=1:numel(entries)
    token=regexpi(entries(k).name,'(?:^|_)t(\d+)_z(\d+)\.tiff?$','tokens','once');
    assert(~isempty(token),'MS2:PlaneName','Cannot parse time/Z: %s',entries(k).name);
    t=str2double(token{1}); z=str2double(token{2});
    assert(t>=1 && t<=info.tvalue && z>=1 && z<=info.zvalue, ...
        'MS2:PlaneRange','Time/Z is outside configured range: %s',entries(k).name);
    index=(t-1)*info.zvalue+z;
    assert(isempty(files{index}),'MS2:DuplicatePlane','Duplicate time %d, Z %d in %s.',t,z,folder);
    files{index}=fullfile(entries(k).folder,entries(k).name);
end
assert(all(~cellfun(@isempty,files)),'MS2:MissingPlane','Incomplete time/Z grid in %s.',folder);
end
