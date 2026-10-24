function Cells=DeNoise(SigInfo)
%DENOISE Preserve the original two-channel subtraction/masking procedure.
% The bundled ND2 reader supports the two-channel acquisition used here;
% its channel assignment depends on the ND2 padding style. Check new data.
p=SigInfo.p;
if isfield(SigInfo,'nd2File') && ~isempty(SigInfo.nd2File)
    file=char(SigInfo.nd2File);
    if ~isfile(file), file=fullfile(p,file); end
else
    [parent,name]=fileparts(p);
    file=fullfile(parent,[name '.nd2']);
end
assert(isfile(file),'MS2:MissingND2','ND2 file does not exist: %s',file);
[Channel_G,Channel_R]=nd2read(file);
assert(isequal(size(Channel_G),size(Channel_R)) && ...
    size(Channel_G,1)==SigInfo.xvalue && size(Channel_G,2)==SigInfo.yvalue && ...
    size(Channel_G,3)==SigInfo.tvalue*SigInfo.zvalue, ...
    'MS2:ND2Dimensions','ND2 dimensions do not match the configured XY, time and Z counts.');
rawfolder=fullfile(p,SigInfo.rawfolder);
if ~isfolder(rawfolder), mkdir(rawfolder); end
Cells=[];
if nargout>0, Cells=zeros(size(Channel_G),'uint16'); end
for frm=1:size(Channel_G,3)
    testnoise=medfilt2(Channel_R(:,:,frm));
    noisebw=imbinarize(testnoise,'adaptive','Sensitivity',0.4);
    noisebwe=imfill(noisebw,'holes');
    noisebwe=imdilate(noisebwe,strel('disk',1));
    noisebwe=imerode(noisebwe,strel('disk',1));
    noisebwe=imerode(noisebwe,strel('disk',1));
    outputSignal=(Channel_G(:,:,frm)-testnoise).*uint16(~noisebwe);
    if nargout>0, Cells(:,:,frm)=outputSignal; end
    z=mod(frm-1,SigInfo.zvalue)+1;
    t=floor((frm-1)/SigInfo.zvalue)+1;
    imwrite(outputSignal,fullfile(rawfolder,sprintf('DenRaw_t%04d_z%04d.tif',t,z)), ...
        'Compression','none');
end
end
