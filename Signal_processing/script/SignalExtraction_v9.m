function SignalExtraction_v9(SigInfo,i,signalsize,bgsize,th)
%This version remove calculation on single zstack
%% input parameters
tvalue = SigInfo.tvalue;
zvalue = SigInfo.zvalue;
xvalue = SigInfo.xvalue;
yvalue = SigInfo.yvalue;
 % cell index
singlesavep = SigInfo.targetfolder;
%maskfolder = SigInfo.trackfolder;
signalfolder = SigInfo.signalfolder;
spotfolder = SigInfo.spotfolder;
MIPfolder = SigInfo.MIPfolder;
p = SigInfo.p;

if ~isfolder(MIPfolder)
    mkdir(MIPfolder)
end

if ~isfolder(signalfolder)
    mkdir(signalfolder)
end
%i:cell index

CellSignal.path = SigInfo.p;
CellSignal.t = tvalue;
CellSignal.z = zvalue;
CellSignal.xyResolution = [xvalue yvalue];
CellSignal.CellIndex = i;
CellSignal.Cycle = SigInfo.Cycle;
CellSignal.T0 = SigInfo.StartTime;% in hour
CellSignal.FrameInterval = SigInfo.FrameInterval; % in seconds
CellSignal.Celltype = '?';
CellSignal.SigInfo = SigInfo;
CellSignal.ci = th;

count=0;
%% 
singleZT = [singlesavep, '\cell', num2str(i)];
MIPsig = [];
MIPci = zeros(1,tvalue);
cellsignalfolder = [p '\' signalfolder '\cell_' num2str(i, '%02d')];
if ~isfolder (cellsignalfolder)
    mkdir(cellsignalfolder)
end
mipfile = fullfile(p,MIPfolder,sprintf('cell_%02d.tiff',i));
assert(~isfile(mipfile), 'MS2:ExistingMIP', 'Cell MIP TIFF already exists: %s',mipfile);
%% generate max projection images
for t = 1:tvalue
    % read all z images for cell i of time t
    t2read = [singleZT, '_t',num2str(t,'%04d'),'*.tif'];
    im2read = dir(t2read);
    cellMIP = []; % never reuse the preceding frame when this frame is missing
    if ~isempty(im2read)
        ImagesZseries = cell(length(im2read),2);
        MIP = [];
        zstck = length(im2read);
        list = [];
        for z = 1:zstck
            imagename = [im2read(z).folder, '\', im2read(z).name];
            cookie = imread(imagename);
            [r0,c0] = find(cookie>0);
            list = [list;r0,c0];
            %[xmin ymin width height]
        end
        assert(~isempty(list), 'MS2:EmptyCellImage', ...
            'Cell %d has no positive image pixels at time %d.',i,t);
        r = (list(:,1));c = (list(:,2));
        for z = 1:zstck
            imagename = [im2read(z).folder, '\', im2read(z).name];
            cookie = imread(imagename);
            try
                cookiecrop = imcrop(cookie,[min(c),min(r),max(c)-min(c), max(r)-min(r)]);
            catch
                cookiecrop = imcrop(cookie,[0,0,1,1]);
            end
            ImagesZseries{z,1} = im2read(z).name;
            ImagesZseries{z,2} = cookiecrop;
            MIP(:,:,z) = cookiecrop;
        end
        CellSignal.RawImage{t} = ImagesZseries;
        cellMIP = max(MIP,[],3);
        CellSignal.MIP{t} = cellMIP;
    end
    assert(~isempty(cellMIP), 'MS2:MissingCellFrame', ...
        'No images for cell %d at time %d.',i,t);
    %% Get Signal on MIP
    
    [sig,row,col] = GetSignal_v4((cellMIP),signalsize,bgsize);
    spot_initial = [];
    civalue = [];
    %% judgment of signal
%     if sig(1)>diffvalue && ((sig(1)/sig(3))*100) >ratiovalue
%         spot_initial = 1;
%     else
%         spot_initial = 0;
%     end
    [spot_initial,civalue] = fitMipSig(cellMIP,[row, col],sig(2),th);
    %% save data
    MIPsig(t,:) = [sig,row,col,spot_initial];
    MIPci(t) = civalue;
    fprintf(repmat('\b',1,count));   
    count=fprintf('current t is : %d',t);

end
fprintf('\n');
CellSignal.MIPsig = MIPsig;
CellSignal.MIPci = MIPci;

%% save mip images
mipsavepath = [p '\' MIPfolder];
imgnum = size(CellSignal.MIP,2);
imagesize = [];
for img = 1:imgnum
    imagesizei = size(CellSignal.MIP{img});
    imagesize = [imagesize;imagesizei];
end
finalimagesize = max(imagesize);

for img = 1:imgnum
    Ioriginal = uint16(CellSignal.MIP{img});
    Iresize = padarray(Ioriginal,finalimagesize-size(Ioriginal),0,'post');
    imwrite(Iresize,[mipsavepath '\cell_' num2str(CellSignal.CellIndex, '%02d') '.tiff'],...
        "WriteMode","append","Compression","none");
end

%% generating time trace
CellSignal.timetrace = generatetimetrace(CellSignal);
%%
save([cellsignalfolder '\cell_' num2str(CellSignal.CellIndex,'%02d') '_tracked_signal.mat'], 'CellSignal');
% SaveMFile([cellsignalfolder '\cell_' num2str(i,'%02d') '_ProofReading.m']);
%SaveMFile_MIP([cellsignalfolder '\cell_' num2str(CellSignal.CellIndex,'%02d') '_ProofReading.m']);
%SaveMFile_MIPsig([cellsignalfolder '\cell_' num2str(CellSignal.CellIndex,'%02d') '_ProofReading_MIPsig.m']);
SaveMFile_MIPsig_ci([cellsignalfolder '\cell_' num2str(CellSignal.CellIndex,'%02d') '_ProofReading_MIPsig_ci.m']);
% savefig(fig,[cellsignalfolder '\cell_' num2str(i,'%02d') '_ProofReading.fig']);
% close all;
disp (['cell ' num2str(i) ' is done.'])
end