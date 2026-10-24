function fig = GenFigure_MIPsig_ci(CellSignal,trange)
% Edit measurements before normalization; the GUI uses the raw 7-column layout.
assert(size(CellSignal.MIPsig,2)==7 && ~isfield(CellSignal,'NormalizedMIPsig'), ...
    'MS2:NormalizedInput','Proofreading must precede normalization.');
if nargin<2 || isempty(trange), trange=1:CellSignal.t; end
validateattributes(trange,{'numeric'},{'row','integer','positive','increasing','<=',CellSignal.t});
assert(all(diff(trange)==1),'MS2:FrameRange','Use a continuous frame range.');
%% input parameters
SigInfo = CellSignal.SigInfo;
tvalue = SigInfo.tvalue;

MipImage = CellSignal.MIP;
tmax = min([length(MipImage),trange(end)]);
trange = (trange(1):1:tmax);
MipSig = CellSignal.MIPsig;
Mipci = CellSignal.MIPci;
 % cell index

cellindex = CellSignal.CellIndex;
% diffvalue = CellSignal.diff;
% ratiovalue = CellSignal.ratio;
ci = CellSignal.ci;
%
%% Create figure
fig = uifigure('Name',['Cell ' num2str(cellindex, '%02d') ' MIPSignalDetection'],'Scrollable','on');
%;
fig.Position = [20 100 1400 700];

% Title label
% indivisual t mip images area
imgnum = size(MipImage,2);
Zt = uipanel('Parent',fig, 'Position',[0 0 400 700],'Scrollable','on','Tag','Zt');
% buttons
% save button
btn = uibutton(fig,'push','Position',[440 25 50 25], 'Text','Save',...
    'ButtonPushedFcn', @(btn,event) GenFig4MipPrvRd_SaveResults(btn,fig),'Tag','btn');
% correction button
Corbtn = uibutton(fig,'push','Position',[500 25 100 25], 'Text','Correction',...
    'ButtonPushedFcn', @(Corbtn,event) GenFig4MipPrvRd_CorCorrection(Corbtn,fig),'Tag','Corbtn');
% threshold buttons
lbl_cithresh = uilabel(fig,'Position',[650 50 100 50],'Text','ci Thresh');
txa_cithresh = uitextarea(fig,'Value', {num2str(ci)},'Position',[650 20 50 25],'Tag','CiThresh');

% lbl_diffthresh = uilabel(fig,'Position',[650 50 100 50],'Text','Diff Thresh');
% txa_diffthresh = uitextarea(fig,'Value', {num2str(diffvalue)},'Position',[650 20 50 25],'Tag','DiffThresh');
% 
% lbl_ratiothresh = uilabel(fig,'Position',[750 50 100 50],'Text','Ratio Thresh');
% txa_ratiothresh = uitextarea(fig,'Value', {num2str(ratiovalue)},'Position',[750 20 50 25],'Tag','RatioThresh');

% lbl_distthresh = uilabel(fig,'Position',[850 50 100 50],'Text','Dist Thresh');
% txa_distthresh = uitextarea(fig,'Value', {num2str(distance4spot)},'Position',[850 20 50 25],'Tag','DistThresh');
% update button

%diffvalue = str2double(txa_diffthresh.Value{1});
%ratiovalue = str2double(txa_ratiothresh.Value{1});

Upbtn = uibutton(fig,'push','Position',[950 25 50 25], 'Text','Update',...
    'ButtonPushedFcn', @(Upbtn,event) GenFig4MipPrvRd_UpdateResults_ci(Upbtn,fig,trange),'Tag','Upbtn');

% Cell signal images
Cellsigpannel = uipanel('Parent',fig, 'Position',[420 150 900 500],'Tag','CellSignalPannel');
%%
imageheight = 80;
% show images and get spot location and signal intensity
for t = trange
    MipImg = MipImage{t};
    ax0 = uiaxes('Parent',Zt,'Position',[0 imageheight*(t-trange(1)) 100 0.8*imageheight],...
            'Tag','SignalImage0'); 
    imshow (MipImg,getThresh(MipImg),'Parent',ax0);
    ax0.UserData{1} = [];
    ax0.UserData{2} = t;
    ax0.UserData{3} = 0;
    ax0.YAxis.Label.String = ['t=' num2str(t)];
    ax0.YAxis.Label.Position(1) = 0;
    
    ax1 = uiaxes('Parent',Zt,'Position',[120 imageheight*(t-trange(1)) 100 0.8*imageheight],...
            'Tag','SignalImage1'); 
    imshow (MipImg,getThresh(MipImg),'Parent',ax1);
    hold on;
    if MipSig(t,end) == 1
        p = scatter(MipSig(t,end-1),MipSig(t,end-2),'Marker','square','MarkerEdgeColor','red','LineWidth',1);
    elseif MipSig(t,end) == 0
        p = scatter(MipSig(t,end-1),MipSig(t,end-2),'Marker','square','MarkerEdgeColor','yellow','LineWidth',1);
    end
    p.Parent = ax1;
    hold on;
    ax1.UserData{1} = MipSig(t,:);
    ax1.UserData{2} = t;
    ax1.UserData{3} = 1;

% sigspot switch button
    sigbtn = uibutton('Parent',Zt,...
        'Position',[300 (20+imageheight*(t-trange(1))) 20 20], 'Text',num2str(MipSig(t,end)),...
        'ButtonPushedFcn', @(sigbtn,event) GenFig4MipPrvRd_SignalCre(sigbtn,fig),...
        'Tag','SignalButton','UserData',t);
end
fig.HandleVisibility = 'on';
%fig.UserData = [MipSig, MipSpot];
fig.UserData = CellSignal;
%% plot signal traces
figure
ax_diff = subplot(3,1,1); 
plot(MipSig(:,1));
axis tight;
ylabel({'signal';'-bg(A.U.)'})
ax_diff.Parent = Cellsigpannel;
ax_diff.Tag = 'DiffAxe';

ax_cisig = subplot(3,1,2); 
plot(MipSig(:,2)); hold on; plot(Mipci);hold off;
axis tight;
ylabel({'signal vs. threshold'})
ax_cisig.Parent = Cellsigpannel;
ax_cisig.Tag = 'CiSigAxe';

ax_spot = subplot(3,1,3); 
plot(MipSig(:,end)); 
axis tight;
ylabel({'signal';'spot'})
ax_spot.Parent = Cellsigpannel;
ax_spot.Tag = 'SpotAxe';

%%
% save([cellsignalfolder '\cell_' num2str(i,'%02d') '_tracked_signal.mat'], 'CellSignal');
% savefig(fig,[cellsignalfolder '\cell_' num2str(i,'%02d') '_ProofReading.fig']);
% close all;
end
