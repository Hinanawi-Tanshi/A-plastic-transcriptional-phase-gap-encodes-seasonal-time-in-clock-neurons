function PlotSingleCellFig(SigInfo)
%% read SigInfo
signalfolder = SigInfo.signalfolder;
p = pwd;
rootpath = [p '\' signalfolder];
signalsubfolder = listSignalFolders(rootpath, SigInfo);
f = figure('Units','normalized', 'Position',[0.1 0.1 0.8 0.8]);
cellnum = length(signalsubfolder);
row = ceil(sqrt(cellnum));
col = ceil(cellnum/row);
rw = (SigInfo.FrameInterval/3600)/2;
cellid = load("Cell_ID.mat");
for cell = 1:cellnum
%% read cellsignal information
    signalfile = dir([signalsubfolder(cell).folder,'\',signalsubfolder(cell).name,'\*tracked_signal.mat']);
    sig = load([signalfile.folder,'\',signalfile.name], 'CellSignal');
    %disp(num2str(c));
    timetraceall = sig.CellSignal.timetrace;
    MipSignalall = sig.CellSignal.MIPsig;
    Signal2Plot = MipSignalall(:,1);
    %Signal2Plot = MipSignalall(:,1)./MipSignalall(:,4);
    %Signal2Plot = (MipSignalall(:,1) -min(MipSignalall(:,1))) *3;

%% plot figure on the same axes
    
    ax_intensity = subplot(row,col,cell);
    plot(timetraceall,Signal2Plot,'Parent',ax_intensity,'Color',[0 0 0]);
    hold on
    ax_intensity.Box = 'off';
    xlabel ('Hour','FontName','Arial'); 
    %ylabel('Intensity / A.U.','FontName','Arial');
    xlim ([timetraceall(1) timetraceall(end)]);
    y1 = ax_intensity.YLim(2);
    y2 = y1 + 0.2*y1;
    
    tvalue = length(timetraceall);
    for fmidx = 1:tvalue
        if MipSignalall(fmidx,end) == 1
            ton = timetraceall(fmidx)-rw;
            toff = timetraceall(fmidx)+rw;
            xloc = [ton,ton,toff,toff];
            yloc = [y1, y2, y2, y1];
            c = patch (xloc, yloc, [0 0.4470 0.7410] ,'EdgeColor', 'none');
    
    %         rectangle('Position',[timetrace(fm)-rw, y0, 2*rw, 1], 'EdgeColor','none', ...
    %             'FaceColor', [0 0.4470 0.7410]);
    %         alpha(.5)
            hold on
    %     elseif testMipSignalall(fmidx,end) == 0
    %         xloc = [timetraceall(fmidx),timetraceall(fmidx), timetraceall(fmidx+1), timetraceall(fmidx+1)];
    %         yloc = [0, 1, 1, 0];
    %         patch (xloc, yloc, 'blue','FaceAlpha',.05, 'EdgeColor', 'none');
    % %         rectangle('Position',[timetrace(fm)-rw, y0, 2*rw, 1], 'EdgeColor','none', ...
    % %             'FaceColor', [0 0.4470 0.7410]);
    % %         alpha(.5)
    %         hold on
        end
    end
    
    %ax_intensity.YLabel.String = [cellid.cellids{cell,1} ': ' cellid.cellids{cell,2}];
    %strmatch(signalsubfolder(cell).name(end-1,end), cell2(cellid.cellids{:,1}));
    cid = str2double(regexp(signalsubfolder(cell).name,'\d+$','match','once'));
    str1 = [(cellid.cellids{cid,1}),': ', (cellid.cellids{cid,2})];
    ax_intensity.YLabel.String = strrep(str1,'_','\_');
%% %plot figure on different axes
% f = figure('Position',[100 100 500 350]);
% ax_intensity = axes('Parent',f,'Units','normalized','Position',[0.1 0.1 0.8 0.6]);
% 
% plot(timetraceall,MipSignalall(:,1),'Parent',ax_intensity,'Color',[0 0 0]);
% ax_intensity.Box = 'off';
% xlabel ('CT / Hour','FontName','Arial'); 
% ylabel('Intensity / A.U.','FontName','Arial');
% xlim ([timetraceall(1) timetraceall(end)]);
% 
% ax_spot = axes('Parent',f,'Units','normalized','Position',[0.1 0.7 0.8 0.1]);
% hold on
% tvalue = length(timetraceall);
% for fmidx = 2:tvalue-1
%     if MipSignalall(fmidx,end) == 1
%         ton = (timetraceall(fmidx)+timetraceall(fmidx-1))/2;
%         toff = (timetraceall(fmidx)+timetraceall(fmidx+1))/2;
%         xloc = [ton,ton,toff,toff];
%         yloc = [0, 1, 1, 0];
%         c = patch (xloc, yloc, [0 0.4470 0.7410] ,'EdgeColor', 'none');
% 
% %         rectangle('Position',[timetrace(fm)-rw, y0, 2*rw, 1], 'EdgeColor','none', ...
% %             'FaceColor', [0 0.4470 0.7410]);
% %         alpha(.5)
%         hold on
% %     elseif testMipSignalall(fmidx,end) == 0
% %         xloc = [timetraceall(fmidx),timetraceall(fmidx), timetraceall(fmidx+1), timetraceall(fmidx+1)];
% %         yloc = [0, 1, 1, 0];
% %         patch (xloc, yloc, 'blue','FaceAlpha',.05, 'EdgeColor', 'none');
% % %         rectangle('Position',[timetrace(fm)-rw, y0, 2*rw, 1], 'EdgeColor','none', ...
% % %             'FaceColor', [0 0.4470 0.7410]);
% % %         alpha(.5)
% %         hold on
%     end
% end
% xlim ([timetraceall(1) timetraceall(end)]);
% ax_spot.XColor = 'none';
% ax_spot.YColor = 'none';
end
savefldr = [p,'\SignalFigs'];
if ~isfolder(savefldr)
    mkdir(savefldr);
end
savefig(f,[savefldr,'\IndividualSignal.fig']);
saveas(f,[savefldr,'\IndividualSignal.png']);
end