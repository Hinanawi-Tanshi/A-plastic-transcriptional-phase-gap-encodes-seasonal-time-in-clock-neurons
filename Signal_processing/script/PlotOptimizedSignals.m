function PlotOptimizedSignals(SigInfo)
%% input parameters
signalfolder = SigInfo.signalfolder;
p = pwd;
rootpath = [p '\' signalfolder];
signalsubfolder = listSignalFolders(rootpath, SigInfo);
tvalue = SigInfo.tvalue;
ytick = load("Cell_ID.mat");

%% timing trace generation
Firstframe = datevec(SigInfo.StartTime); % ZT time for first captured frame
FrameRate = SigInfo.FrameInterval ;% in second
%==========================end of input parameter==========================
FrameInterval = ((FrameRate)/60)/60; % in hour
timevec_hzt = Firstframe(4);
timevec_m = Firstframe(5);
timevec_num = timevec_hzt + timevec_m/60;
timetrace = timevec_num: FrameInterval: (timevec_num + FrameInterval * (tvalue-1));

rw = FrameInterval/2;

AllSignal = {};
%figSig = figure('Units','normalized', 'Position',[0.1 0.1 0.4 0.8]);
cellnum = length(signalsubfolder);
for para_i = 1:5
    figSig = figure('Units','normalized', 'Position',[0.1 0.1 0.4 0.8]);
    for c = 1:cellnum
        signalfile = dir([rootpath,'\',signalsubfolder(c).name,'\*tracked_signal.mat']);
        load([signalfile.folder,'\',signalfile.name], 'CellSignal');
        load([signalfile.folder,'\',signalfile.name], 'para_str');
        load([signalfile.folder,'\',signalfile.name], 'AllMIPsigNewR');
        %disp(num2str(c));
        AllSignal{c}= AllMIPsigNewR{para_i};
        pos1 = [0.1 (1-0.05)*(c-1)/cellnum+0.05  0.4 (1-0.05)/cellnum];
        % 在不够长的序列后补充0
        MIP_sig = padarray(AllMIPsigNewR{para_i}(:,1), length(timetrace)-length(AllMIPsigNewR{para_i}(:,1)), 'post');
        axintensity = subplot('Position',pos1); plot(timetrace, MIP_sig);
        %axintensity.YAxis.Visible = 'off';
        axintensity.Box = 'off';
        if c > 1
            axintensity.XAxis.Visible = 'off';
        end
        axintensity.YTick = (0:50:200);
        axis tight
        pos2 = [0.6 (1-0.05)*(c-1)/cellnum+0.05 0.4 0.5*(1-0.05)/cellnum];
        axspot = subplot('Position',pos2);
        cid = str2double(regexp(signalsubfolder(c).name,'\d+$','match','once'));
        str1 = [(ytick.cellids{cid,1}),': ', (ytick.cellids{cid,2})];
        axspot.Title.String = strrep(str1,'_','\_');
        for fmidx = 1:tvalue
            xloc = [timetrace(fmidx)-rw, timetrace(fmidx)-rw, timetrace(fmidx)+rw, timetrace(fmidx)+rw];
            yloc = [0, 1, 1, 0];
            if CellSignal.MIPsig(fmidx,end) == 1
                patch (xloc, yloc, 'red','FaceAlpha',.3, 'EdgeColor', 'none');
                %         rectangle('Position',[timetrace(fm)-rw, y0, 2*rw, 1], 'EdgeColor','none', ...
                %             'FaceColor', [0 0.4470 0.7410]);
                %         alpha(.5)
                %hold on
            elseif CellSignal.MIPsig(fmidx,end) == 0
                patch (xloc, yloc, 'blue','FaceAlpha',.05, 'EdgeColor', 'none');
                %         rectangle('Position',[timetrace(fm)-rw, y0, 2*rw, 1], 'EdgeColor','none', ...
                %             'FaceColor', [0 0.4470 0.7410]);
                %         alpha(.5)
                %hold on
            elseif isnan(CellSignal.MIPsig(fmidx,end))
                patch (xloc, yloc, 'black','FaceAlpha',.3, 'EdgeColor', 'none');
                %hold on
            end

        end
        axspot.YAxis.Visible = 'off';
        if c > 1
            axspot.XAxis.Visible = 'off';
        end

        axis tight
    end
    savefldr = [p,'\SignalFigs'];
    if ~isfolder(savefldr)
        mkdir(savefldr);
    end
    savefig(figSig,[savefldr,'\SignalFigures_', para_str{para_i},'.fig']);
    saveas(figSig,[savefldr,'\SignalFigures_', para_str{para_i},'.png']);
    save([savefldr '\timetrace.mat'], 'timetrace');
    save([savefldr '\Allsignal_' para_str{para_i} '.mat'],'AllSignal');
end
end