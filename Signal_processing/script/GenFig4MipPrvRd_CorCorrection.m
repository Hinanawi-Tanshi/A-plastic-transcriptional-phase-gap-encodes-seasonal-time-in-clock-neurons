function GenFig4MipPrvRd_CorCorrection (btn,fig)
    %% update Zt
    %[x,y] = ginput(1);
    CurAxe = fig.CurrentAxes;
    ImageData = CurAxe.Children(end).CData;
    sigsave = CurAxe.UserData{1};
    t = CurAxe.UserData{2};
    f = figure('Position',[700 350 100 100]);   
    
    imshow(ImageData,getThresh(ImageData),'InitialMagnification','fit');
    title (['t=' num2str(t)]);
    [x, y, button] = ginput(1);
    close (f);
    
    sigsz = fig.UserData.SigInfo.signalsize;
    bgsz = fig.UserData.SigInfo.bgsize;

    [signalsa,ra,ca] = GetSignal_v4(uint16(ImageData),sigsz,bgsz,[round(y) round(x)]);  % changed
    sigsave(1:end-1) = [signalsa,ra,ca];

    if button == 1
        spotflag = 1;
    elseif button == 3
        spotflag = 0;
    end
    sigsave(end) = spotflag;
    CurAxe.UserData{1} = sigsave;
    fig.UserData.MIPsig(t,:) = sigsave;

    ParentPannel = CurAxe.Parent;
    h = findobj(ParentPannel, 'Tag', 'SignalButton','Userdata',t);
    h.Text = num2str (spotflag);
    
    s = findobj(CurAxe, 'Type','Scatter');
    s.XData = round(x);s.YData = round(y);
    if spotflag == 1
        s.MarkerEdgeColor = [1 0 0];
    elseif spotflag == 0
        s.MarkerEdgeColor = [1 1 0];
    end

%% Update cellsignal pannel

    cellsig = findobj(fig, 'Tag', 'CellSignalPannel');
    % Update cell signal
    ax_diff = findobj(cellsig, 'Tag','DiffAxe');
    ax_cisig = findobj(cellsig, 'Tag','CiSigAxe');
    ax_spot = findobj(cellsig, 'Tag','SpotAxe');

    ax_diff.Children(end).YData = fig.UserData.MIPsig(:,1);
    ax_cisig.Children(end).YData = fig.UserData.MIPsig(:,2);
    ax_spot.Children.YData = fig.UserData.MIPsig(:,end);

end