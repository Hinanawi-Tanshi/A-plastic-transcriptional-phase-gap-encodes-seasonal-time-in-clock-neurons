function GenFig4MipPrvRd_SignalCre (sigbtn,fig)
    %[x,y] = ginput(1);
    %% Update Zt
    t = sigbtn.UserData;
    signalflag = fig.UserData.MIPsig(t,end);
    if signalflag==1
         newflag = 0;
    elseif signalflag == 0
        newflag = 1;
    end
    sigbtn.Text = num2str(newflag);
    
    ax2change = findobj(fig, 'Tag','SignalImage1','UserData',{fig.UserData.MIPsig(t,:), t, 1});
    s2change = findobj(ax2change, 'Type', 'Scatter');
    
    if newflag == 1
        s2change.MarkerEdgeColor = [1,0,0];
    elseif newflag == 0
        s2change.MarkerEdgeColor = [1,1,0];
    end

    fig.UserData.MIPsig(t,end) = newflag;
    ax2change.UserData{1} = fig.UserData.MIPsig(t,:);
    %% Update cellsignal pannel
    ax_spot = findobj(fig, 'Tag', 'SpotAxe');
    ax_spot.Children.YData = fig.UserData.MIPsig(:,end);
end