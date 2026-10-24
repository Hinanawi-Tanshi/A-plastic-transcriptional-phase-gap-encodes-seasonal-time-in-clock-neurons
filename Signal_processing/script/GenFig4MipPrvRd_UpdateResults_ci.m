function GenFig4MipPrvRd_UpdateResults_ci(Upbtn,fig,trange)
%Recalculating the signal spot 0 or 1
    disp('updating...');
    %tvalue = fig.UserData.t;
    oriSigspot = fig.UserData.MIPsig(:,end);
    orisignal = fig.UserData.MIPsig(:,2);
    oriloc = fig.UserData.MIPsig(:,5:6);

    Th_ci = str2double(findobj(fig, 'Tag','CiThresh').Value);
    cellsig = findobj(fig, 'Tag','CellSignalPannel');
    Zt = findobj(fig, 'Tag','Zt');

    cicri = oriSigspot;
    cis = zeros(size(orisignal));

    for t = trange
        cellmip = fig.UserData.MIP{t};
        [c,ci] = fitMipSig(cellmip,[oriloc(t,1), oriloc(t,2)],orisignal(t),Th_ci);
        cicri(t) = c;
        cis(t) = ci;
    end

    newspottag = cicri;
    fig.UserData.MIPsig(trange,end) = newspottag(trange);
    fig.UserData.MIPci(trange) = cis(trange);
    %% Find Changed T and Replot Cellimagepannel
    Tchange = xor (oriSigspot,newspottag);
    loc = find(Tchange == 1);
    %% Replot Cellsignalpannel
    %ax_diff = findobj(cellsig, 'Tag','DiffAxe');
    ax_Cisig = findobj(cellsig, 'Tag','CiSigAxe');
    ax_spot = findobj(cellsig, 'Tag','SpotAxe');


    cla(ax_Cisig); 
    plot(ax_Cisig,(1:1:length(orisignal)),orisignal,(1:1:length(orisignal)),fig.UserData.MIPci);
    axis(ax_Cisig, 'tight');
    ylabel(ax_Cisig,{'signal vs. threshold'})
    ax_Cisig.Tag = 'CiSigAxe';
    
    cla(ax_spot);
    plot(ax_spot,(1:1:length(newspottag))',fig.UserData.MIPsig(:,end)); 
    axis(ax_spot,'tight');
    ylabel(ax_spot,{'signal';'spot'})
    ax_spot.Tag = 'SpotAxe';
    %%  Replot Zt images
    if  ~isempty(loc)
        for ind = 1:length(loc)
            ax2change = findobj(Zt, 'Tag', 'SignalImage1',...
                'UserData',{[fig.UserData.MIPsig(loc(ind),1:6),oriSigspot(loc(ind))],loc(ind),1});
            if ~isempty(ax2change)
                ax2change.UserData{1} = fig.UserData.MIPsig(loc(ind),:);
                s2change = findobj(ax2change,'Type','Scatter');
        
                if newspottag(loc(ind)) == 1
                    s2change.MarkerEdgeColor = [1,0,0];
                elseif newspottag(loc(ind)) == 0
                    s2change.MarkerEdgeColor = [1,1,0];
                end
                
                btn2change = findobj(Zt,'Tag', 'SignalButton','UserData',loc(ind));
                btn2change.Text = num2str(newspottag(loc(ind)));
            end
        end
    end
        disp('updated！')
end 