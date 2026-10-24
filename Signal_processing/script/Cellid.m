function F = Cellid(cellnum)
    F = uifigure('Name','Cell Identification','Scrollable','on');
    F.Position = [1100 300 300 600];
    cellnum = double(cellnum);
    for i = 1:cellnum
        cellid(i)=uilabel(F,'Position',[50 30*i 50 30],'Text',['Cell ' num2str(i)]);
        dd(i) = uidropdown(F,'Items',{'l-LNv','s-LNv','5th s-LNv','LNd',...
            'DN1p', 'DN1a','DN2','DN3'},'Editable','on');
        dd(i).Parent = F;
        dd(i).Position = [100 30*i 100 30];
    end
    
    sbtn = uibutton(F,"ButtonPushedFcn", @(sbtn,event) savecellid(sbtn,dd));
    sbtn.Text = 'Save';
    sbtn.Position = [225 30 50 30];
end

function savecellid(sbtn,dd)
    cellids = {};
    for d = 1:length(dd)
        cellids{d,1} = ['cell_' num2str(d,'%02d')];
        cellids{d,2} = dd(d).Value;
    end
    save('Cell_ID.mat', "cellids");
    disp('CellIDs saved!')
    delete(ancestor(sbtn,'figure'));
end
