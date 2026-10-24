function bboxes = localSearchBbox(bboxSize)
% bboxSize, 局部搜索包围框的宽度，需要-1, WH-1
% bboxSize = [3 4 5];
bboxes = [];
for WH = bboxSize
    if WH==3
        [xx,yy] = meshgrid(-2:0,-2:0);
    elseif WH==4
        [xx,yy] = meshgrid(-2:-1,-2:-1);
    elseif WH==5
        [xx,yy] = meshgrid(-2:-2,-2:-2);
    end
    xx = reshape(xx,[],1);
    yy = reshape(yy,[],1);
    bbox1 = [xx yy ones(length(xx),1)*(WH-1) ones(length(xx),1)*(WH-1)];
    bboxes = [bboxes; bbox1];
end
end % end of func