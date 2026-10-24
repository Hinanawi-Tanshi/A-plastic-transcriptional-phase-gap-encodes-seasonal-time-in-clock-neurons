function OptimizeResults(SigInfo)
%%
[~,~,imFiles] = listSignalFolders(fullfile(SigInfo.p,SigInfo.signalfolder));

dir_SF = fullfile(pwd,'ProcessedBySF_traces');
if isfolder(dir_SF)
    entries = dir(dir_SF);
    assert(all(ismember({entries.name},{'.','..'})), 'MS2:ExistingOptimization', ...
        'Optimization output directory is not empty: %s',dir_SF);
else
    mkdir(dir_SF);
end

% 用hashtable记录每一类细胞最大的亮度值ymaxV
para_str = {'max-all-proj','max-3slice-proj','max-3slice-proj&max-1slice-cal','mean-all-proj','mean-3slice-proj'};
mapping = load('Cell_ID.mat','cellids');
cellids = mapping.cellids;
cellTypes = {};
for i = 1:length(cellids)
    found = 0;
    for j = 1:length(cellTypes)
        if strcmp(cellTypes{j},cellids{i,2}) == 1
            found = 1;
            break;
        end
    end
    if found ==0
        cellTypes{end+1} = cellids{i,2};
    end
end
keys = {};
for i = 1:length(cellTypes)
    for j = 1:length(para_str)
        keys{end+1} = strcat(cellTypes{i}, para_str{j});
    end
end
values = cell(1, length(keys));
for i = 1:length(keys)
    values{i} = 0;
end
hashtable = containers.Map(keys, values);

for fi = 1:length(imFiles)
    fileTar = imFiles{fi};
    load(fullfile(fileTar),'CellSignal')
    % figure; plot(CellSignal.MIPsig(:,1:4)); title('Original processed traces')
    AllMIPsigNewR = {};
    % 提取原始数据信号 dataRaw，x*y*z*t
    dataRaw = single([]);
    sizeIm = [];
    for ti = 1:length(CellSignal.RawImage)
        TR1 = CellSignal.RawImage{ti};
        TR2 = TR1(:,2);
        TR2 = permute(TR2,[2 3 1]);
        TR2 = single(cell2mat(TR2));   % x*y*z
        if isempty(sizeIm)
            sizeIm = [size(TR2,1) size(TR2,2) size(TR2,3)];
        end
        TR3 = imsizeNorm(TR2,sizeIm,'UpprLeft');
        dataRaw = cat(4,dataRaw,TR3);  % x*y*z*t
    end
    % 信号点位置对齐
    roiLocs = round(mean(CellSignal.MIPsig(:,[5 6])));   % (row,col)
    dataAlign = dataRaw;
    for ti = 1:size(dataAlign,4)
        locsTrans = roiLocs-CellSignal.MIPsig(ti,[5 6]); % (row,col)
        xyTrans = locsTrans([2 1]);     % (x,y)
        TR1= imtranslate(dataAlign(:,:,:,ti),[xyTrans 0],'OutputView','same');
        dataAlign(:,:,:,ti) = TR1;
    end
    % visualization
    dataVis = cat(2,dataRaw,dataAlign);
    dataVis_maxP = squeeze(max(dataVis,[],3));
    dataVis_maxP = dataVis_maxP/max(dataVis_maxP,[],'all');
    % implay(dataVis_maxP)

    % 所有层最大投射或max-3 slice 最大投射
    data1 = permute(dataAlign,[3 4 1 2]);
    data2 = mean(mean(mean(data1,4),3),2);
    data2 = smoothdata(data2,'movmean',3);
    [sa,sb] = sort(data2,'descend');
    % settings
    para_ind = {1:size(dataAlign,3),sb(1:3),sb(1),1:size(dataAlign,3),sb(1:3)};
    CellSignal.Methods = para_str;
    figure(1); clf; tiledlayout(2,3);
    for para_i = 1:length(para_str)
        dataAlign1 = dataAlign(:,:,para_ind{para_i},:);
        if contains(para_str{para_i},'cal')
            dataAlign_projZ = squeeze(max(dataAlign1,[],3));
        elseif contains(para_str{para_i},'max')
            dataAlign_projZ = squeeze(max(dataAlign1,[],3));
        elseif contains(para_str{para_i},'mean')
            dataAlign_projZ = squeeze(mean(dataAlign1,3));
        end

        % 扩充边界，zero-padding, 20 pixels，便于后续批量分析
        padsize = 20;
        % row
        TR_size = size(dataAlign_projZ);
        TR1 = ones(padsize,TR_size(2),TR_size(3));
        dataAlign_projZ = cat(1,TR1,dataAlign_projZ,TR1);
        % col
        TR_size = size(dataAlign_projZ);
        TR1 = ones(TR_size(1),padsize,TR_size(3));
        dataAlign_projZ = cat(2,TR1,dataAlign_projZ,TR1);
        % update after extension/zeropading
        roiLocsNew = roiLocs + [1 1]*padsize;
        roiXYNew = roiLocsNew([2 1]);  % (row,col) --> (x,y)
        bboxes = localSearchBbox([4 5]); % [3 4 5]
        bboxes(:,[1 2]) = bboxes(:,[1 2]) + repmat(roiXYNew,size(bboxes,1),1);
        bboxBack = [-7 -7 14 14] + [roiXYNew 0 0];

        MIPsigNew = zeros(size(dataAlign_projZ,3),2);
        signalAll = zeros(size(dataAlign_projZ,3),size(bboxes,1));
        for ti = 1:size(dataAlign_projZ,3)
            im0 = dataAlign_projZ(:,:,ti);
            for bi = 1:size(bboxes,1)
                im1 = imcrop(im0,bboxes(bi,:));
                im1_data = sort(reshape(im1,[],1),'descend');
                im1_data = im1_data(1:round(length(im1_data)*0.75));
                im1_data(im1_data==0) = nan;
                signalAll(ti,bi) = mean(im1_data,'all','omitnan');
            end
            % 用最大的信号bbox，计算背景
            [ma,mb] = max(signalAll(ti,:));
            bboxTar = bboxes(mb,:);
            im0_1 = im0;
            % Legacy expression retained for numerical reproducibility. The X
            % endpoint uses bboxTar(2) in the original implementation.
            im0_1(bboxTar(2):bboxTar(2)+bboxTar(4),bboxTar(1):bboxTar(1)+bboxTar(2)) = 0;
            im2 = imcrop(im0_1,bboxBack);
            im2(im2==0) = nan;
            signalBack = mean(im2,"all",'omitnan');
            % update signal and background
            MIPsigNew(ti,:) = [ma signalBack];
        end
        %% 信号补偿
        dataAlign1 = reshape(dataAlign,[],1);
        dataAlign1(dataAlign1==0)=[];
        TR_data = prctile(dataAlign1,[25 50]);
        MIPsigNew(:,1) = MIPsigNew(:,1)-(TR_data(2)-TR_data(1));
        MIPsigNew(MIPsigNew(:,1)<0,1) = 0;
        MIPsigNew = [MIPsigNew(:,1)-MIPsigNew(:,2) MIPsigNew];
        % 取基准值原则更新：使信号的[0-25th]<0
        MIPsigNewR = MIPsigNew;
        MIPsigNewR(:,1) = MIPsigNewR(:,1) - prctile(MIPsigNewR(:,1),25);
        MIPsigNewR(MIPsigNewR(:,1)<0,1) = 0;
        AllMIPsigNewR{end+1} = MIPsigNewR;
        %% 不进行信号补偿
        %         MIPsigNew = [MIPsigNew(:,1)-MIPsigNew(:,2) MIPsigNew];
        %         MIPsigNewR = MIPsigNew;
        %         AllMIPsigNewR{end+1} = MIPsigNewR;
        %%
        % visualization
        plotInd = 1;
        ymaxV = max([MIPsigNewR(:,plotInd)],[],'all');
        % 按cell类型记录ymaxV
        if ymaxV >  hashtable(strcat(cellids{CellSignal.CellIndex,2},para_str{para_i}))
            hashtable(strcat(cellids{CellSignal.CellIndex,2},para_str{para_i})) = ymaxV;
        end
        ymaxV = ceil(ymaxV/50)*50;
        ymaxV = max([ymaxV 300]);
        ymaxV = ymaxV+30;
        nexttile; title(para_str{para_i})
        plot(MIPsigNewR(:,plotInd))
        y1 = ymaxV-30;
        y2 = ymaxV;
        for ti = 1:length(CellSignal.RawImage)
            if CellSignal.MIPsig(ti, end) == 1
                %text(ti,CellSignal.MIPsig(ti,plotInd),'+');
                ton = ti-0.1;
                toff = ti+0.1;
                xloc = [ton, ton, toff, toff];
                yloc = [y1 ,y2, y2, y1];
                patch (xloc, yloc, [0 0.4470 0.7410] ,'EdgeColor', 'none');
            end
        end
        set(gcf,'Position',[100 100 700 200]);
        save(fileTar,'para_str','AllMIPsigNewR','-append');
    end
    % export
    [cellFolder,matName] = fileparts(fileTar);
    [~,cellName] = fileparts(cellFolder);
    dir1 = fullfile(dir_SF,cellName);
    if ~exist(dir1,'dir')
        mkdir(dir1)
    end
    imname = [matName '.png'];
    saveas(gcf,fullfile(dir1,imname))
end
save('Hash_ymaxV.mat','hashtable','-mat');
end
