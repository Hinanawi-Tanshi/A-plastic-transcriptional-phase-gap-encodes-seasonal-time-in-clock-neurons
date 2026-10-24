function im1 = imsizeNorm(im0,sizeIm,fixpoint)
% 将 im0 的尺寸调整为 sizeIm
% crop/zero-padding
% fixpoint == 'UpprLeft' 零点
switch fixpoint
    case 'UpprLeft'
        im1 = im0;

        % row
        sizeOrigi = size(im1);
        if length(sizeOrigi)==2; sizeOrigi(3)=1; end
        if sizeOrigi(1)>sizeIm(1)
            im1(sizeIm(1)+1:end,:,:)=[];
        elseif sizeOrigi(1)<sizeIm(1)
            TR1 = zeros(sizeIm(1)-sizeOrigi(1),sizeOrigi(2),sizeOrigi(3),'single');
            im1 = cat(1,im1,TR1);
        end
        % col
        sizeOrigi = size(im1);
        if length(sizeOrigi)==2; sizeOrigi(3)=1; end
        if sizeOrigi(2)>sizeIm(2)
            im1(:,sizeIm(2)+1:end,:)=[];
        elseif sizeOrigi(2)<sizeIm(2)
            TR1 = zeros(sizeOrigi(1),sizeIm(2)-sizeOrigi(2),sizeOrigi(3),'single');
            im1 = cat(2,im1,TR1);
        end
        % z-stack
        sizeOrigi = size(im1);
        if length(sizeOrigi)==2; sizeOrigi(3)=1; end
        if sizeOrigi(3)>sizeIm(3)
            im1(:,:,sizeIm(3)+1:end)=[];
        elseif sizeOrigi(3)<sizeIm(3)
            TR1 = zeros(sizeOrigi(1),sizeOrigi(2),sizeIm(3)-sizeOrigi(3),'single');
            im1 = cat(3,im1,TR1);
        end
    case 'Center'
        % no func
end
end % end of func