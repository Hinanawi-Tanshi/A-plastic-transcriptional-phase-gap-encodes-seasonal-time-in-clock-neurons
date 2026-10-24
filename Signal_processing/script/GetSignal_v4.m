function [sig,row,col] = GetSignal_v4(imagename,signalsize,bgsize,...
    searchpoint)

    cookiecrop = imagename;
    cookiebw = imfill(cookiecrop>0,'holes');
    meancell = sum(cookiecrop,"all")/sum(cookiebw,"all");
    %% get mean of cell based on histfit
%     imagevec = cookiecrop(:);imageno0 = imagevec(imagevec>0);
%     h = imghistfit (imageno0,50,'kernel');
%     x = h{1};
%     y = h{2};
% 
%     [~, xmeans] = findpeaks(y,x);
%     try
%         meancell = xmeans(2);
%     catch
%         meancell = xmeans(1);
%     end
    %% assign searchlocation
    if nargin == 3
        [r,c] = find(cookiebw>0);
    else
        newarea = uint16(zeros(size(cookiebw)));
        newarea(searchpoint(1), searchpoint(2))=1;
        newarea_dilate = imdilate(newarea, strel('square',5));
        [r,c] = find(newarea_dilate>0);
    end
    %% calculate
        signalz = [];
        for pxl = 1:length(r)
            pic4search = false(size(cookiebw));
            pic4search(r(pxl),c(pxl)) = 1;
            signalregion = (imdilate(pic4search,strel('square',signalsize))) .* cookiebw;
            bgregion = imdilate(pic4search,strel('square',bgsize)).* ~signalregion .* cookiebw;
            rawsignal = uint16(cookiecrop) .* uint16(signalregion);
            rawbg = uint16(cookiecrop) .* uint16(bgregion);
            signali = sum(rawsignal,"all")/sum(signalregion,"all");
            bgi = sum(rawbg,"all")/sum(bgregion,"all");
            % calculate signals
            % signal = (signali-bgi);
            signal = signali;   % changed
            signalz(pxl,:) = [signal,signali,bgi];
        end
        [~, idx] = max(signalz(:,1));
        row = r(idx);
        col = c(idx);
        sig = [signalz(idx,:),meancell];

end