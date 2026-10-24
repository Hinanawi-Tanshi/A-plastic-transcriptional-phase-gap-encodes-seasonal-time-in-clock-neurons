function [c,ci] = fitMipSig(imagename,sigloc,sig,th)

    cookiecrop = double(imagename);
    
    %meancell = sum(cookiecrop,"all")/sum(cookiebw,"all");
%% get mean of cell
    cellnosig = zeros(size(cookiecrop));
    cellnosig(sigloc(1),sigloc(2)) = 1;
    cellnosig = ~cellnosig;
    cellnosig = cookiecrop.*double(cellnosig);
    %imagenosigvec = double(cookiecrop(:));imagenosigno0 = imagenosigvec(imagenosigvec>0);
    imagenosigvec = cellnosig(:);imagenosigno0 = imagenosigvec(imagenosigvec>0);
    pd = fitdist(imagenosigno0,'Kernel');
    ci = icdf(pd,th);
    c = [];
    if sig>=ci
        c = 1;
    elseif sig<ci
        c = 0;
    end
end