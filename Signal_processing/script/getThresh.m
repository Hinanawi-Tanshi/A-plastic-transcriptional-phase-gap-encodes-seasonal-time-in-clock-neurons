function Th = getThresh(image)
    A = image(:);
    B = A(A>0);
    if isempty(B)
        Th = [0 1];
    else
        Th = double([min(B) max(B)]);
        if Th(2)<=Th(1), Th(2)=Th(1)+1; end
    end
end
