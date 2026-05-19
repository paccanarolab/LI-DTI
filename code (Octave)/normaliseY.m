function [RN] = normaliseY(R,sidebyrow)
    
    [row,col] = size(R);
    if sidebyrow
        D = diag(1./sqrt(sum(R.*R,2)));
        RN = D*R;
        RN(isnan(RN)) = 1/sqrt(col);
    else
        D = diag(1./sqrt(sum(R.*R,1)));
        RN = R*D;
        RN(isnan(RN)) = 1/sqrt(row);
       
    end

     RN(isnan(RN)) = 0;
    
end
