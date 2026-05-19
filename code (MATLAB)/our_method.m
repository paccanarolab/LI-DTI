function [Res] = our_method(R, gamma, SIR, SIC, alphas, bethas)

% function that computes the combined prediction as
%               H*R_r + R_c * W
% where H and W are described in the paper. Alphas are learned in other
% function


    nLeft = length(SIR);
    nRight = length(SIC);
    
    % normalise R
    
    leftNorm = normaliseY(R,1); % by row
    rightNorm = normaliseY(R,0); % by column
    
    sumBetas = 0;
    sumAlphas = 0;
    
    %compute H     
     if nLeft > 0
        leftSideSideInfo = zeros(size(R,1),size(R,1));
        for j = 1:nLeft
            leftSideSideInfo = leftSideSideInfo + bethas(j)*SIR{j};
            sumBetas = sumBetas + bethas(j);
        end
        H = (leftSideSideInfo)./(sumBetas);
    else
        H = zeros(size(R,1),size(R,1));
    end
    %compute W
    
     
    if nRight > 0
        rightSideSideInfo =  zeros(size(R,2),size(R,2));

        for j = 1:nRight
            rightSideSideInfo = rightSideSideInfo + alphas(j)*SIC{j};
            sumAlphas = sumAlphas + alphas(j);
        end
        W = (rightSideSideInfo)./(sumAlphas);
    else
        W =  zeros(size(R,2),size(R,2));

    end
    %compute combination
    HR = H * leftNorm;
    RW = rightNorm * W;
    Res = gamma * RW + (1-gamma) * HR;
end

