function [f,grad] = our_method_gradients(params,R,Rr,Rc,Sc,Sr,lambda1,lambda2)
% returns the value of the following function and its gradients
%
%        sum(alfa(i) * Sc(i))                     sum(beta(i) * Sr(i))
% W = ----------------------------- ;  H = --------------------------------
%            sum(alfa)                               sum(beta)
%
% minimising:
%
%(1/2)*norm2(R - gamma * Rc * ((Rc' * Rc + sum(alfa(i) * Sc(i))/( 1 + sum(alfa)) -
%          (1 - gamma) * ((Rr * Rr' + sum(beta(i) * Sr(i))/( 1 + sum(beta)) * Rr))^2
%           + lambda1 * 1/2 * norm2(beta)^2 + lambda2 * 1/2 * norm2(alpha)^2  
%

    %get number of SI for rows and columns
    nLeft = length(Sr);
    nRight = length(Sc);
    
    %params contains gamma, betas and alphas in that order. 
    %we extract them for clarity
    gamma = params(1);
        
    bethas = params(2:1+nLeft);
    alphas = params(2+nLeft:length(params));
    
    %sum for the denominator
    sumBetas = sum(bethas);
    sumAlphas = sum(alphas);

    %number of rows and columns
    [n_rows,n_cols] = size(R);
    
    %compute current H     
    if nLeft > 0
        leftSideSideInfo = zeros(n_rows,n_rows);
        for j = 1:nLeft
            leftSideSideInfo = leftSideSideInfo + bethas(j)*Sr{j};
            
        end
        
        H = (leftSideSideInfo)./(sumBetas);
    else
        H = zeros(n_rows,n_rows);
    end

    %compute W
    if nRight > 0
        rightSideSideInfo =  zeros(n_cols,n_cols);

        for j = 1:nRight
            rightSideSideInfo = rightSideSideInfo + alphas(j)*Sc{j};
        end
        W = (rightSideSideInfo)./(sumAlphas);
    else
        W =  zeros(n_cols,n_cols);

    end

    %pre compute some terms for faster computation of the gradients

    Rc_times_W = Rc * W;
    H_times_Rr =  H * Rr;

    RcW_HRr = Rc_times_W - H_times_Rr;
    combined_minimise = R - (gamma * RcW_HRr + H_times_Rr);

    lambda_alpha = [lambda1 lambda2*ones(1,length(alphas)-1)];
    lambda_betha = [lambda1 lambda2*ones(1,length(bethas)-1)];
    
    %get current value of the funciton
    f =  (1/2) * norm(combined_minimise,'fro') ^ 2 + (1/2) * lambda1 * norm(bethas(1),'fro')^2 + (1/2) * lambda1 * norm(alphas(1),'fro')^2 ...
        + (1/2) * lambda2 * norm(bethas(2:end),'fro')^2 + (1/2) * lambda2 * norm(alphas(2:end),'fro')^2;
    
    %compute gamma gradient
    grad_gamma = -sum(RcW_HRr(:) .*combined_minimise(:));

    grad_bethas = zeros(1,nLeft);

    %computation of betas

    %precomptation for faster calculation
    Rr_times_combined = Rr * combined_minimise';
    bethas_second_term = sum(leftSideSideInfo(:) .* Rr_times_combined(:))/(sumBetas * sumBetas);

    for j = 1:nLeft

        grad_bethas(j) = -(1 - gamma)* (sum(Sr{j}(:) .* Rr_times_combined(:))/sumBetas -bethas_second_term) + lambda_betha(j) * bethas(j) ;

    end


    %computation of betas

    %precomptation for faster calculation
    grad_alphas = zeros(1,nRight);
    
    %reorder the trace 
    %trace(Rc * rightSideSideInfo * combined_minimise') = trace(rightSideSideInfo * combined_minimise'* Rc)
    
    combined_Rc = combined_minimise'* Rc;
    alphas_second_term = sum( rightSideSideInfo(:) .* combined_Rc(:))/((sumAlphas * sumAlphas));

    for j = 1:nRight
        grad_alphas(j) = - gamma * ( sum(Sc{j}(:) .* combined_Rc(:))/sumAlphas - alphas_second_term) + lambda_alpha(j) * alphas(j) ;
        
    end

    %format output for fmincon
    if nLeft < 1 && nRight <1
        grad = [grad_gamma]';
    elseif nLeft > 0 && nRight <1
         grad = [grad_gamma grad_bethas]';
    elseif nLeft < 1 && nRight >0
         grad = [grad_gamma grad_alphas]';
    else
        grad = [grad_gamma grad_bethas grad_alphas]';
    end

end