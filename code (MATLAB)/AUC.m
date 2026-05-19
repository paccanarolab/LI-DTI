function [AUC_test] = AUC(labels, scores)
      
    % test set evaluation
    %  YTest is an 3 x N matrix with columns: row   column   label                               
    %  where in the original matrix: R(row,column) = label,              

    
    [~,~,~, AUC_test] = perfcurve(labels,scores,1);

%     [~,~,~, AUPR_test] = perfcurve(labels,scores, 1,'xCrit', 'reca', 'yCrit', 'prec');


end