function [AUC_test, AUPR_test, X1,Y1,X2,Y2  ] = AUCDownSample(YTest, Res)
      
    % test set evaluation
    %  YTest is an 3 x N matrix with columns: row   column   label                               
    %  where in the original matrix: R(row,column) = label,              

    % Define the labels
    labels = YTest(:,3)>0;

    % Get the predicted scores.
    idx = sub2ind(size(Res),YTest(:,1),YTest(:,2));
    scores = Res(idx);
    
    [X1,Y1,~, AUC_test] = perfcurve(labels,scores,1);

    [X2,Y2,~, AUPR_test] = perfcurve(labels,scores, 1,'xCrit', 'reca', 'yCrit', 'prec');


end