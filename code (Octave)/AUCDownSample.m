% function [AUC_test, AUPR_test, X1,Y1,X2,Y2  ] = AUCDownSample(YTest, Res)
% 
%     % test set evaluation
%     %  YTest is an 3 x N matrix with columns: row   column   label                               
%     %  where in the original matrix: R(row,column) = label,              
% 
%     % Define the labels
%     labels = YTest(:,3)>0;
% 
%     % Get the predicted scores.
%     idx = sub2ind(size(Res),YTest(:,1),YTest(:,2));
%     scores = Res(idx);
% 
%     [X1,Y1,~, AUC_test] = perfcurve(labels,scores,1);
% 
%     [X2,Y2,~, AUPR_test] = perfcurve(labels,scores, 1,'xCrit', 'reca', 'yCrit', 'prec');
% 
% 
% end

% function [AUC_test, AUPR_test, X1, Y1, X2, Y2] = AUCDownSample(YTest, Res)
% 
%     labels = YTest(:,3) > 0;
%     idx = sub2ind(size(Res), YTest(:,1), YTest(:,2));
%     scores = Res(idx);
% 
%     % Sort scores descending
%     [scores_sorted, order] = sort(scores, 'descend');
%     labels_sorted = labels(order);
% 
%     P = sum(labels_sorted == 1);
%     N = sum(labels_sorted == 0);
% 
%     % Group by unique score thresholds
%     [~, last_idx] = unique(scores_sorted, 'last');
% 
%     TP_all = cumsum(labels_sorted == 1);
%     FP_all = cumsum(labels_sorted == 0);
% 
%     TP = TP_all(last_idx);
%     FP = FP_all(last_idx);
% 
%     % ROC curve
%     TPR = TP / P;
%     FPR = FP / N;
% 
%     X1 = [0; FPR; 1];
%     Y1 = [0; TPR; 1];
% 
%     [X1, sort_idx] = sort(X1);
%     Y1 = Y1(sort_idx);
% 
%     AUC_test = trapz(X1, Y1);
% 
%     % Precision-recall curve
%     recall = TP / P;
%     precision = TP ./ (TP + FP);
% 
%     X2 = [0; recall];
%     Y2 = [1; precision];
% 
%     [X2, sort_idx] = sort(X2);
%     Y2 = Y2(sort_idx);
% 
%     AUPR_test = trapz(X2, Y2);
% 
% end

function [AUC_test, AUPR_test, X1, Y1, X2, Y2] = AUCDownSample(YTest, Res)

    labels = YTest(:,3) > 0;

    idx = sub2ind(size(Res), YTest(:,1), YTest(:,2));
    scores = Res(idx);

    % Sort scores descending
    [scores_sorted, order] = sort(scores, 'descend');
    labels_sorted = labels(order);

    P = sum(labels_sorted == 1);
    N = sum(labels_sorted == 0);

    TP = cumsum(labels_sorted == 1);
    FP = cumsum(labels_sorted == 0);

    % ROC
    TPR = TP / P;
    FPR = FP / N;

    X1 = [0; FPR; 1];
    Y1 = [0; TPR; 1];

    AUC_test = trapz(X1, Y1);

    % Precision-Recall
    recall = TP / P;
    precision = TP ./ (TP + FP);

    X2 = [0; recall];
    Y2 = [1; precision];

    AUPR_test = trapz(X2, Y2);

end