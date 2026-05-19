% function [AUC_test] = AUC(labels, scores)
%     % test set evaluation
%     %  YTest is an 3 x N matrix with columns: row   column   label                               
%     %  where in the original matrix: R(row,column) = label,              
%     [~,~,~, AUC_test] = perfcurve(labels,scores,1);
% 
% %     [~,~,~, AUPR_test] = perfcurve(labels,scores, 1,'xCrit', 'reca', 'yCrit', 'prec');
% 
% end

function [AUC_test] = AUC(labels, scores)

    labels = labels(:) > 0;
    scores = scores(:);

    % Sort scores in descending order
    [scores_sorted, order] = sort(scores, 'descend');
    labels_sorted = labels(order);

    P = sum(labels_sorted == 1);
    N = sum(labels_sorted == 0);

    if P == 0 || N == 0
        AUC_test = NaN;
        return;
    end

    % Group by unique score thresholds, similar to perfcurve behavior
    [~, last_idx] = unique(scores_sorted, 'last');

    TP_all = cumsum(labels_sorted == 1);
    FP_all = cumsum(labels_sorted == 0);

    TP = TP_all(last_idx);
    FP = FP_all(last_idx);

    TPR = TP / P;
    FPR = FP / N;

    X = [0; FPR; 1];
    Y = [0; TPR; 1];

    [X, sort_idx] = sort(X);
    Y = Y(sort_idx);

    AUC_test = trapz(X, Y);

end