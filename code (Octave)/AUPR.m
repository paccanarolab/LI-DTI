% function [AUPR_test] = AUPR(labels, scores)
% 
%     % test set evaluation
%     %  YTest is an 3 x N matrix with columns: row   column   label                               
%     %  where in the original matrix: R(row,column) = label,              
% 
%     [~,~,~, AUPR_test] = perfcurve(labels,scores, 1,'xCrit', 'reca', 'yCrit', 'prec');
% 
% end

function [AUPR_test] = AUPR(labels, scores)

    labels = labels(:) > 0;
    scores = scores(:);

    % Sort scores in descending order
    [scores_sorted, order] = sort(scores, 'descend');
    labels_sorted = labels(order);

    P = sum(labels_sorted == 1);

    if P == 0
        AUPR_test = NaN;
        return;
    end

    % Group by unique score thresholds, similar to perfcurve behavior
    [~, last_idx] = unique(scores_sorted, 'last');

    TP_all = cumsum(labels_sorted == 1);
    FP_all = cumsum(labels_sorted == 0);

    TP = TP_all(last_idx);
    FP = FP_all(last_idx);

    recall = TP / P;
    precision = TP ./ (TP + FP);

    X = [0; recall];
    Y = [1; precision];

    [X, sort_idx] = sort(X);
    Y = Y(sort_idx);

    AUPR_test = trapz(X, Y);

end