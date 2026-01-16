function [average_recall,average_precision] = average_evaluation_at_top_N(test_set,prediction_mat,labels_mat,top_N)

    test_rows = unique(test_set(:,1));
    
    n_tests = length(test_rows);
    recall = zeros(n_tests,length(top_N));
    precision = zeros(n_tests,length(top_N));
    labels_mat = -labels_mat;
    for i = 1:n_tests
        row_to_test = test_rows(i);
        scores = prediction_mat(row_to_test,:);
        labels = labels_mat(row_to_test,:);
        positives_in_testing = test_set(test_set(:,1)==row_to_test,2);
        labels(positives_in_testing) = 1;
        col_to_remove = find(labels==-1);
        labels(col_to_remove) = [];
        scores(col_to_remove) = [];
        [precision(i,:),recall(i,:)] = evaluation_at_top_N(labels, scores, top_N);
    end
    
    average_precision = mean(precision);
    average_recall = mean(recall);

end