function [R_test] = get_test_warm_start_true_ratio(R,nfolds)
%creates nfold test sets from a binary interaction matrix R, checking that
%all rows and columns left in the training set have at list one element
    %get the positives
    [testeable_positives] = get_positives_warm_start(R);
    
    %get the ratio of the original dataset for the negatives
    true_ratio = floor((numel(R)- sum(R(:)))/sum(R(:)));
    n_testing_zeroes_requiered = size(testeable_positives,1) * true_ratio;


    [neg_row,neg_col] = find(~R);
    neg_index_testing = randi(size(neg_row,1),[n_testing_zeroes_requiered 1]);
    test_negatives = [neg_row(neg_index_testing) neg_col(neg_index_testing) zeros(n_testing_zeroes_requiered,1)];
    
    positive_partition = cvpartition([1:size(testeable_positives,1)],"KFold",nfolds);
    negative_partition = cvpartition([1:size(test_negatives,1)],"KFold",nfolds);

    for i = 1:nfolds
        R_test_positives = testeable_positives(test(positive_partition,i),:);
        R_test_negatives = test_negatives(test(negative_partition,i),:);
        R_test{i} = [R_test_positives; R_test_negatives];
    end

end