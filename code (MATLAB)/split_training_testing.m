function [R_train,R_test] = split_training_testing(R,nfolds,train_ratio)
    %function to make a train/test split such that the ratio of negatives
    %to positives in testing is equal to the ratio of positives to
    %negatives in the datset, and the ratio of positives to negatives in
    %
    %the output consist on two cell arrays containing the training and testing for each fold
    %each row in R_train and R_test is a sample of R, where the first element is
    %the row index, the second is the column index and the third column is
    %the value of the sample in R
    % if train_ratio is equal to 
    warning('off','all');

    if train_ratio == 0
        train_ratio = floor((numel(R)- sum(R(:)))/sum(R(:)));
    end

    R_test = get_test_warm_start_true_ratio(R,nfolds);

    for i = 1:nfolds

        %local store of test set
        test_fold = R_test{i};

        %get the test positives and remove them from training pool
        test_fold_positives = test_fold(test_fold(:,3)>0,:);
        tmp_mat = R;
        tmp_mat(sub2ind(size(R),test_fold_positives(:,1),test_fold_positives(:,2))) = 0;

        %get positives in training
        [row_positives,col_positives] = find(tmp_mat);
        train_fold = [row_positives col_positives ones(size(row_positives))];

        %get test negatives and remove them from the training pool
        test_fold_negatives = test_fold(test_fold(:,3)==0,:);
        tmp_mat(sub2ind(size(R),test_fold_negatives(:,1),test_fold_negatives(:,2))) = 1;
        [row_negatives,col_negatives] = find(~tmp_mat);
        
        n_positives_in_training = size(train_fold,1);
        n_negatives_in_training = n_positives_in_training * train_ratio;
        neg_index_training = randi(size(row_negatives,1),[n_negatives_in_training 1]);
        train_fold = [train_fold;row_negatives(neg_index_training) col_negatives(neg_index_training) zeros(size(neg_index_training))];
        R_train{i} = train_fold;
    end

    warning('on','all');


end