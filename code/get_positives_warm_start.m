function [testeable_positives] = get_positives_warm_start(R)
    %R interaction matrix with a valid positive sampling for testing
    testeable_positives = [];
    %positive_pool to check which 
    positive_pool = R;
    
    %get the positives
    [pos_row, pos_col] = find(R);
    positives = [pos_row pos_col];
    %while there are still positives to sample
    while (size([positives]) > 0)
        %get a random positive
        index_pos = randi(size(positives,1));
        row = positives(index_pos,1);
        col = positives(index_pos,2);
        %to check if we can remove the positive from the training set, we
        %remove the positive from a temporary matrix equal to the current state
        tmp_mat = positive_pool;
        tmp_mat(row,col) = 0;
        %if the removal of the positive does not causes a cold start, then
        %we remove it from the training matrix
        if (sum(sum(tmp_mat,1) == 0) == 0) && (sum(sum(tmp_mat,2) == 0) == 0)
            testeable_positives = [testeable_positives;[row col]];
            positive_pool(row,col) = 0;
        end
        %we remove the positive from the list of positive samples
        positives(index_pos,:) = [];
    end
testeable_positives = [testeable_positives ones(size(testeable_positives,1),1)];
end