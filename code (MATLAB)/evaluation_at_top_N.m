function [precision, recall] = evaluation_at_top_N(labels, scores, top_N)
    %given a label vector and a prediction scores vector, compute the
    %recall at top N
    [~, ranking] = sort(scores, 'descend');
    
    true_positives = find(labels);
    
    hits = zeros(1,length(top_N));
    for i = 1:sum(labels)
        rank_of_tp = find(ranking == true_positives(i));
        %if rank_of_tp <= top_N
            hits = hits + (rank_of_tp <= top_N);
        %end
    end
    
    precision = hits./top_N;
    recall = hits./sum(labels);

end