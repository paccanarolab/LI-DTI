function [R_out] = top_k_combined_similarity(R,sim,k)
%function that introduces a prior interaction for empty rows based on its
%similarity to other rows
%   output the input interaction matrix R_out, where the empty rows of the
%   input matrix R have been replaced with the interactions of the top k most
%   similar rows of R, based on the similarity matrix sim
%   The interactions are weighted by the value of the similarity
%   the output matrix R_out should keep the original values of the non zero rows of R
%    
    
    R_out = R;

    %get rows with only 0
    cero_rows = sum(R,2)==0;
    indx_test = find(cero_rows);

    %put the similarity between empty from testing to 0
    %then, we build a big similarity matrix, which is a concatenation of
    %all similarity matrices
    all_sims = [];
    for i = 1:size(sim,2)
        sim{i}(cero_rows,cero_rows) = 0;
        all_sims = [all_sims sim{i}];
    end


    %get the max similarity index for rows with only zeros
    [max_sim_Val, max_sim_Idx] = maxk(all_sims(cero_rows,:),k,2);

    %if the greatest similarity is not zero, then we use the interactions from that row
    no_zero_max = sum(max_sim_Val,2) >0;
    copy_to_row = indx_test(no_zero_max);
    
    %to make things easier, we concatenate copies of the interacion matrix.
    %We then use this big matrix to initialise the interacion matrix
    bigR = repmat(R,size(sim,2),1);

    %get the top k non zero similar rows 
    copy_from_row = max_sim_Idx .* (max_sim_Val >0);
    copy_from_row = copy_from_row(any(copy_from_row,2),:);

    %and their similarities
    sim_weight = max_sim_Val(any(max_sim_Val,2),:);

    %copy to the empty rows the interactions from the most similar row,
    %wheighted by the similarity
    %
    for i = 1:size(copy_from_row,1)
        indexCopy = copy_from_row(i,1:sum(copy_from_row(i,:)>0));
        R_out(copy_to_row(i),:) = sum(bigR(indexCopy,:) .* sim_weight(i,sim_weight(i,:)>0)',1);
    end
   
end
