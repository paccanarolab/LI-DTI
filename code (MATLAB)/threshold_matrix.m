function [filtered_sim] = threshold_matrix(sim,threshold)
    filtered_sim = sim;
    filtered_sim(filtered_sim < threshold) = 0;
    
    [max_sim_Val, max_sim_Idx] = max(sim,[],2);
    
    ind = sub2ind(size(sim), [1:size(sim,1)], max_sim_Idx');
    
    filtered_sim(ind) = max_sim_Val;
    

end