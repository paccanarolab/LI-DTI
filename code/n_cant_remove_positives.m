function n_no_remove = n_cant_remove_positives(R)
    cant_remove = [];
    [pos_row, pos_col] = find(R);
    positives = [pos_row pos_col];
    while (size([positives]) > 0)
        index_pos = randi(size(positives,1));
        tmp_mat = R;
        row = positives(index_pos,1);
        col = positives(index_pos,2);
        tmp_mat(row,col) = 0;
        if (sum(sum(tmp_mat,1) == 0) > 0) || (sum(sum(tmp_mat,2) == 0) > 0)
            cant_remove = [cant_remove;[row col]];
        end
        positives(index_pos,:) = [];
        n_no_remove = size(cant_remove,1);
    end
end