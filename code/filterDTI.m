function [filteredMatrix] = filterDTI(Rtest, drugSim,targetSim,drugTresh,targetTresh,filterSource)
%
%  
    drugIndex = Rtest(:,1);
    targetIndex = Rtest(:,2);
    if strcmpi(filterSource,'drug') || strcmpi(filterSource,'both')
        filterLength = length(drugSim);
        for i = 1:filterLength
            tmp_Sim = drugSim{i} > drugTresh{i};
            tmp_Sim = tmp_Sim - eye(size(tmp_Sim));
            for j = 1 : size(drugIndex)
                if sum(tmp_Sim(drugIndex(j),:)) > 0
                    Rtest(j,1) = 0;
                    Rtest(j,2) = 0;
                    Rtest(j,3) = 0;
                end
            end
        end  
    end
    
    if strcmpi(filterSource,'target') || strcmpi(filterSource,'both')
        filterLength = length(targetSim);
        for i = 1:filterLength
            tmp_Sim = targetSim{i} > targetTresh{i};
            tmp_Sim = tmp_Sim - eye(size(tmp_Sim));
            for j = 1 : size(targetIndex)
                if sum(tmp_Sim(targetIndex(j),:)) > 0
                    Rtest(j,1) = 0;
                    Rtest(j,2) = 0;
                    Rtest(j,3) = 0;
                end
            end
        end  
    end
    indexFilter = sum(Rtest,2) == 0;
    filteredMatrix = Rtest;
    filteredMatrix(indexFilter,:) =[];
    
    


end

