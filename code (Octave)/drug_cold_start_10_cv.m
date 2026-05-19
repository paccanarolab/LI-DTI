%% *LI-DTI: Linear Interpretable Drug-Target Interaction*
% *Drug Cold Start Scenario*
%% 
% * Output: our_method_predictions_10CV_2026.mat file with the predictions of 
% LI-DTI for the DTINet dataset saved in ../repository/cold_start_evaluation_1/drug_cold_start  

%clear all; close all; clc;

if exist('OCTAVE_VERSION', 'builtin')
    pkg load optim
    pkg load statistics
end


%% Read the data
%% 
% Load DTI matrix

R = load('../data/luo_dataset/mat_drug_protein.txt');
%% 
% Similarities for initialisation

chemSim = load('../data/luo_dataset/Similarity_Matrix_Drugs.txt');
protSeqSim = load('../data/luo_dataset/Similarity_Matrix_Proteins.txt');
drugDrug = load('../data/luo_dataset/mat_drug_drug.txt');
sideEffect = load('../data/luo_dataset/mat_drug_se.txt');
drugDisease = load('../data/luo_dataset/mat_drug_disease.txt');
protDisease = load('../data/luo_dataset/mat_protein_disease.txt');
PPI = load('../data/luo_dataset/mat_protein_protein.txt');
%% 
% Keep only Drugs and targets with known DTI

indexProt = sum(R,1) > 0;
indexDrug = sum(R,2) > 0;
R = R(indexDrug,indexProt);
chemSim = chemSim(indexDrug,indexDrug);
chemSim = chemSim - eye(size(chemSim));
protSeqSim = protSeqSim(indexProt,indexProt);
protSeqSimN = protSeqSim./max(protSeqSim);
protSeqSimN = protSeqSimN - eye(size(protSeqSimN));
drugDrug = drugDrug(indexDrug,indexDrug);
sideEffect = sideEffect(indexDrug,:);
sideEffect = sideEffect(:,sum(sideEffect,1) > 0);
drugDisease = drugDisease(indexDrug,:);
protDisease = protDisease(indexProt,:);
PPI = PPI(indexProt,indexProt);
%%
prior_SIM_row = {chemSim};
top_k = 5;
%% 
% Compute the cosine similarity for Side Effects, drug-disease associations,DDI, 
% target-Disease associations and PPI

sideEffectSim = 1 - squareform(pdist(sideEffect,'cosine'));
drugDiseaseSim = 1 - squareform(pdist(drugDisease,'cosine'));
protDiseaseSim = 1 - squareform(pdist(protDisease,'cosine'));
PPISim = 1 - squareform(pdist(PPI,'cosine'));
PPISim(isnan(PPISim)) = 0;
drugDrugSim = 1 - squareform(pdist(drugDrug,'cosine'));
drugDrugSim(isnan(drugDrugSim)) = 0;
%% 
% Remove the Main Diagonal

sideEffectSim = sideEffectSim - eye(size(sideEffectSim));
drugDiseaseSim = drugDiseaseSim - eye(size(drugDiseaseSim));
protDiseaseSim = protDiseaseSim - eye(size(protDiseaseSim));
PPISim = PPISim - eye(size(PPISim));
drugDrugSim = drugDrugSim - eye(size(drugDrugSim));
%% 
% Set setting

nfolds = 10; % number of folds in the cv
bound = 0; % lower limit of the weights of the model
% parameters for matlabs fmincon optimiser
% options = optimoptions(@fmincon,'Algorithm','trust-region-reflective',...
%     'CheckGradients',false,'SpecifyObjectiveGradient',true,...
%     'MaxFunctionEvaluations',10000,'MaxIterations',10000,'Display','off');

if exist('OCTAVE_VERSION', 'builtin')
    options = optimset( ...
        'Algorithm', 'sqp', ...
        'GradObj', 'on', ...
        'MaxIter', 10000, ...
        'MaxFunEvals', 10000, ...
        'Display', 'off' ...
    );
else
    options = optimoptions(@fmincon, ...
        'Algorithm','trust-region-reflective', ...
        'CheckGradients',false, ...
        'SpecifyObjectiveGradient',true, ...
        'MaxFunctionEvaluations',10000, ...
        'MaxIterations',10000, ...
        'Display','off');
end

lambda1 = 1; % lambda for the drug and target covariance similarity
lambda2 = 0.1; %lambda for the remaining similarites
%% 
% Load the 10 CV splits

for k = 1:10
    for j = 1:10
        load(['../data/luo_dataset/drug_cold_start/our_method/cv_' int2str(k) '/test' int2str(j) '.mat']);
        tests{k}{j} = test_fold;
    end
end
%% 
% Set the similarities thresholds (values from the paper) in order: 

optimal_threshold = [0.6 ... chemical similarity
                     0.5 ... side effect similarity
                     0.7 ... Drug drug interaction
                     0 ... drug disease associations
                     0.3 ... sequence similarity
                     0.4 ... protein protein interaction similarity
                     0]; % disease similarity
%% 
% Remove values below threshold value

chemSim(chemSim <optimal_threshold(1)) = 0;
sideEffectSim(sideEffectSim < optimal_threshold(2)) = 0;
drugDrugSim(drugDrugSim < optimal_threshold(3)) = 0;
drugDiseaseSim(drugDiseaseSim<optimal_threshold(4)) = 0;
protSeqSimN(protSeqSimN< optimal_threshold(5)) = 0;
PPISim(PPISim<optimal_threshold(6)) = 0;
protDiseaseSim(protDiseaseSim<optimal_threshold(7)) = 0;
%% 10 runs of 10-fold crossvalidation for our method

SIR = {chemSim,sideEffectSim,drugDrugSim,drugDiseaseSim};
SIC = {protSeqSimN,PPISim,protDiseaseSim};  
initRange = [ones(1,length(SIC)+length(SIR)+3)*0.05]';
%set lower bound of weights for fmincon
low_bound = bound * ones(length(initRange),1);
tic
for j =1:10

    cost_fold = zeros(nfolds,1);
    localX = zeros(nfolds,length(initRange));
    for i = 1:nfolds
        %initialise weights for optimisation
        params =rand(length(initRange),1) .* initRange;
    
        local_R = R;
        test_fold = tests{j}{i};
        idx = sub2ind(size(R),test_fold(:,1),test_fold(:,2));
        local_R(idx) = 0;
        %initialise empty rows with the top k similar drugs
        local_R = top_k_combined_similarity(local_R,prior_SIM_row,top_k);

        SIC_tmp = SIC;
        SIR_tmp = SIR;
        
        
        Rr = normaliseY(local_R,1); % by row
        %compute row covariance similarity
        cov_no_diag = (Rr*Rr' - (Rr*Rr').*eye(size(Rr*Rr',1)));

        % SIR_tmp = [cov_no_diag SIR_tmp(:)'];
        SIR_tmp = [{cov_no_diag}, SIR_tmp(:)'];
        
        Rc = normaliseY(local_R,0); % by column   
        cov_no_diag = (Rc'*Rc - (Rc'*Rc) .* eye(size(Rc'*Rc,1)));
        
        %add covariance to drug similarites
        % SIC_tmp = [cov_no_diag SIC_tmp(:)'];
        SIC_tmp = [{cov_no_diag}, SIC_tmp(:)'];
        
        %get optimal values for gamma, alphas and betas with our method
        [localX(i,:), cost_fold(i)] =  fmincon(@(t)(our_method_gradients(t,local_R,Rr,Rc,SIC_tmp,SIR_tmp,lambda1,lambda2)),params,[],[],[],[],low_bound,[],[],options);
        [Res] = get_predictions( localX(i,:), local_R, SIR_tmp, SIC_tmp);

        our_method_predictions{j}{i} = Res;
    end
end
toc 

outdir = '../repository/cold_start_evaluation_1/drug_cold_start/our_method';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

save('-mat7-binary', fullfile(outdir, 'our_method_predictions_10CV_2026.mat'), 'our_method_predictions');

% save('../repository/cold_start_evaluation_1/drug_cold_start/our_method/our_method_predictions_10CV_2026.mat','our_method_predictions');
%%