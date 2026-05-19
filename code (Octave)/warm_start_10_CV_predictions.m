%% *LI-DTI: Linear Interpretable Drug-Target Interaction* 
% *Warm Start Scenario*
%% 
% * Output: our_method_predictions_10CV_2026.mat file with the predictions of 
% LI-DTI for the DTINet dataset saved in ../repository/warm_start_10_CVs. 

%clear all;
%% Read the data
%% 
% Load DTI matrix

if exist('OCTAVE_VERSION', 'builtin')
    pkg load optim
    pkg load statistics
end

R = load('../data/luo_dataset/mat_drug_protein.txt');
%% 
% Load similarities

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

rng(0);        % For reproducibility 
nfolds = 10;   % Number of folds in the cv
bound = 0;     % Lower limit of the weights of the model
% Parameters for matlabs fmincon optimiser:

% Matlab Implementation
% options = optimoptions(@fmincon,'Algorithm','trust-region-reflective',...
%    'CheckGradients',false,'SpecifyObjectiveGradient',true,...
%    'MaxFunctionEvaluations',10000,'MaxIterations',10000,'Display','off');

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

lambda1 = 1;   % Lambda for the drug and target covariance similarity
lambda2 = 0.1; % Lambda for the remaining similarites
%% 
% Load the 10 CV splits

for k = 1:10
    for j = 1:10
        load(['../data/luo_dataset/warm_start_10_CVs/our_method/cv_' int2str(k) '/test' int2str(j) '.mat']);
        tests{k}{j} = test_fold;
    end
end
%% 
% Set the similarities thresholds (values from the paper) in order: 

optimal_threshold = [0.6 ... Chemical similarity
                     0.5 ... Side effect similarity
                     0.7 ... Drug drug interaction
                     0 ...   Drug disease associations
                     0.3 ... Sequence similarity
                     0.4 ... Protein protein interaction similarity
                     0];     % Disease similarity

%% 
% Remove values below threshold value

chemSim(chemSim <optimal_threshold(1)) = 0;
sideEffectSim(sideEffectSim < optimal_threshold(2)) = 0;
drugDrugSim(drugDrugSim < optimal_threshold(3)) = 0;
drugDiseaseSim(drugDiseaseSim<optimal_threshold(4)) = 0;
protSeqSimN(protSeqSimN< optimal_threshold(5)) = 0;
PPISim(PPISim<optimal_threshold(6)) = 0;
protDiseaseSim(protDiseaseSim<optimal_threshold(7)) = 0;
%% 
% Group similarities 

sideInfoRow_no_diag = {chemSim,sideEffectSim,drugDrugSim,drugDiseaseSim};
sideInfoCol_no_diag = {protSeqSimN,PPISim,protDiseaseSim};
%% 
% Variable to initialise weights 

mParams = zeros(10,10);
%% 10 runs of 10-fold crossvalidation for LI-DTI

tic
for i = 1:10
    %variable to hold cost value and optimal params
    cost_fold = zeros(nfolds,1);
    localX = zeros(1,10);
    for j=1:10
            local_R = R;
            %remove testing DTIs from training
            test_fold = tests{i}{j};
            idx = sub2ind(size(R),test_fold(:,1),test_fold(:,2));
            local_R(idx) = 0;
            % Matlab Implementation
            %SIC_tmp = sideInfoCol_no_diag;
            %SIR_tmp = sideInfoRow_no_diag;

            SIR_tmp = sideInfoRow_no_diag;
            SIC_tmp = sideInfoCol_no_diag;

            %initialise weights for optimisation
            initRange = [ones(1,length(SIC_tmp)+length(SIR_tmp)+3)*0.05]';
            params = rand(length(initRange),1) .* initRange;

            %set lower bound of weights for fmincon
            low_bound = bound * ones(length(params),1);
            
            %normalise DTI matrix by row
            Rr = normaliseY(local_R,1); % by row

            %compute row covariance similarity
            cov_no_diag = (Rr*Rr' - (Rr*Rr').*eye(size(Rr*Rr',1)));

            %add covariance to drug similarites
            %SIR_tmp = [cov_no_diag SIR_tmp(:)'];
            SIR_tmp = [{cov_no_diag}, SIR_tmp(:)'];
                        
            %normalise DTI matrix by column
            Rc = normaliseY(local_R,0); % by column   
            cov_no_diag = (Rc'*Rc - (Rc'*Rc) .* eye(size(Rc'*Rc,1)));
           
            %add covariance to drug similarites
            %SIC_tmp = [cov_no_diag SIC_tmp(:)'];
            SIC_tmp = [{cov_no_diag}, SIC_tmp(:)'];

            %get optimal values for gamma, alphas and betas with our method
            [local_X,~] =  fmincon(@(t)(our_method_gradients(t,local_R,Rr,Rc,SIC_tmp,SIR_tmp,lambda1,lambda2)),params,[],[],[],[],low_bound,[],[],options);
            
            %get predictions for fold with optimal values
            [Res] = get_predictions( local_X, local_R, SIR_tmp, SIC_tmp);
        
        our_method_predictions{i}{j} = Res;
    end
end
toc
%%
%save('../repository/warm_start_10_CVs/our_method/our_method_predictions_10CV_2026.mat','our_method_predictions');

outdir = '../repository/warm_start_10_CVs/our_method';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end
%save(fullfile(outdir, 'our_method_predictions_10CV_2026.mat'), 'our_method_predictions');
save('-mat7-binary', fullfile(outdir, 'our_method_predictions_10CV_2026.mat'), 'our_method_predictions');
