%% *LI-DTI: Linear Interpretable Drug-Target Interaction* 
% *Prospective Evaluation*
%% 
% * Output: our_method_prospective.mat file with the predictions of LI-DTI for 
% the prospective dataset saved in ../repository/prospective_evaluation.
%% Initialization 

% clear all; close all; clc;

if exist('OCTAVE_VERSION', 'builtin')
    pkg load optim
    pkg load statistics
end

%% Read the data

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
%%
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
% Load the optimal threshold values (from the paper)

load('../repository/warm_start_10_CVs/our_method/optimal_threshold.mat');
%%
chemSim(chemSim <optimal_threshold(1)) = 0;
sideEffectSim(sideEffectSim < optimal_threshold(2)) = 0;
drugDrugSim(drugDrugSim < optimal_threshold(3)) = 0;
drugDiseaseSim(drugDiseaseSim<optimal_threshold(4)) = 0;
protSeqSimN(protSeqSimN< optimal_threshold(5)) = 0;
PPISim(PPISim<optimal_threshold(6)) = 0;
protDiseaseSim(protDiseaseSim<optimal_threshold(7)) = 0;
%% 
% Set setting

rng(0);      % For reproducibility 
nfolds = 10; % Number of folds in the cv
bound = 0;   % Lower limit of the weights of the model
% % Parameters for matlabs fmincon optimiser

% MATLAB IMPLEMENTATION
% options = optimoptions(@fmincon,'Algorithm','trust-region-reflective',...
%     'CheckGradients',false,'SpecifyObjectiveGradient',true,'MaxFunctionEvaluations',10000,'MaxIterations',10000);
% %%

% OCTAVE IMPLEMENTATION
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


sideInfoRow_no_diag = {chemSim,sideEffectSim,drugDrugSim,drugDiseaseSim};
sideInfoCol_no_diag = {protSeqSimN,PPISim,protDiseaseSim};
SIC = [sideInfoCol_no_diag(:)'];
SIR = [sideInfoRow_no_diag(:)'];

lambda1 = 1;   % Lambda for the drug and target covariance similarity
lambda2 = 0.1; % Lambda for the remaining similarites
%% 
% Set setting

Rr = normaliseY(R,1); % by row
Rc = normaliseY(R,0); % by column
cov_no_diag = (Rr*Rr' - (Rr*Rr').*eye(size(Rr*Rr',1)));
% cov_no_diag(cov_no_diag > optimal_threshold_cov(1)) = 0;

% SIR = [cov_no_diag SIR(:)'];
SIR = [{cov_no_diag}, SIR(:)'];

cov_no_diag = (Rc'*Rc - (Rc'*Rc) .* eye(size(Rc'*Rc,1)));
% cov_no_diag(cov_no_diag > optimal_threshold_cov(2)) = 0;


% SIC = [cov_no_diag SIC(:)'];
SIC = [{cov_no_diag}, SIC(:)'];

initRange = [ones(1,length(SIC)+length(SIR)+1)*0.05]';
params =rand(length(initRange),1) .* initRange;
low_bound = bound * ones(length(params),1);       
%%
[localX_eq_1, all_costs_eq_1]  = fmincon(@(t)(our_method_gradients(t,R,Rr,Rc,SIC,SIR,lambda1,lambda2)),params,[],[],[],[],low_bound,[],[],options);
[Res] = get_predictions(localX_eq_1, R, SIR, SIC);
%%
localX_eq_1
%%

outdir = '../repository/prospective_evaluation';
if ~exist(outdir, 'dir')
    mkdir(outdir);
end

save('-mat7-binary', fullfile(outdir, 'our_method_prospective.mat'), 'Res');


% save('../repository/prospective_evaluation/our_method_prospective.mat','Res');
