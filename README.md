# LI-DTI: A linear interpretable model for Drug Target Prediction

Official implementation of LI-DTI, a Linear Interpretable Drug-Target Interaction prediction model inspired by recommender systems.
LI-DTI learns from drug-drug and target-target similarity matrices –-chemical, biological, and pharmacological-– and provides interpretable predictions as a linear combination of these interactions.

# LI-DTI
<p align="center">
  <img src="images/model_figure.jpg" alt="LI-DTI model architecture" width="700">
</p>

# Requeriments
- MATLAB (R2018a or later recommended).
- Toolboxes: Statistics and Machine Learning and Optimization.

# Code and Data

# Data Files
All data files are available here: https://zenodo.org/records/XXXXXX

- data/luo_dataset: Contains the data which are used for the prediction tasks in DTINet dataset.
- repository: Contains the replicable predictions for each model in warm and cold start.
- prospective_evaluation: Contains the data for the prospective evaluation in DrugBank.

# Code
Here we describe the code used in LI-DTI.

# Warm and Cold-Start Scenarios
For Warm-Start Scenario:
1. Run warm_start_10_CV_predictions.mlx to get the predictions.
2. Run plot_warm_start.mlx to get the results.
3. Run removing_similarities_warm_start_true_CV.mlx to get the results for the case of removing similarities.

For Drug Cold-Start Scenario:
1. Run drug_cold_10_cv.mlx to get the predictions.
2. Run plot_10CV_drug_cold_start_all_metrics.mlx to get the results.

For Target Cold-Start Scenario:
1. Run target_cold_10_cv.mlx to get the predictions
2. Run plot_10CV_target_cold_start_all_metrics.mlx to get the results.

All predictions are saved in the ./repository/ folder. 

# Web Tool
Web tool to search for drugs and targets and generate sunburst plots that explain predicted scores available here: https://paccanarolab.org/exdtiweb/.

# Prospective Evaluation
1. Run prospective_2022_compute_allSI.mlx for saved the test in ./repository/prospective_evaluation.
2. For Recall-Drugs run corrections_of_prospective_2022_drugswise.mlx.
3. For Recall-Targets run corrections_of_prospective_2022_targetwise.mlx.



# Contacts
If you have any questions or comments, please feel free to contact:
- **Santiago Ferreyra** (`santiago.ferreyra@fgv.br`).
