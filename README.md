# LECA Proteomics

## PPI ML pipeline

```bash
# if setting up a new project directory:
cd $proj_dir
mkdir {data,results,logs,cmds}
mkdir data/{feature_files,feature_subsets,feature_matrix}
mkdir results/{tpot,feature_selection,ppi_predict,walktrap}
tree -d * # confirm directory structure

# ---------------------------------------------------
# build feature matrix
# ---------------------------------------------------

# pickle feats (run time ~5min)
[ ! -e cmds/pickle_feat_files.sh ] || rm cmds/pickle_feat_files.sh; touch cmds/pickle_feat_files.sh
for x in `ls --ignore=*.pkl ${proj_dir}/data/features/`; do echo "python3 ${script_dir}/make_pkl.py -i ${proj_dir}/data/features/${x}"; done >> cmds/pickle_feat_files.sh
cat cmds/pickle_feat_files.sh | parallel -j24

# combine subsets (run time ~45min)
for dir in all chicken dolphin mouse pig rabbit; do echo \
"python3 ${script_dir}/build_featmat.py \
--directory ${proj_dir}/data/feature_files/${dir}/ \
--pickle --outfile_name ${proj_dir}/data/feature_subsets/${dir}_features"; done > cmds/join_subsets.sh
cat cmds/join_subsets.sh | parallel -j6 2>&1 | tee ${proj_dir}/logs/join_subsets.log

# combine all (run time ~20min)
python3 ${script_dir}/build_featmat.py --directory ${proj_dir}/data/feature_subsets/ \
--pickle --outfile_name ${proj_dir}/data/feature_matrix/feature_matrix 2>&1 | tee ${proj_dir}/logs/build_featmat.log

# ---------------------------------------------------
# label feature matrix
# ---------------------------------------------------

# label positive/negative PPIs (run time ~2 hours)
# seed makes negative PPI labels reproducible
python3 ${script_dir}/label_featmat.py \
--featmat ${proj_dir}/data/featmats/featmat.pkl \
--gold_std_file ${proj_dir}/data/gold_stds/all.gold.cmplx.noRibos.merged.txt \
--outfile_name ${proj_dir}/data/featmats/featmat_labeled --seed 13 2>&1 | tee ${proj_dir}/logs/label_featmat.log

# ---------------------------------------------------
# get optimized ML pipeline
# ---------------------------------------------------

# tpot "generations" and "pop_size" set here to default cfmsflow params
# seed makes group splits, train/test splits, & model reproducible
# extremely slow (will probably run overnight)
# to speed things up reduce num_splits, generations, & pop_size
python3 ${script_dir}/run_tpot.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled_traintest.pkl \
--outdir ${proj_dir}/results/tpot/gss/ \
--group_split_method GroupShuffleSplit --num_splits 5 --train_size 0.75 \
--generations 10 --pop_size 20 --seed 17 2>&1 | tee ${proj_dir}/logs/run_tpot_gss.log

python3 ${script_dir}/run_tpot.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled_traintest.pkl \
--outdir ${proj_dir}/results/tpot/gkfold/ \
--group_split_method GroupKFold --num_splits 5 --train_size 0.75 \
--generations 10 --pop_size 20 --seed 17 2>&1 | tee ${proj_dir}/logs/run_tpot_gkfold.log

python3 ${script_dir}/run_tpot.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled_traintest.pkl \
--outdir ${proj_dir}/results/tpot/sgkfold/ \
--group_split_method StratifiedGroupKFold --num_splits 5 --train_size 0.75 \
--generations 10 --pop_size 20 --seed 13 2>&1 | tee ${proj_dir}/logs/run_tpot_sgkfold.log

# get all test scores
grep --with-filename 'Test set score' ../../logs/run_tpot*log > all_test_scores.txt

# get all pipelines
sed -n '/Average/,/Fix random/p' gss/tpot_pipeline_1

[ ! -e all_pipelines.txt ] || rm all_pipelines.txt; touch all_pipelines.txt
for dir in gkfold gss sgkfold; do 
        echo "-------------------------------------------------------"
        echo "-------------------- ${dir} ---------------------------" 
        echo "-------------------------------------------------------"
        sed -n '/Average/,/Fix random/p' ${dir}/tpot_pipeline* |
        sed '/# Fix random state for all the steps in exported pipeline/d' |
        sed '/# Fix random state in exported estimator/d' |
        sed 's/# Average/\n# Average/'
done >> all_pipelines.txt

# best models:
gkfold/tpot_model_3  # Linear SVC (3 steps)
gkfold/tpot_model_5  # SGD (3 steps)
gss/tpot_model_1     # KNN (3 steps)
gss/tpot_model_4     # ExtraTrees (1 step)

# ---------------------------------------------------
# run recursive feature elimination
# ---------------------------------------------------

# choose best tpot models as input
# note: KNN does not provide feature selection logic

# seed makes group splits, train/test splits, & model reproducible
# highly variable run time (will take 0.25-6+ hours depending on # of features)
# to speed things up increase remove_per_step value

# [[ LinearSVC ]]
python3 ${script_dir}/select_features.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled_traintest.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_3.pkl \
--group_split_method GroupKFold \
--outdir ${proj_dir}/results/feature_selection/ \
--threads 12 --remove_per_step 3 --num_splits 5 --seed 17 2>&1 | tee ${proj_dir}/logs/select_features_SVC.log

# [[ SGDClassifier ]]
python3 ${script_dir}/select_features.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled_traintest.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_5.pkl \
--group_split_method GroupKFold \
--outdir ${proj_dir}/results/feature_selection/ \
--threads 12 --remove_per_step 3 --num_splits 5 --seed 17 2>&1 | tee ${proj_dir}/logs/select_features_SGD.log

# [[ ExtraTrees ]]
python3 ${script_dir}/select_features.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled_traintest.pkl \
--model ${proj_dir}/results/tpot/gss/tpot_model_4.pkl \
--group_split_method GroupKFold \
--outdir ${proj_dir}/results/feature_selection/ \
--threads 12 --remove_per_step 3 --num_splits 5 --seed 17 2>&1 | tee ${proj_dir}/logs/select_features_ExtraTrees.log

# ---------------------------------------------------
# predict PPIs
# ---------------------------------------------------

# seed makes group splits, train/test splits, & model reproducible
# fairly quick (run time ~5-15min)

script_dir="/stor/work/Marcotte/project/rmcox/leca/scripts"
proj_dir="/stor/work/Marcotte/project/rmcox/leca/ppi_ml"

# with seed=777 & train_size=0.5, actual train/test split is 75.6/24.4

# [[ LinearSVC ]]
# generate predictions using all features
python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_3.pkl \
--outdir ${proj_dir}/results/ppi_predict/feature_sweep/all/ \
--fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/predict_ppis_LinearSVC_all_feats.log

# generate predictions across a sweep of features
for x in 5 10 25 50 100 250; do echo "python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_3.pkl \
--outdir ${proj_dir}/results/ppi_predict/feature_sweep/${x} \
--feature_selection ${proj_dir}/results/feature_selection/top_feats_LinearSVC_RFECV_all.csv \
--num_feats ${x} --fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/predict_ppis_LinearSVC_${x}_feats.log"; done > cmds/sweep_feats_LinearSVC.sh
cat cmds/sweep_feats_LinearSVC.sh | parallel -j5


# [[ SGDClassifier ]]
# generate predictions using all features
python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_5.pkl \
--outdir ${proj_dir}/results/ppi_predict/feature_sweep/all/ \
--fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/predict_ppis_SGDClassifier_all_feats.log

# generate predictions across a sweep of features
for x in 5 10 25 50 100 250; do echo "python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_5.pkl \
--outdir ${proj_dir}/results/ppi_predict/feature_sweep/${x} \
--feature_selection ${proj_dir}/results/feature_selection/top_feats_SGDClassifier_RFECV_all.csv \
--num_feats ${x} --fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/predict_ppis_SGDClassifier_${x}_feats.log"; done > cmds/sweep_feats_SGDClassifier.sh
cat cmds/sweep_feats_SGDClassifier.sh | parallel -j5


# [[ ExtraTrees ]]
# generate predictions using all features
python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gss/tpot_model_4.pkl \
--outdir ${proj_dir}/results/ppi_predict/feature_sweep/all/ \
--fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/predict_ppis_ExtraTreesClassifier_all_feats.log

# generate predictions across a sweep of features
for x in 5 10 25 50 100 250; do echo "python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gss/tpot_model_4.pkl \
--outdir ${proj_dir}/results/ppi_predict/feature_sweep/${x} \
--feature_selection ${proj_dir}/results/feature_selection/top_feats_ExtraTreesClassifier_RFECV_all.csv \
--num_feats ${x} --fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/predict_ppis_ExtraTreesClassifier_${x}_feats.log"; done > cmds/sweep_feats_ExtraTreesClassifier.sh
cat cmds/sweep_feats_ExtraTreesClassifier.sh | parallel -j5

# ---------------------------------------------------
# run cross-validation
# ---------------------------------------------------

# [[ ExtraTrees ]]
python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gss/tpot_model_4.pkl \
--outdir ${proj_dir}/results/cross_val \
--feature_selection ${proj_dir}/results/feature_selection/top_feats_ExtraTreesClassifier_RFECV_all.csv \
--group_split_method GroupKFold --num_splits 5 \
--num_feats 25 --fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/crossval_ExtraTreesClassifier_25_feats.log

# [[ SGDClassifier ]]
python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_5.pkl \
--outdir ${proj_dir}/results/cross_val \
--feature_selection ${proj_dir}/results/feature_selection/top_feats_SGDClassifier_RFECV_all.csv \
--group_split_method GroupKFold --num_splits 5 \
--num_feats 25 --fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/crossval_SGDClassifier_25_feats.log

# [[ LinearSVC ]]
python3 ${script_dir}/predict_ppis.py \
--featmat ${proj_dir}/data/featmats/featmat_labeled.pkl \
--model ${proj_dir}/results/tpot/gkfold/tpot_model_3.pkl \
--outdir ${proj_dir}/results/cross_val \
--feature_selection ${proj_dir}/results/feature_selection/top_feats_LinearSVC_RFECV_all.csv \
--group_split_method GroupKFold --num_splits 5 \
--num_feats 100 --fdr_cutoff 0.1 --seed 777 --train_size 0.5 2>&1 | tee ${proj_dir}/logs/crossval_LinearSVC_100_feats.log

# ---------------------------------------------------
# cluster PPIs
# ---------------------------------------------------

# best models:
# ExtraTrees (100 feats) -> ppi_predict/feature_sweep/100/scored_interactions_fdr10_ExtraTreesClassifier.csv
# SGD (100 feats) -> ppi_predict/feature_sweep/100/scored_interactions_fdr10_SGDClassifier.csv

# seed makes results reproducible
# very fast (run time ~1min)

# [[ LinearSVC ]]
for x in 3 4 5 6 7; do echo "python3 ${script_dir}/detect_communities.py --scores ${proj_dir}/results/ppi_predict/feature_sweep/100/scored_interactions_fdr10_LinearSVC.csv --steps ${x} --outfile ${proj_dir}/results/walktrap/step_sweep/LinearSVC_100feats_fdr10_${x}steps --annotations ${proj_dir}/annotations/leca_eunog_annots_complete.030721.csv --seed 777 2>&1 | tee ${proj_dir}/logs/walktrap_LinearSVC_100feats_fdr10_${x}steps.log"; done > cmds/sweep_walktrap_steps_LinearSVC_100feats_fdr10.sh
cat cmds/sweep_walktrap_steps_LinearSVC_100feats_fdr10.sh | parallel -j5

# calculate precision-recall per cut
for x in 3 4 5 6 7; do echo "python3 ${script_dir}/calc_cluster_pr.py --cluster_file ${proj_dir}/results/walktrap/step_sweep/LinearSVC_100feats_fdr10_${x}steps.csv --positive_ppi_dict ${proj_dir}/data/featmats/positive_ppi_dict.pkl --negative_ppi_dict ${proj_dir}/data/featmats/negative_ppi_dict.pkl --outfile ${proj_dir}/results/walktrap/step_sweep/LinearSVC_100feats_fdr10_${x}steps_precision_recall.csv 2>&1 | tee ${proj_dir}/logs/walktrap_LinearSVC_100feats_fdr10_${x}steps_pr.log"; done > cmds/walktrap_pr_LinearSVC_100feats_fdr10.sh
cat cmds/walktrap_pr_LinearSVC_100feats_fdr10.sh | parallel -j5

# [[ ExtraTrees ]]
for x in 3 4 5 6 7; do echo "python3 ${script_dir}/detect_communities.py --scores ${proj_dir}/results/ppi_predict/feature_sweep/25/scored_interactions_fdr10_ExtraTreesClassifier.csv --steps ${x} --outfile ${proj_dir}/results/walktrap/step_sweep/ExtraTreesClassifier_25feats_fdr10_${x}steps --annotations ${proj_dir}/annotations/leca_eunog_annots_complete.030721.csv --seed 777 2>&1 | tee ${proj_dir}/logs/walktrap_ExtraTreesClassifier_25feats_fdr10_${x}steps.log"; done > cmds/sweep_walktrap_steps_ExtraTreesClassifier_25feats_fdr10.sh
cat cmds/sweep_walktrap_steps_ExtraTreesClassifier_25feats_fdr10.sh | parallel -j5

# calculate precision-recall per cut
for x in 3 4 5 6 7; do echo "python3 ${script_dir}/calc_cluster_pr.py --cluster_file ${proj_dir}/results/walktrap/step_sweep/ExtraTreesClassifier_25feats_fdr10_${x}steps.csv --positive_ppi_dict ${proj_dir}/data/featmats/positive_ppi_dict.pkl --negative_ppi_dict ${proj_dir}/data/featmats/negative_ppi_dict.pkl --outfile ${proj_dir}/results/walktrap/step_sweep/ExtraTreesClassifier_25feats_fdr10_${x}steps_precision_recall.csv 2>&1 | tee ${proj_dir}/logs/walktrap_ExtraTreesClassifier_25feats_fdr10_${x}steps_pr.log"; done > cmds/walktrap_pr_ExtraTreesClassifier_25feats_fdr10.sh
cat cmds/walktrap_pr_ExtraTreesClassifier_25feats_fdr10.sh | parallel -j5

# [[ SGDClassifier ]]
for x in 3 4 5 6 7; do echo "python3 ${script_dir}/detect_communities.py --scores ${proj_dir}/results/ppi_predict/feature_sweep/25/scored_interactions_fdr10_SGDClassifier.csv --steps ${x} --outfile ${proj_dir}/results/walktrap/step_sweep/SGDClassifier_25feats_fdr10_${x}steps --annotations ${proj_dir}/annotations/leca_eunog_annots_complete.030721.csv --seed 777 2>&1 | tee ${proj_dir}/logs/walktrap_SGDClassifier_25feats_fdr10_${x}steps.log"; done > cmds/sweep_walktrap_steps_SGDClassifier_25feats_fdr10.sh
cat cmds/sweep_walktrap_steps_SGDClassifier_25feats_fdr10.sh | parallel -j5


# calculate precision-recall per cut
for x in 3 4 5 6 7; do echo "python3 ${script_dir}/calc_cluster_pr.py --cluster_file ${proj_dir}/results/walktrap/step_sweep/SGDClassifier_25feats_fdr10_${x}steps.csv  --positive_ppi_dict ${proj_dir}/data/featmats/positive_ppi_dict.pkl --negative_ppi_dict ${proj_dir}/data/featmats/negative_ppi_dict.pkl --outfile ${proj_dir}/results/walktrap/step_sweep/SGDClassifier_25feats_fdr10_${x}steps_precision_recall.csv 2>&1 | tee ${proj_dir}/logs/walktrap_SGDClassifier_25feats_fdr10_${x}steps_pr.log"; done > cmds/walktrap_pr_SGDClassifier_25feats_fdr10.sh
cat cmds/walktrap_pr_SGDClassifier_25feats_fdr10.sh | parallel -j5

# ---------------------------------------------------
# sanity check the pipeline
# ---------------------------------------------------

# pretty quick (run time ~5min)
# will make a nice table in the log file
python3 ${script_dir}/sanity_checks.py \
--featmat_file ${proj_dir}/data/feature_matrix/feature_matrix_labeled.pkl \
--results_file ${proj_dir}/results/ppi_predict/all_feats/scored_interactions_all_ExtraTreesClassifier.csv 2>&1 | tee /project/rmcox/misc/vy.pipeline/logs/sanity_checks.log

# can also check group merge stats
python3 ${script_dir}/label_featmat.py \
--featmat ${proj_dir}/data/feature_matrix/feature_matrix.pkl \
--gold_std_file /project/vyqtdang/MS_brain/nextflow_brain/cfmsflow/gold_std/verNOG/CORUM4.1_core_merged.txt \
--outfile_name ${proj_dir}/data/feature_matrix/feature_matrix_labeled_allgroups \
--keep_cmplx_overlap --seed 13 2>&1 | tee ${proj_dir}/logs/label_featmat_allgroups.log

python3 ${script_dir}/sanity_checks.py \
--featmat_file ${proj_dir}/data/feature_matrix/feature_matrix_labeled_allgroups_traintest.pkl \
--results_file ${proj_dir}/results/ppi_predict/all_feats/scored_interactions_all_ExtraTreesClassifier.csv 2>&1 | tee ${proj_dir}/logs/sanity_checks_groups.log
```
