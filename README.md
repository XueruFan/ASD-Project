# Profiling brain morphology for autism spectrum disorder with two cross-culture large-scale consortia

References
===========
+ Fan X-R, He Y, Wang Y-S, Li L, Lifespan Brain Chart Consortium (LBCC), China Autism Brain Imaging Consortium (CABIC), et al. (2025): Profiling brain morphology for autism spectrum disorder with two cross-culture large-scale consortia. bioRxiv. 2025.02.24.639771.

----

Background
====

We explore neurodevelopmental heterogeneity in Autism Spectrum Disorder (ASD) through normative modeling of cross-cultural cohorts. Leveraging large-scale datasets from Autism Brain Imaging Data Exchange (ABIDE) and China Autism Brain Imaging Consortium (CABIC), the model identifies two ASD subgroups with distinct brain morphological abnormalities: subgroup "L" is characterized by generally smaller brain region volumes and higher rates of abnormality, while subgroup "H" exhibits larger volumes with less pronounced deviations in specific areas. Key areas, such as the isthmus cingulate and transverse temporal gyrus, were identified as critical for subgroup differentiation and ASD trait correlations. In subgroup H, the regional volume of the isthmus cingulate cortex showed a direct correlation with individuals' autistic mannerisms, potentially corresponding to its slower post-peak volumetric declines during development. These findings offer insights into the biological mechanisms underlying ASD and support the advancement of subgroup-driven precision clinical practices.

Code Release and Usage
====

Before starting the analysis, we preprocessed the T1 MRI scans from the ABIDE (I and II) and CABIC datasets with FreeSurfer. The subsequent processing and analysis consist of XXX steps. Detailed explaination is in each script.

**Step 1: Extract brain measurement metrics**

We extracted (and calaulate) 7 global measures and 34 regional gray matter volumes from the preprocessing results.

- `Step1/extract_ABIDE1.R`
- `Step1/extract_ABIDE2.R`
- `Step1/extract_CABIC.R`

**Step 2: ComBat harmonization**

First combine demographic and MRI measurements, format the data based on LBCC standards for subsequent Combat harmonization and individual OoS calculation. Then do Combact.

- `Step2/prepare_Combat_ABIDE.R`

  Note 2-1: During this step, we found that one subject (ABIDE2 ID: 28793) had MRI data but no phenotypic data. This subject was from the GU site. Upon checking the phenotypic file, we noticed that the behavioral data for this subject was misaligned, likely due to a formatting error during data compilation. Since this subject was part of the control group, we decided to exclude it, as the impact would be minimal. Additionally, when merging ABIDE1 and ABIDE2 datasets, we observed inconsistencies in column names. To standardize the data, we aligned all column names with ABIDE2’s naming convention.

  Note 2-2: In the original LBCC work, TCV was calculated as GMV + WMV, so we adopted the same method here.
  
- `Step2/prepare_Combat_CABIC.R`
- `Step2/do_Combat_ABIDE.R`
- `Step2/do_Combat_CABIC.R`

**Step 3: Centile scoring**

We first calculated the OoS centile scores, then visualized the results and tested the normality of the scores with Jarque-Bera test.

- `Step3/calculate_Centile_ABIDE.R`
- `Step3/calculate_Centile_CABIC.R`

Note 3-1: We saved the results separately for different age groups to enable subsequent analyses on narrower age ranges.

- `Step3/Age-13/plot_Centile_ABIDE.R`
- `Step3/Age-13/plot_Centile_CABIC.R`
- `Step3/Age-13/norm_Centile_ABIDE.R`
- `Step3/Age-13/norm_Centile_CABIC.R`

**Step 4: Spectral clustering**

We performed spectral clustering analysis using 34 regional volume OoS scores as classification features to subgroup male ASD individuals (<13 years) from the ABIDE dataset. Tested the normality of the scores with Jarque-Bera test for each cluster. Based on the previous classification results, used SVM-RFECV to select the optimal feature set for classification. Then visualized the results for each cluster. At last, we calculated the extreme percentage of OoS centils for each cluster and project the results on the brain.

- `Step4/Age-13/SpectralCluster/do_Clustering_ABIDE.R`

  Note 4-1: Four participants were removed due to abnormal brain segmentation - three from ABIDE I (IDs: 50752, 51008, 51208) and one from ABIDE II (ID: 29109).
  
- `Step4/Age-13/SpectralCluster/norm_Centile_Clusters_ABIDE.R`
- `Step4/Age-13/SpectralCluster/select_Feature_ABIDE.py`
- `Step4/Age-13/SpectralCluster/plot_Clusters_ABIDE.R`

**Step 5: Analyze differences between clusters**

  We analyzed the differences of population (individuals used for clustering analysis) composition, demographic and cognitive behavioral, and OoS centile scores between clusters.

- `Step5/Age-13/SpectralCluster/statistic_Pheno_ABIDE.R`
- `Step5/Age-13/SpectralCluster/statistic_Difference_Scales_PartA_ABIDE.R`

  Note 5-1: For site difference analysis, only sites with total sample sizes >10 were included. For MRI scanner model/manufacturer comparisons, only those with >30 samples (which is also >10) were analyzed.

- `Step5/Age-13/SpectralCluster/statistic_Difference_Scales_PartB_ABIDE.R`
- `Step5/Age-13/SpectralCluster/arrange_Difference_Scales_ABIDE`
- `Step5/Age-13/SpectralCluster/norm_Scales_Clusters_ABIDE.R`

  Note 5-2: Test whether cognitive scores follow a normal distribution
  
- `Step5/Age-13/SpectralCluster/statistic_Difference_MRI_ABIDE.R`

**Step 6: Validate results using CABIC**

Use the previously trained SVM model to predict clusters in the CABIC dataset. Then did the same ststistic analysis like Step 5.

- `Step6/Age-13/SpectralCluster/predict_Clusters_CABIC.py`
- `Step6/Age-13/SpectralCluster/plot_Clusters_CABIC.R`
- `Step6/Age-13/SpectralCluster/statisic_Pheno_CABIC.R`
- `Step6/Age-13/SpectralCluster/statistic_Difference_Scales_CABIC.R`

**Step 7: Correlation analysis of OoS centile scores and cognitive behaviors**

We performed correlation analysis between brain and behavioral measures involving in the steps before.

- `Step7/Age-13/SpectralCluster/statistic_Correlations_Clusters_ABIDE.R`
- `Step7/Age-13/SpectralCluster/statistic_Correlations_Clusters_CABIC.R`
  
  Note 7-1: Control Site as a fixed effect and TCV as a covariate. Pearson's correlation was applied to most measures, while Spearman's correlation was used for ADOS-2 RRB scores due to their limited range.

- `Step7/Age-13/SpectralCluster/compare_Correlations_Clusters.R`
- `Step7/Age-13/SpectralCluster/plot_Correlations_Clusters_ABIDE.R`
- `Step7/Age-13/SpectralCluster/plot_Correlations_Clusters_CABIC.R`
- `Step7/Age-13/SpectralCluster/compare_Correlations_Clusters.R`

**Step 8: Structural covarianc analysis**

We analyzed structural covariance on OoS scores of each two DK regions. Then we performed correlation analysis on residuals obtained from linear regression models that controlled for site effects and TCV OoS scores.

- `Step7/Age-13/SpectralCluster/stru_Covariance_ABIDE.R`
- `Step7/Age-13/SpectralCluster/stru_Covariance_CABIC.R`
- `Step7/Age-13/SpectralCluster/compare_StruCovariance.R`
- `Step7/Age-13/SpectralCluster/statistic_Correlations_StruCova.R`

**SM-Step 1. Analyze age-related changes in OoS scores**

We use GAMM to explore the relationship between OoS scores of age and spectific brain regions.

- `SM-Step1/Age-13/SpectralCluster/gamm_Centile_Clusters_ABIDE.R`
- `SM-Step1/Age-13/SpectralCluster/gamm_Centile_Clusters_CABIC.R`
- `SM-Step1/Age-13/SpectralCluster/plot_gamm_Centile_ABIDE.R`
- `SM-Step1/Age-13/SpectralCluster/plot_gamm_Centile_CABIC.R`
- `SM-Step1/Age-13/SpectralCluster/plot_gamm_Centile_Both.R`

  Note SM1-1: Only subgroup H were plotted in this code.

**SM-Step 2. Plot with LBCC**

We ploted LBCC-corrected 34 regional volumes with GAM smoothing for clusters.

- `SM-Step1/Age-13/SpectralCluster/plot_LBCC_ABIDE.R`
- `SM-Step1/Age-13/SpectralCluster/plot_LBCC_CABIC.R`
- `SM-Step1/Age-13/SpectralCluster/plot_LBCC_Both.R`

  Note SM2-1: Only subgroup H were plotted in this code.

**SM-Analysis 1. Narrow age range analysis**

After calculated the OoS centile scores, we also performed the same analysis from clustering with a narrow age range 5.0~9.9 years. Related scripts are in the `Age510` folder within each step.

For example, 

- `Step4/Age510/SpectralCluster/do_Clustering_ABIDE.R`

To evaluate the consistency of cluster ID compare to the clustering analysis on broad (<13 years) and narrow (5~9.9 years) age range, use following script,

- `SM-Analysis1/compare_ClusterID.R`

**SM-Analysis 2. GMM clustering**

We performed another clustering method GMM analysis to evaluate the robustness of our results. Related scripts are in the `GmmCluster` folder within each step.

For example, 

- `Step4/Age-13/GmmCluster/do_Clustering_ABIDE.R`

**SM-Analysis 3. Indepentent clustering analysis on CABIC**

We performed spectral clustering on CABIC as did on ABIDE, to evaluate the consistency of cluster ID compare to predicting the cluster with ABIDE classifier.

- `SM-Analysis3/do_Clustering_CABIC.R`
- `SM-Analysis3/plot_Clusters_CABIC.R`
- `SM-Analysis3/compare_ClusterID.R`

**SM-Analysis 4. Clustering analysis on NYU site only**

We performed spectral clustering on ABIDE-NYU site, to evaluate the consistency of cluster ID compare with ABIDE classifier in main analysis.

- `SM-Analysis4/do_Clustering_CABIC.R`
- `SM-Analysis4/compare_ClusterID.R`
