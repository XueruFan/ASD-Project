# Profiling brain morphology for autism spectrum disorder with two cross-culture large-scale consortia

References
===========
+ Fan X-R, He Y, Wang Y-S, Li L, Lifespan Brain Chart Consortium (LBCC), China Autism Brain Imaging Consortium (CABIC), et al. (2025): Profiling brain morphology for autism spectrum disorder with two cross-culture large-scale consortia. bioRxiv. 2025.02.24.639771.

----

Background
====

We explore neurodevelopmental heterogeneity in Autism Spectrum Disorder (ASD) through normative modeling of cross-cultural cohorts. Leveraging large-scale datasets from Autism Brain Imaging Data Exchange (ABIDE) and China Autism Brain Imaging Consortium (CABIC), the model identifies two ASD subgroups with distinct brain morphological abnormalities: subgroup "L" is characterized by generally smaller brain region volumes and higher rates of abnormality, while subgroup "H" exhibits larger volumes with less pronounced deviations in specific areas. Key areas, such as the isthmus cingulate and transverse temporal gyrus, were identified as critical for subgroup differentiation and ASD trait correlations. In subgroup H, the regional volume of the isthmus cingulate cortex showed a direct correlation with individuals' autistic mannerisms, potentially corresponding to its slower post-peak volumetric declines during development. These findings offer insights into the biological mechanisms underlying ASD and support the advancement of subgroup-driven precision clinical practices.

Code Release
====

Before starting the analysis, we preprocessed the T1 MRI scans from the ABIDE (I and II) and CABIC datasets with FreeSurfer. The subsequent processing and analysis consist of XXX steps. Detailed explaination is in each script.

**Step 1: Extract brain measurement metrics**

We extracted (and calaulate) 7 global measures and 34 regional gray matter volumes from the preprocessing results.

- `Step1/extract_ABIDE1.R`
- `Step1/extract_ABIDE2.R`
- `Step1/extract_CABIC.R`

**Step 2: ComBat harmonization**
- `Step2/prepare_Combat_ABIDE.R`

  During this step, we found that one subject (ABIDE2 ID: 28793) had MRI data but no phenotypic data. This subject was from the GU site. Upon checking the phenotypic file, we noticed that the behavioral data for this subject was misaligned, likely due to a formatting error during data compilation. Since this subject was part of the control group, we decided to exclude it, as the impact would be minimal. Additionally, when merging ABIDE1 and ABIDE2 datasets, we observed inconsistencies in column names. To standardize the data, we aligned all column names with ABIDE2’s naming convention.
  
- `Step2/do_Combat_ABIDE.R`
- `Step2/prepare_Combat_CABIC.R`
- `Step2/do_Combat_CABIC.R`



