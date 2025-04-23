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

Before starting the analysis, we preprocessed the T1 MRI scans from the ABIDE (I and II) and CABIC datasets with FreeSurfer.

The subsequent processing and analysis consist of XXX steps:

**Step 1: Extract brain measurement metrics**
We extracted 7 global measures and 34 regional gray matter volumes from the preprocessing results.
