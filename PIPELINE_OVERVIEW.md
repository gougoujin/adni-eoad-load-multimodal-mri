# Pipeline overview

This file gives a practical execution order for the curated scripts.

## 1. Structural connectivity and subject matching
1. `code/analysis/build_subject_structural_connectivity_from_probtrackx.m`
2. `code/analysis/build_cohort_structural_connectivity_from_probtrackx.m`
3. `code/analysis/match_structural_and_functional_subjects.m`

## 2. Structural-functional coupling
4. `code/analysis/compute_global_sc_fc_coupling.m`
5. `code/analysis/compute_regional_sc_fc_coupling.m`
6. `code/analysis/fit_global_coupling_glm.m`
7. `code/analysis/fit_regional_coupling_glm.m`

## 3. Structural decoupling index
8. `code/analysis/compute_sdi_for_subject.m`
9. `code/analysis/compute_sdi_for_cohort.m`
10. `code/analysis/fit_regional_sdi_glm.m`

## 4. Diffusion ROI analyses
11. `code/analysis/extract_roi_diffusion_metrics.m`
12. `code/analysis/fit_roi_diffusion_glm.m`

## 5. Graph summary and multimodal downstream models
13. `code/analysis/summarize_gretna_bonferroni_results.m`
14. `code/analysis/run_coupling_sensitivity_and_clinical_models.m`
15. `code/analysis/run_multimodal_clinical_association_models.m`
16. `code/analysis/run_patient_only_age_sensitivity_analyses.m`

## 6. Figure generation
Run the scripts in `code/figures/` after their required summary tables or matrices are available.

## Optional scripts
The files in `optional_notes/` are retained for transparency and QC but are not required for the core manuscript pipeline.


## Optional Python post-processing
After running the core MATLAB analyses and exporting intermediate figure panels or summary tables,
the Python scripts in `code/python/` may be used for:
- assembling composite figure panels from existing images;
- generating clean manuscript-oriented tables/figures from processed text tables;
- exporting selected supplementary outputs from curated summary files.

These Python scripts are downstream helpers and are not required for the core neuroimaging analyses.
