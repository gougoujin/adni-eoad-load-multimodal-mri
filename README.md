# ADNI EOAD–LOAD Multimodal MRI Code Package

This repository-ready package contains a curated subset of analysis and figure-generation scripts supporting a multimodal MRI study comparing early-onset and late-onset Alzheimer's disease using ADNI-derived data.

## Important note
Raw imaging and clinical data are **not** included here. They were obtained from the Alzheimer's Disease Neuroimaging Initiative (ADNI) and should be accessed directly from ADNI under its registration and data-use policies.

## Directory overview
- `code/analysis/`: core analysis scripts used to build structural connectivity, compute structural-functional coupling, compute structural decoupling index, run ROI-wise diffusion analyses, and fit group models.
- `code/figures/`: scripts used to generate manuscript and supplementary figure panels from processed outputs.
- `code/utils/`: lightweight helper functions reused by analysis or figure scripts.
- `optional_notes/`: optional or alternative scripts retained for transparency, QC, or exploratory workflows.
- `SCRIPT_INDEX.tsv`: compact script inventory with categories and one-line purpose descriptions.
- `PIPELINE_OVERVIEW.md`: suggested execution order and dependencies.
- `environment_or_versions.txt`: software versions used in the study environment.

## Recommended usage
1. Review `PIPELINE_OVERVIEW.md`.
2. Edit local input/output paths in the analysis scripts.
3. Run the core analysis scripts in dependency order.
4. Generate figures from the resulting summary tables and matrices.

## Suggested minimal public release contents
For a lean public repository, prioritize:
- the scripts in `code/analysis/`
- the scripts in `code/figures/`
- a short `README.md`
- `environment_or_versions.txt`
- selected non-identifiable derived summary tables

## Reproducibility reminder
These scripts were lightly sanitized for sharing. You will likely need to:
- update local paths
- confirm toolbox availability
- adjust filenames to match your local outputs
- document any manual preprocessing or ROI selection steps in the repository README


## Programming languages
Core neuroimaging analyses in this repository are primarily implemented in MATLAB.
A small number of Python scripts are also included for image assembly and for
generating manuscript-oriented summary figures and tables.

## License
This repository is distributed under the MIT License. See `LICENSE`.

## Added Python utilities
- `code/python/assemble_main_figure_panel.py`
- `code/python/stack_graph_summary_panels.py`
- `code/python/assemble_alff_falff_node_panels.py`
- `code/python/generate_additional_manuscript_outputs.py`
- `code/python/generate_clean_manuscript_outputs.py`
- `code/python/generate_selected_supplementary_outputs.py`
