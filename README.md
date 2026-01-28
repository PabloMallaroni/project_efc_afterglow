# project_efc_afterglow

Code and resources supporting the manuscript:

**Neurocomputational evidence of sustained Self–Other mergence after psychedelics**  
Mallaroni, P., Mason, N. L., Preller, K. H., Razi, A., Ereira, S., & Ramaekers, J. G.  
medRxiv preprint (doi: 10.1101/2025.10.07.25337510).

## Status

This repository accompanies a manuscript that is **under peer review**.

The analytical outputs are intended to reproduce the reported findings, but the codebase is being **progressively updated for legibility and structure**.

Practical implications:
- Script names, function boundaries, and directory layout may change as refactoring proceeds
- Hard-coded paths are being removed over time; for now, several scripts require editing a `main_path` or `BASE` variable
- Documentation is being expanded incrementally

## What is in this repository

The project integrates:
- Behavioural modelling of a probabilistic false-belief task (pFBT)
- Effective connectivity (spectral DCM, PEB) on peak-effect 7T resting-state fMRI
- Permutation-based multivariate statistics linking computational and neuroimaging predictors to (sub)acute psychosocial outcomes
- NeuroSynth-derived region definitions and peak MNI coordinates for the Theory-of-Mind network used in the DCM analyses

## Repository layout

```
project_efc_afterglow/
  neurosynth_maps/
    pub_getneuro.py
    neurosynth_coord.csv
    vmpfc/ dmpfc/ precuneus/ tpj/     # NeuroSynth association-test maps used to derive coordinates

  scripts/
    pFBT/
      model_fitting/                  # Probabilistic false-belief task model fitting utilities (MATLAB)
      synthetic_data/                 # Synthetic pFBT sessions used for demonstrations
      Parameter Estimates.mat         # Winning-model parameter estimates used by demo scripts
      TaskAccuracy.m                  # Demonstration of task-accuracy estimation (MATLAB)
      VisualiseParameters.m           # Reproduces key behavioural inference/plots (MATLAB)

    dcm/
      pub_dcmfirst.m                  # First-level spectral DCM specification/estimation (MATLAB/SPM)
      pub_dcmpeb.m                    # Within-subject PEB + third-level PEB-of-PEBs (MATLAB/SPM)
      plot_eFC_delta.m                # Plotting helper for effective connectivity contrasts (MATLAB)
      plot_eFC_mat_supplement.m       # Supplementary matrix plotting helper (MATLAB)

    wellbeing/
      pub_manova.py                   # λ predictor: perm MANOVA + canonical analysis + CV (Python)
      pub_manova_efc.py               # eFC predictor variant of the same pipeline (Python)
```

## Getting started

### 1) Behavioural modelling (pFBT)

The `scripts/pFBT/` folder contains self-contained demonstration scripts.

Full model fitting:
- The subfolder `scripts/pFBT/model_fitting/` contains the functions used to fit the nested model family described in the manuscript.

### 2) Effective connectivity (spectral DCM)

The DCM scripts assume an SPM-based workflow and a local project directory with:
- Preprocessed resting-state NIfTI files available per subject/session
- Corresponding BIDS JSON sidecars for acquisition metadata
- ROI definitions derived from the NeuroSynth coordinate file

Key scripts:
- `scripts/dcm/pub_dcmfirst.m` runs first-level GLMs, extracts ROI time series, and estimates a spectral DCM
- `scripts/dcm/pub_dcmpeb.m` aggregates session-level DCMs into within-subject PEBs and runs a third-level PEB-of-PEBs for group effects and associations with behavioural predictors

Important: these scripts currently contain **hard-coded paths** (e.g., `main_path`, `paths.spm`, ROI directories). You will need to edit them to match your local environment.

### 3) Multivariate wellbeing association (Python)

The Python pipelines in `scripts/wellbeing/` implement:
- Drug-aware within-subject residualisation
- Permutation MANOVA with within-subject shuffling
- Univariate follow-ups with Benjamini–Hochberg FDR
- Canonical analysis with structure coefficients
- Repeated-measures correlation
- K-fold cross-validation with permutation p-values

`pub_manova.py` uses `lambda_val` as the predictor.

`pub_manova_efc.py` uses an effective connectivity predictor specified by `PREDICTOR_COL` (default: `rtpjdmpfc`).

Important: these scripts currently assume local subject-level data excel inputs and paths:
- `significant_behaviour.xlsx`
- `significant_outcomes_all.xlsx`

You will need to update `BASE`, `EXCEL_PATH`, and `OUTDIR` to run them.

## Dependencies

### MATLAB
- MATLAB (tested in a modern MATLAB distribution)
- Statistics and Machine Learning Toolbox (for `fitlme`)

For DCM:
- SPM12 (or a compatible SPM build providing spectral DCM and PEB utilities)

### Python
The wellbeing scripts use:
- numpy, pandas
- scipy
- statsmodels
- scikit-learn
- matplotlib
- pingouin

The NeuroSynth utility uses:
- nibabel
- nilearn
- scipy
- seaborn

## Data availability

This repository includes:
- NeuroSynth maps and derived peak coordinates used for ROI definition
- Synthetic pFBT sessions used for demonstrations
- A `.mat` file containing winning-model parameter estimates used by the figure-generation demo scripts

Raw behavioural and neuroimaging data are not included here.

## Contact

For questions, issues, or requests related to analysis details, please use contact the corresponding authors listed in the manuscript.

