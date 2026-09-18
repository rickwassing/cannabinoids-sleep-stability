# `code/` — Analysis pipeline and traceability

This document maps the analysis code to the manuscript (`documents/manuscript.md`) so that any
result can be traced from the raw EEG recording through subject-level processing, first-level
outcomes, and group-level statistics. See `code/REFACTOR_PLAN.md` for the full history of how
this code was reorganised from its original flat layout.

## Pipeline overview

`main.m` is the top-level orchestration script. It:
1. Initializes the environment (`css_init.m`) and adds `code/` (and subfolders) to the MATLAB path.
2. Runs subject-level **PROCESSING** steps via `pipeline/run_processing_section.m`, which dispatches
   to functions in `subject-level/` (preprocessing, segmentation, spectral-feature extraction) and
   `firstlevel-outcomes/` (functions that also compute and save a first-level outcome file) using a
   `Proc`/`Files`/`cfg` pattern.
3. Runs the group-level **ANALYSE** calls (one function per manuscript aim/figure), each living in
   its own subfolder of `group-level/`.

Derivative folders (gitignored, data lives on the `sleep` volume, not in this repo):
`derivatives/EEG-preproc` → `EEG-processed` → `EEG-segmented` → `EEG-output-fstlvl` (first-level
outcomes, `.mat` files) → consumed by `code/group-level/*` together with externally-computed PALM
cluster-statistics in the root-level `/group-level/<name>/` data folder (not to be confused with
the `code/group-level/` code folder documented here).

## First-level outcomes table

Each row is a function in `code/firstlevel-outcomes/` that computes a per-subject/per-session
outcome measure and saves it as a `.mat` file into `derivatives/EEG-output-fstlvl/` (via
`css_createfstlvloutput`), for later use by group-level analyses.

| Function | Derived from | Unit / interpretation | Saved to (derivative folder + filename pattern) | Consumed by |
|---|---|---|---|---|
| `css_infraslowfluctpowerspect.m` | 300-second (or longer) NREM stage-2 bouts of EEG sigma-power time series (`EEG-segmented/*_desc-sigmanrembout*.set`) | Per-channel Gaussian-fit parameters of the infraslow fluctuation (ISF) power spectrum: mean frequency (Hz), peak amplitude (a.u.), bandwidth (Hz), plus binned relative power in adjacent frequency bands (via `css_extractfeatures.m`) | `EEG-output-fstlvl/sub-<sub>/ses-<ses>/sub-<sub>_ses-<ses>_task-psg_desc-a1cnormsigma_fstlvl.mat` (normalized) and `..._desc-a1cabssigma_fstlvl.mat` (absolute) | `group-level/aim1_isf_characterization/analyse_isf_topography.m` and `plot_isf_topography.m` (Figure 2D, F–G); `group-level/aim1b_filter_validation/analyse_filter_edge_artefact.m` (Supplementary Figs S3–S4) |
| `css_crosscorr.m` | Same 300-second NREM sigma-power bouts, cross-correlated against the co-registered instantaneous heart-rate (HR) time series (`EEG-segmented/*_desc-hrnrembout_hr.set`) | Per-channel cross-correlation features: lag of peak correlation (s), peak correlation amplitude (r), and full-width-half-maximum of the cross-correlation peak (s) | `EEG-output-fstlvl/sub-<sub>/ses-<ses>/sub-<sub>_ses-<ses>_task-psg_desc-nremboutxcorr120s_fstlvl.mat` | `group-level/aim1_isf_characterization/plot_isf_topography.m` (Figure 2C, E, H) |
| `css_extractfeatures.m` | Helper called internally by `css_infraslowfluctpowerspect.m` (not run standalone) | Packages the Gaussian-fit outputs (`FIT`) and raw ISF spectrum (`ISF`) into the `features` struct array format expected by `css_createfstlvloutput.m` | n/a — intermediate helper, not a file-writing function itself | `css_infraslowfluctpowerspect.m` |
| `css_createfstlvloutput.m` | Helper called by every first-level-outcome-producing function above (and by `group-level/aim1_isf_characterization/analyse_sigma_spindle_similarity.m` for sigma-power/spindle-density features) | n/a — writes the standard first-level output `.mat` structure (channel locations, features, metadata) to disk | `derivatives/EEG-output-fstlvl/sub-<sub>/ses-<ses>/<filename>.mat` (filename determined by caller) | All group-level analyses that read `EEG-output-fstlvl/*.mat` |

## Group-level analysis → manuscript mapping table

Each row is a subfolder of `code/group-level/`, containing the function(s) that run the
statistical analysis and generate the corresponding manuscript figure/table.

| Folder | Manuscript aim / section | Key figure/table number(s) | Main statistical method |
|---|---|---|---|
| `aim1_isf_characterization/` | Characterize whether ISFσ can be reliably detected and whether THC/CBD alters its frequency, amplitude, and topography | Figure 2 (panels A–E); Supplementary Figure S1 | Paired-sample t-tests (bout number/duration/onset); Fourier/Gaussian-fit spectral analysis; cross-correlation |
| `aim1b_filter_validation/` | Ancillary validation: confirming the FIR bandpass filter's zero-padding negates filter-edge artefacts on the ISFσ phase estimate | Supplementary Figures S3–S4 | Hilbert-transform phase/amplitude comparison between zero-padded ("control") and unpadded ("test") signals; descriptive circular statistics (`circ_mean`, `circ_std`) |
| `aim2_nrem_sleep_stability/` | Whether NREM sleep-stability (awakening vs. continued-sleep arousals) is modulated by ISFσ phase angle/amplitude, and whether THC/CBD alters this | Figure 3; Supplementary Figures S5–S6; Supplementary Tables S1–S2 | Multiple linear mixed-effects models at each discrete time-sample with permutation-based multiple-comparison correction (`prearousal_permutation_models_optimized.m`); circular statistics (Hermans-Rasson test, Dip test, phase-density peak estimation) |
| `aim3_rem_transitions/` | Whether cannabinoid-related delays in REM sleep onset are explained by alterations in ISFσ amplitude/half-period immediately preceding REM transitions | Figure 4 | Linear mixed-effects models (`fitlme`) of ISFσ half-period/amplitude by fluctuation type × treatment condition and by REM latency; permutation-based correction for the pre-REM sigma-power time series (`transrem_permutation_models_optimized.m`) |
| `sleep_macroarchitecture/` | Sleep macroarchitecture and arousal-outcome comparisons between treatment conditions (previously named "descriptives") | Table 1; Supplementary Table (awakening-probability model); Supplementary Figure (arousal descriptives) | Linear mixed-effects models (`fitlme`) per PSG variable; generalized linear mixed-effects (binomial) models (`fitglme`) for awakening probability |
| `supplementary/` | Illustrative-only: explains the ISFσ phase-angle/Hilbert-transform methodology visually (not itself a manuscript result) | Supplementary Figure S5 (phase-angle distribution illustration) | n/a — descriptive/illustrative plot only |

## Notes and known limitations (carried over from the refactor)

- `group-level/supplementary/plot_isf_phase_distribution_illustration.m` loads a file
  `analysis_2a.mat` whose expected location relative to the working directory was not confirmed
  during the refactor — see the `% NOTE:` comment in that file for details.
- `group-level/aim1b_filter_validation/analyse_filter_edge_artefact.m` contains one pre-existing,
  unrelated syntax quirk (a stray `N` token before a cell-divider, line ~248) that was
  intentionally left untouched during the refactor (out of scope; no logic change made).
- `subject-level/spectral-features/css_getsigmapowerusingwavelet.m` and
  `css_getspecpowerusingwavelet.m` are near-duplicates (hardcoded vs. parametrised frequency
  bands) that were flagged but not merged, since their divergence may be intentional.
