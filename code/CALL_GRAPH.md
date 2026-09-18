# `code/` — Full call-graph and dead-code audit

Generated 2026-09-18 as part of `SRP_REFACTOR_PLAN.md` Phase 1. This document is a snapshot of
static analysis, not a maintained-by-hand reference — re-run the audit method described below
after any future rename/move/delete to keep it accurate.

## Method

A Python script parsed every `function <signature>` declaration (handling MATLAB's `...`
line-continuation syntax, the `[a,b] = name(...)`, `a = name(...)`, and bare-name forms, and
treating files with no `function` keyword as top-level scripts) across all 141 `.m` files in
`code/` (`code/toolboxes/` excluded — third-party, untouched). This produced a map of every
function/script name to its defining file.

A second pass grepped every file in the tree for word-boundary (`\bname\b`) references to every
known name, then performed a breadth-first reachability walk starting from the three real entry
points:
- `main.m` — the top-level orchestration script.
- `css_init.m` — environment/path initialization, called by `main.m`.
- `qc/*.m` — manually-run, interactive quality-control scripts (not called from `main.m`, but
  legitimate entry points a user runs directly by hand).

Any non-archived file not reached from one of these three starting points is a dead-code
candidate. Every dead-code candidate was additionally verified with a repo-wide `grep` (not just
within `code/`) against `README.md`, `documents/manuscript.md`, and
`documents/peer-review-comments.md` to rule out being referenced from documentation/manuscript
text before being archived.

**Caveat:** this is regex-based, not a real MATLAB parser — it cannot detect dynamic dispatch
(`feval`, `str2func`, `eval`). A manual grep for these patterns across `code/` (excluding
`toolboxes/`) found exactly one use, in `utils/plotting-generic/css_standard_colors.m` (`eval`
used internally to look up a local variable by constructed name, not to call an external
function) — so dynamic dispatch is not a blind spot for this audit.

---

## Full call table

Each row: file → the function/script it defines → its type → every other file (outside
`code/archive/`) that references its name. A file whose "Called by" column is
**(none — see note)** is the definitive list of currently-unreachable, non-archived code as of
this audit (see "Dead-code findings" below for what was done about each one).

### Top level (`code/`)

| File | Function | Type | Called by |
|---|---|---|---|
| `css_init.m` | `css_init` | FUNC | `./main.m` |
| `main.m` | `main` | SCRIPT | *(entry point — run directly by the user, not called from other code)* |

### `pipeline/`

| File | Function | Type | Called by |
|---|---|---|---|
| `pipeline/run_processing_section.m` | `run_processing_section` | FUNC | `./main.m` |

### `subject-level/`

| File | Function | Type | Called by |
|---|---|---|---|
| `subject-level/preprocessing/css_preproc.m` | `css_preproc` | FUNC | `./pipeline/run_processing_section.m` |
| `subject-level/segmentation/css_extractarousalbouts.m` | `css_extractarousalbouts` | FUNC | `./pipeline/run_processing_section.m` |
| `subject-level/segmentation/css_extractnrembouts.m` | `css_extractnrembouts` | FUNC | `./pipeline/run_processing_section.m` |
| `subject-level/segmentation/css_extractprerembouts.m` | `css_extractprerembouts` | FUNC | `./pipeline/run_processing_section.m` |
| `subject-level/spectral-features/css_extractspindles.m` | `css_extractspindles` | FUNC | `./pipeline/run_processing_section.m` |
| `subject-level/spectral-features/css_getsigmapowerusingwavelet.m` | `css_getsigmapowerusingwavelet` | FUNC | `./pipeline/run_processing_section.m` |
| `subject-level/spectral-features/css_getspecpowerusingwavelet.m` | `css_getspecpowerusingwavelet` | FUNC | `./pipeline/run_processing_section.m` |

### `firstlevel-outcomes/`

| File | Function | Type | Called by |
|---|---|---|---|
| `firstlevel-outcomes/css_createfstlvloutput.m` | `css_createfstlvloutput` | FUNC | `./firstlevel-outcomes/css_crosscorr.m`; `./firstlevel-outcomes/css_infraslowfluctpowerspect.m`; `./group-level/aim1_isf_characterization/analyse_sigma_spindle_similarity.m` |
| `firstlevel-outcomes/css_crosscorr.m` | `css_crosscorr` | FUNC | `./pipeline/run_processing_section.m` |
| `firstlevel-outcomes/css_extractfeatures.m` | `css_extractfeatures` | FUNC | `./firstlevel-outcomes/css_infraslowfluctpowerspect.m`; `./utils/isf-fitting/interp_isfspectrum.m` |
| `firstlevel-outcomes/css_infraslowfluctpowerspect.m` | `css_infraslowfluctpowerspect` | FUNC | `./pipeline/run_processing_section.m` |

### `group-level/`

| File | Function | Type | Called by |
|---|---|---|---|
| `group-level/aim1_isf_characterization/analyse_bout_selection.m` | `analyse_bout_selection` | FUNC | `./main.m` |
| `group-level/aim1_isf_characterization/analyse_isf_topography.m` | `analyse_isf_topography` | FUNC | `./main.m` |
| `group-level/aim1_isf_characterization/analyse_sigma_spindle_similarity.m` | `analyse_sigma_spindle_similarity` | FUNC | `./main.m` |
| `group-level/aim1_isf_characterization/plot_isf_topography.m` | `plot_isf_topography` | FUNC | `./group-level/aim1_isf_characterization/analyse_isf_topography.m` |
| `group-level/aim1b_filter_validation/analyse_filter_edge_artefact.m` | `analyse_filter_edge_artefact` | FUNC | `./main.m` |
| `group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` | `analyse_arousal_isf_phase` | FUNC | `./main.m` |
| `group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_avsigmatrace.m` | `plot_fig2_avsigmatrace` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` |
| `group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_eegtrace.m` | `plot_fig2_eegtrace` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` |
| `group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_instamp.m` | `plot_fig2_instamp` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` |
| `group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_panellabels.m` | `plot_fig2_panellabels` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` |
| `group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_phaseangle.m` | `plot_fig2_phaseangle` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` |
| `group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_probarousal.m` | `plot_fig2_probarousal` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` |
| `group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_sigmatrace.m` | `plot_fig2_sigmatrace` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` |
| `group-level/aim2_nrem_sleep_stability/prearousal_permutation_models_optimized.m` | `prearousal_permutation_models_optimized` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` |
| `group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m` | `analyse_rem_transition_dynamics` | FUNC | `./main.m` |
| `group-level/aim3_rem_transitions/plotting/plot_fig3_ampdur.m` | `plot_fig3_ampdur` | FUNC | `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m` |
| `group-level/aim3_rem_transitions/plotting/plot_fig3_avsigmatrace.m` | `plot_fig3_avsigmatrace` | FUNC | `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m` |
| `group-level/aim3_rem_transitions/plotting/plot_fig3_mdl.m` | `plot_fig3_mdl` | FUNC | `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m` |
| `group-level/aim3_rem_transitions/plotting/plot_fig3_mdlhists.m` | `plot_fig3_mdlhists` | FUNC | `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m` |
| `group-level/aim3_rem_transitions/plotting/plot_fig3_panellabels.m` | `plot_fig3_panellabels` | FUNC | `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m` |
| `group-level/aim3_rem_transitions/plotting/plot_fig3_sigmatrace.m` | `plot_fig3_sigmatrace` | FUNC | `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m` |
| `group-level/aim3_rem_transitions/transrem_permutation_models_optimized.m` | `transrem_permutation_models_optimized` | FUNC | `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m` |
| `group-level/sleep_macroarchitecture/analyse_sleep_macroarchitecture_and_arousal_outcomes.m` | `analyse_sleep_macroarchitecture_and_arousal_outcomes` | FUNC | `./main.m` |
| `group-level/supplementary/plot_isf_phase_distribution_illustration.m` | `plot_isf_phase_distribution_illustration` | FUNC | `./main.m` |

### `qc/`

| File | Function | Type | Called by |
|---|---|---|---|
| `qc/check_isfspectra.m` | `check_isfspectra` | SCRIPT | *(manual QC entry point — run directly by the user)* |
| `qc/check_ultralowfluctuations.m` | `check_ultralowfluctuations` | SCRIPT | *(manual QC entry point — run directly by the user)* |
| `qc/css_inspectspindles.m` | `css_inspectspindles` | FUNC | `./pipeline/run_processing_section.m` |
| `qc/plothrvhypnogram.m` | `plothrvhypnogram` | SCRIPT | *(manual QC entry point — run directly by the user)* |
| `qc/plotisfspect.m` | `plotisfspect` | FUNC | *(manual QC entry point — call directly with an `ISF` struct)* |

### `utils/`

| File | Function | Type | Called by |
|---|---|---|---|
| `utils/circular-stats/averageSigmaPerEvent.m` | `averageSigmaPerEvent` | FUNC | `./group-level/aim2_nrem_sleep_stability/prearousal_permutation_models_optimized.m` |
| `utils/circular-stats/correctPhaseByEmpiricalCDF.m` | `correctPhaseByEmpiricalCDF` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` |
| `utils/circular-stats/getarousalphaseangle.m` | `getarousalphaseangle` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` |
| `utils/circular-stats/getpermpvalue.m` | `getpermpvalue` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m`; `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m` |
| `utils/circular-stats/parsepvalue.m` | `parsepvalue` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` |
| `utils/circular-stats/permuteEventLabels.m` | `permuteEventLabels` | FUNC | `./group-level/aim2_nrem_sleep_stability/prearousal_permutation_models_optimized.m`; `./group-level/aim3_rem_transitions/transrem_permutation_models_optimized.m` |
| `utils/circular-stats/withinChanCircMean.m` | `withinChanCircMean` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m`; `./group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_phaseangle.m` |
| `utils/eeg-io/css_eeglab2hypnogram.m` | `css_eeglab2hypnogram` | FUNC | `./group-level/aim1_isf_characterization/analyse_isf_topography.m`; `./group-level/aim1_isf_characterization/plot_isf_topography.m`; `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m`; `./qc/check_ultralowfluctuations.m`; `./qc/plothrvhypnogram.m`; `./subject-level/segmentation/css_extractarousalbouts.m`; `./subject-level/segmentation/css_extractnrembouts.m`; `./subject-level/segmentation/css_extractprerembouts.m`; `./subject-level/spectral-features/css_extractspindles.m`; `./utils/misc/inspect_n2aros.m` |
| `utils/eeg-io/events2timeseries.m` | `events2timeseries` | FUNC | `./subject-level/spectral-features/css_extractspindles.m` |
| `utils/eeg-io/getarousalbouts.m` | `getarousalbouts` | FUNC | `./subject-level/segmentation/css_extractarousalbouts.m`; `./utils/misc/inspect_n2aros.m` |
| `utils/eeg-io/getnrembouts.m` | `getnrembouts` | FUNC | `./group-level/aim1_isf_characterization/analyse_isf_topography.m`; `./group-level/aim1_isf_characterization/plot_isf_topography.m`; `./qc/check_ultralowfluctuations.m`; `./qc/plothrvhypnogram.m`; `./subject-level/segmentation/css_extractnrembouts.m` |
| `utils/eeg-io/storeoriglatency.m` | `storeoriglatency` | FUNC | `./subject-level/segmentation/css_extractarousalbouts.m`; `./subject-level/segmentation/css_extractnrembouts.m`; `./subject-level/segmentation/css_extractprerembouts.m` |
| `utils/isf-fitting/fitisfspect.m` | `fitisfspect` | FUNC | `./firstlevel-outcomes/css_infraslowfluctpowerspect.m`; `./utils/isf-fitting/interp_isfspectrum.m` |
| `utils/isf-fitting/gaussianfit.m` | `gaussianfit` | FUNC | `./firstlevel-outcomes/css_infraslowfluctpowerspect.m`; `./utils/isf-fitting/fitisfspect.m`; `./utils/isf-fitting/interp_isfspectrum.m` |
| `utils/isf-fitting/interp_isfspectrum.m` | `interp_isfspectrum` | FUNC | `./group-level/aim1_isf_characterization/analyse_isf_topography.m` |
| `utils/isf-fitting/predictISF.m` | `predictISF` | FUNC | `./utils/circular-stats/getarousalphaseangle.m` |
| `utils/isf-fitting/predict_infraslow_eeg.m` | `predict_infraslow_eeg` | FUNC | `./utils/isf-fitting/predictISF.m` |
| `utils/isf-fitting/prerempeakstroughs.m` | `prerempeakstroughs` | FUNC | `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m` |
| `utils/isf-fitting/zerocrosspeakfind.m` | `zerocrosspeakfind` | FUNC | `./utils/circular-stats/getarousalphaseangle.m` |
| `utils/misc/calcarovars.m` | `calcarovars` | FUNC | **(none — see note)** |
| `utils/misc/calcpsgvars.m` | `calcpsgvars` | FUNC | **(none — see note)** |
| `utils/misc/getuuid.m` | `getuuid` | FUNC | `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m`; `./utils/circular-stats/getarousalphaseangle.m`; `./utils/isf-fitting/prerempeakstroughs.m` |
| `utils/misc/inspect_n2aros.m` | `inspect_n2aros` | SCRIPT | **(none — see note)** |
| `utils/misc/loadarousalcsvs.m` | `loadarousalcsvs` | FUNC | **(none — see note)** |
| `utils/misc/loadarousalepochs.m` | `loadarousalepochs` | FUNC | **(none — see note)** |
| `utils/misc/withinConditionNorm.m` | `withinConditionNorm` | FUNC | **(none — see note)** |
| `utils/misc/withinSubMean.m` | `withinSubMean` | FUNC | `./group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_avsigmatrace.m`; `./group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_instamp.m`; `./group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_probarousal.m` |
| `utils/plotting-generic/css_standard_colors.m` | `css_standard_colors` | FUNC | `./group-level/aim1_isf_characterization/plot_isf_topography.m`; `./group-level/aim1b_filter_validation/analyse_filter_edge_artefact.m`; `./group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_avsigmatrace.m`; `./group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_eegtrace.m`; `./group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_instamp.m`; `./group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_phaseangle.m`; `./group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_probarousal.m`; `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m`; `./group-level/aim3_rem_transitions/plotting/plot_fig3_ampdur.m`; `./group-level/aim3_rem_transitions/plotting/plot_fig3_avsigmatrace.m`; `./group-level/aim3_rem_transitions/plotting/plot_fig3_mdl.m`; `./group-level/aim3_rem_transitions/plotting/plot_fig3_mdlhists.m`; `./group-level/aim3_rem_transitions/plotting/plot_fig3_sigmatrace.m`; `./group-level/sleep_macroarchitecture/analyse_sleep_macroarchitecture_and_arousal_outcomes.m`; `./group-level/supplementary/plot_isf_phase_distribution_illustration.m` |
| `utils/plotting-generic/errorpatch.m` | `errorpatch` | FUNC | `./group-level/aim2_nrem_sleep_stability/plotting/plot_fig2_avsigmatrace.m`; `./group-level/aim3_rem_transitions/plotting/plot_fig3_avsigmatrace.m`; `./group-level/sleep_macroarchitecture/analyse_sleep_macroarchitecture_and_arousal_outcomes.m` |
| `utils/plotting-generic/linepatch.m` | `linepatch` | FUNC | `./group-level/aim1b_filter_validation/analyse_filter_edge_artefact.m` |
| `utils/plotting-generic/plotFiltParams.m` | `plotFiltParams` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` |
| `utils/plotting-generic/plotarousalboutselection.m` | `plotarousalboutselection` | FUNC | `./subject-level/segmentation/css_extractarousalbouts.m` |
| `utils/plotting-generic/plotnremboutselection.m` | `plotnremboutselection` | FUNC | `./subject-level/segmentation/css_extractnrembouts.m`; `./subject-level/segmentation/css_extractprerembouts.m` |
| `utils/signal-append/applyappend.m` | `applyappend` | FUNC | `./utils/signal-append/isffilterbout.m` |
| `utils/signal-append/ar_pred.m` | `ar_pred` | FUNC | `./utils/signal-append/signalappend.m` |
| `utils/signal-append/doubleflip.m` | `doubleflip` | FUNC | `./utils/signal-append/signalappend.m` |
| `utils/signal-append/executeappending.m` | `executeappending` | FUNC | `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m` |
| `utils/signal-append/isffilterbout.m` | `isffilterbout` | FUNC | `./utils/circular-stats/getarousalphaseangle.m` |
| `utils/signal-append/rmappend.m` | `rmappend` | FUNC | `./utils/isf-fitting/predictISF.m` |
| `utils/signal-append/signalappend.m` | `signalappend` | FUNC | `./group-level/aim1b_filter_validation/analyse_filter_edge_artefact.m`; `./group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m`; `./utils/signal-append/isffilterbout.m` |
| `utils/signal-processing/calcinsthr.m` | `calcinsthr` | FUNC | `./subject-level/preprocessing/css_preproc.m` |
| `utils/signal-processing/ecg2hr_pantompkin.m` | `ecg2hr_pantompkin` | FUNC | `./utils/signal-processing/calcinsthr.m` |
| `utils/signal-processing/ft_wavelettransform.m` | `ft_wavelettransform` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m`; `./subject-level/spectral-features/css_getsigmapowerusingwavelet.m`; `./subject-level/spectral-features/css_getspecpowerusingwavelet.m` |
| `utils/signal-processing/waveletsettings.m` | `waveletsettings` | FUNC | `./group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m`; `./subject-level/spectral-features/css_getsigmapowerusingwavelet.m`; `./subject-level/spectral-features/css_getspecpowerusingwavelet.m` |
| `utils/spindle-detection/convspindles.m` | `convspindles` | FUNC | **(none — see note)** |
| `utils/spindle-detection/f_SpDetection_Ferrarelli.m` | `f_SpDetection_Ferrarelli` | FUNC | `./subject-level/spectral-features/css_extractspindles.m` |
| `utils/spindle-detection/f_SpDetection_Humans.m` | `f_SpDetection_Humans` | FUNC | `./subject-level/spectral-features/css_extractspindles.m` |
| `utils/spindle-detection/f_SpDetection_Wamsley.m` | `f_SpDetection_Wamsley` | FUNC | `./subject-level/spectral-features/css_extractspindles.m` |

---

## Dead-code findings and disposition

### Archived in this audit (18 functions, zero call sites anywhere in the repo)

Confirmed via both the reachability walk above and a manual repo-wide grep (including
`README.md`, `documents/manuscript.md`, `documents/peer-review-comments.md`) — none were
referenced by name anywhere outside their own file. All moved to `code/archive/scratch/` via
`git mv` (history preserved), each with an `% ARCHIVED 2026-09-18: ...` header comment.

| Archived function | Reason |
|---|---|
| `withinChanCircMedian` | No callers anywhere. |
| `eeglab2arousals` | No callers anywhere. |
| `revisemistakesineventtable` | No callers anywhere. |
| `plothypno` | Superseded by EEG_Processor's `plotHypnogram`, which every active caller in `code/` actually uses. |
| `csapsGCV` | No callers anywhere. |
| `eeg2freqbandpower` | Only used `morletgabortransform`, itself unused — self-contained dead pair. |
| `gethilbert` | No callers anywhere. |
| `magnituderesponse` | No callers anywhere. |
| `morletgabortransform` | Only used by the also-dead `eeg2freqbandpower`. |
| `nonparzscore` | No callers anywhere. |
| `normdistance` | No callers anywhere. |
| `phaseresponse` | No callers anywhere. |
| `shorttimefft` | No callers anywhere. |
| `smooth1q` | No callers anywhere. |
| `zscoreacrosschannels` | No callers anywhere. |
| `f_IFO_Parameters_140324` | Superseded by `css_infraslowfluctpowerspect.m`/`fitisfspect.m`, the ISF-fitting pipeline actually wired into `main.m`. Only referenced internally by the also-dead `f_IFO_WithSpindles`. |
| `f_IFO_WithSpindles` | Same supersession as above. |
| `f_MGT` | Only used by the two `f_IFO_*` functions above — fully self-contained dead cluster of 3. |

### Already flagged and kept (7 functions, unreachable but intentionally retained)

These were already flagged as having no call sites during the original `REFACTOR_PLAN.md` Phase
3 and kept with a `% NOTE: no current call sites found in code/ as of 2026-09-18 (Phase 3...`
comment — most likely because they generate phenotype/data files consumed outside this repo
(e.g. via SPSS or manual analysis), so archiving them risks losing undocumented data-provenance
logic. **Not touched in this audit** — re-confirmed their call-site status is unchanged:

- `utils/misc/calcarovars.m`
- `utils/misc/calcpsgvars.m`
- `utils/misc/inspect_n2aros.m` (manual/interactive QC script)
- `utils/misc/loadarousalcsvs.m`
- `utils/misc/loadarousalepochs.m`
- `utils/misc/withinConditionNorm.m`
- `utils/spindle-detection/convspindles.m`

### Known duplicate — fixed in Phase 2

`utils/eeg-io/eeglab2hypnogram.m` and `utils/eeg-io/css_eeglab2hypnogram.m` contained
byte-for-byte identical logic under two different names. Every active caller used
`css_eeglab2hypnogram` **except** `subject-level/spectral-features/css_extractspindles.m`, which
called the unprefixed `eeglab2hypnogram` — a latent bug (relied on two copies of the same logic
staying in sync). Fixed in `SRP_REFACTOR_PLAN.md` Phase 2: the stray call site was repointed to
`css_eeglab2hypnogram`, and the now-fully-unreferenced `eeglab2hypnogram.m` was archived to
`code/archive/scratch/eeglab2hypnogram.m`.

---

## Bugs found and fixed as a side-effect of this audit

Reading every file in `code/` line-by-line to build the call-graph surfaced two unambiguous
syntax/logic errors, unrelated to the naming/SRP/dead-code task but fixed immediately since they
would break execution and are zero-risk, zero-ambiguity fixes (no numeric/statistical change):

1. **`main.m` line 147:** `cfg.do_parallel = do_parallel;n` — a stray trailing `n` character
   (likely a fat-fingered keystroke during a manual edit) that would cause a MATLAB parse error
   if this section of `main.m` were ever run. Fixed to `cfg.do_parallel = do_parallel;`.
2. **`subject-level/segmentation/css_extractarousalbouts.m` line 69:** a bare `keyboard`
   breakpoint left in the main execution path, immediately after the `getarousalbouts` call and
   before the bout-extraction loop — would halt any real, non-interactive run of this function.
   Removed.

---

## Reachability summary (post Phase 2 state)

- Total non-archive `.m` files in `code/` (excl. `toolboxes/`): **91** (92 after Phase 1, minus
  the `eeglab2hypnogram.m` duplicate archived in Phase 2)
- Reached from `main.m` / `css_init.m` / `qc/*.m`: **84**
- Unreached (all 7 are the intentionally-kept, already-NOTE-flagged files listed above): **7**

This matches expectations exactly — after archiving the 18 Phase-1 dead functions and the Phase-2
duplicate, the only remaining "unreached" files are the ones already known and intentionally kept
from the original `REFACTOR_PLAN.md` Phase 3.

