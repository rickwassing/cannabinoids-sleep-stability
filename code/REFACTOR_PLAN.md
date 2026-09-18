# MATLAB Codebase Refactoring Plan

**Status:** APPROVED — not yet executed. Update the checklist at the bottom as each phase completes.

**Scope:** `code/` directory only (excludes `code/toolboxes/`, which is third-party and untouched).
**Goal:** Reorganise into subject-level processing → first-level outcomes → group-level analysis
(mapped to manuscript aims), without changing any statistical method, model specification,
threshold, filter parameter, or numeric result. Preserve `rickwassing` coding style (spacing,
blank-line block separation, purpose-comments). Use `git mv` for every relocation so history
is preserved.

**Repo context:** BIDS-style layout. `documents/manuscript.md` is the source of truth for how
analyses map to manuscript aims/figures/tables. `code/main.m` is the top-level orchestration
script; `code/section.m` is the generic subject-level job runner it calls into
(`Proc`/`Files`/`cfg` pattern, big switch-case dispatch to `code/processing/*` functions).
Derivative folders (gitignored, on disk only): `derivatives/EEG-preproc` → `EEG-processed` →
`EEG-segmented` → `EEG-output-fstlvl` (first-level outcomes) → consumed by `code/analysis/*`
together with externally-computed PALM cluster-stats in `/group-level/<name>/` (root-level,
data only — NOT the same as the new `code/group-level/` code folder proposed below).

**Approved answers to open questions (from planning phase):**
1. Keep `css_descriptives.m` separate from `css_analyse_2b.m`; rename each to reflect the
   specific manuscript aim/result it addresses (not generic "descriptives").
2. Archive (don't delete) confirmed-superseded files: `css_analyse_2.m`, `css_plot_isffeatures.m`,
   `code/check.m`, `code/tmp/analyze_signal_frequencies.m`.
3. Update `main.m` to call every group-level analysis that actually produced a manuscript
   result (currently missing calls to `css_analyse_2b` and `css_analyse_3`), so the pipeline
   script reflects the full, real analysis chain.
4. Consolidate `code/archive/`, `code/analysis/archive/`, `code/tmp/` into one
   `code/archive/{pipeline,analysis,scratch}/`.
5. Scripts that depend on leftover workspace variables from prior cells
   (`css_analyse_1b.m`, `css_analyse_2.m`, `css_analyse_3.m`, `css_descriptives.m`) are to be
   converted into self-contained functions with explicit inputs — same logic/statistics/output,
   just independently callable. `css_analyse_1b.m`'s undefined `Files` bug will be fixed as part
   of this conversion (bug-fix, not a methodology change — will be called out explicitly in the
   commit/summary).
6. Folder naming: keep `code/group-level/` (accept the semantic overlap with the root-level
   `/group-level/` PALM-output data folder — they sit at different paths and serve different
   purposes: one is code, one is data).

**Tooling notes (general, apply throughout):**
- No MATLAB CLI is available in this environment (`which matlab` fails); `/Applications` has
  MATLAB_R2019b.app and MATLAB_R2021b.app (GUI only, not scriptable headlessly here) — so
  **numeric/behavioural validation of converted functions cannot be run in this environment**.
  Where a validation step calls for "run and confirm output matches", this instead means:
  static review (read old vs. new file side by side, confirm every line of logic is preserved),
  `grep`-based dependency verification, and MATLAB syntax sanity-checking by eye (matching
  `end` keywords, function signatures, etc.). If you have MATLAB available locally, running
  `matlab -batch "addpath(genpath('code')); <function>()"` after each phase and diffing any
  generated `.csv`/`.mat`/figure output against a pre-refactor copy is recommended and is noted
  per-phase below as an optional external validation step you can run yourself.
- No specialised MCP server exists for MATLAB code search/refactoring; the `search_codebase`
  and `read_files` tools (regex + line-range reads) are what's used throughout, which has been
  sufficient. If you want deeper call-graph validation, MATLAB's own `checkcode`/`mlint` or the
  "Find Files" / "Dependency Report" in the MATLAB IDE would be the standard tool — you could
  run `matlab.codetools.requiredFilesAndProducts('code/pipeline/main.m')` locally to get a
  definitive dependency list before/after each phase and diff the two lists.
- No web/documentation lookups are needed for this task (pure internal reorganisation) — Chrome
  web2md / Google are not required unless you want to double check EEGLAB/FieldTrip function
  signatures (`pop_firws`, `pop_select`, `ft_freqanalysis`, etc.), which are not being modified,
  only relocated call sites.

---

## Phase 0 — Setup

**Read first:** `@code/REFACTOR_PLAN.md` (this file, to confirm it matches what was approved).

- [x] Write this plan file (done).
- [x] Create the new directory skeleton (empty folders, so `git mv` targets exist):
  `code/pipeline/`, `code/subject-level/{preprocessing,arousal-detection,spectral-features,segmentation,inspection}/`,
  `code/firstlevel-outcomes/`, `code/group-level/{aim1_isf_characterization,aim1b_filter_validation,aim2_nrem_sleep_stability/plotting,aim3_rem_transitions/plotting,sleep_macroarchitecture,supplementary}/`,
  `code/utils/{eeg-io,signal-processing,spindle-detection,circular-stats,isf-fitting,arousal-bouts,plotting-generic,signal-append,misc}/`,
  `code/archive/{pipeline,analysis,scratch}/`. **Done 2026-09-18.**
- [x] Confirm `code/qc/` and `code/toolboxes/` remain untouched (no action needed). **Confirmed unchanged.**

**Skills/tools:** `run_commands` (`mkdir -p`) only. No MCP server needed.

---

## Phase 1 — Archive consolidation (lowest risk; no logic touched)

**Read first:**
`@code/archive/main_mbpro.m`, `@code/archive/main_supc.m`, `@code/archive/css_figure1.m`,
`@code/archive/spectral_determinants.m`, `@code/archive/tmp_checkN1.m`,
`@code/archive/tmp_inspect_tstat_spectra.m`, `@code/archive/tmp_iso_freqlimits.m`,
`@code/archive/tmp_recalcgausfit.m`, `@code/archive/tmp_recalcismphase.m`,
`@code/analysis/archive/css_arousalawakeningproportion.m`,
`@code/analysis/archive/css_isoampfreqdifference.m`, `@code/analysis/archive/css_isophaseangle.m`,
`@code/analysis/archive/css_prearousalspectrum.m`, `@code/analysis/archive/css_prearousalspindles.m`,
`@code/tmp/analyze_signal_frequencies.m`, `@code/check.m`,
`@code/analysis/css_analyse_2.m`, `@code/analysis/css_plot_isffeatures.m`.
(All already read once during planning — re-read only if resuming in a new session.)

**Actions:**
- [x] `git mv` main_mbpro.m/.asv, main_supc.m/.asv, css_figure1.m, spectral_determinants.m → `code/archive/pipeline/`. (`.asv` files confirmed untracked/gitignored, moved with plain `mv`.)
- [x] `git mv` the 5 `tmp_*.m` files → `code/archive/scratch/`.
- [x] `git mv code/analysis/archive/*.m` → `code/archive/analysis/` (incl. untracked `.asv`).
- [x] `git mv code/analysis/css_analyse_2.m` → `code/archive/analysis/css_analyse_2.m` — header comment added.
- [x] `git mv code/analysis/css_plot_isffeatures.m` (+ untracked `.asv`) → `code/archive/analysis/css_plot_isffeatures.m` — header comment added.
- [x] `git mv code/tmp/analyze_signal_frequencies.m` → `code/archive/scratch/analyze_signal_frequencies.m`
- [x] `git mv code/check.m` → `code/archive/scratch/check.m`
- [x] Removed now-empty `code/analysis/archive/`, `code/tmp/` directories.
- [x] `grep -rn` (word-boundary) for the old archived filenames across `code/` — **zero dangling
      functional references confirmed.** One pre-existing, harmless typo found and left for
      Phase 4 (a comment `%function css_analyse_2` inside the *still-active* `css_analyse_2b.m`,
      mislabelling itself — not a reference to the archived file; will be corrected when
      `css_analyse_2b.m` is renamed in Phase 4f).

**Phase 1 completed 2026-09-18.**

**Skills/tools:** `run_commands` (`git mv`, `grep`), `editor` (header comments). No MCP server needed.

---

## Phase 2 — Subject-level processing + first-level outcomes reorganisation

**Read first:**
`@code/processing/css_preproc.m`, `@code/processing/css_init.m`,
`@code/processing/css_detectarousals.m`, `@code/processing/css_extractarousalbouts.m`,
`@code/processing/css_extractspindles.m`, `@code/processing/css_getsigmapowerusingwavelet.m`,
`@code/processing/css_getspecpowerusingwavelet.m`, `@code/processing/css_extractnrembouts.m`,
`@code/processing/css_extractprerembouts.m`, `@code/processing/css_inspectspindles.m`,
`@code/processing/css_infraslowfluctpowerspect.m`, `@code/processing/css_crosscorr.m`,
`@code/analysis/css_extractfeatures.m`, `@code/analysis/css_createfstlvloutput.m`,
`@code/section.m`, `@code/main.m`.
(All already read during planning; re-read if resuming fresh.)

**File moves (as originally planned):**

| Current path | Planned new path | Notes |
|---|---|---|
| `code/processing/css_preproc.m` | `code/subject-level/preprocessing/css_preproc.m` | pure move |
| `code/processing/css_init.m` | `code/subject-level/preprocessing/css_init.m` | pure move |
| `code/processing/css_detectarousals.m` | `code/subject-level/arousal-detection/css_detectarousals.m` | pure move |
| `code/processing/css_extractarousalbouts.m` | `code/subject-level/arousal-detection/css_extractarousalbouts.m` | pure move |
| `code/processing/css_extractspindles.m` | `code/subject-level/spectral-features/css_extractspindles.m` | pure move |
| `code/processing/css_getsigmapowerusingwavelet.m` | `code/subject-level/spectral-features/css_getsigmapowerusingwavelet.m` | pure move |
| `code/processing/css_getspecpowerusingwavelet.m` | `code/subject-level/spectral-features/css_getspecpowerusingwavelet.m` | pure move; flag near-duplication with sigma version in code comment, do not merge |
| `code/processing/css_extractnrembouts.m` | `code/subject-level/segmentation/css_extractnrembouts.m` | pure move |
| `code/processing/css_extractprerembouts.m` | `code/subject-level/segmentation/css_extractprerembouts.m` | pure move |
| `code/processing/css_inspectspindles.m` | `code/subject-level/inspection/css_inspectspindles.m` | pure move |
| `code/processing/html/*` | `code/subject-level/inspection/html/*` | generated output viewer, pure move |
| `code/processing/css_infraslowfluctpowerspect.m` | `code/firstlevel-outcomes/css_infraslowfluctpowerspect.m` | pure move — see note below |
| `code/processing/css_crosscorr.m` | `code/firstlevel-outcomes/css_crosscorr.m` | pure move |
| `code/analysis/css_extractfeatures.m` | `code/firstlevel-outcomes/css_extractfeatures.m` | pure move |
| `code/analysis/css_createfstlvloutput.m` | `code/firstlevel-outcomes/css_createfstlvloutput.m` | pure move |
| `code/section.m` | `code/pipeline/run_processing_section.m` | pure move + rename; update internal `function [errors] = section(...)` signature name to `run_processing_section` to match filename (MATLAB requires function name == filename) |
| `code/main.m` | `code/pipeline/main.m` | pure move; content updated in Phase 5 |

**⚠️ CORRECTION (added 2026-09-18, after manual edits by repo owner post-Phase-2):**
The table above is the *original plan*. What actually landed on disk (confirmed via
`git show 5bf8aa8 --name-status -M` and a fresh `find code -type f`) diverges from it in seven
places — apparently due to manual edits made outside this planning session. The original Phase 2
completion note below incorrectly asserted "no other discrepancy found"; that assertion is now
corrected. **Actual, current on-disk state (authoritative going forward):**

| File | Plan said | What's actually on disk now |
|---|---|---|
| `code/processing/css_init.m` | → `code/subject-level/preprocessing/css_init.m` | → **`code/css_init.m`** (kept at `code/` root, alongside `main.m`) |
| `code/main.m` | → `code/pipeline/main.m` | **Still at `code/main.m`** (not moved; `code/pipeline/` currently holds only `run_processing_section.m`) |
| `code/processing/css_detectarousals.m` | → `code/subject-level/arousal-detection/...` | **Deleted** (no destination; no remaining references in `code/main.m` or elsewhere) |
| `code/processing/css_extractarousalbouts.m` | → `code/subject-level/arousal-detection/...` | → **`code/subject-level/segmentation/css_extractarousalbouts.m`** (landed in `segmentation/`, not a separate `arousal-detection/` folder) |
| `code/processing/css_inspectspindles.m` | → `code/subject-level/inspection/...` | → **`code/qc/css_inspectspindles.m`** (routed into the pre-existing `qc/` folder instead of a new `inspection/` folder) |
| `code/processing/html/*` | → `code/subject-level/inspection/html/*` | **Deleted** (generated output viewer removed, not relocated) |
| `code/toolboxes/pngquality.py` | not part of Phase 2 (toolboxes out of scope) | **Deleted** |

Net effect: the empty `code/subject-level/arousal-detection/` and `code/subject-level/inspection/`
directories created in Phase 0 are no longer present on disk (nothing was ever moved into them, and
they were not tracked as empty dirs by git), and `code/pipeline/` holds only
`run_processing_section.m`, not `main.m`. All other rows in the table above were executed exactly
as planned. This correction is documentation-only — **no file moves were made or reversed as part
of this correction**; the current on-disk layout (root-level `code/main.m` + `code/css_init.m`,
`code/qc/css_inspectspindles.m`, `code/subject-level/segmentation/css_extractarousalbouts.m`, no
`css_detectarousals.m`) is treated as the accepted, intentional state and is what later phases
(5 and 6) should reference and document, not the original table.

**Note on `css_infraslowfluctpowerspect.m` / `css_crosscorr.m`:** these two functions do subject-level
signal processing AND save the final first-level output (`css_createfstlvloutput`) in one function
body. Per the "avoid premature refactoring" constraint, **do not split the function bodies in this
pass** — moving the whole file to `firstlevel-outcomes/` is sufficient to make the boundary visible
at the folder level. Splitting internals is deferred/optional (see "Deferred / Optional" section).

**Actions:**
- [x] Execute all `git mv` in the table above. **Note:** `code/analysis/css_extractfeatures.m`
      as listed in the table above was a stale path from planning — on disk (confirmed via
      `git log --follow`, present since the initial commit) this file has always lived at
      `code/processing/css_extractfeatures.m`. Moved from its actual location; no other
      discrepancy found.
- [x] In moved `run_processing_section.m`: renamed function declaration line
      `function [errors] = section(Proc, Files, cfg)` → `function [errors] = run_processing_section(Proc, Files, cfg)`.
- [x] `grep -rn "\bsection(" code/` and updated every call site — all 9 were in
      `code/pipeline/main.m` as expected, now `run_processing_section(...)`.
- [x] `grep -rn "'processing/\|processing/css_\|analysis/css_extractfeatures\|analysis/css_createfstlvloutput"` across `code/` and `documents/` — **zero stray path-string references found** (confirmed: file lookups use `dir('derivatives/...')` glob patterns, not `code/` paths).
- [x] Verified `addpath(genpath('code'))` in `main.m`'s init block still present and unchanged — path-based dispatch is folder-agnostic, no code changes needed.

**Phase 2 completed 2026-09-18.** All moves via `git mv` (history preserved). No statistical
method, model specification, threshold, filter parameter, or numeric result was changed — this
was pure relocation plus the one approved function rename (`section` → `run_processing_section`).
`code/processing/` is now drained of all `.m` files (only an empty, pre-existing, untracked
`code/processing/archive/` and `.DS_Store` remain — directory removal deferred to Phase 6 per plan).

**Correction (2026-09-18, later same day):** the "no other discrepancy found" claim above was
inaccurate. Manual edits made after this note was written diverged from the file-moves table in
seven places — see the "⚠️ CORRECTION" block above the table for the authoritative current
on-disk state. Re-verified via `git status --short` (clean) and `find code -maxdepth 2 -type d`
before starting Phase 3.

**Skills/tools:** `run_commands` (`git mv`, `grep`), `editor` (function-name rename in `run_processing_section.m`, call-site updates). No MCP server needed.

---

## Phase 3 — `utils/` (supportfunc) categorisation

**Read first:** all 95 files were listed via `@code/supportfunc/*.m`; the ones whose usage was
directly inspected during planning are listed inline below. If resuming this phase fresh, re-run:
`search_codebase` for each filename to reconfirm call sites before moving, since new files may have
been added/changed since planning.

**Category A — move to `code/utils/eeg-io/`** (EEGLAB struct / hypnogram / event I/O helpers):
`css_eeglab2hypnogram.m`, `eeglab2hypnogram.m`, `eeglab2arousals.m`, `events2timeseries.m`,
`storeoriglatency.m`, `revisemistakesineventtable.m`, `getarousalbouts.m`, `getnrembouts.m`.

**Category B — move to `code/utils/signal-processing/`** (generic spectral/filtering utilities):
`ft_wavelettransform.m`, `waveletsettings.m`, `gethilbert.m`, `shorttimefft.m`, `csapsGCV.m`,
`smooth1q.m`, `morletgabortransform.m`, `eeg2freqbandpower.m`, `magnituderesponse.m`,
`phaseresponse.m`, `ecg2hr_pantompkin.m`, `calcinsthr.m`, `zscoreacrosschannels.m`,
`nonparzscore.m`, `normdistance.m`.

**Category C — move to `code/utils/signal-append/`** (edge-artefact padding helpers used by ISF filtering):
`doubleflip.m`, `ar_pred.m`, `signalappend.m`, `applyappend.m`, `executeappending.m`,
`rmappend.m`, `isffilterbout.m`.

**Category D — move to `code/utils/isf-fitting/`** (ISF Gaussian-fit / phase-extraction helpers):
`fitisfspect.m`, `gaussianfit.m`, `zerocrosspeakfind.m`, `prerempeakstroughs.m`,
`interp_isfspectrum.m`, `predictISF.m`, `predict_infraslow_eeg.m`.
- **Flag for confirmation:** `infraslowmodpowerspect.m`, `infraslowmodpowerspect_old.m` are only
  referenced from `code/archive/scratch/tmp_recalcgausfit.m` (already archived in Phase 1) —
  candidate to move to `code/archive/scratch/` alongside it instead of `utils/`. Propose this;
  confirm before moving.
- **Flag for confirmation:** `extractismphase.m` is only referenced from
  `code/archive/scratch/tmp_recalcismphase.m` — same treatment proposed.

**Category E — move to `code/utils/spindle-detection/`** (third-party-derived spindle algorithms):
`f_IFO_Parameters_140324.m`, `f_IFO_WithSpindles.m`, `f_MGT.m`, `f_SpDetection_Ferrarelli.m`,
`f_SpDetection_Humans.m`, `f_SpDetection_Wamsley.m`.
- **Flag for confirmation:** `convspindles.m`, `tmp_spindle_detection.m` — no call sites found
  anywhere in `code/`. Propose moving `tmp_spindle_detection.m` to `code/archive/scratch/` (name
  itself indicates scratch) and `convspindles.m` to `code/utils/spindle-detection/` with a
  `% NOTE: no current call sites found in code/ as of <date>` comment rather than archiving,
  since it's not obviously scratch-named. Confirm both before moving.

**Category F — move to `code/utils/circular-stats/`:**
`withinChanCircMean.m`, `withinChanCircMedian.m`, `correctPhaseByEmpiricalCDF.m`,
`getarousalphaseangle.m`, `getpermpvalue.m`, `parsepvalue.m`, `permuteEventLabels.m`,
`averageSigmaPerEvent.m`.

**Category G — move to `code/utils/plotting-generic/`** (non analysis-specific plot helpers):
`errorpatch.m`, `linepatch.m`, `css_standard_colors.m`, `plothypno.m`, `plotnremboutselection.m`,
`plotarousalboutselection.m`, `plotFiltParams.m`.
- **Flag for confirmation:** `plothypnoboutselection.m`, `plotarousalersp.m`,
  `plotallsigmachannels.m`, `plot2dspectrum.m`, `plot_psgresults.m`, `plotsubjectlvlismfit.m`,
  `topoplotismfitpeaks.m`, `check_isobw.m` — no call sites found anywhere in current `code/`
  (some only referenced from files already archived in Phase 1, e.g. `plotarousalersp.m` /
  `plot2dspectrum.m` from `archive/spectral_determinants.m`, `plotallsigmachannels.m` from
  `archive/tmp_recalcgausfit.m`). Propose: move the two referenced-only-by-archive files
  (`plotarousalersp.m`, `plot2dspectrum.m`, `plotallsigmachannels.m`) to `code/archive/scratch/`;
  for the remaining unreferenced ones (`plothypnoboutselection.m`, `plot_psgresults.m`,
  `plotsubjectlvlismfit.m`, `topoplotismfitpeaks.m`, `check_isobw.m` — none have a `function`
  keyword at all, i.e. they are bare scripts, likely interactive/manual QC snippets), propose
  moving to `code/archive/scratch/` too. **Confirm before moving — these could be manually
  invoked interactive QC tools you still use.**

**Category H — move to `code/utils/misc/`:**
`getuuid.m`, `withinSubMean.m`.
- **Flag for confirmation:** `withinConditionNorm.m`, `calcarovars.m`, `calcpsgvars.m`,
  `loadarousalcsvs.m`, `loadarousalepochs.m`, `inspect_n2aros.m` — no call sites found in `code/`.
  `calcarovars.m`/`calcpsgvars.m` names strongly suggest they generate the phenotype CSVs consumed
  by `css_descriptives.m` (`phenotype/2024-07-26T171508_psg-variables.csv`,
  `..._arousals.csv`) but that generation may happen outside this repo (e.g., in SPSS or a
  separate script) — **do not archive these two without your confirmation**, since deleting/hiding
  a data-generation function whose output is still used would be a real regression risk. Propose
  keeping in `utils/misc/` with a comment flagging unclear current call site.

**Category I — move to `code/group-level/aim2_nrem_sleep_stability/`** (private helpers of `css_analyse_2b`):
`prearousal_permutation_models_optimized.m` (active), `plot_fig2_avsigmatrace.m`,
`plot_fig2_eegtrace.m`, `plot_fig2_instamp.m`, `plot_fig2_panellabels.m`, `plot_fig2_phaseangle.m`,
`plot_fig2_probarousal.m`, `plot_fig2_sigmatrace.m` (→ into `.../aim2_nrem_sleep_stability/plotting/`).
- **Flag for confirmation:** `prearousal_permutation_models.m` (non-`_optimized` version) — appears
  to be an earlier, slower implementation superseded by `_optimized`; only self-referenced (no
  external call sites found). Propose moving to `code/archive/analysis/` rather than
  `group-level/aim2.../`. Confirm before moving.

**Category J — move to `code/group-level/aim3_rem_transitions/`** (private helpers of `css_analyse_3`):
`transrem_permutation_models_optimized.m` (active), `plot_fig3_ampdur.m`, `plot_fig3_avsigmatrace.m`,
`plot_fig3_mdl.m`, `plot_fig3_mdlhists.m`, `plot_fig3_panellabels.m`, `plot_fig3_sigmatrace.m`
(→ into `.../aim3_rem_transitions/plotting/`).

**Actions:**
- [x] Re-verified each "Flag for confirmation" item via `search_codebase` before moving — all
      still showed zero call sites in `code/` outside this plan doc, matching original planning
      findings. Batch-confirmed with you in one round: **all 5 proposals approved as documented**
      (infraslowmodpowerspect.m/_old.m/extractismphase.m → archive/scratch; tmp_spindle_detection.m
      → archive/scratch, convspindles.m → utils/spindle-detection/ with NOTE; the 8 unreferenced
      plot/QC scripts → archive/scratch; withinConditionNorm.m/calcarovars.m/calcpsgvars.m/
      loadarousalcsvs.m/loadarousalepochs.m/inspect_n2aros.m → kept in utils/misc/ with NOTE
      comments added; prearousal_permutation_models.m (non-optimized) → archive/analysis).
- [x] Executed `git mv` for all 95 `code/supportfunc/*.m` files per the category tables (A–J plus
      all flagged items per the confirmed proposals above): 8 → `utils/eeg-io/`, 15 →
      `utils/signal-processing/`, 7 → `utils/signal-append/`, 7 → `utils/isf-fitting/` (+3 flagged
      → `archive/scratch/`), 7 → `utils/spindle-detection/` (incl. `convspindles.m`; +1 flagged
      `tmp_spindle_detection.m` → `archive/scratch/`), 8 → `utils/circular-stats/`, 7 →
      `utils/plotting-generic/` (+8 flagged → `archive/scratch/`), 8 → `utils/misc/` (incl. the 6
      flagged-but-kept files), 8 → `group-level/aim2_nrem_sleep_stability/` (1 active +7
      plotting), 7 → `group-level/aim3_rem_transitions/` (1 active +6 plotting), 1 flagged →
      `archive/analysis/` (`prearousal_permutation_models.m`). Total: 67 files landed in
      `code/utils/`, 20 in `code/group-level/` group-2/3 helper folders, 8 net-new in
      `code/archive/{scratch,analysis}/`. All via `git mv` (confirmed as renames in
      `git status --short`, history preserved).
- [x] Added `% NOTE:` comments (call-site status + rationale for not archiving) to the 7 files
      kept despite having no confirmed call sites: `convspindles.m`, `withinConditionNorm.m`,
      `calcarovars.m`, `calcpsgvars.m`, `loadarousalcsvs.m`, `loadarousalepochs.m`,
      `inspect_n2aros.m`.
- [x] `grep -rn "\b<oldfilename_without_ext>\b" code/` per category — no broken references found;
      all moved functions still resolve via `addpath(genpath('code'))` (folder-agnostic dispatch).
      Final sweep `grep -rn "supportfunc" --include='*.m' .` across the whole repo returned zero
      matches (only this plan file still mentions the word, in already-completed-phase notes).
- [x] Deleted now-empty `code/supportfunc/` (fully drained; `rmdir` succeeded).

**Phase 3 completed 2026-09-18.** All moves via `git mv` (history preserved, confirmed via
`git status --short` showing `R` for all 95 files). No statistical method, model specification,
threshold, filter parameter, or numeric result was changed — pure categorisation/relocation plus
the approved NOTE-comment additions on 7 kept-but-unreferenced files.

**Skills/tools:** `search_codebase` (re-verify call sites), `run_commands` (`git mv`, `grep`),
`ask_question` (batch-confirmed flagged items with you in one round). No MCP server needed.


---

## Phase 4 — Group-level analysis reorganisation + script→function conversion

**Read first (per sub-item, re-read fully before editing since these are the most complex edits):**
`@code/analysis/css_analyse_1a.m`, `@code/analysis/css_analyse_1b.m`,
`@code/analysis/css_analyse_1c.m`, `@code/analysis/css_plot_1c.m`,
`@code/analysis/css_analyse_2a.m`, `@code/analysis/css_analyse_2b.m`,
`@code/analysis/css_analyse_3.m`, `@code/analysis/css_descriptives.m`,
`@code/analysis/css_analyse_isf_phase_distribution.m`, `@documents/manuscript.md` (for aim/figure
mapping confirmation — read-only reference, not edited).

**4a. `css_analyse_1a.m` → `code/group-level/aim1_isf_characterization/analyse_bout_selection.m`**
- Already a proper function (`function css_analyse_1a()`); pure rename + move, no logic change.
- Update internal `fprintf` report text: none needed (no filename self-reference).

**4b. `css_analyse_1b.m` → `code/group-level/aim1_isf_characterization/analyse_sigma_spindle_similarity.m`**
- Convert script to `function analyse_sigma_spindle_similarity(SigmaFiles, SpdFiles)` (or, to
  preserve current behaviour exactly, `function analyse_sigma_spindle_similarity()` that
  internally does the same `dir(...)` calls the script currently does for `SigmaFiles`/`SpdFiles`
  — **preserve exact current file-discovery glob patterns**).
- **Bug fix (approved per Q5):** replace `for f = 1:length(Files)` (undefined `Files`) with
  `for f = 1:length(SigmaFiles)` — this is what the surrounding code clearly intends (loop bound
  should match `SigmaFiles`/`SpdFiles`, both same length by construction). State this explicitly
  in the phase-completion summary as a bug-fix, not a methodology change.
- No other logic changes.

**4c. `css_analyse_1c.m` → `code/group-level/aim1_isf_characterization/analyse_isf_topography.m`**
- Already a function; pure rename + move.
- Update its call to `css_plot_1c('abs')` / `css_plot_1c('norm')` to the renamed
  `plot_isf_topography('abs')` / `plot_isf_topography('norm')` (see 4d).

**4d. `css_plot_1c.m` → `code/group-level/aim1_isf_characterization/plot_isf_topography.m`**
- Already a function; pure rename + move.

**4e. `css_analyse_2a.m` → `code/group-level/aim1b_filter_validation/analyse_filter_edge_artefact.m`**
- Currently `% function css_analyse_2a` (commented out, run as script with `clear`/`clc` at top).
  Convert to `function analyse_filter_edge_artefact()` — no external inputs are assumed beyond
  what the script itself loads via `dir(...)`, so this is a straightforward wrap, not a signature
  change requiring new arguments.

**4f. `css_analyse_2b.m` → `code/group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m`**
- Currently `%function css_analyse_2` (mislabelled comment — even the commented-out name is wrong,
  says `css_analyse_2` not `css_analyse_2b`; flag this as a pre-existing documentation typo, fix
  when renaming). Convert to `function analyse_arousal_isf_phase()`.
- Update internal calls: `prearousal_permutation_models_optimized` (path-transparent, no rename),
  `plot_fig2_*` functions (path-transparent, no rename — only their *file location* changed in
  Phase 3, not their function names, so call sites need no edits).
- **Rename per Q1 (distinct-from-descriptives labelling):** this is the correct file to carry the
  "NREM sleep stability following cortical arousals" aim label (matches manuscript section
  "Sleep stability following cortical arousals" / Fig 3) — name chosen accordingly.

**4g. `css_analyse_3.m` → `code/group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m`**
- Currently a bare script (no `function` line at all, not even commented). Convert to
  `function analyse_rem_transition_dynamics()`.
- Internal call to `transrem_permutation_models_optimized` and `plot_fig3_*` — path-transparent,
  no edits needed beyond location (Phase 3).


**4h. `css_descriptives.m` → `code/group-level/sleep_macroarchitecture/analyse_sleep_macroarchitecture_and_arousal_outcomes.m`**
- Currently a bare script reading `phenotype/2024-07-26T171508_psg-variables.csv` and
  `..._arousals.csv`. Convert to
  `function analyse_sleep_macroarchitecture_and_arousal_outcomes()` — matches manuscript sections
  "Sleep macroarchitecture" (Table 1 GLMMs) and part of "Sleep stability following cortical
  arousals" (awakening-probability GLMM, Table 2) per Q1's instruction to name for the specific
  aim it addresses rather than "descriptives".
- Two local helper functions already exist at file end (`stringifypvalue`, `calcRiskRatio`) —
  these are MATLAB local functions and will move with the file as-is (valid syntax, no change
  needed); confirm they remain below the main function body per MATLAB's local-function rules.

**4i. `css_analyse_isf_phase_distribution.m` → `code/group-level/supplementary/plot_isf_phase_distribution_illustration.m`**
- Bare script; loads `analysis_2a.mat` (from `permutation-runs/`, gitignored) for `ANGA`.
  **Flag rather than fix:** I could not confirm during planning that `css_analyse_2a.m` actually
  saves a variable named `ANGA` into a file called `analysis_2a.mat` — grep found no `save(...,
  'ANGA'...)` or similar in `css_analyse_2a.m`. This may mean the `.mat` file was produced by an
  even earlier, now-archived version, or created manually. **Do not silently fix or guess** — wrap
  the file load in a function `function plot_isf_phase_distribution_illustration()` (pure
  wrap, no logic change) and add a `% NOTE:` comment flagging the unclear provenance of
  `analysis_2a.mat` for your attention, rather than attempting to regenerate or rename the
  dependency.


**Validation (per file, since no MATLAB CLI is available here):**
- [ ] Side-by-side read of old vs. new file content to confirm 100% logic parity aside from the
      one approved bug-fix (4b) and the function-wrapping itself.
- [ ] Confirm all `end` keywords / brace matching are still valid after adding the wrapping
      `function ... / end` pair (MATLAB scripts don't require function-body `end` for the
      top-level, but do for nested `if`/`for` — check no double-`end` introduced).
- [ ] **Recommended external step (you, locally):** open each converted file in MATLAB, run
      `checkcode('<file>.m')` to confirm no syntax errors, then run the function with a small
      subset of subjects/derivatives available on the `sleep` volume and confirm identical console
      output / figure output to a pre-refactor run.

**Actions:**
- [x] Execute all `git mv` + function-wrapping edits per 4a–4i.
- [x] Update every cross-reference between these files (e.g., 4c calling 4d) to use new names.
- [x] `grep -rn` for every old function name (`css_analyse_1a`, `css_analyse_1b`, ..., `css_descriptives`,
      `css_plot_1c`, `css_analyse_isf_phase_distribution`) across `code/` and confirm the only
      remaining reference is in `code/main.m` (updated in Phase 5) and this plan file.

**Phase 4 completed 2026-09-18.** Notes:
- All 9 `git mv` renames done via `code/analysis/*.m` → `code/group-level/<aim-folder>/*.m`;
  `git status --short` confirms all 9 as `RM` (rename+modify), preserving history.
- 4a (`analyse_bout_selection.m`) and 4c (`analyse_isf_topography.m`) were listed as "pure
  rename + move" in this plan, but since MATLAB requires a file's primary function name to match
  its filename to be callable, their internal `function css_analyse_1a()` / `function
  css_analyse_1c()` declarations were also renamed to `function analyse_bout_selection()` /
  `function analyse_isf_topography()` respectively (beyond just the file move) — noting this as
  a deviation-in-degree from "pure move", though it introduces no logic change, matching the
  spirit of every other 4b–4i conversion in this phase.
- 4b bug-fix applied exactly as specified: `for f = 1:length(Files)` → `for f =
  1:length(SigmaFiles)` in `analyse_sigma_spindle_similarity.m`, with an inline comment marking
  it as a bug-fix.
- 4f's mislabelled `%function css_analyse_2` comment corrected to `function
  analyse_arousal_isf_phase()` as specified.
- The pre-existing `N%%` stray-token typo in `analyse_filter_edge_artefact.m` (line 248, a
  malformed cell-divider) was found again during Phase 4 and **left untouched** — out of scope
  per the plan's "no logic changes beyond the one approved 4b bug-fix" rule.
- 4i's `% NOTE:` comment was **not** copied verbatim from this plan's original wording, because
  re-verification during Phase 4 found `analyse_filter_edge_artefact.m` (formerly
  `css_analyse_2a.m`) *does* contain `save('analysis_2a.mat', ..., 'ANGA', ...)`, contradicting
  this plan's claim that "grep found no `save(..., 'ANGA'...)`". The NOTE actually written flags
  the real ambiguity instead (working-directory-relative save path vs. expected load location),
  rather than repeating a now-inaccurate claim.
- All 9 converted files spot-checked for balanced `function`/`end` pairing (via
  `grep -c '^end$'` / `grep -c '^function'` and manual line-context review) — no double-`end` or
  missing `end` introduced by the wrapping edits.
- `code/analysis/` directory removed (now empty, fully drained).

**Skills/tools:** `read_files` (side-by-side comparison), `editor` (function wrapping, renames),
`run_commands` (`git mv`, `grep`). No MCP server needed. If available to you locally, MATLAB's
built-in `checkcode` / Code Analyzer is the relevant "skill" for validation — no external
download needed, it ships with MATLAB.


---

## Phase 5 — `main.m` update to reflect full pipeline

**Read first:** `@code/pipeline/main.m` (post Phase 2 move), `@documents/manuscript.md` (for aim
ordering/labels), and the final renamed files from Phase 4 to get their exact new names/paths.

**Actions:**
- [ ] In the `% ANALYSE:` section of `main.m`, replace the current 5 calls
      (`css_analyse_1a(); css_analyse_1b(); css_analyse_1c(); css_analyse_2a(); % css_analyse_2b();
      TODO...`) with calls to the renamed functions, in this order (matching manuscript aim order):
      `analyse_bout_selection();` → `analyse_sigma_spindle_similarity();` →
      `analyse_isf_topography();` → `analyse_filter_edge_artefact();` →
      `analyse_arousal_isf_phase();` → `analyse_rem_transition_dynamics();` →
      `analyse_sleep_macroarchitecture_and_arousal_outcomes();`
      (plus, as a clearly-labelled supplementary/illustrative call,
      `plot_isf_phase_distribution_illustration();` at the end).
- [ ] Preserve every existing explanatory comment block above each call (these already describe
      the scientific purpose well) — only update the function name being called, not the prose.
- [ ] Remove the stale `% TODO. Run the statistical analysis...` comment above the former
      `css_analyse_2b` line, since that TODO is now resolved by actually wiring the call in.
- [ ] Update the top-of-file comment referencing `DOI: xxx.x.x.x` — leave as-is (out of scope;
      not a structural issue, and I don't have the actual DOI to fill in).

**Skills/tools:** `editor` only. No MCP server needed.


---

## Phase 6 — README and final verification

**Read first:** `@README.md`, `@code/pipeline/main.m` (final state), all files moved/created in
Phases 1–5 (for the traceability table).

**Actions:**
- [ ] Update `README.md`'s "Project Structure" tree to reflect the new `code/` layout
      (`pipeline/`, `subject-level/`, `firstlevel-outcomes/`, `group-level/`, `utils/`, `qc/`,
      `toolboxes/`, `archive/`), replacing the old `analysis/processing/qc/supportfunc/archive`
      list.
- [ ] Add a new `code/README.md` with two tables required by the original task's traceability
      goal:
      1. **First-level outcomes table**: for each file in `code/firstlevel-outcomes/`, columns =
         {function, what it's derived from, unit/interpretation, where saved (derivative folder +
         filename pattern), which group-level analysis consumes it}.
      2. **Group-level analysis → manuscript mapping table**: for each folder in
         `code/group-level/`, columns = {folder, manuscript aim/section, key figure/table
         number(s), main statistical method used}.
- [ ] Run a final repo-wide `grep -rn` sweep for every old top-level folder name
      (`code/processing/`, `code/analysis/`, `code/supportfunc/`, `code/tmp/`) to confirm zero
      remaining path-string references anywhere in `code/` (comments included).
- [ ] Confirm `code/processing/`, `code/analysis/`, `code/supportfunc/`, `code/tmp/` directories
      no longer exist (fully drained and removed) and `code/archive/` is the single consolidated
      archive location.
- [ ] **Recommended external step (you, locally):** run
      `matlab.codetools.requiredFilesAndProducts('code/pipeline/main.m')` before and after the
      full refactor (keep the "before" list from a pre-refactor git stash/branch) and diff the two
      dependency lists — they should be identical in content (function names / logic), differing
      only in file paths.

**Skills/tools:** `editor`, `run_commands` (`grep`, `git status`). No MCP server needed.


---

## Deferred / Optional (not part of this refactor unless you request it later)

- Splitting `css_infraslowfluctpowerspect.m` / `css_crosscorr.m` into separate "compute" and
  "save first-level output" functions (currently one function does both) — flagged in Phase 2,
  not executed, since it would be a logic restructuring beyond pure reorganisation and the task
  says to avoid premature refactoring.
- Splitting `run_processing_section.m` (formerly `section.m`) into an infrastructure layer
  (parallel host assignment, existence checks) and a step-registry layer (the big switch-case) —
  flagged during planning as a SRP violation, but deferred since it's a larger behavioural
  surface to validate without a MATLAB runtime available in this environment.
- Any change to `css_getsigmapowerusingwavelet.m` vs `css_getspecpowerusingwavelet.m` duplication
  — flagged, not merged, since their divergence (hardcoded vs. parametrised frequency bands) may
  be intentional.

---

## Execution Checklist

- [x] Phase 0 — Setup (2026-09-18)
- [x] Phase 1 — Archive consolidation (2026-09-18)
- [x] Phase 2 — Subject-level + first-level reorganisation (2026-09-18; correction added
      2026-09-18 re: manual edits post-completion — see Phase 2 section)
- [x] Phase 3 — utils/ categorisation (+ batch confirmation of flagged items) (2026-09-18)
- [x] Phase 4 — Group-level reorganisation + script→function conversion (2026-09-18)
- [ ] Phase 5 — main.m update
- [ ] Phase 6 — README + final verification


