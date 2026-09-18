# SRP / Naming / Dead-code Refactoring Plan

**Status:** IN PROGRESS. Update the checklist at the bottom as each phase completes.

**Scope:** `code/` directory only (excludes `code/toolboxes/`, which is third-party and
untouched). This plan is the follow-up to `code/REFACTOR_PLAN.md` (file/folder reorganisation,
completed 2026-09-18). That plan reorganised *where* code lives; this plan addresses three
things it explicitly deferred or didn't cover:

1. A full static call-graph of the codebase, to find dead/unreachable code.
2. Consistent `css_`-prefix usage: prefix only the functions that are directly called as
   top-level entry points from `main.m` (subject-level dispatch via `run_processing_section.m`'s
   switch-case, or group-level `ANALYSE:` calls); everything called *within* those entry points
   should not carry the prefix.
3. Short, clear, `lower_snake_case` function names throughout, and Single Responsibility
   Principle decomposition of large multi-purpose files.

**Goal (same discipline as REFACTOR_PLAN.md):** rename, decompose, and archive without changing
any statistical method, model specification, threshold, filter parameter, or numeric result.
Preserve `rickwassing` coding style (spacing, blank-line block separation, purpose-comments).
Use `git mv` for every relocation so history is preserved. Update every call site whenever a
function is renamed or moved.

**Tooling notes (carried over from REFACTOR_PLAN.md):** No MATLAB CLI is available in this
environment, so numeric/behavioural validation of renamed/decomposed functions cannot be run
here. Validation instead means: static review (read every call site after each rename), a
Python-based regex call-graph tool (built for this plan, see Phase 1) for dependency
verification, and `grep`-based sweeps to confirm zero dangling references. If MATLAB is
available locally, running `matlab -batch "addpath(genpath('code')); main"` (or individual
group-level functions) and diffing generated output against a pre-refactor copy is recommended.

---

## Phase 0 — Setup

- [x] Write this plan file (done).
- [x] Confirm the codebase is unchanged since the audit that informed this plan (re-ran
      `find code -iname '*.m'` and diffed against the cached file list from planning — identical).

**Skills/tools:** `run_commands` (`find`, `diff`). No MCP server needed.

---

## Phase 1 — Call-graph audit, dead-code archival, bug fixes

**Method:** A Python script parsed every `function` declaration (handling `...` line
continuations, `[a,b] = name(...)`, `a = name(...)`, and bare-script files with no `function`
keyword) across all 141 `.m` files in `code/` (excluding `toolboxes/`), building a name → file
map. A second pass grepped every file for word-boundary references to every known function name,
then did a reachability walk from the three real entry points: `main.m`, `css_init.m`, and the
manually-run `qc/*.m` scripts (these are interactive/manual QC tools, not called from `main.m`,
but are legitimate entry points a user runs directly).

**Full results:** see `code/CALL_GRAPH.md`.

**Actions:**
- [x] Built the call-graph and reachability report (`code/CALL_GRAPH.md`).
- [x] Fixed 2 real bugs found as a side-effect of reading every file for the audit (not part of
      the original ask, but unambiguous syntax/logic errors that would break execution):
      1. `main.m` line 147: stray trailing `n` character (`cfg.do_parallel = do_parallel;n`) —
         would throw a MATLAB parse error if that section of `main.m` were run. Removed the `n`.
      2. `subject-level/segmentation/css_extractarousalbouts.m` line 69: a bare `keyboard`
         breakpoint left in the main execution path (after `getarousalbouts`, before the
         bout-extraction loop) — would halt any real run of this function. Removed.
- [x] Archived 18 functions confirmed to have **zero call sites anywhere in `code/`** (verified
      by grep across the whole tree, not just the reachability walk, and confirmed absent from
      `README.md`, `documents/manuscript.md`, and `documents/peer-review-comments.md` too — i.e.
      not referenced by name anywhere in the repo outside their own definitions):
      - `utils/circular-stats/withinChanCircMedian.m`
      - `utils/eeg-io/eeglab2arousals.m`
      - `utils/eeg-io/revisemistakesineventtable.m`
      - `utils/plotting-generic/plothypno.m` (superseded by EEG_Processor's `plotHypnogram`,
        which every active caller in `code/` actually uses)
      - `utils/signal-processing/csapsGCV.m`
      - `utils/signal-processing/eeg2freqbandpower.m`
      - `utils/signal-processing/gethilbert.m`
      - `utils/signal-processing/magnituderesponse.m`
      - `utils/signal-processing/morletgabortransform.m`
      - `utils/signal-processing/nonparzscore.m`
      - `utils/signal-processing/normdistance.m`
      - `utils/signal-processing/phaseresponse.m`
      - `utils/signal-processing/shorttimefft.m`
      - `utils/signal-processing/smooth1q.m`
      - `utils/signal-processing/zscoreacrosschannels.m`
      - `utils/spindle-detection/f_IFO_Parameters_140324.m` (superseded by
        `css_infraslowfluctpowerspect.m`/`fitisfspect.m`, the ISF-fitting pipeline actually wired
        into `main.m`)
      - `utils/spindle-detection/f_IFO_WithSpindles.m` (same supersession as above)
      - `utils/spindle-detection/f_MGT.m` (only used by the two `f_IFO_*` files above, forming a
        fully self-contained dead cluster)
      All moved via `git mv` to `code/archive/scratch/` (confirmed as renames in
      `git status --short`, history preserved), each with an `% ARCHIVED 2026-09-18: ...` header
      comment documenting why, matching the convention established in `REFACTOR_PLAN.md` Phase 1.
- [x] Re-ran the reachability walk after archiving: 92 non-archive files remain, 85 reached, 7
      unreached — and the 7 unreached are **exactly** the files already flagged-and-kept with a
      `% NOTE: no current call sites found...` comment during the original `REFACTOR_PLAN.md`
      Phase 3 (`calcarovars.m`, `calcpsgvars.m`, `inspect_n2aros.m`, `loadarousalcsvs.m`,
      `loadarousalepochs.m`, `withinConditionNorm.m`, `convspindles.m`) — confirming no live code
      was accidentally archived and no new unreferenced files were missed.
- [x] Did **not** archive `utils/eeg-io/eeglab2hypnogram.m` in this phase even though it's a
      near-duplicate of `css_eeglab2hypnogram.m` with one stray caller — that's handled
      separately in Phase 2 below, since it requires a call-site fix (not just archival).

**Phase 1 completed 2026-09-18.**

**Skills/tools:** `run_commands` (`find`, `python3` for the call-graph script, `git mv`, `grep`),
`editor` (bug fixes, archive header comments). No MCP server needed.

---

## Phase 2 — Fix the `eeglab2hypnogram` duplicate

`utils/eeg-io/eeglab2hypnogram.m` and `utils/eeg-io/css_eeglab2hypnogram.m` were byte-for-byte
identical logic under two names. Every active caller in `code/` used `css_eeglab2hypnogram`
**except** `subject-level/spectral-features/css_extractspindles.m`, which called the unprefixed
`eeglab2hypnogram` instead — almost certainly an oversight from an earlier edit.

**Actions:**
- [x] Changed `css_extractspindles.m` line 19 from `eeglab2hypnogram(EEG)` to
      `css_eeglab2hypnogram(EEG)`. Confirmed via `grep` this is now the only reference to either
      name in that file, and both functions had identical logic, so this is a pure clean-up with
      no numeric/behavioural change.
- [x] Archived the now-fully-unreferenced `utils/eeg-io/eeglab2hypnogram.m` to
      `code/archive/scratch/eeglab2hypnogram.m` via `git mv` (history preserved), with a header
      comment explaining it was a duplicate of `css_eeglab2hypnogram.m` and that its only caller
      was repointed rather than the duplicate being kept.
- [x] Re-ran the reachability check to confirm the fix didn't change behaviour: non-archive file
      count dropped by exactly 1 (92→91 total, 85→84 reached), and the 7 unreached files are
      still exactly the same intentionally-kept set as after Phase 1 — no live code broken, no
      new dead code introduced.

**Phase 2 completed 2026-09-18.**

**Skills/tools:** `editor` (call-site fix, archive header comment), `run_commands` (`git mv`,
`grep`, re-running the Phase 1 call-graph script). No MCP server needed.

---

## Phase 3 — Naming convention pass

Renamed all non-`lower_snake_case` function names, updating every call site. Per your explicit
decision, this included the vendor/third-party-derived spindle-detection functions (no exception
for provenance — they're fully ours to maintain now). A re-run of the Phase-1 name-classifier
script during execution turned up one additional non-compliant name not caught during planning
(`averageSigmaPerEvent`) — added to the batch below.

**Renames executed:**

| Old name | New name | File |
|---|---|---|
| `averageSigmaPerEvent` | `average_sigma_per_event` | `utils/circular-stats/averageSigmaPerEvent.m` → `average_sigma_per_event.m` |
| `withinChanCircMean` | `within_chan_circ_mean` | `utils/circular-stats/withinChanCircMean.m` → `within_chan_circ_mean.m` |
| `correctPhaseByEmpiricalCDF` | `correct_phase_by_empirical_cdf` | `utils/circular-stats/correctPhaseByEmpiricalCDF.m` → `correct_phase_by_empirical_cdf.m` |
| `permuteEventLabels` | `permute_event_labels` | `utils/circular-stats/permuteEventLabels.m` → `permute_event_labels.m` |
| `predictISF` | `predict_isf` | `utils/isf-fitting/predictISF.m` → `predict_isf.m` |
| `withinConditionNorm` | `within_condition_norm` | `utils/misc/withinConditionNorm.m` → `within_condition_norm.m` |
| `withinSubMean` | `within_sub_mean` | `utils/misc/withinSubMean.m` → `within_sub_mean.m` |
| `plotFiltParams` | `plot_filt_params` | `utils/plotting-generic/plotFiltParams.m` → `plot_filt_params.m` |
| `f_SpDetection_Ferrarelli` | `detect_spindles_ferrarelli` | `utils/spindle-detection/f_SpDetection_Ferrarelli.m` → `detect_spindles_ferrarelli.m` |
| `f_SpDetection_Humans` | `detect_spindles_fernandez` (implements Fernandez et al. 2018 per its own header comment, despite the filename) | `utils/spindle-detection/f_SpDetection_Humans.m` → `detect_spindles_fernandez.m` |
| `f_SpDetection_Wamsley` | `detect_spindles_wamsley` | `utils/spindle-detection/f_SpDetection_Wamsley.m` → `detect_spindles_wamsley.m` |

(`csapsGCV`, `f_IFO_Parameters_140324`, `f_IFO_WithSpindles`, `f_MGT`, `withinChanCircMedian` were
already archived in Phase 1 as dead code — not renamed in place, since they're not live.)

**Actions:**
- [x] Renamed each file (`git mv`, history preserved) and its internal `function` line
      (including self-referential doc-comment headers, e.g. `%PERMUTEEVENTLABELS` →
      `%PERMUTE_EVENT_LABELS`) to match.
- [x] Updated every call site across `code/` (16 call sites across 8 files, grep-verified before
      and after — used `sed` only for exact, unambiguous, repeated-pattern replacements within a
      single file after visually confirming every occurrence matched, `editor` for all others).
- [x] Re-ran the call-graph/reachability script: 91 non-archive files, 84 reached, same 7
      intentionally-kept files unreached — no regressions.
- [x] Re-ran the naming-convention classifier: **0 non-compliant function names remain** in
      `code/` (excl. `toolboxes/`) — down from 11 flagged at the start of this phase (16 minus 5
      that were already archived dead code from Phase 1, as noted above).
- [x] Full repo-wide `grep` sweep for all 11 old names: zero references remain outside
      `code/archive/`.

**Phase 3 completed 2026-09-18.**

**Skills/tools:** `editor` (function/file renames, call-site updates), `run_commands` (`git mv`,
targeted `sed` for verified-safe bulk replacements, `grep` sweeps, Python call-graph re-run). No
MCP server needed.

---

## Phase 4 — `css_` prefix correction

**Rule:** `css_` prefix = "this is a first-to-call entry point that executes processing/analysis
for the paper" (subject-level dispatch targets called from `run_processing_section.m`'s
switch-case, or group-level functions called directly from `main.m`'s `% ANALYSE:` section).
Everything called *within* those entry points keeps an unprefixed, descriptive name.

**Added `css_` prefix to** (previously unprefixed, called directly from `main.m`):
`analyse_bout_selection` → `css_analyse_bout_selection`,
`analyse_sigma_spindle_similarity` → `css_analyse_sigma_spindle_similarity`,
`analyse_isf_topography` → `css_analyse_isf_topography`,
`analyse_filter_edge_artefact` → `css_analyse_filter_edge_artefact`,
`analyse_arousal_isf_phase` → `css_analyse_arousal_isf_phase`,
`analyse_rem_transition_dynamics` → `css_analyse_rem_transition_dynamics`,
`analyse_sleep_macroarchitecture_and_arousal_outcomes` → `css_analyse_sleep_macroarchitecture`
(shortened per the naming-length rule),
`plot_isf_phase_distribution_illustration` → `css_plot_isf_phase_distribution` (shortened).

**Removed `css_` prefix from** (previously prefixed, but only called as internal helpers, never
as a `main.m`/dispatch entry point):
`css_extractfeatures` → `extract_isf_features` (called only by `css_infraslowfluctpowerspect`),
`css_createfstlvloutput` → `create_fstlvl_output` (helper called by 3 files),
`css_standard_colors` → `standard_colors` (plotting utility, not an entry point).

**Exception discovered during execution — `css_eeglab2hypnogram` keeps its prefix:**
The plan proposed renaming `css_eeglab2hypnogram` → `eeglab2hypnogram` since it's a pure utility,
not a dispatch target. This was attempted, but reverted: `main.m` runs `addpath(genpath('code'))`,
which puts `code/archive/` on the MATLAB path alongside everything else. The exact-duplicate
`eeglab2hypnogram.m` archived in Phase 2 (`code/archive/scratch/eeglab2hypnogram.m`) would then
collide on the path with a *newly un-prefixed* `utils/eeg-io/eeglab2hypnogram.m`, recreating the
ambiguous-shadowing problem Phase 2 specifically fixed. **Decision: `css_eeglab2hypnogram` is kept
as a documented exception to the Phase 4 rule**, purely to avoid this name collision with archived
code — not because it's actually an entry point.

**Borderline cases left as-is (per plan):**
`css_init` (top-level init, called directly by `main.m` — a legitimate entry point) and
`run_processing_section` (generic dispatcher, not itself an analysis step, but not `css_`-style
named either — no change made, consistent with the original refactor's naming for this file).

**Actions:**
- [x] Renamed each file (`git mv`, history preserved) + function declaration for all 8
      "add-prefix" and 3 "remove-prefix" functions (the 4th "remove-prefix" candidate,
      `css_eeglab2hypnogram`, was reverted per the exception above — net-zero change, confirmed
      via `git status` showing no diff for that file).
- [x] Updated every call site: `main.m`'s 8 `% ANALYSE:` calls, plus internal callers of
      `css_extractfeatures`/`css_createfstlvloutput`/`css_standard_colors` (2, 3, and 15 call
      sites respectively — the 15 `css_standard_colors(` → `standard_colors(` call sites were
      replaced via `sed` after visually confirming the pattern was unambiguous and exact across
      all matched files).
- [x] Full repo-wide `grep` sweep for all 11 old names (8 add-prefix + 3 successfully
      remove-prefix): zero references remain outside `code/archive/`.
- [x] Re-ran the call-graph/reachability script: 91 non-archive files, 84 reached, same 7
      intentionally-kept files unreached — no regressions.
- [x] Verified every remaining `css_*` file in `code/` (excl. `archive/`, `toolboxes/`) is either
      a subject-level dispatch target, a `main.m` `% ANALYSE:` entry point, `css_init`, or the
      documented `css_eeglab2hypnogram` exception — no stray/unjustified `css_` prefixes remain.

**Phase 4 completed 2026-09-18.** `code/README.md`'s traceability tables still need updating to
reflect the renamed functions — tracked as a Phase 6 action (final verification), not repeated
here to avoid duplicate/conflicting edits across phases.

**Skills/tools:** `editor` (function/file renames, call-site updates), `run_commands` (`git mv`,
targeted `sed`, `grep` sweeps, Python call-graph re-run). No MCP server needed.

---


## Phase 5 — SRP decomposition of large multi-responsibility files (not yet started)

Files identified as doing multiple distinct jobs in one function body (load data → compute →
fit/model → plot → export, all inline), to be split into single-purpose helper functions with a
slim orchestrating entry function. One file at a time, each reviewed before moving to the next,
since this is the highest-risk phase (real logic extraction, not just renaming/moving).

| File | Lines | Notes |
|---|---|---|
| `group-level/aim1_isf_characterization/plot_isf_topography.m` | 827 | One function, ~16 sequential figure-panel blocks separated by `% ---` comments — natural split points already exist. |
| `group-level/aim1b_filter_validation/analyse_filter_edge_artefact.m` | 428 | Load → filter-comparison compute → Hilbert phase/amplitude compare → plot, all inline. |
| `group-level/aim2_nrem_sleep_stability/analyse_arousal_isf_phase.m` | 417 | Load → wavelet/filter → permutation-model calls → multi-figure plotting, all inline. |
| `group-level/aim3_rem_transitions/analyse_rem_transition_dynamics.m` | 345 | Load → peak/trough extraction → `fitlme` modelling → plotting, all inline. |
| `group-level/sleep_macroarchitecture/analyse_sleep_macroarchitecture_and_arousal_outcomes.m` | 339 | Already has 2 local helper functions (`stringifypvalue`, `calcRiskRatio`) — the main body still mixes `fitlme`/`fitglme` modelling with plotting inline. |
| `pipeline/run_processing_section.m` | 256 | Mixes host-assignment/output-existence-checking infrastructure with the step-dispatch switch-case — flagged as deferred work in the original `REFACTOR_PLAN.md`'s "Deferred / Optional" section. |

**Planned approach per file:** read fully, identify natural seams (existing `% ---` block
comments are often already seam markers), extract each seam into its own unprefixed helper
function with a short descriptive `lower_snake_case` name, leave a slim orchestrating function
(named per Phase 4's `css_` rule) that just calls the helpers in sequence. No statistical method,
threshold, or numeric logic changes — pure decomposition.

**Planned actions (repeated per file above):**
- [ ] Read the full file, propose a decomposition (helper function boundaries + names).
- [ ] Extract helpers via `editor`, update the orchestrating function to call them.
- [ ] Grep-verify no other file called the now-internal logic directly (shouldn't be possible
      since these were single monolithic functions, but verify).
- [ ] Static review: confirm every line of original logic is preserved, just relocated.

---

## Phase 6 — Final verification (not yet started)

- [ ] Repo-wide `grep -rn` sweep for every renamed/removed symbol (old names) to confirm zero
      dangling references anywhere in `code/`.
- [ ] Re-run the full call-graph script one final time; confirm the reachable/unreached sets
      match expectations (same 7 NOTE-flagged files unreached, nothing else).
- [ ] Update `code/README.md`'s first-level-outcomes and group-level-mapping tables to reflect
      any renamed functions.
- [ ] Update this plan's Execution Checklist below.

---

## Execution Checklist

- [x] Phase 0 — Setup (2026-09-18)
- [x] Phase 1 — Call-graph audit, dead-code archival, bug fixes (2026-09-18)
- [x] Phase 2 — Fix `eeglab2hypnogram` duplicate (2026-09-18)
- [x] Phase 3 — Naming convention pass (2026-09-18)
- [x] Phase 4 — `css_` prefix correction (2026-09-18)
- [ ] Phase 5 — SRP decomposition of large files
- [ ] Phase 6 — Final verification

