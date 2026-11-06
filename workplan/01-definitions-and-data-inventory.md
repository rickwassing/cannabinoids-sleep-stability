# Phase 1 — Definitions and data inventory

## Status and progress

| Field | Status |
|---|---|
| Planning | Detailed plan approved |
| Implementation | Complete — inventories and handoffs validated |
| Last planning review | 2026-06-11 |
| Authorized file for this update | This phase file only |
| Analyses/manuscript/figure/response revisions completed | None |

- [x] Read the overarching plan, reviewer comments, both manuscripts, parent papers, and directly relevant project files.
- [x] Decompose relevant reviewer comments and define scope, dependencies, sources, decisions, risks, and acceptance criteria.
- [x] Agree on the detailed plan.
- [x] Build the reviewer traceability, source-of-truth, and analysis-population registers.
- [x] Complete the glossary, variable, metric, timing, subjective, formulation, and REM inventories.
- [x] Resolve or document outcomes/fallbacks for the remaining provenance checks.
- [x] Prepare later-phase handoffs and dual-manuscript location map.
- [x] Validate all deliverables and mark Phase 1 complete.

## Purpose

Establish the factual, operational, and data-provenance foundation needed to address the reviews. Phase 1 defines measurements and populations, inventories available variables and metadata, separates observations from interpretation, and records unresolved facts before analyses or manuscript revisions are designed or implemented.

The main outcome is a source-traceable specification that later phases can use without repeatedly reverse-engineering definitions, samples, timing, or derivations.

## Boundaries

### Included

1. Define arousal, awakening-associated and continued-sleep outcomes, stage assignment, neighboring-event handling, and temporal relationships.
2. Define ISFσ phase and convention, frequency, amplitude measures, half-period, ascending/descending phase, uniformity, and phase timing.
3. Verify formulas, denominators, units, and anchors for clinical sleep metrics and arousal indices.
4. Inventory treatment, sequence, period, washout, recording start/lights off, sleep onset, REM onset, and usable time-of-night variables.
5. Reconcile the parent and current samples and define analysis-specific populations.
6. Inventory parent-trial subjective measures and published findings available when participant-level data are inaccessible.
7. Establish dose/formulation facts and bounded literature questions for later interpretation.
8. Assess candidate healthy-control provenance and eligibility without assuming suitability.
9. Map definitions and variables to Phase 2 and later manuscript, interpretation, figure, synchronization, and response work.

### Excluded

- Running additional or sensitivity analyses (Phase 2).
- Resolving LC–NA mechanisms, pharmacology, clinical efficacy, or therapeutic interpretation (Phase 3).
- Rewriting/restructuring manuscripts (Phase 4), revising figures (Phase 5), final synchronization (Phase 6), or drafting the response letter (Phase 7).
- Treating candidate controls as a randomized placebo group or combining them with the randomized treatment contrast.
- Inferring actual drug concentrations or observed dose-to-sleep intervals from protocol-scheduled timing.

Phase 1 may specify later work but must not execute it.

## Source-of-truth hierarchy

Disagreements must be recorded and resolved, not silently ignored.

1. **Trial procedures/published parent findings:** `documents/parent-trial-protocol-paper.pdf` and `documents/parent-trial-results-paper.pdf`. Use for planned procedures, formulation, scheduled dosing, planned subjective measures, parent sample, and published outcomes. Protocol procedures are intended, not necessarily realized.
2. **Current analytical inputs/derivations:** `spss-analysis/2024-07-26T171508_psg-variables.csv`, `spss-analysis/2024-07-26T171508_arousals.csv`, `group-level/arousal-isf-phaseangle.csv`, and relevant current `code/`. Use for the n=19 sample, variables, transformations, and feature definitions. Confirm provenance where current code may postdate saved outputs or generation functions lack call sites.
3. **Linkage/QC metadata:** `Scoring log and notes_CANSLEEP arousal.xlsx`, `phenotype/recdates.csv`, and `phenotype/groups.csv`. Use for period/treatment linkage, sequence, recording start, washout, candidate controls, and historical notes. Historical notes are not evidence that an issue remained unresolved.
4. **Manuscripts:** `documents/manuscript_imrd.md` and `documents/manuscript_irdm.md` are the content being checked, not independent evidence when conflicting with higher-tier sources.
5. **Candidate controls:** conditional sources only until the eligibility gate is passed.

### Implemented source-of-truth and provenance register

This register records authority for Phase 1 facts. File modification times and abbreviated SHA-256 fingerprints were checked during implementation for core static inputs; fingerprints are provenance aids, not scientific validation.

| Source | Role and authoritative scope | Linkage/granularity | Status and limitations |
|---|---|---|---|
| `documents/peer-review-comments.md` | Complete reviewer wording and sub-question boundaries | Reviewer/comment/sub-comment | Authoritative review input; fingerprint `5bf4a2ec79ba187a` |
| `documents/manuscript_imrd.md` | Current IMRD claims and section locations | Manuscript section/paragraph | Revision target, not independent factual evidence; fingerprint `2da22c740e244e4c` |
| `documents/manuscript_irdm.md` | Current IRDM claims and section locations | Manuscript section/paragraph | Revision target; same scientific content must be retained; fingerprint `654d7274fe30e355` |
| `documents/parent-trial-protocol-paper.pdf` | Planned design, intervention, visit schedule, outcomes, and analysis plan | Trial level | Authoritative for intended procedures, not proof of realized visit-level timing; DOI 10.1136/bmjopen-2019-034421; fingerprint `afbd72c0e11ef10f` |
| `documents/parent-trial-results-paper.pdf` | Published parent n=20 outcomes, including subjective findings | Trial/condition level | Authoritative published parent result; not an n=19 subset analysis; DOI 10.1111/jsr.70124; fingerprint `39014789cb524bd4` |
| `Scoring log and notes_CANSLEEP arousal.xlsx` | Source-folder, visit, blinded/unblinded condition, and historical scoring/QC notes | Source recording/visit | Authoritative linkage aid; historical concerns require resolution status; R008 is confirmed corrected; fingerprint `8aee95a01979d7bf` |
| `phenotype/recdates.csv` | Observed PSG/EEG recording start, effectively lights off | `participant_id + session` | Complete and unique for all 38 analytical trial sessions; includes 18 candidate controls; fingerprint `223b88d1242c2857` |
| `phenotype/groups.csv` | Insomnia/control labels | Participant | Labels 19 trial participants and 22 candidate controls; control eligibility not established; fingerprint `d230ca155a47fc4a` |
| `spss-analysis/2024-07-26T171508_psg-variables.csv` | Current whole-night PSG and arousal-summary values | `participant_id + condition` | 38 rows/19 complete pairs; sleep efficiency and stage-percentage formulas numerically verified; fingerprint `3507a4ff30c0f33a` |
| `spss-analysis/2024-07-26T171508_arousals.csv` | Current all-night event inventory and event-level outcome fields | Event within participant/condition | 8,083 events across all classes; used for duration and stage/outcome inventory; fingerprint `366252df88af5e33` |
| `inspect/selectedaros.csv` | Manual/quality selection list for candidate arousal events | Candidate event ID/file | 1,451 rows; field names include leading spaces and status codes require documentation before independent reuse |
| `group-level/arousal-isf-phaseangle.csv` | Saved channel-expanded N2 phase-analysis event table | Event ID × channel | 22,734 rows = 1,263 unique events × 18 channels; all 19 participants; fingerprint `f771e257044a50cf` |
| `code/utils/eeg-io/getarousalbouts.m` and `code/subject-level/segmentation/css_extractarousalbouts.m` | Stage assignment, awakening classification, bout extraction, and ≥3-second event handling | Event | Authoritative executable definition subject to saved-output provenance checks |
| `code/utils/circular-stats/getarousalphaseangle.m`, `predict_isf.m`, `within_chan_circ_mean.m`, and correction/plot helpers | Neighboring-event selection, offset fields, Hilbert phase, event-level channel averaging, correction, and display convention | Channel/event | Authoritative current implementation; raw and corrected angular interpretations must remain separate |
| `code/utils/misc/calcpsgvars.m` and `calcarovars.m` | Candidate formulas for PSG variables and arousal indices | Participant/condition | No confirmed current generation call site; formulas must be validated against exports before being declared authoritative |
| `code/group-level/aim3_rem_transitions/` | REM extraction, fluctuation features, and pointwise window | Transition/channel/event | Authoritative current implementation; manuscript reports 118 transitions, while a compact saved transition table has not yet been identified |

#### Resolved source conflicts

1. **Parent n=20 versus current n=19:** parent paper governs parent outcomes; current paired HD-EEG inputs govern this analysis. R001's one failed treatment-night HD-EEG explains the difference.
2. **R008 historical note versus final exports:** the note predates review/restaging; final exports are authoritative.
3. **Preceding arousals “not excluded” versus 10-second rule:** the AASM rule requires at least 10 seconds of preceding sleep; events 10–60 seconds earlier remain included. The two statements address different intervals.
4. **Protocol dose time versus observed recording time:** protocol timing is nominal; `recdates.csv` is observed recording start/lights off. Neither supplies actual dose timestamps.
5. **`phase_d0` name versus physiological time:** the current variable name is retained in data/code, but its factual time is −1.5 seconds relative to target onset.

## Implemented reviewer-to-phase traceability matrix

Stable IDs below are the required keys for later analysis, manuscript, figure, and response-letter handoffs. “Phase 1 output” identifies the factual foundation; it does not authorize later-phase implementation.

| ID | Reviewer concern/sub-question | Phase 1 output | Response type | Later owner |
|---|---|---|---|---|
| R1-M1a | Insomnia-only sample may represent altered rather than general physiology | Population register and control eligibility gate | Data inventory; qualification | Phases 2–4, 7 |
| R1-M1b | Can healthy controls show normal versus insomnia sleep composition? | Control provenance/comparability decision | Conditional additional analysis | Phases 2, 7 |
| R1-M1c | Are pharmacological effects general or insomnia-specific? | Inferential boundary: untreated controls cannot test diagnosis-by-treatment interaction | Qualification/direct response | Phases 3, 7 |
| R1-M2a | Continuously timed arousals are paired with 30-second stage labels | Exact operational decision rule | Factual clarification | Phases 2, 4, 7 |
| R1-M2b | Arousal-to-classified-awakening interval varies across events | Epoch-position/temporal-bound variable specification | Sensitivity feasibility | Phase 2 |
| R1-M2c | Variable interval weakens causal interpretation | Epoch-associated terminology and causal limit | Qualification | Phases 3–4, 7 |
| R1-M2d | Shorter scoring or precise interval requested | Minimum clarification plus rescoring feasibility/fallback | Clarification/additional analysis | Phases 2, 7 |
| R1-M3a | Signal-processing terminology is inaccessible | Operational glossary | Factual clarification | Phase 4 |
| R1-M3b | Results prioritize methods over biological meaning | IMRD/IRDM location map | Manuscript restructuring | Phase 4 |
| R1-M4a | How can arousals occur during presumed low-LC/ascending sigma? | Measurement-to-interpretation boundary | Literature/qualification | Phase 3 |
| R1-M4b | Other arousal systems or a more complex phase–LC relationship may contribute | Glossary caveat and Phase 3 research handoff | Literature/qualification | Phase 3 |
| R1-M4c | Alternative noradrenergic interpretations may fit N1/REM/amplitude findings | Distinct observed-variable definitions | Literature/qualification | Phase 3 |
| R1-M5a | No global ISFσ frequency/amplitude/HR effect | Global versus event/state-specific measure definitions | Qualification | Phase 3 |
| R1-M5b | Does this argue against a broad LC mechanism? | Inferential boundary, without resolving mechanism | Direct response/qualification | Phases 3, 7 |
| R1-M6a | Dose time relative to lights off and sleep onset is unclear | Nominal-dose versus observed-lights-off inventory | Factual clarification | Phase 4 |
| R1-M6b | Do effects evolve from early to late night? | Available timing axes and Phase 2 input specification | Additional analysis | Phase 2 |
| R1-M6c | Would decay across night support a drug effect? | Pharmacokinetic limitation: no actual dose times/concentrations | Qualification | Phases 3, 7 |
| R1-m1a | Closely spaced arousals may carry over and be dependent | 0–10-second and 10–60-second rules | Factual clarification | Phases 2, 4 |
| R1-m1b | Repeat analyses with isolated arousals | Fully isolated 60-second sensitivity specification | Sensitivity analysis | Phase 2 |
| R1-m2a | Why use phase 1.5 seconds before onset? | Wavelet-contamination rationale and exact offset map | Factual clarification | Phase 4 |
| R1-m2b | Immediate pre-arousal changes may be meaningful | Distinguish state estimate from evoked contamination | Qualification | Phase 3 |
| R1-m2c | Test sensitivity to phase offset | Prespecified offset fields/times | Sensitivity analysis | Phase 2 |
| R1-m3 | Topographic frequency heterogeneity needs explanation | Regional ISFσ measurement definition | Literature/qualification | Phase 3 |
| R1-m4 | Opposing N1 and N2 effects need explanation | Stage-specific outcome/population definitions | Literature/qualification | Phase 3 |
| R1-m5 | What maximum arousal duration was permitted? | Formal rule, observed range, and 30-second cap inventory | Factual clarification | Phases 1, 4, 7 |
| R1-m6 | Figure 4E–F is difficult to interpret | REM variable definitions and figure handoff | Figure/text revision | Phases 4–5 |
| R1-m7 | Define sleep efficiency and other clinical metrics | Verified clinical-metric dictionary | Factual clarification | Phase 4 |
| R2-Ma | Lack of healthy controls | Control eligibility gate | Qualification/conditional analysis | Phases 2–4, 7 |
| R2-Mb | Add healthy ISFσ phase analysis if data exist | Control phase-variable availability and eligibility | Conditional additional analysis | Phase 2 |
| R2-C1a | Is 10 mg THC high or low? | Route/population/context fact sheet; no universal guideline | Literature/factual context | Phases 3, 7 |
| R2-C1b | Does CBD show comparable concentration-dependent LC/NA effects, and is 200 mg high/low? | CBD evidence question and uncertainty | Literature/qualification | Phases 3, 7 |
| R2-C1c | What is known about the 10/200-mg combination? | Fixed-combination formulation fact sheet and attribution limit | Literature/qualification | Phases 3, 7 |
| R2-C2 | Why n=20 previously and n=19 now? | R001 paired-HD-EEG exclusion fact sheet | Factual clarification/direct response | Phases 4, 7 |
| R2-C3 | Check treatment order | Sequence/period/washout metadata | Additional analysis | Phase 2 |
| R2-C4a | Move FIR details to Methods | Dual-manuscript location map | Manuscript revision | Phase 4 |
| R2-C4b | Explain phase angles/function and uniformity in main text | Glossary and raw/corrected convention | Factual clarification | Phase 4 |
| R2-C5 | Clarify interaction versus final-fluctuation post-hoc result | NREM/final-fluctuation and REM-latency definitions | Manuscript clarification | Phase 4 |
| R2-C6a | Is the approximately −60-second dip significant? | Current tested-window definition | Factual clarification | Phases 2, 7 |
| R2-C6b | Is the +30 to +60-second rise significant? | Untested-window status and new-analysis input | Additional analysis | Phase 2 |
| R2-C6c | Why use 90 seconds? | Pragmatic figure/multiplicity rationale | Direct response | Phases 4, 7 |
| R2-C7 | Improve Figure 4 colors/colorblind accessibility | No Phase 1 data need beyond variable identity | Figure revision | Phase 5 |
| R2-C8 | Add significance indications to figures | Analysis-result/figure handoff | Figure revision | Phase 5 |
| R2-C9 | Reconcile descending-phase continued sleep under placebo | Phase/outcome definitions and indirect LC mapping | Literature/qualification | Phase 3 |
| R2-C10a | TST/SE decreased despite reduced N2 awakening odds | Local-versus-whole-night metric distinction | Qualification | Phase 3 |
| R2-C10b | Why call THC/CBD therapeutic? | Clinical-efficacy boundary for a mechanistic acute study | Qualification/direct response | Phases 3, 7 |
| R2-C10c | What did participants report about sleep quality/amount? | Published n=20 subjective inventory and inaccessible-data limit | Factual clarification | Phases 1, 4, 7 |

### Traceability validation

- All Reviewer 1 major comments, all seven Reviewer 1 minor comments, Reviewer 2's major comment, and all ten numbered Reviewer 2 comments are represented.
- Figure-only concerns are routed to Phase 5 rather than silently omitted.
- Mechanistic questions are supplied with Phase 1 definitions but remain assigned to Phase 3.

## Implemented analysis-population register

Counts below are descriptive source inventories, not inferential analyses.

| Population | Participants/sessions/events | Inclusion and exclusion | Source/linkage | Status/limitation |
|---|---|---|---|---|
| Parent trial | n=20 participants | Parent-trial eligibility/completion | Parent protocol/results | Governs published parent objective/subjective findings |
| Current paired HD-EEG sample | 19 participants; 38 sessions (19/condition) | Complete usable HD-EEG in both crossover conditions; R001 excluded after one night failed | PSG export; participant + condition | Current manuscript-level analytical sample |
| Whole-night PSG macroarchitecture | 19 participants; 38 rows | One row per participant/condition | `...psg-variables.csv` | Complete paired table |
| All scored event export | 19 participants; 8,083 events: 4,506 placebo, 3,577 THC/CBD | Conventional arousal, alpha, and micro classes | `...arousals.csv` | Not all rows are conventional arousals |
| Conventional arousals | 19 participants; 6,090 events | `aro_type = 1arousal`; duration ≥ approximately 3 seconds | Event export | 15 events exceed 30 seconds before model capping |
| Conventional N1 arousals | 19 participants; 1,737 events (951 placebo; 786 THC/CBD) | Conventional events assigned to N1 | Event export | 480 awakening-associated; 1,257 continued-sleep |
| Conventional N2 arousals | 19 participants; 1,892 events (903 placebo; 989 THC/CBD) | Conventional events assigned to N2 | Event export | 1,376 awakening-associated; 516 continued-sleep; broader than phase subset |
| Conventional N3 arousals | 19 participants; 364 events (194 placebo; 170 THC/CBD) | Conventional events assigned to N3 | Event export | 152 awakening-associated; 212 continued-sleep |
| Conventional REM arousals | 19 participants; 2,097 events (1,331 placebo; 766 THC/CBD) | Conventional events assigned to REM (stage code 0) | Event export | 258 awakening-associated; 1,839 continued-sleep |
| Alpha events | 19 participants; 649 events | `aro_type = 2alpha` | Event export | Distinct from conventional arousals |
| Micro-events | 19 participants; 1,344 events | `aro_type = 3micro` | Event export | Distinct from conventional arousals |
| Candidate selected-arousal list | 1,451 rows | Manual/quality status in `inspect/selectedaros.csv` | Filename/event ID | Status-code semantics and attrition to phase table require documentation |
| Saved N2 phase-analysis events | 19 participants; 1,263 unique events; 22,734 channel rows | Selected N2 events surviving phase extraction; 18 channels/event | Event ID in `group-level/arousal-isf-phaseangle.csv` | 695 THC/CBD and 568 placebo events; not interchangeable with all 1,892 N2 arousals |
| N2 phase: awakening-associated | 693 unique events | Phase-table `is_awakening = 1` | Event ID | 368 THC/CBD; 325 placebo |
| N2 phase: continued sleep | 570 unique events | Phase-table `is_awakening = 0` | Event ID | 327 THC/CBD; 243 placebo |
| REM-transition analysis | 19 participants; 118 transitions | Algorithm-identified and visually confirmed transitions in extracted bouts | Manuscript/current REM workflow | 51 THC/CBD; 67 placebo; compact event table not yet located |
| Published subjective outcomes | Parent n=20 | Parent-results analysis population | Parent results paper | Participant-level data inaccessible; no n=19 reanalysis |
| Candidate controls labeled | 22 participants | `group = CNT` | `phenotype/groups.csv` | Eligibility and comparability unresolved |
| Candidate controls with verified processed/timestamp coverage | 18 participants; one recording each | Present in `recdates.csv` and relevant derivatives | Participant ID | Technical presence only; not an approved analysis sample |

### Population distinctions required downstream

1. Do not use counts from the all-event table as conventional-arousal counts.
2. Do not describe the 1,263-event phase subset as all N2 arousals; it is a selected/extracted subset of the 1,892 conventional N2 events.
3. Do not treat 18 processed controls as the eligible control sample until the gate is passed, or the 22 labeled controls as if all have usable recordings.
4. Do not attach parent n=20 subjective findings to the n=19 EEG subset.
5. For REM analyses, distinguish transitions (n=118) from participants (n=19) and account for repeated transitions within participant in later modeling descriptions.

## Agreed decisions and direct observations

### Trial, timing, and crossover

- Randomized, double-blind, placebo-controlled, two-period crossover trial.
- ETC120: plant-derived oral 2-mL formulation containing 10 mg THC and 200 mg CBD (1:20) in MCT oil; matched 2-mL cannabinoid-free placebo.
- Protocol-scheduled administration: one hour before typical sleep onset; intended 8-hour sleep opportunity.
- Actual visit-level dose times are unavailable. Describe dose timing as scheduled/nominal, not observed.
- `phenotype/recdates.csv` is observed PSG/EEG recording start and effectively lights off. It uniquely covers all 38 analytical sessions.
- Recording start, sleep onset, first REM onset, period, sequence, and washout are inventoryable; actual dose-to-lights-off/sleep intervals are not.
- Sequence: 8 placebo→THC/CBD and 11 THC/CBD→placebo. Inter-visit intervals are 7–170 calendar days; none is under one week.

### Samples and QC

- Parent trial: n=20 with analyzable parent outcomes.
- R001 completed the parent trial, but one treatment-night HD-EEG failed. Paired crossover HD-EEG analyses therefore exclude R001 entirely.
- Current sample: 19 complete pairs (`r002`–`r020`). This is analysis-specific data unavailability, not parent-trial attrition or outcome-based exclusion.
- R008's historical staging concern was reviewed/restaged; final exports contain corrected staging and are authoritative.

### Arousal outcome and neighboring events

- Arousals are continuously timed; sleep stages use 30-second epochs.
- An arousal is assigned to its onset epoch, except an onset in the first 15 seconds of a wake/N1 epoch preceded by the requested stage can be assigned to that preceding stage.
- For N2/N3/REM, awakening-associated means the next epoch is wake or N1; for N1, the next epoch must be wake.
- This is an epoch-associated operational outcome, not continuously timed proof that an arousal caused an exact awakening onset.
- Targets with another arousal in the preceding 10 seconds were excluded from phase extraction. Events 10–60 seconds earlier were retained and their pre-target probability modeled. The manuscript statement that events within 60 seconds were not excluded is accurate because the Sleep Scoring manual (The AASM scoring manual) states that arousals must preceed by at least 10 seconds of sleep.

### Arousal duration

- Conventional arousals span approximately 3.00–33.04 seconds in the all-night export.
- Fifteen of 6,090 conventional arousals exceed 30 seconds; six of 1,892 N2 conventional arousals exceed 30 seconds.
- The awakening-probability model caps durations above 30 seconds at 30 seconds. And that is fine.
- ISF event durations are stored in samples at 10 Hz, not seconds.

### ISFσ phase and timing

- Bouts extend to 30 seconds post-target and are cropped by 31.5 seconds before phase estimation. `phase_d0` is therefore −1.5 seconds, not onset. But these slight discrepancies are just in the variable name, and are taken into account with interpreting the results.

| Field | Nominal time relative to onset |
|---|---:|
| `phase_d0` | −1.5 s |
| `phase_d1` | −2.5 s |
| `phase_d2` | −3.5 s |
| `phase_d3` | −4.5 s |
| `phase_d4` | −5.5 s |
| `phase_d8` | −9.5 s |
| `phase_d12` | −13.5 s |
| `phase_d16` | −17.5 s |

- Channel phases are circularly averaged within event; events, not channels, enter circular tests.
- Raw Hilbert convention: 0° peak, 90° descending zero crossing, 180° trough, 270° ascending zero crossing, 360° peak. Raw descending sigma is approximately 0°–180°; ascending sigma is 180°–360°.
- Reported circular analyses apply an empirical-CDF correction for unequal background phase occupancy. Corrected angles form a uniform empirical null but are not automatically literal raw waveform radians. Document raw and corrected conventions separately and validate corrected display labels.

### REM windows

- REM-centered extraction is approximately −300 to +60 seconds for filtering, fluctuation detection, and visualization.
- The current pointwise treatment test is −90 to 0 seconds.
- The 90-second interval was a pragmatic figure-consistency/multiplicity choice, not a formal a priori biological boundary; do not retroactively call it prespecified.
- The reviewer-identified +30 to +60-second interval was not tested, but, this will be persued.

### Subjective outcomes

- Participant-level subjective data are inaccessible; the parent protocol/results are authoritative.
- Parent n=20 findings: no change in subjective sleep quality, a small increase in next-day self-reported sleepiness, and no change in next-day objective alertness.
- Do not present these as a new n=19 subset analysis. Keep perceived quality, perceived quantity, and next-day sleepiness distinct.

### Healthy controls

- Candidate control data exist but are not assumed usable.
- `phenotype/groups.csv` labels 22 controls, while 18 have timestamped/processed coverage; reconcile only if this stream is activated.
- A programmatic `placebo` session label does not make controls a CANSLEEP placebo arm.
- No control analysis may be promised before eligibility review.

## Reviewer issue map

### R1-M1 and Reviewer 2 major — Healthy controls/generalizability

- **Parts:** insomnia-only sample; potentially altered insomnia LC–NA physiology; request for healthy sleep composition and ISFσ phase; question of general versus insomnia-specific response.
- **Issue:** untreated controls can contextualize groups but cannot identify diagnosis-specific treatment modification.
- **Needs:** metadata/feasibility inventory, factual clarification, claim qualification, and conditional analysis.
- **Minimum:** preserve limitation and qualify general claims.
- **Strongest feasible:** exploratory matched/harmonized control versus insomnia-placebo comparison only after eligibility; no diagnosis-by-treatment inference.
- **Handoff:** Phases 2–4 and 7.

### R1-M2 — Variable arousal-to-awakening interval

- **Parts:** continuous arousal versus 30-second staging; variable temporal separation; causal concern; request for shorter scoring or precise definition.
- **Issue:** current label implies more temporal/causal precision than the rule supports.
- **Note:** that the peer reviewer incorrectly indicated that the awakening could occur up to 60 seconds after the start of the arousal. This is incorrect as the awakening could be at most 30 seconds after the start of the arousal.
- **Needs:** extract time delays between arousal and awake epochs to show these are within 0-30 second delays, and that these delays are not different between treatment conditions.
- **Preferred:** document the full rule and feasible interval bounds.

### R1-M3 and R2-C4 — Difficult terminology/results structure

- **Parts:** intuitive definitions for Hilbert phase, corrections, circular statistics, uniformity, and functional significance; technical detail obscures biology.
- **Needs:** glossary and measurement/interpretation boundary in Phase 1; restructuring in Phase 4.

### R1-M4–M5 and R2-C9 — Phase/LC–NA interpretation

- **Parts:** arousal during presumed low-LC states; other arousal systems; unchanged global ISFσ/HR metrics; descending-phase continued sleep under placebo.
- **Issue:** EEG quantities are conflated with direct LC firing/NA concentration.
- **Needs:** cautious measurement definition now; literature/mechanistic reconciliation in Phase 3.
- **Required language principle:** ascending/descending describe filtered sigma; LC–NA mapping is indirect and cannot be made, keep the interpretation at the level of EEG sigma power fluctuations.

### R1-M6 — Administration and time of night

- **Parts:** administration time, relation to lights off/sleep onset, early-versus-late effects, drug-mediated interpretation.
- **Needs:** timing inventory and Phase 2 feasibility; time-of-night analysis for example including time of night as a covariate and testing the interaction between treatment condition and time of night.
- **Minimum:** scheduled one hour before typical sleep onset; actual dose times unavailable.
- **Strongest:** observed recording/sleep-based time axis, without treating it as individual pharmacokinetics.

### R1 minor 1 — Isolated arousals

- **Parts:** preceding arousals may alter sigma and dependence; requested sensitivity.
- **Preferred handoff:** Do not address this minor comment.

### R1 minor 2 — −1.5-second phase estimate

- **Parts:** estimate is not at onset; immediate changes may matter; rationale and offset sensitivity requested.
- **Preferred handoff:** retain −1.5 seconds as contamination-avoidance primary and prespecify this offset was necessary because of the wavelet length.

### R1 minor 5 — Maximum arousal duration

- **Parts:** formal maximum, observed maximum, and missing reporting.
- **Preferred:** the AASM definition does not specify a maximum arousal duration. Raw events above 30 seconds were retained in the export, but duration was analytically capped at 30 seconds in the awakening-probability model; report the raw range and distinguish the model cap from a scoring maximum.

### R1 minor 7 and R2-C10 — Clinical metrics and perceived sleep

- **Parts:** define metrics/denominators; distinguish local stability from whole-night sleep; identify subjective quality/quantity.
- **Needs:** metric dictionary and published subjective inventory. Clinical/therapeutic interpretation belongs to Phase 3.

### R2-C1 — Dose/formulation

- **Parts:** whether THC/CBD doses are high/low; CBD LC/NA effects; evidence for the mixture.
- **Preferred:** contextualize that here are no dosing guidelines. Different by route, population, exposure, tolerance, timing, and outcome. Do not assign fixed-combination effects to THC, CBD, or interaction. Not sure we can address this comment.

### R2-C2 — n=20 versus n=19

- **Answer:** R001 completed the parent trial but one HD-EEG night failed; paired HD-EEG analyses include 19 complete pairs.

### R2-C3 — Treatment order

- **Parts:** order assessment and requested analysis; distinguish sequence, period, and carryover.
- **Preferred handoff:** participant-period table and parsimonious Phase 2 sensitivity. A null sequence result cannot prove no carryover in n=19.

### R2-C5 — NREM versus final pre-REM half-period

- **Need:** separately define ordinary NREM and final transition-to-REM fluctuations. Result rewriting belongs to Phase 4.

### R2-C6 — Figure 4 timing/significance

- **Parts:** pre-REM dip, post-REM rise, 90-second rationale.
- **Current fact:** corrected testing covers −90 to 0 and found no treatment difference; +30 to +60 was not tested; 90 seconds was pragmatic.
- **Preferred handoff:** retain existing interval for continuity and separately test post-REM only if approved, with exploratory label and multiplicity control.

## Manuscript assessment/location map

| Issue | IMRD | IRDM | Problem/type |
|---|---|---|---|
| Arousal outcome | Methods arousal/stability sections; NREM-stability Results | Results before full Methods | Rule, reassignment, and 10-second exclusion incomplete; temporal science limitation plus explanation/placement |
| Phase/offset | Phase Methods; Figure 3 Results/legend | Results before definition | “At onset” inaccurate; raw/corrected phase unclear |
| ISFσ terms | Methods, Results, legends, supplements | Same reordered | Frequency, amplitudes, phase, half-period conflated; explanation plus interpretation risk |
| Clinical metrics | Table 1 and Methods | Table/Results before Methods | Formulas/denominators absent at first use |
| Trial timing/order | Experimental Design/statistics | Later Methods | Scheduled timing, observed lights-off semantics, and order analysis missing |
| Sample | Abstract/Experimental Design | Abstract/later Methods | Paired HD-EEG rationale needs precision |
| Subjective sleep | Absent | Absent | Parent result unavailable to reader |
| Dose context | Introduction/Design | Same | Dose stated; universal high/low framing unsupported |
| REM variables/windows | REM Methods/Results/Figure 4 | Results before Methods | Two REM-latency constructs and windows need distinct names/context |
| Controls | Discussion limitation | Same | Candidate data exist but eligibility unknown |

IRDM Results will eventually need concise outcome/sample/measurement context before Methods, without duplicating full technical detail. Phase 1 identifies facts; Phase 4 drafts prose.

## Required inventories

### Operational glossary

Define each term with operational definition, unit, source, analysis use, biological limit, and manuscript locations. Include: cortical/conventional/alpha/micro arousal classes; awakening-associated and continued-sleep outcomes; stage assignment; recording start/lights off; sleep onset/period; TIB, TST, SOL, WASO, sleep efficiency, whole-night and per-cycle REM latency; stage minutes/percentages; arousal indices; instantaneous sigma power; filtered ISFσ; ISFσ frequency; fitted spectral amplitude; Hilbert amplitude; REM trough-to-peak amplitude; raw and corrected phase; phase landmarks; half/full period; uniform/non-uniform/unimodal/multimodal distributions; NREM versus transition fluctuation; nominal versus observed timing.

### Population register

Separate parent n=20, current n=19, macroarchitecture, all-arousal, conventional-arousal, stage-specific awakening, selected N2 phase, REM-transition, parent subjective, candidate-control, processed-control, and any eligible-control populations. Record participant/session/event counts, inclusion/exclusion, missingness, source, and linkage key.

### Clinical metric dictionary

| Metric | Provisional derivation to verify |
|---|---|
| Recording opportunity/TIB | Number of 30-second epochs × 0.5 min |
| TST | N1 + N2 + N3 + REM minutes |
| SOL | Recording start/lights off to first sleep epoch |
| WASO | Wake between first and last sleep epochs |
| Sleep efficiency | 100 × TST / recording opportunity |
| Whole-night REM latency | First REM relative to first sleep epoch |
| Stage percentage | 100 × stage minutes / TST |
| Overall arousal index | Conventional arousals / TST hours |
| NREM/stage index | Relevant arousals / relevant-stage hours |

Sleep efficiency and stage percentages already reproduce the export numerically; verify all remaining formulas before approval.

### Crossover/timing inventory

Required fields: participant, source visit, treatment, period, sequence, recording start, washout days, protocol minimum, SOL, derived sleep-onset clock time, first-REM latency/time, recording end if available, nominal scheduled dose context, actual-dose-time availability (`unavailable`), source, and validation status.

### Subjective inventory

For each published sleep-related measure record instrument/subscale, construct, assessment time, scoring, parent population, estimate/uncertainty if reported, source page/table, reviewer-question coverage, and unavailable-participant-data limitation.

### Dose/formulation fact sheet

Record product, route, volume, doses, ratio, vehicle, placebo, single-dose context, population/tolerance context, scheduled timing, protocol rationale, and inability to identify constituent-specific effects. Define—but do not conduct—Phase 3 research on oral exposure, dose context, LC/NA evidence, and combination-specific evidence.

### Healthy-control eligibility gate

Verify: source study/permissions; healthy screening; hardware/montage/sampling/channels; recording procedures; scoring/rater comparability; preprocessing; age/sex; medication/substance restrictions; adaptation status; available variables; missingness/QC; matching/harmonization; permissible inference. Outcome: eligible, conditionally eligible, ineligible, or unresolved. No unresolved analysis proceeds.

## Implemented operational glossary and variable inventory

### Arousal and sleep-stability terms

| Term | Operational definition in this project | Unit/granularity | Interpretation limit/source |
|---|---|---|---|
| Conventional cortical arousal | Scored abrupt EEG frequency shift satisfying the applied arousal criteria and lasting at least 3 seconds; export class `1arousal` | Event; seconds | Distinct from alpha and micro-event classes; scoring Methods and event export |
| Alpha event | Separately exported `2alpha` event class | Event; seconds | Do not combine with conventional arousals unless an analysis explicitly specifies it |
| Micro-event | Event shorter than the conventional 3-second threshold, exported as `3micro` | Event; seconds | Not a conventional arousal in the reported awakening models |
| Arousal stage | Stage containing onset, with possible reassignment to the preceding requested stage when onset occurs within the first 15 seconds of a wake/N1 epoch preceded by that stage | 30-second epoch plus continuous onset | This hybrid rule must be disclosed; `getarousalbouts.m` |
| Awakening-associated arousal | For N2/N3/REM-assigned events, the next stage epoch is wake or N1; for N1-assigned events, the next epoch is wake | Binary event outcome | Operational epoch association, not exact wake-onset latency or proof of causation |
| Continued-sleep arousal | Selected event not satisfying the awakening-associated rule | Binary event outcome | Does not mean the EEG remained unchanged or that no later awakening occurred |
| Immediate neighboring-event exclusion | At least 10 seconds of sleep must precede the target; targets with a preceding arousal in that interval are excluded | Event selection | Consistent with the applied AASM arousal requirement |
| Retained preceding event | Another event 10–60 seconds before target onset | Time-varying binary pre-target indicator | Retained in primary phase workflow and modeled; not equivalent to a fully isolated event |
| Fully isolated event | Proposed sensitivity definition: no preceding arousal within the complete 60-second pre-target window | Event selection | Phase 2 sensitivity definition, not the primary reported sample |
| Arousal duration | Scored event duration; raw conventional range approximately 3.00–33.04 seconds | Seconds; ISF table stores samples at 10 Hz | Awakening model caps values above 30 seconds; raw and modeled values must be distinguished |

### Clinical sleep terms

| Term | Verified operational definition | Unit/denominator | Notes |
|---|---|---|---|
| Recording opportunity / TIB | Full scored recording from effective lights off/start to recording end: number of epochs × 0.5 | Minutes | Called TIB in the export; operationally the analyzed recording opportunity |
| Total sleep time (TST) | All N1 + N2 + N3 + REM epochs × 0.5 | Minutes | Whole-recording sleep total |
| Sleep-onset latency (SOL) | Time from recording start/effective lights off to first scored sleep epoch | Minutes | First N1 or deeper sleep epoch |
| Sleep period | From first scored sleep epoch through last scored sleep epoch | Interval | Defines the interval used for WASO and stage minutes in `calcpsgvars.m` |
| Wake after sleep onset (WASO) | Wake epochs within the sleep period × 0.5 | Minutes | Excludes wake before first sleep and after last sleep |
| Snooze | Time from the last scored sleep epoch to recording end | Minutes | Exported but not a standard focal manuscript metric; keep distinct from WASO |
| Sleep efficiency | `100 × TST / recording opportunity` | Percent | Formula reproduced for all 38 rows |
| Stage duration | N1/N2/N3/REM epochs within sleep period × 0.5 | Minutes | Each stage reported separately |
| Stage percentage | `100 × stage minutes / TST` | Percent of TST | All four formulas reproduced for all 38 rows |
| Whole-night first-REM latency | First REM epoch minus first scored sleep epoch | Minutes | Export field `remlat`; not the same as per-cycle transition latency |
| Overall conventional arousal index | `60 × conventional arousal count / TST minutes` | Events/hour TST | Formula reproduced for all rows |
| NREM conventional arousal index | `60 × NREM conventional arousal count / (N1+N2+N3 minutes)` | Events/hour NREM | Formula reproduced for all rows |
| Stage-specific conventional arousal index | `60 × stage-specific conventional arousal count / stage minutes` | Events/hour in stage | N1, N2, N3, and REM formulas reproduced for all rows |

### ISFσ signal and circular-statistics terms

| Term | Operational definition | Unit/granularity | Interpretation limit |
|---|---|---|---|
| Instantaneous sigma power | Time-resolved 11–16-Hz EEG power derived with the specified wavelet | dB-normalized power/time/channel | Observable EEG quantity, not ISFσ phase or direct NA concentration |
| ISFσ signal | Sigma-power time series bandpass filtered in the infraslow range used by the analysis | Filtered dB time series | Oscillatory EEG proxy; indirect relation to LC–NA |
| ISFσ spectral mean frequency | Center (`mu`) of the fitted infraslow spectral peak | Hz/channel/bout | Global/bout spectral property, distinct from event phase |
| ISFσ spectral peak amplitude | Height of the fitted infraslow spectral peak | Normalized or absolute spectral units | Distinct from Hilbert amplitude and trough-to-peak REM-wave amplitude |
| ISFσ spectral bandwidth | Width of the fitted infraslow spectral peak | Hz/channel/bout | Describes spectral spread, not phase duration |
| Instantaneous Hilbert phase | Angle of the analytic filtered ISFσ signal | Radians/degrees/channel/time | Raw convention: peak 0°, descending crossing 90°, trough 180°, ascending crossing 270° |
| Instantaneous Hilbert amplitude | Modulus of the analytic filtered ISFσ signal | Filtered-signal amplitude/channel/time | Event-local envelope magnitude, not spectral peak amplitude |
| Event-level phase | Circular mean of channel-level phase within the selected regional channel set for one event | One angle/event | Event is the circular-statistics observation; channels are not independent observations |
| Empirically corrected phase | Raw phase transformed through a fixed empirical background CDF and mapped back to a circle | Corrected radians/event | Uniform under empirical null; corrected degrees need not equal literal raw waveform landmarks |
| Ascending ISFσ phase | Raw filtered sigma rises from trough toward peak, approximately 180°–360° | Half-cycle | Refers to sigma direction; lower/declining LC–NA is an indirect interpretation, not a measurement |
| Descending ISFσ phase | Raw filtered sigma falls from peak toward trough, approximately 0°–180° | Half-cycle | Refers to sigma direction; rising/higher LC–NA is an indirect interpretation |
| Uniform phase distribution | Event phases are compatible with equal occupancy around the corrected circle | Distribution-level result | Means no evidence of preferred phase, not evidence that phases are biologically irrelevant |
| Non-uniform phase distribution | Event phases depart from equal circular occupancy | Distribution-level result | Supports phase preference but not causation or a unique mechanism |
| Unimodal/multimodal | One versus more than one density mode, evaluated separately from uniformity | Distribution shape | Must not be treated as synonyms for uniform/non-uniform |
| Pre-arousal phase estimate | `phase_d0`, measured nominally 1.5 seconds before target onset to avoid wavelet contamination | Angle/event | Must not be called phase “at onset” |

### REM-transition terms

| Term | Operational definition | Unit/granularity | Notes |
|---|---|---|---|
| REM-centered extraction window | Approximately −300 to +60 seconds around first REM epoch in an identified sleep cycle | Seconds/transition | Supports filtering, feature detection, and visualization |
| Pointwise pre-REM test window | Samples from −90 seconds through immediately before REM onset | Seconds/transition | Pragmatic figure/multiplicity interval, not a unique biological boundary |
| NREM fluctuation | Earlier trough-to-subsequent-peak ISFσ half-wave before the final transition fluctuation | Half-wave/channel, then averaged | Multiple NREM fluctuations may occur per transition |
| Transition-to-REM fluctuation | Final eligible trough-to-subsequent-peak half-wave before REM; trough selected at least 20 seconds before REM onset | One half-wave/channel/transition | Code labels stage `rem`, but this denotes transition-associated feature rather than a fluctuation wholly within REM |
| REM-wave amplitude | Filtered sigma value at peak minus value at preceding trough | dB/channel, averaged across channels | Distinct from spectral and Hilbert amplitude |
| REM-wave half-period | Time from trough to subsequent peak | Seconds/channel, averaged across channels | Code variable `duration`; it is a half-period, not a full oscillatory cycle |
| Per-cycle REM-transition latency | Time from the algorithm-defined cycle start/first preceding non-wake epoch to first REM epoch | Minutes/transition | Code field `remeplat`; distinct from Table 1 whole-night REM latency |
| REM episode duration | Duration of the REM episode associated with the transition | Minutes/transition | Code field `remepdur` |

## Implemented clinical-metric validation record

All formulas that can be reconstructed from the exported columns were checked across all 38 rows. Maximum absolute discrepancies were approximately `4.3×10⁻14` for sleep efficiency, `5.0×10⁻14` for stage percentages, `4.6×10⁻14` for the overall conventional arousal index, `5.0×10⁻14` for NREM/stage NREM indices, and `3.3×10⁻13` for the REM arousal index. These are floating-point rounding differences; the formulas reproduce the export.

**Validation decision:** the displayed formulas are verified for the current export. `calcpsgvars.m` and `calcarovars.m` remain candidate generator functions because their active call sites are unconfirmed, but their formulas match the saved output exactly.

## Implemented crossover and timing inventory

### Available and unavailable fields

| Field | Availability/source | Derivation/interpretation |
|---|---|---|
| Participant and treatment | Complete; scoring workbook and analytical tables | One row per participant/condition |
| Period and sequence | Complete; chronological order of `recdate` joined to treatment | Period 1/2 and one of two sequences |
| Recording start/effective lights off | Complete; `phenotype/recdates.csv` | Observed visit timestamp |
| Inter-visit washout | Complete | Calendar-day difference between observed recording dates |
| SOL and first-REM latency | Complete; PSG export | Minutes from recording start and first sleep, respectively |
| Derived sleep-onset time | Available | Recording start + SOL |
| Derived first-REM time | Available | Recording start + SOL + whole-night REM latency |
| Actual dose-administration time | Unavailable | Do not impute as an observed timestamp |
| Nominal scheduled dose context | Available from protocol | One hour before typical, not actual, sleep onset |
| Drug concentration/exposure | Unavailable | Time-of-night patterns cannot be translated into measured pharmacokinetics |

### Participant-level sequence and washout

| Participant | Sequence | Inter-visit days |
|---|---|---:|
| r002 | placebo→THC/CBD | 7 |
| r003 | THC/CBD→placebo | 29 |
| r004 | placebo→THC/CBD | 8 |
| r005 | THC/CBD→placebo | 16 |
| r006 | THC/CBD→placebo | 8 |
| r007 | placebo→THC/CBD | 24 |
| r008 | THC/CBD→placebo | 7 |
| r009 | THC/CBD→placebo | 52 |
| r010 | THC/CBD→placebo | 37 |
| r011 | placebo→THC/CBD | 7 |
| r012 | THC/CBD→placebo | 9 |
| r013 | THC/CBD→placebo | 63 |
| r014 | THC/CBD→placebo | 8 |
| r015 | THC/CBD→placebo | 8 |
| r016 | placebo→THC/CBD | 35 |
| r017 | placebo→THC/CBD | 144 |
| r018 | placebo→THC/CBD | 7 |
| r019 | placebo→THC/CBD | 170 |
| r020 | THC/CBD→placebo | 118 |

All participants meet the minimum one-week interval. The broad washout range should be retained as descriptive metadata; it does not itself demonstrate or exclude carryover.

### Phase 2 timing-axis feasibility

1. **Observed axes:** time since recording start/lights off, time since scored sleep onset, clock time, and sleep-cycle number.
2. **Nominal-only axis:** time since scheduled dose, because actual dose times are unavailable.
3. **Preferred primary axis:** select one observed axis in Phase 2 based on the scientific question and event support; use at most one complementary sensitivity axis to limit multiplicity.
4. **Interpretive limit:** an early/late pattern can be temporally compatible with acute drug action but cannot establish concentration–effect decay.

## Implemented published subjective-outcome inventory

| Construct/instrument | Protocol assessment time | Published result available to this revision | Applicability/limit |
|---|---|---|---|
| Perceived sleep change/quality — Leeds Sleep Evaluation Questionnaire (LSEQ) | Approximately 07:30 after the treatment night | Parent results report no change in subjective sleep quality | Parent n=20 only; participant-level data unavailable |
| Daytime insomnia symptoms — DISRS | Approximately 07:30 | Instrument was planned; no separate result needed for the reviewer's perceived-sleep question has yet been verified | Do not substitute for perceived sleep amount |
| Subjective sleepiness — Karolinska Sleepiness Scale (KSS) | Before bed, immediately at awakening, then 10:00, 12:00, 14:00, 16:00 | Small increase with THC/CBD reported in the parent abstract: +0.42 points, p=0.02, d=0.22 | Next-day sleepiness is not subjective sleep quality or perceived sleep duration |
| Objective alertness — KDT/MWT-related measures | Before bed/awakening and next-day schedule as specified | Parent abstract reports no observed next-day objective-alertness change | Distinct from self-report |
| Subjective drug effects — DEQ | Baseline, 60 minutes after administration, and approximately 08:00 | Planned context measure | Does not directly answer whether participants thought they slept more |

**Reviewer-answer boundary:** the available publication supports a direct answer about subjective sleep quality and next-day sleepiness. Whether the LSEQ result separately quantifies perceived sleep amount remains a source-text check. If not separately reported, the response must say that perceived amount cannot be answered from accessible participant-level data.

### Exact published LSEQ results

The full parent-results paper reports 40/40 treatment-session observations and no treatment difference in any LSEQ domain:

| LSEQ domain, 0–100 mm (higher is better) | Placebo mean (SD) | THC/CBD mean (SD) | p | Cohen's d [95% CI] |
|---|---:|---:|---:|---:|
| Getting to sleep | 50.5 (19.1) | 48.0 (12.7) | 0.620 | −0.11 [−0.55, 0.33] |
| Quality of sleep | 43.2 (22.4) | 41.1 (25.7) | 0.777 | −0.06 [−0.50, 0.34] |
| Awakening from sleep | 39.7 (19.1) | 34.7 (18.9) | 0.400 | −0.19 [−0.64, 0.25] |
| Behavior following wakefulness | 28.8 (15.9) | 29.3 (14.9) | 0.924 | 0.02 [−0.42, 0.46] |

No direct subjective total-sleep-time or perceived “slept more” outcome is reported. The four LSEQ domains must not be reinterpreted as perceived sleep duration. The defensible later answer is that participants did not report better sleep–wake behavior on any measured LSEQ domain, while whether they perceived greater sleep quantity was not directly measured/reported in an accessible form.

## Implemented dose/formulation fact sheet

| Item | Confirmed fact | Source/limit |
|---|---|---|
| Product | ETC120, plant-derived oral cannabinoid formulation | Parent protocol |
| Active dose | 10 mg Δ9-THC + 200 mg CBD in 2 mL | Parent protocol/results |
| Ratio/vehicle | THC:CBD 1:20 in medium-chain triglyceride oil | Parent protocol |
| Placebo | Matched 2 mL containing no cannabinoids | Parent protocol |
| Administration | Single directly observed oral dose scheduled one hour before typical sleep onset | Actual dose timestamps unavailable |
| Protocol THC rationale | Chosen as the maximum dose expected to produce sleepy/tired subjective effects without substantial cognitive impairment or intoxication in naïve/occasional users | Context-specific rationale, not a universal “high-dose” category |
| Protocol CBD rationale | Included for possible anxiolytic and adverse-effect-mitigating properties | Rationale, not evidence that CBD caused the observed sleep effects |
| Design limitation | Fixed combination with no THC-only or CBD-only arm | Cannot identify THC-specific, CBD-specific, additive, or interactive causation |
| Exposure information | Plasma cannabinoids/metabolites were measured approximately 11 hours post-treatment and published by treatment order | Late single-time-point context only; not an overnight concentration–time profile and participant-level values are unavailable here |

There are no universal dosing guidelines that make “high” or “low” a scientifically sufficient classification across route, formulation, population, tolerance, exposure, and outcome. Phase 3 should answer the reviewer's underlying concern by contextualizing the acute oral doses and transparently stating where human LC–NA evidence for CBD or the exact 1:20 mixture is absent or indirect. It should not promise that every requested high/low categorization has a definitive answer.

## Implemented REM variable and window record

1. **Whole-night REM latency** (`remlat`) is the first REM epoch relative to first scored sleep and is reported in Table 1.
2. **Per-cycle REM-transition latency** (`remeplat`) is the first REM epoch relative to the cycle's algorithm-defined start/first preceding non-wake epoch and is repeated across transitions.
3. The manuscript reports 118 transitions: 51 THC/CBD and 67 placebo.
4. The final transition-associated trough is selected at least 20 seconds before REM onset; amplitude and trough-to-peak time are then extracted.
5. The pointwise model indexes only −90 to <0 seconds. Therefore the pre-REM visual dip is within the tested family and was not significant after correction, while the post-REM visual rise lies outside that model.
6. Phase 2 should add the planned post-REM comparison as an explicitly reviewer-requested/exploratory test and define multiplicity control over the tested interval before execution.

## Implemented arousal-classification decision table

| Onset situation | Assigned stage for analysis | Outcome comparison | Awakening-associated when | Key limitation |
|---|---|---|---|---|
| Onset in an epoch already scored N2, N3, or REM | Current epoch stage | Stage of the following 30-second epoch | Following epoch is wake or N1 | Exact physiological wake/N1 onset within relation to the epoch label is not continuously annotated |
| Onset in an epoch scored N1 | Current epoch N1 | Stage of the following epoch | Following epoch is wake | N1 is treated as sleep for this rule but as an awakening outcome after deeper NREM/REM |
| Onset within first 15 seconds of a wake/N1 epoch preceded by the requested sleep stage | Preceding epoch's requested stage | Stage of the epoch following the onset epoch | Following epoch satisfies the stage-specific rule above | The current wake/N1 epoch motivates reassignment but is not itself the `next_stage` field used for outcome classification |
| Onset after first 15 seconds of a wake/N1 epoch preceded by the requested stage | Not assigned to that requested-stage analysis | Not selected for that stage | Not applicable | Excluded from that stage-specific extraction |
| Onset without at least 10 seconds of preceding stable sleep | Not a valid target under the applied arousal rule | Not selected | Not applicable | Prevents immediate repeat events from entering as targets |

### Supported temporal information

- `onset_epoch` and continuous event latency support calculation of position within the 30-second epoch and time to the next epoch boundary.
- `next_stage` supports only an epoch-level outcome label.
- The project does not currently contain a separate continuously scored physiological wake-onset annotation that would support an exact arousal-to-wake latency.
- Phase 2 can stratify or adjust by target position within the epoch and report the corresponding boundary interval. It must not present that boundary interval as exact physiological wake latency.
- Shorter-timescale rescoring would require a new prespecified rule, blinded review, and validation; it is the strongest but most resource-intensive option rather than a prerequisite for the minimum clarification.

## Implemented >30-second duration audit

- Fifteen of 6,090 conventional arousals (0.25%) exceed 30 seconds; range 30.026–33.040 seconds.
- All 15 are `arousalemg` events. They occur in 10 participants, across N1 (2), N2 (6), N3 (2), and REM (5).
- Four occur under placebo and 11 under THC/CBD; 14 are awakening-associated and one is continued-sleep.
- Six are in the conventional N2 pool relevant to the broader N2 inventory.
- No event exceeds 33.04 seconds, so the cap affects a small tail rather than grossly miscoded multi-minute events.
- The cap is an accepted robustness treatment for the awakening model, but its use and the raw observed range should be disclosed. Phase 2 should decide whether an uncapped/winsorization sensitivity is needed; Phase 1 does not infer a treatment effect from the descriptive imbalance in these 15 events.

## Implemented corrected-phase landmark validation

The fixed empirical background PDF integrates to 1.0. Applying the implemented probability-integral transform to canonical raw landmarks gives:

| Raw landmark | Raw angle | Corrected angle | Shift |
|---|---:|---:|---:|
| Trough | −180°/180° | −180°/180° | 0° |
| Ascending zero crossing | −90° (270°) | −87.29° (272.71°) | +2.71° |
| Peak | 0° | −0.84° (359.16°) | −0.84° |
| Descending zero crossing | 90° | 86.64° | −3.36° |

**Decision:** qualitative landmark labels remain reasonable approximations for the current corrected display because shifts are under 3.4°. Any numerical peak angle must nevertheless be called an empirically corrected phase angle; exact corrected degrees must not be represented as untransformed Hilbert degrees.

## Healthy-control eligibility-gate status

**Current decision: unresolved/not approved for analysis.** Technical presence is established, but source-study eligibility and scientific comparability are not.

| Domain | Current evidence | Status required before analysis |
|---|---|---|
| Source study, consent, and permitted use | Not documented in the reviewed Phase 1 sources | Required |
| Healthy screening/diagnostic definition | Group label only | Required |
| Hardware, montage, channels, sampling | Derivatives exist; equivalence not formally documented | Required |
| Recording opportunity/laboratory procedure | Not established | Required |
| Scoring rules/raters | Not established | Required |
| Preprocessing compatibility | Shared pipeline appears technically possible | Must verify source-level comparability |
| Age/sex comparability | Lightweight participant information exists; matching file is present but not approved as an analysis definition | Must prespecify matching/harmonization |
| Medication/substance restrictions | Not established | Required |
| Adaptation/first-night status | Not established | Required |
| Available PSG/ISFσ variables | Processed derivatives for 18 controls indicate technical availability | Must inventory completeness/QC |
| Missingness/QC | 22 labeled versus 18 processed/timestamped | Must reconcile |
| Permissible inference | Untreated observational comparison only | Must explicitly prohibit diagnosis-by-treatment inference |

**Fallback if documentation is not supplied:** no healthy-control analysis in Phase 2; answer reviewers through transparent limitation and claim moderation. Healthy-control work must not delay unrelated randomized-trial responses.

## Implemented Phase 2 handoff matrix

| Reviewer need | Population/input | Preferred specification | Minimum adequate fallback | Key risk/validation |
|---|---|---|---|---|
| Variable arousal-to-awakening timing | Conventional stage-specific events; onset position, epoch boundary, `next_stage` | Quantify target position/time to boundary and test robustness to epoch position; consider shorter rescoring only if feasible | Exact operational rule plus limitation | Boundary is not physiological wake onset; validate worked cases |
| Fully isolated arousals | N2 phase events with event history | Re-run key phase comparison with no event in prior 60 seconds versus primary no-event-in-10-seconds rule | Report primary 10-second rule and modeled preceding-event probability | Event loss/power and selection toward stable sleep; report attrition |
| Phase-offset sensitivity | Saved/derivable offsets at −1.5, −2.5, −3.5, −4.5, −5.5, −9.5, −13.5, −17.5 s | Prespecify a small near/intermediate/distal set before testing | Rationale plus existing descriptive offset display | Multiplicity/favorable-offset selection; use consistent correction |
| Treatment order/period | 19 participants, 38 sessions, sequence and period table | Parsimonious treatment model adjusted for sequence/order and period as appropriate | Sequence-adjusted sensitivity | n=19 cannot prove no carryover; verify coding and estimand |
| Evolution across night | Event/transition timing relative to recording start, sleep onset, or cycle | Choose one primary observed axis and one sensitivity axis at most | Descriptive temporal context | No actual dose time/concentration; avoid PK claims |
| >30-second arousals | Fifteen conventional events; model-capped duration | If needed, compare capped primary with uncapped or exclusion sensitivity | Disclose raw range and cap | Very small tail; avoid overemphasis; verify transformation in executed model |
| Pre-REM visual dip | 118 transitions; −90 to <0 pointwise window | Report current corrected family result; no additional selective local test unless prespecified | State no corrected difference in tested window | Avoid post hoc sub-window significance fishing |
| Post-REM visual rise | +30 to +60 seconds within extracted data but outside current pointwise model | Reviewer-requested exploratory comparison over prespecified post-REM interval with multiplicity control | State it was not tested and avoid visual inference | REM sigma physiology differs from pre-REM; define model/interval before execution |
| Healthy comparison | Candidate 22/processed 18 controls | Proceed only after eligibility gate and matching/harmonization plan | Limitation and claim moderation | Observational controls cannot test diagnosis-specific drug response |

No handoff authorizes analysis until Phase 2 planning specifies estimands, models, multiplicity, missingness, and outputs.

## Implemented dual-manuscript revision-location map

| Factual content to hand off | IMRD placement | IRDM placement | Synchronization rule |
|---|---|---|---|
| Parent/current sample and R001 exclusion | Abstract and Experimental Design | Abstract; concise Results context if needed; full later Experimental Design | Identical n/reason; order-only difference |
| Formulation and scheduled timing | Experimental Design | Brief first-use Results context only if essential; full Experimental Design | Same facts; do not imply actual dose time |
| Sequence, period, washout | Experimental Design/Statistical analysis; new result in Results after Phase 2 | Concise design context before order result; full later Methods | Same analysis/result/qualification |
| Clinical metric formulas/denominators | First Table 1 use and Methods | Define essential metrics at first Results/Table use; full later Methods | Same formulas; IRDM may use concise parenthetical definitions |
| Arousal stage/outcome decision rule | Sleep scoring/stability Methods | One concise operational sentence before first awakening result; full later Methods | Same rule, reassignment, and causal limitation |
| 10-second versus 10–60-second rule | Phase Methods and sensitivity Results after Phase 2 | Concise primary-rule context before sensitivity result; full later Methods | Same selection and attrition |
| −1.5-second phase timing | Phase Methods; avoid “at onset” in Results/legends | Concise pre-arousal timing at first phase result; full later Methods | Same offset and contamination rationale |
| Raw/corrected phase and intuitive landmarks | Methods plus main-text intuitive explanation; figure legend | Intuitive explanation before result; technical correction later Methods | Same corrected/raw distinction |
| NREM versus transition fluctuation and two REM latencies | REM Methods and first relevant Results paragraph | Concise definitions before REM Results; full later Methods | Use unique names consistently |
| 90-second rationale and post-REM test status | REM Methods/Results/figure legend as appropriate | Results may state tested interval; rationale in later Methods | Same exploratory label and multiplicity result |
| Published subjective findings | Discussion or clinical-context Results only if journal structure supports; source as parent analysis | Corresponding section in reordered manuscript | Explicit parent n=20; no new n=19 analysis |
| Healthy-control limitation/gate outcome | Discussion limitations and response | Same Discussion text | Identical inferential boundary |
| Dose high/low uncertainty | Introduction/Discussion after Phase 3 research | Same scientific text | No version-specific pharmacological claim |

The IRDM exception is limited to concise context required to make Results intelligible before Methods; it must not create a different scientific account or duplicate full technical detail.

## Ordered work breakdown

### 1. Reviewer traceability matrix
- **Implementation status:** Complete; see implemented matrix above.
- **Purpose/concern:** Cover every relevant sub-question without collapsing compound comments.
- **Inputs/action:** Comments, manuscripts, plans; assign stable IDs and map Phase 1 facts, response type, later action, and output.
- **Output/dependency:** Traceability matrix; no dependency.
- **Risk/validation:** Omission/duplication; second pass through every comment line.
- **Acceptance:** Every relevant sub-question has one primary entry and appropriate cross-links.

### 2. Source-of-truth register
- **Implementation status:** Complete for current Phase 1 sources; update if new sources are added.
- **Purpose:** Resolve protocol/data/code/manuscript/history conflicts.
- **Inputs/action:** Source hierarchy; record path, version, authority, scope, limitation, linkage, and conflict.
- **Output/dependency:** Source register; depends on Task 1.
- **Risk/validation:** Current versus executed code; cross-check saved outputs and obtain confirmation where needed.
- **Acceptance:** Every inventory field has a source or unresolved flag.

### 3. Analysis populations
- **Implementation status:** Complete for currently identified analysis families; REM compact-table provenance remains an open source item.
- **Purpose/concern:** Separate all samples; addresses n discrepancy and control feasibility.
- **Inputs/action:** Parent papers, inputs, log, groups; enumerate inclusion/exclusion by analysis.
- **Output/dependency:** Population register/flow; Task 2.
- **Risk/validation:** Assuming identical samples; reproduce counts from inputs.
- **Acceptance:** Every analysis n/exclusion is reproducible.

### 4. Arousal/awakening rules
- **Implementation status:** Complete for the current operational rule and supported timing; Phase 2 rescoring feasibility remains a handoff.
- **Purpose/concern:** Address variable interval and causal precision.
- **Inputs/action:** Scoring rules, code, events/hypnograms; formalize all paths, bounds, and edge cases.
- **Output/dependency:** Decision table/temporal diagram; Tasks 2–3.
- **Risk/validation:** No exact within-epoch onset; worked epoch-start/middle/end examples.
- **Acceptance:** Unambiguous paths and explicit limitation.

### 5. Duration and neighboring events
- **Implementation status:** Complete for inventory and >30-second audit; formal standard wording/cap disclosure remains an editorial handoff.
- **Purpose/concern:** Resolve maximum and isolation questions.
- **Inputs/action:** Exports, code, standard; separate classes/units, audit >30 seconds, document cap and 0–10/10–60 rules.
- **Output/dependency:** Duration/isolation specification; Task 4.
- **Risk/validation:** Samples versus seconds or wrong event class; reproduce counts and inspect events.
- **Acceptance:** Formal, raw, transformed, analyzed, and neighboring-event rules are distinct.

### 6. ISFσ glossary/phase specification
- **Implementation status:** Complete, including corrected-landmark validation.
- **Purpose/concern:** Separate signal measurements from mechanism.
- **Inputs/action:** Methods/code/correction/plots; define features, landmarks, correction, channel averaging, offsets, distribution terms, caveats.
- **Output/dependency:** Glossary and Phase 5 schematic brief; Task 2.
- **Risk/validation:** Raw landmarks assigned to corrected angles; map a known trace through raw/corrected paths.
- **Acceptance:** One definition/unit/source/limit per variable; corrected labels validated or qualified.

### 7. Clinical metrics
- **Implementation status:** Complete; all reconstructable formulas reproduce the current export.
- **Purpose/concern:** Verify formulas and denominators.
- **Inputs/action:** PSG export and derivation functions; dictionary plus numerical checks.
- **Output/dependency:** Verified dictionary; Task 2.
- **Risk/validation:** Inactive function differs from generator; reproduce rows and resolve mismatches.
- **Acceptance:** Exact formula, denominator, unit, anchor, source, and status for all Table 1 metrics.

### 8. Crossover/timing metadata
- **Implementation status:** Complete for available fields; actual dose timing is confirmed unavailable.
- **Purpose/concern:** Supply order, period, washout, and night-time inputs.
- **Inputs/action:** Workbook, recdates, PSG timing; join and derive period, sequence, washout, sleep/REM clock time.
- **Output/dependency:** Analysis-ready metadata specification/table without models; Tasks 2, 3, 7.
- **Risk/validation:** Timezone/date rollover and nominal-dose confusion; one-to-one, chronology, washout, and spot checks.
- **Acceptance:** All 38 sessions link; observed versus nominal timing is explicit.

### 9. n=20→19 fact sheet
- **Implementation status:** Complete.
- **Purpose/concern:** Auditable reviewer answer.
- **Inputs/action:** Parent flow, R001 record, current tables; document completion, one-night failure, paired requirement.
- **Output/dependency:** Population-flow fact sheet; Task 3.
- **Risk/validation:** Calling it withdrawal/outcome exclusion; verify R001 absent and 19 pairs present.
- **Acceptance:** Consistent accurate explanation for later outputs.

### 10. Published subjective outcomes
- **Implementation status:** Complete using the full parent-results Methods, Results, and Table 3; no direct perceived sleep-amount outcome was reported.
- **Purpose/concern:** Answer perceived-sleep questions without participant data.
- **Inputs/action:** Parent papers; record measures, timing, estimates, n, and perceived-amount availability.
- **Output/dependency:** Subjective record; Task 2.
- **Risk/validation:** Conflating quality/quantity/sleepiness; exact page/table citations.
- **Acceptance:** Source-cited n=20 statements with no implied n=19 reanalysis.

### 11. Dose/formulation and research brief
- **Implementation status:** Fact sheet complete; literature interpretation remains assigned to Phase 3.
- **Purpose/concern:** Avoid simplistic exposure and constituent claims.
- **Inputs/action:** Parent papers/manuscript; facts plus bounded Phase 3 questions.
- **Output/dependency:** Fact sheet/research handoff; Task 2.
- **Risk/validation:** Universal high/low or constituent attribution; match protocol.
- **Acceptance:** Supported facts; interpretation deferred.

### 12. Healthy-control gate
- **Implementation status:** Applied; current outcome is unresolved/not approved, with documentation requirements and fallback recorded.
- **Purpose/concern:** Decide defensibility of requested comparison.
- **Inputs/action:** Source documentation and all eligibility domains; assess comparability/QC/inference.
- **Output/dependency:** Eligibility decision and, only if eligible, bounded Phase 2 brief; documentation availability.
- **Risk/validation:** Controls treated as randomized placebo; team-reviewed comparability table.
- **Acceptance:** No analysis without eligible/conditional status and limits.

### 13. REM variables/windows
- **Implementation status:** Complete for current variables/windows and Phase 2 post-REM handoff.
- **Purpose/concern:** Distinguish REM constructs and tested intervals.
- **Inputs/action:** REM code/manuscript; define windows, fluctuation classes, two latency measures, post-REM availability.
- **Output/dependency:** REM specification/Phase 2 handoff; Tasks 2, 6, 7.
- **Risk/validation:** Retroactive biological rationale; code-to-text/axis trace.
- **Acceptance:** Unique names and explicit tested/untested intervals.

### 14. Phase 2 handoff matrix
- **Implementation status:** Complete at input/feasibility level; Phase 2 must still plan models before execution.
- **Purpose:** Enable analyses without redefining inputs.
- **Inputs/action:** Tasks 3–13; specify concern, population, variables, timing, feasibility, minimum/strongest response, fallback, risk, validation.
- **Output/dependency:** Analysis-input matrix; preceding inventories.
- **Risk/validation:** Drifting into execution; one-to-one reviewer mapping.
- **Acceptance:** Phase 2 can start without input ambiguity.

### 15. Dual-manuscript location map
- **Implementation status:** Complete at factual-location level; no manuscript text has been edited.
- **Purpose:** Prepare synchronized IMRD/IRDM revision without editing.
- **Inputs/action:** Both manuscripts and completed inventories; paired locations and IRDM context requirements.
- **Output/dependency:** Location map; completed inventory.
- **Risk/validation:** Scientific divergence or excessive IRDM duplication; side-by-side checklist.
- **Acceptance:** Paired locations or justified structural exception for every factual revision.

## Phase 2 handoffs

Provide inputs—but no results—for: awakening timing/epoch position; fully isolated arousals; phase offsets; sequence/period; observed-axis time of night; >30-second duration sensitivity if needed; pre-/post-REM comparisons; and healthy controls only after eligibility. Each handoff must state the minimum response, strongest feasible response, fallback, validation, and overreaching approaches.

## Risks and mitigation

| Risk | Mitigation |
|---|---|
| Manuscript treated as evidence | Apply source hierarchy/conflict register |
| Current code differs from executed code | Trace outputs/history and confirm provenance |
| Epoch labels treated as exact awakening | Use epoch-associated terms and supported bounds only |
| Corrected phase treated as raw radians | Separate conventions and validate labels |
| Amplitude constructs conflated | Unique names, units, derivations |
| Nominal dose timing treated as observed | Separate timing fields and pharmacokinetic limits |
| Sequence/period/carryover conflated | Define separately; parsimonious Phase 2 plan |
| Parent n=20 presented as n=19 | Label source population for every outcome |
| Controls treated as randomized placebo | Enforce gate and separate streams |
| Historical QC note treated as current failure | Record resolution; R008 is corrected |
| Data-driven sensitivity proliferation | Prespecify a small motivated set |
| Scope drift | Limit Phase 1 to definitions, inventory, feasibility, handoff |

## Deliverables

1. Reviewer traceability matrix and source/provenance register.
2. Operational glossary and analysis-population register.
3. Arousal timing/duration/neighbor specification and ISFσ phase/feature specification.
4. Verified clinical metric dictionary and crossover/timing inventory.
5. n=20→19 fact sheet, published subjective record, and dose/formulation fact sheet.
6. Healthy-control eligibility decision or unresolved status.
7. REM variable/window specification.
8. Phase 2 input/feasibility matrix and later-phase handoffs.
9. IMRD/IRDM location map.
10. Maintained decision, open-item, and progress logs.

## Open non-blocking checks

These items have documented outcomes and therefore do not prevent Phase 1 completion. They remain traceable constraints or later-phase contingencies.

| Item | Phase 1 outcome | Later action/fallback |
|---|---|---|
| Formal arousal maximum | AASM definition supplies a minimum but no maximum; the project retained raw events above 30 seconds but applied a 30-second analytical cap in the awakening model | Report both facts and the observed raw range; do not call the model cap a scoring maximum |
| Fifteen conventional events >30 seconds | Audited: 30.026–33.040 seconds, all `arousalemg`, 10 participants, 14 awakening-associated | Disclose cap; Phase 2 may run robustness handling if needed |
| Corrected phase landmarks | Validated; maximum canonical-landmark shift is 3.36° | Qualitative labels are acceptable approximations; numerical angles remain labeled corrected |
| Perceived sleep amount | Full parent paper reports four LSEQ domains but no subjective TST/“slept more” measure | State that measured LSEQ domains were null and perceived quantity was not directly reported |
| Control documentation | Not supplied/currently unresolved | Gate outcome is not approved; limitation-only response unless documentation is later added |
| 22 labeled vs 18 processed controls | Unreconciled because control stream is not approved | If reactivated, reconcile before defining an analysis sample |
| Rationale for duration cap | Git history records only “Cap extreme durations”; no deeper contemporaneous rationale was found | Describe transparently as project analysis handling; do not invent an a priori rationale |
| REM compact event table | Manuscript and current workflow establish n=118 and feature definitions; no compact saved table identified | Phase 2 must create/validate its analysis table from the workflow before any new REM test |

## Phase 1 completion validation

| Acceptance domain | Evidence in this file | Result |
|---|---|---|
| Reviewer coverage | Stable traceability matrix covers all Reviewer 1 major/minor and Reviewer 2 major/numbered comments | Pass |
| Scientific/methodological defensibility | Epoch limits, indirect LC–NA mapping, fixed-combination limits, and observational-control limits are explicit | Pass |
| Definitions and units | Arousal, clinical sleep, ISFσ, circular, and REM glossaries completed | Pass |
| Population clarity | Parent, current, event-class, phase-subset, REM-transition, subjective, and control populations separated | Pass |
| Data provenance | Source hierarchy/register, conflicts, fingerprints, and unresolved-source fallbacks documented | Pass |
| Crossover/timing | All 38 sessions link; sequence and 7–170-day washout inventoried; nominal versus observed timing separated | Pass |
| Metric validation | All reconstructable PSG/arousal-index formulas reproduce the export to floating-point precision | Pass |
| Arousal reviewer concerns | Decision rule, temporal support, neighboring events, raw/capped duration, and >30-second audit completed | Pass |
| Phase reviewer concerns | Offset map, raw/corrected convention, channel/event unit, and corrected landmarks validated | Pass |
| Subjective outcomes | Full parent Table 3 recorded; n=20 and perceived-amount limitation explicit | Pass |
| Dose/formulation | Factual preparation and inferential limits recorded; research interpretation deferred | Pass |
| Healthy controls | Eligibility gate applied; current not-approved status and fallback explicit | Pass |
| REM definitions/windows | Two latency variables, feature definitions, tested/untested windows, and post-REM handoff recorded | Pass |
| Manuscript equivalence | Paired IMRD/IRDM location map and IRDM concise-context exception documented | Pass |
| Later-phase readiness | Phase 2 handoff matrix plus Phase 3–7 dependencies documented | Pass |
| Scope control | No treatment model, manuscript/figure edit, or response drafting performed | Pass |

## Phase-level acceptance criteria

Phase 1 is complete only when:

1. Every relevant reviewer sub-question has a unique traceable entry.
2. Every operational term has one definition, unit, source, and interpretation limit.
3. Measurements, transformations, statistics, and biological interpretation are separated.
4. Parent n=20, current n=19, analysis-specific, subjective, and control populations are distinct.
5. R001's one-night failure and paired-analysis exclusion are consistently documented.
6. Arousal classification includes epoch timing, 15-second reassignment, and causal limits.
7. The 0–10-second exclusion and retained 10–60-second events are documented.
8. Formal/raw/transformed/analyzed durations and the cap are addressed.
9. `phase_d0` and all offset fields are correctly mapped.
10. Raw/corrected phase conventions and corrected labels are validated or qualified.
11. All ISFσ frequency, amplitude, phase, and period constructs are uniquely defined.
12. Clinical metrics and indices have verified formulas, denominators, units, and anchors.
13. Whole-night and per-cycle REM latency are not conflated.
14. All 38 sessions link uniquely across participant, treatment, period, sequence, and recording start.
15. Washout is verified; observed and nominal timing are separate.
16. Subjective findings specify construct, timing, source, and n=20 limitation.
17. Dose facts avoid universal high/low and constituent-specific attribution.
18. Controls pass the gate or remain explicitly unavailable.
19. Every Phase 2 analysis has inputs, feasibility, fallback, risks, and validation.
20. Both manuscripts have synchronized location mapping and IRDM context requirements.
21. Every unresolved fact has a source/action and fallback.
22. Null, contradictory, and unavailable information is handled transparently.
23. No scientific analysis, manuscript/figure revision, or response drafting has been performed as Phase 1 planning work.

## Dependencies and handoff

- Phase 2 depends on populations, timing, arousal rules, offsets, REM windows, and the control gate.
- Phase 3 depends on measurement/interpretation boundaries, dose facts, global null findings, subjective outcomes, and control limits.
- Phase 4 depends on glossary, metric dictionary, factual clarifications, and location map.
- Phase 5 depends on phase conventions/schematic specification and Phase 2 results.
- Phase 6 depends on paired locations and identical scientific definitions.
- Phase 7 depends on traceability and documented evidence for every answer/action.

## Change control

During Phase 1 implementation, update this file with evidence, decisions, progress, unresolved items, and completion status. Do not implement analyses or manuscript, figure, or response revisions unless explicitly authorized in the appropriate later phase.