# Phase 2 — Reviewer-requested analyses

## Status and progress

| Field | Status |
|---|---|
| Planning | Protocol v0.2 narrowed to four selected reviewer concerns; decisions below are locked unless a documented feasibility check requires amendment |
| Implementation | Not started |
| Results inspected under this protocol | None |
| Last planning review | 2026-06-11 |
| Current authorized scope | Phase 2 specification, implementation, validation, and reporting |

- [x] Read the overarching workplan, Phase 1 inventories/handoffs, reviewer comments, manuscripts, and current analysis code.
- [x] Decompose Phase 2 into reviewer questions and distinguish sensitivity, reviewer-requested, and exploratory work.
- [x] Draft estimands, populations, models, timing definitions, multiplicity families, outputs, and validation gates before inspecting new results.
- [ ] Reconstruct and freeze versioned analysis tables from the existing workflow.
- [ ] Complete pre-model data/QC checks and record any protocol amendments.
- [ ] Implement and validate executable Workstreams 1–3; activate Workstream 4 only after ethics and eligibility approval.
- [ ] Prepare machine-readable summaries, supplementary outputs, and conclusions.
- [ ] Hand results to Phases 3–5 and 7 and mark Phase 2 complete.

## Purpose

Plan, conduct, validate, and report the additional and sensitivity analyses needed to answer the reviewers’ methodological and empirical questions. Phase 2 starts from the definitions and feasibility decisions in `workplan/01-definitions-and-data-inventory.md`; it does not redefine those measurements after results are known.

## Scope and boundaries

### Included

1. Quantify the interval from a continuously marked arousal onset to the relevant 30-second epoch boundary and assess whether awakening results depend on onset position within the epoch, as well as treatment condition.
2. Test whether selected treatment effects change over observed time since sleep onset, with time since recording start as the sole complementary timing axis.
3. Adjust principal randomized treatment contrasts for treatment sequence and explore treatment-order heterogeneity.
4. Integrate healthy-control analyses once ethics have been approved (under review).

### Excluded

- Repeat the principal N2 pre-arousal sigma-power and phase analyses after restricting to targets with no preceding arousal in 60 seconds.
- Assess phase conclusions at a small, prespecified set of pre-arousal offsets.
- Report the already tested pre-REM interval and formally test the reviewer-identified post-REM interval.
- Provide a source-bounded synthesis of published subjective outcomes.
- Inferring exact physiological wake onset from 30-second stage labels.
- Calling observed time of night “time since dose,” estimating overnight exposure, or making pharmacokinetic claims; actual dose timestamps and overnight concentrations are unavailable.
- Selectively testing visually favorable sub-windows or offsets after viewing results.
- Reanalyzing participant-level subjective outcomes, which are not available in this repository.
- Treating healthy controls as a randomized placebo arm or testing diagnosis-by-treatment interactions without an eligible design.
- Revising manuscripts, final figures, interpretation, or the response letter; those are downstream phases, although Phase 2 will supply text-ready inputs.

## Reviewer-question traceability

| ID | Reviewer question | Phase 2 response | Classification |
|---|---|---|---|
| R1-M1 | How variable is the interval between cortical arousal and the awakening-associated stage label? | Workstream 1 | Reviewer-requested methodological analysis |
| R1-M6 | Do treatment effects change across the night/after nominal dosing? | Workstream 2 using observed sleep/recording axes only | Reviewer-requested exploratory analysis |
| R2-m3 | Is there an effect of treatment order? | Workstream 3 | Reviewer-requested sensitivity |
| R1-M1-control | Can healthy controls establish how normal and insomnia sleep differ? | Workstream 4, currently gated | Blocked/unavailable unless ethics and eligibility approval are obtained |

Interpretive, restructuring, and figure-design portions of these comments remain assigned to Phases 3–5 and 7.

## Cross-cutting analysis rules

### Analysis populations and observational units

- The randomized analytical cohort contains 19 participants and 38 treatment sessions unless an outcome-specific availability rule removes observations.
- Every output must identify the participant, session, event/transition, treatment, period, treatment sequence, and analysis-specific inclusion reason where those fields apply.
- Arousal models use one conventional arousal as the event-level observation with participant random effects. Channel values are reduced to the prespecified regional event value; channels are not independent observations.
- Time-of-night analyses use the original endpoint-specific observational unit: arousal event, N2 event, NREM bout, or REM transition. Repeated observations remain clustered within participant.
- If phase is included as a time-of-night outcome, circular summaries use one regional circular mean per event and participant-clustered uncertainty/permutation procedures.
- Healthy controls, if approved, form a separate observational cohort and are never treated as randomized treatment sessions.
- Complete cases are defined separately for each estimand. No outcome or timing value will be imputed.

### Treatment and reference coding

- Treatment: placebo = reference; THC/CBD (`etc120`) = comparison.
- Period: first visit = reference; second visit = comparison.
- Sequence: placebo→THC/CBD = reference; THC/CBD→placebo = comparison.
- Awakening outcome: continued sleep = 0; awakening-associated = 1.
- Continuous timing terms are expressed per 60 minutes and centered at the pooled mean after the analysis population is selected. Both raw and centered values must be saved.
- Categorical coding and all reference levels must be asserted in automated pre-model checks, not inferred from default MATLAB ordering.

### Estimation and reporting

- Report the model coefficient/contrast, 95% confidence interval, and exact or bounded multiplicity-adjusted p-value. For binomial models also report odds ratios and model-based absolute probabilities at meaningful reference values.
- Descriptive tables must report participants, sessions, and events/transitions by treatment and relevant subgroup, plus exclusions and missingness.
- A sensitivity finding is **directionally supported** when its treatment-effect direction agrees with the original estimate; **materially supported** additionally requires an interval that does not include the null under the specified analysis. A wider interval caused by event loss is not by itself evidence of contradiction.
- “No evidence of difference” must not be written as equivalence or proof of no carryover/time trend.
- All stochastic analyses use a recorded seed. Production permutation analyses use at least 5,000 valid permutations; Monte Carlo uncertainty for an adjusted p-value must be reported or bounded. Small development runs cannot be reported as final results.

### Multiplicity

Multiplicity is controlled within the explicitly defined families below; there is no omnibus correction across scientifically distinct reviewer questions.

1. **Arousal-timing family:** the treatment-by-epoch-position interaction in the combined original awakening model is the principal test. Stage-specific interactions are secondary and Holm-corrected across eligible stages; the final-third subset is a sensitivity estimate rather than an additional interaction test.
2. **Time-of-night families:** Holm correction across prespecified treatment-by-time interaction tests within each endpoint group described in Workstream 2. Descriptive early/late plots are not additional tests.
3. **Treatment-order families:** no correction for the sequence-adjusted treatment estimate for each original endpoint; treatment-by-sequence heterogeneity tests are exploratory and Holm-corrected within endpoint group.
4. **Healthy-control families:** to be specified in the required dated amendment before any control result is inspected; no control multiplicity plan is authorized while the workstream remains gated.

Any change to these families requires a dated amendment made before the affected result is inspected.

## Required analysis-table build

Before fitting models, create immutable analysis tables under `group-level/` (gitignored data/output area) and tracked build/analysis code under `code/group-level/`.

| Table | Minimum content | Primary source |
|---|---|---|
| `sessions.csv` | participant, treatment, recording start, sleep onset, period, sequence, washout days | `phenotype/recdates.csv`, PSG export, Phase 1 crossover derivation |
| `arousal_events.csv` | stable event ID, participant/session, onset, onset epoch, position in epoch, time to relevant boundary, assigned/current/next stage, awakening label, raw/capped duration, preceding-event intervals, sleep/recording elapsed time | arousal export and segmented-event workflow |
| `time_of_night_events.csv` / `.mat` | stable event/transition ID, participant/session, treatment, elapsed sleep and recording time, endpoint values, and inclusion flags required by Workstream 2 | validated Aim 1–3 workflows |
| `healthy_controls.csv` / `.mat` | eligible control ID, matching/harmonization variables, approved outcomes, QC and inclusion flags | created only after Workstream 4 activation |

Each table receives a sidecar metadata file containing creation time, git commit, source paths, source modification times or hashes, row counts, variable definitions/units, and code version. Stable IDs—not row positions or newly generated random UUIDs—must link derived outputs.

### Freeze and QC gate

The following must pass before inferential execution:

1. Exactly 38 analytical sessions link uniquely to 19 participants, treatment, period, sequence, and recording start.
2. Treatment/sequence counts reproduce Phase 1: 8 placebo→THC/CBD and 11 THC/CBD→placebo; all washouts are at least seven days.
3. The conventional-arousal and >30-second audits reproduce Phase 1 counts (6,090 events; 15 above 30 seconds) before analysis-specific filtering.
4. Elapsed-time anchors reproduce recording start, scored sleep onset, and event/transition timing for sampled worked cases; no nominal dose time is stored as observed timing.
5. Any endpoint carried into Workstream 2 reproduces its original analysis population and estimate before adding elapsed-time terms.
6. At least ten randomly selected arousals plus all boundary/extreme cases are manually reconciled to source epochs, including the 15-second stage-reassignment rule.
7. Sequence coding reproduces 8 placebo→THC/CBD and 11 THC/CBD→placebo participants and is verified against chronological recording dates.
8. Model matrices are full rank; convergence, singularity, separation, residual behavior, and influential participants are recorded for every fitted model.
9. Healthy-control tables are not created or analyzed until the ethics and eligibility gate in Workstream 4 is passed.

Failure triggers a documented protocol amendment or blocks the affected workstream; it must not be silently repaired after results are viewed.

## Workstream 1 — Arousal-to-awakening timing

### Question and estimand

How much does the interval from arousal onset to the epoch boundary carrying the outcome classification vary, and does the randomized treatment contrast in awakening probability materially change with target position in the onset epoch?

The boundary interval is a property of the scoring rule, **not physiological arousal-to-wake latency**.

### Population and derivation

- Include conventional arousals eligible for the original N1, N2, N3, and REM awakening analyses.
- Derive onset position within the 30-second onset epoch on `[0,30)` and seconds to the relevant next boundary on `(0,30]`, after first verifying the export’s epoch indexing with worked cases.
- Preserve flags for events reassigned from an early wake/N1 epoch to the preceding sleep stage; report these separately because their classification path differs.
- Summarize the interval continuously and in locked thirds: 0–<10, 10–<20, and 20–<30 seconds from epoch start.

### Models

1. Descriptive distributions by assigned stage, awakening label, treatment, reassignment status, participant, and epoch-position third.
2. Main sensitivity model: original awakening GLME plus centered epoch position and `treatment × epoch_position`; retain stage, capped duration, original treatment-by-stage and stage-by-duration terms, and participant random intercept.
3. Stage-specific sensitivity models: `awakening ~ treatment * epoch_position + duration + (1|participant)`.
4. Short-boundary sensitivity: reproduce stage-specific treatment estimates among events in the final epoch third (≤10 seconds to the next boundary), with counts and coverage emphasized.

If a model is unsupported because of sparse cells/separation, report descriptives and the short-boundary sensitivity without escalating model complexity.

### Outputs and acceptance

- Event-flow table and boundary-interval histogram/ECDF.
- Treatment estimates by stage overall and in the short-boundary subset.
- Explicit wording that `next_stage` is epoch-level and no exact wake-onset latency was measured.
- Acceptance requires worked-case verification, complete participant/event attrition, and reproduction of the original model estimate before adding timing terms.

## Workstream 2 — Evolution across the night

### Question and estimand

Does the randomized treatment contrast change linearly with observed elapsed sleep time for selected event-level outcomes?

This is not a measured time-since-dose analysis.

### Timing axes

- Primary: minutes since first scored sleep onset.
- Complementary sensitivity: minutes since observed recording start/effective lights off.
- Clock time and cycle number are descriptive metadata only unless a pre-result amendment establishes that the primary axis is not computable.
- Plot model estimates at 60, 180, 300, and 420 minutes when supported by the data; do not extrapolate beyond observed overlap between treatments.

### Prespecified outcome groups and models

1. **Awakening probability (principal):** stage-specific event GLMEs `awakening ~ treatment * elapsed_time + duration + (1|participant)`. The estimand is the treatment-by-time odds-ratio interaction per hour.
2. **NREM infraslow fluctuation frequency, amplitude and HR-coupling (secondary):** fit `ISF_feature ~ treatment * elapsed_time + (1|participant)`. This tests temporal heterogeneity of the established contrast and is not an independent discovery test. If this is interaction effect is significant, then we can extract Frequency or Amplitude or HR coupling oer time and show how that changes across the night for each of the conditions. 

I don't think we can adjust the analyses for ISF dynamics and NREM sleep stability or for ISF dynamics and REM sleep transitions.

### Acceptance

- Report support and overlap of elapsed time by treatment and outcome.
- Conclusions are limited to evolution over observed sleep/recortime and cannot be translated into drug concentration decay.

## Workstream 3 — Treatment sequence/order

### Question and estimand

Are principal treatment estimates materially altered after accounting for randomized treatment sequence, and is there evidence that the treatment contrast differs between the placebo→THC/CBD and THC/CBD→placebo sequence groups?

### Models

For each principal endpoint, first reproduce the original model and then add treatment sequence as a fixed between-participant design term while retaining the original participant random intercept and endpoint-specific covariates. Include study period as a fixed nuisance term where needed to separate sequence/order from first-versus-second-visit effects; period is not a separate reviewer-requested analysis.

Whilst updating the models listed above, also add in `treatment × sequence` one endpoint at a time as the exploratory treatment-order heterogeneity test; do not fit three-way crossover interactions in 19 participants. Sequence is a between-participant randomized-order term.

### Acceptance

- Sequence must be reconstructed only from chronological observed recording dates and treatment labels and match Phase 1 counts.
- Report sequence-adjusted versus original treatment estimates and CIs side by side, model convergence, and treatment-by-sequence estimates.
- A sequence main effect is contextual; the treatment-by-sequence term directly addresses whether the treatment contrast differs by order.

## Workstream 4 — Healthy-control analyses

**Status: gated; ethics amendment/approval is under review.** Phase 1 found technical evidence for candidate controls but did not establish source-study eligibility, consent/permissions, healthy screening, acquisition/scoring comparability, restrictions, adaptation status, or the reason 22 labeled controls yield only 18 processed records.

No control model will run unless ethics approval is documented, all Phase 1 eligibility-gate domains are resolved, and a dated amendment prespecifies:

- the eligible control population and exclusions;
- matching/weighting or covariate adjustment;
- harmonized PSG/EEG variables and QC;
- missing-data handling and multiplicity; and
- the observational estimand and its limits.

If activated, controls can estimate observational insomnia-versus-control differences under untreated/comparable conditions. They cannot estimate a diagnosis-by-THC/CBD interaction or establish whether the randomized drug response is insomnia-specific. If documentation remains unavailable, the final Phase 2 output is a limitation-only handoff and unrelated work proceeds.

## Implementation order

1. Build and freeze the session and arousal tables; validate timing, treatment, period, and sequence fields.
2. Run the cross-cutting QC gate and reproduce the original endpoint estimates needed by Workstreams 1–3.
3. Implement Workstream 1 because it validates epoch-position derivation and the awakening model.
4. Implement Workstream 2 after elapsed sleep/recording time and each selected endpoint population are validated.
5. Implement Workstream 3 from the verified chronological sequence coding.
6. Keep Workstream 4 gated until ethics approval and all scientific eligibility requirements are documented; then add its prespecified amendment before inspecting control results.
7. Generate synchronized result summaries and handoffs for the four selected reviewer questions.

## Validation strategy

### Code and data validation

- Update existing code
- The human will oversee the code and validate its output

## Phase-level acceptance criteria

Phase 2 is complete only when:

1. Every reviewer question in the traceability table has a result, evidence synthesis, transparent unavailable status, or documented block.
2. Timing is correctly described as epoch boundary, observed time since sleep/recording, or nominal protocol context—never unsupported physiological wake latency or observed time since dose.
3. Sequence, and treatment-by-sequence are interpreted as distinct quantities.
4. Handoffs to Phases 3–5 and 7 identify exact files, estimates, limitations, and synchronized manuscript destinations.

## Immediate next actions

1. Create the tracked `code/group-level/` module and metadata/table-build functions following current MATLAB conventions.
2. Reconstruct `sessions.csv` and `arousal_events.csv`; validate all joins, epoch-position derivation, elapsed-time anchors, period/sequence, and Phase 1 counts.
3. Build deterministic endpoint extracts for the prespecified Workstream 2 outcomes, with stable IDs and no plotting side effects.
5. Record the Workstream 4 ethics decision; if approved, complete the eligibility gate and protocol amendment before creating a control analysis table.
