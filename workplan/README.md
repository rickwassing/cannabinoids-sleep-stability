# Peer-review revision workplan

## Purpose

This directory organizes the work required to address the peer-review comments on the manuscript about THC/CBD, infraslow EEG sigma-power fluctuations (ISFσ), sleep stability, and REM transitions.

The work is divided into seven phases so that each phase can be developed and completed in a separate chat. The files currently provide only a high-level scaffold. Before a phase is implemented, its dedicated chat should expand the relevant file by defining the issues, considering possible solutions, making methodological and editorial decisions, documenting dependencies and risks, and establishing acceptance criteria.

## Manuscript variants

The project contains two versions of the same manuscript:

- `documents/manuscript_imrd.md`: Introduction–Methods–Results–Discussion
- `documents/manuscript_irdm.md`: Introduction–Results–Discussion–Methods

Both versions must ultimately contain the same scientific content, results, qualifications, and conclusions. Their only substantive structural difference should be section order. The Results-before-Methods version may need brief explanatory context so that its Results remain understandable, but it should not duplicate unnecessary technical detail.

## Guiding principles

1. Treat every reviewer comment as a set of distinct concerns or questions rather than as a single item.
2. Address each concern with the appropriate combination of response, clarification, analysis, literature, restructuring, figure revision, or claim moderation.
3. Distinguish observed findings from mechanistic interpretations.
4. Avoid deciding detailed methods or wording before the relevant phase has been fully specified.
5. Keep both manuscript variants synchronized throughout the revision.
6. Maintain traceability from each reviewer concern to the work performed and the final response.
7. Handle the healthy-control work as a related stream that is already underway, while integrating its eventual findings and limitations into the relevant phases.

## Phase sequence

| Phase | File | Purpose |
|---|---|---|
| 1 | [Definitions and data inventory](01-definitions-and-data-inventory.md) | Establish the definitions, available variables, source information, and unresolved factual questions required by later phases. |
| 2 | [Reviewer-requested analyses](02-reviewer-requested-analyses.md) | Plan and conduct the additional and sensitivity analyses needed to answer methodological and empirical comments. |
| 3 | [Interpretation revision](03-interpretation-revision.md) | Reassess the biological, pharmacological, clinical, and therapeutic interpretation of the findings. |
| 4 | [Manuscript restructuring](04-manuscript-restructuring.md) | Improve clarity and move information to the most appropriate sections in both manuscript formats. |
| 5 | [Figure revision](05-figure-revision.md) | Improve interpretability, accessibility, and statistical annotation of the figures. |
| 6 | [Manuscript synchronization](06-manuscript-synchronization.md) | Ensure that the two manuscript-order variants remain scientifically identical. |
| 7 | [Response letter](07-response-letter.md) | Prepare a complete, point-by-point response that accounts for every component of every reviewer comment. |

## Dependencies

The phases are ordered but may overlap:

- Phase 1 supplies definitions and data availability for all later phases.
- Phase 2 supplies new results needed by Phases 3–5.
- Phase 3 determines the claims and qualifications used in Phase 4.
- Phases 4 and 5 should be substantially complete before synchronization in Phase 6.
- Phase 7 should be drafted alongside the work but finalized after Phases 2–6.
- Healthy-control findings should be integrated when available without delaying unrelated reviewer responses.

## How to use a phase file in a separate chat

At the beginning of a phase-specific chat:

1. Read this README, the selected phase file, the peer-review comments, and the relevant manuscript sections.
2. Decompose the phase into individual issues and reviewer sub-questions.
3. Define the scope and boundaries of the phase.
4. Brainstorm and compare possible solutions.
5. Select an approach and document the rationale.
6. Specify inputs, tasks, outputs, dependencies, risks, and acceptance criteria.
7. Only then implement the agreed work.
8. Update the phase file with decisions, progress, evidence, and completion status.

## Current status

These files are initial high-level plans only. Detailed analytical specifications, manuscript wording, implementation decisions, and acceptance criteria are intentionally deferred to the separate phase-specific chats.
