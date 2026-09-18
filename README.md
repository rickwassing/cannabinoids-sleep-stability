# Sleep Stability Analysis

## Overview

This repository contains the analysis code for a research project evaluating sleep stability in insomnia. The study investigates how sleep stability, or rather its instability, can be explained by infraslow fluctuations in sigma power, a proxy measure of LC activity during NREM sleep.

This repository is organized as a [Brain Imaging Data Structure (BIDS)](https://bids.neuroimaging.io/) dataset. Because this is a secondary analysis of an existing dataset, the raw/source data and derivatives are stored separately and are not part of this repository (see `.gitignore`); only the analysis code and lightweight BIDS metadata are tracked here.

## Project Structure

```
├── code/                     # All analysis code (BIDS convention)
│   ├── main.m                # Top-level orchestration script (run this to reproduce all results)
│   ├── css_init.m            # Environment/path initialization
│   ├── pipeline/             # Generic subject-level job runner (Proc/Files/cfg dispatch)
│   ├── subject-level/        # Per-subject EEG processing (preprocessing, segmentation, spectral features)
│   ├── firstlevel-outcomes/  # Functions that compute and save first-level (per-subject) outcome measures
│   ├── group-level/          # Group-level statistical analyses, one subfolder per manuscript aim/figure
│   ├── qc/                   # Quality control / visual inspection helpers
│   ├── utils/                # Re-usable support functions, categorised by purpose
│   ├── toolboxes/            # Third-party dependencies (untouched)
│   └── archive/              # Superseded/scratch code, kept for provenance (pipeline/analysis/scratch)
├── rawdata/                  # BIDS raw data metadata (data itself stored separately)
└── sourcedata/               # BIDS source data placeholder (data itself stored separately)
```

See [`code/README.md`](code/README.md) for a detailed traceability mapping between first-level
outcomes, group-level analyses, and the manuscript's aims/figures/tables.

## License

This project is licensed under a **Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License**.

**You are free to:**

- **Share** — copy and redistribute the material in any medium or format
- **Adapt** — remix, transform, and build upon the material

**Under the following terms:**

- **Attribution** — You must give appropriate credit to the author, provide a link to the license, and indicate if changes were made
- **NonCommercial** — You may not use the material for commercial purposes
- **ShareAlike** — If you remix, transform, or build upon the material, you must distribute your contributions under the same license

For commercial use or licensing inquiries, please contact the author directly.

See the [LICENSE](LICENSE) file for full license details or visit: https://creativecommons.org/licenses/by-nc-sa/4.0/

## Contact

**Rick Wassing**

- GitHub: [@rickwassing](https://github.com/rickwassing)
- Email: rick.wassing@mq.edu.au

## Data Privacy and Ethics

This research follows ethical guidelines for human subjects research. HREC Approval was granted by Bellberry Ltd (2018-04-284). This repository only contains analysis code. All data is stored separately and can be requested.
