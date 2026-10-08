# DENV-GenBank South America: Genomic Curation Pipeline

An automated R pipeline for the retrieval, multi-tier quality control curation, and advanced analytical visualization of Dengue virus (DENV1–4) whole-genome sequences from NCBI GenBank, specifically tailored for South American epidemiological surveillance.

## 📌 Project Overview

This repository contains a reproducible, end-to-end workflow designed to:

1. **Download** Dengue virus sequences and metadata directly from NCBI GenBank up to 2026.
2. **Filter & curate** genomes through a rigorous 4-tier Quality Control (QC) system, removing incomplete sequences, handling ambiguities/Ns, and filtering by length thresholds.
3. **Consolidate** processed metadata and high-quality `.fasta` sequences into a structured `results/` directory for downstream phylogenetic analysis.
4. **Visualize** the curated data through publication-ready, high-resolution TIFF plots, including spatiotemporal heatmaps, sequence length distributions, and QC alluvial flows.

## 📂 Repository Structure

```text
DENV_GENBANK/
├── R/                      # Custom functions and utilities
│   ├── utils.R             # Environment setup and directory management
│   └── qc_utils.R          # Quality control logic and sequence filtering
├── scripts/                # Main pipeline scripts
│   ├── 000_run.R           # Master execution script
│   ├── 01_download.R       # Data retrieval from NCBI
│   ├── 02_qc_tiers.R       # 4-tier sequence curation
│   ├── 03_report.R         # Metadata summary generation
│   ├── 04_plot_results.R   # Basic statistics and geospatial maps
│   └── 05_advanced_plots.R # Advanced analytical visualizations
├── results/                # Automated output directory (git-ignored for heavy files)
│   ├── by_country/         # Curated sequences organized by country
│   ├── plots/              # High-resolution TIFF outputs
│   │   └── advanced/       # Heatmaps, ridgeline, alluvial, and polar plots
│   ├── sequence_curation/  # Tier-by-tier QC fasta results
│   └── summary/            # Processed metadata CSVs
├── .gitignore              # Ignores heavy .fasta/.gb files and R history
└── README.md               # Project documentation
```

## 🚀 How to Use

The entire pipeline is designed to be executed sequentially via the master script.

1. Clone the repository to your local machine:

```bash
git clone https://github.com/your-username/dengue-suramerica-pipeline.git
cd dengue-suramerica-pipeline
```

Replace `your-username` with your GitHub username if needed.

2. Open the R project or set your working directory to the repository root.

3. Install the required dependencies. For example:

```r
install.packages(c(
  "ggplot2", "dplyr", "tidyr", "patchwork",
  "sf", "rnaturalearth", "ggridges", "ggalluvial"
))
```

4. Run the master script to execute the full pipeline:

```r
source("scripts/000_run.R")
```

## 📊 Outputs

All outputs are centralized in the `results/` folder. The visualizations include:

- **Spatiotemporal Heatmaps:** Historical density of DENV genomic surveillance by country.
- **Geospatial Maps:** Sequence distribution across South America using `sf` and `rnaturalearth`.
- **Ridgeline Plots:** Distribution of extracted CDS lengths validating the 9000 nt threshold.
- **Alluvial Diagrams:** Sankey flows demonstrating the sequence attrition through the 4-tier QC methodology.
- **Polar Coordinates:** Seasonal distribution of sample collections.
- **Stacked Waffle Charts:** Proportion of sequences by host origin (Clinical vs. Vector).

## 📚 References

1. R Core Team. *R: A language and environment for statistical computing*. Vienna, Austria: R Foundation for Statistical Computing; 2023. Available from: https://www.R-project.org/.

2. Wickham H. *ggplot2: Elegant Graphics for Data Analysis*. New York: Springer-Verlag; 2016.

3. Wickham H, François R, Henry L, Müller K, Vaughan D. *dplyr: A Grammar of Data Manipulation*. R package version 1.1.4. 2023. Available from: https://CRAN.R-project.org/package=dplyr.

4. Wickham H, Vaughan D, Girlich M. *tidyr: Tidy Messy Data*. R package version 1.3.1. 2024. Available from: https://CRAN.R-project.org/package=tidyr.

5. Pedersen TL. *patchwork: The Composer of Plots*. R package version 1.2.0. 2024. Available from: https://CRAN.R-project.org/package=patchwork.

6. Pebesma E. Simple Features for R: Standardized Support for Spatial Vector Data. *The R Journal*. 2018;10(1):439–446.

7. Massicotte P, South A. *rnaturalearth: World Map Data from Natural Earth*. R package version 1.0.1. 2023. Available from: https://CRAN.R-project.org/package=rnaturalearth.

8. Wilke CO. *ggridges: Ridgeline Plots in 'ggplot2'*. R package version 0.5.6. 2024. Available from: https://CRAN.R-project.org/package=ggridges.

9. Brunson JC. *ggalluvial: Alluvial Plots in 'ggplot2'*. R package version 0.12.5. 2023. Available from: https://CRAN.R-project.org/package=ggalluvial.