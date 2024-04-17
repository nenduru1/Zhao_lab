# PgsRankRannotatR

PgsRankRannotatR is an R wrapper package designed to facilitate the aggregation of polygenic scores from the PGS Catalog for any particular trait, and perform rank aggregation and annotation of effect variants.

## Installation


### Prerequisites

Before installing `PgsRankRannotatR` you must have R installed on your system (version 3.5.0 or higher is recommended, but `PgsRankRannotatR` was developed and tested using R version 4.1.1). You can download it from [CRAN](https://cran.r-project.org).

### Installing Bioconductor and Dependencies

`PgsRankRannotatR` relies on several packages from Bioconductor. To install these dependencies, follow the steps below:

1. **Install Bioconductor Manager**:
   If you do not have `BiocManager` installed, you can install it using the following command in your R console:
   ```R
   if (!requireNamespace("BiocManager", quietly = TRUE))
       install.packages("BiocManager")
   ```

2. **Set Bioconductor version and install dependencies**:
   Set the Bioconductor version to `3.14` and install the required packages:
   ```R
   BiocManager::install(version = '3.14', dependencies = TRUE)
   ```

3. **Install specific Bioconductor packages**:
   You will need to install the following packages. These commands will handle their installation:
   ```R
   BiocManager::install("annotatr")
   BiocManager::install("GenomicRanges")
   BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
   BiocManager::install("org.Hs.eg.db", version = '3.14')
   ```

### Installing PgsRankRannotatR

After installing all necessary dependencies, you can install `PgsRankRannotatR` from GitHub using `devtools`:

```R
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("savannahmwesigwa/Zhao_lab/PgsRankRannotatR")
```

### Getting Started

Once installed, you can load `PgsRankRannotatR` using:

```R
library(PgsRankRannotatR)
```

### Clone repo with test files
```R
install.packages("git2r")
library(git2r)
repo_url <- "https://github.com/savannahmwesigwa/Zhao_lab.git"
dest_path <- "./AlzheimerPGS_Tests"  # Specify your local directory here
repo <- git2r::clone(repo_url, dest_path)
setwd("./AlzheimerPGS_Tests")
```

# Functions
```R
# generate_full_dataset("trait name", "data_path", "AD")
# For example, using Alzheimer's disease
generate_full_dataset("Alzheimer", "../Alzheimer_PGS", "AD")
```
data_path is the path to the local directory where PGSs from the PGS catalog have been downloaded.
This will generate an output file with the specified prefix *_full_dataset.csv
In the above example, the file is AD_full_dataset.csv and can be used for the next step of annotating variants
```R
# process_data("path_to_full_dataset.csv", "file_prefix")
process_data("./AD_full_dataset.csv", "AD")
```
Annotates variants and ranks variants based on absolute values of effect weights.
This produces an output file *_annotated_dataset.csv
In the above example, the file is AD_annotated_dataset.csv
Adds a column named "ranks" that will be used for rank aggregation in the next step.
# Example usage of the perform_rank_aggregation function
```R
# perform_rank_aggregation("path_to_annotated_csv", "File_prefix", "column_to_rank")
perform_rank_aggregation("AD_annotated_dataset.csv," "AD," ranks)

```
This performs rank aggregation of the variants across multiple PGSs on column "ranks"
