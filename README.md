# PgsRankRannotatR

PgsRankRannotatR is an R wrapper package designed to facilitate the aggregation of polygenic scores from the PGS Catalog for any particular trait, and perform rank aggregation and annotation of effect variants.

## Installation

You can install the development version of PgsRankRannotatR from GitHub using the `devtools` package:

```R
devtools::install_github("savannahmwesigwa/PgsRankRannotatR")
```

## Usage
# Load the package
```R
library(PgsRankRannotatR)
```


# Functions
```R
# generate_full_dataset("trait name", "data_path", "AD")
For example, using Alzheimer's disease
generate_full_dataset("Alzheimer", "../Alzheimer_PGS", "AD")
```
data_path is the path to the local directory where PGSs from the PGS catalog have been downloaded.
This will generate an output file with the specified prefix *_full_dataset.csv
In the above example, the file is AD_full_dataset.csv and can be used for the next step of annotating variants
```R
process_data("path_to_full_dataset.csv", "file_prefix")
process_data("./AD_full_dataset.csv", "AD")
```
Annotates variants and ranks variants based on absolute values of effect weights.
This produces an output file *_annotated_dataset.csv
In the above example, the file is AD_annotated_dataset.csv
Adds a column named "ranks" that will be used for rank aggregation in the next step.
# Example usage of the perform_rank_aggregation function
```R
perform_rank_aggregation("path_to_annotated_csv", "File_prefix", "column_to_rank")
perform_rank_aggregation("AD_annotated_dataset.csv," "AD," ranks)

```
This performs rank aggregation of the variants across multiple PGSs on column "ranks"
