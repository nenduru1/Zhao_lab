# PGS Aggregate Annotate and Rank Variants
PgsRankRannotatR is an R wrapper package designed to facilitate rank aggregation and annotation of variants, particularly in the context of polygenic scores from the PGS Catalog.
# PgsRankRannotatR

PgsRankRannotatR is a package designed to perform rank aggregation and annotation tasks in R.

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

# Example usage of the perform_rank_aggregation function
```R
perform_rank_aggregation(data_path, file_prefix, ranks_column)
```
# Functions
## generate_full_dataset
Downloads and aggregates PGS from the PGS catalog.
## process_data
Annotates variants and ranks variants based on absolute values of effect weights.

## perform_rank_aggregation
Performs various algorithms of rank aggregation.
