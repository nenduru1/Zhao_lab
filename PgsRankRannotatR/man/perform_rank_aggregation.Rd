\name{perform_rank_aggregation}
\alias{perform_rank_aggregation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Performs rank aggregation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
perform_rank_aggregation(data_path, file_prefix, ranks_column)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data_path}{
%%     ~~Describe \code{data_path} here~~
}
  \item{file_prefix}{
%%     ~~Describe \code{file_prefix} here~~
}
  \item{ranks_column}{
%%     ~~Describe \code{ranks_column} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (data_path, file_prefix, ranks_column) 
{
    library(dplyr)
    library(tidyr)
    library(RobustRankAggreg)
    library(scales)
    library(readr)
    data <- read.csv(data_path)
    data_InverseRank <- data \%>\% mutate(inv_ranks = 1/{
        {
            ranks_column
        }
    })
    data_wide <- data_InverseRank \%>\% dplyr::select(ID, SNP_coord, 
        inv_ranks) \%>\% tidyr::spread(key = ID, value = inv_ranks, 
        fill = 0)
    data_wide_MRR <- data_wide \%>\% dplyr::rowwise() \%>\% dplyr::mutate(avg_inv_rank = mean(c_across(-SNP_coord)))
    data_MRR <- data_wide_MRR \%>\% dplyr::select(SNP_coord, avg_inv_rank)
    pgs_all_metadata <- read_excel("/data2/gushijimamwesigwa/projects/chat_GPT/pgs_all_metadata.xlsx", 
        sheet = "Score Development Samples")
    column_index <- which(names(pgs_all_metadata) == "Polygenic Score (PGS) ID")
    names(pgs_all_metadata)[column_index] <- "ID"
    pgs_sample_size <- pgs_all_metadata \%>\% group_by(ID) \%>\% 
        summarize(Sum_Number_of_Individuals = sum(`Number of Individuals`, 
            na.rm = TRUE))
    unique_ids <- unique(data$ID)
    sample_size_weight <- data.frame(ID = unique_ids)
    sample_size_weight <- merge(sample_size_weight, pgs_sample_size, 
        by = "ID", all.x = TRUE)
    sample_size_weight$convex_weights <- sample_size_weight$Sum_Number_of_Individuals/sum(sample_size_weight$Sum_Number_of_Individuals)
    weights <- sample_size_weight$convex_weights
    weights <- as.numeric(weights)
    data_wide <- as.data.frame(data_wide)
    data_matrix <- data_wide[, -1]
    data_matrix <- as.matrix(data_matrix)
    rownames(data_matrix) <- data_wide[, 1]
    weighted_values <- t(t(data_matrix) * weights)
    weighted_means <- rowSums(weighted_values)
    weighted_MRR <- data.frame(weighted_means)
    weighted_MRR$SNP_coord <- row.names(weighted_MRR)
    data_wide <- data \%>\% dplyr::select(ID, SNP_coord, {
        {
            ranks_column
        }
    })
    data_wide <- tidyr::spread(data_wide, key = ID, value = {
        {
            ranks_column
        }
    }, fill = NA)
    unique_sorted_groups <- sort(unique(data$ID))
    group_counts <- table(data$ID)
    sorted_counts <- unname(group_counts[match(unique_sorted_groups, 
        names(group_counts))])
    max_rank = nrow(data_wide)
    data_wide <- data_wide \%>\% mutate_all(~ifelse(is.na(.), max_rank, 
        .))
    data_wide[, -1] <- lapply(data_wide[, -1], as.numeric)
    r = as.matrix(data_wide[, -1])
    row.names(r) <- data_wide$SNP_coord
    if (!require(scales)) {
        install.packages("scales")
        library(scales)
    }
    r_norm <- as.data.frame(apply(r, 2, rescale))
    r_norm <- as.matrix(r_norm)
    aggreg_rank_RRA <- aggregateRanks(rmat = r_norm, method = "RRA", 
        exact = TRUE)
    aggreg_rank_stuart <- aggregateRanks(rmat = r_norm, method = "stuart", 
        exact = TRUE, topCutoff = sorted_counts)
    aggreg_rank_min <- aggregateRanks(rmat = r_norm, method = "min", 
        exact = TRUE, topCutoff = sorted_counts)
    aggreg_rank_geo <- aggregateRanks(rmat = r_norm, method = "geom.mean", 
        exact = TRUE, topCutoff = sorted_counts)
    aggreg_rank_mean <- aggregateRanks(rmat = r_norm, method = "mean", 
        exact = TRUE, topCutoff = sorted_counts)
    aggreg_rank_median <- aggregateRanks(rmat = r_norm, method = "median", 
        exact = TRUE, topCutoff = sorted_counts)
    rename_columns <- function(dataframe, substitution) {
        colnames(dataframe) <- sub("^Name$", "SNP_coord", colnames(dataframe))
        colnames(dataframe) <- sub("^Score$", substitution, colnames(dataframe))
        return(dataframe)
    }
    aggreg_rank_RRA <- rename_columns(aggreg_rank_RRA, "RRA_score")
    aggreg_rank_stuart <- rename_columns(aggreg_rank_stuart, 
        "Stuart_score")
    aggreg_rank_min <- rename_columns(aggreg_rank_min, "Min_score")
    aggreg_rank_geo <- rename_columns(aggreg_rank_geo, "Geo_score")
    aggreg_rank_mean <- rename_columns(aggreg_rank_mean, "Mean_score")
    aggreg_rank_median <- rename_columns(aggreg_rank_median, 
        "Med_score")
    data_wide_anno <- data \%>\% dplyr::select(ID, SNP_coord, hm_chr, 
        hm_pos, hm_rsID, annot.symbol, annot.type, ranks) \%>\% 
        tidyr::spread(key = ID, value = ranks, fill = NA)
    data_wide_anno_MRR <- left_join(data_wide_anno, data_MRR, 
        by = "SNP_coord")
    data_wide_anno_MRR <- left_join(data_wide_anno_MRR, weighted_MRR, 
        by = "SNP_coord")
    data_aggreg_rank <- left_join(data_wide_anno_MRR, aggreg_rank_RRA, 
        by = "SNP_coord") \%>\% left_join(aggreg_rank_stuart, by = "SNP_coord") \%>\% 
        left_join(aggreg_rank_min, by = "SNP_coord") \%>\% left_join(aggreg_rank_geo, 
        by = "SNP_coord") \%>\% left_join(aggreg_rank_mean, by = "SNP_coord") \%>\% 
        left_join(aggreg_rank_median, by = "SNP_coord")
    data_aggreg_rank <- data_aggreg_rank \%>\% mutate(RRA_score = rank(RRA_score), 
        Stuart_score = rank(Stuart_score), Min_score = rank(Min_score), 
        Geo_score = rank(Geo_score), Mean_score = rank(Mean_score), 
        Med_score = rank(Med_score), avg_inv_rank = rank(-avg_inv_rank), 
        weighted_means = rank(-weighted_means))
    final_data <- data_aggreg_rank \%>\% dplyr::select(hm_rsID, 
        hm_chr, hm_pos, annot.symbol, annot.type, weighted_means, 
        avg_inv_rank, RRA_score, Stuart_score, Min_score, Geo_score, 
        Mean_score, Med_score)
    file_name <- paste0(file_prefix, "_aggregated_ranks.csv")
    write.csv(final_data, file = file_name, row.names = FALSE)
    return(final_data)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
