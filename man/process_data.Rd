\name{process_data}
\alias{process_data}

\title{
Uses the annotatr package to annotate the variants with nearest gene names and ranks variants based on the absolute value of effect weight.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
process_data(full_dataset_path, output_file_prefix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{full_dataset_path}{
%%     ~~Describe \code{full_dataset_path} here~~
}
  \item{output_file_prefix}{
%%     ~~Describe \code{output_file_prefix} here~~
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
function (full_dataset_path, output_file_prefix) 
{
    options(max.print = 100)
    full_dataset <- read.csv(full_dataset_path)
    full_dataset_wide <- tidyr::spread(full_dataset, key = ID, 
        value = effect_weight, fill = NA)
    chrom_coord <- full_dataset_wide \%>\% dplyr::select(hm_chr, 
        hm_pos, hm_rsID)
    chrom_coord <- chrom_coord \%>\% dplyr::rename(chrom = hm_chr, 
        chromStart = hm_pos) \%>\% dplyr::mutate(chromEnd = chromStart + 
        1) \%>\% dplyr::select(chrom, chromStart, chromEnd, hm_rsID)
    chrom_coord$chrom <- paste("chr", chrom_coord$chrom, sep = "")
    chrom_coord <- subset(chrom_coord, chromEnd != "NA")
    require(dplyr)
    chrom_coord <- mutate_if(chrom_coord, is.numeric, as.integer)
    write.table(chrom_coord, file = "./chrom_coord.bed", sep = "\t", 
        quote = FALSE, row.names = FALSE, col.names = FALSE, 
        na = "")
    dm_regions <- read_regions(con = "./chrom_coord.bed", genome = "hg38", 
        format = "bed", extraCols = c(rsID = "character"))
    annots <- c("hg38_genes_1to5kb", "hg38_genes_promoters", 
        "hg38_genes_cds", "hg38_genes_5UTRs", "hg38_genes_exons", 
        "hg38_genes_firstexons", "hg38_genes_introns", "hg38_genes_intronexonboundaries", 
        "hg38_genes_exonintronboundaries", "hg38_genes_3UTRs", 
        "hg38_genes_intergenic", "hg38_enhancers_fantom", "hg38_basicgenes", 
        "hg38_cpgs")
    annotations <- build_annotations(genome = "hg38", annotations = annots)
    dm_annotated <- annotate_regions(regions = dm_regions, annotations = annotations, 
        minoverlap = 1L, ignore.strand = TRUE, quiet = TRUE)
    df_dm_annotated <- data.frame(dm_annotated)
    df_dm_annotated_select <- df_dm_annotated \%>\% dplyr::select(seqnames, 
        start, rsID, annot.symbol, annot.type, annot.width)
    full_dataset <- full_dataset \%>\% dplyr::mutate(SNP_coord = paste(hm_chr, 
        hm_pos, sep = "_"))
    df_dm_annotated_select <- df_dm_annotated_select \%>\% mutate(SNP_coord = paste(sub("^chr", 
        "", seqnames), (start - 1), sep = "_"))
    df_dm_annotated_filter <- df_dm_annotated_select \%>\% dplyr::group_by(SNP_coord) \%>\% 
        dplyr::filter(!is.na(annot.symbol)) \%>\% dplyr::arrange(annot.width, 
        .by_group = TRUE) \%>\% dplyr::slice(1)
    print(head(df_dm_annotated_select))
    annotated_df <- dplyr::left_join(full_dataset, df_dm_annotated_filter, 
        by = "SNP_coord")
    annotated_df <- annotated_df \%>\% dplyr::select("SNP_coord", 
        "hm_rsID", "hm_chr", "hm_pos", "effect_allele", "effect_weight", 
        "ID", "annot.symbol", "annot.type")
    annotated_df <- annotated_df \%>\% filter(SNP_coord != "NA_NA")
    annotated_df <- annotated_df \%>\% dplyr::group_by(ID) \%>\% 
        dplyr::arrange(desc(effect_weight)) \%>\% dplyr::distinct(SNP_coord, 
        .keep_all = TRUE) \%>\% ungroup()
    annotated_df <- annotated_df \%>\% dplyr::arrange(ID, effect_weight) \%>\% 
        dplyr::group_by(ID) \%>\% dplyr::mutate(ranks = rank(-abs(as.double(effect_weight)), 
        ties.method = "average")) \%>\% ungroup()
    output_file <- paste0(output_file_prefix, "_annotated_dataset.csv")
    write.csv(annotated_df, file = output_file, row.names = FALSE)
    return(annotated_df)
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
