\name{generate_full_dataset}
\alias{generate_full_dataset}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generate_full_dataset(trait_term, directory_path, output_file_prefix)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trait_term}{
%%     ~~Describe \code{trait_term} here~~
}
  \item{directory_path}{
%%     ~~Describe \code{directory_path} here~~
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
function (trait_term, directory_path, output_file_prefix) 
{
    PGS_traits <- get_traits(trait_term = trait_term, exact_term = FALSE)
    trait_vector <- PGS_traits@pgs_ids[4]
    suffix <- "_hmPOS_GRCh38.txt"
    results <- lapply(trait_vector, function(item) {
        full_path <- file.path(directory_path, paste0(item, suffix))
        read_scoring_file(full_path)
    })
    file_names <- names(results$pgs_id)
    data_list <- lapply(file_names, function(file_name) {
        cleaned_name <- basename(file_name)
        cleaned_name <- sub("_hmPOS_GRCh38.txt$", "", cleaned_name)
        data <- results$pgs_id[[file_name]]$data
        data$ID <- cleaned_name
        return(data)
    })
    full_dataset <- dplyr::bind_rows(data_list)
    output_file <- paste0(output_file_prefix, "_full_dataset.csv")
    write.csv(full_dataset, file = output_file, row.names = FALSE)
    cat("Full dataset saved as", output_file, "\n")
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
