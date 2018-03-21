\name{analysis1dat}
\alias{analysis1dat}
\docType{data}
\title{Data package for Analysis 1 type of workflow.
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
This dataset comprises the following:
1. lead_snps:  a vector of snp identifiers used in the analysis
2. populations: a vactor of population identifiers used in the anlysis
3. analysis_intervals:  GRanges objects containing the 650k interval around a lead_snp
4. results: a list of LDAanalysis1 type of objects which comprise:
  - 1000Genomes data as SnpMatrix
  - ...
5. tag_snp_intervals: a GRanges object comprising all tag-snp intervals

}
\usage{data("analysis1dat")}
\format{
  A data frame with 0 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(Analysis 1)
## maybe str(Analysis 1) ; plot(Analysis 1) ...
}
\keyword{datasets}
