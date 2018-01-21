\name{bml}
\alias{bml}
\title{
Creates a BML object
}
\description{
Calculates a BML object
}
\usage{
bml(dataset, ntree, threshold, rep = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataset}{Data frame or matrix of 0's and 1's with genes as columns and observations as rows. Genes has to be named}
  \item{ntree}{Num random restarts, the number of random restarts for searching the tree space}
  \item{threshold}{Threshold for inferring paths, the level of detail in output path file}
  \item{rep}{Number of boostrap he number of bootstrap replicates you want to include, if nrep = 0 (default) no bootstrap will be performed}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Returns an object of class BML which is a list witht following items
%%  If it is a LIST, use
\item{num_edges }{Number of egdes}
\item{num_genes }{Number of genes}
\item{num_unpruned_edges }{Number of unpruned edges}
\item{num_noparent_after_global_pruning }{Number of nodes with no parent after global pruning}
\item{DAG}{Directed aciclic graph}
\item{Bootstrap}{Bootrstrap results}
%% ...
}
\references{
http://bml.molgen.mpg.de/
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
# Normal execution
bml(BML::col_resic, 10, 0.2)

# Normal excuction with bootstrap
bml(BML::col_resic, 10, 0.5, 20)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
