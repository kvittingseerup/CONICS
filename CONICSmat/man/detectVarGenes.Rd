% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Tsne.R
\name{detectVarGenes}
\alias{detectVarGenes}
\title{Identify the most variable genes in a matrix}
\usage{
detectVarGenes(mat, ngenes = 500, FDR = 1e-05)
}
\arguments{
\item{ngenes}{Return the ngenes most variable genes.}

\item{expmat}{A genes X samples expression matrix of log2(CPM/10+1) scaled (single cell) RNA-seq counts.}

\item{ngenes}{FDR threshold}
}
\description{
This function identifies the most variable genes from a matrix of raw gene counts to log2(CPM/10+1).
}
\examples{
detectVarGenes(suva_exp,500)
}
\keyword{Detect}
\keyword{genes}
\keyword{varibale}
