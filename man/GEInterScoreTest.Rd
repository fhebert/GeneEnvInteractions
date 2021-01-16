\name{GEInterScoreTest}
\alias{GEInterScoreTest}
\title{Vector of test statistics for gene-environment interaction effects}
\usage{
GEInterScoreTest(X, E, Y, perm.method = "parametric", N = 1000, U = NULL)
}
\arguments{
\item{X}{a n x p matrix of genotypes (coded as 0/1/2)}

\item{E}{a vector of length n containing the values of the environmental variable}

\item{Y}{a vector of length n containing the values of the binary phenotype (0/1)}

\item{perm.method}{a character giving the chosen resampling method for the phenotype: paramatric bootstrap (\code{perm.method="parametric"}, default value) or simple permutations (any other character)}

\item{N}{the number of permutations of the phenotype (default = 1000)}

\item{U}{an optional matrix of additional covariates}
}
\value{
A list with 5 elements:
\item{Z}{the vector of test statistics}
\item{Z0}{the matrix of permuted vectors of test statistics, displayed as rows}
\item{Sigma}{the correlation matrix corresponding to \code{Z}}
\item{Y0}{the matrix of permuted phenotypes, displayed as columns}
\item{S}{the design matrix for the interaction effects}
\item{W}{the design matrix for the marginal effects}
}
\description{
\code{GEInterScoreTest} is used to compute the vector of score test statistics for testing the
presence of an interaction effect between a gene and an environmental variable on a binary phenotype.
}
