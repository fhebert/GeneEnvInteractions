\name{GeneEnvInteractions-package}
\alias{GeneEnvInteractions-package}
\alias{GeneEnvInteractions}
\docType{package}
\encoding{UTF-8}
\title{
  Global tests for the detection of gene-environment interaction effects on binary phenotypes
}
\description{
  The GeneEnvInteractions package provides functions to compute global p-values for the detection of
  gene-environment interaction effects on a binary phenotype. It also provides a function to generate
  samples of genotype and phenotype values under a gene-environment interaction logistic model.
}
\author{
Florian Hebert, Mathieu Emily, David Causeur

Maintainer: Florian Hebert <florian.hebert@agrocampus-ouest.fr>
}
\keyword{Gene-environment interactions}
\seealso{
  \code{\link{GEInterScoreTest}}
  \code{\link{SampleGEInter}}
}
\examples{
#Matrix of genotypes for one gene, for a population of 100,000 individuals
X = matrix(rbinom(10*100000,2,0.4),ncol=10)
#Environmental factor for the same 100,000 individuals, generated from a normal distribution
E = rnorm(100000)

#Sample of genotype, phenotype and environmental variable values. No main effects, #interaction effect between SNP 5 and the environmental variable:
tmp = SampleGEInter(X,E,-3,0,1,0,0.2,5,1000,1000)

#Score test statistics
res = GEInterScoreTest(tmp$SNP,tmp$Env,tmp$Phenotype)
}
