\name{SampleGEInter}
\alias{SampleGEInter}
\title{Simulation of samples under a GxE interaction model}
\usage{
SampleGEInter(X, E, beta0, betaX, I.betaX, betaE, gamma, I.gamma, n0, n1,
U = NULL, betaU = NULL)
}
\arguments{
\item{X}{a N x p matrix of genotypes (coded as 0/1/2) for a given gene for a population}

\item{E}{a vector of size N containing the values of the environmental variable for the same population}

\item{beta0}{the intercept of the logistic model}

\item{betaX}{the vector of parameters for the main genetic effects}

\item{I.betaX}{the vector of coordinates corresponding to the coefficients in \code{beta}}

\item{betaE}{the parameter for the main effect corresponding to the environmental variable}

\item{gamma}{the vector of parameters for the interaction effects}

\item{I.gamma}{the vector of coordinates corresponding to the coefficients in \code{gamma}}

\item{n0}{the desired number of controls}

\item{n1}{the desired number of cases}

\item{U}{an optional N x q matrix of additional covariates}

\item{betaU}{the vector of coefficients for the covariates}
}
\value{
A list with 3 elements:
\item{SNP}{the sample matrix of genotype values}
\item{Env}{the sample vector of values for the environmental variable}
\item{Y}{the vector of phenotype values}
\item{Covariates}{the sample matrix of covariates}
}
\description{
\code{SampleGEInter} is used to generate a sample under a gene - environment interaction model. The phenotype is generated conditionally to the gene and the environmental variable using a logistic model.
}
