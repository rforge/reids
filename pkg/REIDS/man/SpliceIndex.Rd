% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions.R
\name{SpliceIndex}
\alias{SpliceIndex}
\title{"SpliceIndex"}
\usage{
SpliceIndex(GeneData, ExonData, InformativeExons = NULL,
  groups = list(group1 = NULL, group2 = NULL), paired = FALSE,
  significancelevel = NULL)
}
\arguments{
\item{GeneData}{The microarray data summarized at gene level.}

\item{ExonData}{The microarray data summairzed at exon level.}

\item{InformativeExons}{A character vector of exon IDs. As for the REIDS model probesets are filtered out by I/NI calls model and later on exon score, the remaining exons can be specified here. Only these shall be considered in the FIRMA analysis to make the results between REIDS and FIRMA more comparable}

\item{groups}{The groups of interest in the data. Default two groups are specificied byut more can added as group3, group4,...}

\item{paired}{Logical. Are the groups paired? only used if two groups are present.}

\item{significancelevel}{If specified, filtering is conducted on the p-values.}
}
\value{
A data frame wiith one line per exon. The columns conatin the gene ID, the exon ID, the ratio of the splice indices if two groups are present, a t- or F-statitic, a p-value and an adjusted p-value.
}
\description{
The SpliceIndex function computes the ratio of the splice indices of defined groups. Further, it performs the SI algorithm if the length of the groups is two and the MiDAS algorithm is more groups are specified. Both algorithms are implemented as defined by Affymetrix.
}
\details{
Given the gene level and exon level summarized data, the splice index method for the detection of alternative splicing is performed. The first step is to normalize the exon
data by taking the ratio with the gene level data. These values are referred to as the splice indices. If only two groups are specified, the ratio of their splice indices is taken as a measure for alternative splicing. The more the ratio deviates from zero, the more there is an indication of alternative splicing.
A t-test is conducted on the splice indices of the two groups to test their difference. If more than two groups are specified, an ANOVA model is fitted on the splice indices to discover with an F-test whether there is a difference between the groups somewhere. If a vector of informative exons is given 
to the function, only these are considered for the analysis. Finally, the p-values are adjusted for multiplicity and if a significance level is specified only the significant p-valuesare kept in the data frame.
}

