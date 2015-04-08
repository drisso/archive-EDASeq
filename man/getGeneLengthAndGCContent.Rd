\name{getGeneLengthAndGCContent}
\alias{getGeneLengthAndGCContent}

\title{Get gene length and GC-content}

\description{
Automatically retrieves gene length and GC-content information from Biomart or org.db packages.
}

\usage{
getGeneLengthAndGCContent(id, org, mode=c("biomart", "org.db"))
}

\arguments{
  \item{id}{one or more gene IDs (ensembl or entrez).}
  \item{org}{organism three letter code, e.g. 'hsa' for 'Homo sapiens'.}
  \item{mode}{which mode to use to retrieve the information.}
}

\details{
In 'biomart' mode, it supports all ENSEMBL organisms, but might be time-consuming.

In 'org.db' mode, it is based on Bioconductor's annotation packages. It is much faster,
but only for organisms with a TxDb, BSgenome and OrgDb package. It assumes that the
packages are installed.
}

\value{
A numeric matrix with two columns: gene length and GC-content.
}

\author{
Ludwig Geistlinger
}

\example{
getGeneLengthAndGCContent("ENSG00000012048", "hsa")
}
