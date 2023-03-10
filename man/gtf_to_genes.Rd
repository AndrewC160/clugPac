% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gtf_to_genes.R
\name{gtf_to_genes}
\alias{gtf_to_genes}
\title{GTF to genes}
\usage{
gtf_to_genes(
  gtf_file,
  cache_file = NULL,
  overwrite_cache = FALSE,
  as_grange = TRUE
)
}
\arguments{
\item{gtf_file}{Filename of GFF with gene annotations. If not provided, defaults to packaged "Homo_sapiens.GRCh38.104.chr.tabix.gtf.gz"}

\item{cache_file}{Filename of TSV file in which to store annotations. Defaults to NULL.}

\item{overwrite_cache}{Boolean, defaults to FALSE. Should TSV cache file be rebuilt from scratch?}

\item{as_grange}{Should data be returned as a GRange object? Defaults to TRUE.}
}
\description{
Given a GTF file, extract all genes and produce a GenomicRanges object. If a
cache file is provided, save a TSV copy in that location. If such a file
already exists and <overwrite_cache> == FALSE, load the existing TSV file.
}
