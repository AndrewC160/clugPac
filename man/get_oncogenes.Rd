% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_oncogenes.R
\name{get_oncogenes}
\alias{get_oncogenes}
\title{Get oncogenes}
\usage{
get_oncogenes(gene_source = "asclab")
}
\arguments{
\item{gene_source}{String of oncogene source table name. Defaults to "asclab".}
}
\description{
Retrieve table of oncogenes and their IDs. If "source" is set to "bushman",
this table is retrieved from the Bushman Lab group's compilation of cancer
genes, downloaded from
http://www.bushmanlab.org/assets/doc/allOnco_June2021.tsv (Website:
http://www.bushmanlab.org/links/genelists). If set to
"asclab" (default), this table is from the file
"CancerGeneListSPCG_2023_formatted.tab" provided by Henry Martell. If set to
"oncokb"," this table is from the "OncoKB" database
(https://www.oncokb.org/) from the "oncokb_cancerGeneList.tsv" file.
}
