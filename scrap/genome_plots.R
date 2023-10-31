library(ggplot2)
library(tidyverse)
library(clugPac)
library(hicPac)
library(GenomicRanges)
library(scales)
library(data.table)
library(magrittr)

tb_clusts <- fread("/N/p/asclab/ASC-cutNrun/data/atac_clusters.tsv") %>%
  dplyr::rename(
    cell_line = ATAC.id,
    atac_cluster = Atac.cluster) %>%
  as_tibble %>%
  mutate(atac_cluster = factor(atac_cluster))

tb_bdgs <- tibble(bdg=Sys.glob("/N/p/asclab/ASC-cutNrun/data/pools/*/*treat_pileup.bdg.gz")) %>%
  mutate(sample_name = basename(dirname(bdg)))
tb_npks <- tibble(nPks=Sys.glob("/N/p/asclab/ASC-cutNrun/data/pools/*/narrowPeaks/*_flat.narrowPeak")) %>%
  mutate(sample_name = basename(dirname(dirname(nPks))))
tb_ses  <- tibble(ses=Sys.glob("/N/p/asclab/ASC-cutNrun/data/pools/*/superenhancers/*SuperStitched_REGION_TO_GENE.txt")) %>%
  mutate(sample_name = basename(dirname(dirname(ses))))
tb_cnr  <- tb_bdgs %>%
  full_join(tb_npks,by="sample_name") %>%
  full_join(tb_ses,by="sample_name") %>%
  separate(sample_name,into=c("cell_line","epitope","condition"),remove = FALSE) %>%
  mutate(cell_line = factor(cell_line),
         epitope = factor(epitope,levels=c("K27Ac","K4Me","K4Me3"))) %>%
  left_join(tb_clusts,by="cell_line") %>%
  filter(grepl("^OS",cell_line) & condition == "NONE" & !is.na(epitope)) %>%
  select(sample_name,cell_line,atac_cluster,epitope,bdg,nPks,ses)

tb_atac   <- fread("/N/p/asclab/os_atacseq/tables/os_atacseq_files.tsv") %>%
  as_tibble %>%
  select(cell_line,bdg_file,nPks) %>%
  group_by(cell_line) %>%
  summarize_all(list) %>%
  ungroup %>%
  full_join(tb_clusts,by="cell_line") %>%
  unnest(c(bdg_file,nPks)) %>%
  group_by(cell_line,atac_cluster) %>%
  mutate(rep = row_number()) %>%
  mutate(atac_name = paste0(cell_line,"_",rep)) %>%
  select(atac_name,cell_line,rep,atac_cluster,bdg_file,nPks) %>%
  mutate(cell_line = factor(cell_line))

gr_genes <- gtf_to_genes(cache_file = "/N/p/asclab/ASC-cutNrun/data/genes.tsv",as_grange = TRUE)

gr_window <- gr_genes[gr_genes$gene_name == "MYC"] %>% resize(width=1E6,fix="center")

ses_fls <- tb_cnr %>%
  filter(epitope == "K27Ac") %>%
  vectify(ses,cell_line)


gr_ses <- lapply(names(ses_fls), function(nm) {
  fl_nm <- ses_fls[[nm]]
  fread(fl_nm) %>%
    rename_all(tolower) %>%
    dplyr::rename(seqnames=chrom,
           end = stop,
           se_rank = stitchedpeakrank,
           super=issuper) %>%
    select(seqnames,start,end,region_id,super,se_rank) %>%
    mutate(cell_line = nm)
}) %>% do.call(rbind,.) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE)


plot_base_in <- genome_plot(gr_window = gr_window) %>%
  genome_add_bed_tiles(bed_files =
    filter(tb_cnr,epitope == "K27Ac") %>% vectify(nPks,cell_line),
    bed_fills = "darkgreen",bed_alphas = 0.3,bed_colors = "blue",bed_linewidths=sample(c(0.1,0.25),replace=TRUE,size=11),bed_linetypes='dotted'
  ) %>%
  genome_add_bedgraph(bdg_files =
    filter(tb_cnr,epitope == "K27Ac") %>% vectify(bdg,cell_line),downsample_bins = 2000,
    bdg_fills = "darkgreen",bdg_alphas=1
  ) %>%
  genome_add_gr_tiles(gr_input=gr_ses,gr_name_column = "cell_line",
    gr_fills=NA,gr_colors="red",gr_linewidths=0.5,gr_linetypes="dotted")

plot_base_in$plot
