#' @title Read CMAP file
#'
#' @description
#' Read a BioNano-formatted CMAP file containing OGM tag data (tag positions are
#' relative to contig position, not genome).
#'
#' @param cmap_fl CMAP filename.
#'
#' @import tidyr
#' @import tibble
#' @import stringr
#' @import data.table
#'
#' @export

read_cmap_file <- function(cmap_fl){
  rename <- dplyr::rename
  mutate <- dplyr::mutate

  c_nms <- paste("grep '#h'", cmap_fl) %>%
    system(intern = TRUE) %>%
    str_split("\t",simplify = TRUE) %>%
    gsub("^#. ","",.)

  paste0("grep -v '#' ",cmap_fl) %>%
    fread(cmd = .,col.names=c_nms) %>%
    as_tibble %>%
    rename(contig_id=CMapId,
           site_id = SiteID,
           coverage=Coverage,
           position = Position) %>%
    select(contig_id,site_id,position,coverage) %>%
    return()
}
