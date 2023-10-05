#' @title Get BED saturation
#'
#' @description
#' Given a set of BED files, anticipate the saturation of annotated features by
#' unifying files step-wise in a random order, counting the number of unique
#' features each time. Repeat the process for <iters> iterations, and use a
#' saturation curve model (drc::drm(fct=22.5L())) to extrapolate peaks step-wise
#' up to <extrap_points> samples. Produce a plot summarizing this data, and
#' store it as well as all tables and the saturation model to an RDS object, if
#' <cache_rds> is provided. If <cache_rds> is provided and that file exists,
#' read it directly (if <over_write> is FALSE).
#'
#' Returns a list containing:
#' 1) iterations: Mean and SD peaks added per sample.
#' 2) predictions: Table of modeled mean peaks added per sample up to <extrap_points>.
#' 3) fit: drc model and fit for saturation values.
#' 4) details: Table of details for saturation analysis (steps, iterations, etc.)
#' 5) plot: GGplot summarizing saturation.
#'
#' @param bed_files List of BED files to be analyzed. Assumes columns 1:3 are seqnames,start,end.
#' @param name_out Name column to be appended to each table. Defaults to "Saturation."
#' @param cache_rds RDS file to save results to and/or read cached results from. If not provided, no cache will be saved/loaded.
#' @param iters Number of randomized iterations to perform. Defaults to 10.
#' @param extrap_points Number of samples to extrapolate out in curve. Defaults to 500.
#' @param over_write Should existing <cache_rds> files be overwritten? Defaults to FALSE.
#'
#' @import dplyr
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import magrittr
#' @import data.table
#' @import tibble
#' @import drc
#' @import ggplot2
#' @import tictoc
#' @import scales
#' @import stats
#'
#' @export

get_bed_saturation  <- function(bed_files,name_out="Saturation",cache_rds,iters = 10,extrap_points=500,over_write=FALSE){
  # Aliases.
  mutate <- dplyr::mutate
  arrange<- dplyr::arrange
  select <- dplyr::select
  filter <- dplyr::filter
  rename <- dplyr::rename

  if(missing(cache_rds)) cache_rds <- ""
  read_bed <- function(bed_file_in,){
    fread(bed_file_in) %>%
      select(1:3) %>%
      rename(seqnames=1,start=2,end=3) %>%
      makeGRangesFromDataFrame(keep.extra.columns = FALSE)
  }
  get_union_cnt <- function(fls_in){
    gr_out <- read_bed(fls_in[[1]])
    count_list <- vector(length=length(fls_in),mode="integer")
    names(count_list) <- seq_along(count_list)
    for(i in 1:length(count_list)){
      gr_next<- read_bed(fls_in[[i]])
      gr_out <- suppressWarnings(GenomicRanges::reduce(c(gr_out,gr_next)))
      count_list[i] <- length(gr_out)
    }
    return(count_list)
  }

  if(file.exists(cache_rds) & !over_write){
    lst_out <- readRDS(cache_rds)
  }else{
    tic()
    tb_summ <- enframe(bed_files,name='name',value='bed') %>%
      mutate(idx = row_number(),
             set = name_out)

    tb_iters<- lapply(1:iters,function(i) {
      it_nm <- paste0("iter_",i)
      fl_set<- vectify(tb_summ,bed,idx) %>% sample
      get_union_cnt(fl_set) %>%
        enframe(name="step",value="count") %>%
        mutate(step = as.integer(step)) %>%
        mutate(iter = it_nm)
    }) %>%
      do.call(rbind,.) %>%
      group_by(step) %>%
      summarize(mean_count = mean(count),
                sd_count = sd(count),
                .groups="drop") %>%
      mutate(set = name_out)

    # Fit saturation curve, get predictions.
    fit <- drm(tb_iters$mean_count ~ tb_iters$step,fct = LL2.5())
    tb_pred <- enframe(predict(fit,data.frame(step=c(1:extrap_points))),
                       name = "step",value = "mean_count") %>%
      mutate(set = name_out)
    tb_iters<- mutate(tb_iters,frac = mean_count / max(tb_pred$mean_count))

    # Annotation table.
    tb_scales <- tibble(set = name_out,
                        detected=max(tb_iters$mean_count) %>% as.integer,
                        predicted=max(tb_pred$mean_count) %>% as.integer,
                        steps = nrow(tb_iters),
                        extrap_steps = extrap_points,
                        iters = iters,
                        ymax = 1.1 * (max(tb_pred$mean_count) + max(tb_iters$sd_count))) %>%
      mutate(frac = detected / predicted,
             label_detected = paste0("Detected: ",comma(detected)," (",percent(frac,accuracy=1),")"),
             label_predicted= paste0("Predicted: ",comma(predicted)," (",extrap_steps," samples)")) %>%
      select(set,detected,predicted,frac,steps,iters,extrap_steps,ymax,label_detected,label_predicted)

    # Plot out.
    p <- ggplot(tb_iters,aes(x=step,y=mean_count,
                             ymin=mean_count - sd_count,
                             ymax=mean_count + sd_count)) +
      scale_x_continuous(name="Samples") +
      scale_y_continuous(name="Unique peaks",
                         labels = function(x) case_when(x==0 ~ "0",
                                                        x > 1e3 ~ paste0(x/1e3,"k"),
                                                        x > 1e6 ~ paste0(x/1e6,"M"),
                                                        TRUE ~ as.character(x))) +
      geom_rect(data=tb_scales,mapping=aes(xmax=steps,ymax=ymax,xmin=1,ymin=0),
                fill=NA,color=NA,inherit.aes=FALSE) +
      geom_line(data=filter(tb_pred,step<= max(tb_iters$step)),
                mapping=aes(x=step,y=mean_count),color="red",linewidth=0.5,inherit.aes=FALSE) +
      geom_hline(data=tb_scales,mapping=aes(yintercept=detected),linewidth=0.5,color="gray30",linetype='dashed') +
      geom_hline(data=tb_scales,mapping=aes(yintercept=predicted),linewidth=0.5,color="red",linetype='dashed') +
      geom_label(data=tb_scales,mapping=aes(x=1,y=detected,label=label_detected),
                 color="gray30",hjust=0,vjust=0.9,inherit.aes=FALSE) +
      geom_label(data=tb_scales,mapping=aes(x=1,y=predicted,label=label_predicted),
                 color="red",hjust=0,vjust=0.1,inherit.aes=FALSE) +
      ggtitle(name_out,subtitle = paste0(iters," iterations")) +
      geom_errorbar(linewidth=0.5,width=0.4) +
      theme(plot.background = element_rect(fill="white",color=NA),
            panel.border = element_rect(fill=NA,color="black",linewidth=1),
            panel.background = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(linewidth=0.5,color="gray76"))

    lst_out <- list(iterations = tb_iters,
                    predictions= tb_pred,
                    fit = fit,
                    details = tb_scales,
                    plot = p)

    if(cache_rds != ""){
      saveRDS(lst_out,cache_rds)
    }
    toc()
  }
  return(lst_out)
}
