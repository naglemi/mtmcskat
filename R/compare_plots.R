compare_plots <- function(master_output, plot_out_name, scaffold_ID){

  # Why is this needed? Find out why it is not already numeric
  master_output$position <- as.numeric(as.character(master_output$position))

  master_output$Chr <- rep(scaffold_ID, nrow(master_output))


  graphics::par(mfrow=c(2,1))

  qqman::manhattan(
    master_output[!is.na(master_output$`SKAT_p-val`), ],
    chr = "Chr",
    bp = "position",
    p = "SKAT_p-val",
    xlim = c(min(stats::na.omit(master_output$position)),
             max(stats::na.omit(master_output$position))),
    ylim = c(0,
             ceiling(max(
               stats::na.omit((-log(master_output$`SKAT_p-val`,
                             base = 10)))))),
    suggestiveline = FALSE,
    genomewideline = FALSE)

  if(nrow(master_output[!is.na(master_output$`SKAT_p-val_resampled`), ]) > 0){
    qqman::manhattan(
      master_output[!is.na(master_output$`SKAT_p-val_resampled`), ],
      chr = "Chr",
      bp = "position",
      p = "SKAT_p-val_resampled",
      xlim = c(min(stats::na.omit(master_output$position)),
               max(stats::na.omit(master_output$position))),
      ylim = c(0,
               ceiling(max(
                 stats::na.omit((-log(master_output$`SKAT_p-val`,
                                      base = 10)))))),
      suggestiveline = FALSE,
      genomewideline = FALSE)
  }


  message(paste0("Writing ",
                 plot_out_name))

  p <- grDevices::recordPlot()
  g <- grid::grid.grabExpr(
    gridGraphics::grid.echo(p))

  ggplot2::ggsave(paste0(plot_out_name)
                  ,g, bg = "transparent")

}
