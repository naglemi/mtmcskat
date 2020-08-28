compare_plots <- function(master_output, plot_out_name, scaffold_ID){

  # Why is this needed? Find out why it is not already numeric
  master_output$position <- as.numeric(as.character(master_output$position))



  master_output$Chr <- rep(scaffold_ID, nrow(master_output))

  par(mfrow=c(2,1))

  manhattan(master_output[!is.na(master_output$`SKAT_p-val`), ],
            chr = "Chr",
            bp = "position",
            p = "SKAT_p-val",
            xlim = c(min(master_output$position), max(master_output$position)),
            suggestiveline = FALSE,
            genomewideline = FALSE)

  manhattan(master_output[!is.na(master_output$`SKAT_p-val_resampled`), ],
            chr = "Chr",
            bp = "position",
            p = "SKAT_p-val_resampled",
            xlim = c(min(master_output$position), max(master_output$position)),
            suggestiveline = FALSE,
            genomewideline = FALSE)

  p <- recordPlot()
  g <- grid.grabExpr(grid.echo(p))
  ggsave(paste0(plot_out_name)
         ,g, bg = "transparent")
  message(paste0("Writing ",
                 plot_out_name))

}
