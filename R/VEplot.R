#' @export
#' @import tidyverse ggplot2
VEplot <- function(vaccine, infection_ind) {
  n_type = length(vaccine)
  plot_result = NULL
  for (i in 1:n_type) {
    xlabel = ifelse(i < infection_ind, "Months since vaccination", "Months since infection")
    tmp_dat = as.data.frame(vaccine[[i]])
    xmax = max(tmp_dat$time[which(tmp_dat$VE>=0)],0)
    tmp_plot = ggplot(data = tmp_dat) +
      geom_line(aes(x = time, y = VE*100), size = 0.75) +
      geom_ribbon(aes(x = time, ymin=pmax(`lower .95`*100, 0), ymax=`upper .95`*100), alpha = 0.3, show.legend = F) +
      scale_y_continuous(name = "Effectiveness (%)", breaks=seq(0, 100, 10), limits = c(0,100), expand = c(0, 0))+
      scale_x_continuous(name = xlabel, expand = c(0, 0), limits = c(0,xmax)) +
      theme_light() +
      theme(legend.key = element_rect(colour = NA, fill = NA),
            legend.background = element_blank(),
            aspect.ratio=1/1.5)
    plot_result[[i]] = tmp_plot
  }
  return(plot_result)
}
