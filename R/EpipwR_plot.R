
#' @title Plotting EpipwR Output
#' @description
#' Plots 95% error bars of EpipwR output.
#' Different sample sizes are displayed along the x-axis with lines
#' representing distinct correlation means. This functions should
#' only be used with EpipwR output.
#' @param df A data frame containing output from the `get_power_cont()` function.
#' @examples
#' out <- get_power_cont(100,10000,c(70,80,90,100),.05,c(.3,.4,.5))
#' EpipwR_plot(out)
#' @export
EpipwR_plot<-function(df){
  df$moe_min <- df$avg_power - 1.96*df$se_power
  df$moe_max <- df$avg_power + 1.96*df$se_power
  ggplot2::ggplot(df, aes(x = sample_size, y = avg_power, color = as.character(rho_mu))) +
    ggplot2::geom_line(linewidth = 1.25) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_errorbar(aes(ymin = moe_min, ymax = moe_max), linewidth = 1.25) +
    ggplot2::geom_hline(yintercept=0.8, linetype="dashed", color = "black") +
    ggplot2::geom_hline(yintercept=0.9, linetype="dashed", color = "black") +
    ggplot2::labs(title = "Mean power curve(s) with 95% confidence intervals",
         x = "Sample size",
         y = "Power",
         color=expression(mu[rho])) +
    # theme_bw() +
    ggplot2::ylim(0,1) +
    ggplot2::theme(text=element_text(size=18))
}

