### Converter function for pwrEWAS

#' @title Using pwrEWAS effect-size distributions in EpipwR
#' @description
#' A function that determines the EpipwR inputs (`delta_mu`, `delta_sd`,
#' `det_limit`) for effect-size distribution using pwrEWAS methodology
#' (Centered at 0, truncated around 0 +/- `det_limit` and the 99.99th
#' percentile of the distribution as `maximal_delta`). The output can then
#' be used as inputs to `get_power_cc`.
#'
#' @param maximal_delta The desired 99.99th percentile (unless quantile is
#' specified) of the effect size
#' distribution. Must be between 0 and 1.
#' @param det_limit The minimum effect size to be used in a power calculation.
#' @param quantile The quantile that maximal_delta represents. Note that
#' pwrEWAS uses `quantile=0.9999`.
#' @details
#' As described in Graw, et al. (2019), users specify the "maximal delta" and
#' the detection limit to set the effect size distribution. The purpose of this
#' function is to provide a map between the pwrEWAS inputs and the EpipwR inputs
#' to generate the same effect size distribution.
#'
#' The main purpose of the function is to calculate the standard deviation of
#' the effect size distribution, since the mean is always 0 and the detection
#' limit is user-sepcified. Specifically, the standard deviation is calculated
#' such that the 99.99th percentile of a normal distribution (truncating out
#' 0 +/- `det_limit`) is `maximal_delta`. For further justification, see Graw,
#' et al. (2019). Note that this function does allow users to specify quantiles
#' other than 0.9999. If users wish to use this function for continuous EpipwR,
#' we recommend specifying a smaller `quantile` (i.e., 0.95, 0.99 etc.).
#'
#' @references Graw, S., Henn, R., Thompson, J. A., and Koestler, D. C. (2019).
#' pwrEWAS: A user-friendly tool for comprehensive power estimation for
#' epigenome wide association studies (EWAS). \emph{BMC Bioinformatics},
#' 20(1):218.
#'
#' @return A list with the calculated input values for the `get_power_cc`
#' function.
#'
#'
#' @export


pwrE_to_EpipwR <- function(maximal_delta, det_limit, Quantile=.9999){
  if(maximal_delta > 1 | maximal_delta < 0){stop("maximal_delta must be
                                                 between 0 and 1")}
  if(det_limit < 0 | det_limit > 0.1){stop("det_limit must be
                                                 between 0 and 0.1")}
  if(maximal_delta<=det_limit){stop("maximal_delta must be larger than
                                    the detection limit")}
  sigma <- uniroot(uni_func,c(det_limit + .0001, 1), dmax=maximal_delta,
                   lb=det_limit, q=Quantile,extendInt = "downX")$root
  if(sigma <= 0){stop("This setting requires a standard deviation estimate
                      that is too small for R to recognize. Choose a larger
                      maximal_delta")}
  list(delta_mu=0, delta_sd=sigma, det_limit=det_limit)

}
