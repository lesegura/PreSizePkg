library(tidyverse, quietly = T)

#######################################
#                                     #
#             Parameters              #
#                                     #
#######################################

# BEFORE calculating the sample size based on precision for clustered data, the user needs:
#     rsk_e = the risk in the Exposed,
#     rsk_ne = the risk in the Unexposed
#     ratio_ene = the ratio of Unexposed group size to the Exposed group size
#     c_level = the compatibility level
#     ci_width = the desired width of the compatibility interval for RD measures and ratio of ULCI / LLCI for relative measures
#     deff = the ratio of the variance in the sample with the clustered design to the variance if the sample had been done using
#             simple random sampling


#######################################
#                                     #
#          Calculating DEFF           #
#                                     #
#######################################

deff <- function(icc, n_cluster) {
  deff_table <- 1 + icc * (n_cluster - 1)
  return(deff_table)
}


##############################################
#                                            #
#   Sample size based on precision for RD    #
#                                            #
##############################################


smp_size_rd <- function(rsk_e, rsk_ne, ratio_ene, c_level, ci_width, deff) {
  require(tidyverse, quietly = T)

  z <- qnorm(1 - (1 - c_level) / 2)

  smp_size_tab <- tibble(
    n1 = 4 * z ^ 2 * deff * (ratio_ene * rsk_e * (1 - rsk_e) + rsk_ne * (1 - rsk_ne)) / (ci_width ^ 2 * ratio_ene),
    n0 = n1 * ratio_ene,
    N = n1 + n0) %>%
    round()

  message(paste(
    paste(
    paste(
    paste(
      paste(
        paste(
          paste("Sample size based on precision for a CI width of", ci_width, sep = " "),
                  "a deff of", sep = ", "),
                   deff, sep = " "),
                   "a confidence level of", sep = ", "),
                    c_level, sep = " "), "and an unexposed to exposed group size ratio of", sep = ", "),
                    ratio_ene, sep = " "))


  return(smp_size_tab)
}

smp_size_rd(.4, .3, 3, .90, 0.08, 1) ### checking with Rothman & Greenland's example


##############################################
#                                            #
#   Sample size based on precision for RR    #
#                                            #
##############################################

smp_size_rr <- function(rsk_e,rsk_ne, ratio_ene, c_level, ci_ratio, deff) {
  require(tidyverse, quietly = T)

  z <- qnorm(1 - (1 - c_level) / 2)

  smp_size_tab <- tibble(
    n1 = 4 * z ^ 2 * deff * (ratio_ene * rsk_ne * (1 - rsk_e) + rsk_e * (1 - rsk_ne)) / (ratio_ene * rsk_e * rsk_ne * log(ci_ratio) ^ 2),
    n0 = n1 * ratio_ene,
    N = n1 + n0
  ) %>%
    round()

  message(paste(
    paste(
      paste(
        paste(
          paste(
            paste(
              paste("Sample size (N) based on precision for a ratio of ULCI to LLCI of", ci_ratio, sep = " "),
              "a deff of", sep = ", "),
            deff, sep = " "),
          "a confidence level of", sep = ", "),
        c_level, sep = " "), "and an unexposed to exposed group size ratio of", sep = ", "),
    ratio_ene, sep = " "))


  return(smp_size_tab)
}

smp_size_rr(0.4, 0.3, 3, .95, 2, 1) ### checking with Rothman & Greenland's example

###########################################################
#                                                         #
#   Sample size based on precision for Rate Differences   #
#                                                         #
###########################################################

smp_size_ird <- function(inc_e, inc_ne, ratio_ene, c_level, ci_ratio, deff) {
  require(tidyverse, quietly = T)

  z <- qnorm(1 - (1 - c_level) / 2)

  smp_size_tab <- tibble(
    n1 = 4 * z ^ 2 * deff * (ratio_ene * inc_ne + inc_e) / (ratio_ene * ci_ratio ^ 2),
    n0 = n1 * ratio_ene,
    N = n1 + n0
  ) %>%
    round()

  message(paste(
    paste(
      paste(
        paste(
          paste(
            paste(
              paste("Sample size based on precision for a CI width of", f, sep = " "),
              "a deff of", sep = ", "),
            deff, sep = " "),
          "a confidence level of", sep = ", "),
        cl, sep = " "), "and an unexposed to exposed group size ratio of", sep = ", "),
    r, sep = " "))


  return(smp_size_tab)
}


#############################################
#                                           #
#   Sample size based on precision for IR   #
#                                           #
#############################################

smp_size_irr <- function(inc_e, inc_ne, ratio_ene, c_level, ci_ratio, deff) {
  require(tidyverse, quietly = T)

  z <- qnorm(1 - (1 - c_level) / 2)


  smp_size_tab <- tibble(
    n1 = 4 * z ^ 2 * deff * (ratio_ene * inc_ne + inc_e) / (ratio_ene * inc_e * inc_ne * log(ci_ratio) ^ 2),
    n0 = n1 * ratio_ene,
    N = n1 + n0
  ) %>%
    round()

  message(paste(
    paste(
      paste(
        paste(
          paste(
            paste(
              paste("sample size based on precision for ratio of ULCI to LLCI of", ci_ratio, sep = " "),
              "a deff of", sep = ", "),
            deff, sep = " "),
          "a confidence level of", sep = ", "),
        c_level, sep = " "), "and an unexposed to exposed group size ratio of", sep = ", "),
    ratio_ene, sep = " "))


  return(smp_size_tab)
}

#############################################
#                                           #
#   Sample size based on precision for OR   #
#                                           #
#############################################

### note here m1 is the size of the cases group, m0 is the size of the controls group,
### r is the ratio of controls to cases, p1 is the prevalence of exposure in cases,
### p0 is the prevalence of exposure in controls.

smp_size_or <- function(rsk_e, rsk_ne, ratio_ene, c_level, ci_ratio, deff) {
  require(tidyverse, quietly = T)

  z <- qnorm(1 - (1 - c_level) / 2)


  smp_size_tab <- tibble(
    m1 = 4 * z ^ 2 * deff * (ratio_ene * rsk_ne * (1 - rsk_ne) + rsk_e * (1 - rsk_e)) / ((log(ci_ratio) ^ 2) * (r * rsk_e * rsk_ne * (1 - rsk_e) * (1 - rsk_ne))),
    m0 = m1 * ratio_ene,
    M = m1 + m0
  ) %>%
    round()

  message(paste(
    paste(
      paste(
        paste(
          paste(
            paste(
              paste("sample size based on precision for ratio of ULCI to LLCI of", ci_ratio, sep = " "),
              "a deff of", sep = ", "),
            deff, sep = " "),
          "a confidence level of", sep = ", "),
        c_level, sep = " "), "and an unexposed to exposed group size ratio of", sep = ", "),
    ratio_ene, sep = " "))


  return(smp_size_tab)
}



################################################################
#                                                              #
#   Precision given a fixed sample size for Risk Differences   #
#                                                              #
################################################################


precision_rd <- function(rsk_e, rsk_ne, ratio_ene, c_level, exp_size, deff) {
  require(tidyverse, quietly = T)

  z <- qnorm(1 - (1 - c_level) / 2)

  precision_tab <- tibble(
    f = 2 * z * sqrt(deff) * sqrt((ratio_ene * rsk_e * (1 - rsk_e) + rsk_ne * (1 - rsk_ne)) / (ratio_ene * exp_size)),
    rd = rsk_e - rsk_ne,
    lci = rd - z * sqrt(deff) * sqrt((ratio_ene * rsk_e * (1 - rsk_e) + rsk_ne * (1 - rsk_ne)) / (ratio_ene * exp_size)),
    uci = rd + z * sqrt(deff) * sqrt((ratio_ene * rsk_e * (1 - rsk_e) + rsk_ne * (1 - rsk_ne)) / (ratio_ene * exp_size))
    )

  message(paste(
    paste(
      paste(
        paste(
          paste(
            paste(
              paste("The estimated compatibility interval's precision is of", round(precision_tab$f, 5), sep = " "),
              "for a fixed sample size of", sep = " "),
            exp_size + exp_size*ratio_ene, sep = " "),
          "a design effect of", sep = ", "),
        deff, sep = " "), "and an unexposed to exposed group size ratio of", sep = ", "),
    ratio_ene, sep = " "))

  names(precision_tab) <- c("CI_Width", "Risk_Diff", "Lower_CI", "Upper_CI")

  return(precision_tab)
}

precision_rd(0.026, 0.006, 1, .95, 17478, 1.57) ### checking with NHIS hopelessness
