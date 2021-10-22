#' Title
#'
#' @return
#' @export
#'
#' importFrom magrittr "%>%"
#'
#' @examples
prep_cps03 <- function() {
  
  library(magrittr)
  
  ddi <- ipumsr::read_ipums_ddi(here::here("data", "cps_00009.xml"))
  cps_mort <- ipumsr::read_ipums_micro(ddi)
  
  eitc <- cps_mort %>%
    # limit to households with two children
    dplyr::filter(NCHILD == 2) %>%
    # limit to single and head of household filers
    dplyr::filter(FILESTAT %in% c(4, 5)) %>%
    # drop dependents
    dplyr::filter(DEPSTAT == 0)
  
  # calculate earned income
  eitc <- eitc %>%
    dplyr::mutate(
      earned_income = 
        INCWAGE +
        INCBUS +
        INCFARM
    )
  
  # inflation adjust to 2014 dollars
  eitc <- eitc %>%
    dplyr::mutate(
      earned_income_2014 =
        dplyr::case_when(
          YEAR == 2010 ~ earned_income * 236.715 / 218.076166666666,
          YEAR == 2011 ~ earned_income * 236.715 / 224.923,
          YEAR == 2012 ~ earned_income * 236.715 / 229.586083333333, 
          YEAR == 2013 ~ earned_income * 236.715 / 232.95175,
          YEAR == 2014 ~ earned_income
        )
    ) %>%
    dplyr::select(earned_income, earned_income_2014)
  
  eitc <- eitc %>%
    dplyr::filter(
      earned_income_2014 > 0,
      earned_income_2014 < 30000
    ) 

  return(eitc)
    
}
