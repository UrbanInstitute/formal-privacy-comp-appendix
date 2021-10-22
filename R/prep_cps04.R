#' Title
#'
#' @return
#' @export
#'
#' importFrom magrittr "%>%"
#'
#' @examples
prep_cps04 <- function() {

  library(magrittr)
  
  ddi <- ipumsr::read_ipums_ddi(here::here("data", "cps_00013.xml"))
  asec <- ipumsr::read_ipums_micro(ddi)
  
  # remove labels
  asec <- asec %>%
    dplyr::mutate(
      MONTH = as_factor(ipumsr::lbl_clean(MONTH)),
      AGE = as.numeric(ipumsr::lbl_clean(AGE)),
      SEX = as_factor(ipumsr::lbl_clean(SEX)),
      RACE = as_factor(ipumsr::lbl_clean(RACE)),
      #EDUC = as_factor(lbl_clean(EDUC)),
      UHRSWORKLY = as.numeric(ipumsr::lbl_clean(ipumsr::lbl_na_if(UHRSWORKLY, ~.val == 999))),
      INCWAGE = if_else(condition = as.numeric(ipumsr::lbl_clean(INCWAGE)) == 99999999, 
                        true = as.numeric(NA),
                        false = as.numeric(ipumsr::lbl_clean(INCWAGE)))
    )
  
  asec <- asec %>%
    dplyr::mutate(
      INCWAGE95 = case_when(
        YEAR == 1994 ~ INCWAGE * (151.2 / 147.1),
        YEAR == 1995 ~ INCWAGE,
        YEAR == 1996 ~ INCWAGE * (151.2 / 155.5)
      )
    )
  
  ddi_educ <- ipumsr::read_ipums_ddi(here::here("data", "cps_00012.xml"))
  educ <- ipumsr::read_ipums_micro(ddi_educ)
  
  educ <- educ %>%
    dplyr::select(YEAR, WTFINL, AGE, EDUC, HIGRADE)
  
  
  # 1. map Feb 1990 to ASEC 94 (EDUC)
  educ <- educ %>%
    dplyr::mutate(
      cats = dplyr::case_when(
        EDUC ==   1 ~   1,
        EDUC ==   2 ~   2,
        EDUC ==  11 ~  10, # grade 1 to grades 1-4
        EDUC ==  12 ~  10, # grade 1 to grades 1-4
        EDUC ==  13 ~  10, # grade 1 to grades 1-4
        EDUC ==  14 ~  10,
        EDUC ==  21 ~  20, # grade 5 to grades 5-6
        EDUC ==  22 ~  20, # grade 5 to grades 5-6
        EDUC ==  31 ~  30,
        EDUC ==  32 ~  30,
        EDUC ==  40 ~  40,
        EDUC ==  50 ~  50,
        EDUC ==  60 ~  60,
        EDUC ==  72 ~  71, # 
        EDUC ==  73 ~  73,
        EDUC ==  80 ~  81,
        EDUC ==  90 ~  91,
        EDUC == 100 ~  92,
        EDUC == 110 ~ 111,
        EDUC == 121 ~ 122,
        EDUC == 122 ~ 122,
        EDUC == 123 ~ 122,
        EDUC == 124 ~ 122,
        EDUC == 125 ~ 122
      )
    ) 
  
  educ <- educ %>%
    dplyr::mutate(
      years_of_educ = dplyr::case_when(
        HIGRADE ==  0 ~ 0,
        HIGRADE ==  2 ~ 0,
        HIGRADE <  50 ~ 1,
        HIGRADE <  60 ~ 2,
        HIGRADE <  70 ~ 3,
        HIGRADE <  80 ~ 4,
        HIGRADE <  90 ~ 5,
        HIGRADE < 100 ~ 6,
        HIGRADE < 110 ~ 7,
        HIGRADE < 120 ~ 8,
        HIGRADE < 130 ~ 9,
        HIGRADE < 140 ~ 10,
        HIGRADE < 150 ~ 11,
        HIGRADE < 160 ~ 12,
        HIGRADE < 170 ~ 13,
        HIGRADE < 180 ~ 14,
        HIGRADE < 190 ~ 15,
        HIGRADE < 200 ~ 16,
        HIGRADE < 210 ~ 17,
        HIGRADE == 210 ~ 18,
        TRUE ~ as.numeric(NA)
      )
    )
  
  
  # create lookup table to apply to 1994-1996 ASEC
  # lookup_table <- educ %>%
  #   filter(AGE >= 15) %>%
  #   group_by(cats) %>%
  #   summarize(years_of_educ = weighted.mean(x = years_of_educ, w = WTFINL))
  
  # use library(srvyr) to account for the survey design before calculating the 
  # weighted mean
  educ_design_survey <- educ %>%
    dplyr::select(-EDUC, -HIGRADE) %>%
    dplyr::mutate(fpc = sum(WTFINL)) %>%
    srvyr::as_survey_design(ids = 1, fpc = fpc, weights = WTFINL)
  
  lookup_table <- educ_design_survey %>%
    dplyr::filter(AGE >= 15) %>%  
    dplyr::group_by(cats) %>%
    dplyr::summarize(years_of_educ = srvyr::survey_mean(years_of_educ, vartype = "ci")) %>%
    dplyr::select(cats, years_of_educ)
  
  asec <- asec %>%
    dplyr::mutate(EDUC = ifelse(EDUC %in% c(123, 124, 125), 122, EDUC)) %>%
    dplyr::left_join(lookup_table, by = c("EDUC" = "cats"))
  
  asec <- asec %>%
    dplyr::mutate(potential_experience = pmax(0, AGE - years_of_educ - 6))
  
  asec <- asec %>%
    dplyr::mutate(hourly_wage = case_when(
      INCWAGE95 == 0 ~ 0,
      TRUE ~ INCWAGE95 / (WKSWORK1 * UHRSWORKLY)
    )
    )
  
  asec <- asec %>%
    dplyr::mutate(non_white = as.numeric(!RACE == "White"))
  
  asec_subset <- asec %>%
    dplyr::filter(
      AGE >= 16,
      AGE <= 66
    ) %>%
    dplyr::filter(
      hourly_wage > 2,
      hourly_wage < 150
    ) %>%
    dplyr::filter(INCWAGE > 0)
  
  asec_subset <- asec_subset %>%
    dplyr::mutate(
      log_INCWAGE = log(INCWAGE),
      potential_experience_squared = potential_experience ^ 2,
      potential_experience_cubed = potential_experience ^ 3,
      non_white = factor(non_white)
    ) %>%
    dplyr::select(
      SEX,
      log_INCWAGE,
      years_of_educ,
      potential_experience,
      potential_experience_squared,
      potential_experience_cubed,
      non_white
    )
  
  asec_subset_male <- asec_subset %>%
    dplyr::filter(SEX == "Male") %>%
    dplyr::select(-SEX)
  
  asec_subset_female <- asec_subset %>%
    dplyr::filter(SEX == "Female") %>%
    dplyr::select(-SEX)
  
  return(list(female = asec_subset_female,
              male = asec_subset_male))
    
}
