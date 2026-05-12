# =========================
# BLOCK 2 — FDA recalls -> state-month outcomes -> join NOAA -> join pathogen cost
#   Inputs:
#     - food-enforcement-0001-of-0001.json
#     - pathogen_costs.csv  (your USDA/ERS pathogen cost table)
#     - climate_features (created by Block 1)
#   Output objects:
#     - fda_state_month
#     - ml_panel
#     - ml_panel_cost
#   Output files:
#     - fda_food_recalls_state_month.csv
#     - ml_panel_noaa_x_fda_recalls.csv
#     - ml_panel_noaa_x_fda_recalls_x_cost.csv
# =========================
library(jsonlite)
library(dplyr)
library(lubridate)
library(stringr)
library(readr)
library(tidyr)
library(janitor)

# ---- Load FDA food recalls ----
recall_raw <- fromJSON("food-enforcement-0001-of-0001.json", flatten = TRUE)

fda <- recall_raw$results %>%
  as_tibble() %>%
  filter(country == "United States", product_type == "Food") %>%
  mutate(
    recall_date = ymd(recall_initiation_date),
    year  = year(recall_date),
    month = month(recall_date),
    date  = floor_date(recall_date, "month"),
    state_abbr = toupper(trimws(state)),
    classification = str_squish(classification),
    reason_lc = str_to_lower(reason_for_recall)
  ) %>%
  filter(!is.na(recall_date), !is.na(state_abbr), nchar(state_abbr) == 2)

# ---- Extract pathogen to bridge to cost table ----
fda <- fda %>%
  mutate(
    pathogen = case_when(
      str_detect(reason_lc, regex("listeria\\s+m(ono)?c(ytogenes)?|listeria\\s*m\\.", ignore_case = TRUE)) ~ 
        "Listeria monocytogenes (total)",
      str_detect(reason_lc, regex("salmonella(\\s+spp\\.?|\\s+species|\\s+enteritidis|\\s+typhimurium)?", ignore_case = TRUE)) ~
        "Salmonella spp., nontyphoidal",
      str_detect(reason_lc, regex("campylobacter(\\s+(jejuni|coli))?", ignore_case = TRUE)) ~
        "Campylobacter spp.",
      str_detect(reason_lc, regex("stec\\s*o157|shiga\\s*-?toxin\\s*producing\\s*e\\.?\\s*coli|e\\.?\\s*coli\\s*o157", ignore_case = TRUE)) ~
        "STEC O157",
      str_detect(reason_lc, regex("non-?o157|stec", ignore_case = TRUE)) ~
        "STEC non-O157",
      str_detect(reason_lc, regex("etec", ignore_case = TRUE)) ~
        "ETEC",
      str_detect(reason_lc, regex("clostridium\\s+botulinum|c\\.\\s*botulinum", ignore_case = TRUE)) ~
        "Clostridium botulinum",
      str_detect(reason_lc, regex("clostridium\\s+perfringens|c\\.\\s*perfringens", ignore_case = TRUE)) ~
        "Clostridium perfringens",
      str_detect(reason_lc, regex("staphylococcus\\s+ aureus|s\\.\\s*aureus", ignore_case = TRUE)) ~
        "Staphylococcus aureus",
      str_detect(reason_lc, regex("shigella(\\s+\\w+)?", ignore_case = TRUE)) ~
        "Shigella spp.",
      str_detect(reason_lc, regex("bacillus\\s+cereus|b\\.\\s*cereus", ignore_case = TRUE)) ~
        "Bacillus cereus",
      str_detect(reason_lc, regex("brucella(\\s+spp\\.?|\\s+species)?", ignore_case = TRUE)) ~
        "Brucella spp.",
      str_detect(reason_lc, regex("mycobacterium\\s+bovis|m\\.\\s*bovis", ignore_case = TRUE)) ~
        "Mycobacterium bovis",
      TRUE ~ NA_character_
    )
  )

# ---- Outcomes: recalls per state-month (+ severity) ----
fda_state_month <- fda %>%
  group_by(state_abbr, year, month, date) %>%
  summarise(
    recall_count = n(),
    class_I   = sum(classification == "Class I",   na.rm = TRUE),
    class_II  = sum(classification == "Class II",  na.rm = TRUE),
    class_III = sum(classification == "Class III", na.rm = TRUE),
    .groups = "drop"
  )

write_csv(fda_state_month, "fda_food_recalls_state_month.csv")

# ---- Join to NOAA features (keep zero-recall months) ----
ml_panel <- climate_features %>%
  transmute(
    state_abbr, year, month, date,
    temp_f, pcpn_in, temp_anom, pcpn_anom,
    temp_f_lag1, pcpn_in_lag1, temp_f_roll3, pcpn_in_roll3,
    month_sin, month_cos
  ) %>%
  left_join(fda_state_month, by = c("state_abbr","year","month","date")) %>%
  mutate(
    recall_count = replace_na(recall_count, 0L),
    class_I      = replace_na(class_I, 0L),
    class_II     = replace_na(class_II, 0L),
    class_III    = replace_na(class_III, 0L)
  )

write_csv(ml_panel, "ml_panel_noaa_x_fda_recalls.csv")

# ---- Load pathogen cost table (standardize columns) ----
cost <- read_csv("foodborne-illness-cost-2023.csv", show_col_types = FALSE) %>%
  clean_names() %>%
  transmute(
    pathogen = str_squish(pathogen),
    mean_cases = mean_number_of_cases,
    mean_total_cost_m = mean_total_cost_millions,
    mean_per_case_cost = mean_per_case_cost
  )

# ---- Count recalls by state-month-pathogen -> join cost -> roll up cost index ----
fda_state_month_cost <- fda %>%
  filter(!is.na(pathogen)) %>%
  group_by(state_abbr, year, month, date, pathogen) %>%
  summarise(recall_count_pathogen = n(), .groups = "drop") %>%
  left_join(cost, by = "pathogen") %>%
  mutate(pathogen_weighted_cost_index = recall_count_pathogen * mean_per_case_cost) %>%
  group_by(state_abbr, year, month, date) %>%
  summarise(
    pathogen_recalls = sum(recall_count_pathogen, na.rm = TRUE),
    pathogen_weighted_cost_index = sum(pathogen_weighted_cost_index, na.rm = TRUE),
    .groups = "drop"
  )

# ---- Join cost features into modeling panel ----
ml_panel_cost <- ml_panel %>%
  left_join(fda_state_month_cost, by = c("state_abbr","year","month","date")) %>%
  mutate(
    pathogen_recalls = replace_na(pathogen_recalls, 0L),
    pathogen_weighted_cost_index = replace_na(pathogen_weighted_cost_index, 0)
  )

ml_panel_cost <- ml_panel_cost %>%
  filter(date >= as.Date("2008-02-01"))

# check
ml_panel_cost %>%
  filter(pathogen_weighted_cost_index > 0) %>%
  arrange(desc(pathogen_weighted_cost_index)) %>%
  head(10)

ml_panel_cost %>%
  summarise(
    total_recalls = sum(recall_count),
    total_pathogen_recalls = sum(pathogen_recalls),
    share_pathogen = total_pathogen_recalls / total_recalls
  )

# add log of pathogen cost
ml_panel_cost <- ml_panel_cost %>%
  mutate(
    log_pathogen_cost_index = log1p(pathogen_weighted_cost_index)
  )

write_csv(ml_panel_cost, "ml_panel_noaa_x_fda_recalls_x_cost.csv")
