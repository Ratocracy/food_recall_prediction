# =========================
# BLOCK 1 — Rebuild NOAA climate_features (state-month) from RAW nClimDiv
# Uses NOAA documentation tables to decode ID, then joins to FIPS.
# =========================

library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(lubridate)
library(slider)
library(tigris)

options(tigris_use_cache = TRUE)

tmp_path  <- "climdiv-tmpcst-v1.0.0-20260205.txt"
pcpn_path <- "climdiv-pcpnst-v1.0.0-20260205.txt"

# --- NOAA STATE CODE TABLE (001–050 states only; exclude 101+ regions) ---
noaa_state_lookup <- tibble::tribble(
  ~noaa_state_code, ~state_name, ~state_abbr,
  "001","Alabama","AL",
  "002","Arizona","AZ",
  "003","Arkansas","AR",
  "004","California","CA",
  "005","Colorado","CO",
  "006","Connecticut","CT",
  "007","Delaware","DE",
  "008","Florida","FL",
  "009","Georgia","GA",
  "010","Idaho","ID",
  "011","Illinois","IL",
  "012","Indiana","IN",
  "013","Iowa","IA",
  "014","Kansas","KS",
  "015","Kentucky","KY",
  "016","Louisiana","LA",
  "017","Maine","ME",
  "018","Maryland","MD",
  "019","Massachusetts","MA",
  "020","Michigan","MI",
  "021","Minnesota","MN",
  "022","Mississippi","MS",
  "023","Missouri","MO",
  "024","Montana","MT",
  "025","Nebraska","NE",
  "026","Nevada","NV",
  "027","New Hampshire","NH",
  "028","New Jersey","NJ",
  "029","New Mexico","NM",
  "030","New York","NY",
  "031","North Carolina","NC",
  "032","North Dakota","ND",
  "033","Ohio","OH",
  "034","Oklahoma","OK",
  "035","Oregon","OR",
  "036","Pennsylvania","PA",
  "037","Rhode Island","RI",
  "038","South Carolina","SC",
  "039","South Dakota","SD",
  "040","Tennessee","TN",
  "041","Texas","TX",
  "042","Utah","UT",
  "043","Vermont","VT",
  "044","Virginia","VA",
  "045","Washington","WA",
  "046","West Virginia","WV",
  "047","Wisconsin","WI",
  "048","Wyoming","WY",
  "050","Alaska","AK"
)

# --- FIPS lookup (abbr -> 2-digit FIPS) ---
# This is safe: we force one row per state.
state_fips_lookup <- tigris::fips_codes %>%
  distinct(state_code, state) %>%
  filter(!is.na(state_code), nchar(state_code) == 2) %>%
  transmute(
    state_fips = state_code,
    state_abbr = str_to_upper(state)
  )

read_climdiv_raw <- function(path, value_name, missing_sentinel) {
  lines <- read_lines(path)
  lines <- lines[nzchar(lines)]
  toks  <- str_split(str_trim(lines), "\\s+")
  
  id <- vapply(toks, `[[`, character(1), 1)
  vals_list <- lapply(toks, function(x) suppressWarnings(as.numeric(x[-1])))
  
  # Ensure we have 12 months
  max_len <- max(vapply(vals_list, length, integer(1)))
  mat <- do.call(rbind, lapply(vals_list, function(v) { length(v) <- max_len; v }))
  if (ncol(mat) < 12) stop("Fewer than 12 numeric columns after ID; check file format.")
  
  months12 <- mat[, 1:12, drop = FALSE]
  colnames(months12) <- sprintf("m%02d", 1:12)
  
  as_tibble(months12) %>%
    mutate(id = str_pad(id, 10, side = "left", pad = "0")) %>%
    mutate(
      noaa_state_code = str_sub(id, 1, 3),  # 001–110
      division        = str_sub(id, 4, 4),  # 0 = statewide avg
      element         = str_sub(id, 5, 6),  # 01 precip, 02 temp
      year            = as.integer(str_sub(id, 7, 10))
    ) %>%
    pivot_longer(starts_with("m"), names_to = "month", values_to = value_name) %>%
    mutate(
      month = as.integer(str_remove(month, "^m")),
      date  = as.Date(sprintf("%04d-%02d-01", year, month)),
      !!value_name := na_if(.data[[value_name]], missing_sentinel)
    ) %>%
    left_join(noaa_state_lookup, by = "noaa_state_code") %>%
    left_join(state_fips_lookup, by = "state_abbr") %>%
    filter(!is.na(state_abbr), !is.na(state_fips))
}

temp_df <- read_climdiv_raw(tmp_path,  "temp_f",  missing_sentinel = -99.90) %>%
  filter(division == "0", element == "02") %>%
  select(state_fips, state_abbr, state_name, year, month, date, temp_f)

pcpn_df <- read_climdiv_raw(pcpn_path, "pcpn_in", missing_sentinel = -9.99) %>%
  filter(division == "0", element == "01") %>%
  select(state_fips, state_abbr, state_name, year, month, date, pcpn_in)

# Use FULL JOIN so missing temp/precip doesn't drop the row
climate_state_month <- temp_df %>%
  full_join(pcpn_df, by = c("state_fips","state_abbr","state_name","year","month","date")) %>%
  arrange(state_fips, date)

stopifnot(nrow(climate_state_month) == nrow(distinct(climate_state_month, state_fips, year, month)))

write_csv(climate_state_month, "climdiv_state_month_temp_pcpn.csv")

# ---- Feature engineering ----
climate_features <- climate_state_month %>%
  arrange(state_fips, date) %>%
  group_by(state_fips) %>%
  mutate(
    temp_f_lag1   = lag(temp_f, 1),
    pcpn_in_lag1  = lag(pcpn_in, 1),
    temp_f_roll3  = slide_dbl(temp_f,  mean, .before = 2, .complete = TRUE, na.rm = TRUE),
    pcpn_in_roll3 = slide_dbl(pcpn_in, mean, .before = 2, .complete = TRUE, na.rm = TRUE),
    month_sin = sin(2*pi*month/12),
    month_cos = cos(2*pi*month/12)
  ) %>%
  ungroup()

climatology_1991_2020 <- climate_features %>%
  filter(year >= 1991, year <= 2020) %>%
  group_by(state_fips, month) %>%
  summarise(
    temp_clim_1991_2020 = mean(temp_f,  na.rm = TRUE),
    pcpn_clim_1991_2020 = mean(pcpn_in, na.rm = TRUE),
    .groups = "drop"
  )

climate_features <- climate_features %>%
  left_join(climatology_1991_2020, by = c("state_fips","month")) %>%
  mutate(
    temp_anom = temp_f - temp_clim_1991_2020,
    pcpn_anom = pcpn_in - pcpn_clim_1991_2020
  )
