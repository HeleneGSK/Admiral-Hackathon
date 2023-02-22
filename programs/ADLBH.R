#Libraries
library(admiral)
library(admiraldev)
library(haven)
library(dplyr)
library(lubridate)
library(stringr)
library(xportr)
library(metacore)
library(metatools)

library(prettyR)

# import the lab xpt file (sdtm.lb.xpt)
setwd("/cloud/project/sdtm")
lb<-read_xpt("lb.xpt")

setwd("/cloud/project/adam")
adsl<-read_xpt("adsl.xpt")

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

lb <- convert_blanks_to_na(lb)

unique(lb$LBCAT)

# lbh<-filter(lb,lbcat=="HEMATOLOGY")
lbh <- lb %>%
  filter(LBCAT == "HEMATOLOGY")

# Look-up tables ----

# Assign PARAMCD, PARAM, and PARAMN
param_lookup <- tibble::tribble(
  ~LBTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "ALB", "ALB", "Albumin (g/L)", 1,
  "ALP", "ALKPH", "Alkaline Phosphatase (U/L)", 2,
  "ALT", "ALT", "Alanine Aminotransferase (U/L)", 3,
  "ANISO", "ANISO", "Anisocytes", 4,
  "AST", "AST", "Aspartate Aminotransferase (U/L)", 5,
  "BASO", "BASO", "Basophils (10^9/L)", 6,
  "BASOLE", "BASOLE", "Basophils/Leukocytes (FRACTION)", 7,
  "BILI", "BILI", "Bilirubin (umol/L)", 8,
  "BUN", "BUN", "Blood Urea Nitrogen (mmol/L)", 9,
  "CA", "CA", "Calcium (mmol/L)", 10,
  "CHOL", "CHOLES", "Cholesterol (mmol/L)", 11,
  "CK", "CK", "Creatinine Kinase (U/L)", 12,
  "CL", "CL", "Chloride (mmol/L)", 13,
  "COLOR", "COLOR", "Color", 14,
  "CREAT", "CREAT", "Creatinine (umol/L)", 15,
  "EOS", "EOS", "Eosinophils (10^9/L)", 16,
  "EOSLE", "EOSLE", "Eosinophils/Leukocytes (FRACTION)", 17,
  "GGT", "GGT", "Gamma Glutamyl Transferase (U/L)", 18,
  "GLUC", "GLUC", "Glucose (mmol/L)", 19,
  "HBA1C", "HBA1C", "Hemoglobin A1C (1)", 20,
  "HCT", "HCT", "Hematocrit (1)", 21,
  "HGB", "HGB", "Hemoglobin (mmol/L)", 22,
  "K", "POTAS", "Potassium (mmol/L)", 23,
  "KETONES", "KETON", "Ketones", 24,
  "LYM", "LYMPH", "Lymphocytes (10^9/L)", 25,
  "LYMLE", "LYMPHLE", "Lymphocytes/Leukocytes (FRACTION)", 26,
  "MACROCY", "MACROC", "Macrocytes", 27,
  "MCH", "MCH", "Ery. Mean Corpuscular Hemoglobin (fmol(Fe))", 28,
  "MCHC", "MCHC", "Ery. Mean Corpuscular HGB Concentration (mmol/L)", 29,
  "MCV", "MCV", "Ery. Mean Corpuscular Volume (f/L)", 30,
  "MICROCY", "MICROC", "Microcytes", 31,
  "MONO", "MONO", "Monocytes (10^9/L)", 32,
  "MONOLE", "MONOLE", "Monocytes/Leukocytes (FRACTION)", 33,
  "PH", "PH", "pH", 34,
  "PHOS", "PHOS", "Phosphate (mmol/L)", 35,
  "PLAT", "PLAT", "Platelet (10^9/L)", 36,
  "POIKILO", "POIKIL", "Poikilocytes", 37,
  "POLYCHR", "POLYCH", "Polychromasia", 38,
  "PROT", "PROT", "Protein (g/L)", 39,
  "RBC", "RBC", "Erythrocytes (TI/L)", 40,
  "SODIUM", "SODIUM", "Sodium (mmol/L)", 41,
  "SPGRAV", "SPGRAV", "Specific Gravity", 42,
  "TSH", "TSH", "Thyrotropin (mU/L)", 43,
  "URATE", "URATE", "Urate (umol/L)", 44,
  "UROBIL", "UROBIL", "Urobilinogen", 45,
  "VITB12", "VITB12", "Vitamin B12 (pmol/L)", 46,
  "WBC", "WBC", "Leukocytes (10^9/L)", 47
)

# Derivations ----

# Get list of ADSL vars required for derivations
# adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P, AGE, AGEGR1, AGEGR1N,
#                   COMP24FL, DSRAEFL, RACE, RACEN, SAFFL, SEX, SUBJID,
#                   TRT01AN, TRT01PN)
adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P, AGE,
                  COMP24FL, DSRAEFL, RACE, SAFFL, SEX, SUBJID,
                  TRT01AN, TRT01PN)
adlb <- lbh %>%
  # Join ADSL with LB (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  ) %>%
  ## Calculate ADT, ADY ----
derive_vars_dt(
  new_vars_prefix = "A",
  dtc = LBDTC
) %>%
  derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

adlb <- adlb %>%
  ## Add PARAMCD PARAM and PARAMN - from LOOK-UP table ----
# Replace with PARAMCD lookup function
derive_vars_merged_lookup(
  dataset_add = param_lookup,
  new_vars = vars(PARAMCD, PARAM, PARAMN),
  by_vars = vars(LBTESTCD),
  check_type = "none",
  print_not_mapped = FALSE
) %>%
  ## Calculate PARCAT1 AVAL AVALC ANRLO ANRHI ----
mutate(
  PARCAT1 = 'HEM',
  AVAL = LBSTRESN,
  AVALC = LBSTRESC,
  ANRLO = LBSTNRLO,
  ANRHI = LBSTNRHI
)

# Derive Absolute values from fractional Differentials using WBC
# Only derive where absolute values do not already exist
# Need to populate ANRLO and ANRHI for newly created records

# adlb <- adlb %>%
#   # Derive absolute Basophils
#   derive_param_wbc_abs(
#     by_vars = vars(STUDYID, USUBJID, !!!adsl_vars, DOMAIN, VISIT, VISITNUM, ADT, ADY),
#     set_values_to = vars(
#       PARAMCD = "BASO",
#       PARAM = "Basophils Abs (10^9/L)",
#       PARAMN = 6,
#       DTYPE = "CALCULATION",
#       PARCAT1 = "HEM"
#     ),
#     get_unit_expr = extract_unit(PARAM),
#     diff_code = "BASOLE"
#   ) %>%
#   # Derive absolute Lymphocytes
#   derive_param_wbc_abs(
#     by_vars = vars(STUDYID, USUBJID, !!!adsl_vars, DOMAIN, VISIT, VISITNUM, ADT, ADY),
#     set_values_to = vars(
#       PARAMCD = "LYMPH",
#       PARAM = "Lymphocytes Abs (10^9/L)",
#       PARAMN = 25,
#       DTYPE = "CALCULATION",
#       PARCAT1 = "HEM"
#     ),
#     get_unit_expr = extract_unit(PARAM),
#     diff_code = "LYMPHLE"
#   )


## Get Visit Info ----
# See also the "Visit and Period Variables" vignette
# (https://pharmaverse.github.io/admiral/articles/visits_periods.html#visits)
adlb <- adlb %>%
  # Derive Timing
  mutate(
    AVISIT = case_when(
      str_detect(VISIT, "SCREEN") ~ "Baseline",
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = case_when(
      AVISIT == "Baseline" ~ 0,
      !is.na(VISITNUM) ~ VISITNUM
    )
  )
# Verification
# adlb_t <- adlb %>% select(USUBJID,VISIT,VISITNUM, AVISIT, AVISITN)

## Calculate ANRIND : requires the reference ranges ANRLO, ANRHI ----
adlb <- adlb %>%
  derive_var_anrind()











# Export the HEM lab xpt file (adam.adlbh.xpt)
# setwd("/cloud/project/adam")
# xportr_write(adlbh, "adlbh.xpt")
