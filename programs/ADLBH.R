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

#library(prettyR)

# import the lab xpt file (sdtm.lb.xpt)
# setwd("/cloud/project/sdtm")
# lb<-read_xpt("lb.xpt")
lb <- read_xpt("sdtm/lb.xpt")

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

lb <- convert_blanks_to_na(lb)

# lbh<-filter(lb,lbcat=="HEMATOLOGY")
lbh <- lb %>%
  filter(LBCAT == "HEMATOLOGY")

# import the lab xpt file (adam.adsl.xpt)
#setwd("/cloud/project/adam")
#adsl<-read_xpt("adsl.xpt")
adsl <- read_xpt("adam/adsl.xpt")

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
adlbh <- lbh %>%
  # Join ADSL with LB (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  )
#%>%
#   ## Calculate ADT, ADY ----
# derive_vars_dt(
#   new_vars_prefix = "A",
#   dtc = LBDTC
# ) %>%
#   derive_vars_dy(reference_date = TRTSDT, source_vars = vars(ADT))

adlbh <- adlbh %>%
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
  AVAL  = LBSTRESN,
  AVALC = LBSTRESC,
  ANRLO = LBSTNRLO,
  ANRHI = LBSTNRHI,
  A1LO  = LBSTNRLO,
  A1HI  = LBSTNRHI,
  ABLFL = LBBLFL,
  ADT   = LBDTC,
  ADY   = LBDY
)

# RACEN management
adlbh <- adlbh %>% mutate(
  RACEN = case_when(
    RACE %in% c("AMERICAN INDIAN OR ALASKA NATIVE") ~ 6,
    RACE %in% c("ASIAN") ~ 3,
    RACE %in% c("BLACK OR AFRICAN AMERICAN") ~ 2,
    RACE %in% c("WHITE") ~ 1

  )
)

# Derive Absolute values from fractional Differentials using WBC
# Only derive where absolute values do not already exist
# Need to populate ANRLO and ANRHI for newly created records

# adlbh <- adlbh %>%
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

adlbh <- subset(adlbh,!(VISIT %in% c("AMBUL ECG REMOVAL","RETRIEVAL")))

adlbh <- adlbh %>%
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

## BASE variable
# adlbh_bas <- adlbh %>%
#   mutate(BASE = if_else(ABLFL %in% "Y",adlbh$AVAL,0))
# adlbh_bas_t <- adlbh_bas %>% select(USUBJID,PARAMCD,AVISIT,AVAL,AVALC,ABLFL,BASE)

adlbh <- adlbh %>%
  ## Calculate ONTRTFL ----
derive_var_ontrtfl(
  start_date = ADT,
  ref_start_date = TRTSDT,
  ref_end_date = TRTEDT,
  filter_pre_timepoint = AVISIT == "Baseline"
)


## Calculate ANRIND : requires the reference ranges ANRLO, ANRHI ----
adlbh <- adlbh %>%
  derive_var_anrind()

# Derive baseline flags ----
adlbh <- adlbh %>%
  # Calculate BASETYPE
  mutate(
    BASETYPE = "LAST"
  )
#%>%
#   # Calculate ABLFL
#   restrict_derivation(
#     derivation = derive_var_extreme_flag,
#     args = params(
#       by_vars = vars(STUDYID, USUBJID, BASETYPE, PARAMCD),
#       order = vars(ADT, VISITNUM, LBSEQ),
#       new_var = ABLFL,
#       mode = "last"
#     ),
#     filter = (!is.na(AVAL) & ADT <= TRTSDT & !is.na(BASETYPE))
#   )

## Derive baseline information ----
adlbh <- adlbh %>%
  # Calculate BASE
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate BASEC
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = AVALC,
    new_var = BASEC
  ) %>%
  # Calculate BNRIND
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD, BASETYPE),
    source_var = ANRIND,
    new_var = BNRIND
  ) %>%
  # Calculate CHG
  derive_var_chg() %>%
  # Calculate PCHG
  derive_var_pchg()

## Calculate R2BASE, R2ANRLO and R2ANRHI ----
adlbh <- adlbh %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = BASE,
    new_var = CHG
  ) %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = A1LO,
    new_var = R2A1LO
  ) %>%
  derive_var_analysis_ratio(
    numer_var = AVAL,
    denom_var = A1HI,
    new_var = R2A1HI
  )


# BR2A1H1 - BR2A1LO  ----
adlbh <- adlbh %>%
  mutate(
    BR2A1HI = BASE/A1LO,
    BR2A1LO = BASE/A1HI
  )


## SHIFT derivation ----
adlbh <- adlbh %>%
  # Derive shift from baseline for analysis indicator
  derive_var_shift(
    new_var = SHIFT1,
    from_var = BNRIND,
    to_var = ANRIND
  )

## Flag variables (ANL01FL) ----
# ANL01FL: Flag last result within an AVISIT for post-baseline records
adlbh <- adlbh %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID, PARAMCD, AVISIT),
      order = vars(ADT, AVAL),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVAL) & ONTRTFL == "Y"
  )

## Get extreme values ----
adlbh <- adlbh %>%
  # get MAXIMUM value
  derive_extreme_records(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    order = vars(ADT, AVISITN,LBSEQ),
    mode = "last",
    filter = (!is.na(AVAL) & AVISITN > 1 & AVISITN < 13),
    set_values_to = vars(
      AVISITN = 99,
      AVISIT = "End of Treatment"
    )
  )

# AEMTMTFL
adlbh <- adlbh %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(USUBJID, PARAMCD),
      order = vars(ADT, VISITNUM, LBSEQ),
      new_var = AEMTMTFL,
      mode = "last"
    ),
    filter = !is.na(AVAL) & VISITNUM < 13 & VISITNUM < 1
  )


# Rename variable to have correct ADLBH
adlbh_l <- adlbh %>% rename(TRTA="TRT01A",
                            TRTAN="TRT01AN",
                            TRTP="TRT01P",
                            TRTPN="TRT01PN"
)

# Verification
# adlbh_fin <- adlbh_l %>%
#   select(STUDYID, SUBJID, USUBJID, TRTP, TRTPN, TRTA, TRTAN, TRTSDT, TRTEDT,
#          AGE, AGEGR1, AGEGR1N, RACE, RACEN, SEX, COMP24FL, DSRAEFL, SAFFL,
#          AVISIT, AVISITN, ADY, ADT, VISIT, VISITNUM, PARAM, PARAMCD, PARAMN,
#          PARCAT1, AVAL, BASE, CHG, A1LO, A1HI, R2A1LO, R2A1HI, BR2A1LO,
#          BR2A1HI, ANL01FL, ALBTRVAL, ANRIND, BNRIND, ABLFL, AENTMTFL, LBSEQ,
#          LBNRIND, LBSTRESN)

adlbh_fin <- adlbh_l %>%
  select(STUDYID, SUBJID, USUBJID, TRTP, TRTPN, TRTA, TRTAN, TRTSDT, TRTEDT,
         AGE, RACE, RACEN, SEX, COMP24FL, DSRAEFL, SAFFL,
         AVISIT, AVISITN, ADY, ADT, VISIT, VISITNUM, PARAM, PARAMCD, PARAMN,
         PARCAT1, AVAL, BASE, CHG, A1LO, A1HI, R2A1LO, R2A1HI, BR2A1LO,
         BR2A1HI, ANL01FL, ANRIND, BNRIND, ABLFL, AENTMTFL, LBSEQ,
         LBNRIND, LBSTRESN)

# Sort order
adlbh_fin2 <- adlbh_fin %>%
  arrange(USUBJID,PARAMCD, AVISITN, LBSEQ)


# Apply label
# setwd("/cloud/project/metadata")
# adlbh_spec <-readxl::read_excel("specs.xlsx", sheet = "Variables")%>%
#   filter(Dataset == "ADLBH")

adlbh_spec <-readxl::read_excel("metadata/specs.xlsx", sheet = "Variables")%>%
  filter(Dataset == "ADLBH")

adlbh_spec <- subset(adlbh_spec, select = c(Dataset, Variable, Label))
adlbh_spec <- rename(adlbh_spec, label=Label)
adlbh_spec <- rename(adlbh_spec, dataset=Dataset)
adlbh_spec <- rename(adlbh_spec, variable=Variable)


adlbh_fin2  <- adlbh_fin2 %>%
  xportr_label(adlbh_spec, domain = "ADLBH") # Assigns variable label from metacore specifications


# Export the HEM lab xpt file (adam.adlbh.xpt)
# setwd("/cloud/project/adam")
# xportr_write(adlbh_fin2, "adlbh.xpt")
xportr_write(adlbh_fin2, "adam/adlbh.xpt")

