library(haven)
library(admiral)
library(dplyr)
library(tidyr)
library(metacore)
library(metatools)
library(xportr)
library("stringr")

lb <- read_xpt("sdtm/lb.xpt")
lb <- convert_blanks_to_na(lb)
lb1 <- lb %>%
 filter(LBCAT == "CHEMISTRY")


param_lookup <- tibble::tribble(
  ~LBTESTCD, ~PARAMCD, ~PARAM, ~PARAMN,
  "ALB", "ALB", "Albumin (g/L)", 1,
  "ALP", "ALP", "Alkaline Phosphatase (U/L)", 2,
  "ALT", "ALT", "Alanine Aminotransferase (U/L)", 3,
  "ANISO", "ANISO", "Anisocytes", 4,
  "AST", "AST", "Aspartate Aminotransferase (U/L)", 5,
  "BASO", "BASO", "Basophils (10^9/L)", 6,
  "BASOLE", "BASOLE", "Basophils/Leukocytes (FRACTION)", 7,
  "BILI", "BILI", "Bilirubin (umol/L)", 8,
  "BUN", "BUN", "Blood Urea Nitrogen (mmol/L)", 9,
  "CA", "CA", "Calcium (mmol/L)", 10,
  "CHOL", "CHOL", "Cholesterol (mmol/L)", 11,
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
  "K", "K", "Potassium (mmol/L)", 23,
  "KETONES", "KETONES", "Ketones", 24,
  "LYM", "LYM", "Lymphocytes (10^9/L)", 25,
  "LYMLE", "LYMLE", "Lymphocytes/Leukocytes (FRACTION)", 26,
  "MACROCY", "MACROCY", "Macrocytes", 27,
  "MCH", "MCH", "Ery. Mean Corpuscular Hemoglobin (fmol(Fe))", 28,
  "MCHC", "MCHC", "Ery. Mean Corpuscular HGB Concentration (mmol/L)", 29,
  "MCV", "MCV", "Ery. Mean Corpuscular Volume (f/L)", 30,
  "MICROCY", "MICROCY", "Microcytes", 31,
  "MONO", "MONO", "Monocytes (10^9/L)", 32,
  "MONOLE", "MONOLE", "Monocytes/Leukocytes (FRACTION)", 33,
  "PH", "PH", "pH", 34,
  "PHOS", "PHOS", "Phosphate (mmol/L)", 35,
  "PLAT", "PLAT", "Platelet (10^9/L)", 36,
  "POIKILO", "POIKILO", "Poikilocytes", 37,
  "POLYCHR", "POLYCHR", "Polychromasia", 38,
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

adsl <- read_xpt("adam/adsl.xpt")


# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P, AGE, AGEGR1, AGEGR1N, COMP24FL,
                  DSRAEFL, RACE, RACEN,SAFFL, SEX, SUBJID, TRT01AN, TRT01PN)
# adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P, AGE, COMP24FL,
#                   DSRAEFL, RACE,SAFFL, SEX, SUBJID, TRT01AN, TRT01PN)

adlb <- lb1 %>%
  # Join ADSL with LB (need TRTSDT for ADY derivation)
  derive_vars_merged(
    dataset_add = adsl,
    new_vars = adsl_vars,
    by_vars = vars(STUDYID, USUBJID)
  )%>%
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
)%>%
  ## Calculate PARCAT1 AVAL AVALC ANRLO ANRHI ----
mutate(
  PARCAT1 = 'CHEM',
  AVAL = LBSTRESN,
  AVALC = LBSTRESC,
  A1LO = LBSTNRLO,
  A1HI = LBSTNRHI,
  ANRLO = 0.5*LBSTNRLO,
  ANRHI = 1.5*LBSTNRHI,
  TRTP = TRT01P,
  TRTPN = TRT01PN,
  TRTA = TRT01A,
  TRTAN = TRT01AN
)

adlb <- adlb %>%
  mutate(ALBTRVAL= max(((1.5*A1HI)-LBSTRESN) , (LBSTRESN-(.5*A1LO)))
)

adlb <- adlb %>%
  mutate(
    RACEN = case_when(
      RACE == 'AMERICAN INDIAN OR ALASKA NATIVE' ~ 6,
      RACE == 'ASIAN' ~ 3,
      RACE == 'BLACK OR AFRICAN AMERICAN' ~ 2,
      RACE == 'WHITE' ~ 1
    )
  )

adlb <- adlb %>%
  # Derive Timing
  mutate(
    AVISIT = case_when(
      str_detect(VISIT, "SCREENING 1") ~ "Baseline",
      !is.na(VISIT) ~ str_to_title(VISIT),
      TRUE ~ NA_character_
    ),
    AVISITN = case_when(
      AVISIT == "Baseline" ~ 0,
      !is.na(VISITNUM) ~ VISITNUM
    )
  )
adlb <- adlb %>%
  filter(VISITNUM != 6,
         VISITNUM != 201)


## Calculate ANRIND : requires the reference ranges ANRLO, ANRHI ----
adlb <- adlb %>%
  derive_var_anrind()

## Derive baseline flags ----
# adlb <- adlb %>%
#   # Calculate BASETYPE
#   mutate(
#     BASETYPE = "LAST"
#   ) %>%
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
adlb <- adlb %>%
  mutate(ABLFL = if_else(VISITNUM==1, 'Y', ''))

adlb <- adlb %>%
  mutate(R2A1HI = AVAL/A1LO,
         R2A1LO = AVAL/A1HI)

## Derive baseline information ----
adlb <- adlb %>%
  # Calculate BASE
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    source_var = AVAL,
    new_var = BASE
  ) %>%
  # Calculate BASEC
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    source_var = AVALC,
    new_var = BASEC
  ) %>%
  # Calculate BNRIND
  derive_var_base(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    source_var = ANRIND,
    new_var = BNRIND
  ) %>%
  # Calculate CHG
  derive_var_chg() %>%
  # Calculate PCHG
  derive_var_pchg()

adlb <- adlb %>%
  mutate(BR2A1HI = BASE/A1LO,
         BR2A1LO = BASE/A1HI)

adlb <- adlb %>%

  # Calculate AENTMTFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(STUDYID, USUBJID,  PARAMCD),
      order = vars(ADT, VISITNUM, LBSEQ),
      new_var = AENTMTFL,
      mode = "last"
    ),
    filter = (!is.na(AVAL) & VISITNUM <= 12 &  VISITNUM >= 2)
  )

adlb <- adlb %>%

  # Calculate ANL01FL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(STUDYID, USUBJID,  PARAMCD),
      order = vars(ADT, VISITNUM, LBSEQ),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = (!is.na(AVAL) & VISITNUM <= 12 &  VISITNUM >= 2)
  )

adlb <- adlb %>%

  # get MAXIMUM value
  derive_extreme_records(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    order = vars(desc(ALBTRVAL), ADT, AVISITN),
    mode = "first",
    filter = (!is.na(ALBTRVAL) & VISITNUM >= 2),
    set_values_to = vars(
      ANL01FL = 'Y'
    )
  )

adlb <- adlb %>%
  # get MAXIMUM value
  derive_extreme_records(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    order = vars( ADT, AVISITN),
    mode = "last",
    filter = (!is.na(AVAL) & VISITNUM >= 2 & AVISITN <= 12),
    set_values_to = vars(
      AVISITN = 99,
      AVISIT = "End of Treatment",
    )
  )

adlbc <- select(adlb,
                STUDYID,
                SUBJID,
                USUBJID,
                TRTP,
                TRTPN,
                TRTA,
                TRTAN,
                TRTSDT,
                TRTEDT,
                AGE,
                AGEGR1,
                AGEGR1N,
                RACE,
                RACEN,
                SEX,
                COMP24FL,
                DSRAEFL,
                SAFFL,
                AVISIT,
                AVISITN,
                ADY,
                ADT,
                VISIT,
                VISITNUM,
                PARAM,
                PARAMCD,
                PARAMN,
                PARCAT1,
                AVAL,
                BASE,
                CHG,
                A1LO,
                A1HI,
                R2A1LO,
                BR2A1LO,
                BR2A1HI,
                ANL01FL,
                ALBTRVAL,
                ANRIND,
                BNRIND,
                ABLFL,
                AENTMTFL,
                LBSEQ,
                LBNRIND,
                LBSTRESN)
 xportr_write(adlbc, "adam/adlbc.xpt")
