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
  "ALB", "ALB", "Albumin (g/L)", 33,
  "ALP", "ALP", "Alkaline Phosphatase (U/L)", 22,
  "ALT", "ALT", "Alanine Aminotransferase (U/L)", 24,
  "ANISO", "ANISO", "Anisocytes", 4,
  "AST", "AST", "Aspartate Aminotransferase (U/L)",25,
  "BASO", "BASO", "Basophils (10^9/L)", 6,
  "BASOLE", "BASOLE", "Basophils/Leukocytes (FRACTION)", 7,
  "BILI", "BILI", "Bilirubin (umol/L)", 21,
  "BUN", "BUN", "Blood Urea Nitrogen (mmol/L)", 26,
  "CA", "CA", "Calcium (mmol/L)", 30,
  "CHOL", "CHOL", "Cholesterol (mmol/L)", 34,
  "CK", "CK", "Creatine Kinase (U/L)", 35,
  "CL", "CL", "Chloride (mmol/L)", 20,
  "COLOR", "COLOR", "Color", 14,
  "CREAT", "CREAT", "Creatinine (umol/L)", 27,
  "EOS", "EOS", "Eosinophils (10^9/L)", 16,
  "EOSLE", "EOSLE", "Eosinophils/Leukocytes (FRACTION)", 17,
  "GGT", "GGT", "Gamma Glutamyl Transferase (U/L)", 23,
  "GLUC", "GLUC", "Glucose (mmol/L)", 31,
  "HBA1C", "HBA1C", "Hemoglobin A1C (1)", 13,
  "HCT", "HCT", "Hematocrit (1)", 8,
  "HGB", "HGB", "Hemoglobin (mmol/L)", 2,
  "K", "K", "Potassium (mmol/L)", 19,
  "KETONES", "KETONES", "Ketones", 3,
  "LYM", "LYM", "Lymphocytes (10^9/L)", 5,
  "LYMLE", "LYMLE", "Lymphocytes/Leukocytes (FRACTION)", 9,
  "MACROCY", "MACROCY", "Macrocytes", 15,
  "MCH", "MCH", "Ery. Mean Corpuscular Hemoglobin (fmol(Fe))", 44,
  "MCHC", "MCHC", "Ery. Mean Corpuscular HGB Concentration (mmol/L)", 12,
  "MCV", "MCV", "Ery. Mean Corpuscular Volume (f/L)", 10,
  "MICROCY", "MICROCY", "Microcytes", 41,
  "MONO", "MONO", "Monocytes (10^9/L)", 39,
  "MONOLE", "MONOLE", "Monocytes/Leukocytes (FRACTION)", 1,
  "PH", "PH", "pH", 11,
  "PHOS", "PHOS", "Phosphate (mmol/L)", 29,
  "PLAT", "PLAT", "Platelet (10^9/L)", 36,
  "POIKILO", "POIKILO", "Poikilocytes", 37,
  "POLYCHR", "POLYCHR", "Polychromasia", 38,
  "PROT", "PROT", "Protein (g/L)", 32,
  "RBC", "RBC", "Erythrocytes (TI/L)", 40,
  "SODIUM", "SODIUM", "Sodium (mmol/L)", 18,
  "SPGRAV", "SPGRAV", "Specific Gravity", 42,
  "TSH", "TSH", "Thyrotropin (mU/L)", 43,
  "URATE", "URATE", "Urate (umol/L)", 28,
  "UROBIL", "UROBIL", "Urobilinogen", 45,
  "VITB12", "VITB12", "Vitamin B12 (pmol/L)", 46,
  "WBC", "WBC", "Leukocytes (10^9/L)", 47
)

adsl <- read_xpt("adam/adsl.xpt")


# Derivations ----

# Get list of ADSL vars required for derivations
adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P, AGE, AGEGR1, AGEGR1N, COMP24FL,
                  DSRAEFL, RACE, RACEN,SAFFL, SEX, SUBJID, TRT01AN, TRT01PN)
# adsl_vars <- vars(TRTSDT, TRTEDT, TRT01A, TRT01P, AGE, COMP24FL,  DSRAEFL, RACE,SAFFL, SEX, SUBJID, TRT01AN, TRT01PN)



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
  TRTAN = TRT01AN,
  ADT = as.numeric(ADT-as.Date('1960-01-01') ),
  TRTEDT = as.numeric(TRTEDT-as.Date('1960-01-01') ),
  TRTSDT= as.numeric(TRTSDT-as.Date('1960-01-01') ),
  AVISIT= VISIT,
  AVISITN= VISITNUM
)


adlb <- adlb %>%
  mutate(ALBTRVAL1= ((1.5*A1HI)-LBSTRESN)
  )
adlb <- adlb %>%
  mutate(ALBTRVAL2=  (LBSTRESN-(.5*A1LO))
  )
adlb <- adlb %>%
  mutate(ALBTRVAL = if_else(ALBTRVAL1 >ALBTRVAL2 , ALBTRVAL1, ALBTRVAL2))


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
  mutate(
    AVISIT = case_when(
      VISITNUM == 1 ~ 'Baseline',
      VISITNUM == 9 ~ 'Week 12',
      VISITNUM == 10 ~ 'Week 16',
      VISITNUM == 4 ~ 'Week 2',
      VISITNUM == 5 ~ 'Week 4',
      VISITNUM == 7 ~ 'Week 6',
      VISITNUM == 8 ~ 'Week 8',
      VISITNUM == 11 ~ 'Week 20',
      VISITNUM == 12 ~ 'Week 24',
      VISITNUM == 13 ~ 'Week 26',

    )
  )

adlb <- adlb %>%
  mutate(
    AVISITN = case_when(
      VISITNUM == 1 ~ 0,
      VISITNUM == 9 ~ 12,
      VISITNUM == 10 ~ 16,
      VISITNUM == 4 ~ 2,
      VISITNUM == 5 ~ 4,
      VISITNUM == 7 ~ 6,
      VISITNUM == 8 ~ 8,
      VISITNUM == 11 ~ 20,
      VISITNUM == 12 ~ 24,
      VISITNUM == 13 ~ 26,

    )
  )
adlb <- adlb %>%
  filter(VISITNUM != 6,
         VISITNUM != 201)


## Calculate ANRIND : requires the reference ranges ANRLO, ANRHI ----
# adlb <- adlb %>%
#   derive_var_anrind()

# if the AVAL is > 1.5*A1HI then set to 'H',
# else if AVAL is <.5*A1LO then set to 'L',
# else if AVAL between (.5*A1LO and 1.5*A1HI) then set to 'N',

adlb <- adlb %>%
  mutate(
    ANRIND = case_when(
      AVAL > 1.5*A1HI ~ 'H',
      AVAL < 0.5*A1LO  ~ 'L',
      TRUE ~'N'

    )
  )


# adlb <- adlb %>%
#   mutate(
#     ANRIND = case_when(
#       ANRIND == 'NORMAL' ~ 'N'
#
#     )
#   )
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
  mutate(
    SCHED = case_when(
      VISITNUM == 1 ~ 'Y',
      VISITNUM == 9 ~ 'Y',
      VISITNUM == 10 ~ 'Y',
      VISITNUM == 4 ~ 'Y',
      VISITNUM == 5 ~'Y',
      VISITNUM == 7 ~ 'Y',
      VISITNUM == 8 ~ 'Y',
      VISITNUM == 11 ~ 'Y',
      VISITNUM == 12 ~ 'Y',
      VISITNUM == 13 ~ 'Y',

    )
  )


adlb <- adlb %>%
  mutate(R2A1HI = AVAL/A1HI,
         R2A1LO = AVAL/A1LO)

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
  mutate(CHG = ifelse(VISITNUM==1, NA, CHG))

adlb <- adlb %>%
  mutate(BR2A1HI = BASE/A1HI,
         BR2A1LO = BASE/A1LO)




adlb <- adlb %>%

  # Calculate ANL01FL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(STUDYID, USUBJID,  PARAMCD),
      order = vars(ALBTRVAL,desc(VISITNUM) ),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = (!is.na(ALBTRVAL) & VISITNUM <= 12 &  VISITNUM >= 2  & !is.na(SCHED))
  )

# adlb <- adlb %>%
#
#   # get MAXIMUM value
#   derive_extreme_records(
#     by_vars = vars(STUDYID, USUBJID, PARAMCD),
#     order = vars(desc(ALBTRVAL), ADT, AVISITN, LBSEQ),
#     mode = "first",
#     filter = (!is.na(ALBTRVAL) & VISITNUM >= 2),
#     set_values_to = vars(
#       ANL01FL = 'Y'
#     )
#   )
adlb <- adlb %>%

  # Calculate AENTMTFL
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = vars(STUDYID, USUBJID,  PARAMCD),
      order = vars( VISITNUM, LBSEQ),
      new_var = AENTMTFL,
      mode = "last"
    ),
    filter = ( VISITNUM <= 12 & VISITNUM >= 2  & !is.na(SCHED))
  )


adlb <- adlb %>%
  # get MAXIMUM value
  derive_extreme_records(
    by_vars = vars(STUDYID, USUBJID, PARAMCD),
    order = vars( ADT, VISITNUM, LBSEQ),
    mode = "last",
    filter = ( VISITNUM >= 2 & VISITNUM <= 12 & !is.na(SCHED) ),
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
                R2A1HI,
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



adlbc_spec <-readxl::read_excel("metadata/specs.xlsx", sheet = "Variables")%>%
  filter(Dataset == "ADLBC")
adlbc_spec <- subset(adlbc_spec, select = c(Dataset, Variable, Label))
adlbc_spec <- rename(adlbc_spec, label=Label)
adlbc_spec <- rename(adlbc_spec, dataset=Dataset)
adlbc_spec <- rename(adlbc_spec, variable=Variable)
adlbc <- adlbc %>%
  xportr_label(adlbc_spec, domain = "ADLBC") # Assigns variable label from metacore specifications

 xportr_write(adlbc, "adam/adlbc.xpt")
