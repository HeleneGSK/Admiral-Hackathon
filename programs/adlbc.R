library(haven)
library(admiral)
library(dplyr)
library(tidyr)
library(metacore)
library(metatools)
library(xportr)

lb <- read_xpt("sdtm/lb.xpt")
lb <- convert_blanks_to_na(lb)
lb1 <- lb %>%
 filter(LBCAT == "CHEMISTRY")


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

adsl <- read_xpt("adam/adsl.xpt")



use_ad_template(adlb)
xportr_write(adlbc, "adam/adlbc.xpt")
