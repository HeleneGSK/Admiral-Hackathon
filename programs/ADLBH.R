#Libraries
library(admiral)
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

# When SAS datasets are imported into R using haven::read_sas(), missing
# character values from SAS appear as "" characters in R, instead of appearing
# as NA values. Further details can be obtained via the following link:
# https://pharmaverse.github.io/admiral/articles/admiral.html#handling-of-missing-values

lb <- convert_blanks_to_na(lb)

unique(lb$LBCAT)

# lbh<-filter(lb,lbcat=="HEMATOLOGY")
lbh <- lb %>%
  filter(LBCAT == "HEMATOLOGY")

