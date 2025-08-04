library(tidyverse)

dat <- read.csv("raw.data/WestSib_TA-counts_2008.csv", header = T, sep = ";")

meta <- dat[, 1:12] |> select(sub.zone.code, ecosystem.code, biotope.code) |> 
  rename(zone = sub.zone.code, 
         veg = ecosystem.code,
         top = biotope.code) |>
  mutate(zone = factor(zone),
         zone = fct_relevel(zone, "ftu-n", "ftu", "tai", "stai"))

cdm <- dat[, 13:86]
cdm[is.na(cdm)] <- 0
cdm <- cdm[rowSums(cdm) > 0, colSums(cdm) > 0]

save(cdm, meta, file = "R/clean.data/westsib.rda")

