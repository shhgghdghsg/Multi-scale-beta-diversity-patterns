library(tidyverse)
library(vegan)

load("clean.data/westsib.rda")

# Total number of species and subspecies
ncol(cdm)

# Samples contained 10 to 293 individuals (87.6 ± 51.1, mean ± sd)
cdm |> rowSums() |> range()
cdm |> rowSums() |> mean()
cdm |> rowSums() |> sd()

# Species richness varied from 1 to 26 species (15.5 ± 5.52). 
cdm |> specnumber() |> range()
cdm |> specnumber() |> mean()
cdm |> specnumber() |> sd()


# The proportion of Corythion dubium in the northern forest-tundra samples
sum(cdm$Corythion.dubium[1:18]) / sum(cdm[1:18,])

# The proportion of Corythion dubium in the forest-tundra samples
sum(cdm$Corythion.dubium[19:60]) / sum(cdm[19:60,])

# The proportion of Assulina muscorum in the forest-tundra samples
sum(cdm$Assulina.muscorum[19:60]) / sum(cdm[19:60,])

# The proportion of Assulina seminulum in the northern forest-tundra samples
sum(cdm$Assulina.seminulum[1:18]) / sum(cdm[1:18,])

# The proportion of Centropyxis aerophila in the northern forest-tundra samples
sum(cdm$Centropyxis.aerophila[1:18]) / sum(cdm[1:18,])  

# The proportion of Trinema lineare in the taiga samples
sum(cdm$Trinema.lineare[61:81]) / sum(cdm[61:81,])

# The proportion of Trinema lineare in the sub-taiga samples
sum(cdm$Trinema.lineare[82:111]) / sum(cdm[82:111,])


# The number of species observed per sub-zone
rowsum(cdm, meta$zone) |> specnumber()