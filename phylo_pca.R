##### Phylogenetic PCA #####
# https://cran.r-project.org/web/packages/geomorph/vignettes/geomorph.PCA.html
# source header
source("Header.R")

# load in cleaned presence records
cleaned_pres = readr::read_rds("objects/cleaned_pres.rds")

# load in predictors
predictors = readr::read_rds("objects/predictors.rds")

# extract function
env_vals = lapply(cleaned_pres, get_env_vals, predictors)

env_vals_df = rbindlist(env_vals) |>
  as.data.frame()

# write rds
write_rds(env_vals_df, "objects/env_vals.rds")

##### Phylogenetic tree #####
# read in data
env_vals = readr::read_rds("objects/env_vals.rds")
cleaned_pres = readr::read_rds("objects/cleaned_pres.rds")

# species names
library(stringr)
species = names(cleaned_pres)
species_ = gsub(" ", "_", species)
row.names(env_vals) = species_
genus = stringr::word(species, 1)
family = c("Cyatheaceae", "Cyatheaceae", "Cyatheaceae", 
           "Cyatheaceae", "Dicksoniaceae", "Dicksoniaceae", 
           "Dicksoniaceae", "Loxsomataceae", "Cyatheaceae")
cyatheales = data.frame(species = species, genus = genus,
                        family = family)
# generate phylogenetic tree
library(V.PhyloMaker)
phylo_out = phylo.maker(cyatheales, scenarios=c("S1","S2","S3"))
par(mfrow = c(1, 1))
plot.phylo(phylo_out$scenario.1, cex = 1.5, main = "scenario.1")

phylo = phylo_out[[1]]

# phylogenetic pca
library(phytools)
data(anoletree)
data(anole.data)
## run phylogenetic PCA
pca<-phyl.pca(phylo, env_vals)
print(pca)
## plot results
plot(pca)
biplot(pca)
